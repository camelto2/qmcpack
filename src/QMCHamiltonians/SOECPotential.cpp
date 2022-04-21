//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "Particle/DistanceTable.h"
#include "SOECPotential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
/** constructor
 *\param ionic positions
 *\param els electronic poitions
 *\param psi Trial wave function
*/
SOECPotential::SOECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi, bool compute_forces)
    : ForceBase(ions, els),
      myRNG(nullptr),
      IonConfig(ions),
      Psi(psi),
      compute_forces_(compute_forces),
      Peln(els),
      ElecNeighborIons(els),
      IonNeighborElecs(ions)
{
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(ions, els);
  myTableIndex = els.addTable(ions);
  NumIons      = ions.getTotalNum();
  PP.resize(NumIons, nullptr);
  PPset.resize(IonConfig.getSpeciesSet().getTotalNum());
  pulay_term_.resize(NumIons);
}

void SOECPotential::resetTargetParticleSet(ParticleSet& P) {}

SOECPotential::Return_t SOECPotential::evaluate(ParticleSet& P)
{
  value_ = 0.0;
  for (int ipp = 0; ipp < PPset.size(); ipp++)
    if (PPset[ipp])
      PPset[ipp]->randomize_grid(*myRNG);
  const auto& myTable = P.getDistTableAB(myTableIndex);
  for (int iat = 0; iat < NumIons; iat++)
    IonNeighborElecs.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons.getNeighborList(jel).clear();

  for (int jel = 0; jel < P.getTotalNum(); jel++)
  {
    const auto& dist               = myTable.getDistRow(jel);
    const auto& displ              = myTable.getDisplRow(jel);
    std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(jel);
    for (int iat = 0; iat < NumIons; iat++)
      if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
      {
        RealType pairpot = PP[iat]->evaluateOne(P, iat, Psi, jel, dist[iat], -displ[iat]);
        value_ += pairpot;
        NeighborIons.push_back(iat);
        IonNeighborElecs.getNeighborList(iat).push_back(jel);
      }
  }
  return value_;
}

void SOECPotential::evalIonDerivsImpl(ParticleSet& P,
                                      ParticleSet& ions,
                                      TrialWaveFunction& psi,
                                      ParticleSet::ParticlePos& hf_terms,
                                      ParticleSet::ParticlePos& pulay_terms,
                                      bool keep_grid)
{
  forces      = 0.0;
  pulay_term_ = 0.0;
  value_      = 0.0;
  if (!keep_grid)
  {
    for (int ipp = 0; ipp < PPset.size(); ipp++)
      if (PPset[ipp])
        PPset[ipp]->randomize_grid(*myRNG);
  }
  const auto& myTable = P.getDistTableAB(myTableIndex);
  for (int iat = 0; iat < NumIons; iat++)
    IonNeighborElecs.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons.getNeighborList(jel).clear();

  for (int ig = 0; ig < P.groups(); ig++)
  {
    Psi.prepareGroup(P, ig);
    for (int jel = P.first(ig); jel < P.last(ig); jel++)
    {
      const auto& dist               = myTable.getDistRow(jel);
      const auto& displ              = myTable.getDisplRow(jel);
      std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(jel);
      for (int iat = 0; iat < NumIons; iat++)
        if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
        {
          value_ +=
              PP[iat]->evaluateOneWithForces(P, ions, iat, Psi, jel, dist[iat], -displ[iat], forces[iat], pulay_term_);
          NeighborIons.push_back(iat);
          IonNeighborElecs.getNeighborList(iat).push_back(jel);
        }
    }
  }

  hf_terms -= forces;
  pulay_terms -= pulay_term_;
}

SOECPotential::Return_t SOECPotential::evaluateWithIonDerivs(ParticleSet& P,
                                                             ParticleSet& ions,
                                                             TrialWaveFunction& psi,
                                                             ParticleSet::ParticlePos& hf_terms,
                                                             ParticleSet::ParticlePos& pulay_terms)
{
  evalIonDerivsImpl(P, ions, psi, hf_terms, pulay_terms);
  return value_;
}

SOECPotential::Return_t SOECPotential::evaluateWithIonDerivsDeterministic(ParticleSet& P,
                                                                          ParticleSet& ions,
                                                                          TrialWaveFunction& psi,
                                                                          ParticleSet::ParticlePos& hf_terms,
                                                                          ParticleSet::ParticlePos& pulay_terms)
{
  evalIonDerivsImpl(P, ions, psi, hf_terms, pulay_terms, true);
  return value_;
}

std::unique_ptr<OperatorBase> SOECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<SOECPotential> myclone = std::make_unique<SOECPotential>(IonConfig, qp, psi, compute_forces_);
  for (int ig = 0; ig < PPset.size(); ++ig)
    if (PPset[ig])
      myclone->addComponent(ig, std::unique_ptr<SOECPComponent>(PPset[ig]->makeClone(qp)));
  return myclone;
}

void SOECPotential::addComponent(int groupID, std::unique_ptr<SOECPComponent>&& ppot)
{
  for (int iat = 0; iat < PP.size(); iat++)
    if (IonConfig.GroupID[iat] == groupID)
      PP[iat] = ppot.get();
  PPset[groupID] = std::move(ppot);
}

} // namespace qmcplusplus
