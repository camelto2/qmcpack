//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National
// Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National
// Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "QMCWaveFunctions/Fermion/PfaffianSTU.h"
#include "checkMatrix.hpp"
#include "ConstantSPOSet.h"

namespace qmcplusplus
{

namespace testing
{

class PfaffianSTUTest
{
  using ValueType   = PfaffianSTU::ValueType;
  using RealType    = PfaffianSTU::RealType;
  using ValueVector = PfaffianSTU::ValueVector;
  using ValueMatrix = PfaffianSTU::ValueMatrix;

public:
  std::unique_ptr<ParticleSet> createDummyElec(int nup, int ndn)
  {
    Lattice lattice;
    lattice.BoxBConds = true;
    lattice.R.diagonal(20);
    lattice.LR_dim_cutoff = 15;
    lattice.reset();

    const SimulationCell simulation_cell(lattice);
    auto elec = std::make_unique<ParticleSet>(simulation_cell);

    elec->setName("e");
    elec->create({nup, ndn});

    for (int ie = 0; ie < nup + ndn; ie++)
    {
      elec->R[ie] = {0.3 * ie + 0.1, 0.3 * ie + 0.2, 0.3 * ie + 0.3};
    }

    SpeciesSet& tspecies     = elec->getSpeciesSet();
    int upidx                = tspecies.addSpecies("u");
    int dnidx                = tspecies.addSpecies("d");
    int chgidx               = tspecies.addAttribute("charge");
    int massidx              = tspecies.addAttribute("mass");
    tspecies(chgidx, upidx)  = -1;
    tspecies(chgidx, dnidx)  = -1;
    tspecies(massidx, upidx) = 1;
    tspecies(massidx, dnidx) = 1;

    elec->createSK();
    elec->resetGroups();

    return elec;
  }

  std::unique_ptr<SPOSet> createDummySPO(const std::string name, const int nelec, const int norb) 
  {
    auto spo_ptr = std::make_unique<ConstantSPOSet<ValueType>>(name, nelec, norb);
    return spo_ptr;
  }

  void checkSizes(const int num_elec, const int uporbs, const int dnorbs, const PfaffianSTU& pf)
  {
    CHECK(num_elec == Approx(pf.num_elec_));
    if (num_elec % 2 == 0)
      CHECK(num_elec == Approx(pf.psi_mat_.rows()));
    else
      CHECK(num_elec + 1 == Approx(pf.psi_mat_.rows()));

    CHECK(uporbs*dnorbs == Approx(pf.singlet_mat_.size()));
    CHECK(uporbs*uporbs == Approx(pf.uu_triplet_mat_.size()));
    CHECK(dnorbs*dnorbs == Approx(pf.dd_triplet_mat_.size()));
  }

  void checkEvaluation(PfaffianSTU& pf, ParticleSet& elec)
  {
    //refernce values coming from brute force evaluation of pfaffians in python
    CHECK(pf.psi_mat_.rows() == Approx(6));
    ref_mat_.resize(6, 6);

    ValueVector row0 = {0., 0.18250294, -0.20363673, -0.25720828, -0.02146717, 0.06856466};
    ValueVector row1 = {-0.18250294, 0., 0.2460778, 0.0209037, 0.10986746, 0.3617668};
    ValueVector row2 = {0.20363673, -0.2460778, 0., 0.07461913, -0.01233308, -0.12637421};
    ValueVector row3 = {0.25720828, -0.0209037, -0.07461913, 0., -0.32933761, -0.09674686};
    ValueVector row4 = {0.02146717, -0.10986746, 0.01233308, 0.32933761, 0., 0.06171134};
    ValueVector row5 = {-0.06856466, -0.3617668, 0.12637421, 0.09674686, -0.06171134, 0.};
    for (int i = 0; i < 6; i++)
    {
      ref_mat_(0, i) = row0[i];
      ref_mat_(1, i) = row1[i];
      ref_mat_(2, i) = row2[i];
      ref_mat_(3, i) = row3[i];
      ref_mat_(4, i) = row4[i];
      ref_mat_(5, i) = row5[i];
    }

    pf.psi_mat_  = ref_mat_;
    RealType val = std::real(pf.calculatePfaffian());
    CHECK(val == Approx(-0.024797647365574358));

    //inverse of previous matrix
    row0 = {-2.71503989e-16, 1.81595609e+00, 4.32396207e+00, 9.92372604e-01, -2.35068908e-01, -2.92715884e+00};
    row1 = {-1.81595609e+00, -8.22199236e-17, -1.63444817e+00, 6.50271561e-01, -3.09990794e-01, -2.51198037e+00};
    row2 = {-4.32396207e+00, 1.63444817e+00, -2.22921635e-16, 1.07113524e+00, -3.09811991e+00, -1.30234146e+00};
    row3 = {-9.92372604e-01, -6.50271561e-01, -1.07113524e+00, 4.59511799e-17, 2.72112845e+00, -5.98429146e-01};
    row4 = {2.35068908e-01, 3.09990794e-01, 3.09811991e+00, -2.72112845e+00, -1.03859069e-16, -1.83155586e+00};
    row5 = {2.92715884e+00, 2.51198037e+00, 1.30234146e+00, 5.98429146e-01, 1.83155586e+00, -2.11660363e-17};
    for (int i = 0; i < 6; i++)
    {
      ref_mat_(0, i) = row0[i];
      ref_mat_(1, i) = row1[i];
      ref_mat_(2, i) = row2[i];
      ref_mat_(3, i) = row3[i];
      ref_mat_(4, i) = row4[i];
      ref_mat_(5, i) = row5[i];
    }

    pf.calculateInverse();
    auto check = checkMatrix(ref_mat_, pf.psi_matinv_);
    CHECKED_ELSE(check.result) { FAIL(check.result_message); }

    //new row for update
    //brute forcing update on particle 4
    const int iat = 4;
    pf.active_idx_     = iat;
    ValueVector newrow = {0.56637198, 0.69536227, 0.3739933,  0.69729468, 0.,         0.24455721};

    ValueType ratio = pf.calculateRatio(newrow);
    CHECK(std::real(ratio * val) == Approx(-0.020779936550488032));
  
    //this should update the inverse matrix after accepting proposed move
    pf.acceptMove(elec, iat);

    row0 = {-3.01041281e-15, 3.37758625e+00, 1.51309430e+01,-7.84382153e+00, -2.80518464e-01,-1.03781944e+01};
    row1 = {-3.37758625e+00,-1.78088739e-15,-7.96473187e+00, 7.07501492e+00, -3.69926172e-01,-1.70282163e-01};
    row2 = {-1.51309430e+01, 7.96473187e+00, 8.66730070e-16, 9.71368287e+00, -3.69712798e+00,-1.53008423e+01};
    row3 = { 7.84382153e+00,-7.07501492e+00,-9.71368287e+00,-2.52522350e-15,  3.24724686e+00, 1.68060062e+01};
    row4 = { 2.80518464e-01, 3.69926172e-01, 3.69712798e+00,-3.24724686e+00, -4.24653502e-16,-2.18567926e+00};
    row5 = { 1.03781944e+01, 1.70282163e-01, 1.53008423e+01,-1.68060062e+01,  2.18567926e+00,-8.15626218e-15};
    for (int i = 0; i < 6; i++)
    {
      ref_mat_(0, i) = row0[i];
      ref_mat_(1, i) = row1[i];
      ref_mat_(2, i) = row2[i];
      ref_mat_(3, i) = row3[i];
      ref_mat_(4, i) = row4[i];
      ref_mat_(5, i) = row5[i];
    }
    //ref_mat_ is inverse of updayed move
    check = checkMatrix(ref_mat_, pf.psi_matinv_);
    CHECKED_ELSE(check.result) { FAIL(check.result_message); }
  }

private:
  ValueMatrix ref_mat_;
};

} // namespace testing

TEST_CASE("Pfaffian check sizes", "[wavefunction][fermion]")
{
  Communicate* comm = OHMMS::Controller;
  testing::PfaffianSTUTest pftester;

  int nup = 4;
  int ndn = 2;
  int uporbs = 4;
  int dnorbs = 4;
  std::unique_ptr<ParticleSet> elec1 = pftester.createDummyElec(nup, ndn);
  std::unique_ptr<SPOSet> upspo = pftester.createDummySPO("up", nup, uporbs);
  std::unique_ptr<SPOSet> dnspo = pftester.createDummySPO("dn", ndn, dnorbs);
  std::vector<std::unique_ptr<SPOSet>> spos; 
  spos.emplace_back(std::move(upspo));
  spos.emplace_back(std::move(dnspo));
  PfaffianSTU pf1((*elec1), std::move(spos), "pfaffian1");
  pftester.checkSizes(elec1->getTotalNum(), uporbs, dnorbs, pf1);

  nup = 5;
  ndn = 4;
  uporbs = 7;
  dnorbs = 7;
  std::unique_ptr<ParticleSet> elec2 = pftester.createDummyElec(nup, ndn);
  std::unique_ptr<SPOSet> upspo2 = pftester.createDummySPO("up2", nup, uporbs);
  std::unique_ptr<SPOSet> dnspo2 = pftester.createDummySPO("dn2", ndn, dnorbs);
  std::vector<std::unique_ptr<SPOSet>> spos2;
  spos2.emplace_back(std::move(upspo2));
  spos2.emplace_back(std::move(dnspo2));
  PfaffianSTU pf2((*elec2), std::move(spos2), "pfaffian2");
  pftester.checkSizes(elec2->getTotalNum(), uporbs, dnorbs, pf2);
}

TEST_CASE("Pfaffian checkEvaluations", "[wavefunction][fermion]")
{
  Communicate* comm = OHMMS::Controller;
  testing::PfaffianSTUTest pftester;

  std::unique_ptr<ParticleSet> elec = pftester.createDummyElec(3, 3);
  std::unique_ptr<SPOSet> upspo = pftester.createDummySPO("up", 3, 3);
  std::unique_ptr<SPOSet> dnspo = pftester.createDummySPO("dn", 3, 3);
  std::vector<std::unique_ptr<SPOSet>> spos;
  spos.emplace_back(std::move(upspo));
  spos.emplace_back(std::move(dnspo));
  PfaffianSTU pf((*elec), std::move(spos), "pfaffian");
  pftester.checkEvaluation(pf, (*elec));
}

} // namespace qmcplusplus
