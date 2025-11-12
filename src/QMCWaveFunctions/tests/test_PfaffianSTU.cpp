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

namespace qmcplusplus
{

namespace testing
{

class PfaffianSTUTest
{
  using ValueType   = PfaffianSTU::ValueType;
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

  void checkSizes(const int num_elec, const PfaffianSTU& pf)
  {
    CHECK(num_elec == Approx(pf.num_elec_));
    if (num_elec%2 == 0)
      CHECK(num_elec == Approx(pf.psi_mat_.rows()));
    else
      CHECK(num_elec+1 == Approx(pf.psi_mat_.rows()));
  }

  void checkEvaluation(PfaffianSTU& pf)
  {
    CHECK(pf.psi_mat_.rows() == Approx(6));
    ref_mat_.resize(6,6);
    ValueVector row0 = { 0.,         -0.0476868 , -0.04611147, 0.00385507, -0.04720784,  0.31010055};
    ValueVector row1 = { 0.0476868 ,  0.        , -0.1916138 , 0.07984651,  0.10802211,  0.0899957 };
    ValueVector row2 = { 0.04611147,  0.1916138 ,  0.        , 0.04312142, -0.17709148,  0.44741451};
    ValueVector row3 = {-0.00385507, -0.07984651, -0.04312142, 0.        , -0.24264741, -0.247939  };
    ValueVector row4 = { 0.04720784, -0.10802211,  0.17709148, 0.24264741,  0.        ,  0.09156629};
    ValueVector row5 = {-0.31010055, -0.0899957 , -0.44741451, 0.247939  , -0.09156629,  0.        };
    for (int i = 0; i < 6; i++)
    {
      ref_mat_(0,i) = row0[i];
      ref_mat_(1,i) = row1[i];
      ref_mat_(2,i) = row2[i];
      ref_mat_(3,i) = row3[i];
      ref_mat_(4,i) = row4[i];
      ref_mat_(5,i) = row5[i];
    }

    pf.psi_mat_ = ref_mat_;
    CHECK(std::real(pf.calculatePfaffian()) == Approx(0.028319226149595082));
  }


private:
  ValueMatrix ref_mat_;
};

} // namespace testing

TEST_CASE("Pfaffian check sizes", "[wavefunction][fermion]")
{
  Communicate* comm = OHMMS::Controller;
  testing::PfaffianSTUTest pftester;

  std::unique_ptr<ParticleSet> elec1 = pftester.createDummyElec(4, 2);
  PfaffianSTU pf1((*elec1), "pfaffian1");
  std::unique_ptr<ParticleSet> elec2 = pftester.createDummyElec(5, 4);
  PfaffianSTU pf2((*elec2), "pfaffian2");

  pftester.checkSizes(elec1->getTotalNum(), pf1);
  pftester.checkSizes(elec2->getTotalNum(), pf2);
}

TEST_CASE("Pfaffian evaluate", "[wavefunction][fermion]")
{
  Communicate* comm = OHMMS::Controller;
  testing::PfaffianSTUTest pftester;

  std::unique_ptr<ParticleSet> elec1 = pftester.createDummyElec(3, 3);
  PfaffianSTU pf((*elec1), "pfaffian");
  pftester.checkEvaluation(pf);

}

} // namespace qmcplusplus
