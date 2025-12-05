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

    CHECK(uporbs * dnorbs == Approx(pf.singlet_mat_.size()));
    CHECK(uporbs * uporbs == Approx(pf.uu_triplet_mat_.size()));
    CHECK(dnorbs * dnorbs == Approx(pf.dd_triplet_mat_.size()));
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
    const int iat      = 4;
    pf.active_idx_     = iat;
    ValueVector newrow = {0.56637198, 0.69536227, 0.3739933, 0.69729468, 0., 0.24455721};

    ValueType ratio = pf.calculateRatio(newrow);
    CHECK(std::real(ratio * val) == Approx(-0.020779936550488032));

    //this should update the inverse matrix after accepting proposed move
    pf.acceptMove(elec, iat);

    row0 = {-3.01041281e-15, 3.37758625e+00, 1.51309430e+01, -7.84382153e+00, -2.80518464e-01, -1.03781944e+01};
    row1 = {-3.37758625e+00, -1.78088739e-15, -7.96473187e+00, 7.07501492e+00, -3.69926172e-01, -1.70282163e-01};
    row2 = {-1.51309430e+01, 7.96473187e+00, 8.66730070e-16, 9.71368287e+00, -3.69712798e+00, -1.53008423e+01};
    row3 = {7.84382153e+00, -7.07501492e+00, -9.71368287e+00, -2.52522350e-15, 3.24724686e+00, 1.68060062e+01};
    row4 = {2.80518464e-01, 3.69926172e-01, 3.69712798e+00, -3.24724686e+00, -4.24653502e-16, -2.18567926e+00};
    row5 = {1.03781944e+01, 1.70282163e-01, 1.53008423e+01, -1.68060062e+01, 2.18567926e+00, -8.15626218e-15};
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

  void checkLog(PfaffianSTU& pf, ParticleSet& elec)
  {
    //now going to test how the pfaffian is evaluated given random up/dn SPOSet values and random pairing matrices
    up_ref_mat_.resize(pf.up_psi_mat_.rows(), pf.up_psi_mat_.cols());
    dn_ref_mat_.resize(pf.dn_psi_mat_.rows(), pf.dn_psi_mat_.cols());
    ref_mat_.resize(pf.psi_mat_.rows(), pf.psi_mat_.cols());

    ValueVector row0 = {0.92979771, 0.04481706, 0.36345255, 0.55100725, 0.11389361, 0.83844704};
    ValueVector row1 = {0.8024512, 0.47750522, 0.59375306, 0.12799644, 0.56722092, 0.13476818};
    ValueVector row2 = {0.81711275, 0.82519555, 0.78641738, 0.12231961, 0.52508536, 0.2608769};
    for (int i = 0; i < up_ref_mat_.cols(); i++)
    {
      up_ref_mat_[0][i] = row0[i];
      up_ref_mat_[1][i] = row1[i];
      up_ref_mat_[2][i] = row2[i];
    }
    auto upspo = dynamic_cast<ConstantSPOSet<ValueType>*>(pf.sposets_[0].get());
    upspo->setRefVals(up_ref_mat_);

    row0 = {0.11087974, 0.91875922, 0.70049105, 0.1151218, 0.49077789, 0.36205635};
    row1 = {0.92619697, 0.83850844, 0.38222666, 0.69600158, 0.10640622, 0.1246681};
    for (int i = 0; i < dn_ref_mat_.cols(); i++)
    {
      dn_ref_mat_[0][i] = row0[i];
      dn_ref_mat_[1][i] = row1[i];
    }
    auto dnspo = dynamic_cast<ConstantSPOSet<ValueType>*>(pf.sposets_[1].get());
    dnspo->setRefVals(dn_ref_mat_);

    //singlet
    ValueVector row3, row4, row5;
    row0 = {0.3180414592166366, 0.600104959392495,  0.4642438313743457,
            0.5816172764398225, 0.6037109397256399, 0.2951113484296276};
    row1 = {0.600104959392495,  0.09697842311870486, 0.41937315764161337,
            0.5421635350896339, 0.2548558364689513,  0.4271829154212393};
    row2 = {0.4642438313743457, 0.41937315764161337, 0.6797705794844521,
            0.6307981709790462, 0.28373386258295497, 0.15635159921667535};
    row3 = {0.5816172764398225, 0.5421635350896339,  0.6307981709790462,
            0.8418923779517383, 0.44010957459901245, 0.14809503193049156};
    row4 = {0.6037109397256399,  0.2548558364689513,  0.28373386258295497,
            0.44010957459901245, 0.17414766814937777, 0.6661621789790062};
    row5 = {0.2951113484296276,  0.4271829154212393, 0.15635159921667535,
            0.14809503193049156, 0.6661621789790062, 0.7995336027593793};
    for (int i = 0; i < pf.singlet_mat_.cols(); i++)
    {
      pf.singlet_mat_(0, i) = row0[i];
      pf.singlet_mat_(1, i) = row1[i];
      pf.singlet_mat_(2, i) = row2[i];
      pf.singlet_mat_(3, i) = row3[i];
      pf.singlet_mat_(4, i) = row4[i];
      pf.singlet_mat_(5, i) = row5[i];
    }
    
    //UU triplet
    row0 = { 0.0 , 0.10179712049738304 , 0.04720517011113051 , 0.13767202491931702 , -0.04487759834743943 , 0.1305860299961965 };
    row1 = { -0.10179712049738304 , 0.0 , -0.060719691450885904 , -0.1700387281327197 , 0.1505147613342654 , -0.1258821272345353 };
    row2 = { -0.04720517011113051 , 0.060719691450885904 , 0.0 , -0.17528784342592918 , 0.18023264504075892 , -0.2124931747481953 };
    row3 = { -0.13767202491931702 , 0.1700387281327197 , 0.17528784342592918 , 0.0 , 0.011885828639267404 , -0.2724282490725342 };
    row4 = { 0.04487759834743943 , -0.1505147613342654 , -0.18023264504075892 , -0.011885828639267404 , 0.0 , -0.2527551842515271 };
    row5 = { -0.1305860299961965 , 0.1258821272345353 , 0.2124931747481953 , 0.2724282490725342 , 0.2527551842515271 , 0.0 };

    for (int i = 0; i < pf.uu_triplet_mat_.cols(); i++)
    {
      pf.uu_triplet_mat_(0, i) = row0[i];
      pf.uu_triplet_mat_(1, i) = row1[i];
      pf.uu_triplet_mat_(2, i) = row2[i];
      pf.uu_triplet_mat_(3, i) = row3[i];
      pf.uu_triplet_mat_(4, i) = row4[i];
      pf.uu_triplet_mat_(5, i) = row5[i];
    }

    //DD triplet
    row0 = { 0.0 , -0.2841971014610199 , -0.22729801265851418 , 0.02689802638369143 , -0.014699398554927245 , -0.2408488597164405 };
    row1 = { 0.2841971014610199 , 0.0 , -0.11422553793697526 , -0.09290196456215388 , 0.0845644076456245 , -0.004046219453199551 };
    row2 = { 0.22729801265851418 , 0.11422553793697526 , 0.0 , 0.18361494551113666 , 0.04352167620130676 , 0.19876007777434612 };
    row3 = { -0.02689802638369143 , 0.09290196456215388 , -0.18361494551113666 , 0.0 , 0.4455646986311197 , 0.42441739601566725 };
    row4 = { 0.014699398554927245 , -0.0845644076456245 , -0.04352167620130676 , -0.4455646986311197 , 0.0 , 0.18807375032492452 };
    row5 = { 0.2408488597164405 , 0.004046219453199551 , -0.19876007777434612 , -0.42441739601566725 , -0.18807375032492452 , 0.0 };
    for (int i = 0; i < pf.dd_triplet_mat_.cols(); i++)
    {
      pf.dd_triplet_mat_(0, i) = row0[i];
      pf.dd_triplet_mat_(1, i) = row1[i];
      pf.dd_triplet_mat_(2, i) = row2[i];
      pf.dd_triplet_mat_(3, i) = row3[i];
      pf.dd_triplet_mat_(4, i) = row4[i];
      pf.dd_triplet_mat_(5, i) = row5[i];
    }

    ParticleSet::ParticleGradient G;
    ParticleSet::ParticleLaplacian L;
    pf.evaluateLog(elec, G, L);
    ValueType ref_pfaff = -0.5520438346362099;
    CHECK(std::log(std::abs(ref_pfaff)) == Approx(std::real(pf.log_value_)));
    CHECK(std::arg(ref_pfaff) == Approx(std::imag(pf.log_value_)));
  }

private:
  ValueMatrix ref_mat_;
  ValueMatrix up_ref_mat_;
  ValueMatrix dn_ref_mat_;
  ValueMatrix tmp_pairing_;
};

} // namespace testing

TEST_CASE("Pfaffian check sizes", "[wavefunction][fermion]")
{
  Communicate* comm = OHMMS::Controller;
  testing::PfaffianSTUTest pftester;

  int nup                            = 4;
  int ndn                            = 2;
  int uporbs                         = 4;
  int dnorbs                         = 4;
  std::unique_ptr<ParticleSet> elec1 = pftester.createDummyElec(nup, ndn);
  std::unique_ptr<SPOSet> upspo      = pftester.createDummySPO("up", nup, uporbs);
  std::unique_ptr<SPOSet> dnspo      = pftester.createDummySPO("dn", ndn, dnorbs);
  std::vector<std::unique_ptr<SPOSet>> spos;
  spos.emplace_back(std::move(upspo));
  spos.emplace_back(std::move(dnspo));
  PfaffianSTU pf1((*elec1), std::move(spos), "pfaffian1");
  pftester.checkSizes(elec1->getTotalNum(), uporbs, dnorbs, pf1);

  nup                                = 5;
  ndn                                = 4;
  uporbs                             = 7;
  dnorbs                             = 7;
  std::unique_ptr<ParticleSet> elec2 = pftester.createDummyElec(nup, ndn);
  std::unique_ptr<SPOSet> upspo2     = pftester.createDummySPO("up2", nup, uporbs);
  std::unique_ptr<SPOSet> dnspo2     = pftester.createDummySPO("dn2", ndn, dnorbs);
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
  std::unique_ptr<SPOSet> upspo     = pftester.createDummySPO("up", 3, 3);
  std::unique_ptr<SPOSet> dnspo     = pftester.createDummySPO("dn", 3, 3);
  std::vector<std::unique_ptr<SPOSet>> spos;
  spos.emplace_back(std::move(upspo));
  spos.emplace_back(std::move(dnspo));
  PfaffianSTU pf((*elec), std::move(spos), "pfaffian");
  pftester.checkEvaluation(pf, (*elec));
}

TEST_CASE("Pfaffian evaluateLog", "[wavefunction][fermion]")
{
  Communicate* comm = OHMMS::Controller;
  testing::PfaffianSTUTest pftester;

  int nup                           = 3;
  int ndn                           = 2;
  int norb                          = 6;
  std::unique_ptr<ParticleSet> elec = pftester.createDummyElec(nup, ndn);
  std::unique_ptr<SPOSet> upspo     = pftester.createDummySPO("up", nup, norb);
  std::unique_ptr<SPOSet> dnspo     = pftester.createDummySPO("dn", ndn, norb);
  std::vector<std::unique_ptr<SPOSet>> spos;
  spos.emplace_back(std::move(upspo));
  spos.emplace_back(std::move(dnspo));
  PfaffianSTU pf((*elec), std::move(spos), "pfaffian");
  pftester.checkLog(pf, (*elec));
}

} // namespace qmcplusplus
