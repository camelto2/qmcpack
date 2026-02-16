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
#include "io/hdf/hdf_archive.h"

namespace qmcplusplus
{

namespace testing
{

class PfaffianSTUTest
{
  using ValueType   = PfaffianSTU::ValueType;
  using GradType    = PfaffianSTU::GradType;
  using RealType    = PfaffianSTU::RealType;
  using ValueVector = PfaffianSTU::ValueVector;
  using GradVector  = PfaffianSTU::GradVector;
  using ValueMatrix = PfaffianSTU::ValueMatrix;
  using GradMatrix  = PfaffianSTU::GradMatrix;

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

    pf.psi_mat_   = ref_mat_;
    ValueType val = pf.calculatePfaffian();
    CHECK(std::real(val) == Approx(-0.024797647365574358));

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
    pf.cur_ratio_ =
        ratio; //this normally happens in ratio/ratioGrad, but I'm directly calling  calculateRatio so need this so accceptMove passes since it checks this
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

  void checkVGLEvaluations(PfaffianSTU& pf, ParticleSet& elec)
  {
    //Check evaluateLog for Pfaffian, given a set of up and down orbitals.
    //the orbitals used are sin/cos(k.r) with random k values and electron positions
    //All reference data can be found in pfaffian_testing.ipynb
    up_ref_mat_.resize(pf.up_psi_mat_.rows(), pf.up_psi_mat_.cols());
    dn_ref_mat_.resize(pf.dn_psi_mat_.rows(), pf.dn_psi_mat_.cols());
    up_dref_mat_.resize(pf.up_dpsi_mat_.rows(), pf.up_dpsi_mat_.cols());
    dn_dref_mat_.resize(pf.dn_dpsi_mat_.rows(), pf.dn_dpsi_mat_.cols());
    up_d2ref_mat_.resize(pf.up_d2psi_mat_.rows(), pf.up_d2psi_mat_.cols());
    dn_d2ref_mat_.resize(pf.dn_d2psi_mat_.rows(), pf.dn_d2psi_mat_.cols());

    ref_mat_.resize(pf.psi_mat_.rows(), pf.psi_mat_.cols());

    //up values
    ValueVector row0 = {0.26880981, 0.38297036, 0.67970534, 0.57913668, 0.34154767, 0.68877715};
    ValueVector row1 = {-0.00249844, 0.07922684, 0.51539003, 0.39623433, 0.10879198, 0.5439395};
    ValueVector row2 = {0.38216715, 0.31236072, 0.70069793, 0.55406334, 0.36763894, 0.64983887};
    for (int i = 0; i < up_ref_mat_.cols(); i++)
    {
      up_ref_mat_[0][i] = row0[i];
      up_ref_mat_[1][i] = row1[i];
      up_ref_mat_[2][i] = row2[i];
    }
    auto upspo = dynamic_cast<ConstantSPOSet<ValueType>*>(pf.sposets_[0].get());
    upspo->setRefVals(up_ref_mat_);

    //up grads
    // clang-format off
    GradVector g0 = {{-0.64400252,-0.91610297,-0.48422924},
                     {-0.44093363,-0.45549765,-0.77629632},
                     {-0.37069611,-0.38923311,-0.23519604},
                     {-0.13093672,-0.33648282,-0.6720715 },
                     {-0.19117277,-0.61472653,-0.90141821},
                     {-0.06633773,-0.21230168,-0.56491323}};
    GradVector g1 = {{-0.66860985,-0.95110725,-0.50273163}, 
                     {-0.47582413,-0.49154059,-0.83772365},
                     {-0.43309684,-0.45475425,-0.27478751},
                     {-0.14746686,-0.37896217,-0.75691732},
                     {-0.20219732,-0.65017657,-0.95340118},
                     {-0.076783  ,-0.24572984,-0.65386216}};
    GradVector g2 = {{-0.61785976,-0.8789145 ,-0.46457234},
                     {-0.45344093,-0.46841808,-0.79831636}, 
                     {-0.36057456,-0.37860542,-0.22877421}, 
                     {-0.13370634,-0.34360023,-0.6862874 },
                     {-0.18915991,-0.60825407,-0.89192716},
                     {-0.0695494 ,-0.22258006,-0.59226296}};
    // clang-format on
    for (int i = 0; i < up_dref_mat_.cols(); i++)
    {
      up_dref_mat_[0][i] = g0[i];
      up_dref_mat_[1][i] = g1[i];
      up_dref_mat_[2][i] = g2[i];
    }
    upspo->setRefEGrads(up_dref_mat_);

    //up laps
    // clang-format off
    row0 = {-0.43127665,-0.45082918,-0.43490375,-0.50719738,-0.47441901,-0.4830475 };
    row1 = { 0.00400849,-0.09326511,-0.32976798,-0.34701482,-0.15111502,-0.38147116};
    row2 = {-0.61314639,-0.36770815,-0.44833568,-0.4852386 ,-0.5106605 ,-0.45573962};
    // clang-format on
    for (int i = 0; i < up_d2ref_mat_.cols(); i++)
    {
      up_d2ref_mat_[0][i] = row0[i];
      up_d2ref_mat_[1][i] = row1[i];
      up_d2ref_mat_[2][i] = row2[i];
    }
    upspo->setRefELapls(up_d2ref_mat_);

    //down values
    // clang-format off
    row0 = { 0.6851893 , 0.60400179, 0.89348619, 0.91161367, 0.86358127, 0.94059632};
    row1 = {-0.46062194,-0.66440113, 0.57008769, 0.56464611, 0.02793061, 0.42318932};
    // clang-format on
    for (int i = 0; i < dn_ref_mat_.cols(); i++)
    {
      dn_ref_mat_[0][i] = row0[i];
      dn_ref_mat_[1][i] = row1[i];
    }
    auto dnspo = dynamic_cast<ConstantSPOSet<ValueType>*>(pf.sposets_[1].get());
    dnspo->setRefVals(dn_ref_mat_);

    //down grads
    // clang-format off
    g0 = {{-0.62641093,-0.57186363,-0.54915937}, 
          {-0.75098121,-0.71416503,-0.67599978},
          {-0.06534126,-0.31283256,-0.10731589},
          {-0.14939149,-0.02329714,-0.30623787},
          {-0.40796331,-0.47395767,-0.01600182},
          {-0.27075382,-0.02753156,-0.1706383 }};
    g1 = {{-0.76335365,-0.6968815 ,-0.66921376},
          {-0.70423772,-0.6697131 ,-0.63392338}, 
          {-0.11953783,-0.572308  ,-0.19632784},
          {-0.29995943,-0.04677775,-0.61488735},
          {-0.80879875,-0.93963444,-0.03172407},
          {-0.72251736,-0.07346907,-0.45535511}};
    // clang-format on
    for (int i = 0; i < dn_dref_mat_.cols(); i++)
    {
      dn_dref_mat_[0][i] = g0[i];
      dn_dref_mat_[1][i] = g1[i];
    }
    dnspo->setRefEGrads(dn_dref_mat_);

    //down laps
    // clang-format off
    row0 = {-1.31866924,-1.45582683,-0.50348985,-0.62933418,-1.32928953,-0.84190492};
    row1 = { 0.88648201, 1.60140749,-0.32125104,-0.38980449,-0.0429929 ,-0.37878648};
    // clang-format on
    for (int i = 0; i < dn_d2ref_mat_.cols(); i++)
    {
      dn_d2ref_mat_[0][i] = row0[i];
      dn_d2ref_mat_[1][i] = row1[i];
    }
    dnspo->setRefELapls(dn_d2ref_mat_);

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
    row0 =
        {0.0, 0.10179712049738304, 0.04720517011113051, 0.13767202491931702, -0.04487759834743943, 0.1305860299961965};
    row1 = {-0.10179712049738304, 0.0, -0.060719691450885904, -0.1700387281327197, 0.1505147613342654,
            -0.1258821272345353};
    row2 = {-0.04720517011113051, 0.060719691450885904, 0.0,
            -0.17528784342592918, 0.18023264504075892,  -0.2124931747481953};
    row3 = {-0.13767202491931702, 0.1700387281327197, 0.17528784342592918, 0.0,
            0.011885828639267404, -0.2724282490725342};
    row4 = {0.04487759834743943, -0.1505147613342654, -0.18023264504075892, -0.011885828639267404, 0.0,
            -0.2527551842515271};
    row5 = {-0.1305860299961965, 0.1258821272345353, 0.2124931747481953, 0.2724282490725342, 0.2527551842515271, 0.0};

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
    row0 = {0.000000000000,      -0.2841971014610199,   -0.22729801265851418,
            0.02689802638369143, -0.014699398554927245, -0.2408488597164405};
    row1 = {0.2841971014610199,   0.0, -0.11422553793697526, -0.09290196456215388, 0.0845644076456245,
            -0.004046219453199551};
    row2 = {0.22729801265851418, 0.11422553793697526, 0.0,
            0.18361494551113666, 0.04352167620130676, 0.19876007777434612};
    row3 = {-0.02689802638369143, 0.09290196456215388, -0.18361494551113666, 0.0,
            0.4455646986311197,   0.42441739601566725};
    row4 = {0.014699398554927245, -0.0845644076456245, -0.04352167620130676, -0.4455646986311197, 0.0,
            0.18807375032492452};
    row5 = {0.2408488597164405,   0.004046219453199551, -0.19876007777434612,
            -0.42441739601566725, -0.18807375032492452, 0.0};
    for (int i = 0; i < pf.dd_triplet_mat_.cols(); i++)
    {
      pf.dd_triplet_mat_(0, i) = row0[i];
      pf.dd_triplet_mat_(1, i) = row1[i];
      pf.dd_triplet_mat_(2, i) = row2[i];
      pf.dd_triplet_mat_(3, i) = row3[i];
      pf.dd_triplet_mat_(4, i) = row4[i];
      pf.dd_triplet_mat_(5, i) = row5[i];
    }
    pf.updatePairingNorms();

    ParticleSet::ParticleGradient G(pf.num_elec_);
    ParticleSet::ParticleLaplacian L(pf.num_elec_);
    pf.evaluateLog(elec, G, L);
    ValueType ref_pfaff = 0.010128569344188544;
    CHECK(std::log(std::abs(ref_pfaff)) == Approx(std::real(pf.log_value_)));
    CHECK(std::arg(ref_pfaff) == Approx(std::imag(pf.log_value_)));

    //These reference values are from finite differences
    ParticleSet::ParticleGradient Gref(pf.num_elec_);
    ParticleSet::ParticleLaplacian Lref(pf.num_elec_);
    Gref[0] = {1.69155466, 0.57515117, -6.51474671};
    Lref[0] = -47.03689953825609;
    Gref[1] = {1.39715297, 2.55443726, 3.27057718};
    Lref[1] = -19.380044871245367;
    Gref[2] = {-2.56733439, -2.45232092, 3.30552091};
    Lref[2] = -25.485185886646185;
    Gref[3] = {-0.52971275, -0.47760647, -0.43386478};
    Lref[3] = -2.033374656205245;
    Gref[4] = {0.68774557, 0.2730102, 0.43364048};
    Lref[4] = -2.3526677861768284;
    for (int i = 0; i < pf.num_elec_; i++)
    {
      CHECK(G[i][0] == ValueApprox(Gref[i][0]));
      CHECK(G[i][1] == ValueApprox(Gref[i][1]));
      CHECK(G[i][2] == ValueApprox(Gref[i][2]));
      CHECK(L[i] == ValueApprox(Lref[i]));
    }

    //check evalGrad
    for (int i = 0; i < pf.num_elec_; i++)
    {
      GradType grad = pf.evalGrad(elec, i);
      CHECK(grad[0] == ValueApprox(Gref[i][0]));
      CHECK(grad[1] == ValueApprox(Gref[i][1]));
      CHECK(grad[2] == ValueApprox(Gref[i][2]));
    }

    //now we will move particle 1 to test ratio()
    //note that since I'm working with ConstantSPOSet, I don't actually have to move the electron in
    //the particle set.
    //first need to update the SPOSet under the hood so evaluateValue returns the correct data.
    int iat          = 1;
    ValueVector newv = {0.07497445, 0.21713796, 0.57738733, 0.47357094, 0.18976185, 0.61038528};
    upspo->updateV(elec, iat, newv);

    ValueType ref_ratio = 0.619845148928336;
    ValueType ratio     = pf.ratio(elec, iat);
    CHECK(ratio == ValueApprox(ref_ratio));
    //set values back to original since not accepting move
    newv = {-0.00249844, 0.07922684, 0.51539003, 0.39623433, 0.10879198, 0.5439395};
    upspo->updateV(elec, iat, newv);


    //try new particle and check ratioGrad
    iat = 4;
    // clang-format off
    newv             = {-0.03079311, -0.22518306,  0.67066673,  0.80665308,  0.20567156,  0.6938236};
    GradVector newg  = {{-0.8596155 ,-0.78476096,-0.7536042 },
                        {-0.91807906,-0.87307106,-0.82641381}, 
                        {-0.10792355,-0.51670259,-0.17725264},
                        {-0.21480247,-0.03349778,-0.44032393},
                        {-0.79181642,-0.91990496,-0.03105796},
                        {-0.57427544,-0.05839511,-0.361928  }};
    ValueVector newl = {0.05926236,  0.54275921, -0.3779285,  -0.55687445, -0.3165852,  -0.62102466};
    // clang-format on
    dnspo->updateVGL(elec, iat, newv, newg, newl);
    ref_ratio         = 0.6432811011284547;
    GradType ref_grad = {1.33729694, 0.73665635, 1.02832973};
    GradType grad;
    ratio = pf.ratioGrad(elec, iat, grad);
    CHECK(ratio == ValueApprox(ref_ratio));
    for (int d = 0; d < 3; d++)
      CHECK(grad[d] == ValueApprox(ref_grad[d]));

    //accept Move and lets run evaluateLog to make sure gradient gets updated properly with evalGrad
    pf.acceptMove(elec, iat);
    grad = pf.evalGrad(elec, iat);
    for (int d = 0; d < 3; d++)
      CHECK(grad[d] == ValueApprox(ref_grad[d]));
  }

  void checkReadWrite(PfaffianSTU& pf, ParticleSet& elec)
  {
    pf.buildOptVariables();

    const int num_orbs   = pf.sposets_[0]->size();
    const int num_params = num_orbs * (num_orbs + 1) / 2 + num_orbs * (num_orbs - 1);
    CHECK(num_params == pf.myVars.size());

    std::vector<RealType> parms(num_params);
    for (int i = 0; i < num_params; i++)
      parms[i] = 0.1 * i;

    optimize::VariableSet vs;
    pf.checkInVariablesExclusive(vs);
    for (size_t i = 0; i < vs.size(); i++)
      vs[i] = parms[i];
    pf.resetParametersExclusive(vs);

    hdf_archive hout;
    vs.writeToHDF("pf_vp.h5", hout);

    pf.initializePairingMats(); //reset

    vs.readFromHDF("pf_vp.h5", hout);
    pf.resetParametersExclusive(vs);

    int count = 0;
    for (int i = 0; i < pf.singlet_mat_.rows(); i++)
    {
      CHECK(pf.singlet_mat_(i, i) == ValueApprox(vs[count++]));
      for (int j = i + 1; j < pf.singlet_mat_.rows(); j++, count++)
      {
        CHECK(pf.singlet_mat_(i, j) == ValueApprox(vs[count]));
        CHECK(pf.singlet_mat_(j, i) == ValueApprox(vs[count]));
      }
    }
    for (int i = 0; i < pf.uu_triplet_mat_.rows(); i++)
    {
      for (int j = i + 1; j < pf.uu_triplet_mat_.rows(); j++, count++)
      {
        CHECK(pf.uu_triplet_mat_(i, j) == ValueApprox(vs[count]));
        CHECK(pf.uu_triplet_mat_(j, i) == ValueApprox(-vs[count]));
      }
    }
    for (int i = 0; i < pf.dd_triplet_mat_.rows(); i++)
    {
      for (int j = i + 1; j < pf.dd_triplet_mat_.rows(); j++, count++)
      {
        CHECK(pf.dd_triplet_mat_(i, j) == ValueApprox(vs[count]));
        CHECK(pf.dd_triplet_mat_(j, i) == ValueApprox(-vs[count]));
      }
    }
  }

  void checkParmDerivs(PfaffianSTU& pf, ParticleSet& elec)
  {
    pf.buildOptVariables();

    //set ref vals, for orbitals
    up_ref_mat_.resize(pf.up_psi_mat_.rows(), pf.up_psi_mat_.cols());
    dn_ref_mat_.resize(pf.dn_psi_mat_.rows(), pf.dn_psi_mat_.cols());

    ValueVector row0, row1;
    // clang-format off
    row0 = {0.46170647,  0.64069846,  0.70417798};
    row1 = {0.82174694,  0.88964088,  0.82151254};
    // clang-format on
    for (int i = 0; i < up_ref_mat_.cols(); i++)
    {
      up_ref_mat_[0][i] = row0[i];
      up_ref_mat_[1][i] = row1[i];
    }
    auto upspo = dynamic_cast<ConstantSPOSet<ValueType>*>(pf.sposets_[0].get());
    upspo->setRefVals(up_ref_mat_);

    // clang-format off
    row0 = {0.80904975,  0.86552841,  0.89814096};
    row1 = {0.79386346,  0.91672375,  0.75376951};
    // clang-format on
    for (int i = 0; i < dn_ref_mat_.cols(); i++)
    {
      dn_ref_mat_[0][i] = row0[i];
      dn_ref_mat_[1][i] = row1[i];
    }
    auto dnspo = dynamic_cast<ConstantSPOSet<ValueType>*>(pf.sposets_[1].get());
    dnspo->setRefVals(dn_ref_mat_);

    //parameters to use for S, UU, DD mats
    //set them via resetParametersExclustive


    //from finite differences
    std::vector<ValueType> ref_derivs{-0.15113309767512936, 0.4773609882285376,  -0.7022153622589955,
                                      -0.08069831730561128, -0.7092918379802287, 0.3725852855472653,
                                      -4.615645636943915,   -5.115722397967572,  0.3964642822551843,
                                      -1.318567401160748,   1.2593480402595423,  1.515855381743893};


    //independent values of S, UU, DD matrix. use checkInVariablesExclusive to set matrices
    std::vector<RealType> parms{0.5981391841682085,   0.6942599356718862, 0.04132254367265181,    0.34934851353450946,
                                0.6610527560659689,   0.765125528380152,  0.020775902403176283,   -0.024598622514295176,
                                -0.07553123838606379, 0.2754793798265398, -0.0004650646646693901, 0.2400122192700468};

    optimize::VariableSet vs;
    pf.checkInVariablesExclusive(vs);
    for (size_t i = 0; i < vs.size(); i++)
      vs[i] = parms[i];
    pf.resetParametersExclusive(vs);

    Vector<ValueType> dlogpsi(vs.size());

    pf.evaluateDerivativesWF(elec, vs, dlogpsi);
    for (int i = 0; i < dlogpsi.size(); i++)
      CHECK(dlogpsi[i] == ValueApprox(ref_derivs[i]));
  }

private:
  ValueMatrix ref_mat_;
  ValueMatrix up_ref_mat_;
  ValueMatrix dn_ref_mat_;
  GradMatrix up_dref_mat_;
  GradMatrix dn_dref_mat_;
  ValueMatrix up_d2ref_mat_;
  ValueMatrix dn_d2ref_mat_;
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
  const std::string opt              = "no";
  std::unique_ptr<ParticleSet> elec1 = pftester.createDummyElec(nup, ndn);
  std::unique_ptr<SPOSet> upspo      = pftester.createDummySPO("up", nup, uporbs);
  std::unique_ptr<SPOSet> dnspo      = pftester.createDummySPO("dn", ndn, dnorbs);
  std::vector<std::unique_ptr<SPOSet>> spos;
  spos.emplace_back(std::move(upspo));
  spos.emplace_back(std::move(dnspo));
  PfaffianSTU pf1((*elec1), std::move(spos), "pfaffian1", opt, opt, opt);
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
  PfaffianSTU pf2((*elec2), std::move(spos2), "pfaffian2", opt, opt, opt);
  pftester.checkSizes(elec2->getTotalNum(), uporbs, dnorbs, pf2);
}

TEST_CASE("Pfaffian check pfaffian evaluation", "[wavefunction][fermion]")
{
  Communicate* comm = OHMMS::Controller;
  testing::PfaffianSTUTest pftester;

  std::unique_ptr<ParticleSet> elec = pftester.createDummyElec(3, 3);
  std::unique_ptr<SPOSet> upspo     = pftester.createDummySPO("up", 3, 3);
  std::unique_ptr<SPOSet> dnspo     = pftester.createDummySPO("dn", 3, 3);
  std::vector<std::unique_ptr<SPOSet>> spos;
  spos.emplace_back(std::move(upspo));
  spos.emplace_back(std::move(dnspo));

  const std::string opt = "no";
  PfaffianSTU pf((*elec), std::move(spos), "pfaffian", opt, opt, opt);
  pftester.checkEvaluation(pf, (*elec));
}

TEST_CASE("Pfaffian check VGL evaluations", "[wavefunction][fermion]")
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
  spos.push_back(std::move(upspo));
  spos.push_back(std::move(dnspo));

  const std::string opt = "no";
  PfaffianSTU pf((*elec), std::move(spos), opt, opt, opt);
  pftester.checkVGLEvaluations(pf, (*elec));
}

TEST_CASE("Pfaffian read/write vp", "[wavefunction][fermion]")
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
  spos.push_back(std::move(upspo));
  spos.push_back(std::move(dnspo));

  const std::string opt = "yes";
  PfaffianSTU pf((*elec), std::move(spos), opt, opt, opt);
  pftester.checkReadWrite(pf, (*elec));
}

TEST_CASE("Pfaffian evaluateDerivatives", "[wavefunction][fermion]")
{
  Communicate* comm = OHMMS::Controller;
  testing::PfaffianSTUTest pftester;

  int nup                           = 2;
  int ndn                           = 2;
  int norb                          = 3;
  std::unique_ptr<ParticleSet> elec = pftester.createDummyElec(nup, ndn);
  std::unique_ptr<SPOSet> upspo     = pftester.createDummySPO("up", nup, norb);
  std::unique_ptr<SPOSet> dnspo     = pftester.createDummySPO("dn", ndn, norb);
  std::vector<std::unique_ptr<SPOSet>> spos;
  spos.push_back(std::move(upspo));
  spos.push_back(std::move(dnspo));

  const std::string opt = "yes";
  PfaffianSTU pf((*elec), std::move(spos), opt, opt, opt);
  pftester.checkParmDerivs(pf, (*elec));
}

} // namespace qmcplusplus
