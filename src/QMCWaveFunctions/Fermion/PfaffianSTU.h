//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_PFAFFIANSTU
#define QMCPLUSPLUS_PFAFFIANSTU

#include "QMCWaveFunctions/SPOSet.h"
#include "Utilities/TimerManager.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"

namespace qmcplusplus
{

namespace testing
{
class PfaffianSTUTest;
}

class PfaffianSTU : public WaveFunctionComponent, public OptimizableObject
{
  using ValueVector = SPOSet::ValueVector;
  using GradVector  = SPOSet::GradVector;
  using ValueMatrix = SPOSet::ValueMatrix;
  using GradMatrix  = SPOSet::GradMatrix;

  using mValueType = QMCTraits::QTFull::ValueType;
  using mGradType  = TinyVector<mValueType, DIM>;

public:
  PfaffianSTU(ParticleSet& targetPtcl,
              std::vector<std::unique_ptr<SPOSet>>&& sposets,
              const std::string& opt_singlet, 
              const std::string& opt_uu_triplet, 
              const std::string& opt_dd_triplet,
              const std::string& class_name = "PfaffianSTU");

  ///destructor
  ~PfaffianSTU() override;

  std::string getClassName() const override { return "PfaffianSTU"; }

  bool isFermionic() const final { return true; }
  bool isOptimizable() const override;

  void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override;

  void checkOutVariables(const OptVariables& active) override;

  void checkInVariablesExclusive(OptVariables& active) override;

  void resetParametersExclusive(const OptVariables& active) override; 

  LogValue evaluateLog(const ParticleSet& P,
                       ParticleSet::ParticleGradient& G,
                       ParticleSet::ParticleLaplacian& L) override;

  void updateAfterSweep(const ParticleSet& P,
                        ParticleSet::ParticleGradient& G,
                        ParticleSet::ParticleLaplacian& L);

  void recompute(const ParticleSet& P) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValue updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  PsiValue ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  GradType evalGrad(ParticleSet& P, int iat) override;

  void restore(int iat) override;

  void completeUpdates() override;

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  PsiValue ratio(ParticleSet& P, int iat) override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  void evaluateDerivatives(ParticleSet& P,
                           const OptVariables& active,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi) override;

  void evaluateDerivativesWF(ParticleSet& P, const OptVariables& active, Vector<ValueType>& dlogpsi) override;

  void buildOptVariables();

protected:
  //Timers 
  NewTimer &UpdateTimer, &RatioTimer, &InverseTimer, &BufferTimer, &SPOVTimer, &SPOVGLTimer; 

private:
  void resize();

  //Implements code from M. Bajdich thesis to do row pivot used in calculatePfaffian.
  int rowPivot(ValueMatrix& mat, const int i);

  //Implements code from M. Bajdich thesis to calculate pfaffian from psi_mat_
  ValueType calculatePfaffian();

  void calculateInverse();

  ValueType calculateRatio(const ValueVector& newvals);

  //called by acceptMove, updates the inverse with sherman-morrison-woodbury update
  void updateInverse();

  void initializePairingMats();

  void updatePairingNorms();

  //current wavefunction ratio
  PsiValue cur_ratio_;
  
  //whether to optimize the pairing matrices, allows to enable AGPs with only singlet on
  bool opt_singlet_;
  bool opt_uu_triplet_;
  bool opt_dd_triplet_;

  //current matrix 
  ValueMatrix psi_mat_;

  //row updates only for derivatives. 
  //We don't need the full matrix, only a collection of row updates for each particle
  //We do this since we only ever do dot products of rows with columns of inverse matrix
  GradMatrix dpsi_rows_;
  ValueMatrix d2psi_rows_;

  //current inverse
  ValueMatrix psi_matinv_;

  //change to row/column
  ValueVector psi_delta_;
  //new gradient, lap
  GradVector dpsi_new_;
  ValueVector d2psi_new_;

  //active row/column
  int active_idx_;

  int num_elec_;
  int num_up_;
  int num_dn_;
  int size_;

  //Pairing function coeffs
  ValueMatrix singlet_mat_;
  ValueMatrix uu_triplet_mat_;
  ValueMatrix dd_triplet_mat_;
  RealType singlet_inv_norm_;
  RealType uu_triplet_inv_norm_;
  RealType dd_triplet_inv_norm_;

  //Orbital values, grads, laps
  ValueMatrix up_psi_mat_;
  ValueMatrix dn_psi_mat_;
  GradMatrix  up_dpsi_mat_;
  GradMatrix  dn_dpsi_mat_;
  ValueMatrix up_d2psi_mat_;
  ValueMatrix dn_d2psi_mat_;

  //VGL for particle update
  ValueVector tmp_psi_;
  GradVector  tmp_dpsi_;
  ValueVector tmp_d2psi_;

  std::vector<std::unique_ptr<SPOSet>> sposets_;

  friend class testing::PfaffianSTUTest;
};

} // namespace qmcplusplus

#endif
