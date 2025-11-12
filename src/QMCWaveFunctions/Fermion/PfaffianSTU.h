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

#include "QMCWaveFunctions/WaveFunctionComponent.h"
namespace qmcplusplus
{

namespace testing
{
  class PfaffianSTUTest;
}

class PfaffianSTU : public WaveFunctionComponent
{
  using ValueVector = SPOSet::ValueVector;
  using ValueMatrix = SPOSet::ValueMatrix;

public:
  PfaffianSTU(ParticleSet& targetPtcl, const std::string& class_name = "PfaffianSTU");

  ///destructor
  ~PfaffianSTU() override;

  std::string getClassName() const override { return "PfaffianSTU"; }

  bool isFermionic() const final { return true; }
  bool isOptimizable() const override;

  void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override;

  void checkOutVariables(const opt_variables_type& active) override;

  LogValue evaluateLog(const ParticleSet& P,
                       ParticleSet::ParticleGradient& G,
                       ParticleSet::ParticleLaplacian& L) override;

  void recompute(const ParticleSet& P) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValue updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  PsiValue ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;


  GradType evalGrad(ParticleSet& P, int iat) override;


  void restore(int iat) override;

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  PsiValue ratio(ParticleSet& P, int iat) override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  void evaluateDerivatives(ParticleSet& P, const opt_variables_type& active, Vector<ValueType>& dlogpsi, Vector<ValueType>& dhpsioverpsi) override;
  void evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& active, Vector<ValueType>& dlogpsi) override;

private:

  void resize();

  //Implements code from M. Bajdich thesis to do row pivot used in calculatePfaffian. 
  int rowPivot(ValueMatrix& mat, const int i);

  //Implements code from M. Bajdich thesis to calculate pfaffian from psi_mat_
  ValueType calculatePfaffian();

  //current matrix
  ValueMatrix psi_mat_;

  //current inverse
  ValueMatrix psi_matinv_;

  //values to update
  ValueVector psi_val_;

  //active row/column
  int active_idx_;
  
  int num_elec_;

  friend class testing::PfaffianSTUTest;
};

} // namespace qmcplusplus

#endif
