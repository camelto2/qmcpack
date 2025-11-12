//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "PfaffianSTU.h"

namespace qmcplusplus
{

PfaffianSTU::PfaffianSTU(ParticleSet& targetPtcl, const std::string& class_name) : active_idx_(-1), num_elec_(targetPtcl.getTotalNum())
{
  resize();
}

PfaffianSTU::~PfaffianSTU()  {}

bool PfaffianSTU::isOptimizable() const  { return true; }

void PfaffianSTU::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs)  {}

void PfaffianSTU::checkOutVariables(const opt_variables_type& active)  {}

PfaffianSTU::LogValue PfaffianSTU::evaluateLog(const ParticleSet& P,
                                               ParticleSet::ParticleGradient& G,
                                               ParticleSet::ParticleLaplacian& L) 
{}

void PfaffianSTU::recompute(const ParticleSet& P) 
{}

void PfaffianSTU::registerData(ParticleSet& P, WFBufferType& buf) 
{}

PfaffianSTU::LogValue PfaffianSTU::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch) 
{}

void PfaffianSTU::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{}

PfaffianSTU::PsiValue PfaffianSTU::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) 
{}

PfaffianSTU::GradType PfaffianSTU::evalGrad(ParticleSet& P, int iat) 
{}

void PfaffianSTU::restore(int iat) 
{}

void PfaffianSTU::acceptMove(ParticleSet& P, int iat, bool safe_to_delay) 
{}

PfaffianSTU::PsiValue PfaffianSTU::ratio(ParticleSet& P, int iat) 
{}

std::unique_ptr<WaveFunctionComponent> PfaffianSTU::makeClone(ParticleSet& tqp) const 
{}

void PfaffianSTU::evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& active, Vector<ValueType>& dlogpsi) 
{}

void PfaffianSTU::resize() 
{
  int rowsize = ( num_elec_%2 == 0 ) ? num_elec_ : num_elec_ + 1;
  psi_mat_.resize(rowsize, rowsize);
  psi_matinv_.resize(rowsize, rowsize);
  psi_val_.resize(rowsize);
}

} // namespace qmcplusplus
