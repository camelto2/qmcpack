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

PfaffianSTU::PfaffianSTU(ParticleSet& targetPtcl, const std::string& class_name)
    : active_idx_(-1), num_elec_(targetPtcl.getTotalNum())
{
  resize();
}

PfaffianSTU::~PfaffianSTU() {}

bool PfaffianSTU::isOptimizable() const { return true; }

void PfaffianSTU::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) {}

void PfaffianSTU::checkOutVariables(const opt_variables_type& active) {}

PfaffianSTU::LogValue PfaffianSTU::evaluateLog(const ParticleSet& P,
                                               ParticleSet::ParticleGradient& G,
                                               ParticleSet::ParticleLaplacian& L)
{}

void PfaffianSTU::recompute(const ParticleSet& P) {}

void PfaffianSTU::registerData(ParticleSet& P, WFBufferType& buf) {}

PfaffianSTU::LogValue PfaffianSTU::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch) {}

void PfaffianSTU::copyFromBuffer(ParticleSet& P, WFBufferType& buf) {}

PfaffianSTU::PsiValue PfaffianSTU::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) {}

PfaffianSTU::GradType PfaffianSTU::evalGrad(ParticleSet& P, int iat) {}

void PfaffianSTU::restore(int iat) {}

void PfaffianSTU::acceptMove(ParticleSet& P, int iat, bool safe_to_delay) {}

PfaffianSTU::PsiValue PfaffianSTU::ratio(ParticleSet& P, int iat) {}

std::unique_ptr<WaveFunctionComponent> PfaffianSTU::makeClone(ParticleSet& tqp) const {}

void PfaffianSTU::evaluateDerivatives(ParticleSet& P,
                                      const opt_variables_type& active,
                                      Vector<ValueType>& dlogpsi,
                                      Vector<ValueType>& dhpsioverpsi)
{}

void PfaffianSTU::evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& active, Vector<ValueType>& dlogpsi) {}

void PfaffianSTU::resize()
{
  int rowsize = (num_elec_ % 2 == 0) ? num_elec_ : num_elec_ + 1;
  psi_mat_.resize(rowsize, rowsize);
  psi_matinv_.resize(rowsize, rowsize);
  psi_val_.resize(rowsize);
}

int PfaffianSTU::rowPivot(ValueMatrix& mat, const int i) 
{
  const int size = mat.rows();
  RealType tiny = 1.0e-20;
  ValueType backup;
  int sign = 1;
  RealType big = 0.0;

  int k = 0;
  for (int j = i + 1; j < size; j++)
  {
    RealType temp = std::abs(mat(i,j));
    if (temp > big)
    {
      big = temp;
      k = j;
    }
  }
  if (big < tiny) 
  {
    app_warning() << "Singular row in Pfaffian Matrix" << std::endl;
    mat(i, i+1) = tiny;
  }
  if (k != (i + 1))
  {
    for (int j = i; j < size; j++)
    {
      backup = mat(j, i+1);
      mat(j, i+1) = mat(j, k);
      mat(j,k) = backup;
    }
    for (int j = i; j < size; j++)
    {
      backup = mat(i+1, j);
      mat(i+1,j) = mat(k, j);
      mat(k, j) = backup;
    }
    sign *= -1;
  }
  return sign;
}

PfaffianSTU::ValueType PfaffianSTU::calculatePfaffian()
{
  const int size = psi_mat_.rows();
  ValueMatrix tmp_mat(size, size);
  std::copy(psi_mat_.begin(), psi_mat_.end(), tmp_mat.begin());

  ValueType pf  = 1.0;
  int sign = 1;
  for (int i = 0; i < size; i += 2)
  {
    sign *= rowPivot(tmp_mat, i);
    for (int j = i + 2; j < size; j++)
    {
      ValueType fac = -tmp_mat(i,j) / tmp_mat(i, i+1);
      for (int k = i + 1; k < size; k++)
      {
        tmp_mat(k,j) += fac * tmp_mat(k, i+1);
        tmp_mat(j,k) += fac * tmp_mat(i+1, k);
      }
    }
    pf *= tmp_mat(i, i+1);
  }
  return pf * ValueType(sign);
}

} // namespace qmcplusplus
