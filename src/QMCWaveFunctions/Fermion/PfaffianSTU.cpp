/////////////////////////////////////////////////////////////////////////////////////
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
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus
{

PfaffianSTU::PfaffianSTU(ParticleSet& targetPtcl,
                         std::vector<std::unique_ptr<SPOSet>>&& sposets,
                         const std::string& class_name)
    : active_idx_(-1),
      num_elec_(targetPtcl.getTotalNum()),
      num_up_(targetPtcl.last(0)),
      num_dn_(num_elec_ - num_up_),
      sposets_(std::move(sposets))
{
  resize();
}

PfaffianSTU::~PfaffianSTU() {}

bool PfaffianSTU::isOptimizable() const { return true; }

void PfaffianSTU::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) {}

void PfaffianSTU::checkOutVariables(const OptVariables& active) {}

PfaffianSTU::LogValue PfaffianSTU::evaluateLog(const ParticleSet& P,
                                               ParticleSet::ParticleGradient& G,
                                               ParticleSet::ParticleLaplacian& L)
{
  log_value_ = 0.0;
  recompute(P);
  ValueType val = calculatePfaffian();
  log_value_    = {std::log(std::abs(val)), std::arg(val)};
  calculateInverse();

  for (int ie = 0; ie < num_elec_; ie++)
  {
    mGradType rv   = simd::dot(psi_matinv_[ie], dpsi_mat_[ie], psi_mat_.rows());
    mValueType lap = simd::dot(psi_matinv_[ie], d2psi_mat_[ie], psi_mat_.rows());
    G[ie] += 0.5 * rv;
    L[ie] += 0.5 * (lap - dot(rv, rv));
  }

  return log_value_;
}

void PfaffianSTU::recompute(const ParticleSet& P)
{
  //update up
  sposets_[0]->evaluate_notranspose(P, 0, num_up_, up_psi_mat_, up_dpsi_mat_, up_d2psi_mat_);
  //update dn
  sposets_[1]->evaluate_notranspose(P, 0, num_dn_, dn_psi_mat_, dn_dpsi_mat_, dn_d2psi_mat_);

  const int norb = sposets_[0]->size();
  ValueVector tmpvec(norb);

  psi_mat_ = 0;
  //update upper diagonal of matrix
  for (int i = 0; i < num_elec_; i++)
  {
    for (int j = i + 1; j < num_elec_; j++)
    {
      //triplet uu
      if ((i < num_up_) && (j < num_up_))
      {
        for (int k = 0; k < norb; k++)
        {
          for (int l = 0; l < norb; l++)
          {
            psi_mat_(i, j) += up_psi_mat_(i, k) * uu_triplet_mat_(k, l) * up_psi_mat_(j, l);
            dpsi_mat_(i, j) += up_dpsi_mat_(i, k) * uu_triplet_mat_(k, l) * up_psi_mat_(j, l);
            d2psi_mat_(i, j) += up_d2psi_mat_(i, k) * uu_triplet_mat_(k, l) * up_psi_mat_(j, l);
          }
        }
      }
      // up dn singlet
      else if ((i < num_up_) && (j >= num_up_))
      {
        for (int k = 0; k < norb; k++)
        {
          for (int l = 0; l < norb; l++)
          {
            psi_mat_(i, j) += up_psi_mat_(i, k) * singlet_mat_(k, l) * dn_psi_mat_(j - num_up_, l);
            dpsi_mat_(i, j) += up_dpsi_mat_(i, k) * singlet_mat_(k, l) * dn_psi_mat_(j - num_up_, l);
            d2psi_mat_(i, j) += up_d2psi_mat_(i, k) * singlet_mat_(k, l) * dn_psi_mat_(j - num_up_, l);
          }
        }
      }
      //dn dn triplet
      else if ((i >= num_up_) && (j >= num_up_))
      {
        for (int k = 0; k < norb; k++)
        {
          for (int l = 0; l < norb; l++)
          {
            psi_mat_(i, j) += dn_psi_mat_(i - num_up_, k) * dd_triplet_mat_(k, l) * dn_psi_mat_(j - num_up_, l);
            dpsi_mat_(i, j) += dn_dpsi_mat_(i - num_up_, k) * dd_triplet_mat_(k, l) * dn_psi_mat_(j - num_up_, l);
            d2psi_mat_(i, j) += dn_d2psi_mat_(i - num_up_, k) * dd_triplet_mat_(k, l) * dn_psi_mat_(j - num_up_, l);
          }
        }
      }
      psi_mat_(j, i)   = -psi_mat_(i, j);
      dpsi_mat_(j, i)  = -dpsi_mat_(i, j);
      d2psi_mat_(j, i) = -d2psi_mat_(i, j);
    }
  }
  //Now need to update final col if odd num electrons
  if (psi_mat_.rows() == num_elec_ + 1)
  {
    for (int i = 0; i < num_up_; i++)
    {
      ValueType v              = up_psi_mat_(i, i);
      GradType g               = up_dpsi_mat_(i, i);
      ValueType l              = up_d2psi_mat_(i, i);
      psi_mat_(i, num_elec_)   = v;
      psi_mat_(num_elec_, i)   = -v;
      dpsi_mat_(i, num_elec_)  = g;
      dpsi_mat_(num_elec_, i)  = -g;
      d2psi_mat_(i, num_elec_) = l;
      d2psi_mat_(num_elec_, i) = -l;
    }
    for (int i = 0; i < num_dn_; i++)
    {
      ValueType v                        = dn_psi_mat_(i, i);
      GradType g                         = dn_dpsi_mat_(i, i);
      ValueType l                        = dn_d2psi_mat_(i, i);
      psi_mat_(num_up_ + i, num_elec_)   = v;
      psi_mat_(num_elec_, num_up_ + i)   = -v;
      dpsi_mat_(num_up_ + i, num_elec_)  = g;
      dpsi_mat_(num_elec_, num_up_ + i)  = -g;
      d2psi_mat_(num_up_ + i, num_elec_) = l;
      d2psi_mat_(num_elec_, num_up_ + i) = -l;
    }
  }
}

void PfaffianSTU::registerData(ParticleSet& P, WFBufferType& buf) {}

PfaffianSTU::LogValue PfaffianSTU::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch) {}

void PfaffianSTU::copyFromBuffer(ParticleSet& P, WFBufferType& buf) {}

PfaffianSTU::PsiValue PfaffianSTU::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) {}

PfaffianSTU::GradType PfaffianSTU::evalGrad(ParticleSet& P, int iat) {}

void PfaffianSTU::restore(int iat) {}

void PfaffianSTU::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  assert(iat == active_idx_);
  std::transform(psi_delta_.begin(), psi_delta_.end(), psi_mat_[active_idx_], psi_mat_[active_idx_],
                 [](auto v1, auto v2) { return v1 + v2; });
  std::transform(dpsi_delta_.begin(), dpsi_delta_.end(), dpsi_mat_[active_idx_], dpsi_mat_[active_idx_],
                 [](auto v1, auto v2) { return v1 + v2; });
  std::transform(d2psi_delta_.begin(), d2psi_delta_.end(), d2psi_mat_[active_idx_], d2psi_mat_[active_idx_],
                 [](auto v1, auto v2) { return v1 + v2; });
  for (int i = 0; i < psi_mat_.rows(); i++)
  {
    psi_mat_(i, active_idx_)   = -psi_mat_(active_idx_, i);
    dpsi_mat_(i, active_idx_)  = -dpsi_mat_(active_idx_, i);
    d2psi_mat_(i, active_idx_) = -d2psi_mat_(active_idx_, i);
  }
  updateInverse();
  active_idx_ = -1;
}

PfaffianSTU::PsiValue PfaffianSTU::ratio(ParticleSet& P, int iat) {}

std::unique_ptr<WaveFunctionComponent> PfaffianSTU::makeClone(ParticleSet& tqp) const {}

void PfaffianSTU::evaluateDerivatives(ParticleSet& P,
                                      const OptVariables& active,
                                      Vector<ValueType>& dlogpsi,
                                      Vector<ValueType>& dhpsioverpsi)
{}

void PfaffianSTU::evaluateDerivativesWF(ParticleSet& P, const OptVariables& active, Vector<ValueType>& dlogpsi) {}

void PfaffianSTU::resize()
{
  int rowsize = (num_elec_ % 2 == 0) ? num_elec_ : num_elec_ + 1;
  psi_mat_.resize(rowsize, rowsize);
  dpsi_mat_.resize(rowsize, rowsize);
  d2psi_mat_.resize(rowsize, rowsize);
  psi_matinv_.resize(rowsize, rowsize);
  psi_delta_.resize(rowsize);
  dpsi_delta_.resize(rowsize);
  d2psi_delta_.resize(rowsize);

  //now size the pairing function coefficient matrices matrices
  //up spos must be same size
  assert(sposets_[0]->size() == sposets_[1]->size());
  int norbs = sposets_[0]->size();
  singlet_mat_.resize(norbs, norbs);
  uu_triplet_mat_.resize(norbs, norbs);
  dd_triplet_mat_.resize(norbs, norbs);

  //storage for orbitals
  up_psi_mat_.resize(num_up_, norbs);
  dn_psi_mat_.resize(num_dn_, norbs);
  up_dpsi_mat_.resize(num_up_, norbs);
  dn_dpsi_mat_.resize(num_dn_, norbs);
  up_d2psi_mat_.resize(num_up_, norbs);
  dn_d2psi_mat_.resize(num_dn_, norbs);
}

int PfaffianSTU::rowPivot(ValueMatrix& mat, const int i)
{
  const int size = mat.rows();
  RealType tiny  = 1.0e-20;
  ValueType backup;
  int sign     = 1;
  RealType big = 0.0;

  int k = 0;
  for (int j = i + 1; j < size; j++)
  {
    RealType temp = std::abs(mat(i, j));
    if (temp > big)
    {
      big = temp;
      k   = j;
    }
  }
  if (big < tiny)
  {
    app_warning() << "Singular row in Pfaffian Matrix" << std::endl;
    mat(i, i + 1) = tiny;
  }
  if (k != (i + 1))
  {
    for (int j = i; j < size; j++)
    {
      backup        = mat(j, i + 1);
      mat(j, i + 1) = mat(j, k);
      mat(j, k)     = backup;
    }
    for (int j = i; j < size; j++)
    {
      backup        = mat(i + 1, j);
      mat(i + 1, j) = mat(k, j);
      mat(k, j)     = backup;
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

  ValueType pf = 1.0;
  int sign     = 1;
  for (int i = 0; i < size; i += 2)
  {
    sign *= rowPivot(tmp_mat, i);
    for (int j = i + 2; j < size; j++)
    {
      ValueType fac = -tmp_mat(i, j) / tmp_mat(i, i + 1);
      for (int k = i + 1; k < size; k++)
      {
        tmp_mat(k, j) += fac * tmp_mat(k, i + 1);
        tmp_mat(j, k) += fac * tmp_mat(i + 1, k);
      }
    }
    pf *= tmp_mat(i, i + 1);
  }
  return pf * ValueType(sign);
}

void PfaffianSTU::calculateInverse()
{
  std::copy(psi_mat_.begin(), psi_mat_.end(), psi_matinv_.begin());
  invert_matrix(psi_matinv_, false);
}

PfaffianSTU::ValueType PfaffianSTU::calculateRatio(const ValueVector& newvals)
{
  assert(active_idx_ >= 0);
  assert(newvals.size() == psi_mat_.rows());
  std::transform(newvals.begin(), newvals.end(), psi_mat_[active_idx_], psi_delta_.begin(),
                 [](auto v1, auto v2) { return v1 - v2; });

  //can't use simd::dot since dotting into column of inverse matrix...not contiguous
  ValueType ratio = 0.0;
  for (int i = 0; i < psi_mat_.rows(); i++)
    ratio += psi_delta_[i] * psi_matinv_(i, active_idx_);
  return 1.0 + ratio;
}

void PfaffianSTU::updateInverse()
{
  const int n = psi_matinv_.rows();
  ValueVector u(n, 0.0);
  u[active_idx_] = 1.0;

  ValueMatrix U(n, 2);
  ValueMatrix V(2, n);

  for (int i = 0; i < n; i++)
  {
    U(i, 0) = u[i];
    U(i, 1) = psi_delta_[i];
    V(0, i) = psi_delta_[i];
    V(1, i) = -u[i];
  }

  ValueMatrix tmp(n, 2);
  MatrixOperators::product(psi_matinv_, U, tmp);

  ValueMatrix M(2, 2);
  MatrixOperators::product(V, tmp, M);

  for (int i = 0; i < 2; i++)
    M(i, i) += 1.0;

  invert_matrix(M, false);

  ValueMatrix tmp2(2, n);
  MatrixOperators::product(V, psi_matinv_, tmp2);
  MatrixOperators::product(U, M, tmp);

  ValueMatrix tmp3(n, n);
  ValueMatrix tmp4(n, n);
  MatrixOperators::product(tmp, tmp2, tmp3);
  MatrixOperators::product(psi_matinv_, tmp3, tmp4);

  std::transform(psi_matinv_.begin(), psi_matinv_.end(), tmp4.begin(), psi_matinv_.begin(),
                 [](auto v1, auto v2) { return v1 - v2; });
}

} // namespace qmcplusplus
