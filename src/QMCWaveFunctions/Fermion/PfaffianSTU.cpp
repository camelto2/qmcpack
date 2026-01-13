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
    : RatioTimer(createGlobalTimer(class_name + "::ratio", timer_level_fine)),
      SPOVTimer(createGlobalTimer(class_name + "::spoval", timer_level_fine)),
      active_idx_(-1),
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

  const int size = psi_mat_.rows();
  for (int i = 0; i < num_elec_; i++)
  {
    //exploting symmetry of matrices here. psi_matint_[i] gives a row, but I need to dot with column.
    //Since antisymmetric, add a sign to G and L
    mGradType rv   = simd::dot(psi_matinv_[i], dpsi_rows_[i], size);
    mValueType lap = simd::dot(psi_matinv_[i], d2psi_rows_[i], size);
    G[i] -= rv;
    L[i] -= (lap + dot(rv, rv));
  }

  return log_value_;
}

void PfaffianSTU::recompute(const ParticleSet& P)
{
  //update up
  up_psi_mat_   = 0;
  up_dpsi_mat_  = 0;
  up_d2psi_mat_ = 0;
  sposets_[0]->evaluate_notranspose(P, 0, num_up_, up_psi_mat_, up_dpsi_mat_, up_d2psi_mat_);
  //update dn
  dn_psi_mat_   = 0;
  dn_dpsi_mat_  = 0;
  dn_d2psi_mat_ = 0;
  sposets_[1]->evaluate_notranspose(P, num_up_, num_elec_, dn_psi_mat_, dn_dpsi_mat_, dn_d2psi_mat_);

  const int norb = sposets_[0]->size();

  psi_mat_    = 0;
  dpsi_rows_  = 0;
  d2psi_rows_ = 0;
  //update upper diagonal of matrix
  for (int i = 0; i < num_elec_; i++)
  {
    bool iup = (i < num_up_);
    int ii   = iup ? i : i - num_up_;
    for (int j = i + 1; j < num_elec_; j++)
    {
      bool jup = (j < num_up_);
      int jj   = jup ? j : j - num_up_;
      //now get references and pointers to orbitals and pairing matrix depending on type
      auto& pair_mat = (iup == jup) ? (iup ? uu_triplet_mat_ : dd_triplet_mat_) : singlet_mat_;
      auto* psi_i    = iup ? up_psi_mat_[ii] : dn_psi_mat_[ii];
      auto* psi_j    = jup ? up_psi_mat_[jj] : dn_psi_mat_[jj];
      auto* dpsi_i   = iup ? up_dpsi_mat_[ii] : dn_dpsi_mat_[ii];
      auto* dpsi_j   = jup ? up_dpsi_mat_[jj] : dn_dpsi_mat_[jj];
      auto* d2psi_i  = iup ? up_d2psi_mat_[ii] : dn_d2psi_mat_[ii];
      auto* d2psi_j  = jup ? up_d2psi_mat_[jj] : dn_d2psi_mat_[jj];
      for (int k = 0; k < norb; k++)
      {
        for (int l = 0; l < norb; l++)
        {
          psi_mat_(i, j) += psi_i[k] * pair_mat(k, l) * psi_j[l];
          dpsi_rows_(i, j) += dpsi_i[k] * pair_mat(k, l) * psi_j[l];
          dpsi_rows_(j, i) -= psi_i[k] * pair_mat(k, l) * dpsi_j[l];
          d2psi_rows_(i, j) += d2psi_i[k] * pair_mat(k, l) * psi_j[l];
          d2psi_rows_(j, i) -= psi_i[k] * pair_mat(k, l) * d2psi_j[l];
        }
      }
      psi_mat_(j, i) = -psi_mat_(i, j);
    }
  }
  //Now need to update final col if odd num electrons
  if (psi_mat_.rows() == num_elec_ + 1)
  {
    for (int i = 0; i < num_elec_; i++)
    {
      bool iup                  = (i < num_up_);
      int ii                    = iup ? i : i - num_up_;
      ValueType v               = iup ? up_psi_mat_(ii, ii) : dn_psi_mat_(ii, ii);
      GradType g                = iup ? up_dpsi_mat_(ii, ii) : dn_dpsi_mat_(ii, ii);
      ValueType l               = iup ? up_d2psi_mat_(ii, ii) : dn_d2psi_mat_(ii, ii);
      psi_mat_(i, num_elec_)    = v;
      psi_mat_(num_elec_, i)    = -v;
      dpsi_rows_(i, num_elec_)  = g;
      dpsi_rows_(num_elec_, i)  = -g;
      d2psi_rows_(i, num_elec_) = l;
      d2psi_rows_(num_elec_, i) = -l;
    }
  }
}

void PfaffianSTU::registerData(ParticleSet& P, WFBufferType& buf) {}

//for now just call evaluateLog
PfaffianSTU::LogValue PfaffianSTU::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch) 
{
  ParticleSet::ParticleGradient G(num_elec_);
  ParticleSet::ParticleLaplacian L(num_elec_);
  return evaluateLog(P, G, L);
}

void PfaffianSTU::copyFromBuffer(ParticleSet& P, WFBufferType& buf) {}

PfaffianSTU::PsiValue PfaffianSTU::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) {}

PfaffianSTU::GradType PfaffianSTU::evalGrad(ParticleSet& P, int iat)
{
  ScopedTimer local_timer(RatioTimer);

  const int size = psi_mat_.rows();
  assert((iat >= 0) && (iat < num_elec_));
  //exploting symmetry of matrices here. psi_matint_[i] gives a row, but I need to dot with column.
  //Since antisymmetric, add a sign
  GradType grad = -simd::dot(psi_matinv_[iat], dpsi_rows_[iat], size);
  return grad;
}

void PfaffianSTU::restore(int iat) {}

void PfaffianSTU::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  assert(iat == active_idx_);
  std::transform(psi_delta_.begin(), psi_delta_.end(), psi_mat_[active_idx_], psi_mat_[active_idx_],
                 [](auto v1, auto v2) { return v1 + v2; });
  for (int i = 0; i < psi_mat_.rows(); i++)
  {
    psi_mat_(i, active_idx_) = -psi_mat_(active_idx_, i);
  }
  updateInverse();
  active_idx_ = -1;
}

PfaffianSTU::PsiValue PfaffianSTU::ratio(ParticleSet& P, int iat)
{
  active_idx_ = iat;
  {
    ScopedTimer local_timer(SPOVTimer);
    const int group = P.getGroupID(iat);
    sposets_[group]->evaluateValue(P, iat, tmp_psi_);
  }

  const int norb = sposets_[0]->size();
  ValueVector row_update(psi_mat_.rows());
  bool iup = (iat < num_up_);
  for (int j = 0; j < num_elec_; j++)
  {
    if (j == iat)
      continue;
    bool jup = (j < num_up_);
    int jj   = jup ? j : j - num_up_;

    //for a row update, i need the full pfaffian matrix to be antisymmetric
    //for singlets, the pairing matrix is symmetric
    //for triplets, the pairing matrix is antisymmetric
    //if the pairing matrix is symmetric, I need to add an explicit sign for the lower diagonal
    //but if pairing is antisymmetric, there is no need since pairing(i,j) = -pairing(j,i)
    //look at https://arxiv.org/pdf/1008.2369 and equation 151 for an example
    //therefore, add sign only for singlet and if below diagonal
    ValueType sign = ((j < iat) && (iup != jup)) ? -1.0 : 1.0;

    auto& pair_mat = (iup == jup) ? (iup ? uu_triplet_mat_ : dd_triplet_mat_) : singlet_mat_;
    auto* psi_j    = jup ? up_psi_mat_[jj] : dn_psi_mat_[jj];

    for (int k = 0; k < norb; k++)
    {
      for (int l = 0; l < norb; l++)
      {
        row_update[j] += sign * tmp_psi_[k] * pair_mat(k, l) * psi_j[l];
      }
    }
  }

  if (psi_mat_.rows() == num_elec_ + 1)
  {
    bool iup              = (iat < num_up_);
    int ii                = iup ? iat : iat - num_up_;
    row_update[num_elec_] = tmp_psi_[ii];
  }

  return calculateRatio(row_update);
}

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
  psi_matinv_.resize(rowsize, rowsize);
  psi_delta_.resize(rowsize);
  dpsi_rows_.resize(rowsize, rowsize);
  d2psi_rows_.resize(rowsize, rowsize);

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

  tmp_psi_.resize(norbs);
  tmp_dpsi_.resize(norbs);
  tmp_d2psi_.resize(norbs);
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
