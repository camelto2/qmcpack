/////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "PfaffianSTU.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"
#include "io/hdf/hdf_archive.h"

namespace qmcplusplus
{

PfaffianSTU::PfaffianSTU(ParticleSet& targetPtcl,
                         std::vector<std::unique_ptr<SPOSet>>&& sposets,
                         const std::string& opt_singlet,
                         const std::string& opt_uu_triplet,
                         const std::string& opt_dd_triplet,
                         const std::string& class_name)
    : WaveFunctionComponent(class_name),
      OptimizableObject(class_name),
      UpdateTimer(createGlobalTimer(class_name + "::update", timer_level_fine)),
      RatioTimer(createGlobalTimer(class_name + "::ratio", timer_level_fine)),
      InverseTimer(createGlobalTimer(class_name + "::inverse", timer_level_fine)),
      BufferTimer(createGlobalTimer(class_name + "::buffer", timer_level_fine)),
      SPOVTimer(createGlobalTimer(class_name + "::spoval", timer_level_fine)),
      SPOVGLTimer(createGlobalTimer(class_name + "::spovgl", timer_level_fine)),
      opt_singlet_(opt_singlet == "yes"),
      opt_uu_triplet_(opt_uu_triplet == "yes"),
      opt_dd_triplet_(opt_dd_triplet == "yes"),
      active_idx_(-1),
      num_elec_(targetPtcl.getTotalNum()),
      num_up_(targetPtcl.last(0)),
      num_dn_(num_elec_ - num_up_),
      size_((num_elec_ % 2) == 0 ? num_elec_ : num_elec_ + 1),
      sposets_(std::move(sposets))
{
  resize();
  initializePairingMats();
}

PfaffianSTU::~PfaffianSTU() {}

bool PfaffianSTU::isOptimizable() const { return true; }

void PfaffianSTU::initializePairingMats()
{
  singlet_mat_    = 0.0;
  uu_triplet_mat_ = 0.0;
  dd_triplet_mat_ = 0.0;
  if (opt_singlet_)
  {
    app_log() << "  PfaffianSTU: Initializing singlet" << std::endl;
    const int min = num_up_ >= num_dn_ ? num_dn_ : num_up_;
    for (int i = 0; i < min; i++)
      singlet_mat_(i, i) = 1.0;
  }
  if (opt_uu_triplet_)
  {
    app_log() << "  PfaffianSTU: Initializing uu triplet" << std::endl;
    for (int i = 0; i < num_up_ - 1; i += 2)
    {
      uu_triplet_mat_(i, i + 1) = 1.0;
      uu_triplet_mat_(i + 1, i) = -1.0;
    }
  }
  if (opt_dd_triplet_)
  {
    app_log() << "  PfaffianSTU: Initializing dd triplet" << std::endl;
    for (int i = 0; i < num_dn_ - 1; i += 2)
    {
      dd_triplet_mat_(i, i + 1) = 1.0;
      dd_triplet_mat_(i + 1, i) = -1.0;
    }
  }
}

PfaffianSTU::LogValue PfaffianSTU::evaluateLog(const ParticleSet& P,
                                               ParticleSet::ParticleGradient& G,
                                               ParticleSet::ParticleLaplacian& L)
{
  log_value_ = 0.0;
  recompute(P);
  ValueType val = calculatePfaffian();
  log_value_    = {std::log(std::abs(val)), std::arg(val)};
  calculateInverse();

  for (int i = 0; i < num_elec_; i++)
  {
    //exploting symmetry of matrices here. psi_matinv_[i] gives a row, but I need to dot with column.
    //Since antisymmetric, add a sign to G and L
    mGradType rv   = -simd::dot(psi_matinv_[i], dpsi_rows_[i], size_);
    mValueType lap = -simd::dot(psi_matinv_[i], d2psi_rows_[i], size_);
    G[i] += rv;
    L[i] += (lap - dot(rv, rv));
  }

  return log_value_;
}

void PfaffianSTU::updateAfterSweep(const ParticleSet& P,
                                   ParticleSet::ParticleGradient& G,
                                   ParticleSet::ParticleLaplacian& L)
{
  if (UpdateMode == ORB_PBYP_RATIO)
    recompute(P);

  //no update to the inverse matrix

  for (int i = 0; i < num_elec_; i++)
  {
    //exploting symmetry of matrices here. psi_matinv_[i] gives a row, but I need to dot with column.
    //Since antisymmetric, add a sign to G and L
    mGradType rv   = -simd::dot(psi_matinv_[i], dpsi_rows_[i], size_);
    mValueType lap = -simd::dot(psi_matinv_[i], d2psi_rows_[i], size_);
    G[i] += rv;
    L[i] += (lap - dot(rv, rv));
  }
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
  if (size_ == num_elec_ + 1)
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

  UpdateMode = ORB_WALKER;
}


void PfaffianSTU::registerData(ParticleSet& P, WFBufferType& buf)
{
  const int norbs = sposets_[0]->size();
  if (Bytes_in_WFBuffer == 0)
  {
    Bytes_in_WFBuffer = buf.current();
    buf.add(psi_mat_.first_address(), psi_mat_.last_address());
    buf.add(psi_matinv_.first_address(), psi_matinv_.last_address());
    buf.add(&(dpsi_rows_(0, 0)[0]), &(dpsi_rows_(0, 0)[0]) + size_ * size_ * DIM);
    buf.add(d2psi_rows_.first_address(), d2psi_rows_.last_address());
    buf.add(up_psi_mat_.first_address(), up_psi_mat_.last_address());
    buf.add(&(up_dpsi_mat_(0, 0)[0]), &(up_dpsi_mat_(0, 0)[0]) + num_up_ * norbs * DIM);
    buf.add(up_d2psi_mat_.first_address(), up_d2psi_mat_.last_address());
    buf.add(dn_psi_mat_.first_address(), dn_psi_mat_.last_address());
    buf.add(&(dn_dpsi_mat_(0, 0)[0]), &(dn_dpsi_mat_(0, 0)[0]) + num_dn_ * norbs * DIM);
    buf.add(dn_d2psi_mat_.first_address(), dn_d2psi_mat_.last_address());
    Bytes_in_WFBuffer = buf.current() - Bytes_in_WFBuffer;
    psi_mat_.free();
    psi_matinv_.free();
    dpsi_rows_.free();
    d2psi_rows_.free();
    up_psi_mat_.free();
    up_dpsi_mat_.free();
    up_d2psi_mat_.free();
    dn_psi_mat_.free();
    dn_dpsi_mat_.free();
    dn_d2psi_mat_.free();
  }
  else
  {
    buf.forward(Bytes_in_WFBuffer);
  }
  buf.add(log_value_);
}

//for now just call evaluateLog
PfaffianSTU::LogValue PfaffianSTU::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  if (fromscratch)
    evaluateLog(P, P.G, P.L);
  else
    updateAfterSweep(P, P.G, P.L);
  {
    ScopedTimer local_timer(BufferTimer);
    buf.forward(Bytes_in_WFBuffer);
    buf.put(log_value_);
  }
  return log_value_;
}

void PfaffianSTU::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  const int norbs = sposets_[0]->size();
  ScopedTimer local_timer(BufferTimer);
  psi_mat_.attachReference(buf.lendReference<ValueType>(size_ * size_), size_, size_);
  psi_matinv_.attachReference(buf.lendReference<ValueType>(size_ * size_), size_, size_);
  dpsi_rows_.attachReference(buf.lendReference<GradType>(size_ * size_), size_, size_);
  d2psi_rows_.attachReference(buf.lendReference<ValueType>(size_ * size_), size_, size_);
  up_psi_mat_.attachReference(buf.lendReference<ValueType>(num_up_ * norbs), num_up_, norbs);
  up_dpsi_mat_.attachReference(buf.lendReference<GradType>(num_up_ * norbs), num_up_, norbs);
  up_d2psi_mat_.attachReference(buf.lendReference<ValueType>(num_up_ * norbs), num_up_, norbs);
  dn_psi_mat_.attachReference(buf.lendReference<ValueType>(num_dn_ * norbs), num_dn_, norbs);
  dn_dpsi_mat_.attachReference(buf.lendReference<GradType>(num_dn_ * norbs), num_dn_, norbs);
  dn_d2psi_mat_.attachReference(buf.lendReference<ValueType>(num_dn_ * norbs), num_dn_, norbs);
  buf.get(log_value_);
  active_idx_ = -1;
}

PfaffianSTU::PsiValue PfaffianSTU::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  active_idx_ = iat;
  {
    ScopedTimer local_timer(SPOVGLTimer);
    const int group = P.getGroupID(iat);
    sposets_[group]->evaluateVGL(P, iat, tmp_psi_, tmp_dpsi_, tmp_d2psi_);
  }

  ScopedTimer local_timer(RatioTimer);
  UpdateMode     = ORB_PBYP_PARTIAL;
  const int norb = sposets_[0]->size();
  ValueVector row_update(size_, 0.0);
  dpsi_new_  = 0;
  d2psi_new_ = 0;
  bool iup   = (iat < num_up_);
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
        dpsi_new_[j] += sign * tmp_dpsi_[k] * pair_mat(k, l) * psi_j[l];
        d2psi_new_[j] += sign * tmp_d2psi_[k] * pair_mat(k, l) * psi_j[l];
      }
    }
  }

  if (size_ == num_elec_ + 1)
  {
    bool iup              = (iat < num_up_);
    int ii                = iup ? iat : iat - num_up_;
    row_update[num_elec_] = tmp_psi_[ii];
    dpsi_new_[num_elec_]  = tmp_dpsi_[ii];
    d2psi_new_[num_elec_] = tmp_d2psi_[ii];
  }

  cur_ratio_ = calculateRatio(row_update);
  //exploting symmetry of matrices here. psi_matint_[i] gives a row, but I need to dot with column.
  //Since antisymmetric, add a sign
  grad_iat -= simd::dot(psi_matinv_[iat], dpsi_new_.data(), size_) / cur_ratio_;

  return cur_ratio_;
}

PfaffianSTU::GradType PfaffianSTU::evalGrad(ParticleSet& P, int iat)
{
  ScopedTimer local_timer(RatioTimer);

  assert((iat >= 0) && (iat < num_elec_));
  //exploting symmetry of matrices here. psi_matint_[i] gives a row, but I need to dot with column.
  //Since antisymmetric, add a sign
  GradType grad = -simd::dot(psi_matinv_[iat], dpsi_rows_[iat], size_);
  return grad;
}

void PfaffianSTU::restore(int iat) { cur_ratio_ = 1.0; }

void PfaffianSTU::completeUpdates() { active_idx_ = -1; }

void PfaffianSTU::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  if (cur_ratio_ == PsiValue(0))
  {
    std::ostringstream msg;
    msg << "PfaffianSTU::acceptMove cur_ratio_ is " << cur_ratio_ << "! Report a bug." << std::endl;
    throw std::runtime_error(msg.str());
  }
  ScopedTimer local_timer(UpdateTimer);
  assert(iat == active_idx_);

  log_value_ += convertValueToLog(cur_ratio_);

  std::transform(psi_delta_.begin(), psi_delta_.end(), psi_mat_[active_idx_], psi_mat_[active_idx_],
                 [](auto v1, auto v2) { return v1 + v2; });
  for (int i = 0; i < psi_mat_.rows(); i++)
    psi_mat_(i, active_idx_) = -psi_mat_(active_idx_, i);
  const int ii   = active_idx_ < num_up_ ? active_idx_ : active_idx_ - num_up_;
  const bool iup = (iat < num_up_);
  simd::copy(iup ? up_psi_mat_[ii] : dn_psi_mat_[ii], tmp_psi_.data(), tmp_psi_.size());
  updateInverse();
  if (UpdateMode == ORB_PBYP_PARTIAL)
  {
    simd::copy(dpsi_rows_[active_idx_], dpsi_new_.data(), size_);
    simd::copy(d2psi_rows_[active_idx_], d2psi_new_.data(), size_);

    //active_idx_ column of dpsi_rows_ and d2psi_rows is also needs update
    const int norb = sposets_[0]->size();
    simd::copy(iup ? up_dpsi_mat_[ii] : dn_dpsi_mat_[ii], tmp_dpsi_.data(), norb);
    simd::copy(iup ? up_d2psi_mat_[ii] : dn_d2psi_mat_[ii], tmp_d2psi_.data(), norb);

    for (int j = 0; j < num_elec_; j++)
    {
      if (j == iat)
        continue;

      const bool jup = (j < num_up_);
      const int jj   = jup ? j : j - num_up_;

      auto& pair_mat = (jup == iup) ? (jup ? uu_triplet_mat_ : dd_triplet_mat_) : singlet_mat_;

      auto* dpsi_j  = jup ? up_dpsi_mat_[jj] : dn_dpsi_mat_[jj];
      auto* d2psi_j = jup ? up_d2psi_mat_[jj] : dn_d2psi_mat_[jj];

      dpsi_rows_(j, iat)  = 0.0;
      d2psi_rows_(j, iat) = 0.0;
      for (int k = 0; k < norb; k++)
        for (int l = 0; l < norb; l++)
        {
          dpsi_rows_(j, iat) -= tmp_psi_[k] * pair_mat(k, l) * dpsi_j[l];
          d2psi_rows_(j, iat) -= tmp_psi_[k] * pair_mat(k, l) * d2psi_j[l];
        }

      if (size_ == num_elec_ + 1)
      {
        const int idx         = size_ - 1;
        dpsi_rows_(idx, iat)  = -dpsi_rows_(iat, idx);
        d2psi_rows_(idx, iat) = -d2psi_rows_(iat, idx);
      }
    }
  }
  active_idx_ = -1;
  cur_ratio_  = 1.0;
}

PfaffianSTU::PsiValue PfaffianSTU::ratio(ParticleSet& P, int iat)
{
  UpdateMode  = ORB_PBYP_RATIO;
  active_idx_ = iat;
  {
    ScopedTimer local_timer(SPOVTimer);
    const int group = P.getGroupID(iat);
    sposets_[group]->evaluateValue(P, iat, tmp_psi_);
  }
  ScopedTimer local_timer(RatioTimer);

  const int norb = sposets_[0]->size();
  ValueVector row_update(size_, 0.0);
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
      for (int l = 0; l < norb; l++)
        row_update[j] += sign * tmp_psi_[k] * pair_mat(k, l) * psi_j[l];
  }

  if (size_ == num_elec_ + 1)
  {
    bool iup              = (iat < num_up_);
    int ii                = iup ? iat : iat - num_up_;
    row_update[num_elec_] = tmp_psi_[ii];
  }

  cur_ratio_ = calculateRatio(row_update);
  return cur_ratio_;
}

std::unique_ptr<WaveFunctionComponent> PfaffianSTU::makeClone(ParticleSet& tqp) const
{
  std::vector<std::unique_ptr<SPOSet>> sposet_clones;
  for (const auto& phi : sposets_)
    sposet_clones.emplace_back(phi->makeClone());
  const std::string opt_singlet    = (opt_singlet_) ? "yes" : "no";
  const std::string opt_uu_triplet = (opt_uu_triplet_) ? "yes" : "no";
  const std::string opt_dd_triplet = (opt_dd_triplet_) ? "yes" : "no";
  auto myclone =
      std::make_unique<PfaffianSTU>(tqp, std::move(sposet_clones), opt_singlet, opt_uu_triplet, opt_dd_triplet);

  //need to also copy data needed to actually calculate things. Should only be the
  //pairing matrices

  myclone->singlet_mat_    = this->singlet_mat_;
  myclone->uu_triplet_mat_ = this->uu_triplet_mat_;
  myclone->dd_triplet_mat_ = this->dd_triplet_mat_;
  myclone->myVars          = this->myVars;

  return myclone;
}

void PfaffianSTU::evaluateDerivatives(ParticleSet& P,
                                      const OptVariables& active,
                                      Vector<ValueType>& dlogpsi,
                                      Vector<ValueType>& dhpsioverpsi)
{}

void PfaffianSTU::evaluateDerivativesWF(ParticleSet& P, const OptVariables& active, Vector<ValueType>& dlogpsi)
{
  evaluateLog(P, P.G, P.L); //bring everything up to date

  const int norbs = singlet_mat_.rows();
  int iv          = 0;
  if (opt_singlet_)
  {
    for (int ip = 0; ip < norbs; ip++)
      for (int jp = ip; jp < norbs; jp++, iv++)
      {
        const int loc   = myVars.where(iv);
        ValueType deriv = 0;
        //singlet terms only effected by up, down pairs
        for (int ie = 0; ie < num_up_; ie++)
          for (int je = num_up_; je < num_elec_; je++)
          {
            ValueType val = up_psi_mat_(ie, ip) * dn_psi_mat_(je - num_up_, jp);
            if (ip != jp)
              val += up_psi_mat_(ie, jp) * dn_psi_mat_(je - num_up_, ip);
            deriv += psi_matinv_(je, ie) * val;
          }
        dlogpsi[loc] += deriv;
      }
  }

  if (opt_uu_triplet_)
  {
    for (int ip = 0; ip < norbs; ip++)
      for (int jp = ip + 1; jp < norbs; jp++, iv++)
      {
        const int loc   = myVars.where(iv);
        ValueType deriv = 0;
        for (int ie = 0; ie < num_up_; ie++)
          for (int je = ie + 1; je < num_up_; je++)
          {
            ValueType val = up_psi_mat_(ie, ip) * up_psi_mat_(je, jp) - up_psi_mat_(ie, jp) * up_psi_mat_(je, ip);
            deriv += psi_matinv_(je, ie) * val;
          }
        dlogpsi[loc] += deriv;
      }
  }

  if (opt_dd_triplet_)
  {
    for (int ip = 0; ip < norbs; ip++)
      for (int jp = ip + 1; jp < norbs; jp++, iv++)
      {
        const int loc   = myVars.where(iv);
        ValueType deriv = 0;
        for (int ie = num_up_; ie < num_elec_; ie++)
          for (int je = ie + 1; je < num_elec_; je++)
          {
            int ii        = ie - num_up_;
            int jj        = je - num_up_;
            ValueType val = dn_psi_mat_(ii, ip) * dn_psi_mat_(jj, jp) - dn_psi_mat_(ii, jp) * dn_psi_mat_(jj, ip);
            deriv += psi_matinv_(je, ie) * val;
          }
        dlogpsi[loc] += deriv;
      }
  }
}

void PfaffianSTU::resize()
{
  if (Bytes_in_WFBuffer > 0)
    throw std::runtime_error("PfaffianSTU just went out of sync with buffer");
  psi_mat_.resize(size_, size_);
  psi_matinv_.resize(size_, size_);
  dpsi_rows_.resize(size_, size_);
  d2psi_rows_.resize(size_, size_);
  psi_delta_.resize(size_);
  dpsi_new_.resize(size_);
  d2psi_new_.resize(size_);

  //now size the pairing function coefficient matrices matrices
  //up spos must be same size
  assert(sposets_[0]->size() == sposets_[1]->size());
  const int norbs = sposets_[0]->size();
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
  ValueMatrix tmp_mat(size_, size_);
  std::copy(psi_mat_.begin(), psi_mat_.end(), tmp_mat.begin());

  ValueType pf = 1.0;
  int sign     = 1;
  for (int i = 0; i < size_; i += 2)
  {
    sign *= rowPivot(tmp_mat, i);
    for (int j = i + 2; j < size_; j++)
    {
      ValueType fac = -tmp_mat(i, j) / tmp_mat(i, i + 1);
      for (int k = i + 1; k < size_; k++)
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
  ScopedTimer local_timer(InverseTimer);
  std::copy(psi_mat_.begin(), psi_mat_.end(), psi_matinv_.begin());
  invert_matrix(psi_matinv_, false);
  //active particle becomes invalid after inverse calculation
  active_idx_ = -1;
}

PfaffianSTU::ValueType PfaffianSTU::calculateRatio(const ValueVector& newvals)
{
  assert(active_idx_ >= 0);
  assert(newvals.size() == size_);
  std::transform(newvals.begin(), newvals.end(), psi_mat_[active_idx_], psi_delta_.begin(),
                 [](auto v1, auto v2) { return v1 - v2; });

  //exploting symmetry of matrices here. psi_matint_[i] gives a row, but I need to dot with column.
  //Since antisymmetric, add a sign to G and L
  ValueType ratio = -simd::dot(psi_matinv_[active_idx_], psi_delta_.data(), size_);
  return 1.0 + ratio;
}

void PfaffianSTU::updateInverse()
{
  const int n = size_;
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

void PfaffianSTU::buildOptVariables()
{
  myVars.clear();

  const int num_orbs = singlet_mat_.rows();

  auto registerParam = [this](const int i, const int j, const std::string& label, const ValueMatrix& pair_mat) {
    std::stringstream sstr;
    sstr << my_name_ << "_" << label << "_" << (i < 10 ? "0" : "") << (i < 100 ? "0" : "") << (i < 1000 ? "0" : "") << i
         << "_" << (j < 10 ? "0" : "") << (j < 100 ? "0" : "") << (j < 1000 ? "0" : "") << j << std::endl;
    myVars.insert(sstr.str(), std::real(pair_mat(i, j)));
  };

  if (opt_singlet_)
  {
    const std::string label = "ud";
    for (int i = 0; i < num_orbs; i++)
      for (int j = i; j < num_orbs; j++)
        registerParam(i, j, label, singlet_mat_);
  }

  if (opt_uu_triplet_)
  {
    const std::string label = "uu";
    for (int i = 0; i < num_orbs; i++)
      for (int j = i + 1; j < num_orbs; j++)
        registerParam(i, j, label, uu_triplet_mat_);
  }

  if (opt_dd_triplet_)
  {
    const std::string label = "dd";
    for (int i = 0; i < num_orbs; i++)
      for (int j = i + 1; j < num_orbs; j++)
        registerParam(i, j, label, dd_triplet_mat_);
  }
}

void PfaffianSTU::checkInVariablesExclusive(OptVariables& active)
{
  if (myVars.size())
    active.insertFrom(myVars);
}

void PfaffianSTU::checkOutVariables(const OptVariables& active) { myVars.getIndex(active); }


void PfaffianSTU::resetParametersExclusive(const OptVariables& active)
{
  for (int i = 0; i < myVars.size(); i++)
  {
    int loc   = myVars.where(i);
    myVars[i] = active[loc];
  }

  const int num_orbs = singlet_mat_.rows();
  int idx            = 0;
  if (opt_singlet_)
  {
    for (int i = 0; i < num_orbs; i++)
    {
      for (int j = i; j < num_orbs; j++, idx++)
      {
        singlet_mat_(i, j) = myVars[idx];
        if (j > i)
          singlet_mat_(j, i) = myVars[idx];
      }
    }
  }
  if (opt_uu_triplet_)
  {
    for (int i = 0; i < num_orbs; i++)
    {
      for (int j = i + 1; j < num_orbs; j++, idx++)
      {
        uu_triplet_mat_(i, j) = myVars[idx];
        uu_triplet_mat_(j, i) = -myVars[idx];
      }
    }
  }
  if (opt_dd_triplet_)
  {
    for (int i = 0; i < num_orbs; i++)
    {
      for (int j = i + 1; j < num_orbs; j++, idx++)
      {
        dd_triplet_mat_(i, j) = myVars[idx];
        dd_triplet_mat_(j, i) = -myVars[idx];
      }
    }
  }
  assert(idx == myVars.size());
}

void PfaffianSTU::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) { opt_obj_refs.push_back(*this); }

} // namespace qmcplusplus
