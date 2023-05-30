//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file CrystalLattice.cpp
 *@brief Member function definitions of CrystalLattice<T,D>, included by CrystalLattice.h
 */
#include "LatticeAnalyzer.h"

namespace qmcplusplus
{
/*! \fn CrystalLattice::CrystalLattice()
 *  Default constructor. Initialized to a \p 1x1x1 cubic supercell.
 */
template<class T, unsigned D>
CrystalLattice<T, D>::CrystalLattice()
{
  explicitly_defined = false;
  BoxBConds          = 0;
  VacuumScale        = 1.0;
  R.diagonal(1e10);
  G      = R;
  M      = R;
  Volume = 1;
  reset();
}

template<class T, unsigned D>
template<class TT>
void CrystalLattice<T, D>::set(const Tensor<TT, D>& lat)
{
  explicitly_defined = true;
  R                  = lat;
  reset();
}

template<class T, unsigned D>
void CrystalLattice<T, D>::reset()
{
  G      = inverse(R); //G = transpose(Inverse(R));
  Gt     = transpose(G);
  Volume = std::abs(det(R));
  //M = dot(transpose(R),R);
  M   = dot(R, transpose(R));
  T t = TWOPI * TWOPI;
  Mg  = t * dot(transpose(G), G);
  for (int i = 0; i < D; ++i)
    for (int j = 0; j < D; ++j)
      Rv[i][j] = R(i, j);
  for (int i = 0; i < D; ++i)
    for (int j = 0; j < D; ++j)
      Gv[i][j] = G(j, i);
  for (int i = 0; i < D; ++i)
  {
    Length[i]        = std::sqrt(dot(Rv[i], Rv[i]));
    OneOverLength[i] = 1.0 / Length[i];
  }
  Center = 0.0;
  for (int i = 0; i < D; ++i)
    Center += Rv[i];
  Center *= .5;
  //analysis of a lattice using LatticeAnalyzer
  LatticeAnalyzer<T, D> ldesc;
  SuperCellEnum        = ldesc(BoxBConds);
  DiagonalOnly         = ldesc.isDiagonalOnly(R);
  ABC                  = ldesc.calcSolidAngles(Rv, OneOverLength);
  WignerSeitzRadius    = ldesc.calcWignerSeitzRadius(Rv);
  WignerSeitzRadius_G  = ldesc.calcWignerSeitzRadius(Gv);
  SimulationCellRadius = ldesc.calcSimulationCellRadius(Rv);
  // set equal WignerSeitzRadius and SimulationCellRadius when they are very close.
  if (WignerSeitzRadius > SimulationCellRadius &&
      WignerSeitzRadius - SimulationCellRadius <= WignerSeitzRadius * std::numeric_limits<float>::epsilon() * 2)
    WignerSeitzRadius = SimulationCellRadius;
  CellRadiusSq = SimulationCellRadius * SimulationCellRadius;

  //calclate neighbor cells
  int counter = 0;
  int nsearch = 1;
  int numcell = 2 * nsearch + 1;
  for (int aa = -nsearch; aa <= nsearch; aa++)
    for (int bb = -nsearch; bb <= nsearch; bb++)
      for (int cc = -nsearch; cc <= nsearch; cc++)
      {
        std::vector<T> tmp(D);
        neighbor_cells.push_back(tmp);
        for (int d = 0; d < D; d++)
          neighbor_cells[counter][d] = aa * Rv[0][d] +  bb * Rv[1][d] + cc * Rv[2][d];
        counter++;
      }

  std::vector<std::vector<T>> cp;
  for (int i = 0; i < 3; i++)
  {
    std::vector<T> tmp(3);
    cp.push_back(tmp);
  }
  cp[0][0]=(Rv[1][1]*Rv[2][2]-Rv[1][2]*Rv[2][1]);
  cp[0][1]=(Rv[1][2]*Rv[2][0]-Rv[1][0]*Rv[2][2]); 
  cp[0][2]=(Rv[1][0]*Rv[2][1]-Rv[1][1]*Rv[2][0]);

  cp[1][0]=(Rv[2][1]*Rv[0][2]-Rv[2][2]*Rv[0][1]);
  cp[1][1]=(Rv[2][2]*Rv[0][0]-Rv[2][0]*Rv[0][2]);
  cp[1][2]=(Rv[2][0]*Rv[0][1]-Rv[2][1]*Rv[0][0]);

  cp[2][0]=(Rv[0][1]*Rv[1][2]-Rv[0][2]*Rv[1][1]);
  cp[2][1]=(Rv[0][2]*Rv[1][0]-Rv[0][0]*Rv[1][2]);
  cp[2][2]=(Rv[0][0]*Rv[1][1]-Rv[0][1]*Rv[1][0]);

  smallest_height = 1e99;
  for (int i = 0; i < 3; i++)
  {
    T tmph = 0;
    T l = 0;
    for (int j = 0; j < 3; j++)
    {
        tmph += cp[i][j] * Rv[i][j];
        l += cp[i][j] * cp[i][j];
    }
    tmph = std::abs(tmph)/std::sqrt(l);
    if (tmph < smallest_height)
      smallest_height = tmph;
  }

}

/*! \fn  CrystalLattice<T,D>::operator*=(T sc)
 *  \param sc A scaling factor.
 *  \brief Rescale this supercell by a scalar.
 */
template<class T, unsigned D>
CrystalLattice<T, D>& CrystalLattice<T, D>::operator*=(T sc)
{
  R *= sc;
  reset();
  return *this;
}

template<class T, unsigned D>
void CrystalLattice<T, D>::print(std::ostream& os, int level) const
{
  /*\note level == 0: print only the lattice vectors
   *      level == 1: lattice vectors, boundary conditions, grid
   *      level == 2: + all the internal values
   */
  std::string unit_name = "bohr";

  std::string lattice_name = "  Lattice (" + unit_name + "):";
  std::string pad(lattice_name.length(), ' ');
  os << lattice_name;
  for (int i = 0; i < D; ++i)
  {
    if (i > 0)
    {
      os << pad;
    }
    os << Rv[i] << std::endl;
  }
  if (level > 0)
  {
    os << std::endl;
    os << "  Boundary Conditions: ";
    for (int i = 0; i < D; ++i)
    {
      if (BoxBConds[i])
        os << " p ";
      else
        os << " n ";
    }
    os << std::endl;
    if (VacuumScale != 1.0)
      os << "  Vacuum scale: " << VacuumScale << std::endl;
  }
  if (level > 1)
  {
    os << std::endl;
    os << "  Volume (bohr^3) = " << Volume << std::endl;
    os << std::endl;
    os << "  Reciprocal vectors without 2*pi.\n";
    for (int i = 0; i < D; ++i)
      os << "    g_" << i + 1 << " = " << Gv[i] << std::endl;
    os << std::endl;
    os << "  Metric tensor in real-space.\n";
    for (int i = 0; i < D; ++i)
    {
      os << "    h_" << i + 1 << " = ";
      for (int j = 0; j < D; ++j)
      {
        os << M(i, j) << " ";
      }
      os << std::endl;
    }
    os << std::endl;
    os << "  Metric tensor in g-space.\n";
    for (int i = 0; i < D; ++i)
    {
      os << "    h_" << i + 1 << " = ";
      for (int j = 0; j < D; ++j)
      {
        os << Mg(i, j) << " ";
      }
      os << std::endl;
    }
  }
}

template<class T, unsigned D>
inline bool operator==(const CrystalLattice<T, D>& lhs, const CrystalLattice<T, D>& rhs)
{
  for (int i = 0; i < D * D; ++i)
    if (std::abs(lhs.R[i] - rhs.R[i]) > std::numeric_limits<T>::epsilon())
      return false;
  return true;
}

template<class T, unsigned D>
inline bool operator!=(const CrystalLattice<T, D>& lhs, const CrystalLattice<T, D>& rhs)
{
  return !(lhs == rhs);
}

} // namespace qmcplusplus
