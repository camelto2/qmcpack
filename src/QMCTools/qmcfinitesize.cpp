//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "Message/Communicate.h"
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGenerator.h"
#include "Particle/ParticleSetPool.h"

#include "QMCApp/QMCAppBase.h"
#include "QMCTools/QMCFiniteSize/QMCFiniteSize.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "QMCTools/QMCFiniteSize/SkParserASCII.h"
#include "QMCTools/QMCFiniteSize/SkParserScalarDat.h"
#include "QMCTools/QMCFiniteSize/SkParserHDF5.h"
#include "QMCTools/QMCFiniteSize/NkParserBase.h"
#include "QMCTools/QMCFiniteSize/NkParserHDF5.h"
#include "QMCTools/QMCFiniteSize/NkParserASCII.h"

#include "Numerics/OneDimGridBase.h"

//Purpose of this routine is to compute the finite size effects
//for a given simulation cell from post-processed QMC Data.  For
//the potential, this is done by splining the structure factor
//and performing the integral given in Holzmann et al., PRB 035126 (2016)
//Assuming you have provided a long-ranged jastrow in the input xml, this will also
//calculate the kinetic energy correction from Holzmann et al.
//
//Input:  Cell geometry.  Ion positions, cell geometries, etc, are taken from main.xml file.
//                        Code recognizes ESHDF5 declarations in main.xml.
//                        Will use kspace jastrow to compute kinetic correction if available.
//        S(k):  This is the electron-electron fluctuation structure factor.
//
//Returns: (E(N=infty)-E(N)) for the given simulation cell.

using namespace qmcplusplus;
typedef QMCTraits::RealType RealType;
typedef QMCTraits::PosType PosType;
typedef SkParserBase::Grid_t Grid_t;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc, argv);
  Random.init(0, 1, -1);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right, std::ios::adjustfield);
  std::cout.precision(12);

  std::unique_ptr<SkParserBase> skparser(nullptr);
  std::unique_ptr<NkParserBase> nkparser(nullptr);
  int iargc = 2;

  while (iargc + 1 < argc)
  {
    std::string a(argv[iargc]);
    std::string anxt(argv[iargc + 1]);
    if (a == "--SKascii")
    {
      skparser = std::make_unique<SkParserASCII>();
      skparser->parse(anxt);
    }
    else if (a == "--scalardat")
    {
      skparser = std::make_unique<SkParserScalarDat>();
      skparser->parse(anxt);
    }
    else if (a == "--SKhdf5")
    {
      skparser = std::make_unique<SkParserHDF5>();
      skparser->parse(anxt);
    }
    else if (a == "--NKhdf5")
    {
      nkparser = std::make_unique<NkParserHDF5>();
      nkparser->parse(anxt);
    }
    else if (a == "--NKascii")
    {
      nkparser = std::make_unique<NkParserASCII>();
      nkparser->parse(anxt);
    }
    else if (a == "--help")
    {
      std::cout << "Usage:  qmcfinitesize [main.xml] --[skformat] [SK_FILE]\n";
      std::cout << "  [skformat]\n";
      std::cout << "    --ascii:      S(k) given in kx ky kz sk sk_err format.  Header necessary.\n";
      std::cout << "    --scalardat:  File containing skall elements with energy.pl output format.\n";
      std::cout << "    --hdf5:       stat.h5 file containing skall data.\n";
      return 0;
    }
    iargc++;
  }

  QMCFiniteSize qmcfs(skparser.get(), nkparser.get());
  qmcfs.parse(std::string(argv[1]));
  qmcfs.validateXML();
  qmcfs.execute();

  // Jobs done. Clean up.
  OHMMS::Controller->finalize();
  return 0;
}
