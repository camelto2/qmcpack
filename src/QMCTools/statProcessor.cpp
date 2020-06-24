//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by:  Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

void usage()
{
  std::cout << "Usage: statProcessor [basename] [options]" << std::endl;
  std::cout << "  [basename]: basename for qmc stat file, e.g. basename.g001.s003.stat.h5" << std::endl;
  std::cout << "  [options]" << std::endl;
  std::cout << "    --vmc: vmcIndex" << std::endl;
  std::cout << "    --dmc: dmcIndex" << std::endl;
  std::cout << "    --ntwists: Ntwists" << std::endl;
  std::cout << "    --weights: w_1 w_2 ... w_ntwists. Default: 1" << std::endl;
}

int main(int argc, char** argv)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right, std::ios::adjustfield);
  std::cout.precision(12);

  if (argc == 1)
  {
    usage();
    return 0;
  }

  std::string basename = std::string(argv[1]);
  int vmcIdx;
  bool foundVMC = false;
  int dmcIdx;
  bool foundDMC = false;
  int nTwists;
  std::vector<double> twistWeights;
  bool computeSK = false;
  bool computeNK = false;

  int iargc = 2;
  while (iargc < argc)
  {
    std::string a(argv[iargc]);
    if (a == "--vmc")
    {
      foundVMC = true;
      iargc++;
      vmcIdx = std::stoi(std::string(argv[iargc]));
    }
    else if (a == "--dmc")
    {
      foundDMC = true;
      iargc++;
      dmcIdx = std::stoi(std::string(argv[iargc]));
    }
    else if (a == "--ntwists")
    {
      iargc++;
      nTwists = std::stoi(std::string(argv[iargc]));
    }
    else if (a == "--SK")
    {
      computeSK = true;
    }
    else if (a == "--NK")
    {
      computeNK = true;
    }
    iargc++;
  }
  iargc = 1;
  while (iargc < argc)
  {
    std::string a(argv[iargc]);
    if (a == "--weights")
    {
      iargc++;
      for (int tw = 0; tw < nTwists; tw++)
      {
        twistWeights.push_back(std::stod(argv[iargc]));
        iargc++;
      }
    }
    iargc++;
  }
  if (twistWeights.size() == 0)
  {
    for (int tw = 0; tw < nTwists; tw++)
    {
      twistWeights.push_back(1.0);
    }
  }

  if (foundVMC)
  {
    std::cout << "Processing VMC stat.h5 files: " << std::endl;
    for (int tw = 0; tw < nTwists; tw++)
    {
      std::stringstream ss;
      ss << basename << ".g" << std::setfill('0') << std::setw(3) << tw << ".s" << std::setfill('0') << std::setw(3)
         << vmcIdx << ".stat.h5";
      std::cout << ss.str() << " weight: " << twistWeights[tw] << std::endl;
    }
    if (computeSK)
      std::cout << "Output SK VMC estimator: SK_VMC_twistavg.dat" << std::endl;
    if (computeNK)
      std::cout << "Output NK VMC estimator: NK_VMC_twistavg.dat" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
  }
  if (foundDMC)
  {
    std::cout << "Processing DMC stat.h5 files: " << std::endl;
    for (int tw = 0; tw < nTwists; tw++)
    {
      std::stringstream ss;
      ss << basename << ".g" << std::setfill('0') << std::setw(3) << tw << ".s" << std::setfill('0') << std::setw(3)
         << dmcIdx << ".stat.h5";
      std::cout << ss.str() << " weight: " << twistWeights[tw] << std::endl;
    }
    if (computeSK)
      std::cout << "Output SK mixed estimator: SK_DMC_twistavg.dat" << std::endl;
    if (computeNK)
      std::cout << "Output NK mixed estimator: NK_DMC_twistavg.dat" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
  }
  if (foundVMC && foundDMC)
  {
    if (computeSK)
      std::cout << "Output SK extrapolated estimator: SK_EXTRAP_twistavg.dat" << std::endl;
    if (computeNK)
      std::cout << "Output NK extrapolated estimator: NK_EXTRAP_twistavg.dat" << std::endl;
  }


  return 0;
}
