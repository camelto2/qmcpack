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
#include <fstream>
#include <vector>
#include <iomanip>
#include <assert.h>
#include <algorithm>
#include <cmath>
#include "QMCTools/QMCFiniteSize/SkParserHDF5.h"
#include "QMCTools/QMCFiniteSize/FSUtilities.h"
using namespace qmcplusplus;

const int NumSamples = 1024;

typedef SkParserHDF5::RealType RealType;
typedef SkParserHDF5::FullPrecRealType FullPrecRealType;
typedef SkParserHDF5::PosType PosType;
bool operator==(const PosType& a, const PosType& b)
{
  RealType eps = 1e-8;
  return (std::abs(a[0] - b[0]) < eps) && (std::abs(a[1] - b[1]) < eps) && (std::abs(a[2] - b[2]) < eps);
}

void getExtrapolatedEstimator(const std::vector<PosType>& vmc_kpts,
                              const std::vector<RealType>& vmc_vals,
                              const std::vector<RealType>& vmc_errs,
                              const std::vector<PosType>& dmc_kpts,
                              const std::vector<RealType>& dmc_vals,
                              const std::vector<RealType>& dmc_errs,
                              std::vector<RealType>& ext_vals,
                              std::vector<RealType>& ext_errs)
{
  assert(vmc_kpts.size() == dmc_kpts.size());
  for (int i = 0; i < vmc_kpts.size(); i++)
  {
    assert(vmc_kpts[i] == dmc_kpts[i]);
  }
  for (int i = 0; i < dmc_vals.size(); i++)
  {
    RealType val = 2.0 * dmc_vals[i] - vmc_vals[i];
    RealType err = std::sqrt((2.0 * dmc_errs[i]) * (2.0 * dmc_errs[i]) + vmc_errs[i] * vmc_errs[i]);
    ext_vals.push_back(val);
    ext_errs.push_back(err);
  }
}

void writeToFile(std::string filename,
                 const std::vector<PosType>& kpts,
                 const std::vector<RealType>& vals,
                 const std::vector<RealType>& errs)
{
  ofstream f;
  f.open(filename.c_str());
  f << "# kx ky kz val err" << std::endl;
  for (int i = 0; i < kpts.size(); i++)
  {
    f << kpts[i] << " " << vals[i] << " " << errs[i] << std::endl;
  }
  f.close();
}

//Get SK data for each twist.
//If it is weighted higher than one, it will be appended to the array "Weight" number of times
void getSkDataPerTwist(std::string filename,
                       int weight,
                       int nElec,
                       std::vector<PosType>& kpts_all,
                       std::vector<RealType>& vals_all,
                       std::vector<RealType>& errs_all)
{
  SkParserHDF5 skparser;
  skparser.parse(filename); //performs block averaging and equilibration estimate
  std::vector<PosType> kpts = skparser.get_grid_raw();
  std::vector<RealType> val = skparser.get_sk_raw();
  std::vector<RealType> err = skparser.get_skerr_raw();
  assert(kpts.size() == val.size());
  assert(kpts.size() == err.size());
  for (int j = 0; j < weight; j++)
  {
    for (int i = 0; i < kpts.size(); i++)
    {
      kpts_all.push_back(kpts[i]);
      vals_all.push_back(val[i] / nElec); //makes a normalized S(k)
      errs_all.push_back(err[i] / nElec);
    }
  }
}

void getAveragedData(std::vector<PosType>& kpts_all, std::vector<RealType>& vals_all, std::vector<RealType>& errs_all)
{
  std::vector<PosType> unique_kpts;
  std::vector<RealType> unique_vals;
  std::vector<RealType> unique_errs;
  std::vector<PosType>::iterator it;
  std::vector<std::vector<int>> idxs;
  for (int i = 0; i < kpts_all.size(); i++)
  {
    it = std::find(unique_kpts.begin(), unique_kpts.end(), kpts_all[i]);
    if (it != unique_kpts.end())
    {
      int idx = it - unique_kpts.begin();
      idxs[idx].push_back(i);
    }
    else
    {
      std::vector<int> idx;
      idx.push_back(i);
      unique_kpts.push_back(kpts_all[i]);
      idxs.push_back(idx);
    }
  }

  assert(unique_kpts.size() == idxs.size());
  for (int ik = 0; ik < unique_kpts.size(); ik++)
  {
    std::vector<RealType> vals, errs;
    for (int j = 0; j < idxs[ik].size(); j++)
    {
      vals.push_back(vals_all[idxs[ik][j]]);
      errs.push_back(errs_all[idxs[ik][j]]);
    }
    RealType avg, err;
    RealType avg1, err1;
    getStats(vals, errs, avg, err);
    unique_vals.push_back(avg);
    unique_errs.push_back(err);
  }

  kpts_all.clear();
  vals_all.clear();
  errs_all.clear();
  kpts_all = unique_kpts;
  vals_all = unique_vals;
  errs_all = unique_errs;
}

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
  int nElec = 1;
  std::vector<int> twistWeights;
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
    else if (a == "--nelec")
    {
      iargc++;
      nElec = std::stoi(std::string(argv[iargc]));
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
        twistWeights.push_back(std::stoi(argv[iargc]));
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

  if (computeSK)
  {
    std::vector<PosType> vmc_kpts;
    std::vector<RealType> vmc_sk;
    std::vector<RealType> vmc_skerr;
    std::vector<PosType> dmc_kpts;
    std::vector<RealType> dmc_sk;
    std::vector<RealType> dmc_skerr;
    std::vector<RealType> ext_sk;
    std::vector<RealType> ext_skerr;
    if (foundVMC)
    {
      std::cout << "Processing SK VMC files: " << std::endl;
      for (int tw = 0; tw < nTwists; tw++)
      {
        std::stringstream ss;
        ss << basename << ".g" << std::setfill('0') << std::setw(3) << tw << ".s" << std::setfill('0') << std::setw(3)
           << vmcIdx << ".stat.h5";
        std::cout << "  " << ss.str() << "  with weight: " << twistWeights[tw] << std::endl;
        getSkDataPerTwist(ss.str(), twistWeights[tw], nElec, vmc_kpts, vmc_sk, vmc_skerr);
      }
    }
    if (foundDMC)
    {
      std::cout << "Processing SK DMC files: " << std::endl;
      for (int tw = 0; tw < nTwists; tw++)
      {
        std::stringstream ss;
        ss << basename << ".g" << std::setfill('0') << std::setw(3) << tw << ".s" << std::setfill('0') << std::setw(3)
           << dmcIdx << ".stat.h5";
        std::cout << "  " << ss.str() << "  with weight: " << twistWeights[tw] << std::endl;
        getSkDataPerTwist(ss.str(), twistWeights[tw], nElec, dmc_kpts, dmc_sk, dmc_skerr);
      }
    }

    if (foundVMC)
    {
      std::string filename = "SK_VMC_twistavg.dat";
      std::cout << "Writing twist-averaged S(k) VMC estimator to " << filename << std::endl;
      getAveragedData(vmc_kpts, vmc_sk, vmc_skerr);
      writeToFile(filename, vmc_kpts, vmc_sk, vmc_skerr);
    }
    if (foundDMC)
    {
      std::string filename = "SK_DMC_twistavg.dat";
      std::cout << "Writing twist-averaged S(k) DMC estimator to " << filename << std::endl;
      getAveragedData(dmc_kpts, dmc_sk, dmc_skerr);
      writeToFile(filename, dmc_kpts, dmc_sk, dmc_skerr);
    }
    if (foundVMC && foundDMC)
    {
      std::string filename = "SK_EXTRAP_twistavg.dat";
      std::cout << "Writing twist-averaged S(k) extrapolated estimator to " << filename << std::endl;
      getExtrapolatedEstimator(vmc_kpts, vmc_sk, vmc_skerr, dmc_kpts, dmc_sk, dmc_skerr, ext_sk, ext_skerr);
      writeToFile(filename, dmc_kpts, ext_sk, ext_skerr);
    }
  }

  return 0;
}
