#include "FSUtilities.h"
#include "Utilities/RandomGenerator.h"
#include <algorithm>

namespace qmcplusplus
{
void get_gridinfo_from_posgrid(const std::vector<PosType>& posgridlist,
                               const IndexType& axis,
                               RealType& lx,
                               RealType& rx,
                               RealType& dx,
                               IndexType& Nx)
{
  std::vector<RealType> kx;
  kx.resize(posgridlist.size());

  for (IndexType i = 0; i < posgridlist.size(); i++)
    kx[i] = posgridlist[i][axis];

  std::vector<RealType>::iterator it;

  std::sort(kx.begin(), kx.end());

  it = std::unique(kx.begin(), kx.end());

  lx = *(kx.begin());
  rx = *(it - 1);
  Nx = it - kx.begin();
  dx = (rx - lx) / RealType(Nx - 1);
}

void getStats(const std::vector<RealType>& vals, RealType& avg, RealType& err, int start)
{
  avg   = 0.0;
  int n = 0;
  for (int i = start; i < vals.size(); i++)
  {
    avg += vals[i];
    n++;
  }
  avg /= n;
  err = 0.0;
  for (int i = start; i < vals.size(); i++)
  {
    err += (vals[i] - avg) * (vals[i] - avg);
  }
  err /= n;
  err = std::sqrt(err);
}

void getStats(const std::vector<RealType>& vals,
              const std::vector<RealType>& errs,
              RealType& avg,
              RealType& err,
              int start)
{
  avg             = 0.0;
  err             = 0.0;
  int n           = vals.size() - start;
  RealType varsum = 0.0;
  for (int i = start; i < vals.size(); i++)
  {
    avg += vals[i] / n;
    varsum += errs[i] * errs[i] / n;
  }
  err = std::sqrt(varsum / n);
}

void getStatsWithResampling(const std::vector<RealType>& vals,
                            const std::vector<RealType>& errs,
                            RealType& avg,
                            RealType& err,
                            int start,
                            int NumSamples)
{
  RandomGenerator rng;
  std::vector<RealType> avgs;
  for (int is = 0; is < NumSamples; is++)
  {
    std::vector<RealType> newVals;
    for (int i = start; i < vals.size(); i++)
    {
      FullPrecRealType u1, u2;
      //rng.generate_normal(&chi, 1);
      u1 = rng();
      u2 = rng();
      FullPrecRealType chi = std::sqrt(-2.0 * std::log(u1)) * std::cos(2 * M_PI * u2);
      newVals.push_back(vals[i] + chi * errs[i]);
    }
    RealType tmp1, tmp2;
    getStats(newVals, tmp1, tmp2);
    avgs.push_back(tmp1);
  }
  getStats(avgs, avg, err);
}

int estimateEquilibration(const std::vector<RealType>& vals, RealType frac)
{
  int idx = int(frac * vals.size());
  RealType avg, err;
  getStats(vals, avg, err, idx);
  int c3, c2, c1;
  c3 = vals.size();
  c2 = vals.size();
  c1 = vals.size();
  for (int i = vals.size() - 2; i >= 0; i--)
  {
    if ((vals[i] - avg) * (vals[i + 1] - avg) < 0.0)
    {
      c3 = c2;
      c2 = c1;
      c1 = i;
    }
  }
  if (c3 > frac * vals.size())
    c3 = int(frac * vals.size());
  return c3;
}

} // namespace qmcplusplus
