#include "QMCTools/QMCFiniteSize/NkParserHDF5.h"
#include "QMCTools/QMCFiniteSize/FSUtilities.h"

namespace qmcplusplus
{
void NkParserHDF5::parse(const string& fname)
{
  bool result = statfile.open(fname);
  if (!result)
  {
    cout << "NkParserHDF5::parse could not open " << fname << endl;
    exit(1);
  }

  //read the kpoints
  vector<int> readShape;
  statfile.getShape<int>("nofk/kpoints", readShape);
  assert(readShape[1] == 3);
  int nKpts = readShape[0];
  array<int, 2> kdims{nKpts, 3};
  vector<RealType> ktmp;
  statfile.readSlabReshaped(ktmp, kdims, "nofk/kpoints");
  for (int ik = 0; ik < nKpts; ik++)
  {
    PosType k;
    k[0] = ktmp[3 * ik];
    k[1] = ktmp[3 * ik + 1];
    k[2] = ktmp[3 * ik + 2];
    kgridraw.push_back(k);
  }

  statfile.getShape<int>("nofk/value", readShape);
  assert(readShape[1] == nKpts);
  vector<RealType> nofk_tmp;
  array<int, 2> nofk_dims{readShape[0], nKpts};
  statfile.readSlabReshaped(nofk_tmp, nofk_dims, "nofk/value");
  int nBlocks = readShape[0];

  for (int ik = 0; ik < nKpts; ik++)
  {
    vector<RealType> block_data;
    for (int ib = 0; ib < nBlocks; ib++)
      block_data.push_back(nofk_tmp[nKpts * ib + ik]);
    int ieq = estimateEquilibration(block_data);
    RealType avg, err;
    getStats(block_data, avg, err, ieq);
    nkraw.push_back(avg);
    nkerr_raw.push_back(err);
  }

  statfile.close();

  isParseSuccess = true;

}

} // namespace qmcplusplus
