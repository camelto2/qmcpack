#include "QMCTools/QMCFiniteSize/NkParserASCII.h"
#include <fstream>

namespace qmcplusplus
{
typedef NkParserASCII::RealType RealType;
typedef NkParserASCII::PosType PosType;

vector<vector<RealType>> NkParserASCII::read_nk_file(const string& fname)
{
  vector<vector<RealType>> nkdata(0);

  vector<RealType> tmp(5);

  ifstream f;
  f.open(fname.c_str(), ifstream::in);

  string tmpstring;
  std::getline(f, tmpstring);

  while (!f.eof())
  {
    RealType x = 0, y = 0;
    f >> tmp[KX] >> tmp[KY] >> tmp[KZ] >> tmp[NK] >> tmp[NKERR];
    if (!f.eof())
      nkdata.push_back(tmp);
  }
  return nkdata;
}

vector<PosType> NkParserASCII::get_grid_from_data(vector<vector<RealType>>& filedata)
{
  vector<PosType> grid(filedata.size());
  for (int i = 0; i < filedata.size(); i++)
  {
    grid[i][0] = filedata[i][KX];
    grid[i][1] = filedata[i][KY];
    grid[i][2] = filedata[i][KZ];
  }

  return grid;
}

vector<RealType> NkParserASCII::get_nk_from_data(vector<vector<RealType>>& filedata)
{
  vector<RealType> val(filedata.size());
  for (int i = 0; i < filedata.size(); i++)
    val[i] = filedata[i][NK];

  return val;
}

vector<RealType> NkParserASCII::get_nkerr_from_data(vector<vector<RealType>>& filedata)
{
  vector<RealType> err(filedata.size());
  for (int i = 0; i < filedata.size(); i++)
    err[i] = filedata[i][NKERR];

  return err;
}

void NkParserASCII::parse(const string& fname)
{
  vector<vector<RealType>> rawdata(0);
  rawdata   = read_nk_file(fname);
  kgridraw  = get_grid_from_data(rawdata);
  nkraw     = get_nk_from_data(rawdata);
  nkerr_raw = get_nkerr_from_data(rawdata);

  isParseSuccess = true;
}
} // namespace qmcplusplus
