#ifndef NK_PARSER_ASCII_H
#define NK_PARSER_ASCII_H

#include "QMCTools/QMCFiniteSize/NkParserBase.h"

namespace qmcplusplus
{
class NkParserASCII : public NkParserBase
{
public:
  enum data_layout
  {
    KX,
    KY,
    KZ,
    NK,
    NKERR
  };
  void parse(const string& fname);

private:
  vector<vector<RealType>> read_nk_file(const string& fname);
  vector<PosType> get_grid_from_data(vector<vector<RealType>>& data);
  vector<RealType> get_nk_from_data(vector<vector<RealType>>& data);
  vector<RealType> get_nkerr_from_data(vector<vector<RealType>>& data);
};
} // namespace qmcplusplus

#endif
