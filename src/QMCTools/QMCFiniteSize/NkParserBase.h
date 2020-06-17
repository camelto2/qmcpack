#ifndef NK_PARSER_BASE_H
#define NK_PARSER_BASE_H

#include "Configuration.h"
namespace qmcplusplus
{
using namespace std;

/** Base class for Nk parser
 *
 * parse is the only pure virtual function that must be overridden
 * holds various information about n(k) and access functions used by
 * qmcfinitesize
 */
class NkParserBase : public QMCTraits
{
public:
  virtual void parse(const string& fname) = 0;

protected:
  bool isParseSuccess;
  vector<PosType> kgrid;
  vector<RealType> nk;
  vector<RealType> nk_err;
};

} // namespace qmcplusplus
#endif
