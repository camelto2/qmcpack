#ifndef NK_PARSER_BASE_H
#define NK_PARSER_BASE_H

#include "Configuration.h"
#include "Numerics/OneDimGridBase.h"
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
  typedef LinearGrid<RealType> Grid_t;
  NkParserBase();

  virtual void parse(const string& fname) = 0;
  vector<PosType> get_grid_raw() { return kgridraw; }
  vector<RealType> get_nk_raw() { return nkraw; }
  vector<RealType> get_nkerr_raw() { return nkerr_raw; }

  void compute_grid();
  void get_grid(Grid_t& xgrid, Grid_t& ygrid, Grid_t& zgrid);
  void get_nk(vector<RealType>& nk, vector<RealType>& nkerr);
  void compute_nk();
  inline bool has_grid() { return hasGrid; }
  void set_grid(const vector<PosType>& grid);


protected:
  bool isParseSuccess;
  bool isGridComputed;
  bool isNkComputed;
  bool hasGrid;
  vector<PosType> kgridraw;
  vector<RealType> nkraw;
  vector<RealType> nkerr_raw;

  Grid_t xgrid;
  Grid_t ygrid;
  Grid_t zgrid;

  vector<RealType> nk;
  vector<RealType> nkerr;
  vector<PosType> kgrid;
};

} // namespace qmcplusplus
#endif
