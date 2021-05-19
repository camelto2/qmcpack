#ifndef QMC_FINITE_SIZE_H
#define QMC_FINITE_SIZE_H

#include "QMCApp/QMCAppBase.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "QMCTools/QMCFiniteSize/NkParserBase.h"
#include "Particle/ParticleSetPool.h"
#include "LongRange/LRCoulombSingleton.h"
#include "einspline/bspline_structs.h"
#include "einspline/nubspline_structs.h"

namespace qmcplusplus
{
/** Class to handle FS corrections
 *
 * Implements finite size corrections from Holzmann et al., PRB (2016)
 * Currently implements Eqn. (30), using a long-rang break up of the 
 * Coulomb interaction and a spline representation of S(k). 
 * S(k) is obtained from SkParserBase
 */
class QMCFiniteSize : public QMCAppBase, QMCTraits
{
public:
  typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
  typedef LRCoulombSingleton::GridType GridType;
  typedef LRCoulombSingleton::RadFunctorType RadFunctorType;
  typedef LRHandlerType::mRealType mRealType;
  typedef SkParserBase::Grid_t Grid_t;
  typedef QMCTraits::RealType RealType;
  typedef QMCTraits::FullPrecRealType FullPrecRealType;
  typedef QMCTraits::PosType PosType;
  QMCFiniteSize();
  QMCFiniteSize(SkParserBase* skparser_i, NkParserBase* nkparser_i);
  ~QMCFiniteSize(){};


  inline void setSkParser(SkParserBase* skparser_i) { skparser = skparser_i; };
  inline void setNkParser(NkParserBase* nkparser_i) { nkparser = nkparser_i; };
  bool validateXML();
  bool execute();

  void build_spherical_grid(IndexType mtheta, IndexType mphi);
  void getSkInfo(UBspline_3d_d* spline, vector<RealType>& symmatelem);
  UBspline_3d_d* getSkSpline(vector<RealType> sk, RealType limit = 1.0);
  NUBspline_3d_d* getNkSpline(vector<RealType> nk);

  RealType integrate_spline(NUBspline_1d_d* spline, RealType a, RealType b, IndexType N);
  NUBspline_1d_d* spline_clamped(vector<RealType>& grid, vector<RealType>& vals, RealType lVal, RealType rVal);

  void initialize();
  void initializeSkCorrection();
  void initializeNkCorrection();
  void calcPotentialCorrection();
  void calcKineticCorrection();
  void calcLeadingOrderCorrections();
  void summary();
  RealType calcPotentialDiscrete(vector<RealType> sk);
  RealType calcKineticDiscrete(vector<RealType> nk);
  RealType calcPotentialInt(vector<RealType> sk);
  RealType calcKineticInt(vector<RealType> nk);
  void executeSkCorrection();
  void executeNkCorrection();

private:
  SkParserBase* skparser;
  NkParserBase* nkparser;
  ParticleSetPool ptclPool;
  RealType myRcut;
  RealType myConst;
  ParticleSet* P;
  RealType h; //this is for finite differencing.
  vector<PosType> sphericalgrid;
  GridType* myGrid;
  std::unique_ptr<LRHandlerType> AA;
  RadFunctorType* rVs;
  bool processPWH(xmlNodePtr cur);
  void wfnPut(xmlNodePtr cur);
  void initBreakup();
  Grid_t gridx;
  Grid_t gridy;
  Grid_t gridz;
  NUgrid* nugridx;
  NUgrid* nugridy;
  NUgrid* nugridz;
  void printSkRawSphAvg(const vector<RealType>& sk);
  void printSkSplineSphAvg(UBspline_3d_d* spline);
  void printNkRawSphAvg(const vector<RealType>& nk);
  void printNkSplineSphAvg(NUBspline_3d_d* spline);
  RealType sphericalAvgSk(UBspline_3d_d* spline, RealType k);
  RealType sphericalAvgNk(NUBspline_3d_d* spline, RealType k);
  KContainer Klist;
  vector<TinyVector<int, OHMMS_DIM>> kpts;
  vector<RealType> SK_raw;
  vector<RealType> SKerr_raw;
  vector<RealType> SK;
  vector<RealType> SKerr;
  vector<RealType> NK_raw;
  vector<RealType> NKerr_raw;
  vector<RealType> NK;
  vector<RealType> NKerr;
  vector<PosType> NKkpts;
  vector<PosType> NKkpts_raw;
  map<RealType, vector<int>> kMap;
  RealType NKkmax;
  IndexType mtheta;
  IndexType mphi;
  IndexType NumSamples;
  RealType Ne, Vol, rs, rho;
  RealType tlo, tloerr, vlo, vloerr, Vfs, Vfserr, Tfs, Tfserr;
};
} // namespace qmcplusplus

#endif
