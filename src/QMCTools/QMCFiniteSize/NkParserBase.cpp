#include "QMCTools/QMCFiniteSize/NkParserBase.h"
#include "QMCTools/QMCFiniteSize/FSUtilities.h"

namespace qmcplusplus
{
NkParserBase::NkParserBase() : isParseSuccess(false), isGridComputed(false), isNkComputed(false), hasGrid(false)
{
  nkraw.resize(0);
  nk.resize(0);
  nkerr.resize(0);
  nkerr_raw.resize(0);
  kgridraw.resize(0);
  kgrid.resize(0);
}

void NkParserBase::compute_grid()
{
  if (!isParseSuccess)
    APP_ABORT("NkParserBase::compute_grid(..) : Initial parse failed");

  RealType lx(0), rx(0);
  RealType ly(0), ry(0);
  RealType lz(0), rz(0);

  RealType dx(0), dy(0), dz(0);
  IndexType Nx(0), Ny(0), Nz(0);

  get_gridinfo_from_posgrid(kgridraw, 0, lx, rx, dx, Nx);
  get_gridinfo_from_posgrid(kgridraw, 1, ly, ry, dy, Ny);
  get_gridinfo_from_posgrid(kgridraw, 2, lz, rz, dz, Nz);

  kgrid.resize(Nx * Ny * Nz);

  xgrid.set(lx, rx, Nx);
  ygrid.set(ly, ry, Ny);
  zgrid.set(lz, rz, Nz);

  isGridComputed = true;
}

void NkParserBase::get_grid(Grid_t& xgrid_i, Grid_t& ygrid_i, Grid_t& zgrid_i)
{
  if (!isGridComputed)
    compute_grid();
  xgrid_i.set(xgrid.rmin(), xgrid.rmax(), xgrid.size());
  ygrid_i.set(ygrid.rmin(), ygrid.rmax(), ygrid.size());
  zgrid_i.set(zgrid.rmin(), zgrid.rmax(), zgrid.size());
}

void NkParserBase::set_grid(const vector<PosType>& kgridraw_i)
{
  if (kgridraw_i.size() != nkraw.size())
    APP_ABORT("NkParserBase::set_grid: n(k) and k-grid don't match");
  kgridraw = kgridraw_i;
  compute_grid();
}

void NkParserBase::get_nk(vector<RealType>& nk_i, vector<RealType>& nkerr_i)
{
  if (!isNkComputed)
    compute_nk();
  nk_i    = nk;
  nkerr_i = nkerr;
}

void NkParserBase::compute_nk()
{
  if (!isParseSuccess)
    APP_ABORT("NkParserBase::compute_nk(..) : Initial parse failed");

  if (kgridraw.size() != nkraw.size())
    APP_ABORT("NkParserBase::compute_nk(..) : Kgrid and nofk not the same size");
  if (!isGridComputed)
    compute_grid();

  IndexType nx(0), ny(0), nz(0);
  IndexType Nx(0), Ny(0), Nz(0);
  IndexType newindex(0);

  Nx = xgrid.size();
  Ny = ygrid.size();
  Nz = zgrid.size();

  nk.resize(Nx * Ny * Nz, 0.0);
  nkerr.resize(Nx * Ny * Nz, 0.0);

  for (IndexType i = 0; i < kgridraw.size(); i++)
  {
    nx = xgrid.getIndex(kgridraw[i][0]);
    ny = ygrid.getIndex(kgridraw[i][1]);
    nz = zgrid.getIndex(kgridraw[i][2]);

    newindex        = nx * Ny * Nz + ny * Nz + nz;
    kgrid[newindex] = kgridraw[i];
    nk[newindex]    = nkraw[i];
    nkerr[newindex] = nkerr_raw[i];
  }

  isNkComputed = true;
}


} // namespace qmcplusplus
