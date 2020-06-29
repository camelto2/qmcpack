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

  xgrid = get_NUgrid_from_posgrid(kgridraw, 0);
  ygrid = get_NUgrid_from_posgrid(kgridraw, 1);
  zgrid = get_NUgrid_from_posgrid(kgridraw, 2);

  isGridComputed = true;
}

void NkParserBase::get_grid(NUgrid*& xgrid_i, NUgrid*& ygrid_i, NUgrid*& zgrid_i)
{
  if (!isGridComputed)
    compute_grid();
  xgrid_i = xgrid;
  ygrid_i = ygrid;
  zgrid_i = zgrid;
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

  Nx = xgrid->num_points;
  Ny = ygrid->num_points;
  Nz = zgrid->num_points;

  kgrid.resize(Nx * Ny * Nz, 0.0);
  nk.resize(Nx * Ny * Nz, 0.0);
  nkerr.resize(Nx * Ny * Nz, 0.0);

  for (IndexType i = 0; i < kgridraw.size(); i++)
  {
    nx = xgrid->reverse_map(xgrid, kgridraw[i][0]);
    ny = ygrid->reverse_map(ygrid, kgridraw[i][1]);
    nz = zgrid->reverse_map(zgrid, kgridraw[i][2]);

    newindex        = nx * Ny * Nz + ny * Nz + nz;
    kgrid[newindex] = kgridraw[i];
    nk[newindex]    = nkraw[i];
    nkerr[newindex] = nkerr_raw[i];
  }

  isNkComputed = true;
}


} // namespace qmcplusplus
