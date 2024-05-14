//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ACCELBLAS_H
#define QMCPLUSPLUS_ACCELBLAS_H

#include "PlatformKinds.hpp"

namespace qmcplusplus
{

namespace compute
{

template<PlatformKind PL>
class BLAS;

}

} // namespace qmcplusplus

#endif
