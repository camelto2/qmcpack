//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file RandomGenerator.h
 * @brief Declare a global Random Number Generator
 *
 * Selected among
 * - boost::random
 * - sprng
 * - math::random
 * qmcplusplus::Random() returns a random number [0,1)
 * For OpenMP is enabled, it is important to use thread-safe boost::random. Each
 * thread uses its own random number generator with a distinct seed. This prevents
 * a use of global lock which will slow down the applications significantly.
 */
#ifndef OHMMS_RANDOMGENERATOR
#define OHMMS_RANDOMGENERATOR
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cstdint>
// The definition of the fake RNG should always be available for unit testing
#include "Utilities/FakeRandom.h"
#ifdef HAVE_LIBBOOST
#include "Utilities/BoostRandom.h"
#else
#error -DHAVE_LIBBOOST is missing in the compile line. A cmake dependency fix is needed.
#endif

uint32_t make_seed(int i, int n);

namespace qmcplusplus
{
template<class RNG>
class RNGThreadSafe : public RNG
{
public:
  using result_type = typename RNG::result_type;

  result_type rand();

  /** return a random number [0,1)
   */
  result_type operator()();
};

extern template class RNGThreadSafe<FakeRandom>;
extern template class RNGThreadSafe<BoostRandom<float>>;
extern template class RNGThreadSafe<BoostRandom<double>>;

#ifdef USE_FAKE_RNG
using RandomGenerator_t = FakeRandom;
extern RNGThreadSafe<RandomGenerator_t> fake_random_global;
#define Random fake_random_global
#else
using RandomGenerator_t = BoostRandom<OHMMS_PRECISION_FULL>;
extern RNGThreadSafe<RandomGenerator_t> boost_random_global;
#define Random boost_random_global
#endif
} // namespace qmcplusplus

#endif
