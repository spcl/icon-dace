// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef UTILS_COMMON_H
#define UTILS_COMMON_H

#include "yac_config.h"

#ifdef YAC_OPENMP_ENABLED
#include <omp.h>
#endif

#include <stdlib.h>
#include "ppm/ppm_xfuncs.h"
#include "ppm/core.h"
#include "yac_types.h"
#include "yac_assert.h"

/**
 * remove duplicated entries from a list of integers
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_int(int * array, size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   int prev = array[0];

   for (size_t i = 1; i < N; ++i) {

      if (array[i] == prev) continue;

      prev = array[i];
      ++pos;

      if (pos != i) array[pos] = array[i];
   }

   *n = pos + 1;
}

struct yac_name_type_pair {
  const char * name;
  int type;
};

char const * yac_name_type_pair_get_name(
  struct yac_name_type_pair const * pairs, size_t count, int type);
int yac_name_type_pair_get_type(
  struct yac_name_type_pair const * pairs, size_t count, char const * name);

void yac_qsort_index(
  void * a, size_t count, size_t size,
  int (*compare)(void const *, void const *), size_t * idx);

/* =======================================================================
   Macros
   ======================================================================= */

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifdef YAC_OPENMP_ENABLED
#define YAC_OMP_PARALLEL _Pragma("omp parallel")
#define YAC_OMP_FOR _Pragma("omp for")
#else
#define YAC_OMP_PARALLEL
#define YAC_OMP_FOR
#endif

#endif // UTILS_H

