// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef UTILS_CORE_H
#define UTILS_CORE_H

#include "utils_common.h"

/** \example test_abort_c.c
 * This contains an example of how to use yac_abort_message.
 */

/** \example test_quicksort.c
 * This contains an example of how to use quicksort_index.
 */
void yac_quicksort_index ( int * a, size_t n, int * idx);
void yac_quicksort_index_yac_int_size_t ( yac_int * a, size_t n, size_t * idx);
void yac_quicksort_index_yac_int_yac_int ( yac_int * a, size_t n, yac_int * idx);
void yac_quicksort_index_size_t_yac_int ( size_t * a, size_t n, yac_int * idx);
void yac_quicksort_index_yac_int_uint64_t ( yac_int * a, size_t n, uint64_t * idx);
void yac_quicksort_index_yac_int_int ( yac_int * a, size_t n, int * idx);
void yac_quicksort_index_size_t_size_t ( size_t * a, size_t n, size_t * idx);
void yac_quicksort_index_uint64_t_size_t ( uint64_t * a, size_t n, size_t * idx);
void yac_quicksort_index_int_size_t ( int * a, size_t n, size_t * idx);
void yac_quicksort_index_size_t_int ( size_t * a, size_t n, int * idx);
void yac_quicksort_index_size_t_void_p ( size_t * a, size_t n, void * * idx);
void yac_quicksort_index_int_yac_int ( int * a, size_t n, yac_int * idx);
void yac_quicksort_index_int_double ( int * a, size_t n, double * idx);
void yac_quicksort_index_size_t_size_t_double (
  size_t * a, size_t n, size_t * b, double * c );
void yac_quicksort_index_yac_int_yac_int_double (
  yac_int * a, size_t n, yac_int * b, double * c );
void yac_quicksort_index_yac_int_yac_int_size_t (
  yac_int * a, size_t n, yac_int * b, size_t * c );
void yac_quicksort_index_int_size_t_size_t (
  int * a, size_t n, size_t * b, size_t * c );
void yac_quicksort_index_int_size_t_yac_int (
  int * a, size_t n, size_t * b, yac_int * c );

/** \example test_mergesort.c
 *
 * Natural Merge sort *
 *
 */
void yac_mergesort(void* base, size_t num, size_t size,
                   int (*compar)(const void*,const void*));

/**
 * remove duplicated entries from a list of integers
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_uint(unsigned * array, size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   unsigned prev = array[0];

   for (size_t i = 1; i < N; ++i) {

      if (array[i] == prev) continue;

      prev = array[i];
      ++pos;

      if (pos != i) array[pos] = array[i];
   }

   *n = pos + 1;
}

/**
 * remove duplicated entries from a list of integers
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_double(double * array, size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   double prev = array[0];

   for (size_t i = 1; i < N; ++i) {

      if (array[i] == prev) continue;

      prev = array[i];
      ++pos;

      if (pos != i) array[pos] = array[i];
   }

   *n = pos + 1;
}

/**
 * remove duplicated entries from a list of integers
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_size_t(size_t * array, size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   size_t prev = array[0];

   for (size_t i = 1; i < N; ++i) {

      if (array[i] == prev) continue;

      prev = array[i];
      ++pos;

      if (pos != i) array[pos] = array[i];
   }

   *n = pos + 1;
}

/**
 * remove duplicated entries from a list of integers
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_size_t_2(
  size_t (*array)[2], size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   size_t prev[2] = {array[0][0],
                     array[0][1]};

   for (size_t i = 1; i < N; ++i) {

      if ((array[i][0] == prev[0]) &&
          (array[i][1] == prev[1]))continue;

      prev[0] = array[i][0];
      prev[1] = array[i][1];
      ++pos;

      if (pos != i) {
        array[pos][0] = array[i][0];
        array[pos][1] = array[i][1];
      }
   }

   *n = pos + 1;
}

/**
 * remove duplicated entries from a list of integers
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_size_t_3(
  size_t (*array)[3], size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   size_t prev[3] = {array[0][0],
                     array[0][1],
                     array[0][2]};

   for (size_t i = 1; i < N; ++i) {

      if ((array[i][0] == prev[0]) &&
          (array[i][1] == prev[1]) &&
          (array[i][2] == prev[2]))continue;

      prev[0] = array[i][0];
      prev[1] = array[i][1];
      prev[2] = array[i][2];
      ++pos;

      if (pos != i) {
        array[pos][0] = array[i][0];
        array[pos][1] = array[i][1];
        array[pos][2] = array[i][2];
      }
   }

   *n = pos + 1;
}

/**
 * remove duplicated entries from a list of integers
 * @param[in,out] array array containing a sorted (ascending) list of integers
 * @param[in,out] n     number of entries in array
 */
static inline void yac_remove_duplicates_yac_int(
   yac_int * array, size_t * n) {

   size_t const N = *n;
   size_t pos = 0;

   if (N == 0) return;

   yac_int prev = array[0];

   for (size_t i = 1; i < N; ++i) {

      if (array[i] == prev) continue;

      prev = array[i];
      ++pos;

      if (pos != i) array[pos] = array[i];
   }

   *n = pos + 1;
}

#define ASSERT(c) \
if (!(c)) {\
   fprintf(stderr, "### Assertion violation: %s in %s:%d\n",\
           #c, __FILE__, __LINE__);\
   abort ();\
}

#define COPY_DATA(data, count) \
  (memcpy( \
    xmalloc((size_t)(count) * sizeof(*(data))), \
    (data), (size_t)(count) * sizeof(*(data))))

#endif // UTILS_CORE_H

