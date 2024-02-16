/**
 * @file utils.h
 * @brief Utlity functions
 *
 * Small general utility functions:
 *  - pointer-id conversion
 *  - hash
 *  - sorting
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include "core/ppm_xfuncs.h"
#include "core/core.h"

typedef double(*coordinate_pointer)[3];
typedef double const (* const const_coordinate_pointer)[3] ;

#define YAC_ASSERT(exp, msg) \
  {if(!((exp))) die(((msg)));}

#define YAC_ASSERT_F(exp, format, ...) \
  { \
    if(!((exp))) { \
      char msg_buffer[1024]; \
      int ret = snprintf( \
        msg_buffer, sizeof(msg_buffer), ((format)), __VA_ARGS__); \
      if ((ret >= 0) && ((size_t)ret < sizeof(msg_buffer))) \
        die(((msg_buffer))); \
      else \
        die("an error occured, but error message could not be generated"); \
    } \
  }

int yac_file_exists(const char * filename);

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
void yac_qsort_index(
  void * a, size_t count, size_t size,
  int (*compare)(void const *, void const *), size_t * idx);

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

struct yac_name_type_pair {
  const char * name;
  int type;
};

char const * yac_name_type_pair_get_name(
  struct yac_name_type_pair const * pairs, size_t count, int type);
int yac_name_type_pair_get_type(
  struct yac_name_type_pair const * pairs, size_t count, char const * name);

/* =======================================================================
   Macros
   ======================================================================= */

#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MIN(a,b) ((a) < (b) ? (a) : (b))

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

#endif // UTILS_H

