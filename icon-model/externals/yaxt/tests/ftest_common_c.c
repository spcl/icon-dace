/**
 * @file ftest_common_c.c
 *
 * @copyright Copyright  (C)  2017 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdbool.h>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>


#if defined HAVE_ALLOCA_H && HAVE_ALLOCA_H
#  include <alloca.h>
#elif defined __GNUC__ && !defined alloca
#  define alloca(size) __builtin_alloca(size)
#elif defined _AIX
#  define alloca(size) __alloca(size)
#elif defined _MSC_VER
#  include <malloc.h>
#  define alloca(size) _alloca(size)
#else
#  include <stddef.h>
void *alloca (size_t);
#endif

#define CF_USE_ALLOCA 1
#if defined __clang__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wreserved-id-macro"
#endif
#include "cfortran.h"
#if defined __clang__
#  pragma GCC diagnostic pop
#endif

#include "ftest_common.h"

#ifdef INT64_T_IS_LONG_LONG
// #define INT64 LONGLONG */
#define INT64V LONGLONGV
// #define INT64VV LONGLONGVV
// #define INT64VVV LONGLONGVVV
#elif defined INT64_T_IS_LONG
// #define INT64 LONG
#define INT64V LONGV
// #define INT64VV LONGVV
// #define INT64VVV LONGVVV
#endif


static bool
Cmp_dbl_arrays(int asize, const double *a, const double *b)
{
  assert(asize >= 0);
  return memcmp(a, b, (size_t)asize * sizeof(*a)) != 0;
}

FCALLSCFUN3(INT, Cmp_dbl_arrays, CMP_DBL_ARRAYS, cmp_dbl_arrays,
            INT, DOUBLEV, DOUBLEV)

static bool
Cmp_int16_arrays(int asize, const int16_t *a, const int16_t *b)
{
  assert(asize >= 0);
  return memcmp(a, b, (size_t)asize * sizeof(*a)) != 0;
}

FCALLSCFUN3(INT, Cmp_int16_arrays, CMP_INT16_ARRAYS, cmp_int16_arrays,
            INT, SHORTV, SHORTV)

static bool
Cmp_int32_arrays(int asize, const int32_t *a, const int32_t *b)
{
  assert(asize >= 0);
  return memcmp(a, b, (size_t)asize * sizeof(*a)) != 0;
}

FCALLSCFUN3(INT, Cmp_int32_arrays, CMP_INT32_ARRAYS, cmp_int32_arrays,
            INT, INTV, INTV)

static bool
Cmp_int64_arrays(int asize, const int64_t *a, const int64_t *b)
{
  assert(asize >= 0);
  return memcmp(a, b, (size_t)asize * sizeof(*a)) != 0;
}

FCALLSCFUN3(INT, Cmp_int64_arrays, CMP_INT64_ARRAYS, cmp_int64_arrays,
            INT, INT64V, INT64V)

static bool
pure_Print(const char *s)
{
  fprintf(stderr, "%s", s);
  return 1;
}

FCALLSCFUN1(INT, pure_Print, PURE_PRINT, pure_print, STRING)


/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
