/**
 * @file test_arithmetic_long.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
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
#include <config.h>
#endif

#include <stdlib.h>

#include "tests.h"

#include "xt/xt_core.h"

#include <../src/xt_arithmetic_long.h>

static void
test_xlabs(void);

static void
test_xiimul(void);

static void
test_xlimul(void);

int main(void)
{
  test_xlabs();
  test_xiimul();
  test_xlimul();
  return TEST_EXIT_CODE;
}

#ifdef XT_LONG
#define XL(h,l) ((Xt_long)((Xt_ulong)((Xt_ulong)(h) << xt_int_bits) \
                           | (Xt_ulong)((Xt_ulong)(l)                   \
                                        & (Xt_ulong)(Xt_uint)(~(Xt_uint)0))))
#define XT(h,m,l) { .hi = (Xt_uint)(h), .midlo = (Xt_ulong)XL((m),(l)) }
#else
#define XL(h,l) { .hi = (Xt_uint)(h), .lo = (Xt_uint)(l) }
#define XT(h,m,l) { .hi = (Xt_uint)(h), .mid = (Xt_uint)(m), .lo = (Xt_uint)(l) }
#endif

static void
test_xlabs(void)
{
  struct tstV {
    Xt_long a, r;
  };
  static const struct tstV tst[] = {
    { XL(0, 0), XL(0, 0) }, /* |0| == 0 */
    { XL(0, 1), XL(0, 1) }, /* |-1| == 1 */
    { XL(-1, -1), XL(0, 1) }, /* |-1| == 1 */
    { XL(0, 21), XL(0, 21) }, /* |21| == 21 */
    { XL(-1, -21), XL(0, 21) }, /* |-21| == 21 */
  };
  size_t num_tests = sizeof (tst) / sizeof (tst[0]);
  for (size_t i = 0; i < num_tests; ++i) {
    Xt_long p = xlabs(tst[i].a);
    if (!xllcmp_eq(p, tst[i].r))
      PUT_ERR("error: failed long abs test");
  }
}

static void
test_xiimul(void)
{
  struct tstV {
    Xt_int a, b;
    Xt_long r;
  };
  static const struct tstV tst[] = {
    { 7, 3, XL(0, 21) }, /* 7 * 3 == 21 */
    { 2, -1, XL(-1, -2) }, /* 2 * -1 == -2 */
    { -1, -1, XL(0, 1) }, /* -1 * -1 == 1 */
    { 7, 5, XL(0, 35) },         /*  7 * 5 == 35 */
  };
  size_t num_tests = sizeof (tst) / sizeof (tst[0]);
  for (size_t i = 0; i < num_tests; ++i) {
    Xt_long p = xiimul(tst[i].a, tst[i].b);
    if (!xllcmp_eq(p, tst[i].r)) {
      fprintf(stderr, "i=%zu, p=(%llu,%llu), tst[%zu].r=(%llu,%llu)\n", i,
              (unsigned long long)xlhi(p),
              (unsigned long long)xllo(p),
              i,
              (unsigned long long)xlhi(tst[i].r),
              (unsigned long long)xllo(tst[i].r));
      PUT_ERR("error: failed long multiplication test");
    }
  }
}

static void
test_xlimul(void)
{
  struct tstV {
    Xt_long a;
    Xt_int b;
    Xt_tword r;
  };
  static const struct tstV tst[] = {
    { XL(0,7),  3, XT(0, 0, 21) }, /* 7 * 3 == 21 */
    { XL(0,2), -1, XT(-1, -1, -2) }, /* 2 * -1 == -2 */
    { XL(-1,-1), -1, XT(0, 0, 1) }, /* -1 * -1 == 1 */
    { XL(5, 6), 7, XT(0, 35, 42) },         /*  7 * 5 == 35 */
  };
  size_t num_tests = sizeof (tst) / sizeof (tst[0]);
  for (size_t i = 0; i < num_tests; ++i) {
    Xt_tword p = xlimul(tst[i].a, tst[i].b);
    if (!xttcmp_eq(p, tst[i].r)) {
      fprintf(stderr, "i=%zu, p=(%llu,%llu,%llu), "
              "tst[%zu].r=(%llu,%llu,%llu)\n", i,
              (unsigned long long)xthi(p),
              (unsigned long long)xtmid(p),
              (unsigned long long)xtlo(p),
              i,
              (unsigned long long)xthi(tst[i].r),
              (unsigned long long)xtmid(tst[i].r),
              (unsigned long long)xtlo(tst[i].r));
      PUT_ERR("error: failed long multiplication test");
    }
  }
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
