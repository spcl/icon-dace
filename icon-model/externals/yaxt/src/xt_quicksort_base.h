/**
 * @file xt_quicksort_base.h
 * @brief macros to create quicksort implementations
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
/*
 * xt_quicksort_int was derived from the ScalES-PPM library code,
 * which is derived from genometools, which in turn is derived from
 * the FreeBSD libc.
 */
/*
  Modifications for integration with genometools
  2008 Thomas Jahns <Thomas.Jahns@gmx.net>

  The advertising clause 3. was removed due to the corresponding
  revoke by William Hoskins on July 22, 1999.
  <ftp://ftp.cs.berkeley.edu/pub/4bsd/README.Impt.License.Change>
*/
/*-
 * Copyright (c) 1992, 1993
 *        The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */
#ifndef XT_QUICKSORT_BASE_H
#define XT_QUICKSORT_BASE_H
#define TOKEN_PASTE(a,b) a##_##b
#define NAME_COMPOSE(a,b) TOKEN_PASTE(a,b)
#endif

#ifndef SORT_TYPE
#error "must define type to sort on"
#endif

#ifndef SORT_TYPE_SUFFIX
#error "must define suffix for type to name functions"
#endif

#ifndef SORT_TYPE_CMP_LT
#error "must define macro to compare SORT_TYPE for less than relation"
#endif

#ifndef SORT_TYPE_CMP_LE
#error "must define macro to compare SORT_TYPE for less than or equal relation"
#endif

#ifndef SORT_TYPE_CMP_EQ
#error "must define macro to compare SORT_TYPE for equality"
#endif

#ifndef XT_SORTFUNC_DECL
#define XT_SORTFUNC_DECL
#define XT_SORTFUNC_DECL_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ARGS_DECL
/* these declarations are appended to parameters, defaults to nothing */
#define XT_SORT_EXTRA_ARGS_DECL
#define XT_SORT_EXTRA_ARGS_DECL_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ARGS_PASS
/* determines what is passed to the parameters declared in
 * XT_SORT_EXTRA_ARGS_DECL, defaults to nothing */
#define XT_SORT_EXTRA_ARGS_PASS
#define XT_SORT_EXTRA_ARGS_PASS_UNDEF
#endif

#ifndef XT_SORT_VECSWAP_EXTRA_ARGS_DECL
#define XT_SORT_VECSWAP_EXTRA_ARGS_DECL XT_SORT_EXTRA_ARGS_DECL
#define XT_SORT_VECSWAP_EXTRA_ARGS_DECL_UNDEF
#endif

#ifndef XT_SORT_VECSWAP_EXTRA_ARGS_PASS
#define XT_SORT_VECSWAP_EXTRA_ARGS_PASS XT_SORT_EXTRA_ARGS_PASS
#define XT_SORT_VECSWAP_EXTRA_ARGS_PASS_UNDEF
#endif

#ifndef XT_SORT_MED3_EXTRA_ARGS_DECL
#define XT_SORT_MED3_EXTRA_ARGS_DECL XT_SORT_EXTRA_ARGS_DECL
#define XT_SORT_MED3_EXTRA_ARGS_DECL_UNDEF
#endif

#ifndef XT_SORT_MED3_EXTRA_ARGS_PASS
#define XT_SORT_MED3_EXTRA_ARGS_PASS XT_SORT_EXTRA_ARGS_PASS
#define XT_SORT_MED3_EXTRA_ARGS_PASS_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ARGS_SWAP
#define XT_SORT_EXTRA_ARGS_SWAP(i,j)
#define XT_SORT_EXTRA_ARGS_SWAP_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ARGS_ADVANCE
#define XT_SORT_EXTRA_ARGS_ADVANCE(adv)
#define XT_SORT_EXTRA_ARGS_ADVANCE_UNDEF
#endif

#define MED3 NAME_COMPOSE(med3,SORT_TYPE_SUFFIX)
#define VECSWAP NAME_COMPOSE(vecswap,SORT_TYPE_SUFFIX)
#define XT_QUICKSORT NAME_COMPOSE(xt_quicksort,SORT_TYPE_SUFFIX)
#define XT_QUICKSORT_INNER NAME_COMPOSE(xt_qsort_i,SORT_TYPE_SUFFIX)

static inline size_t
MED3(const SORT_TYPE *a, size_t i, size_t j, size_t k
     XT_SORT_MED3_EXTRA_ARGS_DECL)
{
  return SORT_TYPE_CMP_LT(a[i],a[j],i,j) ?
    (SORT_TYPE_CMP_LT(a[j],a[k],j,k) ?
     j : (SORT_TYPE_CMP_LT(a[i],a[k],i,k) ? k : i ))
    : (SORT_TYPE_CMP_LT(a[k],a[j],k,j) ?
       j : (SORT_TYPE_CMP_LT(a[i],a[k],i,k) ? i : k));
}

#ifndef SWAP
#define SWAP(i,j) do {                                \
    SORT_TYPE t = a[i]; a[i] = a[j]; a[j] = t;        \
    XT_SORT_EXTRA_ARGS_SWAP(i, j);                    \
  } while (0)
#endif

static inline void
VECSWAP(SORT_TYPE *restrict a, size_t ia, size_t ib, size_t n XT_SORT_VECSWAP_EXTRA_ARGS_DECL)
{
  for (size_t i = 0; i < n; ++i)
    SWAP(ia+i, ib+i);
}

static
void XT_QUICKSORT_INNER(SORT_TYPE *restrict a, size_t n XT_SORT_EXTRA_ARGS_DECL)
{
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

  while (1) {
    bool swap_cnt = false;
    if (n < 7) {
      for (size_t m = 1; m < n; ++m)
        for (size_t l = m; l > 0 && SORT_TYPE_CMP_LT(a[l], a[l-1],l,l-1); --l)
          SWAP(l, l-1);
      return;
    }
    {
      size_t m = n/2;
      if (n > 7) {
        size_t l = 0, k = n - 1;
        if (n > 40) {
          size_t d = n / 8;
          l = MED3(a, l,          l + d,  l + 2 * d
                   XT_SORT_MED3_EXTRA_ARGS_PASS);
          m = MED3(a, m - d,      m,      m + d
                   XT_SORT_MED3_EXTRA_ARGS_PASS);
          k = MED3(a, k - 2 * d,  k - d, k
                   XT_SORT_MED3_EXTRA_ARGS_PASS);
        }
        m = MED3(a, l,  m,  k XT_SORT_MED3_EXTRA_ARGS_PASS);
      }
      SWAP(0, m);
    }
    size_t i = 1, b = i;
    size_t c = n - 1, d = c;
    SORT_TYPE pivot = *a;
    for (;;) {
      while (b <= c && SORT_TYPE_CMP_LE(a[b], pivot, b, 0)) {
        if (SORT_TYPE_CMP_EQ(a[b], pivot, b, 0)) {
          swap_cnt = true;
          SWAP(i, b);
          ++i;
        }
        ++b;
      }
      while (b <= c && SORT_TYPE_CMP_LE(pivot, a[c], 0, c)) {
        if (SORT_TYPE_CMP_EQ(pivot, a[c], 0, c)) {
          swap_cnt = true;
          SWAP(c, d);
          --d;
        }
        --c;
      }
      if (b > c)
        break;
      SWAP(b, c);
      swap_cnt = true;
      ++b;
      --c;
    }
    if (n < (16384 / sizeof (SORT_TYPE)) && !swap_cnt) {  /* Switch to insertion sort */
      for (size_t m = 1; m < n; ++m)
        for (size_t l = m; l > 0 && SORT_TYPE_CMP_LT(a[l], a[l-1], l, l-1); --l)
          SWAP(l, l-1);
      return;
    }

    size_t pdiff = MIN(i, b - i);
    VECSWAP(a, 0, b - pdiff, pdiff XT_SORT_VECSWAP_EXTRA_ARGS_PASS);
    pdiff = MIN(d - c, n - d - 1);
    VECSWAP(a, b, n - pdiff, pdiff XT_SORT_VECSWAP_EXTRA_ARGS_PASS);
    if ((pdiff = b - i) > 1U)
      XT_QUICKSORT_INNER(a, pdiff XT_SORT_EXTRA_ARGS_PASS);
    if ((pdiff = d - c) > 1U) {
      /* Iterate rather than recurse to save stack space */
      size_t adv = n - pdiff;
      a += adv;
      n = pdiff;
      XT_SORT_EXTRA_ARGS_ADVANCE(adv);
    }
    else
      break;
  }
#undef MIN
}

XT_SORTFUNC_DECL
void XT_QUICKSORT(SORT_TYPE *restrict a, size_t n XT_SORT_EXTRA_ARGS_DECL)
{
  XT_QUICKSORT_INNER(a, n XT_SORT_EXTRA_ARGS_PASS);
}


#ifdef XT_SORT_MED3_EXTRA_ARGS_PASS_UNDEF
#undef XT_SORT_MED3_EXTRA_ARGS_PASS
#undef XT_SORT_MED3_EXTRA_ARGS_PASS_UNDEF
#endif

#ifdef XT_SORT_MED3_EXTRA_ARGS_DECL_UNDEF
#undef XT_SORT_MED3_EXTRA_ARGS_DECL
#undef XT_SORT_MED3_EXTRA_ARGS_DECL_UNDEF
#endif

#ifdef XT_SORT_VECSWAP_EXTRA_ARGS_PASS_UNDEF
#undef XT_SORT_VECSWAP_EXTRA_ARGS_PASS
#undef XT_SORT_VECSWAP_EXTRA_ARGS_PASS_UNDEF
#endif

#ifdef XT_SORT_VECSWAP_EXTRA_ARGS_DECL_UNDEF
#undef XT_SORT_VECSWAP_EXTRA_ARGS_DECL
#undef XT_SORT_VECSWAP_EXTRA_ARGS_DECL_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_DECL_UNDEF
#undef XT_SORT_EXTRA_ARGS_DECL
#undef XT_SORT_EXTRA_ARGS_DECL_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_DECL_UNDEF
#undef XT_SORT_EXTRA_ARGS_DECL
#undef XT_SORT_EXTRA_ARGS_DECL_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_PASS_UNDEF
#undef XT_SORT_EXTRA_ARGS_PASS
#undef XT_SORT_EXTRA_ARGS_PASS_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_SWAP_UNDEF
#undef XT_SORT_EXTRA_ARGS_SWAP
#undef XT_SORT_EXTRA_ARGS_SWAP_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_ADVANCE_UNDEF
#undef XT_SORT_EXTRA_ARGS_ADVANCE
#undef XT_SORT_EXTRA_ARGS_ADVANCE_UNDEF
#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
