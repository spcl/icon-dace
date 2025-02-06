/**
 * @file qsort_r.c
 * @brief Functions for generic quick sort with extra parameter for
 * improved re-entrancy
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
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
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <inttypes.h>
#include <stddef.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/qsort_r.h"
#include "core/fptr_api.h"
#include "core/minmax.h"
#include "core/swapmacros.h"

static inline void *
med3(void *a, void *b, void *c, PPM_CompareWithData cmp, void *data);

#define vecswap(a, b, n)       if ((n) > 0) swapfunc(a, b, n, swaptype)

static inline void *
med3(void *a, void *b, void *c, PPM_CompareWithData cmp, void *data)
{
  return cmp(a, b, data) < 0 ?
    (cmp(b, c, data) < 0 ? b : (cmp(a, c, data) < 0 ? c : a ))
    :(cmp(b, c, data) > 0 ? b : (cmp(a, c, data) < 0 ? a : c ));
}

void
PPM_qsort_r(void *a, size_t n, size_t es, void *data, PPM_CompareWithData cmp)
{
  unsigned char *pa, *pb, *pc, *pd, *pl, *pm, *pn;
  int r, swaptype, swap_cnt;

  SWAPINIT(a, es);
  while (1)
  {
    swap_cnt = 0;
    if (n < 7) {
      for (pm = (unsigned char *)a + es; pm < (unsigned char *)a + n*es; pm+=es)
        for (pl = pm;
             pl > (unsigned char *)a && cmp(pl - es, pl, data) > 0;
             pl -= es)
          swap(pl, pl - es);
      return;
    }
    pm = (unsigned char *)a + (n / 2) * es;
    if (n > 7) {
      pl = a;
      pn = (unsigned char *)a + (n - 1) * es;
      if (n > 40) {
        size_t d = (n / 8) * es;
        pl = med3(pl, pl + d, pl + 2 * d, cmp, data);
        pm = med3(pm - d, pm, pm + d, cmp, data);
        pn = med3(pn - 2 * d, pn - d, pn, cmp, data);
      }
      pm = med3(pl, pm, pn, cmp, data);
    }
    swap(a, pm);
    pa = pb = (unsigned char *)a + es;

    pc = pd = (unsigned char *)a + (n - 1) * es;
    for (;;) {
      while (pb <= pc && (r = cmp(pb, a, data)) <= 0) {
        if (r == 0) {
          swap_cnt = 1;
          swap(pa, pb);
          pa += es;
        }
        pb += es;
      }
      while (pb <= pc && (r = cmp(pc, a, data)) >= 0) {
        if (r == 0) {
          swap_cnt = 1;
          swap(pc, pd);
          pd -= es;
        }
        pc -= es;
      }
      if (pb > pc)
        break;
      swap(pb, pc);
      swap_cnt = 1;
      pb += es;
      pc -= es;
    }
    if (swap_cnt == 0) {  /* Switch to insertion sort */
      for (pm = (unsigned char *)a + es; pm < (unsigned char *)a + n*es; pm+=es)
        for (pl = pm;
             pl > (unsigned char *)a && cmp(pl - es, pl, data) > 0;
             pl -= es)
          swap(pl, pl - es);
      return;
    }

    pn = (unsigned char *)a + n * es;
    ptrdiff_t pdiff = MIN(pa - (unsigned char *)a, pb - pa);
    vecswap(a, pb - pdiff, (size_t)pdiff);
    pdiff = MIN(pd - pc, pn - pd - (ptrdiff_t)es);
    vecswap(pb, pn - pdiff, (size_t)pdiff);
    if ((pdiff = pb - pa) > (ptrdiff_t)es)
      PPM_qsort_r(a, (size_t)pdiff / es, es, data, cmp);
    if ((pdiff = pd - pc) > (ptrdiff_t)es) {
      /* Iterate rather than recurse to save stack space */
      a = pn - pdiff;
      n = (size_t)pdiff / es;
    }
    else
      break;
  }
/*            qsort(pn - r, r / es, es, cmp);*/
}

enum {
  par_low_limit = 128 * 4,      /* choose something expected to be at
                                 * least 4 cache-lines in size to
                                 * prevent false sharing */
};

static size_t
partition_task(void *arr, size_t n, size_t es,
               void *data, PPM_CompareWithData cmp)
{
  size_t k = 1, l = n - 1;
  unsigned char *pk = (unsigned char *)arr + es * k,
    *pl = (unsigned char *)arr + es * l;
  int swaptype;
  SWAPINIT(arr, es);
  while (1) {
    while (k < n && cmp(pk, arr, data) <= 0)
    {
      pk += es;
      ++k;
    }
    while (cmp(pl, arr, data) > 0)
    {
      pl -= es;
      --l;
    }
    while (pk < pl)
    {
      swap(pk, pl);
      while (k < n )
      {
        pk += es;
        ++k;
        if (cmp(pk, arr, data) > 0)
          break;
      }
      while (l > 0)
      {
        pl -= es;
        --l;
        if (cmp(pl, arr, data) < 0)
          break;
      }
    }
    swap(arr, pl);
    return l;
  }
}

static void
PPM_qsort_r_task(void *a, size_t n, size_t es,
                 void *data, PPM_CompareWithData cmp)
{
  if (n * es <= par_low_limit)
    PPM_qsort_r(a, n, es, data, cmp);
  else
  {
    size_t q = partition_task(a, n, es, data, cmp);
    unsigned char *pq = (unsigned char *)a + es * (q + 1);
    if (q > 1)
#pragma omp task firstprivate(a, n, q, es, data, cmp)
      PPM_qsort_r_task(a, q, es, data, cmp);
    if (n - q > 2)
#pragma omp task firstprivate(pq, n, q, es, data, cmp)
      PPM_qsort_r_task(pq, n - q - 1, es, data, cmp);
#pragma omp taskwait
  }
}

void
PPM_qsort_r_mt(void *a, size_t n, size_t es, void *data, PPM_CompareWithData cmp)
{
  /* divide a into numthreads blocks and run normal qsort on that */
  int numthreads;
#ifdef _OPENMP
  numthreads = omp_get_num_threads();
#else
  numthreads = 1;
#endif
  if (numthreads == 1)
  {
    PPM_qsort_r(a, n, es, data, cmp);
  }
  else
  {
    PPM_qsort_r_task(a, n, es, data, cmp);
  }
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
