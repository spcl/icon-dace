/**
 * @file heap.c
 * @brief routines for sorting arrays into binary heap order
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
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
 *
 */
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>

#include "core/heap.h"
#include "core/ppm_visibility.h"
#include "core/fptr_api.h"
#include "core/swapmacros.h"

static inline size_t
left(size_t i)
{
  return 2*i + 1;
}

static inline size_t
right(size_t i)
{
  return 2*i + 2;
}

static inline size_t
parent(size_t i)
{
  return (i - 1) / 2;
}


#define ep(i) ((unsigned char *)heap + (i) * es)

void PPM_DSO_API_EXPORT
PPM_heapify(void *heap, size_t n, size_t es, size_t i, void *data,
            PPM_CompareWithData cmp)
{
  size_t l, r, largest;
  int swaptype;
  assert(i < n && es > 0);
  SWAPINIT(heap, es);
  do {
    l = left(i);
    r = right(i);
    if (l < n && cmp(ep(l), ep(i), data) > 0)
    {
      largest = l;
    }
    else
      largest = i;
    if (r < n && cmp(ep(r), ep(largest), data) > 0)
      largest = r;
    if (largest == i) break;
    swap(ep(i), ep(largest));
    i = largest;
  } while(1);
}

int PPM_DSO_API_EXPORT
PPM_is_heap(void *heap, size_t n, size_t es, void *data,
            PPM_CompareWithData cmp)
{
  size_t i = n;
  assert(es > 0);
  while (i-- > 1)
    if (cmp(ep(parent(i)), ep(i), data) < 0)
      return 0;
  return 1;
}

void PPM_DSO_API_EXPORT
PPM_build_heap(void *heap, size_t n, size_t es, void *data,
               PPM_CompareWithData cmp)
{
  size_t i = n/2;
  assert(es > 0);
  while (i--)
    PPM_heapify(heap, n, es, i, data, cmp);
}

void PPM_DSO_API_EXPORT
PPM_heap_remove_top(void *heap, size_t n, size_t es, void *data,
                    PPM_CompareWithData cmp)
{
  int swaptype;
  assert(n > 0 && es > 0);
  SWAPINIT(heap, es);
  swap(ep(0), ep(n - 1));
  PPM_heapify(heap, --n, es, 0, data, cmp);
}

void PPM_DSO_API_EXPORT
PPM_heap_elem_increase_sort(void *heap, size_t n, size_t es, size_t i,
                            void *data, PPM_CompareWithData cmp)
{
  int swaptype;
  void *pp, *ip;
  size_t p;
  (void)n;
  assert(es > 0 && i < n);
  SWAPINIT(heap, es);
  while (i && cmp(pp = ep(p = parent(i)), ip = ep(i), data) < 0)
  {
    swap(pp, ip);
    i = p;
  }
}

void PPM_DSO_API_EXPORT
PPM_heap_leaf_minimize(void *heap, size_t n, size_t es,
                       void *data, PPM_CompareWithData cmp)
{
  int swaptype;
  void *last, *p, *minp;
  size_t i, min = 0;
  assert(es > 0);
  SWAPINIT(heap, es);
  last = minp = ep(n - 1);
  if (n < 3)
    return;
  for (i = n - 2; i >= n/2; --i)
    if (cmp(p = ep(i), minp, data) < 0)
    {
      min = i;
      minp = p;
    }
  if (minp != last)
  {
    swap(last, minp);
    PPM_heap_elem_increase_sort(heap, n, es, min, data, cmp);
  }
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
