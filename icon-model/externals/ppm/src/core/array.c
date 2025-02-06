/**
 * @file array.c
 * @brief genometools array class adapted for ScalES-PPM
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
 */
/*
  Copyright (c) 2005-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <assert.h>
#include <limits.h>
#include <string.h>
#include "core/ppm_xfuncs.h"
#include "core/array.h"
#include "core/msort.h"
#include "core/qsort_r.h"
#include "core/ppm_extents.h"

struct PPM_Array {
  void *space;
  unsigned long next_free;
  size_t allocated,
         size_of_elem;
  unsigned int reference_count;
};

PPM_Array* PPM_array_new(size_t size_of_elem)
{
  PPM_Array *a = xcalloc(1, sizeof *a);
  assert(size_of_elem);
  a->size_of_elem = size_of_elem;

  return a;
}

PPM_Array* PPM_array_clone(const PPM_Array *a)
{
  PPM_Array *a_copy;
  assert(a);
  a_copy = xmalloc(sizeof (PPM_Array));
  a_copy->space = xmalloc(a->next_free * a->size_of_elem);
  memcpy(a_copy->space, a->space, a->next_free * a->size_of_elem);
  a_copy->next_free = a_copy->allocated = a->next_free;
  a_copy->size_of_elem = a->size_of_elem;
  a_copy->reference_count = 0;
  return a_copy;
}

PPM_Array* PPM_array_ref(PPM_Array *a)
{
  if (!a) return NULL;
  a->reference_count++;
  return a;
}

void* PPM_array_get(const PPM_Array *a, unsigned long idx)
{
  assert(a && idx < a->next_free);
  return (char*) a->space + idx * a->size_of_elem;
}

void* PPM_array_get_first(const PPM_Array *a)
{
  return PPM_array_get(a, 0);
}

void* PPM_array_get_last(const PPM_Array *a)
{
  assert(a->next_free);
  return PPM_array_get(a, a->next_free-1);
}

void* PPM_array_pop(PPM_Array *a)
{
  assert(a && a->next_free);
  a->next_free--;
  return (char*) a->space + a->next_free * a->size_of_elem;
}

void PPM_array_rem(PPM_Array *a, unsigned long idx)
{
  unsigned long i;
  assert(a && idx < a->next_free);
  /* move elements */
  for (i = idx + 1; i < a->next_free; i++) {
    memcpy((char*) a->space + (i-1) * a->size_of_elem,
           (char*) a->space + i * a->size_of_elem,
           a->size_of_elem);
  }
  /* remove last (now duplicated) element */
  a->next_free--;
}

void PPM_array_rem_span(PPM_Array *a, unsigned long frompos, unsigned long topos)
{
  unsigned long i, len;
  assert(a && frompos <= topos);
  assert(frompos < a->next_free && topos < a->next_free);
  /* move elements */
  len = topos - frompos + 1;
  for (i = topos + 1; i < a->next_free; i++) {
    memcpy((char*) a->space + (i-len) * a->size_of_elem,
           (char*) a->space + i * a->size_of_elem,
           a->size_of_elem);
  }
  /* remove last (now duplicated) elements */
  a->next_free -= len;
}

void PPM_array_reverse(PPM_Array *a)
{
  char *front, *back, *tmp;
  assert(a);
  tmp = xmalloc(a->size_of_elem);
  for (front = a->space,
       back = (char*) a->space + (a->next_free-1) * a->size_of_elem;
       front < back; front += a->size_of_elem, back -= a->size_of_elem) {
    memcpy(tmp, front, a->size_of_elem);
    memcpy(front, back, a->size_of_elem);
    memcpy(back, tmp, a->size_of_elem);
  }
  free(tmp);
}

void* PPM_array_get_space(const PPM_Array *a)
{
  assert(a);
  return a->space;
}

void PPM_array_add_ptr(PPM_Array *a, void *elem)
{
  PPM_array_add(a, elem);
}

void PPM_array_add_elem(PPM_Array *a, void *elem, size_t size_of_elem)
{
  assert(a && elem);
  assert(a->size_of_elem == size_of_elem);
  assert(a->next_free <= a->allocated);
  /* make sure we have enough space */
  if ((a->next_free + 1) * size_of_elem > a->allocated) {
    a->space = xrealloc(a->space, (a->next_free + 1) * size_of_elem);
    a->allocated = (a->next_free + 1) * size_of_elem;
  }
  /* add */
  memcpy((char*) a->space + a->next_free * size_of_elem, elem, size_of_elem);
  a->next_free++;
}

void PPM_array_add_array(PPM_Array *dest, const PPM_Array *src)
{
  unsigned long i;
  assert(dest && src && dest->size_of_elem == src->size_of_elem);
  for (i = 0; i < PPM_array_size(src); i++)
    PPM_array_add_elem(dest, PPM_array_get(src, i), src->size_of_elem);
}

size_t PPM_array_elem_size(const PPM_Array *a)
{
  assert(a);
  return a->size_of_elem;
}

unsigned long PPM_array_size(const PPM_Array *a)
{
  return a ? a->next_free : 0;
}

void PPM_array_set_size(PPM_Array *a, unsigned long size)
{
  assert(a);
  assert(size <= a->next_free);
  a->next_free = size;
}

void PPM_array_reset(PPM_Array *a)
{
  assert(a);
  a->next_free = 0;
}

void PPM_array_sort(PPM_Array *a, PPM_Compare compar)
{
  assert(a && compar);
  qsort(a->space, a->next_free, a->size_of_elem, compar);
}

void PPM_array_sort_stable(PPM_Array *a, PPM_Compare compar)
{
  assert(a && compar);
  PPM_msort(a->space, a->next_free, a->size_of_elem, compar);
}

void PPM_array_sort_with_data(PPM_Array *a, PPM_CompareWithData compar, void *data)
{
  assert(a && compar);
  PPM_qsort_r(a->space, a->next_free, a->size_of_elem, data, compar);
}

void PPM_array_sort_stable_with_data(PPM_Array *a, PPM_CompareWithData compar,
                                    void *data)
{
  assert(a && compar);
  PPM_msort_r(a->space, a->next_free, a->size_of_elem, data, compar);
}

int PPM_array_cmp(const PPM_Array *array_a, const PPM_Array *array_b)
{
  assert(PPM_array_size(array_a) == PPM_array_size(array_b));
  assert(PPM_array_elem_size(array_a) == PPM_array_elem_size(array_b));
  return memcmp(array_a->space, array_b->space,
                array_a->size_of_elem * array_a->next_free);
}

bool PPM_array_equal(const PPM_Array *a, const PPM_Array *b, PPM_Compare cmpfunc)
{
  unsigned long idx, size_a, size_b;
  int cmp;
  assert(PPM_array_elem_size(a) == PPM_array_elem_size(b));
  size_a = PPM_array_size(a);
  size_b = PPM_array_size(b);
  if (size_a < size_b)
    return false;
  if (size_a > size_b)
    return false;
  for (idx = 0; idx < size_a; idx++) {
    cmp = cmpfunc(PPM_array_get(a, idx), PPM_array_get(b, idx));
    if (cmp != 0)
      return false;
  }
  return true;
}

int PPM_array_iterate(PPM_Array *a, PPM_ArrayProcessor array_processor, void *info)
{
  unsigned long idx;
  int rval;
  assert(a && array_processor);
  for (idx = 0; idx < PPM_array_size(a); idx++) {
    if ((rval = array_processor(PPM_array_get(a, idx), info)))
      return rval;
  }
  return 0;
}

int PPM_array_iterate_reverse(PPM_Array *a, PPM_ArrayProcessor array_processor,
                             void *info)
{
  unsigned long idx;
  int rval;
  assert(a && array_processor);
  for (idx = PPM_array_size(a); idx > 0; idx--) {
    if ((rval = array_processor(PPM_array_get(a, idx-1), info)))
      return rval;
  }
  return 0;
}

void PPM_array_prepend_array(PPM_Array *dest, const PPM_Array *src)
{
  unsigned long i;
  assert(dest && src && dest->size_of_elem == src->size_of_elem);
  if (!src->next_free)
    return; /* nothing to do */
  /* make sure <dest> is large enough */
  dest->space = xrealloc(dest->space,
                         (dest->next_free + src->next_free) *
                         dest->size_of_elem);
  dest->allocated = (dest->next_free + src->next_free) * dest->size_of_elem;
  /* move elements in <dest> to the back */
  for (i = dest->next_free; i > 0; i--) {
    memcpy((char*) dest->space + (i-1+src->next_free) * dest->size_of_elem,
           (char*) dest->space + (i-1) * dest->size_of_elem,
           dest->size_of_elem);
  }
  /* copy <src> to the start of <dest> */
  memcpy((char*) dest->space, (char*) src->space,
         src->size_of_elem * src->next_free);
  dest->next_free += src->next_free;
}

void PPM_array_delete(PPM_Array *a)
{
  if (!a) return;
  if (a->reference_count) {
    a->reference_count--;
    return;
  }
  free(a->space);
  free(a);
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
