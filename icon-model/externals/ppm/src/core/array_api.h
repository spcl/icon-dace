/**
 * @file array_api.h
 * @brief genometools array class adapted for ScalES-PPM
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
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

#ifndef ARRAY_API_H
#define ARRAY_API_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdlib.h>
#include "core/fptr_api.h"

/* <PPM_Array> objects are generic arrays for elements of a certain size which
   grow on demand. */
typedef struct PPM_Array PPM_Array;

/* Return a new <PPM_Array> object whose elements have the size <size_of_elem>. */
PPM_Array*      PPM_array_new(size_t size_of_elem);
/* Increase the reference count for <array> and return it.
   If <array> is <NULL>, <NULL> is returned without any side effects. */
PPM_Array*      PPM_array_ref(PPM_Array *array);
/* Return a clone of <array>. */
PPM_Array*      PPM_array_clone(const PPM_Array *array);
/* Return pointer to element number <index> of <array>. <index> has to be
   smaller than <PPM_array_size(array)>. */
void*         PPM_array_get(const PPM_Array *array, unsigned long index);
/* Return pointer to first element of <array>. */
void*         PPM_array_get_first(const PPM_Array *array);
/* Return pointer to last element of <array>. */
void*         PPM_array_get_last(const PPM_Array *array);
/* Return pointer to last element of <array> and remove it from <array>. */
void*         PPM_array_pop(PPM_Array *array);
/* Return pointer to the internal space of <array> where the elements are
   stored.  */
void*         PPM_array_get_space(const PPM_Array *array);
/* Add element <elem> to <array>. The size of <elem> must equal the given
   element size when the <array> was created and is determined automatically
   with the <sizeof> operator. */
#define       PPM_array_add(array, elem) \
              PPM_array_add_elem(array, &(elem), sizeof (elem))
/* Handle most frequent use case: add pointer to array of references */
void PPM_array_add_ptr(PPM_Array *a, void *elem);
/* Add element <elem> with size <size_of_elem> to <array>. <size_of_elem> must
   equal the given element size when the <array> was created. Usually, this
   method is not used directly and the macro <PPM_array_add()> is used
   instead. */
void          PPM_array_add_elem(PPM_Array *array, void *elem,
                                size_t size_of_elem);
/* Add all elements of array <src> to the array <dest>. The element sizes of
   both arrays must be equal. */
void          PPM_array_add_array(PPM_Array *dest, const PPM_Array *src);
/* Remove element with number <index> from <array> in O(<PPM_array_size(array)>)
   time. <index> has to be smaller than <PPM_array_size(array)>. */
void          PPM_array_rem(PPM_Array *array, unsigned long index);
/* Remove elements starting with number <frompos> up to (and including) <topos>
   from <array> in O(<PPM_array_size(array)>) time. <frompos> has to be smaller
   or equal than <topos> and both have to be smaller than
   <PPM_array_size(array)>. */
void          PPM_array_rem_span(PPM_Array *array, unsigned long frompos,
                                unsigned long topos);
/* Reverse the order of the elements in <array>. */
void          PPM_array_reverse(PPM_Array *array);
/* Set the size of <array> to <size>. <size> must be smaller or equal than
   <PPM_array_size(array)>. */
void          PPM_array_set_size(PPM_Array *array, unsigned long size);
/* Reset the <array>. That is, afterwards the array has size 0. */
void          PPM_array_reset(PPM_Array *array);
/* Return the size of the elements stored in <array>. */
size_t        PPM_array_elem_size(const PPM_Array *array);
/* Return the number of elements in <array>. If <array> equals <NULL>, 0 is
   returned. */
unsigned long PPM_array_size(const PPM_Array *array);
/* Sort <array> with the given compare function <compar>. */
void          PPM_array_sort(PPM_Array *array, PPM_Compare compar);
/* Sort <array> in a stable way with the given compare function <compar>. */
void          PPM_array_sort_stable(PPM_Array *array, PPM_Compare compar);
/* Sort <array> with the given compare function <compar>. Passes a pointer with
   userdata <data> to <compar>. */
void          PPM_array_sort_with_data(PPM_Array *array,
                                      PPM_CompareWithData compar,
                                      void *data);
/* Sort <array> in a stable way with the given compare function <compar>. Passes
   a pointer with userdata <data> to <compar>. */
void          PPM_array_sort_stable_with_data(PPM_Array *array,
                                             PPM_CompareWithData compar,
                                             void *data);
/* Compare the content of <array_a> with the content of <array_b>.
   <array_a> and <array_b> must have the same <PPM_array_size()> and
   <PPM_array_elem_size()>. */
int           PPM_array_cmp(const PPM_Array *array_a, const PPM_Array *array_b);
/* Decrease the reference count for <array> or delete it, if this was the last
   reference. */
void          PPM_array_delete(PPM_Array *array);

#endif
/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
