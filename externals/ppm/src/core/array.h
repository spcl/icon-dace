/**
 * @file array.h
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

#ifndef ARRAY_H
#define ARRAY_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdbool.h>

#include "core/array_api.h"

typedef int (*PPM_ArrayProcessor)(void *elem, void *info);

/* Compare the raw content of <array_a> with the content of <array_b>.
   <array_a> and <array_b> must have the same PPM_array_size() and
   PPM_array_elem_size(). */
int           PPM_array_cmp(const PPM_Array  *array_a, const PPM_Array *array_b);
/* Compare the content of <array_a> with the content of <array_b> using the
   comparator function <cmpfunc>. If the elements of both arrays are equal
   w.r.t. <cmpfunc>, true is returned. If the array sizes or content w.r.t.
   <cmpfunc> are different, false is returned. */
bool          PPM_array_equal(const PPM_Array *a, const PPM_Array *b,
                             PPM_Compare cmpfunc);
/* Iterate over all elements in <array> and call <array_processor> with them.
   <info> and <err> are passed to <array_processor>.
   If <array_processor> returns a value != 0, the iteration is stopped and the
   return value of <array_processor> is returned. */
int           PPM_array_iterate(PPM_Array *array,
                               PPM_ArrayProcessor array_processor,
                               void *info);
/* Similar to <array_iterate>, except that the <array> is traversed in reverse
   order. */
int           PPM_array_iterate_reverse(PPM_Array *array,
                                       PPM_ArrayProcessor array_processor,
                                       void *info);
void          PPM_array_prepend_array(PPM_Array *dest, const PPM_Array *src);

#endif

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
