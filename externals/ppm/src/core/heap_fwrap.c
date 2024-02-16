/**
 * @file heap_fwrap.c
 * @brief Fortran wrappers for binary heap routines
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords: binary heap
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

#include <stdlib.h>

#include "cfortran.h"

#include "core/heap.h"

static void
PPM_heapify_wrapper(void *a, int n, int es, int i,
                    void *data, PPM_CompareWithData cmp)
{
  PPM_heapify(a, (size_t)n, (size_t)es, (size_t)i - 1, data, cmp);
}

#undef ROUTINE_6
#define ROUTINE_6 (PPM_CompareWithData)(void (*)(void))
FCALLSCSUB6(PPM_heapify_wrapper, PPM_HEAPIFY, ppm_heapify,
            PVOID, INT, INT, INT, PVOID, ROUTINE)


static int
PPM_is_heap_wrapper(void *a, int n, int es,
                    void *data, PPM_CompareWithData cmp)
{
  return PPM_is_heap(a, (size_t)n, (size_t)es, data, cmp);
}


#undef ROUTINE_5
#define ROUTINE_5 (PPM_CompareWithData)(void (*)(void))
FCALLSCFUN5(LOGICAL,PPM_is_heap_wrapper, PPM_IS_HEAP, ppm_is_heap,
            PVOID, INT, INT, PVOID, ROUTINE)

static void
PPM_build_heap_wrapper(void *a, int n, int es,
                       void *data, PPM_CompareWithData cmp)
{
  PPM_build_heap(a, (size_t)n, (size_t)es, data, cmp);
}

FCALLSCSUB5(PPM_build_heap_wrapper, PPM_BUILD_HEAP, ppm_build_heap,
            PVOID, INT, INT, PVOID, ROUTINE)


static void
PPM_heap_elem_increase_sort_wrapper(void *heap, int n, int es,
                                    int i, void *data, PPM_CompareWithData cmp)
{
  PPM_heap_elem_increase_sort(heap, (size_t)n, (size_t)es, (size_t)i-1, data,
                              cmp);
}

FCALLSCSUB6(PPM_heap_elem_increase_sort_wrapper, PPM_HEAP_ELEM_INCREASE_SORT,
            ppm_heap_elem_increase_sort,
            PVOID, INT, INT, INT, PVOID, ROUTINE)

static void
PPM_heap_remove_top_wrapper(void *heap, int n, int es, void *data,
                            PPM_CompareWithData cmp)
{
  PPM_heap_remove_top(heap, (size_t)n, (size_t)es, data,cmp);
}

FCALLSCSUB5(PPM_heap_remove_top_wrapper, PPM_HEAP_REMOVE_TOP,
            ppm_heap_remove_top,
            PVOID, INT, INT, PVOID, ROUTINE)

static void
PPM_heap_leaf_minimize_wrapper(void *a, int n, int es,
                               void *data, PPM_CompareWithData cmp)
{
  PPM_heap_leaf_minimize(a, (size_t)n, (size_t)es, data, cmp);
}

FCALLSCSUB5(PPM_heap_leaf_minimize_wrapper, PPM_HEAP_LEAF_MINIMIZE,
            ppm_heap_leaf_minimize,
            PVOID, INT, INT, PVOID, ROUTINE)

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
