/**
 * @file heap.h
 * @brief routines for sorting arrays into binary heap order
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

#ifndef PPM_HEAP_BASE_H
#define PPM_HEAP_BASE_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include "core/fptr_api.h"

/**
 * @brief restore binary heap property for sub-tree starting at i
 * This function assumes that the heap property is fulfilled for all
 * child trees of i.
 * @param heap pointer to base address of array
 * @param es size of element
 * @param n number of elements in array
 * @param i start element (multiplied with es to compute offset)
 * @param data pointer to arbitrary data used by cmp
 * @param cmp comparison function
 */
void
PPM_heapify(void *heap, size_t n, size_t es, size_t i, void *data,
            PPM_CompareWithData cmp);

/**
 * @brief restore binary heap property for sub-tree starting at i
 * @param heap pointer to base address of array
 * @param es size of element
 * @param n number of elements in array
 * @param data pointer to arbitrary data used by cmp
 * @param cmp comparison function
 */
void
PPM_build_heap(void *heap, size_t n, size_t es, void *data,
               PPM_CompareWithData cmp);

/**
 * @brief check whether array fulfills binary heap property
 * @param heap pointer to base address of array
 * @param es size of element
 * @param n number of elements in array
 * @param data pointer to arbitrary data used by cmp
 * @param cmp comparison function
 * @return 0 if heap does not fulfill heap property, non-zero otherwise
 */
int
PPM_is_heap(void *heap, size_t n, size_t es, void *data,
            PPM_CompareWithData cmp);


/**
 * @brief remove heap top by replacing with leaf and restore heap property
 * @param heap pointer to base address of array
 * @param es size of element
 * @param n number of elements in array
 * @param data pointer to arbitrary data used by cmp
 * @param cmp comparison function
 */
void
PPM_heap_remove_top(void *heap, size_t n, size_t es, void *data,
                    PPM_CompareWithData cmp);

/**
 * @brief restore heap property after increasing sort of single
 * element i of heap
 * @param heap pointer to base address of array
 * @param es size of element
 * @param n number of elements in array
 * @param i modified element (multiplied with es to compute offset)
 * @param data pointer to arbitrary data used by cmp
 * @param cmp comparison function
 */
void
PPM_heap_elem_increase_sort(void *heap, size_t n, size_t es, size_t i,
                            void *data, PPM_CompareWithData cmp);



/**
 * @brief rearrange leafs such that last element is also lowest sorting
 * @param heap pointer to base address of array
 * @param es size of element
 * @param n number of elements in array
 * @param data pointer to arbitrary data used by cmp
 * @param cmp comparison function
 */
void
PPM_heap_leaf_minimize(void *heap, size_t n, size_t es,
                       void *data, PPM_CompareWithData cmp);


/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */

#endif
