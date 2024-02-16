/*
 * @file bsearch.h
 * @brief generic routines for re-entrant binary search
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * Keywords: binary search
 * @author Thomas Jahns <jahns@dkrz.de>
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

#ifndef PPM_BSEARCH_H
#define PPM_BSEARCH_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include "core/fptr_api.h"

/**
 * @brief Search sorted array a for key.
 *
 * The function requires a to in ascending sort according to cmp.
 * PPM_bsearch_r differs from bsearch(3) in that the comparision
 * function accepts an optional data argument.
 * @param[in] key pointer to element to search for
 * @param[in] a pointer to base address of array
 * @param nmemb number of elements in array
 * @param es size of element
 * @param[in] data pointer to arbitrary data used by cmp
 * @param[in] cmp comparison function
 * @return pointer to element of a matching key or NULL if not found
 */
void PPM_DSO_API_EXPORT*
PPM_bsearch_r(const void *key, const void *a, size_t nmemb, size_t size,
              void *data, PPM_CompareWithData cmp);

/**
 * @brief Search sorted array a for key or next element of lower sort.
 *
 * The function requires a to in ascending sort according to cmp. It
 * differs from PPM_bsearch_r only in the returned element in case key
 * is not found in a.
 * @param[in] key pointer to element to search for
 * @param[in] a pointer to base address of array
 * @param nmemb number of elements in array
 * @param es size of element
 * @param[in] data pointer to arbitrary data used by cmp
 * @param[in] cmp comparison function
 * @return pointer to element of a matching key or greatest element
 * less than key or NULL if cmp(key, a, data) < 0
 */
void PPM_DSO_API_EXPORT*
PPM_bsearch_el_r(const void *key, const void *a, size_t nmemb, size_t size,
                 void *data, PPM_CompareWithData cmp);


/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */

#endif
