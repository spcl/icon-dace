/**
 * @file ensure_array_size.h
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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

#ifndef ENSURE_ARRAY_SIZE_H
#define ENSURE_ARRAY_SIZE_H

#include <stdlib.h>

void
yac_realloc_array(void **array, size_t elem_size, size_t *curr_array_size,
                  size_t requested_size);

#define ENSURE_ARRAY_SIZE(arrayp, curr_array_size, req_size)            \
  do {                                                                  \
    if ((req_size) > (curr_array_size))                                 \
    {                                                                   \
      size_t casize = (curr_array_size);                                \
                                                                        \
      yac_realloc_array((void **)&(arrayp), sizeof(*(arrayp)), &casize, \
                        (req_size));                                    \
      (curr_array_size) = casize;                                       \
    }                                                                   \
  }                                                                     \
  while(0)

#define ENSURE_BYTE_ARRAY_SIZE(arrayp, curr_array_size, req_size) \
  do {                                                            \
    if ((req_size) > (curr_array_size))                           \
    {                                                             \
       size_t casize = (curr_array_size);                         \
                                                                  \
       yac_realloc_array(&(arrayp), 1, &casize, (req_size));      \
                                                                  \
       (curr_array_size) = casize;                                \
    }                                                             \
  }                                                               \
  while(0)

#endif

