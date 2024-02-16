/**
 * @file test_idxlist_utils.h
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

#ifndef TEST_IDXLIST_UTILS_H
#define TEST_IDXLIST_UTILS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>
#include <stdlib.h>

#include <yaxt.h>

void check_idxlist(Xt_idxlist idxlist, Xt_int const * ref_indices,
                   int ref_num_indices);

void
check_stripes(struct Xt_stripe const * stripes, int num_stripes,
              struct Xt_stripe const * ref_stripes, int ref_num_stripes);

void
check_offsets(size_t num_offsets, const int offsets_a[num_offsets],
              const int offsets_b[num_offsets]);

Xt_idxlist
idxlist_pack_unpack_copy(Xt_idxlist idxlist);

void
check_idxlist_copy(Xt_idxlist idxlist, Xt_idxlist idxlist_copy,
                   size_t num_ref_indices, const Xt_int *ref_indices,
                   size_t num_ref_stripes, const struct Xt_stripe *ref_stripes);

static inline int
cmp_idx_arrays(size_t num_indices,
               const Xt_int * restrict a, const Xt_int * restrict b)
{
  bool diff = false;
  for (size_t i = 0; i < num_indices; ++i)
    diff |= a[i] != b[i];
  return diff;
}

#endif // TEST_IDXLIST_UTILS_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
