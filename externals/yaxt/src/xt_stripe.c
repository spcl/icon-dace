/**
 * @file xt_stripe.c
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>
#include <stdlib.h>

#include "xt/xt_core.h"
#include "xt/xt_stripe.h"
#include "xt_stripe_util.h"

#include "xt_stripe_util.h"
#include "core/ppm_xfuncs.h"
#include "instr.h"
#include "ensure_array_size.h"

void xt_convert_indices_to_stripes(const Xt_int *restrict indices,
                                   int num_indices,
                                   struct Xt_stripe **stripes,
                                   int * num_stripes)
{
  *stripes = NULL;
  *num_stripes = 0;
  xt_convert_indices_to_stripes_keep_buf(indices, num_indices,
                                         stripes, num_stripes);
}

void xt_convert_indices_to_stripes_keep_buf(const Xt_int *restrict indices,
                                            int num_indices,
                                            struct Xt_stripe **stripes,
                                            int * num_stripes) {

  INSTR_DEF(instr,"xt_idxstripes_convert_to_stripes")
  INSTR_START(instr);

  struct Xt_stripe *restrict temp_stripes = *stripes;
  size_t temp_stripes_array_size = 0;
  size_t num_temp_stripes = (size_t)*num_stripes * sizeof (*temp_stripes);

  if (num_indices > 0) {
    size_t i = 0;

    while(i < (size_t)num_indices) {
      ++num_temp_stripes;

      ENSURE_ARRAY_SIZE(temp_stripes, temp_stripes_array_size, num_temp_stripes);

      size_t j = 1;

      Xt_int stride = 1;
      if (i + j < (size_t)num_indices) {
        stride = (Xt_int)(indices[i + 1] - indices[i]);
        do {
          ++j;
        } while ((i + j) < (size_t)num_indices
                 && indices[i] + (Xt_int)j * stride == indices[i + j]);
      }
      j-= ((i + j + 1 < (size_t)num_indices)
           && ((indices[i + j - 1] == indices[i + j] - 1)
               & (indices[i + j] == indices[i + j + 1] - 1)));

      temp_stripes[num_temp_stripes-1].start  = indices[i];
      temp_stripes[num_temp_stripes-1].stride = stride;
      temp_stripes[num_temp_stripes-1].nstrides = (int)j;

      i = i + j;
    }
  }

  *stripes = xrealloc(temp_stripes, num_temp_stripes * sizeof(*temp_stripes));
  *num_stripes = (int)num_temp_stripes;
  INSTR_STOP(instr);
}

size_t
xt_stripes_merge_copy(size_t num_stripes,
                      struct Xt_stripe *stripes_dst,
                      const struct Xt_stripe *stripes_src,
                      bool lookback)
{
  size_t skip = 1;
  if (num_stripes) {
    if (lookback) {
      Xt_int stride = stripes_src[0].stride,
        prev_stride = stripes_dst[-1].stride,
        start = stripes_src[0].start,
        prev_start = stripes_dst[-1].start;
      if (stride == prev_stride
          && start == prev_start + stride * (Xt_int)stripes_dst[-1].nstrides) {
        /* merge perfectly aligned stripes */
        stripes_dst[-1].nstrides += stripes_src[0].nstrides;
        ++skip;
        goto copy_loop;
      }
    }
    stripes_dst[0] = stripes_src[0];
  copy_loop:
    if (num_stripes > 1)
      for (size_t i = 1; i < num_stripes; ++i) {
        Xt_int stride = stripes_src[i].stride,
          prev_stride = stripes_dst[i-skip].stride,
          start = stripes_src[i].start,
          prev_start = stripes_dst[i-skip].start;
        if (stride == prev_stride
            && start == prev_start + stride * (Xt_int)stripes_dst[i-skip].nstrides) {
          /* merge perfectly aligned stripes */
          stripes_dst[i-skip].nstrides += stripes_src[i].nstrides;
          ++skip;
        } else
          stripes_dst[i-skip+1] = stripes_src[i];
      }
  }
  return num_stripes - (skip - 1);
}



/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
