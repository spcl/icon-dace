/**
 * @file xt_stripe_util.h
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
#ifndef XT_STRIPE_UTIL_H
#define XT_STRIPE_UTIL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>
#include <stdlib.h>

#include "xt/xt_stripe.h"

#include "core/ppm_visibility.h"
#include "xt_arithmetic_util.h"

struct Xt_stripe_minmax {
  Xt_int min, max;
};

static inline struct Xt_stripe_minmax
xt_stripe2minmax(struct Xt_stripe stripe)
{
  Xt_int min = stripe.start, max = min;
  Xt_int stride = stripe.stride;
  int nstrides = stripe.nstrides;
  Xt_int mask = Xt_isign_mask(stride);
  stride = (Xt_int)((Xt_int)(nstrides - 1) * stride);
  min = (Xt_int)(min + (Xt_int)(mask & stride));
  max = (Xt_int)(max + (Xt_int)((Xt_int)~mask & stride));
  return (struct Xt_stripe_minmax){ .min = min, .max = max};
}


static inline int
xt_stripes_overlap(struct Xt_stripe a, struct Xt_stripe b) {
  struct Xt_stripe_minmax mma = xt_stripe2minmax(a),
    mmb = xt_stripe2minmax(b);
  return (mma.max >= mmb.min) & (mma.min <= mmb.max);
}

PPM_DSO_INTERNAL void
xt_convert_indices_to_stripes_keep_buf(const Xt_int *restrict indices,
                                       int num_indices,
                                       struct Xt_stripe **stripes,
                                       int * num_stripes);

/**
 * copy stripes_src to stripes_dst, fusing trivially adjacent stripes
 * (i.e. having same stride and matching bounds)
 *
 * @param num_stripes number of stripes stored at stripes_src
 * @param stripes_dst target array able to hold at least num_stripes stripes
 * @param stripes_src source array containing \a num_stripes stripes to
 * be copied to \a stripes_dst
 * @param lookback if true, inspects also stripes_dst[-1] for possible
 * fusion with stripes_src[0]
 * @return number of stripes written to stripes_dst
 *
 */
PPM_DSO_INTERNAL size_t
xt_stripes_merge_copy(size_t num_stripes,
                      struct Xt_stripe *stripes_dst,
                      const struct Xt_stripe *stripes_src,
                      bool lookback);


#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
