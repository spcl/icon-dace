/**
 * @file xt_idxstripes_internal.h
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
#ifndef XT_IDXSTRIPES_INTERNAL_H
#define XT_IDXSTRIPES_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>

#include "core/ppm_visibility.h"
#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"
#include "xt_arithmetic_util.h"

PPM_DSO_INTERNAL void
xt_idxstripes_initialize(void);

PPM_DSO_INTERNAL void
xt_idxstripes_finalize(void);

PPM_DSO_INTERNAL Xt_idxlist
xt_idxstripes_get_intersection(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst,
                               Xt_config config);

/**
 * Generates an index list that is built up of stripes of indices.
 * Does not copy the stripes.
 *
 * @param[in] stripes     array defining the stripes
 * @param[in] num_stripes number of stripes
 */
PPM_DSO_INTERNAL Xt_idxlist
xt_idxstripes_prealloc_new(struct Xt_stripe const * stripes, int num_stripes);

/**
 * can b be merged into a, i.e. do a and b form a continuous range?
 *
 * @param a position extent to merge into
 * @param b position extent to merge
 */
static inline bool
xt_can_merge_pos_ext(struct Xt_pos_ext a, struct Xt_pos_ext b)
{
  return (((b.start == a.start + a.size)
           | ((abs(a.size) == 1) & (b.start == a.start - a.size)))
          & ((abs(b.size) == 1) | (isign(a.size) == isign(b.size))));
}

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
