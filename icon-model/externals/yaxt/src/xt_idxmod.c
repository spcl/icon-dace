/**
 * @file xt_idxmod.c
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

#include <assert.h>

#include "xt/xt_core.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"

#include "xt/xt_idxmod.h"


Xt_idxlist xt_idxmod_new(Xt_idxlist patch_idxlist,
                         struct Xt_modifier *modifier,
                         int modifier_num, int *mstate) {
  // ensure that yaxt is initialized
  assert(xt_initialized());

  int patch_size = xt_idxlist_get_num_indices(patch_idxlist);

  // if there is no modifier then we just give back a copy of the original:
  if (modifier_num<1) return xt_idxlist_copy(patch_idxlist);

  Xt_int *workpatch_idx = xmalloc((size_t)patch_size * sizeof(Xt_int));
  Xt_int const* inter_idx;
  int *workpatch_pos = xmalloc((size_t)patch_size * sizeof(int));
  int *extract_pos = xmalloc((size_t)patch_size * sizeof(int));
  Xt_int *subst_idx = xmalloc((size_t)patch_size * sizeof(Xt_int));

  Xt_idxlist workpatch_idxlist = patch_idxlist;

  xt_idxlist_get_indices(workpatch_idxlist, workpatch_idx);

  for (int im = 0; im < modifier_num; im++) {
    struct Xt_modifier *m = &modifier[im];

    // intersection between extract values and workpatch:
    // any multiplicity of workpatch must be repeated in the intersection therefore
    // workpatch must have the target role
    Xt_idxlist intersection_idxlist
      = xt_idxlist_get_intersection(m->extract,workpatch_idxlist);

    // get intersection index array => inter_idx:
    int intersection_size = xt_idxlist_get_num_indices(intersection_idxlist);
    if (intersection_size > patch_size) die("xt_idxmod_new: internal error: (intersection_size > patch_size)");
    inter_idx = xt_idxlist_get_indices_const(intersection_idxlist);

    // get the intersection positions within the extract list
    // m->extract has source role, therefore single_match_only = 0
    // => extract_pos
    int missing = xt_idxlist_get_positions_of_indices(m->extract, inter_idx,
                                                      intersection_size,
                                                      extract_pos, 0);
    if (missing) die("xt_idxmod_new: internal error: cannot locate all intersection positions (1)");

    // get the intersection positions within workpatch
    // we must find each fitting index, so single_match_only = 1
    // => workpatch_pos
    missing = xt_idxlist_get_positions_of_indices(workpatch_idxlist, inter_idx,
                                                  intersection_size,
                                                  workpatch_pos, 1);

    if (missing) die("xt_idxmod_new: internal error: cannot locate all intersection positions (2)");

    // using the positions above, select indices within m->subst:
    // it is an error if we cannot access all positions, so the value of undef_idx does not matter (set to 0)
    int undef_num = xt_idxlist_get_indices_at_positions(m->subst, extract_pos,
                                                        intersection_size,
                                                        subst_idx, 0);
    if (undef_num) die("xt_idxmod_new: internal error: failed access: m->subst is too small");

    // delete workpatch_idxlist
    if (im > 0) {
      xt_idxlist_delete(workpatch_idxlist);
    }

    // substitude indices within workpatch_idx
    int p;
    int mask = m->mask;
    if ( mstate != NULL && mask != 0) {
      // we also update the modification state
      for (int i=0; i<intersection_size; i++) {
        p = workpatch_pos[i];
        workpatch_idx[p] = subst_idx[i];
        mstate[p] |= mask;
      }
    } else {
      for (int i=0; i<intersection_size; i++) {
        p = workpatch_pos[i];
        workpatch_idx[p] = subst_idx[i];
      }
    }
    workpatch_idxlist = xt_idxvec_new(workpatch_idx, patch_size);

    xt_idxlist_delete(intersection_idxlist);
  }

  free(subst_idx);
  free(extract_pos);
  free(workpatch_pos);
  free(workpatch_idx);

  return workpatch_idxlist;
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
