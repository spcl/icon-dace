/**
 * @file xt_idxmod.h
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

#ifndef XT_IDXMOD_H
#define XT_IDXMOD_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"

/** \example test_idxmod.c
 */
/** \example test_idxmod_f.f90
 */

/**
 * A modifier is a pair of compact idxlists that describes an index-to-index
 * stream.
 * Additionally, there is a mask value that can be used to track a sequence of
 * modifications (see \ref xt_idxmod_new).
 */
struct Xt_modifier {
  Xt_idxlist extract; ///<  idx values
  Xt_idxlist subst;   ///<  idx values, must describe at least as many indices
                      ///< as the extract idxlist
  int mask;           ///< modifier mask, can be used to identify modifier
                      ///< usage in the result
};

/**
 * \brief generates a new index list based on an index list and a sequence
 *        of modifiers
 *
 * In the new index list all occurences of indices listed in the
 * \ref Xt_modifier.extract "modifier extract list" are replaced by respective
 * indices of the modifier substitute list. Multiple modifiers are applied
 * consecutively to the input index list. If a modification state array is
 * given then it will be updated using the mask values in the used modifiers
 * (OR operation). An update only affects indices that are modified by the
 * respective modifier.
 *
 * @param[in]     patch_idxlist  input index list
 * @param[in]     modifier       array of modifiers
 * @param[in]     modifier_num   number of modifiers
 * @param[in,out] mstate         modification state
 * @return returns a \ref Xt_idxlist
 * @remark For performance considerations imagine patch_idxlist as a (small)
 *         task-local index list whose transform into a vector scales ideally
 *         and the modifier idxlist-components as compact idxlists of global
 *         extent which should never downgrade to an elementwise index
 *         description.
 * @remark mstate is not initilaized
 * @remark mstate == NULL is valid (disables the modification state handling)
 * @remark if mstate != NULL then it is required to be at least as big as the
 *         input index list (patch_idxlist)
 */

Xt_idxlist xt_idxmod_new(Xt_idxlist patch_idxlist,
                         struct Xt_modifier *modifier,
                         int modifier_num, int *mstate);

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
