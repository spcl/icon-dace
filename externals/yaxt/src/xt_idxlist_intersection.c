/**
 * @file xt_idxlist_intersection.c
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

#include <stdlib.h>

#include "instr.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt_idxlist_unpack.h"
#include "xt/xt_idxempty.h"
#include "xt_idxempty_internal.h"
#include "xt_idxlist_collection_internal.h"
#include "xt/xt_idxvec.h"
#include "xt_config_internal.h"
#include "xt_idxvec_internal.h"
#include "xt_idxsection_internal.h"
#include "xt_idxstripes_internal.h"

typedef Xt_idxlist (*intersection_get)(Xt_idxlist idxlist_src,
                                       Xt_idxlist idxlist_dst,
                                       Xt_config config);

#define empty_isect ((intersection_get)(void (*)(void))xt_idxempty_new)

static const intersection_get
intersection_get_matrix[num_idxlist_classes][num_idxlist_classes] = {
  { empty_isect, empty_isect, empty_isect, empty_isect, empty_isect },
  { empty_isect, xt_idxvec_get_intersection, xt_default_isect, xt_default_isect,
    xt_default_isect },
  { empty_isect, xt_default_isect, /* xt_idxlist_collection_get_intersection, */
    xt_default_isect, xt_default_isect, xt_default_isect },
  { empty_isect, xt_idxsection_get_intersection_with_other_idxlist,
    xt_idxsection_get_intersection_with_other_idxlist,
    xt_idxsection_get_intersection,
    xt_idxsection_get_intersection_with_other_idxlist },
  { empty_isect, xt_default_isect, xt_default_isect, xt_default_isect,
    xt_idxstripes_get_intersection },
};

void xt_idxlist_intersection_init(void)
{
}

Xt_idxlist
xt_idxlist_get_intersection(Xt_idxlist idxlist_src,
                            Xt_idxlist idxlist_dst)
{
  return xt_idxlist_get_intersection_custom(idxlist_src, idxlist_dst,
                                            &xt_default_config);
}

Xt_idxlist
xt_idxlist_get_intersection_custom(Xt_idxlist idxlist_src,
                                   Xt_idxlist idxlist_dst, Xt_config config)
{
  return intersection_get_matrix[idxlist_src->vtable->idxlist_pack_code]
    [idxlist_dst->vtable->idxlist_pack_code](idxlist_src, idxlist_dst, config);
}

Xt_idxlist
xt_default_isect(Xt_idxlist idxlist_src,
                 Xt_idxlist idxlist_dst,
                 Xt_config config)
{

   INSTR_DEF(instr_fallback,"xt_idxlist_get_intersection.fallback")

   // if the get_intersection routine was not able to compute the intersection
   INSTR_START(instr_fallback);

   int num_indices_src = xt_idxlist_get_num_indices(idxlist_src),
     num_indices_dst = xt_idxlist_get_num_indices(idxlist_dst);

   if (num_indices_src == 0 || num_indices_dst == 0)
     return xt_idxempty_new();

   Xt_idxlist intersection;
   if (num_indices_src < config->idxv_cnv_size
       && num_indices_dst < config->idxv_cnv_size) {
     Xt_int *indices_src
       = xmalloc(((size_t)num_indices_src + (size_t)num_indices_dst)
                 * sizeof (indices_src[0])),
       *indices_dst = indices_src + num_indices_src;

     xt_idxlist_get_indices(idxlist_src, indices_src);
     xt_idxlist_get_indices(idxlist_dst, indices_dst);

     Xt_idxlist idxvec_src = xt_idxvec_prealloc_new(indices_src, num_indices_src),
       idxvec_dst = xt_idxvec_prealloc_new(indices_dst, num_indices_dst);

     intersection
       = xt_idxvec_get_intersection(idxvec_src, idxvec_dst, config);

     xt_idxlist_delete(idxvec_src);
     xt_idxlist_delete(idxvec_dst);
     free(indices_src);
   } else {
     int num_stripes_src, num_stripes_dst;
     struct Xt_stripe *stripes_src, *stripes_dst;
     xt_idxlist_get_index_stripes(idxlist_src, &stripes_src, &num_stripes_src);
     xt_idxlist_get_index_stripes(idxlist_dst, &stripes_dst, &num_stripes_dst);
     Xt_idxlist idxstripes_src = xt_idxstripes_prealloc_new(stripes_src,
                                                         num_stripes_src),
       idxstripes_dst = xt_idxstripes_prealloc_new(stripes_dst,
                                                   num_stripes_dst);
     intersection
       = xt_idxstripes_get_intersection(idxstripes_src, idxstripes_dst, config);
     xt_idxlist_delete(idxstripes_dst);
     xt_idxlist_delete(idxstripes_src);
     free(stripes_dst);
     free(stripes_src);
   }

   INSTR_STOP(instr_fallback);
   return intersection;
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

