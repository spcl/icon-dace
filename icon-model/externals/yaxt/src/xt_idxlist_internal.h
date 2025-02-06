/**
 * @brief Provide non-public declarations common to all index lists.
 * @file xt_idxlist_internal.h
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
#ifndef XT_IDXLIST_INTERNAL_H
#define XT_IDXLIST_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "core/ppm_visibility.h"
#include "xt/xt_config.h"
#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"

struct xt_idxlist_vtable {

   void  (*delete)(Xt_idxlist);

   size_t (*get_pack_size)(Xt_idxlist, MPI_Comm);
   void (*pack)(Xt_idxlist, void*, int, int*, MPI_Comm);
   Xt_idxlist (*copy)(Xt_idxlist);
   void (*get_indices)(Xt_idxlist, Xt_int *indices);
   Xt_int const * (*get_indices_const)(Xt_idxlist idxlist);
   void (*get_index_stripes)(Xt_idxlist, struct Xt_stripe **, int *);
   int (*get_index_at_position)(Xt_idxlist, int, Xt_int *);
   int (*get_indices_at_positions)(Xt_idxlist idxlist, const int *positions,
                                   int num, Xt_int *index, Xt_int undef_idx);
   int (*get_position_of_index)(Xt_idxlist, Xt_int, int *);
   size_t (*get_positions_of_indices)(Xt_idxlist, Xt_int const *,
                                      size_t, int *, int);
   int (*get_pos_exts_of_index_stripes)(Xt_idxlist, int,
                                        const struct Xt_stripe *,
                                        int *, struct Xt_pos_ext **, int);
   int (*get_position_of_index_off)(Xt_idxlist, Xt_int, int *, int);
   int (*get_positions_of_indices_off)(Xt_idxlist, Xt_int const *, int,
                                      int *, int *);
   Xt_int (*get_min_index)(Xt_idxlist);
   Xt_int (*get_max_index)(Xt_idxlist);
   void (*get_bounding_box)(Xt_idxlist idxlist, unsigned ndim,
                            const Xt_int global_size[ndim],
                            Xt_int global_start_index,
                            struct Xt_bounds bounds[ndim]);
  int idxlist_pack_code;
};

PPM_DSO_INTERNAL void
xt_idxlist_get_index_stripes_keep_buf(Xt_idxlist idxlist,
                                      struct Xt_stripe ** stripes,
                                      int * num_stripes);

struct Xt_idxlist_ {
  const struct xt_idxlist_vtable *vtable;
  int num_indices;
  Xt_uid uid;
};

PPM_DSO_INTERNAL void xt_idxlist_intersection_init(void);

PPM_DSO_INTERNAL Xt_uid xt_idxlist_new_uid(void);

static inline void Xt_idxlist_init(Xt_idxlist idxlist,
                                   const struct xt_idxlist_vtable *vtable,
                                   int num_indices)
{
  idxlist->vtable = vtable;
  idxlist->uid = xt_idxlist_new_uid();
  idxlist->num_indices = num_indices;
}

PPM_DSO_INTERNAL Xt_idxlist
xt_default_isect(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst,
                 Xt_config config);

/*
 * for index lists this large, stripe representation generally gives
 * better performance for internal computations
 */
enum {
  CHEAP_VECTOR_SIZE = 128,
};

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
