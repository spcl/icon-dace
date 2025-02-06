/**
 * @file xt_idxlist.c
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
#include "config.h"
#endif

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include "xt/xt_core.h"
#include "xt/xt_stripe.h"
#include "xt_stripe_util.h"
#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_idxvec.h"
#include "xt_idxvec_internal.h"
#include "xt/xt_idxstripes.h"
#include "xt_idxstripes_internal.h"
#include "xt/xt_mpi.h"
#include "xt_idxlist_unpack.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "instr.h"

void xt_idxlist_delete(Xt_idxlist idxlist) {

   idxlist->vtable->delete(idxlist);
}

size_t xt_idxlist_get_pack_size(Xt_idxlist idxlist,
                                MPI_Comm comm) {

   return idxlist->vtable->get_pack_size(idxlist, comm);
}

void xt_idxlist_pack(Xt_idxlist idxlist, void *buffer,
                     int buffer_size, int *position,
                     MPI_Comm comm) {

   idxlist->vtable->pack(idxlist, buffer, buffer_size,
                         position, comm);
}

Xt_idxlist xt_idxlist_copy(Xt_idxlist idxlist) {

   return idxlist->vtable->copy(idxlist);
}

int xt_idxlist_get_num_indices(Xt_idxlist idxlist) {
   return idxlist->num_indices;
}

void xt_idxlist_get_indices(Xt_idxlist idxlist, Xt_int *indices) {

   idxlist->vtable->get_indices(idxlist, indices);
}


const Xt_int *xt_idxlist_get_indices_const(Xt_idxlist idxlist) {
  if (idxlist->vtable->get_indices != NULL)
    return idxlist->vtable->get_indices_const (idxlist);

  die("xt_idxlist_get_indices_const: fatal error: "
      "get_indices_const not implemented");
  return NULL;
}


void xt_idxlist_get_index_stripes(Xt_idxlist idxlist,
                                  struct Xt_stripe ** stripes,
                                  int * num_stripes) {

  *stripes = NULL;
  *num_stripes = 0;
  xt_idxlist_get_index_stripes_keep_buf(idxlist, stripes, num_stripes);
}

void xt_idxlist_get_index_stripes_keep_buf(Xt_idxlist idxlist,
                                           struct Xt_stripe ** stripes,
                                           int * num_stripes) {

  INSTR_DEF(instr_fallback,"xt_idxlist_get_index_stripes.fallback")

  if (idxlist->vtable->get_index_stripes != NULL) {
    idxlist->vtable->get_index_stripes(idxlist, stripes, num_stripes);

  } else { // fall-back solution

    INSTR_START(instr_fallback);

    int num_indices = xt_idxlist_get_num_indices(idxlist);

    Xt_int *indices = xmalloc((size_t)num_indices * sizeof (indices[0]));

    xt_idxlist_get_indices(idxlist, indices);

    xt_convert_indices_to_stripes_keep_buf(indices, num_indices,
                                           stripes, num_stripes);

    free(indices);

    INSTR_STOP(instr_fallback);
  }
}




int xt_idxlist_get_index_at_position(Xt_idxlist idxlist, int position,
                                     Xt_int *idx) {

  return idxlist->vtable->get_index_at_position(idxlist, position, idx);
}

int
xt_idxlist_get_indices_at_positions(Xt_idxlist idxlist, const int *positions,
                                    int num_pos, Xt_int *indices,
                                    Xt_int undef_idx) {

  INSTR_DEF(instr_fallback,"xt_idxlist_get_intersection.fallback")

  if ( idxlist->vtable->get_indices_at_positions ) {

    return idxlist->vtable->get_indices_at_positions(idxlist, positions,
                                                     num_pos, indices,
                                                     undef_idx);

  } else {

    INSTR_START(instr_fallback);

    // fallback solution using xt_idxlist_get_index_at_position:
    int undef_count = 0;
    for (int ip=0; ip<num_pos; ip++) {
      if (xt_idxlist_get_index_at_position(idxlist, positions[ip],
                                           &indices[ip]) != 0) {
        indices[ip] = undef_idx;
        undef_count++;
      }
    }

    INSTR_STOP(instr_fallback);
    return undef_count;
  }
}

int xt_idxlist_get_position_of_index(Xt_idxlist idxlist, Xt_int idx,
                                     int * position) {

  return idxlist->vtable->get_position_of_index(idxlist, idx, position);
}

int
xt_idxlist_get_positions_of_indices(Xt_idxlist idxlist, Xt_int const * indices,
                                    int num_indices, int * positions,
                                    int single_match_only) {

  INSTR_DEF(instr_fallback,"xt_idxlist_get_positions_of_indices.fallback")

  size_t num_failed_lookups;
  if (num_indices > 0) {
    if (idxlist->vtable->get_positions_of_indices != NULL)
      num_failed_lookups
        = idxlist->vtable->get_positions_of_indices(idxlist, indices,
                                                    (size_t)num_indices,
                                                    positions,
                                                    single_match_only);
    else {
      INSTR_START(instr_fallback);

      int num_tmp_indices
        = xt_idxlist_get_num_indices(idxlist);
      Xt_int *tmp_indices
        = xmalloc((size_t)num_tmp_indices * sizeof (tmp_indices[0]));
      xt_idxlist_get_indices(idxlist, tmp_indices);
      Xt_idxlist tmp_idxvec
        = xt_idxvec_prealloc_new(tmp_indices, num_tmp_indices);
      num_failed_lookups
        = tmp_idxvec->vtable->get_positions_of_indices(tmp_idxvec, indices,
                                                       (size_t)num_indices,
                                                       positions,
                                                       single_match_only);
      xt_idxlist_delete(tmp_idxvec);
      free(tmp_indices);

      INSTR_STOP(instr_fallback);
    }
  } else
    num_failed_lookups = 0;
  assert(num_failed_lookups <= INT_MAX);
  return (int)num_failed_lookups;
}

int
xt_idxlist_get_pos_exts_of_index_stripes(Xt_idxlist idxlist,
                                         int num_stripes,
                                         const struct Xt_stripe stripes[num_stripes],
                                         int *num_ext,
                                         struct Xt_pos_ext **pos_ext,
                                         int single_match_only)
{
  if (idxlist->vtable->get_pos_exts_of_index_stripes)
    return idxlist->vtable->get_pos_exts_of_index_stripes(
      idxlist,
      num_stripes, stripes, num_ext, pos_ext, single_match_only);
  else {
    int num_tmp_stripes;
    struct Xt_stripe *tmp_stripes;
    xt_idxlist_get_index_stripes(idxlist, &tmp_stripes, &num_tmp_stripes);
    Xt_idxlist idxlist_stripes
      = xt_idxstripes_prealloc_new(tmp_stripes, num_tmp_stripes);
    int retval =
      idxlist_stripes->vtable->get_pos_exts_of_index_stripes(
        idxlist_stripes,
        num_stripes, stripes, num_ext, pos_ext, single_match_only);
    xt_idxlist_delete(idxlist_stripes);
    free(tmp_stripes);
    return retval;
  }
}

int xt_idxlist_get_position_of_index_off(Xt_idxlist idxlist, Xt_int idx,
                                         int * position, int offset) {

  return idxlist->vtable->get_position_of_index_off(idxlist, idx,
                                                    position, offset);
}

int
xt_idxlist_get_positions_of_indices_off(Xt_idxlist idxlist,
                                        const Xt_int *indices,
                                        int num_indices, int * positions,
                                        int * offsets) {

  INSTR_DEF(instr_fallback,"xt_idxlist_get_intersection.fallback")

  if (idxlist->vtable->get_positions_of_indices_off != NULL)
    return idxlist->vtable->get_positions_of_indices_off(idxlist, indices,
                                                         num_indices,
                                                         positions, offsets);
  INSTR_START(instr_fallback);

  int ret_val = 0;
  for (int i = 0; i < num_indices; ++i) {

    int temp_position, offset;

    if (offsets == NULL)
      offset = 0;
    else
      offset = offsets[i];

    ret_val
      = idxlist->vtable->get_position_of_index_off(idxlist, indices[i],
                                                   &temp_position, offset);

    if (ret_val) break;
    positions[i] = temp_position;
  }

  INSTR_STOP(instr_fallback);
  return ret_val;
}

Xt_int xt_idxlist_get_min_index(Xt_idxlist idxlist) {

  return idxlist->vtable->get_min_index(idxlist);
}

Xt_int xt_idxlist_get_max_index(Xt_idxlist idxlist) {

  return idxlist->vtable->get_max_index(idxlist);
}

static void
get_position_in_ndim_space(Xt_int idx, unsigned ndim,
                           Xt_int global_stride[ndim],
                           Xt_int global_start_index, Xt_int position[ndim]) {

  idx = (Xt_int)(idx - global_start_index);

  for (size_t i = 0; i < ndim - 1; ++i) {
    position[i] = (Xt_int)(idx / global_stride[i]);
    idx = (Xt_int)(idx % global_stride[i]);
  }

   position[ndim - 1] = idx;
}

void xt_idxlist_get_bounding_box(Xt_idxlist idxlist, unsigned ndim,
                                 const Xt_int global_size[ndim],
                                 Xt_int global_start_index,
                                 struct Xt_bounds bounds[ndim]) {

   INSTR_DEF(instr_fallback,"xt_idxlist_get_bounding_box.fallback")

   int num_indices = xt_idxlist_get_num_indices(idxlist);

   if (num_indices == 0) {

      for (size_t i = 0; i < ndim; ++i) {
         bounds[i].start = 0;
         bounds[i].size = 0;
      }

      return;

   } else if (idxlist->vtable->get_bounding_box != NULL) {

      idxlist->vtable->get_bounding_box(idxlist, ndim, global_size,
                                        global_start_index, bounds);
      return;
   }

   INSTR_START(instr_fallback);

   Xt_int global_stride[ndim];

   global_stride[ndim - 1] = 1;

   for (size_t i = ndim - 2; i < ndim; --i)
     global_stride[i] = (Xt_int)(global_stride[i+1] * global_size[i+1]);

   Xt_int curr_index;
   Xt_int curr_position[ndim];

   xt_idxlist_get_index_at_position(idxlist, 0, &curr_index);
   get_position_in_ndim_space(curr_index, ndim, global_stride,
                              global_start_index, curr_position);

   for (size_t i = 0; i < ndim; ++i) {
      bounds[i].start = curr_position[i];
      bounds[i].size = 1;
   }

  for (int j = 1; j < num_indices; ++j) {

    xt_idxlist_get_index_at_position(idxlist, j, &curr_index);
    get_position_in_ndim_space(curr_index, ndim, global_stride,
                               global_start_index, curr_position);

    for (size_t i = 0; i < ndim; ++i) {

      if (curr_position[i] < bounds[i].start) {

        bounds[i].size = (Xt_int)(bounds[i].size + curr_position[i] - bounds[i].start);
        bounds[i].start = curr_position[i];

      } else if (curr_position[i] >= bounds[i].start + bounds[i].size) {

        bounds[i].size = (Xt_int)(curr_position[i] - bounds[i].start + 1);
      }
    }
  }

  INSTR_STOP(instr_fallback);
}

static Xt_uid nextId = UINT64_C(1);
#ifdef HAVE_PTHREAD
static pthread_mutex_t nextIdMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

Xt_uid
xt_idxlist_new_uid(void)
{
  Xt_uid thisId;
#ifdef HAVE_PTHREAD
  if (pthread_mutex_lock(&nextIdMutex))
    die("unexpected pthread locking error");
#endif
  thisId = nextId;
  if (!++nextId)
    die("unique ID counter overflow");
#ifdef HAVE_PTHREAD
  if (pthread_mutex_unlock(&nextIdMutex))
    die("unexpected pthread locking error");
#endif
  return thisId;
}

Xt_uid
xt_idxlist_get_uid(Xt_idxlist idxlist)
{
  assert(idxlist->uid);
  return idxlist->uid;
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
