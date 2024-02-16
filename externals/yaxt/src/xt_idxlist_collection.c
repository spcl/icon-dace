/**
 * @file xt_idxlist_collection.c
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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "mpi.h"

#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_idxlist_collection.h"
#include "xt_idxlist_collection_internal.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_mpi.h"
#include "xt_idxlist_unpack.h"
#include "xt_stripe_util.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "ensure_array_size.h"

static void
idxlist_collection_delete(Xt_idxlist data);

static size_t
idxlist_collection_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxlist_collection_pack(Xt_idxlist data, void *buffer, int buffer_size,
                        int *position, MPI_Comm comm);

static Xt_idxlist
idxlist_collection_copy(Xt_idxlist idxlist);

static void
idxlist_collection_get_indices(Xt_idxlist idxlist, Xt_int *indices);

static const Xt_int *
idxlist_collection_get_indices_const(Xt_idxlist idxlist);

static void
idxlist_collection_get_index_stripes(Xt_idxlist idxlist,
                                     struct Xt_stripe ** stripes,
                                     int * num_stripes);

static int
idxlist_collection_get_index_at_position(Xt_idxlist idxlist, int position,
                                         Xt_int * index);

static int
idxlist_collection_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                         int * position);

static int
idxlist_collection_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                             int * position, int offset);

static Xt_int
idxlist_collection_get_min_index(Xt_idxlist idxlist);

static Xt_int
idxlist_collection_get_max_index(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxlist_collection_vtable = {
  .delete                      = idxlist_collection_delete,
  .get_pack_size               = idxlist_collection_get_pack_size,
  .pack                        = idxlist_collection_pack,
  .copy                        = idxlist_collection_copy,
  .get_indices                 = idxlist_collection_get_indices,
  .get_indices_const           = idxlist_collection_get_indices_const,
  .get_index_stripes           = idxlist_collection_get_index_stripes,
  .get_index_at_position       = idxlist_collection_get_index_at_position,
  .get_indices_at_positions    = NULL,
  .get_position_of_index       = idxlist_collection_get_position_of_index,
  .get_positions_of_indices    = NULL,
  .get_position_of_index_off   = idxlist_collection_get_position_of_index_off,
  .get_positions_of_indices_off = NULL,
  .get_min_index               = idxlist_collection_get_min_index,
  .get_max_index               = idxlist_collection_get_max_index,
  .get_bounding_box            = NULL,
  .idxlist_pack_code           = COLLECTION,
};

typedef struct Xt_idxlist_collection_ *Xt_idxlist_collection;

struct Xt_idxlist_collection_ {

  struct Xt_idxlist_ parent;

  int num_idxlists;

  Xt_int *index_array_cache;
  Xt_idxlist idxlists[];
};

Xt_idxlist
xt_idxlist_collection_new(Xt_idxlist *idxlists, int num_idxlists) {
  // ensure that yaxt is initialized
  assert(xt_initialized());

  Xt_idxlist result;
  if (num_idxlists > 1)
  {
    Xt_idxlist_collection collectionlist
      = xmalloc(sizeof (*collectionlist)
                + (size_t)num_idxlists * sizeof (collectionlist->idxlists[0]));

    collectionlist->num_idxlists = num_idxlists;
    collectionlist->index_array_cache = NULL;

    long long num_indices = 0;
    for (int i = 0; i < num_idxlists; ++i)
    {
      collectionlist->idxlists[i] = xt_idxlist_copy(idxlists[i]);
      num_indices += idxlists[i]->num_indices;
    }

    assert(num_indices <= INT_MAX);
    Xt_idxlist_init(&collectionlist->parent, &idxlist_collection_vtable,
                    (int)num_indices);

    result = (Xt_idxlist)collectionlist;
  }
  else if (num_idxlists == 1)
    result = xt_idxlist_copy(idxlists[0]);
  else /* num_idxlists == 0 */
    result = xt_idxempty_new();
  return result;
}

static void
idxlist_collection_delete(Xt_idxlist data) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)data;

   int num_lists = collectionlist->num_idxlists;
   for (int i = 0; i < num_lists; ++i)
      xt_idxlist_delete(collectionlist->idxlists[i]);

   free(collectionlist->index_array_cache);
   free(collectionlist);
}

static size_t
idxlist_collection_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)data;

   int size_header, num_lists = collectionlist->num_idxlists;
   size_t size_idxlists = 0;

   xt_mpi_call(MPI_Pack_size(2, MPI_INT, comm, &size_header), comm);

   for (int i = 0; i < num_lists; ++i)
      size_idxlists
        += xt_idxlist_get_pack_size(collectionlist->idxlists[i], comm);

   return (size_t)size_header + size_idxlists;
}

static void
idxlist_collection_pack(Xt_idxlist data, void *buffer, int buffer_size,
                        int *position, MPI_Comm comm) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)data;
   int num_lists = collectionlist->num_idxlists;
   int header[2] = { COLLECTION, num_lists };

   xt_mpi_call(MPI_Pack(header, 2, MPI_INT, buffer,
                        buffer_size, position, comm), comm);

   for (int i = 0; i < num_lists; ++i)
      xt_idxlist_pack(collectionlist->idxlists[i], buffer, buffer_size,
                      position, comm);
}

Xt_idxlist
xt_idxlist_collection_unpack(void *buffer, int buffer_size, int *position,
                             MPI_Comm comm) {

  int num_lists;
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         &num_lists, 1, MPI_INT, comm), comm);

  Xt_idxlist_collection collectionlist
    = xmalloc(sizeof (*collectionlist)
              + (size_t)num_lists * sizeof (collectionlist->idxlists[0]));

  collectionlist->index_array_cache = NULL;
  collectionlist->num_idxlists = num_lists;

  long long num_indices = 0;
  for (int i = 0; i < num_lists; ++i) {
    collectionlist->idxlists[i] = xt_idxlist_unpack(buffer, buffer_size,
                                                    position, comm);
    num_indices += collectionlist->idxlists[i]->num_indices;
  }

  assert(num_indices <= INT_MAX);
  Xt_idxlist_init(&collectionlist->parent, &idxlist_collection_vtable,
                  (int)num_indices);
  return (Xt_idxlist)collectionlist;
}

Xt_idxlist
xt_idxlist_collection_get_intersection(Xt_idxlist XT_UNUSED(idxlist_src),
                                       Xt_idxlist XT_UNUSED(idxlist_dst),
                                       Xt_config XT_UNUSED(config)) {

  return NULL;
}

static Xt_idxlist
idxlist_collection_copy(Xt_idxlist idxlist) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;

   return xt_idxlist_collection_new(collectionlist->idxlists,
                                    collectionlist->num_idxlists);
}

static void
idxlist_collection_get_indices(Xt_idxlist idxlist, Xt_int *indices) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;
   /// \todo use memcpy with index_array_cache if available
   int offlist = 0, num_lists = collectionlist->num_idxlists;

   for (int i = 0; i < num_lists; ++i) {

      xt_idxlist_get_indices(collectionlist->idxlists[i], indices+offlist);
      offlist += xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
   }
}

static const Xt_int *
idxlist_collection_get_indices_const(Xt_idxlist idxlist) {

  Xt_idxlist_collection collection = (Xt_idxlist_collection)idxlist;

  if (collection->index_array_cache) return collection->index_array_cache;

  unsigned num_indices = (unsigned)idxlist->num_indices;

  Xt_int *tmp_index_array
    = xmalloc(num_indices * sizeof (collection->index_array_cache[0]));

  idxlist_collection_get_indices(idxlist, tmp_index_array);

  collection->index_array_cache = tmp_index_array;

  return collection->index_array_cache;
}


static void
idxlist_collection_get_index_stripes(Xt_idxlist idxlist,
                                     struct Xt_stripe ** stripes,
                                     int * num_stripes) {

  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;

  struct Xt_stripe * temp_stripes = NULL;
  size_t temp_stripes_array_size = 0;
  int num_temp_stripes = 0;

  xt_idxlist_get_index_stripes(collectionlist->idxlists[0], &temp_stripes,
                               &num_temp_stripes);
  temp_stripes_array_size = (size_t)num_temp_stripes;

  struct Xt_stripe * curr_stripes = NULL;
  int curr_num_stripes = 0, num_lists = collectionlist->num_idxlists;
  for (int i = 1; i < num_lists; ++i) {


    xt_idxlist_get_index_stripes_keep_buf(collectionlist->idxlists[i],
                                          &curr_stripes,
                                          &curr_num_stripes);

    ENSURE_ARRAY_SIZE(temp_stripes, temp_stripes_array_size,
                      num_temp_stripes + curr_num_stripes);

    curr_num_stripes
      = (int)xt_stripes_merge_copy((size_t)curr_num_stripes,
                                   temp_stripes + num_temp_stripes,
                                   curr_stripes,
                                   num_temp_stripes > 0);

    num_temp_stripes += curr_num_stripes;
  }

  free(curr_stripes);

  *stripes = xrealloc(temp_stripes,
                      (size_t)num_temp_stripes * sizeof(*temp_stripes));
  *num_stripes = num_temp_stripes;
}

static int
idxlist_collection_get_index_at_position(Xt_idxlist idxlist, int position,
                                         Xt_int * index) {

  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;
  int num_lists = collectionlist->num_idxlists;

  for (int i = 0; i < num_lists; ++i) {
    int n = xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
    if (position >= n)
      position -= n;
    else {
      return xt_idxlist_get_index_at_position(collectionlist->idxlists[i],
                                              position, index);
    }
  }
  return 1;

}

static int
idxlist_collection_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                             int * position, int offset) {

  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;

  int curr_num_indices = 0;

  int idxlist_offsets = 0;

  assert(offset >= 0);

  int i = 0, num_lists = collectionlist->num_idxlists;

  do {
    idxlist_offsets += curr_num_indices;
    curr_num_indices = xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
  } while (idxlist_offsets + curr_num_indices <= offset && ++i < num_lists);

  offset -= idxlist_offsets;

  for (;i < num_lists; ++i)
    if (!xt_idxlist_get_position_of_index_off(collectionlist->idxlists[i],
                                              index, position, offset)) {
      *position += idxlist_offsets;
      return 0;
    } else {
      idxlist_offsets
        += xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
      offset = 0;
    }

  return 1;
}

static int
idxlist_collection_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                         int * position) {

  return idxlist_collection_get_position_of_index_off(idxlist, index,
                                                      position, 0);
}

static Xt_int
idxlist_collection_get_min_index(Xt_idxlist idxlist) {

  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;

  int num_lists = collectionlist->num_idxlists;
  assert(num_lists > 0);

  Xt_int tmp_min, min = xt_idxlist_get_min_index(collectionlist->idxlists[0]);

  for (int i = 1; i < num_lists; ++i)
    if ((tmp_min = xt_idxlist_get_min_index(collectionlist->idxlists[i])) < min)
      min = tmp_min;

  return min;
}

static Xt_int
idxlist_collection_get_max_index(Xt_idxlist idxlist) {

  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;

  int num_lists = collectionlist->num_idxlists;
  assert(num_lists > 0);

  Xt_int tmp_max, max = xt_idxlist_get_max_index(collectionlist->idxlists[0]);

  for (int i = 1; i < num_lists; ++i)
    if ((tmp_max = xt_idxlist_get_max_index(collectionlist->idxlists[i])) > max)
      max = tmp_max;

  return max;
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
