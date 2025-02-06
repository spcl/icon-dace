/**
 * @file xt_idxempty.c
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
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_idxempty.h"
#include "xt_idxempty_internal.h"
#include "xt_idxlist_unpack.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_mpi.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"

static void
idxempty_delete(Xt_idxlist data);

static size_t
idxempty_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxempty_pack(Xt_idxlist data, void *buffer, int buffer_size,
              int *position, MPI_Comm comm);

static Xt_idxlist
idxempty_copy(Xt_idxlist idxlist);

static void
idxempty_get_indices(Xt_idxlist idxlist, Xt_int *indices);

static Xt_int const*
idxempty_get_indices_const(Xt_idxlist idxlist);

static void
idxempty_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                           int * num_stripes);

static int
idxempty_get_index_at_position(Xt_idxlist idxlist, int position, Xt_int * index);

static int
idxempty_get_indices_at_positions(Xt_idxlist idxlist, const int *positions,
                                  int num, Xt_int *index,
                                  Xt_int undef_idx);

static int
idxempty_get_position_of_index(Xt_idxlist idxlist, Xt_int index, int * position);

static int
idxempty_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                   int * position, int offset);

static size_t
idxempty_get_positions_of_indices(Xt_idxlist idxlist, Xt_int const * indices,
                                  size_t num_indices, int *positions,
                                  int single_match_only);

static int
idxempty_get_pos_exts_of_index_stripes(Xt_idxlist idxlist,
                                       int num_stripes,
                                       const struct Xt_stripe *stripes,
                                       int *num_ext,
                                       struct Xt_pos_ext **pos_ext,
                                       int single_match_only);

static int
idxempty_get_positions_of_indices_off(Xt_idxlist idxlist, Xt_int const * indices,
                                      int num_indices, int * positions,
                                      int * offsets);

static Xt_int
idxempty_get_min_index(Xt_idxlist idxlist);

static Xt_int
idxempty_get_max_index(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxempty_vtable = {
  .delete                       = idxempty_delete,
  .get_pack_size                = idxempty_get_pack_size,
  .pack                         = idxempty_pack,
  .copy                         = idxempty_copy,
  .get_indices                  = idxempty_get_indices,
  .get_indices_const            = idxempty_get_indices_const,
  .get_index_stripes            = idxempty_get_index_stripes,
  .get_index_at_position        = idxempty_get_index_at_position,
  .get_indices_at_positions     = idxempty_get_indices_at_positions,
  .get_position_of_index        = idxempty_get_position_of_index,
  .get_positions_of_indices     = idxempty_get_positions_of_indices,
  .get_pos_exts_of_index_stripes = idxempty_get_pos_exts_of_index_stripes,
  .get_position_of_index_off    = idxempty_get_position_of_index_off,
  .get_positions_of_indices_off = idxempty_get_positions_of_indices_off,
  .get_min_index                = idxempty_get_min_index,
  .get_max_index                = idxempty_get_max_index,
  .get_bounding_box             = NULL,
  .idxlist_pack_code            = EMPTY,
};

// index vector data structure
static struct xt_idxempty {
  struct Xt_idxlist_ parent;
} idxempty;

void
xt_idxempty_init(void)
{
  Xt_idxlist_init(&idxempty.parent, &idxempty_vtable, 0);
}

void
xt_idxempty_finalize(void)
{
}

Xt_idxlist xt_idxempty_new(void) {
  // ensure that yaxt is initialized
  assert(xt_initialized());
  return (void*)&idxempty;
}

static void idxempty_delete(Xt_idxlist XT_UNUSED(data)) {
}

static size_t idxempty_get_pack_size(Xt_idxlist XT_UNUSED(data), MPI_Comm comm)
{
  int size_int_type;

  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);

  return (size_t)size_int_type;
}

void idxempty_pack(Xt_idxlist data, void *buffer, int buffer_size,
                 int *position, MPI_Comm comm) {

  (void)data;
  assert(data);
  static const int type = EMPTY;

  xt_mpi_call(MPI_Pack(CAST_MPI_SEND_BUF(&type), 1, MPI_INT, buffer,
                       buffer_size, position, comm), comm);
}

Xt_idxlist xt_idxempty_unpack(void *XT_UNUSED(buffer),
                              int XT_UNUSED(buffer_size),
                              int *XT_UNUSED(position),
                              MPI_Comm XT_UNUSED(comm)) {

   return xt_idxempty_new();
}


static Xt_idxlist
idxempty_copy(Xt_idxlist XT_UNUSED(idxlist)) {

   return xt_idxempty_new();
}

static void
idxempty_get_indices(Xt_idxlist XT_UNUSED(idxlist), Xt_int *XT_UNUSED(indices)) {
  /* do nothing */
}

static Xt_int const*
idxempty_get_indices_const(Xt_idxlist XT_UNUSED(idxlist)) {

  return NULL;
}


static void
idxempty_get_index_stripes(Xt_idxlist XT_UNUSED(idxlist),
                           struct Xt_stripe **stripes,
                           int * num_stripes) {
  (void)stripes;
  *num_stripes = 0;
}

static int
idxempty_get_index_at_position(Xt_idxlist XT_UNUSED(idxlist),
                               int XT_UNUSED(position),
                               Xt_int *XT_UNUSED(index)) {

  return 1;
}

static int
idxempty_get_indices_at_positions(Xt_idxlist XT_UNUSED(idxlist),
                                  const int *XT_UNUSED(positions),
                                  int num_pos,
                                  Xt_int *index, Xt_int undef_idx) {

  assert(num_pos >= 0);

  for (int i = 0; i < num_pos; ++i)
    index[i] = undef_idx;

  return num_pos;
}

static int
idxempty_get_position_of_index_off(Xt_idxlist XT_UNUSED(idxlist),
                                   Xt_int XT_UNUSED(index),
                                   int *XT_UNUSED(position),
                                   int XT_UNUSED(offset)) {
  return 1;
}

static int
idxempty_get_position_of_index(Xt_idxlist XT_UNUSED(idxlist),
                               Xt_int XT_UNUSED(index),
                               int *XT_UNUSED(position)) {

  return 1;
}

size_t idxempty_get_positions_of_indices(Xt_idxlist XT_UNUSED(body_idxlist),
                                         const Xt_int *XT_UNUSED(selection_idx),
                                         size_t num_selection,
                                         int *XT_UNUSED(positions),
                                         int XT_UNUSED(single_match_only)) {

  return num_selection;
}

static int
idxempty_get_pos_exts_of_index_stripes(Xt_idxlist XT_UNUSED(idxlist),
                                       int num_stripes,
                                       const struct Xt_stripe *stripes,
                                       int *num_ext,
                                       struct Xt_pos_ext **pos_ext,
                                       int XT_UNUSED(single_match_only))
{
  *num_ext = 0;
  *pos_ext = NULL;
  unsigned num_idx = 0;
  if (num_stripes > 0)
    for (size_t i = 0; i < (size_t)num_stripes; ++i)
      num_idx += (unsigned)stripes[i].nstrides;
  return (int)num_idx;
}

static int
idxempty_get_positions_of_indices_off(Xt_idxlist XT_UNUSED(idxlist),
                                      const Xt_int *XT_UNUSED(indices),
                                      int XT_UNUSED(num_indices),
                                      int *XT_UNUSED(positions),
                                      int *XT_UNUSED(offsets)) {

  return 1;
}

static Xt_int
idxempty_get_min_index(Xt_idxlist XT_UNUSED(idxlist)) {

  die("idxempty_get_min_index: empty index list");
#ifndef __clang__
  return -1;
#endif
}

static Xt_int
idxempty_get_max_index(Xt_idxlist XT_UNUSED(idxlist)) {

  die("idxempty_get_max_index: empty index list");
#ifndef __clang__
  return -1;
#endif
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
