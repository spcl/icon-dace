/**
 * @file test_idxlist_utils.c
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
#include <string.h>
#include <yaxt.h>
#include "core/ppm_xfuncs.h"
#include "tests.h"
#include "test_idxlist_utils.h"

static void
ref_get_positions_of_indices(Xt_int const * indices, int num_indices,
                             Xt_int const * selection, int selection_size,
                             int * positions, int repeat_first_match) {
  int used[num_indices > 0 ? num_indices : 1];
  int found;

  for (int i = 0; i < num_indices; i++)
    used[i] = 0;

  for (int i = 0; i < selection_size; i++) {

    found = 0;
    for (int j = 0; j < num_indices; j++) {
      if (selection[i] == indices[j] && (repeat_first_match || !used[j])) {
        positions[i] = j;
        used[j] = 1;
        found = 1;
        break;
      }
    }
    if (!found)
      PUT_ERR("ref_get_positions_of_indices: internal test error\n");
  }
}

static Xt_int get_min_index(int num_indices, const Xt_int *indices) {

  Xt_int index = num_indices > 0 ? indices[0] : 0;

  for (int i = 1; i < num_indices; ++i)
    if (index > indices[i])
      index = indices[i];

  return index;
}

static Xt_int get_max_index(int num_indices, const Xt_int *indices) {

  Xt_int index = num_indices > 0 ? indices[0] : 0;

  for (int i = 1; i < num_indices; ++i)
    if (index < indices[i])
      index = indices[i];

  return index;
}

void check_idxlist(Xt_idxlist idxlist, Xt_int const * ref_indices,
                   int ref_num_indices) {

  int num_indices = xt_idxlist_get_num_indices(idxlist);

  if (num_indices != ref_num_indices)
    PUT_ERR("wrong number of indices in idxlist\n");

  size_t ref_num_indices_ = (size_t)(ref_num_indices > 0 ? ref_num_indices : 1);
  Xt_int indices[ref_num_indices_];
  Xt_int index;

  xt_idxlist_get_indices(idxlist, indices);
  const Xt_int *indices_const = xt_idxlist_get_indices_const(idxlist);
  {
    bool mismatch = false;
    for(int i=0; i<ref_num_indices; i++)
      mismatch |= (ref_indices[i] != indices_const[i]);
    if (mismatch)
      PUT_ERR("check_idxlist: ref_indices[i] != indices_const[i]\n");
  }

  {
    int ref_positions[2][ref_num_indices_];
    {

      int position;

      ref_get_positions_of_indices(indices, num_indices, ref_indices,
                                   ref_num_indices, ref_positions[1], 1);

      bool mismatch = false;
      for (int i = 0; i < ref_num_indices; ++i)
        mismatch |= (indices[i] != ref_indices[i]);
      if (mismatch)
        PUT_ERR("index list contains wrong indices\n");

      for (int i = 0; i < ref_num_indices; ++i) {

        if (xt_idxlist_get_index_at_position(idxlist, (int)i, &index)
            || index != ref_indices[i])
          PUT_ERR("error in xt_idxlist_get_index_at_position\n");

        if (xt_idxlist_get_position_of_index(idxlist, ref_indices[i], &position)
            || position != ref_positions[1][i])
          PUT_ERR("error in xt_idxlist_get_position_of_index\n");
      }

      ref_get_positions_of_indices(indices, num_indices, ref_indices,
                                   ref_num_indices, ref_positions[0], 0);
      position = -1;
      for (int i = 0; i < ref_num_indices; ++i) {
        if (xt_idxlist_get_position_of_index_off(idxlist, ref_indices[i],
                                                 &position, position+1))
          PUT_ERR("error in xt_idxlist_get_position_of_index_off\n");

        if (position != ref_positions[0][i])
          PUT_ERR("get_positions_of_indices returned wrong position\n");
      }
    }

    {
      int positions[ref_num_indices_];

      if (xt_idxlist_get_positions_of_indices(idxlist, ref_indices,
                                              ref_num_indices, positions, 1))
        PUT_ERR("get_positions_of_indices failed\n");

      bool mismatch = false;
      for (int i = 0; i < ref_num_indices; ++i)
        mismatch |= (positions[i] != ref_positions[0][i]);
      if (mismatch)
        PUT_ERR("check_idxlist: get_positions_of_indices failed (case 1)\n");
    }

    {
      int positions[ref_num_indices_];

      if (xt_idxlist_get_positions_of_indices(idxlist, ref_indices,
                                              ref_num_indices, positions, 0))
        PUT_ERR("get_positions_of_indices failed\n");

      bool mismatch = false;
      for (int i = 0; i < ref_num_indices; ++i)
        mismatch |= (positions[i] != ref_positions[1][i]);
      if (mismatch)
        PUT_ERR("check_idxlist: get_positions_of_indices failed (case 2)\n");
    }

    {
      Xt_int reverse_ref_indices[ref_num_indices_];

      for (int i = 0; i < ref_num_indices; ++i)
        reverse_ref_indices[i] = ref_indices[ref_num_indices - i - 1];

      {
        int  positions[ref_num_indices_];

        ref_get_positions_of_indices(indices, num_indices, reverse_ref_indices,
                                     ref_num_indices, ref_positions[0], 0);

        if (xt_idxlist_get_positions_of_indices(idxlist, reverse_ref_indices,
                                                ref_num_indices, positions, 1))
          PUT_ERR("get_positions_of_indices failed\n");

        bool mismatch = false;
        for (int i = 0; i < ref_num_indices; ++i)
          mismatch |= (positions[i] != ref_positions[0][i]);
        if (mismatch)
          PUT_ERR("check_idxlist: get_positions_of_indices failed (case 3)\n");
      }

      {
        int positions[ref_num_indices_];

        if (xt_idxlist_get_positions_of_indices(idxlist, reverse_ref_indices,
                                                ref_num_indices, positions, 0))
          PUT_ERR("get_positions_of_indices failed\n");

        bool mismatch = false;
        for (int i = 0; i < ref_num_indices; ++i)
          mismatch |= (positions[i] != ref_positions[1][ref_num_indices-1-i]);
        if (mismatch)
          PUT_ERR("check_idxlist: get_positions_of_indices failed (case 4)\n");
      }
    }
  }

  if (ref_num_indices > 0) {

    if (get_min_index(ref_num_indices, ref_indices) !=
        xt_idxlist_get_min_index(idxlist))
      PUT_ERR("check_idxlist: get_min_index failed\n");

    if (get_max_index(ref_num_indices, ref_indices) !=
        xt_idxlist_get_max_index(idxlist))
      PUT_ERR("check_idxlist: get_max_index failed\n");
  }

  if (!xt_idxlist_get_index_at_position(idxlist, -1, &index))
    PUT_ERR("error in xt_idxlist_get_index_at_position\n");

  if (!xt_idxlist_get_index_at_position(idxlist, (int)ref_num_indices, &index))
    PUT_ERR("error in xt_idxlist_get_index_at_position\n");
}

void
check_stripes(struct Xt_stripe const * stripes, int num_stripes,
              struct Xt_stripe const * ref_stripes, int ref_num_stripes) {

  Xt_int i;

  if (num_stripes != ref_num_stripes)
    PUT_ERR("wrong number of stripes\n");

  bool mismatch = false;
  for (i = 0; i < ref_num_stripes; ++i)
    mismatch |= (stripes[i].start != ref_stripes[i].start ||
                 stripes[i].nstrides != ref_stripes[i].nstrides ||
                 stripes[i].stride != ref_stripes[i].stride);
  if (mismatch)
    PUT_ERR("error in stripe (start, nstrides and/or stride)\n");
}

void
check_offsets(size_t num_offsets, const int offsets_a[num_offsets],
              const int offsets_b[num_offsets])
{
  bool mismatch = false;
  for(size_t i=0; i<num_offsets; i++)
    mismatch |= (offsets_a[i] != offsets_b[i]);
  if (mismatch)
    PUT_ERR("unexpected results for offsets or index positions\n");
}

Xt_idxlist
idxlist_pack_unpack_copy(Xt_idxlist idxlist)
{
  size_t buffer_size
    = xt_idxlist_get_pack_size(idxlist, MPI_COMM_WORLD);

  // allocate send buffers
  void *send_buffer = xmalloc(buffer_size),
    *recv_buffer = xmalloc(buffer_size);

  // pack the index list
  int position = 0;

  assert(buffer_size <= INT_MAX);
  xt_idxlist_pack(idxlist, send_buffer, (int)buffer_size, &position,
                  MPI_COMM_WORLD);
  assert((size_t)position <= buffer_size);

  // send the buffer
  int rank;
  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  MPI_Request request;
  xt_mpi_call(MPI_Isend(send_buffer, (int)buffer_size, MPI_PACKED, rank, 0,
                        MPI_COMM_WORLD, &request), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Request_free(&request), MPI_COMM_WORLD);
  // receive the buffer
  xt_mpi_call(MPI_Recv(recv_buffer, (int)buffer_size, MPI_PACKED, rank, 0,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE), MPI_COMM_WORLD);
  free(send_buffer);
  // unpack the buffer
  position = 0;
  Xt_idxlist idxlist_copy
    = xt_idxlist_unpack(recv_buffer, (int)buffer_size, &position,
                        MPI_COMM_WORLD);
  assert((size_t)position <= buffer_size);
  free(recv_buffer);
  return idxlist_copy;
}


static int
xt_int_compare(const void *a, const void *b)
{
  return (*(const Xt_int *)a > *(const Xt_int *)b)
    - (*(const Xt_int *)a < *(const Xt_int *)b);
}

void
check_idxlist_copy(Xt_idxlist idxlist, Xt_idxlist idxlist_copy,
                   size_t num_ref_indices, const Xt_int *ref_indices,
                   size_t num_ref_stripes, const struct Xt_stripe *ref_stripes)
{
  // check received collection list
  check_idxlist(idxlist_copy, ref_indices, (int)num_ref_indices);

  // compute intersection between the two index lists
  Xt_idxlist intersection
    = xt_idxlist_get_intersection(idxlist, idxlist_copy);

  // check intersection data
  Xt_int *sorted_ref_indices = xmalloc(num_ref_indices
                                       * sizeof (ref_indices[0]));
  memcpy(sorted_ref_indices, ref_indices,
         num_ref_indices * sizeof (ref_indices[0]));
  qsort(sorted_ref_indices, num_ref_indices, sizeof (ref_indices[0]),
        xt_int_compare);
  check_idxlist(intersection, sorted_ref_indices, (int)num_ref_indices);
  free(sorted_ref_indices);
  xt_idxlist_delete(intersection);

  // check the conversion to stripes
  struct Xt_stripe * stripes;
  int num_stripes;
  xt_idxlist_get_index_stripes(idxlist, &stripes, &num_stripes);

  check_stripes(stripes, num_stripes, ref_stripes, (int)num_ref_stripes);
  // clean up
  free(stripes);
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
