/**
 * @file test_xmap_common_parallel.c
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

#include <stdbool.h>
#include <stdlib.h>

#include <mpi.h>
#include <yaxt.h>

#include "tests.h"
#include "ctest_common.h"
#include "test_xmap_common.h"
#include "core/ppm_xfuncs.h"

static void
test_xmap_allgather_analog(xmap_constructor xmap_new,
                           size_t num_indices_per_rank, MPI_Comm comm);
static void
test_pair(xmap_constructor xmap_new, MPI_Comm comm);

static void
test_maxpos(xmap_constructor xmap_new, MPI_Comm comm, int indices_per_rank);

int
xt_xmap_parallel_test_main(xmap_constructor xmap_new)
{
  test_init_mpi(NULL, NULL, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  MPI_Comm comm = MPI_COMM_WORLD;
  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  test_xmap_allgather_analog(xmap_new, 1, comm);
  /* repeat test for large index list that will cause stripifying */
  test_xmap_allgather_analog(xmap_new, 1024, comm);

  if (comm_size > 2) // skip test if there are not enough processes
    test_ring_1d(xmap_new, comm);

  if (comm_size == 2)
    test_pair(xmap_new, comm);

  if (comm_size > 1)
    test_ping_pong(xmap_new, comm, 0, comm_size - 1);

  // test maxpos implementation for xt_xmap_intersection
  test_maxpos(xmap_new, comm, 5);
  // test maxpos implementation for xt_xmap_intersection_ext
  test_maxpos(xmap_new, comm, 501);
  xt_finalize();
  MPI_Finalize();
  return TEST_EXIT_CODE;
}

void
check_xmap_allgather_analog_xmap(Xt_xmap xmap, MPI_Comm comm)
{
  int comm_rank, comm_size, is_inter;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
  int (*get_comm_size)(MPI_Comm, int *)
    = is_inter ? MPI_Comm_remote_size : MPI_Comm_size;
  xt_mpi_call(get_comm_size(comm, &comm_size), comm);

  if (xt_xmap_get_num_destinations(xmap) != comm_size)
    PUT_ERR("error in xt_xmap_get_num_destinations\n");

  if (xt_xmap_get_num_sources(xmap) != comm_size)
    PUT_ERR("error in xt_xmap_get_num_sources\n");

  int *ranks = xmalloc(sizeof (*ranks) * (size_t)comm_size);

  xt_xmap_get_destination_ranks(xmap, ranks);
  bool mismatch = false;
  for (int i = 0; i < comm_size; ++i)
    mismatch |= (ranks[i] != i);
  if (mismatch)
    PUT_ERR("error in xt_xmap_get_destination_ranks\n");

  xt_xmap_get_source_ranks(xmap, ranks);
  mismatch = false;
  for (int i = 0; i < comm_size; ++i)
    mismatch |= (ranks[i] != i);
  if (mismatch)
    PUT_ERR("error in xt_xmap_get_source_ranks\n");
  free(ranks);
}


static void
test_xmap_allgather_analog(xmap_constructor xmap_new,
                           size_t num_indices_per_rank,
                           MPI_Comm comm)
{
  int comm_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);


  // test in which every process requests data from all processes
  // source index list
  Xt_int *src_index_list = xmalloc(num_indices_per_rank * sizeof (Xt_int));
  for (size_t i = 0; i < num_indices_per_rank; ++i)
    src_index_list[i]
      = (Xt_int)((size_t)comm_rank * num_indices_per_rank + i);

  Xt_idxlist src_idxlist
    = xt_idxvec_new(src_index_list, (int)num_indices_per_rank);
  free(src_index_list);

  // destination index list
  size_t num_gathered = (size_t)comm_size * num_indices_per_rank;
  Xt_int *dst_index_list = xmalloc(sizeof(Xt_int) * num_gathered);
  for (size_t i = 0; i < num_gathered; ++i)
    dst_index_list[i] = (Xt_int)i;

  Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index_list, (int)num_gathered);
  free(dst_index_list);

  // test of exchange map
  Xt_xmap xmap = xmap_new(src_idxlist, dst_idxlist, comm);

  // test results
  check_xmap_allgather_analog_xmap(xmap, comm);
  Xt_xmap xmap_copy = xt_xmap_copy(xmap);
  check_xmap_allgather_analog_xmap(xmap_copy, comm);

  // clean up
  xt_xmap_delete(xmap);
  xt_xmap_delete(xmap_copy);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);
}

static void
check_ring_xmap(Xt_xmap xmap, const Xt_int dst_index_list[], bool is_inter)
{
  int num_dst = xt_xmap_get_num_destinations(xmap),
    num_src = xt_xmap_get_num_sources(xmap);
  if (!is_inter && (num_dst > 2 || num_dst < 1))
    PUT_ERR("error in xt_xmap_get_num_destinations\n");

  if (num_src > 2 || num_src < 1)
    PUT_ERR("error in xt_xmap_get_num_sources\n");

  int ranks[2];

  if (!is_inter) {
    xt_xmap_get_destination_ranks(xmap, ranks);
    xt_sort_int(ranks, (size_t)num_dst);
    if (ranks[0] != dst_index_list[0] ||
        ranks[num_dst > 1] != dst_index_list[1])
      PUT_ERR("error in xt_xmap_get_destination_ranks\n");
  }

  xt_xmap_get_source_ranks(xmap, ranks);
  xt_sort_int(ranks, (size_t)num_src);
  if (ranks[0] != dst_index_list[0] ||
      ranks[num_src > 1] != dst_index_list[1])
    PUT_ERR("error in xt_xmap_get_source_ranks\n");
}

void
test_ring_1d(xmap_constructor xmap_new, MPI_Comm comm)
{
  int comm_rank, comm_size, is_inter;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
  int (*get_comm_size)(MPI_Comm, int *)
    = is_inter ? MPI_Comm_remote_size : MPI_Comm_size;
  xt_mpi_call(get_comm_size(comm, &comm_size), comm);
  // test in which each process talks with two other
  // processes

  Xt_int src_index_list[1] = {(Xt_int)(comm_rank)};

  Xt_idxlist src_idxlist = xt_idxvec_new(src_index_list, 1);

  // destination index list
  Xt_int dst_index_list[2] = {(Xt_int)((comm_rank + comm_size - 1)%comm_size),
                              (Xt_int)((comm_rank             + 1)%comm_size)};

  if (dst_index_list[0] > dst_index_list[1]) {
    Xt_int temp = dst_index_list[0];
    dst_index_list[0] = dst_index_list[1];
    dst_index_list[1] = temp;
  }

  Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index_list, 2);

  // test of exchange map
  Xt_xmap xmap = xmap_new(src_idxlist, dst_idxlist, comm);
  xt_idxlist_delete(src_idxlist);

  // test results
  check_ring_xmap(xmap, dst_index_list, (bool)is_inter);
  Xt_xmap xmap_copy = xt_xmap_copy(xmap);
  check_ring_xmap(xmap_copy, dst_index_list, (bool)is_inter);

  // clean up
  xt_xmap_delete(xmap);
  xt_xmap_delete(xmap_copy);
  xt_idxlist_delete(dst_idxlist);
}

void
test_maxpos(xmap_constructor xmap_new, MPI_Comm comm,
            int indices_per_rank)
{
  int comm_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  int world_size = comm_size * indices_per_rank;

  // setup
  Xt_int src_index[indices_per_rank];
  for (int i = 0; i < indices_per_rank; ++i)
    src_index[i] = (Xt_int)(i + comm_rank * indices_per_rank);

  Xt_int dst_index[indices_per_rank];
  for (int i = 0; i < indices_per_rank/2; ++i)
    dst_index[i]
      = (Xt_int)((i - indices_per_rank/2
                  + (comm_rank+comm_size) * indices_per_rank)%world_size);
  for (int i = indices_per_rank/2; i < (indices_per_rank+1)/2; ++i)
    dst_index[i] = (Xt_int)(i + comm_rank * indices_per_rank);
  for (int i = 0; i < indices_per_rank/2; ++i)
    dst_index[(indices_per_rank+1)/2+i]
      = (Xt_int)(((int)i + (comm_rank+1) * indices_per_rank)%world_size);
  Xt_idxlist src_idxlist = xt_idxvec_new(src_index, indices_per_rank);
  Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index, indices_per_rank);

  Xt_xmap xmap = xmap_new(src_idxlist, dst_idxlist, comm);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);

  // test
  // 1. test that initial max positions are in range
  int max_pos_dst = xt_xmap_get_max_dst_pos(xmap),
    max_pos_src = xt_xmap_get_max_src_pos(xmap);
  if (max_pos_dst < indices_per_rank-1)
    PUT_ERR("error in xt_xmap_get_max_dst_pos\n");
  if (max_pos_src < indices_per_rank-1)
    PUT_ERR("error in xt_xmap_get_max_src_pos\n");

  // 2. expand range and verify it is reflected in max pos
  int pos_update1[indices_per_rank];
  for (int i = 0; i < indices_per_rank; ++i)
    pos_update1[i] = 2*i;
  Xt_xmap xmup = xt_xmap_update_positions(xmap, pos_update1, pos_update1);

  int max_pos_dst_u = xt_xmap_get_max_dst_pos(xmup),
    max_pos_src_u = xt_xmap_get_max_src_pos(xmup);
  if (max_pos_dst_u < (indices_per_rank-1)*2)
    PUT_ERR("error in xt_xmap_get_max_dst_pos\n");
  if (max_pos_src_u < (indices_per_rank-1)*2)
    PUT_ERR("error in xt_xmap_get_max_src_pos\n");

  // 3. contract range again and verify max pos is updated
  int pos_update2[2*indices_per_rank];
  for (int i = 0; i < 2*indices_per_rank; ++i)
    pos_update2[i] = i/2;

  Xt_xmap xmup2 = xt_xmap_update_positions(xmup, pos_update2, pos_update2);

  int max_pos_dst_u2 = xt_xmap_get_max_dst_pos(xmup2),
    max_pos_src_u2 = xt_xmap_get_max_src_pos(xmup2);
  if (max_pos_dst_u2 >= indices_per_rank)
    PUT_ERR("error in xt_xmap_get_max_dst_pos\n");
  if (max_pos_src_u2 >= indices_per_rank)
    PUT_ERR("error in xt_xmap_get_max_src_pos\n");

  // 4. apply spread and check max pos range
  int spread[] = { 0, indices_per_rank*3 };
  Xt_xmap xmsp = xt_xmap_spread(xmap, 2, spread, spread);
  int max_pos_dst_s = xt_xmap_get_max_dst_pos(xmsp),
    max_pos_src_s = xt_xmap_get_max_src_pos(xmsp);
  if (max_pos_dst_s < (indices_per_rank-1)*3)
    PUT_ERR("error in xt_xmap_get_max_dst_pos\n");
  if (max_pos_src_s < (indices_per_rank-1)*3)
    PUT_ERR("error in xt_xmap_get_max_src_pos\n");

  // cleanup
  xt_xmap_delete(xmap);
  xt_xmap_delete(xmup);
  xt_xmap_delete(xmup2);
  xt_xmap_delete(xmsp);
}


static void
check_pair_xmap(Xt_xmap xmap)
{
  if (xt_xmap_get_num_destinations(xmap) != 2)
    PUT_ERR("error in xt_xmap_get_num_destinations\n");

  if (xt_xmap_get_num_sources(xmap) != 2)
    PUT_ERR("error in xt_xmap_get_num_sources\n");

  int ranks[2];

  xt_xmap_get_destination_ranks(xmap, ranks);
  if (ranks[0] != 0 || ranks[1] != 1)
    PUT_ERR("error in xt_xmap_get_destination_ranks\n");

  xt_xmap_get_source_ranks(xmap, ranks);
  if (ranks[0] != 0 || ranks[1] != 1)
    PUT_ERR("error in xt_xmap_get_source_ranks\n");
}


static void
test_pair(xmap_constructor xmap_new, MPI_Comm comm)
{
  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  enum { numIdx = 20 };
  static const Xt_int src_index_list[2][numIdx] = { //src_index_list[rank][index];
    {1,2,3,4,5,9,10,11,12,13,17,18,19,20,21,25,26,27,28,29},
    {4,5,6,7,8,12,13,14,15,16,20,21,22,23,24,28,29,30,31,32} };

  Xt_idxlist src_idxlist = xt_idxvec_new(src_index_list[comm_rank], numIdx);

  // destination index list
  static const Xt_int dst_index_list[2][numIdx] = { //dst_index_list[rank][index];
    {10,15,14,13,12,15,10,11,12,13,23,18,19,20,21,31,26,27,28,29},
    {13,12,11,10,15,12,13,14,15,10,20,21,22,23,18,28,29,30,31,26}};

  Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index_list[comm_rank], numIdx);


  // test of exchange map
  Xt_xmap xmap = xmap_new(src_idxlist, dst_idxlist, comm);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);

  // test results

  check_pair_xmap(xmap);
  Xt_xmap xmap_copy = xt_xmap_copy(xmap);
  check_pair_xmap(xmap_copy);
  // clean up
  xt_xmap_delete(xmap);
  xt_xmap_delete(xmap_copy);
}

static void
test_ping_pong_xmap(Xt_xmap xmap, MPI_Comm comm, int ping_rank, int pong_rank)
{
  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);

  if (xt_xmap_get_num_destinations(xmap) != (comm_rank == ping_rank))
    PUT_ERR("error in xt_xmap_get_num_destinations (rank == %d)\n", comm_rank);

  if (xt_xmap_get_num_sources(xmap) != (comm_rank == pong_rank))
    PUT_ERR("error in xt_xmap_get_num_sources (rank == %d)\n", comm_rank);

  if (comm_rank == ping_rank) {
    int dst_rank;
    xt_xmap_get_destination_ranks(xmap, &dst_rank);
    if (dst_rank != pong_rank)
      PUT_ERR("error in xt_xmap_get_destination_ranks\n");
  }

  if (comm_rank == pong_rank) {
    int src_rank;
    xt_xmap_get_source_ranks(xmap, &src_rank);
    if (src_rank != ping_rank)
      PUT_ERR("error in xt_xmap_get_source_ranks\n");
  }
}

void
test_ping_pong(xmap_constructor xmap_new, MPI_Comm comm,
               int ping_rank, int pong_rank)
{
  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);

  enum { numIdx = 5 };
  static const Xt_int index_list[numIdx] = {0,1,2,3,4};

  Xt_idxlist src_idxlist = (comm_rank == ping_rank)?
    xt_idxvec_new(index_list, numIdx) : xt_idxempty_new();

  Xt_idxlist dst_idxlist = (comm_rank == pong_rank)?
    xt_idxvec_new(index_list, numIdx) : xt_idxempty_new();

  // test of exchange map
  Xt_xmap xmap = xmap_new(src_idxlist, dst_idxlist, comm);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);

  // test results
  test_ping_pong_xmap(xmap, comm, ping_rank, pong_rank);
  Xt_xmap xmap_copy = xt_xmap_copy(xmap);
  test_ping_pong_xmap(xmap_copy, comm, ping_rank, pong_rank);

  // clean up
  xt_xmap_delete(xmap);
  xt_xmap_delete(xmap_copy);
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
