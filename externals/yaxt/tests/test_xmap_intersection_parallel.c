/**
 * @file test_xmap_intersection_parallel.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
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
#include <unistd.h>
#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "xt/xt_xmap_intersection.h"
#include "xt/xt_xmap.h"
#include "core/ppm_xfuncs.h"

struct test_message {

  int rank;       // rank of communication partner
  int num_pos; // number of positions
  const int *pos;   // positions to be sent/received
};

static void test_xmap(
  Xt_xmap xmap,
  int num_sends, const struct test_message send_messages[num_sends],
  int num_recvs, const struct test_message recv_messages[num_recvs]);

static Xt_xmap (*xmi_new)(
  int nsrc_com, const struct Xt_com_list src_com[nsrc_com],
  int ndst_com, const struct Xt_com_list dst_com[ndst_com],
  Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist, MPI_Comm comm);

static void
test_strided_block_pos_alltoall(MPI_Comm comm, int nblk, int blksz);

static void
parse_options(int *argc, char ***argv);

int main(int argc, char **argv)
{

  // init mpi
  xt_mpi_call(MPI_Init(&argc, &argv), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  int my_rank, comm_size;

  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  if (comm_size != 3) {

    xt_finalize();
    MPI_Finalize();
    return 77;
  }

  parse_options(&argc, &argv);

  { // simple test (round robin)

    // setup

    Xt_int src_index = (Xt_int)my_rank;
    Xt_int dst_index = (Xt_int)((my_rank + 1)%comm_size);
    Xt_idxlist src_idxlist = xt_idxvec_new(&src_index, 1);
    Xt_idxlist dst_idxlist = xt_idxvec_new(&dst_index, 1);
    int num_src_intersections = 1;
    struct Xt_com_list src_com = {.list = src_idxlist,
                                  .rank = (my_rank+1)%comm_size};
    int num_dst_intersections = 1;
    struct Xt_com_list dst_com = {.list = dst_idxlist,
                                  .rank = (my_rank+comm_size-1)%comm_size};

    Xt_xmap xmap = xmi_new(num_src_intersections, &src_com,
                           num_dst_intersections, &dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test

    int send_pos = 0;
    int num_sends = 1;
    struct test_message send_messages[1] = {{.rank = (my_rank+1)%comm_size,
                                             .pos = &send_pos, .num_pos = 1}};
    int recv_pos = 0;
    int num_recvs = 1;
    struct test_message recv_messages[1] = {{.rank = (my_rank+comm_size-1)%comm_size,
                                             .pos = &recv_pos, .num_pos = 1}};

    test_xmap(xmap, num_sends, send_messages, num_recvs, recv_messages);

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
  }

  { // rank 0 receives the same point from rank 1 and 2

    // setup

    Xt_int src_index = 0;
    Xt_int dst_index = 0;
    Xt_idxlist src_idxlist = (my_rank == 0)?xt_idxempty_new():xt_idxvec_new(&src_index, 1);
    Xt_idxlist dst_idxlist = (my_rank == 0)?xt_idxvec_new(&dst_index, 1):xt_idxempty_new();
    int num_src_intersections = (my_rank == 0)?0:1;
    struct Xt_com_list src_com = {.list = src_idxlist,
                                  .rank = 0};
    int num_dst_intersections = (my_rank == 0)?2:0;
    struct Xt_com_list dst_com[2] = {{.list = dst_idxlist, .rank = 1},
                                     {.list = dst_idxlist, .rank = 2}};

    Xt_xmap xmap = xmi_new(num_src_intersections, &src_com,
                           num_dst_intersections, dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test

    int send_pos = 0;
    int num_sends = (my_rank == 1)?1:0;
    struct test_message send_messages[1] = {{.rank = 0, .pos = &send_pos,
                                             .num_pos = 1}};
    int recv_pos = 0;
    int num_recvs = (my_rank == 0)?1:0;
    struct test_message recv_messages[1] = {{.rank = 1, .pos = &recv_pos,
                                             .num_pos = 1}};

    test_xmap(xmap, num_sends, send_messages, num_recvs, recv_messages);

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
  }

  { // all ranks can receive data from the others
    // rank               |  0  |  1  |  2  |
    // source indices     | 1,2 | 2,0 | 0,1 |
    // destination indice |  0  |  1  |  2  |

    // setup

    Xt_int src_indices[2] = { (Xt_int)((my_rank+1)%comm_size),
                              (Xt_int)((my_rank+2)%comm_size) };
    Xt_int dst_index = (Xt_int)my_rank;
    Xt_idxlist src_idxlist = xt_idxvec_new(src_indices, 2);
    Xt_idxlist dst_idxlist = xt_idxvec_new(&dst_index, 1);
    int num_src_intersections[3] = {2, 1, 0};
    Xt_idxlist src_intersection_idxlist[2] = {xt_idxvec_new(src_indices, 1),
                                              xt_idxvec_new(src_indices+1, 1)};
    struct Xt_com_list src_com[2] = {{.list = src_intersection_idxlist[0],
                                      .rank = 1},
                                     {.list = src_intersection_idxlist[1],
                                      .rank = (my_rank == 0)?2:0}};
    int num_dst_intersections = 1;
    struct Xt_com_list dst_com = {.list = dst_idxlist, .rank = (my_rank == 0)?1:0};

    Xt_xmap xmap = xmi_new(num_src_intersections[my_rank], src_com + my_rank,
                           num_dst_intersections, &dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test

    if (my_rank == 0) {

      int send_pos[2] = {0,1};
      int num_sends = 2;
      struct test_message send_messages[2] = {{.rank = 1, .pos = send_pos+0,
                                               .num_pos = 1},
                                              {.rank = 2, .pos = send_pos+1,
                                               .num_pos = 1}};
      int recv_pos = 0;
      int num_recvs = 1;
      struct test_message recv_messages[1] = {{.rank = 1, .pos = &recv_pos,
                                               .num_pos = 1}};

      test_xmap(xmap, num_sends, send_messages, num_recvs, recv_messages);

    } else if (my_rank == 1) {

      int send_pos = 1;
      int num_sends = 1;
      struct test_message send_messages[1] = {{.rank = 0, .pos = &send_pos,
                                               .num_pos = 1}};
      int recv_pos = 0;
      int num_recvs = 1;
      struct test_message recv_messages[1] = {{.rank = 0, .pos = &recv_pos,
                                               .num_pos = 1}};

      test_xmap(xmap, num_sends, send_messages, num_recvs, recv_messages);

    } else {

      int num_sends = 0;

      int recv_pos = 0;
      int num_recvs = 1;
      struct test_message recv_messages[1] = {{.rank = 0, .pos = &recv_pos,
                                               .num_pos = 1}};

      test_xmap(xmap, num_sends, NULL, num_recvs, recv_messages);
    }

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_intersection_idxlist[1]);
    xt_idxlist_delete(src_intersection_idxlist[0]);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
  }

  { // all ranks can receive data from the others
    // rank               |         0         |         1         |         2         |
    // source indices     |     0,1,2,3,4     |     3,4,5,6,7     |     6,7,8,0,1     |
    // destination indice | 0,1,2,3,4,5,6,7,8 | 0,1,2,3,4,5,6,7,8 | 0,1,2,3,4,5,6,7,8 |

    // setup

    static const Xt_int src_indices[3][5]
      = {{0,1,2,3,4}, {3,4,5,6,7}, {6,7,8,0,1}};
    static const Xt_int dst_indices[9] = {0,1,2,3,4,5,6,7,8};
    Xt_idxlist src_idxlist = xt_idxvec_new(src_indices[my_rank], 5);
    Xt_idxlist dst_idxlist = xt_idxvec_new(dst_indices, 9);

    struct Xt_com_list src_com[3] = {{.list = src_idxlist, .rank = 0},
                                     {.list = src_idxlist, .rank = 1},
                                     {.list = src_idxlist, .rank = 2}};
    struct Xt_com_list dst_com[3] =
      {{.list = xt_idxvec_new(src_indices[0], 5), .rank = 0},
       {.list = xt_idxvec_new(src_indices[1], 5), .rank = 1},
       {.list = xt_idxvec_new(src_indices[2], 5), .rank = 2}};

    Xt_xmap xmap = xmi_new(3, src_com, 3, dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test

    static const int send_pos[3][5] = {{0,1,2,3,4}, {2,3,4}, {2}};
    static const int num_send_pos[3] = {5, 3, 1};
    struct test_message send_messages[3] = {{.rank = 0,
                                             .pos = send_pos[my_rank],
                                             .num_pos = num_send_pos[my_rank]},
                                            {.rank = 1,
                                             .pos = send_pos[my_rank],
                                             .num_pos = num_send_pos[my_rank]},
                                            {.rank = 2,
                                             .pos = send_pos[my_rank],
                                             .num_pos = num_send_pos[my_rank]}};
    static const int recv_pos[3][5] = {{0,1,2,3,4}, {5,6,7}, {8}};
    static const int num_recv_pos[3] = {5, 3, 1};
    struct test_message recv_messages[3] = {{.rank = 0, .pos = recv_pos[0],
                                             .num_pos = num_recv_pos[0]},
                                            {.rank = 1, .pos = recv_pos[1],
                                             .num_pos = num_recv_pos[1]},
                                            {.rank = 2, .pos = recv_pos[2],
                                             .num_pos = num_recv_pos[2]}};

    test_xmap(xmap, 3, send_messages, 3, recv_messages);

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_com[2].list);
    xt_idxlist_delete(dst_com[1].list);
    xt_idxlist_delete(dst_com[0].list);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
  }

  { // one rank receives data from the other two, that have duplicated indices
    // (this provokes a bug found by Joerg Behrens)
    // rank               |  0  |  1  |   2   |
    // source indices     | 0,2 | 1,2 |       |
    // destination indice |     |     | 0,1,2 |

    // setup

    Xt_int src_indices[2][2] = {{0,2}, {1,2}};

    struct Xt_com_list * src_com, * dst_com;
    int num_src_intersections, num_dst_intersections;
    Xt_idxlist src_idxlist, dst_idxlist;

    if (my_rank == 2) {

      src_com = NULL;
      num_src_intersections = 0;

      dst_com = malloc(2 * sizeof(*dst_com));
      num_dst_intersections = 2;
      dst_com[0].list = xt_idxvec_new(src_indices[0], 2);
      dst_com[0].rank = 0;
      dst_com[1].list = xt_idxvec_new(src_indices[1], 2);
      dst_com[1].rank = 1;

      Xt_int dst_indices[3] = {0,1,2};

      src_idxlist = xt_idxempty_new();
      dst_idxlist = xt_idxvec_new(dst_indices, 3);

    } else {

      src_com = malloc(1 * sizeof(*src_com));
      src_com->list = xt_idxvec_new(src_indices[my_rank], 2);
      src_com->rank = 2;
      num_src_intersections = 1;

      dst_com = NULL;
      num_dst_intersections = 0;

      src_idxlist = xt_idxvec_new(src_indices[my_rank], 2);
      dst_idxlist = xt_idxempty_new();
    }

    Xt_xmap xmap = xmi_new(num_src_intersections, src_com,
                           num_dst_intersections, dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test

      if (my_rank == 2) {

        static const int recv_pos[2][2] = {{0,2}, {1}};
        static const struct test_message recv_messages[2]
          = {{.rank = 0, .pos = recv_pos[0], .num_pos = 2},
             {.rank = 1, .pos = recv_pos[1], .num_pos = 1}};

        test_xmap(xmap, 0, NULL, 2, recv_messages);

      } else {

        static const int send_pos[2][2] = {{0,1}, {0}};
        static const int num_send_pos[2] = {2, 1};
        struct test_message send_messages = {.rank = 2,
                                             .pos = send_pos[my_rank],
                                             .num_pos = num_send_pos[my_rank]};

        test_xmap(xmap, 1, &send_messages, 0, NULL);
      }

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
    for (int i = 0; i < num_dst_intersections; ++i)
      xt_idxlist_delete(dst_com[i].list);
    free(dst_com);
    for (int i = 0; i < num_src_intersections; ++i)
      xt_idxlist_delete(src_com[i].list);
    free(src_com);
  }

  { // checks the reorder functionality of exchange maps

    // setup

    enum { count = 6 };
    Xt_int src_indices[count] = {0,5,1,4,2,3};
    Xt_int dst_indices[count] = {5,4,3,2,1,0};
    Xt_int intersection_indices[count] = {0,1,2,3,4,5};

    struct Xt_com_list src_com =
      {.list = xt_idxvec_new(intersection_indices, count),
       .rank = (my_rank + 1) % comm_size};
    struct Xt_com_list dst_com =
      {.list = xt_idxvec_new(intersection_indices, count),
       .rank = (comm_size + my_rank - 1) % comm_size};

    Xt_idxlist src_idxlist = xt_idxvec_new(src_indices, count);
    Xt_idxlist dst_idxlist = xt_idxvec_new(dst_indices, count);

    Xt_xmap xmap = xmi_new(1, &src_com, 1, &dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    enum xt_reorder_type reorder_types[] =
      {XT_REORDER_NONE, XT_REORDER_SEND_UP, XT_REORDER_RECV_UP};

    // test

    for (size_t i = 0; i < sizeof(reorder_types) / sizeof(reorder_types[0]);
         ++i) {

      enum xt_reorder_type curr_reorder_type = reorder_types[i];

      Xt_xmap xmap_reorder = xt_xmap_reorder(xmap, curr_reorder_type);

      struct test_message send_messages =
        {.pos = (int[]){0,2,4,5,3,1}, .num_pos = 6};
      send_messages.rank = (my_rank + 1) % comm_size;
      struct test_message recv_messages =
        {.pos = (int[]){5,4,3,2,1,0}, .num_pos = 6};
      recv_messages.rank = (comm_size + my_rank - 1) % comm_size;

      switch((int)curr_reorder_type) {
        case(XT_REORDER_SEND_UP):
          xt_quicksort_int_permutation(
            (int*)send_messages.pos, (size_t)count, (int*)recv_messages.pos);
          break;
        case(XT_REORDER_RECV_UP):
          xt_quicksort_int_permutation(
            (int*)recv_messages.pos, (size_t)count, (int*)send_messages.pos);
          break;
        default:
          break;
      };

      test_xmap(xmap_reorder, 1, &send_messages, 1, &recv_messages);

      xt_xmap_delete(xmap_reorder);
    }

    // cleanup

    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_com.list);
    xt_idxlist_delete(src_com.list);
  }

  { // checks the spread and update positions functionality of exchange maps

    // setup

    enum { count = 12, nlev = 6, nproma = 4, nblk = 3};
    static const Xt_int src_indices[count] = {0,1,2,3,4,5,6,7,8,9,10,11};
    static const Xt_int dst_indices[count] = {0,1,2,3,4,5,6,7,8,9,10,11};
    static const Xt_int intersection_indices[count] = {0,1,2,3,4,5,6,7,8,9,10,11};
    int displacements[nlev];

    for (int i = 0; i < (int)nlev; ++i) displacements[i] = i * nproma;

    struct Xt_com_list src_com =
      {.list = xt_idxvec_new(intersection_indices, count),
       .rank = (my_rank + 1) % comm_size};
    struct Xt_com_list dst_com =
      {.list = xt_idxvec_new(intersection_indices, count),
       .rank = (comm_size + my_rank - 1) % comm_size};

    Xt_idxlist src_idxlist = xt_idxvec_new(src_indices, count);
    Xt_idxlist dst_idxlist = xt_idxvec_new(dst_indices, count);

    Xt_xmap xmap = xmi_new(1, &src_com, 1, &dst_com,
                           src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    int blocked_positions[(int)count * (int)nlev];
    for (int i = 0, lev = 0; lev < (int)nlev; ++lev)
      for (int blk = 0; blk < nblk; ++blk)
        for (int idx = 0; idx < nproma; ++idx, ++i)
          blocked_positions[i] = idx + (blk * nlev + lev) * nproma;

    Xt_xmap xmap_single_level_blocked =
      xt_xmap_update_positions(xmap, blocked_positions, blocked_positions);
    Xt_xmap xmap_multi_level_blocked =
      xt_xmap_spread(
        xmap_single_level_blocked, (int)nlev, displacements, displacements);

    struct test_message send_messages =
      {.pos = &(blocked_positions[0]), .num_pos = (int)nlev * (int)count};
    send_messages.rank = (my_rank + 1) % comm_size;
    struct test_message recv_messages =
      {.pos = &(blocked_positions[0]), .num_pos = (int)nlev * (int)count};
    recv_messages.rank = (comm_size + my_rank - 1) % comm_size;

    // test

    test_xmap(xmap_multi_level_blocked, 1, &send_messages, 1, &recv_messages);

    // cleanup

    xt_xmap_delete(xmap_multi_level_blocked);
    xt_xmap_delete(xmap_single_level_blocked);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_com.list);
    xt_idxlist_delete(src_com.list);
  }

  test_strided_block_pos_alltoall(MPI_COMM_WORLD, 1, 1);
  test_strided_block_pos_alltoall(MPI_COMM_WORLD, 2, 1);
  test_strided_block_pos_alltoall(MPI_COMM_WORLD, 5, 2000);

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
test_strided_block_pos_alltoall(MPI_Comm comm, int nblk, int blksz)
{
  int comm_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  struct Xt_com_pos *src_com =
    xmalloc(2 * (size_t)comm_size * sizeof(*src_com)),
    *dst_com = src_com + comm_size;
  int (*transfer_pos)[nblk][blksz]
    = xmalloc((size_t)nblk * (size_t)blksz
              * (size_t)comm_size * sizeof(*transfer_pos));
  struct test_message *ref_src_msg =
    xmalloc(2 * (size_t)comm_size * sizeof(*ref_src_msg)),
    *ref_dst_msg = ref_src_msg + comm_size;

  for (int rank = 0; rank < comm_size; ++rank) {
    for (int j = 0; j < nblk; ++j)
      for (int i = 0; i < blksz; ++i)
        transfer_pos[rank][j][i] = i + j*blksz*comm_size;

#if defined __PGI && __PGIC__ <= 15
    ref_src_msg[rank].pos     = src_com[rank].transfer_pos
      = (int *)transfer_pos + rank * nblk * blksz;
#else
    ref_src_msg[rank].pos     = src_com[rank].transfer_pos
      = (int *)(transfer_pos + rank);
#endif
    ref_src_msg[rank].num_pos = src_com[rank].num_transfer_pos = nblk*blksz;
    ref_src_msg[rank].rank    = src_com[rank].rank             = rank;
#if defined __PGI && __PGIC__ <= 15
    ref_dst_msg[rank].pos     = dst_com[rank].transfer_pos
      = (int *)transfer_pos + rank * nblk * blksz;
#else
    ref_dst_msg[rank].pos     = dst_com[rank].transfer_pos
      = (int *)(transfer_pos + rank);
#endif
    ref_dst_msg[rank].num_pos = dst_com[rank].num_transfer_pos = nblk*blksz;
    ref_dst_msg[rank].rank    = dst_com[rank].rank             = rank;
  }

  Xt_xmap xmap =
    xt_xmap_intersection_pos_new(comm_size, src_com, comm_size, dst_com, comm);

  test_xmap(xmap, comm_size, ref_src_msg, comm_size, ref_dst_msg);

  xt_xmap_delete(xmap);

  free(ref_src_msg);
  free(transfer_pos);
  free(src_com);
}

static void test_xmap_iter(Xt_xmap_iter iter, int num_msgs,
                           const struct test_message msgs[num_msgs]) {

  if (num_msgs == 0) {

    if (iter != NULL)
      PUT_ERR("ERROR: xt_xmap_get_*_iterator (iter should be NULL)\n");

  } else if (iter == NULL) {

    PUT_ERR("ERROR: xt_xmap_get_*_iterator (iter should not be NULL)\n");

  } else {

    int i = 0;

    do {

      if (xt_xmap_iterator_get_rank(iter) != msgs[i].rank)
        PUT_ERR("ERROR: xt_xmap_iterator_get_rank\n");

      if (xt_xmap_iterator_get_num_transfer_pos(iter) != msgs[i].num_pos)
        PUT_ERR("ERROR: xt_xmap_iterator_get_num_transfer_pos\n");

      const int *restrict pos = xt_xmap_iterator_get_transfer_pos(iter);

      bool mismatch = false;
      for (int j = 0; j < msgs[i].num_pos; ++j)
        mismatch |= (pos[j] != msgs[i].pos[j]);
      if (mismatch)
        PUT_ERR("ERROR: xt_xmap_iterator_get_transfer_pos\n");

      int num_transfer_pos_ext
        = xt_xmap_iterator_get_num_transfer_pos_ext(iter);
      const struct Xt_pos_ext *restrict pos_ext
        = xt_xmap_iterator_get_transfer_pos_ext(iter);

      mismatch = false;
      size_t ofs = 0;
      for (size_t pe = 0; pe < (size_t)num_transfer_pos_ext; ++pe) {
        int pos_ext_size = abs(pos_ext[pe].size);
        if (pos_ext[pe].size > 0)
          for (int j = 0; j < pos_ext_size; ++j)
            mismatch |= (pos[ofs+(size_t)j] != pos_ext[pe].start + j);
        else
          for (int j = 0; j < pos_ext_size; ++j)
            mismatch |= (pos[ofs+(size_t)j] != pos_ext[pe].start - j);
        ofs += (size_t)pos_ext_size;
      }
      if (mismatch || (int)ofs != msgs[i].num_pos)
        PUT_ERR("ERROR: xt_xmap_iterator_get_transfer_pos_ext\n");

      ++i;

    } while (xt_xmap_iterator_next(iter));

    if (i != num_msgs)
      PUT_ERR("ERROR: xt_xmap_iterator_next (wrong number of message)\n");
  }
}

static void test_xmap(
  Xt_xmap xmap,
  int num_sends, const struct test_message send_messages[num_sends],
  int num_recvs, const struct test_message recv_messages[num_recvs]) {

  enum { numXmaps2Test = 2 };
  Xt_xmap maps[numXmaps2Test] = { xmap, xt_xmap_copy(xmap) };
  for (size_t i = 0; i < numXmaps2Test; ++i) {
    if (xt_xmap_get_num_destinations(maps[i]) != num_sends)
      PUT_ERR("ERROR: xt_xmap_get_num_destinations\n");
    if (xt_xmap_get_num_sources(maps[i]) != num_recvs)
      PUT_ERR("ERROR: xt_xmap_get_num_sources\n");

    Xt_xmap_iter send_iter = xt_xmap_get_out_iterator(maps[i]);
    Xt_xmap_iter recv_iter = xt_xmap_get_in_iterator(maps[i]);

    test_xmap_iter(send_iter, num_sends, send_messages);
    test_xmap_iter(recv_iter, num_recvs, recv_messages);

    if (recv_iter != NULL) xt_xmap_iterator_delete(recv_iter);
    if (send_iter != NULL) xt_xmap_iterator_delete(send_iter);
  }
  xt_xmap_delete(maps[1]);
}

static void
parse_options(int *argc, char ***argv)
{
  xmi_new = xt_xmap_intersection_new;
  int opt;
  while ((opt = getopt(*argc, *argv, "m:")) != -1) {
    switch (opt) {
    case 'm':
      if (!strcmp(optarg, "xt_xmap_intersection_new"))
        xmi_new = xt_xmap_intersection_new;
      else if (!strcmp(optarg, "xt_xmap_intersection_ext_new"))
        xmi_new = xt_xmap_intersection_ext_new;
      else
      {
        fprintf(stderr, "Unknown xmap intersection constructor requested %s\n",
                optarg);
        exit(EXIT_FAILURE);
      }
    }
  }
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
