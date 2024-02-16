/**
 * @file test_xmap_common_intercomm_parallel.c
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

#include <mpi.h>
#include <yaxt.h>

#include "tests.h"
#include "test_xmap_common.h"
#include "core/ppm_xfuncs.h"

static void
test_xmap_allgather_analog(xmap_constructor xmap_new,
                           size_t num_indices_per_rank, MPI_Comm comm);
static void
test_peer(xmap_constructor xmap_new, MPI_Comm comm);

MPI_Comm xt_intra_group_comm;

int
xt_xmap_intercomm_parallel_test_main(xmap_constructor xmap_new,
                                     bool call_initialize, bool call_finalize)
{
  // init mpi
  if (call_initialize) {
    xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
    xt_initialize(MPI_COMM_WORLD);
  }
  MPI_Comm comm = MPI_COMM_WORLD;
  int comm_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int retval;
  if (comm_size > 1) {
    /* split communicator in half to create intercommunicator */
    MPI_Comm inter_comm;
    /* favor inequal splits if the left group does not become too small */
    int splitRank = comm_size/2 - (comm_size > 5);
    int inSecondGroup = comm_rank >= splitRank;
    xt_mpi_call(MPI_Comm_split(comm, inSecondGroup, 0, &xt_intra_group_comm), comm);
    xt_mpi_call(MPI_Intercomm_create(xt_intra_group_comm, 0, comm,
                                     inSecondGroup ? 0 : splitRank,
                                     0, &inter_comm), comm);

    test_xmap_allgather_analog(xmap_new, 1, inter_comm);
    /* repeat test for large index list that will cause stripifying */
    test_xmap_allgather_analog(xmap_new, 1024, inter_comm);

    test_peer(xmap_new, inter_comm);

    test_ping_pong(xmap_new, inter_comm, 0, 0);

    if (splitRank > 2) // skip test if there are not enough processes
      test_ring_1d(xmap_new, inter_comm);

    xt_mpi_call(MPI_Comm_free(&inter_comm), comm);
    xt_mpi_call(MPI_Comm_free(&xt_intra_group_comm), comm);
    retval = TEST_EXIT_CODE;
  } else
    retval = 77;

  if (call_finalize) {
    xt_finalize();
    MPI_Finalize();
  }
  return retval;
}

static void
test_xmap_allgather_analog(xmap_constructor xmap_new,
                           size_t num_indices_per_rank,
                           MPI_Comm comm)
{
  int comm_rank, remote_size;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_remote_size(comm, &remote_size), comm);

  // test in which every process requests data from all processes
  // source index list
  struct Xt_stripe src_index_stripe = {
    .start = (Xt_int)((Xt_int)comm_rank * (Xt_int)num_indices_per_rank),
    .stride = 1,
    .nstrides = (int)num_indices_per_rank
  };
  Xt_idxlist src_idxlist = xt_idxstripes_new(&src_index_stripe, 1);

  // destination index list
  size_t num_gathered = (size_t)remote_size * num_indices_per_rank;
  struct Xt_stripe dst_index_stripe = {
    .start = 0,
    .stride = 1,
    .nstrides = (int)num_gathered
  };
  Xt_idxlist dst_idxlist = xt_idxstripes_new(&dst_index_stripe, 1);

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


typedef int (*xmap_get_num_peers)(Xt_xmap xmap);
typedef void (*xmap_get_peer_ranks)(Xt_xmap xmap, int *ranks);

static void
check_peers(Xt_xmap xmap, size_t num_ref_ranks, const int *ref_ranks,
            int *peer_rank_buf, xmap_get_num_peers get_num_peers,
            xmap_get_peer_ranks get_peer_ranks,
            const char *get_num_peers_name,
            const char *get_peer_ranks_name)
{
  (void)get_num_peers_name;
  (void)get_peer_ranks_name;
  if (get_num_peers(xmap) != (int)num_ref_ranks)
    PUT_ERR("error in %s\n", get_num_peers_name);
  get_peer_ranks(xmap, peer_rank_buf);
  xt_sort_int(peer_rank_buf, num_ref_ranks);
  bool mismatch = false;
  for (size_t i = 0; i < num_ref_ranks; ++i)
    mismatch |= (peer_rank_buf[i] != ref_ranks[i]);
  if (mismatch)
    PUT_ERR("error in %s\n", get_peer_ranks_name);
}

static void
check_peer_xmap(Xt_xmap xmap, struct Xt_stripe stripe_in_local_group,
                int remote_size, long long global_num_idx)
{
  size_t num_indices = (size_t)stripe_in_local_group.nstrides;
  int idx_per_remote_rank = (int)(global_num_idx / remote_size);
  int *ref_ranks = xmalloc(sizeof (*ref_ranks) * num_indices);
  size_t num_remote_peers = 0;
  int last_seen_rank = -1;
  Xt_int start_idx = stripe_in_local_group.start;
  if (num_indices > 0) {
    size_t i = 0;
    for (;;) {
      int remote_rank_i = (int)((start_idx + (Xt_int)i) / idx_per_remote_rank);
      while (i < num_indices && remote_rank_i == last_seen_rank)
        remote_rank_i = (int)((start_idx + (Xt_int)++i) / idx_per_remote_rank);
      if (i >= num_indices)
        goto num_remote_peers_established;
      ref_ranks[num_remote_peers] = remote_rank_i;
      last_seen_rank = remote_rank_i;
      ++num_remote_peers;
      ++i;
    }
  }
num_remote_peers_established: ;
  ref_ranks = xrealloc(ref_ranks, sizeof (*ref_ranks) * num_remote_peers * 2);
  int *rank_buf = ref_ranks + num_remote_peers;
  check_peers(xmap, num_remote_peers, ref_ranks, rank_buf,
              xt_xmap_get_num_destinations, xt_xmap_get_destination_ranks,
              "xt_xmap_get_num_destinations", "xt_xmap_get_destination_ranks");
  check_peers(xmap, num_remote_peers, ref_ranks, rank_buf,
              xt_xmap_get_num_sources, xt_xmap_get_source_ranks,
              "xt_xmap_get_num_sources", "xt_xmap_get_source_ranks");
  free(ref_ranks);
}

static int gcd(int a, int b)
{
  while (b != 0) {
    int t = b;
    b = a % b;
    a = t;
  }
  return a;
}

static long long lcm(int a, int b)
{
  int t = gcd(a, b);
  return (long long)(a / t) * b;
}

static void
test_peer(xmap_constructor xmap_new, MPI_Comm comm)
{
  int comm_rank, comm_size, remote_size;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  xt_mpi_call(MPI_Comm_remote_size(comm, &remote_size), comm);

  long long global_num_idx = lcm(comm_size, remote_size);
  struct Xt_stripe stripe_in_local_group = {
    .start = (Xt_int)(global_num_idx / comm_size * comm_rank),
    .stride = 1,
    .nstrides = (int)(global_num_idx / comm_size)
  };
  Xt_idxlist idxlist = xt_idxstripes_new(&stripe_in_local_group, 1);

  // test of exchange map
  Xt_xmap xmap = xmap_new(idxlist, idxlist, comm);
  xt_idxlist_delete(idxlist);

  // test results

  check_peer_xmap(xmap, stripe_in_local_group, remote_size,
                  global_num_idx);
  Xt_xmap xmap_copy = xt_xmap_copy(xmap);
  check_peer_xmap(xmap_copy, stripe_in_local_group, remote_size,
                  global_num_idx);
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
