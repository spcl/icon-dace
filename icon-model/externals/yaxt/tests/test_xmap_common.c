/**
 * @file test_xmap_common.c
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

#include <mpi.h>

#include <yaxt.h>

#define VERBOSE
#include "tests.h"
#include "ctest_common.h"
#include "test_xmap_common.h"

static void
test_xmap1(xmap_constructor new_xmap,
           int lsize, MPI_Comm comm);
static void
test_xmap2(xmap_constructor new_xmap, MPI_Comm comm);

static int my_rank;

int
xt_xmap_self_test_main(int *argc, char ***argv,
                       xmap_constructor new_xmap)
{
  test_init_mpi(argc, argv, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  MPI_Comm comms[2];
  comms[0] = MPI_COMM_WORLD;
  MPI_Comm_dup(MPI_COMM_WORLD, &comms[1]);
  xt_mpi_comm_mark_exclusive(comms[1]);

  for (size_t i = 0; i < 2; ++i) {
    static const int lsizes[2] = { 7, 1023 };
    for (size_t j = 0; j < 2; ++j)
      test_xmap1(new_xmap, lsizes[j], comms[i]);
    test_xmap2(new_xmap, comms[i]);
  }

  MPI_Comm_free(&comms[1]);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
assert_xmap_is_to_self(Xt_xmap xmap)
{
  // test results
  if (xt_xmap_get_num_destinations(xmap) != 1)
    PUT_ERR("error in xt_xmap_get_num_destinations\n");

  if (xt_xmap_get_num_sources(xmap) != 1)
    PUT_ERR("error in xt_xmap_get_num_sources\n");

  int rank;

  xt_xmap_get_destination_ranks(xmap, &rank);
  if (rank != my_rank)
    PUT_ERR("error in xt_xmap_get_destination_ranks\n");

  xt_xmap_get_source_ranks(xmap, &rank);
  if (rank != my_rank)
    PUT_ERR("error in xt_xmap_get_source_ranks\n");
}


void
test_self_xmap_construct_idxvec(
  const Xt_int *src_index_list, int num_src_indices,
  const Xt_int *dst_index_list, int num_dst_indices,
  xmap_constructor new_xmap, MPI_Comm comm)
{
  Xt_idxlist src_idxlist = xt_idxvec_new(src_index_list, num_src_indices),
    dst_idxlist = xt_idxvec_new(dst_index_list, num_dst_indices);
  test_self_xmap_construct(src_idxlist, dst_idxlist, new_xmap, comm);
}


void
test_self_xmap_construct_idxstripes(
  const struct Xt_stripe *src_indices, int num_src_stripes,
  const struct Xt_stripe *dst_indices, int num_dst_stripes,
  xmap_constructor new_xmap, MPI_Comm comm)
{
  Xt_idxlist
    src_idxlist = xt_idxstripes_new(src_indices, num_src_stripes),
    dst_idxlist = xt_idxstripes_new(dst_indices, num_dst_stripes);
  test_self_xmap_construct(src_idxlist, dst_idxlist, new_xmap, comm);
}

void
test_self_xmap_construct(Xt_idxlist src_idxlist,
                         Xt_idxlist dst_idxlist,
                         xmap_constructor new_xmap, MPI_Comm comm)
{
  // test of exchange map
  Xt_xmap xmap = new_xmap(src_idxlist, dst_idxlist, comm);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);
  assert_xmap_is_to_self(xmap);
  Xt_xmap xmap_copy = xt_xmap_copy(xmap);
  assert_xmap_is_to_self(xmap_copy);
  // clean up
  xt_xmap_delete(xmap);
  xt_xmap_delete(xmap_copy);
}

static inline void shift_idx(Xt_int idx[], int num, int offset)  {
  for (int i=0; i<num; i++) {
    idx[i] = (Xt_int)(idx[i] + my_rank * offset);
  }
}

static void
test_xmap1(xmap_constructor new_xmap,
           int lsize, MPI_Comm comm)
{
  struct Xt_stripe src_stripe, dst_stripe;
  // source index list
  src_stripe.nstrides = lsize;
  src_stripe.start = (Xt_int)(1 + (Xt_int)my_rank * (Xt_int)lsize);
  src_stripe.stride = 1;

  // destination index list
  dst_stripe.nstrides = lsize;
  dst_stripe.start = (Xt_int)(src_stripe.start + lsize - 1);
  dst_stripe.stride = -1;

  test_self_xmap_construct_idxstripes(&src_stripe, 1, &dst_stripe, 1,
                                      new_xmap, comm);
}

static void
test_xmap2(xmap_constructor new_xmap, MPI_Comm comm)
{
  // source index list
  Xt_int src_indices[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
  enum { num_src_indices = sizeof(src_indices) / sizeof(src_indices[0]) };
  shift_idx(src_indices, num_src_indices, 100);

  // destination index list
  Xt_int dst_indices[] = {5,4,3,96,1,5,4,5,4,3,13,2,1};
  enum { num_dst_indices = sizeof(dst_indices) / sizeof(dst_indices[0]) };
  shift_idx(dst_indices, num_dst_indices, 100);

  test_self_xmap_construct_idxvec(src_indices, num_src_indices,
                                  dst_indices, num_dst_indices,
                                  new_xmap, comm);
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
