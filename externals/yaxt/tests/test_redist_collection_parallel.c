/**
 * @file test_redist_collection_parallel.c
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

#include <mpi.h>
#include <yaxt.h>
#include "core/ppm_xfuncs.h"

#include "tests.h"
#include "ctest_common.h"
#include "test_idxlist_utils.h"
#include "test_redist_common.h"


enum {
  list_a = 0,
  list_b = 1,
  list_all = 2,
};

static void
test_4redist(MPI_Comm comm, Xt_config config);

static void
test_rr_redist(MPI_Comm comm, Xt_config config);

int main(int argc, char **argv) {

  test_init_mpi(&argc, &argv, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);
  Xt_config config = redist_exchanger_option(&argc, &argv);
  int comm_size;
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  if (comm_size > 1) {
    test_4redist(MPI_COMM_WORLD, config);

    test_rr_redist(MPI_COMM_WORLD, config);
  }

  xt_config_delete(config);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

enum { num_redists = 4 };

static void
exchange_4redist(Xt_redist redist, MPI_Comm comm,
                 const Xt_int *index_vector_a, const Xt_int *index_vector_b,
                 int sync);

static void
test_4redist(MPI_Comm comm, Xt_config config)
{
  int comm_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // redist test with four different redists
  Xt_idxlist indices_a, indices_b, indices_all;

  assert(comm_size <= XT_INT_MAX / comm_size
         && comm_size <= INT_MAX / comm_size);
  int comm_size_sq = comm_size * comm_size;

  {
    Xt_int start = 0;
    Xt_int global_size[2]
      = { (Xt_int)((Xt_int)2 * comm_size), (Xt_int)((Xt_int)comm_size_sq) };
    int local_size[2] = { comm_size, comm_size };
    Xt_int local_start[2][2] = {
      { 0, (Xt_int)((Xt_int)comm_rank * comm_size) },
      { (Xt_int)comm_size,
        (Xt_int)((Xt_int)comm_size_sq - (Xt_int)(comm_rank+1) * comm_size) }
    };

    Xt_idxlist indices_a_[2];
    for (size_t i = 0; i < 2; ++i)
      indices_a_[i] = xt_idxsection_new(start, 2, global_size, local_size,
                                        local_start[i]);

    indices_a = xt_idxlist_collection_new(indices_a_, 2);

    xt_idxlist_delete(indices_a_[0]);
    xt_idxlist_delete(indices_a_[1]);
  }

  {
    struct Xt_stripe stripe = {
      .start = (Xt_int)((Xt_int)comm_rank*2*comm_size_sq),
      .stride = 1, .nstrides = 2*comm_size_sq
    };
    indices_b = xt_idxstripes_new(&stripe, 1);
  }

  {
    assert(2 <= XT_INT_MAX / comm_size_sq / comm_size
           && 2 <= INT_MAX / comm_size_sq / comm_size);
    struct Xt_stripe stripe = { .start = 0, .stride = 1,
                                .nstrides = 2*comm_size_sq*comm_size };

    indices_all = xt_idxstripes_new(&stripe, 1);
  }

  const int list_sizes[3]
    = { 2*comm_size_sq, 2*comm_size_sq, 2 * comm_size_sq * comm_size };

  Xt_int *index_vector[2];
  for (size_t i = 0; i < 2; ++i)
    index_vector[i] = xmalloc((size_t)list_sizes[i] * sizeof (Xt_int));

  xt_idxlist_get_indices(indices_a, index_vector[list_a]);
  xt_idxlist_get_indices(indices_b, index_vector[list_b]);

  Xt_xmap xmaps[num_redists] = {
    xt_xmap_all2all_new(indices_a, indices_b, comm),
    xt_xmap_all2all_new(indices_b, indices_a, comm),
    xt_xmap_all2all_new(indices_a, indices_all, comm),
    xt_xmap_all2all_new(indices_b, indices_all, comm)
  };

  xt_idxlist_delete(indices_a);
  xt_idxlist_delete(indices_b);
  xt_idxlist_delete(indices_all);

  Xt_redist redists[num_redists] = {xt_redist_p2p_new(xmaps[0], Xt_int_dt),
                                    xt_redist_p2p_new(xmaps[1], Xt_int_dt),
                                    xt_redist_p2p_new(xmaps[2], Xt_int_dt),
                                    xt_redist_p2p_new(xmaps[3], Xt_int_dt)};

  for (size_t i = 0; i < num_redists; ++i)
    xt_xmap_delete(xmaps[i]);

  Xt_redist redist = xt_redist_collection_custom_new(redists, num_redists, -1,
                                                     comm, config);

  // test communicator of redist

  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist),
                                   comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");


  for (size_t i = 0; i < num_redists; ++i)
    xt_redist_delete(redists[i]);

  exchange_4redist(redist, comm, index_vector[list_a], index_vector[list_b], 0);
  exchange_4redist(redist, comm, index_vector[list_a], index_vector[list_b], 1);
  Xt_redist redist_copy = xt_redist_copy(redist);
  xt_redist_delete(redist);
  exchange_4redist(redist_copy, comm,
                   index_vector[list_a], index_vector[list_b], 0);
  exchange_4redist(redist_copy, comm,
                   index_vector[list_a], index_vector[list_b], 1);

  // clean up
  for (size_t i = 0; i < 2; ++i)
    free(index_vector[i]);
  xt_redist_delete(redist_copy);
}


static void
check_4redist_result(int size, void *results[4],
                     const Xt_int *index_vector_a,
                     const Xt_int *index_vector_b);

static void
exchange_4redist(Xt_redist redist, MPI_Comm comm,
                 const Xt_int *index_vector_a, const Xt_int *index_vector_b,
                 int sync)
{
  int comm_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  size_t comm_size_ = (size_t)comm_size;
  assert(comm_size_ <= SIZE_MAX / comm_size_ / comm_size_);
  size_t comm_size_sq_ = comm_size_ * comm_size_,
    comm_size_cb_ = comm_size_sq_ * comm_size_;
  const size_t result_sizes[num_redists] =
    { 2 * comm_size_sq_,
      2 * comm_size_sq_,
      2 * comm_size_cb_,
      2 * comm_size_cb_ };
  const size_t spacing[num_redists] = { 2, 14, 5, 8 };
  size_t buf_size = 0;
  for (size_t i = 0; i < num_redists; ++i)
    buf_size += spacing[i] + result_sizes[i];
  buf_size *= sizeof (Xt_int);
  void *buf = xmalloc(buf_size);
  unsigned char *temp = buf;
  void *results[num_redists];
  results[0] = (void *)(temp += spacing[0] * sizeof (Xt_int));
  for (size_t i = 1; i < num_redists; ++i)
    results[i] = (void *)(
      temp += (result_sizes[i-1] + spacing[i]) * sizeof (Xt_int));

  const void *input[num_redists]
    = { index_vector_a, index_vector_b, index_vector_a, index_vector_b };

  exchange_func_ptr exchange_func
    = sync ? xt_redist_s_exchange : wrap_a_exchange;
  exchange_func(redist, num_redists, input, results);
  check_4redist_result(comm_size, results, index_vector_a,
                       index_vector_b);
  /*
   * create another first buffer, to test successful
   * adaptation to different addresses...
   */
  if (comm_rank == 0)
    results[0] = buf;

  /* ...and repeat exchange */
  exchange_func(redist, num_redists, input, results);

  check_4redist_result(comm_size, results, index_vector_a, index_vector_b);
  free(buf);
}


static void
check_4redist_result(int comm_size, void *results[4],
                     const Xt_int *index_vector_a,
                     const Xt_int *index_vector_b)
{
  size_t comm_size_ = (size_t)comm_size,
    comm_size_sq_ = comm_size_ * comm_size_,
    comm_size_cb_ = comm_size_sq_ * comm_size_;
  if (cmp_idx_arrays(2 * comm_size_sq_, (Xt_int *)results[0], index_vector_b))
    PUT_ERR("error on xt_redist_s_exchange\n");

  if (cmp_idx_arrays(2 * comm_size_sq_, (Xt_int *)results[1], index_vector_a))
    PUT_ERR("error on xt_redist_s_exchange\n");

  bool mismatch = false;
  for (size_t i = 0; i < 2*comm_size_cb_; ++i)
    mismatch |= (((Xt_int *)results[2])[i] != (Xt_int)i);
  if (mismatch)
    PUT_ERR("error on xt_redist_s_exchange\n");

  mismatch = false;
  for (size_t i = 0; i < 2*comm_size_cb_; ++i)
    mismatch |= (((Xt_int *)results[3])[i] != (Xt_int)i);
  if (mismatch)
    PUT_ERR("error on xt_redist_s_exchange\n");
}

enum { elems_per_rank = 5, };

static void
test_rr_redist(MPI_Comm comm, Xt_config config)
{
  int comm_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // redist test with two redists that do a round robin exchange in
  // different directions
  Xt_int src_indices_[elems_per_rank], dst_indices_[2][elems_per_rank];

  for (Xt_int i = 0; i < elems_per_rank; ++i) {
    src_indices_[i] = (Xt_int)(comm_rank * elems_per_rank + i);
    dst_indices_[0][i]
      = (Xt_int)((src_indices_[i] + 1) % (comm_size * elems_per_rank));
    Xt_int temp = (Xt_int)(src_indices_[i] - 1);
    dst_indices_[1][i]
      = (Xt_int)((temp < 0)?(comm_size * elems_per_rank - 1):temp);
  }

  Xt_idxlist src_indices = xt_idxvec_new(src_indices_, elems_per_rank);
  Xt_redist redists[2];
  for (size_t i = 0; i < 2; ++i) {
    Xt_idxlist dst_indices = xt_idxvec_new(dst_indices_[i], elems_per_rank);
    Xt_xmap xmap = xt_xmap_all2all_new(src_indices, dst_indices, comm);
    xt_idxlist_delete(dst_indices);
    redists[i] = xt_redist_p2p_new(xmap, Xt_int_dt);
    xt_xmap_delete(xmap);
  }
  xt_idxlist_delete(src_indices);
  Xt_redist redist = xt_redist_collection_custom_new(redists, 2, -1, comm,
                                                     config);

  // test communicator of redist
  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist), comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  for (size_t i = 0; i < 2; ++i)
    xt_redist_delete(redists[i]);

  Xt_int results_[2][elems_per_rank];
  void *results[2] = {results_[0], results_[1]};
  const void *input[2] = {src_indices_, src_indices_};
  check_redist_coll(redist, sync_mode_test_all,
                    2, input, 2 * elems_per_rank, results,
                    results_, fill_array_xt_int, NULL,
                    dst_indices_, Xt_int_dt, Xt_int_dt);

  Xt_redist redist_copy = xt_redist_copy(redist);
  xt_redist_delete(redist);
  check_redist_coll(redist_copy, sync_mode_test_all,
                    2, input, 2 * elems_per_rank, results,
                    results_, fill_array_xt_int, NULL,
                    dst_indices_, Xt_int_dt, Xt_int_dt);

  // clean up
  xt_redist_delete(redist_copy);
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
