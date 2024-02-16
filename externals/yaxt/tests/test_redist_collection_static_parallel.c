/**
 * @file test_redist_collection_static_parallel.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <yaxt.h>

#include "core/ppm_xfuncs.h"
#include "tests.h"
#include "test_redist_common.h"

static void
test_4redist(MPI_Comm comm, Xt_config config);
static void
test_rr_exchange(MPI_Comm comm, Xt_config config);

int main(int argc, char **argv) {

  // init mpi

  int comm_size;

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);
  Xt_config config = redist_exchanger_option(&argc, &argv);

  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  if (comm_size > 1) {

    test_4redist(MPI_COMM_WORLD, config);
    test_rr_exchange(MPI_COMM_WORLD, config);
  }

  xt_config_delete(config);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
test_transpose_gather(Xt_redist redist,
                      Xt_int *dst, const Xt_int *src,
                      size_t size_a, size_t size_b, size_t size_all,
                      const Xt_int *index_vector_a,
                      const Xt_int *index_vector_b);

static void
test_4redist(MPI_Comm comm, Xt_config config)
{ // redist test with four different redists
  Xt_idxlist indices_a, indices_b, indices_all;
  int comm_size, comm_rank;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  assert(comm_size <= XT_INT_MAX / comm_size);
  int comm_size_sq = comm_size * comm_size;
  {
    Xt_idxlist indices_a_[2];
    Xt_int start = 0;
    Xt_int global_size[2] = {(Xt_int)(2*comm_size), (Xt_int)(comm_size_sq)};
    int local_size[2] = {comm_size,comm_size};
    Xt_int local_start[2][2]
      = {{0, (Xt_int)(comm_rank*comm_size)},
         {(Xt_int)comm_size, (Xt_int)(comm_size_sq-(comm_rank+1)*comm_size)}};

    for (size_t i = 0; i < 2; ++i)
      indices_a_[i] = xt_idxsection_new(start, 2, global_size, local_size,
                                        local_start[i]);

    indices_a = xt_idxlist_collection_new(indices_a_, 2);

    for (size_t i = 0; i < 2; ++i)
      xt_idxlist_delete(indices_a_[i]);
  }

  {
    assert(comm_size - 1 <= INT_MAX / 2 / comm_size_sq);
    struct Xt_stripe stripe = {.start = (Xt_int)(comm_rank*2*comm_size_sq),
                               .stride = 1, .nstrides = 2*comm_size_sq};

    indices_b = xt_idxstripes_new(&stripe, 1);
  }

  {
    assert(comm_size <= INT_MAX / 2 / comm_size_sq);
    struct Xt_stripe stripe = {.start = 0, .stride = 1,
                               .nstrides = 2*comm_size_sq*comm_size};

    indices_all = xt_idxstripes_new(&stripe, 1);
  }

  size_t size_a = 2*(size_t)comm_size*(size_t)comm_size,
    size_b = 2*(size_t)comm_size*(size_t)comm_size,
    size_all = size_a * (size_t)comm_size;
  Xt_int *src = xmalloc(sizeof (*src) * (size_a + size_b + size_all));
  Xt_int *index_vector_a = src, *index_vector_b = src + size_a;
  /* Xt_int *index_vector_all = src + size_a + size_b; */

  xt_idxlist_get_indices(indices_a, src);
  xt_idxlist_get_indices(indices_b, src + size_a);
  xt_idxlist_get_indices(indices_all, src + size_a + size_b);

  Xt_xmap xmaps[4] = {xt_xmap_all2all_new(indices_a, indices_b, comm),
                      xt_xmap_all2all_new(indices_b, indices_a, comm),
                      xt_xmap_all2all_new(indices_a, indices_all, comm),
                      xt_xmap_all2all_new(indices_b, indices_all, comm)};

  xt_idxlist_delete(indices_a);
  xt_idxlist_delete(indices_b);
  xt_idxlist_delete(indices_all);

  Xt_redist redists[4];
  for (size_t i = 0; i < 4; ++i) {
    redists[i] = xt_redist_p2p_new(xmaps[i], Xt_int_dt);
    xt_xmap_delete(xmaps[i]);
  }

  Xt_int *dst = xmalloc(sizeof (*dst) * (size_a + size_b + 2 * size_all));

  Xt_int *results_1 = dst, *results_2 = dst + size_b,
    *results_3 = dst + size_b + size_a,
    *results_4 = dst + size_a + size_b + size_all;

  MPI_Aint src_displacements[4]
    = {0, (MPI_Aint)(size_a * sizeof(Xt_int)),
       0, (MPI_Aint)(size_a * sizeof(Xt_int))};
  MPI_Aint dst_displacements[4]
    = {0, (MPI_Aint)((size_t)(results_2 - results_1) * sizeof(Xt_int)),
       (MPI_Aint)((size_t)(results_3 - results_1) * sizeof(Xt_int)),
       (MPI_Aint)((size_t)(results_4 - results_1) * sizeof(Xt_int))};

  Xt_redist redist = xt_redist_collection_static_custom_new(
    redists, 4, src_displacements, dst_displacements, comm, config);

  // test communicator of redist

  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist), comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  for (size_t i = 0; i < 4; ++i)
    xt_redist_delete(redists[i]);

  test_transpose_gather(redist, dst, src, size_a, size_b, size_all,
                        index_vector_a, index_vector_b);
  Xt_redist redist_copy = xt_redist_copy(redist);
  xt_redist_delete(redist);
  test_transpose_gather(redist_copy, dst, src, size_a, size_b, size_all,
                        index_vector_a, index_vector_b);
  // clean up
  free(src);
  free(dst);
  xt_redist_delete(redist_copy);
}

static void
test_transpose_gather(Xt_redist redist,
                      Xt_int *dst, const Xt_int *src,
                      size_t size_a, size_t size_b, size_t size_all,
                      const Xt_int *index_vector_a,
                      const Xt_int *index_vector_b)
{
  for (int sync_mode = 0; sync_mode < 2; ++sync_mode) {

    fill_array_xt_int(dst, NULL, size_b + size_a + 2 * size_all);

    exchange1_func_ptr exchange1_func
      = sync_mode == 0 ? xt_redist_s_exchange1 : wrap_a_exchange1;
    exchange1_func(redist, src, dst);

    // check results
    Xt_int *results_1 = dst, *results_2 = dst + size_b,
      *results_3 = dst + size_b + size_a,
      *results_4 = dst + size_a + size_b + size_all;
    bool mismatch = false;
    for (size_t i = 0; i < size_b; ++i)
      mismatch |= (results_1[i] != index_vector_b[i]);
    if (mismatch)
      PUT_ERR("error on xt_redist_s_exchange\n");

    mismatch = false;
    for (size_t i = 0; i < size_a; ++i)
      mismatch |= (results_2[i] != index_vector_a[i]);
    if (mismatch)
      PUT_ERR("error on xt_redist_s_exchange\n");

    mismatch = false;
    for (size_t i = 0; i < size_all; ++i)
      mismatch |= (results_3[i] != (int)i);
    if (mismatch)
      PUT_ERR("error on xt_redist_s_exchange\n");

    mismatch = false;
    for (size_t i = 0; i < size_all; ++i)
      mismatch |= (results_4[i] != (int)i);
    if (mismatch)
      PUT_ERR("error on xt_redist_s_exchange\n");
  }
}

static void
test_rr_exchange(MPI_Comm comm, Xt_config config)
{ // redist test with two redists that do a round robin exchange in
  // different directions

  int comm_size, comm_rank;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  enum { numExch = 2, };

  enum { partSize = 5, };
  Xt_int src_indices_[partSize], dst_indices_[numExch][partSize];

  for (Xt_int i = 0; i < partSize; ++i) {
    src_indices_[i] = (Xt_int)(comm_rank * partSize + i);
    dst_indices_[0][i] = (Xt_int)((src_indices_[i] + 1)
                                  % ((Xt_int)comm_size*partSize));
    Xt_int temp = (Xt_int)(src_indices_[i] - 1);
    dst_indices_[1][i] = (Xt_int)((temp < 0)?(comm_size * partSize - 1):temp);
  }

  Xt_xmap xmaps[numExch];
  Xt_idxlist src_indices = xt_idxvec_new(src_indices_, partSize);
  for (size_t i = 0; i < numExch; ++i) {
    Xt_idxlist dst_indices = xt_idxvec_new(dst_indices_[i], partSize);
    xmaps[i] = xt_xmap_all2all_new(src_indices, dst_indices, comm);
    xt_idxlist_delete(dst_indices);
  }
  xt_idxlist_delete(src_indices);

  Xt_redist redists[numExch];
  for (size_t i = 0; i < numExch; ++i) {
    redists[i] = xt_redist_p2p_new(xmaps[i], Xt_int_dt);
    xt_xmap_delete(xmaps[i]);
  }

  Xt_int results[2][partSize];

  MPI_Aint src_displacements[numExch] = {0, 0},
    ofs = (MPI_Aint)((size_t)(results[0]-results[1])*sizeof(Xt_int)),
    dst_displacements[numExch] = {0, ofs};

  Xt_redist redist = xt_redist_collection_static_custom_new(
    redists, numExch, src_displacements, dst_displacements, comm, config);
  xt_redist_delete(redists[0]);
  xt_redist_delete(redists[1]);

  // test communicator of redist
  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist), comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");


  for (int sync_mode = 0; sync_mode < 2; ++sync_mode) {

    fill_array_xt_int(results, NULL, 2 * partSize);

    exchange1_func_ptr exchange1_func
      = sync_mode == 0 ? xt_redist_s_exchange1 : wrap_a_exchange1;
    exchange1_func(redist, (void*)src_indices_, (void*)results[1]);

    // check results
    bool mismatch = false;
    for (int i = 0; i < partSize; ++i)
      mismatch |= (results[1][i] != dst_indices_[0][i]
                   || results[0][i] != dst_indices_[1][i]);
    if (mismatch)
      PUT_ERR("error on xt_redist_s_exchange\n");
  }

  // clean up
  xt_redist_delete(redist);
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
