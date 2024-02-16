/**
 * @file test_redist_repeat_parallel.c
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
#include <stdlib.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "test_redist_common.h"

int main(int argc, char **argv) {

  // init mpi

  int rank, size;

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);
  Xt_config config = redist_exchanger_option(&argc, &argv);

  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  if (size > 1) {

    { // redist test with one redist and 4 levels

      Xt_idxlist indices_a, indices_b;

      {
        Xt_idxlist indices_a_[2];

        Xt_int start = 0;
        assert(size <= XT_INT_MAX / size);
        Xt_int global_size[2] = {(Xt_int)(2*size), (Xt_int)(size*size)};
        int local_size[2] = {size,size};
        Xt_int local_start[2][2]
          = {{ 0, (Xt_int)(rank*size)},
             { (Xt_int)size, (Xt_int)(size*size-(rank+1)*size) }};

        for (size_t i = 0; i < 2; ++i)
          indices_a_[i] = xt_idxsection_new(start, 2, global_size, local_size,
                                            local_start[i]);

        indices_a = xt_idxlist_collection_new(indices_a_, 2);

        xt_idxlist_delete(indices_a_[0]);
        xt_idxlist_delete(indices_a_[1]);
      }

      {
        assert(rank <= XT_INT_MAX / 2 / size / size);
        struct Xt_stripe stripe = {.start = (Xt_int)(rank*2*size*size),
                                   .stride = 1, .nstrides = 2*size*size};

        indices_b = xt_idxstripes_new(&stripe, 1);
      }

      Xt_int index_vector_a[2*size*size], index_vector_b[2*size*size];

      xt_idxlist_get_indices(indices_a, index_vector_a);
      xt_idxlist_get_indices(indices_b, index_vector_b);

      Xt_xmap xmap = xt_xmap_all2all_new(indices_a, indices_b, MPI_COMM_WORLD);

      xt_idxlist_delete(indices_a);
      xt_idxlist_delete(indices_b);

      Xt_redist redist_p2p = xt_redist_p2p_new(xmap, Xt_int_dt);

      xt_xmap_delete(xmap);

      enum { dim1a = 9, rpt_cnt = 4 };
      size_t dim0 = 2*(size_t)size*(size_t)size;
      Xt_int input[dim1a][dim0];
      Xt_int result[rpt_cnt][dim0];
      Xt_int result_2[dim1a][dim0];

      for (int j = 0; j < dim1a; ++j)
        for (size_t i = 0; i < dim0; ++i)
          input[j][i] = (Xt_int)(index_vector_a[i] + j * 2*(Xt_int)dim0);

      static const int displacements[2][rpt_cnt] = { {0,1,2,3}, {1,2,4,8} };
      MPI_Aint extent = (MPI_Aint)(dim0 * sizeof(Xt_int));

      Xt_redist redist_repeat[2];
      for (size_t i = 0; i < 2; ++i)
        redist_repeat[i] = xt_redist_repeat_custom_new(
          redist_p2p, extent, extent, rpt_cnt, displacements[i], config);

      // test communicator of redist_repeat

      for (size_t i = 0; i < 2; ++i)
        if (!communicators_are_congruent(xt_redist_get_MPI_Comm(
                                           redist_repeat[i]), MPI_COMM_WORLD))
          PUT_ERR("error in xt_redist_get_MPI_Comm\n");

      xt_redist_delete(redist_p2p);

      for (int sync_mode = 0; sync_mode < 2; ++sync_mode) {

        for (int j = 0; j < rpt_cnt; ++j)
          for (size_t i = 0; i < dim0; ++i)
            result[j][i] = -1;
        for (int j = 0; j < dim1a; ++j)
          for (size_t i = 0; i < dim0; ++i)
            result_2[j][i] = -1;

        if (sync_mode == 0) {
          xt_redist_s_exchange1(redist_repeat[0], input, result);
          xt_redist_s_exchange1(redist_repeat[1], input, result_2);
        } else {
          Xt_request request[2];
          xt_redist_a_exchange1(redist_repeat[0], input, result, request+0);
          xt_redist_a_exchange1(redist_repeat[1], input, result_2, request+1);
          xt_request_wait(request+0);
          xt_request_wait(request+1);
        }

        // check results

        bool mismatch = false;
        for (int j = 0; j < rpt_cnt; ++j)
          for (size_t i = 0; i < dim0; ++i)
            mismatch |= (result[j][i] != index_vector_b[i] + j * 2*(int)dim0);
        if (mismatch)
          PUT_ERR("ERROR: in first xt_redist_s_exchange1\n");

        mismatch = false;
        for (int j = 1; j <= 8; j<<=1)
          for (size_t i = 0; i < dim0; ++i)
            mismatch |= (result_2[j][i] != index_vector_b[i] + j * 2*(int)dim0);
        if (mismatch)
          PUT_ERR("ERROR: in second xt_redist_s_exchange1\n");
      }

      // clean up
      for (size_t i = 0; i < 2; ++i)
        xt_redist_delete(redist_repeat[i]);
    }
  }

  xt_config_delete(config);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
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
