/**
 * @file test_redist_single_array_base_parallel.c
 *
 * @copyright Copyright  (C)  2017 Jörg Behrens <behrens@dkrz.de>
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
#include "core/ppm_xfuncs.h"

int main(int argc, char **argv) {

  // init mpi

  int rank, size;
  MPI_Comm comm = MPI_COMM_WORLD;

  xt_mpi_call(MPI_Init(NULL, NULL), comm);

  xt_initialize(comm);
  Xt_config config = redist_exchanger_option(&argc, &argv);

  xt_mpi_call(MPI_Comm_rank(comm, &rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &size), comm);

  // simple round robin
  {
    // send and receive messages
    struct Xt_redist_msg send_msgs[] = { {.rank = (rank + 1)%size,
                                          .datatype = MPI_DOUBLE} };
    struct Xt_redist_msg recv_msgs[] = { {.rank = (rank + size - 1)%size,
                                          .datatype = MPI_DOUBLE} };
    enum { nsend = sizeof(send_msgs) / sizeof(send_msgs[0]) };
    enum { nrecv = sizeof(recv_msgs) / sizeof(recv_msgs[0]) };

    // redist_single_array_base
    // test exchange
    double src_data[1];
    enum { num_dst_values = 1 };
    double ref_dst_data[num_dst_values];
    double dst_data[num_dst_values];
    src_data[0] = rank;
    ref_dst_data[0] =  (rank + size - 1)%size;
    test_redist_single_array_base(nsend, send_msgs, nrecv, recv_msgs,
                                  src_data, num_dst_values, dst_data,
                                  fill_array_double, NULL, ref_dst_data,
                                  MPI_DOUBLE, MPI_DOUBLE, comm, config);
  }

  // allgather
  {
    // send and receive messages
    int nsend = size;
    int nrecv = size;
    struct Xt_redist_msg *send_msgs
      = xmalloc((size_t)nsend * sizeof(*send_msgs));
    struct Xt_redist_msg *recv_msgs
      = xmalloc((size_t)nrecv * sizeof(*recv_msgs));

    for (int i = 0; i < nsend; ++i) {
      send_msgs[i].rank = i;
      send_msgs[i].datatype = MPI_DOUBLE;
    }
    for (int i = 0; i < nrecv; ++i) {
      recv_msgs[i].rank = i;
      xt_mpi_call(MPI_Type_create_indexed_block(
                    1, 1, &i, MPI_DOUBLE, &(recv_msgs[i].datatype)), comm);
      xt_mpi_call(MPI_Type_commit(&(recv_msgs[i].datatype)), comm);
    }

    // redist_single_array_base
    // test exchange
    double src_data[1];
    double *ref_dst_data = xmalloc((size_t)size * sizeof(*ref_dst_data));
    double *dst_data = xmalloc((size_t)size * sizeof(*dst_data));
    src_data[0] = rank;
    for (int i = 0; i < size; ++i) ref_dst_data[i] = i;
    test_redist_single_array_base(nsend, send_msgs, nrecv, recv_msgs,
                                  src_data, (size_t)size, dst_data,
                                  fill_array_double, NULL, ref_dst_data,
                                  MPI_DOUBLE, MPI_DOUBLE, comm, config);
    // clean up
    free(dst_data);
    free(ref_dst_data);
    for (int i = 0; i < nrecv; ++i)
      xt_mpi_call(MPI_Type_free(&(recv_msgs[i].datatype)), comm);
    free(recv_msgs);
    free(send_msgs);
  }

  // scatter by rank 0
  {
    // send and receive messages
    int nsend = (rank == 0)?(size):(0);
    int nrecv = 1;
    struct Xt_redist_msg *send_msgs
      = xmalloc((size_t)nsend * sizeof(*send_msgs));
    struct Xt_redist_msg recv_msgs[] = {{.rank = 0, .datatype = MPI_DOUBLE}};

    if (rank == 0) {
      for (int i = 0; i < nsend; ++i) {
        send_msgs[i].rank = i;
        xt_mpi_call(MPI_Type_create_indexed_block(
                      1, 1, &i, MPI_DOUBLE, &(send_msgs[i].datatype)), comm);
        xt_mpi_call(MPI_Type_commit(&(send_msgs[i].datatype)), comm);
      }
    }

    // test exchange
    double *src_data = xmalloc((size_t)nsend * sizeof(*src_data));
    enum { num_dst_values = 1 };
    double ref_dst_data[num_dst_values];
    double dst_data[num_dst_values];
    if (rank == 0) for (int i = 0; i < size; ++i) src_data[i] = i;
    ref_dst_data[0] = rank;

    test_redist_single_array_base(nsend, send_msgs, nrecv, recv_msgs,
                                  src_data, num_dst_values, dst_data,
                                  fill_array_double, NULL, ref_dst_data,
                                  MPI_DOUBLE, MPI_DOUBLE, comm, config);


    // clean up
    free(src_data);
    for (int i = 0; i < nsend; ++i)
      xt_mpi_call(MPI_Type_free(&(send_msgs[i].datatype)), comm);
    free(send_msgs);
  }

  xt_config_delete(config);
  xt_finalize();
  xt_mpi_call(MPI_Finalize(), comm);

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
