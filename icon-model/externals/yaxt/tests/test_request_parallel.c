/**
 * @file test_request_parallel.c
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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <yaxt.h>

#include "core/ppm_xfuncs.h"
#include "tests.h"
#include "ctest_common.h"
#include "test_redist_common.h"

int main(int argc, char **argv) {

  test_init_mpi(&argc, &argv, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  int rank, size;

  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  { // trivial test

    MPI_Request requests[] = {MPI_REQUEST_NULL};
    int n = sizeof(requests) / sizeof(requests[0]);
    MPI_Comm comm = MPI_COMM_WORLD;

    Xt_request request = xt_request_msgs_new(n, requests, comm);

    check_wait_request(&request);
  }

  { // each process receives and sends data from/to every other process

    size_t msg_size = 1 * 1024 * 1024; // 4 * 1 MiB
    int * send_buf = xmalloc(msg_size * sizeof(*send_buf));
    int * recv_buf = xmalloc((size_t)size * msg_size * sizeof(*recv_buf));

    for (size_t i = 0; i < msg_size; ++i) send_buf[i] = rank;
    for (size_t i = 0; i < (size_t)size * msg_size; ++i) recv_buf[i] = -1;

    int n = 2 * size;
    MPI_Request requests[n];
    MPI_Comm comm = MPI_COMM_WORLD;

    for (int i = 0; i < size; ++i)
      xt_mpi_call(MPI_Isend(send_buf, (int)msg_size, MPI_INT, i, 0,
                            MPI_COMM_WORLD, requests + i + size),
                  MPI_COMM_WORLD);
    for (int i = 0; i < size; ++i)
      xt_mpi_call(MPI_Irecv(recv_buf + (size_t)i * msg_size,
                            (int)msg_size, MPI_INT, i, 0,
                            MPI_COMM_WORLD, requests + i), MPI_COMM_WORLD);

    Xt_request request = xt_request_msgs_new(n, requests, comm);

    check_wait_request(&request);

    bool mismatch = false;
    for (int i = 0; i < size; ++i)
      for (size_t j = 0; j < msg_size; ++j)
        mismatch |= (recv_buf[(size_t)i * msg_size + j] != i);
    if (mismatch)
      PUT_ERR("recv_buf[i * msg_size + j] != i\n");
    free(send_buf);
    free(recv_buf);
  }

  { // each process receives and sends data from/to every other process with
    // request_msgs_packed

    int * send_buf = xmalloc(1 * sizeof(*send_buf));
    int * recv_buf = xmalloc(1 * sizeof(*recv_buf));
    void *send_bufp = send_buf, *recv_bufp = recv_buf;

    send_buf[0] = rank;
    recv_buf[0] = -1;

    MPI_Request requests[2];
    MPI_Comm comm = MPI_COMM_WORLD;

    xt_mpi_call(MPI_Isend(send_buf, 1, MPI_INT, (rank + 1)%size, 0,
                          MPI_COMM_WORLD, requests+0), MPI_COMM_WORLD);
    xt_mpi_call(MPI_Irecv(recv_buf, 1, MPI_INT, (rank - 1 + size)%size, 0,
                          MPI_COMM_WORLD, requests+1), MPI_COMM_WORLD);

    MPI_Datatype unpack_datatype[1] = {MPI_INT};
    int neigh_rank;
    Xt_request request =
      xt_request_msgs_packed_new(2, requests, comm, 1, 1, unpack_datatype,
                                 &recv_bufp, &send_bufp, &neigh_rank);

    check_wait_request(&request);

    if (neigh_rank != (rank - 1 + size)%size)
      PUT_ERR("neigh_rank != (rank - 1 + size)%size after xt_request_wait\n");
  }

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
