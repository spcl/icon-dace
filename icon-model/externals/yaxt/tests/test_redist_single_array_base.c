/**
 * @file test_redist_single_array_base.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <yaxt.h>

#include "core/ppm_xfuncs.h"
#include "tests.h"
#include "ctest_common.h"
#include "test_redist_common.h"

int main(int argc, char **argv) {

  MPI_Comm comm = MPI_COMM_WORLD;

  test_init_mpi(&argc, &argv, comm);

  xt_initialize(comm);
  Xt_config config = redist_exchanger_option(&argc, &argv);

  // single double
  {
    // send and receive messages
    static const struct Xt_redist_msg send_msgs[] = { {.rank = 0, .datatype = MPI_DOUBLE} },
      recv_msgs[] = { {.rank = 0, .datatype = MPI_DOUBLE} };
    enum { nsend = sizeof(send_msgs) / sizeof(send_msgs[0]) };
    enum { nrecv = sizeof(recv_msgs) / sizeof(recv_msgs[0]) };

    // test exchange
    static const double src_data[] = {-5};
    static const double ref_dst_data[] = {-5};
    enum { num_ref_values = sizeof(ref_dst_data) / sizeof(ref_dst_data[0]) };
    double dst_data[num_ref_values];
    test_redist_single_array_base(nsend, send_msgs, nrecv, recv_msgs,
                                  src_data, num_ref_values, dst_data,
                                  fill_array_double, NULL, ref_dst_data,
                                  MPI_DOUBLE, MPI_DOUBLE, comm, config);
  }

  // reverse order of some doubles
  {
    // generate datatypes
    MPI_Datatype send_type, recv_type;
    static const int recv_displs[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
    enum { nelem = sizeof(recv_displs) / sizeof(recv_displs[0]) };
    xt_mpi_call(MPI_Type_contiguous(nelem, MPI_FLOAT, &send_type), comm);
    xt_mpi_call(MPI_Type_create_indexed_block(
                  nelem, 1, (int *)recv_displs, MPI_FLOAT, &recv_type), comm);
    xt_mpi_call(MPI_Type_commit(&send_type), comm);
    xt_mpi_call(MPI_Type_commit(&recv_type), comm);

    // send and receive messages
    struct Xt_redist_msg send_msgs[] = { {.rank = 0, .datatype = send_type} };
    struct Xt_redist_msg recv_msgs[] = { {.rank = 0, .datatype = recv_type} };
    enum { nsend = sizeof(send_msgs) / sizeof(send_msgs[0]) };
    enum { nrecv = sizeof(recv_msgs) / sizeof(recv_msgs[0]) };

    // test exchange
    static const float src_data[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    static const float ref_dst_data[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
    enum { num_ref_values = sizeof(ref_dst_data) / sizeof(ref_dst_data[0]) };
    float dst_data[num_ref_values];
    test_redist_single_array_base(nsend, send_msgs, nrecv, recv_msgs,
                                  src_data, num_ref_values, dst_data,
                                  fill_array_float, NULL, ref_dst_data,
                                  MPI_FLOAT, MPI_FLOAT, comm, config);
    // clean up
    xt_mpi_call(MPI_Type_free(&recv_type), comm);
    xt_mpi_call(MPI_Type_free(&send_type), comm);
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
