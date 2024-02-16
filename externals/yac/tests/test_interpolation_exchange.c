/**
 * @file test_interpolation_exchange.c
 *
 * @copyright Copyright  (C)  2023 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

#include "tests.h"
#include "interpolation_exchange.h"
#include "yac_mpi.h"
#include "test_common.h"
#include "yaxt.h"

static void check_exchange(
  struct yac_interpolation_exchange * exchange,
  int is_source, int is_target, int optional_put_get);

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  if (comm_size != 4) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  int is_source, is_target;
  switch (comm_rank) {
    case(0):
      is_source = 1;
      is_target = 1;
      break;
    case(1):
      is_source = 1;
      is_target = 0;
      break;
    case(2):
      is_source = 0;
      is_target = 1;
      break;
    default:
      is_source = 0;
      is_target = 0;
      break;
  }

  // generate redist
  Xt_idxlist src_idxlist =
    is_source?xt_idxvec_new((Xt_int[]){comm_rank}, 1):xt_idxempty_new();
  Xt_idxlist tgt_idxlist =
    is_target?xt_idxvec_new((Xt_int[]){0, 1}, 2):xt_idxempty_new();

  Xt_xmap xmap =
    xt_xmap_dist_dir_new(src_idxlist, tgt_idxlist, MPI_COMM_WORLD);

  Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

  // generate exchange
  size_t num_fields = 1;
  size_t collection_size = 1;
  int with_frac_mask = 0;
  struct yac_interpolation_exchange * exchange =
    yac_interpolation_exchange_new(
      &redist, num_fields, collection_size, with_frac_mask, "simple exchange");
  struct yac_interpolation_exchange * exchange_copy =
    yac_interpolation_exchange_copy(exchange);

  xt_redist_delete(redist);
  xt_xmap_delete(xmap);
  xt_idxlist_delete(tgt_idxlist);
  xt_idxlist_delete(src_idxlist);

  for (int optional_put_get = 0; optional_put_get < 2; ++optional_put_get) {
    check_exchange(exchange, is_source, is_target, optional_put_get);
    check_exchange(exchange_copy, is_source, is_target, optional_put_get);
  }

  yac_interpolation_exchange_delete(exchange, "main cleanup");
  yac_interpolation_exchange_delete(exchange_copy, "main cleanup");

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void check_exchange(
  struct yac_interpolation_exchange * exchange,
  int is_source, int is_target, int optional_put_get) {

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);

  if (yac_interpolation_exchange_is_source(exchange) != is_source)
    PUT_ERR("ERROR in yac_interpolation_exchange_is_source");
  if (yac_interpolation_exchange_is_target(exchange) != is_target)
    PUT_ERR("ERROR in yac_interpolation_exchange_is_target");

  double send_data_[1] = {(double)comm_rank};
  double const * send_data = {is_source?&(send_data_[0]):NULL};
  double recv_data_[2] = {-1.0, -1.0};
  double * recv_data = {is_target?&(recv_data_[0]):NULL};

  // there should be not pending communication
  if (!yac_interpolation_exchange_test(exchange, "check_exchange"))
    PUT_ERR("ERROR in yac_interpolation_exchange_execute");

  // this should not block
  yac_interpolation_exchange_wait(exchange, "check_exchange");

  // simple exchange
  recv_data_[0] = -1.0, recv_data_[1] = -1.0;
  yac_interpolation_exchange_execute(
    exchange, &send_data, &recv_data, "check_exchange");
  if (is_target)
    if ((recv_data_[0] != 0.0) || (recv_data_[1] != 1.0))
      PUT_ERR("ERROR in yac_interpolation_exchange_execute");

  // simple put/get exchange
  if (!optional_put_get || is_source)
    yac_interpolation_exchange_execute_put(
      exchange, &send_data, "check_exchange");
  recv_data_[0] = -1.0, recv_data_[1] = -1.0;
  if (!optional_put_get || is_target)
    yac_interpolation_exchange_execute_get(
      exchange, &recv_data, "check_exchange");
  if (is_target && ((recv_data_[0] != 0.0) || (recv_data_[1] != 1.0)))
    PUT_ERR("ERROR in yac_interpolation_exchange_execute_get");

  if (is_source) yac_interpolation_exchange_wait(exchange, "check_exchange");

  // simple put/get exchange with testing for completion
  if (!optional_put_get || is_source) {
    yac_interpolation_exchange_execute_put(
      exchange, &send_data, "check_exchange");
    if (!is_target)
      // this should return at some point
      while (!yac_interpolation_exchange_test(exchange, "check_exchange"));
  }
  recv_data_[0] = -1.0, recv_data_[1] = -1.0;
  if (!optional_put_get || is_target)
    yac_interpolation_exchange_execute_get(
      exchange, &recv_data, "check_exchange");
  if (is_target && ((recv_data_[0] != 0.0) || (recv_data_[1] != 1.0)))
    PUT_ERR("ERROR in yac_interpolation_exchange_execute_get");
  if (!optional_put_get || is_source)
    yac_interpolation_exchange_wait(exchange, "check_exchange");
}
