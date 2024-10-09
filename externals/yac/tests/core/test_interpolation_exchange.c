// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

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
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

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

  struct yac_interpolation_exchange * exchanges[] = {exchange, exchange_copy};
  enum {NUM_EXCH = sizeof(exchanges) / sizeof(exchanges[0])};

  for (int optional_put_get = 0; optional_put_get < 2; ++optional_put_get) {

// #define NVHPC_BUG
// Starting version 24.3, NVHPC generates invalid code here
#ifdef NVHPC_BUG
    check_exchange(exchanges[0], is_source, is_target, optional_put_get);
    check_exchange(exchanges[1], is_source, is_target, optional_put_get);
#else
    for (size_t exch_idx = 0; exch_idx < NUM_EXCH; ++exch_idx)
      check_exchange(
        exchanges[exch_idx], is_source, is_target, optional_put_get);
#endif
  }

  for (size_t exch_idx = 0; exch_idx < NUM_EXCH; ++exch_idx)
    yac_interpolation_exchange_delete(exchanges[exch_idx], "main cleanup");


  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void check_exchange(
  struct yac_interpolation_exchange * exchange,
  int is_source, int is_target, int optional_put_get) {

  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  if (yac_interpolation_exchange_is_source(exchange) != is_source)
    PUT_ERR("ERROR in yac_interpolation_exchange_is_source");
  if (yac_interpolation_exchange_is_target(exchange) != is_target)
    PUT_ERR("ERROR in yac_interpolation_exchange_is_target");

  double send_data_[1] = {(double)comm_rank};
  double const * send_data = {is_source?&(send_data_[0]):NULL};
  double recv_data_[2] = {-1.0, -1.0};
  double * recv_data = {is_target?&(recv_data_[0]):NULL};
  enum YAC_INTERP_EXCH_STATUS status;

  // there should be not pending communication
  if (!yac_interpolation_exchange_put_test(exchange, "check_exchange"))
    PUT_ERR("ERROR in yac_interpolation_exchange_execute");
  if (!yac_interpolation_exchange_get_test(exchange, "check_exchange"))
    PUT_ERR("ERROR in yac_interpolation_exchange_execute");
  if (yac_interpolation_exchange_status(exchange, "check_exchange") !=
      YAC_INTERP_EXCH_IDLE)
    PUT_ERR("ERROR in yac_interpolation_exchange_status");

  // this should not block
  yac_interpolation_exchange_wait(exchange, "check_exchange");

  // simple exchange
  recv_data_[0] = -1.0, recv_data_[1] = -1.0;
  yac_interpolation_exchange_execute(
    exchange, &send_data, &recv_data, "check_exchange");
  if (is_target)
    if ((recv_data_[0] != 0.0) || (recv_data_[1] != 1.0))
      PUT_ERR("ERROR in yac_interpolation_exchange_execute");
  if (yac_interpolation_exchange_status(exchange, "check_exchange") !=
      YAC_INTERP_EXCH_IDLE)
    PUT_ERR("ERROR in yac_interpolation_exchange_status");

  // simple put/get exchange
  if (!optional_put_get || is_source)
    yac_interpolation_exchange_execute_put(
      exchange, &send_data, "check_exchange");
  status = yac_interpolation_exchange_status(exchange, "check_exchange");
  if (is_target) {
    if (is_source) {
      if (status != YAC_INTERP_EXCH_WAIT_PUT)
        PUT_ERR("ERROR in yac_interpolation_exchange_status");
      if (yac_interpolation_exchange_put_test(exchange, "check_exchange"))
        PUT_ERR("ERROR in yac_interpolation_exchange_execute");
    } else {
      if (status != YAC_INTERP_EXCH_IDLE)
        PUT_ERR("ERROR in yac_interpolation_exchange_status");
    }
  } else {
    if (is_source) {
      if ((status != YAC_INTERP_EXCH_ACTIVE) &&
          (status != YAC_INTERP_EXCH_IDLE))
        PUT_ERR("ERROR in yac_interpolation_exchange_status");
    } else {
      if (status != YAC_INTERP_EXCH_IDLE)
        PUT_ERR("ERROR in yac_interpolation_exchange_status");
    }
  }
  recv_data_[0] = -1.0, recv_data_[1] = -1.0;
  if (!optional_put_get || is_target)
    yac_interpolation_exchange_execute_get(
      exchange, &recv_data, "check_exchange");
  if (is_target && ((recv_data_[0] != 0.0) || (recv_data_[1] != 1.0)))
    PUT_ERR("ERROR in yac_interpolation_exchange_execute_get");
  if (is_source) yac_interpolation_exchange_wait(exchange, "check_exchange");
  if (yac_interpolation_exchange_status(exchange, "check_exchange") !=
      YAC_INTERP_EXCH_IDLE)
    PUT_ERR("ERROR in yac_interpolation_exchange_status");

  // simple put/get exchange with testing for completion
  if (!optional_put_get || is_source) {
    yac_interpolation_exchange_execute_put(
      exchange, &send_data, "check_exchange");
    if (!is_target)
      // this should return at some point
      while (!yac_interpolation_exchange_put_test(exchange, "check_exchange"));
  }
  recv_data_[0] = -1.0, recv_data_[1] = -1.0;
  if (!optional_put_get || is_target)
    yac_interpolation_exchange_execute_get(
      exchange, &recv_data, "check_exchange");
  if (is_target && ((recv_data_[0] != 0.0) || (recv_data_[1] != 1.0)))
    PUT_ERR("ERROR in yac_interpolation_exchange_execute_get");
  if (!optional_put_get || is_source)
    yac_interpolation_exchange_wait(exchange, "check_exchange");

  for (int put_get_order = 0; put_get_order < 2; ++put_get_order) {

    // async get
    for (int i = 0; i < 2; ++i) {
      if (put_get_order == i) {
        if (!optional_put_get || is_source)
          yac_interpolation_exchange_execute_put(
            exchange, &send_data, "check_exchange");
      }
      if (put_get_order != i) {
        recv_data_[0] = -1.0, recv_data_[1] = -1.0;
        if (!optional_put_get || is_target)
          yac_interpolation_exchange_execute_get_async(
            exchange, &recv_data, "check_exchange");
      }
    }
    yac_interpolation_exchange_wait(exchange, "check_exchange");
    if (is_target && ((recv_data_[0] != 0.0) || (recv_data_[1] != 1.0)))
      PUT_ERR("ERROR in yac_interpolation_exchange_execute_get");

    // async get with testing for completion
    for (int i = 0; i < 2; ++i) {
      if (put_get_order == i) {
        if (!optional_put_get || is_source)
          yac_interpolation_exchange_execute_put(
            exchange, &send_data, "check_exchange");
        if ((i == 0) && is_source && is_target)
          if (yac_interpolation_exchange_status(exchange, "check_exchange") !=
              YAC_INTERP_EXCH_WAIT_PUT)
            PUT_ERR("ERROR in yac_interpolation_exchange_status");
      }
      if (put_get_order != i) {
        recv_data_[0] = -1.0, recv_data_[1] = -1.0;
        if (!optional_put_get || is_target)
          yac_interpolation_exchange_execute_get_async(
            exchange, &recv_data, "check_exchange");
        if ((i == 0) && is_source && is_target)
          if (yac_interpolation_exchange_status(exchange, "check_exchange") !=
              YAC_INTERP_EXCH_WAIT_GET)
            PUT_ERR("ERROR in yac_interpolation_exchange_status");
      }
    }
    status = yac_interpolation_exchange_status(exchange, "check_exchange");
    if ((status != YAC_INTERP_EXCH_IDLE) &&
        (status != YAC_INTERP_EXCH_ACTIVE))
      PUT_ERR("ERROR in yac_interpolation_exchange_status");
    // this should return at some point
    while (!(yac_interpolation_exchange_put_test(exchange, "check_exchange") &&
             yac_interpolation_exchange_get_test(exchange, "check_exchange")));
    if (is_target && ((recv_data_[0] != 0.0) || (recv_data_[1] != 1.0)))
      PUT_ERR("ERROR in yac_interpolation_exchange_execute_get");
  }

  // introduce delay to better simulate asynchronous behaviour
  for (int i = 0; i < 4; ++i) {

    for (int j = 0; j < 2; ++j) {
      if ((i & 1) == j) {
        if (!optional_put_get || is_source)
          yac_interpolation_exchange_execute_put(
            exchange, &send_data, "check_exchange");
        if (is_source) sleep(1);
      }
      if ((i & 1) != j) {
        if (!optional_put_get || is_target)
          yac_interpolation_exchange_execute_get_async(
                exchange, &recv_data, "check_exchange");
      }
    }
  }
  yac_interpolation_exchange_wait(exchange, "check_exchange");

  // introduce delay to better simulate asynchronous behaviour
  for (int i = 0; i < 4; ++i) {
    if (is_source && is_target && (i & 1)) {
      yac_interpolation_exchange_execute(
        exchange, &send_data, &recv_data, "check_exchange");
    } else {
      yac_interpolation_exchange_execute_put(
        exchange, &send_data, "check_exchange");
      if (is_source) sleep(1);
      yac_interpolation_exchange_execute_get_async(
            exchange, &recv_data, "check_exchange");
    }
  }
  yac_interpolation_exchange_wait(exchange, "check_exchange");
}
