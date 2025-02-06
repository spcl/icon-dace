// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERPOLATION_EXCHANGE_H
#define INTERPOLATION_EXCHANGE_H

#include "yaxt.h"

/** \example test_interpolation_exchange.c
 * This example show how to use interpolation exchanges.
 */

enum YAC_INTERP_EXCH_STATUS {
  YAC_INTERP_EXCH_IDLE     = 0,
  YAC_INTERP_EXCH_WAIT_PUT = 1,
  YAC_INTERP_EXCH_WAIT_GET = 2,
  YAC_INTERP_EXCH_ACTIVE   = 3,
};

struct yac_interpolation_exchange;

struct yac_interpolation_exchange * yac_interpolation_exchange_new(
  Xt_redist * redists, size_t num_fields, size_t collection_size,
  int with_frac_mask, char const * name);
struct yac_interpolation_exchange * yac_interpolation_exchange_copy(
  struct yac_interpolation_exchange * exchange);

int yac_interpolation_exchange_is_source(
  struct yac_interpolation_exchange * exchange);
int yac_interpolation_exchange_is_target(
  struct yac_interpolation_exchange * exchange);

void yac_interpolation_exchange_execute(
  struct yac_interpolation_exchange * exchange,
  double const ** send_data, double ** recv_data, char const * routine_name);
void yac_interpolation_exchange_execute_put(
  struct yac_interpolation_exchange * exchange, double const ** send_data,
  char const * routine_name);
void yac_interpolation_exchange_wait(
  struct yac_interpolation_exchange * exchange, char const * routine_name);
int yac_interpolation_exchange_put_test(
  struct yac_interpolation_exchange * exchange, char const * routine_name);
int yac_interpolation_exchange_get_test(
  struct yac_interpolation_exchange * exchange, char const * routine_name);
enum YAC_INTERP_EXCH_STATUS yac_interpolation_exchange_status(
  struct yac_interpolation_exchange * exchange, char const * routine_name);
void yac_interpolation_exchange_execute_get(
  struct yac_interpolation_exchange * exchange, double ** recv_data,
  char const * routine_name);
void yac_interpolation_exchange_execute_get_async(
  struct yac_interpolation_exchange * exchange, double ** recv_data,
  char const * routine_name);

void yac_interpolation_exchange_delete(
  struct yac_interpolation_exchange * exchange, char const * routine_name);

#endif // INTERPOLATION_EXCHANGE_H
