/**
 * @file interpolation_exchange.h
 *
 * @copyright Copyright  (C)  2019 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#ifndef INTERPOLATION_EXCHANGE_H
#define INTERPOLATION_EXCHANGE_H

#include "yaxt.h"

/** \example test_interpolation_exchange.c
 * This example show how to use interpolation exchanges.
 */

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
int yac_interpolation_exchange_test(
  struct yac_interpolation_exchange * exchange, char const * routine_name);
void yac_interpolation_exchange_execute_get(
  struct yac_interpolation_exchange * exchange, double ** recv_data,
  char const * routine_name);

void yac_interpolation_exchange_delete(
  struct yac_interpolation_exchange * exchange, char const * routine_name);

#endif // INTERPOLATION_EXCHANGE_H
