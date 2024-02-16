/**
 * @file interpolation_sum_mvp_at_src.c
 *
 * @copyright Copyright  (C)  2023 Moritz Hanke <hanke@dkrz.de>
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

#include <string.h>

#include "interpolation_sum_mvp_at_src.h"
#include "utils.h"
#include "yaxt.h"
#include "interpolation_utils.h"
#include "interpolation_exchange.h"

static int yac_interpolation_sum_mvp_at_src_is_source(
  struct interpolation_type * interp);
static int yac_interpolation_sum_mvp_at_src_is_target(
  struct interpolation_type * interp);
static void yac_interpolation_sum_mvp_at_src_execute(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_sum_mvp_at_src_execute_put(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks,
  int is_target, double frac_mask_fallback_value,
  double scale_factor, double scale_summand);
static void yac_interpolation_sum_mvp_at_src_execute_get(
  struct interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static int yac_interpolation_sum_mvp_at_src_execute_test(
  struct interpolation_type * interp);
static struct interpolation_type * yac_interpolation_sum_mvp_at_src_copy(
  struct interpolation_type * interp);
static void yac_interpolation_sum_mvp_at_src_delete(
  struct interpolation_type * interp);

static struct interpolation_type_vtable const interpolation_sum_mvp_at_src_vtable = {
  .is_source        = yac_interpolation_sum_mvp_at_src_is_source,
  .is_target        = yac_interpolation_sum_mvp_at_src_is_target,
  .execute          = yac_interpolation_sum_mvp_at_src_execute,
  .execute_put      = yac_interpolation_sum_mvp_at_src_execute_put,
  .execute_get      = yac_interpolation_sum_mvp_at_src_execute_get,
  .execute_test     = yac_interpolation_sum_mvp_at_src_execute_test,
  .copy             = yac_interpolation_sum_mvp_at_src_copy,
  .delete           = yac_interpolation_sum_mvp_at_src_delete,
};

struct interpolation_sum_mvp_at_src {

  struct interpolation_type_vtable const * vtable;

  size_t collection_size;
  int with_frac_mask;

  /* data flow:
   *  put:
   *      I. source processes collectively exchange data, such that each process
   *         can process its stencils
   *         - send buffer: src_fields (+src_frac_masks) provided by user
   *         - recv buffer: halo_data
   *         - exchange: src2halo
   *     II. all source processes process their stencils to compute the target
   *         field
   *         - from buffer: src_fields (+src_frac_masks) provided by user
   *                        halo_data
   *         - to buffer: result_data
   *    III. source processes send target field to target processes
   *         - send buffer: result_data
   *         - recv buffer: tgt_field provided by user
   *         - exchange: result2tgt
   *  get:
   *       I. target processes receive target field from source processes
   *         - send buffer: result_data
   *         - recv buffer: tgt_field provided by user
   *         - exchange: result2tgt
   */

  struct yac_interpolation_buffer halo_data;
  struct yac_interpolation_buffer result_data;
  struct yac_interpolation_exchange * src2halo;
  struct yac_interpolation_exchange * result2tgt;

  size_t tgt_count;
  size_t * num_src_per_tgt;
  double * weights;
  size_t * src_field_idx;
  size_t * src_idx;
  size_t num_src_fields;
  double ** src_fields_buffer;

  int is_source;
  int is_target;

  int * ref_count;
};

static struct interpolation_type * yac_interpolation_sum_mvp_at_src_new_(
  size_t collection_size,
  struct yac_interpolation_buffer halo_data,
  struct yac_interpolation_buffer result_data,
  struct yac_interpolation_exchange * src2halo,
  struct yac_interpolation_exchange * result2tgt,
  size_t tgt_count, size_t * num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx, size_t num_src_fields,
  int with_frac_mask, int * ref_count) {

  struct interpolation_sum_mvp_at_src * mvp_at_src =
    xmalloc(1 * sizeof(*mvp_at_src));

  mvp_at_src->vtable = &interpolation_sum_mvp_at_src_vtable;
  mvp_at_src->collection_size = collection_size;
  mvp_at_src->with_frac_mask = with_frac_mask;
  mvp_at_src->halo_data = halo_data;
  mvp_at_src->result_data = result_data;
  mvp_at_src->src2halo = src2halo;
  mvp_at_src->result2tgt = result2tgt;
  mvp_at_src->tgt_count = tgt_count;
  mvp_at_src->num_src_per_tgt = num_src_per_tgt;
  mvp_at_src->weights = weights;
  mvp_at_src->src_field_idx = src_field_idx;
  mvp_at_src->src_idx = src_idx;
  mvp_at_src->num_src_fields = num_src_fields;
  mvp_at_src->src_fields_buffer =
    xmalloc(
      (with_frac_mask?2*collection_size:collection_size) * num_src_fields *
      sizeof(*(mvp_at_src->src_fields_buffer)));
  mvp_at_src->is_source =
    yac_interpolation_exchange_is_source(mvp_at_src->src2halo) ||
    yac_interpolation_exchange_is_target(mvp_at_src->src2halo) ||
    yac_interpolation_exchange_is_source(mvp_at_src->result2tgt);
  mvp_at_src->is_target =
    yac_interpolation_exchange_is_target(mvp_at_src->result2tgt);

  mvp_at_src->ref_count =
    (ref_count == NULL)?xcalloc(1, sizeof(*mvp_at_src->ref_count)):ref_count;
  ++*(mvp_at_src->ref_count);

  return (struct interpolation_type *)mvp_at_src;
}

struct interpolation_type * yac_interpolation_sum_mvp_at_src_new(
  size_t collection_size, Xt_redist * halo_redists,
  size_t tgt_count, size_t * num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, Xt_redist result_redist_, int with_frac_mask) {

  size_t total_num_src = 0;
  for (size_t i = 0; i < tgt_count; ++i) total_num_src += num_src_per_tgt[i];

  return
    yac_interpolation_sum_mvp_at_src_new_(
      collection_size,
      yac_interpolation_buffer_init(
        halo_redists, num_src_fields,
        with_frac_mask?2*collection_size:collection_size, RECV_BUFFER),
      yac_interpolation_buffer_init(
        &result_redist_, 1, collection_size, SEND_BUFFER),
      yac_interpolation_exchange_new(
        halo_redists, num_src_fields,
        collection_size, with_frac_mask, "source to halo"),
      yac_interpolation_exchange_new(
        &result_redist_, 1, collection_size, 0, "result to target"),
      tgt_count,
      COPY_DATA(num_src_per_tgt, tgt_count),
      (weights != NULL)?COPY_DATA(weights, total_num_src):NULL,
      COPY_DATA(src_field_idx, total_num_src),
      COPY_DATA(src_idx, total_num_src),
      num_src_fields, with_frac_mask, NULL);
}

static int yac_interpolation_sum_mvp_at_src_is_source(
  struct interpolation_type * interp) {

  struct interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct interpolation_sum_mvp_at_src *)interp;

  return sum_mvp_at_src->is_source;
}

static int yac_interpolation_sum_mvp_at_src_is_target(
  struct interpolation_type * interp) {

  struct interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct interpolation_sum_mvp_at_src *)interp;

  return sum_mvp_at_src->is_target;
}

static void yac_interpolation_sum_mvp_at_src_execute(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct interpolation_sum_mvp_at_src *)interp;

  double ** results = sum_mvp_at_src->result_data.buffer;

  if (sum_mvp_at_src->is_source) {

    int with_frac_mask = sum_mvp_at_src->with_frac_mask;
    CHECK_WITH_FRAC_MASK("yac_interpolation_sum_mvp_at_src_execute")

    size_t num_src_fields = sum_mvp_at_src->num_src_fields;
    size_t collection_size = sum_mvp_at_src->collection_size;
    double ** temp_src_fields = sum_mvp_at_src->src_fields_buffer;
    double ** halo_buffers = sum_mvp_at_src->halo_data.buffer;
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < num_src_fields; ++j)
        temp_src_fields[i * num_src_fields + j] = src_fields[i][j];
    if (with_frac_mask)
      for (size_t i = 0; i < collection_size; ++i)
        for (size_t j = 0; j < num_src_fields; ++j)
          temp_src_fields[
            i * num_src_fields + j + collection_size * num_src_fields] =
              src_frac_masks[i][j];

    // do halo exchange
    yac_interpolation_exchange_execute(
      sum_mvp_at_src->src2halo, (double const **)temp_src_fields,
      halo_buffers, "yac_interpolation_sum_mvp_at_src_execute");

    compute_tgt_field_wgt(
      (double const * restrict **)src_fields,
      (double const * restrict **)(with_frac_mask?src_frac_masks:NULL),
      (double const * restrict *)halo_buffers,
      (double const * restrict *)(
        with_frac_mask?(halo_buffers + collection_size * num_src_fields):NULL),
      results, NULL, sum_mvp_at_src->tgt_count, sum_mvp_at_src->num_src_per_tgt,
      sum_mvp_at_src->weights, sum_mvp_at_src->src_field_idx, sum_mvp_at_src->src_idx,
      num_src_fields, collection_size, frac_mask_fallback_value,
      scale_factor, scale_summand);

    yac_interpolation_exchange_wait(
      sum_mvp_at_src->result2tgt,
      "yac_interpolation_sum_mvp_at_src_execute");
  }

  // redistribute results
  yac_interpolation_exchange_execute(
    sum_mvp_at_src->result2tgt, (double const **)results, tgt_field,
    "yac_interpolation_sum_mvp_at_src_execute");
}

static void yac_interpolation_sum_mvp_at_src_execute_put(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, int is_target,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct interpolation_sum_mvp_at_src *)interp;

  // wait until previous put is completed
  yac_interpolation_exchange_wait(
    sum_mvp_at_src->result2tgt,
    "yac_interpolation_sum_mvp_at_src_execute_put");

  int with_frac_mask = sum_mvp_at_src->with_frac_mask;
  CHECK_WITH_FRAC_MASK("yac_interpolation_sum_mvp_at_src_execute_put")

  size_t num_src_fields = sum_mvp_at_src->num_src_fields;
  size_t collection_size = sum_mvp_at_src->collection_size;
  double ** temp_src_fields = sum_mvp_at_src->src_fields_buffer;
  double ** halo_buffers = sum_mvp_at_src->halo_data.buffer;
  for (size_t i = 0; i < collection_size; ++i)
    for (size_t j = 0; j < num_src_fields; ++j)
      temp_src_fields[i * num_src_fields + j] = src_fields[i][j];
  if (with_frac_mask) {
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < num_src_fields; ++j)
        temp_src_fields[
          i * num_src_fields + j + collection_size * num_src_fields] =
            src_frac_masks[i][j];
  }

  // do halo exchange
  yac_interpolation_exchange_execute(
    sum_mvp_at_src->src2halo, (double const **)temp_src_fields, halo_buffers,
    "yac_interpolation_sum_mvp_at_src_execute_put");

  double ** results = sum_mvp_at_src->result_data.buffer;

  compute_tgt_field_wgt(
    (double const * restrict **)src_fields,
    (double const * restrict **)(with_frac_mask?src_frac_masks:NULL),
    (double const * restrict *)halo_buffers,
    (double const * restrict *)(
      with_frac_mask?(halo_buffers + collection_size * num_src_fields):NULL),
    results, NULL, sum_mvp_at_src->tgt_count, sum_mvp_at_src->num_src_per_tgt,
    sum_mvp_at_src->weights, sum_mvp_at_src->src_field_idx, sum_mvp_at_src->src_idx,
    num_src_fields, collection_size, frac_mask_fallback_value,
    scale_factor, scale_summand);

  yac_interpolation_exchange_execute_put(
    sum_mvp_at_src->result2tgt, (double const **)results,
    "yac_interpolation_sum_mvp_at_src_execute_put");
}

static void yac_interpolation_sum_mvp_at_src_execute_get(
  struct interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct interpolation_sum_mvp_at_src*)interp;

  yac_interpolation_exchange_execute_get(
    sum_mvp_at_src->result2tgt, tgt_field,
    "yac_interpolation_sum_mvp_at_src_execute_get");
}

static int yac_interpolation_sum_mvp_at_src_execute_test(
  struct interpolation_type * interp) {

  struct interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct interpolation_sum_mvp_at_src*)interp;

  return
    yac_interpolation_exchange_test(
      sum_mvp_at_src->result2tgt,
      "yac_interpolation_sum_mvp_at_src_execute_test");
}

static struct interpolation_type * yac_interpolation_sum_mvp_at_src_copy(
  struct interpolation_type * interp) {

  struct interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct interpolation_sum_mvp_at_src*)interp;

  return
    yac_interpolation_sum_mvp_at_src_new_(
      sum_mvp_at_src->collection_size,
      yac_interpolation_buffer_copy(
        sum_mvp_at_src->halo_data, sum_mvp_at_src->num_src_fields,
        sum_mvp_at_src->with_frac_mask?
          2*sum_mvp_at_src->collection_size:sum_mvp_at_src->collection_size),
      yac_interpolation_buffer_copy(
        sum_mvp_at_src->result_data, 1, sum_mvp_at_src->collection_size),
      yac_interpolation_exchange_copy(sum_mvp_at_src->src2halo),
      yac_interpolation_exchange_copy(sum_mvp_at_src->result2tgt),
      sum_mvp_at_src->tgt_count,
      sum_mvp_at_src->num_src_per_tgt, sum_mvp_at_src->weights,
      sum_mvp_at_src->src_field_idx, sum_mvp_at_src->src_idx,
      sum_mvp_at_src->num_src_fields,
      sum_mvp_at_src->with_frac_mask, sum_mvp_at_src->ref_count);
}

static void yac_interpolation_sum_mvp_at_src_delete(
  struct interpolation_type * interp) {

  if (interp == NULL) return;

  struct interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct interpolation_sum_mvp_at_src*)interp;

  yac_interpolation_exchange_delete(
    sum_mvp_at_src->result2tgt, "yac_interpolation_sum_mvp_at_src_delete");
  yac_interpolation_buffer_free(&(sum_mvp_at_src->result_data));
  yac_interpolation_exchange_delete(
    sum_mvp_at_src->src2halo, "yac_interpolation_sum_mvp_at_src_delete");
  yac_interpolation_buffer_free(&(sum_mvp_at_src->halo_data));
  free(sum_mvp_at_src->src_fields_buffer);

  if (!--(*(sum_mvp_at_src->ref_count))) {
    free(sum_mvp_at_src->num_src_per_tgt);
    free(sum_mvp_at_src->src_idx);
    free(sum_mvp_at_src->src_field_idx);
    free(sum_mvp_at_src->weights);
    free(sum_mvp_at_src->ref_count);
  }

  free(sum_mvp_at_src);
}
