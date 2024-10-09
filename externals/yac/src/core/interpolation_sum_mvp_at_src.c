// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "interpolation_sum_mvp_at_src.h"
#include "utils_core.h"
#include "yaxt.h"
#include "interpolation_utils.h"
#include "interpolation_exchange.h"

static int yac_interpolation_sum_mvp_at_src_is_source(
  struct yac_interpolation_type * interp);
static int yac_interpolation_sum_mvp_at_src_is_target(
  struct yac_interpolation_type * interp);
static void yac_interpolation_sum_mvp_at_src_execute(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_sum_mvp_at_src_execute_put(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks,
  int is_target, double frac_mask_fallback_value,
  double scale_factor, double scale_summand);
static void yac_interpolation_sum_mvp_at_src_execute_get(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_sum_mvp_at_src_execute_get_async(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static int yac_interpolation_sum_mvp_at_src_execute_put_test(
  struct yac_interpolation_type * interp);
static int yac_interpolation_sum_mvp_at_src_execute_get_test(
  struct yac_interpolation_type * interp);
static void yac_interpolation_sum_mvp_at_src_execute_wait(
  struct yac_interpolation_type * interp);
static struct yac_interpolation_type * yac_interpolation_sum_mvp_at_src_copy(
  struct yac_interpolation_type * interp);
static void yac_interpolation_sum_mvp_at_src_delete(
  struct yac_interpolation_type * interp);

static struct yac_interpolation_type_vtable const interpolation_sum_mvp_at_src_vtable = {
  .is_source         = yac_interpolation_sum_mvp_at_src_is_source,
  .is_target         = yac_interpolation_sum_mvp_at_src_is_target,
  .execute           = yac_interpolation_sum_mvp_at_src_execute,
  .execute_put       = yac_interpolation_sum_mvp_at_src_execute_put,
  .execute_get       = yac_interpolation_sum_mvp_at_src_execute_get,
  .execute_get_async = yac_interpolation_sum_mvp_at_src_execute_get_async,
  .execute_put_test  = yac_interpolation_sum_mvp_at_src_execute_put_test,
  .execute_get_test  = yac_interpolation_sum_mvp_at_src_execute_get_test,
  .execute_wait      = yac_interpolation_sum_mvp_at_src_execute_wait,
  .copy              = yac_interpolation_sum_mvp_at_src_copy,
  .delete            = yac_interpolation_sum_mvp_at_src_delete,
};

struct yac_interpolation_sum_mvp_at_src {

  struct yac_interpolation_type_vtable const * vtable;

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
  size_t * prefix_num_src_per_tgt;
  double * weights;
  size_t * src_field_idx;
  size_t * src_idx;
  size_t num_src_fields;
  double ** src_fields_buffer;

  int is_source;
  int is_target;

  int * ref_count;
};

static struct yac_interpolation_type * yac_interpolation_sum_mvp_at_src_new_(
  size_t collection_size,
  struct yac_interpolation_buffer halo_data,
  struct yac_interpolation_buffer result_data,
  struct yac_interpolation_exchange * src2halo,
  struct yac_interpolation_exchange * result2tgt,
  size_t tgt_count, size_t * prefix_num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx, size_t num_src_fields,
  int with_frac_mask, int * ref_count) {

  struct yac_interpolation_sum_mvp_at_src * mvp_at_src =
    xmalloc(1 * sizeof(*mvp_at_src));

  mvp_at_src->vtable = &interpolation_sum_mvp_at_src_vtable;
  mvp_at_src->collection_size = collection_size;
  mvp_at_src->with_frac_mask = with_frac_mask;
  mvp_at_src->halo_data = halo_data;
  mvp_at_src->result_data = result_data;
  mvp_at_src->src2halo = src2halo;
  mvp_at_src->result2tgt = result2tgt;
  mvp_at_src->tgt_count = tgt_count;
  mvp_at_src->prefix_num_src_per_tgt = prefix_num_src_per_tgt;
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

  return (struct yac_interpolation_type *)mvp_at_src;
}

struct yac_interpolation_type * yac_interpolation_sum_mvp_at_src_new(
  size_t collection_size, Xt_redist * halo_redists,
  size_t tgt_count, size_t * num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, Xt_redist result_redist_, int with_frac_mask) {

  size_t total_num_src = 0;
  size_t * prefix_num_src_per_tgt =
    xmalloc((tgt_count + 1) * sizeof(*prefix_num_src_per_tgt));
  for (size_t i = 0; i < tgt_count; ++i) {
    prefix_num_src_per_tgt[i] = total_num_src;
    total_num_src += num_src_per_tgt[i];
  }
  prefix_num_src_per_tgt[tgt_count] = total_num_src;

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
      tgt_count, prefix_num_src_per_tgt,
      (weights != NULL)?COPY_DATA(weights, total_num_src):NULL,
      COPY_DATA(src_field_idx, total_num_src),
      COPY_DATA(src_idx, total_num_src),
      num_src_fields, with_frac_mask, NULL);
}

static int yac_interpolation_sum_mvp_at_src_is_source(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src *)interp;

  return sum_mvp_at_src->is_source;
}

static int yac_interpolation_sum_mvp_at_src_is_target(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src *)interp;

  return sum_mvp_at_src->is_target;
}

static void yac_interpolation_sum_mvp_at_src_execute(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src *)interp;

  yac_interpolation_exchange_wait(
    sum_mvp_at_src->result2tgt,
    "yac_interpolation_sum_mvp_at_src_execute");

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
      results, NULL, sum_mvp_at_src->tgt_count,
      sum_mvp_at_src->prefix_num_src_per_tgt,
      sum_mvp_at_src->weights, sum_mvp_at_src->src_field_idx,
      sum_mvp_at_src->src_idx, num_src_fields, collection_size,
      frac_mask_fallback_value, scale_factor, scale_summand);
  }

  // redistribute results
  yac_interpolation_exchange_execute(
    sum_mvp_at_src->result2tgt, (double const **)results, tgt_field,
    "yac_interpolation_sum_mvp_at_src_execute");
}

static void yac_interpolation_sum_mvp_at_src_execute_put(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, int is_target,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(is_target);

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src *)interp;

  // wait until previous exchange is completed
  if (yac_interpolation_exchange_status(
        sum_mvp_at_src->result2tgt,
        "yac_interpolation_sum_mvp_at_src_execute_put") ==
      YAC_INTERP_EXCH_ACTIVE)
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
    results, NULL, sum_mvp_at_src->tgt_count,
    sum_mvp_at_src->prefix_num_src_per_tgt, sum_mvp_at_src->weights,
    sum_mvp_at_src->src_field_idx, sum_mvp_at_src->src_idx,
    num_src_fields, collection_size, frac_mask_fallback_value,
    scale_factor, scale_summand);

  yac_interpolation_exchange_execute_put(
    sum_mvp_at_src->result2tgt, (double const **)results,
    "yac_interpolation_sum_mvp_at_src_execute_put");
}

static void yac_interpolation_sum_mvp_at_src_execute_get(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(frac_mask_fallback_value);
  UNUSED(scale_factor);
  UNUSED(scale_summand);

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src*)interp;

  yac_interpolation_exchange_execute_get(
    sum_mvp_at_src->result2tgt, tgt_field,
    "yac_interpolation_sum_mvp_at_src_execute_get");
}

static void yac_interpolation_sum_mvp_at_src_execute_get_async(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(frac_mask_fallback_value);
  UNUSED(scale_factor);
  UNUSED(scale_summand);

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src*)interp;

  // wait until previous exchange is completed
  if (yac_interpolation_exchange_status(
        sum_mvp_at_src->result2tgt,
        "yac_interpolation_sum_mvp_at_src_execute_get_async") ==
      YAC_INTERP_EXCH_ACTIVE)
    yac_interpolation_exchange_wait(
      sum_mvp_at_src->result2tgt,
      "yac_interpolation_sum_mvp_at_src_execute_get_async");

  yac_interpolation_exchange_execute_get_async(
    sum_mvp_at_src->result2tgt, tgt_field,
    "yac_interpolation_sum_mvp_at_src_execute_get_async");
}

static int yac_interpolation_sum_mvp_at_src_execute_put_test(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src*)interp;

  return
    yac_interpolation_exchange_put_test(
      sum_mvp_at_src->result2tgt,
      "yac_interpolation_sum_mvp_at_src_execute_put_test");
}

static int yac_interpolation_sum_mvp_at_src_execute_get_test(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src*)interp;

  return
    yac_interpolation_exchange_get_test(
      sum_mvp_at_src->result2tgt,
      "yac_interpolation_sum_mvp_at_src_execute_get_test");
}

static void yac_interpolation_sum_mvp_at_src_execute_wait(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src*)interp;

  yac_interpolation_exchange_wait(
    sum_mvp_at_src->result2tgt,
    "yac_interpolation_sum_mvp_at_src_execute_wait");
}

static struct yac_interpolation_type *
  yac_interpolation_sum_mvp_at_src_copy(
    struct yac_interpolation_type * interp) {

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src*)interp;

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
      sum_mvp_at_src->prefix_num_src_per_tgt, sum_mvp_at_src->weights,
      sum_mvp_at_src->src_field_idx, sum_mvp_at_src->src_idx,
      sum_mvp_at_src->num_src_fields,
      sum_mvp_at_src->with_frac_mask, sum_mvp_at_src->ref_count);
}

static void yac_interpolation_sum_mvp_at_src_delete(
  struct yac_interpolation_type * interp) {

  if (interp == NULL) return;

  struct yac_interpolation_sum_mvp_at_src * sum_mvp_at_src =
    (struct yac_interpolation_sum_mvp_at_src*)interp;

  yac_interpolation_exchange_delete(
    sum_mvp_at_src->result2tgt, "yac_interpolation_sum_mvp_at_src_delete");
  yac_interpolation_buffer_free(&(sum_mvp_at_src->result_data));
  yac_interpolation_exchange_delete(
    sum_mvp_at_src->src2halo, "yac_interpolation_sum_mvp_at_src_delete");
  yac_interpolation_buffer_free(&(sum_mvp_at_src->halo_data));
  free(sum_mvp_at_src->src_fields_buffer);

  if (!--(*(sum_mvp_at_src->ref_count))) {
    free(sum_mvp_at_src->prefix_num_src_per_tgt);
    free(sum_mvp_at_src->src_idx);
    free(sum_mvp_at_src->src_field_idx);
    free(sum_mvp_at_src->weights);
    free(sum_mvp_at_src->ref_count);
  }

  free(sum_mvp_at_src);
}
