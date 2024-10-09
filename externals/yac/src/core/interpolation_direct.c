// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "interpolation_direct.h"
#include "utils_core.h"
#include "yaxt.h"
#include "interpolation_utils.h"
#include "interpolation_exchange.h"

static int yac_interpolation_direct_is_source(
  struct yac_interpolation_type * interp);
static int yac_interpolation_direct_is_target(
  struct yac_interpolation_type * interp);
static void yac_interpolation_direct_execute(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_direct_execute_put(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks,
  int is_target, double frac_mask_fallback_value,
  double scale_factor, double scale_summand);
static void yac_interpolation_direct_execute_get(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_direct_execute_get_async(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static int yac_interpolation_direct_execute_put_test(
  struct yac_interpolation_type * interp);
static int yac_interpolation_direct_execute_get_test(
  struct yac_interpolation_type * interp);
static void yac_interpolation_direct_execute_wait(
  struct yac_interpolation_type * interp);
static struct yac_interpolation_type * yac_interpolation_direct_copy(
  struct yac_interpolation_type * interp);
static void yac_interpolation_direct_delete(
  struct yac_interpolation_type * interp);

static struct yac_interpolation_type_vtable const interpolation_direct_vtable = {
  .is_source         = yac_interpolation_direct_is_source,
  .is_target         = yac_interpolation_direct_is_target,
  .execute           = yac_interpolation_direct_execute,
  .execute_put       = yac_interpolation_direct_execute_put,
  .execute_get       = yac_interpolation_direct_execute_get,
  .execute_get_async = yac_interpolation_direct_execute_get_async,
  .execute_put_test  = yac_interpolation_direct_execute_put_test,
  .execute_get_test  = yac_interpolation_direct_execute_get_test,
  .execute_wait      = yac_interpolation_direct_execute_wait,
  .copy              = yac_interpolation_direct_copy,
  .delete            = yac_interpolation_direct_delete,
};

struct yac_interpolation_direct {

  struct yac_interpolation_type_vtable const * vtable;

  size_t collection_size;

  struct yac_interpolation_buffer src_data;
  struct yac_interpolation_exchange * src2tgt;

  double ** src_field_buffer;
  int is_source;
  int is_target;
};

static struct yac_interpolation_type * yac_interpolation_direct_new_(
  size_t collection_size, struct yac_interpolation_exchange * src2tgt,
  struct yac_interpolation_buffer src_data) {

  struct yac_interpolation_direct * direct = xmalloc(1 * sizeof(*direct));

  direct->vtable = &interpolation_direct_vtable;
  direct->collection_size = collection_size;
  direct->src_data = src_data;
  direct->src2tgt = src2tgt;
  direct->src_field_buffer =
    xcalloc(collection_size, sizeof(*(direct->src_field_buffer)));
  direct->is_source = yac_interpolation_exchange_is_source(direct->src2tgt);
  direct->is_target = yac_interpolation_exchange_is_target(direct->src2tgt);

  return (struct yac_interpolation_type *)direct;
}

struct yac_interpolation_type * yac_interpolation_direct_new(
  size_t collection_size, Xt_redist redist_) {

  return
    yac_interpolation_direct_new_(
      collection_size,
      yac_interpolation_exchange_new(
        &redist_, 1, collection_size, 0, "source to target"),
      yac_interpolation_buffer_init(
        &redist_, 1, collection_size, SEND_BUFFER));
}

static int yac_interpolation_direct_is_source(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct *)interp;

  return direct->is_source;
}

static int yac_interpolation_direct_is_target(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct *)interp;

  return direct->is_target;
}

static void yac_interpolation_direct_execute(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct *)interp;

  double ** src_send_buffer = NULL;

  yac_interpolation_exchange_wait(
    direct->src2tgt, "yac_interpolation_direct_execute");

  if (direct->is_source) {

    if ((frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE) ||
        (scale_factor != 1.0) || (scale_summand != 0.0)) {

      src_send_buffer = direct->src_data.buffer;

      compute_tgt_field(
        (double const * restrict **)src_fields,
        (double const * restrict **)src_frac_masks, src_send_buffer,
        direct->src_data.buffer_sizes, 1,
        direct->collection_size, frac_mask_fallback_value,
        scale_factor, scale_summand);

    } else {

      src_send_buffer = direct->src_field_buffer;
      for (size_t i = 0; i < direct->collection_size; ++i)
        src_send_buffer[i] = src_fields[i][0];
    }
  }

  yac_interpolation_exchange_execute(
    direct->src2tgt, (double const **)src_send_buffer, tgt_field,
    "yac_interpolation_direct_execute");
}

static void yac_interpolation_direct_execute_put(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, int is_target,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(is_target);

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct *)interp;

  double ** src_send_buffer = NULL;

  if (direct->is_source) {

    // wait until previous exchange is completed
    if (yac_interpolation_exchange_status(
          direct->src2tgt, "yac_interpolation_direct_execute_put") ==
        YAC_INTERP_EXCH_ACTIVE)
      yac_interpolation_exchange_wait(
        direct->src2tgt, "yac_interpolation_direct_execute_put");

    src_send_buffer = direct->src_data.buffer;

    compute_tgt_field(
      (double const * restrict **)src_fields,
      (double const * restrict **)src_frac_masks, src_send_buffer,
      direct->src_data.buffer_sizes, 1,
      direct->collection_size, frac_mask_fallback_value,
      scale_factor, scale_summand);
  }

  yac_interpolation_exchange_execute_put(
    direct->src2tgt, (double const **)src_send_buffer,
    "yac_interpolation_direct_execute_put");
}

static void yac_interpolation_direct_execute_get(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(frac_mask_fallback_value);
  UNUSED(scale_factor);
  UNUSED(scale_summand);

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct*)interp;

  yac_interpolation_exchange_execute_get(
    direct->src2tgt, tgt_field, "yac_interpolation_direct_execute_get");
}

static void yac_interpolation_direct_execute_get_async(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(frac_mask_fallback_value);
  UNUSED(scale_factor);
  UNUSED(scale_summand);

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct*)interp;

  // wait until previous exchange is completed
  if (yac_interpolation_exchange_status(
        direct->src2tgt, "yac_interpolation_direct_execute_get_async") ==
      YAC_INTERP_EXCH_ACTIVE)
    yac_interpolation_exchange_wait(
      direct->src2tgt, "yac_interpolation_direct_execute_get_async");

  yac_interpolation_exchange_execute_get_async(
    direct->src2tgt, tgt_field, "yac_interpolation_direct_execute_get_async");
}

static struct yac_interpolation_type * yac_interpolation_direct_copy(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct*)interp;

  return
    yac_interpolation_direct_new_(
      direct->collection_size,
      yac_interpolation_exchange_copy(direct->src2tgt),
      yac_interpolation_buffer_copy(
        direct->src_data, 1, direct->collection_size));
}

static int yac_interpolation_direct_execute_put_test(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct*)interp;

  return
    yac_interpolation_exchange_put_test(
      direct->src2tgt, "yac_interpolation_direct_execute_put_test");
}

static int yac_interpolation_direct_execute_get_test(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct*)interp;

  return
    yac_interpolation_exchange_get_test(
      direct->src2tgt, "yac_interpolation_direct_execute_get_test");
}

static void yac_interpolation_direct_execute_wait(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct*)interp;

  yac_interpolation_exchange_wait(
    direct->src2tgt, "yac_interpolation_direct_execute_wait");
}

static void yac_interpolation_direct_delete(
  struct yac_interpolation_type * interp) {

  if (interp == NULL) return;

  struct yac_interpolation_direct * direct =
    (struct yac_interpolation_direct*)interp;

  yac_interpolation_exchange_delete(
    direct->src2tgt, "yac_interpolation_direct_delete");
  yac_interpolation_buffer_free(&(direct->src_data));
  free(direct->src_field_buffer);
  free(direct);
}
