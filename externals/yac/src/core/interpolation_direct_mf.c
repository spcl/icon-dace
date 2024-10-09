// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "interpolation_direct_mf.h"
#include "utils_core.h"
#include "yaxt.h"
#include "interpolation_utils.h"
#include "interpolation_exchange.h"

static int yac_interpolation_direct_mf_is_source(
  struct yac_interpolation_type * interp);
static int yac_interpolation_direct_mf_is_target(
  struct yac_interpolation_type * interp);
static void yac_interpolation_direct_mf_execute(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_direct_mf_execute_put(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks,
  int is_target, double frac_mask_fallback_value,
  double scale_factor, double scale_summand);
static void yac_interpolation_direct_mf_execute_get(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_direct_mf_execute_get_async(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static int yac_interpolation_direct_mf_execute_put_test(
  struct yac_interpolation_type * interp);
static int yac_interpolation_direct_mf_execute_get_test(
  struct yac_interpolation_type * interp);
static void yac_interpolation_direct_mf_execute_wait(
  struct yac_interpolation_type * interp);
static struct yac_interpolation_type * yac_interpolation_direct_mf_copy(
  struct yac_interpolation_type * interp);
static void yac_interpolation_direct_mf_delete(
  struct yac_interpolation_type * interp);

static struct yac_interpolation_type_vtable const
  interpolation_direct_mf_vtable = {
  .is_source         = yac_interpolation_direct_mf_is_source,
  .is_target         = yac_interpolation_direct_mf_is_target,
  .execute           = yac_interpolation_direct_mf_execute,
  .execute_put       = yac_interpolation_direct_mf_execute_put,
  .execute_get       = yac_interpolation_direct_mf_execute_get,
  .execute_get_async = yac_interpolation_direct_mf_execute_get_async,
  .execute_put_test  = yac_interpolation_direct_mf_execute_put_test,
  .execute_get_test  = yac_interpolation_direct_mf_execute_get_test,
  .execute_wait      = yac_interpolation_direct_mf_execute_wait,
  .copy              = yac_interpolation_direct_mf_copy,
  .delete            = yac_interpolation_direct_mf_delete,
};

struct yac_interpolation_direct_mf {

  struct yac_interpolation_type_vtable const * vtable;

  size_t collection_size;

  struct yac_interpolation_buffer src_data;
  struct yac_interpolation_exchange * src2tgt;

  double ** src_field_buffer;
  double ** tgt_field_buffer;
  size_t num_src_fields;
  int is_source;
  int is_target;
};

static struct yac_interpolation_type * yac_interpolation_direct_mf_new_(
  size_t collection_size, struct yac_interpolation_exchange * src2tgt,
  struct yac_interpolation_buffer src_data, size_t num_src_fields) {

  struct yac_interpolation_direct_mf * direct_mf = xmalloc(1 * sizeof(*direct_mf));

  direct_mf->vtable = &interpolation_direct_mf_vtable;
  direct_mf->collection_size = collection_size;
  direct_mf->src_data = src_data;
  direct_mf->src2tgt = src2tgt;
  direct_mf->src_field_buffer =
    xcalloc(num_src_fields * collection_size,
            sizeof(*(direct_mf->src_field_buffer)));
  direct_mf->tgt_field_buffer =
    xcalloc(num_src_fields * collection_size,
            sizeof(*(direct_mf->tgt_field_buffer)));
  direct_mf->num_src_fields = num_src_fields;
  direct_mf->is_source = yac_interpolation_exchange_is_source(direct_mf->src2tgt);
  direct_mf->is_target = yac_interpolation_exchange_is_target(direct_mf->src2tgt);

  return (struct yac_interpolation_type *)direct_mf;
}

struct yac_interpolation_type * yac_interpolation_direct_mf_new(
  size_t collection_size, Xt_redist * redists, size_t num_src_fields) {

  return
    yac_interpolation_direct_mf_new_(
      collection_size,
      yac_interpolation_exchange_new(
        redists, num_src_fields, collection_size, 0, "source to target"),
      yac_interpolation_buffer_init(
        redists, num_src_fields, collection_size, SEND_BUFFER), num_src_fields);
}

static int yac_interpolation_direct_mf_is_source(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf *)interp;

  return direct_mf->is_source;
}

static int yac_interpolation_direct_mf_is_target(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf *)interp;

  return direct_mf->is_target;
}

static void yac_interpolation_direct_mf_execute(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf *)interp;

  double ** src_send_buffer = NULL;
  size_t collection_size = direct_mf->collection_size;
  size_t num_src_fields = direct_mf->num_src_fields;

  yac_interpolation_exchange_wait(
    direct_mf->src2tgt, "yac_interpolation_direct_mf_execute");

  if (direct_mf->is_source) {

    if ((frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE) ||
        (scale_factor != 1.0) || (scale_summand != 0.0)) {

      src_send_buffer = direct_mf->src_data.buffer;

      compute_tgt_field(
        (double const * restrict **)src_fields,
        (double const * restrict **)src_frac_masks, src_send_buffer,
        direct_mf->src_data.buffer_sizes, num_src_fields,
        collection_size, frac_mask_fallback_value,
        scale_factor, scale_summand);

    } else {

      src_send_buffer = direct_mf->src_field_buffer;
      for (size_t i = 0; i < collection_size; ++i)
        for (size_t j = 0; j < num_src_fields; ++j)
          src_send_buffer[i * num_src_fields + j] = src_fields[i][j];
    }
  }

  double ** tgt_field_buffer = direct_mf->tgt_field_buffer;
  if (direct_mf->is_target)
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < num_src_fields; ++j)
        tgt_field_buffer[i * num_src_fields + j] = tgt_field[i];

  // send source points to the target processes
  yac_interpolation_exchange_execute(
    direct_mf->src2tgt, (double const **)src_send_buffer, tgt_field_buffer,
    "yac_interpolation_direct_mf_execute");
}

static void yac_interpolation_direct_mf_execute_put(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, int is_target,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(is_target);

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf *)interp;

  double ** src_send_buffer = NULL;

  if (direct_mf->is_source) {

    // wait until previous exchange is completed
    if (yac_interpolation_exchange_status(
          direct_mf->src2tgt, "yac_interpolation_direct_mf_execute_put") ==
        YAC_INTERP_EXCH_ACTIVE)
      yac_interpolation_exchange_wait(
        direct_mf->src2tgt, "yac_interpolation_direct_mf_execute_put");

    src_send_buffer = direct_mf->src_data.buffer;

    compute_tgt_field(
      (double const * restrict **)src_fields,
      (double const * restrict **)src_frac_masks, src_send_buffer,
      direct_mf->src_data.buffer_sizes, direct_mf->num_src_fields,
      direct_mf->collection_size, frac_mask_fallback_value,
      scale_factor, scale_summand);
  }

  yac_interpolation_exchange_execute_put(
    direct_mf->src2tgt, (double const **)src_send_buffer,
    "yac_interpolation_direct_mf_execute_put");
}

static void yac_interpolation_direct_mf_execute_get(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(frac_mask_fallback_value);
  UNUSED(scale_factor);
  UNUSED(scale_summand);

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf*)interp;

  double ** tgt_field_buffer = NULL;

  if (direct_mf->is_target) {
    tgt_field_buffer = direct_mf->tgt_field_buffer;
    size_t collection_size = direct_mf->collection_size;
    size_t num_src_fields = direct_mf->num_src_fields;
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < num_src_fields; ++j)
        tgt_field_buffer[i * num_src_fields + j] = tgt_field[i];
  }

  yac_interpolation_exchange_execute_get(
    direct_mf->src2tgt, tgt_field_buffer,
    "yac_interpolation_direct_mf_execute_get");
}

static void yac_interpolation_direct_mf_execute_get_async(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(frac_mask_fallback_value);
  UNUSED(scale_factor);
  UNUSED(scale_summand);

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf*)interp;

  // wait until previous exchange is completed
  if (yac_interpolation_exchange_status(
        direct_mf->src2tgt, "yac_interpolation_direct_mf_execute_get_async") ==
      YAC_INTERP_EXCH_ACTIVE)
    yac_interpolation_exchange_wait(
      direct_mf->src2tgt, "yac_interpolation_direct_mf_execute_get_async");

  double ** tgt_field_buffer = NULL;

  if (direct_mf->is_target) {
    tgt_field_buffer = direct_mf->tgt_field_buffer;
    size_t collection_size = direct_mf->collection_size;
    size_t num_src_fields = direct_mf->num_src_fields;
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < num_src_fields; ++j)
        tgt_field_buffer[i * num_src_fields + j] = tgt_field[i];
  }

  yac_interpolation_exchange_execute_get_async(
    direct_mf->src2tgt, tgt_field_buffer,
    "yac_interpolation_direct_mf_execute_get_async");
}

static int yac_interpolation_direct_mf_execute_put_test(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf*)interp;

  return
    yac_interpolation_exchange_put_test(
      direct_mf->src2tgt, "yac_interpolation_direct_mf_execute_put_test");
}

static int yac_interpolation_direct_mf_execute_get_test(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf*)interp;

  return
    yac_interpolation_exchange_get_test(
      direct_mf->src2tgt, "yac_interpolation_direct_mf_execute_get_test");
}

static void yac_interpolation_direct_mf_execute_wait(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf*)interp;

  yac_interpolation_exchange_wait(
    direct_mf->src2tgt, "yac_interpolation_direct_mf_execute_wait");
}

static struct yac_interpolation_type * yac_interpolation_direct_mf_copy(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf*)interp;

  return
    yac_interpolation_direct_mf_new_(
      direct_mf->collection_size,
      yac_interpolation_exchange_copy(direct_mf->src2tgt),
      yac_interpolation_buffer_copy(
        direct_mf->src_data, direct_mf->num_src_fields,
        direct_mf->collection_size),
      direct_mf->num_src_fields);
}

static void yac_interpolation_direct_mf_delete(
  struct yac_interpolation_type * interp) {

  if (interp == NULL) return;

  struct yac_interpolation_direct_mf * direct_mf =
    (struct yac_interpolation_direct_mf*)interp;

  yac_interpolation_exchange_delete(
    direct_mf->src2tgt, "yac_interpolation_direct_mf_delete");
  yac_interpolation_buffer_free(&(direct_mf->src_data));
  free(direct_mf->src_field_buffer);
  free(direct_mf->tgt_field_buffer);
  free(direct_mf);
}
