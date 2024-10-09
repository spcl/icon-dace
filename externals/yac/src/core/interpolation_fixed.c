// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "interpolation_fixed.h"
#include "utils_core.h"
#include "interpolation_utils.h"

static int yac_interpolation_fixed_is_source(
  struct yac_interpolation_type * interp);
static int yac_interpolation_fixed_is_target(
  struct yac_interpolation_type * interp);
static void yac_interpolation_fixed_execute(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_fixed_execute_put(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks,
  int is_target, double frac_mask_fallback_value,
  double scale_factor, double scale_summand);
static void yac_interpolation_fixed_execute_get(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_fixed_execute_get_async(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static int yac_interpolation_fixed_execute_test(
  struct yac_interpolation_type * interp);
static void yac_interpolation_fixed_execute_wait(
  struct yac_interpolation_type * interp);
static struct yac_interpolation_type * yac_interpolation_fixed_copy(
  struct yac_interpolation_type * interp);
static void yac_interpolation_fixed_delete(
  struct yac_interpolation_type * interp);

static struct yac_interpolation_type_vtable const interpolation_fixed_vtable = {
  .is_source         = yac_interpolation_fixed_is_source,
  .is_target         = yac_interpolation_fixed_is_target,
  .execute           = yac_interpolation_fixed_execute,
  .execute_put       = yac_interpolation_fixed_execute_put,
  .execute_get       = yac_interpolation_fixed_execute_get,
  .execute_get_async = yac_interpolation_fixed_execute_get_async,
  .execute_put_test  = yac_interpolation_fixed_execute_test,
  .execute_get_test  = yac_interpolation_fixed_execute_test,
  .execute_wait      = yac_interpolation_fixed_execute_wait,
  .copy              = yac_interpolation_fixed_copy,
  .delete            = yac_interpolation_fixed_delete,
};

struct yac_interpolation_fixed {

  struct yac_interpolation_type_vtable const * vtable;

  size_t collection_size;

  double value;
  size_t * pos;
  size_t count;

  int ref_count;
};

struct yac_interpolation_type * yac_interpolation_fixed_new(
  size_t collection_size, double value, size_t count, size_t const * pos) {

  struct yac_interpolation_fixed * fixed = xmalloc(1 * sizeof(*fixed));

  fixed->vtable = &interpolation_fixed_vtable;
  fixed->collection_size = collection_size;
  fixed->value = value;
  fixed->count = count;
  fixed->pos = COPY_DATA(pos, count);
  fixed->ref_count = 1;

  return (struct yac_interpolation_type *)fixed;
}

static int yac_interpolation_fixed_is_source(
  struct yac_interpolation_type * interp) {

  UNUSED(interp);

  return 0;
}

static int yac_interpolation_fixed_is_target(
  struct yac_interpolation_type * interp) {

  return ((struct yac_interpolation_fixed*)interp)->count > 0;
}

static void yac_interpolation_fixed_execute_get_(
  struct yac_interpolation_type * interp, double ** tgt_field) {

  struct yac_interpolation_fixed * fixed =
    (struct yac_interpolation_fixed*)interp;

  double const value           = fixed->value;
  size_t * restrict pos    = fixed->pos;
  size_t const count           = fixed->count;
  size_t const collection_size = fixed->collection_size;

  for (size_t l = 0; l < collection_size; ++l)
    for (size_t j = 0; j < count; ++j)
      tgt_field[l][pos[j]] = value;
}

static void yac_interpolation_fixed_execute_get(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(frac_mask_fallback_value);
  UNUSED(scale_factor);
  UNUSED(scale_summand);

  yac_interpolation_fixed_execute_get_(interp, tgt_field);
}

static void yac_interpolation_fixed_execute_get_async(
  struct yac_interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  yac_interpolation_fixed_execute_get(
    interp, tgt_field, frac_mask_fallback_value, scale_factor, scale_summand);
}

static void yac_interpolation_fixed_execute(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(src_fields);
  UNUSED(src_frac_masks);
  UNUSED(frac_mask_fallback_value);
  UNUSED(scale_factor);
  UNUSED(scale_summand);

  yac_interpolation_fixed_execute_get_(interp, tgt_field);
}

static void yac_interpolation_fixed_execute_put(
  struct yac_interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks,
  int is_target, double frac_mask_fallback_value,
  double scale_factor, double scale_summand) {

  UNUSED(interp);
  UNUSED(src_fields);
  UNUSED(src_frac_masks);
  UNUSED(is_target);
  UNUSED(frac_mask_fallback_value);
  UNUSED(scale_factor);
  UNUSED(scale_summand);
  return;
}

static int yac_interpolation_fixed_execute_test(
  struct yac_interpolation_type * interp) {

  UNUSED(interp);
  return 1;
}

static void yac_interpolation_fixed_execute_wait(
  struct yac_interpolation_type * interp) {

  UNUSED(interp);
  return;
}

static struct yac_interpolation_type * yac_interpolation_fixed_copy(
  struct yac_interpolation_type * interp) {

  struct yac_interpolation_fixed * fixed =
    (struct yac_interpolation_fixed*)interp;

  fixed->ref_count++;

  return interp;
}

static void yac_interpolation_fixed_delete(
  struct yac_interpolation_type * interp) {

  if (interp == NULL) return;

  struct yac_interpolation_fixed * fixed =
    (struct yac_interpolation_fixed*)interp;

  if(--(fixed->ref_count)) return;

  free(fixed->pos);
  free(fixed);
}
