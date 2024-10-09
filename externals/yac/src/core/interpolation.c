// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <string.h>
#include <math.h>

#include "ppm/core.h"
#include "utils_core.h"
#include "interpolation_internal.h"
#include "interpolation_fixed.h"
#include "interpolation_direct.h"
#include "interpolation_direct_mf.h"
#include "interpolation_sum_mvp_at_src.h"
#include "interpolation_sum_mvp_at_tgt.h"
#include "yac_mpi_internal.h"

double const YAC_FRAC_MASK_NO_VALUE = 133713371337.0;
double const YAC_FRAC_MASK_UNDEF = -133713371337.0;

struct yac_interpolation {

  struct yac_interpolation_type ** interps;
  size_t interp_count;

  int is_source, is_target;
  size_t collection_size;

  int ref_count;

  double frac_mask_fallback_value;

  double scale_factor;
  double scale_summand;
};

struct yac_interpolation * yac_interpolation_new(
  size_t collection_size, double frac_mask_fallback_value,
  double scale_factor, double scale_summand) {

  struct yac_interpolation * interp = xmalloc(1 * sizeof(*interp));

  interp->interps = NULL;
  interp->interp_count = 0;

  interp->is_source = 0;
  interp->is_target = 0;
  interp->collection_size = collection_size;

  interp->ref_count = 1;

  interp->frac_mask_fallback_value = frac_mask_fallback_value;

  YAC_ASSERT_F(
    isnormal(scale_factor),
    "ERROR(yac_interpolation_new): \"%lf\" is not a valid scaling factor",
    scale_factor);
  YAC_ASSERT_F(
    (scale_summand == 0.0) || isnormal(scale_summand),
    "ERROR(yac_interpolation_new): \"%lf\" is not a valid scaling summand",
    scale_summand);

  interp->scale_factor = scale_factor;
  interp->scale_summand = scale_summand;

  return interp;
}

static void yac_interpolation_add(
  struct yac_interpolation * interp,
  struct yac_interpolation_type * interp_type) {

  if (interp_type == NULL) return;

  interp->interps =
    xrealloc(
      interp->interps,
      (interp->interp_count + 1) * sizeof(*(interp->interps)));
  interp->interps[interp->interp_count] = interp_type;
  interp->interp_count++;

  interp->is_source |= interp_type->vtable->is_source(interp_type);
  interp->is_target |= interp_type->vtable->is_target(interp_type);
}

void yac_interpolation_add_fixed(
  struct yac_interpolation * interp, double value, size_t count,
  size_t * pos) {

  yac_interpolation_add(
    interp,
    yac_interpolation_fixed_new(interp->collection_size, value, count, pos));
}

void yac_interpolation_add_direct(
  struct yac_interpolation * interp, Xt_redist redist) {

  yac_interpolation_add(
    interp, yac_interpolation_direct_new(interp->collection_size, redist));
}

void yac_interpolation_add_direct_mf(
  struct yac_interpolation * interp, Xt_redist * redists,
  size_t num_src_fields) {

  yac_interpolation_add(
    interp,
    yac_interpolation_direct_mf_new(
      interp->collection_size, redists, num_src_fields));
}

void yac_interpolation_add_sum_at_src(
  struct yac_interpolation * interp, Xt_redist * halo_redists,
  size_t tgt_count, size_t * num_src_per_tgt,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, Xt_redist result_redist) {

  yac_interpolation_add(
    interp,
    yac_interpolation_sum_mvp_at_src_new(
      interp->collection_size, halo_redists, tgt_count, num_src_per_tgt, NULL,
      src_field_idx, src_idx, num_src_fields, result_redist,
      interp->frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE));
}

void yac_interpolation_add_sum_at_tgt(
  struct yac_interpolation * interp, Xt_redist * src_redists,
  size_t * tgt_pos, size_t tgt_count, size_t * num_src_per_tgt,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields) {

  yac_interpolation_add(
    interp,
    yac_interpolation_sum_mvp_at_tgt_new(
      interp->collection_size, src_redists, tgt_pos, tgt_count, num_src_per_tgt,
      NULL, src_field_idx, src_idx, num_src_fields,
      interp->frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE));
}

void yac_interpolation_add_weight_sum_mvp_at_src(
  struct yac_interpolation * interp, Xt_redist * halo_redists,
  size_t tgt_count, size_t * num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, Xt_redist result_redist) {

  yac_interpolation_add(
    interp,
    yac_interpolation_sum_mvp_at_src_new(
      interp->collection_size, halo_redists, tgt_count, num_src_per_tgt,
      weights, src_field_idx, src_idx, num_src_fields, result_redist,
      interp->frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE));
}

void yac_interpolation_add_weight_sum_mvp_at_tgt(
  struct yac_interpolation * interp, Xt_redist * src_redists, size_t * tgt_pos,
  size_t tgt_count, size_t * num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields) {

  yac_interpolation_add(
    interp,
    yac_interpolation_sum_mvp_at_tgt_new(
      interp->collection_size, src_redists, tgt_pos, tgt_count, num_src_per_tgt,
      weights, src_field_idx, src_idx, num_src_fields,
      interp->frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE));
}

struct yac_interpolation * yac_interpolation_copy(
  struct yac_interpolation * interp) {

  struct yac_interpolation * interp_copy =
    yac_interpolation_new(
      interp->collection_size, interp->frac_mask_fallback_value,
      interp->scale_factor, interp->scale_summand);

  interp_copy->interps =
    xmalloc(interp->interp_count * sizeof(*interp_copy->interps));
  interp_copy->interp_count = interp->interp_count;

  for (size_t i = 0; i < interp->interp_count; ++i)
    interp_copy->interps[i] =
      interp->interps[i]->vtable->copy(interp->interps[i]);

  interp_copy->is_source = interp->is_source;
  interp_copy->is_target = interp->is_target;

  return interp_copy;
}

// src_fields dimensions [collection_idx]
//                       [field index]
//                       [local_idx]
// src_frac_masks dimensions [collection_idx]
//                           [field index]
//                           [local_idx]
// tgt_field dimensions [collection_idx]
//                      [local_idx]
void yac_interpolation_execute_frac(
  struct yac_interpolation * interp, double *** src_fields,
  double *** src_frac_masks, double ** tgt_field) {

  YAC_ASSERT(
    (interp->frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE) ||
    (src_frac_masks == NULL),
    "ERROR(yac_interpolation_execute_frac): "
    "interpolation was not built for dynamic fractional masking, "
    "use yac_interpolation_execute instead");

  for (size_t i = 0; i < interp->interp_count; ++i)
    interp->interps[i]->vtable->execute(
      interp->interps[i], src_fields, src_frac_masks, tgt_field,
      interp->frac_mask_fallback_value, interp->scale_factor,
      interp->scale_summand);
}

// src_fields dimensions [collection_idx]
//                       [field index]
//                       [local_idx]
// tgt_field dimensions [collection_idx]
//                      [local_idx]
void yac_interpolation_execute(
  struct yac_interpolation * interp, double *** src_fields,
  double ** tgt_field) {

  YAC_ASSERT(
    interp->frac_mask_fallback_value == YAC_FRAC_MASK_NO_VALUE,
    "ERROR(yac_interpolation_execute): "
    "interpolation was built for dynamic fractional masking, "
    "use yac_interpolation_execute_frac instead");

  yac_interpolation_execute_frac(interp, src_fields, NULL, tgt_field);
}

// src_fields dimensions [collection_idx]
//                       [field index]
//                       [local_idx]
// src_frac_masks dimensions [collection_idx]
//                           [field index]
//                           [local_idx]
void yac_interpolation_execute_put_frac(
  struct yac_interpolation * interp, double *** src_fields,
  double *** src_frac_masks) {

  YAC_ASSERT(
    (interp->frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE) ||
    (src_frac_masks == NULL),
    "ERROR(yac_interpolation_execute_put_frac): "
    "interpolation was built for dynamic fractional masking, "
    "use yac_interpolation_execute_put instead");

  if (!interp->is_source) return;

  for (size_t i = 0; i < interp->interp_count; ++i)
    interp->interps[i]->vtable->execute_put(
      interp->interps[i], src_fields, src_frac_masks,
      interp->is_target, interp->frac_mask_fallback_value,
      interp->scale_factor, interp->scale_summand);
}

// src_fields dimensions [collection_idx]
//                       [field index]
//                       [local_idx]
void yac_interpolation_execute_put(
  struct yac_interpolation * interp, double *** src_fields) {

  YAC_ASSERT(
    interp->frac_mask_fallback_value == YAC_FRAC_MASK_NO_VALUE,
    "ERROR(yac_interpolation_execute_put): "
    "interpolation was built for dynamic fractional masking, "
    "use yac_interpolation_execute_put_frac instead");

  yac_interpolation_execute_put_frac(interp, src_fields, NULL);
}

// tgt_field dimensions [collection_idx]
//                      [local_idx]
void yac_interpolation_execute_get(
  struct yac_interpolation * interp, double ** tgt_field) {

  if (!interp->is_target) return;

  for (size_t i = 0; i < interp->interp_count; ++i)
    interp->interps[i]->vtable->execute_get(
      interp->interps[i], tgt_field, interp->frac_mask_fallback_value,
      interp->scale_factor, interp->scale_summand);
}

void yac_interpolation_execute_get_async(
  struct yac_interpolation * interp, double ** tgt_field) {

  if (!interp->is_target) return;

  for (size_t i = 0; i < interp->interp_count; ++i)
    interp->interps[i]->vtable->execute_get_async(
      interp->interps[i], tgt_field, interp->frac_mask_fallback_value,
      interp->scale_factor, interp->scale_summand);
}

int yac_interpolation_execute_put_test(struct yac_interpolation * interp) {

  int test = 0;

  for (size_t i = 0; (i < interp->interp_count) && !test; ++i)
    test |=
      interp->interps[i]->vtable->execute_put_test(
        interp->interps[i]);

  return test;
}

int yac_interpolation_execute_get_test(struct yac_interpolation * interp) {

  int test = 0;

  for (size_t i = 0; (i < interp->interp_count) && !test; ++i)
    test |=
      interp->interps[i]->vtable->execute_get_test(
        interp->interps[i]);

  return test;
}

void yac_interpolation_execute_wait(struct yac_interpolation * interp) {

  for (size_t i = 0; i < interp->interp_count; ++i)
    interp->interps[i]->vtable->execute_wait(interp->interps[i]);
}

void yac_interpolation_inc_ref_count(struct yac_interpolation * interpolation) {

  interpolation->ref_count++;
}

int yac_interpolation_with_frac_mask(struct yac_interpolation * interpolation) {

  return interpolation->frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE;
}

void yac_interpolation_delete(struct yac_interpolation * interp) {

  if (interp == NULL) return;

  if(--(interp->ref_count)) return;

  for (size_t i = 0; i < interp->interp_count; ++i)
    interp->interps[i]->vtable->delete(interp->interps[i]);
  free(interp->interps);
  free(interp);
}

double yac_interpolation_get_const_frac_mask_no_value_c2f() {
  return YAC_FRAC_MASK_NO_VALUE;
}

double yac_interpolation_get_const_frac_mask_undef_c2f() {
  return YAC_FRAC_MASK_UNDEF;
}
