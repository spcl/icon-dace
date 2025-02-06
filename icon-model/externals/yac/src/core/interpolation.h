// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

// YAC PUBLIC HEADER START

/** \example test_interpolation_parallel1_c.c
 * A test for internal interpolation routines.
 */

/** \example test_interpolation_parallel1.F90
 * A test for internal interpolation routines.
 */

/** \example test_interpolation_parallel2.c
 * A test for internal interpolation routines.
 */

/** \example test_interpolation_parallel3.c
 * A test for internal interpolation routines.
 */

/** \example test_interpolation_parallel4.c
 * A test for internal interpolation routines.
 */

/** \example test_interpolation_parallel5.c
 * A test for internal interpolation routines.
 */

/** \example test_interpolation_parallel6.c
 * A test for internal interpolation routines.
 */

#include <yaxt.h>

extern double const YAC_FRAC_MASK_NO_VALUE;
extern double const YAC_FRAC_MASK_UNDEF;

struct yac_interpolation;

struct yac_interpolation * yac_interpolation_new(
  size_t collection_size, double frac_mask_fallback_value,
  double scale_factor, double scale_summand);

void yac_interpolation_inc_ref_count(struct yac_interpolation * interpolation);

int yac_interpolation_with_frac_mask(struct yac_interpolation * interpolation);

void yac_interpolation_add_fixed(
  struct yac_interpolation * interp, double value, size_t count, size_t * pos);

void yac_interpolation_add_direct(
  struct yac_interpolation * interp, Xt_redist redist);

void yac_interpolation_add_direct_mf(
  struct yac_interpolation * interp, Xt_redist * redists, size_t num_src_fields);

/*
 * The target points consist of the sum of a number of source points. The sum
 * can be computed on the source or target processes.
 */
void yac_interpolation_add_sum_at_src(
  struct yac_interpolation * interp, Xt_redist * halo_redists,
  size_t tgt_count, size_t * num_src_per_tgt,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, Xt_redist result_redist);

void yac_interpolation_add_sum_at_tgt(
  struct yac_interpolation * interp, Xt_redist * src_redists,
  size_t * tgt_pos, size_t tgt_count, size_t * num_src_per_tgt,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields);

/*
 * The target points consist of a weighted sum of a number of source points. The
 * operation can be expressed as a distributed Matrix-Vector-Product. This
 * Product can be computed on the source or target processes.
 */
void yac_interpolation_add_weight_sum_mvp_at_src(
  struct yac_interpolation * interp, Xt_redist * halo_redists,
  size_t tgt_count, size_t * num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, Xt_redist result_redist);

void yac_interpolation_add_weight_sum_mvp_at_tgt(
  struct yac_interpolation * interp, Xt_redist * src_redists,
  size_t * tgt_pos, size_t tgt_count, size_t * num_src_per_tgt,
  double * weights, size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields);

struct yac_interpolation * yac_interpolation_copy(
  struct yac_interpolation * interp);

// src_fields dimensions [collection_idx]
//                       [field index]
//                       [local_idx]
// tgt_field dimensions [collection_idx]
//                      [local_idx]
void yac_interpolation_execute(
  struct yac_interpolation * interp, double *** src_fields,
  double ** tgt_field);
void yac_interpolation_execute_frac(
  struct yac_interpolation * interp, double *** src_fields,
  double *** src_frac_masks, double ** tgt_field);

// src_fields dimensions [collection_idx]
//                       [field index]
//                       [local_idx]
void yac_interpolation_execute_put(
  struct yac_interpolation * interp, double *** src_fields);
void yac_interpolation_execute_put_frac(
  struct yac_interpolation * interp, double *** src_fields,
  double *** src_frac_masks);

// tgt_field dimensions [collection_idx]
//                      [local_idx]
void yac_interpolation_execute_get(
  struct yac_interpolation * interp, double ** tgt_field);
void yac_interpolation_execute_get_async(
  struct yac_interpolation * interp, double ** tgt_field);

int yac_interpolation_execute_put_test(struct yac_interpolation * interp);
int yac_interpolation_execute_get_test(struct yac_interpolation * interp);

void yac_interpolation_execute_wait(struct yac_interpolation * interp);

void yac_interpolation_delete(struct yac_interpolation * interp);

// YAC PUBLIC HEADER STOP

#endif // INTERPOLATION_H
