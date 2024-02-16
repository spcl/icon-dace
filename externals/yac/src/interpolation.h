/**
 * @file interpolation.h
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

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "core/core.h"

/** \example test_interpolation_parallel1.c
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

extern double const YAC_FRAC_MASK_NO_VALUE;
extern double const YAC_FRAC_MASK_UNDEF;

struct interpolation_type;
struct interpolation_type_vtable {
  int (*is_source)(struct interpolation_type * interp);
  int (*is_target)(struct interpolation_type * interp);
  void (*execute)(struct interpolation_type * interp,
                  double *** src_fields, double *** src_frac_masks,
                  double ** tgt_field, double frac_mask_fallback_value,
                  double scale_factor, double scale_summand);
  void (*execute_put)(struct interpolation_type * interp,
                      double *** src_fields, double *** src_frac_masks,
                      int is_target, double frac_mask_fallback_value,
                      double scale_factor, double scale_summand);
  void (*execute_get)(struct interpolation_type * interp, double ** tgt_field,
                      double frac_mask_fallback_value, double scale_factor,
                      double scale_summand);
  int (*execute_test)(struct interpolation_type * interp);
  struct interpolation_type * (*copy)(struct interpolation_type * interp);
  void (*delete)(struct interpolation_type * interp);
};
struct interpolation_type {
  struct interpolation_type_vtable const * const vtable;
};

struct interpolation;

struct interpolation * yac_interpolation_new(
  size_t collection_size, double frac_mask_fallback_value,
  double scale_factor, double scale_summand);

void yac_interpolation_inc_ref_count(struct interpolation * interpolation);

int yac_interpolation_with_frac_mask(struct interpolation * interpolation);

void yac_interpolation_add_fixed(
  struct interpolation * interp, double value, size_t count, size_t * pos);

void yac_interpolation_add_direct(
  struct interpolation * interp, Xt_redist redist);

void yac_interpolation_add_direct_mf(
  struct interpolation * interp, Xt_redist * redists, size_t num_src_fields);

/*
 * The target points consist of the sum of a number of source points. The sum
 * can be computed on the source or target processes.
 */
void yac_interpolation_add_sum_at_src(
  struct interpolation * interp, Xt_redist * halo_redists,
  size_t tgt_count, size_t * num_src_per_tgt,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, Xt_redist result_redist);

void yac_interpolation_add_sum_at_tgt(
  struct interpolation * interp, Xt_redist * src_redists,
  size_t * tgt_pos, size_t tgt_count, size_t * num_src_per_tgt,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields);

/*
 * The target points consist of a weighted sum of a number of source points. The
 * operation can be expressed as a distributed Matrix-Vector-Product. This
 * Product can be computed on the source or target processes.
 */
void yac_interpolation_add_weight_sum_mvp_at_src(
  struct interpolation * interp, Xt_redist * halo_redists,
  size_t tgt_count, size_t * num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, Xt_redist result_redist);

void yac_interpolation_add_weight_sum_mvp_at_tgt(
  struct interpolation * interp, Xt_redist * src_redists,
  size_t * tgt_pos, size_t tgt_count, size_t * num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields);

struct interpolation * yac_interpolation_copy(struct interpolation * interp);

// src_fields dimensions [collection_idx]
//                       [field index]
//                       [local_idx]
// tgt_field dimensions [collection_idx]
//                      [local_idx]
void yac_interpolation_execute(
  struct interpolation * interp, double *** src_fields, double ** tgt_field);
void yac_interpolation_execute_frac(
  struct interpolation * interp, double *** src_fields,
  double *** src_frac_masks, double ** tgt_field);

// src_fields dimensions [collection_idx]
//                       [field index]
//                       [local_idx]
void yac_interpolation_execute_put(
  struct interpolation * interp, double *** src_fields);
void yac_interpolation_execute_put_frac(
  struct interpolation * interp, double *** src_fields,
  double *** src_frac_masks);

// tgt_field dimensions [collection_idx]
//                      [local_idx]
void yac_interpolation_execute_get(
  struct interpolation * interp, double ** tgt_field);

int yac_interpolation_execute_test(struct interpolation * interp);

void yac_interpolation_delete(struct interpolation * interp);

#endif // INTERPOLATION_H
