/**
 * @file interpolation_fixed.c
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

#include "interpolation_fixed.h"
#include "utils.h"
#include "interpolation_utils.h"

static int yac_interpolation_fixed_is_source(
  struct interpolation_type * interp);
static int yac_interpolation_fixed_is_target(
  struct interpolation_type * interp);
static void yac_interpolation_fixed_execute(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_fixed_execute_put(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks,
  int is_target, double frac_mask_fallback_value,
  double scale_factor, double scale_summand);
static void yac_interpolation_fixed_execute_get(
  struct interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static int yac_interpolation_fixed_execute_test(
  struct interpolation_type * interp);
static struct interpolation_type * yac_interpolation_fixed_copy(
  struct interpolation_type * interp);
static void yac_interpolation_fixed_delete(
  struct interpolation_type * interp);

static struct interpolation_type_vtable const interpolation_fixed_vtable = {
  .is_source        = yac_interpolation_fixed_is_source,
  .is_target        = yac_interpolation_fixed_is_target,
  .execute          = yac_interpolation_fixed_execute,
  .execute_put      = yac_interpolation_fixed_execute_put,
  .execute_get      = yac_interpolation_fixed_execute_get,
  .execute_test     = yac_interpolation_fixed_execute_test,
  .copy             = yac_interpolation_fixed_copy,
  .delete           = yac_interpolation_fixed_delete,
};

struct interpolation_fixed {

  struct interpolation_type_vtable const * vtable;

  size_t collection_size;

  double value;
  size_t * pos;
  size_t count;

  int ref_count;
};

struct interpolation_type * yac_interpolation_fixed_new(
  size_t collection_size, double value, size_t count, size_t const * pos) {

  struct interpolation_fixed * fixed = xmalloc(1 * sizeof(*fixed));

  fixed->vtable = &interpolation_fixed_vtable;
  fixed->collection_size = collection_size;
  fixed->value = value;
  fixed->count = count;
  fixed->pos = COPY_DATA(pos, count);
  fixed->ref_count = 1;

  return (struct interpolation_type *)fixed;
}

static int yac_interpolation_fixed_is_source(
  struct interpolation_type * interp) {

  return 0;
}

static int yac_interpolation_fixed_is_target(
  struct interpolation_type * interp) {

  return ((struct interpolation_fixed*)interp)->count > 0;
}

static void yac_interpolation_fixed_execute_get_(
  struct interpolation_type * interp, double ** tgt_field) {

  struct interpolation_fixed * fixed = (struct interpolation_fixed*)interp;

  double const value           = fixed->value;
  size_t * restrict pos    = fixed->pos;
  size_t const count           = fixed->count;
  size_t const collection_size = fixed->collection_size;

  for (size_t l = 0; l < collection_size; ++l)
    for (size_t j = 0; j < count; ++j)
      tgt_field[l][pos[j]] = value;
}

static void yac_interpolation_fixed_execute_get(
  struct interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  UNUSED(frac_mask_fallback_value);
  UNUSED(scale_factor);
  UNUSED(scale_summand);

  yac_interpolation_fixed_execute_get_(interp, tgt_field);
}

static void yac_interpolation_fixed_execute(
  struct interpolation_type * interp,
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
  struct interpolation_type * interp,
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
  struct interpolation_type * interp) {

  UNUSED(interp);
  return 1;
}

static struct interpolation_type * yac_interpolation_fixed_copy(
  struct interpolation_type * interp) {

  struct interpolation_fixed * fixed = (struct interpolation_fixed*)interp;

  fixed->ref_count++;

  return interp;
}

static void yac_interpolation_fixed_delete(
  struct interpolation_type * interp) {

  if (interp == NULL) return;

  struct interpolation_fixed * fixed = (struct interpolation_fixed*)interp;

  if(--(fixed->ref_count)) return;

  free(fixed->pos);
  free(fixed);
}
