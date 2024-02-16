/**
 * @file interpolation_direct.c
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

#include "interpolation_direct.h"
#include "utils.h"
#include "yaxt.h"
#include "interpolation_utils.h"
#include "interpolation_exchange.h"

static int yac_interpolation_direct_is_source(
  struct interpolation_type * interp);
static int yac_interpolation_direct_is_target(
  struct interpolation_type * interp);
static void yac_interpolation_direct_execute(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_direct_execute_put(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks,
  int is_target, double frac_mask_fallback_value,
  double scale_factor, double scale_summand);
static void yac_interpolation_direct_execute_get(
  struct interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static int yac_interpolation_direct_execute_test(
  struct interpolation_type * interp);
static struct interpolation_type * yac_interpolation_direct_copy(
  struct interpolation_type * interp);
static void yac_interpolation_direct_delete(
  struct interpolation_type * interp);

static struct interpolation_type_vtable const interpolation_direct_vtable = {
  .is_source        = yac_interpolation_direct_is_source,
  .is_target        = yac_interpolation_direct_is_target,
  .execute          = yac_interpolation_direct_execute,
  .execute_put      = yac_interpolation_direct_execute_put,
  .execute_get      = yac_interpolation_direct_execute_get,
  .execute_test     = yac_interpolation_direct_execute_test,
  .copy             = yac_interpolation_direct_copy,
  .delete           = yac_interpolation_direct_delete,
};

struct interpolation_direct {

  struct interpolation_type_vtable const * vtable;

  size_t collection_size;

  struct yac_interpolation_buffer src_data;
  struct yac_interpolation_exchange * src2tgt;

  double ** src_field_buffer;
  int is_source;
  int is_target;
};

static struct interpolation_type * yac_interpolation_direct_new_(
  size_t collection_size, struct yac_interpolation_exchange * src2tgt,
  struct yac_interpolation_buffer src_data) {

  struct interpolation_direct * direct = xmalloc(1 * sizeof(*direct));

  direct->vtable = &interpolation_direct_vtable;
  direct->collection_size = collection_size;
  direct->src_data = src_data;
  direct->src2tgt = src2tgt;
  direct->src_field_buffer =
    xcalloc(collection_size, sizeof(*(direct->src_field_buffer)));
  direct->is_source = yac_interpolation_exchange_is_source(direct->src2tgt);
  direct->is_target = yac_interpolation_exchange_is_target(direct->src2tgt);

  return (struct interpolation_type *)direct;
}

struct interpolation_type * yac_interpolation_direct_new(
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
  struct interpolation_type * interp) {

  struct interpolation_direct * direct = (struct interpolation_direct *)interp;

  return direct->is_source;
}

static int yac_interpolation_direct_is_target(
  struct interpolation_type * interp) {

  struct interpolation_direct * direct = (struct interpolation_direct *)interp;

  return direct->is_target;
}

static void yac_interpolation_direct_execute(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct interpolation_direct * direct =
    (struct interpolation_direct *)interp;

  double ** src_send_buffer = NULL;

  if (direct->is_source) {

    yac_interpolation_exchange_wait(
      direct->src2tgt, "yac_interpolation_direct_execute");

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
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, int is_target,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct interpolation_direct * direct =
    (struct interpolation_direct *)interp;

  double ** src_send_buffer = NULL;

  if (direct->is_source) {

    // wait until previous put is completed
    yac_interpolation_exchange_wait(
      direct->src2tgt, "yac_interpolation_direct_execute");

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
  struct interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct interpolation_direct * direct = (struct interpolation_direct*)interp;

  yac_interpolation_exchange_execute_get(
    direct->src2tgt, tgt_field, "yac_interpolation_direct_execute_get");
}

static struct interpolation_type * yac_interpolation_direct_copy(
  struct interpolation_type * interp) {

  struct interpolation_direct * direct = (struct interpolation_direct*)interp;

  return
    yac_interpolation_direct_new_(
      direct->collection_size,
      yac_interpolation_exchange_copy(direct->src2tgt),
      yac_interpolation_buffer_copy(
        direct->src_data, 1, direct->collection_size));
}

static int yac_interpolation_direct_execute_test(
  struct interpolation_type * interp) {

  struct interpolation_direct * direct =
    (struct interpolation_direct*)interp;

  return
    yac_interpolation_exchange_test(
      direct->src2tgt, "yac_interpolation_direct_execute_test");
}

static void yac_interpolation_direct_delete(
  struct interpolation_type * interp) {

  if (interp == NULL) return;

  struct interpolation_direct * direct = (struct interpolation_direct*)interp;

  yac_interpolation_exchange_delete(
    direct->src2tgt, "yac_interpolation_direct_delete");
  yac_interpolation_buffer_free(&(direct->src_data));
  free(direct->src_field_buffer);
  free(direct);
}
