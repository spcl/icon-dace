/**
 * @file interpolation_sum_mvp_at_tgt.c
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

#include "interpolation_sum_mvp_at_tgt.h"
#include "utils.h"
#include "yaxt.h"
#include "interpolation_utils.h"
#include "interpolation_exchange.h"

static int yac_interpolation_sum_mvp_at_tgt_is_source(
  struct interpolation_type * interp);
static int yac_interpolation_sum_mvp_at_tgt_is_target(
  struct interpolation_type * interp);
static void yac_interpolation_sum_mvp_at_tgt_execute(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static void yac_interpolation_sum_mvp_at_tgt_execute_put(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks,
  int is_target, double frac_mask_fallback_value,
  double scale_factor, double scale_summand);
static void yac_interpolation_sum_mvp_at_tgt_execute_get(
  struct interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand);
static int yac_interpolation_sum_mvp_at_tgt_execute_test(
  struct interpolation_type * interp);
static struct interpolation_type * yac_interpolation_sum_mvp_at_tgt_copy(
  struct interpolation_type * interp);
static void yac_interpolation_sum_mvp_at_tgt_delete(
  struct interpolation_type * interp);

static struct interpolation_type_vtable const interpolation_sum_mvp_at_tgt_vtable = {
  .is_source        = yac_interpolation_sum_mvp_at_tgt_is_source,
  .is_target        = yac_interpolation_sum_mvp_at_tgt_is_target,
  .execute          = yac_interpolation_sum_mvp_at_tgt_execute,
  .execute_put      = yac_interpolation_sum_mvp_at_tgt_execute_put,
  .execute_get      = yac_interpolation_sum_mvp_at_tgt_execute_get,
  .execute_test     = yac_interpolation_sum_mvp_at_tgt_execute_test,
  .copy             = yac_interpolation_sum_mvp_at_tgt_copy,
  .delete           = yac_interpolation_sum_mvp_at_tgt_delete,
};

struct interpolation_sum_mvp_at_tgt {

  struct interpolation_type_vtable const * vtable;

  size_t collection_size;
  int with_frac_mask;

  /* data flow:
   *  put:
   *      I. source processes pack their source data into a buffer, in order
   *         to be able to asynchronously process the put operation
   *         - from buffer: src_fields (+src_frac_masks) provided by user
   *         - to buffer: src_send_buffer
   *     II. source processes send source data to target processes, such
   *         that each target process can compute its own target field
   *         - send buffer: src_send_buffer
   *         - recv buffer: src_recv_buffer
   *         - exchange: src2tgt
   *  get:
   *      I. target processes receive source data from source processes
   *         - send buffer: src_send_buffer
   *         - recv buffer: src_recv_buffer
   *         - exchange: src2tgt
   *     II. target process comput their target field
   *         - from buffer: src_recv_buffer
   *         - to buffer: tgt_field provied by the user
   */

  struct yac_interpolation_buffer src_send_data;
  struct yac_interpolation_buffer src_recv_data;
  struct yac_interpolation_exchange * src2tgt;

  double *** src_fields;
  double *** src_frac_masks;
  double ** src_fields_buffer;
  size_t * tgt_pos;
  size_t tgt_count;
  size_t * num_src_per_tgt;
  double * weights;
  size_t * src_field_idx;
  size_t * src_idx;
  size_t num_src_fields;
  int is_source;
  int is_target;

  int * ref_count;
};

static struct interpolation_type * yac_interpolation_sum_mvp_at_tgt_new_(
  size_t collection_size,
  struct yac_interpolation_buffer src_send_data,
  struct yac_interpolation_buffer src_recv_data,
  struct yac_interpolation_exchange * src2tgt,
  size_t * tgt_pos,  size_t tgt_count, size_t * num_src_per_tgt,
  double * weights, size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, int with_frac_mask, int * ref_count) {

  struct interpolation_sum_mvp_at_tgt * mvp_at_tgt =
    xmalloc(1 * sizeof(*mvp_at_tgt));

  mvp_at_tgt->vtable = &interpolation_sum_mvp_at_tgt_vtable;
  mvp_at_tgt->collection_size = collection_size;
  mvp_at_tgt->with_frac_mask = with_frac_mask;
  mvp_at_tgt->src_send_data = src_send_data;
  mvp_at_tgt->src_recv_data = src_recv_data;
  mvp_at_tgt->src2tgt = src2tgt;
  mvp_at_tgt->src_fields =
    xmalloc(collection_size * sizeof(*(mvp_at_tgt->src_fields)));
  for (size_t i = 0; i < collection_size; ++i)
    mvp_at_tgt->src_fields[i] = mvp_at_tgt->src_send_data.buffer + i * num_src_fields;
  if (with_frac_mask) {
    mvp_at_tgt->src_frac_masks =
      xmalloc(collection_size * sizeof(*(mvp_at_tgt->src_frac_masks)));
    for (size_t i = 0; i < collection_size; ++i)
      mvp_at_tgt->src_frac_masks[i] =
        mvp_at_tgt->src_send_data.buffer + (collection_size + i) * num_src_fields;
  } else {
    mvp_at_tgt->src_frac_masks = NULL;
  }
  mvp_at_tgt->src_fields_buffer =
    xcalloc((with_frac_mask?2*collection_size:collection_size) *
            num_src_fields, sizeof(*(mvp_at_tgt->src_fields_buffer)));
  mvp_at_tgt->tgt_pos = tgt_pos;
  mvp_at_tgt->tgt_count = tgt_count;
  mvp_at_tgt->num_src_per_tgt = num_src_per_tgt;
  mvp_at_tgt->weights = weights;
  mvp_at_tgt->src_field_idx = src_field_idx;
  mvp_at_tgt->src_idx = src_idx;
  mvp_at_tgt->num_src_fields = num_src_fields;
  mvp_at_tgt->is_source =
    yac_interpolation_exchange_is_source(mvp_at_tgt->src2tgt);
  mvp_at_tgt->is_target = tgt_count > 0;

  mvp_at_tgt->ref_count =
    (ref_count == NULL)?xcalloc(1, sizeof(*mvp_at_tgt->ref_count)):ref_count;
  ++*(mvp_at_tgt->ref_count);

  return (struct interpolation_type *)mvp_at_tgt;
}

struct interpolation_type * yac_interpolation_sum_mvp_at_tgt_new(
  size_t collection_size, Xt_redist * src_redists, size_t * tgt_pos,
  size_t tgt_count, size_t * num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, int with_frac_mask) {

  size_t total_num_src = 0;
  for (size_t i = 0; i < tgt_count; ++i) total_num_src += num_src_per_tgt[i];

  // in case a process is both source and target, the send buffer is also used
  // to copy the send field (and send_frac_mask) in the put for the get
  // --> we have to make sure that the buffer contains all locally required
  //     data
  size_t * min_buffer_sizes =
    xcalloc(num_src_fields, sizeof(*min_buffer_sizes));
  for (size_t i = 0, offset = 0; i < tgt_count; ++i) {
    size_t curr_num_src = num_src_per_tgt[i];
    for (size_t j = 0; j < curr_num_src; ++j, ++offset) {
      size_t curr_src_field_idx = src_field_idx[offset];
      if (curr_src_field_idx == SIZE_MAX) continue;
      size_t curr_src_extent = (src_idx[offset] + 1) * sizeof(double);
      if (min_buffer_sizes[curr_src_field_idx] < curr_src_extent)
        min_buffer_sizes[curr_src_field_idx] = curr_src_extent;
    }
  }

  struct interpolation_type * interp =
    yac_interpolation_sum_mvp_at_tgt_new_(
      collection_size,
      yac_interpolation_buffer_init_2(
        src_redists, min_buffer_sizes, num_src_fields,
        with_frac_mask?2*collection_size:collection_size, SEND_BUFFER),
      yac_interpolation_buffer_init(
        src_redists, num_src_fields,
        with_frac_mask?2*collection_size:collection_size, RECV_BUFFER),
      yac_interpolation_exchange_new(
        src_redists, num_src_fields,
        collection_size, with_frac_mask, "source to target"),
      COPY_DATA(tgt_pos, tgt_count), tgt_count,
      COPY_DATA(num_src_per_tgt, tgt_count),
      (weights != NULL)?COPY_DATA(weights, total_num_src):NULL,
      COPY_DATA(src_field_idx, total_num_src),
      COPY_DATA(src_idx, total_num_src),
      num_src_fields, with_frac_mask, NULL);

  free(min_buffer_sizes);
  return interp;
}

static int yac_interpolation_sum_mvp_at_tgt_is_source(
  struct interpolation_type * interp) {

  struct interpolation_sum_mvp_at_tgt * sum_mvp_at_tgt =
    (struct interpolation_sum_mvp_at_tgt *)interp;

  return sum_mvp_at_tgt->is_source;
}

static int yac_interpolation_sum_mvp_at_tgt_is_target(
  struct interpolation_type * interp) {

  struct interpolation_sum_mvp_at_tgt * sum_mvp_at_tgt =
    (struct interpolation_sum_mvp_at_tgt *)interp;

  return sum_mvp_at_tgt->is_target;
}

static void yac_interpolation_sum_mvp_at_tgt_execute(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct interpolation_sum_mvp_at_tgt * sum_mvp_at_tgt =
    (struct interpolation_sum_mvp_at_tgt *)interp;

  int with_frac_mask = sum_mvp_at_tgt->with_frac_mask;
  CHECK_WITH_FRAC_MASK("yac_interpolation_sum_mvp_at_tgt_execute")

  size_t collection_size = sum_mvp_at_tgt->collection_size;
  size_t num_src_fields = sum_mvp_at_tgt->num_src_fields;
  double ** temp_src_fields = sum_mvp_at_tgt->src_fields_buffer;
  double ** src_recv_buffer = sum_mvp_at_tgt->src_recv_data.buffer;
  if (sum_mvp_at_tgt->is_source) {
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < num_src_fields; ++j)
        temp_src_fields[i * num_src_fields + j] = src_fields[i][j];
    if (with_frac_mask)
      for (size_t i = 0; i < collection_size; ++i)
        for (size_t j = 0; j < num_src_fields; ++j)
          temp_src_fields
            [i * num_src_fields + j + collection_size * num_src_fields] =
            src_frac_masks[i][j];
  }

  yac_interpolation_exchange_wait(
    sum_mvp_at_tgt->src2tgt,
    "yac_interpolation_sum_mvp_at_tgt_execute");

  // send source points to targets
  yac_interpolation_exchange_execute(
    sum_mvp_at_tgt->src2tgt, (double const **)temp_src_fields,
    src_recv_buffer, "yac_interpolation_sum_mvp_at_tgt_execute");

  compute_tgt_field_wgt(
    (double const * restrict **)src_fields,
    (double const * restrict **)(with_frac_mask?src_frac_masks:NULL),
    (double const * restrict *)src_recv_buffer,
    (double const * restrict *)(
      with_frac_mask?(src_recv_buffer + collection_size * num_src_fields):NULL),
    tgt_field, sum_mvp_at_tgt->tgt_pos, sum_mvp_at_tgt->tgt_count,
    sum_mvp_at_tgt->num_src_per_tgt, sum_mvp_at_tgt->weights,
    sum_mvp_at_tgt->src_field_idx, sum_mvp_at_tgt->src_idx,
    num_src_fields, collection_size, frac_mask_fallback_value,
    scale_factor, scale_summand);
}

static void yac_interpolation_sum_mvp_at_tgt_execute_put(
  struct interpolation_type * interp,
  double *** src_fields, double *** src_frac_masks, int is_target,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct interpolation_sum_mvp_at_tgt * sum_mvp_at_tgt =
    (struct interpolation_sum_mvp_at_tgt *)interp;

  double ** src_send_buffer = NULL;

  if (sum_mvp_at_tgt->is_source) {

    // wait until previous put is completed
    yac_interpolation_exchange_wait(
      sum_mvp_at_tgt->src2tgt,
      "yac_interpolation_sum_mvp_at_tgt_execute_put");

    int with_frac_mask = sum_mvp_at_tgt->with_frac_mask;
    CHECK_WITH_FRAC_MASK("yac_interpolation_sum_mvp_at_tgt_execute_put")

    size_t collection_size = sum_mvp_at_tgt->collection_size;
    size_t num_src_fields = sum_mvp_at_tgt->num_src_fields;
    size_t * src_send_buffer_sizes = sum_mvp_at_tgt->src_send_data.buffer_sizes;
    src_send_buffer = sum_mvp_at_tgt->src_send_data.buffer;
    if (sum_mvp_at_tgt->is_source)
      for (size_t i = 0; i < collection_size; ++i)
        for (size_t j = 0; j < num_src_fields; ++j)
          memcpy(src_send_buffer[i * num_src_fields + j], src_fields[i][j],
                 src_send_buffer_sizes[j]);
    if (with_frac_mask) {
        for (size_t i = 0; i < collection_size; ++i)
          for (size_t j = 0; j < num_src_fields; ++j)
            memcpy(
              src_send_buffer
                [i * num_src_fields + j + collection_size * num_src_fields],
                src_frac_masks[i][j], src_send_buffer_sizes[j]);
    }
  }

  yac_interpolation_exchange_execute_put(
    sum_mvp_at_tgt->src2tgt, (double const **)src_send_buffer,
    "yac_interpolation_sum_mvp_at_tgt_execute_put");
}

static void yac_interpolation_sum_mvp_at_tgt_execute_get(
  struct interpolation_type * interp, double ** tgt_field,
  double frac_mask_fallback_value, double scale_factor, double scale_summand) {

  struct interpolation_sum_mvp_at_tgt * sum_mvp_at_tgt =
    (struct interpolation_sum_mvp_at_tgt*)interp;

  int with_frac_mask = sum_mvp_at_tgt->with_frac_mask;
  CHECK_WITH_FRAC_MASK("yac_interpolation_sum_mvp_at_tgt_execute_get")

  size_t collection_size = sum_mvp_at_tgt->collection_size;
  size_t num_src_fields = sum_mvp_at_tgt->num_src_fields;
  double ** src_recv_buffer = sum_mvp_at_tgt->src_recv_data.buffer;

  // receive source field data
  yac_interpolation_exchange_execute_get(
    sum_mvp_at_tgt->src2tgt, src_recv_buffer,
    "yac_interpolation_sum_mvp_at_tgt_execute_get");

  compute_tgt_field_wgt(
    (double const * restrict **)(sum_mvp_at_tgt->src_fields),
    (double const * restrict **)(
      with_frac_mask?sum_mvp_at_tgt->src_frac_masks:NULL),
    (double const * restrict *)src_recv_buffer,
    (double const * restrict *)(
      with_frac_mask?(src_recv_buffer + collection_size * num_src_fields):NULL),
    tgt_field, sum_mvp_at_tgt->tgt_pos, sum_mvp_at_tgt->tgt_count,
    sum_mvp_at_tgt->num_src_per_tgt,
    sum_mvp_at_tgt->weights, sum_mvp_at_tgt->src_field_idx,
    sum_mvp_at_tgt->src_idx, num_src_fields, collection_size,
    frac_mask_fallback_value, scale_factor, scale_summand);
}

static int yac_interpolation_sum_mvp_at_tgt_execute_test(
  struct interpolation_type * interp) {

  struct interpolation_sum_mvp_at_tgt * sum_mvp_at_tgt =
    (struct interpolation_sum_mvp_at_tgt*)interp;

  return
    yac_interpolation_exchange_test(
      sum_mvp_at_tgt->src2tgt,
      "yac_interpolation_sum_mvp_at_tgt_execute_test");
}

static struct interpolation_type * yac_interpolation_sum_mvp_at_tgt_copy(
  struct interpolation_type * interp) {

  struct interpolation_sum_mvp_at_tgt * sum_mvp_at_tgt =
    (struct interpolation_sum_mvp_at_tgt*)interp;

  return
    yac_interpolation_sum_mvp_at_tgt_new_(
      sum_mvp_at_tgt->collection_size,
      yac_interpolation_buffer_copy(
        sum_mvp_at_tgt->src_send_data, sum_mvp_at_tgt->num_src_fields,
        sum_mvp_at_tgt->collection_size),
      yac_interpolation_buffer_copy(
        sum_mvp_at_tgt->src_recv_data, sum_mvp_at_tgt->num_src_fields,
        sum_mvp_at_tgt->collection_size),
      yac_interpolation_exchange_copy(sum_mvp_at_tgt->src2tgt),
      sum_mvp_at_tgt->tgt_pos,  sum_mvp_at_tgt->tgt_count,
      sum_mvp_at_tgt->num_src_per_tgt, sum_mvp_at_tgt->weights,
      sum_mvp_at_tgt->src_field_idx, sum_mvp_at_tgt->src_idx,
      sum_mvp_at_tgt->num_src_fields, sum_mvp_at_tgt->with_frac_mask,
      sum_mvp_at_tgt->ref_count);
}

static void yac_interpolation_sum_mvp_at_tgt_delete(
  struct interpolation_type * interp) {

  if (interp == NULL) return;

  struct interpolation_sum_mvp_at_tgt * sum_mvp_at_tgt =
    (struct interpolation_sum_mvp_at_tgt*)interp;

  free(sum_mvp_at_tgt->src_fields);
  if (sum_mvp_at_tgt->src_frac_masks != NULL)
    free(sum_mvp_at_tgt->src_frac_masks);
  yac_interpolation_exchange_delete(
    sum_mvp_at_tgt->src2tgt, "yac_interpolation_sum_mvp_at_tgt_delete");
  yac_interpolation_buffer_free(&(sum_mvp_at_tgt->src_send_data));
  yac_interpolation_buffer_free(&(sum_mvp_at_tgt->src_recv_data));
  free(sum_mvp_at_tgt->src_fields_buffer);

  if (!--(*(sum_mvp_at_tgt->ref_count))) {
    free(sum_mvp_at_tgt->tgt_pos);
    free(sum_mvp_at_tgt->src_idx);
    free(sum_mvp_at_tgt->src_field_idx);
    free(sum_mvp_at_tgt->weights);
    free(sum_mvp_at_tgt->num_src_per_tgt);
    free(sum_mvp_at_tgt->ref_count);
  }
  free(sum_mvp_at_tgt);
}
