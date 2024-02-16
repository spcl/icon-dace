/**
 * @file interpolation_utils.h
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

#ifndef INTERPOLATION_UTILS_H
#define INTERPOLATION_UTILS_H

#include <math.h>

#include "core/ppm_xfuncs.h"
#include "yac_mpi.h"

static inline void compute_tgt_field_wgt(
  double const * restrict ** src_fields,
  double const * restrict ** src_frac_masks,
  double const * restrict * remote_src_fields,
  double const * restrict * remote_src_frac_masks,
  double * restrict * tgt_field,
  size_t const * restrict tgt_pos,
  size_t tgt_count, size_t const * restrict num_src_per_tgt,
  double const * restrict weights,
  size_t const * restrict src_field_idx,
  size_t const * restrict src_idx,
  size_t num_src_fields, size_t collection_size,
  double frac_mask_fallback_value,
  double scale_factor, double scale_summand) {

  // return if there is nothing to do
  if (tgt_count == 0) return;

#define FRAC_MASK_TOL (1e-12)

#define COMPUTE_TGT_FIELD_WGT_FRAC_(TGT_POS, WEIGHT, SCALE) \
  { \
    for (size_t l = 0; l < collection_size; ++l) { \
      double const * restrict * curr_local_field_data = \
        src_fields?src_fields[l]:NULL; \
      double const * restrict * curr_local_frac_mask_data = \
        src_frac_masks?src_frac_masks[l]:NULL; \
      double const * restrict curr_remote_field_data = \
        remote_src_fields[l * num_src_fields]; \
      double const * restrict curr_remote_frac_mask_data = \
        remote_src_frac_masks[l * num_src_fields]; \
      double * restrict curr_tgt_field = tgt_field[l]; \
      for (size_t i = 0, k = 0; i < tgt_count; ++i) { \
        size_t curr_num_src_per_tgt = num_src_per_tgt[i]; \
        double result = 0.0; \
        double frac_weight_sum = 0.0; \
        for (size_t j = 0; j < curr_num_src_per_tgt; ++j, ++k) { \
          double const * restrict frac_mask_data; \
          double const * restrict src_field_data; \
          if (src_field_idx[k] == SIZE_MAX) { \
            frac_mask_data = curr_remote_frac_mask_data; \
            src_field_data = curr_remote_field_data; \
          } else { \
            frac_mask_data = \
              curr_local_frac_mask_data[src_field_idx[k]]; \
            src_field_data = curr_local_field_data[src_field_idx[k]]; \
          } \
          result += src_field_data[src_idx[k]] * (WEIGHT); \
          frac_weight_sum += frac_mask_data[src_idx[k]] * (WEIGHT); \
        } \
        curr_tgt_field[(TGT_POS)] = \
          (fabs(frac_weight_sum) > FRAC_MASK_TOL)? \
            (SCALE(result / frac_weight_sum)): \
            frac_mask_fallback_value; \
      } \
    } \
  }

#define COMPUTE_TGT_FIELD_WGT_NOFRAC_(TGT_POS, WEIGHT, SCALE) \
  { \
    for (size_t l = 0; l < collection_size; ++l) { \
      double const * restrict * curr_local_field_data = \
        src_fields?src_fields[l]:NULL; \
      double const * restrict curr_remote_field_data = \
        remote_src_fields[l * num_src_fields]; \
      double * restrict curr_tgt_field = tgt_field[l]; \
      for (size_t i = 0, k = 0; i < tgt_count; ++i) { \
        size_t curr_num_src_per_tgt = num_src_per_tgt[i]; \
        double result = 0.0; \
        for (size_t j = 0; j < curr_num_src_per_tgt; ++j, ++k) { \
          double const * restrict src_field_data; \
          if (src_field_idx[k] == SIZE_MAX) { \
            src_field_data = curr_remote_field_data; \
          } else { \
            src_field_data = curr_local_field_data[src_field_idx[k]]; \
          } \
          result += src_field_data[src_idx[k]] * (WEIGHT); \
        } \
        curr_tgt_field[(TGT_POS)] = SCALE(result); \
      } \
    } \
  }

#define COMPUTE_TGT_FIELD_WGT(COMPUTE, SCALE) \
  { \
    if (weights != NULL) { \
      if (tgt_pos != NULL) COMPUTE(tgt_pos[i], weights[k], SCALE) \
      else                 COMPUTE(i, weights[k], SCALE) \
    } else { \
      if (tgt_pos != NULL) COMPUTE(tgt_pos[i], 1.0, SCALE) \
      else                 COMPUTE(i, 1.0, SCALE) \
    } \
  }

#define COMPUTE_TGT_FIELD_WGT_FRAC(SCALE) \
  COMPUTE_TGT_FIELD_WGT(COMPUTE_TGT_FIELD_WGT_FRAC_, SCALE)

#define COMPUTE_TGT_FIELD_WGT_NOFRAC(SCALE) \
  COMPUTE_TGT_FIELD_WGT(COMPUTE_TGT_FIELD_WGT_NOFRAC_, SCALE)

#define NO_SCALING(RESULT) (RESULT)
#define MULT(RESULT) ((RESULT) * scale_factor)
#define ADD(RESULT) ((RESULT) + scale_summand)
#define MULT_ADD(RESULT) ((RESULT) * scale_factor + scale_summand)

#define COMPUTE_FIELD(COMPUTE_FIELD_FRAC, COMPUTE_FIELD_NOFRAC) \
  { \
    if (frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE) { \
      if (scale_factor == 1.0) { \
        if (scale_summand == 0.0) COMPUTE_FIELD_FRAC(NO_SCALING) \
        else                      COMPUTE_FIELD_FRAC(ADD) \
      } else { \
        if (scale_summand == 0.0) COMPUTE_FIELD_FRAC(MULT) \
        else                      COMPUTE_FIELD_FRAC(MULT_ADD) \
      } \
    } else { \
      if (scale_factor == 1.0) { \
        if (scale_summand == 0.0) COMPUTE_FIELD_NOFRAC(NO_SCALING) \
        else                      COMPUTE_FIELD_NOFRAC(ADD) \
      } else { \
        if (scale_summand == 0.0) COMPUTE_FIELD_NOFRAC(MULT) \
        else                      COMPUTE_FIELD_NOFRAC(MULT_ADD) \
      } \
    } \
  }

  COMPUTE_FIELD(COMPUTE_TGT_FIELD_WGT_FRAC, COMPUTE_TGT_FIELD_WGT_NOFRAC)

#undef COMPUTE_TGT_FIELD_WGT
#undef COMPUTE_TGT_FIELD_WGT_FRAC
#undef COMPUTE_TGT_FIELD_WGT_NOFRAC
#undef COMPUTE_TGT_FIELD_WGT_FRAC_
#undef COMPUTE_TGT_FIELD_WGT_NOFRAC_
}

static inline void compute_tgt_field(
  double const * restrict ** src_fields,
  double const * restrict ** src_frac_masks,
  double * restrict * tgt_field,
  size_t * restrict tgt_buffer_sizes,
  size_t num_src_fields, size_t collection_size,
  double frac_mask_fallback_value,
  double scale_factor, double scale_summand)  {

#define COMPUTE_TGT_FIELD_FRAC(SCALE) \
  { \
    for (size_t i = 0; i < collection_size; ++i) { \
      for (size_t j = 0; j < num_src_fields; ++j) { \
        memcpy(tgt_field[i * num_src_fields + j],  \
                src_fields[i][j], tgt_buffer_sizes[j]); \
        for (size_t k = 0, offset = 0; offset < tgt_buffer_sizes[j]; \
            ++k, offset += sizeof(***src_fields)) { \
          if (src_frac_masks[i][j][k] != 0.0) \
            tgt_field[i * num_src_fields + j][k] = \
              SCALE( \
                tgt_field[i * num_src_fields + j][k] / \
                src_frac_masks[i][j][k]); \
          else \
            tgt_field[i * num_src_fields + j][k] = \
              frac_mask_fallback_value; \
        } \
      } \
    } \
  }

#define COMPUTE_TGT_FIELD_NOFRAC(SCALE) \
  { \
    for (size_t i = 0; i < collection_size; ++i) { \
      for (size_t j = 0; j < num_src_fields; ++j) { \
        memcpy(tgt_field[i * num_src_fields + j],  \
                src_fields[i][j], tgt_buffer_sizes[j]); \
        for (size_t k = 0, offset = 0; offset < tgt_buffer_sizes[j]; \
            ++k, offset += sizeof(***src_fields)) \
          tgt_field[i * num_src_fields + j][k] = \
            SCALE(tgt_field[i * num_src_fields + j][k]); \
      } \
    } \
  }

  COMPUTE_FIELD(COMPUTE_TGT_FIELD_FRAC, COMPUTE_TGT_FIELD_NOFRAC)

#undef COMPUTE_TGT_FIELD_NOFRAC
#undef COMPUTE_TGT_FIELD_FRAC
}

#undef MULT_ADD
#undef ADD
#undef MULT
#undef NO_SCALING
#undef COMPUTE_FIELD
#undef FRAC_MASK_TOL

#define CHECK_WITH_FRAC_MASK(ROUTINE) \
  YAC_ASSERT_F( \
    with_frac_mask == (frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE), \
    "ERROR(%s) " \
    "with_frac_mask does not match value provided to constructor\n" \
    "(frac_mask_fallback_value = %lf with_frac_mask = %d " \
    "with_frac_mask(constructor) %d)", (ROUTINE), \
    frac_mask_fallback_value, \
    frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE, with_frac_mask)

enum yac_interpolation_buffer_type {
  SEND_BUFFER,
  RECV_BUFFER,
};

struct yac_interpolation_buffer {

  double ** buffer;
  size_t * buffer_sizes;
};

struct yac_interpolation_buffer yac_interpolation_buffer_init(
  Xt_redist * redists, size_t num_fields, size_t collection_size,
  enum yac_interpolation_buffer_type type);

struct yac_interpolation_buffer yac_interpolation_buffer_init_2(
  Xt_redist * redists, size_t * min_buffer_sizes, size_t num_fields,
  size_t collection_size, enum yac_interpolation_buffer_type type);

struct yac_interpolation_buffer yac_interpolation_buffer_copy(
  struct yac_interpolation_buffer buffer, size_t num_fields,
  size_t collection_size);

void yac_interpolation_buffer_free(struct yac_interpolation_buffer * buffer);

#endif // INTERPOLATION_UTILS_H
