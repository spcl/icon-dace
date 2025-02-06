// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERPOLATION_SUM_MVP_AT_TGT_H
#define INTERPOLATION_SUM_MVP_AT_TGT_H

#include "interpolation_internal.h"

struct yac_interpolation_type * yac_interpolation_sum_mvp_at_tgt_new(
  size_t collection_size, Xt_redist * src_redists, size_t * tgt_pos,
  size_t tgt_count, size_t * num_src_per_tgt, double * weights,
  size_t * src_field_idx, size_t * src_idx,
  size_t num_src_fields, int with_frac_mask);

#endif // INTERPOLATION_SUM_MVP_AT_TGT_H
