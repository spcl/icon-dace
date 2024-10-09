// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DUPLICATE_STENCILS_H
#define DUPLICATE_STENCILS_H

#include "interp_weights.h"

// YAC PUBLIC HEADER START

/** \example test_duplicate_stencils_parallel.c
 * Test for stencil duplication.
 */

/**
 * Duplicates stencils for specified target points in a yac_interp_weights
 * data structure.
 * @param[in,out] weights weights data structure
 * @param[in]     tgt_grid target grid
 * @param[in]     tgt_orig_global_id global ids of target points, whose
 *                                   interpolation stencils are to be
 *                                   duplicated
 * @param[in]     tgt_duplicated_idx local ids of target points, who are to
 *                                   interpolated using the duplicated stencils
 * @param[in]     nbr_duplicated     number of stencils that are to be
 *                                   duplicated
 * @param[in]     location           location of target points
 * @remark this operation is collective for all processes involved in the
 *         creation of the weights data structure
 */
void yac_duplicate_stencils(
  struct yac_interp_weights * weights, struct yac_basic_grid * tgt_grid,
  yac_int * tgt_orig_global_id, size_t * tgt_duplicated_idx,
  size_t nbr_duplicated, enum yac_location location);

// YAC PUBLIC HEADER STOP

#endif // DUPLICATE_STENCILS_H
