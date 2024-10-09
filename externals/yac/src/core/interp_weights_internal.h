// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_WEIGHTS_INTERNAL_H
#define INTERP_WEIGHTS_INTERNAL_H

#include "interp_weights.h"
#include "interpolation_internal.h"
#include "dist_grid_internal.h"

/**
 * adds targets that are to get a fixed value
 * @param[in] weights     interpolation weights
 * @param[in] tgts        targets that get a fixed value
 * @param[in] fixed_value fixed value that is to be assigned to
 *                        the provided targets
 */
void yac_interp_weights_add_fixed(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  double fixed_value);

/**
 * adds targets that are to get a weighted sum of source values
 * @param[in] weights         interpolation weights
 * @param[in] tgts            targets that get the sum
 * @param[in] num_src_per_tgt number of sources per target
 * @param[in] srcs            sources
 * @param[in] w               weights
 */
void yac_interp_weights_add_wsum(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * num_src_per_tgt, struct remote_point * srcs, double * w);

/**
 * adds targets that are to get a sum of source values
 * @param[in] weights         interpolation weights
 * @param[in] tgts            targets that get the weighted sum
 * @param[in] num_src_per_tgt number of sources per target
 * @param[in] srcs            sources
 */
void yac_interp_weights_add_sum(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * num_src_per_tgt, struct remote_point * srcs);

/**
 * adds targets that are to get single source values
 * @param[in] weights interpolation weights
 * @param[in] tgts    targets that get the value of a selected source point
 * @param[in] srcs    sources
 */
void yac_interp_weights_add_direct(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  struct remote_point * srcs);

/**
 * adds targets that are to get a weighted sum of source values
 * @param[in] weights                   interpolation weights
 * @param[in] tgts                      targets that get the weighted sum
 * @param[in] num_src_per_field_per_tgt number of sources per target per
 *                                      source field
 * @param[in] srcs_per_field            sources per source field
 * @param[in] w                         weights
 * @param[in] num_src_fields            number of input source fields
 * @remark num_src_per_field_per_tgt is a 1-D Array with the
 *         following data layout:\n
 *         num_src_per_field_per_tgt[tgts->count][num_src_fields]
 */
void yac_interp_weights_add_wsum_mf(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * num_src_per_field_per_tgt, struct remote_point ** srcs_per_field,
  double * w, size_t num_src_fields);

/**
 * adds targets that are to get a sum of source values
 * @param[in] weights                   interpolation weights
 * @param[in] tgts                      targets that get the sum
 * @param[in] num_src_per_field_per_tgt number of sources per target per
 *                                      source field
 * @param[in] srcs_per_field            sources per source field
 * @param[in] num_src_fields            number of input source fields
 * @remark num_src_per_field_per_tgt is a 1-D Array with the
 *         following data layout:\n
 *         num_src_per_field_per_tgt[tgts->count][num_src_fields]
 */
void yac_interp_weights_add_sum_mf(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * num_src_per_field_per_tgt, struct remote_point ** srcs_per_field,
  size_t num_src_fields);

/**
 * adds targets that are to get single source values
 * @param[in] weights           interpolation weights
 * @param[in] tgts              targets that get the value of a selected
 *                              source point
 * @param[in] src_field_indices source field indices of selected source point
 * @param[in] srcs_per_field    sources per source field
 * @param[in] num_src_fields    number of input source fields
 */
void yac_interp_weights_add_direct_mf(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * src_field_indices, struct remote_point ** srcs_per_field,
  size_t num_src_fields);

/**
 * adds targets whose stencil is the weighted sum of the copies of
 * existing stencils
 * @param[in] weights              interpolation weights
 * @param[in] tgts                 targets that get the weightes sum
 * @param[in] num_stencils_per_tgt number of stencils per target
 * @param[in] stencil_indices      indices of the stencils to be copied
 * @param[in] stencil_ranks        ranks of the processs owning the respective
 *                                 stencil
 * @param[in] w                    weights
 * @remark this call is collective
 */
void yac_interp_weights_wcopy_weights(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * num_stencils_per_tgt, size_t * stencil_indices,
  int * stencil_ranks, double * w);

/**
 * returns the MPI communicator used by weights
 * @param[in] weights interpololation weights
 * @return MPI communicator used by weights
 */
MPI_Comm yac_interp_weights_get_comm(struct yac_interp_weights * weights);

#endif // INTERP_WEIGHTS_INTERNAL_H

