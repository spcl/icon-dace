// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_WEIGHTS_H
#define INTERP_WEIGHTS_H

#include "yac_types.h"
#include "location.h"
#include "interpolation.h"

// YAC PUBLIC HEADER START

/** \example test_interp_weights_parallel.c
 * This contains some examples on how to use interp_weights.
 */

enum yac_interp_weights_reorder_type {
  YAC_MAPPING_ON_SRC, //!< weights will be appied at source processes
  YAC_MAPPING_ON_TGT, //!< weights will be applied at target processes
};

struct yac_interp_weights;

/**
 * Constructor for interpolation weights.
 * @param[in] comm           MPI communicator
 * @param[in] tgt_location   location of target field
 * @param[in] src_locations  locations of source fields
 * @param[in] num_src_fields number of source fields
 * @return interpolation weights
 */
struct yac_interp_weights * yac_interp_weights_new(
  MPI_Comm comm, enum yac_location tgt_location,
  enum yac_location * src_locations, size_t num_src_fields);

/**
 * writes interpolation weights to file
 * @param[in] weights       interpolation weights
 * @param[in] filename      file name
 * @param[in] src_grid_name name of the source grid
 * @param[in] tgt_grid_name name of the target grid
 * @param[in] src_grid_size global size of the source grid
 * @param[in] tgt_grid_size global size of the target grid
 * @remark this call is collective
 * @remark Global grid size argument can be either the global grid size or
 *         zero. If a valid global grid size was provided by at least one
 *         process, it will be added as a dimension to the weight file.
 */
void yac_interp_weights_write_to_file(
  struct yac_interp_weights * weights, char const * filename,
  char const * src_grid_name, char const * tgt_grid_name,
  size_t src_grid_size, size_t tgt_grid_size);

/**
 * generates an interpolation from interpolation weights
 * @param[in] weights                  interpolation weights
 * @param[in] reorder                  determines at which processes the
 *                                     weights are
 *                                     to be applied
 * @param[in] collection_size          collection size
 * @param[in] frac_mask_fallback_value fallback value for dynamic
 *                                     fractional masking
 * @param[in] scaling_factor           scaling factor
 * @param[in] scaling_summand          scaling summand
 * @return interpolation
 * @remark if frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE, dynamic
 *         fractional masking will be used
 * @remark all target field values, whose source points are not masked by
 *         the fractional mask, that receive a interpolation value, which is
 *         not a fixed value will by scaled by the following formula:\n
 *         y = scaling_factor * x + scaling_summand
 */
struct yac_interpolation * yac_interp_weights_get_interpolation(
  struct yac_interp_weights * weights,
  enum yac_interp_weights_reorder_type reorder,
  size_t collection_size, double frac_mask_fallback_value,
  double scaling_factor, double scaling_summand);

/**
 * returns the count of all target for which the weights contain a stencil
 * @param[in] weights interpolation weights
 * @return count of all targets in weights with a stencil
 */
size_t yac_interp_weights_get_interp_count(
  struct yac_interp_weights * weights);

/**
 * returns the global ids of all targets for which the weights contain a
 * stencil
 * @param[in] weights interpolation weights
 * @return global ids of all targets in weights with a stencil
 */
yac_int * yac_interp_weights_get_interp_tgt(
  struct yac_interp_weights * weights);

/**
 * Destructor for interpolation weights.
 * @param[inout] weights interpolation weights
 */
void yac_interp_weights_delete(struct yac_interp_weights * weights);

// YAC PUBLIC HEADER STOP

#endif // INTERP_WEIGHTS_H
