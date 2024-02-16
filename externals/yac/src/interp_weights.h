/**
 * @file interp_weights.h
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

#ifndef INTERP_WEIGHTS_H
#define INTERP_WEIGHTS_H

#include "core/core.h"
#include "interp_grid.h"
#include "interpolation.h"

/** \example test_interp_weights_parallel.c
 * This contains some examples on how to use interp_weights.
 */

enum interp_weights_reorder_type {
  MAPPING_ON_SRC, //!< weights will be appied at source processes
  MAPPING_ON_TGT, //!< weights will be applied at target processes
};

struct interp_weights;

/**
 * Constructor for interpolation weights.
 * @param[in] comm           MPI communicator
 * @param[in] tgt_location   location of target field
 * @param[in] src_locations  locations of source fields
 * @param[in] num_src_fields number of source fields
 * @return interpolation weights
 */
struct interp_weights * yac_interp_weights_new(
  MPI_Comm comm, enum yac_location tgt_location,
  enum yac_location * src_locations, size_t num_src_fields);

/**
 * adds targets that are to get a fixed value
 * @param[in] weights     interpolation weights
 * @param[in] tgts        targets that get a fixed value
 * @param[in] fixed_value fixed value that is to be assigned to
 *                        the provided targets
 */
void yac_interp_weights_add_fixed(
  struct interp_weights * weights, struct remote_points * tgts,
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
  struct interp_weights * weights, struct remote_points * tgts,
  size_t * num_src_per_tgt, struct remote_point * srcs, double * w);

/**
 * adds targets that are to get a sum of source values
 * @param[in] weights         interpolation weights
 * @param[in] tgts            targets that get the weighted sum
 * @param[in] num_src_per_tgt number of sources per target
 * @param[in] srcs            sources
 */
void yac_interp_weights_add_sum(
  struct interp_weights * weights, struct remote_points * tgts,
  size_t * num_src_per_tgt, struct remote_point * srcs);

/**
 * adds targets that are to get single source values
 * @param[in] weights interpolation weights
 * @param[in] tgts    targets that get the value of a selected source point
 * @param[in] srcs    sources
 */
void yac_interp_weights_add_direct(
  struct interp_weights * weights, struct remote_points * tgts,
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
  struct interp_weights * weights, struct remote_points * tgts,
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
  struct interp_weights * weights, struct remote_points * tgts,
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
  struct interp_weights * weights, struct remote_points * tgts,
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
  struct interp_weights * weights, struct remote_points * tgts,
  size_t * num_stencils_per_tgt, size_t * stencil_indices,
  int * stencil_ranks, double * w);

/**
 * writes interpolation weights to file
 * @param[in] weights       interpolation weights
 * @param[in] filename      file name
 * @param[in] src_grid_name name of the source grid
 * @param[in] tgt_grid_name name of the target grid
 * @remark this call is collective
 */
void yac_interp_weights_write_to_file(
  struct interp_weights * weights, char const * filename,
  char const * src_grid_name, char const * tgt_grid_name);

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
struct interpolation * yac_interp_weights_get_interpolation(
  struct interp_weights * weights, enum interp_weights_reorder_type reorder,
  size_t collection_size, double frac_mask_fallback_value,
  double scaling_factor, double scaling_summand);

/**
 * returns the count of all target for which the weights contain a stencil
 * @param[in] weights interpolation weights
 * @return count of all targets in weights with a stencil
 */
size_t yac_interp_weights_get_interp_count(struct interp_weights * weights);

/**
 * returns the global ids of all targets for which the weights contain a
 * stencil
 * @param[in] weights interpolation weights
 * @return global ids of all targets in weights with a stencil
 */
yac_int * yac_interp_weights_get_interp_tgt(struct interp_weights * weights);

/**
 * Destructor for interpolation weights.
 * @param[inout] weights interpolation weights
 */
void yac_interp_weights_delete(struct interp_weights * weights);

#endif // INTERP_WEIGHTS_H
