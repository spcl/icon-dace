// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_CALLBACK_H
#define INTERP_METHOD_CALLBACK_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

#ifndef TYPEDEF_YAC_FUNC_COMPUTE_WEIGHTS
#define TYPEDEF_YAC_FUNC_COMPUTE_WEIGHTS

// Remark: make sure that this typedef is consistent with the one in yac.h
typedef void (*yac_func_compute_weights)(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);
#endif

#define YAC_INTERP_CALLBACK_COMPUTE_WEIGHTS_KEY_DEFAULT (NULL)

/** \example test_interp_method_callback_parallel.c
 * A test for the parallel callback interpolation method.
 */

/**
 * constructor for a interpolation method of type interp_method_callback
 * @param[in] compute_weights_callback pointer to routine that is used by YAC to
 *                                     compute the weights
 * @param[in] user_data                data pointer associated to the function
 *                                     pointer
 * @return returns a pointer to an interpolation method
 */
struct interp_method * yac_interp_method_callback_new(
  yac_func_compute_weights compute_weights_callback, void * user_data);

/**
 * sets a weight computation routine that can afterwards be retrieved by
 * \ref yac_interp_method_callback_get_compute_weights_callback
 * @param[in] compute_weights_callback pointer to a weight computation routine
 * @param[in] user_data                data pointer that will be passed to
 *                                     compute_weights_callback
 * @param[in] key                      key that can afterwards be used to
 *                                     retrieve the provided pointers
 */
void yac_interp_method_callback_add_compute_weights_callback(
  yac_func_compute_weights compute_weights_callback,
  void * user_data, char const * key);

/**
 * retrieves a compute_weights_callback pointer that was set by
 * \ref yac_interp_method_callback_add_compute_weights_callback
 * @param[in]  key                      key that identifies the pointers that
 *                                      are to be retrieved
 * @param[out] compute_weights_callback pointer to a weight computation routine
 * @param[out] user_data                data pointer that was provided together
 *                                      with the function pointer
 */
void yac_interp_method_callback_get_compute_weights_callback(
  char const * key, yac_func_compute_weights * compute_weights_callback,
  void ** user_data);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_CALLBACK_H
