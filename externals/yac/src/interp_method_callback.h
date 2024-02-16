/**
 * @file interp_method_callback.h
 *
 * @copyright Copyright  (C)  2021 Moritz Hanke <hanke@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
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

#ifndef INTERP_METHOD_CALLBACK_H
#define INTERP_METHOD_CALLBACK_H

#include "interp_method.h"
#include "yac_interface.h"

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

#endif // INTERP_METHOD_CALLBACK_H
