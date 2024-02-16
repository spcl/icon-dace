/**
 * @file interp_method_check.h
 *
 * @copyright Copyright  (C)  2020 Moritz Hanke <hanke@dkrz.de>
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

#ifndef INTERP_METHOD_CHECK_H
#define INTERP_METHOD_CHECK_H

#include "interp_method.h"

typedef void (*func_constructor)(void * user_data);
typedef void (*func_do_search)(yac_int const * global_ids,
  double const (*coordinates_xyz)[3], size_t count, void * user_data);

/**
 * constructor for a interpolation method of type interp_method_check
 * @param[in] constructor_callback  pointer to routine that is called in
 *                                  \ref yac_interp_method_check_new
 * @param[in] constructor_user_data pointer passed to constructor_callback when
 *                                  it is called
 * @param[in] do_search_callback    pointer to routine that is to be called when
 *                                  the do_search routine of this
 *                                  interp_method is called
 * @param[in] do_search_user_data   pointer passed to do_search_callback when it
 *                                  is called
 * @return returns a pointer to an interpolation method
 */
struct interp_method * yac_interp_method_check_new(
  func_constructor constructor_callback,
  void * constructor_user_data,
  func_do_search do_search_callback,
  void * do_search_user_data);

/**
 * sets a constructor_callback that can afterwards be retrieved by
 * \ref yac_interp_method_check_get_constructor_callback
 * @param[in] constructor_callback pointer to a constructor_callback routine
 * @param[in] user_data            pointer to user data
 * @param[in] key                  string that can afterwards be used to
 *                                 retrieve the pointer
 */
void yac_interp_method_check_add_constructor_callback(
  func_constructor constructor_callback, void * user_data, char const * key);

/**
 * retrieves a constructor_callback pointer that was set by
 * \ref yac_interp_method_check_add_constructor_callback
 * @param[in] key                 string to identify a pointer that was
 *                                previously set
 * @param[out] constructor_callback pointer to constructor_callback routine
 * @param[out] user_data            pointer to user data
 */
void yac_interp_method_check_get_constructor_callback(
  char const * key, func_constructor * constructor_callback, void ** user_data);

/**
 * sets a do_search_callback that can afterwards be retrieved by
 * \ref yac_interp_method_check_get_do_search_callback
 * @param[in] do_search_callback pointer to a do_search_callback routine
 * @param[in] user_data            pointer to user data
 * @param[in] key                string that can afterwards be used to retrieve
 *                               the pointer
 */
void yac_interp_method_check_add_do_search_callback(
  func_do_search do_search_callback, void * user_data, char const * key);

/**
 * retrieves a do_search_callback pointer that was set by
 * \ref yac_interp_method_check_add_do_search_callback
 * @param[in] key                 string to identify a pointer that was
 *                                previously set
 * @param[out] do_search_callback pointer to do_search_callback routine
 * @param[out] user_data          pointer to user data
 */
void yac_interp_method_check_get_do_search_callback(
  char const * key, func_do_search * do_search_callback, void ** user_data);

#endif // INTERP_METHOD_CHECK_H
