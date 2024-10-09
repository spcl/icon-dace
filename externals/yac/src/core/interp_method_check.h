// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_CHECK_H
#define INTERP_METHOD_CHECK_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

#define YAC_INTERP_CHECK_CONSTRUCTOR_KEY_DEFAULT ("")
#define YAC_INTERP_CHECK_DO_SEARCH_KEY_DEFAULT ("")

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

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_CHECK_H
