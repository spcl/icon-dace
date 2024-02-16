/**
 * @file interp_method_check.c
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

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include "string.h"
#include "utils.h"
#include "ensure_array_size.h"
#include "interp_method_check.h"

static size_t do_search_check(struct interp_method * method,
                              struct interp_grid * interp_grid,
                              size_t * tgt_points, size_t count,
                              struct interp_weights * weights);
static void delete_check(struct interp_method * method);

static struct interp_method_vtable
  interp_method_check_vtable = {
    .do_search = do_search_check,
    .delete = delete_check};

struct interp_method_check {

  struct interp_method_vtable * vtable;
  func_do_search do_search_callback;
  void * do_search_user_data;
};

static size_t do_search_check (struct interp_method * method,
                               struct interp_grid * interp_grid,
                               size_t * tgt_points, size_t count,
                               struct interp_weights * weights) {

  struct interp_method_check * method_check =
    (struct interp_method_check *)method;

  yac_int * tgt_global_ids = xmalloc(count * sizeof(*tgt_global_ids));
  coordinate_pointer tgt_coordinates =
    xmalloc(count * sizeof(*tgt_coordinates));

  yac_interp_grid_get_tgt_global_ids(
    interp_grid, tgt_points, count, tgt_global_ids);
  yac_interp_grid_get_tgt_coordinates(
    interp_grid, tgt_points, count, tgt_coordinates);

  if (method_check->do_search_callback != NULL)
    method_check->do_search_callback(
      tgt_global_ids, (const_coordinate_pointer)tgt_coordinates, count,
      method_check->do_search_user_data);

  free(tgt_coordinates);
  free(tgt_global_ids);

  return 0;
}

struct interp_method * yac_interp_method_check_new(
  func_constructor constructor_callback,
  void * constructor_user_data,
  func_do_search do_search_callback,
  void * do_search_user_data) {

  struct interp_method_check * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_check_vtable;
  method->do_search_callback = do_search_callback;
  method->do_search_user_data = do_search_user_data;

  if (constructor_callback != NULL)
    constructor_callback(constructor_user_data);

  return (struct interp_method*)method;
}

static void delete_check(struct interp_method * method) {
  free(method);
}

enum callback_type {
  CONSTRUCTOR,
  DO_SEARCH,
};
typedef void (*func_dummy)(void);
static struct {
  struct {
    func_dummy callback;
    void * user_data;
  } value;
  enum callback_type type;
  char * key;
} * callback_lookup_table = NULL;
static size_t callback_lookup_table_array_size = 0;
static size_t callback_lookup_table_size = 0;

static void interp_method_check_add_callback(
  func_dummy callback, enum callback_type type, void * user_data,
  char const * key) {

  // ensure that the lookup table does already contain an entry with the same
  // key and type
  for (size_t i = 0; i < callback_lookup_table_size; ++i)
    YAC_ASSERT_F(
      strcmp(callback_lookup_table[i].key, key) ||
      (callback_lookup_table[i].type != type),
      "ERROR(interp_method_check_add_callback): "
      "key \"%s\" has been set for the same callback type (%d)",
      key, (int)type)

  ENSURE_ARRAY_SIZE(callback_lookup_table, callback_lookup_table_array_size,
                    callback_lookup_table_size + 1);

  callback_lookup_table[callback_lookup_table_size].key =
    xmalloc((strlen(key)+1) * sizeof(*key));
  strcpy(callback_lookup_table[callback_lookup_table_size].key, key);
  callback_lookup_table[callback_lookup_table_size].type = type;
  callback_lookup_table[callback_lookup_table_size].value.callback = callback;
  callback_lookup_table[callback_lookup_table_size].value.user_data = user_data;
  callback_lookup_table_size++;
}

static void interp_method_get_callback(
  char const * key, func_dummy * callback, enum callback_type type,
  void ** user_data) {

  *callback = NULL;
  *user_data = NULL;
  for (size_t i = 0;
      (i < callback_lookup_table_size) &&
      (*callback == NULL) && (*user_data == NULL); ++i) {
    if ((!strcmp(callback_lookup_table[i].key, key)) &&
        (callback_lookup_table[i].type == type)) {
      *callback = callback_lookup_table[i].value.callback;
      *user_data = callback_lookup_table[i].value.user_data;
      return;
    }
  }
}

void yac_interp_method_check_add_constructor_callback(
  func_constructor constructor_callback, void * user_data, char const * key) {

  interp_method_check_add_callback(
    (func_dummy)constructor_callback, CONSTRUCTOR, user_data, key);
}

void yac_interp_method_check_get_constructor_callback(
  char const * key, func_constructor * constructor_callback,
  void ** user_data) {

  interp_method_get_callback(
    key, (func_dummy*)constructor_callback, CONSTRUCTOR, user_data);
}

void yac_interp_method_check_add_do_search_callback(
  func_do_search do_search_callback, void * user_data, char const * key) {

  interp_method_check_add_callback(
    (func_dummy)do_search_callback, DO_SEARCH, user_data, key);
}

void yac_interp_method_check_get_do_search_callback(
  char const * key, func_do_search * do_search_callback, void ** user_data) {

  interp_method_get_callback(
    key, (func_dummy*)do_search_callback, DO_SEARCH, user_data);
}
