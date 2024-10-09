// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include "interp_method_internal.h"
#include "interp_method_fixed.h"

static size_t do_search_fixed(struct interp_method * method,
                              struct yac_interp_grid * interp_grid,
                              size_t * tgt_points, size_t count,
                              struct yac_interp_weights * weights);
static void delete_fixed(struct interp_method * method);

static struct interp_method_vtable
  interp_method_fixed_vtable = {
    .do_search = do_search_fixed,
    .delete = delete_fixed};

struct interp_method_fixed {

  struct interp_method_vtable * vtable;
  double value;
};

static size_t do_search_fixed (struct interp_method * method,
                               struct yac_interp_grid * interp_grid,
                               size_t * tgt_points, size_t count,
                               struct yac_interp_weights * weights) {

  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(interp_grid, tgt_points, count),
    .count = count};

  yac_interp_weights_add_fixed(
    weights, &tgts, ((struct interp_method_fixed *)method)->value);

  free(tgts.data);

  return count;
}

struct interp_method * yac_interp_method_fixed_new(double value) {

  struct interp_method_fixed * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_fixed_vtable;
  method->value = value;

  return (struct interp_method*)method;
}

static void delete_fixed(struct interp_method * method) {
  free(method);
}
