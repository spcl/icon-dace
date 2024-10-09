// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_INTERNAL_H
#define INTERP_METHOD_INTERNAL_H

#include "interp_method.h"
#include "dist_grid_internal.h"
#include "interp_grid_internal.h"
#include "interp_weights_internal.h"

struct interp_method_vtable {

   size_t (*do_search)(struct interp_method * method,
                       struct yac_interp_grid * grid,
                       size_t * tgt_points, size_t count,
                       struct yac_interp_weights * weights);
   void (*delete)(struct interp_method * method);
};

struct interp_method {

   struct interp_method_vtable *vtable;
};

#endif // INTERP_METHOD_INTERNAL_H
