// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERPOLATION_INTERNAL_H
#define INTERPOLATION_INTERNAL_H

#include "interpolation.h"

struct yac_interpolation_type;
struct yac_interpolation_type_vtable {
  int (*is_source)(struct yac_interpolation_type * interp);
  int (*is_target)(struct yac_interpolation_type * interp);
  void (*execute)(struct yac_interpolation_type * interp,
                  double *** src_fields, double *** src_frac_masks,
                  double ** tgt_field, double frac_mask_fallback_value,
                  double scale_factor, double scale_summand);
  void (*execute_put)(struct yac_interpolation_type * interp,
                      double *** src_fields, double *** src_frac_masks,
                      int is_target, double frac_mask_fallback_value,
                      double scale_factor, double scale_summand);
  void (*execute_get)(struct yac_interpolation_type * interp,
                      double ** tgt_field, double frac_mask_fallback_value,
                      double scale_factor, double scale_summand);
  void (*execute_get_async)(struct yac_interpolation_type * interp,
                            double ** tgt_field, double frac_mask_fallback_value,
                            double scale_factor, double scale_summand);
  int (*execute_put_test)(struct yac_interpolation_type * interp);
  int (*execute_get_test)(struct yac_interpolation_type * interp);
  void (*execute_wait)(struct yac_interpolation_type * interp);
  struct yac_interpolation_type * (*copy)(
    struct yac_interpolation_type * interp);
  void (*delete)(struct yac_interpolation_type * interp);
};
struct yac_interpolation_type {
  struct yac_interpolation_type_vtable const * const vtable;
};

#endif // INTERPOLATION_INTERNAL_H
