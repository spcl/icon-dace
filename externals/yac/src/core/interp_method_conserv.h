// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_CONSERV_H
#define INTERP_METHOD_CONSERV_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

/** \example test_interp_method_conserv_parallel.c
 * A test for the parallel conservative‚ interpolation method.
 */

/**
 * normalisation options for conservative interpolation\n
 * (see SCRIP user manual for a more detailed description of the options)
 */
enum yac_interp_method_conserv_normalisation {
   YAC_INTERP_CONSERV_DESTAREA = 0, //!< interpolation values are normalised by the area of
                                    //!< the respective target cell
                                    //!< (this is the default option)
   YAC_INTERP_CONSERV_FRACAREA = 1, //!< interpolation values will be normalised by the area
                                    //!< of the respective target cell that is covered by
                                    //!< non-masked target cells
};

#define YAC_INTERP_CONSERV_ORDER_DEFAULT (1)
#define YAC_INTERP_CONSERV_ENFORCED_CONSERV_DEFAULT (0)
#define YAC_INTERP_CONSERV_PARTIAL_COVERAGE_DEFAULT (0)
#define YAC_INTERP_CONSERV_NORMALISATION_DEFAULT (0)

/**
 * constructor for a interpolation method of type interp_method_conserv
 * @param[in] order            1st or 2nd order remapping
 * @param[in] enforced_conserv enforce conservation by correcting truncation errors
 * @param[in] partial_coverage if == 0, target cells that are not completely
 *                             covered by non-masked source cells will not be
 *                             interpolated (this is the default option)\n
 *                             if != 0, for target cells that are only
 *                             partially covered by non-mask source cells the
 *                             interpolation value will be determined by the
 *                             contributions of the respective source cells
 * @param[in] normalisation    specifies how to do the normalisation
 */
struct interp_method * yac_interp_method_conserv_new(
  int order, int enforced_conserv, int partial_coverage,
  enum yac_interp_method_conserv_normalisation normalisation);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_CONSERV_H
