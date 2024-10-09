// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_SPMAP_H
#define INTERP_METHOD_SPMAP_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

/** \example test_interp_method_spmap_parallel.c
 * A test for the parallel source point mapping interpolation method.
 */

enum yac_interp_spmap_weight_type {
  YAC_INTERP_SPMAP_AVG  = 0, // simple average
  YAC_INTERP_SPMAP_DIST = 1, // distance weighted
};

enum yac_interp_spmap_scale_type {
  YAC_INTERP_SPMAP_NONE       = 0, //!< weights are not scaled
  YAC_INTERP_SPMAP_SRCAREA    = 1, //!< weights are multiplied by
                                   //!< the area of the associated source cell
  YAC_INTERP_SPMAP_INVTGTAREA = 2, //!< weights are muliplied by
                                   //!< the inverse of the area of the
                                   //!< associated target cell
  YAC_INTERP_SPMAP_FRACAREA  = 3,  //!< weights are multiplied by
                                   //!< the area of the associated source cell
                                   //!< and the inverse of the area of the
                                   //!< associated target cell
};

#define YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT (0.0)
#define YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT (0.0)
#define YAC_INTERP_SPMAP_WEIGHTED_DEFAULT (YAC_INTERP_SPMAP_AVG)
#define YAC_INTERP_SPMAP_SCALE_DEFAULT (YAC_INTERP_SPMAP_NONE)
#define YAC_INTERP_SPMAP_SRC_SPHERE_RADIUS_DEFAULT (1.0)
#define YAC_INTERP_SPMAP_TGT_SPHERE_RADIUS_DEFAULT (1.0)

/**
 * Constructor for a interpolation method of type interp_method_spmap\n
 * This method searches for each unmasked source point the closest unmasked
 * target point.\n
 * If the maximum search distance is > 0.0, only target points that are within
 * this distance from the source points are being considered.
 * If spread_distance is > 0.0, the method uses the previously found target
 * points as a starting point. Around each starting point, a bounding circle is
 * generated. Afterwards for each starting point all target cells whose bounding
 * circles intersect with the generated one are put into a list. Out of this
 * list all target cells connected directly or indirectly through other cells
 * from this list to the target cell of the starting are selected for the
 * interpolation. Then a weighting method is applied to the selected target
 * cells to generate the weights. Afterwards the weights are scaled.
 * @param[in] spread_distance     spreading distance
 * @param[in] max_search_distance maximum search distance
 * @param[in] weight_type         weighting type
 * @param[in] scale_type          scaling type
 * @param[in] src_sphere_radius   sphere radius used for
 *                                source cell area computation
 * @param[in] tgt_sphere_radius   sphere radius used for
 *                                target cell area computation
 * @remark the unit for the spread and maximum search distance is Radian
 */
struct interp_method * yac_interp_method_spmap_new(
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type,
  enum yac_interp_spmap_scale_type scale_type,
  double src_sphere_radius, double tgt_sphere_radius);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_SPMAP_H
