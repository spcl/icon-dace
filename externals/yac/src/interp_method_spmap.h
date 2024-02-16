/**
 * @file interp_method_spmap.h
 *
 * @copyright Copyright  (C)  2020 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
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

#ifndef INTERP_METHOD_SPMAP_H
#define INTERP_METHOD_SPMAP_H

#include "interp_method.h"

/** \example test_interp_method_spmap_parallel.c
 * A test for the parallel source point mapping interpolation method.
 */

#define YAC_DEFAULT_SPREAD_DISTANCE (0.0)
#define YAC_DEFAULT_MAX_SEARCH_DISTANCE (0.0)

enum yac_interp_spmap_weight_type {
  SPMAP_AVG  = 0, // simple average
  SPMAP_DIST = 1, // distance weighted
};

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
 * cells to generate the weights.
 * @param[in] spread_distance     spreading distance
 * @param[in] max_search_distance maximum search distance
 * @param[in] weight_type         weighting type
 * @remark the unit for the spread and maximum search distance is Radian
 */
struct interp_method * yac_interp_method_spmap_new(
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type);

#endif // INTERP_METHOD_SPMAP_H
