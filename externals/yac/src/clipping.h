/**
 * @file clipping.h
 * @brief Structs and interfaces for cell clipping
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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

#ifndef CLIPPING_H
#define CLIPPING_H

#include "grid_cell.h"

/** \example test_clipping.c
 * This contains some examples on how to use the \ref yac_cell_clipping
 * routine.
 */
/** \example test_lat_clipping.c
 * This contains some examples on how to use the yac_cell_lat_clipping
 * routine.
 */

struct yac_circle;

/**
  * \brief cell clipping to get the cells describing the intersections
  *
  * The routine takes (a list of) source cells and a target cell. It sets the
  * target cell data and does some further initialisation. Thus it needs to be
  * called for each new target cell intersection calculation
  *
  * The vertices of source and target cells can be either provided in a clockwise
  * or anticlockwise sense.
  *
  * @param[in] N              number of source cells
  * @param[in] source_cell    list of source cells
  * @param[in] target_cell    target cell
  * @param[in] overlap_buffer buffer for the overlaps between the target and
  *                           the source cells
  *
  * \remark source cells can be either convex or concave
  * \remark target cell has to be convex
  * \remark cells in overlap_buffer can be concave (even if source cell was convex)
  * \remark overlap_buffer must contain valid grid_cells (have to be initialised
  *         using \ref yac_init_grid_cell; initialisation have to be done only once,
  *         in consecutive calls, the cells can be reused with have to be
  *         reinitialised)
  *
 **/
void yac_cell_clipping (size_t N,
                        struct grid_cell * source_cell,
                        struct grid_cell target_cell,
                        struct grid_cell * overlap_buffer);

/**
  * \brief cell clipping to get the cells describing the intersections
  *
  * The routine takes (a list of) cells and two latitude bounds.
  *
  * @param[in] N              number of cells
  * @param[in] cells          list of cells
  * @param[in] lat_bounds     latitude bounds in radiant
  * @param[in] overlap_buffer buffer for the overlaps between the cells and
  *                           latitude band
  *
  * \remark cells in overlap_buffer can be concave
  * \remark overlap_buffer must contain valid grid_cells (have to be initialised
  *         using \ref yac_init_grid_cell; initialisation have to be done only once,
  *         in consecutive calls, the cells can be reused with have to be
  *         reinitialised)
  * \remark this routine is currently not being used within YAC but potentially
  *         used within the CDOs
  *
 **/
void yac_cell_lat_clipping (size_t N,
                            struct grid_cell * cells,
                            double lat_bounds[2],
                            struct grid_cell * overlap_buffer);

/** \example test_partial_areas.c
 * This contains examples on how to use \ref yac_compute_overlap_areas.
 */

/** \example test_compute_overlap_area.c
 * This contains examples on how to use \ref yac_compute_overlap_areas.
 */

/**
  * \brief calculates partial areas for all overlapping parts of the source
  *        cells with arbitrary target cells, this is required for conservative
  *        remapping.
  *
  * Some of the underlying concepts can be found in
  *
  * See e.g. Joseph O'Rourke, Computational Geometry in C, 2nd Ed., 1998
  *          Sec. 7.6 Intersections of Convex Polygons, page 253.
  *
  * The routine takes (a list of) source cells and a target cell. It determines the
  * clipping points of the intersection between a source and the target cells using
  * cell_clipping internally. In a second step areas are calculated for each
  * intersection described in the overlap cells. If a target cell is fully
  * covered by N source cells the N partial areas should add up to the area of
  * the target cell.
  *
  * @param[in]  N               number of source cells
  * @param[in]  source_cell     list of source cells
  * @param[in]  target_cell     target cell
  * @param[out] partial_areas   list of N partial weights, one weight for each
  *                             source-target intersection
  *
 **/
void yac_compute_overlap_areas (size_t N,
                                struct grid_cell * source_cell,
                                struct grid_cell target_cell,
                                double * partial_areas);

/**
  * \brief calculates partial areas for all overlapping parts of the source
  *        cells with arbitrary target cells, this is required for conservative
  *        remapping. In addition, the barycenter of each overlap is calculated.
  *
  * Some of the underlying concepts can be found in
  *
  * See e.g. Joseph O'Rourke, Computational Geometry in C, 2nd Ed., 1998
  *          Sec. 7.6 Intersections of Convex Polygons, page 253.
  *
  * The routine takes (a list of) source cells and a target cell. It determines the
  * clipping points of the intersection between a source and the target cells using
  * cell_clipping internally. In a second step areas are calculated for each
  * intersection described in the overlap cells. If a target cell is fully
  * covered by N source cells the N partial areas should add up to the area of
  * the target cell.
  *
  * @param[in]  N                   number of source cells
  * @param[in]  source_cell         list of source cells
  * @param[in]  target_cell         target cell
  * @param[out] overlap_areas       list of N partial weights, one weight for
  *                                 each source-target intersection
  * @param[out] overlap_barycenters coordinates of the barycenters of the
  *                                 overlap cell
  *
 **/
void yac_compute_overlap_info (size_t N,
                               struct grid_cell * source_cell,
                               struct grid_cell target_cell,
                               double * overlap_areas,
                               double (*overlap_barycenters)[3]);
/**
  * \brief correct interpolation weights
  *
  * Returns weights with a sum close to 1.0
  *
  * @param[in]  N                 number of source cells
  * @param[out] weight            list of N partial weights
  *
 **/
void yac_correct_weights (size_t N, double * weight);

/**
 * Generate a yac circle
 *
 * @param[in]  a             point a
 * @param[in]  b             point b
 * @param[in]  type          typ of edge connecting points a and b
 * @param[in]  edge_ordering ordering of the points a and b
 * @param[out] circle        yac circle
 */
void yac_circle_generate(
  double const * a, double const * b, enum yac_edge_type type,
  int edge_ordering, struct yac_circle * circle);

/**
 * Compare routine for yac circles (first by type and second parameters
 *                                  of the circle, if types are identical)
 *
 * @param[in] a yac circle a
 * @param[in] b yac circle b
 * @returns -1, 0, or 1 depending on the circles (lat circles always come last)
 */
int yac_circle_compare(void const * a, void const * b);

/**
 * Determines whether the provided yac circle contains the north pole
 * @param[in] circle yac circle
 * @return 1 if the circle contains the north pole, 0 otherwise
 */
int yac_circle_contains_north_pole(struct yac_circle * circle);

/**
 * Determines whether a given point is on the "inside"-side of a plane
 * defined by the given yac circle
 *
 * @param[in] point  point to be checked
 * @param[in] circle yac circle
 * @returns  0 if the point is not inside\n
 *           1 if the point is inside\n
 *           2 if the point is on the plane
 */
int yac_circle_point_is_inside(
  double const point[3], struct yac_circle * circle);

/** Compares the distances of two points to a given yac circle
 *
 * @param[in] a      point a
 * @param[in] b      point b
 * @param[in] circle yac circle
 * @returns -1 if point a is closer to the circle than point b\n
 *           0 if point a and b have the same distance to the circle\n
 *           1 if point b is closer to the circle than point a
 */
int yac_circle_compare_distances(
  double const a[3], double const b[3], struct yac_circle * circle);

#endif // CLIPPING_H
