/**
 * @file sphere_part.h
 *
 * @copyright Copyright  (C)  2014 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
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

#ifndef SPHERE_PART_H
#define SPHERE_PART_H

/**
 * \file sphere_part.h
 * \brief algorithm for searching cells and points on a grid
 *
 * \ref yac_point_sphere_part_search_new generates a tree structure, which
 * makes it easy to look for points. A documentation of the respective
 * algorithm can be found at \ref sphere_part_docu.
 */
 
/**
 * \page sphere_part_docu Sphere Partitioning Algorithm
 *
 * The following describes how the Sphere Partitioning algorithm generates a
 * tree data structure for a given set of polygons on a sphere.
 * This tree structure allows to easily search all cells in the given data
 * set that overlaps with another given cell or point (set of cells or points).
 *
 * \section sphere_part_tree_gen Generation of tree structure
 * to partition a set of polygons on the sphere the following data structure can be constructed:
 * - c = center of sphere
 * - S = set of all points
 *
 * \subsection recurspart_sec Description for a routine that generates the tree
 * recurspart(P, v, t)
 * -# Set p to balance point of P
 * -# Compute great circle L with base plane collinear to |c p| and orthogonal to v
 * -# Partition data
 *    -# Set I to subset of S intersecting L
 *    -# Set T to subset of S/I with positive scalar product to normal vector of L
 *       - Ti = recurspart(T, norm(L)) if |T| > threshold t else Ti = list(T)
 *    -# Set U to subset of S/I with negative scalar product to normal vector of L
 *       - Ui = recurspart(U, norm(L)) if |U| > threshold t else Ui = list(U)
 * -# return node(list(I), Ti, Ui, L, alpha=angle over L containing I)
 *
 * \subsection recuspart_observ_sec Observations
 * - I, T, U and form a partition of S
 * - first L is orthogonal to equator, i.e. a longitude circle
 * - initial call: recurspart(S, (0,0,1), t)
 *
 * \section sphere_part_search Searching in the tree structure
 *
 * \subsection search_sec Description for a routine that searches for a list of cells
 * search(n, p)
 * -# if n is leaf
 *    -# P = search_list(n, p)
 * -# else
 *    -# P = {}
 *    -# if p in n.alpha
 *       - P = P united search_list(n.I, p)
 *    -# if p * norm(n.L) > 0
 *       - P = P united search(n.Ti, p)
 *    -# else if p * norm(n.L) < 0
 *       - P = P united search(n.Ui, p)
 * -# return P
 *
 * \subsection search_observ_sec Observations
 * - returns list of matching polygons
 */

#include "grid.h"
#include "geometry.h"

/** \example test_bnd_sphere_part.c
 * This contains a test of the bounding circle sphere part algorithm.
 */

struct point_sphere_part_search;

struct point_sphere_part_search * yac_point_sphere_part_search_new (
  size_t num_points, coordinate_pointer coordinates_xyz, yac_int const * ids);

struct point_sphere_part_search * yac_point_sphere_part_search_mask_new (
  size_t num_points, coordinate_pointer coordinates_xyz,
  yac_int const * ids, int const * mask);

void yac_delete_point_sphere_part_search(
  struct point_sphere_part_search * search);

/**
 * This routine does a nearest neighbour search between the points provided to
 * this routine and the matching yac_point_sphere_part_search_new call.
 */
void yac_point_sphere_part_search_NN(struct point_sphere_part_search * search,
                                     size_t num_points,
                                     double (*coordinates_xyz)[3],
                                     double * cos_angles,
                                     double (**result_coordinates_xyz)[3],
                                     size_t * result_coordinates_xyz_array_size,
                                     size_t ** local_point_ids,
                                     size_t * local_point_ids_array_size,
                                     size_t * num_local_point_ids);

/**
 * This routine does a n nearest neighbour search between the points provided to
 * this routine and the matching yac_point_sphere_part_search_new call.
 */
void yac_point_sphere_part_search_NNN(struct point_sphere_part_search * search,
                                      size_t num_points,
                                      double (*coordinates_xyz)[3], size_t n,
                                      double ** cos_angles,
                                      size_t * cos_angles_array_size,
                                      double (**result_coordinates_xyz)[3],
                                      size_t * result_coordinates_xyz_array_size,
                                      size_t ** local_point_ids,
                                      size_t * local_point_ids_array_size,
                                      size_t * num_local_point_ids);

/**
 * This routine returns true if the provided point_sphere_part_search contains
 * a point that is within the provided bounding circle.
 */
int yac_point_sphere_part_search_bnd_circle_contains_points(
  struct point_sphere_part_search * search, struct bounding_circle circle);

struct bnd_sphere_part_search;
struct bnd_sphere_part_search * yac_bnd_sphere_part_search_new(
  struct bounding_circle * circles, size_t num_circles);
void yac_bnd_sphere_part_search_delete(struct bnd_sphere_part_search * search);
void yac_bnd_sphere_part_search_do_point_search(
  struct bnd_sphere_part_search * search, coordinate_pointer coordinates_xyz,
  size_t count, size_t ** cells, size_t * num_cells_per_coordinate);
void yac_bnd_sphere_part_search_do_bnd_circle_search(
  struct bnd_sphere_part_search * search, struct bounding_circle * bnd_circles,
  size_t count, size_t ** cells, size_t * num_cells_per_bnd_circle);

#endif // SPHERE_PART_H

