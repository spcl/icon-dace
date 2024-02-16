/**
 * @file check_overlap.c
 * @brief Set of functions to determine an overlap between any two cells 
 *
 * Compared to \ref check_overlap.c these routine take care of how vertex
 * points of a cell are connected, either by great circles or along the loxodrome.
 *
 * very interesting literature:
 * - http://geospatialmethods.org/spheres/GCIntersect.html
 *
 * Note: Not all functions are documented by Doxygen. See the source code
 * and \ref geometry.h for further details.
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

#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "geometry.h"
#include "grid_cell.h"

static double const tol = 1.0e-12;

// returns the direction of the given edge along the equator
// switching the points a and b of the edge give the following result:
// a->b == -(b-a)
// (exeption: if the edge is a meridian the result in always 0)
static inline int edge_direction(double * a, double * b) {

   double cross_ab_2 = a[0] * b[1] - a[1] * b[0];

   return (cross_ab_2 > - tol) - (cross_ab_2 < tol);
}

static int point_in_cell_generic(
  double point_coords[3], struct grid_cell cell) {

   double second_point[3];

   // if the point is on the pole
   if (fabs(fabs(point_coords[2]) - 1.0) < tol) {
      second_point[0] = 1;
      second_point[1] = 0;
      second_point[2] = 0;
   } else if (point_coords[2] > 0) {
      second_point[0] = 0;
      second_point[1] = 0;
      second_point[2] = -1;
   } else {
      second_point[0] = 0;
      second_point[1] = 0;
      second_point[2] = 1;
   }

   int through_cell_corner;
   int edge_crossings;
   through_cell_corner = 0;
   edge_crossings = 0;

   // for all edges of cell
   for (size_t i = 0; i < cell.num_corners; ++i) {

      double * a = cell.coordinates_xyz[i];
      double * b = cell.coordinates_xyz[(i+1)%cell.num_corners];

      int ret_value;
      double p[3], q[3];

      ret_value = yac_intersect_vec(cell.edge_type[i], a, b, LON_CIRCLE_EDGE,
                                    point_coords, second_point, p, q);

      // if both edges do not intersect
      if ((ret_value == -1) || (ret_value == 0)) continue;

      // if p is not the intersection point of both edges
      if ((ret_value & ((1 << 0) + (1 << 2))) != (1 << 0) + (1 << 2)) {

         // if q is the intersection point of both edges
         if ((ret_value & ((1 << 1) + (1 << 3))) == (1 << 1) + (1 << 3))
            p[0] = q[0], p[1] = q[1], p[2] = q[2];
         else
            continue;
      }

      // if the intersection is the point itself
      if (points_are_identically(point_coords, p))
         return 1;

      //----------
      // in case the point edge goes through a endpoint of a cell edge
      // there are three special cases that need to be taken care of:
      // 1: the cell is concave an the point edge only touched an inner corner of the cell -> does not switch in/outside
      // 2: the point edge only touched an outer corner of the cell -> does not switch in/outside
      // 3: the point edge entered/left the cell through a corner -> switch in/outside
      //---------
      if ((points_are_identically(p, a)) ||
          (points_are_identically(p, b))) {

         // check the direction of the cell edge (edge can have a positive or
         // negative direction in terms of longitudes)
         // each cell endpoint is checked twice (for each adjacent edge)
         // If the direction for two adjacent edges sharing an endpoint that
         // is crossed by the point edge, then the related endpoint
         // intersection can be classified as case 3.
         // case 1 and 2 result in no change to through_cell_corner
         // case 3 result in an in/decrease of through_cell_corner by 2

         through_cell_corner += edge_direction(a, b);
      } else {

         //standard case -> we crossed an edge
         edge_crossings++;
      }
   } // (j = 0; j < cell.num_corners; ++j)

   edge_crossings += through_cell_corner / 2;

   // if we have an odd number of edge crossings -> point was inside cell
   // or the point was directly on an longitude edge
   return ( (edge_crossings & 1) || (through_cell_corner & 1) );
}

static int point_in_cell_unstruct_triangle(
  double point_coords[3], struct grid_cell cell) {

  double * a = cell.coordinates_xyz[0];
  double * b = cell.coordinates_xyz[1];
  double * c = cell.coordinates_xyz[2];

  double cross_edges[3][3];

  crossproduct_kahan(a, b, &(cross_edges[0][0]));
  crossproduct_kahan(b, c, &(cross_edges[1][0]));
  crossproduct_kahan(c, a, &(cross_edges[2][0]));

  double sq_sin_edge[3] = {
    cross_edges[0][0] * cross_edges[0][0] +
    cross_edges[0][1] * cross_edges[0][1] +
    cross_edges[0][2] * cross_edges[0][2],
    cross_edges[1][0] * cross_edges[1][0] +
    cross_edges[1][1] * cross_edges[1][1] +
    cross_edges[1][2] * cross_edges[1][2],
    cross_edges[2][0] * cross_edges[2][0] +
    cross_edges[2][1] * cross_edges[2][1] +
    cross_edges[2][2] * cross_edges[2][2]};

  // if all edges of the triangle have a length of at least yac_angle_tol
  if ((sq_sin_edge[0] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[1] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[2] > (yac_angle_tol * yac_angle_tol))) {

    normalise_vector(&(cross_edges[0][0]));
    normalise_vector(&(cross_edges[1][0]));
    normalise_vector(&(cross_edges[2][0]));

    double sin_point_edge[3] = {
      cross_edges[0][0] * point_coords[0] +
      cross_edges[0][1] * point_coords[1] +
      cross_edges[0][2] * point_coords[2],
      cross_edges[1][0] * point_coords[0] +
      cross_edges[1][1] * point_coords[1] +
      cross_edges[1][2] * point_coords[2],
      cross_edges[2][0] * point_coords[0] +
      cross_edges[2][1] * point_coords[1] +
      cross_edges[2][2] * point_coords[2]};

    int cell_corner_order =
      cross_edges[0][0] * c[0] + cross_edges[0][1] * c[1] + cross_edges[0][2] * c[2] >= 0;

    if (cell_corner_order)
      return (sin_point_edge[0] > - yac_angle_tol) &&
             (sin_point_edge[1] > - yac_angle_tol) &&
             (sin_point_edge[2] > - yac_angle_tol);
    else
      return (sin_point_edge[0] < yac_angle_tol) &&
             (sin_point_edge[1] < yac_angle_tol) &&
             (sin_point_edge[2] < yac_angle_tol);

  // if the length of the first edge is shorter than yac_angle_tol
  } else if (sq_sin_edge[0] <= (yac_angle_tol * yac_angle_tol)) {

    // if the length of the second edge is also shorter than yac_angle_tol
    // -> triangle is a point
    if (sq_sin_edge[1] <= (yac_angle_tol * yac_angle_tol))
      return
        compare_angles(get_vector_angle_2(a, point_coords), SIN_COS_TOL) <= 0;
    // -> triangle is an edge
    else {

      normalise_vector(&(cross_edges[1][0]));
      double sin_point_edge_1 =
        cross_edges[1][0] * point_coords[0] +
        cross_edges[1][1] * point_coords[1] +
        cross_edges[1][2] * point_coords[2];

      struct sin_cos_angle edge_angle =
        sum_angles_no_check(get_vector_angle_2(a, c), SIN_COS_TOL);
      struct sin_cos_angle ap_angle = get_vector_angle_2(a, point_coords);
      struct sin_cos_angle cp_angle = get_vector_angle_2(c, point_coords);
      // point is on the plane of the edge and between b and c
      return (fabs(sin_point_edge_1) <= yac_angle_tol) &&
             (compare_angles(ap_angle, edge_angle) <= 0) &&
             (compare_angles(cp_angle, edge_angle) <= 0);
    }

  // -> triangle is an edge
  } else {

    normalise_vector(&(cross_edges[0][0]));

    double sin_point_edge_0 =
      cross_edges[0][0] * point_coords[0] +
      cross_edges[0][1] * point_coords[1] +
      cross_edges[0][2] * point_coords[2];

    struct sin_cos_angle edge_angle =
      sum_angles_no_check(get_vector_angle_2(a, b), SIN_COS_TOL);
    struct sin_cos_angle ap_angle = get_vector_angle_2(a, point_coords);
    struct sin_cos_angle bp_angle = get_vector_angle_2(b, point_coords);
    // point is on the plane of the edge and between b and c
    return (fabs(sin_point_edge_0) <= yac_angle_tol) &&
           (compare_angles(ap_angle, edge_angle) <= 0) &&
           (compare_angles(bp_angle, edge_angle) <= 0);
  }
}

// static int any_point_in_cell_unstruct_triangle(
  // double * point_coords, size_t num_points, struct grid_cell cell) {

  // double * a = cell.coordinates_xyz[0];
  // double * b = cell.coordinates_xyz[1];
  // double * c = cell.coordinates_xyz[2];

  // double cross_edges[3][3];

  // crossproduct_kahan(a, b, &(cross_edges[0][0]));
  // crossproduct_kahan(b, c, &(cross_edges[1][0]));
  // crossproduct_kahan(c, a, &(cross_edges[2][0]));

  // double sq_sin_edge[3] = {
    // cross_edges[0][0] * cross_edges[0][0] +
    // cross_edges[0][1] * cross_edges[0][1] +
    // cross_edges[0][2] * cross_edges[0][2],
    // cross_edges[1][0] * cross_edges[1][0] +
    // cross_edges[1][1] * cross_edges[1][1] +
    // cross_edges[1][2] * cross_edges[1][2],
    // cross_edges[2][0] * cross_edges[2][0] +
    // cross_edges[2][1] * cross_edges[2][1] +
    // cross_edges[2][2] * cross_edges[2][2]};

  // normalise_vector(&(cross_edges[0][0]));
  // normalise_vector(&(cross_edges[1][0]));
  // normalise_vector(&(cross_edges[2][0]));

  // // if all edges of the triangle have a length of at least yac_angle_tol
  // if ((sq_sin_edge[0] > (yac_angle_tol * yac_angle_tol)) &&
      // (sq_sin_edge[1] > (yac_angle_tol * yac_angle_tol)) &&
      // (sq_sin_edge[2] > (yac_angle_tol * yac_angle_tol))) {

    // int cell_corner_order =
      // cross_edges[0][0] * c[0] + cross_edges[0][1] * c[1] + cross_edges[0][2] * c[2] >= 0;

    // if (cell_corner_order) {
      // for (size_t i = 0; i < num_points; ++i) {
        // double sin_point_edge[3] = {
          // cross_edges[0][0] * point_coords[0+3*i] +
          // cross_edges[0][1] * point_coords[1+3*i] +
          // cross_edges[0][2] * point_coords[2+3*i],
          // cross_edges[1][0] * point_coords[0+3*i] +
          // cross_edges[1][1] * point_coords[1+3*i] +
          // cross_edges[1][2] * point_coords[2+3*i],
          // cross_edges[2][0] * point_coords[0+3*i] +
          // cross_edges[2][1] * point_coords[1+3*i] +
          // cross_edges[2][2] * point_coords[2+3*i]};
        // if ((sin_point_edge[0] > - yac_angle_tol) &&
            // (sin_point_edge[1] > - yac_angle_tol) &&
            // (sin_point_edge[2] > - yac_angle_tol)) return 1;
      // }
    // } else {
      // for (size_t i = 0; i < num_points; ++i) {
        // double sin_point_edge[3] = {
          // cross_edges[0][0] * point_coords[0+3*i] +
          // cross_edges[0][1] * point_coords[1+3*i] +
          // cross_edges[0][2] * point_coords[2+3*i],
          // cross_edges[1][0] * point_coords[0+3*i] +
          // cross_edges[1][1] * point_coords[1+3*i] +
          // cross_edges[1][2] * point_coords[2+3*i],
          // cross_edges[2][0] * point_coords[0+3*i] +
          // cross_edges[2][1] * point_coords[1+3*i] +
          // cross_edges[2][2] * point_coords[2+3*i]};
        // if ((sin_point_edge[0] < yac_angle_tol) &&
            // (sin_point_edge[1] < yac_angle_tol) &&
            // (sin_point_edge[2] < yac_angle_tol)) return 1;
      // }
    // }

  // // if the length of the first edge is shorter than yac_angle_tol
  // } else if (sq_sin_edge[0] <= (yac_angle_tol * yac_angle_tol)) {

    // // if the length of the second edge is also shorter than yac_angle_tol
    // // -> triangle is a point
    // if (sq_sin_edge[1] <= (yac_angle_tol * yac_angle_tol))
      // return
        // compare_angles(get_vector_angle_2(a, point_coords), SIN_COS_TOL) <= 0;
    // // -> triangle is an edge
    // else {
      // struct sin_cos_angle edge_angle =
        // sum_angles_no_check(get_vector_angle_2(a, c), SIN_COS_TOL);

      // for (size_t i = 0; i < num_points; ++i) {
        // struct sin_cos_angle ap_angle = get_vector_angle_2(a, point_coords+3*i);
        // struct sin_cos_angle cp_angle = get_vector_angle_2(c, point_coords+3*i);
        // double sin_point_edge = cross_edges[1][0] * point_coords[0+3*i] +
                                // cross_edges[1][1] * point_coords[1+3*i] +
                                // cross_edges[1][2] * point_coords[2+3*i];
        // // point is on the plane of the edge and between b and c
        // if ((fabs(sin_point_edge) <= yac_angle_tol) &&
            // (compare_angles(ap_angle, edge_angle) <= 0) &&
            // (compare_angles(cp_angle, edge_angle) <= 0)) return 1;
      // }
    // }

  // // -> triangle is an edge
  // } else {
    // struct sin_cos_angle edge_angle =
      // sum_angles_no_check(get_vector_angle_2(a, b), SIN_COS_TOL);

    // for (size_t i = 0; i < num_points; ++i) {
      // struct sin_cos_angle ap_angle = get_vector_angle_2(a, point_coords+3*i);
      // struct sin_cos_angle bp_angle = get_vector_angle_2(b, point_coords+3*i);
      // double sin_point_edge = cross_edges[0][0] * point_coords[0+3*i] +
                              // cross_edges[0][1] * point_coords[1+3*i] +
                              // cross_edges[0][2] * point_coords[2+3*i];
      // // point is on the plane of the edge and between b and c
      // if  ((fabs(sin_point_edge) <= yac_angle_tol) &&
           // (compare_angles(ap_angle, edge_angle) <= 0) &&
           // (compare_angles(bp_angle, edge_angle) <= 0)) return 1;
    // }
  // }

  // return 0;
// }

static inline double dotproduct(double a[], double b[]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static int point_in_cell_unstruct_quad(
  double point_coords[3], struct grid_cell cell) {

  double * a = cell.coordinates_xyz[0];
  double * b = cell.coordinates_xyz[1];
  double * c = cell.coordinates_xyz[2];
  double * d = cell.coordinates_xyz[3];

  double cross[4][3];

  crossproduct_kahan(a, b, &(cross[0][0]));
  crossproduct_kahan(b, c, &(cross[1][0]));
  crossproduct_kahan(c, d, &(cross[2][0]));
  crossproduct_kahan(d, a, &(cross[3][0]));

  double sq_sin_edge[4] = {
    cross[0][0]*cross[0][0]+cross[0][1]*cross[0][1]+cross[0][2]*cross[0][2],
    cross[1][0]*cross[1][0]+cross[1][1]*cross[1][1]+cross[1][2]*cross[1][2],
    cross[2][0]*cross[2][0]+cross[2][1]*cross[2][1]+cross[2][2]*cross[2][2],
    cross[3][0]*cross[3][0]+cross[3][1]*cross[3][1]+cross[3][2]*cross[3][2]};

  // if all edges of the quad have a length of at least yac_angle_tol
  if ((sq_sin_edge[0] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[1] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[2] > (yac_angle_tol * yac_angle_tol)) &&
      (sq_sin_edge[3] > (yac_angle_tol * yac_angle_tol))) {

    normalise_vector(&(cross[0][0]));
    normalise_vector(&(cross[1][0]));
    normalise_vector(&(cross[2][0]));
    normalise_vector(&(cross[3][0]));

    int tmp = dotproduct(c, cross[0]) >= 0.0;
    tmp += dotproduct(d, cross[0]) >= 0.0;
    tmp += dotproduct(a, cross[1]) >= 0.0;
    tmp += dotproduct(d, cross[1]) >= 0.0;
    tmp += dotproduct(a, cross[2]) >= 0.0;
    tmp += dotproduct(b, cross[2]) >= 0.0;
    tmp += dotproduct(b, cross[3]) >= 0.0;
    tmp += dotproduct(c, cross[3]) >= 0.0;
    // check whether the quad is convex
    if (!(tmp & 7)) {
      if (tmp)
        return
          (dotproduct(point_coords, cross[0]) + yac_angle_tol >=  0.0) &&
          (dotproduct(point_coords, cross[1]) + yac_angle_tol >=  0.0) &&
          (dotproduct(point_coords, cross[2]) + yac_angle_tol >=  0.0) &&
          (dotproduct(point_coords, cross[3]) + yac_angle_tol >=  0.0);
      else
        return
          (dotproduct(point_coords, cross[0]) - yac_angle_tol <= -0.0) &&
          (dotproduct(point_coords, cross[1]) - yac_angle_tol <= -0.0) &&
          (dotproduct(point_coords, cross[2]) - yac_angle_tol <= -0.0) &&
          (dotproduct(point_coords, cross[3]) - yac_angle_tol <= -0.0);
    }
  }

  return point_in_cell_generic(point_coords, cell);
}

// static int any_point_in_cell_unstruct_quad(
  // double * point_coords, size_t num_points, struct grid_cell cell) {

  // double * a = cell.coordinates_xyz[0];
  // double * b = cell.coordinates_xyz[1];
  // double * c = cell.coordinates_xyz[2];
  // double * d = cell.coordinates_xyz[3];

  // double cross[4][3];

  // crossproduct_kahan(a, b, &(cross[0][0]));
  // crossproduct_kahan(b, c, &(cross[1][0]));
  // crossproduct_kahan(c, d, &(cross[2][0]));
  // crossproduct_kahan(d, a, &(cross[3][0]));

  // double sq_sin_edge[4] = {
    // cross[0][0]*cross[0][0]+cross[0][1]*cross[0][1]+cross[0][2]*cross[0][2],
    // cross[1][0]*cross[1][0]+cross[1][1]*cross[1][1]+cross[1][2]*cross[1][2],
    // cross[2][0]*cross[2][0]+cross[2][1]*cross[2][1]+cross[2][2]*cross[2][2],
    // cross[3][0]*cross[3][0]+cross[3][1]*cross[3][1]+cross[3][2]*cross[3][2]};

  // normalise_vector(&(cross[0][0]));
  // normalise_vector(&(cross[1][0]));
  // normalise_vector(&(cross[2][0]));
  // normalise_vector(&(cross[3][0]));

  // // if all edges of the quad have a length of at least yac_angle_tol
  // if ((sq_sin_edge[0] > (yac_angle_tol * yac_angle_tol)) &&
      // (sq_sin_edge[1] > (yac_angle_tol * yac_angle_tol)) &&
      // (sq_sin_edge[2] > (yac_angle_tol * yac_angle_tol)) &&
      // (sq_sin_edge[3] > (yac_angle_tol * yac_angle_tol))) {

    // int tmp = dotproduct(c, cross[0]) >= 0.0;
    // tmp += dotproduct(d, cross[0]) >= 0.0;
    // tmp += dotproduct(a, cross[1]) >= 0.0;
    // tmp += dotproduct(d, cross[1]) >= 0.0;
    // tmp += dotproduct(a, cross[2]) >= 0.0;
    // tmp += dotproduct(b, cross[2]) >= 0.0;
    // tmp += dotproduct(b, cross[3]) >= 0.0;
    // tmp += dotproduct(c, cross[3]) >= 0.0;
    // // check whether the quad is convex
    // if (!(tmp & 7)) {
      // if (tmp) {
        // for (size_t i = 0; i < num_points; ++i) {
          // if ((dotproduct(point_coords+3*i, cross[0]) + yac_angle_tol >=  0.0) &&
              // (dotproduct(point_coords+3*i, cross[1]) + yac_angle_tol >=  0.0) &&
              // (dotproduct(point_coords+3*i, cross[2]) + yac_angle_tol >=  0.0) &&
              // (dotproduct(point_coords+3*i, cross[3]) + yac_angle_tol >=  0.0))
            // return 1;
        // }
      // } else {
        // for (size_t i = 0; i < num_points; ++i) {
          // if ((dotproduct(point_coords+3*i, cross[0]) - yac_angle_tol <= -0.0) &&
              // (dotproduct(point_coords+3*i, cross[1]) - yac_angle_tol <= -0.0) &&
              // (dotproduct(point_coords+3*i, cross[2]) - yac_angle_tol <= -0.0) &&
              // (dotproduct(point_coords+3*i, cross[3]) - yac_angle_tol <= -0.0))
            // return 1;
        // }
      // }
    // }
  // } else {
    // for (size_t i = 0; i < num_points; ++i)
      // return point_in_cell_generic(point_coords+3*i, cell);
  // }
  // return 0;
// }

// checks whether a point is within a cell
int yac_point_in_cell (double point_coords[3], struct grid_cell cell) {

   struct bounding_circle bnd_circle;
   yac_get_cell_bounding_circle(cell, &bnd_circle);

   return yac_point_in_cell2(point_coords, cell, bnd_circle);
}

// checks whether a point is within a cell
int yac_point_in_cell2 (double point_coords[3], struct grid_cell cell,
                        struct bounding_circle bnd_circle) {

   // check whether the point is within the bounding circle of the cell;

   if (!yac_point_in_bounding_circle_vec(point_coords, &bnd_circle))
      return 1 == 0;

   // check if the point matches any of cell corners
   for (size_t i = 0; i < cell.num_corners; ++i)
     if (points_are_identically(point_coords, cell.coordinates_xyz[i]))
       return 1;

   if ((cell.num_corners == 3) &&
       (cell.edge_type[0] == GREAT_CIRCLE_EDGE) &&
       (cell.edge_type[1] == GREAT_CIRCLE_EDGE) &&
       (cell.edge_type[2] == GREAT_CIRCLE_EDGE))
      return point_in_cell_unstruct_triangle(point_coords, cell);
   else if ((cell.num_corners == 4) &&
       (cell.edge_type[0] == GREAT_CIRCLE_EDGE) &&
       (cell.edge_type[1] == GREAT_CIRCLE_EDGE) &&
       (cell.edge_type[2] == GREAT_CIRCLE_EDGE) &&
       (cell.edge_type[3] == GREAT_CIRCLE_EDGE))
      return point_in_cell_unstruct_quad(point_coords, cell);
   else
      return point_in_cell_generic(point_coords, cell);
}
