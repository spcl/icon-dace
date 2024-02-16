/**
 * @file bnd_circle.c
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

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "geometry.h"
#include "utils.h"
#include "grid.h"
#include "ensure_array_size.h"

/** \file bnd_circle.c
 *  \brief Set of functions to calculate a bounding circle around a certain set of points.
 *
 **/

// computes the circumscribe circle of a quad on the sphere
// it is assumed that the edges are circles of longitude and latitude
void yac_get_cell_circumscribe_circle_reg_quad(
   double a[3], double b[3], double c[3],
   struct bounding_circle * bnd_circle) {

   yac_get_cell_circumscribe_circle_unstruct_triangle(a, b, c, bnd_circle);
}

// computes the bounding circle of a quad on the sphere
// it is assumed that the edges are circles of longitude and latitude
void yac_get_cell_bounding_circle_reg_quad(
   double a[3], double b[3], double c[3],
   struct bounding_circle * bnd_circle) {

   yac_get_cell_circumscribe_circle_unstruct_triangle(a, b, c, bnd_circle);
}

// computes the circumscribe circle of a triangle on the sphere
// it is assumed that all edges are great circles
void yac_get_cell_circumscribe_circle_unstruct_triangle(
   double a[3], double b[3], double c[3],
   struct bounding_circle * bnd_circle) {

   double ab[3] = {a[0]-b[0], a[1]-b[1], a[2]-b[2]},
          ac[3] = {b[0]-c[0], b[1]-c[1], b[2]-c[2]};

   // it is assumed that the angles of a triangle do not get too small...
   // crossproduct_kahan(ab, ac, bnd_circle->base_vector);
   crossproduct_d(ab, ac, bnd_circle->base_vector);
   normalise_vector(bnd_circle->base_vector);

   // find biggest component of base_vector
   double fabs_base_vector[3] = {fabs(bnd_circle->base_vector[0]),
                                 fabs(bnd_circle->base_vector[1]),
                                 fabs(bnd_circle->base_vector[2])};
   int biggest_component_index =
    ((fabs_base_vector[0] < fabs_base_vector[1]) |
     (fabs_base_vector[0] < fabs_base_vector[2])) <<
    (fabs_base_vector[1] < fabs_base_vector[2]);

   if ((bnd_circle->base_vector[biggest_component_index] > 0) ^
       (a[biggest_component_index] > 0)) {
      bnd_circle->base_vector[0] = -bnd_circle->base_vector[0];
      bnd_circle->base_vector[1] = -bnd_circle->base_vector[1];
      bnd_circle->base_vector[2] = -bnd_circle->base_vector[2];
   }

   bnd_circle->inc_angle =
       sum_angles_no_check(
          get_vector_angle_2(
             bnd_circle->base_vector, a), SIN_COS_TOL);
   bnd_circle->sq_crd = DBL_MAX;
}

static inline double get_sin_vector_angle(
  double a[3], double b[3]) {

  double cross_ab[3];
  crossproduct_kahan(a, b, cross_ab);

  return sqrt(cross_ab[0]*cross_ab[0] +
              cross_ab[1]*cross_ab[1] +
              cross_ab[2]*cross_ab[2]);
}

// computes the bounding circle of a triangle on the sphere
// it is assumed that all edges are great circles
void yac_get_cell_bounding_circle_unstruct_triangle(
  double a[3], double b[3], double c[3],
  struct bounding_circle * bnd_circle) {

  double middle_point[3];

  middle_point[0] = a[0] + b[0] + c[0];
  middle_point[1] = a[1] + b[1] + c[1];
  middle_point[2] = a[2] + b[2] + c[2];

  normalise_vector(middle_point);

  double cos_angles[3] = {middle_point[0] * a[0] +
                          middle_point[1] * a[1] +
                          middle_point[2] * a[2],
                          middle_point[0] * b[0] +
                          middle_point[1] * b[1] +
                          middle_point[2] * b[2],
                          middle_point[0] * c[0] +
                          middle_point[1] * c[1] +
                          middle_point[2] * c[2]};

  struct sin_cos_angle inc_angle;

  // find the biggest angle

  if (cos_angles[0] < cos_angles[1]) {
    if (cos_angles[0] < cos_angles[2]) {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, a), cos_angles[0]);
    } else {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, c), cos_angles[2]);
    }
  } else {
    if (cos_angles[1] < cos_angles[2]) {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, b), cos_angles[1]);
    } else {
      inc_angle =
        sin_cos_angle_new(get_sin_vector_angle(middle_point, c), cos_angles[2]);
    }
  }

  bnd_circle->base_vector[0] = middle_point[0];
  bnd_circle->base_vector[1] = middle_point[1];
  bnd_circle->base_vector[2] = middle_point[2];
  bnd_circle->inc_angle = sum_angles_no_check(inc_angle, SIN_COS_TOL);
  bnd_circle->sq_crd = DBL_MAX;
}

// computes a) the angle between the middle point of the edge with the middle
//             point of the polygon
//          b) half of the angle between between the two vertices of the edge
// returns the sum of both angles
static inline struct sin_cos_angle compute_edge_inc_angle(
  double * restrict a, double * restrict b, double * restrict middle_point) {

  double edge_middle_point[3] = {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
  normalise_vector(edge_middle_point);

  struct sin_cos_angle t1 = get_vector_angle_2(edge_middle_point, a);
  struct sin_cos_angle t2 = get_vector_angle_2(edge_middle_point, middle_point);

  return sum_angles_no_check(t1, t2);
}

void yac_get_cell_bounding_circle(struct grid_cell cell,
                                  struct bounding_circle * bnd_circle) {

   double middle_point[3];

   middle_point[0] = 0.0;
   middle_point[1] = 0.0;
   middle_point[2] = 0.0;

   size_t num_corners = cell.num_corners;

   // compute the coordinates in rad and 3d
   for (size_t i = 0; i < num_corners; ++i) {

      middle_point[0] += cell.coordinates_xyz[i][0];
      middle_point[1] += cell.coordinates_xyz[i][1];
      middle_point[2] += cell.coordinates_xyz[i][2];
   }

   normalise_vector(middle_point);

   // compute the angle required for the bounding circle
   struct sin_cos_angle max_angle =
      compute_edge_inc_angle(
         cell.coordinates_xyz[num_corners-1],
         cell.coordinates_xyz[0], middle_point);

   for (size_t i = 0; i < num_corners-1; ++i) {

      struct sin_cos_angle curr_angle =
         compute_edge_inc_angle(
            cell.coordinates_xyz[i], cell.coordinates_xyz[i+1], middle_point);

      if (compare_angles(max_angle, curr_angle) < 0) max_angle = curr_angle;
   }

   bnd_circle->base_vector[0] = middle_point[0];
   bnd_circle->base_vector[1] = middle_point[1];
   bnd_circle->base_vector[2] = middle_point[2];

   bnd_circle->inc_angle = sum_angles_no_check(max_angle, SIN_COS_TOL);
   bnd_circle->sq_crd = DBL_MAX;
}

int yac_extents_overlap(struct bounding_circle * extent_a,
                        struct bounding_circle * extent_b) {

  struct sin_cos_angle base_vector_angle =
    get_vector_angle_2(extent_a->base_vector, extent_b->base_vector);

  struct sin_cos_angle tmp_angle, inc_angle_sum;
  int big_sum =
    sum_angles(extent_a->inc_angle, extent_b->inc_angle, &tmp_angle);

  if (big_sum ||
      sum_angles(tmp_angle, SIN_COS_TOL, &inc_angle_sum))
    return 1;

  return compare_angles(base_vector_angle, inc_angle_sum) <= 0;
}
