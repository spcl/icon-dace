/**
 * @file intersection.c
 * @brief Set of functions to determine the intersection between two edges
 *
 * very interesting literature:
 * - http://geospatialmethods.org/spheres/GCIntersect.html
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
#include <stdio.h>

#include "utils.h"
#include "geometry.h"

// angle tolerance
static double const tol = 1.0e-12;

enum {
  p_on_a = 1,
  q_on_a = 2,
  p_on_b = 4,
  q_on_b = 8,
  circles_are_identical = 16,
};

enum {
  no_intersection = 0,
  a_between_cd = 1,
  b_between_cd = 2,
  c_between_ab = 4,
  d_between_ab = 8,
};

/**
 * computes the intersection of two great circles
 * @param[in]  norm_vector_a norm vector of circle a
 * @param[in]  norm_vector_b norm vector of circle b
 * @param[out] p             first intersection point
 * @param[out] q             second intersection point
 * @return -1 if the two circles are identical; 2 otherwise
 */
static int gcxgc(
  double const norm_vector_a[3], double const norm_vector_b[3],
  double p[3], double q[3]) {

  // compute p and q
  // p' = (a X b) X (c X d)
  // p  = p' / length(p')
  // q = -p
  double temp_p[3];
  crossproduct_kahan(norm_vector_a, norm_vector_b, temp_p);

  double sq_sin_plane_angle =
    temp_p[0] * temp_p[0] + temp_p[1] * temp_p[1] + temp_p[2] * temp_p[2];

  // if both circles are on the same plane
  if (sq_sin_plane_angle <= yac_sq_angle_tol) return -1;

  double scale = 1.0 / sqrt(sq_sin_plane_angle);
  p[0] = (temp_p[0] * scale);
  p[1] = (temp_p[1] * scale);
  p[2] = (temp_p[2] * scale);
  q[0] = -p[0], q[1] = -p[1], q[2] = -p[2];

  return 2;
}

#if defined __INTEL_COMPILER
#pragma intel optimization_level 0
#elif defined _CRAYC
#pragma _CRI noopt
#endif

static int loncxlatc(
  double const gc_norm_vector[3], double const z,
  double p[3], double q[3]) {

  // based on gcxlatc with some optimisations for circles of longitude

  // square of z of lat circle
  double sq_z = z * z;

  double b = sqrt(MAX(1.0 - sq_z,0.0));
  double b_n0 = b * gc_norm_vector[0];
  double b_n1 = b * gc_norm_vector[1];

  p[0] = b_n1;
  p[1] = -b_n0;
  p[2] = z;
  q[0] = -p[0];
  q[1] = -p[1];
  q[2] = z;

  return 2;
}

static int gcxlatc(
  double const gc_norm_vector[3], double const z,
  double p[3], double q[3]) {

/* How to compute intersection point(s) between great circle (given by its
 * norm vector) and a circle of latitude (give by z coordinate):
 *
 * n = norm vector of great circle
 * x = vector to intersection point(s)
 * z = cos(lat)   // where lat is the latitude of the circle of latitude
 *
 * I.   n * x = 0 // intersection is orthogonal to the norm vector
 * II.  x * x = 1 // intersection vector has a length of 1
 * III. x_2   = z // z-coordinate of intersection point is z
 *
 * from I we can derive:
 * x_1=-(n_0*x_0+n_2*z)/n_1
 *
 * x_1 and x_2 in II:
 * x_0^2+(n_0*x_0+n_2*z)^2/n_1^2+z^2=1
 * x_0^2+(2*n_0*n_2*z)/(n_0^2+n_1^2)*x_0+
 *   (n_2^2*z^2+n_1^2*z^2-n_1^2)/(n_0^2+n_1^2)=0
 *
 * with M = n_0^2+n_1^2:
 * x_0^2 + (2*n_0*n_2*z)/M*x_0 + (n_2^2*z^2+n_1^2*z^2-n_1^2)/M = 0
 *
 * x_0 = -(n_0*n_2*z)/M
 *       +- sqrt((n_0*n_2*z)^2/M^2 - (n_2^2*z^2+n_1^2*z^2-n_1^2)/M)
 *     = -(*n_2*z/M)*n_0 +- sqrt(M-z^2)/M*n_1
 *
 * with a = -n_2*z/M and b = sqrt(M-z^2)/M:
 * x_0 = a * n_0 +- b * n_1
 * x_1 = a * n_1 -+ b * n_0
 */

  // if the great circle is a circle of longitude
  if (fabsl(gc_norm_vector[2]) < yac_angle_low_tol)
     return loncxlatc(gc_norm_vector, z, p, q);

  double sq_n[3] = {gc_norm_vector[0] * gc_norm_vector[0],
                    gc_norm_vector[1] * gc_norm_vector[1],
                    gc_norm_vector[2] * gc_norm_vector[2]};

  // square of highest z-value of great circle
  double max_sq_gc_z = sq_n[0] + sq_n[1];

  // square of z of lat circle
  double sq_z = z * z;

  // if both circles are the equator -> identical circles
  if ((max_sq_gc_z < yac_sq_angle_tol) && (sq_z < yac_sq_angle_tol)) return -1;

  // determine number of intersections
  int num_intersections =
    ((max_sq_gc_z > (sq_z - tol)) - (max_sq_gc_z < (sq_z + tol))) + 1;

  switch (num_intersections) {

    // if no intersection is possible
    case (0): break;

    // the great circles touches the circle of latitude in a single point
    // ==> (b = 0)
    case (1): {

      double a = -(gc_norm_vector[2] * z) / max_sq_gc_z;

      p[0] = ((q[0] = gc_norm_vector[0] * a));
      p[1] = ((q[1] = gc_norm_vector[1] * a));
      p[2] = z, q[2] = z;
      break;
    }

    default:
    case (2): {

      double a = -(gc_norm_vector[2] * z) / max_sq_gc_z;
      double b = sqrt(max_sq_gc_z - sq_z)/max_sq_gc_z;
      double a_n0 = a * gc_norm_vector[0];
      double a_n1 = a * gc_norm_vector[1];
      double b_n0 = b * gc_norm_vector[0];
      double b_n1 = b * gc_norm_vector[1];

      p[0] = a_n0 + b_n1;
      p[1] = a_n1 - b_n0;
      p[2] = z;
      q[0] = a_n0 - b_n1;
      q[1] = a_n1 + b_n0;
      q[2] = z;
      break;
    }
  }

  return num_intersections;
}

#if defined _CRAYC
#pragma _CRI opt
#endif

static int loncxlonc(
  double const norm_vector_a[3], double const norm_vector_b[3],
  double p[3], double q[3]) {

  double sq_plane_diff =
    MIN((norm_vector_a[0] - norm_vector_b[0]) *
        (norm_vector_a[0] - norm_vector_b[0]) +
        (norm_vector_a[1] - norm_vector_b[1]) *
        (norm_vector_a[1] - norm_vector_b[1]),
        (norm_vector_a[0] + norm_vector_b[0]) *
        (norm_vector_a[0] + norm_vector_b[0]) +
        (norm_vector_a[1] + norm_vector_b[1]) *
        (norm_vector_a[1] + norm_vector_b[1]));

  // if both circles are on the same circle of longitude
  if (sq_plane_diff < yac_sq_angle_tol) return -1;

  // the intersection points of the two circles of longitude are the poles
  p[0] = 0, p[1] = 0; p[2] = 1;
  q[0] = 0, q[1] = 0; q[2] = -1;

  return 2;
}

/** \example test_pxgc.c
 * This example contains some tests for the pxgc routine
 */

static int pxgc(
  double const vec[3], double const norm_vector[3], double p[3], double q[3]) {

  // compute the angle between the point vector and the norm vector
  // we use kahan summation to increase precision
  double dot = vec[0] * norm_vector[0];
  double err = fma(vec[0], norm_vector[0], -dot);
  dot = fma(vec[1], norm_vector[1], dot);
  err += fma(vec[1], norm_vector[1], -vec[1]*norm_vector[1]);
  dot = fma(vec[2], norm_vector[2], dot);
  err += fma(vec[2], norm_vector[2], -vec[2]*norm_vector[2]);
  dot += err;

  // if the point is not on the plane
  if (fabsl(dot) > yac_angle_tol) return 0;

  q[0] = -((p[0] = vec[0]));
  q[1] = -((p[1] = vec[1]));
  q[2] = -((p[2] = vec[2]));
  return 1;
}

static int pxlatc(
  double const vec[3], double z, double p[3], double q[3]) {

  // if the point is not on the plane of the circle of latitude
  if (fabs(vec[2] - z) > tol) return 0;

  q[0] = -((p[0] = vec[0]));
  q[1] = -((p[1] = vec[1]));
  q[2] =  ((p[2] = vec[2]));
  return 1;
}

static int pxp(
  double const vec_a[3], double const vec_b[3], double p[3], double q[3]) {

  // if the points are identical
  if (points_are_identically(vec_a, vec_b)) {
    q[0] = -((p[0] = vec_a[0]));
    q[1] = -((p[1] = vec_a[1]));
    q[2] = -((p[2] = vec_a[2]));
    return 1;
  } else {
    return 0;
  }
}

int yac_circle_intersect(
  struct yac_circle a, struct yac_circle b, double p[3], double q[3]) {

  YAC_ASSERT(
    ((a.type == GREAT_CIRCLE) ||
     (a.type == LON_CIRCLE) ||
     (a.type == LAT_CIRCLE) ||
     (a.type == POINT)) &&
    ((b.type == LON_CIRCLE) ||
     (b.type == GREAT_CIRCLE) ||
     (b.type == LAT_CIRCLE) ||
     (b.type == POINT)),
    "ERROR(yac_circle_intersect): unknown circle type.")

  if (((a.type == GREAT_CIRCLE) && (b.type == GREAT_CIRCLE)) ||
      ((a.type == GREAT_CIRCLE) && (b.type == LON_CIRCLE)) ||
      ((a.type == LON_CIRCLE) && (b.type == GREAT_CIRCLE))) {
    return gcxgc(a.data.gc.norm_vector, b.data.gc.norm_vector, p, q);
  } else if (((a.type == GREAT_CIRCLE) && (b.type == LAT_CIRCLE))) {
    return gcxlatc(a.data.gc.norm_vector, b.data.lat.z, p, q);
  } else if (((a.type == LAT_CIRCLE) && (b.type == GREAT_CIRCLE))) {
    return gcxlatc(b.data.gc.norm_vector, a.data.lat.z, p, q);
  } else if ((a.type == LON_CIRCLE) && (b.type == LON_CIRCLE)) {
    return loncxlonc(a.data.lon.norm_vector, b.data.lon.norm_vector, p, q);
  } else if ((a.type == LON_CIRCLE) && (b.type == LAT_CIRCLE)) {
    return loncxlatc(a.data.lon.norm_vector, b.data.lat.z, p, q);
  } else if ((a.type == LAT_CIRCLE) && (b.type == LON_CIRCLE)) {
    return loncxlatc(b.data.lon.norm_vector, a.data.lat.z, p, q);
  } else if ((a.type == POINT) && (b.type == GREAT_CIRCLE)) {
    return pxgc(a.data.p.vec, b.data.gc.norm_vector, p, q);
  } else if ((a.type == POINT) && (b.type == LON_CIRCLE)) {
    return pxgc(a.data.p.vec, b.data.gc.norm_vector, p, q);
  } else if ((a.type == POINT) && (b.type == LAT_CIRCLE)) {
    return pxlatc(a.data.p.vec, b.data.lat.z, p, q);
  } else if ((a.type == POINT) && (b.type == POINT)) {
    return pxp(a.data.p.vec, b.data.p.vec, p, q);
  } else if ((a.type == GREAT_CIRCLE) && (b.type == POINT)) {
    return pxgc(b.data.p.vec, a.data.gc.norm_vector, p, q);
  } else if ((a.type == LON_CIRCLE) && (b.type == POINT)) {
    return pxgc(b.data.p.vec, a.data.gc.norm_vector, p, q);
  } else /*if ((a.type == LAT_CIRCLE) && (b.type == POINT))*/ {
    return pxlatc(b.data.p.vec, a.data.lat.z, p, q);
  }
}

/**
 * determines whether vector p is between vectors a and b
 * (it is assumed, that a, b, and p are in the same plane)
 */
static inline int vector_is_between (
  double const a[], double const b[], double const p[],
  double sq_len_diff_ab) {

  // In case the angle between the vectors a and b is 180 degree, the angle
  // between the vectors pa and pb is exactly 90 degree (Pythagorean theorem).
  // In all other cases the angle between pa and pb must be bigger than 90
  // degree for p to be between a and b
  // from this we can deduce:
  // ||ab||^2 >= ||ap||^2 + ||bp||^2 => (Pythagorean theorem)

  // the tol is required in case p is very close to a or b
  double sq_len_diff_vec_ap = sq_len_diff_vec(a,p);
  double sq_len_diff_vec_bp = sq_len_diff_vec(b,p);
  return sq_len_diff_ab + tol >= sq_len_diff_vec_ap + sq_len_diff_vec_bp;
}

static int vector_is_between_lat (
  double const a[], double const b[], double const p[]) {

/* determines whether p is between a and b
   (a, b, p have the same latitude)*/

/*  I. a_0 * alpha + b_0 * beta = p_0
   II. a_1 * alpha + b_1 * beta = p_1

   if alpha > 0 and beta > 0 -> p is between a and b */

   if (fabs(fabs(a[2]) - 1.0) < tol) return 1;

   double fabs_a[2] = {fabs(a[0]), fabs(a[1])};
   int flag = fabs_a[0] > fabs_a[1];
   double a_0 = a[flag], a_1 = a[flag^1];
   double b_0 = b[flag], b_1 = b[flag^1];
   double p_0 = p[flag], p_1 = p[flag^1];

   double temp = b_0 - (b_1 * a_0) / a_1;

   YAC_ASSERT(
      (fabs(temp) > tol),
      "ERROR(vector_is_between_lat): "
      "routine does not support zero length edges")
   /*// if a and b are nearly identical
   if (fabs(temp) < tol)
     return (fabs(a_0 - p_0) < tol) && (fabs(a_1 - p_1) < tol);*/

   double beta = (p_0 - (p_1 * a_0) / a_1) / temp;

   if (beta < -tol) return 0;

   double alpha = (p_1 - b_1 * beta) / a_1;

   return alpha > -tol;
}

static void compute_edge_middle_point_vec(
  enum yac_edge_type edge_type, double const a[3], double const b[3],
  double middle[3]) {

  YAC_ASSERT(
    (edge_type == GREAT_CIRCLE_EDGE) ||
    (edge_type == LAT_CIRCLE_EDGE) ||
    (edge_type == LON_CIRCLE_EDGE),
    "ERROR(compute_edge_middle_point_vec): invalide edge_type")

  if (edge_type == LAT_CIRCLE_EDGE) {
    middle[0] = a[0]+b[0], middle[1] = a[1]+b[1], middle[2] = a[2];
    double temp = middle[0]*middle[0] + middle[1]*middle[1];
    double scale = sqrt((1.0 - middle[2] * middle[2])/temp);
    middle[0] *= scale;
    middle[1] *= scale;
  } else {
    middle[0] = a[0]+b[0], middle[1] = a[1]+b[1], middle[2] = a[2]+b[2];
    normalise_vector(middle);
  }
}

static int yac_identical_circles_vec(
  enum yac_edge_type edge_type, double const a[3], double const b[3],
  double const c[3], double const d[3], double p[3], double q[3],
  int intersection_type) {

  YAC_ASSERT(
    (intersection_type == no_intersection) ||
    (intersection_type == a_between_cd + b_between_cd) ||
    (intersection_type == a_between_cd + b_between_cd + c_between_ab) ||
    (intersection_type == a_between_cd + b_between_cd + d_between_ab) ||
    (intersection_type == a_between_cd + c_between_ab) ||
    (intersection_type == a_between_cd + c_between_ab + d_between_ab) ||
    (intersection_type == a_between_cd + d_between_ab) ||
    (intersection_type == b_between_cd + c_between_ab) ||
    (intersection_type == b_between_cd + c_between_ab + d_between_ab) ||
    (intersection_type == b_between_cd + d_between_ab) ||
    (intersection_type == c_between_ab + d_between_ab) ||
    (intersection_type == a_between_cd + b_between_cd +
                          c_between_ab + d_between_ab),
    "internal error")

  switch (intersection_type) {
    default:

    // edges do not intersect
    case (no_intersection):
      // use middle points of both edges as intersection points
      compute_edge_middle_point_vec(edge_type, a, b, p);
      compute_edge_middle_point_vec(edge_type, c, d, q);
      return p_on_a + q_on_b + circles_are_identical;

    // both vertices of first edge overlap with the second edge or
    // edges are identical
    case (a_between_cd + b_between_cd + c_between_ab + d_between_ab):
    case (a_between_cd + b_between_cd):
    case (a_between_cd + b_between_cd + c_between_ab):
    case (a_between_cd + b_between_cd + d_between_ab):
      p[0] = a[0], p[1] = a[1], p[2] = a[2];
      q[0] = b[0], q[1] = b[1], q[2] = b[2];
      break;

    // both edges of the second edge overlap with the first edge
    case (c_between_ab + d_between_ab):
    case (a_between_cd + c_between_ab + d_between_ab):
    case (b_between_cd + c_between_ab + d_between_ab):
      p[0] = c[0], p[1] = c[1], p[2] = c[2];
      q[0] = d[0], q[1] = d[1], q[2] = d[2];
      break;

    // only one vertex of each edge overlaps with the other edge
    case (a_between_cd + c_between_ab):
      p[0] = a[0], p[1] = a[1], p[2] = a[2];
      q[0] = c[0], q[1] = c[1], q[2] = c[2];
      break;
    case (a_between_cd + d_between_ab):
      p[0] = a[0], p[1] = a[1], p[2] = a[2];
      q[0] = d[0], q[1] = d[1], q[2] = d[2];
      break;
    case (b_between_cd + c_between_ab):
      p[0] = b[0], p[1] = b[1], p[2] = b[2];
      q[0] = c[0], q[1] = c[1], q[2] = c[2];
      break;
    case (b_between_cd + d_between_ab):
      p[0] = b[0], p[1] = b[1], p[2] = b[2];
      q[0] = d[0], q[1] = d[1], q[2] = d[2];
      break;
  }
  return p_on_a + q_on_a + p_on_b + q_on_b + circles_are_identical;
}

// computes intersection of great circle edges that are on
// the same great circle
static int yac_identical_gcxgc_vec(
  double const a[3], double const b[3], double const c[3], double const d[3],
  double p[3], double q[3]) {

  double sq_len_diff_ab = sq_len_diff_vec(a, b);
  double sq_len_diff_cd = sq_len_diff_vec(c, d);

  return
    yac_identical_circles_vec(
      GREAT_CIRCLE_EDGE, a, b, c, d, p, q,
      (vector_is_between(c, d, a, sq_len_diff_cd) << 0) +
      (vector_is_between(c, d, b, sq_len_diff_cd) << 1) +
      (vector_is_between(a, b, c, sq_len_diff_ab) << 2) +
      (vector_is_between(a, b, d, sq_len_diff_ab) << 3));
}

// determines the return value for the provided intersection points p and q
static int yac_check_pq_gcxgc(
  double const a[3], double const b[3], double const c[3], double const d[3],
  double const p[3], double const q[3]) {

  // square of the Euclidean distance between the vertices of the two edges
  double sq_len_diff_ab = sq_len_diff_vec(a, b);
  double sq_len_diff_cd = sq_len_diff_vec(c, d);

  return (vector_is_between(a, b, p, sq_len_diff_ab) << 0) +
         (vector_is_between(a, b, q, sq_len_diff_ab) << 1) +
         (vector_is_between(c, d, p, sq_len_diff_cd) << 2) +
         (vector_is_between(c, d, q, sq_len_diff_cd) << 3);
}

static void compute_norm_vector(
  double const a[3], double const b[3], double norm_vector[3]) {

  norm_vector[0] = a[1] * b[2] - a[2] * b[1];
  norm_vector[1] = a[2] * b[0] - a[0] * b[2];
  norm_vector[2] = a[0] * b[1] - a[1] * b[0];

  double scale = 1.0 / sqrt(norm_vector[0] * norm_vector[0] +
                            norm_vector[1] * norm_vector[1] +
                            norm_vector[2] * norm_vector[2]);
  norm_vector[0] *= scale;
  norm_vector[1] *= scale;
  norm_vector[2] *= scale;
}

/** \example test_gcxgc.c
 * This contains examples on how to use \ref yac_gcxgc_vec
 */
/**
 * computes the intersection points of two great circles \n
 * based on http://www.geoclub.de/viewtopic.php?f=54&t=29689
 * @param[in]  a first point of edge a the is on a great circle
 * @param[in]  b second point of edge a the is on a great circle
 * @param[in]  c first point of edge b the is on a great circle
 * @param[in]  d second point of edge b the is on a great circle
 * @param[out] p intersection point
 * @param[out] q intersection point
 * @return  0 if the intersection points are neither on edge a or b \n
 *          -1 if an error occurred \n
 *          1st bit will be set if p is on edge a \n
 *          2nd bit will be set if q is on edge a \n
 *          3rd bit will be set if p is on edge b \n
 *          4th bit will be set if q is on edge b
 *          5th bit will be set if both great circles are identically
 * @remarks both edges need to have length of at least yac_angle_tol
 */
static int yac_gcxgc_vec(
  double const a[3], double const b[3], double const c[3], double const d[3],
  double p[3], double q[3]) {

  double norm_vector_ab[3], norm_vector_cd[3];

  compute_norm_vector(a, b, norm_vector_ab);
  compute_norm_vector(c, d, norm_vector_cd);

  switch (gcxgc(norm_vector_ab, norm_vector_cd, p, q)) {
    default:
    case(2):
      return yac_check_pq_gcxgc(a, b, c, d, p, q);
    case(-1):
      return yac_identical_gcxgc_vec(a, b, c, d, p, q);
  };
}

/** \example test_latcxlatc.c
 * This contains examples on \ref yac_latcxlatc_vec
 */
/** \brief compute the intersection point two circles of latitude
 *
 * compute the intersection points of two circle of latitude
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 *      - 5th bit will be set if both edges are on the same circle of latitude
 * @remarks both edges need to have length of at least yac_angle_tol
 **/
static int yac_latcxlatc_vec(
  double const a[3], double const b[3], double const c[3], double const d[3],
  double p[3], double q[3]) {

  // two circles of latitude can only intersect if they are on the same latitude
  if (fabs(a[2] - c[2]) > tol) return -1;
  else
    return
      yac_identical_circles_vec(
        LAT_CIRCLE_EDGE, a, b, c, d, p, q,
        (vector_is_between_lat(c, d, a) << 0) +
        (vector_is_between_lat(c, d, b) << 1) +
        (vector_is_between_lat(a, b, c) << 2) +
        (vector_is_between_lat(a, b, d) << 3));
}

static void compute_lon_norm_vector(
  double const a[3], double const b[3], double norm_vector[3]) {

  norm_vector[0] = a[1] * b[2] - a[2] * b[1];
  norm_vector[1] = a[2] * b[0] - a[0] * b[2];
  norm_vector[2] = 0.0;

  double scale = 1.0 / sqrt(norm_vector[0] * norm_vector[0] +
                            norm_vector[1] * norm_vector[1]);
  norm_vector[0] *= scale;
  norm_vector[1] *= scale;
}

/** \example test_loncxlonc.c
 * This contains examples on \ref yac_loncxlonc_vec
 */
/** \brief compute the intersection point two circles of longitude
 *
 * compute the intersection points of two circle of longitude
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 *      - 5th bit will be set if both edges are on the same circle of longitude
 * @remarks both edges need to have length of at least yac_angle_tol
 **/
static int yac_loncxlonc_vec(
  double const a[3], double const b[3], double const c[3], double const d[3],
  double p[3], double q[3]) {

  double norm_vector_ab[3], norm_vector_cd[3];
  compute_lon_norm_vector(a, b, norm_vector_ab);
  compute_lon_norm_vector(c, d, norm_vector_cd);

  switch(loncxlonc(norm_vector_ab, norm_vector_cd, p, q)) {

    default:
    case(2):
      return yac_check_pq_gcxgc(a, b, c, d, p, q);
    case(-1):
      return yac_identical_gcxgc_vec(a, b, c, d, p, q);
  };
}

// determines the return value for the provided intersection points p and q
static int yac_check_pq_gcxlatc(
  double const a[3], double const b[3], double const c[3], double const d[3],
  double const p[3], double const q[3]) {

  // square of the Euclidean distance between the vertices of the two edges
  double sq_len_diff_ab = sq_len_diff_vec(a, b);

  return
    (vector_is_between(a, b, p, sq_len_diff_ab) << 0) +
    (vector_is_between(a, b, q, sq_len_diff_ab) << 1) +
    (vector_is_between_lat(c, d, p)             << 2) +
    (vector_is_between_lat(c, d, q)             << 3);
}

/** \example test_loncxlatc.c
 * This contains examples on \ref yac_loncxlatc_vec
 */
/** \brief compute the intersection point of a meridian and a parallel
 *
 * compute the intersection points of a circle of longitude (defined by a and b)
 * and a circle of latitude (defined by c and d)
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *      - 0 if the intersection points are neither between (a and b) or (c and d)
 *      - -1 if an error occurred
 *      - 1st bit will be set if p is between a and b
 *      - 2nd bit will be set if q is between a and b
 *      - 3rd bit will be set if p is between c and d
 *      - 4th bit will be set if q is between c and d
 * @remarks both edges need to have length of at least yac_angle_tol
 **/
static int yac_loncxlatc_vec (
  double const a[3], double const b[3], double const c[3], double const d[3],
  double p[3], double q[3]) {

  double norm_vector_ab[3];
  compute_lon_norm_vector(a, b, norm_vector_ab);

  int num_intersections = loncxlatc(norm_vector_ab, c[2], p, q);
  YAC_ASSERT(num_intersections == 2, "ERROR(yac_loncxlatc_vec): internal error")

  return yac_check_pq_gcxlatc(a, b, c, d, p, q);
}

/** \example test_gcxlatc.c
 * This contains examples on gcxlatc_vec.
 */
/** \brief compute the intersection of a great circle with the parallel
 *
 * compute the intersection points of a great circle (defined by a and b)
 * and a circle of latitude (defined by c and d)
 * (both circles need to have a length of at least yac_angle_tol)
 * if p and q are != NULL they contain the intersection points
 * the return value is:
 *    - 0 if the intersection points are neither between (a and b) or (c and d)
 *    - -1 if the two circles do not intersect or an error occurred
 *    - 1st bit will be set if p is between a and b
 *    - 2nd bit will be set if q is between a and b
 *    - 3rd bit will be set if p is between c and d
 *    - 4th bit will be set if q is between c and d
 *    - 5th bit will be set if both circles are identically
 * \remarks if -1 is returned neither p or q is set
 * \remarks if the two circles only have one intersection point,
 *          p and q will be identically, but only the p bits will be set
 * @remarks both edges need to have length of at least yac_angle_tol
 **/

#if defined __INTEL_COMPILER
#pragma intel optimization_level 0
#elif defined _CRAYC
#pragma _CRI noopt
#endif

static int yac_gcxlatc_vec(
  double const a[3], double const b[3], double const c[3], double const d[3],
  double p[3], double q[3]) {

  double norm_vector_ab[3];
  compute_norm_vector(a, b, norm_vector_ab);

  int num_intersections = gcxlatc(norm_vector_ab, c[2], p, q);

  switch (num_intersections) {
    default:
    case(2):
    case(1):
      return yac_check_pq_gcxlatc(a, b, c, d, p, q);
    case(0):
      return -1;
    case(-1):
      return yac_latcxlatc_vec(a, b, c, d, p, q);
  }
}

#if defined _CRAYC
#pragma _CRI opt
#endif

static int yac_pxp_vec(
  double const a[3], double const b[3], double p[3], double q[3]) {

  return (pxp(a, b, p, q))?(p_on_a + p_on_b):-1;
}

static int yac_pxgc_vec(
  double const point[3], double const a[3], double const b[3],
  double p[3], double q[3]) {

  // compute norm vector for the plane of ab
  double norm_ab[3];
  compute_norm_vector(a, b, norm_ab);

  // if the point is on the plane of ab
  if (pxgc(point, norm_ab, p, q)) {

    int ret_value = p_on_a;
    double sq_len_diff_ab = sq_len_diff_vec(a, b);

    if      (vector_is_between(a, b, p, sq_len_diff_ab)) ret_value |= p_on_b;
    else if (vector_is_between(a, b, q, sq_len_diff_ab)) ret_value |= q_on_b;

    return ret_value;

  } else {
    return -1;
  }
}

static int yac_pxlat_vec(
  double const point[3], double const a[3], double const b[3],
  double p[3], double q[3]) {

  if (pxlatc(point, a[2], p, q)) {

    int ret_value = p_on_a;

    if      (vector_is_between_lat(a, b, p)) ret_value |= p_on_b;
    else if (vector_is_between_lat(a, b, q)) ret_value |= q_on_b;

    return ret_value;
  } else {
    return -1;
  }
}

static inline int adjust_ret_value(int ret_value) {
  return (ret_value & (~(1 + 2 + 4 + 8))) +
         ((ret_value & (1 + 2)) << 2) +
         ((ret_value & (4 + 8)) >> 2);
}

int yac_intersect_vec (
  enum yac_edge_type edge_type_a, double const a[3], double const b[3],
  enum yac_edge_type edge_type_b, double const c[3], double const d[3],
  double p[3], double q[3]) {

  YAC_ASSERT(
    ((edge_type_a == LON_CIRCLE_EDGE) ||
     (edge_type_a == LAT_CIRCLE_EDGE) ||
     (edge_type_a == GREAT_CIRCLE_EDGE)) &&
    ((edge_type_b == LON_CIRCLE_EDGE) ||
     (edge_type_b == LAT_CIRCLE_EDGE) ||
     (edge_type_b == GREAT_CIRCLE_EDGE)), "ERROR: unknown edge type.")

  int edge_a_is_point = points_are_identically(a, b);
  int edge_b_is_point = points_are_identically(c, d);

  // if both edges are points
  if (edge_a_is_point && edge_b_is_point) {

    return yac_pxp_vec(a, c, p, q);

  // if one edge is a point
  } else if (edge_a_is_point || edge_b_is_point) {

    enum yac_edge_type edge_type[2] = {edge_type_a, edge_type_b};
    double const * edge[2][2] = {{a,b}, {c,d}};
    double const * point[2] = {c,a};

    int ret_value;
    switch (edge_type[edge_a_is_point]) {
      case (LAT_CIRCLE_EDGE):
        ret_value =
          yac_pxlat_vec(
            point[edge_a_is_point], 
            edge[edge_a_is_point][0], edge[edge_a_is_point][1], p, q);
        break;
      default:
      case (LON_CIRCLE_EDGE):
      case (GREAT_CIRCLE_EDGE):
        ret_value =
          yac_pxgc_vec(
            point[edge_a_is_point], 
            edge[edge_a_is_point][0], edge[edge_a_is_point][1], p, q);
        break;
    }

    if (edge_a_is_point) return ret_value;
    else                 return adjust_ret_value(ret_value);

   // if both edges are on circles of latitude
  } else if (edge_type_a == LAT_CIRCLE_EDGE &&
             edge_type_b == LAT_CIRCLE_EDGE) {

    return yac_latcxlatc_vec(a, b, c, d, p, q);

  // if both edges are on circle of longitude
  } else if (edge_type_a == LON_CIRCLE_EDGE &&
             edge_type_b == LON_CIRCLE_EDGE) {

    return yac_loncxlonc_vec(a, b, c, d, p, q);

  // if both edges are on great circles
  } else if ((edge_type_a == GREAT_CIRCLE_EDGE &&
              edge_type_b == GREAT_CIRCLE_EDGE) ||
             (edge_type_a == LON_CIRCLE_EDGE   &&
              edge_type_b == GREAT_CIRCLE_EDGE) ||
             (edge_type_a == GREAT_CIRCLE_EDGE &&
              edge_type_b == LON_CIRCLE_EDGE)) {

    return yac_gcxgc_vec(a, b, c, d, p, q);

  // if one edge a is on a great circle and edge b on a circle of latitude
  } else if (edge_type_a == GREAT_CIRCLE_EDGE &&
             edge_type_b == LAT_CIRCLE_EDGE) {

    return yac_gcxlatc_vec(a, b, c, d, p, q);

  // if one edge a is on a circle of latitude and edge b on a great circle
  } else if (edge_type_a == LAT_CIRCLE_EDGE &&
             edge_type_b == GREAT_CIRCLE_EDGE ) {

    return adjust_ret_value(yac_gcxlatc_vec(c, d, a, b, p, q));

  // if one edge a is on a circle of longitude and edge b on a circle of latitude
  } else if (edge_type_a == LON_CIRCLE_EDGE &&
             edge_type_b == LAT_CIRCLE_EDGE) {

    return yac_loncxlatc_vec(a, b, c, d, p, q);

  // if one edge a is on a circle of latitude and edge b on a circle of longitude
  } else /* if (edge_type_a == LAT_CIRCLE_EDGE &&
             edge_type_b == LON_CIRCLE_EDGE )*/ {

    return adjust_ret_value(yac_loncxlatc_vec(c, d, a, b, p, q));
  }
}

int yac_point_on_edge(
  double p[3], double const a[3], double const b[3],
  enum yac_circle_type circle_type) {

  YAC_ASSERT(
    (circle_type == LON_CIRCLE) ||
    (circle_type == LAT_CIRCLE) ||
    (circle_type == GREAT_CIRCLE) ||
    (circle_type == POINT),
    "ERROR(yac_point_on_edge): unknown edge type.")

  switch (circle_type) {
    default:
    case (GREAT_CIRCLE):
    case (LON_CIRCLE):
      return vector_is_between(a, b, p, sq_len_diff_vec(a, b));
    case (LAT_CIRCLE):
      return vector_is_between_lat(a, b, p);
    case(POINT):
      return points_are_identically(p, a);
  }
}
