/**
 * @file bnd_triangle.c
 * @brief implementation of algorithm for the computation of a bounding triangle
 *
 * @copyright Copyright  (C)  2018 Moritz Hanke <hanke@dkrz.de>
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
#include <string.h>

#include "utils.h"
#include "geometry.h"
#include "yac_lapack_interface.h"

static inline void compute_middle_point(
  double (* restrict vertices)[3], size_t num_vertices,
  double * restrict middle_point) {

  middle_point[0] = 0;
  middle_point[1] = 0;
  middle_point[2] = 0;

  for (size_t i = 0; i < num_vertices; ++i) {

    middle_point[0] += vertices[i][0];
    middle_point[1] += vertices[i][1];
    middle_point[2] += vertices[i][2];
  }

  double scale = 1.0 / sqrt(middle_point[0] * middle_point[0] +
                            middle_point[1] * middle_point[1] +
                            middle_point[2] * middle_point[2]);

  middle_point[0] *= scale;
  middle_point[1] *= scale;
  middle_point[2] *= scale;
}

static inline double * find_furthest_vertex(
  double (* restrict vertices)[3], size_t num_vertices,
  double * restrict middle_point) {

  double * x_furthest = NULL;
  double min_dot = 2.0; // actual co-domain is [-1;1]

#ifdef __NEC__
// vectorization of the following loop leads to failue of
// test_interp_method_hcsbb_parallel
// with NEC compiler when 'CFLAGS='-O1 -finline-functions'
#pragma _NEC novector
#endif
  for (size_t i = 0; i < num_vertices; ++i) {

    double curr_dot = vertices[i][0] * middle_point[0] +
                      vertices[i][1] * middle_point[1] +
                      vertices[i][2] * middle_point[2];

    // the dot product is defined by dot = cos(a,b)*|a|*|b|
    // since all our vectors have a length of 1, we directly get the cosine of
    // the angle between the vectors
    // we are only interested in determining which vertex is the furtherest,
    // therefore we do not need the actual angle
    if (curr_dot < min_dot) {
      x_furthest = &(vertices[i][0]);
      min_dot = curr_dot;
    }
  }

  return x_furthest;
}

// computes the triangle
static void compute_triangle(
  double * a, double * x_middle, double * x_t, double (*x_triangle)[3]) {

  // compute the three norm vectors for the planes defined by the edges and the
  // origin
  double edge_norm_vector[3][3];
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j)
      edge_norm_vector[i][j] =
        -1.0 * a[1] * x_middle[j] + a[0] * x_t[3*i+j];
    normalise_vector(edge_norm_vector[i]);
  }

  // the corners of the triangle are the intersections of planes
  // they are perpendicular to the norm vectors of the adjacent edges, hence
  // we can use the cross product to determine them
  // we just have to be carefull with the ordering of the norm vector or else we
  // might get a vector opposite of the actual corner of the triangle
  crossproduct_d(edge_norm_vector[0], edge_norm_vector[1], &(x_triangle[0][0]));
  crossproduct_d(edge_norm_vector[1], edge_norm_vector[2], &(x_triangle[1][0]));
  crossproduct_d(edge_norm_vector[2], edge_norm_vector[0], &(x_triangle[2][0]));
  for (size_t i = 0; i < 3; ++i) normalise_vector(&(x_triangle[i][0]));
}

void yac_compute_bnd_triangle(
  double * vertices, size_t num_vertices, double triangle[][3],
  size_t num_tests) {

  // the middle point of the triangle that is to contain all vertices
  double x_middle[3];
  compute_middle_point((double(*)[3])vertices, num_vertices, x_middle);

  // find the vertex that is furthest away from the middle point
  double * x_furthest =
    find_furthest_vertex((double(*)[3])vertices, num_vertices, x_middle);

  // a triangle whose edge goes through the furthest index and whose edge
  // through the furthest index is perpendicular to the plane defined by
  // x_middle, x_furthest and x_origin will always contain all points

  // now we compute three vectors x_tj (j = 0..2), with the following properties
  // I. all x_tj are orthogonal to x_middle
  // II. x_t0 is in the plane (x_middle, x_furthest, x_origin)
  // III. the angle between x_t0 and x_t1/x_t2 is 120�/-120�
  double x_t[3][3];
  {
    double t[3];
    crossproduct_d(x_middle, x_furthest, t);
    crossproduct_d(t, x_middle, x_t[0]);
  }
  normalise_vector(x_t[0]);
  rotate_vector(x_middle,  (2.0 * M_PI) / 3.0, x_t[0], x_t[1]);
  rotate_vector(x_middle, -(2.0 * M_PI) / 3.0, x_t[0], x_t[2]);

  // we define the opening angle of the triangle to be the angle between
  // x_middle and the vector to the middle of the edge of the triangle

  double min_cos_angle = -2.0; // is an invalid value
  double best_x_t[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  double best_a[2] = {0.0, 0.0};

  double * a = xmalloc(3 * num_vertices * 3 * sizeof(*a));

  // we have one solution that should work. now we check whether different
  // rotations of the current solution give better results
  // we only need to check a rotation range of ]0�,60�[
  double d_rot_angle = M_PI / (double)(3 * num_tests);
  for (size_t test_idx = 0; test_idx < num_tests; ++test_idx) {

    // current rotation angle of triangle to be checked
    double curr_rot_angle = (double)test_idx * d_rot_angle;

    double x_t_[3][3];

    // for the three edges
    for (size_t i = 0; i < 3; ++i) {

      memcpy(a + i*3*num_vertices, vertices, 3*num_vertices * sizeof(*a));

      // rotation of x_ti gives x_ti_
      rotate_vector(x_middle, curr_rot_angle, x_t[i], x_t_[i]);

      double A[3][3];
      // x_o is orthogonal to x_ti_ and x_middle
      crossproduct_d(x_middle, x_t_[i], &A[2][0]/*x_o*/);

      // the vectors x_middle, x_ti_ and x_o define a coordinate system
      // each vector can be define as
      // v = a0*x_middle + a1*x_ti_ + a2*x_o

      // the sine and cosine of the opening angle required for the current
      // triangle to contain a vertex can be computed by:
      // cos_angle = a0 / sqrt(a0*a0 + a1*a1)
      // sin_angle = a1 / sqrt(a0*a0 + a1*a1)

      // to get a0 and a1 we have to solve the following linear system:
      // A*a = v
      // where: A is a matrix consisting of (x_middle, x_ti_, x_o) and
      //        a is the vector (a0, a1, a2)

      lapack_int n = 3, nrhs = (lapack_int) num_vertices, lda = n, ldx = n, ipiv[3];
      A[0][0] = x_middle[0], A[0][1] = x_middle[1], A[0][2] = x_middle[2];
      A[1][0] = x_t_[i][0], A[1][1] = x_t_[i][1], A[1][2] = x_t_[i][2];

      // we use LAPACK to solve the linear system for all vertices at once
      // initially a contains all vertices, after the call to LAPACKE_dgesv it
      // contains the results
      YAC_ASSERT(
        !LAPACKE_dgesv(
          LAPACK_COL_MAJOR, n, nrhs, &A[0][0], lda,
          &ipiv[0], a + i * 3 * num_vertices, ldx),
        "ERROR: internal error (could not solve linear 3x3 system)")
    }

    // for each test_idx we actually check two triangles, which have a rotation
    // angle of 180� between themselves
    // each edge of the two triangles has an opposite counterpart. a vertex will
    // influence the required opening angle for one of the two edges. the
    // effected edge is determined base on sign of the sine of the required
    // opening angle (see flag)
    double curr_min_cos_angle[2] = {2.0, 2.0};
    double * curr_x_furthest[2] = {NULL, NULL};
    double curr_best_a[2][2];

    // for all three edges and each vertex
    for (size_t i = 0; i < 3 * num_vertices; ++i) {

      int flag = a[1+3*i] < 0.0;
      double scale = 1.0 / sqrt(a[0+3*i]*a[0+3*i] + a[1+3*i]*a[1+3*i]);
      double cos_angle = a[0+3*i] * scale;

      if (cos_angle < curr_min_cos_angle[flag]) {
        curr_min_cos_angle[flag] = cos_angle;
        curr_x_furthest[flag] = vertices + 3 * (i%num_vertices);
        curr_best_a[flag][0] = cos_angle;
        curr_best_a[flag][1] = a[1+3*i]*scale;
      }
    }

    int flag = curr_min_cos_angle[0] < curr_min_cos_angle[1];

#ifdef DEBUG
    // debug output
    if (curr_x_furthest[flag] != NULL) {
      double x_triangle[3][3];
      compute_triangle(&curr_best_a[flag][0], &x_middle[0], &x_t_[0][0],
                       &x_triangle[0][0]);
      print_bnd_triangle(vertices, num_vertices, &x_triangle[0][0], test_idx);
    }
#endif // DEBUG

    // if the current triangle is smaller than the best solution till now
    if ((curr_x_furthest[flag] != NULL) &&
        (curr_min_cos_angle[flag] > min_cos_angle)) {
      memcpy(&best_x_t[0][0], &x_t_[0][0], 9 * sizeof(x_t_[0][0]));
      min_cos_angle = curr_min_cos_angle[flag];
      x_furthest = curr_x_furthest[flag];
      best_a[0] = curr_best_a[flag][0];
      best_a[1] = curr_best_a[flag][1];
    }
  }

  free(a);

  YAC_ASSERT(
    min_cos_angle > -2.0,
    "ERROR: internal error (could not fine any matching triangle)")

  // to get the actual triangle we have to compute the intersections of the
  // edges
  compute_triangle(&best_a[0], &x_middle[0], &best_x_t[0][0], triangle);
}

