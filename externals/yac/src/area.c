/**
 * @file area.c
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
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "area.h"
#include "grid.h"
#include "clipping.h"
#include "geometry.h"
#include "utils.h"
#include "ensure_array_size.h"

static inline double scalar_product(double a[], double b[]);

/* ----------------------------------- */

double yac_triangle_area ( struct grid_cell cell ) {

  /* taken from the ICON code, mo_base_geometry.f90
     provided by Luis Kornblueh, MPI-M. */

  const double tol = 1e-18;

  double s01, s12, s20;
  double ca1, ca2, ca3;
  double a1, a2, a3;

  double * triangle[3];
  double u01[3], u12[3], u20[3];

  YAC_ASSERT(
    cell.num_corners == 3, "ERROR(yac_triangle_area): cell is not a triangle")

  triangle[0] = cell.coordinates_xyz[0];
  triangle[1] = cell.coordinates_xyz[1];
  triangle[2] = cell.coordinates_xyz[2];

  /* First, compute cross products Uij = Vi x Vj. */

  crossproduct_kahan(triangle[0], triangle[1], u01);
  crossproduct_kahan(triangle[1], triangle[2], u12);
  crossproduct_kahan(triangle[2], triangle[0], u20);

  /*  Normalize Uij to unit vectors. */

  s01 = scalar_product(u01, u01);
  s12 = scalar_product(u12, u12);
  s20 = scalar_product(u20, u20);

  /* Test for a degenerated triangle associated with collinear vertices. */

  if ( fabs(s01) < tol ||
       fabs(s12) < tol ||
       fabs(s20) < tol )
    return 0.0;

  s01 = sqrt(s01);
  s12 = sqrt(s12);
  s20 = sqrt(s20);

  for (int m = 0; m < 3; m++ ) {
    u01[m] = u01[m]/s01;
    u12[m] = u12[m]/s12;
    u20[m] = u20[m]/s20;
  }

  /*  Compute interior angles Ai as the dihedral angles between planes:

      CA1 = cos(A1) = -<U01,U20>
      CA2 = cos(A2) = -<U12,U01>
      CA3 = cos(A3) = -<U20,U12>

  */

  ca1 = -u01[0]*u20[0]-u01[1]*u20[1]-u01[2]*u20[2];
  ca2 = -u12[0]*u01[0]-u12[1]*u01[1]-u12[2]*u01[2];
  ca3 = -u20[0]*u12[0]-u20[1]*u12[1]-u20[2]*u12[2];

  if ( ca1 < -1.0 ) ca1 = -1.0;
  if ( ca1 >  1.0 ) ca1 =  1.0;
  if ( ca2 < -1.0 ) ca2 = -1.0;
  if ( ca2 >  1.0 ) ca2 =  1.0;
  if ( ca3 < -1.0 ) ca3 = -1.0;
  if ( ca3 >  1.0 ) ca3 =  1.0;

  a1 = acos(ca1);
  a2 = acos(ca2);
  a3 = acos(ca3);

  /*  Compute areas = a1 + a2 + a3 - pi.

      here for a unit sphere: */

  // return MAX(( a1+a2+a3-M_PI ) * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS, 0.0);
  return MAX( (a1+a2+a3-M_PI) , 0.0 );
}

/* ----------------------------------- */

static inline struct sin_cos_angle
tri_area_quarter_angle(struct sin_cos_angle angle) {

  // the angle passed to this routine should be in the range [0;3*PI/4],
  // in case the angle is outside this, we have to assume that this is due to
  // numerical inaccuracy or
  // tri_area was called for a too big triangle...has to be really big...

  if (compare_angles(angle, SIN_COS_7_M_PI_4) >= 0) return SIN_COS_ZERO;
  else return quarter_angle(angle);
}

/** area of a spherical triangle based on L'Huilier's Theorem
  *
  * source code is taken from code by Robert Oehmke of Earth System Modeling
  * Framework (www.earthsystemmodeling.org)
  *
  * it has been extended by a more accurate computation of vector angles
  *
  * the license statement for this routine is as follows:
  * Earth System Modeling Framework
  * Copyright 2002-2013, University Corporation for Atmospheric Research,
  * Massachusetts Institute of Technology, Geophysical Fluid Dynamics
  * Laboratory, University of Michigan, National Centers for Environmental
  * Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
  * NASA Goddard Space Flight Center.
  * Licensed under the University of Illinois-NCSA License.
  *
  * \remark all edges are on great circle
  */
static double tri_area_(
  struct sin_cos_angle angle_a,
  struct sin_cos_angle angle_b,
  struct sin_cos_angle angle_c) {

  if (compare_angles(angle_a, SIN_COS_LOW_TOL) < 0) return 0.0;
  if (compare_angles(angle_b, SIN_COS_LOW_TOL) < 0) return 0.0;
  if (compare_angles(angle_c, SIN_COS_LOW_TOL) < 0) return 0.0;

  double sin_sin = angle_a.sin * angle_b.sin;
  double sin_cos = angle_a.sin * angle_b.cos;
  double cos_sin = angle_a.cos * angle_b.sin;
  double cos_cos = angle_a.cos * angle_b.cos;

  double sin_sin_sin = sin_sin * angle_c.sin;
  double sin_sin_cos = sin_sin * angle_c.cos;
  double sin_cos_sin = sin_cos * angle_c.sin;
  double sin_cos_cos = sin_cos * angle_c.cos;
  double cos_sin_sin = cos_sin * angle_c.sin;
  double cos_sin_cos = cos_sin * angle_c.cos;
  double cos_cos_sin = cos_cos * angle_c.sin;
  double cos_cos_cos = cos_cos * angle_c.cos;

  double t_sin_a = sin_sin_sin - sin_cos_cos;
  double t_sin_b = cos_sin_cos + cos_cos_sin;
  double t_sin_c = sin_sin_sin + sin_cos_cos;
  double t_sin_d = cos_sin_cos - cos_cos_sin;
  double t_cos_a = cos_cos_cos - cos_sin_sin;
  double t_cos_b = sin_sin_cos + sin_cos_sin;
  double t_cos_c = cos_cos_cos + cos_sin_sin;
  double t_cos_d = sin_sin_cos - sin_cos_sin;

  struct sin_cos_angle t_angle[4] = {
    tri_area_quarter_angle(
      (struct sin_cos_angle){.sin = - t_sin_a + t_sin_b,
                             .cos = + t_cos_a - t_cos_b}),
    tri_area_quarter_angle(
      (struct sin_cos_angle){.sin = + t_sin_a + t_sin_b,
                             .cos = + t_cos_a + t_cos_b}),
    tri_area_quarter_angle(
      (struct sin_cos_angle){.sin = + t_sin_c - t_sin_d,
                             .cos = + t_cos_c + t_cos_d}),
    tri_area_quarter_angle(
      (struct sin_cos_angle){.sin = + t_sin_c + t_sin_d,
                             .cos = + t_cos_c - t_cos_d})};

  double t = (t_angle[0].sin*t_angle[1].sin*t_angle[2].sin*t_angle[3].sin) /
             (t_angle[0].cos*t_angle[1].cos*t_angle[2].cos*t_angle[3].cos);

  return 4.0 * atan(sqrt(fabs(t)));
}

static double tri_area(double u[3], double v[3], double w[3]) {

  return
    tri_area_(get_vector_angle_2(u,v),
              get_vector_angle_2(u,w),
              get_vector_angle_2(w,v));
}

/* ----------------------------------- */

static inline int compute_norm_vector(double a[], double b[], double norm[]) {

  crossproduct_kahan(a, b, norm);

  double scale = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);

  if (scale <= yac_angle_tol) return 1;

  scale = 1.0 / scale;

  norm[0] *= scale;
  norm[1] *= scale;
  norm[2] *= scale;

  return 0;
}

static double inline XYZtoLon(double a[3]) {
  return atan2(a[1] , a[0]);
}

static double
lat_edge_correction(double base_vec[3], double a[3], double b[3]) {

  double const tol = 1e-8;

  YAC_ASSERT(
    fabs(a[2]-b[2]) <= tol, "ERROR: latitude of both corners is not identical")

  double h = fabs(a[2]);

  // if we are at the equator or at a pole
  if (h < tol || fabs(1.0 - h) < tol)
    return 0.0;

  double lat_area = fabs((1.0 - h) * get_angle(XYZtoLon(a), XYZtoLon(b)));

  double pole[3] = {0, 0, (a[2] >= 0.0)?1.0:-1.0};
  double gc_area = tri_area(a, b, pole);

  double correction = MAX(lat_area - gc_area, 0.0);

  double middle_lat[3] = {a[0]+b[0], a[1]+b[1], a[2]};
  double scale = sqrt(middle_lat[0]*middle_lat[0]+middle_lat[1]*middle_lat[1]);

  YAC_ASSERT(fabs(scale) >= 1e-18, "internal error")

  scale = sqrt(1.0 - a[2]*a[2]) / scale;

  middle_lat[0] *= scale;
  middle_lat[1] *= scale;

  double norm_ab[3];

  // if the angle between a and b is to small to compute a norm vector
  if (compute_norm_vector(a, b, norm_ab)) return 0.0;

  double scalar_base = scalar_product(norm_ab, base_vec);
  double scalar_middle_lat = scalar_product(norm_ab, middle_lat);

  // if the base point is on the same plane as a and b
  if (fabs(scalar_base) < 1e-11) {

    double norm_middle[3];

    if (compute_norm_vector(middle_lat, pole, norm_middle)) return 0.0;

    double scalar_a = scalar_product(norm_middle, a);

    if (scalar_a > 0)
      return correction;
    else
      return - correction;

  } else {

    if (scalar_middle_lat < 0)
      return correction;
    else
      return - correction;
  }
}



double yac_pole_area(struct grid_cell cell) {

  size_t const M = cell.num_corners;

  double coordinates_x[M];
  double coordinates_y[M];

  for (size_t i = 0; i < M; ++i)
    XYZtoLL(cell.coordinates_xyz[i], &(coordinates_x[i]), &(coordinates_y[i]));

  double area = 0.0;

  if (M < 2) return 0.0;

  int closer_to_south_pole = coordinates_y[0] < 0;

  double pole_vec[3] = {0, 0, (closer_to_south_pole)?-1.0:1.0};

  // it would also be possible to use the equator instead
  // of the poles as the baseline
  // this could be used as special case for cell close
  // the equator (were the other method is probably most
  // inaccurate)

  for (size_t i = 0; i < M; ++i) {

    // if one of the points it at the pole
    if (fabs(fabs(coordinates_y[i]) - M_PI_2) < 1e-12) continue;
    if (fabs(fabs(coordinates_y[(i+1)%M]) - M_PI_2) < 1e-12) {
      ++i; // we can skip the next edge
      continue;
    }

    YAC_ASSERT(
      (cell.edge_type[i] == GREAT_CIRCLE_EDGE) ||
      (cell.edge_type[i] == LON_CIRCLE_EDGE) ||
      (cell.edge_type[i] == LAT_CIRCLE_EDGE),
      "ERROR: unsupported edge type")

    if (cell.edge_type[i] == GREAT_CIRCLE_EDGE ||
        cell.edge_type[i] == LON_CIRCLE_EDGE) {

      double * a;
      double * b;

      a = cell.coordinates_xyz[i];
      b = cell.coordinates_xyz[(i+1)%M];

      double edge_direction = a[0]*b[1]-a[1]*b[0]; // 3. component of cross product

      // if the edge is nearly on a circle of longitude
      if (fabs(edge_direction) < 1e-12) continue;

      double tmp_area = tri_area(a, b, pole_vec);

      // or the other way round
      if (edge_direction > 0) area -= tmp_area;
      else                    area += tmp_area;

    } else {

      // the area of a sphere cap is:
      // A = 2 * PI * r * h (where r == 1 and h == 1 - cos(d_lat))
      // scaled with the longitude angle this is:
      // A' = (d_lon / (2 * PI)) * A
      // => A' = d_lon * (1 - cos(d_lat))

      double d_lon = get_angle(coordinates_x[i], coordinates_x[(i+1)%M]);
      double d_lat = M_PI_2;

      if (closer_to_south_pole)
        d_lat += coordinates_y[i];
      else
        d_lat -= coordinates_y[i];

      double h = 1 - cos(d_lat);

      area += d_lon * h;

    }
  }
  // return fabs(area) * YAC_EARTH_RADIUS * YAC_EARTH_RADIUS;
  return fabs(area);
}

double yac_planar_3dcell_area (struct grid_cell cell) {

 /*
  * source code is originally based on http://gaim.umbc.edu/2010/06/03/polygon-area/
  *
  */

  double norm[3] = {0,0,0};
  size_t M = cell.num_corners;

  if (M < 3) return 0.0;

  for (size_t i0 = M - 1, i1 = 0; i1 < M; i0 = i1, ++i1) {
    norm[0] += cell.coordinates_xyz[i0][1]*cell.coordinates_xyz[i1][2] -
               cell.coordinates_xyz[i1][1]*cell.coordinates_xyz[i0][2];
    norm[1] += cell.coordinates_xyz[i0][2]*cell.coordinates_xyz[i1][0] -
               cell.coordinates_xyz[i1][2]*cell.coordinates_xyz[i0][0];
    norm[2] += cell.coordinates_xyz[i0][0]*cell.coordinates_xyz[i1][1] -
               cell.coordinates_xyz[i1][0]*cell.coordinates_xyz[i0][1];
  };

  return 0.5 * sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
}

 /*
  * source code is originally based on code by Robert Oehmke of Earth System Modeling
  * Framework (www.earthsystemmodeling.org)
  *
  * it has be extended to support YAC data structures and two types of
  * grid cell edges (great circle and circle of latitude)
  *
  */
double yac_huiliers_area (struct grid_cell cell) {

  size_t M = cell.num_corners;

  if (M < 2) return 0.0;

  int lat_flag = 0;

  for (size_t i = 0; i < M; i++)
    lat_flag |= cell.edge_type[i] == LAT_CIRCLE_EDGE;

  if (M == 3 && !lat_flag)
    return tri_area(cell.coordinates_xyz[0],
                    cell.coordinates_xyz[1],
                    cell.coordinates_xyz[2]);

  // sum areas around cell
  double area = 0.0;

  for (size_t i = 2; i < M; ++i) {

    double tmp_area = tri_area(cell.coordinates_xyz[0],
                               cell.coordinates_xyz[i-1],
                               cell.coordinates_xyz[i]);

    double norm[3];

    if (compute_norm_vector(cell.coordinates_xyz[i-1],
                            cell.coordinates_xyz[i], norm)) continue;

    double scalar_base = scalar_product(norm, cell.coordinates_xyz[0]);

    if (scalar_base > 0) area += tmp_area;
    else area -= tmp_area;
  }

  // if there is at least one latitude circle edge
  if (lat_flag) {

    for (size_t i = 0; i < M; ++i) {

      if (cell.edge_type[i] == LAT_CIRCLE_EDGE) {

        size_t i_ = (i+1)%cell.num_corners;

        area += lat_edge_correction(cell.coordinates_xyz[0],
                                    cell.coordinates_xyz[i],
                                    cell.coordinates_xyz[i_]);
      }
    }
  }

  return fabs(area);
}

static double tri_area_info(
  double ref[3], double a[3], double b[3],
  double * barycenter, double sign) {

  double * corners[3] = {b, ref, a};
  double cross[3][3];
  struct sin_cos_angle angles[3];
  for (int i = 0, j = 2; i < 3; j = i++) {
    crossproduct_kahan(corners[j], corners[i], cross[i]);
    angles[i] =
      sin_cos_angle_new(sqrt(cross[i][0]*cross[i][0] +
                             cross[i][1]*cross[i][1] +
                             cross[i][2]*cross[i][2]),
                       scalar_product(corners[j], corners[i]));
  }

  double area = tri_area_(angles[0], angles[1], angles[2]);

  // the barycenter of the triangle is given by the sum edge norm vector
  // scaled by half of the associated edge length
  for (int i = 0; i < 3; ++i) {
    double scale = 0.5 * compute_angle(angles[i]) / angles[i].sin * sign;
    for (int j = 0; j < 3; ++j)
      barycenter[j] += cross[i][j] * scale;
  }

  if (scalar_product(ref, cross[0]) < 0.0) area = -area;

  return area * sign;
}

double yac_huiliers_area_info(
  struct grid_cell cell, double * barycenter, double sign) {

  size_t M = cell.num_corners;

  if (M < 2) return 0.0;

  int lat_flag = 0;

  for (size_t i = 0; i < M; i++)
    lat_flag |= cell.edge_type[i] == LAT_CIRCLE_EDGE;

  if (M == 3 && !lat_flag)
    return tri_area_info(cell.coordinates_xyz[0],
                         cell.coordinates_xyz[1],
                         cell.coordinates_xyz[2],
                         barycenter, sign);

  // sum areas around cell
  double area = 0.0;

  for (size_t i = 2; i < M; ++i) {
    area +=
      tri_area_info(
        cell.coordinates_xyz[0],
        cell.coordinates_xyz[i-1],
        cell.coordinates_xyz[i],
        barycenter, sign);
  }

  // if there is at least one latitude circle edge
  if (lat_flag) {

    for (size_t i = 0; i < M; ++i) {

      if (cell.edge_type[i] == LAT_CIRCLE_EDGE) {

        size_t i_ = (i+1)%cell.num_corners;

        area += lat_edge_correction(cell.coordinates_xyz[0],
                                    cell.coordinates_xyz[i],
                                    cell.coordinates_xyz[i_]) * sign;
      }
    }
  }

  return area;
}

/* ----------------------------------- */

static inline double scalar_product(double a[], double b[]) {
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

/* ----------------------------------- */

