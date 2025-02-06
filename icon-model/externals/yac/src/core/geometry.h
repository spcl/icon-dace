// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

/** \example test_geometry.c
 * This contains some tests for basic routines of \ref geometry.h.
 */

/** \example test_angle.c
 * This contains some tests for angle calculations of \ref geometry.h.
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#ifdef HAVE_CONFIG_H
// We include 'config.h' in this header file (which is a bad practice)
// because we need the definition of the 'restrict' keyword for the
// inlined functions.
#include "config.h"
#endif

#include <math.h>
#include <float.h>

#include "basic_grid.h"
#include "grid_cell.h"
#include "utils_core.h"

#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

// angle tolerance

#define yac_angle_tol         (1e-9)
#define yac_sq_angle_tol      (1e-9 * 1e-9)
#define yac_cos_angle_tol     (0.9999999999999999995)
#define yac_angle_low_tol     (1e-11)
#define yac_cos_angle_low_tol (1.0)

struct sin_cos_angle {
  double sin, cos;
};

static const struct sin_cos_angle SIN_COS_ZERO = {0.0, 1.0};
static const struct sin_cos_angle SIN_COS_TOL = {yac_angle_tol, yac_cos_angle_tol};
static const struct sin_cos_angle SIN_COS_LOW_TOL = {yac_angle_low_tol, yac_cos_angle_low_tol};

static const struct sin_cos_angle SIN_COS_M_PI_2 = {1.0, 0.0}; /* PI/2 */
static const struct sin_cos_angle SIN_COS_M_PI = {0.0, -1.0}; /* PI */
static const struct sin_cos_angle SIN_COS_7_M_PI_4 = {-0.707106781186547524401, +0.707106781186547524401}; /* (7 * PI)/4 */

/**
 * defines a spherical cap
 */
struct bounding_circle {

   //! the middle point of the spherical cap in cartesian coordinates
   //! (is a unit vector)
   double base_vector[3];
   //! angle between the middle point and the boundary of the spherical cap
   struct sin_cos_angle inc_angle;
   //! squared Euclidean distance between the base vector and the edge of the
   //! bounding circle (chord of inc_angle)
   double sq_crd; // sq_crd = (chord)^2
};

/** \example test_circle.c
 * This contains some examples on how to use the yac_circle interfaces
 * routines.
 */

enum yac_circle_type {
   GREAT_CIRCLE = YAC_GREAT_CIRCLE_EDGE,
   LAT_CIRCLE   = YAC_LAT_CIRCLE_EDGE,
   LON_CIRCLE   = YAC_LON_CIRCLE_EDGE,
   POINT,
};

// each circle partitions the sphere into an inside and an outside part
// (the point is the exception; everything except the point is outside)
struct yac_circle {
  enum yac_circle_type type;
  union {
    struct {
      // the norm vector points into the inside direction
      double norm_vector[3];
    } gc, lon;
    struct {
      int north_is_out;
      double z;
    } lat;
    struct {
      double vec[3];
    } p;
  } data;
};

/** \example test_point_in_cell.c
 * This contains examples on how to use point_in_cell.
 */

/**
 * checks whether a given point is within a given cell \n
 * @param[in] point_coords
 * @param[in] cell
 * @return 0 if the point is not in the cell
 */
int yac_point_in_cell (double point_coords[3], struct yac_grid_cell cell);

/**
 * checks whether a given point is within a given cell \n
 * @param[in] point_coords
 * @param[in] cell
 * @param[in] bnd_circle
 * @return 0 if the point is not in the cell
 */
int yac_point_in_cell2 (double point_coords[3], struct yac_grid_cell cell,
                        struct bounding_circle bnd_circle);

/**
 * computes the angle between two longitude coordinates (in rad) \n
 * takes into account that longitude have a period of 2 pi
 * @param[in] a_lon
 * @param[in] b_lon
 * @return angle between both coordinates (in rad)
 */
static inline double get_angle (double a_lon, double b_lon) {
   double diff = a_lon - b_lon;
#if defined(CDO)
   return diff - lround(diff / (2.0 * M_PI)) * (2.0 * M_PI);
#else
   return diff - round(diff / (2.0 * M_PI)) * (2.0 * M_PI);
#endif
}

/**
 * computes the angle between two longitude coordinates (in deg) \n
 * takes into account that longitude have a period of 360
 * @param[in] a_lon
 * @param[in] b_lon
 * @return angle between both coordinates (in rad)
 */
static inline double get_angle_deg (double a_lon, double b_lon) {
   double diff = a_lon - b_lon;
#if defined(CDO)
   return diff - lround(diff / 360.0) * (360.0);
#else
   return diff - round(diff / 360.0) * (360.0);
#endif
}

/**
 * computes the intersection points of two edges
 * @param[in]  edge_type_a type of edge a
 * @param[in]  a           first point of edge a
 * @param[in]  b           second point of edge a
 * @param[in]  edge_type_b type of edge b
 * @param[in]  c           first point of edge b
 * @param[in]  d           second point of edge b
 * @param[out] p           first intersection point
 * @param[out] q           second intersection point
 * @return  0 if the intersection points are neither on edge a or b \n
 *         -1 if the two circles do not intersect \n
 *          1st bit will be set if p is on edge a \n
 *          2nd bit will be set if q is on edge a \n
 *          3rd bit will be set if p is on edge b \n
 *          4th bit will be set if q is on edge b \n
 *          5th bit will be set if both circles are identically
 *
 * \remarks if -1 is returned neither p or q is set
 * \remarks if the two circles only have one intersection point,
 *          p and q will be identically, but only the p bits will be set
 */
int yac_intersect_vec (
  enum yac_edge_type edge_type_a, double const a[3], double const b[3],
  enum yac_edge_type edge_type_b, double const c[3], double const d[3],
  double p[3], double q[3]);

/** \example test_cell_bnd_circle.c
 * These are some examples on how to use \ref yac_get_cell_bounding_circle.
 */

/**
 * gets the bounding circle for a grid cell
 * @param[in] cell grid cell (coordinates have to be in radian)
 * @param[out] bnd_circle bounding circle of the grid cell
 */
void yac_get_cell_bounding_circle(struct yac_grid_cell cell,
                                  struct bounding_circle * bnd_circle);

/**
 * computes the circumscribe circle for a triangle on the sphere
 * @param[in]  a          coordinates of first point (xyz)
 * @param[in]  b          coordinates of second point (xyz)
 * @param[in]  c          coordinates of thrid point (xyz)
 * @param[out] bnd_circle circumscribe circle
 * @remark it is assumed that all three edges of the triangle are great circles
 */
void yac_get_cell_circumscribe_circle_unstruct_triangle(
   double a[3], double b[3], double c[3], struct bounding_circle * bnd_circle);

/**
 * computes the smallest bounding circle for a triangle on the sphere
 * @param[in]  a          coordinates of first point (xyz)
 * @param[in]  b          coordinates of second point (xyz)
 * @param[in]  c          coordinates of thrid point (xyz)
 * @param[out] bnd_circle bounding circle
 * @remark it is assumed that all three edges of the triangle are great circles
 */
void yac_get_cell_bounding_circle_unstruct_triangle(
   double a[3], double b[3], double c[3], struct bounding_circle * bnd_circle);

/**
 * computes the circumscribe circle for a quad on the sphere
 * @param[in]  a          coordinates of first point (xyz)
 * @param[in]  b          coordinates of second point (xyz)
 * @param[in]  c          coordinates of thrid point (xyz)
 * @param[out] bnd_circle circumscribe circle
 * @remark it is assumed that all edges of the quad are either circles of
 *         longitude or latitude
 */
void yac_get_cell_circumscribe_circle_reg_quad(
   double a[3], double b[3], double c[3],
   struct bounding_circle * bnd_circle);

/**
 * computes the smallest bounding circle for a triangle on the sphere
 * @param[in]  a          coordinates of first point (xyz)
 * @param[in]  b          coordinates of second point (xyz)
 * @param[in]  c          coordinates of thrid point (xyz)
 * @param[out] bnd_circle bounding circle
 * @remark it is assumed that all edges of the quad are either circles of
 *         longitude or latitude
 */
void yac_get_cell_bounding_circle_reg_quad(
   double a[3], double b[3], double c[3],
   struct bounding_circle * bnd_circle);

/**
 * checks whether two extents overlap
 * @param[in] extent_a bounding circle
 * @param[in] extent_b bounding circle
 * @return 0 if the bounding circles do not overlap
 */
int yac_extents_overlap(struct bounding_circle * extent_a,
                        struct bounding_circle * extent_b);

//! computes square of the lenght of the vector ab
static inline double sq_len_diff_vec(double const a[3], double const b[3]) {
  double ab[3] = {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
  return ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
}

static inline double compute_sq_crd(struct sin_cos_angle angle) {

  double temp = 1.0 - angle.cos;
  return temp * temp + angle.sin * angle.sin;
}

/**
 * checks whether a point is within a bounding circle
 * @param[in] point_vector point to be checked
 * @param[in] bnd_circle bounding circle
 * @return 0 if point is not within the bounding circle
 */
static inline int yac_point_in_bounding_circle_vec(
  double point_vector[3], struct bounding_circle * bnd_circle) {

  if (bnd_circle->sq_crd == DBL_MAX)
    bnd_circle->sq_crd = compute_sq_crd(bnd_circle->inc_angle);

  return
    bnd_circle->sq_crd >=
    sq_len_diff_vec(bnd_circle->base_vector, point_vector);
}

/**
 * converts lon-lat coordinates into xyz ones
 *
 * Further information:
 * http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations
 *
 * @param[in]  lon   longitude coordinates in radian
 * @param[in]  lat   latitude coordinates in radian
 * @param[out] p_out xyz coordinates
 */
static inline void LLtoXYZ(double lon, double lat, double p_out[]) {

   while (lon < -M_PI) lon += 2.0 * M_PI;
   while (lon >= M_PI) lon -= 2.0 * M_PI;

   double cos_lat = cos(lat);
   p_out[0] = cos_lat * cos(lon);
   p_out[1] = cos_lat * sin(lon);
   p_out[2] = sin(lat);
}

/**
 * converts lon-lat coordinates into xyz ones
 * @param[in]  lon   longitude coordinates in deg
 * @param[in]  lat   latitude coordinates in deg
 * @param[out] p_out xyz coordinates
 */
static inline void LLtoXYZ_deg(double lon, double lat, double p_out[]) {
   LLtoXYZ(lon*YAC_RAD, lat*YAC_RAD, p_out);
}

/**
 * converts lon-lat coordinates into xyz ones
 *
 * Further information:
 * http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations
 *
 * @param[in]  p_in xyz coordinates
 * @param[out] lon  longitude coordinate in radian
 * @param[out] lat  latitude coordinate in radian
 */
static inline void XYZtoLL (double const p_in[], double * lon, double * lat) {

   *lon = atan2(p_in[1] , p_in[0]);
   *lat = M_PI_2 - acos(p_in[2]);
}

// computation of the determinant of a 2x2 matrix using Kahan summation
static inline double det2_kahan(double a, double b, double c, double d) {
  double bc = b*c;
  double e_bc = fma(b,c,-bc); // the rounding error of the multiplication
  double det = fma(a,d,-bc);
  return det + e_bc;
}

static inline void crossproduct_kahan (
  double const a[], double const b[], double cross[]) {

  /* crossproduct in Cartesian coordinates */

  cross[0] = det2_kahan(a[1], a[2], b[1], b[2]);
  cross[1] = det2_kahan(a[2], a[0], b[2], b[0]);
  cross[2] = det2_kahan(a[0], a[1], b[0], b[1]);

}

/**
 * for small angles <= 1e-?8? the crossproduct is inaccurate\n
 * use \ref crossproduct_kahan for these cases
 */
static inline void crossproduct_d (
  const double a[], const double b[], double cross[]) {

/* crossproduct in Cartesian coordinates */

   cross[0] = a[1] * b[2] - a[2] * b[1];
   cross[1] = a[2] * b[0] - a[0] * b[2];
   cross[2] = a[0] * b[1] - a[1] * b[0];
}

/**
 * computes the great circle distance in rad for two points given in xyz coordinates
 * taken from http://johnblackburne.blogspot.de/2012/05/angle-between-two-3d-vectors.html
 * @param[in] a point coordinates of point a
 * @param[in] b point coordinates of point b
 * @return great circle distance in rad between both points
 */
static inline double get_vector_angle(double const a[3], double const b[3]) {

   double dot_product = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

   // the acos most accurate in the range [-M_SQRT1_2;M_SQRT1_2]
   if (fabs(dot_product) > M_SQRT1_2) {

      double cross_ab[3];

      crossproduct_kahan(a, b, cross_ab);

      double asin_tmp = asin(sqrt(cross_ab[0]*cross_ab[0]+
                                  cross_ab[1]*cross_ab[1]+
                                  cross_ab[2]*cross_ab[2]));

      if (dot_product > 0.0) // if the angle is smaller than (PI / 2)
         return MAX(asin_tmp,0.0);
      else
         return MIN(M_PI - asin_tmp, M_PI);

   } else {

      return acos(dot_product);
   }

   /*
   // this solution is simpler, but has a much worse performance
   double cross[3], dot, cross_abs;

   crossproduct_kahan(a, b, cross);
   dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
   cross_abs = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

   return fabs(atan2(cross_abs, dot));
   */
}

static inline double clamp_abs_one(double val) {
  return val > -1.0 ? (val < 1.0 ? val : 1.0) : -1.0;
}

static inline struct sin_cos_angle sin_cos_angle_new(double sin, double cos) {

  struct sin_cos_angle angle;

  angle.sin = clamp_abs_one(sin);
  angle.cos = clamp_abs_one(cos);

  return angle;
}

static inline struct sin_cos_angle get_vector_angle_2(
  double const a[3], double const b[3]) {

  double cross_ab[3];
  crossproduct_kahan(a, b, cross_ab);

  return sin_cos_angle_new(sqrt(cross_ab[0]*cross_ab[0] +
                                cross_ab[1]*cross_ab[1] +
                                cross_ab[2]*cross_ab[2]),
                           a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

// works for angles in the range [0;2*PI[
static inline int compare_angles(
  struct sin_cos_angle a, struct sin_cos_angle b) {

  // there are 5 sections:
  // 0: 0      <= angle < PI/4
  // 1: PI/4   <= angle < 3*PI/4
  // 2: 3*PI/4 <= angle < 5*PI/4
  // 3: 5*PI/4 <= angle < 7*PI/4
  // 4: 7*PI/4 <= angle < 2*PI
  int t_a = fabs(a.cos) <= M_SQRT1_2;
  int t_b = fabs(b.cos) <= M_SQRT1_2;
  int a_section = t_a | ((((a.sin < 0.0) & t_a) |
                          ((a.cos < 0.0) & (fabs(a.sin) < M_SQRT1_2))) << 1);
  int b_section = t_b | ((((b.sin < 0.0) & t_b) |
                          ((b.cos < 0.0) & (fabs(b.sin) < M_SQRT1_2))) << 1);
  // if current section is 0, then it could be actually section 4
  if (!a_section) a_section = (a.sin < 0.0) << 2;
  if (!b_section) b_section = (b.sin < 0.0) << 2;

  if (a_section != b_section)
    return (a_section > b_section) - (a_section < b_section);

  int ret;

  switch (a_section) {
    case(0):
    case(4):
    default:
      ret = (b.sin < a.sin + yac_angle_low_tol) -
            (a.sin < b.sin + yac_angle_low_tol);
      if (ret) return ret;
      else {
        ret = (a.cos < b.cos + yac_angle_low_tol) -
              (b.cos < a.cos + yac_angle_low_tol);
        if (a.sin >= 0.0) return ret;
        else return -ret;
      }
    case(1):
      ret = (a.cos < b.cos + yac_angle_low_tol) -
            (b.cos < a.cos + yac_angle_low_tol);
      if (ret) return ret;
      else {
        ret = (b.sin < a.sin + yac_angle_low_tol) -
              (a.sin < b.sin + yac_angle_low_tol);
        if (a.cos >= 0.0) return ret;
        else return -ret;
      }
    case(2):
      ret = (a.sin < b.sin + yac_angle_low_tol) -
            (b.sin < a.sin + yac_angle_low_tol);
      if (ret) return ret;
      else {
        ret = (a.cos < b.cos + yac_angle_low_tol) -
              (b.cos < a.cos + yac_angle_low_tol);
        if (a.sin >= 0.0) return ret;
        else return -ret;
      }
    case(3):
      ret = (b.cos < a.cos + yac_angle_low_tol) -
            (a.cos < b.cos + yac_angle_low_tol);
      if (ret) return ret;
      else {
        ret = (a.sin < b.sin + yac_angle_low_tol) -
              (b.sin < a.sin + yac_angle_low_tol);
        if (a.cos <= 0.0) return ret;
        else return -ret;
      }
  }
}

//! this routine converts a struct sin_cos_angle into a double, which allows to
//! easily compare angles
//! 0: 0      <= angle < PI/4   => ret = 0 * M_SQRT1_2 + sin(angle)
//! 1: PI/4   <= angle < 3*PI/4 => ret = 2 * M_SQRT1_2 - cos(angle)
//! 2: 3*PI/4 <= angle < 5*PI/4 => ret = 4 * M_SQRT1_2 - sin(angle)
//! 3: 5*PI/4 <= angle < 7*PI/4 => ret = 6 * M_SQRT1_2 + cos(angle)
//! 4: 7*PI/4 <= angle < 2*PI   => ret = 8 * M_SQRT1_2 + sin(angle)
static inline double sin_cos_angle_to_dble(struct sin_cos_angle angle) {

  int t = fabs(angle.cos) <= M_SQRT1_2;
  int section = t | ((((angle.sin < 0.0) & t) |
                      ((angle.cos < 0.0) & (fabs(angle.sin) < M_SQRT1_2))) << 1);
  if (!section) section = (angle.sin < 0.0) << 2;

  switch (section) {
    default:
    case(0): return 0.0 * M_SQRT1_2 + angle.sin;
    case(1): return 2.0 * M_SQRT1_2 - angle.cos;
    case(2): return 4.0 * M_SQRT1_2 - angle.sin;
    case(3): return 6.0 * M_SQRT1_2 + angle.cos;
    case(4): return 8.0 * M_SQRT1_2 + angle.sin;
  }
}

//! computes (a+b), if (a+b) >= 2*PI, the result is (a+b-2*PI)
//! a and be have to be in the range [0;2*PI[
static inline struct sin_cos_angle sum_angles_no_check(
  struct sin_cos_angle a, struct sin_cos_angle b) {

  return sin_cos_angle_new(a.sin * b.cos + a.cos * b.sin,
                           a.cos * b.cos - a.sin * b.sin);
}

//! computes (a+b), if (a+b) >= 2*PI, the result is (a+b-2*PI)
//! a and be have to be in the range [0;2*PI[
//! \returns 1 if (a+b) is >= 2*PI
static inline int sum_angles(
  struct sin_cos_angle a, struct sin_cos_angle b,
  struct sin_cos_angle * restrict sum) {

  struct sin_cos_angle sum_ = sum_angles_no_check(a, b);

  *sum = sum_;

  // if a or b is smaller than the result
  return (compare_angles(sum_, a) < 0) || (compare_angles(sum_, b) < 0);
}

//! computes (a-b), if (a-b) < 0, the result is (a-b+2*PI)
//! a and be have to be in the range [0;2*PI[
static inline struct sin_cos_angle sub_angles_no_check(
  struct sin_cos_angle a, struct sin_cos_angle b) {

  return sin_cos_angle_new(a.sin * b.cos - a.cos * b.sin,
                           a.cos * b.cos + a.sin * b.sin);
}

//! computes (a-b), if (a-b) <= 0, the result is (a-b+2*PI)
//! a and be have to be in the range [0;2*PI[
//! \returns 1 if (a-b) < 0
static inline int sub_angles(
  struct sin_cos_angle a, struct sin_cos_angle b,
  struct sin_cos_angle * restrict sub) {

  int compare_result = compare_angles(a, b);

  // return sin=0.0 and cos=1.0 if the angles are equal to each other,
  // i.e. compare_result == 0
  *sub = compare_result ? sub_angles_no_check(a, b) : SIN_COS_ZERO;

  return compare_result < 0;
}

//! return angles in the range of [0;2PI[
static inline double compute_angle(struct sin_cos_angle angle) {

  // the acos and asin are most accurate in the range [-M_SQRT1_2;M_SQRT1_2]
  if (angle.cos > M_SQRT1_2) {

    double angle_ = asin(angle.sin);

    if (angle_ < 0.0) angle_ += 2.0 * M_PI;
    return angle_;

  } else if (angle.cos > - M_SQRT1_2) {

    double angle_ = acos(angle.cos);

    if (angle.sin > 0.0) return angle_;
    else return 2.0 * M_PI - angle_;

  } else {
    return M_PI - asin(angle.sin);
  }
}

//! computes angle / 2
//! The basic idea is to imagin angle being a vector v with (cos(x); sin(x)) and
//! w being a vector with (1.0; 0.0). The normalised sum of v and w gives us v_h
//! with (cos(x/2); sin(x/2)).
//! For quadrants 2, 3, and 4 we have to apply some little tricks for higher
//! accuracy.
static inline struct sin_cos_angle half_angle(struct sin_cos_angle angle) {

  double x = (1.0 + fabs(angle.cos));

  double scale = 1.0 / sqrt(x * x + angle.sin * angle.sin);

  // first or fourth quadrant
  if (angle.cos >= 0) {
    scale = copysign(scale, angle.sin);
    return sin_cos_angle_new(angle.sin * scale, x * scale);

  // second and third quadrant
  } else return sin_cos_angle_new(x * scale, angle.sin * scale);
}

//! computes angle / 4
//! I derived it based on \ref half_angle
static inline struct sin_cos_angle quarter_angle(struct sin_cos_angle angle) {

  double tan = fabs(angle.sin / (1.0 + fabs(angle.cos)));
  double one_plus_sq_tan = 1.0 + tan * tan;
  double sqrt_one_plus_sq_tan = sqrt(one_plus_sq_tan);

  double a = 1.0;
  double b = tan;

  // second and third quadrant
  if (angle.cos < 0.0) {
    a = tan;
    b = 1.0;
  }

  double scale = M_SQRT1_2 / sqrt(one_plus_sq_tan + a * sqrt_one_plus_sq_tan);
  double x = b * scale;
  double y = (a + sqrt_one_plus_sq_tan) * scale;

  // first and second quadrant
  if (angle.sin >= 0.0) return sin_cos_angle_new(x, y);
  // third and fourth quadrant
  else return sin_cos_angle_new(y, x);
}

/**
 * determines whether two given points are (nearly) identically
 * @param[in] a point coordinates of point a
 * @param[in] b point coordinates of point b
 * @return true if both points are (nearly) identically
 */
static inline int points_are_identically(double const * a, double const * b) {

  // for small angles the Euclidean distance is nearly identical to
  // great circle distance between vectors a and b
  return sq_len_diff_vec(a, b) <= yac_sq_angle_tol;
}

/**
 * compares two 3d coordinates (vectors need to have a length of 1.0)
 * @param[in] a coordinates of point a
 * @param[in] b coordinates of point b
 * @return 0 if points are identical, otherwise -1 or 1 depending on their
 *           relation
 */
static inline int compare_coords(double const * a, double const * b) {

   double dot_product = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

   // if the angle is smaller than ~0.81 degree
   // (in this range the acos is still rather accurate)
   if (dot_product > 0.9999) { // (acos(0.9999) = ~0.81 degree)

      // both points are close to each other -> use cross product for higher
     // accuracy

      double cross_ab[3];
      crossproduct_kahan(a, b, cross_ab);

      // for very small angles: asin(alpha) = ~alpha   (alpha in rad)
      if (sqrt(cross_ab[0]*cross_ab[0] +
               cross_ab[1]*cross_ab[1] +
               cross_ab[2]*cross_ab[2]) < yac_angle_tol) return 0;
   }

   int ret;
   if ((ret = (a[0] > b[0]) - (a[0] < b[0]))) return ret;
   if ((ret = (a[1] > b[1]) - (a[1] < b[1]))) return ret;
   return (a[2] > b[2]) - (a[2] < b[2]);
}

/** normalises a vector */
static inline void normalise_vector(double v[]) {

   double norm = 1.0 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

   v[0] *= norm;
   v[1] *= norm;
   v[2] *= norm;
}

/** rotate vector v_in around the given axis by a given angle
 * @param[in] axis axis around which v_in is to rotated
 * @param[in] angle rotation angle
 * @param[in] v_in vector to be rotated
 * @param[out] v_out rotated vector
 */
static inline void rotate_vector2(
  double axis[], struct sin_cos_angle angle, double v_in[], double v_out[]) {

  // using Rodrigues' rotation formula
  // v_out = v_in * cos(angle) +
  //         (axis x v_in) * sin(angle) +
  //         axis * (axis * v_in) * (1 - cos(angle))

  double cross_axis_v_in[3];
  crossproduct_d(axis, v_in, cross_axis_v_in);

  double dot_axis_v_in = axis[0]*v_in[0] + axis[1]*v_in[1] + axis[2]*v_in[2];
  double temp = dot_axis_v_in * (1.0 - angle.cos);

  v_out[0] =
    v_in[0] * angle.cos + cross_axis_v_in[0] * angle.sin + axis[0] * temp;
  v_out[1] =
    v_in[1] * angle.cos + cross_axis_v_in[1] * angle.sin + axis[1] * temp;
  v_out[2] =
    v_in[2] * angle.cos + cross_axis_v_in[2] * angle.sin + axis[2] * temp;
}

/** rotate vector v_in around the given axis by a given angle
 * @param[in] axis axis around which v_in is to rotated
 * @param[in] angle rotation angle
 * @param[in] v_in vector to be rotated
 * @param[out] v_out rotated vector
 */
static inline void rotate_vector(
  double axis[], double angle, double v_in[], double v_out[]) {

  double sin_angle = sin(angle);
  double cos_angle = cos(angle);
  rotate_vector2(
    axis, sin_cos_angle_new(sin_angle, cos_angle), v_in, v_out);
}

/**
 * splits given cell into triangles
 * @param[in] cell cell to be triangulated
 * @param[in] start_corner start algorithm at corner with this index
 *                         (0 <= start_corner < n; with n being number of cell
 *                          corners)
 * @param[out] triangles triangles that are the result of the triangulation
 * @remark the user needs to provide (n-2) initialised grid_cells for the
 *         argument triangles (n being the number of cell corners of the cell)
 * @remark this routine currently only works for convex grid cells
 * @remark the algorithm ignores that actual edge types for cells with more than
 *         3 corners and assumes them to be great circle edges
 */
void yac_triangulate_cell(
  struct yac_grid_cell cell, size_t start_corner, struct yac_grid_cell * triangles);

/**
 * splits given indices of a cell into triangles
 * @param[in] corner_indices cell corner indices to be triangulated
 * @param[in] num_corners    number of corners of the cell
 * @param[in] start_corner   start algorithm at corner with this index
 *                           (0 <= start_corner < n; with n being number of cell
 *                            corners)
 * @param[out] triangle_indices triangle indices that are the result of the
 *                              triangulation
 * @remark the user needs to provide provide an array of size 3*(num_corners-2)
 *         for the argument triangle_indices
 */
void yac_triangulate_cell_indices(
  size_t const * corner_indices, size_t num_corners, size_t start_corner,
  size_t (*triangle_indices)[3]);

/**
 * computes a spherical triangle that contains all given vertices,
 * the algorithm tries to minimise the size of the triangle
 * @param[in] vertices vertices that are supposed to be within the triangle
 * @param[in] num_vertices number of vertices in vertices
 * @param[out] triangle bounding triangle
 * @param[in] num_tests with increasing number of tests computation time and
 *                      chance of smaller bounding triangle increases
 * @remark num_tests has to be > 0
 */
void yac_compute_bnd_triangle(
  double * vertices, size_t num_vertices, double triangle[][3],
  size_t num_tests);

int yac_circle_intersect(
  struct yac_circle a, struct yac_circle b, double p[3], double q[3]);

int yac_point_on_edge(
  double p[3], double const a[3], double const b[3],
  enum yac_circle_type circle_type);

#endif // GEOMETRY_H

