// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <math.h>

#include "geometry.h"
#include "test_cxc.h"

enum {
  error = -1,
  p_on_a = 1,
  q_on_a = 2,
  p_on_b = 4,
  q_on_b = 8,
  circles_are_identical = 16,
};

static void test_results_vec(double p[3], double q[3], int ret_val,
                             double ref_p[3], double ref_q[3],
                             int ref_ret_value) {

  int const mask_p = p_on_a + p_on_b;
  int const mask_q = q_on_a + q_on_b;

  // if intersection could be computed or no intersection is possible
  if (ref_ret_value == error) {

    if (ret_val != error)
      PUT_ERR("error in ret_val (no error indicated)\n")

  } else {

    if ((ret_val & circles_are_identical) ^
        (ref_ret_value & circles_are_identical))
      PUT_ERR("error in ret_val (identically circles)\n")

    // if there is only one intersection point
    if (points_are_identically(ref_q, ref_p)) {

      if (ret_val != ref_ret_value)
          PUT_ERR("error in ret_val (a)");

      if (!points_are_identically(p, ref_p))
          PUT_ERR("p does not match\n")

      if (!points_are_identically(q, ref_p))
          PUT_ERR("q does not match\n")

    } else {

      if (points_are_identically(p, ref_p)) {

        if ((ret_val & mask_p) != (ref_ret_value & mask_p))
          PUT_ERR("error in ret_val (b)\n")

        if (!points_are_identically(q, ref_q))
          PUT_ERR("q does not match\n")

        if ((ret_val & mask_q) != (ref_ret_value & mask_q))
          PUT_ERR("error in ret_val (c)\n")

      } else {

        if (!points_are_identically(p, ref_q))
          PUT_ERR("p does not match\n")

        if ((ret_val & mask_p) << 1 != (ref_ret_value & mask_q))
          PUT_ERR("error in ret_val (d)\n")

        if (!points_are_identically(q, ref_p))
          PUT_ERR("q does not match\n")

        if ((ret_val & mask_q) >> 1 != (ref_ret_value & mask_p))
          PUT_ERR("error in ret_val (e)\n")
      }

      if (ref_ret_value == 0 && ret_val != 0)
        PUT_ERR("error in return value (f)\n")
    }
  }
}

void test_cxc_rad(double lon_a, double lat_a, double lon_b, double lat_b,
                  double lon_c, double lat_c, double lon_d, double lat_d,
                  enum yac_edge_type edge_type_a, enum yac_edge_type edge_type_b,
                  double lon_ref_p, double lat_ref_p,
                  double lon_ref_q, double lat_ref_q, int ref_ret_val) {

  double edges[2][2][3];
  LLtoXYZ(lon_a, lat_a, edges[0][0]);
  LLtoXYZ(lon_b, lat_b, edges[0][1]);
  LLtoXYZ(lon_c, lat_c, edges[1][0]);
  LLtoXYZ(lon_d, lat_d, edges[1][1]);

  enum yac_edge_type edge_types[2] = {edge_type_a, edge_type_b};

  double ref_pq[2][3];
  LLtoXYZ(lon_ref_p, lat_ref_p, ref_pq[0]);
  LLtoXYZ(lon_ref_q, lat_ref_q, ref_pq[1]);

  int ref_ret_vals[2] =
    {ref_ret_val,
     // remove p/q bits
     (ref_ret_val & (~(1+2+4+8))) |
     // move p_on_a bit to q_on_b
     ((ref_ret_val & p_on_a)?(q_on_b):(0)) |
     // move q_on_a bit to p_on_b
     ((ref_ret_val & q_on_a)?(p_on_b):(0)) |
     // move p_on_b bit to q_on_a
     ((ref_ret_val & p_on_b)?(q_on_a):(0)) |
     // move q_on_b bit to p_on_a
     ((ref_ret_val & q_on_b)?(p_on_a):(0))};

  // edge order (results have to be independent of the order of both edge)
  for (int i = 0; i < 2; ++i) {
    // vertex order of first edge
    for (int j = 0; j < 2; ++j) {
      // vertex order of second edge
      for (int k = 0; k < 2; ++k) {

        double p[3], q[3];
        int ret_val =
          yac_intersect_vec(
            edge_types[i], edges[i][j], edges[i][j^1],
            edge_types[i^1], edges[i^1][k], edges[i^1][k^1],
            p, q);
        test_results_vec(
          p, q, ret_val, ref_pq[i], ref_pq[i^1], ref_ret_vals[i]);
      }
    }
  }
}

void test_cxc(double lon_a, double lat_a, double lon_b, double lat_b,
              double lon_c, double lat_c, double lon_d, double lat_d,
              enum yac_edge_type edge_type_a, enum yac_edge_type edge_type_b,
              double lon_ref_p, double lat_ref_p,
              double lon_ref_q, double lat_ref_q, int ref_ret_val) {

   test_cxc_rad(lon_a*YAC_RAD, lat_a*YAC_RAD, lon_b*YAC_RAD, lat_b*YAC_RAD,
                lon_c*YAC_RAD, lat_c*YAC_RAD, lon_d*YAC_RAD, lat_d*YAC_RAD,
                edge_type_a, edge_type_b, lon_ref_p*YAC_RAD, lat_ref_p*YAC_RAD,
                lon_ref_q*YAC_RAD, lat_ref_q*YAC_RAD, ref_ret_val);
}

static double adjust_lon(double lon, double ref) {
  while (lon < ref - 180) lon += 360;
  while (lon > ref + 180) lon -= 360;
  return lon;
}

void get_edge_middle_point(
  enum yac_edge_type edge_type,
  double lon_a, double lat_a, double lon_b, double lat_b,
  double * lon_middle, double * lat_middle) {

  double a[3], b[3];
  LLtoXYZ_deg(lon_a, lat_a, a);
  LLtoXYZ_deg(lon_b, lat_b, b);

  YAC_ASSERT(
    (edge_type == YAC_GREAT_CIRCLE_EDGE) ||
    (edge_type == YAC_LON_CIRCLE_EDGE) ||
    (edge_type == YAC_LAT_CIRCLE_EDGE),
    "ERROR(get_edge_middle_point): invalid edge_type")

  switch (edge_type) {
    default:
    case (YAC_LON_CIRCLE_EDGE):
    case (YAC_GREAT_CIRCLE_EDGE): {
      double middle[3];
      middle[0] = a[0] + b[0];
      middle[1] = a[1] + b[1];
      middle[2] = a[2] + b[2];

      normalise_vector(middle);

      XYZtoLL(middle, lon_middle, lat_middle);
      *lon_middle /= YAC_RAD;
      *lat_middle /= YAC_RAD;
      break;
    }
    case (YAC_LAT_CIRCLE_EDGE): {
      *lon_middle = (lon_a + adjust_lon(lon_b, lon_a))*0.5;
      *lat_middle = lat_a;
    }
  }
}
