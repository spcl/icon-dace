// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TEST_CXC_H
#define TEST_CXC_H

#include "geometry.h"
#include "tests.h"

/** \example test_cxc.c
 * A test for cell intersections.
 */

void test_cxc(double lon_a, double lat_a, double lon_b, double lat_b,
              double lon_c, double lat_c, double lon_d, double lat_d,
              enum yac_edge_type edge_type_a, enum yac_edge_type edge_type_b,
              double lon_ref_p, double lat_ref_p,
              double lon_ref_q, double lat_ref_q, int ref_ret_val);

void test_cxc_rad(double lon_a, double lat_a, double lon_b, double lat_b,
                  double lon_c, double lat_c, double lon_d, double lat_d,
                  enum yac_edge_type edge_type_a, enum yac_edge_type edge_type_b,
                  double lon_ref_p, double lat_ref_p,
                  double lon_ref_q, double lat_ref_q, int ref_ret_val);
void get_edge_middle_point(
  enum yac_edge_type edge_type,
  double lon_a, double lat_a, double lon_b, double lat_b,
  double * lon_middle, double * lat_middle);

#endif // TEST_CXC_H

