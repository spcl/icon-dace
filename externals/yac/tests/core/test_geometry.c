// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tests.h"

#include "geometry.h"

int main (void) {

  {
    double base_point[2] = {2.253 * YAC_RAD, 3.25 * YAC_RAD};
    double point_a[2] = {2 * YAC_RAD, 2 * YAC_RAD};
    double point_b[2] = {1 * YAC_RAD, 3 * YAC_RAD};

    double base_vector[3], vector_a[3], vector_b[3];

    LLtoXYZ(base_point[0], base_point[1], base_vector);
    LLtoXYZ(point_a[0], point_a[1], vector_a);
    LLtoXYZ(point_b[0], point_b[1], vector_b);

    double distance_a, distance_b;

    distance_a = get_vector_angle(base_vector, vector_a);
    distance_b = get_vector_angle(base_vector, vector_b);

    if (distance_a >= distance_b)
      PUT_ERR("error in distance calculation\n");
  }

  {
    double const tol = 1e-10;
    double angles[] = {-M_PI_2 + tol, -M_PI_4, -tol, 0,
                         tol, M_PI_4, M_PI_2 - tol};
    double offsets[] = {-6.0*M_PI, -4.0*M_PI, -2.0*M_PI, 0,
                         2.0*M_PI, 4.0*M_PI, 6.0*M_PI};
    double start[] = {-M_PI, -M_PI_2, -M_PI_4, -M_1_PI, -tol, 0,
                       tol, M_1_PI, M_PI_4, M_PI_2, M_PI};

    for (size_t i = 0; i < sizeof(angles) / sizeof(angles[0]); ++i) {
      for (size_t j = 0; j < sizeof(offsets) / sizeof(offsets[0]); ++j) {
        for (size_t k = 0; k < sizeof(offsets) / sizeof(offsets[0]); ++k) {
          for (size_t l = 0; l < sizeof(start) / sizeof(start[0]); ++l) {

            double a_lon = start[l] + offsets[j] + angles[i];
            double b_lon = start[l] + offsets[k];

            if (fabs(get_angle(a_lon, b_lon) - angles[i]) > tol)
              PUT_ERR("error in get_angle\n");
          }
        }
      }
    }
  }

  {
    struct bounding_circle a =
      {.base_vector = {0.0,0.0,1.0},
       .inc_angle   = SIN_COS_M_PI,
       .sq_crd      = DBL_MAX};
    struct bounding_circle b =
      {.base_vector = {0.0,0.0,-1.0},
       .inc_angle   = SIN_COS_M_PI,
       .sq_crd      = DBL_MAX};

    if (!yac_extents_overlap(&a, &b))
      PUT_ERR("error in yac_extents_overlap\n");
  }

  {
    struct bounding_circle a =
      {.base_vector = {0.0,0.0,1.0},
       .inc_angle   = SIN_COS_M_PI,
       .sq_crd      = DBL_MAX};
    struct bounding_circle b =
      {.base_vector = {0.0,0.0,-1.0},
       .inc_angle   = SIN_COS_M_PI,
       .sq_crd      = DBL_MAX};
    sub_angles(SIN_COS_M_PI, SIN_COS_LOW_TOL, &b.inc_angle);

    if (!yac_extents_overlap(&a, &b))
      PUT_ERR("error in yac_extents_overlap\n");
  }
  return TEST_EXIT_CODE;
}

