/**
 * @file test_angle.c
 *
 * @copyright Copyright  (C)  2018 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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
#include "geometry.h"
#include "tests.h"

#define NUM_TESTS (1000)

static void generate_rand_angle(
  double * angle, struct sin_cos_angle * sin_cos_angle);

static struct sin_cos_angle compute_sin_cos_angle(double angle);
static int compare_double(double a, double b);
static int compare_size_t(size_t a, size_t b);
static void check_angle_sum_sub(
  struct sin_cos_angle a, struct sin_cos_angle b,
  double ref_sum, double ref_sub);
static void check_compare_angles(double dble_angle_a, double dble_angle_b,
                                 struct sin_cos_angle sin_cos_angle_a,
                                 struct sin_cos_angle sin_cos_angle_b);

int main (void) {

  srand(1337);

  struct sin_cos_angle sin_cos_angle;
  double angle;

  // multiples of PI
  double angles[] = {0.0, 0.24, 0.26, 0.5, 0.8, 1.0, 1.1, 1.5, 1.6, 1.9};
  for (size_t i = 0; i < sizeof(angles) / sizeof(angles[0]); ++i)
    angles[i] *= M_PI;

  struct sin_cos_angle sin_cos_angles[] =
    {(struct sin_cos_angle){.sin = 0.0, .cos = 1.0},               // 0PI
     (struct sin_cos_angle){.sin = yac_angle_tol, .cos = 1.0},     // 0PI + tol
     (struct sin_cos_angle){.sin = 0.5, .cos = M_SQRT3_4},         // PI/6
     (struct sin_cos_angle){.sin =  M_SQRT1_2-yac_angle_tol,
                            .cos =  M_SQRT1_2+yac_angle_tol},      // PI/4 - tol
     (struct sin_cos_angle){.sin =  M_SQRT1_2, .cos = M_SQRT1_2},  // PI/4
     (struct sin_cos_angle){.sin =  M_SQRT1_2+yac_angle_tol,
                            .cos =  M_SQRT1_2-yac_angle_tol},      // PI/4 + tol
     (struct sin_cos_angle){.sin =  M_SQRT3_4, .cos = 0.5},        // PI/3
     (struct sin_cos_angle){.sin = 1.0, .cos = yac_angle_tol},     // PI/2 - tol
     (struct sin_cos_angle){.sin = 1.0, .cos = 0.0},               // PI/2
     (struct sin_cos_angle){.sin = 1.0, .cos = - yac_angle_tol},   // PI/2 + tol
     (struct sin_cos_angle){.sin =  M_SQRT3_4, .cos = -0.5},       // 2PI/3
     (struct sin_cos_angle){.sin =  M_SQRT1_2+yac_angle_tol,
                            .cos = -M_SQRT1_2+yac_angle_tol},      // 3PI/4 - tol
     (struct sin_cos_angle){.sin =  M_SQRT1_2, .cos = -M_SQRT1_2}, // 3PI/4
     (struct sin_cos_angle){.sin =  M_SQRT1_2-yac_angle_tol,
                            .cos = -M_SQRT1_2-yac_angle_tol},      // 3PI/4 + tol
     (struct sin_cos_angle){.sin = 0.5, .cos = -M_SQRT3_4},        // 5PI/6
     (struct sin_cos_angle){.sin = yac_angle_tol, .cos = -1.0},    // PI/2 - tol
     (struct sin_cos_angle){.sin = 0.0, .cos = -1.0},              // PI/2
     (struct sin_cos_angle){.sin = -yac_angle_tol, .cos = -1.0},   // PI/2 + tol
     (struct sin_cos_angle){.sin = -0.5, .cos = -M_SQRT3_4},       // 7PI/6
     (struct sin_cos_angle){.sin = -M_SQRT1_2+yac_angle_tol,
                            .cos = -M_SQRT1_2-yac_angle_tol},      // 5PI/4 - tol
     (struct sin_cos_angle){.sin = -M_SQRT1_2, .cos = -M_SQRT1_2}, // 5PI/4
     (struct sin_cos_angle){.sin = -M_SQRT1_2-yac_angle_tol,
                            .cos = -M_SQRT1_2+yac_angle_tol},      // 5PI/4 + tol
     (struct sin_cos_angle){.sin = -M_SQRT3_4, .cos = -0.5},       // 4PI/3
     (struct sin_cos_angle){.sin = -1.0, .cos = -yac_angle_tol},   // 3PI/2 - tol
     (struct sin_cos_angle){.sin = -1.0, .cos = 0.0},              // 3PI/2
     (struct sin_cos_angle){.sin = -1.0, .cos = yac_angle_tol},    // 3PI/2 + tol

     (struct sin_cos_angle){.sin = -M_SQRT3_4, .cos = 0.5},        // 5PI/3
     (struct sin_cos_angle){.sin = -M_SQRT1_2-yac_angle_tol,
                            .cos =  M_SQRT1_2-yac_angle_tol},      // 7PI/4 - tol
     (struct sin_cos_angle){.sin = -M_SQRT1_2, .cos = M_SQRT1_2},  // 7PI/4
     (struct sin_cos_angle){.sin = -M_SQRT1_2+yac_angle_tol,
                            .cos =  M_SQRT1_2+yac_angle_tol},      // 7PI/4 + tol
     (struct sin_cos_angle){.sin = -0.5, .cos = M_SQRT3_4},        // 11PI/6
     (struct sin_cos_angle){.sin = -yac_angle_tol, .cos = 1.0}     // 2PI - tol
     };

  { // test routine compute_angle
    for (size_t i = 0; i < sizeof(angles) / sizeof(angles[0]); ++i)
      if (fabs(compute_angle(compute_sin_cos_angle(angles[i])) - angles[i]) >
          yac_angle_tol)
        PUT_ERR("error in compute_angle");

    for (size_t i = 0; i < NUM_TESTS; ++i) {
      double temp_angle;
      generate_rand_angle(&temp_angle, &sin_cos_angle);
      if (fabs(temp_angle - compute_angle(sin_cos_angle)) > yac_angle_tol)
        PUT_ERR("error in compute_angle");
    }
  }

  { // test routine compare_angles and sin_cos_angle_to_dble
    for (size_t i = 0; i < sizeof(angles)/sizeof(angles[0]); ++i) {
      for (size_t j = 0; j < sizeof(angles)/sizeof(angles[0]); ++j) {
        struct sin_cos_angle temp_sin_cos_angles[2] =
          {compute_sin_cos_angle(angles[i]), compute_sin_cos_angle(angles[j])};
        check_compare_angles(
          angles[i], angles[j], temp_sin_cos_angles[0], temp_sin_cos_angles[1]);
      }
    }

    for (size_t i = 0; i < NUM_TESTS; ++i) {
      double temp_angles[2];
      struct sin_cos_angle temp_sin_cos_angles[2];
      for (size_t j = 0; j < 2; ++j)
        generate_rand_angle(&(temp_angles[j]), &(temp_sin_cos_angles[j]));
      for (size_t j = 0; j < 2; ++j)
        for (size_t k = 0; k < 2; ++k)
          check_compare_angles(temp_angles[j], temp_angles[k],
                               temp_sin_cos_angles[j], temp_sin_cos_angles[k]);
    }

    for (size_t i = 0; i < sizeof(sin_cos_angles)/sizeof(sin_cos_angles[0]); ++i) {
      for (size_t j = 0; j < sizeof(sin_cos_angles)/sizeof(sin_cos_angles[0]); ++j) {
        double temp[2] = {sin_cos_angle_to_dble(sin_cos_angles[i]),
                          sin_cos_angle_to_dble(sin_cos_angles[j])};
        int compare_ij = compare_size_t(i, j);
        if ((compare_ij != compare_angles(sin_cos_angles[i], sin_cos_angles[j])) ||
            (compare_ij != compare_double(temp[0], temp[1])))
          PUT_ERR("error in compare_angles");
      }
    }
  }

  { // test routine sum_angles, sub_angles
    for (size_t i = 0; i < sizeof(angles)/sizeof(angles[0]); ++i) {
      for (size_t j = 0; j < sizeof(angles)/sizeof(angles[0]); ++j) {
        struct sin_cos_angle temp_sin_cos_angles[2] =
        {compute_sin_cos_angle(angles[i]), compute_sin_cos_angle(angles[j])};
        check_angle_sum_sub(temp_sin_cos_angles[0], temp_sin_cos_angles[1],
                            angles[i] + angles[j], angles[i] - angles[j]);
      }
    }

    for (size_t i = 0; i < NUM_TESTS; ++i) {
      double temp_angles[2];
      struct sin_cos_angle temp_sin_cos_angles[2];
      for (size_t j = 0; j < 2; ++j)
        generate_rand_angle(&(temp_angles[j]), &(temp_sin_cos_angles[j]));
      for (size_t j = 0; j < 2; ++j)
        for (size_t k = 0; k < 2; ++k)
          check_angle_sum_sub(
            temp_sin_cos_angles[j], temp_sin_cos_angles[k],
            temp_angles[j] + temp_angles[k], temp_angles[j] - temp_angles[k]);
    }
  }

  { // test routine half_angle
    for (size_t i = 0; i < sizeof(angles)/sizeof(angles[0]); ++i)
      if (fabs(angles[i] * 0.5 -
               compute_angle(half_angle(compute_sin_cos_angle(angles[i])))) >
          yac_angle_tol)
        PUT_ERR("error in half_angle");

    for (size_t i = 0; i < NUM_TESTS; ++i) {
      generate_rand_angle(&angle, &sin_cos_angle);
      if (fabs(angle * 0.5 -
               compute_angle(half_angle(sin_cos_angle))) > yac_angle_tol)
        PUT_ERR("error in half_angle");
    }
  }

  { // test routine quarter_angle
    for (size_t i = 0; i < sizeof(angles)/sizeof(angles[0]); ++i)
      if (fabs(angles[i] * 0.25 -
               compute_angle(quarter_angle(compute_sin_cos_angle(angles[i])))) >
          yac_angle_tol)
        PUT_ERR("error in quarter_angle");

    if (fabs(M_PI * 0.25 - compute_angle(quarter_angle(SIN_COS_M_PI))) >
        yac_angle_tol) PUT_ERR("error in quarter_angle");

    for (size_t i = 0; i < NUM_TESTS; ++i) {
      generate_rand_angle(&angle, &sin_cos_angle);
      if (fabs(angle * 0.25 -
               compute_angle(quarter_angle(sin_cos_angle))) > yac_angle_tol)
        PUT_ERR("error in quarter_angle");
    }
  }

  return TEST_EXIT_CODE;
}

static void generate_rand_angle(
  double * angle, struct sin_cos_angle * sin_cos_angle) {

  *angle = (2.0 * M_PI) * ((double)rand() / (double)RAND_MAX);
  *sin_cos_angle = compute_sin_cos_angle(*angle);
}

static struct sin_cos_angle compute_sin_cos_angle(double angle) {

  double sin_angle = sin(angle);
  double cos_angle = cos(angle);
  return sin_cos_angle_new(sin_angle, cos_angle);
}

static int compare_double(double a, double b) {

  return (a > b) - (a < b);
}

static int compare_size_t(size_t a, size_t b) {

  return (a > b) - (a < b);
}

static void check_angle_sum_sub(
  struct sin_cos_angle a, struct sin_cos_angle b,
  double ref_sum, double ref_sub) {

  struct sin_cos_angle sum;
  int big_sum = sum_angles(a, b, &sum);

  double angle_sum = compute_angle(sum) + ((big_sum)?(2.0*M_PI):0.0);

  if (fabs(angle_sum - ref_sum) > yac_angle_tol) PUT_ERR("error sum_angles");

  struct sin_cos_angle sub;
  int neg = sub_angles(a, b, &sub);

  double angle_sub = compute_angle(sub) - ((neg)?(2.0*M_PI):0.0);

  if (fabs(angle_sub - ref_sub) > yac_angle_tol) PUT_ERR("error sub_angles");
}

static void check_compare_angles(double dble_angle_a, double dble_angle_b,
                                 struct sin_cos_angle sin_cos_angle_a,
                                 struct sin_cos_angle sin_cos_angle_b) {

  int ret = compare_double(dble_angle_a, dble_angle_b);
  if (ret != compare_angles(sin_cos_angle_a, sin_cos_angle_b))
    PUT_ERR("error in compare_angles");
  double temp[2] = {sin_cos_angle_to_dble(sin_cos_angle_a),
                    sin_cos_angle_to_dble(sin_cos_angle_b)};
  if (ret != compare_double(temp[0], temp[1]))
    PUT_ERR("error in sin_cos_angle_to_dble");
}

