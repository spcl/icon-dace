/**
 * @file test_bnd_sphere_part.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *
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
#include <float.h>

#include "tests.h"
#include "test_common.h"
#include "grid.h"
#include "sphere_part.h"
#include "geometry.h"

static void random_test(size_t count);
static void random_test_big(size_t count);
static void empty_test();
static void zero_balance_point();
static void mean_point_test();

int main(void) {

  srand(1365441);

  for (size_t i = 0; i < 10; ++i) {
    random_test(1);
    random_test(10);
    random_test(100);
    random_test(650);
    random_test_big(50);
  }
  empty_test();
  zero_balance_point();
  mean_point_test();

  return TEST_EXIT_CODE;
}

static void gen_random_point(double * point) {

  double x_coordinate = 2.0 * M_PI * (((double)rand()) / ((double)RAND_MAX));
  double y_coordinate = M_PI * (((double)rand()) / ((double)RAND_MAX)) - M_PI_2;
  LLtoXYZ(x_coordinate, y_coordinate, point);
}

static void gen_random_bnd_circle(struct bounding_circle * circle) {

  gen_random_point(circle->base_vector);
  double angle = M_PI_2 * (((double)rand()) / ((double)RAND_MAX));
  circle->inc_angle = sin_cos_angle_new(sin(angle), cos(angle));
  circle->sq_crd = DBL_MAX;
}

static void gen_random_bnd_circle_big(struct bounding_circle * circle) {

  gen_random_point(circle->base_vector);
  double angle = 3.0 * M_PI_2 * (((double)rand()) / ((double)RAND_MAX));
  circle->inc_angle = sin_cos_angle_new(sin(angle), cos(angle));
  circle->sq_crd = DBL_MAX;
}

static int compare_size_t(void const * a, void const * b) {
  return
    (*(size_t const *)a > *(size_t const *)b) -
    (*(size_t const *)a < *(size_t const *)b);
}

static int check_results(
  size_t * list, size_t list_size, size_t * sorted_indices, size_t count) {

  qsort(list, list_size, sizeof(size_t), compare_size_t);

  size_t match_count = 0;
  for (size_t i = 0, j = 0; i < count; ++i) {
    while ((j < list_size) && (list[j] < sorted_indices[i])) ++j;
    if ((j < list_size) && (list[j] == sorted_indices[i])) ++match_count;
  }
  return match_count != count;
}

static void random_test(size_t count) {

  struct bounding_circle * bnd_circles =
    xmalloc(2 * count * sizeof(*bnd_circles));
  struct bounding_circle * bnd_circles_a = bnd_circles;
  struct bounding_circle * bnd_circles_b = bnd_circles + count;
  double (*point_coords)[3] = xmalloc(count * sizeof(*point_coords));

  for (size_t i = 0; i < count; ++i) {

    gen_random_bnd_circle(bnd_circles_a + i);
    gen_random_bnd_circle(bnd_circles_b + i);
    gen_random_point(point_coords[i]);
  }

  struct bnd_sphere_part_search * search =
    yac_bnd_sphere_part_search_new(bnd_circles_a, count);

  { // check yac_bnd_sphere_part_search_do_point_search
    size_t * results;
    size_t * num_results = xmalloc(count * sizeof(*num_results));
    yac_bnd_sphere_part_search_do_point_search(
      search, point_coords, count, &results, num_results);

    size_t * ref_results = xmalloc(count * sizeof(*ref_results));
    for (size_t i = 0, offset = 0; i < count; ++i) {
      size_t ref_num_results = 0;
      for (size_t j = 0; j < count; ++j)
        if (yac_point_in_bounding_circle_vec(
              point_coords[i], bnd_circles_a + j))
          ref_results[ref_num_results++] = j;

      if ((ref_num_results > num_results[i]) ||
          check_results(
            results + offset, num_results[i], ref_results, ref_num_results))
        PUT_ERR("ERROR in yac_bnd_sphere_part_search_do_point_search");

      offset += num_results[i];
    }
    free(ref_results);
    free(results);
    free(num_results);
  }

  { // yac_bnd_sphere_part_search_do_bnd_circle_search
    size_t * results;
    size_t * num_results = xmalloc(count * sizeof(*num_results));
    yac_bnd_sphere_part_search_do_bnd_circle_search(
      search, bnd_circles_b, count, &results, num_results);

    size_t * ref_results = xmalloc(count * sizeof(*ref_results));
    for (size_t i = 0, offset = 0; i < count; ++i) {
      size_t ref_num_results = 0;
      for (size_t j = 0; j < count; ++j)
        if (yac_extents_overlap(
              bnd_circles_b + i, bnd_circles_a + j))
          ref_results[ref_num_results++] = j;

      if ((ref_num_results > num_results[i]) ||
          check_results(
            results + offset, num_results[i], ref_results, ref_num_results))
        PUT_ERR("ERROR in yac_bnd_sphere_part_search_do_bnd_circle_search");

      offset += num_results[i];
    }
    free(ref_results);
    free(results);
    free(num_results);
  }

  yac_bnd_sphere_part_search_delete(search);
  free(point_coords);
  free(bnd_circles);
}

static void random_test_big(size_t count) {

  struct bounding_circle * bnd_circles =
    xmalloc(2 * count * sizeof(*bnd_circles));
  struct bounding_circle * bnd_circles_a = bnd_circles;
  struct bounding_circle * bnd_circles_b = bnd_circles + count;
  double (*point_coords)[3] = xmalloc(count * sizeof(*point_coords));

  for (size_t i = 0; i < count; ++i) {

    gen_random_bnd_circle(bnd_circles_a + i);
    gen_random_bnd_circle_big(bnd_circles_b + i);
    gen_random_point(point_coords[i]);
  }

  struct bnd_sphere_part_search * search =
    yac_bnd_sphere_part_search_new(bnd_circles_a, count);

  { // yac_bnd_sphere_part_search_do_bnd_circle_search
    size_t * results;
    size_t * num_results = xmalloc(count * sizeof(*num_results));
    yac_bnd_sphere_part_search_do_bnd_circle_search(
      search, bnd_circles_b, count, &results, num_results);

    size_t * ref_results = xmalloc(count * sizeof(*ref_results));
    for (size_t i = 0, offset = 0; i < count; ++i) {
      size_t ref_num_results = 0;
      for (size_t j = 0; j < count; ++j)
        if (yac_extents_overlap(
              bnd_circles_b + i, bnd_circles_a + j))
          ref_results[ref_num_results++] = j;

      if ((ref_num_results > num_results[i]) ||
          check_results(
            results + offset, num_results[i], ref_results, ref_num_results))
        PUT_ERR("ERROR in yac_bnd_sphere_part_search_do_bnd_circle_search");

      offset += num_results[i];
    }
    free(ref_results);
    free(results);
    free(num_results);
  }

  yac_bnd_sphere_part_search_delete(search);
  free(point_coords);
  free(bnd_circles);
}

static void empty_test() {

  size_t count = 32;
  struct bounding_circle * bnd_circles =
    xmalloc(count * sizeof(*bnd_circles));
  double (*point_coords)[3] = xmalloc(count * sizeof(*point_coords));

  for (size_t i = 0; i < count; ++i) {

    gen_random_bnd_circle(bnd_circles + i);
    gen_random_point(point_coords[i]);
  }

  struct bnd_sphere_part_search * search =
    yac_bnd_sphere_part_search_new(NULL, 0);

  { // check yac_bnd_sphere_part_search_do_point_search
    size_t * results;
    size_t * num_results = xmalloc(count * sizeof(*num_results));
    yac_bnd_sphere_part_search_do_point_search(
      search, point_coords, count, &results, num_results);

    for (size_t i = 0; i < count; ++i)
      if (num_results[i] != 0)
        PUT_ERR("ERROR in yac_bnd_sphere_part_search_do_point_search");

    free(results);
    free(num_results);
  }

  { // yac_bnd_sphere_part_search_do_bnd_circle_search
    size_t * results;
    size_t * num_results = xmalloc(count * sizeof(*num_results));
    yac_bnd_sphere_part_search_do_bnd_circle_search(
      search, bnd_circles, count, &results, num_results);

    for (size_t i = 0; i < count; ++i)
      if (num_results[i] != 0)
        PUT_ERR("ERROR in yac_bnd_sphere_part_search_do_bnd_circle_search");

    free(results);
    free(num_results);
  }

  yac_bnd_sphere_part_search_delete(search);
  free(point_coords);
  free(bnd_circles);
}

static void zero_balance_point() {

  struct bounding_circle zero_balance_point_bnd_circles[6] =
    {{.base_vector = {1,0,0}, .inc_angle = SIN_COS_ZERO, .sq_crd = DBL_MAX},
     {.base_vector = {-1,0,0}, .inc_angle = SIN_COS_ZERO, .sq_crd = DBL_MAX},
     {.base_vector = {0,1,0}, .inc_angle = SIN_COS_ZERO, .sq_crd = DBL_MAX},
     {.base_vector = {0,-1,0}, .inc_angle = SIN_COS_ZERO, .sq_crd = DBL_MAX},
     {.base_vector = {0,0,1}, .inc_angle = SIN_COS_ZERO, .sq_crd = DBL_MAX},
     {.base_vector = {0,0,-1}, .inc_angle = SIN_COS_ZERO, .sq_crd = DBL_MAX}};

  struct bnd_sphere_part_search * search =
    yac_bnd_sphere_part_search_new(
      zero_balance_point_bnd_circles,
      sizeof(zero_balance_point_bnd_circles)/
        sizeof(zero_balance_point_bnd_circles[0]));

  yac_bnd_sphere_part_search_delete(search);
}

static void mean_point_test() {

  struct bounding_circle pole_balance_point_bnd_circles[5] =
    {{.base_vector = {1,0,0}, .inc_angle = SIN_COS_M_PI_2, .sq_crd = DBL_MAX},
     {.base_vector = {-1,0,0}, .inc_angle = SIN_COS_M_PI_2, .sq_crd = DBL_MAX},
     {.base_vector = {0,1,0}, .inc_angle = SIN_COS_M_PI_2, .sq_crd = DBL_MAX},
     {.base_vector = {0,-1,0}, .inc_angle = SIN_COS_M_PI_2, .sq_crd = DBL_MAX},
     {.base_vector = {0,0,1}, .inc_angle = SIN_COS_M_PI_2, .sq_crd = DBL_MAX}};

  double search_points[6][3] =
    {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

  struct bnd_sphere_part_search * search =
    yac_bnd_sphere_part_search_new(
      pole_balance_point_bnd_circles,
      sizeof(pole_balance_point_bnd_circles)/
        sizeof(pole_balance_point_bnd_circles[0]));

  size_t num_points = sizeof(search_points)/sizeof(search_points[0]);
  size_t * results;
  size_t num_results[num_points];
  yac_bnd_sphere_part_search_do_point_search(
    search, search_points, num_points, &results, num_results);

  size_t ref_results[6][5] =
    {{0,2,3,4},{1,2,3,4},{0,1,3,4},{0,1,2,4},{0,1,2,3,4},{0,1,2,3}};
  size_t ref_num_results[6] = {4,4,4,4,5,4};

  for (size_t i = 0, offset = 0; i < num_points; offset += num_results[i++])
    if ((ref_num_results[i] > num_results[i]) ||
        check_results(
          results + offset, num_results[i], ref_results[i], ref_num_results[i]))
      PUT_ERR("ERROR in yac_bnd_sphere_part_search_do_point_search");

  free(results);
  yac_bnd_sphere_part_search_delete(search);
}
