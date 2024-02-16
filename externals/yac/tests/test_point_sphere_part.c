/**
 * @file test_point_sphere_part.c
 *
 * @copyright Copyright  (C)  2016 DKRZ, MPI-M
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

static void random_test(size_t const num_points);
static void zero_point_test();
static void single_point_test();
static void empty_branch_test(size_t num_points);
static void zero_balance_point();
static void duplicated_point_test();

int main(void) {

  srand(1365441);

  for (size_t i = 0; i < 10; ++i) {
    random_test(1);
    random_test(10);
    random_test(100);
    random_test(650);
  }
  zero_point_test();
  single_point_test();
  empty_branch_test(4);
  empty_branch_test(128);
  zero_balance_point();
  duplicated_point_test();

  return TEST_EXIT_CODE;
}

struct distance_index {
  double distance;
  size_t index;
};

static int compare_distance_index(const void * a, const void * b) {

  int ret = (((struct distance_index*)a)->distance >
             ((struct distance_index*)b)->distance) -
            (((struct distance_index*)a)->distance <
             ((struct distance_index*)b)->distance);
  if (ret == 0)
    ret = (((struct distance_index*)a)->index >
           ((struct distance_index*)b)->index) -
          (((struct distance_index*)a)->index <
           ((struct distance_index*)b)->index);

  return ret;
}

static void point_sphere_test_NNN(
  struct point_sphere_part_search * search, size_t n,
  double (*xyz_coordinates)[3], size_t num_points, size_t ** local_point_ids,
  size_t * local_point_ids_array_size,
  size_t * num_local_point_ids, struct distance_index * distances) {

  yac_point_sphere_part_search_NNN(
    search, num_points, xyz_coordinates, n, NULL, NULL, NULL, NULL,
    local_point_ids, local_point_ids_array_size, num_local_point_ids);

  // check results
  for (size_t j = 0, offset = 0; j < num_points; ++j) {

    struct distance_index * curr_distances = distances + j * num_points;

    size_t ref_num_local_points = 1, m = 0;
    for (; (ref_num_local_points < num_points); ++ref_num_local_points)
      if (double_are_unequal(
            curr_distances[ref_num_local_points-1].distance,
            curr_distances[ref_num_local_points].distance))
        if (++m == n) break;

    if (num_local_point_ids[j] != ref_num_local_points) {
      PUT_ERR("wrong number of local points found\n");
      continue;
    }

    size_t num_matching_points = 0;
    for (size_t k = 0; k < ref_num_local_points; ++k)
      for (size_t l = 0; l < ref_num_local_points; ++l)
        if (curr_distances[k].index == (*local_point_ids)[offset+l])
          ++num_matching_points;

    if (num_matching_points != ref_num_local_points)
      PUT_ERR("wrong results\n");

    offset += ref_num_local_points;
  }
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

static void random_test(size_t const num_points) {

  double (*coordinates_xyz)[3];
  double (*coordinates_xyz_a)[3];
  double (*coordinates_xyz_b)[3];
  struct bounding_circle * bnd_circles =
    xmalloc(num_points * sizeof(*bnd_circles));
  yac_int * global_ids = xmalloc(num_points * sizeof(*global_ids));

  coordinates_xyz = xmalloc(2 * num_points * sizeof(*coordinates_xyz));
  coordinates_xyz_a = coordinates_xyz;
  coordinates_xyz_b = coordinates_xyz + num_points;

  size_t * local_point_ids = NULL;
  size_t local_point_ids_array_size = 0;
  size_t * num_local_point_ids =
    xmalloc(num_points * sizeof(*num_local_point_ids));
  struct distance_index * distances =
    xmalloc(num_points * num_points * sizeof(*distances));

  for (size_t i = 0; i < num_points; ++i) {
    gen_random_bnd_circle(bnd_circles + i);
    global_ids[i] = (yac_int)i;
  }
  for (size_t k = 0; k < 2 * num_points; ++k)
    gen_random_point(coordinates_xyz[k]);

  for (size_t j = 0; j < num_points; ++j) {
    for (size_t k = 0; k < num_points; ++k) {

      distances[j * num_points + k].distance =
      get_vector_angle(coordinates_xyz_a[k], coordinates_xyz_b[j]);
      distances[j * num_points + k].index = k;
    }
    qsort(distances + j * num_points, num_points, sizeof(*distances),
          compare_distance_index);
  }

  struct point_sphere_part_search * search =
    yac_point_sphere_part_search_new(
      num_points, coordinates_xyz_a, global_ids);

  for (size_t n = 1; n < 16; ++n)
    point_sphere_test_NNN(
      search, n, coordinates_xyz_b, num_points, &local_point_ids,
      &local_point_ids_array_size, num_local_point_ids,  distances);


  { // yac_point_sphere_part_search_bnd_circle_contains_points
    for (size_t i = 0; i < num_points; ++i) {
      int ref_flag = 0;
      for (size_t j = 0; (j < num_points) && !ref_flag; ++j)
        ref_flag =
          yac_point_in_bounding_circle_vec(
            coordinates_xyz_a[j], bnd_circles + i);
      if (ref_flag !=
          yac_point_sphere_part_search_bnd_circle_contains_points(
            search, bnd_circles[i]))
        PUT_ERR(
          "ERROR in yac_point_sphere_part_search_bnd_circle_contains_points");
    }
  }

  yac_delete_point_sphere_part_search(search);

  free(distances);
  free(local_point_ids);
  free(num_local_point_ids);
  free(global_ids);
  free(bnd_circles);
  free(coordinates_xyz);
}

static void single_point_test() {

  struct point_sphere_part_search * search;
  {
    double x_coordinate = 0, y_coordinate = 0;
    double xyz_coordinate[1][3];
    yac_int global_id[1] = {0};
    LLtoXYZ(x_coordinate, y_coordinate, xyz_coordinate[0]);
    search = yac_point_sphere_part_search_new(1, xyz_coordinate, global_id);
  }

  double x_coordinates[] = {0.0,90.0,180.0,-90.0, 0.0,  0.0};
  double y_coordinates[] = {0.0, 0.0,  0.0,  0.0,90.0,-90.0};
  size_t num_points = sizeof(x_coordinates) / sizeof(x_coordinates[0]);
  double xyz_coordinates[num_points][3];
  for (size_t i = 0; i < num_points; ++i)
    LLtoXYZ_deg(x_coordinates[i], y_coordinates[i], xyz_coordinates[i]);

  size_t * local_point_ids = NULL;
  size_t local_point_ids_array_size = 0;
  size_t * num_local_point_ids =
    xmalloc(num_points * sizeof(*num_local_point_ids));

  yac_point_sphere_part_search_NN(
    search, num_points, xyz_coordinates, NULL, NULL, NULL, &local_point_ids,
    &local_point_ids_array_size, num_local_point_ids);

  for (size_t i = 0; i < num_points; ++i)
    if ((num_local_point_ids[i] != 1) || (local_point_ids[i] != 0))
      PUT_ERR("single_point_test: wrong result\n");

  free(local_point_ids);
  free(num_local_point_ids);
  yac_delete_point_sphere_part_search(search);
}

static void zero_point_test() {

  struct point_sphere_part_search * search =
    yac_point_sphere_part_search_new(0, NULL, NULL);

  double x_coordinates[] = {0.0,90.0,180.0,-90.0, 0.0,  0.0};
  double y_coordinates[] = {0.0, 0.0,  0.0,  0.0,90.0,-90.0};
  size_t num_points = sizeof(x_coordinates) / sizeof(x_coordinates[0]);
  double xyz_coordinates[num_points][3];
  for (size_t i = 0; i < num_points; ++i)
    LLtoXYZ_deg(x_coordinates[i], y_coordinates[i], xyz_coordinates[i]);

  size_t * local_point_ids = NULL;
  size_t local_point_ids_array_size = 0;
  size_t * num_local_point_ids =
    xmalloc(num_points * sizeof(*num_local_point_ids));

  for (size_t i = 0; i < num_points; ++i) num_local_point_ids[i] = -1;
  yac_point_sphere_part_search_NN(
    search, num_points, xyz_coordinates, NULL, NULL, NULL, &local_point_ids,
    &local_point_ids_array_size, num_local_point_ids);

  for (size_t i = 0; i < num_points; ++i)
    if (num_local_point_ids[i] != 0) PUT_ERR("zero_point_test: wrong result\n");

  for (size_t i = 0; i < num_points; ++i) num_local_point_ids[i] = -1;
  yac_point_sphere_part_search_NNN(
    search, num_points, xyz_coordinates, 16, NULL, NULL, NULL, NULL,
    &local_point_ids, &local_point_ids_array_size, num_local_point_ids);

  for (size_t i = 0; i < num_points; ++i)
    if (num_local_point_ids[i] != 0) PUT_ERR("zero_point_test: wrong result\n");

  free(local_point_ids);
  free(num_local_point_ids);
  yac_delete_point_sphere_part_search(search);
}

static void empty_branch_test(size_t num_points_a) {

  coordinate_pointer xyz_coordinates_a =
    malloc(num_points_a * sizeof(*xyz_coordinates_a));
  yac_int * global_ids = malloc(num_points_a * sizeof(*global_ids));
  for (size_t i = 0; i < num_points_a; ++i) {
    LLtoXYZ_deg(0.0, -5.0 + (10.0 / (double)num_points_a) * (double)i,
                xyz_coordinates_a[i]);
    global_ids[i] = (yac_int)i;
  }

  struct point_sphere_part_search * search =
    yac_point_sphere_part_search_new(
      num_points_a, xyz_coordinates_a, global_ids);
  free(global_ids);

  double x_coordinates_b[] = {-1.0,-1.0, 1.0,1.0};
  double y_coordinates_b[] = {-1.0, 1.0,-1.0,1.0};
  size_t num_points_b = sizeof(x_coordinates_b) / sizeof(x_coordinates_b[0]);
  double xyz_coordinates_b[num_points_b][3];
  for (size_t i = 0; i < num_points_b; ++i)
    LLtoXYZ_deg(x_coordinates_b[i], y_coordinates_b[i], xyz_coordinates_b[i]);

  size_t * local_point_ids = NULL;
  size_t local_point_ids_array_size = 0;
  size_t * num_local_point_ids =
    xmalloc(num_points_b * sizeof(*num_local_point_ids));

  struct distance_index distances[num_points_a * num_points_b];

  for (size_t i = 0; i < num_points_b; ++i) {

    for (size_t j = 0; j < num_points_a; ++j) {

      distances[i * num_points_a + j].distance =
        get_vector_angle(xyz_coordinates_a[j], xyz_coordinates_b[i]);
      distances[i * num_points_a + j].index = j;
    }
    qsort(distances + i * num_points_a, num_points_a, sizeof(*distances),
          compare_distance_index);
  }

  for (size_t i = 0; i < num_points_b; ++i) num_local_point_ids[i] = -1;
  yac_point_sphere_part_search_NN(
    search, num_points_b, xyz_coordinates_b, NULL, NULL, NULL,
    &local_point_ids, &local_point_ids_array_size, num_local_point_ids);

    // check results
  for (size_t j = 0, offset = 0; j < num_points_b; ++j) {

    struct distance_index * curr_distances = distances + j * num_points_a;

    size_t ref_num_local_points;
    for (ref_num_local_points = 1;
         (ref_num_local_points < num_points_a); ++ref_num_local_points)
      if (double_are_unequal(
            curr_distances[ref_num_local_points-1].distance,
            curr_distances[ref_num_local_points].distance)) break;

    if (num_local_point_ids[j] != ref_num_local_points) {
      PUT_ERR("wrong number of local points found\n");
      continue;
    }

    size_t num_matching_points = 0;
    for (size_t k = 0; k < ref_num_local_points; ++k)
      for (size_t l = 0; l < ref_num_local_points; ++l)
        if (curr_distances[k].index == local_point_ids[offset+l])
          ++num_matching_points;

    if (num_matching_points != ref_num_local_points)
      PUT_ERR("wrong results\n");

    offset += ref_num_local_points;
  }

  free(xyz_coordinates_a);
  free(local_point_ids);
  free(num_local_point_ids);
  yac_delete_point_sphere_part_search(search);
}

static void zero_balance_point() {

  struct point_sphere_part_search * search;
  {
    double x_coordinate = 0, y_coordinate = 0;
    double xyz_coordinate[6][3] =
      {{1,0,0}, {-1,0,0},
       {0,1,0}, {0,-1,0},
       {0,0,1}, {0,0,-1}};
    yac_int global_ids[6] = {0,1,2,3,4,5};
    LLtoXYZ(x_coordinate, y_coordinate, xyz_coordinate[0]);
    search = yac_point_sphere_part_search_new(6, xyz_coordinate, global_ids);
  }
  yac_delete_point_sphere_part_search(search);
}

static void duplicated_point_test() {

  struct point_sphere_part_search * search;
  {
    double x_coordinates[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                              1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                              2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
                              3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,
                              4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,
                              2.0,2.0,2.0};
    double y_coordinates[] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,
                              0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,
                              0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,
                              0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,
                              0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,
                              4.0,4.0,4.0};
    size_t num_points = sizeof(x_coordinates) / sizeof(x_coordinates[0]);
    double xyz_coordinates[num_points][3];
    yac_int global_ids[] = {
      0,1,2,3,4,5,6,7,8,9,
      10,11,12,13,14,15,16,17,18,19,
      20,21,22,23,24,25,26,27,28,29,
      30,31,32,33,34,35,36,37,38,39,
      40,41,42,43,44,45,46,47,48,49,
      100,100,101};
    for (size_t i = 0; i < num_points; ++i) {
      LLtoXYZ_deg(x_coordinates[i], y_coordinates[i], xyz_coordinates[i]);
      global_ids[i] = (yac_int)i;
    }
    search =
      yac_point_sphere_part_search_new(
        num_points, xyz_coordinates, global_ids);
  }

  double x_coordinates[] = {1.0,1.0,2.0,2.0,2.0,3.0,3.0,2.1,0.1};
  double y_coordinates[] = {1.0,8.0,3.0,4.0,5.0,1.0,8.0,4.0,0.0};
  size_t num_points = sizeof(x_coordinates) / sizeof(x_coordinates[0]);
  double xyz_coordinates[num_points][3];
  for (size_t i = 0; i < num_points; ++i)
    LLtoXYZ_deg(x_coordinates[i], y_coordinates[i], xyz_coordinates[i]);

  size_t * local_point_ids = NULL;
  size_t local_point_ids_array_size = 0;
  size_t num_local_point_ids[num_points];

  double (*result_coordinates_xyz)[3] = NULL;
  size_t result_coordinates_xyz_array_size = 0;

  yac_point_sphere_part_search_NN(
    search, num_points, xyz_coordinates, NULL,
    &result_coordinates_xyz, &result_coordinates_xyz_array_size,
    &local_point_ids, &local_point_ids_array_size, num_local_point_ids);

  size_t ref_local_point_ids[9][3] =
    {{11}, {18}, {23}, {52}, {25}, {31}, {38}, {52}, {0}};
  size_t ref_num_local_point_ids[9] = {1,1,1,1,1,1,1,1,1};
  double ref_x_coordinates[9] = {1.0,1.0,2.0,2.0,2.0,3.0,3.0,2.0,0.0};
  double ref_y_coordinates[9] = {1.0,8.0,3.0,4.0,5.0,1.0,8.0,4.0,0.0};
  double ref_result_coordinates_xyz[9][3][3];

  for (int i = 0; i < num_points; ++i)
    for (size_t j = 0; j < ref_num_local_point_ids[i]; ++j)
      LLtoXYZ_deg(
        ref_x_coordinates[i], ref_y_coordinates[i],
        ref_result_coordinates_xyz[i][j]);

  for (size_t i = 0, offset = 0; i < num_points;
       offset += num_local_point_ids[i++]) {

    if (num_local_point_ids[i] != ref_num_local_point_ids[i]) {
      PUT_ERR("ERROR in yac_point_sphere_part_search_NN (num_local_point_ids)");
    } else {

      for (size_t j = 0; j < ref_num_local_point_ids[i]; ++j) {

        if (!points_are_identically(
               result_coordinates_xyz[offset+j],
               ref_result_coordinates_xyz[i][j]))
          PUT_ERR("ERROR in yac_point_sphere_part_search_NN "
                  "(result_coordinates_xyz)");
        if (local_point_ids[offset+j] != ref_local_point_ids[i][j])
          PUT_ERR("ERROR in yac_point_sphere_part_search_NN (local_point_ids)");
      }
    }
  }

  free(local_point_ids);
  free(result_coordinates_xyz);
  yac_delete_point_sphere_part_search(search);
}
