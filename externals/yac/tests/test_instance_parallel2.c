/**
 * @file test_instance_parallel2.c
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
#include <unistd.h>
#include <string.h>

#include <mpi.h>

#include "tests.h"
#include "test_common.h"
#include "instance.h"
#include "utils.h"
#include "interp_method_check.h"
#include "yac_interface.h"
#include "yac_mpi.h"
#include "event.h"

static void check_constructor(void * user_data);

static void check_do_search(
  yac_int const * global_ids, double const (*coordinates_xyz)[3],
  size_t count, void * user_data);

struct field_config {
  size_t grid_idx;
  size_t point_idx;
};

struct field_config_src_tgt {
  struct field_config src, tgt;
};

static int get_num_unique_configurations(
  struct field_config field_config_a[4],
  struct field_config field_config_b[4]);

enum DIRECTION {
  IN = 0,
  OUT = 1,
};

static void run_test_configuration(struct field_config test_config[2][3][4]);

int main (void) {

  // This test checks whether weights are reused for different fields
  // if grid and field data is identical. In case any of these two differs
  // weights cannot be reused.

  yac_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
  xt_initialize(MPI_COMM_WORLD);
  int rank, size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);

  if (size != 9) {

    PUT_ERR("ERROR: wrong number of processes\n");
    return TEST_EXIT_CODE;
  }

  int local_comp_idx = rank % 3;

  // dimensions: [test_idx][field_idx]
  struct field_config test_cases[][4] =
    {{{.grid_idx = 0, .point_idx = 0},
      {.grid_idx = 1, .point_idx = 1},
      {.grid_idx = 2, .point_idx = 2},
      {.grid_idx = 3, .point_idx = 3}},
     {{.grid_idx = 1, .point_idx = 1},
      {.grid_idx = 1, .point_idx = 2},
      {.grid_idx = 1, .point_idx = 2},
      {.grid_idx = 3, .point_idx = 2}},
     {{.grid_idx = 3, .point_idx = 2},
      {.grid_idx = 3, .point_idx = 2},
      {.grid_idx = 3, .point_idx = 2},
      {.grid_idx = 3, .point_idx = 2}},
     {{.grid_idx = 0, .point_idx = 0},
      {.grid_idx = 0, .point_idx = 1},
      {.grid_idx = 0, .point_idx = 2},
      {.grid_idx = 0, .point_idx = 3}},
     {{.grid_idx = 0, .point_idx = 0},
      {.grid_idx = 1, .point_idx = 0},
      {.grid_idx = 2, .point_idx = 0},
      {.grid_idx = 3, .point_idx = 0}},
     {{.grid_idx = 0, .point_idx = 2},
      {.grid_idx = 1, .point_idx = 3},
      {.grid_idx = 0, .point_idx = 2},
      {.grid_idx = 1, .point_idx = 3}}};
  size_t num_test_cases = sizeof(test_cases) / sizeof(test_cases[0]);

  //[direction][other_comp]
  int check_constructor_call_count[2][3];
  int check_do_search_call_count[2][3];

  for (int comp_offset = -1; comp_offset <= 1; comp_offset += 2) {
    int other_comp_idx = (local_comp_idx  + comp_offset + 3)%3;
    for (int direction = 0; direction < 2; ++direction) {

      char check_constructor_key[256];
      sprintf(check_constructor_key, "check_constructor_%u_to_%u",
        (direction == IN)?(other_comp_idx + 1):(local_comp_idx  + 1),
        (direction == IN)?(local_comp_idx  + 1):(other_comp_idx + 1));
      yac_interp_method_check_add_constructor_callback(
        check_constructor,
        &(check_constructor_call_count[direction][other_comp_idx]),
        check_constructor_key);

      char check_do_search_key[256];
      sprintf(check_do_search_key, "check_do_search_%u_to_%u",
        (direction == IN)?(other_comp_idx + 1):(local_comp_idx  + 1),
        (direction == IN)?(local_comp_idx  + 1):(other_comp_idx + 1));
      yac_interp_method_check_add_do_search_callback(
        check_do_search,
        &(check_do_search_call_count[direction][other_comp_idx]),
        check_do_search_key);
    }
  }

  struct field_config test_config[2][3][4];

  // [direction][comp_idx]
  int test_idx[2][3];

  // for all in field configurations
  for (int test_idx_in = 0; test_idx_in < num_test_cases; ++test_idx_in) {

    for (int comp_idx = 0; comp_idx < 3; ++comp_idx)
      test_idx[IN][comp_idx] = test_idx_in;

    // check all out configurations agains the current in config
    for (int test_idx_out = 0, test_run_out = 0;
          test_run_out < (num_test_cases + 2) / 3; ++test_run_out) {

      for (int comp_idx = 0; comp_idx < 3; ++comp_idx,
            test_idx_out = (test_idx_out + 1)%num_test_cases)
        test_idx[OUT][comp_idx] = test_idx_out;

      for (int direction = 0; direction < 2; ++direction)
        for (int comp_idx = 0; comp_idx < 3; ++comp_idx)
          for (int field_idx = 0; field_idx < 4; ++field_idx)
            test_config[direction][comp_idx][field_idx] =
              test_cases[test_idx[direction][comp_idx]][field_idx];

      memset(check_constructor_call_count, 0,
            sizeof(check_constructor_call_count));
      memset(check_do_search_call_count, 0,
            sizeof(check_do_search_call_count));

      // run the current test configuration
      run_test_configuration(test_config);

      // check results
      for (int comp_offset = -1; comp_offset <= 1; comp_offset += 2) {
        int other_comp_idx = (local_comp_idx  + comp_offset + 3)%3;
        for (int direction = 0; direction < 2; ++direction)
          if ((check_constructor_call_count[direction][other_comp_idx] != 
                check_do_search_call_count[direction][other_comp_idx]) ||
              (check_constructor_call_count[direction][other_comp_idx] !=
               get_num_unique_configurations(
                 test_config[direction][local_comp_idx ],
                 test_config[direction^1][other_comp_idx])))
            PUT_ERR("invalid call count");
      }
    }
  }

  xt_finalize();
  yac_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);

  return TEST_EXIT_CODE;
}

static void check_constructor(void * user_data) {

  ++*(int*)user_data;

  return;
}

static void check_do_search(
  yac_int const * global_ids, double const (*coordinates_xyz)[3],
  size_t num_points, void * user_data) {

  UNUSED(global_ids);
  UNUSED(coordinates_xyz);
  UNUSED(num_points);

  ++*(int*)user_data;

  return;
}

static int compare_size_t(const void * a, const void * b) {

  return (*((const size_t*)a) > *((const size_t*)b)) -
         (*((const size_t*)a) < *((const size_t*)b));
}

static int compare_field_config(const void * a, const void * b) {

  int ret = compare_size_t(&(((const struct field_config*)a)->grid_idx),
                           &(((const struct field_config*)b)->grid_idx));
  if (ret) return ret;
  else return compare_size_t(&(((const struct field_config*)a)->point_idx),
                             &(((const struct field_config*)b)->point_idx));
}

static int compare_field_config_src_tgt(const void * a, const void * b) {

  int ret =
    compare_field_config(&(((const struct field_config_src_tgt*)a)->src),
                         &(((const struct field_config_src_tgt*)b)->src));
  if (ret) return ret;
  else
    return compare_field_config(&(((const struct field_config_src_tgt*)a)->tgt),
                                &(((const struct field_config_src_tgt*)b)->tgt));
}

static int get_num_unique_configurations(
  struct field_config field_config_a[4],
  struct field_config field_config_b[4]) {

  struct field_config_src_tgt field_configs[4];

  for (int i = 0; i < 4; ++i) {
    field_configs[i].src = field_config_a[i];
    field_configs[i].tgt = field_config_b[i];
  }

  qsort(
    field_configs, 4, sizeof(field_configs[0]), compare_field_config_src_tgt);

  int num_unique_config = 1;
  for (int i = 1; i < 4; ++i)
    if (compare_field_config_src_tgt(&field_configs[i-1], &field_configs[i]))
      ++num_unique_config;

  return num_unique_config;
}

static void run_test_configuration(struct field_config test_config[2][3][4]) {

  int rank;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);

  int comp_idx = rank % 3;

  struct yac_instance * instance =
    yac_instance_new(MPI_COMM_WORLD);

  // define couples
  yac_instance_def_datetime(
    instance, "2008-03-09T16:05:07", "2008-03-10T16:05:07");
  char * coupling_period =
    strdup(yac_time_to_ISO("60", YAC_TIME_UNIT_SECOND));
  for (size_t src_comp_idx = 0; src_comp_idx < 3; ++src_comp_idx) {
    char src_comp_id = '1' + (char)src_comp_idx;
    for (size_t tgt_comp_idx = 0; tgt_comp_idx < 3; ++tgt_comp_idx) {
      char tgt_comp_id = '1' + (char)tgt_comp_idx;
      if (src_comp_id == tgt_comp_id) continue;
      char src_comp_name[32], tgt_comp_name[32];
      sprintf(src_comp_name, "comp_%c", src_comp_id);
      sprintf(tgt_comp_name, "comp_%c", tgt_comp_id);
      char constructor_key[32], do_search_key[32];
      sprintf(constructor_key, "check_constructor_%c_to_%c",
              src_comp_id, tgt_comp_id);
      sprintf(do_search_key, "check_do_search_%c_to_%c",
              src_comp_id, tgt_comp_id);
      struct yac_interp_stack_config * interp_stack =
        yac_interp_stack_config_new();
      yac_interp_stack_config_add_check(
        interp_stack, constructor_key, do_search_key);
      yac_interp_stack_config_add_fixed(interp_stack, -1.0);
      for (size_t field_idx = 0; field_idx < 4; ++field_idx) {
        char field_id = '1' + (char)field_idx;
        char src_grid_name[32], tgt_grid_name[32];
        sprintf(src_grid_name, "grid_%c_%zu_out",
                src_comp_id, test_config[OUT][src_comp_idx][field_idx].grid_idx+1);
        sprintf(tgt_grid_name, "grid_%c_%zu_in",
                tgt_comp_id, test_config[IN][tgt_comp_idx][field_idx].grid_idx+1);
        char field_name[32];
        sprintf(field_name, "field_from_%c_to_%c_%c",
                src_comp_id, tgt_comp_id, field_id);
        yac_instance_def_couple(
          instance, src_comp_name, src_grid_name, field_name,
          tgt_comp_name, tgt_grid_name, field_name,
          coupling_period, YAC_REDUCTION_TIME_NONE, interp_stack, 60, 60,
          NULL, 0, 1.0, 0.0, 0, NULL, NULL);
      }
      yac_interp_stack_config_delete(interp_stack);
    }
  }
  free(coupling_period);

  char const * component_names[3] = {"comp_1", "comp_2", "comp_3"};
  yac_instance_def_components(
    instance, &(component_names[comp_idx]), 1);

  int comp_rank;
  MPI_Comm component_comm =
    yac_component_get_communicator(
      yac_instance_get_component(instance, 0));
  MPI_Comm_rank(component_comm, &comp_rank);

  size_t num_vertices[2] = {2,2};
  int cyclic[2] = {0,0};
  double coordinates_x[2];
  double coordinates_y[4][2] = {{0,1},{1,2},{2,3},{3,4}};

  coordinates_x[0] = comp_rank;
  coordinates_x[1] = comp_rank + 1;

  yac_int global_cell_ids[1] = {comp_rank};
  yac_int global_corner_ids[4] =
    {4 * comp_rank + 0, 4 * comp_rank + 1,
     4 * comp_rank + 2, 4 * comp_rank + 3};
  int cell_core_mask[1] = {1};
  int corner_core_mask[4] = {1,1,1,1};

  // [direction][grid_idx]
  struct yac_basic_grid * grid[2][4];
  // [direction][grid_idx][field_idx]
  struct interp_field interp_fields[2][4][4];

  char * timestep[2] = {strdup(yac_time_to_ISO("20", YAC_TIME_UNIT_SECOND)),
                        strdup(yac_time_to_ISO("10", YAC_TIME_UNIT_SECOND))};
  size_t collection_size = 1;
  char const * str_direction[2] = {"in", "out"};

  for (int grid_idx = 0; grid_idx < 4; ++grid_idx) {
    for (int direction = 0; direction < 2; ++direction) {
      char grid_name[256];
      sprintf(
        grid_name, "grid_%d_%d_%s", comp_idx+1,
        grid_idx+1, str_direction[direction]);
      struct basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg_2d_deg(
          num_vertices, cyclic, coordinates_x, coordinates_y[grid_idx]);
      grid_data.cell_ids = TO_POINTER(global_cell_ids);
      grid_data.vertex_ids = TO_POINTER(global_corner_ids);
      grid_data.core_cell_mask = TO_POINTER(cell_core_mask);
      grid_data.core_vertex_mask = TO_POINTER(corner_core_mask);
      grid[direction][grid_idx] = yac_basic_grid_new(grid_name, grid_data);
      for (int point_idx = 0; point_idx < 4; ++point_idx) {
        size_t coordinates_idx =
          yac_basic_grid_add_coordinates(
            grid[direction][grid_idx], CORNER,
            grid_data.vertex_coordinates, grid_data.num_vertices);
        interp_fields[direction]
                     [grid_idx]
                     [point_idx].location = CORNER;
        interp_fields[direction]
                     [grid_idx]
                     [point_idx].coordinates_idx = coordinates_idx;
        interp_fields[direction]
                     [grid_idx]
                     [point_idx].masks_idx = SIZE_MAX;
      }
    }
  }

  for (int comp_offset = -1; comp_offset <= 1; comp_offset += 2) {
    int other_comp_idx = (comp_idx + comp_offset + 3) % 3;
    for (int field_idx = 0; field_idx < 4; ++field_idx) {
      for (int direction = 0; direction < 2; ++direction) {
        char field_name[256];
        sprintf(field_name, "field_from_%u_to_%u_%u",
          ((direction == IN)?other_comp_idx:comp_idx)+1,
          ((direction == IN)?comp_idx:other_comp_idx)+1,
          field_idx + 1);
        struct field_config transient_config =
          test_config[direction][comp_idx][field_idx];
        int grid_idx = transient_config.grid_idx;
        int point_idx = transient_config.point_idx;
        yac_instance_add_field(
          instance, field_name, 0, grid[direction][grid_idx],
          &interp_fields[direction][grid_idx][point_idx], 1,
          collection_size, timestep[direction]);
      }
    }
  }
  free(timestep[1]);
  free(timestep[0]);

  yac_instance_setup(instance);

  yac_instance_delete(instance);

  for (int grid_idx = 0; grid_idx < 4; ++grid_idx)
    for (int direction = 0; direction < 2; ++direction)
      yac_basic_grid_delete(grid[direction][grid_idx]);
}
