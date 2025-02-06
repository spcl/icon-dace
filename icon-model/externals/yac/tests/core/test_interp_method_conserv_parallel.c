// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tests.h"
#include "interp_method.h"
#include "interp_method_conserv.h"
#include "interp_method_fixed.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"
#include "geometry.h"

#include <mpi.h>
#include <yaxt.h>

char const src_grid_name[] = "src_grid";
char const tgt_grid_name[] = "tgt_grid";

// conservative remapping with one source and one target process
static void test1();

// conservative remapping with two source and one target process
static void test2();

// conservative remapping with three source and two target processes
static void test3();

// conservative remapping with three source and two target processes
// this test checks all possible configurations for the conservative
// remapping
static void test4();

// conservative remapping with three source and two target processes
// this test checks how the interpolation method handles the case in which
// it does not receive any target points because they are masked
static void test5();

// 2nd order conservative remapping with one source and one target process
static void test6();

// 2nd order conservative remapping with three sources and one target process
static void test7();

// 2nd order conservative remapping with three sources and one target process
static void test8();

// conservative remapping with 4 source and 4 target processes
// contain partially overlapping cells and non-covered tgt cells
static void test9();

// target cell do not overlap with any non-masked source cell
static void test10();
static void test11();

// 2nd order conservative remapping with three sources and one target process
static void test12();

// only half of the target cells overlap with source grid
// (replicates a bug found in the packing of the final results)
static void test13();

// check a bug found in 2nd order conservative
static void test14();

int main (void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
  test8();
  test9();
  test10();
  test11();
  test12();
  test13();
  test14();

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void test1() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 2) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 2;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank == 1;
  int is_target = comm_rank == 0;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 2x2 grid:
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02

    double coordinates_x[] = {-1,0,1};
    double coordinates_y[] = {-1,0,1};
    size_t const num_cells[2] = {2,2};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{2,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;


    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    double coordinates_x[] = {-0.5,0.5};
    double coordinates_y[] = {-0.5,0.5};
    size_t const num_cells[2] = {1,1};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{1,1}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  int order = 1;
  int enforced_conserv = 0;
  int partial_coverage = 0;
  enum yac_interp_method_conserv_normalisation normalisation =
    YAC_INTERP_CONSERV_DESTAREA;

  struct interp_method * method_stack[] = {
    yac_interp_method_conserv_new(
      order, enforced_conserv, partial_coverage, normalisation), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  enum yac_interp_weights_reorder_type reorder_type[2] =
    {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

  for (size_t j = 0; j < 2; ++j) {
    for (size_t collection_size = 1; collection_size < 4;
         collection_size += 2) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      if (is_source) {
        double *** src_data = xmalloc(collection_size * sizeof(*src_data));
        for (size_t i = 0; i < collection_size; ++i)
          src_data[i] = xmalloc(1 * sizeof(**src_data));

        struct yac_basic_grid_data * src_grid_data =
          yac_basic_grid_get_data(src_grid);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          double * src_field =
            xmalloc(src_grid_data->num_cells * sizeof(*src_field));
          for (size_t i = 0; i < src_grid_data->num_cells; ++i)
            src_field[i] =
              (src_grid_data->core_cell_mask[i])?
                (double)(src_grid_data->cell_ids[i] + 1) +
                (double)(collection_idx * 4):
                -1.0;
          src_data[collection_idx][0] = src_field;
        }

        yac_interpolation_execute_put(interpolation, src_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          free(src_data[collection_idx][0]);
          free(src_data[collection_idx]);
        }
        free(src_data);
      }
      if (is_target) {
        double ref_tgt_field[1] = {2.5};

        struct yac_basic_grid_data * tgt_grid_data =
          yac_basic_grid_get_data(tgt_grid);

        double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          tgt_data[collection_idx] =
            xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
          for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
            tgt_data[collection_idx][k] = -1;
        }

        yac_interpolation_execute_get(interpolation, tgt_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          for (size_t j = 0, offset = collection_idx * 4;
               j < tgt_grid_data->num_cells; ++j) {
            if (tgt_grid_data->core_cell_mask[j] &&
                (ref_tgt_field[tgt_grid_data->cell_ids[j]] != -1.0)) {
              if (fabs((ref_tgt_field[tgt_grid_data->cell_ids[j]] +
                        (double)offset) - tgt_data[collection_idx][j]) > 1e-3)
                PUT_ERR("wrong interpolation result");
            } else {
              if (tgt_data[collection_idx][j] != -1.0)
                PUT_ERR("wrong interpolation result");
            }
          }
        }

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          free(tgt_data[collection_idx]);
        free(tgt_data);
      }

      yac_interpolation_delete(interpolation);
    }
  }

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test2() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 3) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 3;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 1;
  int is_target = comm_rank < 1;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 2x2 grid:
    //
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02

    double coordinates_x[] = {-1,0,1};
    double coordinates_y[] = {-1,0,1};
    size_t const num_cells[2] = {2,2};
    size_t local_start[2][2] = {{0,0},{1,0}};
    size_t local_count[2][2] = {{1,2},{1,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    double coordinates_x[] = {-0.5,0.5};
    double coordinates_y[] = {-0.5,0.5};
    size_t const num_cells[2] = {1,1};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{1,1}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  int order = 1;
  int enforced_conserv = 1;
  int partial_coverage = 0;
  enum yac_interp_method_conserv_normalisation normalisation =
    YAC_INTERP_CONSERV_DESTAREA;

  struct interp_method * method_stack[] = {
    yac_interp_method_conserv_new(
      order, enforced_conserv, partial_coverage, normalisation), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  enum yac_interp_weights_reorder_type reorder_type[2] =
    {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

  for (size_t j = 0; j < 2; ++j) {
    for (size_t collection_size = 1; collection_size < 4;
         collection_size += 2) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      if (is_source) {
        double *** src_data = xmalloc(collection_size * sizeof(*src_data));
        for (size_t i = 0; i < collection_size; ++i)
          src_data[i] = xmalloc(1 * sizeof(**src_data));

        struct yac_basic_grid_data * src_grid_data =
          yac_basic_grid_get_data(src_grid);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          src_data[collection_idx][0] =
            xmalloc(src_grid_data->num_cells * sizeof(***src_data));
          for (size_t i = 0; i < src_grid_data->num_cells; ++i)
            src_data[collection_idx][0][i] =
              (src_grid_data->core_cell_mask[i])?
                (double)(src_grid_data->cell_ids[i] + 1) +
                (double)(collection_idx * 4):
                -1.0;
        }

        yac_interpolation_execute_put(interpolation, src_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          free(src_data[collection_idx][0]);
          free(src_data[collection_idx]);
        }
        free(src_data);
      }
      if (is_target) {
        double ref_tgt_field[1] = {2.5};

        struct yac_basic_grid_data * tgt_grid_data =
          yac_basic_grid_get_data(tgt_grid);

        double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          tgt_data[collection_idx] =
            xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
          for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
            tgt_data[collection_idx][k] = -1;
        }

        yac_interpolation_execute_get(interpolation, tgt_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          for (size_t j = 0, offset = collection_idx * 4;
               j < tgt_grid_data->num_cells; ++j) {
            if (tgt_grid_data->core_cell_mask[j] &&
                (ref_tgt_field[tgt_grid_data->cell_ids[j]] != -1.0)) {
              if (fabs((ref_tgt_field[tgt_grid_data->cell_ids[j]] +
                        (double)offset) - tgt_data[collection_idx][j]) > 1e-3)
                PUT_ERR("wrong interpolation result");
            } else {
              if (tgt_data[collection_idx][j] != -1.0)
                PUT_ERR("wrong interpolation result");
            }
          }
        }

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          free(tgt_data[collection_idx]);
        free(tgt_data);
      }

      yac_interpolation_delete(interpolation);
    }
  }

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test3() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 5) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 5;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 2;
  int is_target = comm_rank < 2;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 3x3 grid:
    //
    //       12--21--13--22--14--23--15
    //       |       |       |       |
    //       15  06  17  07  19  08  20
    //       |       |       |       |
    //       08--14--09--16--10--18--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03

    double coordinates_x[] = {-1.5,-0.5,0.5,1.5};
    double coordinates_y[] = {-1.5,-0.5,0.5,1.5};
    size_t const num_cells[2] = {3,3};
    size_t local_start[3][2] = {{0,0},{1,0},{2,0}};
    size_t local_count[3][2] = {{1,3},{1,3},{1,3}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    // the global target grid is a 2x2 grid:
    //
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02

    double coordinates_x[] = {-1,0,1};
    double coordinates_y[] = {-1,0,1};
    size_t const num_cells[2] = {2,2};
    size_t local_start[2][2] = {{0,0},{0,1}};
    size_t local_count[2][2] = {{2,1},{2,1}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  int order = 1;
  int enforced_conserv = 0;
  int partial_coverage = 0;
  enum yac_interp_method_conserv_normalisation normalisation =
    YAC_INTERP_CONSERV_DESTAREA;

  struct interp_method * method_stack[] = {
    yac_interp_method_conserv_new(
      order, enforced_conserv, partial_coverage, normalisation), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  enum yac_interp_weights_reorder_type reorder_type[2] =
    {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

  for (size_t j = 0; j < 2; ++j) {
    for (size_t collection_size = 1; collection_size < 4;
         collection_size += 2) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      if (is_source) {
        double *** src_data = xmalloc(collection_size * sizeof(*src_data));
        for (size_t i = 0; i < collection_size; ++i)
          src_data[i] = xmalloc(1 * sizeof(**src_data));

        struct yac_basic_grid_data * src_grid_data =
          yac_basic_grid_get_data(src_grid);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          src_data[collection_idx][0] =
            xmalloc(src_grid_data->num_cells * sizeof(***src_data));
          for (size_t i = 0; i < src_grid_data->num_cells; ++i)
            src_data[collection_idx][0][i] =
              (src_grid_data->core_cell_mask[i])?
                (double)(src_grid_data->cell_ids[i] + 1) +
                (double)(collection_idx * 9):
                -1.0;
        }

        yac_interpolation_execute_put(interpolation, src_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          free(src_data[collection_idx][0]);
          free(src_data[collection_idx]);
        }
        free(src_data);
      }
      if (is_target) {
        double ref_tgt_field[4] = {3,4,6,7};

        struct yac_basic_grid_data * tgt_grid_data =
          yac_basic_grid_get_data(tgt_grid);

        double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          tgt_data[collection_idx] =
            xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
          for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
            tgt_data[collection_idx][k] = -1;
        }

        yac_interpolation_execute_get(interpolation, tgt_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          for (size_t j = 0, offset = collection_idx * 9;
               j < tgt_grid_data->num_cells; ++j) {
            if (tgt_grid_data->core_cell_mask[j] &&
                (ref_tgt_field[tgt_grid_data->cell_ids[j]] != -1.0)) {
              if (fabs((ref_tgt_field[tgt_grid_data->cell_ids[j]] +
                        (double)offset) - tgt_data[collection_idx][j]) > 1e-3)
                PUT_ERR("wrong interpolation result");
            } else {
              if (tgt_data[collection_idx][j] != -1.0)
                PUT_ERR("wrong interpolation result");
            }
          }
        }

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          free(tgt_data[collection_idx]);
        free(tgt_data);
      }

      yac_interpolation_delete(interpolation);
    }
  }

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test4() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 8) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 5;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 2;
  int is_target = comm_rank < 2;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 4x4 grid:
    //
    //       20--36--21--37--22--38--23--39--24
    //       |       |       |       |       |
    //       28  12  30  13  32  14  34  15  35
    //       |       |       |       |       |
    //       15--27--16--29--17--31--18--33--19
    //       |       |       |       |       |
    //       19  08  21  09  23  10  25  11  26
    //       |       |       |       |       |
    //       10--18--11--20--12--22--13--24--14
    //       |       |       |       |       |
    //       10  04  12  05  14  06  16  07  17
    //       |       |       |       |       |
    //       05--09--06--11--07--13--08--15--09
    //       |       |       |       |       |
    //       01  00  03  01  05  02  07  03  08
    //       |       |       |       |       |
    //       00--00--01--02--02--04--03--06--04
    //
    // mask for the source grid (cells with "x" are masked)
    //
    //       +---+---+---+---+
    //       | x | x | x |   |
    //       +---+---+---+---+
    //       | x | x |   |   |
    //       +---+---+---+---+
    //       | x |   |   |   |
    //       +---+---+---+---+
    //       |   |   |   |   |
    //       +---+---+---+---+

    double coordinates_x[] = {0,1,2,3,4};
    double coordinates_y[] = {0,1,2,3,4};
    size_t const num_cells[2] = {4,4};
    size_t local_start[3][2] = {{0,0},{0,1},{0,3}};
    size_t local_count[3][2] = {{4,1},{4,2},{4,1}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {
    // the global target grid is a 3x3 grid:
    //
    //       12--21--13--22--14--23--15
    //       |       |       |       |
    //       15  06  17  07  19  08  20
    //       |       |       |       |
    //       08--14--09--16--10--18--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03

    double coordinates_x[] = {0.5,1.5,2.5,3.5};
    double coordinates_y[] = {0.5,1.5,2.5,3.5};
    size_t const num_cells[2] = {3,3};
    size_t local_start[2][2] = {{0,0},{0,1}};
    size_t local_count[2][2] = {{3,1},{3,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  int * src_cell_mask = NULL;
  if (is_source) {
    int global_src_mask [] = {1,1,1,1,
                              0,1,1,1,
                              0,0,1,1,
                              0,0,0,1};

    struct yac_basic_grid_data * src_grid_data =
      yac_basic_grid_get_data(src_grid);

    src_cell_mask = xmalloc(src_grid_data->num_cells * sizeof(*src_cell_mask));
    for (size_t i = 0; i < src_grid_data->num_cells; ++i)
      src_cell_mask[i] = global_src_mask[src_grid_data->cell_ids[i]];
    yac_basic_grid_add_mask_nocpy(src_grid, YAC_LOC_CELL, src_cell_mask, NULL);
  }

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  enum yac_interp_method_conserv_normalisation
     normalisation[2] = {YAC_INTERP_CONSERV_DESTAREA,
                         YAC_INTERP_CONSERV_FRACAREA};

  for (int partial_coverage = 0; partial_coverage <= 1; ++partial_coverage) {
    for (int norm_idx = 0; norm_idx <= 1; ++norm_idx) {

      int order = 1;
      int enforced_conserv = 0;

      struct interp_method * method_stack[] = {
        yac_interp_method_conserv_new(
          order, enforced_conserv, partial_coverage,
          normalisation[norm_idx]), NULL};

      struct yac_interp_weights * weights =
        yac_interp_method_do_search(method_stack, interp_grid);

      yac_interp_method_delete(method_stack);

      enum yac_interp_weights_reorder_type reorder_type[2] =
        {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

      for (size_t j = 0; j < 2; ++j) {

        size_t collection_size = 1;
        struct yac_interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        if (is_source) {
          double *** src_data = xmalloc(collection_size * sizeof(*src_data));
          for (size_t i = 0; i < collection_size; ++i)
            src_data[i] = xmalloc(1 * sizeof(**src_data));

          struct yac_basic_grid_data * src_grid_data =
            yac_basic_grid_get_data(src_grid);

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            src_data[collection_idx][0] =
              xmalloc(src_grid_data->num_cells * sizeof(***src_data));
            for (size_t i = 0; i < src_grid_data->num_cells; ++i)
              src_data[collection_idx][0][i] =
                (src_grid_data->core_cell_mask[i])?
                  (double)(src_grid_data->cell_ids[i]) +
                  (double)(collection_idx * 16):
                  -1.0;
          }

          yac_interpolation_execute_put(interpolation, src_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            free(src_data[collection_idx][0]);
            free(src_data[collection_idx]);
          }
          free(src_data);
        }
        if (is_target) {
          // [norm_idx][partial_coverage][tgt_idx]
          double ref_tgt_field[2][2][9] =
            {{{-1,3.5,4.5,-1,-1,8.5,-1,-1,-1},
              {1.5,3.5,4.5,1.25,5.25,8.5,-1,2.5,9.0}},
             {{-1,3.5,4.5,-1,-1,8.5,-1,-1,-1},
              {2.0,3.5,4.5,5.0,7.0,8.5,-1,10.0,12.0}}};

          struct yac_basic_grid_data * tgt_grid_data =
            yac_basic_grid_get_data(tgt_grid);

          double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            tgt_data[collection_idx] =
              xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
            for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
              tgt_data[collection_idx][k] = -1;
          }

          yac_interpolation_execute_get(interpolation, tgt_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            for (size_t j = 0, offset = collection_idx * 16;
                 j < tgt_grid_data->num_cells; ++j) {
              if (tgt_grid_data->core_cell_mask[j] &&
                  (ref_tgt_field[norm_idx]
                                [partial_coverage]
                                [tgt_grid_data->cell_ids[j]] != -1.0)) {
                if (fabs((ref_tgt_field[norm_idx]
                                       [partial_coverage]
                                       [tgt_grid_data->cell_ids[j]] +
                          (double)offset) -
                         tgt_data[collection_idx][j]) > 1e-3)
                  PUT_ERR("wrong interpolation result");
              } else {
                if (tgt_data[collection_idx][j] != -1.0)
                  PUT_ERR("wrong interpolation result");
              }
            }
          }

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx)
            free(tgt_data[collection_idx]);
          free(tgt_data);
        }

        yac_interpolation_delete(interpolation);
      }
      yac_interp_weights_delete(weights);
    }
  }

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test5() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 5) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 5;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 2;
  int is_target = comm_rank < 2;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 4x4 grid:
    //
    //       20--36--21--37--22--38--23--39--24
    //       |       |       |       |       |
    //       28  12  30  13  32  14  34  15  35
    //       |       |       |       |       |
    //       15--27--16--29--17--31--18--33--19
    //       |       |       |       |       |
    //       19  08  21  09  23  10  25  11  26
    //       |       |       |       |       |
    //       10--18--11--20--12--22--13--24--14
    //       |       |       |       |       |
    //       10  04  12  05  14  06  16  07  17
    //       |       |       |       |       |
    //       05--09--06--11--07--13--08--15--09
    //       |       |       |       |       |
    //       01  00  03  01  05  02  07  03  08
    //       |       |       |       |       |
    //       00--00--01--02--02--04--03--06--04
    //
    // mask for the source grid (cells with "x" are masked)
    //
    //       +---+---+---+---+
    //       | x | x | x |   |
    //       +---+---+---+---+
    //       | x | x |   |   |
    //       +---+---+---+---+
    //       | x |   |   |   |
    //       +---+---+---+---+
    //       |   |   |   |   |
    //       +---+---+---+---+

    double coordinates_x[] = {0,1,2,3,4};
    double coordinates_y[] = {0,1,2,3,4};
    size_t const num_cells[2] = {4,4};
    size_t local_start[3][2] = {{0,0},{0,1},{0,3}};
    size_t local_count[3][2] = {{4,1},{4,2},{4,1}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {
    // the global target grid is a 3x3 grid:
    //
    //       12--21--13--22--14--23--15
    //       |       |       |       |
    //       15  06  17  07  19  08  20
    //       |       |       |       |
    //       08--14--09--16--10--18--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03
    //
    // target mask:
    //
    //      +---+---+---+
    //      | x | x | x |
    //      +---+---+---+
    //      | x | x | x |
    //      +---+---+---+
    //      |   |   |   |
    //      +---+---+---+

    double coordinates_x[] = {0.5,1.5,2.5,3.5};
    double coordinates_y[] = {0.5,1.5,2.5,3.5};
    size_t const num_cells[2] = {3,3};
    size_t local_start[2][2] = {{0,0},{0,1}};
    size_t local_count[2][2] = {{3,1},{3,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  int * src_cell_mask = NULL, * tgt_cell_mask = NULL;
  if (is_source) {
    int global_src_mask [] = {1,1,1,1,
                              0,1,1,1,
                              0,0,1,1,
                              0,0,0,1};

    struct yac_basic_grid_data * src_grid_data =
      yac_basic_grid_get_data(src_grid);

    src_cell_mask = xmalloc(src_grid_data->num_cells * sizeof(*src_cell_mask));
    for (size_t i = 0; i < src_grid_data->num_cells; ++i)
      src_cell_mask[i] = global_src_mask[src_grid_data->cell_ids[i]];
    yac_basic_grid_add_mask_nocpy(src_grid, YAC_LOC_CELL, src_cell_mask, NULL);
  }
  if (is_target) {
    int global_tgt_mask [] = {1,1,1,
                              0,0,0,
                              0,0,0};
    double coordinates_x [] = {1.0,2.0,3.0};
    double coordinates_y [] = {1.0,2.0,3.0};

    struct yac_basic_grid_data * tgt_grid_data =
      yac_basic_grid_get_data(tgt_grid);

    tgt_cell_mask = xmalloc(tgt_grid_data->num_cells * sizeof(*tgt_cell_mask));
    yac_coordinate_pointer tgt_cell_coords =
      xmalloc(tgt_grid_data->num_cells * sizeof(*tgt_cell_coords));
    for (size_t i = 0; i < tgt_grid_data->num_cells; ++i) {
      tgt_cell_mask[i] = global_tgt_mask[tgt_grid_data->cell_ids[i]];
      LLtoXYZ_deg(
        coordinates_x[tgt_grid_data->cell_ids[i]%3],
        coordinates_y[tgt_grid_data->cell_ids[i]/3],
        tgt_cell_coords[i]);
    }
    yac_basic_grid_add_mask_nocpy(
      tgt_grid, YAC_LOC_CELL, tgt_cell_mask, NULL);
    yac_basic_grid_add_coordinates_nocpy(
      tgt_grid, YAC_LOC_CELL, tgt_cell_coords);
  }

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  enum yac_interp_method_conserv_normalisation
     normalisation[2] = {YAC_INTERP_CONSERV_DESTAREA,
                         YAC_INTERP_CONSERV_FRACAREA};

  for (int partial_coverage = 0; partial_coverage <= 1; ++partial_coverage) {
    for (int norm_idx = 0; norm_idx <= 1; ++norm_idx) {

      int order = 1;
      int enforced_conserv = 0;

      struct interp_method * method_stack[] = {
        yac_interp_method_conserv_new(
          order, enforced_conserv, partial_coverage,
          normalisation[norm_idx]), NULL};

      struct yac_interp_weights * weights =
        yac_interp_method_do_search(method_stack, interp_grid);

      yac_interp_method_delete(method_stack);

      enum yac_interp_weights_reorder_type reorder_type[2] =
        {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

      for (size_t j = 0; j < 2; ++j) {

        size_t collection_size = 1;
        struct yac_interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        if (is_source) {
          double *** src_data = xmalloc(collection_size * sizeof(*src_data));
          for (size_t i = 0; i < collection_size; ++i)
            src_data[i] = xmalloc(1 * sizeof(**src_data));

          struct yac_basic_grid_data * src_grid_data =
            yac_basic_grid_get_data(src_grid);

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            src_data[collection_idx][0] =
              xmalloc(src_grid_data->num_cells * sizeof(***src_data));
            for (size_t i = 0; i < src_grid_data->num_cells; ++i)
              src_data[collection_idx][0][i] =
                (src_grid_data->core_cell_mask[i])?
                  (double)(src_grid_data->cell_ids[i]) +
                  (double)(collection_idx * 16):
                  -1.0;
          }

          yac_interpolation_execute_put(interpolation, src_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            free(src_data[collection_idx][0]);
            free(src_data[collection_idx]);
          }
          free(src_data);
        }
        if (is_target) {
          // [norm_idx][partial_coverage][tgt_idx]
          double ref_tgt_field[2][2][9] =
            {{{-1,3.5,4.5,-1,-1,-1,-1,-1,-1},
              {1.5,3.5,4.5,-1,-1,-1,-1,-1,-1}},
             {{-1,3.5,4.5,-1,-1,-1,-1,-1,-1},
              {2.0,3.5,4.5,-1,-1,-1,-1,-1,-1}}};

          struct yac_basic_grid_data * tgt_grid_data =
            yac_basic_grid_get_data(tgt_grid);

          double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            tgt_data[collection_idx] =
              xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
            for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
              tgt_data[collection_idx][k] = -1;
          }

          yac_interpolation_execute_get(interpolation, tgt_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx) {
            for (size_t j = 0, offset = collection_idx * 16;
                 j < tgt_grid_data->num_cells; ++j) {
              if (tgt_grid_data->core_cell_mask[j] &&
                  (ref_tgt_field[norm_idx]
                                [partial_coverage]
                                [tgt_grid_data->cell_ids[j]] != -1.0)) {
                if (fabs((ref_tgt_field[norm_idx]
                                       [partial_coverage]
                                       [tgt_grid_data->cell_ids[j]] +
                          (double)offset) -
                         tgt_data[collection_idx][j]) > 1e-3)
                  PUT_ERR("wrong interpolation result");
              } else {
                if (tgt_data[collection_idx][j] != -1.0)
                  PUT_ERR("wrong interpolation result");
              }
            }
          }

          for (size_t collection_idx = 0; collection_idx < collection_size;
               ++collection_idx)
            free(tgt_data[collection_idx]);
          free(tgt_data);
        }

        yac_interpolation_delete(interpolation);
      }
      yac_interp_weights_delete(weights);
    }
  }

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test6() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 2) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 2;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 1;
  int is_target = comm_rank < 1;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 3x3 grid:
    //
    //       12--21--13--22--14--23--15
    //       |       |       |       |
    //       15  06  17  07  19  08  20
    //       |       |       |       |
    //       08--14--09--16--10--18--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03

    double coordinates_x[] = {-1.5,-0.5,0.5,1.5};
    double coordinates_y[] = {-1.5,-0.5,0.5,1.5};
    size_t const num_cells[2] = {3,3};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{3,3}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    // the global target grid is a 2x2 grid:
    //
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02

    double coordinates_x[] = {-0.5,0,0.5};
    double coordinates_y[] = {-0.5,0,0.5};
    size_t const num_cells[2] = {2,2};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{2,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  int order = 2;
  int enforced_conserv = 0;
  int partial_coverage = 0;
  enum yac_interp_method_conserv_normalisation normalisation =
    YAC_INTERP_CONSERV_DESTAREA;

  struct interp_method * method_stack[] = {
    yac_interp_method_conserv_new(
      order, enforced_conserv, partial_coverage, normalisation), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  enum yac_interp_weights_reorder_type reorder_type[2] =
    {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

#define NUM_TESTS (4)

  for (size_t t = 0; t < NUM_TESTS; ++t) {
    for (size_t j = 0; j < 2; ++j) {

      size_t collection_size = 1;

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      if (is_source) {
        double global_field_data[NUM_TESTS][9] = {
          {1,1,1, 1,1,1, 1,1,1},
          {0,0,0, 1,1,1, 2,2,2},
          {0,1,2, 0,1,2, 0,1,2},
          {0,1,2, 1,2,3, 2,3,4}};

        double *** src_data = xmalloc(collection_size * sizeof(*src_data));
        for (size_t i = 0; i < collection_size; ++i)
          src_data[i] = xmalloc(1 * sizeof(**src_data));

        struct yac_basic_grid_data * src_grid_data =
          yac_basic_grid_get_data(src_grid);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          double * src_field =
            xmalloc(src_grid_data->num_cells * sizeof(*src_field));
          for (size_t i = 0; i < src_grid_data->num_cells; ++i)
            src_field[i] =
              (src_grid_data->core_cell_mask[i])?
                global_field_data[t][src_grid_data->cell_ids[i]]:-1.0;
          src_data[collection_idx][0] = src_field;
        }

        yac_interpolation_execute_put(interpolation, src_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          free(src_data[collection_idx][0]);
          free(src_data[collection_idx]);
        }
        free(src_data);
      }
      if (is_target) {
        double ref_tgt_field[NUM_TESTS][4] = {
          {1,1,
           1,1},
          {0.75,0.75,
           1.25,1.25},
          {0.75,1.25,
           0.75,1.25},
          {1.5,2.0,
           2.0,2.5}};

        struct yac_basic_grid_data * tgt_grid_data =
          yac_basic_grid_get_data(tgt_grid);

        double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          tgt_data[collection_idx] =
            xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
          for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
            tgt_data[collection_idx][k] = -1;
        }

        yac_interpolation_execute_get(interpolation, tgt_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          for (size_t j = 0; j < tgt_grid_data->num_cells; ++j) {
            if (tgt_grid_data->core_cell_mask[j] &&
                (ref_tgt_field[t][tgt_grid_data->cell_ids[j]] != -1.0)) {
              if (fabs((ref_tgt_field[t][tgt_grid_data->cell_ids[j]]) -
                       tgt_data[collection_idx][j]) > 1e-3)
                PUT_ERR("wrong interpolation result");
            } else {
              if (tgt_data[collection_idx][j] != -1.0)
                PUT_ERR("wrong interpolation result");
            }
          }
        }

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          free(tgt_data[collection_idx]);
        free(tgt_data);
      }

      yac_interpolation_delete(interpolation);
    }
  }
#undef NUM_TESTS

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test7() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 4) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 4;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 1;
  int is_target = comm_rank < 1;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 3x3 grid:
    //
    //       12--21--13--22--14--23--15
    //       |       |       |       |
    //       15  06  17  07  19  08  20
    //       |       |       |       |
    //       08--14--09--16--10--18--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03

    double coordinates_x[] = {-1.5,-0.5,0.5,1.5};
    double coordinates_y[] = {-1.5,-0.5,0.5,1.5};
    size_t const num_cells[2] = {3,3};
    size_t local_start[3][2] = {{0,0},{0,1},{0,2}};
    size_t local_count[3][2] = {{3,1},{3,1},{3,1}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    // the global target grid is a 2x2 grid:
    //
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02

    double coordinates_x[] = {-0.5,0,0.5};
    double coordinates_y[] = {-0.5,0,0.5};
    size_t const num_cells[2] = {2,2};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{2,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  int order = 2;
  int enforced_conserv = 0;
  int partial_coverage = 0;
  enum yac_interp_method_conserv_normalisation normalisation =
    YAC_INTERP_CONSERV_DESTAREA;

  struct interp_method * method_stack[] = {
    yac_interp_method_conserv_new(
      order, enforced_conserv, partial_coverage, normalisation), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  enum yac_interp_weights_reorder_type reorder_type[2] =
    {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

#define NUM_TESTS (4)

  for (size_t t = 0; t < NUM_TESTS; ++t) {
    for (size_t j = 0; j < 2; ++j) {

      size_t collection_size = 1;

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      if (is_source) {
        double global_field_data[NUM_TESTS][9] = {
          {1,1,1, 1,1,1, 1,1,1},
          {0,0,0, 1,1,1, 2,2,2},
          {0,1,2, 0,1,2, 0,1,2},
          {0,1,2, 1,2,3, 2,3,4}};
        double *** src_data = xmalloc(collection_size * sizeof(*src_data));
        for (size_t i = 0; i < collection_size; ++i)
          src_data[i] = xmalloc(1 * sizeof(**src_data));

        struct yac_basic_grid_data * src_grid_data =
          yac_basic_grid_get_data(src_grid);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          src_data[collection_idx][0] =
            xmalloc(src_grid_data->num_cells * sizeof(***src_data));
          for (size_t i = 0; i < src_grid_data->num_cells; ++i)
            src_data[collection_idx][0][i] =
              (src_grid_data->core_cell_mask[i])?
                global_field_data[t][src_grid_data->cell_ids[i]]:-1.0;
        }

        yac_interpolation_execute_put(interpolation, src_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          free(src_data[collection_idx][0]);
          free(src_data[collection_idx]);
        }
        free(src_data);
      }
      if (is_target) {
        double ref_tgt_field[NUM_TESTS][4] = {
          {1,1,
           1,1},
          {0.75,0.75,
           1.25,1.25},
          {0.75,1.25,
           0.75,1.25},
          {1.5,2.0,
           2.0,2.5}};

        struct yac_basic_grid_data * tgt_grid_data =
          yac_basic_grid_get_data(tgt_grid);

        double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          tgt_data[collection_idx] =
            xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
          for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
            tgt_data[collection_idx][k] = -1;
        }

        yac_interpolation_execute_get(interpolation, tgt_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          for (size_t j = 0; j < tgt_grid_data->num_cells; ++j) {
            if (tgt_grid_data->core_cell_mask[j] &&
                (ref_tgt_field[t][tgt_grid_data->cell_ids[j]] != -1.0)) {
              if (fabs((ref_tgt_field[t][tgt_grid_data->cell_ids[j]]) -
                       tgt_data[collection_idx][j]) > 1e-3)
                PUT_ERR("wrong interpolation result");
            } else {
              if (tgt_data[collection_idx][j] != -1.0)
                PUT_ERR("wrong interpolation result");
            }
          }
        }

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          free(tgt_data[collection_idx]);
        free(tgt_data);
      }

      yac_interpolation_delete(interpolation);
    }
  }
#undef NUM_TESTS

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test8() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 4) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 4;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 1;
  int is_target = comm_rank < 1;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 5x5 grid:
    //
    //       30--55--31--56--32--57--33-58--34-59--35
    //       |       |       |       |      |      |
    //       45  20  47  21  49  22  51 23  53 24  54
    //       |       |       |       |      |      |
    //       24--44--25--46--26--48--27-50--28-52--29
    //       |       |       |       |      |      |
    //       34  15  36  16  38  17  40 18  42 19  43
    //       |       |       |       |      |      |
    //       18--33--19--35--20--37--21-39--22-41--23
    //       |       |       |       |      |      |
    //       23  10  25  11  27  12  29 13  31 14  32
    //       |       |       |       |      |      |
    //       12--22--13--24--14--26--15-28--16-30--17
    //       |       |       |       |      |      |
    //       12  05  14  06  16  07  18 08  20 09  21
    //       |       |       |       |      |      |
    //       06--11--07--13--08--15--09-17--10-19--11
    //       |       |       |       |      |      |
    //       01  00  03  01  05  02  07 03  09 04  10
    //       |       |       |       |      |      |
    //       00--00--01--02--02--04--03-06--04-08--05

    double coordinates_x[] = {-2.5,-1.5,-0.5,0.5,1.5,2.5};
    double coordinates_y[] = {-2.5,-1.5,-0.5,0.5,1.5,2.5};
    size_t const num_cells[2] = {5,5};
    size_t local_start[3][2] = {{0,0},{0,2},{3,2}};
    size_t local_count[3][2] = {{5,2},{3,3},{2,3}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    // the global target grid is a 5x5 grid and is slightly shifted to the
    // source grid

    double coordinates_x[] = {-2,-1,0,1,2,3};
    double coordinates_y[] = {-3,-2,-1,0,1,2};
    size_t const num_cells[2] = {5,5};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{5,5}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  int order = 2;
  int enforced_conserv = 0;
  int partial_coverage = 1;
  enum yac_interp_method_conserv_normalisation normalisation =
    YAC_INTERP_CONSERV_DESTAREA;

  struct interp_method * method_stack[] = {
    yac_interp_method_conserv_new(
      order, enforced_conserv, partial_coverage, normalisation), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  enum yac_interp_weights_reorder_type reorder_type[2] =
    {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

#define NUM_TESTS (4)

  for (size_t t = 0; t < NUM_TESTS; ++t) {
    for (size_t j = 0; j < 2; ++j) {

      size_t collection_size = 1;

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      double global_src_field_data[NUM_TESTS][25] =
        {{1,1,1,1,1,
          1,1,1,1,1,
          1,1,1,1,1,
          1,1,1,1,1,
          1,1,1,1,1},
          {0.00,0.25,0.50,0.75,1.00,
          0.00,0.25,0.50,0.75,1.00,
          0.00,0.25,0.50,0.75,1.00,
          0.00,0.25,0.50,0.75,1.00,
          0.00,0.25,0.50,0.75,1.00},
          {0.00,0.00,0.00,0.00,0.00,
          0.25,0.25,0.25,0.25,0.25,
          0.50,0.50,0.50,0.50,0.50,
          0.75,0.75,0.75,0.75,0.75,
          1.00,1.00,1.00,1.00,1.00},
          {0.00,0.25,0.50,0.75,1.00,
          0.25,0.50,0.75,1.00,1.25,
          0.50,0.75,1.00,1.25,1.50,
          0.75,1.00,1.25,1.50,1.75,
          1.00,1.25,1.50,1.75,2.00}};

      size_t src_neighs[25][4] =
        {{0,1,5,0},{1,2,6,0},{2,3,7,1},{3,4,8,2},{4,4,9,3},
         {0,6,10,5},{1,7,11,5},{2,8,12,6},{3,9,13,7},{4,9,14,8},
         {5,11,15,10},{6,12,16,10},{7,13,17,11},{8,14,18,12},{9,14,19,13},
         {10,16,20,15},{11,17,21,15},{12,18,22,16},{13,19,23,17},{14,19,24,18},
         {15,21,20,20},{16,22,21,20},{17,23,22,21},{18,24,23,22},{19,24,24,23}};
      double src_gradient_weights_signs[4][4] =
        {{-1.0,1.0,1.0,-1.0},{-1.0,-1.0,1.0,1.0},
         {1.0,1.0,-1.0,-1.0},{1.0,-1.0,-1.0,1.0}};
      double src_gradient_weights[25][4] =
        {{0.25,0,0.25,0},{0.25,0.125,0.25,0.125},{0.25,0.125,0.25,0.125},
           {0.25,0.125,0.25,0.125},{0.25,0.25,0.25,0.25},
         {0.125,0,0.125,0},{0.125,0.125,0.125,0.125},{0.125,0.125,0.125,0.125},
           {0.125,0.125,0.125,0.125},{0.125,0.25,0.125,0.25},
         {0.125,0,0.125,0},{0.125,0.125,0.125,0.125},{0.125,0.125,0.125,0.125},
           {0.125,0.125,0.125,0.125},{0.125,0.25,0.125,0.25},
         {0.125,0,0.125,0},{0.125,0.125,0.125,0.125},{0.125,0.125,0.125,0.125},
           {0.125,0.125,0.125,0.125},{0.125,0.25,0.125,0.25},
         {0,0,0,0},{0,0.125,0,0.125},{0,0.125,0,0.125},
           {0,0.125,0,0.125},{0,0.25,0,0.25}};
      size_t tgt_to_src[25][4] =
        {{SIZE_MAX,SIZE_MAX,0,1}, {SIZE_MAX,SIZE_MAX,1,2},
           {SIZE_MAX,SIZE_MAX,2,3},{SIZE_MAX,SIZE_MAX,3,4},
           {SIZE_MAX,SIZE_MAX,4,SIZE_MAX},
         {0,1,5,6},{1,2,6,7},{2,3,7,8},{3,4,8,9},{4,SIZE_MAX,9,SIZE_MAX},
         {5,6,10,11},{6,7,11,12},{7,8,12,13},{8,9,13,14},
           {9,SIZE_MAX,14,SIZE_MAX},
         {10,11,15,16},{11,12,16,17},{12,13,17,18},{13,14,18,19},
           {14,SIZE_MAX,19,SIZE_MAX},
         {15,16,20,21},{16,17,21,22},{17,18,22,23},{18,19,23,24},
           {19,SIZE_MAX,24,SIZE_MAX}};

      // check generated interpolation
      if (is_source) {
        double *** src_data = xmalloc(collection_size * sizeof(*src_data));
        for (size_t i = 0; i < collection_size; ++i)
          src_data[i] = xmalloc(1 * sizeof(**src_data));

        struct yac_basic_grid_data * src_grid_data =
          yac_basic_grid_get_data(src_grid);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          src_data[collection_idx][0] =
            xmalloc(src_grid_data->num_cells * sizeof(***src_data));
          for (size_t i = 0; i < src_grid_data->num_cells; ++i)
            src_data[collection_idx][0][i] =
              (src_grid_data->core_cell_mask[i])?
                global_src_field_data[t][src_grid_data->cell_ids[i]]:-1.0;
        }

        yac_interpolation_execute_put(interpolation, src_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          free(src_data[collection_idx][0]);
          free(src_data[collection_idx]);
        }
        free(src_data);
      }
      if (is_target) {

        double ref_tgt_field[25];

        for (size_t i = 0; i < 25; ++i) {
          double tgt_field_value = 0.0;
          for (size_t j = 0; j < 4; ++j) {
            size_t src_global_id = tgt_to_src[i][j];
            if (src_global_id != SIZE_MAX) {
              tgt_field_value +=
                global_src_field_data[t][src_global_id];
              for (size_t k = 0; k < 4; ++k)
                tgt_field_value +=
                  global_src_field_data[t][src_neighs[src_global_id][k]] *
                  src_gradient_weights_signs[j][k] *
                  src_gradient_weights[src_global_id][k];
            }
          }
          ref_tgt_field[i] = tgt_field_value / 4.0;
        }

        struct yac_basic_grid_data * tgt_grid_data =
          yac_basic_grid_get_data(tgt_grid);

        double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          tgt_data[collection_idx] =
            xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
          for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
            tgt_data[collection_idx][k] = -1;
        }

        yac_interpolation_execute_get(interpolation, tgt_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          for (size_t j = 0; j < tgt_grid_data->num_cells; ++j) {
            if (tgt_grid_data->core_cell_mask[j]) {
              if (fabs((ref_tgt_field[tgt_grid_data->cell_ids[j]]) -
                       tgt_data[collection_idx][j]) > 1e-3)
                PUT_ERR("wrong interpolation result");
            } else {
              if (tgt_data[collection_idx][j] != -1.0)
                PUT_ERR("wrong interpolation result");
            }
          }
        }

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          free(tgt_data[collection_idx]);
        free(tgt_data);
      }

      yac_interpolation_delete(interpolation);
    }
  }
#undef NUM_TESTS

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test9() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 8) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 8;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 4;
  int is_target = comm_rank < 4;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // corner and cell ids for a 6 x 5 grid
    //
    // 35-----36-----37-----38-----39-----40-----41
    //  |      |      |      |      |      |      |
    //  |  24  |  25  |  26  |  27  |  28  |  29  |
    //  |      |      |      |      |      |      |
    // 28-----29-----30-----31-----32-----33-----34
    //  |      |      |      |      |      |      |
    //  |  18  |  19  |  20  |  21  |  22  |  23  |
    //  |      |      |      |      |      |      |
    // 21-----22-----23-----24-----25-----26-----27
    //  |      |      |      |      |      |      |
    //  |  12  |  13  |  14  |  15  |  16  |  17  |
    //  |      |      |      |      |      |      |
    // 14-----15-----16-----17-----18-----19-----20
    //  |      |      |      |      |      |      |
    //  |  06  |  07  |  08  |  09  |  10  |  11  |
    //  |      |      |      |      |      |      |
    // 07-----08-----09-----10-----11-----12-----13
    //  |      |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |  05  |
    //  |      |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05-----06
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // +---+---+---+---+---+---+
    // | 2 | 2 | 2 | 3 | 3 | 3 |
    // +---+---+---+---+---+---+
    // | 2 | 2 | 2 | 3 | 3 | 3 |
    // +---+---+---+---+---+---+
    // | 2 | 2 | 2 | 3 | 3 | 3 |
    // +---+---+---+---+---+---+
    // | 0 | 0 | 0 | 1 | 1 | 1 |
    // +---+---+---+---+---+---+
    // | 0 | 0 | 0 | 1 | 1 | 1 |
    // +---+---+---+---+---+---+

    double coordinates_x[] = {-3,-2,-1,0,1,2,3};
    double coordinates_y[] = {-3,-2,-1,0,1,2};
    size_t const num_cells[2] = {6,5};
    size_t local_start[4][2] = {{0,0},{3,0},{0,2},{3,2}};
    size_t local_count[4][2] = {{3,2},{3,2},{3,3},{3,3}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {
    // corner and cell ids for a 5 x 6 grid
    //
    // 36-----37-----38-----39-----40-----41
    //  |      |      |      |      |      |
    //  |  25  |  26  |  27  |  28  |  29  |
    //  |      |      |      |      |      |
    // 30-----31-----32-----33-----34-----35
    //  |      |      |      |      |      |
    //  |  20  |  21  |  22  |  23  |  24  |
    //  |      |      |      |      |      |
    // 24-----25-----26-----27-----28-----29
    //  |      |      |      |      |      |
    //  |  15  |  16  |  17  |  18  |  19  |
    //  |      |      |      |      |      |
    // 18-----19-----20-----21-----22-----23
    //  |      |      |      |      |      |
    //  |  10  |  11  |  12  |  13  |  14  |
    //  |      |      |      |      |      |
    // 12-----13-----14-----15-----16-----17
    //  |      |      |      |      |      |
    //  |  05  |  06  |  07  |  08  |  09  |
    //  |      |      |      |      |      |
    // 06-----07-----08-----09-----10-----11
    //  |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |
    //  |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // +---+---+---+---+---+
    // | 2 | 2 | 3 | 3 | 3 |
    // +---+---+---+---+---+
    // | 2 | 2 | 3 | 3 | 3 |
    // +---+---+---+---+---+
    // | 2 | 2 | 3 | 3 | 3 |
    // +---+---+---+---+---+
    // | 0 | 0 | 0 | 1 | 1 |
    // +---+---+---+---+---+
    // | 0 | 0 | 0 | 1 | 1 |
    // +---+---+---+---+---+
    // | 0 | 0 | 0 | 1 | 1 |
    // +---+---+---+---+---+

    double coordinates_x[] = {-2.5,-1.5,-0.5,0.5,1.5,2.5};
    double coordinates_y[] = {-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5};
    size_t const num_cells[2] = {5,6};
    size_t local_start[4][2] = {{0,0},{3,0},{0,3},{2,3}};
    size_t local_count[4][2] = {{3,2},{3,2},{3,3},{3,3}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name, yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  int order = 1;
  int enforced_conserv = 0;
  int partial_coverage = 0;
  enum yac_interp_method_conserv_normalisation normalisation =
    YAC_INTERP_CONSERV_DESTAREA;

  struct interp_method * method_stack[] = {
    yac_interp_method_conserv_new(
      order, enforced_conserv, partial_coverage, normalisation), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  enum yac_interp_weights_reorder_type reorder_type[2] =
    {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

  for (size_t j = 0; j < 2; ++j) {
    for (size_t collection_size = 1; collection_size < 4;
         collection_size += 2) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      if (is_source) {
        double *** src_data = xmalloc(collection_size * sizeof(*src_data));
        for (size_t i = 0; i < collection_size; ++i)
          src_data[i] = xmalloc(1 * sizeof(**src_data));

        struct yac_basic_grid_data * src_grid_data =
          yac_basic_grid_get_data(src_grid);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          src_data[collection_idx][0] =
            xmalloc(src_grid_data->num_cells * sizeof(***src_data));
          for (size_t i = 0; i < src_grid_data->num_cells; ++i)
            src_data[collection_idx][0][i] =
              (src_grid_data->core_cell_mask[i])?
                (double)(src_grid_data->cell_ids[i]) +
                (double)(collection_idx * 30):
                -1.0;
        }

        yac_interpolation_execute_put(interpolation, src_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          free(src_data[collection_idx][0]);
          free(src_data[collection_idx]);
        }
        free(src_data);
      }
      if (is_target) {
        double ref_tgt_field[30] = {
          (0+1+6+7)/4.0,(1+2+7+8)/4.0,(2+3+8+9)/4.0,
            (3+4+9+10)/4.0,(4+5+10+11)/4.0,
          (6+7+12+13)/4.0,(7+8+13+14)/4.0,(8+9+14+15)/4.0,
            (9+10+15+16)/4.0,(10+11+16+17)/4.0,
          (12+13+18+19)/4.0,(13+14+19+20)/4.0,(14+15+20+21)/4.0,
            (15+16+21+22)/4.0,(16+17+22+23)/4.0,
          (18+19+24+25)/4.0,(19+20+25+26)/4.0,(20+21+26+27)/4.0,
            (21+22+27+28)/4.0,(22+23+28+29)/4.0,
          -1,-1,-1,-1,-1,
          -1,-1,-1,-1,-1,};

        struct yac_basic_grid_data * tgt_grid_data =
          yac_basic_grid_get_data(tgt_grid);

        double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          tgt_data[collection_idx] =
            xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
          for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
            tgt_data[collection_idx][k] = -1;
        }

        yac_interpolation_execute_get(interpolation, tgt_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          for (size_t j = 0, offset = collection_idx * 30;
               j < tgt_grid_data->num_cells; ++j) {
            if (tgt_grid_data->core_cell_mask[j] &&
                (ref_tgt_field[tgt_grid_data->cell_ids[j]] != -1.0)) {
              if (fabs((ref_tgt_field[tgt_grid_data->cell_ids[j]] +
                        (double)offset) - tgt_data[collection_idx][j]) > 1e-3)
                PUT_ERR("wrong interpolation result");
            } else {
              if (tgt_data[collection_idx][j] != -1.0)
                PUT_ERR("wrong interpolation result");
            }
          }
        }

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          free(tgt_data[collection_idx]);
        free(tgt_data);
      }

      yac_interpolation_delete(interpolation);
    }
  }

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test10() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 2) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 2;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank == 1;
  int is_target = comm_rank == 0;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 3x3 grid:
    //
    //       12--21--13--22--14--23--15
    //       |       |       |       |
    //       15  06  17  07  19  08  20
    //       |       |       |       |
    //       08--14--09--16--10--18--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03

    double coordinates_x[] = {0.0,1.0,2.0,3.0};
    double coordinates_y[] = {0.0,1.0,2.0,3.0};
    size_t const num_cells[2] = {3,3};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{3,3}};
    int cell_mask[] = {0,0,0, 0,0,0, 0,0,1};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
    yac_basic_grid_add_mask(src_grid, YAC_LOC_CELL, cell_mask, 9, NULL);
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    double coordinates_x[] = {0.0,0.5,1.0};
    double coordinates_y[] = {0.0,0.5,1.0};
    size_t const num_cells[2] = {2,2};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{2,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  for (int config = 0; config < (1 << 4); ++config) {

    int order = ((config >> 0) & 1) + 1;
    int enforced_conserv = (config >> 1) & 1;
    int partial_coverage = (config >> 2) & 1;
    enum yac_interp_method_conserv_normalisation normalisation =
      (enum yac_interp_method_conserv_normalisation)((config >> 3) & 1);

    if ((order == 2) && (enforced_conserv == 1)) continue;

    struct interp_method * method_stack[] = {
      yac_interp_method_conserv_new(
        order, enforced_conserv, partial_coverage, normalisation),
      yac_interp_method_fixed_new(-2.0), NULL};

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum yac_interp_weights_reorder_type reorder_type[2] =
      {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

    for (size_t j = 0; j < 2; ++j) {
      for (size_t collection_size = 1; collection_size < 4;
          collection_size += 2) {

        struct yac_interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        if (is_source) {
          double *** src_data = xmalloc(collection_size * sizeof(*src_data));
          for (size_t i = 0; i < collection_size; ++i)
            src_data[i] = xmalloc(1 * sizeof(**src_data));

          struct yac_basic_grid_data * src_grid_data =
            yac_basic_grid_get_data(src_grid);

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            double * src_field =
              xmalloc(src_grid_data->num_cells * sizeof(*src_field));
            for (size_t i = 0; i < src_grid_data->num_cells; ++i)
              src_field[i] =
                (src_grid_data->core_cell_mask[i])?
                  (double)(src_grid_data->cell_ids[i] + 1) +
                  (double)(collection_idx * 4):
                  -1.0;
            src_data[collection_idx][0] = src_field;
          }

          yac_interpolation_execute_put(interpolation, src_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            free(src_data[collection_idx][0]);
            free(src_data[collection_idx]);
          }
          free(src_data);
        }
        if (is_target) {
          double ref_tgt_field[4] = {-2.0, -2.0, -2.0, -2.0};

          struct yac_basic_grid_data * tgt_grid_data =
            yac_basic_grid_get_data(tgt_grid);

          double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            tgt_data[collection_idx] =
              xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
            for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
              tgt_data[collection_idx][k] = -1;
          }

          yac_interpolation_execute_get(interpolation, tgt_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            for (size_t j = 0; j < tgt_grid_data->num_cells; ++j) {
              if (tgt_grid_data->core_cell_mask[j] &&
                  (ref_tgt_field[tgt_grid_data->cell_ids[j]] != -1.0)) {
                if (fabs(ref_tgt_field[tgt_grid_data->cell_ids[j]] -
                         tgt_data[collection_idx][j]) > 1e-3)
                  PUT_ERR("wrong interpolation result");
              } else {
                if (tgt_data[collection_idx][j] != -1.0)
                  PUT_ERR("wrong interpolation result");
              }
            }
          }

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx)
            free(tgt_data[collection_idx]);
          free(tgt_data);
        }

        yac_interpolation_delete(interpolation);
      }
    }

    yac_interp_weights_delete(weights);
  }
  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test11() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 2) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 2;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank == 1;
  int is_target = comm_rank == 0;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 3x3 grid:
    //
    //       12--21--13--22--14--23--15
    //       |       |       |       |
    //       15  06  17  07  19  08  20
    //       |       |       |       |
    //       08--14--09--16--10--18--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03

    double coordinates_x[] = {-1.5,-0.5,0.5,1.5};
    double coordinates_y[] = {-1.5,-0.5,0.5,1.5};
    size_t const num_cells[2] = {3,3};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{3,3}};
    int cell_mask[] = {1,1,1, 1,0,1, 1,1,1};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
    yac_basic_grid_add_mask(src_grid, YAC_LOC_CELL, cell_mask, 9, NULL);
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    double coordinates_x[] = {-0.49,0.49};
    double coordinates_y[] = {-0.49,0.49};
    size_t const num_cells[2] = {1,1};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{1,1}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  for (int config = 0; config < (1 << 4); ++config) {

    int order = ((config >> 0) & 1) + 1;
    int enforced_conserv = (config >> 1) & 1;
    int partial_coverage = (config >> 2) & 1;
    enum yac_interp_method_conserv_normalisation normalisation =
      (enum yac_interp_method_conserv_normalisation)((config >> 3) & 1);

    if ((order == 2) && (enforced_conserv == 1)) continue;

    struct interp_method * method_stack[] = {
      yac_interp_method_conserv_new(
        order, enforced_conserv, partial_coverage, normalisation),
      yac_interp_method_fixed_new(-2.0), NULL};

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum yac_interp_weights_reorder_type reorder_type[2] =
      {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

    for (size_t j = 0; j < 2; ++j) {
      for (size_t collection_size = 1; collection_size < 4;
          collection_size += 2) {

        struct yac_interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        if (is_source) {
          double *** src_data = xmalloc(collection_size * sizeof(*src_data));
          for (size_t i = 0; i < collection_size; ++i)
            src_data[i] = xmalloc(1 * sizeof(**src_data));

          struct yac_basic_grid_data * src_grid_data =
            yac_basic_grid_get_data(src_grid);

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            double * src_field =
              xmalloc(src_grid_data->num_cells * sizeof(*src_field));
            for (size_t i = 0; i < src_grid_data->num_cells; ++i)
              src_field[i] =
                (src_grid_data->core_cell_mask[i])?
                  (double)(src_grid_data->cell_ids[i] + 1) +
                  (double)(collection_idx * 4):
                  -1.0;
            src_data[collection_idx][0] = src_field;
          }

          yac_interpolation_execute_put(interpolation, src_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            free(src_data[collection_idx][0]);
            free(src_data[collection_idx]);
          }
          free(src_data);
        }
        if (is_target) {
          double ref_tgt_field[1] = {-2.0};

          struct yac_basic_grid_data * tgt_grid_data =
            yac_basic_grid_get_data(tgt_grid);

          double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            tgt_data[collection_idx] =
              xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
            for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
              tgt_data[collection_idx][k] = -1;
          }

          yac_interpolation_execute_get(interpolation, tgt_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            for (size_t j = 0; j < tgt_grid_data->num_cells; ++j) {
              if (tgt_grid_data->core_cell_mask[j] &&
                  (ref_tgt_field[tgt_grid_data->cell_ids[j]] != -2.0)) {
                if (fabs(ref_tgt_field[tgt_grid_data->cell_ids[j]] -
                         tgt_data[collection_idx][j]) > 1e-3)
                  PUT_ERR("wrong interpolation result");
              } else {
                if (tgt_data[collection_idx][j] != -2.0)
                  PUT_ERR("wrong interpolation result");
              }
            }
          }

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx)
            free(tgt_data[collection_idx]);
          free(tgt_data);
        }

        yac_interpolation_delete(interpolation);
      }
    }

    yac_interp_weights_delete(weights);
  }
  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test12() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 4) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 4;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 1;
  int is_target = comm_rank < 1;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 3x3 grid:
    //
    //       12--21--13--22--14--23--15
    //       |       |       |       |
    //       15  06  17  07  19  08  20
    //       |       |       |       |
    //       08--14--09--16--10--18--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03

    double coordinates_x[] = {-1.0,0.0,1.5,2.0};
    double coordinates_y[] = {-1.0,0.0,1.5,2.0};
    size_t const num_cells[2] = {3,3};
    size_t local_start[3][2] = {{0,0},{0,1},{0,2}};
    size_t local_count[3][2] = {{3,1},{3,1},{3,1}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    // the global target grid is a 4x4 grid:
    //
    //       20--36--21--37--22--38--23--39--24
    //       |       |       |       |       |
    //       28  12  30  13  32  14  34  15  35
    //       |       |       |       |       |
    //       15--27--16--29--17--31--18--33--19
    //       |       |       |       |       |
    //       19  08  21  09  23  10  25  11  26
    //       |       |       |       |       |
    //       10--18--11--20--12--22--13--24--14
    //       |       |       |       |       |
    //       10  04  12  05  14  06  16  07  17
    //       |       |       |       |       |
    //       05--09--06--11--07--13--08--15--09
    //       |       |       |       |       |
    //       01  00  03  01  05  02  07  03  08
    //       |       |       |       |       |
    //       00--00--01--02--02--04--03--06--04

    double coordinates_x[] = {-3.0,-2.0,-1.0,1.0,1.5};
    double coordinates_y[] = {-3.0,-2.0,-1.0,1.0,1.5};
    size_t const num_cells[2] = {4,4};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{4,4}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  if (is_source) {
    int global_src_mask [] = {1,1,1,
                              1,0,1,
                              1,1,1};

    struct yac_basic_grid_data * src_grid_data =
      yac_basic_grid_get_data(src_grid);

    int * src_cell_mask =
      xmalloc(src_grid_data->num_cells * sizeof(*src_cell_mask));
    for (size_t i = 0; i < src_grid_data->num_cells; ++i)
      src_cell_mask[i] = global_src_mask[src_grid_data->cell_ids[i]];
    yac_basic_grid_add_mask_nocpy(
      src_grid, YAC_LOC_CELL, src_cell_mask, NULL);
  }

  if (is_target) {
    int global_tgt_mask [] = {0,0,1,1,
                              0,1,1,1,
                              1,1,1,1,
                              1,1,1,0};
    double coordinates_x[] = {-2.5,-1.75,-0.75,0.25};
    double coordinates_y[] = {-2.5,-1.75,-0.75,0.25};

    struct yac_basic_grid_data * tgt_grid_data =
      yac_basic_grid_get_data(tgt_grid);

    int * tgt_cell_mask =
      xmalloc(tgt_grid_data->num_cells * sizeof(*tgt_cell_mask));
    yac_coordinate_pointer tgt_cell_coords =
      xmalloc(tgt_grid_data->num_cells * sizeof(*tgt_cell_coords));
    for (size_t i = 0; i < tgt_grid_data->num_cells; ++i) {
      tgt_cell_mask[i] = global_tgt_mask[tgt_grid_data->cell_ids[i]];
      LLtoXYZ_deg(
        coordinates_x[tgt_grid_data->cell_ids[i]%4],
        coordinates_y[tgt_grid_data->cell_ids[i]/4],
        tgt_cell_coords[i]);
    }
    yac_basic_grid_add_mask_nocpy(
      tgt_grid, YAC_LOC_CELL, tgt_cell_mask, NULL);
    yac_basic_grid_add_coordinates_nocpy(
      tgt_grid, YAC_LOC_CELL, tgt_cell_coords);
  }

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  int order = 2;
  int enforced_conserv = 0;
  int partial_coverage = 1;
  enum yac_interp_method_conserv_normalisation normalisation =
    YAC_INTERP_CONSERV_FRACAREA;

  struct interp_method * method_stack[] = {
    yac_interp_method_conserv_new(
      order, enforced_conserv, partial_coverage, normalisation),
    yac_interp_method_fixed_new(-2.0), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  enum yac_interp_weights_reorder_type reorder_type[2] =
    {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

#define NUM_TESTS (4)

  for (size_t t = 0; t < NUM_TESTS; ++t) {
    for (size_t j = 0; j < 2; ++j) {

      size_t collection_size = 1;

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      if (is_source) {
        double global_field_data[NUM_TESTS][9] = {
          {1,1,1, 1,1,1, 1,1,1},
          {0,0,0, 1,1,1, 2,2,2},
          {0,1,2, 0,1,2, 0,1,2},
          {0,1,2, 1,2,3, 2,3,4}};
        double *** src_data = xmalloc(collection_size * sizeof(*src_data));
        for (size_t i = 0; i < collection_size; ++i)
          src_data[i] = xmalloc(1 * sizeof(**src_data));

        struct yac_basic_grid_data * src_grid_data =
          yac_basic_grid_get_data(src_grid);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          src_data[collection_idx][0] =
            xmalloc(src_grid_data->num_cells * sizeof(***src_data));
          for (size_t i = 0; i < src_grid_data->num_cells; ++i)
            src_data[collection_idx][0][i] =
              (src_grid_data->core_cell_mask[i])?
                global_field_data[t][src_grid_data->cell_ids[i]]:-1.0;
        }

        yac_interpolation_execute_put(interpolation, src_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          free(src_data[collection_idx][0]);
          free(src_data[collection_idx]);
        }
        free(src_data);
      }
      if (is_target) {
        double ref_tgt_field[NUM_TESTS][16] = {
          {-1.0,-1.0,-2.0,-2.0,
           -1.0,-2.0,-2.0,-2.0,
           -2.0,-2.0,1.0,1.0,
           -2.0,-2.0,1.0,-1.0},
          {-1.0,-1.0,-2.0,-2.0,
           -1.0,-2.0,-2.0,-2.0,
           -2.0,-2.0,(0.0+0.0+1.0)/3.0,0.0,
           -2.0,-2.0,1.0,-1.0},
          {-1.0,-1.0,-2.0,-2.0,
           -1.0,-2.0,-2.0,-2.0,
           -2.0,-2.0,(0.0+1.0+0.0)/3.0,1.0,
           -2.0,-2.0,0.0,-1.0},
          {-1.0,-1.0,-2.0,-2.0,
           -1.0,-2.0,-2.0,-2.0,
           -2.0,-2.0,(0.0+1.0+1.0)/3.0,1.0,
           -2.0,-2.0,1.0,-1.0}};

        struct yac_basic_grid_data * tgt_grid_data =
          yac_basic_grid_get_data(tgt_grid);

        double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          tgt_data[collection_idx] =
            xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
          for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
            tgt_data[collection_idx][k] = -1;
        }

        yac_interpolation_execute_get(interpolation, tgt_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          for (size_t j = 0; j < tgt_grid_data->num_cells; ++j) {
            if (tgt_grid_data->core_cell_mask[j] &&
                (ref_tgt_field[t][tgt_grid_data->cell_ids[j]] != -1.0)) {
              if (fabs((ref_tgt_field[t][tgt_grid_data->cell_ids[j]]) -
                       tgt_data[collection_idx][j]) > 1e-3)
                PUT_ERR("wrong interpolation result");
            } else {
              if (tgt_data[collection_idx][j] != -1.0)
                PUT_ERR("wrong interpolation result");
            }
          }
        }

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          free(tgt_data[collection_idx]);
        free(tgt_data);
      }

      yac_interpolation_delete(interpolation);
    }
  }
#undef NUM_TESTS

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test13() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 2) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 2;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 1;
  int is_target = comm_rank < 1;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 1x1 grid:
    //
    //       02--03--03
    //       |       |
    //       01  00  02
    //       |       |
    //       00--00--01

    double coordinates_x[] = {-0.1,0.1};
    double coordinates_y[] = {-0.1,0.1};
    size_t const num_cells[2] = {1,1};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{1,1}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    src_grid =
      yac_basic_grid_new(
        src_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    // the global target grid:
    //
    //       10--18--11--19--12--20--13--21--14
    //       |       |       |       |       |
    //       10  04  12  05  14  06  16  07  17
    //       |       |       |       |       |
    //       05--09--06--11--07--13--08--15--09
    //       |       |       |       |       |
    //       01  00  03  01  05  02  07  03  08
    //       |       |       |       |       |
    //       00--00--01--02--02--04--03--06--04

    enum {NBR_VERTICES = 15, NBR_CELLS = 8, NUM_VERTICES_PER_CELL = 4};
    int num_vertices_per_cell[NBR_CELLS] = {NUM_VERTICES_PER_CELL,
                                            NUM_VERTICES_PER_CELL,
                                            NUM_VERTICES_PER_CELL,
                                            NUM_VERTICES_PER_CELL,
                                            NUM_VERTICES_PER_CELL,
                                            NUM_VERTICES_PER_CELL,
                                            NUM_VERTICES_PER_CELL,
                                            NUM_VERTICES_PER_CELL};
    double coordinates_x[NBR_VERTICES] = {-3,-2,0,2,3, -3,-2,0,2,3, -3,-2,0,2,3};
    double coordinates_y[NBR_VERTICES] = {-1,-1,-1,-1,-1, 0,0,0,0,0, 1,1,1,1,1};
    int cell_to_vertex[NBR_CELLS][NUM_VERTICES_PER_CELL] =
      {{0,1,6,5}, {1,2,7,6}, {2,3,8,7}, {3,4,9,8},
       {5,6,11,10}, {6,7,12,11}, {7,8,13,12}, {8,9,14,13}};

    struct yac_basic_grid_data tgt_grid_data =
      yac_generate_basic_grid_data_unstruct_deg(
        NBR_VERTICES, NBR_CELLS, num_vertices_per_cell,
        coordinates_x, coordinates_y, &cell_to_vertex[0][0]);
    tgt_grid_data.cell_ids =
      xmalloc(NBR_CELLS * sizeof(*tgt_grid_data.cell_ids));
    for (size_t i = 0; i < NBR_CELLS; ++i)
      tgt_grid_data.cell_ids[i] = (yac_int)i;
    tgt_grid =
      yac_basic_grid_new(tgt_grid_name, tgt_grid_data);
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  int order = 1;
  int enforced_conserv = 0;
  int partial_coverage = 1;
  enum yac_interp_method_conserv_normalisation normalisation =
    YAC_INTERP_CONSERV_FRACAREA;

  struct interp_method * method_stack[] = {
    yac_interp_method_conserv_new(
      order, enforced_conserv, partial_coverage, normalisation), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  size_t collection_size = 1;
  struct yac_interpolation * interpolation =
    yac_interp_weights_get_interpolation(
      weights, YAC_MAPPING_ON_SRC, collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

  // check generated interpolation
  if (is_source) {
    double *** src_data = xmalloc(collection_size * sizeof(*src_data));
    for (size_t i = 0; i < collection_size; ++i)
      src_data[i] = xmalloc(1 * sizeof(**src_data));

    struct yac_basic_grid_data * src_grid_data =
      yac_basic_grid_get_data(src_grid);

    for (size_t collection_idx = 0; collection_idx < collection_size;
          ++collection_idx) {
      src_data[collection_idx][0] =
        xmalloc(src_grid_data->num_cells * sizeof(***src_data));
      for (size_t i = 0; i < src_grid_data->num_cells; ++i)
        src_data[collection_idx][0][i] =
          (src_grid_data->core_cell_mask[i])?
            (double)(src_grid_data->cell_ids[i] + 1) +
            (double)(collection_idx * 9):
            -1.0;
    }

    yac_interpolation_execute_put(interpolation, src_data);

    for (size_t collection_idx = 0; collection_idx < collection_size;
          ++collection_idx) {
      free(src_data[collection_idx][0]);
      free(src_data[collection_idx]);
    }
    free(src_data);
  }
  if (is_target) {
    double ref_tgt_field[8] = {-1,1,1,-1,-1,1,1,-1};

    struct yac_basic_grid_data * tgt_grid_data =
      yac_basic_grid_get_data(tgt_grid);

    double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
    for (size_t collection_idx = 0; collection_idx < collection_size;
          ++collection_idx) {
      tgt_data[collection_idx] =
        xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
      for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
        tgt_data[collection_idx][k] = -1;
    }

    yac_interpolation_execute_get(interpolation, tgt_data);

    for (size_t collection_idx = 0; collection_idx < collection_size;
          ++collection_idx) {
      for (size_t j = 0, offset = collection_idx * 9;
            j < tgt_grid_data->num_cells; ++j) {
        if (((tgt_grid_data->core_cell_mask == NULL) ||
             tgt_grid_data->core_cell_mask[j]) &&
            (ref_tgt_field[tgt_grid_data->cell_ids[j]] != -1.0)) {
          if (fabs((ref_tgt_field[tgt_grid_data->cell_ids[j]] +
                    (double)offset) - tgt_data[collection_idx][j]) > 1e-3)
            PUT_ERR("wrong interpolation result");
        } else {
          if (tgt_data[collection_idx][j] != -1.0)
            PUT_ERR("wrong interpolation result");
        }
      }
    }

    for (size_t collection_idx = 0; collection_idx < collection_size;
          ++collection_idx)
      free(tgt_data[collection_idx]);
    free(tgt_data);
  }

  yac_interpolation_delete(interpolation);

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}

static void test14() {

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size < 2) {
    PUT_ERR("ERROR: too few processes");
    xt_finalize();
    MPI_Finalize();
    exit(TEST_EXIT_CODE);
  }

  int is_active = comm_rank < 2;

  MPI_Comm global_comm;
  MPI_Comm_split(MPI_COMM_WORLD, is_active, 0, &global_comm);

  if (!is_active) {
    MPI_Comm_free(&global_comm);
    return;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int is_source = comm_rank >= 1;
  int is_target = comm_rank < 1;

  MPI_Comm split_comm;
  MPI_Comm_split(global_comm, is_target, 0, &split_comm);

  MPI_Comm_rank(split_comm, &comm_rank);
  MPI_Comm_size(split_comm, &comm_size);

  struct yac_basic_grid * src_grid = NULL;
  struct yac_basic_grid * tgt_grid = NULL;

  if (is_source) {

    // the global source grid is a 3x3 grid:
    //
    //       12--21--13--22--14--23--15
    //       |       |       |       |
    //       15  06  17  07  19  08  20
    //       |       |       |       |
    //       08--14--09--16--10--18--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03
    //
    // mask for the source grid (cells with "x" are masked)
    //
    //       +---+---+---+
    //       | x |   | x |
    //       +---+---+---+
    //       |   |   |   |
    //       +---+---+---+
    //       |   |   |   |
    //       +---+---+---+

    double coordinates_x[] = {-1.5,-0.5,0.5,1.5};
    double coordinates_y[] = {-1.5,-0.5,0.5,1.5};
    size_t const num_cells[2] = {3,3};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{3,3}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    int cell_core_mask[] = {1,1,1,
                            1,1,1,
                            1,0,1};
    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo);
    for (size_t i = 0; i < grid_data.num_cells; ++i)
      grid_data.core_cell_mask[i] =
        grid_data.core_cell_mask[i] &&
        cell_core_mask[grid_data.cell_ids[i]];
    src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
  } else src_grid = yac_basic_grid_empty_new(src_grid_name);

  if (is_target) {

    // the global target grid is a 2x2 grid:
    //
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02

    double coordinates_x[] = {-2,0,2};
    double coordinates_y[] = {-2,0,2};
    size_t const num_cells[2] = {2,2};
    size_t local_start[1][2] = {{0,0}};
    size_t local_count[1][2] = {{2,2}};
    int with_halo = 1;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    tgt_grid =
      yac_basic_grid_new(
        tgt_grid_name,
        yac_generate_basic_grid_data_reg2d(
          coordinates_x, coordinates_y, num_cells,
          local_start[comm_rank], local_count[comm_rank], with_halo));
  } else tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, global_comm);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  int order = 2;
  int enforced_conserv = 0;
  int partial_coverage = 1;
  enum yac_interp_method_conserv_normalisation normalisation =
    YAC_INTERP_CONSERV_FRACAREA;

  struct interp_method * method_stack[] = {
    yac_interp_method_conserv_new(
      order, enforced_conserv, partial_coverage, normalisation), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  size_t collection_size = 1;
  struct yac_interpolation * interpolation =
    yac_interp_weights_get_interpolation(
      weights, YAC_MAPPING_ON_SRC, collection_size, YAC_FRAC_MASK_NO_VALUE,
      1.0, 0.0);

  // check generated interpolation
  if (is_source) {
    double *** src_data = xmalloc(collection_size * sizeof(*src_data));
    for (size_t i = 0; i < collection_size; ++i)
      src_data[i] = xmalloc(1 * sizeof(**src_data));

    struct yac_basic_grid_data * src_grid_data =
      yac_basic_grid_get_data(src_grid);

    for (size_t collection_idx = 0; collection_idx < collection_size;
          ++collection_idx) {
      src_data[collection_idx][0] =
        xmalloc(src_grid_data->num_cells * sizeof(***src_data));
      for (size_t i = 0; i < src_grid_data->num_cells; ++i)
        src_data[collection_idx][0][i] =
          (src_grid_data->core_cell_mask[i])?
            (1.0 + (double)(collection_idx * 4)):-1.0;
    }

    yac_interpolation_execute_put(interpolation, src_data);

    for (size_t collection_idx = 0; collection_idx < collection_size;
          ++collection_idx) {
      free(src_data[collection_idx][0]);
      free(src_data[collection_idx]);
    }
    free(src_data);
  }
  if (is_target) {
    double ref_tgt_field[4] = {1,1,1,1};

    struct yac_basic_grid_data * tgt_grid_data =
      yac_basic_grid_get_data(tgt_grid);

    double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
    for (size_t collection_idx = 0; collection_idx < collection_size;
          ++collection_idx) {
      tgt_data[collection_idx] =
        xmalloc(tgt_grid_data->num_cells * sizeof(**tgt_data));
      for (size_t k = 0; k < tgt_grid_data->num_cells; ++k)
        tgt_data[collection_idx][k] = -1;
    }

    yac_interpolation_execute_get(interpolation, tgt_data);

    for (size_t collection_idx = 0; collection_idx < collection_size;
          ++collection_idx) {
      for (size_t j = 0, offset = collection_idx * 4;
            j < tgt_grid_data->num_cells; ++j) {
        if (tgt_grid_data->core_cell_mask[j] &&
            (ref_tgt_field[tgt_grid_data->cell_ids[j]] != -1.0)) {
          if (fabs((ref_tgt_field[tgt_grid_data->cell_ids[j]] +
                    (double)offset) - tgt_data[collection_idx][j]) > 1e-3)
            PUT_ERR("wrong interpolation result");
        } else {
          if (tgt_data[collection_idx][j] != -1.0)
            PUT_ERR("wrong interpolation result");
        }
      }
    }

    for (size_t collection_idx = 0; collection_idx < collection_size;
          ++collection_idx)
      free(tgt_data[collection_idx]);
    free(tgt_data);
  }

  yac_interpolation_delete(interpolation);

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  MPI_Comm_free(&split_comm);
  MPI_Comm_free(&global_comm);
}
