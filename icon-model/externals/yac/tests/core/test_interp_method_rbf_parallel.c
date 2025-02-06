// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tests.h"
#include "interp_method.h"
#include "interp_method_nnn.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"
#include "geometry.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

char const src_grid_name[] = "src_grid";
char const tgt_grid_name[] = "tgt_grid";

static void target_main(MPI_Comm target_comm);
static void source_main(MPI_Comm source_comm);

int main (void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size != 2) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  int tgt_flag = comm_rank < 1;

  MPI_Comm split_comm;
  MPI_Comm_split(MPI_COMM_WORLD, tgt_flag, 0, &split_comm);

  if (tgt_flag) target_main(split_comm);
  else          source_main(split_comm);

  MPI_Comm_free(&split_comm);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void source_main(MPI_Comm source_comm) {

  // 1 and 5 nearest neighbours per target point with fixed value fallback

  // corner and cell ids for a 7 x 7 grid (x == target point position)

  // 56-----57-----58-----59-----60-----61-----62-----63
  //  |     x|     x|     x|     x|     x|     x|     x|
  //  |  42  |  43  |  44  |  45  |  46  |  47  |  48  |
  //  |      |      |      |      |      |      |      |
  // 48-----49-----50-----51-----52-----53-----54-----55
  //  |     x|     x|     x|     x|     x|     x|     x|
  //  |  35  |  36  |  37  |  38  |  39  |  40  |  41  |
  //  |      |      |      |      |      |      |      |
  // 40-----41-----42-----43-----44-----45-----46-----47
  //  |     x|     x|     x|     x|     x|     x|     x|
  //  |  28  |  29  |  30  |  31  |  32  |  33  |  34  |
  //  |      |      |      |      |      |      |      |
  // 32-----33-----34-----35-----36-----37-----38-----39
  //  |     x|     x|     x|     x|     x|     x|     x|
  //  |  21  |  22  |  23  |  24  |  25  |  26  |  27  |
  //  |      |      |      |      |      |      |      |
  // 24-----25-----26-----27-----28-----29-----30-----31
  //  |     x|     x|     x|     x|     x|     x|     x|
  //  |  14  |  15  |  16  |  17  |  18  |  19  |  20  |
  //  |      |      |      |      |      |      |      |
  // 16-----17-----18-----19-----20-----21-----22-----23
  //  |     x|     x|     x|     x|     x|     x|     x|
  //  |  07  |  08  |  09  |  10  |  11  |  12  |  13  |
  //  |      |      |      |      |      |      |      |
  // 08-----09-----10-----11-----12-----13-----14-----15
  //  |     x|     x|     x|     x|     x|     x|     x|
  //  |  00  |  01  |  02  |  03  |  04  |  05  |  06  |
  //  |      |      |      |      |      |      |      |
  // 00-----01-----02-----03-----04-----05-----06-----07
  //
  // the grid is distributed among the processes as follows:
  // (index == process)
  //
  // 3---3---3---3---3---3---3---3
  // | 3 | 3 | 3 | 3 | 3 | 3 | 3 |
  // 3---3---3---3---3---3---3---3
  // | 3 | 3 | 3 | 3 | 3 | 3 | 3 |
  // 3---3---3---1---2---2---3---3
  // | 1 | 1 | 1 | 2 | 2 | 2 | 2 |
  // 1---1---1---2---2---2---2---2
  // | 1 | 1 | 1 | 2 | 2 | 2 | 2 |
  // 1---1---1---1---2---2---2---2
  // | 1 | 1 | 1 | 2 | 2 | 2 | 2 |
  // 1---1---1---0---0---0---2---2
  // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
  // 0---0---0---0---0---0---0---0
  // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
  // 0---0---0---0---0---0---0---0
  //
  // the source mask looks as follows (# == masked point)
  //
  // +---+---+---+---+---+---+---#
  // |   |   |   |   |   |   |   |
  // +---+---+---+---+---+---#---+
  // |   |   |   |   |   |   |   |
  // +---+---+---+---+---#---+---+
  // |   |   |   |   |   |   |   |
  // +---+---+---+---#---+---+---+
  // |   |   |   |   |   |   |   |
  // +---+---+---#---+---+---+---+
  // |   |   |   |   |   |   |   |
  // +---+---#---+---+---+---+---+
  // |   |   |   |   |   |   |   |
  // +---#---+---+---+---+---+---+
  // |   |   |   |   |   |   |   |
  // #---+---+---+---+---+---+---+

  int my_source_rank;
  MPI_Comm_rank(source_comm, &my_source_rank);

  double coordinates_x[] = {-2.0, -1.0, 0.0, 1.0, 2.0};
  double coordinates_y[] = {-2.0, -1.0, 0.0, 1.0, 2.0};
  size_t const num_cells[2] = {4,4};
  size_t local_start[2] = {0,0};
  size_t local_count[2] = {4,4};
  int with_halo = 0;
  for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
  for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

  struct yac_basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg2d(
      coordinates_x, coordinates_y, num_cells,
      local_start, local_count, with_halo);

  struct yac_basic_grid * grids[2] =
    {yac_basic_grid_new(src_grid_name, grid_data),
     yac_basic_grid_empty_new(tgt_grid_name)};
  yac_basic_grid_add_coordinates_nocpy(grids[0], YAC_LOC_CELL, NULL);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  struct interp_method * method_stack[2] =
    {yac_interp_method_nnn_new(
       (struct yac_nnn_config){
         .type = YAC_INTERP_NNN_RBF, .n = 9,
          .max_search_distance = YAC_INTERP_RBF_MAX_SEARCH_DISTANCE_DEFAULT,
         .data.rbf_scale = YAC_INTERP_RBF_SCALE_DEFAULT}),
     NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  enum yac_interp_weights_reorder_type reorder_type[2] =
    {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

  for (size_t i = 0; i < 2; ++i) {
    for (size_t collection_size = 1; collection_size < 4;
         collection_size += 2) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      {
        double *** src_data = xmalloc(collection_size * sizeof(*src_data));
        double ref_src_data[25] = {0,1,2,3,4,
                                   1,2,3,4,5,
                                   2,3,4,5,6,
                                   3,4,5,6,7,
                                   4,5,6,7,8};

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          // only one field
          src_data[collection_idx] = xmalloc(1 * sizeof(**src_data));
          src_data[collection_idx][0] =
            xmalloc(grid_data.num_vertices * sizeof(***src_data));
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
            src_data[collection_idx][0][i] =
              ref_src_data[i] + (double)(collection_idx * 9);
        }

        yac_interpolation_execute(interpolation, src_data, NULL);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          free(src_data[collection_idx][0]);
          free(src_data[collection_idx]);
        }
        free(src_data);
      }

      yac_interpolation_delete(interpolation);
    }
  }

  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);
  yac_basic_grid_delete(grids[1]);
  yac_basic_grid_delete(grids[0]);
}

static void target_main(MPI_Comm target_comm) {

  UNUSED(target_comm);

  double coordinates_x[] = {-1.5,-0.5,0.5,1.5};
  double coordinates_y[] = {-1.5,-0.5,0.5,1.5};
  double cell_coordinates_x[] = {-1.0,0.0,1.0};
  double cell_coordinates_y[] = {-1.0,0.0,1.0};
  yac_coordinate_pointer cell_coordinates = xmalloc(9 * sizeof(*cell_coordinates));
  size_t const num_cells[2] = {3,3};
  size_t local_start[2] = {0,0};
  size_t local_count[2] = {3,3};
  int with_halo = 0;
  for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
  for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;
  for (size_t i = 0, k = 0; i < num_cells[1]; ++i)
    for (size_t j = 0; j < num_cells[0]; ++j, ++k)
      LLtoXYZ_deg(
        cell_coordinates_x[j], cell_coordinates_y[i], cell_coordinates[k]);

  struct yac_basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg2d(
      coordinates_x, coordinates_y, num_cells,
      local_start, local_count, with_halo);

  struct yac_basic_grid * grids[2] =
    {yac_basic_grid_new(tgt_grid_name, grid_data),
     yac_basic_grid_empty_new(src_grid_name)};
  yac_basic_grid_add_coordinates_nocpy(
    grids[0], YAC_LOC_CELL, cell_coordinates);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                        num_src_fields, src_fields, tgt_field);

  struct interp_method * method_stack[2] =
    {yac_interp_method_nnn_new(
       (struct yac_nnn_config){
         .type = YAC_INTERP_NNN_RBF, .n = 9,
          .max_search_distance = YAC_INTERP_RBF_MAX_SEARCH_DISTANCE_DEFAULT,
         .data.rbf_scale = YAC_INTERP_RBF_SCALE_DEFAULT}),
     NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  enum yac_interp_weights_reorder_type reorder_type[2] =
    {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

  for (size_t i = 0; i < 2; ++i) {
    for (size_t collection_size = 1; collection_size < 4;
         collection_size += 2) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      {
        double ref_tgt_field[9] = {2,3,4, 3,4,5, 4,5,6};

        double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx) {
          tgt_data[collection_idx] =
            xmalloc(grid_data.num_cells * sizeof(**tgt_data));
          for (size_t j = 0; j < 9; ++j) tgt_data[collection_idx][j] = -1.0;
        }

        yac_interpolation_execute(interpolation, NULL, tgt_data);

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          for (size_t j = 0, offset = collection_idx * 9;
               j < grid_data.num_cells; ++j)
            if (fabs((ref_tgt_field[j] + (double)offset) -
                     tgt_data[collection_idx][j]) > 1e-3)
              PUT_ERR("wrong interpolation result");

        for (size_t collection_idx = 0; collection_idx < collection_size;
             ++collection_idx)
          free(tgt_data[collection_idx]);
        free(tgt_data);
      }

      yac_interpolation_delete(interpolation);
    }
  }

  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);
  yac_basic_grid_delete(grids[1]);
  yac_basic_grid_delete(grids[0]);
}
