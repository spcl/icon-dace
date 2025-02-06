// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tests.h"
#include "interp_method.h"
#include "interp_method_nnn.h"
#include "interp_method_fixed.h"
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

  if (comm_size != 5) {
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

  { // 1 and 5 nearest neighbours per target point with fixed value fallback

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

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    double coordinates_y[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[4][2] = {{0,0},{0,2},{3,2},{0,5}};
    size_t local_count[4][2] = {{7,2},{3,3},{4,3},{7,2}};
    int with_halo = 0;
    int global_corner_mask[8][8] = {
      {0,1,1,1,1,1,1,1},
      {1,0,1,1,1,1,1,1},
      {1,1,0,1,1,1,1,1},
      {1,1,1,0,1,1,1,1},
      {1,1,1,1,0,1,1,1},
      {1,1,1,1,1,0,1,1},
      {1,1,1,1,1,1,0,1},
      {1,1,1,1,1,1,1,0}};
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[my_source_rank], local_count[my_source_rank], with_halo);

    int * src_corner_mask =
      xmalloc(grid_data.num_vertices * sizeof(*src_corner_mask));
    for (size_t i = 0; i < grid_data.num_vertices; ++i)
      src_corner_mask[i] =
        ((int*)(&(global_corner_mask[0][0])))[grid_data.vertex_ids[i]];

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
    yac_basic_grid_add_mask_nocpy(
      src_grid, YAC_LOC_CORNER, src_corner_mask, NULL);
    struct yac_basic_grid * tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    size_t num_stacks = 7;
    struct interp_method * method_stack[7][3] = {
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_AVG, .n = 1,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_AVG, .n = 5,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_DIST, .n = 5,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_GAUSS, .n = 5,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
          .data.gauss_scale = YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_ZERO, .n = 1,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_AVG, .n = 1,
          .max_search_distance = 0.4 * YAC_RAD}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_AVG, .n = 5, .max_search_distance = 1.3 * YAC_RAD}),
       yac_interp_method_fixed_new(-1), NULL}};

    struct yac_interp_weights * weights[7];
    for (size_t i = 0; i < num_stacks; ++i) {
      weights[i] = yac_interp_method_do_search(method_stack[i], interp_grid);
      yac_interp_method_delete(method_stack[i]);
    }

    enum yac_interp_weights_reorder_type reorder_type[2] =
      {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

    for (size_t i = 0; i < num_stacks; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        for (size_t collection_size = 1; collection_size < 4;
            collection_size += 2) {

          struct yac_interpolation * interpolation =
            yac_interp_weights_get_interpolation(
              weights[i], reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

          // check generated interpolation
          {
            double *** src_data = xmalloc(collection_size * sizeof(*src_data));

            for (size_t collection_idx = 0; collection_idx < collection_size;
                ++collection_idx) {
              // only one field
              src_data[collection_idx] = xmalloc(1 * sizeof(**src_data));
              src_data[collection_idx][0] =
                xmalloc(grid_data.num_vertices * sizeof(***src_data));
              for (size_t i = 0; i < grid_data.num_vertices; ++i)
                src_data[collection_idx][0][i] =
                  (double)(grid_data.vertex_ids[i]) +
                  (double)(collection_idx * 49);
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
    }

    for (size_t i = 0; i < num_stacks; ++i)
      yac_interp_weights_delete(weights[i]);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);
  }

  { // 1 and 5 nearest neighbours per target point with fixed value fallback

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
    // +---+---#---#---#---#---#---#
    // |   |   |   |   |   |   |   |
    // +---+---#---#---#---#---#---#
    // |   |   |   |   |   |   |   |
    // #---#---#---#---#---#---#---#
    // |   |   |   |   |   |   |   |
    // #---#---#---#---#---#---#---#
    // |   |   |   |   |   |   |   |
    // #---#---#---#---#---#---#---#
    // |   |   |   |   |   |   |   |
    // #---#---#---#---#---#---#---#
    // |   |   |   |   |   |   |   |
    // +---#---#---#---#---#---#---+
    // |   |   |   |   |   |   |   |
    // +---+---#---#---#---#---+---+
    //
    // the target mask looks as follows (# == masked point)
    //
    // +---+---+---+---+---+---+---+
    // |   |   |   |   |   |   |   |
    // +---+---+---+---+---+---+---+
    // |   |   |   |   |   |   |   |
    // +---+---#---#---#---#---+---+
    // |   |   |   |   |   |   |   |
    // +---+---#---#---#---#---+---+
    // |   |   |   |   |   |   |   |
    // +---+---#---#---#---#---+---+
    // |   |   |   |   |   |   |   |
    // +---+---#---#---#---#---+---+
    // |   |   |   |   |   |   |   |
    // +---+---+---+---+---+---+---+
    // |   |   |   |   |   |   |   |
    // +---+---+---+---+---+---+---+
    //--------------------
    // set up global search
    //--------------------

    int my_source_rank;
    MPI_Comm_rank(source_comm, &my_source_rank);

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    double coordinates_y[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[4][2] = {{0,0},{0,2},{3,2},{0,5}};
    size_t local_count[4][2] = {{7,2},{3,3},{4,3},{7,2}};
    int with_halo = 0;
    int global_corner_mask[8][8] = {
      {1,1,0,0,0,0,1,1},
      {1,0,0,0,0,0,0,1},
      {0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0},
      {1,1,0,0,0,0,0,0},
      {1,1,0,0,0,0,0,0}};
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[my_source_rank], local_count[my_source_rank], with_halo);

    int * src_corner_mask =
      xmalloc(grid_data.num_vertices * sizeof(*src_corner_mask));
    for (size_t i = 0; i < grid_data.num_vertices; ++i)
      src_corner_mask[i] =
        ((int*)(&(global_corner_mask[0][0])))[grid_data.vertex_ids[i]];

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
    yac_basic_grid_add_mask_nocpy(
      src_grid, YAC_LOC_CORNER, src_corner_mask, NULL);
    struct yac_basic_grid * tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[2] = {
      yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_AVG, .n = 4,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
      NULL};

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum yac_interp_weights_reorder_type reorder_type[2] =
      {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

    for (size_t i = 0; i < 2; ++i) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      {
        double * src_field = NULL;
        double ** src_fields = &src_field;

        src_field = xmalloc(grid_data.num_vertices * sizeof(*src_field));
        for (size_t i = 0; i < grid_data.num_vertices; ++i)
          src_field[i] = (double)(grid_data.vertex_ids[i]);

        yac_interpolation_execute(interpolation, &src_fields, NULL);

        free(src_field);
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);
  }

  { // 4 nearest neighbours per target point with fixed value fallback
    // with matching vertices between source and target

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0};
    double coordinates_y[] = {0.0, 1.0, 2.0, 3.0};
    size_t const num_cells[2] = {3,3};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {3,3};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start, local_count, with_halo);

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
    struct yac_basic_grid * tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    enum {NUM_STACKS = 2};
    struct interp_method * method_stack[NUM_STACKS][3] = {
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_DIST, .n = 4,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
      yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_GAUSS, .n = 4,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
          .data.gauss_scale = YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT}),
      yac_interp_method_fixed_new(-1), NULL}};

    struct yac_interp_weights * weights[NUM_STACKS];
    for (size_t i = 0; i < NUM_STACKS; ++i) {
      weights[i] = yac_interp_method_do_search(method_stack[i], interp_grid);
      yac_interp_method_delete(method_stack[i]);
    }

    for (size_t i = 0; i < NUM_STACKS; ++i) {
      for (size_t collection_size = 1; collection_size < 4;
          collection_size += 2) {

        struct yac_interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights[i], YAC_MAPPING_ON_SRC, collection_size,
            YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        {
          double *** src_data = xmalloc(collection_size * sizeof(*src_data));

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            // only one field
            src_data[collection_idx] = xmalloc(1 * sizeof(**src_data));
            src_data[collection_idx][0] =
              xmalloc(grid_data.num_vertices * sizeof(***src_data));
            for (size_t i = 0; i < grid_data.num_vertices; ++i)
              src_data[collection_idx][0][i] =
                (double)(grid_data.vertex_ids[i]) +
                (double)(collection_idx * 16);
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

    for (size_t i = 0; i < NUM_STACKS; ++i)
      yac_interp_weights_delete(weights[i]);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);
  }

  { // 4 nearest neighbours per target point
    // where the 4 source points are far away from the target point

    double coordinates_x[] = {0.0, 1.0, 2.0, 3.0};
    double coordinates_y[] = {0.0, 1.0, 2.0, 3.0};
    size_t const num_cells[2] = {3,3};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {3,3};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start, local_count, with_halo);

    struct yac_basic_grid * src_grid =
      yac_basic_grid_new(src_grid_name, grid_data);
    struct yac_basic_grid * tgt_grid = yac_basic_grid_empty_new(tgt_grid_name);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    enum {NUM_STACKS = 2};
    struct interp_method * method_stack[NUM_STACKS][3] = {
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_DIST, .n = 4,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
      yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_GAUSS, .n = 4,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
          .data.gauss_scale = YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT}),
      yac_interp_method_fixed_new(-1), NULL}};

    struct yac_interp_weights * weights[NUM_STACKS];
    for (size_t i = 0; i < NUM_STACKS; ++i) {
      weights[i] = yac_interp_method_do_search(method_stack[i], interp_grid);
      yac_interp_method_delete(method_stack[i]);
    }

    for (size_t i = 0; i < NUM_STACKS; ++i) {
      for (size_t collection_size = 1; collection_size < 4;
          collection_size += 2) {

        struct yac_interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights[i], YAC_MAPPING_ON_SRC, collection_size,
            YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        {
          double *** src_data = xmalloc(collection_size * sizeof(*src_data));

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            // only one field
            src_data[collection_idx] = xmalloc(1 * sizeof(**src_data));
            src_data[collection_idx][0] =
              xmalloc(grid_data.num_vertices * sizeof(***src_data));
            for (size_t i = 0; i < grid_data.num_vertices; ++i)
              src_data[collection_idx][0][i] =
                (double)(grid_data.vertex_ids[i]) +
                (double)(collection_idx * 16);
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

    for (size_t i = 0; i < NUM_STACKS; ++i)
      yac_interp_weights_delete(weights[i]);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);
  }
}

static double compute_dist_result_(
  size_t n, size_t * src_indices, size_t tgt_index,
  double * src_coordinates_x, double * src_coordinates_y, size_t src_size_x,
  double * tgt_coordinates_x, double * tgt_coordinates_y, size_t tgt_size_x) {

  double tgt_coord[3];
  LLtoXYZ_deg(
    tgt_coordinates_x[tgt_index%tgt_size_x],
    tgt_coordinates_y[tgt_index/tgt_size_x], tgt_coord);

  double result = 0.0;
  double weights_sum = 0.0;

  for (size_t i = 0; i < n; ++i) {

    double src_coord[3];
    LLtoXYZ_deg(src_coordinates_x[src_indices[i]%src_size_x],
                src_coordinates_y[src_indices[i]/src_size_x], src_coord);

    double weight = 1.0 / get_vector_angle(tgt_coord, src_coord);
    result += weight * (double)src_indices[i];
    weights_sum += weight;
  }
  return result / weights_sum;
}

static double compute_dist_result(
  size_t n, size_t * src_indices, size_t tgt_index) {

  return
    compute_dist_result_(
      n, src_indices, tgt_index,
      (double[]){0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},
      (double[]){0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0}, 8,
      (double[]){0.75,1.75,2.75,3.75,4.75,5.75,6.75},
      (double[]){0.75,1.75,2.75,3.75,4.75,5.75,6.75}, 7);
}

static double compute_dist_result_2(
  size_t n, size_t * src_indices, size_t tgt_index) {

  return
    compute_dist_result_(
      n, src_indices, tgt_index,
      (double[]){0.0,1.0,2.0,3.0},
      (double[]){0.0,1.0,2.0,3.0}, 4,
      (double[]){1.0,1.5,2.0},
      (double[]){1.0,1.5,2.0}, 3);
}

static double compute_dist_result_3(
  size_t n, size_t * src_indices, size_t tgt_index) {

  return
    compute_dist_result_(
      n, src_indices, tgt_index,
      (double[]){0.0,1.0,2.0,3.0},
      (double[]){0.0,1.0,2.0,3.0}, 4,
      (double[]){1.0+90.0,1.5+90.0,2.0+90.0},
      (double[]){1.0,1.5,2.0}, 3);
}

static double compute_gauss_result_(
  size_t n, size_t * src_indices, size_t tgt_index,
  double * src_coordinates_x, double * src_coordinates_y, size_t src_size_x,
  double * tgt_coordinates_x, double * tgt_coordinates_y, size_t tgt_size_x) {

  double src_coords[n][3];

  for (size_t i = 0; i < n; ++i)
    LLtoXYZ_deg(src_coordinates_x[src_indices[i]%src_size_x],
                src_coordinates_y[src_indices[i]/src_size_x], src_coords[i]);

  // "- n" because we do not count the diagonal
  double src_distances_sum = 0.0;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      src_distances_sum += get_vector_angle(src_coords[i], src_coords[j]);

  size_t src_distance_count = n * n - n;

  // compute mean distance
  double src_distance_mean = src_distances_sum / (double)src_distance_count;
  double src_distance_mean_sq = src_distance_mean * src_distance_mean;

  // compute weights
  double tgt_coord[3];
  LLtoXYZ_deg(tgt_coordinates_x[tgt_index%tgt_size_x],
              tgt_coordinates_y[tgt_index/tgt_size_x], tgt_coord);
  double weights[n];
  for (size_t i = 0; i < n; ++i) {
    double tgt_distance = get_vector_angle(src_coords[i], tgt_coord);
    weights[i] =
      exp(-1.0 * (tgt_distance * tgt_distance) /
                 (YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT * src_distance_mean_sq));
  }

  // compute sum of weights
  double weights_sum = 0.0;
  for (size_t i = 0; i < n; ++i) weights_sum += weights[i];
  for (size_t i = 0; i < n; ++i) weights[i] /= weights_sum;

  // compute interpolation results
  double result = 0.0;
  for (size_t i = 0; i < n; ++i)
    result += weights[i] * (double)src_indices[i];

  return result;
}

static double compute_gauss_result(
  size_t n, size_t * src_indices, size_t tgt_index) {

  return
    compute_gauss_result_(
      n, src_indices, tgt_index,
      (double[]){0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},
      (double[]){0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0}, 8,
      (double[]){0.75,1.75,2.75,3.75,4.75,5.75,6.75},
      (double[]){0.75,1.75,2.75,3.75,4.75,5.75,6.75}, 7);
}

static double compute_gauss_result_2(
  size_t n, size_t * src_indices, size_t tgt_index) {

  return
    compute_gauss_result_(
      n, src_indices, tgt_index,
      (double[]){0.0,1.0,2.0,3.0},
      (double[]){0.0,1.0,2.0,3.0}, 4,
      (double[]){1.0,1.5,2.0},
      (double[]){1.0,1.5,2.0}, 3);
}

static void target_main(MPI_Comm target_comm) {

  UNUSED(target_comm);

  {
    double coordinates_x[] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double cell_coordinates_x[] = {0.75,1.75,2.75,3.75,4.75,5.75,6.75};
    double cell_coordinates_y[] = {0.75,1.75,2.75,3.75,4.75,5.75,6.75};
    yac_coordinate_pointer cell_coordinates = xmalloc(49 * sizeof(*cell_coordinates));
    size_t const num_cells[2] = {7,7};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {7,7};
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

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    yac_basic_grid_add_coordinates_nocpy(
      tgt_grid, YAC_LOC_CELL, cell_coordinates);
    struct yac_basic_grid * src_grid = yac_basic_grid_empty_new(src_grid_name);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    size_t num_stacks = 7;
    struct interp_method * method_stack[7][3] = {
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_AVG, .n = 1,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_AVG, .n = 5,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_DIST, .n = 5,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_GAUSS, .n = 5,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
          .data.gauss_scale = YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_ZERO, .n = 1,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_AVG, .n = 1,
          .max_search_distance = 0.4 * YAC_RAD}),
       yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_AVG, .n = 5,
          .max_search_distance = 1.3 * YAC_RAD}),
       yac_interp_method_fixed_new(-1), NULL}};


    struct yac_interp_weights * weights[7];
    for (size_t i = 0; i < num_stacks; ++i) {
      weights[i] = yac_interp_method_do_search(method_stack[i], interp_grid);
      yac_interp_method_delete(method_stack[i]);
    }

    enum yac_interp_weights_reorder_type reorder_type[2] =
      {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

    for (size_t i = 0; i < num_stacks; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        for (size_t collection_size = 1; collection_size < 4;
            collection_size += 2) {

          struct yac_interpolation * interpolation =
            yac_interp_weights_get_interpolation(
              weights[i], reorder_type[j], collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

          // check generated interpolation
          {
            double ref_tgt_field[][49] =
              {{8,10,11,12,13,14,15,
                17,17,19,20,21,22,23,
                25,26,26,28,29,30,31,
                33,34,35,35,37,38,39,
                41,42,43,44,44,46,47,
                49,50,51,52,53,53,55,
                57,58,59,60,61,62,62},
              {(1+8+10+16+17)/5.0,(1+2+10+11+17)/5.0,
                (2+3+10+11+12)/5.0,(3+4+11+12+13)/5.0,
                (4+5+12+13+14)/5.0,(5+6+13+14+15)/5.0,
                (6+7+14+15+23)/5.0,(8+16+17+24+25)/5.0,
                (10+17+19+25+26)/5.0,(10+11+19+20+26)/5.0,
                (11+12+19+20+21)/5.0,(12+13+20+21+22)/5.0,
                (13+14+21+22+23)/5.0,(14+15+22+23+31)/5.0,
                (16+17+24+25+26)/5.0,(17+19+25+26+34)/5.0,
                (19+20+26+28+35)/5.0,(19+20+21+28+29)/5.0,
                (20+21+28+29+30)/5.0,(21+22+29+30+31)/5.0,
                (22+23+30+31+39)/5.0,(24+25+32+33+34)/5.0,
                (25+26+33+34+35)/5.0,(26+28+34+35+43)/5.0,
                (28+29+35+37+44)/5.0,(28+29+30+37+38)/5.0,
                (29+30+37+38+39)/5.0,(30+31+38+39+47)/5.0,
                (32+33+40+41+42)/5.0,(33+34+41+42+43)/5.0,
                (34+35+42+43+44)/5.0,(35+37+43+44+52)/5.0,
                (37+38+44+46+53)/5.0,(37+38+39+46+47)/5.0,
                (38+39+46+47+55)/5.0,(40+41+48+49+50)/5.0,
                (41+42+49+50+51)/5.0,(42+43+50+51+52)/5.0,
                (43+44+51+52+53)/5.0,(44+46+52+53+61)/5.0,
                (46+47+53+55+62)/5.0,(46+47+53+55+62)/5.0,
                (48+49+56+57+58)/5.0,(49+50+57+58+59)/5.0,
                (50+51+58+59+60)/5.0,(51+52+59+60+61)/5.0,
                (52+53+60+61+62)/5.0,(53+55+60+61+62)/5.0,
                (47+53+55+61+62)/5.0},
              {compute_dist_result(5, (size_t[]){1,8,10,16,17}, 0),
                compute_dist_result(5, (size_t[]){1,2,10,11,17}, 1),
                compute_dist_result(5, (size_t[]){2,3,10,11,12}, 2),
                compute_dist_result(5, (size_t[]){3,4,11,12,13}, 3),
                compute_dist_result(5, (size_t[]){4,5,12,13,14}, 4),
                compute_dist_result(5, (size_t[]){5,6,13,14,15}, 5),
                compute_dist_result(5, (size_t[]){6,7,14,15,23}, 6),
                compute_dist_result(5, (size_t[]){8,16,17,24,25}, 7),
                compute_dist_result(5, (size_t[]){10,17,19,25,26}, 8),
                compute_dist_result(5, (size_t[]){10,11,19,20,26}, 9),
                compute_dist_result(5, (size_t[]){11,12,19,20,21}, 10),
                compute_dist_result(5, (size_t[]){12,13,20,21,22}, 11),
                compute_dist_result(5, (size_t[]){13,14,21,22,23}, 12),
                compute_dist_result(5, (size_t[]){14,15,22,23,31}, 13),
                compute_dist_result(5, (size_t[]){16,17,24,25,26}, 14),
                compute_dist_result(5, (size_t[]){17,19,25,26,34}, 15),
                compute_dist_result(5, (size_t[]){19,20,26,28,35}, 16),
                compute_dist_result(5, (size_t[]){19,20,21,28,29}, 17),
                compute_dist_result(5, (size_t[]){20,21,28,29,30}, 18),
                compute_dist_result(5, (size_t[]){21,22,29,30,31}, 19),
                compute_dist_result(5, (size_t[]){22,23,30,31,39}, 20),
                compute_dist_result(5, (size_t[]){24,25,32,33,34}, 21),
                compute_dist_result(5, (size_t[]){25,26,33,34,35}, 22),
                compute_dist_result(5, (size_t[]){26,28,34,35,43}, 23),
                compute_dist_result(5, (size_t[]){28,29,35,37,44}, 24),
                compute_dist_result(5, (size_t[]){28,29,30,37,38}, 25),
                compute_dist_result(5, (size_t[]){29,30,37,38,39}, 26),
                compute_dist_result(5, (size_t[]){30,31,38,39,47}, 27),
                compute_dist_result(5, (size_t[]){32,33,40,41,42}, 28),
                compute_dist_result(5, (size_t[]){33,34,41,42,43}, 29),
                compute_dist_result(5, (size_t[]){34,35,42,43,44}, 30),
                compute_dist_result(5, (size_t[]){35,37,43,44,52}, 31),
                compute_dist_result(5, (size_t[]){37,38,44,46,53}, 32),
                compute_dist_result(5, (size_t[]){37,38,39,46,47}, 33),
                compute_dist_result(5, (size_t[]){38,39,46,47,55}, 34),
                compute_dist_result(5, (size_t[]){40,41,48,49,50}, 35),
                compute_dist_result(5, (size_t[]){41,42,49,50,51}, 36),
                compute_dist_result(5, (size_t[]){42,43,50,51,52}, 37),
                compute_dist_result(5, (size_t[]){43,44,51,52,53}, 38),
                compute_dist_result(5, (size_t[]){44,46,52,53,61}, 39),
                compute_dist_result(5, (size_t[]){46,47,53,55,62}, 40),
                compute_dist_result(5, (size_t[]){46,47,53,55,62}, 41),
                compute_dist_result(5, (size_t[]){48,49,56,57,58}, 42),
                compute_dist_result(5, (size_t[]){49,50,57,58,59}, 43),
                compute_dist_result(5, (size_t[]){50,51,58,59,60}, 44),
                compute_dist_result(5, (size_t[]){51,52,59,60,61}, 45),
                compute_dist_result(5, (size_t[]){52,53,60,61,62}, 46),
                compute_dist_result(5, (size_t[]){53,55,60,61,62}, 47),
                compute_dist_result(5, (size_t[]){47,53,55,61,62}, 48)},
              {compute_gauss_result(5, (size_t[]){1,8,10,16,17}, 0),
                compute_gauss_result(5, (size_t[]){1,2,10,11,17}, 1),
                compute_gauss_result(5, (size_t[]){2,3,10,11,12}, 2),
                compute_gauss_result(5, (size_t[]){3,4,11,12,13}, 3),
                compute_gauss_result(5, (size_t[]){4,5,12,13,14}, 4),
                compute_gauss_result(5, (size_t[]){5,6,13,14,15}, 5),
                compute_gauss_result(5, (size_t[]){6,7,14,15,23}, 6),
                compute_gauss_result(5, (size_t[]){8,16,17,24,25}, 7),
                compute_gauss_result(5, (size_t[]){10,17,19,25,26}, 8),
                compute_gauss_result(5, (size_t[]){10,11,19,20,26}, 9),
                compute_gauss_result(5, (size_t[]){11,12,19,20,21}, 10),
                compute_gauss_result(5, (size_t[]){12,13,20,21,22}, 11),
                compute_gauss_result(5, (size_t[]){13,14,21,22,23}, 12),
                compute_gauss_result(5, (size_t[]){14,15,22,23,31}, 13),
                compute_gauss_result(5, (size_t[]){16,17,24,25,26}, 14),
                compute_gauss_result(5, (size_t[]){17,19,25,26,34}, 15),
                compute_gauss_result(5, (size_t[]){19,20,26,28,35}, 16),
                compute_gauss_result(5, (size_t[]){19,20,21,28,29}, 17),
                compute_gauss_result(5, (size_t[]){20,21,28,29,30}, 18),
                compute_gauss_result(5, (size_t[]){21,22,29,30,31}, 19),
                compute_gauss_result(5, (size_t[]){22,23,30,31,39}, 20),
                compute_gauss_result(5, (size_t[]){24,25,32,33,34}, 21),
                compute_gauss_result(5, (size_t[]){25,26,33,34,35}, 22),
                compute_gauss_result(5, (size_t[]){26,28,34,35,43}, 23),
                compute_gauss_result(5, (size_t[]){28,29,35,37,44}, 24),
                compute_gauss_result(5, (size_t[]){28,29,30,37,38}, 25),
                compute_gauss_result(5, (size_t[]){29,30,37,38,39}, 26),
                compute_gauss_result(5, (size_t[]){30,31,38,39,47}, 27),
                compute_gauss_result(5, (size_t[]){32,33,40,41,42}, 28),
                compute_gauss_result(5, (size_t[]){33,34,41,42,43}, 29),
                compute_gauss_result(5, (size_t[]){34,35,42,43,44}, 30),
                compute_gauss_result(5, (size_t[]){35,37,43,44,52}, 31),
                compute_gauss_result(5, (size_t[]){37,38,44,46,53}, 32),
                compute_gauss_result(5, (size_t[]){37,38,39,46,47}, 33),
                compute_gauss_result(5, (size_t[]){38,39,46,47,55}, 34),
                compute_gauss_result(5, (size_t[]){40,41,48,49,50}, 35),
                compute_gauss_result(5, (size_t[]){41,42,49,50,51}, 36),
                compute_gauss_result(5, (size_t[]){42,43,50,51,52}, 37),
                compute_gauss_result(5, (size_t[]){43,44,51,52,53}, 38),
                compute_gauss_result(5, (size_t[]){44,46,52,53,61}, 39),
                compute_gauss_result(5, (size_t[]){46,47,53,55,62}, 40),
                compute_gauss_result(5, (size_t[]){46,47,53,55,62}, 41),
                compute_gauss_result(5, (size_t[]){48,49,56,57,58}, 42),
                compute_gauss_result(5, (size_t[]){49,50,57,58,59}, 43),
                compute_gauss_result(5, (size_t[]){50,51,58,59,60}, 44),
                compute_gauss_result(5, (size_t[]){51,52,59,60,61}, 45),
                compute_gauss_result(5, (size_t[]){52,53,60,61,62}, 46),
                compute_gauss_result(5, (size_t[]){53,55,60,61,62}, 47),
                compute_gauss_result(5, (size_t[]){47,53,55,61,62}, 48)},
              {0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,
                0,0,0,0,0,0,0},
              {-1,10,11,12,13,14,15,
                17,-1,19,20,21,22,23,
                25,26,-1,28,29,30,31,
                33,34,35,-1,37,38,39,
                41,42,43,44,-1,46,47,
                49,50,51,52,53,-1,55,
                57,58,59,60,61,62,-1},
              {-1.0,-1.0,
                (2+3+10+11+12)/5.0,(3+4+11+12+13)/5.0,
                (4+5+12+13+14)/5.0,(5+6+13+14+15)/5.0,
                (6+7+14+15+23)/5.0,-1.0,
                -1.0,-1.0,
                (11+12+19+20+21)/5.0,(12+13+20+21+22)/5.0,
                (13+14+21+22+23)/5.0,(14+15+22+23+31)/5.0,
                (16+17+24+25+26)/5.0,-1.0,
                -1.0,-1.0,
                (20+21+28+29+30)/5.0,(21+22+29+30+31)/5.0,
                (22+23+30+31+39)/5.0,(24+25+32+33+34)/5.0,
                (25+26+33+34+35)/5.0,-1.0,
                -1.0,-1.0,
                (29+30+37+38+39)/5.0,(30+31+38+39+47)/5.0,
                (32+33+40+41+42)/5.0,(33+34+41+42+43)/5.0,
                (34+35+42+43+44)/5.0,-1.0,
                -1.0,-1.0,
                (38+39+46+47+55)/5.0,(40+41+48+49+50)/5.0,
                (41+42+49+50+51)/5.0,(42+43+50+51+52)/5.0,
                (43+44+51+52+53)/5.0,-1.0,
                -1.0,-1.0,
                (48+49+56+57+58)/5.0,(49+50+57+58+59)/5.0,
                (50+51+58+59+60)/5.0,(51+52+59+60+61)/5.0,
                (52+53+60+61+62)/5.0,-1.0,
                -1.0}};

            double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
            for (size_t collection_idx = 0; collection_idx < collection_size;
                ++collection_idx)
              tgt_data[collection_idx] =
                xmalloc(grid_data.num_cells * sizeof(**tgt_data));

            yac_interpolation_execute(interpolation, NULL, tgt_data);

            for (size_t collection_idx = 0; collection_idx < collection_size;
                ++collection_idx)
              for (size_t j = 0, offset = collection_idx * 49;
                  j < grid_data.num_cells; ++j) {
                if (ref_tgt_field[i][j] == -1.0) {
                  if (tgt_data[collection_idx][j] != -1.0)
                    PUT_ERR("wrong interpolation result (fallback)");
                } else {
                  if ((tgt_data[collection_idx][j] != 0.0) &&
                      (fabs((ref_tgt_field[i][j] + (double)offset) -
                            tgt_data[collection_idx][j]) > 1e-9))
                    PUT_ERR("wrong interpolation result");
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
    }

    for (size_t i = 0; i < num_stacks; ++i)
      yac_interp_weights_delete(weights[i]);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }

  {
    double coordinates_x[] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double cell_coordinates_x[] = {0.75,1.75,2.75,3.75,4.75,5.75,6.75};
    double cell_coordinates_y[] = {0.75,1.75,2.75,3.75,4.75,5.75,6.75};
    yac_coordinate_pointer cell_coordinates = xmalloc(49 * sizeof(*cell_coordinates));
    size_t const num_cells[2] = {7,7};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {7,7};
    int with_halo = 0;
    int global_cell_mask[7][7] = {
      {1,1,1,1,1,1,1},
      {1,1,1,1,1,1,1},
      {1,1,0,0,0,1,1},
      {1,1,0,0,0,1,1},
      {1,1,0,0,0,1,1},
      {1,1,1,1,1,1,1},
      {1,1,1,1,1,1,1}};
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

    int * tgt_cell_mask =
      xmalloc(grid_data.num_cells * sizeof(*tgt_cell_mask));
    for (size_t i = 0; i < grid_data.num_cells; ++i)
      tgt_cell_mask[i] =
        ((int*)(&(global_cell_mask[0][0])))[grid_data.cell_ids[i]];

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    yac_basic_grid_add_mask_nocpy(
      tgt_grid, YAC_LOC_CELL, tgt_cell_mask, NULL);
    yac_basic_grid_add_coordinates_nocpy(
      tgt_grid, YAC_LOC_CELL, cell_coordinates);
    struct yac_basic_grid * src_grid = yac_basic_grid_empty_new(src_grid_name);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[2] = {
      yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_AVG, .n = 4,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}), NULL};

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    yac_interp_method_delete(method_stack);

    enum yac_interp_weights_reorder_type reorder_type[2] =
      {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};

    for (size_t i = 0; i < 2; ++i) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      {
        double * tgt_field = NULL;
        double ref_tgt_field[49] =
          {(0+1+8+49)/4.0,(0+1+8+6)/4.0,(0+1+8+6)/4.0,(1+6+7+15)/4.0,
            (1+6+7+15)/4.0,(1+6+7+15)/4.0,(1+6+7+15)/4.0,
          (0+1+8+49)/4.0,(0+1+8+49)/4.0,(0+1+8+6)/4.0,(1+6+7+15)/4.0,
            (1+6+7+15)/4.0,(1+6+7+15)/4.0,(1+6+7+15)/4.0,
          (8+1+0+49)/4.0,(8+1+0+49)/4.0,-1,-1,
            -1,(1+6+7+15)/4.0,(1+6+7+15)/4.0,
          (48+49+8+57)/4.0,(48+49+8+57)/4.0,-1,-1,
            -1,(6+7+15+49)/4.0,(6+7+15+49)/4.0,
          (48+49+56+57)/4.0,(48+49+56+57)/4.0,-1,-1,
            -1,(6+7+15+49)/4.0,(6+7+15+49)/4.0,
          (48+49+56+57)/4.0,(48+49+56+57)/4.0,(48+49+56+57)/4.0,(48+49+56+57)/4.0,
            (48+49+56+57)/4.0,(48+49+15+57)/4.0,(6+7+15+49)/4.0,
          (48+49+56+57)/4.0,(48+49+56+57)/4.0,(48+49+56+57)/4.0,(48+49+56+57)/4.0,
            (48+49+56+57)/4.0,(48+49+56+57)/4.0,(49+57+56+15)/4.0};

        tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));

        for (size_t j = 0; j < grid_data.num_cells; ++j) tgt_field[j] = -1.0;

        yac_interpolation_execute(interpolation, NULL, &tgt_field);

        for (size_t j = 0; j < grid_data.num_cells; ++j)
          if (fabs(ref_tgt_field[j] - tgt_field[j]) > 1e-9)
            PUT_ERR("wrong interpolation result");

        free(tgt_field);
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }

  {
    double coordinates_x[] = {1.0,1.5,2.0};
    double coordinates_y[] = {1.0,1.5,2.0};
    size_t const num_cells[2] = {2,2};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {2,2};
    int with_halo = 0;
    int global_corner_mask[3*3] = {1,0,1, 0,1,0, 1,0,1};
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start, local_count, with_halo);

    int * tgt_corner_mask =
      xmalloc(grid_data.num_vertices * sizeof(*tgt_corner_mask));
    for (size_t i = 0; i < grid_data.num_vertices; ++i)
      tgt_corner_mask[i] = global_corner_mask[grid_data.vertex_ids[i]];

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    yac_basic_grid_add_mask_nocpy(
      tgt_grid, YAC_LOC_CORNER, tgt_corner_mask, NULL);
    struct yac_basic_grid * src_grid = yac_basic_grid_empty_new(src_grid_name);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    enum {NUM_STACKS = 2};
    struct interp_method * method_stack[NUM_STACKS][3] = {
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_DIST, .n = 4,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
      yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_GAUSS, .n = 4,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
          .data.gauss_scale = YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT}),
      yac_interp_method_fixed_new(-1), NULL}};


    struct yac_interp_weights * weights[NUM_STACKS];
    for (size_t i = 0; i < NUM_STACKS; ++i) {
      weights[i] = yac_interp_method_do_search(method_stack[i], interp_grid);
      yac_interp_method_delete(method_stack[i]);
    }

    for (size_t i = 0; i < NUM_STACKS; ++i) {
      for (size_t collection_size = 1; collection_size < 4;
          collection_size += 2) {

        struct yac_interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights[i], YAC_MAPPING_ON_SRC, collection_size,
            YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        {
          double ref_tgt_field[][3*3] =
            {{5.0,-2.0,6.0,
              -2.0, compute_dist_result_2(4, (size_t[]){5,6,9,10}, 4),-2.0,
              9.0,-2.0,10.0},
             {5.0,-2.0,6.0,
              -2.0, compute_gauss_result_2(4, (size_t[]){5,6,9,10}, 4),-2.0,
              9.0,-2.0,10.0}};

          double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            tgt_data[collection_idx] =
              xmalloc(grid_data.num_vertices * sizeof(**tgt_data));
            for (size_t i = 0; i < grid_data.num_vertices; ++i)
              tgt_data[collection_idx][i] = -2.0;
          }

          yac_interpolation_execute(interpolation, NULL, tgt_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx)
            for (size_t j = 0, offset = collection_idx * 16;
                j < grid_data.num_vertices; ++j) {
              if (ref_tgt_field[i][j] != -2.0) {
                if ((fabs(
                       (ref_tgt_field[i][j] + (double)offset) -
                       tgt_data[collection_idx][j]) > 1e-9))
                  PUT_ERR("wrong interpolation result");
              } else {
                if (tgt_data[collection_idx][j] != -2.0)
                  PUT_ERR("wrong interpolation result");
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

    for (size_t i = 0; i < NUM_STACKS; ++i)
      yac_interp_weights_delete(weights[i]);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }

  {
    double coordinates_x[] = {1.0+90.0,1.5+90.0,2.0+90.0};
    double coordinates_y[] = {1.0,1.5,2.0};
    size_t const num_cells[2] = {2,2};
    size_t local_start[2] = {0,0};
    size_t local_count[2] = {2,2};
    int with_halo = 0;
    int global_corner_mask[3*3] = {1,0,1, 0,1,0, 1,0,1};
    for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start, local_count, with_halo);

    int * tgt_corner_mask =
      xmalloc(grid_data.num_vertices * sizeof(*tgt_corner_mask));
    for (size_t i = 0; i < grid_data.num_vertices; ++i)
      tgt_corner_mask[i] = global_corner_mask[grid_data.vertex_ids[i]];

    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_new(tgt_grid_name, grid_data);
    yac_basic_grid_add_mask_nocpy(
      tgt_grid, YAC_LOC_CORNER, tgt_corner_mask, NULL);
    struct yac_basic_grid * src_grid = yac_basic_grid_empty_new(src_grid_name);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    enum {NUM_STACKS = 2};
    struct interp_method * method_stack[NUM_STACKS][3] = {
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_DIST, .n = 4,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
      yac_interp_method_fixed_new(-1), NULL},
      {yac_interp_method_nnn_new(
        (struct yac_nnn_config){
          .type = YAC_INTERP_NNN_GAUSS, .n = 4,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
          .data.gauss_scale = YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT}),
      yac_interp_method_fixed_new(-1), NULL}};


    struct yac_interp_weights * weights[NUM_STACKS];
    for (size_t i = 0; i < NUM_STACKS; ++i) {
      weights[i] = yac_interp_method_do_search(method_stack[i], interp_grid);
      yac_interp_method_delete(method_stack[i]);
    }

    for (size_t i = 0; i < NUM_STACKS; ++i) {
      for (size_t collection_size = 1; collection_size < 4;
          collection_size += 2) {

        struct yac_interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights[i], YAC_MAPPING_ON_SRC, collection_size,
            YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        {
          double ref_tgt_field[][3*3] =
            {{compute_dist_result_3(4, (size_t[]){3,7,11,15}, 0),-2.0,
                compute_dist_result_3(4, (size_t[]){3,7,11,15}, 2),
              -2.0, compute_dist_result_3(4, (size_t[]){3,7,11,15}, 4),-2.0,
              compute_dist_result_3(4, (size_t[]){3,7,11,15}, 6),-2.0,
                compute_dist_result_3(4, (size_t[]){3,7,11,15}, 8)},
             {compute_dist_result_3(4, (size_t[]){3,7,11,15}, 0),-2.0,
                compute_dist_result_3(4, (size_t[]){3,7,11,15}, 2),
              -2.0, compute_dist_result_3(4, (size_t[]){3,7,11,15}, 4),-2.0,
              compute_dist_result_3(4, (size_t[]){3,7,11,15}, 6),-2.0,
                compute_dist_result_3(4, (size_t[]){3,7,11,15}, 8)}};

          double ** tgt_data = xmalloc(collection_size * sizeof(*tgt_data));
          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            tgt_data[collection_idx] =
              xmalloc(grid_data.num_vertices * sizeof(**tgt_data));
            for (size_t i = 0; i < grid_data.num_vertices; ++i)
              tgt_data[collection_idx][i] = -2.0;
          }

          yac_interpolation_execute(interpolation, NULL, tgt_data);

          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx)
            for (size_t j = 0, offset = collection_idx * 16;
                j < grid_data.num_vertices; ++j) {
              if (ref_tgt_field[i][j] != -2.0) {
                if ((fabs(
                       (ref_tgt_field[i][j] + (double)offset) -
                       tgt_data[collection_idx][j]) > 1e-9) ||
                    (tgt_data[collection_idx][j] !=
                       tgt_data[collection_idx][j]))
                  PUT_ERR("wrong interpolation result");
              } else {
                if (tgt_data[collection_idx][j] != -2.0)
                  PUT_ERR("wrong interpolation result");
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

    for (size_t i = 0; i < NUM_STACKS; ++i)
      yac_interp_weights_delete(weights[i]);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(src_grid);
    yac_basic_grid_delete(tgt_grid);
  }
}
