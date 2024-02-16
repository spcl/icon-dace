/**
 * @file test_interp_method_spmap_parallel.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "interp_method.h"
#include "interp_method_spmap.h"
#include "interp_method_fixed.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

enum yac_interp_spmap_weight_type weight_types[] = {
  SPMAP_AVG,
  SPMAP_DIST
};
size_t num_weight_types = sizeof(weight_types) / sizeof(weight_types[0]);

enum interp_weights_reorder_type reorder_types[] = {
  MAPPING_ON_SRC,
  MAPPING_ON_TGT
};
size_t num_reorder_types = sizeof(reorder_types) / sizeof(reorder_types[0]);

static char const * grid_names[2] = {"grid1", "grid2"};

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size != 6) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  MPI_Comm split_comm;
  yac_mpi_call(
    MPI_Comm_split(
      MPI_COMM_WORLD, comm_rank < 1, 0, &split_comm), MPI_COMM_WORLD);

  int split_comm_rank, split_comm_size;
  yac_mpi_call(MPI_Comm_rank(split_comm, &split_comm_rank), split_comm);
  yac_mpi_call(MPI_Comm_size(split_comm, &split_comm_size), split_comm);

  { // parallel interpolation process
    // corner and cell ids for a 7 x 7 grid
    // 56-----57-----58-----59-----60-----61-----62-----63
    //  |      |      |      |      |      |      |      |
    //  |  42  |  43  |  44  |  45  |  46  |  47  |  48  |
    //  |      |      |      |      |      |      |      |
    // 48-----49-----50-----51-----52-----53-----54-----55
    //  |      |      |      |      |      |      |      |
    //  |  35  |  36  |  37  |  38  |  39  |  40  |  41  |
    //  |      |      |      |      |      |      |      |
    // 40-----41-----42-----43-----44-----45-----46-----47
    //  |      |      |      |      |      |      |      |
    //  |  28  |  29  |  30  |  31  |  32  |  33  |  34  |
    //  |      |      |      |      |      |      |      |
    // 32-----33-----34-----35-----36-----37-----38-----39
    //  |      |      |      |      |      |      |      |
    //  |  21  |  22  |  23  |  24  |  25  |  26  |  27  |
    //  |      |      |      |      |      |      |      |
    // 24-----25-----26-----27-----28-----29-----30-----31
    //  |      |      |      |      |      |      |      |
    //  |  14  |  15  |  16  |  17  |  18  |  19  |  20  |
    //  |      |      |      |      |      |      |      |
    // 16-----17-----18-----19-----20-----21-----22-----23
    //  |      |      |      |      |      |      |      |
    //  |  07  |  08  |  09  |  10  |  11  |  12  |  13  |
    //  |      |      |      |      |      |      |      |
    // 08-----09-----10-----11-----12-----13-----14-----15
    //  |      |      |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |  05  |  06  |
    //  |      |      |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05-----06-----07
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 3---3---3---3---4---4---4---4
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 2---2---2---2---3---3---3---3
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 1---1---1---1---2---2---2---2
    // | 1 | 1 | 1 | 1 | 3 | 3 | 3 |
    // 0---0---0---0---1---1---1---1
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    //
    // mask
    // land = 0
    // coast = 1
    // ocean = 2
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 1 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 0 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 1 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 1 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 1 | 1 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 1 | 1 | 1 | 1 | 1 |
    // +---+---+---+---+---+---+---+
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{7,7}}};
    int global_mask[7*7] = {
      0,0,1,1,1,1,1,
      0,1,1,2,2,2,2,
      1,1,2,2,2,2,2,
      1,2,2,1,1,1,2,
      1,2,2,1,0,1,2,
      1,2,2,1,1,1,2,
      1,2,2,2,2,2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    {
      int valid_mask_value = (is_tgt)?2:1;
      coordinate_pointer point_coordinates =
        xmalloc(grid_data.num_cells * sizeof(*point_coordinates));
      int * mask = xmalloc(grid_data.num_cells * sizeof(*mask));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = point_coordinates[i];
        for (size_t k = 0; k < 3; ++k) middle_point[k] = 0.0;
        size_t * curr_vertices =
          grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
        size_t curr_num_vertices = grid_data.num_vertices_per_cell[i];
        for (size_t j = 0; j < curr_num_vertices; ++j) {
          double * curr_vertex_coord =
            grid_data.vertex_coordinates[curr_vertices[j]];
          for (size_t k = 0; k < 3; ++k)
            middle_point[k] += curr_vertex_coord[k];
        }
        normalise_vector(middle_point);
        mask[i] = global_mask[grid_data.cell_ids[i]] == valid_mask_value;
      }
      yac_basic_grid_add_coordinates_nocpy(grids[0], CELL, point_coordinates);
      yac_basic_grid_add_mask_nocpy(grids[0], CELL, mask, NULL);
    }

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CELL, .coordinates_idx = 0, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_spmap_new(
         YAC_DEFAULT_SPREAD_DISTANCE,
         YAC_DEFAULT_MAX_SEARCH_DISTANCE,
         SPMAP_AVG),
       yac_interp_method_fixed_new(-2.0), NULL};

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    for (size_t i = 0; i < num_reorder_types; ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_types[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      {
        double * src_field = NULL;
        double ** src_fields = &src_field;
        double * tgt_field = NULL;

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i) tgt_field[i] = -1;
        } else {
          src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            src_field[i] = (double)(grid_data.cell_ids[i]);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        double ref_tgt_field[7*7] = {
          -1,-1,-1,-1,-1,-1,-1,
          -1,-1,-1,2+3+9,4,5,6,
          -1,-1,8+15,-2,25,-2,-2,
          -1,14+21,24,-1,-1,-1,26,
          -1,28,31,-1,-1,-1,33,
          -1,35,38,-1,-1,-1,40,
          -1,42,-2,-2,39,-2,-2};

        if (is_tgt)
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            if (ref_tgt_field[grid_data.cell_ids[i]] != tgt_field[i])
              PUT_ERR("wrong interpolation result");

        free(tgt_field);
        free(src_field);
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
    yac_interp_method_delete(method_stack);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  { // parallel interpolation process
    // corner and cell ids for a 7 x 7 grid
    // 56-----57-----58-----59-----60-----61-----62-----63
    //  |      |      |      |      |      |      |      |
    //  |  42  |  43  |  44  |  45  |  46  |  47  |  48  |
    //  |      |      |      |      |      |      |      |
    // 48-----49-----50-----51-----52-----53-----54-----55
    //  |      |      |      |      |      |      |      |
    //  |  35  |  36  |  37  |  38  |  39  |  40  |  41  |
    //  |      |      |      |      |      |      |      |
    // 40-----41-----42-----43-----44-----45-----46-----47
    //  |      |      |      |      |      |      |      |
    //  |  28  |  29  |  30  |  31  |  32  |  33  |  34  |
    //  |      |      |      |      |      |      |      |
    // 32-----33-----34-----35-----36-----37-----38-----39
    //  |      |      |      |      |      |      |      |
    //  |  21  |  22  |  23  |  24  |  25  |  26  |  27  |
    //  |      |      |      |      |      |      |      |
    // 24-----25-----26-----27-----28-----29-----30-----31
    //  |      |      |      |      |      |      |      |
    //  |  14  |  15  |  16  |  17  |  18  |  19  |  20  |
    //  |      |      |      |      |      |      |      |
    // 16-----17-----18-----19-----20-----21-----22-----23
    //  |      |      |      |      |      |      |      |
    //  |  07  |  08  |  09  |  10  |  11  |  12  |  13  |
    //  |      |      |      |      |      |      |      |
    // 08-----09-----10-----11-----12-----13-----14-----15
    //  |      |      |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |  05  |  06  |
    //  |      |      |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05-----06-----07
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 3---3---3---3---4---4---4---4
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 2---2---2---2---3---3---3---3
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 1---1---1---1---2---2---2---2
    // | 1 | 1 | 1 | 1 | 3 | 3 | 3 |
    // 0---0---0---0---1---1---1---1
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    //
    // mask
    // land = 0
    // coast = 1
    // ocean = 2
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 1 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 0 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 1 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 1 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 1 | 1 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 1 | 1 | 1 | 1 | 1 |
    // +---+---+---+---+---+---+---+
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{7,7}}};
    int global_mask[7*7] = {
      0,0,1,1,1,1,1,
      0,1,1,2,2,2,2,
      1,1,2,2,2,2,2,
      1,2,2,1,1,1,2,
      1,2,2,1,0,1,2,
      1,2,2,1,1,1,2,
      1,2,2,2,2,2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    {
      int valid_mask_value = (is_tgt)?2:1;
      coordinate_pointer point_coordinates =
        xmalloc(grid_data.num_cells * sizeof(*point_coordinates));
      int * mask = xmalloc(grid_data.num_cells * sizeof(*mask));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = point_coordinates[i];
        for (size_t k = 0; k < 3; ++k) middle_point[k] = 0.0;
        size_t * curr_vertices =
          grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
        size_t curr_num_vertices = grid_data.num_vertices_per_cell[i];
        for (size_t j = 0; j < curr_num_vertices; ++j) {
          double * curr_vertex_coord =
            grid_data.vertex_coordinates[curr_vertices[j]];
          for (size_t k = 0; k < 3; ++k)
            middle_point[k] += curr_vertex_coord[k];
        }
        normalise_vector(middle_point);
        mask[i] = global_mask[grid_data.cell_ids[i]] == valid_mask_value;
      }
      yac_basic_grid_add_coordinates_nocpy(grids[0], CELL, point_coordinates);
      yac_basic_grid_add_mask_nocpy(grids[0], CELL, mask, NULL);
    }

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CELL, .coordinates_idx = 0, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_spmap_new(
         YAC_RAD * 1.1, YAC_DEFAULT_MAX_SEARCH_DISTANCE, SPMAP_AVG),
       yac_interp_method_fixed_new(-2.0), NULL};

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    for (size_t i = 0; i < num_reorder_types; ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_types[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      {
        double * src_field = NULL;
        double ** src_fields = &src_field;
        double * tgt_field = NULL;

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i) tgt_field[i] = -1;
        } else {
          src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            src_field[i] = (double)(grid_data.cell_ids[i]);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        double ref_tgt_field[7*7] = {
          -1,-1,-1,-1,-1,-1,-1,
          -1,-1,-1, 0, 0, 0, 0,
          -1,-1, 0, 0, 0, 0, 0,
          -1, 0, 0,-1,-1,-1, 0,
          -1, 0, 0,-1,-1,-1, 0,
          -1, 0, 0,-1,-1,-1, 0,
          -1, 0, 0, 0, 0, 0, 0};
        size_t coast_point[] = {2,3,4,5,6,
                                8,9,
                                14,15,
                                21,24,25,26,
                                28,31,33,
                                35,38,39,40,
                                42};
        size_t num_coast_points = sizeof(coast_point)/sizeof(coast_point[0]);
        size_t num_tgt_per_coast[] = {3,3,4,4,3,
                                      3,3,
                                      3,3,
                                      3,4,4,3,
                                      4,4,3,
                                      4,4,3,3,
                                      3};
        size_t tgts[] = {10,11,17, 10,11,17, 10,11,12,18, 11,12,13,19, 12,13,20,
                         16,17,23, 10,11,17,
                         22,23,29, 16,17,23,
                         22,23,29, 16,22,23,30, 11,17,18,19, 20,27,34,
                         22,29,30,36, 23,29,30,37, 27,34,41,
                         29,36,37,43, 30,36,37,44, 45,46,47, 34,41,48,
                         36,43,44};

        for (size_t i = 0, k = 0; i < num_coast_points; ++i) {
          size_t curr_num_tgt = num_tgt_per_coast[i];
          double curr_data = (double)(coast_point[i]) / (double)curr_num_tgt;
          for (size_t j = 0; j < curr_num_tgt; ++j, ++k)
            ref_tgt_field[tgts[k]] += curr_data;
        }

        if (is_tgt)
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            if (fabs(ref_tgt_field[grid_data.cell_ids[i]] -
                     tgt_field[i]) > 1e-6)
              PUT_ERR("wrong interpolation result");

        free(tgt_field);
        free(src_field);
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
    yac_interp_method_delete(method_stack);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  { // parallel interpolation process
    // corner and cell ids for a 7 x 7 grid
    // 56-----57-----58-----59-----60-----61-----62-----63
    //  |      |      |      |      |      |      |      |
    //  |  42  |  43  |  44  |  45  |  46  |  47  |  48  |
    //  |      |      |      |      |      |      |      |
    // 48-----49-----50-----51-----52-----53-----54-----55
    //  |      |      |      |      |      |      |      |
    //  |  35  |  36  |  37  |  38  |  39  |  40  |  41  |
    //  |      |      |      |      |      |      |      |
    // 40-----41-----42-----43-----44-----45-----46-----47
    //  |      |      |      |      |      |      |      |
    //  |  28  |  29  |  30  |  31  |  32  |  33  |  34  |
    //  |      |      |      |      |      |      |      |
    // 32-----33-----34-----35-----36-----37-----38-----39
    //  |      |      |      |      |      |      |      |
    //  |  21  |  22  |  23  |  24  |  25  |  26  |  27  |
    //  |      |      |      |      |      |      |      |
    // 24-----25-----26-----27-----28-----29-----30-----31
    //  |      |      |      |      |      |      |      |
    //  |  14  |  15  |  16  |  17  |  18  |  19  |  20  |
    //  |      |      |      |      |      |      |      |
    // 16-----17-----18-----19-----20-----21-----22-----23
    //  |      |      |      |      |      |      |      |
    //  |  07  |  08  |  09  |  10  |  11  |  12  |  13  |
    //  |      |      |      |      |      |      |      |
    // 08-----09-----10-----11-----12-----13-----14-----15
    //  |      |      |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |  05  |  06  |
    //  |      |      |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05-----06-----07
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 3---3---3---3---4---4---4---4
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 2---2---2---2---3---3---3---3
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 1---1---1---1---2---2---2---2
    // | 1 | 1 | 1 | 1 | 3 | 3 | 3 |
    // 0---0---0---0---1---1---1---1
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    //
    // mask
    // land = 0
    // coast = 1
    // ocean = 2
    // +---+---+---+---+---+---+---+
    // | 1 | 0 | 0 | 0 | 0 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 2 | 1 | 0 | 0 | 1 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 2 | 2 | 1 | 1 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 2 | 2 | 1 | 1 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 2 | 1 | 0 | 0 | 1 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 0 | 0 | 0 | 0 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
    // +---+---+---+---+---+---+---+
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{7,7}}};
    int global_mask[7*7] = {
      0,0,0,0,0,0,1,
      1,0,0,0,0,1,2,
      2,1,0,0,1,2,2,
      2,2,1,1,2,2,2,
      2,2,1,1,2,2,2,
      2,1,0,0,1,2,2,
      1,0,0,0,0,1,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    coordinate_pointer point_coordinates =
      xmalloc(grid_data.num_cells * sizeof(*point_coordinates));
    {
      int valid_mask_value = (is_tgt)?2:1;
      int * mask = xmalloc(grid_data.num_cells * sizeof(*mask));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = point_coordinates[i];
        for (size_t k = 0; k < 3; ++k) middle_point[k] = 0.0;
        size_t * curr_vertices =
          grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
        size_t curr_num_vertices = grid_data.num_vertices_per_cell[i];
        for (size_t j = 0; j < curr_num_vertices; ++j) {
          double * curr_vertex_coord =
            grid_data.vertex_coordinates[curr_vertices[j]];
          for (size_t k = 0; k < 3; ++k)
            middle_point[k] += curr_vertex_coord[k];
        }
        normalise_vector(middle_point);
        mask[i] = global_mask[grid_data.cell_ids[i]] == valid_mask_value;
      }
      yac_basic_grid_add_coordinates_nocpy(grids[0], CELL, point_coordinates);
      yac_basic_grid_add_mask_nocpy(grids[0], CELL, mask, NULL);
    }

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CELL, .coordinates_idx = 0, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    for (size_t weight_type_idx = 0; weight_type_idx < num_weight_types;
         ++weight_type_idx) {

      struct interp_method * method_stack[] =
        {yac_interp_method_spmap_new(
           YAC_RAD * 3.7, YAC_DEFAULT_MAX_SEARCH_DISTANCE,
           weight_types[weight_type_idx]),
         yac_interp_method_fixed_new(-2.0), NULL};

      struct interp_weights * weights =
        yac_interp_method_do_search(method_stack, interp_grid);

      for (size_t i = 0; i < num_reorder_types; ++i) {

        struct interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, reorder_types[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        {
          double * src_field = NULL;
          double ** src_fields = &src_field;
          double * tgt_field = NULL;

          if (is_tgt) {
            tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));
            for (size_t i = 0; i < grid_data.num_cells; ++i) tgt_field[i] = -1;
          } else {
            src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
            for (size_t i = 0; i < grid_data.num_cells; ++i)
              src_field[i] = (double)(grid_data.cell_ids[i]);
          }

          yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

          if (is_tgt) {

            double ref_tgt_field[7*7] = {
               -1,-1,-1,-1,-1,-1,-1,
               -1,-1,-1,-1,-1,-1, 0,
                0,-1,-1,-1,-1, 0, 0,
                0, 0,-1,-1, 0, 0, 0,
                0, 0,-1,-1, 0, 0, 0,
                0,-1,-1,-1,-1, 0, 0,
               -1,-1,-1,-1,-1,-1, 0};
            size_t coast_point[] = {7,15,23,30,36,42,
                                    6,12,18,24,31,39,47};
            size_t num_coast_points = sizeof(coast_point)/sizeof(coast_point[0]);
            size_t num_tgt_per_coast[] = {6,6,6,6,6,6,
                                          9,9,11,12,12,11,9};
            size_t tgts[] = {14,21,22,28,29,35,
                             14,21,22,28,29,35,
                             14,21,22,28,29,35,
                             14,21,22,28,29,35,
                             14,21,22,28,29,35,
                             14,21,22,28,29,35,

                             13,19,20,25,26,27,32,33,34,
                             13,19,20,25,26,27,32,33,34,
                             13,19,20,25,26,27,32,33,34,40,41,
                             13,19,20,25,26,27,32,33,34,40,41,48,
                             13,19,20,25,26,27,32,33,34,40,41,48,
                             19,20,25,26,27,32,33,34,40,41,48,
                             25,26,27,32,33,34,40,41,48};

            switch (weight_types[weight_type_idx]) {
              default:
              case(SPMAP_AVG): {
                for (size_t i = 0, k = 0; i < num_coast_points; ++i) {
                  size_t curr_num_tgt = num_tgt_per_coast[i];
                  double curr_data = (double)(coast_point[i]) / (double)curr_num_tgt;
                  for (size_t j = 0; j < curr_num_tgt; ++j, ++k)
                    ref_tgt_field[tgts[k]] += curr_data;
                }
                break;
              }
              case(SPMAP_DIST): {
                size_t * curr_tgts = tgts;
                for (size_t i = 0; i < num_coast_points; ++i) {
                  size_t curr_num_tgt = num_tgt_per_coast[i];
                  size_t curr_coast_point = coast_point[i];
                  double curr_src_data = (double)(curr_coast_point);
                  double inv_distances[curr_num_tgt];
                  double inv_distances_sum = 0.0;
                  for (size_t j = 0; j < curr_num_tgt; ++j)
                    inv_distances_sum +=
                      ((inv_distances[j] =
                          1.0 / get_vector_angle(
                                  point_coordinates[curr_coast_point],
                                  point_coordinates[curr_tgts[j]])));
                  for (size_t j = 0; j < curr_num_tgt; ++j)
                    ref_tgt_field[curr_tgts[j]] +=
                      curr_src_data * (inv_distances[j] / inv_distances_sum);
                  curr_tgts += curr_num_tgt;
                }
              }
            }

            for (size_t i = 0; i < grid_data.num_cells; ++i)
              if (fabs(ref_tgt_field[grid_data.cell_ids[i]] -
                       tgt_field[i]) > 1e-6)
                PUT_ERR("wrong interpolation result");
          }

          free(tgt_field);
          free(src_field);
        }

        yac_interpolation_delete(interpolation);
      }

      yac_interp_weights_delete(weights);
      yac_interp_method_delete(method_stack);
    }
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  { // parallel interpolation process
    // corner and cell ids for a 7 x 7 grid
    // 56-----57-----58-----59-----60-----61-----62-----63
    //  |      |      |      |      |      |      |      |
    //  |  42  |  43  |  44  |  45  |  46  |  47  |  48  |
    //  |      |      |      |      |      |      |      |
    // 48-----49-----50-----51-----52-----53-----54-----55
    //  |      |      |      |      |      |      |      |
    //  |  35  |  36  |  37  |  38  |  39  |  40  |  41  |
    //  |      |      |      |      |      |      |      |
    // 40-----41-----42-----43-----44-----45-----46-----47
    //  |      |      |      |      |      |      |      |
    //  |  28  |  29  |  30  |  31  |  32  |  33  |  34  |
    //  |      |      |      |      |      |      |      |
    // 32-----33-----34-----35-----36-----37-----38-----39
    //  |      |      |      |      |      |      |      |
    //  |  21  |  22  |  23  |  24  |  25  |  26  |  27  |
    //  |      |      |      |      |      |      |      |
    // 24-----25-----26-----27-----28-----29-----30-----31
    //  |      |      |      |      |      |      |      |
    //  |  14  |  15  |  16  |  17  |  18  |  19  |  20  |
    //  |      |      |      |      |      |      |      |
    // 16-----17-----18-----19-----20-----21-----22-----23
    //  |      |      |      |      |      |      |      |
    //  |  07  |  08  |  09  |  10  |  11  |  12  |  13  |
    //  |      |      |      |      |      |      |      |
    // 08-----09-----10-----11-----12-----13-----14-----15
    //  |      |      |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |  05  |  06  |
    //  |      |      |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05-----06-----07
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 3---3---3---3---4---4---4---4
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 2---2---2---2---3---3---3---3
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 1---1---1---1---2---2---2---2
    // | 1 | 1 | 1 | 1 | 3 | 3 | 3 |
    // 0---0---0---0---1---1---1---1
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    //
    // mask
    // land = 0
    // coast = 1
    // ocean = 2
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 1 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 0 | 1 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 0 | 1 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 1 | 0 | 0 | 1 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
    // +---+---+---+---+---+---+---+
    // | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
    // +---+---+---+---+---+---+---+
    // | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
    // +---+---+---+---+---+---+---+
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{7,7}}};
    int global_mask[7*7] = {
      1,0,0,0,0,0,0,
      0,1,0,0,0,0,0,
      0,0,1,0,0,0,0,
      0,0,0,1,0,0,1,
      0,0,0,0,1,1,2,
      0,0,0,0,1,2,2,
      0,0,0,1,2,2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    coordinate_pointer point_coordinates =
      xmalloc(grid_data.num_cells * sizeof(*point_coordinates));
    {
      int valid_mask_value = (is_tgt)?2:1;
      int * mask = xmalloc(grid_data.num_cells * sizeof(*mask));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = point_coordinates[i];
        for (size_t k = 0; k < 3; ++k) middle_point[k] = 0.0;
        size_t * curr_vertices =
          grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
        size_t curr_num_vertices = grid_data.num_vertices_per_cell[i];
        for (size_t j = 0; j < curr_num_vertices; ++j) {
          double * curr_vertex_coord =
            grid_data.vertex_coordinates[curr_vertices[j]];
          for (size_t k = 0; k < 3; ++k)
            middle_point[k] += curr_vertex_coord[k];
        }
        normalise_vector(middle_point);
        mask[i] = global_mask[grid_data.cell_ids[i]] == valid_mask_value;
      }
      yac_basic_grid_add_coordinates_nocpy(grids[0], CELL, point_coordinates);
      yac_basic_grid_add_mask_nocpy(grids[0], CELL, mask, NULL);
    }

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CELL, .coordinates_idx = 0, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_spmap_new(
         YAC_DEFAULT_SPREAD_DISTANCE, YAC_RAD * 3.0, SPMAP_AVG),
       yac_interp_method_fixed_new(-2.0), NULL};

    struct interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    for (size_t i = 0; i < num_reorder_types; ++i) {

      struct interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_types[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      // check generated interpolation
      {
        double * src_field = NULL;
        double ** src_fields = &src_field;
        double * tgt_field = NULL;

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i) tgt_field[i] = -1;
        } else {
          src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            src_field[i] = (double)(grid_data.cell_ids[i]);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        if (is_tgt) {

          double ref_tgt_field[7*7] = {
             -1,-1,-1,-1,-1,-1,-1,
             -1,-1,-1,-1,-1,-1,-1,
             -1,-1,-1,-1,-1,-1,-1,
             -1,-1,-1,-1,-1,-1,-1,
             -1,-1,-1,-1,-1,-1, 0,
             -1,-1,-1,-1,-1, 0,-2,
             -1,-1,-1,-1, 0,-2,-2};
          size_t coast_point[] = {0,8,16,24,27,32,33,39,45};
          size_t num_coast_points = sizeof(coast_point)/sizeof(coast_point[0]);
          size_t num_tgt_per_coast[] = {0,0,0,1,1,1,1,1,1};
          size_t tgts[] = {40,34,40,34,40,46};

          for (size_t i = 0, k = 0; i < num_coast_points; ++i) {
            size_t curr_num_tgt = num_tgt_per_coast[i];
            if (curr_num_tgt == 0) continue;
            double curr_data = (double)(coast_point[i]) / (double)curr_num_tgt;
            for (size_t j = 0; j < curr_num_tgt; ++j, ++k)
              ref_tgt_field[tgts[k]] += curr_data;
          }

          for (size_t i = 0; i < grid_data.num_cells; ++i)
            if (fabs(ref_tgt_field[grid_data.cell_ids[i]] -
                     tgt_field[i]) > 1e-6)
              PUT_ERR("wrong interpolation result");
        }

        free(tgt_field);
        free(src_field);
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
    yac_interp_method_delete(method_stack);

    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  yac_mpi_call(MPI_Comm_free(&split_comm), MPI_COMM_WORLD);

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
