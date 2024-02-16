/**
 * @file test_interp_method_creep_parallel.c
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
#include "interp_method_avg.h"
#include "interp_method_conserv.h"
#include "interp_method_creep.h"
#include "interp_method_fixed.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

static char const * grid_names[2] = {"src_grid", "tgt_grid"};

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size != 3) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  MPI_Comm split_comm;
  yac_mpi_call(
    MPI_Comm_split(
      MPI_COMM_WORLD, comm_rank < 2, 0, &split_comm), MPI_COMM_WORLD);

  int split_comm_rank, split_comm_size;
  yac_mpi_call(MPI_Comm_rank(split_comm, &split_comm_rank), split_comm);
  yac_mpi_call(MPI_Comm_size(split_comm, &split_comm_size), split_comm);

  {// stack that starts with linear point interpolation, followed by creep
   // fill and fixed as backup with two source processes and a single target
   // process

    // the global source grid is a 3x2 grid:
    //       08--14--09--15--10--16--11
    //       |       |       |       |
    //       08  03  10  04  12  05  13
    //       |       |       |       |
    //       04--07--05--09--06--11--07
    //       |       |       |       |
    //       01  00  03  01  05  02  06
    //       |       |       |       |
    //       00--00--01--02--02--04--03
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[] = {0.0,1.0,2.0,3.0};
    double coordinates_y[] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0};
    double tgt_cell_coordinates_x[] = {0.5,1.5,2.5};
    double tgt_cell_coordinates_y[] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5};
    double tgt_cell_coordinates[24][3];
    int tgt_cell_mask[24] = {1, 1, 1,
                             1, 1, 1,
                             1, 1, 1,
                             1, 0, 1,
                             1, 0, 0,
                             1, 1, 1,
                             0, 0, 0,
                             1, 1, 1};
    size_t const num_cells[2][2] = {{3,2}, {3,8}};
    size_t local_start[2][2][2] = {{{0,0},{1,0}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{2,2},{2,2}}, {{3,8}}};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[is_tgt][0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[is_tgt][1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells[is_tgt],
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    if (is_tgt) {
      for (size_t i = 0, k = 0; i < num_cells[is_tgt][1]; ++i)
        for (size_t j = 0; j < num_cells[is_tgt][0]; ++j, ++k)
          LLtoXYZ_deg(
            tgt_cell_coordinates_x[j], tgt_cell_coordinates_y[i],
            tgt_cell_coordinates[k]);
      yac_basic_grid_add_coordinates(
        grids[0], CELL, tgt_cell_coordinates, grid_data.num_cells);
      yac_basic_grid_add_mask(
        grids[0], CELL, tgt_cell_mask, grid_data.num_cells, NULL);
    }

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[is_tgt], grids[is_tgt^1], MPI_COMM_WORLD);
    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    for (int creep_distance = -1; creep_distance <= 7; ++creep_distance) {

      struct interp_method * method_stack[] =
        {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
         yac_interp_method_creep_new(creep_distance),
         yac_interp_method_fixed_new(-1.0), NULL};

      struct interp_weights * weights =
        yac_interp_method_do_search(method_stack, interp_grid);

      enum interp_weights_reorder_type reorder_type[2] =
        {MAPPING_ON_SRC, MAPPING_ON_TGT};

      for (size_t i = 0; i < 2; ++i) {

        struct interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        {
          double * src_field = NULL;
          double ** src_fields = &src_field;
          double * tgt_field = NULL;
          double * ref_tgt_field = NULL;
          double ref_global_tgt_field[9][24] =
          {{ 2.5,  3.5,  4.5, // creep_distance = -1
             6.5,  7.5,  8.5,
             6.5,  7.5,  8.5,
             6.5, -2.0,  8.5,
             6.5, -2.0, -2.0,
             6.5,  6.5,  6.5,
            -2.0, -2.0, -2.0,
            -1.0, -1.0, -1.0},
           { 2.5,  3.5,  4.5, // creep_distance = 0
             6.5,  7.5,  8.5,
            -1.0, -1.0, -1.0,
            -1.0, -2.0, -1.0,
            -1.0, -2.0, -2.0,
            -1.0, -1.0, -1.0,
            -2.0, -2.0, -2.0,
            -1.0, -1.0, -1.0},
           { 2.5,  3.5,  4.5, // creep_distance = 1
             6.5,  7.5,  8.5,
             6.5,  7.5,  8.5,
            -1.0, -2.0, -1.0,
            -1.0, -2.0, -2.0,
            -1.0, -1.0, -1.0,
            -2.0, -2.0, -2.0,
            -1.0, -1.0, -1.0},
           { 2.5,  3.5,  4.5, // creep_distance = 2
             6.5,  7.5,  8.5,
             6.5,  7.5,  8.5,
             6.5, -2.0,  8.5,
            -1.0, -2.0, -2.0,
            -1.0, -1.0, -1.0,
            -2.0, -2.0, -2.0,
            -1.0, -1.0, -1.0},
           { 2.5,  3.5,  4.5, // creep_distance = 3
             6.5,  7.5,  8.5,
             6.5,  7.5,  8.5,
             6.5, -2.0,  8.5,
             6.5, -2.0, -2.0,
            -1.0, -1.0, -1.0,
            -2.0, -2.0, -2.0,
            -1.0, -1.0, -1.0},
           { 2.5,  3.5,  4.5, // creep_distance = 4
             6.5,  7.5,  8.5,
             6.5,  7.5,  8.5,
             6.5, -2.0,  8.5,
             6.5, -2.0, -2.0,
             6.5, -1.0, -1.0,
            -2.0, -2.0, -2.0,
            -1.0, -1.0, -1.0},
           { 2.5,  3.5,  4.5, // creep_distance = 5
             6.5,  7.5,  8.5,
             6.5,  7.5,  8.5,
             6.5, -2.0,  8.5,
             6.5, -2.0, -2.0,
             6.5,  6.5, -1.0,
            -2.0, -2.0, -2.0,
            -1.0, -1.0, -1.0},
           { 2.5,  3.5,  4.5, // creep_distance = 6
             6.5,  7.5,  8.5,
             6.5,  7.5,  8.5,
             6.5, -2.0,  8.5,
             6.5, -2.0, -2.0,
             6.5,  6.5,  6.5,
            -2.0, -2.0, -2.0,
            -1.0, -1.0, -1.0},
           { 2.5,  3.5,  4.5, // creep_distance = 7
             6.5,  7.5,  8.5,
             6.5,  7.5,  8.5,
             6.5, -2.0,  8.5,
             6.5, -2.0, -2.0,
             6.5,  6.5,  6.5,
            -2.0, -2.0, -2.0,
            -1.0, -1.0, -1.0}};

          if (is_tgt) {
            tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));
            ref_tgt_field = xmalloc(grid_data.num_cells * sizeof(*ref_tgt_field));
            for (size_t i = 0; i < grid_data.num_cells; ++i) {
              tgt_field[i] = -2.0;
              ref_tgt_field[i] =
                ref_global_tgt_field[creep_distance+1][grid_data.cell_ids[i]];
            }
          } else {
            src_field = xmalloc(grid_data.num_vertices * sizeof(*src_field));
            for (size_t i = 0; i < grid_data.num_vertices; ++i)
              src_field[i] = (double)(grid_data.vertex_ids[i]);
          }

          yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

          if (is_tgt)
            for (size_t i = 0; i < grid_data.num_cells; ++i)
              if (fabs(ref_tgt_field[i] - tgt_field[i]) > 1e-9)
                PUT_ERR("wrong interpolation result");

          free(src_field);
          free(tgt_field);
          free(ref_tgt_field);
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

  {// stack that starts with linear point interpolation, followed by creep
   // fill and fixed as backup with two source processes and a single target
   // process

    // the global source grid is a 7x7 grid
    // the global target grid is a 6x6 grid
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][8] = {{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},
                                  {0.5,1.5,2.5,3.5,4.5,5.5,6.5}};
    double coordinates_y[2][8] = {{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},
                                  {0.5,1.5,2.5,3.5,4.5,5.5,6.5}};
    int vertex_mask[2][8*8] =
      {{1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,
        1,1,0,0,0,0,1,1,
        1,1,0,0,0,0,1,1,
        1,1,0,0,0,0,1,1,
        1,1,0,0,0,0,1,1,
        1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,
         1,1,1,1,1,1,1,
         1,1,1,0,1,1,1,
         1,1,0,1,0,1,1,
         1,1,1,0,1,1,1,
         1,1,1,1,1,1,1,
         1,1,1,1,1,1,1}};
    size_t const num_cells[2][2] = {{7,7}, {6,6}};
    size_t local_start[2][2][2] = {{{0,0},{3,0}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{3,7},{4,7}}, {{6,6}}};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[is_tgt][0]; ++i)
      coordinates_x[is_tgt][i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[is_tgt][1]; ++i)
      coordinates_y[is_tgt][i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x[is_tgt], coordinates_y[is_tgt], num_cells[is_tgt],
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    int temp_vertex_mask[8*8];
    for (size_t i = 0; i < grid_data.num_vertices; ++i)
      temp_vertex_mask[i] = vertex_mask[is_tgt][grid_data.vertex_ids[i]];
    yac_basic_grid_add_mask(
      grids[0], CORNER, temp_vertex_mask, grid_data.num_vertices, NULL);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[is_tgt], grids[is_tgt^1], MPI_COMM_WORLD);
    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    for (int creep_distance = -1; creep_distance <= 3; ++creep_distance) {

      struct interp_method * method_stack[] =
        {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
         yac_interp_method_creep_new(creep_distance),
         yac_interp_method_fixed_new(-1.0), NULL};

      struct interp_weights * weights =
        yac_interp_method_do_search(method_stack, interp_grid);

      enum interp_weights_reorder_type reorder_type[2] =
        {MAPPING_ON_SRC, MAPPING_ON_TGT};

      for (size_t i = 0; i < 2; ++i) {

        struct interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        {
          double * src_field = NULL;
          double ** src_fields = &src_field;
          double * tgt_field = NULL;
          double * ref_tgt_field = NULL;
          double ref_global_tgt_field[5][7*7] =
          {{ 4.5,  5.5,  6.5,  7.5,  8.5,  9.5, 10.5, // creep_distance = -1
            12.5,  9.0,  6.5,  7.5,  8.5, 14.0, 18.5,
            20.5, 20.5, 13.5, -2.0, 17.5, 26.5, 26.5,
            28.5, 28.5, -2.0, -1.0, -2.0, 34.5, 34.5,
            36.5, 36.5, 45.5, -2.0, 49.5, 42.5, 42.5,
            44.5, 49.0, 54.5, 55.5, 56.5, 54.0, 50.5,
            52.5, 53.5, 54.5, 55.5, 56.5, 57.5, 58.5},
           { 4.5,  5.5,  6.5,  7.5,  8.5,  9.5, 10.5, // creep_distance = 0
            12.5, -1.0, -1.0, -1.0, -1.0, -1.0, 18.5,
            20.5, -1.0, -1.0, -2.0, -1.0, -1.0, 26.5,
            28.5, -1.0, -2.0, -1.0, -2.0, -1.0, 34.5,
            36.5, -1.0, -1.0, -2.0, -1.0, -1.0, 42.5,
            44.5, -1.0, -1.0, -1.0, -1.0, -1.0, 50.5,
            52.5, 53.5, 54.5, 55.5, 56.5, 57.5, 58.5},
           { 4.5,  5.5,  6.5,  7.5,  8.5,  9.5, 10.5, // creep_distance = 1
            12.5,  9.0,  6.5,  7.5,  8.5, 14.0, 18.5,
            20.5, 20.5, -1.0, -2.0, -1.0, 26.5, 26.5,
            28.5, 28.5, -2.0, -1.0, -2.0, 34.5, 34.5,
            36.5, 36.5, -1.0, -2.0, -1.0, 42.5, 42.5,
            44.5, 49.0, 54.5, 55.5, 56.5, 54.0, 50.5,
            52.5, 53.5, 54.5, 55.5, 56.5, 57.5, 58.5},
           { 4.5,  5.5,  6.5,  7.5,  8.5,  9.5, 10.5, // creep_distance = 2
            12.5,  9.0,  6.5,  7.5,  8.5, 14.0, 18.5,
            20.5, 20.5, 13.5, -2.0, 17.5, 26.5, 26.5,
            28.5, 28.5, -2.0, -1.0, -2.0, 34.5, 34.5,
            36.5, 36.5, 45.5, -2.0, 49.5, 42.5, 42.5,
            44.5, 49.0, 54.5, 55.5, 56.5, 54.0, 50.5,
            52.5, 53.5, 54.5, 55.5, 56.5, 57.5, 58.5},
           { 4.5,  5.5,  6.5,  7.5,  8.5,  9.5, 10.5, // creep_distance = 3
            12.5,  9.0,  6.5,  7.5,  8.5, 14.0, 18.5,
            20.5, 20.5, 13.5, -2.0, 17.5, 26.5, 26.5,
            28.5, 28.5, -2.0, -1.0, -2.0, 34.5, 34.5,
            36.5, 36.5, 45.5, -2.0, 49.5, 42.5, 42.5,
            44.5, 49.0, 54.5, 55.5, 56.5, 54.0, 50.5,
            52.5, 53.5, 54.5, 55.5, 56.5, 57.5, 58.5}};

          if (is_tgt) {
            tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
            ref_tgt_field = xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
            for (size_t i = 0; i < grid_data.num_vertices; ++i) {
              tgt_field[i] = -2.0;
              ref_tgt_field[i] =
                ref_global_tgt_field[creep_distance+1][grid_data.vertex_ids[i]];
            }
          } else {
            src_field = xmalloc(grid_data.num_vertices * sizeof(*src_field));
            for (size_t i = 0; i < grid_data.num_vertices; ++i)
              src_field[i] = (double)(grid_data.vertex_ids[i]);
          }

          yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

          if (is_tgt)
            for (size_t i = 0; i < grid_data.num_vertices; ++i)
              if (fabs(ref_tgt_field[i] - tgt_field[i]) > 1e-9)
                PUT_ERR("wrong interpolation result");

          free(src_field);
          free(tgt_field);
          free(ref_tgt_field);
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

  {// stack that starts with linear point interpolation, followed by creep
   // fill and fixed as backup with two source processes and a single target
   // process

    // the global source grid is a 5x1 grid
    // the global target grid is a 5x3 grid
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[6] = {0.0,1.0,2.0,3.0,4.0,5.0};
    double coordinates_y[4] = {0.0,1.0,2.0,3.0};
    int cell_mask[2][5*3] =
      {{1,0,1,0,1,},
       {1,1,1,1,1,
        1,1,0,1,1,
        1,1,1,1,1}};
    size_t const num_cells[2][2] = {{5,1}, {5,3}};
    size_t local_start[2][2][2] = {{{0,0},{3,0}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{3,1},{2,1}}, {{5,3}}};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[is_tgt][0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[is_tgt][1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells[is_tgt],
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    int temp_cell_mask[5*3];
    for (size_t i = 0; i < grid_data.num_cells; ++i)
      temp_cell_mask[i] = cell_mask[is_tgt][grid_data.cell_ids[i]];
    yac_basic_grid_add_mask(
      grids[0], CELL, temp_cell_mask, grid_data.num_cells, NULL);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[is_tgt], grids[is_tgt^1], MPI_COMM_WORLD);
    struct interp_field src_fields[] =
      {{.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = SIZE_MAX, .masks_idx = 0};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    for (int creep_distance = -1; creep_distance <= 5; ++creep_distance) {

      struct interp_method * method_stack[] =
        {yac_interp_method_conserv_new(1, 0, 1, CONSERV_DESTAREA),
         yac_interp_method_creep_new(creep_distance),
         yac_interp_method_fixed_new(-1.0), NULL};

      struct interp_weights * weights =
        yac_interp_method_do_search(method_stack, interp_grid);

      enum interp_weights_reorder_type reorder_type[2] =
        {MAPPING_ON_SRC, MAPPING_ON_TGT};

      for (size_t i = 0; i < 2; ++i) {

        struct interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, reorder_type[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        {
          double * src_field = NULL;
          double ** src_fields = &src_field;
          double * tgt_field = NULL;
          double * ref_tgt_field = NULL;
          double ref_global_tgt_field[7][7*7] =
          {{ 0.0,  1.0 ,  2.0 ,  3.0 ,  4.0, // creep_distance = -1
             0.0,  0.5 , -2.0 ,  3.5 ,  4.0,
             0.0,  0.25,  2.0 ,  3.75,  4.0},
           { 0.0, -1.0 ,  2.0 , -1.0 ,  4.0, // creep_distance = 0
            -1.0, -1.0 , -2.0 , -1.0 , -1.0,
            -1.0, -1.0 , -1.0 , -1.0 , -1.0},
           { 0.0,  1.0 ,  2.0 ,  3.0 ,  4.0, // creep_distance = 1
             0.0, -1.0 , -2.0 , -1.0 ,  4.0,
            -1.0, -1.0 , -1.0 , -1.0 , -1.0},
           { 0.0,  1.0 ,  2.0 ,  3.0 ,  4.0, // creep_distance = 2
             0.0,  0.5 , -2.0 ,  3.5 ,  4.0,
             0.0, -1.0 , -1.0 , -1.0 ,  4.0},
           { 0.0,  1.0 ,  2.0 ,  3.0 ,  4.0, // creep_distance = 3
             0.0,  0.5 , -2.0 ,  3.5 ,  4.0,
             0.0,  0.25, -1.0 ,  3.75,  4.0},
           { 0.0,  1.0 ,  2.0 ,  3.0 ,  4.0, // creep_distance = 4
             0.0,  0.5 , -2.0 ,  3.5 ,  4.0,
             0.0,  0.25,  2.0 ,  3.75,  4.0},
           { 0.0,  1.0 ,  2.0 ,  3.0 ,  4.0, // creep_distance = 5
             0.0,  0.5 , -2.0 ,  3.5 ,  4.0,
             0.0,  0.25,  2.0 ,  3.75,  4.0}};

          if (is_tgt) {
            tgt_field = xmalloc(grid_data.num_cells* sizeof(*tgt_field));
            ref_tgt_field = xmalloc(grid_data.num_cells* sizeof(*ref_tgt_field));
            for (size_t i = 0; i < grid_data.num_cells; ++i) {
              tgt_field[i] = -2.0;
              ref_tgt_field[i] =
                ref_global_tgt_field[creep_distance+1][grid_data.cell_ids[i]];
            }
          } else {
            src_field = xmalloc(grid_data.num_cells* sizeof(*src_field));
            for (size_t i = 0; i < grid_data.num_cells; ++i)
              src_field[i] = (double)(grid_data.cell_ids[i]);
          }

          yac_interpolation_execute(interpolation, &src_fields, &tgt_field);
 
          if (is_tgt)
            for (size_t i = 0; i < grid_data.num_cells; ++i)
              if (fabs(ref_tgt_field[i] - tgt_field[i]) > 1e-9)
                PUT_ERR("wrong interpolation result");

          free(src_field);
          free(tgt_field);
          free(ref_tgt_field);
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

  yac_mpi_call(MPI_Comm_free(&split_comm), MPI_COMM_WORLD);

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
