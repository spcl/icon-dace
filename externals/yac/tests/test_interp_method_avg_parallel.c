/**
 * @file test_interp_method_avg_parallel.c
 *
 * @copyright Copyright  (C)  2019 DKRZ, MPI-M
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
#include "interp_method_fixed.h"
#include "dist_grid_utils.h"
#include "clipping.h"
#include "yac_mpi.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

enum interp_weights_reorder_type reorder_types[] =
  {MAPPING_ON_SRC, MAPPING_ON_TGT};
size_t num_reorder_types = sizeof(reorder_types) / sizeof(reorder_types[0]);
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

  {// Test 1a
   // linear point interpolation (fixed is backup) with two source processes and a
   // single target process the source fields have two struct points per grid
   // (the test checks how the interpolation methods handle NULL target coordinates )

    // the global grid is a 4x2 grid:
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
    double coordinates_x[2][4] = {{0.0,1.0,2.0,3.0}, {0.5,1.5,2.5}};
    double coordinates_y[2][3] = {{0.0,1.0,2.0}, {0.5,1.5,2.5}};
    size_t const num_cells[2][2] = {{3,2}, {2,2}};
    size_t local_start[2][2][2] = {{{0,0},{1,0}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{2,2},{2,2}}, {{2,2}}};
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

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    for (int partial_coverage = 0; partial_coverage <= 1; ++partial_coverage) {

      struct interp_method * method_stack[] =
        {yac_interp_method_avg_new(AVG_ARITHMETIC, partial_coverage),
        yac_interp_method_fixed_new(-1.0), NULL};

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
          double * ref_tgt_field = NULL;
          double ref_global_tgt_field[9] =
            {2.5, 3.5, 4.5, 6.5, 7.5, 8.5, -1, -1, -1};

          if (is_tgt) {
            tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
            ref_tgt_field = xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
            for (size_t i = 0; i < grid_data.num_vertices; ++i)
              ref_tgt_field[i] = ref_global_tgt_field[grid_data.vertex_ids[i]];
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
    } // partial_coverage
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  {// Test 2a
   // linear point interpolation (fixed is backup) with two source processes and a
   // single target process the source fields have two struct points per grid
   // (the test check how the interpolation methods handle NULL target coordinates )

    // the global grid is a 4x2 grid:
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
    // source mask:
    //       0-------1-------1-------1
    //       |       |       |       |
    //       |       |       |       |
    //       |       |       |       |
    //       1-------1-------1-------1
    //       |       |       |       |
    //       |       |       |       |
    //       |       |       |       |
    //       1-------1-------1-------0
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][4] = {{0.0,1.0,2.0,3.0}, {0.5,1.5,2.5}};
    double coordinates_y[2][3] = {{0.0,1.0,2.0}, {0.5,1.5,2.5}};
    size_t const num_cells[2][2] = {{3,2}, {2,2}};
    size_t local_start[2][2][2] = {{{0,0},{1,0}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{2,2},{2,2}}, {{2,2}}};
    int global_src_mask[12] = {1, 1, 1, 0,
                               1, 1, 1, 1,
                               0, 1, 1, 1};
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
    int * src_mask = NULL;
    if (!is_tgt) {
      src_mask = xmalloc(grid_data.num_vertices * sizeof(*src_mask));
      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        src_mask[i] = global_src_mask[grid_data.vertex_ids[i]];
      yac_basic_grid_add_mask_nocpy(grids[0], CORNER, src_mask, NULL);
    }

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);
    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       yac_interp_method_fixed_new(-1.0), NULL};

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
        double * ref_tgt_field = NULL;
        double ref_global_tgt_field[9] =
          {2.5, 3.5, -1, -1, 7.5, 8.5, -1, -1, -1};

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
          ref_tgt_field = xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
            ref_tgt_field[i] = ref_global_tgt_field[grid_data.vertex_ids[i]];
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
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  {// Test 3
   // linear distance weighted (no partial coverage) point interpolation
   // with two source processes and a single target process

    // the global grid is a 3x3 grid:
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
    // source mask:
    //       0-------1-------1-------0
    //       |       |       |       |
    //       |       |       |       |
    //       |       |       |       |
    //       1-------1-------1-------1
    //       |       |       |       |
    //       |       |       |       |
    //       |       |       |       |
    //       1-------1-------1-------1
    //       |       |       |       |
    //       |       |       |       |
    //       |       |       |       |
    //       0-------1-------1-------0
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][4] = {{0.0,1.0,2.0,3.0}, {0.0,1.0,2.0,3.0}};
    double coordinates_y[2][4] = {{0.0,1.0,2.0,3.0}, {0.0,1.0,2.0,3.0}};
    size_t const num_cells[2][2] = {{3,3}, {3,3}};
    size_t local_start[2][2][2] = {{{0,0},{1,0}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{2,3},{2,3}}, {{3,3}}};
    int global_src_mask[16] = {0, 1, 1, 0,
                               1, 1, 1, 1,
                               1, 1, 1, 1,
                               0, 1, 1, 0};
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

    int * src_mask = NULL;
    if (!is_tgt) {
      src_mask = xmalloc(grid_data.num_vertices * sizeof(*src_mask));
      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        src_mask[i] = global_src_mask[grid_data.vertex_ids[i]];
      yac_basic_grid_add_mask_nocpy(grids[0], CORNER, src_mask, NULL);
    }

    coordinate_pointer tgt_cell_coordinates = NULL;
    if (is_tgt) {
      double cell_coordinates_x[3] = {0.75,1.75,2.75};
      double cell_coordinates_y[3] = {0.75,1.75,2.75};
      tgt_cell_coordinates = xmalloc(9 * sizeof(*tgt_cell_coordinates));
      for (int i = 0, k = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j, ++k)
          LLtoXYZ_deg(
            cell_coordinates_x[j], cell_coordinates_y[i],
            tgt_cell_coordinates[k]);
      yac_basic_grid_add_coordinates_nocpy(grids[0], CELL, tgt_cell_coordinates);
    }

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_avg_new(AVG_DIST, 0),
       yac_interp_method_fixed_new(-1.0), NULL};

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
        double * ref_tgt_field = NULL;
        double ref_global_tgt_field[9] =
          {-1.0,//0.14963* 0.0 + 0.20075*( 1.0+ 4.0) + 0.44888* 5.0,
           0.14963* 1.0 + 0.20075*( 2.0+ 5.0) + 0.44888* 6.0,
           -1.0,//0.14963* 2.0 + 0.20075*( 3.0+ 6.0) + 0.44888* 7.0,
           0.14963* 4.0 + 0.20075*( 5.0+ 8.0) + 0.44888* 9.0,
           0.14963* 5.0 + 0.20075*( 6.0+ 9.0) + 0.44888*10.0,
           0.14963* 6.0 + 0.20075*( 7.0+10.0) + 0.44888*11.0,
           -1.0,//0.14963* 8.0 + 0.20075*( 9.0+12.0) + 0.44888*13.0,
           0.14963* 9.0 + 0.20075*(10.0+13.0) + 0.44888*14.0,
          -1.0};//0.14963*10.0 + 0.20075*(11.0+14.0) + 0.44888*15.0};

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
          ref_tgt_field = xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
            ref_tgt_field[i] = ref_global_tgt_field[grid_data.vertex_ids[i]];
        } else {
          src_field = xmalloc(grid_data.num_vertices * sizeof(*src_field));
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
            src_field[i] = (double)(grid_data.vertex_ids[i]);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        if (is_tgt)
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            if (fabs(ref_tgt_field[i] - tgt_field[i]) > 1e-3)
              PUT_ERR("wrong interpolation result");

        free(src_field);
        free(tgt_field);
        free(ref_tgt_field);
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

  {// Test 4
   // linear distance weighted (with partial coverage) point interpolation
   // with two source processes and a single target process

    // the global grid is a 3x3 grid:
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
    // source mask:
    //       0-------1-------1-------0
    //       |       |       |       |
    //       |       |       |       |
    //       |       |       |       |
    //       1-------0-------1-------1
    //       |       |       |       |
    //       |       |       |       |
    //       |       |       |       |
    //       1-------1-------0-------0
    //       |       |       |       |
    //       |       |       |       |
    //       |       |       |       |
    //       0-------1-------0-------0
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][4] = {{0.0,1.0,2.0,3.0}, {0.0,1.0,2.0,3.0}};
    double coordinates_y[2][4] = {{0.0,1.0,2.0,3.0}, {0.0,1.0,2.0,3.0}};
    size_t const num_cells[2][2] = {{3,3}, {3,3}};
    size_t local_start[2][2][2] = {{{0,0},{1,0}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{2,3},{2,3}}, {{3,3}}};
    int global_src_mask[16] = {0, 1, 0, 0,
                               1, 1, 0, 0,
                               1, 0, 1, 1,
                               0, 1, 1, 0};
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

    int * src_mask = NULL;
    if (!is_tgt) {
      src_mask = xmalloc(grid_data.num_vertices * sizeof(*src_mask));
      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        src_mask[i] = global_src_mask[grid_data.vertex_ids[i]];
      yac_basic_grid_add_mask_nocpy(grids[0], CORNER, src_mask, NULL);
    }

    coordinate_pointer tgt_cell_coordinates = NULL;
    if (is_tgt) {
      double cell_coordinates_x[3] = {0.75,1.75,2.75};
      double cell_coordinates_y[3] = {0.75,1.75,2.75};
      tgt_cell_coordinates = xmalloc(9 * sizeof(*tgt_cell_coordinates));
      for (int i = 0, k = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j, ++k)
          LLtoXYZ_deg(
            cell_coordinates_x[j], cell_coordinates_y[i],
            tgt_cell_coordinates[k]);
      yac_basic_grid_add_coordinates_nocpy(grids[0], CELL, tgt_cell_coordinates);
    }

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_avg_new(AVG_DIST, 1),
       yac_interp_method_fixed_new(-1.0), NULL};

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
        double * ref_tgt_field = NULL;
        double inv_dist_a = sqrt(8.0);
        double inv_dist_b = sqrt(8.0/5.0);
        double inv_dist_c = sqrt(8.0/9.0);
        double ref_global_tgt_field[9] =
          {(inv_dist_b*( 1.0+ 4.0) + inv_dist_a* 5.0)/(inv_dist_a + 2.0 * inv_dist_b),
           (inv_dist_c* 1.0 + inv_dist_b* 5.0)/(inv_dist_b + inv_dist_c),
           -1.0,
           (inv_dist_c* 4.0 + inv_dist_b*( 5.0+ 8.0))/(2.0 * inv_dist_b + inv_dist_c),
           (inv_dist_c* 5.0 + inv_dist_a*10.0)/(inv_dist_a + inv_dist_c),
           (inv_dist_b*10.0 + inv_dist_a*11.0)/(inv_dist_a + inv_dist_b),
           (inv_dist_c* 8.0 + inv_dist_a*13.0)/(inv_dist_a + inv_dist_c),
           (inv_dist_b*(10.0+13.0) + inv_dist_a*14.0)/(inv_dist_a + 2.0 * inv_dist_b),
           (inv_dist_c*10.0 + inv_dist_b*(11.0+14.0))/(2.0 * inv_dist_b + inv_dist_c)};

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
          ref_tgt_field = xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
            ref_tgt_field[i] = ref_global_tgt_field[grid_data.vertex_ids[i]];
        } else {
          src_field = xmalloc(grid_data.num_vertices * sizeof(*src_field));
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
            src_field[i] = (double)(grid_data.vertex_ids[i]);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        if (is_tgt)
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            if (fabs(ref_tgt_field[i] - tgt_field[i]) > 1e-3)
              PUT_ERR("wrong interpolation result");

        free(src_field);
        free(tgt_field);
        free(ref_tgt_field);
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

  {// Test 5
   // linear distance weighted (with partial coverage) point interpolation
   // with two source processes and a single target process
   // (target points are on vertices of the source grid

    // the global grid is a 3x3 grid:
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
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][4] = {{0.0,1.0,2.0,3.0}, {0.0,1.0,2.0,3.0}};
    double coordinates_y[2][4] = {{0.0,1.0,2.0,3.0}, {0.0,1.0,2.0,3.0}};
    size_t const num_cells[2][2] = {{3,3}, {3,3}};
    size_t local_start[2][2][2] = {{{0,0},{1,0}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{2,3},{2,3}}, {{3,3}}};
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

    coordinate_pointer tgt_cell_coordinates = NULL;
    if (is_tgt) {
      double cell_coordinates_x[3] = {1.0,2.0,3.0};
      double cell_coordinates_y[3] = {1.0,2.0,3.0};
      tgt_cell_coordinates = xmalloc(9 * sizeof(*tgt_cell_coordinates));
      for (int i = 0, k = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j, ++k)
          LLtoXYZ_deg(
            cell_coordinates_x[j], cell_coordinates_y[i],
            tgt_cell_coordinates[k]);
      yac_basic_grid_add_coordinates_nocpy(grids[0], CELL, tgt_cell_coordinates);
    }

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct interp_field src_fields[] =
      {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_avg_new(AVG_DIST, 1),
       yac_interp_method_fixed_new(-1.0), NULL};

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
        double * ref_tgt_field = NULL;
        double ref_global_tgt_field[9] =
          {5.0, 6.0, 7.0, 9.0, 10.0, 11.0, 13.0, 14.0, 15.0};

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
          ref_tgt_field = xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
            ref_tgt_field[i] = ref_global_tgt_field[grid_data.vertex_ids[i]];
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
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  {// Test 6
   // tests cell based source field

    // the global grid is a 7x7 grid:
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
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[2][2][2] = {{{0,0},{0,4}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{7,4},{7,3}}, {{7,7}}};
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

    coordinate_pointer src_point_coordinates = NULL;
    if (!is_tgt) {
      src_point_coordinates =
        xmalloc(grid_data.num_cells * sizeof(*src_point_coordinates));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = src_point_coordinates[i];
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
      }
      yac_basic_grid_add_coordinates_nocpy(grids[0], CELL, src_point_coordinates);
    }

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);
    struct interp_field src_fields[] =
      {{.location = CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct interp_field tgt_field =
      {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_avg_new(AVG_ARITHMETIC, 0),
       yac_interp_method_fixed_new(-1.0), NULL};

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
        double * ref_tgt_field = NULL;
        double ref_global_tgt_field[64] =
          {-4,-4,-4,-4,-4,-4,-4,-4,
           -4, 0+ 1+ 7+ 8, 1+ 2+ 8+ 9, 2+ 3+ 9+10, 3+ 4+10+11, 4+ 5+11+12, 5+ 6+12+13,-4,
           -4, 7+ 8+14+15, 8+ 9+15+16, 9+10+16+17,10+11+17+18,11+12+18+19,12+13+19+20,-4,
           -4,14+15+21+22,15+16+22+23,16+17+23+24,17+18+24+25,18+19+25+26,19+20+26+27,-4,
           -4,21+22+28+29,22+23+29+30,23+24+30+31,24+25+31+32,25+26+32+33,26+27+33+34,-4,
           -4,28+29+35+36,29+30+36+37,30+31+37+38,31+32+38+39,32+33+39+40,33+34+40+41,-4,
           -4,35+36+42+43,36+37+43+44,37+38+44+45,38+39+45+46,39+40+46+47,40+41+47+48,-4,
           -4,-4,-4,-4,-4,-4,-4,-4};
        for (size_t i = 0;
             i < sizeof(ref_global_tgt_field) / sizeof(ref_global_tgt_field[0]);
             ++i)
          ref_global_tgt_field[i] /= 4.0;

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
          ref_tgt_field = xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
            ref_tgt_field[i] = ref_global_tgt_field[grid_data.vertex_ids[i]];
        } else {
          src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            src_field[i] = (double)(grid_data.cell_ids[i]);
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
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  {// Test 7
   // barycentric coordinates weighted (with partial coverage) point
   // interpolation with two source processes and a single target process

    // the global source grid is a 1x1 grid:
    //       02--03--03
    //       |       |
    //       01  00  02
    //       |       |
    //       00--00--01
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double src_coordinates_x[2][2] = {{-0.5,0.5}, {0.5,-0.5}};
    double src_coordinates_y[2][2][2] =
      {{{-0.5,0.5}, {0.5,-0.5}}, {{89,90}, {90,89}}};
    double src_coordinates_x_unstruct[4];
    double src_coordinates_y_unstruct[4];
    size_t const src_num_cells[2] = {1,1};
    size_t src_local_start[2] = {0,0};
    size_t src_local_count[2] = {1,1};
    int src_cell_mask[4] = {1,1, 0,1};
    yac_int src_reorder[2][2][4] =
      {{{0,1,2,3},{2,3,0,1}},{{1,0,3,2},{3,2,1,0}}};
    double tgt_coordinates_x[5] = {-0.5,-0.25,0.0,0.25,0.5};
    double tgt_coordinates_y[2][5] =
      {{-0.5,-0.25,0.0,0.25,0.5}, {89.00,89.25,89.50,89.75,90.00}};
    size_t const tgt_num_cells[2] = {4,4};
    size_t tgt_local_start[2] = {0,0};
    size_t tgt_local_count[2] = {4,4};

    for (size_t i = 0; i < 2 * 2; ++i)
      (&(src_coordinates_x[0][0]))[i] *= YAC_RAD;
    for (size_t i = 0; i < 2 * 2 * 2; ++i)
      (&(src_coordinates_y[0][0][0]))[i] *= YAC_RAD;
    for (size_t i = 0; i < 5; ++i) tgt_coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i < 2 * 5; ++i)
      (&(tgt_coordinates_y[0][0]))[i] *= YAC_RAD;

    for (int partial_coverage = 0; partial_coverage < 2; ++partial_coverage) {
      for (int is_reg = 0; is_reg < 2; ++is_reg) {
        for (int at_pole = 0; at_pole <= is_reg; ++at_pole) {
          for (int with_mask = 0; with_mask < 2; ++with_mask) {
            for (int order_x = 0; order_x < 2; ++order_x) {
              for (int order_y = 0; order_y < 2; ++order_y) {

                struct basic_grid_data grid_data;

                if (is_tgt) {
                  grid_data =
                    yac_generate_basic_grid_data_reg2d(
                      tgt_coordinates_x, tgt_coordinates_y[at_pole],
                      tgt_num_cells, tgt_local_start, tgt_local_count, 0);
                } else {
                  if (is_reg) {
                    grid_data =
                      yac_generate_basic_grid_data_reg2d(
                        src_coordinates_x[order_x],
                        src_coordinates_y[at_pole][order_y], src_num_cells,
                        src_local_start, src_local_count, 0);
                  } else {
                    for (int i = 0; i < 4; ++i) {
                      src_coordinates_x_unstruct[i] =
                        src_coordinates_x[order_x][i&1];
                      src_coordinates_y_unstruct[i] =
                        src_coordinates_y[at_pole][order_y][i>>1];
                    }
                    grid_data =
                      yac_generate_basic_grid_data_unstruct(
                        4, 1, (int[]){4},
                        src_coordinates_x_unstruct,
                        src_coordinates_y_unstruct, (int[]){0,1,3,2});
                  }
                }

                struct yac_basic_grid * grids[2] =
                  {yac_basic_grid_new(grid_names[is_tgt], grid_data),
                   yac_basic_grid_empty_new(grid_names[is_tgt^1])};

                if (with_mask && !is_tgt) {
                  int * mask = xmalloc(grid_data.num_vertices * sizeof(*mask));
                  for (int i = 0; i < grid_data.num_vertices; ++i)
                    mask[i] =
                      src_cell_mask[src_reorder[order_x][order_y][i]];
                  yac_basic_grid_add_mask_nocpy(grids[0], CORNER, mask, NULL);
                }

                struct dist_grid_pair * grid_pair =
                  yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

                struct interp_field src_fields[] =
                  {{.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = (with_mask)?0:SIZE_MAX}};
                size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
                struct interp_field tgt_field =
                  {.location = CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

                struct interp_grid * interp_grid =
                  yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                                      num_src_fields, src_fields, tgt_field);

                struct interp_method * method_stack[] =
                  {yac_interp_method_avg_new(AVG_BARY, partial_coverage),
                   yac_interp_method_fixed_new(-1.0), NULL};

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
                    double * ref_tgt_field = NULL;
                    double ref_global_tgt_field[25];
                    double ref_weights[2][2][25][4] =
                      {{{{1.00, 0.00, 0.00, 0.00},
                         {0.75, 0.25, 0.00, 0.00},
                         {0.50, 0.50, 0.00, 0.00},
                         {0.25, 0.75, 0.00, 0.00},
                         {0.00, 1.00, 0.00, 0.00},
                         {0.75, 0.00, 0.25, 0.00},
                         {0.75, 0.00, 0.00, 0.25},
                         {0.50, 0.25, 0.00, 0.25},
                         {0.25, 0.50, 0.00, 0.25},
                         {0.00, 0.75, 0.00, 0.25},
                         {0.50, 0.00, 0.50, 0.00},
                         {0.50, 0.00, 0.25, 0.25},
                         {0.50, 0.00, 0.00, 0.50},
                         {0.25, 0.25, 0.00, 0.50},
                         {0.00, 0.50, 0.00, 0.50},
                         {0.25, 0.00, 0.75, 0.00},
                         {0.25, 0.00, 0.50, 0.25},
                         {0.25, 0.00, 0.25, 0.50},
                         {0.25, 0.00, 0.00, 0.75},
                         {0.00, 0.25, 0.00, 0.75},
                         {0.00, 0.00, 1.00, 0.00},
                         {0.00, 0.00, 0.75, 0.25},
                         {0.00, 0.00, 0.50, 0.50},
                         {0.00, 0.00, 0.25, 0.75},
                         {0.00, 0.00, 0.00, 1.00}},
                        {{1.00, 0.00, 0.00, 0.00},
                         {0.75, 0.25, 0.00, 0.00},
                         {0.50, 0.50, 0.00, 0.00},
                         {0.25, 0.75, 0.00, 0.00},
                         {0.00, 1.00, 0.00, 0.00},
                         {0.00, 0.00, 0.00, 0.00},//{0.75, 0.00, 0.25, 0.00}
                         {0.75, 0.00, 0.00, 0.25},
                         {0.50, 0.25, 0.00, 0.25},
                         {0.25, 0.50, 0.00, 0.25},
                         {0.00, 0.75, 0.00, 0.25},
                         {0.00, 0.00, 0.00, 0.00},//{0.50, 0.00, 0.50, 0.00}
                         {0.00, 0.00, 0.00, 0.00},//{0.50, 0.00, 0.25, 0.25}
                         {0.50, 0.00, 0.00, 0.50},
                         {0.25, 0.25, 0.00, 0.50},
                         {0.00, 0.50, 0.00, 0.50},
                         {0.00, 0.00, 0.00, 0.00},//{0.25, 0.00, 0.75, 0.00}
                         {0.00, 0.00, 0.00, 0.00},//{0.25, 0.00, 0.50, 0.25}
                         {0.00, 0.00, 0.00, 0.00},//{0.25, 0.00, 0.25, 0.50}
                         {0.25, 0.00, 0.00, 0.75},
                         {0.00, 0.25, 0.00, 0.75},
                         {0.00, 0.00, 0.00, 0.00},//{0.00, 0.00, 1.00, 0.00}
                         {0.00, 0.00, 0.00, 0.00},//{0.00, 0.00, 0.75, 0.25}
                         {0.00, 0.00, 0.00, 0.00},//{0.00, 0.00, 0.50, 0.50}
                         {0.00, 0.00, 0.00, 0.00},//{0.00, 0.00, 0.25, 0.75}
                         {0.00, 0.00, 0.00, 1.00}}},
                       {{{1.00, 0.00, 0.00, 0.00},
                         {0.75, 0.25, 0.00, 0.00},
                         {0.50, 0.50, 0.00, 0.00},
                         {0.25, 0.75, 0.00, 0.00},
                         {0.00, 1.00, 0.00, 0.00},
                         {0.75, 0.00, 0.25, 0.00},
                         {0.75, 0.00, 0.00, 0.25},
                         {0.50, 0.25, 0.00, 0.25},
                         {0.25, 0.50, 0.00, 0.25},
                         {0.00, 0.75, 0.00, 0.25},
                         {0.50, 0.00, 0.50, 0.00},
                         {0.50, 0.00, 0.25, 0.25},
                         {0.50, 0.00, 0.00, 0.50},
                         {0.25, 0.25, 0.00, 0.50},
                         {0.00, 0.50, 0.00, 0.50},
                         {0.25, 0.00, 0.75, 0.00},
                         {0.25, 0.00, 0.50, 0.25},
                         {0.25, 0.00, 0.25, 0.50},
                         {0.25, 0.00, 0.00, 0.75},
                         {0.00, 0.25, 0.00, 0.75},
                         {0.00, 0.00, 1.00, 0.00},
                         {0.00, 0.00, 0.75, 0.25},
                         {0.00, 0.00, 0.50, 0.50},
                         {0.00, 0.00, 0.25, 0.75},
                         {0.00, 0.00, 0.00, 1.00}},
                        {{1.00, 0.00, 0.00, 0.00},
                         {0.75, 0.25, 0.00, 0.00},
                         {0.50, 0.50, 0.00, 0.00},
                         {0.25, 0.75, 0.00, 0.00},
                         {0.00, 1.00, 0.00, 0.00},
                         {1.00, 0.00, 0.00, 0.00},//{0.75, 0.00, 0.25, 0.00}
                         {0.75, 0.00, 0.00, 0.25},
                         {0.50, 0.25, 0.00, 0.25},
                         {0.25, 0.50, 0.00, 0.25},
                         {0.00, 0.75, 0.00, 0.25},
                         {1.00, 0.00, 0.00, 0.00},//{0.50, 0.00, 0.50, 0.00}
                         {0.50/0.75, 0.00, 0.00, 0.25/0.75},//{0.50, 0.00, 0.25, 0.25}
                         {0.50, 0.00, 0.00, 0.50},
                         {0.25, 0.25, 0.00, 0.50},
                         {0.00, 0.50, 0.00, 0.50},
                         {1.00, 0.00, 0.00, 0.00},//{0.25, 0.00, 0.75, 0.00}
                         {0.50, 0.00, 0.00, 0.50},//{0.25, 0.00, 0.50, 0.25}
                         {0.25/0.75, 0.00, 0.00, 0.50/0.75},//{0.25, 0.00, 0.25, 0.50}
                         {0.25, 0.00, 0.00, 0.75},
                         {0.00, 0.25, 0.00, 0.75},
                         {0.00, 0.00, 0.00, 0.00},//{0.00, 0.00, 1.00, 0.00}
                         {0.00, 0.00, 0.00, 1.00},//{0.00, 0.00, 0.75, 0.25}
                         {0.00, 0.00, 0.00, 1.00},//{0.00, 0.00, 0.50, 0.50}
                         {0.00, 0.00, 0.00, 1.00},//{0.00, 0.00, 0.25, 0.75}
                         {0.00, 0.00, 0.00, 1.00}}}};

                    for (size_t j = 0; j < 25; ++j) {
                      ref_global_tgt_field[j] =
                        (double)(src_reorder[order_x][order_y][0] + 1) *
                        ref_weights[partial_coverage][with_mask][j][0] +
                        (double)(src_reorder[order_x][order_y][1] + 1) *
                        ref_weights[partial_coverage][with_mask][j][1] +
                        (double)(src_reorder[order_x][order_y][2] + 1) *
                        ref_weights[partial_coverage][with_mask][j][2] +
                        (double)(src_reorder[order_x][order_y][3] + 1) *
                        ref_weights[partial_coverage][with_mask][j][3];
                      if (ref_global_tgt_field[j] == 0.0)
                        ref_global_tgt_field[j] = -1.0;
                    }

                    if (is_tgt) {
                      tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
                      ref_tgt_field = xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
                      for (size_t i = 0; i < grid_data.num_vertices; ++i)
                        ref_tgt_field[i] = ref_global_tgt_field[i];
                    } else {
                      src_field = xmalloc(grid_data.num_vertices * sizeof(*src_field));
                      for (size_t i = 0; i < grid_data.num_vertices; ++i)
                        src_field[i] = (double)(i) + 1.0;
                    }

                    yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

                    if (is_tgt)
                      for (size_t i = 0; i < grid_data.num_vertices; ++i)
                        if (fabs(ref_tgt_field[i] - tgt_field[i]) > 1e-3)
                          PUT_ERR("wrong interpolation result");

                    free(src_field);
                    free(tgt_field);
                    free(ref_tgt_field);
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
            }
          }
        }
      }
    }
  }

  {// Test 8
   // barycentric coordinates weighted (with partial coverage) point
   // interpolation with two source processes and a single target process

    // the global source grid:
    //       02------03
    //       |     / |
    //       |   /   |
    //       | /     |
    //       00------01
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    size_t src_num_vertices = 4;
    size_t src_num_cells = 2;
    int src_num_vertices_per_cell[] = {3, 3};
    int src_cell_to_vertex[] = {0,1,3, 0,2,3};
    double src_coordinates_x[4] = {-0.5,0.5,-0.5,0.5};
    double src_coordinates_y[4] = {-0.5,-0.5,0.5,0.5};
    size_t tgt_num_vertices[2] = {5,5};
    int tgt_cyclic[2] = {0,0};
    double tgt_coordinates_x[5] = {-0.5,-0.25,0.0,0.25,0.5};
    double tgt_coordinates_y[5] = {-0.5,-0.25,0.0,0.25,0.5};
    int src_cell_mask[4] = {1,1, 0,1};

    for (int with_mask = 0; with_mask < 2; ++with_mask) {

      struct basic_grid_data grid_data;
      if (is_tgt) {
        grid_data =
          yac_generate_basic_grid_data_reg_2d_deg(
            tgt_num_vertices, tgt_cyclic,
            tgt_coordinates_x, tgt_coordinates_y);
      } else {
        grid_data =
          yac_generate_basic_grid_data_unstruct_deg(
            src_num_vertices, src_num_cells, src_num_vertices_per_cell,
            src_coordinates_x, src_coordinates_y, src_cell_to_vertex);
      }

      struct yac_basic_grid * grids[2] =
        {yac_basic_grid_new(grid_names[is_tgt], grid_data),
         yac_basic_grid_empty_new(grid_names[is_tgt^1])};

      if (with_mask && !is_tgt) {
        int * mask = xmalloc(grid_data.num_vertices * sizeof(*mask));
        for (int i = 0; i < grid_data.num_vertices; ++i)
          mask[i] = src_cell_mask[i];
        yac_basic_grid_add_mask_nocpy(grids[0], CORNER, mask, NULL);
      }

      struct dist_grid_pair * grid_pair =
        yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

      struct interp_field src_fields[] =
        {{.location = CORNER, .coordinates_idx = SIZE_MAX,
          .masks_idx = (with_mask)?0:SIZE_MAX}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct interp_field tgt_field =
        {.location = CORNER, .coordinates_idx = SIZE_MAX,
         .masks_idx = SIZE_MAX};

      struct interp_grid * interp_grid =
        yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                            num_src_fields, src_fields, tgt_field);

      for (int partial_coverage = 0; partial_coverage < 2; ++partial_coverage) {

        struct interp_method * method_stack[] =
          {yac_interp_method_avg_new(AVG_BARY, partial_coverage),
           yac_interp_method_fixed_new(-1.0), NULL};

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
            double * ref_tgt_field = NULL;
            double ref_global_tgt_field[25];
            double ref_weights[2][2][25][4] =
             {{{{1.00, 0.00, 0.00, 0.00},
                 {0.75, 0.25, 0.00, 0.00},
                 {0.50, 0.50, 0.00, 0.00},
                 {0.25, 0.75, 0.00, 0.00},
                 {0.00, 1.00, 0.00, 0.00},
                 {0.75, 0.00, 0.25, 0.00},
                 {0.75, 0.00, 0.00, 0.25},
                 {0.50, 0.25, 0.00, 0.25},
                 {0.25, 0.50, 0.00, 0.25},
                 {0.00, 0.75, 0.00, 0.25},
                 {0.50, 0.00, 0.50, 0.00},
                 {0.50, 0.00, 0.25, 0.25},
                 {0.50, 0.00, 0.00, 0.50},
                 {0.25, 0.25, 0.00, 0.50},
                 {0.00, 0.50, 0.00, 0.50},
                 {0.25, 0.00, 0.75, 0.00},
                 {0.25, 0.00, 0.50, 0.25},
                 {0.25, 0.00, 0.25, 0.50},
                 {0.25, 0.00, 0.00, 0.75},
                 {0.00, 0.25, 0.00, 0.75},
                 {0.00, 0.00, 1.00, 0.00},
                 {0.00, 0.00, 0.75, 0.25},
                 {0.00, 0.00, 0.50, 0.50},
                 {0.00, 0.00, 0.25, 0.75},
                 {0.00, 0.00, 0.00, 1.00}},
                {{1.00, 0.00, 0.00, 0.00},
                 {0.75, 0.25, 0.00, 0.00},
                 {0.50, 0.50, 0.00, 0.00},
                 {0.25, 0.75, 0.00, 0.00},
                 {0.00, 1.00, 0.00, 0.00},
                 {0.00, 0.00, 0.00, 0.00},//{0.75, 0.00, 0.25, 0.00}
                 {0.75, 0.00, 0.00, 0.25},
                 {0.50, 0.25, 0.00, 0.25},
                 {0.25, 0.50, 0.00, 0.25},
                 {0.00, 0.75, 0.00, 0.25},
                 {0.00, 0.00, 0.00, 0.00},//{0.50, 0.00, 0.50, 0.00}
                 {0.00, 0.00, 0.00, 0.00},//{0.50, 0.00, 0.25, 0.25}
                 {0.50, 0.00, 0.00, 0.50},
                 {0.25, 0.25, 0.00, 0.50},
                 {0.00, 0.50, 0.00, 0.50},
                 {0.00, 0.00, 0.00, 0.00},//{0.25, 0.00, 0.75, 0.00}
                 {0.00, 0.00, 0.00, 0.00},//{0.25, 0.00, 0.50, 0.25}
                 {0.00, 0.00, 0.00, 0.00},//{0.25, 0.00, 0.25, 0.50}
                 {0.25, 0.00, 0.00, 0.75},
                 {0.00, 0.25, 0.00, 0.75},
                 {0.00, 0.00, 0.00, 0.00},//{0.00, 0.00, 1.00, 0.00}
                 {0.00, 0.00, 0.00, 0.00},//{0.00, 0.00, 0.75, 0.25}
                 {0.00, 0.00, 0.00, 0.00},//{0.00, 0.00, 0.50, 0.50}
                 {0.00, 0.00, 0.00, 0.00},//{0.00, 0.00, 0.25, 0.75}
                 {0.00, 0.00, 0.00, 1.00}}},
               {{{1.00, 0.00, 0.00, 0.00},
                 {0.75, 0.25, 0.00, 0.00},
                 {0.50, 0.50, 0.00, 0.00},
                 {0.25, 0.75, 0.00, 0.00},
                 {0.00, 1.00, 0.00, 0.00},
                 {0.75, 0.00, 0.25, 0.00},
                 {0.75, 0.00, 0.00, 0.25},
                 {0.50, 0.25, 0.00, 0.25},
                 {0.25, 0.50, 0.00, 0.25},
                 {0.00, 0.75, 0.00, 0.25},
                 {0.50, 0.00, 0.50, 0.00},
                 {0.50, 0.00, 0.25, 0.25},
                 {0.50, 0.00, 0.00, 0.50},
                 {0.25, 0.25, 0.00, 0.50},
                 {0.00, 0.50, 0.00, 0.50},
                 {0.25, 0.00, 0.75, 0.00},
                 {0.25, 0.00, 0.50, 0.25},
                 {0.25, 0.00, 0.25, 0.50},
                 {0.25, 0.00, 0.00, 0.75},
                 {0.00, 0.25, 0.00, 0.75},
                 {0.00, 0.00, 1.00, 0.00},
                 {0.00, 0.00, 0.75, 0.25},
                 {0.00, 0.00, 0.50, 0.50},
                 {0.00, 0.00, 0.25, 0.75},
                 {0.00, 0.00, 0.00, 1.00}},
                {{1.00, 0.00, 0.00, 0.00},
                 {0.75, 0.25, 0.00, 0.00},
                 {0.50, 0.50, 0.00, 0.00},
                 {0.25, 0.75, 0.00, 0.00},
                 {0.00, 1.00, 0.00, 0.00},
                 {1.00, 0.00, 0.00, 0.00},//{0.75, 0.00, 0.25, 0.00}
                 {0.75, 0.00, 0.00, 0.25},
                 {0.50, 0.25, 0.00, 0.25},
                 {0.25, 0.50, 0.00, 0.25},
                 {0.00, 0.75, 0.00, 0.25},
                 {1.00, 0.00, 0.00, 0.00},//{0.50, 0.00, 0.50, 0.00}
                 {0.50/0.75, 0.00, 0.00, 0.25/0.75},//{0.50, 0.00, 0.25, 0.25}
                 {0.50, 0.00, 0.00, 0.50},
                 {0.25, 0.25, 0.00, 0.50},
                 {0.00, 0.50, 0.00, 0.50},
                 {1.00, 0.00, 0.00, 0.00},//{0.25, 0.00, 0.75, 0.00}
                 {0.50, 0.00, 0.00, 0.50},//{0.25, 0.00, 0.50, 0.25}
                 {0.25/0.75, 0.00, 0.00, 0.50/0.75},//{0.25, 0.00, 0.25, 0.50}
                 {0.25, 0.00, 0.00, 0.75},
                 {0.00, 0.25, 0.00, 0.75},
                 {0.00, 0.00, 0.00, 0.00},//{0.00, 0.00, 1.00, 0.00}
                 {0.00, 0.00, 0.00, 1.00},//{0.00, 0.00, 0.75, 0.25}
                 {0.00, 0.00, 0.00, 1.00},//{0.00, 0.00, 0.50, 0.50}
                 {0.00, 0.00, 0.00, 1.00},//{0.00, 0.00, 0.25, 0.75}
                 {0.00, 0.00, 0.00, 1.00}}}};

            for (size_t j = 0; j < 25; ++j) {
              ref_global_tgt_field[j] =
                1.0 * ref_weights[partial_coverage][with_mask][j][0] +
                2.0 * ref_weights[partial_coverage][with_mask][j][1] +
                3.0 * ref_weights[partial_coverage][with_mask][j][2] +
                4.0 * ref_weights[partial_coverage][with_mask][j][3];
              if (ref_global_tgt_field[j] == 0.0)
                ref_global_tgt_field[j] = -1.0;
            }

            if (is_tgt) {
              tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
              ref_tgt_field = xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
              for (size_t i = 0; i < grid_data.num_vertices; ++i)
                ref_tgt_field[i] = ref_global_tgt_field[i];
            } else {
              src_field = xmalloc(grid_data.num_vertices * sizeof(*src_field));
              for (size_t i = 0; i < grid_data.num_vertices; ++i)
                src_field[i] = (double)(i) + 1.0;
            }

            yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

            if (is_tgt)
              for (size_t i = 0; i < grid_data.num_vertices; ++i)
                if (fabs(ref_tgt_field[i] - tgt_field[i]) > 1e-3)
                  PUT_ERR("wrong interpolation result");

            free(src_field);
            free(tgt_field);
            free(ref_tgt_field);
          }

          yac_interpolation_delete(interpolation);
        }

        yac_interp_weights_delete(weights);
        yac_interp_method_delete(method_stack);
      } // partial_coverage
      yac_interp_grid_delete(interp_grid);
      yac_dist_grid_pair_delete(grid_pair);
      yac_basic_grid_delete(grids[1]);
      yac_basic_grid_delete(grids[0]);
    }
  }



  {// Test 9
   // simple interpolation test with various configuration options

    // the global source grid:
    //       06--10--07--11--08
    //       |       |       |
    //       06  02  08  03  09
    //       |       |       |
    //       03--05--04--07--05
    //       |       |       |
    //       01  00  03  01  04
    //       |       |       |
    //       00--00--01--02--02
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][7] = {{0.0,0.1,0.2}, {0.0,0.05,0.1,0.15,0.2}};
    double coordinates_y[2][7] = {{0.0,0.1,0.2}, {0.0,0.05,0.1,0.15,0.2}};
    size_t const num_cells[2][2] = {{2,2}, {4,4}};
    size_t local_start[2][2][2] = {{{0,0},{0,1}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{2,1},{2,1}}, {{4,4}}};
    int global_mask[2][5*5] = {{1,1,1,
                                0,1,1,
                                0,0,1},
                               {1,1,1,1,1,
                                1,1,1,1,0,
                                1,1,1,0,0,
                                1,1,0,0,0,
                                1,0,0,0,0}};
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

    int field_mask[5*5];
    for (size_t i = 0; i < grid_data.num_vertices; ++i)
      field_mask[i] = global_mask[is_tgt][grid_data.vertex_ids[i]];
    yac_basic_grid_add_mask(
      grids[0], CORNER, field_mask, grid_data.num_vertices, NULL);

    struct dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    for (int with_src_mask = 0; with_src_mask < 2; ++with_src_mask) {
      for (int with_tgt_mask = 0; with_tgt_mask < 2; ++with_tgt_mask) {

        struct interp_field src_fields[] =
          {{.location = CORNER,
            .coordinates_idx = SIZE_MAX,
            .masks_idx = (with_src_mask)?0:SIZE_MAX}};
        size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
        struct interp_field tgt_field =
          {.location = CORNER,
           .coordinates_idx = SIZE_MAX,
           .masks_idx = (with_tgt_mask)?0:SIZE_MAX};

        struct interp_grid * interp_grid =
          yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                              num_src_fields, src_fields, tgt_field);

        enum yac_interp_avg_weight_type weight_types[] =
          {AVG_ARITHMETIC, AVG_DIST, AVG_BARY};
        size_t const num_weight_types =
          sizeof(weight_types) / sizeof(weight_types[0]);

        for (size_t weight_types_idx = 0; weight_types_idx < num_weight_types;
             ++weight_types_idx) {
          for (int partial_coverage = 0; partial_coverage < 2;
               ++partial_coverage) {

            struct interp_method * method_stack[] =
              {yac_interp_method_avg_new(
                 weight_types[weight_types_idx], partial_coverage),
               yac_interp_method_fixed_new(-1.0), NULL};

            struct interp_weights * weights =
              yac_interp_method_do_search(method_stack, interp_grid);

            for (size_t i = 0; i < num_reorder_types; ++i) {

              struct interpolation * interpolation =
                yac_interp_weights_get_interpolation(
                  weights, reorder_types[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

              // check generated interpolation
              {
                double src_field_data[3*3];
                double * src_field = src_field_data;
                double ** src_fields = &src_field;
                double tgt_field_data[5*5] =
                  {-2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0};
                double * tgt_field = tgt_field_data;

                if (!is_tgt)
                  for (size_t i = 0; i < grid_data.num_vertices; ++i)
                    src_field[i] = (double)(grid_data.vertex_ids[i]);

                yac_interpolation_execute(
                  interpolation, &src_fields, &tgt_field);

                if (is_tgt) {
                  size_t src_per_tgt[5*5][4] =
                    {{0,1,3,4},{0,1,3,4},{0,1,3,4},{1,2,4,5},{1,2,4,5},
                     {0,1,3,4},{0,1,3,4},{0,1,3,4},{1,2,4,5},{1,2,4,5},
                     {0,1,3,4},{0,1,3,4},{0,1,3,4},{1,2,4,5},{1,2,4,5},
                     {3,4,6,7},{3,4,6,7},{3,4,6,7},{4,5,7,8},{4,5,7,8},
                     {3,4,6,7},{3,4,6,7},{3,4,6,7},{4,5,7,8},{4,5,7,8}};
                  double A = 0.345491502812526;
                  double B = 0.154508497187473;
                  double C = 0.25;
                  double D = 0.5;
                  double weights[3][5*5][4] =
                    {{{C,C,C,C},{C,C,C,C},{C,C,C,C},{C,C,C,C},{C,C,C,C},
                      {C,C,C,C},{C,C,C,C},{C,C,C,C},{C,C,C,C},{C,C,C,C},
                      {C,C,C,C},{C,C,C,C},{C,C,C,C},{C,C,C,C},{C,C,C,C},
                      {C,C,C,C},{C,C,C,C},{C,C,C,C},{C,C,C,C},{C,C,C,C},
                      {C,C,C,C},{C,C,C,C},{C,C,C,C},{C,C,C,C},{C,C,C,C}},
                     {{1,0,0,0},{A,A,B,B},{0,1,0,0},{A,A,B,B},{0,1,0,0},
                      {A,B,A,B},{C,C,C,C},{B,A,B,A},{C,C,C,C},{B,A,B,A},
                      {0,0,1,0},{B,B,A,A},{0,0,0,1},{B,B,A,A},{0,0,0,1},
                      {A,B,A,B},{C,C,C,C},{B,A,B,A},{C,C,C,C},{B,A,B,A},
                      {0,0,1,0},{B,B,A,A},{0,0,0,1},{B,B,A,A},{0,0,0,1}},
                     {{1,0,0,0},{D,D,0,0},{0,1,0,0},{D,D,0,0},{0,1,0,0},
                      {D,0,D,0},{D,0,0,D},{0,D,0,D},{D,0,0,D},{0,D,0,D},
                      {0,0,1,0},{0,0,D,D},{0,0,0,1},{0,0,D,D},{0,0,0,1},
                      {D,0,D,0},{D,0,0,D},{0,D,0,D},{D,0,0,D},{0,D,0,D},
                      {0,0,1,0},{0,0,D,D},{0,0,0,1},{0,0,D,D},{0,0,0,1}}};

                  for (size_t i = 0; i < grid_data.num_vertices; ++i) {

                    size_t tgt_idx = grid_data.vertex_ids[i];
                    double src_values[4];
                    double w[4];

                    double ref_tgt_value = 0.0;

                    int contains_masked_src = 0;
                    int contains_unmasked_src = 0;

                    for (int j = 0; j < 4; ++j) {
                      size_t src_idx = src_per_tgt[tgt_idx][j];
                      w[j] = weights[weight_types_idx][tgt_idx][j];
                      if (with_src_mask &&
                          !global_mask[0][src_idx] &&
                          (w[j] != 0.0)) {
                        src_values[j] = -1.0;
                        w[j] = 0.0;
                        contains_masked_src = 1;
                      } else {
                        src_values[j] = (double)src_idx;
                        contains_unmasked_src |= w[j] != 0.0;
                      }
                    }

                    if (with_tgt_mask && !global_mask[1][tgt_idx]) {
                      ref_tgt_value = -2.0;
                    } else if (contains_masked_src) {
                      if (partial_coverage && contains_unmasked_src) {
                        yac_correct_weights(4, w);
                        for (int j = 0; j < 4; ++j)
                          ref_tgt_value += w[j] * src_values[j];
                      } else {
                        ref_tgt_value = -1.0;
                      }
                    } else {
                      for (int j = 0; j < 4; ++j)
                        ref_tgt_value += w[j] * src_values[j];
                    }

                    // some special cases
                    if (with_src_mask && partial_coverage &&
                        (weight_types[weight_types_idx] == AVG_DIST) &&
                        (!with_tgt_mask || global_mask[1][tgt_idx])) {
                      switch (tgt_idx) {
                        case(10):
                          ref_tgt_value = 1.73879612503;
                          break;
                        case(20):
                        case(22):
                          ref_tgt_value = 4.0;
                        default:
                          break;
                      }
                    }

                    if (fabs(ref_tgt_value - tgt_field[i]) > 1e-3)
                      PUT_ERR("wrong interpolation result");
                  }
                }
              }

              yac_interpolation_delete(interpolation);
            } // reorder type
            yac_interp_weights_delete(weights);
            yac_interp_method_delete(method_stack);
          } // partial coverage
        } // weight type
        yac_interp_grid_delete(interp_grid);
      } // with tgt mask
    } // with src mask
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  yac_mpi_call(MPI_Comm_free(&split_comm), MPI_COMM_WORLD);

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
