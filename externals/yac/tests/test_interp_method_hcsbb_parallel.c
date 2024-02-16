/**
 * @file test_interp_method_hcsbb_parallel.c
 *
 * @copyright Copyright  (C)  2020 DKRZ, MPI-M
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
#include "interp_method_hcsbb.h"
#include "interp_method_fixed.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

enum interp_weights_reorder_type reorder_types[] =
  {MAPPING_ON_SRC, MAPPING_ON_TGT};
size_t num_reorder_types = sizeof(reorder_types) / sizeof(reorder_types[0]);

static void compute_field_data_XYZ(
  double (*vertices)[3], size_t num_vertices, double * f_v);
static char const * grid_names[2] = {"src_grid", "tgt_grid"};

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
    // case 1)
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
    // case 2)
    //
    // 0---0---1---1---3---3---4---4
    // | 0 | 1 | 1 | 3 | 3 | 4 | 4 |
    // 0---0---1---1---3---3---4---4
    // | 0 | 1 | 1 | 3 | 3 | 4 | 4 |
    // 0---0---1---1---3---3---4---4
    // | 0 | 1 | 1 | 3 | 3 | 4 | 4 |
    // 0---0---1---1---3---3---4---4
    // | 0 | 0 | 1 | 3 | 3 | 4 | 4 |
    // 0---0---1---1---2---2---4---4
    // | 0 | 0 | 1 | 2 | 2 | 4 | 4 |
    // 0---0---0---1---2---2---4---4
    // | 0 | 0 | 1 | 2 | 2 | 4 | 4 |
    // 0---0---0---1---2---2---4---4
    // | 0 | 0 | 1 | 2 | 2 | 4 | 4 |
    // 0---0---0---1---2---2---4---4
    //
    //---------------
    // setup
    //---------------

    enum {NUM_SRC_DECOMP = 2};
    size_t num_points[NUM_SRC_DECOMP][num_reorder_types];
    double * tgt_fields[NUM_SRC_DECOMP][num_reorder_types];
    int is_tgt = split_comm_size == 1;
    for (int case_idx = 0; case_idx < NUM_SRC_DECOMP; ++case_idx) {

      double coordinates_x[NUM_SRC_DECOMP][2][10] =
        {{{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},{-1.0}},
         {{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},{-1.0}}};
      double coordinates_y[NUM_SRC_DECOMP][2][10] =
        {{{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},{-1.0}},
         {{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},{-1.0}}};
      size_t const num_cells[2][2] = {{7,7}, {9,9}};
      size_t local_start[NUM_SRC_DECOMP][2][5][2] =
        {{{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}},
         {{{0,0},{1,0},{3,0},{3,3},{5,0}}, {{0,0}}}};
      size_t local_count[NUM_SRC_DECOMP][2][5][2] =
        {{{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{9,9}}},
         {{{2,7},{2,7},{2,3},{2,4},{2,7}}, {{9,9}}}};
      int with_halo = 0;
      for (size_t i = 0; i < 10; ++i) {
        coordinates_x[case_idx][1][i] = 2.5 + 2.0 * ((double)i / 9.0);
        coordinates_y[case_idx][1][i] = 2.5 + 2.0 * ((double)i / 9.0);
      }
      for (size_t i = 0; i <= num_cells[is_tgt][0]; ++i)
        coordinates_x[case_idx][is_tgt][i] *= YAC_RAD;
      for (size_t i = 0; i <= num_cells[is_tgt][1]; ++i)
        coordinates_y[case_idx][is_tgt][i] *= YAC_RAD;

      struct basic_grid_data grid_data =
        yac_generate_basic_grid_data_reg2d(
          coordinates_x[case_idx][is_tgt], coordinates_y[case_idx][is_tgt],
          num_cells[is_tgt], local_start[case_idx][is_tgt][split_comm_rank],
          local_count[case_idx][is_tgt][split_comm_rank], with_halo);

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

      struct interp_method * method_stack[] =
        {yac_interp_method_hcsbb_new(),
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

          if (is_tgt) {
            tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
            ref_tgt_field =
              xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
            compute_field_data_XYZ(
              grid_data.vertex_coordinates, 100, ref_tgt_field);
            for (size_t i = 0; i < 100; ++i) tgt_field[i] = -1;
            num_points[case_idx][i] = grid_data.num_vertices;
            tgt_fields[case_idx][i] = tgt_field;
          } else {
            src_field = xmalloc(grid_data.num_vertices * sizeof(*src_field));
            compute_field_data_XYZ(
              grid_data.vertex_coordinates, grid_data.num_vertices, src_field);
          }

          yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

          if (is_tgt)
            for (size_t i = 0; i < 100; ++i)
              if (fabs(1.0 - ref_tgt_field[i] / tgt_field[i]) > 1e-2)
                PUT_ERR("wrong interpolation result");

          free(ref_tgt_field);
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

    if (is_tgt) {
      for (int case_idx = 0; case_idx < NUM_SRC_DECOMP; ++case_idx) {
        for (int reorder_type = 0; reorder_type < num_reorder_types;
            ++reorder_type) {

          if (num_points[0][0] != num_points[case_idx][reorder_type])
            PUT_ERR("error in test setup");
          for (size_t i = 0; i < num_points[0][0]; ++i)
            if (tgt_fields[0][0][i] != tgt_fields[case_idx][reorder_type][i])
              PUT_ERR("interpolation results are not "
                      "decomposition-independant");
        }
      }
      for (int case_idx = 0; case_idx < NUM_SRC_DECOMP; ++case_idx)
        for (int reorder_type = 0; reorder_type < num_reorder_types;
            ++reorder_type) free(tgt_fields[case_idx][reorder_type]);
    }
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
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][10] = {{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},{-1.0}};
    double coordinates_y[2][10] = {{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},{-1.0}};
    size_t const num_cells[2][2] = {{7,7}, {9,9}};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{9,9}}};
    int with_halo = 0;
    for (size_t i = 0; i < 10; ++i) {
      coordinates_x[1][i] = 2.0 + 3.0 * ((double)i / 9.0);
      coordinates_y[1][i] = 2.0 + 3.0 * ((double)i / 9.0);
    }
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

    struct interp_method * method_stack[] =
      {yac_interp_method_hcsbb_new(),
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

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
          ref_tgt_field =
            xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
          compute_field_data_XYZ(
            grid_data.vertex_coordinates, 100, ref_tgt_field);
          for (size_t i = 0; i < 100; ++i) tgt_field[i] = -1;
        } else {
          src_field = xmalloc(grid_data.num_vertices * sizeof(*src_field));
          compute_field_data_XYZ(
            grid_data.vertex_coordinates, grid_data.num_vertices, src_field);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        if (is_tgt)
          for (size_t i = 0; i < 100; ++i)
            if (fabs(1.0 - ref_tgt_field[i] / tgt_field[i]) > 1e-2)
              PUT_ERR("wrong interpolation result");

        free(ref_tgt_field);
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
    // (source point location is at the cell centers)
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
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][10] = {{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},{-1.0}};
    double coordinates_y[2][10] = {{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},{-1.0}};
    size_t const num_cells[2][2] = {{7,7}, {9,9}};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{9,9}}};
    int with_halo = 0;
    for (size_t i = 0; i < 10; ++i) {
      coordinates_x[1][i] = 2.0 + 3.0 * ((double)i / 9.0);
      coordinates_y[1][i] = 2.0 + 3.0 * ((double)i / 9.0);
    }
    for (size_t i = 0; i <= num_cells[is_tgt][0]; ++i)
      coordinates_x[is_tgt][i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[is_tgt][1]; ++i)
      coordinates_y[is_tgt][i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x[is_tgt], coordinates_y[is_tgt], num_cells[is_tgt],
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

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
    }

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};
    if (!is_tgt)
      yac_basic_grid_add_coordinates_nocpy(
        grids[0], CELL, src_point_coordinates);

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
      {yac_interp_method_hcsbb_new(),
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

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
          ref_tgt_field =
            xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
          compute_field_data_XYZ(
            grid_data.vertex_coordinates, 100, ref_tgt_field);
          for (size_t i = 0; i < 100; ++i) tgt_field[i] = -1;
        } else {
          src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
          compute_field_data_XYZ(
            src_point_coordinates, grid_data.num_cells, src_field);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        if (is_tgt)
          for (size_t i = 0; i < 100; ++i)
            if (fabs(1.0 - ref_tgt_field[i] / tgt_field[i]) > 1e-2)
              PUT_ERR("wrong interpolation result");

        free(ref_tgt_field);
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
    // (source point location is at the cell centers)
    // (target grid contains cells that are not covered by any source cell)
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
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][10] = {{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},{-1.0}};
    double coordinates_y[2][10] = {{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0},{-1.0}};
    size_t const num_cells[2][2] = {{7,7}, {9,9}};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{9,9}}};
    int with_halo = 0;
    for (size_t i = 0; i < 10; ++i) {
      coordinates_x[1][i] = 4.0 + 4.0 * ((double)i / 9.0);
      coordinates_y[1][i] = 4.0 + 4.0 * ((double)i / 9.0);
    }
    for (size_t i = 0; i <= num_cells[is_tgt][0]; ++i)
      coordinates_x[is_tgt][i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[is_tgt][1]; ++i)
      coordinates_y[is_tgt][i] *= YAC_RAD;

    struct basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x[is_tgt], coordinates_y[is_tgt], num_cells[is_tgt],
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

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
    }

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};
    if (!is_tgt)
      yac_basic_grid_add_coordinates_nocpy(
        grids[0], CELL, src_point_coordinates);

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
      {yac_interp_method_hcsbb_new(),
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

        if (is_tgt) {
          int is_outside_src_grid[] = {
            0,0,0,0,0,0,1,1,1,1,
            0,0,0,0,0,0,1,1,1,1,
            0,0,0,0,0,0,1,1,1,1,
            0,0,0,0,0,0,1,1,1,1,
            0,0,0,0,0,0,1,1,1,1,
            0,0,0,0,0,0,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1};
          tgt_field = xmalloc(grid_data.num_vertices * sizeof(*tgt_field));
          ref_tgt_field =
            xmalloc(grid_data.num_vertices * sizeof(*ref_tgt_field));
          compute_field_data_XYZ(
            grid_data.vertex_coordinates, 100, ref_tgt_field);
          for (size_t i = 0; i < grid_data.num_vertices; ++i)
          if (is_outside_src_grid[grid_data.vertex_ids[i]])
              ref_tgt_field[i] = -1.0;
          for (size_t i = 0; i < 100; ++i) tgt_field[i] = -2.0;
        } else {
          src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
          compute_field_data_XYZ(
            src_point_coordinates, grid_data.num_cells, src_field);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        if (is_tgt)
          for (size_t i = 0; i < 100; ++i)
            if (fabs(1.0 - ref_tgt_field[i] / tgt_field[i]) > 1e-2)
              PUT_ERR("wrong interpolation result");

        free(ref_tgt_field);
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

#if defined __INTEL_COMPILER
#pragma intel optimization_level 0
#elif defined _CRAYC
#pragma _CRI noopt
#endif
static void compute_field_data_XYZ(
  double (*vertices)[3], size_t num_vertices, double * f_v) {

  for (size_t i = 0; i < num_vertices; ++i) {

    double x = vertices[i][0];
    double y = vertices[i][1];
    double z = vertices[i][2];

    // f = 1 + x^8 + e^(2y^3) + e^(2x^2) + 10xyz
    f_v[i] = 1.0 + pow(x, 8.0);
    f_v[i] += exp(2.0*y*y*y);
    f_v[i] += exp(2.0*x*x) + 10.0*x*y*z;

    // f = x^3 + 3.0 * x*y^2 + 5 * x*y*z
    // f_v[i] = x*x*x + 3.0*x*y*y + 5.0*x*y*z + 1.0;

    // f_v[i] = 1.0;
  }
}
