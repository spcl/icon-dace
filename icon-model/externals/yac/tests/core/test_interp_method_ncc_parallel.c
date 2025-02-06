// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "interp_method.h"
#include "interp_method_ncc.h"
#include "interp_method_fixed.h"
#include "dist_grid_utils.h"
#include "clipping.h"
#include "yac_mpi.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

enum yac_interp_weights_reorder_type reorder_types[] =
  {YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT};
size_t num_reorder_types = sizeof(reorder_types) / sizeof(reorder_types[0]);
static char const * grid_names[2] = {"src_grid", "tgt_grid"};

static double compute_ref_tgt_value(
  size_t tgt_corner, int * tgt_global_core_mask, int * tgt_global_field_mask,
  int * src_global_core_mask, int * src_global_field_maske, int with_src_mask,
  int with_tgt_mask, enum yac_interp_ncc_weight_type weight_type,
  int partial_coverage);

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size != 3) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  MPI_Comm split_comm;
  MPI_Comm_split(
    MPI_COMM_WORLD, comm_rank < 2, 0, &split_comm);

  int split_comm_rank, split_comm_size;
  MPI_Comm_rank(split_comm, &split_comm_rank);
  MPI_Comm_size(split_comm, &split_comm_size);

  {// test with various configuration options

    // The global source grid has 5x4 cells:
    //
    //    24--44--25--45--26--46--27--47--28--48--29
    //    |       |       |       |       |       |
    //    34  15  36  16  38  17  40  18  42  19  43
    //    |       |       |       |       |       |
    //    18--33--19--35--20--37--21--39--22--41--23
    //    |       |       |       |       |       |
    //    23  10  25  11  27  12  29  13  31  14  32
    //    |       |       |       |       |       |
    //    12--22--13--24--14--26--15--28--16--30--17
    //    |       |       |       |       |       |
    //    12  05  14  06  16  07  18  08  20  09  21
    //    |       |       |       |       |       |
    //    06--11--07--13--08--15--09--17--10--19--11
    //    |       |       |       |       |       |
    //    01  00  03  01  05  02  07  03  09  04  10
    //    |       |       |       |       |       |
    //    00--01--01--02--02--04--03--06--04--08--05
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][10] =
      {{0,1,2,3,4,5}, {0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75}};
    double coordinates_y[2][8] =
      {{-2,-1,0,1,2}, {-1.75,-1.25,-0.75,-0.25,0.25,0.75,1.25,1.75}};
    double cell_coordiantes_x[5 * 4] =
      {0.5,1.5,2.5,3.5,4.5,
       0.5,1.5,2.5,3.5,4.5,
       0.5,1.5,2.5,3.5,4.5,
       0.5,1.5,2.5,3.5,4.5};
    double cell_coordiantes_y[5 * 4] =
      {-1.5,-1.5,-1.5,-1.5,-1.5,
       -0.5,-0.5,-0.5,-0.5,-0.5,
       0.5,0.5,0.5,0.5,0.5,
       1.5,1.5,1.5,1.5,1.5};
    size_t const num_cells[2][2] = {{5,4}, {9,7}};
    size_t local_start[2][2][2] = {{{0,0},{3,0}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{3,4},{2,4}}, {{9,7}}};
    int global_core_mask[2][10*8] =
      {{1,0,0,0,0,
        1,1,0,0,0,
        1,1,1,1,0,
        1,1,1,1,0},
       {1,1,1,0,0,0,0,0,0,0,
        1,1,1,0,0,0,0,0,0,0,
        1,1,1,1,1,1,0,0,0,0,
        1,1,1,1,1,1,0,0,0,0,
        1,1,1,1,1,1,1,1,1,0,
        1,1,1,1,1,1,1,1,1,0,
        1,1,1,1,1,1,1,1,1,0,
        1,1,1,1,1,1,1,1,1,0}};
    int global_field_mask[2][10*8] =
      {{0,1,1,1,1,
        0,1,1,1,1,
        0,1,1,1,1,
        0,1,1,1,1},
       {1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        0,0,0,0,0,0,0,0,0,0}};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[is_tgt][0]; ++i)
      coordinates_x[is_tgt][i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[is_tgt][1]; ++i)
      coordinates_y[is_tgt][i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x[is_tgt], coordinates_y[is_tgt], num_cells[is_tgt],
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    if (is_tgt) {
      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        grid_data.core_vertex_mask[i] =
          global_core_mask[is_tgt][grid_data.vertex_ids[i]];
      free(grid_data.core_edge_mask);
      free(grid_data.core_cell_mask);
      grid_data.core_edge_mask = NULL;
      grid_data.core_cell_mask = NULL;
    } else {
      for (size_t i = 0; i < grid_data.num_cells; ++i)
        grid_data.core_cell_mask[i] =
          global_core_mask[is_tgt][grid_data.cell_ids[i]];
    }

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    int field_mask[10*8];
    double field_coords[10*8][3];
    if (is_tgt) {
      for (size_t i = 0; i < grid_data.num_vertices; ++i)
        field_mask[i] = global_field_mask[is_tgt][grid_data.vertex_ids[i]];
      yac_basic_grid_add_mask(
        grids[0], YAC_LOC_CORNER, field_mask, grid_data.num_vertices, NULL);
    } else {
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        field_mask[i] = global_field_mask[is_tgt][grid_data.cell_ids[i]];
        LLtoXYZ_deg(
          cell_coordiantes_x[grid_data.cell_ids[i]],
          cell_coordiantes_y[grid_data.cell_ids[i]], field_coords[i]);
      }
      yac_basic_grid_add_mask(
        grids[0], YAC_LOC_CELL, field_mask, grid_data.num_cells, NULL);
      yac_basic_grid_add_coordinates(
        grids[0], YAC_LOC_CELL, field_coords, grid_data.num_cells);
    }

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    for (int with_src_mask = 0; with_src_mask < 2; ++with_src_mask) {
      for (int with_tgt_mask = 0; with_tgt_mask < 2; ++with_tgt_mask) {

        struct yac_interp_field src_fields[] =
          {{.location = YAC_LOC_CELL,
            .coordinates_idx = 0,
            .masks_idx = (with_src_mask)?0:SIZE_MAX}};
        size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
        struct yac_interp_field tgt_field =
          {.location = YAC_LOC_CORNER,
           .coordinates_idx = SIZE_MAX,
           .masks_idx = (with_tgt_mask)?0:SIZE_MAX};

        struct yac_interp_grid * interp_grid =
          yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                              num_src_fields, src_fields, tgt_field);

        enum yac_interp_ncc_weight_type weight_types[] =
          {YAC_INTERP_NCC_AVG,
           YAC_INTERP_NCC_DIST};
        size_t const num_weight_types =
          sizeof(weight_types) / sizeof(weight_types[0]);

        for (size_t weight_types_idx = 0; weight_types_idx < num_weight_types;
             ++weight_types_idx) {
          for (int partial_coverage = 0; partial_coverage < 2;
               ++partial_coverage) {

            struct interp_method * method_stack[] =
              {yac_interp_method_ncc_new(
                 weight_types[weight_types_idx], partial_coverage),
               yac_interp_method_fixed_new(-1.0), NULL};

            struct yac_interp_weights * weights =
              yac_interp_method_do_search(method_stack, interp_grid);

            for (size_t i = 0; i < num_reorder_types; ++i) {

              struct yac_interpolation * interpolation =
                yac_interp_weights_get_interpolation(
                  weights, reorder_types[i], 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

              // check generated interpolation
              {
                double src_field_data[5*4];
                double * src_field = src_field_data;
                double ** src_fields = &src_field;
                double tgt_field_data[10*8] =
                  {-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,
                   -2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0};
                double * tgt_field = tgt_field_data;

                if (!is_tgt)
                  for (size_t i = 0; i < grid_data.num_cells; ++i)
                    src_field[i] = (double)(grid_data.cell_ids[i] + 1);

                yac_interpolation_execute(
                  interpolation, &src_fields, &tgt_field);

                if (is_tgt) {

                  for (size_t i = 0; i < grid_data.num_vertices; ++i) {

                    double ref_tgt_value =
                      compute_ref_tgt_value(
                        i, global_core_mask[1], global_field_mask[1],
                        global_core_mask[0], global_field_mask[0],
                        with_src_mask, with_tgt_mask,
                        weight_types[weight_types_idx], partial_coverage);

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

  {// small test

    // The global source grid has 5x4 cells:
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[2][4] =
      {{-2,-1,1,2}, {-1.0,0.0,1.0}};
    double coordinates_y[2][4] =
      {{-2,-1,1,2}, {-1.0,0.0,1.0}};
    double field_coordiantes_x[2][3*3] =
      {{-1.5,0.0,1.5,
        -1.5,0.0,1.5,
        -1.5,0.0,1.5},
       {-1.0,0.0,1.0,
        -1.0,0.0,1.0,
        -1.0,0.0,1.0}};
    double field_coordiantes_y[2][3*3] =
      {{-1.5,-1.5,-1.5,
        0.0,0.0,0.0,
        1.5,1.5,1.5},
       {-1.0,-1.0,-1.0,
        0.0,0.0,0.0,
        1.0,1.0,1.0}};
    int field_mask[3*3] = {0,1,1, 1,1,1, 1,1,0};
    size_t const num_cells[2][2] = {{3,3}, {2,2}};
    size_t local_start[2][2][2] = {{{0,0},{0,1}}, {{0,0}}};
    size_t local_count[2][2][2] = {{{3,2},{3,2}}, {{2,2}}};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[is_tgt][0]; ++i)
      coordinates_x[is_tgt][i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[is_tgt][1]; ++i)
      coordinates_y[is_tgt][i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x[is_tgt], coordinates_y[is_tgt], num_cells[is_tgt],
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    double field_coords[3*3][3];
    {
      const yac_int * ids = (is_tgt)?grid_data.vertex_ids:grid_data.cell_ids;
      size_t count = (is_tgt)?grid_data.num_vertices:grid_data.num_cells;
      enum yac_location field_location[2] = {YAC_LOC_CELL, YAC_LOC_CORNER};
      int src_field_mask[3*3];
      for (size_t i = 0; i < count; ++i) {
        LLtoXYZ_deg(
          field_coordiantes_x[is_tgt][ids[i]],
          field_coordiantes_y[is_tgt][ids[i]], field_coords[i]);
        if (!is_tgt) src_field_mask[i] = field_mask[ids[i]];
      }
      yac_basic_grid_add_coordinates(
        grids[0], field_location[is_tgt], field_coords, count);
      if (!is_tgt)
        yac_basic_grid_add_mask(
          grids[0], field_location[is_tgt], src_field_mask, count, NULL);
    }

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    for (int with_src_mask = 0; with_src_mask < 2; ++with_src_mask) {

      struct yac_interp_field src_fields[] =
        {{.location = YAC_LOC_CELL,
          .coordinates_idx = 0,
          .masks_idx = with_src_mask?0:SIZE_MAX}};
      size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
      struct yac_interp_field tgt_field =
        {.location = YAC_LOC_CORNER,
          .coordinates_idx = 0,
          .masks_idx = SIZE_MAX};

      struct yac_interp_grid * interp_grid =
        yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                            num_src_fields, src_fields, tgt_field);

      for (int partial_coverage = 0; partial_coverage < 2; ++partial_coverage) {

        struct interp_method * method_stack[] =
          {yac_interp_method_ncc_new(
              YAC_INTERP_NCC_DIST, partial_coverage),
            yac_interp_method_fixed_new(-1.0), NULL};

        struct yac_interp_weights * weights =
          yac_interp_method_do_search(method_stack, interp_grid);

        struct yac_interpolation * interpolation =
          yac_interp_weights_get_interpolation(
            weights, YAC_MAPPING_ON_SRC, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

        // check generated interpolation
        {
          double src_field_data[3*3];
          double * src_field = src_field_data;
          double ** src_fields = &src_field;
          double tgt_field_data[3*3] =
            {-2.0,-2.0,-2.0,
            -2.0,-2.0,-2.0,
            -2.0,-2.0,-2.0};
          double * tgt_field = tgt_field_data;

          if (!is_tgt)
            for (size_t i = 0; i < grid_data.num_cells; ++i)
              src_field[i] = (double)(grid_data.cell_ids[i] + 1);

          yac_interpolation_execute(
            interpolation, &src_fields, &tgt_field);

          if (is_tgt) {

            yac_int tgt_to_src[2][2][9][4] =
            {{{{0,1,3,4},{0,1,3,4},{1,2,4,5},
              {0,1,3,4},{4},{1,2,4,5},
              {3,4,6,7},{3,4,6,7},{4,5,7,8}},
              {{0,1,3,4},{1,2,4,5},{1,2,4,5},
              {3,4,6,7},{4},{4,5,7,8},
              {3,4,6,7},{4,5,7,8},{4,5,7,8}}},
            {{{1,3,4},{1,3,4},{1,2,4,5},
              {1,3,4},{4},{1,2,4,5},
              {3,4,6,7},{3,4,6,7},{4,5,7}},
              {{1,3,4},{1,2,4,5},{1,2,4,5},
              {3,4,6,7},{4},{4,5,7},
              {3,4,6,7},{4,5,7},{4,5,7}}}};
            int num_src_per_tgt[2][2][2][9] =
              {{{{4,4,4, 4,1,4, 4,4,4},{4,4,4, 4,1,4, 4,4,4}},
                {{0,0,4, 0,1,4, 4,4,0},{0,4,4, 4,1,0, 4,0,0}}},
               {{{4,4,4, 4,1,4, 4,4,4},{4,4,4, 4,1,4, 4,4,4}},
                {{3,3,4, 3,1,4, 4,4,3},{3,4,4, 4,1,3, 4,3,3}}}};

            for (size_t i = 0; i < grid_data.num_vertices; ++i) {

              yac_int tgt_id = grid_data.vertex_ids[i];
              double ref_tgt_value[2] = {0.0, 0.0};

              double tgt_coord[3];
              LLtoXYZ_deg(
                field_coordiantes_x[1][tgt_id], field_coordiantes_y[1][tgt_id],
                tgt_coord);

              for (int j = 0; j < 2; ++j) {

                double inv_distance_sum = 0.0;
                int num_src =
                  num_src_per_tgt[partial_coverage][with_src_mask][j][tgt_id];

                if (num_src > 0) {
                  for (int k = 0; k < num_src; ++k) {
                    yac_int src_id = tgt_to_src[with_src_mask][j][tgt_id][k];
                    double src_coord[3];
                    LLtoXYZ_deg(
                      field_coordiantes_x[0][src_id],
                      field_coordiantes_y[0][src_id], src_coord);
                    double angle = get_vector_angle(src_coord, tgt_coord);
                    if (angle < yac_angle_tol) {
                      ref_tgt_value[j] = (double)(src_id + 1);
                      inv_distance_sum = 1.0;
                      break;
                    }
                    double inv_distance =
                      1.0 / get_vector_angle(src_coord, tgt_coord);
                    inv_distance_sum += inv_distance;
                    ref_tgt_value[j] += (double)(src_id + 1) * inv_distance;
                  }

                  ref_tgt_value[j] /= inv_distance_sum;
                } else {
                  ref_tgt_value[j] = -1.0;
                }
              }

              if ((fabs(ref_tgt_value[0] - tgt_field[i]) > 1e-3) &&
                  (fabs(ref_tgt_value[1] - tgt_field[i]) > 1e-3))
                PUT_ERR("wrong interpolation result");
            }
          }
        }

        yac_interpolation_delete(interpolation);
        yac_interp_weights_delete(weights);
        yac_interp_method_delete(method_stack);
      } // partial_coverage
      yac_interp_grid_delete(interp_grid);
    } // with_src_mask
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  MPI_Comm_free(&split_comm);

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static double compute_ref_tgt_value(
  size_t tgt_corner, int * tgt_global_core_mask, int * tgt_global_field_mask,
  int * src_global_core_mask, int * src_global_field_mask, int with_src_mask,
  int with_tgt_mask, enum yac_interp_ncc_weight_type weight_type,
  int partial_coverage) {

  // if target corner is masked
  if (!tgt_global_core_mask[tgt_corner] ||
      (with_tgt_mask && !tgt_global_field_mask[tgt_corner])) return -2.0;

  // determine matching source cell
  int source_grid_x_idx = ((int)tgt_corner % 10) / 2;
  int source_grid_y_idx = ((int)tgt_corner / 10) / 2;
  int source_cell = source_grid_x_idx + source_grid_y_idx * 5;

  // check source cell core mask
  if (!src_global_core_mask[source_cell]) return -1.0;

  // determine closest source cell corner
  int is_right_corner = (int)tgt_corner & 1;
  int is_upper_corner = ((int)tgt_corner / 10) & 1;
  int source_corner =
    source_grid_x_idx + is_right_corner +
    6 * (source_grid_y_idx + is_upper_corner);

  // determine source cells surrounding the source vertex
  int source_cells[4];
  source_cells[0] = (5 * (source_corner / 6) + source_corner % 6) - 5 - 1;
  source_cells[1] = (5 * (source_corner / 6) + source_corner % 6) - 5;
  source_cells[2] = (5 * (source_corner / 6) + source_corner % 6) - 1;
  source_cells[3] = (5 * (source_corner / 6) + source_corner % 6);

  // special handling for grid edge corners
  if (source_corner <= 5) {
    source_cells[0] = -1;
    source_cells[1] = -1;
  }
  if (source_corner >= 24) {
    source_cells[2] = -1;
    source_cells[3] = -1;
  }
  if ((source_corner % 6) == 0) {
    source_cells[0] = -1;
    source_cells[2] = -1;
  }
  if (((source_corner + 1) % 6) == 0) {
    source_cells[3] = -1;
    source_cells[1] = -1;
  }

  // apply core mask to source cells
  for (int i = 0; i < 4; ++i)
    if ((source_cells[i] >= 0) &&
        (!src_global_core_mask[source_cells[i]]))
      source_cells[i] = -1;

  // apply field mask to source cells
  if (with_src_mask) {
    for (int i = 0; i < 4; ++i) {
      if ((source_cells[i] >= 0) &&
          (!src_global_field_mask[source_cells[i]])) {
        if (!partial_coverage) return -1.0;
        source_cells[i] = -1;
      }
    }
  }

  double const ref_distances[2][2] =
    {{sqrt(0.25*0.25 + 0.25*0.25),
      sqrt(0.75*0.75 + 0.25*0.25)},
     {sqrt(0.25*0.25 + 0.75*0.75),
      sqrt(0.75*0.75 + 0.75*0.75)}};

  double distances[4];
  int num_src_cells = 0;
  for (int i = 0, k = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j, ++k) {
      if (source_cells[k] >= 0) {
        distances[num_src_cells] =
          ref_distances[i^is_upper_corner^1][j^is_right_corner^1];
        source_cells[num_src_cells] = source_cells[k];
        ++num_src_cells;
      }
    }
  }

  if (num_src_cells == 0) return -1.0;

  double tgt_value = 0.0;
  switch(weight_type) {
    default:
    case (YAC_INTERP_NCC_AVG): {
      double weight = 1.0 / (double)num_src_cells;
      for (int i = 0; i < num_src_cells; ++i)
        tgt_value += (double)(source_cells[i] + 1) * weight;
      break;
    }
    case (YAC_INTERP_NCC_DIST): {
      double inv_distance_sum = 0.0;
      for (int i = 0; i < num_src_cells; ++i)
        inv_distance_sum += 1.0 / distances[i];
      for (int i = 0; i < num_src_cells; ++i)
        tgt_value +=
          (double)(source_cells[i] + 1) / (distances[i] * inv_distance_sum);
      break;
    }
  };
  return tgt_value;
}
