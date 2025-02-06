// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <mpi.h>

#include "tests.h"
#include "test_common.h"
#include "weight_file_common.h"
#include "dist_grid_utils.h"
#include "dist_grid.h"
#include "geometry.h"
#include "interp_weights.h"
#include "interp_method_avg.h"
#include "interp_method_file.h"
#include "interp_method_fixed.h"
#include "interp_method_nnn.h"

#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

int main(int argc, char *argv[]) {

  if (argc != 2) {
    PUT_ERR("wrong number of arguments\n");
    return TEST_EXIT_CODE;
  }

  enum yac_interp_weights_reorder_type reorder_type =
   (strcmp(argv[1], "src") == 0)?YAC_MAPPING_ON_SRC:YAC_MAPPING_ON_TGT;

  if ((reorder_type != YAC_MAPPING_ON_SRC) && strcmp(argv[1], "tgt")) {
    PUT_ERR("invalid argument (has to be either \"src\" or \"tgt\")\n");
    return TEST_EXIT_CODE;
  }

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

  unsigned is_source = comm_rank > 0;
  unsigned is_target = comm_rank < 2;

  /*
   * The source grid is distributed among 2 processes.
   *
   * The global source grid has 5x4 cells:
   *
   *    24--44--25--45--26--46--27--47--28--48--29
   *    |       |       |       |       |       |
   *    34  15  36  16  38  17  40  18  42  19  43
   *    |       |       |       |       |       |
   *    18--33--19--35--20--37--21--39--22--41--23
   *    |       |       |       |       |       |
   *    23  10  25  11  27  12  29  13  31  14  32
   *    |       |       |       |       |       |
   *    12--22--13--24--14--26--15--28--16--30--17
   *    |       |       |       |       |       |
   *    12  05  14  06  16  07  18  08  20  09  21
   *    |       |       |       |       |       |
   *    06--11--07--13--08--15--09--17--10--19--11
   *    |       |       |       |       |       |
   *    01  00  03  01  05  02  07  03  09  04  10
   *    |       |       |       |       |       |
   *    00--01--01--02--02--04--03--06--04--08--05
   */
  /*
   * The target grid is distributed among 2 processes
   *
   * The global target grid has 6x3 cells:
   *
   *    21--39--22--40--23--41--24--42--25--43--26--44--27
   *    |       |       |       |       |       |       |
   *    27  12  29  13  31  14  33  15  35  16  37  17  38
   *    |       |       |       |       |       |       |
   *    14--26--15--28--16--30--17--32--18--34--19--36--20
   *    |       |       |       |       |       |       |
   *    14  06  16  07  18  08  20  09  22  10  24  11  25
   *    |       |       |       |       |       |       |
   *    07--13--08--15--09--17--10--19--11--21--12--23--13
   *    |       |       |       |       |       |       |
   *    01  00  03  01  05  02  07  03  09  04  11  05  12
   *    |       |       |       |       |       |       |
   *    00--00--01--02--02--04--03--06--04--08--05--10--06
   */

  char const source_grid_name[] = "source_grid";
  char const target_grid_name[] = "target_grid";

  double src_coordinates_x[] = {0,1,2,3,4,5};
  double src_coordinates_y[] = {0,1,2,3,4};
  size_t src_global_num_cells[] = {5,4};
  int src_with_halo = 1;
  size_t src_local_start[2][2] = {{0,0}, {0,2}};
  size_t src_local_count[2] = {5,2};

  for (size_t i = 0; i <= src_global_num_cells[0]; ++i)
    src_coordinates_x[i] *= YAC_RAD;
  for (size_t i = 0; i <= src_global_num_cells[1]; ++i)
    src_coordinates_y[i] *= YAC_RAD;

  struct yac_basic_grid * src_grid =
    (is_source)?
      yac_basic_grid_new(
        source_grid_name,
        yac_generate_basic_grid_data_reg2d(
          src_coordinates_x, src_coordinates_y, src_global_num_cells,
          src_local_start[comm_rank-1], src_local_count, src_with_halo)):
      yac_basic_grid_empty_new(source_grid_name);

  double tgt_coordinates_x[] = {0.7,1.7,2.7,3.7,4.7,5.7,6.7};
  double tgt_coordinates_y[] = {0.7,1.7,2.7,3.7};
  size_t tgt_global_num_cells[] = {6,3};
  int tgt_with_halo = 1;
  size_t tgt_local_start[2][2] = {{3,0}, {0,0}};
  size_t tgt_local_count[2] = {3,3};

  for (size_t i = 0; i <= tgt_global_num_cells[0]; ++i)
    tgt_coordinates_x[i] *= YAC_RAD;
  for (size_t i = 0; i <= tgt_global_num_cells[1]; ++i)
    tgt_coordinates_y[i] *= YAC_RAD;

  struct yac_basic_grid * tgt_grid =
    (is_target)?
      yac_basic_grid_new(
        target_grid_name,
        yac_generate_basic_grid_data_reg2d(
          tgt_coordinates_x, tgt_coordinates_y, tgt_global_num_cells,
          tgt_local_start[comm_rank], tgt_local_count, tgt_with_halo)):
      yac_basic_grid_empty_new(target_grid_name);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

  { // test interpolation

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, source_grid_name, target_grid_name,
                          num_src_fields, src_fields, tgt_field);

    // generate a weight file in which all weights are one
    char const * sum_weight_file_name =
      "test_interpolation_parallel5_weight_file.nc";
    if (comm_rank == 0) {
      int src_indices[4*5*4] =
        {0,1,6,7,     1,2,7,8,     2,3,8,9,     3,4,9,10,    4,5,10,11,
         6,7,12,13,   7,8,13,14,   8,9,14,15,   9,10,15,16,  10,11,16,17,
         12,13,18,19, 13,14,19,20, 14,15,20,21, 15,16,21,22, 16,17,22,23,
         18,19,24,25, 19,20,25,26, 20,21,26,27, 21,22,27,28, 22,23,28,29};
      int tgt_indices[4*5*4] =
        {0,0,0,0,     1,1,1,1,     2,2,2,2,     3,3,3,3,     4,4,4,4,
         7,7,7,7,     8,8,8,8,     9,9,9,9,     10,10,10,10, 11,11,11,11,
         14,14,14,14, 15,15,15,15, 16,16,16,16, 17,17,17,17, 18,18,18,18,
         21,21,21,21, 22,22,22,22, 23,23,23,23, 24,24,24,24, 25,25,25,25};
      double weights[4*5*4];
      for (size_t i = 0; i < 4*5*4; ++i) weights[i] = 1.0;
      size_t num_links = 4*5*4;
      enum yac_location src_locations[1] = {YAC_LOC_CORNER};
      enum yac_location tgt_location = YAC_LOC_CORNER;
      unsigned num_src_fields = 1;
      int num_links_per_field[1] = {num_links};
      int * tgt_id_fixed = NULL;
      size_t num_fixed_tgt = 0;
      double * fixed_values = NULL;
      int * num_tgt_per_fixed_value = NULL;
      size_t num_fixed_values = 0;

      write_weight_file(
        sum_weight_file_name, src_indices, tgt_indices, weights, num_links,
        src_locations, num_src_fields, num_links_per_field, tgt_id_fixed,
        num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
        num_fixed_values, tgt_location,
        source_grid_name, target_grid_name);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    size_t const num_method_stacks = 3;
    struct interp_method * method_stacks[3][3] =
      {{yac_interp_method_avg_new(YAC_INTERP_AVG_ARITHMETIC, 0),
        yac_interp_method_fixed_new(1337),
        NULL},
       {yac_interp_method_nnn_new(
          (struct yac_nnn_config){
            .type = YAC_INTERP_NNN_AVG, .n = 1,
            .max_search_distance =
              YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT}),
        yac_interp_method_fixed_new(1337),
        NULL},
       {yac_interp_method_file_new(sum_weight_file_name),
        yac_interp_method_fixed_new(1337),
        NULL}};

    for (size_t method_stack_index = 0;
         method_stack_index < num_method_stacks; ++method_stack_index) {

      struct yac_interp_weights * weights =
        yac_interp_method_do_search(
          method_stacks[method_stack_index], interp_grid);

      yac_interp_method_delete(method_stacks[method_stack_index]);

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_type, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

      char const weight_file_out[] =
        "test_interpolation_parallel5_weight_file_out.nc";
      yac_interp_weights_write_to_file(
        weights, weight_file_out, source_grid_name, target_grid_name, 0, 0);

      //---------------------
      // do the interpolation
      //---------------------

      struct yac_basic_grid_data * src_grid_data =
        is_source?yac_basic_grid_get_data(src_grid):NULL;
      double * source_data_field =
        is_source?
          xmalloc(src_grid_data->num_vertices * sizeof(*source_data_field)):NULL;
      double * source_data_pointset[1] = {source_data_field};  // num_pointset == 1
      double ** source_data[1] = {source_data_pointset}; // collection_size == 1

      struct yac_basic_grid_data * tgt_grid_data =
        is_target?yac_basic_grid_get_data(tgt_grid):NULL;
      double * target_data_field =
        is_target?
          xmalloc(tgt_grid_data->num_vertices * sizeof(*target_data_field)):NULL;
      double * target_data[1] = {target_data_field}; // collection_size == 1

      if (is_source) {

        // source_data dimensions [collection_idx]
        //                        [pointset_idx]
        //                        [local_idx]
        for (size_t i = 0; i < src_grid_data->num_vertices; ++i)
          source_data_field[i] =
            (src_grid_data->core_vertex_mask[i])?
              ((double)(src_grid_data->vertex_ids[i])):(-1.0);
      }

      for (int use_put_get = 0; use_put_get < 2; ++use_put_get) {

        if (is_target) {

          // target_data dimensions [collection_idx]
          //                        [local_idx]
          for (size_t i = 0; i < tgt_grid_data->num_vertices; ++i)
            target_data_field[i] = -1.0;
        }

        if (use_put_get) {
          if (is_source)
            yac_interpolation_execute_put(interpolation, source_data);
          if (is_target)
            yac_interpolation_execute_get(interpolation, target_data);
        } else {
          yac_interpolation_execute(interpolation, source_data, target_data);
        }

        if (is_target) {

          //----------------------------
          // check interpolation results
          //----------------------------

          double ref_target_data[3][4*7] =
            {{3.5,4.5,5.5,6.5,7.5,1337,1337,
              9.5,10.5,11.5,12.5,13.5,1337,1337,
              15.5,16.5,17.5,18.5,19.5,1337,1337,
              21.5,22.5,23.5,24.5,25.5,1337,1337},
             {7,8,9,10,11,11,11,
              13,14,15,16,17,17,17,
              19,20,21,22,23,23,23,
              25,26,27,28,29,29,29},
             {14,18,22,26,30,1337,1337,
              38,42,46,50,54,1337,1337,
              62,66,70,74,78,1337,1337,
              86,90,94,98,102,1337,1337}};

          for (size_t i = 0; i < tgt_grid_data->num_vertices; ++i) {
            if (tgt_grid_data->core_vertex_mask[i]) {
              if (double_are_unequal(
                    target_data[0][i],
                    ref_target_data
                      [method_stack_index][tgt_grid_data->vertex_ids[i]]))
                PUT_ERR("error in interpolated data on target side\n");
            } else {
              if (target_data[0][i] != -1.0)
                PUT_ERR("error in interpolated data on target side\n");
            }
          }
        }
      }

      if (is_target) free(target_data_field);
      if (is_source) free(source_data_field);

      //--------
      // cleanup
      //--------

      if (is_source || is_target) yac_interpolation_delete(interpolation);
      yac_interp_weights_delete(weights);

      if (comm_rank == 0) unlink(weight_file_out);
    }

    if (comm_rank == 0) unlink(sum_weight_file_name);

    yac_interp_grid_delete(interp_grid);
  }

  //--------
  // cleanup
  //--------

  yac_dist_grid_pair_delete(grid_pair);
  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

