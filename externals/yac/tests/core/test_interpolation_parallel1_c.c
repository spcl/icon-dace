// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

#include "tests.h"
#include "test_common.h"
#include "dist_grid_utils.h"
#include "interp_method_avg.h"
#include "interp_method_fixed.h"

#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

#define MAX_COLLECTION_SIZE (10)

char const grid_name_src[] = "src_grid";
char const grid_name_tgt[] = "tgt_grid";

static void target_main(
  MPI_Comm target_comm, enum yac_interp_weights_reorder_type reorder_type);

static void source_main(
  MPI_Comm source_comm, enum yac_interp_weights_reorder_type reorder_type);

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

  if (comm_size != 4) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  // split processes into source and target

  int tgt_flag = comm_rank < 2;

  MPI_Comm split_comm;
  MPI_Comm_split(MPI_COMM_WORLD, tgt_flag, 0, &split_comm);

  if (tgt_flag) target_main(split_comm, reorder_type);
  else          source_main(split_comm, reorder_type);

  MPI_Comm_free(&split_comm);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

/*
 * The source grid is distributed among 2 processes
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
static void source_main(
  MPI_Comm source_comm, enum yac_interp_weights_reorder_type reorder_type) {

  int my_source_rank;
  MPI_Comm_rank(source_comm, &my_source_rank);

  double coordinates_x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  double coordinates_y[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  size_t const num_cells[2] = {5,4};
  size_t local_start[2][2] = {{0,0},{3,0}};
  size_t local_count[2][2] = {{3,4},{2,4}};
  int with_halo = 1;
  for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
  for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

  struct yac_basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg2d(
      coordinates_x, coordinates_y, num_cells,
      local_start[my_source_rank], local_count[my_source_rank], with_halo);
  struct yac_basic_grid * src_grid =
    yac_basic_grid_new(grid_name_src, grid_data);
  struct yac_basic_grid * tgt_grid =
    yac_basic_grid_empty_new(grid_name_tgt);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, grid_name_src, grid_name_tgt,
                        num_src_fields, src_fields, tgt_field);

  struct interp_method * method_stack[] =
    {yac_interp_method_avg_new(YAC_INTERP_AVG_ARITHMETIC, 0),
     yac_interp_method_fixed_new(1337),
     NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

  // -------------------
  // set up source data
  // -------------------

  // src_data dimensions [collection_idx]
  //                     [pointset_idx]
  //                     [local_idx]
  double *** src_data = xmalloc(MAX_COLLECTION_SIZE * sizeof(*src_data));
  for (size_t collection_idx = 0; collection_idx < MAX_COLLECTION_SIZE;
       ++collection_idx) {
    src_data[collection_idx] = xmalloc(1 * sizeof(**src_data));
    src_data[collection_idx][0] =
      xmalloc(grid_data.num_vertices * sizeof(***src_data));
    for (size_t i = 0; i < grid_data.num_vertices; ++i)
      src_data[collection_idx][0][i] =
        (grid_data.core_vertex_mask[i])?
          ((double)(grid_data.vertex_ids[i]) + (double)(collection_idx * 30)):
          (-1.0);
  }

  for (size_t collection_size = 1; collection_size <= MAX_COLLECTION_SIZE;
       ++collection_size) {

    //---------------------------
    // set up field interpolation
    //---------------------------

    struct yac_interpolation * interpolation =
      yac_interp_weights_get_interpolation(
        weights, reorder_type, collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    struct yac_interpolation * interpolation_copy =
      yac_interpolation_copy(interpolation);

    //---------------------
    // do the interpolation
    //---------------------

    yac_interpolation_execute_put(interpolation, src_data);
    yac_interpolation_execute_put(interpolation_copy, src_data);

    yac_interpolation_delete(interpolation_copy);
    yac_interpolation_delete(interpolation);
  }

  //--------
  // cleanup
  //--------


  for (size_t collection_idx = 0; collection_idx < MAX_COLLECTION_SIZE;
       ++collection_idx) {
    free(src_data[collection_idx][0]);
    free(src_data[collection_idx]);
  }
  free(src_data);

  yac_basic_grid_delete(tgt_grid);
  yac_basic_grid_delete(src_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);
}

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
static void target_main(
  MPI_Comm target_comm, enum yac_interp_weights_reorder_type reorder_type) {

  int my_target_rank;
  MPI_Comm_rank(target_comm, &my_target_rank);

  double coordinates_x[] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5};
  double coordinates_y[] = {0.5,1.5,2.5,3.5};
  size_t const num_cells[2] = {6,3};
  size_t local_start[2][2] = {{0,0},{3,0}};
  size_t local_count[2][2] = {{3,3},{3,3}};
  int with_halo = 0;
  for (size_t i = 0; i <= num_cells[0]; ++i) coordinates_x[i] *= YAC_RAD;
  for (size_t i = 0; i <= num_cells[1]; ++i) coordinates_y[i] *= YAC_RAD;

  struct yac_basic_grid_data grid_data =
    yac_generate_basic_grid_data_reg2d(
      coordinates_x, coordinates_y, num_cells,
      local_start[my_target_rank], local_count[my_target_rank], with_halo);
  struct yac_basic_grid * tgt_grid =
    yac_basic_grid_new(grid_name_tgt, grid_data);
  struct yac_basic_grid * src_grid =
    yac_basic_grid_empty_new(grid_name_src);

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(tgt_grid, src_grid, MPI_COMM_WORLD);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CORNER, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, grid_name_src, grid_name_tgt,
                        num_src_fields, src_fields, tgt_field);

  struct interp_method * method_stack[] =
    {yac_interp_method_avg_new(YAC_INTERP_AVG_ARITHMETIC, 0),
     yac_interp_method_fixed_new(1337),
     NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);

   //---------------------
   // do the interpolation
   //---------------------

   double target_field[MAX_COLLECTION_SIZE][16];
   double * target_data[MAX_COLLECTION_SIZE];
   for (unsigned i = 0; i < MAX_COLLECTION_SIZE; ++i)
      target_data[i] = target_field[i];

  for (unsigned collection_size = 1; collection_size <= MAX_COLLECTION_SIZE;
       ++collection_size) {

    //---------------------------
    // set up field interpolation
    //---------------------------

    struct yac_interpolation * interpolation[2];
    interpolation[0] =
      yac_interp_weights_get_interpolation(
        weights, reorder_type, collection_size, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
    interpolation[1] = yac_interpolation_copy(interpolation[0]);

    for (size_t k = 0; k < 2; ++k) {

      for (unsigned i = 0; i < collection_size; ++i)
        for (unsigned j = 0; j < 16; ++j)
          target_field[i][j] = -1;

      yac_interpolation_execute_get(interpolation[k], target_data);

      //----------------------------
      // check interpolation results
      //----------------------------

      double ref_target_data[2][16] = {{3.5,4.5,5.5,6.5,
                                        9.5,10.5,11.5,12.5,
                                        15.5,16.5,17.5,18.5,
                                        21.5,22.5,23.5,24.5},
                                      {6.5,7.5,1337,1337,
                                        12.5,13.5,1337,1337,
                                        18.5,19.5,1337,1337,
                                        24.5,25.5,1337,1337}};

      for (unsigned i = 0; i < collection_size; ++i) {
        for (unsigned j = 0; j < 16; ++j) {
          if ((double_are_equal(ref_target_data[my_target_rank][j], -1.0)) ||
              (double_are_equal(ref_target_data[my_target_rank][j], 1337.0))) {
            if (double_are_unequal(target_data[i][j],
                                  ref_target_data[my_target_rank][j]))
              PUT_ERR("error in interpolated data on target side\n");
          } else {
            if (double_are_unequal(
                  target_data[i][j],
                  ref_target_data[my_target_rank][j] + (double)(i * 30)))
              PUT_ERR("error in interpolated data on target side\n");
          }
        }
      }
    }

    for (int i = 0; i < 2; ++i)
      yac_interpolation_delete(interpolation[i]);
  }

  //---------
  // clean up
  //---------

  yac_basic_grid_delete(src_grid);
  yac_basic_grid_delete(tgt_grid);
  yac_interp_weights_delete(weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);
}

