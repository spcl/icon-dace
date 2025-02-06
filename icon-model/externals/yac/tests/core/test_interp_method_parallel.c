// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>
#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "read_icon_grid.h"
#include "interp_method.h"
#include "interp_method_fixed.h"
#include "yac_mpi.h"

int main(int argc, char** argv) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  set_even_io_rank_list(MPI_COMM_WORLD);

  if (argc != 2) {
    PUT_ERR("ERROR: missing grid file directory");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  char * filenames[2];
  char * grid_filenames[] ={"icon_grid_0030_R02B03_G.nc", "icon_grid_0043_R02B04_G.nc"};
  for (int i = 0; i < 2; ++i)
    filenames[i] =
      strcat(
        strcpy(
          malloc(strlen(argv[1]) + strlen(grid_filenames[i]) + 2), argv[1]),
        grid_filenames[i]);

  struct yac_basic_grid_data grid_data[2];

  for (int i = 0; i < 2; ++i)
    grid_data[i] =
      yac_read_icon_basic_grid_data_parallel(
        filenames[i], MPI_COMM_WORLD);

  struct yac_basic_grid * grids[2] =
    {yac_basic_grid_new(filenames[0], grid_data[0]),
     yac_basic_grid_new(filenames[1], grid_data[1])};

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, filenames[0], filenames[1],
                        num_src_fields, src_fields, tgt_field);

  struct interp_method * method_stack[] =
    {yac_interp_method_fixed_new(-1.0), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  struct yac_interpolation * interpolation_src =
    yac_interp_weights_get_interpolation(
      weights, YAC_MAPPING_ON_SRC, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);
  struct yac_interpolation * interpolation_tgt =
    yac_interp_weights_get_interpolation(
      weights, YAC_MAPPING_ON_TGT, 1, YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0);

  yac_interpolation_delete(interpolation_src);
  yac_interpolation_delete(interpolation_tgt);

  yac_interp_weights_delete(weights);

  yac_interp_method_delete(method_stack);

  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  for (int i = 0; i < 2; ++i) {
    yac_basic_grid_delete(grids[i]);
    free(filenames[i]);
  }

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
