// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tests.h"
#include "read_cube_csv_grid.h"
#include "grid2vtk.h"

int main(int argc, char** argv) {

  if (argc != 2) {
    PUT_ERR("ERROR: missing grid file directory");
    return TEST_EXIT_CODE;
  }

  char * grid_filename =
    strcat(
      strcpy(
        malloc(strlen(argv[1]) + 32), argv[1]), "cube_10x10.csv");

  struct yac_basic_grid_data cube_grid =
    yac_read_cube_csv_grid(grid_filename);

  free(grid_filename);

  if (cube_grid.num_cells != 600)
    PUT_ERR("ERROR: wrong number of cells\n");

  if (cube_grid.num_vertices != 602)
    PUT_ERR("ERROR: wrong number of grid vertices\n")

// #define WRITE_VTK_GRID_FILE
#ifdef WRITE_VTK_GRID_FILE
  yac_write_basic_grid_data_to_file(&cube_grid, "cube");
#endif // WRITE_VTK_GRID_FILE

  yac_basic_grid_data_free(cube_grid);

  return TEST_EXIT_CODE;
}
