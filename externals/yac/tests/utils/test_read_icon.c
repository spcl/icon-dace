// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "read_icon_grid.h"
#include "grid2vtk.h"
#include "io_utils.h"

int main(int argc, char** argv) {

   if (argc != 2) {
      PUT_ERR("ERROR: missing grid file directory");
      return TEST_EXIT_CODE;
   }

   char * grid_filename =
      strcat(
         strcpy(
           malloc(strlen(argv[1]) + 32), argv[1]), "icon_grid_0030_R02B03_G.nc");

   if (!yac_file_exists(grid_filename))
      return EXIT_SKIP_TEST;

   struct yac_basic_grid_data icon_grid =
      yac_read_icon_basic_grid_data(grid_filename);

   free(grid_filename);

   if (icon_grid.num_cells != 5120)
      PUT_ERR("wrong number of grid cells");

// #define WRITE_VTK_GRID_FILE
#ifdef WRITE_VTK_GRID_FILE
   yac_write_basic_grid_data_to_file(&icon_grid, "icon");
#endif // WRITE_VTK_GRID_FILE

   yac_basic_grid_data_free(icon_grid);

   return TEST_EXIT_CODE;
}

