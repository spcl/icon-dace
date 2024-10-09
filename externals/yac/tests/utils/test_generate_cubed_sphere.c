// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <unistd.h>

#include "tests.h"
#include "generate_cubed_sphere.h"
#include "grid2vtk.h"
#include "read_icon_grid.h"

int main(void) {

  {
    unsigned n = 20;

    struct yac_basic_grid_data cube_grid =
      yac_generate_cubed_sphere_grid(n);

    if (cube_grid.num_cells != n * n * 6)
      PUT_ERR("ERROR: wrong number of cells\n");

    if (cube_grid.num_vertices != n * n * 6 + 2)
      PUT_ERR("ERROR: wrong number of grid vertices\n")

// #define WRITE_VTK_GRID_FILE
#ifdef WRITE_VTK_GRID_FILE
    yac_write_basic_grid_data_to_file(&cube_grid, "cubed_sphere");
#endif // WRITE_VTK_GRID_FILE

    yac_basic_grid_data_free(cube_grid);
  }

  {
    unsigned n = 32;
    unsigned ref_num_core_cells = n * n * 6;
    int sizes[] = {1,3,6,11};
    for (size_t i = 0; i < (sizeof(sizes)/sizeof(sizes[0])); ++i) {

      unsigned num_core_cells = 0;

      for (int rank = 0, size = sizes[i]; rank < size; ++rank) {

        unsigned nbr_vertices;
        unsigned nbr_cells;
        unsigned * num_vertices_per_cell;
        unsigned * cell_to_vertex;
        double * x_vertices;
        double * y_vertices;
        double * x_cells;
        double * y_cells;
        int * global_cell_id;
        int * cell_core_mask;
        int * global_corner_id;
        int * corner_core_mask;

        yac_generate_part_cube_grid_information(
          n, &nbr_vertices, &nbr_cells, &num_vertices_per_cell, &cell_to_vertex,
          &x_vertices, &y_vertices, &x_cells, &y_cells,
          &global_cell_id, &cell_core_mask, &global_corner_id, &corner_core_mask,
          rank, size);

        for (unsigned j = 0; j < nbr_cells; ++j)
          if (cell_core_mask[j]) ++num_core_cells;

        free(num_vertices_per_cell);
        free(cell_to_vertex);
        free(x_vertices);
        free(y_vertices);
        free(x_cells);
        free(y_cells);
        free(global_cell_id);
        free(cell_core_mask);
        free(global_corner_id);
        free(corner_core_mask);
      }
      if (ref_num_core_cells != num_core_cells)
        PUT_ERR("wrong number of cells");
    }
  }

  {
    char const * grid_filename = "test_generate_cubed_sphere.nc";

    unsigned n = 20;

    yac_write_cubed_sphere_grid(n, grid_filename);

    struct yac_basic_grid_data cube_grid =
      yac_read_icon_basic_grid_data(grid_filename);

    if (cube_grid.num_cells != n * n * 6)
      PUT_ERR("ERROR: wrong number of cells\n");

    if (cube_grid.num_vertices != n * n * 6 + 2)
      PUT_ERR("ERROR: wrong number of grid vertices\n")

// #define WRITE_VTK_GRID_FILE
#ifdef WRITE_VTK_GRID_FILE
    yac_write_basic_grid_data_to_file(&cube_grid, "cube");
#endif // WRITE_VTK_GRID_FILE

    yac_basic_grid_data_free(cube_grid);

    unlink(grid_filename);
  }

  return TEST_EXIT_CODE;
}

