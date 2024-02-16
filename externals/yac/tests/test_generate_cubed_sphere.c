/**
 * @file test_generate_cubed_sphere.c
 *
 * @copyright Copyright  (C)  2015 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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

#include "tests.h"
#include "generate_cubed_sphere.h"
#include "grid2vtk.h"

int main(void) {

  {
    unsigned n = 20;

    struct basic_grid_data cube_grid = generate_cubed_sphere_grid(n);

    if (cube_grid.num_cells != n * n * 6)
      PUT_ERR("ERROR: wrong number of cells\n");

    if (cube_grid.num_vertices != n * n * 6 + 2)
      PUT_ERR("ERROR: wrong number of grid vertices\n")

#define WRITE_VTK_GRID_FILE
#ifdef WRITE_VTK_GRID_FILE
    yac_write_basic_grid_to_file(&cube_grid, "cubed_sphere");
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

        generate_part_cube_grid_information(
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

  return TEST_EXIT_CODE;
}

