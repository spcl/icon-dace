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
#include "grid.h"

/** \example test_read_mpiom.c
 * This contains examples for read_mpiom_grid.
 */

/**
 * reads in an mpiom grid netcdf file and generates a struct %grid from it
 * @param[in] filename name of the mpiom grid netcdf file
 * @returns %grid structure that contains the mpiom grid
 */
struct basic_grid_data read_mpiom_grid(char * filename);

/**
 * reads in an mpiom grid netcdf file and return the grid information in
 * a format that is supported by the YAC user interface.
 * @param[in]  filename              name of the mpiom grid netcdf file
 * @param[out] num_vertices          number of vertices in the grid
 * @param[out] num_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] cellmask              mask value for cells
 *
 */
void read_mpiom_grid_information(const char * filename, int * num_vertices,
                                 int * num_cells, int ** num_vertices_per_cell,
                                 int ** cell_to_vertex, double ** x_vertices,
                                 double ** y_vertices, double ** x_cells,
                                 double ** y_cells, int ** cellmask);

/**
 * reads in an partition mpiom grid netcdf file and return the grid information in
 * a format that is supported by the YAC user interface.
 * @param[in]  filename              name of the mpiom grid netcdf file
 * @param[out] num_vertices          number of vertices in the grid
 * @param[out] num_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] global_cell_id        global cell IDs
 * @param[out] cell_mask             mask value for cells
 * @param[out] cell_core_mask        cell core mask
 * @param[out] global_corner_id      global corner IDs
 * @param[out] corner_core_mask      corner core mask
 * @param[out] rank                  local MPI rank
 * @param[out] size                  number of MPI ranks
 *
 */
void read_part_mpiom_grid_information(
  const char * filename, int * num_vertices, int * num_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex, double ** x_vertices,
  double ** y_vertices, double ** x_cells,
  double ** y_cells, int ** global_cell_id,
  int ** cell_mask, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size);

