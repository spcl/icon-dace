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

/** \example test_read_scrip.c
 * This contains examples for read_scrip_grid.
 */

/**
 * reads in an grid in SCRIP format
 * @param[in] grid_filename    name of the SCRIP grid netcdf file
 * @param[in] mask_filename    name of the SCRIP mask netcdf file
 * @param[in] grid_name        name of the grid in the file
 * @param[in] valid_mask_value value that marks cells as valid
 * @returns %basic_grid_data structure that contains the grid
 */
struct basic_grid_data read_scrip_grid(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value);

/**
 * reads in an grid in SCRIP format
 * @param[in]  grid_filename         name of the SCRIP grid netcdf file
 * @param[in]  mask_filename         name of the SCRIP mask netcdf file
 * @param[in]  grid_name             name of the grid in the file
 * @param[in]  valid_mask_value      value that marks cells as valid
 * @param[out] num_vertices          number of vertices in the grid
 * @param[out] num_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] x_vertices            longitude coordinates of the vertices
 * @param[out] y_vertices            latitude coordinates of the vertices
 * @param[out] x_cells               longitude coordinates of cell points
 * @param[out] y_cells               latitude coordinates of cell points
 * @param[out] cell_to_vertex        vertices indices per cell
 * @param[out] cell_core_mask        cell core mask
 */
void read_scrip_grid_information(
  char const * grid_filename, char const * mask_filename,
  char const * grid_name, int valid_mask_value,
  size_t * num_vertices, size_t * num_cells, int ** num_vertices_per_cell,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells,
  int ** cell_to_vertex, int ** cell_core_mask);
