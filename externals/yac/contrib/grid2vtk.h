/**
 * @file grid2vtk.h
 *
 * @copyright Copyright  (C)  2018 DKRZ, MPI-M
 *
 * @version 1.0
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

#ifndef GRID2VTK_H
#define GRID2VTK_H

#include "grid.h"

/** \example test_grid2vtk.c
 * This contains some examples on how to use the
 * \ref yac_write_basic_grid_to_file and
 * \ref yac_write_grid_cells_to_file routine.
 */


/**
 * writes a \ref basic_grid_data "basic grid" to file
 * @param[in] grid basic grid data
 * @param[in] name file name (".vtk" will be added to this)
 * @remark as a reference, the equator will be included in the file
 */
void yac_write_basic_grid_to_file(struct basic_grid_data * grid, char * name);

/**
 * writes a list of cells into a vtk file, which can be visualised by paraview
 * @param[in] cells               list of cells
 * @param[in] num_cells           number of entries in cells
 * @param[in] name                file name (".vtk" will be added to this)
 * @param[in] num_points_per_edge each cell edge will be approximated by N
 *                                straight lines in 3D space, where N is
 *                                (num_points_per_edge - 1)
 * @remark the edge type will be taken into account
 */
void yac_write_grid_cells_to_file(
  struct grid_cell * cells, size_t num_cells, char * name,
  size_t num_points_per_edge);

#endif // GRID2VTK_H

