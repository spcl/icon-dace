/**
 * @file grid_cell.h
 * @brief Structs and interfaces to handle grid cells
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
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

#include <stdio.h>

#ifndef GRID_CELL_H
#define GRID_CELL_H

enum yac_edge_type {
   GREAT_CIRCLE_EDGE = 0, //!< great circle
   LAT_CIRCLE_EDGE   = 1, //!< latitude circle
   LON_CIRCLE_EDGE   = 2, //!< longitude circle
};

struct grid_cell {
   double (*coordinates_xyz)[3];
   enum yac_edge_type * edge_type;
   size_t num_corners;
   size_t array_size; //!< size in elements of the arrays: coordinates_x,
                      //!< coordinates_y, edge_type and 1/3 of coordinates_xyz
};

enum yac_cell_type {
  LON_LAT_CELL,
  LAT_CELL,
  GREAT_CIRCLE_CELL,
  MIXED_CELL
};

/**
 * initiates a grid_cell object
 * before the first being used a grid_cell object has to be initialised
 * @param[in] cell object to be initialised
 * @see free_grid_cell
 * @see get_grid_cell
 */
void yac_init_grid_cell(struct grid_cell * cell);

/**
 * copies a given grid cell
 * @param[in]  in_cell  cell to be copied
 * @param[out] out_cell copied cell
 * @remarks out_cell needs to be a cell that has previously been
 *          initialised or a cell that already contains valid data
 */
void yac_copy_grid_cell(struct grid_cell in_cell, struct grid_cell * out_cell);

/**
 * frees all memory associated with a grid_cell object and reinitialised
 * the cell
 * @param[in,out] cell
 */
void yac_free_grid_cell(struct grid_cell * cell);

#ifdef YAC_DEBUG_GRIC_CELL
/**
 * prints out info about a grid_cell object and reinitialised
 * the cell
 * @param[in] stream
 * @param[in] cell
 * @param[in] name
 */
void print_grid_cell(FILE * stream, struct grid_cell cell, char * name);
#endif

#endif // GRID_CELL_H

