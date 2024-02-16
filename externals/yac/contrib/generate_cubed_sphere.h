/**
 * @file generate_cubed_sphere.h
 * @brief contains routines to generate cubed sphere grids
 *
 * @copyright Copyright  (C)  2015 Moritz Hanke <hanke@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
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


#ifndef GENERATE_CUBED_SPHERE_H
#define GENERATE_CUBED_SPHERE_H

#include "grid.h"

/** \example test_generate_cubed_sphere.c
 * A test for the generation of cubed sphere grids used in several toy models.
 */

void generate_cubed_sphere_grid_information(
  unsigned n, unsigned * num_cells, unsigned * num_vertices,
  double ** x_vertices, double ** y_vertices, double ** z_vertices,
  unsigned ** vertices_of_cell, unsigned ** face_id);

void generate_part_cube_grid_information(
  unsigned n, unsigned * nbr_vertices, unsigned * nbr_cells,
  unsigned ** num_vertices_per_cell, unsigned ** cell_to_vertex,
  double ** x_vertices, double ** y_vertices, double ** x_cells,
  double ** y_cells, int ** global_cell_id, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size);

struct basic_grid_data generate_cubed_sphere_grid(unsigned n);

#endif // GENERATE_CUBED_SPHERE_H

