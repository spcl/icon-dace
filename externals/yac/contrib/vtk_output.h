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
#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

#include <stdio.h>

/** \example test_vtk_output.c
 * Test for vtk output utility functions.
 */

/** \file vtk_output.h
  * \brief general routines for writing vtk files
  *
  * To create a vtk file you have to execute the following
  * steps in the specified order:
  *
  * -# generate a vtk file (\ref vtk_open)
  * -# define the grid
  *   -# write the point data (\ref vtk_write_point_data)
  *   -# write the cell data (\ref vtk_write_cell_data)
  * -# provide the field data (these routines can be called in any order)
  *   -# field cell data
  *     -# \ref vtk_write_cell_scalars_uint
  *     -# \ref vtk_write_cell_scalars_int
  *     -# \ref vtk_write_cell_scalars_float
  *     -# \ref vtk_write_cell_scalars_double
  *   -# field point data
  *     -# \ref vtk_write_point_scalars_uint
  *     -# \ref vtk_write_point_scalars_int
  *     -# \ref vtk_write_point_scalars_float
  *     -# \ref vtk_write_point_scalars_double
  * -# close the vtk file (\ref vtk_close)
 **/

typedef struct VTK_FILE_ VTK_FILE;

/**
 * initialises a vtk file
 * @param[in] filename name of the vtk file
 * @param[in] title title for the data inside the file
 * @return handle the the vtk file
 */
VTK_FILE * vtk_open(const char * filename, const char * title);

/**
 * writes the 3d coordinates of all points into the vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] point_data array containing the 3d coordinates of all points
 * @param[in] num_points number of points to be written
 *
 * \remark the array associated to point_data should have the size 3 * num_points
 * \remark points[i*3+0], points[i*3+1] and points[i*3+2] should contain the x, y and z coordinate of the i'th point
 */
void vtk_write_point_data(VTK_FILE * vtk_file, double * point_data, unsigned num_points);

/**
 * writes the cell data (which cell consists of which points) to the vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] cell_corners contains for all cells the indices of points the respective cells are made up of
 * @param[in] num_points_per_cell contains contains for each cell the number of corners it is made up of
 * @param[in] num_cells number of cells
 *
 * \remark the indices in cell_corners refere to the index in points array passed to \ref vtk_write_point_data
 * \remark the size of num_points_per_cell should be num_cells
 * \remark the size of cell_corners is the sum of all num_points_per_cell[i] for 0<=i<num_cells
 */
void vtk_write_cell_data(VTK_FILE * vtk_file, unsigned * cell_corners,
                         unsigned * num_points_per_cell, unsigned num_cells);

/**
 * writes an array of unsigned integer scalar cell data values for each cell to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each cell
 * @param[in] num_cells number of cells
 * @param[in] name name of the written field
 *
 * /remark num_cells should be identically to the number of cells passed to a previous call to vtk_write_cell_data
 * /remark the size of scalars should be num_cells
 */
void vtk_write_cell_scalars_uint(VTK_FILE * vtk_file, unsigned * scalars,
                                 unsigned num_cells, char * name);

/**
 * writes an array of unsigned integer scalar point data values for each point to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each point
 * @param[in] num_points number of points
 * @param[in] name name of the written field
 *
 * /remark num_points should be identically to the number of points passed to a previous call to vtk_write_point_data
 * /remark the size of scalars should be num_points
 */
void vtk_write_point_scalars_uint(VTK_FILE * vtk_file, unsigned * scalars,
                                  unsigned num_points, char * name);

/**
 * writes an array of integer scalar cell data values for each cell to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each cell
 * @param[in] num_cells number of cells
 * @param[in] name name of the written field
 *
 * /remark num_cells should be identically to the number of cells passed to a previous call to vtk_write_cell_data
 * /remark the size of scalars should be num_cells
 */
void vtk_write_cell_scalars_int(VTK_FILE * vtk_file, int * scalars, unsigned num_cells,
                                char * name);

/**
 * writes an array of integer scalar point data values for each point to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each point
 * @param[in] num_points number of points
 * @param[in] name name of the written field
 *
 * /remark num_points should be identically to the number of points passed to a previous call to vtk_write_point_data
 * /remark the size of scalars should be num_points
 */
void vtk_write_point_scalars_int(VTK_FILE * vtk_file, int * scalars,
                                 unsigned num_points, char * name);

/**
 * writes an array of float scalar cell data values for each cell to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each cell
 * @param[in] num_cells number of cells
 * @param[in] name name of the written field
 *
 * /remark num_cells should be identically to the number of cells passed to a previous call to vtk_write_cell_data
 * /remark the size of scalars should be num_cells
 */
void vtk_write_cell_scalars_float(VTK_FILE * vtk_file, float * scalars,
                                  unsigned num_cells, char * name);

/**
 * writes an array of float scalar point data values for each point to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each point
 * @param[in] num_points number of points
 * @param[in] name name of the written field
 *
 * /remark num_points should be identically to the number of points passed to a previous call to vtk_write_point_data
 * /remark the size of scalars should be num_points
 */
void vtk_write_point_scalars_float(VTK_FILE * vtk_file, float * scalars,
                                   unsigned num_points, char * name);

/**
 * writes an array of double scalar cell data values for each cell to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each cell
 * @param[in] num_cells number of cells
 * @param[in] name name of the written field
 *
 * /remark num_cells should be identically to the number of cells passed to a previous call to vtk_write_cell_data
 * /remark the size of scalars should be num_cells
 */

void vtk_write_cell_scalars_double(VTK_FILE * vtk_file, double * scalars,
                                   unsigned num_cells, char * name);

/**
 * writes an array of double scalar point data values for each point to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each point
 * @param[in] num_points number of points
 * @param[in] name name of the written field
 *
 * /remark num_points should be identically to the number of points passed to a previous call to vtk_write_point_data
 * /remark the size of scalars should be num_points
 */
void vtk_write_point_scalars_double(VTK_FILE * vtk_file, double * scalars,
                                    unsigned num_points, char * name);

/**
 * closes a vtk file
 * @param[in] vtk_file file pointer to an already open file
 */
void vtk_close(VTK_FILE * vtk_file);

#endif // VTK_OUTPUT_H

