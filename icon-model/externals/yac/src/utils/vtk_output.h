// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

// YAC PUBLIC HEADER START

/** \example test_vtk_output.c
 * Test for vtk output utility functions.
 */

/** \file vtk_output.h
  * \brief general routines for writing vtk files
  *
  * To create a vtk file you have to execute the following
  * steps in the specified order:
  *
  * -# generate a vtk file (\ref yac_vtk_open)
  * -# define the grid
  *   -# write the point data (\ref yac_vtk_write_point_data)
  *   -# write the cell data (\ref yac_vtk_write_cell_data)
  * -# provide the field data (these routines can be called in any order)
  *   -# field cell data
  *     -# \ref yac_vtk_write_cell_scalars_uint
  *     -# \ref yac_vtk_write_cell_scalars_int
  *     -# \ref yac_vtk_write_cell_scalars_float
  *     -# \ref yac_vtk_write_cell_scalars_double
  *   -# field point data
  *     -# \ref yac_vtk_write_point_scalars_uint
  *     -# \ref yac_vtk_write_point_scalars_int
  *     -# \ref yac_vtk_write_point_scalars_float
  *     -# \ref yac_vtk_write_point_scalars_double
  * -# close the vtk file (\ref yac_vtk_close)
 **/

typedef struct YAC_VTK_FILE_ YAC_VTK_FILE;

/**
 * initialises a vtk file
 * @param[in] filename name of the vtk file
 * @param[in] title title for the data inside the file
 * @return handle the the vtk file
 */
YAC_VTK_FILE * yac_vtk_open(const char * filename, const char * title);

/**
 * writes the 3d coordinates of all points into the vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] point_data array containing the 3d coordinates of all points
 * @param[in] num_points number of points to be written
 *
 * \remark the array associated to point_data should have the size 3 * num_points
 * \remark points[i*3+0], points[i*3+1] and points[i*3+2] should contain the x, y and z coordinate of the i'th point
 */
void yac_vtk_write_point_data(
  YAC_VTK_FILE * vtk_file, double * point_data, unsigned num_points);

/**
 * writes the 3d coordinates of all points into the vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] point_data_lon array containing the longitude coordinates of all points (in rad)
 * @param[in] point_data_lat array containing the latitude coordinates of all points (in rad)
 * @param[in] num_points number of points to be written
 *
 * \remark the array associated to point_data should have the size 3 * num_points
 * \remark points[i*3+0], points[i*3+1] and points[i*3+2] should contain the x, y and z coordinate of the i'th point
 */
void yac_vtk_write_point_data_ll(
  YAC_VTK_FILE * vtk_file, double * point_data_lon, double * point_data_lat,
  unsigned num_points);

/**
 * writes the cell data (which cell consists of which points) to the vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] cell_corners contains for all cells the indices of points the respective cells are made up of
 * @param[in] num_points_per_cell contains contains for each cell the number of corners it is made up of
 * @param[in] num_cells number of cells
 *
 * \remark the indices in cell_corners refere to the index in points array passed to
 *         \ref yac_vtk_write_point_data
 * \remark the size of num_points_per_cell should be num_cells
 * \remark the size of cell_corners is the sum of all num_points_per_cell[i] for 0<=i<num_cells
 */
void yac_vtk_write_cell_data(
  YAC_VTK_FILE * vtk_file, unsigned * cell_corners,
  unsigned * num_points_per_cell, unsigned num_cells);

/**
 * writes an array of unsigned integer scalar cell data values for each cell to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each cell
 * @param[in] num_cells number of cells
 * @param[in] name name of the written field
 *
 * /remark num_cells should be identically to the number of cells passed to a
 *         previous call to yac_vtk_write_cell_data
 * /remark the size of scalars should be num_cells
 */
void yac_vtk_write_cell_scalars_uint(
  YAC_VTK_FILE * vtk_file, unsigned * scalars,
  unsigned num_cells, char * name);

/**
 * writes an array of unsigned integer scalar point data values for each point to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each point
 * @param[in] num_points number of points
 * @param[in] name name of the written field
 *
 * /remark num_points should be identically to the number of points passed to a
 *         previous call to yac_vtk_write_point_data
 * /remark the size of scalars should be num_points
 */
void yac_vtk_write_point_scalars_uint(
  YAC_VTK_FILE * vtk_file, unsigned * scalars,
  unsigned num_points, char * name);

/**
 * writes an array of integer scalar cell data values for each cell to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each cell
 * @param[in] num_cells number of cells
 * @param[in] name name of the written field
 *
 * /remark num_cells should be identically to the number of cells passed to a
 *         previous call to yac_vtk_write_cell_data
 * /remark the size of scalars should be num_cells
 */
void yac_vtk_write_cell_scalars_int(
  YAC_VTK_FILE * vtk_file, int * scalars, unsigned num_cells, char * name);

/**
 * writes an array of integer scalar point data values for each point to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each point
 * @param[in] num_points number of points
 * @param[in] name name of the written field
 *
 * /remark num_points should be identically to the number of points passed to
 *         a previous call to yac_vtk_write_point_data
 * /remark the size of scalars should be num_points
 */
void yac_vtk_write_point_scalars_int(
  YAC_VTK_FILE * vtk_file, int * scalars,
  unsigned num_points, char * name);

/**
 * writes an array of float scalar cell data values for each cell to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each cell
 * @param[in] num_cells number of cells
 * @param[in] name name of the written field
 *
 * /remark num_cells should be identically to the number of cells passed to
 *         a previous call to yac_vtk_write_cell_data
 * /remark the size of scalars should be num_cells
 */
void yac_vtk_write_cell_scalars_float(
  YAC_VTK_FILE * vtk_file, float * scalars,
  unsigned num_cells, char * name);

/**
 * writes an array of float scalar point data values for each point to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each point
 * @param[in] num_points number of points
 * @param[in] name name of the written field
 *
 * /remark num_points should be identically to the number of points passed to a
 *         previous call to yac_vtk_write_point_data
 * /remark the size of scalars should be num_points
 */
void yac_vtk_write_point_scalars_float(
  YAC_VTK_FILE * vtk_file, float * scalars,
  unsigned num_points, char * name);

/**
 * writes an array of double scalar cell data values for each cell to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each cell
 * @param[in] num_cells number of cells
 * @param[in] name name of the written field
 *
 * /remark num_cells should be identically to the number of cells passed to a
 *         previous call to yac_vtk_write_cell_data
 * /remark the size of scalars should be num_cells
 */

void yac_vtk_write_cell_scalars_double(
  YAC_VTK_FILE * vtk_file, double * scalars,
  unsigned num_cells, char * name);

/**
 * writes an array of double scalar point data values for each point to a vtk file
 * @param[in] vtk_file file pointer to an already open file
 * @param[in] scalars and array containing a scalar value for each point
 * @param[in] num_points number of points
 * @param[in] name name of the written field
 *
 * /remark num_points should be identically to the number of points passed to a
 *         previous call to yac_vtk_write_point_data
 * /remark the size of scalars should be num_points
 */
void yac_vtk_write_point_scalars_double(
  YAC_VTK_FILE * vtk_file, double * scalars,
  unsigned num_points, char * name);

/**
 * closes a vtk file
 * @param[in] vtk_file file pointer to an already open file
 */
void yac_vtk_close(YAC_VTK_FILE * vtk_file);

// YAC PUBLIC HEADER STOP

#endif // VTK_OUTPUT_H

