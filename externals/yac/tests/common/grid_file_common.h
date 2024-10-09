// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef GRID_FILE_COMMON_H
#define GRID_FILE_COMMON_H

/**
 * checks whether data in a scrip formated grid file is identical to a
 * reference data provided
 * @param[in] filename                 name of the grid file to be checked
 * @param[in] grid_name                name of the grid
 * @param[in] ref_num_cells            reference number of cells
 * @param[in] ref_num_corners_per_cell reference number of corners per cell
 * @param[in] ref_cla                  reference latitude coordinates in degree
 * @param[in] ref_clo                  reference longitude coordinates in degree
 * @param[in] ref_lat                  reference latitude cell center coordinates in degree
 * @param[in] ref_lon                  reference longitude cell center coordinates in degree
 * @param[in] ref_cell_global_ids      reference global cell ids
 * @param[in] ref_core_cell_mask       reference core cell mask
 * @param[in] ref_vertex_global_ids    reference global vertex ids
 * @param[in] ref_core_vertex_mask     reference core vertex mask
 * @param[in] ref_edge_global_ids      reference global edge ids
 * @param[in] ref_core_edge_mask       reference core edge mask
 * @remark arguments ref_lat, ref_lon, ref_global_ids, ref_core_cell_mask can be NULL
 */
void check_grid_file(
  char const * filename, char const * grid_name,
  size_t ref_num_cells, size_t ref_num_corners_per_cell,
  double * ref_cla, double * ref_clo, double * ref_lat, double * ref_lon,
  int * ref_cell_global_ids, int * ref_core_cell_mask,
  int * ref_vertex_global_ids, int * ref_core_vertex_mask,
  int * ref_edge_global_ids, int * ref_core_edge_mask);

#endif // GRID_FILE_COMMON_H
