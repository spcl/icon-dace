/**
 * @file interp_grid.h
 *
 * @copyright Copyright  (C)  2019 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#ifndef INTERP_GRID_H
#define INTERP_GRID_H

#include "dist_grid.h"

/** \example test_interp_grid_parallel.c
 * A test for parallel grid interpolation.
 */

/**
 * generate a interpolation grid
 * @param[in] grid_pair      distributed grid pair
 * @param[in] src_grid_name  name of the source grid
 * @param[in] tgt_grid_name  name of the target grid
 * @param[in] num_src_fields number of source fields
 * @param[in] src_fields     specifies the source fields
 * @param[in] tgt_field      specifies the target field
 * @return interpolation grid
 */
struct interp_grid * yac_interp_grid_new(
  struct dist_grid_pair * grid_pair,
  char const * src_grid_name, char const * tgt_grid_name,
  size_t num_src_fields, struct interp_field const * src_fields,
  struct interp_field const tgt_field);

/**
 * gets all unmasked points available in the local part of the
 * distributed source grids
 * @param[in]  interp_grid   interpolation grid
 * @param[in]  src_field_idx index of the source field
 * @param[out] src_indices   local indices of all unmasked local source cells
 * @param[out] count         number of entries in src_indices
 * @remark each unmask point of the global source grid is returned only by
 *         a single process
 */
void yac_interp_grid_get_src_points(
  struct interp_grid * interp_grid, size_t src_field_idx,
  size_t ** src_indices, size_t * count);

/**
 * gets all unmasked points available in the local part of the
 * distributed target grids
 * @param[in]  interp_grid   interpolation grid
 * @param[out] tgt_indices   local indices of all unmasked local target cells
 * @param[out] count         number of entries in tgt_indices
 * @remark each unmask point of the global target grid is returned only by
 *         a single process
 */
void yac_interp_grid_get_tgt_points(
  struct interp_grid * interp_grid, size_t ** tgt_indices, size_t * count);

/**
 * allocates and returns remote_point informations for a list of selected
 * source points
 * @param[in] interp_grid   interpolation grid
 * @param[in] src_field_idx index of the source field
 * @param[in] src_points    local ids of selected source points
 * @param[in] count         number of entries in src_points
 * @return remote_point information for the selected source points
 * @remark the user is responsible for freeing the memory allocated for the
 *         remote_point information
 */
struct remote_point * yac_interp_grid_get_src_remote_points(
  struct interp_grid * interp_grid, size_t src_field_idx,
  size_t * src_points, size_t count);

/**
 * allocates and returns remote_point informations for a list of selected
 * source points
 * @param[in] interp_grid interpolation grid
 * @param[in] location    location for which the remote point info is requested
 * @param[in] src_points  local ids of selected source points
 * @param[in] count       number of entries in src_points
 * @return remote_point information for the selected source points
 * @remark the user is responsible for freeing the memory allocated for the
 *         remote_point information
 */
struct remote_point * yac_interp_grid_get_src_remote_points2(
  struct interp_grid * interp_grid, enum yac_location location,
  size_t * src_points, size_t count);

/**
 * allocates and returns remote_point informations for a list of selected
 * target points
 * @param[in] interp_grid   interpolation grid
 * @param[in] tgt_points    local ids of selected target points
 * @param[in] count         number of entries in tgt_points
 * @return remote_point information for the selected target points
 * @remark the user is responsible for freeing the memory allocated for the
 *         remote_point information
 */
struct remote_point * yac_interp_grid_get_tgt_remote_points(
  struct interp_grid * interp_grid, size_t * tgt_points, size_t count);

/**
 * returns the local ids for the provided global source ids
 * @param[in]  interp_grid    interpolation grid
 * @param[in]  src_field_idx  index of the source field
 * @param[in]  src_global_ids global ids of the selected source points
 * @param[in]  count          number of entries in src_global_ids
 * @param[out] src_local_ids  local ids of the selected source points
 * @remark the user has to ensure that the array associated to
 *         src_local_ids is big enough to hold enough elements
 * @remark in case one or more selected global source ids are not available
 *         in the local data, the local source data will be extended
 *         accordingly
 */
void yac_interp_grid_src_global_to_local(
  struct interp_grid * interp_grid, size_t src_field_idx,
  yac_int * src_global_ids, size_t count, size_t * src_local_ids);

/**
 * returns the local ids for the provided global target ids
 * @param[in]  interp_grid    interpolation grid
 * @param[in]  tgt_global_ids global ids of the selected target points
 * @param[in]  count          number of entries in tgt_global_ids
 * @param[out] tgt_local_ids  local ids of the selected target points
 * @remark the user has to ensure that the array associated to
 *         tgt_local_ids is big enough to hold enough elements
 * @remark in case one or more selected global target ids are not available
 *         in the local data, the local target data will be extended
 *         accordingly
 */
void yac_interp_grid_tgt_global_to_local(
  struct interp_grid * interp_grid, yac_int * tgt_global_ids,
  size_t count, size_t * tgt_local_ids);

/**
 * gets the location of the target field
 * @param[in] interp_grid interpolation grid
 * @return location of the target field
 */
enum yac_location yac_interp_grid_get_tgt_field_location(
  struct interp_grid * interp_grid);

/**
 * gets the global ids of the selected source points
 * @param[in]  interp_grid    interpolation grid
 * @param[in]  src_points     local ids of the selected source points
 * @param[in]  count          number of entries in src_points
 * @param[in]  src_field_idx  index of the source field
 * @param[in]  src_global_ids global ids of selected source points
 * @remark the user has to ensure that the array associated to
 *         src_global_ids is big enough to hold enough elements
 */
void yac_interp_grid_get_src_global_ids(
  struct interp_grid * interp_grid, size_t * src_points, size_t count,
  size_t src_field_idx, yac_int * src_global_ids);

/**
 * gets the global ids of the selected target points
 * @param[in]  interp_grid    interpolation grid
 * @param[in]  tgt_points     local ids of the selected target points
 * @param[in]  count          number of entries in tgt_points
 * @param[in]  tgt_global_ids global ids of selected target points
 * @remark the user has to ensure that the array associated to
 *         tgt_global_ids is big enough to hold enough elements
 */
void yac_interp_grid_get_tgt_global_ids(
  struct interp_grid * interp_grid, size_t * tgt_points, size_t count,
  yac_int * tgt_global_ids);

/**
 * gets the coordinates of the selected source points
 * @param[in]  interp_grid     interpolation grid
 * @param[in]  src_points      local ids of the selected source points
 * @param[in]  count           number of entries in src_points
 * @param[in]  src_field_idx   index of the source field
 * @param[in]  src_coordinates coordinates of selected source points
 * @remark the user has to ensure that the array associated to
 *         src_coordinates is big enough to hold enough elements
 * @remark if the source field location is CORNER and no field coordinates
 *         are available/defined, this routine will return a pointer to the
 *         vertex coordinates of the source grid
 */
void yac_interp_grid_get_src_coordinates(
  struct interp_grid * interp_grid, size_t * src_points, size_t count,
  size_t src_field_idx, coordinate_pointer src_coordinates);

/**
 * gets the coordinates of the selected target points
 * @param[in]  interp_grid     interpolation grid
 * @param[in]  tgt_points      local ids of the selected target points
 * @param[in]  count           number of entries in tgt_points
 * @param[in]  tgt_coordinates coordinates of selected target points
 * @remark the user has to ensure that the array associated to
 *         tgt_coordinates is big enough to hold enough elements
 * @remark if the target field location is CORNER and no field coordinates
 *         are available/defined, this routine will return a pointer to the
 *         vertex coordinates of the source grid
 */
void yac_interp_grid_get_tgt_coordinates(
  struct interp_grid * interp_grid, size_t * tgt_points, size_t count,
  coordinate_pointer tgt_coordinates);

/**
 * gets the number of source fields in the provided interpolation grid
 * @param[in] interp_grid interpolation grid
 * @return number of source fields
 */
size_t yac_interp_grid_get_num_src_fields(struct interp_grid * interp_grid);

/**
 * gets the location of a source field
 * @param[in] interp_grid   interpolation grid
 * @param[in] src_field_idx index of the source field
 * @return location of selected source field
 */
enum yac_location yac_interp_grid_get_src_field_location(
  struct interp_grid * interp_grid, size_t src_field_idx);

/**
 * gets the global ids of the selected field of all points in the local part
 * of the distributed source grid
 * @param[in] interp_grid   interpolation grid
 * @param[in] src_field_idx index of the source field
 * @return global ids for the selected source field
 */
const_yac_int_pointer yac_interp_grid_get_src_field_global_ids(
  struct interp_grid * interp_grid, size_t src_field_idx);

/**
 * gets the coordinates of the selected field of all points in the local part
 * of the distributed source grid
 * @param[in] interp_grid   interpolation grid
 * @param[in] src_field_idx index of the source field
 * @return coordinates for the selected source field
 * @remark if the source field location is CORNER and no field coordinates
 *         are available/defined, this routine will return a pointer to the
 *         vertex coordinates of the source grid
 * @remark if the source field location is not CORNER and no field coordinates
 *         are available/defined, this routine will return NULL
 */
const_coordinate_pointer yac_interp_grid_get_src_field_coords(
  struct interp_grid * interp_grid, size_t src_field_idx);

/**
 * gets the coordinates of the target field of all points in the local part
 * of the distributed target grid
 * @param[in] interp_grid   interpolation grid
 * @return coordinates for the target field
 * @remark if the target field location is CORNER and no field coordinates
 *         are available/defined, this routine will return a pointer to the
 *         vertex coordinates of the target grid
 * @remark if the target field location is not CORNER and no field coordinates
 *         are available/defined, this routine will return NULL
 */
const_coordinate_pointer yac_interp_grid_get_tgt_field_coords(
  struct interp_grid * interp_grid);

/**
 * gets the mask of the selected field of all points in the local part
 * of the distributed source grid
 * @param[in] interp_grid   interpolation grid
 * @param[in] src_field_idx index of the source field
 * @return mask for the selected source field
 * @remark if for the selected field no mask is available/defined,
 *         NULL is returned
 */
const_int_pointer yac_interp_grid_get_src_field_mask(
  struct interp_grid * interp_grid, size_t src_field_idx);

/**
 * search for source cells that map to the provided search coordinates
 * @param[in]  interp_grid   interpolation grid
 * @param[in]  search_coords search coordinates
 * @param[in]  count         number of search coordinates
 * @param[out] src_cells     local ids of matching source cells
 * @remark in case for a search coordinate no matching source cell was found,
 *         the respective entry in src_cells is SIZE_MAX
 */
void yac_interp_grid_do_points_search(
  struct interp_grid * interp_grid, coordinate_pointer search_coords,
  size_t count, size_t * src_cells);

/**
 * search for source cells that map to the provided search coordinates
 * @param[in]  interp_grid   interpolation grid
 * @param[in]  search_coords search coordinates
 * @param[in]  count         number of search coordinates
 * @param[out] src_cells     local ids of matching source cells
 * @remark in case for a search coordinate no matching source cell was found,
 *         the respective entry in src_cells is SIZE_MAX
 * @remark assumes that all grid edges are on great circle edges
 */
void yac_interp_grid_do_points_search_gc(
  struct interp_grid * interp_grid, coordinate_pointer search_coords,
  size_t count, size_t * src_cells);

/**
 * does a n-nearest-neighbour search on the source grid
 * @param[in]  interp_grid   interpolation grid
 * @param[in]  search_coords search coordinates
 * @param[in]  count         number of search coordinates
 * @param[in]  n             number of source points per search coordinate
 * @param[out] src_points    local ids of result source points
 */
void yac_interp_grid_do_nnn_search_src(
  struct interp_grid * interp_grid, coordinate_pointer search_coords,
  size_t count, size_t n, size_t * src_points);

/**
 * does a n-nearest-neighbour search on the target grid
 * @param[in]  interp_grid   interpolation grid
 * @param[in]  search_coords search coordinates
 * @param[in]  count         number of search coordinates
 * @param[in]  n             number of target points per search coordinate
 * @param[out] tgt_points    local ids of result target points
 */
void yac_interp_grid_do_nnn_search_tgt(
  struct interp_grid * interp_grid, coordinate_pointer search_coords,
  size_t count, size_t n, size_t * tgt_points);

/**
 * search for all source cells matching the provided bounding circles
 * @param[in]  interp_grid            interpolation grid
 * @param[in]  bnd_circles            search bounding circles
 * @param[in]  count                  number of search bounding circles
 * @param[in]  src_field_idx          index of the source field
 * @param[out] src_cells              local ids of result source points
 * @param[out] num_src_per_bnd_circle number of result source points per
 *                                    bounding circle
 */
void yac_interp_grid_do_bnd_circle_search_src(
  struct interp_grid * interp_grid, const_bounding_circle_pointer bnd_circles,
  size_t count, size_t src_field_idx, size_t ** src_cells,
  size_t * num_src_per_bnd_circle);

/**
 * search for all target cells matching the provided bounding circles
 * @param[in]  interp_grid            interpolation grid
 * @param[in]  bnd_circles            search bounding circles
 * @param[in]  count                  number of search bounding circles
 * @param[out] tgt_cells              local ids of result target points
 * @param[out] num_tgt_per_bnd_circle number of result target points per
 *                                    bounding circle
 */
void yac_interp_grid_do_bnd_circle_search_tgt(
  struct interp_grid * interp_grid, const_bounding_circle_pointer bnd_circles,
  size_t count, size_t ** tgt_cells, size_t * num_tgt_per_bnd_circle);

/**
 * search for all source cells matching the provided target cells
 * @param[in]  interp_grid     interpolation grid
 * @param[in]  tgt_cells       local ids of target cells
 * @param[in]  count           number of target cells
 * @param[out] src_cells       local ids of result source points
 * @param[out] num_src_per_tgt number of result source points per target cell
 */
void yac_interp_grid_do_cell_search_src(
  struct interp_grid * interp_grid, size_t * tgt_cells, size_t count,
  size_t ** src_cells, size_t * num_src_per_tgt);

/**
 * search for all target cells matching the provided source cells
 * @param[in]  interp_grid     interpolation grid
 * @param[in]  src_cells       local ids of source cells
 * @param[in]  count           number of source cells
 * @param[out] tgt_cells       local ids of result target points
 * @param[out] num_tgt_per_src number of result target points per source cell
 */
void yac_interp_grid_do_cell_search_tgt(
  struct interp_grid * interp_grid, size_t * src_cells, size_t count,
  size_t ** tgt_cells, size_t * num_tgt_per_src);

/**
 * gets a communicator containing all ranks of the interpolation grid
 * @param[in] interp_grid interpolation grid
 * @return comm communicator
 */
MPI_Comm yac_interp_grid_get_MPI_Comm(struct interp_grid * interp_grid);

/**
 * gets the basic grid data of the source grid
 * @param[in] interp_grid interpolation grid
 * @return basic grid data
 * @remark the contents of the basic grid data may change by some calls
 *         (for example by \ref yac_interp_grid_do_bnd_circle_search_src)
 */
struct const_basic_grid_data * yac_interp_grid_get_basic_grid_data_src(
  struct interp_grid * interp_grid);

/**
 * gets the basic grid data of the target grid
 * @param[in] interp_grid interpolation grid
 * @return basic grid data
 * @remark the contents of the basic grid data may change by some calls
 *         (for example by \ref yac_interp_grid_do_bnd_circle_search_src)
 */
struct const_basic_grid_data * yac_interp_grid_get_basic_grid_data_tgt(
  struct interp_grid * interp_grid);

/**
 * gets neighbouring source grid cells for all provided source cells
 * @param[in]  interp_grid interpolation grid
 * @param[in]  src_cells   local ids of source cells for which the neighbours
 *                         are to be determined
 * @param[in]  count       number for entries in src_cells
 * @param[out] neighbours  neighbours for all edges of all provided source cells
 * @remark the user has to ensure that the array neighbours has the
 *         appropriate size, which is the sum of the number of edges for each
 *         provided cell
 * @remark for each provided source cell, this routine sets the
 *         neighbours for all edges of the cells
 * @remark in case there is no neighbour of an edge, the respective entry
 *         contains SIZE_MAX
 * @remark the result cells are sorted in clock- or counter
 *         clockwise order
 */
void yac_interp_grid_get_src_cell_neighbours(
  struct interp_grid * interp_grid, size_t * src_cells, size_t count,
  size_t * neighbours);

/**
 * gets neighbouring target grid cells for all provided target cells
 * @param[in]  interp_grid interpolation grid
 * @param[in]  tgt_cells   local ids of target cells for which the neighbours
 *                         are to be determined
 * @param[in]  count       number for entries in tgt_cells
 * @param[out] neighbours  neighbours for all edges of all provided target cells
 * @remark the user has to ensure that the array neighbours has the
 *         appropriate size, which is the sum of the number of edges for each
 *         provided cell
 * @remark for each provided source cell, this routine sets the
 *         neighbours for all edges of the cells
 * @remark in case there is no neighbour of an edge, the respective entry
 *         contains SIZE_MAX
 * @remark the result cells are sorted in clock- or counter
 *         clockwise order
 */
void yac_interp_grid_get_tgt_cell_neighbours(
  struct interp_grid * interp_grid, size_t * tgt_cells, size_t count,
  size_t * neighbours);

/**
 * generates for all vertices of the source grid the list of all cells
 * surrounding the vertices\n
 * vertices not belonging to the selected cells will have a zero entry in
 * num_cells_per_vertex
 * @param[in]  interp_grid            interpolation grid
 * @param[in]  cells                  selected source cells
 * @param[in]  count                  number of entries in cells
 * @param[out] vertex_to_cell         list of result cells for all vertices
 *                                    associated to the selected cells
 * @param[out] vertex_to_cell_offsets the results cells for vertex i can be
 *                                    found at:\n
 *                                    vertex_to_cell + vertex_to_cell_offsets[i]
 * @param[out] num_cells_per_vertex   number of result cells for each vertex
 * @remark the result cells for a vertex are sorted in clock- or counter
 *         clockwise order
 * @remark in case a vertex is not fully surrounded by unmasked cells,
 *         the respective entry in num_cells_per_vertex is 0
 * @remark the user is responsible to free the memory returned through the
 *         out arrays
 */
void yac_interp_grid_get_aux_grid_src(
  struct interp_grid * interp_grid, size_t * cells, size_t count,
  size_t ** vertex_to_cell, size_t ** vertex_to_cell_offsets,
  int ** num_cells_per_vertex);

/**
 * Relocates source-target-point-pairs. The flag "to_tgt_owner" determines
 * whether the pairs go to the dist owner process of the target or source
 * points for each pair.
 * @param[in]    interp_grid   interpolation grid
 * @param[in]    to_tgt_owner  determines whether the pairs go to the dist
 *                             owner process of the target or source points
 * @param[in]    src_field_idx index of the source field
 * @param[inout] src_points    local ids of the source-part of the point pairs
 * @param[inout] tgt_points    local ids of the target-part of the point pairs
 * @param[inout] weights       weights for all point pairs
 * @param[inout] count         number of point pairs before and after
 *                             this routine
 * @remark if the pointer weights points to NULL on all processes, no weights
 *         will be exchanged
 */
void yac_interp_grid_relocate_src_tgt_pairs(
  struct interp_grid * interp_grid, int to_tgt_owner,
  size_t src_field_idx, size_t ** src_points, size_t ** tgt_points,
  double ** weights, size_t * count);

/**
 * Gets the distributed owners of the given target indices
 * @param[in]  interp_grid interpolation grid
 * @param[in]  tgt_indices local indices of target points
 * @param[in]  count       number of indices in tgt_indices
 * @param[out] owners      ranks of processes to which the target points in
 *                         tgt_indices are assigned to by the distributed grid
 *                         decomposition
 * @remark the user has to ensure that owners is big enough to hold the results
 *         for all target points
 */
void yac_interp_grid_determine_dist_tgt_owners(
  struct interp_grid * interp_grid, size_t * tgt_indices, size_t count,
  int * owners);

/**
 * Determines all non-masked neighour vertices for the selected
 * target vertices
 * @param[in]  interp_grid           interpolation grid
 * @param[in]  vertices              local ids of selected vertices
 * @param[in]  count                 number of entries in vertices
 * @param[out] neigh_vertices        neighbour vertices
 * @param[out] num_neighs_per_vertex number of neighbours per vertex
 * @remark the user has to ensure that num_neighs_per_vertex is allocated
 *         big enough
 * @remark the user is responsible for freeing memory associated to
 *         neigh_vertices
 */
void yac_interp_grid_get_tgt_vertex_neighbours(
  struct interp_grid * interp_grid, size_t * vertices, size_t count,
  size_t ** neigh_vertices, int * num_neighs_per_vertex);

/**
 * Relocates source-target-point-pairs. The flag "to_tgt_owner" determines
 * whether the pairs go to the orig owner process of the target or source
 * points for each pair.
 * @param[in]    interp_grid  interpolation grid
 * @param[in]    to_tgt_owner determines whether the pairs go to the dist
 *                            owner process of the target or source points
 * @param[in]    src_location location of the source points
 * @param[inout] src_points   local ids of the source-part of the point pairs
 * @param[inout] tgt_points   local ids of the target-part of the point pairs
 * @param[inout] weights      weights for all point pairs
 * @param[inout] count        number of point pairs before and after
 *                            this routine
 * @remark if the pointer weights points to NULL on all processes, no weights
 *         will be exchanged
 */
void yac_interp_grid_relocate_src_tgt_pairs_orig(
  struct interp_grid * interp_grid, int to_tgt_owner,
  enum yac_location src_location, size_t ** src_points, size_t ** tgt_points,
  double ** weights, size_t * count);

/**
 * deletes all memory associated with the provided interpolation grid
 * @param[inout] interp_grid interpolation grid
 */
void yac_interp_grid_delete(struct interp_grid * interp_grid);

#endif // INTERP_GRID_H
