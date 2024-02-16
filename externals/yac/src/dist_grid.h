/**
 * @file dist_grid.h
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

#ifndef DIST_GRID_H
#define DIST_GRID_H

#include <yaxt.h>

#include "geometry.h"
#include "proc_sphere_part.h"
#include "grid.h"

/** \example test_dist_grid_pair_parallel.c
 * A test for distributed grid pairs.
 */

struct proc_sphere_part_node;
struct dist_grid_pair;
struct dist_grid;
typedef size_t const (* const const_size_t_2_pointer)[2];
typedef int const * const const_int_pointer;
typedef size_t const * const const_size_t_pointer;
typedef yac_int const * const const_yac_int_pointer;
typedef enum yac_edge_type const * const const_yac_edge_type_pointer;
typedef struct bounding_circle const * const const_bounding_circle_pointer;
typedef struct remote_point_infos const * const const_remote_point_infos_pointer;

struct const_basic_grid_data {
  const_coordinate_pointer vertex_coordinates;
  const_yac_int_pointer ids[3];
  const_int_pointer num_vertices_per_cell;
  const_size_t_pointer cell_to_vertex;
  const_size_t_pointer cell_to_vertex_offsets;
  const_size_t_pointer cell_to_edge;
  const_size_t_pointer cell_to_edge_offsets;
  const_size_t_2_pointer edge_to_vertex;
  const_bounding_circle_pointer cell_bnd_circles;
  const_yac_edge_type_pointer edge_type;
  const_remote_point_infos_pointer cell_owners;
  const_remote_point_infos_pointer vertex_owners;
  const_remote_point_infos_pointer edge_owners;
  const size_t count[3];
};

// information about a point on a remote process
struct remote_point_info {
  int rank;          // MPI rank
  uint64_t orig_pos; // local id on the remote process
};

struct remote_point_infos {
  int count;
  union {
    struct remote_point_info single; // valid if count == 1;
    struct remote_point_info * multi; // valid if count != -1;
  } data;
};

struct remote_point {
  yac_int global_id;
  struct remote_point_infos data;
};

struct remote_points {
  struct remote_point * data;
  size_t count;
  struct remote_point_info buffer[];
};

/**
 * generate a distributed grid pair
 * @param[in] grid_a      first basic grid
 * @param[in] grid_b      second basic grid
 * @param[in] comm        communicator containing the ranks of all processes
 *                        participating in this call
 * @return distributed grid pair
 * @remark this routine is collective for all processes in comm
 */
struct dist_grid_pair * yac_dist_grid_pair_new(
  struct yac_basic_grid * grid_a, struct yac_basic_grid * grid_b,
  MPI_Comm comm);

/**
 * delete all memory associated to the provided grid pair
 * @param[inout] grid_pair
 */
void yac_dist_grid_pair_delete(struct dist_grid_pair * grid_pair);

/**
 * gets a communicator containing all ranks of the distributed grid pair
 * @param[in] grid_pair distributed grid pair
 * @return comm communicator
 */
MPI_Comm yac_dist_grid_pair_get_MPI_Comm(struct dist_grid_pair * grid_pair);

/**
 * gets one of the two grids from the distributed grid pairs
 * @param[in] grid_pair distributed grid pair
 * @param[in] grid_name grid name of the selected grid
 * @return distributed grid
 */
struct dist_grid * yac_dist_grid_pair_get_dist_grid(
  struct dist_grid_pair * grid_pair, char const * grid_name);

/**
 * search for cells that map to the provided search coordinates
 * @param[in]  grid_pair     distributed grid pair
 * @param[in]  grid_name     grid name of the grid for which the search is to
 *                           be performed
 * @param[in]  search_coords search coordinates
 * @param[in]  count         number of search coordinates
 * @param[out] cells         local ids of matching cells
 * @remark in case for a search coordinate no matching cell was found,
 *         the respective entry in cells is SIZE_MAX
 */
void yac_dist_grid_pair_do_point_search(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  coordinate_pointer search_coords, size_t count, size_t * cells);

/**
 * search for cells that map to the provided search coordinates
 * @param[in]  grid_pair     distributed grid pair
 * @param[in]  grid_name     grid name of the grid for which the search is to
 *                           be performed
 * @param[in]  search_coords search coordinates
 * @param[in]  count         number of search coordinates
 * @param[out] cells         local ids of matching cells
 * @remark SIZE_MAX is the returned result in case no cell was found
 * @remark assumes that all grid edges are on great circle edges
 */
void yac_dist_grid_pair_do_point_search_gc(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  coordinate_pointer search_coords, size_t count, size_t * cells);

/**
 * does a n-nearest-neighbour search
 * @param[in]  grid_pair     distributed grid pair
 * @param[in]  grid_name     grid name of the grid for which the search is to
 *                           be performed
 * @param[in]  search_coords search coordinates
 * @param[in]  count         number of search coordinates
 * @param[out] local_ids     local ids of results points
 * @param[in]  n             number of points per search coordinate
 * @param[in]  field         field description
 */
void yac_dist_grid_pair_do_nnn_search(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  coordinate_pointer search_coords, size_t count, size_t * local_ids,
  size_t n, struct interp_field field);

/**
 * search for all cells matching the provided bounding circles
 * @param[in]  grid_pair                  distributed grid pair
 * @param[in]  grid_name                  grid name of the grid for which the
 *                                        search is to be performed
 * @param[in]  bnd_circles                search bounding circles
 * @param[in]  count                      number of search bounding circles
 * @param[out] cells                      local ids of results cells
 * @param[out] num_results_per_bnd_circle number of results per bounding circle
 * @param[in]  field                      field description
 */
void yac_dist_grid_pair_do_bnd_circle_search(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  const_bounding_circle_pointer bnd_circles, size_t count, size_t ** cells,
  size_t * num_results_per_bnd_circle, struct interp_field field);

/**
 * search for matches between cells of the two grids within the grid pair
 * @param[in]  grid_pair                   distributed grid pair
 * @param[in]  search_grid_name            grid name of the grid from which to
 *                                         take the search cells
 * @param[in]  result_grid_name            grid name of the grid from which to
 *                                         take the matching cells
 * @param[in]  search_cells                local ids of search cells
 * @param[in]  count                       number of search cells
 * @param[out] result_cells                local ids of results cells
 * @param[out] num_results_per_search_cell number of results per search cell
 * @param[in]  result_field                field description
 */
void yac_dist_grid_pair_do_cell_search(
  struct dist_grid_pair * grid_pair,
  char const * search_grid_name, char const * result_grid_name,
  size_t * search_cells, size_t count, size_t ** result_cells,
  size_t * num_results_per_search_cell, struct interp_field result_field);

/**
 * gets neighbouring grid cells for all provided cells
 * @param[in]  grid_pair  distributed grid pair
 * @param[in]  grid_name  grid name of the grid for which the
 *                        search is to be performed
 * @param[in]  cells      local ids of cells for which the neighbours are to
 *                        be determined
 * @param[in]  count      number for entries in cells
 * @param[out] neighbours neighbours for all edges of all provided cells
 * @remark the user has to ensure that the array neighbours has the
 *         appropriate size, which is the sum of the number of edges for each
 *         provided cell
 * @remark for each provided cell, this routine sets the
 *         neighbours for all edges of the cells
 * @remark in case there is no neighbour of an edge, the respective entry
 *         contains SIZE_MAX
 * @remark the result cells are sorted in clock- or counter
 *         clockwise order
 */
void yac_dist_grid_pair_get_cell_neighbours(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  size_t * cells, size_t count, size_t * neighbours);

/**
 * gets the basic grid data of a distributed grid
 * @param[in] dist_grid distributed grid
 * @return basic grid data
 * @remark the contents of the basic grid data may change by some calls
 *         (for example by \ref yac_dist_grid_pair_do_bnd_circle_search)
 */
struct const_basic_grid_data *
yac_dist_grid_get_basic_grid_data(struct dist_grid * dist_grid);

/**
 * gets the number of points in the local part of the distributed grid
 * @param[in] dist_grid distributed grid
 * @param[in] location  location of the points
 * @remark each global point is only counted once by a single process
 */
size_t yac_dist_grid_get_local_count(
  struct dist_grid * dist_grid, enum yac_location location);

/**
 * gets all unmasked points available in the local part of the distributed grid
 * @param[in]  grid    distributed grid
 * @param[in]  field   field description
 * @param[out] indices local indices of all unmasked local cells
 * @param[out] count   number of entries in indices
 * @remark each unmask point of the global grid is returned only by a single
 *         process
 */
void yac_dist_grid_get_local_unmasked_points(
  struct dist_grid * grid, struct interp_field field,
  size_t ** indices, size_t * count);

/**
 * generates an MPI Datatype for struct remote_point_info
 * @param[in] comm communicator
 * @return MPI Datatype for struct remote_point_info
 * @remark the user has to free the returned MPI Datatype using MPI_Type_free
 */
MPI_Datatype yac_get_remote_point_info_mpi_datatype(MPI_Comm comm);

/**
 * computes the maximum size required by MPI to pack the provided point
 * of type struct remote_point
 * @param[in] point         point for which the pack size is to be determined
 * @param[in] point_info_dt MPI Datatype for packing struct point_info
 * @param[in] comm          communicator
 * @return maximum packing size
 */
int yac_remote_point_get_pack_size(
  struct remote_point * point, MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * packs a provided remote_point; this is simlar to MPI_Pack
 * @param[in]    point         remote_point to be packed
 * @param[inout] buffer        packing buffer
 * @param[in]    buffer_size   size of packing buffer
 * @param[inout] position      packing position
 * @param[in]    point_info_dt MPI Datatype for packing struct point_info
 * @param[in]    comm          communicator
 */
void yac_remote_point_pack(
  struct remote_point * point, void * buffer, int buffer_size, int * position,
  MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * unpacks and allocates a remote_point from a buffer; this is similar to
 * MPI_Unpack
 * @param[in]    buffer        packing buffer
 * @param[in]    buffer_size   size of packing buffer
 * @param[inout] position      unpacking position
 * @param[out]   point         unpacked point
 * @param[in]    point_info_dt MPI Datatype for unpacking struct point_info
 * @param[in]    comm          communicator
 */
void yac_remote_point_unpack(
  void * buffer, int buffer_size, int * position, struct remote_point * point,
  MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * computes the maximum size required by MPI to pack the provided points
 * of type struct remote_points
 * @param[in] points        points for which the pack size is to be determined
 * @param[in] point_info_dt MPI Datatype for unpacking struct point_info
 * @param[in] comm          communicator
 * @return maximum packing size
 */
int yac_remote_points_get_pack_size(
  struct remote_points * points, MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * packs provided remote_points; this is simlar to MPI_Pack
 * @param[in]    points        remote_points to be packed
 * @param[inout] buffer        packing buffer
 * @param[in]    buffer_size   size of packing buffer
 * @param[inout] position      packing position
 * @param[in]    point_info_dt MPI Datatype for packing struct point_info
 * @param[in]    comm          communicator
 */
void yac_remote_points_pack(
  struct remote_points * points, void * buffer, int buffer_size, int * position,
  MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * unpacks and allocates remote_points from a buffer; this is similar to
 * MPI_Unpack
 * @param[in]    buffer        packing buffer
 * @param[in]    buffer_size   size of packing buffer
 * @param[inout] position      unpacking position
 * @param[out]   points        unpacked points
 * @param[in]    point_info_dt MPI Datatype for unpacking struct point_info
 * @param[in]    comm          communicator
 */
void yac_remote_points_unpack(
  void * buffer, int buffer_size, int * position,
  struct remote_points ** points, MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * allocates and returns remote_point informations for a list of selected points
 * @param[in] dist_grid distributed grid
 * @param[in] location  location of the requested points
 * @param[in] points    local ids of selected points
 * @param[in] count     number of entries in points
 * @return remote_point information for the selected points
 * @remark the user is responsible for freeing the memory allocated for the
 *         remote_point information
 */
struct remote_point * yac_dist_grid_get_remote_points(
  struct dist_grid * dist_grid, enum yac_location location,
  size_t * points, size_t count);

/**
 * returns the local ids for the provided global ids
 * @param[in]  dist_grid  distributed grid
 * @param[in]  location   location of the selected points
 * @param[in]  global_ids global ids of the selected points
 * @param[in]  count      number of entries in global_ids
 * @param[out] local_ids  local ids of the selected points
 * @remark the user has to ensure that the array associated to local_ids is
 *         big enough to hold enough elements
 * @remark in case one or more selected global ids are not available in the
 *         local data, the local data will be extended accordingly
 */
void yac_dist_grid_global_to_local(
  struct dist_grid * dist_grid, enum yac_location location,
  yac_int * global_ids, size_t count, size_t * local_ids);

/**
 * generates for all vertices of the grid the list of all cells
 * surrounding the vertices\n
 * vertices not belonging to the selected cells will have a zero entry in
 * num_cells_per_vertex
 * @param[in]  grid_pair              distributed grid pair
 * @param[in]  grid_name              grid name of the grid on which this
 *                                    routine is supposed to work on
 * @param[in]  cells                  selected cells
 * @param[in]  count                  number of entries in cells
 * @param[out] vertex_to_cell         list of result cells for all vertices
 *                                    associated to the selected cells
 * @param[out] vertex_to_cell_offsets the results cells for vertex i can be
 *                                    found at:\n
 *                                    vertex_to_cell + vertex_to_cell_offsets[i]
 * @param[out] num_cells_per_vertex   number of result cells for each vertex
 * @param[in]  field                  if the provided field contains a mask
 *                                    for cells, it will be used
 * @remark the result cells for a vertex are sorted in clock- or counter
 *         clockwise order
 * @remark in case a vertex is not fully surrounded by unmasked cells,
 *         the respective entry in num_cells_per_vertex is 0
 * @remark the user is responsible to free the memory returned through the
 *         out arrays
 */
void yac_dist_grid_pair_get_aux_grid(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  size_t * cells, size_t count,
  size_t ** vertex_to_cell, size_t ** vertex_to_cell_offsets,
  int ** num_cells_per_vertex, struct interp_field field);

/**
 * Each point in a distributed grid has a unique owner rank.\n
 * This routine relocates point pairs of two grids. If (a_is_ref != 0), this
 * routine will relocate the point pairs such that each process only has pairs
 * with points_a being owned by the respective process.
 * @param[in]    grid_pair     distributed grid pair
 * @param[in]    a_is_ref      determines whether points_a or points_b are
 *                             used as the sorting reference
 * @param[in]    to_dist_owner determines whether the sorting is based on
 *                             the dist or orig owner of the respective points
 * @param[in]    grid_name_a   grid name for the a-part of point pairs
 * @param[inout] points_a      local ids of the a-part of the point pairs
 * @param[in]    location_a    location of the a-part of the point pairs
 * @param[in]    grid_name_b   grid name for the b-part of point pairs
 * @param[inout] points_b      local ids of the b-part of the point pairs
 * @param[in]    location_b    location of the b-part of the point pairs
 * @param[inout] weights       weights for all point pairs
 * @param[inout] count         number of point pairs before and after this
 *                             routine
 */
void yac_dist_grid_pair_relocate_point_pairs(
  struct dist_grid_pair * grid_pair, int a_is_ref, int to_dist_owner,
  char const * grid_name_a, size_t ** points_a, enum yac_location location_a,
  char const * grid_name_b, size_t ** points_b, enum yac_location location_b,
  double ** weights, size_t * count);

/**
 * returns a pointer to a field mask
 * @param[in] dist_grid distributed grid
 * @param[in] field     field for which the mask is to be returned
 * @return field mask (NULL if no mask is available/defined)
 */
int const * yac_dist_grid_get_field_mask(
  struct dist_grid * dist_grid, struct interp_field field);

/**
 * returns a pointer to a field coordinates
 * @param[in] dist_grid distributed grid
 * @param[in] field     field for which the coordinates are to be returned
 * @return field coordinates (NULL if no mask is available)
 * @remark if the field location is CORNER and no field coordinates are
 *         available/defined, this routine will return a pointer to the
 *         vertex coordinates of the grid
 */
coordinate_pointer yac_dist_grid_get_field_coords(
  struct dist_grid * dist_grid, struct interp_field field);

/**
 * returns the data for a selected grid cell
 * @param[in]  grid_data basic grid data
 * @param[in]  cell_idx  local id of the selected cell
 * @param[out] cell      cell data
 * @remark the user has to ensure that the cell was probably initialised
 * @remark the user has to free the memory associated to cell
 *         (see \ref yac_free_grid_cell)
 */
void yac_const_basic_grid_data_get_grid_cell(
  struct const_basic_grid_data * grid_data, size_t cell_idx,
  struct grid_cell * cell);

/**
 * Determines the ranks of the distributed owners for the selected points
 * @param[in]  grid_pair distributed grid pair
 * @param[in]  grid_name grid name of the grid on which this
 *                       routine is supposed to work on
 * @param[in]  points    local ids of selected points
 * @param[in]  count     number of entries in points
 * @param[in]  location  location of the requested points
 * @param[out] ranks     distributed owners of selected points
 */
void yac_dist_grid_pair_determine_dist_owner(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  size_t * points, size_t count, enum yac_location location, int * ranks);

/**
 * Determines all non-masked neighour vertices for the selected vertices
 * @param[in]  grid_pair             distributed grid pair
 * @param[in]  grid_name             grid name of the grid on which this
 *                                   routine is supposed to work on
 * @param[in]  vertices              local ids of selected vertices
 * @param[in]  count                 number of entries in vertices
 * @param[out] neigh_vertices        neighbour vertices
 * @param[out] num_neighs_per_vertex number of neighbours per vertex
 * @param[in]  field                 if the provided field contains a mask for
 *                                   corners, it will be used
 * @remark the user has to ensure that num_neighs_per_vertex is allocated
 *         big enough
 * @remark the user is responsible for freeing memory associated to
 *         neigh_vertices
 */
void yac_dist_grid_pair_get_vertex_neighbours(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  size_t * vertices, size_t count, size_t ** neigh_vertices,
  int * num_neighs_per_vertex, struct interp_field field);

#endif // DIST_GRID_H
