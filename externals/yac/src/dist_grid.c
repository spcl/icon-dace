/**
 * @file dist_grid.c
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

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <yaxt.h>

#include "grid.h"
#include "dist_grid.h"
#include "geometry.h"
#include "yac_mpi.h"
#include "utils.h"
#include "sphere_part.h"
#include "ensure_array_size.h"
#include "interp_grid.h"
#include "yac_interface.h"

struct single_remote_point {
  yac_int global_id;
  struct remote_point_info data;
};

struct single_remote_point_reorder {
  struct single_remote_point data;
  size_t reorder_idx;
};

struct remote_point_reorder {
  struct remote_point data;
  size_t reorder_idx;
};

// reorder_idx has to be first
struct global_vertex_reorder {
  size_t reorder_idx;
  yac_int global_id;
  double coord[3];
  struct remote_point_infos owners;
};

// reorder_idx has to be first
struct global_edge_reorder {
  size_t reorder_idx;
  yac_int global_id;
  enum yac_edge_type edge_type;
  yac_int edge_to_vertex[2];
  struct remote_point_infos owners;
};

struct id_pos {
  yac_int global_id;
  uint64_t orig_pos;
};

struct n_ids_reorder {

  int count;
  yac_int * ids;
  size_t reorder_idx;
  yac_int global_id;
};

struct id_reorder_coord {
  yac_int global_id;
  size_t reorder_idx;
  double coord[3];
};

struct temp_field_data {
  int ** masks;
  size_t * masks_array_sizes;
  size_t masks_count;
  coordinate_pointer * coordinates;
  size_t * coordinates_array_sizes;
  size_t coordinates_count;
};

// warning: when changing this, ensure that struct const_basic_grid_data is
// changed accordingly
struct dist_grid {
  coordinate_pointer vertex_coordinates;
  yac_int * ids[3];
  int * num_vertices_per_cell;
  size_t * cell_to_vertex;
  size_t * cell_to_vertex_offsets;
  size_t * cell_to_edge;
  size_t * cell_to_edge_offsets;
  size_t_2_pointer edge_to_vertex;
  struct bounding_circle * cell_bnd_circles;
  enum yac_edge_type * edge_type;
  struct remote_point_infos * owners[3];   // original owner(s)
  size_t total_count[3];
  size_t count[3];
  int * core_mask[3];
  yac_int * sorted_ids[3];
  size_t * sorted_reorder_idx[3];
  struct yac_field_data_set field_data;
  MPI_Comm comm;
};

struct dist_grid_pair {
  struct dist_grid dist_grid[2];
  char const grid_names[2][YAC_MAX_CHARLEN];
  struct proc_sphere_part_node * proc_sphere_part;
  struct point_sphere_part_search * vertex_sphere_part[2];
  struct bnd_sphere_part_search * cell_sphere_part[2];
  MPI_Comm comm;
};

struct nnn_search_result {
  struct single_remote_point point_info;
  double cos_angle;
};

struct missing_edge_neighbour {
  struct {
    size_t local_id;
    yac_int global_id;
  } cell, edge;
  size_t neigh_idx;
};

// looks up positions of ids in an array of sorted ids
static void id2idx(
  yac_int * ids, size_t * idx, size_t num_ids,
  yac_int * ref_sorted_ids, size_t num_sorted_ids) {

  size_t * reorder = xmalloc(num_ids * sizeof(*reorder));
  for (size_t i = 0; i < num_ids; ++i) reorder[i] = i;

  yac_quicksort_index_yac_int_size_t(ids, num_ids, reorder);

  for (size_t i = 0, j = 0; i < num_ids; ++i) {

    yac_int curr_id = ids[i];
    while ((j < num_sorted_ids) && (ref_sorted_ids[j] < curr_id)) ++j;
    YAC_ASSERT(
      (j < num_sorted_ids) && (ref_sorted_ids[j] == curr_id),
      "ERROR(id2idx): id not found")
    idx[reorder[i]] = j;
  }

  free(reorder);
}

// returns the index of the cell vertex with the lowest global id
static size_t get_cell_reference_vertex(
  struct dist_grid * dist_grid, size_t cell_idx) {

  int num_vertices = dist_grid->num_vertices_per_cell[cell_idx];
  if (num_vertices == 0) return SIZE_MAX;
  yac_int * grid_vertex_ids = dist_grid->ids[CORNER];
  size_t * vertices =
    dist_grid->cell_to_vertex + dist_grid->cell_to_vertex_offsets[cell_idx];
  // get the cell corner with the smallest global id
  size_t min_idx = vertices[0];
  yac_int min_global_id = grid_vertex_ids[min_idx];
  for (int j = 1; j < num_vertices; ++j) {
    size_t curr_vertex = vertices[j];
    yac_int curr_global_id = grid_vertex_ids[curr_vertex];
    if (min_global_id > curr_global_id) {
      min_global_id = curr_global_id;
      min_idx = curr_vertex;
    }
  }
  return min_idx;
}

// generate cell core mask (true for all cells belonging to the local part of
// the distributed directory and the have a valid original owner)
static int * determine_cell_core_mask(
  struct dist_grid * dist_grid, int is_root, int * vertex_owner_mask) {

  size_t num_cells = dist_grid->count[CELL];
  int * core_cell_mask = xmalloc(num_cells * sizeof(*core_cell_mask));

  //-------------------------
  // determine cell core mask
  //-------------------------
  for (size_t i = 0; i < num_cells; ++i) {
    size_t ref_vertex = get_cell_reference_vertex(dist_grid, i);
    core_cell_mask[i] =
      (dist_grid->owners[CELL][i].count > 0) &&
      ((ref_vertex != SIZE_MAX)?(vertex_owner_mask[ref_vertex]):is_root);
  }
  return core_cell_mask;
}

// returns the index of the edge vertex with the lowest global id
static inline size_t get_edge_reference_vertex(
  struct dist_grid * dist_grid, size_t edge_idx) {

  yac_int * vertex_ids = dist_grid->ids[CORNER];
  size_t * edge_vertices = &(dist_grid->edge_to_vertex[edge_idx][0]);
  // get the edge corner with the smallest global id
  return edge_vertices[
          (vertex_ids[edge_vertices[0]] > vertex_ids[edge_vertices[1]])?1:0];
}

// generate edge core mask (true for all edge belonging to the local part of
// the distributed directory and the have a valid original owner)
static int * determine_edge_core_mask(
  struct dist_grid * dist_grid, int * vertex_owner_mask) {

  size_t num_edges = dist_grid->count[EDGE];
  int * core_edge_mask = xmalloc(num_edges * sizeof(*core_edge_mask));

  //-------------------------
  // determine edge core mask
  //-------------------------
  for (size_t i = 0; i < num_edges; ++i)
    core_edge_mask[i] =
      (dist_grid->owners[EDGE][i].count > 0) &&
      vertex_owner_mask[get_edge_reference_vertex(dist_grid, i)];
  return core_edge_mask;
}

// generate core masks for cells, vertices and edges
static void generate_core_masks(
  struct dist_grid * dist_grid,
  struct proc_sphere_part_node * proc_sphere_part, int comm_rank) {

  // determine distributed owner for vertices in the local part of the
  // distributed grid
  int * vertex_owner_mask =
    xmalloc(dist_grid->count[CORNER] * sizeof(*vertex_owner_mask));
  yac_proc_sphere_part_do_point_search(
    proc_sphere_part, dist_grid->vertex_coordinates,
    dist_grid->count[CORNER], vertex_owner_mask);
  for (size_t i = 0; i < dist_grid->count[CORNER]; ++i)
    vertex_owner_mask[i] = vertex_owner_mask[i] == comm_rank;

  // generate core mask for cells
  dist_grid->core_mask[CELL] =
    determine_cell_core_mask(dist_grid, comm_rank == 0, vertex_owner_mask);

  // generate core mask for edges
  dist_grid->core_mask[EDGE] =
    determine_edge_core_mask(dist_grid, vertex_owner_mask);

  // generate vertex core mask (true for vertices belonging to the local
  //                            process and have an original owner)
  for (size_t i = 0; i < dist_grid->count[CORNER]; ++i)
    vertex_owner_mask[i] =
      vertex_owner_mask[i] && (dist_grid->owners[CORNER][i].count > 0);
  dist_grid->core_mask[CORNER] = vertex_owner_mask;
}

static MPI_Datatype yac_get_single_remote_point_mpi_datatype(MPI_Comm comm) {

  struct single_remote_point dummy;
  MPI_Datatype single_id_owner_dt;
  int array_of_blocklengths[] = {1, 1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.global_id) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.data.rank) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.data.orig_pos) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {yac_int_dt, MPI_INT, MPI_UINT64_T};
  yac_mpi_call(
    MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements,
                           array_of_types, &single_id_owner_dt), comm);
  return yac_create_resized(single_id_owner_dt, sizeof(dummy), comm);
}

static MPI_Datatype yac_get_n_single_remote_point_reorder_mpi_datatype(
  int n, MPI_Comm comm) {

  MPI_Datatype dt_single_remote_point_ =
    yac_get_single_remote_point_mpi_datatype(comm);
  MPI_Datatype dt_single_remote_point__ =
    yac_create_resized(
      dt_single_remote_point_, sizeof(struct single_remote_point_reorder), comm);
  MPI_Datatype dt_single_remote_point;
  yac_mpi_call(
    MPI_Type_contiguous(
      n, dt_single_remote_point__, &dt_single_remote_point), comm);
  yac_mpi_call(MPI_Type_free(&dt_single_remote_point__), comm);
  yac_mpi_call(MPI_Type_commit(&dt_single_remote_point), comm);
  return dt_single_remote_point;
}

static int compare_single_remote_point_global_id (
  const void * a, const void * b) {

  return (((const struct single_remote_point *)a)->global_id >
          ((const struct single_remote_point *)b)->global_id) -
         (((const struct single_remote_point *)a)->global_id <
          ((const struct single_remote_point *)b)->global_id);
}

static MPI_Datatype yac_get_id_pos_mpi_datatype(MPI_Comm comm) {

  struct id_pos dummy;
  MPI_Datatype id_pos_dt;
  int array_of_blocklengths[] = {1,1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.global_id) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.orig_pos) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] = {yac_int_dt, MPI_UINT64_T};
  yac_mpi_call(
    MPI_Type_create_struct(
      2, array_of_blocklengths, array_of_displacements,
      array_of_types, &id_pos_dt), comm);
  return yac_create_resized(id_pos_dt, sizeof(dummy), comm);
}

MPI_Datatype yac_get_remote_point_info_mpi_datatype(MPI_Comm comm) {

  struct remote_point_info dummy;
  MPI_Datatype remote_point_info_dt;
  int array_of_blocklengths[] = {1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.rank) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.orig_pos) -
     (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {MPI_INT, MPI_UINT64_T};
  yac_mpi_call(
    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           array_of_types, &remote_point_info_dt), comm);
  return yac_create_resized(remote_point_info_dt, sizeof(dummy), comm);
}

static struct remote_point_infos copy_remote_point_infos(
  struct remote_point_infos point_infos) {

  int count = point_infos.count;
  if (count > 1) {
    struct remote_point_info * temp = xmalloc((size_t)count * sizeof(*temp));
    memcpy(temp, point_infos.data.multi, (size_t)count * sizeof(*temp));
    point_infos.data.multi = temp;
  }
  return point_infos;
}

static struct remote_point_infos * generate_remote_point_infos(
  yac_int * sorted_ids, size_t * reorder_idx, size_t count,
  struct remote_points * points) {

  struct remote_point_infos * point_infos = xmalloc(count * sizeof(*point_infos));

  size_t id_owner_count = points->count;
  struct remote_point * points_ = points->data;

  for (size_t i = 0, j = 0; i < count; ++i) {

    yac_int curr_id = sorted_ids[i];

    while ((j < id_owner_count) && (points_[j].global_id < curr_id)) ++j;

    if ((j >= id_owner_count) || (points_[j].global_id != curr_id))
      point_infos[reorder_idx[i]].count = -1;
    else
      point_infos[reorder_idx[i]] = copy_remote_point_infos(points_[j].data);
  }

  return point_infos;
}

static void lookup_remote_points(
  yac_int * requested_ids, struct remote_point_infos ** results,
  size_t * reorder_buffer, size_t request_count,
  struct remote_points * local_remote_data) {

  struct remote_point * local_remote_points = local_remote_data->data;
  size_t local_count = local_remote_data->count;
  static struct remote_point_infos dummy_data =
    {.count = 0, .data.single = {.rank = INT_MAX, .orig_pos = UINT64_MAX}};

  for (size_t j = 0; j < request_count; ++j) reorder_buffer[j] = j;

  // sort ids for easier owner lookup
  yac_quicksort_index_yac_int_size_t(
    requested_ids, request_count, reorder_buffer);

  for (size_t i = 0, j = 0; i < request_count; ++i) {
    yac_int curr_id = requested_ids[i];
    while ((j < local_count) && (curr_id > local_remote_points[j].global_id))
      ++j;
    results[reorder_buffer[i]] =
      ((j < local_count) && (curr_id == local_remote_points[j].global_id))?
        &(local_remote_points[j].data):&dummy_data;
  }
}

static int yac_remote_point_infos_get_pack_size(
  struct remote_point_infos const * infos, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  int pack_size_count,
      pack_size_data;

  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &pack_size_count), comm);
  yac_mpi_call(
    MPI_Pack_size(infos->count, point_info_dt, comm, &pack_size_data), comm);

  return pack_size_count + pack_size_data;
}

static void yac_remote_point_infos_get_pack_sizes(
  struct remote_point_infos ** infos, size_t count, int * pack_sizes,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  int pack_size_count;

  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &pack_size_count), comm);

  for (size_t i = 0; i < count; ++i) {
    yac_mpi_call(
      MPI_Pack_size(
        infos[i]->count, point_info_dt, comm, pack_sizes + i), comm);
    pack_sizes[i] += pack_size_count;
  }
}

static void yac_remote_point_infos_pack(
  struct remote_point_infos ** infos, int * pack_sizes, size_t count,
  void * buffer, MPI_Datatype point_info_dt, MPI_Comm comm) {

  for (size_t i = 0; i < count; ++i) {

    int curr_count = infos[i]->count;
    struct remote_point_info * curr_point_infos =
      (curr_count == 1)?
        (&(infos[i]->data.single)):(infos[i]->data.multi);

    int position = 0;
    int buffer_size = pack_sizes[i];
    yac_mpi_call(
        MPI_Pack(&curr_count, 1, MPI_INT, buffer,
                 buffer_size, &position, comm), comm);
    yac_mpi_call(
        MPI_Pack(curr_point_infos, curr_count, point_info_dt, buffer,
                 buffer_size, &position, comm), comm);
    buffer = (void*)(((unsigned char*)buffer) + buffer_size);
  }
}

static void yac_remote_point_infos_unpack(
  void * buffer, int buffer_size, int * position,
  struct remote_point_infos * infos, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int count;
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, &count, 1, MPI_INT, comm), comm);

  infos->count = count;

  if (count > 0) {
    struct remote_point_info * point_infos;
    if (count == 1)
      point_infos = &(infos->data.single);
    else
      point_infos =
        (infos->data.multi = xmalloc((size_t)count * sizeof(*point_infos)));

    yac_mpi_call(
      MPI_Unpack(buffer, buffer_size, position, point_infos, count,
                 point_info_dt, comm), comm);
  } else {
    infos->data.single.rank = INT_MAX;
    infos->data.single.orig_pos = UINT64_MAX;
  }
}

static void yac_remote_point_infos_unpack_info_buffer(
  void * buffer, int buffer_size, int * position,
  struct remote_point_info * info_buffer, size_t * info_buffer_position,
  struct remote_point_infos * infos, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int count;
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, &count, 1, MPI_INT, comm), comm);

  infos->count = count;

  struct remote_point_info * point_infos;
  if (count == 1)
    point_infos = &(infos->data.single);
  else {
    point_infos =
      (infos->data.multi = info_buffer + *info_buffer_position);
    *info_buffer_position += count;
  }

  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, point_infos, count,
               point_info_dt, comm), comm);
}

int yac_remote_point_get_pack_size(
  struct remote_point * point, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int pack_size_id,
      pack_size_remote_point_infos;

  yac_mpi_call(MPI_Pack_size(1, yac_int_dt, comm, &pack_size_id), comm);
  pack_size_remote_point_infos =
    yac_remote_point_infos_get_pack_size(&(point->data), point_info_dt, comm);

  return pack_size_id + pack_size_remote_point_infos;
}

static void yac_remote_point_infos_pack_single(
  struct remote_point_infos const * infos, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int count = infos->count;

  struct remote_point_info const * info =
    (count == 1)?(&(infos->data.single)):(infos->data.multi);

  yac_mpi_call(
      MPI_Pack(&count, 1, MPI_INT, buffer, buffer_size, position, comm), comm);
  yac_mpi_call(
      MPI_Pack(info, count, point_info_dt, buffer, buffer_size, position, comm),
      comm);
}

void yac_remote_point_pack(
  struct remote_point * point, void * buffer, int buffer_size, int * position,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  yac_mpi_call(
      MPI_Pack(&(point->global_id), 1, yac_int_dt, buffer,
               buffer_size, position, comm), comm);

  yac_remote_point_infos_pack_single(
    &(point->data), buffer, buffer_size, position, point_info_dt, comm);
}

void yac_remote_point_unpack(
  void * buffer, int buffer_size, int * position, struct remote_point * point,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &(point->global_id), 1, yac_int_dt, comm),
    comm);
  yac_remote_point_infos_unpack(
    buffer, buffer_size, position, &(point->data), point_info_dt, comm);
}

void yac_remote_point_unpack_info_buffer(
  void * buffer, int buffer_size, int * position,
  struct remote_point_info * info_buffer, size_t * info_buffer_position,
  struct remote_point * point, MPI_Datatype point_info_dt, MPI_Comm comm) {

  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &(point->global_id), 1, yac_int_dt, comm),
    comm);
  yac_remote_point_infos_unpack_info_buffer(
    buffer, buffer_size, position, info_buffer, info_buffer_position,
    &(point->data), point_info_dt, comm);
}

int yac_remote_points_get_pack_size(
  struct remote_points * points, MPI_Datatype point_info_dt, MPI_Comm comm) {

  size_t count = points->count;
  struct remote_point * points_data = points->data;

  int count_pack_size,
      remote_points_pack_size;

  yac_mpi_call(MPI_Pack_size(2, MPI_UINT64_T, comm, &count_pack_size), comm);
  remote_points_pack_size = 0;
  for (size_t i = 0; i < count; ++i)
    remote_points_pack_size +=
      yac_remote_point_get_pack_size(points_data + i, point_info_dt, comm);

  return count_pack_size + remote_points_pack_size;
}

void yac_remote_points_pack(
  struct remote_points * points, void * buffer, int buffer_size, int * position,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  size_t count = points->count;
  struct remote_point * point_data = points->data;
  uint64_t counts[2] = {(uint64_t)count, 0};
  for (size_t i = 0; i < count; ++i)
    if (point_data[i].data.count > 1)
      counts[1] += (uint64_t)(point_data[i].data.count);

  yac_mpi_call(
      MPI_Pack(counts, 2, MPI_UINT64_T, buffer,
               buffer_size, position, comm), comm);
  for (size_t i = 0; i < count; ++i)
    yac_remote_point_pack(
      point_data + i, buffer, buffer_size, position, point_info_dt, comm);
}

void yac_remote_points_unpack(
  void * buffer, int buffer_size, int * position,
  struct remote_points ** points, MPI_Datatype point_info_dt, MPI_Comm comm) {

  uint64_t counts[2];

  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, counts, 2, MPI_UINT64_T, comm), comm);

  *points = xmalloc(((size_t)counts[1]) * sizeof(struct remote_point_infos) +
                    sizeof(**points));

  size_t count = ((*points)->count = (size_t)(counts[0]));
  struct remote_point * point_data =
    ((*points)->data = xmalloc(count * sizeof(*((*points)->data))));
  struct remote_point_info * remote_point_info_buffer =
    &((*points)->buffer[0]);

  for (size_t i = 0, offset = 0; i < count; ++i) {

    yac_remote_point_unpack_info_buffer(
      buffer, buffer_size, position, remote_point_info_buffer, &offset,
      point_data + i, point_info_dt, comm);
  }
}

static void generate_dist_grid_remote_point_infos(
  struct dist_grid * dist_grid,
  struct proc_sphere_part_node * proc_sphere_part,
  struct remote_points ** dist_owner) {

  MPI_Comm comm = dist_grid->comm;

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  for (int location = 0; location < 3; ++location)
    dist_grid->owners[location] =
      generate_remote_point_infos(
        dist_grid->sorted_ids[location],
        dist_grid->sorted_reorder_idx[location],
        dist_grid->count[location], dist_owner[location]);

  size_t num_missing_vertices = 0;
  size_t num_missing_edges = 0;

  for (size_t i = 0; i < dist_grid->count[CELL]; ++i)
    YAC_ASSERT(
      dist_grid->owners[CELL][i].count != -1,
      "ERROR(generate_dist_grid_remote_point_infos):"
      "no owner information for a cell available")
  for (size_t i = 0; i < dist_grid->count[CORNER]; ++i)
    if (dist_grid->owners[CORNER][i].count == -1) ++num_missing_vertices;
  for (size_t i = 0; i < dist_grid->count[EDGE]; ++i)
    if (dist_grid->owners[EDGE][i].count == -1) ++num_missing_edges;

  size_t total_num_missing = num_missing_vertices + num_missing_edges;
  int * ranks_buffer = xmalloc(total_num_missing * sizeof(*ranks_buffer));
  int * vertex_ranks = ranks_buffer;
  int * edge_ranks = ranks_buffer + num_missing_vertices;
  size_t * size_t_buffer =
    xmalloc((3 * total_num_missing +
             4 * (size_t)comm_size) * sizeof(*size_t_buffer));
  size_t * vertex_reorder = size_t_buffer;
  size_t * edge_reorder = size_t_buffer + num_missing_vertices;
  size_t * coord_indices = size_t_buffer + total_num_missing;
  size_t * ranks_idxs = coord_indices;
  size_t * cve_reorder = size_t_buffer + 2 * total_num_missing;
  size_t * total_sendcounts = size_t_buffer + 3 * total_num_missing + 0 * comm_size;
  size_t * total_recvcounts = size_t_buffer + 3 * total_num_missing + 1 * comm_size;
  size_t * total_sdispls =    size_t_buffer + 3 * total_num_missing + 2 * comm_size;
  size_t * total_rdispls =    size_t_buffer + 3 * total_num_missing + 3 * comm_size;

  size_t offset = 0;
  for (size_t i = 0, j = 0; i < dist_grid->count[CORNER]; ++i) {
    if (dist_grid->owners[CORNER][i].count == -1) {

      vertex_reorder[j] = i;
      ++j;
      coord_indices[offset] = i;
      cve_reorder[offset] = offset;
      ++offset;
    }
  }
  for (size_t i = 0, j = 0; i < dist_grid->count[EDGE]; ++i) {
    if (dist_grid->owners[EDGE][i].count == -1) {
      // compute middle point of edge
      size_t * vertex_idxes = &(dist_grid->edge_to_vertex[i][0]);
      yac_int vertex_ids[2] = {dist_grid->ids[CORNER][vertex_idxes[0]],
                               dist_grid->ids[CORNER][vertex_idxes[1]]};
      edge_reorder[j] = i;
      ++j;
      coord_indices[offset] =
        vertex_idxes[(vertex_ids[0] > vertex_ids[1])?1:0];
      cve_reorder[offset] = offset;
      ++offset;
    }
  }

  yac_quicksort_index_size_t_size_t(
    coord_indices, total_num_missing, cve_reorder);

  size_t num_unique_idxs = 0;
  for (size_t i = 0, prev_idx = SIZE_MAX; i < total_num_missing; ++i) {
    size_t curr_idx = coord_indices[i];
    if (curr_idx != prev_idx) {
      prev_idx = curr_idx;
      ++num_unique_idxs;
    }
  }

  coordinate_pointer search_coords =
    xmalloc(num_unique_idxs * sizeof(*search_coords));
  int * search_ranks = xmalloc(num_unique_idxs * sizeof(*search_ranks));

  for (size_t i = 0, j = 0, prev_idx = SIZE_MAX; i < total_num_missing; ++i) {
    size_t curr_idx = coord_indices[i];
    if (curr_idx != prev_idx) {
      memcpy(search_coords[j], dist_grid->vertex_coordinates[curr_idx],
             3 * sizeof(double));
      prev_idx = curr_idx;
      ++j;
    }
    ranks_idxs[i] = j-1;
  }

  yac_proc_sphere_part_do_point_search(
    proc_sphere_part, search_coords, num_unique_idxs, search_ranks);
  free(search_coords);

  yac_quicksort_index_size_t_size_t(
    cve_reorder, total_num_missing, ranks_idxs);

  offset = 0;
  for (size_t i = 0; i < num_missing_vertices; ++i, ++offset)
      vertex_ranks[i] = search_ranks[ranks_idxs[offset]];
  for (size_t i = 0; i < num_missing_edges; ++i, ++offset)
      edge_ranks[i] = search_ranks[ranks_idxs[offset]];
  free(search_ranks);

  yac_quicksort_index_int_size_t(
    vertex_ranks, num_missing_vertices, vertex_reorder);
  yac_quicksort_index_int_size_t(
    edge_ranks, num_missing_edges, edge_reorder);

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    2, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0; i < num_missing_vertices; ++i)
      sendcounts[2 * vertex_ranks[i] + 0]++;
  for (size_t i = 0; i < num_missing_edges; ++i)
      sendcounts[2 * edge_ranks[i] + 1]++;

  yac_generate_alltoallv_args(
    2, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t num_requested_ids[2] = {0,0};
  size_t saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {
    total_sdispls[i] = saccu;
    total_rdispls[i] = raccu;
    total_sendcounts[i] = 0;
    total_recvcounts[i] = 0;
    for (int j = 0; j < 2; ++j) {
      total_sendcounts[i] += sendcounts[2 * i + j];
      total_recvcounts[i] += recvcounts[2 * i + j];
      num_requested_ids[j] += recvcounts[2 * i + j];
    }
    saccu += total_sendcounts[i];
    raccu += total_recvcounts[i];
  }
  size_t request_count = total_recvcounts[comm_size - 1] +
                         total_rdispls[comm_size - 1];

  yac_int * exchange_id_buffer =
    xmalloc(
      (total_num_missing + 2 * request_count) * sizeof(*exchange_id_buffer));
  yac_int * id_send_buffer = exchange_id_buffer + request_count;
  yac_int * id_recv_buffer = exchange_id_buffer + request_count + total_num_missing;

  // pack the ids for which the local process is missing the owner information
  for (size_t i = 0; i < num_missing_vertices; ++i) {
    size_t pos = sdispls[2 * vertex_ranks[i] + 0 + 1]++;
    id_send_buffer[pos] = dist_grid->ids[CORNER][vertex_reorder[i]];
  }
  for (size_t i = 0; i < num_missing_edges; ++i) {
    size_t pos = sdispls[2 * edge_ranks[i] + 1 + 1]++;
    id_send_buffer[pos] = dist_grid->ids[EDGE][edge_reorder[i]];
  }

  // exchange the ids for which the local process is missing the owner information
  yac_alltoallv_yac_int_p2p(
    id_send_buffer, total_sendcounts, total_sdispls,
    id_recv_buffer, total_recvcounts, total_rdispls, comm);

  yac_int * requested_ids[2] =
    {exchange_id_buffer, exchange_id_buffer + num_requested_ids[0]};

  // unpack data
  {
    size_t ids_offset[2] = {0, 0}, offset = 0;
    for (int i = 0; i < comm_size; ++i) {
      for (int j = 0; j < 2; ++j) {
        size_t curr_count = recvcounts[2 * i + j];
        memcpy(requested_ids[j] + ids_offset[j], id_recv_buffer + offset,
               curr_count * sizeof(*(requested_ids[0])));
        ids_offset[j] += curr_count;
        offset += curr_count;
      }
    }
  }

  struct remote_point_infos ** remote_point_infos_buffer =
    xmalloc(request_count * sizeof(*remote_point_infos_buffer));
  struct remote_point_infos ** found_points[2] =
    {remote_point_infos_buffer,
     remote_point_infos_buffer + num_requested_ids[0]};
  struct remote_points * local_remote_points[2] =
    {dist_owner[CORNER], dist_owner[EDGE]};

  size_t max_num_requested_ids = MAX(num_requested_ids[0], num_requested_ids[1]);
  size_t * reorder_buffer = xmalloc(max_num_requested_ids * sizeof(*reorder_buffer));
  // lookup owners for requested ids
  for (int i = 0; i < 2; ++i)
    lookup_remote_points(
      requested_ids[i], found_points[i], reorder_buffer, num_requested_ids[i],
      local_remote_points[i]);
  free(reorder_buffer);
  free(exchange_id_buffer);

  int * pack_size_buffer =
    xmalloc((total_num_missing + request_count) * sizeof(*pack_size_buffer));
  int * send_pack_sizes = pack_size_buffer;
  int * recv_pack_sizes = pack_size_buffer + request_count;
  MPI_Datatype point_info_dt = yac_get_remote_point_info_mpi_datatype(comm);

  {
    size_t offset = 0, offsets[2] = {0, 0};
    for (int i = 0; i < comm_size; ++i) {
      for (int j = 0; j < 2; ++j) {
        size_t curr_count = (size_t)(recvcounts[2 * i + j]);
        // get pack sizes
        yac_remote_point_infos_get_pack_sizes(
          found_points[j] + offsets[j], curr_count,
          send_pack_sizes + offset, point_info_dt, comm);
        offsets[j] += curr_count;
        offset += curr_count;
      }
    }
  }

  // exchange pack sizes
  yac_alltoallv_int_p2p(
    send_pack_sizes, total_recvcounts, total_rdispls,
    recv_pack_sizes, total_sendcounts, total_sdispls, comm);

  size_t send_buffer_size = 0, recv_buffer_size = 0;
  for (size_t i = 0; i < request_count; ++i)
    send_buffer_size += (size_t)(send_pack_sizes[i]);
  for (size_t i = 0; i < total_num_missing; ++i)
    recv_buffer_size += (size_t)(recv_pack_sizes[i]);

  void * result_buffer = xmalloc(send_buffer_size + recv_buffer_size);
  void * send_result_buffer = result_buffer;
  void * recv_result_buffer = (unsigned char *)result_buffer + send_buffer_size;

  // pack data
  {
    size_t send_buffer_size = 0, offset = 0, offsets[2] = {0,0};
    for (int i = 0; i < comm_size; ++i) {
      for (int j = 0; j < 2; ++j) {
        size_t curr_count = recvcounts[2 * i + j];
        // get pack sizes
        yac_remote_point_infos_pack(
          found_points[j] + offsets[j],
          send_pack_sizes + offset, curr_count,
          (void*)((unsigned char*)send_result_buffer + send_buffer_size),
          point_info_dt, comm);
        for (size_t k = 0; k < curr_count; ++k)
          send_buffer_size += send_pack_sizes[offset + k];
        offsets[j] += curr_count;
        offset += curr_count;
      }
    }
  }

  {
    size_t send_offset = 0, recv_offset = 0;
    size_t saccu = 0, raccu = 0;
    for (int i = 0; i < comm_size; ++i) {
      total_sdispls[i] = saccu;
      total_rdispls[i] = raccu;
      size_t curr_sendcount = 0, curr_recvcount = 0;
      for (int k = 0; k < 2; ++k) {
        for (size_t j = 0; j < recvcounts[2 * i + k]; ++j, ++send_offset)
          curr_sendcount += (size_t)(send_pack_sizes[send_offset]);
        for (size_t j = 0; j < sendcounts[2 * i + k]; ++j, ++recv_offset)
          curr_recvcount += (size_t)(recv_pack_sizes[recv_offset]);
      }
      total_sendcounts[i] = curr_sendcount;
      total_recvcounts[i] = curr_recvcount;
      saccu += curr_sendcount;
      raccu += curr_recvcount;
    }
  }

  // exchange results
  yac_alltoallv_packed_p2p(
    send_result_buffer, total_sendcounts, total_sdispls,
    recv_result_buffer, total_recvcounts, total_rdispls, comm);

  // unpack results
  {
    size_t counts[2] = {0, 0}, count = 0, offset = 0;
    struct remote_point_infos * remote_points[2] =
      {dist_grid->owners[CORNER], dist_grid->owners[EDGE]};
    size_t * reorder[2] = {vertex_reorder, edge_reorder};
    for (int i = 0; i < comm_size; ++i) {
      for (int j = 0; j < 2; ++j) {
        size_t curr_count = sendcounts[2 * i + j];
        for (size_t k = 0; k < curr_count; ++k, ++counts[j], ++count) {
          int curr_buffer_size = recv_pack_sizes[sdispls[2 * i + j] + k];
          int position = 0;
          yac_remote_point_infos_unpack(
            (void*)((unsigned char*)recv_result_buffer + offset),
            curr_buffer_size, &position,
            remote_points[j] + reorder[j][counts[j]], point_info_dt, comm);
          offset += (size_t)(curr_buffer_size);
        }
      }
    }
  }

  yac_mpi_call(MPI_Type_free(&point_info_dt), comm);

  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  free(result_buffer);
  free(pack_size_buffer);
  free(remote_point_infos_buffer);
  free(size_t_buffer);
  free(ranks_buffer);
}

static int compare_remote_point(const void * a, const void * b) {

  return ((const struct remote_point*)a)->global_id -
         ((const struct remote_point*)b)->global_id;
}

// inserts an element into an array and increases the corresponding size
static void insert_global_id(yac_int * ids, size_t n, yac_int id) {

  size_t i;
  for (i = 0; i < n; ++i) if (ids[i] >= id) break;
  // copy new id into array and move bigger elements one position up
  if (n != i) memmove(ids + i + 1, ids + i, (n - i) * sizeof(*ids));
  ids[i] = id;
}

// inserts an element into an array and increases the corresponding size
// if the element already exists in the array nothing is done
static void insert_rank(int * ranks, int * count, int rank) {

  int i;
  int n = *count;

  for (i = 0; i < n; ++i) if (ranks[i] >= rank) break;

  // if the rank is already in the array
  if (i != n) {
    if (ranks[i] == rank) return;
    else memmove(ranks + i + 1, ranks + i, ((size_t)(n - i)) * sizeof(*ranks));
  }
  ranks[i] = rank;
  *count = n + 1;
}

static int compare_n_ids_reorder_ids (
  const void * a, const void * b) {

  int count_a = ((struct n_ids_reorder const *)a)->count;
  int count_b = ((struct n_ids_reorder const *)b)->count;
  yac_int * a_ids = ((struct n_ids_reorder const *)a)->ids;
  yac_int * b_ids = ((struct n_ids_reorder const *)b)->ids;
  int ret = count_a - count_b;
  for (int i = 0; !ret && (i < count_a); ++i)
    ret = (a_ids[i] > b_ids[i]) - (a_ids[i] < b_ids[i]);
  return ret;
}

static int compare_n_ids_reorder_reorder (
  const void * a, const void * b) {

  return (((struct n_ids_reorder const *)a)->reorder_idx >
          ((struct n_ids_reorder const *)b)->reorder_idx) -
         (((struct n_ids_reorder const *)a)->reorder_idx <
          ((struct n_ids_reorder const *)b)->reorder_idx);
}

// determines for all cells, in the grid data provided by the user, to
// which processes they belong according to the decomposition of the
// distributed grid
static void determine_dist_cell_ranks(
  struct proc_sphere_part_node * proc_sphere_part,
  struct basic_grid_data * grid, MPI_Comm comm, int ** dist_cell_ranks,
  int * dist_cell_rank_counts, size_t * dist_cell_rank_offsets) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t num_cells = grid->num_cells;
  int * ranks_buffer = xmalloc((size_t)comm_size * sizeof(*ranks_buffer));
  size_t dist_cell_ranks_array_size = num_cells;
  int * dist_cell_ranks_ = xmalloc(num_cells * sizeof(*dist_cell_ranks_));
  size_t offset = 0;

  // determine maximum number of vertices per cell
  int max_num_vertices_per_cell = 0;
  for (size_t i = 0; i < num_cells; ++i)
    if (grid->num_vertices_per_cell[i] > max_num_vertices_per_cell)
      max_num_vertices_per_cell = grid->num_vertices_per_cell[i];

  struct grid_cell cell;
  cell.coordinates_xyz =
    xmalloc(
      (size_t)max_num_vertices_per_cell * sizeof(*(cell.coordinates_xyz)));
  cell.edge_type =
    xmalloc((size_t)max_num_vertices_per_cell * sizeof(*(cell.edge_type)));

  size_t * cell_to_vertex = grid->cell_to_vertex;
  size_t * cell_to_vertex_offsets = grid->cell_to_vertex_offsets;
  size_t * cell_to_edge = grid->cell_to_edge;
  size_t * cell_to_edge_offsets = grid->cell_to_edge_offsets;
  coordinate_pointer vertex_coordinates = grid->vertex_coordinates;
  enum yac_edge_type * edge_type = grid->edge_type;
  int * num_vertices_per_cell = grid->num_vertices_per_cell;

  // generate a bounding circle for each cell and use it to determine
  // ranks of the processes that require this cell
  for (size_t i = 0; i < num_cells; ++i) {

    // get bounding circle for current cell
    cell.num_corners = num_vertices_per_cell[i];
    size_t * curr_cell_to_vertex =
      cell_to_vertex + cell_to_vertex_offsets[i];
    size_t * curr_cell_to_edge =
      cell_to_edge + cell_to_edge_offsets[i];
    for (int j = 0; j < num_vertices_per_cell[i]; ++j) {
      coordinate_pointer curr_vertex_coords =
        vertex_coordinates + curr_cell_to_vertex[j];
      for (int k = 0; k < 3; ++k)
        cell.coordinates_xyz[j][k] = (*curr_vertex_coords)[k];
      cell.edge_type[j] = edge_type[curr_cell_to_edge[j]];
    }
    int rank_count;
    if (cell.num_corners > 0) {

      struct bounding_circle bnd_circle;
      yac_get_cell_bounding_circle(cell, &bnd_circle);

      yac_proc_sphere_part_do_bnd_circle_search(
        proc_sphere_part, bnd_circle, ranks_buffer, &rank_count);
    } else {
      ranks_buffer[0] = 0;
      rank_count = 1;
    }

    ENSURE_ARRAY_SIZE(dist_cell_ranks_, dist_cell_ranks_array_size,
                      offset + (size_t)rank_count);
    memcpy(dist_cell_ranks_ + offset, ranks_buffer,
           (size_t)rank_count * sizeof(*ranks_buffer));
    dist_cell_rank_counts[i] = rank_count;
    dist_cell_rank_offsets[i] = offset;
    offset += (size_t)rank_count;
  }

  dist_cell_rank_offsets[num_cells] = offset;

  free(cell.edge_type);
  free(cell.coordinates_xyz);
  free(ranks_buffer);

  *dist_cell_ranks = dist_cell_ranks_;
}

static MPI_Datatype yac_get_id_reorder_coord_coord_mpi_datatype(MPI_Comm comm) {

  struct id_reorder_coord dummy;
  MPI_Datatype coord_dt;
  int array_of_blocklengths[] = {3};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.coord) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] = {MPI_DOUBLE};
  yac_mpi_call(
    MPI_Type_create_struct(1, array_of_blocklengths, array_of_displacements,
                           array_of_types, &coord_dt), comm);
  return yac_create_resized(coord_dt, sizeof(dummy), comm);
}

static MPI_Datatype yac_get_id_reorder_coord_id_mpi_datatype(MPI_Comm comm) {

  struct id_reorder_coord dummy;
  MPI_Datatype global_id_dt;
  int array_of_blocklengths[] = {1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.global_id) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] = {yac_int_dt};
  yac_mpi_call(
    MPI_Type_create_struct(1, array_of_blocklengths, array_of_displacements,
                           array_of_types, &global_id_dt), comm);
  return yac_create_resized(global_id_dt, sizeof(dummy), comm);
}

static int compare_id_reorder_coord_coord (
  const void * a, const void * b) {

  return compare_coords(((const struct id_reorder_coord *)a)->coord,
                        ((const struct id_reorder_coord *)b)->coord);
}

static int compare_id_reorder_coord_reorder_idx (
  const void * a, const void * b) {

  return (((const struct id_reorder_coord *)a)->reorder_idx >
          ((const struct id_reorder_coord *)b)->reorder_idx) -
         (((const struct id_reorder_coord *)a)->reorder_idx <
          ((const struct id_reorder_coord *)b)->reorder_idx);
}

// generates global vertex ids if they are not provied by the user
// and returns the distributed owners for each vertex in the provided
// grid data
static int * generate_vertex_ids(
  struct proc_sphere_part_node * proc_sphere_part,
  struct basic_grid_data * grid, MPI_Comm comm) {

  int * vertex_ranks = xmalloc(grid->num_vertices * sizeof(*vertex_ranks));

  // determine the dist grid ranks for all vertices
  yac_proc_sphere_part_do_point_search(
    proc_sphere_part, grid->vertex_coordinates, grid->num_vertices,
    vertex_ranks);

  // check whether only of subset of the processes have defined
  // their global ids, which is not supported
  int ids_available_local =
    (grid->num_vertices > 0) && (grid->vertex_ids != NULL);
  int ids_available_global;
  yac_mpi_call(MPI_Allreduce(&ids_available_local, &ids_available_global, 1,
                             MPI_INT, MPI_MAX, comm), comm);

  // if there are vertices defined locally for the current grid
  if (grid->num_vertices != 0) {
    YAC_ASSERT(
      ids_available_local || !ids_available_global,
      "ERROR(generate_vertex_ids): inconsistent global ids")

    // generate core mask if not available
    if (grid->core_vertex_mask == NULL) {
      int * core_mask =
        (grid->core_vertex_mask =
           xmalloc(grid->num_vertices * sizeof(*core_mask)));
      for (size_t j = 0; j < grid->num_vertices; ++j) core_mask[j] = 1;
    }
  }

  // if we do not need to generate the global ids
  if (ids_available_global) return vertex_ranks;

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0; i < grid->num_vertices; ++i)
    sendcounts[vertex_ranks[i]]++;

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t orig_vertex_count = grid->num_vertices;
  size_t dist_vertex_count =
    recvcounts[comm_size - 1] + rdispls[comm_size - 1];

  struct id_reorder_coord * id_reorder_coord_buffer =
    xmalloc((orig_vertex_count + dist_vertex_count) *
            sizeof(*id_reorder_coord_buffer));
  struct id_reorder_coord *
    id_reorder_coord_send_buffer = id_reorder_coord_buffer;
  struct id_reorder_coord *
    id_reorder_coord_recv_buffer = id_reorder_coord_buffer + orig_vertex_count;

  // pack send buffer
  for (size_t i = 0; i < orig_vertex_count; ++i) {
    size_t pos = sdispls[vertex_ranks[i] + 1]++;
    id_reorder_coord_send_buffer[pos].reorder_idx = i;
    for (int j = 0; j < 3; ++j)
      id_reorder_coord_send_buffer[pos].coord[j] =
        grid->vertex_coordinates[i][j];
  }

  MPI_Datatype id_reorder_coord_coord_dt =
    yac_get_id_reorder_coord_coord_mpi_datatype(comm);

  // exchange data
  yac_alltoallv_p2p(
    id_reorder_coord_send_buffer, sendcounts, sdispls,
    id_reorder_coord_recv_buffer, recvcounts, rdispls,
    sizeof(*id_reorder_coord_send_buffer), id_reorder_coord_coord_dt, comm);

  yac_mpi_call(MPI_Type_free(&id_reorder_coord_coord_dt), comm);

  for (size_t i = 0; i < dist_vertex_count; ++i)
    id_reorder_coord_recv_buffer[i].reorder_idx = i;

  // sort received vertices based on coordinates
  qsort(id_reorder_coord_recv_buffer, dist_vertex_count,
        sizeof(*id_reorder_coord_recv_buffer), compare_id_reorder_coord_coord);

  size_t unique_count = dist_vertex_count > 0;
  struct id_reorder_coord * prev = id_reorder_coord_recv_buffer;

  // determine number of unique coordinates
  for (size_t i = 0; i < dist_vertex_count; ++i) {
    struct id_reorder_coord * curr = id_reorder_coord_recv_buffer + i;
    if (compare_id_reorder_coord_coord(prev, curr)) {
      ++unique_count;
      prev = curr;
    }
    curr->global_id = (yac_int)unique_count - 1;
  }

  YAC_ASSERT(
    unique_count <= (size_t)XT_INT_MAX,
    "ERROR(generate_vertex_ids): global_id out of bounds")

  yac_int yac_int_unique_count = (yac_int)unique_count;
  yac_int id_offset;

  // determine exclusive scan of sum of numbers of unique
  // coordinates on all ranks
  yac_mpi_call(MPI_Exscan(&yac_int_unique_count, &id_offset, 1, yac_int_dt,
                          MPI_SUM, comm), comm);
  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  if (comm_rank == 0) id_offset = 0;

  YAC_ASSERT(
    ((size_t)id_offset + unique_count) <= (size_t)XT_INT_MAX,
    "ERROR(generate_vertex_ids): global_id out of bounds")

  // adjust global ids
  for (size_t i = 0; i < dist_vertex_count; ++i)
    id_reorder_coord_recv_buffer[i].global_id += id_offset;

  // return received vertices into original order
  qsort(id_reorder_coord_recv_buffer, dist_vertex_count,
        sizeof(*id_reorder_coord_recv_buffer),
        compare_id_reorder_coord_reorder_idx);

  MPI_Datatype id_reorder_coord_id_dt =
    yac_get_id_reorder_coord_id_mpi_datatype(comm);

  // return generated global ids data
  yac_alltoallv_p2p(
    id_reorder_coord_recv_buffer, recvcounts, rdispls,
    id_reorder_coord_send_buffer, sendcounts, sdispls,
    sizeof(*id_reorder_coord_send_buffer), id_reorder_coord_id_dt, comm);

  yac_mpi_call(MPI_Type_free(&id_reorder_coord_id_dt), comm);

  yac_int * vertex_ids =
    ((grid->vertex_ids =
      (grid->num_vertices > 0)?
        xmalloc(grid->num_vertices * sizeof(*(grid->vertex_ids))):NULL));

  for (size_t i = 0; i < grid->num_vertices; ++i)
    vertex_ids[id_reorder_coord_send_buffer[i].reorder_idx] =
      id_reorder_coord_send_buffer[i].global_id;

  free(id_reorder_coord_buffer);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  return vertex_ranks;
}

// generates global cell/edge ids if they are not provied by the user
static void generate_ce_ids(
  struct proc_sphere_part_node * proc_sphere_part,
  struct basic_grid_data * grid, int * vertex_ranks, MPI_Comm comm) {

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // check whether only of subset of the processes have defined
  // their global ids, which is not supported
  int ids_available_local[2], ids_available_global[2];
  ids_available_local[0] = (grid->num_cells > 0) && (grid->cell_ids != NULL);
  ids_available_local[1] = (grid->num_edges > 0) && (grid->edge_ids != NULL);
  yac_mpi_call(MPI_Allreduce(ids_available_local, ids_available_global, 2,
                             MPI_INT, MPI_MAX, comm), comm);

  if (grid->num_cells > 0) {
    YAC_ASSERT(
      ids_available_local[0] == ids_available_global[0],
      "ERROR(generate_ce_ids): inconsistent global ids")
    // generate core mask if not available
    if (grid->core_cell_mask == NULL) {
      int * core_cell_mask =
        (grid->core_cell_mask =
           xmalloc(grid->num_cells * sizeof(*core_cell_mask)));
      for (size_t j = 0; j < grid->num_cells; ++j) core_cell_mask[j] = 1;
    }
  }
  if (grid->num_edges > 0) {
    YAC_ASSERT(
      ids_available_local[1] == ids_available_global[1],
      "ERROR(generate_ce_ids): inconsistent global ids")
    // generate core mask if not available
    if (grid->core_edge_mask == NULL) {
      int * core_edge_mask =
        (grid->core_edge_mask =
           xmalloc(grid->num_edges * sizeof(*core_edge_mask)));
      for (size_t j = 0; j < grid->num_edges; ++j) core_edge_mask[j] = 1;
    }
  }

  // if no ids need to be generated
  if (ids_available_global[0] && ids_available_global[1]) return;

  int * rank_buffer =
    xmalloc(
      (((ids_available_global[0])?0:(grid->num_cells)) +
       ((ids_available_global[1])?0:(grid->num_edges))) * sizeof(*rank_buffer));
  int * cell_ranks = rank_buffer;
  int * edge_ranks =
    rank_buffer + ((ids_available_global[0])?0:(grid->num_cells));

  size_t * size_t_buffer =
    xmalloc((8 * (size_t)comm_size + 1) * sizeof(*size_t_buffer));
  size_t * sendcounts = size_t_buffer + 0 * comm_size;
  size_t * recvcounts = size_t_buffer + 2 * comm_size;
  size_t * total_sendcounts = size_t_buffer + 4 * comm_size;
  size_t * total_recvcounts = size_t_buffer + 5 * comm_size;
  size_t * total_sdispls = size_t_buffer + 6 * comm_size;
  size_t * total_rdispls = size_t_buffer + 7 * comm_size + 1;
  memset(sendcounts, 0, 2 * (size_t)comm_size * sizeof(*sendcounts));

  yac_int * cell_to_vertex_ids = NULL;
  int max_num_vertices_per_cell = 0;

  if (!ids_available_global[0]) {

    { // compute the global maximum number of vertices per cell
      for (size_t i = 0; i < grid->num_cells; ++i)
        if (grid->num_vertices_per_cell[i] > max_num_vertices_per_cell)
          max_num_vertices_per_cell = grid->num_vertices_per_cell[i];
      yac_mpi_call(MPI_Allreduce(MPI_IN_PLACE, &max_num_vertices_per_cell, 1,
                                 MPI_INT, MPI_MAX, comm), comm);
    }

    cell_to_vertex_ids =
      xmalloc(grid->num_cells * max_num_vertices_per_cell *
              sizeof(*cell_to_vertex_ids));

    for (size_t i = 0; i < grid->num_cells; ++i) {

      size_t curr_num_vertices = grid->num_vertices_per_cell[i];
      yac_int * curr_cell_to_vertex_ids =
        cell_to_vertex_ids + i * max_num_vertices_per_cell;

      int cell_rank;
      if (curr_num_vertices > 0) {
        size_t * curr_cell_vertices =
          grid->cell_to_vertex + grid->cell_to_vertex_offsets[i];
        size_t min_vertex = curr_cell_vertices[0];
        curr_cell_to_vertex_ids[0] = grid->vertex_ids[min_vertex];
        for (size_t j = 1; j < curr_num_vertices; ++j) {
          size_t curr_vertex_idx = curr_cell_vertices[j];
          yac_int curr_vertex_id = grid->vertex_ids[curr_vertex_idx];
          insert_global_id(curr_cell_to_vertex_ids, j, curr_vertex_id);
          if (curr_cell_to_vertex_ids[0] == curr_vertex_id)
            min_vertex = curr_vertex_idx;
        }
        cell_rank = vertex_ranks[min_vertex];
      } else {
        cell_rank = 0;
      }
      for (size_t j = curr_num_vertices; j < max_num_vertices_per_cell; ++j)
        curr_cell_to_vertex_ids[j] = XT_INT_MAX;

      sendcounts[2 * ((cell_ranks[i] = cell_rank)) + 0]++;
    }
  }

  yac_int * edge_to_vertex_ids = NULL;

  if (!ids_available_global[1]) {

    edge_to_vertex_ids =
      xmalloc(2 * grid->num_edges * sizeof(*edge_to_vertex_ids));

    for (size_t i = 0; i < grid->num_edges; ++i) {

      size_t * curr_edge_to_vertex = &(grid->edge_to_vertex[i][0]);
      yac_int * curr_edge_vertex_ids = edge_to_vertex_ids + 2 * i;
      curr_edge_vertex_ids[0] = grid->vertex_ids[curr_edge_to_vertex[0]];
      curr_edge_vertex_ids[1] = grid->vertex_ids[curr_edge_to_vertex[1]];

      if (curr_edge_vertex_ids[0] > curr_edge_vertex_ids[1]) {
        yac_int temp = curr_edge_vertex_ids[0];
        curr_edge_vertex_ids[0] = curr_edge_vertex_ids[1];
        curr_edge_vertex_ids[1] = temp;
        sendcounts[
          2 * ((edge_ranks[i] = vertex_ranks[curr_edge_to_vertex[1]])) + 1]++;
      } else {
        sendcounts[
          2 * ((edge_ranks[i] = vertex_ranks[curr_edge_to_vertex[0]])) + 1]++;
      }
    }
  }

  // exchange the number of cells and edges
  yac_mpi_call(MPI_Alltoall(sendcounts, 2, YAC_MPI_SIZE_T,
                            recvcounts, 2, YAC_MPI_SIZE_T, comm), comm);

  total_sdispls[0] = 0;
  size_t recv_counts[2] = {0,0};
  size_t saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {
    total_sdispls[i+1] = saccu;
    total_rdispls[i] = raccu;
    recv_counts[0] += recvcounts[2 * i + 0];
    recv_counts[1] += recvcounts[2 * i + 1];
    total_sendcounts[i] = sendcounts[2 * i + 0] *
                          (size_t)max_num_vertices_per_cell +
                          sendcounts[2 * i + 1] * 2;
    total_recvcounts[i] = recvcounts[2 * i + 0] *
                          (size_t)max_num_vertices_per_cell +
                          recvcounts[2 * i + 1] * 2;
    saccu += total_sendcounts[i];
    raccu += total_recvcounts[i];
  }
  size_t local_data_count = total_sendcounts[comm_size - 1] +
                            total_sdispls[comm_size];
  size_t recv_count = total_recvcounts[comm_size - 1] +
                      total_rdispls[comm_size - 1];

  yac_int * yac_int_buffer =
    xmalloc((local_data_count + recv_count) * sizeof(*yac_int_buffer));
  yac_int * send_buffer = yac_int_buffer;
  yac_int * recv_buffer = yac_int_buffer + local_data_count;

  // pack send buffer
  if (!ids_available_global[0])
    for (size_t i = 0; i < grid->num_cells; ++i)
      for (int j = 0; j < max_num_vertices_per_cell; ++j)
        send_buffer[total_sdispls[cell_ranks[i] + 1]++] =
          cell_to_vertex_ids[i * max_num_vertices_per_cell + j];
  if (!ids_available_global[1])
    for (size_t i = 0; i < grid->num_edges; ++i)
      for (int j = 0; j < 2; ++j)
        send_buffer[total_sdispls[edge_ranks[i] + 1]++] =
          edge_to_vertex_ids[2 * i + j];

  free(edge_to_vertex_ids);
  free(cell_to_vertex_ids);

  // exchange data
  yac_alltoallv_yac_int_p2p(
    send_buffer, total_sendcounts, total_sdispls,
    recv_buffer, total_recvcounts, total_rdispls, comm);

  struct n_ids_reorder * n_ids_reorder_buffer =
    xmalloc((recv_counts[0] + recv_counts[1]) * sizeof(*n_ids_reorder_buffer));
  struct n_ids_reorder * n_ids_reorder[2] =
    {n_ids_reorder_buffer, n_ids_reorder_buffer + recv_counts[0]};

  size_t offset = 0;
  int index_counts[2] = {max_num_vertices_per_cell, 2};
  size_t reorder_idx = 0;
  recv_counts[0] = 0;
  recv_counts[1] = 0;
  for (int i = 0; i < comm_size; ++i) {
    for (int j = 0; j < 2; ++j) {
      size_t curr_count = recvcounts[2 * i + j];
      for (size_t k = 0; k < curr_count; ++k, ++reorder_idx, ++recv_counts[j]) {
        n_ids_reorder[j][recv_counts[j]].count = index_counts[j];
        n_ids_reorder[j][recv_counts[j]].ids = recv_buffer + offset;
        n_ids_reorder[j][recv_counts[j]].reorder_idx = reorder_idx;
        offset += index_counts[j];
      }
    }
  }

  for (int i = 0; i < 2; ++i) {

    if (ids_available_global[i]) continue;

    qsort(n_ids_reorder[i], recv_counts[i], sizeof(*(n_ids_reorder[i])),
          compare_n_ids_reorder_ids);

    size_t unique_count = recv_counts[i] > 0;
    struct n_ids_reorder * prev = n_ids_reorder[i];
    struct n_ids_reorder * curr = n_ids_reorder[i];

    for (size_t j = 0; j < recv_counts[i]; ++j, ++curr) {
      if (compare_n_ids_reorder_ids(prev, curr)) {
        ++unique_count;
        prev = curr;
      }
      curr->global_id = (yac_int)(unique_count - 1);
    }

    YAC_ASSERT(
      unique_count <= (size_t)XT_INT_MAX,
      "ERROR(generate_global_ce_ids): global_id out of bounds")

    yac_int yac_int_unique_count = (yac_int)unique_count;
    yac_int id_offset;

    // determine exclusive scan of sum of numbers of unique ids on all ranks
    yac_mpi_call(MPI_Exscan(&yac_int_unique_count, &id_offset, 1, yac_int_dt,
                            MPI_SUM, comm), comm);
    if (comm_rank == 0) id_offset = 0;

    YAC_ASSERT(
      ((size_t)id_offset + unique_count) <= (size_t)XT_INT_MAX,
      "ERROR(generate_global_ce_ids): global_id out of bounds")

    // adjust global ids
    for (size_t j = 0; j < recv_counts[i]; ++j)
      n_ids_reorder[i][j].global_id += id_offset;
  }
  free(yac_int_buffer);

  qsort(n_ids_reorder_buffer, recv_counts[0] + recv_counts[1],
        sizeof(*n_ids_reorder_buffer), compare_n_ids_reorder_reorder);

  yac_int * global_ids_buffer =
    xmalloc((recv_counts[0] + recv_counts[1] +
             ((ids_available_global[0])?0:(grid->num_cells)) +
             ((ids_available_global[1])?0:(grid->num_edges))) *
             sizeof(*global_ids_buffer));
  yac_int * send_global_ids = global_ids_buffer;
  yac_int * recv_global_ids =
    global_ids_buffer + recv_counts[0] + recv_counts[1];

  for (size_t i = 0; i < recv_counts[0] + recv_counts[1]; ++i)
    send_global_ids[i] = n_ids_reorder_buffer[i].global_id;
  free(n_ids_reorder_buffer);

  // generate count and displs data
  saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {
    total_sdispls[i] = saccu;
    total_rdispls[i] = raccu;
    saccu +=
      ((total_sendcounts[i] = recvcounts[2 * i + 0] + recvcounts[2 * i + 1]));
    raccu +=
      ((total_recvcounts[i] = sendcounts[2 * i + 0] + sendcounts[2 * i + 1]));
  }

  // exchange generated global ids data
  yac_alltoallv_yac_int_p2p(
    send_global_ids, total_sendcounts, total_sdispls,
    recv_global_ids, total_recvcounts, total_rdispls, comm);

  if ((!ids_available_global[0]) && (grid->num_cells > 0))
    grid->cell_ids = xmalloc(grid->num_cells * sizeof(*grid->cell_ids));
  if ((!ids_available_global[1]) && (grid->num_edges > 0))
    grid->edge_ids = xmalloc(grid->num_edges * sizeof(grid->edge_ids));

  // unpack generated global ids
  if (!ids_available_global[0])
    for (size_t i = 0; i < grid->num_cells; ++i)
      grid->cell_ids[i] = recv_global_ids[total_rdispls[cell_ranks[i]]++];
  if (!ids_available_global[1])
    for (size_t i = 0; i < grid->num_edges; ++i)
      grid->edge_ids[i] = recv_global_ids[total_rdispls[edge_ranks[i]]++];

  free(rank_buffer);
  free(size_t_buffer);
  free(global_ids_buffer);
}

// generates global cell/vertex/edge ids if they are not provided by the user
static void generate_global_ids(
  struct proc_sphere_part_node * proc_sphere_part,
  struct basic_grid_data * grid, MPI_Comm comm) {

  // generate global ids for all vertices (if non-existent)
  int * vertex_ranks = generate_vertex_ids(proc_sphere_part, grid,  comm);

  // generate global ids for all cell and edge (if non-existent)
  generate_ce_ids(proc_sphere_part, grid, vertex_ranks, comm);

  free(vertex_ranks);
}

// generate edge to cell mapping
static size_t_2_pointer generate_edge_to_cell(
  const_size_t_pointer cell_to_edge, const_int_pointer num_edges_per_cell,
  int * core_cell_mask, size_t num_cells, size_t num_edges) {

  if (num_cells == 0) return NULL;

  size_t_2_pointer edge_to_cell = xmalloc(num_edges * sizeof(*edge_to_cell));

  for (size_t i = 0; i < num_edges; ++i) {
    edge_to_cell[i][0] = SIZE_MAX;
    edge_to_cell[i][1] = SIZE_MAX;
  }

  for (size_t i = 0, offset = 0; i < num_cells; ++i) {

    size_t curr_num_edges = num_edges_per_cell[i];
    const_size_t_pointer curr_cell_to_edge = cell_to_edge + offset;
    offset += curr_num_edges;

    if ((core_cell_mask == NULL) || core_cell_mask[i]) {

      for (size_t j = 0; j < curr_num_edges; ++j) {

        size_t curr_edge = curr_cell_to_edge[j];
        size_t * curr_edge_to_cell = edge_to_cell[curr_edge];
        curr_edge_to_cell += *curr_edge_to_cell != SIZE_MAX;
        YAC_ASSERT(
          *curr_edge_to_cell == SIZE_MAX,
          "ERROR(generate_edge_to_cell): "
          "more than two cells point to a single edge "
          "(does the grid contain degenrated cells (less than 3 corners); "
          "these can be masked out using the core mask)")
        *curr_edge_to_cell = i;
      }
    }
  }

  return edge_to_cell;
}

static void get_vertex_cells(
  void * grid, size_t vertex_idx, size_t ** cells, int * num_cells) {

  *cells =
    ((struct basic_grid_data *)grid)->vertex_to_cell +
    ((struct basic_grid_data *)grid)->vertex_to_cell_offsets[vertex_idx];
  *num_cells =
    ((struct basic_grid_data *)grid)->num_cells_per_vertex[vertex_idx];
}

static void get_edge_cells(
  void * edge_to_cell, size_t edge_idx, size_t ** cells, int * num_cells) {

  *cells = ((size_t_2_pointer)edge_to_cell)[edge_idx];
  *num_cells = 2;
}

static void get_ve_ranks(
  int ** rank_buffer, size_t * rank_buffer_size,
  size_t * rank_buffer_array_size, int * num_ve_ranks,
  int * core_mask, size_t count,
  void(*get_cells)(
    void * ve_to_cell_data, size_t idx, size_t ** cells, int * num_cells),
  void * ve_to_cell_data, int * num_cell_ranks,
  size_t * dist_cell_rank_offsets) {

  // determine dist ranks for all vertices/edges based on the dist ranks of
  // adjacent cells
  for (size_t i = 0; i < count; ++i) {

    if (core_mask[i]) {

      size_t * cells;
      int num_cells;
      get_cells(ve_to_cell_data, i, &cells, &num_cells);

      size_t max_num_ranks = 0;
      for (int j = 0; j < num_cells; ++j)
        if (cells[j] != SIZE_MAX)
          max_num_ranks += (size_t)(num_cell_ranks[cells[j]]);

      ENSURE_ARRAY_SIZE(*rank_buffer, *rank_buffer_array_size,
                        *rank_buffer_size + (size_t)max_num_ranks);

      int * curr_ranks = *rank_buffer + *rank_buffer_size;
      int curr_num_ranks = 0;

      for (int j = 0; j < num_cells; ++j) {

        size_t curr_cell_idx = cells[j];
        if (curr_cell_idx == SIZE_MAX) continue;

        int curr_num_cell_ranks = num_cell_ranks[curr_cell_idx];
        int * curr_cell_ranks =
          *rank_buffer + dist_cell_rank_offsets[curr_cell_idx];

        for (int k = 0; k < curr_num_cell_ranks; ++k)
          insert_rank(curr_ranks, &curr_num_ranks, curr_cell_ranks[k]);
      }
      *rank_buffer_size += (size_t)(num_ve_ranks[i] = curr_num_ranks);
    } else {
      num_ve_ranks[i] = 0;
    }
  }
}

// generates mapping information from the decomposition provided by the user to
// one of the distributed grid
static struct remote_points ** generate_dist_remote_points(
  struct proc_sphere_part_node * proc_sphere_part,
  struct basic_grid_data * grid, MPI_Comm comm) {

  // generate global ids if they are missing
  generate_global_ids(proc_sphere_part, grid, comm);

  // determine dist ranks for all cells based on their bounding circles
  int * rank_buffer;
  int * num_ranks_buffer =
    xmalloc((grid->num_cells + grid->num_vertices + grid->num_edges) *
             sizeof(*num_ranks_buffer));
  int * num_cell_ranks = num_ranks_buffer;
  int * num_vertex_ranks = num_ranks_buffer + grid->num_cells;
  int * num_edge_ranks =
    num_ranks_buffer + grid->num_cells + grid->num_vertices;
  size_t * dist_cell_rank_offsets =
    xmalloc((grid->num_cells + 1) * sizeof(*dist_cell_rank_offsets));
  determine_dist_cell_ranks(
    proc_sphere_part, grid, comm,
    &rank_buffer, num_cell_ranks, dist_cell_rank_offsets);

  size_t rank_buffer_size = dist_cell_rank_offsets[grid->num_cells];
  size_t rank_buffer_array_size = dist_cell_rank_offsets[grid->num_cells];

  // determine dist ranks for all vertices based on the dist ranks of
  // adjacent cells
  get_ve_ranks(
    &rank_buffer, &rank_buffer_size, &rank_buffer_array_size, num_vertex_ranks,
    grid->core_vertex_mask, grid->num_vertices, get_vertex_cells,
    (void*)grid, num_cell_ranks, dist_cell_rank_offsets);

  // determine dist ranks for all edges based on the dist ranks of
  // adjacent cells
  size_t_2_pointer edge_to_cell =
    generate_edge_to_cell(
      grid->cell_to_edge, grid->num_vertices_per_cell,
      grid->core_cell_mask, grid->num_cells, grid->num_edges);
  get_ve_ranks(
    &rank_buffer, &rank_buffer_size, &rank_buffer_array_size, num_edge_ranks,
    grid->core_edge_mask, grid->num_edges, get_edge_cells,
    (void*)edge_to_cell, num_cell_ranks, dist_cell_rank_offsets);
  free(edge_to_cell);

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    3, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);
  size_t * size_t_buffer =
    xmalloc(4 * (size_t)comm_size * sizeof(*size_t_buffer));
  size_t * total_sendcounts = size_t_buffer + 0 * comm_size;
  size_t * total_recvcounts = size_t_buffer + 1 * comm_size;
  size_t * total_sdispls =    size_t_buffer + 2 * comm_size;
  size_t * total_rdispls =    size_t_buffer + 3 * comm_size;

  struct {
    size_t count;
    int * num_ranks;
    int * core_mask;
    yac_int * ids;
  } cve_data[3] =
    {{.count = grid->num_cells,
      .num_ranks = num_cell_ranks,
      .core_mask = grid->core_cell_mask,
      .ids = grid->cell_ids},
     {.count = grid->num_vertices,
      .num_ranks = num_vertex_ranks,
      .core_mask = grid->core_vertex_mask,
      .ids = grid->vertex_ids},
     {.count = grid->num_edges,
      .num_ranks = num_edge_ranks,
      .core_mask = grid->core_edge_mask,
      .ids = grid->edge_ids}};
  size_t offset = 0;

  // determine number of cells/vertices/edges that have to be
  // sent to other processes
  for (int location = 0; location < 3; ++location) {
    size_t count = cve_data[location].count;
    int * num_ranks = cve_data[location].num_ranks;
    int * core_mask = cve_data[location].core_mask;
    for (size_t i = 0; i < count; ++i) {
      int curr_num_ranks = num_ranks[i];
      int * ranks = rank_buffer + offset;
      offset += (size_t)curr_num_ranks;
      if (core_mask[i])
        for (int j = 0; j < curr_num_ranks; ++j)
          sendcounts[3 * ranks[j] + location]++;
    }
  }

  yac_generate_alltoallv_args(
    3, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t receive_counts[3] = {0,0,0};
  size_t saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {
    total_sdispls[i] = saccu;
    total_rdispls[i] = raccu;
    total_sendcounts[i] = 0;
    total_recvcounts[i] = 0;
    for (int location = 0; location < 3; ++location) {
      total_sendcounts[i] += sendcounts[3 * i + location];
      total_recvcounts[i] += recvcounts[3 * i + location];
      receive_counts[location] += recvcounts[3 * i + location];
    }
    saccu += total_sendcounts[i];
    raccu += total_recvcounts[i];
  }
  size_t local_data_count = total_sendcounts[comm_size - 1] +
                            total_sdispls[comm_size - 1];
  size_t recv_count = total_recvcounts[comm_size - 1] +
                      total_rdispls[comm_size - 1];

  struct id_pos * id_pos_buffer =
    xmalloc((local_data_count + recv_count) * sizeof(*id_pos_buffer));
  struct id_pos * id_pos_send_buffer = id_pos_buffer;
  struct id_pos * id_pos_recv_buffer =
    id_pos_buffer + local_data_count;

  // pack cell/vertex/edge information for distributed owners
  offset = 0;
  for (int location = 0; location < 3; ++location) {
    size_t count = cve_data[location].count;
    int * num_ranks = cve_data[location].num_ranks;
    int * core_mask = cve_data[location].core_mask;
    yac_int * ids = cve_data[location].ids;
    for (size_t i = 0; i < count; ++i) {
      int curr_num_ranks = num_ranks[i];
      int * ranks = rank_buffer + offset;
      offset += (size_t)curr_num_ranks;
      if (core_mask[i]) {
        yac_int global_id = ids[i];
        for (int j = 0; j < curr_num_ranks; ++j) {
          size_t pos = sdispls[3 * ranks[j] + location + 1]++;
          id_pos_send_buffer[pos].global_id = global_id;
          id_pos_send_buffer[pos].orig_pos = i;
        }
      }
    }
  }
  free(num_ranks_buffer);
  free(rank_buffer);
  free(dist_cell_rank_offsets);

  MPI_Datatype id_pos_dt = yac_get_id_pos_mpi_datatype(comm);

  // exchange cell/vertex/edge information for distributed owners
  yac_alltoallv_p2p(
    id_pos_send_buffer, total_sendcounts, total_sdispls,
    id_pos_recv_buffer, total_recvcounts, total_rdispls,
    sizeof(*id_pos_send_buffer), id_pos_dt, comm);

  yac_mpi_call(MPI_Type_free(&id_pos_dt), comm);

  size_t dist_owner_counts[3] = {0, 0, 0};
  for (int i = 0; i < comm_size; ++i)
    for (int location = 0; location < 3; ++location)
      dist_owner_counts[location] += recvcounts[3 * i + location];
  size_t max_dist_owner_count =
    MAX(MAX(dist_owner_counts[0], dist_owner_counts[1]), dist_owner_counts[2]);
  struct single_remote_point * temp_buffer =
    xmalloc(max_dist_owner_count * sizeof(*temp_buffer));

  struct remote_points ** dist_owners = xmalloc(3 * sizeof(*dist_owners));

  // unpack data
  for (int location = 0; location < 3; ++location) {

    size_t count = 0;
    for (int i = 0; i < comm_size; ++i) {
      size_t curr_recvcount = recvcounts[3 * i + location];
      struct id_pos * curr_id_pos =
        id_pos_recv_buffer + rdispls[3 * i + location];
      for (size_t k = 0; k < curr_recvcount; ++k, ++count) {
        temp_buffer[count].global_id = curr_id_pos[k].global_id;
        temp_buffer[count].data.orig_pos = curr_id_pos[k].orig_pos;
        temp_buffer[count].data.rank = i;
      }
    }

    // sort received global ids
    qsort(temp_buffer, count, sizeof(*temp_buffer),
          compare_single_remote_point_global_id);

    struct remote_point * unique_ids = xmalloc(count * sizeof(*unique_ids));
    size_t num_unique_ids = 0;

    // determine unique global ids
    yac_int prev_id = (count > 0)?temp_buffer[0].global_id - 1:-1;
    for (size_t i = 0; i < count; ++i) {

      yac_int curr_id = temp_buffer[i].global_id;
      if (curr_id != prev_id) {
        prev_id = curr_id;
        unique_ids[num_unique_ids].global_id = curr_id;
        unique_ids[num_unique_ids].data.count = 1;
        num_unique_ids++;
      } else {
        unique_ids[num_unique_ids-1].data.count++;
      }
    }

    size_t remote_point_info_buffer_size = 0;
    for (size_t i = 0; i < num_unique_ids; ++i)
      if (unique_ids[i].data.count > 1)
        remote_point_info_buffer_size +=
          (size_t)(unique_ids[i].data.count);

    dist_owners[location] =
      xmalloc(remote_point_info_buffer_size * sizeof(struct remote_point_info) +
              sizeof(**dist_owners));
    dist_owners[location]->data =
      (unique_ids = xrealloc(unique_ids, num_unique_ids * sizeof(*unique_ids)));
    dist_owners[location]->count = num_unique_ids;

    for (size_t i = 0, l = 0, offset = 0; i < num_unique_ids; ++i) {
      int curr_count = unique_ids[i].data.count;
      if (curr_count == 1) {
        unique_ids[i].data.data.single = temp_buffer[l].data;
        ++l;
      } else {
        unique_ids[i].data.data.multi = &(dist_owners[location]->buffer[offset]);
        for (int k = 0; k < curr_count; ++k, ++l, ++offset) {
          dist_owners[location]->buffer[offset] = temp_buffer[l].data;
        }
      }
    }

    qsort(
      unique_ids, num_unique_ids, sizeof(*unique_ids), compare_remote_point);
  }

  free(id_pos_buffer);
  free(temp_buffer);
  free(size_t_buffer);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  return dist_owners;
}

static void generate_sorted_ids(
  yac_int * ids, size_t count, yac_int ** sorted_ids, size_t ** reorder_idx) {

  yac_int * sorted_ids_ = xrealloc(*sorted_ids, count * sizeof(*sorted_ids_));
  size_t * reorder_idx_ = xrealloc(*reorder_idx, count * sizeof(*reorder_idx_));

  memcpy(sorted_ids_, ids, count * sizeof(*sorted_ids_));
  for (size_t i = 0; i < count; ++i) reorder_idx_[i] = i;

  yac_quicksort_index_yac_int_size_t(sorted_ids_, count, reorder_idx_);

  *sorted_ids = sorted_ids_;
  *reorder_idx = reorder_idx_;
}

static Xt_xmap generate_xmap_data(
  void * data, size_t count, MPI_Comm comm,
  void(*set_sendcounts)(void*,size_t,size_t*),
  void(*pack)(void*,size_t,size_t*,int*,int*)) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  set_sendcounts(data, count, sendcounts);

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);
  size_t num_src_msg = 0, num_dst_msg = 0;
  for (int i = 0; i < comm_size; ++i) {
    num_src_msg += (recvcounts[i] > 0);
    num_dst_msg += (sendcounts[i] > 0);
  }

  size_t recv_count =
    rdispls[comm_size-1] + recvcounts[comm_size-1];

  int * pos_buffer =
    xmalloc((recv_count + 2 * count) * sizeof(*pos_buffer));
  int * src_pos_buffer = pos_buffer;
  int * dst_pos_buffer = pos_buffer + recv_count;
  int * send_pos_buffer = pos_buffer + recv_count + count;

  // pack send buffer
  pack(data, count, sdispls, dst_pos_buffer, send_pos_buffer);

  // redistribute positions of requested data
  yac_alltoallv_int_p2p(
    send_pos_buffer, sendcounts, sdispls,
    src_pos_buffer, recvcounts, rdispls, comm);

  struct Xt_com_pos * com_pos =
    xmalloc(((size_t)num_src_msg + (size_t)num_dst_msg) * sizeof(*com_pos));
  struct Xt_com_pos * src_com = com_pos;
  struct Xt_com_pos * dst_com = com_pos + num_src_msg;

  // set transfer_pos pointers and transfer_pos counts in com_pos's
  num_src_msg = 0;
  num_dst_msg = 0;
  for (int i = 0; i < comm_size; ++i) {
    if (recvcounts[i] > 0) {
      src_com[num_src_msg].transfer_pos = src_pos_buffer;
      src_com[num_src_msg].num_transfer_pos = recvcounts[i];
      src_com[num_src_msg].rank = i;
      src_pos_buffer += recvcounts[i];
      ++num_src_msg;
    }
    if (sendcounts[i] > 0) {
      dst_com[num_dst_msg].transfer_pos = dst_pos_buffer;
      dst_com[num_dst_msg].num_transfer_pos = sendcounts[i];
      dst_com[num_dst_msg].rank = i;
      dst_pos_buffer += sendcounts[i];
      ++num_dst_msg;
    }
  }
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  Xt_xmap xmap =
    xt_xmap_intersection_pos_new(
      num_src_msg, src_com, num_dst_msg, dst_com, comm);

  free(com_pos);
  free(pos_buffer);

  return xmap;
}

static void generate_xmap_set_sendcounts_remote_points(
  void * data, size_t count, size_t * sendcounts) {

  struct remote_point * remote_point_data = (struct remote_point *)data;

  for (size_t i = 0; i < count; ++i) {
    struct remote_point_info * curr_info =
      (remote_point_data[i].data.count > 1)?
        (remote_point_data[i].data.data.multi):
        (&(remote_point_data[i].data.data.single));
    sendcounts[curr_info->rank]++;
  }
}

static void generate_xmap_pack_remote_points(
  void * data, size_t count, size_t * sdispls,
  int * dst_pos_buffer, int * send_pos_buffer) {

  struct remote_point * remote_point_data = (struct remote_point *)data;

  YAC_ASSERT(
    count <= INT_MAX,
    "ERROR(generate_xmap_pack_remote_points): count exceeds INT_MAX");

  for (size_t i = 0; i < count; ++i) {
    struct remote_point_info * curr_info =
      (remote_point_data[i].data.count > 1)?
        (remote_point_data[i].data.data.multi):
        (&(remote_point_data[i].data.data.single));
    size_t pos = sdispls[curr_info->rank+1]++;
    dst_pos_buffer[pos] = i;
    send_pos_buffer[pos] = (int)(curr_info->orig_pos);
  }
}

static Xt_xmap generate_xmap_from_remote_points(
  struct remote_point * remote_point_data, size_t count, MPI_Comm comm) {

  return
    generate_xmap_data(remote_point_data, count, comm,
                       generate_xmap_set_sendcounts_remote_points,
                       generate_xmap_pack_remote_points);
}

static void generate_xmap_set_sendcounts_remote_point_info(
  void * data, size_t count, size_t * sendcounts) {

  struct remote_point_info * point_infos = (struct remote_point_info *)data;

  for (size_t i = 0; i < count; ++i) sendcounts[point_infos[i].rank]++;
}

static void generate_xmap_pack_remote_point_info(
  void * data, size_t count, size_t * sdispls,
  int * dst_pos_buffer, int * send_pos_buffer) {

  struct remote_point_info * point_infos = (struct remote_point_info *)data;

  YAC_ASSERT(
    count <= INT_MAX,
    "ERROR(generate_xmap_pack_remote_point_info): count exceeds INT_MAX");

  for (size_t i = 0; i < count; ++i) {
    size_t pos = sdispls[point_infos[i].rank+1]++;
    dst_pos_buffer[pos] = i;
    send_pos_buffer[pos] = (int)(point_infos[i].orig_pos);
  }
}

static Xt_xmap generate_xmap_from_remote_point_info(
  struct remote_point_info * point_infos, size_t count, MPI_Comm comm) {

  return
    generate_xmap_data(point_infos, count, comm,
                       generate_xmap_set_sendcounts_remote_point_info,
                       generate_xmap_pack_remote_point_info);
}

static int compare_single_remote_point_reorder_global_id(
  const void * a, const void * b) {

  return (((const struct single_remote_point_reorder *)a)->data.global_id >
          ((const struct single_remote_point_reorder *)b)->data.global_id) -
         (((const struct single_remote_point_reorder *)a)->data.global_id <
          ((const struct single_remote_point_reorder *)b)->data.global_id);
}

static int compare_single_remote_point_reorder_reorder_idx(
  const void * a, const void * b) {

  return (((const struct single_remote_point_reorder *)a)->reorder_idx >
          ((const struct single_remote_point_reorder *)b)->reorder_idx) -
         (((const struct single_remote_point_reorder *)a)->reorder_idx <
          ((const struct single_remote_point_reorder *)b)->reorder_idx);
}

static void compute_basic_ve_data(
  size_t max_num_cell_ve, size_t num_cells, int * num_ve_per_cell,
  struct single_remote_point_reorder * cell_ve_point_data,
  yac_int ** ids, size_t * num_ids, size_t ** cell_to_ve, Xt_xmap * xmap,
  MPI_Comm comm) {

  size_t num_total_ve = num_cells * max_num_cell_ve;

  size_t total_num_cell_ve = 0;
  for (size_t i = 0, k = 0, l = 0; i < num_cells; ++i) {
    size_t num_vertices = (size_t)(num_ve_per_cell[i]);
    total_num_cell_ve += num_vertices;
    size_t j = 0;
    for (; j < num_vertices; ++j, ++k, ++l) {
      cell_ve_point_data[k].reorder_idx = l;
    }
    for (; j < max_num_cell_ve; ++j, ++k) {
      cell_ve_point_data[k].data.global_id = XT_INT_MAX;
      cell_ve_point_data[k].reorder_idx = SIZE_MAX;
    }
  }
  qsort(cell_ve_point_data, num_total_ve, sizeof(*cell_ve_point_data),
        compare_single_remote_point_reorder_global_id);
  yac_int * ids_ = xmalloc(total_num_cell_ve * sizeof(*ids_));
  size_t num_ids_ = 0;
  size_t * cell_to_ve_ =
    xmalloc(total_num_cell_ve * sizeof(*cell_to_ve_));
  struct remote_point_info * ve_point_info =
    xmalloc(total_num_cell_ve * sizeof(*ve_point_info));
  yac_int prev_global_id = XT_INT_MAX;
  for (size_t i = 0; i < num_total_ve; ++i) {
    yac_int curr_global_id = cell_ve_point_data[i].data.global_id;
    size_t curr_reorder_idx = cell_ve_point_data[i].reorder_idx;
    if (curr_global_id == XT_INT_MAX) break;
    if (curr_global_id != prev_global_id) {
      ids_[num_ids_] = (prev_global_id = curr_global_id);
      ve_point_info[num_ids_] = cell_ve_point_data[i].data.data;
      num_ids_++;
    }
    cell_to_ve_[curr_reorder_idx] = num_ids_ - 1;
  }

  *ids = xrealloc(ids_, num_ids_ * sizeof(*ids_));
  *num_ids = num_ids_;
  *cell_to_ve = cell_to_ve_;
  *xmap = generate_xmap_from_remote_point_info(ve_point_info, num_ids_, comm);
  free(ve_point_info);
}

static MPI_Datatype yac_get_coordinate_mpi_datatype(MPI_Comm comm) {

  MPI_Datatype coord_dt;
  yac_mpi_call(MPI_Type_contiguous(3, MPI_DOUBLE, &coord_dt), comm);
  yac_mpi_call(MPI_Type_commit(&coord_dt), comm);
  return coord_dt;
}

static struct yac_field_data field_data_init(
  struct yac_field_data * orig_field_data, size_t dist_size,
  Xt_redist redist_mask, Xt_redist redist_coords, MPI_Comm comm) {

  struct yac_field_data dist_field_data = yac_field_data_init();

  uint64_t counts[2], max_counts[2];
  if (orig_field_data != NULL) {
    counts[0] = orig_field_data->masks_count;
    counts[1] = orig_field_data->coordinates_count;
  } else {
    counts[0] = 0;
    counts[1] = 0;
  }
  yac_mpi_call(
    MPI_Allreduce(
      counts, max_counts, 2, MPI_UINT64_T, MPI_MAX, comm), comm);
  YAC_ASSERT(
    (orig_field_data == NULL) ||
    ((counts[0] == max_counts[0]) && (counts[1] == max_counts[1])),
    "ERROR(field_data_init): inconsistent number of masks or coordinates")

  int * data_available_flag =
    xcalloc(2 * max_counts[0] + max_counts[1], sizeof(*data_available_flag));

  if (orig_field_data != NULL) {
    for (size_t i = 0; i < orig_field_data->masks_count; ++i) {
      data_available_flag[i] = orig_field_data->masks[i].data != NULL;
      data_available_flag[i + orig_field_data->masks_count] =
        (orig_field_data->masks[i].name != NULL)?
          ((int)strlen(orig_field_data->masks[i].name)+1):0;
    }
    for (size_t i = 0; i < orig_field_data->coordinates_count; ++i)
      data_available_flag[i + 2 * orig_field_data->masks_count] =
        orig_field_data->coordinates[i] != NULL;
  }

  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, data_available_flag,
      (int)(2 * max_counts[0] + max_counts[1]), MPI_INT, MPI_MAX, comm), comm);

  if (orig_field_data != NULL) {
    for (size_t i = 0; i < orig_field_data->masks_count; ++i) {
      YAC_ASSERT(
        data_available_flag[i] == (orig_field_data->masks[i].data != NULL),
        "ERROR(field_data_init): inconsistent availability of masks")
      int mask_name_len =
        (orig_field_data->masks[i].name != NULL)?
          ((int)strlen(orig_field_data->masks[i].name)+1):0;
      YAC_ASSERT(
        data_available_flag[i + orig_field_data->masks_count] ==
          mask_name_len,
        "ERROR(field_data_init): inconsistent mask names")
    }

    for (size_t i = 0; i < orig_field_data->coordinates_count; ++i)
      YAC_ASSERT(
        data_available_flag[i + 2 * orig_field_data->masks_count] ==
        (orig_field_data->coordinates[i] != NULL),
        "ERROR(field_data_init): inconsistent availability of coordinates")
  }

  for (uint64_t i = 0; i < max_counts[0]; ++i) {
    int * dist_mask = NULL;
    if (data_available_flag[i]) {
      dist_mask = xmalloc(dist_size * sizeof(*dist_mask));
      int * orig_mask =
        (orig_field_data != NULL)?orig_field_data->masks[i].data:NULL;
      xt_redist_s_exchange1(redist_mask, orig_mask, dist_mask);
    }

    int mask_name_len = data_available_flag[i + max_counts[0]];
    char * mask_name = NULL;
    if (mask_name_len > 0) {
      mask_name = xmalloc((size_t)mask_name_len * sizeof(*mask_name));
      if ((orig_field_data != NULL) &&
          (orig_field_data->masks[i].name != NULL))
        memcpy(
          mask_name, orig_field_data->masks[i].name, (size_t)mask_name_len);
      else
        memset(mask_name, 0, (size_t)mask_name_len * sizeof(*mask_name));
      yac_mpi_call(
        MPI_Allreduce(
          MPI_IN_PLACE, mask_name, mask_name_len, MPI_CHAR, MPI_MAX, comm),
        comm);
      YAC_ASSERT(
        (orig_field_data == NULL) ||
        (orig_field_data->masks[i].name == NULL) ||
        !memcmp(
           orig_field_data->masks[i].name, mask_name,
           (size_t)mask_name_len * sizeof(*mask_name)),
        "ERROR(field_data_init): inconsistent mask names")
    }
    yac_field_data_add_mask_nocpy(&dist_field_data, dist_mask, mask_name);
  }

  for (uint64_t i = 0; i < max_counts[1]; ++i) {
    coordinate_pointer dist_coordinates = NULL;
    if (data_available_flag[i + 2 * max_counts[0]]) {
      dist_coordinates = xmalloc(dist_size * sizeof(*dist_coordinates));
      coordinate_pointer orig_coordinates =
        (orig_field_data != NULL)?orig_field_data->coordinates[i]:NULL;
      xt_redist_s_exchange1(redist_coords, orig_coordinates, dist_coordinates);
    }
    yac_field_data_add_coordinates_nocpy(&dist_field_data, dist_coordinates);
  }

  free(data_available_flag);

  return dist_field_data;
}

static struct dist_grid generate_dist_grid(
  struct proc_sphere_part_node * proc_sphere_part,
  struct yac_basic_grid * grid, MPI_Comm comm) {

  struct basic_grid_data * grid_data = yac_basic_grid_get_data(grid);

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  int max_num_vertices_per_cell = 0;
  { // compute the global maximum number of vertices per cell
    for (size_t i = 0; i < grid_data->num_cells; ++i)
      if (grid_data->num_vertices_per_cell[i] > max_num_vertices_per_cell)
        max_num_vertices_per_cell = grid_data->num_vertices_per_cell[i];
    yac_mpi_call(MPI_Allreduce(MPI_IN_PLACE, &max_num_vertices_per_cell, 1,
                               MPI_INT, MPI_MAX, comm), comm);
  }

  struct remote_points ** dist_owner =
    generate_dist_remote_points(proc_sphere_part, grid_data, comm);

  Xt_xmap xmap_cell_orig_to_dist =
    generate_xmap_from_remote_points(
      dist_owner[CELL]->data, dist_owner[CELL]->count, comm);

  MPI_Datatype dt_single_remote_point =
    yac_get_n_single_remote_point_reorder_mpi_datatype(
      max_num_vertices_per_cell, comm);
  MPI_Datatype dt_coord = yac_get_coordinate_mpi_datatype(comm);

  Xt_redist redist_cell_yac_int,
            redist_cell_int,
            redist_cell_corner_point_data,
            redist_cell_coords;
  redist_cell_yac_int = xt_redist_p2p_new(xmap_cell_orig_to_dist, yac_int_dt);
  redist_cell_int = xt_redist_p2p_new(xmap_cell_orig_to_dist, MPI_INT);
  redist_cell_corner_point_data =
    xt_redist_p2p_new(xmap_cell_orig_to_dist, dt_single_remote_point);
  yac_mpi_call(MPI_Type_free(&dt_single_remote_point), comm);
  redist_cell_coords = xt_redist_p2p_new(xmap_cell_orig_to_dist, dt_coord);

  xt_xmap_delete(xmap_cell_orig_to_dist);

  size_t num_cells = dist_owner[CELL]->count;
  size_t max_num_cell_corners =
    (size_t)max_num_vertices_per_cell * num_cells;
  size_t max_num_local_cell_corners =
    (size_t)max_num_vertices_per_cell * grid_data->num_cells;
  yac_int * cell_ids = xmalloc(num_cells * sizeof(*cell_ids));
  int * num_vertices_per_cell =
    xmalloc(num_cells * sizeof(*num_vertices_per_cell));
  struct single_remote_point_reorder *
    cell_corner_point_data =
    xmalloc(2 * (max_num_cell_corners + max_num_local_cell_corners) *
            sizeof(*cell_corner_point_data));
  struct single_remote_point_reorder * cell_edge_point_data =
    cell_corner_point_data + max_num_cell_corners;
  struct single_remote_point_reorder * local_cell_corner_point_data =
    cell_corner_point_data + 2 * max_num_cell_corners;
  struct single_remote_point_reorder * local_cell_edge_point_data =
    cell_corner_point_data + 2 * max_num_cell_corners +
    max_num_local_cell_corners;
  for (size_t i = 0; i < grid_data->num_cells; ++i){
    int num_corners = (size_t)(grid_data->num_vertices_per_cell[i]);
    size_t * corners = grid_data->cell_to_vertex + grid_data->cell_to_vertex_offsets[i];
    size_t * edges = grid_data->cell_to_edge + grid_data->cell_to_edge_offsets[i];
    struct single_remote_point_reorder * curr_local_cell_corner_point_data =
      local_cell_corner_point_data + i * (size_t)max_num_vertices_per_cell;
    struct single_remote_point_reorder * curr_local_cell_edge_point_data =
      local_cell_edge_point_data + i * (size_t)max_num_vertices_per_cell;
    for (size_t j = 0; j < num_corners; ++j) {
      curr_local_cell_corner_point_data[j].data.global_id = grid_data->vertex_ids[corners[j]];
      curr_local_cell_corner_point_data[j].data.data.rank = comm_rank;
      curr_local_cell_corner_point_data[j].data.data.orig_pos = corners[j];
      curr_local_cell_edge_point_data[j].data.global_id = grid_data->edge_ids[edges[j]];
      curr_local_cell_edge_point_data[j].data.data.rank = comm_rank;
      curr_local_cell_edge_point_data[j].data.data.orig_pos = edges[j];
    }
  }

  xt_redist_s_exchange1(redist_cell_yac_int, grid_data->cell_ids, cell_ids);
  xt_redist_s_exchange1(
    redist_cell_int, grid_data->num_vertices_per_cell, num_vertices_per_cell);
  xt_redist_s_exchange1(
    redist_cell_corner_point_data,
    local_cell_corner_point_data, cell_corner_point_data);
  xt_redist_s_exchange1(
    redist_cell_corner_point_data,
    local_cell_edge_point_data, cell_edge_point_data);

  struct yac_field_data cell_field_data =
    field_data_init(
      yac_basic_grid_get_field_data(grid, CELL), num_cells,
      redist_cell_int, redist_cell_coords, comm);

  xt_redist_delete(redist_cell_coords);
  xt_redist_delete(redist_cell_corner_point_data);
  xt_redist_delete(redist_cell_yac_int);
  xt_redist_delete(redist_cell_int);

  yac_int * vertex_ids;
  size_t num_vertices;
  size_t * cell_to_vertex;
  size_t * cell_to_vertex_offsets =
    xmalloc(num_cells * sizeof(*cell_to_vertex_offsets));
  Xt_xmap xmap_vertex_orig_to_dist;
  for (size_t i = 0, accu = 0; i < num_cells; ++i) {
    cell_to_vertex_offsets[i] = accu;
    accu += (size_t)(num_vertices_per_cell[i]);
  }
  compute_basic_ve_data(
    max_num_vertices_per_cell, num_cells, num_vertices_per_cell,
    cell_corner_point_data, &vertex_ids, &num_vertices,
    &cell_to_vertex, &xmap_vertex_orig_to_dist, comm);

  yac_int * edge_ids;
  size_t num_edges;
  size_t * cell_to_edge;
  size_t * cell_to_edge_offsets = cell_to_vertex_offsets;
  Xt_xmap xmap_edge_orig_to_dist;
  compute_basic_ve_data(
    max_num_vertices_per_cell, num_cells, num_vertices_per_cell,
    cell_edge_point_data, &edge_ids, &num_edges,
    &cell_to_edge, &xmap_edge_orig_to_dist, comm);

  free(cell_corner_point_data);

  // generate redists for redistribution of vertex data
  Xt_redist redist_vertex_coords =
    xt_redist_p2p_new(xmap_vertex_orig_to_dist, dt_coord);
  Xt_redist redist_vertex_int =
    xt_redist_p2p_new(xmap_vertex_orig_to_dist, MPI_INT);
  xt_xmap_delete(xmap_vertex_orig_to_dist);

  // get vertex coordinates
  coordinate_pointer vertex_coordinates =
    xmalloc(num_vertices * sizeof(*vertex_coordinates));
  xt_redist_s_exchange1(
    redist_vertex_coords, grid_data->vertex_coordinates, vertex_coordinates);

  struct yac_field_data vertex_field_data =
    field_data_init(
      yac_basic_grid_get_field_data(grid, CORNER), num_vertices,
      redist_vertex_int, redist_vertex_coords, comm);

  xt_redist_delete(redist_vertex_int);
  xt_redist_delete(redist_vertex_coords);

  Xt_redist redist_edge_int,
            redist_edge_coords,
            redist_edge_2yac_int;
  { // generate redists for redistribution of edge data

    MPI_Datatype dt_2yac_int;
    yac_mpi_call(MPI_Type_contiguous(2, yac_int_dt, &dt_2yac_int), comm);
    yac_mpi_call(MPI_Type_commit(&dt_2yac_int), comm);

    redist_edge_int = xt_redist_p2p_new(xmap_edge_orig_to_dist, MPI_INT);
    redist_edge_coords = xt_redist_p2p_new(xmap_edge_orig_to_dist, dt_coord);
    redist_edge_2yac_int = xt_redist_p2p_new(xmap_edge_orig_to_dist, dt_2yac_int);

    yac_mpi_call(MPI_Type_free(&dt_2yac_int), comm);

  }
  xt_xmap_delete(xmap_edge_orig_to_dist);

  enum yac_edge_type * edge_type = xmalloc(num_edges * sizeof(*edge_type));

  { // get edge types
    int * temp_edge_type_src, * temp_edge_type_dst;
    if (sizeof(*edge_type) == sizeof(int)) {
      temp_edge_type_src = (int*)(grid_data->edge_type);
      temp_edge_type_dst = (int*)edge_type;
    } else {
        temp_edge_type_src =
          xmalloc(grid_data->num_edges * sizeof(*temp_edge_type_src));
        for (size_t i = 0; i < grid_data->num_edges; ++i)
          temp_edge_type_src[i] = (int)(grid_data->edge_type[i]);
        temp_edge_type_dst = xmalloc(num_edges * sizeof(*temp_edge_type_dst));
        for (size_t i = 0; i < num_edges; ++i)
          temp_edge_type_dst[i] = (int)(edge_type[i]);
    }

    xt_redist_s_exchange1(
      redist_edge_int, temp_edge_type_src, temp_edge_type_dst);

    if (sizeof(*edge_type) != sizeof(int)) {

      for (size_t i = 0; i < num_edges; ++i)
        edge_type[i] = (enum yac_edge_type)(temp_edge_type_dst[i]);

      free(temp_edge_type_src);
      free(temp_edge_type_dst);
    }
  }

  struct yac_field_data edge_field_data =
    field_data_init(
      yac_basic_grid_get_field_data(grid, EDGE), num_edges,
      redist_edge_int, redist_edge_coords, comm);

  xt_redist_delete(redist_edge_int);
  xt_redist_delete(redist_edge_coords);

  struct bounding_circle * cell_bnd_circles =
    xmalloc(num_cells * sizeof(*cell_bnd_circles));
  {
    struct grid_cell cell;
    cell.coordinates_xyz = xmalloc((size_t)max_num_vertices_per_cell *
                                   sizeof(*(cell.coordinates_xyz)));
    cell.edge_type = xmalloc((size_t)max_num_vertices_per_cell *
                             sizeof(*(cell.edge_type)));

    for (size_t i = 0; i < num_cells; ++i) {
      size_t * curr_cell_to_vertex =
        cell_to_vertex + cell_to_vertex_offsets[i];
      size_t * curr_cell_to_edge =
        cell_to_edge + cell_to_edge_offsets[i];
      for (int j = 0; j < num_vertices_per_cell[i]; ++j) {
        coordinate_pointer curr_vertex_coords =
          vertex_coordinates + curr_cell_to_vertex[j];
        for (int k = 0; k < 3; ++k)
          cell.coordinates_xyz[j][k] = (*curr_vertex_coords)[k];
        cell.edge_type[j] = edge_type[curr_cell_to_edge[j]];
      }
      cell.num_corners = num_vertices_per_cell[i];
      if (cell.num_corners > 0)
        yac_get_cell_bounding_circle(cell, cell_bnd_circles + i);
      else
        cell_bnd_circles[i] =
          (struct bounding_circle) {
            .base_vector = {1.0, 0.0, 0.0},
            .inc_angle = SIN_COS_ZERO,
            .sq_crd = DBL_MAX};
    }
    free(cell.edge_type);
    free(cell.coordinates_xyz);
  }

  size_t_2_pointer edge_to_vertex =
    xmalloc(num_edges * sizeof(*edge_to_vertex));

  // get edge to vertex
  {
    yac_int * vertex_id_buffer =
      xmalloc(
        2 * (grid_data->num_edges + num_edges) * sizeof(*vertex_id_buffer));
    yac_int * grid_edge_vertex_ids = vertex_id_buffer;
    yac_int * edge_vertex_ids = vertex_id_buffer + 2 * grid_data->num_edges;

    size_t * grid_edge_to_vertex = &(grid_data->edge_to_vertex[0][0]);

    for (size_t i = 0; i < 2 * grid_data->num_edges; ++i)
      grid_edge_vertex_ids[i] =
        grid_data->vertex_ids[grid_edge_to_vertex[i]];

    xt_redist_s_exchange1(
      redist_edge_2yac_int, grid_edge_vertex_ids, edge_vertex_ids);

    id2idx(edge_vertex_ids, &(edge_to_vertex[0][0]), 2 * num_edges,
           vertex_ids, num_vertices);

    free(vertex_id_buffer);
  }

  xt_redist_delete(redist_edge_2yac_int);

  yac_mpi_call(MPI_Type_free(&dt_coord), comm);


  struct dist_grid dist_grid =
    {.comm = comm,
     .vertex_coordinates = vertex_coordinates,
     .ids = {cell_ids, vertex_ids, edge_ids},
     .total_count = {num_cells, num_vertices, num_edges},
     .count = {num_cells, num_vertices, num_edges},
     .num_vertices_per_cell = num_vertices_per_cell,
     .cell_to_vertex = cell_to_vertex,
     .cell_to_vertex_offsets = cell_to_vertex_offsets,
     .cell_to_edge = cell_to_edge,
     .cell_to_edge_offsets = cell_to_edge_offsets,
     .edge_to_vertex = edge_to_vertex,
     .cell_bnd_circles = cell_bnd_circles,
     .edge_type = edge_type,
     .sorted_ids = {NULL, NULL, NULL},
     .sorted_reorder_idx = {NULL, NULL, NULL}};
  generate_sorted_ids(cell_ids, num_cells, &(dist_grid.sorted_ids[CELL]),
                      &(dist_grid.sorted_reorder_idx[CELL]));
  generate_sorted_ids(vertex_ids, num_vertices, &(dist_grid.sorted_ids[CORNER]),
                      &(dist_grid.sorted_reorder_idx[CORNER]));
  generate_sorted_ids(edge_ids, num_edges, &(dist_grid.sorted_ids[EDGE]),
                      &(dist_grid.sorted_reorder_idx[EDGE]));
  dist_grid.field_data.cell = cell_field_data;
  dist_grid.field_data.vertex = vertex_field_data;
  dist_grid.field_data.edge = edge_field_data;

  generate_dist_grid_remote_point_infos(
    &dist_grid, proc_sphere_part, dist_owner);

  generate_core_masks(&dist_grid, proc_sphere_part, comm_rank);

  for (int i = 0; i < 3; ++i) {
    free(dist_owner[i]->data);
    free(dist_owner[i]);
  }
  free(dist_owner);

  return dist_grid;
}

static void setup_search_data(struct dist_grid_pair * dist_grid_pair) {

  // build sphere part for vertices
  for (int i = 0; i < 2; ++i)
    dist_grid_pair->vertex_sphere_part[i] =
      yac_point_sphere_part_search_new(
        dist_grid_pair->dist_grid[i].count[CORNER],
        dist_grid_pair->dist_grid[i].vertex_coordinates,
        dist_grid_pair->dist_grid[i].ids[CORNER]);

  // build sphere part for bounding circle of cells
  for (int i = 0; i < 2; ++i) {
    dist_grid_pair->cell_sphere_part[i] =
      yac_bnd_sphere_part_search_new(
        dist_grid_pair->dist_grid[i].cell_bnd_circles,
        dist_grid_pair->dist_grid[i].count[CELL]);
  }
}

static void get_grid_cells(
  struct yac_basic_grid * grid, struct dist_cell * cells, MPI_Comm comm) {

  struct basic_grid_data * grid_data = yac_basic_grid_get_data(grid);

  size_t num_cells = grid_data->num_cells;
  for (size_t i = 0; i < num_cells; ++i) {

    double * cell_coord = cells[i].coord;
    for (int j = 0; j < 3; ++j) cell_coord[j] = 0.0;

    if (grid_data->num_vertices_per_cell[i] == 0) continue;

    size_t * curr_cell_to_vertex =
      grid_data->cell_to_vertex + grid_data->cell_to_vertex_offsets[i];
    // compute the middle point of the cell
    for (int j = 0; j < grid_data->num_vertices_per_cell[i]; ++j) {
      coordinate_pointer curr_vertex_coords =
        grid_data->vertex_coordinates + curr_cell_to_vertex[j];
      for (int k = 0; k < 3; ++k)
        cell_coord[k] += (*curr_vertex_coords)[k];
    }
    normalise_vector(cell_coord);
  }
}

static struct proc_sphere_part_node *
generate_dist_grid_decomposition(
  struct yac_basic_grid * grid_a, struct yac_basic_grid * grid_b,
  MPI_Comm comm) {

  size_t num_cells_a = yac_basic_grid_get_data(grid_a)->num_cells;
  size_t num_cells_b = yac_basic_grid_get_data(grid_b)->num_cells;

  // gather cell centers from both grids into single array
  size_t num_cells = num_cells_a + num_cells_b;
  struct dist_cell * cells = xmalloc(num_cells * sizeof(*cells));
  get_grid_cells(grid_a, cells, comm);
  get_grid_cells(grid_b, cells + num_cells_a, comm);

  // redistribute all cell centers and build parallel sphere
  // part for vertices
  struct proc_sphere_part_node * proc_sphere_part =
    yac_redistribute_cells(&cells, &num_cells, comm);

  free(cells);

  return proc_sphere_part;
}

struct dist_grid_pair * yac_dist_grid_pair_new(
  struct yac_basic_grid * grid_a, struct yac_basic_grid * grid_b,
  MPI_Comm comm) {

  YAC_ASSERT(
    grid_a,
    "ERROR(yac_dist_grid_pair_new): "
    "NULL is not a valid value for parameter grid_a")
  YAC_ASSERT(
    grid_b,
    "ERROR(yac_dist_grid_pair_new): "
    "NULL is not a valid value for parameter grid_b")

  YAC_ASSERT(
    strcmp(yac_basic_grid_get_name(grid_a),
           yac_basic_grid_get_name(grid_b)),
    "ERROR(yac_dist_grid_pair_new): identical grid names")

  MPI_Comm comm_copy;
  yac_mpi_call(MPI_Comm_dup(comm, &comm_copy), comm);
  comm = comm_copy;

  // ensure same grid ordering on all processes
  if (strcmp(yac_basic_grid_get_name(grid_a),
             yac_basic_grid_get_name(grid_b)) > 0) {
    struct yac_basic_grid * grid_swap = grid_a;
    grid_a = grid_b;
    grid_b = grid_swap;
  }

  struct dist_grid_pair * grid_pair = xmalloc(1 * sizeof(*grid_pair));

  strncpy((char*)grid_pair->grid_names[0],
          yac_basic_grid_get_name(grid_a), YAC_MAX_CHARLEN-1);
  strncpy((char*)grid_pair->grid_names[1],
          yac_basic_grid_get_name(grid_b), YAC_MAX_CHARLEN-1);
  grid_pair->comm = comm;

  // generate a decomposition for the distributed grid pair
  // with the following properties:
  //   * each process covers a unique area on the sphere
  //   * number of cells (from grid_a and/or grid_b) per area
  //     is roughly the same
  grid_pair->proc_sphere_part =
    generate_dist_grid_decomposition(grid_a, grid_b, comm);

  // redistribute grid_a and grid_b according to decomposition
  grid_pair->dist_grid[0] =
    generate_dist_grid(grid_pair->proc_sphere_part, grid_a, comm);
  grid_pair->dist_grid[1] =
    generate_dist_grid(grid_pair->proc_sphere_part, grid_b, comm);

  // build search data structures for local cells and vertices in
  // distributed grid
  setup_search_data(grid_pair);

  return grid_pair;
}

MPI_Comm yac_dist_grid_pair_get_MPI_Comm(struct dist_grid_pair * grid_pair) {
  return grid_pair->comm;
}

struct dist_grid * yac_dist_grid_pair_get_dist_grid(
  struct dist_grid_pair * grid_pair, char const * grid_name) {

  YAC_ASSERT_F(
    strlen(grid_name) < YAC_MAX_CHARLEN,
    "ERROR(yac_dist_grid_pair_get_dist_grid): grid name \"%s\"is too long "
    "(maximum length is YAC_MAX_CHARLEN)", grid_name)

  struct dist_grid * dist_grid = NULL;
  for (int i = 0; (i < 2) && (dist_grid == NULL); ++i)
    if (!strcmp(grid_name, grid_pair->grid_names[i]))
      dist_grid = &(grid_pair->dist_grid[i]);
  YAC_ASSERT(
    dist_grid, "ERROR(yac_dist_grid_pair_get_dist_grid): invalid grid_name")
  return dist_grid;
}

struct const_basic_grid_data * yac_dist_grid_get_basic_grid_data(
  struct dist_grid * dist_grid) {

  return (struct const_basic_grid_data *)dist_grid;
}

size_t yac_dist_grid_get_local_count(
  struct dist_grid * dist_grid, enum yac_location location) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_dist_grid_get_local_count): invalid location")
  size_t local_count = 0;
  size_t count = dist_grid->count[location];
  int * core_mask = dist_grid->core_mask[location];
  for (size_t i = 0; i < count; ++i) if (core_mask[i]) ++local_count;
  return local_count;
}

static size_t yac_dist_grid_get_count(
  struct dist_grid * dist_grid, enum yac_location location) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_dist_grid_get_count): invalid location")
  return dist_grid->total_count[location];
}

static int const * yac_dist_grid_get_core_mask(
  struct dist_grid * dist_grid, enum yac_location location) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_dist_grid_get_core_mask): invalid location")
  return dist_grid->core_mask[location];
}

static yac_int const * yac_dist_grid_get_global_ids(
  struct dist_grid * dist_grid, enum yac_location location) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_dist_grid_get_global_ids): invalid location")
  return dist_grid->ids[location];
}

void yac_dist_grid_get_local_unmasked_points(
  struct dist_grid * dist_grid, struct interp_field field,
  size_t ** indices, size_t * count) {

  int const * mask = yac_dist_grid_get_field_mask(dist_grid, field);

  size_t total_count = yac_dist_grid_get_count(dist_grid, field.location);
  int const * core_mask =
    yac_dist_grid_get_core_mask(dist_grid, field.location);

  size_t * temp_indices = xmalloc(total_count * sizeof(*temp_indices));

  size_t count_ = 0;

  if (mask != NULL) {
    for (size_t i = 0; i < total_count; ++i)
      if (core_mask[i] && mask[i]) temp_indices[count_++] = i;
  } else {
    for (size_t i = 0; i < total_count; ++i)
      if (core_mask[i]) temp_indices[count_++] = i;
  }

  *indices = xrealloc(temp_indices, count_ * sizeof(**indices));
  *count = count_;
}

static struct yac_field_data * yac_dist_grid_get_field_data(
  struct dist_grid * dist_grid, enum yac_location location) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_dist_grid_get_field_data): invalid location")

  switch (location) {
    default:
    case(CELL): return &(dist_grid->field_data.cell);
    case(CORNER): return &(dist_grid->field_data.vertex);
    case(EDGE): return &(dist_grid->field_data.edge);
  };
}

int const * yac_dist_grid_get_field_mask(
  struct dist_grid * dist_grid, struct interp_field field) {

  if (field.masks_idx == SIZE_MAX) return NULL;

  struct yac_field_data * data =
    yac_dist_grid_get_field_data(dist_grid, field.location);

  YAC_ASSERT(
    field.masks_idx < data->masks_count,
    "ERROR(yac_dist_grid_get_field_mask): invalid mask_idx")

  return data->masks[field.masks_idx].data;
}

coordinate_pointer yac_dist_grid_get_field_coords(
  struct dist_grid * dist_grid, struct interp_field field) {

  coordinate_pointer coords = NULL;

  if (field.coordinates_idx != SIZE_MAX) {

    struct yac_field_data * data =
      yac_dist_grid_get_field_data(dist_grid, field.location);

    YAC_ASSERT(
      field.coordinates_idx < data->coordinates_count,
      "ERROR(yac_dist_grid_get_field_coords): invalid coord_idx")

    coords = data->coordinates[field.coordinates_idx];
  }

  // if no field coordinates are defined, but the location is at the corners of
  // of the grid cells, return coordinates of them
  if ((coords == NULL) && (field.location == CORNER))
    coords = dist_grid->vertex_coordinates;

  return coords;
}

size_t yac_dist_grid_get_unmasked_local_count(
  struct dist_grid * dist_grid, struct interp_field field) {

  size_t total_count = yac_dist_grid_get_count(dist_grid, field.location);

  int const * mask = yac_dist_grid_get_field_mask(dist_grid, field);
  if (mask == NULL)
    return yac_dist_grid_get_local_count(dist_grid, field.location);

  int const * core_mask =
    yac_dist_grid_get_core_mask(dist_grid, field.location);

  size_t unmasked_local_count = 0;
  for (size_t i = 0; i < total_count; ++i)
    if (core_mask[i] && mask[i]) ++unmasked_local_count;
  return unmasked_local_count;
}

static void yac_remote_point_infos_free(
  struct remote_point_infos * point_infos, size_t count) {

  for (size_t i = 0; i < count; ++i)
    if (point_infos[i].count > 1) free(point_infos[i].data.multi);
  free(point_infos);
}

static void yac_dist_grid_free(struct dist_grid grid) {

  free(grid.vertex_coordinates);
  free(grid.num_vertices_per_cell);
  free(grid.cell_to_vertex);
  free(grid.cell_to_vertex_offsets);
  free(grid.cell_to_edge);
  free(grid.cell_bnd_circles);
  free(grid.edge_type);
  free(grid.edge_to_vertex);
  for (int i = 0; i < 3; ++i) {
    free(grid.ids[i]);
    free(grid.core_mask[i]);
    free(grid.sorted_ids[i]);
    free(grid.sorted_reorder_idx[i]);
    yac_remote_point_infos_free(grid.owners[i], grid.total_count[i]);
  }
  yac_field_data_set_free(grid.field_data);
}

void yac_dist_grid_pair_delete(struct dist_grid_pair * grid_pair) {

  if (grid_pair == NULL) return;
  yac_mpi_call(MPI_Comm_free(&(grid_pair->comm)), MPI_COMM_WORLD);
  yac_proc_sphere_part_node_delete(grid_pair->proc_sphere_part);
  for (int i = 0; i < 2; ++i) {
    yac_dist_grid_free(grid_pair->dist_grid[i]);
    yac_delete_point_sphere_part_search(grid_pair->vertex_sphere_part[i]);
    yac_bnd_sphere_part_search_delete(grid_pair->cell_sphere_part[i]);
  }
  free(grid_pair);
}

static int coord_in_cell(double coord[3], struct dist_grid * dist_grid,
                         size_t cell_idx, struct grid_cell * buffer_cell) {

  yac_const_basic_grid_data_get_grid_cell(
    (struct const_basic_grid_data *)dist_grid, cell_idx, buffer_cell);

  return
    yac_point_in_cell2(
      coord, *buffer_cell, dist_grid->cell_bnd_circles[cell_idx]);
}

static int coord_in_cell_gc(double coord[3], struct dist_grid * dist_grid,
                         size_t cell_idx, struct grid_cell * buffer_cell) {

  yac_const_basic_grid_data_get_grid_cell(
    (struct const_basic_grid_data *)dist_grid, cell_idx, buffer_cell);
  for (size_t i = 0; i < buffer_cell->num_corners; ++i)
    buffer_cell->edge_type[i] = GREAT_CIRCLE_EDGE;

  return
    yac_point_in_cell2(
      coord, *buffer_cell, dist_grid->cell_bnd_circles[cell_idx]);
}

void yac_const_basic_grid_data_get_grid_cell(
  struct const_basic_grid_data * grid_data, size_t cell_idx,
  struct grid_cell * buffer_cell) {

  size_t num_vertices = (size_t)(grid_data->num_vertices_per_cell[cell_idx]);

  struct grid_cell cell = *buffer_cell;

  if (cell.array_size < num_vertices) {
    cell.coordinates_xyz =
      xrealloc(cell.coordinates_xyz, num_vertices *
               sizeof(*(cell.coordinates_xyz)));
    cell.edge_type = xrealloc(cell.edge_type, num_vertices *
                                      sizeof(*(cell.edge_type)));
    cell.array_size = num_vertices;
    *buffer_cell = cell;
  }

  for (size_t i = 0; i < num_vertices; ++i) {
    size_t vertex_idx =
      grid_data->cell_to_vertex[
        grid_data->cell_to_vertex_offsets[cell_idx] + i];
    cell.coordinates_xyz[i][0] = grid_data->vertex_coordinates[vertex_idx][0];
    cell.coordinates_xyz[i][1] = grid_data->vertex_coordinates[vertex_idx][1];
    cell.coordinates_xyz[i][2] = grid_data->vertex_coordinates[vertex_idx][2];
    size_t edge_idx =
      grid_data->cell_to_edge[grid_data->cell_to_edge_offsets[cell_idx] + i];
    cell.edge_type[i] = grid_data->edge_type[edge_idx];
  }
  buffer_cell->num_corners = num_vertices;
}

static struct bnd_sphere_part_search * dist_grid_pair_get_cell_sphere_part(
  struct dist_grid_pair * grid_pair, char const * grid_name) {

  YAC_ASSERT_F(
    strlen(grid_name) < YAC_MAX_CHARLEN,
    "ERROR(dist_grid_pair_get_cell_sphere_part): grid name \"%s\"is too long "
    "(maximum length is YAC_MAX_CHARLEN)", grid_name)

  struct bnd_sphere_part_search * search = NULL;

  for (int i = 0; (i < 2) && (search == NULL); ++i)
    if (!strcmp(grid_name, grid_pair->grid_names[i]))
      search = grid_pair->cell_sphere_part[i];
  YAC_ASSERT(
    search != NULL,
    "ERROR(yac_dist_grid_pair_get_cell_sphere_part): invalid grid_name")
  return search;
}

static void dist_grid_pair_do_point_search_local(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  coordinate_pointer search_coords, size_t count, size_t * cells,
  int (*coord_in_cell)(
    double coord[3], struct dist_grid * dist_grid, size_t cell_idx,
    struct grid_cell * buffer_cell)) {

  struct bnd_sphere_part_search * cell_sphere_part =
    dist_grid_pair_get_cell_sphere_part(grid_pair, grid_name);
  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name);

  size_t * temp_cells;
  size_t * num_cells_per_coord =
    xmalloc(count * sizeof(*num_cells_per_coord));

  // search for all matching source cells
  yac_bnd_sphere_part_search_do_point_search(
    cell_sphere_part, search_coords, count, &temp_cells,
    num_cells_per_coord);

  struct grid_cell buffer_cell;
  yac_init_grid_cell(&buffer_cell);

  // if we have multiple source cells for a single search coordinate, get the
  // source cell with the lowest global id
  for (size_t i = 0, k = 0; i < count; ++i) {
    size_t curr_num_cells = num_cells_per_coord[i];
    if (curr_num_cells == 0) {
      cells[i] = SIZE_MAX;
    } else if (curr_num_cells == 1) {
      if (coord_in_cell(
            search_coords[i], dist_grid, temp_cells[k], &buffer_cell))
        cells[i] = temp_cells[k];
      else
        cells[i] = SIZE_MAX;
      ++k;
    } else {
      size_t cell_idx = SIZE_MAX;
      yac_int cell_id = XT_INT_MAX;
      for (size_t j = 0; j < curr_num_cells; ++j, ++k) {
        size_t curr_cell_idx = temp_cells[k];
        yac_int curr_cell_id = dist_grid->ids[CELL][curr_cell_idx];
        if (!coord_in_cell(
               search_coords[i], dist_grid, curr_cell_idx, &buffer_cell))
          continue;
        if (curr_cell_id < cell_id) {
          cell_idx = curr_cell_idx;
          cell_id = curr_cell_id;
        }
      }
      cells[i] = cell_idx;
    }
  }

  yac_free_grid_cell(&buffer_cell);
  free(num_cells_per_coord);
  free(temp_cells);
}

static void lookup_single_remote_point_reorder_locally(
  struct dist_grid * dist_grid, enum yac_location location,
  struct single_remote_point_reorder * ids, size_t * count, size_t * idx) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(lookup_single_remote_point_reorder_locally): invalid location")

  yac_int * sorted_ids = dist_grid->sorted_ids[location];
  size_t * reorder_idx = dist_grid->sorted_reorder_idx[location];
  size_t num_ids = dist_grid->total_count[location];

  size_t count_ = *count;
  size_t new_count = 0;

  // sort ids by global ids
  qsort(ids, count_, sizeof(*ids),
    compare_single_remote_point_reorder_global_id);

  for (size_t i = 0, j = 0; i < count_; ++i) {
    yac_int curr_id = ids[i].data.global_id;
    while ((j < num_ids) && (sorted_ids[j] < curr_id)) ++j;
    if ((j < num_ids) && (sorted_ids[j] == curr_id)) {
      idx[ids[i].reorder_idx] = reorder_idx[j];
    } else {
      if (i != new_count) ids[new_count] = ids[i];
      ++new_count;
    }
  }

  *count = new_count;
}

static int get_pack_size_field_data(
  struct yac_field_data field_data, MPI_Comm comm) {

  int pack_size_field_coord, pack_size_field_mask;

  yac_mpi_call(
    MPI_Pack_size(3, MPI_DOUBLE, comm, &pack_size_field_coord), comm);
  pack_size_field_coord *= (int)field_data.coordinates_count;

  yac_mpi_call(
    MPI_Pack_size(1, MPI_INT, comm, &pack_size_field_mask), comm);
  pack_size_field_mask *= (int)field_data.masks_count;

  return pack_size_field_coord + pack_size_field_mask;
}

static int get_pack_size_base_cell(
  struct yac_field_data cell_field_data,
  MPI_Datatype bnd_circle_dt, MPI_Comm comm) {

  int pack_size_id,
      pack_size_num_vertices,
      pack_size_bnd_circle;

  // id
  yac_mpi_call(MPI_Pack_size(1, yac_int_dt, comm, &pack_size_id), comm);
  // num_vertices
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &pack_size_num_vertices), comm);
  // bounding circle
  yac_mpi_call(
    MPI_Pack_size(1, bnd_circle_dt, comm, &pack_size_bnd_circle), comm);

  return pack_size_id + pack_size_num_vertices + pack_size_bnd_circle +
         get_pack_size_field_data(cell_field_data, comm);
}

static int get_pack_size_base_vertex(
  struct yac_field_data vertex_field_data, MPI_Comm comm) {

  int pack_size_id,
      pack_size_vertex_coords;
  // id
  yac_mpi_call(MPI_Pack_size(1, yac_int_dt, comm, &pack_size_id), comm);
  // vertex coordinates
  yac_mpi_call(
    MPI_Pack_size(3, MPI_DOUBLE, comm, &pack_size_vertex_coords), comm);

  return pack_size_id + pack_size_vertex_coords +
         get_pack_size_field_data(vertex_field_data, comm);
}

static int get_pack_size_base_edge(
  struct yac_field_data edge_field_data, MPI_Comm comm) {

  int pack_size_id,
      pack_size_edge_to_vertex,
      pack_size_edge_type;
  // id
  yac_mpi_call(MPI_Pack_size(1, yac_int_dt, comm, &pack_size_id), comm);
  // edge type
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &pack_size_edge_type), comm);
  // edge vertex ids
  yac_mpi_call(
    MPI_Pack_size(2, yac_int_dt, comm, &pack_size_edge_to_vertex), comm);

  return pack_size_id + pack_size_edge_type + pack_size_edge_to_vertex +
         get_pack_size_field_data(edge_field_data, comm);
}

static void get_pack_sizes_cell(
  struct dist_grid * dist_grid, uint64_t * pos, size_t count, int * pack_sizes,
  MPI_Datatype bnd_circle_dt, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int pack_size_base_cell =
    get_pack_size_base_cell(dist_grid->field_data.cell, bnd_circle_dt, comm);
  int pack_size_base_vertex =
    get_pack_size_base_vertex(dist_grid->field_data.vertex, comm);
  int pack_size_base_edge =
    get_pack_size_base_edge(dist_grid->field_data.edge, comm);

  for (size_t i = 0; i < count; ++i) {
    size_t idx = (size_t)(pos[i]);
    int num_vertices = dist_grid->num_vertices_per_cell[idx];
    size_t * curr_vertices =
      dist_grid->cell_to_vertex + dist_grid->cell_to_vertex_offsets[idx];
    size_t * curr_edges =
      dist_grid->cell_to_edge + dist_grid->cell_to_edge_offsets[idx];
    int pack_size =
      pack_size_base_cell +
      num_vertices * (pack_size_base_vertex + pack_size_base_edge) +
      yac_remote_point_infos_get_pack_size(
        dist_grid->owners[CELL] + idx, point_info_dt, comm);
    for (int j = 0; j < num_vertices; ++j) {
      pack_size +=
        yac_remote_point_infos_get_pack_size(
          dist_grid->owners[CORNER] + curr_vertices[j], point_info_dt, comm) +
        yac_remote_point_infos_get_pack_size(
          dist_grid->owners[EDGE] + curr_edges[j], point_info_dt, comm);
    }
    pack_sizes[i] = pack_size;
  }
}

static void get_pack_sizes_vertex(
  struct dist_grid * dist_grid, uint64_t * pos, size_t count, int * pack_sizes,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  int pack_size_base_vertex =
    get_pack_size_base_vertex(dist_grid->field_data.vertex, comm);
  for (size_t i = 0; i < count; ++i)
    pack_sizes[i] =
      pack_size_base_vertex +
      yac_remote_point_infos_get_pack_size(
        dist_grid->owners[CORNER] + pos[i], point_info_dt, comm);
}

static void get_pack_sizes_edge(
  struct dist_grid * dist_grid, uint64_t * pos, size_t count, int * pack_sizes,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  int pack_size_base_vertex =
    get_pack_size_base_vertex(dist_grid->field_data.vertex, comm);
  int pack_size_base_edge =
    get_pack_size_base_edge(dist_grid->field_data.edge, comm);
  for (size_t i = 0; i < count; ++i) {
    size_t * curr_vertices = dist_grid->edge_to_vertex[pos[i]];
    pack_sizes[i] =
      pack_size_base_edge +
      2 * pack_size_base_vertex +
      yac_remote_point_infos_get_pack_size(
        dist_grid->owners[EDGE] + pos[i], point_info_dt, comm) +
      yac_remote_point_infos_get_pack_size(
        dist_grid->owners[CORNER] + curr_vertices[0], point_info_dt, comm) +
      yac_remote_point_infos_get_pack_size(
        dist_grid->owners[CORNER] + curr_vertices[1], point_info_dt, comm);
  }
}

static void get_pack_sizes(
  struct dist_grid * dist_grid, enum yac_location location, uint64_t * pos,
  size_t count, int * pack_sizes, MPI_Datatype bnd_circle_dt,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(get_pack_sizes): invalid location")

  switch(location) {
    default:
    case(CELL):
      get_pack_sizes_cell(
        dist_grid, pos, count, pack_sizes, bnd_circle_dt, point_info_dt, comm);
      break;
    case(CORNER):
      get_pack_sizes_vertex(
        dist_grid, pos, count, pack_sizes, point_info_dt, comm);
      break;
    case(EDGE):
      get_pack_sizes_edge(
        dist_grid, pos, count, pack_sizes, point_info_dt, comm);
      break;
  };
}

static void pack_field_data(
  size_t idx, void * buffer, int buffer_size, int * position,
  struct yac_field_data field_data, MPI_Comm comm) {

  // coordinates
  for (size_t i = 0; i < field_data.coordinates_count; ++i)
    yac_mpi_call(
      MPI_Pack(field_data.coordinates[i][idx], 3, MPI_DOUBLE, buffer,
               buffer_size, position, comm), comm);

  // masks
  for (size_t i = 0; i < field_data.masks_count; ++i)
    yac_mpi_call(
      MPI_Pack(field_data.masks[i].data + idx, 1, MPI_INT, buffer,
               buffer_size, position, comm), comm);
}

static void pack_grid_data_vertex(
  struct dist_grid * dist_grid, size_t idx, void * buffer, int buffer_size,
  int * position, MPI_Datatype bnd_circle_dt, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  UNUSED(bnd_circle_dt);

  // id
  yac_mpi_call(
    MPI_Pack(dist_grid->ids[CORNER] + idx, 1, yac_int_dt, buffer, buffer_size,
             position, comm), comm);
  // vertex coordinates
  yac_mpi_call(
    MPI_Pack(&(dist_grid->vertex_coordinates[idx][0]), 3, MPI_DOUBLE, buffer,
             buffer_size, position, comm), comm);
  // vertex owner
  yac_remote_point_infos_pack_single(
    dist_grid->owners[CORNER] + idx, buffer, buffer_size, position,
    point_info_dt, comm);
  // pack field data
  pack_field_data(
    idx, buffer, buffer_size, position, dist_grid->field_data.vertex, comm);
}

static void pack_grid_data_edge_(
  struct dist_grid * dist_grid, size_t idx, void * buffer, int buffer_size,
  int * position, MPI_Datatype bnd_circle_dt, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  UNUSED(bnd_circle_dt);

  int edge_type = (int)(dist_grid->edge_type[idx]);

  // id
  yac_mpi_call(
    MPI_Pack(
      dist_grid->ids[EDGE] + idx, 1, yac_int_dt, buffer, buffer_size, position,
      comm), comm);
  // edge type
  yac_mpi_call(
    MPI_Pack(
      &edge_type, 1, MPI_INT, buffer, buffer_size, position, comm), comm);
  // edge to vertex
  yac_int edge_to_vertex[2] = {
    dist_grid->ids[CORNER][dist_grid->edge_to_vertex[idx][0]],
    dist_grid->ids[CORNER][dist_grid->edge_to_vertex[idx][1]]};
  yac_mpi_call(
    MPI_Pack(
      edge_to_vertex, 2, yac_int_dt, buffer, buffer_size, position, comm),
    comm);
  // edge owner
  yac_remote_point_infos_pack_single(
    dist_grid->owners[EDGE] + idx, buffer, buffer_size, position,
    point_info_dt, comm);
  // pack field data
  pack_field_data(
    idx, buffer, buffer_size, position, dist_grid->field_data.edge, comm);
}

static void pack_grid_data_edge(
  struct dist_grid * dist_grid, size_t idx, void * buffer, int buffer_size,
  int * position, MPI_Datatype bnd_circle_dt, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  pack_grid_data_edge_(
    dist_grid, idx, buffer, buffer_size, position,
    bnd_circle_dt, point_info_dt, comm);

  // pack edge vertices
  for (int i = 0; i < 2; ++i)
    pack_grid_data_vertex(
      dist_grid, dist_grid->edge_to_vertex[idx][i],
      buffer, buffer_size, position, bnd_circle_dt, point_info_dt, comm);
}

static void pack_grid_data_cell(
  struct dist_grid * dist_grid, size_t idx, void * buffer, int buffer_size,
  int * position, MPI_Datatype bnd_circle_dt, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  int num_vertices = dist_grid->num_vertices_per_cell[idx];

  // id
  yac_mpi_call(
    MPI_Pack(
      dist_grid->ids[CELL] + idx, 1, yac_int_dt, buffer, buffer_size, position,
      comm), comm);
  // pack field data
  pack_field_data(
    idx, buffer, buffer_size, position, dist_grid->field_data.cell, comm);
  // num_vertices
  yac_mpi_call(
    MPI_Pack(&num_vertices, 1, MPI_INT, buffer,
             buffer_size, position, comm), comm);
  // bounding_circle
  yac_mpi_call(
    MPI_Pack(dist_grid->cell_bnd_circles + idx, 1, bnd_circle_dt, buffer,
             buffer_size, position, comm), comm);
  // cell owner
  yac_remote_point_infos_pack_single(
    dist_grid->owners[CELL] + idx, buffer, buffer_size, position,
    point_info_dt, comm);

  for (int i = 0; i < num_vertices; ++i) {
    pack_grid_data_vertex(
      dist_grid,
      dist_grid->cell_to_vertex[dist_grid->cell_to_vertex_offsets[idx] + i],
      buffer, buffer_size, position, bnd_circle_dt, point_info_dt, comm);
    pack_grid_data_edge_(
      dist_grid,
      dist_grid->cell_to_edge[dist_grid->cell_to_edge_offsets[idx] + i],
      buffer, buffer_size, position, bnd_circle_dt, point_info_dt, comm);
  }
}

static void pack_grid_data(
  struct dist_grid * dist_grid, enum yac_location location, uint64_t * pos,
  size_t count, void ** pack_data, int * pack_sizes,
  MPI_Datatype bnd_circle_dt, MPI_Datatype point_info_dt, MPI_Comm comm) {

  get_pack_sizes(dist_grid, location, pos, count, pack_sizes,
                 bnd_circle_dt, point_info_dt, comm);

  size_t pack_size = 0;
  for (size_t i = 0; i < count; ++i) pack_size += (size_t)(pack_sizes[i]);

  void * pack_data_ = xmalloc(pack_size);

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(pack_grid_data): invalid location")

  void (*func_pack[3])(
    struct dist_grid * dist_grid, size_t idx, void * buffer, int buffer_size,
    int * position, MPI_Datatype bnd_circle_dt, MPI_Datatype point_info_dt,
    MPI_Comm comm) =
    {pack_grid_data_cell, pack_grid_data_vertex, pack_grid_data_edge};

  for (size_t i = 0, offset = 0; i < count; ++i) {
    int position = 0;
    func_pack[location](
      dist_grid, pos[i], (char*)pack_data_ + offset, pack_sizes[i],
      &position, bnd_circle_dt, point_info_dt, comm);
    pack_sizes[i] = position;
    offset += (size_t)position;
  }

  *pack_data = pack_data_;
}

static void unpack_field_data(
  void * buffer, int buffer_size, int * position, size_t idx,
  struct temp_field_data temp_field_data, MPI_Comm comm) {

  for (size_t i = 0; i < temp_field_data.coordinates_count; ++i)
    yac_mpi_call(
      MPI_Unpack(
        buffer, buffer_size, position, temp_field_data.coordinates[i][idx],
        3, MPI_DOUBLE, comm), comm);

  for (size_t i = 0; i < temp_field_data.masks_count; ++i)
    yac_mpi_call(
      MPI_Unpack(buffer, buffer_size, position,
                 temp_field_data.masks[i] + idx, 1, MPI_INT, comm), comm);
}

static void unpack_grid_data_vertex(
  struct global_vertex_reorder * vertex, size_t idx, void * buffer,
  int buffer_size, int * position,
  struct temp_field_data temp_vertex_field_data,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  // id
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, &(vertex[idx].global_id), 1,
               yac_int_dt, comm), comm);
  // vertex coordinates
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, &(vertex[idx].coord[0]), 3,
               MPI_DOUBLE, comm), comm);
  // vertex owners
  yac_remote_point_infos_unpack(
    buffer, buffer_size, position, &(vertex[idx].owners), point_info_dt, comm);
  // unpack field data
  unpack_field_data(
    buffer, buffer_size, position, idx, temp_vertex_field_data, comm);
}

static void unpack_grid_data_edge(
  struct global_edge_reorder * edge, size_t idx, void * buffer,
  int buffer_size, int * position,
  struct temp_field_data temp_edge_field_data,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  int edge_type;

  // id
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, &(edge[idx].global_id), 1,
               yac_int_dt, comm), comm);
  // edge type
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, &edge_type, 1,
               MPI_INT, comm), comm);
  edge[idx].edge_type = (enum yac_edge_type)edge_type;
  // edge to vertex
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, edge[idx].edge_to_vertex, 2,
               yac_int_dt, comm), comm);
  // edge owners
  yac_remote_point_infos_unpack(buffer, buffer_size, position,
                                &(edge[idx].owners), point_info_dt, comm);
  // unpack field data
  unpack_field_data(
    buffer, buffer_size, position, idx, temp_edge_field_data, comm);
}

static int compare_global_vertex_reorder_global_id(
  const void * a, const void * b) {

  return (((const struct global_vertex_reorder *)a)->global_id >
          ((const struct global_vertex_reorder *)b)->global_id) -
         (((const struct global_vertex_reorder *)a)->global_id <
          ((const struct global_vertex_reorder *)b)->global_id);
}

static void add_field_data(
  struct yac_field_data * field_data, struct temp_field_data temp_field_data,
  void * reorder_idx, size_t reorder_idx_size,
  size_t old_count, size_t new_count) {

  size_t add_count = new_count - old_count;

  for (size_t i = 0; i < temp_field_data.masks_count; ++i) {
    int * temp_mask = temp_field_data.masks[i];
    int * mask =
      ((field_data->masks[i].data =
         xrealloc(field_data->masks[i].data, new_count * sizeof(*mask))));
    for (size_t i = 0, j = old_count; i < add_count; ++i, ++j) {
      size_t idx =
        *(size_t*)((unsigned char*)reorder_idx + i * reorder_idx_size);
      mask[j] = temp_mask[idx];
    }
  }

  for (size_t i = 0; i < temp_field_data.coordinates_count; ++i) {
    coordinate_pointer temp_coordinates = temp_field_data.coordinates[i];
    coordinate_pointer coordinates =
      ((field_data->coordinates[i] =
         xrealloc(
          field_data->coordinates[i], new_count * sizeof(*coordinates))));
    for (size_t i = 0, j = old_count; i < add_count; ++i, ++j) {
      size_t idx =
        *(size_t*)((unsigned char*)reorder_idx + i * reorder_idx_size);
      coordinates[j][0] = temp_coordinates[idx][0];
      coordinates[j][1] = temp_coordinates[idx][1];
      coordinates[j][2] = temp_coordinates[idx][2];
    }
  }
}

static void yac_remote_point_infos_single_free(
  struct remote_point_infos * point_infos) {

  if (point_infos->count > 1) free(point_infos->data.multi);
}

static void yac_dist_grid_add_vertices(
  struct dist_grid * dist_grid, struct global_vertex_reorder * vertices,
  size_t count, size_t * idx,
  struct temp_field_data temp_vertex_field_data) {

  if (count == 0) return;

  // sort vertices global ids
  qsort(vertices, count, sizeof(*vertices),
        compare_global_vertex_reorder_global_id);

  yac_int * sorted_vertex_ids = dist_grid->sorted_ids[CORNER];
  size_t * sorted_vertex_reorder_idx = dist_grid->sorted_reorder_idx[CORNER];

  yac_int prev_global_id = vertices[0].global_id - 1;
  size_t prev_idx = 0;
  size_t add_count = 0;
  size_t num_total_vertices = dist_grid->total_count[CORNER];

  // determine which vertices need to be added to local data
  for (size_t i = 0, j = 0; i < count; ++i) {

    yac_int curr_global_id = vertices[i].global_id;
    size_t curr_reorder_idx = vertices[i].reorder_idx;

    // if the current global id is a duplicate
    if (prev_global_id == curr_global_id) {
      if (idx != NULL) idx[curr_reorder_idx] = prev_idx;
      yac_remote_point_infos_single_free(&(vertices[i].owners));
      continue;
    }
    prev_global_id = curr_global_id;

    // check whether the current global id is already part of the local
    // grid data
    while ((j < num_total_vertices) && (sorted_vertex_ids[j] < curr_global_id))
      ++j;

    // if we found a match in the local data
    if ((j < num_total_vertices) && (sorted_vertex_ids[j] == curr_global_id)) {

      if (idx != NULL) idx[curr_reorder_idx] = sorted_vertex_reorder_idx[j];
      prev_idx = sorted_vertex_reorder_idx[j];
      yac_remote_point_infos_single_free(&(vertices[i].owners));

    // if we need to add the current vertex to the local data
    } else {

      if (idx != NULL) idx[curr_reorder_idx] = num_total_vertices + add_count;
      prev_idx = num_total_vertices + add_count;
      if (add_count != i) vertices[add_count] = vertices[i];
      ++add_count;
    }
  }

  size_t new_num_total_vertices = num_total_vertices + add_count;
  coordinate_pointer vertex_coordinates =
    xrealloc(dist_grid->vertex_coordinates, new_num_total_vertices *
             sizeof(*vertex_coordinates));
  yac_int * vertex_ids =
    xrealloc(dist_grid->ids[CORNER], new_num_total_vertices * sizeof(*vertex_ids));
  int * core_vertex_mask =
    xrealloc(dist_grid->core_mask[CORNER], new_num_total_vertices *
             sizeof(*core_vertex_mask));
  struct remote_point_infos * vertex_owners =
    xrealloc(dist_grid->owners[CORNER], new_num_total_vertices *
             sizeof(*vertex_owners));
  sorted_vertex_ids =
    xrealloc(
      sorted_vertex_ids, new_num_total_vertices * sizeof(*sorted_vertex_ids));
  sorted_vertex_reorder_idx =
    xrealloc(
      sorted_vertex_reorder_idx, new_num_total_vertices *
      sizeof(*sorted_vertex_reorder_idx));

  // add the selected vertices to the local grid data
  for (size_t i = 0, j = num_total_vertices; i < add_count; ++i, ++j) {

    vertex_coordinates[j][0] = vertices[i].coord[0];
    vertex_coordinates[j][1] = vertices[i].coord[1];
    vertex_coordinates[j][2] = vertices[i].coord[2];
    vertex_ids[j] = vertices[i].global_id;
    core_vertex_mask[j] = 0;
    vertex_owners[j] = vertices[i].owners;
    sorted_vertex_ids[j] = vertices[i].global_id;
    sorted_vertex_reorder_idx[j] = j;
  }
  // add field data
  add_field_data(
    &(dist_grid->field_data.vertex), temp_vertex_field_data,
    vertices, sizeof(*vertices),
    num_total_vertices, new_num_total_vertices);
  yac_quicksort_index_yac_int_size_t(
    sorted_vertex_ids, new_num_total_vertices, sorted_vertex_reorder_idx);

  dist_grid->vertex_coordinates = vertex_coordinates;
  dist_grid->ids[CORNER] = vertex_ids;
  dist_grid->core_mask[CORNER] = core_vertex_mask;
  dist_grid->owners[CORNER] = vertex_owners;
  dist_grid->sorted_ids[CORNER] = sorted_vertex_ids;
  dist_grid->sorted_reorder_idx[CORNER] = sorted_vertex_reorder_idx;
  dist_grid->total_count[CORNER] = new_num_total_vertices;
}

static int compare_global_edge_reorder_global_id(
  const void * a, const void * b) {

  return (((const struct global_edge_reorder *)a)->global_id >
          ((const struct global_edge_reorder *)b)->global_id) -
         (((const struct global_edge_reorder *)a)->global_id <
          ((const struct global_edge_reorder *)b)->global_id);
}

static void yac_dist_grid_add_edges(
  struct dist_grid * dist_grid, struct global_edge_reorder * edges,
  size_t count, size_t * idx, struct temp_field_data temp_edge_field_data) {

  if (count == 0) return;

  // sort edges global ids
  qsort(edges, count, sizeof(*edges), compare_global_edge_reorder_global_id);

  yac_int * sorted_edge_ids = dist_grid->sorted_ids[EDGE];
  size_t * sorted_edge_reorder_idx = dist_grid->sorted_reorder_idx[EDGE];

  yac_int prev_global_id = edges[0].global_id - 1;
  size_t prev_idx = 0;
  size_t add_count = 0;
  size_t num_total_edges = dist_grid->total_count[EDGE];

  // determine which edges need to be added to local data
  for (size_t i = 0, j = 0; i < count; ++i) {

    yac_int curr_global_id = edges[i].global_id;
    size_t curr_reorder_idx = edges[i].reorder_idx;

    // if the current global id is a duplicate
    if (prev_global_id == curr_global_id) {
      if (idx != NULL) idx[curr_reorder_idx] = prev_idx;
      yac_remote_point_infos_single_free(&(edges[i].owners));
      continue;
    }
    prev_global_id = curr_global_id;

    // check whether the current global id is already part of the local
    // grid data
    while ((j < num_total_edges) && (sorted_edge_ids[j] < curr_global_id)) ++j;

    // if we found a match in the local data
    if ((j < num_total_edges) && (sorted_edge_ids[j] == curr_global_id)) {

      if (idx != NULL) idx[curr_reorder_idx] = sorted_edge_reorder_idx[j];
      prev_idx = sorted_edge_reorder_idx[j];
      yac_remote_point_infos_single_free(&(edges[i].owners));

    // if we need to add the current edge to the local data
    } else {

      if (idx != NULL) idx[curr_reorder_idx] = num_total_edges + add_count;
      prev_idx = num_total_edges + add_count;
      if (add_count != i) edges[add_count] = edges[i];
      ++add_count;
    }
  }

  size_t new_num_total_edges = num_total_edges + add_count;
  yac_int * edge_ids =
    xrealloc(dist_grid->ids[EDGE], new_num_total_edges * sizeof(*edge_ids));
  enum yac_edge_type * edge_type =
    xrealloc(dist_grid->edge_type, new_num_total_edges * sizeof(*edge_type));
  size_t_2_pointer edge_to_vertex =
    xrealloc(dist_grid->edge_to_vertex, new_num_total_edges * sizeof(*edge_to_vertex));
  struct remote_point_infos * edge_owners =
    xrealloc(dist_grid->owners[EDGE], new_num_total_edges * sizeof(*edge_owners));
  int * core_edge_mask =
    xrealloc(dist_grid->core_mask[EDGE], new_num_total_edges *
             sizeof(*core_edge_mask));
  sorted_edge_ids =
    xrealloc(
      sorted_edge_ids, new_num_total_edges * sizeof(*sorted_edge_ids));
  sorted_edge_reorder_idx =
    xrealloc(
      sorted_edge_reorder_idx, new_num_total_edges *
      sizeof(*sorted_edge_reorder_idx));

  yac_int * vertex_ids = xmalloc(2 * add_count * sizeof(*vertex_ids));
  size_t * reorder = xmalloc(2 * add_count * sizeof(*reorder));

  // add the selected edges to the local grid data
  for (size_t i = 0, j = num_total_edges; i < add_count; ++i, ++j) {

    edge_ids[j] = edges[i].global_id;
    edge_type[j] = edges[i].edge_type;
    core_edge_mask[j] = 0;
    edge_owners[j] = edges[i].owners;
    sorted_edge_ids[j] = edges[i].global_id;
    sorted_edge_reorder_idx[j] = j;

    vertex_ids[2 * i + 0] = edges[i].edge_to_vertex[0];
    vertex_ids[2 * i + 1] = edges[i].edge_to_vertex[1];
    reorder[2 * i + 0] = 2 * num_total_edges + 2 * i + 0;
    reorder[2 * i + 1] = 2 * num_total_edges + 2 * i + 1;
  }
  // add field data
  add_field_data(
    &(dist_grid->field_data.edge), temp_edge_field_data, edges, sizeof(*edges),
    num_total_edges, new_num_total_edges);
  yac_quicksort_index_yac_int_size_t(
    sorted_edge_ids, new_num_total_edges, sorted_edge_reorder_idx);

  { // determine vertex indices for edge_to_vertex
    yac_quicksort_index_yac_int_size_t(vertex_ids, 2 * add_count, reorder);
    yac_int * sorted_vertex_ids = dist_grid->sorted_ids[CORNER];
    size_t * sorted_vertex_reorder_idx = dist_grid->sorted_reorder_idx[CORNER];
    size_t total_num_vertices = dist_grid->total_count[CORNER];
    size_t * edge_to_vertex_ = (size_t*)&(edge_to_vertex[0][0]);
    // lookup global ids
    for (size_t i = 0, j = 0; i < 2 * add_count; ++i) {
      yac_int curr_id = vertex_ids[i];
      while ((j < total_num_vertices) && (sorted_vertex_ids[j] < curr_id)) ++j;
      YAC_ASSERT(
        (j < total_num_vertices) && (sorted_vertex_ids[j] == curr_id),
        "ERROR(yac_dist_grid_add_edges): vertex id not found")
      edge_to_vertex_[reorder[i]] = sorted_vertex_reorder_idx[j];
    }
  }

  free(vertex_ids);
  free(reorder);

  dist_grid->ids[EDGE] = edge_ids;
  dist_grid->edge_type = edge_type;
  dist_grid->edge_to_vertex = edge_to_vertex;
  dist_grid->owners[EDGE] = edge_owners;
  dist_grid->core_mask[EDGE] = core_edge_mask;
  dist_grid->sorted_ids[EDGE] = sorted_edge_ids;
  dist_grid->sorted_reorder_idx[EDGE] = sorted_edge_reorder_idx;
  dist_grid->total_count[EDGE] = new_num_total_edges;
}

static void yac_dist_grid_add_cells(
  struct dist_grid * dist_grid, yac_int * cell_ids, int * num_vertices_per_cell,
  struct bounding_circle * cell_bnd_circles, size_t count,
  size_t * cell_to_vertex, size_t * cell_to_edge,
  struct remote_point_infos * cell_owners,
  struct temp_field_data temp_cell_field_data) {

  if (count == 0) return;

  size_t * reorder_idx = xmalloc(count * sizeof(reorder_idx));
  for (size_t i = 0; i < count; ++i) reorder_idx[i] = i;

  size_t * prescan = xmalloc(count * sizeof(*prescan));
  for (size_t i = 0, accu = 0; i < count;
       accu += (size_t)(num_vertices_per_cell[i++])) prescan[i] = accu;

  // sort cells global ids
  yac_quicksort_index_yac_int_size_t(cell_ids, count, reorder_idx);

  yac_int * sorted_cell_ids = dist_grid->sorted_ids[CELL];
  size_t * sorted_cell_reorder_idx = dist_grid->sorted_reorder_idx[CELL];

  yac_int prev_global_id = cell_ids[0] - 1;
  size_t cell_add_count = 0;
  size_t relations_add_count = 0;
  size_t num_total_cells = dist_grid->total_count[CELL];

  // determine which cells need to be added to local data
  for (size_t i = 0, j = 0; i < count; ++i) {

    yac_int curr_global_id = cell_ids[i];
    size_t curr_reorder_idx = reorder_idx[i];

    // if the current global id is a duplicate
    if (prev_global_id == curr_global_id) {
      yac_remote_point_infos_single_free(cell_owners + curr_reorder_idx);
      continue;
    }
    prev_global_id = curr_global_id;

    // check whether the current global id is already part of the local
    // grid data
    while ((j < num_total_cells) && (sorted_cell_ids[j] < curr_global_id)) ++j;

    // if we did not find a match in the local data
    if ((j >= num_total_cells) || (sorted_cell_ids[j] != curr_global_id)) {

      if (cell_add_count != i) {
        cell_ids[cell_add_count] = curr_global_id;
        reorder_idx[cell_add_count] = curr_reorder_idx;
      }
      ++cell_add_count;
      relations_add_count += (size_t)(num_vertices_per_cell[curr_reorder_idx]);
    }
  }

  size_t new_num_total_cells = num_total_cells + cell_add_count;
  size_t num_total_relations =
    (num_total_cells > 0)?
      (dist_grid->cell_to_vertex_offsets[num_total_cells-1] +
       (size_t)(dist_grid->num_vertices_per_cell[num_total_cells-1])):0;
  size_t new_num_total_relations = num_total_relations + relations_add_count;
  yac_int * new_cell_ids =
    xrealloc(dist_grid->ids[CELL], new_num_total_cells * sizeof(*new_cell_ids));
  int * new_num_vertices_per_cell =
    xrealloc(dist_grid->num_vertices_per_cell, new_num_total_cells *
             sizeof(*new_num_vertices_per_cell));
  size_t * new_cell_to_vertex =
    xrealloc(dist_grid->cell_to_vertex, new_num_total_relations *
             sizeof(*new_cell_to_vertex));
  size_t * cell_to_vertex_offsets =
    xrealloc(dist_grid->cell_to_vertex_offsets, new_num_total_cells *
             sizeof(*cell_to_vertex_offsets));
  size_t * new_cell_to_edge =
    xrealloc(dist_grid->cell_to_edge, new_num_total_relations *
             sizeof(*new_cell_to_edge));
  struct bounding_circle * new_cell_bnd_circles =
    xrealloc(dist_grid->cell_bnd_circles, new_num_total_cells *
             sizeof(*new_cell_bnd_circles));
  int * core_cell_mask =
    xrealloc(dist_grid->core_mask[CELL], new_num_total_cells * sizeof(*core_cell_mask));
  struct remote_point_infos * new_cell_owners =
    xrealloc(dist_grid->owners[CELL], new_num_total_cells * sizeof(*cell_owners));
  sorted_cell_ids =
    xrealloc(
      sorted_cell_ids, new_num_total_cells * sizeof(*sorted_cell_ids));
  sorted_cell_reorder_idx =
    xrealloc(
      sorted_cell_reorder_idx, new_num_total_cells *
      sizeof(*sorted_cell_reorder_idx));

  // add the selected cells to the local grid data
  for (size_t i = 0, j = num_total_cells; i < cell_add_count;
       ++i, ++j) {

    size_t curr_reorder_idx = reorder_idx[i];
    int curr_num_vertices = num_vertices_per_cell[curr_reorder_idx];
    size_t curr_relation_idx = prescan[curr_reorder_idx];

    new_cell_ids[j] = cell_ids[i];
    new_num_vertices_per_cell[j] = curr_num_vertices;
    cell_to_vertex_offsets[j] = num_total_relations;
    for (int j = 0; j < curr_num_vertices;
         ++j, ++num_total_relations, ++curr_relation_idx) {
      new_cell_to_vertex[num_total_relations] =
        cell_to_vertex[curr_relation_idx];
      new_cell_to_edge[num_total_relations] = cell_to_edge[curr_relation_idx];
    }
    core_cell_mask[j] = 0;
    sorted_cell_ids[j] = cell_ids[i];
    sorted_cell_reorder_idx[j] = j;
    new_cell_bnd_circles[j] = cell_bnd_circles[curr_reorder_idx];
    new_cell_owners[j] = cell_owners[curr_reorder_idx];
  }
  // add field data
  add_field_data(
    &(dist_grid->field_data.cell), temp_cell_field_data,
    reorder_idx, sizeof(*reorder_idx),
    num_total_cells, new_num_total_cells);
  yac_quicksort_index_yac_int_size_t(
    sorted_cell_ids, new_num_total_cells, sorted_cell_reorder_idx);

  dist_grid->ids[CELL] = new_cell_ids;
  dist_grid->num_vertices_per_cell = new_num_vertices_per_cell;
  dist_grid->cell_to_vertex = new_cell_to_vertex;
  dist_grid->cell_to_vertex_offsets = cell_to_vertex_offsets;
  dist_grid->cell_to_edge = new_cell_to_edge;
  dist_grid->cell_to_edge_offsets = cell_to_vertex_offsets;
  dist_grid->core_mask[CELL] = core_cell_mask;
  dist_grid->owners[CELL] = new_cell_owners;
  dist_grid->sorted_ids[CELL] = sorted_cell_ids;
  dist_grid->sorted_reorder_idx[CELL] = sorted_cell_reorder_idx;
  dist_grid->cell_bnd_circles = new_cell_bnd_circles;
  dist_grid->total_count[CELL] = new_num_total_cells;

  free(prescan);
  free(reorder_idx);
}

static void ensure_temp_field_data_sizes(
  struct temp_field_data * temp_field_data, size_t size) {

  for (size_t i = 0; i < temp_field_data->masks_count; ++i)
    ENSURE_ARRAY_SIZE(
      temp_field_data->masks[i], temp_field_data->masks_array_sizes[i], size);
  for (size_t i = 0; i < temp_field_data->coordinates_count; ++i)
    ENSURE_ARRAY_SIZE(
      temp_field_data->coordinates[i],
      temp_field_data->coordinates_array_sizes[i], size);
}

static struct temp_field_data temp_field_data_init(
  struct yac_field_data field_data, size_t count) {

  struct temp_field_data temp_field_data;

  temp_field_data.masks =
    xmalloc(field_data.masks_count * sizeof(*temp_field_data.masks));
  temp_field_data.masks_array_sizes =
    xmalloc(
      field_data.masks_count * sizeof(*temp_field_data.masks_array_sizes));
  temp_field_data.masks_count = field_data.masks_count;
  for (size_t i = 0; i < field_data.masks_count; ++i) {
    temp_field_data.masks[i] =
      xmalloc(count * sizeof(**temp_field_data.masks));
    temp_field_data.masks_array_sizes[i] = count;
  }

  temp_field_data.coordinates =
    xmalloc(
      field_data.coordinates_count * sizeof(*temp_field_data.coordinates));
  temp_field_data.coordinates_array_sizes =
    xmalloc(
      field_data.coordinates_count *
      sizeof(*temp_field_data.coordinates_array_sizes));
  temp_field_data.coordinates_count = field_data.coordinates_count;
  for (size_t i = 0; i < field_data.coordinates_count; ++i) {
    temp_field_data.coordinates[i] =
      xmalloc(count * sizeof(**temp_field_data.coordinates));
    temp_field_data.coordinates_array_sizes[i] = count;
  }

  return temp_field_data;
}

static void temp_field_data_free(struct temp_field_data temp_field_data) {

  for (size_t i = 0; i < temp_field_data.masks_count; ++i)
    free(temp_field_data.masks[i]);
  free(temp_field_data.masks);
  free(temp_field_data.masks_array_sizes);
  for (size_t i = 0; i < temp_field_data.coordinates_count; ++i)
    free(temp_field_data.coordinates[i]);
  free(temp_field_data.coordinates);
  free(temp_field_data.coordinates_array_sizes);
}

static void unpack_grid_data_cells(
  struct dist_grid * dist_grid, size_t count, void * buffer, int buffer_size,
  MPI_Datatype bnd_circle_dt, MPI_Datatype point_info_dt, MPI_Comm comm) {

  yac_int * cell_ids = xmalloc(count * sizeof(*cell_ids));
  int * num_vertices_per_cell = xmalloc(count * sizeof(*num_vertices_per_cell));
  struct bounding_circle * cell_bnd_circles =
    xmalloc(count * sizeof(*cell_bnd_circles));
  struct remote_point_infos * cell_owners =
    xmalloc(count * sizeof(*cell_owners));

  struct global_vertex_reorder * vertices = NULL;
  size_t vertices_array_size = 0;
  size_t total_num_vertices = 0;

  struct global_edge_reorder * edges = NULL;
  size_t edges_array_size = 0;

  struct temp_field_data temp_cell_field_data =
    temp_field_data_init(dist_grid->field_data.cell, count);
  struct temp_field_data temp_vertex_field_data =
    temp_field_data_init(dist_grid->field_data.vertex, 3 * count);
  struct temp_field_data temp_edge_field_data =
    temp_field_data_init(dist_grid->field_data.edge, 3 * count);

  for (size_t i = 0, buffer_offset = 0; i < count; ++i) {

    int position = 0;
    void * curr_buffer = (char*)buffer + buffer_offset;
    int num_vertices;

    // cell id
    yac_mpi_call(
      MPI_Unpack(curr_buffer, buffer_size, &position, cell_ids + i, 1,
                 yac_int_dt, comm), comm);
    // unpack field data
    unpack_field_data(
      curr_buffer, buffer_size, &position, i, temp_cell_field_data, comm);
    // num vertices
    yac_mpi_call(
      MPI_Unpack(curr_buffer, buffer_size, &position, &num_vertices, 1,
                 MPI_INT, comm), comm);
    // bounding circle
    yac_mpi_call(
      MPI_Unpack(curr_buffer, buffer_size, &position, cell_bnd_circles + i, 1,
                 bnd_circle_dt, comm), comm);
    // cell owners
    yac_remote_point_infos_unpack(
      curr_buffer, buffer_size, &position, cell_owners + i,
      point_info_dt, comm);

    num_vertices_per_cell[i] = num_vertices;

    ENSURE_ARRAY_SIZE(
      vertices, vertices_array_size, total_num_vertices + (size_t)num_vertices);
    ENSURE_ARRAY_SIZE(
      edges, edges_array_size, total_num_vertices + (size_t)num_vertices);
    ensure_temp_field_data_sizes(
      &temp_vertex_field_data, total_num_vertices + (size_t)num_vertices);
    ensure_temp_field_data_sizes(
      &temp_edge_field_data, total_num_vertices + (size_t)num_vertices);

    for (int j = 0; j < num_vertices; ++j, ++total_num_vertices) {
      unpack_grid_data_vertex(
        vertices, total_num_vertices, curr_buffer, buffer_size, &position,
        temp_vertex_field_data, point_info_dt, comm);
      unpack_grid_data_edge(
        edges, total_num_vertices, curr_buffer, buffer_size, &position,
        temp_edge_field_data, point_info_dt, comm);
      vertices[total_num_vertices].reorder_idx = total_num_vertices;
      edges[total_num_vertices].reorder_idx = total_num_vertices;
    }

    buffer_offset += (size_t)position;
    buffer_size -= position;
  }

  size_t * cell_to_vertex = xmalloc(total_num_vertices * sizeof(*cell_to_vertex));
  size_t * cell_to_edge = xmalloc(total_num_vertices * sizeof(*cell_to_edge));

  yac_dist_grid_add_vertices(
    dist_grid, vertices, total_num_vertices, cell_to_vertex,
    temp_vertex_field_data);
  yac_dist_grid_add_edges(
    dist_grid, edges, total_num_vertices, cell_to_edge,
    temp_edge_field_data);
  yac_dist_grid_add_cells(
    dist_grid, cell_ids, num_vertices_per_cell, cell_bnd_circles, count,
    cell_to_vertex, cell_to_edge, cell_owners, temp_cell_field_data);

  temp_field_data_free(temp_cell_field_data);
  temp_field_data_free(temp_vertex_field_data);
  temp_field_data_free(temp_edge_field_data);
  free(cell_to_edge);
  free(cell_to_vertex);
  free(vertices);
  free(edges);
  free(cell_owners);
  free(cell_bnd_circles);
  free(num_vertices_per_cell);
  free(cell_ids);
}

static void unpack_grid_data_vertices(
  struct dist_grid * dist_grid, size_t count, void * buffer, int buffer_size,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  struct global_vertex_reorder * vertices = xmalloc(count * sizeof(*vertices));

  struct temp_field_data temp_vertex_field_data =
    temp_field_data_init(dist_grid->field_data.vertex, count);

  for (size_t i = 0, buffer_offset = 0; i < count; ++i) {

    int position = 0;
    void * curr_buffer = (char*)buffer + buffer_offset;

    unpack_grid_data_vertex(
      vertices, i, curr_buffer, buffer_size, &position,
      temp_vertex_field_data, point_info_dt, comm);
    vertices[i].reorder_idx = i;

    buffer_offset += (size_t)position;
    buffer_size -= position;
  }

  yac_dist_grid_add_vertices(
    dist_grid, vertices, count, NULL, temp_vertex_field_data);

  temp_field_data_free(temp_vertex_field_data);

  free(vertices);
}

static void unpack_grid_data_edges(
  struct dist_grid * dist_grid, size_t count, void * buffer, int buffer_size,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  struct global_edge_reorder * edges = xmalloc(count * sizeof(*edges));
  struct global_vertex_reorder * vertices =
    xmalloc(2 * count * sizeof(*vertices));

  struct temp_field_data temp_edge_field_data =
    temp_field_data_init(dist_grid->field_data.edge, count);
  struct temp_field_data temp_vertex_field_data =
    temp_field_data_init(dist_grid->field_data.vertex, 2 * count);

  for (size_t i = 0, buffer_offset = 0; i < count; ++i) {

    int position = 0;
    void * curr_buffer = (char*)buffer + buffer_offset;

    unpack_grid_data_edge(
      edges, i, curr_buffer, buffer_size, &position,
      temp_edge_field_data, point_info_dt, comm);
    edges[i].reorder_idx = i;

    for (size_t j = 0; j < 2; ++j)
      unpack_grid_data_vertex(
        vertices, 2 * i + j, curr_buffer, buffer_size, &position,
        temp_vertex_field_data, point_info_dt, comm);

    buffer_offset += (size_t)position;
    buffer_size -= position;
  }

  yac_dist_grid_add_vertices(
    dist_grid, vertices, 2 * count, NULL, temp_vertex_field_data);
  yac_dist_grid_add_edges(
    dist_grid, edges, count, NULL, temp_edge_field_data);

  temp_field_data_free(temp_vertex_field_data);
  temp_field_data_free(temp_edge_field_data);

  free(vertices);
  free(edges);
}

static void unpack_grid_data(
  struct dist_grid * dist_grid, enum yac_location location, size_t count,
  void * buffer, int buffer_size, MPI_Datatype bnd_circle_dt,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(unpack_grid_data): invalid location")

  switch(location) {
    default:
    case(CELL):
      unpack_grid_data_cells(
        dist_grid, count, buffer, buffer_size, bnd_circle_dt,
        point_info_dt, comm);
      break;
    case(CORNER):
      unpack_grid_data_vertices(
        dist_grid, count, buffer, buffer_size, point_info_dt, comm);
      break;
    case(EDGE):
      unpack_grid_data_edges(
        dist_grid, count, buffer, buffer_size, point_info_dt, comm);
      break;
  };
}

static int compare_single_remote_point_reorder_owner(
  const void * a, const void * b) {

  return ((const struct single_remote_point_reorder *)a)->data.data.rank -
         ((const struct single_remote_point_reorder *)b)->data.data.rank;
}

static void yac_dist_grid_single_remote_point_to_local(
  struct dist_grid * dist_grid, struct single_remote_point * ids, size_t count,
  enum yac_location location, size_t * idx) {

  MPI_Comm comm = dist_grid->comm;
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t remote_count = 0;

  for (size_t i = 0; i < count; ++i) {
    if (ids[i].global_id == XT_INT_MAX) idx[i] = SIZE_MAX;
    else if (ids[i].data.rank != comm_rank) ++remote_count;
    else idx[i] = ids[i].data.orig_pos;
  }

  struct single_remote_point_reorder * missing_ids =
    xmalloc(remote_count * sizeof(*missing_ids));

  for (size_t i = 0, j = 0; i < count; ++i) {
    if ((ids[i].data.rank != comm_rank) &&
        (ids[i].global_id != XT_INT_MAX)) {
      missing_ids[j].data = ids[i];
      missing_ids[j].reorder_idx = i;
      ++j;
    }
  }

  // check whether we already have some of the missing ids locally
  lookup_single_remote_point_reorder_locally(
    dist_grid, location, missing_ids, &remote_count, idx);

  // sort data by owner
  qsort(missing_ids, remote_count, sizeof(*missing_ids),
        compare_single_remote_point_reorder_owner);

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0; i < remote_count; ++i)
    sendcounts[missing_ids[i].data.data.rank]++;

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t recv_count = rdispls[comm_size-1] + recvcounts[comm_size-1];

  uint64_t * uint64_t_buffer =
    xmalloc((remote_count + recv_count) * sizeof(*uint64_t_buffer));
  uint64_t * orig_pos_send_buffer = uint64_t_buffer;
  uint64_t * orig_pos_recv_buffer = uint64_t_buffer + remote_count;

  // pack send buffer
  for (size_t i = 0; i < remote_count; ++i) {
    int rank = missing_ids[i].data.data.rank;
    if (rank != comm_rank)
      orig_pos_send_buffer[sdispls[rank+1]++] =
        (uint64_t)(missing_ids[i].data.data.orig_pos);
  }

  // redistribute ids
  yac_alltoallv_uint64_p2p(
    orig_pos_send_buffer, sendcounts, sdispls,
    orig_pos_recv_buffer, recvcounts, rdispls, comm);

  MPI_Datatype bnd_circle_dt = yac_get_bounding_circle_mpi_datatype(comm);
  yac_mpi_call(MPI_Type_commit(&bnd_circle_dt), comm);
  MPI_Datatype point_info_dt = yac_get_remote_point_info_mpi_datatype(comm);
  yac_mpi_call(MPI_Type_commit(&point_info_dt), comm);

  void * packed_send_data = NULL;
  int * pack_sizes = xmalloc(recv_count * sizeof(*pack_sizes));

  // pack all requested grid data
  pack_grid_data(
    dist_grid, location, orig_pos_recv_buffer, recv_count,
    &packed_send_data, pack_sizes,
    bnd_circle_dt, point_info_dt, comm);
  free(uint64_t_buffer);

  memset(sendcounts, 0, (size_t)comm_size * sizeof(*sendcounts));
  for (int i = 0, k = 0; i < comm_size; ++i)
    for (size_t j = 0; j < recvcounts[i]; ++j, ++k)
      sendcounts[i] += (size_t)(pack_sizes[k]);

  free(pack_sizes);

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  recv_count = rdispls[comm_size-1] + recvcounts[comm_size-1];

  void * packed_recv_data = xmalloc(recv_count);

  // redistribute packed grid data
  yac_alltoallv_packed_p2p(
    packed_send_data, sendcounts, sdispls+1,
    packed_recv_data, recvcounts, rdispls, comm);

  // unpack requested grid data
  unpack_grid_data(
    dist_grid, location, remote_count, packed_recv_data, (int)recv_count,
    bnd_circle_dt, point_info_dt, comm);

  yac_mpi_call(MPI_Type_free(&point_info_dt), comm);
  yac_mpi_call(MPI_Type_free(&bnd_circle_dt), comm);

  // get the local ids for the remaining missing ids
  lookup_single_remote_point_reorder_locally(
    dist_grid, location, missing_ids, &remote_count, idx);

  free(missing_ids);
  free(packed_recv_data);
  free(packed_send_data);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
}

void yac_dist_grid_pair_do_point_search_(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  coordinate_pointer search_coords, size_t count, size_t * cells,
  int (*coord_in_cell)(
    double coord[3], struct dist_grid * dist_grid, size_t cell_idx,
    struct grid_cell * buffer_cell)) {

  MPI_Comm comm = grid_pair->comm;
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int * ranks = xmalloc(count * sizeof(ranks));

  // search for every point in the process tree
  yac_proc_sphere_part_do_point_search(
    grid_pair->proc_sphere_part, search_coords, count, ranks);

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);
  for (size_t i = 0; i < count; ++i) sendcounts[ranks[i]]++;

  size_t local_count = sendcounts[comm_rank];
  sendcounts[comm_rank] = 0;

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t remote_count = sdispls[comm_size] + sendcounts[comm_size-1];
  size_t request_count = rdispls[comm_size-1] + recvcounts[comm_size-1];

  coordinate_pointer coord_buffer =
    xmalloc((remote_count + request_count + local_count) *
            sizeof(*coord_buffer));
  coordinate_pointer coord_send_buffer = coord_buffer + 0;
  coordinate_pointer coord_recv_buffer = coord_buffer + remote_count;
  coordinate_pointer coord_local_buffer =
    coord_buffer + remote_count + request_count;

  // pack search coordinates
  for (size_t i = 0, k = 0; i < count; ++i) {
    if (ranks[i] == comm_rank) {
      coord_local_buffer[k][0] = search_coords[i][0];
      coord_local_buffer[k][1] = search_coords[i][1];
      coord_local_buffer[k][2] = search_coords[i][2];
      ++k;
    } else {
      size_t displ = sdispls[ranks[i]+1]++;
      coord_send_buffer[displ][0] = search_coords[i][0];
      coord_send_buffer[displ][1] = search_coords[i][1];
      coord_send_buffer[displ][2] = search_coords[i][2];
    }
  }

  MPI_Datatype dt_coord;
  yac_mpi_call(MPI_Type_contiguous(3, MPI_DOUBLE, &dt_coord), comm);
  yac_mpi_call(MPI_Type_commit(&dt_coord), comm);

  // redistribute search coordinates
  yac_alltoallv_p2p(
    coord_send_buffer, sendcounts, sdispls,
    coord_recv_buffer, recvcounts, rdispls,
    sizeof(*coord_send_buffer), dt_coord, comm);

  yac_mpi_call(MPI_Type_free(&dt_coord), comm);

  size_t * local_cells =
    xmalloc((request_count + local_count) * sizeof(*local_cells));

  // do local search
  dist_grid_pair_do_point_search_local(
    grid_pair, grid_name, coord_recv_buffer, request_count + local_count,
    local_cells, coord_in_cell);

  struct single_remote_point * single_remote_point_buffer =
    xmalloc((remote_count + request_count) *
            sizeof(*single_remote_point_buffer));
  struct single_remote_point * id_send_buffer = single_remote_point_buffer;
  struct single_remote_point * id_recv_buffer = single_remote_point_buffer +
                                                request_count;

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name);

  // pack global ids of found source cells
  for (size_t i = 0; i < request_count; ++i) {
    size_t cell_idx = local_cells[i];
    id_send_buffer[i].data.rank = comm_rank;
    if (cell_idx != SIZE_MAX) {
      id_send_buffer[i].global_id = dist_grid->ids[CELL][cell_idx];
      id_send_buffer[i].data.orig_pos = (uint64_t)cell_idx;
    } else {
      id_send_buffer[i].global_id = XT_INT_MAX;
      id_send_buffer[i].data.orig_pos = UINT64_MAX;
    }
  }

  MPI_Datatype single_remote_point_dt =
    yac_get_single_remote_point_mpi_datatype(comm);

  // redistribute results (global ids of found source cells)
  yac_alltoallv_p2p(
    id_send_buffer, recvcounts, rdispls, id_recv_buffer, sendcounts, sdispls,
    sizeof(*id_send_buffer), single_remote_point_dt, comm);

  yac_mpi_call(MPI_Type_free(&single_remote_point_dt), comm);

  size_t * new_local_cells =
    xmalloc(remote_count * sizeof(*new_local_cells));

  // convert all remote ids to local ones, extend local dist_grid data,
  // if necessary
  yac_dist_grid_single_remote_point_to_local(
    dist_grid, id_recv_buffer, remote_count, CELL, new_local_cells);

  // extract results from local and remote search
  for (size_t i = 0, k = 0; i < count; ++i) {
    if (ranks[i] == comm_rank) {
      cells[i] = local_cells[request_count + k];
      ++k;
    } else {
      size_t displ = sdispls[ranks[i]]++;
      cells[i] = new_local_cells[displ];
    }
  }

  free(new_local_cells);
  free(single_remote_point_buffer);
  free(local_cells);
  free(coord_buffer);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  free(ranks);
}

void yac_dist_grid_pair_do_point_search(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  coordinate_pointer search_coords, size_t count, size_t * cells) {

  yac_dist_grid_pair_do_point_search_(
    grid_pair, grid_name, search_coords, count, cells, coord_in_cell);
}

void yac_dist_grid_pair_do_point_search_gc(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  coordinate_pointer search_coords, size_t count, size_t * cells) {

  yac_dist_grid_pair_do_point_search_(
    grid_pair, grid_name, search_coords, count, cells, coord_in_cell_gc);
}

static struct point_sphere_part_search * yac_dist_grid_get_field_sphere_part(
  struct dist_grid * dist_grid, struct interp_field field) {

  coordinate_pointer field_coords =
    yac_dist_grid_get_field_coords(dist_grid, field);
  yac_int const * global_ids =
    yac_dist_grid_get_global_ids(dist_grid, field.location);
  int const * mask = yac_dist_grid_get_field_mask(dist_grid, field);
  size_t total_count = yac_dist_grid_get_count(dist_grid, field.location);

  if (mask == NULL)
    return
      yac_point_sphere_part_search_new(
        total_count, field_coords, global_ids);
  else
    return
      yac_point_sphere_part_search_mask_new(
        total_count, field_coords, global_ids, mask);
}

static void yac_dist_grid_get_n_unmasked_local_points(
  struct dist_grid * dist_grid, struct interp_field field, int comm_rank,
  size_t n, struct single_remote_point * points) {

  size_t total_count = yac_dist_grid_get_count(dist_grid, field.location);
  int const * mask = yac_dist_grid_get_field_mask(dist_grid, field);
  int const * core_mask =
    yac_dist_grid_get_core_mask(dist_grid, field.location);
  yac_int const * global_ids =
    yac_dist_grid_get_global_ids(dist_grid, field.location);

  if (mask == NULL) {

    for (size_t i = 0, j = 0; i < total_count; ++i) {
      if (core_mask[i]) {
        points[j].global_id = global_ids[i];
        points[j].data.rank = comm_rank;
        points[j].data.orig_pos = i;
        if (n == ++j) return;
      }
    }

  } else {

    for (size_t i = 0, j = 0; i < total_count; ++i) {
      if (core_mask[i] && mask[i]) {
        points[j].global_id = global_ids[i];
        points[j].data.rank = comm_rank;
        points[j].data.orig_pos = i;
        if (n == ++j) return;
      }
    }
  }
}

static MPI_Datatype yac_get_nnn_search_result_mpi_datatype(
  size_t n, MPI_Comm comm) {

  struct nnn_search_result dummy;
  MPI_Datatype single_remote_point_dt =
                 yac_get_single_remote_point_mpi_datatype(comm),
               single_nnn_search_result_dt_,
               nnn_search_result_dt;
  int array_of_blocklengths[] = {1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.point_info) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.cos_angle) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] = {single_remote_point_dt, MPI_DOUBLE};
  yac_mpi_call(
    MPI_Type_create_struct(
      2, array_of_blocklengths, array_of_displacements,
      array_of_types, &single_nnn_search_result_dt_), comm);
  MPI_Datatype single_nnn_search_result_dt =
    yac_create_resized(single_nnn_search_result_dt_, sizeof(dummy), comm);
  yac_mpi_call(
    MPI_Type_contiguous(
      (int)n, single_nnn_search_result_dt, &nnn_search_result_dt), comm);
  yac_mpi_call(MPI_Type_commit(&nnn_search_result_dt), comm);
  yac_mpi_call(MPI_Type_free(&single_nnn_search_result_dt), comm);
  yac_mpi_call(MPI_Type_free(&single_remote_point_dt), comm);
  return nnn_search_result_dt;
}

static int compare_nnn_search_results_cos_angle(
  void const * a, void const * b) {

  struct nnn_search_result const * result_a =
    (struct nnn_search_result const *)a;
  struct nnn_search_result const * result_b =
    (struct nnn_search_result const *)b;

  int ret = (result_a->cos_angle < result_b->cos_angle) -
            (result_a->cos_angle > result_b->cos_angle);
  if (ret) return ret;
  return (result_a->point_info.global_id > result_b->point_info.global_id) -
         (result_a->point_info.global_id < result_b->point_info.global_id);
}

static void insert_nnn_result_points(
  struct nnn_search_result * array, size_t array_size,
  struct nnn_search_result * insert, size_t insert_size) {

  if (insert_size == 0) return;

  // sort arrays first by cos_angle and second by global id
  qsort(array, array_size, sizeof(*array),
        compare_nnn_search_results_cos_angle);
  qsort(insert, insert_size, sizeof(*insert),
        compare_nnn_search_results_cos_angle);

  struct nnn_search_result temp_results[array_size];
  memcpy(temp_results, array, array_size * sizeof(*array));

  for (size_t i = 0, j = 0, k = 0; k < array_size; ++k) {

    int cmp =
      (j < insert_size)?
        compare_nnn_search_results_cos_angle(
          &temp_results[i], insert + j):(-1);

    if (cmp <= 0) {
      if (k != i) array[k] = temp_results[i];
      ++i;
    } else {
      array[k] = insert[j];
    }

    j += cmp >= 0;
  }
}

void yac_dist_grid_pair_do_nnn_search(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  coordinate_pointer search_coords, size_t count, size_t * local_ids,
  size_t n, struct interp_field field) {

  MPI_Comm comm = grid_pair->comm;
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name);

  uint64_t unmasked_local_count =
    yac_dist_grid_get_unmasked_local_count(dist_grid, field);

  uint64_t * unmasked_local_counts =
    xmalloc((size_t)comm_size * sizeof(*unmasked_local_counts));

  // exchange number of local source points and target points
  yac_mpi_call(
    MPI_Allgather(
      &unmasked_local_count, 1, MPI_UINT64_T,
      unmasked_local_counts, 1, MPI_UINT64_T, comm), comm);

  // check whether there is a rank with too few source points
  int flag = 0;
  for (int i = 0; i < comm_size; ++i)
    flag |= unmasked_local_counts[i] < (uint64_t)n;

  // if ranks with insufficient number of local source points
  if (flag) {

    uint64_t global_num_unmasked_count = 0;
    for (int i = 0; i < comm_size; ++i)
      global_num_unmasked_count += unmasked_local_counts[i];

    YAC_ASSERT(
      (size_t)global_num_unmasked_count >= n,
      "ERROR(yac_dist_grid_pair_do_nnn_search): insufficient number of "
      "unmasked points")

    size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
    yac_get_comm_buffers(
      1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

    // get ranks of processes that have additional data or require data in
    // order to have enough source points to do a initial nnn search
    int * flag_buffer = xcalloc(2 * (size_t)comm_size, sizeof(*flag_buffer));
    int * send_flags = flag_buffer;
    int * recv_flags = flag_buffer + comm_size;
    yac_proc_sphere_part_get_neigh_ranks(
      grid_pair->proc_sphere_part, unmasked_local_counts, (uint64_t)n,
      send_flags, recv_flags, comm_rank, comm_size);
    for (int i = 0; i < comm_size; ++i) {
      sendcounts[i] = (size_t)send_flags[i];
      recvcounts[i] = (size_t)recv_flags[i];
    }
    free(flag_buffer);

    size_t local_send_count = (size_t)(MIN(unmasked_local_count, n));

    size_t raccu = 0;
    for (int i = 0; i < comm_size; ++i) {
      sdispls[i] = 0;
      rdispls[i] = raccu;
      sendcounts[i] *= local_send_count;
      raccu += (recvcounts[i] *= (int)(MIN(unmasked_local_counts[i], n)));
    }

    size_t recv_count = recvcounts[comm_size-1] + rdispls[comm_size-1];

    struct single_remote_point * single_remote_point_buffer =
      xmalloc(
        (local_send_count + recv_count) * sizeof(*single_remote_point_buffer));
    struct single_remote_point * local_send_ids = single_remote_point_buffer;
    struct single_remote_point * recv_ids =
      single_remote_point_buffer + local_send_count;

    // get local source points that can be sent to other processes
    yac_dist_grid_get_n_unmasked_local_points(
      dist_grid, field, comm_rank, local_send_count, local_send_ids);

    MPI_Datatype single_remote_point_dt =
      yac_get_single_remote_point_mpi_datatype(comm);

    // exchange source points (integrate points into local data)
    yac_alltoallv_p2p(
      local_send_ids, sendcounts, sdispls, recv_ids, recvcounts, rdispls,
      sizeof(*local_send_ids), single_remote_point_dt, comm);

    yac_mpi_call(MPI_Type_free(&single_remote_point_dt), comm);

    size_t * dummy = xmalloc(recv_count * sizeof(*dummy));

    // convert all remote ids to local ones, extend local dist_grid data,
    // if necessary
    yac_dist_grid_single_remote_point_to_local(
      dist_grid, recv_ids, recv_count, field.location, dummy);

    free(dummy);
    free(single_remote_point_buffer);
    yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  }
  free(unmasked_local_counts);

  struct point_sphere_part_search * sphere_part =
    yac_dist_grid_get_field_sphere_part(dist_grid, field);

  double * cos_angles = NULL;
  size_t cos_angles_array_size = 0;
  coordinate_pointer result_coords = NULL;
  size_t result_coords_array_size = 0;
  size_t * result_points = NULL;
  size_t result_points_array_size = 0;
  size_t * num_results = xcalloc(count, sizeof(*num_results));

  // do local search
  yac_point_sphere_part_search_NNN(
    sphere_part, count, search_coords, n, &cos_angles, &cos_angles_array_size,
    &result_coords, &result_coords_array_size, &result_points,
    &result_points_array_size, num_results);

  size_t max_num_results = 0;
  for (size_t i = 0; i < count; ++i)
    if (max_num_results < num_results[i]) max_num_results = num_results[i];

  // extract results (only the local core points)
  struct nnn_search_result * results = xmalloc(n * count * sizeof(*results));
  yac_int const * global_ids =
    yac_dist_grid_get_global_ids(dist_grid, field.location);

  struct nnn_search_result * temp_results =
    xmalloc(max_num_results * sizeof(*temp_results));

  for (size_t i = 0, k = 0; i < count; ++i) {

    struct nnn_search_result * curr_results = results + i * n;
    for (size_t j = 0; j < n; ++j, ++k) {
      size_t curr_local_id = result_points[k];
      curr_results[j].point_info.global_id = global_ids[curr_local_id];
      curr_results[j].point_info.data.rank = comm_rank;
      curr_results[j].point_info.data.orig_pos = (uint64_t)curr_local_id;
      curr_results[j].cos_angle = cos_angles[k];
    }

    for (size_t j = n; j < num_results[i]; ++j, ++k) {
      size_t curr_local_id = result_points[k];
      temp_results[j-n].point_info.global_id = global_ids[curr_local_id];
      temp_results[j-n].point_info.data.rank = comm_rank;
      temp_results[j-n].point_info.data.orig_pos = (uint64_t)curr_local_id;
      temp_results[j-n].cos_angle = cos_angles[k];
    }
    insert_nnn_result_points(curr_results, n, temp_results, num_results[i] - n);
  }

  coordinate_pointer field_coords =
    yac_dist_grid_get_field_coords(dist_grid, field);
  int * request_ranks = NULL;
  size_t request_ranks_array_size = 0;
  size_t num_request_ranks = 0;
  int * num_requests = xmalloc(count * sizeof(*num_requests));

  // for all search points
  for (size_t i = 0; i < count; ++i) {

    // generate bounding circles for all search points
    struct bounding_circle result_bnd_circle;
    memcpy(result_bnd_circle.base_vector, search_coords[i],
           3 * sizeof(search_coords[0][0]));
    result_bnd_circle.inc_angle =
      get_vector_angle_2(
        search_coords[i],
        field_coords[results[(i + 1) * n - 1].point_info.data.orig_pos]);

    ENSURE_ARRAY_SIZE(request_ranks, request_ranks_array_size,
                      num_request_ranks + (size_t)comm_size);

    // search for processes that might be able to contribute to the search
    // results
    int * curr_request_ranks = request_ranks + num_request_ranks;
    yac_proc_sphere_part_do_bnd_circle_search(
      grid_pair->proc_sphere_part, result_bnd_circle,
      curr_request_ranks, num_requests + i);

    // remove requests for local process
    int new_num_requests = 0;
    for (int j = 0; j < num_requests[i]; ++j) {
      if (curr_request_ranks[j] == comm_rank) continue;
      if (new_num_requests != j)
        curr_request_ranks[new_num_requests] = curr_request_ranks[j];
      ++new_num_requests;
    }

    num_request_ranks += (size_t)(num_requests[i] = new_num_requests);
  }

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0; i < num_request_ranks; ++i) sendcounts[request_ranks[i]]++;

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t recv_count = rdispls[comm_size-1] + recvcounts[comm_size-1];

  coordinate_pointer coord_buffer =
    xmalloc((num_request_ranks + recv_count) * sizeof(*coord_buffer));
  coordinate_pointer send_request_coords = coord_buffer;
  coordinate_pointer recv_request_coords = coord_buffer + num_request_ranks;

  // pack bounding circles
  for (size_t i = 0, k = 0; i < count; ++i)
    for (size_t j = 0; j < num_requests[i]; ++j, ++k)
      memcpy(send_request_coords[sdispls[request_ranks[k]+1]++],
             search_coords[i], 3 * sizeof(double));

  MPI_Datatype coord_dt = yac_get_coordinate_mpi_datatype(comm);

  // exchange requests to other processes
  yac_alltoallv_p2p(
    send_request_coords, sendcounts, sdispls,
    recv_request_coords, recvcounts, rdispls,
    sizeof(*send_request_coords), coord_dt, comm);

  yac_mpi_call(MPI_Type_free(&coord_dt), comm);

  num_results = xrealloc(num_results, recv_count * sizeof(*num_results));
  memset(num_results, 0, recv_count * sizeof(*num_results));

  // do local search for requests
  yac_point_sphere_part_search_NNN(
    sphere_part, recv_count, recv_request_coords, n, &cos_angles,
    &cos_angles_array_size, &result_coords, &result_coords_array_size,
    &result_points, &result_points_array_size, num_results);
  yac_delete_point_sphere_part_search(sphere_part);
  free(coord_buffer);

  for (size_t i = 0; i < recv_count; ++i)
    if (max_num_results < num_results[i]) max_num_results = num_results[i];
  temp_results =
    xrealloc(temp_results, max_num_results * sizeof(*temp_results));

  // extract results (only the ones, that are local core points)
  struct nnn_search_result * request_results =
    xmalloc(n * (num_request_ranks + recv_count) * sizeof(*request_results));
  struct nnn_search_result * send_request_results = request_results;
  struct nnn_search_result * recv_request_results =
    request_results + n * recv_count;

  for (size_t i = 0, k = 0; i < recv_count; ++i) {

    struct nnn_search_result * curr_results = send_request_results + i * n;
    for (size_t j = 0; j < n; ++j, ++k) {
      size_t curr_local_id = result_points[k];
      curr_results[j].point_info.global_id = global_ids[curr_local_id];
      curr_results[j].point_info.data.rank = comm_rank;
      curr_results[j].point_info.data.orig_pos = (uint64_t)curr_local_id;
      curr_results[j].cos_angle = cos_angles[k];
    }

    for (size_t j = n; j < num_results[i]; ++j, ++k) {
      size_t curr_local_id = result_points[k];
      temp_results[j-n].point_info.global_id = global_ids[curr_local_id];
      temp_results[j-n].point_info.data.rank = comm_rank;
      temp_results[j-n].point_info.data.orig_pos = (uint64_t)curr_local_id;
      temp_results[j-n].cos_angle = cos_angles[k];
    }
    insert_nnn_result_points(curr_results, n, temp_results, num_results[i]-n);
  }

  free(cos_angles);
  free(result_coords);
  free(result_points);
  free(num_results);

  MPI_Datatype nnn_search_result_dt =
    yac_get_nnn_search_result_mpi_datatype(n, comm);

  // return results
  yac_alltoallv_p2p(
    send_request_results, recvcounts, rdispls,
    recv_request_results, sendcounts, sdispls,
    n * sizeof(*send_request_results), nnn_search_result_dt, comm);

  yac_mpi_call(MPI_Type_free(&nnn_search_result_dt), comm);

  // unpack results of remote search and integrate them into local results
  for (size_t i = 0, k = 0; i < count; ++i) {
    struct nnn_search_result * curr_results = results + n * i;
    for (size_t j = 0; j < num_requests[i]; ++k, ++j) {
      size_t pos = (size_t)(sdispls[request_ranks[k]]++);
      struct nnn_search_result * curr_request_results =
        recv_request_results + n * pos;
      insert_nnn_result_points(curr_results, n, curr_request_results, n);
    }
  }
  free(request_ranks);
  free(request_results);
  free(num_requests);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  size_t remote_count = 0;
  for (size_t i = 0; i < n * count; ++i)
    if (results[i].point_info.data.rank != comm_rank) ++remote_count;

  // gather remote points from search results
  struct single_remote_point * remote_points =
    xmalloc(remote_count * sizeof(*remote_points));
  for (size_t i = 0, k = 0; i < n * count; ++i)
    if (results[i].point_info.data.rank != comm_rank)
      remote_points[k++] = results[i].point_info;

  // convert all remote ids to local ones, extend local dist_grid data,
  // if necessary
  yac_dist_grid_single_remote_point_to_local(
    dist_grid, remote_points, remote_count, field.location,
    local_ids + n * count - remote_count);

  free(remote_points);

  // gather final result
  for (size_t i = 0, offset = n * count - remote_count; i < n * count; ++i) {
    if (results[i].point_info.data.rank == comm_rank)
      local_ids[i] = (size_t)(results[i].point_info.data.orig_pos);
    else local_ids[i] = local_ids[offset++];
  }
  free(temp_results);
  free(results);
}

void yac_dist_grid_pair_do_bnd_circle_search(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  const_bounding_circle_pointer bnd_circles, size_t count, size_t ** cells,
  size_t * num_results_per_bnd_circle, struct interp_field field) {

  MPI_Comm comm = grid_pair->comm;
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name);

  int * num_ranks = xcalloc(count, sizeof(*num_ranks));
  int * rank_buffer = NULL;
  size_t rank_buffer_size = 0;
  size_t rank_buffer_array_size = 0;

  for (size_t i = 0; i < count; ++i) {

    ENSURE_ARRAY_SIZE(rank_buffer, rank_buffer_array_size,
                      rank_buffer_size + (size_t)comm_size);

    // finds all processes whose core area overlaps with the bounding circle
    // of the current cell
    // beware: even if the core are of a process does not overlap with the
    // bounding circle, it may have core cells that overlap nevertheless
    yac_proc_sphere_part_do_bnd_circle_search(
      grid_pair->proc_sphere_part, bnd_circles[i],
      rank_buffer + rank_buffer_size, num_ranks + i);
    rank_buffer_size += (size_t)(num_ranks[i]);
  }

  size_t * size_t_buffer =
    xmalloc(4 * (size_t)comm_size * sizeof(*size_t_buffer));
  size_t * result_sendcounts = size_t_buffer + 0 * comm_size;
  size_t * result_recvcounts = size_t_buffer + 1 * comm_size;
  size_t * result_sdispls =    size_t_buffer + 2 * comm_size;
  size_t * result_rdispls =    size_t_buffer + 3 * comm_size;

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0, offset = 0; i < count; ++i) {
    int curr_num_ranks = num_ranks[i];
    int * ranks = rank_buffer + offset;
    offset += (size_t)curr_num_ranks;
    for (int j = 0; j < curr_num_ranks; ++j) sendcounts[ranks[j]]++;
  }

  // local overlaps do not need to be send around
  size_t local_count = sendcounts[comm_rank];
  sendcounts[comm_rank] = 0;

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t send_count = sdispls[comm_size] + sendcounts[comm_size-1];
  size_t recv_count = rdispls[comm_size-1] + recvcounts[comm_size-1];

  struct bounding_circle * bnd_circle_buffer =
    xmalloc((send_count + recv_count + local_count) *
            sizeof(*bnd_circle_buffer));
  struct bounding_circle * send_buffer = bnd_circle_buffer;
  struct bounding_circle * recv_buffer = bnd_circle_buffer + send_count;
  struct bounding_circle * local_buffer =
    bnd_circle_buffer + send_count + recv_count;

  // pack bounding circles
  for (size_t i = 0, offset = 0, local_offset = 0; i < count; ++i) {
    int curr_num_ranks = num_ranks[i];
    int * ranks = rank_buffer + offset;
    offset += (size_t)curr_num_ranks;
    for (int j = 0; j < curr_num_ranks; ++j) {
      int rank = ranks[j];
      if (rank == comm_rank)
        local_buffer[local_offset++] = bnd_circles[i];
      else
        send_buffer[sdispls[rank + 1]++] = bnd_circles[i];
    }
  }

  MPI_Datatype bnd_circle_dt = yac_get_bounding_circle_mpi_datatype(comm);
  yac_mpi_call(MPI_Type_commit(&bnd_circle_dt), comm);

  // exchange bounding circles
  yac_alltoallv_p2p(
    send_buffer, sendcounts, sdispls, recv_buffer, recvcounts, rdispls,
    sizeof(*send_buffer), bnd_circle_dt, comm);

  yac_mpi_call(MPI_Type_free(&bnd_circle_dt), comm);

  struct bnd_sphere_part_search * cell_sphere_part =
    dist_grid_pair_get_cell_sphere_part(grid_pair, grid_name);

  size_t * local_cells = NULL;
  size_t * num_local_cells_per_bnd_circle =
    xmalloc((recv_count + local_count) *
            sizeof(*num_local_cells_per_bnd_circle));

  uint64_t * uint64_t_buffer =
    xmalloc((send_count + recv_count + local_count) *
            sizeof(*uint64_t_buffer));
  uint64_t * num_local_cells_per_bnd_circle_uint64_t = uint64_t_buffer;
  uint64_t * num_remote_cells_per_bnd_circle =
    uint64_t_buffer + recv_count + local_count;

  // search for received bounding circles in local data
  yac_bnd_sphere_part_search_do_bnd_circle_search(
    cell_sphere_part, recv_buffer, recv_count + local_count, &local_cells,
    num_local_cells_per_bnd_circle);

  int const * field_mask =
    (field.location == CELL)?
      yac_dist_grid_get_field_mask(dist_grid, field):NULL;

  struct bounding_circle * cell_bnd_circles = dist_grid->cell_bnd_circles;

  // check field mask and actual overlap of bounding circles
  for (size_t i = 0, offset = 0, new_offset = 0;
       i < recv_count + local_count; ++i) {

    struct bounding_circle * curr_bnd_circle = recv_buffer + i;
    size_t curr_num_results = num_local_cells_per_bnd_circle[i];

    // remove cells whose bounding circle do not overlap with the current one
    uint64_t new_num_results = 0;
    for (size_t j = 0; j < curr_num_results; ++j, ++offset) {
      size_t local_cell_id = local_cells[offset];
      if (!yac_extents_overlap(curr_bnd_circle,
                               cell_bnd_circles + local_cell_id)) continue;
      if ((field_mask == NULL) || (field_mask[local_cell_id])) {
        if (offset != new_offset) local_cells[new_offset] = local_cell_id;
        new_num_results++;
        new_offset++;
      }
    }
    num_local_cells_per_bnd_circle_uint64_t[i] = new_num_results;
  }
  free(num_local_cells_per_bnd_circle);
  free(bnd_circle_buffer);

  // exchange number of results per bounding circle
  yac_alltoallv_p2p(
    num_local_cells_per_bnd_circle_uint64_t, recvcounts, rdispls,
    num_remote_cells_per_bnd_circle, sendcounts, sdispls,
    sizeof(*num_local_cells_per_bnd_circle_uint64_t), MPI_UINT64_T, comm);

  size_t saccu = 0, raccu = 0, soffset = 0, roffset = 0;
  for (int i = 0; i < comm_size; ++i) {

    result_sdispls[i] = saccu;
    result_rdispls[i] = raccu;

    size_t sendcount = recvcounts[i];
    size_t recvcount = sendcounts[i];

    result_sendcounts[i] = 0;
    result_recvcounts[i] = 0;
    for (size_t j = 0; j < sendcount; ++j, ++soffset)
      result_sendcounts[i] +=
        (size_t)(num_local_cells_per_bnd_circle_uint64_t[soffset]);
    for (size_t j = 0; j < recvcount; ++j, ++roffset)
      result_recvcounts[i] +=
        (size_t)(num_remote_cells_per_bnd_circle[roffset]);

    saccu += result_sendcounts[i];
    raccu += result_recvcounts[i];
  }

  // count the number of results for bounding circles, which had a match with#
  // their original process
  size_t result_local_count = 0;
  for (size_t i = recv_count; i < recv_count + local_count; ++i)
    result_local_count += (size_t)(num_local_cells_per_bnd_circle_uint64_t[i]);

  size_t result_send_count = (size_t)(result_sdispls[comm_size-1]) +
                             (size_t)(result_sendcounts[comm_size-1]);
  size_t result_recv_count = (size_t)(result_rdispls[comm_size-1]) +
                             (size_t)(result_recvcounts[comm_size-1]);

  struct single_remote_point * single_remote_point_buffer =
    xmalloc((result_recv_count + result_send_count) *
            sizeof(*single_remote_point_buffer));
  struct single_remote_point * id_send_buffer = single_remote_point_buffer;
  struct single_remote_point * id_recv_buffer = single_remote_point_buffer +
                                                result_send_count;

  yac_int * cell_ids = dist_grid->ids[CELL];

  for (size_t i = 0; i < result_send_count; ++i) {
    size_t local_cell_id = local_cells[i];
    id_send_buffer[i].global_id = cell_ids[local_cell_id];
    id_send_buffer[i].data.rank = comm_rank;
    id_send_buffer[i].data.orig_pos = local_cell_id;
  }

  MPI_Datatype single_remote_point_dt =
    yac_get_single_remote_point_mpi_datatype(comm);

  // redistribute results (global ids of found source cells)
  yac_alltoallv_p2p(
    id_send_buffer, result_sendcounts, result_sdispls,
    id_recv_buffer, result_recvcounts, result_rdispls,
    sizeof(*id_send_buffer), single_remote_point_dt, comm);

  yac_mpi_call(MPI_Type_free(&single_remote_point_dt), comm);

  size_t * new_local_cells =
    xmalloc((result_recv_count + result_local_count) *
            sizeof(*new_local_cells));

  memcpy(new_local_cells + result_recv_count,
         local_cells + result_send_count,
         result_local_count * sizeof(*new_local_cells));
  free(local_cells);

  // convert all remote ids to local ones, extend local dist_grid data,
  // if necessary
  yac_dist_grid_single_remote_point_to_local(
    dist_grid, id_recv_buffer, result_recv_count, CELL, new_local_cells);

  free(single_remote_point_buffer);

  size_t * reorder_idx =
    xmalloc((result_recv_count + result_local_count) * sizeof(*reorder_idx));

  memset(
    num_results_per_bnd_circle, 0, count * sizeof(*num_results_per_bnd_circle));

  for (size_t i = 0, offset = 0, reorder = 0, local_search_idx = recv_count,
              local_offset = result_recv_count; i < count; ++i) {
    int curr_num_ranks = num_ranks[i];
    int * ranks = rank_buffer + offset;
    offset += (size_t)curr_num_ranks;
    for (int j = 0; j < curr_num_ranks; ++j) {
      int rank = ranks[j];
      if (rank == comm_rank) {
        uint64_t curr_num_results =
          num_local_cells_per_bnd_circle_uint64_t[local_search_idx++];
        num_results_per_bnd_circle[i] += (size_t)curr_num_results;
        for (uint64_t k = 0; k < curr_num_results; ++k, ++reorder)
          reorder_idx[local_offset++] = reorder;
      } else {
        size_t rank_pos = sdispls[rank]++;
        uint64_t curr_num_results = num_remote_cells_per_bnd_circle[rank_pos];
        num_results_per_bnd_circle[i] += (size_t)curr_num_results;
        for (uint64_t k = 0; k < curr_num_results; ++k, ++reorder)
          reorder_idx[result_rdispls[rank]++] = reorder;
      }
    }
  }
  free(uint64_t_buffer);
  free(num_ranks);
  free(rank_buffer);
  free(size_t_buffer);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  yac_quicksort_index_size_t_size_t(
    reorder_idx, result_recv_count + result_local_count, new_local_cells);
  free(reorder_idx);

  // remove duplicated results
  for (size_t i = 0, offset = 0, new_offset = 0; i < count; ++i) {

    size_t * curr_local_cells = new_local_cells + offset;
    size_t curr_num_results_per_bnd_circle = num_results_per_bnd_circle[i];
    size_t new_num_results_per_bnd_circle = 0;
    size_t prev_cell = SIZE_MAX;
    offset += curr_num_results_per_bnd_circle;

    yac_quicksort_index_size_t_size_t(
      curr_local_cells, curr_num_results_per_bnd_circle, NULL);

    for (size_t j = 0; j < curr_num_results_per_bnd_circle; ++j) {
      size_t curr_cell = curr_local_cells[j];
      if (curr_cell != prev_cell) {
        new_local_cells[new_offset++] = (prev_cell = curr_cell);
        ++new_num_results_per_bnd_circle;
      }
    }
    num_results_per_bnd_circle[i] = new_num_results_per_bnd_circle;
  }

  *cells = new_local_cells;
}

void yac_dist_grid_pair_do_cell_search(
  struct dist_grid_pair * grid_pair,
  char const * search_grid_name, char const * result_grid_name,
  size_t * search_cells, size_t count, size_t ** result_cells,
  size_t * num_results_per_search_cell, struct interp_field result_field) {

  struct bounding_circle * search_bnd_circles =
    xmalloc(count * sizeof(*search_bnd_circles));

  struct dist_grid * search_dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, search_grid_name);
  struct dist_grid * result_dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, result_grid_name);

  const_bounding_circle_pointer search_grid_cell_bnd_circles =
    search_dist_grid->cell_bnd_circles;

  for (size_t i = 0; i < count; ++i)
    search_bnd_circles[i] = search_grid_cell_bnd_circles[search_cells[i]];

  yac_dist_grid_pair_do_bnd_circle_search(
    grid_pair, result_grid_name, search_bnd_circles, count, result_cells,
    num_results_per_search_cell, result_field);

  size_t total_num_result_cells = 0;

  // struct grid_cell search_cell, result_cell;
  // yac_init_grid_cell(&search_cell);
  // yac_init_grid_cell(&result_cell);

  // filter out obvious mismachtes
  // (currently only the bounding circles are checked)
  for (size_t i = 0, offset = 0; i < count; ++i) {

    size_t curr_num_results_per_bnd_circle = num_results_per_search_cell[i];
    size_t new_num_results_per_search_cell = 0;
    size_t * curr_result_cells = *result_cells + offset;
    offset += curr_num_results_per_bnd_circle;

    // yac_const_basic_grid_data_get_grid_cell(
      // (struct const_basic_grid_data *)search_dist_grid,
      // search_cells[i], &search_cell);

    for (size_t j = 0; j < curr_num_results_per_bnd_circle; ++j) {

      size_t curr_result_cell = curr_result_cells[j];

      // yac_const_basic_grid_data_get_grid_cell(
        // (struct const_basic_grid_data *)result_dist_grid,
        // curr_result_cell, &result_cell);

      // if (yac_check_overlap_cells2(
            // search_cell, search_dist_grid->cell_bnd_circles[search_cells[i]],
            // result_cell, result_dist_grid->cell_bnd_circles[curr_result_cell])) {
      if (yac_extents_overlap(
            search_bnd_circles + i,
            result_dist_grid->cell_bnd_circles + curr_result_cell)) {

        (*result_cells)[total_num_result_cells++] = curr_result_cell;
        ++new_num_results_per_search_cell;
      }
    }

    num_results_per_search_cell[i] = new_num_results_per_search_cell;
  }

  // yac_free_grid_cell(&result_cell);
  // yac_free_grid_cell(&search_cell);

  // *result_cells =
    // xrealloc(*result_cells, total_num_result_cells * sizeof(**result_cells));
  free(search_bnd_circles);
}

static struct bounding_circle compute_edge_bnd_circle(
  struct dist_grid * dist_grid, size_t edge_id) {

  struct bounding_circle bnd_circle;

  size_t * edge_to_vertex = &(dist_grid->edge_to_vertex[edge_id][0]);
  double * vertices[2] =
    {&(dist_grid->vertex_coordinates[edge_to_vertex[0]][0]),
     &(dist_grid->vertex_coordinates[edge_to_vertex[1]][0])};

  bnd_circle.base_vector[0] = vertices[0][0] + vertices[1][0];
  bnd_circle.base_vector[1] = vertices[0][1] + vertices[1][1];
  bnd_circle.base_vector[2] = vertices[0][2] + vertices[1][2];
  normalise_vector(bnd_circle.base_vector);
  bnd_circle.inc_angle =
    half_angle(get_vector_angle_2(vertices[0], vertices[1]));
  bnd_circle.sq_crd = DBL_MAX;

  return bnd_circle;
}

static void yac_dist_grid_get_cell_neighbours(
  struct dist_grid * dist_grid,
  struct proc_sphere_part_node * proc_sphere_part,
  size_t * cells, size_t count, size_t * neighbours) {

  // generate edge to cell
  size_t_2_pointer edge_to_cell =
    generate_edge_to_cell(
      dist_grid->cell_to_edge, dist_grid->num_vertices_per_cell,
      NULL, dist_grid->total_count[CELL], dist_grid->total_count[EDGE]);

  // get maximum number of edges per cell
  size_t max_num_edges_per_cell = 0;
  for (size_t i = 0; i < dist_grid->total_count[CELL]; ++i)
    if (max_num_edges_per_cell < dist_grid->num_vertices_per_cell[i])
      max_num_edges_per_cell = dist_grid->num_vertices_per_cell[i];

  size_t_2_pointer edge_vertices =
    xmalloc(max_num_edges_per_cell * sizeof(*edge_vertices));

  size_t neigh_idx = 0;

  struct missing_edge_neighbour * missing_edge_neighbour = NULL;
  size_t missing_edge_neighbour_array_size = 0;
  size_t num_missing_neighbours = 0;

  // for each cell
  for (size_t i = 0; i < count; ++i) {

    size_t curr_cell = cells[i];

    // get all edges
    size_t curr_num_edges = dist_grid->num_vertices_per_cell[curr_cell];
    size_t const * cell_edges =
      dist_grid->cell_to_edge + dist_grid->cell_to_edge_offsets[curr_cell];
    for (size_t j = 0; j < curr_num_edges; ++j) {
      size_t const * curr_edge_to_vertex =
        dist_grid->edge_to_vertex[cell_edges[j]];
      edge_vertices[j][0] = curr_edge_to_vertex[0];
      edge_vertices[j][1] = curr_edge_to_vertex[1];
    }

    ENSURE_ARRAY_SIZE(
      missing_edge_neighbour, missing_edge_neighbour_array_size,
      num_missing_neighbours + curr_num_edges);

    // get the neighbour cells by following the edges and vertices around
    // the cell
    size_t prev_vertex = edge_vertices[0][0];
    for (size_t j = 0, edge_idx = 0; j < curr_num_edges; ++j, ++neigh_idx) {

      // get the neighbour cell associated with the current edge
      size_t curr_edge = cell_edges[edge_idx];
      size_t * curr_edge_cells = edge_to_cell[curr_edge];
      size_t other_cell = curr_edge_cells[curr_edge_cells[0] == curr_cell];
      neighbours[neigh_idx] = other_cell;

      if (other_cell == SIZE_MAX) {
        struct missing_edge_neighbour * curr_miss_neigh =
          missing_edge_neighbour + num_missing_neighbours++;
        curr_miss_neigh->edge.local_id = curr_edge;
        curr_miss_neigh->edge.global_id = dist_grid->ids[EDGE][curr_edge];
        curr_miss_neigh->cell.local_id = curr_cell;
        curr_miss_neigh->cell.global_id = dist_grid->ids[CELL][curr_cell];
        curr_miss_neigh->neigh_idx = neigh_idx;
      }

      // get an edge that shares a vertex with the current edge
      size_t new_edge_idx = SIZE_MAX;
      for (size_t k = 0; k < curr_num_edges; ++k) {
        if (k == edge_idx) continue;
        else if (edge_vertices[k][0] == prev_vertex) {
          new_edge_idx = k;
          prev_vertex = edge_vertices[k][1];
          break;
        } else if (edge_vertices[k][1] == prev_vertex) {
          new_edge_idx = k;
          prev_vertex = edge_vertices[k][0];
          break;
        }
      }
      YAC_ASSERT(
        new_edge_idx < SIZE_MAX,
        "ERROR(yac_basic_grid_data_get_cell_neighbours): "
        "inconsistent cell_to_edge/edge_to_vertex data")
      edge_idx = new_edge_idx;
    }

    // check whether we went once completely around the cell
    YAC_ASSERT(
      prev_vertex == edge_vertices[0][0],
      "ERROR(yac_basic_grid_data_get_cell_neighbours): "
      "inconsistent cell_to_edge/edge_to_vertex data")
  }

  { // get the missing neighbours
    MPI_Comm comm = dist_grid->comm;
    int comm_rank, comm_size;
    yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
    yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

    size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
    yac_get_comm_buffers(
      1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);
    int * int_buffer =
      xmalloc(
        ((size_t)comm_size + num_missing_neighbours) * sizeof(*int_buffer));
    int * temp_ranks = int_buffer;
    int * num_ranks =  int_buffer + comm_size;
    memset(num_ranks, 0, num_missing_neighbours * sizeof(*num_ranks));

    int * rank_buffer = NULL;
    size_t rank_buffer_array_size = 0;
    size_t rank_buffer_size = 0;

    for (size_t i = 0; i < num_missing_neighbours; ++i) {

      // alternatively get the dist owner ranks for all cell vertices
      int curr_num_ranks;
      yac_proc_sphere_part_do_bnd_circle_search(
        proc_sphere_part,
        compute_edge_bnd_circle(dist_grid, missing_edge_neighbour[i].edge.local_id),
        temp_ranks, &curr_num_ranks);

      ENSURE_ARRAY_SIZE(rank_buffer, rank_buffer_array_size,
                        rank_buffer_size + (size_t)curr_num_ranks);

      for (int j = 0; j < curr_num_ranks; ++j) {
        int curr_rank = temp_ranks[j];
        if (curr_rank != comm_rank) {
          sendcounts[curr_rank] += 2;
          num_ranks[i]++;
          rank_buffer[rank_buffer_size++] = curr_rank;
        }
      }
    }

    yac_generate_alltoallv_args(
      1, sendcounts, recvcounts, sdispls, rdispls, comm);

    size_t send_count =
      (sdispls[comm_size] + sendcounts[comm_size-1])/2;
    size_t recv_count =
      (rdispls[comm_size-1] + recvcounts[comm_size-1])/2;

    yac_int * yac_int_buffer =
      xmalloc(2 * (send_count + recv_count) * sizeof(*yac_int_buffer));
    yac_int * send_buffer = yac_int_buffer;
    yac_int * recv_buffer = yac_int_buffer + 2 * send_count;
    size_t * result_reorder_idx =
      xmalloc(send_count * sizeof(*result_reorder_idx));

    // pack send buffer
    for (size_t i = 0, k = 0, rank_offset = 0; i < num_missing_neighbours; ++i) {

      int * curr_rank_buffer = rank_buffer + rank_offset;
      int curr_num_ranks = num_ranks[i];
      rank_offset += (size_t)curr_num_ranks;

      for (int j = 0; j < curr_num_ranks; ++j, ++k) {

        int rank = curr_rank_buffer[j];

        if (rank == comm_rank) continue;

        size_t pos = sdispls[rank + 1];
        sdispls[rank + 1] += 2;

        send_buffer[pos + 0] = missing_edge_neighbour[i].edge.global_id;
        send_buffer[pos + 1] = missing_edge_neighbour[i].cell.global_id;
        result_reorder_idx[pos / 2] = missing_edge_neighbour[i].neigh_idx;
      }
    }
    free(rank_buffer);
    free(missing_edge_neighbour);

    // redistribute requested global ids and remote points of core points
    yac_alltoallv_yac_int_p2p(
      send_buffer, sendcounts, sdispls,
      recv_buffer, recvcounts, rdispls, comm);

    yac_int * request_edge_ids =
      xmalloc(recv_count * sizeof(*request_edge_ids));
    size_t * reorder_idx = xmalloc(recv_count * sizeof(*reorder_idx));
    struct single_remote_point * point_info_buffer =
      xmalloc(
        (send_count + recv_count) * sizeof(*point_info_buffer));
    struct single_remote_point * point_send_buffer = point_info_buffer;
    struct single_remote_point * point_recv_buffer =
      point_info_buffer + recv_count;
    for (size_t i = 0; i < recv_count; ++i) {
      request_edge_ids[i] = recv_buffer[2 * i + 0];
      reorder_idx[i] = i;
    }

    yac_quicksort_index_yac_int_size_t(
      request_edge_ids, recv_count, reorder_idx);
    yac_int * sorted_edge_ids = dist_grid->sorted_ids[EDGE];
    size_t * sorted_edge_reorder_idx = dist_grid->sorted_reorder_idx[EDGE];
    size_t num_edges = dist_grid->total_count[EDGE];
    yac_int * cell_ids = dist_grid->ids[CELL];

    // lookup the requested edge locally
    for (size_t i = 0, j = 0; i < recv_count; ++i) {

      yac_int curr_edge_id = request_edge_ids[i];
      size_t curr_reorder_idx = reorder_idx[i];

      while ((j < num_edges) && (sorted_edge_ids[j] < curr_edge_id)) ++j;

      // if the edge is not available locally
      if ((j >= num_edges) || (sorted_edge_ids[j] != curr_edge_id)) {
        point_send_buffer[curr_reorder_idx] =
        (struct single_remote_point)
          {.global_id = XT_INT_MAX,
           .data = {.rank = comm_rank, .orig_pos = UINT64_MAX}};
        continue;
      }

      // id of the cell that is available on the other process
      yac_int available_edge_cell_id = recv_buffer[2 * curr_reorder_idx + 1];

      size_t * local_edge_cell_ids = edge_to_cell[sorted_edge_reorder_idx[j]];
      yac_int global_edge_cell_ids[2];
      for (int k = 0; k < 2; ++k)
        global_edge_cell_ids[k] =
          (local_edge_cell_ids[k] == SIZE_MAX)?
            XT_INT_MAX:cell_ids[local_edge_cell_ids[k]];

      int missing_idx = global_edge_cell_ids[0] == available_edge_cell_id;

      // consistency check
      YAC_ASSERT(
        (global_edge_cell_ids[missing_idx^1] == available_edge_cell_id) ||
        (global_edge_cell_ids[missing_idx^1] == XT_INT_MAX),
        "ERROR(yac_dist_grid_get_cell_neighbours): "
        "inconsistent cell edge grid data")

      point_send_buffer[curr_reorder_idx].global_id =
        global_edge_cell_ids[missing_idx];
      point_send_buffer[curr_reorder_idx].data.rank = comm_rank;
      point_send_buffer[curr_reorder_idx].data.orig_pos =
        (local_edge_cell_ids[missing_idx] == SIZE_MAX)?
          (uint64_t)UINT64_MAX:local_edge_cell_ids[missing_idx];
    }
    free(reorder_idx);
    free(request_edge_ids);
    free(yac_int_buffer);

    for (int i = 0; i < comm_size; ++i) {
      sdispls[i] /= 2;
      rdispls[i] /= 2;
      sendcounts[i] /= 2;
      recvcounts[i] /= 2;
    }

    MPI_Datatype single_remote_point_dt =
      yac_get_single_remote_point_mpi_datatype(comm);

    yac_alltoallv_p2p(
      point_send_buffer, recvcounts, rdispls,
      point_recv_buffer, sendcounts, sdispls,
      sizeof(*point_send_buffer), single_remote_point_dt, comm);

    yac_mpi_call(MPI_Type_free(&single_remote_point_dt), comm);

    struct single_remote_point_reorder * results =
      xmalloc(send_count * sizeof(*results));

    for (size_t i = 0; i < send_count; ++i) {
      results[i].data = point_recv_buffer[i];
      results[i].reorder_idx = result_reorder_idx[i];
    }

    qsort(results, send_count, sizeof(*results),
          compare_single_remote_point_reorder_reorder_idx);

    // remove duplicated results
    size_t result_count = 0;
    yac_int prev_global_id = XT_INT_MAX;
    for (size_t i = 0, prev_reorder_idx = SIZE_MAX; i < send_count; ++i) {

      // if the current result does not contain useful data
      yac_int curr_global_id = results[i].data.global_id;
      if (curr_global_id == XT_INT_MAX) continue;

      size_t curr_reorder_idx = results[i].reorder_idx;
      if (curr_reorder_idx != prev_reorder_idx){

        results[result_count++] = results[i];
        prev_reorder_idx = curr_reorder_idx;
        prev_global_id = curr_global_id;

      } else {

        YAC_ASSERT(
          prev_global_id == curr_global_id,
          "ERROR(yac_dist_grid_get_cell_neighbours): "
          "inconsistent cell edge data")
      }
    }
    for (size_t i = 0; i < result_count; ++i) {
      point_send_buffer[i] = results[i].data;
      result_reorder_idx[i] = results[i].reorder_idx;
    }
    free(results);

    size_t * local_ids = xmalloc(result_count * sizeof(*local_ids));

    yac_dist_grid_single_remote_point_to_local(
      dist_grid, point_send_buffer, result_count, CELL, local_ids);

    for (size_t i = 0; i < result_count; ++i)
      neighbours[result_reorder_idx[i]] = local_ids[i];

    free(local_ids);
    free(result_reorder_idx);
    free(point_send_buffer);
    free(int_buffer);
    yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  }

  free(edge_vertices);
  free(edge_to_cell);
}

void yac_dist_grid_pair_get_cell_neighbours(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  size_t * cells, size_t count, size_t * neighbours) {

  yac_dist_grid_get_cell_neighbours(
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name),
    grid_pair->proc_sphere_part, cells, count, neighbours);
}

struct remote_point * yac_dist_grid_get_remote_points(
  struct dist_grid * dist_grid, enum yac_location location,
  size_t * points, size_t count) {

  struct remote_point * remote_points = xmalloc(count * sizeof(*remote_points));

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) ||(location == EDGE),
    "ERROR(yac_dist_grid_get_remote_points): invalid location")
  yac_int * global_ids = dist_grid->ids[location];
  struct remote_point_infos * point_infos = dist_grid->owners[location];

  for (size_t i = 0; i < count; ++i) {
    remote_points[i].global_id = global_ids[points[i]];
    remote_points[i].data = point_infos[points[i]];
  }

  return remote_points;
}

static inline int
compute_bucket(yac_int value, int comm_size) {
  return (int)(value / 128) % comm_size;
}

static int get_global_id_pack_size(MPI_Comm comm) {

  int global_id_pack_size;

  yac_mpi_call(
    MPI_Pack_size(1, yac_int_dt, comm, &global_id_pack_size), comm);

  return global_id_pack_size;
}

static void pack_global_id(
  yac_int global_id, void * buffer, int buffer_size, int * position,
  MPI_Comm comm) {

  yac_mpi_call(
      MPI_Pack(&global_id, 1, yac_int_dt, buffer,
               buffer_size, position, comm), comm);
}

static void unpack_global_id(
  void * buffer, int buffer_size, int * position, yac_int * global_id,
  MPI_Comm comm) {

  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, global_id, 1,
               yac_int_dt, comm), comm);
}

static int get_single_remote_point_pack_size(
  MPI_Datatype single_remote_point_dt, MPI_Comm comm) {

  int pack_size;

  yac_mpi_call(
    MPI_Pack_size(1, single_remote_point_dt, comm, &pack_size), comm);

  return pack_size;
}

static void yac_single_remote_point_pack(
  struct single_remote_point * point,
  void * buffer, int buffer_size, int * position,
  MPI_Datatype single_remote_point_dt, MPI_Comm comm) {

  yac_mpi_call(
      MPI_Pack(point, 1, single_remote_point_dt, buffer,
               buffer_size, position, comm), comm);
}

static void yac_single_remote_point_unpack(
  void * buffer, int buffer_size, int * position,
  struct single_remote_point * point, MPI_Datatype single_remote_point_dt,
  MPI_Comm comm) {

  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, point, 1,
               single_remote_point_dt, comm), comm);
}

void yac_dist_grid_global_to_local(
  struct dist_grid * dist_grid, enum yac_location location,
  yac_int * global_ids, size_t count, size_t * local_ids) {

  MPI_Comm comm = dist_grid->comm;
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * size_t_buffer =
    xmalloc((8 * (size_t)comm_size + 1) * sizeof(*size_t_buffer));
  size_t * sendcounts           = size_t_buffer + 0 * comm_size;
  size_t * recvcounts           = size_t_buffer + 2 * comm_size;
  size_t * sdispls              = size_t_buffer + 4 * comm_size;
  size_t * rdispls              = size_t_buffer + 5 * comm_size + 1;
  size_t * total_sendcounts     = size_t_buffer + 6 * comm_size + 1;
  size_t * total_recvcounts     = size_t_buffer + 7 * comm_size + 1;
  memset(sendcounts, 0, 2 * (size_t)comm_size * sizeof(*sendcounts));

  size_t * core_points, core_count;
  struct interp_field dummy_interp_field =
    {.location = location,
     .coordinates_idx = SIZE_MAX,
     .masks_idx = SIZE_MAX};
  yac_dist_grid_get_local_unmasked_points(
    dist_grid, dummy_interp_field, &core_points, &core_count);
  yac_int const * grid_global_ids =
    yac_dist_grid_get_global_ids(dist_grid, location);

  int * rank_buffer = xmalloc((count + core_count) * sizeof(*rank_buffer));
  int * core_point_ranks = rank_buffer;
  int * global_id_ranks = rank_buffer + core_count;

  size_t * global_id_reorder_idx =
    xmalloc(count * sizeof(*global_id_reorder_idx));

  for (size_t i = 0; i < count; ++i) {
    int rank =
      (global_id_ranks[i] = compute_bucket(global_ids[i], comm_size));
    sendcounts[2 * rank + 0]++;
    global_id_reorder_idx[i] = i;
  }

  yac_quicksort_index_int_size_t(global_id_ranks, count, global_id_reorder_idx);

  for (size_t i = 0; i < core_count; ++i) {
    size_t point_idx = core_points[i];
    int rank =
      (core_point_ranks[i] =
         compute_bucket(grid_global_ids[point_idx], comm_size));
    sendcounts[2 * rank + 1]++;
  }

  // sort core points by rank
  yac_quicksort_index_int_size_t(core_point_ranks, core_count, core_points);

  // exchange number of requsted global ids, number of core points in dist_grid
  // and total pack size
  yac_mpi_call(MPI_Alltoall(sendcounts, 2, YAC_MPI_SIZE_T,
                            recvcounts, 2, YAC_MPI_SIZE_T, comm), comm);

  size_t recv_core_count = 0;
  for (int i = 0; i < comm_size; ++i)
    recv_core_count += recvcounts[2*i+1];

  MPI_Datatype single_remote_point_dt =
    yac_get_single_remote_point_mpi_datatype(comm);
  int global_id_pack_size = get_global_id_pack_size(comm);
  int core_point_pack_size =
    get_single_remote_point_pack_size(single_remote_point_dt, comm);

  for (int i = 0; i < comm_size; ++i) {
    total_sendcounts[i] = sendcounts[2*i+0] * (size_t)global_id_pack_size +
                          sendcounts[2*i+1] * (size_t)core_point_pack_size;
    total_recvcounts[i] = recvcounts[2*i+0] * (size_t)global_id_pack_size +
                          recvcounts[2*i+1] * (size_t)core_point_pack_size;
  }

  size_t saccu = 0, raccu = 0;
  sdispls[0] = 0;
  for (int i = 0; i < comm_size; ++i) {
    sdispls[i+1] = saccu;
    rdispls[i] = raccu;
    saccu += total_sendcounts[i];
    raccu += total_recvcounts[i];
  }

  size_t send_size = sdispls[comm_size] + total_sendcounts[comm_size-1];
  size_t recv_size = rdispls[comm_size-1] + total_recvcounts[comm_size-1];
  void * pack_buffer = xmalloc(send_size + recv_size);
  void * send_buffer = pack_buffer;
  void * recv_buffer = (void*)((unsigned char *)pack_buffer + send_size);

  // pack data
  for (size_t i = 0; i < count; ++i) {
    yac_int curr_global_id = global_ids[global_id_reorder_idx[i]];
    int rank = global_id_ranks[i];
    size_t pos = sdispls[rank + 1];
    int position = 0;
    pack_global_id(
      curr_global_id, (void*)((unsigned char*)send_buffer + pos),
      global_id_pack_size, &position, comm);
    sdispls[rank + 1] += global_id_pack_size;
  }
  for (size_t i = 0; i < core_count; ++i) {
    size_t point_idx = core_points[i];
    int rank = core_point_ranks[i];
    size_t pos = sdispls[rank + 1];
    struct single_remote_point curr_core_point;
    curr_core_point.global_id = grid_global_ids[point_idx];
    curr_core_point.data.rank = comm_rank;
    curr_core_point.data.orig_pos = point_idx;
    int position = 0;
    yac_single_remote_point_pack(
      &curr_core_point,
      (void*)((unsigned char*)send_buffer + pos),
      core_point_pack_size, &position, single_remote_point_dt, comm);
    sdispls[rank + 1] += core_point_pack_size;
  }
  free(rank_buffer);
  free(core_points);

  // redistribute requested global ids and remote points of core points
  yac_alltoallv_packed_p2p(
    send_buffer, total_sendcounts, sdispls,
    recv_buffer, total_recvcounts, rdispls, comm);

  size_t num_requested_ids = 0;
  for(int i = 0; i < comm_size; ++i)
    num_requested_ids += recvcounts[2*i+0];

  yac_int * request_global_ids =
    xmalloc(num_requested_ids * sizeof(*request_global_ids));
  size_t * reorder_idx =
    xmalloc(num_requested_ids * sizeof(*reorder_idx));
  struct single_remote_point * point_info_buffer =
    xmalloc((recv_core_count + num_requested_ids + count) *
            sizeof(*point_info_buffer));
  struct single_remote_point * recv_core_points = point_info_buffer;
  struct single_remote_point * send_point_info =
    point_info_buffer + recv_core_count;
  struct single_remote_point * recv_point_info =
    point_info_buffer + recv_core_count + num_requested_ids;

  // unpack data
  num_requested_ids = 0;
  recv_core_count = 0;
  for (int i = 0; i < comm_size; ++i) {

    size_t curr_num_requested_ids = recvcounts[2*i+0];
    size_t curr_num_core_points = recvcounts[2*i+1];

    for (size_t j = 0; j < curr_num_requested_ids; ++j, ++num_requested_ids) {

      int position = 0;
      unpack_global_id(
        recv_buffer, global_id_pack_size, &position,
        request_global_ids + num_requested_ids, comm);
      reorder_idx[num_requested_ids] = num_requested_ids;
      recv_buffer = (void*)((unsigned char*)recv_buffer + global_id_pack_size);
    }
    for (size_t j = 0; j < curr_num_core_points; ++j, ++recv_core_count) {

      int position = 0;
      yac_single_remote_point_unpack(
        recv_buffer, core_point_pack_size, &position,
        recv_core_points + recv_core_count, single_remote_point_dt, comm);
      recv_buffer = (void*)((unsigned char*)recv_buffer + core_point_pack_size);
    }
  }
  free(pack_buffer);

  // sort the requested global ids
  yac_quicksort_index_yac_int_size_t(
    request_global_ids, num_requested_ids, reorder_idx);

  // sort the received core remote_points
  qsort(recv_core_points, recv_core_count, sizeof(*recv_core_points),
        compare_single_remote_point_global_id);

  // match request global ids with remote_points
  for (size_t i = 0, j = 0; i < num_requested_ids; ++i) {

    yac_int curr_global_id = request_global_ids[i];
    while ((j < recv_core_count) &&
           (recv_core_points[j].global_id < curr_global_id)) ++j;

    YAC_ASSERT_F(
      (j < recv_core_count) &&
      (recv_core_points[j].global_id == curr_global_id),
      "ERROR(yac_dist_grid_global_to_local): "
      "no matching core point found for global id %zu",
      (size_t)curr_global_id)

    send_point_info[reorder_idx[i]] = recv_core_points[j];
  }
  free(reorder_idx);
  free(request_global_ids);

  saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {

    sdispls[i] = saccu;
    rdispls[i] = raccu;
    int recvcount = sendcounts[2*i+0];
    int sendcount = recvcounts[2*i+0];
    saccu += (sendcounts[i] = sendcount);
    raccu += (recvcounts[i] = recvcount);
  }

  // redistribute remote_point_infos
  yac_alltoallv_p2p(
    send_point_info, sendcounts, sdispls,
    recv_point_info, recvcounts, rdispls,
    sizeof(*send_point_info), single_remote_point_dt, comm);

  free(size_t_buffer);
  yac_mpi_call(MPI_Type_free(&single_remote_point_dt), comm);

  yac_dist_grid_single_remote_point_to_local(
    dist_grid, recv_point_info, count, location, local_ids);

  yac_quicksort_index_size_t_size_t(global_id_reorder_idx, count, local_ids);

  free(global_id_reorder_idx);
  free(point_info_buffer);
}

/**
 * returns for each selected vertex the list of all cells surrounding
 * this vertex
 * @param[in]  grid_pair            distributed grid pair
 * @param[in]  grid_name            grid name of the grid on which this
 *                                  routine is supposed to work on
 * @param[in]  vertices             local ids of selected vertices
 * @param[in]  count                number of entries in vertices
 * @param[out] cells                list of result cells for the selected
 *                                  vertices
 * @param[out] num_cells_per_vertex number of cells per selected vertex
 * @param[in]  field                if the provided field contains a mask for
 *                                  cells, it will be used
 * @remark the result cells for a vertex are sorted in clock- or counter
 *         clockwise order
 * @remark in case a vertex is not fully surrounded by unmasked cells,
 *         the respective entry in num_cells_per_vertex is 0
 * @remark the user is responsible to free the memory returned through cells
 * @remark the user has to ensure that the array associated to
 *         num_vertices_per_cell is big enough to hold enough elements
 */
static void yac_dist_grid_pair_get_aux_grid_cells(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  size_t * vertices, size_t count, size_t ** cells,
  int * num_cells_per_vertex, struct interp_field field) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name);

  // generate small bounding circle for all vertices
  struct bounding_circle * vertex_bnd_circles =
    xmalloc(count * sizeof(*vertex_bnd_circles));

  struct sin_cos_angle sin_cos_high_tol =
    {yac_angle_tol*100.0, cos(yac_angle_tol*100.0)};

  for (size_t i = 0; i < count; ++i) {
    memcpy(vertex_bnd_circles[i].base_vector,
           dist_grid->vertex_coordinates[vertices[i]], 3 * sizeof(double));
    vertex_bnd_circles[i].inc_angle = sin_cos_high_tol;
    vertex_bnd_circles[i].sq_crd = DBL_MAX;
  }

  size_t * bnd_result_cells;
  size_t * num_result_per_vertex =
    xmalloc(count * sizeof(*num_result_per_vertex));

  // do bounding circle search
  yac_dist_grid_pair_do_bnd_circle_search(
    grid_pair, grid_name, vertex_bnd_circles, count, &bnd_result_cells,
    num_result_per_vertex, field);
  free(vertex_bnd_circles);

  size_t max_num_cell_per_vertex = 0;
  for (size_t i = 0; i < count; ++i)
    if (num_result_per_vertex[i] > max_num_cell_per_vertex)
      max_num_cell_per_vertex = num_result_per_vertex[i];
  size_t (*neigh_vertices)[2] =
    xmalloc(max_num_cell_per_vertex * sizeof(*neigh_vertices));
  size_t * temp_vertex_cell =
    xmalloc(max_num_cell_per_vertex * sizeof(*temp_vertex_cell));
  yac_int * global_cell_ids =
    xmalloc(max_num_cell_per_vertex * sizeof(*global_cell_ids));

  // for all vertices
  for (size_t i = 0, offset = 0, new_offset = 0; i < count; ++i) {

    size_t curr_vertex = vertices[i];
    size_t * curr_cells = bnd_result_cells + offset;
    size_t curr_num_cells_per_vertex = num_result_per_vertex[i];
    size_t new_num_cell_per_vertex = 0;

    // remove all cells that do not contain the current vertex
    for (size_t j = 0; j < curr_num_cells_per_vertex; ++j) {

      size_t curr_cell = curr_cells[j];
      size_t * curr_cell_vertices =
        dist_grid->cell_to_vertex + dist_grid->cell_to_vertex_offsets[curr_cell];
      size_t curr_cell_size = dist_grid->num_vertices_per_cell[curr_cell];

      size_t vertex_idx = 0;
      for (; vertex_idx < curr_cell_size; ++vertex_idx)
        if (curr_cell_vertices[vertex_idx] == curr_vertex) break;

      // if the current cell is not linked to the current vertex
      if (vertex_idx == curr_cell_size) continue;

      // get the vertex adjacent to the current vertex
      neigh_vertices[new_num_cell_per_vertex][0] =
        curr_cell_vertices[((vertex_idx + curr_cell_size) - 1)%curr_cell_size];
      neigh_vertices[new_num_cell_per_vertex][1] =
        curr_cell_vertices[(vertex_idx + 1)%curr_cell_size];

      if (new_offset + new_num_cell_per_vertex != offset + j)
        bnd_result_cells[new_offset + new_num_cell_per_vertex] = curr_cell;
      ++new_num_cell_per_vertex;
    }

    offset += curr_num_cells_per_vertex;
    curr_cells = bnd_result_cells + new_offset;

    if (new_num_cell_per_vertex == 0) {
      num_cells_per_vertex[i] = 0;
      continue;
    }

    // determin order of cell around the current vertex
    temp_vertex_cell[0] = curr_cells[0];
    size_t start_neigh_vertex = neigh_vertices[0][0];
    size_t prev_vertex = neigh_vertices[0][1];
    for (size_t j = 1, prev_cell_idx = 0; j < new_num_cell_per_vertex; ++j) {

      size_t k;
      for (k = 0; k < new_num_cell_per_vertex; ++k) {

        if (k == prev_cell_idx) continue;

        int flag = neigh_vertices[k][0] == prev_vertex;
        if (flag || (neigh_vertices[k][1] == prev_vertex)) {
          temp_vertex_cell[j] = curr_cells[k];
          prev_cell_idx = k;
          prev_vertex = neigh_vertices[k][flag];
          break;
        }
      }

      // if we could not find the next neighbour
      // (vertex is at an edge of the grid)
      if (k == new_num_cell_per_vertex) {
        new_num_cell_per_vertex = 0;
        break;
      }
    }
    if ((prev_vertex != start_neigh_vertex) || (new_num_cell_per_vertex < 3))
      new_num_cell_per_vertex = 0;

    new_offset += new_num_cell_per_vertex;
    num_cells_per_vertex[i] = (int)new_num_cell_per_vertex;

    if (new_num_cell_per_vertex == 0) continue;

    // get global ids of all cells
    yac_int min_global_cell_id = XT_INT_MAX;
    size_t min_global_cell_id_idx = SIZE_MAX;
    for (size_t j = 0; j < new_num_cell_per_vertex; ++j) {
      yac_int curr_global_cell_id =
        ((global_cell_ids[j] = dist_grid->ids[CELL][temp_vertex_cell[j]]));
      if (curr_global_cell_id < min_global_cell_id) {
        min_global_cell_id = curr_global_cell_id;
        min_global_cell_id_idx = j;
      }
    }

    // determine order in which to store the cells
    int order =
      (global_cell_ids[
        ((min_global_cell_id_idx + new_num_cell_per_vertex) - 1)%
          new_num_cell_per_vertex] >
       global_cell_ids[(min_global_cell_id_idx + 1)%new_num_cell_per_vertex])?-1:1;

    // store cells in a partition-independent order
    for (size_t j = 0; j < new_num_cell_per_vertex; ++j)
      curr_cells[j] =
        temp_vertex_cell[
          ((int)(min_global_cell_id_idx + new_num_cell_per_vertex) +
           (int)j * order)%(int)new_num_cell_per_vertex];
  }

  *cells = bnd_result_cells;
  free(num_result_per_vertex);
  free(global_cell_ids);
  free(temp_vertex_cell);
  free(neigh_vertices);
}

static inline int compare_size_t(const void * a, const void * b) {

  size_t const * a_ = a, * b_ = b;

  return (*a_ > *b_) - (*b_ > *a_);
}

void yac_dist_grid_pair_get_vertex_neighbours(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  size_t * vertices, size_t count, size_t ** neigh_vertices_,
  int * num_neighs_per_vertex, struct interp_field field) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name);

  // generate small bounding circle for all vertices
  struct bounding_circle * vertex_bnd_circles =
    xmalloc(count * sizeof(*vertex_bnd_circles));

  struct sin_cos_angle sin_cos_high_tol =
    {yac_angle_tol*100.0, cos(yac_angle_tol*100.0)};

  for (size_t i = 0; i < count; ++i) {
    memcpy(vertex_bnd_circles[i].base_vector,
           dist_grid->vertex_coordinates[vertices[i]], 3 * sizeof(double));
    vertex_bnd_circles[i].inc_angle = sin_cos_high_tol;
    vertex_bnd_circles[i].sq_crd = DBL_MAX;
  }

  size_t * bnd_result_cells;
  size_t * num_result_per_vertex =
    xmalloc(count * sizeof(*num_result_per_vertex));

  // do bounding circle search
  struct interp_field dummy_interp_field =
    {.location = CELL,
     .coordinates_idx = SIZE_MAX,
     .masks_idx = SIZE_MAX};
  yac_dist_grid_pair_do_bnd_circle_search(
    grid_pair, grid_name, vertex_bnd_circles, count,
    &bnd_result_cells, num_result_per_vertex, dummy_interp_field);
  free(vertex_bnd_circles);

  size_t total_num_neigh = 0;
  size_t max_num_neigh = 0;
  for (size_t i = 0; i < count; ++i) {
    total_num_neigh += num_result_per_vertex[i];
    if (num_result_per_vertex[i] > max_num_neigh)
      max_num_neigh = num_result_per_vertex[i];
  }

  int const * vertex_mask = NULL;
  if (field.location == CORNER)
    vertex_mask = yac_dist_grid_get_field_mask(dist_grid, field);

  size_t * neigh_vertices =
    xmalloc(total_num_neigh * sizeof(*neigh_vertices));
  size_t * temp_neigh_vertices =
    xmalloc(2 * max_num_neigh * sizeof(*temp_neigh_vertices));
  total_num_neigh = 0;

  // for all vertices
  for (size_t i = 0, offset = 0; i < count; ++i) {

    size_t curr_vertex = vertices[i];
    size_t * curr_cells = bnd_result_cells + offset;
    size_t curr_num_cells_per_vertex = num_result_per_vertex[i];

    size_t curr_num_neigh_vertices = 0;

    // remove all cells that do not contain the current vertex
    for (size_t j = 0; j < curr_num_cells_per_vertex; ++j) {

      size_t curr_cell = curr_cells[j];
      size_t * curr_cell_vertices =
        dist_grid->cell_to_vertex +
        dist_grid->cell_to_vertex_offsets[curr_cell];
      size_t curr_cell_size = dist_grid->num_vertices_per_cell[curr_cell];

      size_t vertex_idx = 0;
      for (; vertex_idx < curr_cell_size; ++vertex_idx)
        if (curr_cell_vertices[vertex_idx] == curr_vertex) break;

      // if the current cell is not linked to the current vertex
      if (vertex_idx == curr_cell_size) continue;


      // get the vertex adjacent to the current vertex
      size_t neigh_vertex_idx =
        curr_cell_vertices[((vertex_idx + curr_cell_size) - 1)%curr_cell_size];
      if ((vertex_mask != NULL) && (vertex_mask[neigh_vertex_idx]))
        temp_neigh_vertices[curr_num_neigh_vertices++] = neigh_vertex_idx;
      neigh_vertex_idx =
        curr_cell_vertices[(vertex_idx + 1)%curr_cell_size];
      if ((vertex_mask != NULL) && (vertex_mask[neigh_vertex_idx]))
        temp_neigh_vertices[curr_num_neigh_vertices++] = neigh_vertex_idx;
    }

    qsort(temp_neigh_vertices, curr_num_neigh_vertices,
          sizeof(*temp_neigh_vertices), compare_size_t);
    yac_remove_duplicates_size_t(
      temp_neigh_vertices, &curr_num_neigh_vertices);
    memcpy(neigh_vertices + total_num_neigh, temp_neigh_vertices,
           curr_num_neigh_vertices * sizeof(*neigh_vertices));

    total_num_neigh += curr_num_neigh_vertices;
    num_neighs_per_vertex[i] = curr_num_neigh_vertices;

    offset += curr_num_cells_per_vertex;
  }
  free(temp_neigh_vertices);

  free(bnd_result_cells);
  free(num_result_per_vertex);

  *neigh_vertices_ = neigh_vertices;
}

void yac_dist_grid_pair_get_aux_grid(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  size_t * cells, size_t count,
  size_t ** vertex_to_cell, size_t ** vertex_to_cell_offsets_,
  int ** num_cells_per_vertex_, struct interp_field field) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name);

  // determine required vertices
  size_t * temp_cells = xmalloc(count * sizeof(*temp_cells));
  int * required_vertices =
    xcalloc(dist_grid->total_count[CORNER], sizeof(*required_vertices));
  memcpy(temp_cells, cells, count * sizeof(*cells));
  qsort(temp_cells, count, sizeof(*temp_cells), compare_size_t);
  for (size_t i = 0, prev_cell = SIZE_MAX; i < count; ++i) {
    size_t curr_cell = temp_cells[i];
    if (curr_cell == SIZE_MAX) break;
    if (curr_cell != prev_cell) {
      prev_cell = curr_cell;
      size_t curr_num_vertices = dist_grid->num_vertices_per_cell[curr_cell];
      size_t const * curr_vertices =
        dist_grid->cell_to_vertex + dist_grid->cell_to_vertex_offsets[curr_cell];
      for (size_t j = 0; j < curr_num_vertices; ++j)
        required_vertices[curr_vertices[j]] = 1;
    }
  }
  free(temp_cells);

  // generate list of all required vertices
  size_t num_unique_vertices = 0;
  for (size_t i = 0; i < dist_grid->total_count[CORNER]; ++i)
    if (required_vertices[i]) ++num_unique_vertices;
  size_t * unique_vertices =
    xmalloc(num_unique_vertices * sizeof(*unique_vertices));
  for (size_t i = 0, j = 0; i < dist_grid->total_count[CORNER]; ++i)
    if (required_vertices[i]) unique_vertices[j++] = i;
  free(required_vertices);

  // get aux cells for all unique vertices
  int * num_cells_per_vertex =
    xcalloc(num_unique_vertices, sizeof(*num_cells_per_vertex));
  yac_dist_grid_pair_get_aux_grid_cells(
    grid_pair, grid_name, unique_vertices, num_unique_vertices,
    vertex_to_cell, num_cells_per_vertex, field);

  int * grid_num_cells_per_vertex =
    xcalloc(
      dist_grid->total_count[CORNER], sizeof(*grid_num_cells_per_vertex));

  for (size_t i = 0; i < num_unique_vertices; ++i)
    grid_num_cells_per_vertex[unique_vertices[i]] =
      num_cells_per_vertex[i];
  free(num_cells_per_vertex);
  free(unique_vertices);

  size_t * vertex_to_cell_offsets =
    xmalloc(dist_grid->total_count[CORNER] * sizeof(*vertex_to_cell_offsets));
  for (size_t i = 0, offset = 0; i < dist_grid->total_count[CORNER]; ++i) {
    vertex_to_cell_offsets[i] = offset;
    offset += grid_num_cells_per_vertex[i];
  }

  *vertex_to_cell_offsets_ = vertex_to_cell_offsets;
  *num_cells_per_vertex_ = grid_num_cells_per_vertex;
}

static void relocate_points(
  size_t ** points, size_t * reorder_idx, int * ranks, size_t count,
  enum yac_location location, size_t local_count, size_t recv_count,
  struct single_remote_point * id_send_buffer,
  struct single_remote_point * id_recv_buffer,
  size_t * sendcounts, size_t * sdispls, size_t * recvcounts, size_t * rdispls,
  MPI_Datatype single_remote_point_dt, MPI_Comm comm,
  struct dist_grid * dist_grid) {

  yac_int const * global_ids =
    yac_dist_grid_get_global_ids(dist_grid, location);

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  size_t * old_points = *points;
  size_t * new_points =
    xmalloc((local_count + recv_count) * sizeof(*new_points));

  for (size_t i = 0, j = 0, k = 0; i < count; ++i) {
    size_t idx = reorder_idx[i];
    size_t curr_point = old_points[idx];
    int curr_rank = ranks[i];
    if (curr_rank == comm_rank) {
      new_points[j++] = curr_point;
    } else {
      id_send_buffer[k].global_id = global_ids[curr_point];
      id_send_buffer[k].data.rank = comm_rank;
      id_send_buffer[k].data.orig_pos = curr_point;
      ++k;
    }
  }

  // redistribute points
  yac_alltoallv_p2p(
    id_send_buffer, sendcounts, sdispls,
    id_recv_buffer, recvcounts, rdispls,
    sizeof(*id_send_buffer), single_remote_point_dt, comm);

  // convert all remote ids to local ones, extend local dist_grid data,
  // if necessary
  yac_dist_grid_single_remote_point_to_local(
    dist_grid, id_recv_buffer, recv_count, location,
    new_points + local_count);

  *points = new_points;
  free(old_points);
}

static void relocate_weights(
  double ** weights, size_t * reorder_idx, int * ranks, size_t count,
  size_t send_count, size_t local_count, size_t recv_count,
  size_t * sendcounts, size_t * sdispls, size_t * recvcounts, size_t * rdispls,
  MPI_Comm comm) {

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  double * old_weights = *weights;
  double * send_buffer = xmalloc(send_count * sizeof(*send_buffer));
  double * recv_buffer =
    xmalloc((local_count + recv_count) * sizeof(*recv_buffer));

  for (size_t i = 0, j = 0, k = 0; i < count; ++i) {
    size_t idx = reorder_idx[i];
    double curr_weight = old_weights[idx];
    int curr_rank = ranks[i];
    if (curr_rank == comm_rank) recv_buffer[j++] = curr_weight;
    else send_buffer[k++] = curr_weight;
  }

  // redistribute points
  yac_alltoallv_p2p(
    send_buffer, sendcounts, sdispls,
    recv_buffer + local_count, recvcounts, rdispls,
    sizeof(*send_buffer), MPI_DOUBLE, comm);

  free(old_weights);
  *weights = recv_buffer;
  free(send_buffer);
}

void yac_dist_grid_determine_dist_vertex_owner(
  struct dist_grid * dist_grid,
  struct proc_sphere_part_node * proc_sphere_part,
  size_t * vertices, size_t count, int * ranks) {

  size_t * reorder_idx = xmalloc(count * sizeof(*reorder_idx));
  for (size_t i = 0; i < count; ++i) reorder_idx[i] = i;
  yac_quicksort_index_size_t_size_t(vertices, count, reorder_idx);

  size_t valid_count;
  for (valid_count = 0;
       (valid_count < count) && (vertices[valid_count] != SIZE_MAX);
       ++valid_count);

  for (size_t i = valid_count; i < count; ++i) ranks[reorder_idx[i]] = 0;

  size_t unique_count = 0;
  size_t prev_vertex = SIZE_MAX;
  for (size_t i = 0; i < valid_count; ++i) {
    size_t curr_vertex = vertices[i];
    if (curr_vertex != prev_vertex) {
      prev_vertex = curr_vertex;
      ++unique_count;
    }
  }

  coordinate_pointer grid_coords = dist_grid->vertex_coordinates;
  coordinate_pointer search_coords =
    xmalloc(unique_count * sizeof(*search_coords));

  prev_vertex = SIZE_MAX;
  for (size_t i = 0, j = 0; i < valid_count; ++i) {
    size_t curr_vertex = vertices[i];
    if (curr_vertex != prev_vertex) {
      prev_vertex = curr_vertex;
      memcpy(
        search_coords[j++], grid_coords[curr_vertex], 3 * sizeof(double));
    }
  }

  int * temp_ranks = xmalloc(unique_count * sizeof(*temp_ranks));
  yac_proc_sphere_part_do_point_search(
    proc_sphere_part, search_coords, unique_count, temp_ranks);

  prev_vertex = SIZE_MAX;
  for (size_t i = 0, j = 0; i < valid_count; ++i) {
    size_t curr_vertex = vertices[i];
    if (curr_vertex != prev_vertex) {
      prev_vertex = curr_vertex;
      ++j;
    }
    ranks[reorder_idx[i]] = temp_ranks[j-1];
  }

  yac_quicksort_index_size_t_size_t(reorder_idx, count, vertices);
  free(temp_ranks);
  free(search_coords);
  free(reorder_idx);
}

void yac_dist_grid_determine_dist_ce_owner(
  struct dist_grid * dist_grid,
  struct proc_sphere_part_node * proc_sphere_part,
  size_t * indices, size_t count, int * ranks,
  size_t (*get_ce_reference_vertex)(struct dist_grid *, size_t)) {

  size_t * reorder_idx = xmalloc(count * sizeof(*reorder_idx));
  for (size_t i = 0; i < count; ++i) reorder_idx[i] = i;
  yac_quicksort_index_size_t_size_t(indices, count, reorder_idx);

  size_t unique_count = 0;
  size_t prev_index = SIZE_MAX;
  for (size_t i = 0; i < count; ++i) {
    size_t curr_index = indices[i];
    if (curr_index != prev_index) {
      prev_index = curr_index;
      ++unique_count;
    }
  }

  size_t * ref_vertices = xmalloc(unique_count * sizeof(*ref_vertices));
  prev_index = SIZE_MAX;
  for (size_t i = 0, j = 0; i < count; ++i) {
    size_t curr_index = indices[i];
    if (curr_index != prev_index) {
      prev_index = curr_index;
      ref_vertices[j++] =
        get_ce_reference_vertex(dist_grid, curr_index);
    }
  }

  int * temp_ranks = xmalloc(unique_count * sizeof(*temp_ranks));
  yac_dist_grid_determine_dist_vertex_owner(
    dist_grid, proc_sphere_part, ref_vertices, unique_count, temp_ranks);
  free(ref_vertices);

  prev_index = SIZE_MAX;
  for (size_t i = 0, j = 0; i < count; ++i) {
    size_t curr_index = indices[i];
    if (curr_index != prev_index) {
      prev_index = curr_index;
      ++j;
    }
    ranks[reorder_idx[i]] = temp_ranks[j-1];
  }

  yac_quicksort_index_size_t_size_t(reorder_idx, count, indices);
  free(temp_ranks);
  free(reorder_idx);
}

void yac_dist_grid_determine_dist_cell_owner(
  struct dist_grid * dist_grid,
  struct proc_sphere_part_node * proc_sphere_part,
  size_t * cells, size_t count, int * ranks) {

  yac_dist_grid_determine_dist_ce_owner(
    dist_grid, proc_sphere_part, cells, count, ranks,
    get_cell_reference_vertex);
}

void yac_dist_grid_determine_dist_edge_owner(
  struct dist_grid * dist_grid, struct proc_sphere_part_node * proc_sphere_part,
  size_t * edges, size_t count, int * ranks) {

  yac_dist_grid_determine_dist_ce_owner(
    dist_grid, proc_sphere_part, edges, count, ranks,
    get_edge_reference_vertex);
}

void yac_dist_grid_pair_determine_dist_owner(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  size_t * points, size_t count, enum yac_location location, int * ranks) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name);

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_dist_grid_pair_determine_dist_owner): invalid location")

  void (*determine_dist_owner[3])(
    struct dist_grid * dist_grid,
    struct proc_sphere_part_node * proc_sphere_part,
    size_t * cells, size_t count, int * ranks) =
    {yac_dist_grid_determine_dist_cell_owner,
     yac_dist_grid_determine_dist_vertex_owner,
     yac_dist_grid_determine_dist_edge_owner};
  determine_dist_owner[location](
    dist_grid, grid_pair->proc_sphere_part, points, count, ranks);
}

void yac_dist_grid_pair_determine_orig_owner(
  struct dist_grid_pair * grid_pair, char const * grid_name,
  size_t * points, size_t count, enum yac_location location, int * ranks) {

  struct dist_grid * dist_grid =
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name);

  struct remote_point_infos * orig_owners;

  YAC_ASSERT(
    (location == CELL) || (location == CORNER) || (location == EDGE),
    "ERROR(yac_dist_grid_pair_determine_orig_owner): invalid location");
  orig_owners = dist_grid->owners[location];

  for (size_t i = 0; i < count; ++i) {

    struct remote_point_infos * curr_orig_owner = orig_owners + points[i];

    int min_rank;
    int curr_count = curr_orig_owner->count;
    if (curr_count == 1) {
      min_rank = curr_orig_owner->data.single.rank;
    } else {
      min_rank = curr_orig_owner->data.multi[0].rank;
      for (int j = 1; j < curr_count; ++j) {
        int curr_rank = curr_orig_owner->data.multi[j].rank;
        if (min_rank > curr_rank) min_rank = curr_rank;
      }
    }
    ranks[i] = min_rank;
  }
}

void yac_dist_grid_pair_relocate_point_pairs(
  struct dist_grid_pair * grid_pair, int a_is_ref, int to_dist_owner,
  char const * grid_name_a, size_t ** points_a, enum yac_location location_a,
  char const * grid_name_b, size_t ** points_b, enum yac_location location_b,
  double ** weights, size_t * count) {

  size_t count_ = *count;

  MPI_Comm comm = grid_pair->comm;
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // check whether we have to exchange weights
  int weight_flag_local =
    (count_ > 0) && (weights != NULL) && (*weights != NULL);
  int weight_flag;
  yac_mpi_call(MPI_Allreduce(&weight_flag_local, &weight_flag, 1,
                             MPI_INT, MPI_MAX, comm), comm);

  // if there are points defined locally for the current grid
  YAC_ASSERT(
    (count_ <= 0) || weight_flag_local || !weight_flag,
    "ERROR(yac_dist_grid_pair_relocate_point_pairs): weights")

  // get the owner ranks for the reference points
  int * ranks = xmalloc(count_ * sizeof(ranks));
  size_t * reorder_idx = xmalloc(count_ * sizeof(reorder_idx));
  {
    char const * grid_name = (a_is_ref)?grid_name_a:grid_name_b;
    size_t * points = (a_is_ref)?*points_a:*points_b;
    enum yac_location location = (a_is_ref)?location_a:location_b;
    if (to_dist_owner)
      yac_dist_grid_pair_determine_dist_owner(
        grid_pair, grid_name, points, count_, location, ranks);
    else
      yac_dist_grid_pair_determine_orig_owner(
        grid_pair, grid_name, points, count_, location, ranks);
  }
  for (size_t i = 0; i < count_; ++i) reorder_idx[i] = i;
  yac_quicksort_index_int_size_t(ranks, count_, reorder_idx);

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);
  for (size_t i = 0; i < count_; ++i) sendcounts[ranks[i]]++;

  size_t local_count = sendcounts[comm_rank];
  sendcounts[comm_rank] = 0;

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t send_count = sdispls[comm_size] + sendcounts[comm_size-1];
  size_t recv_count = rdispls[comm_size-1] + recvcounts[comm_size-1];

  struct single_remote_point * single_remote_point_buffer =
    xmalloc((send_count + recv_count) *
            sizeof(*single_remote_point_buffer));
  struct single_remote_point * id_send_buffer = single_remote_point_buffer;
  struct single_remote_point * id_recv_buffer = single_remote_point_buffer +
                                                send_count;

  MPI_Datatype single_remote_point_dt =
    yac_get_single_remote_point_mpi_datatype(comm);

  relocate_points(
    points_a, reorder_idx, ranks, count_, location_a,
    local_count, recv_count, id_send_buffer, id_recv_buffer,
    sendcounts, sdispls + 1, recvcounts, rdispls,
    single_remote_point_dt, comm,
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name_a));
  relocate_points(
    points_b, reorder_idx, ranks, count_, location_b,
    local_count, recv_count, id_send_buffer, id_recv_buffer,
    sendcounts, sdispls + 1, recvcounts, rdispls,
    single_remote_point_dt, comm,
    yac_dist_grid_pair_get_dist_grid(grid_pair, grid_name_b));
  yac_mpi_call(MPI_Type_free(&single_remote_point_dt), comm);

  if (weight_flag)
    relocate_weights(
      weights, reorder_idx, ranks, count_,
      send_count, local_count, recv_count,
      sendcounts, sdispls + 1, recvcounts, rdispls, comm);

  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  free(single_remote_point_buffer);
  free(reorder_idx);
  free(ranks);

  *count = local_count + recv_count;
}
