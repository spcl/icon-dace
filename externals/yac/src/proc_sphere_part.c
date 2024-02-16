/**
 * @file proc_sphere_part.c
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

#include "string.h"
#include "proc_sphere_part.h"
#include "yac_mpi.h"
#ifdef SCOREP_USER_ENABLE
#include "scorep/SCOREP_User.h"
#endif // SCOREP_USER_ENABLE

enum splicomm_tags {
  DATA_SIZE_TAG,
  DATA_TAG,
};

struct proc_sphere_part_node;

struct proc_sphere_part_node_data {
  union {
    struct proc_sphere_part_node * node;
    int rank;
  } data;
  int is_leaf;
};
struct proc_sphere_part_node {
  struct proc_sphere_part_node_data U, T;
  double gc_norm_vector[3];
};

struct neigh_search_data {
  int * ranks;
  int num_ranks;
  struct proc_sphere_part_node * node;
};

static int get_proc_sphere_part_node_data_pack_size(
  struct proc_sphere_part_node_data node_data, MPI_Comm comm);

static int get_proc_sphere_part_node_pack_size(
  struct proc_sphere_part_node node, MPI_Comm comm) {

  int vec_pack_size;
  yac_mpi_call(
    MPI_Pack_size(3, MPI_DOUBLE, comm, &vec_pack_size), comm);

  return vec_pack_size +
         get_proc_sphere_part_node_data_pack_size(node.U, comm) +
         get_proc_sphere_part_node_data_pack_size(node.T, comm);
}

static int get_proc_sphere_part_node_data_pack_size(
  struct proc_sphere_part_node_data node_data, MPI_Comm comm) {

  int int_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &int_pack_size), comm);

  int data_size = int_pack_size;
  if (node_data.is_leaf)
    data_size += int_pack_size;
  else
    data_size +=
      get_proc_sphere_part_node_pack_size((*node_data.data.node), comm);

  return data_size;
}

static void pack_proc_sphere_part_node_data(
  struct proc_sphere_part_node_data node_data, void * pack_buffer,
  int pack_buffer_size, int * position, MPI_Comm comm);

static void pack_proc_sphere_part_node(
  struct proc_sphere_part_node * node, void * pack_buffer,
  int pack_buffer_size, int * position, MPI_Comm comm) {

  yac_mpi_call(MPI_Pack(&(node->gc_norm_vector[0]), 3, MPI_DOUBLE, pack_buffer,
                        pack_buffer_size, position, comm), comm);
  pack_proc_sphere_part_node_data(
    node->U, pack_buffer, pack_buffer_size, position, comm);
  pack_proc_sphere_part_node_data(
    node->T, pack_buffer, pack_buffer_size, position, comm);
}

static void pack_proc_sphere_part_node_data(
  struct proc_sphere_part_node_data node_data, void * pack_buffer,
  int pack_buffer_size, int * position, MPI_Comm comm) {

  yac_mpi_call(MPI_Pack(&(node_data.is_leaf), 1, MPI_INT, pack_buffer,
                        pack_buffer_size, position, comm), comm);

  if (node_data.is_leaf)
    yac_mpi_call(MPI_Pack(&(node_data.data.rank), 1, MPI_INT, pack_buffer,
                          pack_buffer_size, position, comm), comm);
  else
    pack_proc_sphere_part_node(node_data.data.node, pack_buffer,
                               pack_buffer_size, position, comm);
}

static struct proc_sphere_part_node_data unpack_proc_sphere_part_node_data(
  void * pack_buffer, int pack_buffer_size, int * position,
  MPI_Comm comm);

static struct proc_sphere_part_node * unpack_proc_sphere_part_node(
  void * pack_buffer, int pack_buffer_size, int * position,
  MPI_Comm comm) {

  struct proc_sphere_part_node * node = xmalloc(1 * sizeof(*node));

  yac_mpi_call(
    MPI_Unpack(pack_buffer, pack_buffer_size, position,
               &(node->gc_norm_vector[0]), 3, MPI_DOUBLE, comm), comm);

  node->U = unpack_proc_sphere_part_node_data(pack_buffer, pack_buffer_size,
                                              position, comm);
  node->T = unpack_proc_sphere_part_node_data(pack_buffer, pack_buffer_size,
                                              position, comm);

  return node;
}

static struct proc_sphere_part_node_data unpack_proc_sphere_part_node_data(
  void * pack_buffer, int pack_buffer_size, int * position,
  MPI_Comm comm) {

  struct proc_sphere_part_node_data node_data;

  yac_mpi_call(
    MPI_Unpack(pack_buffer, pack_buffer_size, position,
               &(node_data.is_leaf), 1, MPI_INT, comm), comm);

  if (node_data.is_leaf)
    yac_mpi_call(
      MPI_Unpack(pack_buffer, pack_buffer_size, position,
                 &(node_data.data.rank), 1, MPI_INT, comm), comm);
  else
    node_data.data.node =
      unpack_proc_sphere_part_node(
        pack_buffer, pack_buffer_size, position, comm);

  return node_data;
}

static struct proc_sphere_part_node_data get_remote_data(
  struct proc_sphere_part_node_data local_data,
  struct yac_group_comm local_group_comm,
  struct yac_group_comm remote_group_comm) {

  int comm_rank;
  MPI_Comm comm = local_group_comm.comm;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  int data_size;
  void * recv_buffer;

  int order = local_group_comm.start < remote_group_comm.start;

  for (int i = 0; i < 2; ++i) {

    if (order == i) {

      yac_bcast_group(
        &data_size, 1, MPI_INT, remote_group_comm.start, local_group_comm);
      recv_buffer = xmalloc((size_t)data_size);
      yac_bcast_group(recv_buffer, data_size, MPI_PACKED,
                      remote_group_comm.start, local_group_comm);

    } else {
      if (comm_rank == local_group_comm.start) {

        // pack local_data
        int pack_buffer_size =
          get_proc_sphere_part_node_data_pack_size(local_data, comm);
        void * pack_buffer = xmalloc((size_t)pack_buffer_size);
        int position = 0;
        pack_proc_sphere_part_node_data(
          local_data, pack_buffer, pack_buffer_size, &position, comm);

        // broadcast data size to other group
        yac_bcast_group(
          &position, 1, MPI_INT, local_group_comm.start, remote_group_comm);

        // broadcast remote_data to other root
        yac_bcast_group(pack_buffer, position, MPI_PACKED,
                        local_group_comm.start, remote_group_comm);
        free(pack_buffer);
      }
    }
  }

  // unpack node data
  int position = 0;
  struct proc_sphere_part_node_data other_data =
    unpack_proc_sphere_part_node_data(
      recv_buffer, data_size, &position, comm);
  free(recv_buffer);

  return other_data;
}

static void compute_redist_recvcounts_rdispls(
  uint64_t comm_rank, uint64_t split_rank, uint64_t comm_size,
  uint64_t global_sizes[2], uint64_t (*all_bucket_sizes)[2],
  int * counts, int * displs, size_t * recv_count) {

  int color = comm_rank >= split_rank;

  uint64_t split_comm_size = split_rank;
  uint64_t split_comm_rank = comm_rank;
  if (color) {
    split_comm_rank -= split_comm_size;
    split_comm_size = comm_size - split_comm_size;
  }

  uint64_t global_size = global_sizes[color];
  uint64_t local_interval_start =
    (global_size * split_comm_rank + split_comm_size - 1) /
    split_comm_size;
  uint64_t local_interval_end =
    (global_size * (split_comm_rank+1) + split_comm_size - 1) /
    split_comm_size;

  *recv_count = (size_t)(local_interval_end - local_interval_start);

  for (uint64_t i = 0, start_idx = 0; i < comm_size; ++i) {

    uint64_t next_start_idx = start_idx + all_bucket_sizes[i][color];
    uint64_t interval_start = MAX(start_idx, local_interval_start);
    uint64_t interval_end = MIN(next_start_idx, local_interval_end);

    if (interval_start < interval_end) {

      uint64_t count = interval_end - interval_start;
      uint64_t disp = interval_start - local_interval_start;

      YAC_ASSERT(
        (count <= INT_MAX) && (disp <= INT_MAX),
        "ERROR(compute_redist_recvcounts_rdispls): invalid interval")

      counts[i] = (int)count;
      displs[i] = (int)disp;
    } else {
      counts[i] = 0;
      displs[i] = 0;
    }

    start_idx = next_start_idx;
  }
}

static void compute_redist_sendcounts_sdispls(
  uint64_t comm_rank, uint64_t split_rank, uint64_t comm_size,
  uint64_t global_sizes[2], uint64_t (*all_bucket_sizes)[2],
  uint64_t U_size, int * counts, int * displs) {

  uint64_t local_interval_start[2] = {0, 0};
  uint64_t local_interval_end[2];
  for (uint64_t i = 0; i < comm_rank; ++i)
    for (int j = 0; j < 2; ++j)
      local_interval_start[j] += all_bucket_sizes[i][j];
  for (int j = 0; j < 2; ++j)
    local_interval_end[j] =
      local_interval_start[j] + all_bucket_sizes[comm_rank][j];

  uint64_t comm_sizes[2] = {split_rank, comm_size - split_rank};

  for (uint64_t i = 0, start_idx[2] = {0,0}; i < comm_size; ++i) {

    int color = i >= split_rank;
    uint64_t global_size = global_sizes[color];
    uint64_t split_comm_rank = i - (color?(split_rank):0);
    uint64_t split_comm_size = comm_sizes[color];
    uint64_t next_start_idx =
      (global_size * (split_comm_rank + 1) + split_comm_size - 1) /
      split_comm_size;
    uint64_t interval_start = MAX(start_idx[color], local_interval_start[color]);
    uint64_t interval_end = MIN(next_start_idx, local_interval_end[color]);

    if (interval_start < interval_end) {

      uint64_t count = interval_end - interval_start;
      uint64_t disp = interval_start - local_interval_start[color] +
                      ((color)?(U_size):(0));

      YAC_ASSERT(
        (count <= INT_MAX) && (disp <= INT_MAX),
        "ERROR(compute_redist_sendcounts_sdispls): invalid interval")

      counts[i] = (int)count;
      displs[i] = (int)disp;
    } else {
      counts[i] = 0;
      displs[i] = 0;
    }

    start_idx[color] = next_start_idx;
  }
}

static int compare_dist_cells_coord(const void * a, const void * b) {

  return memcmp(&(((const struct dist_cell *)a)->coord[0]),
                &(((const struct dist_cell *)b)->coord[0]), 3 * sizeof(double));
}

static void remove_duplicated_cells(
  struct dist_cell * cells, size_t * num_cells) {

  if (*num_cells == 0) return;

  //------------------------------------
  // sort cells by grid id and global id
  //------------------------------------

  qsort(cells, *num_cells, sizeof(*cells), compare_dist_cells_coord);

  size_t old_num_cells = *num_cells;
  struct dist_cell * last_new_cell = cells;
  for (size_t i = 1; i < old_num_cells; ++i) {
    struct dist_cell * curr_cell = cells + i;
    if (compare_dist_cells_coord(last_new_cell, curr_cell)) {
      ++last_new_cell;
      if (last_new_cell != curr_cell) *last_new_cell = *curr_cell;
    }
  }
  *num_cells = (size_t)(last_new_cell - cells) + 1;
}

static struct proc_sphere_part_node * redistribute_cells_recursive(
  struct dist_cell ** cells, size_t * num_cells,
  struct yac_group_comm group_comm, double prev_gc_norm_vector[3],
  MPI_Datatype dist_cell_dt) {

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_DEFINE( local_balance_point_region )
SCOREP_USER_REGION_DEFINE( global_balance_point_region )
SCOREP_USER_REGION_DEFINE( splitting_region )
SCOREP_USER_REGION_DEFINE( comm_split_region )
SCOREP_USER_REGION_DEFINE( redist_data_region )
#endif

  int group_rank = yac_group_comm_get_rank(group_comm);
  int group_size = yac_group_comm_get_size(group_comm);

  remove_duplicated_cells(*cells, num_cells);

  //--------------------
  // compute split plane
  //--------------------
#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_BEGIN(
  local_balance_point_region, "local balance point",
  SCOREP_USER_REGION_TYPE_COMMON )
#endif
  // compute local balance point
  double balance_point[3] = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < *num_cells; ++i) {
    double * cell_coord = (*cells)[i].coord;
    for (int j = 0; j < 3; ++j) balance_point[j] += cell_coord[j];
  }
#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_END( local_balance_point_region )
SCOREP_USER_REGION_BEGIN(
  global_balance_point_region, "global balance point",
  SCOREP_USER_REGION_TYPE_COMMON )
#endif

  // in case *num_cells is zero or something else happens where the coordinates
  // of the balance point are zero, the normalisation would generate NaN's
  if ((balance_point[0] != 0.0) ||
      (balance_point[1] != 0.0) ||
      (balance_point[2] != 0.0))
    normalise_vector(balance_point);

  // compute global balance point (make sure that the allreduce operation
  // generates bit-identical results on all processes)
  yac_allreduce_sum_dble(&(balance_point[0]), 3, group_comm);

  if ((fabs(balance_point[0]) > 1e-9) ||
      (fabs(balance_point[1]) > 1e-9) ||
      (fabs(balance_point[2]) > 1e-9)) {
     normalise_vector(balance_point);
  } else {
     balance_point[0] = prev_gc_norm_vector[2];
     balance_point[1] = prev_gc_norm_vector[0];
     balance_point[2] = prev_gc_norm_vector[1];
  }

  double gc_norm_vector[3];
  crossproduct_kahan(balance_point, prev_gc_norm_vector, gc_norm_vector);

  if ((fabs(gc_norm_vector[0]) > 1e-9) ||
      (fabs(gc_norm_vector[1]) > 1e-9) ||
      (fabs(gc_norm_vector[2]) > 1e-9)) {
     normalise_vector(gc_norm_vector);
  } else {
     gc_norm_vector[0] = prev_gc_norm_vector[2];
     gc_norm_vector[1] = prev_gc_norm_vector[0];
     gc_norm_vector[2] = prev_gc_norm_vector[1];
  }

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_END( global_balance_point_region )
#endif

  //-----------------
  // split local data
  //-----------------

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_BEGIN( splitting_region, "splitting data", SCOREP_USER_REGION_TYPE_COMMON )
#endif

  // angle between a cell coord and the great circle plane:
  // acos(dot(gc_norm_vector, cell_coord)) = angle(gc_norm_vector, cell_coord)
  // acos(dot(gc_norm_vector, cell_coord)) - PI/2 = angle(gc_plane, cell_coord)
  // dot <= 0.0    -> U list
  // dot >  0.0    -> T list

  struct dist_cell * left = *cells, * right = *cells + *num_cells - 1;

  // sort such that all cells for the U list come first, followed by the
  // cells of the T list
  while (1) {
    // find a cell that does not belong into U-list
    while (left <= right) {
      double * curr_coordinates_xyz = left->coord;
      double dot = curr_coordinates_xyz[0] * gc_norm_vector[0] +
                   curr_coordinates_xyz[1] * gc_norm_vector[1] +
                   curr_coordinates_xyz[2] * gc_norm_vector[2];

      // if (angle < M_PI_2)
      if (dot > 0.0) break;
      ++left;
    };

    // find a cell that does not belong into T-list
    while (left < right) {
      double * curr_coordinates_xyz = right->coord;
      double dot = curr_coordinates_xyz[0] * gc_norm_vector[0] +
                   curr_coordinates_xyz[1] * gc_norm_vector[1] +
                   curr_coordinates_xyz[2] * gc_norm_vector[2];

      // if (angle >= M_PI_2)
      if (dot <= 0.0) break;
      --right;
    }

    if (left < right) {
      struct dist_cell tmp_cell = *left;
      *left = *right;
      *right = tmp_cell;
      ++left;
      --right;
    } else {
      break;
    }
  }


  uint64_t U_size = (uint64_t)(left - *cells);
  uint64_t T_size = (uint64_t)(*num_cells) - U_size;

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_END( splitting_region )
#endif

  // there is no MPI datatype for size_t -> we use uint64_t here
  uint64_t bucket_sizes[2] = {U_size, T_size};

  // exchange local U/T sizes between all processes
  uint64_t (*all_bucket_sizes)[2] =
    xmalloc((size_t)group_size * sizeof(*all_bucket_sizes));
  yac_allgather_uint64(
    &(bucket_sizes[0]), &(all_bucket_sizes[0][0]), 2, group_comm);

  // determine global U/T sizes
  uint64_t global_bucket_sizes[2] = {0, 0};
  for (int i = 0; i < group_size; ++i)
    for (int j = 0; j < 2; ++j)
      global_bucket_sizes[j] += all_bucket_sizes[i][j];
  uint64_t global_num_cells = global_bucket_sizes[0] + global_bucket_sizes[1];

  //----------------------
  // split into two groups
  //----------------------
#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_BEGIN(
  comm_split_region, "creating splitcomm", SCOREP_USER_REGION_TYPE_COMMON )
#endif
  // determine processor groups
  int split_rank =
    MIN(
      (int)MAX(
        ((global_bucket_sizes[0] * (uint64_t)group_size + global_num_cells/2) /
         global_num_cells), 1),
      group_size - 1);

  // generate processor groups
  struct yac_group_comm local_group_comm, remote_group_comm;
  yac_group_comm_split(
    group_comm, split_rank, &local_group_comm, &remote_group_comm);

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_END( comm_split_region )
#endif
  //------------------
  // redistribute data
  //------------------

  int * int_buffer = xmalloc(4 * (size_t)group_size * sizeof(*int_buffer));
  int * sendcounts = int_buffer + 0 * group_size;
  int * recvcounts = int_buffer + 1 * group_size;
  int * sdispls =    int_buffer + 2 * group_size;
  int * rdispls =    int_buffer + 3 * group_size;

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_BEGIN(
  redist_data_region, "data redistribution", SCOREP_USER_REGION_TYPE_COMMON )
#endif
  // compute send and receive counts and respective ranks for data
  // redistribution
  compute_redist_sendcounts_sdispls(
    (uint64_t)group_rank, (uint64_t)split_rank, (uint64_t)group_size,
    global_bucket_sizes, all_bucket_sizes, (uint64_t)U_size, sendcounts, sdispls);
  size_t new_num_cells;
  compute_redist_recvcounts_rdispls(
    (uint64_t)group_rank, (uint64_t)split_rank, (uint64_t)group_size,
    global_bucket_sizes, all_bucket_sizes, recvcounts, rdispls, &new_num_cells);

  // exchange data between both groups
  struct dist_cell * new_cells = xmalloc(new_num_cells * sizeof(*new_cells));
  yac_alltoallv_p2p_group(
    *cells, sendcounts, sdispls, new_cells, recvcounts, rdispls,
    sizeof(**cells), dist_cell_dt, group_comm);
#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_END( redist_data_region )
#endif
  free(*cells);
  *cells = new_cells;
  *num_cells = new_num_cells;

  free(int_buffer);
  free(all_bucket_sizes);

  //----------
  // recursion
  //----------

  // redistribute data within local group
  struct proc_sphere_part_node_data local_data;

  if (yac_group_comm_get_size(local_group_comm) > 1) {
    local_data.data.node =
      redistribute_cells_recursive(
        cells, num_cells, local_group_comm, gc_norm_vector, dist_cell_dt);
    local_data.is_leaf = 0;
  } else {

    remove_duplicated_cells(*cells, num_cells);
    local_data.data.rank = yac_group_comm_get_global_rank(group_comm);
    local_data.is_leaf = 1;
  }

  // get proc_sphere_part_node_data from remote group
  struct proc_sphere_part_node_data remote_data =
    get_remote_data(local_data, local_group_comm, remote_group_comm);

  // generate node
  struct proc_sphere_part_node * node = xmalloc(1 * sizeof(*node));
  if (group_rank < split_rank) {
    node->U = local_data;
    node->T = remote_data;
  } else {
    node->U = remote_data;
    node->T = local_data;
  }
  for (int i = 0; i < 3; ++i) node->gc_norm_vector[i] = gc_norm_vector[i];

  return node;
}

static MPI_Datatype yac_get_dist_cell_mpi_datatype(MPI_Comm comm) {

  MPI_Datatype coord_dt;
  yac_mpi_call(MPI_Type_contiguous(3, MPI_DOUBLE, &coord_dt), comm);
  return yac_create_resized(coord_dt, sizeof(struct dist_cell), comm);
}

struct proc_sphere_part_node * yac_redistribute_cells(
  struct dist_cell ** cells, size_t * num_cells, MPI_Comm comm) {

  double base_gc_norm_vector[3] = {0.0,0.0,1.0};

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  uint64_t global_num_cells = (uint64_t)*num_cells;
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &global_num_cells, 1, MPI_UINT64_T, MPI_SUM, comm), comm);

  struct proc_sphere_part_node * proc_sphere_part;

  if ((comm_size > 1) && (global_num_cells > 0)) {

    MPI_Datatype dist_cell_dt = yac_get_dist_cell_mpi_datatype(comm);
    struct yac_group_comm group_comm = yac_group_comm_new(comm);
    proc_sphere_part =
      redistribute_cells_recursive(cells, num_cells, group_comm,
                                   base_gc_norm_vector, dist_cell_dt);
    yac_group_comm_delete(group_comm);

    yac_mpi_call(MPI_Type_free(&dist_cell_dt), comm);
  } else {
    proc_sphere_part = xmalloc(1 * sizeof(*proc_sphere_part));
    proc_sphere_part->U.data.rank = 0;
    proc_sphere_part->U.is_leaf = 1;
    proc_sphere_part->T.data.rank = 0;
    proc_sphere_part->T.is_leaf = 1;
    proc_sphere_part->gc_norm_vector[0] = base_gc_norm_vector[0];
    proc_sphere_part->gc_norm_vector[1] = base_gc_norm_vector[1];
    proc_sphere_part->gc_norm_vector[2] = base_gc_norm_vector[2];
  }

  return proc_sphere_part;
}

static int is_serial_node(struct proc_sphere_part_node * node) {
  return (node->U.is_leaf) && (node->T.is_leaf) &&
         (node->U.data.rank == 0) && (node->T.data.rank == 0);
}

void yac_proc_sphere_part_do_point_search(
  struct proc_sphere_part_node * node, coordinate_pointer search_coords,
  size_t count, int * ranks) {

  if (is_serial_node(node)) {
    for (size_t i = 0; i < count; ++i) ranks[i] = 0;
    return;
  }

  for (size_t i = 0; i < count; ++i) {

    double curr_coord[3] =
      {search_coords[i][0], search_coords[i][1], search_coords[i][2]};

    struct proc_sphere_part_node * curr_node = node;

    while (1) {

      // get angle between the norm vector of the split plane and the base
      // point of the bounding circle
      struct sin_cos_angle angle =
        get_vector_angle_2(curr_coord, curr_node->gc_norm_vector);

      if (angle.cos <= 0.0) {
        if (curr_node->U.is_leaf) {
          ranks[i] = curr_node->U.data.rank;
          break;
        } else {
          curr_node = curr_node->U.data.node;
          continue;
        }
      } else {
        if (curr_node->T.is_leaf) {
          ranks[i] = curr_node->T.data.rank;
          break;
        } else {
          curr_node = curr_node->T.data.node;
          continue;
        }
      }
    }
  }
}

static void bnd_circle_search(
  struct proc_sphere_part_node * node, struct bounding_circle bnd_circle,
  int * ranks, int * rank_count) {

  double dot = bnd_circle.base_vector[0] * node->gc_norm_vector[0] +
               bnd_circle.base_vector[1] * node->gc_norm_vector[1] +
               bnd_circle.base_vector[2] * node->gc_norm_vector[2];

  // angle < M_PI_2 + bnd_circle.inc_angle
  if (dot > - bnd_circle.inc_angle.sin) {

    if (node->T.is_leaf) {

      ranks[*rank_count] = node->T.data.rank;
      ++*rank_count;

    } else {
      bnd_circle_search(node->T.data.node, bnd_circle, ranks, rank_count);
    }
  }

  // angle > M_PI_2 - bnd_circle.inc_angle
  if (dot < bnd_circle.inc_angle.sin) {

    if (node->U.is_leaf) {

      ranks[*rank_count] = node->U.data.rank;
      ++*rank_count;

    } else {
      bnd_circle_search(node->U.data.node, bnd_circle, ranks, rank_count);
    }
  }
}

void yac_proc_sphere_part_do_bnd_circle_search(
  struct proc_sphere_part_node * node, struct bounding_circle bnd_circle,
  int * ranks, int * rank_count) {

  YAC_ASSERT(
    compare_angles(bnd_circle.inc_angle, SIN_COS_M_PI) == -1,
    "ERROR(yac_proc_sphere_part_do_bnd_circle_search): angle is >= PI")

  // special case in which the proc_sphere_part_node only contains a single
  // rank
  if (is_serial_node(node)) {
    ranks[0] = 0;
    *rank_count = 1;
  } else {
    *rank_count = 0;
    bnd_circle_search(node, bnd_circle, ranks, rank_count);
  }
}

static int get_leaf_ranks(struct proc_sphere_part_node * node, int * ranks) {

  int curr_size;

  if (node->U.is_leaf) {
    *ranks = node->U.data.rank;
    curr_size = 1;
  } else {
    curr_size = get_leaf_ranks(node->U.data.node, ranks);
  }
  if (node->T.is_leaf) {
    ranks[curr_size++] = node->T.data.rank;
  } else {
    curr_size += get_leaf_ranks(node->T.data.node, ranks + curr_size);
  }

  return curr_size;
}

static void get_neigh_ranks(
  struct proc_sphere_part_node * node, uint64_t * leaf_sizes, uint64_t min_size,
  uint64_t ** inner_node_sizes, int * send_flags, int * recv_flags,
  int comm_rank, struct neigh_search_data * last_valid_node) {

  int curr_node_is_valid = (**inner_node_sizes >= min_size);

  struct neigh_search_data * curr_valid_node;
  struct neigh_search_data temp_valid_node;

  int * ranks = last_valid_node->ranks;

  if (curr_node_is_valid) {
    temp_valid_node.ranks = ranks;
    temp_valid_node.num_ranks = 0;
    temp_valid_node.node = node;
    curr_valid_node = &temp_valid_node;
    last_valid_node->num_ranks = 0;
  } else {
    curr_valid_node = last_valid_node;
  }

  for (int j = 0; j < 2; ++j) {

    struct proc_sphere_part_node_data node_data = (j == 0)?(node->U):(node->T);

    if (node_data.is_leaf) {

      int rank = node_data.data.rank;

      if (leaf_sizes[rank] < min_size) {

        if (curr_valid_node->num_ranks == 0)
          curr_valid_node->num_ranks =
            get_leaf_ranks(curr_valid_node->node, ranks);

        // if the current leaf is the local process
        if (rank == comm_rank)
          for (int i = 0; i < curr_valid_node->num_ranks; ++i)
            recv_flags[ranks[i]] = 1;

        // if the process of the current leaf required data from the
        // local process
        for (int i = 0; i < curr_valid_node->num_ranks; ++i) {
          if (ranks[i] == comm_rank) {
            send_flags[rank] = 1;
            break;
          }
        }
      }

    } else {
      ++*inner_node_sizes;
      get_neigh_ranks(
        node_data.data.node, leaf_sizes, min_size, inner_node_sizes,
        send_flags, recv_flags, comm_rank, curr_valid_node);
    }
  }
}

static uint64_t determine_node_sizes(
  struct proc_sphere_part_node * node, uint64_t * leaf_sizes,
  uint64_t ** inner_node_sizes) {

  uint64_t * curr_inner_node_size = *inner_node_sizes;
  uint64_t node_size;
  if (node->U.is_leaf) {
    node_size = leaf_sizes[node->U.data.rank];
  } else {
    ++*inner_node_sizes;
    node_size =
      determine_node_sizes(node->U.data.node, leaf_sizes, inner_node_sizes);
  }
  if (node->T.is_leaf) {
    node_size += leaf_sizes[node->T.data.rank];
  } else {
    ++*inner_node_sizes;
    node_size +=
      determine_node_sizes(node->T.data.node, leaf_sizes, inner_node_sizes);
  }
  return (*curr_inner_node_size = node_size);
}

void yac_proc_sphere_part_get_neigh_ranks(
  struct proc_sphere_part_node * node, uint64_t * leaf_sizes,
  uint64_t min_size, int * send_flags, int * recv_flags,
  int comm_rank, int comm_size) {

  uint64_t * inner_node_sizes =
    xcalloc((size_t)comm_size, sizeof(*inner_node_sizes));

  uint64_t * temp_inner_node_sizes = inner_node_sizes;
  determine_node_sizes(node, leaf_sizes, &temp_inner_node_sizes);

  YAC_ASSERT(
    *inner_node_sizes >= min_size,
    "ERROR(yac_proc_sphere_part_get_neigh_ranks): sum of global leaf sizes "
    "is < min_size")

  struct neigh_search_data search_data = {
    .ranks = xmalloc((size_t)comm_size * sizeof(int)),
    .num_ranks = 0,
    .node = node
  };

  temp_inner_node_sizes = inner_node_sizes;
  get_neigh_ranks(node, leaf_sizes, min_size, &temp_inner_node_sizes,
                  send_flags, recv_flags, comm_rank, &search_data);

  send_flags[comm_rank] = 0;
  recv_flags[comm_rank] = 0;

  free(search_data.ranks);
  free(inner_node_sizes);
}

void yac_proc_sphere_part_node_delete(struct proc_sphere_part_node * node) {

  if (!(node->U.is_leaf)) yac_proc_sphere_part_node_delete(node->U.data.node);
  if (!(node->T.is_leaf)) yac_proc_sphere_part_node_delete(node->T.data.node);
  free(node);
}
