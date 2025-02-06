// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>
#include <mpi.h>

#include "interp_weights_internal.h"
#include "basic_grid.h"
#include "utils_common.h"
#include "yac_mpi_internal.h"

struct dup_info {
  yac_int orig_global_id;
  yac_int duplicated_global_id;
  size_t duplicated_orig_pos;
  int duplicated_rank;
};

struct stencil_info {
  yac_int tgt_global_id;
  size_t idx;
  int rank;
};

static MPI_Datatype get_dup_info_mpi_datatype(MPI_Comm comm) {

  struct dup_info dummy;
  MPI_Datatype dup_info_dt;
  int array_of_blocklengths[] = {1, 1, 1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.orig_global_id) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.duplicated_global_id) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.duplicated_orig_pos) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.duplicated_rank) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {yac_int_dt, yac_int_dt, YAC_MPI_SIZE_T, MPI_INT};
  yac_mpi_call(
    MPI_Type_create_struct(4, array_of_blocklengths, array_of_displacements,
                           array_of_types, &dup_info_dt), comm);
  return yac_create_resized(dup_info_dt, sizeof(dummy), comm);
}

static MPI_Datatype get_stencil_info_mpi_datatype(MPI_Comm comm) {

  struct stencil_info dummy;
  MPI_Datatype stencil_info_dt;
  int array_of_blocklengths[] = {1, 1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.tgt_global_id) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.idx) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.rank) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {yac_int_dt, YAC_MPI_SIZE_T, MPI_INT};
  yac_mpi_call(
    MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements,
                           array_of_types, &stencil_info_dt), comm);
  return yac_create_resized(stencil_info_dt, sizeof(dummy), comm);
}

static int compare_dup_info_orig_global_ids(
  const void * a, const void * b) {

  return (((const struct dup_info *)a)->orig_global_id >
          ((const struct dup_info *)b)->orig_global_id) -
         (((const struct dup_info *)a)->orig_global_id <
          ((const struct dup_info *)b)->orig_global_id);
}

static int compare_stencil_info_tgt_global_id(
  const void * a, const void * b) {

  return (((const struct stencil_info *)a)->tgt_global_id >
          ((const struct stencil_info *)b)->tgt_global_id) -
         (((const struct stencil_info *)a)->tgt_global_id <
          ((const struct stencil_info *)b)->tgt_global_id);
}

static inline int compute_bucket(yac_int value, int comm_size) {
  return (int)(value / 128) % comm_size;
}

void yac_duplicate_stencils(
  struct yac_interp_weights * weights, struct yac_basic_grid * tgt_grid,
  yac_int * tgt_orig_global_id, size_t * tgt_duplicated_idx,
  size_t nbr_duplicated, enum yac_location location) {

  YAC_ASSERT(
    location == YAC_LOC_CELL,
    "ERROR(yac_duplicate_stencils): only supported for cells");

  yac_int * duplicated_global_ids =
    xmalloc(nbr_duplicated * sizeof(*duplicated_global_ids));

  struct yac_basic_grid_data * grid_data =
    yac_basic_grid_get_data(tgt_grid);
  yac_int const * cell_ids = grid_data->cell_ids;

  // get global ids of all duplicated target cells
  for (size_t i = 0; i < nbr_duplicated; ++i)
    duplicated_global_ids[i] = cell_ids[tgt_duplicated_idx[i]];

  // get target global ids of all local stencils
  size_t stencil_count = yac_interp_weights_get_interp_count(weights);
  yac_int * stencil_tgt_ids =
    yac_interp_weights_get_interp_tgt(weights);

  MPI_Comm comm = yac_interp_weights_get_comm(weights);
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int * int_buffer =
    xmalloc((nbr_duplicated + stencil_count) * sizeof(*int_buffer));
  int * dup_dist_ranks = int_buffer;
  int * stencil_dist_ranks = int_buffer + nbr_duplicated;

  // determine the bucket ranks
  for (size_t i = 0; i < nbr_duplicated; ++i)
    dup_dist_ranks[i] = compute_bucket(tgt_orig_global_id[i], comm_size);
  for (size_t i = 0; i < stencil_count; ++i)
    stencil_dist_ranks[i] = compute_bucket(stencil_tgt_ids[i], comm_size);

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  // set up counts and displs for redistribution of information about
  // duplicated points
  for (size_t i = 0; i < nbr_duplicated; ++i)
    sendcounts[dup_dist_ranks[i]]++;
  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t dup_info_count = recvcounts[comm_size - 1] + rdispls[comm_size - 1];
  struct dup_info * dup_info_buffer =
    xmalloc((nbr_duplicated + dup_info_count) * sizeof(*dup_info_buffer));
  struct dup_info * dup_send_buffer = dup_info_buffer + dup_info_count;
  struct dup_info * dup_recv_buffer = dup_info_buffer;

  // pack information about duplicated points
  for (size_t i = 0; i < nbr_duplicated; ++i) {
    size_t pos = sdispls[dup_dist_ranks[i]+1]++;
    dup_send_buffer[pos].orig_global_id = tgt_orig_global_id[i];
    dup_send_buffer[pos].duplicated_global_id = duplicated_global_ids[i];
    dup_send_buffer[pos].duplicated_orig_pos = tgt_duplicated_idx[i];
    dup_send_buffer[pos].duplicated_rank = comm_rank;
  }

  // redistribute information about duplicated points
  MPI_Datatype dup_info_dt = get_dup_info_mpi_datatype(comm);
  yac_mpi_call(MPI_Type_commit(&dup_info_dt), comm);
  yac_alltoallv_p2p(
    dup_send_buffer, sendcounts, sdispls,
    dup_recv_buffer, recvcounts, rdispls,
    sizeof(*dup_send_buffer), dup_info_dt, comm);
  yac_mpi_call(MPI_Type_free(&dup_info_dt), comm);
  dup_recv_buffer =
    xrealloc(dup_recv_buffer, dup_info_count * sizeof(*dup_recv_buffer));
  qsort(
    dup_recv_buffer, dup_info_count, sizeof(*dup_recv_buffer),
    compare_dup_info_orig_global_ids);

  // set up counts and displs for redistribution of stencil information
  memset(sendcounts, 0, (size_t)comm_size * sizeof(*sendcounts));
  for (size_t i = 0; i < stencil_count; ++i)
    sendcounts[stencil_dist_ranks[i]]++;
  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t stencil_info_count = recvcounts[comm_size - 1] + rdispls[comm_size - 1];
  struct stencil_info * stencil_info_buffer =
    xmalloc((stencil_count + stencil_info_count) * sizeof(*stencil_info_buffer));
  struct stencil_info * stencil_send_buffer =
    stencil_info_buffer + stencil_info_count;
  struct stencil_info * stencil_recv_buffer = stencil_info_buffer;

  // pack stencil information
  for (size_t i = 0; i < stencil_count; ++i) {
    size_t pos = sdispls[stencil_dist_ranks[i]+1]++;
    stencil_send_buffer[pos].tgt_global_id = stencil_tgt_ids[i];
    stencil_send_buffer[pos].idx = i;
    stencil_send_buffer[pos].rank = comm_rank;
  }
  free(stencil_tgt_ids);

  // redistribute stencil information
  MPI_Datatype stencil_info_dt = get_stencil_info_mpi_datatype(comm);
  yac_mpi_call(MPI_Type_commit(&stencil_info_dt), comm);
  yac_alltoallv_p2p(
    stencil_send_buffer, sendcounts, sdispls,
    stencil_recv_buffer, recvcounts, rdispls,
    sizeof(*stencil_send_buffer), stencil_info_dt, comm);
  yac_mpi_call(MPI_Type_free(&stencil_info_dt), comm);
  stencil_recv_buffer =
    xrealloc(
      stencil_recv_buffer, stencil_info_count * sizeof(*stencil_recv_buffer));
  qsort(
    stencil_recv_buffer, stencil_info_count, sizeof(*stencil_recv_buffer),
    compare_stencil_info_tgt_global_id);

  // count number of duplicated targets for which at least one stencil exists
  // and the number of stencils that match with a duplicated points
  size_t dup_match_count = 0;
  size_t stencil_match_count = 0;
  for (size_t i = 0, j = 0; i < dup_info_count; ++i) {
    yac_int curr_orig_global_id = dup_recv_buffer[i].orig_global_id;
    while ((j < stencil_info_count) &&
           (stencil_recv_buffer[j].tgt_global_id < curr_orig_global_id)) ++j;
    if ((j < stencil_info_count) &&
        (stencil_recv_buffer[j].tgt_global_id == curr_orig_global_id))
      ++dup_match_count;
    size_t k = j;
    while ((k < stencil_info_count) &&
           (stencil_recv_buffer[k].tgt_global_id == curr_orig_global_id)) {
      ++k;
      ++stencil_match_count;
    }
  }

  struct remote_points dup_remote_points;
  dup_remote_points.data =
    xmalloc(dup_match_count * sizeof(*(dup_remote_points.data)));
  dup_remote_points.count = dup_match_count;
  size_t * num_stencils_per_tgt =
    xmalloc(stencil_match_count * sizeof(*num_stencils_per_tgt));
  size_t * stencil_indices =
    xmalloc(stencil_match_count * sizeof(*stencil_indices));
  int * stencil_ranks =
    xmalloc(stencil_match_count * sizeof(*stencil_ranks));
  double * w = xmalloc(stencil_match_count * sizeof(*w));

  // match duplicated targets with stencils
  dup_match_count = 0;
  stencil_match_count = 0;
  for (size_t i = 0, j = 0; i < dup_info_count; ++i) {
    yac_int curr_orig_global_id = dup_recv_buffer[i].orig_global_id;

    // search for a matching stencil
    while ((j < stencil_info_count) &&
           (stencil_recv_buffer[j].tgt_global_id < curr_orig_global_id)) ++j;

    // if we found a matching stencil
    if ((j < stencil_info_count) &&
        (stencil_recv_buffer[j].tgt_global_id == curr_orig_global_id)) {

      size_t curr_num_stencils_per_tgt = 0;

      // for the current duplicated point -> get all matching stencils
      size_t k = j;
      while ((k < stencil_info_count) &&
            (stencil_recv_buffer[k].tgt_global_id == curr_orig_global_id)) {
        stencil_indices[stencil_match_count] = stencil_recv_buffer[k].idx;
        stencil_ranks[stencil_match_count] = stencil_recv_buffer[k].rank;
        w[stencil_match_count] = 1.0;
        ++k;
        ++stencil_match_count;
        ++curr_num_stencils_per_tgt;
      }

      dup_remote_points.data[dup_match_count].global_id =
        dup_recv_buffer[i].duplicated_global_id;
      dup_remote_points.data[dup_match_count].data.count = 1;
      dup_remote_points.data[dup_match_count].data.data.single.rank =
        dup_recv_buffer[i].duplicated_rank;
      dup_remote_points.data[dup_match_count].data.data.single.orig_pos =
        (uint64_t)(dup_recv_buffer[i].duplicated_orig_pos);
      num_stencils_per_tgt[dup_match_count] = curr_num_stencils_per_tgt;
      ++dup_match_count;
    }
  }

  // duplicate stencils
  yac_interp_weights_wcopy_weights(
    weights, &dup_remote_points, num_stencils_per_tgt,
    stencil_indices, stencil_ranks, w);

  // clean up
  free(w);
  free(stencil_ranks);
  free(stencil_indices);
  free(num_stencils_per_tgt);
  free(stencil_recv_buffer);
  free(dup_remote_points.data);
  free(dup_recv_buffer);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  free(int_buffer);
  free(duplicated_global_ids);
}

void yac_duplicate_stencils_f2c(
  struct yac_interp_weights * weights, struct yac_basic_grid * tgt_grid,
  yac_int * tgt_orig_global_id, size_t * tgt_duplicated_idx,
  size_t nbr_duplicated, int location) {

  yac_duplicate_stencils(
    weights, tgt_grid, tgt_orig_global_id, tgt_duplicated_idx,
    nbr_duplicated, yac_get_location(location));
}
