// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include "tests.h"
#include "yac_mpi_internal.h"

int main (void) {

  MPI_Init(NULL, NULL);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  if (comm_size != 9) {

    PUT_ERR("ERROR: wrong number of processes\n");
    return TEST_EXIT_CODE;
  }

  struct yac_group_comm world_group_comm = yac_group_comm_new(MPI_COMM_WORLD);

  if (yac_group_comm_get_global_rank(world_group_comm) != comm_rank)
    PUT_ERR("ERROR in yac_group_comm_get_global_rank");
  if (yac_group_comm_get_global_size(world_group_comm) != comm_size)
    PUT_ERR("ERROR in yac_group_comm_get_global_size");
  if (yac_group_comm_get_rank(world_group_comm) != comm_rank)
    PUT_ERR("ERROR in yac_group_comm_get_rank");
  if (yac_group_comm_get_size(world_group_comm) != comm_size)
    PUT_ERR("ERROR in yac_group_comm_get_size");

  struct yac_group_comm local_group_comm, remote_group_comm;
  int split_rank = 3;
  yac_group_comm_split(
    world_group_comm, split_rank, &local_group_comm, &remote_group_comm);

  int group_idx = comm_rank >= split_rank;

  int ref_group_size[2] = {3, 6};
  int group_rank_offset[2] = {0, 3};

  if (yac_group_comm_get_global_rank(local_group_comm) != comm_rank)
    PUT_ERR("ERROR in yac_group_comm_get_global_rank");
  if (yac_group_comm_get_global_size(local_group_comm) != comm_size)
    PUT_ERR("ERROR in yac_group_comm_get_global_size");
  if (yac_group_comm_get_global_rank(remote_group_comm) != comm_rank)
    PUT_ERR("ERROR in yac_group_comm_get_global_rank");
  if (yac_group_comm_get_global_size(remote_group_comm) != comm_size)
    PUT_ERR("ERROR in yac_group_comm_get_global_size");

  if (yac_group_comm_get_rank(local_group_comm) !=
      (comm_rank - group_rank_offset[group_idx]))
    PUT_ERR("ERROR in yac_group_comm_get_rank");
  if (yac_group_comm_get_size(local_group_comm) != ref_group_size[group_idx])
    PUT_ERR("ERROR in yac_group_comm_get_size");
  if (yac_group_comm_get_rank(remote_group_comm) !=
      (comm_rank - group_rank_offset[group_idx^1]))
    PUT_ERR("ERROR in yac_group_comm_get_rank");
  if (yac_group_comm_get_size(remote_group_comm) != ref_group_size[group_idx^1])
    PUT_ERR("ERROR in yac_group_comm_get_size");

  double dble_proc_data[9][3] = {{1,2,3},
                                 {4,5,6},
                                 {7,8,9},
                                 {10,11,12},
                                 {13,14,15},
                                 {16,17,18},
                                 {19,20,21},
                                 {22,23,24},
                                 {25,26,27}};
  double ref_dble_sum[2][3] = {{1+4+7,2+5+8,3+6+9},
                               {10+13+16+19+22+25,
                                11+14+17+20+23+26,
                                12+15+18+21+24+27}};

  double local_dble_proc_data[3];
  for (int i = 0; i < 3; ++i)
    local_dble_proc_data[i] = dble_proc_data[comm_rank][i];
  yac_allreduce_sum_dble(local_dble_proc_data, 3, local_group_comm);
  for (int i = 0; i < 3; ++i)
    if (fabs(local_dble_proc_data[i] - ref_dble_sum[group_idx][i]) > 1e-9)
      PUT_ERR("ERROR in yac_allreduce_sum_dble");

  uint64_t uint64_proc_data[9][2] = {{1,2},
                                     {3,4},
                                     {5,6},
                                     {7,8},
                                     {9,10},
                                     {11,12},
                                     {13,14},
                                     {15,16},
                                     {17,18}};
  uint64_t uint64_recv_buffer[12];
  uint64_t ref_uint64_allgather[2][12] = {{1,2,3,4,5,6},
                                          {7,8,9,10,11,12,13,14,15,16,17,18}};
  size_t ref_uint64_allgather_size[2] = {2*3, 2*6};

  yac_allgather_uint64(
    uint64_proc_data[comm_rank], uint64_recv_buffer, 2, local_group_comm);
  for (size_t i = 0; i < ref_uint64_allgather_size[group_idx]; ++i)
    if (uint64_recv_buffer[i] != ref_uint64_allgather[group_idx][i])
      PUT_ERR("ERROR in yac_allgather_uint64");

  for (int i = 0; i < comm_size; ++i) {

    // if current rank is root
    if (i == comm_rank) {

      // bcast to remote group
      yac_bcast_group(
        dble_proc_data[comm_rank], 3, MPI_DOUBLE, i, remote_group_comm);

      // bcast to local group
      yac_bcast_group(
        dble_proc_data[comm_rank], 3, MPI_DOUBLE, i, local_group_comm);

    } else {

      // bcast within local group
      double bcast_data[3];
      yac_bcast_group(bcast_data, 3, MPI_DOUBLE, i, local_group_comm);
      for (int j = 0; j < 3; ++j)
        if (fabs(bcast_data[j] - dble_proc_data[i][j]) > 1e-9)
          PUT_ERR("ERROR in yac_bcast_group");
    }
  }

  yac_group_comm_delete(world_group_comm);

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
