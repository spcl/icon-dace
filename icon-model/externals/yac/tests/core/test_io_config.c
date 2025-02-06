// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "tests.h"
#include "test_common.h"
#include "io_utils.h"

static void check_io_config(
  MPI_Comm comm, int ref_local_is_io, int * ref_io_ranks, int ref_num_io_ranks);
static void generate_node_wise_ref_data(
  MPI_Comm comm, int num_io_ranks_per_node,
  int ** ref_io_ranks, int * ref_num_io_ranks, int * ref_local_is_io);

int main(void) {

  MPI_Init(NULL, NULL);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  if (comm_size != 12) {
    PUT_ERR("ERROR: wrong number of processes");
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  { // rank list

    // clear environment
    clear_yac_io_env();

    setenv("YAC_IO_RANK_LIST", "1,3,5,7,9,11", 1);
    setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", "12", 1);

    int ref_local_is_io = comm_rank & 1;
    int ref_io_ranks[] = {1,3,5,7,9,11};
    int ref_num_io_ranks = (int)(sizeof(ref_io_ranks) / sizeof(ref_io_ranks[0]));

    check_io_config(
      MPI_COMM_WORLD, ref_local_is_io, ref_io_ranks, ref_num_io_ranks);
  }

  { // rank list (order of processes in reversed compared to MPI_COMM_WORLD)

    // clear environment
    clear_yac_io_env();

    setenv("YAC_IO_RANK_LIST", "0,1,4,7", 1);
    setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", "12", 1);

    MPI_Comm reverse_world_comm;
    MPI_Comm_split(
      MPI_COMM_WORLD, 1, comm_size - comm_rank, &reverse_world_comm);

    int ref_local_is_io = (comm_rank == 0) || (comm_rank == 1) ||
                          (comm_rank == 4) || (comm_rank == 7);
    int ref_io_ranks[] = {4,7,10,11};
    int ref_num_io_ranks = (int)(sizeof(ref_io_ranks) / sizeof(ref_io_ranks[0]));

    check_io_config(
      reverse_world_comm, ref_local_is_io, ref_io_ranks, ref_num_io_ranks);

    MPI_Comm_free(&reverse_world_comm);
  }

  { // rank list (comm contains only a subset of MPI_COMM_WORLD)

    // clear environment
    clear_yac_io_env();

    setenv("YAC_IO_RANK_LIST", "1,4,7,10", 1);
    setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", "12", 1);

    int is_in_subgroup = (comm_rank >= 4) && (comm_rank <= 7);
    MPI_Comm sub_world_comm;
    MPI_Comm_split(
      MPI_COMM_WORLD, is_in_subgroup, comm_rank, &sub_world_comm);

    int ref_local_is_io = (comm_rank == 4) || (comm_rank == 7);
    int ref_io_ranks[] = {0,3};
    int ref_num_io_ranks = (int)(sizeof(ref_io_ranks) / sizeof(ref_io_ranks[0]));

    if (is_in_subgroup)
      check_io_config(
        sub_world_comm, ref_local_is_io, ref_io_ranks, ref_num_io_ranks);

    MPI_Comm_free(&sub_world_comm);
  }

  { // rank list + maximum number of io ranks

    // clear environment
    clear_yac_io_env();

    setenv("YAC_IO_RANK_LIST", "1,3,5,7,9,11", 1);
    setenv("YAC_IO_MAX_NUM_RANKS", "3", 1);
    setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", "12", 1);

    int ref_local_is_io = (comm_rank == 1) || (comm_rank == 3) ||
                          (comm_rank == 5);
    int ref_io_ranks[] = {1,3,5};
    int ref_num_io_ranks = (int)(sizeof(ref_io_ranks) / sizeof(ref_io_ranks[0]));

    check_io_config(
      MPI_COMM_WORLD, ref_local_is_io, ref_io_ranks, ref_num_io_ranks);
  }

  { // rank list + rank exclude list

    // clear environment
    clear_yac_io_env();

    setenv("YAC_IO_RANK_LIST", "0,1,2,3,4,5,6,7,8,9,10,11", 1);
    setenv("YAC_IO_RANK_EXCLUDE_LIST", "0,2,4,6,8,10", 1);
    setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", "12", 1);

    int ref_local_is_io = comm_rank & 1;
    int ref_io_ranks[] = {1,3,5,7,9,11};
    int ref_num_io_ranks = (int)(sizeof(ref_io_ranks) / sizeof(ref_io_ranks[0]));

    check_io_config(
      MPI_COMM_WORLD, ref_local_is_io, ref_io_ranks, ref_num_io_ranks);
  }

  { // rank list + rank exclude list + sub-communicator

    // clear environment
    clear_yac_io_env();

    setenv("YAC_IO_RANK_LIST", "0,1,2,3,4,5,6,7,8,9,10,11", 1);
    setenv("YAC_IO_RANK_EXCLUDE_LIST", "0,2,4,6,8,10", 1);
    setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", "12", 1);

    int is_in_subgroup = (comm_rank >= 4) && (comm_rank <= 7);
    MPI_Comm sub_world_comm;
    MPI_Comm_split(
      MPI_COMM_WORLD, is_in_subgroup, comm_rank, &sub_world_comm);

    int ref_local_is_io = (comm_rank & 1) && is_in_subgroup;
    int ref_io_ranks[] = {1,3};
    int ref_num_io_ranks = (int)(sizeof(ref_io_ranks) / sizeof(ref_io_ranks[0]));

    if (is_in_subgroup)
      check_io_config(
        sub_world_comm, ref_local_is_io, ref_io_ranks, ref_num_io_ranks);

    MPI_Comm_free(&sub_world_comm);
  }

  { // limited number of ranks per node

    // clear environment
    clear_yac_io_env();

    setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", "1", 1);

    int ref_local_is_io;
    int * ref_io_ranks;
    int ref_num_io_ranks;

    generate_node_wise_ref_data(
      MPI_COMM_WORLD, 1, &ref_io_ranks, &ref_num_io_ranks, &ref_local_is_io);

    check_io_config(
      MPI_COMM_WORLD, ref_local_is_io, ref_io_ranks, ref_num_io_ranks);

    free(ref_io_ranks);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void check_io_config(
  MPI_Comm comm, int ref_local_is_io, int * ref_io_ranks, int ref_num_io_ranks) {

  int local_is_io;
  int * io_ranks;
  int num_io_ranks;
  yac_get_io_ranks(comm, &local_is_io, &io_ranks, &num_io_ranks);

  int comm_size;
  MPI_Comm_size(comm, &comm_size);

  int * flag = calloc((size_t)comm_size, sizeof(*flag));

  if (ref_local_is_io != local_is_io)
    PUT_ERR("wrong local_is_io");
  if (ref_num_io_ranks != num_io_ranks)
    PUT_ERR("wrong num io ranks");

  for (int i = 0; i < num_io_ranks; ++i) {
    if (flag[io_ranks[i]])
      PUT_ERR("duplicated entry in io_ranks");
    flag[io_ranks[i]] = 1;
  }

  for (int i = 0; i < ref_num_io_ranks; ++i)
    if (!flag[ref_io_ranks[i]])
      PUT_ERR("wrong io_ranks");

  free(flag);
  free(io_ranks);
}

static void generate_node_wise_ref_data(
  MPI_Comm comm, int num_io_ranks_per_node,
  int ** ref_io_ranks, int * ref_num_io_ranks, int * ref_local_is_io) {

  int comm_rank, comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  MPI_Comm split_comm;
  MPI_Comm_split_type(
    comm, MPI_COMM_TYPE_SHARED, comm_rank, MPI_INFO_NULL, &split_comm);

  int split_comm_rank;
  MPI_Comm_rank(split_comm, &split_comm_rank);

  *ref_local_is_io = split_comm_rank < num_io_ranks_per_node;

  int * is_io_rank = xmalloc((size_t)comm_size * sizeof(*is_io_rank));
  MPI_Allgather(
    ref_local_is_io, 1, MPI_INT, is_io_rank, 1, MPI_INT, comm);

  int num_io_ranks = 0;
  for (int i = 0; i < comm_size; ++i)
    if (is_io_rank[i])
      is_io_rank[num_io_ranks++] = i;

  *ref_io_ranks =
    xrealloc(is_io_rank, (size_t)num_io_ranks * sizeof(**ref_io_ranks));
  *ref_num_io_ranks = num_io_ranks;
}
