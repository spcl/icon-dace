// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <mpi.h>
#include <string.h>

#include "yac_mpi_common.h"
#include "utils_common.h"

void yac_mpi_handshake(MPI_Comm comm, size_t n, char const** group_names,
  MPI_Comm * group_comms) {
  int is_intercomm;
  yac_mpi_call(MPI_Comm_test_inter(comm, &is_intercomm), comm);
  YAC_ASSERT(!is_intercomm,
             "ERROR(yac_mpi_handshake): inter-communicators are not supported");

  // STEP 1: Version exchange
  enum {MPI_HANDSHAKE_VERSION = 1};
  int version = MPI_HANDSHAKE_VERSION;
  yac_mpi_call(
    MPI_Allreduce(MPI_IN_PLACE, &version, 1, MPI_INT, MPI_MIN, comm),
    comm);
  YAC_ASSERT_F(
    version == MPI_HANDSHAKE_VERSION,
    "ERROR(yac_mpi_handshake): "
    "Version check failed. YAC only supports MPI handshake version %d",
    MPI_HANDSHAKE_VERSION);

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  for (size_t i = 0; i < n; ++i) group_comms[i] = MPI_COMM_NULL;

  while(1){
    // STEP 2: determine broadcasting rank
    size_t group_idx = SIZE_MAX;
    for (size_t i = 0; (i < n) && (group_idx == SIZE_MAX); ++i)
      if (group_comms[i] == MPI_COMM_NULL) group_idx = i;
    int broadcasting_rank = group_idx != SIZE_MAX ? rank : size;
    yac_mpi_call(
      MPI_Allreduce(MPI_IN_PLACE, &broadcasting_rank, 1, MPI_INT, MPI_MIN, comm),
      comm);
    YAC_ASSERT(broadcasting_rank >= 0 && broadcasting_rank <= size,
      "ERROR(yac_mpi_handshake): "
      "broadcasting rank cannot be negativ or greater than communicator size.");
    if(broadcasting_rank == size) break;

    // STEP 3: broadcast group name
    int groupnamelen = 0;
    if(broadcasting_rank == rank){
      size_t len = strlen(group_names[group_idx]);
      YAC_ASSERT(len <= INT_MAX,
        "ERROR(yac_mpi_handshake): group name is too long");
      groupnamelen = (int)len;
    }
    yac_mpi_call(
      MPI_Bcast(&groupnamelen, 1, MPI_INT, broadcasting_rank, comm),
      comm);
    char * groupname = xmalloc((size_t)(groupnamelen + 1) * sizeof(*groupname));
    if(broadcasting_rank == rank){
      strcpy(groupname, group_names[group_idx]);
    }
    yac_mpi_call(
      MPI_Bcast(groupname, groupnamelen, MPI_CHAR, broadcasting_rank, comm),
      comm);
    groupname[groupnamelen] = '\0';

    // STEP 4: split communicator
    group_idx = SIZE_MAX;
    for (size_t i = 0; (i < n) && (group_idx == SIZE_MAX); ++i)
      if (!strcmp(groupname, group_names[i])){
        YAC_ASSERT_F(group_comms[i] == MPI_COMM_NULL,
          "ERROR(yac_mpi_handshake): "
          "Group communicator for group '%s' was already created, "
          "but was broadcasted again.",
          groupname);
        group_idx = i;
      }
    free(groupname);
    MPI_Comm new_comm;
    yac_mpi_call(
      MPI_Comm_split(comm, (group_idx != SIZE_MAX)?0:MPI_UNDEFINED, rank, &new_comm),
      comm);
    if(group_idx != SIZE_MAX)
      group_comms[group_idx] = new_comm;
  }
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
