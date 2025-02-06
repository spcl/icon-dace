// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef YAC_MPI_INTERNAL_H
#define YAC_MPI_INTERNAL_H

#include <stdint.h>
#include <mpi.h>
#include "utils_core.h"
#include "yac_mpi_common.h"
#include "yac_mpi.h"

#if SIZE_MAX == UINT8_MAX
#define YAC_MPI_SIZE_T MPI_UINT8_T
#define YAC_MPI_SIZE_T_TYPE uint8_t
#elif SIZE_MAX == UINT16_MAX
#define YAC_MPI_SIZE_T MPI_UINT16_T
#define YAC_MPI_SIZE_T_TYPE uint16_t
#elif SIZE_MAX == UINT32_MAX
#define YAC_MPI_SIZE_T MPI_UINT32_T
#define YAC_MPI_SIZE_T_TYPE uint32_t
#elif SIZE_MAX == UINT64_MAX
#define YAC_MPI_SIZE_T MPI_UINT64_T
#define YAC_MPI_SIZE_T_TYPE uint64_t
#else
#error "YAC: unable to determine MPI data type for size_t"
#endif

/** \example test_group_comm.c
 * This example tests the usage of the YAC group communicator interfaces.
 */

struct yac_group_comm {
  int start;
  int size;
  MPI_Comm comm;
};

#define YAC_ALLTOALLV_P2P_DEC(N, T) \
void yac_alltoallv_ ## N ## _p2p( \
  T const * send_buffer, size_t const * sendcounts, size_t const * sdispls, \
  T * recv_buffer, size_t const * recvcounts, size_t const * rdispls, \
  MPI_Comm comm); \

YAC_ALLTOALLV_P2P_DEC(int, int)
YAC_ALLTOALLV_P2P_DEC(yac_int, yac_int)
YAC_ALLTOALLV_P2P_DEC(uint64, uint64_t)
YAC_ALLTOALLV_P2P_DEC(packed, void)
YAC_ALLTOALLV_P2P_DEC(dble, double)
YAC_ALLTOALLV_P2P_DEC(size_t, size_t)

#undef YAC_ALLTOALLV_P2P_DEC

void yac_alltoallv_p2p(
  void const * send_buffer, size_t const * sendcounts, size_t const * sdispls,
  void * recv_buffer, size_t const * recvcounts, size_t const * rdispls,
  size_t dt_size, MPI_Datatype dt, MPI_Comm comm);

void yac_alltoallv_p2p_group(
  void const * send_buffer, int const * sendcounts, int const * sdispls,
  void * recv_buffer, int const * recvcounts, int const * rdispls,
  size_t dt_size, MPI_Datatype dt, struct yac_group_comm group_comm);

void yac_allreduce_sum_dble(
  double * buffer, int count, struct yac_group_comm group_comm);

void yac_allgather_uint64(
  const uint64_t * sendbuf, uint64_t * recvbuf, int count,
  struct yac_group_comm group_comm);

void yac_bcast_group(
  void * buffer, int count, MPI_Datatype datatype, int root,
  struct yac_group_comm group_comm);

struct yac_group_comm yac_group_comm_new(MPI_Comm comm);
void yac_group_comm_delete(struct yac_group_comm group_comm);
int yac_group_comm_get_rank(struct yac_group_comm group_comm);
int yac_group_comm_get_size(struct yac_group_comm group_comm);
int yac_group_comm_get_global_rank(struct yac_group_comm group_comm);
int yac_group_comm_get_global_size(struct yac_group_comm group_comm);

void yac_group_comm_split(
  struct yac_group_comm group_comm, int split_rank,
  struct yac_group_comm * local_group_comm,
  struct yac_group_comm * remote_group_comm);

MPI_Datatype yac_get_bounding_circle_mpi_datatype(MPI_Comm comm);
MPI_Datatype yac_create_resized(
  MPI_Datatype dt, size_t new_size, MPI_Comm comm);

void yac_generate_alltoallv_args(
  int count, size_t const * sendcounts, size_t * recvcounts,
  size_t * sdispls, size_t * rdispls, MPI_Comm comm);

void yac_get_comm_buffers(
  int count, size_t ** sendcounts, size_t ** recvcounts,
  size_t ** sdispls, size_t ** rdispls, MPI_Comm comm);
void yac_free_comm_buffers(
  size_t * sendcounts, size_t * recvcounts,
  size_t * sdispls, size_t * rdispls);

#endif // YAC_MPI_INTERNAL_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
