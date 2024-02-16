/**
 * @file yac_mpi.h
 * @brief utility routines for MPI
 *
 * contains utility routines for handling MPI
 *
 * @copyright Copyright  (C)  2020 Moritz Hanke <hanke@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
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

#ifndef YAC_MPI_H
#define YAC_MPI_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdint.h>
#include <mpi.h>
#include "utils.h"

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

/**
 * check return code of MPI call and call abort function if needed
 */
#define yac_mpi_call(call, comm)                \
  do {                                          \
    int error_code = (call);                    \
    if (error_code != MPI_SUCCESS)              \
      yac_mpi_error(error_code, comm);          \
  } while(0)

/**
 * report error return of MPI call
 *
 * @param[in] error_code return code of an MPI call
 * @param[in] comm       communicator which was used for the respective MPI call
 */
void yac_mpi_error(int error_code, MPI_Comm comm);

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

void yac_yaxt_init(MPI_Comm comm);
int yac_mpi_is_initialised();
void yac_mpi_init();
void yac_mpi_cleanup();
void yac_mpi_finalize();
void yac_mpi_handshake(MPI_Comm comm, size_t n, char const** group_names,
  MPI_Comm * group_comms);

#endif // YAC_MPI_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
