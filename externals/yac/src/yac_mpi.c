/**
 * @file yac_mpi.c
 *
 * copied from yaxt
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>
#include "yac_mpi.h"
#include "geometry.h"
#include "ensure_array_size.h"
#include "core/core.h"

static int mpi_initialised_by_yac = 0;
static int yaxt_initialised_by_yac = 0;
static size_t * comm_buffer = NULL;
static size_t comm_buffer_array_size = 0;
static int comm_buffer_in_use = 0;
static int init_count = 0;
static int yaxt_init_count = 0;

//taken from http://beige.ucs.indiana.edu/I590/node85.html
void yac_mpi_error(int error_code, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  char error_string[MPI_MAX_ERROR_STRING];
  int length_of_error_string, error_class;

  MPI_Error_class(error_code, &error_class);
  MPI_Error_string(error_class, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Abort(comm, error_code);
}

void yac_alltoallv_p2p(
  void const * send_buffer, size_t const * sendcounts, size_t const * sdispls,
  void * recv_buffer, size_t const * recvcounts, size_t const * rdispls,
  size_t dt_size, MPI_Datatype dt, MPI_Comm comm) {

#define USE_P2P_ALLTOALLV
#ifdef USE_P2P_ALLTOALLV
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int req_count = 0;
  for (int i = 0; i < comm_size; ++i)
    req_count += (sendcounts[i] > 0) + (recvcounts[i] > 0);
  MPI_Request * req = xmalloc((size_t)req_count * sizeof(*req));

  req_count = 0;
  for (int j = 0, lb = comm_rank, ub = comm_size; j < 2;
       ++j, lb = 0, ub = comm_rank) {
    for (int i = lb; i < ub; ++i) {
      if (sendcounts[i] > 0) {
        YAC_ASSERT_F(
          sendcounts[i] <= INT_MAX,
          "ERROR(yac_alltoallv_p2p): sendcounts[%d] = %zu exceeds INT_MAX (%d)",
          i, sendcounts[i], (int)INT_MAX)
        yac_mpi_call(
          MPI_Isend(
            (void const *)((unsigned char *)send_buffer +
                           dt_size * sdispls[i]),
            (int)(sendcounts[i]), dt, i, 0,
            comm, req + req_count), comm);
        ++req_count;
      }
      if (recvcounts[i] > 0) {
        YAC_ASSERT_F(
          recvcounts[i] <= INT_MAX,
          "ERROR(yac_alltoallv_p2p): recvcounts[%d] = %zu exceeds INT_MAX (%d)",
          i, recvcounts[i], (int)INT_MAX)
        yac_mpi_call(
          MPI_Irecv(
            (void *)((unsigned char *)recv_buffer +
                     dt_size * rdispls[i]),
            (int)(recvcounts[i]), dt, i, 0,
            comm, req + req_count), comm);
        ++req_count;
      }
    }
  }
  yac_mpi_call(MPI_Waitall(req_count, req, MPI_STATUSES_IGNORE), comm);
  free(req);
#else // USE_P2P_ALLTOALLV
  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  int * int_buffer = xmalloc(4 * comm_size * sizeof(*int_buffer));
  int * int_sendcounts = int_buffer + 0 * comm_size;
  int * int_sdispls    = int_buffer + 1 * comm_size;
  int * int_recvcounts = int_buffer + 2 * comm_size;
  int * int_rdispls    = int_buffer + 3 * comm_size;
  for (int i = 0; i < comm_size; ++i) {
    YAC_ASSERT_F(
      sendcounts[i] <= INT_MAX,
      "ERROR(yac_alltoallv_p2p): sendcounts[%d] = %zu exceeds INT_MAX (%d)",
      i, sendcounts[i], (int)INT_MAX)
    YAC_ASSERT_F(
      sdispls[i] <= INT_MAX,
      "ERROR(yac_alltoallv_p2p): sdispls[%d] = %zu exceeds INT_MAX (%d)",
      i, sdispls[i], (int)INT_MAX)
    YAC_ASSERT_F(
      recvcounts[i] <= INT_MAX,
      "ERROR(yac_alltoallv_p2p): recvcounts[%d] = %zu exceeds INT_MAX (%d)",
      i, recvcounts[i], (int)INT_MAX)
    YAC_ASSERT_F(
      rdispls[i] <= INT_MAX,
      "ERROR(yac_alltoallv_p2p): rdispls[%d] = %zu exceeds INT_MAX (%d)",
      i, rdispls[i], (int)INT_MAX)
    int_sendcounts[i] = (int)(sendcounts[i]);
    int_sdispls[i] = (int)(sdispls[i]);
    int_recvcounts[i] = (int)(recvcounts[i]);
    int_rdispls[i] = (int)(rdispls[i]);
  }
  yac_mpi_call(
    MPI_Alltoallv(send_buffer, int_sendcounts, int_sdispls, dt,
                  recv_buffer, int_recvcounts, int_rdispls, dt, comm), comm);
  free(int_buffer);
#endif // USE_P2P_ALLTOALLV
}

#define YAC_ALLTOALL_P2P_TYPE(NAME, TYPE, TYPE_SIZE, MPI_TYPE) \
  void yac_alltoallv_ ## NAME ## _p2p( \
    TYPE const * send_buffer, size_t const * sendcounts, size_t const * sdispls, \
    TYPE * recv_buffer, size_t const * recvcounts, size_t const * rdispls, \
    MPI_Comm comm) { \
    yac_alltoallv_p2p( \
      (void const *)send_buffer, sendcounts, sdispls, \
      (void *)recv_buffer, recvcounts, rdispls, \
      TYPE_SIZE, MPI_TYPE, comm); \
  }

YAC_ALLTOALL_P2P_TYPE(int, int, sizeof(int), MPI_INT)
YAC_ALLTOALL_P2P_TYPE(yac_int, yac_int, sizeof(yac_int), yac_int_dt)
YAC_ALLTOALL_P2P_TYPE(uint64, uint64_t, sizeof(uint64_t), MPI_UINT64_T)
YAC_ALLTOALL_P2P_TYPE(packed, void, 1, MPI_PACKED)
YAC_ALLTOALL_P2P_TYPE(dble, double, sizeof(double), MPI_DOUBLE)
YAC_ALLTOALL_P2P_TYPE(size_t, size_t, sizeof(size_t), YAC_MPI_SIZE_T)

void yac_alltoallv_p2p_group(
  void const * send_buffer, int const * sendcounts, int const * sdispls,
  void * recv_buffer, int const * recvcounts, int const * rdispls,
  size_t dt_size, MPI_Datatype dt, struct yac_group_comm group_comm) {

  MPI_Comm comm = group_comm.comm;
  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  int rank = comm_rank - group_comm.start;

  int req_count = 0;
  for (int i = 0; i < group_comm.size; ++i)
    req_count += (sendcounts[i] > 0) + (recvcounts[i] > 0);
  MPI_Request * req = xmalloc((size_t)req_count * sizeof(*req));

  req_count = 0;
  for (int j = 0, lb = rank, ub = group_comm.size; j < 2;
       ++j, lb = 0, ub = rank) {
    for (int i = lb; i < ub; ++i) {
      if (sendcounts[i] > 0) {

        yac_mpi_call(
          MPI_Isend(
            (void const *)((unsigned char *)send_buffer +
                           dt_size * (size_t)(sdispls[i])),
            sendcounts[i], dt, i + group_comm.start, 0,
            comm, req + req_count), comm);
        ++req_count;
      }
      if (recvcounts[i] > 0) {
        yac_mpi_call(
          MPI_Irecv(
            (void *)((unsigned char *)recv_buffer +
                     dt_size * (size_t)(rdispls[i])),
            recvcounts[i], dt, i + group_comm.start, 0,
            comm, req + req_count), comm);
        ++req_count;
      }
    }
  }
  yac_mpi_call(MPI_Waitall(req_count, req, MPI_STATUSES_IGNORE), comm);
  free(req);
}

static int nearest_power_of_two(int x) {
  int power = 1;

  while(power < x) power *= 2;

  return power / 2;
}

// based on https://doi.org/10.1016/j.parco.2017.08.004
void yac_allreduce_sum_dble(
  double * buffer, int count, struct yac_group_comm group_comm) {

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(group_comm.comm, &comm_rank), group_comm.comm);

  int rank = comm_rank - group_comm.start;
  int pof2 = nearest_power_of_two(group_comm.size);
  int rem = group_comm.size - pof2;
  int my_rank;
  double * recv_buffer = xmalloc((size_t)count * sizeof(*recv_buffer));

  if (rank < 2 * rem) {

    if (rank & 1) {
      yac_mpi_call(
        MPI_Recv(
          (void*)recv_buffer, count, MPI_DOUBLE, rank - 1 + group_comm.start, 0,
          group_comm.comm, MPI_STATUS_IGNORE), group_comm.comm);
      for (int i = 0; i < count; ++i) buffer[i] += recv_buffer[i];
      my_rank = rank / 2;
    } else {
      yac_mpi_call(
        MPI_Send(
          (void const *)buffer, count, MPI_DOUBLE, rank + 1 + group_comm.start,
          0, group_comm.comm), group_comm.comm);
      my_rank = -1;
    }
  } else {
    my_rank = rank - rem;
  }
  if (my_rank != -1) {
    int mask = 1;
    while (mask < pof2) {
      int newdst = my_rank ^ mask;
      int dst;
      if (newdst < rem) dst = newdst * 2 + 1;
      else              dst = newdst + rem;
      yac_mpi_call(
        MPI_Sendrecv(
          (void const*)buffer, count, MPI_DOUBLE, dst + group_comm.start, 0,
          (void*)recv_buffer, count, MPI_DOUBLE, dst + group_comm.start, 0,
          group_comm.comm, MPI_STATUS_IGNORE),
        group_comm.comm);
      for (int i = 0; i < count; ++i) buffer[i] += recv_buffer[i];
      mask <<= 1;
    }
  }
  free(recv_buffer);
  if (rank < 2 * rem) {
    if (rank & 1) {
      yac_mpi_call(
        MPI_Send(
          (void const*)buffer, count, MPI_DOUBLE, rank - 1 + group_comm.start,
          0, group_comm.comm), group_comm.comm);
    } else {
      yac_mpi_call(
        MPI_Recv(
          (void*)buffer, count, MPI_DOUBLE, rank + 1 + group_comm.start, 0,
          group_comm.comm, MPI_STATUS_IGNORE), group_comm.comm);
    }
  }
}

static int log2_(int x) {
  if (x <= 1) return 0;
  int l2 = 0;
  while (x >>= 1) ++l2;
  return l2;
}

// based on https://doi.org/10.1109/71.642949
void yac_allgather_uint64(
  const uint64_t * sendbuf, uint64_t * recvbuf, int count,
  struct yac_group_comm group_comm) {

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(group_comm.comm, &comm_rank), group_comm.comm);
  int rank = comm_rank - group_comm.start;

  uint64_t * temp = xmalloc((size_t)group_comm.size * (size_t)count * sizeof(*temp));

  int lg2 = log2_(group_comm.size);
  memcpy(temp, sendbuf, (size_t)count * sizeof(*temp));
  int nblk = 1;
  int curr_len = count;

  for (int r = 0; r < lg2; ++r) {
    int dst = (rank - nblk + group_comm.size) % group_comm.size;
    int src = (rank + nblk) % group_comm.size;
    yac_mpi_call(
      MPI_Sendrecv(
        (void const*)temp, curr_len, MPI_UINT64_T, dst + group_comm.start, 0,
        (void *)(temp + (size_t)curr_len), curr_len, MPI_UINT64_T,
        src + group_comm.start, 0, group_comm.comm, MPI_STATUS_IGNORE),
      group_comm.comm);
    nblk *= 2;
    curr_len *= 2;
  }
  int rest = count * group_comm.size - curr_len;
  int dst = (rank - nblk + group_comm.size) % group_comm.size;
  int src = (rank + nblk) % group_comm.size;
  yac_mpi_call(
    MPI_Sendrecv(
      (void const*)temp, rest, MPI_UINT64_T, dst + group_comm.start, 0,
      (void*)(temp + (size_t)curr_len), rest, MPI_UINT64_T,
      src + group_comm.start, 0, group_comm.comm, MPI_STATUS_IGNORE),
    group_comm.comm);
  memcpy(recvbuf + (size_t)count * (size_t)rank,
         temp, (size_t)count * (size_t)(group_comm.size - rank) * sizeof(*temp));
  memcpy(recvbuf, temp + (size_t)count * (size_t)(group_comm.size - rank),
         (size_t)count * (size_t)rank * sizeof(*temp));

  free(temp);
}

void yac_bcast_group(
  void * buffer, int count, MPI_Datatype datatype, int root,
  struct yac_group_comm group_comm) {

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(group_comm.comm, &comm_rank), group_comm.comm);
  int rank = comm_rank - group_comm.start;

  // if root is not part of the group
  if ((root < group_comm.start) ||
      (root >= group_comm.start + group_comm.size)) {

    if (comm_rank == root) {
      yac_mpi_call(
        MPI_Send(
          (void const*)buffer, count, datatype, group_comm.start, 0,
          group_comm.comm), group_comm.comm);
      return;
    } else if (comm_rank == group_comm.start) {
      yac_mpi_call(
        MPI_Recv(
          buffer, count, datatype, root, 0, group_comm.comm,
          MPI_STATUS_IGNORE), group_comm.comm);
    }
    root = 0;
  } else {
    root -= group_comm.start;
  }

  // if not root, receive data
  if (rank != root) {

    int temp_rank = (group_comm.size + rank - root) % group_comm.size;
    int bit = 1;

    while (bit <= temp_rank) bit <<= 1;
    bit >>= 1;

    int src_rank =
      (((temp_rank ^ bit) + root) % group_comm.size) + group_comm.start;

    yac_mpi_call(
      MPI_Recv(buffer, count, datatype, src_rank, 0, group_comm.comm,
               MPI_STATUS_IGNORE), group_comm.comm);
  }

  // relative rank in respect to root
  int temp_rank = (group_comm.size + rank - root) % group_comm.size;
  int bit = 1, send_rank;

  while(bit <= temp_rank) bit <<= 1;

  while ((send_rank = temp_rank | bit) < group_comm.size) {

    bit <<= 1;

    send_rank = ((send_rank + root) % group_comm.size) + group_comm.start;

    yac_mpi_call(
      MPI_Send(
        (void const*)buffer, count, datatype, send_rank, 0, group_comm.comm),
      group_comm.comm);
  }
}

struct yac_group_comm yac_group_comm_new(MPI_Comm comm) {

  struct yac_group_comm group_comm;
  group_comm.start = 0;
  yac_mpi_call(MPI_Comm_size(comm, &(group_comm.size)), comm);
  yac_mpi_call(MPI_Comm_dup(comm, &(group_comm.comm)), comm);

  return group_comm;
}

void yac_group_comm_delete(struct yac_group_comm group_comm) {

  yac_mpi_call(MPI_Comm_free(&(group_comm.comm)), group_comm.comm);
}

int yac_group_comm_get_rank(struct yac_group_comm group_comm) {
  return yac_group_comm_get_global_rank(group_comm) - group_comm.start;
}

int yac_group_comm_get_size(struct yac_group_comm group_comm) {
  return group_comm.size;
}

int yac_group_comm_get_global_rank(struct yac_group_comm group_comm) {
  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(group_comm.comm, &comm_rank), group_comm.comm);
  return comm_rank;
}

int yac_group_comm_get_global_size(struct yac_group_comm group_comm) {
  int comm_size;
  yac_mpi_call(MPI_Comm_size(group_comm.comm, &comm_size), group_comm.comm);
  return comm_size;
}

void yac_group_comm_split(
  struct yac_group_comm group_comm, int split_rank,
  struct yac_group_comm * local_group_comm,
  struct yac_group_comm * remote_group_comm) {

  int comm_rank;
  MPI_Comm comm = group_comm.comm;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  YAC_ASSERT(
    (split_rank >= 0) && (split_rank < group_comm.size),
    "ERROR(yac_group_comm_split): invalid split rank")

  int start[2] = {group_comm.start, group_comm.start + split_rank};
  int size[2] = {split_rank, group_comm.size - split_rank};
  int local_idx = (comm_rank - group_comm.start) >= split_rank;

  local_group_comm->start = start[local_idx];
  local_group_comm->size = size[local_idx];
  local_group_comm->comm = comm;
  remote_group_comm->start = start[local_idx^1];
  remote_group_comm->size = size[local_idx^1];
  remote_group_comm->comm = comm;
}

MPI_Datatype yac_get_bounding_circle_mpi_datatype(MPI_Comm comm) {

  struct bounding_circle dummy;
  MPI_Datatype bnd_circle_dt;
  int array_of_blocklengths[] = {3, 1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.base_vector[0]) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.inc_angle.sin) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.inc_angle.cos) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  yac_mpi_call(
    MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements,
                           array_of_types, &bnd_circle_dt), comm);
  return yac_create_resized(bnd_circle_dt, sizeof(dummy), comm);
}

MPI_Datatype yac_create_resized(
  MPI_Datatype dt, size_t new_size, MPI_Comm comm) {

  MPI_Datatype resized_dt;

#define OPENMPI_WORKAROUND
#ifdef OPENMPI_WORKAROUND
  MPI_Aint lb, extent;
  MPI_Type_get_extent(dt, &lb, &extent);
  yac_mpi_call(
    MPI_Type_create_resized(dt, lb, (MPI_Aint)new_size, &resized_dt), comm);
#else
  yac_mpi_call(
    MPI_Type_create_resized(dt, 0, (MPI_Aint)new_size, &resized_dt), comm);
#endif
#undef OPENMPI_WORKAROUND
  yac_mpi_call(MPI_Type_free(&dt), comm);
  yac_mpi_call(MPI_Type_commit(&resized_dt), comm);
  return resized_dt;
}

void yac_generate_alltoallv_args(
  int count, size_t const * sendcounts, size_t * recvcounts,
  size_t * sdispls, size_t * rdispls, MPI_Comm comm) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // exchange the number of requested points
  yac_mpi_call(
    MPI_Alltoall(
      (void const*)sendcounts, count, YAC_MPI_SIZE_T,
      (void*)recvcounts, count, YAC_MPI_SIZE_T, comm), comm);

  // sdispls are offset by one position, this is intentional because this is
  // usefull for packing of data
  sdispls[0] = 0;
  for (size_t i = 0, saccu = 0, raccu = 0; i < count * comm_size; ++i) {
    sdispls[i+1] = saccu;
    rdispls[i] = raccu;
    saccu += sendcounts[i];
    raccu += recvcounts[i];
  }
}

void yac_get_comm_buffers(
  int count, size_t ** sendcounts, size_t ** recvcounts,
  size_t ** sdispls, size_t ** rdispls, MPI_Comm comm) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * comm_buffer_;
  if (!comm_buffer_in_use) {
    ENSURE_ARRAY_SIZE(
      comm_buffer, comm_buffer_array_size,
      4 * (size_t)count * (size_t)comm_size + 1);
    comm_buffer_ = comm_buffer;
    comm_buffer_in_use = 1;
  } else {
    comm_buffer_ =
      xmalloc(
        (4 * (size_t)count * (size_t)comm_size + 1) * sizeof(*comm_buffer_));
  }

  size_t offset = (size_t)count * (size_t)comm_size;
  *sendcounts = comm_buffer_ + 0 * offset;
  *recvcounts = comm_buffer_ + 1 * offset;
  *rdispls =    comm_buffer_ + 2 * offset;
  *sdispls =    comm_buffer_ + 3 * offset; // sdispls is bigger by one element,
                                       // which is usefull when packing data
                                       // for alltoallv operation
  memset(
    comm_buffer_, 0, (size_t)count * (size_t)comm_size * sizeof(*comm_buffer_));
}

void yac_free_comm_buffers(
  size_t * sendcounts, size_t * recvcounts,
  size_t * sdispls, size_t * rdispls) {

  if (sendcounts != comm_buffer) free(sendcounts);
  else comm_buffer_in_use = 0;
}

int yac_mpi_is_initialised() {

  int mpi_initialized;
  MPI_Initialized(&mpi_initialized);

  return mpi_initialized;
}

void yac_yaxt_init(MPI_Comm comm) {

  YAC_ASSERT(
    yac_mpi_is_initialised(),
    "ERROR(yac_yaxt_init): MPI has not yet been initialised");

  YAC_ASSERT(
    (yaxt_init_count == 0) || !yaxt_initialised_by_yac,
    "ERROR(yac_yaxt_init): YAXT was initialised by YAC. \n"
    "In case there are multiple instances of YAC in parallel, the user has "
    "to initialise YAXT such that it is available on all processes that "
    "use YAC.")

  if ((yaxt_init_count == 0) && (!xt_initialized() || xt_finalized())) {
    xt_initialize(comm);
    yaxt_initialised_by_yac = 1;
  }
  yaxt_init_count++;
}

static void yac_yaxt_cleanup() {

  if ((yaxt_init_count == 1) && yaxt_initialised_by_yac) {
    xt_finalize();
    yaxt_initialised_by_yac = 0;
  }
  yaxt_init_count--;
}

#define XSTR(s) STR(s)
#define STR(s) #s
#define CHECK_YAC_MPI_SIZE_T \
  YAC_ASSERT_F( \
    sizeof(size_t) == sizeof(YAC_MPI_SIZE_T_TYPE), \
    "ERROR(CHECK_YAC_MPI_SIZE_T): " \
    "could not determine MPI data type for size_t " \
    "(sizeof(size_t): %zu; sizeof(%s): %zu)", \
    sizeof(size_t), XSTR(YAC_MPI_SIZE_T_TYPE), \
    sizeof(YAC_MPI_SIZE_T_TYPE))

void yac_mpi_init() {

  CHECK_YAC_MPI_SIZE_T

  if ((init_count == 0) && (!yac_mpi_is_initialised())) {
    MPI_Init(NULL, NULL);
    mpi_initialised_by_yac = 1;
  }

  init_count++;
}

void yac_mpi_cleanup() {

  if (init_count == 1) {
    YAC_ASSERT(
      !comm_buffer_in_use, "ERROR(yac_mpi_finalize): comm_buffer still in use")
    free(comm_buffer);
    comm_buffer = NULL;
    comm_buffer_array_size = 0;
  }
  init_count--;
  yac_yaxt_cleanup();
}

void yac_mpi_finalize() {

  yac_mpi_cleanup();
  if (mpi_initialised_by_yac)
    yac_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);
}

void yac_mpi_handshake(MPI_Comm comm, size_t n, char const** group_names,
  MPI_Comm * group_comms){
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
