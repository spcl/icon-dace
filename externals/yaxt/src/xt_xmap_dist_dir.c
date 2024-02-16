/**
 * @file xt_xmap_dist_dir.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#include <config.h>
#endif

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <mpi.h>

#include "xt/xt_idxlist.h"
#include "xt/xt_idxlist_collection.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_xmap.h"
#include "xt/xt_xmap_dist_dir.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_xmap_intersection.h"
#include "xt_config_internal.h"
#include "xt_idxlist_internal.h"
#include "instr.h"
#include "xt_xmap_dist_dir_common.h"

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ == 11 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static inline Xt_int
get_dist_dir_global_interval_size(Xt_idxlist src, Xt_idxlist dst,
                                  bool *stripify, MPI_Comm comm, int comm_size,
                                  Xt_config config)
{
  unsigned long long local_vals[2], global_sums[2];

  unsigned num_indices_src = (unsigned)xt_idxlist_get_num_indices(src);
  local_vals[0] = num_indices_src;
  local_vals[1] = (num_indices_src >= (unsigned)config->idxv_cnv_size)
    || (xt_idxlist_get_num_indices(dst) >= config->idxv_cnv_size);

  xt_mpi_call(MPI_Allreduce(local_vals, global_sums, 2,
                            MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm), comm);

  *stripify = global_sums[1] > 0;
  return (Xt_int)(MAX(((global_sums[0] + (unsigned)comm_size - 1)
                      / (unsigned)comm_size), 1) * (unsigned)comm_size);
}

static inline Xt_int get_min_idxlist_index(Xt_idxlist a, Xt_idxlist b) {

  int num_a = xt_idxlist_get_num_indices(a),
    num_b = xt_idxlist_get_num_indices(b);
  Xt_int min_index_a = num_a ? xt_idxlist_get_min_index(a) : XT_INT_MAX,
    min_index_b = num_b ? xt_idxlist_get_min_index(b) : XT_INT_MAX,
    min_index = (Xt_int)MIN(min_index_a, min_index_b);
  return min_index;
}

static inline Xt_int get_max_idxlist_index(Xt_idxlist a, Xt_idxlist b) {

  int num_a = xt_idxlist_get_num_indices(a),
    num_b = xt_idxlist_get_num_indices(b);
  Xt_int max_index_a = num_a ? xt_idxlist_get_max_index(a) : XT_INT_MIN,
    max_index_b = num_b ? xt_idxlist_get_max_index(b) : XT_INT_MIN,
    max_index = MAX(max_index_a, max_index_b);
  return max_index;
}

static struct bucket_params
get_bucket_params(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                  bool *stripify, MPI_Comm comm, int comm_size,
                  Xt_config config)
{
  /* global_interval is a multiple of comm_size and has a size of at least
     comm_size */
  Xt_int global_interval
    = get_dist_dir_global_interval_size(src_idxlist, dst_idxlist, stripify,
                                        comm, comm_size, config);
  Xt_int local_interval = (Xt_int)(global_interval / comm_size);
  Xt_int local_index_range_lbound
    = get_min_idxlist_index(src_idxlist, dst_idxlist);
  Xt_int local_index_range_ubound
    = get_max_idxlist_index(src_idxlist, dst_idxlist);
  return (struct bucket_params){
    .global_interval = global_interval,
    .local_interval = local_interval,
    .local_index_range_lbound = local_index_range_lbound,
    .local_index_range_ubound = local_index_range_ubound,
  };
}

enum {
  SEND_SIZE_SRC = 0,
  SEND_SIZE_DST = 1,
  SEND_NUM_SRC = 2,
  SEND_NUM_DST = 3,
  SEND_SIZE_ASIZE,
};

static inline void
rank_no_send(size_t rank, int (*restrict send_size)[SEND_SIZE_ASIZE])
{
  send_size[rank][SEND_SIZE_SRC] = 0;
  send_size[rank][SEND_NUM_SRC] = 0;
  send_size[rank][SEND_SIZE_DST] = 0;
  send_size[rank][SEND_NUM_DST] = 0;
}

static int
compute_and_pack_bucket_intersections(struct bucket_params *bucket_params,
                                      Xt_idxlist src_idxlist,
                                      Xt_idxlist dst_idxlist,
                                      int (*send_size)[SEND_SIZE_ASIZE],
                                      void **send_buffer_,
                                      MPI_Comm comm, int comm_size) {


  int num_msg = 0;
  Xt_int local_index_range_lbound = bucket_params->local_index_range_lbound;
  Xt_int local_index_range_ubound = bucket_params->local_index_range_ubound;
  if (local_index_range_lbound <= local_index_range_ubound) {
    size_t send_buffer_size = 0;
    struct Xt_stripe *stripes = NULL;
    size_t stripes_array_size = 0;
    Xt_int global_interval = bucket_params->global_interval;
    Xt_int local_interval = bucket_params->local_interval;
    size_t first_overlapping_bucket = 0;
    /* is it impossible for early buckets to overlap our lists? */
    if (local_index_range_lbound >= 0
        && (local_index_range_ubound < global_interval)) {
      first_overlapping_bucket
        = (size_t)(local_index_range_lbound / local_interval);
      for (size_t i = 0; i < first_overlapping_bucket; ++i)
        rank_no_send(i, send_size);
    }
    /* is it impossible for later ranks to overlap our lists? */
    size_t start_of_non_overlapping_bucket_suffix
      = (size_t)(((long long)local_index_range_ubound + local_interval - 1)
                 / local_interval) + 1;
    if (local_index_range_lbound < 0
        || start_of_non_overlapping_bucket_suffix > (size_t)comm_size)
      start_of_non_overlapping_bucket_suffix = (size_t)comm_size;
    size_t max_num_intersect
      = start_of_non_overlapping_bucket_suffix - first_overlapping_bucket;
    struct {
      Xt_idxlist list;
      int rank;
    } *restrict sends_dst
      = xmalloc(2 * max_num_intersect * sizeof(*sends_dst)),
      *restrict sends_src = sends_dst + max_num_intersect;
    size_t i = first_overlapping_bucket;
    size_t num_src_msg = 0, num_dst_msg = 0;
    for (; i < (size_t)start_of_non_overlapping_bucket_suffix; ++i) {
      Xt_idxlist bucket
        = xt_xmap_dist_dir_get_bucket(bucket_params,
                                      &stripes, &stripes_array_size, (int)i);
      if (bucket) {
        Xt_idxlist send4src = xt_idxlist_get_intersection(src_idxlist, bucket);
        if (xt_idxlist_get_num_indices(send4src) > 0) {
          sends_src[num_src_msg].list = send4src;
          sends_src[num_src_msg].rank = (int)i;
          send_buffer_size += xt_idxlist_get_pack_size(send4src, comm);
          send_size[i][SEND_NUM_SRC] = 1;
          ++num_src_msg;
        } else {
          send_size[i][SEND_SIZE_SRC] = 0;
          send_size[i][SEND_NUM_SRC] = 0;
          xt_idxlist_delete(send4src);
        }

        Xt_idxlist send4dst = xt_idxlist_get_intersection(bucket, dst_idxlist);
        xt_idxlist_delete(bucket);
        if (xt_idxlist_get_num_indices(send4dst) > 0) {
          sends_dst[num_dst_msg].list = send4dst;
          sends_dst[num_dst_msg].rank = (int)i;
          send_buffer_size += xt_idxlist_get_pack_size(send4dst, comm);
          send_size[i][SEND_NUM_DST] = 1;
          ++num_dst_msg;
        } else {
          send_size[i][SEND_SIZE_DST] = 0;
          send_size[i][SEND_NUM_DST] = 0;
          xt_idxlist_delete(send4dst);
        }
      } else
        rank_no_send(i, send_size);
    }
    for (; i < (size_t)comm_size; ++i)
      rank_no_send(i, send_size);

    unsigned char *send_buffer
      = *send_buffer_ = xrealloc(stripes, send_buffer_size);
    size_t ofs = 0;
    for (i = 0; i < num_src_msg; ++i) {
      int position = 0;
      xt_idxlist_pack(sends_src[i].list, send_buffer+ofs,
                      (int)(send_buffer_size-ofs), &position, comm);
      send_size[sends_src[i].rank][SEND_SIZE_SRC] = position;
      ofs += (size_t)position;
      xt_idxlist_delete(sends_src[i].list);
    }

    for (i = 0; i < num_dst_msg; ++i) {
      int position = 0;
      xt_idxlist_pack(sends_dst[i].list, send_buffer+ofs,
                      (int)(send_buffer_size-ofs), &position, comm);
      send_size[sends_dst[i].rank][SEND_SIZE_DST] = position;
      ofs += (size_t)position;
      xt_idxlist_delete(sends_dst[i].list);
    }
    free(sends_dst);
    num_msg = (int)(num_src_msg + num_dst_msg);
  } else
    memset(send_size, 0, (size_t)comm_size * sizeof (*send_size));
  return num_msg;
}

static void
recv_and_unpack_intersection(struct dist_dir *dist_dir, int recv_size,
                             int recv_count, void * recv_buffer, int tag,
                             MPI_Comm comm) {

  // initialize distributed directories
  int total_recv_size = 0;

  for (int i = 0; i < recv_count; ++i)
  {
    MPI_Status status;

    xt_mpi_call(MPI_Recv(recv_buffer, recv_size, MPI_PACKED, MPI_ANY_SOURCE,
                         tag, comm, &status), comm);

    int received_count;
    xt_mpi_call(MPI_Get_count(&status, MPI_PACKED, &received_count), comm);

    int position = 0;

    dist_dir->entries[i].rank = status.MPI_SOURCE;
    dist_dir->entries[i].list =
      xt_idxlist_unpack(recv_buffer, received_count, &position, comm);

    total_recv_size += received_count;
  }

  if (total_recv_size != recv_size)
    Xt_abort(comm, "ERROR: recv_intersections received wrong number of bytes",
             __FILE__, __LINE__);
  dist_dir->num_entries = recv_count;
}

static void send_intersections(void *send_buffer,
                               const int (*send_size)[SEND_SIZE_ASIZE],
                               MPI_Request *dir_init_send_requests,
                               int tag_offset, MPI_Comm comm, int comm_size) {
  int src_tag = tag_offset + xt_mpi_tag_xmap_dist_dir_src_send;
  struct Xt_xmdd_txstat txstat
    = xt_xmap_dist_dir_send_intersections(send_buffer,
                                          SEND_SIZE_ASIZE, SEND_SIZE_SRC,
                                          src_tag, comm, comm_size,
                                          dir_init_send_requests,
                                          send_size);
  int dst_tag = tag_offset + xt_mpi_tag_xmap_dist_dir_dst_send;
  xt_xmap_dist_dir_send_intersections((unsigned char *)send_buffer
                                      + txstat.bytes,
                                      SEND_SIZE_ASIZE, SEND_SIZE_DST,
                                      dst_tag, comm, comm_size,
                                      dir_init_send_requests + txstat.num_msg,
                                      send_size);
}

static void
recv_and_unpack_intersections(int recv_size[SEND_SIZE_ASIZE],
                              struct dist_dir **src_dist_dir,
                              struct dist_dir **dst_dist_dir,
                              int tag_offset, MPI_Comm comm) {

  *src_dist_dir = xmalloc(sizeof (struct dist_dir)
                          + (sizeof (struct Xt_com_list)
                             * (size_t)recv_size[SEND_NUM_SRC]));
  *dst_dist_dir = xmalloc(sizeof (struct dist_dir)
                          + (sizeof (struct Xt_com_list)
                             * (size_t)recv_size[SEND_NUM_DST]));

  void * recv_buffer = xmalloc((size_t)MAX(recv_size[SEND_SIZE_SRC],
                                           recv_size[SEND_SIZE_DST]));

  recv_and_unpack_intersection(*src_dist_dir, recv_size[SEND_SIZE_SRC],
                               recv_size[SEND_NUM_SRC], recv_buffer,
                               tag_offset + xt_mpi_tag_xmap_dist_dir_src_send,
                               comm);
  recv_and_unpack_intersection(*dst_dist_dir, recv_size[SEND_SIZE_DST],
                               recv_size[SEND_NUM_DST], recv_buffer,
                               tag_offset + xt_mpi_tag_xmap_dist_dir_dst_send,
                               comm);

  free(recv_buffer);
}


static size_t
buf_size_from_intersections(size_t num_intersections,
                            const struct isect *restrict src_dst_intersections,
                            MPI_Comm comm, int comm_size,
                            int (*restrict send_size)[SEND_SIZE_ASIZE])
{
  size_t total_send_size = 0;
  for (int i = 0; i < comm_size; ++i)
    (void)(send_size[i][SEND_SIZE_SRC] = 0),
      (void)(send_size[i][SEND_SIZE_DST] = 0),
      (void)(send_size[i][SEND_NUM_SRC] = 0),
      (void)(send_size[i][SEND_NUM_DST] = 0);

  int rank_pack_size;
  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &rank_pack_size), comm);

  for (size_t i = 0; i < num_intersections; ++i)
  {
    int msg_size = rank_pack_size
      + (int)xt_idxlist_get_pack_size(src_dst_intersections[i].idxlist,
                                      comm);
    size_t src_rank
      = (size_t)src_dst_intersections[i].rank[xt_xmdd_direction_src],
      dst_rank = (size_t)src_dst_intersections[i].rank[xt_xmdd_direction_dst];
    /* send_size[i][SEND_SIZE_(SRC|DST)] are set when actually
     * packing, because that provides a potentially tighter bound,
     * see xt_xmap_dist_dir_pack_intersections */
    ++(send_size[src_rank][SEND_NUM_SRC]);
    ++(send_size[dst_rank][SEND_NUM_DST]);
    total_send_size += 2*(size_t)msg_size;
  }
  assert(total_send_size <= INT_MAX);
  return total_send_size;
}


static int
pack_src_dst_dist_dirs(size_t num_intersections,
                       struct isect *restrict src_dst_intersections,
                       int (*send_size)[SEND_SIZE_ASIZE],
                       void **send_buffer_,
                       MPI_Comm comm, int comm_size) {

  size_t total_send_size
    = buf_size_from_intersections(num_intersections,
                                  src_dst_intersections,
                                  comm, comm_size, send_size);

  unsigned char *send_buffer = (*send_buffer_)
    = xmalloc((size_t)total_send_size);
  size_t ofs = 0;
  if (num_intersections > 1)
    qsort(src_dst_intersections, num_intersections,
          sizeof (src_dst_intersections[0]), xt_xmdd_cmp_isect_src_rank);
  size_t num_send_indices_requests
    = xt_xmap_dist_dir_pack_intersections(
      xt_xmdd_direction_src, num_intersections, src_dst_intersections, false,
      SEND_SIZE_ASIZE, SEND_SIZE_SRC, send_size,
      send_buffer, total_send_size, &ofs, comm);

  if (num_intersections > 1)
    qsort(src_dst_intersections, num_intersections,
          sizeof (src_dst_intersections[0]), xt_xmdd_cmp_isect_dst_rank);
  num_send_indices_requests
    += xt_xmap_dist_dir_pack_intersections(
      xt_xmdd_direction_dst, num_intersections, src_dst_intersections, true,
      SEND_SIZE_ASIZE, SEND_SIZE_DST, send_size,
      send_buffer, total_send_size, &ofs, comm);
  assert(num_send_indices_requests <= INT_MAX);
  return (int)num_send_indices_requests;
}

/**
 * @brief wrapper for MPI_Reduce_scatter_block if available or
 * MPI_Reduce_scatter if not
 * @param num_sizes number of size entries to reduce over and to be
 * received via @a recv_size
 * @param recv_size array to hold result of reduction
 * @param send_size sizes to sum over, array size must correspond to
 * corresponding size of comm times @a num_sizes
 * @param comm MPI communicator to use
 */
static void
xt_xmap_dist_dir_reduce_scatter_sizes(int num_sizes,
                                      int recv_size[num_sizes],
                                      int (*send_size)[num_sizes],
                                      MPI_Comm comm) {

#if MPI_VERSION > 2 || ( MPI_VERSION == 2 && MPI_SUBVERSION >= 2)
  xt_mpi_call(MPI_Reduce_scatter_block((int *)send_size, (int *)recv_size,
                                       num_sizes, MPI_INT, MPI_SUM,
                                       comm), comm);
#else
  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int *recv_count = xmalloc((size_t)comm_size * sizeof(*recv_count));
  for (int i = 0; i < comm_size; ++i) recv_count[i] = num_sizes;

  xt_mpi_call(MPI_Reduce_scatter(send_size, recv_size, recv_count, MPI_INT,
                                 MPI_SUM, comm), comm);

  free(recv_count);
#endif
}

static void generate_distributed_directories(struct dist_dir **src_dist_dir,
                                             struct dist_dir **dst_dist_dir,
                                             bool *stripify,
                                             Xt_idxlist src_idxlist,
                                             Xt_idxlist dst_idxlist,
                                             int tag_offset,
                                             MPI_Comm comm, int comm_size,
                                             Xt_config config) {

  struct bucket_params bucket_params
    = get_bucket_params(src_idxlist, dst_idxlist, stripify, comm, comm_size,
                        config);

  void *send_buffer = NULL;

  int (*send_size)[SEND_SIZE_ASIZE]
    = xmalloc((size_t)comm_size * sizeof(*send_size));

  int num_msg
    = compute_and_pack_bucket_intersections(&bucket_params,
                                            src_idxlist, dst_idxlist,
                                            send_size, &send_buffer,
                                            comm, comm_size);

  int recv_size[SEND_SIZE_ASIZE]; // for src and dst

  /* get packed intersection sizes to be sent from other ranks */
  xt_xmap_dist_dir_reduce_scatter_sizes(SEND_SIZE_ASIZE, recv_size, send_size, comm);

  MPI_Request *dir_init_send_requests
    = xmalloc((size_t)num_msg * sizeof(*dir_init_send_requests));
  send_intersections(send_buffer, (const int (*)[SEND_SIZE_ASIZE])send_size,
                     dir_init_send_requests, tag_offset, comm, comm_size);

  free(send_size);

  recv_and_unpack_intersections(recv_size, src_dist_dir, dst_dist_dir,
                                tag_offset, comm);

  // wait for the sends to be completed
  xt_mpi_call(MPI_Waitall(num_msg, dir_init_send_requests,
                          MPI_STATUSES_IGNORE), comm);
  free(dir_init_send_requests);
  free(send_buffer);
}

static void
recv_and_unpack_dist_dir_result(struct dist_dir *dist_dir, int recv_size,
                                void *restrict recv_buffer, int tag,
                                MPI_Comm comm)
{

  // initiate distributed directories
  int num_entries = 0;
  while (recv_size > 0) {

    MPI_Status status;

    xt_mpi_call(MPI_Recv(recv_buffer, recv_size, MPI_PACKED,
                         MPI_ANY_SOURCE, tag, comm, &status), comm);

    int received_count;
    xt_mpi_call(MPI_Get_count(&status, MPI_PACKED, &received_count), comm);

    recv_size -= received_count;

    int position = 0;

    while (received_count > position) {

      xt_mpi_call(MPI_Unpack(recv_buffer, received_count, &position,
                             &dist_dir->entries[num_entries].rank,
                             1, MPI_INT, comm), comm);

      dist_dir->entries[num_entries].list =
        xt_idxlist_unpack(recv_buffer, received_count, &position, comm);

      ++num_entries;
    }
  }
  qsort(dist_dir->entries, (size_t)num_entries, sizeof(*dist_dir->entries),
        xt_com_list_rank_cmp);

  if (0 != recv_size)
    Xt_abort(comm, "ERROR: recv_and_unpack_dist_dir_result"
             " received wrong number of bytes", __FILE__, __LINE__);

  dist_dir->num_entries = num_entries;

}


static void
recv_and_unpack_dist_dir_results(int recv_size[SEND_SIZE_ASIZE],
                                 struct dist_dir **src_intersections,
                                 struct dist_dir **dst_intersections,
                                 int *num_send_indices_requests,
                                 MPI_Request *send_indices_requests,
                                 int tag_offset,
                                 MPI_Comm comm) {

  struct dist_dir *src_dist_dir_results
    = xmalloc(sizeof (struct dist_dir)
              + (sizeof (struct Xt_com_list)
                 * (size_t)recv_size[SEND_NUM_SRC])),
    *dst_dist_dir_results
    = xmalloc(sizeof (struct dist_dir)
              + (sizeof (struct Xt_com_list)
                 * (size_t)recv_size[SEND_NUM_DST]));

  void *recv_buffer = xmalloc((size_t)MAX(recv_size[SEND_SIZE_SRC],
                                          recv_size[SEND_SIZE_DST]));

  recv_and_unpack_dist_dir_result(src_dist_dir_results,
                                  recv_size[SEND_SIZE_SRC],
                                  recv_buffer, tag_offset
                                  + xt_mpi_tag_xmap_dist_dir_src_send, comm);
  assert(src_dist_dir_results->num_entries
         == recv_size[SEND_NUM_SRC]);

  enum { ops_completed_auto_size = 16 };
  int ops_completed_auto[ops_completed_auto_size];
  int *ops_completed
    = *num_send_indices_requests > ops_completed_auto_size
    ? xmalloc((size_t)*num_send_indices_requests * sizeof (*ops_completed))
    : ops_completed_auto;
  bool all_sends_done
    = xt_mpi_test_some(num_send_indices_requests, send_indices_requests,
                       ops_completed, comm);

  recv_and_unpack_dist_dir_result(dst_dist_dir_results,
                                  recv_size[SEND_SIZE_DST],
                                  recv_buffer, tag_offset
                                  + xt_mpi_tag_xmap_dist_dir_dst_send, comm);
  assert(dst_dist_dir_results->num_entries
         == recv_size[SEND_NUM_DST]);

  if (!all_sends_done)
    all_sends_done
      = xt_mpi_test_some(num_send_indices_requests, send_indices_requests,
                         ops_completed, comm);
  free(recv_buffer);

  xt_xmap_dist_dir_same_rank_merge(&src_dist_dir_results);
  *src_intersections = src_dist_dir_results;

  if (!all_sends_done)
    all_sends_done
      = xt_mpi_test_some(num_send_indices_requests, send_indices_requests,
                         ops_completed, comm);

  xt_xmap_dist_dir_same_rank_merge(&dst_dist_dir_results);
  *dst_intersections = dst_dist_dir_results;
  if (ops_completed != ops_completed_auto) free(ops_completed);
}

static void exchange_idxlists(struct dist_dir **src_intersections,
                              struct dist_dir **dst_intersections,
                              bool *stripify,
                              Xt_idxlist src_idxlist,
                              Xt_idxlist dst_idxlist,
                              int tag_offset,
                              MPI_Comm comm,
                              Xt_config config) {

  int comm_size;

  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  struct dist_dir *src_dist_dir, *dst_dist_dir;

  generate_distributed_directories(&src_dist_dir, &dst_dist_dir, stripify,
                                   src_idxlist, dst_idxlist,
                                   tag_offset, comm, comm_size, config);

  void * send_buffer;

  int recv_size[SEND_SIZE_ASIZE], (*send_size)[SEND_SIZE_ASIZE]
    = xmalloc((size_t)comm_size * sizeof(*send_size));

  /* match the source and destination entries in the local distributed
   * directories... */
  struct isect *src_dst_intersections;
  size_t num_intersections
    = xt_xmap_dist_dir_match_src_dst(src_dist_dir, dst_dist_dir,
                                     &src_dst_intersections);
  xt_xmdd_free_dist_dir(src_dist_dir);
  xt_xmdd_free_dist_dir(dst_dist_dir);
  /* ... and pack the results into a sendable format */
  int num_send_indices_requests
    = pack_src_dst_dist_dirs(num_intersections, src_dst_intersections,
                             send_size, &send_buffer, comm, comm_size);
  free(src_dst_intersections);

  // get the data size the local process will receive from other processes
  xt_xmap_dist_dir_reduce_scatter_sizes(SEND_SIZE_ASIZE, recv_size,
                                        send_size, comm);

  MPI_Request *send_indices_requests
    = xmalloc((size_t)num_send_indices_requests
              * sizeof(*send_indices_requests));

  send_intersections(send_buffer, (const int (*)[SEND_SIZE_ASIZE])send_size,
                     send_indices_requests, tag_offset, comm, comm_size);

  recv_and_unpack_dist_dir_results(recv_size,
                                   src_intersections, dst_intersections,
                                   &num_send_indices_requests,
                                   send_indices_requests,
                                   tag_offset, comm);

  xt_mpi_call(MPI_Waitall(num_send_indices_requests, send_indices_requests,
                          MPI_STATUSES_IGNORE), comm);

  free(send_buffer);
  free(send_size);
  free(send_indices_requests);
}

Xt_xmap
xt_xmap_dist_dir_intracomm_custom_new(
  Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
  MPI_Comm comm, Xt_config config)
{
  INSTR_DEF(this_instr,"xt_xmap_all2all_new")
  INSTR_START(this_instr);

  // ensure that yaxt is initialized
  assert(xt_initialized());

  int tag_offset;
  MPI_Comm newcomm = xt_mpi_comm_smart_dup(comm, &tag_offset);

  struct dist_dir *src_intersections, *dst_intersections;

  bool stripify;
  exchange_idxlists(&src_intersections, &dst_intersections, &stripify,
                    src_idxlist, dst_idxlist, tag_offset, newcomm, config);

  Xt_xmap (*xmap_new)(int num_src_intersections,
                      const struct Xt_com_list *src_com,
                      int num_dst_intersections,
                      const struct Xt_com_list *dst_com,
                      Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                      MPI_Comm comm)
    = stripify ? xt_xmap_intersection_ext_new : xt_xmap_intersection_new;

  Xt_xmap xmap
    = xmap_new(src_intersections->num_entries, src_intersections->entries,
               dst_intersections->num_entries, dst_intersections->entries,
               src_idxlist, dst_idxlist, newcomm);

  xt_mpi_comm_smart_dedup(&newcomm, tag_offset);

  xt_xmdd_free_dist_dir(src_intersections);
  xt_xmdd_free_dist_dir(dst_intersections);
  INSTR_STOP(this_instr);
  return xmap;
}

Xt_xmap
xt_xmap_dist_dir_intracomm_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                               MPI_Comm comm)
{
  return xt_xmap_dist_dir_intracomm_custom_new(src_idxlist, dst_idxlist, comm,
                                               &xt_default_config);
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
