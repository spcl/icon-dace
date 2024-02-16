/**
 * @file xt_xmap_dist_dir_intercomm.c
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
#include "xt/xt_xmap_dist_dir_intercomm.h"
#include "xt/xt_mpi.h"
#include "xt_arithmetic_util.h"
#include "xt_idxstripes_internal.h"
#include "xt_mpi_internal.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_xmap_intersection.h"
#include "xt_idxlist_internal.h"
#include "xt_xmap_dist_dir_common.h"
#include "xt_config_internal.h"
#include "instr.h"

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ == 11 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static inline void
get_dist_dir_global_interval_size(Xt_idxlist src, Xt_idxlist dst,
                                  bool *stripify, Xt_int interval_size[2],
                                  MPI_Comm intra_comm, MPI_Comm inter_comm,
                                  int comm_size, int remote_size,
                                  int comm_rank, int tag_offset_inter,
                                  Xt_config config)
{
  /* global_sums[0] and [1] refer to the local and remote group of
   * intercommunicator inter_comm */
  unsigned long long local_vals[2], global_sums[2][2];

  unsigned num_indices_src = (unsigned)xt_idxlist_get_num_indices(src);
  local_vals[0] = num_indices_src;
  local_vals[1] = (num_indices_src >= (unsigned)config->idxv_cnv_size)
    || (xt_idxlist_get_num_indices(dst) >= config->idxv_cnv_size);

  xt_mpi_call(MPI_Allreduce(local_vals, global_sums[0], 2,
                            MPI_UNSIGNED_LONG_LONG, MPI_SUM, intra_comm),
              intra_comm);
  /* instead of sendrecv one might use hand-programmed multi-casts
   * sending to each rank in a range from the remote group and
   * receiving from the first rank in that group,
   * the better choice probably depends on the asymmetry of the group
   * sizes, i.e. use bcast from a very small to a very large group
   * and few sends from a large to a small group */
  if (comm_rank == 0) {
    int tag = tag_offset_inter + xt_mpi_tag_xmap_dist_dir_src_send;
    xt_mpi_call(MPI_Sendrecv(global_sums[0], 2, MPI_UNSIGNED_LONG_LONG, 0, tag,
                             global_sums[1], 2, MPI_UNSIGNED_LONG_LONG, 0, tag,
                             inter_comm, MPI_STATUS_IGNORE), inter_comm);
  }
  xt_mpi_call(MPI_Bcast(global_sums[1], 2, MPI_UNSIGNED_LONG_LONG,
                        0, intra_comm), intra_comm);
  *stripify = (global_sums[0][1] > 0 || global_sums[1][1] > 0);
  interval_size[0]
    = (Xt_int)(((global_sums[0][0] + (unsigned)comm_size - 1)
                / (unsigned)comm_size) * (unsigned)comm_size);
  interval_size[1]
    = (Xt_int)(((global_sums[1][0] + (unsigned)remote_size - 1)
                / (unsigned)remote_size) * (unsigned)remote_size);
}

enum {
  SEND_SIZE = 0,
  SEND_NUM = 1,
  SEND_SIZE_ASIZE,
};

static inline void
rank_no_send(size_t rank, int (*restrict send_size)[SEND_SIZE_ASIZE])
{
  send_size[rank][SEND_SIZE] = 0;
  send_size[rank][SEND_NUM] = 0;
}

static size_t
compute_and_pack_bucket_intersections(struct bucket_params *bucket_params,
                                      Xt_idxlist idxlist,
                                      int (*send_size)[SEND_SIZE_ASIZE],
                                      void **send_buffer_,
                                      MPI_Comm comm, int comm_size)
{
  size_t num_msg = 0;
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
    } *restrict sends
      = xmalloc(max_num_intersect * sizeof(*sends));
    size_t i = first_overlapping_bucket;
    for (; i < (size_t)start_of_non_overlapping_bucket_suffix; ++i) {
      Xt_idxlist bucket
        = xt_xmap_dist_dir_get_bucket(bucket_params,
                                      &stripes, &stripes_array_size, (int)i);
      if (bucket) {
        Xt_idxlist isect2send = xt_idxlist_get_intersection(idxlist, bucket);
        xt_idxlist_delete(bucket);
        if (xt_idxlist_get_num_indices(isect2send) > 0) {
          sends[num_msg].list = isect2send;
          sends[num_msg].rank = (int)i;
          send_buffer_size += xt_idxlist_get_pack_size(isect2send, comm);
          /* send_size[i][SEND_SIZE] is set below after the actual
           * pack, because MPI_Pack_size only gives an upper bound,
           * not the actually needed size */
          send_size[i][SEND_NUM] = 1;
          ++num_msg;
        } else {
          rank_no_send(i, send_size);
          xt_idxlist_delete(isect2send);
        }
      } else
        rank_no_send(i, send_size);

    }
    for (; i < (size_t)comm_size; ++i)
      rank_no_send(i, send_size);

    unsigned char *send_buffer
      = *send_buffer_ = xrealloc(stripes, send_buffer_size);
    size_t ofs = 0;
    for (i = 0; i < num_msg; ++i) {
      int position = 0;
      xt_idxlist_pack(sends[i].list, send_buffer + ofs,
                      (int)(send_buffer_size-ofs), &position, comm);
      send_size[sends[i].rank][SEND_SIZE] = position;
      ofs += (size_t)position;
      xt_idxlist_delete(sends[i].list);
    }

    free(sends);
  } else
    memset(send_size, 0, (size_t)comm_size * sizeof (*send_size));
  return num_msg;
}

static inline Xt_int
get_min_idxlist_index(Xt_idxlist l)
{
  int num_idx = xt_idxlist_get_num_indices(l);
  Xt_int min_index = num_idx ? xt_idxlist_get_min_index(l) : XT_INT_MAX;
  return min_index;
}

static inline Xt_int
get_max_idxlist_index(Xt_idxlist l)
{
  int num_idx = xt_idxlist_get_num_indices(l);
  Xt_int max_index = num_idx ? xt_idxlist_get_max_index(l) : XT_INT_MIN;
  return max_index;
}

static struct bucket_params
get_bucket_params(Xt_idxlist idxlist,
                  Xt_int global_interval, int comm_size)
{
  /* guard vs. comm_size being larger than number of indices */
  Xt_int local_interval = MAX((Xt_int)1, (Xt_int)(global_interval / comm_size));
  Xt_int local_index_range_lbound = get_min_idxlist_index(idxlist);
  Xt_int local_index_range_ubound = get_max_idxlist_index(idxlist);
  return (struct bucket_params){
    .global_interval = global_interval,
    .local_interval = local_interval,
    .local_index_range_lbound = local_index_range_lbound,
    .local_index_range_ubound = local_index_range_ubound,
  };
}

static void
compress_sizes(int (*restrict sizes)[SEND_SIZE_ASIZE], int comm_size,
               struct Xt_xmdd_txstat *tx_stat, int *counts)
{
  size_t tx_num = 0, size_sum = 0;
  for (size_t i = 0; i < (size_t)comm_size; ++i)
    if (sizes[i][SEND_SIZE]) {
      int tx_size = sizes[i][SEND_SIZE];
      size_sum += (size_t)tx_size;
      sizes[tx_num][SEND_SIZE] = tx_size;
      if (counts) counts[tx_num] = sizes[i][SEND_NUM];
      sizes[tx_num][SEND_NUM] = (int)i;
      ++tx_num;
    }
  *tx_stat = (struct Xt_xmdd_txstat){ .bytes = size_sum, .num_msg = tx_num };
}

static void
create_intersections(struct Xt_xmdd_txstat tx_stat[2],
                     int recv_size[][SEND_SIZE_ASIZE],
                     void **send_buffer, int send_size[][SEND_SIZE_ASIZE],
                     Xt_idxlist idxlist, Xt_int interval_size,
                     MPI_Comm comm, int comm_size)
{
  struct bucket_params bucket_params
    = get_bucket_params(idxlist, interval_size, comm_size);
  *send_buffer = NULL;
#ifndef NDEBUG
  size_t num_msg =
#endif
    compute_and_pack_bucket_intersections(
      &bucket_params, idxlist, send_size, send_buffer, comm, comm_size);
  xt_mpi_call(MPI_Alltoall(send_size, SEND_SIZE_ASIZE, MPI_INT,
                           recv_size, SEND_SIZE_ASIZE, MPI_INT, comm), comm);
  compress_sizes(recv_size, comm_size, tx_stat + 0, NULL);
  compress_sizes(send_size, comm_size, tx_stat + 1, NULL);
  assert(num_msg == tx_stat[1].num_msg);
}

typedef int (*tx_fp)(void *, int, MPI_Datatype, int, int,
                     MPI_Comm, MPI_Request *);
static void
tx_intersections(size_t num_msg,
                 const int (*sizes)[SEND_SIZE_ASIZE],
                 unsigned char *buffer, MPI_Request *requests,
                 int tag, MPI_Comm comm, tx_fp tx_op)
{
  size_t ofs = 0;
  for (size_t i = 0; i < num_msg; ++i)
  {
    int rank = sizes[i][SEND_NUM], count = sizes[i][SEND_SIZE];
    xt_mpi_call(tx_op(buffer + ofs,
                      count, MPI_PACKED, rank, tag, comm, requests + i), comm);
    ofs += (size_t)count;
  }
}

static void
irecv_intersections(size_t num_msg,
                    const int (*recv_size)[SEND_SIZE_ASIZE],
                    void *recv_buffer, MPI_Request *requests,
                    int tag, MPI_Comm comm)
{
  tx_intersections(num_msg, recv_size, recv_buffer, requests, tag, comm,
                   (tx_fp)MPI_Irecv);
}

static void
isend_intersections(size_t num_msg,
                    const int (*send_size)[SEND_SIZE_ASIZE],
                    void *send_buffer, MPI_Request *requests,
                    int tag, MPI_Comm comm)
{
  tx_intersections(num_msg, send_size, send_buffer, requests, tag, comm,
                   (tx_fp)MPI_Isend);
}


static void
unpack_dist_dir(struct Xt_xmdd_txstat tx_stat,
                const int (*sizes)[SEND_SIZE_ASIZE],
                void *buffer,
                struct dist_dir **dist_dir,
                MPI_Comm comm)
{
  size_t num_msg = tx_stat.num_msg, buf_size = tx_stat.bytes;
  struct dist_dir *restrict dist_dir_
    = *dist_dir = xmalloc(sizeof (*dist_dir_)
                          + sizeof (*dist_dir_->entries) * num_msg);
  dist_dir_->num_entries = (int)num_msg;
  int position = 0;
  for (size_t i = 0; i < num_msg; ++i)
  {
    int rank = sizes[i][SEND_NUM];
    dist_dir_->entries[i].rank = rank;
    dist_dir_->entries[i].list
      = xt_idxlist_unpack(buffer, (int)buf_size, &position, comm);
  }
}

static void
generate_distributed_directories(struct dist_dir **src_dist_dir,
                                 struct dist_dir **dst_dist_dir,
                                 bool *stripify,
                                 Xt_idxlist src_idxlist,
                                 Xt_idxlist dst_idxlist,
                                 int tag_offset_inter, int tag_offset_intra,
                                 MPI_Comm inter_comm, MPI_Comm intra_comm,
                                 int remote_size, int comm_size,
                                 int comm_rank, Xt_config config) {

  /* interval_size[0] and interval_size[1] are the global interval
   * size for the local and remote group */
  Xt_int interval_size[2];
  get_dist_dir_global_interval_size(src_idxlist, dst_idxlist,
                                    stripify, interval_size,
                                    intra_comm, inter_comm,
                                    comm_size, remote_size,
                                    comm_rank, tag_offset_inter, config);
  void *send_buffer_local, *send_buffer_remote;
  int (*send_size_local)[SEND_SIZE_ASIZE]
    = xmalloc(((size_t)comm_size + (size_t)remote_size)
              * 2 * sizeof(*send_size_local)),
    (*send_size_remote)[SEND_SIZE_ASIZE] = send_size_local + comm_size,
    (*recv_size_local)[SEND_SIZE_ASIZE] = send_size_remote + remote_size,
    (*recv_size_remote)[SEND_SIZE_ASIZE] = recv_size_local + comm_size;
  struct Xt_xmdd_txstat tx_stat_local[2], tx_stat_remote[2];
  create_intersections(tx_stat_local, recv_size_local, &send_buffer_local,
                       send_size_local, src_idxlist, interval_size[0],
                       intra_comm, comm_size);
  create_intersections(tx_stat_remote, recv_size_remote, &send_buffer_remote,
                       send_size_remote, dst_idxlist, interval_size[1],
                       inter_comm, remote_size);

  size_t num_req = tx_stat_local[0].num_msg + tx_stat_remote[0].num_msg
    + tx_stat_local[1].num_msg + tx_stat_remote[1].num_msg;
  MPI_Request *dir_init_requests
    = xmalloc(num_req * sizeof(*dir_init_requests)
              + tx_stat_local[0].bytes + tx_stat_remote[0].bytes);
  void *recv_buffer_local = dir_init_requests + num_req,
    *recv_buffer_remote = ((unsigned char *)recv_buffer_local
                           + tx_stat_local[0].bytes);
  int tag_intra = tag_offset_intra + xt_mpi_tag_xmap_dist_dir_src_send;
  size_t req_ofs = tx_stat_local[0].num_msg;
  irecv_intersections(tx_stat_local[0].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])recv_size_local,
                      recv_buffer_local, dir_init_requests,
                      tag_intra, intra_comm);
  int tag_inter = tag_offset_inter + xt_mpi_tag_xmap_dist_dir_src_send;
  irecv_intersections(tx_stat_remote[0].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])recv_size_remote,
                      recv_buffer_remote, dir_init_requests + req_ofs,
                      tag_inter, inter_comm);
  req_ofs += tx_stat_remote[0].num_msg;
  isend_intersections(tx_stat_local[1].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])send_size_local,
                      send_buffer_local, dir_init_requests + req_ofs,
                      tag_intra, intra_comm);
  req_ofs += tx_stat_local[1].num_msg;
  isend_intersections(tx_stat_remote[1].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])send_size_remote,
                      send_buffer_remote, dir_init_requests + req_ofs,
                      tag_inter, inter_comm);
  // wait for data transfers to complete
  xt_mpi_call(MPI_Waitall((int)num_req, dir_init_requests,
                          MPI_STATUSES_IGNORE), inter_comm);
  free(send_buffer_local);
  free(send_buffer_remote);
  unpack_dist_dir(tx_stat_local[0],
                  (const int (*)[SEND_SIZE_ASIZE])recv_size_local,
                  recv_buffer_local, src_dist_dir, intra_comm);
  unpack_dist_dir(tx_stat_remote[0],
                  (const int (*)[SEND_SIZE_ASIZE])recv_size_remote,
                  recv_buffer_remote, dst_dist_dir, inter_comm);
  free(send_size_local);
  free(dir_init_requests);
}


static size_t
send_size_from_intersections(size_t num_intersections,
                             const struct isect *restrict src_dst_intersections,
                             enum xt_xmdd_direction target,
                             MPI_Comm comm, int comm_size,
                             int (*restrict send_size_target)[SEND_SIZE_ASIZE])
{
  size_t total_send_size = 0;
  for (int i = 0; i < comm_size; ++i)
    (void)(send_size_target[i][SEND_SIZE] = 0),
      (void)(send_size_target[i][SEND_NUM] = 0);

  int rank_pack_size;
  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &rank_pack_size), comm);

  for (size_t i = 0; i < num_intersections; ++i)
  {
    int msg_size = rank_pack_size
      + (int)xt_idxlist_get_pack_size(src_dst_intersections[i].idxlist, comm);
    size_t target_rank = (size_t)src_dst_intersections[i].rank[target];
    /* send_size_target[target_rank][SEND_SIZE] += msg_size; */
    ++(send_size_target[target_rank][SEND_NUM]);
    total_send_size += (size_t)msg_size;
  }
  assert(total_send_size <= INT_MAX);
  return total_send_size;
}


static size_t
pack_dist_dirs(size_t num_intersections,
               struct isect *restrict src_dst_intersections,
               int (*send_size)[SEND_SIZE_ASIZE],
               void **send_buffer_, enum xt_xmdd_direction target,
               bool isect_idxlist_delete, MPI_Comm comm, int comm_size) {

  size_t total_send_size
    = send_size_from_intersections(num_intersections,
                                   src_dst_intersections,
                                   target,
                                   comm, comm_size, send_size);

  unsigned char *send_buffer = (*send_buffer_)
    = xmalloc((size_t)total_send_size);
  qsort(src_dst_intersections, num_intersections,
        sizeof (src_dst_intersections[0]),
        target == xt_xmdd_direction_src
        ? xt_xmdd_cmp_isect_src_rank : xt_xmdd_cmp_isect_dst_rank);
  size_t ofs = 0;
  size_t num_requests
    = xt_xmap_dist_dir_pack_intersections(
      target, num_intersections, src_dst_intersections,
      isect_idxlist_delete,
      SEND_SIZE_ASIZE, SEND_SIZE, send_size,
      send_buffer, total_send_size, &ofs, comm);
  return num_requests;
}

static void
unpack_dist_dir_results(struct Xt_xmdd_txstat tx_stat,
                        struct dist_dir **dist_dir,
                        void *restrict recv_buffer,
                        int *restrict entry_counts,
                        MPI_Comm comm)
{
  size_t num_msg = tx_stat.num_msg;
  int buf_size = (int)tx_stat.bytes;
  int position = 0;
  size_t num_entries_sent = 0;
  for (size_t i = 0; i < num_msg; ++i)
    num_entries_sent += (size_t)entry_counts[i];
  *dist_dir = xmalloc(sizeof (struct dist_dir)
                      + (sizeof (struct Xt_com_list) * num_entries_sent));
  (*dist_dir)->num_entries = (int)num_entries_sent;
  struct Xt_com_list *restrict entries = (*dist_dir)->entries;
  size_t num_entries = 0;
  for (size_t i = 0; i < num_msg; ++i) {
    size_t num_entries_from_rank = (size_t)entry_counts[i];
    for (size_t j = 0; j < num_entries_from_rank; ++j) {
      xt_mpi_call(MPI_Unpack(recv_buffer, buf_size, &position,
                             &entries[num_entries].rank,
                             1, MPI_INT, comm), comm);
      entries[num_entries].list =
        xt_idxlist_unpack(recv_buffer, buf_size, &position, comm);
      ++num_entries;
    }
  }
  assert(num_entries == num_entries_sent);
  qsort(entries, num_entries_sent, sizeof(*entries), xt_com_list_rank_cmp);
  xt_xmap_dist_dir_same_rank_merge(dist_dir);
}


static void
exchange_idxlists(struct dist_dir **src_intersections,
                  struct dist_dir **dst_intersections,
                  bool *stripify,
                  Xt_idxlist src_idxlist,
                  Xt_idxlist dst_idxlist,
                  int tag_offset_inter, int tag_offset_intra,
                  MPI_Comm inter_comm, MPI_Comm intra_comm,
                  Xt_config config) {

  int comm_size, remote_size, comm_rank;
  xt_mpi_call(MPI_Comm_size(inter_comm, &comm_size), inter_comm);
  xt_mpi_call(MPI_Comm_rank(inter_comm, &comm_rank), inter_comm);
  xt_mpi_call(MPI_Comm_remote_size(inter_comm, &remote_size), inter_comm);

  struct dist_dir *src_dist_dir, *dst_dist_dir;

  generate_distributed_directories(&src_dist_dir, &dst_dist_dir, stripify,
                                   src_idxlist, dst_idxlist,
                                   tag_offset_inter, tag_offset_intra,
                                   inter_comm, intra_comm,
                                   remote_size, comm_size,
                                   comm_rank, config);


  int (*send_size_local)[SEND_SIZE_ASIZE]
    = xmalloc(((size_t)comm_size + (size_t)remote_size)
              * 2U * sizeof(*send_size_local)),
    (*recv_size_local)[SEND_SIZE_ASIZE] = send_size_local + comm_size,
    (*send_size_remote)[SEND_SIZE_ASIZE] = recv_size_local + comm_size,
    (*recv_size_remote)[SEND_SIZE_ASIZE] = send_size_remote + remote_size;

  /* match the source and destination entries in the local distributed
   * directories... */
  struct isect *src_dst_intersections;
  size_t num_intersections
    = xt_xmap_dist_dir_match_src_dst(src_dist_dir, dst_dist_dir,
                                     &src_dst_intersections);
  xt_xmdd_free_dist_dir(src_dist_dir);
  xt_xmdd_free_dist_dir(dst_dist_dir);
  /* ... and pack the results into a sendable format */
  void *send_buffer_local, *send_buffer_remote;
  size_t num_send_requests_local
    = pack_dist_dirs(num_intersections, src_dst_intersections,
                     send_size_local, &send_buffer_local,
                     xt_xmdd_direction_src, false, intra_comm, comm_size),
    num_send_requests_remote
    = pack_dist_dirs(num_intersections, src_dst_intersections,
                     send_size_remote, &send_buffer_remote,
                     xt_xmdd_direction_dst, true, inter_comm, remote_size);
  free(src_dst_intersections);

  // get the data size the local process will receive from other processes
  xt_mpi_call(MPI_Alltoall(send_size_local, SEND_SIZE_ASIZE, MPI_INT,
                           recv_size_local, SEND_SIZE_ASIZE, MPI_INT,
                           intra_comm), intra_comm);
  xt_mpi_call(MPI_Alltoall(send_size_remote, SEND_SIZE_ASIZE, MPI_INT,
                           recv_size_remote, SEND_SIZE_ASIZE, MPI_INT,
                           inter_comm), inter_comm);

  struct Xt_xmdd_txstat tx_stat_local[2], tx_stat_remote[2];
  int *isect_counts_recv_local
    = xmalloc(((size_t)comm_size + (size_t)remote_size) * sizeof (int)),
    *isect_counts_recv_remote = isect_counts_recv_local + comm_size;
  compress_sizes(send_size_local, comm_size, tx_stat_local+1, NULL);
  compress_sizes(recv_size_local, comm_size, tx_stat_local+0,
                 isect_counts_recv_local);
  compress_sizes(send_size_remote, remote_size, tx_stat_remote+1, NULL);
  compress_sizes(recv_size_remote, remote_size, tx_stat_remote+0,
                 isect_counts_recv_remote);
  assert(tx_stat_local[1].num_msg == num_send_requests_local
         && tx_stat_remote[1].num_msg == num_send_requests_remote);
  size_t num_requests
    = num_send_requests_local + num_send_requests_remote
    + tx_stat_local[0].num_msg + tx_stat_remote[0].num_msg;
  assert(num_requests <= INT_MAX);
  MPI_Request *requests
    = xmalloc(num_requests * sizeof(*requests)
              + tx_stat_local[0].bytes + tx_stat_remote[0].bytes);
  void *recv_buf_local = requests + num_requests,
    *recv_buf_remote = (unsigned char *)recv_buf_local + tx_stat_local[0].bytes;
  size_t req_ofs = tx_stat_local[0].num_msg;
  int tag_intra = tag_offset_intra + xt_mpi_tag_xmap_dist_dir_src_send;
  irecv_intersections(tx_stat_local[0].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])recv_size_local,
                      recv_buf_local, requests, tag_intra, intra_comm);
  int tag_inter = tag_offset_inter + xt_mpi_tag_xmap_dist_dir_src_send;
  irecv_intersections(tx_stat_remote[0].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])recv_size_remote,
                      recv_buf_remote, requests+req_ofs, tag_inter, inter_comm);
  req_ofs += tx_stat_remote[0].num_msg;
  isend_intersections(tx_stat_local[1].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])send_size_local,
                      send_buffer_local, requests+req_ofs, tag_intra,
                      intra_comm);
  req_ofs += tx_stat_local[1].num_msg;
  isend_intersections(tx_stat_remote[1].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])send_size_remote,
                      send_buffer_remote, requests+req_ofs, tag_inter,
                      inter_comm);
  xt_mpi_call(MPI_Waitall((int)num_requests, requests, MPI_STATUSES_IGNORE),
              inter_comm);
  free(send_buffer_local);
  free(send_buffer_remote);
  free(send_size_local);

  unpack_dist_dir_results(tx_stat_local[0], src_intersections, recv_buf_local,
                          isect_counts_recv_local, intra_comm);
  unpack_dist_dir_results(tx_stat_remote[0], dst_intersections, recv_buf_remote,
                          isect_counts_recv_remote, inter_comm);
  free(requests);
  free(isect_counts_recv_local);
}



Xt_xmap
xt_xmap_dist_dir_intercomm_custom_new(Xt_idxlist src_idxlist,
                                      Xt_idxlist dst_idxlist,
                                      MPI_Comm inter_comm_,
                                      MPI_Comm intra_comm_,
                                      Xt_config config)
{
  INSTR_DEF(this_instr,"xt_xmap_dist_dir_intercomm_new")
  INSTR_START(this_instr);

  // ensure that yaxt is initialized
  assert(xt_initialized());

  int tag_offset_inter, tag_offset_intra;
  MPI_Comm inter_comm = xt_mpi_comm_smart_dup(inter_comm_, &tag_offset_inter),
    intra_comm = xt_mpi_comm_smart_dup(intra_comm_, &tag_offset_intra);

  struct dist_dir *src_intersections, *dst_intersections;

  bool stripify;
  exchange_idxlists(&src_intersections, &dst_intersections, &stripify,
                    src_idxlist, dst_idxlist,
                    tag_offset_inter, tag_offset_intra,
                    inter_comm, intra_comm, config);

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
               src_idxlist, dst_idxlist, inter_comm);

  xt_mpi_comm_smart_dedup(&inter_comm, tag_offset_inter);
  xt_mpi_comm_smart_dedup(&intra_comm, tag_offset_intra);

  xt_xmdd_free_dist_dir(src_intersections);
  xt_xmdd_free_dist_dir(dst_intersections);
  INSTR_STOP(this_instr);
  return xmap;
}

Xt_xmap
xt_xmap_dist_dir_intercomm_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                               MPI_Comm inter_comm, MPI_Comm intra_comm)
{
  return xt_xmap_dist_dir_intercomm_custom_new(
    src_idxlist, dst_idxlist, inter_comm, intra_comm, &xt_default_config);
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
