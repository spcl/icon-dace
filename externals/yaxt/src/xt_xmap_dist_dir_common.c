/**
 * @file xt_xmap_dist_dir_common.c
 *
 * @brief Implementation of uitlity functions for creation of
 * distributed directories.
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

#include <string.h>
#include <assert.h>

#include "core/ppm_xfuncs.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_mpi.h"
#include "xt/xt_xmap_dist_dir.h"
#include "xt/xt_xmap_dist_dir_intercomm.h"
#include "xt_arithmetic_util.h"
#include "xt_xmap_dist_dir_common.h"
#include "ensure_array_size.h"
#include "xt_idxstripes_internal.h"
#include "xt_config_internal.h"

void xt_xmdd_free_dist_dir(struct dist_dir *dist_dir) {
  size_t num_entries
    = dist_dir->num_entries > 0
    ? (size_t)dist_dir->num_entries : (size_t)0;
  struct Xt_com_list *entries = dist_dir->entries;
  for (size_t i = 0; i < num_entries; ++i)
    xt_idxlist_delete(entries[i].list);
  free(dist_dir);
}

/**
 * \brief generates the buckets of the distributed directory
 *
 * The buckets of the distributed directory are computed as follows:
 * - compute sum of the sizes of all src_idxlist's and round to next
 *   multiple of comm_size -> global_interval
 * - local_interval_size = global_interval / comm_size
 * - local_interval_start = rank * local_interval_size
 * - local_interval_end = (rank + 1) * local_interval_size
 * - local_interval = $[local_interval_start,local_interval_end)$
 * - bucket of each rank is the set of stripes that is defined as:
 *      stripe[i].start = start + rank*local_interval_start + i*global_interval
 *      stripe[i].stride = 1;
 *      stripe[i].nstrides = (int)local_interval;
 *      with i in [0,num_stripes)
 *      bucket[rank] = xt_idxstripes_new(stripes, num_stripes);
 * - num_stripes and start are choosen such that:
 *      start + local_interval
 *         > get_min_idxlist_index(src_idxlist, dst_idxlist)
 *      (which implies for rank 0
 *       stripe[0].start <= get_min_idxlist_index(src_idxlist, dst_idxlist))
 *   where
 *       start = global_interval * k,
 *       k is minimal
 *   and
 *      start + num_stripes*global_interval
 *         > get_max_idxlist_index(src_idxlist, dst_idxlist);
 *   i.e. local_interval is replicated modulo global_interval for the
 *   whole range of indices.
 *
 * \param[in] bucket_params the parameters
 * \param[out] stripes_ pointer to pointer of array of the stripes
 * describing the bucket for rank \a dist_dir_rank. Must point to a
 * pointer allocated by malloc or related and will be resized if necessary.
 * \param[in,out] stripes_array_size number of stripes stored in \a stripes_
 * \param[in] dist_dir_rank rank for which to compute bucket
 * \return newly created bucket. The bucket must not be longer lived
 * than the dynamic memory referenced by \a stripes_ .
 */
Xt_idxlist
xt_xmap_dist_dir_get_bucket(const struct bucket_params *bucket_params,
                            struct Xt_stripe **stripes_,
                            size_t *stripes_array_size,
                            int dist_dir_rank)
{
  Xt_int global_interval = bucket_params->global_interval;
  Xt_int local_interval = bucket_params->local_interval;
  Xt_int local_index_range_lbound = bucket_params->local_index_range_lbound;
  Xt_int local_index_range_ubound = bucket_params->local_index_range_ubound;
  int num_stripes = 0;

  /* find first index in bucket of dist_dir_rank
     <= local_index_range_lbound */
  Xt_int start = (Xt_int)(0 + dist_dir_rank * local_interval);
  {
    long long start_correction
      = (long long)local_index_range_lbound - (long long)start;
    Xt_int corr_steps
      = (Xt_int)((start_correction
                  - (llsign_mask(start_correction)
                     & (long long)(global_interval - 1)))
                 / (long long)global_interval);
    start = (Xt_int)(start + corr_steps * global_interval);
  }
  /* next find last stripe in bucket of dist_dir_rank
   * <= local_index_range_ubound */
  Xt_int end
    = (Xt_int)(start
               + (((long long)local_index_range_ubound - (long long)start)
                  / global_interval) * global_interval);
  Xt_int use_start_stripe
    = (Xt_int)(start + local_interval > local_index_range_lbound);
  num_stripes = (int)(((long long)end - (long long)start)/global_interval)
    + (int)use_start_stripe;
  start = (Xt_int)(start
                   + ((Xt_int)((Xt_uint)use_start_stripe - (Xt_uint)1)
                      & global_interval));
  if (!num_stripes)
    return NULL;

  struct Xt_stripe *restrict stripes = *stripes_;
  ENSURE_ARRAY_SIZE(stripes, *stripes_array_size, (size_t)num_stripes);

  for (int j = 0; j < num_stripes; ++j) {

    stripes[j].start = (Xt_int)(start + j * global_interval);
    stripes[j].stride = 1;
    stripes[j].nstrides = (int)local_interval;
  }

  *stripes_ = stripes;
  return xt_idxstripes_prealloc_new(stripes, num_stripes);
}

struct Xt_xmdd_txstat
xt_xmap_dist_dir_send_intersections(
  void *restrict send_buffer,
  size_t send_size_asize, size_t send_size_entry,
  int tag, MPI_Comm comm, int rank_lim,
  MPI_Request *restrict requests,
  const int send_size[rank_lim][send_size_asize])
{
  size_t offset = 0;
  size_t reqOfs = 0;

  // pack the intersections into the send buffer
  for (size_t rank = 0; rank < (size_t)rank_lim; ++rank)
    if (send_size[rank][send_size_entry] > 0) {
      xt_mpi_call(MPI_Isend((char *)send_buffer + offset,
                            send_size[rank][send_size_entry],
                            MPI_PACKED, (int)rank, tag,
                            comm, requests + reqOfs),
                  comm);
      ++reqOfs;
      offset += (size_t)send_size[rank][send_size_entry];
    }
  return (struct Xt_xmdd_txstat){ .bytes = offset, .num_msg = reqOfs };
}

size_t
xt_xmap_dist_dir_match_src_dst(const struct dist_dir *src_dist_dir,
                               const struct dist_dir *dst_dist_dir,
                               struct isect **src_dst_intersections)
{
  struct isect *src_dst_intersections_ = *src_dst_intersections
    = xmalloc((size_t)src_dist_dir->num_entries
              * (size_t)dst_dist_dir->num_entries
              * sizeof(**src_dst_intersections));
  size_t isect_fill = 0;
  const struct Xt_com_list *restrict entries_src = src_dist_dir->entries,
    *restrict entries_dst = dst_dist_dir->entries;
  size_t num_entries_src = (size_t)src_dist_dir->num_entries,
    num_entries_dst = (size_t)dst_dist_dir->num_entries;
  for (size_t i = 0; i < num_entries_src; ++i)
    for (size_t j = 0; j < num_entries_dst; ++j)
    {
      Xt_idxlist intersection
        = xt_idxlist_get_intersection(entries_src[i].list, entries_dst[j].list);
      if (xt_idxlist_get_num_indices(intersection) > 0) {
        src_dst_intersections_[isect_fill]
          = (struct isect){
          .rank = { [xt_xmdd_direction_src]=entries_src[i].rank,
                    [xt_xmdd_direction_dst]=entries_dst[j].rank},
          .idxlist = intersection };
        ++isect_fill;
      } else
        xt_idxlist_delete(intersection);
    }
  *src_dst_intersections
    = xrealloc(src_dst_intersections_,
               isect_fill * sizeof (*src_dst_intersections_));
  return isect_fill;
}

#if __GNUC__ == 11 && __GNUC_MINOR__ <= 2
/* gcc 11.1 has a bug in the -fsanitize=undefined functionality
 * which creates a bogus warning without the below suppression, bug
 * report at https://gcc.gnu.org/bugzilla/show_bug.cgi?id=101585 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvla-parameter"
#endif

size_t
xt_xmap_dist_dir_pack_intersections(
  enum xt_xmdd_direction target,
  size_t num_intersections,
  const struct isect *restrict src_dst_intersections,
  bool isect_idxlist_delete,
  size_t send_size_asize, size_t send_size_idx,
  int (*send_size)[send_size_asize],
  unsigned char *buffer, size_t buf_size, size_t *ofs, MPI_Comm comm)
{
  int prev_send_rank = -1;
  size_t num_send_indices_requests = 0;
  size_t origin = 1 ^ target, ofs_ = *ofs;
  int position = 0;
  for (size_t i = 0; i < num_intersections; ++i)
  {
    /* see if this generates a new request? */
    int send_rank = src_dst_intersections[i].rank[target];
    num_send_indices_requests += send_rank != prev_send_rank;

    // pack rank
    XT_MPI_SEND_BUF_CONST int *prank
      = CAST_MPI_SEND_BUF(src_dst_intersections[i].rank + origin);
    if (send_rank != prev_send_rank && prev_send_rank != -1) {
      send_size[prev_send_rank][send_size_idx] = position;
      ofs_ += (size_t)position;
      position = 0;
    }
    prev_send_rank = send_rank;
    xt_mpi_call(MPI_Pack(prank, 1, MPI_INT, buffer+ofs_, (int)(buf_size-ofs_),
                         &position, comm), comm);
    // pack intersection
    xt_idxlist_pack(src_dst_intersections[i].idxlist, buffer+ofs_,
                    (int)(buf_size-ofs_), &position, comm);

    if (isect_idxlist_delete)
      xt_idxlist_delete(src_dst_intersections[i].idxlist);
  }
  if (prev_send_rank != -1)
    send_size[prev_send_rank][send_size_idx] = position;

  *ofs = ofs_ + (size_t)position;
  return num_send_indices_requests;
}

#if __GNUC__ == 11 && __GNUC_MINOR__ <= 2
#pragma GCC diagnostic pop
#endif


static int
stripe_cmp(const void *a, const void *b)
{
  typedef const struct Xt_stripe *csx;
  return (((csx)a)->start > ((csx)b)->start)
    - (((csx)b)->start > ((csx)a)->start);
}

/*
 * @param dist_dir_results contains the intersections of this ranks
 * dst or src idxlist with other ranks in bucket-sized chunks, the
 * chunks belonging to the same communication partner are merged
 * in-place here, i.e. on return (*dist_dir_results)->num_entries is
 * less than or equal to the previous count and *dist_dir_results
 * might point somewhere else
 */
void
xt_xmap_dist_dir_same_rank_merge(struct dist_dir **dist_dir_results) {

  struct Xt_com_list *restrict entries = (*dist_dir_results)->entries;
  size_t num_isect_agg = 0;

  size_t i = 0, num_shards = (size_t)(*dist_dir_results)->num_entries;
  while (i < num_shards) {
    int rank = entries[i].rank;
    size_t j = i;
    /* find all entries matching the currently considered rank */
    do
      ++j;
    while (j < num_shards && entries[j].rank == rank);

    struct Xt_stripe *restrict stripes = NULL;
    size_t num_stripes = 0;
    for (; i < j; ++i) {
      struct Xt_stripe *stripes_of_intersection;
      int num_stripes_of_intersection;
      xt_idxlist_get_index_stripes(entries[i].list,
                                   &stripes_of_intersection,
                                   &num_stripes_of_intersection);
      xt_idxlist_delete(entries[i].list);
      if (stripes) {
        stripes = xrealloc(stripes,
                           (num_stripes + (size_t)num_stripes_of_intersection)
                           * sizeof (*stripes));
        memcpy(stripes + num_stripes, stripes_of_intersection,
               (size_t)num_stripes_of_intersection * sizeof (*stripes));
        free(stripes_of_intersection);
      } else
        stripes = stripes_of_intersection;
      num_stripes += (size_t)num_stripes_of_intersection;
    }
    qsort(stripes, num_stripes, sizeof (*stripes), stripe_cmp);
    entries[num_isect_agg].list = xt_idxstripes_new(stripes, (int)num_stripes);
    free(stripes);
    entries[num_isect_agg].rank = rank;
    ++num_isect_agg;
  }
  (*dist_dir_results)->num_entries = (int)num_isect_agg;
  *dist_dir_results = xrealloc(*dist_dir_results, sizeof (struct dist_dir)
                               + (size_t)num_isect_agg
                               * sizeof(struct Xt_com_list));
}


int
xt_xmdd_cmp_isect_src_rank(const void *a_, const void *b_)
{
  const struct isect *a = a_, *b = b_;
  /* this is safe vs. overflow because ranks are in [0..MAX_INT) */
  return a->rank[xt_xmdd_direction_src] - b->rank[xt_xmdd_direction_src];
}

int
xt_xmdd_cmp_isect_dst_rank(const void *a_, const void *b_)
{
  const struct isect *a = a_, *b = b_;
  /* this is safe vs. overflow because ranks are in [0..MAX_INT) */
  return a->rank[xt_xmdd_direction_dst] - b->rank[xt_xmdd_direction_dst];
}

int
xt_com_list_rank_cmp(const void *a_, const void *b_)
{
  const struct Xt_com_list *a = a_, *b = b_;
  /* this is overflow-safe because rank's are non-negative ints */
  return a->rank - b->rank;
}

Xt_xmap
xt_xmap_dist_dir_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                     MPI_Comm comm)
{
  return xt_xmap_dist_dir_custom_new(src_idxlist, dst_idxlist, comm,
                                     &xt_default_config);
}

Xt_xmap
xt_xmap_dist_dir_custom_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                            MPI_Comm comm, Xt_config config)
{
  // ensure that yaxt is initialized
  assert(xt_initialized());

  int is_inter;
  xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
  Xt_xmap xmap;
  if (!is_inter)
    xmap = xt_xmap_dist_dir_intracomm_custom_new(src_idxlist, dst_idxlist, comm,
                                                 config);
  else {
    MPI_Comm merge_comm, local_intra_comm;
    MPI_Group local_group;
    xt_mpi_call(MPI_Comm_group(comm, &local_group), comm);
    xt_mpi_call(MPI_Intercomm_merge(comm, 0, &merge_comm), comm);
    xt_mpi_call(MPI_Comm_create(merge_comm, local_group, &local_intra_comm),
                comm);
    xt_mpi_call(MPI_Group_free(&local_group), comm);
    xt_mpi_call(MPI_Comm_free(&merge_comm), comm);
    xt_mpi_comm_mark_exclusive(local_intra_comm);

    xmap
      = xt_xmap_dist_dir_intercomm_custom_new(src_idxlist, dst_idxlist,
                                              comm, local_intra_comm, config);

    xt_mpi_call(MPI_Comm_free(&local_intra_comm), local_intra_comm);
  }
  return xmap;
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
