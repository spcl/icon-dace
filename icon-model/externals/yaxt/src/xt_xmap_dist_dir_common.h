/**
 * @file xt_xmap_dist_dir_common.h
 *
 * @brief Uitlity functions for creation of distributed directories.
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
#ifndef XT_XMAP_DIST_DIR_COMMON_H
#define XT_XMAP_DIST_DIR_COMMON_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>

#include "core/ppm_visibility.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_xmap_intersection.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct dist_dir {
  int num_entries;
  struct Xt_com_list entries[];
};

PPM_DSO_INTERNAL void
xt_xmdd_free_dist_dir(struct dist_dir *dist_dir);

PPM_DSO_INTERNAL int
xt_com_list_rank_cmp(const void *a_, const void *b_);

struct bucket_params {
  Xt_int global_interval, local_interval,
    local_index_range_lbound, local_index_range_ubound;
};

PPM_DSO_INTERNAL Xt_idxlist
xt_xmap_dist_dir_get_bucket(const struct bucket_params *bucket_params,
                            struct Xt_stripe **stripes_,
                            size_t *stripes_array_size,
                            int dist_dir_rank);

struct Xt_xmdd_txstat {
  size_t bytes, num_msg;
};

/**
 * @brief initiate sends of packed index lists
 * @param send_buffer pointer to buffer of packed index lists
 * @param send_size_asize number of per-rank entries in send_size
 * @param send_size_entry index of @a send_size to act on
 * @param tag tag to use in MPI_Isend
 * @param comm communicator to use
 * @param rank_lim iterate over rank from 0 to rank_lim - 1 when
 * investigating send_size[rank][send_size_entry]
 * @param requests the resulting MPI requests from the issued
 * MPI_Isend calls are stored here, size must be at least number of
 * greater-than-zero entries in @a send_size column @a send_size_entry
 * @param send_size 2-dimensional array of individual message sizes
 */
PPM_DSO_INTERNAL struct Xt_xmdd_txstat
xt_xmap_dist_dir_send_intersections(
  void *restrict send_buffer,
  size_t send_size_asize, size_t send_size_entry,
  int tag, MPI_Comm comm, int rank_lim,
  MPI_Request requests[],
  const int send_size[rank_lim][send_size_asize]);


enum xt_xmdd_direction {
  xt_xmdd_direction_src = 0,
  xt_xmdd_direction_dst = 1,
};

struct isect
{
  int rank[2];
  Xt_idxlist idxlist;
};

PPM_DSO_INTERNAL int
xt_xmdd_cmp_isect_src_rank(const void *a_, const void *b_);
PPM_DSO_INTERNAL int
xt_xmdd_cmp_isect_dst_rank(const void *a_, const void *b_);


PPM_DSO_INTERNAL size_t
xt_xmap_dist_dir_match_src_dst(const struct dist_dir *src_dist_dir,
                               const struct dist_dir *dst_dist_dir,
                               struct isect **src_dst_intersections);

PPM_DSO_INTERNAL size_t
xt_xmap_dist_dir_pack_intersections(
  enum xt_xmdd_direction target,
  size_t num_intersections,
  const struct isect *restrict src_dst_intersections,
  bool isect_idxlist_delete,
  size_t send_size_asize, size_t send_size_idx,
  int (*send_size)[send_size_asize],
  unsigned char *buffer, size_t buf_size, size_t *ofs, MPI_Comm comm);


PPM_DSO_INTERNAL void
xt_xmap_dist_dir_same_rank_merge(struct dist_dir **dist_dir_results);


#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
