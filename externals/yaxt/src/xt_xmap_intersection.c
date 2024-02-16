/**
 * @file xt_xmap_intersection.c
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <mpi.h>

#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_xmap.h"
#include "xt_xmap_internal.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_xmap_intersection.h"
#include "xt_xmap_intersection_common.h"
#include "ensure_array_size.h"
#include "xt_arithmetic_util.h"
#include "xt/quicksort.h"

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ == 11 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static MPI_Comm     xmap_intersection_get_communicator(Xt_xmap xmap);
static int          xmap_intersection_get_num_destinations(Xt_xmap xmap);
static int          xmap_intersection_get_num_sources(Xt_xmap xmap);
static void
xmap_intersection_get_destination_ranks(Xt_xmap xmap, int * ranks);
static void
xmap_intersection_get_source_ranks(Xt_xmap xmap, int * ranks);
static Xt_xmap_iter xmap_intersection_get_in_iterator(Xt_xmap xmap);
static Xt_xmap_iter xmap_intersection_get_out_iterator(Xt_xmap xmap);
static Xt_xmap      xmap_intersection_copy(Xt_xmap xmap);
static void         xmap_intersection_delete(Xt_xmap xmap);
static int          xmap_intersection_iterator_next(Xt_xmap_iter iter);
static int          xmap_intersection_iterator_get_rank(Xt_xmap_iter iter);
static int const *
xmap_intersection_iterator_get_transfer_pos(Xt_xmap_iter iter);
static int
xmap_intersection_iterator_get_num_transfer_pos(Xt_xmap_iter iter);
static const struct Xt_pos_ext *
xmap_intersection_iterator_get_transfer_pos_ext(Xt_xmap_iter iter);
static int
xmap_intersection_iterator_get_num_transfer_pos_ext(Xt_xmap_iter iter);
static void         xmap_intersection_iterator_delete(Xt_xmap_iter iter);
static int          xmap_intersection_get_max_src_pos(Xt_xmap xmap);
static int          xmap_intersection_get_max_dst_pos(Xt_xmap xmap);
static Xt_xmap
xmap_intersection_reorder(Xt_xmap xmap, enum xt_reorder_type type);
static Xt_xmap
xmap_intersection_update_positions(Xt_xmap xmap,
                                   const int * src_positions,
                                   const int * dst_positions);
static Xt_xmap
xmap_intersection_spread(Xt_xmap xmap, int num_repetitions,
                         const int src_displacements[num_repetitions],
                         const int dst_displacements[num_repetitions]);


static const struct Xt_xmap_iter_vtable
xmap_iterator_intersection_vtable = {
  .next                 = xmap_intersection_iterator_next,
  .get_rank             = xmap_intersection_iterator_get_rank,
  .get_transfer_pos     = xmap_intersection_iterator_get_transfer_pos,
  .get_num_transfer_pos = xmap_intersection_iterator_get_num_transfer_pos,
  .get_transfer_pos_ext = xmap_intersection_iterator_get_transfer_pos_ext,
  .get_num_transfer_pos_ext
  = xmap_intersection_iterator_get_num_transfer_pos_ext,
  .delete               = xmap_intersection_iterator_delete};

typedef struct Xt_xmap_iter_intersection_ *Xt_xmap_iter_intersection;

struct Xt_xmap_iter_intersection_ {

  const struct Xt_xmap_iter_vtable * vtable;

  struct exchange_data * msg;
  int msgs_left;
};

static inline Xt_xmap_iter_intersection
xmii(void *iter)
{
  return (Xt_xmap_iter_intersection)iter;
}


static const struct Xt_xmap_vtable xmap_intersection_vtable = {
        .get_communicator      = xmap_intersection_get_communicator,
        .get_num_destinations  = xmap_intersection_get_num_destinations,
        .get_num_sources       = xmap_intersection_get_num_sources,
        .get_destination_ranks = xmap_intersection_get_destination_ranks,
        .get_source_ranks      = xmap_intersection_get_source_ranks,
        .get_out_iterator      = xmap_intersection_get_out_iterator,
        .get_in_iterator       = xmap_intersection_get_in_iterator,
        .copy                  = xmap_intersection_copy,
        .delete                = xmap_intersection_delete,
        .get_max_src_pos       = xmap_intersection_get_max_src_pos,
        .get_max_dst_pos       = xmap_intersection_get_max_dst_pos,
        .reorder               = xmap_intersection_reorder,
        .update_pos            = xmap_intersection_update_positions,
        .spread                = xmap_intersection_spread};

struct exchange_data {
  // list of relative positions in memory to send or receive
  int * transfer_pos;
  struct Xt_pos_ext *transfer_pos_ext_cache;
  int num_transfer_pos, num_transfer_pos_ext;
  int rank;
};

struct Xt_xmap_intersection_ {

  const struct Xt_xmap_vtable * vtable;

  int n_in, n_out;

  // we need the max position in order to enable quick range-checks
  // for xmap-users like redist
  int max_src_pos; // max possible pos over all src transfer_pos (always >= 0)
  int max_dst_pos; // same for dst
  int tag_offset;  /* add to make tags on same communicator non-overlapping */

  MPI_Comm comm;
  struct exchange_data msg[];
};

typedef struct Xt_xmap_intersection_ *Xt_xmap_intersection;

static inline Xt_xmap_intersection
xmi(void *xmap)
{
  return (Xt_xmap_intersection)xmap;
}

static MPI_Comm xmap_intersection_get_communicator(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  return xmap_intersection->comm;
}

static int xmap_intersection_get_num_destinations(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  // the number of destinations equals the number of source messages
  return xmap_intersection->n_out;
}

static int xmap_intersection_get_num_sources(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  // the number of sources equals the number of destination messages
  return xmap_intersection->n_in;
}

static void xmap_intersection_get_destination_ranks(Xt_xmap xmap, int * ranks) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);
  struct exchange_data *restrict out_msg
    = xmap_intersection->msg + xmap_intersection->n_in;
  for (int i = 0; i < xmap_intersection->n_out; ++i)
    ranks[i] = out_msg[i].rank;
}

static void xmap_intersection_get_source_ranks(Xt_xmap xmap, int * ranks) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);
  struct exchange_data *restrict in_msg = xmap_intersection->msg;
  for (int i = 0; i < xmap_intersection->n_in; ++i)
    ranks[i] = in_msg[i].rank;
}

enum {
  bitsPerCoverageElement = sizeof (unsigned long) * CHAR_BIT,
};

struct tpd_result {
  Xt_int *indices_to_remove;
  int resCount;
  bool all_dst_covered;
};

/* compute list positions for recv direction */
static struct tpd_result
generate_dir_transfer_pos_dst(
  int num_intersections,
  const struct Xt_com_list intersections[num_intersections],
  Xt_idxlist mypart_idxlist,
  struct exchange_data *restrict resSets,
  int *restrict num_indices_to_remove_per_intersection)
{
  int mypart_num_indices = xt_idxlist_get_num_indices(mypart_idxlist);
  size_t coverage_size = (size_t)mypart_num_indices;
  coverage_size = (coverage_size + bitsPerCoverageElement - 1)
    /bitsPerCoverageElement;
  unsigned long *restrict coverage = xcalloc(coverage_size, sizeof(*coverage));
  /* set uncovered top-most bits to ease later comparison */
  if (mypart_num_indices%bitsPerCoverageElement)
    coverage[coverage_size-1]
      = ~((1UL << (mypart_num_indices%bitsPerCoverageElement)) - 1UL);

  int new_num_intersections = 0;
  size_t total_num_indices_to_remove = 0;
  size_t curr_indices_to_remove_size = 0;
  Xt_int *restrict indices_to_remove = NULL;
  int *restrict intersection_pos = NULL;

  for (int i = 0; i < num_intersections; ++i) {

    const Xt_int *restrict intersection_idxvec
      = xt_idxlist_get_indices_const(intersections[i].list);
    int max_intersection_size
      = xt_idxlist_get_num_indices(intersections[i].list);
    intersection_pos
      = xrealloc(intersection_pos,
                 (size_t)max_intersection_size * sizeof(*intersection_pos));

#ifndef NDEBUG
    int retval =
#endif
      xt_idxlist_get_positions_of_indices(
        mypart_idxlist, intersection_idxvec, max_intersection_size,
        intersection_pos, 1);
    assert(retval == 0);

    // we have to enforce single_match_only not only within a single
    // intersection, but also between all intersections

    int intersection_size = 0;
    int num_indices_to_remove_isect = 0;

    /* at most max_intersection_size many indices need to be removed */
    ENSURE_ARRAY_SIZE(indices_to_remove, curr_indices_to_remove_size,
                      total_num_indices_to_remove
                      + (size_t)max_intersection_size);

    for (int j = 0; j < max_intersection_size; ++j) {

      int pos = intersection_pos[j];
      /* the predicate effectively conditionalizes adding of either
       * the position to intersection_pos
       *   if the current value was NOT already in another intersection
       * or
       * the index to indices_to_remove_
       *   if the current value was already in another intersection
       */
      unsigned long mask = 1UL << (pos % bitsPerCoverageElement);
      int is_duplicate = (coverage[pos/bitsPerCoverageElement] & mask) != 0UL;
      intersection_pos[intersection_size] = pos;
      indices_to_remove[total_num_indices_to_remove
                        + (size_t)num_indices_to_remove_isect]
        = intersection_idxvec[j];
      intersection_size += is_duplicate ^ 1;
      num_indices_to_remove_isect += is_duplicate;
      coverage[pos/bitsPerCoverageElement] |= mask;
    }

    total_num_indices_to_remove += (size_t)num_indices_to_remove_isect;
    num_indices_to_remove_per_intersection[i] = num_indices_to_remove_isect;

    if (intersection_size > 0) {
      resSets[new_num_intersections].transfer_pos = intersection_pos;
      resSets[new_num_intersections].num_transfer_pos = intersection_size;
      resSets[new_num_intersections].transfer_pos_ext_cache = NULL;
      resSets[new_num_intersections].num_transfer_pos_ext
        = (int)count_pos_ext((size_t)intersection_size, intersection_pos);
      resSets[new_num_intersections].rank = intersections[i].rank;
      new_num_intersections++;
      intersection_pos = NULL;
    }
  }
  free(intersection_pos);

  indices_to_remove
    = xrealloc(indices_to_remove, total_num_indices_to_remove
               * sizeof (*indices_to_remove));

  // check resulting bit map
  unsigned long all_bits_set = ~0UL;
  for (size_t i = 0; i < coverage_size; ++i)
    all_bits_set &= coverage[i];

  free(coverage);
  return (struct tpd_result){
    .indices_to_remove = indices_to_remove,
    .resCount = new_num_intersections,
    .all_dst_covered = all_bits_set == ~0UL };
}


struct pos_count_max {
  int count, max_pos;
};

/* compute list positions for send direction */
static struct pos_count_max
generate_dir_transfer_pos_src(int num_intersections,
                              const struct Xt_com_list
                              intersections[num_intersections],
                              Xt_idxlist mypart_idxlist,
                              struct exchange_data *restrict resSets,
                              const Xt_int *indices_to_remove,
                              const int *num_indices_to_remove_per_intersection)
{
  int new_num_intersections = 0;
  int offset = 0;

  Xt_int * new_intersection_idxvec = NULL;
  size_t curr_new_intersection_idxvec_size = 0;
  int *restrict intersection_pos = NULL;
  int max_pos_ = -1;

  for (int i = 0; i < num_intersections; ++i) {

    const Xt_int *restrict intersection_idxvec
      = xt_idxlist_get_indices_const(intersections[i].list);
    int intersection_size
      = xt_idxlist_get_num_indices(intersections[i].list);
    intersection_pos = xrealloc(intersection_pos,
                                (size_t)intersection_size
                                * sizeof(*intersection_pos));

    int num_indices_to_remove = num_indices_to_remove_per_intersection[i];
    if (num_indices_to_remove > 0) {

      ENSURE_ARRAY_SIZE(
        new_intersection_idxvec, curr_new_intersection_idxvec_size,
        intersection_size - num_indices_to_remove + 1);
      int new_intersection_size = 0;

      for (int j = 0; j < intersection_size; ++j) {

        int discard = 0;

        Xt_int idx = intersection_idxvec[j];
        /* could be improved with BLOOM-filter if
         * num_indices_to_remove was sufficiently large */
        for (int k = 0; k < num_indices_to_remove; ++k)
          discard |= (idx == indices_to_remove[offset + k]);

        new_intersection_idxvec[new_intersection_size] = idx;
        new_intersection_size += !discard;
      }

      intersection_idxvec = new_intersection_idxvec;
      intersection_size = new_intersection_size;
      offset = offset + num_indices_to_remove;
    }

#ifndef NDEBUG
    int retval =
#endif
      xt_idxlist_get_positions_of_indices(
        mypart_idxlist, intersection_idxvec, intersection_size,
        intersection_pos, 0);
    assert(retval == 0);

    if (intersection_size > 0) {
      resSets[new_num_intersections].transfer_pos = intersection_pos;
      resSets[new_num_intersections].num_transfer_pos = intersection_size;
      for (int j = 0; j < intersection_size; ++j)
        if (intersection_pos[j] > max_pos_) max_pos_ = intersection_pos[j];
      resSets[new_num_intersections].transfer_pos_ext_cache = NULL;
      resSets[new_num_intersections].num_transfer_pos_ext
        = (int)count_pos_ext((size_t)intersection_size, intersection_pos);
      resSets[new_num_intersections].rank = intersections[i].rank;
      new_num_intersections++;
      intersection_pos = NULL;
    }
  }

  free(new_intersection_idxvec);
  free(intersection_pos);

  return (struct pos_count_max){ .max_pos = max_pos_,
      .count = new_num_intersections };
}

static Xt_int *
exchange_points_to_remove(int num_src_intersections,
                          const struct Xt_com_list
                          src_com[num_src_intersections],
                          int num_dst_intersections,
                          const struct Xt_com_list
                          dst_com[num_dst_intersections],
                          int *restrict num_src_indices_to_remove_per_intersection,
                          Xt_int *dst_indices_to_remove,
                          const int *restrict
                          num_dst_indices_to_remove_per_intersection,
                          int tag_offset,
                          MPI_Comm comm) {

  MPI_Request * requests
    = xmalloc((size_t)(num_src_intersections + 2 * num_dst_intersections) *
              sizeof(*requests));
  MPI_Request *send_header_requests = requests,
    *recv_requests = requests + num_dst_intersections,
    *send_data_requests = recv_requests + num_src_intersections;

  // set up receives for indices that need to be removed from the send messages
  for (int i = 0; i < num_src_intersections; ++i)
    xt_mpi_call(MPI_Irecv(
                  num_src_indices_to_remove_per_intersection + i, 1, MPI_INT,
                  src_com[i].rank,
                  tag_offset + xt_mpi_tag_xmap_intersection_header_exchange,
                  comm, recv_requests+i), comm);

  // send indices that need to be removed on the target side due to duplicated
  // receives
  int offset = 0;
  unsigned num_nonempty_dst_intersections = 0;
  for (int i = 0; i < num_dst_intersections; ++i) {
    xt_mpi_call(MPI_Isend(
                  CAST_MPI_SEND_BUF(
                    num_dst_indices_to_remove_per_intersection + i), 1, MPI_INT,
                  dst_com[i].rank,
                  tag_offset + xt_mpi_tag_xmap_intersection_header_exchange,
                  comm, send_header_requests + i), comm);

    if (num_dst_indices_to_remove_per_intersection[i] > 0) {

      xt_mpi_call(MPI_Isend(
                    dst_indices_to_remove + offset,
                    num_dst_indices_to_remove_per_intersection[i], Xt_int_dt,
                    dst_com[i].rank,
                    tag_offset + xt_mpi_tag_xmap_intersection_data_exchange,
                    comm, send_data_requests + num_nonempty_dst_intersections),
                  comm);
      offset += num_dst_indices_to_remove_per_intersection[i];
      ++num_nonempty_dst_intersections;
    }
  }

  // wait for the receiving of headers to complete
  xt_mpi_call(MPI_Waitall(num_src_intersections + num_dst_intersections,
                          send_header_requests, MPI_STATUSES_IGNORE), comm);

  size_t total_num_src_indices_to_recv = 0;

  for (int i = 0; i < num_src_intersections; ++i)
    total_num_src_indices_to_recv
      += (size_t)num_src_indices_to_remove_per_intersection[i];

  unsigned num_nonempty_src_intersections = 0;
  Xt_int *src_indices_to_remove;
  if (total_num_src_indices_to_recv > 0) {

    src_indices_to_remove = xmalloc(total_num_src_indices_to_recv
                                    * sizeof(*src_indices_to_remove));

    // set up receive for indices that need to be removed
    offset = 0;
    for (int i = 0; i < num_src_intersections; ++i)
      if (num_src_indices_to_remove_per_intersection[i] > 0) {
        ++num_nonempty_src_intersections;
        xt_mpi_call(MPI_Irecv(
                      src_indices_to_remove + offset,
                      num_src_indices_to_remove_per_intersection[i],
                      Xt_int_dt, src_com[i].rank,
                      tag_offset + xt_mpi_tag_xmap_intersection_data_exchange,
                      comm,
                      send_data_requests - num_nonempty_src_intersections),
                    comm);

        offset += num_src_indices_to_remove_per_intersection[i];
      }

  } else {
    src_indices_to_remove = NULL;
  }

  // wait until all communication is completed
  xt_mpi_call(MPI_Waitall((int)num_nonempty_src_intersections
                          + (int)num_nonempty_dst_intersections,
                          send_data_requests-num_nonempty_src_intersections,
                          MPI_STATUSES_IGNORE), comm);
  free(requests);
  return src_indices_to_remove;
}

static int
generate_transfer_pos(struct Xt_xmap_intersection_ *xmap,
                      int num_src_intersections,
                      const struct Xt_com_list src_com[num_src_intersections],
                      int num_dst_intersections,
                      const struct Xt_com_list dst_com[num_dst_intersections],
                      Xt_idxlist src_idxlist_local,
                      Xt_idxlist dst_idxlist_local,
                      MPI_Comm comm) {

  int *num_src_indices_to_remove_per_intersection =
    xmalloc(((size_t)num_src_intersections + (size_t)num_dst_intersections)
            * sizeof(int)),
    *num_dst_indices_to_remove_per_intersection =
    num_src_indices_to_remove_per_intersection + num_src_intersections;

  struct tpd_result tpdr
    = generate_dir_transfer_pos_dst(
      num_dst_intersections, dst_com, dst_idxlist_local, xmap->msg,
      num_dst_indices_to_remove_per_intersection);
  int all_dst_covered = tpdr.all_dst_covered;
  xmap->n_in = tpdr.resCount;
  Xt_int *dst_indices_to_remove = tpdr.indices_to_remove;
  // exchange the points that need to be removed
  Xt_int *src_indices_to_remove
    = exchange_points_to_remove(
      num_src_intersections, src_com, num_dst_intersections, dst_com,
      num_src_indices_to_remove_per_intersection,
      dst_indices_to_remove, num_dst_indices_to_remove_per_intersection,
      xmap->tag_offset, comm);

  free(dst_indices_to_remove);
  num_src_indices_to_remove_per_intersection
    = xrealloc(num_src_indices_to_remove_per_intersection,
               (size_t)num_src_intersections * sizeof(int));

  struct pos_count_max tpsr
    = generate_dir_transfer_pos_src(
      num_src_intersections, src_com, src_idxlist_local, xmap->msg + xmap->n_in,
      src_indices_to_remove, num_src_indices_to_remove_per_intersection);
  xmap->max_src_pos = tpsr.max_pos;
  xmap->n_out = tpsr.count;

  free(src_indices_to_remove);
  free(num_src_indices_to_remove_per_intersection);
  return all_dst_covered;
}

Xt_xmap
xt_xmap_intersection_new(int num_src_intersections,
                         const struct Xt_com_list
                         src_com[num_src_intersections],
                         int num_dst_intersections,
                         const struct Xt_com_list
                         dst_com[num_dst_intersections],
                         Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                         MPI_Comm comm) {

  // ensure that yaxt is initialized
  assert(xt_initialized());

  size_t num_isect
    = (size_t)num_dst_intersections + (size_t)num_src_intersections;
  Xt_xmap_intersection xmap
    = xmalloc(sizeof (*xmap) + num_isect * sizeof (struct exchange_data));

  xmap->vtable = &xmap_intersection_vtable;

  xmap->comm = comm = xt_mpi_comm_smart_dup(comm, &xmap->tag_offset);

  // generate exchange lists
  if (!generate_transfer_pos(xmap,
                             num_src_intersections, src_com,
                             num_dst_intersections, dst_com,
                             src_idxlist, dst_idxlist, comm)) {

    int num_dst_indices = xt_idxlist_get_num_indices(dst_idxlist);
    const Xt_int *dst_indices = xt_idxlist_get_indices_const(dst_idxlist);
    int * found_index_mask = xcalloc((size_t)num_dst_indices,
                                     sizeof(*found_index_mask));
    int * index_positions = xmalloc((size_t)num_dst_indices
                                    * sizeof(*index_positions));

    for (size_t i = 0; i < (size_t)num_dst_intersections; ++i) {
      xt_idxlist_get_positions_of_indices(dst_com[i].list, dst_indices,
                                          num_dst_indices, index_positions, 0);
      for (size_t j = 0; j < (size_t)num_dst_indices; ++j)
        found_index_mask[j] |= index_positions[j] != -1;
    }

    int first_missing_pos = 0;
    while ((first_missing_pos < (num_dst_indices - 1)) &&
           (found_index_mask[first_missing_pos])) ++first_missing_pos;
    free(found_index_mask);
    free(index_positions);
    xmap_intersection_delete((Xt_xmap)xmap);
    print_miss_msg(dst_idxlist, first_missing_pos, comm, __FILE__, __LINE__);
  }

  int n_in = xmap->n_in, n_out = xmap->n_out;
  size_t new_num_isect = (size_t)n_in + (size_t)n_out;
  if (new_num_isect != num_isect)
    xmap = xrealloc(xmap, sizeof (*xmap) + (new_num_isect
                                            * sizeof(struct exchange_data)));

  xmap->max_dst_pos = xt_idxlist_get_num_indices(dst_idxlist) - 1;

  return (Xt_xmap)xmap;
}

static int xmap_intersection_get_max_src_pos(Xt_xmap xmap) {
  return xmi(xmap)->max_src_pos;
}

static int xmap_intersection_get_max_dst_pos(Xt_xmap xmap) {
  return xmi(xmap)->max_dst_pos;
}


typedef int (*Xt_pos_copy)(size_t num_pos,
                           int *pos, const int *orig_pos,
                           void *state);

static int
pos_copy_verbatim(size_t num_pos,
                  int *pos, const int *orig_pos, void *state)
{
  (void)state;
  memcpy(pos, orig_pos, sizeof (*pos) * num_pos);
  return -1;
}

typedef void (*Xt_pos_ncopy)(size_t num_pos, int *pos, const int *orig_pos,
                             void *state, int num_repetitions,
                             const int displacements[num_repetitions]);

static void
xmap_intersection_msg_copy(size_t nmsg,
                           const struct exchange_data *restrict msg,
                           int *nmsg_copy,
                           struct exchange_data *restrict msg_copy,
                           int *max_pos_, int num_repetitions,
                           Xt_pos_copy pos_copy, void *pos_copy_state) {
  *nmsg_copy = (int)nmsg;
  int max_pos = 0;
  for (size_t i = 0; i < nmsg; ++i) {
    size_t num_transfer_pos
      = (size_t)(msg_copy[i].num_transfer_pos
                 = num_repetitions * msg[i].num_transfer_pos);
    msg_copy[i].rank = msg[i].rank;
    msg_copy[i].transfer_pos_ext_cache = NULL;
    size_t size_transfer_pos
      = num_transfer_pos * sizeof (*(msg[i].transfer_pos));
    msg_copy[i].transfer_pos = xmalloc(size_transfer_pos);
    int new_max_pos
      = pos_copy(num_transfer_pos,
                 msg_copy[i].transfer_pos, msg[i].transfer_pos,
                 pos_copy_state);
    if (new_max_pos > max_pos)
      max_pos = new_max_pos;
    if (pos_copy == pos_copy_verbatim)
      msg_copy[i].num_transfer_pos_ext = msg[i].num_transfer_pos_ext;
    else
      msg_copy[i].num_transfer_pos_ext =
        (int)(count_pos_ext(num_transfer_pos, msg_copy[i].transfer_pos));
  }
  if (pos_copy != pos_copy_verbatim)
    *max_pos_ = max_pos;
}

static Xt_xmap
xmap_intersection_copy_(Xt_xmap xmap,
                        int num_repetitions,
                        Xt_pos_copy pos_copy_in, void *pci_state,
                        Xt_pos_copy pos_copy_out, void *pco_state)
{
  Xt_xmap_intersection xmap_intersection = xmi(xmap), xmap_intersection_new;
  size_t n_in = (size_t)xmap_intersection->n_in,
    n_out = (size_t)xmap_intersection->n_out,
    num_isect = n_in + n_out;
  xmap_intersection_new
    = xmalloc(sizeof (*xmap_intersection_new)
              + num_isect * sizeof (struct exchange_data));
  xmap_intersection_new->vtable = xmap_intersection->vtable;
  xmap_intersection_new->n_in = (int)n_in;
  xmap_intersection_new->n_out = (int)n_out;
  xmap_intersection_new->max_src_pos = xmap_intersection->max_src_pos;
  xmap_intersection_new->max_dst_pos = xmap_intersection->max_dst_pos;
  xmap_intersection_msg_copy(n_in, xmap_intersection->msg,
                             &xmap_intersection_new->n_in,
                             xmap_intersection_new->msg,
                             &xmap_intersection_new->max_dst_pos,
                             num_repetitions,
                             pos_copy_in, pci_state);
  xmap_intersection_msg_copy(n_out, xmap_intersection->msg + n_in,
                             &xmap_intersection_new->n_out,
                             xmap_intersection_new->msg+n_in,
                             &xmap_intersection_new->max_src_pos,
                             num_repetitions,
                             pos_copy_out, pco_state);
  xmap_intersection_new->comm
    = xt_mpi_comm_smart_dup(xmap_intersection->comm,
                            &xmap_intersection_new->tag_offset);
  return (Xt_xmap)xmap_intersection_new;
}

static Xt_xmap
xmap_intersection_copy(Xt_xmap xmap)
{
  return xmap_intersection_copy_(xmap, 1,
                                 pos_copy_verbatim, NULL,
                                 pos_copy_verbatim, NULL);
}

static void
xmap_intersection_msg_delete(int nmsg, struct exchange_data *msg) {
  for (int i = 0; i < nmsg; ++i) {
    free(msg[i].transfer_pos_ext_cache);
    free(msg[i].transfer_pos);
  }
}

static void xmap_intersection_delete(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  int num_isect = xmap_intersection->n_in + xmap_intersection->n_out;
  xmap_intersection_msg_delete(num_isect, xmap_intersection->msg);
  xt_mpi_comm_smart_dedup(&xmap_intersection->comm,
                          xmap_intersection->tag_offset);
  free(xmap_intersection);
}

static Xt_xmap_iter xmap_intersection_get_in_iterator(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  if (xmap_intersection->n_in == 0)
    return NULL;

  Xt_xmap_iter_intersection iter = xmalloc(sizeof (*iter));

  iter->vtable = &xmap_iterator_intersection_vtable;
  iter->msg = xmap_intersection->msg;
  iter->msgs_left = xmap_intersection->n_in - 1;

  return (Xt_xmap_iter)iter;
}

static Xt_xmap_iter xmap_intersection_get_out_iterator(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  if (xmap_intersection->n_out == 0)
    return NULL;

  Xt_xmap_iter_intersection iter = xmalloc(sizeof (*iter));

  iter->vtable = &xmap_iterator_intersection_vtable;
  iter->msg = xmap_intersection->msg + xmap_intersection->n_in;
  iter->msgs_left = xmap_intersection->n_out - 1;

  return (Xt_xmap_iter)iter;
}

struct a2abuf
{
  int *buffer, *counts, *displs;
};

static inline struct a2abuf
setup_buffer(int comm_size, int n, const struct exchange_data *msgs)
{
  size_t buffer_size = 0;
  for (int i = 0; i < n; ++i)
    buffer_size += (size_t)msgs[i].num_transfer_pos;

  int *restrict counts = xmalloc((2 * (size_t)comm_size + buffer_size)
                                 * sizeof(*counts)),
    *restrict displs = counts + comm_size,
    *restrict buffer = counts + 2 * comm_size;
  for (int i = 0; i < comm_size; ++i)
    counts[i] = 0;

  for (int i = 0; i < n; ++i)
    counts[msgs[i].rank] = msgs[i].num_transfer_pos;

  for (int i = 0, accu = 0; i < comm_size; ++i) {
    displs[i] = accu;
    accu += counts[i];
  }
  return (struct a2abuf){ .buffer=buffer, .counts=counts, .displs=displs };
}

static void reorder_transfer_pos(
  int n_out, int n_in, struct exchange_data * out_msg,
  struct exchange_data * in_msg, MPI_Comm comm) {

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  struct a2abuf
    out = setup_buffer(comm_size, n_out, out_msg),
    in = setup_buffer(comm_size, n_in, in_msg);

  for (int i = 0; i < n_out; ++i) {
    const struct exchange_data * curr_msg = out_msg + i;
    int count = curr_msg->num_transfer_pos;
    if (count > 0) {
      // pack local out transfer_pos into out buffer
      memcpy(out.buffer + out.displs[curr_msg->rank],
             curr_msg->transfer_pos, (size_t)count * sizeof(*out.buffer));
      // sort local out transfer_pos
      xt_quicksort_int(curr_msg->transfer_pos, (size_t)count);
    }
  }

  // exchange local out transfer_pos
  xt_mpi_call(MPI_Alltoallv(out.buffer, out.counts, out.displs, MPI_INT,
                            in.buffer, in.counts, in.displs, MPI_INT,
                            comm), comm);

  for (int i = 0; i < n_in; ++i) {
    const struct exchange_data * curr_msg = in_msg + i;
    int count = curr_msg->num_transfer_pos;
    if (count > 0)
      // reorder local receive transfer_pos (curr_msg->transfer_pos) based
      // on remote out transfer_pos (in_buffer + in_displs[curr_msg->rank])
      xt_quicksort_int_permutation(
        in.buffer + in.displs[curr_msg->rank], (size_t)count,
        curr_msg->transfer_pos);
  }

  free(out.counts);
  free(in.counts);

  free(out_msg->transfer_pos_ext_cache);
  free(in_msg->transfer_pos_ext_cache);
  out_msg->transfer_pos_ext_cache = NULL;
  in_msg->transfer_pos_ext_cache = NULL;
  out_msg->num_transfer_pos_ext =
    (int)count_pos_ext(
      (size_t)out_msg->num_transfer_pos, out_msg->transfer_pos);
  in_msg->num_transfer_pos_ext =
    (int)count_pos_ext(
      (size_t)in_msg->num_transfer_pos, in_msg->transfer_pos);
}

static Xt_xmap
xmap_intersection_reorder(Xt_xmap xmap, enum xt_reorder_type type) {

  Xt_xmap_intersection xmap_intersection_new =
    xmi(xmap_intersection_copy(xmap));

  int n_out = xmap_intersection_new->n_out;
  int n_in = xmap_intersection_new->n_in;
  struct exchange_data *in_msg = xmap_intersection_new->msg,
    *out_msg = in_msg+n_in;
  MPI_Comm comm = xmap_intersection_new->comm;

  switch ((int)type) {
    case (XT_REORDER_NONE):
      break;
    case (XT_REORDER_SEND_UP):
      reorder_transfer_pos(n_out, n_in, out_msg, in_msg, comm);
      break;
    case (XT_REORDER_RECV_UP):
      reorder_transfer_pos(n_in, n_out, in_msg, out_msg, comm);
      break;
    default:
      Xt_abort(comm, "ERROR(xmap_intersection_reorder):invalid reorder type",
               __FILE__, __LINE__);
  };

  return (Xt_xmap)xmap_intersection_new;
}

static int
subst_positions(size_t num_pos,
                int *restrict pos, const int *restrict orig_pos, void *new_pos_)
{
  const int *restrict new_pos = new_pos_;
  int max_pos = 0;
  for (size_t i = 0; i < num_pos; ++i) {
    int np = new_pos[orig_pos[i]];
    pos[i] = np;
    if (np > max_pos)
      max_pos = np;
  }
  return max_pos;
}

static Xt_xmap
xmap_intersection_update_positions(Xt_xmap xmap,
                                   const int *src_positions,
                                   const int *dst_positions) {

  return xmap_intersection_copy_(xmap, 1,
                                 subst_positions, (void *)dst_positions,
                                 subst_positions, (void *)src_positions);
}

struct spread_state
{
  int num_repetitions;
  const int *restrict displacements;
};

static int
pos_copy_spread(size_t num_pos, int *restrict pos,
                const int *restrict orig_pos, void *state)
{
  struct spread_state *sp = state;
  int num_repetitions = sp->num_repetitions;
  const int *restrict displacements = sp->displacements;
  num_pos = num_pos / (size_t)num_repetitions;
  int max_pos = -1;
  for (int i = 0, k = 0; i < num_repetitions; ++i) {
    int curr_disp = displacements[i];
    for (size_t j = 0; j < num_pos; ++j, ++k) {
      int np = orig_pos[j] + curr_disp;
      pos[k] = np;
      if (np > max_pos)
        max_pos = np;
    }
  }
  return max_pos;
}


static Xt_xmap
xmap_intersection_spread(Xt_xmap xmap, int num_repetitions,
                         const int src_displacements[num_repetitions],
                         const int dst_displacements[num_repetitions])
{
  return xmap_intersection_copy_(xmap, num_repetitions,
                                 pos_copy_spread,
                                 &(struct spread_state){
                                   .num_repetitions = num_repetitions,
                                   .displacements = src_displacements },
                                 pos_copy_spread,
                                 &(struct spread_state){
                                   .num_repetitions = num_repetitions,
                                   .displacements = dst_displacements });
}

/* how many pos values have monotonically either positively or
 * negatively consecutive values and copy to pos_copy */
static inline struct pos_run copy_get_pos_run_len(
  size_t num_pos, const int *restrict pos,
  int *restrict pos_copy)
{
  size_t i = 0, j = 1;
  int direction = 0;
  int start = pos_copy[0] = pos[0];
  if (j < num_pos) {
    direction = isign_mask(pos[1] - pos[0]);
    while (j < num_pos
           && (pos_copy[j] = pos[j]) == start + (~direction & (int)(j - i)) +
           (direction & -(int)(j - i))) {
      pos_copy[j] = pos[j];
      ++j;
    }
    direction = direction & ((j == 1) - 1);
  }
  return (struct pos_run){ .start = start, .len = j, .direction = direction };
}

/* compute number of position extents that would be required
   to represent positions array and copy to pos_copy */
static struct pos_count_max
max_count_pos_ext_and_copy(int max_pos, size_t num_pos, const int *restrict pos,
                           int *restrict pos_copy)
{
  size_t i = 0, num_pos_ext = 0;
  while (i < num_pos) {
    struct pos_run run = copy_get_pos_run_len(num_pos - i, pos + i, pos_copy + i);
    i += run.len;
    int max_of_run = (run.start & run.direction) | ((run.start + (int)run.len - 1) & ~run.direction);
    if (max_of_run > max_pos) max_pos = max_of_run;
    ++num_pos_ext;
  }
  return (struct pos_count_max){ .count = (int)num_pos_ext, .max_pos = max_pos };
}

static void init_exchange_data_from_com_pos(
  int count, struct exchange_data *restrict msgs,
  const struct Xt_com_pos *restrict com, int *max_pos) {

  int max_pos_ = -1;
  for (int i = 0; i < count; ++i) {
    int num_transfer_pos = com[i].num_transfer_pos;
    int *restrict transfer_pos =
      xmalloc((size_t)num_transfer_pos * sizeof(*transfer_pos));
    msgs[i].transfer_pos = transfer_pos;
    msgs[i].transfer_pos_ext_cache = NULL;
    msgs[i].num_transfer_pos = num_transfer_pos;
    msgs[i].rank = com[i].rank;
    struct pos_count_max max_count
      = max_count_pos_ext_and_copy(max_pos_, (size_t)num_transfer_pos,
                                   com[i].transfer_pos, transfer_pos);
    msgs[i].num_transfer_pos_ext = max_count.count;
    if (max_count.max_pos > max_pos_) max_pos_ = max_count.max_pos;
  }
  *max_pos = max_pos_;
}

Xt_xmap
xt_xmap_intersection_pos_new(
  int num_src_msg, const struct Xt_com_pos src_com[num_src_msg],
  int num_dst_msg, const struct Xt_com_pos dst_com[num_dst_msg],
  MPI_Comm comm) {

  // ensure that yaxt is initialized
  assert(xt_initialized());

  size_t num_msg = (size_t)num_src_msg + (size_t)num_dst_msg;
  Xt_xmap_intersection xmap
    = xmalloc(sizeof (*xmap) + num_msg * sizeof (struct exchange_data));

  xmap->vtable = &xmap_intersection_vtable;
  xmap->n_in = num_dst_msg;
  xmap->n_out = num_src_msg;
  xmap->comm = comm = xt_mpi_comm_smart_dup(comm, &xmap->tag_offset);
  init_exchange_data_from_com_pos(
    num_dst_msg, xmap->msg, dst_com, &xmap->max_dst_pos);
  init_exchange_data_from_com_pos(
    num_src_msg, xmap->msg + num_dst_msg, src_com, &xmap->max_src_pos);

  return (Xt_xmap)xmap;
}

static int xmap_intersection_iterator_next(Xt_xmap_iter iter) {

  Xt_xmap_iter_intersection iter_intersection = xmii(iter);

  if (iter_intersection == NULL || iter_intersection->msgs_left == 0)
    return 0;

  iter_intersection->msg++;
  iter_intersection->msgs_left--;

  return 1;
}

static int xmap_intersection_iterator_get_rank(Xt_xmap_iter iter) {

  assert(iter != NULL);
  return xmii(iter)->msg->rank;
}

static int const *
xmap_intersection_iterator_get_transfer_pos(Xt_xmap_iter iter) {

  assert(iter != NULL);
  return xmii(iter)->msg->transfer_pos;
}

static const struct Xt_pos_ext *
xmap_intersection_iterator_get_transfer_pos_ext(Xt_xmap_iter iter) {

  assert(iter != NULL);
  if (!xmii(iter)->msg->transfer_pos_ext_cache) {
    size_t num_transfer_pos_ext = (size_t)((size_t)xmii(iter)->msg->num_transfer_pos_ext);
    struct Xt_pos_ext * transfer_pos_ext =
      (xmii(iter)->msg->transfer_pos_ext_cache =
         xmalloc(num_transfer_pos_ext * sizeof(*transfer_pos_ext)));
    generate_pos_ext((size_t)xmii(iter)->msg->num_transfer_pos,
                     xmii(iter)->msg->transfer_pos,
                     num_transfer_pos_ext, transfer_pos_ext);
  }
  return xmii(iter)->msg->transfer_pos_ext_cache;
}

static int
xmap_intersection_iterator_get_num_transfer_pos_ext(Xt_xmap_iter iter) {
  assert(iter != NULL);
  return xmii(iter)->msg->num_transfer_pos_ext;
}

static int
xmap_intersection_iterator_get_num_transfer_pos(Xt_xmap_iter iter) {

  assert(iter != NULL);
  return xmii(iter)->msg->num_transfer_pos;
}

static void xmap_intersection_iterator_delete(Xt_xmap_iter iter) {

  free(iter);
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
