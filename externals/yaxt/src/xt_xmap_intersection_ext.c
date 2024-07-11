/**
 * @file xt_xmap_intersection_ext.c
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

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
#include "xt_arithmetic_util.h"
#include "ensure_array_size.h"
#include "xt_cover.h"
#include "xt/quicksort.h"

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static MPI_Comm     xmap_intersection_ext_get_communicator(Xt_xmap xmap);
static int          xmap_intersection_ext_get_num_destinations(Xt_xmap xmap);
static int          xmap_intersection_ext_get_num_sources(Xt_xmap xmap);
static void
xmap_intersection_ext_get_destination_ranks(Xt_xmap xmap, int * ranks);
static void
xmap_intersection_ext_get_source_ranks(Xt_xmap xmap, int * ranks);
static Xt_xmap_iter xmap_intersection_ext_get_in_iterator(Xt_xmap xmap);
static Xt_xmap_iter xmap_intersection_ext_get_out_iterator(Xt_xmap xmap);
static Xt_xmap xmap_intersection_ext_copy(Xt_xmap xmap);
static void xmap_intersection_ext_delete(Xt_xmap xmap);
static int xmap_intersection_ext_get_max_src_pos(Xt_xmap xmap);
static int xmap_intersection_ext_get_max_dst_pos(Xt_xmap xmap);
static Xt_xmap
xmap_intersection_ext_reorder(Xt_xmap xmap, enum xt_reorder_type type);
static Xt_xmap
xmap_intersection_ext_update_positions(Xt_xmap xmap,
                                       const int * src_positions,
                                       const int * dst_positions);
static Xt_xmap
xmap_intersection_ext_spread(Xt_xmap xmap, int num_repetitions,
                             const int src_displacements[num_repetitions],
                             const int dst_displacements[num_repetitions]);


static const struct Xt_xmap_vtable xmap_intersection_vtable = {
        .get_communicator      = xmap_intersection_ext_get_communicator,
        .get_num_destinations  = xmap_intersection_ext_get_num_destinations,
        .get_num_sources       = xmap_intersection_ext_get_num_sources,
        .get_destination_ranks = xmap_intersection_ext_get_destination_ranks,
        .get_source_ranks      = xmap_intersection_ext_get_source_ranks,
        .get_out_iterator      = xmap_intersection_ext_get_out_iterator,
        .get_in_iterator       = xmap_intersection_ext_get_in_iterator,
        .copy                  = xmap_intersection_ext_copy,
        .delete                = xmap_intersection_ext_delete,
        .get_max_src_pos       = xmap_intersection_ext_get_max_src_pos,
        .get_max_dst_pos       = xmap_intersection_ext_get_max_dst_pos,
        .reorder               = xmap_intersection_ext_reorder,
        .update_pos            = xmap_intersection_ext_update_positions,
        .spread                = xmap_intersection_ext_spread};

struct exchange_ext {
  // list of relative position extents in index list to send or receive
  struct Xt_pos_ext *transfer_pos_ext;
  /* generated on-demand */
  int *transfer_pos;
  int num_transfer_pos, num_transfer_pos_ext;
  int rank;
};

struct Xt_xmap_intersection_ext_ {

  const struct Xt_xmap_vtable * vtable;

  int n_in, n_out;

  // we need the max position in order to enable quick range-checks
  // for xmap-users like redist
  int max_src_pos; // max possible pos over all src transfer_pos (always >= 0)
  int max_dst_pos; // same for dst
  int tag_offset;  /* offset to add to message tags for uniqueness */
  MPI_Comm comm;
  struct exchange_ext msg[];
};

typedef struct Xt_xmap_intersection_ext_ *Xt_xmap_intersection_ext;

static inline Xt_xmap_intersection_ext xmie(void *xmap)
{
  return (Xt_xmap_intersection_ext)xmap;
}

static MPI_Comm xmap_intersection_ext_get_communicator(Xt_xmap xmap)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  return xmap_intersection_ext->comm;
}

static int xmap_intersection_ext_get_num_destinations(Xt_xmap xmap)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  // the number of destination equals the number of source messages
  return xmap_intersection_ext->n_out;
}

static int xmap_intersection_ext_get_num_sources(Xt_xmap xmap)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  // the number of sources equals the number of destination messages
  return xmap_intersection_ext->n_in;
}

static void
xmap_intersection_ext_get_destination_ranks(Xt_xmap xmap, int *restrict ranks)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  size_t n_out = (size_t)xmap_intersection_ext->n_out;
  const struct exchange_ext *restrict out_msg
    = xmap_intersection_ext->msg + xmap_intersection_ext->n_in;
  for (size_t i = 0; i < n_out; ++i)
    ranks[i] = out_msg[i].rank;
}

static void
xmap_intersection_ext_get_source_ranks(Xt_xmap xmap, int *restrict ranks)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  size_t n_in = (size_t)xmap_intersection_ext->n_in;
  const struct exchange_ext *restrict in_msg = xmap_intersection_ext->msg;
  for (size_t i = 0; i < n_in; ++i)
    ranks[i] = in_msg[i].rank;
}

static int xmap_intersection_ext_get_max_src_pos(Xt_xmap xmap) {
  return xmie(xmap)->max_src_pos;
}

static int xmap_intersection_ext_get_max_dst_pos(Xt_xmap xmap) {
  return xmie(xmap)->max_dst_pos;
}

typedef int (*Xt_pos_ext_copy)(size_t num_orig_pos_ext,
                               size_t *num_pos_ext,
                               struct Xt_pos_ext **pos_ext,
                               const struct Xt_pos_ext *orig_pos_ext,
                               size_t num_orig_pos, const int *orig_pos,
                               void *state);

static int
pos_ext_copy_verbatim(size_t num_orig_pos_ext,
                      size_t *num_pos_ext,
                      struct Xt_pos_ext **pos_ext,
                      const struct Xt_pos_ext *orig_pos_ext,
                      size_t num_orig_pos, const int *orig_pos,
                      void *state)
{
  (void)state; (void)num_orig_pos; (void)orig_pos;
  size_t size_pos_ext = num_orig_pos_ext * sizeof (**pos_ext);
  struct Xt_pos_ext *pos_ext_ = *pos_ext = xmalloc(size_pos_ext);
  memcpy(pos_ext_, orig_pos_ext, size_pos_ext);
  *num_pos_ext = num_orig_pos_ext;
  return -1;
}

static void
xmap_intersection_ext_msg_copy(size_t nmsg,
                               const struct exchange_ext *restrict msg,
                               int *nmsg_copy,
                               struct exchange_ext *restrict msg_copy,
                               int *max_pos_, int num_repetitions,
                               Xt_pos_ext_copy pos_ext_copy, void *pec_state)
{
  *nmsg_copy = (int)nmsg;
  int max_pos = 0;
  for (size_t i = 0; i < nmsg; ++i) {
    msg_copy[i].num_transfer_pos = num_repetitions * msg[i].num_transfer_pos;
    msg_copy[i].rank = msg[i].rank;
    msg_copy[i].transfer_pos = NULL;
    size_t num_transfer_pos_ext;
    int new_max_pos
      = pos_ext_copy((size_t)msg[i].num_transfer_pos_ext, &num_transfer_pos_ext,
                     &msg_copy[i].transfer_pos_ext, msg[i].transfer_pos_ext,
                     (size_t)msg[i].num_transfer_pos, msg[i].transfer_pos,
                     pec_state);
    if (new_max_pos > max_pos)
      max_pos = new_max_pos;
    msg_copy[i].num_transfer_pos_ext = (int)num_transfer_pos_ext;
  }
  if (pos_ext_copy != pos_ext_copy_verbatim)
    *max_pos_ = max_pos;
}

static Xt_xmap
xmap_intersection_ext_copy_(Xt_xmap xmap, int num_repetitions,
                            Xt_pos_ext_copy pe_cpy_in, void *peci_state,
                            Xt_pos_ext_copy pe_cpy_out, void *peco_state)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap),
    xmap_intersection_ext_new;
  size_t n_in = (size_t)xmap_intersection_ext->n_in,
    n_out = (size_t)xmap_intersection_ext->n_out,
    num_isect = n_in + n_out;
  xmap_intersection_ext_new
    = xmalloc(sizeof (*xmap_intersection_ext_new)
              + num_isect * sizeof (struct exchange_ext));
  xmap_intersection_ext_new->vtable = xmap_intersection_ext->vtable;
  xmap_intersection_ext_new->n_in = (int)n_in;
  xmap_intersection_ext_new->n_out = (int)n_out;
  xmap_intersection_ext_new->max_src_pos = xmap_intersection_ext->max_src_pos;
  xmap_intersection_ext_new->max_dst_pos = xmap_intersection_ext->max_dst_pos;
  xmap_intersection_ext_msg_copy(n_in, xmap_intersection_ext->msg,
                                 &xmap_intersection_ext_new->n_in,
                                 xmap_intersection_ext_new->msg,
                                 &xmap_intersection_ext_new->max_dst_pos,
                                 num_repetitions, pe_cpy_in, peci_state);
  xmap_intersection_ext_msg_copy(n_out, xmap_intersection_ext->msg+n_in,
                                 &xmap_intersection_ext_new->n_out,
                                 xmap_intersection_ext_new->msg+n_in,
                                 &xmap_intersection_ext_new->max_src_pos,
                                 num_repetitions, pe_cpy_out, peco_state);
  xmap_intersection_ext_new->comm
    = xt_mpi_comm_smart_dup(xmap_intersection_ext->comm,
                            &xmap_intersection_ext_new->tag_offset);
  return (Xt_xmap)xmap_intersection_ext_new;
}

static Xt_xmap
xmap_intersection_ext_copy(Xt_xmap xmap)
{
  return xmap_intersection_ext_copy_(xmap, 1,
                                     pos_ext_copy_verbatim, NULL,
                                     pos_ext_copy_verbatim, NULL);
}


static void
xt_free_exchange_ext(size_t num_msg, struct exchange_ext *restrict msg)
{
  for (size_t i = 0; i < num_msg; ++i) {
    free(msg[i].transfer_pos);
    free(msg[i].transfer_pos_ext);
  }
}

static void xmap_intersection_ext_delete(Xt_xmap xmap) {

  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  size_t num_isect = (size_t)xmap_intersection_ext->n_in
    + (size_t)xmap_intersection_ext->n_out;
  xt_free_exchange_ext(num_isect, xmap_intersection_ext->msg);
  xt_mpi_comm_smart_dedup(&xmap_intersection_ext->comm,
                          xmap_intersection_ext->tag_offset);
  free(xmap_intersection_ext);
}

static void
generate_transfer_ext(struct Xt_xmap_intersection_ext_ *xmap,
                      int num_src_intersections,
                      const struct Xt_com_list src_com[num_src_intersections],
                      int num_dst_intersections,
                      const struct Xt_com_list dst_com[num_dst_intersections],
                      Xt_idxlist src_idxlist_local,
                      Xt_idxlist dst_idxlist_local,
                      MPI_Comm comm);

Xt_xmap
xt_xmap_intersection_ext_new(int num_src_intersections,
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
  Xt_xmap_intersection_ext xmap
    = xmalloc(sizeof (*xmap) + num_isect * sizeof (struct exchange_ext));

  xmap->vtable = &xmap_intersection_vtable;

  xmap->comm = comm = xt_mpi_comm_smart_dup(comm, &xmap->tag_offset);

  // generate exchange lists
  generate_transfer_ext(xmap,
                        num_src_intersections, src_com,
                        num_dst_intersections, dst_com,
                        src_idxlist, dst_idxlist, comm);

  size_t new_num_isect = (size_t)xmap->n_in + (size_t)xmap->n_out;
  if (new_num_isect != num_isect)
    xmap = xrealloc(xmap, sizeof (*xmap) + (new_num_isect
                                            * sizeof(struct exchange_ext)));

  xmap->max_dst_pos = xt_idxlist_get_num_indices(dst_idxlist) - 1;

  return (Xt_xmap)xmap;
}

struct ted_result {
  struct Xt_pos_ext_vec cover;
  int resCount;
};

static struct ted_result
generate_dir_transfer_pos_ext_dst(
  int num_intersections,
  const struct Xt_com_list intersections[num_intersections],
  Xt_idxlist mypart_idxlist,
  struct exchange_ext *resSets,
  int (*restrict dst_removals_per_intersection)[2]);


static struct Xt_pos_ext *
exchange_pos_ext_modifications(
  int num_src_intersections,
  const struct Xt_com_list src_com[num_src_intersections],
  int num_dst_intersections,
  const struct Xt_com_list dst_com[num_dst_intersections],
  struct exchange_ext dst_ext[num_dst_intersections],
  int (*restrict src_removals_per_intersection)[2],
  const int (*restrict dst_removals_per_intersection)[2],
  int tag_offset,
  MPI_Comm comm);

static void
remap_dst_intersections(int num_dst_intersections,
                        const struct Xt_com_list dst_com[num_dst_intersections],
                        Xt_idxlist mypart_idxlist,
                        int resCount,
                        struct exchange_ext resSets[resCount],
                        const int (*removals_per_intersection)[2]);

struct tes_result {
  int resCount;
  int max_pos;
};

static struct tes_result
generate_dir_transfer_pos_ext_src(
  int num_intersections,
  const struct Xt_com_list intersections[num_intersections],
  Xt_idxlist mypart_idxlist,
  struct exchange_ext *resSets,
  const int (*restrict removals_per_intersection)[2],
  const struct Xt_pos_ext *pos_updates);

static void
generate_transfer_ext(struct Xt_xmap_intersection_ext_ *xmap,
                      int num_src_intersections,
                      const struct Xt_com_list src_com[num_src_intersections],
                      int num_dst_intersections,
                      const struct Xt_com_list dst_com[num_dst_intersections],
                      Xt_idxlist src_idxlist,
                      Xt_idxlist dst_idxlist,
                      MPI_Comm comm) {

  /* {dst|src}_removals_per_intersection[i][0] denotes the number of
   * indices to be removed from the intersection with {src|dst}_com[i].rank.
   * {dst|src}_removals_per_intersection[rank][1] denotes the number of
   * pos_ext needed to represent this change (0 if either none or all
   * indices got removed).
   */
  int (*src_removals_per_intersection)[2] =
    xmalloc(((size_t)num_dst_intersections + (size_t)num_src_intersections)
            * sizeof(*src_removals_per_intersection)),
    (*dst_removals_per_intersection)[2]
    = src_removals_per_intersection + num_src_intersections;

  {
    struct ted_result tedr
      = generate_dir_transfer_pos_ext_dst(
        num_dst_intersections, dst_com, dst_idxlist,
        xmap->msg, dst_removals_per_intersection);
    struct Xt_pos_ext_vec cover = tedr.cover;
    if (!xt_idxlist_pos_ext_is_full_cover(dst_idxlist, cover)) {
      if (xt_idxlist_get_num_indices(dst_idxlist) == 0)
        Xt_abort(comm, "ERROR: ups...this should not have happend...", __FILE__,
                 __LINE__);
      int first_missing_pos
        = ((cover.num_pos_ext > 0) && (cover.pos_ext[0].start == 0))
        ? cover.pos_ext[0].start + cover.pos_ext[0].size : 0;
      print_miss_msg(dst_idxlist, first_missing_pos, comm, __FILE__, __LINE__);
    }
    xt_cover_finish(&cover);
    xmap->n_in = tedr.resCount;
  }

  // exchange pos_ext of lists where additional indices need to be removed
  struct Xt_pos_ext *pos_updates
    = exchange_pos_ext_modifications(
      num_src_intersections, src_com, num_dst_intersections, dst_com, xmap->msg,
      src_removals_per_intersection,
      (const int (*)[2])dst_removals_per_intersection, xmap->tag_offset, comm);

  remap_dst_intersections(num_dst_intersections, dst_com, dst_idxlist,
                          xmap->n_in, xmap->msg,
                          (const int (*)[2])dst_removals_per_intersection);

  src_removals_per_intersection =
    xrealloc(src_removals_per_intersection, (size_t)num_src_intersections
             * sizeof(*src_removals_per_intersection));

  struct tes_result tesr
    = generate_dir_transfer_pos_ext_src(
      num_src_intersections, src_com, src_idxlist, xmap->msg+xmap->n_in,
      (const int (*)[2])src_removals_per_intersection, pos_updates);
  xmap->n_out = tesr.resCount;
  xmap->max_src_pos = tesr.max_pos;
  free(src_removals_per_intersection);
  free(pos_updates);
}

struct Xt_pos_ext_overlap {
  int skip, overlap, tail;
};

static struct Xt_pos_ext_overlap
Xt_get_pos_ext_overlap(struct Xt_pos_ext a, struct Xt_pos_ext b)
{
  /* == 0 if a.size >= 0 ; == ~0 if a.size < 0 */
  int aSizeMaskNeg = isign_mask(a.size),
    /* compute start and end indices of ranges */
    a_s = a.start +  (aSizeMaskNeg & (a.size + 1)),
    a_e   = a.start + (~aSizeMaskNeg & (a.size - 1)),
    bSizeMaskNeg = isign_mask(b.size),
    b_s = b.start +  (bSizeMaskNeg & (b.size + 1)),
    b_e = b.start + (~bSizeMaskNeg & (b.size - 1));
  /* does overlap exist? */
  if ((b_s > a_e) | (a_s > b_e))
    return (struct Xt_pos_ext_overlap){ a.size, 0, 0};
  else {
    /* determine length of overlap parts */
    int lowSkipA = b_s - a_s;
    int lowSkipB = -lowSkipA;
    lowSkipA = (int)((unsigned)(lowSkipA + abs(lowSkipA))/2U);
    lowSkipB = (int)((unsigned)(lowSkipB + abs(lowSkipB))/2U);
    int overlapLen = imin(b_e - b_s - lowSkipB + 1,
                          abs(a.size) - lowSkipA);
    int highSkipA = abs(a.size) - lowSkipA - overlapLen;
    /* then adjust lengths to direction of overlap (from
     * perspective of a */
    int aSkipLen = (~aSizeMaskNeg & lowSkipA)
      | (aSizeMaskNeg & -highSkipA),
      aTailLen = (aSizeMaskNeg & -lowSkipA)
      | (~aSizeMaskNeg & highSkipA);
    return (struct Xt_pos_ext_overlap){ aSkipLen, overlapLen, aTailLen };
  }
}



static void
cut_pos_ext_from_pos_exts(struct Xt_pos_ext pos_ext,
                          struct Xt_pos_ext_vec *pos_exts);

static struct Xt_pos_ext *
get_pos_exts_of_index_stripes(Xt_idxlist idxlist,
                              int num_stripes,
                              const struct Xt_stripe stripes[num_stripes],
                              int *num_ext,
                              int single_match_only)
{
  struct Xt_pos_ext *pos_ext;
#ifndef NDEBUG
  int retval =
#endif
    xt_idxlist_get_pos_exts_of_index_stripes(
      idxlist, num_stripes, stripes,
      num_ext, &pos_ext, single_match_only);
  assert(retval == 0);
  return pos_ext;
}

static struct ted_result
generate_dir_transfer_pos_ext_dst(
  int num_intersections,
  const struct Xt_com_list intersections[num_intersections],
  Xt_idxlist mypart_idxlist,
  struct exchange_ext *restrict resSets,
  int (*restrict dst_removals_per_intersection)[2])
{
  int new_num_intersections = 0;

  // we have to enforce single_match_only not only within a single
  // intersection, but also between all intersections
  /* ranges already covered from previous intersections, i.e. which
   * must not be transmitted twice */
  enum { initial_vec_size = 8 };
  struct Xt_pos_ext_vec cover;
  xt_cover_start(&cover, initial_vec_size);

  for (int i = 0; i < num_intersections; ++i) {

    int num_stripes, num_indices_to_remove = 0;
    struct Xt_stripe *intersection_idxstripes;
    xt_idxlist_get_index_stripes(intersections[i].list,
                                 &intersection_idxstripes, &num_stripes);
    int num_isect_pos_exts;
    struct Xt_pos_ext *restrict isect_pos_exts
      = get_pos_exts_of_index_stripes(
        mypart_idxlist, num_stripes, intersection_idxstripes,
        &num_isect_pos_exts, 1);
    int isect_pos_exts_size_psum = 0;
    int intersection_size = xt_idxlist_get_num_indices(intersections[i].list);
    /* start with all indices from intersection as used,
       later split ranges, if overlaps are found */
    struct Xt_pos_ext_vec transferable
      = (struct Xt_pos_ext_vec){ .num_pos_ext = 1,
          .size_pos_ext = initial_vec_size,
          .pos_ext = xrealloc(intersection_idxstripes,
                              sizeof (struct Xt_pos_ext) * initial_vec_size) };
    intersection_idxstripes = NULL;
    transferable.pos_ext[0]
      = (struct Xt_pos_ext){ .start = 0, .size = intersection_size };
    /* find overlap(s) with previously found ranges for all
     * stripes mapped to position extents */
    for (size_t j = 0; j < (size_t)num_isect_pos_exts; ++j) {
      struct Xt_pos_ext isect_pos_ext = isect_pos_exts[j];
      /* ensure isect_pos_ext is oriented with ascending positions */
      int isign_mask_isect_pos_ext_size = isign_mask(isect_pos_ext.size);
      isect_pos_ext.start
        += isign_mask_isect_pos_ext_size & (isect_pos_ext.size + 1);
      int isect_pos_ext_orig_size = isect_pos_ext.size;
      isect_pos_ext.size = abs(isect_pos_ext.size);
      isect_pos_exts_size_psum += isect_pos_ext.size;
      /* keep progress as inverse of change to psum to compensate for
       * eventual correction later */
      int progress = -isect_pos_ext.size;
      size_t search_start_pos = 0, insert_pos;
      do {
        struct Xt_pos_range query = (struct Xt_pos_range){
          .start = isect_pos_ext.start,
          .end = isect_pos_ext.start + isect_pos_ext.size - 1 };
        insert_pos
          = xt_cover_insert_or_overlap(&cover, query, true, search_start_pos);
        if (insert_pos == SIZE_MAX)
          goto next_isect_pos_ext;
        struct Xt_pos_ext_overlap overlap_desc
          = Xt_get_pos_ext_overlap(isect_pos_ext, cover.pos_ext[insert_pos]);
        /* insert overlap into updates
         * by ...*/
        /* ...first inserting the skipped part into
         * cover.pos_ext, since that is sorted
         * and obviously precedes cover.pos_ext[insert_pos],
         * cover.pos_ext[insert_pos] can be seemlessly extended...
         */
        cover.pos_ext[insert_pos].start -= overlap_desc.skip;
        cover.pos_ext[insert_pos].size += overlap_desc.skip;
        /* ...and optionally merged with its predecessor, if the
         * intervening range becomes zero, ... */
        if (insert_pos > 0
            && (cover.pos_ext[insert_pos].start
                == (cover.pos_ext[insert_pos - 1].start
                    + cover.pos_ext[insert_pos - 1].size)))
        {
          cover.pos_ext[insert_pos - 1].size
            += cover.pos_ext[insert_pos].size;
          memmove(cover.pos_ext + insert_pos, cover.pos_ext + insert_pos + 1,
                  (--cover.num_pos_ext - insert_pos)
                  * sizeof (*cover.pos_ext));
          --insert_pos;
        }
        progress = (~isign_mask_isect_pos_ext_size
                    & (progress + overlap_desc.skip))
          | (isign_mask_isect_pos_ext_size
             & (isect_pos_ext_orig_size + overlap_desc.tail));
        /* ... then splitting transferable accordingly, ... */
        num_indices_to_remove += overlap_desc.overlap;
        cut_pos_ext_from_pos_exts((struct Xt_pos_ext){
            .start = isect_pos_exts_size_psum + progress,
              .size = overlap_desc.overlap }, &transferable);
        progress += overlap_desc.overlap;
        /* ... lastly the search can continue with the tail ... */
        isect_pos_ext.start += overlap_desc.skip + overlap_desc.overlap;
        /* ... if there is any */
        isect_pos_ext.size = overlap_desc.tail;
        search_start_pos = ++insert_pos;
      } while ((isect_pos_ext.size != 0)
               & (search_start_pos != cover.num_pos_ext));
      if (isect_pos_ext.size)
        /* already at end of list -> append ... */
        xt_cover_range_append(&cover, isect_pos_ext);
      /* ... and start the next intersection range */
    next_isect_pos_ext:
      ;
    }

    if (intersection_size > num_indices_to_remove) {
      resSets[new_num_intersections].transfer_pos_ext
        = xrealloc(transferable.pos_ext, sizeof (struct Xt_pos_ext)
                   * transferable.num_pos_ext);
      /* start with empty cache of positions to transfer */
      resSets[new_num_intersections].transfer_pos = NULL;
      resSets[new_num_intersections].num_transfer_pos
        = intersection_size - num_indices_to_remove;
      resSets[new_num_intersections].num_transfer_pos_ext
        = (int)transferable.num_pos_ext;
      resSets[new_num_intersections].rank = intersections[i].rank;
      ++new_num_intersections;
    } else
      free(transferable.pos_ext);
    dst_removals_per_intersection[i][0] = num_indices_to_remove;
    dst_removals_per_intersection[i][1]
      = ((num_indices_to_remove == intersection_size)
         | (num_indices_to_remove == 0))?0:(int)transferable.num_pos_ext;
    free(isect_pos_exts);
  }
  /* since cover is a struct, at least pgcc 11-13 cannot compile this with a
   * compound literal or initialize r directly */
#if defined __PGI && __PGIC__ <= 13
  struct ted_result r;
  r.cover = cover;
  r.resCount = new_num_intersections;
  return r;
#else
  return (struct ted_result){
    .cover = cover,
    .resCount = new_num_intersections };
#endif
}

static void
cut_pos_ext_from_pos_exts(struct Xt_pos_ext pos_ext,
                          struct Xt_pos_ext_vec *pos_exts)
{
  struct Xt_pos_ext *restrict pos_exts_ = pos_exts->pos_ext;
  size_t num_pos_exts_ = pos_exts->num_pos_ext;
  size_t i = num_pos_exts_;
  while (pos_exts_[--i].start > pos_ext.start)
    ;
  int db_skip = pos_ext.start - pos_exts_[i].start;
  if ((!db_skip) & (pos_ext.size == pos_exts_[i].size))
  {
    /* delete fully overlapped transfer part */
    memmove(pos_exts_ + i, pos_exts_ + i + 1,
            sizeof (*pos_exts_) * (num_pos_exts_ - i - 1));
    pos_exts->num_pos_ext = --num_pos_exts_;
  }
  else if (db_skip + pos_ext.size == pos_exts_[i].size)
  {
    /* pos_ext overlaps end of pos_exts_[i] */
    pos_exts_[i].size -= pos_ext.size;
  }
  else if (db_skip == 0)
  {
    /* pos_ext overlaps start of pos_exts_[i] */
    pos_exts_[i].start = pos_ext.start + pos_ext.size;
    pos_exts_[i].size -= pos_ext.size;
  }
  else
  {
    struct Xt_pos_ext orig = pos_exts_[i];
    ENSURE_ARRAY_SIZE(pos_exts->pos_ext, pos_exts->size_pos_ext,
                      num_pos_exts_ + 1);
    pos_exts_ = pos_exts->pos_ext;
    memmove(pos_exts_ + i + 1, pos_exts_ + i,
            (num_pos_exts_ - i) * sizeof (*pos_exts_));
    pos_exts_[i] = (struct Xt_pos_ext){.start = orig.start,
                                       .size = db_skip };
    pos_exts_[i + 1] = (struct Xt_pos_ext){
      .start = pos_ext.start + pos_ext.size,
      .size = orig.size - db_skip - pos_ext.size };
    pos_exts->num_pos_ext = ++num_pos_exts_;
  }
}

static struct Xt_pos_ext *
exchange_pos_ext_modifications(
  int num_src_intersections,
  const struct Xt_com_list src_com[num_src_intersections],
  int num_dst_intersections,
  const struct Xt_com_list dst_com[num_dst_intersections],
  struct exchange_ext dst_ext[num_dst_intersections],
  int (*restrict src_removals_per_intersection)[2],
  const int (*restrict dst_removals_per_intersection)[2],
  int tag_offset,
  MPI_Comm comm)
{
  MPI_Request * requests
    = xmalloc((size_t)(num_src_intersections + 2 * num_dst_intersections) *
              sizeof(*requests));
  MPI_Request *restrict send_header_requests = requests,
    *restrict recv_requests = requests + num_dst_intersections,
    *restrict send_data_requests = recv_requests + num_src_intersections;

  // set up receives for indices that need to be removed from the send messages
  for (int i = 0; i < num_src_intersections; ++i)
    xt_mpi_call(MPI_Irecv(
                  src_removals_per_intersection[i], 2, MPI_INT, src_com[i].rank,
                  tag_offset + xt_mpi_tag_xmap_intersection_header_exchange,
                  comm, recv_requests + i), comm);

  /* send rebuilt pos_ext vectors that needed to be modified on the
   * target side due to duplicated receives
   */
  unsigned num_active_dst = 0, num_dst_changes = 0;
  for (int i = 0; i < num_dst_intersections; ++i) {
    xt_mpi_call(MPI_Isend(
                  CAST_MPI_SEND_BUF(dst_removals_per_intersection[i]),
                  2, MPI_INT, dst_com[i].rank,
                  tag_offset + xt_mpi_tag_xmap_intersection_header_exchange,
                  comm, send_header_requests + i), comm);

    if (dst_removals_per_intersection[i][1] > 0) {

      assert(dst_removals_per_intersection[i][1]
             == dst_ext[num_active_dst].num_transfer_pos_ext
             && dst_com[i].rank == dst_ext[num_active_dst].rank);
      xt_mpi_call(MPI_Isend(
                    dst_ext[num_active_dst].transfer_pos_ext,
                    dst_removals_per_intersection[i][1],
                    MPI_2INT, dst_com[i].rank,
                    tag_offset + xt_mpi_tag_xmap_intersection_data_exchange,
                    comm, send_data_requests + num_dst_changes),
                  comm);
      ++num_dst_changes;
    }
    num_active_dst += (unsigned)((dst_removals_per_intersection[i][0] == 0)
                                 | (dst_removals_per_intersection[i][1] != 0));
  }

  // wait for the receiving of headers to complete
  xt_mpi_call(MPI_Waitall(num_src_intersections + num_dst_intersections,
                          requests, MPI_STATUSES_IGNORE), comm);

  size_t total_num_pos_ext_to_recv = 0;

  for (size_t i = 0; i < (size_t)num_src_intersections; ++i)
    total_num_pos_ext_to_recv += (size_t)src_removals_per_intersection[i][1];

  struct Xt_pos_ext *src_updated_pos_ext;
  unsigned num_src_changes = 0;
  if (total_num_pos_ext_to_recv > 0) {

    src_updated_pos_ext
      = xmalloc(total_num_pos_ext_to_recv * sizeof(*src_updated_pos_ext));

    /* set up receive for pos_ext that need to be modified because
     * indices needed to be removed from the intersection */
    size_t offset = 0;
    for (int i = 0; i < num_src_intersections; ++i)
      if (src_removals_per_intersection[i][1] > 0) {
        ++num_src_changes;
        xt_mpi_call(MPI_Irecv(
                      src_updated_pos_ext + offset,
                      src_removals_per_intersection[i][1], MPI_2INT,
                      src_com[i].rank,
                      tag_offset + xt_mpi_tag_xmap_intersection_data_exchange,
                      comm, send_data_requests - num_src_changes), comm);

        offset += (size_t)src_removals_per_intersection[i][1];
      }
  } else
    src_updated_pos_ext = NULL;

  // wait until all communication is completed
  xt_mpi_call(MPI_Waitall((int)num_src_changes + (int)num_dst_changes,
                          send_data_requests - num_src_changes,
                          MPI_STATUSES_IGNORE), comm);

  free(requests);
  return src_updated_pos_ext;
}

static void
remap_intersection(Xt_idxlist mypart_idxlist,
                   Xt_idxlist intersection,
                   size_t num_pos_updates,
                   const struct Xt_pos_ext pos_updates[num_pos_updates],
                   struct exchange_ext *resSet,
                   int single_match_only);

static void
remap_dst_intersections(
  int num_dst_intersections,
  const struct Xt_com_list intersections[num_dst_intersections],
  Xt_idxlist mypart_idxlist,
  int resCount,
  struct exchange_ext resSets[resCount],
  const int (*removals_per_intersection)[2])
{
  size_t resIdx = 0;
  for (size_t i = 0; i < (size_t)num_dst_intersections; ++i)
  {
    int intersection_size
      = xt_idxlist_get_num_indices(intersections[i].list);

    int num_indices_to_remove = removals_per_intersection[i][0];

    if (num_indices_to_remove != intersection_size) {} else
      /* intersection is made redundant */
      continue;

    struct Xt_pos_ext *pos_updates = resSets[resIdx].transfer_pos_ext;
    remap_intersection(mypart_idxlist, intersections[i].list,
                       (size_t)removals_per_intersection[i][1],
                       pos_updates, resSets + resIdx, 1);
    free(pos_updates);
    ++resIdx;
  }
  assert(resIdx == (size_t)resCount);
}

static inline int
pos_ext_find_max_pos(int num_pos_ext,
                     const struct Xt_pos_ext *restrict pos_ext)
{
  int max_pos = -1;
  for (size_t i = 0; i < (size_t)num_pos_ext; ++i) {
    int start = pos_ext[i].start,
      size = pos_ext[i].size,
      max = size > 0 ? start + size - 1 : start;
    if (max > max_pos) max_pos = max;
  }
  return max_pos;
}

/* compute updated positions for send direction */
static struct tes_result
generate_dir_transfer_pos_ext_src(
  int num_intersections,
  const struct Xt_com_list intersections[num_intersections],
  Xt_idxlist mypart_idxlist,
  struct exchange_ext *resSets,
  const int (*restrict removals_per_intersection)[2],
  const struct Xt_pos_ext *pos_updates)
{
  /* count total number of intersections that results */
  int new_num_intersections = 0;
  /* indexes into pos_updates */
  size_t intersection_pos_ext = 0;

  int max_pos = -1;
  for (int i = 0; i < num_intersections; ++i) {

    int intersection_size
      = xt_idxlist_get_num_indices(intersections[i].list);

    int num_indices_to_remove = removals_per_intersection[i][0];

    /* intersection might be redundant */
    if (num_indices_to_remove != intersection_size) {

      remap_intersection(mypart_idxlist, intersections[i].list,
                         (size_t)removals_per_intersection[i][1],
                         pos_updates + intersection_pos_ext,
                         resSets + new_num_intersections, 0);

      int max = pos_ext_find_max_pos(
        resSets[new_num_intersections].num_transfer_pos_ext,
        resSets[new_num_intersections].transfer_pos_ext);
      if (max > max_pos) max_pos = max;
      /* evaluate cache lazily */
      resSets[new_num_intersections].transfer_pos = NULL;
      resSets[new_num_intersections].num_transfer_pos
        = intersection_size - num_indices_to_remove;
      resSets[new_num_intersections].rank = intersections[i].rank;
      new_num_intersections++;
      intersection_pos_ext += (size_t)removals_per_intersection[i][1];
    }
  }

  return (struct tes_result) {
    .resCount = new_num_intersections,
    .max_pos = max_pos,
  };
}

static struct Xt_stripe *
refine_stripes(int *num_stripes_,
               struct Xt_stripe *restrict intersection_idxstripes,
               size_t num_pos_updates,
               const struct Xt_pos_ext *restrict pos_updates)
{
  /* trim intersection_idxstripes to those actually used */
  size_t num_refined_intersection_idxstripes = 0,
    size_refined_intersection_idxstripes = num_pos_updates;
  struct Xt_stripe *restrict refined_intersection_idxstripes
    = xmalloc(size_refined_intersection_idxstripes
              * sizeof (*refined_intersection_idxstripes));
  size_t i_stripe = 0;
  int nstrides_psum = 0;
  for (size_t i_pos_ext = 0; i_pos_ext < num_pos_updates; ++i_pos_ext)
  {
    int pos = pos_updates[i_pos_ext].start;
    int size = pos_updates[i_pos_ext].size;
    while (nstrides_psum + intersection_idxstripes[i_stripe].nstrides <= pos)
    {
      nstrides_psum += intersection_idxstripes[i_stripe].nstrides;
      ++i_stripe;
    }
    do {
      int instripe_pos = pos - nstrides_psum;
      ENSURE_ARRAY_SIZE(refined_intersection_idxstripes,
                        size_refined_intersection_idxstripes,
                        num_refined_intersection_idxstripes + 1);
      struct Xt_stripe cur_stripe = intersection_idxstripes[i_stripe];
      int cur_stripe_nstrides = cur_stripe.nstrides;
      int overlap = imin(cur_stripe_nstrides - instripe_pos, size);
      cur_stripe.start
        = (Xt_int)(cur_stripe.start
                   + (Xt_int)instripe_pos * cur_stripe.stride);
      cur_stripe.nstrides = overlap;
      refined_intersection_idxstripes[num_refined_intersection_idxstripes]
        = cur_stripe;
      ++num_refined_intersection_idxstripes;
      i_stripe += (instripe_pos + overlap == cur_stripe_nstrides);
      nstrides_psum += (instripe_pos + overlap == cur_stripe_nstrides)
        ? cur_stripe_nstrides : 0;
      pos += overlap;
      size -= overlap;
    } while (size);
  }
  free(intersection_idxstripes);
  *num_stripes_ = (int)num_refined_intersection_idxstripes;
  return refined_intersection_idxstripes;
}


/* match index stripes of intersection to corresponding positions in
 * part list, optionally updating the stripes
 * @param num_pos_updates number of position extents describing the
 *                        subset of positions from intersections to use
 * @param pos_updates list of position extents to use from \a intersection
 */
static void
remap_intersection(Xt_idxlist mypart_idxlist,
                   Xt_idxlist intersection,
                   size_t num_pos_updates,
                   const struct Xt_pos_ext pos_updates[num_pos_updates],
                   struct exchange_ext *resSet,
                   int single_match_only)
{
  struct Xt_stripe *intersection_idxstripes;
  int num_stripes;
  xt_idxlist_get_index_stripes(intersection,
                               &intersection_idxstripes,
                               &num_stripes);
  if (num_pos_updates)
    intersection_idxstripes
      = refine_stripes(&num_stripes, intersection_idxstripes,
                       num_pos_updates, pos_updates);

  /* match back intersection_idxstripes to positions in mypart */
  resSet->transfer_pos_ext = NULL;
#ifndef NDEBUG
  int retval =
#endif
    xt_idxlist_get_pos_exts_of_index_stripes(
      mypart_idxlist, num_stripes, intersection_idxstripes,
      &resSet->num_transfer_pos_ext,
      &resSet->transfer_pos_ext, single_match_only);
  assert(retval == 0);
  free(intersection_idxstripes);
}

static struct Xt_pos_ext * exchange_transfer_pos_ext(
  int n_out, const struct exchange_ext *restrict out_msg,
  int n_in, const struct exchange_ext *restrict in_msg,
  struct exchange_ext *restrict remote_out_msg,
  int tag_offset, MPI_Comm comm) {

  MPI_Request * requests
    = xmalloc((size_t)(n_in + 2 * n_out) * sizeof(*requests));
  MPI_Request *send_header_requests = requests,
    *recv_requests = requests + n_out,
    *send_data_requests = recv_requests + n_in;

  // set up receives for number of transfer_pos_ext per remote
  // for (int i = 0; i < n_in; ++i)
  for (int i = 0; i < n_in; ++i)
    xt_mpi_call(MPI_Irecv(
                  &(remote_out_msg[i].num_transfer_pos_ext), 1, MPI_INT,
                  in_msg[i].rank,
                  tag_offset + xt_mpi_tag_xmap_intersection_header_exchange,
                  comm, recv_requests + i), comm);

  // send number of transfer_pos_ext
  for (int i = 0; i < n_out; ++i) {
    xt_mpi_call(MPI_Isend(
                  CAST_MPI_SEND_BUF(&(out_msg[i].num_transfer_pos_ext)),
                  1, MPI_INT, out_msg[i].rank,
                  tag_offset + xt_mpi_tag_xmap_intersection_header_exchange,
                  comm, send_header_requests + i), comm);

    xt_mpi_call(MPI_Isend(
                  CAST_MPI_SEND_BUF(out_msg[i].transfer_pos_ext),
                  out_msg[i].num_transfer_pos_ext,
                  MPI_2INT, out_msg[i].rank,
                  tag_offset + xt_mpi_tag_xmap_intersection_data_exchange,
                  comm, send_data_requests + i),
                comm);
  }

  // wait for the receiving of headers to complete
  xt_mpi_call(
    MPI_Waitall(n_out + n_in, send_header_requests, MPI_STATUSES_IGNORE), comm);

  size_t total_num_pos_ext_to_recv = 0;

  for (size_t i = 0; i < (size_t)n_in; ++i)
    total_num_pos_ext_to_recv +=
      (size_t)(remote_out_msg[i].num_transfer_pos_ext);

  struct Xt_pos_ext *transfer_pos_ext_buffer;
  if (total_num_pos_ext_to_recv > 0) {

    transfer_pos_ext_buffer
      = xmalloc(total_num_pos_ext_to_recv * sizeof(*transfer_pos_ext_buffer));

    // set up receive for transfer_pos_ext
    struct Xt_pos_ext *curr_transfer_pos_ext = transfer_pos_ext_buffer;
    for (int i = 0; i < n_in; ++i) {
      xt_mpi_call(MPI_Irecv(
                    curr_transfer_pos_ext,
                    remote_out_msg[i].num_transfer_pos_ext, MPI_2INT,
                    in_msg[i].rank,
                    tag_offset + xt_mpi_tag_xmap_intersection_data_exchange,
                    comm, recv_requests + i), comm);

      remote_out_msg[i].transfer_pos_ext = curr_transfer_pos_ext;
      curr_transfer_pos_ext += remote_out_msg[i].num_transfer_pos_ext;
    }
  } else
    transfer_pos_ext_buffer = NULL;

  // wait until all communication is completed
  xt_mpi_call(
    MPI_Waitall(n_in + n_out, recv_requests, MPI_STATUSES_IGNORE), comm);

  free(requests);
  return transfer_pos_ext_buffer;
}

static void
sort_transfer_pos_ext(int n, struct exchange_ext *msg) {

  int buffer_size = 0;
  for (int i = 0; i < n; ++i)
    if (msg[i].transfer_pos == NULL && msg[i].num_transfer_pos > buffer_size)
      buffer_size = msg[i].num_transfer_pos;

  int *transfer_pos_buffer
    = buffer_size > 0
    ? xmalloc((size_t)buffer_size * sizeof(*transfer_pos_buffer))
    : NULL;

  for (int i = 0; i < n; ++i) {

    // get the positions array
    int *restrict transfer_pos;
    size_t num_transfer_pos = (size_t)(msg[i].num_transfer_pos);
    if (msg[i].transfer_pos != NULL) {
      transfer_pos = msg[i].transfer_pos;
    } else {
      transfer_pos = transfer_pos_buffer;
      generate_pos(
        (size_t)(msg[i].num_transfer_pos_ext), msg[i].transfer_pos_ext,
        num_transfer_pos, transfer_pos);
    }

    // sort the positions array
    xt_quicksort_int(transfer_pos, num_transfer_pos);

    // convert the positions array to position extents array
    size_t num_transfer_pos_ext = count_pos_ext(num_transfer_pos, transfer_pos);
    struct Xt_pos_ext * transfer_pos_ext;
    if (num_transfer_pos_ext != (size_t)(msg[i].num_transfer_pos_ext)) {
      msg[i].num_transfer_pos_ext = (int)num_transfer_pos_ext;
      transfer_pos_ext =
        (msg[i].transfer_pos_ext =
           xrealloc(msg[i].transfer_pos_ext,
                    num_transfer_pos_ext * sizeof(*transfer_pos_ext)));
    } else {
      transfer_pos_ext = msg[i].transfer_pos_ext;
    }
    generate_pos_ext(
      num_transfer_pos, transfer_pos, num_transfer_pos_ext, transfer_pos_ext);
  }

  if (buffer_size > 0) free(transfer_pos_buffer);
}

static void
sort_transfer_pos_ext_permutation(int n,
                                  struct exchange_ext *msg,
                                  struct exchange_ext *permutation_msg) {

  size_t buffer_size = 0;
  for (int i = 0; i < n; ++i) {
    assert(msg[i].transfer_pos == NULL
           && permutation_msg[i].transfer_pos == NULL);
    size_t curr_buffer_size
      = (size_t)(msg[i].num_transfer_pos)
      + (size_t)(permutation_msg[i].num_transfer_pos);
    if (curr_buffer_size > buffer_size) buffer_size = curr_buffer_size;
  }

  int *transfer_pos_buffer
    = buffer_size > 0
    ? xmalloc(buffer_size * sizeof(*transfer_pos_buffer))
    : NULL;

  for (int i = 0; i < n; ++i) {

    // get the positions array
    size_t num_transfer_pos = (size_t)(msg[i].num_transfer_pos);

    int *restrict transfer_pos = transfer_pos_buffer;
    generate_pos(
      (size_t)(msg[i].num_transfer_pos_ext), msg[i].transfer_pos_ext,
      num_transfer_pos, transfer_pos);

    // get the permutation array
    int *permutation = transfer_pos_buffer + num_transfer_pos;
    generate_pos(
      (size_t)(permutation_msg[i].num_transfer_pos_ext),
      permutation_msg[i].transfer_pos_ext, num_transfer_pos, permutation);

    // sort the positions array
    xt_quicksort_int_permutation(transfer_pos, num_transfer_pos, permutation);

    // convert the permutation array to position extents array
    size_t num_transfer_pos_ext = count_pos_ext(num_transfer_pos, permutation);
    struct Xt_pos_ext * transfer_pos_ext;
    if (num_transfer_pos_ext !=
        (size_t)(permutation_msg[i].num_transfer_pos_ext)) {
      permutation_msg[i].num_transfer_pos_ext = (int)num_transfer_pos_ext;
      permutation_msg[i].transfer_pos_ext
        = transfer_pos_ext
        = xrealloc(permutation_msg[i].transfer_pos_ext,
                   num_transfer_pos_ext * sizeof(*transfer_pos_ext));
    } else {
      transfer_pos_ext = permutation_msg[i].transfer_pos_ext;
    }
    generate_pos_ext(
      num_transfer_pos, permutation, num_transfer_pos_ext, transfer_pos_ext);
  }

  if (buffer_size > 0) free(transfer_pos_buffer);
}

static void reorder_transfer_pos_ext(
  int n_out, int n_in, struct exchange_ext * out_msg,
  struct exchange_ext * in_msg, int tag_offset, MPI_Comm comm) {

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  struct exchange_ext *restrict remote_out_msg =
    xmalloc((size_t)n_in * sizeof(*remote_out_msg));

  for (int i = 0; i < n_in; ++i) {
    remote_out_msg[i].rank = in_msg[i].rank;
    remote_out_msg[i].num_transfer_pos = in_msg[i].num_transfer_pos;
    remote_out_msg[i].transfer_pos = NULL;
  }

  // exchange transfer_pos_ext of out messages to respective receivers
  struct Xt_pos_ext *transfer_pos_ext_buffer =
    exchange_transfer_pos_ext(
      n_out, out_msg, n_in, in_msg, remote_out_msg, tag_offset, comm);

  // sort the transfer positions in all out messages
  sort_transfer_pos_ext(n_out, out_msg);

  // sort the transfer positions in all in messages based on the remote out
  // messages
  sort_transfer_pos_ext_permutation(n_in, remote_out_msg, in_msg);

  free(transfer_pos_ext_buffer);
  free(remote_out_msg);
}

static Xt_xmap
xmap_intersection_ext_reorder(Xt_xmap xmap, enum xt_reorder_type type) {

  Xt_xmap_intersection_ext xmap_intersection_ext_new =
    xmie(xmap_intersection_ext_copy(xmap));

  int n_out = xmap_intersection_ext_new->n_out;
  int n_in = xmap_intersection_ext_new->n_in;
  struct exchange_ext *in_msg = xmap_intersection_ext_new->msg,
    *out_msg = in_msg + n_in;
  MPI_Comm comm = xmap_intersection_ext_new->comm;
  int tag_offset = xmap_intersection_ext_new->tag_offset;

  switch ((int)type) {
    case (XT_REORDER_NONE):
      break;
    case (XT_REORDER_SEND_UP):
      reorder_transfer_pos_ext(n_out, n_in, out_msg, in_msg, tag_offset, comm);
      break;
    case (XT_REORDER_RECV_UP):
      reorder_transfer_pos_ext(n_in, n_out, in_msg, out_msg, tag_offset, comm);
      break;
    default:
      Xt_abort(comm, "ERROR(xmap_intersection_ext_reorder):invalid reorder "
               "type", __FILE__, __LINE__);
  };

  return (Xt_xmap)xmap_intersection_ext_new;
}

struct up_state
{
  int *pos_buffer;
  const int *new_pos;
};

static int
update_positions(size_t num_orig_pos_ext,
                 size_t *num_pos_ext,
                 struct Xt_pos_ext **pos_ext,
                 const struct Xt_pos_ext *orig_pos_ext,
                 size_t num_orig_pos, const int *orig_pos,
                 void *state_)
{
  (void)num_orig_pos_ext;
  struct up_state *state = state_;
  int *pos_buffer = state->pos_buffer;
  const int *restrict new_pos = state->new_pos;
  const int *pos;
  if (orig_pos)
    pos = orig_pos;
  else {
    generate_pos(num_orig_pos_ext, orig_pos_ext, num_orig_pos, pos_buffer);
    pos = pos_buffer;
  }
  // update positions
  int max_pos = 0;
  for (size_t j = 0; j < num_orig_pos; ++j) {
    int np = new_pos[pos[j]];
    pos_buffer[j] = np;
    if (np > max_pos)
      max_pos = np;
  }

  // convert the array of substituted positions into position extents array
  size_t num_pos_ext_ = *num_pos_ext = count_pos_ext(num_orig_pos, pos_buffer);
  struct Xt_pos_ext *pos_ext_
    = *pos_ext = xmalloc(num_pos_ext_ * sizeof (*pos_ext_));
  generate_pos_ext(num_orig_pos, pos_buffer, num_pos_ext_, pos_ext_);
  return max_pos;
}

static Xt_xmap
xmap_intersection_ext_update_positions(Xt_xmap xmap,
                                       const int *src_positions,
                                       const int *dst_positions) {

  Xt_xmap_intersection_ext xmie_orig = xmie(xmap);
  size_t max_num_pos = 0;
  size_t n = (size_t)xmie_orig->n_in + (size_t)xmie_orig->n_out;
  const struct exchange_ext *restrict msg = xmie_orig->msg;
  for (size_t i = 0; i < n; ++i)
    if ((size_t)msg[i].num_transfer_pos > max_num_pos)
      max_num_pos = (size_t)msg[i].num_transfer_pos;
  int *pos_buffer
    = max_num_pos > 0
    ? xmalloc((size_t)max_num_pos * sizeof(*pos_buffer))
    : NULL;
  struct up_state ups_in = { pos_buffer, dst_positions },
    ups_out = { pos_buffer, src_positions };
  Xt_xmap xmap_new
    = xmap_intersection_ext_copy_(xmap, 1,
                                  update_positions, &ups_in,
                                  update_positions, &ups_out);
  free(pos_buffer);
  return xmap_new;
}

struct spread_state
{
  int num_repetitions;
  const int *displacements;
};

static int
pos_ext_copy_spread(size_t num_orig_pos_ext,
                    size_t *num_pos_ext,
                    struct Xt_pos_ext **pos_ext,
                    const struct Xt_pos_ext *orig_pos_ext,
                    size_t num_orig_pos, const int *orig_pos,
                    void *state)
{
  (void)num_orig_pos; (void)orig_pos;
  struct spread_state *sp = state;
  int num_repetitions = sp->num_repetitions;
  const int *restrict displacements = sp->displacements;
  size_t new_num_pos_ext = num_orig_pos_ext * (size_t)num_repetitions;
  size_t size_pos_ext = new_num_pos_ext * sizeof (**pos_ext);
  struct Xt_pos_ext *pos_ext_ = *pos_ext = xmalloc(size_pos_ext);
  int max_pos = -1;
  for (int i = 0; i < num_repetitions; ++i) {
    struct Xt_pos_ext *restrict curr_pos_ext =
      pos_ext_ + (size_t)i * num_orig_pos_ext;
    const int curr_displacement = displacements[i];
    for (size_t j = 0; j < num_orig_pos_ext; ++j) {
      int start = orig_pos_ext[j].start + curr_displacement,
        size = orig_pos_ext[j].size,
        max = size > 0 ? start + size - 1 : start;
      if (max > max_pos)
        max_pos = max;
      curr_pos_ext[j] = (struct Xt_pos_ext){ .start = start, .size = size };
    }
  }
  *num_pos_ext = new_num_pos_ext;
  return max_pos;
}


static Xt_xmap
xmap_intersection_ext_spread(Xt_xmap xmap, int num_repetitions,
                             const int src_displacements[num_repetitions],
                             const int dst_displacements[num_repetitions]) {


  return xmap_intersection_ext_copy_(xmap, num_repetitions,
                                     pos_ext_copy_spread,
                                     &(struct spread_state){
                                       .num_repetitions = num_repetitions,
                                       .displacements = src_displacements },
                                     pos_ext_copy_spread,
                                     &(struct spread_state){
                                       .num_repetitions = num_repetitions,
                                       .displacements = dst_displacements });
}


/* iterator operations */

static int xmap_intersection_ext_iterator_next(Xt_xmap_iter iter);
static int xmap_intersection_ext_iterator_get_rank(Xt_xmap_iter iter);
static int const *
xmap_intersection_ext_iterator_get_transfer_pos(Xt_xmap_iter iter);
static int
xmap_intersection_ext_iterator_get_num_transfer_pos(Xt_xmap_iter iter);
static const struct Xt_pos_ext *
xmap_intersection_ext_iterator_get_transfer_pos_ext(Xt_xmap_iter iter);
static int
xmap_intersection_ext_iterator_get_num_transfer_pos_ext(Xt_xmap_iter iter);
static void xmap_intersection_ext_iterator_delete(Xt_xmap_iter iter);

static const struct Xt_xmap_iter_vtable
xmap_iterator_intersection_ext_vtable = {
  .next                 = xmap_intersection_ext_iterator_next,
  .get_rank             = xmap_intersection_ext_iterator_get_rank,
  .get_transfer_pos     = xmap_intersection_ext_iterator_get_transfer_pos,
  .get_num_transfer_pos = xmap_intersection_ext_iterator_get_num_transfer_pos,
  .get_transfer_pos_ext = xmap_intersection_ext_iterator_get_transfer_pos_ext,
  .get_num_transfer_pos_ext
  = xmap_intersection_ext_iterator_get_num_transfer_pos_ext,
  .delete               = xmap_intersection_ext_iterator_delete};

typedef struct Xt_xmap_iter_intersection_ext_ *Xt_xmap_iter_intersection_ext;

struct Xt_xmap_iter_intersection_ext_ {

  const struct Xt_xmap_iter_vtable *vtable;

  struct exchange_ext *msg;
  int msgs_left;
};

static Xt_xmap_iter xmap_intersection_ext_get_in_iterator(Xt_xmap xmap) {

  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);

  if (xmap_intersection_ext->n_in == 0)
    return NULL;

  Xt_xmap_iter_intersection_ext iter = xmalloc(sizeof (*iter));

  iter->vtable = &xmap_iterator_intersection_ext_vtable;
  iter->msg = xmap_intersection_ext->msg;
  iter->msgs_left = xmap_intersection_ext->n_in - 1;

  return (Xt_xmap_iter)iter;
}

static Xt_xmap_iter xmap_intersection_ext_get_out_iterator(Xt_xmap xmap) {

  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);

  if (xmap_intersection_ext->n_out == 0)
    return NULL;

  Xt_xmap_iter_intersection_ext iter = xmalloc(sizeof (*iter));

  iter->vtable = &xmap_iterator_intersection_ext_vtable;
  iter->msg = xmap_intersection_ext->msg + xmap_intersection_ext->n_in;
  iter->msgs_left = xmap_intersection_ext->n_out - 1;

  return (Xt_xmap_iter)iter;
}

static inline Xt_xmap_iter_intersection_ext
xmiei(void *iter)
{
  return (Xt_xmap_iter_intersection_ext)iter;
}

static int xmap_intersection_ext_iterator_next(Xt_xmap_iter iter) {

  Xt_xmap_iter_intersection_ext iter_intersection = xmiei(iter);

  if (iter_intersection == NULL || iter_intersection->msgs_left == 0)
    return 0;

  iter_intersection->msg++;
  iter_intersection->msgs_left--;

  return 1;
}

static int xmap_intersection_ext_iterator_get_rank(Xt_xmap_iter iter) {

  assert(iter != NULL);
  return xmiei(iter)->msg->rank;
}

static int const *
xmap_intersection_ext_iterator_get_transfer_pos(Xt_xmap_iter iter) {

  assert(iter != NULL);
  struct exchange_ext *restrict msg = xmiei(iter)->msg;
  if ((!msg->num_transfer_pos) | (msg->transfer_pos != NULL)) { } else {
    size_t num_transfer_pos = (size_t)msg->num_transfer_pos;
    int * transfer_pos =
      (msg->transfer_pos = xmalloc(num_transfer_pos * sizeof(*transfer_pos)));
    generate_pos(
      (size_t)msg->num_transfer_pos_ext, msg->transfer_pos_ext,
      num_transfer_pos, transfer_pos);
  }
  return msg->transfer_pos;
}

static int
xmap_intersection_ext_iterator_get_num_transfer_pos(Xt_xmap_iter iter) {
  assert(iter != NULL);
  return xmiei(iter)->msg->num_transfer_pos;
}

static const struct Xt_pos_ext *
xmap_intersection_ext_iterator_get_transfer_pos_ext(Xt_xmap_iter iter) {
  assert(iter != NULL);
  return xmiei(iter)->msg->transfer_pos_ext;
}

static int
xmap_intersection_ext_iterator_get_num_transfer_pos_ext(Xt_xmap_iter iter) {
  assert(iter != NULL);
  return xmiei(iter)->msg->num_transfer_pos_ext;
}

static void xmap_intersection_ext_iterator_delete(Xt_xmap_iter iter) {

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
