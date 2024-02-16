/**
 * @file xt_redist_p2p.c
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
#include <assert.h>

#include <mpi.h>

#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_p2p.h"
#include "xt_redist_internal.h"
#include "xt/xt_redist_single_array_base.h"
#include "xt/xt_xmap.h"
#include "xt/xt_idxlist.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "xt_config_internal.h"

#include "xt_arithmetic_util.h"

static MPI_Datatype
generate_datatype(const int *transfer_pos, int num_transfer_pos,
                  const int *offsets, MPI_Datatype base_datatype, MPI_Comm comm)
{
  int const * displ;
  int * tmp_displ = NULL;

  if (offsets != NULL) {

    tmp_displ = xmalloc((size_t)num_transfer_pos * sizeof(int));

    for (int i = 0; i < num_transfer_pos; ++i)
      tmp_displ[i] = offsets[transfer_pos[i]];

    displ = tmp_displ;

  } else
    displ = transfer_pos;


  MPI_Datatype type
    = xt_mpi_generate_datatype(displ, num_transfer_pos, base_datatype, comm);

  free(tmp_displ);

  return type;
}

static void
generate_msg_infos(int num_msgs, Xt_xmap_iter iter, const int *offsets,
                   MPI_Datatype base_datatype, struct Xt_redist_msg *msgs,
                   MPI_Comm comm) {

  if (num_msgs > 0) {
    struct Xt_redist_msg *restrict curr_msg = msgs;

    do {

      const int *curr_transfer_pos = xt_xmap_iterator_get_transfer_pos(iter);
      int curr_num_transfer_pos = xt_xmap_iterator_get_num_transfer_pos(iter);

      curr_msg->datatype
        = generate_datatype(curr_transfer_pos, curr_num_transfer_pos,
                            offsets, base_datatype, comm);
      curr_msg->rank = xt_xmap_iterator_get_rank(iter);

      curr_msg++;
    } while (xt_xmap_iterator_next(iter));
  }
}

Xt_redist
xt_redist_p2p_off_new(Xt_xmap xmap, const int *src_offsets,
                      const int *dst_offsets, MPI_Datatype datatype)
{
  return xt_redist_p2p_off_custom_new(xmap, src_offsets, dst_offsets, datatype,
                                      (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_p2p_off_custom_new(Xt_xmap xmap, const int *src_offsets,
                             const int *dst_offsets, MPI_Datatype datatype,
                             Xt_config config) {
  // ensure that yaxt is initialized
  assert(xt_initialized());

  int nsend = xt_xmap_get_num_destinations(xmap),
    nrecv = xt_xmap_get_num_sources(xmap);
  size_t nmsg = (size_t)nsend + (size_t)nrecv;
  struct Xt_redist_msg *msgs = xmalloc(nmsg * sizeof (*msgs)),
    *send_msgs = msgs, *recv_msgs = msgs + nsend;
  int tag_offset;
  MPI_Comm comm
    = xt_mpi_comm_smart_dup(xt_xmap_get_communicator(xmap), &tag_offset);

  Xt_xmap_iter dst_iter = xt_xmap_get_in_iterator(xmap);
  generate_msg_infos(nrecv, dst_iter, dst_offsets, datatype, recv_msgs,
                     comm);
  if (dst_iter) xt_xmap_iterator_delete(dst_iter);

  Xt_xmap_iter src_iter = xt_xmap_get_out_iterator(xmap);
  generate_msg_infos(nsend, src_iter, src_offsets, datatype, send_msgs,
                     comm);
  if (src_iter) xt_xmap_iterator_delete(src_iter);

  struct Xt_config_ config_ = *config;
  config_.flags |= exch_no_dt_dup;

  Xt_redist result = xt_redist_single_array_base_custom_new(
    nsend, nrecv, send_msgs, recv_msgs, comm, &config_);

  free(msgs);
  xt_mpi_comm_smart_dedup(&comm, tag_offset);
  return result;
}

/* ====================================================================== */

struct ext_disp
{
  int disp, ext_idx;
};

static inline struct ext_disp
pos2disp(int pos, int num_ext, const struct Xt_offset_ext extents[],
         const int psum_ext_size[])
{
  int j = 0;
  /* FIXME: use bsearch if linear search is too slow, i.e. num_ext >> 1000 */
  /* what extent covers the pos'th position? */
  while (j < num_ext && pos >= psum_ext_size[j + 1])
    ++j;
  int disp = extents[j].start + (pos - psum_ext_size[j]) * extents[j].stride;
  return (struct ext_disp){ .disp = disp, .ext_idx = j };
}

static inline struct ext_disp
pos2disp2(int pos, int num_ext, const struct Xt_offset_ext extents[],
          const int psum_ext_size[], int start_ext)
{
  int j = start_ext;
  if (pos < psum_ext_size[j + 1] && pos >= psum_ext_size[j])
    ;
  else if (pos < psum_ext_size[j + 1])
  {
    j = 0;
    while (j < start_ext && pos >= psum_ext_size[j + 1])
      ++j;
  }
  else
    while (j < num_ext && pos >= psum_ext_size[j + 1])
      ++j;
  int disp = extents[j].start + (pos - psum_ext_size[j]) * extents[j].stride;
  return (struct ext_disp){ .disp = disp, .ext_idx = j };
}

static MPI_Datatype
generate_ext_datatype(int num_transfer_pos_ext,
                      const struct Xt_pos_ext transfer_pos_ext[],
                      int num_ext, const struct Xt_offset_ext extents[],
                      const int psum_ext_size[],
                      void **work_buf, size_t *work_buf_size,
                      MPI_Datatype base_datatype, MPI_Comm comm)
{
  if (num_transfer_pos_ext > 0)
  {
    struct Xt_offset_ext *dt_stripes;
    size_t size_dt_stripes, num_dt_stripes = 0;
    enum
    {
      dt_stripes_init_size = 8,
      dt_stripes_init_alloc = dt_stripes_init_size * sizeof (*dt_stripes),
    };
    if (*work_buf_size < dt_stripes_init_alloc) {
      dt_stripes = xrealloc(*work_buf, dt_stripes_init_alloc);
      size_dt_stripes = dt_stripes_init_size;
    } else {
      dt_stripes = *work_buf;
      size_dt_stripes = *work_buf_size / sizeof (*dt_stripes);
    }
    int i = 0,
      search_start_ext
      = pos2disp(transfer_pos_ext[0].start,
                 num_ext, extents, psum_ext_size).ext_idx;
    do
    {
      struct Xt_pos_ext current_pos_ext = transfer_pos_ext[i];
      if (num_dt_stripes >= size_dt_stripes)
      {
      more_stripes:
        size_dt_stripes *= 2;
        dt_stripes = xrealloc(dt_stripes,
                              size_dt_stripes * sizeof (*dt_stripes));
      }
      do {
        /* find extent containing start position of current range */
        struct ext_disp pos = pos2disp2(current_pos_ext.start,
                                        num_ext, extents, psum_ext_size,
                                        search_start_ext);
        search_start_ext = pos.ext_idx;
        struct Xt_offset_ext base_ext = extents[pos.ext_idx],
          derived_ext;
        int preceding = psum_ext_size[pos.ext_idx];
        derived_ext.start
          = base_ext.start + ((current_pos_ext.start - preceding)
                              * base_ext.stride);
        int isign_mask_current_pos_ext_size = isign_mask(current_pos_ext.size);
        /* find number of positions in containing extent,
         * which precede current_pos_ext.start
         * if (current_pos_ext.size < 0)
         * or follow current_pos_ext.start
         * if (current_pos_ext.size > 0) */
        derived_ext.size = imin(abs(current_pos_ext.size),
                                (~isign_mask_current_pos_ext_size
                                 & (base_ext.size
                                    - (current_pos_ext.start - preceding)))
                                | (isign_mask_current_pos_ext_size
                                   & (current_pos_ext.start - preceding + 1)));
        derived_ext.stride
          = (~isign_mask_current_pos_ext_size & base_ext.stride)
          | (isign_mask_current_pos_ext_size & -base_ext.stride);
        dt_stripes[num_dt_stripes++] = derived_ext;
        current_pos_ext.size
          += (~isign_mask_current_pos_ext_size & -derived_ext.size)
          | (isign_mask_current_pos_ext_size & derived_ext.size);
        current_pos_ext.start += derived_ext.size;
      } while ((abs(current_pos_ext.size) > 0)
               & (num_dt_stripes < size_dt_stripes));
      /* current_pos_ext hasn't been mapped completely, get more
       * stripe memory */
      if (abs(current_pos_ext.size) > 0)
        goto more_stripes;
      /* only advance current_pos_ext after it has been mapped completely */
    } while (++i < num_transfer_pos_ext);
    MPI_Datatype type
      = xt_mpi_generate_datatype_stripe(dt_stripes, (int)num_dt_stripes,
                                        base_datatype, comm);
    *work_buf = dt_stripes;
    *work_buf_size = size_dt_stripes * sizeof (*dt_stripes);
    return type;
  }
  else
    return MPI_DATATYPE_NULL;
}

static void
generate_ext_msg_infos(int num_msgs, Xt_xmap_iter iter,
                       int num_ext,
                       const struct Xt_offset_ext extents[],
                       MPI_Datatype base_datatype,
                       struct Xt_redist_msg *msgs,
                       MPI_Comm comm)
{
  if (num_msgs > 0) {
    /* partial sums of ext sizes */
    int *restrict psum_ext_size
      = xmalloc(((size_t)num_ext + 1) * sizeof (psum_ext_size[0]));
    int accum = 0;
    for (size_t i = 0; i < (size_t)num_ext; ++i) {
      psum_ext_size[i] = accum;
      accum += extents[i].size;
    }
    psum_ext_size[num_ext] = accum;

    void *buf = NULL;
    size_t buf_size = 0;

    struct Xt_redist_msg *curr_msg = msgs;
    do {

      const struct Xt_pos_ext *curr_transfer_pos_ext
        = xt_xmap_iterator_get_transfer_pos_ext(iter);
      int curr_num_transfer_pos_ext
        = xt_xmap_iterator_get_num_transfer_pos_ext(iter);

      curr_msg->datatype
        = generate_ext_datatype(curr_num_transfer_pos_ext,
                                curr_transfer_pos_ext,
                                num_ext, extents, psum_ext_size,
                                &buf, &buf_size,
                                base_datatype, comm);
      curr_msg->rank = xt_xmap_iterator_get_rank(iter);

      curr_msg++;

    } while (xt_xmap_iterator_next(iter));
    free(psum_ext_size);
    free(buf);
  }
}

Xt_redist
xt_redist_p2p_ext_new(Xt_xmap xmap,
                      int num_src_ext,
                      const struct Xt_offset_ext src_extents[],
                      int num_dst_ext,
                      const struct Xt_offset_ext dst_extents[],
                      MPI_Datatype datatype)
{
  return xt_redist_p2p_ext_custom_new(xmap, num_src_ext, src_extents,
                                      num_dst_ext, dst_extents, datatype,
                                      (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_p2p_ext_custom_new(Xt_xmap xmap,
                             int num_src_ext,
                             const struct Xt_offset_ext src_extents[],
                             int num_dst_ext,
                             const struct Xt_offset_ext dst_extents[],
                             MPI_Datatype datatype,
                             Xt_config config)
{
  // ensure that yaxt is initialized
  assert(xt_initialized());
  int tag_offset;
  MPI_Comm comm
    = xt_mpi_comm_smart_dup(xt_xmap_get_communicator(xmap), &tag_offset);

  int nrecv = xt_xmap_get_num_sources(xmap),
    nsend = xt_xmap_get_num_destinations(xmap);
  size_t nmsg = (size_t)nrecv + (size_t)nsend;
  struct Xt_redist_msg *msgs = xmalloc(nmsg * sizeof (*msgs)),
    *send_msgs = msgs, *recv_msgs = msgs + nsend;
  Xt_xmap_iter dst_iter = xt_xmap_get_in_iterator(xmap);
  generate_ext_msg_infos(nrecv, dst_iter, num_dst_ext, dst_extents,
                         datatype, recv_msgs, comm);
  if (dst_iter) xt_xmap_iterator_delete(dst_iter);

  Xt_xmap_iter src_iter = xt_xmap_get_out_iterator(xmap);
  generate_ext_msg_infos(nsend, src_iter, num_src_ext, src_extents,
                         datatype, send_msgs, comm);
  if (src_iter) xt_xmap_iterator_delete(src_iter);

  struct Xt_config_ config_ = *config;
  config_.flags |= exch_no_dt_dup;

  Xt_redist result = xt_redist_single_array_base_custom_new(
    nsend, nrecv, send_msgs, recv_msgs, comm, &config_);

  free(msgs);
  xt_mpi_comm_smart_dedup(&comm, tag_offset);
  return result;
}

/* ====================================================================== */

static inline void
aux_gen_simple_block_offsets(int block_offsets[], const int block_sizes[],
                             size_t num_blocks) {

  if (num_blocks > 0) {
    int accum = 0;
    for (size_t i = 0; i < num_blocks; ++i) {
      block_offsets[i] = accum;
      accum += block_sizes[i];
    }
  }
}

static MPI_Datatype
generate_block_datatype(const int *transfer_pos, int num_transfer_pos,
                        const int *block_offsets, const int *block_sizes,
                        MPI_Datatype base_datatype, MPI_Comm comm) {

  assert(block_sizes && block_offsets);

  int *bdispl_vec = xmalloc(2 * (size_t)num_transfer_pos * sizeof(*bdispl_vec)),
    *blen_vec = bdispl_vec + num_transfer_pos;

  for (int i = 0; i < num_transfer_pos; ++i) {
    int j = transfer_pos[i];
    bdispl_vec[i] = block_offsets[j];
    blen_vec[i] = block_sizes[j];
  }

  MPI_Datatype type
    = xt_mpi_generate_datatype_block(bdispl_vec, blen_vec,
                                     num_transfer_pos, base_datatype, comm);

  free(bdispl_vec);

  return type;
}

static void
generate_block_msg_infos(int num_msgs, Xt_xmap_iter iter,
                         const int *block_offsets,
                         const int *block_sizes, int **aux_offsets,
                         size_t num_blocks,
                         MPI_Datatype base_datatype,
                         struct Xt_redist_msg *msgs, MPI_Comm comm) {

  if (num_msgs > 0) {

    const int *block_offsets_;
    if (block_offsets)
      block_offsets_ = block_offsets;
    else {
      block_offsets_ = *aux_offsets
        = xrealloc(*aux_offsets, num_blocks * sizeof(*block_offsets_));
      aux_gen_simple_block_offsets(*aux_offsets, block_sizes, num_blocks);
    }

    size_t ofs = 0;
    do {
      const int *curr_transfer_pos = xt_xmap_iterator_get_transfer_pos(iter);
      int curr_num_transfer_pos = xt_xmap_iterator_get_num_transfer_pos(iter);
      msgs[ofs].datatype
        = generate_block_datatype(curr_transfer_pos, curr_num_transfer_pos,
                                  block_offsets_, block_sizes, base_datatype,
                                  comm);
      msgs[ofs].rank = xt_xmap_iterator_get_rank(iter);

      ofs++;
    } while (xt_xmap_iterator_next(iter));

  }
}

Xt_redist
xt_redist_p2p_blocks_off_new(Xt_xmap xmap,
                             const int *src_block_offsets,
                             const int *src_block_sizes,
                             int src_block_num,
                             const int *dst_block_offsets,
                             const int *dst_block_sizes,
                             int dst_block_num,
                             MPI_Datatype datatype)
{
  return xt_redist_p2p_blocks_off_custom_new(
    xmap, src_block_offsets, src_block_sizes, src_block_num, dst_block_offsets,
    dst_block_sizes, dst_block_num, datatype, (Xt_config)&xt_default_config);
}


Xt_redist
xt_redist_p2p_blocks_off_custom_new(Xt_xmap xmap,
                                    const int *src_block_offsets,
                                    const int *src_block_sizes,
                                    int src_block_num,
                                    const int *dst_block_offsets,
                                    const int *dst_block_sizes,
                                    int dst_block_num,
                                    MPI_Datatype datatype,
                                    Xt_config config)
{
  // ensure that yaxt is initialized
  assert(xt_initialized() && src_block_sizes && dst_block_sizes);

  int tag_offset;
  MPI_Comm comm
    = xt_mpi_comm_smart_dup(xt_xmap_get_communicator(xmap), &tag_offset);


  int nsend = xt_xmap_get_num_destinations(xmap),
    nrecv = xt_xmap_get_num_sources(xmap);

  size_t nmsg = ((size_t)nsend + (size_t)nrecv);
  struct Xt_redist_msg *msgs = xmalloc(nmsg * sizeof (*msgs));

  int *aux_offsets = NULL;

  Xt_xmap_iter dst_iter = xt_xmap_get_in_iterator(xmap),
    src_iter = xt_xmap_get_out_iterator(xmap);

  // dst part:
#ifndef NDEBUG
  int max_dst_pos = xt_xmap_get_max_dst_pos(xmap);
  if (dst_block_num < max_dst_pos)
    die("xt_redist_p2p_blocks_off_new: dst_block_num too small");
#endif
  generate_block_msg_infos(nrecv, dst_iter, dst_block_offsets, dst_block_sizes,
                           &aux_offsets, (size_t)dst_block_num,
                           datatype, msgs, comm);
  if (dst_iter) xt_xmap_iterator_delete(dst_iter);

  // src part:
#ifndef NDEBUG
  int max_src_pos = xt_xmap_get_max_src_pos(xmap);
  if (src_block_num < max_src_pos)
    die("xt_redist_p2p_blocks_off_new: src_block_num too small");
#endif
  generate_block_msg_infos(nsend, src_iter, src_block_offsets, src_block_sizes,
                           &aux_offsets, (size_t)src_block_num,
                           datatype, msgs+nrecv, comm);
  free(aux_offsets);

  if (src_iter) xt_xmap_iterator_delete(src_iter);

  struct Xt_config_ config_ = *config;
  config_.flags |= exch_no_dt_dup;

  Xt_redist result
    = xt_redist_single_array_base_custom_new(
      nsend, nrecv, msgs+nrecv, msgs, comm, &config_);

  free(msgs);
  xt_mpi_comm_smart_dedup(&comm, tag_offset);
  return result;
}

Xt_redist xt_redist_p2p_blocks_new(Xt_xmap xmap,
                                   const int *src_block_sizes,
                                   int src_block_num,
                                   const int *dst_block_sizes,
                                   int dst_block_num,
                                   MPI_Datatype datatype)
{
  return xt_redist_p2p_blocks_custom_new(
    xmap, src_block_sizes, src_block_num, dst_block_sizes, dst_block_num,
    datatype, (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_p2p_blocks_custom_new(Xt_xmap xmap,
                                const int *src_block_sizes, int src_block_num,
                                const int *dst_block_sizes, int dst_block_num,
                                MPI_Datatype datatype,
                                Xt_config config)
{
  return xt_redist_p2p_blocks_off_custom_new(
    xmap, NULL, src_block_sizes, src_block_num,
    NULL, dst_block_sizes, dst_block_num, datatype, config);
}


Xt_redist xt_redist_p2p_new(Xt_xmap xmap, MPI_Datatype datatype)
{
  return xt_redist_p2p_custom_new(xmap, datatype,
                                  (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_p2p_custom_new(Xt_xmap xmap, MPI_Datatype datatype, Xt_config config)
{
  return xt_redist_p2p_off_custom_new(xmap, NULL, NULL, datatype, config);
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
