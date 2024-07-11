/**
 * @file xt_redist_p2p_ext.c
 *
 * @copyright Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
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

#define TOKEN_PASTE2_(a,b) a##b
#define TOKEN_PASTE2(a,b) TOKEN_PASTE2_(a,b)
#define TOKEN_PASTE3_(a,b,c) a##b##c
#define TOKEN_PASTE3(a,b,c) TOKEN_PASTE3_(a,b,c)


#define XT_GENERATE_EXT_DATATYPE \
  TOKEN_PASTE3(generate_,XT_EXT_TAG,_datatype)
#define XT_GENERATE_EXT_MSG_INFOS \
  TOKEN_PASTE3(generate_,XT_EXT_TAG,_msg_infos)
#define XT_REDIST_P2P_EXT_NEW \
  TOKEN_PASTE3(xt_redist_p2p_,XT_EXT_TAG,_new)
#define XT_REDIST_P2P_EXT_CUSTOM_NEW \
  TOKEN_PASTE3(xt_redist_p2p_,XT_EXT_TAG,_custom_new)

static MPI_Datatype
XT_GENERATE_EXT_DATATYPE(int num_transfer_pos_ext,
                      const struct Xt_pos_ext transfer_pos_ext[],
                      int num_ext, const XT_EXT_TYPE extents[],
                      const int psum_ext_size[],
                      void **work_buf, size_t *work_buf_size,
                      MPI_Datatype base_datatype, MPI_Comm comm)
{
  if (num_transfer_pos_ext > 0)
  {
    XT_EXT_TYPE *dt_stripes;
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
                 num_ext, psum_ext_size);
    do
    {
      struct Xt_pos_ext current_pos_ext = transfer_pos_ext[i];
      if (num_dt_stripes >= size_dt_stripes) {
      more_stripes:
        size_dt_stripes *= 2;
        dt_stripes = xrealloc(dt_stripes,
                              size_dt_stripes * sizeof (*dt_stripes));
      }
      do {
        /* find extent containing start position of current range */
        search_start_ext = pos2disp2(current_pos_ext.start,
                                     num_ext, psum_ext_size,
                                     search_start_ext);
        XT_EXT_TYPE base_ext = extents[search_start_ext],
          derived_ext;
        int pos_remaining = current_pos_ext.start
          - psum_ext_size[search_start_ext];
        derived_ext.start
          = base_ext.start + pos_remaining * base_ext.stride;
        int isign_mask_current_pos_ext_size = isign_mask(current_pos_ext.size);
        XT_EXT_STRIDE_MASK_PREP;
        /* find number of positions in containing extent,
         * which precede current_pos_ext.start
         * if (current_pos_ext.size < 0)
         * or follow current_pos_ext.start
         * if (current_pos_ext.size > 0) */
        derived_ext.size = imin(abs(current_pos_ext.size),
                                (~isign_mask_current_pos_ext_size
                                 & (base_ext.size - pos_remaining))
                                | (isign_mask_current_pos_ext_size
                                   & (pos_remaining + 1)));
        derived_ext.stride
          = (~XT_EXT_STRIDE_MASK & base_ext.stride)
          | (XT_EXT_STRIDE_MASK & -base_ext.stride);
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
      = XT_MPI_GENERATE_DATATYPE(dt_stripes, (int)num_dt_stripes,
                                 base_datatype, comm);
    *work_buf = dt_stripes;
    *work_buf_size = size_dt_stripes * sizeof (*dt_stripes);
    return type;
  }
  else
    return MPI_DATATYPE_NULL;
}

static void
XT_GENERATE_EXT_MSG_INFOS(int num_msgs, Xt_xmap_iter iter,
                       int num_ext,
                       const XT_EXT_TYPE extents[],
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
        = XT_GENERATE_EXT_DATATYPE(curr_num_transfer_pos_ext,
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
XT_REDIST_P2P_EXT_NEW(Xt_xmap xmap,
                      int num_src_ext,
                      const XT_EXT_TYPE src_extents[],
                      int num_dst_ext,
                      const XT_EXT_TYPE dst_extents[],
                      MPI_Datatype datatype)
{
  return XT_REDIST_P2P_EXT_CUSTOM_NEW(xmap, num_src_ext, src_extents,
                                      num_dst_ext, dst_extents, datatype,
                                      (Xt_config)&xt_default_config);
}

Xt_redist
XT_REDIST_P2P_EXT_CUSTOM_NEW(Xt_xmap xmap,
                             int num_src_ext,
                             const XT_EXT_TYPE src_extents[],
                             int num_dst_ext,
                             const XT_EXT_TYPE dst_extents[],
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
  XT_GENERATE_EXT_MSG_INFOS(nrecv, dst_iter, num_dst_ext, dst_extents,
                            datatype, recv_msgs, comm);
  if (dst_iter) xt_xmap_iterator_delete(dst_iter);

  Xt_xmap_iter src_iter = xt_xmap_get_out_iterator(xmap);
  XT_GENERATE_EXT_MSG_INFOS(nsend, src_iter, num_src_ext, src_extents,
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

#undef TOKEN_PASTE2
#undef TOKEN_PASTE2_
#undef TOKEN_PASTE3
#undef TOKEN_PASTE3_

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
