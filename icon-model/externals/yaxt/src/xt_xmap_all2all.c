/**
 * @file xt_xmap_all2all.c
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
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_xmap_all2all.h"
#include "xt/xt_xmap_intersection.h"
#include "xt_config_internal.h"
#include "xt_idxlist_internal.h"
#include "instr.h"

static void exchange_idxlists(struct Xt_com_list **src_intersections,
                              size_t *num_src_intersections,
                              struct Xt_com_list **dst_intersections,
                              size_t *num_dst_intersections,
                              int * stripify,
                              Xt_idxlist src_idxlist_local,
                              Xt_idxlist dst_idxlist_local,
                              MPI_Comm comm, Xt_config config)
{

  /*
    Note: The meaning of source (src) and destination (dst) points can already be understood by
    looking at the serial case, where it is just a transformation of sequences of integers
    (called indices). The starting state (source sequence) is transformed into an end state
    (dst sequence). The transformation does not have to be bijective. For each position dpos of
    the dst sequence we have at least one position spos in the src sequence with:
    dst(dpos) = src(spos)
   */

  int comm_size, rank, is_inter;
  xt_mpi_call(MPI_Comm_rank(comm, &rank), comm);
  xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
  int (*get_comm_size)(MPI_Comm, int *)
    = is_inter ? MPI_Comm_remote_size : MPI_Comm_size;
  xt_mpi_call(get_comm_size(comm, &comm_size), comm);

  // allocate memory for intersections
  struct Xt_com_list *dsti = xmalloc((size_t)comm_size * sizeof (*dsti));
  struct Xt_com_list *srci = xmalloc((size_t)comm_size * sizeof (*srci));

  // compute size of local index lists
  size_t src_pack_size = xt_idxlist_get_pack_size(src_idxlist_local, comm);
  size_t dst_pack_size = xt_idxlist_get_pack_size(dst_idxlist_local, comm);
  size_t size_sum = src_pack_size + dst_pack_size;

  if (size_sum >= INT_MAX || size_sum < src_pack_size
      || size_sum < dst_pack_size)
    die("local src+dst index lists are too large");

  int send_buffer_size = (int)size_sum;

  // exchange buffer sizes
  int *restrict pack_sizes
    = xmalloc((size_t)comm_size * sizeof(*pack_sizes) * 2);

  xt_mpi_call(MPI_Allgather(&send_buffer_size, 1, MPI_INT,
                            pack_sizes, 1, MPI_INT, comm), comm);

  int *restrict displ = pack_sizes + comm_size;
  unsigned recv_buffer_size = 0, size_overflow = 0;
  for (size_t i = 0; i < (size_t)comm_size; ++i) {
    displ[i] = (int)recv_buffer_size;
    recv_buffer_size += (unsigned)pack_sizes[i];
    size_overflow |= recv_buffer_size & (1U << (sizeof(int) * CHAR_BIT - 1));
  }
  if (size_overflow)
    die("accumulated buffer sizes too big,"
        " use distributed directory (xt_xmap_dist_dir_new)!");
  void *recv_buffer = xmalloc((size_t)recv_buffer_size + size_sum),
    *send_buffer = (unsigned char *)recv_buffer + (size_t)recv_buffer_size;
  // pack local index lists
  {
    int position = 0;
    xt_idxlist_pack(src_idxlist_local, send_buffer, send_buffer_size,
                    &position, comm);
    xt_idxlist_pack(dst_idxlist_local, send_buffer, send_buffer_size,
                    &position, comm);
  }
  // exchange buffers
  xt_mpi_call(MPI_Allgatherv(send_buffer, send_buffer_size, MPI_PACKED,
                             recv_buffer, pack_sizes, displ, MPI_PACKED,
                             comm), comm);

  size_t dst_isect_count = 0, src_isect_count = 0;
  int large_list_seen = 0;
  // compute intersections
  for (int i = 0; i < comm_size; ++i) {

    int position = 0;
    // unpack buffers unless local
    Xt_idxlist src, dst;
    if (is_inter || i != rank) {
      src = xt_idxlist_unpack((unsigned char *)recv_buffer + displ[i],
                              pack_sizes[i], &position, comm);
      dst = xt_idxlist_unpack((unsigned char *)recv_buffer + displ[i],
                              pack_sizes[i], &position, comm);
    } else {
      src = src_idxlist_local;
      dst = dst_idxlist_local;
    }
    large_list_seen
      |= (xt_idxlist_get_num_indices(src) >= config->idxv_cnv_size)
      | (xt_idxlist_get_num_indices(dst) >= config->idxv_cnv_size);
    Xt_idxlist intersect
      = xt_idxlist_get_intersection_custom(src, dst_idxlist_local, config);
    if (xt_idxlist_get_num_indices(intersect) > 0) {

      dsti[dst_isect_count].list = intersect;
      dsti[dst_isect_count].rank = i;
      ++dst_isect_count;
    }
    else
      xt_idxlist_delete(intersect);

    intersect
      = xt_idxlist_get_intersection_custom(src_idxlist_local, dst, config);
    if (xt_idxlist_get_num_indices(intersect) > 0) {

      srci[src_isect_count].list = intersect;
      srci[src_isect_count].rank = i;
      ++src_isect_count;
    }
    else
      xt_idxlist_delete(intersect);
    if (is_inter || i != rank) {
      xt_idxlist_delete(src);
      xt_idxlist_delete(dst);
    }
  }

  *stripify = large_list_seen;

  free(recv_buffer);
  free(pack_sizes);

  /* minimize memory use of tables */
  *num_src_intersections = src_isect_count;
  srci = xrealloc(srci, src_isect_count * sizeof (**src_intersections));
  *src_intersections = srci;

  *num_dst_intersections = dst_isect_count;
  dsti = xrealloc(dsti, dst_isect_count * sizeof (**dst_intersections));
  *dst_intersections = dsti;
}

Xt_xmap
xt_xmap_all2all_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                    MPI_Comm comm)
{
  return xt_xmap_all2all_custom_new(src_idxlist, dst_idxlist,
                                    comm, &xt_default_config);
}

Xt_xmap
xt_xmap_all2all_custom_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                           MPI_Comm comm, Xt_config config) {
  INSTR_DEF(t_xt_xmap_all2all_new,"xt_xmap_all2all_new")
  INSTR_START(t_xt_xmap_all2all_new);

  // ensure that yaxt is initialized
  assert(xt_initialized());

  int tag_offset;
  MPI_Comm newcomm = xt_mpi_comm_smart_dup(comm, &tag_offset);

  struct Xt_com_list * src_intersections = NULL, * dst_intersections = NULL;
  size_t num_src_intersections, num_dst_intersections;

  int stripify;
  // exchange index lists between all processes in comm
  exchange_idxlists(&src_intersections, &num_src_intersections,
                    &dst_intersections, &num_dst_intersections,
                    &stripify, src_idxlist, dst_idxlist, newcomm,
                    config);

  Xt_xmap (*xmap_new)(int num_src_intersections,
                      const struct Xt_com_list *src_com,
                      int num_dst_intersections,
                      const struct Xt_com_list *dst_com,
                      Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                      MPI_Comm comm)
    = stripify ? xt_xmap_intersection_ext_new : xt_xmap_intersection_new;

  Xt_xmap xmap = xmap_new((int)num_src_intersections, src_intersections,
                          (int)num_dst_intersections, dst_intersections,
                          src_idxlist, dst_idxlist, newcomm);

  xt_mpi_comm_smart_dedup(&newcomm, tag_offset);
  for (size_t i = 0; i < num_src_intersections; ++i)
    if (src_intersections[i].list != NULL)
      xt_idxlist_delete(src_intersections[i].list);
  for (size_t i = 0; i < num_dst_intersections; ++i)
    if (dst_intersections[i].list != NULL)
      xt_idxlist_delete(dst_intersections[i].list);
  free(src_intersections);
  free(dst_intersections);
  INSTR_STOP(t_xt_xmap_all2all_new);
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
