/**
 * @file xt_xmap_intersection.h
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

#ifndef XT_XMAP_INTERSECTION_H
#define XT_XMAP_INTERSECTION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xt/xt_core.h"

/** \example test_xmap_intersection_parallel.c
 */

struct Xt_com_list {
  Xt_idxlist list;
  int rank;
};

struct Xt_com_pos {
  // list of relative positions in memory to send or receive
  int * transfer_pos;
  int num_transfer_pos;
  int rank;
};

/**
 * constructor for an exchange map \n
 * this operation is collective over all processes in comm \n
 * it uses the provided intersection information to generate the exchange map
 *
 * @param[in] num_src_intersections number of source intersections
 * @param[in] src_com               array containing the source intersections
 *                                  with the other processes and their rank
 * @param[in] num_dst_intersections number of destination intersections
 * @param[in] dst_com               array containing the destination
 *                                  intersections with the other
 *                                  processes and their rank
 * @param[in] src_idxlist           source index list
 * @param[in] dst_idxlist           destination index list
 * @param[in] comm                  MPI communicator that contains all processes
 *                                  that part in the exchange
 */
Xt_xmap
xt_xmap_intersection_new(int num_src_intersections,
                         const struct Xt_com_list
                         src_com[num_src_intersections],
                         int num_dst_intersections,
                         const struct Xt_com_list
                         dst_com[num_dst_intersections],
                         Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                         MPI_Comm comm);

/**
 * constructor for an exchange map \n
 * This operation is collective over all processes in comm \n
 * it uses the provided intersection information to generate the
 * exchange map.
 * Internally this function uses ranges to represent index list
 * positions and is therefore conserving space for somewhat
 * contiguous index lists. Depending on the size and shape of
 * intersections this can impact performance drastically..
 *
 * @param[in] num_src_intersections number of source intersections
 * @param[in] src_com               array containing the source intersections
 *                                  with the other processes and their rank
 * @param[in] num_dst_intersections number of destination intersections
 * @param[in] dst_com               array containing the destination
 *                                  intersections with the other
 *                                  processes and their rank
 * @param[in] src_idxlist           source index list
 * @param[in] dst_idxlist           destination index list
 * @param[in] comm                  MPI communicator that contains all processes
 *                                  that part in the exchange
 */
Xt_xmap
xt_xmap_intersection_ext_new(
  int num_src_intersections,
  const struct Xt_com_list src_com[num_src_intersections],
  int num_dst_intersections,
  const struct Xt_com_list dst_com[num_dst_intersections],
  Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
  MPI_Comm comm);

/**
 * constructor for an exchange map \n
 * this operation is collective over all processes in comm \n
 * it uses the provided intersection information to generate the exchange map
 *
 * @param[in] num_src_msg number of source messages
 * @param[in] src_com     array containing relative positions for all
 *                        source messages and the destination rank
 * @param[in] num_dst_msg number of destination messages
 * @param[in] dst_com     array containing relative positions for all
 *                        destination messages and the source rank
 * @param[in] comm        MPI communicator that contains all processes
 *                        that part in the exchange
 */
Xt_xmap
xt_xmap_intersection_pos_new(
  int num_src_msg, const struct Xt_com_pos src_com[num_src_msg],
  int num_dst_msg, const struct Xt_com_pos dst_com[num_dst_msg],
  MPI_Comm comm);

#endif // XT_XMAP_INTERSECTION_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
