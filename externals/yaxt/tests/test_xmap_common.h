/**
 * @file test_xmap_common.h
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
#ifndef TEST_XMAP_COMMON_H
#define TEST_XMAP_COMMON_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>

#include <mpi.h>
#include <yaxt.h>

typedef Xt_xmap (*xmap_constructor)(Xt_idxlist src_idxlist,
                                    Xt_idxlist dst_idxlist,
                                    MPI_Comm comm);

enum test_idxlist_size {
  SMALL,
  BIG,
};

int
xt_xmap_self_test_main(xmap_constructor xmap_constructor);

/* note: calls xt_idxlist_delete for src_idxlist and dst_idxlist */
void
test_self_xmap_construct(Xt_idxlist src_idxlist,
                         Xt_idxlist dst_idxlist,
                         xmap_constructor new_xmap, MPI_Comm comm);
void
test_self_xmap_construct_idxvec(
  const Xt_int *src_index_list, int num_src_indices,
  const Xt_int *dst_index_list, int num_dst_indices,
  xmap_constructor new_xmap, MPI_Comm comm);

void
test_self_xmap_construct_idxstripes(
  const struct Xt_stripe *src_index_list, int num_src_stripes,
  const struct Xt_stripe *dst_index_list, int num_dst_stripes,
  xmap_constructor new_xmap, MPI_Comm comm);

int
xt_xmap_parallel_test_main(xmap_constructor xmap_constructor);

extern MPI_Comm xt_intra_group_comm;

int
xt_xmap_intercomm_parallel_test_main(xmap_constructor xmap_constructor,
                                     bool call_initialize, bool call_finalize);

void
test_ping_pong(xmap_constructor xmap_new, MPI_Comm comm,
               int ping_rank, int pong_rank);

void
check_xmap_allgather_analog_xmap(Xt_xmap xmap, MPI_Comm comm);

void
test_ring_1d(xmap_constructor xmap_new, MPI_Comm comm);

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
