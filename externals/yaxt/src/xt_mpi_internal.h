/**
 * @file xt_mpi_internal.h
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
#ifndef XT_MPI_INTERNAL_H
#define XT_MPI_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>
#include <mpi.h>

#include "core/ppm_visibility.h"

/** \example test_mpi_smartdedup.c
 */

enum xt_mpi_tags {
  xt_mpi_tag_exchange_msg,
  xt_mpi_tag_xmap_dist_dir_src_send,
  xt_mpi_tag_xmap_dist_dir_dst_send,
  xt_mpi_tag_xmap_intersection_header_exchange,
  xt_mpi_tag_xmap_intersection_data_exchange,
  xt_mpi_num_tags,
};

PPM_DSO_INTERNAL size_t
xt_disp2ext_count(size_t disp_len, const int *disp);

PPM_DSO_INTERNAL size_t
xt_disp2ext(size_t disp_len, const int *disp,
            struct Xt_offset_ext *restrict v);


PPM_DSO_INTERNAL void
xt_mpi_init(void);
PPM_DSO_INTERNAL void
xt_mpi_finalize(void);

PPM_DSO_INTERNAL MPI_Comm
xt_mpi_comm_smart_dup(MPI_Comm comm, int *tag_offset);

PPM_DSO_INTERNAL void
xt_mpi_comm_smart_dedup(MPI_Comm *comm, int tag_offset);

/**
 * Given an array of MPI requests, call MPI_Test_some and
 *
 * 1. return if no requests are left unfinished
 * 2. sort non-finished requests into remaining leading part of array
 * 3. update count of remaining requests
 *
 * @param[in,out] num_req pointer to count of requests
 * @param[in,out] req array of requests to test (size *num_req)
 * @param[out] ops_completed array of size at least matching req, this
 * is used as temporary scratch space and overwritten
 * @param[in] comm communicator to use for failure notifications
 */
PPM_DSO_INTERNAL bool
xt_mpi_test_some(int *restrict num_req,
                 MPI_Request req[],
                 int ops_completed[], MPI_Comm comm);

#ifdef _OPENMP
/**
 * multi-thread version, meant to be called by all OpenMP threads,
 * contains omp barrier
 */
PPM_DSO_INTERNAL bool
xt_mpi_test_some_mt(int *restrict num_req,
                    MPI_Request *restrict req,
                    int *restrict ops_completed, MPI_Comm comm);
#endif
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
