/**
 * @file rrobin.c
 *
 * @copyright Copyright  (C)  2012 Anna Fuchs   <anna.fuchsmakaeva@googlemail.com>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Anna Fuchs   <anna.fuchsmakaeva@googlemail.com>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Anna Fuchs   <anna.fuchsmakaeva@googlemail.com>
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

#include <yaxt.h>

int main(void) {

  //init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize (MPI_COMM_WORLD);

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  {
    // init source and destination data array, local data 5 elements length
    enum { len = 5 };
    int src_array[len];
    int dst_array[len];

    for (int i = 0; i < len; i++) {
      // first position = 1, second position = rank,
      // third position element_index; for rank < 10, len < 10
      src_array[i] = 100 + rank*10 + i;
      dst_array[i] = -1;

      //print source
      printf("SOURCEvalue: %d, element_index: %d, rank: %d \n", src_array[i],
             i, rank);
    }

    // output barrier
    MPI_Barrier (MPI_COMM_WORLD);


    struct Xt_stripe src_stripes = {.start = (Xt_int)(rank*len), .stride = 1,  .nstrides = len};
    // source index list by stripes
    Xt_idxlist src_idxlist = xt_idxstripes_new(&src_stripes, 1);

    struct Xt_stripe dst_stripes = {.start = (Xt_int)(((rank+1)*len)%(size*len)),
                                    .stride = 1,  .nstrides = len};
    // destination index list by stripes
    Xt_idxlist dst_idxlist = xt_idxstripes_new(&dst_stripes, 1);

    // xmap
    Xt_xmap xmap
      = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_INTEGER);

    // array pointer, especially necessary for data array number > 1
    int* src_array_p = &src_array[0];
    int* dst_array_p = &dst_array[0];

    //Exchange
    xt_redist_s_exchange1(redist, src_array_p, dst_array_p);

    // print result
    // for emelent_index < 10, len < 10:
    // first and last positions: should stay;
    // second position: should be +1, except last rank should have 0
    for (int p = 0; p < len; p++)
      printf("DESTvalue: %d, element_index: %d, rank: %d \n",
             dst_array[p], p, rank);


    xt_redist_delete(redist);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
  }

  MPI_Finalize();

  return 0;
}
/*
 * Local Variables:
 * coding: utf-8
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
 * license-default: "bsd"
 * End:
 */
