/**
 * @file row2col.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

enum {
  global_size_x = 1024,
  global_size_y = 512,
};

int main(void) {

  // init

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize (MPI_COMM_WORLD);

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  int part_size_x = (global_size_x + size - 1) / size;
  int part_size_y = (global_size_y + size - 1) / size;
  int local_part_size_x = MIN(part_size_x, global_size_x - part_size_x);
  int local_part_size_y = MIN(part_size_y, global_size_y - part_size_y);

  int src_array_size[2] = {local_part_size_y, global_size_x};
  int tgt_array_size[2] = {global_size_y, local_part_size_x};
  Xt_int global_size[2] = {(Xt_int)global_size_y, (Xt_int)global_size_x};
  Xt_int src_local_start[2] = {(Xt_int)(rank * part_size_y), 0};
  Xt_int tgt_local_start[2] = {0, (Xt_int)(rank * part_size_x)};

  // define decomposition

  Xt_idxlist src_idxlist = xt_idxsection_new(0, 2, global_size, src_array_size,
                                             src_local_start);
  Xt_idxlist tgt_idxlist = xt_idxsection_new(0, 2, global_size, tgt_array_size,
                                             tgt_local_start);

  // generate exchange map

  Xt_xmap xmap = xt_xmap_all2all_new(src_idxlist, tgt_idxlist, MPI_COMM_WORLD);

  // generate redistribution object

  Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

  // do the exchange

  double * src_array = malloc((size_t)src_array_size[0]
                              * (size_t)src_array_size[1] * sizeof(*src_array));
  double * tgt_array = malloc((size_t)tgt_array_size[0]
                              * (size_t)tgt_array_size[1] * sizeof(*tgt_array));

  for (Xt_int i = 0; i < src_array_size[0] * src_array_size[1]; ++i)
    src_array[i] = (double)rank;
  for (Xt_int i = 0; i < tgt_array_size[0] * tgt_array_size[1]; ++i)
    tgt_array[i] = (double)rank;

  xt_redist_s_exchange1(redist, src_array, tgt_array);

  // clean up

  free(tgt_array);
  free(src_array);
  xt_redist_delete(redist);
  xt_xmap_delete(xmap);
  xt_idxlist_delete(tgt_idxlist);
  xt_idxlist_delete(src_idxlist);

  // finalise

  xt_finalize();

  MPI_Finalize();

  return 0;
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * coding: utf-8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
 * license-default: "bsd"
 * End:
 */
