/**
 * @file mpich_zero_size_contiguous.c
 * @brief test for a defect in MPICH releases 3.4.3, 4.0, 4.0.1 and 4.0.2
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

/* run-time settings/expectations for configure-time checks
 * acx_mpirun_num_tasks=1
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  // contiguous data type with count zero
  MPI_Datatype cont_zero_ddt;
  MPI_Type_contiguous(0, MPI_CHAR, &cont_zero_ddt);
  MPI_Type_commit(&cont_zero_ddt);

  double inbuf[1];
  double outbuf[1];

  int pack_size;
  MPI_Pack_size(1, cont_zero_ddt, MPI_COMM_WORLD, &pack_size);

  void *pack_buf = malloc(pack_size);
  if (!pack_buf) {
    fprintf(stderr, "error: failed to allocate buffer of size %d!\n",
            pack_size);
    exit(1);
  }

  int position = 0;
  inbuf[0] = 3.1415;
  outbuf[0] = 23;

  MPI_Pack(
    inbuf, 1, cont_zero_ddt, pack_buf, pack_size, &position,
    MPI_COMM_WORLD);
  position = 0;
  MPI_Unpack(
    pack_buf, 1, &position, outbuf, 1,
    cont_zero_ddt, MPI_COMM_WORLD);
  /* this test usually fails in one of the above two functions with a
   * segfault or internal error but just to make sure outbuf wasn't
   * modified inspect its contents here */
  if (outbuf[0] < 23 || outbuf[0] > 23) {
    fputs("Unexpected modification of output buffer!\n", stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  free(pack_buf);

  MPI_Type_free(&cont_zero_ddt);

  MPI_Finalize();

  return EXIT_SUCCESS;
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
 * license-default: "bsd"
 * license-markup: "doxygen"
 * End:
 */
