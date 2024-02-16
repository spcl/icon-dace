/**
 * @file mpich_3.4.3_mpi_waitall.c
 *
 * @copyright Copyright  (C)  2023 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Nils-Arne Dreier <dreier@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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
 * acx_mpi_job_count=1
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main (void) {

  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  int const count = 64800;
  MPI_Datatype dt;
  MPI_Type_contiguous(count, MPI_DOUBLE, &dt);
  MPI_Type_commit(&dt);

  double * send_data = calloc(2 * count, sizeof(send_data));
  double * recv_data = send_data + count;
  MPI_Request requests[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
  MPI_Isend(send_data, 1, dt, rank, 0, comm, &requests[0]);
  MPI_Irecv(recv_data, 1, dt, rank, 0, comm, &requests[1]);
#define MPI_BUG
#ifdef MPI_BUG
  MPI_Comm_free(&comm);
  MPI_Type_free(&dt);
  MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
#else
  MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
  MPI_Comm_free(&comm);
  MPI_Type_free(&dt);
#endif

  MPI_Finalize();

  return EXIT_SUCCESS;
}
