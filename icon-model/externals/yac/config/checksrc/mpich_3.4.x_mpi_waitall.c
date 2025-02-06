// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

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
