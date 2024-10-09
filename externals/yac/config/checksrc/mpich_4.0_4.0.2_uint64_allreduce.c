// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <yaxt.h>

// acx_mpi_job_count=2

int main(int argc, char **argv) {

  // init mpi
  MPI_Init(&argc, &argv);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  if (comm_size != 2) {
    fputs("wrong number of processes\n", stderr);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  {
    uint64_t send_data[2] = {UINT64_MAX, UINT64_MAX};
    send_data[comm_rank] = 0;

    uint64_t recv_data[2] = {13, 13};

    MPI_Allreduce(send_data, recv_data, 2, MPI_UINT64_T, MPI_MIN,
                  MPI_COMM_WORLD);

    if ((recv_data[0] != 0) || (recv_data[1] != 0)) {
      fprintf(stderr,
              "ERROR in MPI_Allreduce (MPI_UINT64_T): "
              "recv_data = [%" PRIu64 "; %" PRIu64 "] (expected [0,0])\n",
              recv_data[0], recv_data[1]);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  {
    int send_data[2] = {INT_MAX, INT_MAX};
    send_data[comm_rank] = 0;

    int recv_data[2] = {13, 13};

    MPI_Allreduce(send_data, recv_data, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if ((recv_data[0] != 0) || (recv_data[1] != 0)) {
      fprintf(stderr,
              "ERROR in MPI_Allreduce (MPI_INT): "
              "recv_data = [%i; %i] (expected [0,0])\n",
              recv_data[0], recv_data[1]);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  {
    Xt_int send_data[2] = {XT_INT_MAX, XT_INT_MAX};
    send_data[comm_rank] = 0;

    Xt_int recv_data[2] = {13, 13};

    MPI_Allreduce(send_data, recv_data, 2, Xt_int_dt, MPI_MIN, MPI_COMM_WORLD);

    if ((recv_data[0] != 0) || (recv_data[1] != 0)) {
      fprintf(stderr,
              "ERROR in MPI_Allreduce (Xt_int_dt): "
              "recv_data = [%" XT_INT_FMT "; %" XT_INT_FMT
              "] (expected [0,0])\n",
              recv_data[0], recv_data[1]);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  // finalise mpi
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return EXIT_SUCCESS;
}
