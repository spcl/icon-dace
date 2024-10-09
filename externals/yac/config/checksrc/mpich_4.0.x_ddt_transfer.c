// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

/* run-time settings/expectations for configure-time checks
 * acx_mpi_job_count=2
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int size, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  enum {
    DDT_COUNT = 1024,
    EXCH_COUNT = 128,
    BASE_COUNT = 4,
    RANK_COUNT = 2};

  if (size != RANK_COUNT) {
    fprintf(stderr, "wrong number of processes (has to be %d)", RANK_COUNT);
    exit(EXIT_FAILURE);
  }

  MPI_Datatype * ddts = malloc(DDT_COUNT * sizeof(*ddts));
  for (int i = 0; i < DDT_COUNT; ++i) {
    MPI_Type_contiguous(i + 1, MPI_DOUBLE, ddts + i);
    MPI_Type_commit(ddts + i);
  }

  {
    MPI_Datatype ddt;
    MPI_Type_contiguous(BASE_COUNT, MPI_DOUBLE, &ddt);
    MPI_Type_commit(&ddt);

    double * send_buffer =
      malloc(EXCH_COUNT * BASE_COUNT * sizeof(*send_buffer));
    double * recv_buffer =
      malloc(RANK_COUNT * EXCH_COUNT * BASE_COUNT * sizeof(*recv_buffer));

    for (int i = 0; i < EXCH_COUNT; ++i)
      for (int j = 0; j < BASE_COUNT; ++j)
        send_buffer[i * BASE_COUNT + j] = (double)rank;
    for (int i = 0; i < RANK_COUNT; ++i)
      for (int j = 0; j < EXCH_COUNT; ++j)
        for (int k = 0; k < BASE_COUNT; ++k)
          recv_buffer[i * EXCH_COUNT * BASE_COUNT +
                      j * BASE_COUNT + k] = -1.0;

    MPI_Request requests[2][RANK_COUNT];

    for (int i = 0; i < RANK_COUNT; ++i)
      MPI_Irecv(
        recv_buffer + i * EXCH_COUNT * BASE_COUNT, EXCH_COUNT,
        ddt, i, 0, MPI_COMM_WORLD, requests[0] + i);
    for (int i = 0; i < RANK_COUNT; ++i)
      MPI_Isend(
        send_buffer, EXCH_COUNT, ddt, i, 0, MPI_COMM_WORLD,
        requests[1] + i);
    MPI_Waitall(2 * RANK_COUNT, requests[0], MPI_STATUSES_IGNORE);

    free(recv_buffer);
    free(send_buffer);
    MPI_Type_free(&ddt);
  }

  for (int i = 0; i < DDT_COUNT; ++i)
    MPI_Type_free(ddts + i);
  free(ddts);

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
