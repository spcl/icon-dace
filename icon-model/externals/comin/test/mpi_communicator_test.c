/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#include <stdbool.h>
#include <stdio.h>

#include <comin.h>
#include <mpi.h>

MPI_Comm host_comm, plugin_comm;
int rank;

void mpi_communicator_test_setup(){
  host_comm = MPI_Comm_f2c(comin_parallel_get_host_mpi_comm());
  plugin_comm = MPI_Comm_f2c(comin_parallel_get_plugin_mpi_comm());

  int host_comm_size, plugin_comm_size;
  MPI_Comm_size(host_comm, &host_comm_size);
  MPI_Comm_size(plugin_comm, &plugin_comm_size);
  MPI_Comm_rank(host_comm, &rank);

  if (rank == 0){
    fprintf(stderr, "host_comm has %d processes\n", host_comm_size);
    fprintf(stderr, "plugin_comm has %d processes\n", plugin_comm_size);
  }
}
