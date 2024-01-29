/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#include <stdio.h>

#include <mpi.h>

#include "../include/mpi_handshake.h"

int main(int argc, char** argv){
  MPI_Init(&argc, &argv);

  // first handshake
  const char* group_name = "comin";
  MPI_Comm comin_comm;
  mpi_handshake(&group_name, &comin_comm, 1, MPI_COMM_WORLD);

  // second handshake the name for the comin communicator is read from the command line
  char external_name[255];
  sprintf(external_name, "external_%s", argv[1]);
  const char* group_names[3] = {argv[1], external_name, "all_external"};
  MPI_Comm comms[3];
  mpi_handshake(group_names, comms, 3, comin_comm);

  int plugin_size, external_size, all_external_size;
  MPI_Comm_size(comms[0], &plugin_size);
  MPI_Comm_size(comms[1], &external_size);
  MPI_Comm_size(comms[2], &all_external_size);
  int rank;
  MPI_Comm_rank(comms[1], &rank);

  if(rank == 0){
    printf("plugin(%s) size: %d\n", argv[1], plugin_size);
    printf("this external size: %d\n", external_size);
    printf("all external size: %d\n", all_external_size);
  }

  MPI_Comm_free(&comms[0]);
  MPI_Comm_free(&comms[1]);
  MPI_Comm_free(&comin_comm);
  MPI_Finalize();
}
