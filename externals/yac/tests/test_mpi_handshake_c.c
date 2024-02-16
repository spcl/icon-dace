/**
 * @file test_mpi_handshake_c.c
 *
 * @copyright Copyright  (C)  2022 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
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

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <yaxt.h>

#include "yac_interface.h"
#include "tests.h"

int main (void) {

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);

  int global_rank, global_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &global_size);

  if (global_size != 6) {
    fprintf(stderr, "Wrong number of processes (should be 6)\n");
    exit(EXIT_FAILURE);
  }

  {
    int color = MPI_UNDEFINED;
    MPI_Comm group_comm = MPI_COMM_NULL, main_comm = MPI_COMM_NULL;
    char const * groupMain = "main";
    char const * groups[2];
    MPI_Comm comms_out[2];
    switch (global_rank) {
      default:
      case (0):
      case (1):
        // first two ranks are part of group A
        groups[0] = groupMain;
        groups[1] = "group a";
        yac_cmpi_handshake(MPI_COMM_WORLD, 2, groups, comms_out);
        main_comm = comms_out[0];
        group_comm = comms_out[1];
        color = 0;
        break;
      case (2):
      case (3):
        // next two ranks are part of group B
        groups[0] = "group b";
        groups[1] = groupMain;
        yac_cmpi_handshake(MPI_COMM_WORLD, 2, groups, comms_out);
        group_comm = comms_out[0];
        main_comm = comms_out[1];
        color = 1;
        break;
      case (4):
        // rank four does only provide the main group
        yac_cmpi_handshake(MPI_COMM_WORLD, 1, &groupMain, &main_comm);
        break;
      case (5):
        // rank five does a dummy initialisation
        yac_cmpi_handshake(MPI_COMM_WORLD, 0, NULL, NULL);
      break;
    }

    MPI_Comm ref_group_comm, ref_main_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &ref_group_comm);
    MPI_Comm_split(MPI_COMM_WORLD, (global_rank == 5)?MPI_UNDEFINED:0, 0,
                   &ref_main_comm);

    if (ref_group_comm != MPI_COMM_NULL) {

      int result;
      MPI_Comm_compare(ref_group_comm, group_comm, &result);

      if ((result != MPI_CONGRUENT) && (result != MPI_IDENT))
        PUT_ERR("ERROR in yac_cmpi_handshake");

      MPI_Comm_free(&ref_group_comm);
      MPI_Comm_free(&group_comm);

    } else if (group_comm != MPI_COMM_NULL) {
      PUT_ERR("ERROR in yac_cmpi_handshake");
    }

    if (ref_main_comm != MPI_COMM_NULL) {

      int result;
      MPI_Comm_compare(ref_main_comm, main_comm, &result);

      if ((result != MPI_CONGRUENT) && (result != MPI_IDENT))
        PUT_ERR("ERROR in yac_cmpi_handshake");

      MPI_Comm_free(&ref_main_comm);
      MPI_Comm_free(&main_comm);

    } else if (group_comm != MPI_COMM_NULL) {
      PUT_ERR("ERROR in yac_cmpi_handshake");
    }
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}
