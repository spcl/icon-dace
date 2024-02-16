/**
 * @file test_instance_parallel4.c
 *
 * @copyright Copyright  (C)  2020 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <mpi.h>

#include "tests.h"
#include "test_common.h"
#include "instance.h"
#include "utils.h"
#include "yac_interface.h"
#include "yac_mpi.h"
#include "yaxt.h"
#include "yac_mpi.h"

static char const * component_names[3] = {"comp_1", "comp_2", "comp_3"};

void check_comm(
  struct yac_instance * instance,
  size_t * comp_idxs, size_t num_comps, int ref_size) ;

int main (void) {

  yac_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  int rank, size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  if (size != 9) {

    PUT_ERR("ERROR: wrong number of processes\n");
    return TEST_EXIT_CODE;
  }

  { // three component, defined on disjoint sets of processes
    struct yac_instance * instance =
      yac_instance_new(MPI_COMM_WORLD);

    yac_instance_def_components(
      instance, &(component_names[rank / 3]), 1);

    // comm between comp_1 and comp_2
    if (rank < 6) check_comm(instance, (size_t[]){0, 1}, 2, 6);

    // comm between comp_1 and comp_3
    if ((rank < 3) || (rank >= 6))
      check_comm(instance, (size_t[]){0, 2}, 2, 6);

    // comm between comp_2 and comp_3
    if (rank >= 3) check_comm(instance, (size_t[]){1, 2}, 2, 6);

    // comm between comp_1, comp_2, and comp_3
    check_comm(instance, (size_t[]){0, 1, 2}, 3, 9);

    yac_instance_delete(instance);
  }

  { // two components sharing all processes
    struct yac_instance * instance =
      yac_instance_new(MPI_COMM_WORLD);

    char const * comp_names[2] = {component_names[0], component_names[1]};
    yac_instance_def_components(instance, comp_names, 2);

    check_comm(instance, (size_t[]){0, 1}, 2, 9);

    yac_instance_delete(instance);
  }

  { // three components one is defined on all processes, the other two are
    // a subset
    struct yac_instance * instance =
      yac_instance_new(MPI_COMM_WORLD);

    int num_comps;
    char const * comp_names[2];

    comp_names[0] = component_names[0];
    if (rank < 3) {
      comp_names[1] = component_names[1];
      num_comps = 2;
    } else if (rank < 6) {
      num_comps = 1;
    } else {
      comp_names[1] = component_names[2];
      num_comps = 2;
    }

    yac_instance_def_components(instance, comp_names, num_comps);

    // comm between comp_1 and comp_2
    check_comm(instance, (size_t[]){0, 1}, 2, 9);

    // comm between comp_1 and comp_3
    check_comm(instance, (size_t[]){0, 2}, 2, 9);

    // comm between comp_2 and comp_3
    if ((rank < 3) || (rank >= 6))
      check_comm(instance, (size_t[]){1, 2}, 2, 6);

    // comm between comp_1, comp_2, and comp_3
    check_comm(instance, (size_t[]){0, 1, 2}, 3, 9);

    yac_instance_delete(instance);
  }

  { // three components, one having its one processes, the other two sharing
    // some
    struct yac_instance * instance =
      yac_instance_new(MPI_COMM_WORLD);

    int num_comps = 0;
    char const * comp_names[2];

    if (rank < 3)      comp_names[num_comps++] = component_names[0];
    else if (rank < 8) comp_names[num_comps++] = component_names[1];
    if (rank >= 6) comp_names[num_comps++] = component_names[2];

    yac_instance_def_components(instance, comp_names, num_comps);

    // comm between comp_1 and comp_2
    if (rank < 8) check_comm(instance, (size_t[]){0, 1}, 2, 8);

    // comm between comp_1 and comp_3
    if ((rank < 3) || (rank >= 6))
      check_comm(instance, (size_t[]){0, 2}, 2, 6);

    // comm between comp_2 and comp_3
    if (rank >= 3) check_comm(instance, (size_t[]){1, 2}, 2, 6);

    // comm between comp_1, comp_2, and comp_3
    check_comm(instance, (size_t[]){0, 1, 2}, 3, 9);

    yac_instance_delete(instance);
  }

  { // two components and a some processes that do not define any YAC instance

    if ((rank % 3) == 2) {

      yac_cmpi_handshake(MPI_COMM_WORLD, 0, NULL, NULL);

    } else {

      MPI_Comm yac_comm;
      char const * yac_groupname = "yac";

      yac_cmpi_handshake(MPI_COMM_WORLD, 1, &yac_groupname, &yac_comm);

      struct yac_instance * instance =
        yac_instance_new(yac_comm);

      yac_mpi_call(MPI_Comm_free(&yac_comm), MPI_COMM_WORLD);

      char const * comp_names[] = {component_names[(rank & 1)]};

      yac_instance_def_components(instance, comp_names, 1);

      // comm between comp_1 and comp_2
      check_comm(instance, (size_t[]){0, 1}, 2, 6);

      yac_instance_delete(instance);
    }
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

void check_comm(
  struct yac_instance * instance,
  size_t * comp_idxs, size_t num_comps, int ref_size) {

  char const ** comp_names = xmalloc(num_comps * sizeof(*comp_names));
  for (size_t i = 0; i < num_comps; ++i)
    comp_names[i] = component_names[comp_idxs[i]];

  MPI_Comm comps_comm =
    yac_instance_get_comps_comm(instance, comp_names, num_comps);

  free(comp_names);

  int comps_comm_size;
  yac_mpi_call(MPI_Comm_size(comps_comm, &comps_comm_size), MPI_COMM_WORLD);

  if (comps_comm_size != ref_size) PUT_ERR("invalid size of comps_comm");

  yac_mpi_call(MPI_Comm_free(&comps_comm), MPI_COMM_WORLD);
}
