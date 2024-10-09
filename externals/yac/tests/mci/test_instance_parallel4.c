// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <mpi.h>

#include "tests.h"
#include "test_common.h"
#include "instance.h"
#include "yac.h"
#include "yac_mpi.h"
#include "yaxt.h"
#include "yac_mpi.h"

static char const * component_names[3] = {"comp_1", "comp_2", "comp_3"};

void check_comm(
  struct yac_instance * instance,
  size_t * comp_idxs, size_t num_comps, int ref_size) ;

int main (void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

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

      MPI_Comm_free(&yac_comm);

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
  MPI_Comm_size(comps_comm, &comps_comm_size);

  if (comps_comm_size != ref_size) PUT_ERR("invalid size of comps_comm");

  MPI_Comm_free(&comps_comm);
}
