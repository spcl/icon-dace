/**
 * @file test_component_config.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "component.h"
#include "couple_config.h"

static struct yac_couple_config * generate_couple_config(
  char ** comp_names, size_t count);

int main (void) {

  // initialise MPI
  yac_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
  int rank, size;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  if (size != 8) {

    PUT_ERR("ERROR: wrong number of processes\n");
    return TEST_EXIT_CODE;
  }

  //----------------------------------------------------------------------------
  // setup
  //----------------------------------------------------------------------------

  // names of all components
  char const * all_comp_names[] = {"comp_a", "comp_b", "comp_c"};
  // ranks for each component
  int comp_ranks[3][4] = {{0, 1, 2, 4}, {0, 1, 3, 5}, {0, 2, 3, 6}};
  // number of ranks per component
  int num_comp_ranks[3] = {4, 4, 4};

  // get the names of all components that are to be defined on the local process
  char const * local_comp_names[3];
  size_t num_local_comps = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < num_comp_ranks[i]; ++j)
      if (rank == comp_ranks[i][j])
        local_comp_names[num_local_comps++] = all_comp_names[i];

  // generate dummy couple_config
  struct yac_couple_config * couple_config =
    generate_couple_config((char**)all_comp_names, 3);

  // generate component configuration
  struct yac_component_config * comp_config =
    yac_component_config_new(
      couple_config, (char const **)local_comp_names, num_local_comps,
      MPI_COMM_WORLD);

  // free memory in dummy couple_config
  yac_couple_config_delete(couple_config);

  //----------------------------------------------------------------------------
  // generate reference data
  //----------------------------------------------------------------------------

  // generate reference communicators for each component and for each
  // component pair
  MPI_Comm ref_comp_comm[3], ref_comp_pair_comm[3], ref_all_comps_comm;
  {
    MPI_Group world_group, comp_group[3];
    yac_mpi_call(
      MPI_Comm_group(MPI_COMM_WORLD, &world_group), MPI_COMM_WORLD);
    for (int i = 0; i < 3; ++i) {
      yac_mpi_call(
        MPI_Group_incl(
          world_group, num_comp_ranks[i], comp_ranks[i], &comp_group[i]),
        MPI_COMM_WORLD);
      yac_mpi_call(
        MPI_Comm_create(MPI_COMM_WORLD, comp_group[i], &ref_comp_comm[i]),
        MPI_COMM_WORLD);
    }
    for (int i = 0, k = 0; i < 3; ++i) {
      for (int j = i + 1; j < 3; ++j, ++k) {
        MPI_Group comp_pair_group;
        yac_mpi_call(
          MPI_Group_union(comp_group[i], comp_group[j], &comp_pair_group),
          MPI_COMM_WORLD);
        yac_mpi_call(
          MPI_Comm_create(
            MPI_COMM_WORLD, comp_pair_group, &ref_comp_pair_comm[k]),
          MPI_COMM_WORLD);
        yac_mpi_call(MPI_Group_free(&comp_pair_group), MPI_COMM_WORLD);
      }
    }
    {
      MPI_Group all_comps_group = MPI_GROUP_EMPTY;
      for (int i = 0; i < 3; ++i) {
        MPI_Group merge_group;
        yac_mpi_call(
          MPI_Group_union(all_comps_group, comp_group[i], &merge_group),
          MPI_COMM_WORLD);
        if (all_comps_group != MPI_GROUP_EMPTY)
          yac_mpi_call(MPI_Group_free(&all_comps_group), MPI_COMM_WORLD);
        all_comps_group = merge_group;
      }
      yac_mpi_call(
        MPI_Comm_create(
          MPI_COMM_WORLD, all_comps_group, &ref_all_comps_comm),
        MPI_COMM_WORLD);
      yac_mpi_call(MPI_Group_free(&all_comps_group), MPI_COMM_WORLD);
    }
    for (int i = 0; i < 3; ++i)
      yac_mpi_call(MPI_Group_free(&comp_group[i]), MPI_COMM_WORLD);
    yac_mpi_call(MPI_Group_free(&world_group), MPI_COMM_WORLD);
  }

  //----------------------------------------------------------------------------
  // testing
  //----------------------------------------------------------------------------

  if (yac_component_config_get_num_components(comp_config) != num_local_comps)
    PUT_ERR("error in yac_component_config_get_num_components");

  struct component ** components =
    yac_component_config_get_components(comp_config);

  for (size_t i = 0; i < num_local_comps; ++i) {

    if (strcmp(local_comp_names[i],  yac_component_get_name(components[i])))
      PUT_ERR("error in yac_component_get_name");

    int compare_result;
    yac_mpi_call(
      MPI_Comm_compare(
        ref_comp_comm[local_comp_names[i][5]-'a'],
        yac_component_get_communicator(components[i]), &compare_result),
      MPI_COMM_WORLD);
    if (compare_result != MPI_CONGRUENT)
      PUT_ERR("error in yac_component_get_communicator");
  }

  // encode into a single int which components are defined locally
  int comp_flags = 0;
  for (size_t i = 0; i < num_local_comps; ++i)
    comp_flags |= (1 << (local_comp_names[i][5] - 'a'));

  // check component pair communicators
  for (int i = 0, k = 0; i < 3; ++i) { // for all component
    for (int j = i + 1; j < 3; ++j, ++k) { // for all remaining components
      if ((comp_flags & (1 << i)) || (comp_flags & (1 << j))) {
        MPI_Comm comp_pair_comm =
          yac_component_config_get_comps_comm(
            comp_config,
            (char const *[]){all_comp_names[i], all_comp_names[j]}, 2);
        int compare_result;
        yac_mpi_call(
          MPI_Comm_compare(
            ref_comp_pair_comm[k], comp_pair_comm, &compare_result),
          MPI_COMM_WORLD);
        if (compare_result != MPI_CONGRUENT)
          PUT_ERR("error in yac_component_config_get_comps_comm");
        yac_mpi_call(MPI_Comm_free(&comp_pair_comm), MPI_COMM_WORLD);
      }
    }
  }

  // check communicators containing all processes
  if (ref_all_comps_comm != MPI_COMM_NULL) {
    MPI_Comm all_comps_comm =
      yac_component_config_get_comps_comm(comp_config, all_comp_names, 3);
    int compare_result;
    yac_mpi_call(
      MPI_Comm_compare(ref_all_comps_comm, all_comps_comm, &compare_result),
      MPI_COMM_WORLD);
    if (compare_result != MPI_CONGRUENT)
      PUT_ERR("error in yac_component_config_get_comps_comm");
    yac_mpi_call(MPI_Comm_free(&all_comps_comm), MPI_COMM_WORLD);
  }

  //----------------------------------------------------------------------------
  // clean-up
  //----------------------------------------------------------------------------

  yac_component_config_delete(comp_config);

  if (ref_all_comps_comm != MPI_COMM_NULL)
    yac_mpi_call(MPI_Comm_free(&ref_all_comps_comm), MPI_COMM_WORLD);
  for (int i = 0; i < 3; ++i) {
    if (ref_comp_comm[i] != MPI_COMM_NULL)
      yac_mpi_call(MPI_Comm_free(&ref_comp_comm[i]), MPI_COMM_WORLD);
    if (ref_comp_pair_comm[i] != MPI_COMM_NULL)
      yac_mpi_call(MPI_Comm_free(&ref_comp_pair_comm[i]), MPI_COMM_WORLD);
  }

  // finalize MPI
  yac_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);

  return TEST_EXIT_CODE;
}

static struct yac_couple_config * generate_couple_config(
  char ** comp_names, size_t count) {

  struct yac_couple_config * couple_config = yac_couple_config_new();
  for (size_t i = 0; i < count; ++i)
    yac_couple_config_add_component(couple_config, comp_names[i]);
  return couple_config;
}
