#define VERBOSE

/**
 * @file dynamic_config.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Nils-Arne Dreier <dreier@dkrz.de>
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
#include "yac_interface.h"

int main () {
  // Initialisation of MPI

  MPI_Init (0, NULL);

  /* The initialisation phase */
#ifdef VERBOSE
  printf (". main: calling yac_cinit\n");
#endif

  // start with an empty config
  yac_cinit ();
  yac_cdef_datetime("2022-11-11T13:18:00", "2022-11-12T13:18:00");
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);

  int global_rank, global_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&global_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&global_size);

  /* Dynamic component definition. */

#ifdef VERBOSE
  printf (". main: calling yac_cdef_comp\n");
#endif

  int comp_id;
  char * comp_name;
  if(global_rank == 0)
    comp_name = "ICON-atmosphere";
  else if(global_rank == 1)
    comp_name = "ICON-ocean";
  else if(global_rank == 2)
    comp_name = "VIZ";
  else
    comp_name = "OUTPUT";

  yac_cdef_comp ( comp_name, &comp_id );
#ifdef VERBOSE
  printf ( ". main: defined %s with local comp ID %i \n", comp_name, comp_id );
#endif

  MPI_Comm localcomm;
  yac_cget_comp_comm(comp_id, &localcomm);
  int comp_size;
  MPI_Comm_size(localcomm, &comp_size);
#ifdef VERBOSE
  printf ( ". main: component %s has %i processes\n", comp_name, comp_size );
#endif

  MPI_Comm_free(&localcomm);

  // register grids (for simplicity one grid per comp)
  int num_vert[2] = {3,3};
  int cyclic[2] = {0,0};
  double vertices[3] = {0, 1., 2.};
  int grid_id;
  yac_cdef_grid_reg2d(comp_name, num_vert, cyclic,
    vertices, vertices, &grid_id);

  int point_id;
  yac_cdef_points_reg2d(grid_id,
    num_vert,YAC_LOCATION_CORNER, vertices, vertices, &point_id);

  int t_id, t2_id, s_id;
  yac_cdef_field(
    "t", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &t_id);
  yac_cdef_field(
    "t2", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &t2_id);
  if(global_rank < 2)
    yac_cdef_field(
      "s", comp_id, &point_id, 1, 1, "1", YAC_TIME_UNIT_SECOND, &s_id);

  if(global_rank != 1)
    yac_csync_def();
  printf("synced!\n");

  if(global_rank == 0){
    int interp_stack_config_id;
    yac_cget_interp_stack_config(&interp_stack_config_id);
    yac_cadd_interp_stack_config_nnn(interp_stack_config_id,
		YAC_NNN_AVG, 1,1.);
    if(global_size > 1){
      yac_cdef_couple("ICON-atmosphere", "ICON-atmosphere", "t", // source (comp, grid, field)
        "ICON-ocean", "ICON-ocean", "t", // target (comp, grid, field)
        "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_AVERAGE,
        interp_stack_config_id, 0, 0);
      yac_cdef_couple("ICON-atmosphere", "ICON-atmosphere", "t2", // source (comp, grid, field)
        "ICON-ocean", "ICON-ocean", "t2", // target (comp, grid, field)
        "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_AVERAGE,
        interp_stack_config_id,  0, 0);
      yac_cdef_couple("ICON-ocean", "ICON-ocean", "s", // source (comp, grid, field)
        "ICON-atmosphere", "ICON-atmosphere", "s", // target (comp, grid, field)
        "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_AVERAGE,
        interp_stack_config_id,  0, 0);
      yac_cfree_interp_stack_config(interp_stack_config_id);
    }
  }else if (global_rank == 2){
    int interp_stack_config_id;
    yac_cget_interp_stack_config(&interp_stack_config_id);
    yac_cadd_interp_stack_config_nnn(interp_stack_config_id,
		YAC_NNN_AVG, 1,1.);
    yac_cdef_couple("ICON-atmosphere", "ICON-atmosphere", "t", // source (comp, grid, field)
      "VIZ", "VIZ", "t", // target (comp, grid, field)
      "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_AVERAGE,
      interp_stack_config_id,  0, 0);
    if(global_size > 3) // add couple for two other components
      yac_cdef_couple("ICON-atmosphere", "ICON-atmosphere", "t", // source (comp, grid, field)
        "OUTPUT", "OUTPUT", "t", // target (comp, grid, field)
        "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_AVERAGE,
        interp_stack_config_id,  0, 0);
    yac_cfree_interp_stack_config(interp_stack_config_id);
  }

  int ierr;
  yac_cenddef( );

  double * data = malloc(sizeof(double) * num_vert[0]*num_vert[1]);
  int info;
  if(global_rank == 0){
    for(size_t i = 0; i < num_vert[0]*num_vert[1]; ++i)
      data[i] = 42.+i;
    double ** data_ptr = &data;
    yac_cput(t_id, 1, &data_ptr, &info, &ierr);
    for(size_t i = 0; i < num_vert[0]*num_vert[1]; ++i)
      data[i] = -42.+i;
    yac_cput(t2_id, 1, &data_ptr, &info, &ierr);
  }else{
    yac_cget(t2_id, 1, &data, &info, &ierr);
    printf("info: %d\nerr: %d\n%f\n", info, ierr, data[2]);
    yac_cget(t_id, 1, &data, &info, &ierr);
    printf("info: %d\nerr: %d\n%f\n", info, ierr, data[2]);
  }

  yac_cfinalize();

  MPI_Finalize();

  return EXIT_SUCCESS;
}
