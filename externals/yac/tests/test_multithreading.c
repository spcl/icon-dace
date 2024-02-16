/**
 * @file test_multithreading.c
 *
 * @copyright Copyright  (C)  2023 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *         Nils-Arne Dreier <dreier@dkrz.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 *             Nils-Arne Dreier <dreier@dkrz.de>
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

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include <mpi.h>

#include "yac_interface.h"

int comp_id, grid_id, point_id, field_id1, field_id2, interp_stack_id;

void defs(const char* comp_name){

  yac_cdef_comp(comp_name, &comp_id);

  int nbr_vertices[2] = {3, 3};
  int cyclic[2] = {0,0};
  double x_vertices[3] = {-1, 0, 1};
  double y_vertices[3] = {-1, 0, 1};
  yac_cdef_grid_reg2d(comp_name, nbr_vertices, cyclic, x_vertices, y_vertices, &grid_id);

  yac_cdef_points_reg2d(grid_id, nbr_vertices, YAC_LOCATION_CORNER,
    x_vertices, y_vertices, &point_id );

  yac_cdef_field ( "field1",
    comp_id,
    &point_id,
    1,
    1,
    "1",
    YAC_TIME_UNIT_SECOND,
    &field_id1);
  yac_cdef_field ( "field2",
    comp_id,
    &point_id,
    1,
    1,
    "1",
    YAC_TIME_UNIT_SECOND,
    &field_id2);

  yac_csync_def();

  yac_cget_interp_stack_config(&interp_stack_id);
  yac_cadd_interp_stack_config_nnn(
    interp_stack_id, YAC_NNN_AVG, 1, 1.0);

  yac_cdef_couple("compA", "compA", "field1",
    "compB", "compB", "field1",
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE,
    interp_stack_id, 0, 0);
  yac_cdef_couple("compA", "compA", "field2",
    "compB", "compB", "field2",
    "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE,
    interp_stack_id, 0, 0);
  yac_cenddef();
}

void put_field_loop(int* ptr_field_id){
  int info = 0, ierror;
  double buffer[9];
  for(int i = 0; i<60; ++i){
    printf("put: field_id: %d, i: %d\n", *ptr_field_id, i);
    buffer[5] = *ptr_field_id + i;
    yac_cput_(*ptr_field_id, 1, buffer, &info, &ierror);
  }
}

void compA(){
  defs("compA");
  pthread_t thread1, thread2;
  pthread_create( &thread1, NULL, (void*(*)(void*)) put_field_loop, (void*) &field_id1);
  pthread_create( &thread2, NULL, (void*(*)(void*)) put_field_loop, (void*) &field_id2);

  pthread_join( thread1, NULL);
  pthread_join( thread2, NULL);
}

void get_field_loop(int* ptr_field_id){
  int info = 0, ierror;
  double buffer[9];
  for(int i = 0; i<60; ++i){
    printf("get: field_id: %d, i: %d\n", *ptr_field_id, i);
    yac_cget_(*ptr_field_id, 1, buffer, &info, &ierror);
    if(buffer[5] != *ptr_field_id + i)
      exit(1);
  }
}

void compB(){
  defs("compB");
  pthread_t thread1, thread2;
  pthread_create( &thread1, NULL, (void*(*)(void*)) get_field_loop, (void*) &field_id1);
  pthread_create( &thread2, NULL, (void*(*)(void*)) get_field_loop, (void*) &field_id2);

  pthread_join( thread1, NULL);
  pthread_join( thread2, NULL);
}

int main(void) {
  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);

  if(provided < MPI_THREAD_MULTIPLE){
    printf("Skip test due to not compatible MPI (MPI_Query_thread: %d)\n", provided);
    return 77;
  }

  yac_cinit();
  yac_cdef_datetime (
    "2023-01-09T14:20:00",
    "2023-01-09T14:21:00" );
  yac_cdef_calendar (YAC_PROLEPTIC_GREGORIAN);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank%2 == 0)
    compA();
  else
    compB();

  yac_cfinalize();
  MPI_Finalize();
  return 0;
}
