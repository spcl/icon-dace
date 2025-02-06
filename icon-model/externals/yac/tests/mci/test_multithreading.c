// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include <mpi.h>

#include "yac.h"

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
    interp_stack_id, YAC_NNN_AVG, 1, 0.0, 1.0);

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

void * put_field_loop(void* ptr_field_id_){
  int info = 0, ierror;
  double buffer[9];
  int* ptr_field_id = ptr_field_id_;
  for(int i = 0; i<60; ++i){
    printf("put: field_id: %d, i: %d\n", *ptr_field_id, i);
    buffer[5] = *ptr_field_id + i;
    yac_cput_(*ptr_field_id, 1, buffer, &info, &ierror);
  }
  return NULL;
}

void compA(){
  defs("compA");
  pthread_t thread1, thread2;
  pthread_create( &thread1, NULL, &put_field_loop, (void*) &field_id1);
  pthread_create( &thread2, NULL, &put_field_loop, (void*) &field_id2);

  pthread_join( thread1, NULL);
  pthread_join( thread2, NULL);
}

void * get_field_loop(void* ptr_field_id_){
  int info = 0, ierror;
  double buffer[9];
  int* ptr_field_id = ptr_field_id_;
  for(int i = 0; i<60; ++i){
    printf("get: field_id: %d, i: %d\n", *ptr_field_id, i);
    yac_cget_(*ptr_field_id, 1, buffer, &info, &ierror);
    if(buffer[5] != *ptr_field_id + i)
      exit(1);
  }
  return NULL;
}

void compB(){
  defs("compB");
  pthread_t thread1, thread2;
  pthread_create( &thread1, NULL, &get_field_loop, (void*) &field_id1);
  pthread_create( &thread2, NULL, &get_field_loop, (void*) &field_id2);

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
  yac_cdef_calendar (YAC_PROLEPTIC_GREGORIAN);
  yac_cdef_datetime (
    "2023-01-09T14:20:00",
    "2023-01-09T14:21:00" );
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
