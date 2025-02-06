! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

program test_def_points

  use utest
  use yac

  implicit none

  integer :: comp_id
  integer :: grid_id
  integer :: point_id
  integer :: location

  integer, parameter :: nbr_points = 2

  integer, parameter :: nbr_vertices = 6
  integer, parameter :: nbr_cells    = 2
  integer, parameter :: nbr_vertices_per_cell = 4
  integer            :: cell_to_vertex(nbr_vertices_per_cell,nbr_cells)

#ifdef __use_REAL
  real :: x_vertices(nbr_vertices)
  real :: y_vertices(nbr_vertices)
#else
  double precision :: x_vertices(nbr_vertices)
  double precision :: y_vertices(nbr_vertices)
#endif

#ifdef __use_REAL
  real :: x_points(nbr_points)
  real :: y_points(nbr_points)
#else
  double precision :: x_points(nbr_points)
  double precision :: y_points(nbr_points)
#endif

  call yac_finit ( )

  ! set up dummy component
  call yac_fdef_comp ( "ICON-ocean", comp_id )

  print *, ' def_comp returned comp_id ', comp_id
  call test ( comp_id /= -99 )

  ! uniform unstructured grid

  !  1-------2
  !  |       |
  !  |   1   |
  !  |       |
  !  3-------4
  !  |       |
  !  |   2   |
  !  |       |
  !  5-------6

  x_vertices(1) = -0.5; y_vertices(1) =  1.0
  x_vertices(2) =  0.5; y_vertices(2) =  1.0
  x_vertices(3) = -0.5; y_vertices(3) =  0.0
  x_vertices(4) =  0.5; y_vertices(4) =  0.0
  x_vertices(5) = -0.5; y_vertices(5) = -1.0
  x_vertices(6) =  0.5; y_vertices(6) = -1.0

  cell_to_vertex(1,1) = 1
  cell_to_vertex(2,1) = 3
  cell_to_vertex(3,1) = 4
  cell_to_vertex(4,1) = 2
  cell_to_vertex(1,2) = 3
  cell_to_vertex(2,2) = 5
  cell_to_vertex(3,2) = 6
  cell_to_vertex(4,2) = 4

  grid_id = -99

  call yac_fdef_grid ( 'grid1',               &
                       nbr_vertices,          &
                       nbr_cells,             &
                       nbr_vertices_per_cell, &
                       x_vertices,            &
                       y_vertices,            &
                       cell_to_vertex,        &
                       grid_id )

  print *, ' def_grid returned grid_id ', grid_id
  call test ( grid_id /= -99 )

  call start_test("def_points")

  location = YAC_LOCATION_CELL ! one point per cell

  x_points(1) =   0.0
  y_points(1) =   0.5
  x_points(2) =   0.0
  y_points(2) =  -0.5

  point_id = -99

  call yac_fdef_points ( grid_id, &
                         nbr_points,   &
                         location,     &
                         x_points,     &
                         y_points,     &
                         point_id )

  print *, ' def_points returned point_id ', point_id
  call test ( point_id /= -99 )

  call test ( 2 == yac_fget_points_size( point_id ))

  call yac_ffinalize ( )

  call stop_test

  call exit_tests

end program test_def_points
