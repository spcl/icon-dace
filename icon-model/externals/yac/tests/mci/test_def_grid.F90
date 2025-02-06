! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

program test_def_grid

  use utest
  use yac

  implicit none

  integer :: comp_id
  integer :: grid_id

  integer, parameter :: nbr_vertices          = 6
  integer, parameter :: nbr_vertices_2d(2)    = (/3,2/)
  integer, parameter :: cyclic(2)             = (/0,0/)
  integer, parameter :: nbr_cells             = 2
  integer, parameter :: nbr_vertices_per_cell = 4
  integer, parameter :: nbr_connections       = 7
  integer, parameter :: nbr_connections_ll    = 8

  integer            :: nbr_vertices_per_cell_nu(nbr_cells)
  integer            :: cell_to_vertex_nu(nbr_connections)
  integer            :: cell_to_vertex_nu_ll(nbr_connections_ll)

#ifdef __use_REAL
  real :: x_vertices(nbr_vertices)
  real :: y_vertices(nbr_vertices)
#else
  double precision :: x_vertices(nbr_vertices)
  double precision :: y_vertices(nbr_vertices)
#endif

  integer :: cell_to_vertex(nbr_vertices_per_cell,nbr_cells)

  call start_test("def_grid")

  call yac_finit ( )

  comp_id = -99
  call yac_fdef_comp ( 'ICON-ocean', comp_id )

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

  call test ( 2 == yac_fget_grid_size(yac_location_cell, grid_id) )
  call test ( 6 == yac_fget_grid_size(yac_location_corner, grid_id) )
  call test ( 7 == yac_fget_grid_size(yac_location_edge, grid_id) )

  ! uniform unstructured lonlat grid

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

  call yac_fdef_grid ( 'grid1_ll',            &
                       nbr_vertices,          &
                       nbr_cells,             &
                       nbr_vertices_per_cell, &
                       x_vertices,            &
                       y_vertices,            &
                       cell_to_vertex,        &
                       grid_id,               &
                       .true. )

  print *, ' def_grid returned grid_id ', grid_id
  call test ( grid_id /= -99 )

  ! non-uniform unstructured grid

  !      1
  !     / \
  !    / 1 \
  !   /     \
  !  2-------3
  !  |       |
  !  |   2   |
  !  |       |
  !  4-------5

  x_vertices(1) =  0.0; y_vertices(1) =  1.0
  x_vertices(2) = -0.5; y_vertices(2) =  0.0
  x_vertices(3) =  0.5; y_vertices(3) =  0.0
  x_vertices(4) = -0.5; y_vertices(4) = -1.0
  x_vertices(5) =  0.5; y_vertices(5) = -1.0

  cell_to_vertex_nu(1) = 1
  cell_to_vertex_nu(2) = 2
  cell_to_vertex_nu(3) = 3
  cell_to_vertex_nu(4) = 2
  cell_to_vertex_nu(5) = 4
  cell_to_vertex_nu(6) = 5
  cell_to_vertex_nu(7) = 3

  nbr_vertices_per_cell_nu(1) = 3
  nbr_vertices_per_cell_nu(2) = 4

  call yac_fdef_grid ( 'grid2',                  &
                       nbr_vertices,             &
                       nbr_cells,                &
                       nbr_connections,          &
                       nbr_vertices_per_cell_nu, &
                       x_vertices,               &
                       y_vertices,               &
                       cell_to_vertex_nu,        &
                       grid_id )

  print *, ' def_grid returned grid_id ', grid_id
  call test ( grid_id /= -99 )

  ! non-uniform unstructured lonlat grid

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

  cell_to_vertex_nu_ll(1) = 1
  cell_to_vertex_nu_ll(2) = 3
  cell_to_vertex_nu_ll(3) = 4
  cell_to_vertex_nu_ll(4) = 2
  cell_to_vertex_nu_ll(5) = 3
  cell_to_vertex_nu_ll(6) = 5
  cell_to_vertex_nu_ll(7) = 6
  cell_to_vertex_nu_ll(8) = 4

  nbr_vertices_per_cell_nu(1) = 4
  nbr_vertices_per_cell_nu(2) = 4

  call yac_fdef_grid ( 'grid2_ll',               &
                       nbr_vertices,             &
                       nbr_cells,                &
                       nbr_connections_ll,       &
                       nbr_vertices_per_cell_nu, &
                       x_vertices,               &
                       y_vertices,               &
                       cell_to_vertex_nu_ll,     &
                       grid_id,                  &
                       .true. )

  print *, ' def_grid returned grid_id ', grid_id
  call test ( grid_id /= -99 )

  ! regular 2d grid

  !  3-------4-------5
  !  |       |       |
  !  |   0   |   1   |
  !  |       |       |
  !  0-------1-------2

  x_vertices(1) = 0.0
  x_vertices(2) = 1.0
  x_vertices(3) = 2.0
  y_vertices(1) = 0.0
  y_vertices(2) = 1.0

  call yac_fdef_grid ( 'grid3',         &
                       nbr_vertices_2d, &
                       cyclic,          &
                       x_vertices,      &
                       y_vertices,      &
                       grid_id )

  print *, ' def_grid returned grid_id ', grid_id
  call test ( grid_id /= -99 )

  call yac_ffinalize

  call stop_test

  call exit_tests

end program test_def_grid
