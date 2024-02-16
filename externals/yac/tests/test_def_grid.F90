!!
!! @file test_def_grid.F90
!!
!! @copyright Copyright  (C)  2013 DKRZ, MPI-M
!!
!! @author Moritz Hanke <hanke@dkrz.de>
!!         Rene Redler  <rene.redler@mpimet.mpg.de>
!!
!
! Keywords:
! Maintainer: Moritz Hanke <hanke@dkrz.de>
!             Rene Redler <rene.redler@mpimet.mpg.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
!
! This file is part of YAC.
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!


#include "test_macros.inc"

program test_def_grid

  use utest
  use mo_yac_finterface

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

#ifdef VERBOSE
  print *, ' def_comp returned comp_id ', comp_id
#endif
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

#ifdef VERBOSE
  print *, ' def_grid returned grid_id ', grid_id
#endif
  call test ( grid_id /= -99 )

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

#ifdef VERBOSE
  print *, ' def_grid returned grid_id ', grid_id
#endif
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

#ifdef VERBOSE
  print *, ' def_grid returned grid_id ', grid_id
#endif
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

#ifdef VERBOSE
  print *, ' def_grid returned grid_id ', grid_id
#endif
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

#ifdef VERBOSE
  print *, ' def_grid returned grid_id ', grid_id
#endif
  call test ( grid_id /= -99 )

  call yac_ffinalize

  call stop_test

  call exit_tests

end program test_def_grid
