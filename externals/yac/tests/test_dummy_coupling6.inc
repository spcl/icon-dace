!
! @file test_dummy_coupling6.f90
!
! @copyright Copyright  (C)  2022 DKRZ, MPI-M
!
! @author Moritz Hanke <hanke@dkrz.de>
!         Rene Redler  <rene.redler@mpimet.mpg.de>
!
!
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

#define DUMMY_VALUE (-1337.0)

#include "test_macros.inc"

#ifndef TEST_PRECISION
#error "TEST_PRECISION is not defined"
#endif

#ifndef YAC_PTR_TYPE
#error "YAC_PTR_TYPE is not defined"
#endif

PROGRAM main

  USE utest
  USE mo_yac_finterface
  USE yaxt
  USE mpi

  IMPLICIT NONE

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)   !< single precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307) !< double precision
  INTEGER, PARAMETER :: wp = TEST_PRECISION              !< selected working precision

  ! M_PI / 180
  REAL(wp), PARAMETER :: YAC_RAD = 0.017453292519943295769_wp

  CHARACTER (LEN=YAC_MAX_CHARLEN), PARAMETER :: comp_A = 'comp_A'
  CHARACTER (LEN=YAC_MAX_CHARLEN), PARAMETER :: comp_B = 'comp_B'

  INTEGER :: rank, ierror

  INTEGER :: group_comms(1)
  CHARACTER(len=YAC_MAX_CHARLEN) :: group_names(1)

  ! ===================================================================

  CALL MPI_Init(ierror)
  CALL xt_initialize(MPI_COMM_WORLD)
  CALL start_test('dummy_coupling6')

  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

  SELECT CASE (rank)
    CASE (0)
      group_names(1) = 'yac'
      CALL yac_fmpi_handshake(MPI_COMM_WORLD, group_names, group_comms)
      CALL yac_finit_comm_dummy(group_comms(1))
      CALL MPI_Comm_free(group_comms(1), ierror)
    CASE (1)
      CALL run_comp_a()
    CASE (2)
      CALL run_comp_b()
    CASE DEFAULT
      CALL test(.FALSE.)
  END SELECT

  CALL xt_finalize()
  CALL MPI_Finalize(ierror)

  CALL stop_test
  CALL exit_tests

CONTAINS

  SUBROUTINE check_pair_comm(pair_comm, proc_name)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: pair_comm

    CHARACTER, INTENT(IN) :: proc_name

    CHARACTER(LEN=1) :: sendbuf
    CHARACTER(LEN=2) :: recvbuf
    INTEGER :: ierror

    sendbuf(1:1) = proc_name

    CALL MPI_Allgather( &
      sendbuf, 1, MPI_CHARACTER, recvbuf, 1, MPI_CHARACTER, pair_comm, ierror)

    CALL test(recvbuf == "AB")

  END SUBROUTINE

  SUBROUTINE run_comp_a()

    IMPLICIT NONE

    INTEGER :: comp_id, grid_id, cell_point_id, vertex_point_id, &
               cell_mask_id, vertex_mask_id, field_id, multi_field_id
    INTEGER :: pair_comm

    INTEGER, PARAMETER :: nbr_vertices = 4
    INTEGER, PARAMETER :: nbr_cells = 2
    INTEGER, DIMENSION(2), PARAMETER :: nbr_vertices_per_cell = (/3,3/)
    REAL(wp), DIMENSION(4), PARAMETER :: &
      x_vertices =                               &
        (/0.0 * YAC_RAD, -1.0 * YAC_RAD, 1.0 * YAC_RAD, 0.0 * YAC_RAD/)
    REAL(wp), DIMENSION(4), PARAMETER :: &
      y_vertices =                               &
        (/1.0 * YAC_RAD, 0.0 * YAC_RAD, 0.0 * YAC_RAD, -1.0 * YAC_RAD/)
    INTEGER, DIMENSION(6), PARAMETER :: cell_to_vertex = (/1,2,3, 2,4,3/)
    REAL(wp), DIMENSION(2), PARAMETER :: &
      x_cells = (/0.0 * YAC_RAD, 0.0 * YAC_RAD/)
    REAL(wp), DIMENSION(2), PARAMETER :: &
      y_cells = (/0.5 * YAC_RAD, -0.5 * YAC_RAD/)
    REAL(wp), TARGET :: cell_data(nbr_cells), vertex_data(nbr_vertices)

    INTEGER, PARAMETER :: nbr_hor_points = 2
    INTEGER, PARAMETER :: collection_size = 1
    INTEGER, PARAMETER :: nbr_pointsets = 1
    INTEGER :: info

    TYPE(YAC_PTR_TYPE) :: send_field_ptr(2, collection_size)

    ! initialise YAC default instance
    CALL yac_finit()
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fread_config_json('coupling_test6.json')

    ! define local component
    CALL yac_fdef_comp(comp_A, comp_id)

    ! get communicator that contains both components
    CALL yac_fget_comps_comm((/comp_A, comp_B/), 2, pair_comm)

    ! check the pair_comm
    CALL check_pair_comm(pair_comm, 'A')
    CALL MPI_Comm_free(pair_comm, ierror)

    !   define local grid
    !
    !      0
    !     / \
    !    / o \
    !   /     \
    !  1-------2   Eq.
    !   \     /
    !    \ o /
    !     \ /
    !      3
    CALL yac_fdef_grid(                                              &
      'grid_A', nbr_vertices, nbr_cells, SUM(nbr_vertices_per_cell), &
      nbr_vertices_per_cell, x_vertices, y_vertices, cell_to_vertex, grid_id)

    ! center points in cells (needed e.g. for nearest neighbour)
    CALL yac_fdef_points( &
      grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, cell_point_id)
    ! vertex points
    CALL yac_fdef_points( &
      grid_id, nbr_vertices, YAC_LOCATION_CORNER, x_vertices, y_vertices, &
      vertex_point_id)

    ! masks
    CALL yac_fdef_mask( &
      grid_id, nbr_cells, YAC_LOCATION_CELL, (/.TRUE.,.TRUE./), cell_mask_id)
    CALL yac_fdef_mask( &
      grid_id, nbr_vertices, YAC_LOCATION_CORNER, (/1,1,1,1/), vertex_mask_id)

    ! define field
    CALL yac_fdef_field( &
      'A_to_B_src', comp_id, (/cell_point_id/), 1, 1, "1", YAC_TIME_UNIT_SECOND, &
      field_id)
    CALL yac_fdef_field_mask( &
      'A_to_B_multi', comp_id, (/cell_point_id, vertex_point_id/), &
      (/cell_mask_id, vertex_mask_id/), 2, collection_size, "1", &
      YAC_TIME_UNIT_SECOND, multi_field_id)

    ! generate coupling
    CALL yac_fenddef()

    ! move data from comp_A to comp_B
    cell_data = (/3.0, 4.0/)
    vertex_data = (/3.0, 4.0, 5.0, 6.0/)
    CALL yac_fput(                                                            &
      field_id, nbr_hor_points, nbr_pointsets, collection_size,               &
      RESHAPE(cell_data, (/nbr_hor_points, nbr_pointsets, collection_size/)), &
      info, ierror)
    send_field_ptr(1,1)%p => cell_data
    send_field_ptr(2,1)%p => vertex_data
    CALL yac_fput( &
      multi_field_id, 2, collection_size, send_field_ptr, info, ierror)

    CALL yac_ffinalize()

  END SUBROUTINE run_comp_a

  SUBROUTINE run_comp_b()

    IMPLICIT NONE

    INTEGER :: default_instance_id, comp_id, grid_id_B, point_id_B, &
               field_id, multi_field_id
    INTEGER :: pair_comm

    INTEGER :: internal_comm
    INTEGER :: internal_instance_id, internal_comp_id, grid_id_C, point_id_C
    INTEGER :: field_id_src, field_id_tgt

    INTEGER, DIMENSION(2), PARAMETER :: nbr_vertices = (/2,3/)
    INTEGER, DIMENSION(2), PARAMETER :: cyclic = (/0,0/)
    INTEGER, DIMENSION(2), PARAMETER :: nbr_cells = (/1,2/)
    REAL(wp), DIMENSION(2), PARAMETER :: &
      x_vertices = (/-0.5 * YAC_RAD, 0.5 * YAC_RAD/)
    REAL(wp), DIMENSION(3), PARAMETER :: &
      y_vertices = (/-1.0 * YAC_RAD, 0.0 * YAC_RAD, 1.0 * YAC_RAD/)
    REAL(wp), DIMENSION(1), PARAMETER :: &
      x_cells = (/0.0 * YAC_RAD/)
    REAL(wp), DIMENSION(2), PARAMETER :: &
      y_cells = (/-0.5 * YAC_RAD, 0.5 * YAC_RAD/)

    REAL(wp), DIMENSION(2,3), PARAMETER :: &
      internal_x_vertices = &
        RESHAPE((/-0.4 * YAC_RAD, 0.4 * YAC_RAD, &
                  -0.4 * YAC_RAD, 0.4 * YAC_RAD, &
                  -0.4 * YAC_RAD, 0.4 * YAC_RAD/), (/2,3/))
    REAL(wp), DIMENSION(2,3), PARAMETER :: &
      internal_y_vertices = &
        RESHAPE((/-0.7 * YAC_RAD, -0.7 * YAC_RAD, &
                   0.0 * YAC_RAD,  0.0 * YAC_RAD, &
                   0.8 * YAC_RAD,  0.8 * YAC_RAD/), (/2,3/))

    INTEGER, PARAMETER :: nbr_hor_points = 2
    INTEGER, PARAMETER :: collection_size = 1
    REAL(wp), DIMENSION(nbr_hor_points,collection_size), TARGET :: recv_field_B
    INTEGER :: info

    REAL(wp), DIMENSION(6,collection_size) :: recv_field_C
    INTEGER :: send_info, recv_info
    TYPE(YAC_PTR_TYPE) :: recv_field_ptr(collection_size)

    ! initialise YAC default instance
    CALL yac_finit(default_instance_id)
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fread_config_json(default_instance_id, 'coupling_test6.json')

    ! define local component
    CALL yac_fdef_comp(default_instance_id, comp_B, comp_id)

    ! get communicator that contains both components
    CALL yac_fget_comps_comm( &
      default_instance_id, (/comp_A, comp_B/), 2, pair_comm)

    ! check the pair_comm
    CALL check_pair_comm(pair_comm, 'B')
    CALL MPI_Comm_free(pair_comm, ierror)

    !   define local grid
    !
    !   4-------5
    !   |       |
    !   |   o   |
    !   |       |
    !   2-------3
    !   |       |
    !   |   o   |
    !   |       |
    !   0-------1
    CALL yac_fdef_grid( &
      'grid_B', nbr_vertices, cyclic, x_vertices, y_vertices, grid_id_B)

    ! center points in cells (needed e.g. for nearest neighbour)
    CALL yac_fdef_points( &
      grid_id_B, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, point_id_B)

    ! define field
    CALL yac_fdef_field( &
      'A_to_B_tgt', comp_id, (/point_id_B/), 1, 1, "1", YAC_TIME_UNIT_SECOND, &
      field_id)
    CALL yac_fdef_field( &
      'A_to_B_multi', comp_id, (/point_id_B/), 1, 1, "1", YAC_TIME_UNIT_SECOND, &
      multi_field_id)

    ! generate coupling
    CALL yac_fenddef(default_instance_id)

    ! move data from comp_A to comp_B
    recv_field_B = DUMMY_VALUE
    recv_field_ptr(1)%p => recv_field_B(:,1)
    CALL yac_fget( &
      multi_field_id, collection_size, recv_field_ptr, info, ierror)
    CALL test(ALL(recv_field_B == -2.0))
    CALL yac_fget( &
      field_id, nbr_hor_points, collection_size, recv_field_B, info, ierror)

    !-----------------------------------------------------
    ! run internal instance
    !-----------------------------------------------------

      CALL yac_fget_comp_comm(comp_id, internal_comm)

      ! initialise internal YAC instance
      CALL yac_finit_comm(internal_comm, internal_instance_id)
      CALL yac_fread_config_yaml( &
        internal_instance_id, 'coupling_test6_local.yaml')

      ! define internal component
      CALL yac_fdef_comp(internal_instance_id, comp_B, internal_comp_id)

      !   define local grid
      !
      !   4-------5
      !   |       |
      !   |   o   |
      !   |       |
      !   2-------3
      !   |       |
      !   |   o   |
      !   |       |
      !   0-------1
      CALL yac_fdef_grid(               &
        'grid_C', nbr_vertices, cyclic, &
        internal_x_vertices, internal_y_vertices, grid_id_C)

      ! vertices points (needed e.g. for nearest neighbour)
      CALL yac_fdef_points(                           &
        grid_id_C, nbr_vertices, YAC_LOCATION_CORNER, &
        internal_x_vertices, internal_y_vertices, point_id_C)

      ! define field
      CALL yac_fdef_field( &
        'B_to_C', internal_comp_id, (/point_id_B/), 1, 1, "1", &
        YAC_TIME_UNIT_SECOND, field_id_src)
      CALL yac_fdef_field( &
        'B_to_C', internal_comp_id, (/point_id_C/), 1, 1, "1", &
        YAC_TIME_UNIT_SECOND, field_id_tgt)

      ! generate internal coupling
      CALL yac_fenddef(internal_instance_id)

      ! move data from grid B to grid C
      recv_field_C = DUMMY_VALUE
      CALL yac_fexchange( &
        field_id_src, field_id_tgt, nbr_hor_points, 6, collection_size,&
        RESHAPE( &
          recv_field_B, (/nbr_hor_points, collection_size/)), &
        recv_field_C, send_info, recv_info, ierror)

      CALL test(ALL(RESHAPE(recv_field_C, (/6/)) == (/4.0, 4.0, 4.0, 4.0, 3.0, 3.0/)))

      CALL yac_fcleanup(internal_instance_id)

      CALL MPI_Comm_free(internal_comm, ierror)

    !-----------------------------------------------------

    CALL yac_ffinalize(default_instance_id)

  END SUBROUTINE

END PROGRAM
