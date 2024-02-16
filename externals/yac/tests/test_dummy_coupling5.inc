!
! @file test_dummy_coupling5_dble.f90
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
  USE mpi
  USE, INTRINSIC :: ieee_arithmetic

  IMPLICIT NONE

  INTEGER :: instance_id, comp_id, grid_ids(2), point_ids(2), field_ids(2,2), &
             dummy_field_ids(2)
  INTEGER :: global_rank, global_size
  LOGICAL :: is_target, is_source

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)   !< single precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307) !< double precision
  INTEGER, PARAMETER :: wp = TEST_PRECISION              !< selected working precision

  REAL(wp), PARAMETER :: YAC_RAD = 0.017453292519943 ! M_PI / 180.0
  REAL(wp) :: x_vertices(3,2)
  REAL(wp) :: y_vertices(3)

  INTEGER :: cell_global_index(4,2)
  INTEGER :: corner_global_index(9,2)
  INTEGER :: edge_global_index(12,2)
  INTEGER :: cell_core_mask(4,2)
  INTEGER :: corner_core_mask(9,2)
  INTEGER :: edge_core_mask(12,2)
  INTEGER :: corner_field_mask(9,2)

  REAL(wp), TARGET :: source_data(9,2,10)
  REAL(wp) :: ref_recv_field(9,2,7,2)
  INTEGER :: ref_send_info(10,2)
  INTEGER :: ref_recv_info(7)

  INTEGER, PARAMETER :: source_period = 2
  INTEGER, PARAMETER :: target_period = 3
  INTEGER, PARAMETER :: test_time = 18
  INTEGER, PARAMETER :: collection_size = 1

  INTEGER :: i, t, ierror

  LOGICAL :: is_source_timestep, is_target_timestep
  LOGICAL :: is_even_timestep
  TYPE(YAC_PTR_TYPE) :: send_field(1, collection_size)
  TYPE(YAC_PTR_TYPE) :: recv_field(collection_size)
  REAL(wp), TARGET :: recv_field_data(9)
  INTEGER :: send_info, recv_info
  INTEGER :: send_action, recv_action
  INTEGER :: dummy_send_info, dummy_recv_info
  INTEGER :: dummy_send_action, dummy_recv_action
  INTEGER :: dummy_send_ierror, dummy_recv_ierror

  INTEGER :: put_idx, get_idx

  INTEGER :: with_core_mask, with_field_mask

  ! ===================================================================

  CALL start_test('dummy_coupling5')

  CALL MPI_Init (ierror)

  DO with_core_mask = 0, 1
    DO with_field_mask = 0, 1

      CALL yac_finit_instance(instance_id)
      CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
      CALL yac_fread_config_yaml(instance_id, "coupling_test5.yaml")

      CALL MPI_Comm_rank(MPI_COMM_WORLD, global_rank, ierror)
      CALL MPI_Comm_size(MPI_COMM_WORLD, global_size, ierror)

      IF (global_size /= 3) THEN
        WRITE ( * , * ) "Wrong number of processes (should be 3)"
        CALL error_exit
      ENDIF

      is_target = global_rank < 2
      is_source = global_rank > 0

      ! define local component
      CALL yac_fdef_comp(instance_id, 'main_comp', comp_id)

      x_vertices(:,1) = (/0.0, 1.0, 2.0/)
      x_vertices(:,2) = (/1.0, 2.0, 3.0/)
      y_vertices(:) =   (/0.0, 1.0, 2.0/)

      x_vertices = x_vertices * YAC_RAD
      y_vertices = y_vertices * YAC_RAD

      ! define grid (both components use an identical grids)
      IF (is_target) THEN
        CALL yac_fdef_grid(                &
          'target_grid', (/3,3/), (/0,0/), &
          x_vertices(:,global_rank+1), y_vertices, grid_ids(1))
      END IF
      IF (is_source) THEN
        CALL yac_fdef_grid(                &
          'source_grid', (/3,3/), (/0,0/), &
          x_vertices(:,global_rank), y_vertices, grid_ids(2))
      END IF

      cell_global_index = RESHAPE((/0,1,3,4, 1,2,4,5/), SHAPE(cell_global_index))
      corner_global_index =                                &
        RESHAPE(                                           &
          (/0,1,2, 4,5,6, 8,9,10, 1,2,3, 5,6,7, 9,10,11/), &
          SHAPE(corner_global_index))
      edge_global_index =                                                     &
        RESHAPE(                                                              &
          (/0,1,2,3,5, 7,8,9,10,12, 14,15, 2,3,4,5,6, 9,10,11,12,13, 15,16/), &
          SHAPE(edge_global_index))
      cell_core_mask = RESHAPE((/1,1,1,0, 0,1,1,1/), SHAPE(cell_core_mask))
      corner_core_mask = &
        RESHAPE(         &
          (/1,1,0, 0,1,0, 1,1,0, 0,1,1, 0,1,1, 0,1,1/), SHAPE(corner_core_mask))
      edge_core_mask = &
        RESHAPE(       &
          (/1,1,1,1,0, 1,1,1,1,0, 1,0, 0,0,1,1,1, 0,0,1,1,1, 1,1/), &
          SHAPE(edge_core_mask))
      corner_field_mask = &
        RESHAPE(         &
          (/1,1,1, 0,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1/), SHAPE(corner_field_mask))

      IF (is_target) THEN
        CALL yac_fset_global_index( &
          cell_global_index(:,global_rank+1), YAC_LOCATION_CELL, grid_ids(1))
        CALL yac_fset_global_index( &
          corner_global_index(:,global_rank+1), YAC_LOCATION_CORNER, grid_ids(1))
        CALL yac_fset_global_index( &
          edge_global_index(:,global_rank+1), YAC_LOCATION_EDGE, grid_ids(1))
        IF (with_core_mask == 1) THEN
          CALL yac_fset_core_mask( &
            cell_core_mask(:,global_rank+1), YAC_LOCATION_CELL, grid_ids(1))
          CALL yac_fset_core_mask( &
            corner_core_mask(:,global_rank+1), YAC_LOCATION_CORNER, grid_ids(1))
          CALL yac_fset_core_mask( &
            edge_core_mask(:,global_rank+1), YAC_LOCATION_EDGE, grid_ids(1))
        END IF
      END IF
      IF (is_source) THEN
        CALL yac_fset_global_index( &
          cell_global_index(:,global_rank), YAC_LOCATION_CELL, grid_ids(2))
        CALL yac_fset_global_index( &
          corner_global_index(:,global_rank), YAC_LOCATION_CORNER, grid_ids(2))
        CALL yac_fset_global_index( &
          edge_global_index(:,global_rank), YAC_LOCATION_EDGE, grid_ids(2))
        IF (with_core_mask == 1) THEN
          CALL yac_fset_core_mask( &
            cell_core_mask(:,global_rank), YAC_LOCATION_CELL, grid_ids(2))
          CALL yac_fset_core_mask( &
            corner_core_mask(:,global_rank), YAC_LOCATION_CORNER, grid_ids(2))
          CALL yac_fset_core_mask( &
            edge_core_mask(:,global_rank), YAC_LOCATION_EDGE, grid_ids(2))
        END IF
      END IF

      ! define points at the vertices of the grid
      IF (is_target) THEN
        CALL yac_fdef_points( &
          grid_ids(1), (/3,3/), YAC_LOCATION_CORNER, &
          x_vertices(:,global_rank+1), y_vertices, point_ids(1))
        IF (with_field_mask == 1) THEN
          CALL yac_fset_mask(corner_field_mask(:,global_rank+1), point_ids(1))
        END IF
      END IF
      IF (is_source) THEN
        CALL yac_fdef_points( &
          grid_ids(2), (/3,3/), YAC_LOCATION_CORNER, &
          x_vertices(:,global_rank), y_vertices, point_ids(2))
        IF (with_field_mask == 1) THEN
          CALL yac_fset_mask(corner_field_mask(:,global_rank), point_ids(2))
        END IF
      END IF

      ! define fields
      IF (is_target) THEN
        CALL yac_fdef_field( &
          "field_a", comp_id, point_ids(1), 1, 1, "3", YAC_TIME_UNIT_SECOND, &
          field_ids(1,1))
        CALL yac_fdef_field( &
          "field_b", comp_id, point_ids(1), 1, 1, "3", YAC_TIME_UNIT_SECOND, &
          field_ids(1,2))
        CALL yac_fdef_field( &
          "dummy_field", comp_id, point_ids(1), 1, 1, "3", YAC_TIME_UNIT_SECOND, &
          dummy_field_ids(1))
      END IF
      IF (is_source) THEN
        CALL yac_fdef_field( &
          "field_a", comp_id, point_ids(2), 1, 1, "2", YAC_TIME_UNIT_SECOND, &
          field_ids(2,1))
        CALL yac_fdef_field( &
          "field_b", comp_id, point_ids(2), 1, 1, "2", YAC_TIME_UNIT_SECOND, &
          field_ids(2,2))
        CALL yac_fdef_field( &
          "dummy_field", comp_id, point_ids(2), 1, 1, "2", YAC_TIME_UNIT_SECOND, &
          dummy_field_ids(2))
      END IF

      CALL yac_fenddef(instance_id)

      source_data =                                                  &
        RESHAPE(                                                     &
          (/0, 1, 2, 4, 5, 6, 8, 9,10,  1, 2, 3, 5, 6, 7, 9,10,11,   &
            1, 2, 3, 5, 6, 7, 9,10,11,  2, 3, 4, 6, 7, 8,10,11,12,   &
            2, 3, 4, 6, 7, 8,10,11,12,  3, 4, 5, 7, 8, 9,11,12,13,   &
            3, 4, 5, 7, 8, 9,11,12,13,  4, 5, 6, 8, 9,10,12,13,14,   &
            4, 5, 6, 8, 9,10,12,13,14,  5, 6, 7, 9,10,11,13,14,15,   &
            5, 6, 7, 9,10,11,13,14,15,  6, 7, 8,10,11,12,14,15,16,   &
            6, 7, 8,10,11,12,14,15,16,  7, 8, 9,11,12,13,15,16,17,   &
           -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,   &
           -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,   &
           -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1/), &
          SHAPE(source_data))
      ref_recv_field =                                                &
        RESHAPE(                                                      &
          (/ 0, 1, 2, 4, 5, 6, 8, 9,10,  1, 2, 3, 5, 6, 7, 9,10,11,   &
            -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,   &
             6, 9,12,18,21,24,30,33,36,  9,12,15,21,24,27,33,36,39,   &
            -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,   &
            15,18,21,27,30,33,39,42,45, 18,21,24,30,33,36,42,45,48,   &
            -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,   &
            -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,   &
             0, 1, 2, 4, 5, 6, 8, 9,10,  1, 2, 3, 5, 6, 7, 9,10,11,   &
            -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,   &
             3, 4, 5, 7, 8, 9,11,12,13,  4, 5, 6, 8, 9,10,12,13,14,   &
            -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,   &
             6, 7, 8,10,11,12,14,15,16,  7, 8, 9,11,12,13,15,16,17,   &
            -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,   &
            -1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1/), &
          SHAPE(ref_recv_field))

      IF (is_source) THEN
        DO i = 1, 9
          IF (((with_core_mask == 1) .AND. &
               (corner_core_mask(i,global_rank) == 0)) .OR. &
              ((with_field_mask == 1) .AND. &
               (corner_field_mask(i,global_rank) == 0))) THEN
            source_data(i, global_rank, :) = &
              ieee_value(source_data(i, global_rank, :), ieee_signaling_nan)
          END IF
        END DO
      END IF
      IF (is_target) THEN
        DO i = 1, 9
          IF (((with_core_mask == 1) .AND. &
               (corner_core_mask(i,global_rank+1) == 0)) .OR. &
              ((with_field_mask == 1) .AND. &
               (corner_field_mask(i,global_rank+1) == 0))) THEN
            ref_recv_field(i, global_rank+1, :,:) = -1
          END IF
        END DO
      END IF

      ref_send_info =                   &
        RESHAPE(                        &
          (/YAC_ACTION_COUPLING,        &
            YAC_ACTION_REDUCTION,       &
            YAC_ACTION_REDUCTION,       &
            YAC_ACTION_COUPLING,        &
            YAC_ACTION_REDUCTION,       &
            YAC_ACTION_REDUCTION,       &
            YAC_ACTION_PUT_FOR_RESTART, &
            YAC_ACTION_OUT_OF_BOUND,    &
            YAC_ACTION_OUT_OF_BOUND,    &
            YAC_ACTION_OUT_OF_BOUND,    &
            YAC_ACTION_COUPLING,        &
            YAC_ACTION_NONE,            &
            YAC_ACTION_NONE,            &
            YAC_ACTION_COUPLING,        &
            YAC_ACTION_NONE,            &
            YAC_ACTION_NONE,            &
            YAC_ACTION_PUT_FOR_RESTART, &
            YAC_ACTION_OUT_OF_BOUND,    &
            YAC_ACTION_OUT_OF_BOUND,    &
            YAC_ACTION_OUT_OF_BOUND/),  &
          SHAPE(ref_send_info))
      ref_recv_info =      &
        (/YAC_ACTION_COUPLING,        &
          YAC_ACTION_NONE,            &
          YAC_ACTION_COUPLING,        &
          YAC_ACTION_NONE,            &
          YAC_ACTION_GET_FOR_RESTART, &
          YAC_ACTION_OUT_OF_BOUND,    &
          YAC_ACTION_OUT_OF_BOUND/)

      put_idx = 1
      get_idx = 1

      ! do time steps
      DO t = 1, test_time

        is_source_timestep = is_source .AND. (MOD(t, source_period) == 0)
        is_target_timestep = is_target .AND. (MOD(t, target_period) == 0)
        is_even_timestep = MOD(t, 2 * source_period) == 0

        IF (is_source_timestep) THEN
          send_field(1,1)%p => source_data(:, global_rank, put_idx)
        ELSE
          send_field(1,1)%p => NULL()
        END IF
        recv_field_data = -1
        IF (is_target_timestep) THEN
          recv_field(1)%p => recv_field_data
        ELSE
          recv_field(1)%p => NULL()
        END IF

        IF (is_source_timestep) THEN
          CALL yac_fget_action(dummy_field_ids(2), dummy_send_action)
          IF (is_even_timestep) THEN
            CALL yac_fput(                                        &
              dummy_field_ids(2), 1, collection_size, send_field, &
              dummy_send_info, dummy_send_ierror)
            CALL test(dummy_send_ierror == 0)
            CALL test(dummy_send_info == YAC_ACTION_NONE)
            CALL test(dummy_send_action == YAC_ACTION_NONE)
          ELSE
            CALL yac_fupdate(dummy_field_ids(2))
          END IF
        END IF

        IF (is_target_timestep) THEN
          CALL yac_fget_action(dummy_field_ids(1), dummy_recv_action)
          IF (is_even_timestep) THEN
            CALL yac_fget(                                     &
              dummy_field_ids(1), collection_size, recv_field, &
              dummy_recv_info, dummy_recv_ierror)
            CALL test(dummy_recv_ierror == 0)
            CALL test(dummy_recv_info == YAC_ACTION_NONE)
            CALL test(dummy_recv_action == YAC_ACTION_NONE)
          ELSE
            CALL yac_fupdate(dummy_field_ids(1))
          END IF
        END IF

        DO i = 1, 2
          ierror = 0
          IF (is_source_timestep .AND. .NOT. is_target_timestep) THEN
            CALL yac_fget_action(field_ids(2,i), send_action)
            IF (is_even_timestep .AND. send_action == YAC_ACTION_NONE) THEN
              CALL yac_fupdate(field_ids(2,i))
              send_info = YAC_ACTION_NONE
            ELSE
              CALL yac_fput( &
                field_ids(2,i), 1, collection_size, send_field, send_info, ierror)
            END IF
          ELSE IF (.NOT. is_source_timestep .AND. is_target_timestep) THEN
            CALL yac_fget_action(field_ids(1,i), recv_action)
            IF (is_even_timestep .AND. recv_action == YAC_ACTION_NONE) THEN
              CALL yac_fupdate(field_ids(1,i))
              recv_info = YAC_ACTION_NONE
            ELSE
              CALL yac_fget( &
                field_ids(1,i), collection_size, recv_field, recv_info, ierror)
            END IF
          ELSE IF (is_source_timestep .AND. is_target_timestep) THEN
            CALL yac_fget_action(field_ids(2,i), send_action)
            CALL yac_fget_action(field_ids(1,i), recv_action)
            CALL yac_fexchange(                               &
              field_ids(2,i), field_ids(1,i), 1, collection_size, &
              send_field, recv_field, send_info, recv_info, ierror)
          END IF

          CALL test(ierror == 0)
          IF (is_source_timestep) THEN
            CALL test(send_info == ref_send_info(put_idx,i))
            CALL test(send_action == ref_send_info(put_idx,i))
          END IF
          IF (is_target_timestep) THEN
            CALL test(recv_info == ref_recv_info(get_idx))
            CALL test(recv_action == ref_recv_info(get_idx))
            CALL test(ALL(recv_field_data == ref_recv_field(:,global_rank+1,get_idx,i)))
          END IF
        END DO

        IF (is_source_timestep) THEN
          put_idx = put_idx + 1
        END IF
        IF (is_target_timestep) THEN
          get_idx = get_idx + 1
        END IF
      END DO

      CALL yac_ffinalize(instance_id)

    END DO ! with_field_mask
  END DO ! with_core_mask

  CALL MPI_Finalize(ierror)

  CALL stop_test
  CALL exit_tests

! ----------------------------------------------------------

CONTAINS

  SUBROUTINE error_exit ()

    USE mpi, ONLY : mpi_abort, MPI_COMM_WORLD
    USE utest

    INTEGER :: ierror

    CALL test ( .FALSE. )
    CALL stop_test
    CALL exit_tests
    CALL mpi_abort ( MPI_COMM_WORLD, 999, ierror )

  END SUBROUTINE error_exit

END PROGRAM main
