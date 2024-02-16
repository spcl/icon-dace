!
! @file test_dummy_coupling3.inc
!
! @copyright Copyright  (C)  2021 DKRZ, MPI-M
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

#ifndef TEST_PRECISION
#error "TEST_PRECISION is not defined"
#endif

#define DUMMY_VALUE (-1337.0)

#include "test_macros.inc"

! module name has to be specific to each TEST_PRECISION in order
! to avoid conflicts, because the compiler may generate object files
! for each module using the module name as file name
MODULE MOD_NAME

  PUBLIC
  INTEGER, PARAMETER :: wp = TEST_PRECISION              !< selected working precision

  TYPE grid_data

    CHARACTER (LEN=128)       :: name

    INTEGER                   :: nbr_vertices
    INTEGER                   :: nbr_cells
    INTEGER, POINTER          :: num_vertices_per_cell(:)
    INTEGER, POINTER          :: cell_to_vertex(:)

    REAL(wp), ALLOCATABLE :: x_vertices(:)
    REAL(wp), ALLOCATABLE :: y_vertices(:)
    REAL(wp), ALLOCATABLE :: x_cells(:)
    REAL(wp), ALLOCATABLE :: y_cells(:)

    INTEGER, POINTER          :: cell_mask (:)
    INTEGER, POINTER          :: global_cell_id (:)
    INTEGER, POINTER          :: cell_core_mask (:)
    INTEGER, POINTER          :: global_corner_id (:)
    INTEGER, POINTER          :: corner_core_mask (:)

  END TYPE grid_data

  TYPE  intra_comp_yac_info

    LOGICAL                   :: is_comp_root
    INTEGER                   :: yac_id
    INTEGER                   :: out_field_id
    INTEGER                   :: in_field_id

  END TYPE intra_comp_yac_info

  TYPE comp_pair_yac_info

    INTEGER                   :: yac_id
    INTEGER                   :: nbr_cells(2)
    INTEGER                   :: out_field_ids(2)
    INTEGER                   :: in_field_ids(2)

  END TYPE comp_pair_yac_info

END MODULE MOD_NAME

!
! This tests checks the support for multiple parallel YAC instances.
! There are three main components and one dummy component. The dummy component
! does not take part in any coupling.
!
!  All four components are registered in the default YAC instance.
!
!  The three active components couple to each other.
!
!  Each active component shares one process with each of the two other active
!  components.
!
!  Each active component pair generates their own YAC instance.
!
!  Each active component has its own YAC instance in which one of its processes
!  generates its own component.
!
!  The whole setup is run twice in order to simulate the restarting of YAC
!  instances.
!

PROGRAM main

  USE mpi
  USE MOD_NAME
  USE mo_yac_finterface
  USE yaxt
  USE utest

  IMPLICIT NONE

  CHARACTER (LEN=6)               :: global_comp_names(4)
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: local_comp_names(2)
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: comp_pair_names(2)
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: timestep_string

  INTEGER                         :: i, j, jj, k, t
  INTEGER                         :: ierror
  INTEGER                         :: info
  INTEGER                         :: comp_idx
  INTEGER                         :: field_idx
  INTEGER                         :: out_info, in_info
  INTEGER                         :: global_rank
  INTEGER                         :: global_size

  INTEGER                         :: num_local_comps = 0

  INTEGER                         :: default_comp_ids(2)
  INTEGER                         :: default_num_comps = 0
  INTEGER                         :: default_comp_comms(2)

  INTEGER                         :: grid_ids(2)
  INTEGER                         :: cell_point_ids(2)
  INTEGER                         :: regional_grid_id      = -1
  INTEGER                         :: regional_point_id     = -1
  INTEGER                         :: nbr_vertices_regional = 0

  INTEGER                         :: num_vertices_regional(2)
  INTEGER                         :: cyclic(2)

  INTEGER                         :: default_in_field_ids(2,2)
  INTEGER                         :: default_out_field_ids(2)
  INTEGER                         :: local_comp_idx(2)
  INTEGER                         :: nbr_cells(2)
  INTEGER                         :: id
  INTEGER                         :: num_cells
  INTEGER                         :: grid_id, out_id, in_id

  LOGICAL                         :: is_dummy
  REAL(wp), PARAMETER             :: YAC_RAD = 0.017453292519943 ! M_PI / 180.0
  REAL(wp)                        :: regional_coord_x(21)
  REAL(wp)                        :: regional_coord_y(21)

  TYPE(grid_data)                 :: comp_grid_data(2)
  TYPE(intra_comp_yac_info)       :: yac_intra_infos(2)
  TYPE(comp_pair_yac_info)        :: yac_comp_pair_infos(3)

  REAL(wp), ALLOCATABLE           :: out_field_data(:,:,:)
  REAL(wp), ALLOCATABLE           :: in_field_data(:,:)
  REAL(wp), ALLOCATABLE           :: regional_data(:,:)

  INTEGER, PARAMETER              :: collection_size = 4
  INTEGER                         :: temp_collection_size

  TYPE field_data
    REAL(wp), ALLOCATABLE :: data(:,:,:)
  END TYPE field_data

  TYPE(field_data) :: out_field(2)
  TYPE(field_data) :: in_field(2)

  ! ===================================================================

  CALL start_test('dummy_coupling3')
  CALL MPI_Init(ierror)
  CALL xt_initialize(MPI_COMM_WORLD)

  CALL MPI_Comm_rank(MPI_COMM_WORLD, global_rank, ierror)
  CALL MPI_Comm_size(MPI_COMM_WORLD, global_size, ierror)

  IF (global_size /= 10) THEN
    WRITE ( * , * ) "Wrong number of processes (should be 10)"
    CALL error_exit
  ENDIF

  ! ranks 8 and 9 do not initialise YAC, but still need to call the
  ! respective dummy initilialisation and finalise routines
  IF (global_rank > 7) THEN

    CALL yac_finit_dummy()
    CALL xt_finalize()
    CALL MPI_Finalize(ierror)
    CALL stop_test
    CALL exit_tests
    STOP
  END IF

  global_comp_names = (/ 'comp_a', 'comp_b', 'comp_c', 'dummy ' /)

  ! determine the local components

  is_dummy = (global_rank == 3) .OR. (global_rank == 7)

  num_local_comps = 0

  IF ((global_rank == 0) .OR. (global_rank == 4) .OR. (global_rank == 5)) THEN
    num_local_comps = num_local_comps + 1
    local_comp_names(num_local_comps) = global_comp_names(1)
  ENDIF

  IF ((global_rank == 1) .OR. (global_rank == 4) .OR. (global_rank == 6)) THEN
    num_local_comps = num_local_comps + 1
    local_comp_names(num_local_comps) = global_comp_names(2)
  ENDIF

  IF ((global_rank == 2) .OR. (global_rank == 5) .OR. (global_rank == 6)) THEN
    num_local_comps = num_local_comps + 1
    local_comp_names(num_local_comps) = global_comp_names(3)
  ENDIF

  IF ((global_rank == 3) .OR. (global_rank == 7))                         THEN
    num_local_comps = num_local_comps + 1
    local_comp_names(num_local_comps) = global_comp_names(4)
  ENDIF

  ! initialise default YAC instance
  CALL yac_finit ()
  CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
  CALL yac_fread_config_yaml("coupling_test3_default.yaml")

  ! register local components in default YAC instance
  default_num_comps = num_local_comps
  CALL yac_fdef_comps(local_comp_names, default_num_comps, default_comp_ids)

  DO comp_idx = 1, num_local_comps
    CALL yac_fget_comp_comm ( default_comp_ids(comp_idx), default_comp_comms(comp_idx) )
  END DO

  ! generate and register grid, global ids, core mask, and points
  IF ( .NOT. is_dummy ) THEN

    DO comp_idx = 1, num_local_comps

      SELECT CASE (local_comp_names(comp_idx)(6:6))

      CASE ("a")
        CALL generate_icon_grid_data ( default_comp_comms(comp_idx), &
                                       comp_grid_data(comp_idx) )
      CASE ("b")
        CALL generate_cube_grid_data ( default_comp_comms(comp_idx), &
                                       comp_grid_data(comp_idx) )
      CASE ("c")
        CALL generate_mpiom_grid_data ( default_comp_comms(comp_idx), &
                                        comp_grid_data(comp_idx) )
      CASE DEFAULT
        WRITE ( * , * ) ' Unknown component!'
        CALL error_exit ()

      END SELECT

      CALL yac_fdef_grid (                                      &
           TRIM(comp_grid_data(comp_idx)%name),                 &
           comp_grid_data(comp_idx)%nbr_vertices,               &
           comp_grid_data(comp_idx)%nbr_cells,                  &
           SUM(comp_grid_data(comp_idx)%num_vertices_per_cell), &
           comp_grid_data(comp_idx)%num_vertices_per_cell,      &
           comp_grid_data(comp_idx)%x_vertices,                 &
           comp_grid_data(comp_idx)%y_vertices,                 &
           comp_grid_data(comp_idx)%cell_to_vertex,             &
           grid_id )

      grid_ids(comp_idx) = grid_id

      CALL yac_fset_global_index(                          &
           comp_grid_data(comp_idx)%global_corner_id,      &
           YAC_LOCATION_CORNER,                            &
           grid_ids(comp_idx) )

      CALL yac_fset_core_mask (                            &
           comp_grid_data(comp_idx)%corner_core_mask,      &
           YAC_LOCATION_CORNER,                            &
           grid_ids(comp_idx) )

      CALL yac_fset_global_index(                          &
           comp_grid_data(comp_idx)%global_cell_id,        &
           YAC_LOCATION_CELL,                              &
           grid_ids(comp_idx) )

      CALL yac_fset_core_mask(                             &
           comp_grid_data(comp_idx)%cell_core_mask,        &
           YAC_LOCATION_CELL,                              &
           grid_ids(comp_idx) )

      CALL yac_fdef_points (                               &
           grid_ids(comp_idx),                             &
           comp_grid_data(comp_idx)%nbr_cells,             &
           YAC_LOCATION_CELL,                              &
           comp_grid_data(comp_idx)%x_cells,               &
           comp_grid_data(comp_idx)%y_cells,               &
           cell_point_ids(comp_idx) )

      IF ( ASSOCIATED(comp_grid_data(comp_idx)%cell_mask) ) &
           CALL yac_fset_mask (                             &
           comp_grid_data(comp_idx)%cell_mask,              &
           cell_point_ids(comp_idx) )

    END DO ! comp_idx

    ! regional grid on root processes of each component
    IF ( (global_rank == 0) .OR. (global_rank == 1) .OR. (global_rank == 2) ) THEN

      ! register regional grid on rank 0
      regional_coord_x = (/ -20,-18,-16,-14,-12,-10, -8, -6, -4, -2, &
           0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20 /)
      regional_coord_y = (/ -20,-18,-16,-14,-12,-10, -8, -6, -4, -2, &
           0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20 /)

      num_vertices_regional(1) = SIZE(regional_coord_x)
      num_vertices_regional(2) = SIZE(regional_coord_y)

      nbr_vertices_regional = num_vertices_regional(1) * num_vertices_regional(2)

      cyclic = (/ 0, 0 /)
      DO i = 1, num_vertices_regional(1)
        regional_coord_x(i) = regional_coord_x(i) * YAC_RAD
      ENDDO

      DO i = 1, num_vertices_regional(2)
        regional_coord_y(i) = regional_coord_y(i) * YAC_RAD;
      ENDDO

      CALL yac_fdef_grid ( "regional_grid",       &
                           num_vertices_regional, &
                           cyclic,                &
                           regional_coord_x,      &
                           regional_coord_y,      &
                           regional_grid_id )

      CALL yac_fdef_points ( regional_grid_id,      &
                             num_vertices_regional, &
                             YAC_LOCATION_CORNER,   &
                             regional_coord_x,      &
                             regional_coord_y,      &
                             regional_point_id )
    END IF ! global_rank

    ! setup yac instance which couples a component to one of its own ranks

    DO comp_idx = 1, num_local_comps
      yac_intra_infos(comp_idx) =        &
           setup_intra_comp_yac (        &
           default_comp_comms(comp_idx), &
           local_comp_names(comp_idx),   &
           cell_point_ids(comp_idx),     &
           regional_point_id,            &
           collection_size )
    END DO

    ! register in and out fields

    DO comp_idx = 1, num_local_comps
      CALL def_default_fields (              &
           default_comp_ids(comp_idx),       &
           local_comp_names(comp_idx),       &
           cell_point_ids(comp_idx),         &
           default_in_field_ids(:,comp_idx), &
           default_out_field_ids(comp_idx),  &
           collection_size )
    END DO

  ENDIF ! .NOT. is_dummy

  CALL yac_fsync_def()

  CALL test(yac_fget_start_datetime() == '+1800-01-01T00:00:00.000')
  CALL test(yac_fget_end_datetime() == '+2100-01-01T00:00:00.000')

  ! generate interpolations for the default YAC instance

  CALL yac_fenddef ( )

  ! setup yac instances that only couple all component combinations

  k = 1
  DO i = 1, 3

    local_comp_idx(1) = HUGE(i)

    DO comp_idx = 1, num_local_comps
      IF ( TRIM(global_comp_names(i)) == TRIM(local_comp_names(comp_idx)) ) &
           local_comp_idx(1) = comp_idx
    END DO

    DO j = i + 1, 3
      local_comp_idx(2) = HUGE(i)
      DO comp_idx = 1, num_local_comps
        IF ( TRIM(global_comp_names(j)) == TRIM(local_comp_names(comp_idx)) ) &
             local_comp_idx(2) = comp_idx
      END DO

        comp_pair_names(1) = TRIM(global_comp_names(i))
        comp_pair_names(2) = TRIM(global_comp_names(j))

        nbr_cells = (/ comp_grid_data(1)%nbr_cells, comp_grid_data(2)%nbr_cells /)

        yac_comp_pair_infos(k) =                  &
           setup_comp_pair_yac ( comp_pair_names, &
                                 cell_point_ids,  &
                                 nbr_cells,       &
                                 local_comp_idx,  &
                                 collection_size )

      k = k + 1

    END DO ! j-loop

  END DO ! i-loop

  ! if this is not the dummy process, which is not involved in any coupling

  IF ( .NOT. is_dummy ) THEN

    DO comp_idx = 1, num_local_comps

      timestep_string = yac_fget_timestep_from_field_id( &
           default_out_field_ids(comp_idx));
      CALL test(TRIM(timestep_string) == "PT01.000S")
      temp_collection_size = &
        yac_fget_collection_size_from_field_id(default_out_field_ids(comp_idx))
      CALL test(temp_collection_size == collection_size)

    END DO

    DO comp_idx = 1, num_local_comps
      ALLOCATE(                                                  &
        in_field(comp_idx)%data(                                 &
          comp_grid_data(comp_idx)%nbr_cells,1,collection_size), &
        out_field(comp_idx)%data(                                &
          comp_grid_data(comp_idx)%nbr_cells,1,collection_size))
      DO i = 1, comp_grid_data(comp_idx)%nbr_cells
        IF ( comp_grid_data(comp_idx)%cell_core_mask(i) == 1 ) THEN
          out_field(comp_idx)%data(i,:,:) = &
            test_harmonic(comp_grid_data(comp_idx)%x_cells(i), &
                          comp_grid_data(comp_idx)%y_cells(i))
        ELSE
          out_field(comp_idx)%data(i,:,:) = DUMMY_VALUE
        END IF
      END DO
    END DO

    ! do some ping-pongs in the default yac instance

    DO t = 1, 10

      DO comp_idx = 1, num_local_comps
        id = default_out_field_ids(comp_idx)
        CALL yac_fput(                                                &
          id, comp_grid_data(comp_idx)%nbr_cells, 1, collection_size, &
          out_field(comp_idx)%data, info, ierror )
      END DO

      DO comp_idx = 1, num_local_comps
        DO field_idx = 1, 2
          id = default_in_field_ids(field_idx, comp_idx)
          in_field(comp_idx)%data(1:comp_grid_data(comp_idx)%nbr_cells,:,:) = &
            DUMMY_VALUE
          CALL yac_fget(                                             &
            id, comp_grid_data(comp_idx)%nbr_cells, collection_size, &
            in_field(comp_idx)%data(:,1,:), info, ierror )
        END DO
      END DO

    END DO ! t-loop

    ! interpolate received results to the regional grid defined on
    ! rank 0 of the current component

    DO comp_idx = 1, num_local_comps

      IF ( yac_intra_infos(comp_idx)%is_comp_root ) THEN

        out_id = yac_intra_infos(comp_idx)%out_field_id
        in_id = yac_intra_infos(comp_idx)%in_field_id

        ALLOCATE ( regional_data(nbr_vertices_regional,collection_size) )
        regional_data = DUMMY_VALUE;

        CALL yac_fexchange (                            &
                    out_id,                             &
                    in_id,                              &
                    comp_grid_data(comp_idx)%nbr_cells, &
                    1,                                  &
                    nbr_vertices_regional,              &
                    collection_size,                    &
                    in_field(comp_idx)%data,            &
                    regional_data,                      &
                    out_info, in_info, ierror )

       DEALLOCATE (regional_data)

      ELSE

        id = yac_intra_infos(comp_idx)%out_field_id;
        CALL yac_fput( &
          id, comp_grid_data(comp_idx)%nbr_cells, 1, collection_size, &
          in_field(comp_idx)%data, info, ierror )

      END IF

    END DO

    ! clean-up

    DO comp_idx = 1, num_local_comps
      DEALLOCATE(out_field(comp_idx)%data)
      DEALLOCATE(in_field(comp_idx)%data)
      CALL grid_data_free ( comp_grid_data(comp_idx) )
    END DO

    ! do some exchanges between within the component pair YAC instances
    ! for all available component pairs

    DO i = 1,  3;

      ! if the local process is not a part of the current component pair
      IF ( yac_comp_pair_infos(i)%yac_id == HUGE(i) ) CYCLE

      ! if the local process has both components of the current pair
      IF ( (yac_comp_pair_infos(i)%nbr_cells(1) > 0) .AND. &
           (yac_comp_pair_infos(i)%nbr_cells(2) > 0) )  THEN

        DO j = 1, 2
          num_cells = yac_comp_pair_infos(i)%nbr_cells(j)
          ALLOCATE(out_field(j)%data(num_cells,1,collection_size), &
                   in_field(j)%data(num_cells,1,collection_size))
          out_field(j)%data(1:num_cells,:,:) = 0.0;
          in_field(j)%data(1:num_cells,:,:) = DUMMY_VALUE
        END DO

        DO j = 1, 2

          IF ( j == 1 ) jj = 2
          IF ( j == 2 ) jj = 1
          in_id = yac_comp_pair_infos(i)%in_field_ids(jj)
          out_id = yac_comp_pair_infos(i)%out_field_ids(j)

          CALL yac_fexchange ( out_id,                               &
                               in_id,                                &
                               yac_comp_pair_infos(i)%nbr_cells(j),  &
                               1,                                    &
                               yac_comp_pair_infos(i)%nbr_cells(jj), &
                               collection_size,                      &
                               out_field(j)%data,                    &
                               in_field(jj)%data(:,1,:),             &
                               out_info, in_info, ierror )
        END DO

        DO j = 1, 2
          DEALLOCATE(out_field(j)%data)
          DEALLOCATE(in_field(j)%data)
        END DO

      ELSE

        DO j = 1, 2

          ! if the source component is available locally
          IF (yac_comp_pair_infos(i)%nbr_cells(j) > 0) THEN

            num_cells = yac_comp_pair_infos(i)%nbr_cells(j)

            ! send some dummy data
            ALLOCATE(out_field_data(num_cells,1,collection_size))
            out_field_data = 0.0
            id = yac_comp_pair_infos(i)%out_field_ids(j)
            CALL yac_fput ( id, num_cells, 1, collection_size, &
                            out_field_data, info, ierror )
            DEALLOCATE(out_field_data)

          ELSE

            num_cells = yac_comp_pair_infos(i)%nbr_cells(MERGE(1, 2, j == 2))

            ! receive data
            ALLOCATE(in_field_data(num_cells,collection_size))
            id = yac_comp_pair_infos(i)%in_field_ids(MERGE(1, 2, j == 2))
            in_field_data = DUMMY_VALUE
            CALL yac_fget ( id, num_cells, collection_size, in_field_data, &
                            info, ierror )
            DEALLOCATE(in_field_data)

          END IF

        END DO

      END IF

    END DO ! i-loop

  END IF ! dummy

  DO comp_idx = 1, num_local_comps
    CALL MPI_Comm_free(default_comp_comms(comp_idx), ierror)
  END DO

  IF ( .NOT. is_dummy ) THEN
    DO i = 1, 3
      IF ( yac_comp_pair_infos(i)%yac_id /= HUGE(i) ) THEN
        CALL yac_fcleanup_instance ( yac_comp_pair_infos(i)%yac_id )
      END IF
    END DO
    DO comp_idx = 1, num_local_comps
      CALL yac_fcleanup_instance ( yac_intra_infos(comp_idx)%yac_id )
    END DO
  END IF

  CALL yac_ffinalize ()

  CALL xt_finalize ()
  CALL MPI_Finalize (ierror)
  CALL stop_test
  CALL exit_tests

CONTAINS

  FUNCTION test_harmonic ( lon, lat )
    USE MOD_NAME, ONLY: wp
    real(wp) :: lon, lat
    real(wp) :: test_harmonic
    test_harmonic = 2.0 + ((sin(2.0 * lat))**16.0  * cos(16.0 * lon))
  END FUNCTION test_harmonic

  FUNCTION setup_intra_comp_yac ( comp_comm,         &
                                comp_name,         &
                                base_point_id,     &
                                regional_point_id, &
                                collection_size) RESULT ( yac_intra_info )

  USE mo_yac_finterface, ONLY  : yac_finit_comm_instance,    &
                                 yac_fdef_comps_instance,    &
                                 yac_fdef_field,             &
                                 yac_fget_start_datetime,    &
                                 yac_fget_end_datetime,      &
                                 yac_fenddef_instance,       &
                                 yac_fget_timestep_from_field_id,    &
                                 yac_fget_collection_size_from_field_id,   &
                                 yac_fsync_def,              &
                                 yac_fread_config_yaml,      &
                                 YAC_MAX_CHARLEN,            &
                                 YAC_TIME_UNIT_SECOND

  USE MOD_NAME, ONLY : intra_comp_yac_info, wp
  USE utest

  IMPLICIT NONE

  INTEGER           :: comp_comm
  CHARACTER (LEN=*) :: comp_name
  INTEGER           :: base_point_id
  INTEGER           :: regional_point_id
  INTEGER           :: collection_size

  ! Return Type
  TYPE(intra_comp_yac_info) :: yac_intra_info

  ! initialise intra component YAC instance
  INTEGER yac_id
  CHARACTER (LEN=128)             :: yaml_file_name
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: comp_names(2)
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: timestep_string
  INTEGER                         :: comp_ids(2)
  INTEGER                         :: field_ids(2)

  INTEGER                         :: point_ids(1)
  INTEGER                         :: temp_collection_size

  WRITE ( yaml_file_name , '(A15,A6,A11)' ) 'coupling_test3_', TRIM(comp_name), '_intra.yaml'
  CALL yac_finit_comm_instance ( comp_comm, yac_id)
  CALL yac_fread_config_yaml(yac_id, TRIM(yaml_file_name))

  ! register component(s) (all processes are part of the base component
  !                        process 0 has an additional component)

  comp_names(1) = TRIM(comp_name)
  comp_names(2) = "comp_regional"

  IF ( regional_point_id /= -1 ) THEN
    CALL yac_fdef_comps_instance ( yac_id, comp_names, 2, comp_ids )
  ELSE
    CALL yac_fdef_comp_instance ( yac_id, comp_names(1), comp_ids(1) )
  ENDIF


  ! register field (from base component to regional one)

  point_ids(1) = base_point_id
  CALL yac_fdef_field ( "base_to_regional",   &
                        comp_ids(1),          &
                        point_ids,            &
                        1,                    &
                        collection_size,      &
                        "1",                  &
                        YAC_TIME_UNIT_SECOND, &
                        field_ids(1) )

  IF (regional_point_id /= -1) THEN
    point_ids(1) = regional_point_id
    CALL yac_fdef_field( "base_to_regional",   &
                         comp_ids(2),          &
                         point_ids,            &
                         1,                    &
                         collection_size,      &
                         "1",                  &
                         YAC_TIME_UNIT_SECOND, &
                         field_ids(2) )
  ENDIF

  CALL yac_fsync_def(yac_id)

  CALL test(yac_fget_start_datetime(yac_id) == '+1800-01-01T00:00:00.000')
  CALL test(yac_fget_end_datetime(yac_id) == '+2100-01-01T00:00:00.000')

  ! generate interpolations for the default YAC instance

  CALL yac_fenddef_instance ( yac_id )

  yac_intra_info%yac_id       = yac_id
  yac_intra_info%is_comp_root = regional_point_id /= -1
  yac_intra_info%out_field_id = field_ids(1)
  yac_intra_info%in_field_id  = field_ids(2)

  timestep_string = yac_fget_timestep_from_field_id(field_ids(1));
  CALL test(TRIM(timestep_string) == "PT01.000S")
  temp_collection_size = yac_fget_collection_size_from_field_id(field_ids(1))
  CALL test(temp_collection_size == collection_size)

END FUNCTION setup_intra_comp_yac

! ----------------------------------------------------------

FUNCTION setup_comp_pair_yac ( comp_pair_names,    &
                               point_ids,          &
                               nbr_cells,          &
                               local_comp_indices, &
                               collection_size ) RESULT ( yac_info )

  USE mo_yac_finterface, ONLY  : yac_fget_comps_comm,         &
                                 yac_finit_comm_instance,     &
                                 yac_fdef_field,              &
                                 yac_fdef_comps_instance,     &
                                 yac_fenddef_instance,        &
                                 yac_fget_role_from_field_id, &
                                 yac_fread_config_yaml,       &
                                 YAC_MAX_CHARLEN,             &
                                 YAC_TIME_UNIT_SECOND

  USE MOD_NAME, ONLY : comp_pair_yac_info, wp

  IMPLICIT NONE

  CHARACTER (LEN=*)        :: comp_pair_names(2)
  INTEGER                  :: point_ids(2)
  INTEGER                  :: nbr_cells(2)
  INTEGER                  :: local_comp_indices(2)
  INTEGER                  :: collection_size

  ! Return Type
  TYPE(comp_pair_yac_info) :: yac_info

  ! local
  INTEGER                         :: i, j, k
  INTEGER                         :: yac_id
  INTEGER                         :: num_comps
  INTEGER                         :: comp_ids(2)
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: comp_names(2)
  CHARACTER (LEN=16)              :: field_name
  CHARACTER (LEN=256)             :: yaml_file_name
  INTEGER                         :: comp_pair_comm
  INTEGER                         :: field_ids(2,2)
  INTEGER                         :: comp_idx
  INTEGER                         :: ierror
  INTEGER                         :: roles(2)
  INTEGER                         :: imax

  IF ((local_comp_indices(1) == HUGE(imax)) .AND.  &
      (local_comp_indices(2) == HUGE(imax))) THEN
    yac_info%yac_id = HUGE(imax)
    RETURN
  END IF

  ! get the communicator for the component pair
  comp_names(1) = comp_pair_names(1)
  comp_names(2) = comp_pair_names(2)
  CALL yac_fget_comps_comm(comp_names, 2, comp_pair_comm)

  WRITE (yaml_file_name, '(A15,A,A1,A,A5)' ) 'coupling_test3_', &
       TRIM(comp_pair_names(1)), '_', TRIM(comp_pair_names(2)), '.yaml'

  ! initialise component pair YAC instance
  CALL yac_finit_comm_instance( comp_pair_comm, yac_id )
  CALL yac_fread_config_yaml(yac_id, TRIM(yaml_file_name))

  ! register component(s) the local process can have one or two components
  ! associated to the YAC instance
  num_comps = 0
  DO i = 1, 2
    IF ( local_comp_indices(i) /= HUGE(imax) ) THEN
      num_comps = num_comps+1
      comp_names(num_comps) = comp_pair_names(i)
    ENDIF
  END DO
  CALL yac_fdef_comps_instance ( yac_id,     &
                                 comp_names, &
                                 num_comps,  &
                                 comp_ids )

  k = 1

  DO i = 1, 2

    comp_idx = local_comp_indices(i);
    if ( comp_idx == HUGE(i) ) cycle

    ! register field
    DO j = 1, 2
      WRITE ( field_name, '(A5,A1,A4)' ) 'comp_', comp_pair_names(j)(6:6), '_out'
      CALL yac_fdef_field ( TRIM(field_name),     &
                            comp_ids(k),          &
                            point_ids(comp_idx),  &
                            1,                    &
                            collection_size,      &
                            "1",                  &
                            YAC_TIME_UNIT_SECOND, &
                            field_ids(i,j) )
    END DO

    k = k + 1

  END DO

  ! generate interpolations for the component pair YAC instance
  CALL yac_fenddef_instance ( yac_id )

  CALL MPI_Comm_free ( comp_pair_comm, ierror )

  yac_info%yac_id = yac_id

  DO i = 1, 2
    comp_idx = local_comp_indices(i)
    IF ( comp_idx == HUGE(i) ) THEN

      yac_info%nbr_cells(i)     = 0
      yac_info%out_field_ids(i) = HUGE(i)
      yac_info%in_field_ids(i)  = HUGE(i)

    ELSE

      yac_info%nbr_cells(i) = nbr_cells(comp_idx)

      DO j = 1, 2

         roles(j) = yac_fget_role_from_field_id( field_ids(i,j) )

        SELECT CASE (roles(j) )

        CASE(0)

        CASE(1)
          yac_info%out_field_ids(i) = field_ids(i,j)

        CASE(2)
          yac_info%in_field_ids(i) = field_ids(i,j)

        CASE DEFAULT
          WRITE ( * , * ) "invalid role"
          CALL error_exit ()

        END SELECT

      END DO

      IF ( roles(1) == roles (2) ) THEN
        WRITE ( * , * ) "roles should not match"
        CALL error_exit ()
      END IF

    END IF

  END DO

END FUNCTION setup_comp_pair_yac
! ----------------------------------------------------------

SUBROUTINE def_default_fields ( comp_id,       &
                                comp_name,     &
                                cell_point_id, &
                                in_field_ids,  &
                                out_field_id,  &
                                collection_size)

  USE mo_yac_finterface, ONLY : yac_fdef_field, YAC_TIME_UNIT_SECOND

  IMPLICIT NONE

  INTEGER             :: comp_id
  INTEGER             :: cell_point_id
  INTEGER             :: in_field_ids(2)
  INTEGER             :: out_field_id
  CHARACTER (LEN=*)   :: comp_name
  INTEGER             :: collection_size

  INTEGER             :: i
  INTEGER             :: num_in_fields
  INTEGER             :: point_ids(1)

  CHARACTER (LEN=10)  :: field_names(3)
  CHARACTER (LEN=10)  :: curr_field_name

  INTEGER :: field_id

  field_names = (/ "comp_a_out", "comp_b_out", "comp_c_out" /)

  num_in_fields = 0
  DO i = 1, 3

    curr_field_name = TRIM(field_names(i))
    point_ids(1) = cell_point_id
    CALL yac_fdef_field ( TRIM(curr_field_name), &
         comp_id,               &
         point_ids,             &
         1,                     &
         collection_size,       &
         "1",                   &
         YAC_TIME_UNIT_SECOND,  &
         field_id )

    IF ( comp_name(6:6) == curr_field_name(6:6) ) THEN
      out_field_id = field_id
    ELSE
      num_in_fields = num_in_fields+1
      in_field_ids(num_in_fields) = field_id
    END IF

  END DO

END SUBROUTINE def_default_fields

! ----------------------------------------------------------

SUBROUTINE error_exit ()

  USE mpi, ONLY : mpi_abort, MPI_COMM_WORLD
  USE utest

  IMPLICIT NONE

  INTEGER :: ierror

  CALL test ( .FALSE. )
  CALL stop_test
  CALL exit_tests
  CALL mpi_abort ( MPI_COMM_WORLD, 999, ierror )

END SUBROUTINE error_exit

! ----------------------------------------------------------

SUBROUTINE grid_data_free ( grid )

  USE MOD_NAME, ONLY : grid_data, wp
  USE iso_c_binding, ONLY : c_loc

  IMPLICIT NONE

  TYPE(grid_data) :: grid

  INTERFACE

    SUBROUTINE C_FREE ( ptr ) BIND ( c, name='free' )

      USE, INTRINSIC :: iso_c_binding, only : c_ptr

      TYPE ( c_ptr ), INTENT(IN), VALUE :: ptr

    END SUBROUTINE C_FREE

  END INTERFACE

  IF (ASSOCIATED(grid%num_vertices_per_cell)) &
    CALL C_FREE(c_loc(grid%num_vertices_per_cell(1)))
  IF (ASSOCIATED(grid%cell_to_vertex)) &
    CALL C_FREE(c_loc(grid%cell_to_vertex(1)))
  IF (ALLOCATED(grid%x_vertices)) DEALLOCATE(grid%x_vertices)
  IF (ALLOCATED(grid%y_vertices)) DEALLOCATE(grid%y_vertices)
  IF (ALLOCATED(grid%x_cells)) DEALLOCATE(grid%x_cells)
  IF (ALLOCATED(grid%y_cells)) DEALLOCATE(grid%y_cells)
  IF (ASSOCIATED(grid%cell_mask)) &
    CALL C_FREE(c_loc(grid%cell_mask(1)))
  IF (ASSOCIATED(grid%global_cell_id)) &
    CALL C_FREE(c_loc(grid%global_cell_id(1)))
  IF (ASSOCIATED(grid%cell_core_mask)) &
    CALL C_FREE(c_loc(grid%cell_core_mask(1)))
  IF (ASSOCIATED(grid%global_corner_id)) &
    CALL C_FREE(c_loc(grid%global_corner_id(1)))
  IF (ASSOCIATED(grid%corner_core_mask)) &
    CALL C_FREE(c_loc(grid%corner_core_mask(1)))

END SUBROUTINE grid_data_free

! ----------------------------------------------------------

SUBROUTINE generate_icon_grid_data ( comm, grid )

  USE mpi
  USE MOD_NAME, ONLY : grid_data, wp
  USE iso_c_binding, ONLY : c_loc

  IMPLICIT NONE

  INTEGER          :: comm
  TYPE (grid_data) :: grid

  INTEGER          :: comm_rank, comm_size, ierror
  INTEGER          :: i

  DOUBLE PRECISION, POINTER :: x_vertices(:)
  DOUBLE PRECISION, POINTER :: y_vertices(:)
  DOUBLE PRECISION, POINTER :: x_cells(:)
  DOUBLE PRECISION, POINTER :: y_cells(:)

  INTERFACE

    SUBROUTINE C_FREE ( ptr ) BIND ( c, name='free' )

      USE, INTRINSIC :: iso_c_binding, only : c_ptr

      TYPE ( c_ptr ), INTENT(IN), VALUE :: ptr

    END SUBROUTINE C_FREE

  END INTERFACE

  CALL MPI_Comm_rank( comm, comm_rank, ierror )
  CALL MPI_Comm_size( comm, comm_size, ierror )

  grid%name = "icon_grid";

  CALL fread_part_icon_grid_information ( "icon_grid_R02B02.nc", &
    grid%nbr_vertices, grid%nbr_cells,                           &
    grid%num_vertices_per_cell, grid%cell_to_vertex,             &
    x_vertices, y_vertices, x_cells, y_cells,                    &
    grid%global_cell_id, grid%cell_mask,                         &
    grid%cell_core_mask, grid%global_corner_id,                  &
    grid%corner_core_mask, comm_rank, comm_size )

  IF (ASSOCIATED(x_vertices)) THEN
    ALLOCATE(grid%x_vertices(SIZE(x_vertices)))
    grid%x_vertices = REAL(x_vertices, kind=wp)
    CALL C_FREE(c_loc(x_vertices(1)))
  END IF
  IF (ASSOCIATED(y_vertices)) THEN
    ALLOCATE(grid%y_vertices(SIZE(y_vertices)))
    grid%y_vertices = REAL(y_vertices, kind=wp)
    CALL C_FREE(c_loc(y_vertices(1)))
  END IF
  IF (ASSOCIATED(x_cells)) THEN
    ALLOCATE(grid%x_cells(SIZE(x_cells)))
    grid%x_cells = REAL(x_cells, kind=wp)
    CALL C_FREE(c_loc(x_cells(1)))
  END IF
  IF (ASSOCIATED(y_cells)) THEN
    ALLOCATE(grid%y_cells(SIZE(y_cells)))
    grid%y_cells = REAL(y_cells, kind=wp)
    CALL C_FREE(c_loc(y_cells(1)))
  END IF

  DO i = 1, grid%nbr_cells
     IF ( grid%cell_mask(i) == 0 ) THEN
        grid%cell_mask(i) = 1
     END IF
  END DO

END SUBROUTINE generate_icon_grid_data

SUBROUTINE fread_part_icon_grid_information ( grid_file,             &
                                              nbr_vertices,          &
                                              nbr_cells,             &
                                              num_vertices_per_cell, &
                                              cell_to_vertex,        &
                                              x_vertices,            &
                                              y_vertices,            &
                                              x_cells,               &
                                              y_cells,               &
                                              global_cell_ids,       &
                                              cell_mask,             &
                                              cell_core_mask,        &
                                              global_corner_ids,     &
                                              corner_core_mask,      &
                                              comm_rank,             &
                                              comm_size )

  USE iso_c_binding, ONLY: c_int, c_ptr, c_f_pointer, c_null_char

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE read_part_icon_grid_information ( grid_file,              &
                                                  nbr_vertices,           &
                                                  nbr_cells,              &
                                                  num_vertices_per_cell,  &
                                                  cell_to_vertex,         &
                                                  x_vertices, y_vertices, &
                                                  x_cells, y_cells,       &
                                                  global_cell_ids,        &
                                                  cell_mask,              &
                                                  cell_core_mask,         &
                                                  global_corner_ids,      &
                                                  corner_core_mask,       &
                                                  comm_rank,              &
                                                  comm_size )             &
      BIND ( c, name='read_part_icon_grid_information' )

      USE, INTRINSIC :: iso_c_binding, only : c_int, c_char, c_ptr

      CHARACTER ( kind=c_char), DIMENSION(*) :: grid_file
      INTEGER ( kind=c_int )                 :: nbr_vertices
      INTEGER ( kind=c_int )                 :: nbr_cells
      TYPE (c_ptr )                          :: num_vertices_per_cell
      TYPE (c_ptr )                          :: cell_to_vertex
      TYPE (c_ptr )                          :: x_vertices
      TYPE (c_ptr )                          :: y_vertices
      TYPE (c_ptr )                          :: x_cells
      TYPE (c_ptr )                          :: y_cells
      TYPE (c_ptr )                          :: global_cell_ids
      TYPE (c_ptr )                          :: cell_mask
      TYPE (c_ptr )                          :: cell_core_mask
      TYPE (c_ptr )                          :: global_corner_ids
      TYPE (c_ptr )                          :: corner_core_mask
      INTEGER ( kind=c_int ), VALUE          :: comm_rank
      INTEGER ( kind=c_int ), VALUE          :: comm_size

    END SUBROUTINE read_part_icon_grid_information

  END INTERFACE

  CHARACTER (LEN=*)         :: grid_file
  INTEGER ( kind=c_int )    :: nbr_vertices
  INTEGER ( kind=c_int )    :: nbr_cells
  INTEGER, POINTER          :: num_vertices_per_cell(:)
  INTEGER, POINTER          :: cell_to_vertex(:)
  DOUBLE PRECISION, POINTER :: x_vertices(:)
  DOUBLE PRECISION, POINTER :: y_vertices(:)
  DOUBLE PRECISION, POINTER :: x_cells(:)
  DOUBLE PRECISION, POINTER :: y_cells(:)
  INTEGER, POINTER          :: global_cell_ids(:)
  INTEGER, POINTER          :: cell_mask(:)
  INTEGER, POINTER          :: cell_core_mask(:)
  INTEGER, POINTER          :: global_corner_ids(:)
  INTEGER, POINTER          :: corner_core_mask(:)
  INTEGER ( kind=c_int )    :: comm_rank
  INTEGER ( kind=c_int )    :: comm_size

  TYPE(c_ptr)  :: num_vertices_per_cell_
  TYPE(c_ptr)  :: cell_to_vertex_
  TYPE(c_ptr)  :: x_vertices_
  TYPE(c_ptr)  :: y_vertices_
  TYPE(c_ptr)  :: x_cells_
  TYPE(c_ptr)  :: y_cells_
  TYPE(c_ptr)  :: global_cell_ids_
  TYPE(c_ptr)  :: cell_mask_
  TYPE(c_ptr)  :: cell_core_mask_
  TYPE(c_ptr)  :: global_corner_ids_
  TYPE(c_ptr)  :: corner_core_mask_

  CALL read_part_icon_grid_information ( TRIM(grid_file) // c_null_char, &
                                          nbr_vertices,                   &
                                          nbr_cells,                      &
                                          num_vertices_per_cell_,         &
                                          cell_to_vertex_,                &
                                          x_vertices_, y_vertices_,       &
                                          x_cells_, y_cells_,             &
                                          global_cell_ids_,               &
                                          cell_mask_,                     &
                                          cell_core_mask_,                &
                                          global_corner_ids_,             &
                                          corner_core_mask_,              &
                                          comm_rank,                      &
                                          comm_size )

  CALL C_F_POINTER(num_vertices_per_cell_, num_vertices_per_cell, (/nbr_cells/))
  CALL C_F_POINTER(cell_to_vertex_, cell_to_vertex, (/SUM(num_vertices_per_cell)/))
  CALL C_F_POINTER(x_vertices_, x_vertices, (/nbr_vertices/))
  CALL C_F_POINTER(y_vertices_, y_vertices, (/nbr_vertices/))
  CALL C_F_POINTER(x_cells_, x_cells, (/nbr_cells/))
  CALL C_F_POINTER(y_cells_, y_cells, (/nbr_cells/))
  CALL C_F_POINTER(global_cell_ids_, global_cell_ids, (/nbr_cells/))
  CALL C_F_POINTER(cell_mask_, cell_mask, (/nbr_cells/))
  CALL C_F_POINTER(cell_core_mask_, cell_core_mask, (/nbr_cells/))
  CALL C_F_POINTER(global_corner_ids_, global_corner_ids, (/nbr_vertices/))
  CALL C_F_POINTER(corner_core_mask_, corner_core_mask, (/nbr_vertices/))

  cell_to_vertex = cell_to_vertex + 1

END SUBROUTINE fread_part_icon_grid_information

! ----------------------------------------------------------

SUBROUTINE generate_cube_grid_data ( comm, grid )

  USE mpi
  USE MOD_NAME, ONLY : grid_data, wp
  USE iso_c_binding, ONLY : c_loc

  IMPLICIT NONE

  INTEGER          :: comm
  TYPE (grid_data) :: grid

  INTEGER          :: n = 50
  INTEGER          :: comm_rank, comm_size, ierror

  DOUBLE PRECISION, POINTER :: x_vertices(:)
  DOUBLE PRECISION, POINTER :: y_vertices(:)
  DOUBLE PRECISION, POINTER :: x_cells(:)
  DOUBLE PRECISION, POINTER :: y_cells(:)

  INTERFACE

    SUBROUTINE C_FREE ( ptr ) BIND ( c, name='free' )

      USE, INTRINSIC :: iso_c_binding, only : c_ptr

      TYPE ( c_ptr ), INTENT(IN), VALUE :: ptr

    END SUBROUTINE C_FREE

  END INTERFACE

  CALL MPI_Comm_rank( comm, comm_rank, ierror )
  CALL MPI_Comm_size( comm, comm_size, ierror )

  grid%name = "cube_grid"

  CALL fgenerate_part_cube_grid_information ( n,     &
    grid%nbr_vertices, grid%nbr_cells,               &
    grid%num_vertices_per_cell, grid%cell_to_vertex, &
    x_vertices, y_vertices, x_cells, y_cells,        &
    grid%global_cell_id,                             &
    grid%cell_core_mask, grid%global_corner_id,      &
    grid%corner_core_mask, comm_rank, comm_size )

  IF (ASSOCIATED(x_vertices)) THEN
    ALLOCATE(grid%x_vertices(SIZE(x_vertices)))
    grid%x_vertices = REAL(x_vertices, kind=wp)
    CALL C_FREE(c_loc(x_vertices(1)))
  END IF
  IF (ASSOCIATED(y_vertices)) THEN
    ALLOCATE(grid%y_vertices(SIZE(y_vertices)))
    grid%y_vertices = REAL(y_vertices, kind=wp)
    CALL C_FREE(c_loc(y_vertices(1)))
  END IF
  IF (ASSOCIATED(x_cells)) THEN
    ALLOCATE(grid%x_cells(SIZE(x_cells)))
    grid%x_cells = REAL(x_cells, kind=wp)
    CALL C_FREE(c_loc(x_cells(1)))
  END IF
  IF (ASSOCIATED(y_cells)) THEN
    ALLOCATE(grid%y_cells(SIZE(y_cells)))
    grid%y_cells = REAL(y_cells, kind=wp)
    CALL C_FREE(c_loc(y_cells(1)))
  END IF

  grid%cell_mask => NULL()

END SUBROUTINE generate_cube_grid_data

SUBROUTINE fgenerate_part_cube_grid_information ( n,                     &
                                                  nbr_vertices,          &
                                                  nbr_cells,             &
                                                  num_vertices_per_cell, &
                                                  cell_to_vertex,        &
                                                  x_vertices,            &
                                                  y_vertices,            &
                                                  x_cells,               &
                                                  y_cells,               &
                                                  global_cell_ids,       &
                                                  cell_core_mask,        &
                                                  global_corner_ids,     &
                                                  corner_core_mask,      &
                                                  comm_rank,             &
                                                  comm_size )

  USE iso_c_binding, ONLY: c_int, c_ptr, c_f_pointer

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE generate_part_cube_grid_information ( n,                      &
                                                      nbr_vertices,           &
                                                      nbr_cells,              &
                                                      num_vertices_per_cell,  &
                                                      cell_to_vertex,         &
                                                      x_vertices, y_vertices, &
                                                      x_cells, y_cells,       &
                                                      global_cell_ids,        &
                                                      cell_core_mask,         &
                                                      global_corner_ids,      &
                                                      corner_core_mask,       &
                                                      comm_rank,              &
                                                      comm_size )             &
      BIND ( c, name='generate_part_cube_grid_information' )

      USE, INTRINSIC :: iso_c_binding, only : c_int, c_ptr

      INTEGER ( kind=c_int ), VALUE :: n
      INTEGER ( kind=c_int )        :: nbr_vertices
      INTEGER ( kind=c_int )        :: nbr_cells
      TYPE(c_ptr)                   :: num_vertices_per_cell
      TYPE(c_ptr)                   :: cell_to_vertex
      TYPE(c_ptr)                   :: x_vertices
      TYPE(c_ptr)                   :: y_vertices
      TYPE(c_ptr)                   :: x_cells
      TYPE(c_ptr)                   :: y_cells
      TYPE(c_ptr)                   :: global_cell_ids
      TYPE(c_ptr)                   :: cell_core_mask
      TYPE(c_ptr)                   :: global_corner_ids
      TYPE(c_ptr)                   :: corner_core_mask
      INTEGER ( kind=c_int ), VALUE :: comm_rank
      INTEGER ( kind=c_int ), VALUE :: comm_size

    END SUBROUTINE generate_part_cube_grid_information

  END INTERFACE

  INTEGER ( kind=c_int )          :: n
  INTEGER ( kind=c_int )          :: nbr_vertices
  INTEGER ( kind=c_int )          :: nbr_cells
  INTEGER ( kind=c_int ), POINTER :: num_vertices_per_cell(:)
  INTEGER ( kind=c_int ), POINTER :: cell_to_vertex(:)
  DOUBLE PRECISION, POINTER       :: x_vertices(:)
  DOUBLE PRECISION, POINTER       :: y_vertices(:)
  DOUBLE PRECISION, POINTER       :: x_cells(:)
  DOUBLE PRECISION, POINTER       :: y_cells(:)
  INTEGER ( kind=c_int ), POINTER :: global_cell_ids(:)
  INTEGER ( kind=c_int ), POINTER :: cell_core_mask(:)
  INTEGER ( kind=c_int ), POINTER :: global_corner_ids(:)
  INTEGER ( kind=c_int ), POINTER :: corner_core_mask(:)
  INTEGER ( kind=c_int )          :: comm_rank
  INTEGER ( kind=c_int )          :: comm_size

  TYPE(c_ptr) :: num_vertices_per_cell_
  TYPE(c_ptr) :: cell_to_vertex_
  TYPE(c_ptr) :: x_vertices_
  TYPE(c_ptr) :: y_vertices_
  TYPE(c_ptr) :: x_cells_
  TYPE(c_ptr) :: y_cells_
  TYPE(c_ptr) :: global_cell_ids_
  TYPE(c_ptr) :: cell_core_mask_
  TYPE(c_ptr) :: global_corner_ids_
  TYPE(c_ptr) :: corner_core_mask_

  CALL generate_part_cube_grid_information ( n,                        &
                                              nbr_vertices,             &
                                              nbr_cells,                &
                                              num_vertices_per_cell_,   &
                                              cell_to_vertex_,          &
                                              x_vertices_, y_vertices_, &
                                              x_cells_, y_cells_,       &
                                              global_cell_ids_,         &
                                              cell_core_mask_,          &
                                              global_corner_ids_,       &
                                              corner_core_mask_,        &
                                              comm_rank,                &
                                              comm_size )

  CALL C_F_POINTER(num_vertices_per_cell_, num_vertices_per_cell, (/nbr_cells/))
  CALL C_F_POINTER(cell_to_vertex_, cell_to_vertex, (/SUM(num_vertices_per_cell)/))
  CALL C_F_POINTER(x_vertices_, x_vertices, (/nbr_vertices/))
  CALL C_F_POINTER(y_vertices_, y_vertices, (/nbr_vertices/))
  CALL C_F_POINTER(x_cells_, x_cells, (/nbr_cells/))
  CALL C_F_POINTER(y_cells_, y_cells, (/nbr_cells/))
  CALL C_F_POINTER(global_cell_ids_, global_cell_ids, (/nbr_cells/))
  CALL C_F_POINTER(cell_core_mask_, cell_core_mask, (/nbr_cells/))
  CALL C_F_POINTER(global_corner_ids_, global_corner_ids, (/nbr_vertices/))
  CALL C_F_POINTER(corner_core_mask_, corner_core_mask, (/nbr_vertices/))

  cell_to_vertex = cell_to_vertex + 1

END SUBROUTINE fgenerate_part_cube_grid_information

! ----------------------------------------------------------

SUBROUTINE generate_mpiom_grid_data ( comm, grid )

  USE mpi
  USE MOD_NAME, ONLY : grid_data, wp
  USE iso_c_binding, ONLY : c_loc

  IMPLICIT NONE

  INTEGER          :: comm
  TYPE (grid_data) :: grid

  INTEGER          :: i
  INTEGER          :: comm_rank, comm_size, ierror

  DOUBLE PRECISION, POINTER :: x_vertices(:)
  DOUBLE PRECISION, POINTER :: y_vertices(:)
  DOUBLE PRECISION, POINTER :: x_cells(:)
  DOUBLE PRECISION, POINTER :: y_cells(:)

  INTERFACE

    SUBROUTINE C_FREE ( ptr ) BIND ( c, name='free' )

      USE, INTRINSIC :: iso_c_binding, only : c_ptr

      TYPE ( c_ptr ), INTENT(IN), VALUE :: ptr

    END SUBROUTINE C_FREE

  END INTERFACE

  CALL MPI_Comm_rank( comm, comm_rank, ierror )
  CALL MPI_Comm_size( comm, comm_size, ierror )

  grid%name = "mpiom_grid"

  CALL fread_part_mpiom_grid_information ( "GR30_lsm.nc", &
    grid%nbr_vertices, grid%nbr_cells,                    &
    grid%num_vertices_per_cell, grid%cell_to_vertex,      &
    x_vertices, y_vertices, x_cells, y_cells,             &
    grid%global_cell_id, grid%cell_mask,                  &
    grid%cell_core_mask, grid%global_corner_id,           &
    grid%corner_core_mask, comm_rank, comm_size )

  IF (ASSOCIATED(x_vertices)) THEN
    ALLOCATE(grid%x_vertices(SIZE(x_vertices)))
    grid%x_vertices = REAL(x_vertices, kind=wp)
    CALL C_FREE(c_loc(x_vertices(1)))
  END IF
  IF (ASSOCIATED(y_vertices)) THEN
    ALLOCATE(grid%y_vertices(SIZE(y_vertices)))
    grid%y_vertices = REAL(y_vertices, kind=wp)
    CALL C_FREE(c_loc(y_vertices(1)))
  END IF
  IF (ASSOCIATED(x_cells)) THEN
    ALLOCATE(grid%x_cells(SIZE(x_cells)))
    grid%x_cells = REAL(x_cells, kind=wp)
    CALL C_FREE(c_loc(x_cells(1)))
  END IF
  IF (ASSOCIATED(y_cells)) THEN
    ALLOCATE(grid%y_cells(SIZE(y_cells)))
    grid%y_cells = REAL(y_cells, kind=wp)
    CALL C_FREE(c_loc(y_cells(1)))
  END IF

  DO i = 1, grid%nbr_cells
    IF ( grid%cell_mask(i) == 0 ) THEN
      grid%cell_mask(i) = 1
    ELSE
      grid%cell_mask(i) = 0
    END IF
  END DO

END SUBROUTINE generate_mpiom_grid_data

SUBROUTINE fread_part_mpiom_grid_information ( grid_file,              &
                                                nbr_vertices,           &
                                                nbr_cells,              &
                                                num_vertices_per_cell,  &
                                                cell_to_vertex,         &
                                                x_vertices, y_vertices, &
                                                x_cells, y_cells,       &
                                                global_cell_ids,        &
                                                cell_mask,              &
                                                cell_core_mask,         &
                                                global_corner_ids,      &
                                                corner_core_mask,       &
                                                comm_rank,              &
                                                comm_size )

  USE iso_c_binding, ONLY: c_int, c_ptr, c_f_pointer, c_null_char

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE read_part_mpiom_grid_information ( grid_file,              &
                                                  nbr_vertices,           &
                                                  nbr_cells,              &
                                                  num_vertices_per_cell,  &
                                                  cell_to_vertex,         &
                                                  x_vertices, y_vertices, &
                                                  x_cells, y_cells,       &
                                                  global_cell_ids,        &
                                                  cell_mask,              &
                                                  cell_core_mask,         &
                                                  global_corner_ids,      &
                                                  corner_core_mask,       &
                                                  comm_rank,              &
                                                  comm_size )             &
      BIND ( c, name='read_part_mpiom_grid_information' )

      USE, INTRINSIC :: iso_c_binding, only : c_int, c_char, c_ptr

      CHARACTER ( kind=c_char), DIMENSION(*) :: grid_file
      INTEGER ( kind=c_int )                 :: nbr_vertices
      INTEGER ( kind=c_int )                 :: nbr_cells
      TYPE (c_ptr )                          :: num_vertices_per_cell
      TYPE (c_ptr )                          :: cell_to_vertex
      TYPE (c_ptr )                          :: x_vertices
      TYPE (c_ptr )                          :: y_vertices
      TYPE (c_ptr )                          :: x_cells
      TYPE (c_ptr )                          :: y_cells
      TYPE (c_ptr )                          :: global_cell_ids
      TYPE (c_ptr )                          :: cell_mask
      TYPE (c_ptr )                          :: cell_core_mask
      TYPE (c_ptr )                          :: global_corner_ids
      TYPE (c_ptr )                          :: corner_core_mask
      INTEGER ( kind=c_int ), VALUE          :: comm_rank
      INTEGER ( kind=c_int ), VALUE          :: comm_size

    END SUBROUTINE read_part_mpiom_grid_information

  END INTERFACE

  CHARACTER (LEN=*)         :: grid_file
  INTEGER ( kind=c_int )    :: nbr_vertices
  INTEGER ( kind=c_int )    :: nbr_cells
  INTEGER, POINTER          :: num_vertices_per_cell(:)
  INTEGER, POINTER          :: cell_to_vertex(:)
  DOUBLE PRECISION, POINTER :: x_vertices(:)
  DOUBLE PRECISION, POINTER :: y_vertices(:)
  DOUBLE PRECISION, POINTER :: x_cells(:)
  DOUBLE PRECISION, POINTER :: y_cells(:)
  INTEGER, POINTER          :: global_cell_ids(:)
  INTEGER, POINTER          :: cell_mask(:)
  INTEGER, POINTER          :: cell_core_mask(:)
  INTEGER, POINTER          :: global_corner_ids(:)
  INTEGER, POINTER          :: corner_core_mask(:)
  INTEGER ( kind=c_int )    :: comm_rank
  INTEGER ( kind=c_int )    :: comm_size

  TYPE(c_ptr)  :: num_vertices_per_cell_
  TYPE(c_ptr)  :: cell_to_vertex_
  TYPE(c_ptr)  :: x_vertices_
  TYPE(c_ptr)  :: y_vertices_
  TYPE(c_ptr)  :: x_cells_
  TYPE(c_ptr)  :: y_cells_
  TYPE(c_ptr)  :: global_cell_ids_
  TYPE(c_ptr)  :: cell_mask_
  TYPE(c_ptr)  :: cell_core_mask_
  TYPE(c_ptr)  :: global_corner_ids_
  TYPE(c_ptr)  :: corner_core_mask_

  CALL read_part_mpiom_grid_information ( TRIM(grid_file) // c_null_char, &
                                          nbr_vertices,                   &
                                          nbr_cells,                      &
                                          num_vertices_per_cell_,         &
                                          cell_to_vertex_,                &
                                          x_vertices_, y_vertices_,       &
                                          x_cells_, y_cells_,             &
                                          global_cell_ids_,               &
                                          cell_mask_,                     &
                                          cell_core_mask_,                &
                                          global_corner_ids_,             &
                                          corner_core_mask_,              &
                                          comm_rank,                      &
                                          comm_size )

  CALL C_F_POINTER(num_vertices_per_cell_, num_vertices_per_cell, (/nbr_cells/))
  CALL C_F_POINTER(cell_to_vertex_, cell_to_vertex, (/SUM(num_vertices_per_cell)/))
  CALL C_F_POINTER(x_vertices_, x_vertices, (/nbr_vertices/))
  CALL C_F_POINTER(y_vertices_, y_vertices, (/nbr_vertices/))
  CALL C_F_POINTER(x_cells_, x_cells, (/nbr_cells/))
  CALL C_F_POINTER(y_cells_, y_cells, (/nbr_cells/))
  CALL C_F_POINTER(global_cell_ids_, global_cell_ids, (/nbr_cells/))
  CALL C_F_POINTER(cell_mask_, cell_mask, (/nbr_cells/))
  CALL C_F_POINTER(cell_core_mask_, cell_core_mask, (/nbr_cells/))
  CALL C_F_POINTER(global_corner_ids_, global_corner_ids, (/nbr_vertices/))
  CALL C_F_POINTER(corner_core_mask_, corner_core_mask, (/nbr_vertices/))

  cell_to_vertex = cell_to_vertex + 1

END SUBROUTINE fread_part_mpiom_grid_information

END PROGRAM main
