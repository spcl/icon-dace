! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#define DUMMY_VALUE (-1337.0d0)

#include "test_macros.inc"

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

MODULE field_data
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: field_out_data(:,:)
  DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: field_in_data(:,:)
  INTEGER, PUBLIC :: field_data_size
  INTEGER, PUBLIC :: point_id = HUGE(0)
  CHARACTER(LEN=64), PUBLIC :: config_filename
END MODULE

PROGRAM main

  USE mpi
  USE utest
  USE yac
  USE field_data
  USE, INTRINSIC :: iso_c_binding, only : c_null_char

  IMPLICIT NONE

  INTEGER :: global_rank, ierror

  INTEGER, PARAMETER :: max_opt_arg_len = 1024
  CHARACTER(max_opt_arg_len) :: grid_dir
  INTEGER :: arg_len

  INTERFACE

    SUBROUTINE C_UNLINK ( path ) BIND ( c, name='unlink' )

      USE, INTRINSIC :: iso_c_binding, only : c_char

      CHARACTER(KIND=c_char), DIMENSION(*) :: path

    END SUBROUTINE C_UNLINK

  END INTERFACE

  ! ===================================================================

  CALL start_test('dummy_restart_dble')

  CALL test(COMMAND_ARGUMENT_COUNT() == 1)
  CALL GET_COMMAND_ARGUMENT(1, grid_dir, arg_len)

  config_filename = 'test_restart_dble.yaml'

  ! run the models multiple time to check whether yac can be restarted

  ! no coupling
  CALL run_model_config(.FALSE., .FALSE., .FALSE., grid_dir)
  CALL run_model_config(.FALSE., .FALSE., .TRUE., grid_dir)
  ! coupling only in one direction
  CALL run_model_config(.TRUE., .FALSE., .FALSE., grid_dir)
  CALL run_model_config(.TRUE., .FALSE., .TRUE., grid_dir)
  ! no coupling
  CALL run_model_config(.FALSE., .TRUE., .FALSE., grid_dir)
  CALL run_model_config(.FALSE., .TRUE., .TRUE., grid_dir)
  ! coupling in both direction
  CALL run_model_config(.TRUE., .TRUE., .FALSE., grid_dir)
  CALL run_model_config(.TRUE., .TRUE., .TRUE., grid_dir)

  DEALLOCATE(field_out_data)
  DEALLOCATE(field_in_data)

  CALL MPI_Comm_rank(MPI_COMM_WORLD, global_rank, ierror)
  IF (global_rank == 0) CALL C_UNLINK(TRIM(config_filename) // c_null_char)

  CALL yac_ffinalize()

  CALL stop_test
  CALL exit_tests

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

  SUBROUTINE run_model_config( &
    icon_to_cube, cube_to_icon, config_from_file, grid_dir)

    USE mpi

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: icon_to_cube
    LOGICAL, INTENT(IN) :: cube_to_icon
    LOGICAL, INTENT(IN) :: config_from_file
    CHARACTER (LEN=*)   :: grid_dir

    INTEGER :: ierror, info
    INTEGER :: global_rank, global_size
    INTEGER :: local_comm, local_rank, local_size
    LOGICAL :: is_icon

    INTEGER :: comp_id
    INTEGER :: field_out_id, field_in_id
    LOGICAL :: receives_data
    INTEGER :: t, i
    INTEGER :: interp_stack_id
    INTEGER, PARAMETER :: collection_size = 1

    CHARACTER (LEN=:), ALLOCATABLE :: config

    INTERFACE
      SUBROUTINE generate_icon_grid(comm_rank, comm_size, grid_dir)
        INTEGER, INTENT(IN) :: comm_rank
        INTEGER, INTENT(IN) :: comm_size
        CHARACTER (LEN=*)   :: grid_dir
      END SUBROUTINE generate_icon_grid
      SUBROUTINE generate_cube_grid(comm_rank, comm_size)
        INTEGER, INTENT(IN) :: comm_rank
        INTEGER, INTENT(IN) :: comm_size
      END SUBROUTINE generate_cube_grid
    END INTERFACE

    CALL yac_finit()
    CALL MPI_Comm_rank(MPI_COMM_WORLD, global_rank, ierror)
    CALL MPI_Comm_size(MPI_COMM_WORLD, global_size, ierror)

    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fdef_datetime( &
      "+1800-01-01T00:00:00.000", "+2100-01-01T00:00:00.000")

    IF (global_size < 3) THEN
      WRITE ( * , * ) "Wrong number of processes (should be at least 3)"
      CALL error_exit
    ENDIF

    is_icon = global_rank < (global_size / 2)

    ! register component
    CALL yac_fdef_comp(MERGE('icon', 'cube', is_icon), comp_id)

    CALL yac_fget_comp_comm(comp_id, local_comm)
    CALL MPI_Comm_rank(local_comm, local_rank, ierror)
    CALL MPI_Comm_size(local_comm, local_size, ierror)
    CALL MPI_Comm_free(local_comm, ierror)

    ! if no points have been registered yet
    IF (point_id == HUGE(point_id)) THEN
      IF (is_icon) THEN
        CALL generate_icon_grid(local_rank, local_size, grid_dir)
      ELSE
        CALL generate_cube_grid(local_rank, local_size)
      END IF
    END IF

    ! register grid
    IF (is_icon) THEN
      receives_data = cube_to_icon
      CALL yac_fdef_field( &
        'icon_to_cube', comp_id, (/point_id/), 1, 1, "1", &
        YAC_TIME_UNIT_SECOND, field_out_id)
      CALL yac_fdef_field( &
        'cube_to_icon', comp_id, (/point_id/), 1, 1, "1", &
        YAC_TIME_UNIT_SECOND, field_in_id)
    ELSE
      receives_data = icon_to_cube
      CALL yac_fdef_field( &
        'cube_to_icon', comp_id, (/point_id/), 1, 1, "1", &
        YAC_TIME_UNIT_SECOND, field_out_id)
      CALL yac_fdef_field( &
        'icon_to_cube', comp_id, (/point_id/), 1, 1, "1", &
        YAC_TIME_UNIT_SECOND, field_in_id)
    END IF

    ! end definition phase
    CALL yac_fsync_def()

    IF (config_from_file) THEN

      CALL yac_fread_config_yaml(config_filename)
      CALL yac_fenddef()

    ELSE

      ! generate an interpolation interpolation stack
      CALL yac_fget_interp_stack_config(interp_stack_id)
      CALL yac_fadd_interp_stack_config_conservative( &
        interp_stack_id, 1, 0, 0, 0)

      ! generate couplings
      IF (icon_to_cube .AND. .NOT. is_icon) THEN
        CALL yac_fdef_couple( &
          'icon', 'icon', 'icon_to_cube', &
          'cube', 'cube', 'icon_to_cube', &
          '1', YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
          interp_stack_id)
      END IF
      IF (cube_to_icon .AND. is_icon) THEN
        CALL yac_fdef_couple( &
          'cube', 'cube', 'cube_to_icon', &
          'icon', 'icon', 'cube_to_icon', &
          '1', YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
          interp_stack_id)
      END IF

      CALL yac_ffree_interp_stack_config(interp_stack_id)

      CALL yac_fenddef(YAC_YAML_EMITTER_DEFAULT_F, config)

      IF (global_rank == 0) THEN
        open(unit=999, file=config_filename, status='replace', &
             form='formatted', action='write')
        write(999,"(A)") config
        close(999)
      END IF

      DEALLOCATE(config)
    END IF

    DO t = 1, 10

      CALL yac_fput( &
        field_out_id, SIZE(field_out_data, 1), collection_size, &
        field_out_data, info, ierror)

      field_in_data = DUMMY_VALUE
      CALL yac_fget( &
        field_in_id, SIZE(field_in_data, 1), collection_size, &
        field_in_data, info, ierror)

      DO i = 1, field_data_size

        IF (receives_data) THEN
          CALL test(ABS(field_in_data(i,1) - field_out_data(i,1)) < 0.001)
        ELSE
          CALL test(field_in_data(i,1) == DUMMY_VALUE)
        END IF
      END DO
    END DO

    CALL yac_fcleanup()

  END SUBROUTINE run_model_config

END PROGRAM

SUBROUTINE generate_icon_grid(comm_rank, comm_size, grid_dir)

  USE yac
  USE field_data, ONLY: field_out_data, field_in_data, field_data_size, point_id
  USE, INTRINSIC :: iso_c_binding, ONLY : c_loc

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE C_FREE ( ptr ) BIND ( c, NAME='free' )

      USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr

      TYPE ( c_ptr ), INTENT(IN), VALUE :: ptr

    END SUBROUTINE C_FREE

  END INTERFACE

  INTEGER, INTENT(IN) :: comm_rank
  INTEGER, INTENT(IN) :: comm_size
  CHARACTER (LEN=*)   :: grid_dir

  INTEGER :: grid_id

  INTEGER :: nbr_vertices
  INTEGER :: nbr_cells
  INTEGER, POINTER :: num_vertices_per_cell(:)
  INTEGER, POINTER :: cell_to_vertex(:)
  DOUBLE PRECISION, POINTER :: x_vertices(:)
  DOUBLE PRECISION, POINTER :: y_vertices(:)
  DOUBLE PRECISION, POINTER :: x_cells(:)
  DOUBLE PRECISION, POINTER :: y_cells(:)

  INTEGER, POINTER :: cell_mask_int(:)
  INTEGER, POINTER :: global_cell_id(:)
  INTEGER, POINTER :: cell_core_mask(:)
  INTEGER, POINTER :: global_corner_id(:)
  INTEGER, POINTER :: corner_core_mask(:)

  LOGICAL, ALLOCATABLE :: cell_mask(:)

  CALL yac_read_part_icon_grid_information(                                  &
    TRIM(grid_dir) // 'icon_grid_0030_R02B03_G.nc', nbr_vertices, nbr_cells, &
    num_vertices_per_cell, cell_to_vertex, x_vertices, y_vertices,           &
    x_cells, y_cells, global_cell_id, cell_mask_int, cell_core_mask,         &
    global_corner_id, corner_core_mask, comm_rank, comm_size)

  ALLOCATE(cell_mask(nbr_cells))
  cell_mask = cell_mask_int == 0

  CALL yac_fdef_grid(                                            &
    'icon', nbr_vertices, nbr_cells, SUM(num_vertices_per_cell), &
    num_vertices_per_cell, x_vertices, y_vertices, cell_to_vertex, grid_id)

  CALL yac_fset_global_index(global_corner_id, YAC_LOCATION_CORNER, grid_id)
  CALL yac_fset_core_mask(corner_core_mask == 1, YAC_LOCATION_CORNER, grid_id)
  CALL yac_fset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id)
  CALL yac_fset_core_mask(cell_core_mask == 1, YAC_LOCATION_CELL, grid_id)

  CALL yac_fdef_points( &
    grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, point_id)

  CALL yac_fset_mask(cell_mask, point_id)

  ALLOCATE(field_out_data(nbr_cells,1))
  field_data_size = nbr_cells
  field_out_data(:,1) =                          &
    MERGE(test_harmonic(x_cells(:), y_cells(:)), &
          DUMMY_VALUE, cell_core_mask == 1)
  ALLOCATE(field_in_data(nbr_cells,1))

  IF (ASSOCIATED(num_vertices_per_cell)) CALL C_FREE(c_loc(num_vertices_per_cell(1)))
  IF (ASSOCIATED(cell_to_vertex)) CALL C_FREE(c_loc(cell_to_vertex(1)))
  IF (ASSOCIATED(x_vertices)) CALL C_FREE(c_loc(x_vertices(1)))
  IF (ASSOCIATED(y_vertices)) CALL C_FREE(c_loc(y_vertices(1)))
  IF (ASSOCIATED(x_cells)) CALL C_FREE(c_loc(x_cells(1)))
  IF (ASSOCIATED(y_cells)) CALL C_FREE(c_loc(y_cells(1)))

  IF (ASSOCIATED(cell_mask_int)) CALL C_FREE(c_loc(cell_mask_int(1)))
  IF (ASSOCIATED(global_cell_id)) CALL C_FREE(c_loc(global_cell_id(1)))
  IF (ASSOCIATED(cell_core_mask)) CALL C_FREE(c_loc(cell_core_mask(1)))
  IF (ASSOCIATED(global_corner_id)) CALL C_FREE(c_loc(global_corner_id(1)))
  IF (ASSOCIATED(corner_core_mask)) CALL C_FREE(c_loc(corner_core_mask(1)))

CONTAINS

  ELEMENTAL DOUBLE PRECISION FUNCTION test_harmonic ( lon, lat )
    DOUBLE PRECISION, INTENT(IN) :: lon, lat

    DOUBLE PRECISION, PARAMETER :: YAC_RAD = 0.017453292519943 ! M_PI / 180.0
    test_harmonic = &
      2.0 + ((sin(2.0 * lat * YAC_RAD))**16.0  * cos(16.0 * lon * YAC_RAD))
  END FUNCTION test_harmonic

  SUBROUTINE yac_read_part_icon_grid_information ( grid_file,             &
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

      SUBROUTINE yac_read_part_icon_grid_information_c ( grid_file,              &
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
        BIND ( c, name='yac_read_part_icon_grid_information' )

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

      END SUBROUTINE yac_read_part_icon_grid_information_c

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

    CALL yac_read_part_icon_grid_information_c ( TRIM(grid_file) // c_null_char, &
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

  END SUBROUTINE yac_read_part_icon_grid_information

END SUBROUTINE generate_icon_grid

! ----------------------------------------------------------

SUBROUTINE generate_cube_grid(comm_rank, comm_size)

  USE yac
  USE field_data, ONLY: field_out_data, field_in_data, field_data_size, point_id
  USE, INTRINSIC :: iso_c_binding, ONLY : c_loc

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE C_FREE ( ptr ) BIND ( c, NAME='free' )

      USE, INTRINSIC :: iso_c_binding, ONLY : c_ptr

      TYPE ( c_ptr ), INTENT(IN), VALUE :: ptr

    END SUBROUTINE C_FREE

  END INTERFACE

  INTEGER, INTENT(IN) :: comm_rank
  INTEGER, INTENT(IN) :: comm_size

  INTEGER :: grid_id

  INTEGER, PARAMETER :: n = 50

  INTEGER :: nbr_vertices
  INTEGER :: nbr_cells
  INTEGER, POINTER :: num_vertices_per_cell(:)
  INTEGER, POINTER :: cell_to_vertex(:)
  DOUBLE PRECISION, POINTER :: x_vertices(:)
  DOUBLE PRECISION, POINTER :: y_vertices(:)
  DOUBLE PRECISION, POINTER :: x_cells(:)
  DOUBLE PRECISION, POINTER :: y_cells(:)

  INTEGER, POINTER :: cell_core_mask(:)
  INTEGER, POINTER :: corner_core_mask(:)
  INTEGER, POINTER :: global_cell_id(:)
  INTEGER, POINTER :: global_corner_id(:)

  CALL yac_generate_part_cube_grid_information(                        &
    n, nbr_vertices, nbr_cells, num_vertices_per_cell, cell_to_vertex, &
    x_vertices, y_vertices, x_cells, y_cells, global_cell_id,          &
    cell_core_mask, global_corner_id, corner_core_mask,                &
    comm_rank, comm_size)

  CALL yac_fdef_grid(                                            &
    'cube', nbr_vertices, nbr_cells, SUM(num_vertices_per_cell), &
    num_vertices_per_cell, x_vertices, y_vertices, cell_to_vertex, grid_id)

  CALL yac_fset_global_index(global_corner_id, YAC_LOCATION_CORNER, grid_id)
  CALL yac_fset_core_mask(corner_core_mask == 1, YAC_LOCATION_CORNER, grid_id)
  CALL yac_fset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id)
  CALL yac_fset_core_mask(cell_core_mask == 1, YAC_LOCATION_CELL, grid_id)

  CALL yac_fdef_points( &
    grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, point_id)

  ALLOCATE(field_out_data(nbr_cells,1))
  field_data_size = nbr_cells
  field_out_data(:,1) =                          &
    MERGE(test_harmonic(x_cells(:), y_cells(:)), &
          DUMMY_VALUE, cell_core_mask == 1)
  ALLOCATE(field_in_data(nbr_cells,1))

  IF (ASSOCIATED(num_vertices_per_cell)) CALL C_FREE(c_loc(num_vertices_per_cell(1)))
  IF (ASSOCIATED(cell_to_vertex)) CALL C_FREE(c_loc(cell_to_vertex(1)))
  IF (ASSOCIATED(x_vertices)) CALL C_FREE(c_loc(x_vertices(1)))
  IF (ASSOCIATED(y_vertices)) CALL C_FREE(c_loc(y_vertices(1)))
  IF (ASSOCIATED(x_cells)) CALL C_FREE(c_loc(x_cells(1)))
  IF (ASSOCIATED(y_cells)) CALL C_FREE(c_loc(y_cells(1)))

  IF (ASSOCIATED(global_cell_id)) CALL C_FREE(c_loc(global_cell_id(1)))
  IF (ASSOCIATED(cell_core_mask)) CALL C_FREE(c_loc(cell_core_mask(1)))
  IF (ASSOCIATED(global_corner_id)) CALL C_FREE(c_loc(global_corner_id(1)))
  IF (ASSOCIATED(corner_core_mask)) CALL C_FREE(c_loc(corner_core_mask(1)))

CONTAINS

  ELEMENTAL DOUBLE PRECISION FUNCTION test_harmonic ( lon, lat )
    DOUBLE PRECISION, INTENT(IN) :: lon, lat

    DOUBLE PRECISION, PARAMETER :: YAC_RAD = 0.017453292519943 ! M_PI / 180.0
    test_harmonic = &
      2.0 + ((sin(2.0 * lat * YAC_RAD))**16.0  * cos(16.0 * lon * YAC_RAD))
  END FUNCTION test_harmonic

  SUBROUTINE yac_generate_part_cube_grid_information ( n,                     &
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

      SUBROUTINE yac_generate_part_cube_grid_information_c ( n,                      &
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
        BIND ( c, name='yac_generate_part_cube_grid_information' )

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

      END SUBROUTINE yac_generate_part_cube_grid_information_c

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

    CALL yac_generate_part_cube_grid_information_c ( n,                        &
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

  END SUBROUTINE yac_generate_part_cube_grid_information

END SUBROUTINE generate_cube_grid
