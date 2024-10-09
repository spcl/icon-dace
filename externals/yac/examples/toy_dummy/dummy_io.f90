! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

PROGRAM dummy_io

  USE mpi
  USE yac

  IMPLICIT NONE

  INTEGER, PARAMETER :: pd =  12
  INTEGER, PARAMETER :: rd = 307

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd) !< double precision
  INTEGER, PARAMETER :: wp = dp                        !< selected working precision

  INTEGER, PARAMETER :: no_of_fields = 8
  INTEGER, PARAMETER :: max_char_length = 132

  INTEGER, PARAMETER :: nbr_cells = 2
  INTEGER, PARAMETER :: nbr_vertices = 4

  REAL(wp), PARAMETER :: YAC_RAD = 0.017453292519943295769_wp ! M_PI / 180

  CHARACTER(LEN=max_char_length) :: field_name(no_of_fields)
  CHARACTER(LEN=max_char_length) :: yaml_filename
  CHARACTER(LEN=max_char_length) :: grid_name
  CHARACTER(LEN=max_char_length) :: comp_name

  INTEGER :: i, igrid, ierror

  INTEGER :: comp_id
  INTEGER :: cell_point_ids(2)
  INTEGER :: grid_ids(2)

  INTEGER :: glb_index(nbr_cells)
  INTEGER :: cell_core_mask(nbr_cells)
  INTEGER :: nbr_vertices_per_cell

  REAL(wp), ALLOCATABLE :: buffer_lon(:)
  REAL(wp), ALLOCATABLE :: buffer_lat(:)
  INTEGER, ALLOCATABLE  :: cell_to_vertex(:,:)

  INTEGER, ALLOCATABLE  :: cell_mask(:)
  INTEGER, ALLOCATABLE  :: field_id(:)

  INTEGER :: local_comm, npes, rank

  CALL mpi_init (ierror)

  ! Initialise the coupler
  CALL yac_finit ( )
  yaml_filename = "toy_dummy.yaml" ! default configuration file name
  CALL parse_arguments(yaml_filename)
  CALL yac_fread_config_yaml(yaml_filename)

  ! Inform the coupler about what we are
  comp_name = "dummy_io"
  CALL yac_fdef_comp ( comp_name, comp_id )

  CALL yac_fget_comp_comm ( comp_id, local_comm )

  CALL mpi_comm_rank ( local_comm, rank, ierror )
  CALL mpi_comm_size ( local_comm, npes, ierror )

  WRITE ( 6 , * ) TRIM(comp_name), " rank ", rank, ": local size is ", npes

  ! Here we define two grids (although numrically identical) to mimick
  ! two different Output server groups

  DO igrid = 1, 2

     IF ( igrid == 1 ) THEN
        grid_name = "ocean_grid"
     ELSE
        grid_name = "atmos_grid"
     ENDIF

     ALLOCATE(buffer_lon(nbr_vertices))
     ALLOCATE(buffer_lat(nbr_vertices))
     ALLOCATE(cell_to_vertex(3,nbr_cells))

     nbr_vertices_per_cell = 3

     ! Define vertices

     !      1
     !     / \
     !    / o \
     !   /     \
     !  2-------3   Eq.
     !   \     /
     !    \ o /
     !     \ /
     !      4

     buffer_lon(1) =  0.0 * YAC_RAD
     buffer_lon(2) = -1.0 * YAC_RAD
     buffer_lon(3) =  1.0 * YAC_RAD
     buffer_lon(4) =  0.0 * YAC_RAD
     buffer_lat(1) =  1.0 * YAC_RAD
     buffer_lat(2) =  0.0 * YAC_RAD
     buffer_lat(3) =  0.0 * YAC_RAD
     buffer_lat(4) = -1.0 * YAC_RAD

     ! Connectivity
     cell_to_vertex(1,1) = 1
     cell_to_vertex(2,1) = 2
     cell_to_vertex(3,1) = 3 ! cell 1
     cell_to_vertex(1,2) = 2
     cell_to_vertex(2,2) = 4
     cell_to_vertex(3,2) = 3 ! cell 2

     ! Definition of an unstructured grid
     CALL yac_fdef_grid (          &
          & grid_name,             &
          & nbr_vertices,          &
          & nbr_cells,             &
          & nbr_vertices_per_cell, &
          & buffer_lon,            &
          & buffer_lat,            &
          & cell_to_vertex,        &
          & grid_ids(igrid) )

     ! Decomposition information

     DO i = 1, nbr_cells
        glb_index(i) = i
        cell_core_mask(i) = 1
     ENDDO

     CALL yac_fset_global_index ( &
          & glb_index,            &
          & YAC_LOCATION_CELL,    &
          & grid_ids(igrid) )
     CALL yac_fset_core_mask ( &
          & cell_core_mask,    &
          & YAC_LOCATION_CELL, &
          & grid_ids(igrid) )

     ! Center points in cells (needed e.g. for nearest neighbour)

     buffer_lon(1) =  0.0 * YAC_RAD
     buffer_lon(2) =  0.0 * YAC_RAD
     buffer_lat(1) =  0.5 * YAC_RAD
     buffer_lat(2) = -0.5 * YAC_RAD

     CALL yac_fdef_points (           &
          & grid_ids(igrid),          &
          & nbr_cells,                &
          & YAC_LOCATION_CELL,        &
          & buffer_lon,               &
          & buffer_lat,               &
          & cell_point_ids(igrid) )

     DEALLOCATE (buffer_lon, buffer_lat, cell_to_vertex)

     ! Mask generation
     ALLOCATE(cell_mask(nbr_cells))
     DO i = 1, nbr_cells
        cell_mask(i) = 1
     ENDDO

     CALL yac_fset_mask (          &
          & cell_mask,             &
          & cell_point_ids(igrid) )

     DEALLOCATE (cell_mask)

  ENDDO

  field_name(1) = "atmos_out1" ! output field
  field_name(2) = "atmos_out2" ! output field
  field_name(3) = "atmos_out3" ! output field
  field_name(4) = "atmos_out4" ! output field

  field_name(5) = "ocean_out1" ! output field
  field_name(6) = "ocean_out2" ! output field
  field_name(7) = "ocean_out3" ! output field
  field_name(8) = "ocean_out4" ! output field

  ALLOCATE(field_id(no_of_fields))

  ! atmosphere output

  DO i = 1, no_of_fields-4
     CALL yac_fdef_field (        &
          & field_name(i),        &
          & comp_id,              &
          & cell_point_ids(2),    &
          & 1,                    &
          & 1,                    &
          & "1",                  &
          & YAC_TIME_UNIT_SECOND, &
          & field_id(i) )
  ENDDO

  ! ocean output

  DO i = no_of_fields-3, no_of_fields
     CALL yac_fdef_field (        &
          & field_name(i),        &
          & comp_id,              &
          & cell_point_ids(1),    &
          & 1,                    &
          & 1,                    &
          & "1",                  &
          & YAC_TIME_UNIT_SECOND, &
          & field_id(i) )
  ENDDO

  CALL yac_fenddef ( )

  CALL yac_ffinalize

  CALL mpi_finalize (ierror)

CONTAINS

  SUBROUTINE parse_arguments(configFilename)

    CHARACTER(LEN=max_char_length) :: configFilename

    CHARACTER(LEN=max_char_length) :: arg
    INTEGER :: i
    LOGICAL :: skip_arg = .false.

    DO i = 1, command_argument_count()

      IF (.NOT. skip_arg) THEN

        CALL get_command_argument(i, arg)

        SELECT CASE (arg)

          CASE ('-c')
            IF (i == command_argument_count()) THEN
              print '(2a, /)', 'missing parameter for command-line option: ', arg
              print '(a, /)', 'command-line options:'
              print '(a)',    '  -c configFilename'
              STOP
            ELSE
              CALL get_command_argument(i+1, configFilename)
              skip_arg = .true.
            END IF

          CASE default
            print '(2a, /)', 'unrecognised command-line option: ', arg
            print '(a, /)', 'command-line options:'
            print '(a)',    '  -c configFilename'
            STOP
        END SELECT
      ELSE
        skip_arg = .false.
      END IF
    END DO

  END SUBROUTINE parse_arguments

END PROGRAM dummy_io
