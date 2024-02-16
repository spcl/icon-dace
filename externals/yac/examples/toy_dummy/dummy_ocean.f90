PROGRAM dummy_ocean

  USE mpi
  USE mo_yac_finterface

  IMPLICIT NONE

  INTEGER, PARAMETER :: pd =  12
  INTEGER, PARAMETER :: rd = 307

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd) !< double precision
  INTEGER, PARAMETER :: wp = dp                        !< selected working precision

  INTEGER, PARAMETER :: no_of_fields = 14
  INTEGER, PARAMETER :: max_char_length = 132

  INTEGER, PARAMETER :: nbr_cells = 2
  INTEGER, PARAMETER :: nbr_vertices = 4

  REAL(wp), PARAMETER :: YAC_RAD = 0.017453292519943295769_wp ! M_PI / 180

  CHARACTER(LEN=max_char_length) :: dummy_name
  CHARACTER(LEN=max_char_length) :: field_name(no_of_fields)
  INTEGER                        :: field_collection_size(no_of_fields)
  CHARACTER(LEN=max_char_length) :: yaml_filename
  CHARACTER(LEN=max_char_length) :: grid_name
  CHARACTER(LEN=max_char_length) :: comp_name
  CHARACTER(LEN=max_char_length) :: timestep_string

  INTEGER :: role
  INTEGER :: i, info, ierror
  INTEGER :: comp_id
  INTEGER :: cell_point_ids(1)
  INTEGER :: grid_id

  INTEGER :: glb_index(nbr_cells)
  INTEGER :: cell_core_mask(nbr_cells)
  INTEGER :: nbr_vertices_per_cell

  REAL(wp), ALLOCATABLE :: buffer(:,:)
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
  comp_name = "dummy_ocean"
  grid_name = "dummy_ocean_grid"
  CALL yac_fdef_comp ( comp_name, comp_id )

  print *, "YAC Version: ", TRIM(yac_fget_version())

  CALL yac_fget_comp_comm ( comp_id, local_comm )
  print *, 'Local Comm', local_comm

  CALL mpi_comm_rank ( local_comm, rank, ierror )
  CALL mpi_comm_size ( local_comm, npes, ierror )

  WRITE ( 6 , * ) TRIM(comp_name), " rank ", rank, ": local size is ", npes

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
       & grid_id )

  ! Decomposition information

  DO i = 1, nbr_cells
     glb_index(i) = i
     cell_core_mask(i) = 1
  ENDDO

  CALL yac_fset_global_index ( &
       & glb_index,            &
       & YAC_LOCATION_CELL,    &
       & grid_id )
  CALL yac_fset_core_mask ( &
       & cell_core_mask,    &
       & YAC_LOCATION_CELL, &
       & grid_id )

  ! Center points in cells (needed e.g. for nearest neighbour)

  buffer_lon(1) =  0.0 * YAC_RAD
  buffer_lon(2) =  0.0 * YAC_RAD
  buffer_lat(1) =  0.5 * YAC_RAD
  buffer_lat(2) = -0.5 * YAC_RAD

  CALL yac_fdef_points (    &
       & grid_id,           &
       & nbr_cells,         &
       & YAC_LOCATION_CELL, &
       & buffer_lon,        &
       & buffer_lat,        &
       & cell_point_ids(1) )

  DEALLOCATE (buffer_lon, buffer_lat, cell_to_vertex)

  ! Mask generation
  ALLOCATE(cell_mask(nbr_cells))
  DO i = 1, nbr_cells
     cell_mask(i) = 1
  ENDDO

  CALL yac_fset_mask ( &
       & cell_mask,    &
       & cell_point_ids(1) )

  DEALLOCATE (cell_mask)

  field_name(1) =  "surface_downward_eastward_stress"  ! bundled field containing two components
  field_name(2) =  "surface_downward_northward_stress" ! bundled field containing two components
  field_name(3) =  "surface_fresh_water_flux"          ! bundled field containing three components
  field_name(4) =  "surface_temperature"
  field_name(5) =  "total_heat_flux"                   ! bundled field containing four components
  field_name(6) =  "atmosphere_sea_ice_bundle"         ! bundled field containing four components
  field_name(7) =  "sea_surface_temperature"
  field_name(8) =  "eastward_sea_water_velocity"
  field_name(9) =  "northward_sea_water_velocity"
  field_name(10) = "ocean_sea_ice_bundle"              ! bundled field containing four components
  field_name(11) = "ocean_out1"                        ! output field
  field_name(12) = "ocean_out2"                        ! output field
  field_name(13) = "ocean_out3"                        ! output field
  field_name(14) = "ocean_out4"                        ! output field

  field_collection_size(1) =  2
  field_collection_size(2) =  2
  field_collection_size(3) =  3
  field_collection_size(4) =  1
  field_collection_size(5) =  4
  field_collection_size(6) =  4
  field_collection_size(7) =  1
  field_collection_size(8) =  1
  field_collection_size(9) =  1
  field_collection_size(10) = 5
  field_collection_size(11) = 1
  field_collection_size(12) = 1
  field_collection_size(13) = 1
  field_collection_size(14) = 1

  ALLOCATE(field_id(no_of_fields))

  ! fields for coupling

  DO i = 1, no_of_fields-4
     CALL yac_fdef_field (            &
          & field_name(i),            &
          & comp_id,                  &
          & cell_point_ids,           &
          & 1,                        &
          & field_collection_size(i), &
          & "1",                      &
          & YAC_TIME_UNIT_SECOND,     &
          & field_id(i) )
  ENDDO

  ! fields for output server

  DO i = no_of_fields-3, no_of_fields
     CALL yac_fdef_field (            &
          & field_name(i),            &
          & comp_id,                  &
          & cell_point_ids,           &
          & 1,                        &
          & field_collection_size(i), &
          & "1",                      &
          & YAC_TIME_UNIT_SECOND,     &
          & field_id(i) )
  ENDDO

  CALL yac_fenddef ( )

  ! Data exchange

  ALLOCATE(buffer(nbr_cells,5))
  buffer(:,:) = 0.0_wp

  !   field_id(1) represents "TAUX"   wind stress component
  !   field_id(2) represents "TAUY"   wind stress component
  !   field_id(3) represents "SFWFLX" surface fresh water flux
  !   field_id(4) represents "SFTEMP" surface temperature
  !   field_id(5) represents "THFLX"  total heat flux
  !   field_id(6) represents "ICEATM" ice temperatures and melt potential
  !
  !   field_id(7) represents "SST"    sea surface temperature
  !   field_id(8) represents "OCEANU" u component of ocean surface current
  !   field_id(9) represents "OCEANV" v component of ocean surface current
  !   field_id(10)represents "ICEOCE" ice thickness, concentration and temperatures
  !
  !   field_id(11) - field_id(14) represent output fields
  !
  !   Get some info back from the coupling configuration
  !
  DO i = 1, no_of_fields
    timestep_string = yac_fget_field_timestep ( field_id(i) )
    WRITE ( 6 , * ) "Field ID ", field_id(i), TRIM(timestep_string)
  ENDDO

  DO i = 1, no_of_fields
    role = yac_fget_field_role ( field_id(i) )
    dummy_name = yac_fget_field_name_from_field_id ( field_id(i) )
    WRITE ( 6 , * ) "Requested role for ", TRIM(dummy_name) , " is ", role
  ENDDO

  !
  ! Send fields from ocean to atmosphere
  ! ------------------------------------
  !
  ! SST
  buffer(:,1) = 110.0_wp
  CALL yac_fput ( field_id(7), nbr_cells, 1, buffer(1:nbr_cells,1:1), info, ierror )
  IF ( info > 0 ) &
       WRITE ( 6 , * ) "ocean CPL SST", minval(buffer(1:nbr_cells,1:1)), maxval(buffer(1:nbr_cells,1:1))

  ! zonal velocity
  buffer(:,1) = 120.0_wp
  CALL yac_fput ( field_id(8), nbr_cells, 1, buffer(1:nbr_cells,1:1), info, ierror )
  IF ( info > 0 ) &
       WRITE ( 6 , * ) "ocean CPL U", minval(buffer(1:nbr_cells,1:1)), maxval(buffer(1:nbr_cells,1:1))

  ! meridional velocity
  buffer(:,1) = 130.0_wp
  CALL yac_fput ( field_id(9), nbr_cells, 1, buffer(1:nbr_cells,1:1), info, ierror )
  IF ( info > 0 ) &
       WRITE ( 6 , * ) "ocean CPL V", minval(buffer(1:nbr_cells,1:1)), maxval(buffer(1:nbr_cells,1:1))

  ! Ice thickness, concentration, T1 and T2
  buffer(:,1) = 140.1_wp
  buffer(:,2) = 140.2_wp
  buffer(:,3) = 140.3_wp
  buffer(:,4) = 140.4_wp
  buffer(:,5) = 140.5_wp

  CALL yac_fput ( field_id(10), nbr_cells, 5, buffer(1:nbr_cells,1:5), info, ierror )
  IF ( info > 0 ) THEN
     WRITE ( 6 , * ) "ocean CPL ice 1", minval(buffer(1:nbr_cells,1:1)), maxval(buffer(1:nbr_cells,1:1))
     WRITE ( 6 , * ) "ocean CPL ice 2", minval(buffer(1:nbr_cells,2:2)), maxval(buffer(1:nbr_cells,2:2))
     WRITE ( 6 , * ) "ocean CPL ice 3", minval(buffer(1:nbr_cells,3:3)), maxval(buffer(1:nbr_cells,3:3))
     WRITE ( 6 , * ) "ocean CPL ice 4", minval(buffer(1:nbr_cells,4:4)), maxval(buffer(1:nbr_cells,4:4))
     WRITE ( 6 , * ) "ocean CPL ice 5", minval(buffer(1:nbr_cells,5:5)), maxval(buffer(1:nbr_cells,5:5))
  ENDIF

  !
  ! Receive fields from atmosphere
  ! ------------------------------

  ! zonal wind stress

  CALL yac_fget ( field_id(1), nbr_cells, 2, buffer, info, ierror )

  ! meridional wind stress

  CALL yac_fget ( field_id(2), nbr_cells, 2, buffer, info, ierror )

  ! freshwater flux

  CALL yac_fget ( field_id(3), nbr_cells, 3, buffer, info, ierror )

  ! surface air temperature
  CALL yac_fget ( field_id(4), nbr_cells, 1, buffer, info, ierror )

  ! total heat flux - 4 parts - record 5

  CALL yac_fget ( field_id(5), nbr_cells, 4, buffer, info, ierror )

  ! ice parameter

  CALL yac_fget ( field_id(6), nbr_cells, 4, buffer, info, ierror )

  DEALLOCATE(buffer)
  DEALLOCATE(field_id)

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

END PROGRAM dummy_ocean
