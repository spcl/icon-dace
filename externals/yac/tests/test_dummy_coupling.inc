!!
!! @file test_dummy_coupling.inc
!!
!! @copyright Copyright  (C)  2013 DKRZ, MPI-M
!!
!! @author Rene Redler <rene.redler@mpimet.mpg.de>
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

#ifndef TEST_PRECISION
#error "TEST_PRECISION is not defined"
#endif

! module name has to be specific to each TEST_PRECISION in order
! to avoid conflicts, because the compiler may generate object files
! for each module using the module name as file name
MODULE MOD_NAME

  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: wp = TEST_PRECISION              !< selected working precision

  INTEGER, PARAMETER :: no_of_fields = 10
  INTEGER, PARAMETER :: nbr_cells = 2
  INTEGER, PARAMETER :: nbr_vertices = 4

  INTEGER, PARAMETER :: max_char_length = 132

  CHARACTER(LEN=max_char_length) :: dummy_name
  CHARACTER(LEN=max_char_length) :: field_name(no_of_fields)
  CHARACTER(LEN=max_char_length) :: yaml_filename
  CHARACTER(LEN=max_char_length) :: grid_name
  CHARACTER(LEN=max_char_length) :: comp_name
  CHARACTER(LEN=max_char_length) :: timestep_string
  INTEGER                        :: collection_size
  INTEGER                        :: ref_collection_size(no_of_fields)

  INTEGER :: info
  INTEGER :: comp_id
  INTEGER :: comp_ids(1)
  INTEGER :: grid_id
  INTEGER :: cell_point_ids(1)
  INTEGER :: cell_mask_ids(1)

  INTEGER :: glb_index(nbr_cells)
  INTEGER :: core_mask(nbr_cells)
  INTEGER :: nbr_vertices_per_cell

  REAL(wp), ALLOCATABLE :: buffer(:,:)
  REAL(wp), ALLOCATABLE :: buffer_lon(:)
  REAL(wp), ALLOCATABLE :: buffer_lat(:)
  INTEGER, ALLOCATABLE  :: cell_to_vertex(:,:)

  INTEGER, ALLOCATABLE  :: cell_mask(:)
  INTEGER, ALLOCATABLE  :: field_id(:)

END MODULE MOD_NAME

! ------------------------------------------------------------------------

#include "test_macros.inc"

! ------------------------------------------------------------------------

PROGRAM dummy_coupling

  USE mpi
  USE MOD_NAME
  USE mo_yac_finterface
  USE utest

  IMPLICIT NONE

  INTEGER :: ierror, npes, rank

  CALL start_test("dummy_coupling")

  ! Initialise the coupler
  yaml_filename = "coupling_test.yaml"
  CALL yac_finit ( )

  CALL mpi_comm_rank ( MPI_COMM_WORLD, rank, ierror )
  CALL mpi_comm_size ( MPI_COMM_WORLD, npes, ierror )

  SELECT CASE ( rank )

      CASE ( 0 )
          CALL dummy_atmosphere ( )

      CASE ( 1 )
          CALL dummy_ocean ( )

      CASE ( 2 )
          CALL dummy_io ( )

      CASE DEFAULT
          WRITE ( * , * ) "Too many processes have been launched"
          CALL test ( .FALSE. )
          CALL stop_test
          CALL exit_tests
          CALL mpi_abort (MPI_COMM_WORLD, 999, ierror)

  END SELECT

  CALL yac_ffinalize

  CALL stop_test
  CALL exit_tests

CONTAINS

! ------------------------------------------------------------------------

SUBROUTINE dummy_atmosphere ( )

  USE mpi
  USE MOD_NAME
  USE mo_yac_finterface
  USE utest

  IMPLICIT NONE

  INTEGER :: role
  INTEGER :: i, ierror
  INTEGER :: local_comm, npes, rank

  INTEGER :: iflag
  LOGICAL :: lflag

  ! Inform the coupler about what we are
  comp_name = "dummy_atmosphere"
  CALL yac_fdef_comp ( comp_name, comp_id )
  comp_ids(1) = comp_id

  CALL yac_fget_comp_comm ( comp_id, local_comm )

  CALL mpi_comm_rank ( local_comm, rank, ierror )
  CALL mpi_comm_size ( local_comm, npes, ierror )

  ! WRITE ( 6 , * ) TRIM(comp_name), " rank ", rank, ": local size is ", npes

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

  grid_name = "dummy_atmosphere_grid"

  buffer_lon(1) =  0.0; buffer_lat(1) =  1.0
  buffer_lon(2) = -1.0; buffer_lat(2) =  0.0
  buffer_lon(3) =  1.0; buffer_lat(3) =  0.0
  buffer_lon(4) =  0.0; buffer_lat(4) = -1.0

  ! Connectivity
  cell_to_vertex(1,1) = 1; cell_to_vertex(2,1) = 2; cell_to_vertex(3,1) = 3 ! cell 1
  cell_to_vertex(1,2) = 2; cell_to_vertex(2,2) = 4; cell_to_vertex(3,2) = 3 ! cell 2

  ! Define unstructured grid
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
     core_mask(i) = -1
  ENDDO

  CALL yac_fset_global_index ( &
       & glb_index,            &
       & YAC_LOCATION_CELL,    &
       & grid_id )

  CALL yac_fset_core_mask ( &
       & glb_index,          &
       & YAC_LOCATION_CELL,  &
       & grid_id )

  ! Center points in cells (needed e.g. for nearest neighbour)

  buffer_lon(1) = 0.0; buffer_lat(1) =  0.5
  buffer_lon(2) = 0.0; buffer_lat(2) = -0.5

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

  CALL yac_fdef_mask(       &
       & grid_id,           &
       & nbr_cells,         &
       & YAC_LOCATION_CELL, &
       & cell_mask,         &
       cell_mask_ids(1) );

  CALL yac_fset_mask ( &
       & cell_mask,    &
       & cell_point_ids(1) )

  DEALLOCATE (cell_mask)

  CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)

  field_name( 1) = "TAUX"   ! bundled field containing two components
  field_name( 2) = "TAUY"   ! bundled field containing two components
  field_name( 3) = "SFWFLX" ! bundled field containing three components
  field_name( 4) = "SFTEMP"
  field_name( 5) = "THFLX"  ! bundled field containing four components
  field_name( 6) = "ICEATM" ! bundled field containing four components
  field_name( 7) = "SST"
  field_name( 8) = "OCEANU"
  field_name( 9) = "OCEANV"
  field_name(10) = "ICEOCE" ! bundled field containing four components

  ref_collection_size( 1) = 2
  ref_collection_size( 2) = 2
  ref_collection_size( 3) = 3
  ref_collection_size( 4) = 1
  ref_collection_size( 5) = 4
  ref_collection_size( 6) = 4
  ref_collection_size( 7) = 1
  ref_collection_size( 8) = 1
  ref_collection_size( 9) = 1
  ref_collection_size(10) = 5

  ALLOCATE(field_id(no_of_fields))

  DO i = 1, no_of_fields/2
     CALL yac_fdef_field (          &
          & field_name(i),          &
          & comp_id,                &
          & cell_point_ids,         &
          & 1,                      &
          & ref_collection_size(i), &
          & "600",                  &
          & YAC_TIME_UNIT_SECOND,   &
          & field_id(i) )
  ENDDO

  DO i = no_of_fields/2+1, no_of_fields
     CALL yac_fdef_field_mask (     &
          & field_name(i),          &
          & comp_id,                &
          & cell_point_ids,         &
          & cell_mask_ids,          &
          & 1,                      &
          & ref_collection_size(i), &
          & "600",                  &
          & YAC_TIME_UNIT_SECOND,   &
          & field_id(i) )
  ENDDO

  ! read couplings from YAML configuration file
  CALL yac_fread_config_yaml(yaml_filename)

  CALL yac_fenddef ( )

  !these queries can only be fulfilled once the search has been completed.
  DO i = 1, no_of_fields
    timestep_string = yac_fget_timestep_from_field_id ( field_id(i) )
    IF ( TRIM(timestep_string) .NE. "PT10M" ) THEN
       WRITE ( 6 , * ) "Wrong time step ", TRIM(field_name(i)), " ", TRIM(timestep_string)
       CALL test (.FALSE. )
    ENDIF
  ENDDO

  DO i = 1, no_of_fields
    collection_size = yac_fget_collection_size_from_field_id ( field_id(i) )
    IF ( collection_size .NE. ref_collection_size(i) ) THEN
       WRITE ( 6 , * ) "Wrong collection size ", TRIM(field_name(i)), " ", collection_size
       CALL test (.FALSE. )
    ENDIF
  ENDDO

  DO i = 1, no_of_fields
    role = yac_fget_role_from_field_id ( field_id(i) )
    dummy_name = yac_fget_field_name_from_field_id ( field_id(i) )
    IF ( i < 7 .AND. role .NE. 1 ) THEN
      WRITE ( 6 , * ) "Wrong requested role for ", TRIM(dummy_name) , " is ", role
      CALL test (.FALSE. )
    ELSE IF ( i > 6 .AND. role .NE. 2 ) THEN
      WRITE ( 6 , * ) "Wrong requested role for ", TRIM(dummy_name) , " is ", role
      CALL test (.FALSE. )
    ENDIF
  ENDDO

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
  ! Send fields to oceanid
  ! --------------------

  ! meridional wind stress
  buffer(:,1) = 10.1_wp
  buffer(:,2) = 10.2_wp

  CALL yac_fput ( field_id(1), nbr_cells, 2, buffer(1:nbr_cells,1:2), info, ierror )
  ! IF ( info > 0 ) &
  !   WRITE ( 6 , * ) "atmosphere CPL TAUX 1", &
  !                   MINVAL(buffer(1:nbr_cells,1:1)), &
  !                   MAXVAL(buffer(1:nbr_cells,1:1))
  ! IF ( info > 0 ) &
  !   WRITE ( 6 , * ) "atmosphere CPL TAUX 2", &
  !                   MINVAL(buffer(1:nbr_cells,2:2)), &
  !                   MAXVAL(buffer(1:nbr_cells,2:2))


  ! zonal  wind stress
  buffer(:,1) = 20.1_wp
  buffer(:,2) = 20.2_wp

  CALL yac_fput ( field_id(2), nbr_cells, 2, buffer(1:nbr_cells,1:2), info, ierror )
  ! IF ( info > 0 ) &
  !   WRITE ( 6 , * ) "atmosphere CPL TAUY 1", &
  !                   MINVAL(buffer(1:nbr_cells,1:1)), &
  !                   MAXVAL(buffer(1:nbr_cells,1:1))
  ! IF ( info > 0 ) &
  !   WRITE ( 6 , * ) "atmosphere CPL TAUY 2", &
  !                   MINVAL(buffer(1:nbr_cells,2:2)), &
  !                   MAXVAL(buffer(1:nbr_cells,2:2))

  ! surface fresh water flux
  buffer(:,1) = 30.1_wp
  buffer(:,2) = 30.2_wp
  buffer(:,3) = 30.3_wp

  CALL yac_fput ( field_id(3), nbr_cells, 3, buffer(1:nbr_cells,1:3), info, ierror )

  ! surface temperature
  buffer(:,1) = 40.1_wp
  CALL yac_fput ( field_id(4), nbr_cells, 1, buffer(1:nbr_cells,1:1), info, ierror )

  ! total heat flux
  buffer(:,1) = 50.1_wp
  buffer(:,2) = 50.2_wp
  buffer(:,3) = 50.3_wp
  buffer(:,4) = 50.4_wp

  CALL yac_fput ( field_id(5), nbr_cells, 4, buffer(1:nbr_cells,1:4), info, ierror )

  ! ice temperatures and melt potential
  buffer(:,1) = 60.1_wp
  buffer(:,2) = 60.2_wp
  buffer(:,3) = 60.3_wp
  buffer(:,4) = 60.4_wp

  CALL yac_fput ( field_id(6), nbr_cells, 4, buffer(1:nbr_cells,1:4), info, ierror )

   ! check whether yac_ftest does not crash here
   CALL yac_ftest ( field_id(6), lflag )

  !
  ! Receive fields from ocean
  ! -------------------------

  ! SST
  CALL yac_fget ( field_id(7), nbr_cells, 1, buffer(1:nbr_cells,1:1), info, ierror )
  IF ( info > 0 .AND. buffer(1,1) /= 110.0_wp ) THEN
     WRITE ( 6 , * ) "atmosphere CPL SST", buffer(1,1)
     CALL test (.FALSE. )
  ENDIF

  ! zonal velocity
  CALL yac_fget ( field_id(8), nbr_cells, 1, buffer(1:nbr_cells,1:1), info, ierror )
  IF ( info > 0 .AND. buffer(1,1) /= 120.0_wp ) THEN
     WRITE ( 6 , * ) "atmosphere CPL OCEANU ", buffer(1,1)
     CALL test (.FALSE. )
  ENDIF

  ! meridional velocity
  CALL yac_fget ( field_id(9), nbr_cells, 1, buffer(1:nbr_cells,1:1), info, ierror )
  IF ( info > 0 .AND. buffer(1,1) /= 130.0_wp ) THEN
     WRITE ( 6 , * ) "atmosphere CPL OCEANV ", buffer(1,1)
     CALL test (.FALSE. )
  ENDIF

  ! Ice thickness, concentration, T1 and T2

  CALL yac_fget ( field_id(10), nbr_cells, 5, buffer(1:nbr_cells,1:5), info, ierror )

  IF ( info > 0 ) THEN
     IF ( buffer(1,1) /= 140.1_wp ) THEN
        WRITE ( 6 , * ) "atmosphere CPL ice 1 ", buffer(1,1)
        CALL test (.FALSE. )
     ENDIF

     IF ( buffer(1,2) /= 140.2_wp ) THEN
        WRITE ( 6 , * ) "atmosphere CPL ice 2 ", buffer(1,1)
        CALL test (.FALSE. )
     ENDIF
     IF ( buffer(1,3) /= 140.3_wp ) THEN
        WRITE ( 6 , * ) "atmosphere CPL ice 3 ", buffer(1,1)
        CALL test (.FALSE. )
     ENDIF
     IF ( buffer(1,4) /= 140.4_wp ) THEN
        WRITE ( 6 , * ) "atmosphere CPL ice 4 ", buffer(1,1)
        CALL test (.FALSE. )
     ENDIF
     IF ( buffer(1,5) /= 140.5_wp ) THEN
        WRITE ( 6 , * ) "atmosphere CPL ice 5 ", buffer(1,1)
        CALL test (.FALSE. )
     ENDIF
  ENDIF

  ! for target fields, yac_ftest should always return true/1
  CALL yac_ftest ( field_id(10), lflag )
  CALL test ( lflag )
  CALL yac_ftest ( field_id(10), iflag )
  CALL test ( iflag == 1 )

  DEALLOCATE(buffer)
  DEALLOCATE(field_id)

END SUBROUTINE dummy_atmosphere

! ------------------------------------------------------------------------

SUBROUTINE dummy_ocean ( )

  USE mpi
  USE MOD_NAME
  USE mo_yac_finterface

  IMPLICIT NONE

  INTEGER :: i, ierror
  INTEGER :: local_comm, npes, rank

  ! Inform the coupler about what we are
  comp_name = "dummy_ocean"
  CALL yac_fdef_comp ( comp_name, comp_id )
  comp_ids(1) = comp_id

  CALL yac_fget_comp_comm ( comp_id, local_comm )

  CALL mpi_comm_rank ( local_comm, rank, ierror )
  CALL mpi_comm_size ( local_comm, npes, ierror )

  ! WRITE ( 6 , * ) TRIM(comp_name), " rank ", rank, ": local size is ", npes

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

  grid_name = "dummy_ocean_grid"

  buffer_lon(1) =  0.0; buffer_lat(1) =  1.0
  buffer_lon(2) = -1.0; buffer_lat(2) =  0.0
  buffer_lon(3) =  1.0; buffer_lat(3) =  0.0
  buffer_lon(4) =  0.0; buffer_lat(4) = -1.0

  ! Connectivity
  cell_to_vertex(1,1) = 1; cell_to_vertex(2,1) = 2; cell_to_vertex(3,1) = 3 ! cell 1
  cell_to_vertex(1,2) = 2; cell_to_vertex(2,2) = 4; cell_to_vertex(3,2) = 3 ! cell 2

  ! Define unstructured grid
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
     core_mask(i) = 1
  ENDDO

  CALL yac_fset_global_index ( &
       & glb_index,            &
       & YAC_LOCATION_CELL,    &
       & grid_id )

  CALL yac_fset_core_mask ( &
       & glb_index,          &
       & YAC_LOCATION_CELL,  &
       & grid_id )

  ! Center points in cells (needed e.g. for nearest neighbour)

  buffer_lon(1) = 0.0; buffer_lat(1) =  0.5
  buffer_lon(2) = 0.0; buffer_lat(2) = -0.5

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

  CALL yac_fdef_mask(       &
       & grid_id,           &
       & nbr_cells,         &
       & YAC_LOCATION_CELL, &
       & cell_mask,         &
       & cell_mask_ids(1) );

  CALL yac_fset_mask ( &
       & cell_mask,    &
       & cell_point_ids(1) )

  DEALLOCATE (cell_mask)

  CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)

  field_name(1) =  "TAUX"   ! bundled field containing two components
  field_name(2) =  "TAUY"   ! bundled field containing two components
  field_name(3) =  "SFWFLX" ! bundled field containing three components
  field_name(4) =  "SFTEMP"
  field_name(5) =  "THFLX"  ! bundled field containing four components
  field_name(6) =  "ICEATM" ! bundled field containing four components
  field_name(7) =  "SST"
  field_name(8) =  "OCEANU"
  field_name(9) =  "OCEANV"
  field_name(10) = "ICEOCE" ! bundled field containing four components

  ref_collection_size( 1) = 2
  ref_collection_size( 2) = 2
  ref_collection_size( 3) = 3
  ref_collection_size( 4) = 1
  ref_collection_size( 5) = 4
  ref_collection_size( 6) = 4
  ref_collection_size( 7) = 1
  ref_collection_size( 8) = 1
  ref_collection_size( 9) = 1
  ref_collection_size(10) = 5

  ALLOCATE(field_id(no_of_fields))

  DO i = 1, no_of_fields/2
     CALL yac_fdef_field (          &
          & field_name(i),          &
          & comp_id,                &
          & cell_point_ids,         &
          & 1,                      &
          & ref_collection_size(i), &
          & "3600",                 &
          & YAC_TIME_UNIT_SECOND,   &
          & field_id(i) )
  ENDDO

  DO i = no_of_fields/2+1, no_of_fields
     CALL yac_fdef_field_mask (     &
          & field_name(i),          &
          & comp_id,                &
          & cell_point_ids,         &
          & cell_mask_ids,          &
          & 1,                      &
          & ref_collection_size(i), &
          & "3600",                 &
          & YAC_TIME_UNIT_SECOND,   &
          & field_id(i) )
  ENDDO

  ! read couplings from YAML configuration file
  CALL yac_fread_config_yaml(yaml_filename)

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
  !
  ! Send fields from ocean to atmosphere
  ! ------------------------------------
  !
  ! SST
  buffer(:,1) = 110.0_wp
  CALL yac_fput ( field_id(7), nbr_cells, 1, buffer(1:nbr_cells,1:1), info, ierror )

  ! zonal velocity
  buffer(:,1) = 120.0_wp
  CALL yac_fput ( field_id(8), nbr_cells, 1, buffer(1:nbr_cells,1:1), info, ierror )

  ! meridional velocity
  buffer(:,1) = 130.0_wp
  CALL yac_fput ( field_id(9), nbr_cells, 1, buffer(1:nbr_cells,1:1), info, ierror )

  ! Ice thickness, concentration, T1 and T2
  buffer(:,1) = 140.1_wp
  buffer(:,2) = 140.2_wp
  buffer(:,3) = 140.3_wp
  buffer(:,4) = 140.4_wp
  buffer(:,5) = 140.5_wp

  CALL yac_fput ( field_id(10), nbr_cells, 5, buffer(1:nbr_cells,1:5), info, ierror )

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

END SUBROUTINE dummy_ocean

! ------------------------------------------------------------------------

SUBROUTINE dummy_io ( )

  USE mpi
  USE mo_yac_finterface
  USE MOD_NAME

  IMPLICIT NONE

  INTEGER :: ierror
  INTEGER :: local_comm, npes, rank

  ! Inform the coupler about what we are
  comp_name = "dummy_io"
  CALL yac_fdef_comp ( comp_name, comp_id )

  CALL yac_fget_comp_comm ( comp_id, local_comm )

  CALL mpi_comm_rank ( local_comm, rank, ierror )
  CALL mpi_comm_size ( local_comm, npes, ierror )

  ! WRITE ( 6 , * ) TRIM(comp_name), " rank ", rank, ": local size is ", npes

  CALL yac_fenddef ( )

END SUBROUTINE dummy_io

END PROGRAM dummy_coupling
