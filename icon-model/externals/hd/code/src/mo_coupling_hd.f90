! mo_coupling_hd.f90 - Register HD to the coupler and access basic coupler functionality
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Author: Moritz Hanke (DKRZ)
! Contact: <stefan.hagemann@hereon.de>, <moritz.hanke@dkrz.de>
!_________________________________________
!----------------------------------
!>
!! Routines to access basic coupler functionality and routines for registering
!! the HD model to the coupler
!!
!! @par Revision History
!! First version by Moritz Hanke,  DKRZ, May 2022.
!!
!----------------------------

#ifdef NOMPI
#undef COUP_OAS
#undef COUP_YAC
#endif

#if defined(COUP_OAS) && defined(COUP_YAC)
#error "HD was configured with OASIS and YAC support. You can only use one of them."
#endif

#if defined(COUP_OAS) || defined(COUP_YAC)
#define COUPLED_RUN
#endif

MODULE mo_coupling_hd

#ifdef COUP_OAS
  USE oas_hd, ONLY: oas_hd_define, hd_receive_fld, hd_send_fld, hd_send_fld_nemo
#endif

  
#ifdef COUP_YAC
  USE mo_yac_finterface
  USE mpi
  USE mo_coupling, ONLY: runoff_s, runoff_dr, fdir_hd, get_grid_dimensions, &
                         lcoupling_atm, lcoupling_oce, icpl_mask_tohd
  USE mo_exception,ONLY: finish, message, message_text
  USE mo_io,       ONLY: io_open, io_close, io_read
  USE mo_netcdf,   ONLY: file_info, io_inq_varid, io_get_var_int, io_get_var_double
  USE mo_time_control, ONLY: delta_time
  USE mo_util_string,  ONLY: int2string
#else
  USE mo_coupling, ONLY: get_grid_dimensions, lcoupling_oce
#endif

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: coupling_hd_init
  PUBLIC :: coupling_hd_recv_from_land
  PUBLIC :: coupling_hd_send_to_ocean
  PUBLIC :: coupling_hd_send_to_ocean_direct

#ifdef COUP_YAC
  CHARACTER(LEN=*), PARAMETER :: comp_name = "HD"
  CHARACTER(LEN=*), PARAMETER :: grid_name = "HD_GRID"
  CHARACTER(LEN=*), PARAMETER :: land_comp_name = "LAND"
  CHARACTER(LEN=*), PARAMETER :: ocean_comp_name = "OCEAN"

  INTEGER :: nbr_cells(2) ! number of cells in lon and lat direction of
                          ! the local grid part assigned to this process

  INTEGER :: runoff_s_field_id ! surface water runoff; sum over forecast [kg/m2]
  INTEGER :: runoff_g_field_id ! soil water runoff; sum over forecast [kg/m2]
  INTEGER :: rdc2ocn_field_id ! river discharge to ocean [m3/s]

  INTEGER, PARAMETER :: nlev = 1
#endif

  REAL(dp), ALLOCATABLE, PUBLIC :: hd_receive_mask(:,:) ! mask on which HD may receive input used for icpl_mask_tohd=2
  REAL(dp), ALLOCATABLE, PUBLIC :: hd_outflow(:,:)      ! Discharge field send to ocean model, i.e. at river mouths

  CHARACTER(*), PARAMETER :: modname = "mo_coupling_hd"

CONTAINS

#ifdef COUP_YAC

  SUBROUTINE send_HD_domain_to_land(nlon, nlat, lon_org, lat_org, res)

    INTEGER, INTENT(IN) :: nlon, nlat
    REAL(dp), INTENT(IN) :: lon_org, lat_org, res

    INTEGER :: comm, root, comm_rank, comm_size, ierror
    REAL(dp) :: send_buffer(4)
    CHARACTER (LEN=YAC_MAX_CHARLEN) :: comp_names(2)

    comp_names(1) = TRIM(comp_name)
    comp_names(2) = TRIM(land_comp_name)
    CALL yac_fget_comps_comm(comp_names, 2, comm)

    send_buffer(1) = lon_org
    send_buffer(2) = lat_org
    send_buffer(3) = lon_org + REAL(nlon, dp) * res
    send_buffer(4) = lat_org - REAL(nlat, dp) * res

    CALL MPI_Comm_rank(comm, comm_rank, ierror)
    CALL MPI_Allreduce(comm_rank, root, 1, MPI_INTEGER, MPI_MIN, comm, ierror)
    CALL MPI_Bcast(send_buffer, 4, MPI_DOUBLE_PRECISION, root, comm, ierror)

  END SUBROUTINE send_HD_domain_to_land

#endif

  SUBROUTINE coupling_hd_init()

    INTEGER :: nlon
    INTEGER :: nlat
    REAL(dp) :: lon_org ! north-west corner longitude of gridbox(1,1)
    REAL(dp) :: lat_org ! north-west corner latitude of gridbox(1,1)
    REAL(dp) :: res     ! grid resolution

#ifdef COUP_YAC
    ! PI / 180
    REAL(dp), PARAMETER :: deg_to_rad = 0.01745329251994329576923690768489_dp

    INTEGER :: comp_id
    INTEGER :: grid_id
    INTEGER :: cell_point_id
    INTEGER :: land_mask_id
    INTEGER :: river_mouth_mask_id

    INTEGER :: nbr_vertices(2) ! number of vertices in lon and lat direction of
                               ! the local grid part assigned to this process

    INTEGER :: cyclic(2) ! ignore for now
    REAL, ALLOCATABLE :: x_vertices(:) ! longitude coordinates of the local grid
                                       ! part assigned to this process
    REAL, ALLOCATABLE :: y_vertices(:) ! latitude coordinates of the local grid
                                       ! part assigned to this process
    REAL, ALLOCATABLE :: x_cell_centers(:) ! longitude coordinates of the cell
                                           ! centers of the local grid part
                                           ! assigned to this process
    REAL, ALLOCATABLE :: y_cell_centers(:) ! latitude coordinates of the cell
                                           ! centers of the local grid part
                                           ! assigned to this process
    LOGICAL, ALLOCATABLE :: land_mask(:) ! mask for all HD cells (e.g. land) that receive water from atmosphere/land
    LOGICAL, ALLOCATABLE :: river_mouth_mask(:) ! mask for all river mouth cells
    INTEGER, ALLOCATABLE :: glb_cell_index(:)
    INTEGER, ALLOCATABLE :: glb_corner_index(:)

    ! get the timestep from HD configuration
    CHARACTER(LEN=80) :: timestep  	!delta_time

    ! TODO: set correct time unit, so that YAC can correctly interpret
    !       the timestep value
    !       (yac_time_unit_millisecond, yac_time_unit_second,
    !        yac_time_unit_minute, ..., or yac_time_unit_iso_format)
    INTEGER :: time_unit = yac_time_unit_second

    INTEGER :: ierror, i, j, k
#endif

    CALL get_grid_dimensions(nlon, nlat, lon_org, lat_org, res)

    IF (lcoupling_oce) THEN
      ALLOCATE (hd_outflow(nlon,nlat))
    ENDIF

#ifdef COUP_OAS

    CALL oas_hd_define()

#endif

#ifdef COUP_YAC
    timestep = int2string(int(delta_time))
    timestep = TRIM(ADJUSTL(timestep))

    lon_org = lon_org * deg_to_rad
    lat_org = lat_org * deg_to_rad
    res = res * deg_to_rad

    ! register HD component in YAC
    CALL yac_fdef_comp ( comp_name, comp_id )

    ! send HD domain to other component
    ! CALL send_HD_domain_to_land(nlon ,nlat, lon_org, lat_org, res)

    nbr_vertices(1) = nlon + 1
    nbr_vertices(2) = nlat + 1

    nbr_cells(1) = nlon
    nbr_cells(2) = nlat

    cyclic(1) = 0
    cyclic(2) = 0

    ALLOCATE(x_vertices(nbr_vertices(1)))
    ALLOCATE(y_vertices(nbr_vertices(2)))

    DO i = 1, nbr_vertices(1)
      x_vertices(i) = REAL(i - 1, dp) * res + lon_org
    END DO
    DO i = 1, nbr_vertices(2)
      y_vertices(i) = lat_org - REAL(i - 1, dp) * res
    END DO

    ! register HD grid in YAC
    CALL yac_fdef_grid ( &
         & grid_name,    &
         & nbr_vertices, &
         & cyclic,       &
         & x_vertices,   &
         & y_vertices,   &
         & grid_id )

    DEALLOCATE(y_vertices)
    DEALLOCATE(x_vertices)

    ALLOCATE(glb_cell_index(nbr_cells(1) * nbr_cells(2)))
    ALLOCATE(glb_corner_index(nbr_vertices(1) * nbr_vertices(2)))

    DO i = 1, nbr_cells(1) * nbr_cells(2)
      glb_cell_index(i) = i - 1
    END DO
    DO i = 1, nbr_vertices(1) * nbr_vertices(2)
      glb_corner_index(i) = i - 1
    END DO

    ! set global cell/corner ids
    CALL yac_fset_global_index ( &
         & glb_cell_index,       &
         & YAC_LOCATION_CELL,    &
         & grid_id )
    CALL yac_fset_global_index ( &
         & glb_corner_index,     &
         & YAC_LOCATION_CORNER,  &
         & grid_id )

    DEALLOCATE(glb_cell_index)
    DEALLOCATE(glb_corner_index)

! MoHa: the following lines are required in case the local field data
!       contains halo cells which do not contain valid data

!    ! set cell/corner core masks
!    CALL yac_fset_core_mask ( &
!         & cell_core_mask,    &
!         & YAC_LOCATION_CELL, &
!         & grid_id )
!    CALL yac_fset_core_mask (   &
!         & corner_core_mask,    &
!         & YAC_LOCATION_CORNER, &
!         & grid_id )

    ALLOCATE(x_cell_centers(nbr_cells(1)))
    ALLOCATE(y_cell_centers(nbr_cells(2)))

    DO i = 1, nbr_cells(1)
      x_cell_centers(i) = (REAL(i, dp) - 0.5_dp) * res + lon_org
    END DO
    DO i = 1, nbr_cells(2)
      y_cell_centers(i) = lat_org - (REAL(i - 1, dp) + 0.5) * res
    END DO

    ! define coordinates for the field data
    CALL yac_fdef_points (    &
         & grid_id,           &
         & nbr_cells,         &
         & YAC_LOCATION_CELL, &
         & x_cell_centers,    &
         & y_cell_centers,    &
         & cell_point_id )

    DEALLOCATE(y_cell_centers)
    DEALLOCATE(x_cell_centers)

    ALLOCATE(land_mask(nbr_cells(1) * nbr_cells(2)))
    ALLOCATE(river_mouth_mask(nbr_cells(1) * nbr_cells(2)))

    ! Read hd_receive_mask if icpl_mask_tohd = 2
    IF (icpl_mask_tohd.EQ.2) CALL read_receiving_info(nlon, nlat, 'hd_receive.nc')

    k = 1
    DO i = 1, nlat
      DO j = 1, nlon
        IF (icpl_mask_tohd.EQ.0) land_mask(k) = fdir_hd(j,i) > 0
        IF (icpl_mask_tohd.EQ.2) THEN
          land_mask(k) = hd_receive_mask(j,i) > 0._dp .OR. fdir_hd(j,i) > 0
          river_mouth_mask(k) = fdir_hd(j,i) == 0 .OR.  &
                                (hd_receive_mask(j,i) > 0._dp .AND. fdir_hd(j,i) == -1)
        ELSE
          river_mouth_mask(k) = fdir_hd(j,i) == 0
        ENDIF
        k = k + 1
      END DO
    END DO
    IF (icpl_mask_tohd.EQ.1) land_mask(:) = .TRUE.

    CALL yac_fdef_mask( &
      grid_id, nbr_cells(1) * nbr_cells(2), YAC_LOCATION_CELL, &
      land_mask, land_mask_id)
    CALL yac_fdef_mask( &
      grid_id, nbr_cells(1) * nbr_cells(2), YAC_LOCATION_CELL, &
      river_mouth_mask, river_mouth_mask_id)

    !---------------------
    ! define output fields
    !---------------------
    ! river discharge to ocean
    IF (lcoupling_oce) THEN
      CALL yac_fdef_field_mask (      &
         & 'river_runoff',          &
         & comp_id,                 &
         & (/cell_point_id/),       &
         & (/river_mouth_mask_id/), &
         & 1,                       &
         & nlev,                    &
         & timestep,                &
         & time_unit,               &
         & rdc2ocn_field_id )
    ENDIF

    !---------------------
    ! define input fields
    !---------------------

    IF (lcoupling_atm) THEN
      ! surface water runoff; sum over forecast
      CALL yac_fdef_field_mask ( &
         & 'surface_water_runoff',         &
         & comp_id,            &
         & (/cell_point_id/),  &
         & (/land_mask_id/),   &
         & 1,                  &
         & nlev,               &
         & timestep,           &
         & time_unit,          &
         & runoff_s_field_id )

      ! soil water runoff = subsurface runoff = drainage; sum over forecast
      CALL yac_fdef_field_mask ( &
         & 'soil_water_runoff',  &
         & comp_id,            &
         & (/cell_point_id/),  &
         & (/land_mask_id/),   &
         & 1,                  &
         & nlev,               &
         & timestep,           &
         & time_unit,          &
         & runoff_g_field_id )
    ENDIF

    ! setup coupling (incl. weight computation)
    ! this call is collective for all processes
    CALL yac_fenddef ( )
#endif

  END SUBROUTINE coupling_hd_init

  SUBROUTINE coupling_hd_recv_from_land()

#ifdef COUP_OAS

    CALL hd_receive_fld()

#endif

#ifdef COUP_YAC
    ! Note that the unit conversion to m/s for both runoff fluxes is done afterwards using
    ! the unit factor ufakru from the namelist HD_CTL

    ! 1: RUNOFF_S (kg/m2/forecast time): surface water runoff; sum over forecast or average flux
    CALL get_field(runoff_s, runoff_s_field_id)

    ! 2: RUNOFF_G (kg/m2/forecast time): soil water runoff; sum over forecast or average flux
    CALL get_field(runoff_dr, runoff_g_field_id)
#endif

  END SUBROUTINE coupling_hd_recv_from_land

#ifdef COUP_YAC
  SUBROUTINE get_field(field, field_id)

    REAL(dp), TARGET :: field(*)
    INTEGER :: field_id

    TYPE(yac_dble_ptr) :: recv_field(nlev)
    INTEGER :: i
    INTEGER :: nbr_hor_points, level_start, level_end
    INTEGER :: info, ierror

    nbr_hor_points = nbr_cells(1) * nbr_cells(2)

    DO i = 1, nlev
      level_start = (i-1) * nbr_hor_points + 1
      level_end = i * nbr_hor_points
      recv_field(i)%p => field(level_start:level_end)
    END DO

    CALL yac_fget(field_id, nlev, recv_field, info, ierror)

  END SUBROUTINE
#endif

  SUBROUTINE coupling_hd_send_to_ocean(friv_hd)

    REAL(dp), INTENT(IN) :: friv_hd(:,:) ! river flow

#ifdef COUP_OAS
    CALL hd_send_fld(friv_hd)
#endif

#ifdef COUP_YAC
    ! RDC2NEMO [m3/s]: river discharge to ocean
    CALL put_field(friv_hd, rdc2ocn_field_id)
#endif

  END SUBROUTINE coupling_hd_send_to_ocean

  SUBROUTINE coupling_hd_send_to_ocean_direct()

    CHARACTER(*), PARAMETER :: routine = modname//":coupling_hd_send_to_ocean_direct"

#ifdef COUP_OAS
    CALL hd_send_fld_nemo()
#endif

#ifdef COUP_YAC
    CALL yac_abort_message( &
        routine // " : this routine is not implented for YAC", &
        __FILE__, __LINE__)
#endif

  END SUBROUTINE coupling_hd_send_to_ocean_direct

#ifdef COUP_YAC
  SUBROUTINE put_field(field, field_id)

    REAL(dp), TARGET :: field(*)
    INTEGER :: field_id

    TYPE(yac_dble_ptr) :: send_field(1, nlev)
    INTEGER :: i
    INTEGER :: nbr_hor_points, level_start, level_end
    INTEGER :: info, ierror

    nbr_hor_points = nbr_cells(1) * nbr_cells(2)

    DO i = 1, nlev
      level_start = (i-1) * nbr_hor_points + 1
      level_end = i * nbr_hor_points
      send_field(1,i)%p => field(level_start:level_end)
    END DO

    CALL yac_fput(field_id, 1, nlev, send_field, info, ierror)

  END SUBROUTINE put_field
#endif

#ifdef COUP_YAC
  SUBROUTINE read_receiving_info(nlon, nlat, receiving_file)
    !
    ! Reads mask on HD grid with comprise boxes that may receive input from atmosphere/land model.
    ! Due to different grids and interpolation, these boxes may also comprise ocean boxes on HD grid.
    !
    !
    INTEGER, INTENT(in) :: nlon
    INTEGER, INTENT(in) :: nlat
    CHARACTER(len=*), INTENT(in) :: receiving_file

    TYPE (FILE_INFO)  :: fileinfo

    INTEGER nvarid, fileid
    LOGICAL :: lex

    INQUIRE (file=TRIM(receiving_file), exist=lex)
    IF (.NOT. lex) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(receiving_file),'>'
      CALL message('read_receiving_info', message_text)
      CALL finish ('read_receiving_info', 'run terminated.')
    ENDIF

    fileinfo%opened = .FALSE.
    CALL IO_open (TRIM(receiving_file), fileinfo, IO_READ)
    CALL message('', '')
    WRITE (message_text,*) 'Reading receiving info from file ', TRIM(receiving_file)
    CALL message('read_receiving_info', message_text)

    fileID = fileinfo%file_id
    !
    ! Allocate and read hd_receive_mask 
    ALLOCATE (hd_receive_mask(nlon,nlat))
    CALL IO_inq_varid (fileID, 'hd_receive_mask', nvarid)
    CALL IO_get_var_double (fileID, nvarid, hd_receive_mask)    

    CALL IO_close(fileinfo)

  END SUBROUTINE read_receiving_info
#endif

! **************************************************************************
  SUBROUTINE halo_around_mouth(nlon, nlat, lglobe, fdir_hd, mouth_with_halo)
! **************************************************************************
    ! Routine is not tested up to now.

    INTEGER, INTENT(in) :: nlon
    INTEGER, INTENT(in) :: nlat
    LOGICAL, INTENT(in) :: lglobe                ! Global/cyclic data for longitudes No/Yes=F/T
    INTEGER, INTENT(in) :: fdir_hd(:,:)          ! river flow directions
    INTEGER, INTENT(out) :: mouth_with_halo(:,:) ! mouth mask plus halo of one ocean (fdir=-1) grid box

    INTEGER :: jl, jb

    mouth_with_halo(:,:) = 0
    WHERE (fdir_hd(:,:) == 0)
      mouth_with_halo(:,:) = 1
    END WHERE    

    DO jb = 2, nlat-1
      DO jl = 2, nlon-1
        IF (fdir_hd(jl,jb) == 0) THEN
          WHERE (fdir_hd(jl-1:jl+1, jb-1:jb+1) == -1)  
            mouth_with_halo(jl-1:jl+1, jb-1:jb+1) = 1
          END WHERE    
        ENDIF
!!        river_mouth_mask(k) = fdir_hd(j,i) == 0
      END DO
      ! Western Border
      IF (fdir_hd(1,jb) == 0) THEN
        WHERE (fdir_hd(1:2, jb-1:jb+1) == -1)  
          mouth_with_halo(1:2, jb-1:jb+1) = 1
        END WHERE    
        IF (lglobe) THEN
          WHERE (fdir_hd(nlon, jb-1:jb+1) == -1)  
            mouth_with_halo(nlon, jb-1:jb+1) = 1
          END WHERE    
        ENDIF
      ENDIF
      ! Eastern Border
      IF (fdir_hd(nlon,jb) == 0) THEN
        WHERE (fdir_hd(nlon-1:nlon, jb-1:jb+1) == -1)  
          mouth_with_halo(nlon-1:nlon, jb-1:jb+1) = 1
        END WHERE    
        IF (lglobe) THEN
          WHERE (fdir_hd(1, jb-1:jb+1) == -1)  
            mouth_with_halo(1, jb-1:jb+1) = 1
          END WHERE    
        ENDIF
      ENDIF
    END DO
!
!   *** Northern and Southern latitude
    DO jl = 2, nlon-1
      IF (fdir_hd(jl,1) == 0) THEN
        WHERE (fdir_hd(jl-1:jl+1, 1:2) == -1)  
          mouth_with_halo(jl-1:jl+1, 1:2) = 1
        END WHERE    
      ENDIF
      IF (fdir_hd(jl,nlon) == 0) THEN
        WHERE (fdir_hd(jl-1:jl+1, nlat-1:nlat) == -1)  
          mouth_with_halo(jl-1:jl+1, nlat-1:nlat) = 1
        END WHERE    
      ENDIF
    END DO
   ! Western Border
    IF (fdir_hd(1,1) == 0) THEN    ! North
      WHERE (fdir_hd(1:2, 1:2) == -1)  
        mouth_with_halo(1:2, 1:2) = 1
      END WHERE    
      IF (lglobe) THEN
        IF (fdir_hd(nlon, 1) == -1) mouth_with_halo(nlon, 1) = 1
        IF (fdir_hd(nlon, 2) == -1) mouth_with_halo(nlon, 2) = 1
      ENDIF
    ENDIF
    IF (fdir_hd(1,nlat) == 0) THEN  ! South
      WHERE (fdir_hd(1:2, nlat-1:nlat) == -1)  
        mouth_with_halo(1:2, nlat-1:nlat) = 1
      END WHERE    
      IF (lglobe) THEN
        IF (fdir_hd(nlon, nlat-1) == -1) mouth_with_halo(nlon, nlat-1) = 1
        IF (fdir_hd(nlon, nlat) == -1) mouth_with_halo(nlon, nlat) = 1
      ENDIF
    ENDIF

    ! Eastern Border
    IF (fdir_hd(nlon,1) == 0) THEN   ! North
      WHERE (fdir_hd(nlon-1:nlon, 1:2) == -1)  
        mouth_with_halo(nlon-1:nlon, 1:2) = 1
      END WHERE    
      IF (lglobe) THEN
        IF (fdir_hd(1, 1) == -1) mouth_with_halo(1, 1) = 1
        IF (fdir_hd(1, 2) == -1) mouth_with_halo(1, 2) = 1
      ENDIF
    ENDIF
    IF (fdir_hd(nlon,nlat) == 0) THEN  ! South
      WHERE (fdir_hd(nlon-1:nlon, nlat-1:nlat) == -1)  
        mouth_with_halo(nlon-1:nlon, nlat-1:nlat) = 1
      END WHERE    
      IF (lglobe) THEN
        IF (fdir_hd(1, nlat-1) == -1) mouth_with_halo(1, nlat-1) = 1
        IF (fdir_hd(1, nlat) == -1) mouth_with_halo(1, nlat) = 1
      ENDIF
    ENDIF

  END SUBROUTINE halo_around_mouth

END MODULE mo_coupling_hd
