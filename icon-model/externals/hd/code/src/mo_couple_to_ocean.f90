! mo_couple_to_ocean.f90 - Routine for direct coupling to an ocean model
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Author: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

MODULE mo_couple_to_ocean

  !
  ! *** Module that includes routine for the direct coupling to the ocean, where
  !     where target gridboxes on the ocean grid are directly associated with 
  !     HD river mouth point within the HD model. For the coupling to the ocean model, 
  !     the related ocean inflow on the target grid can be passed without interpolation
  !     to the ocean model (e.g. NEMO). e.g. via OASIS.
  !
  ! Vs. 1.0 - S. Hagemann, HZG-IfK, November 2017, original source
  !


  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_grid,          ONLY: domain
  USE mo_time_control,  ONLY: get_date_components, current_date, initial_date 
  USE mo_time_conversion,  ONLY: TC_get 
  USE mo_hydrology,     ONLY: grid_hd
  USE mo_io,            ONLY: io_open, io_close, io_read, io_write
  USE mo_netcdf,        ONLY: nf_global, nf_double, nf_real, file_info, nf_enddef, &
                              io_inq_dimid, io_inq_varid, io_inq_dimlen, &
                              nf_def_dim, nf_def_var, nf_put_att_text, &
                              io_get_var_double, io_get_var_int,     &
                              nf_put_var_double, nf_put_vara_double, &
                              nf_get_att_int, nf_noerr, nf_close, nf_def_var_deflate,  &
                              nf_check, nf_netcdf4, nf_create, nf_unlimited

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_coupling_info, dis_to_ocean, hd_on_ocean_open, hd_on_ocean_write, hd_on_ocean_close

  ! On HD grid
  INTEGER, ALLOCATABLE, PUBLIC :: fmou_hd_to_ocean(:,:)     ! HD mouth points that have a target on ocean grid
                                                   ! 0= no target, 1= target grid box has only one source
                                                   ! 2 or larger = target grid boxes has several sources
  INTEGER, ALLOCATABLE :: indexx(:,:)              ! Target x-indices (longitude) of HD mouths on ocean grid
  INTEGER, ALLOCATABLE :: indexy(:,:)              ! Target y-indices (latitude) of HD mouths on ocean grid

  ! On ocean grid, e.g. NEMO
  INTEGER, PUBLIC       :: nxocean                 ! x dimension of ocean grid (Longitude)
  INTEGER, PUBLIC       :: nyocean                 ! y dimension of ocean grid (Latitude)
  INTEGER, ALLOCATABLE  :: fmou_hd_on_ocean(:,:)    ! Target points of HD mouths on ocean grid that specify number of sources
  REAL(dp), ALLOCATABLE, PUBLIC :: discharge_on_ocean(:,:)    ! Inflow into the ocean on ocean grid
  REAL(dp), ALLOCATABLE :: olon(:,:)               ! Longitudes on ocean grid
  REAL(dp), ALLOCATABLE :: olat(:,:)               ! Latitudes on ocean grid

  ! netCDF id of discharge_on_ocean used in the hd_on_ocean routines
  INTEGER ::  ncid_hdoo            ! File ID
  INTEGER ::  time_id              ! Time varable ID
  INTEGER ::  dis_on_ocean_id      ! Variable ID of output on ocean grid
  INTEGER, SAVE :: icount = 0

CONTAINS

  SUBROUTINE read_coupling_info(coupling_file)
    !
    ! reads necessary fields and ocean grid info for coupling to ocean using coupling_type 2
    !
    ! Default input file:  hdcouple.nc: 
    !
    TYPE (FILE_INFO)  :: fileinfo

    INTEGER nvarid, fileid, dimid, i, status
    CHARACTER(len=80) :: coupling_file
    LOGICAL :: lex

    INQUIRE (file=coupling_file, exist=lex)
    IF (.NOT. lex) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(coupling_file),'>'
      CALL message('read_coupling_info', message_text)
      CALL finish ('read_coupling_info', 'run terminated.')
    ENDIF

    fileinfo%opened = .FALSE.
    CALL IO_open (coupling_file, fileinfo, IO_READ)
    CALL message('', '')
    WRITE (message_text,*) 'Reading coupling info from file ', TRIM(coupling_file)
    CALL message('read_coupling_info', message_text)

    fileID = fileinfo%file_id
    !
    ! Allocate HD fields
    ALLOCATE (fmou_hd_to_ocean(grid_hd%nlon,grid_hd%nlat))
    ALLOCATE (indexx(grid_hd%nlon,grid_hd%nlat))
    ALLOCATE (indexy(grid_hd%nlon,grid_hd%nlat))
    !
    ! Read HD fields
    CALL IO_inq_varid (fileID, 'FMOU_HD_TO_NEMO', nvarid)
    CALL IO_get_var_int (fileID, nvarid, fmou_hd_to_ocean)    
    CALL IO_inq_varid (fileID, 'INDEXX', nvarid)
    CALL IO_get_var_int (fileID, nvarid, indexx)
    CALL IO_inq_varid (fileID, 'INDEXY', nvarid)
    CALL IO_get_var_int (fileID, nvarid, indexy)
    
    ! Read ocean dimensions
    CALL IO_inq_dimid (fileID, 'x', dimid)
    CALL IO_inq_dimlen (fileID, dimid, nxocean)
    CALL IO_inq_dimid (fileID, 'y', dimid)
    CALL IO_inq_dimlen (fileID, dimid, nyocean)
    ALLOCATE (fmou_hd_on_ocean(nxocean,nyocean))
    ALLOCATE (olon(nxocean,nyocean))
    ALLOCATE (olat(nxocean,nyocean))
    ALLOCATE (discharge_on_ocean(nxocean,nyocean))

    ! Read ocean fields
    CALL IO_inq_varid (fileID, 'FMOU_HD_ON_NEMO', nvarid)
    CALL IO_get_var_int (fileID, nvarid, fmou_hd_on_ocean)    
    CALL IO_inq_varid (fileID, 'olon', nvarid)
    CALL IO_get_var_double (fileID, nvarid, olon)    
    CALL IO_inq_varid (fileID, 'olat', nvarid)
    CALL IO_get_var_double (fileID, nvarid, olat)

    CALL IO_close(fileinfo)

  END SUBROUTINE read_coupling_info

!
! *********************************************************************
  SUBROUTINE dis_to_ocean(hd_outflow)
! *********************************************************************
!
!     Transfer HD discharge to ocean grid, e.g. NEMO


      REAL(dp), INTENT(in) :: hd_outflow(:,:)
      INTEGER :: jl, jb
!
      discharge_on_ocean(:,:) = 0._dp
      DO jb=1, grid_hd%nlat
      DO jl=1, grid_hd%nlon 
      IF (indexx(jl,jb).GT.0.5) THEN
        discharge_on_ocean(indexx(jl,jb), indexy(jl,jb)) = &
          discharge_on_ocean(indexx(jl,jb), indexy(jl,jb)) + hd_outflow(jl,jb)
      ENDIF
      ENDDO	  
      ENDDO	  
!
  END SUBROUTINE
!
! *********************************************************************
  SUBROUTINE hd_on_ocean_open
! *********************************************************************

    USE mo_namelist,         ONLY: open_nml, position_nml, POSITIONED
    USE mo_io_units,         ONLY: nout

    ! dimension ids

    INTEGER ::  lon_dim
    INTEGER ::  lat_dim
    INTEGER ::  time_dim

    ! dimension lengths

    INTEGER, PARAMETER ::  time_len = NF_UNLIMITED

    ! variable ids

    INTEGER ::  lon_id
    INTEGER ::  lat_id

    ! variable shapes

    INTEGER :: dim1(1)
    INTEGER :: dim2(2)
    INTEGER :: dim3(3)

    INTEGER :: read_status, inml, iunit
    INTEGER :: year, month, day
    CHARACTER(len=80) :: dnoo
    CHARACTER(8) :: cdate     ! Date in reality = date of creating output file
    CHARACTER(30) :: cstart

    ! Namelist with user attributes
    CHARACTER(200) :: hd_user, hd_cont, hd_inst, hd_instid
    NAMELIST /HDUSER_CTL/ &
    hd_user, hd_cont, hd_inst, hd_instid

    dnoo = "discharge_on_ocean.nc"

    ! enter define mode

    CALL nf_check(nf_create(TRIM(dnoo), NF_NETCDF4, ncid_hdoo),fname=TRIM(dnoo))

    WRITE (message_text,*)                                    &
         'HD model river discharge output on ocean grid: ', TRIM(dnoo)
    CALL message('hd_on_ocean_open', message_text)

    ! read namelist hduser_ctl

    inml = open_nml ('namelist.hduser')
    iunit = position_nml ('HDUSER_CTL', inml, status=read_status)
    SELECT CASE (read_status)
    CASE (POSITIONED)
       READ (iunit, hduser_ctl)
       CALL message('hd_on_ocean_open', 'Namelist HDUSER_CTL: ')
       WRITE(nout, hduser_ctl)
    END SELECT

    ! define dimensions

    CALL nf_check(nf_def_dim(ncid_hdoo, 'x', nxocean, lon_dim))
    CALL nf_check(nf_def_dim(ncid_hdoo, 'y', nyocean, lat_dim))
    CALL nf_check(nf_def_dim(ncid_hdoo, 'time', NF_UNLIMITED, time_dim))

    ! define variables

    dim1(1) = time_dim
    CALL nf_check(nf_def_var(ncid_hdoo, 'time', NF_DOUBLE, 1, dim1, time_id))
    CALL nf_check(nf_def_var_deflate(ncid_hdoo, time_id, 0, 1, 2))  ! shuffle, deflate, dflate_level

    dim2(1:2) = (/ lon_dim, lat_dim /)
    CALL nf_check(nf_def_var(ncid_hdoo, 'lon', NF_DOUBLE, 2, dim2, lon_id))
    CALL nf_check(nf_def_var_deflate(ncid_hdoo, lon_id, 0, 1, 2))  ! shuffle, deflate, dflate_level
    CALL nf_check(nf_def_var(ncid_hdoo, 'lat', NF_DOUBLE, 2, dim2, lat_id))
    CALL nf_check(nf_def_var_deflate(ncid_hdoo, lat_id, 0, 1, 2))  ! shuffle, deflate, dflate_level
            
    dim3(1:3) = (/ lon_dim, lat_dim, time_dim /)
    CALL nf_check(nf_def_var(ncid_hdoo, 'discharge_on_ocean', NF_REAL, 3, dim3, dis_on_ocean_id))
    CALL nf_check(nf_def_var_deflate(ncid_hdoo, dis_on_ocean_id, 0, 1, 2))  ! shuffle, deflate, dflate_level

    ! assign attributes
    CALL nf_check(nf_put_att_text(ncid_hdoo, lon_id, 'long_name', 9, 'longitude'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, lon_id, 'units', 12, 'degrees_east'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, lon_id, 'standard_name', 9, 'longitude'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, lon_id, 'axis',          1, 'X'))

    CALL nf_check(nf_put_att_text(ncid_hdoo, lat_id, 'long_name', 8, 'latitude'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, lat_id, 'units', 13, 'degrees_north'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, lat_id, 'standard_name', 8, 'latitude'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, lat_id, 'axis',          1, 'Y'))

    ! dt_start = (yr, mo, dy, hr, mi, se) starting date of actual run, i.e. restart date.
    CALL get_date_components(initial_date, year, month, day)
    WRITE(cstart, '(A11,I4.4, 2(A1,I2.2),A9 )') 'days since ', year, &
          '-',month, '-',day, ' 00:00:00'   
    CALL nf_check(nf_put_att_text(ncid_hdoo, time_id, 'units',       30, cstart))
    CALL nf_check(nf_put_att_text(ncid_hdoo, time_id, 'calendar',19, 'proleptic_gregorian'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, time_id, 'axis',         1, 'T'))

    CALL nf_check(nf_put_att_text(ncid_hdoo, dis_on_ocean_id, 'long_name', 18, 'discharge on ocean'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, dis_on_ocean_id, 'standard_name', 45, 'water_volume_transport_into_ocean_from_rivers'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, dis_on_ocean_id, 'units', 6, 'm3 s-1'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, dis_on_ocean_id, 'coordinates', 9, 'lon lat'))
        
    CALL DATE_AND_TIME(cdate)          ! returns YYYYMMDD

    CALL nf_check(nf_put_att_text(ncid_hdoo, NF_GLOBAL, 'Conventions',    6, 'CF-1.6'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, NF_GLOBAL, 'source',  42, 'HD model Vs. 5, doi:10.5281/zenodo.4893099'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, NF_GLOBAL, 'title',         32, 'HD model river inflow into ocean'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, NF_GLOBAL, 'licence',        9, 'CC-BY 4.0'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, NF_GLOBAL, 'crs',           16, 'WGS84, EPSG:4326'))
    CALL nf_check(nf_put_att_text(ncid_hdoo, NF_GLOBAL, 'creation_date', 10, cdate(7:8)//'.'//cdate(5:6)//'.'//cdate(1:4) ))

    ! Institution and person related attributes
    CALL nf_check(nf_put_att_text(ncid_hdoo, NF_GLOBAL, 'institution', len(TRIM(hd_inst)), &
         TRIM(hd_inst) ))
    CALL nf_check(nf_put_att_text(ncid_hdoo, NF_GLOBAL, 'institution-ID', len(TRIM(hd_instid)), &
         TRIM(hd_instid) ))
    CALL nf_check(nf_put_att_text(ncid_hdoo, NF_GLOBAL, 'contact', len(TRIM(hd_cont)), TRIM(hd_cont) ))
    CALL nf_check(nf_put_att_text(ncid_hdoo, NF_GLOBAL, 'creator', len(TRIM(hd_user)), TRIM(hd_user) ))

    ! leave define mode
    CALL nf_check(nf_enddef(ncid_hdoo))

    CALL nf_check(nf_put_var_double(ncid_hdoo, lat_id, olat))
    CALL nf_check(nf_put_var_double(ncid_hdoo, lon_id, olon))

  END SUBROUTINE hd_on_ocean_open

!===========================================================================
  SUBROUTINE hd_on_ocean_close

    CALL nf_check(nf_close(ncid_hdoo))

  END SUBROUTINE hd_on_ocean_close

!===========================================================================

  SUBROUTINE hd_on_ocean_write(friv_on_ocean)

    REAL(dp), INTENT(in) :: friv_on_ocean(:,:)

    ! starts and counts for array sections of record variables

    INTEGER ::  tstart(1), tcount(1)
    INTEGER ::  fstart(3), fcount(3)
    INTEGER  :: day, second, day1, sec1
    REAL(dp) :: days_since_start(1)

    icount = icount+1

    ! Calculate days since start 
    CALL TC_get(initial_date, day1, sec1)
    CALL TC_get(current_date, day, second)
    days_since_start(1) = day - day1 + (second-sec1)/86400._dp

    tstart(1) = icount
    tcount(1) = 1
    CALL nf_check(nf_put_vara_double(ncid_hdoo, time_id, tstart, tcount, days_since_start))

    fstart(1:3) = (/ 1, 1, icount /)
    fcount(1:3) = (/ nxocean, nyocean, 1 /)
    CALL nf_check(nf_put_vara_double(ncid_hdoo, dis_on_ocean_id, fstart, fcount, friv_on_ocean))

  END SUBROUTINE hd_on_ocean_write
  
END MODULE mo_couple_to_ocean

