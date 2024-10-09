! mo_hd_highres_io.f90 - HD model IO routines on HD grid
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann and Ha Ho-Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

MODULE mo_hd_highres_io

! Modifications
!
! S. Hagemann, MPI-M, Oct 1999 : original source
! MPI-M members: Restructuring for the use in MPI-ESM, the Earth System Model of the 
!    Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
!    Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
!    version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
!    doi: 10.1029/2018MS001400.
! H. Ho-Hagemann, HZG, June 2014, introducing OASIS coupling in offline version.
! S. Hagemann, MPI-M, 2015, further refinements for offline version.
! H. Ho-Hagemann, HZG, 2018, Further Adaptations for OASIS coupling.
! S. Hagemann, HZG, 2018, Small refinements for HD 5 Min. version.

  USE mo_kind,         ONLY: dp
  USE mo_exception,    ONLY: message_text, message
  USE mo_time_control, ONLY: get_date_components, next_date, current_date, initial_date 
  USE mo_time_conversion,  ONLY: TC_get
  USE mo_filename,     ONLY: compose_filenames, standard_output_file
  USE mo_grid,         ONLY: domain
  USE mo_netcdf,       ONLY: nf_check, nf_create, nf_close,          &
                             nf_unlimited, nf_put_vara_double,       & 
                             nf_def_dim, nf_put_att_text, nf_enddef, & 
                             nf_put_var_double, nf_def_var, nf_def_var_deflate,  &
                             nf_global, nf_double, nf_real, nf_clobber, nf_max_name, nf_netcdf4

  IMPLICIT NONE

  PRIVATE
  
  ! netCDF id
  INTEGER ::  ncid

  ! dimension ids

  INTEGER ::  lon_dim
  INTEGER ::  lat_dim
  INTEGER ::  time_dim

  ! HD grid definition
  TYPE(domain)  :: grid_hd    

  ! dimension lengths

  INTEGER, PARAMETER ::  time_len = NF_UNLIMITED

  ! variable ids

  INTEGER ::  lon_id, nlon_id
  INTEGER ::  lat_id, nlat_id
  INTEGER ::  time_id

!! #ifdef COUP_OAS
  INTEGER ::  hdarea_id
!! #endif
  INTEGER ::  friv_id

  ! variable shapes

  INTEGER :: dim1(1)
  INTEGER :: dim2(2)
  INTEGER :: dim3(3)
  
  INTEGER, SAVE :: icount = 0

  PUBLIC :: hd_highres_init
  PUBLIC :: hd_highres_open
  PUBLIC :: hd_highres_close
  PUBLIC :: hd_highres_write

  ! Namelist with user attributes
  CHARACTER(200) :: hd_user, hd_cont, hd_inst, hd_instid
  NAMELIST /HDUSER_CTL/ &
     hd_user, hd_cont, hd_inst, hd_instid

CONTAINS

  SUBROUTINE hd_highres_init (grid_in)

    TYPE(domain), INTENT(in)  :: grid_in

    grid_hd = grid_in

  END SUBROUTINE hd_highres_init

  SUBROUTINE hd_highres_open (filename, varname, ncid, klon, klat)

    USE mo_namelist,         ONLY: open_nml, position_nml, POSITIONED
    USE mo_io_units,         ONLY: nout

    CHARACTER(nf_max_name), INTENT(in)  ::  filename
    CHARACTER(nf_max_name), INTENT(in)  ::  varname
    INTEGER,                INTENT(out) ::  ncid
    INTEGER, OPTIONAL,      INTENT(in)  ::  klon
    INTEGER, OPTIONAL,      INTENT(in)  ::  klat

    ! data variables

    REAL(dp), ALLOCATABLE :: lon(:)
    REAL(dp), ALLOCATABLE :: lat(:)

    INTEGER :: jb, jl
    INTEGER :: year, month, day
    INTEGER :: read_status, inml, iunit

    CHARACTER(nf_max_name) :: total_filename
    CHARACTER(30) :: cstart
    CHARACTER(8) :: cdate     ! Date in reality = date of creating output file
    CHARACTER(120) :: ctitle

    CALL compose_filenames
    total_filename = TRIM(standard_output_file)//'_'//TRIM(filename)

    ! enter define mode

    CALL nf_check(nf_create(TRIM(total_filename), NF_NETCDF4, ncid),fname=TRIM(total_filename))

    WRITE (message_text,*)                                    &
         'HD model high resolution river discharge output: ', &
         TRIM(total_filename)
    CALL message('hd_highres_open', message_text)

    ! read namelist hduser_ctl

    inml = open_nml ('namelist.hduser')
    iunit = position_nml ('HDUSER_CTL', inml, status=read_status)
    SELECT CASE (read_status)
    CASE (POSITIONED)
       READ (iunit, hduser_ctl)
       CALL message('hd_highres_open', 'Namelist HDUSER_CTL: ')
       WRITE(nout, hduser_ctl)
    END SELECT

    ! define dimensions
    IF (PRESENT(klon)) THEN
      CALL nf_check(nf_def_dim(ncid, 'lon', klon, lon_dim))
      CALL nf_check(nf_def_dim(ncid, 'lat', klat, lat_dim))
    ELSE
      CALL nf_check(nf_def_dim(ncid, 'lon', grid_hd%nlon, lon_dim))
      CALL nf_check(nf_def_dim(ncid, 'lat', grid_hd%nlat, lat_dim))
    ENDIF
    CALL nf_check(nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim))

    ! define variables

    dim1(1) = time_dim
    CALL nf_check(nf_def_var(ncid, 'time', NF_DOUBLE, 1, dim1, time_id))
    CALL nf_check(nf_def_var_deflate(ncid, time_id, 0, 1, 2))  ! shuffle, deflate, dflate_level

    dim2(1:2) = (/ lon_dim, lat_dim /)
    ! ncview cannot deal with 2-D coordinates  
    CALL nf_check(nf_def_var(ncid, 'lon', NF_DOUBLE, 1, dim2(1), lon_id))
    CALL nf_check(nf_def_var(ncid, 'lat', NF_DOUBLE, 1, dim2(2), lat_id))

!! #ifdef COUP_OAS
    CALL nf_check(nf_def_var(ncid, 'hdarea', NF_DOUBLE, 2, dim2, hdarea_id))
    CALL nf_check(nf_def_var_deflate(ncid, hdarea_id, 0, 1, 2))  ! shuffle, deflate, dflate_level
!!  #endif
            
    dim3(1:3) = (/ lon_dim, lat_dim, time_dim /)
    CALL nf_check(nf_def_var(ncid, TRIM(varname), NF_REAL, 3, dim3, friv_id))
    CALL nf_check(nf_def_var_deflate(ncid, friv_id, 0, 1, 2))  ! shuffle, deflate, dflate_level

    ! assign attributes
    CALL nf_check(nf_put_att_text(ncid, lon_id, 'long_name',     9, 'longitude'))
    CALL nf_check(nf_put_att_text(ncid, lon_id, 'units',        12, 'degrees_east'))
    CALL nf_check(nf_put_att_text(ncid, lon_id, 'standard_name', 9, 'longitude'))
    CALL nf_check(nf_put_att_text(ncid, lon_id, 'axis',          1, 'X'))

    CALL nf_check(nf_put_att_text(ncid, lat_id, 'long_name',     8, 'latitude'))
    CALL nf_check(nf_put_att_text(ncid, lat_id, 'units',        13, 'degrees_north'))
    CALL nf_check(nf_put_att_text(ncid, lat_id, 'standard_name', 8, 'latitude'))
    CALL nf_check(nf_put_att_text(ncid, lat_id, 'axis',          1, 'Y'))

    ! dt_start = (yr, mo, dy, hr, mi, se) starting date of actual run, i.e. restart date.
    CALL get_date_components(initial_date, year, month, day)
    WRITE(cstart, '(A11,I4.4, 2(A1,I2.2),A9 )') 'days since ', year, &
          '-',month, '-',day, ' 00:00:00'   
    CALL nf_check(nf_put_att_text(ncid, time_id, 'units',       30, cstart))
    CALL nf_check(nf_put_att_text(ncid, time_id, 'calendar',    19, 'proleptic_gregorian'))
    CALL nf_check(nf_put_att_text(ncid, time_id, 'axis',         1, 'T'))

!! #ifdef COUP_OAS
    CALL nf_check(nf_put_att_text(ncid, hdarea_id, 'long_name', 11, 'HDgrid area'))
    CALL nf_check(nf_put_att_text(ncid, hdarea_id, 'standard_name', 9, 'cell_area'))
    CALL nf_check(nf_put_att_text(ncid, hdarea_id, 'units',      2, 'm2'))
!! #endif

    !Standard settings
    ctitle = "Simulated HD model discharge"

    SELECT CASE (TRIM(varname))
       CASE ('friv')
          CALL nf_check(nf_put_att_text(ncid, friv_id, 'long_name',     15, 'river discharge'))
          CALL nf_check(nf_put_att_text(ncid, friv_id, 'standard_name', 39, 'water_volume_transport_in_river_channel'))
          CALL nf_check(nf_put_att_text(ncid, friv_id, 'units',          6, 'm3 s-1'))
       CASE ('disch')
          CALL nf_check(nf_put_att_text(ncid, friv_id, 'long_name',     15, 'discharge to ocean'))
          CALL nf_check(nf_put_att_text(ncid, friv_id, 'standard_name', 45, 'water_volume_transport_into_ocean_from_rivers'))
          CALL nf_check(nf_put_att_text(ncid, friv_id, 'units',          6, 'm s-1'))
       CASE ('friv_bc')
          CALL nf_check(nf_put_att_text(ncid, friv_id, 'long_name',     30, 'bias corrected river discharge'))
          CALL nf_check(nf_put_att_text(ncid, friv_id, 'standard_name', 39, 'water_volume_transport_in_river_channel'))
          CALL nf_check(nf_put_att_text(ncid, friv_id, 'units',          6, 'm3 s-1'))
          ctitle = "Bias corrected HD model discharge"
    END SELECT

    CALL DATE_AND_TIME(cdate)          ! returns YYYYMMDD

    CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'Conventions',    6, 'CF-1.6'))
    CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'source',  42, 'HD model Vs. 5, doi:10.5281/zenodo.4893099'))
    CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'title',   len(TRIM(ctitle)), ctitle))
    CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'licence',        9, 'CC-BY 4.0'))
    CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'crs',           16, 'WGS84, EPSG:4326'))
    CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'creation_date', 10, cdate(7:8)//'.'//cdate(5:6)//'.'//cdate(1:4) ))

    ! Institution and person related attributes
    CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'institution', len(TRIM(hd_inst)), &
         TRIM(hd_inst) ))
    CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'institution-ID', len(TRIM(hd_instid)), &
         TRIM(hd_instid) ))
    CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'contact', len(TRIM(hd_cont)), TRIM(hd_cont) ))
    CALL nf_check(nf_put_att_text(ncid, NF_GLOBAL, 'creator', len(TRIM(hd_user)), TRIM(hd_user) ))

    ! leave define mode
    CALL nf_check(nf_enddef(ncid))

    IF (PRESENT(klon)) THEN
       WRITE (message_text,*)                                    &
         'WARNING: Lon/Lat calculation not adapted for output on non-HD grid'
       CALL message('hd_highres_open', message_text)
       ALLOCATE (lon(klon))
       ALLOCATE (lat(klat))
       DO jb = 1, klat
         lat(jb) = grid_hd%origin_lat-REAL(jb,dp)*grid_hd%resolution+0.5_dp*grid_hd%resolution
       ENDDO
       DO jl = 1, klon
         lon(jl) = REAL(jl,dp)*grid_hd%resolution+grid_hd%origin_lon-0.5_dp*grid_hd%resolution
         IF (lon(jl) >= 180.0_dp) lon(jl) = lon(jl)-360.0_dp
       ENDDO

    ELSE
       ALLOCATE (lon(grid_hd%nlon))
       ALLOCATE (lat(grid_hd%nlat))
       DO jb = 1, grid_hd%nlat
         lat(jb) = grid_hd%origin_lat-REAL(jb,dp)*grid_hd%resolution+0.5_dp*grid_hd%resolution
       ENDDO
       DO jl = 1, grid_hd%nlon
         lon(jl) = REAL(jl,dp)*grid_hd%resolution+grid_hd%origin_lon-0.5_dp*grid_hd%resolution
         IF (lon(jl) >= 180.0_dp) lon(jl) = lon(jl)-360.0_dp
       ENDDO
    ENDIF

    CALL nf_check(nf_put_var_double(ncid, lat_id, lat(:) ))
    CALL nf_check(nf_put_var_double(ncid, lon_id, lon(:) ))

    DEALLOCATE (lon, lat)

  END SUBROUTINE hd_highres_open

  SUBROUTINE hd_highres_close(ncid)

    INTEGER, INTENT(in) :: ncid

    CALL nf_check(nf_close(ncid))

  END SUBROUTINE hd_highres_close

!===========================================================================
  SUBROUTINE hd_highres_write(ncid, friv, icount, klon, klat, hd_area)

    INTEGER,  INTENT(in)    :: ncid
    REAL(dp), INTENT(in)    :: friv(:,:)
    INTEGER,  INTENT(inout) :: icount
    INTEGER,  INTENT(in)    :: klon, klat
    REAL(dp), OPTIONAL, INTENT(in) :: hd_area(:)
    INTEGER :: jl
    REAL(dp), ALLOCATABLE :: hdarea(:,:)

    ! starts and counts for array sections of record variables

    INTEGER ::  tstart(1), tcount(1)
    INTEGER ::  fstart(3), fcount(3)

    INTEGER  :: day, second, day1, sec1
    REAL(dp) :: days_since_start(1)

    IF (PRESENT(hd_area)) THEN
      ALLOCATE (hdarea(klon,klat)) ; hdarea(:,:)   = 0.0_dp
      do jl = 1, klon
        hdarea(jl,:) = hd_area(:)
      enddo
    ENDIF

    icount = icount+1

    ! Current date includes days since start   
!!    CALL get_date_components(current_date, year, month, day, hour, minute, second)
    ! Calculate yyyymmdd
!!    yyyymmdd = ABS(year)*10000+month*100+day              &
!!             + (hour*3600+minute*60+second)/86400._dp
!!    IF (year < 0) yyyymmdd = -yyyymmdd

    ! Calculate days since start 
    CALL TC_get(initial_date, day1, sec1)
    CALL TC_get(current_date, day, second)
    days_since_start(1) = day - day1 + (second-sec1)/86400._dp

    tstart(1) = icount
    tcount(1) = 1
    CALL nf_check(nf_put_vara_double(ncid, time_id, tstart, tcount, days_since_start))

    IF (PRESENT(hd_area)) THEN
      fstart(1:2) = (/ 1, 1 /)
      fcount(1:2) = (/ klon, klat /)
      CALL nf_check(nf_put_vara_double(ncid, hdarea_id, fstart, fcount, hdarea))
    ENDIF

    fstart(1:3) = (/ 1, 1, icount /)
    fcount(1:3) = (/ klon, klat, 1 /)
    CALL nf_check(nf_put_vara_double(ncid, friv_id, fstart, fcount, friv))

    IF (PRESENT(hd_area)) THEN
      DEALLOCATE  (hdarea)
    ENDIF

  END SUBROUTINE hd_highres_write
  
END MODULE mo_hd_highres_io
