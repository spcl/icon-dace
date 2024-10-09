! mo_flow_inout.f90 - Utilities for Netcdf IN and OUTput
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

MODULE mo_flow_inout

  !
  ! Authors:
  !
  ! S. Hagemann, HZG-IfK, June 2020- , original source

  use netcdf

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: open_inflow, read_inflow, open_outflow, write_inflow_on_ocean, &
            close_inflow, close_outflow, read_ute_inflow, set_variable, &
            open_output, define_output_variable, write_output_variable, close_output, &
            set_hereon_attributes

  INTEGER, PUBLIC      :: ncid_in     ! ID of Netcdf input file
  INTEGER, PUBLIC      :: ncid_out    ! ID of Netcdf output file
  INTEGER, PUBLIC      :: time_id_in  ! Time ID of input (0 if none) 
  INTEGER, PUBLIC      :: time_id_out ! Time ID of output 

  INTEGER, PUBLIC      :: nvar_in = 1     ! No. of Input Variables
  INTEGER, PUBLIC      :: nvar_out = 1    ! No. of Output Variables - initially set to nvar_in
                                          ! They can be overwritten in routine set_variable

  INTEGER, PUBLIC      :: nvar        ! No. of Variables
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: idvar_in  ! Variables IDs input
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: idvar_out ! Variables IDs output
!
! *** Variables read from input and initially set for output - They can be overwritten in routine set_variable
  CHARACTER (len=40), PUBLIC, DIMENSION(:), ALLOCATABLE :: cvar  ! Variable names
  CHARACTER (len=40), PUBLIC, DIMENSION(:), ALLOCATABLE :: cunit ! Units of Variables
  CHARACTER (len=60), PUBLIC, DIMENSION(:), ALLOCATABLE :: clong ! Long names

  INTEGER :: ierr
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)


CONTAINS

! ******************************************************************
  SUBROUTINE open_inflow(dninp)
! ******************************************************************
!
!   *** Open inflow data and check time

    CHARACTER (LEN=*), INTENT(in)  :: dninp       ! Input file name

    INTEGER :: ierr
    INTEGER, DIMENSION(:), ALLOCATABLE :: idvar
    INTEGER :: ndum, i, ndel, j, idfix
    CHARACTER (len=40), DIMENSION(5) :: cvar2
    CHARACTER (len=40), DIMENSION(:), ALLOCATABLE :: cvarname

    INTEGER, PARAMETER :: nfix = 13 
    CHARACTER (len=40), DIMENSION(nfix) :: cvfix
    DATA cvfix / 'time', 'lat', 'lon', 'northing', 'easting', 'plev', &
                 'lat_2', 'lon_2', 'lat_2_bnds', 'lon_2_bnds', 'depth', 'hdarea', &
                 'time_bnds'  /
!
!   ******** Open infow (e.g. discharge) input file
    ierr = nf90_open(dninp, NF90_NOWRITE, ncid_in)      ! input
    IF (ierr.NE.0) THEN
       WRITE(*,*) 'Error opening ', TRIM(dninp)
       STOP ' ---> Program is terminated !!!!!'
    ENDIF
    ierr = nf90_inquire(ncid_in, nvariables = ndum)

    WRITE(*,*) 'Input file = ', TRIM(dninp)
    WRITE(*,*) 'No. of Variables in Input file = ', ndum
    ALLOCATE(idvar(ndum))
    ALLOCATE(cvarname(ndum))
    ierr = nf90_inq_varids(ncid_in, ndum, idvar)
    DO i=1, ndum
      ierr = nf90_inquire_variable(ncid_in, idvar(i), name=cvarname(i))
      WRITE(*,*) idvar(i), cvarname(i)
    ENDDO
!
!   *** Exclude non time dependent variables
    DO j=1, nfix
      ierr = nf90_inq_varid(ncid_in,cvfix(j),idfix)
!     *** Is time included in variables (J=1)
      IF (j.EQ.1) THEN
        IF (ierr.NE.0) THEN
          time_id_in = 0
          WRITE(*,*) 'No time variable found in input data'
        ELSE
          time_id_in = idfix
        ENDIF
      ENDIF
      IF (ierr.EQ.0) THEN
        ndel=0
        DO i=1, ndum
          IF (idvar(i).EQ.idfix) THEN
            WRITE(*,*) TRIM(cvarname(i)), ' variable will not be regarded in inflow'
            ndel = 1          !
          ELSE IF (ndel.EQ.1) THEN
            idvar(i-1) = idvar(i)
            cvarname(i-1) = cvarname(i)
          ENDIF
        ENDDO
        ndum = ndum - ndel
      ENDIF
    ENDDO

    nvar = ndum
    ALLOCATE(idvar_in(nvar))
    if (ALLOCATED(cvar)) DEALLOCATE(cvar)
    ALLOCATE(cvar(nvar))
    DO i=1,nvar
      idvar_in(i) = idvar(i)
      cvar(i) = cvarname(i)
      IF (TRIM(cvar(i)).EQ.'var0') cvar(i)='discharge'
      IF (TRIM(cvar(i)).EQ.'friv') cvar(i)='discharge'
    ENDDO

    DEALLOCATE(idvar)
    DEALLOCATE(cvarname)

  END SUBROUTINE open_inflow
!
! ******************************************************************
  SUBROUTINE read_inflow(nlon, nlat, istep, ivar, finflow, yyyymmdd, ierr)
! ******************************************************************
    INTEGER, INTENT(in) :: nlon 
    INTEGER, INTENT(in) :: nlat 
    INTEGER, INTENT(in) :: istep  ! Time step number
    INTEGER, INTENT(in) :: ivar   ! Variable number
    REAL, DIMENSION(nlon,nlat), INTENT(out) :: finflow 
    DOUBLE PRECISION, INTENT(out) :: yyyymmdd     ! date + time in unit [day]
    INTEGER, INTENT(out) :: ierr 

    IF (time_id_in.EQ.0) THEN
      ierr = nf90_get_var(ncid_in, idvar_in(ivar), finflow)
      yyyymmdd = 0.
    ELSE
      ierr = nf90_get_var(ncid_in, time_id_in, yyyymmdd, start = (/ istep /) )   
      IF (yyyymmdd .LT. 10000.) yyyymmdd = yyyymmdd * 10000. 
      ierr = nf90_get_var(ncid_in, idvar_in(ivar), finflow, &
           start = (/ 1, 1, istep /), count = (/ nlon, nlat, 1 /))
    ENDIF

  END SUBROUTINE read_inflow
!
! ******************************************************************
  SUBROUTINE read_ute_inflow(nlon, istep, ivar, xmiss, finflow, yyyymmdd, ierr)
! ******************************************************************
    INTEGER, INTENT(in) :: nlon 
    INTEGER, INTENT(in) :: istep  ! Time step number
    INTEGER, INTENT(in) :: ivar   ! Variable number
    DOUBLE PRECISION, INTENT(in)    :: xmiss 
    REAL, DIMENSION(nlon), INTENT(out) :: finflow 
    DOUBLE PRECISION, INTENT(out) :: yyyymmdd     ! date + time in unit [day]
    INTEGER, INTENT(out) :: ierr 
!
    DOUBLE PRECISION, DIMENSION(nlon, 28125)  :: flow_ute
!
    ierr = nf90_get_var(ncid_in, idvar_in(ivar), flow_ute)
    WHERE (flow_ute(:,istep).NE.flow_ute(:,istep))
      flow_ute(:,istep) = 0.
    END WHERE
    finflow(:) = REAL(flow_ute(:, istep))
    yyyymmdd = 0.

  END SUBROUTINE read_ute_inflow
!
! ******************************************************************
  SUBROUTINE set_variable(id_var, nset)
! ******************************************************************
!
!   *** Set variable names cvar(:) for output - Currently fixed for nvar_out = 1
!   *** Useful if variable names were not read from NetCDF file
!
    INTEGER, DIMENSION(:), INTENT(in)  :: id_var   ! Own variable No.
    INTEGER, INTENT(in), OPTIONAL      :: nset     ! No. of variables to be set 
    INTEGER :: i
!
    IF (PRESENT(nset)) THEN
      nvar_out = nset
    ELSE
      nvar_out = size(id_var)
    ENDIF
!
    IF (ALLOCATED(cvar)) DEALLOCATE(cvar)
    IF (ALLOCATED(cunit)) DEALLOCATE(cunit)
    IF (ALLOCATED(clong)) DEALLOCATE(clong)
    ALLOCATE(cvar(nvar_out))
    ALLOCATE(cunit(nvar_out))
    ALLOCATE(clong(nvar_out))
    DO i = 1, nvar_out
      SELECT CASE (id_var(i))
        CASE(70) ; cvar(i) = 'Discharge' ; cunit(i) = '[m³/s]' ; clong(i) = "water_volume_transport_in_river_channel"
        CASE(71) ; cvar(i) = 'friv_nitrogen' ; cunit(i) = '[t/d]' ; clong(i) = "nitrogen_mass_transport_in_river_channel"
        CASE(72) ; cvar(i) = 'friv_phosphorus' ; cunit(i) = '[t/d]' ; clong(i) = "phosphorus_mass_transport_in_river_channel"
        CASE(73) ; cvar(i) = 'friv_silicate' ; cunit(i) = '[t/d]' ; clong(i) = "silicate_mass_transport_in_river_channel"
        CASE(74) ; cvar(i) = 'friv_doc' ; cunit(i) = '[t/d]' ; clong(i) = "DOC_mass_transport_in_river_channel"
        CASE DEFAULT ; WRITE(*,*) 'SET_VARIABLE: Variable ID ', id_var(i), ' not defined --> ABBRUCH' ; STOP
      END SELECT
    ENDDO

  END SUBROUTINE set_variable
!
! ******************************************************************
  SUBROUTINE open_outflow(dnout, mask_target, iconv)
! ******************************************************************
!
!     *** Open discharge on ocean grid output file and write coordinates
    use mo_grid,         ONLY: model
	  
    CHARACTER (LEN=*), INTENT(in)  :: dnout       ! Output file name
    TYPE(model), INTENT(in)        :: mask_target ! Target info, e.g.
    INTEGER, INTENT(in) :: iconv ! Method of conversion: 0 = no, 1=Daily data (def.),
                                 ! 2=monthly clim., 3=5 bgc sequence
                                 ! 4=General sequence of bgc variables and timesteps
!
!     *** NETCDF variables
    INTEGER :: ierr, dimids(2), lon_dim, lat_dim, time_dim
    INTEGER :: dim1(1), dim2(2), dim3(3), varid
    INTEGER :: i
    REAL :: xmiss = -9999.
!
    ierr = nf90_create(dnout, NF90_HDF5, ncid_out)  ! NF90_CLOBBER was used 
    ierr = nf90_def_dim(ncid_out, 'x', mask_target%nlon, lon_dim)
    ierr = nf90_def_dim(ncid_out, 'y', mask_target%nlat, lat_dim)
    dimids(1) = lon_dim
    dimids(2) = lat_dim
!   *** Time
    ierr = nf90_def_dim(ncid_out, 'time', NF90_UNLIMITED, time_dim)
    dim1(1) = time_dim
    ierr = nf90_def_var(ncid_out, 'time', NF90_DOUBLE, dim1, time_id_out)
    ierr = nf90_def_var_deflate(ncid_out, time_id_out, shuffle = 0, &
                     deflate = 1, deflate_level = 2)
    ierr = nf90_put_att(ncid_out,time_id_out,'units','day as %Y%m%d.%f')
    ierr = nf90_put_att(ncid_out,time_id_out,'standard_name','time')
    ierr = nf90_put_att(ncid_out,time_id_out,'calendar','proleptic_gregorian')
    ierr = nf90_put_att(ncid_out,time_id_out,'axis','T')

    ierr = nf90_def_var(ncid_out, 'lon', NF90_DOUBLE, dimids, varid)
    ierr = nf90_def_var_deflate(ncid_out, varid, shuffle = 0, &
                     deflate = 1, deflate_level = 2)
    ierr = nf90_put_att(ncid_out,varid,'units','degrees_E')
    ierr = nf90_put_att(ncid_out,varid,'long_name','longitude')
    ierr = nf90_put_att(ncid_out,varid,'standard_name','longitude')
    ierr = nf90_enddef(ncid_out)
    ierr = nf90_put_var(ncid_out, varid, mask_target%xlon)

    ierr = nf90_redef(ncid_out)
    ierr = nf90_def_var(ncid_out, 'lat', NF90_DOUBLE, dimids, varid)
    ierr = nf90_def_var_deflate(ncid_out, varid, shuffle = 0, &
                     deflate = 1, deflate_level = 2)
    ierr = nf90_put_att(ncid_out,varid,'units','degrees_N')
    ierr = nf90_put_att(ncid_out,varid,'long_name','latitude')
    ierr = nf90_put_att(ncid_out,varid,'standard_name','latitude')
    ierr = nf90_enddef(ncid_out)
    ierr = nf90_put_var(ncid_out, varid, mask_target%xlat)

    ierr = nf90_redef(ncid_out)

    ALLOCATE(idvar_out(nvar))
    IF (iconv.LE.2) THEN
      dim3(1:3) = (/ lon_dim, lat_dim, time_dim /)
      DO i=1, nvar
!!        IF (nvar.EQ.1 .OR. TRIM(cvar(i)).EQ.'discharge') THEN
        IF (TRIM(cvar(i)).EQ.'discharge') THEN
          ierr = nf90_def_var(ncid_out, 'discharge_on_ocean', NF90_float, dim3, idvar_out(i))
          ierr = nf90_def_var_deflate(ncid_out, idvar_out(i), shuffle = 0, &
                     deflate = 1, deflate_level = 2)
          ierr = nf90_put_att(ncid_out,idvar_out(i),'units','[m3 s-1]')
          ierr = nf90_put_att(ncid_out,idvar_out(i),'code', 219)
        ELSE
          ierr = nf90_def_var(ncid_out, TRIM(cvar(i)), NF90_float, dim3, idvar_out(i))
          ierr = nf90_def_var_deflate(ncid_out, idvar_out(i), shuffle = 0, &
                     deflate = 1, deflate_level = 2)
          ierr = nf90_put_att(ncid_out,idvar_out(i),'units','[t/d]')
          ierr = nf90_put_att(ncid_out,idvar_out(i),'code', 900+i)
        ENDIF
        ierr = nf90_put_att(ncid_out,idvar_out(i),'long_name', TRIM(cvar(i)) // ' inflow')
        ierr = nf90_put_att(ncid_out,idvar_out(i),'standard_name',TRIM(cvar(i))  &
               // "_transport_into_ocean_from_rivers")
        ierr = nf90_put_att(ncid_out,idvar_out(i),'coordinates','lon lat')
        ierr = nf90_put_att(ncid_out,idvar_out(i),'missing_value', xmiss)
        WRITE(*,*) ' Defining output for ', TRIM(cvar(i)) // '_inflow', &
                ' time-id in:', time_id_in
      ENDDO
      ierr = nf90_enddef(ncid_out)
!
    ELSE IF (iconv.GE.3) THEN
      IF (time_id_in.EQ.0) time_id_out = 0 
      IF (time_id_out.EQ.0) THEN
        dim2(1:2) = (/ lon_dim, lat_dim /)
      ELSE
        dim3(1:3) = (/ lon_dim, lat_dim, time_dim /)
      ENDIF
      DO i=1, nvar
        IF (time_id_out.EQ.0) THEN
          ierr = nf90_def_var(ncid_out, TRIM(cvar(i)), NF90_float, dim2, idvar_out(i))
        ELSE
          ierr = nf90_def_var(ncid_out, TRIM(cvar(i)), NF90_float, dim3, idvar_out(i))
        ENDIF
        ierr = nf90_put_att(ncid_out,idvar_out(i),'code', 900+i)
        ierr = nf90_put_att(ncid_out,idvar_out(i),'units','[Mg ' // TRIM(cvar(i)) // 'yr-1]')
        ierr = nf90_put_att(ncid_out,idvar_out(i),'long_name', TRIM(cvar(i)) // ' inflow')
        ierr = nf90_put_att(ncid_out,idvar_out(i),'standard_name',TRIM(cvar(i))  &
               // "_transport_into_${COC}_ocean_from_rivers")
        ierr = nf90_put_att(ncid_out,idvar_out(i),'coordinates','lon lat')
        WRITE(*,*) ' Defining output for ', TRIM(cvar(i)) // '_inflow', &
                ' time-id in:', time_id_in
      ENDDO
      ierr = nf90_enddef(ncid_out)
    ENDIF

  END SUBROUTINE open_outflow

! ******************************************************************
  SUBROUTINE write_inflow_on_ocean(nxocean, nyocean, icount, ivar, &
             friv_on_ocean, yyyymmdd)
! ******************************************************************

    INTEGER,  INTENT(in) :: nxocean
    INTEGER,  INTENT(in) :: nyocean
    INTEGER,  INTENT(in) :: icount
    INTEGER,  INTENT(in) :: ivar
    REAL,     INTENT(in) :: friv_on_ocean(:,:)
    DOUBLE PRECISION, INTENT(in) :: yyyymmdd     ! date + time in unit [day]

    ! starts and counts for array sections of record variables

    INTEGER ::  tstart(1), tcount(1)
    INTEGER ::  fstart(3), fcount(3)
    INTEGER :: ierr
    INTEGER :: isav = 0
    SAVE isav

    IF (time_id_out.EQ.0) THEN
      ierr = nf90_put_var(ncid_out, idvar_out(ivar), friv_on_ocean)
    ELSE
      IF (isav.NE.icount) THEN
        tstart(1) = icount
        ierr = nf90_put_var(ncid_out, time_id_out, yyyymmdd, start=tstart)
        IF (ierr.NE.0) WRITE(*,*) icount, '. error writing time: ', ierr
      ENDIF

      fstart(1:3) = (/ 1, 1, icount /)
      fcount(1:3) = (/ nxocean, nyocean, 1 /)

      ierr = nf90_put_var(ncid_out, idvar_out(ivar), friv_on_ocean, start=fstart, count=fcount)
      isav = icount
    ENDIF

  END SUBROUTINE write_inflow_on_ocean

  SUBROUTINE close_inflow

    INTEGER :: ierr
    ierr = nf90_close(ncid_in)
    DEALLOCATE(idvar_in)

  END SUBROUTINE close_inflow

  SUBROUTINE close_outflow

    INTEGER :: ierr
    ierr = nf90_close(ncid_out)
    DEALLOCATE(idvar_out)

  END SUBROUTINE close_outflow
!
! ******************************************************************
  SUBROUTINE open_output(dnout, mask_target, ncid_output, dimids, idim)
! ******************************************************************
!
!   *** Open output file for general data with lat/lon info but without time dependence
!
    use mo_grid,         ONLY: model
	  
    CHARACTER (LEN=*), INTENT(in)  :: dnout       ! Output file name
    TYPE(model), INTENT(in)        :: mask_target ! Target info, e.g.
    INTEGER, INTENT(out) :: ncid_output           ! Output file ID
    INTEGER, INTENT(out) :: dimids(2)             ! Dimension IDs of 2D array     
    INTEGER, INTENT(in), OPTIONAL  :: idim    ! Number of coordinate dimensions: 1 or 2 (default)
!
!     *** NETCDF variables
    INTEGER :: lon_dim, lat_dim, varid
    LOGICAL :: log2
!
    log2 = .True.
    IF (PRESENT(idim)) THEN
      IF (idim.EQ.1) log2 = .FALSE.
    ENDIF

    ierr = nf90_create(dnout, NF90_HDF5, ncid_output) 
    ierr = nf90_def_dim(ncid_output, 'lon', mask_target%nlon, lon_dim)
    ierr = nf90_def_dim(ncid_output, 'lat', mask_target%nlat, lat_dim)
    dimids(1) = lon_dim
    dimids(2) = lat_dim

    IF (log2) THEN
      ierr = nf90_def_var(ncid_output, 'lon', NF90_DOUBLE, dimids, varid)
    ELSE
      ierr = nf90_def_var(ncid_output, 'lon', NF90_DOUBLE, lon_dim, varid)
    ENDIF
    ierr = nf90_def_var_deflate(ncid_output, varid, shuffle = 0, &
                     deflate = 1, deflate_level = 2)
    ierr = nf90_put_att(ncid_output,varid,'units','degrees_E')
    ierr = nf90_put_att(ncid_output,varid,'long_name','longitude')
    ierr = nf90_put_att(ncid_output,varid,'standard_name','longitude')
    ierr = nf90_put_att(ncid_output,varid,'axis','X')
    ierr = nf90_enddef(ncid_output)
    IF (log2) THEN
      ierr = nf90_put_var(ncid_output, varid, mask_target%xlon)
    ELSE
      ierr = nf90_put_var(ncid_output, varid, mask_target%xlon(:,1))
    ENDIF

    ierr = nf90_redef(ncid_output)
    IF (log2) THEN
      ierr = nf90_def_var(ncid_output, 'lat', NF90_DOUBLE, dimids, varid)
    ELSE
      ierr = nf90_def_var(ncid_output, 'lat', NF90_DOUBLE, lat_dim, varid)
    ENDIF
    ierr = nf90_def_var_deflate(ncid_output, varid, shuffle = 0, &
                     deflate = 1, deflate_level = 2)
    ierr = nf90_put_att(ncid_output,varid,'units','degrees_N')
    ierr = nf90_put_att(ncid_output,varid,'long_name','latitude')
    ierr = nf90_put_att(ncid_output,varid,'standard_name','latitude')
    ierr = nf90_put_att(ncid_output,varid,'axis','Y')
    ierr = nf90_enddef(ncid_output)
    IF (log2) THEN
      ierr = nf90_put_var(ncid_output, varid, mask_target%xlat)
    ELSE
      ierr = nf90_put_var(ncid_output, varid, mask_target%xlat(1,:))
    ENDIF
!
  END SUBROUTINE open_output

! ******************************************************************
  SUBROUTINE define_output_variable(ncid, cvar, ctype, dimids, varid, xmiss, cunit, clong, icode)
! ******************************************************************
!
    INTEGER, INTENT(in)           :: ncid   ! Output file ID
    CHARACTER (len=*), INTENT(in) :: cvar          ! Variable name
    CHARACTER (len=*), INTENT(in) :: ctype         ! Variable type: 'int', 'float, 'dp'
    INTEGER, INTENT(in)           :: dimids(2)     ! Dimensional IDs
    INTEGER, INTENT(out)          :: varid
    REAL(dp), INTENT(in), OPTIONAL  :: xmiss         ! Missing Value
    CHARACTER (LEN=*), INTENT(in), OPTIONAL  :: cunit       !  Unit
    CHARACTER (LEN=*), INTENT(in), OPTIONAL  :: clong       !  Long name
    INTEGER, INTENT(in), OPTIONAL :: icode         ! Code no. 

    ierr = nf90_redef(ncid)
    IF (ctype.EQ.'int') ierr = nf90_def_var(ncid, cvar, NF90_int, dimids, varid)
    IF (ctype.EQ.'float') ierr = nf90_def_var(ncid, cvar, NF90_float, dimids, varid)
    IF (ctype.EQ.'dp') ierr = nf90_def_var(ncid, cvar, NF90_DOUBLE, dimids, varid)
    ierr = nf90_def_var_deflate(ncid, varid, shuffle = 0, &
                     deflate = 1, deflate_level = 2)
    IF (PRESENT(xmiss)) THEN
      IF (ctype.EQ.'int') ierr = nf90_put_att(ncid,varid,'missing_value', NINT(xmiss))
      IF (ctype.EQ.'float') ierr = nf90_put_att(ncid,varid,'missing_value', REAL(xmiss))
      IF (ctype.EQ.'dp') ierr = nf90_put_att(ncid,varid,'missing_value', xmiss)
    ENDIF
    IF (PRESENT(cunit)) THEN
      ierr = nf90_put_att(ncid,varid,'units', cunit)
    ELSE
      ierr = nf90_put_att(ncid,varid,'units','[-]')
    ENDIF
    IF (PRESENT(clong)) ierr = nf90_put_att(ncid,varid,'long_name',cvar )
    IF (PRESENT(icode)) THEN
      ierr = nf90_put_att(ncid,varid,'code', icode)
    ENDIF
!!    ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
    ierr = nf90_enddef(ncid)

  END SUBROUTINE define_output_variable

! ******************************************************************
  SUBROUTINE write_output_variable(ncid, varid, fdat)
! ******************************************************************
!
    INTEGER, INTENT(in)              :: ncid   ! Output file ID
    INTEGER, INTENT(in)              :: varid  ! Variable ID
    REAL, DIMENSION(:,:), INTENT(in) :: fdat   ! 2D Data array 

    ierr = nf90_put_var(ncid, varid, fdat )

  END SUBROUTINE write_output_variable

! ******************************************************************
  SUBROUTINE close_output(ncid)
! ******************************************************************

    INTEGER, INTENT(in)              :: ncid   ! Output file ID

    ierr = nf90_close(ncid)

  END SUBROUTINE close_output

! ******************************************************************
  SUBROUTINE set_hereon_attributes(ctitle, ncid_output, cexpid, csource, cforcing, cref, corg)
! ******************************************************************
!
!!    INTEGER, INTENT(in) :: ncid_out           ! Output file ID
    CHARACTER (LEN=*), INTENT(in)  :: ctitle              ! Title of output data
    INTEGER, INTENT(in), OPTIONAL  :: ncid_output           ! Output file ID
    CHARACTER (LEN=*), INTENT(in), OPTIONAL  :: cexpid     !  Experiment ID
    CHARACTER (LEN=*), INTENT(in), OPTIONAL  :: csource    !  Source of model data: Model Vs. and DOI
    CHARACTER (LEN=*), INTENT(in), OPTIONAL  :: cforcing   !  Forcing used
    CHARACTER (LEN=*), INTENT(in), OPTIONAL  :: cref       !  References
    CHARACTER (LEN=*), INTENT(in), OPTIONAL  :: corg       !  Creator (formerly originator)

    CHARACTER(8) :: cdate     ! Date in reality = date of creating output file
    INTEGER :: ierr
    INTEGER :: ncid

    ! Namelist with user attributes
    CHARACTER(200) :: hd_user, hd_cont, hd_inst, hd_instid
!!    NAMELIST /HDUSER_CTL/ &
!!       hd_user, hd_cont, hd_inst, hd_instid
    hd_user = "Stefan Hagemann"
    hd_cont = "stefan.hagemann@hereon.de, https://coastmod.hereon.de"
    hd_inst = "Helmholtz-Zentrum Hereon, Institute of Coastal Systems, Germany"
    hd_instid = "ROR: 03qjp1d79"

    IF (PRESENT(ncid_output)) THEN
      ncid = ncid_output
    ELSE
      ncid = ncid_out
    ENDIF

    CALL DATE_AND_TIME(cdate)          ! returns YYYYMMDD

    ierr = nf90_redef(ncid) 
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions',   'CF-1.6')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'title',         TRIM(ctitle))
    IF (PRESENT(cexpid)) THEN
      ierr = nf90_put_att(ncid, NF90_GLOBAL, 'experiment_id', TRIM(cexpid) )
    ENDIF
    IF (PRESENT(csource)) THEN
      ierr = nf90_put_att(ncid, NF90_GLOBAL, 'source', TRIM(csource) )
    ENDIF
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'licence',       'CC-BY 4.0')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'crs',           'WGS84, EPSG:4326')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', cdate(7:8)//'.'//cdate(5:6)//'.'//cdate(1:4) )
    IF (PRESENT(cforcing)) THEN
      ierr = nf90_put_att(ncid, NF90_GLOBAL, 'forcing', TRIM(cforcing) )
    ENDIF
    IF (PRESENT(cref)) THEN
      ierr = nf90_put_att(ncid, NF90_GLOBAL, 'references', TRIM(cref) )
    ENDIF

    ! Institution and person related attributes
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'institution',   TRIM(hd_inst) )
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'institution-ID', TRIM(hd_instid) )
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'contact',       TRIM(hd_cont) )
    IF (PRESENT(corg)) THEN
      ierr = nf90_put_att(ncid, NF90_GLOBAL, 'creator',  TRIM(corg) )
    ELSE
      ierr = nf90_put_att(ncid, NF90_GLOBAL, 'creator',  TRIM(hd_user) )
    ENDIF

    ! leave define mode
    ierr = nf90_enddef(ncid)

  END SUBROUTINE set_hereon_attributes
 
END MODULE mo_flow_inout

