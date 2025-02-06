! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Chemical heating rate
!
! This module allows to compute chemical heating by interpolating from zonal mean climatology (12-month)
!
!  Modification by Guidi Zhou, MPI-M, 2017-02-28:
!  - added support for using chemical heating only above a certain level
!  - added output of efficiency factor for standard shortwave radiation
!  Modification by Guidi Zhou, MPI-M (2017-03-06)
!  - added the ability to compute chemical heating only above a certain altitude for performance

MODULE mo_upatmo_phy_chemheat

  USE mo_kind,                   ONLY: wp, i8
  USE mo_parallel_config,        ONLY: p_test_run
  USE mo_mpi,                    ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                                  p_comm_work_test, p_comm_work
  USE mo_util_string,            ONLY: tolower, t_keyword_list, associate_keyword, with_keywords
  USE mo_netcdf_errhandler,      ONLY: nf
  USE mo_netcdf
  USE mo_master_config,          ONLY: getModelBaseDir
  USE mo_impl_constants,         ONLY: max_char_length, SUCCESS
  USE mo_math_constants,         ONLY: deg2rad
  USE mo_exception,              ONLY: finish, message
  USE mtime,                     ONLY: datetime, julianday, newDatetime, newJulianday, getJulianDayFromDatetime, &
    &                                  deallocateDatetime, deallocateJulianday

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: chem_heat_check      ! subroutine to check existence and consistency of chemical heating input file
  PUBLIC :: chem_heat_init       ! subroutine to initialize/read chemical heating climatology data from input file
  PUBLIC :: chem_heat            ! subroutine to compute (interpolate) chemical heating at current time step
  PUBLIC :: chem_heat_clean      ! subroutine to clean up

  ! private
  LOGICAL :: lchemheat_checked     = .FALSE.
  LOGICAL :: lchemheat_initialized = .FALSE.
  LOGICAL :: lchemheat_finalized   = .FALSE.

  ! name of chemical heating climatology data input file
  CHARACTER(len=*), PARAMETER :: filename_template = "<path>chemheat.nc" ! name consisting of template
  CHARACTER(len=max_char_length) :: filename                             ! actual name

  INTEGER, PARAMETER :: ntime = 12                ! number of timesteps of chemical heating climatology read from input file
  INTEGER :: nlev                                 ! number of levels    of chemical heating climatology read from input file
  INTEGER :: nlat                                 ! number of latitudes of chemical heating climatology read from input file
  INTEGER :: varid_ch, varid_lev, varid_lat       ! variable IDs        of chemical heating climatology read from input file
  INTEGER :: varid_time                           ! variable ID for time
  REAL(wp) :: minlev, maxlev, minlat, maxlat      ! min/max of level and latitude
  LOGICAL :: levinc, latinc                       ! is level/latitude in ascending order?

  REAL(wp), PARAMETER, PUBLIC :: zeroz = 70000.0_wp  ! altitude below which chemical heating is zero
  REAL(wp), PARAMETER, PUBLIC :: onez  = 80000.0_wp  ! altitude above which full chemical heating is used
  REAL(wp), PARAMETER, PUBLIC :: effrswmin = 0.23_wp ! efficiency factor for standard shortwave radiation (>=200nm) 
                                                     ! when full chemical heating is used

  REAL(wp) :: scl_ch                              ! scale factor to convert chemical heating data to K/s
  REAL(wp) :: scl_lev                             ! scale factor to convert level data to Pa or m, 
                                                  ! depending on type of vertical level

  REAL(wp), ALLOCATABLE :: levclim(:)             ! level    data for chemical heating climatology read from input file
  REAL(wp), ALLOCATABLE :: latclim(:)             ! latitude data for chemical heating climatology read from input file
  REAL(wp), ALLOCATABLE :: timeclim(:)            ! time     data for chemical heating climatology read from input file
  REAL(wp), ALLOCATABLE :: chclim(:, :, :)        ! chemical heating climatology data read from input file

  CHARACTER(len=*), PARAMETER :: acceptable_names_lev(12) = ['lev     ', &
    &                                                        'levs    ', &
    &                                                        'level   ', &
    &                                                        'levels  ', &
    &                                                        'z       ', &
    &                                                        'zlev    ', &
    &                                                        'zlevs   ', &
    &                                                        'zlevel  ', &
    &                                                        'zlevels ', &
    &                                                        'alt     ', &
    &                                                        'altitude', &
    &                                                        'height  '  &
    &                                                       ]

  CHARACTER(len=*), PARAMETER :: acceptable_names_lat(2) = ['lat     ', &
    &                                                       'latitude'  &
    &                                                      ]

  CHARACTER(len=*), PARAMETER :: acceptable_names_ch(8) =  ['dt_cheat         ', &
    &                                                       'dt_chem          ', &
    &                                                       'dt_chemheat      ', &
    &                                                       'tend_ta_cheat    ', &
    &                                                       'tend_ta_chem     ', &
    &                                                       'tend_ta_chemheat ', &
    &                                                       'chem_heat        ', &
    &                                                       'ddt_temp_chemheat'  &
    &                                                      ]

  CHARACTER(len=*), PARAMETER :: modname = "mo_upatmo_phy_chemheat"

  TYPE t_intp_weight
    INTEGER :: prev, next
    REAL(wp) :: factor
  END TYPE t_intp_weight

CONTAINS

  SUBROUTINE chem_heat_check(opt_filename, opt_nlat, opt_nlev, opt_ntime, &
    &                        opt_unitchemheat, opt_unitlat, opt_unitlev, opt_unittime)

    ! IN/OUT
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: opt_filename
    INTEGER,          OPTIONAL, INTENT(OUT) :: opt_nlat
    INTEGER,          OPTIONAL, INTENT(OUT) :: opt_nlev
    INTEGER,          OPTIONAL, INTENT(OUT) :: opt_ntime
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: opt_unitchemheat
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: opt_unitlat
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: opt_unitlev
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: opt_unittime

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: routine = modname // ":chem_heat_check"
    LOGICAL :: file_exists
    CHARACTER(len=max_char_length) :: vn, temp, levname, varname, units, message_text
    INTEGER :: mpi_comm, ncid, nvars, id, ian, dimid_time, dimid_lat, dimid_lev, ntimestep, ndims, dimlen
    INTEGER, ALLOCATABLE :: dimids(:)
    CHARACTER(LEN=NF90_MAX_NAME) :: attname
    INTEGER                    :: natts, iatt

    !---------------------------------------------------------
    
    IF (.NOT. lchemheat_checked) THEN
      
      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF
      
      IF (my_process_is_stdio()) THEN
        ! generate filename
        IF (PRESENT(opt_filename)) THEN
          filename = TRIM(opt_filename)
        ELSE
          filename = generate_filename(TRIM(ADJUSTL(filename_template)), getModelBaseDir())
        ENDIF
        CALL message(routine, 'chemheat_filename=' // TRIM(filename))
        
        ! check file existence
        INQUIRE(FILE=TRIM(filename), EXIST=file_exists)  
        IF (.NOT. file_exists) CALL finish(routine, "input file " // TRIM(filename) // " does not exist")
        
        CALL nf(nf90_open(TRIM(filename), NF90_NOWRITE, ncid), routine)
        
        ! check dimensions - time
        CALL nf(nf90_inq_dimid(ncid, 'time', dimid_time), routine) 
        IF (dimid_time < 0) CALL finish(routine, "input file " // TRIM(filename) // " does not have a time dimension")
        CALL nf(nf90_inquire_dimension(ncid, dimid_time, len = ntimestep), routine)
        IF (ntimestep /= ntime) CALL finish(routine, "input file " // TRIM(filename) // " must have 12 timesteps")
        CALL nf(nf90_inq_varid(ncid, 'time', varid_time), routine)
        
        ! check dimensions - latitude
        CALL nf(nf90_inq_dimid(ncid, 'lat', dimid_lat), routine) 
        IF (dimid_lat < 0) THEN
          CALL nf(nf90_inq_dimid(ncid, 'latitude', dimid_lat), routine) 
          IF (dimid_lat < 0) CALL finish(routine, "input file " // TRIM(filename) // " does not have a latitude dimension")
        END IF
        CALL nf(nf90_inquire_dimension(ncid, dimid_lat, len = nlat), routine)
        WRITE(temp, *)nlat
        CALL message(routine, "input file " // TRIM(filename) // " has " // TRIM(ADJUSTL(temp)) // " latitudes")
        
        ! check dimensions - level
        CALL nf(nf90_inq_dimid(ncid, 'lev', dimid_lev), routine) 
        IF (dimid_lev < 0) THEN
          CALL nf(nf90_inq_dimid(ncid, 'level', dimid_lev), routine) 
          IF (dimid_lev < 0) CALL finish(routine, "input file " // TRIM(filename) // " does not have a level dimension")
        END IF
        CALL nf(nf90_inquire_dimension(ncid, dimid_lev, len = nlev), routine)
        WRITE(temp, *)nlev
        CALL message(routine, "input file " // TRIM(filename) // " has " // TRIM(ADJUSTL(temp)) // " levels")
        
        ! check variables
        CALL nf(nf90_inquire(ncid, nVariables = nvars), routine)
        
        varid_lev = -1
        varid_lat = -1
        varid_ch = -1
        DO id = 1, nvars
          CALL nf(nf90_inquire_variable(ncid, id, name = vn), routine)
          
          IF (varid_lev == -1) THEN
            DO ian = 1, SIZE(acceptable_names_lev)
              IF ( TRIM(tolower(vn)) == TRIM(tolower(acceptable_names_lev(ian))) ) THEN
                varid_lev = id
                levname = vn
                EXIT
              ENDIF
            END DO
          END IF
          
          IF (varid_lat == -1) THEN
            DO ian = 1, SIZE(acceptable_names_lat)
              IF ( TRIM(tolower(vn)) == TRIM(tolower(acceptable_names_lat(ian))) ) THEN
                varid_lat = id
                EXIT
              ENDIF
            END DO
          END IF
          
          IF (varid_ch == -1) THEN
            DO ian = 1, SIZE(acceptable_names_ch)
              IF ( TRIM(tolower(vn)) == TRIM(tolower(acceptable_names_ch(ian))) ) THEN
                varid_ch = id
                varname = vn
                EXIT
              ENDIF
            END DO
          END IF
        ENDDO
        
        IF (varid_ch == -1) CALL finish(routine, "input file " // TRIM(filename) // &
          & " does not contain a recognizable variable " // &
          & "for chemical heating")
        CALL message(routine, "variable " // TRIM(vn) // " found in input file " // TRIM(filename))
        
        ! check variable dimensionality
        CALL nf(nf90_inquire_variable(ncid, varid_ch, ndims = ndims), routine)
        IF (ndims < 3) CALL finish(routine, "variable " // TRIM(varname) // " in input file " // TRIM(filename) // &
          & " is defined on less than 3 dimensions")
        ALLOCATE(dimids(ndims))
        CALL nf(nf90_inquire_variable(ncid, varid_ch, dimids = dimids), routine)
        IF (dimids(ndims) /= dimid_time) CALL finish(routine, "the 1st dimension of variable " // TRIM(varname) // &
          & " in input file " // TRIM(filename) // " must be time")
        IF (dimids(ndims - 1) /= dimid_lev) CALL finish(routine, "the 2nd dimension of variable " // TRIM(varname) // &
          & " in input file " // TRIM(filename) // " must be level")
        IF (dimids(ndims - 2) /= dimid_lat) CALL finish(routine, "the 3rd dimension of variable " // TRIM(varname) // &
          & " in input file " // TRIM(filename) // " must be latitude")
        DO id = ndims - 3, 1, -1
          CALL nf(nf90_inquire_dimension(ncid, dimids(id), len = dimlen), routine)
          IF (dimlen /= 1) CALL finish(routine, "variable " // TRIM(varname) // " in input file " // TRIM(filename) // &
            & " has non-singleton dimensions other than time, level, and latitude")
        END DO
        DEALLOCATE(dimids)
        
        ! get unit of chemical heating variable
        CALL nf(nf90_inquire_variable(ncid, varid_ch, nAtts = natts), routine)
        units = ''
        DO iatt = 1, natts
          CALL nf(nf90_inq_attname(ncid, varid_ch, iatt, attname), routine)
          IF (TRIM(attname) == 'units') THEN
            CALL nf(nf90_get_att(ncid, varid_ch, 'units', units), routine)
          END IF
        END DO
        scl_ch = heating2kps_scl(units)
        
        IF (scl_ch < 0) THEN
          scl_ch = 1.0_wp
          CALL message(routine, "variable " // TRIM(varname) // " in input file " // TRIM(filename) // &
            & " doesn't have a recognizable unit attribute: assuming K/s")
        ELSE
          WRITE(message_text, *)scl_ch
          CALL message(routine, "variable " // TRIM(varname) // " in input file " // TRIM(filename) // &
            & " has unit " // TRIM(units) // "; conversion factor to K/s is " // TRIM(ADJUSTL(message_text)))
        END IF
        
        ! get type and unit of level variable
        CALL nf(nf90_inquire_variable(ncid, varid_lev, nAtts = natts), routine)
        units = ''
        DO iatt = 1, natts
          CALL nf(nf90_inq_attname(ncid, varid_lev, iatt, attname), routine)
          IF (TRIM(attname) == 'units') THEN
            CALL nf(nf90_get_att(ncid, varid_lev, 'units', units), routine)
          END IF
        END DO
        scl_lev = length2meter_scl(units)      ! try height level (convert to m)
        IF (scl_lev > 0) THEN
          WRITE(message_text, *)scl_lev
          CALL message(routine, "dimension " // TRIM(levname) // " in input file " // TRIM(filename) // &
            & " has unit " // TRIM(units) // "; conversion factor to m is " // TRIM(ADJUSTL(message_text)))
        ELSE
          scl_lev = 1.0_wp
          CALL message(routine, "level variable " // TRIM(levname) // " in input file " // TRIM(filename) // &
            & " doesn't have a recognizable unit attribute: assuming m")
        END IF
        
        CALL nf(nf90_close(ncid), routine)
      END IF
      
      CALL p_bcast(filename, p_io, mpi_comm)
      CALL p_bcast(nlat, p_io, mpi_comm)
      CALL p_bcast(nlev, p_io, mpi_comm)
      CALL p_bcast(varid_ch, p_io, mpi_comm)
      CALL p_bcast(varid_lev, p_io, mpi_comm)
      CALL p_bcast(varid_lat, p_io, mpi_comm)
      CALL p_bcast(varid_time, p_io, mpi_comm)
      CALL p_bcast(scl_ch, p_io, mpi_comm)
      CALL p_bcast(scl_lev, p_io, mpi_comm)
      
    ENDIF
    
    IF (PRESENT(opt_nlat))         opt_nlat         = nlat
    IF (PRESENT(opt_nlev))         opt_nlev         = nlev
    IF (PRESENT(opt_ntime))        opt_ntime        = ntime
    IF (PRESENT(opt_unitchemheat)) opt_unitchemheat = 'K s-1'
    IF (PRESENT(opt_unitlat))      opt_unitlat      = 'rad'
    IF (PRESENT(opt_unitlev))      opt_unitlev      = 'm'
    IF (PRESENT(opt_unittime))     opt_unittime     = 'month'
    
    lchemheat_checked = .TRUE.
  END SUBROUTINE chem_heat_check

  SUBROUTINE chem_heat_init(opt_chemheat, opt_lat, opt_lev, opt_time)

    ! IN/OUT
    REAL(wp), OPTIONAL, INTENT(OUT) :: opt_chemheat(:,:,:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: opt_lat(:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: opt_lev(:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: opt_time(:)

    ! LOCAL
    INTEGER :: mpi_comm, ncid
    CHARACTER(len=*), PARAMETER :: routine = modname // ":chem_heat_init"

    !---------------------------------------------------------

    ! apart from the determination of the start level 
    ! everything in this subroutine should be run through only once
    IF (.NOT. lchemheat_initialized) THEN

      IF (.NOT. lchemheat_checked) CALL chem_heat_check()
      
      ALLOCATE(levclim(nlev))
      ALLOCATE(latclim(nlat))
      ALLOCATE(timeclim(ntime))
      ALLOCATE(chclim(nlat, nlev, ntime))
      
      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF
      
      IF (my_process_is_stdio()) THEN
        ! open
        CALL nf(nf90_open(TRIM(filename), NF90_NOWRITE, ncid), routine)
        
        ! get level
        CALL nf(nf90_get_var(ncid, varid_lev, levclim), routine)
        levclim = levclim * scl_lev       ! unit conversion to m
        minlev = MINVAL(levclim)
        maxlev = MAXVAL(levclim)
        levinc = levclim(1) < levclim(2)
        
        ! get latitude
        CALL nf(nf90_get_var(ncid, varid_lat, latclim), routine)
        latclim = latclim * deg2rad       ! convert from degree to radian
        minlat = MINVAL(latclim)
        maxlat = MAXVAL(latclim)
        latinc = latclim(1) < latclim(2)

        ! get time
        CALL nf(nf90_get_var(ncid, varid_time, timeclim), routine)
        
        ! get chemical heating
        CALL nf(nf90_get_var(ncid, varid_ch, chclim), routine)
        chclim = chclim * scl_ch          ! unit conversion to K/s
        
        ! close
        CALL nf(nf90_close(ncid), routine)
      END IF
      
      CALL p_bcast(minlev, p_io, mpi_comm)
      CALL p_bcast(maxlev, p_io, mpi_comm)
      CALL p_bcast(minlat, p_io, mpi_comm)
      CALL p_bcast(maxlat, p_io, mpi_comm)
      CALL p_bcast(levinc, p_io, mpi_comm)
      CALL p_bcast(latinc, p_io, mpi_comm)
      CALL p_bcast(levclim, p_io, mpi_comm)
      CALL p_bcast(latclim, p_io, mpi_comm)
      CALL p_bcast(timeclim, p_io, mpi_comm)
      CALL p_bcast(chclim, p_io, mpi_comm)

    ENDIF  !lchemheat_initialized

    IF (PRESENT(opt_chemheat)) opt_chemheat = chclim
    IF (PRESENT(opt_lat))      opt_lat      = latclim
    IF (PRESENT(opt_lev))      opt_lev      = levclim
    IF (PRESENT(opt_time))     opt_time     = timeclim
    
    lchemheat_initialized = .TRUE.
  END SUBROUTINE chem_heat_init

  SUBROUTINE chem_heat_clean
    IF (lchemheat_finalized) RETURN
    IF (ALLOCATED(levclim))  DEALLOCATE(levclim)
    IF (ALLOCATED(latclim))  DEALLOCATE(latclim)
    IF (ALLOCATED(timeclim)) DEALLOCATE(timeclim)
    IF (ALLOCATED(chclim))   DEALLOCATE(chclim)
    lchemheat_finalized = .TRUE.
  END SUBROUTINE chem_heat_clean

  SUBROUTINE chem_heat(jcs, jce, kbdim, klev, lat, zh, this_datetime, cheat, effrsw,    &
    &                  opt_istartlev, opt_iendlev, opt_loffline)

    ! IN/OUT
    INTEGER , INTENT(IN) :: jcs, jce, kbdim, klev   ! dimensions
    REAL(wp), INTENT(IN) :: lat(kbdim)              ! latitude                          [rad]
    REAL(wp), INTENT(IN) :: zh(kbdim, klev)         ! geopotential height at full level [m]
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    REAL(wp), INTENT(OUT):: cheat(kbdim,klev)       ! tendency dT/dt (K/s)
    REAL(wp), INTENT(OUT):: effrsw(kbdim,klev)      ! efficiency factor for standard shartwave radiation
    INTEGER,  OPTIONAL, INTENT(IN)  :: opt_istartlev, opt_iendlev ! optional vertical start and end indices
    LOGICAL,  OPTIONAL, INTENT(IN)  :: opt_loffline  ! optional offline mode

    ! LOCAL
    REAL(wp), ALLOCATABLE :: lev(:, :), sclfac(:, :)
    REAL(wp), ALLOCATABLE :: chclim_time(:, :)      ! time-interpolated climatology
    REAL(wp), ALLOCATABLE :: chclim_lat(:, :)       ! latitude-interpolated climatology
    REAL(wp), ALLOCATABLE :: chclim_lev(:, :)       ! level-interpolated climatology
    INTEGER :: jk, jks, istartlev, iendlev, nactivelev
    LOGICAL :: loffline

    !---------------------------------------------------------

    ! determine start and end indices of vertical grid layers,
    ! for which tendencies should be computed
    IF (PRESENT(opt_istartlev)) THEN
      istartlev = MIN(MAX(1, opt_istartlev), klev)
    ELSE
      istartlev = 1
    ENDIF

    IF (PRESENT(opt_iendlev)) THEN
      iendlev = MIN(MAX(1, opt_iendlev), klev)
    ELSE
      iendlev = klev
    ENDIF

    ! in off-line mode tendencies are computed, 
    ! but there is no feedback on the dynamics; 
    ! here, this information is required for the efficiency factor, 
    ! which has to be set to unity in this case
    IF (PRESENT(opt_loffline)) THEN
      loffline = opt_loffline
    ELSE
      loffline = .FALSE.
    ENDIF

    IF (.NOT. lchemheat_initialized) CALL chem_heat_init()

    ! please do not limit range of assignment 
    ! (e.g., cheat(jcs:jce,istartlev:iendlev) = 0._wp)), 
    ! since tendencies have attribute INTENT(OUT)
    cheat(:, :)  = 0._wp
    effrsw(:, :) = 1._wp

    IF (istartlev > iendlev) RETURN
 
    nactivelev = iendlev - istartlev + 1

    ALLOCATE(lev(kbdim, nactivelev), sclfac(kbdim, nactivelev), &
      &      chclim_time(nlat, nlev), chclim_lat(kbdim, nlev),  &
      &      chclim_lev(kbdim, nactivelev))

    lev = zh(:, istartlev : iendlev)

    !******* interpolate in time *******
    chclim_time = intp_anncy(nlat, nlev, this_datetime, chclim)    ! chclim_time is on old lat, old lev, and new time (scalar)

    !******** interpolate in lat ********
    chclim_lat = intp_lat(nlev, nlat, latclim, chclim_time, kbdim, jcs, jce, lat)  ! chclim_lat is on old lev and new lat

    !******** interpolate in lev ********
    ! scale factor: 0 at zeroz and below, 1 at onez and above
    sclfac = MIN(1._wp, MAX(0._wp, (lev - zeroz) / (onez - zeroz)))
    ! cheat is on new lev and new lat
    chclim_lev(:, 1 : nactivelev) = intp_lev(jcs, jce, kbdim, nlev, levclim, chclim_lat, nactivelev, lev)
    ! scale cheat
    DO jk = istartlev, iendlev
      ! shifted index should start counting from 1
      jks = jk - istartlev + 1
      cheat(:, jk) = sclfac(:, jks) * chclim_lev(:, jks)
    ENDDO
    ! in the offline mode we keep the above initialization with 1
    IF (.NOT. loffline) THEN
      DO jk = istartlev, iendlev
        jks = jk - istartlev + 1
        effrsw(:, jk) = 1._wp - sclfac(:, jks) * ( 1._wp - effrswmin)
      ENDDO
    ENDIF

    !******** finalize ********
    IF (ALLOCATED(lev)) DEALLOCATE(lev)
    IF (ALLOCATED(sclfac)) DEALLOCATE(sclfac)
    IF (ALLOCATED(chclim_time)) DEALLOCATE(chclim_time)
    IF (ALLOCATED(chclim_lat)) DEALLOCATE(chclim_lat)
    IF (ALLOCATED(chclim_lev)) DEALLOCATE(chclim_lev)
  END SUBROUTINE chem_heat

  FUNCTION intp_time_weight(curr_datetime, nmonth) RESULT (wi)

    ! IN/OUT
    TYPE(datetime), POINTER, INTENT(IN) :: curr_datetime
    INTEGER, INTENT(IN) :: nmonth               ! length of data: 12 for climatology (annual cycle), 14 for amip (0:13)
    TYPE(t_intp_weight) :: wi

    ! LOCAL
    TYPE(datetime),  POINTER :: prev15, next15, this15
    TYPE(julianday), POINTER :: curr_jul, this15_jul, prev15_jul, next15_jul
    REAL(wp) :: curr_jd, this15_jd, prev15_jd, next15_jd
    REAL(wp) :: diff_this15, diff_prev15, diff_next15_prev15

    CHARACTER(LEN=*), PARAMETER :: routine = modname // ":intp_time_weight"

    !---------------------------------------------------------

    curr_jul => newJulianday(0_i8, 0_i8)
    this15_jul => newJulianday(0_i8, 0_i8)
    prev15_jul => newJulianday(0_i8, 0_i8)
    next15_jul => newJulianday(0_i8, 0_i8)

    ! current julian day
    CALL getJulianDayFromDatetime(curr_datetime, curr_jul)
    curr_jd = curr_jul%day + curr_jul%ms / 86400000

    ! mid-month point of this month and julian day
    
    !write(*,*)'intp_time_weight:',curr_datetime%date%year,',',curr_datetime%date%month

    this15 => newDatetime(curr_datetime%date%year, curr_datetime%date%month, 15, 12, 0, 0, 0)
    CALL getJulianDayFromDatetime(this15, this15_jul)
    this15_jd = this15_jul%day + this15_jul%ms / 86400000

    ! time-difference between now and mid-month point of this month
    diff_this15 = curr_jd - this15_jd

    IF (diff_this15 >= 0) THEN     ! mid-month point or second half of a month
      prev15 => newDatetime(curr_datetime%date%year, curr_datetime%date%month, 15, 12, 0, 0, 0)  ! same as this15

      IF (curr_datetime%date%month == 12) THEN
        next15 => newDateTime(curr_datetime%date%year + 1, 1, 15, 12, 0, 0, 0)
      ELSE
        next15 => newDateTime(curr_datetime%date%year, curr_datetime%date%month + 1, 15, 12, 0, 0, 0)
      END IF
    ELSE                          ! first half of a month
      IF (curr_datetime%date%month == 1) THEN
        prev15 => newDatetime(curr_datetime%date%year - 1, 12, 15, 12, 0, 0, 0)
      ELSE
        prev15 => newDatetime(curr_datetime%date%year, curr_datetime%date%month - 1, 15, 12, 0, 0, 0)
      END IF

      next15 => newDatetime(curr_datetime%date%year, curr_datetime%date%month, 15, 12, 0, 0, 0)  ! same as this15
    END IF

    CALL getJulianDayFromDatetime(prev15, prev15_jul)
    prev15_jd = prev15_jul%day + prev15_jul%ms / 86400000
    CALL getJulianDayFromDatetime(next15, next15_jul)
    next15_jd = next15_jul%day + next15_jul%ms / 86400000

    diff_prev15 = curr_jd - prev15_jd
    diff_next15_prev15 = next15_jd - prev15_jd

    IF (nmonth == 12) THEN
      wi%prev = prev15%date%month
      wi%next = next15%date%month
    ELSEIF (nmonth == 14) THEN
      IF (prev15%date%year == curr_datetime%date%year - 1) THEN
        wi%prev = 1
      ELSE
        wi%prev = prev15%date%month + 1
      END IF
      IF (next15%date%year == curr_datetime%date%year + 1) THEN
        wi%next = 14
      ELSE
        wi%next = next15%date%month + 1
      END IF
    ELSE
      CALL finish(routine, "number of months must be 12 or 14")
    END IF

    wi%factor = diff_prev15 / diff_next15_prev15

    CALL deallocateDatetime(this15)
    CALL deallocateDatetime(prev15)
    CALL deallocateDatetime(next15)
    CALL deallocateJulianday(curr_jul)
    CALL deallocateJulianday(this15_jul)
    CALL deallocateJulianday(prev15_jul)
    CALL deallocateJulianday(next15_jul)
  END FUNCTION intp_time_weight

  FUNCTION intp_lev_weight(kilev, ilev, olev) RESULT(wi)

    ! IN/OUT
    INTEGER,  INTENT(IN) :: kilev
    REAL(wp), INTENT(IN) :: ilev(kilev)   ! input level
    REAL(wp), INTENT(IN) :: olev          ! output level
    TYPE(t_intp_weight)  :: wi

    ! LOCAL
    REAL(wp) :: minilev, maxilev
    LOGICAL :: incilev

    INTEGER :: poslev, prevlev, nextlev
    REAL(wp) :: difflev(kilev)

    !---------------------------------------------------------

    minilev = MINVAL(ilev)
    maxilev = MAXVAL(ilev)
    incilev = ilev(1) < ilev(2)

    IF (olev <= minilev) THEN
      IF (incilev) THEN
        wi%prev = 1
        wi%next = 1
        wi%factor = 0._wp
      ELSE
        wi%prev = kilev
        wi%next = kilev
        wi%factor = 0._wp
      END IF
    ELSEIF (olev >= maxilev) THEN
      IF (incilev) THEN
        wi%prev = kilev
        wi%next = kilev
        wi%factor = 0._wp
      ELSE
        wi%prev = 1
        wi%next = 1
        wi%factor = 0._wp
      END IF
    ELSE
      difflev = ilev - olev

      ! position of the 1st element in ilev greater than olev
      poslev = MINLOC(difflev, 1, difflev >= 0.0_wp)   ! poslev belongs to [2, kilev], 
                                                       ! if incilev == .true. or [1, kilev - 1], 
                                                       ! if incilev == .false.)
      IF (difflev(poslev) == 0.0_wp) THEN
        wi%prev = poslev
        wi%next = poslev
        wi%factor = 0._wp
      ELSE
        nextlev = poslev                               ! position of element of ilev greater than olev(
        IF (incilev) THEN
          prevlev = poslev - 1                         ! position of element of ilev less than olev
        ELSE
          prevlev = poslev + 1
        END IF

        wi%prev = prevlev
        wi%next = nextlev
        wi%factor = (olev - ilev(prevlev)) / (ilev(nextlev) - ilev(prevlev))
      END IF
    END IF
  END FUNCTION intp_lev_weight

  FUNCTION intp_lat_weight(kilat, ilat, olat) RESULT(wi)

    ! IN/OUT
    INTEGER,  INTENT(IN) :: kilat
    REAL(wp), INTENT(IN) :: ilat(kilat)         ! input latitude
    REAL(wp), INTENT(IN) :: olat                ! output latitude
    TYPE(t_intp_weight)  :: wi

    ! LOCAL
    REAL(wp) :: minilat, maxilat
    LOGICAL :: incilat

    INTEGER :: poslat, prevlat, nextlat
    REAL(wp) :: difflat(kilat)

    !---------------------------------------------------------

    minilat = MINVAL(ilat)
    maxilat = MAXVAL(ilat)
    incilat = ilat(1) < ilat(2)

      IF (olat <= minilat) THEN
        IF (incilat) THEN
          wi%prev = 1
          wi%next = 1
          wi%factor = 0._wp
        ELSE
          wi%prev = kilat
          wi%next = kilat
          wi%factor = 0._wp
        END IF
      ELSEIF (olat >= maxilat) THEN
        IF (incilat) THEN
          wi%prev = kilat
          wi%next = kilat
          wi%factor = 0._wp
        ELSE
          wi%prev = 1
          wi%next = 1
          wi%factor = 0._wp
        END IF
      ELSE
        difflat = ilat - olat

        ! position of the 1st element in ilat greater than olat
        poslat = MINLOC(difflat, 1, difflat >= 0.0_wp)   ! poslat belongs to [2, kilat], 
                                                         ! if incilat == .true. or [1, kilat - 1], 
                                                         ! if incilat == .false.)
        IF (difflat(poslat) == 0.0_wp) THEN
          wi%prev = poslat
          wi%next = poslat
          wi%factor = 0._wp
        ELSE
          nextlat = poslat                               ! position of element of ilat greater than olat
          IF (incilat) THEN
            prevlat = poslat - 1                         ! position of element of ilat less than olat
          ELSE
            prevlat = poslat + 1
          END IF

          wi%prev = prevlat
          wi%next = nextlat
          wi%factor = (olat - ilat(prevlat)) / (ilat(nextlat) - ilat(prevlat))
        END IF
      END IF
  END FUNCTION intp_lat_weight

  FUNCTION intp_anncy(klat, klev, this_datetime, anncy) RESULT (intp)

    ! IN/OUT
    INTEGER, INTENT(IN) :: klat, klev
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
    REAL(wp), INTENT(IN) :: anncy(klat, klev, 12)   ! monthly-mean annual cycle
    REAL(wp) :: intp(klat, klev)                    ! interpolated value at current time step

    ! LOCAL
    TYPE(t_intp_weight) :: wi

    !---------------------------------------------------------

    wi = intp_time_weight(this_datetime, 12)

    intp(:, :) = anncy(:, :, wi%prev) + wi%factor * (anncy(:, :, wi%next) - anncy(:, :, wi%prev))
  END FUNCTION intp_anncy

  FUNCTION intp_lev(jcs, jce, kbdim, kilev, ilev, iprof, kolev, olev) RESULT(oprof)

    ! IN/OUT
    INTEGER,  INTENT(IN) :: jcs, jce, kbdim, kilev, kolev
    REAL(wp), INTENT(IN) :: ilev(kilev)          ! input level
    REAL(wp), INTENT(IN) :: iprof(kbdim, kilev)  ! input profile
    REAL(wp), INTENT(IN) :: olev(kbdim, kolev)   ! output level
    REAL(wp) :: oprof(kbdim, kolev)              ! output profile

    ! LOCAL
    TYPE(t_intp_weight) :: wi
    INTEGER :: jl, jk

    !---------------------------------------------------------

    oprof(:,:) = 0._wp

    DO jl = jcs, jce
      DO jk = 1, kolev
        wi = intp_lev_weight(kilev, ilev, olev(jl, jk))

        oprof(jl, jk) = iprof(jl, wi%prev) + wi%factor * (iprof(jl, wi%next) - iprof(jl, wi%prev))
      END DO
    END DO
  END FUNCTION intp_lev

  FUNCTION intp_lat(klev, kilat, ilat, iprof, kolat, solat, eolat, olat) RESULT(oprof)

    ! IN/OUT
    INTEGER,  INTENT(IN) :: klev, kilat, kolat
    REAL(wp), INTENT(IN) :: ilat(kilat)            ! input latitude
    REAL(wp), INTENT(IN) :: iprof(kilat, klev)     ! input profile
    INTEGER,  INTENT(IN) :: solat, eolat           ! start and end index of valid input latitude
    REAL(wp), INTENT(IN) :: olat(kolat)            ! output latitude
    REAL(wp) :: oprof(kolat, klev)                 ! output profile

    ! LOCAL
    INTEGER :: jl
    TYPE(t_intp_weight) :: wi

    !---------------------------------------------------------

    oprof(:, :) = 0._wp

    DO jl = solat, eolat
      wi = intp_lat_weight(kilat, ilat, olat(jl))

      oprof(jl, :) = iprof(wi%prev, :) + wi%factor * (iprof(wi%next, :) - iprof(wi%prev, :))
    END DO
  END FUNCTION intp_lat

  FUNCTION heating2kps_scl(units, panic) RESULT(scl)
    
    ! IN/OUT
    CHARACTER(LEN=*),  INTENT(IN) :: units
    LOGICAL, OPTIONAL, INTENT(IN) :: panic

    REAL(wp) :: scl

    !---------------------------------------------------------

    SELECT CASE (TRIM(tolower(units)))
    CASE ('k/s', 'k s-1', 'k*s-1', 'k s^-1', 'k*s^-1')
      scl = 1.0_wp
    CASE ('k/day', 'k day-1', 'k*day-1', 'k day^-1', 'k*day^-1')
      scl = 1.0_wp / 86400.0_wp
    CASE DEFAULT
      scl = -1
      IF (present(panic) .AND. panic) CALL finish('heating2kps_scl: invalid heating rate unit: "'//TRIM(units)//'"')
    END SELECT
  END FUNCTION heating2kps_scl

  FUNCTION length2meter_scl(units, panic) RESULT (scl)

    ! IN/OUT
    CHARACTER(LEN=*), INTENT(IN) :: units
    LOGICAL, OPTIONAL, INTENT(IN) :: panic

    REAL(wp) :: scl

    !---------------------------------------------------------

    SELECT CASE (TRIM(units))
    CASE('Ym', 'yottameter', 'yottametre')
      scl = 1e24_wp
    CASE('Zm', 'zettameter', 'zettametre')
      scl = 1e21_wp
    CASE('Em', 'exameter', 'exametre')
      scl = 1e18_wp
    CASE('Pm', 'petameter', 'petametre')
      scl = 1e15_wp
    CASE('Tm', 'terameter', 'terametre')
      scl = 1e12_wp
    CASE('Gm', 'gigameter', 'gigametre')
      scl = 1e9_wp
    CASE('Mm', 'megameter', 'megametre')
      scl = 1e6_wp
    CASE('league')
      scl = 3 * 1609.344_wp
    CASE ('nmi', 'nautical mile')
      scl = 1852.0_wp
    CASE('mile')
      scl = 1609.344_wp
    CASE('km', 'kilometer', 'kilometre')
      scl = 1000.0_wp
    CASE('hm', 'hectometer', 'hectometre')
      scl = 100.0_wp
    CASE('dam', 'decameter', 'decametre')
      scl = 10.0_wp
    CASE ('m', 'meter', 'metre')
      scl = 1.0_wp
    CASE ('yd', 'yard')
      scl = 0.9144_wp
    CASE ('ft', 'foot', 'feet')
      scl = 0.3048_wp
    CASE ('dm', 'decimeter', 'decimetre')
      scl = 0.1_wp
    CASE ('in', 'inch')
      scl = 0.0254_wp
    CASE ('cm', 'centimeter', 'centimetre')
      scl = 0.01_wp
    CASE ('mm', 'millimeter', 'millimetre')
      scl = 0.001_wp
    CASE ('thou', 'mil')
      scl = 2.54e-5_wp
    CASE ('um', 'micrometer', 'micrometre', 'micron')
      scl = 1e-6_wp
    CASE ('nm', 'nanometer', 'nanometre')
      scl = 1e-9_wp
    CASE ('A', 'angstrom')
      scl = 1e-10_wp
    CASE ('pm', 'picometer', 'picometre')
      scl = 1e-12_wp
    CASE ('xu')
      scl = 1e-13_wp
    CASE ('fm', 'femtometer', 'femtometre', 'fermi')
      scl = 1e-15_wp
    CASE ('am', 'attometer', 'attometre')
      scl = 1e-18_wp
    CASE ('zm', 'zeptometer', 'zeptometre')
      scl = 1e-21_wp
    CASE ('ym', 'yoctometer', 'yoctometre')
      scl = 1e-24_wp
    CASE DEFAULT
      scl = -1
      IF (present(panic) .AND. panic) CALL finish('length2meter_scl: invalid length unit: "'//TRIM(units)//'"')
    END SELECT
  END FUNCTION length2meter_scl

  FUNCTION generate_filename(filename, model_base_dir) RESULT(result_str)

    ! IN/OUT
    CHARACTER(len=*), INTENT(IN)    :: filename, &
      &                                model_base_dir
    CHARACTER(len=MAX_CHAR_LENGTH)  :: result_str

    ! LOCAL
    TYPE (t_keyword_list), POINTER  :: keywords => NULL()

    !---------------------------------------------------------

    CALL associate_keyword("<path>",     TRIM(model_base_dir), keywords)
    ! replace keywords in "filename", which is by default
    ! filename = "<path>chemheat"
    result_str = TRIM(with_keywords(keywords, TRIM(filename)))

  END FUNCTION generate_filename

END MODULE mo_upatmo_phy_chemheat
