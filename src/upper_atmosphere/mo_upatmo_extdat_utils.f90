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

! Auxiliary subroutines for the processing
! of the external upper-atmosphere data.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_upatmo_extdat_utils

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message_text, message
  USE mo_impl_constants,         ONLY: SUCCESS, MAX_CHAR_LENGTH, &
    &                                  min_rlcell_int
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c
  USE mo_math_constants,         ONLY: deg2rad, pi_2
  USE mo_upatmo_impl_const,      ONLY: iUpatmoExtdatLatId, iUpatmoExtdatLevId, &
    &                                  iUpatmoExtdatTimeId
  USE mo_model_domain,           ONLY: t_patch
  USE mo_upatmo_types,           ONLY: t_extdat_latlevtime
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_mpi,                    ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                                  p_comm_work
  USE mo_netcdf_errhandler,      ONLY: nf
  USE mo_netcdf
  USE mo_upatmo_phy_chemheat,    ONLY: chem_heat_check, chem_heat_init
  USE mo_util_string,            ONLY: int2string, real2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_extdat_gas
  PUBLIC :: read_extdat_chemheat
  PUBLIC :: construct_interpolation_lat
  PUBLIC :: construct_interpolation_lev

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_upatmo_extdat_utils'

CONTAINS

  !>
  !! Read external gas data for the upper atmosphere 
  !! under NWP forcing.
  !!
  SUBROUTINE read_extdat_gas( gas,         &  !inout
    &                         gasname,     &  !in
    &                         vmr2mmr,     &  !in
    &                         filename,    &  !in
    &                         opt_lmessage )  !optin

    ! In/out variables
    TYPE(t_extdat_latlevtime), INTENT(INOUT) :: gas
    CHARACTER(LEN=*),          INTENT(IN)    :: gasname
    REAL(wp),                  INTENT(IN)    :: vmr2mmr
    CHARACTER(LEN=*),          INTENT(IN)    :: filename
    LOGICAL, OPTIONAL,         INTENT(IN)    :: opt_lmessage

    ! Local variables
    INTEGER, ALLOCATABLE :: dimids(:)
    INTEGER  :: ndims(3)
    INTEGER  :: ncid
    INTEGER  :: dimid_time, dimid_lev, dimid_lat
    INTEGER  :: ntime, nlev, nlat, ndim
    INTEGER  :: varid_time, varid_lev, varid_lat, varid_gas
    INTEGER  :: i, istart, iend, istep, jlev, jlat, jtime
    INTEGER  :: istat
    INTEGER  :: mpi_comm

    LOGICAL  :: lstdioproc, lexist, lmessage

    CHARACTER(LEN = :), ALLOCATABLE :: varunit, dimunit_lat, dimunit_lev, dimunit_time

    CHARACTER(LEN=4),  PARAMETER :: levname  = 'plev'
    CHARACTER(LEN=3),  PARAMETER :: latname  = 'lat'
    CHARACTER(LEN=4),  PARAMETER :: timename = 'time'
    CHARACTER(LEN=2),  PARAMETER :: levunit  = 'Pa' 
    CHARACTER(LEN=13), PARAMETER :: latunit  = 'degrees_north'
    CHARACTER(LEN=15), PARAMETER :: timeunit = 'month_of_a_year'
    CHARACTER(LEN=9),  PARAMETER :: gasunit  = 'mol mol-1' 
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':read_extdat_gas'
    
    !-------------------------------------------------------------- 

    !--------
    ! Checks
    !--------

    IF (gas%lat_id /= iUpatmoExtdatLatId%deg) THEN
      CALL finish(TRIM(routine), 'Dimension of latitude has to be: degree north.')
    ELSEIF (gas%lev_id /= iUpatmoExtdatLevId%p) THEN
      CALL finish(TRIM(routine), 'Levels have to be pressure levels.')
    ELSEIF (gas%time_id /= iUpatmoExtdatTimeId%month) THEN
      CALL finish(TRIM(routine), 'Dimension of time has to be: month .')
    ENDIF

    IF (PRESENT(opt_lmessage)) THEN
      lmessage = opt_lmessage
    ELSE
      lmessage = .FALSE.
    ENDIF

    IF (lmessage) CALL message(TRIM(routine), &
      & 'Processing external data for gas '//TRIM(gasname))

    ! Initializations
    ntime = -999
    nlev  = -999
    nlat  = -999

    ! Please do not remove this!
    ! (There are some problems with regard to the read-in of strings 
    ! from the default file with the external data. 
    ! The reasons are not yet know (maybe it is related to numbers in strings). 
    ! This is a workaround for the time being.)
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(gasunit)) :: varunit, STAT=istat)
    IF (istat/=SUCCESS) CALL finish (TRIM(routine), 'Allocation of varunit failed')
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(latunit)) :: dimunit_lat, STAT=istat)
    IF (istat/=SUCCESS) CALL finish (TRIM(routine), 'Allocation of dimunit_lat failed')
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(levunit)) :: dimunit_lev, STAT=istat)
    IF (istat/=SUCCESS) CALL finish (TRIM(routine), 'Allocation of dimunit_lev failed')
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(timeunit)) :: dimunit_time, STAT=istat)
    IF (istat/=SUCCESS) CALL finish (TRIM(routine), 'Allocation of dimunit_time failed')

    ! Is this the standard I/O-process?
    lstdioproc = my_process_is_stdio()

    ! MPI communicator
    mpi_comm = p_comm_work

    IF (lstdioproc) THEN

      !------------------
      ! Does file exist?
      !------------------

      INQUIRE(file=TRIM(filename), EXIST=lexist)
      IF (.NOT. lexist) THEN
        message_text = 'Gas file '//TRIM(filename)//' does not exist.'
        CALL finish(TRIM(routine), TRIM(message_text))
      ENDIF

      !-----------
      ! Open file
      !-----------

      CALL nf(nf90_open(TRIM(filename), NF90_NOWRITE, ncid), routine)

      !---------------------
      ! Evaluate dimensions
      !---------------------

      ! Time
      CALL get_dim( ncid    = ncid,         & !in
        &           dimname = timename,     & !in
        &           dimunit = dimunit_time, & !inout
        &           dimid   = dimid_time,   & !inout
        &           ndim    = ntime,        & !inout
        &           varid   = varid_time    ) !inout

      ! Levels
      CALL get_dim( ncid    = ncid,        & !in
        &           dimname = levname,     & !in
        &           dimunit = dimunit_lev, & !inout
        &           dimid   = dimid_lev,   & !inout
        &           ndim    = nlev,        & !inout
        &           varid   = varid_lev    ) !inout

      ! Latitudes
      CALL get_dim( ncid    = ncid,        & !in
        &           dimname = latname,     & !in
        &           dimunit = dimunit_lat, & !inout
        &           dimid   = dimid_lat,   & !inout
        &           ndim    = nlat,        & !inout
        &           varid   = varid_lat    ) !inout

      ! Some checks
      IF (ntime < 1) THEN
        message_text = 'Size of time dimension in gas file '//TRIM(filename) &
          & //' has to be > 0, but it is '// TRIM(int2string(ntime))
        CALL finish(TRIM(routine), TRIM(message_text))
      ELSEIF (nlev < 1) THEN
        message_text = 'Size of level dimension in gas file '//TRIM(filename) &
          & //' has to be > 0, but it is '// TRIM(int2string(nlev))
        CALL finish(TRIM(routine), TRIM(message_text))
      ELSEIF (nlat < 1) THEN
        message_text = 'Size of latitude dimension in gas file '//TRIM(filename) &
          & //' has to be > 0, but it is '// TRIM(int2string(nlat))
        CALL finish(TRIM(routine), TRIM(message_text))
      ELSEIF (ntime /= 12) THEN
        message_text = 'Size of time dimension in gas file '//TRIM(filename) &
          & //' has to be 12, but it is '// TRIM(int2string(ntime))
        CALL finish(TRIM(routine), TRIM(message_text))
      ELSEIF (TRIM(dimunit_time) /= timeunit) THEN
        message_text = 'Exclusively supportet time unit: '//timeunit &
          & //' but unit in file is: '//TRIM(dimunit_time)
        CALL finish(TRIM(routine), TRIM(message_text))
      ELSEIF (TRIM(dimunit_lev) /= levunit) THEN
        message_text = 'Exclusively supportet level unit: '//levunit &
          & //' but unit in file is: '//TRIM(dimunit_lev)
        CALL finish(TRIM(routine), TRIM(message_text))
      ELSEIF (TRIM(dimunit_lat) /= latunit) THEN
        message_text = 'Exclusively supportet latitude unit: '//latunit &
          & //' but unit in file is: '//TRIM(dimunit_lat)
        CALL finish(TRIM(routine), TRIM(message_text))
      ENDIF

      ! Gas:
      ! If the variable could not be found, nf would call 'finish'
      CALL nf(nf90_inq_varid(ncid, TRIM(gasname), varid_gas), routine)
      ! Get number of variable dimensions
      CALL nf(nf90_inquire_variable(ncid, varid_gas, ndims = ndim), routine)
      ! Number of dimensions should be 3 (or maybe larger)
      IF (ndim < 3) THEN
        message_text = 'Gas '//TRIM(gasname)//' in gas file '//TRIM(filename) &
          & //' has less than 3 dimensions.' 
        CALL finish(TRIM(routine), TRIM(message_text))
      ENDIF  !IF (ndim < 3)
      ALLOCATE(dimids(ndim), STAT=istat)
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of dimids failed.')
      CALL nf(nf90_inquire_variable(ncid, varid_gas, dimids = dimids), routine)
      ! Check, if variable varies in correct dimensions
      IF (dimids(ndim) /= dimid_time) THEN
        message_text = 'First dimension of gas '//TRIM(gasname) &
          & //' in gas file '//TRIM(filename)//' needs to be time.' 
        CALL finish(TRIM(routine), TRIM(message_text))
      ELSEIF (dimids(ndim-1) /= dimid_lev) THEN
        message_text = 'Second dimension of gas '//TRIM(gasname) &
          & //' in gas file '//TRIM(filename)//' needs to be level.' 
        CALL finish(TRIM(routine), TRIM(message_text))
      ELSEIF (dimids(ndim-2) /= dimid_lat) THEN
        message_text = 'Third dimension of gas '//TRIM(gasname) &
          & //' in gas file '//TRIM(filename)//' needs to be latitude.' 
        CALL finish(TRIM(routine), TRIM(message_text))
      ENDIF
      DEALLOCATE(dimids, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of dimids failed.')     
      ! Check variable unit
      CALL nf(nf90_get_att(ncid, varid_gas, 'units', varunit), routine)
      IF (TRIM(varunit) /= gasunit) THEN
        message_text = 'Exclusively supportet gas unit: '//gasunit &
          & //' but unit in file is: '//TRIM(varunit)
        CALL finish(TRIM(routine), TRIM(message_text))
      ENDIF
 
    ENDIF  !IF (lstdioproc)

    !---------------------------
    ! Broadcast dimension sizes
    !---------------------------

    ndims = (/nlat, nlev, ntime/)
    CALL p_bcast(ndims, p_io, mpi_comm)
    nlat  = ndims(1)
    nlev  = ndims(2)
    ntime = ndims(3)

    ! Store sizes
    gas%nlat   = nlat
    gas%nlev   = nlev
    gas%ntime  = ntime

    ! Allocate fields in external data type
    ALLOCATE( gas%data(nlat,nlev,ntime), &
      &       gas%lat(nlat),             &
      &       gas%lev(nlev),             &
      &       gas%lev_half(nlev+1),      &
      &       gas%time(ntime),           &
      &       STAT=istat                 )
    IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of gas failed.')

    !-----------------
    ! Broadcast units
    !-----------------

    ! Below, we will convert the gas unit 
    ! from volume mixing ratio (mole fraction) 
    ! into mass mixing ratio
    gas%unit_data = 'kg kg-1'
    ! Below, we will convert latitudes to rad
    gas%unit_lat  = 'rad'
    gas%unit_lev  = levunit
    gas%unit_time = timeunit

    !--------------------------
    ! Read dimensions and data
    !--------------------------

    IF (lstdioproc) THEN

      ! Read gas data
      CALL nf(nf90_get_var(ncid, varid_gas, gas%data), routine)
      ! Read levels
      CALL nf(nf90_get_var(ncid, varid_lev, gas%lev), routine)
      ! Read latitudes
      CALL nf(nf90_get_var(ncid, varid_lat, gas%lat), routine)
      ! Read times
      CALL nf(nf90_get_var(ncid, varid_time, gas%time), routine)

      !------------
      ! Close file
      !------------

      CALL nf(nf90_close(ncid), routine)

    ENDIF  !IF (lstdioproc)

    !----------------
    ! Broadcast data
    !----------------

    CALL p_bcast(gas%data, p_io, mpi_comm)
    CALL p_bcast(gas%lev, p_io, mpi_comm)
    CALL p_bcast(gas%lat, p_io, mpi_comm)
    CALL p_bcast(gas%time, p_io, mpi_comm)

    !----------------
    ! Postprocessing
    !----------------

    ! Determine start and end indices and steps for loops:
    ! Times: Jan., Feb., ..., Dec.
    IF (gas%time(1) < gas%time(ntime)) THEN
      gas%istarttime = 1
      gas%iendtime   = ntime
      gas%isteptime  = 1
    ELSE
      gas%istarttime = ntime
      gas%iendtime   = 1
      gas%isteptime  = -1
    ENDIF
    ! Whatever the read-in times are, 
    ! we overwrite it with 1, 2, 3, ..., 12
    istart = gas%istarttime
    iend   = gas%iendtime
    istep  = gas%isteptime
    gas%time = (/ (REAL(i,wp),i=istart,iend,istep) /)
    ! Latitudes from south to north
    IF (gas%lat(1) < gas%lat(nlat)) THEN
      gas%istartlat = 1
      gas%iendlat   = nlat
      gas%isteplat  = 1
    ELSE
      gas%istartlat = nlat
      gas%iendlat   = 1
      gas%isteplat  = -1
    ENDIF
    ! Boundary check
    IF (MAXVAL(gas%lat) > 90._wp) &
      & CALL finish(TRIM(routine), 'Max(lat) = '//TRIM(real2string(MAXVAL(gas%lat))))
    IF (MINVAL(gas%lat) < -90._wp) &
      & CALL finish(TRIM(routine), 'Min(lat) = '//TRIM(real2string(MINVAL(gas%lat))))
    ! Convert latitudes from degree north to radian
    gas%lat = gas%lat * deg2rad
    ! Pressure levels from model top to model bottom
    IF (gas%lev(1) < gas%lev(nlev)) THEN
      gas%istartlev = 1
      gas%iendlev   = nlev
      gas%isteplev  = 1
    ELSE
      gas%istartlev = nlev
      gas%iendlev   = 1
      gas%isteplev  = -1
    ENDIF

    ! The half-level pressures, p(jk+-1/2), 
    ! are defined:
    !
    !   p(jk+1/2) = [p(jk) + p(jk+1)] / 2, 
    !
    ! and halve the mass per unit area between p(jk+1) and p(jk). 
    ! The uppermost and lowermost half-level pressures are defined:
    !
    !   p(1-1/2)    = 0 Pa, 
    !   p(nlev+1/2) = 125000 Pa,
    !
    ! so we start with a corresponding bounary check of the read in full pressure levels.
    IF (MAXVAL(gas%lev) > 125000._wp) &
      & CALL finish(TRIM(routine), 'Max(lev) = '//TRIM(real2string(MAXVAL(gas%lev))))
    IF (MINVAL(gas%lev) < 0._wp) &
      & CALL finish(TRIM(routine), 'Min(lev) = '//TRIM(real2string(MINVAL(gas%lev))))

    istep  = gas%isteplev
    istart = gas%istartlev + ( 1 - istep ) / 2
    iend   = gas%iendlev + ( 1 + istep ) / 2

    gas%lev_half(istart) = 0._wp
    gas%lev_half(iend)   = 125000._wp
    DO jlev = istart + istep, iend - istep, istep
      ! Half-level pressure
      gas%lev_half(jlev) = 0.5_wp *( gas%lev(jlev-1) + gas%lev(jlev) )
    ENDDO  !jlev

    ! Convert gas unit from volume mixing ratio (mole fraction) [mol mol-1]
    ! into mass mixing ratio [kg kg-1]
    DO jtime = gas%istarttime, gas%iendtime, gas%isteptime
      DO jlev = gas%istartlev, gas%iendlev, gas%isteplev
        DO jlat = gas%istartlat, gas%iendlat, gas%isteplat
          gas%data(jlat,jlev,jtime) = vmr2mmr * gas%data(jlat,jlev,jtime)
        ENDDO  !jlat
      ENDDO  !jlev
    ENDDO  !jtime

    !------
    ! Info
    !------

    IF (lmessage) THEN 
      CALL message(TRIM(routine), 'istarttime, iendtime, isteptime: ' &
        & //TRIM(int2string(gas%istarttime))//', '                    &
        & //TRIM(int2string(gas%iendtime))//', '                      &
        & //TRIM(int2string(gas%isteptime))                           )
      CALL message(TRIM(routine), 'istartlat, iendlat, isteplat: ' &
        & //TRIM(int2string(gas%istartlat))//', '                  &
        & //TRIM(int2string(gas%iendlat))//', '                    &
        & //TRIM(int2string(gas%isteplat))                         )
      CALL message(TRIM(routine), 'istartlev, iendlev, isteplev: ' &
        & //TRIM(int2string(gas%istartlev))//', '                  &
        & //TRIM(int2string(gas%iendlev))//', '                    &
        & //TRIM(int2string(gas%isteplev))                         )
    ENDIF

    !----------
    ! Clean-up
    !----------

    DEALLOCATE(varunit, dimunit_lat, dimunit_lev, dimunit_time, STAT=istat)
    IF (istat/=SUCCESS) CALL finish (TRIM(routine), 'Deallocation of characters failed')

  END SUBROUTINE read_extdat_gas

  !====================================================================================

  !>
  !! Read external gas data for the upper atmosphere 
  !! under NWP forcing.
  !!
  SUBROUTINE read_extdat_chemheat( chemheat,    &  !inout
    &                              filename,    &  !in
    &                              opt_lmessage )  !optin

    ! In/out variables
    TYPE(t_extdat_latlevtime), INTENT(INOUT) :: chemheat
    CHARACTER(LEN=*),          INTENT(IN)    :: filename
    LOGICAL, OPTIONAL,         INTENT(IN)    :: opt_lmessage

    ! Local variables
    INTEGER  :: istat, i, istart, iend, istep

    LOGICAL :: lmessage

    CHARACTER(LEN = :), ALLOCATABLE :: varunit, dimunit_lat, dimunit_lev, dimunit_time

    CHARACTER(LEN=1), PARAMETER :: levunit      = 'm' 
    CHARACTER(LEN=3), PARAMETER :: latunit      = 'rad'
    CHARACTER(LEN=5), PARAMETER :: timeunit     = 'month'
    CHARACTER(LEN=5), PARAMETER :: chemheatunit = 'K s-1' 
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':read_extdat_chemheat'
    
    !-------------------------------------------------------------- 

    !--------
    ! Checks
    !--------

    IF (chemheat%lat_id /= iUpatmoExtdatLatId%deg) THEN
      CALL finish(TRIM(routine), 'Dimension of latitude has to be: degree north.')
    ELSEIF (chemheat%lev_id /= iUpatmoExtdatLevId%z) THEN
      CALL finish(TRIM(routine), 'Levels have to be geometric heights.')
    ELSEIF (chemheat%time_id /= iUpatmoExtdatTimeId%month) THEN
      CALL finish(TRIM(routine), 'Dimension of time has to be: month .')
    ENDIF

    IF (PRESENT(opt_lmessage)) THEN
      lmessage = opt_lmessage
    ELSE
      lmessage = .FALSE.
    ENDIF

    IF (lmessage) CALL message(TRIM(routine), &
      & 'Processing external data for chemical heating tendencies')

    !--------------------------
    ! Get dimensions and units
    !--------------------------

    ! Please do not remove this!
    ! (There are some problems with regard to the read-in of strings 
    ! from the default file with the external data. 
    ! The reasons are not yet know. This is a workaround for the time being.)
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(chemheatunit)) :: varunit, STAT=istat)
    IF (istat/=SUCCESS) CALL finish (TRIM(routine), 'Allocation of varunit failed')
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(latunit)) :: dimunit_lat, STAT=istat)
    IF (istat/=SUCCESS) CALL finish (TRIM(routine), 'Allocation of dimunit_lat failed')
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(levunit)) :: dimunit_lev, STAT=istat)
    IF (istat/=SUCCESS) CALL finish (TRIM(routine), 'Allocation of dimunit_lev failed')
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(timeunit)) :: dimunit_time, STAT=istat)
    IF (istat/=SUCCESS) CALL finish (TRIM(routine), 'Allocation of dimunit_time failed')

    CALL chem_heat_check( opt_filename     = filename,       & !optin
      &                   opt_nlat         = chemheat%nlat,  & !optout
      &                   opt_nlev         = chemheat%nlev,  & !optout
      &                   opt_ntime        = chemheat%ntime, & !optout
      &                   opt_unitchemheat = varunit,        & !optout
      &                   opt_unitlat      = dimunit_lat,    & !optout
      &                   opt_unitlev      = dimunit_lev,    & !optout
      &                   opt_unittime     = dimunit_time    ) !optout

    ! Some checks
    IF (chemheat%ntime < 1) THEN
      message_text = 'Size of time dimension in chemheat file '//TRIM(filename) &
        & //' has to be > 0, but it is '// TRIM(int2string(chemheat%ntime))
      CALL finish(TRIM(routine), TRIM(message_text))
    ELSEIF (chemheat%nlev < 1) THEN
      message_text = 'Size of level dimension in chemheat file '//TRIM(filename) &
        & //' has to be > 0, but it is '// TRIM(int2string(chemheat%nlev))
      CALL finish(TRIM(routine), TRIM(message_text))
    ELSEIF (chemheat%nlat < 1) THEN
      message_text = 'Size of latitude dimension in chemheat file '//TRIM(filename) &
        & //' has to be > 0, but it is '// TRIM(int2string(chemheat%nlat))
      CALL finish(TRIM(routine), TRIM(message_text))
    ELSEIF (chemheat%ntime /= 12) THEN
      message_text = 'Size of time dimension in chemheat file '//TRIM(filename) &
        & //' has to be 12, but it is '// TRIM(int2string(chemheat%ntime))
      CALL finish(TRIM(routine), TRIM(message_text))
    ELSEIF (dimunit_time /= timeunit) THEN
      message_text = 'Exclusively supportet time unit: '//timeunit &
        & //' but unit in file is: '//TRIM(chemheat%unit_time)
      CALL finish(TRIM(routine), TRIM(message_text))
    ELSEIF (dimunit_lev /= levunit) THEN
      message_text = 'Exclusively supportet level unit: '//levunit &
        & //' but unit in file is: '//TRIM(chemheat%unit_lev)
      CALL finish(TRIM(routine), TRIM(message_text))
    ELSEIF (dimunit_lat /= latunit) THEN
      message_text = 'Exclusively supportet latitude unit: '//latunit &
        & //' but unit in file is: '//TRIM(chemheat%unit_lat)
      CALL finish(TRIM(routine), TRIM(message_text))
    ELSEIF (varunit /= chemheatunit) THEN
      message_text = 'Exclusively supportet gas unit: '//chemheatunit &
        & //' but unit in file is: '//TRIM(chemheat%unit_data)
      CALL finish(TRIM(routine), TRIM(message_text))
    ENDIF

    chemheat%unit_time = dimunit_time
    chemheat%unit_lev  = dimunit_lev
    chemheat%unit_lat  = dimunit_lat
    chemheat%unit_data = varunit

    !----------------------
    ! Allocate data fields
    !----------------------

    ALLOCATE( chemheat%data(chemheat%nlat,chemheat%nlev,chemheat%ntime), &
      &       chemheat%lat(chemheat%nlat),                               &
      &       chemheat%lev(chemheat%nlev),                               &
      &       chemheat%time(chemheat%ntime),                             &
      &       STAT=istat                                                 )
    IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of gas failed.')
    
    !----------
    ! Get data
    !----------

    CALL chem_heat_init( opt_chemheat = chemheat%data, & !optout
      &                  opt_lat      = chemheat%lat,  & !optout
      &                  opt_lev      = chemheat%lev,  & !optout
      &                  opt_time     = chemheat%time  ) !optout

    !----------------
    ! Postprocessing
    !----------------

    ! Determine start and end indices and steps for loops:
    ! Times: Jan., Feb., ..., Dec.
    IF (chemheat%time(1) < chemheat%time(chemheat%ntime)) THEN
      chemheat%istarttime = 1
      chemheat%iendtime   = chemheat%ntime
      chemheat%isteptime  = 1
    ELSE
      chemheat%istarttime = chemheat%ntime
      chemheat%iendtime   = 1
      chemheat%isteptime  = -1
    ENDIF
    ! Whatever the read-in times are, 
    ! we overwrite it with 1, 2, 3, ..., 12
    istart = chemheat%istarttime
    iend   = chemheat%iendtime
    istep  = chemheat%isteptime
    chemheat%time = (/ (REAL(i,wp),i=istart,iend,istep) /)
    ! Latitudes from south to north
    IF (chemheat%lat(1) < chemheat%lat(chemheat%nlat)) THEN
      chemheat%istartlat = 1
      chemheat%iendlat   = chemheat%nlat
      chemheat%isteplat  = 1
    ELSE
      chemheat%istartlat = chemheat%nlat
      chemheat%iendlat   = 1
      chemheat%isteplat  = -1
    ENDIF
    ! Boundary check (latitudes are already in rad)
    IF (MAXVAL(chemheat%lat) > pi_2) &
      & CALL finish(TRIM(routine), 'Max(lat) = '//TRIM(real2string(MAXVAL(chemheat%lat))))
    IF (MINVAL(chemheat%lat) < -pi_2) &
      & CALL finish(TRIM(routine), 'Min(lat) = '//TRIM(real2string(MINVAL(chemheat%lat))))
    ! Geometric height from top to bottom
    IF (chemheat%lev(1) > chemheat%lev(chemheat%nlev)) THEN
      chemheat%istartlev = 1
      chemheat%iendlev   = chemheat%nlev
      chemheat%isteplev  = 1
    ELSE
      chemheat%istartlev = chemheat%nlev
      chemheat%iendlev   = 1
      chemheat%isteplev  = -1
    ENDIF
    ! Boundary check
    IF (MINVAL(chemheat%lev) < 0._wp) &
      & CALL finish(TRIM(routine), 'Min(lev) = '//TRIM(real2string(MINVAL(chemheat%lev))))

    !------
    ! Info
    !------

    IF (lmessage) THEN 
      CALL message(TRIM(routine), 'istarttime, iendtime, isteptime: ' &
        & //TRIM(int2string(chemheat%istarttime))//', '               &
        & //TRIM(int2string(chemheat%iendtime))//', '                 &
        & //TRIM(int2string(chemheat%isteptime))                      )
      CALL message(TRIM(routine), 'istartlat, iendlat, isteplat: ' &
        & //TRIM(int2string(chemheat%istartlat))//', '             &
        & //TRIM(int2string(chemheat%iendlat))//', '               &
        & //TRIM(int2string(chemheat%isteplat))                    )
      CALL message(TRIM(routine), 'istartlev, iendlev, isteplev: ' &
        & //TRIM(int2string(chemheat%istartlev))//', '             &
        & //TRIM(int2string(chemheat%iendlev))//', '               &
        & //TRIM(int2string(chemheat%isteplev))                    )
    ENDIF

    !----------
    ! Clean-up
    !----------

    DEALLOCATE(varunit, dimunit_lat, dimunit_lev, dimunit_time, STAT=istat)
    IF (istat/=SUCCESS) CALL finish (TRIM(routine), 'Deallocation of characters failed')

  END SUBROUTINE read_extdat_chemheat

  !====================================================================================

  SUBROUTINE get_dim( ncid,    & !in
    &                 dimname, & !in
    &                 dimunit, & !inout
    &                 dimid,   & !inout
    &                 ndim,    & !inout
    &                 varid    ) !inout

    ! In/out variables
    INTEGER,          INTENT(IN)    :: ncid
    CHARACTER(LEN=*), INTENT(IN)    :: dimname
    CHARACTER(LEN=*), INTENT(INOUT) :: dimunit
    INTEGER,          INTENT(INOUT) :: dimid
    INTEGER,          INTENT(INOUT) :: ndim
    INTEGER,          INTENT(INOUT) :: varid

    ! Local variables
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':get_dim'
    
    !-------------------------------------------------------------- 

    ! Initialize out variables
    dimunit = " "
    dimid   = -999
    ndim    = -999
    varid   = -999
    ! Get id of dimension
    CALL nf(nf90_inq_dimid(ncid, dimname, dimid), routine)
    ! Get size of dimension
    CALL nf(nf90_inquire_dimension(ncid, dimid, len = ndim), routine)
    ! Get dimension unit
    CALL nf(nf90_inq_varid(ncid, dimname, varid), routine)
    CALL nf(nf90_get_att(ncid, varid, 'units', dimunit), routine)

  END SUBROUTINE get_dim

  !==================================================================================== 

  !>
  !! Determine auxiliary variables for meridional interpolation.
  !! Please note: we cannot use 
  !! * src/shared/mo_latitude_interpolation: latitude_weights_li
  !! here, because the latitudes from which we interpolate 
  !! (stored in lat_stzstl) are not equidistant in general.
  !!
  SUBROUTINE construct_interpolation_lat( p_patch,    &  !in
    &                                     lat_stzstl, &  !in
    &                                     istart,     &  !in
    &                                     iend,       &  !in
    &                                     istep,      &  !in
    &                                     intrpl_idx, &  !out
    &                                     intrpl_wgt  )  !out

    ! In/out variables
    TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch
    REAL(wp),              INTENT(IN)  :: lat_stzstl(:)       ! (nlat)
    INTEGER,               INTENT(IN)  :: istart, iend, istep
    INTEGER,               INTENT(OUT) :: intrpl_idx(:,:,:)   ! (2,nproma,nblks_c)
    REAL(wp),              INTENT(OUT) :: intrpl_wgt(:,:,:)   ! (2,nproma,nblks_c)

    ! Local variables
    REAL(wp) :: lat, scale

    INTEGER  :: jc, jb, jlat
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx

    LOGICAL  :: lfound

    REAL(wp), PARAMETER :: eps_dlat = 1.e-10_wp

    !--------------------------------------------------------------

    ! Initialization (because of INTENT(OUT)!)
    intrpl_idx(:,:,:) = -999
    intrpl_wgt(:,:,:) = 0._wp

    ! Loop boundaries for prognostic domain.
    rl_start   = grf_bdywidth_c + 1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, i_startidx, i_endidx, jlat, lat, scale, lfound) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      
      DO jc = i_startidx, i_endidx
        
        ! Get meridional position of grid cell
        lat = p_patch%cells%center(jc,jb)%lat

        ! Search latitude among latitudes, 
        ! from which to interpolate
        lfound = .FALSE.
        LAT_LOOP: DO jlat = istart, iend - istep, istep
          IF (lat >= lat_stzstl(jlat) .AND. lat < lat_stzstl(jlat + istep)) THEN
            intrpl_idx(1,jc,jb) = jlat
            intrpl_idx(2,jc,jb) = jlat + istep
            ! Linear interpolation
            scale = 1._wp / MAX(eps_dlat, ABS(lat_stzstl(jlat + istep) - lat_stzstl(jlat)))
            intrpl_wgt(1,jc,jb) = scale * ABS(lat_stzstl(jlat + istep) - lat)
            intrpl_wgt(2,jc,jb) = scale * ABS(lat - lat_stzstl(jlat))       
            lfound = .TRUE.
            EXIT LAT_LOOP
          ENDIF
        ENDDO  LAT_LOOP

        IF (.NOT. lfound) THEN
          IF (lat < lat_stzstl(istart)) THEN
            intrpl_idx(1,jc,jb) = istart
            intrpl_idx(2,jc,jb) = istart
            intrpl_wgt(1,jc,jb) = 0.5_wp
            intrpl_wgt(2,jc,jb) = 0.5_wp
            lfound = .TRUE.            
          ELSEIF (lat >= lat_stzstl(iend)) THEN
            intrpl_idx(1,jc,jb) = iend
            intrpl_idx(2,jc,jb) = iend
            intrpl_wgt(1,jc,jb) = 0.5_wp
            intrpl_wgt(2,jc,jb) = 0.5_wp
            lfound = .TRUE.            
          ENDIF
        ENDIF  !IF (.NOT. lfound)

      ENDDO  !jc

    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE construct_interpolation_lat

  !==================================================================================== 

  !>
  !! Determine auxiliary variables for vertical interpolation.
  !!
  SUBROUTINE construct_interpolation_lev( p_patch,    &  !in
    &                                     lev_stzstl, &  !in
    &                                     istart,     &  !in
    &                                     iend,       &  !in
    &                                     istep,      &  !in
    &                                     intrpl_idx, &  !out
    &                                     intrpl_wgt, &  !out
    &                                     vct_a       )  !(opt)in
    ! In/out variables
    TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch
    REAL(wp),              INTENT(IN)  :: lev_stzstl(:)       ! (nlev_extdat)
    INTEGER,               INTENT(IN)  :: istart, iend, istep
    INTEGER,               INTENT(OUT) :: intrpl_idx(:,:)     ! (2,nlev)
    REAL(wp),              INTENT(OUT) :: intrpl_wgt(:,:)     ! (2,nlev)
    REAL(wp),    OPTIONAL, INTENT(IN)  :: vct_a(:)            ! (nlev+1)

    ! Local variables
    REAL(wp) :: lev, scale

    INTEGER  :: nlev, nshift_total
    INTEGER  :: jk, jks, jlev, istart_dyn, istart_aux

    LOGICAL  :: lfound

    REAL(wp), PARAMETER :: eps_dz = 1.e-10_wp
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':construct_interpolation_lev'
    
    !--------------------------------------------------------------

    IF (.NOT. PRESENT(vct_a)) CALL finish(TRIM(routine), 'vct_a has to be present.')

    ! Initialization (because of INTENT(OUT)!)
    intrpl_idx(:,:) = -999
    intrpl_wgt(:,:) = 0._wp

    ! Number of grid layers on domain
    nlev = p_patch%nlev

    ! Shift of grid layer index, 
    ! to account for vertical nesting
    nshift_total = p_patch%nshift_total

    ! Initialize start level for search
    istart_dyn = istart
    istart_aux = istart

    DO jk = 1, nlev

      jks = jk + nshift_total

      ! Height of grid layer center
      lev = 0.5_wp * ( vct_a(jks) + vct_a(jks + 1) )

      ! Search height among heights, 
      ! from which to interpolate
      lfound = .FALSE.
      LEV_EXT_LOOP: DO jlev = istart_dyn, iend - istep, istep
        IF (lev <= lev_stzstl(jlev) .AND. lev > lev_stzstl(jlev + istep)) THEN
          intrpl_idx(1,jk) = jlev
          intrpl_idx(2,jk) = jlev + istep
          ! Linear interpolation
          scale = 1._wp / MAX(eps_dz, ABS(lev_stzstl(jlev + istep) - lev_stzstl(jlev)))
          intrpl_wgt(1,jk) = scale * ABS(lev_stzstl(jlev + istep) - lev)
          intrpl_wgt(2,jk) = scale * ABS(lev - lev_stzstl(jlev))
          lfound = .TRUE.
          ! For the next search, LEV_EXT_LOOP can start 
          ! at the current position
          istart_aux = jlev
          EXIT LEV_EXT_LOOP
        ENDIF
      ENDDO LEV_EXT_LOOP

      istart_dyn = istart_aux

      IF (.NOT. lfound) THEN
        IF (lev > lev_stzstl(istart)) THEN
          intrpl_idx(1,jk) = istart
          intrpl_idx(2,jk) = istart
          intrpl_wgt(1,jk) = 0.5_wp
          intrpl_wgt(2,jk) = 0.5_wp
          lfound = .TRUE.            
        ELSEIF (lev <= lev_stzstl(iend)) THEN
          intrpl_idx(1,jk) = iend
          intrpl_idx(2,jk) = iend
          intrpl_wgt(1,jk) = 0.5_wp
          intrpl_wgt(2,jk) = 0.5_wp
          lfound = .TRUE.            
        ENDIF
      ENDIF  !IF (.NOT. lfound)

      IF (.NOT. lfound) CALL finish(TRIM(routine), 'Height level '//TRIM(real2string(lev)) &
        & //' not found among stuetzstellen')

    ENDDO  !jk
    
  END SUBROUTINE construct_interpolation_lev

END MODULE mo_upatmo_extdat_utils
