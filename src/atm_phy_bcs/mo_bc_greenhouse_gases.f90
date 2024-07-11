!
! Preliminary read and time interpolation of greenhouse gases data
!
! This is  a clone of the respective ECHAM routine
!
! Time series of various greenhouse gases are read from
! file bc_greenhouse_gases.nc (CO2, CH4, N2O, and CFC's).
! Provides interpolation in time and conversion from volume mixing ratio
! to mass mixing ratio of CO2, CH4, and N2O - not for CFC's!
!
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

MODULE mo_bc_greenhouse_gases

  USE mo_kind,               ONLY: wp, dp, i8
  USE mo_exception,          ONLY: finish, message, message_text, warning
  USE mo_physical_constants, ONLY: vmr_to_mmr_co2, vmr_to_mmr_ch4, vmr_to_mmr_n2o, vmr_to_mmr_c11, vmr_to_mmr_c12
  USE mo_netcdf,             ONLY: nf90_noerr, nf90_nowrite, nf90_max_var_dims
  USE mo_netcdf_parallel,    ONLY: p_nf90_open, p_nf90_inq_dimid, p_nf90_inquire_dimension, &
       &                           p_nf90_inq_varid, p_nf90_get_var, p_nf90_close, &
       &                           p_nf90_inquire_variable
  USE mtime,                 ONLY: datetime, no_of_sec_in_a_day, &
       &                           getNoOfDaysInYearDateTime, &
       &                           getdayofyearfromdatetime,  &
       &                           getnoofsecondselapsedindaydatetime

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_bc_greenhouse_gases
  PUBLIC :: bc_greenhouse_gases_time_interpolation
  PUBLIC :: cleanup_greenhouse_gases

  PUBLIC :: ghg_no_cfc

  PUBLIC :: bc_greenhouse_gases_file_read
  PUBLIC :: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcmmr
  PUBLIC :: ghg_co2vmr, ghg_ch4vmr, ghg_n2ovmr, ghg_cfcvmr

  INTEGER, PARAMETER :: ghg_no_cfc = 2
  CHARACTER(len=*), PARAMETER :: ghg_cfc_names(ghg_no_cfc) = (/ "CFC_11", "CFC_12" /)

  REAL(wp) :: ghg_base_year

  INTEGER :: ghg_no_years

  REAL(wp), ALLOCATABLE :: ghg_years(:)
  REAL(wp), ALLOCATABLE :: ghg_co2(:)
  REAL(wp), ALLOCATABLE :: ghg_ch4(:)
  REAL(wp), ALLOCATABLE :: ghg_n2o(:)
  REAL(wp), ALLOCATABLE :: ghg_cfc(:,:)

  REAL(wp), PROTECTED :: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr
  REAL(wp), PROTECTED :: ghg_cfcmmr(ghg_no_cfc)

  REAL(wp), PROTECTED :: ghg_co2vmr, ghg_ch4vmr, ghg_n2ovmr
  REAL(wp), PROTECTED :: ghg_cfcvmr(ghg_no_cfc)

  LOGICAL, SAVE :: bc_greenhouse_gases_file_read = .FALSE.

CONTAINS

  SUBROUTINE read_bc_greenhouse_gases(ghg_filename)

    INTEGER :: ncid, nvarid, time_dimid, var_ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: var_dimids, var_count
    INTEGER :: i
    CHARACTER(LEN=*) :: ghg_filename

    IF (bc_greenhouse_gases_file_read) THEN
      CALL message('','Greenhouse gases already read ...')
      RETURN
    ENDIF

    CALL message('','Use transient, annually resolved greenhouse gases secenario based on CMIP5')
    CALL nf_check(p_nf90_open(ghg_filename, nf90_nowrite, ncid))
    CALL nf_check(p_nf90_inq_dimid (ncid, 'time', time_dimid))
    CALL nf_check(p_nf90_inquire_dimension (ncid, time_dimid, len = ghg_no_years))

    ALLOCATE (ghg_years(ghg_no_years))
    ALLOCATE (ghg_co2(ghg_no_years))
    ALLOCATE (ghg_ch4(ghg_no_years))
    ALLOCATE (ghg_n2o(ghg_no_years))
    ALLOCATE (ghg_cfc(ghg_no_years,ghg_no_cfc))
    !$ACC ENTER DATA PCREATE(ghg_years, ghg_co2, ghg_ch4, ghg_n2o, ghg_cfc, ghg_cfcmmr, ghg_cfcvmr)

    CALL nf_check(p_nf90_inq_varid (ncid, 'time', nvarid))
    CALL nf_check(p_nf90_get_var (ncid, nvarid, ghg_years))
      
    CALL nf_check(p_nf90_inq_varid (ncid, 'CO2', nvarid))
    ! Find the time dimension in the variable and initialize 'var_count' as the
    ! 'count' argument of the 'p_nf90_get_var' function:
    CALL nf_check(p_nf90_inquire_variable(ncid, nvarid, &
                                        & ndims = var_ndims, &
                                        & dimids = var_dimids))
    var_count(:) = 1
    DO i = 1, var_ndims
      IF (var_dimids(i) == time_dimid) THEN
        var_count(i) = ghg_no_years
        EXIT
      ENDIF
    ENDDO
    CALL nf_check(p_nf90_get_var (ncid, nvarid, ghg_co2, count = var_count))
      
    ! Assume that the rest of the variables have the same dimensions:
    CALL nf_check(p_nf90_inq_varid (ncid, 'CH4', nvarid))
    CALL nf_check(p_nf90_get_var (ncid, nvarid, ghg_ch4, count = var_count))
      
    CALL nf_check(p_nf90_inq_varid (ncid, 'N2O', nvarid))
    CALL nf_check(p_nf90_get_var (ncid, nvarid, ghg_n2o, count = var_count))
      
    DO i = 1, ghg_no_cfc
      CALL nf_check(p_nf90_inq_varid (ncid, TRIM(ghg_cfc_names(i)), nvarid))
      CALL nf_check(p_nf90_get_var (ncid, nvarid, ghg_cfc(:,i), count = var_count))
    ENDDO

    bc_greenhouse_gases_file_read = .TRUE.

    CALL nf_check(p_nf90_close(ncid))

    ghg_base_year = ghg_years(1)

    !$ACC UPDATE DEVICE(ghg_years, ghg_co2, ghg_ch4, ghg_n2o, ghg_cfc) ASYNC(1)

  END SUBROUTINE read_bc_greenhouse_gases

  SUBROUTINE bc_greenhouse_gases_time_interpolation(radiation_date, print_report)

    TYPE(datetime), POINTER, INTENT(in) :: radiation_date
    LOGICAL, INTENT(IN), OPTIONAL :: print_report

    REAL(dp) :: zsecref, zsecnow
    REAL(dp) :: zw1, zw2
    REAL(wp) :: zco2int, zch4int, zn2oint
    REAL(wp) :: zcfc(ghg_no_cfc)
    INTEGER(i8) :: yearlen, yearday
    INTEGER :: iyear, iyearm, iyearp

    ! interpolation in time

    yearlen = getNoOfDaysInYearDateTime(radiation_date)*no_of_sec_in_a_day
    yearday = (getdayofyearfromdatetime(radiation_date)-1)*no_of_sec_in_a_day &
         &   +getnoofsecondselapsedindaydatetime(radiation_date)    
    zsecref = REAL(yearlen, dp)
    zsecnow = REAL(yearday, dp)

    iyear  = radiation_date%date%year - INT(ghg_base_year) + 1   ! set right index to access in ghg fields
    iyearm = iyear - 1
    iyearp = iyear + 1

    ! Data are allocated from 1 to ghg_no_years, thus
    ! iyear, iyearm and iyearp shall stay within this
    ! range

    IF (radiation_date%date%month <= 6) THEN     ! first half of year

      IF ( iyear  < 1 .OR. iyear  > ghg_no_years .OR. &
   &       iyearm < 1 .OR. iyearm > ghg_no_years ) THEN

        WRITE (message_text,'(a,i8,a,i8,a,i8)') 'iyear ', iyear,    &
   &                                         ' or iyearm ', iyearm, &
   &                                         ' are out of range 1 - ', ghg_no_years
        CALL finish('mo_bc_greenhouse_gases', message_text)
 
      ENDIF

      zw1 = zsecnow/zsecref + 0.5_dp
      zw2 = 1.0_dp - zw1

      zco2int   = 1.0e-06_wp * ( zw1*ghg_co2(iyear)   + zw2*ghg_co2(iyearm)   )
      zch4int   = 1.0e-09_wp * ( zw1*ghg_ch4(iyear)   + zw2*ghg_ch4(iyearm)   )
      zn2oint   = 1.0e-09_wp * ( zw1*ghg_n2o(iyear)   + zw2*ghg_n2o(iyearm)   )
      zcfc(:)   = 1.0e-12_wp * ( zw1*ghg_cfc(iyear,:) + zw2*ghg_cfc(iyearm,:) )

    ELSE                                         ! second half of year

      IF ( iyear  < 1 .OR. iyear  > ghg_no_years .OR. &
   &       iyearp < 1 .OR. iyearp > ghg_no_years ) THEN

        WRITE (message_text,'(a,i8,a,i8,a,i8)') 'iyear ', iyear,    &
   &                                         ' or iyearp ', iyearp, &
   &                                         ' are out of range 1 - ', ghg_no_years
        CALL finish('mo_bc_greenhouse_gases', message_text)

      ENDIF

      zw2= zsecnow/zsecref - 0.5_dp
      zw1= 1.0_dp - zw2

      zco2int   = 1.0e-06_wp * ( zw1*ghg_co2(iyear)   + zw2*ghg_co2(iyearp)   )
      zch4int   = 1.0e-09_wp * ( zw1*ghg_ch4(iyear)   + zw2*ghg_ch4(iyearp)   )
      zn2oint   = 1.0e-09_wp * ( zw1*ghg_n2o(iyear)   + zw2*ghg_n2o(iyearp)   )
      zcfc(:)   = 1.0e-12_wp * ( zw1*ghg_cfc(iyear,:) + zw2*ghg_cfc(iyearp,:) )
    END IF

    ghg_co2vmr    = zco2int
    ghg_ch4vmr    = zch4int
    ghg_n2ovmr    = zn2oint

    ghg_cfcvmr(1) = zcfc(1)
    ghg_cfcvmr(2) = zcfc(2)

    ghg_co2mmr = ghg_co2vmr * vmr_to_mmr_co2
    ghg_ch4mmr = ghg_ch4vmr * vmr_to_mmr_ch4
    ghg_n2ommr = ghg_n2ovmr * vmr_to_mmr_n2o

    ghg_cfcmmr(1) = ghg_cfcvmr(1) * vmr_to_mmr_c11
    ghg_cfcmmr(2) = ghg_cfcvmr(2) * vmr_to_mmr_c12

    IF (PRESENT(print_report)) THEN
      IF (print_report) THEN
        WRITE (message_text,'(a,5(a,"=",f0.2," ",a,:,", "))') &
            & 'Interpolated VMRs: ', &
            & 'CO2', 1e6_wp * ghg_co2vmr, 'ppm', &
            & 'CH4', 1e9_wp * ghg_ch4vmr, 'ppb', &
            & 'N2O', 1e9_wp * ghg_n2ovmr, 'ppb', &
            & 'CFC11', 1e12_wp * ghg_cfcvmr(1), 'ppt', &
            & 'CFC12', 1e12_wp * ghg_cfcvmr(2), 'ppt'
        CALL message('mo_bc_greenhouse_gases', message_text)
      END IF
    END IF

    !$ACC UPDATE DEVICE(ghg_cfcvmr, ghg_cfcmmr) ASYNC(1)

  END SUBROUTINE bc_greenhouse_gases_time_interpolation

  SUBROUTINE cleanup_greenhouse_gases
    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(ghg_years) IF(ALLOCATED(ghg_years))
    !$ACC EXIT DATA DELETE(ghg_co2) IF(ALLOCATED(ghg_co2))
    !$ACC EXIT DATA DELETE(ghg_ch4) IF(ALLOCATED(ghg_ch4))
    !$ACC EXIT DATA DELETE(ghg_n2o) IF(ALLOCATED(ghg_n2o))
    !$ACC EXIT DATA DELETE(ghg_cfc) IF(ALLOCATED(ghg_cfc))
    !$ACC EXIT DATA DELETE(ghg_cfcmmr)
    IF (ALLOCATED(ghg_years)) DEALLOCATE(ghg_years)
    IF (ALLOCATED(ghg_co2))   DEALLOCATE(ghg_co2)
    IF (ALLOCATED(ghg_ch4))   DEALLOCATE(ghg_ch4)
    IF (ALLOCATED(ghg_n2o))   DEALLOCATE(ghg_n2o)
    IF (ALLOCATED(ghg_cfc))   DEALLOCATE(ghg_cfc)
  END SUBROUTINE cleanup_greenhouse_gases

  SUBROUTINE nf_check(iret)
    USE mo_netcdf_errhandler, ONLY: nf
    INTEGER, INTENT(in) :: iret

    CALL nf( iret, 'mo_bc_greenhouse_gases')
  END SUBROUTINE nf_check

END MODULE mo_bc_greenhouse_gases
