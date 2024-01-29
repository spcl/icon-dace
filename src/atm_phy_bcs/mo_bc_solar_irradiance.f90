! This is  a clone of the respective ECHAM routine
!
! Read spectrally resolved solar irradiance yearly, apply primitive time
! interpolation to monthly mean values and apply
!
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

MODULE mo_bc_solar_irradiance

  USE mo_kind,            ONLY: dp, i8
  USE mo_exception,       ONLY: finish, message, warning, message_text
  USE mo_netcdf_parallel, ONLY: p_nf_open, p_nf_inq_dimid, p_nf_inq_dimlen, &
       &                        p_nf_inq_varid, p_nf_get_vara_double, p_nf_close, &
       &                        nf_read, nf_noerr, p_nf_get_var_int
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights
  USE mo_run_config,      ONLY: msg_level


  IMPLICIT NONE
  PRIVATE

  REAL(dp), POINTER :: tsi_radt_m(:) => NULL(), tsi_m(:) => NULL()
  REAL(dp), POINTER :: ssi_radt_m(:,:) => NULL()

  INTEGER, ALLOCATABLE :: ssi_years(:)
  INTEGER, ALLOCATABLE :: ssi_months(:)

  PUBLIC :: read_bc_solar_irradiance, ssi_time_interpolation

  LOGICAL, SAVE :: lread_solar = .TRUE., lread_solar_radt = .TRUE.
  INTEGER(i8), SAVE :: last_year = -999999, last_year_radt = -999999

CONTAINS

  SUBROUTINE read_bc_solar_irradiance(year, lradt)
    INTEGER(i8), INTENT(in) :: year
    LOGICAL, INTENT(in) :: lradt ! lradt=.true.: read data for radiation time step
                                 ! the radiation time step may be in the future
                                 ! lradt=.false.: read data for heating rates only
                                 ! data needed at every integration time step

    INTEGER :: ncid, ndimid, nvarid
    INTEGER :: ssi_time_entries
    INTEGER :: ssi_numwl
    INTEGER :: start(2), cnt(2)

    INTEGER :: idx, first_year

    IF (lradt) THEN
       IF (last_year_radt /= year ) lread_solar_radt=.TRUE.
    ELSE
       IF (last_year /= year) lread_solar=.TRUE.
    ENDIF

    IF ((lradt .AND. .NOT. lread_solar_radt) .OR. (.NOT. lradt .AND. .NOT. lread_solar)) THEN
      RETURN
    ENDIF

    WRITE (message_text,'(A,I4,A)') 'reading bc_solar_irradiance_sw_b14.nc (year ', year, ')'
    CALL message('mo_bc_solar_irradiance:read_bc_solar_irradiance', message_text)

    CALL nf_check(p_nf_open('bc_solar_irradiance_sw_b14.nc', nf_read, ncid))

    CALL nf_check(p_nf_inq_dimid(ncid, 'time', ndimid))
    CALL nf_check(p_nf_inq_dimlen(ncid, ndimid, ssi_time_entries))
    CALL nf_check(p_nf_inq_dimid(ncid, 'numwl', ndimid))
    CALL nf_check(p_nf_inq_dimlen(ncid, ndimid, ssi_numwl))

    ALLOCATE (ssi_years(ssi_time_entries))
    ALLOCATE (ssi_months(ssi_time_entries))

    IF (lradt) THEN
       IF (.NOT.(ASSOCIATED(tsi_radt_m))) ALLOCATE(tsi_radt_m(0:13))
       IF (.NOT.(ASSOCIATED(ssi_radt_m))) ALLOCATE(ssi_radt_m(ssi_numwl,0:13))
       !$ACC ENTER DATA PCREATE(tsi_radt_m, ssi_radt_m)
    ELSE
       IF (.NOT.(ASSOCIATED(tsi_m)))      ALLOCATE(tsi_m(0:13))
       !$ACC ENTER DATA PCREATE(tsi_m)
    END IF

    CALL nf_check(p_nf_inq_varid(ncid, 'year', nvarid))
    CALL nf_check(p_nf_get_var_int (ncid, nvarid, ssi_years))
    CALL nf_check(p_nf_inq_varid(ncid, 'month', nvarid))
    CALL nf_check(p_nf_get_var_int (ncid, nvarid, ssi_months))

    first_year = ssi_years(1)

    ! not adding 1 in calculating the offset leads to an index to December of year-1
    idx = 12*INT(year - INT(first_year,i8))
    IF (idx < 1) THEN
      CALL finish('','No solar irradiance data available for the requested year')
    END IF

    CALL nf_check(p_nf_inq_varid (ncid, 'TSI', nvarid))
    start(1) = idx
    cnt(1) = 14
    IF (lradt) THEN
       CALL nf_check(p_nf_get_vara_double(ncid, nvarid, start, cnt, tsi_radt_m))
       CALL nf_check(p_nf_inq_varid (ncid, 'SSI', nvarid))
       start(1) = 1;   cnt(1) = ssi_numwl;
       start(2) = idx; cnt(2) = 14;
       CALL nf_check(p_nf_get_vara_double (ncid, nvarid, start, cnt, ssi_radt_m))
       lread_solar_radt=.FALSE.
       last_year_radt=year
       !$ACC UPDATE DEVICE(tsi_radt_m, ssi_radt_m) ASYNC(1)
    ELSE
       CALL nf_check(p_nf_get_vara_double(ncid, nvarid, start, cnt, tsi_m))
       lread_solar=.FALSE.
       last_year=year
       !$ACC UPDATE DEVICE(tsi_m) ASYNC(1)
    END IF

    CALL nf_check(p_nf_close(ncid))

    DEALLOCATE(ssi_years)
    DEALLOCATE(ssi_months)

  END SUBROUTINE read_bc_solar_irradiance

  SUBROUTINE ssi_time_interpolation(tiw, lradt, tsi, ssi)
    TYPE( t_time_interpolation_weights), INTENT(in) :: tiw
    LOGICAL, INTENT(in)             :: lradt
    REAL(dp), INTENT(out)           :: tsi
    REAL(dp), INTENT(out), OPTIONAL :: ssi(:)
    CHARACTER(len=14)               :: ctsi

    IF (lradt) THEN
      IF (.NOT. (tiw%year1 == last_year_radt .OR. tiw%year2 == last_year_radt)) THEN
        WRITE (message_text,'(A,I4,A,I4,A,I4)') 'Stale data: requested years are ', tiw%year1, ' and ', &
            & tiw%year2, ' but data is for ', last_year_radt
        CALL finish('mo_bc_solar_irradiance:ssi_time_interpolation', message_text)
      END IF
      IF (.NOT.PRESENT(ssi)) THEN
        CALL finish ('ssi_time_interpolation of mo_bc_solar_irradiance', &
                     'Interpolation to radiation time step needs ssi')
      ELSE
        tsi    = tiw%weight1 * tsi_radt_m(tiw%month1_index) + tiw%weight2 * tsi_radt_m(tiw%month2_index)
        ssi(:) = tiw%weight1*ssi_radt_m(:,tiw%month1_index) + tiw%weight2*ssi_radt_m(:,tiw%month2_index)
        WRITE(ctsi,'(F14.8)') tsi
        IF (msg_level >= 11) CALL message('','Interpolated total solar irradiance and spectral ' &
          &          //'bands for radiation transfer, tsi= '//ctsi)
        !$ACC ENTER DATA PCREATE(ssi)
        !$ACC UPDATE DEVICE(ssi) ASYNC(1)
      END IF
    ELSE
      IF (PRESENT(ssi)) THEN
        IF (msg_level >= 11) CALL message ('ssi_time_interpolation of mo_bc_solar_irradiance', &
                     'Interpolation of ssi not necessary')
      END IF
      tsi    = tiw%weight1 * tsi_m(tiw%month1_index) + tiw%weight2 * tsi_m(tiw%month2_index)
    END IF

  END SUBROUTINE ssi_time_interpolation


  SUBROUTINE nf_check(iret)
    USE mo_netcdf_errhandler, ONLY: nf
    INTEGER, INTENT(in) :: iret

    CALL nf(iret, 'mo_bc_solar_irradiance')
  END SUBROUTINE nf_check

END MODULE mo_bc_solar_irradiance
