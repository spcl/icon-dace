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
MODULE mo_netcdf_errhandler

  USE mo_exception,          ONLY: finish, message_text, warning
  USE mo_netcdf

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: nf

CONTAINS

  SUBROUTINE nf(errstat, routine, warnonly, silent)

    INTEGER, INTENT(IN)           :: errstat
    CHARACTER(*), INTENT(IN) :: routine
    LOGICAL, INTENT(IN), OPTIONAL :: warnonly, silent
    LOGICAL :: lwarnonly

    IF(PRESENT(silent)) THEN
      IF (silent) RETURN
    END IF
    lwarnonly = .FALSE.
    IF(PRESENT(warnonly)) lwarnonly = .TRUE.
    IF (errstat .NE. nf90_noerr) THEN
      IF (lwarnonly) THEN
        CALL warning(TRIM(routine)//' netCDF error', nf90_strerror(errstat))
      ELSE
        CALL finish(TRIM(routine)//' netCDF error', nf90_strerror(errstat))
      ENDIF
    ENDIF
  END SUBROUTINE nf

END MODULE mo_netcdf_errhandler
