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

! Adapter module for NetCDF. Re-exports NetCDF library functions and adds some
! convenience functions.
MODULE mo_netcdf

USE netcdf

IMPLICIT NONE

PUBLIC

INTERFACE nf90x_put_att_converted
  MODULE PROCEDURE nf90x_put_att_converted_logical
END INTERFACE nf90x_put_att_converted

PRIVATE :: nf90x_put_att_converted_logical

INTERFACE nf90x_get_att_converted
  MODULE PROCEDURE nf90x_get_att_converted_logical
END INTERFACE nf90x_get_att_converted

PRIVATE :: nf90x_get_att_converted_logical

CONTAINS

INTEGER FUNCTION nf90x_put_att_converted_logical(ncid, varid, name, values)
  INTEGER, INTENT(in) :: ncid, varid
  CHARACTER(len=*), INTENT(in) :: name
  logical, INTENT(in) :: values

  INTEGER :: buf

  IF (values) THEN; buf = 1; ELSE; buf = 0; ENDIF

  nf90x_put_att_converted_logical = nf90_put_att(ncid, varid, name, buf)
END FUNCTION nf90x_put_att_converted_logical

INTEGER FUNCTION nf90x_get_att_converted_logical(ncid, varid, name, values)
    INTEGER, INTENT(in) :: ncid, varid
    CHARACTER(len=*), INTENT(in) :: name
    logical, INTENT(out) :: values

    INTEGER :: buf

    nf90x_get_att_converted_logical = nf90_get_att(ncid, varid, name, buf)
    values = (buf /= 0)
END FUNCTION nf90x_get_att_converted_logical

END MODULE mo_netcdf
