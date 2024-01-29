!> Adapter module for NetCDF.
!! Re-exports NetCDF library functions and adds some convenience
!! functions.
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
MODULE mo_netcdf

  USE netcdf4_nf_interfaces
  USE netcdf_nf_interfaces
  USE netcdf_nf_data

  IMPLICIT NONE

  ! Public default, to reexport NetCDF routines and variables
  PUBLIC

  PUBLIC :: nfx_put_att
  PUBLIC :: nfx_get_att

  !> Convenience function to put scalar or array attributes into a NetCDF file.
  !! Logicals are converted to integers 1 (`.TRUE.`) or 0 (`.FALSE.`).
  INTERFACE nfx_put_att
    MODULE PROCEDURE nfx_put_att_i2_0
    MODULE PROCEDURE nfx_put_att_i2_1
    MODULE PROCEDURE nfx_put_att_i4_0
    MODULE PROCEDURE nfx_put_att_i4_1
    MODULE PROCEDURE nfx_put_att_r4_0
    MODULE PROCEDURE nfx_put_att_r4_1
    MODULE PROCEDURE nfx_put_att_r8_0
    MODULE PROCEDURE nfx_put_att_r8_1
    MODULE PROCEDURE nfx_put_att_l_0
    MODULE PROCEDURE nfx_put_att_l_1
  END INTERFACE

  PRIVATE :: nfx_put_att_i2_0, nfx_put_att_i2_1
  PRIVATE :: nfx_put_att_i4_0, nfx_put_att_i4_1
  PRIVATE :: nfx_put_att_r4_0, nfx_put_att_r4_1
  PRIVATE :: nfx_put_att_r8_0, nfx_put_att_r8_1
  PRIVATE :: nfx_put_att_l_0, nfx_put_att_l_1

  !> Convenience function to get scalar or array attributes from a NetCDF file.
  !! Logicals are converted from integers, with nonzero values considered `.TRUE.`.
  INTERFACE nfx_get_att
    MODULE PROCEDURE nfx_get_att_i2_0
    MODULE PROCEDURE nfx_get_att_i2_1
    MODULE PROCEDURE nfx_get_att_i4_0
    MODULE PROCEDURE nfx_get_att_i4_1
    MODULE PROCEDURE nfx_get_att_r4_0
    MODULE PROCEDURE nfx_get_att_r4_1
    MODULE PROCEDURE nfx_get_att_r8_0
    MODULE PROCEDURE nfx_get_att_r8_1
    MODULE PROCEDURE nfx_get_att_l_0
    MODULE PROCEDURE nfx_get_att_l_1
  END INTERFACE

  PRIVATE :: nfx_get_att_i2_0, nfx_get_att_i2_1
  PRIVATE :: nfx_get_att_i4_0, nfx_get_att_i4_1
  PRIVATE :: nfx_get_att_r4_0, nfx_get_att_r4_1
  PRIVATE :: nfx_get_att_r8_0, nfx_get_att_r8_1
  PRIVATE :: nfx_get_att_l_0, nfx_get_att_l_1

CONTAINS

Function nfx_put_att_i2_0 (ncid, varid, name, xtype, i2vals) &
  RESULT(status)

  Integer,          Intent(IN) :: ncid, varid, xtype
  Character(LEN=*), Intent(IN) :: name
  Integer(NFINT2),   Intent(IN) :: i2vals
  Integer                      :: status

  Integer(NFINT2) :: ival(1)

    ival(1) = i2vals
    status = nf_put_att_int2(ncid, varid, name, xtype, 1, ival)

End Function

Function nfx_put_att_i2_1 (ncid, varid, name, xtype, i2vals) &
  RESULT(status)

  Integer,          Intent(IN) :: ncid, varid, xtype
  Character(LEN=*), Intent(IN) :: name
  Integer(NFINT2), Intent(IN), Contiguous :: i2vals(:)
  Integer                      :: status

  status = nf_put_att_int2(ncid, varid, name, xtype, SIZE(i2vals), i2vals)

End Function


Function nfx_put_att_i4_0 (ncid, varid, name, xtype, ivals) &
    RESULT(status)

    Integer,          Intent(IN) :: ncid, varid, xtype
    Character(LEN=*), Intent(IN) :: name
    Integer(NFINT),   Intent(IN) :: ivals
    Integer                      :: status

    Integer(NFINT) :: ival(1)

      ival(1) = ivals
      status = nf_put_att_int(ncid, varid, name, xtype, 1, ival)

  End Function

  Function nfx_put_att_i4_1 (ncid, varid, name, xtype, ivals) &
    RESULT(status)

    Integer,          Intent(IN) :: ncid, varid, xtype
    Character(LEN=*), Intent(IN) :: name
    Integer(NFINT), Intent(IN), Contiguous :: ivals(:)
    Integer                      :: status

    status = nf_put_att_int(ncid, varid, name, xtype, SIZE(ivals), ivals)

  End Function


  Function nfx_put_att_r4_0 (ncid, varid, name, xtype, rvals) &
    RESULT(status)

    Integer,          Intent(IN) :: ncid, varid, xtype
    Character(LEN=*), Intent(IN) :: name
    Real(RK4),        Intent(IN) :: rvals
    Integer                      :: status

    Real(RK4) :: dval(1)

    dval(1) = rvals
    status = nf_put_att_real(ncid, varid, name, xtype, 1, dval)

  End Function

  Function nfx_put_att_r4_1 (ncid, varid, name, xtype, rvals) &
    RESULT(status)

    Integer,          Intent(IN) :: ncid, varid, xtype
    Character(LEN=*), Intent(IN) :: name
    Real(RK4), Intent(IN), Contiguous :: rvals(:)
    Integer                      :: status

    status = nf_put_att_real(ncid, varid, name, xtype, SIZE(rvals), rvals)

  End Function


  Function nfx_put_att_r8_0 (ncid, varid, name, xtype, dvals) &
    RESULT(status)

    Integer,          Intent(IN) :: ncid, varid, xtype
    Character(LEN=*), Intent(IN) :: name
    Real(RK8),        Intent(IN) :: dvals
    Integer                      :: status

    Real(RK8) :: dval(1)

    dval(1) = dvals
    status = nf_put_att_double(ncid, varid, name, xtype, 1, dval)

  End Function

  Function nfx_put_att_r8_1 (ncid, varid, name, xtype, dvals) &
    RESULT(status)

    Integer,          Intent(IN) :: ncid, varid, xtype
    Character(LEN=*), Intent(IN) :: name
    Real(RK8), Intent(IN), Contiguous :: dvals(:)
    Integer                      :: status

    status = nf_put_att_double(ncid, varid, name, xtype, SIZE(dvals), dvals)

  End Function


  Function nfx_put_att_l_0 (ncid, varid, name, xtype, lvals) &
    RESULT(status)

    Integer,          Intent(IN) :: ncid, varid, xtype
    Character(LEN=*), Intent(IN) :: name
    Logical,          Intent(IN) :: lvals
    Integer                      :: status

    Integer(NFINT) :: lval(1)

    lval(1) = MERGE(1, 0, lvals)
    status = nf_put_att_int(ncid, varid, name, xtype, 1, lval)

  End Function

  Function nfx_put_att_l_1 (ncid, varid, name, xtype, lvals) &
    RESULT(status)

    Integer,          Intent(IN) :: ncid, varid, xtype
    Character(LEN=*), Intent(IN) :: name
    Logical, Intent(IN), Contiguous :: lvals(:)
    Integer                      :: status

    Integer(NFINT) :: ivals(SIZE(lvals))

    ivals(:) = MERGE(1, 0, lvals)

    status = nf_put_att_int(ncid, varid, name, xtype, SIZE(ivals), ivals)

  End Function


  Function nfx_get_att_i2_0(ncid, varid, name, i2vals) RESULT(status)

    Integer,          Intent(IN)  :: ncid, varid
    Character(LEN=*), Intent(IN)  :: name
    Integer(NFINT2),  Intent(OUT) :: i2vals
    Integer                       :: status

    Integer(NFINT2) :: ival(1)

    status = nf_get_att_int2(ncid, varid, name, ival)
    i2vals = ival(1)

  End Function

  Function nfx_get_att_i2_1(ncid, varid, name, i2vals) RESULT(status)

    Integer,          Intent(IN)  :: ncid, varid
    Character(LEN=*), Intent(IN)  :: name
    Integer(NFINT2),  Intent(OUT) :: i2vals(:)
    Integer                       :: status

    status = nf_get_att_int2(ncid, varid, name, i2vals)

  End Function


  Function nfx_get_att_i4_0(ncid, varid, name, ivals) RESULT(status)

    Integer,          Intent(IN)  :: ncid, varid
    Character(LEN=*), Intent(IN)  :: name
    Integer(NFINT),   Intent(OUT) :: ivals
    Integer                       :: status

    Integer(NFINT) :: ival(1)

    status = nf_get_att_int(ncid, varid, name, ival)
    ivals = ival(1)

  End Function

  Function nfx_get_att_i4_1(ncid, varid, name, ivals) RESULT(status)

    Integer,          Intent(IN)  :: ncid, varid
    Character(LEN=*), Intent(IN)  :: name
    Integer(NFINT),   Intent(OUT) :: ivals(:)
    Integer                       :: status

    status = nf_get_att_int(ncid, varid, name, ivals)

  End Function


  Function nfx_get_att_r4_0(ncid, varid, name, rvals) RESULT(status)

    Integer,          Intent(IN)  :: ncid, varid
    Character(LEN=*), Intent(IN)  :: name
    Real(RK4),        Intent(OUT) :: rvals
    Integer                       :: status

    Real(RK4) :: dval(1)

    status = nf_get_att_real(ncid, varid, name, dval)
    rvals = dval(1)

  End Function

  Function nfx_get_att_r4_1(ncid, varid, name, rvals) RESULT(status)

    Integer,          Intent(IN)  :: ncid, varid
    Character(LEN=*), Intent(IN)  :: name
    Real(RK4),        Intent(OUT) :: rvals(:)
    Integer                       :: status

    status = nf_get_att_real(ncid, varid, name, rvals)

  End Function


  Function nfx_get_att_r8_0(ncid, varid, name, dvals) RESULT(status)

    Integer,          Intent(IN)  :: ncid, varid
    Character(LEN=*), Intent(IN)  :: name
    Real(RK8),        Intent(OUT) :: dvals
    Integer                       :: status

    Real(RK8) :: dval(1)

    status = nf_get_att_double(ncid, varid, name, dval)
    dvals = dval(1)

  End Function

  Function nfx_get_att_r8_1(ncid, varid, name, dvals) RESULT(status)

    Integer,          Intent(IN)  :: ncid, varid
    Character(LEN=*), Intent(IN)  :: name
    Real(RK8),        Intent(OUT) :: dvals(:)
    Integer                       :: status

    status = nf_get_att_double(ncid, varid, name, dvals)

  End Function


  Function nfx_get_att_l_0(ncid, varid, name, lvals) RESULT(status)

    Integer,          Intent(IN)  :: ncid, varid
    Character(LEN=*), Intent(IN)  :: name
    Logical,          Intent(OUT) :: lvals
    Integer                       :: status

    Integer(NFINT) :: ival(1)

    status = nf_get_att_int(ncid, varid, name, ival)
    lvals = (ival(1) /= 0)

  End Function

  Function nfx_get_att_l_1(ncid, varid, name, lvals) RESULT(status)

    Integer,          Intent(IN)  :: ncid, varid
    Character(LEN=*), Intent(IN)  :: name
    Logical,          Intent(OUT) :: lvals(:)
    Integer                       :: status

    Integer(NFINT) :: ivals(SIZE(lvals))

    status = nf_get_att_int(ncid, varid, name, ivals)
    lvals(:) = (ivals /= 0)

  End Function

END MODULE mo_netcdf
