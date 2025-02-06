!> Contains basic definitions and functions for I/O
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_io
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp, dp
  USE mo_jsb_io_iface, ONLY:                                           &
    Get_netcdf_precision,                                              &
    DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24,  &
    ZAXIS_SURFACE, ZAXIS_DEPTH_BELOW_LAND, ZAXIS_GENERIC,              &
    FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2,           &
    GRID_CELL, GRID_UNSTRUCTURED, GRID_UNSTRUCTURED_CELL,              &
    GRID_LONLAT, TSTEP_CONSTANT, TSTEP_INSTANT,                        &
    cdiDefMissval, read_jsb_io_namelist

  USE mo_exception, ONLY: finish

  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: grib_bits     = DATATYPE_PACK16, &
                        grib_extbits  = DATATYPE_PACK24, &
                        nc_bits       = DATATYPE_FLT32,  &
                        nc_extbits    = DATATYPE_FLT64

  TYPE t_cf
    CHARACTER(len=128) :: standard_name = ''
    CHARACTER(len=128) :: units         = ''
    CHARACTER(len=128) :: long_name     = ''
    INTEGER            :: datatype      = -1
  END TYPE t_cf

  TYPE t_grib1
    INTEGER :: table
    INTEGER :: parameter
    INTEGER :: bits
  END TYPE t_grib1

  TYPE t_grib2
    INTEGER :: discipline
    INTEGER :: category
    INTEGER :: number
    INTEGER :: bits
  END TYPE t_grib2

  ! Derived-type constructor for t_cf
  INTERFACE t_cf
    PROCEDURE Construct_cf
  END INTERFACE t_cf

  INTEGER, PARAMETER :: max_tables = 5

  INTEGER :: tables(max_tables) = (/186, 187, 188, 189, 190/)

  ! Remark on usage of "cdiDefMissval" (taken from icon:mo_initicon.f90)

  ! Inside the GRIB_API (v.1.9.18) the missing value is converted into
  ! LONG INT for a test, but the default CDI missing value is outside
  ! of the valid range for LONG INT (U. Schulzweida, bug report
  ! SUP-277). This causes a crash with INVALID OPERATION.

  ! As a workaround we can choose a different missing value in the
  ! calling subroutine (here). For the SX-9 this must lie within 53
  ! bits, because "the SX compiler generates codes using HW
  ! instructions for floating-point data instead of instructions for
  ! integers. Therefore, the operation result is not guaranteed if the
  ! value cannot be represented as integer within 53 bits."
  REAL(dp), PARAMETER :: missval = -9.E+15_dp

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io'

CONTAINS

  FUNCTION Construct_cf(standard_name, units, long_name, datatype) RESULT(return_value)

    CHARACTER(len=*), INTENT(in) :: standard_name
    CHARACTER(len=*), INTENT(in) :: units
    CHARACTER(len=*), INTENT(in) :: long_name
    INTEGER, OPTIONAL,  INTENT(in) :: datatype
    TYPE(t_cf)         :: return_value

    return_value%standard_name = TRIM(standard_name)
    return_value%units         = TRIM(units)
    return_value%long_name     = TRIM(long_name)
    IF (PRESENT(datatype)) THEN
      return_value%datatype = datatype
    ELSE
      return_value%datatype = Get_netcdf_precision()
    END IF

  END FUNCTION Construct_cf

  !>
  !> Initialize io control for JSBACH
  !!
  SUBROUTINE init_jsb_io(namelist_filename)

    CHARACTER(len=*), INTENT(in) :: namelist_filename

    CALL read_jsb_io_namelist(TRIM(namelist_filename))

  END SUBROUTINE init_jsb_io


#endif
END MODULE mo_jsb_io
