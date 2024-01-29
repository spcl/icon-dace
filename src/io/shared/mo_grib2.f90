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

MODULE mo_grib2

  USE mo_kind,                  ONLY: dp

  IMPLICIT NONE

  PRIVATE

  ! max. number of additional GRIB2 integer keys per variable
  INTEGER, PARAMETER :: MAX_INT_KEYS = 10

  ! max. number of additional GRIB2 double keys per variable
  INTEGER, PARAMETER :: MAX_DBL_KEYS = 5


  TYPE t_grib2_global
    INTEGER :: centre
    INTEGER :: subcentre
    INTEGER :: generating_process
  END TYPE t_grib2_global

  TYPE t_grib2_int_key
    CHARACTER(len=50) :: key
    INTEGER           :: val
  END TYPE t_grib2_int_key

  TYPE t_grib2_dbl_key
    CHARACTER(len=50) :: key
    REAL(dp)          :: val
  END TYPE t_grib2_dbl_key

  TYPE t_grib2_key_list
    INTEGER                :: nint_keys  ! no. of integer keys
    INTEGER                :: ndbl_keys  ! no. of double keys
    TYPE (t_grib2_int_key) :: int_key(MAX_INT_KEYS)    
    TYPE (t_grib2_dbl_key) :: dbl_key(MAX_DBL_KEYS)
  END TYPE t_grib2_key_list

  TYPE t_grib2_var
    INTEGER :: discipline
    INTEGER :: category
    INTEGER :: number
    INTEGER :: bits
    INTEGER :: gridtype
    INTEGER :: subgridtype
    
    ! list of additional GRIB2 key/value pairs
    TYPE (t_grib2_key_list) :: additional_keys
  END TYPE t_grib2_var



  INTERFACE OPERATOR(+)
    MODULE PROCEDURE grib2_key_list_plus_int
    MODULE PROCEDURE grib2_key_list_plus_dbl
  END INTERFACE OPERATOR(+)

  PUBLIC :: t_grib2_global
  PUBLIC :: t_grib2_var
  PUBLIC :: t_grib2_int_key
  PUBLIC :: t_grib2_dbl_key
  PUBLIC :: t_grib2_key_list
  PUBLIC :: OPERATOR(+)

  ! constructor
  PUBLIC :: grib2_var


CONTAINS


  ! constructor for GRIB2 derived data type
  FUNCTION grib2_var(discipline, category, number, bits, gridtype, subgridtype)
    INTEGER, INTENT(IN) :: discipline
    INTEGER, INTENT(IN) :: category
    INTEGER, INTENT(IN) :: number
    INTEGER, INTENT(IN) :: bits
    INTEGER, INTENT(IN) :: gridtype
    INTEGER, INTENT(IN) :: subgridtype
    TYPE(t_grib2_var) :: grib2_var

    grib2_var%discipline  = discipline 
    grib2_var%category    = category   
    grib2_var%number      = number     
    grib2_var%bits        = bits       
    grib2_var%gridtype    = gridtype   
    grib2_var%subgridtype = subgridtype

    grib2_var%additional_keys%nint_keys = 0
    grib2_var%additional_keys%ndbl_keys = 0
  END FUNCTION grib2_var

  FUNCTION grib2_key_list_plus_int(a, b)
    TYPE(t_grib2_var) :: grib2_key_list_plus_int
    TYPE(t_grib2_var),     INTENT(IN) :: a
    TYPE(t_grib2_int_key), INTENT(IN) :: b
    ! local variables
    INTEGER :: i

    grib2_key_list_plus_int = a
    IF (a%additional_keys%nint_keys < MAX_INT_KEYS) THEN
      i = grib2_key_list_plus_int%additional_keys%nint_keys
      grib2_key_list_plus_int%additional_keys%nint_keys = i + 1
      grib2_key_list_plus_int%additional_keys%int_key(i+1) = b
    END IF
  END FUNCTION grib2_key_list_plus_int

  FUNCTION grib2_key_list_plus_dbl(a, b)
    TYPE(t_grib2_var) :: grib2_key_list_plus_dbl
    TYPE(t_grib2_var),     INTENT(IN) :: a
    TYPE(t_grib2_dbl_key), INTENT(IN) :: b
    ! local variables
    INTEGER :: i

    grib2_key_list_plus_dbl = a
    IF (a%additional_keys%ndbl_keys < MAX_DBL_KEYS) THEN
      i = grib2_key_list_plus_dbl%additional_keys%ndbl_keys
      grib2_key_list_plus_dbl%additional_keys%ndbl_keys = i + 1
      grib2_key_list_plus_dbl%additional_keys%dbl_key(i+1) = b
    END IF
  END FUNCTION grib2_key_list_plus_dbl

END MODULE mo_grib2
