!> @file comin_c_utils.F90
!! @brief ComIn C utilities, contains helper routines for converting data
!!        between C and Fortran.
!
!  @authors 01/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_c_utils

  USE iso_c_binding, ONLY: C_CHAR, C_NULL_CHAR
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: convert_c_string, convert_f_string

  INTERFACE convert_c_string
    MODULE PROCEDURE convert_c_string_cptr
    MODULE PROCEDURE convert_c_string_arr
  END INTERFACE convert_c_string

CONTAINS

  !> Convert c-style character array into Fortran string.
  !
  FUNCTION convert_c_string_cptr( cptr ) RESULT (string)

    USE, intrinsic :: iso_c_binding, only: c_ptr, c_char, &
                                                                              c_f_pointer,c_size_t

    TYPE(c_ptr), intent(in) :: cptr
    CHARACTER(len=:), allocatable :: string
    CHARACTER(kind=c_char), dimension(:), pointer :: chars
    INTEGER(kind=c_size_t) :: strlen

    INTERFACE
      FUNCTION c_strlen(str_ptr) BIND ( C, name = "strlen" ) RESULT(len)
        USE, INTRINSIC :: iso_c_binding
        TYPE(c_ptr), VALUE      :: str_ptr
        INTEGER(kind=c_size_t)  :: len
      END FUNCTION c_strlen
    END INTERFACE

    strlen = c_strlen(cptr)
    CALL c_f_pointer(cptr, chars, [ strlen ])
    string = convert_c_string_arr(chars)
  END FUNCTION convert_c_string_cptr

  FUNCTION convert_c_string_arr( arr ) RESULT (string)
    CHARACTER(kind=C_CHAR), INTENT(IN) :: arr(:)
    CHARACTER(len=:), ALLOCATABLE :: string
    INTEGER :: i, strlen

    DO strlen=1,SIZE(arr)
      IF (arr(strlen) .EQ. c_null_char) EXIT
    END DO

    ALLOCATE(CHARACTER(len=strlen-1) :: string)
    DO i=1,strlen-1
      string(i:i) = arr(i)
    END DO
  END FUNCTION convert_c_string_arr

  !> Convert Fortran string into C-style character array.
  !
  subroutine convert_f_string( string, arr)

    CHARACTER(len=*), INTENT(IN) :: string
    CHARACTER(len=1, kind=c_char), INTENT(INOUT) :: arr(:)
    INTEGER :: i

    DO i=1,LEN_TRIM(string)
      arr(i) = string(i:i)
    END DO
    arr(i) = c_null_char
  end subroutine  convert_f_string

END MODULE comin_c_utils
