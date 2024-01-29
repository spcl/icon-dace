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

!> This module contains bindings to the functions of C standard library.

MODULE mo_util_libc

  USE ISO_C_BINDING, ONLY: c_ptr, c_char, c_int, c_size_t, C_ASSOCIATED, C_F_POINTER
  USE mo_exception, ONLY: finish
  USE mo_util_string, ONLY: int2string

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: memset_f, memcmp_f, memcpy_f, strerror

  CHARACTER(*), PARAMETER :: modname = "mo_util_libc"
  INTEGER, PARAMETER :: SUCCESS = 0

CONTAINS

  FUNCTION strerror(errorNumber) RESULT(resultVar)
    CHARACTER(:), ALLOCATABLE :: resultVar
    INTEGER, VALUE :: errorNumber

    TYPE(c_ptr) :: c_pointer
    INTEGER :: charPointerShape(1), error, i
    CHARACTER(KIND=c_char), DIMENSION(:), POINTER :: f_pointer
    CHARACTER(*), PARAMETER :: routine = modname//":strerror"

    INTERFACE
      FUNCTION c_strerror(c_error) BIND(C, NAME="strerror") RESULT(c_result)
        IMPORT c_int, c_ptr
        TYPE(c_ptr) :: c_result
        INTEGER(c_int), VALUE :: c_error
      END FUNCTION c_strerror

      INTEGER(c_size_t) FUNCTION c_strlen(charPtr) BIND(C, NAME="strlen")
        IMPORT c_size_t, c_ptr
        TYPE(c_ptr), VALUE :: charPtr
      END FUNCTION c_strlen
    END INTERFACE

    c_pointer = c_strerror(errorNumber)
    IF (C_ASSOCIATED(c_pointer)) THEN
      charPointerShape(1) = INT(c_strlen(c_pointer))
      CALL C_F_POINTER(c_pointer, f_pointer, charPointerShape)
      ALLOCATE (CHARACTER(LEN=charPointerShape(1)) :: resultVar, STAT=error)
      IF (error /= SUCCESS) CALL finish(routine, "memory allocation failure")
      DO i = 1, charPointerShape(1)
        resultVar(i:i) = f_pointer(i)
      END DO
    ELSE
      CALL finish(routine, "strerror("//TRIM(int2string(errorNumber))//") returned NULL")
    END IF
  END FUNCTION strerror

  SUBROUTINE memset_f(str, c, n)
    TYPE(c_ptr), VALUE :: str
    INTEGER(c_int), VALUE :: c
    INTEGER(c_size_t), VALUE :: n

    INTERFACE
      TYPE(c_ptr) FUNCTION c_memset(str, c, n) BIND(C, NAME='memset')
        IMPORT c_ptr, c_int, c_size_t
        TYPE(c_ptr), VALUE :: str
        INTEGER(c_int), VALUE :: c
        INTEGER(c_size_t), VALUE :: n
      END FUNCTION c_memset
    END INTERFACE

    str = c_memset(str, c, n)
  END SUBROUTINE memset_f

  FUNCTION memcmp_f(str1, str2, n) RESULT(f_result)
    TYPE(c_ptr), VALUE, INTENT(IN) :: str1, str2
    INTEGER(c_size_t), VALUE :: n
    LOGICAL :: f_result

    INTERFACE
      INTEGER(c_int) FUNCTION c_memcmp(a, b, n) BIND(C, NAME='memcmp')
        IMPORT c_ptr, c_int, c_size_t
        TYPE(c_ptr), VALUE, INTENT(IN) :: a, b
        INTEGER(c_size_t), VALUE :: n
      END FUNCTION c_memcmp
    END INTERFACE

    f_result = (c_memcmp(str1, str2, n) /= 0)
  END FUNCTION memcmp_f

  TYPE(c_ptr) FUNCTION memcpy_f(dest, src, bsize) RESULT(ret)
    TYPE(c_ptr), VALUE :: dest
    TYPE(c_ptr), INTENT(IN), VALUE :: src
    INTEGER(c_size_t), INTENT(IN) :: bsize

    INTERFACE
      TYPE(c_ptr) FUNCTION c_memcpy(a, b, s) BIND(C, NAME='memcpy')
        IMPORT c_size_t, c_ptr
        TYPE(c_ptr), VALUE :: a
        TYPE(c_ptr), INTENT(IN), VALUE :: b
        INTEGER(c_size_t), INTENT(IN), VALUE :: s
      END FUNCTION c_memcpy
    END INTERFACE

    ret = c_memcpy(dest, src, bsize)
  END FUNCTION memcpy_f

END MODULE mo_util_libc

