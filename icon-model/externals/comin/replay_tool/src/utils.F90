!> @file utils.F90
!! @brief Utility function to make programming with FORTRAN bearable
!
!  @authors 01/2024 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.

MODULE utils

  IMPLICIT NONE

  PUBLIC :: int2string

CONTAINS

  FUNCTION int2string(n, opt_fmt)
    CHARACTER(:), ALLOCATABLE :: int2string ! result
    CHARACTER(len=128) :: res
    INTEGER, INTENT(in) :: n
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=128) :: fmt

    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(i0)'
    END IF
    WRITE(res,fmt) n
    res = ADJUSTL(res)
    int2string = TRIM(res)
  END FUNCTION int2string

END MODULE utils
