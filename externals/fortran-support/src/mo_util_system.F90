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

MODULE mo_util_system

  USE ISO_C_BINDING, ONLY: c_int, c_char, c_null_char

  IMPLICIT NONE

  PRIVATE

#ifdef __XT3__
  PUBLIC :: util_base_iobuf
#endif
  PUBLIC :: util_exit
  PUBLIC :: util_abort
  PUBLIC :: util_system

  INTERFACE
#ifdef __XT3__
    SUBROUTINE util_base_iobuf() BIND(C)
    END SUBROUTINE util_base_iobuf
#endif
    SUBROUTINE util_exit(exit_no) BIND(C)
      IMPORT c_int
      INTEGER(c_int), VALUE :: exit_no
    END SUBROUTINE util_exit
    SUBROUTINE util_abort() BIND(C)
    END SUBROUTINE util_abort
  END INTERFACE

CONTAINS

  FUNCTION util_system(f_s) RESULT(f_result)
    INTEGER(c_int) :: f_result
    CHARACTER(KIND=c_char, LEN=*), INTENT(IN) :: f_s

    CHARACTER(KIND=c_char) :: c_s(LEN(f_s) + 1)
    INTEGER :: i

    INTERFACE
      FUNCTION c_util_system(c_s) BIND(C, NAME='util_system') RESULT(c_result)
        IMPORT c_int, c_char
        INTEGER(c_int) :: c_result
        CHARACTER(KIND=c_char), INTENT(IN) :: c_s(*)
      END FUNCTION c_util_system
    END INTERFACE
    DO i = 1, LEN(f_s)
      c_s(i) = f_s(i:i)
    END DO
    c_s(LEN(f_s) + 1) = c_null_char
    f_result = c_util_system(c_s)
  END FUNCTION util_system

END MODULE mo_util_system

