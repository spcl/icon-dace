!>
!! @file test_misc_f.f90
!! @brief Test minor functionality needed for higher level parts of YAXT.
!!
!! @copyright Copyright  (C)  2021 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
#include "fc_feature_defs.inc"
PROGRAM test_misc
  USE xt_core, ONLY: i2, i4, i8
  USE ftest_common, ONLY: icbrt, test_abort, run_randomized_tests, &
       init_fortran_random
  IMPLICIT NONE
  LOGICAL :: fully_random_tests
  CHARACTER(len=*), PARAMETER :: filename = 'test_misc_f.f90'

  fully_random_tests = run_randomized_tests()
  CALL init_fortran_random(fully_random_tests)
  CALL test_icbrt
CONTAINS
  SUBROUTINE test_icbrt
    CALL test_icbrt_i2
    CALL test_icbrt_i4
    CALL test_icbrt_i8
  END SUBROUTINE test_icbrt

  SUBROUTINE test_icbrt_i2
    INTEGER(i2), PARAMETER :: ulim = 31_i2
    INTEGER(i2) :: i, cubed, cbrt
    CHARACTER(len=132) :: msg
    DO i = -ulim, ulim, 1_i2
      cubed = i**3_i2
      cbrt = icbrt(cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", i, &
             "**3 = ", cubed, ", but icbrt(", cubed, ") = ", cbrt
        CALL test_abort(msg, filename, __LINE__)
      END IF
    END DO
  END SUBROUTINE test_icbrt_i2

  SUBROUTINE test_icbrt_i4
    INTEGER(i4), PARAMETER :: ulim = 1290
    INTEGER(i4) :: i, cubed, cbrt, prev_cubed, other_cubed
    CHARACTER(len=132) :: msg
    DOUBLE PRECISION :: rnd
    i = icbrt(-8)
    prev_cubed = -HUGE(1_i4)
    DO i = -ulim, -1_i4, 1_i4
      cubed = i**3_i4
      cbrt = icbrt(cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", i, &
             "**3 = ", cubed, ", but icbrt(", cubed, ") = ", cbrt
        CALL test_abort(msg, filename, __LINE__)
      END IF
      CALL random_number(rnd)
      other_cubed = prev_cubed + MAX(1_i4, INT(rnd * DBLE(cubed - prev_cubed), i4))
      cbrt = icbrt(other_cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", other_cubed, &
             ", expected ", i, &
             ", but icbrt(", other_cubed, ") = ", cbrt
        CALL test_abort(msg, filename, __LINE__)
      END IF
      prev_cubed = cubed
    END DO
    prev_cubed = HUGE(1_i4)
    DO i = ulim, 0_i4, -1_i4
      cubed = i**3_i4
      cbrt = icbrt(cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", i, &
             "**3 = ", cubed, ", but icbrt(", cubed, ") = ", cbrt
        CALL test_abort(msg, filename, __LINE__)
      END IF
      CALL random_number(rnd)
      other_cubed = cubed + MIN(INT(rnd * DBLE(prev_cubed - cubed), i4), -1_i4)
      cbrt = icbrt(cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", other_cubed, &
             ", expected ", i, &
             ", but icbrt(", other_cubed, ") = ", cbrt
        CALL test_abort(msg, filename, __LINE__)
      END IF
      prev_cubed = cubed
    END DO
  END SUBROUTINE test_icbrt_i4

  SUBROUTINE test_icbrt_i8
    INTEGER(i8), PARAMETER :: ulim = 2097151_i8
    INTEGER(i8) :: i, cubed, cbrt
    CHARACTER(len=132) :: msg
    DO i = -ulim, ulim, 1_i8
      ! need to spell that out explcitely for crayftn 8.5.5 at least
#ifndef HAVE_FC_PRECISE_INTEGER_EXPONENTIATION
      cubed = i * i * i
#else
      cubed = i**3
#endif
      cbrt = icbrt(cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", i, &
             "**3 = ", cubed, ", but icbrt(", cubed, ") = ", cbrt
        CALL test_abort(msg, filename, __LINE__)
      END IF
    END DO
  END SUBROUTINE test_icbrt_i8
END PROGRAM test_misc
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
! license-default: "bsd"
! End:
!
