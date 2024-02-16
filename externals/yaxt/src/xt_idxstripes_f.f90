!>
!! @file xt_idxstripes_f.f90
!! @brief Fortran interface to yaxt implementation
!!
!! @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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

!>
!! @example test_perf_stripes.f90
MODULE xt_idxstripes
  USE xt_core, ONLY: xt_stripe, i2, i4, i8, xt_abort
  USE xt_idxlist_abstract, ONLY: xt_idxlist, xt_idxlist_c2f, xt_idxlist_f2c
  USE iso_c_binding, ONLY: c_ptr, c_int
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xt_idxstripes_new
  INTERFACE xt_idxstripes_new
    MODULE PROCEDURE xt_idxstripes_new_s
    MODULE PROCEDURE xt_idxstripes_new_a1d
    MODULE PROCEDURE xt_idxstripes_new_a1d_i2
    MODULE PROCEDURE xt_idxstripes_new_a1d_i4
    MODULE PROCEDURE xt_idxstripes_new_a1d_i8
    MODULE PROCEDURE xt_idxstripes_new_a2d
    MODULE PROCEDURE xt_idxstripes_new_a2d_i2
    MODULE PROCEDURE xt_idxstripes_new_a2d_i4
    MODULE PROCEDURE xt_idxstripes_new_a2d_i8
    MODULE PROCEDURE xt_idxstripes_new_a3d
    MODULE PROCEDURE xt_idxstripes_new_a3d_i2
    MODULE PROCEDURE xt_idxstripes_new_a3d_i4
    MODULE PROCEDURE xt_idxstripes_new_a3d_i8
    MODULE PROCEDURE xt_idxstripes_new_a4d
    MODULE PROCEDURE xt_idxstripes_new_a4d_i2
    MODULE PROCEDURE xt_idxstripes_new_a4d_i4
    MODULE PROCEDURE xt_idxstripes_new_a4d_i8
    MODULE PROCEDURE xt_idxstripes_new_a5d
    MODULE PROCEDURE xt_idxstripes_new_a5d_i2
    MODULE PROCEDURE xt_idxstripes_new_a5d_i4
    MODULE PROCEDURE xt_idxstripes_new_a5d_i8
    MODULE PROCEDURE xt_idxstripes_new_a6d
    MODULE PROCEDURE xt_idxstripes_new_a6d_i2
    MODULE PROCEDURE xt_idxstripes_new_a6d_i4
    MODULE PROCEDURE xt_idxstripes_new_a6d_i8
    MODULE PROCEDURE xt_idxstripes_new_a7d
    MODULE PROCEDURE xt_idxstripes_new_a7d_i2
    MODULE PROCEDURE xt_idxstripes_new_a7d_i4
    MODULE PROCEDURE xt_idxstripes_new_a7d_i8
  END INTERFACE xt_idxstripes_new

  INTERFACE
    FUNCTION xt_idxstripes_new_c(stripes, num_stripes) &
         BIND(C, name='xt_idxstripes_new') RESULT(res_ptr)
      IMPORT:: xt_idxlist, xt_stripe, c_ptr, c_int
      IMPLICIT NONE
      TYPE(xt_stripe), INTENT(in) :: stripes(*)
      INTEGER(c_int), VALUE, INTENT(in) :: num_stripes
      TYPE(c_ptr) :: res_ptr
    END FUNCTION xt_idxstripes_new_c
  END INTERFACE

  PUBLIC :: xt_idxstripes_from_idxlist_new

  CHARACTER(len=*), PARAMETER :: filename = 'xt_idxstripes_f.f90'
CONTAINS

  FUNCTION xt_idxstripes_new_s(stripe) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripe
    TYPE(xt_idxlist) :: res

    res = xt_idxlist_c2f(xt_idxstripes_new_c((/ stripe /), 1_c_int))
  END FUNCTION xt_idxstripes_new_s

  FUNCTION xt_idxstripes_new_a1d(stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (SIZE(stripes) > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(SIZE(stripes), c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a1d

  FUNCTION xt_idxstripes_new_a1d_i2(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(*)
    INTEGER(i2), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, &
         INT(num_stripes, c_int)))
  END FUNCTION xt_idxstripes_new_a1d_i2

  FUNCTION xt_idxstripes_new_a1d_i4(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(*)
    INTEGER(i4), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a1d_i4

  FUNCTION xt_idxstripes_new_a1d_i8(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(*)
    INTEGER(i8), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a1d_i8

  FUNCTION xt_idxstripes_new_a2d(stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(:,:)
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (SIZE(stripes) > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(SIZE(stripes), c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a2d

  FUNCTION xt_idxstripes_new_a2d_i2(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, &
         INT(num_stripes, c_int)))
  END FUNCTION xt_idxstripes_new_a2d_i2

  FUNCTION xt_idxstripes_new_a2d_i4(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res
    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a2d_i4

  FUNCTION xt_idxstripes_new_a2d_i8(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a2d_i8

  FUNCTION xt_idxstripes_new_a3d(stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(:,:,:)
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (SIZE(stripes) > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(SIZE(stripes), c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a3d

  FUNCTION xt_idxstripes_new_a3d_i2(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, &
         INT(num_stripes, c_int)))
  END FUNCTION xt_idxstripes_new_a3d_i2

  FUNCTION xt_idxstripes_new_a3d_i4(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a3d_i4

  FUNCTION xt_idxstripes_new_a3d_i8(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a3d_i8

  FUNCTION xt_idxstripes_new_a4d(stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(:,:,:,:)
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (SIZE(stripes) > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(SIZE(stripes), c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a4d

  FUNCTION xt_idxstripes_new_a4d_i2(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, &
         INT(num_stripes, c_int)))
  END FUNCTION xt_idxstripes_new_a4d_i2

  FUNCTION xt_idxstripes_new_a4d_i4(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a4d_i4

  FUNCTION xt_idxstripes_new_a4d_i8(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a4d_i8

  FUNCTION xt_idxstripes_new_a5d(stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(:,:,:,:,:)
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (SIZE(stripes) > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(SIZE(stripes), c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a5d

  FUNCTION xt_idxstripes_new_a5d_i2(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, &
         INT(num_stripes, c_int)))
  END FUNCTION xt_idxstripes_new_a5d_i2

  FUNCTION xt_idxstripes_new_a5d_i4(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a5d_i4

  FUNCTION xt_idxstripes_new_a5d_i8(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a5d_i8

  FUNCTION xt_idxstripes_new_a6d(stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(:,:,:,:,:,:)
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (SIZE(stripes) > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(SIZE(stripes), c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a6d

  FUNCTION xt_idxstripes_new_a6d_i2(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,1,1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, &
         INT(num_stripes, c_int)))
  END FUNCTION xt_idxstripes_new_a6d_i2

  FUNCTION xt_idxstripes_new_a6d_i4(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,1,1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a6d_i4

  FUNCTION xt_idxstripes_new_a6d_i8(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,1,1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a6d_i8

  FUNCTION xt_idxstripes_new_a7d(stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(:,:,:,:,:,:,:)
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (SIZE(stripes) > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(SIZE(stripes), c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a7d

  FUNCTION xt_idxstripes_new_a7d_i2(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,1,1,1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, &
         INT(num_stripes, c_int)))
  END FUNCTION xt_idxstripes_new_a7d_i2

  FUNCTION xt_idxstripes_new_a7d_i4(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,1,1,1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a7d_i4

  FUNCTION xt_idxstripes_new_a7d_i8(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(1,1,1,1,1,1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTEGER(c_int) :: num_stripes_c

    IF (num_stripes > HUGE(num_stripes_c)) &
         CALL xt_abort("too many idxstripes elements", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxstripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxstripes_new_a7d_i8

  FUNCTION xt_idxstripes_from_idxlist_new(idxlist_src) RESULT(idxstripes)
    TYPE(xt_idxlist), INTENT(in) :: idxlist_src
    TYPE(xt_idxlist) :: idxstripes
    INTERFACE
      FUNCTION xt_idxstripes_from_idxlist_new_c(idxlist_src) &
           BIND(c, name='xt_idxstripes_from_idxlist_new') RESULT(idxstripes)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist_src
        TYPE(c_ptr) :: idxstripes
      END FUNCTION xt_idxstripes_from_idxlist_new_c
    END INTERFACE
    idxstripes = xt_idxlist_c2f(&
         xt_idxstripes_from_idxlist_new_c(xt_idxlist_f2c(idxlist_src)))
  END FUNCTION xt_idxstripes_from_idxlist_new

END MODULE xt_idxstripes
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
