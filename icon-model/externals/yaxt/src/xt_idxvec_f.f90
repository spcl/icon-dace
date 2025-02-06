!>
!! @file xt_idxvec_f.f90
!! @brief Fortran interface to yaxt xt_idxvec functions
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
MODULE xt_idxvec
  USE xt_core, ONLY: i2, i4, i8, xt_int_kind, xt_abort, xt_stripe
  USE xt_idxlist_abstract, ONLY: xt_idxlist, xt_idxlist_c2f
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_int
  IMPLICIT NONE
  PRIVATE
  INTERFACE xt_idxvec_new
    MODULE PROCEDURE xt_idxvec_new_a1d
    MODULE PROCEDURE xt_idxvec_new_a1d_i2
    MODULE PROCEDURE xt_idxvec_new_a1d_i4
    MODULE PROCEDURE xt_idxvec_new_a1d_i8
    MODULE PROCEDURE xt_idxvec_new_a2d
    MODULE PROCEDURE xt_idxvec_new_a2d_i2
    MODULE PROCEDURE xt_idxvec_new_a2d_i4
    MODULE PROCEDURE xt_idxvec_new_a2d_i8
    MODULE PROCEDURE xt_idxvec_new_a3d
    MODULE PROCEDURE xt_idxvec_new_a3d_i2
    MODULE PROCEDURE xt_idxvec_new_a3d_i4
    MODULE PROCEDURE xt_idxvec_new_a3d_i8
    MODULE PROCEDURE xt_idxvec_new_a4d
    MODULE PROCEDURE xt_idxvec_new_a4d_i2
    MODULE PROCEDURE xt_idxvec_new_a4d_i4
    MODULE PROCEDURE xt_idxvec_new_a4d_i8
    MODULE PROCEDURE xt_idxvec_new_a5d
    MODULE PROCEDURE xt_idxvec_new_a5d_i2
    MODULE PROCEDURE xt_idxvec_new_a5d_i4
    MODULE PROCEDURE xt_idxvec_new_a5d_i8
    MODULE PROCEDURE xt_idxvec_new_a6d
    MODULE PROCEDURE xt_idxvec_new_a6d_i2
    MODULE PROCEDURE xt_idxvec_new_a6d_i4
    MODULE PROCEDURE xt_idxvec_new_a6d_i8
    MODULE PROCEDURE xt_idxvec_new_a7d
    MODULE PROCEDURE xt_idxvec_new_a7d_i2
    MODULE PROCEDURE xt_idxvec_new_a7d_i4
    MODULE PROCEDURE xt_idxvec_new_a7d_i8
  END INTERFACE xt_idxvec_new

  INTERFACE xt_idxvec_from_stripes_new
    MODULE PROCEDURE xt_idxvec_from_stripes_new_a
    MODULE PROCEDURE xt_idxvec_from_stripes_new_a_i2
    MODULE PROCEDURE xt_idxvec_from_stripes_new_a_i4
    MODULE PROCEDURE xt_idxvec_from_stripes_new_a_i8
  END INTERFACE xt_idxvec_from_stripes_new

  PUBLIC :: xt_idxvec_from_stripes_new, xt_idxvec_new

  INTERFACE
    FUNCTION xt_idxvec_new_c(idxvec, num_indices) &
         BIND(C, name='xt_idxvec_new') RESULT(res_ptr)
      IMPORT :: xt_int_kind, c_ptr, c_int
      IMPLICIT NONE
      INTEGER(xt_int_kind), INTENT(in) :: idxvec(*)
      INTEGER(c_int), VALUE, INTENT(in) :: num_indices
      TYPE(c_ptr) :: res_ptr
    END FUNCTION xt_idxvec_new_c

    FUNCTION xt_idxvec_from_stripes_new_c(stripes, num_stripes) &
         BIND(C, name='xt_idxvec_from_stripes_new') RESULT(res_ptr)
      IMPORT :: xt_stripe, c_int, c_ptr
      IMPLICIT NONE
      TYPE(xt_stripe), INTENT(in) :: stripes(*)
      INTEGER(c_int), VALUE, INTENT(in) :: num_stripes
      TYPE(c_ptr) :: res_ptr
    END FUNCTION xt_idxvec_from_stripes_new_c
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: filename = 'xt_idxvec_f.f90'
CONTAINS

  FUNCTION xt_idxvec_new_a1d(idxvec) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(:)
    TYPE(Xt_idxlist) :: res

    INTEGER(xt_int_kind) :: idxvec_dummy(1)
    INTEGER(c_int) :: num_indices_c
    IF (SIZE(idxvec) > HUGE(num_indices_c)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(SIZE(idxvec), c_int)
    IF (num_indices_c > 0_c_int) THEN
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
    ELSE
      idxvec_dummy(1) = HUGE(idxvec_dummy)
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec_dummy, num_indices_c))
    END IF
  END FUNCTION xt_idxvec_new_a1d

  FUNCTION xt_idxvec_new_a1d_i2(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(*)
    INTEGER(i2), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a1d_i2

  FUNCTION xt_idxvec_new_a1d_i4(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(*)
    INTEGER(i4), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i4 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a1d_i4

  FUNCTION xt_idxvec_new_a1d_i8(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(*)
    INTEGER(i8), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i8 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a1d_i8

  FUNCTION xt_idxvec_new_a2d(idxvec) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(:,:)
    TYPE(Xt_idxlist) :: res
    INTEGER(xt_int_kind) :: idxvec_dummy(1)
    INTEGER(c_int) :: num_indices_c
    IF (SIZE(idxvec) > HUGE(num_indices_c)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(SIZE(idxvec), c_int)
    IF (num_indices_c > 0_c_int) THEN
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
    ELSE
      idxvec_dummy(1) = HUGE(idxvec_dummy)
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec_dummy, num_indices_c))
    END IF
  END FUNCTION xt_idxvec_new_a2d

  FUNCTION xt_idxvec_new_a2d_i2(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a2d_i2

  FUNCTION xt_idxvec_new_a2d_i4(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i4 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a2d_i4

  FUNCTION xt_idxvec_new_a2d_i8(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i8 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a2d_i8

  FUNCTION xt_idxvec_new_a3d(idxvec) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(:,:,:)
    TYPE(Xt_idxlist) :: res

    INTEGER(xt_int_kind) :: idxvec_dummy(1)
    INTEGER(c_int) :: num_indices_c
    IF (SIZE(idxvec) > HUGE(num_indices_c)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(SIZE(idxvec), c_int)
    IF (num_indices_c > 0_c_int) THEN
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
    ELSE
      idxvec_dummy(1) = HUGE(idxvec_dummy)
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec_dummy, num_indices_c))
    END IF
  END FUNCTION xt_idxvec_new_a3d

  FUNCTION xt_idxvec_new_a3d_i2(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a3d_i2

  FUNCTION xt_idxvec_new_a3d_i4(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c
    IF (i4 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a3d_i4

  FUNCTION xt_idxvec_new_a3d_i8(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c
    IF (i8 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a3d_i8

  FUNCTION xt_idxvec_new_a4d(idxvec) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(:,:,:,:)
    TYPE(Xt_idxlist) :: res

    INTEGER(xt_int_kind) :: idxvec_dummy(1)
    INTEGER(c_int) :: num_indices
    IF (SIZE(idxvec) > HUGE(num_indices)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices = INT(SIZE(idxvec), c_int)
    IF (num_indices > 0) THEN
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices))
    ELSE
      idxvec_dummy(1) = HUGE(idxvec_dummy)
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec_dummy, num_indices))
    END IF
  END FUNCTION xt_idxvec_new_a4d

  FUNCTION xt_idxvec_new_a4d_i2(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a4d_i2

  FUNCTION xt_idxvec_new_a4d_i4(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i4 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a4d_i4

  FUNCTION xt_idxvec_new_a4d_i8(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i8 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a4d_i8

  FUNCTION xt_idxvec_new_a5d(idxvec) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(:,:,:,:,:)
    TYPE(Xt_idxlist) :: res

    INTEGER(xt_int_kind) :: idxvec_dummy(1)
    INTEGER(c_int) :: num_indices_c
    IF (SIZE(idxvec) > HUGE(num_indices_c)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(SIZE(idxvec), c_int)
    IF (num_indices_c > 0_c_int) THEN
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
    ELSE
      idxvec_dummy(1) = HUGE(idxvec_dummy)
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec_dummy, num_indices_c))
    END IF
  END FUNCTION xt_idxvec_new_a5d

  FUNCTION xt_idxvec_new_a5d_i2(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))

  END FUNCTION xt_idxvec_new_a5d_i2

  FUNCTION xt_idxvec_new_a5d_i4(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i4 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a5d_i4

  FUNCTION xt_idxvec_new_a5d_i8(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i8 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a5d_i8

  FUNCTION xt_idxvec_new_a6d(idxvec) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(:,:,:,:,:,:)
    TYPE(Xt_idxlist) :: res

    INTEGER(xt_int_kind) :: idxvec_dummy(1)
    INTEGER(c_int) :: num_indices_c

    IF (SIZE(idxvec) > HUGE(num_indices_c)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(SIZE(idxvec), c_int)
    IF (num_indices_c > 0_c_int) THEN
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
    ELSE
      idxvec_dummy(1) = HUGE(idxvec_dummy)
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec_dummy, num_indices_c))
    END IF
  END FUNCTION xt_idxvec_new_a6d

  FUNCTION xt_idxvec_new_a6d_i2(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,1,1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_indices

    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a6d_i2

  FUNCTION xt_idxvec_new_a6d_i4(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,1,1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i4 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a6d_i4

  FUNCTION xt_idxvec_new_a6d_i8(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,1,1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i8 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a6d_i8

  FUNCTION xt_idxvec_new_a7d(idxvec) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(:,:,:,:,:,:,:)
    TYPE(Xt_idxlist) :: res

    INTEGER(xt_int_kind) :: idxvec_dummy(1)
    INTEGER(c_int) :: num_indices_c
    IF (SIZE(idxvec) > HUGE(num_indices_c)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(SIZE(idxvec), c_int)
    IF (num_indices_c > 0_c_int) THEN
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
    ELSE
      idxvec_dummy(1) = HUGE(idxvec_dummy)
      res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec_dummy, num_indices_c))
    END IF
  END FUNCTION xt_idxvec_new_a7d

  FUNCTION xt_idxvec_new_a7d_i2(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,1,1,1,*)
    INTEGER(i2), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a7d_i2

  FUNCTION xt_idxvec_new_a7d_i4(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,1,1,1,*)
    INTEGER(i4), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i4 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a7d_i4

  FUNCTION xt_idxvec_new_a7d_i8(idxvec, num_indices) RESULT(res)
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(1,1,1,1,1,1,*)
    INTEGER(i8), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_indices_c

    IF (i8 /= c_int .AND. num_indices > HUGE(1_c_int)) &
         CALL xt_abort("too many idxvec elements", filename, __LINE__)
    num_indices_c = INT(num_indices, c_int)
    res = xt_idxlist_c2f(xt_idxvec_new_c(idxvec, num_indices_c))
  END FUNCTION xt_idxvec_new_a7d_i8

  FUNCTION xt_idxvec_from_stripes_new_a(stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_stripes_c
    num_stripes_c = INT(SIZE(stripes), c_int)
    res = xt_idxlist_c2f(xt_idxvec_from_stripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxvec_from_stripes_new_a

  FUNCTION xt_idxvec_from_stripes_new_a_i2(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(*)
    INTEGER(i2), INTENT(in) :: num_stripes
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_stripes_c
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxvec_from_stripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxvec_from_stripes_new_a_i2

  FUNCTION xt_idxvec_from_stripes_new_a_i4(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(*)
    INTEGER(i4), INTENT(in) :: num_stripes
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_stripes_c

    IF (i4 /= c_int .AND. num_stripes > HUGE(1_c_int)) &
         CALL xt_abort("too many stripes", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxvec_from_stripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxvec_from_stripes_new_a_i4

  FUNCTION xt_idxvec_from_stripes_new_a_i8(stripes, num_stripes) RESULT(res)
    TYPE(xt_stripe), INTENT(in) :: stripes(*)
    INTEGER(i8), INTENT(in) :: num_stripes
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_stripes_c

    IF (i8 /= c_int .AND. num_stripes > HUGE(1_c_int)) &
         CALL xt_abort("too many stripes", filename, __LINE__)
    num_stripes_c = INT(num_stripes, c_int)
    res = xt_idxlist_c2f(xt_idxvec_from_stripes_new_c(stripes, num_stripes_c))
  END FUNCTION xt_idxvec_from_stripes_new_a_i8

END MODULE xt_idxvec
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
