!>
!! @file test_idxlist_utils_f.f90
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
MODULE test_idxlist_utils
  USE yaxt, ONLY: xt_int_kind, xt_idxlist, xt_idxlist_c2f, xt_idxlist_f2c, &
       xt_stripe
  USE ftest_common, ONLY: test_abort
  USE iso_c_binding, ONLY: c_ptr, c_int, c_size_t
  IMPLICIT NONE
  PRIVATE
  INTERFACE
    FUNCTION test_err_count() BIND(c, name='test_err_count') RESULT(code)
      IMPORT :: c_int
      INTEGER(c_int) :: code
    END FUNCTION test_err_count
  END INTERFACE
  PUBLIC :: test_err_count
  PUBLIC :: check_idxlist, check_stripes, check_offsets, &
       idxlist_pack_unpack_copy, check_idxlist_copy

  CHARACTER(len=*), PARAMETER :: filename = 'test_idxlist_utils_f.f90'
CONTAINS
  SUBROUTINE check_idxlist(idxlist, ref_indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)

    INTEGER :: num_ref_indices
    INTEGER(xt_int_kind) :: dummy(1)
    INTEGER(c_int) :: num_ref_indices_c

    INTERFACE
      SUBROUTINE check_idxlist_c(idxlist, ref_indices, ref_num_indices) &
           BIND(c, name='check_idxlist')
        IMPORT :: xt_int_kind, c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(xt_int_kind), INTENT(in) :: ref_indices(*)
        INTEGER(c_int), VALUE, INTENT(in) :: ref_num_indices
      END SUBROUTINE check_idxlist_c
    END INTERFACE

    num_ref_indices = SIZE(ref_indices)
    IF (num_ref_indices > 0) THEN
      num_ref_indices_c = INT(num_ref_indices, c_int)
      CALL check_idxlist_c(xt_idxlist_f2c(idxlist), ref_indices, &
           num_ref_indices_c)
    ELSE
      dummy(1) = 0_xt_int_kind
      CALL check_idxlist_c(xt_idxlist_f2c(idxlist), dummy, &
           0_c_int)
    END IF
  END SUBROUTINE check_idxlist

  SUBROUTINE check_stripes(stripes, ref_stripes)
    TYPE(xt_stripe), INTENT(in) :: stripes(:), ref_stripes(:)
    INTEGER(c_int) :: num_stripes_c, ref_num_stripes_c
    INTERFACE
      SUBROUTINE check_stripes_c(stripes, num_stripes, ref_stripes, &
           ref_num_stripes) BIND(c, name='check_stripes')
        IMPORT :: xt_stripe, c_int
        IMPLICIT NONE
        TYPE(xt_stripe), INTENT(in) :: stripes(*), ref_stripes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: num_stripes, ref_num_stripes
      END SUBROUTINE check_stripes_c
    END INTERFACE

    num_stripes_c = INT(SIZE(stripes), c_int)
    ref_num_stripes_c = INT(SIZE(ref_stripes), c_int)
    CALL check_stripes_c(stripes, num_stripes_c, ref_stripes, ref_num_stripes_c)

  END SUBROUTINE check_stripes

  SUBROUTINE check_offsets(offsets_a, offsets_b)
    INTEGER(c_int), INTENT(in) :: offsets_a(:), offsets_b(:)
    INTEGER(c_size_t) :: num_offsets_c
    INTERFACE
      SUBROUTINE check_offsets_c(num_offsets, offsets_a, offsets_b) &
           BIND(c, name='check_offsets')
        IMPORT :: c_size_t, c_int
        IMPLICIT NONE
        INTEGER(c_size_t), VALUE, INTENT(in) :: num_offsets
        INTEGER(c_int), INTENT(IN) :: offsets_a(num_offsets), &
             offsets_b(num_offsets)
      END SUBROUTINE check_offsets_c
    END INTERFACE

    IF (SIZE(offsets_a) /= SIZE(offsets_b)) &
         CALL test_abort("inequal number of array elements in eq test", &
         filename, __LINE__)

    num_offsets_c = INT(SIZE(offsets_a), c_size_t)
    CALL check_offsets_c(num_offsets_c, offsets_a, offsets_b)

  END SUBROUTINE check_offsets

  FUNCTION idxlist_pack_unpack_copy(idxlist) RESULT(idxlist_copy)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    TYPE(xt_idxlist) :: idxlist_copy

    INTERFACE
      FUNCTION idxlist_pack_unpack_copy_c(idxlist) RESULT(idxlist_copy) &
           BIND(c, name='idxlist_pack_unpack_copy')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        TYPE(c_ptr) :: idxlist_copy
      END FUNCTION idxlist_pack_unpack_copy_c
    END INTERFACE

    idxlist_copy &
         = xt_idxlist_c2f(idxlist_pack_unpack_copy_c(xt_idxlist_f2c(idxlist)))

  END FUNCTION idxlist_pack_unpack_copy

  SUBROUTINE check_idxlist_copy(idxlist, idxlist_copy, ref_indices, ref_stripes)
    TYPE(xt_idxlist), INTENT(in) :: idxlist, idxlist_copy
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)
    TYPE(xt_stripe), INTENT(in) :: ref_stripes(:)
    INTEGER(c_size_t) :: num_ref_indices_c, num_ref_stripes_c
    INTERFACE
      SUBROUTINE check_idxlist_copy_c(idxlist, idxlist_copy, &
           num_ref_indices, ref_indices, &
           num_ref_stripes, ref_stripes) BIND(c, name='check_idxlist_copy')
        IMPORT :: c_size_t, c_ptr, xt_int_kind, xt_stripe
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist, idxlist_copy
        INTEGER(c_size_t), VALUE, INTENT(in) :: num_ref_indices, num_ref_stripes
        INTEGER(xt_int_kind), INTENT(in) :: ref_indices(*)
        TYPE(xt_stripe), INTENT(in) :: ref_stripes(*)
      END SUBROUTINE check_idxlist_copy_c
    END INTERFACE
    num_ref_indices_c = INT(SIZE(ref_indices), c_size_t)
    num_ref_stripes_c = INT(SIZE(ref_stripes), c_size_t)
    CALL check_idxlist_copy_c(xt_idxlist_f2c(idxlist), &
         xt_idxlist_f2c(idxlist_copy), num_ref_indices_c, ref_indices, &
         num_ref_stripes_c, ref_stripes)
  END SUBROUTINE check_idxlist_copy

END MODULE test_idxlist_utils
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
