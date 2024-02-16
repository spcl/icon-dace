!>
!! @file test_idxstripes_f.f90
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
PROGRAM test_idxstripes_f
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, &
       run_randomized_tests, init_fortran_random
  USE mpi
  USE test_idxlist_utils, ONLY: check_idxlist, test_err_count, &
       idxlist_pack_unpack_copy
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_bounds, xt_pos_ext, &
       xt_stripe, xt_idxlist, xt_idxlist_delete, xt_idxstripes_new, &
       xt_idxvec_from_stripes_new, xt_int_kind, xt_idxlist_copy, &
       xt_idxvec_from_stripes_new, &
       xt_idxlist_get_index_stripes, xt_idxlist_get_intersection, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_bounding_box, OPERATOR(/=), &
       xt_idxlist_get_pos_exts_of_index_stripes, &
       xt_idxlist_get_num_indices, xt_idxvec_new, &
       xt_idxstripes_from_idxlist_new
  USE iso_c_binding, ONLY: c_int
  IMPLICIT NONE
  INTEGER, PARAMETER :: xi = xt_int_kind
  LOGICAL :: fully_random_tests
  CHARACTER(len=*), PARAMETER :: filename = 'test_idxstripes_f.f90'

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  CALL stripe_test_general1
  CALL stripe_test_general2
  CALL stripe_test_general3
  CALL stripe_test_general4
  CALL stripe_test_general5
  CALL stripe_test_general6
  CALL test_intersection0
  CALL test_intersection1
  CALL test_intersection2
  CALL test_intersection3
  CALL test_intersection4
  CALL test_intersection5
  CALL test_intersection6
  CALL test_intersection7
  CALL test_intersection8
  CALL test_intersection9
  CALL test_intersection10
  CALL test_intersection11
  CALL test_intersection12
  CALL test_intersection13
  CALL test_intersection14
  CALL test_intersection15
  CALL test_intersection_stripe2vec
  CALL test_idxlist_stripes_pos_ext1
  CALL test_idxlist_stripes_pos_ext2
  CALL test_idxlist_stripes_pos_ext3
#if SIZEOF_XT_INT > 2
  CALL test_idxlist_stripes_pos_ext4
  CALL test_idxlist_stripes_pos_ext5
#endif
  CALL test_idxlist_stripes_pos_ext_randomized1(.FALSE.)
  fully_random_tests = run_randomized_tests()
  IF (fully_random_tests) &
       CALL test_idxlist_stripes_pos_ext_randomized1(.TRUE.)
  CALL test_get_pos1
  CALL test_get_pos2
  CALL test_get_pos3
  CALL test_get_pos4
  CALL test_stripe_overlap
  CALL test_stripe_bb1
  CALL test_stripe_bb2
  CALL check_pos_ext1
  CALL check_pos_ext2
  CALL check_pos_ext3
  CALL check_pos_ext4
  CALL check_pos_ext5
  CALL check_pos_ext6
  CALL check_pos_ext7
  CALL check_pos_ext8
  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_finalize
  CALL finish_mpi

CONTAINS
  SUBROUTINE stripe_test_general(stripes, ref_indices)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)

    TYPE(xt_idxlist) :: idxstripes, idxvec
    INTEGER :: num_ext, num_unmatched, num_pos, i
    INTEGER(c_int) :: ext_size
    TYPE(xt_pos_ext), ALLOCATABLE :: pos_ext(:)

    idxstripes = xt_idxstripes_new(stripes, SIZE(stripes))
    CALL do_tests(idxstripes, ref_indices)

    num_unmatched = xt_idxlist_get_pos_exts_of_index_stripes(idxstripes, &
         stripes, pos_ext, .TRUE.)
    IF (num_unmatched /= 0) &
         CALL test_abort("stripes not found", filename, __LINE__)

    num_pos = 0
    IF (ALLOCATED(pos_ext)) THEN
      num_ext = SIZE(pos_ext)
    ELSE
      num_ext = 0
    END IF
    DO i = 1, num_ext
      ext_size = pos_ext(i)%size
      IF (num_pos /= pos_ext(i)%start) &
           CALL test_abort("position/start mismatch", filename, __LINE__)
      num_pos = num_pos + ext_size
    END DO
    IF (num_pos /= xt_idxlist_get_num_indices(idxstripes)) &
         CALL test_abort("index list length/positions overlap mismatch", &
         filename, __LINE__)

    IF (ALLOCATED(pos_ext)) DEALLOCATE(pos_ext)
    CALL xt_idxlist_delete(idxstripes)

    ! test recreation of stripes from reference vector
    idxvec = xt_idxvec_new(ref_indices)
    idxstripes = xt_idxstripes_from_idxlist_new(idxvec)
    CALL check_idxlist(idxstripes, ref_indices)
    CALL xt_idxlist_delete(idxvec)
    CALL xt_idxlist_delete(idxstripes)
  END SUBROUTINE stripe_test_general

  SUBROUTINE stripe_test_general1
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 1, 5), xt_stripe(20, 1, 5) /);
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0_xi,  1_xi,  2_xi,  3_xi,  4_xi, &
         &   10_xi, 11_xi, 12_xi, 13_xi, 14_xi, &
         &   20_xi, 21_xi, 22_xi, 23_xi, 24_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general1

  SUBROUTINE stripe_test_general2
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 2, 5), xt_stripe(20, 3, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0_xi,  1_xi,  2_xi,  3_xi,  4_xi, &
         &   10_xi, 12_xi, 14_xi, 16_xi, 18_xi, &
         &   20_xi, 23_xi, 26_xi, 29_xi, 32_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general2

  SUBROUTINE stripe_test_general3
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(0, 6, 5), &
         xt_stripe(1, 3, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, 6_xi, 12_xi, 18_xi, 24_xi, &
         &    1_xi, 4_xi,  7_xi, 10_xi, 13_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general3

  SUBROUTINE stripe_test_general4
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(0, -1, 5), &
         xt_stripe(1, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, -1_xi, -2_xi, -3_xi, -4_xi, &
         &    1_xi,  2_xi,  3_xi,  4_xi,  5_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general4

  SUBROUTINE stripe_test_general5
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(9, -2, 5), &
         xt_stripe(0, 2, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 9_xi,  7_xi,  5_xi,  3_xi,  1_xi, &
         &    0_xi,  2_xi,  4_xi,  6_xi,  8_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general5

  SUBROUTINE stripe_test_general6
    TYPE(xt_stripe), PARAMETER :: stripes(1) = (/ xt_stripe(179, -2, 0) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(1) = (/ 0_xi /)
    CALL stripe_test_general(stripes, ref_indices(1:0))
  END SUBROUTINE stripe_test_general6

  SUBROUTINE test_intersection(stripes_a, stripes_b, ref_indices_a, ref_indices_b)
    TYPE(xt_stripe), INTENT(in) :: stripes_a(:), stripes_b(:)
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices_a(:)
    INTEGER(xt_int_kind), OPTIONAL, INTENT(in) :: ref_indices_b(:)
    TYPE(xt_idxlist) :: idxstripes_a, idxstripes_b, intersection(2)

    idxstripes_a = xt_idxstripes_new(stripes_a)
    idxstripes_b = xt_idxstripes_new(stripes_b)
    intersection(1) = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b)
    intersection(2) = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a)
    CALL do_tests(intersection(1), ref_indices_a)
    IF (PRESENT(ref_indices_b)) THEN
      CALL do_tests(intersection(2), ref_indices_b)
    ELSE
      CALL do_tests(intersection(2), ref_indices_a)
    END IF
    CALL xt_idxlist_delete(intersection(2))
    CALL xt_idxlist_delete(intersection(1))
    CALL xt_idxlist_delete(idxstripes_a)
    CALL xt_idxlist_delete(idxstripes_b)
  END SUBROUTINE test_intersection

  SUBROUTINE test_intersection0
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 0, 0) /), &
         stripes_b(1) = (/ xt_stripe(1, 1, 8) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(1) = (/ 0_xi /)
    CALL test_intersection(stripes_a(1:0), stripes_b, ref_indices(1:0))
  END SUBROUTINE test_intersection0

  SUBROUTINE test_intersection1
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 1, 4), &
         xt_stripe(6, 1, 4) /), &
         stripes_b(1) = (/ xt_stripe(1, 1, 8) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(6) &
         = (/ 1_xi, 2_xi, 3_xi, 6_xi, 7_xi, 8_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection1

  SUBROUTINE test_intersection2
    TYPE(xt_stripe), PARAMETER :: stripes_a(3) = (/ xt_stripe(0, 1, 4), &
         xt_stripe(6, 1, 4), xt_stripe(11, 1, 4) /), &
         stripes_b(2) = (/ xt_stripe(1, 1, 7), xt_stripe(9, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(9) &
         = (/ 1_xi, 2_xi, 3_xi, 6_xi, 7_xi, 9_xi, 11_xi, 12_xi, 13_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection2

  SUBROUTINE test_intersection3
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 1, 3), &
         xt_stripe(8, 1, 3) /), &
         stripes_b(2) = (/ xt_stripe(3, 1, 5), xt_stripe(11, 1, 3) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(1) = (/ -1_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices(1:0))
  END SUBROUTINE test_intersection3

  SUBROUTINE test_intersection4
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 1, 10) /), &
         stripes_b(2) = (/ xt_stripe(0, 2, 5), xt_stripe(9, -2, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, 8_xi, 9_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection4

  SUBROUTINE test_intersection5
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 3, 5), &
         xt_stripe(1, 7, 5) /), &
         stripes_b(2) = (/ xt_stripe(0, 2, 7), xt_stripe(24, -1, 10) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(6) &
         = (/ 0_xi, 6_xi, 8_xi, 12_xi, 15_xi, 22_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection5

  SUBROUTINE test_intersection6
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 1, 10) /), &
         stripes_b(2) = (/ xt_stripe(5, 1, 5), xt_stripe(4, -1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, 8_xi, 9_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection6

  SUBROUTINE test_intersection7
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 1, 10) , &
            xt_stripe(20, 1, 5) /), &
         stripes_b(2) = (/ xt_stripe(3, 1, 5), xt_stripe(17, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(7) &
         = (/ 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, 20_xi, 21_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection7

  SUBROUTINE test_intersection8
    TYPE(xt_stripe), PARAMETER :: stripes_a(10) = (/ xt_stripe(0, 1, 2), &
         xt_stripe(3, 1, 2), xt_stripe(5, 1, 2), xt_stripe(8, 1, 2), &
         xt_stripe(10, 1, 2), xt_stripe(14, 1, 2), xt_stripe(17, 1, 2), &
         xt_stripe(20, 1, 2), xt_stripe(23, 1, 2), xt_stripe(25, 1, 2) /), &
         stripes_b(5) = (/ xt_stripe(5, 1, 3), xt_stripe(8, 1, 2), &
         xt_stripe(19, 1, 1), xt_stripe(20, 1, 2), xt_stripe(30, 1, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(6) &
         = (/ 5_xi, 6_xi, 8_xi, 9_xi, 20_xi, 21_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection8

  SUBROUTINE test_intersection9
    TYPE(xt_stripe), PARAMETER :: stripes_a(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(1, 1, 5), xt_stripe(2, 1, 5) /), &
         stripes_b(1) = (/ xt_stripe(-2, 1, 10) /)
#ifndef __G95__
    INTEGER(xi) :: i
    INTEGER(xt_int_kind), PARAMETER :: ref_indices_a(7) &
         = (/ (i, i=0_xi,6_xi) /), &
#else
    INTEGER :: i
    INTEGER(xt_int_kind), PARAMETER :: ref_indices_a(7) &
         = (/ (INT(i, xi), i=0_xi,6_xi) /), &
#endif
         ref_indices_b(15) = (/ 0_xi, 1_xi, 1_xi, 2_xi, 2_xi, 2_xi, 3_xi, &
         &                      3_xi, 3_xi, 4_xi, 4_xi, 4_xi, 5_xi, 5_xi, 6_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices_a, ref_indices_b)
  END SUBROUTINE test_intersection9

  SUBROUTINE test_intersection10
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 2, 5) /), &
         stripes_b(1) = (/ xt_stripe(1, 2, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: dummy(1) = (/ -1_xi /)
    CALL test_intersection(stripes_a, stripes_b, dummy(1:0))
  END SUBROUTINE test_intersection10

  SUBROUTINE test_intersection11
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 5, 20) /), &
         stripes_b(1) = (/ xt_stripe(1, 7, 15) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(3) = (/ 15_xi, 50_xi, 85_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection11

  ! both ranges overlap in range but have no
  ! indices in common because of stride
  SUBROUTINE test_intersection12
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(34, 29, 12) /), &
         stripes_b(1) = (/ xt_stripe(36, 7, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: dummy(1) = (/ -1_xi /)

    CALL test_intersection(stripes_a, stripes_b, dummy(1:0))
  END SUBROUTINE test_intersection12

  ! same as test_intersection12 but with negative stride
  SUBROUTINE test_intersection13
    TYPE(xt_stripe), PARAMETER :: &
         stripes_a(1) = (/ xt_stripe(353, -29, 12) /), &
         stripes_b(1) = (/ xt_stripe(36, 7, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: dummy(1) = (/ -1_xi /)

    CALL test_intersection(stripes_a, stripes_b, dummy(1:0))
  END SUBROUTINE test_intersection13

  SUBROUTINE test_intersection14
    TYPE(xt_stripe), PARAMETER :: &
         stripes_a(1) = (/ xt_stripe(95, -29, 2) /), &
         stripes_b(1) = (/ xt_stripe(81, 14, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(1) = (/ 95_xi /)

    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection14

  SUBROUTINE test_intersection15
    TYPE(xt_stripe), PARAMETER :: &
         stripes_a(1) = (/ xt_stripe(546, 14, 2) /), &
         stripes_b(1) = (/ xt_stripe(354, 206, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(1) = (/ 560_xi /)

    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection15

  SUBROUTINE test_intersection_stripe2vec
    INTEGER, PARAMETER :: num_stripes = 3
    TYPE(xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(4, 1, 1), xt_stripe(5, 1, 1), xt_stripe(10, -10, 2) /)
    TYPE(xt_idxlist) :: idxvec_a, idxvec_b, intersection
    INTEGER(xt_int_kind), PARAMETER :: index_vector(1) = (/ 5_xi /)
    INTEGER(xt_int_kind) :: intersection_idx
    LOGICAL :: not_found
    idxvec_a = xt_idxvec_from_stripes_new(stripes)
    idxvec_b = xt_idxvec_new(index_vector)
    intersection = xt_idxlist_get_intersection(idxvec_a, idxvec_b)
    IF (xt_idxlist_get_num_indices(intersection) /= 1) &
         CALL test_abort("unexpected number of indices in intersection!", &
         filename, __LINE__)
    not_found = xt_idxlist_get_index_at_position(intersection, 0, &
         intersection_idx)
    IF (not_found .OR. intersection_idx /= index_vector(1)) &
         CALL test_abort("unexpected index in intersection!", &
         filename, __LINE__)
    CALL xt_idxlist_delete(intersection)
    CALL xt_idxlist_delete(idxvec_a)
    CALL xt_idxlist_delete(idxvec_b)
  END SUBROUTINE test_intersection_stripe2vec

  SUBROUTINE test_idxlist_stripes_pos_ext1
    INTEGER, PARAMETER :: num_indices = 223
    INTEGER(xt_int_kind), PARAMETER :: index_vector(num_indices) = (/ &
         3375_xi, 3376_xi, 3379_xi, 3380_xi, 3381_xi, 3387_xi, 3388_xi, &
         3389_xi, 3390_xi, 3391_xi, 3392_xi, 3393_xi, 3421_xi, 3422_xi, &
         3423_xi, 3424_xi, 3425_xi, 3426_xi, 3427_xi, 3444_xi, 3458_xi, &
         3459_xi, 3461_xi, 3462_xi, 3463_xi, 3464_xi, 3465_xi, 3466_xi, &
         3467_xi, 3468_xi, 3469_xi, 3470_xi, 3471_xi, 3472_xi, 3473_xi, &
         3474_xi, 3475_xi, 3476_xi, 3477_xi, 3478_xi, 3479_xi, 3480_xi, &
         3529_xi, 3606_xi, 3607_xi, 3608_xi, 3611_xi, 3612_xi, 3613_xi, &
         3614_xi, 3617_xi, 3620_xi, 3621_xi, 3622_xi, 3623_xi, 3624_xi, &
         3625_xi, 3626_xi, 3627_xi, 3628_xi, 3629_xi, 3630_xi, 3631_xi, &
         3684_xi, 3685_xi, 3686_xi, 3687_xi, 3688_xi, 3689_xi, 3690_xi, &
         3691_xi, 3692_xi, 3693_xi, 3694_xi, 3695_xi, 3696_xi, 3697_xi, &
         3698_xi, 3699_xi, 3700_xi, 3701_xi, 3702_xi, 3703_xi, 3704_xi, &
         3705_xi, 3706_xi, 3707_xi, 3708_xi, 3709_xi, 3713_xi, 3714_xi, &
         3715_xi, 3716_xi, 3717_xi, 3718_xi, 3719_xi, 3720_xi, 3721_xi, &
         3722_xi, 3723_xi, 3724_xi, 3725_xi, 3726_xi, 3727_xi, 3728_xi, &
         3729_xi, 3730_xi, 3731_xi, 3741_xi, 3742_xi, 3931_xi, 3932_xi, &
         3374_xi, 3382_xi, 3385_xi, 3394_xi, 3404_xi, 3408_xi, 3412_xi, &
         3440_xi, 3443_xi, 3457_xi, 3481_xi, 3483_xi, 3527_xi, 3619_xi, &
         3735_xi, 3743_xi, 3925_xi, 3930_xi, 3377_xi, 3378_xi, 3383_xi, &
         3384_xi, 3386_xi, 3395_xi, 3397_xi, 3398_xi, 3400_xi, 3402_xi, &
         3403_xi, 3407_xi, 3409_xi, 3410_xi, 3413_xi, 3420_xi, 3441_xi, &
         3442_xi, 3445_xi, 3448_xi, 3449_xi, 3451_xi, 3460_xi, 3482_xi, &
         3519_xi, 3520_xi, 3526_xi, 3528_xi, 3530_xi, 3592_xi, 3593_xi, &
         3595_xi, 3596_xi, 3597_xi, 3609_xi, 3610_xi, 3615_xi, 3616_xi, &
         3618_xi, 3644_xi, 3710_xi, 3711_xi, 3712_xi, 3732_xi, 3733_xi, &
         3736_xi, 3737_xi, 3748_xi, 3749_xi, 3753_xi, 3754_xi, 3759_xi, &
         3760_xi, 3766_xi, 3767_xi, 3919_xi, 3920_xi, 3924_xi, 3926_xi, &
         3933_xi, 3934_xi, 2589_xi, 2602_xi, 2680_xi, 3326_xi, 3340_xi, &
         3341_xi, 3396_xi, 3401_xi, 3411_xi, 3414_xi, 3418_xi, 3446_xi, &
         3447_xi, 3450_xi, 3515_xi, 3521_xi, 3525_xi, 3582_xi, 3590_xi, &
         3591_xi, 3594_xi, 3642_xi, 3734_xi, 3738_xi, 3747_xi, 3750_xi, &
         3761_xi, 3765_xi, 3865_xi, 3918_xi, 3923_xi, 3935_xi /)
    INTEGER, PARAMETER :: num_stripes = 26
    TYPE(xt_stripe), PARAMETER :: stripes(num_stripes) = (/ &
         xt_stripe(3326, 14,  2), xt_stripe(3341, 33,  1), &
         xt_stripe(3374,  1, 25), xt_stripe(3400,  1,  5), &
         xt_stripe(3407,  1,  8), xt_stripe(3418,  2,  1), &
         xt_stripe(3420,  1,  8), xt_stripe(3440,  1, 12), &
         xt_stripe(3457,  1, 27), xt_stripe(3515,  4,  1), &
         xt_stripe(3519,  1,  3), xt_stripe(3525,  1,  6), &
         xt_stripe(3582,  8,  1), xt_stripe(3590,  1,  8), &
         xt_stripe(3606,  1, 26), xt_stripe(3642,  2,  2), &
         xt_stripe(3684,  1, 55), xt_stripe(3741,  1,  3), &
         xt_stripe(3747,  1,  4), xt_stripe(3753,  1,  2), &
         xt_stripe(3759,  1,  3), xt_stripe(3765,  1,  3), &
         xt_stripe(3865, 53,  1), xt_stripe(3918,  1,  3), &
         xt_stripe(3923,  1,  4), xt_stripe(3930,  1,  6) /)
    TYPE(xt_idxlist) :: idxlist

    idxlist = xt_idxvec_new(index_vector, num_indices)
    CALL check_idxlist_stripes_pos_ext(idxlist, stripes)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_idxlist_stripes_pos_ext1

  SUBROUTINE test_idxlist_stripes_pos_ext2
    INTEGER, PARAMETER :: num_indices = 201
    INTEGER(xt_int_kind), PARAMETER :: index_vector(num_indices) = (/ &
         &  178_xi,  179_xi,  180_xi,  181_xi,  182_xi,  183_xi,  184_xi, &
         &  186_xi,  187_xi,  188_xi,  189_xi,  190_xi,  194_xi,  195_xi, &
         &  196_xi,  197_xi,  198_xi,  199_xi,  200_xi,  201_xi,  202_xi, &
         &  203_xi,  204_xi,  205_xi,  206_xi,  207_xi,  208_xi,  209_xi, &
         &  210_xi,  211_xi,  212_xi,  217_xi,  223_xi,  426_xi,  428_xi, &
         &  429_xi,  430_xi,  434_xi,  435_xi,  436_xi,  437_xi,  438_xi, &
         &  439_xi,  440_xi,  442_xi,  443_xi,  444_xi,  445_xi,  446_xi, &
         &  447_xi,  448_xi,  449_xi,  450_xi,  451_xi,  452_xi,  453_xi, &
         &  454_xi,  455_xi,  456_xi,  457_xi,  458_xi,  670_xi,  671_xi, &
         &  672_xi,  673_xi,  674_xi,  675_xi,  676_xi,  677_xi,  682_xi, &
         &  684_xi,  685_xi,  686_xi,  687_xi,  688_xi,  689_xi,  690_xi, &
         &  692_xi,  695_xi,  703_xi,  704_xi,  705_xi,  706_xi,  707_xi, &
         &  894_xi,  895_xi,  896_xi,  897_xi,  898_xi,  899_xi,  900_xi, &
         &  901_xi,  906_xi,  907_xi,  908_xi,  913_xi,  915_xi,  921_xi, &
         &  922_xi,  923_xi,  924_xi,  925_xi,  926_xi,  927_xi, 1096_xi, &
         & 1097_xi, 1098_xi, 1099_xi, 1100_xi, 1101_xi, 1102_xi, 1103_xi, &
         & 1107_xi, 1108_xi, 1109_xi, 1110_xi, 1111_xi, 1113_xi, 1114_xi, &
         & 1119_xi, 1120_xi, 1121_xi, 2095_xi, 2096_xi, 2097_xi, 2098_xi, &
         & 2100_xi, 2102_xi, 2103_xi, 2104_xi, 2105_xi, 2107_xi, 2108_xi, &
         & 2109_xi, 2110_xi, 2112_xi, 2118_xi, 2120_xi, 2121_xi, 2122_xi, &
         & 2123_xi, 2124_xi, 2125_xi, 2127_xi, 2128_xi, 2129_xi, 2130_xi, &
         & 2134_xi, 2140_xi, 2141_xi, 2142_xi, 2143_xi, 2145_xi, 2148_xi, &
         & 2149_xi, 2151_xi, 2152_xi, 2153_xi, 2154_xi, 2155_xi, 2156_xi, &
         &  683_xi,  691_xi,  903_xi,  914_xi, 1105_xi, 1115_xi, 2099_xi, &
         & 2106_xi, 2111_xi, 2115_xi, 2126_xi, 2132_xi, 2139_xi, 2144_xi, &
         & 2147_xi, 2150_xi, 2305_xi,  427_xi,  465_xi,  466_xi,  678_xi, &
         &  693_xi,  902_xi,  909_xi, 1104_xi, 1112_xi, 2101_xi, 2113_xi, &
         & 2114_xi, 2116_xi, 2117_xi, 2119_xi, 2131_xi, 2136_xi, 2138_xi, &
         & 2146_xi, 2297_xi, 2302_xi, 2304_xi, 2307_xi /)
    INTEGER, PARAMETER :: num_stripes = 8
    TYPE(xt_stripe), PARAMETER :: stripes(num_stripes) = (/ &
      xt_stripe(670, 1,  9), xt_stripe(682, 1, 12), &
      xt_stripe(695, 8,  1), xt_stripe(703, 1,  5), &
      xt_stripe(894, 1, 10), xt_stripe(906, 1,  4), &
      xt_stripe(913, 1,  3), xt_stripe(921, 1,  7) /)
    TYPE(xt_idxlist) :: idxlist

    idxlist = xt_idxvec_new(index_vector)
    CALL check_idxlist_stripes_pos_ext(idxlist, stripes)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_idxlist_stripes_pos_ext2

  SUBROUTINE test_idxlist_stripes_pos_ext3
    INTEGER, PARAMETER :: num_indices = 1144
    INTEGER(xt_int_kind), PARAMETER :: index_vector(num_indices) = (/ &
      2055_xi, 2056_xi, 2060_xi, 2193_xi, 2199_xi, 2203_xi, 2211_xi, 2212_xi, &
      2278_xi, 2281_xi, 2311_xi, 2312_xi, 2316_xi, 2317_xi, 2322_xi, 2332_xi, &
      2447_xi, 2448_xi, 2452_xi, 2585_xi, 2591_xi, 2595_xi, 2603_xi, 2604_xi, &
      2670_xi, 2673_xi, 2703_xi, 2704_xi, 2708_xi, 2709_xi, 2714_xi, 2724_xi, &
      2839_xi, 2840_xi, 2844_xi, 2977_xi, 2983_xi, 2987_xi, 2995_xi, 2996_xi, &
      3062_xi, 3065_xi, 3095_xi, 3096_xi, 3100_xi, 3101_xi, 3106_xi, 3116_xi, &
      3231_xi, 3232_xi, 3236_xi, 3369_xi, 3375_xi, 3379_xi, 3387_xi, 3388_xi, &
      3454_xi, 3457_xi, 3487_xi, 3488_xi, 3492_xi, 3493_xi, 3498_xi, 3508_xi, &
      3623_xi, 3624_xi, 3628_xi, 3761_xi, 3767_xi, 3771_xi, 3779_xi, 3780_xi, &
      3846_xi, 3849_xi, 3879_xi, 3880_xi, 3884_xi, 3885_xi, 3890_xi, 3900_xi, &
      3997_xi, 4001_xi, 4002_xi, 4053_xi, 4057_xi, 4084_xi, 4085_xi, 4092_xi, &
      4102_xi, 4188_xi, 4192_xi, 4201_xi, 4373_xi, 4377_xi, 4378_xi, 4429_xi, &
      4433_xi, 4460_xi, 4461_xi, 4468_xi, 4478_xi, 4564_xi, 4568_xi, 4577_xi, &
      4749_xi, 4753_xi, 4754_xi, 4805_xi, 4809_xi, 4836_xi, 4837_xi, 4844_xi, &
      4854_xi, 4945_xi, 4953_xi, 5125_xi, 5129_xi, 5130_xi, 5181_xi, 5185_xi, &
      5212_xi, 5213_xi, 5220_xi, 5230_xi, 5321_xi, 5329_xi, 5501_xi, 5505_xi, &
      5506_xi, 5557_xi, 5561_xi, 5588_xi, 5589_xi, 5596_xi, 5606_xi, 5697_xi, &
      5705_xi,  162_xi,  163_xi,  166_xi,  168_xi,  171_xi,  172_xi,  173_xi, &
       177_xi,  181_xi,  362_xi,  363_xi,  367_xi,  369_xi,  375_xi,  378_xi, &
       382_xi,  383_xi,  386_xi,  570_xi,  571_xi,  574_xi,  576_xi,  579_xi, &
       580_xi,  581_xi,  585_xi,  589_xi,  758_xi,  759_xi,  763_xi,  765_xi, &
       769_xi,  774_xi,  775_xi,  778_xi,  962_xi,  963_xi,  966_xi,  968_xi, &
       971_xi,  972_xi,  973_xi,  977_xi,  981_xi, 1150_xi, 1151_xi, 1155_xi, &
      1157_xi, 1161_xi, 1166_xi, 1167_xi, 1170_xi, 1354_xi, 1355_xi, 1358_xi, &
      1360_xi, 1363_xi, 1364_xi, 1365_xi, 1369_xi, 1373_xi, 1542_xi, 1543_xi, &
      1547_xi, 1549_xi, 1553_xi, 1558_xi, 1559_xi, 1562_xi, 1746_xi, 1747_xi, &
      1750_xi, 1752_xi, 1755_xi, 1756_xi, 1757_xi, 1761_xi, 1918_xi, 1919_xi, &
      1923_xi, 1925_xi, 1929_xi, 1934_xi, 1935_xi, 1938_xi, 1988_xi, 1989_xi, &
      2024_xi, 2025_xi, 2032_xi, 2033_xi, 2036_xi, 2038_xi, 2039_xi, 2048_xi, &
      2049_xi, 2053_xi, 2054_xi, 2057_xi, 2058_xi, 2061_xi, 2076_xi, 2077_xi, &
      2091_xi, 2092_xi, 2093_xi, 2095_xi, 2097_xi, 2126_xi, 2127_xi, 2144_xi, &
      2145_xi, 2149_xi, 2150_xi, 2156_xi, 2198_xi, 2204_xi, 2205_xi, 2207_xi, &
      2245_xi, 2253_xi, 2254_xi, 2256_xi, 2268_xi, 2269_xi, 2277_xi, 2279_xi, &
      2280_xi, 2283_xi, 2287_xi, 2298_xi, 2299_xi, 2307_xi, 2308_xi, 2309_xi, &
      2310_xi, 2333_xi, 2334_xi, 2380_xi, 2381_xi, 2416_xi, 2417_xi, 2424_xi, &
      2425_xi, 2428_xi, 2430_xi, 2431_xi, 2440_xi, 2441_xi, 2445_xi, 2446_xi, &
      2449_xi, 2450_xi, 2453_xi, 2468_xi, 2469_xi, 2483_xi, 2484_xi, 2485_xi, &
      2487_xi, 2489_xi, 2518_xi, 2519_xi, 2536_xi, 2537_xi, 2541_xi, 2542_xi, &
      2548_xi, 2590_xi, 2596_xi, 2597_xi, 2599_xi, 2637_xi, 2645_xi, 2646_xi, &
      2648_xi, 2660_xi, 2661_xi, 2669_xi, 2671_xi, 2672_xi, 2675_xi, 2679_xi, &
      2690_xi, 2691_xi, 2699_xi, 2700_xi, 2701_xi, 2702_xi, 2725_xi, 2726_xi, &
      2772_xi, 2773_xi, 2808_xi, 2809_xi, 2816_xi, 2817_xi, 2820_xi, 2822_xi, &
      2823_xi, 2832_xi, 2833_xi, 2837_xi, 2838_xi, 2841_xi, 2842_xi, 2845_xi, &
      2860_xi, 2861_xi, 2875_xi, 2876_xi, 2877_xi, 2879_xi, 2881_xi, 2910_xi, &
      2911_xi, 2928_xi, 2929_xi, 2933_xi, 2934_xi, 2940_xi, 2982_xi, 2988_xi, &
      2989_xi, 2991_xi, 3029_xi, 3037_xi, 3038_xi, 3040_xi, 3052_xi, 3053_xi, &
      3061_xi, 3063_xi, 3064_xi, 3067_xi, 3071_xi, 3082_xi, 3083_xi, 3091_xi, &
      3092_xi, 3093_xi, 3094_xi, 3117_xi, 3118_xi, 3164_xi, 3165_xi, 3200_xi, &
      3201_xi, 3208_xi, 3209_xi, 3212_xi, 3214_xi, 3215_xi, 3224_xi, 3225_xi, &
      3229_xi, 3230_xi, 3233_xi, 3234_xi, 3237_xi, 3252_xi, 3253_xi, 3267_xi, &
      3268_xi, 3269_xi, 3271_xi, 3273_xi, 3302_xi, 3303_xi, 3320_xi, 3321_xi, &
      3325_xi, 3326_xi, 3332_xi, 3374_xi, 3380_xi, 3381_xi, 3383_xi, 3421_xi, &
      3429_xi, 3430_xi, 3432_xi, 3444_xi, 3445_xi, 3453_xi, 3455_xi, 3456_xi, &
      3459_xi, 3463_xi, 3474_xi, 3475_xi, 3483_xi, 3484_xi, 3485_xi, 3486_xi, &
      3509_xi, 3510_xi, 3556_xi, 3557_xi, 3592_xi, 3593_xi, 3600_xi, 3601_xi, &
      3604_xi, 3606_xi, 3607_xi, 3616_xi, 3617_xi, 3621_xi, 3622_xi, 3625_xi, &
      3626_xi, 3629_xi, 3644_xi, 3645_xi, 3659_xi, 3660_xi, 3661_xi, 3663_xi, &
      3665_xi, 3694_xi, 3695_xi, 3712_xi, 3713_xi, 3717_xi, 3718_xi, 3724_xi, &
      3766_xi, 3772_xi, 3773_xi, 3775_xi, 3813_xi, 3821_xi, 3822_xi, 3824_xi, &
      3836_xi, 3837_xi, 3845_xi, 3847_xi, 3848_xi, 3851_xi, 3855_xi, 3866_xi, &
      3867_xi, 3875_xi, 3876_xi, 3877_xi, 3878_xi, 3901_xi, 3902_xi, 3948_xi, &
      3949_xi, 3984_xi, 3985_xi, 3992_xi, 3993_xi, 3996_xi, 3998_xi, 3999_xi, &
      4008_xi, 4009_xi, 4013_xi, 4014_xi, 4017_xi, 4018_xi, 4021_xi, 4036_xi, &
      4037_xi, 4051_xi, 4052_xi, 4054_xi, 4055_xi, 4058_xi, 4090_xi, 4091_xi, &
      4093_xi, 4108_xi, 4109_xi, 4112_xi, 4113_xi, 4114_xi, 4158_xi, 4164_xi, &
      4165_xi, 4193_xi, 4199_xi, 4200_xi, 4212_xi, 4213_xi, 4222_xi, 4223_xi, &
      4225_xi, 4227_xi, 4231_xi, 4242_xi, 4243_xi, 4250_xi, 4251_xi, 4271_xi, &
      4272_xi, 4274_xi, 4324_xi, 4325_xi, 4360_xi, 4361_xi, 4368_xi, 4369_xi, &
      4372_xi, 4374_xi, 4375_xi, 4384_xi, 4385_xi, 4389_xi, 4390_xi, 4393_xi, &
      4394_xi, 4397_xi, 4412_xi, 4413_xi, 4427_xi, 4428_xi, 4430_xi, 4431_xi, &
      4434_xi, 4466_xi, 4467_xi, 4469_xi, 4484_xi, 4485_xi, 4488_xi, 4489_xi, &
      4490_xi, 4534_xi, 4540_xi, 4541_xi, 4569_xi, 4575_xi, 4576_xi, 4588_xi, &
      4589_xi, 4598_xi, 4599_xi, 4601_xi, 4603_xi, 4607_xi, 4618_xi, 4619_xi, &
      4626_xi, 4627_xi, 4647_xi, 4648_xi, 4650_xi, 4700_xi, 4701_xi, 4736_xi, &
      4737_xi, 4744_xi, 4745_xi, 4748_xi, 4750_xi, 4751_xi, 4760_xi, 4761_xi, &
      4765_xi, 4766_xi, 4769_xi, 4770_xi, 4773_xi, 4788_xi, 4789_xi, 4803_xi, &
      4804_xi, 4806_xi, 4807_xi, 4810_xi, 4842_xi, 4843_xi, 4845_xi, 4860_xi, &
      4861_xi, 4864_xi, 4865_xi, 4866_xi, 4910_xi, 4916_xi, 4917_xi, 4951_xi, &
      4952_xi, 4964_xi, 4965_xi, 4974_xi, 4975_xi, 4977_xi, 4979_xi, 4983_xi, &
      4994_xi, 4995_xi, 5002_xi, 5003_xi, 5023_xi, 5024_xi, 5026_xi, 5076_xi, &
      5077_xi, 5112_xi, 5113_xi, 5120_xi, 5121_xi, 5124_xi, 5126_xi, 5127_xi, &
      5136_xi, 5137_xi, 5141_xi, 5142_xi, 5145_xi, 5146_xi, 5149_xi, 5164_xi, &
      5165_xi, 5179_xi, 5180_xi, 5182_xi, 5183_xi, 5186_xi, 5218_xi, 5219_xi, &
      5221_xi, 5236_xi, 5237_xi, 5240_xi, 5241_xi, 5242_xi, 5286_xi, 5292_xi, &
      5293_xi, 5327_xi, 5328_xi, 5340_xi, 5341_xi, 5350_xi, 5351_xi, 5353_xi, &
      5355_xi, 5359_xi, 5370_xi, 5371_xi, 5378_xi, 5379_xi, 5399_xi, 5400_xi, &
      5402_xi, 5452_xi, 5453_xi, 5488_xi, 5489_xi, 5496_xi, 5497_xi, 5500_xi, &
      5502_xi, 5503_xi, 5512_xi, 5513_xi, 5517_xi, 5518_xi, 5521_xi, 5522_xi, &
      5525_xi, 5540_xi, 5541_xi, 5555_xi, 5556_xi, 5558_xi, 5559_xi, 5562_xi, &
      5594_xi, 5595_xi, 5597_xi, 5612_xi, 5613_xi, 5616_xi, 5617_xi, 5618_xi, &
      5662_xi, 5668_xi, 5669_xi, 5703_xi, 5704_xi, 5716_xi, 5717_xi, 5726_xi, &
      5727_xi, 5729_xi, 5731_xi, 5735_xi, 5746_xi, 5747_xi, 5754_xi, 5755_xi, &
      5775_xi, 5776_xi, 5778_xi, 5958_xi, 5959_xi, 5962_xi, 5964_xi, 5967_xi, &
      5968_xi, 5971_xi, 5973_xi, 6154_xi, 6155_xi, 6159_xi, 6161_xi, 6167_xi, &
      6170_xi, 6172_xi, 6173_xi, 6350_xi, 6351_xi, 6354_xi, 6356_xi, 6359_xi, &
      6360_xi, 6363_xi, 6530_xi, 6531_xi, 6535_xi, 6537_xi, 6543_xi, 6546_xi, &
      6548_xi, 6549_xi, 6726_xi, 6727_xi, 6730_xi, 6732_xi, 6735_xi, 6736_xi, &
      6739_xi, 6906_xi, 6907_xi, 6911_xi, 6913_xi, 6919_xi, 6922_xi, 6924_xi, &
      6925_xi, 7102_xi, 7103_xi, 7106_xi, 7108_xi, 7111_xi, 7112_xi, 7115_xi, &
      7282_xi, 7283_xi, 7287_xi, 7289_xi, 7295_xi, 7298_xi, 7300_xi, 7301_xi, &
      7478_xi, 7479_xi, 7482_xi, 7484_xi, 7487_xi, 7488_xi, 7491_xi, 7646_xi, &
      7647_xi, 7651_xi, 7653_xi, 7657_xi, 7660_xi, 7661_xi,  130_xi,  161_xi, &
       169_xi,  170_xi,  336_xi,  361_xi,  366_xi,  384_xi,  538_xi,  569_xi, &
       577_xi,  578_xi,  736_xi,  757_xi,  762_xi,  776_xi,  930_xi,  961_xi, &
       969_xi,  970_xi, 1128_xi, 1149_xi, 1154_xi, 1168_xi, 1322_xi, 1353_xi, &
      1361_xi, 1362_xi, 1520_xi, 1541_xi, 1546_xi, 1560_xi, 1714_xi, 1745_xi, &
      1753_xi, 1754_xi, 1896_xi, 1917_xi, 1922_xi, 1936_xi, 1985_xi, 2019_xi, &
      2031_xi, 2035_xi, 2040_xi, 2044_xi, 2052_xi, 2059_xi, 2062_xi, 2071_xi, &
      2087_xi, 2090_xi, 2094_xi, 2140_xi, 2148_xi, 2153_xi, 2157_xi, 2206_xi, &
      2257_xi, 2263_xi, 2267_xi, 2284_xi, 2288_xi, 2293_xi, 2295_xi, 2305_xi, &
      2306_xi, 2377_xi, 2411_xi, 2423_xi, 2427_xi, 2432_xi, 2436_xi, 2444_xi, &
      2451_xi, 2454_xi, 2463_xi, 2479_xi, 2482_xi, 2486_xi, 2532_xi, 2540_xi, &
      2545_xi, 2549_xi, 2598_xi, 2649_xi, 2655_xi, 2659_xi, 2676_xi, 2680_xi, &
      2685_xi, 2687_xi, 2697_xi, 2698_xi, 2769_xi, 2803_xi, 2815_xi, 2819_xi, &
      2824_xi, 2828_xi, 2836_xi, 2843_xi, 2846_xi, 2855_xi, 2871_xi, 2874_xi, &
      2878_xi, 2924_xi, 2932_xi, 2937_xi, 2941_xi, 2990_xi, 3041_xi, 3047_xi, &
      3051_xi, 3068_xi, 3072_xi, 3077_xi, 3079_xi, 3089_xi, 3090_xi, 3161_xi, &
      3195_xi, 3207_xi, 3211_xi, 3216_xi, 3220_xi, 3228_xi, 3235_xi, 3238_xi, &
      3247_xi, 3263_xi, 3266_xi, 3270_xi, 3316_xi, 3324_xi, 3329_xi, 3333_xi, &
      3382_xi, 3433_xi, 3439_xi, 3443_xi, 3460_xi, 3464_xi, 3469_xi, 3471_xi, &
      3481_xi, 3482_xi, 3553_xi, 3587_xi, 3599_xi, 3603_xi, 3608_xi, 3612_xi, &
      3620_xi, 3627_xi, 3630_xi, 3639_xi, 3655_xi, 3658_xi, 3662_xi, 3708_xi, &
      3716_xi, 3721_xi, 3725_xi, 3774_xi, 3825_xi, 3831_xi, 3835_xi, 3852_xi, &
      3856_xi, 3861_xi, 3863_xi, 3873_xi, 3874_xi, 3945_xi, 3979_xi, 3991_xi, &
      3995_xi, 4000_xi, 4004_xi, 4012_xi, 4019_xi, 4022_xi, 4031_xi, 4033_xi, &
      4047_xi, 4050_xi, 4104_xi, 4106_xi, 4115_xi, 4207_xi, 4221_xi, 4228_xi, &
      4232_xi, 4237_xi, 4249_xi, 4252_xi, 4321_xi, 4355_xi, 4367_xi, 4371_xi, &
      4376_xi, 4380_xi, 4388_xi, 4395_xi, 4398_xi, 4407_xi, 4409_xi, 4423_xi, &
      4426_xi, 4480_xi, 4482_xi, 4491_xi, 4583_xi, 4597_xi, 4604_xi, 4608_xi, &
      4613_xi, 4625_xi, 4628_xi, 4697_xi, 4731_xi, 4743_xi, 4747_xi, 4752_xi, &
      4756_xi, 4764_xi, 4771_xi, 4774_xi, 4783_xi, 4785_xi, 4799_xi, 4802_xi, &
      4856_xi, 4858_xi, 4867_xi, 4959_xi, 4973_xi, 4980_xi, 4984_xi, 4989_xi, &
      5001_xi, 5004_xi, 5073_xi, 5107_xi, 5119_xi, 5123_xi, 5128_xi, 5132_xi, &
      5140_xi, 5147_xi, 5150_xi, 5159_xi, 5161_xi, 5175_xi, 5178_xi, 5232_xi, &
      5234_xi, 5243_xi, 5335_xi, 5349_xi, 5356_xi, 5360_xi, 5365_xi, 5377_xi, &
      5380_xi, 5449_xi, 5483_xi, 5495_xi, 5499_xi, 5504_xi, 5508_xi, 5516_xi, &
      5523_xi, 5526_xi, 5535_xi, 5537_xi, 5551_xi, 5554_xi, 5608_xi, 5610_xi, &
      5619_xi, 5711_xi, 5725_xi, 5732_xi, 5736_xi, 5741_xi, 5753_xi, 5756_xi, &
      5930_xi, 5957_xi, 5965_xi, 5966_xi, 6128_xi, 6153_xi, 6158_xi, 6174_xi, &
      6322_xi, 6349_xi, 6357_xi, 6358_xi, 6504_xi, 6529_xi, 6534_xi, 6550_xi, &
      6698_xi, 6725_xi, 6733_xi, 6734_xi, 6880_xi, 6905_xi, 6910_xi, 6926_xi, &
      7074_xi, 7101_xi, 7109_xi, 7110_xi, 7256_xi, 7281_xi, 7286_xi, 7302_xi, &
      7450_xi, 7477_xi, 7485_xi, 7486_xi, 7624_xi, 7645_xi, 7650_xi, 7662_xi /)
    INTEGER, PARAMETER :: num_stripes = 187
    TYPE(xt_stripe), PARAMETER :: stripes(num_stripes) = (/ &
      xt_stripe(173, 408, 2), xt_stripe(973, 392, 3), xt_stripe(1985, 4, 2), &
      xt_stripe(2044, 4, 2), xt_stripe(2049, 3, 1), xt_stripe(2052, 1, 9), &
      xt_stripe(2062, 131, 2), xt_stripe(2198, 1, 2), xt_stripe(2203, 1, 5), &
      xt_stripe(2211, 1, 2), xt_stripe(2263, 4, 1), xt_stripe(2267, 1, 3), &
      xt_stripe(2277, 1, 5), xt_stripe(2283, 1, 2), xt_stripe(2287, 1, 2), &
      xt_stripe(2293, 2, 2), xt_stripe(2298, 1, 2), xt_stripe(2305, 1, 8), &
      xt_stripe(2316, 1, 2), xt_stripe(2322, 10, 1), xt_stripe(2332, 1, 3), &
      xt_stripe(2377, 4, 2), xt_stripe(2436, 4, 2), xt_stripe(2441, 3, 1), &
      xt_stripe(2444, 1, 9), xt_stripe(2454, 131, 2), xt_stripe(2590, 1, 2), &
      xt_stripe(2595, 1, 5), xt_stripe(2603, 1, 2), xt_stripe(2655, 4, 1), &
      xt_stripe(2659, 1, 3), xt_stripe(2669, 1, 5), xt_stripe(2675, 1, 2), &
      xt_stripe(2679, 1, 2), xt_stripe(2685, 2, 2), xt_stripe(2690, 1, 2), &
      xt_stripe(2697, 1, 8), xt_stripe(2708, 1, 2), xt_stripe(2714, 10, 1), &
      xt_stripe(2724, 1, 3), xt_stripe(2769, 4, 2), xt_stripe(2828, 4, 2), &
      xt_stripe(2833, 3, 1), xt_stripe(2836, 1, 9), xt_stripe(2846, 131, 2), &
      xt_stripe(2982, 1, 2), xt_stripe(2987, 1, 5), xt_stripe(2995, 1, 2), &
      xt_stripe(3047, 4, 1), xt_stripe(3051, 1, 3), xt_stripe(3061, 1, 5), &
      xt_stripe(3067, 1, 2), xt_stripe(3071, 1, 2), xt_stripe(3077, 2, 2), &
      xt_stripe(3082, 1, 2), xt_stripe(3089, 1, 8), xt_stripe(3100, 1, 2), &
      xt_stripe(3106, 10, 1), xt_stripe(3116, 1, 3), xt_stripe(3161, 4, 2), &
      xt_stripe(3220, 4, 2), xt_stripe(3225, 3, 1), xt_stripe(3228, 1, 9), &
      xt_stripe(3238, 131, 2), xt_stripe(3374, 1, 2), xt_stripe(3379, 1, 5), &
      xt_stripe(3387, 1, 2), xt_stripe(3439, 4, 1), xt_stripe(3443, 1, 3), &
      xt_stripe(3453, 1, 5), xt_stripe(3459, 1, 2), xt_stripe(3463, 1, 2), &
      xt_stripe(3469, 2, 2), xt_stripe(3474, 1, 2), xt_stripe(3481, 1, 8), &
      xt_stripe(3492, 1, 2), xt_stripe(3498, 10, 1), xt_stripe(3508, 1, 3), &
      xt_stripe(3553, 4, 2), xt_stripe(3612, 4, 2), xt_stripe(3617, 3, 1), &
      xt_stripe(3620, 1, 9), xt_stripe(3630, 131, 2), xt_stripe(3766, 1, 2), &
      xt_stripe(3771, 1, 5), xt_stripe(3779, 1, 2), xt_stripe(3831, 4, 1), &
      xt_stripe(3835, 1, 3), xt_stripe(3845, 1, 5), xt_stripe(3851, 1, 2), &
      xt_stripe(3855, 1, 2), xt_stripe(3861, 2, 2), xt_stripe(3866, 1, 2), &
      xt_stripe(3873, 1, 8), xt_stripe(3884, 1, 2), xt_stripe(3890, 10, 1), &
      xt_stripe(3900, 1, 3), xt_stripe(3945, 3, 2), xt_stripe(3979, 5, 2), &
      xt_stripe(3985, 6, 1), xt_stripe(3991, 1, 3), xt_stripe(3995, 2, 1), &
      xt_stripe(3997, 1, 6), xt_stripe(4031, 2, 2), xt_stripe(4036, 1, 2), &
      xt_stripe(4047, 3, 1), xt_stripe(4050, 1, 6), xt_stripe(4057, 1, 2), &
      xt_stripe(4084, 1, 2), xt_stripe(4090, 1, 4), xt_stripe(4102, 2, 4), &
      xt_stripe(4109, 3, 1), xt_stripe(4112, 1, 4), xt_stripe(4188, 4, 2), &
      xt_stripe(4193, 6, 1), xt_stripe(4199, 1, 3), xt_stripe(4321, 3, 2), &
      xt_stripe(4355, 5, 2), xt_stripe(4361, 6, 1), xt_stripe(4367, 1, 3), &
      xt_stripe(4371, 2, 1), xt_stripe(4373, 1, 6), xt_stripe(4407, 2, 2), &
      xt_stripe(4412, 1, 2), xt_stripe(4423, 3, 1), xt_stripe(4426, 1, 6), &
      xt_stripe(4433, 1, 2), xt_stripe(4460, 1, 2), xt_stripe(4466, 1, 4), &
      xt_stripe(4478, 2, 4), xt_stripe(4485, 3, 1), xt_stripe(4488, 1, 4), &
      xt_stripe(4564, 4, 2), xt_stripe(4569, 6, 1), xt_stripe(4575, 1, 3), &
      xt_stripe(4697, 3, 2), xt_stripe(4731, 5, 2), xt_stripe(4737, 6, 1), &
      xt_stripe(4743, 1, 3), xt_stripe(4747, 2, 1), xt_stripe(4749, 1, 6), &
      xt_stripe(4783, 2, 2), xt_stripe(4788, 1, 2), xt_stripe(4799, 3, 1), &
      xt_stripe(4802, 1, 6), xt_stripe(4809, 1, 2), xt_stripe(4836, 1, 2), &
      xt_stripe(4842, 1, 4), xt_stripe(4854, 2, 4), xt_stripe(4861, 3, 1), &
      xt_stripe(4864, 1, 4), xt_stripe(4945, 6, 1), xt_stripe(4951, 1, 3), &
      xt_stripe(5107, 5, 2), xt_stripe(5113, 6, 1), xt_stripe(5119, 1, 3), &
      xt_stripe(5123, 2, 1), xt_stripe(5125, 1, 6), xt_stripe(5159, 2, 2), &
      xt_stripe(5164, 1, 2), xt_stripe(5175, 3, 1), xt_stripe(5178, 1, 6), &
      xt_stripe(5185, 1, 2), xt_stripe(5212, 1, 2), xt_stripe(5218, 1, 4), &
      xt_stripe(5230, 2, 4), xt_stripe(5237, 3, 1), xt_stripe(5240, 1, 4), &
      xt_stripe(5321, 6, 1), xt_stripe(5327, 1, 3), xt_stripe(5483, 5, 2), &
      xt_stripe(5489, 6, 1), xt_stripe(5495, 1, 3), xt_stripe(5499, 2, 1), &
      xt_stripe(5501, 1, 6), xt_stripe(5535, 2, 2), xt_stripe(5540, 1, 2), &
      xt_stripe(5551, 3, 1), xt_stripe(5554, 1, 6), xt_stripe(5561, 1, 2), &
      xt_stripe(5588, 1, 2), xt_stripe(5594, 1, 4), xt_stripe(5606, 2, 4), &
      xt_stripe(5613, 3, 1), xt_stripe(5616, 1, 4), xt_stripe(5697, 6, 1), &
      xt_stripe(5703, 1, 3) /)
    TYPE(xt_idxlist) :: idxlist

    idxlist = xt_idxvec_new(index_vector)
    CALL check_idxlist_stripes_pos_ext(idxlist, stripes)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_idxlist_stripes_pos_ext3

#if SIZEOF_XT_INT > 2
  SUBROUTINE test_idxlist_stripes_pos_ext4
    INTEGER, PARAMETER :: num_indices = 3
    INTEGER(xt_int_kind), PARAMETER :: index_vector(num_indices) &
         = (/ 328669_xi, 30608_xi, 38403_xi /)
    INTEGER, PARAMETER :: num_stripes = 1
    TYPE(xt_stripe), PARAMETER :: stripes(num_stripes) = (/ &
      xt_stripe(30608_xi, 7795_xi, 2)/)
    TYPE(xt_idxlist) :: idxlist

    idxlist = xt_idxvec_new(index_vector)
    CALL check_idxlist_stripes_pos_ext(idxlist, stripes)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_idxlist_stripes_pos_ext4

  SUBROUTINE test_idxlist_stripes_pos_ext5
    INTEGER, PARAMETER :: num_indices = 3
    INTEGER(xt_int_kind), PARAMETER :: index_vector(num_indices) &
         = (/ 679605_xi, 726349_xi, 726346_xi /)
    INTEGER, PARAMETER :: num_stripes = 1
    TYPE(xt_stripe), PARAMETER :: stripes(num_stripes) = (/ &
      xt_stripe(679605_xi, 46741_xi, 2)/)
    TYPE(xt_idxlist) :: idxlist

    idxlist = xt_idxvec_new(index_vector)
    CALL check_idxlist_stripes_pos_ext(idxlist, stripes)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_idxlist_stripes_pos_ext5
#endif

  SUBROUTINE test_idxlist_stripes_pos_ext_randomized1(full_random)
    LOGICAL, INTENT(in) :: full_random
    INTEGER, PARAMETER :: num_iterations=128, &
         max_num_indices=1024, max_index=1024

    INTEGER :: i, iteration, num_indices
    INTEGER(xt_int_kind), ALLOCATABLE :: indices(:)
    REAL, ALLOCATABLE :: rvals(:)
    TYPE(xt_idxlist) :: idxlist
    TYPE(xt_stripe), ALLOCATABLE :: stripes(:)
    TYPE(xt_stripe) :: stripes_dummy(1)

    CALL init_fortran_random(full_random)
    ALLOCATE(indices(max_num_indices), rvals(max_num_indices))
    DO iteration = 1, num_iterations
      CALL random_number(rvals(1))
      num_indices = NINT(rvals(1) * REAL(max_num_indices))

      CALL random_number(rvals(1:num_indices))
      DO i = 1, num_indices
        indices(i) = NINT(rvals(i)*REAL((2*max_index)-max_index), xt_int_kind)
      END DO
      idxlist = xt_idxvec_new(indices(1:num_indices))

      CALL xt_idxlist_get_index_stripes(idxlist, stripes)
      IF (ALLOCATED(stripes) .EQV. num_indices == 0) &
         CALL test_abort("get index stripes returned values for empty list", &
         filename, __LINE__)
      IF (num_indices > 0) THEN
        CALL check_idxlist_stripes_pos_ext(idxlist, stripes)
      ELSE
        CALL check_idxlist_stripes_pos_ext(idxlist, stripes_dummy(1:0))
      END IF

      CALL xt_idxlist_delete(idxlist)
    END DO
  END SUBROUTINE test_idxlist_stripes_pos_ext_randomized1

  SUBROUTINE check_idxlist_stripes_pos_ext(idxlist, stripes)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    TYPE(xt_stripe), INTENT(in) :: stripes(:)

    TYPE(xt_pos_ext), ALLOCATABLE :: pos_ext(:)
    INTEGER :: num_stripes, num_ext, num_unmatched
    INTEGER :: abs_pos_ext_size, jsign, i, j, k, send_pos
    INTEGER(xt_int_kind) :: intersection_index, orig_index
    LOGICAL, PARAMETER :: single_match_only = .TRUE.
    LOGICAL :: unmatched_in_intersection, unmatched_in_idxlist
    TYPE(xt_idxlist) :: intersection
    num_stripes = SIZE(stripes)

    num_unmatched = xt_idxlist_get_pos_exts_of_index_stripes( &
         idxlist, num_stripes, stripes, num_ext, pos_ext, single_match_only)

    ! testing of results
    IF (num_unmatched /= 0) &
         CALL test_abort("error in xt_idxlist_get_pos_exts_of_index_stripes", &
         filename, __LINE__)
    intersection = xt_idxvec_from_stripes_new(stripes)
    k = 0
    DO i = 1, num_ext
      abs_pos_ext_size = INT(ABS(pos_ext(i)%size))
      jsign = MERGE(1, -1, pos_ext(i)%size >= 0)
      DO j = 0, abs_pos_ext_size-1
        unmatched_in_intersection &
             = xt_idxlist_get_index_at_position(intersection, k, &
             intersection_index)
        send_pos = pos_ext(i)%start + jsign * j
        unmatched_in_idxlist &
             = xt_idxlist_get_index_at_position(idxlist, send_pos, orig_index)
        IF (unmatched_in_intersection .OR. unmatched_in_idxlist &
             .OR. intersection_index /= orig_index) THEN
          WRITE (0, '(4(a,i0))') "intersection pos ", k, &
               " index ", intersection_index, &
               " orig pos ", send_pos, &
               " index ", orig_index
          CALL test_abort("error in xt_idxlist_get_pos_exts_of_index_stripes", &
               filename, __LINE__)
        END IF
        k = k + 1
      END DO
    END DO
    CALL xt_idxlist_delete(intersection)
  END SUBROUTINE check_idxlist_stripes_pos_ext

  SUBROUTINE test_get_pos(stripes, pos)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER, INTENT(in) :: pos(:)
    INTEGER(xt_int_kind), PARAMETER :: dummy = 1_xi
    INTEGER(xt_int_kind) :: ref_sel_idx(SIZE(pos)), sel_idx(SIZE(pos))
    INTEGER(xt_int_kind), PARAMETER :: undef_idx = -HUGE(dummy)
    INTEGER :: num_pos, ip, p, ref_undef_count, undef_count
    TYPE(xt_idxlist) :: idxlist
    idxlist = xt_idxstripes_new(stripes)
    num_pos = SIZE(pos)
    ref_undef_count = 0
    DO ip = 1, num_pos
      p = pos(ip)
      IF (xt_idxlist_get_index_at_position(idxlist, p, ref_sel_idx(ip))) THEN
        ref_sel_idx(ip) = undef_idx
        ref_undef_count = ref_undef_count + 1
      END IF
    END DO
    undef_count = xt_idxlist_get_indices_at_positions(idxlist, pos, sel_idx, &
         undef_idx)
    IF (undef_count /= ref_undef_count) &
         CALL test_abort("inequal undef count!", filename, __LINE__)
    IF (ANY(sel_idx /= ref_sel_idx)) &
         CALL test_abort("incorrect index returned for position!", &
         filename, __LINE__)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_get_pos

  SUBROUTINE test_get_pos1
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 1, 5), xt_stripe(20, -1, 5) /)
    INTEGER, PARAMETER :: pos(13) = &
         (/   0,   2,   7,   9,  11, &
         &  100,  11, 200,   9, 300, &
         &   18, 400,   5 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos1

  SUBROUTINE test_get_pos2
    TYPE(xt_stripe), PARAMETER :: stripes(4) = (/ xt_stripe(0, 1, 3), &
         xt_stripe(10, 1, 2), xt_stripe(20, -1, 6), xt_stripe(30, -1, 7) /)
    INTEGER, PARAMETER :: pos(19) = &
         (/   -1,    0,    1,    2,    3,    4,   23,    5,    6,    7, &
         &     8,    9,   10,   11,   12,    0,    2,  100, 2000 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos2

  SUBROUTINE test_get_pos3
    TYPE(xt_stripe), PARAMETER :: stripes(4) = (/ xt_stripe(0, 1, 3), &
         xt_stripe(10, 1, 2), xt_stripe(20, -1, 6), xt_stripe(30, -1, 7) /)
    INTEGER, PARAMETER :: pos(13) = &
         (/    4,    7,    2,    5,    9,    0,   10,    6,   11,    8, &
         &    12,    1,    3 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos3

  SUBROUTINE test_get_pos4
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 1, 5), xt_stripe(20, -1, 5) /)
    INTEGER, PARAMETER :: pos(7) = &
         (/  -10,  200,  700,   90,   90,   18,  141 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos4

  SUBROUTINE test_stripe_overlap
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(1, 1, 5) /)
#ifndef __G95__
    INTEGER(xi) :: i, j
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ ((i + j, i=0,4), j = 0, 1) /)
#else
    INTEGER :: i, j
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ ((INT(i + j, xi), i=0,4), j = 0, 1) /)
#endif
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE test_stripe_overlap

  SUBROUTINE test_stripe_bb(stripes, global_size, global_start_index, bounds_ref)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER(xt_int_kind), INTENT(in) :: global_size(:), global_start_index
    TYPE(xt_bounds), INTENT(in) :: bounds_ref(:)

    TYPE(xt_bounds) :: bounds(SIZE(global_size))
    TYPE(xt_idxlist) :: idxstripes

    IF (SIZE(global_size) /= SIZE(bounds_ref)) &
         CALL test_abort("size mismatch for bounding-box", filename, __LINE__)
    idxstripes = xt_idxstripes_new(stripes, SIZE(stripes))

    bounds = xt_idxlist_get_bounding_box(idxstripes, global_size, &
         global_start_index)
    IF (ANY(bounds /= bounds_ref)) &
         CALL test_abort("boundary box doesn't match reference", &
         filename, __LINE__)
    CALL xt_idxlist_delete(idxstripes)
  END SUBROUTINE test_stripe_bb

  SUBROUTINE test_stripe_bb1
    TYPE(xt_stripe), PARAMETER :: stripes(1) = (/ xt_stripe(-1, -1, -1) /)
    INTEGER(xt_int_kind), PARAMETER :: global_size(3) = 4_xi, &
         global_start_index = 0
    TYPE(xt_bounds), PARAMETER :: bounds_ref(3) = xt_bounds(0, 0)
    CALL test_stripe_bb(stripes(1:0), global_size, global_start_index, bounds_ref)
  END SUBROUTINE test_stripe_bb1

  SUBROUTINE test_stripe_bb2
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(47, -12, 2), &
         xt_stripe(32, 12, 2), xt_stripe(36, 12, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: global_size(3) = (/ 5_xi, 4_xi, 3_xi /), &
         global_start_index = 1
    TYPE(xt_bounds), PARAMETER :: bounds_ref(3) = (/ xt_bounds(2, 2), &
         xt_bounds(2, 2), xt_bounds(1, 2) /)
    CALL test_stripe_bb(stripes, global_size, global_start_index, bounds_ref)
  END SUBROUTINE test_stripe_bb2

  SUBROUTINE do_tests(idxlist, ref_indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)

    TYPE(xt_stripe), ALLOCATABLE :: stripes(:)
    TYPE(xt_stripe), PARAMETER :: dummy(1) = (/ xt_stripe(0,0,0) /)
    INTEGER :: num_stripes
    TYPE(xt_idxlist) :: temp_idxlist, idxlist_copy

    CALL check_idxlist(idxlist, ref_indices)
    CALL xt_idxlist_get_index_stripes(idxlist, stripes)
    IF (ALLOCATED(stripes)) THEN
      num_stripes = SIZE(stripes)
      temp_idxlist = xt_idxvec_from_stripes_new(stripes, num_stripes)
    ELSE
      num_stripes = 0
      temp_idxlist = xt_idxvec_from_stripes_new(dummy, num_stripes)
    END IF
    CALL check_idxlist(temp_idxlist, ref_indices)

    CALL xt_idxlist_delete(temp_idxlist)

    IF (ALLOCATED(stripes)) DEALLOCATE(stripes)

    ! test packing and unpacking
    idxlist_copy = idxlist_pack_unpack_copy(idxlist)

    ! check copy
    CALL check_idxlist(idxlist_copy, ref_indices)

    CALL xt_idxlist_delete(idxlist_copy)

    ! test copying
    idxlist_copy = xt_idxlist_copy(idxlist)

    ! check copy
    CALL check_idxlist(idxlist_copy, ref_indices)

    ! clean up
    CALL xt_idxlist_delete(idxlist_copy)
  END SUBROUTINE do_tests

  SUBROUTINE check_pos_ext(stripes, search_stripes, ref_pos_ext, &
       single_match_only, ref_unmatched, test_desc)
    TYPE(xt_stripe), INTENT(in) :: stripes(:), search_stripes(:)
    TYPE(xt_pos_ext), intent(in) :: ref_pos_ext(:)
    LOGICAL, INTENT(in) :: single_match_only
    INTEGER, INTENT(in) :: ref_unmatched
    CHARACTER(len=*) :: test_desc

    INTEGER :: num_search_stripes, num_ref_pos_ext, num_ext, &
         unmatched
    TYPE(xt_idxlist) :: idxstripes
    TYPE(xt_pos_ext), ALLOCATABLE :: pos_ext(:)

    num_search_stripes = SIZE(search_stripes)
    num_ref_pos_ext = SIZE(ref_pos_ext)

    idxstripes = xt_idxstripes_new(stripes)
    unmatched = xt_idxlist_get_pos_exts_of_index_stripes(idxstripes, &
         num_search_stripes, search_stripes, &
         num_ext, pos_ext, single_match_only)
    IF (unmatched /= ref_unmatched) &
         CALL test_abort("error in number of unmatched indices for " &
         // test_desc, filename, __LINE__)
    IF (num_ext < 0 .OR. num_ext /= num_ref_pos_ext) &
         CALL test_abort("error finding " // test_desc, filename, __LINE__)
    IF (ANY(pos_ext /= ref_pos_ext)) &
         CALL test_abort("incorrect position extent length found in "&
         // test_desc, filename, __LINE__)
    DEALLOCATE(pos_ext)
    CALL xt_idxlist_delete(idxstripes)
  END SUBROUTINE check_pos_ext

  SUBROUTINE check_pos_ext1
    INTEGER, PARAMETER :: num_stripes = 1, num_ref_pos_ext = 1, &
         num_ref_unmatched = 0

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10) /), &
         search_stripes(1) = (/ xt_stripe(10_xi, -1_xi, 5) /)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(9, -5) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "simple inverted stripe")
  END SUBROUTINE check_pos_ext1

  SUBROUTINE check_pos_ext2
    INTEGER, PARAMETER :: num_stripes = 1, num_ref_pos_ext = 1, &
         num_ref_unmatched = 5

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10) /), &
         search_stripes(2) = xt_stripe(10_xi, -1_xi, 5)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(9, -5) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "simple inverted stripe")
  END SUBROUTINE check_pos_ext2

  SUBROUTINE check_pos_ext3
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10) /), &
         search_stripes(1) = xt_stripe(10_xi, 1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(9, 2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search inc stripe over inc gap")
  END SUBROUTINE check_pos_ext3

  SUBROUTINE check_pos_ext4
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(25_xi, -1_xi, 11), xt_stripe(10_xi, -1_xi, 10) /), &
         search_stripes(1) = xt_stripe(10_xi, 1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(11, -2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search inc stripe over dec gap")
  END SUBROUTINE check_pos_ext4

  SUBROUTINE check_pos_ext5
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(25_xi, -1_xi, 11), xt_stripe(10_xi, -1_xi, 10) /), &
         search_stripes(1) = xt_stripe(15_xi, -1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(10, 2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over dec gap")
  END SUBROUTINE check_pos_ext5

  SUBROUTINE check_pos_ext6
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10) /), &
         search_stripes(1) = xt_stripe(15_xi, -1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(10, -2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over inc gap")
  END SUBROUTINE check_pos_ext6

  SUBROUTINE check_pos_ext7
    INTEGER, PARAMETER :: num_stripes = 3, num_ref_pos_ext = 1, &
         num_ref_unmatched = 8

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10), &
         &    xt_stripe(29_xi, 1_xi, 10) /), &
         search_stripes(1) = xt_stripe(32_xi, -1_xi, 30)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(23, -22) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over 2 inc gap")
  END SUBROUTINE check_pos_ext7

  SUBROUTINE check_pos_ext8
    INTEGER, PARAMETER :: num_stripes = 5, num_ref_pos_ext = 5, &
         num_ref_unmatched = 0

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10), &
         &    xt_stripe(29_xi, 1_xi, 10), xt_stripe(14_xi, -1_xi, 4), &
         &    xt_stripe(28_xi, -1_xi, 4) /), &
         search_stripes(1) = xt_stripe(32_xi, -1_xi, 30)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(23, -4), xt_pos_ext(34, 4), xt_pos_ext(19, -10), &
         &    xt_pos_ext(30, 4), xt_pos_ext(9, -8) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over jumbled stripes")
  END SUBROUTINE check_pos_ext8

END PROGRAM test_idxstripes_f
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
