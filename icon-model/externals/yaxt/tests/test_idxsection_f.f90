!>
!! @file test_idxsection_f.f90
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
PROGRAM test_idxsection
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_idxlist, &
       xt_int_kind, xt_bounds, xt_stripe, &
       xt_idxlist_copy, xt_idxlist_get_index_stripes, xt_idxlist_delete, &
       xt_idxlist_get_intersection, xt_idxlist_get_positions_of_indices, &
       xt_idxlist_get_position_of_index, xt_idxlist_get_indices_const, &
       xt_idxlist_get_bounding_box, &
       xt_idxsection_new, xt_idxvec_new, OPERATOR(/=), OPERATOR(==)
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: check_idxlist, test_err_count, check_stripes, &
       idxlist_pack_unpack_copy
  IMPLICIT NONE
  INTEGER, PARAMETER :: xi = xt_int_kind
  CHARACTER(len=*), PARAMETER :: filename = 'test_idxsection_f.f90'

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL test_1d_section
  CALL test_2d_section
  CALL test_3d_section
  CALL test_4d_section
  CALL test_2d_simple
  CALL test_1d_intersection1
  CALL test_1d_intersection2
  CALL test_2d_intersection1
  CALL test_2d_1
  CALL test_2d_2
  CALL test_get_positions1
  CALL test_get_positions2
  CALL test_other_intersection
  CALL test_signed_sizes1
  CALL test_signed_sizes2
  CALL test_signed_size_positions
  CALL test_signed_size_intersections
  CALL test_section_with_stride1
  CALL test_section_with_stride2
  CALL test_bb1
  CALL test_bb2
  CALL test_bb3

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_finalize
  CALL finish_mpi

CONTAINS
  SUBROUTINE test_1d_section
    INTEGER(xt_int_kind), PARAMETER :: start = 0_xi
    INTEGER, PARAMETER :: num_dimensions = 1
    INTEGER(xt_int_kind), PARAMETER :: global_size(num_dimensions) = 10_xi, &
         local_start(num_dimensions) = 3_xi, &
         ref_indices(5) = (/ 3_xi, 4_xi, 5_xi, 6_xi, 7_xi /)
    INTEGER, PARAMETER :: local_size(num_dimensions) = 5
    TYPE(xt_stripe), PARAMETER :: ref_stripes(1) = xt_stripe(3, 1, 5)
    TYPE(xt_idxlist) :: idxsection

    ! create index section
    idxsection = xt_idxsection_new(start, global_size, local_size, local_start)

    ! testing
    CALL do_tests(idxsection, ref_indices, ref_stripes)

    ! clean up
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_1d_section

  SUBROUTINE test_2d_section
    INTEGER(xt_int_kind), PARAMETER :: start = 0_xi
    INTEGER, PARAMETER :: num_dimensions = 2
    INTEGER(xt_int_kind), PARAMETER :: global_size(num_dimensions) &
         = (/ 5_xi, 6_xi /), &
         local_start(num_dimensions) = (/ 1_xi, 2_xi /), &
         ref_indices(6) = (/ 8_xi, 9_xi, 14_xi, 15_xi, 20_xi, 21_xi /)
    INTEGER, PARAMETER :: local_size(num_dimensions) = (/ 3, 2 /)
    TYPE(xt_stripe), PARAMETER :: ref_stripes(3) = (/ xt_stripe(8, 1, 2), &
         xt_stripe(14, 1, 2), xt_stripe(20, 1, 2) /)
    TYPE(xt_idxlist) :: idxsection

    ! create index section
    idxsection = xt_idxsection_new(start, global_size, local_size, local_start)

    ! testing
    CALL do_tests(idxsection, ref_indices, ref_stripes)

    ! clean up
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_2d_section

  SUBROUTINE test_3d_section
    INTEGER(xt_int_kind), PARAMETER :: start = 0_xi
    INTEGER, PARAMETER :: num_dimensions = 3
    INTEGER(xt_int_kind), PARAMETER :: global_size(num_dimensions) = 4_xi, &
         local_start(num_dimensions) = (/ 0_xi, 1_xi, 1_xi /), &
         ref_indices(16) = (/ 5_xi, 6_xi, 9_xi, 10_xi, 21_xi, 22_xi, 25_xi, &
         26_xi, 37_xi, 38_xi, 41_xi, 42_xi, 53_xi, 54_xi, 57_xi, 58_xi /)
    INTEGER, PARAMETER :: local_size(num_dimensions) = (/ 4, 2, 2 /)
    TYPE(xt_stripe), PARAMETER :: ref_stripes(8) = (/ xt_stripe(5, 1, 2), &
         xt_stripe(9, 1, 2), xt_stripe(21, 1, 2), xt_stripe(25, 1, 2), &
         xt_stripe(37, 1, 2), xt_stripe(41, 1, 2), xt_stripe(53, 1, 2), &
         xt_stripe(57, 1, 2) /)
    TYPE(xt_idxlist) :: idxsection

    ! create index section
    idxsection = xt_idxsection_new(start, global_size, local_size, local_start)

    ! testing
    CALL do_tests(idxsection, ref_indices, ref_stripes)

    ! clean up
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_3d_section

  SUBROUTINE test_4d_section
    INTEGER(xt_int_kind), PARAMETER :: start = 0_xi
    INTEGER, PARAMETER :: num_dimensions = 4
    INTEGER(xt_int_kind) :: i, j, k, l
    INTEGER(xt_int_kind), PARAMETER :: global_size(num_dimensions) &
         = (/ 3_xi, 4_xi, 4_xi, 3_xi /), &
         local_start(num_dimensions) &
         = (/ 0_xi, 1_xi, 1_xi, 1_xi /), &
#ifdef __xlC__
         ref_indices(36) = &
              (/ 16_xi,17_xi,19_xi,20_xi,22_xi,23_xi, &
              28_xi,29_xi,31_xi,32_xi,34_xi,35_xi, &
              40_xi,41_xi,43_xi,44_xi,46_xi,47_xi, &
              64_xi,65_xi,67_xi,68_xi,70_xi,71_xi, &
              76_xi,77_xi,79_xi,80_xi,82_xi,83_xi, &
              88_xi,89_xi,91_xi,92_xi,94_xi,95_xi /)
#else
         ref_indices(36) &
         = (/ ((((16_xi + i + j*3_xi + k*12_xi + l*48_xi, &
         &        i=0_xi,1_xi), j=0_xi,2_xi), k=0_xi,2_xi), l=0_xi,1_xi) /)
#endif
    INTEGER, PARAMETER :: local_size(num_dimensions) = (/ 2, 3, 3, 2 /)
    TYPE(xt_stripe), PARAMETER :: ref_stripes(18) &
         = (/ (((xt_stripe(16_xi + j*3_xi + k*12_xi + l*48_xi, 1_xi, 2), &
         &       j = 0_xi,2_xi), k=0_xi,2_xi), l=0_xi,1_xi) /)
    TYPE(xt_idxlist) :: idxsection

    ! create index section
    idxsection = xt_idxsection_new(start, global_size, local_size, local_start)

    ! testing
    CALL do_tests(idxsection, ref_indices, ref_stripes)

    ! clean up
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_4d_section

  SUBROUTINE test_2d_simple
    INTEGER(xt_int_kind), PARAMETER :: start = 0_xi
    INTEGER, PARAMETER :: num_dimensions = 2
    INTEGER(xt_int_kind) :: i, j
    INTEGER(xt_int_kind), PARAMETER :: global_size(num_dimensions) &
         = (/ 5_xi, 10_xi /), &
         local_start(num_dimensions) = (/ 1_xi, 2_xi /), &
         ref_indices(12) &
         = (/ ((12_xi + i + j*10_xi, i=0_xi,3_xi), j=0_xi,2_xi) /)
    INTEGER, PARAMETER :: local_size(num_dimensions) = (/ 3, 4 /)
    TYPE(xt_idxlist) :: idxsection

    ! create index section
    idxsection = xt_idxsection_new(start, global_size, local_size, local_start)

    ! testing
    CALL check_idxlist(idxsection, ref_indices)

    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_2d_simple

  SUBROUTINE test_intersection(&
       start_a, global_size_a, local_size_a, local_start_a, &
       start_b, global_size_b, local_size_b, local_start_b, &
       ref_indices, ref_stripes)
    INTEGER(xt_int_kind), INTENT(in) :: start_a, start_b, global_size_a(:), &
         global_size_b(:), local_start_a(:), local_start_b(:), ref_indices(:)
    INTEGER, INTENT(in) :: local_size_a(:), local_size_b(:)
    TYPE(xt_stripe), INTENT(in) :: ref_stripes(:)
    TYPE(xt_idxlist) :: idxsection(2), intersection
    idxsection(1) = xt_idxsection_new(start_a, global_size_a, local_size_a, &
         local_start_a)
    idxsection(2) = xt_idxsection_new(start_b, global_size_b, local_size_b, &
         local_start_b)
    intersection = xt_idxlist_get_intersection(idxsection(1), idxsection(2))
    CALL xt_idxlist_delete(idxsection(2))
    CALL xt_idxlist_delete(idxsection(1))
    CALL do_tests(intersection, ref_indices, ref_stripes)
    CALL xt_idxlist_delete(intersection)
  END SUBROUTINE test_intersection

  SUBROUTINE test_1d_intersection1
    INTEGER(xt_int_kind), PARAMETER :: start=0_xi, &
         global_size_a(1) = 10_xi, global_size_b(1) = 15_xi, &
         local_start_a(1) = 4_xi, local_start_b(1) = 7_xi, &
         ref_indices(2) = (/ 7_xi, 8_xi /)
    INTEGER, PARAMETER :: local_size_a(1) = 5, local_size_b(1) = 6
    TYPE(xt_stripe), PARAMETER :: ref_stripes(1) = (/ xt_stripe(7, 1, 2) /)
    CALL test_intersection(&
         start, global_size_a, local_size_a, local_start_a, &
         start, global_size_b, local_size_b, local_start_b, &
         ref_indices, ref_stripes)
  END SUBROUTINE test_1d_intersection1

  SUBROUTINE test_1d_intersection2
    INTEGER(xt_int_kind), PARAMETER :: start=0_xi, global_size_a(1) = 10_xi, &
         global_size_b(1) = 10_xi, local_start_a(1) = 3_xi, &
         local_start_b(1) = 4_xi, ref_indices(1) = (/ -1_xi /)
    INTEGER, PARAMETER :: local_size_a(1) = 1, local_size_b(1) = 5
    TYPE(xt_stripe), PARAMETER :: ref_stripes(1) = (/ xt_stripe(-1, -1, -1) /)
    CALL test_intersection(&
         start, global_size_a, local_size_a, local_start_a, &
         start, global_size_b, local_size_b, local_start_b, &
         ref_indices(1:0), ref_stripes(1:0))
  END SUBROUTINE test_1d_intersection2

  SUBROUTINE test_2d_intersection1
    INTEGER, PARAMETER :: n = 2
    INTEGER(xt_int_kind), PARAMETER :: start=0_xi, &
         global_size_a(n) = 6_xi, global_size_b(n) = 6_xi, &
         local_start_a(n) = 1_xi, local_start_b(n) = (/ 3_xi, 2_xi /), &
         ref_indices(2) = (/ 20_xi, 26_xi /)
    INTEGER, PARAMETER :: local_size_a(n) = (/ 4, 2 /), local_size_b(n) = 3
    TYPE(xt_stripe), PARAMETER :: ref_stripes(2) = (/ xt_stripe(20, 1, 1), &
         xt_stripe(26, 1, 1) /)
    CALL test_intersection(&
         start, global_size_a, local_size_a, local_start_a, &
         start, global_size_b, local_size_b, local_start_b, &
         ref_indices, ref_stripes)
  END SUBROUTINE test_2d_intersection1

  SUBROUTINE test_2d_1
    INTEGER, PARAMETER :: n = 2
    INTEGER(xt_int_kind), PARAMETER :: start=0_xi, &
         global_size(n) = 4, local_start(n) = (/ 0_xi, 2_xi /), &
         ref_indices(4) = (/ 2_xi, 3_xi, 6_xi, 7_xi /)
    INTEGER, PARAMETER :: local_size(n) = 2
    TYPE(xt_idxlist) :: idxsection
    idxsection = xt_idxsection_new(start, n, global_size, local_size, &
         local_start)
    CALL check_idxlist(idxsection, ref_indices)
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_2d_1

  SUBROUTINE test_2d_2
    INTEGER, PARAMETER :: n = 2
    INTEGER(xt_int_kind), PARAMETER :: start=1_xi, &
         global_size(n) = 4_xi, local_start(n) = (/ 0_xi, 2_xi /), &
         ref_indices(4) = (/ 3_xi, 4_xi, 7_xi, 8_xi /)
    INTEGER, PARAMETER :: local_size(n) = 2
    TYPE(xt_idxlist) :: idxsection
    idxsection = xt_idxsection_new(start, n, global_size, local_size, &
         local_start)
    CALL check_idxlist(idxsection, ref_indices)
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_2d_2

  SUBROUTINE test_get_positions1
    INTEGER, PARAMETER :: n = 2, num_selection = 6
    INTEGER(xt_int_kind), PARAMETER :: start=0_xi, &
         global_size(n) = 4_xi, local_start(n) = (/ 0_xi, 2_xi /)
    INTEGER, PARAMETER :: local_size(n) = 2
    INTEGER(xt_int_kind), PARAMETER :: selection(num_selection) &
         = (/ 1_xi, 2_xi, 5_xi, 6_xi, 7_xi, 8_xi /)
    INTEGER, PARAMETER :: ref_positions(num_selection) &
         = (/ 1*0 - 1, 2*0 + 0, 5*0 - 1, 6*0 + 2, 7*0 + 3, 8*0 - 1 /)
    INTEGER :: positions(num_selection), num_found
    TYPE(xt_idxlist) :: idxsection
    idxsection = xt_idxsection_new(start, n, global_size, local_size, &
         local_start)
    num_found = xt_idxlist_get_positions_of_indices(idxsection, selection, &
         positions, .FALSE.)
    IF (num_found /= 3) &
         CALL test_abort("xt_idxlist_get_positions_of_indices &
         &returned incorrect num_unmatched", &
         filename, __LINE__)
    IF (ANY(positions /= ref_positions)) &
         CALL test_abort("xt_idxlist_get_positions_of_indices &
         &returned incorrect position", &
         filename, __LINE__)
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_get_positions1

  SUBROUTINE test_get_positions2
    INTEGER, PARAMETER :: n = 2, num_selection = 9
    INTEGER(xt_int_kind), PARAMETER :: start=0_xi, &
         global_size(n) = 4_xi, local_start(n) = (/ 0_xi, 2_xi /)
    INTEGER, PARAMETER :: local_size(n) = 2
    INTEGER(xt_int_kind), PARAMETER :: selection(num_selection) &
         = (/ 2_xi, 1_xi, 5_xi, 7_xi, 6_xi, 7_xi, 7_xi, 6_xi, 8_xi /)
    INTEGER, PARAMETER :: ref_positions(num_selection) &
         = (/ 2*0 + 0, 1*0 - 1, 5*0 - 1, 7*0 + 3, 6*0 + 2, 7*0 + 3, 7*0 + 3, &
         &    6*0 + 2, 8*0 - 1 /)
    INTEGER :: positions(num_selection), num_found, i, p
    LOGICAL :: notfound
    TYPE(xt_idxlist) :: idxsection
    idxsection = xt_idxsection_new(start, n, global_size, local_size, &
         local_start)
    num_found = xt_idxlist_get_positions_of_indices(idxsection, selection, &
         positions, .FALSE.)
    IF (num_found /= 3) &
         CALL test_abort("xt_idxlist_get_position_of_indices &
         &returned incorrect num_unmatched", &
         filename, __LINE__)
    IF (ANY(positions /= ref_positions)) &
         CALL test_abort("xt_idxlist_get_position_of_indices &
         &returned incorrect position", &
         filename, __LINE__)
    DO i = 1, num_selection
      notfound = xt_idxlist_get_position_of_index(idxsection, selection(i), p)
      IF (p /= ref_positions(i) &
           .OR. (notfound .AND. ref_positions(i) /= -1)) &
           CALL test_abort("xt_idxlist_get_position_of_index &
           &returned incorrect position", &
           filename, __LINE__)
    END DO
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_get_positions2

  SUBROUTINE test_other_intersection
    INTEGER, PARAMETER :: n = 2, num_sel_idx = 9
    INTEGER(xt_int_kind), PARAMETER :: start=0_xi, &
         global_size(n) = 4_xi, local_start(n) = (/ 0_xi, 2_xi /), &
         sel_idx(num_sel_idx) &
         = (/ 2_xi, 1_xi, 5_xi, 7_xi, 6_xi, 7_xi, 7_xi, 6_xi, 8_xi /), &
         ref_inter_idx(6) = (/ 2_xi, 6_xi, 6_xi, 7_xi, 7_xi, 7_xi /)
    INTEGER, PARAMETER :: local_size(n) = 2
    TYPE(xt_idxlist) :: idxsection, sel_idxlist, inter_idxlist

    idxsection = xt_idxsection_new(start, n, global_size, local_size, &
         local_start)
    sel_idxlist = xt_idxvec_new(sel_idx)
    inter_idxlist = xt_idxlist_get_intersection(idxsection, sel_idxlist)
    CALL xt_idxlist_delete(sel_idxlist)
    CALL xt_idxlist_delete(idxsection)
    CALL check_idxlist(inter_idxlist, ref_inter_idx)
    CALL xt_idxlist_delete(inter_idxlist)
  END SUBROUTINE test_other_intersection

  ! test 2D section with arbitrary size signs
  SUBROUTINE test_signed_sizes1
    INTEGER :: i
    TYPE(xt_idxlist) :: idxsection
    INTEGER, PARAMETER :: n = 2
    INTEGER(xt_int_kind), PARAMETER :: start = 0, &
         global_size(n, 4) = RESHAPE( &
         (/ 5_xi, 10_xi, 5_xi,-10_xi, -5_xi, 10_xi, -5_xi, -10_xi /), &
         (/ n, 4 /) ), &
         local_start(2) = (/ 1_xi, 2_xi /), &
         ref_indices(12, 16) = RESHAPE( &
         (/ 12_xi, 13_xi, 14_xi, 15_xi, 22_xi, 23_xi, 24_xi, 25_xi, 32_xi, &
         &  33_xi, 34_xi, 35_xi, &
         &  15_xi, 14_xi, 13_xi, 12_xi, 25_xi, 24_xi, 23_xi, 22_xi, 35_xi, &
         &  34_xi, 33_xi, 32_xi, &
         &  32_xi, 33_xi, 34_xi, 35_xi, 22_xi, 23_xi, 24_xi, 25_xi, 12_xi, &
         &  13_xi, 14_xi, 15_xi, &
         &  35_xi, 34_xi, 33_xi, 32_xi, 25_xi, 24_xi, 23_xi, 22_xi, 15_xi, &
         &  14_xi, 13_xi, 12_xi, &
         &  17_xi, 16_xi, 15_xi, 14_xi, 27_xi, 26_xi, 25_xi, 24_xi, 37_xi, &
         &  36_xi, 35_xi, 34_xi, &
         &  14_xi, 15_xi, 16_xi, 17_xi, 24_xi, 25_xi, 26_xi, 27_xi, 34_xi, &
         &  35_xi, 36_xi, 37_xi, &
         &  37_xi, 36_xi, 35_xi, 34_xi, 27_xi, 26_xi, 25_xi, 24_xi, 17_xi, &
         &  16_xi, 15_xi, 14_xi, &
         &  34_xi, 35_xi, 36_xi, 37_xi, 24_xi, 25_xi, 26_xi, 27_xi, 14_xi, &
         &  15_xi, 16_xi, 17_xi, &
         &  32_xi, 33_xi, 34_xi, 35_xi, 22_xi, 23_xi, 24_xi, 25_xi, 12_xi, &
         &  13_xi, 14_xi, 15_xi, &
         &  35_xi, 34_xi, 33_xi, 32_xi, 25_xi, 24_xi, 23_xi, 22_xi, 15_xi, &
         &  14_xi, 13_xi, 12_xi, &
         &  12_xi, 13_xi, 14_xi, 15_xi, 22_xi, 23_xi, 24_xi, 25_xi, 32_xi, &
         &  33_xi, 34_xi, 35_xi, &
         &  15_xi, 14_xi, 13_xi, 12_xi, 25_xi, 24_xi, 23_xi, 22_xi, 35_xi, &
         &  34_xi, 33_xi, 32_xi, &
         &  37_xi, 36_xi, 35_xi, 34_xi, 27_xi, 26_xi, 25_xi, 24_xi, 17_xi, &
         &  16_xi, 15_xi, 14_xi, &
         &  34_xi, 35_xi, 36_xi, 37_xi, 24_xi, 25_xi, 26_xi, 27_xi, 14_xi, &
         &  15_xi, 16_xi, 17_xi, &
         &  17_xi, 16_xi, 15_xi, 14_xi, 27_xi, 26_xi, 25_xi, 24_xi, 37_xi, &
         &  36_xi, 35_xi, 34_xi, &
         &  14_xi, 15_xi, 16_xi, 17_xi, 24_xi, 25_xi, 26_xi, 27_xi, 34_xi, &
         &  35_xi, 36_xi, 37_xi /), (/ 12, 16 /) )
    INTEGER, PARAMETER :: local_size(n, 4) = RESHAPE( &
         (/ 3, 4, 3, -4, -3, 4, -3, -4 /), (/ n, 4 /) )
    ! iterate through all sign combinations of -/+ for local and global
    ! and for x and y, giving 2^2^2 combinations
    DO i = 0, 15
      ! create index section

      idxsection = xt_idxsection_new(start, n, global_size(:, i/4 + 1), &
           local_size(:, MOD(i, 4) + 1), local_start)

      ! testing
      CALL check_idxlist(idxsection, ref_indices(:, i + 1))

      ! clean up
      CALL xt_idxlist_delete(idxsection)
    END DO
  END SUBROUTINE test_signed_sizes1

  ! test 2D section with arbitrary size signs
  SUBROUTINE test_signed_sizes2
    INTEGER :: i
    TYPE(xt_idxlist) :: idxsection
    INTEGER, PARAMETER :: n = 2
    INTEGER(xt_int_kind), PARAMETER :: start = 0, &
         global_size(n, 4) = RESHAPE( &
         (/ 5_xi, 6_xi, 5_xi,-6_xi, -5_xi, 6_xi, -5_xi, -6_xi /), &
         (/ n, 4 /) ), &
         local_start(2) = (/ 1_xi, 2_xi /), &
         ref_indices(6, 16) = RESHAPE( &
         (/  8_xi,  9_xi, 10_xi, 14_xi, 15_xi, 16_xi, &
         &  10_xi,  9_xi,  8_xi, 16_xi, 15_xi, 14_xi, &
         &  14_xi, 15_xi, 16_xi,  8_xi,  9_xi, 10_xi, &
         &  16_xi, 15_xi, 14_xi, 10_xi,  9_xi,  8_xi, &
         &   9_xi,  8_xi,  7_xi, 15_xi, 14_xi, 13_xi, &
         &   7_xi,  8_xi,  9_xi, 13_xi, 14_xi, 15_xi, &
         &  15_xi, 14_xi, 13_xi,  9_xi,  8_xi,  7_xi, &
         &  13_xi, 14_xi, 15_xi,  7_xi,  8_xi,  9_xi, &
         &  20_xi, 21_xi, 22_xi, 14_xi, 15_xi, 16_xi, &
         &  22_xi, 21_xi, 20_xi, 16_xi, 15_xi, 14_xi, &
         &  14_xi, 15_xi, 16_xi, 20_xi, 21_xi, 22_xi, &
         &  16_xi, 15_xi, 14_xi, 22_xi, 21_xi, 20_xi, &
         &  21_xi, 20_xi, 19_xi, 15_xi, 14_xi, 13_xi, &
         &  19_xi, 20_xi, 21_xi, 13_xi, 14_xi, 15_xi, &
         &  15_xi, 14_xi, 13_xi, 21_xi, 20_xi, 19_xi, &
         &  13_xi, 14_xi, 15_xi, 19_xi, 20_xi, 21_xi /), (/ 6, 16 /) )
    ! iterate through all sign combinations of -/+ for local and global
    INTEGER, PARAMETER :: local_size(n, 4) = RESHAPE( &
         (/ 2, 3, 2, -3, -2, 3, -2, -3 /), (/ n, 4 /) )
    ! iterate through all sign combinations of -/+ for local and global
    ! and for x and y, giving 2^2^2 combinations
    DO i = 0, 15
      ! create index section

      idxsection = xt_idxsection_new(start, n, global_size(:, i/4 + 1), &
           local_size(:, MOD(i, 4) + 1), local_start)

      ! testing
      CALL check_idxlist(idxsection, ref_indices(:, i + 1))

      ! clean up
      CALL xt_idxlist_delete(idxsection)
    END DO
  END SUBROUTINE test_signed_sizes2

  SUBROUTINE test_signed_size_intersections
    INTEGER, PARAMETER :: n = 2
    INTEGER(xt_int_kind), PARAMETER :: start = 0, &
         global_size(n, 4) = RESHAPE( &
         (/ 5_xi, 10_xi, 5_xi,-10_xi, -5_xi, 10_xi, -5_xi, -10_xi /), &
         (/ n, 4 /) ), &
         local_start(2) = (/ 1_xi, 2_xi /), &
         indices(12, 16) = RESHAPE( &
         (/ 12_xi, 13_xi, 14_xi, 15_xi, 22_xi, 23_xi, 24_xi, 25_xi, 32_xi, &
         &  33_xi, 34_xi, 35_xi, &
         &  15_xi, 14_xi, 13_xi, 12_xi, 25_xi, 24_xi, 23_xi, 22_xi, 35_xi, &
         &  34_xi, 33_xi, 32_xi, &
         &  32_xi, 33_xi, 34_xi, 35_xi, 22_xi, 23_xi, 24_xi, 25_xi, 12_xi, &
         &  13_xi, 14_xi, 15_xi, &
         &  35_xi, 34_xi, 33_xi, 32_xi, 25_xi, 24_xi, 23_xi, 22_xi, 15_xi, &
         &  14_xi, 13_xi, 12_xi, &
         &  17_xi, 16_xi, 15_xi, 14_xi, 27_xi, 26_xi, 25_xi, 24_xi, 37_xi, &
         &  36_xi, 35_xi, 34_xi, &
         &  14_xi, 15_xi, 16_xi, 17_xi, 24_xi, 25_xi, 26_xi, 27_xi, 34_xi, &
         &  35_xi, 36_xi, 37_xi, &
         &  37_xi, 36_xi, 35_xi, 34_xi, 27_xi, 26_xi, 25_xi, 24_xi, 17_xi, &
         &  16_xi, 15_xi, 14_xi, &
         &  34_xi, 35_xi, 36_xi, 37_xi, 24_xi, 25_xi, 26_xi, 27_xi, 14_xi, &
         &  15_xi, 16_xi, 17_xi, &
         &  32_xi, 33_xi, 34_xi, 35_xi, 22_xi, 23_xi, 24_xi, 25_xi, 12_xi, &
         &  13_xi, 14_xi, 15_xi, &
         &  35_xi, 34_xi, 33_xi, 32_xi, 25_xi, 24_xi, 23_xi, 22_xi, 15_xi, &
         &  14_xi, 13_xi, 12_xi, &
         &  12_xi, 13_xi, 14_xi, 15_xi, 22_xi, 23_xi, 24_xi, 25_xi, 32_xi, &
         &  33_xi, 34_xi, 35_xi, &
         &  15_xi, 14_xi, 13_xi, 12_xi, 25_xi, 24_xi, 23_xi, 22_xi, 35_xi, &
         &  34_xi, 33_xi, 32_xi, &
         &  37_xi, 36_xi, 35_xi, 34_xi, 27_xi, 26_xi, 25_xi, 24_xi, 17_xi, &
         &  16_xi, 15_xi, 14_xi, &
         &  34_xi, 35_xi, 36_xi, 37_xi, 24_xi, 25_xi, 26_xi, 27_xi, 14_xi, &
         &  15_xi, 16_xi, 17_xi, &
         &  17_xi, 16_xi, 15_xi, 14_xi, 27_xi, 26_xi, 25_xi, 24_xi, 37_xi, &
         &  36_xi, 35_xi, 34_xi, &
         &  14_xi, 15_xi, 16_xi, 17_xi, 24_xi, 25_xi, 26_xi, 27_xi, 34_xi, &
         &  35_xi, 36_xi, 37_xi /), (/ 12, 16 /) )
    INTEGER, PARAMETER :: local_size(n, 4) = RESHAPE( &
         (/ 3, 4, 3, -4, -3, 4, -3, -4 /), (/ n, 4 /) )
    INTEGER :: i, j
    TYPE(xt_idxlist) :: idxsection_a, idxsection_b, &
         idxvec_a, idxvec_b, idxsection_intersection, &
         idxsection_intersection_other, idxvec_intersection

    DO i = 0, 15
      DO j = 0, 15
        ! create index section
        idxsection_a = xt_idxsection_new(start, n, global_size(:, i/4 + 1), &
             local_size(:, MOD(i, 4) + 1), local_start)
        idxsection_b = xt_idxsection_new(start, n, global_size(:, j/4 + 1), &
             local_size(:, MOD(j, 4) + 1), local_start)
        ! create reference index vectors
        idxvec_a = xt_idxvec_new(indices(:, i+1))
        idxvec_b = xt_idxvec_new(indices(:, j+1))

        ! testing
        idxsection_intersection = xt_idxlist_get_intersection(idxsection_a, &
             idxsection_b)
        idxsection_intersection_other &
             = xt_idxlist_get_intersection(idxsection_a, idxvec_b)
        idxvec_intersection = xt_idxlist_get_intersection(idxvec_a, idxvec_b)

        CALL check_idxlist(idxsection_intersection, &
             xt_idxlist_get_indices_const(idxvec_intersection))
        CALL check_idxlist(idxsection_intersection_other, &
             xt_idxlist_get_indices_const(idxvec_intersection))

        ! clean up

        CALL xt_idxlist_delete(idxsection_a)
        CALL xt_idxlist_delete(idxsection_b)
        CALL xt_idxlist_delete(idxvec_a)
        CALL xt_idxlist_delete(idxvec_b)
        CALL xt_idxlist_delete(idxsection_intersection)
        CALL xt_idxlist_delete(idxsection_intersection_other)
        CALL xt_idxlist_delete(idxvec_intersection)
      END DO
    END DO
  END SUBROUTINE test_signed_size_intersections

  SUBROUTINE test_signed_size_positions
    TYPE(xt_idxlist) :: idxsection
    INTEGER, PARAMETER :: n = 2, num_pos = 34
    INTEGER :: positions(num_pos)

    INTEGER(xt_int_kind), PARAMETER :: start = 0, &
         global_size(n) = (/ -5_xi, 6_xi /), &
         local_start(n) = (/ 1_xi, 2_xi /), &
         ref_indices(6) = (/ 16_xi, 15_xi, 14_xi, 22_xi, 21_xi, 20_xi /), &
         indices(num_pos) = &
         (/ -1_xi,  0_xi,  1_xi,  2_xi,  3_xi, &
         &   4_xi,  5_xi,  6_xi,  7_xi,  8_xi, &
         &   9_xi, 10_xi, 11_xi, 12_xi, 13_xi, &
         &  14_xi, 15_xi, 14_xi, 16_xi, 17_xi, &
         &  18_xi, 19_xi, 20_xi, 20_xi, 21_xi, &
         &  22_xi, 23_xi, 24_xi, 25_xi, 26_xi, &
         &  27_xi, 28_xi, 29_xi, 30_xi /)
    INTEGER, PARAMETER :: local_size(n) = (/ -2, -3 /)
    INTEGER, PARAMETER :: ref_positions(num_pos) = &
         (/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
         &  -1, -1, -1, -1, -1,  2,  1, -1,  0, -1, &
         &  -1, -1,  5, -1,  4,  3, -1, -1, -1, -1, &
         &  -1, -1, -1, -1 /)

    ! create index section
    idxsection = xt_idxsection_new(start, n, global_size, &
         local_size, local_start)

    ! testing
    CALL check_idxlist(idxsection, ref_indices)

    ! check get_positions_of_indices
    if (xt_idxlist_get_positions_of_indices(idxsection, indices, positions, &
         .TRUE.) /= 28) &
         CALL test_abort("error in xt_idxlist_get_positions_of_indices &
         &(wrong number of unmatched indices)", &
         filename, __LINE__)

    IF (ANY(ref_positions /= positions)) &
         call test_abort("error in xt_idxlist_get_positions_of_indices &
         &(wrong position)", &
         filename, __LINE__)

    ! clean up
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_signed_size_positions

  SUBROUTINE test_section_with_stride(start, global_size, local_size, &
       local_start, ref_indices)
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(:), &
         local_start(:), ref_indices(:)
    INTEGER, INTENT(in) :: local_size(:)
    TYPE(xt_idxlist) :: idxsection
    idxsection = xt_idxsection_new(start, SIZE(global_size), global_size, &
         local_size, local_start)
    CALL check_idxlist(idxsection, ref_indices)
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE test_section_with_stride

  SUBROUTINE test_section_with_stride1
    INTEGER, PARAMETER :: num_dimensions = 3
    INTEGER(xt_int_kind), PARAMETER :: start = 0_xi, &
         global_size(num_dimensions) = (/ 5_xi, 5_xi, 2_xi /), &
         local_start(num_dimensions) = (/ 2_xi, 0_xi, 1_xi /), &
         ref_indices(12) = &
         (/ 21_xi, 23_xi, 25_xi, 27_xi, &
         &  31_xi, 33_xi, 35_xi, 37_xi, &
         &  41_xi, 43_xi, 45_xi, 47_xi /)
    INTEGER, PARAMETER :: local_size(num_dimensions) = (/ 3, 4, 1 /)
    CALL test_section_with_stride(start, global_size, local_size, local_start, &
         ref_indices)
  END SUBROUTINE test_section_with_stride1

  SUBROUTINE test_section_with_stride2
    INTEGER, PARAMETER :: num_dimensions = 4
    INTEGER(xt_int_kind), PARAMETER :: start = 0_xi, &
         global_size(num_dimensions) = (/ 3_xi, 2_xi, 5_xi, 2_xi /), &
         local_start(num_dimensions) = (/ 0_xi, 1_xi, 1_xi, 0_xi /), &
         ref_indices(12) = &
         (/ 12_xi, 14_xi, 16_xi, 18_xi, &
         &  32_xi, 34_xi, 36_xi, 38_xi, &
         &  52_xi, 54_xi, 56_xi, 58_xi /)
    INTEGER, PARAMETER :: local_size(num_dimensions) = (/ 3, 1, 4, 1 /)
    CALL test_section_with_stride(start, global_size, local_size, local_start, &
         ref_indices)
  END SUBROUTINE test_section_with_stride2

  SUBROUTINE check_bb(start, global_size, local_size, &
       local_start, bb_start, global_bb_size, ref_bb)
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(:), local_start(:), &
         bb_start, global_bb_size(:)
    INTEGER, INTENT(in) :: local_size(:)
    TYPE(xt_bounds), INTENT(in) :: ref_bb(:)
    TYPE(xt_idxlist) :: idxsection
    TYPE(xt_bounds) :: bounds(SIZE(global_bb_size))
    idxsection = xt_idxsection_new(start, SIZE(global_size), global_size, &
         local_size, local_start)
    bounds = xt_idxlist_get_bounding_box(idxsection, global_bb_size, bb_start)
    IF (ANY(bounds /= ref_bb)) &
         CALL test_abort("bounding box mismatch", filename, __LINE__)
    CALL xt_idxlist_delete(idxsection)
  END SUBROUTINE check_bb

  SUBROUTINE test_bb1
    INTEGER, PARAMETER :: num_dimensions = 3
    INTEGER(xt_int_kind), PARAMETER :: start = 0_xi, &
         global_size(num_dimensions) = 4_xi, &
         local_start(num_dimensions) = (/ 2_xi, 0_xi, 1_xi /)
    INTEGER, PARAMETER :: local_size(num_dimensions) = 0
    TYPE(xt_bounds), PARAMETER :: ref_bb(num_dimensions) = xt_bounds(0, 0)
    CALL check_bb(start, global_size, local_size, local_start, start, &
         INT(global_size, xt_int_kind), ref_bb)
  END SUBROUTINE test_bb1

  SUBROUTINE test_bb2
    INTEGER, PARAMETER :: num_dimensions = 3
    INTEGER(xt_int_kind), PARAMETER :: start = 1_xi, &
         global_size(num_dimensions) = (/ 5_xi, 4_xi, 3_xi /), &
         local_start(num_dimensions) = (/ 2_xi, 2_xi, 1_xi /)
    INTEGER, PARAMETER :: local_size(num_dimensions) = 2
    TYPE(xt_bounds), PARAMETER :: ref_bb(num_dimensions) = &
         (/ xt_bounds(2, 2), xt_bounds(2, 2), xt_bounds(1, 2) /)
    CALL check_bb(start, global_size, local_size, local_start, start, &
         INT(global_size, xt_int_kind), ref_bb)
  END SUBROUTINE test_bb2

  SUBROUTINE test_bb3
    INTEGER, PARAMETER :: num_dimensions = 4, bb_ndim = 3
    INTEGER(xt_int_kind), PARAMETER :: start = 1_xi, &
         global_size(num_dimensions) = (/ 5_xi, 2_xi, 2_xi, 3_xi /), &
         local_start(num_dimensions) = (/ 2_xi, 0_xi, 1_xi, 1_xi /), &
         global_bb_size(bb_ndim) = (/ 5_xi, 4_xi, 3_xi /)
    INTEGER, PARAMETER :: local_size(num_dimensions) = (/ 2, 2, 1, 2 /)
    TYPE(xt_bounds), PARAMETER :: ref_bb(bb_ndim) = &
         (/ xt_bounds(2, 2), xt_bounds(1, 3), xt_bounds(1, 2) /)
    CALL check_bb(start, global_size, local_size, local_start, start, &
         global_bb_size, ref_bb)
  END SUBROUTINE test_bb3

  SUBROUTINE do_tests(idxlist, ref_indices, ref_stripes)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)
    TYPE(xt_stripe), OPTIONAL, INTENT(in) :: ref_stripes(:)

    TYPE(xt_stripe), ALLOCATABLE :: stripes(:)
    TYPE(xt_idxlist) :: idxlist_copy

    CALL check_idxlist(idxlist, ref_indices)
    IF (PRESENT(ref_stripes)) THEN
      CALL xt_idxlist_get_index_stripes(idxlist, stripes)
      IF (ALLOCATED(stripes)) THEN
        CALL check_stripes(stripes, ref_stripes)
        DEALLOCATE(stripes)
      ELSE
        IF (SIZE(ref_stripes) /= 0) &
             CALL test_abort("failed to reproduce stripes", filename, __LINE__)
      END IF
    END IF

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

END PROGRAM test_idxsection
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
