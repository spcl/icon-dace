!>
!! @file test_idxvec_f.f90
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
PROGRAM test_idxvec
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_idxlist, xt_idxvec_new, &
       xt_idxvec_from_stripes_new, &
       xt_int_kind, xt_stripe, xt_idxlist_get_intersection, &
       xt_idxlist_get_index_stripes, xt_idxlist_delete, &
       xt_idxlist_get_positions_of_indices, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_position_of_index, xt_idxlist_get_position_of_index_off, &
       xt_bounds, xt_idxlist_get_bounding_box, xt_idxlist_copy
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: check_idxlist, test_err_count, check_stripes, &
       check_offsets, idxlist_pack_unpack_copy
  IMPLICIT NONE
  INTEGER, PARAMETER :: xi = xt_int_kind
  INTEGER(xt_int_kind), PARAMETER :: index_vector(7) &
       = (/ 1_xi, 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, 7_xi /)
  CHARACTER(len=*), PARAMETER :: filename = 'test_idxvec_f.f90'

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL test_idxvec_pack_unpack
  CALL test_copying
  CALL test_repeated_equal_indices
  CALL test_positions
  CALL test_intersection_surjective
  CALL test_intersection_partial
  CALL test_intersection_inverse_partial
  CALL test_intersection_unsort_partial
  CALL test_intersection_unsort_inverse_partial
  CALL test_stripes1
  CALL test_stripes2
  CALL test_stripes3
  CALL test_stripes4
  CALL test_stripes5
  CALL test_stripes6
  CALL test_stripes7
  CALL test_get_indices_at_positions1
  CALL test_get_indices_at_positions2
  CALL test_get_indices_at_positions3
  CALL test_get_positions_of_indices
  CALL test_bounding_box1
  CALL test_bounding_box2
  CALL test_bounding_box3

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE test_idxvec_pack_unpack
    TYPE(xt_idxlist) :: idxvector, idxvector_copy, intersection
    TYPE(xt_stripe), PARAMETER :: ref_stripes(1) = (/ xt_stripe(1, 1, 7) /)
    CALL setup_idxvec(idxvector, index_vector)

    idxvector_copy = idxlist_pack_unpack_copy(idxvector)

    CALL check_idxlist(idxvector_copy, index_vector)

    intersection = xt_idxlist_get_intersection(idxvector, idxvector_copy)

    CALL check_idxlist(intersection, index_vector)

    CALL compare_stripes(idxvector, ref_stripes)

    CALL xt_idxlist_delete(idxvector)
    CALL xt_idxlist_delete(idxvector_copy)
    CALL xt_idxlist_delete(intersection)

  END SUBROUTINE test_idxvec_pack_unpack

  SUBROUTINE test_copying
    TYPE(xt_idxlist) :: idxvector, idxvector_copy
    CALL setup_idxvec(idxvector, index_vector)
    idxvector_copy = xt_idxlist_copy(idxvector)
    CALL check_idxlist(idxvector_copy, index_vector)
    CALL xt_idxlist_delete(idxvector)
    CALL xt_idxlist_delete(idxvector_copy)
  END SUBROUTINE test_copying

  SUBROUTINE test_repeated_equal_indices
    INTEGER(xt_int_kind), PARAMETER :: index_vector(8) &
         = (/ 1_xi, 2_xi, 3_xi, 7_xi, 5_xi, 6_xi, 7_xi, 7_xi /)
    TYPE(xt_idxlist) :: idxvector
    CALL setup_idxvec(idxvector, index_vector)
    CALL xt_idxlist_delete(idxvector)
  END SUBROUTINE test_repeated_equal_indices

  SUBROUTINE test_positions
    LOGICAL, PARAMETER :: single_match_only = .TRUE.
    INTEGER(xt_int_kind), PARAMETER :: index_vector(20) &
         = (/ 10_xi, 15_xi, 14_xi, 13_xi, 12_xi, &
         &    15_xi, 10_xi, 11_xi, 12_xi, 13_xi, &
         &    23_xi, 18_xi, 19_xi, 20_xi, 21_xi, &
         &    31_xi, 26_xi, 27_xi, 28_xi, 29_xi /), &
         intersection_vector(13) &
         = (/ 12_xi, 12_xi, 13_xi, 13_xi, 14_xi, &
         &    15_xi, 15_xi, 20_xi, 21_xi, 23_xi, &
         &    28_xi, 29_xi, 31_xi /)
    INTEGER :: intersection_pos(SIZE(intersection_vector))
    INTEGER, PARAMETER :: ref_intersection_pos(SIZE(intersection_vector)) &
         = (/ 4, 8, 3, 9, 2, 1, 5, 13, 14, 10, 18, 19, 15 /)
    TYPE(xt_idxlist) :: idxvector
    INTEGER :: notfound
    CALL setup_idxvec(idxvector, index_vector)
    notfound = xt_idxlist_get_positions_of_indices(idxvector, &
         intersection_vector, intersection_pos, single_match_only)
    IF (notfound /= 0) &
         CALL test_abort('expected indices not found!', filename, __LINE__)
    CALL check_offsets(intersection_pos, ref_intersection_pos)
    CALL xt_idxlist_delete(idxvector)
  END SUBROUTINE test_positions

  SUBROUTINE test_intersection(index_vector_a, index_vector_b, &
       ref_intersection_indices)
    INTEGER(xt_int_kind), INTENT(in) :: index_vector_a(:), index_vector_b(:), &
         ref_intersection_indices(:)
    ! note: instead of declaring this as the move intuitive idxvector(2),
    ! two distinct variables are used to remain compatible with NAG 5.2
    TYPE(xt_idxlist) :: idxvector1, idxvector2, intersection
    CALL setup_idxvec(idxvector1, index_vector_a)
    CALL setup_idxvec(idxvector2, index_vector_b)
    intersection = xt_idxlist_get_intersection(idxvector1, idxvector2)
    CALL check_idxlist(intersection, ref_intersection_indices)
    CALL xt_idxlist_delete(intersection)
    CALL xt_idxlist_delete(idxvector2)
    CALL xt_idxlist_delete(idxvector1)
  END SUBROUTINE test_intersection

  SUBROUTINE test_intersection_surjective
    INTEGER(xt_int_kind), PARAMETER :: index_vector(3, 2) &
         = RESHAPE((/ 1_xi, 2_xi, 3_xi, 1_xi, 2_xi, 3_xi /), &
         &         SHAPE(index_vector)), &
         ref_intersection_indices(3) = (/ 1_xi, 2_xi, 3_xi /)
    CALL test_intersection(index_vector(:, 1), index_vector(:, 2), &
         ref_intersection_indices)
  END SUBROUTINE test_intersection_surjective

  SUBROUTINE test_intersection_partial
    INTEGER(xt_int_kind), PARAMETER :: index_vector(3, 2) &
         = RESHAPE((/ 1_xi, 2_xi, 3_xi, 2_xi, 3_xi, 4_xi /), &
         &         SHAPE(index_vector)), &
         ref_intersection_indices(2) = (/ 2_xi, 3_xi /)
    CALL test_intersection(index_vector(:, 1), index_vector(:, 2), &
         ref_intersection_indices)
  END SUBROUTINE test_intersection_partial

  SUBROUTINE test_intersection_inverse_partial
    INTEGER(xt_int_kind), PARAMETER :: index_vector(3, 2) &
         = RESHAPE((/ 2_xi, 3_xi, 4_xi, 1_xi, 2_xi, 3_xi /), &
         &         SHAPE(index_vector)), &
         ref_intersection_indices(2) = (/ 2_xi, 3_xi /)
    CALL test_intersection(index_vector(:, 1), index_vector(:, 2), &
         ref_intersection_indices)
  END SUBROUTINE test_intersection_inverse_partial

  SUBROUTINE test_intersection_unsort_partial
    INTEGER(xt_int_kind), PARAMETER :: index_vector(3, 2) &
         = RESHAPE((/ 4_xi, 2_xi, 3_xi, 3_xi, 1_xi, 2_xi /), &
         &         SHAPE(index_vector)), &
         ref_intersection_indices(2) = (/ 2_xi, 3_xi /)
    CALL test_intersection(index_vector(:, 1), index_vector(:, 2), &
         ref_intersection_indices)
  END SUBROUTINE test_intersection_unsort_partial

  SUBROUTINE test_intersection_unsort_inverse_partial
    INTEGER(xt_int_kind), PARAMETER :: index_vector(3, 2) &
         = RESHAPE((/ 3_xi, 1_xi, 2_xi, 4_xi, 2_xi, 3_xi /), &
         &         SHAPE(index_vector)), &
         ref_intersection_indices(2) = (/ 2_xi, 3_xi /)
    CALL test_intersection(index_vector(:, 1), index_vector(:, 2), &
         ref_intersection_indices)
  END SUBROUTINE test_intersection_unsort_inverse_partial

  SUBROUTINE test_idxvec_from_stripes(stripes, ref_indices)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)

    TYPE(xt_idxlist) :: idxvec
    idxvec = xt_idxvec_from_stripes_new(stripes)
    CALL check_idxlist(idxvec, ref_indices)
    CALL xt_idxlist_delete(idxvec)
  END SUBROUTINE test_idxvec_from_stripes

  SUBROUTINE test_stripes1
    TYPE(xt_stripe), PARAMETER :: stripes(2) = &
         (/ xt_stripe(5, 1, 5), xt_stripe(4, -1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 5_xi, 6_xi, 7_xi, 8_xi, 9_xi, 4_xi, 3_xi, 2_xi, 1_xi, 0_xi /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes1

  SUBROUTINE test_stripes2
    TYPE(xt_stripe), PARAMETER :: stripes(3) &
         = (/ xt_stripe(0, 1, 5), xt_stripe(2, 1, 5), xt_stripe(4, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi, &
         &    2_xi, 3_xi, 4_xi, 5_xi, 6_xi, &
         &    4_xi, 5_xi, 6_xi, 7_xi, 8_xi /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes2

  SUBROUTINE test_stripes3
    TYPE(xt_stripe), PARAMETER :: stripes(3) &
         = (/ xt_stripe(2, 1, 5), xt_stripe(0, 1, 5), xt_stripe(4, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, &
         &    0_xi, 1_xi, 2_xi, 3_xi, 4_xi, &
         &    4_xi, 5_xi, 6_xi, 7_xi, 8_xi /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes3

  SUBROUTINE test_stripes4
    TYPE(xt_stripe), PARAMETER :: stripes(3) &
         = (/ xt_stripe(2, 1, 5), xt_stripe(4, -1, 5), xt_stripe(4, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, &
         &    4_xi, 3_xi, 2_xi, 1_xi, 0_xi, &
         &    4_xi, 5_xi, 6_xi, 7_xi, 8_xi /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes4

  SUBROUTINE test_stripes5
    TYPE(xt_stripe), PARAMETER :: stripes(3) &
         = (/ xt_stripe(0, 3, 5), xt_stripe(1, 3, 5), xt_stripe(14, -3, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0_xi,  3_xi, 6_xi,  9_xi, 12_xi, &
         &    1_xi,  4_xi, 7_xi, 10_xi, 13_xi, &
         &   14_xi, 11_xi, 8_xi,  5_xi,  2_xi /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes5

  SUBROUTINE test_stripes6
    TYPE(xt_stripe), PARAMETER :: stripes(3) &
         = (/ xt_stripe(0, 3, 5), xt_stripe(2, 3, 5), xt_stripe(14, -3, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0_xi,  3_xi, 6_xi,  9_xi, 12_xi, &
         &    2_xi,  5_xi, 8_xi, 11_xi, 14_xi, &
         &   14_xi, 11_xi, 8_xi,  5_xi,  2_xi /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes6

  SUBROUTINE test_stripes7
    TYPE(xt_stripe), PARAMETER :: stripes(4) = (/ xt_stripe(0, -1, 5), &
         xt_stripe(1, 1, 5), xt_stripe(-5, -1, 5), xt_stripe(6, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(20) &
         = (/ 0_xi, -1_xi, -2_xi, -3_xi, -4_xi, &
         &    1_xi,  2_xi,  3_xi,  4_xi,  5_xi, &
         &   -5_xi, -6_xi, -7_xi, -8_xi, -9_xi, &
         &    6_xi,  7_xi,  8_xi,  9_xi, 10_xi /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes7

  SUBROUTINE test_get_indices_at_positions(indices, undef_idx, pos)
    INTEGER(xt_int_kind), INTENT(in) :: indices(:), undef_idx
    INTEGER, INTENT(in) :: pos(:)
    INTEGER(xt_int_kind) :: ref_sel_idx(SIZE(pos)), sel_idx(SIZE(pos))
    INTEGER :: i, num_pos, undef_count, ref_undef_count
    TYPE(xt_idxlist) :: idxvec
    CHARACTER(80) :: msg
    num_pos = SIZE(pos)
    idxvec = xt_idxvec_new(indices, SIZE(indices))
    ref_undef_count = 0
    DO i = 1, num_pos
      IF (xt_idxlist_get_index_at_position(idxvec, pos(i), ref_sel_idx(i))) &
           THEN
        ref_sel_idx(i) = undef_idx
        ref_undef_count = ref_undef_count + 1
      END IF
    END DO
    undef_count = xt_idxlist_get_indices_at_positions(idxvec, pos, sel_idx, &
         undef_idx)
    IF (undef_count /= ref_undef_count) THEN
      CALL test_abort("test_idxvec_f.f90: (undef_count /= ref_undef_count)", &
           filename, __LINE__)
    END IF
    DO i = 1, num_pos
      IF (sel_idx(i) /= ref_sel_idx(i)) THEN
        WRITE (msg, '(2(a,i0),a)') "test_idxvec_f.f90: (sel_idx(", i, &
             ") /= ref_sel_idx(", i, "))"
        CALL test_abort(msg, filename, __LINE__)
      END IF
    END DO
    CALL xt_idxlist_delete(idxvec)
  END SUBROUTINE test_get_indices_at_positions

  SUBROUTINE test_get_indices_at_positions1
    INTEGER(xt_int_kind), PARAMETER :: indices(16) &
         = (/ 0_xi,  3_xi,  6_xi,  9_xi, 12_xi, 1_xi,  4_xi, 7_xi, &
         &   10_xi, 13_xi, 14_xi, 11_xi,  8_xi, 5_xi,  2_xi, 1_xi /)
    INTEGER(xt_int_kind), PARAMETER :: undef_idx = -HUGE(undef_idx)
    INTEGER, PARAMETER :: pos(13) &
         = (/ 0, 2, 7, 9, 11, 100, 11, 200, 9, 300, 7, 400, 5 /)
    CALL test_get_indices_at_positions(indices, undef_idx, pos)
  END SUBROUTINE test_get_indices_at_positions1

  SUBROUTINE test_get_indices_at_positions2
    INTEGER(xt_int_kind), PARAMETER :: indices(16) &
         = (/ 0_xi,  3_xi,  6_xi,  9_xi, 12_xi, 1_xi, 4_xi, 7_xi, &
         &   10_xi, 13_xi, 14_xi, 11_xi,  8_xi, 5_xi, 2_xi, 1_xi /)
    INTEGER(xt_int_kind), PARAMETER :: undef_idx = -HUGE(undef_idx)
    INTEGER, PARAMETER :: pos(9) &
         = (/ 0, 2, 7, 9, 11, 11, 9, 7, 5 /)
    CALL test_get_indices_at_positions(indices, undef_idx, pos)
  END SUBROUTINE test_get_indices_at_positions2

  SUBROUTINE test_get_indices_at_positions3
    INTEGER(xt_int_kind), PARAMETER :: indices(16) &
         = (/ 0_xi,  3_xi,  6_xi,  9_xi, 12_xi, 1_xi, 4_xi, 7_xi, &
         &   10_xi, 13_xi, 14_xi, 11_xi,  8_xi, 5_xi, 2_xi, 1_xi /)
    INTEGER(xt_int_kind), PARAMETER :: undef_idx = -HUGE(indices(1))
    INTEGER, PARAMETER :: pos(9) &
         = (/ 100, 102, 107, 109, 1011, 1011, 109, 107, 105 /)
    CALL test_get_indices_at_positions(indices, undef_idx, pos)
  END SUBROUTINE test_get_indices_at_positions3

  SUBROUTINE test_get_positions_of_indices
    INTEGER(xt_int_kind), PARAMETER :: indices(2) = (/ 0_xi, 2_xi /), &
         selection(3) = (/ 1_xi, 2_xi, 3_xi /)
    INTEGER :: i, position, positions(SIZE(selection))
    INTEGER, PARAMETER :: ref_positions(3) = (/ -1, 1, -1 /)
    TYPE(xt_idxlist) :: idxvec

    idxvec = xt_idxvec_new(indices, SIZE(indices))
    IF (.NOT. xt_idxlist_get_position_of_index(idxvec, 1_xt_int_kind, &
         position)) &
         CALL test_abort('xt_idxlist_get_position_of_index did not return &
         &an error', &
         filename, __LINE__)
    IF (.NOT. xt_idxlist_get_position_of_index_off(idxvec, 1_xt_int_kind, &
         position, 0)) &
         CALL test_abort('xt_idxlist_get_position_of_index_off did not return &
         &an error', &
         filename, __LINE__)
    IF (.NOT. xt_idxlist_get_position_of_index_off(idxvec, 0_xt_int_kind, &
         position, 1)) &
         CALL test_abort('xt_idxlist_get_position_of_index_off did not return &
         &an error', &
         filename, __LINE__)
    IF (xt_idxlist_get_positions_of_indices(idxvec, selection, positions, &
         .FALSE.) /= 2) THEN
      CALL test_abort('xt_idxlist_get_positions_of_indices did not return &
           &correct number of matches', &
           filename, __LINE__)
    END IF
    DO i = 1, SIZE(selection)
      IF (positions(i) /= ref_positions(i)) &
           CALL test_abort('xt_idxlist_get_positions_of_indices returned &
           &incorrect position', &
           filename, __LINE__)
    END DO
    CALL xt_idxlist_delete(idxvec)
  END SUBROUTINE test_get_positions_of_indices

  SUBROUTINE test_bounding_box(indices, global_size, global_start_index, &
       ref_bounds)
    INTEGER(xt_int_kind), INTENT(in) :: indices(:), global_size(:), &
         global_start_index
    TYPE(xt_bounds), INTENT(in) :: ref_bounds(:)

    TYPE(xt_bounds) :: bounds(SIZE(global_size))
    INTEGER :: i, ndim
    TYPE(xt_idxlist) :: idxvec
    CHARACTER(80) :: msg

    ndim = SIZE(global_size)
    IF (SIZE(ref_bounds) /= ndim) &
         CALL test_abort('ERROR: inequal dimensions', filename, __LINE__)

    idxvec = xt_idxvec_new(indices, SIZE(indices))

    bounds = xt_idxlist_get_bounding_box(idxvec, global_size, &
         global_start_index)

    DO i = 1, ndim
      IF (bounds(i)%start /= ref_bounds(i)%start) THEN
        WRITE (0, '(2(a,i0))') 'bounds(', i, ')%start=', bounds(i)%start
        WRITE (0, '(2(a,i0))') 'ref_bounds(', i, ')%start=', ref_bounds(i)%start
        WRITE (msg, '(a,i0)') "ERROR: xt_idxlist_get_bounding_box inequal &
             &starts at i=", i
        CALL test_abort(msg, filename, __LINE__)
      END IF
      IF (bounds(i)%size /= ref_bounds(i)%size) THEN
        WRITE (0, '(2(a,i0))') 'bounds(', i, ')%size=', bounds(i)%size
        WRITE (0, '(2(a,i0))') 'ref_bounds(', i, ')%size=', ref_bounds(i)%size
        WRITE (msg, '(a,i0)') "ERROR: xt_idxlist_get_bounding_box inequal &
             &size at i=", i
        CALL test_abort(msg, filename, __LINE__)
      END IF
    END DO
    CALL xt_idxlist_delete(idxvec)
  END SUBROUTINE test_bounding_box

  SUBROUTINE test_bounding_box1
    INTEGER(xt_int_kind), PARAMETER :: indices(2) = (/ 21_xi, 42_xi /), &
         global_size(3) = 4_xi, global_start_index = 0_xi
    TYPE(xt_bounds), PARAMETER :: ref_bounds(3) = xt_bounds(1_xi, 2_xi)
    CALL test_bounding_box(indices, global_size, global_start_index, ref_bounds)
  END SUBROUTINE test_bounding_box1

  SUBROUTINE test_bounding_box2
    INTEGER(xt_int_kind), PARAMETER :: indices(5) &
         = (/ 45_xi, 35_xi, 32_xi, 48_xi, 33_xi /), &
         global_size(3) = (/ 5_xi, 4_xi, 3_xi /), global_start_index = 1_xi
    TYPE(xt_bounds), PARAMETER :: ref_bounds(3) &
         = (/ xt_bounds(2, 2), xt_bounds(2, 2), xt_bounds(1, 2) /)
    CALL test_bounding_box(indices, global_size, global_start_index, ref_bounds)
  END SUBROUTINE test_bounding_box2

  SUBROUTINE test_bounding_box3
    INTEGER(xt_int_kind), PARAMETER :: indices(1) = (/ -1_xi /), &
         global_size(3) = 4_xi, global_start_index = 0_xi
    TYPE(xt_bounds), PARAMETER :: ref_bounds(3) = xt_bounds(0, 0)
    CALL test_bounding_box(indices(1:0), global_size, global_start_index, &
         ref_bounds)
  END SUBROUTINE test_bounding_box3

  SUBROUTINE setup_idxvec(idxlist, index_vector)
    TYPE(xt_idxlist), INTENT(out) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: index_vector(:)

    idxlist = xt_idxvec_new(index_vector, SIZE(index_vector))
    CALL check_idxlist(idxlist, index_vector)
  END SUBROUTINE setup_idxvec

  SUBROUTINE compare_stripes(idxlist, ref_stripes)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    TYPE(xt_stripe), INTENT(in) :: ref_stripes(:)

    TYPE(xt_stripe), ALLOCATABLE :: stripes(:)

    CALL xt_idxlist_get_index_stripes(idxlist, stripes)

    CALL check_stripes(stripes, ref_stripes)

  END SUBROUTINE compare_stripes

END PROGRAM test_idxvec
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
