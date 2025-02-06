!>
!! @file test_idxlist_collection_f.f90
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
PROGRAM test_idxlist_collection_f
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE mpi
  USE test_idxlist_utils, ONLY: check_idxlist, test_err_count, &
       idxlist_pack_unpack_copy, check_idxlist_copy
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, &
       xt_idxlist, xt_idxvec_new, &
       xt_idxlist_collection_new, xt_idxlist_delete, xt_stripe, &
       xt_idxlist_get_intersection, xt_idxsection_new, xt_idxstripes_new, &
       xt_bounds, xt_idxempty_new, xt_idxlist_get_bounding_box, OPERATOR(/=)
  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER :: filename = 'test_idxlist_collection_f.f90'
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL test_idxlist_collection_pack_unpack
  CALL test_idxlist_collection_copy
  CALL test_idxlist_collection_intersection
  CALL test_idxlist_collection_heterogeneous
  CALL test_bounding_box1
  CALL test_bounding_box2

  CALL xt_finalize
  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL finish_mpi

CONTAINS
  SUBROUTINE test_idxlist_collection_pack_unpack
    INTEGER, PARAMETER :: num_indices = 7, num_vec = 2
    INTEGER(xt_int_kind) :: i, j
    INTEGER(xt_int_kind), PARAMETER :: index_list(num_indices, num_vec) = &
         RESHAPE((/ ((INT(i, xt_int_kind), i = 1, num_indices), &
         &                                 j = 1, num_vec) /), &
         &       shape = (/ num_indices, num_vec /))
    TYPE(xt_idxlist) :: idxlists(num_vec), collectionlist, collectionlist_copy
    TYPE(xt_stripe), PARAMETER :: ref_stripes(num_vec) = xt_stripe(1, 1, 7)
    INTEGER :: k
    DO k = 1, num_vec
      idxlists(k) = xt_idxvec_new(index_list(:, k), num_indices)
    END DO
    collectionlist = xt_idxlist_collection_new(idxlists)
    CALL xt_idxlist_delete(idxlists)
    CALL check_idxlist(collectionlist, &
         RESHAPE(index_list, (/ SIZE(index_list) /)))
    collectionlist_copy = idxlist_pack_unpack_copy(collectionlist)
    CALL check_idxlist_copy(collectionlist, collectionlist_copy, &
         RESHAPE(index_list, (/ SIZE(index_list) /)), ref_stripes)
    CALL xt_idxlist_delete(collectionlist_copy)
    CALL xt_idxlist_delete(collectionlist)
  END SUBROUTINE test_idxlist_collection_pack_unpack

  SUBROUTINE test_idxlist_collection_copy
    INTEGER, PARAMETER :: num_indices = 7, num_vec = 2
    INTEGER(xt_int_kind) :: i, j
    INTEGER(xt_int_kind), PARAMETER :: index_list(num_indices, num_vec) = &
         RESHAPE((/ ((INT(num_indices - (j * num_indices + 1 - j - i) &
         &                * (2*j - 1), xt_int_kind), &
         &           i=1, num_indices), j=1,0,-1) /), &
         &       (/ num_indices, num_vec /))
    TYPE(xt_idxlist) :: idxlists(num_vec), collectionlist, collectionlist_copy
    TYPE(xt_stripe), PARAMETER :: ref_stripes(num_vec) &
         = (/ xt_stripe(1, 1, 7), xt_stripe(7, -1, 7) /)
    INTEGER :: k
    DO k = 1, num_vec
      idxlists(k) = xt_idxvec_new(index_list(:, k), num_indices)
    END DO
    collectionlist = xt_idxlist_collection_new(idxlists)
    CALL xt_idxlist_delete(idxlists)
    CALL check_idxlist(collectionlist, &
         RESHAPE(index_list, (/ SIZE(index_list) /)))
    collectionlist_copy = idxlist_pack_unpack_copy(collectionlist)
    CALL check_idxlist_copy(collectionlist, collectionlist_copy, &
         RESHAPE(index_list, (/ SIZE(index_list) /)), ref_stripes)
    CALL xt_idxlist_delete(collectionlist_copy)
    CALL xt_idxlist_delete(collectionlist)
  END SUBROUTINE test_idxlist_collection_copy

  SUBROUTINE test_idxlist_collection_intersection
    INTEGER, PARAMETER :: num_indices = 7, num_lists = 3
    INTEGER, PARAMETER :: xi = xt_int_kind
    INTEGER(xt_int_kind), PARAMETER :: index_list(num_indices, num_lists) &
         = RESHAPE((/ 1_xi, 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, &
         &            7_xi, 6_xi, 5_xi, 4_xi, 3_xi, 2_xi, 1_xi, &
         &            2_xi, 6_xi, 1_xi, 4_xi, 7_xi, 3_xi, 0_xi /),  &
         &         (/ num_indices,  num_lists /)), &
         sorted_index_list(SIZE(index_list)) &
         = (/ 0_xi, 1_xi, 1_xi, 1_xi, 2_xi, 2_xi, 2_xi, &
         &    3_xi, 3_xi, 3_xi, 4_xi, 4_xi, 4_xi, 5_xi, &
         &    5_xi, 6_xi, 6_xi, 6_xi, 7_xi, 7_xi, 7_xi /)
    TYPE(xt_idxlist) :: idxlists(num_lists), collectionlist, intersection, &
         ref_idxvec
    INTEGER :: i

    DO i = 1, 3
      idxlists(i) = xt_idxvec_new(index_list(:, i), num_indices)
    END DO
    collectionlist = xt_idxlist_collection_new(idxlists)
    DO i = 1, 3
      CALL xt_idxlist_delete(idxlists(i))
    END DO
    CALL check_idxlist(collectionlist, &
         RESHAPE(index_list, (/ SIZE(index_list) /)))
    ref_idxvec = xt_idxvec_new(RESHAPE(index_list, (/ SIZE(index_list) /)), &
         SIZE(index_list))
    intersection = xt_idxlist_get_intersection(ref_idxvec, collectionlist)
    CALL check_idxlist(intersection, sorted_index_list)
    CALL xt_idxlist_delete(intersection)
    intersection = xt_idxlist_get_intersection(collectionlist, ref_idxvec)
    CALL check_idxlist(intersection, sorted_index_list)
    CALL xt_idxlist_delete(intersection)
    CALL xt_idxlist_delete(ref_idxvec)
    CALL xt_idxlist_delete(collectionlist)

  END SUBROUTINE test_idxlist_collection_intersection

  SUBROUTINE test_idxlist_collection_heterogeneous
    INTEGER, PARAMETER :: num_indices = 6, num_lists = 3
    INTEGER, PARAMETER :: xi = xt_int_kind
    INTEGER(xt_int_kind), PARAMETER :: &
         index_list(num_indices) = (/ 1_xi, 3_xi, 5_xi, 7_xi, 9_xi, 11_xi /)
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(0, 2, 5), &
         xt_stripe(1, 2, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: local_start(2) = 2
    INTEGER(xt_int_kind), PARAMETER :: global_size(2) = (/ 10_xi, 10_xi /)
    INTEGER, PARAMETER :: local_size(2) = 5
    INTEGER, PARAMETER :: ref_size = num_indices + stripes(1)%nstrides &
         + stripes(2)%nstrides + local_size(1) * local_size(2)
    INTEGER(xt_int_kind), PARAMETER :: ref_index_list(ref_size) &
         = (/ 1_xi, 3_xi, 5_xi, 7_xi, 9_xi, 11_xi, &
         &    0_xi, 2_xi, 4_xi, 6_xi, 8_xi, 1_xi, 3_xi, 5_xi, 7_xi, 9_xi, &
         &    22_xi, 23_xi, 24_xi, 25_xi, 26_xi, &
         &    32_xi, 33_xi, 34_xi, 35_xi, 36_xi, &
         &    42_xi, 43_xi, 44_xi, 45_xi, 46_xi, &
         &    52_xi, 53_xi, 54_xi, 55_xi, 56_xi, &
         &    62_xi, 63_xi, 64_xi, 65_xi, 66_xi /)
    TYPE(xt_idxlist) :: idxlists(num_lists), collectionlist

    idxlists(1) = xt_idxvec_new(index_list, SIZE(index_list))
    idxlists(2) = xt_idxstripes_new(stripes, SIZE(stripes))
    idxlists(3) = xt_idxsection_new(0_xt_int_kind, global_size, local_size, &
         local_start)

    ! generate a collection index list
    collectionlist = xt_idxlist_collection_new(idxlists)

    CALL xt_idxlist_delete(idxlists)

    ! test generated collection list
    CALL check_idxlist(collectionlist, ref_index_list)

    CALL xt_idxlist_delete(collectionlist)

  END SUBROUTINE test_idxlist_collection_heterogeneous

  SUBROUTINE test_bounding_box1
    INTEGER, PARAMETER :: ndim=3, num_lists = 2
    INTEGER(xt_int_kind), PARAMETER :: global_size_bb(ndim) = 4, &
         global_start_index = 0
    TYPE(xt_idxlist) :: idxlists(num_lists), collectionlist
    TYPE(xt_bounds) :: bounds(ndim)
    INTEGER :: i

    DO i = 1, num_lists
      idxlists(i) = xt_idxempty_new()
    END DO
    collectionlist = xt_idxlist_collection_new(idxlists)
    CALL xt_idxlist_delete(idxlists)
    bounds = xt_idxlist_get_bounding_box(collectionlist, global_size_bb, &
         global_start_index)
    IF (ANY(bounds%size /= 0)) &
         CALL test_abort("ERROR: non-zero bounding box size", &
         filename, __LINE__)
    CALL xt_idxlist_delete(collectionlist)
  END SUBROUTINE test_bounding_box1

  SUBROUTINE test_bounding_box2
    INTEGER, PARAMETER :: ndim = 3, num_lists = 2, num_indices = 3
    INTEGER, PARAMETER :: xi = xt_int_kind
    INTEGER(xt_int_kind), PARAMETER :: indices(num_indices, num_lists) &
         = RESHAPE( (/ 45_xi, 35_xi, 32_xi, 32_xi, 48_xi, 33_xi /), &
         &          (/ num_indices, num_lists /)), &
         global_size(ndim) = (/ 5_xi, 4_xi, 3_xi /), &
         global_start_index = 1
    TYPE(xt_idxlist) :: idxlists(num_lists), collectionlist
    TYPE(xt_bounds) :: bounds(ndim)
    TYPE(xt_bounds), PARAMETER :: bounds_ref(ndim) = (/ xt_bounds(2, 2), &
         xt_bounds(2, 2), xt_bounds(1, 2) /)
    INTEGER :: i

    DO i = 1, num_lists
      idxlists(i) = xt_idxvec_new(indices(:, i), SIZE(indices, 1))
    END DO
    collectionlist = xt_idxlist_collection_new(idxlists)
    CALL xt_idxlist_delete(idxlists)

    bounds = xt_idxlist_get_bounding_box(collectionlist, global_size, &
         global_start_index)
    CALL xt_idxlist_delete(collectionlist)
    IF (ANY(bounds /= bounds_ref)) &
         CALL test_abort("ERROR: unexpected boundaries", filename, __LINE__)
  END SUBROUTINE test_bounding_box2

END PROGRAM test_idxlist_collection_f
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
