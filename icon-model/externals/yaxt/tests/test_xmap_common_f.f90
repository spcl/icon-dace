!>
!! @file test_xmap_common_f.f90
!! @brief generic Fortran xmap test procedures
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
MODULE test_xmap_common
  USE iso_c_binding, ONLY: c_int
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, &
       xt_stripe, xt_idxlist, xt_idxlist_delete, xt_idxvec_new, &
       xt_idxstripes_new, xi => xt_int_kind, &
       xt_xmap, xt_xmap_copy, xt_xmap_delete, &
       xt_xmap_get_num_destinations, xt_xmap_get_num_sources, &
       xt_xmap_get_destination_ranks, xt_xmap_get_source_ranks, &
       xt_mpi_comm_mark_exclusive
  IMPLICIT NONE
  PRIVATE
  INTEGER :: my_rank
  PUBLIC :: xmap_self_test_main, test_self_xmap_construct
  INTERFACE test_self_xmap_construct
    MODULE PROCEDURE test_self_xmap_construct_idxlist
    MODULE PROCEDURE test_self_xmap_construct_indices
    MODULE PROCEDURE test_self_xmap_construct_stripes
  END INTERFACE test_self_xmap_construct
  CHARACTER(len=*), PARAMETER :: filename = 'test_xmap_common_f.f90'
CONTAINS
  SUBROUTINE xmap_self_test_main(xmap_new)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER :: ierror, i, j
    INTEGER :: comms(2)
    INTEGER(xi), PARAMETER :: lsize(2) = (/ 7_xi, 1023_xi /)

    CALL init_mpi
    CALL xt_initialize(mpi_comm_world)
    CALL mpi_comm_rank(mpi_comm_world, my_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", filename, __LINE__)

    comms(1) = mpi_comm_world
    CALL mpi_comm_dup(mpi_comm_world, comms(2), ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", filename, __LINE__)
    CALL xt_mpi_comm_mark_exclusive(comms(2))

    DO i = 1, SIZE(comms)
      DO j = 1, 2
        CALL test_xmap1a(xmap_new, lsize(j), comms(i))
        CALL test_xmap1b(xmap_new, lsize(j), comms(i))
      END DO
      CALL test_xmap2(xmap_new, comms(i))
    END DO

    CALL mpi_comm_free(comms(2), ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", filename, __LINE__)

    IF (test_err_count() /= 0) &
         CALL test_abort("non-zero error count!", filename, __LINE__)
    CALL xt_finalize
    CALL finish_mpi
  END SUBROUTINE xmap_self_test_main

  SUBROUTINE shift_idx(idx, offset)
    INTEGER(xi), INTENT(inout) :: idx(:)
    INTEGER(xi), INTENT(in) :: offset
    INTEGER :: i
    DO i = 1, SIZE(idx)
      idx(i) = idx(i) + INT(my_rank, xi) * offset
    END DO
  END SUBROUTINE shift_idx

  SUBROUTINE assert_xmap_is_to_self(xmap)
    TYPE(xt_xmap) :: xmap
    INTEGER :: rank(1)
    IF (xt_xmap_get_num_destinations(xmap) /= 1) &
         CALL test_abort("error in xmap construction", filename, __LINE__)

    IF (xt_xmap_get_num_sources(xmap) /= 1) &
         CALL test_abort("error in xt_xmap_get_num_sources", filename, __LINE__)
    CALL xt_xmap_get_destination_ranks(xmap, rank)
    IF (rank(1) /= my_rank) &
         CALL test_abort("error in xt_xmap_get_destination_ranks", &
         filename, __LINE__)

    CALL xt_xmap_get_source_ranks(xmap, rank)
    IF (rank(1) /= my_rank) &
         CALL test_abort("error in xt_xmap_get_source_ranks", &
         filename, __LINE__)

  END SUBROUTINE assert_xmap_is_to_self

  SUBROUTINE test_self_xmap_construct_stripes(src_stripes, dst_stripes, &
       xmap_new, comm)
    TYPE(xt_stripe), INTENT(in) :: src_stripes(:), dst_stripes(:)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER, INTENT(in) :: comm
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    src_idxlist = xt_idxstripes_new(src_stripes)
    dst_idxlist = xt_idxstripes_new(dst_stripes)
    CALL test_self_xmap_construct(src_idxlist, dst_idxlist, xmap_new, comm)
  END SUBROUTINE test_self_xmap_construct_stripes

  SUBROUTINE test_self_xmap_construct_indices(src_indices, dst_indices, &
       xmap_new, comm)
    INTEGER(xi), INTENT(in) :: src_indices(:), dst_indices(:)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER, INTENT(in) :: comm
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    src_idxlist = xt_idxvec_new(src_indices)
    dst_idxlist = xt_idxvec_new(dst_indices)
    CALL test_self_xmap_construct(src_idxlist, dst_idxlist, xmap_new, comm)
  END SUBROUTINE test_self_xmap_construct_indices

  SUBROUTINE test_self_xmap_construct_idxlist(src_idxlist, dst_idxlist, &
       xmap_new, comm)
    TYPE(xt_idxlist), INTENT(inout) :: src_idxlist, dst_idxlist
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER, INTENT(in) :: comm

    TYPE(xt_xmap) :: xmap, xmap_copy

    xmap = xmap_new(src_idxlist, dst_idxlist, comm)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    CALL assert_xmap_is_to_self(xmap)
    xmap_copy = xt_xmap_copy(xmap)
    CALL assert_xmap_is_to_self(xmap_copy)

    CALL xt_xmap_delete(xmap)
    CALL xt_xmap_delete(xmap_copy)
  END SUBROUTINE test_self_xmap_construct_idxlist

  SUBROUTINE test_xmap1a(xmap_new, lsize, comm)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER(xi), INTENT(in) :: lsize
    INTEGER, INTENT(in) :: comm

    INTEGER(xt_int_kind) :: src_indices(lsize), dst_indices(lsize), i
    DO i = 1_xi, lsize
      src_indices(i) = i
    END DO
    CALL shift_idx(src_indices, lsize)
    DO i = 1_xi, lsize
      dst_indices(i) = lsize - i + 1_xi
    END DO
    CALL shift_idx(dst_indices, lsize)

    CALL test_self_xmap_construct(src_indices, dst_indices, &
         xmap_new, comm)
  END SUBROUTINE test_xmap1a

  SUBROUTINE test_xmap1b(xmap_new, lsize, comm)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER(xi), INTENT(in) :: lsize
    INTEGER, INTENT(in) :: comm

    TYPE(xt_stripe) :: src_stripe(1), dst_stripe(1)
    src_stripe(1) = xt_stripe(1_xi + INT(my_rank, xi) * lsize, &
         1_xi, INT(lsize, c_int))
    dst_stripe(1) = xt_stripe(INT(my_rank+1, xi) * lsize, &
         -1_xi, INT(lsize, c_int))
    CALL test_self_xmap_construct(src_stripe, dst_stripe, &
         xmap_new, comm)
  END SUBROUTINE test_xmap1b

  SUBROUTINE test_xmap2(xmap_new, comm)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER, INTENT(in) :: comm

    INTEGER(xi) :: src_index_list(14), dst_index_list(13)
    src_index_list = &
         (/ 5_xi, 67_xi, 4_xi, 5_xi, 13_xi, &
         &  9_xi,  2_xi, 1_xi, 0_xi, 96_xi, &
         & 13_xi, 12_xi, 1_xi, 3_xi /)
    dst_index_list = &
         (/ 5_xi, 4_xi, 3_xi, 96_xi, 1_xi, &
         &  5_xi, 4_xi, 5_xi,  4_xi, 3_xi, &
         & 13_xi, 2_xi, 1_xi /)
    CALL test_self_xmap_construct(src_index_list, dst_index_list, xmap_new, comm)
  END SUBROUTINE test_xmap2

END MODULE test_xmap_common
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
