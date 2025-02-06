!>
!! @file test_xmap_common_parallel_f.f90
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
MODULE test_xmap_common_parallel
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, xt_stripe, &
       xi => xt_int_kind, xt_sort_int, &
       xt_idxlist, xt_idxlist_delete, xt_idxvec_new, &
       xt_idxstripes_new, xt_idxempty_new, &
       xt_xmap, xt_xmap_copy, xt_xmap_delete, &
       xt_xmap_get_num_destinations, xt_xmap_get_num_sources, &
       xt_xmap_get_destination_ranks, xt_xmap_get_source_ranks, &
       xt_xmap_get_max_dst_pos, xt_xmap_get_max_src_pos, &
       xt_xmap_update_positions, xt_xmap_spread
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xmap_parallel_test_main
  PUBLIC :: get_rank_range
  PUBLIC :: check_allgather_analog_xmap
  PUBLIC :: test_ring_1d
  PUBLIC :: test_ping_pong
  CHARACTER(len=*), PARAMETER :: filename = 'test_xmap_common_parallel_f.f90'
CONTAINS
  SUBROUTINE xmap_parallel_test_main(xmap_new)
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
    INTEGER :: comm, comm_rank, comm_size
    INTEGER :: ierror
    CALL init_mpi
    comm = mpi_comm_world
    CALL xt_initialize(comm)
    CALL mpi_comm_rank(comm, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_rank", filename, __LINE__)
    CALL mpi_comm_size(comm, comm_size, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_size", filename, __LINE__)
    IF (comm_size > HUGE(1_xi)) &
         CALL test_abort("number of ranks exceeds test limit", &
         filename, __LINE__)

    CALL test_allgather_analog(xmap_new, 1, comm)
    ! repeat test for large index list that will cause stripifying
    CALL test_allgather_analog(xmap_new, 1024, comm)
    IF (comm_size > 2) CALL test_ring_1d(xmap_new, comm)
    IF (comm_size == 2) CALL test_pair(xmap_new, comm)
    IF (comm_size > 1) CALL test_ping_pong(xmap_new, comm, 0, comm_size - 1)

    ! test maxpos implementation for xt_xmap_intersection
    CALL test_maxpos(xmap_new, comm, 5)
    ! test maxpos implementation for xt_xmap_intersection_ext
    CALL test_maxpos(xmap_new, comm, 501)

    IF (test_err_count() /= 0) &
         CALL test_abort("non-zero error count!", filename, __LINE__)
    CALL xt_finalize
    CALL finish_mpi
  END SUBROUTINE xmap_parallel_test_main

  SUBROUTINE get_rank_range(comm, is_inter, comm_rank, comm_size)
    INTEGER, INTENT(inout) :: comm
    INTEGER, INTENT(out) :: comm_rank, comm_size
    LOGICAL, INTENT(out) :: is_inter
    INTEGER :: ierror

    CALL mpi_comm_rank(comm, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_rank", filename, __LINE__)
    CALL mpi_comm_test_inter(comm, is_inter, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_test_inter", &
         filename, __LINE__)
    IF (is_inter) THEN
      CALL mpi_comm_remote_size(comm, comm_size, ierror)
    ELSE
      CALL mpi_comm_size(comm, comm_size, ierror)
    END IF
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_(remote)_size", &
         filename, __LINE__)
  END SUBROUTINE get_rank_range

  SUBROUTINE check_allgather_analog_xmap(xmap, comm)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(inout) :: comm
    INTEGER, ALLOCATABLE :: ranks(:)
    INTEGER(xt_int_kind) :: i
    INTEGER :: comm_rank, comm_size
    LOGICAL :: is_inter

    CALL get_rank_range(comm, is_inter, comm_rank, comm_size)
    IF (xt_xmap_get_num_destinations(xmap) /= INT(comm_size, xi)) &
         CALL test_abort("error in xmap construction", filename, __LINE__)

    IF (xt_xmap_get_num_sources(xmap) /= INT(comm_size, xi)) &
         CALL test_abort("error in xt_xmap_get_num_sources", &
         filename, __LINE__)

    ALLOCATE(ranks(comm_size))

    CALL xt_xmap_get_destination_ranks(xmap, ranks)
    IF (ANY(ranks /= (/ (i, i=0_xi,INT(comm_size-1, xi)) /))) &
         CALL test_abort("error in xt_xmap_get_destination_ranks", &
         filename, __LINE__)

    CALL xt_xmap_get_source_ranks(xmap, ranks)
    IF (ANY(ranks /= (/ (i, i=0_xi,INT(comm_size-1, xi)) /))) &
         CALL test_abort("error in xt_xmap_get_source_ranks", &
         filename, __LINE__)
    DEALLOCATE(ranks)
  END SUBROUTINE check_allgather_analog_xmap

  SUBROUTINE test_allgather_analog(xmap_new, num_indices_per_rank, comm)
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
    INTEGER, INTENT(inout) :: comm
    INTEGER, INTENT(in) :: num_indices_per_rank
    INTEGER(xi), ALLOCATABLE :: src_index_list(:)
    INTEGER(xi) :: i
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap, xmap_copy
    TYPE(xt_stripe) :: dst_index_stripe(1)
    INTEGER :: comm_size, comm_rank
    INTEGER(xi) :: comm_rank_xi, num_indices_per_rank_xi
    LOGICAL :: is_inter

    CALL get_rank_range(comm, is_inter, comm_rank, comm_size)
    comm_rank_xi = INT(comm_rank, xi)
    num_indices_per_rank_xi = INT(num_indices_per_rank, xi)
    ! setup
    ALLOCATE(src_index_list(num_indices_per_rank))
    DO i = 1_xi, num_indices_per_rank
      src_index_list(i) = comm_rank_xi * num_indices_per_rank_xi + i - 1_xi
    END DO
    src_idxlist = xt_idxvec_new(src_index_list)
    dst_index_stripe(1) = xt_stripe(0, 1, comm_size * num_indices_per_rank)
    dst_idxlist = xt_idxstripes_new(dst_index_stripe)
    xmap = xmap_new(src_idxlist, dst_idxlist, comm)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    ! verify expected results
    CALL check_allgather_analog_xmap(xmap, comm)
    xmap_copy = xt_xmap_copy(xmap)
    CALL check_allgather_analog_xmap(xmap, comm)

    ! clean up
    CALL xt_xmap_delete(xmap)
    CALL xt_xmap_delete(xmap_copy)
  END SUBROUTINE test_allgather_analog

  SUBROUTINE check_ring_xmap(xmap, dst_index_list, is_inter)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(xt_int_kind), INTENT(in) :: dst_index_list(2)
    LOGICAL, INTENT(in) :: is_inter
    INTEGER :: ranks(2), num_dst, num_src
    num_dst = xt_xmap_get_num_destinations(xmap)
    IF (.NOT. is_inter .AND. (num_dst > 2 .OR. num_dst < 1)) &
         CALL test_abort("error in xt_xmap_get_num_destinations", &
         filename, __LINE__)

    num_src = xt_xmap_get_num_sources(xmap)
    IF (num_src > 2 .OR. num_src < 1) &
         CALL test_abort("error in xt_xmap_get_num_sources", filename, __LINE__)

    IF (.NOT. is_inter) THEN
      CALL xt_xmap_get_destination_ranks(xmap, ranks)
      CALL xt_sort_int(ranks(1:num_dst))
      IF (ANY(ranks /= dst_index_list)) &
           CALL test_abort("error in xt_xmap_get_destination_ranks", &
           filename, __LINE__)
    END IF

    CALL xt_xmap_get_source_ranks(xmap, ranks)
    CALL xt_sort_int(ranks(1:num_src))
    IF (ANY(ranks /= dst_index_list)) &
         CALL test_abort("error in xt_xmap_get_source_ranks", &
         filename, __LINE__)
  END SUBROUTINE check_ring_xmap

  SUBROUTINE test_ring_1d(xmap_new, comm)
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
    INTEGER, INTENT(inout) :: comm
    ! test in which each process talks WITH two other processes
    INTEGER(xt_int_kind) :: src_index_list(1), dst_index_list(2), temp
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap, xmap_copy
    INTEGER :: comm_size, comm_rank
    LOGICAL :: is_inter

    CALL get_rank_range(comm, is_inter, comm_rank, comm_size)
    src_index_list(1) = INT(comm_rank, xi)
    src_idxlist = xt_idxvec_new(src_index_list)

    ! destination index list
    dst_index_list(1) = INT(MOD(comm_rank + comm_size - 1, comm_size), xi)
    dst_index_list(2) = INT(MOD(comm_rank             + 1, comm_size), xi)
    IF (dst_index_list(1) > dst_index_list(2)) THEN
      temp = dst_index_list(1)
      dst_index_list(1) = dst_index_list(2)
      dst_index_list(2) = temp
    END IF
    dst_idxlist = xt_idxvec_new(dst_index_list, 2)

    ! test of exchange map
    xmap = xmap_new(src_idxlist, dst_idxlist, comm)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    ! test results
    CALL check_ring_xmap(xmap, dst_index_list, is_inter)
    xmap_copy = xt_xmap_copy(xmap)
    CALL check_ring_xmap(xmap_copy, dst_index_list, is_inter)

    ! clean up
    CALL xt_xmap_delete(xmap)
    CALL xt_xmap_delete(xmap_copy)

  END SUBROUTINE test_ring_1d

  SUBROUTINE test_maxpos(xmap_new, comm, indices_per_rank)
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
    INTEGER, INTENT(in) :: indices_per_rank
    ! first setup simple pattern of boundary exchange
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap, xmup, xmup2, xmsp
    INTEGER :: indices_for_exch
    INTEGER(xt_int_kind) :: src_index(indices_per_rank), &
         dst_index(indices_per_rank)
    INTEGER :: max_pos_src, max_pos_dst, max_pos_src_u, max_pos_dst_u, &
         max_pos_src_u2, max_pos_dst_u2, max_pos_src_s, max_pos_dst_s
    INTEGER :: comm_rank, comm_size, world_size
    INTEGER :: ierror
    INTEGER :: i, xmspread(2)
    INTEGER :: pos_update1(indices_per_rank), pos_update2(2*indices_per_rank)

    CALL mpi_comm_rank(comm, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_rank", filename, __LINE__)
    CALL mpi_comm_size(comm, comm_size, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_size", filename, __LINE__)

    world_size = comm_size * indices_per_rank

    ! setup
    indices_for_exch = indices_per_rank/2
    DO i = 1, indices_per_rank
      src_index(i) = INT(i-1 + comm_rank * indices_per_rank, xi)
    END DO
    DO i = 1, indices_for_exch
      dst_index(i) = INT(MOD(i - 1 - indices_for_exch &
           + (comm_rank+comm_size) * indices_per_rank, world_size), xi)
    END DO
    DO i = indices_for_exch+1, indices_per_rank-indices_for_exch
      dst_index(i) = INT(i-1 + comm_rank * indices_per_rank, xi)
    END DO
    DO i = 1, indices_for_exch
      dst_index(indices_per_rank-indices_for_exch+i) &
           = INT(MOD(i + (comm_rank+1) * indices_per_rank, world_size), xi)
    END DO
    src_idxlist = xt_idxvec_new(src_index)
    dst_idxlist = xt_idxvec_new(dst_index)

    xmap = xmap_new(src_idxlist, dst_idxlist, comm)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    ! test
    ! 1. test that initial max positions are in range
    max_pos_dst = xt_xmap_get_max_dst_pos(xmap)
    max_pos_src = xt_xmap_get_max_src_pos(xmap)
    IF (max_pos_src < indices_per_rank-1) &
         CALL test_abort("error in xt_xmap_get_max_src_pos", filename, __LINE__)
    IF (max_pos_dst < indices_per_rank-1) &
         CALL test_abort("error in xt_xmap_get_max_dst_pos", filename, __LINE__)

    ! 2. expand range and verify it is reflected in max pos
    DO i = 1,indices_per_rank
      pos_update1(i) = (i-1)*2
    END DO

    xmup = xt_xmap_update_positions(xmap, pos_update1, pos_update1)

    max_pos_dst_u = xt_xmap_get_max_dst_pos(xmup)
    max_pos_src_u = xt_xmap_get_max_src_pos(xmup)
    IF (max_pos_src_u < (indices_per_rank-1)*2) &
         CALL test_abort("error in xt_xmap_get_max_src_pos", filename, __LINE__)
    IF (max_pos_dst_u < (indices_per_rank-1)*2) &
         CALL test_abort("error in xt_xmap_get_max_dst_pos", filename, __LINE__)

    ! 3. contract range again and verify max pos is updated
    DO i = 1, indices_per_rank*2
      pos_update2(i) = (i-1)/2
    END DO
    xmup2 = xt_xmap_update_positions(xmap, pos_update2, pos_update2)

    max_pos_dst_u2 = xt_xmap_get_max_dst_pos(xmup2)
    max_pos_src_u2 = xt_xmap_get_max_src_pos(xmup2)
    IF (max_pos_src_u2 >= indices_per_rank) &
         CALL test_abort("error in xt_xmap_get_max_src_pos", filename, __LINE__)
    IF (max_pos_dst_u2 >= indices_per_rank) &
         CALL test_abort("error in xt_xmap_get_max_dst_pos", filename, __LINE__)

    ! 4. apply spread and check max pos range
    xmspread(1) = 0
    xmspread(2) = indices_per_rank*3
    xmsp = xt_xmap_spread(xmap, 2, xmspread, xmspread)
    max_pos_dst_s = xt_xmap_get_max_dst_pos(xmsp)
    max_pos_src_s = xt_xmap_get_max_src_pos(xmsp)
    IF (max_pos_dst_s < (indices_per_rank-1)*3) &
         CALL test_abort("error in xt_xmap_get_max_dst_pos", filename, __LINE__)
    IF (max_pos_src_s < (indices_per_rank-1)*3) &
         CALL test_abort("error in xt_xmap_get_max_src_pos", filename, __LINE__)

    ! cleanup
    CALL xt_xmap_delete(xmap)
    CALL xt_xmap_delete(xmup)
    CALL xt_xmap_delete(xmup2)
    CALL xt_xmap_delete(xmsp)
  END SUBROUTINE test_maxpos

  SUBROUTINE check_pair_xmap(xmap)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER :: ranks(2)
    ! test results
    IF (xt_xmap_get_num_destinations(xmap) /= 2) &
         CALL test_abort("error in xt_xmap_get_num_destinations", &
         filename, __LINE__)

    IF (xt_xmap_get_num_sources(xmap) /= 2) &
         CALL test_abort("error in xt_xmap_get_num_sources", filename, __LINE__)

    CALL xt_xmap_get_destination_ranks(xmap, ranks)
    IF (ranks(1) /= 0 .OR. ranks(2) /= 1) &
         CALL test_abort("error in xt_xmap_get_destination_ranks", &
         filename, __LINE__)

    CALL xt_xmap_get_source_ranks(xmap, ranks)
    IF (ranks(1) /= 0 .OR. ranks(2) /= 1) &
         CALL test_abort("error in xt_xmap_get_source_ranks", &
         filename, __LINE__)
  END SUBROUTINE check_pair_xmap

  SUBROUTINE test_pair(xmap_new, comm)
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
    !src_index_list(index, rank)
    INTEGER(xt_int_kind) :: i, j, k
#ifdef __xlC__
    INTEGER(xt_int_kind), PARAMETER :: src_index_list(20, 0:1) = RESHAPE((/ &
         &  1_xi,  2_xi,  3_xi,  4_xi,  5_xi, &
         &  9_xi, 10_xi, 11_xi, 12_xi, 13_xi, &
         & 17_xi, 18_xi, 19_xi, 20_xi, 21_xi, &
         & 25_xi, 26_xi, 27_xi, 28_xi, 29_xi, &
         &  4_xi,  5_xi,  6_xi,  7_xi,  8_xi, &
         & 12_xi, 13_xi, 14_xi, 15_xi, 16_xi, &
         & 20_xi, 21_xi, 22_xi, 23_xi, 24_xi, &
         & 28_xi, 29_xi, 30_xi, 31_xi, 32_xi /),  &
         (/ 20, 2 /))
#else
    INTEGER(xt_int_kind), PARAMETER :: src_index_list(20, 0:1) = RESHAPE((/ &
         (((i + j * 8_xi + k * 3_xi, i = 1_xi, 5_xi), j = 0_xi,3_xi), &
         k = 0_xi,1_xi) /), (/ 20, 2 /))
#endif
    ! dst_index_list(index,rank)
    INTEGER(xt_int_kind), PARAMETER :: dst_index_list(20, 0:1) = RESHAPE((/ &
         10_xi, 15_xi, 14_xi, 13_xi, 12_xi, &
         15_xi, 10_xi, 11_xi, 12_xi, 13_xi, &
         23_xi, 18_xi, 19_xi, 20_xi, 21_xi, &
         31_xi, 26_xi, 27_xi, 28_xi, 29_xi, &
         13_xi, 12_xi, 11_xi, 10_xi, 15_xi, &
         12_xi, 13_xi, 14_xi, 15_xi, 10_xi, &
         20_xi, 21_xi, 22_xi, 23_xi, 18_xi, &
         28_xi, 29_xi, 30_xi, 31_xi, 26_xi /),  &
         (/ 20, 2 /))
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap, xmap_copy
    INTEGER :: comm_rank, ierror

    CALL mpi_comm_rank(comm, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_rank", filename, __LINE__)

    src_idxlist = xt_idxvec_new(src_index_list(:, comm_rank))

    ! destination index list
    dst_idxlist = xt_idxvec_new(dst_index_list(:, comm_rank))

    ! test of exchange map
    xmap = xmap_new(src_idxlist, dst_idxlist, comm)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    CALL check_pair_xmap(xmap)
    xmap_copy = xt_xmap_copy(xmap)
    CALL check_pair_xmap(xmap_copy)

    ! clean up
    CALL xt_xmap_delete(xmap)
    CALL xt_xmap_delete(xmap_copy)
  END SUBROUTINE test_pair

  SUBROUTINE check_ping_pong_xmap(xmap, comm, ping_rank, pong_rank)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: comm, ping_rank, pong_rank
    INTEGER :: expect, dst_rank(1), src_rank(1), comm_rank, ierror
    CHARACTER(len=80) :: msg

    CALL mpi_comm_rank(comm, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort('error calling mpi_comm_rank', filename, __LINE__)
    WRITE (msg, '(a,i0,a)') "error in xt_xmap_get_num_destinations (rank == ", &
         comm_rank, ")"
    expect = MERGE(1, 0, comm_rank == ping_rank)
    IF (xt_xmap_get_num_destinations(xmap) /= expect) &
         CALL test_abort(msg, filename, __LINE__)

    expect = MERGE(1, 0, comm_rank == pong_rank)
    IF (xt_xmap_get_num_sources(xmap) /= expect) &
         CALL test_abort(msg, filename, __LINE__)

    IF (comm_rank == ping_rank) THEN
      CALL xt_xmap_get_destination_ranks(xmap, dst_rank)
      IF (dst_rank(1) /= pong_rank) &
           CALL test_abort("error in xt_xmap_get_destination_ranks", &
           filename, __LINE__)
    END IF
    IF (comm_rank == pong_rank) THEN
      CALL xt_xmap_get_source_ranks(xmap, src_rank)
      IF (src_rank(1) /= ping_rank) &
           CALL test_abort("error in xt_xmap_get_source_ranks", &
           filename, __LINE__)
    END IF
  END SUBROUTINE check_ping_pong_xmap

  SUBROUTINE test_ping_pong(xmap_new, comm, ping_rank, pong_rank)
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
    INTEGER, INTENT(in) :: ping_rank, pong_rank
    INTEGER, INTENT(inout) :: comm
    INTEGER(xt_int_kind), PARAMETER :: &
         index_list(5) = (/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi /)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap, xmap_copy
    INTEGER :: comm_rank, comm_size
    LOGICAL :: is_inter
    CALL get_rank_range(comm, is_inter, comm_rank, comm_size)
    IF (comm_rank == ping_rank) THEN
      src_idxlist = xt_idxvec_new(index_list)
    ELSE
      src_idxlist = xt_idxempty_new()
    END IF


    IF (comm_rank == pong_rank) THEN
      dst_idxlist = xt_idxvec_new(index_list)
    ELSE
      dst_idxlist = xt_idxempty_new()
    END IF

    ! test of exchange map

    xmap = xmap_new(src_idxlist, dst_idxlist, comm)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    ! test results
    CALL check_ping_pong_xmap(xmap, comm, ping_rank, pong_rank)
    xmap_copy = xt_xmap_copy(xmap)
    CALL check_ping_pong_xmap(xmap_copy, comm, ping_rank, pong_rank)
    ! clean up
    CALL xt_xmap_delete(xmap)
    CALL xt_xmap_delete(xmap_copy)
  END SUBROUTINE test_ping_pong
END MODULE test_xmap_common_parallel
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
