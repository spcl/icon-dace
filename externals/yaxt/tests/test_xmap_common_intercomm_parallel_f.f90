!>
!! @file test_xmap_common_intercomm_parallel_f.f90
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
MODULE test_xmap_common_intercomm_parallel
  USE iso_c_binding, ONLY: c_int
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, posix_exit
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, xt_stripe, &
       xi => xt_int_kind, &
       xt_idxlist, xt_idxlist_delete, xt_idxstripes_new, &
       xt_xmap, xt_xmap_copy, xt_xmap_delete, &
       xt_xmap_get_num_destinations, xt_xmap_get_num_sources, &
       xt_xmap_get_destination_ranks, xt_xmap_get_source_ranks, &
       xt_sort_int
  USE xt_core, ONLY: i8
  USE test_xmap_common_parallel, ONLY: get_rank_range, &
       check_allgather_analog_xmap, test_ping_pong, test_ring_1d
  IMPLICIT NONE
  PRIVATE
  INTEGER :: intra_group_comm
  PUBLIC :: xmap_intercomm_parallel_test_main, intra_group_comm
  CHARACTER(len=*), PARAMETER :: &
       filename = 'test_xmap_common_intercomm_parallel_f.f90'
CONTAINS
  SUBROUTINE xmap_intercomm_parallel_test_main(xmap_new, call_initialize, &
       call_finalize)
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
    LOGICAL, OPTIONAL, INTENT(in) :: call_initialize, call_finalize
    INTEGER :: comm, comm_rank, comm_size, ierror, inter_comm, split_rank, &
         retval
    LOGICAL :: in_second_group, call_finalize_, call_initialize_

    IF (PRESENT(call_initialize)) THEN
      call_initialize_ = call_initialize
    ELSE
      call_initialize_ = .TRUE.
    END IF
    IF (PRESENT(call_finalize)) THEN
      call_finalize_ = call_finalize
    ELSE
      call_finalize_ = .TRUE.
    END IF
    IF (call_initialize_) THEN
      CALL init_mpi
      CALL xt_initialize(mpi_comm_world)
    END IF
    comm = mpi_comm_world
    CALL mpi_comm_rank(comm, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_rank", filename, __LINE__)
    CALL mpi_comm_size(comm, comm_size, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_size", filename, __LINE__)
    IF (comm_size > HUGE(1_xi)) &
         CALL test_abort("number of ranks exceeds test limit", &
         filename, __LINE__)

    IF (comm_size > 1) THEN
      retval = 0
      split_rank = comm_size/2 - MERGE(1, 0, comm_size > 5)
      in_second_group = comm_rank >= split_rank
      CALL mpi_comm_split(comm, MERGE(1, 0, in_second_group), 0, &
           intra_group_comm, ierror)
      IF (ierror /= mpi_success) &
           CALL test_abort("error calling mpi_comm_split", filename, __LINE__)
      CALL mpi_intercomm_create(intra_group_comm, 0, comm, &
           MERGE(0, split_rank, in_second_group), 0, inter_comm, ierror)
      IF (ierror /= mpi_success) &
           CALL test_abort("error calling mpi_intercomm_create", &
           filename, __LINE__)
      CALL test_allgather_analog(xmap_new, 1_xi, inter_comm)
      ! repeat test for large index list that will cause stripifying
      CALL test_allgather_analog(xmap_new, 1024_xi, inter_comm)

      CALL test_peer(xmap_new, inter_comm)

      CALL test_ping_pong(xmap_new, inter_comm, 0, 0)

      IF (split_rank > 2) CALL test_ring_1d(xmap_new, inter_comm)

      CALL mpi_comm_free(inter_comm, ierror)
      IF (ierror /= mpi_success) &
           CALL test_abort("error calling mpi_comm_free", filename, __LINE__)
      CALL mpi_comm_free(intra_group_comm, ierror)
      IF (ierror /= mpi_success) &
           CALL test_abort("error calling mpi_comm_free", filename, __LINE__)
    ELSE
      retval = 77
    END IF
    IF (test_err_count() /= 0) &
         CALL test_abort("non-zero error count!", filename, __LINE__)
    IF (call_finalize_) THEN
      CALL xt_finalize
      CALL finish_mpi
    END IF
    IF (retval /= 0) CALL posix_exit(retval)
  END SUBROUTINE xmap_intercomm_parallel_test_main

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
    INTEGER(xt_int_kind), INTENT(in) :: num_indices_per_rank
    INTEGER, INTENT(inout) :: comm
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap, xmap_copy
    TYPE(xt_stripe) :: src_index_stripe(1), dst_index_stripe(1)
    INTEGER :: comm_rank, remote_size
    LOGICAL :: is_inter
    CALL get_rank_range(comm, is_inter, comm_rank, remote_size)
    ! setup
    src_index_stripe(1) = xt_stripe(INT(comm_rank, xi) * num_indices_per_rank, &
         1_xi, INT(num_indices_per_rank, c_int))
    src_idxlist = xt_idxstripes_new(src_index_stripe)
    dst_index_stripe(1) = xt_stripe(0, 1, &
         INT(INT(remote_size, xi) * num_indices_per_rank, c_int))
    dst_idxlist = xt_idxstripes_new(dst_index_stripe)
    xmap = xmap_new(src_idxlist, dst_idxlist, comm)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    ! verify expected results
    CALL check_allgather_analog_xmap(xmap, comm)
    xmap_copy = xt_xmap_copy(xmap)
    CALL check_allgather_analog_xmap(xmap_copy, comm)

    ! clean up
    CALL xt_xmap_delete(xmap)
    CALL xt_xmap_delete(xmap_copy)
  END SUBROUTINE test_allgather_analog


  SUBROUTINE check_peers(xmap, num_ref_ranks, ref_ranks, peer_rank_buf, &
       get_num_peers, get_peer_ranks, get_num_peers_name, get_peer_ranks_name)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: num_ref_ranks
    INTEGER(c_int), INTENT(in) :: ref_ranks(num_ref_ranks)
    INTEGER(c_int), INTENT(out) :: peer_rank_buf(num_ref_ranks)
    INTERFACE
      FUNCTION get_num_peers(xmap) RESULT(num)
        IMPORT :: xt_xmap
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER :: num
      END FUNCTION get_num_peers
      SUBROUTINE get_peer_ranks(xmap, ranks)
        IMPORT :: xt_xmap, c_int
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), INTENT(out) :: ranks(*)
      END SUBROUTINE get_peer_ranks
    END INTERFACE
    CHARACTER(len=*), INTENT(in) :: get_num_peers_name, get_peer_ranks_name
    CHARACTER(len=80) :: msg
    IF (get_num_peers(xmap) /= num_ref_ranks) THEN
      WRITE (msg, '(2a)') "error in ", get_num_peers_name
      CALL test_abort(msg, filename, __LINE__)
    END IF
    CALL get_peer_ranks(xmap, peer_rank_buf)
    CALL xt_sort_int(peer_rank_buf)
    IF (ANY(peer_rank_buf /= ref_ranks)) THEN
      WRITE (msg, '(2a)') "error in ", get_peer_ranks_name
      CALL test_abort(msg, filename, __LINE__)
    END IF
  END SUBROUTINE check_peers

  SUBROUTINE check_peer_xmap(xmap, stripe_in_local_group, remote_size, &
       global_num_idx)
    TYPE(xt_xmap), INTENT(in) :: xmap
    TYPE(xt_stripe), INTENT(in) :: stripe_in_local_group
    INTEGER, INTENT(in) :: remote_size
    INTEGER(i8), INTENT(in) :: global_num_idx
    INTEGER :: num_indices, idx_per_remote_rank, num_remote_peers, &
         last_seen_rank, i, remote_rank_i
    INTEGER(c_int), ALLOCATABLE :: ref_ranks(:), rank_buf(:)
    num_indices = INT(stripe_in_local_group%nstrides)
    idx_per_remote_rank = INT(global_num_idx / INT(remote_size, i8))
    ALLOCATE(ref_ranks(num_indices))
    num_remote_peers = 0
    last_seen_rank = -1
    DO i = 1, num_indices
      remote_rank_i = INT((stripe_in_local_group%start + INT(i-1, xi)) &
           &              / idx_per_remote_rank)
      IF (remote_rank_i /= last_seen_rank) THEN
        num_remote_peers = num_remote_peers + 1
        ref_ranks(num_remote_peers) = INT(remote_rank_i, c_int)
        last_seen_rank = remote_rank_i
      END IF
    END DO
    ALLOCATE(rank_buf(num_remote_peers))
    CALL check_peers(xmap, num_remote_peers, ref_ranks, rank_buf, &
         xt_xmap_get_num_destinations, xt_xmap_get_destination_ranks, &
         "xt_xmap_get_num_destinations", "xt_xmap_get_destination_ranks")
    CALL check_peers(xmap, num_remote_peers, ref_ranks, rank_buf, &
         xt_xmap_get_num_sources, xt_xmap_get_source_ranks, &
         "xt_xmap_get_num_sources", "xt_xmap_get_source_ranks")
  END SUBROUTINE check_peer_xmap

  FUNCTION gcd(a, b)
    INTEGER, INTENT(in) :: a, b
    INTEGER :: a_, b_, t_, gcd
    a_ = a ; b_ = b
    DO WHILE (b_ /= 0)
      t_ = b_
      b_ = MOD(a_, b_)
      a_ = t_
    END DO
    gcd = a_
  END FUNCTION gcd

  FUNCTION lcm(a, b)
    INTEGER(i8) :: lcm
    INTEGER, INTENT(in) :: a, b
    INTEGER :: t
    t = gcd(a, b)
    lcm = INT(a / t, i8) * INT(b, i8)
  END FUNCTION lcm

  SUBROUTINE test_peer(xmap_new, comm)
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
    INTEGER :: comm_rank, comm_size, remote_size, ierror
    INTEGER(i8) :: global_num_idx
    TYPE(xt_stripe) :: stripe_in_local_group(1)
    TYPE(xt_idxlist) :: idxlist
    TYPE(xt_xmap) :: xmap, xmap_copy
    CALL mpi_comm_rank(comm, comm_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_rank", filename, __LINE__)
    CALL mpi_comm_size(comm, comm_size, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_size", filename, __LINE__)
    CALL mpi_comm_remote_size(comm, remote_size, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("error calling mpi_comm_remote_size", &
         filename, __LINE__)
    global_num_idx = lcm(comm_size, remote_size)
    stripe_in_local_group(1) &
         = xt_stripe(INT(global_num_idx / INT(comm_size, i8) &
         &               * INT(comm_rank, i8), xi), &
         &           1, INT(global_num_idx / INT(comm_size, i8), c_int))
    idxlist = xt_idxstripes_new(stripe_in_local_group)
    xmap = xmap_new(idxlist, idxlist, comm)
    CALL xt_idxlist_delete(idxlist)
    CALL check_peer_xmap(xmap, stripe_in_local_group(1), remote_size, &
         global_num_idx)
    xmap_copy = xt_xmap_copy(xmap)
    CALL check_peer_xmap(xmap_copy, stripe_in_local_group(1), remote_size, &
         global_num_idx)
    CALL xt_xmap_delete(xmap)
    CALL xt_xmap_delete(xmap_copy)
  END SUBROUTINE test_peer

END MODULE test_xmap_common_intercomm_parallel
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
