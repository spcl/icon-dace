!>
!! @file test_redist_p2p_f.f90
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
!  Keywords:
!  Maintainer: Jörg Behrens <behrens@dkrz.de>
!              Moritz Hanke <hanke@dkrz.de>
!              Thomas Jahns <jahns@dkrz.de>
!  URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are  permitted provided that the following conditions are
!  met:
!
!  Redistributions of source code must retain the above copyright notice,
!  this list of conditions and the following disclaimer.
!
!  Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.
!
!  Neither the name of the DKRZ GmbH nor the names of its contributors
!  may be used to endorse or promote products derived from this software
!  without specific prior written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
#include "fc_feature_defs.inc"
PROGRAM test_redist_p2p_f
  USE mpi
  USE yaxt, ONLY: xt_int_kind, xt_xmap, xt_idxlist, xt_redist, xt_offset_ext, &
       xi => xt_int_kind, xt_int_mpidt, xt_initialize, xt_finalize, &
       xt_idxvec_new, xt_idxlist_delete, &
       xt_redist_p2p_new, xt_redist_p2p_off_custom_new, xt_redist_p2p_ext_new, &
       xt_redist_s_exchange, &
       xt_redist_copy, xt_redist_delete,  xt_redist_get_mpi_comm, &
       xt_xmap_all2all_new, xt_xmap_delete, &
       xt_config, xt_config_delete
  ! pgfortran is in most versions well incapable of handling multiply extended
  ! generic interfaces
#ifdef __PGI
  USE xt_redist_logical, ONLY: xt_redist_s_exchange
#endif
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_redist_common, ONLY: check_redist, communicators_are_congruent, &
       redist_exchanger_option
  USE test_idxlist_utils, ONLY: test_err_count
  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER :: filename = 'test_redist_p2p_f.f90'
  TYPE(xt_config) :: config

  CALL init_mpi

  CALL xt_initialize(mpi_comm_world)
  config = redist_exchanger_option()

  ! offset-free test:
  ! source index list
  CALL test_without_offsets
  CALL test_with_offsets
  CALL test_offset_extents

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)

  CALL xt_config_delete(config)
  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  SUBROUTINE test_without_offsets
    INTEGER, PARAMETER :: src_num_indices = 14, dst_num_indices = 13
    INTEGER(xt_int_kind), PARAMETER :: src_index_list(src_num_indices) &
         = (/  5_xi, 67_xi,  4_xi,  5_xi, 13_xi, &
         &     9_xi,  2_xi,  1_xi,  0_xi, 96_xi, &
         &    13_xi, 12_xi,  1_xi,  3_xi /), &
         dst_index_list(dst_num_indices) = &
         & (/  5_xi,  4_xi,  3_xi, 96_xi,  1_xi, &
         &     5_xi,  4_xi,  5_xi,  4_xi,  3_xi, &
         &    13_xi,  2_xi,  1_xi /)
    INTEGER :: i
#ifndef __PGI
    DOUBLE PRECISION, PARAMETER :: src_data(src_num_indices) = &
         (/ (DBLE(i), i=0,src_num_indices-1) /)
#else
    ! for PGI Fortran DBLE must be evaluated at run-time
    DOUBLE PRECISION :: src_data(src_num_indices)
#endif
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(dst_num_indices) &
         = (/ 0.0d0,  2.0d0, 13.0d0,  9.0d0,  7.0d0, &
         &    0.0d0,  2.0d0,  0.0d0,  2.0d0, 13.0d0, &
         &    4.0d0,  6.0d0,  7.0d0 /)
    LOGICAL :: src_l(src_num_indices), &
         dst_l(dst_num_indices), ref_dst_l(dst_num_indices)
    DOUBLE PRECISION :: dst_data(dst_num_indices)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist_dp, redist_copy, redist_l

#ifdef __PGI
    DO i = 1, src_num_indices
      src_data(i) = DBLE(i - 1)
    END DO
#endif

    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices)

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices)

    ! xmap
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)

    ! redist_p2p
    redist_dp = xt_redist_p2p_new(xmap, mpi_double_precision, config)

    ! test communicator of redist
    IF (.NOT. communicators_are_congruent(xt_redist_get_mpi_comm(redist_dp), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_Comm", filename, __LINE__)

    ! test exchange
    CALL check_redist(redist_dp, src_data, dst_data, ref_dst_data)

    ! repeat for logicals
    src_l = NINT(MOD(src_data, 2.0d0)) == 1
    dst_l = .FALSE.
    ref_dst_l = NINT(MOD(ref_dst_data, 2.0d0)) == 1
    redist_l = xt_redist_p2p_new(xmap, mpi_logical, config)
    IF (.NOT. communicators_are_congruent(xt_redist_get_mpi_comm(redist_l), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_Comm", filename, __LINE__)
    CALL xt_redist_s_exchange(redist_l, src_l, dst_l)
    IF (ANY(dst_l .NEQV. ref_dst_l)) &
         CALL test_abort("error in xt_redist_s_exchange for 1D logical array", &
         filename, __LINE__)
    redist_copy = xt_redist_copy(redist_dp)
    CALL xt_redist_delete(redist_dp)
    CALL check_redist(redist_copy, src_data, dst_data, ref_dst_data)

    ! clean up
    CALL xt_redist_delete(redist_copy)
    CALL xt_redist_delete(redist_l)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END SUBROUTINE test_without_offsets

  SUBROUTINE test_with_offsets
    ! source index list
    INTEGER, PARAMETER :: src_num = 14, dst_num = 13
    INTEGER(xt_int_kind), PARAMETER :: src_index_list(src_num) = &
         (/  5_xi, 67_xi,  4_xi,  5_xi, 13_xi, &
         &   9_xi,  2_xi,  1_xi,  0_xi, 96_xi, &
         &  13_xi, 12_xi,  1_xi,  3_xi /), &
         dst_index_list(dst_num) = &
         (/  5_xi,  4_xi,  3_xi, 96_xi,  1_xi, &
         &   5_xi,  4_xi,  5_xi,  4_xi,  3_xi, &
         &  13_xi,  2_xi,  1_xi /)
    INTEGER :: i
    INTEGER, PARAMETER :: src_pos(src_num) = (/ (i, i = 0, src_num - 1) /), &
         dst_pos(dst_num) = (/ ( dst_num - i, i = 1, dst_num ) /)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_copy
#ifndef __PGI
    DOUBLE PRECISION, PARAMETER :: src_data(src_num) = &
         (/ (DBLE(i), i=0,src_num-1) /)
#else
    ! for PGI Fortran DBLE must be evaluated at run-time
    DOUBLE PRECISION :: src_data(src_num)
#endif
    DOUBLE PRECISION :: dst_data(dst_num)
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(dst_num) = &
         (/  0.0d0,  2.0d0, 13.0d0,  9.0d0,  7.0d0, &
         &   0.0d0,  2.0d0,  0.0d0,  2.0d0, 13.0d0, &
         &   4.0d0,  6.0d0,  7.0d0 /)

#ifdef __PGI
    DO i = 1, src_num
      src_data(i) = DBLE(i - 1)
    END DO
#endif

    src_idxlist = xt_idxvec_new(src_index_list)

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num)

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)

    ! redist_p2p with offsets
    redist = xt_redist_p2p_off_custom_new(xmap, src_pos, dst_pos, &
         mpi_double_precision, config)

    ! test communicator of redist
    IF (.NOT. communicators_are_congruent(xt_redist_get_mpi_comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_MPI_Comm", filename, __LINE__)

    ! test exchange
    CALL check_redist(redist, src_data, dst_data, ref_dst_data(dst_num:1:-1))

    redist_copy = xt_redist_copy(redist)
    CALL xt_redist_delete(redist)
    CALL check_redist(redist_copy, src_data, dst_data, &
         ref_dst_data(dst_num:1:-1))

    ! clean up
    CALL xt_redist_delete(redist_copy)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END SUBROUTINE test_with_offsets

  SUBROUTINE test_offset_extents
    ! source/destination index lists
    INTEGER, PARAMETER :: src_num = 14, dst_num = 13
    INTEGER(xt_int_kind), PARAMETER :: src_index_list(src_num) = &
         (/  5_xi, 67_xi,  4_xi,  5_xi, 13_xi, &
         &   9_xi,  2_xi,  1_xi,  0_xi, 96_xi, &
         &  13_xi, 12_xi,  1_xi,  3_xi /), &
         dst_index_list(dst_num) = &
         (/  5_xi,  4_xi,  3_xi, 96_xi,  1_xi, &
         &   5_xi,  4_xi,  5_xi,  4_xi,  3_xi, &
         &  13_xi, 2_xi, 1_xi /)
#ifdef __G95__
    INTEGER :: i
    INTEGER(xt_int_kind), PARAMETER :: src_data(src_num) &
         = (/ (INT(i, xi), i = 0, 13) /)
#else
    INTEGER(xt_int_kind) :: i
    INTEGER(xt_int_kind), PARAMETER :: src_data(src_num) &
         = (/ (i, i = 0_xi, 13_xi) /)
#endif
    INTEGER(xt_int_kind) :: dst_data(dst_num)
    INTEGER(xt_int_kind), PARAMETER :: ref_dst_data(dst_num) = &
         (/  7_xi,  6_xi,  4_xi, 13_xi,  2_xi, &
         &   0_xi,  2_xi,  0_xi,  7_xi,  9_xi, &
         &  13_xi,  2_xi,  0_xi /)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_copy
    TYPE(xt_offset_ext), PARAMETER :: &
         src_pos(1) = (/ xt_offset_ext(0, src_num, 1) /), &
         dst_pos(1) = (/ xt_offset_ext(dst_num - 1, dst_num, -1) /)

    src_idxlist = xt_idxvec_new(src_index_list)
    dst_idxlist = xt_idxvec_new(dst_index_list)

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)

    ! redist_p2p with extents of offsets
    redist = xt_redist_p2p_ext_new(xmap, src_pos, dst_pos, xt_int_mpidt, config)
    ! test communicator of redist
    IF (.NOT. communicators_are_congruent(xt_redist_get_MPI_Comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_MPI_Comm", &
         filename, __LINE__)

    ! test exchange
    CALL check_redist(redist, src_data, dst_data, ref_dst_data)

    redist_copy = xt_redist_copy(redist)
    CALL xt_redist_delete(redist)
    CALL check_redist(redist_copy, src_data, dst_data, ref_dst_data)

    ! clean up
    CALL xt_redist_delete(redist_copy)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END SUBROUTINE test_offset_extents

END PROGRAM test_redist_p2p_f
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
