!>
!! @file unstruct_halo_exchange_f.f90
!! @brief Fortran example of unstructured halo exchange
!!
!! @copyright Copyright  (C)  2015 Moritz Hanke <hanke@dkrz.de>
!!
!! @author Moritz Hanke <hanke@dkrz.de>
!!

!
! Keywords:
! Maintainer: Moritz Hanke <hanke@dkrz.de>
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

  ! This program demonstrates the halo exchange for unstructured data that is
  ! distributed among three processes. The "global" grid looks as follows
  ! (The numbers represent the indices of the corners of the grid.)
  !
  !      01------04------08------13
  !     /  \    /  \    /  \    /  \
  !    /    \  /    \  /    \  /    \
  !  00------03------07------12------17
  !    \    /  \    /  \    /  \    /
  !     \  /    \  /    \  /    \  /
  !      02------06------11------16
  !        \    /  \    /  \    /
  !         \  /    \  /    \  /
  !          05------10------15
  !            \    /  \    /
  !             \  /    \  /
  !              09------14
  !
  ! Each process has a part of this grid plus a halo around it (marked with '*').
  !
  ! Process 0 (halo points: 8, 12, 11 10, 5):
  !
  !      01------04******08
  !     /  \    /  \    *  *
  !    /    \  /    \  *    *
  !  00------03------07******12
  !    \    /  \    /  *    *
  !     \  /    \  /    *  *
  !      02------06******11
  !        *    *  *    *
  !         *  *    *  *
  !          05******10
  !
  ! Process 1 (halo points: 4, 3, 6, 10, 15):
  !
  !      04******08------13
  !     *  *    /  \    /  \
  !    *    *  /    \  /    \
  !  03******07------12------17
  !    *    *  \    /  \    /
  !     *  *    \  /    \  /
  !      06******11------16
  !        *    *  *    *
  !         *  *    *  *
  !          10******15
  !
  ! Process 2 (halo points: 3, 7, 12, 2, 16)
  !
  !      03******07******12
  !     *  *    *  *    *  *
  !    *    *  *    *  *    *
  !  02******06------11******16
  !    *    /  \    /  \    *
  !     *  /    \  /    \  *
  !      05------10------15
  !        \    /  \    /
  !         \  /    \  /
  !          09------14
#include "fc_feature_defs.inc"
PROGRAM unstruct_halo_exchang
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, &
       xt_idxlist, xt_idxvec_new, xt_idxlist_delete, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_off_new, xt_redist_delete, &
       xt_redist_s_exchange, xi => xt_int_kind
  ! PGI compilers do not handle generic interface correctly
#if defined __PGI && ( __PGIC__ == 15 || __PGIC__ == 14 )
  USE xt_redist_real_sp, ONLY: xt_redist_s_exchange
  USE xt_redist_real_dp, ONLY: xt_redist_s_exchange
#endif
  IMPLICIT NONE

  INTEGER :: comm_rank, comm_size
  INTEGER :: i, j, rank, ierror

  INTEGER(xi) :: local_indices(12, 0:2), idx
  INTEGER(xi) :: halo_indices(5, 0:2)
  INTEGER(xi) :: src_indices(12, 0:2), tgt_indices(5, 0:2)
  INTEGER :: src_offsets(12)
  INTEGER :: tgt_offsets(5)

  TYPE(xt_idxlist) :: src_idxlist, tgt_idxlist
  TYPE(xt_xmap) :: xmap
  TYPE(xt_redist) :: redist
  DOUBLE PRECISION :: src_data(12), tgt_data(12)

  CHARACTER(64) :: fmt

  CALL mpi_init(ierror)
  IF (ierror /= mpi_success) STOP 1
  CALL xt_initialize(mpi_comm_world)
  CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
  IF (ierror /= mpi_success) STOP 1
  CALL mpi_comm_size (mpi_comm_world, comm_size, ierror)
  IF (ierror /= mpi_success) STOP 1

  IF (comm_size /= 3) STOP 3

  ! create description of locally stored indices...
  local_indices(:,0) = (/    1_xi,  4_xi,  8_xi, &
       &                  0_xi,  3_xi,  7_xi, 12_xi, &
       &                     2_xi,  6_xi, 11_xi, &
       &                         5_xi, 10_xi                 /)
  local_indices(:,1) = (/    4_xi,  8_xi, 13_xi, &
       &                  3_xi,  7_xi, 12_xi, 17_xi, &
       &                     6_xi, 11_xi, 16_xi, &
       &                        10_xi, 15_xi                 /)
  local_indices(:,2) = (/    3_xi,  7_xi, 12_xi, &
       &                  2_xi,  6_xi, 11_xi, 16_xi, &
       &                     5_xi, 10_xi, 15_xi, &
       &                         9_xi, 14_xi                 /)

  ! ...and list of which of those are shadows of data "owned" by another rank
  halo_indices(:,0) = (/  8_xi, 12_xi, 11_xi, 10_xi,  5_xi /)
  halo_indices(:,1) = (/  4_xi,  3_xi,  6_xi, 10_xi, 15_xi /)
  halo_indices(:,2) = (/  3_xi,  7_xi, 12_xi,  2_xi, 16_xi /)

  DO rank = 0, comm_size - 1
    DO i = 1, SIZE(local_indices, 1)
      idx = local_indices(i, rank)
      src_indices(i, rank) = MERGE(-1_xi, idx, ANY(halo_indices(:, rank) == idx))
    END DO
    tgt_indices(:, rank) = halo_indices(:, rank)
  END DO

  ! create decomposition descriptors
  src_idxlist = xt_idxvec_new(src_indices(:,comm_rank), SIZE(src_indices, 1))
  tgt_idxlist = xt_idxvec_new(tgt_indices(:,comm_rank), SIZE(tgt_indices, 1))

  ! generate exchange map
  xmap = xt_xmap_all2all_new(src_idxlist, tgt_idxlist, mpi_comm_world)

  ! generate redistribution object
  src_offsets = (/(i, i = 0, SIZE(src_indices, 1) - 1)/)
  DO i = 1, SIZE(halo_indices, 1)
    DO j = 1, SIZE(local_indices, 1)
      IF (halo_indices(i, comm_rank) == local_indices(j, comm_rank)) &
        tgt_offsets(i) = j - 1
    END DO
  END DO
  redist = xt_redist_p2p_off_new(xmap, src_offsets, tgt_offsets, &
    &                            mpi_double_precision)

  ! prepare arrays
  src_data(:) = DBLE(src_indices(:,comm_rank))
  tgt_data(:) = DBLE(src_indices(:,comm_rank))

  fmt = '(I1, " ", A, " exchange: ", xx(I3))'
  WRITE(fmt(29:30), '(I2)') SIZE(src_data, 1)
  PRINT fmt, comm_rank, 'before', INT(src_data(:))

  ! do the exchange
  CALL xt_redist_s_exchange(redist, src_data, tgt_data)

  PRINT fmt, comm_rank, 'after ', INT(tgt_data(:))

  ! clean up
  CALL xt_redist_delete(redist)
  CALL xt_xmap_delete(xmap)
  CALL xt_idxlist_delete(tgt_idxlist)
  CALL xt_idxlist_delete(src_idxlist)

  ! finalise
  CALL xt_finalize()
  CALL mpi_finalize(ierror)
  IF (ierror /= mpi_success) STOP 1
END PROGRAM unstruct_halo_exchang
!
! Local Variables:
! coding: utf-8
! f90-continuation-indent: 5
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
! license-default: "bsd"
! End:
!
