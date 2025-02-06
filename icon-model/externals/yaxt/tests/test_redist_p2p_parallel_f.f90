!>
!! @file test_redist_p2p_parallel_f.f90
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
PROGRAM test_redist_p2p_parallel
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, &
       xt_int_kind, xi => xt_int_kind, &
       xt_idxlist, xt_idxvec_new, xt_idxlist_delete, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_get_mpi_comm, &
       xt_redist_delete, &
       xt_redist_p2p_blocks_off_custom_new, xt_redist_p2p_blocks_custom_new, &
       xt_config, xt_config_delete
  USE test_idxlist_utils, ONLY: test_err_count
  USE test_redist_common, ONLY: communicators_are_congruent, &
       check_redist, redist_exchanger_option
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER :: filename = 'test_redist_p2p_parallel_f.f90'
  TYPE(xt_config) :: config
  INTEGER :: comm_rank, comm_size, ierror

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  config = redist_exchanger_option()

  CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("MPI error!", filename, __LINE__)

  CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("MPI error!", filename, __LINE__)

  CALL simple_test
  CALL nonuniform_test
  CALL block_redist_test

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_config_delete(config)
  CALL xt_finalize
  CALL finish_mpi

CONTAINS
  SUBROUTINE simple_test
    INTEGER, PARAMETER :: data_size = 10
    INTEGER, PARAMETER :: src_num_indices = data_size, &
         dst_num_indices = data_size
    INTEGER(xt_int_kind) :: src_index_list(data_size), &
         dst_index_list(data_size)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist
    DOUBLE PRECISION :: src_data(data_size), dst_data(data_size)
    INTEGER :: i

    ! source index list
    DO i = 1, src_num_indices
      src_index_list(i) = INT(comm_rank * data_size + (i - 1), xi)
    END DO

    src_idxlist = xt_idxvec_new(src_index_list)
    ! destination index list
    DO i = 1, dst_num_indices
      dst_index_list(i) &
           = INT(MOD(comm_rank * data_size + i + 1, comm_size * data_size), xi)
    END DO
    dst_idxlist = xt_idxvec_new(dst_index_list)
    ! xmap
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)
    ! redist_p2p
    redist = xt_redist_p2p_new(xmap, mpi_double_precision, config)

    ! test communicator of redist
    IF (.NOT. communicators_are_congruent(xt_redist_get_mpi_comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", filename, __LINE__)

    DO i = 1, src_num_indices
      src_data(i) = DBLE(comm_rank * data_size + i - 1)
    END DO

    CALL check_redist(redist, src_data, dst_data, dst_index_list)

    ! clean up
    CALL xt_redist_delete(redist)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END SUBROUTINE simple_test

  ! test nonuniform numbers of send and receive partners
  SUBROUTINE nonuniform_test
    ! source index list
    INTEGER(xt_int_kind), ALLOCATABLE :: src_index_list(:), dst_index_list(:)
    DOUBLE PRECISION, ALLOCATABLE :: src_data(:), dst_data(:), ref_dst_data(:)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(Xt_redist) :: redist
    INTEGER :: i, src_num_indices, dst_num_indices

    ALLOCATE(src_index_list(comm_size), dst_index_list(comm_size), &
         src_data(comm_size), dst_data(comm_size), ref_dst_data(comm_size))
    src_num_indices = MERGE(comm_size, 0, comm_rank == 0)
    DO i = 1, src_num_indices
      src_index_list(i) = INT(i - 1, xi)
    END DO

    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices)

    ! destination index list
    dst_num_indices = comm_size
    DO i = 1, dst_num_indices
      dst_index_list(i) = INT(i - 1, xi)
    END DO

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices)

    ! xmap
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)

    ! redist_p2p
    redist = xt_redist_p2p_new(xmap, mpi_double_precision, config)

    ! test communicator of redist
    IF (.NOT. communicators_are_congruent(xt_redist_get_mpi_comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", filename, __LINE__)

    ! test exchange
    IF (comm_rank == 0) THEN
      DO i = 1, comm_size
        src_data(i) = DBLE(i - 1)
      END DO
    ELSE
      src_data(:) = -2.0d0
    END IF

    DO i = 1, comm_size
      ref_dst_data(i) = DBLE(i-1)
    END DO
    CALL check_redist(redist, src_data, dst_data, ref_dst_data)

    ! clean up
    CALL xt_redist_delete(redist)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END SUBROUTINE nonuniform_test

  ! test redist with blocks
  SUBROUTINE block_redist_test
    ! gvol_size: volume of deep ocean
    INTEGER :: ngdom, gvol_size, i, nwin, ig0, ig, j, p, qa, qb
    ! gdepth: ocean depth of an one dim. ocean
    INTEGER, ALLOCATABLE :: gdoma(:), gdomb(:), gsurfdata(:), &
         gdepth(:), ig2col_off(:), b_surfdata_ref(:), gvoldata(:), &
         src_block_offsets(:), src_block_sizes(:), dst_block_offsets(:), &
         dst_block_sizes(:), b_voldata_ref(:)
    INTEGER, ALLOCATABLE :: a_surfdata(:), b_surfdata(:), &
         a_voldata(:), b_voldata(:)
    INTEGER(xi), ALLOCATABLE :: iveca(:), ivecb(:)
    INTEGER :: ia, ib, blk_ofs_accum, gdepth_i
    TYPE(Xt_idxlist) :: idxlist_a, idxlist_b
    TYPE(xt_xmap) :: xmap
    TYPE(Xt_redist) :: redist, block_redist, block_redist2

    IF (2 * comm_size > HUGE(1_xt_int_kind)) &
         CALL test_abort('too large number of tasks', filename, __LINE__)
    ! the global index domain (1dim problem):
    ngdom = 2 * comm_size
    ! start state (index distribution) of global domain
    ALLOCATE(gdoma(ngdom), gdomb(ngdom))
    ! end state ""
    ALLOCATE(gsurfdata(ngdom), gdepth(ngdom))
    ALLOCATE(ig2col_off(ngdom)) ! offset of surface DATA within vol
    gvol_size = 0
    DO i = 1, ngdom
      gdoma(i) = i - 1
      gdomb(i) = ngdom - i
      gsurfdata(i) = 99 + i
      gdepth(i) = i
      ig2col_off(i) = gvol_size
      gvol_size = gvol_size + gdepth(i)
    END DO

    nwin = ngdom / comm_size ! my local window size of the global surface domain
    ! start of my window within global index domain (== global offset)
    ig0 = comm_rank * nwin
    IF (nwin * comm_size /= ngdom) &
         CALL test_abort("internal error", filename, __LINE__)

    ! local index
    ALLOCATE(iveca(nwin), ivecb(nwin))
    DO i = 1, nwin
      ig = ig0 + i
      iveca(i) = INT(gdoma(ig), xi)
      ivecb(i) = INT(gdomb(ig), xi)
    END DO

    idxlist_a = xt_idxvec_new(iveca, nwin)
    idxlist_b = xt_idxvec_new(ivecb, nwin)

    xmap = xt_xmap_all2all_new(idxlist_a, idxlist_b, mpi_comm_world)

    ! simple redist
    redist = xt_redist_p2p_new(xmap, mpi_integer, config)

    ! test communicator of redist
    IF (.NOT. communicators_are_congruent(xt_redist_get_mpi_comm(redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", filename, __LINE__)

    ALLOCATE(a_surfdata(nwin), b_surfdata(nwin), b_surfdata_ref(nwin))
    DO i = 1, nwin
      a_surfdata(i) = gsurfdata(iveca(i) + 1)
      b_surfdata(i) = -1
      b_surfdata_ref(i) = gsurfdata(ivecb(i) + 1)
    END DO

    CALL check_redist(redist, a_surfdata, b_surfdata, b_surfdata_ref)
    CALL xt_redist_delete(redist)

    ! generate global volume data
    ALLOCATE(gvoldata(gvol_size))
    DO i = 1, ngdom
      DO j = 1, gdepth(i)
        p = ig2col_off(i) + j
        gvoldata(p) = (i - 1) * 100 + j - 1
      END DO
    END DO

    ! generate blocks
    ALLOCATE(src_block_offsets(nwin), src_block_sizes(nwin), &
         dst_block_offsets(nwin), dst_block_sizes(nwin))
    ! we only need local size but simply oversize here
    ALLOCATE(a_voldata(gvol_size), b_voldata(gvol_size), &
         b_voldata_ref(gvol_size))
    a_voldata(:) = -1
    b_voldata_ref(:) = -1

    qa = 0
    blk_ofs_accum = 0
    DO i = 1, nwin
      ia = INT(iveca(i)) + 1
      gdepth_i = gdepth(ia)
      src_block_offsets(i) = blk_ofs_accum
      blk_ofs_accum = blk_ofs_accum + gdepth_i
      src_block_sizes(i) = gdepth_i
      p = ig2col_off(ia)
      DO j = 1, gdepth_i
        a_voldata(qa + j) = gvoldata(p + j)
      END DO
      qa = qa + gdepth_i
    END DO

    qb = 0
    blk_ofs_accum = 0
    DO i = 1, nwin
      ib = INT(ivecb(i)) + 1
      gdepth_i = gdepth(ib)
      dst_block_offsets(i) = blk_ofs_accum
      blk_ofs_accum = blk_ofs_accum + gdepth_i
      dst_block_sizes(i) = gdepth_i
      p = ig2col_off(ib)
      DO j = 1, gdepth_i
        b_voldata_ref(qb + j) = gvoldata(p + j)
      END DO
      qb = qb + gdepth_i
    END DO

    ! redist with blocks
    block_redist = xt_redist_p2p_blocks_off_custom_new(xmap, &
         src_block_offsets, src_block_sizes, nwin, &
         dst_block_offsets, dst_block_sizes, nwin, mpi_integer, config)
    ! test communicator of redist
    IF (.NOT. communicators_are_congruent(xt_redist_get_mpi_comm(block_redist), &
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", filename, __LINE__)

    CALL check_redist(block_redist, a_voldata, b_voldata, b_voldata_ref)

    ! redist with blocks but without explicit offsets:
    block_redist2 = xt_redist_p2p_blocks_custom_new(xmap, &
         src_block_sizes, nwin, dst_block_sizes, nwin, mpi_integer, config)
    ! test communicator of redist

    IF (.NOT. communicators_are_congruent(xt_redist_get_mpi_comm(block_redist2),&
         mpi_comm_world)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", filename, __LINE__)

    CALL check_redist(block_redist2, a_voldata, b_voldata, b_voldata_ref)

    ! cleanup
    CALL xt_redist_delete(block_redist2)
    CALL xt_redist_delete(block_redist)
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(idxlist_a)
    CALL xt_idxlist_delete(idxlist_b)
  END SUBROUTINE block_redist_test

END PROGRAM test_redist_p2p_parallel
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
