!>
!! @file test_redist_collection_static_f.f90
!! @brief Fortran test of redist_collection_static class
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
PROGRAM test_redist_collection_static
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, &
       xt_xmap, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_collection_static_new, &
       xt_redist_delete, xt_config, xt_config_delete
  USE test_redist_common, ONLY: build_odd_selection_xmap, check_redist, &
       redist_exchanger_option
  IMPLICIT NONE
  TYPE(xt_config) :: config
  CHARACTER(len=*), PARAMETER :: &
       filename = 'test_redist_collection_static_f.f90'
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  config = redist_exchanger_option()

  CALL simple_test(mpi_comm_world, config)
  CALL test_repeated_redist(mpi_comm_world, config)

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_config_delete(config)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE simple_test(comm, config)
    ! general test with one redist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    ! set up data
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_coll
    INTEGER, PARAMETER :: src_slice_len = 5, dst_slice_len = 3
    DOUBLE PRECISION, PARAMETER :: &
         src_data(src_slice_len) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /), &
         ref_dst_data(dst_slice_len) = (/ 1.0d0, 3.0d0, 5.0d0 /)
    DOUBLE PRECISION :: dst_data(dst_slice_len)
    INTEGER(mpi_address_kind), PARAMETER :: &
         displacements(1) = 0_mpi_address_kind

    xmap = build_odd_selection_xmap(src_slice_len, comm)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_collection
    redist_coll = xt_redist_collection_static_new((/ redist /), 1, &
         displacements, displacements, comm, config)

    CALL xt_redist_delete(redist)

    ! test exchange
    CALL check_redist(redist_coll, src_data, dst_data, ref_dst_data)

    ! clean up
    CALL xt_redist_delete(redist_coll)
  END SUBROUTINE simple_test

  SUBROUTINE test_repeated_redist_ds1(redist_coll)
    TYPE(xt_redist), INTENT(in) :: redist_coll
    INTEGER :: i, j
    DOUBLE PRECISION, PARAMETER :: src_data(5, 3) &
         = RESHAPE((/ (DBLE(i), i = 1, 15)/), (/ 5, 3 /)), &
         ref_dst_data(3, 3) &
         = RESHAPE((/ ((DBLE(i + j), i = 1,5,2), j = 0,10,5) /), (/ 3, 3 /))
    DOUBLE PRECISION :: dst_data(3, 3)

    CALL check_redist(redist_coll, src_data, dst_data, ref_dst_data)
  END SUBROUTINE test_repeated_redist_ds1

  SUBROUTINE test_repeated_redist_ds2(redist_coll)
    TYPE(xt_redist), INTENT(in) :: redist_coll
    INTEGER :: i, j
    DOUBLE PRECISION, PARAMETER :: src_data(5, 3) = RESHAPE((/&
         (DBLE(i), i = 20, 34)/), (/ 5, 3 /)), &
         ref_dst_data(3, 3) &
         = RESHAPE((/ ((DBLE(i + j), i = 1,5,2), j = 19,33,5) /), (/ 3, 3 /))
    DOUBLE PRECISION, SAVE :: dst_data(3, 3)

    CALL check_redist(redist_coll, src_data, dst_data, ref_dst_data)
  END SUBROUTINE test_repeated_redist_ds2

  SUBROUTINE test_repeated_redist(comm, config)
    ! test with one redist used three times (with two different input data
    ! displacements -> test of cache) (with default cache size)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    ! set up data
    INTEGER, PARAMETER :: num_slice = 3
    INTEGER, PARAMETER :: src_slice_len = 5
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redists(num_slice), redist_coll
    INTEGER(mpi_address_kind) :: src_displacements(num_slice), &
         dst_displacements(num_slice), src_base, dst_base, temp
    DOUBLE PRECISION, TARGET :: src_template(5, 3), dst_template(3, 3)
    INTEGER :: i, ierror

    xmap = build_odd_selection_xmap(src_slice_len, comm)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_collection
    redists = redist
    src_displacements(1) = 0_mpi_address_kind
    dst_displacements(1) = 0_mpi_address_kind
    CALL mpi_get_address(src_template(1, 1), src_base, ierror)
    CALL mpi_get_address(dst_template(1, 1), dst_base, ierror)
    DO i = 2, num_slice
      CALL mpi_get_address(src_template(1, i), temp, ierror)
      src_displacements(i) = temp - src_base
      CALL mpi_get_address(dst_template(1, i), temp, ierror)
      dst_displacements(i) = temp - dst_base
    END DO

    redist_coll = xt_redist_collection_static_new(redists, num_slice, &
         src_displacements, dst_displacements, comm, config)
    CALL xt_redist_delete(redist)

    ! test exchange
    CALL test_repeated_redist_ds1(redist_coll)
    ! test exchange with changed displacements
    CALL test_repeated_redist_ds2(redist_coll)
    ! test exchange with original displacements
    CALL test_repeated_redist_ds1(redist_coll)
    ! clean up
    CALL xt_redist_delete(redist_coll)
  END SUBROUTINE test_repeated_redist

END PROGRAM test_redist_collection_static
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
