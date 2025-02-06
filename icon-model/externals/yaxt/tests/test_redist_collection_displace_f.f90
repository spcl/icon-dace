!>
!! @file test_redist_collection_displace_f.f90
!! @brief Fortran cache displacement test of redist_collection class
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
MODULE redist_collection_displace
  USE mpi
  USE ftest_common, ONLY: test_abort, cmp_arrays
  USE yaxt, ONLY: &
       xt_xmap, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_collection_new, &
       Xt_redist_copy, xt_redist_delete, xt_redist_s_exchange, &
       xt_request, xt_redist_a_exchange, xt_config
  ! older PGI compilers do not handle generic interface correctly
#if defined __PGI && (__PGIC__ < 12 || (__PGIC__ ==  12 && __PGIC_MINOR__ <= 10))
  USE xt_redist_base, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
#endif
  ! and when taking the slice address and the optimizer is on,
  ! some other random failure occurs even with very recent compiler versions
#if defined __PGI && (__PGIC__ <= 21 || __PGIC__ == 22 && __PGIC_MINOR__ <= 5)
#undef HAVE_FC_C_LOC_OF_SLICE
#endif
  USE test_redist_common, ONLY: build_odd_selection_xmap, &
       check_wait_request
  USE iso_c_binding, ONLY: c_ptr
#include "xt_slice_c_loc.inc"
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: cache_size = 16, cache_overrun = 2
  INTEGER, PARAMETER :: num_slice = 3, dst_step = 2
  INTEGER, PARAMETER :: src_slice_len = 5
  INTEGER, PARAMETER :: dst_slice_len &
       = (src_slice_len + dst_step - 1)/dst_step
  CHARACTER(len=*), PARAMETER :: &
       filename = 'test_redist_collection_displace_f.f90'
  CHARACTER(len=*), PARAMETER :: err_msg(2) = &
       (/ "error on xt_redist_s_exchange", "error on xt_redist_a_exchange" /)
  PUBLIC :: test_displacement_variations
CONTAINS
  ! test with one redist used three times (with different input
  ! data displacements until the cache is full)
  ! set up data
  SUBROUTINE test_displacement_variations(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redists(num_slice), redist_coll, &
         redist_coll_copy

    xmap = build_odd_selection_xmap(src_slice_len, comm)
    redist = xt_redist_p2p_new(xmap, mpi_double_precision, config)

    CALL xt_xmap_delete(xmap)

    ! generate redist_collection
    redists = redist

    redist_coll = xt_redist_collection_new(redists, num_slice, &
         cache_size, comm, config)

    CALL xt_redist_delete(redist)

    CALL run_displacement_check(redist_coll)
    redist_coll_copy = xt_redist_copy(redist_coll)
    CALL run_displacement_check(redist_coll_copy)

    ! clean up
    CALL xt_redist_delete(redist_coll)
    CALL xt_redist_delete(redist_coll_copy)
  END SUBROUTINE test_displacement_variations

#ifndef HAVE_FC_PTR_BOUND_REMAP
  SUBROUTINE ptr_bind(p, ub, a)
    DOUBLE PRECISION, POINTER :: p(:,:)
    INTEGER, INTENT(in) :: ub(2)
    DOUBLE PRECISION, TARGET :: a(ub(1), ub(2))
    p => a
  END SUBROUTINE ptr_bind
#endif

  SUBROUTINE run_displacement_check(redist_coll)
    TYPE(xt_redist), INTENT(in) :: redist_coll
    INTEGER :: i, j, k
    DOUBLE PRECISION, TARGET :: &
         src(src_slice_len * (num_slice+1) + cache_size + cache_overrun), &
         dst(dst_slice_len * (num_slice+1) + cache_size + cache_overrun)
    DOUBLE PRECISION, POINTER :: src_data(:, :), dst_data(:, :)
    DOUBLE PRECISION, POINTER :: src_data_(:), dst_data_(:)
    TYPE(c_ptr) :: src_data_p(num_slice), dst_data_p(num_slice)
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(dst_slice_len, num_slice) = &
         RESHAPE((/ ((DBLE(i + j * src_slice_len), &
         &            i = 1, src_slice_len, dst_step), &
         &           j = 0, num_slice - 1) /), &
         &       (/ dst_slice_len, num_slice /))
    TYPE(xt_request) :: request
    INTEGER :: iexch

#ifdef HAVE_FC_PTR_BOUND_REMAP
    src_data(1:src_slice_len, 1:num_slice) => src(1:src_slice_len*num_slice)
#else
    CALL ptr_bind(src_data, (/ src_slice_len, num_slice /), src)
#endif
    DO j = 1, num_slice
      DO i = 1, src_slice_len
        src_data(i, j) = DBLE(i + (j - 1) * src_slice_len)
      END DO
    END DO

    src_data_ => src(src_slice_len*num_slice+1:)

#ifdef HAVE_FC_PTR_BOUND_REMAP
    dst_data(1:dst_slice_len, 1:num_slice) => dst(1:dst_slice_len*num_slice)
#else
    CALL ptr_bind(dst_data, (/ dst_slice_len, num_slice /), dst)
#endif

    dst_data_ => dst(dst_slice_len*num_slice+1:)

    DO i = 1, num_slice - 1
      XT_SLICE_C_LOC(src_data(:, i), src_data_p(i))
      XT_SLICE_C_LOC(dst_data(:, i), dst_data_p(i))
    END DO

    ! test exchange
    DO k = 1, cache_size + cache_overrun
      src_data_(k:k+src_slice_len-1) = src_data(:,num_slice)


      XT_SLICE_C_LOC(src_data_(k:k+src_slice_len-1), src_data_p(3))
      XT_SLICE_C_LOC(dst_data_(k:k+dst_slice_len-1), dst_data_p(3))

      DO iexch = 1, 2
        dst = -1.0d0
        IF (iexch == 1) THEN
          CALL xt_redist_s_exchange(redist_coll, num_slice, src_data_p, &
               dst_data_p)
        ELSE
          CALL xt_redist_a_exchange(redist_coll, num_slice, src_data_p, &
               dst_data_p, request)
          CALL check_wait_request(request, filename, __LINE__)
        ENDIF
        IF (cmp_arrays(ref_dst_data(:, 1:num_slice-1), &
             &         dst_data(:, 1:num_slice-1)) &
             .OR. cmp_arrays(ref_dst_data(:,num_slice), &
             &               dst_data_(k:k+dst_slice_len-1))) &
           CALL test_abort(err_msg(iexch), filename, __LINE__)
      ENDDO
    END DO
  END SUBROUTINE run_displacement_check

END MODULE redist_collection_displace
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
