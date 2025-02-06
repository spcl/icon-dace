!>
!! @file test_redist_repeat_f.f90
!! @brief Fortran test of redist_repeatection_static class
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
PROGRAM test_redist_repeat
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, cmp_arrays
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_xmap, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_p2p_off_new, &
       xt_redist_repeat_new, xt_redist_delete, &
       xt_redist_s_exchange, xt_redist_a_exchange, &
       xt_request, xt_config, xt_config_delete
#if defined __PGI && ( __PGIC__ == 15 || __PGIC__ == 14 )
  USE xt_redist_real_sp, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
  USE xt_redist_real_dp, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
#endif
  USE test_redist_common, ONLY: build_odd_selection_xmap, check_redist, &
       check_wait_request, redist_exchanger_option
  USE iso_c_binding, ONLY: c_int
  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER :: filename = 'test_redist_repeat_f.f90'
  CHARACTER(len=*), PARAMETER :: exch1name(2) = &
       (/ "xt_redist_s_exchange1", "xt_redist_a_exchange1" /)
  TYPE(xt_config) :: config

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  config = redist_exchanger_option()

  CALL simple_test(mpi_comm_world, config)
  CALL test_repeated_redist(mpi_comm_world, config)
  CALL test_repeated_redist_with_gap(mpi_comm_world, config)
  CALL test_repeated_overlapping_redist(mpi_comm_world, config)
  CALL test_repeated_redist_asym(mpi_comm_world, config)

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_config_delete(config)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE simple_test(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    ! general test with one redist
    ! set up data
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_repeat
    INTEGER, PARAMETER :: src_slice_len = 5, dst_slice_len = 3
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(dst_slice_len) &
         = (/ 1.0d0, 3.0d0, 5.0d0 /), &
         src_data(src_slice_len) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)
    DOUBLE PRECISION :: dst_data(dst_slice_len)
    INTEGER(mpi_address_kind) :: src_extent, dst_extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int) :: displacements(1) = 0
    INTEGER :: ierror

    xmap = build_odd_selection_xmap(src_slice_len, comm)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    CALL mpi_get_address(src_data(1), base_address, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort('mpi_get_address failed', filename, __LINE__)
    CALL mpi_get_address(src_data(2), temp_address, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort('mpi_get_address failed', filename, __LINE__)
    src_extent = (temp_address - base_address) * src_slice_len
    CALL mpi_get_address(dst_data(1), base_address, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort('mpi_get_address failed', filename, __LINE__)
    CALL mpi_get_address(dst_data(2), temp_address, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort('mpi_get_address failed', filename, __LINE__)
    dst_extent = (temp_address - base_address) * dst_slice_len

    ! generate redist_repeat
    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, 1, &
                                         displacements, config)

    CALL xt_redist_delete(redist)

    ! test exchange
    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)

    ! clean up
    CALL xt_redist_delete(redist_repeat)
  END SUBROUTINE simple_test

  SUBROUTINE test_repeated_redist_ds1(redist_repeat)
    TYPE(xt_redist), INTENT(in) :: redist_repeat
    INTEGER :: i, j
    DOUBLE PRECISION, PARAMETER :: src_data(5, 3) = RESHAPE((/&
         (DBLE(i), i = 1, 15)/), (/ 5, 3 /))
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(3, 3) &
         = RESHAPE((/ ((DBLE(i + j), i = 1,5,2), j = 0,10,5) /), (/ 3, 3 /))
    DOUBLE PRECISION :: dst_data(3, 3)

    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)
  END SUBROUTINE test_repeated_redist_ds1

#ifdef __PGI
#  define NO_2D_PARAM
#elif defined(_CRAYFTN)
#  if _RELEASE_MAJOR < 8 || (_RELEASE_MAJOR == 8 && _RELEASE_MINOR < 7)
#    define NO_2D_PARAM
#  endif
#endif

  SUBROUTINE test_repeated_redist_ds1_with_gap(redist_repeat)
    TYPE(xt_redist), INTENT(in) :: redist_repeat
    INTEGER :: i, j
    DOUBLE PRECISION, PARAMETER :: src_data(5, 5) = RESHAPE((/&
         (DBLE(i), i = 1, 25)/), (/ 5, 5 /))
    DOUBLE PRECISION :: dst_data(3, 5)
#ifdef NO_2D_PARAM
    DOUBLE PRECISION :: ref_dst_data(3, 5)
    ref_dst_data &
         = RESHAPE((/ ((DBLE((i + j)*MOD(j+1,2)-MOD(j,2)), i = 1,5,2), &
         j = 0,20,5) /), (/ 3, 5 /))
#else
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(3, 5) &
         = RESHAPE((/ ((DBLE((i + j)*MOD(j+1,2)-MOD(j,2)), i = 1,5,2), j = 0,20,5) /), (/ 3, 5 /))
#endif
    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)
  END SUBROUTINE test_repeated_redist_ds1_with_gap

  SUBROUTINE test_repeated_redist_ds2(redist_repeat)
    TYPE(xt_redist), INTENT(in) :: redist_repeat
    INTEGER :: i, j
    DOUBLE PRECISION, PARAMETER :: src_data(5, 3) = RESHAPE((/&
         (DBLE(i), i = 20, 34)/), (/ 5, 3 /))
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(3, 3) &
         = RESHAPE((/ ((DBLE(i + j), i = 1,5,2), j = 19,33,5) /), (/ 3, 3 /))
    DOUBLE PRECISION :: dst_data(3, 3)

    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)
  END SUBROUTINE test_repeated_redist_ds2

  SUBROUTINE test_repeated_redist(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    ! test with one redist used three times (with two different input data
    ! displacements -> test of cache) (with default cache size)
    ! set up data
    INTEGER, PARAMETER :: num_slice = 3
    INTEGER, PARAMETER :: src_slice_len = 5
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_repeat
    INTEGER(mpi_address_kind) :: src_extent, dst_extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int), PARAMETER :: &
         displacements(3)= (/ 0_c_int, 1_c_int, 2_c_int /)
    DOUBLE PRECISION, TARGET :: src_template(5, 3), dst_template(3, 3)
    INTEGER :: ierror

    xmap = build_odd_selection_xmap(src_slice_len, comm)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_repeat
    CALL mpi_get_address(src_template(1,1), base_address, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort('mpi_get_address failed', filename, __LINE__)
    CALL mpi_get_address(src_template(1,2), temp_address, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort('mpi_get_address failed', filename, __LINE__)
    src_extent = temp_address - base_address
    CALL mpi_get_address(dst_template(1,1), base_address, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort('mpi_get_address failed', filename, __LINE__)
    CALL mpi_get_address(dst_template(1,2), temp_address, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort('mpi_get_address failed', filename, __LINE__)
    dst_extent = temp_address - base_address

    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, &
         num_slice, displacements, config)
    CALL xt_redist_delete(redist)

    ! test exchange
    CALL test_repeated_redist_ds1(redist_repeat)
    ! test exchange
    CALL test_repeated_redist_ds2(redist_repeat)
    ! clean up
    CALL xt_redist_delete(redist_repeat)
  END SUBROUTINE test_repeated_redist

  SUBROUTINE test_repeated_redist_asym(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    ! test asymmetric variant of redist_repeat

    INTEGER, PARAMETER :: num_slice = 3
    INTEGER, PARAMETER :: src_slice_len = 5
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_repeat
    INTEGER(mpi_address_kind) :: src_extent, dst_extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int) :: src_displacements(3), dst_displacements(3)
    INTEGER :: i, ierror
    DOUBLE PRECISION, PARAMETER :: &
         ref_dst_data(3, 3) = RESHAPE([ 6.0d0, 8.0d0, 10.0d0, 11.0d0, 13.0d0, &
         &                             15.0d0, 1.0d0,  3.0d0,  5.0d0 ], [3,3] )
    DOUBLE PRECISION, TARGET :: dst_data(3, 3)
    DOUBLE PRECISION, TARGET, SAVE :: &
         src_data(5, 3) = RESHAPE([(DBLE(i), i = 1, 15)], [5,3])
    INTEGER, PARAMETER :: dp = KIND(src_data)

    ! xmap: [1,2,3,4,5] -> [1,3,5]
    xmap = build_odd_selection_xmap(src_slice_len, comm)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_repeat:
    CALL mpi_get_address(src_data(1,1), base_address, ierror)
    CALL mpi_get_address(src_data(1,2), temp_address, ierror)
    src_extent = temp_address - base_address
    CALL mpi_get_address(dst_data(1,1), base_address, ierror)
    CALL mpi_get_address(dst_data(1,2), temp_address, ierror)
    dst_extent = temp_address - base_address

    ! repeated redist parameters:
    src_displacements = [0,1,2]
    dst_displacements = [2,0,1]

    ! connect to explicit shape interface:
    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, &
         num_slice, src_displacements, dst_displacements, config)
    dst_data = -1.0_dp
    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)
    CALL xt_redist_delete(redist_repeat)

    ! connect to assumed shape interface:
    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, &
         src_displacements, dst_displacements, config)
    dst_data = -1.0_dp
    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)
    CALL xt_redist_delete(redist_repeat)

    CALL xt_redist_delete(redist)
  END SUBROUTINE test_repeated_redist_asym

  SUBROUTINE test_repeated_redist_with_gap(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    ! test with one redist used three times (with two different input data
    ! displacements -> test of cache) (with default cache size)
    ! set up data
    INTEGER, PARAMETER :: num_slice = 3
    INTEGER, PARAMETER :: src_slice_len = 5
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_repeat
    INTEGER(mpi_address_kind) :: src_extent, dst_extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int), PARAMETER :: displacements(3) = (/0,2,4/)
    DOUBLE PRECISION, TARGET :: src_template(5, 3), dst_template(3, 3)
    INTEGER :: ierror

    xmap = build_odd_selection_xmap(src_slice_len, comm)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_repeat
    CALL mpi_get_address(src_template(1,1), base_address, ierror)
    CALL mpi_get_address(src_template(1,2), temp_address, ierror)
    src_extent = temp_address - base_address
    CALL mpi_get_address(dst_template(1,1), base_address, ierror)
    CALL mpi_get_address(dst_template(1,2), temp_address, ierror)
    dst_extent = temp_address - base_address

    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, &
         num_slice, displacements, config)
    CALL xt_redist_delete(redist)

    ! test exchange
    CALL test_repeated_redist_ds1_with_gap(redist_repeat)
    ! clean up
    CALL xt_redist_delete(redist_repeat)
  END SUBROUTINE test_repeated_redist_with_gap

  SUBROUTINE test_repeated_overlapping_redist(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    ! test with one redist used three times (with two different input data
    ! displacements -> test of cache) (with default cache size)
    ! set up data
    INTEGER, PARAMETER :: npt = 9, selection_len = 6
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_repeat
    INTEGER(mpi_address_kind) :: src_extent, dst_extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int), PARAMETER :: displacements(2) = (/ 0_c_int, 1_c_int /)
    INTEGER :: i, j, ierror
    INTEGER, PARAMETER :: src_pos(npt) = (/ (i, i=1,npt) /), &
         dst_pos(npt) = (/ (2*i, i = 0, npt-1) /)
    DOUBLE PRECISION :: src_data(npt), dst_data(npt)
#if __INTEL_COMPILER >= 1600 && __INTEL_COMPILER <= 1602 || defined __PGI
    DOUBLE PRECISION :: ref_dst_data(npt)
#else
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(npt) &
         = (/ ((DBLE(((2-j)*3+i+101)*((ABS(j)+j)/ABS(j+MAX(1,j))) &
         &           +(j-1-ABS(j-1))/2), &
         &     i=1,3 ),j=2,0,-1) /)
#endif
    DOUBLE PRECISION, TARGET :: src_template(2), dst_template(2)
    INTEGER :: iexch
    TYPE(xt_request) :: request(2)

    xmap = build_odd_selection_xmap(selection_len, comm)

    redist = xt_redist_p2p_off_new(xmap, src_pos, dst_pos, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! init data
#if __INTEL_COMPILER >= 1600 && __INTEL_COMPILER <= 1602 || defined __PGI
    DO j = 2, 0, -1
      DO i = 1, 3
        ref_dst_data(i + (2-j)*3) = DBLE(((2-j)*3+i+101)*((ABS(j)+j)/ABS(j+1)) &
             &                           +(j-1-ABS(j-1))/2)
      END DO
    END DO
#endif
    DO i = 1, npt
      src_data(i) = 1.0d2 + DBLE(i)
    END DO

    DO iexch = 1, 2

      dst_data = -1.0d0

      ! test individual redists
      IF (iexch == 1) THEN
        CALL xt_redist_s_exchange(redist, src_data, dst_data)
        CALL xt_redist_s_exchange(redist, src_data(2:), dst_data(2:))
      ELSE
        CALL xt_redist_a_exchange(redist, src_data, dst_data, request(1))
        CALL xt_redist_a_exchange(redist, src_data(2:), dst_data(2:), request(2))
        CALL check_wait_request(request(1), filename, __LINE__)
        CALL check_wait_request(request(2), filename, __LINE__)
      ENDIF
      ! check individual redists to have desired effect
      IF (cmp_arrays(dst_data, ref_dst_data)) &
           CALL test_abort("error in "//exch1name(iexch), filename,__LINE__)
    ENDDO
    dst_data = -1.0d0
    ! generate redist_repeat
    CALL mpi_get_address(src_template(1), base_address, ierror)
    CALL mpi_get_address(src_template(2), temp_address, ierror)
    src_extent = temp_address - base_address
    CALL mpi_get_address(dst_template(1), base_address, ierror)
    CALL mpi_get_address(dst_template(2), temp_address, ierror)
    dst_extent = temp_address - base_address

    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, &
         displacements, config)
    CALL xt_redist_delete(redist)

    ! test exchange
    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)
    ! clean up
    CALL xt_redist_delete(redist_repeat)
  END SUBROUTINE test_repeated_overlapping_redist

END PROGRAM test_redist_repeat
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
