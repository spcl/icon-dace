!>
!! @file test_sort_f.f90
!! @brief Test Fortran interface to yaxt sorting functions
!!
!! @copyright Copyright  (C)  2021 Jörg Behrens <behrens@dkrz.de>
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
PROGRAM test_sort
  USE iso_c_binding, ONLY: c_int, c_double, c_int32_t
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, &
       xt_sort_index, xt_sort_int, xt_sort_permutation, &
       xi => xt_int_kind
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, &
       random_fill, crc32, permute
  IMPLICIT NONE
  INTEGER, PARAMETER :: range_size=500000
  CHARACTER(len=*), PARAMETER :: filename = 'test_sort_f.f90'
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  CALL test_int_sort()
  CALL test_xt_int_sort()
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE test_int_sort()
    INTEGER(c_int), ALLOCATABLE :: a(:), permutation(:)
    REAL(c_double) :: rnd
    INTEGER :: n, i
    INTEGER(c_int32_t) :: crc_before, crc_after
    LOGICAL :: err_found
    CALL RANDOM_NUMBER(rnd)
    n = NINT(rnd * range_size)
    ALLOCATE(a(n), permutation(n))
    CALL random_fill(a)
    CALL xt_sort_int(a)
    err_found = .FALSE.
    DO i = 2, n
      err_found = err_found .OR. a(i-1) > a(i)
    END DO
    IF (err_found) CALL test_abort("discontinuity found!", filename, __LINE__)
    CALL random_fill(a)
    DO i = 1, n
      permutation(i) = i
    END DO
    crc_before = crc32(a)
    CALL xt_sort_permutation(a, permutation)
    DO i = 2, n
      err_found = err_found .OR. a(i-1) > a(i)
    END DO
    IF (err_found) CALL test_abort("discontinuity found!", filename, __LINE__)
    CALL permute(a, permutation)
    crc_after = crc32(a)
    IF (crc_after /= crc_before) &
         CALL test_abort("error in data reconstruction!", filename, __LINE__)
  END SUBROUTINE test_int_sort

  SUBROUTINE test_xt_int_sort()
    INTEGER(c_int), ALLOCATABLE :: permutation(:)
    INTEGER(xi), ALLOCATABLE :: b(:)
    REAL(c_double) :: rnd
    INTEGER :: n, i
    INTEGER(c_int32_t) :: crc_before, crc_after
    LOGICAL :: err_found
    CALL RANDOM_NUMBER(rnd)
    n = NINT(rnd * range_size)
    ALLOCATE(b(n), permutation(n))
    CALL random_fill(b)
    crc_before = crc32(b)
    DO i = 1, n
      permutation(i) = i
    END DO
    CALL xt_sort_index(b, permutation, reset_positions=.FALSE.)
    err_found = .FALSE.
    DO i = 2, n
      err_found = err_found .OR. b(i-1) > b(i)
    END DO
    IF (err_found) CALL test_abort("discontinuity found!", filename, __LINE__)
    CALL permute(b, permutation)
    crc_after = crc32(b)
    IF (crc_after /= crc_before) &
         CALL test_abort("error in data reconstruction!", filename, __LINE__)
  END SUBROUTINE test_xt_int_sort

END PROGRAM test_sort
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
! license-default: "bsd"
! End:
!
