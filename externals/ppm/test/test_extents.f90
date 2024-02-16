!>
!! @file test_extents.f90
!! @brief unit tests for extents data type
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!!
!! @author Thomas Jahns <jahns@dkrz.de>
!! @version 1.0
!
! Maintainer: Thomas Jahns <jahns@dkrz.de>
! URL: https://www.dkrz.de/redmine/projects/scales-ppm
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
MODULE extent_tests
  USE ppm_std_type_kinds, ONLY: i4, i8
  USE ppm_extents, ONLY: ASSIGNMENT(=), extent, extent_size, iinterval, &
       OPERATOR(==), OPERATOR(/=), extents_do_intersect, extent_intersect, &
       char, ext2s_len, extent_i8, iinterval_i8, ext2s_len_i8
  USE ppm_base, ONLY: abort_ppm
  USE ppm_random, ONLY: irand, irand8, irandr
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_extent_i4, test_extent_i8

  CHARACTER(len=*), PARAMETER :: filename = 'test_extents.f90'
CONTAINS
  SUBROUTINE setup_extent_test_errmsg(eq_warn, overlap_pred_warn, &
       overlap_size_warn)
    CHARACTER(len=*), INTENT(out) :: eq_warn, overlap_pred_warn, &
         overlap_size_warn
    INTEGER :: tid
#ifdef _OPENMP
    tid = omp_get_thread_num()
#else
    tid = 0
#endif
#if defined __GNUC__ && __GNUC__ == 7
    ! work-around for GCC internal write race bug
!$omp critical
#endif
    WRITE (eq_warn, '(a,i0)') 'equality comparison failed, tid=', tid
    WRITE (overlap_pred_warn, '(a,i0)') 'overlap assertion failed, tid=', tid
    WRITE (overlap_size_warn, '(a,i0)') 'overlap size mismatch, tid=', tid
#if defined __GNUC__ && __GNUC__ == 7
!$omp end critical
#endif
  END SUBROUTINE setup_extent_test_errmsg

  SUBROUTINE test_extent_i4
    INTEGER(i4), PARAMETER :: rmod=1000, rstart_mod=10000, num_tests=10000
    TYPE(extent) :: a(2)
    TYPE(iinterval) :: b(2)
    CHARACTER(len=ext2s_len) :: b_str(2), b_itsect_str
    INTEGER :: i, tid
    INTEGER(i4) :: l, m
    CHARACTER(len=64) :: eq_warn, overlap_pred_warn, overlap_size_warn, &
         min_intersect_eq_warn

#ifdef _OPENMP
    tid = omp_get_thread_num()
#else
    tid = 0
#endif

    CALL setup_extent_test_errmsg(eq_warn, overlap_pred_warn, overlap_size_warn)
!$omp do
    DO i = 1, num_tests
      a(1) = extent(0_i4, 0_i4)
      DO WHILE(extent_size(a(1)) == 0)
        a(1)%first = irand()
        l = irand()
        ! prevent overflow
        IF (a(1)%first < 0 .AND. l < 0) l = MOD(l, HUGE(l) + a(1)%first + 1)
        IF (a(1)%first > 0 .AND. l > 0) l = MOD(l, HUGE(l) - a(1)%first + 1)
        a(1)%size = l
      END DO
      b(1) = a(1)
      a(2) = b(1)
      IF (a(1) /= a(2)) THEN
        CALL abort_ppm(TRIM(eq_warn), filename, __LINE__)
      END IF
      l = irandr(iinterval(-rstart_mod,rstart_mod))
      m = irandr(iinterval(0, 2*rmod-1))
      b(1) = iinterval(l, l + m)
      b(2) = iinterval(l, l + irandr(iinterval(0, 2*rmod-1)))
      IF (.NOT. extents_do_intersect(b(1), b(2))) THEN
        CALL abort_ppm(TRIM(overlap_pred_warn), filename, __LINE__)
      END IF
      b_str = CHAR(b)
      b_itsect_str = CHAR(extent_intersect(b(1), b(2)))
      IF (extent_intersect(b(1), b(2)) &
           /= MERGE(b(1), b(2), b(1)%last < b(2)%last)) THEN
        WRITE (min_intersect_eq_warn, '(7a,i0)') &
             'extent_intersect(', TRIM(b_str(1)), ',', TRIM(b_str(2)), &
             ') = ', TRIM(b_itsect_str), '; tid=', tid
        CALL abort_ppm(TRIM(min_intersect_eq_warn), filename, __LINE__)
      END IF
      a = b
      IF (.NOT. extents_do_intersect(a(1), a(2))) THEN
        CALL abort_ppm(TRIM(overlap_pred_warn), filename, __LINE__)
      END IF
      IF (   MIN(extent_size(a(1)), extent_size(a(2))) &
        & /= extent_size(extent_intersect(a(1), a(2)))) THEN
        CALL abort_ppm(TRIM(overlap_size_warn), filename, __LINE__)
      END IF
    END DO
!$omp end do
  END SUBROUTINE test_extent_i4

  SUBROUTINE test_extent_i8
    INTEGER(i8), PARAMETER :: rmod=1000_i8, rstart_mod=10000_i8
    INTEGER, PARAMETER :: num_tests=10000
    TYPE(extent_i8) :: a(2)
    TYPE(iinterval_i8) :: b(2)
    CHARACTER(len=ext2s_len_i8) :: b_str(2), b_itsect_str
    INTEGER :: i, tid
    INTEGER(i8) :: l, m
    CHARACTER(len=96) :: eq_warn, overlap_pred_warn, overlap_size_warn, &
         min_intersect_eq_warn

#ifdef _OPENMP
    tid = omp_get_thread_num()
#else
    tid = 0
#endif

    CALL setup_extent_test_errmsg(eq_warn, overlap_pred_warn, overlap_size_warn)
!$omp do
    DO i = 1, num_tests
      a(1) = extent_i8(0_i8, 0_i8)
      DO WHILE(extent_size(a(1)) == 0_i8)
        a(1)%first = irand8()
        l = irand8()
        ! prevent overflow
        IF (a(1)%first < 0_i8 .AND. l < 0_i8) &
             l = MOD(l, HUGE(l) + a(1)%first + 1_i8)
        IF (a(1)%first > 0_i8 .AND. l > 0_i8) &
             l = MOD(l, HUGE(l) - a(1)%first + 1_i8)
        a(1)%size = l
      END DO
      b(1) = a(1)
      a(2) = b(1)
      IF (a(1) /= a(2)) THEN
        CALL abort_ppm(TRIM(eq_warn), filename, __LINE__)
      END IF
      l = irandr(iinterval_i8(-rstart_mod,rstart_mod))
      m = irandr(iinterval_i8(0_i8, 2_i8 * rmod - 1_i8))
      b(1) = iinterval_i8(l, l + m)
      b(2) = iinterval_i8(l, l + irandr(iinterval_i8(0_i8, 2_i8 * rmod - 1_i8)))
      IF (.NOT. extents_do_intersect(b(1), b(2))) THEN
        CALL abort_ppm(TRIM(overlap_pred_warn), filename, __LINE__)
      END IF
      b_str = CHAR(b)
      b_itsect_str = CHAR(extent_intersect(b(1), b(2)))
      IF (extent_intersect(b(1), b(2)) &
           /= MERGE(b(1), b(2), b(1)%last < b(2)%last)) THEN
        WRITE (min_intersect_eq_warn, '(7a,i0)') &
             'extent_intersect(', TRIM(b_str(1)), ',', TRIM(b_str(2)), &
             ') = ', TRIM(b_itsect_str), '; tid=', tid
        CALL abort_ppm(TRIM(min_intersect_eq_warn), filename, __LINE__)
      END IF
      a = b
      IF (.NOT. extents_do_intersect(a(1), a(2))) THEN
        CALL abort_ppm(TRIM(overlap_pred_warn), filename, __LINE__)
      END IF
      IF (   MIN(extent_size(a(1)), extent_size(a(2))) &
        & /= extent_size(extent_intersect(a(1), a(2)))) THEN
        CALL abort_ppm(TRIM(overlap_size_warn), filename, __LINE__)
      END IF
    END DO
!$omp end do
  END SUBROUTINE test_extent_i8

END MODULE extent_tests

PROGRAM test_extents
  USE extent_tests, ONLY: test_extent_i4, test_extent_i8
  USE ppm_base, ONLY: ppm_default_comm
  USE ppm_random, ONLY: initialize_irand
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
  IMPLICIT NONE
  INTEGER :: rand_seed, tid
!$omp parallel private(tid, rand_seed)
  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
#ifdef _OPENMP
  tid = omp_get_thread_num()
#else
  tid = 0
#endif
  PRINT '(2(a,i0))', 'thread id=', tid, ', random seed=', rand_seed

  CALL test_extent_i4

  CALL test_extent_i8

!$omp end parallel
END PROGRAM test_extents
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
