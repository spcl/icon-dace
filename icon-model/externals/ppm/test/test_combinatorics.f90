!>
!! @file test_combinatorics.f90
!! @brief test combinatorial routines
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
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
PROGRAM test_combinatorics
  USE ppm_base, ONLY: ppm_default_comm, assertion, abort_ppm
  USE ppm_std_type_kinds, ONLY: i8, i4
  USE ppm_compare, ONLY: cmp_i4
  USE ppm_sort, ONLY: qsort_r_mt
  USE ppm_extents, ONLY: iinterval
  USE ppm_random, ONLY: initialize_irand, irand, irandr, lrand, a_rand_mt
  USE ppm_combinatorics, ONLY: permute, is_permutation, combination
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
  IMPLICIT NONE
  INTEGER(i4), ALLOCATABLE :: a(:), a_ref(:), a_selection(:), a_selection_i(:)
  INTEGER, PARAMETER :: max_a_size=1000000
  INTEGER :: rand_seed, i, n, tid, a_selection_size
  LOGICAL :: p
  LOGICAL, ALLOCATABLE :: seen(:)
  TYPE(iinterval) :: bound
  INTEGER(i8) :: sum_total, sum_selected, sum_not_selected
  CHARACTER(len=*), PARAMETER :: filename = 'test_combinatorics.f90'

!$omp parallel private(i, p, rand_seed, tid) &
!$omp shared(a, a_ref, bound, n, ppm_default_comm, &
!$omp sum_total, sum_selected, sum_not_selected)
#ifdef _OPENMP
  tid = omp_get_thread_num()
#else
  tid = 0
#endif
  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(2(a,i0))', 'thread id=', tid, ', random seed=', rand_seed
!$omp master
  n = irandr(iinterval(2, max_a_size))

  ALLOCATE(a(n), a_ref(n))
!$omp end master
!$omp barrier
  ! test permutation of random integers
  CALL a_rand_mt(a)
!$omp do
  DO i = 1, n
    a_ref(i) = a(i)
  END DO
!$omp end do
!$omp single
  CALL permute(a)
!$omp end single
!$omp sections
!$omp section
  CALL qsort_r_mt(a, SIZE(a), 4, 0, cmp_i4)
!$omp section
  CALL qsort_r_mt(a_ref, SIZE(a_ref), 4, 0, cmp_i4)
!$omp end sections
  p = .TRUE.
!$omp do
  DO i = 1, n
    p = p .AND. a(i) == a_ref(i)
  END DO
!$omp end do

  CALL assertion(p, filename, __LINE__, &
       'permutation does not contain same values as original')

  ! next permute sequence of integer range
!$omp single
  bound%first = MIN(irand(), HUGE(n) - n + 1)
  bound%last = bound%first + n - 1
!$omp end single
!$omp do
  DO i = 1, n
    a(i) = bound%first + i - 1
    a_ref(i) = a(i)
  END DO
!$omp end do
!$omp sections
!$omp section
  CALL assertion(is_permutation(a_ref, bound), filename, __LINE__, &
       'failure in permutation generation')
!$omp section
  CALL permute(a)
#ifdef HAVE_FC_OPENMP_TASK
!$omp task
#else
!$omp end sections
!$omp sections
!$omp section
#endif
  CALL assertion(is_permutation(a, bound), filename, __LINE__, &
       'failure in permutation check')
#ifdef HAVE_FC_OPENMP_TASK
!$omp end task
!$omp task
#else
!$omp section
#endif
  CALL assertion(is_permutation(a, a_ref), filename, __LINE__, &
       'failure in permutation check')
#ifdef HAVE_FC_OPENMP_TASK
!$omp end task
#endif
!$omp end sections
!$omp single
  CALL qsort_r_mt(a, SIZE(a), 4, 0, cmp_i4)
!$omp end single
  p = .TRUE.
!$omp do
  DO i = 1, n
    p = p .AND. a(i) == a_ref(i)
  END DO
!$omp end do
  CALL assertion(p, filename, __LINE__, &
       'permutation of integers does not sort into sequence of same integers')
!$omp master
  a_selection_size = irandr(iinterval(0, n))
  a_selection_size = MERGE(0, n, lrand())
  ALLOCATE(a_selection(a_selection_size), a_selection_i(n - a_selection_size))
  CALL combination(a_selection, a_selection_i, iinterval(1, n))
  ALLOCATE(seen(n))
  sum_total = 0_i8
  sum_selected = 0_i8
  sum_not_selected = 0_i8
!$omp end master
!$omp barrier
!$omp do
  DO i = 1, n
    seen(i) = .FALSE.
  END DO
!$omp end do
!$omp do
  DO i = 1, a_selection_size
    IF (a_selection(i) > n .OR. a_selection(i) < 1) &
         CALL abort_ppm("invalid element in selection!", filename, __LINE__)
    IF (seen(a_selection(i))) &
         CALL abort_ppm("duplicate element in selection!", filename, __LINE__)
    seen(a_selection(i)) = .TRUE.
  END DO
!$omp end do nowait
!$omp do
  DO i = 1, n - a_selection_size
    IF (a_selection_i(i) > n .OR. a_selection_i(i) < 1) &
         CALL abort_ppm("invalid element in selection!", filename, __LINE__)
    IF (seen(a_selection_i(i))) &
         CALL abort_ppm("duplicate element in selection!", filename, __LINE__)
    seen(a_selection_i(i)) = .TRUE.
  END DO
!$omp end do
#ifdef __INTEL_COMPILER
  IF (n > 0) THEN
#endif
!$omp do reduction(+: sum_total)
    DO i = 1, n
      sum_total = sum_total + INT(a(i), i8)
    END DO
!$omp end do
#ifdef __INTEL_COMPILER
  END IF
  IF (a_selection_size > 0) THEN
#endif
!$omp do reduction(+: sum_selected)
    DO i = 1, a_selection_size
      sum_selected = sum_selected + INT(a(a_selection(i)), i8)
    END DO
!$omp end do
#ifdef __INTEL_COMPILER
  END IF
  IF (n - a_selection_size > 0) THEN
#endif
!$omp do reduction(+: sum_not_selected)
    DO i = 1, n - a_selection_size
      sum_not_selected = sum_not_selected + INT(a(a_selection_i(i)), i8)
    END DO
!$omp end do
#ifdef __INTEL_COMPILER
  END IF
#endif
!$omp master
  CALL assertion(sum_total == sum_selected + sum_not_selected, &
       filename, __LINE__, "summation of combination error")
!$omp end master
!$omp end parallel
END PROGRAM test_combinatorics
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!
