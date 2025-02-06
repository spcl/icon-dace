!>
!! @file test_qsort.f90
!! @brief test of qsort implementation in library
!!
!! @copyright Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: test quicksort
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
PROGRAM test_qsort
  USE ppm_random, ONLY: irand, initialize_irand, a_rand
  USE ppm_std_type_kinds, ONLY: i4, i8
  USE ppm_compare, ONLY: cmp_i4, cmp_i8
  USE ppm_base, ONLY: ppm_default_comm, abort => abort_ppm
  USE ppm_sort, ONLY: qsort_r, insertion_sort
  IMPLICIT NONE
  INTEGER, PARAMETER :: maxasize = 30000, minasize = 2
  INTEGER, PARAMETER :: maxbsize = 3000, minbsize = 2
  CHARACTER(len=132) :: msg
  INTEGER(i4), ALLOCATABLE :: a(:)
  INTEGER(i8), ALLOCATABLE :: b(:)
  INTEGER :: asize, bsize, i, rand_seed
  CHARACTER(len=14), PARAMETER :: filename = 'test_qsort.f90'

  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(a,i0)', 'random seed=', rand_seed

  asize = MODULO(irand(), maxasize - minasize + 1) + minasize
  ALLOCATE(a(asize))
  CALL a_rand(a)
  CALL qsort_r(a, SIZE(a), 4, 0, cmp_i4)
  DO i = 2, asize
    IF (a(i) < a(i - 1)) THEN
      WRITE (0, '(2(a,i0))') "discontinuity at index ", i - 1, " vs ", i
      WRITE (msg, '(4(a,i0))') "value a(", i - 1, ")=", a(i - 1), &
           "value a(", i, ")=", a(i)
      CALL abort(msg, filename, __LINE__)
    END IF
  END DO
  DEALLOCATE(a)
  bsize = MODULO(irand(), maxbsize - minbsize + 1) + minbsize
  ALLOCATE(b(bsize))
  CALL a_rand(b)
  CALL insertion_sort(b, SIZE(b), 8, 0, cmp_i8)
  DO i = 2, bsize
    IF (b(i) < b(i - 1)) THEN
      WRITE (0, '(2(a,i0))') "discontinuity at index ", i - 1, " vs ", i
      WRITE (msg, '(4(a,i0))') "value b(", i - 1, ")=", b(i - 1), &
           "value b(", i, ")=", b(i)
      CALL abort(msg, filename, __LINE__)
    END IF
  END DO
END PROGRAM test_qsort
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
