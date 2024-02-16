!>
!! @file test_bsearch.f90
!! @brief test binary search functions
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
PROGRAM test_bsearch
  USE ppm_base, ONLY: assertion, ppm_default_comm
  USE ppm_compare, ONLY: cmp_i4
  USE ppm_extents, ONLY: iinterval
  USE ppm_random, ONLY: initialize_irand, irandr, a_rand
  USE ppm_search, ONLY: bsearch_r, bsearch_el_r
  USE ppm_sort, ONLY: qsort_r
  USE ppm_std_type_kinds, ONLY: i4
  IMPLICIT NONE
  INTEGER, PARAMETER :: maxasize = 30000, minasize = 2
  TYPE(iinterval), PARAMETER :: asize_range = iinterval(minasize, maxasize)
  INTEGER(i4), ALLOCATABLE :: a(:)
  INTEGER :: asize, i, rand_seed
  CHARACTER(len=*), PARAMETER :: filename = 'test_bsearch.f90'

  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(a,i0)', 'random seed=', rand_seed

  asize = irandr(asize_range)
  ALLOCATE(a(asize))

  DO i=1, asize
    a(i) = i * 10
  END DO

  i = irandr(iinterval(1, asize))
  CALL assertion(bsearch_r(i * 10, a, asize, 4, 0, cmp_i4) == i, &
       filename, __LINE__, 'expected value not found in array')
  CALL assertion(bsearch_el_r(i * 10, a, asize, 4, 0, cmp_i4) == i, &
       filename, __LINE__, 'expected value not found in array')
  a(i) = a(i) + 5
  CALL assertion(bsearch_r(i * 10, a, asize, 4, 0, cmp_i4) == 0, &
       filename, __LINE__, 'unexpected value found in array')
  CALL assertion(bsearch_el_r(i * 10, a, asize, 4, 0, cmp_i4) == i - 1, &
       filename, __LINE__, 'expected element not found in array')
  CALL a_rand(a)
  CALL qsort_r(a, asize, 4, 0, cmp_i4)
  i = irandr(iinterval(1, asize))
  CALL assertion(a(bsearch_r(a(i), a, asize, 4, 0, cmp_i4)) == a(i), &
       filename, __LINE__, 'expected value not found in array')
  CALL assertion(a(bsearch_el_r(a(i), a, asize, 4, 0, cmp_i4)) == a(i), &
       filename, __LINE__, 'expected value not found in array')
END PROGRAM test_bsearch
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
