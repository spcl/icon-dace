!>
!! @file test_heap.f90
!! @brief test heap data structure
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
PROGRAM test_heap
  USE ppm_std_type_kinds, ONLY: i4
  USE ppm_heap, ONLY: build_heap, heapify, is_heap, heap_elem_increase_sort
  USE ppm_compare, ONLY: cmp_i4
  USE ppm_base, ONLY: assertion, ppm_default_comm
  USE ppm_random, ONLY: initialize_irand, finalize_irand, a_rand, irandp, lrand
  IMPLICIT NONE
  INTEGER, PARAMETER :: max_size=100000
  INTEGER(i4), ALLOCATABLE :: a(:)
  INTEGER :: adjust, i, max_adjust, n, rand_seed
  LOGICAL :: l
  CHARACTER(len=*), PARAMETER :: filename = 'test_heap.f90'

  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(a,i0)', 'random seed=', rand_seed

  n = MOD(irandp(), max_size) + 1

  ALLOCATE(a(n))

  CALL a_rand(a)

  CALL build_heap(a, n, 4, 0, cmp_i4)

  CALL assertion(is_heap(a, n, 4, 0, cmp_i4), filename, __LINE__, &
       'sorted array is not heap?')

  i = MOD(irandp(), n) + 1
  l = lrand()

  IF (a(i) < HUGE(a(i)) .AND. l) THEN
    max_adjust = HUGE(a(i)) - MAX(a(i), 0)
    adjust = MOD(irandp(), max_adjust) + 1
    a(i) = a(i) + adjust
    CALL heap_elem_increase_sort(a, n, 4, i, 0, cmp_i4)
  ELSE
    max_adjust = HUGE(a(i)) - MAX(-a(i), 0)
    adjust = MOD(irandp(), max_adjust) + 1
    a(i) = a(i) - adjust
    CALL heapify(a, n, 4, i, 0, cmp_i4)
  END IF

  CALL assertion(is_heap(a, n, 4, 0, cmp_i4), filename, __LINE__, &
       'array lost heap property after adjustment')

  CALL finalize_irand

END PROGRAM test_heap
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!
