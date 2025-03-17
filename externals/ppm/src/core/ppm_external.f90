!! @file ppm_external.f90 --- gathers declarations for polymorphic C routines
!!                            and introduces appropriate renames
!
! Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!
! Version: 1.0
! Author: Thomas Jahns <jahns@dkrz.de>
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
!> do not use ppm_external directly, instead only the wrapper modules
!! ppm_sort and ppm_heap should be USE'd.

MODULE ppm_external
  IMPLICIT NONE
  PRIVATE
  ! generic quicksort implementation
  EXTERNAL :: ppm_qsort_r, ppm_qsort_r_mt
  PUBLIC :: ppm_qsort_r, ppm_qsort_r_mt
  ! generic binary heap data structure
  EXTERNAL :: ppm_build_heap, ppm_heapify, &
       ppm_heap_elem_increase_sort, ppm_heap_remove_top, &
       ppm_heap_leaf_minimize
  LOGICAL, EXTERNAL :: ppm_is_heap
  PUBLIC :: ppm_build_heap, ppm_heapify, ppm_is_heap, &
       ppm_heap_elem_increase_sort, ppm_heap_remove_top, &
       ppm_heap_leaf_minimize
  EXTERNAL :: ppm_insertion_sort, ppm_sorted_insertion, ppm_insertion_sort_once
  PUBLIC :: ppm_insertion_sort, ppm_sorted_insertion, ppm_insertion_sort_once
  INTEGER, EXTERNAL :: ppm_bsearch_r, ppm_bsearch_el_r
  PUBLIC :: ppm_bsearch_r, ppm_bsearch_el_r
END MODULE ppm_external

!> generic sort implementations
MODULE ppm_sort
  ! quick sort, good, all-purpose sort, but weak on specific sets
  USE ppm_external, ONLY: qsort_r => ppm_qsort_r, &
       qsort_r_mt => ppm_qsort_r_mt
  ! simpler sort routine which performs well on already sorted data
  USE ppm_external, ONLY: insertion_sort => ppm_insertion_sort, &
       sorted_insertion => ppm_sorted_insertion, &
       insertion_sort_once => ppm_insertion_sort_once
  USE ppm_std_type_kinds, ONLY: i4
  IMPLICIT NONE
  PRIVATE
  INTERFACE
    PURE SUBROUTINE ppm_qsort_i4_f(a, n, cmp)
      USE ppm_std_type_kinds, ONLY: i4
      INTEGER(i4), INTENT(in) :: n
      INTEGER(i4), INTENT(inout) :: a(n)
      INTERFACE
        PURE FUNCTION cmp(a, b) RESULT(t)
          USE ppm_std_type_kinds, ONLY: i4
          INTEGER(i4), INTENT(in) :: a, b
          INTEGER(i4) :: t
        END FUNCTION cmp
      END INTERFACE
    END SUBROUTINE ppm_qsort_i4_f
  END INTERFACE
  PUBLIC :: qsort_r, qsort_r_mt
  PUBLIC :: qsort_i4
  PUBLIC :: insertion_sort, sorted_insertion, insertion_sort_once
CONTAINS
  PURE SUBROUTINE qsort_i4(a, cmp)
    INTEGER(i4), INTENT(inout) :: a(:)
    INTERFACE
      PURE FUNCTION cmp(a, b) RESULT(t)
        IMPORT :: i4
        INTEGER(i4), INTENT(in) :: a, b
        INTEGER(i4) :: t
      END FUNCTION cmp
    END INTERFACE
    CALL ppm_qsort_i4_f(a, SIZE(a), cmp)
  END SUBROUTINE qsort_i4
END MODULE ppm_sort

!> generic binary heap data structure
MODULE ppm_heap
  USE ppm_external, ONLY: build_heap => ppm_build_heap, &
       heapify => ppm_heapify, is_heap => ppm_is_heap, &
       heap_elem_increase_sort => ppm_heap_elem_increase_sort, &
       heap_remove_top => ppm_heap_remove_top, &
       heap_leaf_minimize => ppm_heap_leaf_minimize
  IMPLICIT NONE
  PUBLIC
END MODULE ppm_heap

!> generic binary search routine
MODULE ppm_search
  USE ppm_external, ONLY: bsearch_r => ppm_bsearch_r, &
       bsearch_el_r => ppm_bsearch_el_r
  IMPLICIT NONE
  PUBLIC
END MODULE ppm_search
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
