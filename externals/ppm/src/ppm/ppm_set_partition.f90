!>
!! @file ppm_set_partition.f90
!! @brief compute set partitions
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
MODULE ppm_set_partition
  USE ppm_std_type_kinds, ONLY: i4, i8, sizeof_integer
  USE ppm_base, ONLY: assertion
  USE ppm_compare, ONLY: cmp_i4, cmp_i4_indirect_i8
  USE ppm_set_partition_base, ONLY : set_i4, partition_vec
  USE ppm_sort, ONLY: qsort_r
  IMPLICIT NONE
  PRIVATE
  INTERFACE greedy_partitioning
    MODULE PROCEDURE greedy_partitioning_set_i4
    MODULE PROCEDURE greedy_partitioning_pv
  END INTERFACE greedy_partitioning
  PUBLIC :: greedy_partitioning

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_set_partition.f90'
CONTAINS
  !> partition a set, simple heuristics
  SUBROUTINE greedy_partitioning_set_i4(partitioning, weight)
    TYPE(set_i4), INTENT(out) :: partitioning(:)
    INTEGER(i4), INTENT(in) :: weight(:)

    INTEGER :: i, j, m, n, nparts
    INTEGER(i4), ALLOCATABLE :: part_assign(:,:)
    INTEGER :: part_size(SIZE(partitioning)), part_heap(SIZE(partitioning))
    INTEGER(i8) :: part_weight(SIZE(partitioning))
    !> low water partition, i.e. a partition with minimal weight-sum
    INTEGER :: lwp

    nparts = SIZE(partitioning)
    CALL assertion(nparts > 0, line=__LINE__, source=filename, &
         msg="number of partitions must be positive integer")

    n = SIZE(weight)

    ALLOCATE(part_assign(2, n))

    DO j = 1, n
      part_assign(1, j) = weight(j)
      part_assign(2, j) = j
    END DO
    CALL qsort_r(part_assign, n, 8, 0, cmp_i4)

    m = MIN(nparts, n)
    DO i = 1, m
      !> at this point: part_assign(1, i) == weight(part_assign(2, i))
      part_size(i) = 1
      part_weight(i) = INT(part_assign(1, i), i8)
      part_heap(nparts-i+1) = i
      part_assign(1, i) = i
    END DO
    part_size(m+1:nparts) = 0
    part_weight(m+1:nparts) = 0_i8
    DO j = 1, nparts - m
      part_heap(j) = j+m
    END DO
    DO j = nparts+1, n
      lwp = part_heap(1)
      part_weight(lwp) = part_weight(lwp) + INT(part_assign(1, j), i8)
      part_size(lwp) = part_size(lwp) + 1
      part_assign(1, j) = lwp
      !> workaround for unavailable heap routines
      CALL qsort_r(part_heap, nparts, sizeof_integer, part_weight, &
           cmp_i4_indirect_i8)
    END DO
    DO i = 1, nparts
      ALLOCATE(partitioning(i)%elem(part_size(i)))
    END DO
    !> reuse part_size as part_fill
    part_size = 0
    DO j = 1, n
      i = part_assign(1, j)
      m = part_assign(2, j)
      part_size(i) = part_size(i) + 1
      partitioning(i)%elem(part_size(i)) = m
    END DO
    DEALLOCATE(part_assign)
  END SUBROUTINE greedy_partitioning_set_i4

  SUBROUTINE greedy_partitioning_pv(partitioning, weight)
    TYPE(partition_vec), INTENT(inout) :: partitioning
    INTEGER(i4), INTENT(in) :: weight(:)

    INTEGER :: i, j, m, n, nparts
    INTEGER, ALLOCATABLE :: part_assign(:,:)
    INTEGER :: part_size(SIZE(partitioning%start) - 1), &
         part_heap(SIZE(partitioning%start) - 1)
    INTEGER(i8) :: part_weight(SIZE(partitioning%start) - 1)
    !> low water partition, i.e. a partition with minimal weight-sum
    INTEGER :: lwp

    nparts = SIZE(partitioning%start) - 1
    CALL assertion(nparts > 0, line=__LINE__, source=filename, &
         msg="number of partitions must be positive integer")

    n = SIZE(weight)
    CALL assertion(n == SIZE(partitioning%elements), line=__LINE__, &
         source=filename, &
         msg="partitioning must be allocated to correct size")

    ALLOCATE(part_assign(2, n))

    DO j = 1, n
      part_assign(1, j) = weight(j)
      part_assign(2, j) = j
    END DO
    CALL qsort_r(part_assign, n, sizeof_integer * 2, 0, cmp_i4)

    m = MIN(nparts, n)
    DO i = 1, m
      !> at this point: part_assign(1, i) == weight(part_assign(2, i))
      part_size(i) = 1
      part_weight(i) = INT(part_assign(1, i), i8)
      part_heap(nparts-i+1) = i
      part_assign(1, i) = i
    END DO
    part_size(m+1:nparts) = 0
    part_weight(m+1:nparts) = 0_i8
    DO j = 1, nparts - m
      part_heap(j) = j+m
    END DO
    DO j = nparts+1, n
      lwp = part_heap(1)
      part_weight(lwp) = part_weight(lwp) + INT(part_assign(1, j), i8)
      part_size(lwp) = part_size(lwp) + 1
      part_assign(1, j) = lwp
      !> workaround for unavailable heap routines
      CALL qsort_r(part_heap, nparts, sizeof_integer, part_weight, &
           cmp_i4_indirect_i8)
    END DO
    partitioning%start(1) = 1
    DO i = 1, nparts
      partitioning%start(i + 1) = partitioning%start(i) + part_size(i)
    END DO
    !> reuse part_size as part_fill
    part_size = 0
    DO j = 1, n
      i = part_assign(1, j)
      m = part_assign(2, j)
      partitioning%elements(partitioning%start(i) + part_size(i)) = m
      part_size(i) = part_size(i) + 1
    END DO
    DEALLOCATE(part_assign)
  END SUBROUTINE greedy_partitioning_pv

END MODULE ppm_set_partition
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
