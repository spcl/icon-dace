!>
!! @file test_set_partition_base.f90
!! @brief test equivalence of basic partition data structures
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
PROGRAM test_set_partition_base
  USE ppm_base, ONLY: assertion, ppm_default_comm
  USE ppm_random, ONLY: initialize_irand, a_randr, irandp, irandr
  USE ppm_extents, ONLY: iinterval, ASSIGNMENT(=)
  USE ppm_combinatorics, ONLY: combination
  USE ppm_std_type_kinds, ONLY: i4
  USE ppm_set_partition_base, ONLY: OPERATOR(==), ASSIGNMENT(=), &
       partition_assignment, set_i4, partition_vec, assign_set_i4_2_pv, &
       OPERATOR(/=)
  IMPLICIT NONE
  TYPE(partition_assignment) :: parts_a, parts_d
  TYPE(set_i4), ALLOCATABLE :: parts_b(:)
  TYPE(partition_vec) :: parts_c
  INTEGER, PARAMETER :: max_num_parts=1024, max_num_elements=1024**2
!  INTEGER, PARAMETER :: max_num_parts=5, max_num_elements=10
  INTEGER :: num_parts, num_elements, rand_seed
  LOGICAL :: test_changed
  CHARACTER(len=*), PARAMETER :: filename = 'test_set_partition_base.f90'

  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(a,i0)', 'random seed=', rand_seed

  num_elements = MOD(irandp(), max_num_elements)
  num_parts = MOD(irandp(), max_num_parts) + 1

  ALLOCATE(parts_a%assigned(num_elements))
  parts_a%part_range = iinterval(1, num_parts)
  CALL a_randr(parts_a%assigned, iinterval(1, num_parts))

  ALLOCATE(parts_b(num_parts))
  parts_b = parts_a

  CALL assertion(parts_a == parts_b, filename, __LINE__, &
       "representation inequal after assignment")

  CALL assign_set_i4_2_pv(parts_c, parts_b)

  CALL assertion(parts_a == parts_c .AND. parts_c == parts_b, &
       filename, __LINE__, "representation inequal after assignment")

  parts_d = parts_c

  CALL assertion(parts_d == parts_c .AND. parts_d == parts_b &
       .AND. parts_d == parts_a, filename, &
       __LINE__, "representation inequal after assignment")

  CALL change_part(parts_b, test_changed)

  IF (test_changed) THEN
    CALL assertion(parts_b /= parts_a .AND. parts_b /= parts_c &
         .AND. parts_b /= parts_d, filename, &
         __LINE__, "representation equal after modification")
  ELSE
    CALL assertion(parts_b == parts_a .AND. parts_b == parts_c &
         .AND. parts_b == parts_d, filename, &
         __LINE__, "representation equal after modification")
  END IF

CONTAINS
  SUBROUTINE change_part(sets, changed)
    TYPE(set_i4), ALLOCATABLE, INTENT(inout) :: sets(:)
    LOGICAL, INTENT(out) :: changed
    INTEGER(i4), ALLOCATABLE :: buf(:), selected(:), not_selected(:)
    INTEGER :: i, j, m, n, src, dest, src_size, dest_size
    INTEGER :: nchanges

    IF (SIZE(sets) < 2) THEN
      changed = .FALSE.
      RETURN
    END IF
    i = irandr(iinterval(1, num_parts))
    j = irandr(iinterval(1, num_parts-1))
    j = j + MERGE(1, 0, j >= i)
    m = SIZE(sets(i)%elem)
    n = SIZE(sets(j)%elem)

    nchanges = irandr(iinterval(-m, n))

    IF (nchanges < 0) THEN
      ! move -nchanges elements from partition i to partition j
      src = i
      dest = j
      nchanges = -nchanges
    ELSE IF (nchanges > 0) THEN
      ! move nchanges elements from partition j to partition i
      src = j
      dest = i
    ELSE ! nchanges == 0
      changed = .FALSE.
      RETURN
    END IF
    src_size = SIZE(sets(src)%elem)
    dest_size = SIZE(sets(dest)%elem)

    ALLOCATE(buf(dest_size))

    buf = sets(dest)%elem
    DEALLOCATE(sets(dest)%elem)
    ALLOCATE(sets(dest)%elem(dest_size + nchanges))
    sets(dest)%elem(1:dest_size) = buf

    ALLOCATE(selected(nchanges), not_selected(src_size-nchanges))
    CALL combination(selected, not_selected, iinterval(1, src_size))
    DO i = 1, nchanges
      sets(dest)%elem(dest_size+i) = sets(src)%elem(selected(i))
    END DO
    DEALLOCATE(buf)
    ALLOCATE(buf(src_size - nchanges))
    buf = sets(src)%elem(not_selected)
    DEALLOCATE(sets(src)%elem)
    ALLOCATE(sets(src)%elem(src_size - nchanges))
    sets(src)%elem = buf
    DEALLOCATE(buf)
    changed = .TRUE.

  END SUBROUTINE change_part
END PROGRAM test_set_partition_base
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!
