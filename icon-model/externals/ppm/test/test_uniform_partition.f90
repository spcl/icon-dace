!>
!! @file test_uniform_partition.f90
!! @brief test uniform partition generator functions
!!
!! @copyright Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: test uniform partitioning
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
PROGRAM test_uniform_partition
  USE ppm_base, ONLY: assertion, ppm_default_comm
  USE ppm_extents, ONLY: extent_size, extent_start, extent_end, extent, &
       OPERATOR(==)
  USE ppm_random, ONLY: irandp, lrand, initialize_irand
  USE ppm_set_partition_base, ONLY: block_decomposition
  USE ppm_uniform_partition, ONLY: uniform_partition, uniform_decomposition, &
       partidx_of_elem
  IMPLICIT NONE
  INTEGER, PARAMETER :: maxsize = 1000000, maxparts=maxsize
  TYPE(extent), ALLOCATABLE :: domain(:)
  INTEGER, ALLOCATABLE :: nparts(:)
  LOGICAL, ALLOCATABLE :: symmetry(:)
  TYPE(extent), ALLOCATABLE :: refpart(:)
  TYPE(block_decomposition), ALLOCATABLE :: test_nd(:)
  INTEGER :: i, n, sbase, ssize, rand_seed
  CHARACTER(len=*), PARAMETER :: filename = 'test_uniform_partition.f90'

  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(a,i0)', 'random seed=', rand_seed

  ! generate number of dimensions to partition in
  n = MOD(irandp(), 9) + 2

  ALLOCATE(domain(n), nparts(n), symmetry(n), test_nd(n))
  DO i = 1, n
    sbase = MOD(irandp(), maxsize)
    ssize = MOD(irandp(), MIN(maxparts, HUGE(i)-sbase))
    nparts(i) = MOD(irandp(), ssize) + 1
    symmetry(i) = lrand()
    IF (symmetry(i)) THEN
      ! can't divide odd number of items into even number of evenly sized parts
      IF (MOD(nparts(i), 2) /= MOD(ssize, 2)) ssize = ssize + 1
    END IF
    domain(i) = extent(sbase, ssize)
  END DO
  CALL uniform_decomposition(test_nd, domain, nparts, symmetry)
  DO i = 1, n
    IF (ALLOCATED(refpart)) DEALLOCATE(refpart)
    ALLOCATE(refpart(nparts(i)))
    CALL test_1d_part(refpart, domain(i), nparts(i), symmetry(i))
    CALL assertion(ALL(refpart == test_nd(i)%partition), filename, __LINE__, &
         msg='inequal partitioning results')
  END DO

CONTAINS

  SUBROUTINE test_1d_part(parts1d, part_range, nparts, symmetric)
    TYPE(extent), INTENT(out) :: parts1d(:)
    TYPE(extent), INTENT(in) :: part_range
    INTEGER, INTENT(in) :: nparts
    LOGICAL, INTENT(in) :: symmetric

    INTEGER :: i, imax, j, je, part_idx
    imax = MERGE(nparts, (nparts + 1)/2, .NOT. symmetric)
    parts1d(1) = uniform_partition(part_range, &
         nparts, 1, symmetric)
    CALL assertion(extent_start(parts1d(1)) == extent_start(part_range), &
         line=__LINE__, source=filename)
    IF (symmetric) THEN
      parts1d(nparts) = uniform_partition(part_range, nparts, nparts, symmetric)
      CALL assertion(extent_size(parts1d(1)) &
           &         == extent_size(parts1d(nparts)), filename, &
           __LINE__, "symmetric size equality assertion")
    END IF
    DO i = 2, imax
      parts1d(i) = uniform_partition(part_range, &
           nparts, i, symmetric)
      CALL assertion(extent_end(parts1d(i - 1)) &
           == extent_start(parts1d(i)) - 1, filename, __LINE__)
      CALL assertion(ABS(extent_size(parts1d(i - 1)) &
           - extent_size(parts1d(i))) &
           < MERGE(3, 2, symmetric), filename, __LINE__, &
           msg=TRIM(MERGE("symmetric    ","non-symmetric",symmetric)) &
           &   // " partitions size difference too big")
      IF (symmetric) THEN
        parts1d(nparts - i + 1) = &
             uniform_partition(part_range, nparts, nparts - i + 1, symmetric)
        CALL assertion(extent_size(parts1d(i)) &
             == extent_size(parts1d(nparts - i + 1)), filename, &
             __LINE__, "symmetric size equality assertion")
      END IF
    END DO
    IF (.NOT. symmetric) THEN
      DO i = 1, imax
        je = extent_end(parts1d(i))
        DO j = extent_start(parts1d(i)), je
          part_idx = partidx_of_elem(part_range, nparts, j)
          IF (.NOT. part_idx == i) THEN
            CALL assertion(part_idx == i, filename, __LINE__, &
                 "part assignment mismatch!")
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE test_1d_part
END PROGRAM test_uniform_partition
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
