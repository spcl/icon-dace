!> @file ppm_uniform_partition.f90
!! @brief compute uniform partitions
!!
!! @copyright Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
! Keywords: uniform distribution partition
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
!
!> compute uniform partitioning of n-dimensional rectilinear
#include "fc_feature_defs.inc"
MODULE ppm_uniform_partition
  USE ppm_base, ONLY: assertion
  USE ppm_extents, ONLY: extent, extent_size, extent_start
  USE ppm_std_type_kinds, ONLY: i8
  USE ppm_set_partition_base, ONLY: block_decomposition
  IMPLICIT NONE
  PRIVATE

  INTERFACE uniform_decomposition
    MODULE PROCEDURE uniform_decomposition_1d
    MODULE PROCEDURE uniform_decomposition_nd
  END INTERFACE uniform_decomposition

  INTERFACE partidx_of_elem
    MODULE PROCEDURE partidx_of_elem_uniform_deco
  END INTERFACE partidx_of_elem

  PUBLIC :: uniform_partition, uniform_decomposition
  PUBLIC :: partidx_of_elem, uniform_partition_start

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_uniform_partition.f90'
CONTAINS
  !> compute start integer of uniform interval partition
  ELEMENTAL FUNCTION uniform_partition_start(set_interval, nparts, &
       part_idx, symmetric) RESULT(start)
    INTEGER, INTENT(in) :: nparts
    TYPE(extent), INTENT(in) :: set_interval
    INTEGER, INTENT(in) :: part_idx
    LOGICAL, INTENT(in) :: symmetric

    INTEGER :: start, part_offset, sym_part_idx
    INTEGER(i8) :: sym_size

    IF (symmetric) THEN
      sym_size = (INT(extent_size(set_interval), i8) * INT(nparts/2, i8)) &
           / INT(nparts, i8)
      sym_part_idx = MERGE(part_idx - 1, nparts - part_idx + 1, &
           part_idx <= (nparts+1)/2)
      part_offset = INT((sym_size * INT(sym_part_idx, i8)) / INT(nparts/2, i8))
      IF (part_idx > (nparts + 1)/2) THEN
        part_offset = extent_size(set_interval) - part_offset
      END IF
    ELSE
      part_offset = INT((INT(extent_size(set_interval), i8) &
           &             * INT(part_idx - 1, i8)) / INT(nparts, i8))
    END IF
    start = extent_start(set_interval) + part_offset
  END FUNCTION uniform_partition_start

  !> compute nth part of integer set interval
  !!
  !! The interval is divided into roughly same sized sub-intervals
  !! forming a uniform partition.
  !! @param set_interval global domain
  !! @param nparts number of parts to decompose into
  !! @param part_idx number of sub-interval to compute
  !! @param symmetric if .true. ensure that
  !! <tt>SIZE(uniform_partition(i)) == SIZE(uniform_partition(nparts - i + 1))</tt>
  !! @return part range corresponding to part_idx
  ELEMENTAL FUNCTION uniform_partition(set_interval, nparts, part_idx, &
       symmetric) RESULT(interval)
    INTEGER, INTENT(in) :: nparts
    TYPE(extent), INTENT(in) :: set_interval
    INTEGER, INTENT(in) :: part_idx
    LOGICAL, OPTIONAL, INTENT(in) :: symmetric

    TYPE(extent) :: interval
    LOGICAL :: symmetric_pass
    INTEGER :: start_part, start_next_part
    IF (PRESENT(symmetric)) THEN
      symmetric_pass = symmetric
    ELSE
      symmetric_pass = .FALSE.
    END IF
    start_part = uniform_partition_start(set_interval, nparts, part_idx, &
         symmetric_pass)
    start_next_part = uniform_partition_start(set_interval, nparts, &
         part_idx + 1, symmetric_pass)
    interval = extent(start_part, start_next_part - start_part)
  END FUNCTION uniform_partition

  !> divide integer set interval into evenly sized sub-intervals
  !! forming a partition
  !!
  !! @param set_interval interval to partition
  !! @param nparts number of parts to compute
  !! @param parts \a nparts sub-intervals
  !! @param symmetric if .true. partitions(i) will be same-size as
  !! partitions(nparts-i+1)
  !!
  SUBROUTINE uniform_decomposition_1d(set_interval, nparts, &
       parts, symmetric)
    INTEGER, INTENT(in) :: nparts
    TYPE(extent), INTENT(in) :: set_interval
    TYPE(extent), INTENT(out) :: parts(nparts)
    LOGICAL, OPTIONAL, INTENT(in) :: symmetric

    INTEGER :: i

    DO i = 1, nparts
      parts(i) = uniform_partition(set_interval, nparts, i, symmetric=symmetric)
    END DO

  END SUBROUTINE uniform_decomposition_1d


  !> Compute uniform divisions of rectilinear structure.
  !!
  !! @param pgrid same dimension as rect and nparts, holds list of
  !! partition extents
  !! @param rect intervals for indices of rectilinear to partition
  !! @param nparts number of parts to compute in each dimension
  !! @param symmetric symmetry requirements for each dimension
  SUBROUTINE uniform_decomposition_nd(pgrid, rect, nparts, symmetric)
    TYPE(extent), INTENT(in) :: rect(:)
    INTEGER, INTENT(in) :: nparts(:)
    TYPE(block_decomposition), INTENT(out) :: pgrid(:)
    LOGICAL, OPTIONAL, INTENT(in) :: symmetric(:)

    INTEGER :: i, n
    LOGICAL :: symmetry(SIZE(rect))

    n = SIZE(rect)
    CALL assertion(SIZE(nparts) == n, filename, __LINE__, &
         msg='argument parts must be same size as argument rect')
    CALL assertion(SIZE(pgrid) == n, filename, __LINE__, &
         msg='argument pgrids must be same size as argument rect')
    symmetry = .FALSE.
    IF (PRESENT(symmetric)) THEN
      CALL assertion(SIZE(symmetric) == n, filename, __LINE__, &
           msg='argument symmetry must be same size as argument rect')
      symmetry = symmetric
    END IF
    DO i = 1, n
      ALLOCATE(pgrid(i)%partition(nparts(i)))
      CALL uniform_decomposition_1d(rect(i), nparts(i), &
           parts=pgrid(i)%partition, symmetric=symmetry(i))
    END DO
  END SUBROUTINE uniform_decomposition_nd

  ELEMENTAL FUNCTION partidx_of_elem_uniform_deco(set_interval, nparts, &
       elem_idx) RESULT(part_idx)
    TYPE(extent), INTENT(in) :: set_interval
    INTEGER, INTENT(in) :: nparts, elem_idx
    INTEGER :: part_idx
    part_idx = INT((INT(elem_idx - extent_start(set_interval), i8) &
         * INT(nparts, i8) + INT(nparts, i8) - 1_i8) &
         / INT(extent_size(set_interval), i8)) + 1
  END FUNCTION partidx_of_elem_uniform_deco

END MODULE ppm_uniform_partition
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
