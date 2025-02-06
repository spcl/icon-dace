!>
!! @file ppm_general_block_decomposition.f90
!! @brief build general BLOCK distribution
!!
!! @copyright Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: ScalES PPM general block distribution
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
MODULE ppm_general_block_decomposition
  USE ppm_extents, ONLY: extent, extent_start, extent_end, &
       extent_size
  USE ppm_base, ONLY: assertion
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: gbd_refine_2d

  CHARACTER(len=*), PARAMETER :: &
       filename = 'ppm_general_block_decomposition.f90'
CONTAINS
  !> serial implementation of iterative refinement of a partitioning
  !! by alternating optimization of horizontal and vertical divisions
  !! yields a general block decomposition
  SUBROUTINE gbd_refine_2d(global_range, weights, partitions_x, partitions_y)
    TYPE(extent), INTENT(in) :: global_range(2)
    REAL, INTENT(in) :: weights(global_range(1)%first:, &
         global_range(2)%first:)
    TYPE(extent), INTENT(inout) :: partitions_x(:), &
         partitions_y(:)

    INTEGER :: size_x, size_y, nparts_x, nparts_y
    REAL :: max_weight
!    REAL :: part_weight(SIZE(partitions_x),SIZE(partitions_y))

    nparts_x = SIZE(partitions_x)
    nparts_y = SIZE(partitions_y)
    size_x = SIZE(weights, 1)
    size_y = SIZE(weights, 2)
    CALL assertion(size_x == extent_size(global_range(1)) &
         .AND. size_y == extent_size(global_range(2)), line=__LINE__, &
         source=filename)

    CALL compute_weights

  CONTAINS
    SUBROUTINE compute_weights
      INTEGER :: i, j
      REAL :: weight
      weight = -HUGE(weight)
      DO j = 1, nparts_y
        DO i = 1, nparts_x
          weight = SUM(weights(extent_start(partitions_x(i)): &
               extent_end(partitions_x(i)), extent_start(partitions_y(j)): &
               extent_end(partitions_y(j))))
!          part_weight(i, j) = weight
          IF (weight > max_weight) max_weight = weight
        END DO
      END DO
    END SUBROUTINE compute_weights
  END SUBROUTINE gbd_refine_2d

END MODULE ppm_general_block_decomposition
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
