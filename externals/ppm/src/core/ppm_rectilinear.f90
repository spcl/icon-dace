!>
!! @file ppm_rectilinear.f90
!! @brief utility routines for handling rectilinear data
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
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
!> Utility routines to use array of extents as rectilinear descriptor
#include "fc_feature_defs.inc"
MODULE ppm_rectilinear
  USE ppm_base, ONLY: assertion, abort_ppm
  USE ppm_extents, ONLY: extent, extent_end, extent_size, extent_start, &
       iinterval
  USE ppm_std_type_kinds, ONLY: i4
  IMPLICIT NONE
  PRIVATE
  INTERFACE rlcoord2lidx
    MODULE PROCEDURE rlcoord2lidx_e
    MODULE PROCEDURE rlcoord2lidx_i
  END INTERFACE rlcoord2lidx
  INTERFACE lidx2rlcoord
    MODULE PROCEDURE lidx2rlcoord_e
    MODULE PROCEDURE lidx2rlcoord_i
  END INTERFACE lidx2rlcoord
  INTERFACE num_neighbours_of_rect_elem
    MODULE PROCEDURE num_neighbours_of_rect_elem_e
    MODULE PROCEDURE num_neighbours_of_rect_elem_i
  END INTERFACE num_neighbours_of_rect_elem
  INTERFACE lidx_nb_coords
    MODULE PROCEDURE lidx_nb_coords_e
    MODULE PROCEDURE lidx_nb_coords_i
  END INTERFACE lidx_nb_coords
  INTERFACE lidx_nb_indices
    MODULE PROCEDURE lidx_nb_indices_e
    MODULE PROCEDURE lidx_nb_indices_i
  END INTERFACE lidx_nb_indices
  PUBLIC :: rlcoord2lidx, lidx2rlcoord
  PUBLIC :: num_neighbours_of_rect_elem
  PUBLIC :: lidx_nb_coords, lidx_nb_indices

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_rectilinear.f90'
CONTAINS
  FUNCTION rlcoord2lidx_e(rect, coord) RESULT(idx)
    TYPE(extent), INTENT(in) :: rect(:)
    INTEGER(i4), INTENT(in) :: coord(:)
    INTEGER(i4) :: idx

    INTEGER :: i, n

    n = SIZE(rect)
    CALL assertion(n > 0 .AND. SIZE(coord) == n, filename, __LINE__, &
         'rect and coordinate must have identical, non-zero dimensions')

    idx = 1
    IF (n > 0) THEN
      idx = idx + coord(1) - extent_start(rect(1))
    END IF
    DO i = 2, n
      idx = idx + extent_size(rect(1:i-1)) * (coord(i) - extent_start(rect(i)))
    END DO
  END FUNCTION rlcoord2lidx_e

  FUNCTION rlcoord2lidx_i(rect, coord) RESULT(idx)
    TYPE(iinterval), INTENT(in) :: rect(:)
    INTEGER(i4), INTENT(in) :: coord(:)
    INTEGER(i4) :: idx

    INTEGER :: i, n

    n = SIZE(rect)
    CALL assertion(n > 0 .AND. SIZE(coord) == n, filename, __LINE__, &
         'rect and coordinate must have identical, non-zero dimensions')

    idx = 1
    IF (n > 0) THEN
      idx = idx + coord(1) - extent_start(rect(1))
    END IF
    DO i = 2, n
      idx = idx + extent_size(rect(1:i-1)) * (coord(i) - extent_start(rect(i)))
    END DO
  END FUNCTION rlcoord2lidx_i

  FUNCTION lidx2rlcoord_e(rect, idx) RESULT(coord)
    TYPE(extent), INTENT(in) :: rect(:)
    INTEGER(i4), INTENT(in) :: idx
    INTEGER(i4) :: coord(SIZE(rect))

    INTEGER :: base, i, ofs, rest

    rest = idx - 1
    DO i = SIZE(coord), 2, -1
      base = extent_size(rect(1:i-1))
      ofs = rest / base
      coord(i) =  ofs + extent_start(rect(i))
      rest = rest - ofs * base
    END DO
    coord(1) = rest + extent_start(rect(1))
  END FUNCTION lidx2rlcoord_e

  FUNCTION lidx2rlcoord_i(rect, idx) RESULT(coord)
    TYPE(iinterval), INTENT(in) :: rect(:)
    INTEGER(i4), INTENT(in) :: idx
    INTEGER(i4) :: coord(SIZE(rect))

    INTEGER :: base, i, ofs, rest

    rest = idx - 1
    DO i = SIZE(coord), 2, -1
      base = extent_size(rect(1:i-1))
      ofs = rest / base
      coord(i) =  ofs + extent_start(rect(i))
      rest = rest - ofs * base
    END DO
    coord(1) = rest + extent_start(rect(1))
  END FUNCTION lidx2rlcoord_i

  FUNCTION num_neighbours_of_rect_elem_e(rect, coord) RESULT(nnb)
    TYPE(extent), INTENT(in) :: rect(:)
    INTEGER(i4), INTENT(in) :: coord(:)
    INTEGER :: nnb

    INTEGER :: i, n, s, e
    LOGICAL :: is_in_bounds

    n = SIZE(rect)
    CALL assertion(SIZE(coord) == n, filename, __LINE__, &
         'rect and coordinate must have identical, non-zero dimensions')
    nnb = 0
    is_in_bounds = .TRUE.
    DO i = 1, n
      s = extent_start(rect(i))
      e = extent_end(rect(i))
      is_in_bounds = is_in_bounds .AND. coord(i) <= e .AND. coord(i) >= s
      nnb = nnb + MERGE(1, 0, coord(i) > s) &
           &    + MERGE(1, 0, coord(i) < e)
    END DO
    IF (.NOT. is_in_bounds) &
         CALL abort_ppm('coordinate must be contained in rect', filename, &
         __LINE__)
  END FUNCTION num_neighbours_of_rect_elem_e

  FUNCTION num_neighbours_of_rect_elem_i(rect, coord) RESULT(nnb)
    TYPE(iinterval), INTENT(in) :: rect(:)
    INTEGER(i4), INTENT(in) :: coord(:)
    INTEGER :: nnb

    INTEGER :: i, n, s, e
    LOGICAL :: is_in_bounds

    n = SIZE(rect)
    CALL assertion(SIZE(coord) == n, filename, __LINE__, &
         'rect and coordinate must have identical, non-zero dimensions')
    nnb = 0
    is_in_bounds = .TRUE.
    DO i = 1, n
      s = extent_start(rect(i))
      e = extent_end(rect(i))
      is_in_bounds = is_in_bounds .AND. coord(i) <= e .AND. coord(i) >= s
      nnb = nnb + MERGE(1, 0, coord(i) > s) &
           &    + MERGE(1, 0, coord(i) < e)
    END DO
    IF (.NOT. is_in_bounds) &
         CALL abort_ppm('coordinate must be contained in rect', filename, &
         __LINE__)
  END FUNCTION num_neighbours_of_rect_elem_i

  SUBROUTINE lidx_nb_coords_e(rect, idx, nbcoord)
    TYPE(extent), INTENT(in) :: rect(:)
    INTEGER(i4), INTENT(in) :: idx
    INTEGER(i4), INTENT(out) :: nbcoord(:, :)

    INTEGER :: i, k, n
    INTEGER(i4) :: coord(SIZE(rect))

    n = SIZE(rect)
    CALL assertion(SIZE(nbcoord, 1) == n, filename, __LINE__, &
         'rect and neighbour coordinate(s) must have identical dimensions')
    coord = lidx2rlcoord(rect, idx)
    k = 1
    DO i = 1, n
      IF (coord(i) > extent_start(rect(i))) THEN
        nbcoord(:, k) = coord
        nbcoord(i, k) = coord(i) - 1
        k = k + 1
      END IF
      IF (coord(i) < extent_end(rect(i))) THEN
        nbcoord(:, k) = coord
        nbcoord(i, k) = coord(i) + 1
        k = k + 1
      END IF
    END DO
  END SUBROUTINE lidx_nb_coords_e

  SUBROUTINE lidx_nb_coords_i(rect, idx, nbcoord)
    TYPE(iinterval), INTENT(in) :: rect(:)
    INTEGER(i4), INTENT(in) :: idx
    INTEGER(i4), INTENT(out) :: nbcoord(:, :)

    INTEGER :: i, k, n
    INTEGER(i4) :: coord(SIZE(rect))

    n = SIZE(rect)
    CALL assertion(SIZE(nbcoord, 1) == n, filename, __LINE__, &
         'rect and neighbour coordinate(s) must have identical dimensions')
    coord = lidx2rlcoord(rect, idx)
    k = 1
    DO i = 1, n
      IF (coord(i) > extent_start(rect(i))) THEN
        nbcoord(:, k) = coord
        nbcoord(i, k) = coord(i) - 1
        k = k + 1
      END IF
      IF (coord(i) < extent_end(rect(i))) THEN
        nbcoord(:, k) = coord
        nbcoord(i, k) = coord(i) + 1
        k = k + 1
      END IF
    END DO
  END SUBROUTINE lidx_nb_coords_i

  SUBROUTINE lidx_nb_indices_e(rect, idx, nbidx)
    TYPE(extent), INTENT(in) :: rect(:)
    INTEGER(i4), INTENT(in) :: idx
    INTEGER(i4), INTENT(out) :: nbidx(:)

    INTEGER :: i, k, n
    INTEGER(i4) :: coord(SIZE(rect))

    n = SIZE(rect)
    coord = lidx2rlcoord(rect, idx)
    k = 1
    DO i = 1, n
      IF (coord(i) > extent_start(rect(i))) THEN
        coord(i) = coord(i) - 1
        nbidx(k) = rlcoord2lidx(rect, coord)
        coord(i) = coord(i) + 1
        k = k + 1
      END IF
      IF (coord(i) < extent_end(rect(i))) THEN
        coord(i) = coord(i) + 1
        nbidx(k) = rlcoord2lidx(rect, coord)
        coord(i) = coord(i) - 1
        k = k + 1
      END IF
    END DO
  END SUBROUTINE lidx_nb_indices_e

  SUBROUTINE lidx_nb_indices_i(rect, idx, nbidx)
    TYPE(iinterval), INTENT(in) :: rect(:)
    INTEGER(i4), INTENT(in) :: idx
    INTEGER(i4), INTENT(out) :: nbidx(:)

    INTEGER :: i, k, n
    INTEGER(i4) :: coord(SIZE(rect))

    n = SIZE(rect)
    coord = lidx2rlcoord(rect, idx)
    k = 1
    DO i = 1, n
      IF (coord(i) > extent_start(rect(i))) THEN
        coord(i) = coord(i) - 1
        nbidx(k) = rlcoord2lidx(rect, coord)
        coord(i) = coord(i) + 1
        k = k + 1
      END IF
      IF (coord(i) < extent_end(rect(i))) THEN
        coord(i) = coord(i) + 1
        nbidx(k) = rlcoord2lidx(rect, coord)
        coord(i) = coord(i) - 1
        k = k + 1
      END IF
    END DO
  END SUBROUTINE lidx_nb_indices_i
END MODULE ppm_rectilinear
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
