! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

PROGRAM test_divide_cell

  USE mo_decomposition_tools, ONLY: t_cell_info, divide_cells_by_location

  IMPLICIT NONE

  INTEGER, PARAMETER :: ncells = 10
  INTEGER, PARAMETER :: nparts = 3, firstpart = 0
  INTEGER :: i, owner
  LOGICAL, ALLOCATABLE :: seen(:)

  TYPE(t_cell_info) :: cell_desc(ncells) &
       = (/ t_cell_info(-5000, -1000, 1, HUGE(1)), &
       &    t_cell_info(5000, -1000, 2, -HUGE(1)), &
       &    t_cell_info(-5000, 1000, 3, HUGE(1)), &
       &    t_cell_info(5000, 1000, 4, -HUGE(1)), &
       &    t_cell_info(-5000, 0, 5, HUGE(1)), &
       &    t_cell_info(-5000, -2000, 6, HUGE(1)), &
       &    t_cell_info(5000, -2000, 7, -HUGE(1)), &
       &    t_cell_info(-5000, 2000, 8, HUGE(1)), &
       &    t_cell_info(5000, 2000, 9, -HUGE(1)), &
       &    t_cell_info(-5000, 0, 10, HUGE(1)) /)

  CALL divide_cells_by_location(ncells, cell_desc, &
       firstpart, nparts - 1 + firstpart)

  ALLOCATE(seen(0:MAX(ncells, nparts)))
  seen = .FALSE.
  seen(firstpart:nparts+firstpart-1) = .FALSE.
  DO i = 1, ncells
    owner = cell_desc(i)%owner
    IF (owner < firstpart .OR. owner >= nparts + firstpart) THEN
      PRINT *, ' cell ', i, 'found invalid owner assignment!'
      PRINT *, TRIM(cell_desc2str(cell_desc(i)))
      STOP 1
    END IF
    seen(owner) = .TRUE.
  END DO
  DO i = firstpart, nparts+firstpart-1
    IF (.NOT. seen(i)) THEN
      PRINT *, 'no cells assigned to rank ', i, '!'
      STOP 1
    END IF
  END DO

#if 0
  PRINT *, ''
  DO i = 1, ncells
    PRINT *, TRIM(cell_desc2str(cell_desc(i)))
  END DO
#endif

CONTAINS
  FUNCTION cell_desc2str(c) RESULT(s)
    TYPE(t_cell_info), INTENT(in) :: c
    INTEGER, PARAMETER :: idig = RANGE(1) + 1 + 1
    CHARACTER(4 * idig + (4 - 1) * 2) :: s
    WRITE (s, '(3(i0,", "),i0)') c%lat, c%lon, c%cell_number, c%owner
  END FUNCTION cell_desc2str
END PROGRAM test_divide_cell
