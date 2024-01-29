! Module defining the derived types for the selection of vertical
! levels for the output module.
!
! Details of the implementation: see "mo_level_selection".
!
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

MODULE mo_level_selection_types

  USE mo_impl_constants,    ONLY: SUCCESS
  USE mo_exception,         ONLY: finish
  IMPLICIT NONE

  PRIVATE

  ! subroutines
  PUBLIC :: t_level_selection

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_level_selection_types'


  TYPE t_level_selection
    !> number of selected levels
    INTEGER :: n_selected

    !> level selection as input in the form of a LOGICAL array
    !  s(1...N), where "s(i)=.TRUE." means that level "i" is selected:
    LOGICAL, ALLOCATABLE :: s(:)

    !> integer list idx(1...n_selected) containing the selected level
    !  indices.
    INTEGER, ALLOCATABLE :: global_idx(:)

    !> local index in the list of selected level indices (i.e. an
    !  integer number in the range 1...n_selected).
    INTEGER, ALLOCATABLE :: local_idx(:)

  CONTAINS
    PROCEDURE :: finalize => t_level_selection_finalize        !< destructor
  END TYPE t_level_selection

CONTAINS


  ! --------------------------------------------------------------------------------------
  !> Frees a level selection data object.
  !
  SUBROUTINE t_level_selection_finalize(selection)
    CLASS(t_level_selection), INTENT(INOUT) :: selection !< level selection object
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = '::deallocate_level_selection'
    INTEGER :: ierrstat

    IF (ALLOCATED(selection%s))  DEALLOCATE(selection%s, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    IF (ALLOCATED(selection%global_idx))  DEALLOCATE(selection%global_idx, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    IF (ALLOCATED(selection%local_idx))  DEALLOCATE(selection%local_idx, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE t_level_selection_finalize

END MODULE mo_level_selection_types

