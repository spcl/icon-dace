!> Contains functions for grid instances
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_grid
#ifndef __NO_JSBACH__

  USE mo_exception,      ONLY: message, message_text, finish
  USE mo_jsb_grid_class, ONLY: t_jsb_grid, t_jsb_grid_ptr, t_jsb_vgrid, t_jsb_vgrid_ptr

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_grid, Register_vgrid, Get_grid, Get_vgrid

  TYPE(t_jsb_grid_ptr),  ALLOCATABLE, SAVE :: jsbach_grids(:)

CONTAINS

  SUBROUTINE Register_grid(grid)

    TYPE(t_jsb_grid),  POINTER, INTENT(in)    :: grid

    INTEGER :: i, n
    TYPE(t_jsb_grid_ptr), ALLOCATABLE :: grids_temp(:)

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Register_grid'

    IF (ALLOCATED(jsbach_grids)) THEN
      DO i=1,SIZE(jsbach_grids)
        IF (ASSOCIATED(jsbach_grids(i)%p)) THEN
          IF (jsbach_grids(i)%p%id == grid%id) THEN
            RETURN ! grid already exists and is already registered
          END IF
        END IF
      END DO
      ! grid is new
      n = SIZE(jsbach_grids)
      ALLOCATE(grids_temp(n+1))
      grids_temp(1:n) = jsbach_grids
      grids_temp(n+1) = t_jsb_grid_ptr(grid)
      CALL MOVE_ALLOC(grids_temp, jsbach_grids)
    ELSE
      ALLOCATE(jsbach_grids(1))
      jsbach_grids(1) = t_jsb_grid_ptr(grid)
    END IF

    CALL message(TRIM(routine), 'registered new JSBACH grid instance '//TRIM(grid%name))

    !$ACC ENTER DATA PCREATE(grid, grid%lat, grid%lon, grid%area, grid%lsm, grid%lsf)
    !$ACC UPDATE DEVICE(grid%lat, grid%lon, grid%area, grid%lsm, grid%lsf)

  END SUBROUTINE Register_grid

  SUBROUTINE Register_vgrid(vgrid)

    USE mo_jsb_io_iface, ONLY: Create_zaxis
    USE mo_jsb_vertical_axes, ONLY: jsbach_vgrids

    TYPE(t_jsb_vgrid), POINTER, INTENT(inout) :: vgrid

    TYPE(t_jsb_vgrid_ptr), ALLOCATABLE :: vgrids_temp(:)
    INTEGER :: i, n, length

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Register_vgrid'

    IF (ALLOCATED(jsbach_vgrids)) THEN
      DO i=1,SIZE(jsbach_vgrids)
        IF (ASSOCIATED(jsbach_vgrids(i)%p)) THEN
          IF (jsbach_vgrids(i)%p%ZaxisID == vgrid%ZaxisID) THEN
            RETURN ! axis already exists and is already registered
          END IF
        END IF
      END DO
      ! Axis is new
      n = SIZE(jsbach_vgrids)
      ALLOCATE(vgrids_temp(n+1))
      vgrids_temp(1:n) = jsbach_vgrids
      vgrids_temp(n+1) = t_jsb_vgrid_ptr(vgrid)
      CALL MOVE_ALLOC(vgrids_temp, jsbach_vgrids)
    ELSE
      ALLOCATE(jsbach_vgrids(1))
      jsbach_vgrids(1) = t_jsb_vgrid_ptr(vgrid)
    END IF
    vgrid%id = SIZE(jsbach_vgrids)

    length = vgrid%n_levels
    CALL Create_zaxis(vgrid%ZaxisID, vgrid%echamZaxisIdx, &
      & vgrid%cdi_axis_type, length, vgrid%name, vgrid%levels(1:length), vgrid%units, vgrid%longname, &
      & vgrid%lbounds(1:length), vgrid%ubounds(1:length))

    WRITE(message_text, '(A,A,A,I10,I3)') 'Registered new vgrid: name, ID, length= ',TRIM(vgrid%name), ' ', vgrid%ZaxisID, length
    CALL message(TRIM(routine), message_text)

    !$ACC ENTER DATA CREATE(vgrid, vgrid%levels, vgrid%lbounds, vgrid%ubounds, vgrid%dz)
    !$ACC UPDATE DEVICE(vgrid%levels, vgrid%lbounds, vgrid%ubounds, vgrid%dz)
  END SUBROUTINE Register_vgrid

  FUNCTION Get_grid(id) RESULT(grid)

    INTEGER, INTENT(in)        :: id
    TYPE(t_jsb_grid), POINTER  :: grid

    INTEGER :: i
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Get_grid'

    grid => NULL()
    DO i=1,SIZE(jsbach_grids)
      IF (ASSOCIATED(jsbach_grids(i)%p)) THEN
        IF (jsbach_grids(i)%p%id == id) THEN
          grid => jsbach_grids(i)%p
          EXIT
        END IF
      END IF
    END DO

    IF (.NOT. ASSOCIATED(grid)) THEN
      CALL finish(TRIM(routine), 'Grid not found.')
    END IF

  END FUNCTION Get_grid

  FUNCTION Get_vgrid(name) RESULT(vgrid)

    USE mo_jsb_vertical_axes, ONLY: jsbach_vgrids

    CHARACTER(len=*),           INTENT(in)    :: name
    TYPE(t_jsb_vgrid), POINTER                :: vgrid

    INTEGER :: i
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Get_vgrid'

    vgrid => NULL()
    DO i=1,SIZE(jsbach_vgrids)
      IF (ASSOCIATED(jsbach_vgrids(i)%p)) THEN
        IF (jsbach_vgrids(i)%p%id > 0 .AND. TRIM(jsbach_vgrids(i)%p%name) == TRIM(name)) THEN
          vgrid => jsbach_vgrids(i)%p
          EXIT
        END IF
      END IF
    END DO

    IF (.NOT. ASSOCIATED(vgrid)) THEN
      CALL finish(TRIM(routine), 'Vgrid '//TRIM(name)//' not found.')
    END IF

  END FUNCTION Get_vgrid

#endif
END MODULE mo_jsb_grid
