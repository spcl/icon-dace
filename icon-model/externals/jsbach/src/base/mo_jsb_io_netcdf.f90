!> Contains basic functions for NetCDF I/O
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
MODULE mo_jsb_io_netcdf
#ifndef __NO_JSBACH__

  USE mo_jsb_io_netcdf_iface, ONLY: t_input_file, netcdf_open_input
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_grid,           ONLY: Get_grid

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_input_file, jsb_netcdf_open_input

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_netcdf'

CONTAINS

  TYPE(t_input_file) FUNCTION jsb_netcdf_open_input(filename, grid_id)

    CHARACTER(LEN=*), INTENT(in)           :: filename
    INTEGER,          INTENT(in), OPTIONAL :: grid_id

    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_netcdf_open_input'

    IF (PRESENT(grid_id)) THEN
      grid => Get_grid(grid_id)
      jsb_netcdf_open_input = netcdf_open_input(filename, grid%patch)
    ELSE
      jsb_netcdf_open_input = netcdf_open_input(filename)
    END IF

  END FUNCTION jsb_netcdf_open_input

#endif
END MODULE mo_jsb_io_netcdf
