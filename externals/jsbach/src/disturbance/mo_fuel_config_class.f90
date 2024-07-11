!> Contains structures and methods for fuel config
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
MODULE mo_fuel_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_fuel_config

  TYPE, EXTENDS(t_jsb_config) :: t_fuel_config
    INTEGER :: fuel_algorithm  ! Select fuel algorithm: 1 : jsbach (default), 2 : Arora and Boer (2005),
                               !                        3 : Thonicke et al.(2001), 4 : Read/GFED
  CONTAINS
     PROCEDURE :: Init => Init_fuel_config
  END type t_fuel_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize fuel process
  !!
  !! @param[inout]     config     Configuration type of process (t_fuel_config)
  !!
  SUBROUTINE Init_fuel_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_fuel_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    INTEGER                     :: fuel_algorithm

    NAMELIST /jsb_fuel_nml/     &
        & active,               &
        & ic_filename,          &
        & bc_filename,          &
        & fuel_algorithm

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_fuel_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting fuel configuration')

    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_fuel.nc'
    ic_filename            = 'ic_land_fuel.nc'
    fuel_algorithm         =  1

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_fuel_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_fuel_nml)

    CALL close_nml(nml_handler)

    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%fuel_algorithm         = fuel_algorithm

  END SUBROUTINE Init_fuel_config

#endif
END MODULE mo_fuel_config_class
