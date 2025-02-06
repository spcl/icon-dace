!> Contains structures and methods for carbon config
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
MODULE mo_carbon_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_carbon_config

  TYPE, EXTENDS(t_jsb_config) :: t_carbon_config
    LOGICAL  :: with_nitrogen
    LOGICAL  :: read_cpools ! Read Cpools from external file when experiment is newly started (cpools are
                            ! initialized with values from file) or restarted (cpools from restart file are
                            ! overwritten by values from file).
                            ! ATTENTION: The run script has to take care that read_cpools is set to .FALSE.
                            !            in the experiment jobs following the first one! Or re-set it manually!
    LOGICAL  :: diag_humus_fluxes ! Activate diagnostic carbon fluxes to and from the yasso humus pools
                                  ! required for analytical humus pool equilibration (compare jsb3.2 svn rev 9934)
    REAL(wp) :: fire_fract_wood_2_atmos

   CONTAINS
     PROCEDURE :: Init => Init_carbon_config
  END type t_carbon_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize carbon process
  !!
  !! @param[inout]     config     Configuration type of process (t_carbon_config)
  !!
  SUBROUTINE Init_carbon_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_carbon_constants,   ONLY: FireFracWood2Atmos

    CLASS(t_carbon_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: with_nitrogen
    LOGICAL                     :: read_cpools, diag_humus_fluxes
    REAL(wp)                    :: fire_fract_wood_2_atmos
    ! LOGICAL           :: lemissions

    NAMELIST /jsb_carbon_nml/      &
         active,                   &
         ic_filename,              &
         bc_filename,              &
         with_nitrogen,            &
         read_cpools,              &
         diag_humus_fluxes,        &
         fire_fract_wood_2_atmos

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_carbon_config'

    CALL message(TRIM(routine), 'Starting carbon configuration')

    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_carbon.nc'
    ic_filename            = 'ic_land_carbon.nc'
    with_nitrogen          = .FALSE.
    read_cpools            = .FALSE.
    diag_humus_fluxes      = .FALSE.
    fire_fract_wood_2_atmos = FireFracWood2Atmos

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_carbon_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_carbon_nml)

    CALL close_nml(nml_handler)

    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%with_nitrogen          = with_nitrogen
    config%read_cpools            = read_cpools
    config%diag_humus_fluxes      = diag_humus_fluxes
    config%fire_fract_wood_2_atmos = fire_fract_wood_2_atmos

  END SUBROUTINE Init_carbon_config

#endif
END MODULE mo_carbon_config_class
