!> Contains structures and methods for turbulence config
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
MODULE mo_turb_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_turb_config

  TYPE, EXTENDS(t_jsb_config) :: t_turb_config
    REAL(wp) :: &
      & blending_height,             & !< Blending height [m]
      & roughness_snow,              & !< Roughness length of snow-covered surfaces [m]
      & roughness_bare,              & !< Roughness length of bare land surfaces [m]
      & roughness_water,             & !< Roughness length of lake open water [m]
      & roughness_lai_saturation,    & !< factor in roughness length dependence on LAI (should be between 0.1 and 1)
      & roughness_momentum_to_heat,  & !< factor by which roughness for heat is smaller than roughness for momentum
      & max_ini_rough_m                !< Limit roughness length read from ini file (<= 0: no limit)
    LOGICAL :: &
      & l_roughness_lai  !< true:  calculation of roughness length depends on LAI
                         !  false: PFT-specific values from lctlib are used
   CONTAINS
     PROCEDURE :: Init => Init_turb_config
  END type t_turb_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_turb_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize turbulence process
  !!
  !! @param[inout]     config     Configuration type of process (t_turb_config)
  !!
  SUBROUTINE Init_turb_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_turb_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    REAL(wp) :: &
      & blending_height, roughness_snow, roughness_bare, roughness_water, &
      & roughness_lai_saturation, roughness_momentum_to_heat, max_ini_rough_m
    LOGICAL :: l_roughness_lai

    NAMELIST /jsb_turb_nml/           &
      &   active,                     &
      &   ic_filename,                &
      &   bc_filename,                &
      &   blending_height,            &
      &   roughness_snow,             &
      &   roughness_bare,             &
      &   roughness_water,            &
      &   roughness_lai_saturation,   &
      &   roughness_momentum_to_heat, &
      &   max_ini_rough_m,            &
      &   l_roughness_lai

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_turb_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting turb configuration')

    ! Set defaults
    active                     = .TRUE.
    bc_filename                = 'bc_land_phys.nc'
    ic_filename                = 'ic_land_phys.nc'
    blending_height            = 100._wp
    roughness_snow             = 0.001_wp
    roughness_bare             = 0.005_wp
    roughness_water            = 0.0002_wp   ! Claussen, 1991: Estimation of areally-averaged surface fluxes
    roughness_lai_saturation   = 0.4_wp
    roughness_momentum_to_heat = 3._wp
    max_ini_rough_m            = -1._wp      ! No limit
    l_roughness_lai            = .TRUE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_turb_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_turb_nml)

    CALL close_nml(nml_handler)

    config%active                     = active
    config%ic_filename                = ic_filename
    config%bc_filename                = bc_filename
    config%blending_height            = blending_height
    config%roughness_snow             = roughness_snow
    config%roughness_bare             = roughness_bare
    config%roughness_water            = roughness_water
    config%roughness_lai_saturation   = roughness_lai_saturation
    config%roughness_momentum_to_heat = roughness_momentum_to_heat
    config%max_ini_rough_m            = max_ini_rough_m
    config%l_roughness_lai            = l_roughness_lai

  END SUBROUTINE Init_turb_config

#endif
END MODULE mo_turb_config_class
