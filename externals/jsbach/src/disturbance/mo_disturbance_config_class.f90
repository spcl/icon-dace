!> Contains structures and methods for disturbance config
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
MODULE mo_disturb_config_class
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
  PUBLIC :: t_disturb_config

  TYPE, EXTENDS(t_jsb_config) :: t_disturb_config
    INTEGER           :: fire_algorithm      ! Select fire algorithm: 0 : none, 1 : jsbach (default),
                                             ! 2 : Arora and Boer (2005) (not yet), 3 : Thonicke et al.(2001)
    INTEGER           :: windbreak_algorithm ! Select windbreak algorithm: 0 : none, 1 : jsbach (default)
    LOGICAL           :: lburn_pasture
    REAL(wp)          :: wnd_threshold, &    ! factor by which the previous day maximum wind speed
                                             ! must be larger than the climatological daily maximum
                                             ! wind speed to allow any wind damage
                         wnd_damage_scale    ! scaling factor for wind damage

   REAL(wp)           :: fire_litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                         fire_rel_hum_threshold, & ! maximal relative humidity for fire
                         fire_minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                         fire_minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                         fire_tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                         fire_tau_grass            ! return period of fire for grass PFT [year] at 0% relative humidity

   CONTAINS
     PROCEDURE :: Init => Init_disturb_config
  END type t_disturb_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize disturb process
  !!
  !! @param[inout]     config     Configuration type of process (t_disturb_config)
  !!
  SUBROUTINE Init_disturb_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_disturb_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    INTEGER                     :: fire_algorithm
    INTEGER                     :: windbreak_algorithm
    LOGICAL                     :: lburn_pasture
    ! Fire namelist parameters:
    REAL(wp)                    :: fire_litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                                   fire_rel_hum_threshold, & ! maximal relative humidity for fire
                                   fire_minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                                   fire_minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                                   fire_tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                                   fire_tau_grass            ! return period of fire for grass PFT [year] at 0% relative humidity
    ! Wind break namelist parameters:
    REAL(wp)                    :: wnd_threshold, & ! factor by which the previous day maximum wind speed
                                                    ! must be larger than the climatological daily maximum
                                                    ! wind speed to allow any wind damage
                                   wnd_damage_scale ! scaling factor for wind damage

    NAMELIST /jsb_disturb_nml/     &
         active,                   &
         ic_filename,              &
         bc_filename,              &
         fire_algorithm,           &
         windbreak_algorithm,      &
         lburn_pasture,            &
         fire_litter_threshold,    &
         fire_rel_hum_threshold,   &
         fire_minimum_woody,       &
         fire_minimum_grass,       &
         fire_tau_woody,           &
         fire_tau_grass,           &
         wnd_threshold,            &
         wnd_damage_scale

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_disturb_config'

    CALL message(TRIM(routine), 'Starting disturb configuration')

    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_disturb.nc'
    ic_filename            = 'ic_land_disturb.nc'
    fire_algorithm         =  1
    windbreak_algorithm    =  1
    lburn_pasture          = .FALSE.
    !
    fire_litter_threshold  = 27.87_wp   ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                                        ! In JS3 this value was 16.67 because it was related to fuel per gridbox area.
                                        ! In JS4 it is related to fuel per potential vegetated area of the veg tile,
                                        ! which is on global average 0.5980 of the grid box area.
    fire_rel_hum_threshold = 70._wp     ! maximal relative humidity for fire
    fire_minimum_woody     =  0.002_wp  ! minimal fraction of act_fpc of woody PFT to be burned each year
    fire_minimum_grass     =  0.006_wp  ! minimal fraction of act_fpc of grass PFT to be burned each year
    fire_tau_woody         =  6.0_wp    ! return period of fire for woody PFT [year] at 0% relative humidity
    fire_tau_grass         =  2.0_wp    ! return period of fire for grass PFT [year] at 0% relative humidity
    !
    wnd_threshold          =  2.25_wp   ! factor by which the previous day maximum wind speed
                                        ! must be larger than the climatological daily maximum
                                        ! wind speed to allow any wind damage
    wnd_damage_scale       =  5.e-03_wp ! scaling factor for wind damage

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_disturb_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_disturb_nml)

    CALL close_nml(nml_handler)

    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%fire_algorithm         = fire_algorithm
    config%windbreak_algorithm    = windbreak_algorithm
    config%lburn_pasture          = lburn_pasture
    config%fire_litter_threshold  = fire_litter_threshold
    config%fire_rel_hum_threshold = fire_rel_hum_threshold
    config%fire_minimum_woody     = fire_minimum_woody
    config%fire_minimum_grass     = fire_minimum_grass
    config%fire_tau_woody         = fire_tau_woody
    config%fire_tau_grass         = fire_tau_grass
    config%wnd_threshold          = wnd_threshold
    config%wnd_damage_scale       = wnd_damage_scale

  END SUBROUTINE Init_disturb_config

#endif
END MODULE mo_disturb_config_class
