!> Contains structures and methods for radiation config
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
MODULE mo_rad_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_rad_config

  ! ======================================================================================================= !
  !> RAD_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_rad_config
    LOGICAL ::                     &
      & use_alb_veg_simple,        &  !< Use the simplified albedo calculations for vegetation tiles
      & use_alb_soil_scheme,       &  !< For vegetation and bare soil tiles:
                                      !   Use the complex calculations of the soil albedo from carbon in and
                                      !   on the soil (accounts for the soil under vegetation and for bare soil)
                                      !   if vegetation tile: only used if use_alb_veg_simple=.FALSE.
      & use_alb_mineralsoil_const, &  !< For vegetation and bare soil tiles:
                                      !   if TRUE constant parameter values are used for albedo of the mineral soil
                                      !   if FALSE values for albedo of the mineral soil are taken from
                                      !   bc_land_phys.nc file.
      & use_alb_soil_litter,       &  !< For vegetation and bare soil tiles:
                                      !   TRUE: soil albedo includes the albedo of the litter cover above it
                                      !   FALSE soil albedo ignores the albedo of the litter cover above it

      & use_alb_canopy                !< Only for vegetation tiles:
                                      !   TRUE:  canopy albedo is taken from bc_land_phys.nc file
                                      !   FALSE: canopy albedo is taken from lctlib file
    CHARACTER(len=10) :: &
      & use_alb_soil_organic_C        !< For vegetation and bare soil tiles:
                                      !   if linear: the albedo of soil without covering litter but with carbon in it
                                      !     is calculated from the mineral soil albedo minus a reduction calculated
                                      !     linearily from the carbon in the soil
                                      !   if log: the albedo of soil without covering litter but with carbon in it
                                      !     is calculated from a logarithmically decreasing function shaped by
                                      !     the mineral soil albedo and the carbon content within the soil
                                      !     linearily from the carbon in the soil
    REAL(wp) :: &
      & albedo_age_weight, &          !< 0: snow age scheme with zero weight: only temperature scheme for snow albedo
                                      !      is used
                                      ! 1: only snow age scheme is used, temperature scheme with zero weight
                                      ! 0 < albedo_age_weight < 1: snow albedo is calculated by linearly weighting
                                      !                            the snow albedo resulting from both schemes
      & AlbedoCanopySnow              ! Snow albedo for snow on canopy. This configuration value is not read in from
                                      ! the namelist. It is calculated in subroutine rad_init_ic, but should nevertheless
                                      ! placed here!
  CONTAINS
    PROCEDURE :: Init => Init_rad_config
  END type t_rad_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read RAD_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_rad_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_rad_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: use_alb_veg_simple
    LOGICAL                     :: use_alb_mineralsoil_const
    CHARACTER(len=10)           :: use_alb_soil_organic_C
    LOGICAL                     :: use_alb_soil_litter
    LOGICAL                     :: use_alb_soil_scheme
    LOGICAL                     :: use_alb_canopy
    REAL(wp)                    :: albedo_age_weight

    NAMELIST /jsb_rad_nml/            &
         active,                      &
         ic_filename,                 &
         bc_filename,                 &
         use_alb_veg_simple,          &
         use_alb_mineralsoil_const,   &
         use_alb_soil_organic_C,      &
         use_alb_soil_litter,         &
         use_alb_soil_scheme,         &
         use_alb_canopy,              &
         albedo_age_weight

    INTEGER :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_rad_config'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Starting radiation configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active                    = .TRUE.
    bc_filename               = 'bc_land_rad.nc'
    ic_filename               = 'ic_land_rad.nc'
    use_alb_veg_simple        = .FALSE.
    use_alb_mineralsoil_const = .FALSE.
    use_alb_soil_organic_C    = ''
    use_alb_soil_litter       = .FALSE.
    use_alb_soil_scheme       = .FALSE.
    use_alb_canopy            = .FALSE.
    albedo_age_weight         = 0.5_wp
    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_rad_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_rad_nml)

    CALL close_nml(nml_handler)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%use_alb_veg_simple        = use_alb_veg_simple
    config%use_alb_mineralsoil_const = use_alb_mineralsoil_const
    config%use_alb_soil_organic_C    = use_alb_soil_organic_C
    config%use_alb_soil_litter       = use_alb_soil_litter
    config%use_alb_soil_scheme       = use_alb_soil_scheme
    config%use_alb_canopy            = use_alb_canopy
    config%albedo_age_weight         = albedo_age_weight

  END SUBROUTINE Init_rad_config

#endif
END MODULE mo_rad_config_class
