!> Contains structures and methods for hydrology config
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
MODULE mo_hydro_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_jsb_impl_constants,ONLY: WB_IGNORE, WB_LOGGING, WB_ERROR
  USE mo_util,              ONLY: real2string
  USE mo_util_string,       ONLY: tolower
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hydro_config

  TYPE, EXTENDS(t_jsb_config) :: t_hydro_config
    CHARACTER(len=filename_max) :: bc_sso_filename   !< For elevation and oro_stddev
    CHARACTER(len=filename_max) :: ic_moist_filename !< For `l_read_initial_moist=.TRUE.`
    REAL(wp)               :: &
      & w_skin_max,           & !< Maximum water holding capacity of the soil and canopy skins, as well as
                                !  maximum snow depth on leaves  [m water equivalent]
      & w_soil_limit,         & !< Upper limit for maximum soil moisture content
      & w_soil_crit_fract,    & !< Fraction of field capacity at which soil starts to dry,
                                !  i.e. at which water stress starts to be above zero
      & w_soil_wilt_fract,    & !< Fraction of field capacity at which soil is at wilting point,
                                !  i.e. when no transpiration occurs anymore
      & snow_depth_max          !< Limit for snow depth [m water equivalent] to avoid grid cells with
                                !  infinitely growing snow depth; -1. for no limitation

    CHARACTER(len=10)     ::  &
      & scheme_stom_cond        !< Scheme to use for computation of stomatal conductance (echam5 only for now)
    LOGICAL               ::  &
      & l_read_initial_moist, & !< Force reading initial soil moisture for layers from `ic_moist_filename`
                                !  even if `init_from_ifs=.TRUE.`
      & sanitize_restart,     & !< Sanitize water budget in case a restart file does not exacly match
      & l_organic,            & !< Calculate soil hydrological parameters depending on organic matter fractions
      & l_socmap,             & !< Read organic matter fractions from 3D map
      & l_soil_texture,       & !< Deduce soil hydrological parameters from soil texture
      & l_glac_in_soil,       & !< Treat glacier as part of soil, not separate glacier tile
      & l_infil_subzero,      & !< Allow infiltration at soil temperatures below 0degC
      & l_ponds                 !< Additional surface water storage
    INTEGER               ::  &
      & soilhydmodel,         & !< Choice of soil hydrological model to calc soil hydraulic properties:
                                !! 'VanGenuchten' (default), 'BrooksCorey', or 'Campbell'
      & interpol_mean,        & !< Interpolation scheme, to optain K & D at layer interface:
                                !! 'upstream' (default) or 'arithmetic' mean
      & pond_dynamics,        & !< Choice of pond dynamics scheme (i.e. depth-area scaling)
                                !! 'tanh' (default), 'quadratic'
      & hydro_scale,          & !< 'semi_distributed': use the ARNO scheme assuming sub-grid scale variations in infiltration
                                !! and subsurface drainage (default)
                                !! 'uniform': assume uniform soil moisture distribution
                                !! (e.g. for applications like HydroTiles, HighRes or sitelevel)
      & enforce_water_budget    !< Set water balance check level: 'unset', 'ignore', 'logging', 'error'


  CONTAINS
    PROCEDURE             :: Init => Init_hydro_config
  END type t_hydro_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_config_class'

CONTAINS

  SUBROUTINE Init_hydro_config(config)

    USE mo_jsb_model_class,    ONLY: MODEL_JSBACH, MODEL_QUINCY
    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_hydro_constants,    ONLY: BrooksCorey_, Campbell_, VanGenuchten_, Upstream_, Arithmetic_, &
      &                              Quad_, Tanh_, Semi_Distributed_, Uniform_

    CLASS(t_hydro_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename, bc_sso_filename
    REAL(wp)                    :: w_skin_max, snow_depth_max
    REAL(wp)                    :: w_soil_limit, w_soil_crit_fract, w_soil_wilt_fract
    CHARACTER(len=16)           :: scheme_stom_cond, soilhydmodel, interpol_mean, pond_dynamics, hydro_scale
    CHARACTER(len=10)           :: enforce_water_budget
    LOGICAL                     :: l_read_initial_moist, sanitize_restart, l_organic, &
      &                            l_socmap, l_soil_texture, l_glac_in_soil, l_infil_subzero, l_ponds

    NAMELIST /jsb_hydro_nml/        &
      & active,                     &
      & ic_filename,                &
      & bc_filename,                &
      & bc_sso_filename,            &
      & w_skin_max,                 &
      & w_soil_limit,               &
      & w_soil_crit_fract,          &
      & w_soil_wilt_fract,          &
      & snow_depth_max,             &
      & scheme_stom_cond,           &
      & l_read_initial_moist,       &
      & sanitize_restart,           &
      & enforce_water_budget,       &
      & l_organic,                  &
      & l_socmap,                   &
      & l_soil_texture,             &
      & l_glac_in_soil,             &
      & l_infil_subzero,            &
      & l_ponds,                    &
      & soilhydmodel,               &
      & interpol_mean,              &
      & pond_dynamics,              &
      & hydro_scale

    INTEGER :: nml_handler, nml_unit, istat
    TYPE(t_input_file) :: input_file
    REAL(wp), POINTER :: ptr_1D(:)
    REAL(wp), ALLOCATABLE :: depths(:), mids(:)
    INTEGER :: nsoil

    TYPE(t_jsb_vgrid), POINTER :: soil_w

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hydro_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting hydro configuration')

    ! Set defaults
    SELECT CASE (config%model_config%model_scheme)
    CASE (MODEL_QUINCY)
      active                 = .FALSE.
    CASE (MODEL_JSBACH)
      active                 = .TRUE.
    END SELECT
    bc_filename            = 'bc_land_hydro.nc'
    ic_filename            = 'ic_land_hydro.nc'
    bc_sso_filename        = 'bc_land_sso.nc'
    w_skin_max             = 2.E-4_wp            ! 0.2 mm corresponds to Roesch et al. 2001
    w_soil_limit           = -1.0_wp
    w_soil_crit_fract      = 0.75_wp
    w_soil_wilt_fract      = 0.35_wp
    snow_depth_max         = -1.0_wp             ! -1: no limitation
    scheme_stom_cond       = 'echam5'
    l_read_initial_moist   = .FALSE.
    sanitize_restart       = .FALSE.
    enforce_water_budget   = 'unset'
    l_organic              = .TRUE.
    l_socmap               = .TRUE.
    l_soil_texture         = .FALSE.
    l_glac_in_soil         = .FALSE.
    l_infil_subzero        = .TRUE.
    l_ponds                = .FALSE.
    pond_dynamics          = "tanh"
    soilhydmodel           = "VanGenuchten"
    interpol_mean          = "upstream"
    hydro_scale            = "semi_distributed"

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_hydro_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_hydro_nml)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%bc_sso_filename     = bc_sso_filename

    config%w_skin_max          = w_skin_max
    CALL message(TRIM(routine), 'Maximum water holding capacity of soil or canopy skin: '//TRIM(real2string(w_skin_max)))

    config%w_soil_limit        = w_soil_limit
    CALL message(TRIM(routine), 'Upper limit for maximum soil moisture content: '//TRIM(real2string(w_soil_limit)))

    config%w_soil_crit_fract   = w_soil_crit_fract
    CALL message(TRIM(routine), 'Fraction of field capacity at critical point: '//TRIM(real2string(w_soil_crit_fract)))

    config%w_soil_wilt_fract   = w_soil_wilt_fract
    CALL message(TRIM(routine), 'Fraction of field capacity at permanent wilting point: '//TRIM(real2string(w_soil_wilt_fract)))

    config%snow_depth_max      = snow_depth_max
    IF (snow_depth_max >= 0._wp) THEN
      CALL message(TRIM(routine), 'Snow depth limitation: '//TRIM(real2string(snow_depth_max)))
    END IF

    config%scheme_stom_cond       = scheme_stom_cond

    config%l_read_initial_moist   = l_read_initial_moist

    config%sanitize_restart       = sanitize_restart
    IF (config%sanitize_restart) THEN
      CALL message(TRIM(routine), 'Soil hydrology of the restart file is sanitized.')
    END IF

    SELECT CASE (tolower(TRIM(enforce_water_budget)))
    CASE ("unset")
      config%enforce_water_budget = config%model_config%enforce_water_budget
    CASE ("ignore")
      config%enforce_water_budget = WB_IGNORE
      CALL message(TRIM(routine), 'WARNING: Land surface water balance will not be checked during simulation.')
    CASE ("logging")
      config%enforce_water_budget = WB_LOGGING
      CALL message(TRIM(routine), 'WARNING: Simulation will not stop due to any land surface water balance violation '// &
        & 'but information will be added to the log file')
    CASE ("error")
      config%enforce_water_budget = WB_ERROR
      CALL message(TRIM(routine), 'WARNING: Simulation will stop due to any land surface water balance violation.')
    CASE DEFAULT
      CALL finish(TRIM(routine), 'enforce_water_budget == '//tolower(TRIM(enforce_water_budget))//' not available.')
    END SELECT
    IF (config%enforce_water_budget == WB_ERROR .AND. sanitize_restart) THEN
      CALL message(TRIM(routine), 'WARNING: Sanitizing restart files causes water balance violations.' &
        &                       //' You should switch off enforce_water_budget in case sanitizing is needed.')
    END IF

    config%l_soil_texture = l_soil_texture
    IF (l_soil_texture) THEN
      CALL message(TRIM(routine), 'Calculate hydrological parameters for mineral soils from soil texture.')
    END IF

    config%l_organic = l_organic
    IF (l_organic) THEN
      CALL message(TRIM(routine), 'Include organic matter in calculation of hydrological soil parameters.')
    ELSE
      IF (l_soil_texture) THEN
        CALL message (TRIM(routine), 'WARNING: l_organic=.F. and l_soil_texture=.T.' &
                                   //' => Hydrological soil parameters represent purely mineral soils!')
      ELSE
        CALL message (TRIM(routine), 'Soil hydrological parameters read from bc file.')
      END IF
    END IF

    config%l_socmap = l_socmap
    IF (l_socmap .AND. l_organic) THEN
      CALL message(TRIM(routine), 'Read map with soil organic carbon fractions from bc file.')
    ELSE IF (.NOT. l_socmap .AND. l_organic) THEN
      CALL message(TRIM(routine), 'Calculate soil organic carbon fractions from forest fractions.')
    ELSE IF (l_socmap .AND. .NOT. l_organic) THEN
      CALL message(TRIM(routine), 'Setting l_socmap=false since l_organic=false and organic matter fractions are not used.')
      config%l_socmap = .FALSE.
    END IF

    config%l_glac_in_soil      = l_glac_in_soil
    config%l_infil_subzero     = l_infil_subzero
    IF (l_infil_subzero) CALL message(TRIM(routine), 'Allow infiltration at sub-zero temperatures')

    config%l_ponds = l_ponds
    IF (l_ponds) THEN
      CALL message(TRIM(routine), '*** WARNING: Use of ponds is still experimental and not thoroughly validated!!')
    ELSE
      CALL message(TRIM(routine), 'Experimental pond scheme is disabled.')
    END IF

    SELECT CASE(TRIM(tolower(soilhydmodel)))
    CASE("vangenuchten")
      config%soilhydmodel = VanGenuchten_
      CALL message(TRIM(routine), 'Soil hydrology uses the Van Genuchten model')
    CASE("brookscorey")
      config%soilhydmodel = BrooksCorey_
      CALL message(TRIM(routine), 'Soil hydrology uses the Brooks & Corey model')
    CASE("campbell")
      config%soilhydmodel = Campbell_
      CALL message(TRIM(routine), 'Soil hydrology uses the Campbell model')
    CASE default
      CALL finish(TRIM(routine), 'Soil hydrology model '//TRIM(soilhydmodel)//' not implemented')
    END SELECT

    SELECT CASE(TRIM(tolower(interpol_mean)))
    CASE("upstream")
      config%interpol_mean = Upstream_
      CALL message(TRIM(routine), 'Soil diffusivity/conductivity interpolation uses upstream means')
    CASE("arithmetic")
      config%interpol_mean = Arithmetic_
      CALL message(TRIM(routine), 'Soil diffusivity/conductivity interpolation uses arithmetic means')
    CASE default
      CALL finish(TRIM(routine), 'Soil diffusivity/conductivity interpolation method '//TRIM(interpol_mean)//' not implemented')
    END SELECT

    SELECT CASE(TRIM(tolower(pond_dynamics)))
    CASE("quadratic")
      config%pond_dynamics = Quad_
      CALL message(TRIM(routine), 'Pond dynamics scheme uses quadratic scaling')
    CASE("tanh")
      config%pond_dynamics = Tanh_
      CALL message(TRIM(routine), 'Pond dynamics scheme uses tanh scaling')
    CASE default
      CALL finish(TRIM(routine), 'Pond dynamics scheme '//TRIM(pond_dynamics)//' not implemented')
    END SELECT

    SELECT CASE(TRIM(tolower(hydro_scale)))
    CASE("semi_distributed")
      config%hydro_scale = Semi_Distributed_
      CALL message(TRIM(routine), 'Using semi-distributed scale parametrization (ARNO scheme) for soil hydrological fluxes')
    CASE("uniform")
      config%hydro_scale = Uniform_
      CALL message(TRIM(routine), 'Using uniform scale process implementation for hydrological fluxes (still experimental)')
    CASE default
      CALL finish(TRIM(routine), 'Soil hydrology scale '//TRIM(hydro_scale)//' not implemented')
    END SELECT

    CALL close_nml(nml_handler)

    IF (.NOT. active) RETURN

    input_file = jsb_netcdf_open_input(ic_filename)

    ptr_1D => input_file%Read_1d(variable_name='soillev')    ! Depth of layer bottom
    nsoil = SIZE(ptr_1D)

    IF (l_organic .AND. nsoil < 3) &
      & CALL finish(TRIM(routine), 'At least three soil layers required for using organic soil component')

    ALLOCATE(depths(nsoil+1))
    ALLOCATE(mids(nsoil))
    depths(1) = 0._wp
    depths(2:nsoil+1) = ptr_1D(1:nsoil)
    mids(1:nsoil) = (depths(1:nsoil) + depths(2:nsoil+1)) / 2._wp
    DEALLOCATE(ptr_1D)

    CALL input_file%Close()

    soil_w  => new_vgrid('soil_depth_water', ZAXIS_DEPTH_BELOW_LAND, nsoil, &
      & levels  = mids                 (1:nsoil  ),                         &
      & lbounds = depths               (1:nsoil  ),                         &
      & ubounds = depths               (2:nsoil+1),                         &
      & units='m')
    CALL register_vgrid(soil_w)
    WRITE(message_text, *) 'Soil levels in hydrology (upper) [m]: ', soil_w%lbounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in hydrology (mid)   [m]: ', soil_w%levels
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in hydrology (lower) [m]: ', soil_w%ubounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil level depths in hydrology   [m]: ', soil_w%dz
    CALL message(TRIM(routine), message_text)
    DEALLOCATE(depths, mids)

  END SUBROUTINE Init_hydro_config

#endif
END MODULE mo_hydro_config_class
