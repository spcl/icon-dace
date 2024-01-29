!> Contains the memory class for the hydrology process.
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
MODULE mo_hydro_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_control,      ONLY: jsbach_runs_standalone
  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,    ONLY: LAND_TYPE, BARE_TYPE, VEG_TYPE, LAKE_TYPE, GLACIER_TYPE
  USE mo_jsb_varlist,      ONLY: BASIC !, MEDIUM, FULL
  USE mo_jsb_process_class, ONLY: ON_TILE_, AGGREGATE_

  ! Use of prcesses in this module
  dsl4jsb_Use_processes HYDRO_, SSE_, SEB_, ASSIMI_

  ! Use process configurations
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(ASSIMI_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hydro_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 90

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_hydro_memory

    TYPE(t_jsb_var_real2d) :: &
      & fract_snow,           & !< Snow fraction of tile                                                []
      & w_snow,               & !< Water content of snow reservoir                                      [m water equivalent]
      & evapo_snow,           & !< Evaporation from snow                                                [kg m-2 s-1]
      & snowmelt,             & !< Snow melt                                                            [kg m-2 s-1]
      & evapopot,             & !< Potential evaporation                                                [kg m-2 s-1]
      & evapotrans,           & !< Evapotranspiration independent of the type of tile                   [kg m-2 s-1]
                                !< e.g. for the lake tile it is evapopot and for the land tile it is evapotrans_lnd.
      & water_flux,           & !< Total water flux                                                      [m3/s]
      & water_content,        & !< Total water content                                                   [m3]
      & water_budget            !< Water budget error                                                    [m3]

    TYPE(t_jsb_var_real2d) :: &
      & runoff,               & !< surface runoff
      & drainage,             & !< drainage
      & ws_negative,          & !< sum of negative soil moisture (diagnosis of errors) [m]
      & discharge,            & !< discharge (local)
      & discharge_ocean,      & !< discharge to the ocean
      & internal_drain          !< internal drain

    TYPE(t_jsb_var_real2d) :: &
      & elevation,            & !< Elevation
      & oro_stddev              !< Standard deviation of orography

    ! Additional variables for land type
    TYPE(t_jsb_var_real2d) :: &
      & fract_snow_soil,      & !< Snow fraction on soil                                               []
      & w_snow_soil,          & !< Water content of snow reservoir of soil                             [m water equivalent]
      ! & w_snow_soil_old,      & !< Water content of snow reservoir of soil (old)                     [m water equivalent]
      & snow_soil_dens,       & !< Density of snow on soil                                             [kg m-3]
      & evapotrans_lnd,       & !< Evapotranspiration over vegetation and bare soil type               [kg m-2 s-1]
      & q_snocpymlt             !< Heating by melting of snow on canopy                                [W m-2]

    ! Additional variables for soil

    ! Parameters
    TYPE(t_jsb_var_real3d) :: &
      & soil_depth_sl,        & !< Depth of each soil layer that can be saturated with water (until bedrock) [m]
      & fract_org_sl            !< Fractions of organic material in soil layers                        []
    TYPE(t_jsb_var_real2d) :: &
      & soil_depth,           & !< Soil depth derived from textures (bedrock)                          [m]
      & max_moist,            & !< Maximum root zone moisture content                                  [m]
      & vol_field_cap,        & !< Volumetric soil field capacity                                      [m/m]
      & vol_p_wilt,           & !< Volumetric permanent wilting point                                  [m/m]
      & vol_porosity,         & !< Volumetric porosity of mineral soil                                 [m/m]
      & pore_size_index,      & !< Soil pore size distribution index                                   []
      & bclapp,               & !< Exponent B in Clapp and Hornberger
      & matrix_pot,           & !< Matrix potential [m]
      & hyd_cond_sat            !< Saturated hydraulic conductivity [m/s]

    TYPE(t_jsb_var_real2d) :: &
      & fract_water,          & !< Wet (skin reservoir) fraction of tile (soil and canopy)                  []
      & w_skin,               & !< Water content in skin reservoir (soil and canopy)                        [m]
      & snow_accum,           & !< Snow accumulation at non-glacier/non-lake land points                    [m water equivalent]
      & evapotrans_soil,      & !< Evapotranspiration from soil w/o snow and skin evaporation               [kg m-2 s-1]
      !& evapo_soil,           & !< Evaporation from ground                                                  [kg m-2 s-1]
      & evapo_skin,           & !< Evaporation from skin reservoir                                          [kg m-2 s-1]
      & evapo_deficit,        & !< Evaporation deficit flux due to inconsistent treatment of snow evap.     [m water equivalent]
      & water_excess,         & !< Water available for infiltration into the soil                           [m water equivalent]
      & w_soil_overflow,      & !< Water content of reservoir for w_soil_overflown soil water (terraplanet) [m]
      & w_soil_column,        & !< Water content in the whole soil column (root zone)                       [m]
      & w_soil_rel              !< Filling of the soil water bucket relative to MaxMoisture (=field capacity?) TODO

    TYPE(t_jsb_var_real3d) :: &
      & vol_field_cap_sl,     & !< Volumetric soil field capacity                                           [m/m]
      & vol_p_wilt_sl,        & !< Volumetric permanent wilting point                                       [m/m]
      & vol_porosity_sl,      & !< Volumetric porosity of mineral                                           [m/m]
      & pore_size_index_sl,   & !< Soil pore size distribution index                                        []
      & bclapp_sl,            & !< Exponent B in Clapp and Hornberger
      & matrix_pot_sl,        & !< Matrix potential                                                         [m]
      & hyd_cond_sat_sl,      & !< Saturated hydraulic conductivity                                         [m/s]
      & w_soil_sl,            & !< Water content in soil layers                                             [m]
      & w_ice_sl,             & !< Ice content in soil layers                                               [m]
      & w_soil_freeze_sl,     & !< Flux from freezing water in soil layers                                  [kg m-2 s-1]
      & w_ice_melt_sl,        & !< Flux from melting ice in soil layers                                     [kg m-2 s-1]
      & w_soil_sat_sl,        & !< Water content of soil layers at saturation (derived from soil porosity)  [m]
      & w_soil_fc_sl,         & !< Water content of soil layers at field capacity
                                !! (derived from volumetric field capacity)                                 [m]
      & w_soil_pwp_sl           !< Water content of soil layers at permanent wilting point
                                !! (derived from vol. perm wilt point)                                      [m]

    ! Additional variables for PFT lct_type

    ! Parameters
    TYPE(t_jsb_var_real3d) :: &
      & root_depth_sl           !< Rooted depth per soil layer (until rooting depth) [m]
    TYPE(t_jsb_var_real2d) :: &
      & root_depth,           & !< Rooting depth [m]
      & fract_snow_can,       & !< Snow fraction on canopy                                                  []
      & w_snow_can,           & !< Water content of snow reservoir on canopy                                [m water equivalent]
      & w_soil_root,          & !< Water content in root zone of the soil                                   [m]
      & w_soil_root_fc,       & !< Water content at field capacity in root zone of the soil                 [m]
      & w_soil_root_pwp,      & !< Water content at permanent wilting point in root zone of the soil        [m]
      & water_stress,         & !< Water stress factor of canopy                        (1: no water stress, 0: infinite stress)
      & canopy_cond_unlimited,& !< Canopy conductance without water limit  [m/s] (mean value)
      & canopy_cond_limited,  & !< Canopy conductance with water limit [m/s] (mean value)
      & transpiration           !< Transpiration

    ! Additional variables for GLACIER lct_type
    TYPE(t_jsb_var_real2d) :: &
      & fract_snow_glac,      & !< Snow fraction on glacier                                                 []
      & w_glac,               & !< Glacier depth including snow on glacier                                  [m water equivalent]
      & runoff_glac             !< Runoff from glacier (rain+melt, no calving)                              [m water equivalent]

    ! Additional variables if no separate glacier lct_type is used and glaciers are treated as part of SOIL lct_type
    TYPE(t_jsb_var_real2d) :: &
      & fract_glac          !< Glacier fraction of tile                                                 []

    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & evapo_wtr,            & !< Evaporation over lake water
      & evapo_ice,            & !< Evaporation over lake ice
      & fract_snow_lice,      & !< Fraction of swow on lake ice                                             []
      & w_snow_lice             !< Water content of snow on lake ice                                        [m water equivalent]

    !! quincy
    !Soil Water addons for QUINCY in HYDRO
    TYPE(t_jsb_var_real2d) :: &
      & w_soil_root_theta,    & !< root zone integrated plant available water, fraction of maximum 
      & w_soil_root_pot         !< soil water potential [MPa]

    ! Diagnostic global land means/sums e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) :: &
      trans_gmean,            & !< Global mean transpiration                [kg m-2 s-1]
      evapotrans_gmean,       & !< Global land mean evapotranspiration      [kg m-2 s-1]
      water_content_gsum,     & !< Global land water content                [km3]
      discharge_ocean_gsum,   & !< Global water discharge to the oceans     [Sv]
      w_soil_rel_gmean,       & !< Global mean relative soil moisture       [1]
      fract_snow_gsum,        & !< Global snow area on non-glacier land     [Mio m2]
      w_snow_gsum               !< Global snow amount on non-glacier land   [Gt]

  CONTAINS
    PROCEDURE :: Init => Init_hydro_memory
  END TYPE t_hydro_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_memory_class'

CONTAINS

  SUBROUTINE Init_hydro_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    & TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    USE mo_jsb_physical_constants, ONLY: &
      & dens_snow,                       & ! Density of snow (l_dynsnow=.FALSE.)
      & dens_snow_min                      ! Density of fresh snow (l_dynsnow=.TRUE.)

    CLASS(t_hydro_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)    :: prefix
    CHARACTER(len=*),      INTENT(in)    :: suffix
    INTEGER,               INTENT(in)    :: lct_ids(:)
    INTEGER,               INTENT(in)    :: model_id

    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_config(ASSIMI_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid            ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_w  ! Vertical grids
    INTEGER :: table
    TYPE(t_grib2) :: grib2_desc

    REAL(wp) :: ini_snow_dens

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hydro_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')
    soil_w       => Get_vgrid('soil_depth_water')

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_config(ASSIMI_)

    IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
      ini_snow_dens = dens_snow_min
    ELSE
      ini_snow_dens = dens_snow
    END IF

    ! Common variables

    IF (TRIM(mem%owner_tile_name) == 'box') THEN
      grib2_desc = t_grib2(1,0,202, grib_bits)
    ELSE
      grib2_desc = t_grib2(255, 255, 255, grib_bits)
    END IF
    CALL mem%Add_var('fract_snow', mem%fract_snow,                             &
      & hgrid, surface,                                                        &
      & t_cf('fraction of snow on surface', '-', ''),                          &
      & t_grib1(table, 255, grib_bits), grib2_desc,                            &
      & prefix, suffix,                                                        &
      & output_level=BASIC,                                                    &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    IF (TRIM(mem%owner_tile_name) == 'box') THEN
      grib2_desc = t_grib2(1,0,212, grib_bits)
    ELSE
      grib2_desc = t_grib2(255, 255, 255, grib_bits)
    END IF
    CALL mem%Add_var( 'w_snow', mem%w_snow,                                           &
      & hgrid, surface,                                                               &
      & t_cf('Water content of snow reservoir on surface', 'm water equivalent', ''), &
      & t_grib1(table, 255, grib_bits), grib2_desc,                                   &
      & prefix, suffix,                                                               &
      & output_level=BASIC, &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! Total evapotranspiration, including sublimation
    CALL mem%Add_var( 'evapotrans', mem%evapotrans,                                          &
      & hgrid, surface,                                                                      &
      & t_cf('surface_evapotranspiration', 'kg m-2 s-1', 'Evapotranspiration from surface'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
      & prefix, suffix,                                                                      &
      & lrestart=.TRUE.,                                                                     &
      & output_level=BASIC,                                                                  &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! Potential evapotranspiration/sublimation
    CALL mem%Add_var( 'evapopot', mem%evapopot,                                                    &
      & hgrid, surface,                                                                            &
      & t_cf('surface_evaporation_potential', 'kg m-2 s-1', 'Potential evaporation from surface'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                         &
      & prefix, suffix,                                                                            &
      & lrestart=.TRUE.,                                                                           &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'evapo_snow', mem%evapo_snow,                            &
      & hgrid, surface,                                                        &
      & t_cf('evapo_snow', 'kg m-2 s-1', 'Evaporation from snow'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.FALSE.,                                                      &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'snowmelt', mem%snowmelt,                                &
      & hgrid, surface,                                                        &
      & t_cf('snow_melt', 'kg m-2 s-1', 'Snow melt on surface'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.FALSE.,                                                      &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'runoff', mem%runoff,                                   &
      & hgrid, surface,                                                       &
      & t_cf('surface_runoff', 'kg m-2 s-1', 'Surface runoff'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & lrestart=.FALSE.,                                                     &
      & loutput=.TRUE., output_level=BASIC,                                   &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'drainage', mem%drainage,                              &
      & hgrid, surface,                                                      &
      & t_cf('drainage', 'kg m-2 s-1', 'drainage'),                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & loutput=.TRUE., output_level=BASIC,                                  &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'ws_negative', mem%ws_negative,                        &
      & hgrid, surface,                                                      &
      & t_cf('ws_negative', 'kg m-2', 'negative soil water'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & loutput=.TRUE., output_level=BASIC,                                  &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'discharge', mem%discharge,                            &
      & hgrid, surface,                                                      &
      & t_cf('discharge', 'kg m-2 s-1', 'local discharge'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'discharge_ocean', mem%discharge_ocean,                &
      & hgrid, surface,                                                      &
      & t_cf('discharge_ocean', 'm3 s-1', 'discharge to the ocean'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'internal_drain', mem%internal_drain,                &
      & hgrid, surface,                                                    &
      & t_cf('internal_drain', 'kg m-2 s-1', 'internal drainage'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,                                                    &
      & lrestart=.FALSE.,                                                  &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'water_flux', mem%water_flux,                                  &
      & hgrid, surface,                                                              &
      & t_cf('water_flux', 'm3 s-1', 'total water flux'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_content', mem%water_content,                            &
      & hgrid, surface,                                                              &
      & t_cf('water_content', 'm3', 'total water content'),                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_budget', mem%water_budget,                              &
      & hgrid, surface,                                                              &
      & t_cf('water_budget', 'm3', 'water budget error'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp )

    IF (jsbach_runs_standalone()) THEN
       CALL mem%Add_var( 'elevation', mem%elevation,                                    &
         & hgrid, surface,                                                              &
         & t_cf('elevation', '', ''),                                                   &
         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
         & prefix, suffix,                                                              &
         & lrestart=.FALSE.,                                                            &
         & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )
    END IF

    CALL mem%Add_var( 'oro_stddev', mem%oro_stddev,                                  &
      & hgrid, surface,                                                              &
      & t_cf('orography_standard_deviation', '', ''),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    IF (model%config%use_tmx) THEN
      CALL mem%Add_var( 'q_snocpymlt', mem%q_snocpymlt,                          &
        & hgrid, surface,                                                        &
        & t_cf('heating_snow_cpy_melt', 'W m-2', ''),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE.,                                                      &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    END IF

      ! Additional variables for land type
    IF ( (     One_of(VEG_TYPE,     lct_ids(:)) > 0 &
      &   .OR. One_of(BARE_TYPE,    lct_ids(:)) > 0 &
      &   .OR. One_of(GLACIER_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE,    lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'fract_snow_soil', mem%fract_snow_soil,                 &
        & hgrid, surface,                                                       &
        & t_cf('fraction of snow on soil', '-', ''),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
        & prefix, suffix,                                                       &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_snow_soil', mem%w_snow_soil,                              &
        & hgrid, surface,                                                            &
        & t_cf('Water content of snow reservoir on soil', 'm water equivalent', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! CALL mem%Add_var( 'w_snow_soil_old', mem%w_snow_soil_old,                            &
      !   & hgrid, surface,                                                                  &
      !   & t_cf('Water content of snow reservoir on soil (old)', 'm water equivalent', ''), &
      !   & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      !   & prefix, suffix,                                                                  &
      !   & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'snow_soil_dens', mem%snow_soil_dens,                &
        & hgrid, surface,                                                    &
        & t_cf('snow_soil_dens', 'kg m-3', 'Density of snow on soil'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & initval_r=ini_snow_dens, l_aggregate_all=.TRUE. )

      ! Total evapotranspiration, including sublimation
      CALL mem%Add_var( 'evapotrans_lnd', mem%evapotrans_lnd,                                                     &
        & hgrid, surface,                                                                                         &
        & t_cf('surface_evapotranspiration_lnd', 'kg m-2 s-1', 'Evapotranspiration from land surface w/o lakes'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & lrestart=.FALSE.,                                                                                       &
        & loutput=.TRUE., output_level=BASIC,                                                                     &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (.NOT. model%config%use_tmx) THEN
        CALL mem%Add_var( 'q_snocpymlt', mem%q_snocpymlt,                          &
          & hgrid, surface,                                                        &
          & t_cf('heating_snow_cpy_melt', 'W m-2', ''),                            &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
          & prefix, suffix,                                                        &
          & lrestart=.FALSE.,                                                      &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
      END IF

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(1,0,201, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'fract_water', mem%fract_water,                          &
        & hgrid, surface,                                                        &
        & t_cf('fraction of water on surface surface_wet_fraction', '-', ''),    &
        & t_grib1(table, 255, grib_bits), grib2_desc,                            &
        & prefix, suffix,                                                        &
        & output_level=BASIC,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    END IF

    ! Additional variables for soil
    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(BARE_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'soil_depth', mem%soil_depth,                          &
        & hgrid, surface,                                                      &
        & t_cf('soil_depth', 'm', 'Depth of soil until bedrock'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE.,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'soil_depth_sl', mem%soil_depth_sl,                    &
        & hgrid, soil_w,                                                       &
        & t_cf('soil_depth_layer', 'm', ''),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput=.TRUE.,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
        CALL mem%Add_var( 'fract_org_sl', mem%fract_org_sl,                         &
          & hgrid, soil_w,                                                          &
          & t_cf('fract_org_sl', '', 'Fraction of organic material in soil layer'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
          & prefix, suffix,                                                         &
          & lrestart=.TRUE.,                                                        &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
      END IF

      CALL mem%Add_var( 'max_moist', mem%max_moist,                            &
        & hgrid, surface,                                                      &
        & t_cf('max_moist', 'm', 'Maximum root zone moisture'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_field_cap', mem%vol_field_cap,                    &
        & hgrid, surface,                                                      &
        & t_cf('vol_field_capacity', 'm/m', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_field_cap_sl', mem%vol_field_cap_sl,              &
        & hgrid, soil_w,                                                       &
        & t_cf('vol_field_capacity_sl', 'm/m', ''),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'vol_p_wilt', mem%vol_p_wilt,                          &
        & hgrid, surface,                                                      &
        & t_cf('vol_p_wilt', 'm/m', ''),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_p_wilt_sl', mem%vol_p_wilt_sl,                    &
        & hgrid, soil_w,                                                       &
        & t_cf('vol_p_wilt_sl', 'm/m', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'vol_porosity', mem%vol_porosity,                      &
        & hgrid, surface,                                                      &
        & t_cf('vol_porosity', '', ''),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'vol_porosity_sl', mem%vol_porosity_sl,                &
        & hgrid, soil_w,                                                       &
        & t_cf('vol_porosity_sl', 'm/m', ''),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'pore_size_index', mem%pore_size_index,                &
        & hgrid, surface,                                                      &
        & t_cf('pore_size_index', '', ''),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'pore_size_index_sl', mem%pore_size_index_sl,          &
        & hgrid, soil_w,                                                       &
        & t_cf('pore_size_index_sl', '', ''),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'bclapp', mem%bclapp,                                  &
        & hgrid, surface,                                                      &
        & t_cf('bclapp_sl', '', ''),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'bclapp_sl', mem%bclapp_sl,                            &
        & hgrid, soil_w,                                                       &
        & t_cf('bclapp_sl', '', ''),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'matrix_pot', mem%matrix_pot,                          &
        & hgrid, surface,                                                      &
        & t_cf('matrix_potential', 'm', ''),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'matrix_pot_sl', mem%matrix_pot_sl,                    &
        & hgrid, soil_w,                                                       &
        & t_cf('matrix_potential_sl', 'm', ''),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'hyd_cond_sat', mem%hyd_cond_sat,                      &
        & hgrid, surface,                                                      &
        & t_cf('sat_hydraulic_cond', 'm s-1', ''),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'hyd_cond_sat_sl', mem%hyd_cond_sat_sl,                &
        & hgrid, soil_w,                                                       &
        & t_cf('sat_hydraulic_cond_sl', 'm s-1', ''),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(1,0,211, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'w_skin', mem%w_skin,                                                         &
        & hgrid, surface,                                                                             &
        & t_cf('skin_reservoir', 'm water equivalent', 'Water content in skin reservoir of surface'), &
        & t_grib1(table, 255, grib_bits), grib2_desc,                                                 &
        & prefix, suffix,                                                                             &
        & output_level=BASIC,                                                                         &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'snow_accum', mem%snow_accum,                                                        &
        & hgrid, surface,                                                                                    &
        & t_cf('snow_accum', 'm water equivalent', 'Snow accumulation at non-glacier/non-lake land points'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
        & prefix, suffix,                                                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    ! END IF

    ! IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
    !   &   .OR. One_of(BARE_TYPE, lct_ids(:)) > 0 &
    !   &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
    !   &  )                                       &
    !   &  .AND. .NOT. ASSOCIATED(mem%fract_water%ptr) ) THEN

    CALL mem%Add_var( 'evapotrans_soil', mem%evapotrans_soil,                                            &
        & hgrid, surface,                                                                                &
        & t_cf('evapotrans_soil', 'kg m-2 s-1', 'Evapotranspiration from soil (w/o snow and skin res.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
        & prefix, suffix,                                                                                &
        & loutput=.TRUE.,                                                                                &
        & lrestart=.FALSE.,                                                                              &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'evapo_skin', mem%evapo_skin,                               &
        & hgrid, surface,                                                         &
        & t_cf('evapo_skin', 'kg m-2 s-1', 'Evaporation from skin reservoir'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
        & prefix, suffix,                                                         &
        & lrestart=.FALSE.,                                                       &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'evapo_deficit', mem%evapo_deficit,                      &
        & hgrid, surface,                                                      &
        & t_cf('evaporation_deficit', 'm water equivalent', ''),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'water_excess', mem%water_excess,                        &
        & hgrid, surface,                                                      &
        & t_cf('water_excess', 'm water equivalent', ''),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_overflow', mem%w_soil_overflow,                        &
        & hgrid, surface,                                                              &
        & t_cf('Water content of reservoir for soil water', 'm water equivalent', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(1,0,213, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'w_soil_column', mem%w_soil_column,                          &
        & hgrid, surface,                                                            &
        & t_cf('Water content in the whole soil column', 'm water equivalent', ''),  &
        & t_grib1(table, 255, grib_bits), grib2_desc,                                &
        & prefix, suffix,                                                            &
        & loutput=.TRUE., output_level=BASIC,                                        &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_sl', mem%w_soil_sl,                          &
        & hgrid, soil_w,                                                     &
        & t_cf('Water content in soil layers', 'm water equivalent', ''),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_ice_sl', mem%w_ice_sl,                            &
        & hgrid, soil_w,                                                     &
        & t_cf('Ice content in soil layers', 'm water equivalent', ''),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_freeze_sl', mem%w_soil_freeze_sl,            &
        & hgrid, soil_w,                                                     &
        & t_cf('Freezing water flux in soil layers', 'kg m-2 s-1', ''),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE.,                                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_ice_melt_sl', mem%w_ice_melt_sl,                  &
        & hgrid, soil_w,                                                     &
        & t_cf('Melting ice flux in soil layers', 'kg m-2 s-1', ''),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE.,                                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_sat_sl', mem%w_soil_sat_sl,                             &
        & hgrid, soil_w,                                                                &
        & t_cf('Water content in soil layers at saturation', 'm water equivalent', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
        & prefix, suffix,                                                               &
        & loutput=.TRUE., &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_fc_sl', mem%w_soil_fc_sl,                                   &
        & hgrid, soil_w,                                                                    &
        & t_cf('Water content in soil layers at field capacity', 'm water equivalent', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & loutput=.TRUE., &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_pwp_sl', mem%w_soil_pwp_sl,                                          &
        & hgrid, soil_w,                                                                             &
        & t_cf('Water content in soil layers at permanent wilting point', 'm water equivalent', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                         &
        & prefix, suffix,                                                                            &
        & loutput=.TRUE., &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var('w_soil_rel', mem%w_soil_rel,                                &
        &  hgrid, surface,                                                          &
        &  t_cf('water_soil_relative', '', ''),                                     &
        &  t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp )

    END IF

    ! Additional variables for vegetation
    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'root_depth', mem%root_depth,                        &
        & hgrid, surface,                                                    &
        & t_cf('rooting_depth', '', ''),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'root_depth_sl', mem%root_depth_sl,                  &
        & hgrid, soil_w,                                                     &
        & t_cf('root_depth_layer', 'm', ''),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL mem%Add_var( 'fract_snow_can', mem%fract_snow_can,                &
        & hgrid, surface,                                                    &
        & t_cf('fraction of snow on canopy', '-', ''),                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & output_level=BASIC,                                                &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_snow_can', mem%w_snow_can,                                  &
        & hgrid, surface,                                                              &
        & t_cf('Water content of snow reservoir on canopy', 'm water equivalent', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC,                                                          &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_root', mem%w_soil_root,                              &
        & hgrid, surface,                                                            &
        & t_cf('Water content in root zone of the soil', 'm water equivalent', ''),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & lrestart=.FALSE.,                                                          &
        & loutput=.TRUE., &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_root_fc', mem%w_soil_root_fc,                                         &
        & hgrid, surface,                                                                             &
        & t_cf('Water content at field capacity in root zone of the soil', 'm water equivalent', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
        & prefix, suffix,                                                                             &
        & lrestart=.FALSE.,                                                                           &
        & loutput=.TRUE., &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_soil_root_pwp', mem%w_soil_root_pwp,                                                &
        & hgrid, surface,                                                                                      &
        & t_cf('Water content at permanent wilting point in root zone of the soil', 'm water equivalent', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                   &
        & prefix, suffix,                                                                                      &
        & lrestart=.FALSE.,                                                                                    &
        & loutput=.TRUE., &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'water_stress', mem%water_stress,                     &
        & hgrid, surface,                                                     &
        & t_cf('Water stress factor of canopy', '-', ''),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & loutput=.TRUE., &
        & lrestart=.FALSE., initval_r=0.0_wp )

      IF (.NOT. dsl4jsb_Config(ASSIMI_)%active) THEN
        CALL mem%Add_var( 'canopy_cond_unlimited', mem%canopy_cond_unlimited,  &
          & hgrid, surface,                                                    &
          & t_cf('canopy_conductance_unlimited', '', ''),                      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & loutput=.TRUE., &
          & initval_r=0.0_wp )

        CALL mem%Add_var( 'canopy_cond_limited', mem%canopy_cond_limited,      &
          & hgrid, surface,                                                    &
          & t_cf('canopy_conductance_limited', '', ''),                        &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          & loutput=.TRUE., &
          & initval_r=0.0_wp )
      END IF

      CALL mem%Add_var( 'transpiration', mem%transpiration,                          &
        & hgrid, surface,                                                            &
        & t_cf('surface_transpiration', 'kg m-2 s-1', 'Transpiration from surface'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & lrestart=.FALSE.,                                                          &
        & loutput=.TRUE., &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    END IF

    ! Additional variables for GLACIER lct
    IF ( (     One_of(GLACIER_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE,    lct_ids(:)) > 0 &
      &  ) ) THEN

      CALL mem%Add_var( 'fract_snow_glac', mem%fract_snow_glac,              &
        & hgrid, surface,                                                    &
        & t_cf('fraction of snow on glacier', '-', ''),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & output_level=BASIC,                                                &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'w_glac', mem%w_glac,                                                  &
        & hgrid, surface,                                                                      &
        & t_cf('water_glac', 'm water equivalent', 'Water content of glacier including snow'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix,                                                                      &
        & output_level=BASIC,                                                                  &
        & initval_r=20.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'runoff_glac', mem%runoff_glac,                            &
        & hgrid, surface,                                                          &
        & t_cf('runoff_glac', 'kg m-2 s-1', 'Runoff from glacier (rain + melt)'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
        & prefix, suffix,                                                          &
        & initval_r=20.0_wp, l_aggregate_all=.TRUE. )

    END IF

    ! Additional variables for LAKE lct
    IF (      One_of(LAKE_TYPE, lct_ids(:)) > 0 ) THEN

      CALL mem%Add_var( 'evapo_wtr', mem%evapo_wtr,                                                   &
        & hgrid, surface,                                                                             &
        & t_cf('surface_evapotranspiration_wtr', 'kg m-2 s-1', 'Evapotranspiration from lake water'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
        & prefix, suffix,                                                                             &
        & lrestart=.FALSE.,                                                                           &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        CALL mem%Add_var( 'evapo_ice', mem%evapo_ice,                                                 &
          & hgrid, surface,                                                                           &
          & t_cf('surface_evapotranspiration_ice', 'kg m-2 s-1', 'Evapotranspiration from lake ice'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
          & prefix, suffix,                                                                           &
          & lrestart=.FALSE.,                                                                         &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        !Fraction of snow on lake ice, rel. to ice fraction
        CALL mem%Add_var( 'fract_snow_lice', mem%fract_snow_lice,                                     &
          & hgrid, surface,                                                                           &
          & t_cf('fract_snow_lice', '-', 'fraction of snow on lake ice (rel. to ice fraction)'),      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
          & prefix, suffix,                                                                           &
          & output_level=BASIC,                                                                       &
          & initval_r=0.0_wp )

        CALL mem%Add_var( 'w_snow_lice', mem%w_snow_lice,                                             &
          & hgrid, surface,                                                                           &
          & t_cf('w_snow_lice', 'm water equivalent', 'Water content of snow reservoir on lake ice'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
          & prefix, suffix,                                                                           &
          & output_level=BASIC,                                                                       &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
      END IF

    END IF

    !! quincy
    IF (model%config%use_quincy) THEN
      CALL mem%Add_var('w_soil_root_theta', mem%w_soil_root_theta,                                           &
        & hgrid, surface,                                                                                    &
        & t_cf('w_soil_root_theta', '', 'root zone integrated plant available water, fraction of maximum '), &
        & t_grib1(table, 255, grib_bits),                                                                    &
        & t_grib2(255, 255, 255, grib_bits),                                                                 &
        & prefix, suffix,                                                                                    &
        & loutput=.TRUE., lrestart=.TRUE.,                                                                   &
        & initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('w_soil_root_pot', mem%w_soil_root_pot,                                               &
        & hgrid, surface,                                                                                    &
        & t_cf('w_soil_root_pot', 'MPa', 'soil water potential'),                                            &
        & t_grib1(table, 255, grib_bits),                                                                    &
        & t_grib2(255, 255, 255, grib_bits),                                                                 &
        & prefix, suffix,                                                                                    &
        & loutput=.TRUE., lrestart=.TRUE.,                                                                   &
        & initval_r=0.0_wp) ! initvalue
    ENDIF

#ifdef __ICON__
    ! Diagnostic 1d global land variables for experiment monitoring
    !       (1d stream variables are not supported with echam)
    !
    IF ( TRIM(suffix) == 'box' ) THEN
      CALL mem%Add_var('trans_gmean', mem%trans_gmean,                                   &
        & hgrid, surface,                                                                &
        & t_cf('trans_gmean', 'kg m-2 s-1', 'Global mean transpiration'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('evapotrans_gmean', mem%evapotrans_gmean,                         &
        & hgrid, surface,                                                                &
        & t_cf('evapotrans_gmean', 'kg m-2 s-1', 'Global land evapotranspiration'),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('water_content_gsum', mem%water_content_gsum,                     &
        & hgrid, surface,                                                                &
        & t_cf('water_content_gsum', 'km3', 'Global land water content'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('discharge_ocean_gsum', mem%discharge_ocean_gsum,                 &
        & hgrid, surface,                                                                &
        & t_cf('discharge_ocean_gsum', 'Sv', 'Global water discharge to the oceans'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('w_soil_rel_gmean', mem%w_soil_rel_gmean,                         &
        & hgrid, surface,                                                                &
        & t_cf('w_soil_rel_gmean', '1', 'Global mean relative soil moisture'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('fract_snow_gsum', mem%fract_snow_gsum,                           &
        & hgrid, surface,                                                                &
        & t_cf('fract_snow_gsum', 'Mio km2', 'Global snow area (including glaciers)'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
      CALL mem%Add_var('w_snow_gsum', mem%w_snow_gsum,                                   &
        & hgrid, surface,                                                                &
        & t_cf('w_snow_gsum', 'Gt', 'Global snow amount on non-glacier land'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )
    END IF
#endif

    END SUBROUTINE Init_hydro_memory
#endif
END MODULE mo_hydro_memory_class
