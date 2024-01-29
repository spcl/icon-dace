!> Contains types for the atmosphere-land interface memory.
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
MODULE mo_a2l_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  !USE mo_exception,         ONLY: message
  USE mo_util,              ONLY: One_of

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d
  USE mo_jsb_lct_class,     ONLY: LAKE_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_a2l_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 35

  TYPE, EXTENDS(t_jsb_memory) :: t_a2l_memory
    TYPE(t_jsb_var_real2d) :: &
      DEBUG_VAR          , &
      t_air              , &
      q_air              , &
      press_air          , &
      rain               , &
      snow               , &
      wind_air           , &
      wind_10m           , &
      lw_srf_down        , &
      swvis_srf_down     , &
      swnir_srf_down     , &
      swpar_srf_down     , &
      fract_par_diffuse  , &
      dz_srf             , &
      press_srf          , &
      rho_srf            , &
      drag_srf           , &
      t_acoef            , &
      t_bcoef            , &
      q_acoef            , &
      q_bcoef            , &
      pch                , &
      cos_zenith_angle   , &
      CO2_air            , &
      CO2_air_mol        , &
      ! For lakes:
      drag_wtr           , &
      drag_ice           , &
      t_acoef_wtr        , &
      t_bcoef_wtr        , &
      q_acoef_wtr        , &
      q_bcoef_wtr        , &
      t_acoef_ice        , &
      t_bcoef_ice        , &
      q_acoef_ice        , &
      q_bcoef_ice
    
    !! quincy
    !TYPE(t_jsb_var_real2d) :: &
      ! & t_air                 , &  !> 2m air temperature [K] 
      ! & q_air                 , &  !> 2m air specific humidity [g g-1] 
      ! & rain                  , &  !> liquid water precipitation [kg m-2 s-1]
      ! & snow                  , &  !> solid water precipitation [kg m-2 s-1]
      ! & wind_air              , &  !> 2m air wind speed [m s-1] 
      ! & wind_10m              , &  !> 10m air wind speed [m s-1] 
      ! & lw_srf_down           , &  !> longwave surface downward flux [W m-2] 
      ! & swvis_srf_down        , &  !> shortwave surface downward flux in the visible band [W m-2]
      ! & swnir_srf_down        , &  !> shortwave surface downward flux in the near infrared band [W m-2] 
      ! & swpar_srf_down        , &  !> shortwave surface downward flux in the PAR band [W m-2] 
      ! & swpar_fdiffuse        , &  !> fraction of shortwave surface downward flux that is diffuse [fraction]
      ! & press_srf             , &  !> surface air pressure [Pa]
      ! & drag_srf              , &  !> total surface drag [unitless] 
      ! & t_acoef               , &  !> Richtmeyer-Morten coefficient A for temperature [unitless] 
      ! & t_bcoef               , &  !> Richtmeyer-Morten coefficient B for temperature [unitless]
      ! & q_acoef               , &  !> Richtmeyer-Morten coefficient A for specific humidity [unitless] 
      ! & q_bcoef               , &  !> Richtmeyer-Morten coefficient B for specific humidity [unitless]
      ! & pch                   , &  !> total surface drag [unitless] 
      ! & cos_zenith_angle           !> Cosine of solar zenith angle [unitless]

    ! CO2
    TYPE(t_jsb_var_real2d) :: &
      ! jsbach4 variables, in quincy routines replaced by co2_mixing_ratio  _c13  _c14
      ! & CO2_air               , &  !> molar ratio (mass)
      ! & CO2_air_mol           , &  !> molar ratio (volume)
      & CO2_mixing_ratio      , &  !< surface air CO2 mixing ratio [ppmv]   >> mo_atmland_interface: CO2_mixing_ratio = CO2_air_mol * 1000000._wp
      & CO2_mixing_ratio_C13  , &  !< surface air 13CO2 mixing ratio [ppmv]
      & CO2_mixing_ratio_C14       !< surface air 14CO2 mixing ratio [scaled]
    
    
  CONTAINS
    PROCEDURE :: Init => Init_a2l_memory
  END TYPE t_a2l_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_a2l_memory_class'

CONTAINS

  SUBROUTINE Init_a2l_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_a2l_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_grid),    POINTER :: hgrid
    TYPE(t_jsb_vgrid),   POINTER :: surface
    TYPE(t_jsb_model),   POINTER :: model
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_a2l_memory'
    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used
    table = tables(1)

    model   => Get_model(model_id)

    hgrid => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL mem%Add_var( 'DEBUG_VAR', mem%DEBUG_VAR,                            &
      & hgrid, surface,                                                      &
      & t_cf('DEBUG_VAR', '', ''),                                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=777.7_wp )

    CALL mem%Add_var( 't_air', mem%t_air,                                    &
      & hgrid, surface,                                                      &
      & t_cf('air_temperature', 'K', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=283.15_wp )  ! R: I just invented a value to start with...

    CALL mem%Add_var( 'q_air', mem%q_air,                                                 &
      & hgrid, surface,                                                                   &
      & t_cf('air_specific_humidity', 'kg kg-1', 'Specific humidity of air at surface.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.FALSE.,                                                                 &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'rain', mem%rain,                                      &
      & hgrid, surface,                                                      &
      & t_cf('rain_fall', 'kg/m2/s', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'snow', mem%snow,                                      &
      & hgrid, surface,                                                      &
      & t_cf('snow_fall', 'kg/m2/s', ''),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'wind_air', mem%wind_air,                              &
      & hgrid, surface,                                                      &
      & t_cf('air_wind_speed', 'm/s', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'wind_10m', mem%wind_10m,                              &
      & hgrid, surface,                                                      &
      & t_cf('10m_wind_speed', 'm/s', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'lw_srf_down', mem%lw_srf_down,                        &
      & hgrid, surface,                                                      &
      & t_cf('downward_longwave_radiation', 'W/m2', ''),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=300.0_wp )

    CALL mem%Add_var( 'swvis_srf_down', mem%swvis_srf_down,                  &
      & hgrid, surface,                                                      &
      & t_cf('downward_visible_shortwave_radiation', 'W/m2', ''),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=100.0_wp )

    CALL mem%Add_var( 'swnir_srf_down', mem%swnir_srf_down,                  &
      & hgrid, surface,                                                      &
      & t_cf('downward_nir_shortwave_radiation', 'W/m2', ''),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=90.0_wp )

    CALL mem%Add_var( 'swpar_srf_down', mem%swpar_srf_down,                  &
      & hgrid, surface,                                                      &
      & t_cf('downward_par_shortwave_radiation', 'W/m2', ''),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'fract_par_diffuse', mem%fract_par_diffuse,            &
      & hgrid, surface,                                                      &
      & t_cf('fract_par_diffuse', '', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'press_srf', mem%press_srf,                            &
      & hgrid, surface,                                                      &
      & t_cf('surface_pressure', 'N/m2', ''),                                & ! N/m2 =Pa
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=101325.0_wp )  ! R: I just took average atmospheric pressure
                                 !    at sea level as value to start with...

    IF (model%config%use_tmx) THEN
      CALL mem%Add_var( 'press_air', mem%press_air,                                         &
        & hgrid, surface,                                                                   &
        & t_cf('air_pressure', 'N m-2', 'Pressure at lowest level'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & lrestart=.FALSE.,                                                                 &
        & output_level=BASIC, &
        & initval_r=98300.0_wp )

      CALL mem%Add_var( 'dz_srf', mem%dz_srf,                                  &
        & hgrid, surface,                                                      &
        & t_cf('reference height in surface layer', 'm', ''),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=FULL,                                                   &
        & initval_r=10.0_wp )

      CALL mem%Add_var( 'rho_srf', mem%rho_srf,                                &
        & hgrid, surface,                                                      &
        & t_cf('air density at surface', '', ''),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=1.0_wp )
    ELSE
      CALL mem%Add_var( 'drag_srf', mem%drag_srf,                              &
        & hgrid, surface,                                                      &
        & t_cf('surface_drag', '', ''),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )
      CALL mem%Add_var( 'pch', mem%pch,                                        &
        & hgrid, surface,                                                      &
        & t_cf('pch', '', ''),                                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )  
      END IF

    CALL mem%Add_var( 't_acoef', mem%t_acoef,                                &
      & hgrid, surface,                                                      &
      & t_cf('temperature_acoef', '', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 't_bcoef', mem%t_bcoef,                                &
      & hgrid, surface,                                                      &
      & t_cf('temperature_acoef', '', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'q_acoef', mem%q_acoef,                                &
      & hgrid, surface,                                                      &
      & t_cf('specific_humidity_acoef', '', ''),                             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'q_bcoef', mem%q_bcoef,                                &
      & hgrid, surface,                                                      &
      & t_cf('specific_humidity_bcoef', '', ''),                             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'cos_zenith_angle', mem%cos_zenith_angle,              &
      & hgrid, surface,                                                      &
      & t_cf('cosine_zenith_angle', '', ''),                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'co2_air', mem%CO2_air,                                                      &
      & hgrid, surface,                                                                            &
      & t_cf('CO2_air', '[kg(CO2)/kg(air)]', 'CO2 mass mixing ratio of lowest atmosphere level'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                         &
      & prefix, suffix,                                                                            &
      & output_level=BASIC, &
      & lrestart=.FALSE., initval_r=6.077E-4_wp ) ! R: I just freely choosed this value. Corresponds to 400 ppm:
                                                  !    CO2:44.011 g/mol; dry air: 28,97 g/mol
                                                  !    400 ppm => 6.077*10^-4

    CALL mem%Add_var( 'co2_air_mol', mem%CO2_air_mol,                                                  &
      & hgrid, surface,                                                                                &
      & t_cf('CO2_air_mol', '[mol(CO2)/mol(air)]', 'CO2 mol mixing ratio of lowest atmosphere level'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=BASIC, &
      & lrestart=.FALSE., initval_r=4.0E-4_wp ) ! R: I just choosed this value for 400 ppm

    IF (.NOT. model%config%use_tmx .AND. One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%t_acoef_wtr%ptr)) THEN

      CALL mem%Add_var( 'drag_wtr', mem%drag_wtr,                              &
        & hgrid, surface,                                                      &
        & t_cf('surface_drag_wtr', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

      CALL  mem%Add_var( 'drag_ice', mem%drag_ice,                             &
        & hgrid, surface,                                                      &
        & t_cf('surface_drag_ice', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 't_acoef_wtr', mem%t_acoef_wtr,                        &
        & hgrid, surface,                                                      &
        & t_cf('temperature_acoef_wtr', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 't_bcoef_wtr', mem%t_bcoef_wtr,                        &
        & hgrid, surface,                                                      &
        & t_cf('temperature_acoef_wtr', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'q_acoef_wtr', mem%q_acoef_wtr,                        &
        & hgrid, surface,                                                      &
        & t_cf('specific_humidity_acoef_wtr', '', ''),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'q_bcoef_wtr', mem%q_bcoef_wtr,                        &
        & hgrid, surface,                                                      &
        & t_cf('specific_humidity_bcoef_wtr', '', ''),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 't_acoef_ice', mem%t_acoef_ice,                        &
        & hgrid, surface,                                                      &
        & t_cf('temperature_acoef_ice', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 't_bcoef_ice', mem%t_bcoef_ice,                        &
        & hgrid, surface,                                                      &
        & t_cf('temperature_acoef_ice', '', ''),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'q_acoef_ice', mem%q_acoef_ice,                        &
        & hgrid, surface,                                                      &
        & t_cf('specific_humidity_acoef_ice', '', ''),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'q_bcoef_ice', mem%q_bcoef_ice,                        &
        & hgrid, surface,                                                      &
        & t_cf('specific_humidity_bcoef_ice', '', ''),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    &
        & output_level=BASIC,                                                  &
        & initval_r=0.0_wp )

    END IF

    IF (.NOT. model%config%use_quincy) RETURN

    !! quincy
    ! CALL mem%Add_var('swpar_fdiffuse', mem%swpar_fdiffuse, &
    !  & hgrid, surface, &
    !  & t_cf('swpar_fdiffuse', 'fraction', 'fraction of shortwave surface downward flux that is diffuse'), &
    !  & t_grib1(table, 255, grib_bits), &
    !  & t_grib2(255, 255, 255, grib_bits), &
    !  & prefix, suffix, &
    !  & loutput=.TRUE., lrestart=.TRUE., &
    !  & initval_r=0.0_wp)
    
    CALL mem%Add_var('CO2_mixing_ratio', mem%co2_mixing_ratio,                    &
      & hgrid, surface,                                                           &
      & t_cf('CO2_mixing_ratio', 'ppm', 'surface air CO2 mixing ratio'),          &
      & t_grib1(table, 255, grib_bits),                                           &
      & t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                           &
      & loutput=.TRUE., lrestart=.TRUE.,                                          &
      & initval_r=0.0_wp)

    CALL mem%Add_var('CO2_mixing_ratio_C13', mem%CO2_mixing_ratio_C13,            &
      & hgrid, surface,                                                           &
      & t_cf('CO2_mixing_ratio_C13', 'ppm', 'surface air 13CO2 mixing ratio'),    &
      & t_grib1(table, 255, grib_bits),                                           &
      & t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                           &
      & loutput=.TRUE., lrestart=.TRUE.,                                          &
      & initval_r=0.0_wp)

    CALL mem%Add_var('CO2_mixing_ratio_C14', mem%CO2_mixing_ratio_C14,            &
      & hgrid, surface,                                                           &
      & t_cf('CO2_mixing_ratio_C14', 'scaled', 'surface air 14CO2 mixing ratio'), &
      & t_grib1(table, 255, grib_bits),                                           &
      & t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                           &
      & loutput=.TRUE., lrestart=.TRUE.,                                          &
      & initval_r=0.0_wp)

  END SUBROUTINE Init_a2l_memory

#endif
END MODULE mo_a2l_memory_class
