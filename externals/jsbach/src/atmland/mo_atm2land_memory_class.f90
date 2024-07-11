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
  USE mo_util,              ONLY: One_of

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_control,       ONLY: jsbach_runs_standalone
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d
  USE mo_jsb_lct_class,     ONLY: LAKE_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_a2l_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 60

  TYPE, EXTENDS(t_jsb_memory) :: t_a2l_memory

    !------------------------------------------------------------------------------------------------------ !
    !>1.0 jsbach4
    !>

    TYPE(t_jsb_var_real2d) :: &
      & DEBUG_VAR               , &  !<
      & CO2_air                 , &  !<
      & CO2_air_mol                  !<
    TYPE(t_jsb_var_real2d) :: &
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

#ifndef __NO_QUINCY__
    !------------------------------------------------------------------------------------------------------ !
    !>3.0 quincy
    !>

    ! CO2
    TYPE(t_jsb_var_real2d) :: &
      ! jsbach4 variables, in quincy routines replaced by co2_mixing_ratio  _c13  _c14
      ! & CO2_air               , &  !< molar ratio (mass)
      ! & CO2_air_mol           , &  !< molar ratio (volume)
      & CO2_mixing_ratio      , &  !< surface air CO2 mixing ratio [ppmv]   >> mo_atmland_interface: CO2_mixing_ratio = CO2_air_mol * 1000000._wp
      & CO2_mixing_ratio_C13  , &  !< surface air 13CO2 mixing ratio [ppmv]
      & CO2_mixing_ratio_C14       !< surface air 14CO2 mixing ratio [scaled]

    ! ga   ( @TODO may be moved elsewhere in the future )
    TYPE(t_jsb_var_real2d) :: &
      & ga                         !< aerodynamic conductance [m s-1]

    ! N & P deposition
    TYPE(t_jsb_var_real2d) ::   &
      & nhx_deposition        , &  !< surface downward reduced nitrogen deposition velocity [mumol m-2 sec-1]
      & noy_deposition        , &  !< surface downward oxidised nitrogen deposition velocity [mumol m-2 sec-1]
      & nhx_n15_deposition    , &  !< surface downward reduced nitrogen-15 deposition velocity [mumol m-2 sec-1]
      & noy_n15_deposition    , &  !< surface downward oxidised nitrogen-15 deposition velocity [mumol m-2 sec-1]
      & p_deposition               !< surface downward reduced phosphoruns deposition velocity [mumol m-2 sec-1]

    ! For calculation of the timestep after local midnight (1st timestep of the new day)
    TYPE(t_jsb_var_real2d) ::   &
      & daytime_counter       , &  !< number of timesteps of current day with light available [#]
      & local_time_day_seconds     !< seconds passed on this day (local cell time)

#ifdef __QUINCY_STANDALONE__
    ! For calculation of constant input values (add values over the 1st year, calc avrg, and use this as input for teh rest of the simulation)
    TYPE(t_jsb_var_real2d) :: &
      & nhx_deposition_acc, &             !< accumulate (simply add) input value of nhx_deposition [mumol m-2 sec-1]
      & nhy_deposition_acc, &             !< accumulate (simply add) input value of nhy_deposition [mumol m-2 sec-1]
      & p_deposition_acc, &               !< accumulate (simply add) input value of p_deposition [mumol m-2 sec-1]
      & co2_mixing_ratio_acc, &           !< accumulate (simply add) input value of co2_mixing_ratio [ppm]
      & co2_dC13_acc, &                   !< accumulate (simply add) input value of co2_dC13 [ppm]
      & co2_DC14_acc, &                   !< accumulate (simply add) input value of co2_DC14 [scaled]
      & const_input_timestep_counter      !< count number of timesteps since simulation start [unitless]

    TYPE(t_jsb_var_real2d) :: &
      & nhx_deposition_const, &           !< constant value (avrg of 1st year of input value) of nhx_deposition [mumol m-2 sec-1]
      & nhy_deposition_const, &           !< constant value (avrg of 1st year of input value) of nhy_deposition [mumol m-2 sec-1]
      & p_deposition_const, &             !< constant value (avrg of 1st year of input value) of p_deposition [mumol m-2 sec-1]
      & co2_mixing_ratio_const, &         !< constant value (avrg of 1st year of input value) of co2_mixing_ratio [ppm]
      & co2_dC13_const, &                 !< constant value (avrg of 1st year of input value) of co2_dC13 [ppm]
      & co2_DC14_const                    !< constant value (avrg of 1st year of input value) of co2_DC14 [scaled]
#endif
#endif

  CONTAINS
    PROCEDURE :: Init => Init_a2l_memory
  END TYPE t_a2l_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_a2l_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !>initialize memory (variables) for the process: atm2land
  !>
  SUBROUTINE Init_a2l_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
    USE mo_jsb_class,         ONLY: Get_model
#ifndef __NO_QUINCY__
    ! quincy - needed for init of variables @TODO move to actual init routine
    USE mo_atmland_constants, ONLY: def_wind, def_co2_deltaC13, def_co2_deltaC14, def_t_air, standard_press_srf, def_vpd, &
                                    eps_vpd, def_co2_mixing_ratio, def_swdown, frac_vis_swtot_srf, frac_par_swvis_srf, &
                                    def_fdiffuse, def_cos_angle
    USE mo_isotope_util,      ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C
    USE mo_atmland_util,      ONLY: calc_spec_humidity_sat
#endif

    CLASS(t_a2l_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)            :: prefix          !< process name
    CHARACTER(len=*),    INTENT(in)            :: suffix          !< tile name
    INTEGER,             INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,             INTENT(in)            :: model_id        !< model ID model\%id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid    ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface  ! Vertical grid
    INTEGER                    :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_a2l_memory'

#ifndef __NO_QUINCY__
    ! quincy - local variables needed for init of memory variables @TODO move to actual init routine
    REAL(wp) :: C13C12_mixing_ratio, &
                C14C_mixing_ratio
    REAL(wp) :: wind_air_init, &
                wind_10m_init, &
                drag_srf_init, &
                ga_init, &
                q_air_init, &
                t_air_init, &
                press_srf_init, &
                swvis_srf_down_init, &
                swnir_srf_down_init, &
                swpar_srf_down_init, &
                fract_par_diffuse_init, &
                cos_zenith_angle_init, &
                co2_mixing_ratio_init, &
                co2_mixing_ratio_c13_init, &
                co2_mixing_ratio_c14_init
#endif

    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model        => Get_model(model_id)
    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')

    table        = tables(1)

#ifndef __NO_QUINCY__
    ! @TODO move to atm2land_init and init a2l variables there
    IF (model%config%model_scheme == MODEL_QUINCY) THEN
      C13C12_mixing_ratio         = calc_mixing_ratio_C13C12(def_co2_deltaC13)
      C14C_mixing_ratio           = calc_mixing_ratio_C14C(def_co2_deltaC13,def_co2_deltaC14)
      t_air_init                  = def_t_air
      press_srf_init              = standard_press_srf
      q_air_init                  = calc_spec_humidity_sat(t_air_init,press_srf_init) - def_vpd / press_srf_init * eps_vpd
      wind_air_init               = def_wind
      wind_10m_init               = wind_air_init
      drag_srf_init               = 1000.0_wp
      ga_init                     = 40._wp / drag_srf_init * wind_10m_init
      swvis_srf_down_init         = def_swdown * frac_vis_swtot_srf
      swnir_srf_down_init         = def_swdown * (1._wp - frac_vis_swtot_srf)
      swpar_srf_down_init         = swvis_srf_down_init * frac_par_swvis_srf
      fract_par_diffuse_init      = def_fdiffuse !TODO: is this still required? Probably for QS?!
      cos_zenith_angle_init       = def_cos_angle
      co2_mixing_ratio_init       = def_co2_mixing_ratio
      co2_mixing_ratio_c13_init   = co2_mixing_ratio_init / (1._wp + 1._wp / C13C12_mixing_ratio)
      co2_mixing_ratio_c14_init   = C14C_mixing_ratio * co2_mixing_ratio_init
    END IF
#endif

    ! ----------------------------------------------------------------------------------------------------- !
    ! default jsbach4 memory
    !
    CALL mem%Add_var( 'DEBUG_VAR', mem%DEBUG_VAR,                            &
      & hgrid, surface,                                                      &
      & t_cf('DEBUG_VAR', '', ''),                                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=777.7_wp )

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

    ! --------------------------------------------------------------------------------------------------- !
    ! variables used with jsbach4 & quincy
    ! jsbach4 settings
    !
    IF (model%config%model_scheme == MODEL_JSBACH) THEN

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
        & lrestart=jsbach_runs_standalone(),                                   &
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

    ! lakes
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
        & output_level=BASIC,                                                  &
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

    END IF ! lakes
    END IF ! IF (model%config%model_scheme == MODEL_JSBACH) THEN

#ifndef __NO_QUINCY__
    ! ----------------------------------------------------------------------------------------------------- !
    ! quincy
    IF (model%config%model_scheme == MODEL_QUINCY) THEN

      ! ------------------------------------------------------------------------------------------------- !
      ! variables used with jsbach4 & quincy
      ! quincy settings
      !
      CALL mem%Add_var('t_air', mem%t_air, &
        & hgrid, surface, &
        & t_cf('t_air', 'K', '2m air temperature'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & output_level=BASIC, &
        & initval_r = t_air_init)

      CALL mem%Add_var('q_air', mem%q_air, &
        & hgrid, surface, &
        & t_cf('q_air', 'g g-1', '2m air specific humidity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & output_level=BASIC, &
        & initval_r = q_air_init)

      CALL mem%Add_var('rain', mem%rain, &
        & hgrid, surface, &
        & t_cf('rain', 'kg m-2 s-1', 'liquid water precipitation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & output_level=BASIC, &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('snow', mem%snow, &
        & hgrid, surface, &
        & t_cf('snow', 'kg m-2 s-1', 'solid water precipitation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & output_level=BASIC, &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('wind_air', mem%wind_air, &
        & hgrid, surface, &
        & t_cf('wind_air', 'm s-1', '2m air wind speed'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = wind_air_init)

      CALL mem%Add_var('wind_10m', mem%wind_10m, &
        & hgrid, surface, &
        & t_cf('wind_10m', 'm s-1', '10m air wind speed'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = wind_10m_init)

      CALL mem%Add_var('lw_srf_down', mem%lw_srf_down, &
        & hgrid, surface, &
        & t_cf('lw_srf_down', 'W m-2', 'longwave surface downward flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('swvis_srf_down', mem%swvis_srf_down, &
        & hgrid, surface, &
        & t_cf('swvis_srf_down', 'W m-2', 'shortwave surface downward flux in the visible band'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & output_level=BASIC, &
        & initval_r = swvis_srf_down_init)

      CALL mem%Add_var('swnir_srf_down', mem%swnir_srf_down, &
        & hgrid, surface, &
        & t_cf('swnir_srf_down', 'W m-2', 'shortwave surface downward flux in the near infrared band'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & output_level=BASIC, &
        & initval_r = swnir_srf_down_init)

      CALL mem%Add_var('swpar_srf_down', mem%swpar_srf_down, &
        & hgrid, surface, &
        & t_cf('swpar_srf_down', 'W m-2', 'shortwave surface downward flux in the PAR band'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = swpar_srf_down_init)

      CALL mem%Add_var('fract_par_diffuse', mem%fract_par_diffuse, &
        & hgrid, surface, &
        & t_cf('fract_par_diffuse', 'fraction', 'fraction of shortwave surface downward flux that is diffuse'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = fract_par_diffuse_init)

      CALL mem%Add_var('press_srf', mem%press_srf, &
        & hgrid, surface, &
        & t_cf('press_srf', 'Pa', 'surface air pressure'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & output_level=BASIC, &
        & initval_r = press_srf_init)

      CALL mem%Add_var('drag_srf', mem%drag_srf, &
        & hgrid, surface, &
        & t_cf('drag_srf', 'unitless', 'total surface drag'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = drag_srf_init)

      CALL mem%Add_var('t_acoef', mem%t_acoef, &
        & hgrid, surface, &
        & t_cf('t_acoef', 'unitless', 'Richtmeyer-Morten coefficient A for temperature'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_bcoef', mem%t_bcoef, &
        & hgrid, surface, &
        & t_cf('t_bcoef', 'unitless', 'Richtmeyer-Morten coefficient B for temperature'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('q_acoef', mem%q_acoef, &
        & hgrid, surface, &
        & t_cf('q_acoef', 'unitless', 'Richtmeyer-Morten coefficient A for specific humidity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('q_bcoef', mem%q_bcoef, &
        & hgrid, surface, &
        & t_cf('q_bcoef', 'unitless', 'Richtmeyer-Morten coefficient B for specific humidity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('pch', mem%pch, &
        & hgrid, surface, &
        & t_cf('pch', 'unitless', 'total surface drag'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('cos_zenith_angle', mem%cos_zenith_angle, &
        & hgrid, surface, &
        & t_cf('cos_zenith_angle', 'unitless', 'Cosine of solar zenith angle'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = cos_zenith_angle_init)

      ! lakes
      IF (One_of(LAKE_TYPE, lct_ids(:)) > 0) THEN

        CALL mem%Add_var('drag_wtr', mem%drag_wtr, &
          & hgrid, surface, &
          & t_cf('drag_wtr', 'unitless', 'drag coefficient over water'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('drag_ice', mem%drag_ice, &
          & hgrid, surface, &
          & t_cf('drag_ice', 'unitless', 'drag coefficient over ice'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('t_acoef_wtr', mem%t_acoef_wtr, &
          & hgrid, surface, &
          & t_cf('t_acoef_wtr', 'unitless', 'Richtmeyer-Morten coefficient A for temperature over water'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('t_bcoef_wtr', mem%t_bcoef_wtr, &
          & hgrid, surface, &
          & t_cf('t_bcoef_wtr', 'unitless', 'Richtmeyer-Morten coefficient A for temperature over water'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('q_acoef_wtr', mem%q_acoef_wtr, &
          & hgrid, surface, &
          & t_cf('q_acoef_wtr', 'unitless', 'Richtmeyer-Morten coefficient A for specific humidity over water'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('q_bcoef_wtr', mem%q_bcoef_wtr, &
          & hgrid, surface, &
          & t_cf('q_bcoef_wtr', 'unitless', 'Richtmeyer-Morten coefficient A for specific humidity over water'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('t_acoef_ice', mem%t_acoef_ice, &
          & hgrid, surface, &
          & t_cf('t_acoef_ice', 'unitless', 'Richtmeyer-Morten coefficient A for temperature over ice'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('t_bcoef_ice', mem%t_bcoef_ice, &
          & hgrid, surface, &
          & t_cf('t_bcoef_ice', 'unitless', 'Richtmeyer-Morten coefficient A for temperature over ice'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('q_acoef_ice', mem%q_acoef_ice, &
          & hgrid, surface, &
          & t_cf('q_acoef_ice', 'unitless', 'Richtmeyer-Morten coefficient A for specific humidity over ice'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('q_bcoef_ice', mem%q_bcoef_ice, &
          & hgrid, surface, &
          & t_cf('q_bcoef_ice', 'unitless', 'Richtmeyer-Morten coefficient A for specific humidity over ice'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

      END IF ! lakes

      ! ------------------------------------------------------------------------------------------------- !
      ! variables used only with quincy
      !
      CALL mem%Add_var('co2_mixing_ratio', mem%co2_mixing_ratio, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio', 'ppm', 'surface air co2 mixing ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = co2_mixing_ratio_init)

      CALL mem%Add_var('co2_mixing_ratio_c13', mem%co2_mixing_ratio_c13, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio_C13', 'ppm', 'surface air 13co2 mixing ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = co2_mixing_ratio_c13_init)

      CALL mem%Add_var('co2_mixing_ratio_C14', mem%co2_mixing_ratio_c14, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio_c14', 'scaled', 'surface air 14co2 mixing ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = co2_mixing_ratio_c14_init)

      CALL mem%Add_var('ga', mem%ga, &
        & hgrid, surface, &
        & t_cf('ga', 'm s-1', 'aerodynamic conductance'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = ga_init)

      CALL mem%Add_var('nhx_deposition', mem%nhx_deposition, &
        & hgrid, surface, &
        & t_cf('nhx_deposition', 'mumol m-2 sec-1', 'surface downward reduced nitrogen deposition velocity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & output_level=BASIC, &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('noy_deposition', mem%noy_deposition, &
        & hgrid, surface, &
        & t_cf('noy_deposition', 'mumol m-2 sec-1', 'surface downward oxidised nitrogen deposition velocity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & output_level=BASIC, &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nhx_n15_deposition', mem%nhx_n15_deposition, &
        & hgrid, surface, &
        & t_cf('nhx_n15_deposition', 'mumol m-2 sec-1', 'surface downward reduced nitrogen-15 deposition velocity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('noy_n15_deposition', mem%noy_n15_deposition, &
        & hgrid, surface, &
        & t_cf('noy_n15_deposition', 'mumol m-2 sec-1', 'surface downward oxidised nitrogen-15 deposition velocity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('p_deposition', mem%p_deposition, &
        & hgrid, surface, &
        & t_cf('p_deposition', 'mumol m-2 sec-1', 'surface downward reduced phosphoruns deposition velocity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & output_level=BASIC, &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('local_time_day_seconds', mem%local_time_day_seconds, &
        & hgrid, surface, &
        & t_cf('local_time_day_seconds', 's', 'local time of this grid-cell'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('daytime_counter', mem%daytime_counter, &
        & hgrid, surface, &
        & t_cf('daytime_counter', '', 'number of timesteps of current day with light available'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    END IF ! IF (model%config%model_scheme == MODEL_QUINCY) THEN
#endif

  END SUBROUTINE Init_a2l_memory

#endif
END MODULE mo_a2l_memory_class
