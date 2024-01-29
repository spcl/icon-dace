!> Contains the memory class for the surface energy balance process.
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

MODULE mo_seb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, LAKE_TYPE, GLACIER_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real1d, t_jsb_var_real2d
  USE mo_jsb_physical_constants, ONLY: tmelt

  dsl4jsb_Use_processes SEB_
  dsl4jsb_Use_config(SEB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_seb_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 50

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_seb_memory
    !
    ! Common variables
    TYPE(t_jsb_var_real2d) ::   &
      & t,                      & !< Surface temperature                                             [K]
      & t_old,                  & !< Surface temperature at previous timestep                        [K]
      & t_unfilt,               & !< Surface temperature, unfiltered (before Asselin)                [K]
      & t_rad4,                 & !< Surface radiative temperature **4                               [K**4]
      & t_eff4,                 & !< Effective surface temperature **4                               [K**4]
                                  !! This is the surface temperature used in the atmospheric
                                  !! radiative heating
      & qsat_star,              & !< Saturation specific humidity at surface (qsat^star, see manual) [kg kg-1]
                                  !! This is the qsat_srf that is used in the vertical diffusion
                                  !! scheme in the atmosphere
      & heat_cap,               & !< Heat capacity of surface layer                                  [J m-2 K-1]
      & heat_cap_old,           & !< Heat capacity of surface layer at previous timestep             [J m-2 K-1]
      & latent_hflx,            & !< Latent heat flux                                                [W m-2]
      & sensible_hflx,          & !< Sensible heat flux                                              [W m-2]
      & forc_hflx                 !< Additional heat flux at surface                                 [W m-2]
                                  !! (additional to radiative fluxes and latent and sensible
                                  !!  heat fluxes, e.g. ground heat flux)
    !
    ! Additional variables for LAND lct_type
    TYPE(t_jsb_var_real2d) ::   &
      & s_star,                 & !< Surface dry static energy (s^star, see manual)                  [m2 s-2]
                                  !! This is the s that is used in the vertical diffusion
                                  !! scheme in the atmosphere
      & t_unfilt_old,           & !< Surface temperature at previous timestep, unfiltered            [K]
                                  !! (before Asselin)
      & t_filt,                 & !< Surface temperature, filtered (after Asselin)                   [K]
      & latent_hflx_lnd,        & !< Latent heat flux                                                [W m-2]
      & sensible_hflx_lnd,      & !< Sensible heat flux                                              [W m-2]
      & previous_day_temp_mean, & ! mean day air-temperature [Celsius] of previous day,
                                  ! because not known for ..
                                  ! current day; updated at the beginning of each new day (stream)
      & previous_day_temp_min,  & ! minimum air-temperature [Celsius] of previous day, because not known for  ..
                                  ! current day; updated at the beginning of each new day (stream)
      & previous_day_temp_max,  & ! maximum air-temperature [Celsius] of previous day, because not known for   ..
                                  ! current day; updated at the beginning of each new day (stream)
      & pseudo_soil_temp,       & ! Pseudo-soil-temperature [Celsius]: a weighted running mean of the land ...
                                  ! surface temperature, to simulate soil temperature. Updated every time step ..
                                  ! by the subroutine "pseudo_soil_temp()". (stream)
      & day_temp_sum,           & ! sum of air temperatures since midnight (needed to compute mean air
                                  ! temperature) [Celsius * 1] (stream)
      & day_temp_min,           & ! minimum of air temperatures since midnight (needed to compute growing
                                  ! degree days GDD) [Celsius * 1] (stream)
      & day_temp_max,           & ! maximum of air temperatures since midnight (needed to compute growing
                                  ! degree days GDD) [Celsius * 1] (stream)
      & N_pseudo_soil_temp,     & ! Normalization for computing the pseudo soil temperature
      & F_pseudo_soil_temp        ! exponential factor used for updating the pseudo soil temperature
    !
    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & t_lwtr,               & !< Lake surface temperature (water)                                [K]
      & t_lice,               & !< Lake surface temperature (ice)                                  [K]
!!$      & t_lwtr_old,            & !< Lake surface temperature (water) at previous timestep           [K]
!!$      & t_lice_old,            & !< Lake surface temperature (ice) at previous timestep             [K]
      & s_lwtr,               & !< Surface dry static energy over lake water                       [m2 s-2]
      & qsat_lwtr,            & !< Saturation specific humidity at surface over lake water        [kg kg-1]
      & s_lice,               & !< Surface dry static energy over lake ice                         [m2 s-2]
      & qsat_lice,            & !< Saturation specific humidity at surface over lake ice           [kg kg-1]
      & latent_hflx_wtr,      &
      & sensible_hflx_wtr,    &
      & latent_hflx_ice,      &
      & sensible_hflx_ice,    &
      & depth_lice,           & !< Thickness of ice on lake                                        [m]
      & fract_lice,           & !< Fraction of lake ice                                            []
      & hflx_form_ice,        & !< Heat flux for ice formation                                     [W m-2]
      & hflx_cond_ice,        & !< Conductive heat flux for ice on lakes
      & hflx_melt_lice,       & !< Heat flux from melting ice on lakes
      & hflx_melt_snow          !< Heat flux from melting snow on lake ice

    ! Diagnostic global land mean variables e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) ::        &
      t_gmean                   !< Global land mean surface temperature (box tile)                 [K]


  CONTAINS
    PROCEDURE :: Init => Init_seb_memory
  END TYPE t_seb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_memory_class'

CONTAINS

  SUBROUTINE Init_seb_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC !, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables !, TSTEP_CONSTANT
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_seb_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    dsl4jsb_Def_config(SEB_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_seb_memory'
    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    dsl4jsb_Get_config(SEB_)

    CALL mem%Add_var( 't', mem%t,                                                    &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp', 'K', 'surface temperature'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC, &
      & initval_r=280.0_wp )

    IF (model%config%use_tmx) THEN
      CALL mem%Add_var( 't_rad**4', mem%t_rad4,                                      &
        & hgrid, surface,                                                            &
        & t_cf('sfc_temp_rad**4', 'K', 'surface radiative temperature **4'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & output_level=BASIC,                                                        &
        & initval_r=280.0_wp**4 )
    END IF

    CALL mem%Add_var( 't_old', mem%t_old,                                            &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp_old', 'K', 'surface temperature at previous timestep'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=280.0_wp )

    CALL mem%Add_var( 't_unfilt', mem%t_unfilt,                                      &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp_unfilt', 'K', 'unfiltered surface temperature'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=280.0_wp )

    CALL mem%Add_var( 'heat_cap', mem%heat_cap,                                      &
      & hgrid, surface,                                                              &
      & t_cf('sfc_heat_cap', 'J m-2 K-1)', 'surface layer heat capacity'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'heat_cap_old', mem%heat_cap_old,                              &
      & hgrid, surface,                                                              &
      & t_cf('sfc_heat_cap_old', 'J m-2 K-1)', 'srf layer heat cap. prev timestep'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'latent_hflx', mem%latent_hflx,                                &
      & hgrid, surface,                                                              &
      & t_cf('latent_hflx', 'W m-2', 'latent heat flux density at surface'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'sensible_hflx', mem%sensible_hflx,                            &
      & hgrid, surface,                                                              &
      & t_cf('sensible_hflx', 'W m-2', 'sensible heat flux at surface'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'forc_hflx', mem%forc_hflx,                                    &
      & hgrid, surface,                                                              &
      & t_cf('forc_hflx', 'W m-2)', 'additional heat flux at surface'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    IF (model%config%use_tmx .OR. One_of(LAND_TYPE, lct_ids(:)) > 0) THEN
      CALL mem%Add_var( 't_eff**4', mem%t_eff4,                                              &
        & hgrid, surface,                                                                    &
        & t_cf('surface_temperature_eff**4', 'K', 'effective surface temperature**4'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=280.0_wp**4 )
    END IF

    ! Additional variables for LAND lct
    IF ( model%config%use_tmx .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 ) THEN

      CALL mem%Add_var( 's_star', mem%s_star,                                                &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_dry_static_energy', 'm2 s-2', 'surface dry static energy'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & output_level=BASIC, &
        & initval_r=2.9E5_wp )

      CALL mem%Add_var( 'qsat_star', mem%qsat_star,                                          &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_qsat', 'm2 s-2', 'surface specific humidity at saturation'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.TRUE.,                                                                   &
        & output_level=BASIC,                                                                &
        & initval_r=0.0075_wp )

      CALL mem%Add_var( 't_unfilt_old', mem%t_unfilt_old,                                    &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_temp_unfilt_old', 'K', 'unf. surface temperature at previous timestep'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=280.0_wp )

      CALL mem%Add_var( 't_filt', mem%t_filt,                                                &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_temp_filt', 'K', 'filtered surface temperature'),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=280.0_wp )

      CALL mem%Add_var( 'latent_hflx_lnd', mem%latent_hflx_lnd,                              &
        & hgrid, surface,                                                                    &
        & t_cf('latent_hflx_lnd', 'W m-2', 'latent heat flux density at surface (lnd)'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'sensible_hflx_lnd', mem%sensible_hflx_lnd,                          &
        & hgrid, surface,                                                                    &
        & t_cf('sensible_hflx_lnd', 'W m-2', 'sensible heat flux at surface (lnd)'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'previous_day_temp_mean', mem%previous_day_temp_mean,      &
        & hgrid, surface,                                                          &
        & t_cf('previous_day_temp_mean', '', ''),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
        & prefix, suffix,                                                          &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10.0_wp )

      CALL mem%Add_var( 'previous_day_temp_min', mem%previous_day_temp_min,      &
        & hgrid, surface,                                                        &
        & t_cf('previous_day_temp_min', '', ''),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=5.0_wp )

      CALL mem%Add_var( 'previous_day_temp_max', mem%previous_day_temp_max,      &
        & hgrid, surface,                                                        &
        & t_cf('albedo_nir_previous_day_temp_maxsnow', '', ''),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=15.0_wp )

      CALL mem%Add_var( 'pseudo_soil_temp', mem%pseudo_soil_temp,                &
        & hgrid, surface,                                                        &
        & t_cf('pseudo_soil_temp', '', ''),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10.0_wp )

      CALL mem%Add_var('day_temp_sum', mem%day_temp_sum,                         &
        & hgrid, surface,                                                        &
        & t_cf('day_temp_sum', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0._wp )

      CALL mem%Add_var('day_temp_min', mem%day_temp_min,                         &
        & hgrid, surface,                                                        &
        & t_cf('day_temp_min', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0._wp )

      CALL mem%Add_var('day_temp_max', mem%day_temp_max,                         &
        & hgrid, surface,                                                        &
        & t_cf('day_temp_max', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0._wp)

      CALL mem%Add_var('N_pseudo_soil_temp', mem%N_pseudo_soil_temp,             &
        & hgrid, surface,                                                        &
        & t_cf('N_pseudo_soil_temp', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10._wp )

      CALL mem%Add_var('F_pseudo_soil_temp', mem%F_pseudo_soil_temp,             &
        & hgrid, surface,                                                        &
        & t_cf('F_pseudo_soil_temp', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10._wp )

    END IF

    ! Additional variables for lakes
    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 ) THEN

      CALL mem%Add_var( 't_lwtr', mem%t_lwtr,                                          &
        & hgrid, surface,                                                              &
        & t_cf('sfc_temp_wtr', 'K', 'lake surface temperature wtr '),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=280.0_wp )

      CALL mem%Add_var( 's_lwtr', mem%s_lwtr,                                                  &
        & hgrid, surface,                                                                      &
        & t_cf('sfc_dry_static_energy_wtr', 'm2 s-2', 'surface dry static energy (lake wtr)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix,                                                                      &
        & lrestart=.FALSE.,                                                                    &
        & initval_r=2.9E5_wp )

      CALL mem%Add_var( 'qsat_lwtr', mem%qsat_lwtr,                                              &
        & hgrid, surface,                                                                        &
        & t_cf('sfc_qsat_lwtr', 'm2 s-2', 'surface specific humidity at saturation (lake wtr)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE.,                                                                      &
        & initval_r=0.0075_wp )

      CALL mem%Add_var( 'latent_hflx_wtr', mem%latent_hflx_wtr,                        &
        & hgrid, surface,                                                              &
        & t_cf('latent_hflx_wtr', 'W m-2', 'latent heat flux over lake water'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.TRUE.,                                                             &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'sensible_hflx_wtr', mem%sensible_hflx_wtr,                    &
        & hgrid, surface,                                                              &
        & t_cf('sensible_hflx_wtr', 'W m-2', 'sensible heat flux over lake water'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.TRUE.,                                                             &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! Fraction of ice on lake
      CALL mem%Add_var( 'fract_lice', mem%fract_lice,                                  &
        & hgrid, surface,                                                              &
        & t_cf('fract_lice', '', 'fraction of ice on lake'),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )
    
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      
        CALL mem%Add_var( 't_lice', mem%t_lice,                                          &
          & hgrid, surface,                                                              &
          & t_cf('sfc_temp_ice', 'K', 'lake surface temperature ice '),                  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=tmelt )
    
        CALL mem%Add_var( 's_lice', mem%s_lice,                                                  &
          & hgrid, surface,                                                                      &
          & t_cf('sfc_dry_static_energy_ice', 'm2 s-2', 'surface dry static energy (lake ice)'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
          & prefix, suffix,                                                                      &
          & lrestart=.FALSE.,                                                                    &
          & initval_r=2.9E5_wp )
    
        CALL mem%Add_var( 'qsat_lice', mem%qsat_lice,                                              &
          & hgrid, surface,                                                                        &
          & t_cf('sfc_qsat_lice', 'm2 s-2', 'surface specific humidity at saturation (lake ice)'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix,                                                                        &
          & lrestart=.FALSE.,                                                                      &
          & initval_r=0.0075_wp )
    
        CALL mem%Add_var( 'latent_hflx_ice', mem%latent_hflx_ice,                        &
          & hgrid, surface,                                                              &
          & t_cf('latent_hflx_ice', 'W m-2', 'latent heat flux over lake ice'),          &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & lrestart=.TRUE.,                                                             &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    
        CALL mem%Add_var( 'sensible_hflx_ice', mem%sensible_hflx_ice,                    &
          & hgrid, surface,                                                              &
          & t_cf('sensible_hflx_ice', 'W m-2', 'sensible heat flux over lake ice'),      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & lrestart=.TRUE.,                                                             &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    
        CALL mem%Add_var( 'depth_lice', mem%depth_lice,                                  &
          & hgrid, surface,                                                              &
          & t_cf('ice_thickness', 'm', 'thickness of ice on lake'),                      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp )
    
        CALL mem%Add_var( 'hflx_form_ice', mem%hflx_form_ice,                            &
          & hgrid, surface,                                                              &
          & t_cf('hflx_form_ice', 'W m-2)', 'heat flux for ice formation on lake ice'),  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    
        CALL mem%Add_var( 'hflx_cond_ice', mem%hflx_cond_ice,                            &
          & hgrid, surface,                                                              &
          & t_cf('hflx_cond_ice', 'W m-2)', 'conductive heat flux on lake ice'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    
        CALL mem%Add_var( 'hflx_melt_lice', mem%hflx_melt_lice,                          &
          & hgrid, surface,                                                              &
          & t_cf('hflx_melt_lice', 'W m-2)', 'heat flux from ice melt on lakes'),        &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    
        CALL mem%Add_var( 'hflx_melt_snow', mem%hflx_melt_snow,                          &
          & hgrid, surface,                                                              &
          & t_cf('hflx_melt_snow', 'W m-2)', 'heat_flux from snow melt on lakes'),       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        END IF

    END IF

#ifdef __ICON__
    ! Diagnostic global carbon sums for experiment monitoring
    !       (1d stream variables are not supported with echam)
    !
    IF (TRIM(suffix) == 'box') THEN

      CALL mem%Add_var( 't_gmean', mem%t_gmean,                                 &
        & hgrid, surface,                                                       &
        & t_cf('t_gmean', 'K', 'Global land mean surface temperature'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
        & prefix, suffix,                                                       &
        & lrestart=.FALSE., initval_r=280.0_wp )

    END IF
#endif

  END SUBROUTINE Init_seb_memory

#endif
END MODULE mo_seb_memory_class
