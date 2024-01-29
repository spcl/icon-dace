!> Interface to run JSBACH for one time step
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
MODULE mo_atmland_interface
#ifndef __NO_JSBACH__

  USE mo_exception,       ONLY: finish
  USE mo_kind,            ONLY: wp

  USE mo_jsb_control,        ONLY: jsbach_runs_standalone
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_lct_class,      ONLY: LAKE_TYPE
  !USE mo_jsb_config_class,   ONLY: t_jsb_config
  !USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options

  dsl4jsb_Use_processes A2L_, SEB_, TURB_, HYDRO_, RAD_, CARBON_
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(CARBON_)


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_atm2land, update_land2atm

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_interface'

CONTAINS

  SUBROUTINE update_atm2land( &
    & tile, options,          &
    & t_air,                  &
    & q_air,                  &
    & press_air,              &
    & rain,                   &
    & snow,                   &
    & wind_air,               &
    & wind_10m,               &
    & lw_srf_down,            &
    & swvis_srf_down,         &
    & swnir_srf_down,         &
    & swpar_srf_down,         &
    & fract_par_diffuse,      &
    & dz_srf,                 &
    & press_srf,              &
    & rho_srf,                &
    & drag_srf,               &
    & t_acoef,                &
    & t_bcoef,                &
    & q_acoef,                &
    & q_bcoef,                &
    & pch,                    &
    & cos_zenith_angle,       &
    & CO2_air,                &
    & DEBUG_VAR,              &
    ! For lakes:
    & drag_wtr,               &
    & drag_ice,               &
    & t_acoef_wtr,            &
    & t_bcoef_wtr,            &
    & q_acoef_wtr,            &
    & q_bcoef_wtr,            &
    & t_acoef_ice,            &
    & t_bcoef_ice,            &
    & q_acoef_ice,            &
    & q_bcoef_ice             &
    & )

    USE mo_jsb_physical_constants, ONLY: molarMassDryAir, molarMassCO2

    CLASS(t_jsb_tile_abstract), INTENT(inout)             :: tile
    TYPE(t_jsb_task_options),   INTENT(in)                :: options
    ! TODO This line doesn't work with intel on mistral and standalone JSBACH with ECHAM infrastructure ... why?
    ! REAL(wp), OPTIONAL, DIMENSION(options%nc), INTENT(in) ::                                                       &
    REAL(wp), OPTIONAL, DIMENSION(:), INTENT(in) ::                                                                &
      & t_air, q_air, press_air, rain, snow, wind_air, wind_10m, lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, &
      & fract_par_diffuse, dz_srf, press_srf, rho_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle, &
      & CO2_air, DEBUG_VAR, drag_wtr, drag_ice,                                                                    &
      & t_acoef_wtr, t_bcoef_wtr, q_acoef_wtr, q_bcoef_wtr,                                                        &
      & t_acoef_ice, t_bcoef_ice, q_acoef_ice, q_bcoef_ice

    dsl4jsb_Def_memory(A2L_)

    INTEGER ::  iblk, ics, ice, nc, i
    TYPE(t_jsb_model), POINTER :: model
    LOGICAL :: use_quincy
    LOGICAL :: tile_contains_lake

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_atm2land'

    dsl4jsb_Real2D_onChunk ::  &
      & t_air_ptr,             &
      & q_air_ptr,             &
      & press_air_ptr,         &
      & rain_ptr,              &
      & snow_ptr,              &
      & wind_air_ptr,          &
      & wind_10m_ptr,          &
      & lw_srf_down_ptr,       &
      & swvis_srf_down_ptr,    &
      & swnir_srf_down_ptr,    &
      & swpar_srf_down_ptr,    &
      & fract_par_diffuse_ptr, &
      & dz_srf_ptr,            &
      & press_srf_ptr,         &
      & rho_srf_ptr,           &
      & drag_srf_ptr,          &
      & t_acoef_ptr,           &
      & t_bcoef_ptr,           &
      & q_acoef_ptr,           &
      & q_bcoef_ptr,           &
      & pch_ptr,               &
      & cos_zenith_angle_ptr,  &
      & CO2_air_ptr,           &
      & CO2_air_mol_ptr,       &
      & CO2_mixing_ratio_ptr,  &
      & DEBUG_VAR_ptr,         &
      & drag_wtr_ptr,          &
      & drag_ice_ptr,          &
      & t_acoef_wtr_ptr,       &
      & t_bcoef_wtr_ptr,       &
      & q_acoef_wtr_ptr,       &
      & q_bcoef_wtr_ptr,       &
      & t_acoef_ice_ptr,       &
      & t_bcoef_ice_ptr,       &
      & q_acoef_ice_ptr,       &
      & q_bcoef_ice_ptr

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    model => Get_model(tile%owner_model_id)
    use_quincy = model%config%use_quincy

    ! make this logical accessible in the OpenACC code directly
    tile_contains_lake = tile%contains_lake

    IF (nc /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(A2L_)
    IF (PRESENT(DEBUG_VAR)) THEN
      DEBUG_VAR_ptr         => dsl4jsb_var2D_onChunk(A2L_, DEBUG_VAR)
    END IF

    IF (PRESENT(t_air)) THEN
      t_air_ptr             => dsl4jsb_var2D_onChunk(A2L_, t_air)
    END IF
    IF (PRESENT(q_air)) THEN
      q_air_ptr             => dsl4jsb_var2D_onChunk(A2L_, q_air)
    END IF
    IF (PRESENT(press_air)) THEN
      press_air_ptr         => dsl4jsb_var2D_onChunk(A2L_, press_air)
    END IF
    IF(PRESENT(rain)) THEN
      rain_ptr              => dsl4jsb_var2D_onChunk(A2L_, rain)
    END IF
    IF (PRESENT(snow)) THEN
      snow_ptr              => dsl4jsb_var2D_onChunk(A2L_, snow)
    END IF
    IF (PRESENT(wind_air)) THEN
      wind_air_ptr          => dsl4jsb_var2D_onChunk(A2L_, wind_air)
    END IF
    IF (PRESENT(wind_10m)) THEN
      wind_10m_ptr          => dsl4jsb_var2D_onChunk(A2L_, wind_10m)
    END IF
    IF (PRESENT(lw_srf_down)) THEN
      lw_srf_down_ptr       => dsl4jsb_var2D_onChunk(A2L_, lw_srf_down)
    END IF
    IF (PRESENT(swvis_srf_down)) THEN
      swvis_srf_down_ptr    => dsl4jsb_var2D_onChunk(A2L_, swvis_srf_down)
    END IF
    IF(PRESENT(swnir_srf_down)) THEN
      swnir_srf_down_ptr    => dsl4jsb_var2D_onChunk(A2L_, swnir_srf_down)
    END IF
    IF (PRESENT(swpar_srf_down)) THEN
      swpar_srf_down_ptr    => dsl4jsb_var2D_onChunk(A2L_, swpar_srf_down)
    END IF
    IF (PRESENT(fract_par_diffuse)) THEN
      fract_par_diffuse_ptr => dsl4jsb_var2D_onChunk(A2L_, fract_par_diffuse)
    END IF
    IF (PRESENT(dz_srf)) THEN
      dz_srf_ptr            => dsl4jsb_var2D_onChunk(A2L_, dz_srf)
    END IF
    IF (PRESENT(press_srf)) THEN
      press_srf_ptr         => dsl4jsb_var2D_onChunk(A2L_, press_srf)
    END IF
    IF (PRESENT(rho_srf)) THEN
      rho_srf_ptr           => dsl4jsb_var2D_onChunk(A2L_, rho_srf)
    END IF
    IF(PRESENT(drag_srf)) THEN
      drag_srf_ptr          => dsl4jsb_var2D_onChunk(A2L_, drag_srf)
    END IF
    IF (PRESENT(t_acoef)) THEN
      t_acoef_ptr           => dsl4jsb_var2D_onChunk(A2L_, t_acoef)
    END IF
    IF (PRESENT(t_bcoef)) THEN
      t_bcoef_ptr           => dsl4jsb_var2D_onChunk(A2L_, t_bcoef)
    END IF
    IF (PRESENT(q_acoef)) THEN
      q_acoef_ptr           => dsl4jsb_var2D_onChunk(A2L_, q_acoef)
    END IF
    IF (PRESENT(q_bcoef)) THEN
      q_bcoef_ptr           => dsl4jsb_var2D_onChunk(A2L_, q_bcoef)
    END IF
    IF (PRESENT(pch)) THEN
      pch_ptr               => dsl4jsb_var2D_onChunk(A2L_, pch)
    END IF
    IF (PRESENT(cos_zenith_angle)) THEN
      cos_zenith_angle_ptr  => dsl4jsb_var2D_onChunk(A2L_, cos_zenith_angle)
    END IF
    IF (PRESENT(CO2_air)) THEN
      CO2_air_ptr           => dsl4jsb_var2D_onChunk(A2L_, CO2_air)
      CO2_air_mol_ptr       => dsl4jsb_var2D_onChunk(A2L_, CO2_air_mol)
      IF (use_quincy) CO2_mixing_ratio_ptr  => dsl4jsb_var2D_onChunk(A2L_, CO2_mixing_ratio)
    END IF
    IF (PRESENT(drag_wtr)) THEN
      drag_wtr_ptr          => dsl4jsb_var2D_onChunk(A2L_, drag_wtr)
    END IF
    IF (PRESENT(drag_ice)) THEN
      drag_ice_ptr          => dsl4jsb_var2D_onChunk(A2L_, drag_ice)
    END IF
    IF (PRESENT(t_acoef_wtr)) THEN
      t_acoef_wtr_ptr       => dsl4jsb_var2D_onChunk(A2L_, t_acoef_wtr)
    END IF
    IF (PRESENT(t_bcoef_wtr)) THEN
      t_bcoef_wtr_ptr       => dsl4jsb_var2D_onChunk(A2L_, t_bcoef_wtr)
    END IF
    IF (PRESENT(q_acoef_wtr)) THEN
      q_acoef_wtr_ptr       => dsl4jsb_var2D_onChunk(A2L_, q_acoef_wtr)
    END IF
    IF (PRESENT(q_bcoef_wtr)) THEN
      q_bcoef_wtr_ptr       => dsl4jsb_var2D_onChunk(A2L_, q_bcoef_wtr)
    END IF
    IF (PRESENT(t_acoef_ice)) THEN
      t_acoef_ice_ptr       => dsl4jsb_var2D_onChunk(A2L_, t_acoef_ice)
    END IF
    IF (PRESENT(t_bcoef_ice)) THEN
      t_bcoef_ice_ptr       => dsl4jsb_var2D_onChunk(A2L_, t_bcoef_ice)
    END IF
    IF (PRESENT(q_acoef_ice)) THEN
      q_acoef_ice_ptr       => dsl4jsb_var2D_onChunk(A2L_, q_acoef_ice)
    END IF
    IF (PRESENT(q_bcoef_ice)) THEN
      q_bcoef_ice_ptr       => dsl4jsb_var2D_onChunk(A2L_, q_bcoef_ice)
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT)

    IF (PRESENT(DEBUG_VAR)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        DEBUG_VAR_ptr(i) = DEBUG_VAR(i)
      END DO
    END IF

    IF (PRESENT(t_air)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        t_air_ptr(i) = t_air(i)
      END DO
    END IF
    IF (PRESENT(q_air)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        q_air_ptr(i) = q_air(i)
      END DO
    END IF
    IF (PRESENT(press_air)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        press_air_ptr(i) = press_air(i)
      END DO
    END IF
    IF (PRESENT(rain)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        rain_ptr(i) = rain(i)
      END DO
    END IF
    IF (PRESENT(snow)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        snow_ptr(i) = snow(i)
      END DO
    END IF
    IF (PRESENT(wind_air)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        wind_air_ptr(i) = wind_air(i)
      END DO
    END IF
    IF (PRESENT(wind_10m)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        wind_10m_ptr(i) = wind_10m(i)
      END DO
    END IF
    IF (PRESENT(lw_srf_down)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        lw_srf_down_ptr(i) = lw_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swvis_srf_down)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        swvis_srf_down_ptr(i) = swvis_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swnir_srf_down)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        swnir_srf_down_ptr(i) = swnir_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swpar_srf_down)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        swpar_srf_down_ptr(i) = swpar_srf_down(i)
      END DO
    END IF
    IF (PRESENT(fract_par_diffuse)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        fract_par_diffuse_ptr(i)  = fract_par_diffuse(i)
      END DO
    END IF
    IF (PRESENT(dz_srf)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        dz_srf_ptr(i) = dz_srf(i)
      END DO
    END IF
    IF (PRESENT(press_srf)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        press_srf_ptr(i) = press_srf(i)
      END DO
    END IF
    IF (PRESENT(drag_srf)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        drag_srf_ptr(i) = drag_srf(i)
      END DO
    END IF
    IF (PRESENT(rho_srf)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        rho_srf_ptr(i) = rho_srf(i)
      END DO
    END IF
    IF (PRESENT(t_acoef)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        t_acoef_ptr(i) = t_acoef(i)
      END DO
    END IF
    IF (PRESENT(t_bcoef)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        t_bcoef_ptr(i) = t_bcoef(i)
      END DO
    END IF
    IF (PRESENT(q_acoef)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        q_acoef_ptr(i) = q_acoef(i)
      END DO
    END IF
    IF (PRESENT(q_bcoef)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        q_bcoef_ptr(i) = q_bcoef(i)
      END DO
    END IF
    IF (PRESENT(pch)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        pch_ptr(i) = pch(i)
      END DO
    END IF
    IF (PRESENT(cos_zenith_angle)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        cos_zenith_angle_ptr(i) = cos_zenith_angle(i)
      END DO
    END IF
    IF (PRESENT(CO2_air)) THEN
      !$ACC LOOP GANG VECTOR
      DO i=1,nc
        CO2_air_ptr(i) = CO2_air(i)
        ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
        CO2_air_mol_ptr(i) = CO2_air(i) * molarMassDryAir / molarMassCO2
        ! convert CO2 from "molar ratio (volume)" to "co2 mixing ratio ppmv" (quincy | used in e.g. update_canopy_fluxes)
        IF (use_quincy) CO2_mixing_ratio_ptr(i) = CO2_air_mol_ptr(i) * 1000000._wp
      END DO
    END IF

    IF (tile_contains_lake) THEN

      IF (PRESENT(drag_wtr)) THEN
      !$ACC LOOP GANG VECTOR
        DO i=1,nc
          drag_wtr_ptr(i) = drag_wtr(i)
        END DO
      END IF
      IF (PRESENT(drag_ice)) THEN
      !$ACC LOOP GANG VECTOR
        DO i=1,nc
          drag_ice_ptr(i) = drag_ice(i)
        END DO
      END IF
      IF (PRESENT(t_acoef_wtr)) THEN
      !$ACC LOOP GANG VECTOR
        DO i=1,nc
          t_acoef_wtr_ptr(i) = t_acoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(t_bcoef_wtr)) THEN
      !$ACC LOOP GANG VECTOR
        DO i=1,nc
          t_bcoef_wtr_ptr(i) = t_bcoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(q_acoef_wtr)) THEN
      !$ACC LOOP GANG VECTOR
        DO i=1,nc
          q_acoef_wtr_ptr(i) = q_acoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(q_bcoef_wtr)) THEN
      !$ACC LOOP GANG VECTOR
        DO i=1,nc
          q_bcoef_wtr_ptr(i) = q_bcoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(t_acoef_ice)) THEN
      !$ACC LOOP GANG VECTOR
        DO i=1,nc
          t_acoef_ice_ptr(i) = t_acoef_ice(i)
        END DO
      END IF
      IF (PRESENT(t_bcoef_ice)) THEN
      !$ACC LOOP GANG VECTOR
        DO i=1,nc
          t_bcoef_ice_ptr(i) = t_bcoef_ice(i)
        END DO
      END IF
      IF (PRESENT(q_acoef_ice)) THEN
      !$ACC LOOP GANG VECTOR
        DO i=1,nc
          q_acoef_ice_ptr(i) = q_acoef_ice(i)
        END DO
      END IF
      IF (PRESENT(q_bcoef_ice)) THEN
      !$ACC LOOP GANG VECTOR
        DO i=1,nc
          q_bcoef_ice_ptr(i) = q_bcoef_ice(i)
        END DO
      END IF

    END IF

    !$ACC END PARALLEL

  END SUBROUTINE update_atm2land

  SUBROUTINE update_land2atm(tile, options,                                          &
    & t_srf, t_srf_rad, t_eff_srf, qsat_srf, s_srf,                                  &
    & fact_q_air, fact_qsat_srf, evapopot,                                           &
    & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,                  &
    & rough_h_srf, rough_m_srf, q_snocpymlt,                                         &
    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,                            &
    & kh, km, kh_neutral, km_neutral, CO2_flux,                                      &
    & t_lwtr, t_lice, qsat_lwtr, qsat_lice, s_lwtr, s_lice,                          &
    & evapo_wtr, latent_hflx_wtr, sensible_hflx_wtr,                                 &
    & evapo_ice, latent_hflx_ice, sensible_hflx_ice,                                 &
    & ice_fract_lake,                                                                &
    & alb_vis_dir_wtr, alb_vis_dif_wtr, alb_nir_dir_wtr, alb_nir_dif_wtr,            &
    & albedo_lwtr, albedo_lice)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(wp), OPTIONAL, INTENT(out) ::                                               &
      & t_srf(:), t_srf_rad(:), t_eff_srf(:), qsat_srf(:), s_srf(:),                 &
      & fact_q_air(:), fact_qsat_srf(:), evapopot(:),                                &
      & evapotrans(:), latent_hflx(:), sensible_hflx(:), grnd_hflx(:), grnd_hcap(:), &
      & rough_h_srf(:), rough_m_srf(:), q_snocpymlt(:),                              &
      & t_lwtr(:), t_lice(:), qsat_lwtr(:), qsat_lice(:), s_lwtr(:), s_lice(:),      &
      & evapo_wtr(:), latent_hflx_wtr(:), sensible_hflx_wtr(:),                      &
      & evapo_ice(:), latent_hflx_ice(:), sensible_hflx_ice(:),                      &
      & ice_fract_lake(:),                                                           &
      & alb_vis_dir_wtr(:), alb_vis_dif_wtr(:), alb_nir_dir_wtr(:),                  &
      & alb_nir_dif_wtr(:), albedo_lwtr(:), albedo_lice(:),                          &
      & alb_vis_dir(:), alb_nir_dir(:), alb_vis_dif(:), alb_nir_dif(:),              &
      & CO2_flux(:), kh(:), km(:), kh_neutral(:), km_neutral(:)

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(CARBON_)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), POINTER :: land_tile
    dsl4jsb_Def_memory_tile(TURB_, land_tile)

    INTEGER ::  iblk, ics, ice, nc, i, j, it, nt
    REAL(wp) :: fact_lake(options%nc), fact_land(options%nc)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_land2atm'

    dsl4jsb_Real2D_onChunk :: &
      & t_filt_ptr, &
      & t_eff_ptr, &
      & qsat_star_ptr, &
      & s_star_ptr, &
      & fact_q_air_ptr, &
      & fact_q_air_land_tile, &
      & fact_qsat_srf_ptr, &
      & fact_qsat_srf_land_tile, &
      & evapopot_ptr, &
      & evapotrans_lnd_ptr, &
      & latent_hflx_lnd_ptr, &
      & sensible_hflx_lnd_ptr, &
      & forc_hflx_ptr, &
      & heat_cap_old_ptr, &
      & rough_h_ptr, &
      & rough_m_ptr, &
      & q_snocpymlt_ptr, &
      & alb_vis_lnd_ptr, &
      & alb_nir_lnd_ptr, &
      & kh_ptr, &
      & km_ptr, &
      & kh_neutral_ptr, &
      & km_neutral_ptr, &
      & co2flux_npp_2_atm_ta_ptr, &
      & co2flux_soilresp_2_atm_ta_ptr, &
      & co2flux_herb_2_atm_ta_ptr, &
      & co2flux_fire_all_2_atm_ta_ptr, &
      & t_lwtr_ptr, &
      & t_lice_ptr, &
      & qsat_lwtr_ptr, &
      & qsat_lice_ptr, &
      & s_lwtr_ptr, &
      & s_lice_ptr, &
      & evapo_wtr_ptr, &
      & latent_hflx_wtr_ptr, &
      & sensible_hflx_wtr_ptr, &
      & evapo_ice_ptr, &
      & latent_hflx_ice_ptr, &
      & fract_lice_ptr, &
      & albedo_lwtr_ptr, &
      & albedo_lice_ptr, &
      & sensible_hflx_ice_ptr

    ! avoid compiler warnings about dummy arguments not being used
    IF (PRESENT(alb_nir_dif_wtr)) CONTINUE
    IF (PRESENT(alb_nir_dir_wtr)) CONTINUE
    IF (PRESENT(alb_vis_dif_wtr)) CONTINUE
    IF (PRESENT(alb_vis_dir_wtr)) CONTINUE

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    ! IF (nc /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(RAD_)
    IF (tile%Is_process_active(CARBON_)) THEN
      dsl4jsb_Get_memory(CARBON_)
    END IF

    !$ACC DATA CREATE(fact_land, fact_lake)
    IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN

      ! With lakes, re-scale _lnd, _lice and _lwtr fluxes given back to the atmosphere so that they are relative to the
      ! jsbach grid box, not/only including lakes (JSBACH grid boxes do not include ocean fractions.). Re-scaling is
      ! necessary, as in the atmosphere fluxes are regarded as land-only (or lake-only) fluxes, and lake, land and ocean
      ! fractions are considered separately.
      ! Note that, if lakes are considered, the lake tile must be a direct child of the root tile.

      DO j=1,SIZE(tile%lcts)
        IF (tile%lcts(j)%id == LAKE_TYPE) EXIT         ! We get index j of the lake tile
      END DO

      !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(tile%lcts)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        fact_lake(i) = tile%lcts(j)%fract(ics+i-1,iblk)      ! Lake fraction
        fact_land(i) = 1._wp - fact_lake(i)            ! Land fraction
        IF (fact_lake(i) > 0._wp) THEN
          fact_lake(i) = 1._wp / fact_lake(i)          ! Lake factor: 1/lake_fraction
        ELSE
          fact_lake(i) = 1._wp
        END IF
        IF (fact_land(i) > 0._wp) THEN
          fact_land(i) = 1._wp / fact_land(i)          ! Land factor: 1/land_fraction
        ELSE
          fact_land(i) = 1._wp
        END IF
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

    ELSE    ! without lakes
 
      !$ACC PARALLEL LOOP DEFAULT(PRESENT)
      DO i=1,nc
        fact_land(i) = 1._wp
        fact_lake(i) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP

    END IF

    ! Exchange fields on the land tile
    ! --------------------------------

    ! For some variables the atmosphere expects values for the land tile (excluding lakes).
    ! Thus we define a pointer to the land tile.
    nt = tile%Get_no_of_children()
    DO it=1,nt
      IF (it == 1) THEN
        land_tile => tile%Get_first_child_tile()
      ELSE
        land_tile => land_tile%Get_next_sibling_tile()
      END IF
      ! Exit if the land tile was found
      IF (TRIM(land_tile%name) == 'land') EXIT
    END DO
    dsl4jsb_Get_memory_tile(TURB_, land_tile)

    IF (model%config%use_tmx) THEN
      IF (PRESENT(t_srf)) THEN
        t_filt_ptr => dsl4jsb_var2D_onChunk(SEB_, t)
        !$ACC DATA PRESENT(t_srf, t_filt_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          t_srf(i) = t_filt_ptr(i)
        END DO
        !$ACC END DATA
      END IF
      IF (PRESENT(t_srf_rad)) THEN
        t_filt_ptr => dsl4jsb_var2D_onChunk(SEB_, t_rad4)
        !$ACC DATA PRESENT(t_srf_rad, t_filt_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          t_srf_rad(i) = t_filt_ptr(i)**0.25_wp
        END DO
        !$ACC END DATA
      END IF
    ELSE
      IF (PRESENT(t_srf)) THEN
        t_filt_ptr => dsl4jsb_var2D_onChunk(SEB_, t_filt)
        !$ACC DATA PRESENT(t_srf, t_filt_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          t_srf(i) = t_filt_ptr(i)
        END DO
        !$ACC END DATA
      END IF
    END IF
    IF (PRESENT(t_eff_srf)) THEN
      t_eff_ptr => dsl4jsb_var2D_onChunk(SEB_,   t_eff4)
      !$ACC DATA PRESENT(t_eff_srf, t_eff_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        t_eff_srf(i) = t_eff_ptr(i)**0.25_wp
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(qsat_srf)) THEN
      qsat_star_ptr => dsl4jsb_var2D_onChunk(SEB_,   qsat_star)
      !$ACC DATA PRESENT(qsat_srf, qsat_star_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        qsat_srf(i) = qsat_star_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(s_srf)) THEN
      s_star_ptr => dsl4jsb_var2D_onChunk(SEB_,   s_star)
      !$ACC DATA PRESENT(s_srf, s_star_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        s_srf(i) = s_star_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(fact_q_air)) THEN
      dsl4jsb_Get_var2d_onChunk_tile_name(TURB_, fact_q_air, land_tile)
      !$ACC DATA PRESENT(fact_q_air, fact_q_air_land_tile)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        fact_q_air(i) = fact_q_air_land_tile(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(fact_qsat_srf)) THEN
      dsl4jsb_Get_var2d_onChunk_tile_name(TURB_, fact_qsat_srf, land_tile)
      !$ACC DATA PRESENT(fact_qsat_srf, fact_qsat_srf_land_tile)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        fact_qsat_srf(i) = fact_qsat_srf_land_tile(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(evapopot)) THEN
      evapopot_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapopot) ! For offline only
      !$ACC DATA PRESENT(evapopot, evapopot_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        evapopot(i) = evapopot_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(evapotrans)) THEN
      IF (model%config%use_tmx) THEN
        evapotrans_lnd_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapotrans)
      ELSE
        evapotrans_lnd_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapotrans_lnd)
      END IF
      !$ACC DATA PRESENT(evapotrans, fact_land, evapotrans_lnd_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        evapotrans(i) = fact_land(i) * evapotrans_lnd_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(latent_hflx)) THEN
      IF (model%config%use_tmx) THEN
        latent_hflx_lnd_ptr => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx)
      ELSE
        latent_hflx_lnd_ptr => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_lnd)
      END IF
      !$ACC DATA PRESENT(latent_hflx, fact_land, latent_hflx_lnd_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        latent_hflx(i) = fact_land(i) * latent_hflx_lnd_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(sensible_hflx)) THEN
      IF (model%config%use_tmx) THEN
        sensible_hflx_lnd_ptr => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx)
      ELSE
        sensible_hflx_lnd_ptr => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_lnd)
      END IF
      !$ACC DATA PRESENT(sensible_hflx, fact_land, sensible_hflx_lnd_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        sensible_hflx(i) = fact_land(i) * sensible_hflx_lnd_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(grnd_hflx)) THEN
      forc_hflx_ptr => dsl4jsb_var2D_onChunk(SEB_,   forc_hflx)    ! TODO Not used
      !$ACC DATA PRESENT(grnd_hflx, forc_hflx_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        grnd_hflx(i) = forc_hflx_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(grnd_hcap)) THEN
      heat_cap_old_ptr => dsl4jsb_var2D_onChunk(SEB_,   heat_cap_old)  ! not computed for lake
      !$ACC DATA PRESENT(grnd_hcap, heat_cap_old_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        grnd_hcap(i) = heat_cap_old_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(rough_h_srf)) THEN
      rough_h_ptr => dsl4jsb_var2D_onChunk(TURB_,  rough_h)
      !$ACC DATA PRESENT(rough_h_srf, rough_h_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        rough_h_srf(i) = rough_h_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(rough_m_srf)) THEN
      rough_m_ptr => dsl4jsb_var2D_onChunk(TURB_,  rough_m)
      !$ACC DATA PRESENT(rough_m_srf, rough_m_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        rough_m_srf(i) = rough_m_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(q_snocpymlt)) THEN
      q_snocpymlt_ptr => dsl4jsb_var2D_onChunk(HYDRO_, q_snocpymlt)
      !$ACC DATA PRESENT(q_snocpymlt, q_snocpymlt_ptr, fact_land)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        q_snocpymlt(i) = fact_land(i) * q_snocpymlt_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(alb_vis_dir)) THEN
      IF (model%config%use_tmx) THEN
        alb_vis_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_, alb_vis)      ! TODO
      ELSE
        alb_vis_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_, alb_vis_lnd)  ! TODO
      END IF
      !$ACC DATA PRESENT(alb_vis_lnd_ptr, alb_vis_dir)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        alb_vis_dir(i) = alb_vis_lnd_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(alb_vis_dif)) THEN
      IF (model%config%use_tmx) THEN
        alb_vis_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_, alb_vis)      ! TODO
      ELSE
        alb_vis_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_, alb_vis_lnd)  ! TODO
      END IF
      !$ACC DATA PRESENT(alb_vis_dif, alb_vis_lnd_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        alb_vis_dif(i) = alb_vis_lnd_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(alb_nir_dir)) THEN
      IF (model%config%use_tmx) THEN
        alb_nir_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_, alb_nir)      ! TODO
      ELSE
        alb_nir_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_, alb_nir_lnd)  ! TODO
      END IF
      !$ACC DATA PRESENT(alb_nir_dir, alb_nir_lnd_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        alb_nir_dir(i) = alb_nir_lnd_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(alb_nir_dif)) THEN
      IF (model%config%use_tmx) THEN
        alb_nir_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_, alb_nir)
      ELSE
        alb_nir_lnd_ptr => dsl4jsb_var2D_onChunk(RAD_, alb_nir_lnd)
      END IF
      !$ACC DATA PRESENT(alb_nir_dif, alb_nir_lnd_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        alb_nir_dif(i) = alb_nir_lnd_ptr(i)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC END DATA
    END IF
    IF (PRESENT(CO2_flux)) THEN
      IF (tile%Is_process_active(CARBON_)) THEN
        co2flux_npp_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,co2flux_npp_2_atm_ta)
        co2flux_soilresp_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,co2flux_soilresp_2_atm_ta)
        co2flux_herb_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,co2flux_herb_2_atm_ta)
        co2flux_fire_all_2_atm_ta_ptr => dsl4jsb_var2D_onChunk(CARBON_,co2flux_fire_all_2_atm_ta)

        !$ACC DATA PRESENT(CO2_flux, fact_land, co2flux_npp_2_atm_ta_ptr, co2flux_soilresp_2_atm_ta_ptr) &
        !$ACC   PRESENT(co2flux_herb_2_atm_ta_ptr, co2flux_fire_all_2_atm_ta_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          CO2_flux(i) = fact_land(i) *              &
            & (   co2flux_npp_2_atm_ta_ptr(i)       &
            &   + co2flux_soilresp_2_atm_ta_ptr(i)  &
            &   + co2flux_herb_2_atm_ta_ptr(i)      &
            &   + co2flux_fire_all_2_atm_ta_ptr(i)  &
            & )
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        CO2_flux(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(kh)) THEN
      kh_ptr => dsl4jsb_var2D_onChunk(TURB_, kh)
      !$ACC DATA PRESENT(kh, kh_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        kh(i) = kh_ptr(i)
      END DO
      !$ACC END DATA
    END IF
    IF (PRESENT(km)) THEN
      km_ptr => dsl4jsb_var2D_onChunk(TURB_, km)
      !$ACC DATA PRESENT(km, km_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        km(i) = km_ptr(i)
      END DO
      !$ACC END DATA
    END IF
    IF (PRESENT(kh_neutral)) THEN
      kh_neutral_ptr => dsl4jsb_var2D_onChunk(TURB_, kh_neutral)
      !$ACC DATA PRESENT(kh_neutral, kh_neutral_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        kh_neutral(i) = kh_neutral_ptr(i)
      END DO
      !$ACC END DATA
    END IF
    IF (PRESENT(km_neutral)) THEN
      km_neutral_ptr => dsl4jsb_var2D_onChunk(TURB_, km_neutral)
      !$ACC DATA PRESENT(km_neutral, km_neutral_ptr)
      !$ACC PARALLEL LOOP
      DO i=1,nc
        km_neutral(i) = km_neutral_ptr(i)
      END DO
      !$ACC END DATA
    END IF

    ! Exchange fields on lake water fractions
    ! ---------------------------------------

    IF (PRESENT(t_lwtr)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN
        t_lwtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   t_lwtr)
        !$ACC DATA PRESENT(t_lwtr_ptr, t_lwtr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          t_lwtr(i) = t_lwtr_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
      !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        t_lwtr(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(qsat_lwtr)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN
        qsat_lwtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   qsat_lwtr)
        !$ACC DATA PRESENT(qsat_lwtr, qsat_lwtr_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          qsat_lwtr(i) = qsat_lwtr_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        qsat_lwtr(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(s_lwtr)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN
        s_lwtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   s_lwtr)
        !$ACC DATA PRESENT(s_lwtr, s_lwtr_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          s_lwtr(i) = s_lwtr_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
          s_lwtr(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(evapo_wtr)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN
        evapo_wtr_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapo_wtr)
        !$ACC DATA PRESENT(evapo_wtr, evapo_wtr_ptr, fact_lake)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          evapo_wtr(i) = fact_lake(i) * evapo_wtr_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        evapo_wtr(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(latent_hflx_wtr)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN
        latent_hflx_wtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_wtr)
        !$ACC DATA PRESENT(latent_hflx_wtr, latent_hflx_wtr_ptr, fact_lake)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          latent_hflx_wtr(i) = fact_lake(i) * latent_hflx_wtr_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        latent_hflx_wtr(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(sensible_hflx_wtr)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN
        sensible_hflx_wtr_ptr => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_wtr)
        !$ACC DATA PRESENT(sensible_hflx_wtr, sensible_hflx_wtr_ptr, fact_lake)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          sensible_hflx_wtr(i) = fact_lake(i) * sensible_hflx_wtr_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        sensible_hflx_wtr(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(albedo_lwtr)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN
        albedo_lwtr_ptr => dsl4jsb_var2D_onChunk(RAD_,   albedo_lwtr)
        !$ACC DATA PRESENT(albedo_lwtr, albedo_lwtr_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          albedo_lwtr(i) = albedo_lwtr_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        albedo_lwtr(:) = 0.07_wp
        !$ACC END KERNELS
      END IF
    END IF

    ! Exchange fields on the ice fraction of lakes
    ! --------------------------------------------

    IF (PRESENT(ice_fract_lake)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone()) THEN
        fract_lice_ptr => dsl4jsb_var2D_onChunk(SEB_,   fract_lice)
        !$ACC DATA PRESENT(ice_fract_lake, fract_lice_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          ice_fract_lake(i) = fract_lice_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        ice_fract_lake(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF

    IF (PRESENT(t_lice)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone() .AND. dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        t_lice_ptr => dsl4jsb_var2D_onChunk(SEB_,   t_lice)
        !$ACC DATA PRESENT(t_lice, t_lice_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          t_lice(i) = t_lice_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        t_lice(:) = 273.15_wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(qsat_lice)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone() .AND. dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        qsat_lice_ptr => dsl4jsb_var2D_onChunk(SEB_,   qsat_lice)
        !$ACC DATA PRESENT(qsat_lice, qsat_lice_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          qsat_lice(i) = qsat_lice_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        qsat_lice(:) = 0.0075_wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(s_lice)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone() .AND. dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        s_lice_ptr => dsl4jsb_var2D_onChunk(SEB_,   s_lice)
        !$ACC DATA PRESENT(s_lice, s_lice_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          s_lice(i) = s_lice_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        s_lice(:) = 2.9E5_wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(evapo_ice)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone() .AND. dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        evapo_ice_ptr => dsl4jsb_var2D_onChunk(HYDRO_, evapo_ice)
        !$ACC DATA PRESENT(evapo_ice_ptr, evapo_ice, fact_lake)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          evapo_ice(i) = fact_lake(i) * evapo_ice_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        evapo_ice(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(latent_hflx_ice)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone() .AND. dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        latent_hflx_ice_ptr => dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_ice)
        !$ACC DATA PRESENT(latent_hflx_ice, latent_hflx_ice_ptr, fact_lake)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          latent_hflx_ice(i) = fact_lake(i) * latent_hflx_ice_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        latent_hflx_ice(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(sensible_hflx_ice)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone() .AND. dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        sensible_hflx_ice_ptr  => dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_ice)
        !$ACC DATA PRESENT(sensible_hflx_ice_ptr, sensible_hflx_ice, fact_lake)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          sensible_hflx_ice(i) = fact_lake(i) * sensible_hflx_ice_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        sensible_hflx_ice(:) = 0._wp
        !$ACC END KERNELS
      END IF
    END IF
    IF (PRESENT(albedo_lice)) THEN
      IF (tile%contains_lake .AND. .NOT. jsbach_runs_standalone() .AND. dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        albedo_lice_ptr => dsl4jsb_var2D_onChunk(RAD_,   albedo_lice)
        !$ACC DATA PRESENT(albedo_lice, albedo_lice_ptr)
        !$ACC PARALLEL LOOP
        DO i=1,nc
          albedo_lice(i) = albedo_lice_ptr(i)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT)
        albedo_lice(:) = 0.55_wp
        !$ACC END KERNELS
      END IF
    END IF

    !$ACC END DATA

  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
