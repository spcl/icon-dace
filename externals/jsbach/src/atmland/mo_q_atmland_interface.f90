!> Interface to run ICON-Land with QUINCY for one time step
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
MODULE mo_q_atmland_interface
#ifndef __NO_QUINCY__

  USE mo_exception,       ONLY: finish
  USE mo_kind,            ONLY: wp

  USE mo_jsb_control,        ONLY: jsbach_runs_standalone
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid,           ONLY: Get_grid
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_lct_class,      ONLY: LAKE_TYPE
  !USE mo_jsb_config_class,   ONLY: t_jsb_config
  !USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options

  dsl4jsb_Use_processes A2L_, Q_RAD_, SPQ_, VEG_
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(Q_RAD_)
  dsl4jsb_Use_memory(SPQ_)
  dsl4jsb_Use_memory(VEG_)


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_atm2land_quincy, update_land2atm_quincy

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_atmland_interface'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Pass atm forcing
  !
  SUBROUTINE update_atm2land_quincy( &
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
    & nhx_deposition,         &
    & noy_deposition,         &
    & nhx_n15_deposition,     &
    & noy_n15_deposition,     &
    & p_deposition,           &
    ! For lakes:
    & DEBUG_VAR,              &
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
    USE mtime,                     ONLY: datetime
    USE mo_time_config,            ONLY: time_config
    USE mo_jsb_time,               ONLY: get_secs_of_day, get_year, is_newyear
    USE mo_iq_atm2land_process,    ONLY: update_local_time_and_daytime_counter, update_slow_sb_pool_accelerator_bookkeeping
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)             :: tile
    TYPE(t_jsb_task_options),   INTENT(in)                :: options
    REAL(wp), OPTIONAL, DIMENSION(:), INTENT(in)          ::  &
      & t_air, &
      & q_air, &
      & press_air, &
      & rain, &
      & snow, &
      & wind_air, &
      & wind_10m, &
      & lw_srf_down, &
      & swvis_srf_down, &
      & swnir_srf_down, &
      & swpar_srf_down, &
      & fract_par_diffuse, &
      & dz_srf, &
      & press_srf, &
      & rho_srf, &
      & drag_srf, &
      & t_acoef, &
      & t_bcoef, &
      & q_acoef, &
      & q_bcoef, &
      & pch, &
      & cos_zenith_angle, &
      & CO2_air, &
      & nhx_deposition, &
      & noy_deposition, &
      & nhx_n15_deposition, &
      & noy_n15_deposition, &
      & p_deposition, &
      & DEBUG_VAR, &
      & drag_wtr, &
      & drag_ice, &
      & t_acoef_wtr, &
      & t_bcoef_wtr, &
      & q_acoef_wtr, &
      & q_bcoef_wtr, &
      & t_acoef_ice, &
      & t_bcoef_ice, &
      & q_acoef_ice, &
      & q_bcoef_ice
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(SPQ_)

    INTEGER  :: iblk, ics, ice, nc, i
    INTEGER  :: global_seconds_day, current_year
    REAL(wp) :: dtime
    TYPE(datetime),    POINTER :: mtime_current !< elapsed simulation time
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    REAL(wp), POINTER :: lon(:)
    INTEGER :: model_scheme
    LOGICAL :: run_spinup_accelerator                    !< model configuration: if running with slow sb pool spin-up accelerator
    INTEGER :: sb_pool_spinup_accelerator_max_executions !< bookkeeping configurations for the spin-up accelerator: max number of executions
    INTEGER :: sb_pool_spinup_accelerator_frequency      !< bookkeeping configurations for the spin-up accelerator: frequency of executions
    INTEGER :: sb_pool_spinup_accelerator_start_year     !< bookkeeping configurations for the spin-up accelerator: start year of executions

    LOGICAL :: tile_contains_lake

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_atm2land_quincy'

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
      & nhx_deposition_ptr,    &
      & noy_deposition_ptr,    &
      & nhx_n15_deposition_ptr,&
      & noy_n15_deposition_ptr,&
      & p_deposition_ptr,      &
      ! --
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

    dsl4jsb_Real2D_onChunk :: daytime_counter
    dsl4jsb_Real2D_onChunk :: local_time_day_seconds
    dsl4jsb_Real2D_onChunk :: slow_sb_pool_accelerator_execute
    dsl4jsb_Real2D_onChunk :: slow_sb_pool_accelerator_execution_counter
    ! ----------------------------------------------------------------------------------------------------- !
    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    model => Get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    lon   => grid%lon(ics:ice, iblk)
    mtime_current => time_config%tc_current_date
    global_seconds_day = get_secs_of_day(mtime_current)
    current_year = get_year(mtime_current)

    model_scheme = model%config%model_scheme
    run_spinup_accelerator = model%config%flag_slow_sb_pool_spinup_accelerator
    sb_pool_spinup_accelerator_max_executions = model%config%slow_sb_pool_spinup_accelerator_max_executions
    sb_pool_spinup_accelerator_frequency = model%config%slow_sb_pool_spinup_accelerator_frequency
    sb_pool_spinup_accelerator_start_year = model%config%slow_sb_pool_spinup_accelerator_start_year


    ! make this logical accessible in the OpenACC code directly
    tile_contains_lake = tile%contains_lake

    IF (nc /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SPQ_)
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
      drag_srf_ptr          => dsl4jsb_var2D_onChunk(SPQ_, spq_drag_srf)
    END IF
    IF (PRESENT(t_acoef)) THEN
      t_acoef_ptr           => dsl4jsb_var2D_onChunk(SPQ_, spq_t_acoef)
    END IF
    IF (PRESENT(t_bcoef)) THEN
      t_bcoef_ptr           => dsl4jsb_var2D_onChunk(SPQ_, spq_t_bcoef)
    END IF
    IF (PRESENT(q_acoef)) THEN
      q_acoef_ptr           => dsl4jsb_var2D_onChunk(SPQ_, spq_q_acoef)
    END IF
    IF (PRESENT(q_bcoef)) THEN
      q_bcoef_ptr           => dsl4jsb_var2D_onChunk(SPQ_, spq_q_bcoef)
    END IF
    IF (PRESENT(pch)) THEN
      pch_ptr               => dsl4jsb_var2D_onChunk(SPQ_, spq_pch)
    END IF
    IF (PRESENT(cos_zenith_angle)) THEN
      cos_zenith_angle_ptr  => dsl4jsb_var2D_onChunk(A2L_, cos_zenith_angle)
    END IF
    IF (PRESENT(CO2_air)) THEN
      CO2_air_ptr           => dsl4jsb_var2D_onChunk(A2L_, CO2_air)
      CO2_air_mol_ptr       => dsl4jsb_var2D_onChunk(A2L_, CO2_air_mol)
      CO2_mixing_ratio_ptr  => dsl4jsb_var2D_onChunk(A2L_, CO2_mixing_ratio)
    END IF
    IF (PRESENT(nhx_deposition)) THEN
      nhx_deposition_ptr      => dsl4jsb_var2D_onChunk(A2L_, nhx_deposition)
    END IF
    IF (PRESENT(noy_deposition)) THEN
      noy_deposition_ptr      => dsl4jsb_var2D_onChunk(A2L_, noy_deposition)
    END IF
    IF (PRESENT(nhx_n15_deposition)) THEN
      nhx_n15_deposition_ptr  => dsl4jsb_var2D_onChunk(A2L_, nhx_n15_deposition)
    END IF
    IF (PRESENT(noy_n15_deposition)) THEN
      noy_n15_deposition_ptr  => dsl4jsb_var2D_onChunk(A2L_, noy_n15_deposition)
    END IF
    IF (PRESENT(p_deposition)) THEN
      p_deposition_ptr        => dsl4jsb_var2D_onChunk(A2L_, p_deposition)
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

    daytime_counter         => dsl4jsb_var2D_onChunk(A2L_, daytime_counter)
    local_time_day_seconds  => dsl4jsb_var2D_onChunk(A2L_, local_time_day_seconds)
    IF (run_spinup_accelerator) THEN
      slow_sb_pool_accelerator_execute           => dsl4jsb_var2D_onChunk(A2L_, slow_sb_pool_accelerator_execute)
      slow_sb_pool_accelerator_execution_counter => dsl4jsb_var2D_onChunk(A2L_, slow_sb_pool_accelerator_execution_counter)
    END IF
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

    IF (PRESENT(DEBUG_VAR)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        DEBUG_VAR_ptr(i) = DEBUG_VAR(i)
      END DO
    END IF

    IF (PRESENT(t_air)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        t_air_ptr(i) = t_air(i)
      END DO
    END IF
    IF (PRESENT(q_air)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        q_air_ptr(i) = q_air(i)
      END DO
    END IF
    IF (PRESENT(press_air)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        press_air_ptr(i) = press_air(i)
      END DO
    END IF
    IF (PRESENT(rain)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        rain_ptr(i) = rain(i)
      END DO
    END IF
    IF (PRESENT(snow)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        snow_ptr(i) = snow(i)
      END DO
    END IF
    IF (PRESENT(wind_air)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        wind_air_ptr(i) = wind_air(i)
      END DO
    END IF
    IF (PRESENT(wind_10m)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        wind_10m_ptr(i) = wind_10m(i)
      END DO
    END IF
    IF (PRESENT(lw_srf_down)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        lw_srf_down_ptr(i) = lw_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swvis_srf_down)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        swvis_srf_down_ptr(i) = swvis_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swnir_srf_down)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        swnir_srf_down_ptr(i) = swnir_srf_down(i)
      END DO
    END IF
    IF (PRESENT(swpar_srf_down)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        swpar_srf_down_ptr(i) = swpar_srf_down(i)
      END DO
    END IF
    IF (PRESENT(fract_par_diffuse)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        fract_par_diffuse_ptr(i)  = fract_par_diffuse(i)
      END DO
    END IF
    IF (PRESENT(dz_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        dz_srf_ptr(i) = dz_srf(i)
      END DO
    END IF
    IF (PRESENT(press_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        press_srf_ptr(i) = press_srf(i)
      END DO
    END IF
    IF (PRESENT(drag_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        drag_srf_ptr(i) = drag_srf(i)
      END DO
    END IF
    IF (PRESENT(rho_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        rho_srf_ptr(i) = rho_srf(i)
      END DO
    END IF
    IF (PRESENT(t_acoef)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        t_acoef_ptr(i) = t_acoef(i)
      END DO
    END IF
    IF (PRESENT(t_bcoef)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        t_bcoef_ptr(i) = t_bcoef(i)
      END DO
    END IF
    IF (PRESENT(q_acoef)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        q_acoef_ptr(i) = q_acoef(i)
      END DO
    END IF
    IF (PRESENT(q_bcoef)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        q_bcoef_ptr(i) = q_bcoef(i)
      END DO
    END IF
    IF (PRESENT(pch)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        pch_ptr(i) = pch(i)
      END DO
    END IF
    IF (PRESENT(cos_zenith_angle)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        cos_zenith_angle_ptr(i) = cos_zenith_angle(i)
      END DO
    END IF
    IF (PRESENT(CO2_air)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        CO2_air_ptr(i) = CO2_air(i)
        ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
        CO2_air_mol_ptr(i) = CO2_air(i) * molarMassDryAir / molarMassCO2
        ! convert CO2 from "molar ratio (volume)" to "co2 mixing ratio ppmv"
        CO2_mixing_ratio_ptr(i) = CO2_air_mol_ptr(i) * 1000000._wp
      END DO
    END IF
    IF (PRESENT(nhx_deposition)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        nhx_deposition_ptr(i) = nhx_deposition(i)
      END DO
    END IF
    IF (PRESENT(noy_deposition)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        noy_deposition_ptr(i) = noy_deposition(i)
      END DO
    END IF
    IF (PRESENT(nhx_n15_deposition)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        nhx_n15_deposition_ptr(i) = nhx_n15_deposition(i)
      END DO
    END IF
    IF (PRESENT(noy_n15_deposition)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        noy_n15_deposition_ptr(i) = noy_n15_deposition(i)
      END DO
    END IF
    IF (PRESENT(p_deposition)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        p_deposition_ptr(i) = p_deposition(i)
      END DO
    END IF

    !> lakes
    !>
    IF (tile_contains_lake) THEN
      IF (PRESENT(drag_wtr)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          drag_wtr_ptr(i) = drag_wtr(i)
        END DO
      END IF
      IF (PRESENT(drag_ice)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          drag_ice_ptr(i) = drag_ice(i)
        END DO
      END IF
      IF (PRESENT(t_acoef_wtr)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          t_acoef_wtr_ptr(i) = t_acoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(t_bcoef_wtr)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          t_bcoef_wtr_ptr(i) = t_bcoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(q_acoef_wtr)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          q_acoef_wtr_ptr(i) = q_acoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(q_bcoef_wtr)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          q_bcoef_wtr_ptr(i) = q_bcoef_wtr(i)
        END DO
      END IF
      IF (PRESENT(t_acoef_ice)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          t_acoef_ice_ptr(i) = t_acoef_ice(i)
        END DO
      END IF
      IF (PRESENT(t_bcoef_ice)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          t_bcoef_ice_ptr(i) = t_bcoef_ice(i)
        END DO
      END IF
      IF (PRESENT(q_acoef_ice)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          q_acoef_ice_ptr(i) = q_acoef_ice(i)
        END DO
      END IF
      IF (PRESENT(q_bcoef_ice)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=1,nc
          q_bcoef_ice_ptr(i) = q_bcoef_ice(i)
        END DO
      END IF
    END IF

    !$ACC END PARALLEL

    !$ACC WAIT(1)

    IF (.NOT. PRESENT(cos_zenith_angle)) THEN
      CALL finish(TRIM(routine), 'cos_zenith_angle not present but needed for QUINCY model')
    END IF

    ! Update local time and the daytime counter
    CALL update_local_time_and_daytime_counter( &
      &     global_seconds_day, dtime, lon, swpar_srf_down, daytime_counter, local_time_day_seconds)

    ! In case that we run with spin-up acceleration we need to update the accelerator bookkeeping
    IF (run_spinup_accelerator) THEN
      ! in most years the spin-up will not be accelerated
      slow_sb_pool_accelerator_execute(:) = 0.0_wp

      IF(is_newyear(mtime_current, dtime)) THEN
        CALL finish(TRIM(routine), 'This functionality has so far not been tested within a coupled run, but should be!')

        CALL update_slow_sb_pool_accelerator_bookkeeping( dtime, current_year, &
          & sb_pool_spinup_accelerator_max_executions, sb_pool_spinup_accelerator_frequency, sb_pool_spinup_accelerator_start_year, &
          & slow_sb_pool_accelerator_execution_counter, slow_sb_pool_accelerator_execute)
      ENDIF
    END IF
  END SUBROUTINE update_atm2land_quincy

  SUBROUTINE update_land2atm_quincy(tile, options, &
    & t_srf, &
    & t_srf_rad, &
    & t_eff_srf, &
    & qsat_srf, &
    & s_srf, &
    & fact_q_air, &
    & fact_qsat_srf, &
    & evapopot, &
    & evapotrans, &
    & latent_hflx, &
    & sensible_hflx, &
    & grnd_hflx, &
    & grnd_hcap, &
    & rough_h_srf, &
    & rough_m_srf, &
    & q_snocpymlt, &
    & alb_vis_dir, &
    & alb_nir_dir, &
    & alb_vis_dif, &
    & alb_nir_dif, &
    & kh, &
    & km, &
    & kh_neutral, &
    & km_neutral, &
    & CO2_flux, &
    & t_lwtr, &
    & t_lice, &
    & qsat_lwtr, &
    & qsat_lice, &
    & s_lwtr, &
    & s_lice, &
    & evapo_wtr, &
    & latent_hflx_wtr, &
    & sensible_hflx_wtr, &
    & evapo_ice, &
    & latent_hflx_ice, &
    & sensible_hflx_ice, &
    & ice_fract_lake, &
    & alb_vis_dir_wtr, &
    & alb_vis_dif_wtr, &
    & alb_nir_dir_wtr, &
    & alb_nir_dif_wtr, &
    & albedo_lwtr, &
    & albedo_lice)

    USE mo_jsb_physical_constants, ONLY: molarMassCO2

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    REAL(wp), OPTIONAL, INTENT(out) :: &
      & t_srf(:), &
      & t_srf_rad(:), &
      & t_eff_srf(:), &
      & qsat_srf(:), &
      & s_srf(:), &
      & fact_q_air(:), &
      & fact_qsat_srf(:), &
      & evapopot(:), &
      & evapotrans(:), &
      & latent_hflx(:), &
      & sensible_hflx(:), &
      & grnd_hflx(:), &
      & grnd_hcap(:), &
      & rough_h_srf(:), &
      & rough_m_srf(:), &
      & q_snocpymlt(:), &
      & alb_vis_dir(:), &
      & alb_nir_dir(:), &
      & alb_vis_dif(:), &
      & alb_nir_dif(:), &
      & kh(:), &
      & km(:), &
      & kh_neutral(:), &
      & km_neutral(:), &
      & CO2_flux(:), &
      & t_lwtr(:), &
      & t_lice(:), &
      & qsat_lwtr(:), &
      & qsat_lice(:), &
      & s_lwtr(:), &
      & s_lice(:), &
      & evapo_wtr(:), &
      & latent_hflx_wtr(:), &
      & sensible_hflx_wtr(:), &
      & evapo_ice(:), &
      & latent_hflx_ice(:), &
      & sensible_hflx_ice(:), &
      & ice_fract_lake(:), &
      & alb_vis_dir_wtr(:), &
      & alb_vis_dif_wtr(:), &
      & alb_nir_dir_wtr(:), &
      & alb_nir_dif_wtr(:), &
      & albedo_lwtr(:), &
      & albedo_lice(:)

    dsl4jsb_Def_memory(Q_RAD_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(VEG_)

    TYPE(t_jsb_model), POINTER :: model

    INTEGER  :: iblk, ics, ice, nc, i, j, it, nt

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_land2atm_quincy'

    dsl4jsb_Real2D_onChunk :: &
      & t_srf_ptr, &
      ! & t_srf_rad_ptr, &
      & t_eff_srf_ptr, &
      & qsat_srf_ptr, &
      & s_srf_ptr, &
      & fact_q_air_ptr, &
      & fact_qsat_srf_ptr, &
      & evapopot_ptr, &
      & evapotrans_ptr, &
      & latent_hflx_ptr, &
      & sensible_hflx_ptr, &
      & grnd_hflx_ptr
    dsl4jsb_Real3D_onChunk :: &
      & grnd_hcap_ptr
    dsl4jsb_Real2D_onChunk :: &
      & rough_h_srf_ptr, &
      & rough_m_srf_ptr, &
      & q_snocpymlt_ptr, &
      & alb_vis_dir_ptr, &
      & alb_nir_dir_ptr, &
      & alb_vis_dif_ptr, &
      & alb_nir_dif_ptr, &
      ! & kh_ptr, &
      ! & km_ptr, &
      ! & kh_neutral_ptr, &
      ! & km_neutral_ptr, &
      & net_biosphere_production_ptr !, & ! the JSBACH interface uses the CO2_flux here (different unit&sign)
      ! & t_lwtr_ptr, &
      ! & t_lice_ptr, &
      ! & qsat_lwtr_ptr, &
      ! & qsat_lice_ptr, &
      ! & s_lwtr_ptr, &
      ! & s_lice_ptr, &
      ! & evapo_wtr_ptr, &
      ! & latent_hflx_wtr_ptr, &
      ! & sensible_hflx_wtr_ptr, &
      ! & evapo_ice_ptr, &
      ! & latent_hflx_ice_ptr, &
      ! & sensible_hflx_ice_ptr, &
      ! & ice_fract_lake_ptr, &
      ! & alb_vis_dir_wtr_ptr, &
      ! & alb_vis_dif_wtr_ptr, &
      ! & alb_nir_dir_wtr_ptr, &
      ! & alb_nir_dif_wtr_ptr, &
      ! & albedo_lwtr_ptr, &
      ! & albedo_lice

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

    IF (nc /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    model => Get_model(tile%owner_model_id)
    ! use_tmx =  model%config%use_tmx

    dsl4jsb_Get_memory(Q_RAD_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)

    ! set all output var to zero
    !   because at the moment QUINCY does not provide all output var, but only the ones needed without tmx and without lakes
    !   this is a temporary solution until quincy does provide also lake and ice variables
    IF (PRESENT(t_srf)) t_srf(:)                  = 0.0_wp
    IF (PRESENT(t_srf_rad)) t_srf_rad(:)              = 0.0_wp
    IF (PRESENT(t_eff_srf)) t_eff_srf(:)              = 0.0_wp
    IF (PRESENT(qsat_srf)) qsat_srf(:)               = 0.0_wp
    IF (PRESENT(s_srf)) s_srf(:)                  = 0.0_wp
    IF (PRESENT(fact_q_air)) fact_q_air(:)             = 0.0_wp
    IF (PRESENT(fact_qsat_srf)) fact_qsat_srf(:)          = 0.0_wp
    IF (PRESENT(evapopot)) evapopot(:)               = 0.0_wp
    IF (PRESENT(evapotrans)) evapotrans(:)             = 0.0_wp
    IF (PRESENT(latent_hflx)) latent_hflx(:)            = 0.0_wp
    IF (PRESENT(sensible_hflx)) sensible_hflx(:)          = 0.0_wp
    IF (PRESENT(grnd_hflx)) grnd_hflx(:)              = 0.0_wp
    IF (PRESENT(grnd_hcap)) grnd_hcap(:)              = 0.0_wp
    IF (PRESENT(rough_h_srf)) rough_h_srf(:)            = 0.0_wp
    IF (PRESENT(rough_m_srf)) rough_m_srf(:)            = 0.0_wp
    IF (PRESENT(q_snocpymlt)) q_snocpymlt(:)            = 0.0_wp
    IF (PRESENT(alb_vis_dir)) alb_vis_dir(:)            = 0.0_wp
    IF (PRESENT(alb_nir_dir)) alb_nir_dir(:)            = 0.0_wp
    IF (PRESENT(alb_vis_dif)) alb_vis_dif(:)            = 0.0_wp
    IF (PRESENT(alb_nir_dif)) alb_nir_dif(:)            = 0.0_wp
    IF (PRESENT(kh)) kh(:)                     = 0.0_wp
    IF (PRESENT(km)) km(:)                     = 0.0_wp
    IF (PRESENT(kh_neutral)) kh_neutral(:)             = 0.0_wp
    IF (PRESENT(km_neutral)) km_neutral(:)             = 0.0_wp
    IF (PRESENT(CO2_flux)) CO2_flux(:)               = 0.0_wp
    IF (PRESENT(t_lwtr)) t_lwtr(:)                 = 0.0_wp
    IF (PRESENT(t_lice)) t_lice(:)                 = 0.0_wp
    IF (PRESENT(qsat_lwtr)) qsat_lwtr(:)              = 0.0_wp
    IF (PRESENT(qsat_lice)) qsat_lice(:)              = 0.0_wp
    IF (PRESENT(s_lwtr)) s_lwtr(:)                 = 0.0_wp
    IF (PRESENT(s_lice)) s_lice(:)                 = 0.0_wp
    IF (PRESENT(evapo_wtr)) evapo_wtr(:)              = 0.0_wp
    IF (PRESENT(latent_hflx_wtr)) latent_hflx_wtr(:)        = 0.0_wp
    IF (PRESENT(sensible_hflx_wtr)) sensible_hflx_wtr(:)      = 0.0_wp
    IF (PRESENT(evapo_ice)) evapo_ice(:)              = 0.0_wp
    IF (PRESENT(latent_hflx_ice)) latent_hflx_ice(:)        = 0.0_wp
    IF (PRESENT(sensible_hflx_ice)) sensible_hflx_ice(:)      = 0.0_wp
    IF (PRESENT(ice_fract_lake)) ice_fract_lake(:)         = 0.0_wp
    IF (PRESENT(alb_vis_dir_wtr)) alb_vis_dir_wtr(:)        = 0.0_wp
    IF (PRESENT(alb_vis_dif_wtr)) alb_vis_dif_wtr(:)        = 0.0_wp
    IF (PRESENT(alb_nir_dir_wtr)) alb_nir_dir_wtr(:)        = 0.0_wp
    IF (PRESENT(alb_nir_dif_wtr)) alb_nir_dif_wtr(:)        = 0.0_wp
    IF (PRESENT(albedo_lwtr)) albedo_lwtr(:)            = 0.0_wp
    IF (PRESENT(albedo_lice)) albedo_lice(:)            = 0.0_wp


    ! Exchange fields on the land tile
    ! --------------------------------

    t_srf_ptr          => dsl4jsb_var2D_onChunk(SPQ_, t_srf_new)
    t_eff_srf_ptr      => dsl4jsb_var2D_onChunk(SPQ_, temp_srf_eff_4)
    qsat_srf_ptr       => dsl4jsb_var2D_onChunk(SPQ_, qsat_star)
    s_srf_ptr          => dsl4jsb_var2D_onChunk(SPQ_, s_star)
    fact_q_air_ptr     => dsl4jsb_var2D_onChunk(SPQ_, fact_q_air)
    fact_qsat_srf_ptr  => dsl4jsb_var2D_onChunk(SPQ_, fact_qsat_srf)
    evapopot_ptr       => dsl4jsb_var2D_onChunk(SPQ_, evapopot)
    evapotrans_ptr     => dsl4jsb_var2D_onChunk(SPQ_, evapotranspiration)
    latent_hflx_ptr    => dsl4jsb_var2D_onChunk(SPQ_, latent_heat_flx)
    sensible_hflx_ptr  => dsl4jsb_var2D_onChunk(SPQ_, sensible_heat_flx)
    grnd_hflx_ptr      => dsl4jsb_var2D_onChunk(SPQ_, ground_heat_flx)
    grnd_hcap_ptr      => dsl4jsb_var3D_onChunk(SPQ_, heat_capa_sl)           ! 3D var, use first soil layer
    rough_h_srf_ptr    => dsl4jsb_var2D_onChunk(SPQ_, z0h)
    rough_m_srf_ptr    => dsl4jsb_var2D_onChunk(SPQ_, z0m)
    !q_snocpymlt_ptr    => dsl4jsb_var2D_onChunk(SPQ_, ) ! just set to zero, no variable available
    alb_vis_dir_ptr    => dsl4jsb_var2D_onChunk(Q_RAD_, alb_vis)
    alb_nir_dir_ptr    => dsl4jsb_var2D_onChunk(Q_RAD_, alb_nir)
    alb_vis_dif_ptr    => dsl4jsb_var2D_onChunk(Q_RAD_, alb_vis)
    alb_nir_dif_ptr    => dsl4jsb_var2D_onChunk(Q_RAD_, alb_nir)
    ! the JSBACH interface uses the CO2_flux here (different unit and different sign)
    net_biosphere_production_ptr => dsl4jsb_var2D_onChunk(VEG_, net_biosphere_production)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

    IF (PRESENT(t_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        t_srf(i) = t_srf_ptr(i)
      END DO
    END IF

    IF (PRESENT(t_eff_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        t_eff_srf(i) = t_eff_srf_ptr(i) ** 0.25_wp
      END DO
    END IF

    IF (PRESENT(qsat_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        qsat_srf(i) = qsat_srf_ptr(i)
      END DO
    END IF

    IF (PRESENT(s_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        s_srf(i) = s_srf_ptr(i)
      END DO
    END IF

    IF (PRESENT(fact_q_air)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        fact_q_air(i) = fact_q_air_ptr(i)
      END DO
    END IF

    IF (PRESENT(fact_qsat_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        fact_qsat_srf(i) = fact_qsat_srf_ptr(i)
      END DO
    END IF

    IF (PRESENT(evapopot)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        evapopot(i) = evapopot_ptr(i)
      END DO
    END IF

    IF (PRESENT(evapotrans)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        evapotrans(i) = evapotrans_ptr(i)
      END DO
    END IF

    IF (PRESENT(latent_hflx)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        latent_hflx(i) = latent_hflx_ptr(i)
      END DO
    END IF

    IF (PRESENT(sensible_hflx)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        sensible_hflx(i) = sensible_hflx_ptr(i)
      END DO
    END IF

    IF (PRESENT(grnd_hflx)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        grnd_hflx(i) = grnd_hflx_ptr(i)
      END DO
    END IF

    IF (PRESENT(grnd_hcap)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        grnd_hcap(i) = grnd_hcap_ptr(i,1)   ! use value form 1st soil layer
      END DO
    END IF

    IF (PRESENT(rough_h_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        rough_h_srf(i) = rough_h_srf_ptr(i)
      END DO
    END IF

    IF (PRESENT(rough_m_srf)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        rough_m_srf(i) = rough_m_srf_ptr(i)
      END DO
    END IF

    ! IF (PRESENT(q_snocpymlt)) THEN
    !   !$ACC LOOP GANG(STATIC: 1) VECTOR
    !   DO i=1,nc
    !     q_snocpymlt(i) = q_snocpymlt_ptr(i)
    !   END DO
    ! END IF

    IF (PRESENT(alb_vis_dir)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        alb_vis_dir(i) = alb_vis_dir_ptr(i)
      END DO
    END IF

    IF (PRESENT(alb_nir_dir)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        alb_nir_dir(i) = alb_nir_dir_ptr(i)
      END DO
    END IF

    IF (PRESENT(alb_vis_dif)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        alb_vis_dif(i) = alb_vis_dif_ptr(i)
      END DO
    END IF

    IF (PRESENT(alb_nir_dif)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        alb_nir_dif(i) = alb_nir_dif_ptr(i)
      END DO
    END IF

    IF (PRESENT(CO2_flux)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=1,nc
        ! micro-mol CO2 m-2 s-1 -> kg(CO2) m-2 s-1
        ! * -1.0_wp because the atmosphere expects CO2 sources as positive and sinks as negative values
        !   ... differs to the convention for fluxes in ICON which would be "fluxes point downwards"
        !   ... note: in jsbach the negation of the fluxes is already done in the carbon interface
        CO2_flux(i) = net_biosphere_production_ptr(i) * molarMassCO2 * 0.001_wp * 0.000001_wp * (-1.0_wp)
      END DO
    END IF

    !$ACC END PARALLEL
    !$ACC WAIT(1)

  END SUBROUTINE update_land2atm_quincy

#endif
END MODULE mo_q_atmland_interface
