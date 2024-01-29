!> Contains the routines for the surface energy balance on LAND lct_type.
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_seb_land
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, message_text, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_lct_class,      ONLY: GLACIER_TYPE, Contains_lct
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_jsb_time,           ONLY: is_time_experiment_start

  dsl4jsb_Use_processes SEB_, RAD_, TURB_, A2L_, SSE_, HYDRO_
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(SSE_)

  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_energy_land, update_asselin_land, update_surface_fluxes_land

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_land'

CONTAINS
  !
  ! ================================================================================================================================
  !
  SUBROUTINE update_surface_energy_land(tile, options)

    USE mo_sse_process,            ONLY: Get_liquid_max
    USE mo_phy_schemes,            ONLY: qsat_water, qsat_ice, q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, tpfac1, tpfac2, tpfac3, cpd
    USE mo_sse_constants,          ONLY: snow_depth_min
    USE mo_jsb_time,               ONLY: is_newday, timesteps_per_day ! get_asselin_coef, timeStep_in_days
    USE mo_jsb_grid,               ONLY: Get_grid
    USE mo_turb_interface,         ONLY: update_exchange_coefficients

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(A2L_)

    LOGICAL                :: lstart

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: day_temp_min
    dsl4jsb_Real2D_onChunk :: day_temp_max
    dsl4jsb_Real2D_onChunk :: day_temp_sum
    dsl4jsb_Real2D_onChunk :: previous_day_temp_mean
    dsl4jsb_Real2D_onChunk :: previous_day_temp_min
    dsl4jsb_Real2D_onChunk :: previous_day_temp_max
    dsl4jsb_Real2D_onChunk :: F_pseudo_soil_temp
    dsl4jsb_Real2D_onChunk :: N_pseudo_soil_temp
    dsl4jsb_Real2D_onChunk :: pseudo_soil_temp

    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: t_old
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: t_unfilt_old
    dsl4jsb_Real2D_onChunk :: t_eff4
    dsl4jsb_Real2D_onChunk :: qsat_star
    dsl4jsb_Real2D_onChunk :: s_star
    dsl4jsb_Real2D_onChunk :: heat_cap
    dsl4jsb_Real2D_onChunk :: forc_hflx
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: fact_qsat_srf
    dsl4jsb_Real2D_onChunk :: fact_q_air
    dsl4jsb_Real2D_onChunk :: drag_srf      ! old turb
    dsl4jsb_Real2D_onChunk :: ch            ! new turb
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: grnd_hflx
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: w_snow
    dsl4jsb_Real2D_onChunk :: fract_snow
    dsl4jsb_Real2D_onChunk :: snow_soil_dens
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: w_ice_sl
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl
    dsl4jsb_Real3D_onChunk :: matrix_pot_sl
    dsl4jsb_Real3D_onChunk :: bclapp_sl
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: q_acoef
    dsl4jsb_Real2D_onChunk :: q_bcoef

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc)  :: &
     & qsat_srf_old,                    &
     & s_old,                           &
     & dQdT,                            & !< Sensitivity of saturated surface specific humidity to temperature
     & t2s_conv,                        & !< Conversion factor from temperature to dry static energy (C_pd * (1+(delta-1)*q_v))
     & heat_tcoef,                      & !< Heat transfer coefficient (rho*C_h*|v|)
     & t_star,                          &
     & liquid_max
    REAL(wp) :: t_air_in_Celcius

    LOGICAL, DIMENSION(options%nc) :: &
     & is_glacier,                    &
     & is_lake,                       &
     & is_ocean,                      &
     & has_snow_layer

    REAL(wp) :: &
      & w_soil_critical_tmp, &
      & soil_depth_sl_vol_porosity_sl_tmp
    LOGICAL  :: l_snow_tmp, l_freeze_tmp, l_supercool_tmp, use_tmx

    INTEGER  :: iblk, ics, ice, nc, ic, i, iter
    REAL(wp) :: steplen
    INTEGER  :: tmp_timesteps_per_day
    LOGICAL  :: l_is_newday
    INTEGER  :: ilct

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_energy_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    steplen = options%steplen

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    use_tmx = model%config%use_tmx
    l_is_newday = is_newday(options%current_datetime, options%dtime)
    tmp_timesteps_per_day = timesteps_per_day(options%dtime)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_config(SSE_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)

    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)          ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_acoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_bcoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   q_acoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   q_bcoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   press_srf)      ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   q_air)          ! IN

    dsl4jsb_Get_var2D_onChunk(SEB_,   previous_day_temp_mean) ! IN
    dsl4jsb_Get_var2D_onChunk(SEB_,   day_temp_sum)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   day_temp_min)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   day_temp_max)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   previous_day_temp_min)  ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   previous_day_temp_max)  ! OUT

    dsl4jsb_Get_var2D_onChunk(SEB_,   N_pseudo_soil_temp)     ! IN
    dsl4jsb_Get_var2D_onChunk(SEB_,   F_pseudo_soil_temp)     ! IN
    dsl4jsb_Get_var2D_onChunk(SEB_,   pseudo_soil_temp)       ! INOUT

    dsl4jsb_Get_var2D_onChunk(SEB_,   t)              ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_old)          ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)       ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt_old)   ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_eff4)         ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   qsat_star)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   s_star)         ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   heat_cap)       ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   forc_hflx)      ! OUT
    dsl4jsb_Get_var2D_onChunk(TURB_,  fact_q_air)     ! in
    dsl4jsb_Get_var2D_onChunk(TURB_,  fact_qsat_srf)  ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   rad_srf_net)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, w_snow)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)     ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens) ! in
    IF (.NOT. tile%is_glacier) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)   ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, matrix_pot_sl)   ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)       ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl) ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, w_ice_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, w_soil_sl)       ! in
    END IF

    dsl4jsb_Get_var2D_onChunk(SSE_,   grnd_hflx)      ! IN
    dsl4jsb_Get_var2D_onChunk(SSE_,   hcap_grnd)      ! IN

    IF (use_tmx) THEN
      dsl4jsb_Get_var2D_onChunk(TURB_,  ch)           ! in
    ELSE
      dsl4jsb_Get_var2D_onChunk(A2L_,   drag_srf)     ! IN
    END IF

    !$ACC DATA &
    !$ACC   CREATE(qsat_srf_old, s_old, dQdT, t2s_conv, heat_tcoef, t_star, liquid_max) &
    !$ACC   CREATE(is_glacier, is_ocean, is_lake, has_snow_layer)

    ! @todo: currently, fraction of glacier has to be either zero or one!
    IF (use_tmx) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        is_glacier(ic) = tile%is_glacier
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      IF (Contains_lct(tile%lcts, GLACIER_TYPE)) THEN
        DO ilct=1,SIZE(tile%lcts)
          IF (tile%lcts(ilct)%id == GLACIER_TYPE) THEN
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
            DO ic=1,nc
              is_glacier(ic) = tile%lcts(ilct)%fract(ics+ic-1,iblk) > 0._wp
            END DO
            EXIT
          END IF
        END DO
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO ic=1,nc
          is_glacier(ic) = .FALSE.
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END IF

    IF (use_tmx) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        is_lake(ic) = .FALSE.
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(i)
      DO ic=1,nc
        i = ics + ic - 1
        is_lake(ic) = tile%fract(i,iblk) <= 0._wp
      END DO
      !$ACC END PARALLEL LOOP
    END IF
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      is_ocean(ic) = .NOT. grid%lsm(ic,iblk)
    END DO
    !$ACC END PARALLEL LOOP

    IF (is_time_experiment_start(options%current_datetime)) THEN            ! Start of experiment
      lstart = .TRUE.
    ELSE
      lstart = .FALSE.
    END IF

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc

      ! In the moment, the calls of calc_previous_day_variables and calc_pseudo_soil_temp are only
      ! necessary for the phenology process
      t_air_in_Celcius = t_air(ic) - tmelt  ! convert Kelvin in Celcius
      
      ! Sum up of the previous day temperatures
      CALL calc_previous_day_variables(l_is_newday,               & ! Input
                                      lstart,                    & ! Input
                                      t_air_in_Celcius,          & ! Input
                                      !dtime,                     & ! Input
                                      tmp_timesteps_per_day,     & ! Input
                                      previous_day_temp_mean(ic), & ! InOut (for summer- and evergreen)
                                      day_temp_sum(ic),           & ! InOut
                                      previous_day_temp_min(ic),  & ! InOut (for crop)
                                      day_temp_min(ic),           & ! InOut
                                      previous_day_temp_max(ic),  & ! InOut (for crop)
                                      day_temp_max(ic) )            ! InOut

      ! Update of pseudo-soil temperature for each time step
      CALL calc_pseudo_soil_temp(t_air_in_Celcius    , & ! Input
                                N_pseudo_soil_temp(ic), & ! Input
                                F_pseudo_soil_temp(ic), & ! Input
                                pseudo_soil_temp(ic)  )   ! InOut (for summer-, evergreen and crop)

    END DO
    !$ACC END PARALLEL LOOP

    w_soil_critical_tmp = dsl4jsb_Config(SSE_)%w_soil_critical
    l_snow_tmp = dsl4jsb_Config(SSE_)%l_snow
    l_freeze_tmp = dsl4jsb_Config(SSE_)%l_freeze
    l_supercool_tmp = dsl4jsb_Config(SSE_)%l_supercool
                          
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc

      IF (is_ocean(ic) .OR. is_lake(ic)) THEN
        ! Temporary fix for temperatures over complete ocean or lake boxes getting too warm or too cold,
        ! leading to a crash in convect_tables
        t_unfilt_old(ic) = 280._wp
        t_old(ic)        = 280._wp
      ELSE
        ! Save old surface temperature and saturation spec. humidity
        t_unfilt_old(ic) = t_unfilt(ic)
        t_old       (ic) = t(ic)
      END IF

      ! Save old surface saturation specific humidity and compute its sensitivity to temperature
      IF (tile%is_glacier .AND. use_tmx) THEN
        qsat_srf_old(ic) = qsat_ice(t_old(ic), press_srf(ic), use_convect_tables=.NOT. use_tmx)
        dQdT(ic) = ( qsat_ice(t(ic) + 0.001_wp, press_srf(ic), use_convect_tables=.NOT. use_tmx) - qsat_srf_old(ic) ) * 1000._wp
      ELSE
        qsat_srf_old(ic) = qsat_water(t_old(ic), press_srf(ic), use_convect_tables=.NOT. use_tmx)
        dQdT(ic) = ( qsat_water(t(ic) + 0.001_wp, press_srf(ic), use_convect_tables=.NOT. use_tmx) - qsat_srf_old(ic) ) * 1000._wp
      END IF

      ! The heat capacity of the surface layer is currently taken as the heat capacity of the upper soil layer
      heat_cap(ic) = hcap_grnd(ic)
      ! The last term of the rhs of the energy balance equation is the conductive heat flux from the ground below
      forc_hflx(ic) = grnd_hflx(ic)

      ! Old dry static energy
      s_old(ic) = surface_dry_static_energy( &
        & t_old(ic), &
        & q_effective(qsat_srf_old(ic), q_air(ic), fact_qsat_srf(ic), fact_q_air(ic)), cpd)

      ! Transfer coefficient
      IF (use_tmx) THEN
        heat_tcoef(ic) = ch(ic)
      ELSE
        heat_tcoef(ic) = heat_transfer_coef(drag_srf(ic), steplen)
      END IF

      t2s_conv(ic) = s_old(ic) / t_old(ic)

      ! Compute new surface temperature and saturation humidity
      IF (lstart) THEN
        s_star(ic) = s_old(ic)
      ELSE
        CALL surface_temp_implicit(nc,                    & ! in
          MERGE(1._wp, tpfac1, use_tmx), & ! in
          steplen,      & ! in
          t2s_conv(ic),                           & ! in
          t_acoef(ic), t_bcoef(ic), q_acoef(ic), q_bcoef(ic), &  ! in
          s_old(ic), qsat_srf_old(ic), dQdT(ic),    & ! in
          rad_srf_net(ic), forc_hflx(ic), heat_tcoef(ic),    &  ! in
          fact_q_air(ic), fact_qsat_srf(ic),       & ! in
          fract_snow(ic), heat_cap(ic),            & ! in
          & s_star(ic))                             ! out
      END IF

    END DO
    !$ACC END PARALLEL LOOP

    IF (use_tmx .AND. .NOT. lstart) THEN
      DO iter=1,dsl4jsb_Config(SEB_)%niter_tmx
        !$ACC WAIT(1)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO ic=1,nc
          t_unfilt(ic) = s_star(ic) / t2s_conv(ic)
        END DO
        !$ACC END PARALLEL LOOP
        CALL update_exchange_coefficients(tile, options)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO ic=1,nc
          CALL surface_temp_implicit( nc,                    & ! in
            MERGE(1._wp, tpfac1, use_tmx),                   & ! in
            steplen,                                         & ! in
            t2s_conv(ic),                                     & ! in
            t_acoef(ic), t_bcoef(ic), q_acoef(ic), q_bcoef(ic),  & ! in
            s_old(ic), qsat_srf_old(ic), dQdT(ic),              & ! in
            rad_srf_net(ic), forc_hflx(ic), ch(ic),             & ! in
            fact_q_air(ic), fact_qsat_srf(ic),                 & ! in
            fract_snow(ic), heat_cap(ic),                      & ! in
            & s_star(ic))                                       ! out
        END DO
        !$ACC END PARALLEL LOOP
      END DO
    END IF
    !$ACC WAIT(1)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(soil_depth_sl_vol_porosity_sl_tmp)
    DO ic=1,nc

      IF (is_ocean(ic) .OR. is_lake(ic)) THEN
        ! Todo: Temporary fix for temperatures over complete ocean or lake boxes getting too warm or too cold,
        ! leading to a crash in convect_tables (lake and land tiles need to be direct children of root tile!)
        ! s_star(:) = s_old(:)
        s_star(ic) = 280._wp * t2s_conv(ic)
      END IF

      ! New unfiltered surface temperature (\tilde(X)^(t+1)) (see below Eq. 3.3.1.5 in ECHAM3 manual, alpha = 1/tpfac2)
      IF (use_tmx) THEN
        t_unfilt(ic) = s_star(ic) / t2s_conv(ic)
      ELSE
        t_unfilt(ic) = ( tpfac2 * s_star(ic) + tpfac3 * s_old(ic) ) / t2s_conv(ic)
      END IF

      ! Correct new dry static energy for melting of snow on the surface and (if l_freeze=.TRUE.) melting/freezing of
      ! ice/water in the top soil layer that will happen during this time step
      ! @todo: shouldn't this also consider wether there is actually enough energy to melt water/ice or freeze water?
      t_star(ic) = s_star(ic) / t2s_conv(ic)

      IF (l_snow_tmp) THEN
        has_snow_layer(ic) = w_snow(ic) * rhoh2o / snow_soil_dens(ic) > snow_depth_min
      ELSE
        has_snow_layer(ic) = w_snow(ic) > w_soil_critical_tmp
      END IF

      ! Actual snowmelt will be computed later in the hydrology process (update_surface_hydrology)
      IF (has_snow_layer(ic) .OR. is_glacier(ic)) THEN
        ! Snow melt where glacier or snow layer
        t_star(ic) = MIN(t_star(ic), tmelt)
      END IF

      IF (l_freeze_tmp .AND. .NOT. tile%is_glacier) THEN

        IF (l_supercool_tmp) THEN
          soil_depth_sl_vol_porosity_sl_tmp = soil_depth_sl(ic,1) * vol_porosity_sl(ic,1)

          liquid_max(ic) = Get_liquid_max( &
            & t_star(ic)                                   , &
            & soil_depth_sl_vol_porosity_sl_tmp            , & ! Maximum water storage
            & matrix_pot_sl  (ic,1)                        , &
            & bclapp_sl      (ic,1)                          &
            & )
        ELSE
          liquid_max(ic) = 0._wp
        END IF

        IF (.NOT. has_snow_layer(ic) .AND. .NOT. is_glacier(ic)) THEN
          IF (w_soil_sl(ic,1) > MAX(liquid_max(ic), w_soil_critical_tmp)) THEN
            ! Freezing of water in uppermost soil layer where no snow layers above surface and considerable amount of moisture
            t_star(ic) = MAX(t_star(ic), tmelt)
          END IF
          IF (w_ice_sl(ic,1) > w_soil_critical_tmp ) THEN
            ! Ice melt in uppermost soil layer where no snow layers above surface and considerable amount of ice
            t_star(ic) = MIN(t_star(ic), tmelt)
          END IF
        END IF

      END IF

      s_star(ic) = t_star(ic) * t2s_conv(ic)

      ! Compute new surface qsat (after correction of dry static energy for snow melt!)
      IF (lstart) THEN
        qsat_star(ic) = qsat_srf_old(ic)
      ELSE
        qsat_star(ic) = qsat_srf_old(ic) + dQdT(ic) * (s_star(ic) - s_old(ic)) / t2s_conv(ic)
        ! qsat_star(ic) = qsat_srf_old(ic) + dQdT(ic) * (t_star(ic) - t_old(ic))
      END IF

    END DO
    !$ACC END PARALLEL LOOP

#ifndef _OPENACC
    IF (ANY(4._wp * t_star(:) - 3._wp * t_old(:) <= 0._wp .AND. tile%fract(ics:ice,iblk) > 0._wp)) THEN
       i = MINLOC(4._wp * t_star(:) - 3._wp * t_old(:), DIM=1)
       WRITE (message_text,*) 'Instability: Extreme temperature difference from one time step to the next at ',     &
            & '(', grid%lon(i,iblk), ';', grid%lat(i,iblk), '): t_star: ', t_star(i), 'K,  t_old: ' , t_old(i), 'K'
       !CALL finish(TRIM(routine), message_text)
       CALL message(TRIM(routine), message_text)
    END IF
#endif

    ! 'Effective temperature' used in radheat of the atmosphere model (see Eq. 6.3 in ECHAM5 manual)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
    DO ic=1,nc
      t_eff4(ic) = t_old(ic)**3 * (4._wp * t_star(ic) - 3._wp * t_old(ic))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT(1)

    !$ACC END DATA

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_energy_land 

  SUBROUTINE surface_temp_implicit(klon, alpha, steplen,    &
    &            pcp,                                         &
    &            pescoe, pfscoe, peqcoe, pfqcoe,              &
    &            psold, pqsold, pdqsold,                      &
    &            pnetrad, pgrdfl,                             &
    &            pcfh, pcair, pcsat, pfracsu, pgrdcap,        &
    &            psnew)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: &
      !& cpd, &
      !& cpvd1, &
      & zemiss_def, &
      & stbo, &
      & alv,  &
      & als

    INTEGER,  INTENT(in)    :: klon
    REAL(wp), INTENT(in)    :: alpha
    REAL(wp), INTENT(in)    :: steplen
    ! REAL(wp), INTENT(in)    :: pcp(:), pfscoe(:), pescoe(:), pfqcoe(:), peqcoe(:)
    ! REAL(wp), INTENT(in)    :: psold(:), pqsold(:), pdqsold(:)
    ! REAL(wp), INTENT(in)    :: pnetrad(:), pgrdfl(:)
    ! REAL(wp), INTENT(in)    :: pcfh(:), pcair(:), pcsat(:), pfracsu(:)
    ! REAL(wp), INTENT(in)    :: pgrdcap(:)
    ! REAL(wp), INTENT(out)   :: psnew(:)
    REAL(wp), INTENT(in)    :: pcp, pfscoe, pescoe, pfqcoe, peqcoe
    REAL(wp), INTENT(in)    :: psold, pqsold, pdqsold
    REAL(wp), INTENT(in)    :: pnetrad, pgrdfl
    REAL(wp), INTENT(in)    :: pcfh, pcair, pcsat, pfracsu
    REAL(wp), INTENT(in)    :: pgrdcap
    REAL(wp), INTENT(out)   :: psnew
        !
    REAL(wp) :: zcolin, zcohfl, zcoind, zicp, zca, zcs
    REAL(wp) :: pdt
    INTEGER :: nc, ic
    !REAL(wp) :: pc16

    ! nc = SIZE(psold)

    pdt     = alpha * steplen
    !pc16    = cpd * cpvd1

    !$noACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$noACC LOOP GANG VECTOR &
    !$noACC   PRIVATE(zcolin, zcohfl, zcoind, zicp, zca, zcs)
    ! DO ic=1,nc
      ! zicp = 1._wp / pcp(ic)
      ! !
      ! zca = als * pfracsu(ic) +  alv * (pcair(ic) - pfracsu(ic))
      ! zcs = als * pfracsu(ic) +  alv * (pcsat(ic) - pfracsu(ic))
      ! !
      ! zcolin = pgrdcap(ic)*zicp +                                                      &
      !   &      pdt * (zicp * 4._wp * zemiss_def * stbo * ((zicp * psold(ic))**3)  -    &
      !   &      pcfh(ic) * (zca * peqcoe(ic) - zcs) *                                    &
      !   &      zicp * pdqsold(ic))
      ! !
      ! zcohfl = -pdt * pcfh(ic) * (pescoe(ic) - 1._wp)
      ! !
      ! zcoind = pdt * (pnetrad(ic) + pcfh(ic) * pfscoe(ic) +  pcfh(ic) *                              &
      !           &      ((zca * peqcoe(ic) - zcs) * pqsold(ic) + zca * pfqcoe(ic)) + pgrdfl(ic))
      ! !
      ! psnew(ic)  = (zcolin * psold(ic) + zcoind) / (zcolin + zcohfl)
      zicp = 1._wp / pcp
      !
      zca = als * pfracsu +  alv * (pcair - pfracsu)
      zcs = als * pfracsu +  alv * (pcsat - pfracsu)
      !
      zcolin = pgrdcap*zicp +                                                      &
        &      pdt * (zicp * 4._wp * zemiss_def * stbo * ((zicp * psold)**3)  -    &
        &      pcfh * (zca * peqcoe - zcs) *                                    &
        &      zicp * pdqsold)
      !
      zcohfl = -pdt * pcfh * (pescoe - 1._wp)
      !
      zcoind = pdt * (pnetrad + pcfh * pfscoe +  pcfh *                              &
                &      ((zca * peqcoe - zcs) * pqsold + zca * pfqcoe) + pgrdfl)
      !
      psnew  = (zcolin * psold + zcoind) / (zcolin + zcohfl)
    ! END DO
    !$noACC END PARALLEL

  END SUBROUTINE surface_temp_implicit
  !
  ! ================================================================================================================================
  !
  SUBROUTINE update_asselin_land(tile, options)

    USE mo_jsb_time,               ONLY: get_asselin_coef

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    !dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & t,           &
      & t_rad4,      &
      & t_filt,      &
      & t_old,       &
      & t_unfilt,    &
      & t_unfilt_old

    ! Locally allocated variables
    !
    INTEGER :: iblk, ics, ice, nc, ic
    REAL(wp) :: eps
    LOGICAL  :: use_tmx

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_asselin_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    use_tmx = model%config%use_tmx

    ! Get reference to variables for current block
    !
    !dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    dsl4jsb_Get_var2D_onChunk(SEB_,      t)              ! out
    IF (use_tmx) THEN
      dsl4jsb_Get_var2D_onChunk(SEB_,    t_rad4)         ! out
    END IF
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_filt)         ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_old)          ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_unfilt)       ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_unfilt_old)   ! in

    ! Asselin time filter, if applicable
    eps = get_asselin_coef()
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
    DO ic=1,nc
      IF (eps > 0._wp) THEN
        t_filt(ic) = t_unfilt_old(ic) + eps * (t_old(ic) - 2._wp * t_unfilt_old(ic) + t_unfilt(ic))
      ELSE
        t_filt(ic) = t_unfilt(ic)
      END IF
      t(ic) = t_filt(ic)
      IF (use_tmx) THEN
        t_rad4(ic) = t(ic)**4
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_asselin_land

  SUBROUTINE update_surface_fluxes_land(tile, options)

    USE mo_phy_schemes,            ONLY: heat_transfer_coef
    USE mo_jsb_physical_constants, ONLY: alv, als, cpd, cvd

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(A2L_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: s_star
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: ch
    dsl4jsb_Real2D_onChunk :: latent_hflx
    dsl4jsb_Real2D_onChunk :: sensible_hflx
    dsl4jsb_Real2D_onChunk :: latent_hflx_lnd
    dsl4jsb_Real2D_onChunk :: sensible_hflx_lnd
    dsl4jsb_Real2D_onChunk :: evapopot
    dsl4jsb_Real2D_onChunk :: evapotrans
    dsl4jsb_Real2D_onChunk :: fract_snow

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) ::                      &
      & s_air      , &  !< Dry static energy at lowest atmospheric level
      & heat_tcoef      !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc, ic
    REAL(wp) :: steplen
    LOGICAL  :: use_tmx

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_fluxes_land'

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (.NOT. tile%Is_process_active(SEB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    use_tmx = model%config%use_tmx

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)

    dsl4jsb_Get_var2D_onChunk(A2L_,   t_acoef)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_bcoef)             ! in

    dsl4jsb_Get_var2D_onChunk(SEB_,   s_star)              ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)  ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapopot)            ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)          ! in

    dsl4jsb_Get_var2D_onChunk(SEB_,   sensible_hflx)       ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   latent_hflx)         ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   sensible_hflx_lnd)   ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   latent_hflx_lnd)     ! out

    IF (use_tmx) THEN
      dsl4jsb_Get_memory(TURB_)
      dsl4jsb_Get_var2D_onChunk(TURB_, ch)                 ! in
    ELSE
      dsl4jsb_Get_var2D_onChunk(A2L_, drag_srf)            ! in
    END IF

    ! Compute new dry static energy at lowest atmospheric level by back-substitution
    !
    !$ACC DATA CREATE(s_air, heat_tcoef)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc

      s_air(ic) = t_acoef(ic) * s_star   (ic) + t_bcoef(ic)

      IF (use_tmx) THEN
        sensible_hflx(ic) = ch(ic) * (s_air(ic) - s_star(ic)) * cvd / cpd ! Sensible heat flux
      ELSE
        heat_tcoef(ic) = heat_transfer_coef(drag_srf(ic), steplen)    ! Transfer coefficient
        sensible_hflx(ic) = heat_tcoef(ic) * (s_air(ic) - s_star(ic)) ! Sensible heat flux
      END IF

      ! Compute latent heat flux
      !
      latent_hflx(ic) = alv * evapotrans(ic) + (als - alv) * fract_snow(ic) * evapopot(ic)

      ! These two variables are not aggregated together with lake and are the fluxes given back to the atmosphere
      sensible_hflx_lnd(ic) = sensible_hflx(ic)
      latent_hflx_lnd  (ic) = latent_hflx  (ic)

    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
    !$ACC END DATA

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_fluxes_land


    ! Updates the mean day temperature ("previous_day_temp_mean(:)") from air temperature and the mean day NPP-rate
  ! (field previous_day_NPP(:,:)
  SUBROUTINE calc_previous_day_variables( &
    & is_newday,                                    & ! Input
    & l_start,                                      & ! Input
    & t_air_in_Celcius,                             & ! Input
    ! & dtime,                                        & ! Input
    & time_steps_per_day,                           & ! Input
    & previous_day_temp_mean,                       & ! InOut
    & day_temp_sum,                                 & ! InOut
    & previous_day_temp_min,                        & ! InOut
    & day_temp_min,                                 & ! InOut
    & previous_day_temp_max,                        & ! InOut
    &  day_temp_max )                                 ! InOut

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    !  DECLARATIONS

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    LOGICAL,   intent(in)    :: is_newday
    LOGICAL,   intent(in)    :: l_start

    REAL(wp),  intent(in)    :: t_air_in_Celcius!,  & ! air temperature at current time step in lowest atmospheric layer in Celcius
                                !dtime

    INTEGER,  intent(in)     :: time_steps_per_day

    REAL(wp),  intent(inout) :: previous_day_temp_mean, & ! "in" because if it is not calculated new, it should remain
                                previous_day_temp_min,  & ! as before. Without "in" it would be NaN in the output.
                                previous_day_temp_max

    REAL(wp),  intent(inout) :: day_temp_sum,           &
                                day_temp_min,           &
                                day_temp_max
    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES


    !-----------------------------------------------------------------------
    ! CONTENT

    ! --- update mean day values

    ! Note that updating is done globally at the same time step, i.e. for different longitudes the updating happens at different
    ! local times.

    IF (.NOT. l_start .AND. is_newday) THEN  ! day has changed --> recompute mean, min, max day temperature of previous day
                                             ! and reinitialize day_temp_sum(), day_temp_min(), day_temp_max()
       previous_day_temp_mean = day_temp_sum/time_steps_per_day
       day_temp_sum = t_air_in_Celcius

       previous_day_temp_min = day_temp_min
       day_temp_min = t_air_in_Celcius

       previous_day_temp_max = day_temp_max
       day_temp_max = t_air_in_Celcius

    ELSE  ! day has not changed or start of experiment (day_temp_sum is initialized with zero!)
       day_temp_sum = day_temp_sum  + t_air_in_Celcius

       IF (day_temp_min > t_air_in_Celcius)   day_temp_min = t_air_in_Celcius
       IF (day_temp_max < t_air_in_Celcius)   day_temp_max = t_air_in_Celcius
    END IF

  END SUBROUTINE calc_previous_day_variables


  ! --- update_pseudo_soil_temp() --------------------------------------------------------------------------------------------------
  !
  ! This routine computes a weighted running mean of the air temperature, which is interpreted here as a pseudo soil temperature:
  ! (1)   T_ps(t) = N^(-1) * SUM(t'=-infty,t) T(t')*exp(-(t-t')*delta/tau_soil),
  ! where "T(t)" is the air temperature at time step with number "t", "delta" the length of the time step (in days) and
  ! "tau_soil" is the characteristic time for loosing the memory of temperature in the soil (also in days; this is a tuning
  ! parameter! The normalization "N" is
  ! (2)   N = SUM(t'=-infty,t) exp(-(t-t')*delta/tau_soil) = 1/(1 - exp(-delta/tau_soil)).
  ! This normalization constant (called "N_pseudo_soil_temp") is computed during initialization of this phenology module.
  ! Computation of T_ps(t) is performed iteratively: it follows from (1)
  !                    T(t+1)          delta
  ! (3)   T_ps(t+1) =  ------ + exp(- --------) * T_ps(t).
  !                      N            tau_soil
  ! The exponential factor (called F_pseudo_soil_temp) is computed during initialization of this phenology module.
  !
  ! Technically the only effect of this routine is an update of the field "pseudo_soil_temp(:)". The routine has to be called
  ! once every time step for every grid point.
  !
  ! Remark: Instead of air-temperature one could try to use the bottom temperature to compute the pseudo_soil_temp.
  !
  SUBROUTINE calc_pseudo_soil_temp( &
    & t_air_in_Celcius,                       & ! Input
    & N_pseudo_soil_temp,                     & ! Input
    & F_pseudo_soil_temp,                     & ! Input
    & pseudo_soil_temp                        & ! InOut
    & )

    !$ACC ROUTINE SEQ

    REAL(wp), intent(in) :: &
     & t_air_in_Celcius,    &
     & N_pseudo_soil_temp,  &
     & F_pseudo_soil_temp

    REAL(wp),  intent(inout) :: pseudo_soil_temp

    pseudo_soil_temp= t_air_in_Celcius / N_pseudo_soil_temp  +  F_pseudo_soil_temp * pseudo_soil_temp

  END SUBROUTINE calc_pseudo_soil_temp

#endif
END MODULE mo_seb_land
