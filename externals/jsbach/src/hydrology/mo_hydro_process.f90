!> Contains the routines for the hydro processes
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

!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"

MODULE mo_hydro_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish, message_text

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_surface_hydrology_land, calc_surface_hydrology_glacier, calc_soil_hydrology, &
    & calc_wskin_fractions_lice, calc_wet_fractions_veg, calc_wet_fractions_bare,             &
    & get_canopy_conductance, get_water_stress_factor, calc_orographic_features

  INTERFACE get_canopy_conductance
    MODULE PROCEDURE get_canopy_cond_unstressed_simple
    MODULE PROCEDURE get_canopy_cond_stressed_simple
  END INTERFACE get_canopy_conductance

  REAL(wp), PARAMETER ::  zwdmin = 0.05_wp
  !$ACC DECLARE COPYIN(zwdmin)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_process'

CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !! Compute surface hydrology on glacier-free land
  !!
  !! @param [in]     lstart                T: Start of experiment
  !! @param [in]     dtime                 Time step
  !! @param [in]     l_dynsnow             T: Compute density of snow dynamically
  !! @param [in]     snow_depth_max        Snow depth limit [m water equivalent]; -1. for no limit
  !! @param [in]     t_unfilt              Surface tempemperature [K] (unfiltered)
  !! @param [in]     wind_10m              wind speed at 10m height [m/s]
  !! @param [in]     t_air                 lowest layer atmosphere temperature [K]
  !! @param [in]     skinres_canopy_max    Capacity of the canopy skin reservoir (also used to limit snow on canopy) [m]
  !! @param [in]     skinres_max           Total capacity of the skin reservoirs, i.e. soil and canopy [m]
  !! @param [in]     fract_snow            surface snow fraction []
  !! @param [in]     fract_skin            skin reservoir fraction []
  !! @param [in]     hcap_grnd             heat capacity of the uppermost soil layer [J m-2 K-1]
  !! @param [in]     evapotrans            evapotranspiration (including sublimation) [kg m-2 s-1]
  !! @param [in]     evapopot              potential evaporation (if there was enough water/ice) [kg m-2 s-1]
  !! @param [in]     transpiration         transpiration [kg m-2 s-1]
  !! @param [in]     rain                  liquid precipitation [kg m-2 s-1]
  !! @param [in]     snow                  solid precipitation [kg m-2 s-1]
  !! @param [in,out] wtr_skin              water content of the skin reservoir (vegetation and bare soil) [m]
  !! @param [in,out] weq_snow_soil         snow depth at the ground [m water equivalent]
  !! @param [in,out] weq_snow_can          snow depth on the canopy [m water equivalent]
  !! @param [in,out] snow_soil_dens        snow density on soil at non-glacier points [m water equivalent]
  !! @param [out]    q_snocpymlt           Heating due to snow melt on canopy [W m-2]
  !! @param [out]    snow_accum            snow budget change within time step [m water equivalent]
  !! @param [out]    snowmelt_soil         snow/ice melt at land points (excluding canopy) [kg m-2 s-1]
  !! @param [out]    evapotrans_soil       evapotranspiration from soil w/o skin and snow reservoirs [kg m-2 s-1]
  !! @param [out]    evapo_skin            evaporation from skin reservoir [kg m-2 s-1]
  !! @param [out]    evapo_snow            evaporation from snow [kg m-2 s-1]
  !! @param [out]    water_to_soil         Water available for infiltration into the soil [m water equivalent]
  !! @param [out]    evapo_deficit         amount of evaporation extracted from different storage then intented [m]
  !
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE calc_surface_hydrology_land (                          &
    & lstart, dtime, ltpe_closed, l_dynsnow, l_infil_subzero,       &
    & snow_depth_max, steepness, t_soil_sl1, t_unfilt, wind_10m,    &
    & t_air, skinres_canopy_max, skinres_max, weq_pond_max,         &
    & fract_snow, fract_skin, fract_pond,                           &
    & hcap_grnd, evapotrans, evapopot, transpiration, rain, snow,   &
    & weq_rootzone, weq_rootzone_max, weq_soil,                     &
    & wtr_skin, weq_snow_soil, weq_snow_can, snow_soil_dens,        &
    & wtr_pond, ice_pond, q_snocpymlt,                              &
    & snow_accum, snowmelt_soil, pond_freeze, pond_melt,            &
    & evapotrans_soil, evapo_skin, evapo_snow, evapo_pond,          &
    & wtr_pond_net_flx, water_to_soil, evapo_deficit, infilt, runoff)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, alf, dens_snow_min, dens_snow_max, dens_snow
    USE mo_hydro_constants,        ONLY: InterceptionEfficiency
    USE mo_sse_constants,          ONLY: snow_depth_min

    LOGICAL,  INTENT(in) :: lstart, ltpe_closed, l_dynsnow, l_infil_subzero
    REAL(wp), INTENT(in) :: dtime, snow_depth_max
    REAL(wp), INTENT(in) ::                                                 &
      & steepness(:), t_soil_sl1(:), t_unfilt(:), wind_10m(:), t_air(:),    &
      & skinres_canopy_max(:), skinres_max(:), weq_pond_max(:),             &
      & fract_snow(:), fract_skin(:), fract_pond(:), hcap_grnd(:), evapotrans(:),  &
      & evapopot(:), transpiration(:), rain(:), snow(:), weq_rootzone(:),   &
      & weq_rootzone_max(:), weq_soil(:)
    REAL(wp), INTENT(inout) :: &
      & wtr_skin(:), weq_snow_soil(:), weq_snow_can(:), snow_soil_dens(:),  &
      & wtr_pond(:), ice_pond(:)
    REAL(wp), INTENT(out)   ::                                              &
      & q_snocpymlt(:), snow_accum(:), snowmelt_soil(:), pond_freeze(:),    &
      & pond_melt(:), evapotrans_soil(:), evapo_skin(:), evapo_snow(:),     &
      & evapo_pond(:), wtr_pond_net_flx(:), water_to_soil(:),               &
      & evapo_deficit(:), infilt(:), runoff(:)
    !
    !  local variables
    !
    REAL(wp) ::               &
      & rain_in_m,            & !< amount of rainfall within time step [m]
      & new_snow,             & !< amount of snowfall within time step [m]
      & weq_snow_soil_old,    & !< amount of snow prior to update [m]
      & evapotrans_in_m,      & !< amount of evapotranspiration within time step [m]
      & evapotrans_soil_in_m, & !< amount of evapotranspiration from soil (w/o snow, skin and ponds) within time step [m]
      & evapo_soil_in_m,      & !< amount of ground evaporation within time step [m]
      & evapo_skin_in_m,      & !< amount of evaporation from skin reservoir within time step [m]
      & evapo_pond_in_m,      & !< amount of evaporation from pond reservoir within time step [m]
      & transpiration_in_m,   & !< amount of transpiration within time step [m]
      & evapo_snow_in_m,      & !< amount of sublimation from snow within time step [m]
      & snowmelt_can,         & !< snow melt from canopy [m]
      & snowmelt_pot_in_m,    & !< potential amount of snow melt [m]
      & new_snow_can,         & !< amount of snowfall on canopy within time step [m]
      & new_snow_soil,        & !< amount of snowfall to soil within time step [m]
      & exp_t, exp_w,         & !< exponents needed for unloading of snow due to temperature / wind
      & snow_blown,           & !< amount of snow blown from canopy to the ground
      & canopy_pre_snow,      & !< snow depth on the canopy prior to its update
      & evapotrans_no_snow,   & !< evapotranspiration without snow evaporation [m]
      & evapo_snow_pot_in_m,  & !< potential snow evaporation
      & evapo_skin_pot_in_m,  & !< potential evaporation from wet surface (skin reservoir)
      & evapo_pond_pot_in_m,  & !< potential evaporation from pond reservoir
      & rain_to_skinres,      & !< amount of rain going into the skin reservoir
      & wtr_skin_pre_rain,    & !< amount of water in the skin reservoir prior to the update [m]
      & t_pond,               & !< surface pond water temperature [K]
      & wtr_pond_inflow,      & !< pond inflow [m]
      & wtr_pond_outflow,     & !< pond outflow [m]
      & wtr_pond_overflow,    & !< water exceeding the maximum allowed pond water content (liquid + ice) [m]
      & pond_change_pot         !< amount of water that could potentially be frozen or thawed [m]

    INTEGER :: ic, nc
    !
    !  Parameters  - compare chapter 2 "Land Physics" of the JSBACH3 documentation
    !
    REAL(wp), PARAMETER :: zc1 = tmelt - 3._wp ! [K]
    REAL(wp), PARAMETER :: zc2 = 1.87E5_wp     ! [Ks]
    REAL(wp), PARAMETER :: zc3 = 1.56E5_wp     ! [m]
    REAL(wp), PARAMETER :: k_sat_clay = 8.50E-8_wp ! saturated hydraulic conductivity for clay (Terra model) [m s-1]

    nc = SIZE(snowmelt_soil)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR &
    !$ACC   PRIVATE(rain_in_m, new_snow, weq_snow_soil_old)                                  &
    !$ACC   PRIVATE(evapotrans_in_m, evapotrans_soil_in_m, evapo_soil_in_m, evapo_skin_in_m) &
    !$ACC   PRIVATE(transpiration_in_m, evapo_snow_in_m, snowmelt_can, snowmelt_pot_in_m)    &
    !$ACC   PRIVATE(new_snow_can, new_snow_soil, exp_t, exp_w, snow_blown, canopy_pre_snow)  &
    !$ACC   PRIVATE(evapotrans_no_snow, evapo_snow_pot_in_m, evapo_skin_pot_in_m)            &
    !$ACC   PRIVATE(evapo_pond_pot_in_m, rain_to_skinres, wtr_skin_pre_rain)

    DO ic = 1, nc
      snow_accum(ic)        = 0._wp
      weq_snow_soil_old     = weq_snow_soil(ic)
      snowmelt_soil(ic)     = 0._wp
      evapotrans_soil(ic)   = 0._wp
      evapo_skin(ic)        = 0._wp
      evapo_snow(ic)        = 0._wp
      pond_freeze(ic)       = 0._wp
      pond_melt(ic)         = 0._wp
      water_to_soil(ic)     = 0._wp
      evapo_deficit(ic)     = 0._wp
      wtr_pond_net_flx(ic)  = 0._wp

      !-----------------------------------------------------------------
      ! Convert water fluxes to m water equivalent within the time step
      !-----------------------------------------------------------------
      rain_in_m           = rain(ic)          * dtime/rhoh2o
      new_snow            = snow(ic)          * dtime/rhoh2o   !in_m       ! snow_fall
      evapotrans_in_m     = evapotrans(ic)    * dtime/rhoh2o
      transpiration_in_m  = transpiration(ic) * dtime/rhoh2o
      evapo_snow_pot_in_m = fract_snow(ic)    * evapopot(ic) * dtime/rhoh2o
      evapotrans_no_snow  = evapotrans_in_m - evapo_snow_pot_in_m
      evapo_skin_pot_in_m = (1._wp - fract_snow(ic)) * fract_skin(ic) * evapopot(ic) * dtime/rhoh2o
      evapo_pond_pot_in_m = (1._wp - fract_snow(ic)) * fract_pond(ic) * evapopot(ic) * dtime/rhoh2o

      !---------------------------------------------
      !  Budgets of snow (canopy, ground)
      !---------------------------------------------

      ! amount of snow remaining on the canopy
      new_snow_can = MIN(new_snow * InterceptionEfficiency, MAX(skinres_canopy_max(ic) - weq_snow_can(ic), 0._wp))

      ! remaining snow falls on the ground
      new_snow_soil = new_snow - new_snow_can

      ! update snow on the canopy
      ! Note: evaporation happens from canopy, as long as there is enough snow
      canopy_pre_snow  = weq_snow_can(ic)
      weq_snow_can(ic) = MIN(MAX(0._wp, canopy_pre_snow + new_snow_can + evapo_snow_pot_in_m), skinres_canopy_max(ic))
      evapo_snow_in_m  = evapo_snow_pot_in_m - (weq_snow_can(ic) - canopy_pre_snow - new_snow_can)
      ! unloading of snow from the canopy due to melting
      exp_t            = MAX(0._wp, t_air(ic) - zc1) / zc2 * dtime
      snowmelt_can     = weq_snow_can(ic) * (1._wp-EXP(-exp_t))
      weq_snow_can(ic) = weq_snow_can(ic) - snowmelt_can

      ! unloading of snow from the canopy due to wind
      exp_w            = wind_10m(ic) / zc3 * dtime
      snow_blown       = weq_snow_can(ic) * (1._wp-EXP(-exp_w))
      weq_snow_can(ic) = weq_snow_can(ic)  - snow_blown
      new_snow_soil    = new_snow_soil + snow_blown

      ! Heating due to snow melt on canopy
      ! @Todo Why is heating only computed for canopy?
      q_snocpymlt(ic) = snowmelt_can * rhoh2o * alf / dtime

      !  Snowfall and sublimation
      !-----------------------------------------------
      weq_snow_soil(ic) = weq_snow_soil(ic) + new_snow_soil + evapo_snow_in_m

      ! Correction if there was too much snow evaporation
      IF (weq_snow_soil(ic) < 0._wp) THEN
        evapotrans_no_snow = evapotrans_no_snow + weq_snow_soil(ic)
        evapo_deficit(ic)  = weq_snow_soil(ic)
        weq_snow_soil(ic)  = 0._wp
      ELSE
        evapo_deficit(ic)  = 0._wp
      END IF

      !  Snow melt
      !------------------------
      IF (.NOT. lstart) THEN      ! hcap_ground = 0. at model start
        IF (t_unfilt(ic) > tmelt .AND. weq_snow_soil(ic) > 0._wp) THEN
          snowmelt_pot_in_m = hcap_grnd(ic) * (t_unfilt(ic)-tmelt)/(alf*rhoh2o)      ! potential snow and ice melt [m]
          snowmelt_soil(ic) = MAX(MIN(snowmelt_pot_in_m, weq_snow_soil(ic)), 0._wp)  ! snow melt limited by actual snow depth
          weq_snow_soil(ic) = weq_snow_soil(ic) - snowmelt_soil(ic)                  ! reduce snow depth according to melting
          ! surface temperature correction for snow/ice melt is applied in task "snowmelt_correction"
        END IF
      END IF

      !  Snow budget and meltwater
      !-----------------------------------------------------
      wtr_skin(ic)      = wtr_skin(ic) + snowmelt_can                              ! Add melt water from canopy to skin reservoir
      water_to_soil(ic) = snowmelt_soil(ic) + MAX(0._wp, wtr_skin(ic) - skinres_max(ic)) ! Excess water enters the soil, and
      wtr_skin(ic)      = MIN(skinres_max(ic), wtr_skin(ic))                       ! skin reservoir is limited to the maximum
      snow_accum(ic)    = new_snow + evapo_snow_pot_in_m - snowmelt_soil(ic) - snowmelt_can

      IF (snow_depth_max >= 0._wp) THEN
        ! Limit snow depth to avoid grid cells with infinitely growing snow depths in cooler climates
        water_to_soil(ic) = water_to_soil(ic) + MAX(0._wp, weq_snow_soil(ic) - snow_depth_max)
        weq_snow_soil(ic) = MIN(snow_depth_max, weq_snow_soil(ic))
      END IF

      !  Freezing and melting of surface pond reservoir
      ! @todo add t_pond as persistent variable, provide specific heat capacity
      !       and move to sse routine. Also use ice density for ice fluxes
      !------------------------
      ! 1.) Potential melting:
      t_pond = t_unfilt(ic) - snowmelt_soil(ic) / hcap_grnd(ic) * (alf*rhoh2o)
      IF (t_pond > tmelt .AND. ice_pond(ic) > 0._wp) THEN
        pond_change_pot      = hcap_grnd(ic) * (t_pond-tmelt)/(alf*rhoh2o)
        pond_melt(ic)        = MAX(MIN(pond_change_pot, ice_pond(ic)), 0._wp)
        ice_pond(ic)         = ice_pond(ic) - pond_melt(ic)
        wtr_pond(ic)         = wtr_pond(ic) + pond_melt(ic)
        wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) + pond_melt(ic)
      ELSE IF (t_pond < tmelt .AND. wtr_pond(ic) > 0._wp) THEN
        pond_change_pot      = hcap_grnd(ic) * (tmelt-t_pond)/(alf*rhoh2o)
        pond_freeze(ic)      = MAX(MIN(pond_change_pot, wtr_pond(ic)), 0._wp)
        wtr_pond(ic)         = wtr_pond(ic) - pond_freeze(ic)
        ice_pond(ic)         = ice_pond(ic) + pond_freeze(ic)
        wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) - pond_freeze(ic)
      END IF

      !-----------------------------------------------------------
      !   Budget of water in skin and pond reservoir (canopy and ground)
      !-----------------------------------------------------------
      ! Interception of rain (part of rain goes into skin reservoir)
      rain_to_skinres = MIN(rain_in_m * InterceptionEfficiency, MAX(skinres_max(ic) - wtr_skin(ic), 0._wp))
      ! Rest of rain goes to soil
      ! Note: the MAX in following line is only to account for apparent precision error (small negative values)
      water_to_soil(ic) = water_to_soil(ic) + MAX(rain_in_m - rain_to_skinres, 0._wp)
      wtr_skin_pre_rain = wtr_skin(ic)
      wtr_skin(ic)      = MIN(skinres_max(ic), MAX(0._wp, wtr_skin_pre_rain + rain_to_skinres + evapo_skin_pot_in_m))

      ! Note: at this point ignoring the amount of dew that exceeds the skin reservoir
      !       as dew is calculated later based on evapo_soil_in_m (f(evapotrans_soil))
      evapo_skin_in_m = wtr_skin(ic) - (wtr_skin_pre_rain + rain_to_skinres)

      ! Limit pond evaporation to the available pond weq
      evapo_pond_in_m = MAX(-(wtr_pond(ic) + ice_pond(ic)), evapo_pond_pot_in_m)
      ! Only evaporate from ponds if demand cannot be fullfilled from skin reservoir and transpiration
      evapo_pond_in_m = MAX(evapo_pond_in_m, &
        & MIN(0._wp, evapotrans_no_snow - evapo_skin_in_m - transpiration_in_m))

      ! Update evaporation required from the soil
      evapotrans_soil_in_m = evapotrans_no_snow - evapo_skin_in_m - evapo_pond_in_m
      IF (weq_soil(ic) + evapotrans_soil_in_m < 0._wp) THEN
        ! If the soil does not contain enough moisture, evaporate from ponds instead if possible
        evapo_pond_in_m  = MAX(-(wtr_pond(ic) + ice_pond(ic)), &
          &                    evapo_pond_in_m + (weq_soil(ic) + evapotrans_soil_in_m))
        evapotrans_soil_in_m = evapotrans_no_snow - evapo_skin_in_m - evapo_pond_in_m
      END IF
      evapo_soil_in_m = evapotrans_soil_in_m - transpiration_in_m

      ! Add dew to liquid pond storage or evaporate pond water
      IF (evapo_pond_in_m > 0._wp .OR. -evapo_pond_in_m <= wtr_pond(ic)) THEN
        wtr_pond(ic) = wtr_pond(ic) + evapo_pond_in_m
      ELSE IF (-evapo_pond_in_m < wtr_pond(ic) + ice_pond(ic)) THEN
        ice_pond(ic) = ice_pond(ic) + (evapo_pond_in_m + wtr_pond(ic))
        wtr_pond(ic) = 0._wp
      ELSE
        evapo_deficit(ic) = evapo_deficit(ic) + (evapo_pond_in_m + wtr_pond(ic) + ice_pond(ic))
        wtr_pond(ic) = 0._wp
        ice_pond(ic) = 0._wp
      END IF
      wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) + evapo_pond_in_m

      IF (evapo_skin_pot_in_m < 0._wp) THEN
        evapo_deficit(ic) = evapo_deficit(ic) + evapo_skin_pot_in_m - evapo_skin_in_m
      END IF

      ! Positive values of evaporation and transpiration (dew) are added to water_to_soil.
      ! Negative fluxes change soil moisture later in calc_soil_hydrology.
      IF (evapo_soil_in_m > 0._wp)    water_to_soil(ic) = water_to_soil(ic) + evapo_soil_in_m
      IF (transpiration_in_m > 0._wp) water_to_soil(ic) = water_to_soil(ic) + transpiration_in_m

      ! Compute snow density
      IF (l_dynsnow) THEN
        IF (weq_snow_soil_old * rhoh2o / snow_soil_dens(ic) > snow_depth_min) THEN
          ! Update density of old snow following eq. 31 of Verseghy, 1991
          snow_soil_dens(ic) = dens_snow_max - (dens_snow_max - snow_soil_dens(ic)) * EXP(-0.01_wp * dtime/3600._wp)
        ELSE
          snow_soil_dens(ic) = dens_snow_min  ! fresh snow
        END IF
        ! Calculate weighted mean density from old and fresh snow if any compaction happened already
        IF (weq_snow_soil(ic) > weq_snow_soil_old .AND. snow_soil_dens(ic) > dens_snow_min) THEN
          snow_soil_dens(ic) = (snow_soil_dens(ic) * weq_snow_soil_old                    &
            &                  + dens_snow_min * (weq_snow_soil(ic) - weq_snow_soil_old)) &
            &                 / weq_snow_soil(ic)
        END IF
      ELSE
        snow_soil_dens(ic) = dens_snow
      END IF

      !---------------------
      !> Infiltration and surface runoff
      !---------------------
      !> Runoff is either set to zero for Terraplanet setup or is computed via the
      !! ARNO Scheme (default). Note that both, infiltration and runoff, are further
      !! modified in the soil hydrology routine.

      IF (ltpe_closed) THEN
        ! Terra planet setup without runoff
        runoff(ic) = 0._wp
        infilt(ic) = water_to_soil(ic)

      ELSE
        ! Lateral water flows from ponds are connected to the choice of the surface runoff
        ! and infiltration scheme. The following pond flow scheme is based on the WEED
        ! scheme (see https://doi.org/10.5194/tc-15-1097-2021) and is modified to work
        ! with the ARNO runoff scheme.

        ! Compute lateral inflow into ponds, expected to originate
        ! from the whole grid cell and flow towards the depressions. Thus, most of the
        ! available water is expected to end up in ponds
        wtr_pond_inflow = water_to_soil(ic) * fract_pond(ic)**(1.0_wp/3.0_wp)

        ! As ponds form in depressions, no lateral outflow is assumed but
        ! rather infiltration into the soil. Here we use several assumptions:
        ! - ponds usually correspond to clay rich soils --> sat_hyd_con for clay is used as infiltration
        !   limitation and is scaled with pond fraction as only a part of the cell is covered by ponds
        ! - infiltration is only allowed for ponds without ice (used as proxy for frozen ground)
        IF (ice_pond(ic) > 0._wp) THEN
          wtr_pond_outflow = 0._wp
        ELSE
          wtr_pond_outflow = MIN(wtr_pond(ic) + 0.5_wp * wtr_pond_inflow, &
          & k_sat_clay * dtime * fract_pond(ic))
        END IF
        ! Update pond storages
        wtr_pond(ic)         = wtr_pond(ic) + wtr_pond_inflow - wtr_pond_outflow

        ! If pond water (+ ice storage) exceeds the maximum allowed pond content, compute
        ! overflow to be added to water_to_soil
        IF (wtr_pond(ic) + ice_pond(ic) > weq_pond_max(ic)) THEN
          wtr_pond_overflow = MAX(0._wp, wtr_pond(ic) + ice_pond(ic) - weq_pond_max(ic))
          wtr_pond(ic)      = MAX(0._wp, weq_pond_max(ic) - ice_pond(ic))
        ELSE
          wtr_pond_overflow = 0._wp
        END IF
        wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) + wtr_pond_inflow  &
          &                  - wtr_pond_outflow - wtr_pond_overflow

        ! Update amount of water available for runoff and infiltration
        water_to_soil(ic) = water_to_soil(ic) - wtr_pond_inflow + wtr_pond_outflow &
          &              + wtr_pond_overflow

        ! Compute surface runoff and infiltration
        CALL arno_scheme(l_infil_subzero, t_soil_sl1(ic),        &
          &              water_to_soil(ic), weq_rootzone(ic),    &
          &              steepness(ic), weq_rootzone_max(ic),    &
          &              runoff(ic), infilt(ic))

      END IF

      ! Transform fluxes in output from [m] to [kg m-2 s-1]
      snowmelt_soil(ic) = snowmelt_soil(ic) * rhoh2o / dtime
      pond_freeze(ic)   = pond_freeze(ic)   * rhoh2o / dtime
      pond_melt(ic)     = pond_melt(ic)     * rhoh2o / dtime

      evapotrans_soil(ic)  = evapotrans_soil_in_m * rhoh2o / dtime
      evapo_skin(ic)       = evapo_skin_in_m      * rhoh2o / dtime
      evapo_snow(ic)       = evapo_snow_in_m      * rhoh2o / dtime
      evapo_pond(ic)       = evapo_pond_in_m      * rhoh2o / dtime
      wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) * rhoh2o / dtime

    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    !$ACC WAIT(1)

  END SUBROUTINE calc_surface_hydrology_land

  ! ===============================================================================================================================
  !>
  !! calc surface hydrology on glaciers
  !!
  !! @param [in]     lstart                T: Start of experiment
  !! @param [in]     dtime                 Time step
  !! @param [in]     t_unfilt              surface tempemperature [K] (unfiltered)
  !! @param [in]     fract_snow            surface snow fraction []
  !! @param [in]     hcap_grnd             heat capacity of the uppermost soil layer [J m-2 K-1]
  !! @param [in]     evapotrans            evapotranspiration (including sublimation) [kg m-2 s-1]
  !! @param [in]     evapopot              potential evaporation [kg m-2 s-1]
  !! @param [in]     rain                  liquid precipitation [kg m-2 s-1]
  !! @param [in]     snow                  solid precipitation [kg m-2 s-1]
  !! @param [in,out] weq_glac              glacier depth (snow and ice) [m water equivalent]
  !! @param [out]    q_snocpymlt           Heating due to snow melt on canopy [W m-2]
  !! @param [out]    snowmelt              Snow/ice melt at glacier points [kg m-2 s-1]
  !! @param [out]    runoff_glac           glacier runoff (rain+snow/ice melt, no calving) [kg m-2 s-1]
  !! @param [out]    pme_glacier           precipitation minus sublimation on glacier [kg m-2 s-1]
  !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_surface_hydrology_glacier ( &
    & lstart, dtime,                                         &
    & t_unfilt,                                              &
    & fract_snow, hcap_grnd,                                 &
    & evapotrans, evapopot, rain, snow,                      &
    & weq_glac, q_snocpymlt,                                 &
    & snowmelt, runoff_glac,                                 &
    & pme_glacier)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, alf

    LOGICAL,  INTENT(in) :: lstart
    REAL(wp), INTENT(IN) :: dtime
    REAL(wp), INTENT(in) ::       &
      & t_unfilt,                 &
      & fract_snow, hcap_grnd, evapotrans, evapopot, rain, snow
    REAL(wp), INTENT(inout) :: &
      & weq_glac
    REAL(wp), INTENT(out)   :: &
      & q_snocpymlt, snowmelt, runoff_glac, pme_glacier

    !
    !  local variables
    !
    REAL(wp) ::              &
      & rain_in_m,           & !< amount of rainfall within time step [m]
      & new_snow,            & !< amount of snowfall within time step [m]
      & evapotrans_in_m,     & !< amount of evapotranspiration within time step [m]
      & snowmelt_pot_in_m,   & !< potential amount of snow melt [m]
      & evapo_snow_pot_in_m, & !< potential snow evaporation
      & pme_glacier_in_m,    & !< P-E in m
      & runoff_glac_in_m,    & !< Runoff in m
      & snowmelt_in_m          !< Snow/ice melt in m

    !
    !  Parameters  - compare chapter 2 "Land Physics" of the JSBACH3 documentation
    !
    REAL(wp), PARAMETER :: zc1 = tmelt - 3._wp !< [K]
    REAL(wp), PARAMETER :: zc2=1.87E5_wp       !< [Ks]
    REAL(wp), PARAMETER :: zc3=1.56E5_wp       !< [m]
    !$ACC DECLARE COPYIN(zc1, zc2, zc3)

    !----------------------------------------------------------------------------------------------

    ! Convert water fluxes to m water equivalent within the time step
    !-----------------------------------------------------------------

    rain_in_m          = rain       * dtime/rhoh2o
    new_snow           = snow       * dtime/rhoh2o   !in_m       ! snow_fall
    evapotrans_in_m    = evapotrans * dtime/rhoh2o
    evapo_snow_pot_in_m = fract_snow * evapopot * dtime/rhoh2o

    !  Snowfall and sublimation on glaciers
    !---------------------------------------

    weq_glac         = weq_glac  + new_snow + evapo_snow_pot_in_m   ! glacier depth [m]
    pme_glacier_in_m = rain_in_m + new_snow + evapotrans_in_m       ! P-E on glaciers [m]
    runoff_glac_in_m = rain_in_m                                    ! there is no infiltration on glaciers


    !  Snow and glacier melt
    !------------------------

    snowmelt_in_m = 0._wp
    IF (.NOT. lstart) THEN      ! hcap_ground = 0. at model start
      IF (t_unfilt > tmelt) THEN
        snowmelt_pot_in_m = hcap_grnd * (t_unfilt-tmelt)/(alf*rhoh2o) !  potential snow and ice melt [m]
        snowmelt_in_m     = snowmelt_pot_in_m                       ! there is an unlimited amount of snow
        weq_glac          = weq_glac - snowmelt_in_m                ! reduce glacier depth according to melting
        runoff_glac_in_m  = runoff_glac_in_m + snowmelt_in_m        ! add melt water to the runoff
        ! surface temperature correction for snow/ice melt is applied in task "snowmelt_correction"
      END IF
    END IF

    q_snocpymlt = 0._wp

    snowmelt    = snowmelt_in_m    * rhoh2o / dtime
    runoff_glac = runoff_glac_in_m * rhoh2o / dtime
    pme_glacier = pme_glacier_in_m * rhoh2o / dtime

  END SUBROUTINE calc_surface_hydrology_glacier

  ! ===============================================================================================================================
  !>
  !! Calculates vertical water movement through the soil column
  !
  ! @param [in]      nc                  Vector lenth
  ! @param [in]      nsoil               Number of soil layers
  ! @param [in]      dtime               Time step [s]
  ! @param [in]      l_infil_subzero     Allow infiltration at soil temperatures below 0degC
  ! @param [in]      soil_depth_sl       Thicknesses of soil layers until bedrocj [m]
  ! @param [in]      root_depth_sl       Thicknesses of soil layers until rooting depth [m]
  ! @param [in]      hyd_cond_sat_sl     Hydraulic conductivity of saturated soil [m/s]
  ! @param [in]      matrix_pot_sl       Soil matrix potential [m]
  ! @param [in]      bclapp_sl           Exponent B in Clapp and Hornberger
  ! @param [in]      pore_size_index_sl  Soil pore size distribution index used in Van Genuchten method
  ! @param [in]      vol_porosity_sl     Soil porosity                      - ice not considered [m/m]
  ! @param [in]      vol_field_cap_sl    Volumetric field capacity          - ice not considered [m/m]
  ! @param [in]      vol_p_wilt_sl       Volumetric permanent wilting point - ice not considered [m/m]
  ! @param [in]      evapotrans_soil     Evapotranspiration from the soil (w/o snow or canopy) [kg m-2 s-1]
  ! @param [in]      transpiration       Evapotranspiration [kg m-2 s-1]
  ! @param [in]      wtr_soil_pot_scool_sl Amount of potentially supercooled water [m]
  ! @param [in,out]  wtr_soil_sl         Content un-frozen water of soil layers in rooting zone  [m]
  ! @param [in,out]  ice_soil_sl         Content of frozen water of soil layers in rooting zone  [m]
  ! @param [out]     wsat_sl             Maximum possible amount of water in soil layer          [m]
  ! @param [out]     field_cap_sl        Water content at field capacity in soil layer           [m]
  ! @param [out]     p_wilt_sl           Water content at permanent wilting point in soil layers [m]
  ! @param [out]     runoff              Runoff
  ! @param [out]     drainage            Drainage
  !
  SUBROUTINE calc_soil_hydrology(                                                        &
    ! in
    & nc, l_fract, lat, lon, nsoil, dtime, ltpe_closed, ltpe_open,                       &
    & soilhydmodel, interpol_mean,                                                       &
    & w_soil_wilt_fract, soil_depth_sl, root_depth_sl,                                   &
    & hyd_cond_sat_sl, matrix_pot_sl, bclapp_sl, pore_size_index_sl, vol_porosity_sl,    &
    & vol_field_cap_sl, vol_p_wilt_sl, wtr_soil_pot_scool_sl,                            &
    & fract_pond_max, evapotrans_soil, transpiration, ice_pond, weq_pond_max,            &
    ! inout
    & infilt, runoff, wtr_soil_sl, ice_soil_sl, wtr_pond, wtr_pond_net_flx,              &
    & tpe_overflow, evapo_deficit,                                                       &
    ! out
    & wtr_wsat_sl, wtr_field_cap_sl, wtr_p_wilt_sl, infilt_over, drainage, drain_rock,   &
      drain_lat, wtr_soilhyd_res                                                         &
    & )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    INTEGER,  INTENT(in) :: &
      & nc, nsoil
    INTEGER, INTENT(in) :: &
      & soilhydmodel, interpol_mean
    REAL(wp), INTENT(in) :: &
      & dtime, w_soil_wilt_fract
    LOGICAL, INTENT(in) :: ltpe_closed
    LOGICAL, INTENT(in) :: ltpe_open
    LOGICAL, INTENT(in), DIMENSION(:) :: l_fract      ! Tile fraction > 0
    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & soil_depth_sl, root_depth_sl, hyd_cond_sat_sl, matrix_pot_sl, bclapp_sl, pore_size_index_sl, &
      & vol_porosity_sl, vol_field_cap_sl, vol_p_wilt_sl, wtr_soil_pot_scool_sl
    REAL(wp), INTENT(in), DIMENSION(:) ::               &
      & lat, lon,                                       &
      & fract_pond_max, evapotrans_soil, transpiration, &
      & ice_pond, weq_pond_max
    REAL(wp), INTENT(inout), DIMENSION(:) :: &
      & tpe_overflow, infilt, runoff, evapo_deficit, wtr_pond, wtr_pond_net_flx
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: &
      & wtr_soil_sl, ice_soil_sl
    REAL(wp), INTENT(out), DIMENSION(:) :: &
      & infilt_over, drainage,             &
      & drain_rock, drain_lat, wtr_soilhyd_res
    REAL(wp), INTENT(out), DIMENSION(:,:) :: &
      & wtr_wsat_sl, wtr_field_cap_sl, wtr_p_wilt_sl

    ! Local variables
    REAL(wp) ::                       &
      & evapo_soil_in_m(nc),          & !< evaporation from the soil (without transpiration) [m]
      & transpiration_in_m(nc),       & !< transpiration within time step [m]
      & wtr_pond_corr(nc)               !< correction term for ponds fluxes due to saturated soils [m]

    REAL(wp) ::              &
      & drain_lat_sl(nc,nsoil)

    INTEGER :: ic, is

    INTEGER, PARAMETER :: ilog = 0 ! Switch for debugging output
    !$ACC DECLARE COPYIN(ilog)

    !$ACC DATA &
    !$ACC   CREATE(evapo_soil_in_m, transpiration_in_m, wtr_pond_corr, drain_lat_sl)

    !
    !> Calculation of soil saturation, field capacity and wilting point
    !

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc
        ! Note: here we are calculating soil saturation (this does not correspond to weq_rootzone_max, which was
        !    calculated from field capacity)
        wtr_wsat_sl     (ic,is) = MAX(vol_porosity_sl (ic,is) * soil_depth_sl(ic,is) - ice_soil_sl(ic,is), 0._wp)  ! -> diffusion
        wtr_field_cap_sl(ic,is) = MAX(vol_field_cap_sl(ic,is) * soil_depth_sl(ic,is) - ice_soil_sl(ic,is), 0._wp)  ! -> percolation
        ! @todo: to be discussed: Does it make sense, to reduce water content at wilting point relative to the ice fraction?
        IF (vol_porosity_sl(ic,is) > 0._wp) THEN
          wtr_p_wilt_sl(ic,is) = wtr_wsat_sl(ic,is) * vol_p_wilt_sl(ic,is) / vol_porosity_sl(ic,is)
        ELSE
          wtr_p_wilt_sl(ic,is) = 0._wp
        END IF
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      ! Initialize diagnostic runoff and drainage variables
      infilt_over(ic)    = 0._wp
      drain_rock(ic)     = 0._wp
      ! Note: any dew is removed from these fluxes, as dew is already applied
      !       during the surface hydrology computation
      evapo_soil_in_m(ic)    = MIN(0._wp, evapotrans_soil(ic) - transpiration(ic)) * dtime / rhoh2o
      transpiration_in_m(ic) = MIN(0._wp, transpiration(ic))  * dtime / rhoh2o

      !---------------------
      !> 2. Vertical water transport within the soil column
      !---------------------

      ! Initilalization of drainage
      drain_lat(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc
        drain_lat_sl(ic,is) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    CALL soilhyd(nc, l_fract, lat, lon, nsoil, dtime, soilhydmodel, interpol_mean,   &
      &     w_soil_wilt_fract, soil_depth_sl, wtr_field_cap_sl, wtr_wsat_sl,         &
      &     hyd_cond_sat_sl, vol_porosity_sl, vol_field_cap_sl, vol_p_wilt_sl,       &
      &     bclapp_sl, matrix_pot_sl, pore_size_index_sl,                            &
      &     wtr_soil_pot_scool_sl, ice_soil_sl, wtr_soil_sl,                         &
      &     transpiration_in_m, evapo_soil_in_m, root_depth_sl,                      &
      &     infilt, infilt_over, drain_rock, drain_lat_sl,                           &
      &     tpe_overflow, evapo_deficit, wtr_soilhyd_res, ltpe_closed, ltpe_open)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

    ! Aggregate or distribute the different runoff flux components
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        drain_lat(ic) = drain_lat(ic) + drain_lat_sl(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      drainage(ic)  = drain_rock(ic) + drain_lat(ic)

      ! Infiltration overflow goes into ponds water storage
      IF (fract_pond_max(ic) > EPSILON(1._wp) .AND. infilt_over(ic) > 0._wp) THEN
        ! Pond state need to be corrected accordingly
        ! Add infiltration overflow to pond storage
        wtr_pond_corr(ic)    = MIN(MAX(0._wp, weq_pond_max(ic) - wtr_pond(ic) - ice_pond(ic)), infilt_over(ic))
        wtr_pond(ic)         = wtr_pond(ic) + wtr_pond_corr(ic)
        wtr_pond_net_flx(ic) = wtr_pond_net_flx(ic) + wtr_pond_corr(ic) * rhoh2o / dtime
        ! Check for remaining overflow and add to surface runoff
        wtr_pond_corr(ic)    = MAX(0._wp, infilt_over(ic) - wtr_pond_corr(ic))
        runoff(ic)           = runoff(ic) + wtr_pond_corr(ic)

      ! Without ponds, all infitration overflow goes directly into surface runoff
      ELSE
        runoff(ic)    = runoff(ic) + infilt_over(ic)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      infilt(ic)        = infilt       (ic) * rhoh2o / dtime
      runoff(ic)        = runoff       (ic) * rhoh2o / dtime
      infilt_over(ic)   = infilt_over  (ic) * rhoh2o / dtime
      drainage(ic)      = drainage     (ic) * rhoh2o / dtime
      drain_rock(ic)    = drain_rock   (ic) * rhoh2o / dtime
      drain_lat(ic)     = drain_lat    (ic) * rhoh2o / dtime
      evapo_deficit(ic) = evapo_deficit(ic) * rhoh2o / dtime
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE calc_soil_hydrology


  !------------------------------------------------------------------------------------------------
  !! @brief   Calculation of subgrid orographic features
  !!
  !! @par Revision History
  !!   Written by Tobias Stacke (MPI-M) in January 2023 based on the
  !!   ARNO Scheme infiltration routine

  SUBROUTINE calc_orographic_features(nc, nlat, oro_stddev, steepness)

    USE mo_hydro_constants,        ONLY: oro_var_min, oro_var_max

    INTEGER,  INTENT(in) :: &
      & nc, nlat

    REAL(wp), INTENT(in)  :: oro_stddev(:)  !< standard deviation of orography                   [m]
    REAL(wp), INTENT(out) :: steepness(:)   !< parameter defining the subgrid slope distribution [/]

    INTEGER :: ic

    REAL(wp) ::                       &
      & sigma_0,                      & !< minimum value (100 m); below b = 0.01
      & sigma_max                       !< resolution-dependent maximum
    !
    ! steepness b: shape parameter, defining the sub-gid scale steepness of the orographie:
    !
    ! b = (oro_stddev - sigma_0) / (oro_stddev + sigma_max)
    !
    sigma_0   = oro_var_min
    sigma_max = oro_var_max * 64._wp / REAL(nlat, wp) ! @todo: site level

    ! Compute steepness parameter

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      steepness(ic) = MAX(0._wp, oro_stddev(ic) - sigma_0) / (oro_stddev(ic) + sigma_max)
      ! Limit parameter to realistic bounds
      steepness(ic) = MAX(MIN(steepness(ic), 0.5_wp), 0.01_wp)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE calc_orographic_features

  !------------------------------------------------------------------------------------------------
  !! @brief   Calculation of surface runoff based on the ARNO Scheme
  !! (E. Todini, The ARNO rainfall-runoff model, doi:10.1016/S0022-1694(96)80016-3.)
  !!
  !! @par Revision History
  !!   Revised by Tobias Stacke (MPI-M) in January 2023

#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE arno_scheme(l_infil_subzero, t_soil_sl1, water_to_soil, weq_rootzone, &
    &                    steepness, weq_rootzone_max,                              &
    &                    runoff, infilt)

  !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: tmelt

    LOGICAL, INTENT(in)  :: l_infil_subzero

    REAL(wp), INTENT(in) ::        &
      & t_soil_sl1,                & !< temperature of uppermost soil layer               [K]
      & water_to_soil,             & !< amount of water going into the soil               [m]
      & weq_rootzone,              & !< liquid water + ice within the root zone           [m]
      & steepness,                 & !< parameter defining the subgrid slope distribution [/]
      & weq_rootzone_max             !< maximum water holding capacity (liquid + ice) of the root zone [m]

    REAL(wp), INTENT(out) ::       &
      & runoff,                    & !< surface runoff                                    [m]
      & infilt                       !< infiltration of water into the soil               [m]

    ! Local variables
    REAL(wp) ::                    &
      & ws_min_drain,              & !< minimum amount of root zone soil water for drainage
      & ws_rel,                    & !< relative rootzone soil moisture (water + ice)
      & zb1, zbm, zconw1, zvol

    ! todo: parameters should be defined as parameters
    zb1    = 1._wp + steepness
    zbm    = 1._wp / zb1
    zconw1 = weq_rootzone_max * zb1

    !>  Surface runoff and infiltration
    !   -----------------------------------
    !   f(w) = 1 - (1 - w/w_max)**b

    ! No drainage is possible below a certain minimum amount of water (zwdmin). This corresponds
    ! to the layered bucket scheme, and has nothing to do with field capacity or wilting point.
    ws_min_drain = zwdmin * weq_rootzone_max

    ! Temperature dependence of infiltration should not be included when using the 5-layer snow scheme because:
    !     when tsurf == 0.C and T_soil < 0.C --> all snow melt is going to surface runoff
    IF (.NOT. l_infil_subzero .AND. (t_soil_sl1 < tmelt)) THEN
      runoff = water_to_soil       !  no infiltration as uppermost soil layer is frozen -> runoff
      infilt = 0._wp

    ELSE
      IF (water_to_soil > 0._wp .AND. weq_rootzone > ws_min_drain) THEN    ! soilwater content above limit for drainage
        ! compare jsbach3 documentation section 2.3.2.2
        ws_rel = MIN(1._wp, weq_rootzone / weq_rootzone_max)               ! relative rootzone soil moisture
        zvol   = (1._wp - ws_rel)**zbm - water_to_soil / zconw1            ! factor taking account subgrid scale
                                                                           ! inhomogeneity and steepness

        runoff = water_to_soil - (weq_rootzone_max - weq_rootzone)         ! > 0: water exceeding soil capacity
        IF (zvol > 0._wp) runoff = runoff + weq_rootzone_max * zvol**zb1
        runoff = MAX(MIN(runoff, water_to_soil), 0._wp)                    ! no runoff < 0, limited by available water
        infilt = water_to_soil - runoff                                    ! water not going into runoff is infiltrated

      ELSE
        runoff = 0._wp                                                     ! no runoff as soil is too dry
        infilt = water_to_soil                                             ! all water infiltrated

      END IF
    END IF

  END SUBROUTINE arno_scheme

  !------------------------------------------------------------------------------------------------
  !! @brief   Calculation of vertical soil water movement for multiple soil layers
  !!
  !! @par Revision History
  !!   Routine based on the soil hydrology proposed by Christian Reick (MPI-M) and written
  !!   by Philipp de Vrese (MPI-M)
  !!   Revised by Tobias Stacke (MPI-M) in September 2022

  SUBROUTINE soilhyd(nc, l_fract, lat, lon, nsoil, dtime,          &
    & soilhydmodel, interpol_mean,                                 &
    & w_soil_wilt_fract, dsoil,                                    &
    & wtr_field_cap_sl, wtr_wsat_sl, hyd_cond_sat_sl,              &
    & vol_porosity_sl, vol_field_cap_sl, vol_p_wilt_sl,            &
    & bclapp_sl, matrix_pot_sl, pore_size_index_sl,                &
    & wtr_soil_pot_scool_sl, ice_soil_sl, wtr_soil_sl,             &
    & transpiration_in_m, evapo_soil_in_m, root_depth_sl,          &
    & infilt, infilt_over, drain_rock, drain_lat,                  &
    & tpe_overflow, evapo_deficit, wtr_soilhyd_res,                &
    & ltpe_closed, ltpe_open)

    INTEGER,  INTENT(IN)            :: nc                         ! vector length
    LOGICAL,  INTENT(IN)            :: l_fract(:)                 ! Tile fraction > 0
    REAL(wp), INTENT(IN)            :: lat(:)                     ! Grid cell latitudes
    REAL(wp), INTENT(IN)            :: lon(:)                     ! Grid cell longitudes
    INTEGER,  INTENT(IN)            :: nsoil                      ! number of soil layers
    INTEGER,  INTENT(in)            :: soilhydmodel               ! hydrology model scheme
    INTEGER,  INTENT(in)            :: interpol_mean              ! lateral interpolation scheme
    REAL(wp), INTENT(IN)            :: dtime                      ! model time step legth [s]
    REAL(wp), INTENT(in)            :: w_soil_wilt_fract       ! fractional plant wilting point [/]
    REAL(wp), INTENT(IN)            :: dsoil(:,:)                 ! soil depth until bedrock within each layer [m]
    REAL(wp), INTENT(IN)            :: wtr_field_cap_sl(:,:)      ! field capacity of the soil layer (reduced by ice) [m]
    REAL(wp), INTENT(IN)            :: wtr_wsat_sl(:,:)           ! saturation capacity of the soil layer (reduced by ice) [m]
    REAL(wp), INTENT(IN)            :: hyd_cond_sat_sl(:,:)       ! hydraulic conductivity of saturated soil [m/s]
    REAL(wp), INTENT(IN)            :: vol_porosity_sl(:,:)       ! Soil porosity                      - ice not considered [m/m]
    REAL(wp), INTENT(IN)            :: vol_field_cap_sl(:,:)      ! Volumetric field capacity          - ice not considered [m/m]
    REAL(wp), INTENT(IN)            :: vol_p_wilt_sl(:,:)         ! Volumetric permanent wilting point - ice not considered [m/m]
    REAL(wp), INTENT(IN)            :: bclapp_sl(:,:)             ! exponent B in Clapp and Hornberger
    REAL(wp), INTENT(IN)            :: matrix_pot_sl(:,:)         ! matrix potential [m]
    REAL(wp), INTENT(IN)            :: pore_size_index_sl(:,:)    ! soil pore size distribution index used in Van Genuchten
    REAL(wp), INTENT(IN)            :: wtr_soil_pot_scool_sl(:,:) ! potentially supercooled water [m]
    REAL(wp), INTENT(IN)            :: root_depth_sl(:,:)         ! thicknesses of soil layers until rooting depth [m]
    REAL(wp), INTENT(IN)            :: transpiration_in_m(:)      ! amount of transpiration within timestep [m]
    REAL(wp), INTENT(IN)            :: evapo_soil_in_m(:)         ! amount of soil evaporation within timestep [m]
    REAL(wp), INTENT(INOUT)         :: infilt(:)                  ! amount of infiltration within timestep [m]
    REAL(wp), INTENT(INOUT)         :: ice_soil_sl(:,:)           ! amount of ice in the soil [m]
    REAL(wp), INTENT(INOUT)         :: wtr_soil_sl(:,:)           ! soil moisture of each layer [m]
    REAL(wp), INTENT(INOUT)         :: infilt_over(:)             ! amount of infiltration overflow (Horton + Dunne) [m]
    REAL(wp), INTENT(INOUT)         :: drain_rock(:)              ! drainage from lowest soil layer [m]
    REAL(wp), INTENT(INOUT)         :: drain_lat(:,:)             ! overflow from subsurface soil layers [m]
    REAL(wp), INTENT(INOUT)         :: tpe_overflow(:)
    REAL(wp), INTENT(INOUT)         :: evapo_deficit(:)           ! Water evaporated from unintended sources [m]
    REAL(wp), INTENT(OUT)           :: wtr_soilhyd_res(:)         ! Residual error of the vertical transport scheme

    LOGICAL, INTENT(IN) :: ltpe_closed
    LOGICAL, INTENT(IN) :: ltpe_open

    ! Local variables
    INTEGER  :: ic, is                       ! looping index (cells index and soil layer))
    INTEGER  :: last_soil_layer(nc)          ! deepest soil layer index (above the bedrock boundary)
    CHARACTER(len=4096) :: message_text_long ! long string for soil hydrology checks

    REAL(wp) :: weq_soil(nc)                 ! total column soil moisture (water + ice)
    REAL(wp) :: ws_inter_sl(nc,nsoil)        ! soil moisture with infiltration added to first layer
    REAL(wp) :: ws_vol_sl(nc,nsoil)          ! volumetric soil mositure
    REAL(wp) :: ice_impeda_sl(nc,nsoil)      ! ice impedance factor on all soil layers
    REAL(wp) :: hyd_cond_sl(nc,nsoil)        ! hydraulic conductivity of a layer [m/s]
    REAL(wp) :: diffus_sl(nc,nsoil)          ! diffusivity at mid-depth of the layer[m2/s]
    REAL(wp) :: hyd_cond_li(nc,nsoil)        ! hydraulic conductivity at upper layer interface
    REAL(wp) :: diffus_li(nc,nsoil)          ! diffusivity at upper layer interface
    REAL(wp) :: drain_sl(nc)                 ! drainage from the current layer
    REAL(wp) :: evapo_soil_sl1(nc)           ! bare soil evaporation extracted from top layer (mobile part)
    REAL(wp) :: transpiration_sl(nc,nsoil)   ! transpiration extracted from every soil layer
    REAL(wp) :: deficit_evapotrans(nc)       ! water that could not be extracted from the soil [m]
    REAL(wp) :: deficit_trans(nc)            ! transpiration that could not be extracted from the soil [m]
    REAL(wp) :: deficit_sevap(nc)            ! soil evaporation that could not be extracted from the soil [m]
    REAL(wp) :: remoist(nc)                  ! for TPE
    REAL(wp) :: wtr_rsdl_sl(nc,nsoil)        ! residual liquid water that can not move through the soil (supercooled or < pwilt_sl)
    REAL(wp) :: weq_wsat_sl(nc,nsoil)        ! Maximum storage in each soil layer (porosity) (not reduced by ice) [m]
    REAL(wp) :: weq_field_cap_sl(nc,nsoil)   ! field capacity of the soil layer (not reduced by ice) [m]
    REAL(wp) :: weq_p_wilt_sl(nc,nsoil)      ! wilting point - not reduced relative to soil ice fraction [m]
    REAL(wp) :: wtr_transp_residual(nc)      ! numerical error caused by the vertical transport scheme [m]
    REAL(wp) :: wtr_transp_corr(nc)          ! uncorrected part of the transport residual [m]
    REAL(wp) :: wtr_exinfcap(nc)             ! surface runoff due to exceeding the soil infiltration capacity [m]
    REAL(wp) :: wtr_exwsat(nc)               ! surface runoff due to exceeding the soil moisture saturation [m]
    REAL(wp) :: wtr_soil_mobile_sl(nc,nsoil) ! mobile soil water (needed with ACC)
    REAL(wp) :: hlp1(nc,nsoil), hlp4(nc,nsoil)
    REAL(wp) :: hlp2(nc), hlp3(nc), hlp5

    !$ACC DATA &
    !$ACC   CREATE(weq_soil, ws_inter_sl, ws_vol_sl, ice_impeda_sl, hyd_cond_sl, diffus_sl) &
    !$ACC   CREATE(hyd_cond_li, diffus_li, drain_sl, evapo_soil_sl1, transpiration_sl) &
    !$ACC   CREATE(deficit_evapotrans, deficit_trans, deficit_sevap, remoist) &
    !$ACC   CREATE(wtr_rsdl_sl, weq_wsat_sl, weq_field_cap_sl, weq_p_wilt_sl, wtr_transp_residual) &
    !$ACC   CREATE(wtr_transp_corr, wtr_exinfcap, wtr_exwsat, wtr_soil_mobile_sl, last_soil_layer, hlp1, hlp2, hlp3, hlp4)

    !> Finish in case of negative soil moisture
    !
#ifndef _OPENACC
    DO is = 1, nsoil
      DO ic = 1, nc
        IF (wtr_soil_sl(ic,is) < 0._wp .AND. l_fract(ic)) THEN
          WRITE (message_text,*) 'negative soil moisture (routine start) at ', &
            & lat(ic),'N and ',lon(ic),'E, layer ', is, ':', NEW_LINE('a'), &
            & 'Soil moisture:    ', wtr_soil_sl(ic,is)
          CALL finish ('soilhyd', message_text)
        END IF
      END DO
    END DO
#endif

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO is = 1, nsoil
      DO ic = 1, nc

        !> Calculation of soil saturation, field capacity and wilting point taking into account soil ice
        !
        weq_wsat_sl     (ic,is) = vol_porosity_sl (ic,is) * dsoil(ic,is)
        weq_field_cap_sl(ic,is) = vol_field_cap_sl(ic,is) * dsoil(ic,is)
        weq_p_wilt_sl   (ic,is) = vol_p_wilt_sl(ic,is)    * dsoil(ic,is)

        ! supercooled water in the soil should not be able to move vertically
        ! @todo at some point p_wilt has to be renamed to something more appropriate, as it
        !       basically describes immobile adhesion water which is not really the same
        !       as the plant wilting point.
        wtr_rsdl_sl(ic,is) = MIN(MAX(wtr_soil_pot_scool_sl(ic,is), weq_p_wilt_sl(ic,is)), wtr_soil_sl(ic,is))

      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      !> Initialize surface runoff components
      wtr_exinfcap(ic) = infilt(ic)
      wtr_exwsat(ic)   = 0._wp
      !> Initialize active layers
      last_soil_layer(ic) = 0
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (dsoil(ic,is) > 0._wp) last_soil_layer(ic) = is
      END DO
    END DO
    !$ACC END PARALLEL

#ifndef _OPENACC
    IF (ANY(last_soil_layer(:) < 1)) CALL finish('soilhyd', 'Problem with no. of active soil layers (=0)')
#endif

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      ws_inter_sl(ic,1) = wtr_soil_sl(ic,1) + infilt(ic)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO is = 2, nsoil
      DO ic = 1, nc
        ws_inter_sl(ic,is) = wtr_soil_sl(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! Determining hydraulic conductivity and diffusivity on layers and at layer interfaces
    !   Note that the indexing uses the upper interface -> interface Ii = I(i,i-1)
    !   For numerical reasons, soil loops always cover all soil layers, even if they are not active (bedrock)
    CALL get_soilhyd_cond_properties(                                       &
          & soilhydmodel, interpol_mean,                                    &
          & nc, nsoil, dsoil(:,:),                                          &
          & ws_inter_sl(:,:), ice_soil_sl(:,:), wtr_soil_pot_scool_sl(:,:), &
          & weq_wsat_sl(:,:), weq_field_cap_sl(:,:), weq_p_wilt_sl(:,:),    &
          & hyd_cond_sat_sl(:,:), matrix_pot_sl(:,:),                       &
          & bclapp_sl(:,:), pore_size_index_sl(:,:),                        &
          & dt=dtime, last_soil_layer=last_soil_layer(:),                   &
          & ice_impedance=ice_impeda_sl(:,:),                               &
          & K=hyd_cond_sl(:,:),       D=diffus_sl(:,:),                     &
          & K_inter=hyd_cond_li(:,:), D_inter=diffus_li(:,:)                &
          )

    !< Diagnose evapotranspiration fluxes from the different soil layer based on the available
    !    mobile water
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO is = 1, nsoil
      DO ic = 1, nc
        wtr_soil_mobile_sl(ic,is) = wtr_soil_sl(ic,is) - wtr_rsdl_sl(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    CALL diagnose_evapotrans(nc, l_fract, lat, lon, nsoil, w_soil_wilt_fract, root_depth_sl, &
      &                      dsoil, wtr_field_cap_sl, transpiration_in_m, evapo_soil_in_m,   &
      &                      wtr_soil_mobile_sl, evapo_deficit, evapo_soil_sl1,              &
      &                      transpiration_sl, deficit_sevap, deficit_trans)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      ! Diagnose infiltration from the hydraulic conductivity of the top soil layer [m s-1 -->] [m]
      ! @todo this is implicitely done already within the ARNO scheme, thus it cannot be re-applied here.
      !       As soon as the Arno Scheme is retired and the different soil properties associated with
      !       orography are explizitly resolved, infiltration could be computed based on the conductivity instead.
      !       However, that with the current model configuration this approach either results in too high
      !       infiltration (if top layer conductivity is not constraint with the amout of water + infiltration
      !       of the top layer) or too low if the soil is not fully saturated. Smart ideas are needed here!
      ! infilt(ic)       = MIN(infilt(ic), hyd_cond_sl(ic,1) * dtime)

      ! Any water exceeding the infiltration capacity is given to the surface runoff (Horton runoff)
      !   Currently zero, as it is already included in the surface runoff from the ARNO scheme
      wtr_exinfcap(ic) = MAX(wtr_exinfcap(ic) - infilt(ic), 0._wp)

      ! Prescribe drainage from hydraulic conductivity of the lowest layer [m s-1] --> [m]
      drain_rock(ic) = hyd_cond_sl(ic,last_soil_layer(ic)) * dtime

      !< store actual soil water state, infiltration (upper boundary condition) and evaporation
      !  to compute the residual transport error
      wtr_transp_residual(ic) = SUM(wtr_soil_sl(ic,:)) + infilt(ic) + evapo_soil_sl1(ic) + SUM(transpiration_sl(ic,:))
    END DO
    !$ACC END PARALLEL LOOP

    ! Determination of Vertical movement of mobile water through the soil.
    !   Note that normalized storage values (volumetric soil moisture and storage change) are
    !   needed for the transport and immobile water is temporarily removed from the soil.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc
        IF (dsoil(ic,is) > 0._wp) THEN
          ws_vol_sl(ic,is) = (wtr_soil_sl(ic,is) - wtr_rsdl_sl(ic,is)) / dsoil(ic,is)
          transpiration_sl(ic,is) = transpiration_sl(ic,is) / dsoil(ic,is)
        END IF
        hlp1(ic,is) = 1._wp
        hlp4(ic,is) = transpiration_sl(ic,is) / dtime ! Normalized storage change term (Transpiration)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      hlp2(ic) = (infilt(ic) + evapo_soil_sl1(ic)) / dtime     ! Upper boundary condition
      hlp3(ic) = hyd_cond_sl(ic,last_soil_layer(ic))           ! Lower boundary condition
    END DO
    !$ACC END PARALLEL LOOP

    hlp5 = 1._wp
    ! CALL calc_vertical_transport(nc, n, d, alpha, dzf, top_bound, bot_bound, S, C, P, Q, X_old)
    !   *-bound, S are fluxes i.e. x/sec
    CALL calc_vertical_transport(           &
                        & dtime, hlp5,                    &
                        & last_soil_layer(1:nc), &
                        & dsoil(1:nc,1:nsoil),                     & ! Layer thickness
                        & hlp2(1:nc),                                          & ! Upper boundary condition
                        & hlp3(1:nc),                                          & ! Lower boundary condition
                        & hlp4(1:nc,1:nsoil),                         & ! Normalized storage change term (Transpiration)
                        & hlp1(1:nc,1:nsoil),                         &
                        & diffus_li(1:nc,1:nsoil),                 & ! Diffusivity
                        & hyd_cond_li(1:nc,1:nsoil),               & ! Hydraulic conductivity
                        & ws_vol_sl(1:nc,1:nsoil)                  & ! Normalized soil moisture
                        & )

    ! Converting back to absolute water content
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc
        IF (dsoil(ic,is) > 0._wp) THEN
          wtr_soil_sl(ic,is) = ws_vol_sl(ic,is) * dsoil(ic,is) + wtr_rsdl_sl(ic,is)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !< If lowest soil layer becomes negative, reduce bottom drainage
    !    In rare cases - when the bottom layer is very thin - also above layers might get negative.
    !    Generally, this only happens right after initialization.
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = nsoil, 1, -1
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (wtr_soil_sl(ic, last_soil_layer(ic)) < 0._wp .AND. l_fract(ic)) THEN
          IF (is <= last_soil_layer(ic) .AND. wtr_soil_sl(ic,is) < 0._wp) THEN
#ifndef _OPENACC
            WRITE (message_text_long,*) 'Soil moisture correction needed at ',   &
              &  lat(ic),'N and ',lon(ic),'E:' ,         NEW_LINE('a'),          &
              &  'Soil moisture:   ', wtr_soil_sl(ic,:), NEW_LINE('a'),          &
              &  'Bottom drainage: ', drain_rock(ic)
            CALL message('soilhyd', message_text_long, all_print=.TRUE.)
#endif
            drain_rock(ic) = MAX((drain_rock(ic) + wtr_soil_sl(ic,is)), 0._wp)
            wtr_soil_sl(ic,is) = 0._wp
#ifndef _OPENACC
            WRITE (message_text_long,*) 'After correction: wtr_soil_sl(ic,:): ', &
              &  lat(ic),'N and ',lon(ic),'E:' ,         NEW_LINE('a'),          &
              &  'Soil moisture:   ', wtr_soil_sl(ic,:), NEW_LINE('a'),          &
              &  'Bottom drainage: ', drain_rock(ic)
            CALL message('soilhyd', message_text_long, all_print=.TRUE.)
#endif
          END IF
        END IF
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

    !< Subtract new soil state and bottom drainage (lower boundary condition) from water balance
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      wtr_transp_residual(ic) = wtr_transp_residual(ic) - (SUM(wtr_soil_sl(ic,:)) + drain_rock(ic))
      wtr_soilhyd_res(ic) = wtr_transp_residual(ic)
    END DO

    !< Correct for the vertical transport error.

#ifndef _OPENACC
    !< 1. Complain if the overall transport error is too large
    DO ic = 1, nc
      IF (ABS(wtr_transp_residual(ic)) > 1.0e-10_wp .AND. l_fract(ic)) THEN
        WRITE (message_text_long,*) 'Soil water transport residual too large at ', &
          & lat(ic),'N and ',lon(ic),'E:',                   NEW_LINE('a'), &
          & 'Transport residual: ', wtr_transp_residual(ic), NEW_LINE('a'), &
          & 'Infiltration:       ', infilt(ic),              NEW_LINE('a'), &
          & 'Soil evaporation:   ', evapo_soil_sl1(ic),      NEW_LINE('a'), &
          & 'Bottom drainage:    ', drain_rock(ic),          NEW_LINE('a'), &
          & 'Transpiration:      ', transpiration_sl(ic,:),  NEW_LINE('a'), &
          & 'Soil moisture:      ', wtr_soil_sl(ic,:)
        CALL finish ('soilhyd', message_text_long)
      END IF
    END DO
#endif

    !< 2. Modify lower boundary condition (drainage) as much as possible
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      IF (ABS(wtr_transp_residual(ic)) > 0._wp) THEN
        wtr_transp_corr(ic)     = MIN(wtr_transp_residual(ic) + drain_rock(ic), 0._wp)
        drain_rock(ic)          = MAX(drain_rock(ic) + wtr_transp_residual(ic), 0._wp)
        wtr_transp_residual(ic) = wtr_transp_corr(ic)
      END IF
    END DO

    !$ACC END PARALLEL

    !< 3. Modify soil layers if necessary
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = nsoil, 1, -1
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (ABS(wtr_transp_residual(ic)) > 0._wp) THEN
          ! Make sure to correct the residual and not to put negative soil moisture onto the residual
          IF (wtr_soil_sl(ic,is) >= 0._wp .OR. wtr_soil_sl(ic,is) + wtr_transp_residual(ic) > 0._wp) THEN
            wtr_transp_corr(ic)     = MIN(wtr_transp_residual(ic) + wtr_soil_sl(ic,is), 0._wp)
            wtr_soil_sl(ic,is)       = MAX(wtr_soil_sl(ic,is) + wtr_transp_residual(ic), 0._wp)
            wtr_transp_residual(ic) = wtr_transp_corr(ic)
          END IF
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

#ifndef _OPENACC
    !> 4. Complain if correction is not sufficient
    DO ic = 1, nc
      IF (ABS(wtr_transp_residual(ic)) > 1.0e-10_wp .AND. l_fract(ic)) THEN
        WRITE (message_text_long,*) 'Cannot correct soil water transport residual at ', &
          & lat(ic),'N and ',lon(ic),'E:',                   NEW_LINE('a'), &
          & 'Transport residual: ', wtr_transp_residual(ic), NEW_LINE('a'), &
          & 'Infiltration:       ', infilt(ic),              NEW_LINE('a'), &
          & 'Soil evaporation:   ', evapo_soil_sl1(ic),      NEW_LINE('a'), &
          & 'Bottom drainage:    ', drain_rock(ic),          NEW_LINE('a'), &
          & 'Transpiration:      ', transpiration_sl(ic,:),  NEW_LINE('a'), &
          & 'Soil moisture:      ', wtr_soil_sl(ic,:)
        CALL finish ('soilhyd', message_text_long)
      END IF
    END DO
#endif

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

    ! Extract remaining bare soil evaporation from top layer (mobile and immoble soil water,
    ! because this has to be consistent with the computation of the upper soil layer relative humidity
    ! in sse_processes:relative_humidity_soil)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      IF (deficit_sevap(ic) < 0._wp) THEN                              ! remaining soil evaporation demand
        wtr_soil_sl(ic,1) = wtr_soil_sl(ic,1) + deficit_sevap(ic)      ! take water from uppermost soil layer
        IF (wtr_soil_sl(ic,1) < 0._wp) THEN                            ! if there is not enough liquid soil water
          deficit_sevap(ic) = wtr_soil_sl(ic,1)                        ! this is the amount of water missing
          wtr_soil_sl(ic,1) = 0._wp                                    ! soil moisture cannot be negative
        ELSE                                                           ! no deficit remains
          deficit_sevap(ic) = 0._wp
        END IF
      END IF
      evapo_deficit(ic) = evapo_deficit(ic) + deficit_sevap(ic)        ! remaining demand is added to deficit
    END DO

    ! Correct remaining evapotranspiration deficit
    ! @todo this is an unfortunate necessity as sometimes the evaporative demand is not fully
    !       backed by the water available in the land surface reservoirs.
    !       Changing the soil ice content here obviously violates the energy balance. A correction
    !       term is needed for this, if no better solution can be found to avoid the evaporation deficit.
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO ic = 1, nc
      deficit_evapotrans(ic) = deficit_sevap(ic) + deficit_trans(ic)
      weq_soil(ic) = SUM(wtr_soil_sl(ic,:)) + SUM(ice_soil_sl(ic,:))  ! water and ice within the column
#ifndef _OPENACC
      IF (weq_soil(ic) + deficit_evapotrans(ic) < -EPSILON(1.0_wp) .AND. l_fract(ic)) THEN
        WRITE (message_text,*) 'ET deficit cannot be compensated with the soil moisture storage'
        CALL message ('soilhyd', message_text, all_print=.TRUE.)
      END IF
#endif
    END DO

    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        IF (deficit_evapotrans(ic) < 0._wp .AND. weq_soil(ic) > 0._wp) THEN
          ! reduce soil moisture and ice relative to the water content of each layer
          wtr_soil_sl(ic,is) = MAX(0._wp, wtr_soil_sl(ic,is) + deficit_evapotrans(ic) * wtr_soil_sl(ic,is) / weq_soil(ic))
          ice_soil_sl(ic,is) = MAX(0._wp, ice_soil_sl(ic,is) + deficit_evapotrans(ic) * ice_soil_sl(ic,is) / weq_soil(ic))
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    ! Additional corrections for excess soil moisture
    ! In case of over-saturation water would not have infiltrated in the first place --> pile
    ! water upwards
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = nsoil, 2, -1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        ! If saturation moisture is exceeded, add the excess water to next higher soil layer
        IF (wtr_soil_sl(ic,is) > wtr_wsat_sl(ic,is)) THEN
          wtr_soil_sl(ic,is-1) = wtr_soil_sl(ic,is-1) + (wtr_soil_sl(ic,is) - wtr_wsat_sl(ic,is))
          wtr_soil_sl(ic,is)   = wtr_wsat_sl(ic,is)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      ! Top layer saturation excess is added to the surface runoff (Dunne runoff)
      IF (wtr_soil_sl(ic,1) > wtr_wsat_sl(ic,1)) THEN
        wtr_exwsat(ic)      = wtr_soil_sl(ic,1) - wtr_wsat_sl(ic,1)
        wtr_soil_sl(ic,1)   = wtr_wsat_sl(ic,1)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !
    ! Lateral Drainage (i.e. drainage from upper soil layers)
    ! ---- In supersaturated (or frozen) soils the water above the field capacity is assumed to
    !      move lateraly until it reaches the stream (or the neighbouring topo-tile).
    !      In this implementation all water above field capacity is added to the drainage. However,
    !      this is a very simplistic assumption. For e.g. high resolution simulations it might
    !      be worthwile to consider lateral drainage as a reservoir outflow with a predefined
    !      lagtime to allow for oversaturation of soils.
    ! @todo In the current model configuration lateral drainage seems to deteriorate river
    !      discharge while it has little effect on evapotranspiration. Therefore it is
    !      disabled, but this should be revised, e.g. for high resolution simulations.
    ! DO is = nsoil, 1, -1
    !   ! If maximum moisture is exceeded, add the excess water to lateral drainage
    !   IF (wtr_soil_sl(ic,is) > wtr_field_cap_sl(ic,is)) THEN
    !     drain_sl(ic)      = wtr_soil_sl(ic,is) - wtr_field_cap_sl(ic,is) ! Put excess water in lateral drainage
    !     wtr_soil_sl(ic,is) = wtr_soil_sl(ic,is) - drain_sl(ic)      ! remove lateral drainage of the layer from soil moisture
    !     IF (.NOT. ltpe_closed) THEN
    !       ! Add excess water to lateral drainage
    !       drain_lat(ic,is)   = drain_lat(ic,is)  + drain_sl(ic)
    !     ELSE
    !       ! no drainage for TPE closed case --> redirect overflow from drainage to overflow pool
    !       tpe_overflow(ic) = tpe_overflow(ic) + drain_sl(ic)
    !     END IF
    !   ELSE IF (ltpe_closed) THEN
    !     ! refill soil moisture from overflow pool
    !     remoist(ic)       = MIN(tpe_overflow(ic), wtr_field_cap_sl(ic,is) - wtr_soil_sl(ic,is))
    !     wtr_soil_sl(ic,is) = wtr_soil_sl(ic,is) + remoist(ic)
    !     tpe_overflow(ic)  = tpe_overflow(ic) - remoist(ic)
    !   END IF
    ! END DO

    ! Modification for TPE (open)
    IF (ltpe_open) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO is = 1, nsoil
        !$ACC LOOP GANG VECTOR
        DO ic = 1, nc
          !IF (i == nsoil - 2 .OR. i == nsoil - 1 .OR. i == nsoil) THEN   ! special adaptation of Shirisha
          IF (is == nsoil-1 .OR. is == nsoil) THEN
            IF (dsoil(ic,is) > 0._wp) THEN
              ! keep the two lowest soil layers always very wet
              wtr_soil_sl(ic,is) = MAX(wtr_soil_sl(ic,is),wtr_field_cap_sl(ic,is) * 0.9_wp)
            END IF
          END IF
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

    ! Add surface runoff components as infiltration overflow
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      infilt_over(ic) = wtr_exinfcap(ic) + wtr_exwsat(ic)
    END DO
    !$ACC END PARALLEL LOOP

#ifndef _OPENACC
    ! Check value range for soil moisture
    DO ic = 1, nc
      IF (ANY(wtr_soil_sl(ic,:) < -EPSILON(1.0_wp)) .AND. l_fract(ic)) THEN
        WRITE (message_text_long,*) 'negative soil moisture (routine end) at ', &
          & lat(ic),'N and ',lon(ic),'E:',                NEW_LINE('a'), &
          & 'Infiltration:     ', infilt(ic),             NEW_LINE('a'), &
          & 'Soil evaporation: ', evapo_soil_sl1(ic),     NEW_LINE('a'), &
          & 'Bottom drainage:  ', drain_rock(ic),         NEW_LINE('a'), &
          & 'Transpiration:    ', transpiration_sl(ic,:), NEW_LINE('a'), &
          & 'Soil moisture:    ', wtr_soil_sl(ic,:)
        CALL finish ('soilhyd', message_text_long)
      ELSE IF (ANY(wtr_soil_sl(ic,:) < 0.0_wp) .AND. l_fract(ic)) THEN
        ! WRITE (message_text,*) 'setting negative soil moisture (routine end) to zero at ', &
        !   & lat(ic),'N and ',lon(ic),'E:', NEW_LINE('a'), &
        !   & 'Soil moisture:    ', wtr_soil_sl(ic,:)
        ! Enable for debugging purpose, otherwise it fills up the logfile
        ! CALL message ('soilhyd', message_text)
        wtr_soil_sl(ic,:) = MAX(wtr_soil_sl(ic,:), 0.0_wp)
      END IF

      IF (ANY(wtr_soil_sl(ic,:) > wtr_wsat_sl(ic,:) + EPSILON(1.0_wp)) .AND. l_fract(ic)) THEN
        WRITE (message_text_long,*) 'soil moisture exceeds saturation capacity(routine end) at ', &
          &  lat(ic),'N and ',lon(ic),'E:',                 NEW_LINE('a'), &
          & 'Infiltration:       ', infilt(ic),             NEW_LINE('a'), &
          & 'Soil evaporation:   ', evapo_soil_sl1(ic),     NEW_LINE('a'), &
          & 'Bottom drainage:    ', drain_rock(ic),         NEW_LINE('a'), &
          & 'Transpiration:      ', transpiration_sl(ic,:), NEW_LINE('a'), &
          & 'Soil moisture:      ', wtr_soil_sl(ic,:),      NEW_LINE('a'), &
          & 'Saturation capcity: ', wtr_wsat_sl(ic,:)
        CALL finish ('soilhyd', message_text_long)
      END IF
    END DO
#else
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO is=1,nsoil
      DO ic=1,nc
        IF (wtr_soil_sl(ic,is) < 0._wp) THEN
          wtr_soil_sl(ic,is) = 0._wp
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP
#endif

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE soilhyd

  !----------------------------------------------------------------------------------------------

  SUBROUTINE diagnose_evapotrans( &
    & nc,                         &
    & l_fract,                    &
    & lat,                        &
    & lon,                        &
    & nsoil,                      &
    & w_soil_wilt_fract,          &
    & droot,                      &
    & dsoil,                      &
    & wtr_field_cap_sl,           &
    & transpiration_in_m,         &
    & evapo_soil_in_m,            &
    & wtr_soil_sl,                &
    & evapo_deficit,              &
    & evapo_soil_sl1,             &
    & transpiration_sl,           &
    & deficit_sevap,              &
    & deficit_trans               &
    & )

    !----------------------------------------------------------------------------------------------
    !< This routine distributes the evaporative demand of soil evaporation and transpiration
    !! between the soil layers of the root zone. Transpiration that cannot be satisfied
    !! from the root zone is subtracted from deeper layers. Also, only that part of soil evaporation
    !! that can be satisfied with mobile soil water is considered here.
    !! Note, that the routine only diagnoses the layered evaporation fluxes but the actual
    !! soil moisture update is done in the vertical transport routine. Also, it only consideres
    !! negative evapotranspiration fluxes (actual evapotranspiration) as positive fluxes (dew)
    !! are already considered in calc_surface_hydrology.
    !! The routine is based on routine soilchange by Stefan Hagemann.
    !----------------------------------------------------------------------------------------------

    ! arguments

    INTEGER,  INTENT(in)     :: nc                      ! vector length
    LOGICAL,  INTENT(IN)     :: l_fract(:)              ! Tile fraction > 0
    REAL(wp), INTENT(IN)     :: lat(:)                  ! Grid cell latitudes
    REAL(wp), INTENT(IN)     :: lon(:)                  ! Grid cell longitudes
    INTEGER,  INTENT(in)     :: nsoil                   ! number of soil layers
    REAL(wp), INTENT(in)     :: w_soil_wilt_fract       ! fractional plant wilting point [/]
    REAL(wp), INTENT(in)     :: droot(:,:)              ! root depth within layer [m]
    REAL(wp), INTENT(in)     :: dsoil(:,:)              ! soil depth until bedrock within layer [m]
    REAL(wp), INTENT(in)     :: wtr_field_cap_sl(:,:)   ! field capacity of the layer (reduced by ice)
    REAL(wp), INTENT(in)     :: transpiration_in_m(:)   ! transpiration [m]
    REAL(wp), INTENT(in)     :: evapo_soil_in_m(:)      ! bare soil evaporation [m]
    REAL(wp), INTENT(in)     :: wtr_soil_sl(:,:)        ! water content of the soil layer [m]
    REAL(wp), INTENT(inout)  :: evapo_deficit(:)        ! evaporation from unintented sources [m]
    REAL(wp), INTENT(OUT)    :: evapo_soil_sl1(:)       ! Soil evaporation from the top soil layer [m]
    REAL(wp), INTENT(OUT)    :: transpiration_sl(:,:)   ! Transpiration from soil layers [m]
    REAL(wp), INTENT(OUT)    :: deficit_trans(:)        ! Unaccounted transpiration [m]
    REAL(wp), INTENT(OUT)    :: deficit_sevap(:)        ! Unaccounted soil evaporation [m]

    ! local variables

    INTEGER  :: ic, is
    CHARACTER(len=4096) :: message_text_long ! long string for soil hydrology checks
    REAL(wp) :: dummy_wtr_soil(nc, nsoil)    ! dummy soil moisture storage [m]
    REAL(wp) :: deficit                      ! water deficit within the actual soil layer [m]
    REAL(wp) :: remaining                    ! water remaining in the soil layer [m]
    REAL(wp) :: rootfract                    ! fraction of the soil layer within the root zone [/]
    REAL(wp) :: fixed                        ! amount of water below the wilting point, not available for plants [m]
    REAL(wp) :: root_depth(nc)               ! Total depth of rootzone [m]

    !----------------------------------------------------------------------------------------------
    !  diagnose evaporation and transpiration fluxes
    !----------------------------------------------------------------------------------------------

    !$ACC DATA &
    !$ACC   CREATE(dummy_wtr_soil, root_depth)

#ifndef _OPENACC
    DO ic = 1, nc
      IF (ANY(wtr_soil_sl(ic,:) < 0._wp) .AND. l_fract(ic)) THEN
        WRITE (message_text,*) 'negative soil moisture (routine start) at ', &
          & lat(ic),'N and ',lon(ic),'E:', NEW_LINE('a'), &
          & 'Soil moisture: ', wtr_soil_sl(ic,:)
        CALL finish ('diagnose_evapotrans', message_text)
      END IF
    END DO
#endif

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      root_depth(ic) = SUM(droot(ic,:))
    END DO
    !$ACC END PARALLEL LOOP

    !  bare soil evaporation
    ! -----------------------

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO is = 1, nsoil
      DO ic = 1, nc
        ! Use dummy soil moisture storage for diagnostics, because the real soil moisture storage
        !   is only to be changed during the vertical transport
        dummy_wtr_soil(ic,is)   = wtr_soil_sl(ic,is)
        transpiration_sl(ic,is) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc

      evapo_soil_sl1(ic) = 0._wp
      deficit_sevap(ic)  = 0._wp     ! water missing for soil evaporation due to limited soil water
      deficit_trans(ic)  = 0._wp     ! water missing for transpiration due to limited soil water

      ! 1. Soil evaporation (mobile water only)
      ! ---------------

      IF (evapo_soil_in_m(ic) < 0._wp) THEN                             ! actual evaporation flux (not dew)
        dummy_wtr_soil(ic,1) = dummy_wtr_soil(ic,1) + evapo_soil_in_m(ic) ! take water from uppermost soil layer
        IF (dummy_wtr_soil(ic,1) < 0._wp) THEN                          ! if there is not enough mobile soil water
          deficit_sevap(ic) = dummy_wtr_soil(ic,1)                       ! this is the amount of water missing
          dummy_wtr_soil(ic,1) = 0._wp                                  ! soil moisture cannot be negative
        END IF
        evapo_soil_sl1(ic) = evapo_soil_in_m(ic) - deficit_sevap(ic)      ! reduce soil evaporation if top layer soil is too dry
      END IF

      !  2. Transpiration
      ! ---------------

      ! sealed grid cells (rooting_depth zero) and actual transpiration (which should not happen)
      !    => no water reachable for transpiration, add flux to deficit term
      ! @todo what does rooting depth=zero for bare soil tile?
      IF(root_depth(ic) <= 0._wp .AND. transpiration_in_m(ic) < 0._wp) THEN
        deficit_trans(ic) = transpiration_in_m(ic)                       ! add transpiration to deficit term
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR PRIVATE(deficit, rootfract, fixed, remaining)
      DO ic = 1, nc
        ! grid cell not sealed and actual transpiration (not dew) takes place
        !    => take water equally from all over the root zone
        IF (root_depth(ic) > 0._wp .AND. transpiration_in_m(ic) < 0._wp) THEN

          ! water needed for transpiration from the current layer
          transpiration_sl(ic,is) = transpiration_in_m(ic) * droot(ic,is) / root_depth(ic)

          IF (droot(ic,is) > 0._wp) THEN

            rootfract = droot(ic,is) / dsoil(ic,is) ! fraction of the soil layer that is within the root zone
            ! @todo at some point wtr_field_cap * w_soil_wilt_fract needs to be changed to the real wilting point moisture of the
            !       specific plant growing on this tile
            fixed = wtr_field_cap_sl(ic,is) * w_soil_wilt_fract * rootfract     ! amount of unavailable root zone water

            IF (dummy_wtr_soil(ic,is) * rootfract >= fixed) THEN                ! there is water in the root zone available
              remaining = dummy_wtr_soil(ic,is) * rootfract + transpiration_sl(ic,is) ! remaining water after transpiration
              IF (remaining < fixed) THEN                                       ! transpiration exceeds available water
                deficit = remaining - fixed
                dummy_wtr_soil(ic,is) = fixed + dummy_wtr_soil(ic,is) * (1._wp-rootfract)
              ELSE                                                                      ! enough soil water
                deficit = 0._wp
                dummy_wtr_soil(ic,is) = remaining + dummy_wtr_soil(ic,is) * (1._wp-rootfract)
              END IF
            ELSE                                                                        ! not enough water for transpiration
              deficit = transpiration_sl(ic,is)
            END IF
            transpiration_sl(ic,is) = transpiration_sl(ic,is) - deficit             ! reduce transpiration if necessary

            ! sum up deficits from the different layers
            deficit_trans(ic) = deficit_trans(ic) + deficit

          END IF
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      evapo_deficit(ic) = evapo_deficit(ic) + deficit_trans(ic)
    END DO
    !$ACC END PARALLEL LOOP

    ! If the transpiration deficit can be reduced by redistributing transpiration extraction from different layers,
    !   it should be done already here. The deficit should be taken from layers starting at the
    !   surface to accelerate the feedback on evapotranspiration.
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG VECTOR PRIVATE(remaining)
      DO ic = 1, nc
        IF (deficit_trans(ic) < 0._wp .AND. dummy_wtr_soil(ic,is) > wtr_field_cap_sl(ic,is) * w_soil_wilt_fract) THEN
          ! add water that can be still transpired to the transpiration to reduce the overall deficit
          remaining          = dummy_wtr_soil(ic,is) - &
            & MAX(dummy_wtr_soil(ic,is) + deficit_trans(ic), wtr_field_cap_sl(ic,is) * w_soil_wilt_fract)
          deficit_trans(ic)      = deficit_trans(ic)      + remaining
          dummy_wtr_soil(ic,is)   = dummy_wtr_soil(ic,is)   - remaining
          transpiration_sl(ic,is) = transpiration_sl(ic,is) - remaining
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

#ifndef _OPENACC
    DO ic = 1, nc
      ! make sure there would be no negative soil moisture after applying ET
      IF (ANY(dummy_wtr_soil(ic,:) < 0._wp) .AND. l_fract(ic)) THEN
        WRITE (message_text_long,*) 'negative soil moisture (routine end) at ', &
          & lat(ic),'N and ',lon(ic),'E:',                        NEW_LINE('a'), &
          & 'Soil evaporation:         ', evapo_soil_sl1(ic),     NEW_LINE('a'), &
          & 'Soil evaporation deficit: ', deficit_sevap(ic),      NEW_LINE('a'), &
          & 'Transpiration:            ', transpiration_sl(ic,:), NEW_LINE('a'), &
          & 'Transpiration deficit:    ', deficit_trans(ic),      NEW_LINE('a'), &
          & 'Soil moisture:            ', dummy_wtr_soil(ic,:)
        CALL finish ('diagnose_evapotrans', message_text_long)
      END IF
      ! check that diagnosed ET and expected ET is still identical
      IF ( ABS((evapo_soil_sl1(ic) + deficit_sevap(ic) + SUM(transpiration_sl(ic,:)) + deficit_trans(ic)) &
        &    - (evapo_soil_in_m(ic) + transpiration_in_m(ic)) ) > 1.0e-13_wp &
        &  .AND. l_fract(ic)) THEN
        WRITE (message_text_long,*) 'ET fluxes mismatch at ', &
          & lat(ic),'N and ',lon(ic),'E:',                                  NEW_LINE('a'), &
          & 'Initial mobile soil water:          ', wtr_soil_sl(ic,:),      NEW_LINE('a'), &
          & 'Expected soil evaporation:          ', evapo_soil_in_m(ic),    NEW_LINE('a'), &
          & 'Expected transpiration:             ', transpiration_in_m(ic), NEW_LINE('a'), &
          & 'Diagnosed soil evaporation :        ', evapo_soil_sl1(ic),     NEW_LINE('a'), &
          & 'Diagnosed transpiration:            ', transpiration_sl(ic,:), NEW_LINE('a'), &
          & 'Diagnosed soil moisture:            ', dummy_wtr_soil(ic,:),   NEW_LINE('a'), &
          & 'Remaining soil evaporation deficit: ', deficit_sevap(ic),      NEW_LINE('a'), &
          & 'Remaining transpiration deficit:    ', deficit_trans(ic)
        CALL finish ('diagnose_evapotrans', message_text_long)
      END IF
    END DO
#endif

  END SUBROUTINE diagnose_evapotrans

  !----------------------------------------------------------------------------------------------
  ! Calculates the snow fraction on lake ice.
  !
  ! @param [in]     weq_snow_lice    Snow depth on lake ice [m water equivalent]
  ! @param [out]    fract_snow_lice  Fraction of snow on lake ice
  !!
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE calc_wskin_fractions_lice( &
    & weq_snow_lice,                    & ! in
    & fract_snow_lice                   & ! out
    & )

    !$ACC ROUTINE SEQ

    REAL(wp),  INTENT(in)    :: weq_snow_lice
    REAL(wp),  INTENT(out)   :: fract_snow_lice

    fract_snow_lice = TANH(weq_snow_lice * 100._wp)

  END SUBROUTINE calc_wskin_fractions_lice

  SUBROUTINE calc_wet_fractions_veg(   &
    & dtime,                           & ! in
    & use_tmx,                         & ! in
    & skinres_max,                     & ! in
    & weq_pond_max,                    & ! in
    & fract_pond_max,                  & ! in
    & pond_dynamics_scheme,            & ! in
    & oro_stddev,                      & ! in
    & t_srf_old,                       & ! in
    & press_srf,                       & ! in
    & heat_tcoef,                      & ! in
    & q_air,                           & ! in
    & wtr_skin,                        & ! in
    & weq_pond,                        & ! in
    & weq_snow_soil,                   & ! in
    & weq_snow_can,                    & ! in
    & fract_snow_can,                  & ! inout
    & fract_skin,                      & ! out
    & fract_pond,                      & ! out
    & fract_wet,                       & ! out
    & fract_snow_soil                  & ! out
    & )

    USE mo_phy_schemes,            ONLY: qsat_water
    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_jsb_math_constants,     ONLY: pi
    USE mo_hydro_constants,        ONLY: Quad_, Tanh_, oro_crit

    REAL(wp), INTENT(in) :: &
      & dtime

    INTEGER, INTENT(in)  :: pond_dynamics_scheme

    LOGICAL, INTENT(in) :: use_tmx

    REAL(wp), INTENT(in), DIMENSION(:) :: &
      & oro_stddev,         &
      & t_srf_old,          &
      & press_srf,          &
      & heat_tcoef,         &
      & q_air,              &
      & skinres_max,        &
      & weq_pond_max,       &
      & fract_pond_max,     &
      & wtr_skin,           &
      & weq_pond,           &
      & weq_snow_soil,      &
      & weq_snow_can

    REAL(wp), INTENT(out), DIMENSION(:) :: &
      & fract_snow_can,      &
      & fract_skin,          &
      & fract_pond,          &
      & fract_wet,           &
      & fract_snow_soil

    REAL(wp) ::       &
      & qsat_srf_old, &
      & evapo_pot

    INTEGER :: nc, ic

    nc = SIZE(press_srf)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(qsat_srf_old, evapo_pot)
    DO ic=1,nc

      ! Snow fraction on canopy
      ! Wet skin reservoir fraction on soil and canopy (between 0 and 1)
      IF (skinres_max(ic) > EPSILON(1._wp)) THEN
        fract_snow_can(ic) = MIN(1._wp, weq_snow_can(ic) / skinres_max(ic))
        fract_skin(ic)      = MIN(1._wp, wtr_skin(ic)     / skinres_max(ic))
      ELSE
        fract_snow_can(ic) = 0._wp
        fract_skin(ic)      = 0._wp
      END IF

      ! Snow fraction on soil
      fract_snow_soil(ic) = Get_snow_fract_noforest(weq_snow_soil(ic), oro_stddev(ic)) ! snow on soil below forest
      fract_snow_soil(ic) = MERGE(fract_snow_can(ic), fract_snow_soil(ic), &
                               fract_snow_soil(ic) < EPSILON(1._wp) .AND. fract_snow_can(ic) > EPSILON(1._wp))

      ! Update pond fraction
      IF (weq_pond_max(ic) > EPSILON(1.0_wp) .AND. weq_pond(ic) > EPSILON(1.0_wp) .AND. &
        & pond_dynamics_scheme == Quad_) THEN
        fract_pond(ic) = MIN(1._wp, (weq_pond(ic) / weq_pond_max(ic))**0.5_wp) * fract_pond_max(ic)
      ELSE IF (weq_pond_max(ic) > EPSILON(1.0_wp) .AND. weq_pond(ic) > EPSILON(1.0_wp) .AND. &
        &      pond_dynamics_scheme == Tanh_) THEN
        fract_pond(ic) = MIN(1._wp, (TANH(pi * (weq_pond(ic) / weq_pond_max(ic))))**(oro_stddev(ic)/oro_crit)) &
          &              * fract_pond_max(ic)
      ELSE
        fract_pond(ic) = 0._wp
      END IF
      ! Potential evaporation using old values of air and surface humidity
      qsat_srf_old = qsat_water(t_srf_old(ic), press_srf(ic), use_convect_tables=.NOT. use_tmx)
      evapo_pot    = -1._wp * heat_tcoef(ic) * (q_air(ic) - qsat_srf_old)  ! Positive upwards

      ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than
      ! snow water content from soil and canopy; same for skin and pond reservoir
      IF (fract_snow_soil(ic) > 0._wp) THEN
        ! @todo Shouldn't one take rhoice here insteady of rhoh2o?
        fract_snow_soil(ic) = fract_snow_soil(ic) / MAX(1._wp, fract_snow_soil(ic) * evapo_pot * dtime             &
          &                                                    / (rhoh2o * (weq_snow_soil(ic) + weq_snow_can(ic))) &
          &                                            )
      END IF
      IF (fract_skin(ic) > 0._wp ) THEN
        fract_skin(ic) = fract_skin(ic) / MAX(1._wp, (1._wp - fract_snow_soil(ic)) * evapo_pot * dtime      &
          &                                          / (rhoh2o * MAX(EPSILON(1._wp), wtr_skin(ic)))  &
          &                                )
      END IF
      IF (fract_pond(ic) > 0._wp ) THEN
        fract_pond(ic) = fract_pond(ic) / MAX(1._wp, fract_pond(ic) * evapo_pot * dtime                &
        &                                          / (rhoh2o * MAX(EPSILON(1._wp), weq_pond(ic)))    &
        &                                  )
      END IF
      ! Combine the different wet surface fractions assuming that skin and pond locations
      ! would always coincide.
      fract_wet(ic) = MAX(fract_skin(ic), fract_pond(ic))

    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

  END SUBROUTINE calc_wet_fractions_veg

  SUBROUTINE calc_wet_fractions_bare( &
    & dtime,                            & ! in
    & use_tmx,                          & ! in
    & skinres_max,                      & ! in
    & weq_pond_max,                     & ! in
    & fract_pond_max,                   & ! in
    & pond_dynamics_scheme,             & ! in
    & oro_stddev,                       & ! in
    & t_srf_old,                        & ! in
    & press_srf,                        & ! in
    & heat_tcoef,                       & ! in
    & q_air,                            & ! in
    & wtr_skin,                         & ! in
    & weq_pond,                         & ! in
    & weq_snow_soil,                    & ! in
    & fract_skin,                       & ! out
    & fract_pond,                       & ! out
    & fract_wet,                        & ! out
    & fract_snow_soil                   & ! out
    & )

    USE mo_phy_schemes,            ONLY: qsat_water
    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_jsb_math_constants,     ONLY: pi
    USE mo_hydro_constants,        ONLY: Quad_, Tanh_, oro_crit

    REAL(wp), INTENT(in) :: &
      & dtime

    INTEGER, INTENT(in) :: pond_dynamics_scheme

    LOGICAL, INTENT(in) :: use_tmx

    REAL(wp), INTENT(in), DIMENSION(:) :: &
      & oro_stddev,         &
      & t_srf_old,          &
      & press_srf,          &
      & heat_tcoef,         &
      & q_air,              &
      & skinres_max,        &
      & weq_pond_max,       &
      & fract_pond_max,     &
      & wtr_skin,           &
      & weq_pond,           &
      & weq_snow_soil


    REAL(wp), INTENT(out), DIMENSION(:) :: &
      & fract_skin,          &
      & fract_pond,          &
      & fract_wet,           &
      & fract_snow_soil

    REAL(wp) ::       &
      & qsat_srf_old, &
      & evapo_pot

    INTEGER :: nc, ic

    nc = SIZE(press_srf)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(qsat_srf_old, evapo_pot)
    DO ic=1,nc

      IF (skinres_max(ic) > EPSILON(1._wp)) THEN
        fract_skin(ic) = MIN(1._wp, wtr_skin(ic) / skinres_max(ic))
      ELSE
        fract_skin(ic) = 0._wp
      END IF

      ! Snow fraction on soil
      fract_snow_soil(ic) = Get_snow_fract_noforest(weq_snow_soil(ic), oro_stddev(ic)) ! snow on soil below forest

      ! Update pond fraction
      IF (weq_pond_max(ic) > EPSILON(1.0_wp) .AND. weq_pond(ic) > EPSILON(1.0_wp) .AND. &
        & pond_dynamics_scheme == Quad_) THEN
        fract_pond(ic) = MIN(1._wp, (weq_pond(ic) / weq_pond_max(ic))**0.5_wp) * fract_pond_max(ic)
      ELSE IF (weq_pond_max(ic) > EPSILON(1.0_wp) .AND. weq_pond(ic) > EPSILON(1.0_wp) .AND. &
        &      pond_dynamics_scheme == Tanh_) THEN
        fract_pond(ic) = MIN(1._wp, (TANH(pi * (weq_pond(ic) / weq_pond_max(ic))))**(oro_stddev(ic)/oro_crit)) &
          &              * fract_pond_max(ic)
      ELSE
        fract_pond(ic) = 0._wp
      END IF

      ! Potential evaporation using old values of air and surface humidity
      qsat_srf_old = qsat_water(t_srf_old(ic), press_srf(ic), use_convect_tables=.NOT. use_tmx)
      evapo_pot    = heat_tcoef(ic) * (qsat_srf_old - q_air(ic)) ! Positive upwards

      ! Reduce snow cover on soil if snow loss during the time step due to potential evaporation is larger than
      ! snow water content from soil and canopy together
      IF (fract_snow_soil(ic) > 0._wp .AND. weq_snow_soil(ic) > EPSILON(1._wp)) THEN
        fract_snow_soil(ic) = fract_snow_soil(ic) / MAX(1._wp, fract_snow_soil(ic) * evapo_pot * dtime            &
          &                                                    / (rhoh2o * MAX(EPSILON(1._wp), weq_snow_soil(ic))) )
      END IF
      ! Same for water cover on soil
      IF (fract_skin(ic) > 0._wp .AND. wtr_skin(ic) > EPSILON(1._wp)) THEN
        fract_skin(ic) = fract_skin(ic) / MAX(1._wp, (1._wp - fract_snow_soil(ic)) * evapo_pot * dtime      &
          &                                                 / (rhoh2o * MAX(EPSILON(1._wp), wtr_skin(ic))) )
      END IF
      ! Same for surface water ponds
      IF (fract_pond(ic) > 0._wp ) THEN
        fract_pond(ic) = fract_pond(ic) / MAX(1._wp, fract_pond(ic) * evapo_pot * dtime                &
          &                                          / (rhoh2o * MAX(EPSILON(1._wp), weq_pond(ic)))    &
          &                                  )
      END IF

      ! Combine the different wet surface fractions assuming that skin and pond locations
      ! would always coincide.
      fract_wet(ic) = MAX(fract_skin(ic), fract_pond(ic))

    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT(1)

  END SUBROUTINE calc_wet_fractions_bare

  !< Roesch et al. 2002, Climate Dynamics
  !
#ifndef _OPENACC
  ELEMENTAL &
#endif
  REAL(wp) FUNCTION Get_snow_fract_noforest(snow, orodev)

    !$ACC ROUTINE SEQ

  USE mo_hydro_constants, ONLY: wsn2fract_const, wsn2fract_eps, wsn2fract_sigfac

    REAL(wp), INTENT(in) :: &
      & snow,   &
      & orodev

    Get_snow_fract_noforest = wsn2fract_const * TANH(snow * 100._wp) &
                               & * SQRT(snow * 1000._wp / (snow * 1000._wp + wsn2fract_eps + wsn2fract_sigfac * orodev))

  END FUNCTION Get_snow_fract_noforest

#ifndef _OPENACC
  ELEMENTAL &
#endif
  REAL(wp) FUNCTION get_canopy_cond_unstressed_simple(lai, par) RESULT(conductance)

    USE mo_hydro_constants, ONLY: k => conductance_k, a => conductance_a, b => conductance_b, c => conductance_c

    REAL(wp), INTENT(in) ::     &
                           lai, &
                           par

    ! Local variables
    !
!!$    REAL(wp) :: d(SIZE(par))
!!$    REAL(wp) :: zpar(SIZE(par))
    REAL(wp) :: d, zpar

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_canopy_cond_unstressed_simple'

    !$ACC ROUTINE SEQ

!!$    CALL message(TRIM(routine), 'Computing unstressed canopy conductance')

    zpar = MAX(1.E-10_wp, par)

!!$    WHERE (lai > EPSILON(1._wp))
    IF (lai > EPSILON(1._wp)) THEN
      d = (a + b*c) / (c * zpar)
      conductance = ( LOG((d * EXP(k*lai) + 1._wp) / (d + 1._wp)) * b / (d * zpar) - &
                    & LOG((d + EXP(-k*lai)) / (d + 1._wp))                           &
                    ) / (k * c)
    ELSE
      conductance = EPSILON(1._wp)
    END IF

  END FUNCTION get_canopy_cond_unstressed_simple

#ifndef _OPENACC
  ELEMENTAL &
#endif
  REAL(wp) FUNCTION get_canopy_cond_stressed_simple(cond_unstressed, water_stress, air_is_saturated) RESULT(conductance)

    REAL(wp), INTENT(in) :: &
      & cond_unstressed, &
      & water_stress
    LOGICAL,  INTENT(in) :: air_is_saturated

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_canopy_cond_stressed_simple'

!!$    CALL message(TRIM(routine), 'Computing unstressed canopy conductance')

    !$ACC ROUTINE SEQ

    IF (air_is_saturated) THEN
      conductance = EPSILON(1._wp)
    ELSE
      conductance = cond_unstressed * water_stress
    END IF

  END FUNCTION get_canopy_cond_stressed_simple

  !
  ! !>      Gets the water stress factor.
  !
  ! @param      w_soil               The w soil sl
  ! @param      w_soil_max           The w soil fc sl
  ! @param      w_soil_crit_fract    The w soil crit fract
  ! @param      w_soil_wilt_fract    The w soil pwp sl
  ! @param      water_stress_factor  The water stress factor
  !
  ! @return     The water stress factor. !
  !
#ifndef _OPENACC
  ELEMENTAL &
#endif
  FUNCTION get_water_stress_factor ( &
    & w_soil, w_soil_max, w_soil_crit_fract, w_soil_wilt_fract) RESULT(water_stress_factor)

    !$ACC ROUTINE SEQ

    ! TODO
    ! TBD: Maybe change later so that it uses geographically dependend critical values (like wilting point)

    REAL(wp), INTENT(in) :: w_soil            !< Soil water content
    REAL(wp), INTENT(in) :: w_soil_max        !< Soil water content at field capacity
    REAL(wp), INTENT(in) :: w_soil_crit_fract !< Fraction of max. soil water content at critical point
    REAL(wp), INTENT(in) :: w_soil_wilt_fract !< Fraction of max. soil water content at wilting point
    REAL(wp)             :: water_stress_factor

    REAL(wp) :: w_crit, w_wilt                       !< Soil water content at critical/wilting point

    w_crit = w_soil_max * w_soil_crit_fract
    w_wilt = w_soil_max * w_soil_wilt_fract

    IF (w_crit - w_wilt > 0._wp) THEN
      water_stress_factor = MAX (0._wp, &
        &                        MIN (1._wp, (w_soil - w_wilt) / (w_crit - w_wilt) ))
    ELSE
      water_stress_factor = 0._wp
    END IF

  END FUNCTION get_water_stress_factor

  SUBROUTINE calc_vertical_transport(dt, alpha, n_act, dzf, top_bound, bot_bound, S, C, P, Q, X)

    USE mo_util, ONLY: tdma_solver_vec ! OpenACC-enabled solver for tridiagonal matrix (Thomas algorithm)

  ! INPUT
  INTEGER, INTENT(in) :: n_act(:)
  REAL(wp), Intent(in):: dt, alpha,            &    ! time step
                       & top_bound(:), bot_bound(:), &    ! fluxes at upper and lower boundary
                       & dzf(:,:),               &    ! soil layer thickness
                       & C(:,:),                 &    ! heat capacity if used for heat transport,
                                                    ! set to 1 for water transport
                       & S(:,:),                 &    ! prescribed water/heat change flux per layer
                       & Q(:,:),                 &    ! Diffusivity at (upper) layer boundaries
                       & P(:,:)                       ! Conductivity at (upper) layer boundaries
  REAL(wp), INTENT(inout) :: &
                       & X(:,:)                       ! Normalized soil state for every layer

  ! LOCAL
  INTEGER            :: j, ic, nc
  INTEGER            :: nsoil
  REAL(wp) :: dzh(SIZE(S,1),SIZE(S,2)), &
    &         matrix(SIZE(S,1),SIZE(S,2),4) ! Matrix containing sub-diagonal (1), diagonal (2),
                                            ! super-diagonal (3) of triangular matrix and rhs (4) of
                                            ! linear equation system

  nc    = SIZE(S,1)
  nsoil = SIZE(S,2)

  !$ACC DATA &
  !$ACC   CREATE(dzh, matrix)

  !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
  DO j=1,nsoil
    DO ic=1,nc
      dzh(ic,j) = 0.0                             ! thickness of half levels
    END DO
  END DO
  !$ACC END PARALLEL LOOP

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP SEQ
  DO j=2,nsoil
    !$ACC LOOP GANG VECTOR
    DO ic=1,nc
      IF (j <= n_act(ic)) THEN
        dzh(ic,j) = (dzf(ic,j-1) + dzf(ic,j)) / 2._wp
      END IF
    END DO
  END DO
  !$ACC END PARALLEL

  !$ACC WAIT(1)

  ! Note that P & Q refer to the inter face above a given lvl not below.
  ! Hence for j == 1, there is no P(j) and Q(j)
  ! Make sure that this is consistent with the definition of input params.

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP SEQ
  DO j = 1, nsoil
    !$ACC LOOP GANG VECTOR
    DO ic=1,nc
      IF(j == 1) THEN
        ! Top soil layer
        matrix(ic,j,1) =  0.0_wp
        matrix(ic,j,2) =  P(ic,j+1) / dzh(ic,j+1) + dzf(ic,j) * C(ic,j) / (alpha * dt)
        matrix(ic,j,3) = -P(ic,j+1) / dzh(ic,j+1)
        matrix(ic,j,4) = -Q(ic,j+1) + dzf(ic,j) * S(ic,j) + dzf(ic,j) * C(ic,j) * X(ic,j) / (alpha * dt) + top_bound(ic)
      ELSE IF (j == n_act(ic)) THEN
        ! Bottom soil layers
        matrix(ic,j,1) = -P(ic,j) / dzh(ic,j)
        matrix(ic,j,2) =  P(ic,j) / dzh(ic,j) + dzf(ic,j) * C(ic,j) / (alpha * dt)
        matrix(ic,j,3) =  0.0_wp
        matrix(ic,j,4) =  Q(ic,j)          + dzf(ic,j) * S(ic,j) + dzf(ic,j) * C(ic,j) * X(ic,j) / (alpha * dt) - bot_bound(ic)
      ELSE IF (j < n_act(ic)) THEN
        ! Middle soil layers
        matrix(ic,j,1) = -P(ic,j) / dzh(ic,j)
        matrix(ic,j,2) =  P(ic,j) / dzh(ic,j) +  P(ic,j+1) / dzh(ic,j+1) + dzf(ic,j) * C(ic,j) / (alpha * dt)
        matrix(ic,j,3) = -P(ic,j+1) / dzh(ic,j+1)
        matrix(ic,j,4) =  Q(ic,j) - Q(ic,j+1) + dzf(ic,j) * S(ic,j) + dzf(ic,j) * C(ic,j) * X(ic,j) / (alpha * dt)
      ELSE
        matrix(ic,j,1) = 0._wp
        matrix(ic,j,2) = 1._wp
        matrix(ic,j,3) = 0._wp
        matrix(ic,j,4) = 1._wp
      END IF
    END DO
  END DO
  !$ACC END PARALLEL

#ifndef _OPENACC
  DO ic=1,nc
    IF( matrix(ic,1,2) .EQ. 0.0) THEN
      WRITE (message_text,*) 'something went terribly wrong ... rewrite equations!', C(ic,1), dzf(ic,1)
      CALL finish('calc_vertical_transport', message_text)
    END IF
  END DO
#endif

  CALL tdma_solver_vec(matrix(:,:,1), matrix(:,:,2), matrix(:,:,3), matrix(:,:,4), &
    &                  1, nsoil, 1, nc, X(:,:))

  !$ACC WAIT(1)
  !$ACC END DATA

  END SUBROUTINE calc_vertical_transport

  SUBROUTINE get_soilhyd_cond_properties(           &
            & soilhydmodel, interpol_mean,          &
            & nc, nsoil, dsoil,                     &
            & ws, ice, wscools,                     &
            & wsat, fc, pwp,                        &
            & k_sat, mpot_sat, bclapp, ps_index,    &
            & dt, last_soil_layer,                  &
            & ice_impedance, K, D, K_inter,         &
            & D_inter, mpot_act                     &
            )

    USE mo_hydro_constants,    ONLY: BrooksCorey_, ClappHornberger_, VanGenuchten_, &
      &                              Upstream_, Arithmetic_

    INTEGER,  INTENT(in) :: soilhydmodel         !< Model to determine conductivites & diffusivities (KD)
    INTEGER,  INTENT(in) :: interpol_mean        !< and how to derive values at layer interface (AVG scheme)
    INTEGER,  INTENT(in) :: nc                   !< vector length
    INTEGER,  INTENT(in) :: nsoil                !< number of below ground layers
    REAL(wp), INTENT(in) :: dsoil(:,:)           !< soil depth (until bedrock) per layer
    REAL(wp), INTENT(in) :: ws(:,:)              !< water content
    REAL(wp), INTENT(in) :: ice(:,:)             !< ice content
    REAL(wp), INTENT(in) :: wscools(:,:)         !< potentially supercooled water
    REAL(wp), INTENT(in) :: wsat(:,:)            !< soil depth (until bedrock) per layer
    REAL(wp), INTENT(in) :: fc(:,:)              !< field capacity [m]
    REAL(wp), INTENT(in) :: pwp(:,:)             !< residual water (until bedrock) per layer
    REAL(wp), INTENT(in) :: k_sat(:,:)           !< hydraaulic conductivity
    REAL(wp), INTENT(in) :: mpot_sat(:,:)        !< matric potential for saturated soils
    REAL(wp), INTENT(in) :: bclapp(:,:)          !< clapp & Hornberger exponent
    REAL(wp), INTENT(in) :: ps_index(:,:)        !< Pore size index
    REAL(wp), INTENT(in),  OPTIONAL :: dt                 !< time step length
    INTEGER,  INTENT(in),  OPTIONAL :: last_soil_layer(:) !< index of deepest soil layers (above bedrock)
    REAL(wp), INTENT(out), OPTIONAL :: ice_impedance(:,:) !< Impedance factor to account for ice blocking flowpaths
    REAL(wp), INTENT(out), OPTIONAL :: K(:,:)             !< hydraulic condctivity of the soil layer [m]
    REAL(wp), INTENT(out), OPTIONAL :: D(:,:)             !< diffusivity of the soil layer [m]
    REAL(wp), INTENT(out), OPTIONAL :: K_inter(:,:)       !< hydraulic condctivity at (upper) layer interface [m]
    REAL(wp), INTENT(out), OPTIONAL :: D_inter(:,:)       !< diffusivity at (upper) layer interface [m]
    REAL(wp), INTENT(out), OPTIONAL :: mpot_act(:,:)      !< matric potential at actual soil state

    ! local variables
    INTEGER  :: ic, is ! cell index and soil layer index
    REAL(wp) :: ck, local_dt
    REAL(wp) :: ws_vol(nc,nsoil), ice_vol(nc,nsoil), scool_vol(nc,nsoil),     &
             &  nvgn(nc,nsoil), mvgn(nc,nsoil), ws_rel_vol(nc,nsoil),         &
             &  ws_range_vol(nc,nsoil), wsat_vol(nc,nsoil), fc_vol(nc,nsoil), &
             &  pwp_vol(nc,nsoil), ws_free_vol(nc,nsoil), ice_imp(nc,nsoil)

    !$ACC DATA &
    !$ACC   CREATE(ws_vol, ice_vol, scool_vol, nvgn, mvgn, ws_rel_vol, ws_range_vol) &
    !$ACC   CREATE(wsat_vol, fc_vol, pwp_vol, ws_free_vol, ice_imp)

    ck = 8.0_wp ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc

        ! Set defaults
        ws_rel_vol(ic,is) = 0._wp
        ice_imp(ic,is)    = 1._wp
        nvgn(ic,is)       = ps_index(ic,is) + 1._wp
        mvgn(ic,is)       = ps_index(ic,is) / nvgn(ic,is)

        ! Convert to volumetric quantities
        IF (dsoil(ic,is) > 0._wp) THEN
          ws_vol(ic,is)    =  ws(ic,is)                      / dsoil(ic,is)
          ice_vol(ic,is)   =  ice(ic,is)                     / dsoil(ic,is)
          scool_vol(ic,is) =  MIN(wscools(ic,is), ws(ic,is)) / dsoil(ic,is)
          wsat_vol(ic,is)  =  wsat(ic,is)                    / dsoil(ic,is)
          fc_vol(ic,is)    =  fc(ic,is)                      / dsoil(ic,is)
          pwp_vol(ic,is)   =  pwp(ic,is)                     / dsoil(ic,is)
        ELSE
          ws_vol(ic,is)    =  0._wp
          ice_vol(ic,is)   =  0._wp
          scool_vol(ic,is) =  0._wp
          wsat_vol(ic,is)  =  0._wp
          fc_vol(ic,is)    =  0._wp
          pwp_vol(ic,is)   =  0._wp
        END IF
      END DO
    END DO
    !$ACC END LOOP

    IF (PRESENT(K)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          K(ic,is)       = 0._wp
        END DO
      END DO
      !$ACC END LOOP
    END IF
    IF (PRESENT(D)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          D(ic,is)       = 0._wp
        END DO
      END DO
      !$ACC END LOOP
    END IF
    IF (PRESENT(K_inter)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          K_inter(ic,is) = 0._wp
        END DO
      END DO
      !$ACC END LOOP
    END IF
    IF (PRESENT(D_inter)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          D_inter(ic,is) = 0._wp
        END DO
      END DO
      !$ACC END LOOP
    END IF
    IF (PRESENT(mpot_act)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          mpot_act(ic,is) = 0._wp
        END DO
      END DO
      !$ACC END LOOP
    END IF
    !$ACC END PARALLEL

    ! Compute impedance factor due to soil ice
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO is = 1, nsoil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = 1, nc
        IF (fc_vol(ic,is) + ice_vol(ic,is) > 1.0e-10_wp) THEN
          ice_imp(ic,is)  = 10._wp**(-6._wp * (ice_vol(ic,is) + scool_vol(ic,is) ) &
            & / (fc_vol(ic,is) + ice_vol(ic,is)))
        ELSE
          ice_imp(ic,is) = 1._wp
        END IF
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    IF (PRESENT(ice_impedance)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          ice_impedance(ic,is) = ice_imp(ic,is)
        END DO
      END DO
      !$ACC END LOOP
    END IF

    ! Compute relative soil moisture depending on soil hydrological model
    IF (soilhydmodel == ClappHornberger_) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          ws_free_vol(ic,is)   = MAX(0._wp, ws_vol(ic,is) - scool_vol(ic,is))
          ws_range_vol(ic,is)  = wsat_vol(ic,is)
        END DO
      END DO
      !$ACC END LOOP
    ELSE  ! Brooks & Corey; Van Genuchten
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          ws_free_vol(ic,is)   = MAX(0._wp, ws_vol(ic,is) - MAX(scool_vol(ic,is), pwp_vol(ic,is)))
          ws_range_vol(ic,is)  = fc_vol(ic,is) - pwp_vol(ic,is)
        END DO
      END DO
      !$ACC END LOOP
    END IF
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO is = 1, nsoil
      DO ic = 1, nc
        IF (ws_range_vol(ic,is) > 0._wp) THEN
          ws_rel_vol(ic,is)    = ws_free_vol(ic,is) / ws_range_vol(ic,is)
        END IF
      END DO
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    ! Making sure ws_rel is in the physical range
    ! This is especially important when infiltration is added to the first layer
    IF (soilhydmodel == VanGenuchten_) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          ws_rel_vol(ic,is) = MIN(0.9999_wp, MAX(0._wp, ws_rel_vol(ic,is))) ! doesn't work for wrel == 1;
                                                                            ! (1-ws_rel(ic,is)**(1/mvgn(ic,is)))**(-mvgn)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    ELSE  ! Brooks & Corey; Van Genuchten
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          ws_rel_vol(ic,is) = MIN(1._wp, MAX(0._wp, ws_rel_vol(ic,is)))
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ! Compute hydrological conductivity (K) for every layer
    IF (PRESENT(K)) THEN

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

      IF (soilhydmodel == BrooksCorey_ .OR. soilhydmodel == ClappHornberger_) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO is = 1, nsoil
          DO ic = 1, nc
            K(ic,is) =  k_sat(ic,is)                                                   &
              & * ws_rel_vol(ic,is)**(2._wp * bclapp(ic,is) + 3._wp)
          END DO
        END DO
        !$ACC END LOOP

      ELSE IF (soilhydmodel == VanGenuchten_) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO is = 1, nsoil
          DO ic = 1, nc
            IF (mvgn(ic,is) > 0._wp) THEN
              K(ic,is) = k_sat(ic,is) * ws_rel_vol(ic,is)**(0.5_wp) &
                & *(1._wp-(1._wp-ws_rel_vol(ic,is)**(1._wp/mvgn(ic,is)))**mvgn(ic,is))**2._wp
            ELSE
              K(ic,is) = 0._wp
            END IF
          END DO
        END DO
        !$ACC END LOOP
      END IF

      ! Reduce transport velocity in presence of ice
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is = 1, nsoil
        DO ic = 1, nc
          K(ic,is) = K(ic,is) * ice_imp(ic,is)
        END DO
      END DO
      !$ACC END LOOP

      !$ACC END PARALLEL
      !$ACC WAIT(1)

      ! Make sure percolation from top layer is limited to water content !!!
      IF (PRESENT(dt)) THEN
        local_dt = dt
        IF (local_dt > 0._wp) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
          DO ic = 1, nc
            K(ic,1) = MIN(K(ic,1), MAX(0._wp,(ws(ic,1)-pwp(ic,1)) / local_dt))
          END DO
          !$ACC END PARALLEL LOOP
        END IF
      END IF

    END IF

    ! Compute hydrological diffusivity (D) for every layer
    IF (PRESENT(D)) THEN

      IF (soilhydmodel == BrooksCorey_ .OR. soilhydmodel == ClappHornberger_) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO is = 1, nsoil
          DO ic = 1, nc
            IF (ws_range_vol(ic,is) > 0._wp) THEN
              D(ic,is) = -k_sat(ic,is) * mpot_sat(ic,is) * bclapp(ic,is) / ws_range_vol(ic,is) &
                & * ws_rel_vol(ic,is)**(1._wp * bclapp(ic,is) + 2._wp)
            ELSE
              D(ic,is) = 0._wp
            END IF
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      ELSE IF (soilhydmodel == VanGenuchten_) THEN
#ifndef __OPENACC__
        IF (.NOT. PRESENT(K)) THEN
          WRITE (message_text,*) 'VanGenuchten diffusivity computation requires hydrological conductivity'
          CALL finish ('get_soilhyd_cond_properties', message_text)
        END IF
#endif
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO is = 1, nsoil
          DO ic = 1, nc
            IF (ws_range_vol(ic,is) > 0._wp) THEN
              D(ic,is) = -K(ic,is)* mpot_sat(ic,is) * bclapp(ic,is) / ws_range_vol(ic,is) &
                & * ws_rel_vol(ic,is)**(1._wp/mvgn(ic,is)) * (1._wp-ws_rel_vol(ic,is)**(1._wp/mvgn(ic,is)))**(-mvgn(ic,is))
            ELSE
              D(ic,is) = 0._wp
            END IF
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      ! Reduce transport velocity in presence of ice
      ! However, ice seems to attract water --> is this the right approach for diffusivity?
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          D(ic,is) = D(ic,is) * ice_imp(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF
    !$ACC WAIT(1)

    ! Determining K & D on layer interface
    ! Possibility
    ! 1) Upstream mean [K_lev; Dmax(x_lev, x_lev-1)}
    ! 2) Adjusted arithmetic mean A (of k, D)
    ! Note that in the vertical transport scheme the notation refers to the upper interface of a level.
    ! Hence the loop goes from last_soil_layer to 2

    IF (PRESENT(K_inter)) THEN

#ifndef __OPENACC__
      IF (.NOT. PRESENT(K) .OR. .NOT. PRESENT(last_soil_layer)) THEN
        WRITE (message_text,*) 'K at layer interface cannot be computed without '//&
          & 'information about actual soil depth and computing K on layers'
        CALL finish ('get_soilhyd_cond_properties', message_text)
      END IF
#endif

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO is = nsoil, 2, -1
        !$ACC LOOP GANG VECTOR
        DO ic = 1, nc
          IF (is <= last_soil_layer(ic)) THEN
            IF (interpol_mean == Upstream_) THEN
              K_inter(ic,is)  =   K(ic,is-1)
            ELSE IF (interpol_mean == Arithmetic_) THEN
              IF (K(ic,is-1) > K(ic,is)) THEN ! MOVEMENT OF WETTING FRONTS LIMITED BY SATURATION OF SUBJACENT LAYER. ELSE ...
                K_inter(ic,is)  =   (K(ic,is)*dsoil(ic,is)+K(ic,is-1)*dsoil(ic,is-1))/(dsoil(ic,is)+dsoil(ic,is-1))
              ELSE ! ... PREVENT FLUX FROM COMPLETLY DRYING OUT UPPER LAYER
                K_inter(ic,is)  = K(ic,is-1)
              END IF
            END IF
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

    END IF

    IF (PRESENT(D_inter)) THEN

#ifndef __OPENACC__
      IF (.NOT. PRESENT(D) .OR. .NOT. PRESENT(last_soil_layer)) THEN
        WRITE (message_text,*) 'D at layer interface cannot be computed without '//&
          & 'information about actual soil depth and computing D on layers'
        CALL finish ('get_soilhyd_cond_properties', message_text)
      END IF
#endif
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO is = nsoil, 2, -1
        !$ACC LOOP GANG VECTOR
        DO ic = 1, nc
          IF (is <= last_soil_layer(ic)) THEN
            IF (interpol_mean == Upstream_) THEN
              D_inter(ic,is)  = MAX(D(ic,is),D(ic,is-1))
            ELSE IF (interpol_mean == Arithmetic_) THEN
              D_inter(ic,is)  = (D(ic,is)*dsoil(ic,is)+D(ic,is-1)*dsoil(ic,is-1))/(dsoil(ic,is)+dsoil(ic,is-1))
            END IF
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

    END IF

    ! Determine matric potential for actual soil moisture state
    ! Formula from Stuurop 2021 (https://doi.org/10.1016/j.coldregions.2021.103456),
    ! adjusted from Zhang 2007 (https://doi.org/10.1175/JHM605.1)
    IF (PRESENT(mpot_act)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          IF (ws_rel_vol(ic,is) > 0._wp) THEN
            IF (soilhydmodel == BrooksCorey_ .OR. soilhydmodel == ClappHornberger_) THEN
              mpot_act(ic,is) = mpot_sat(ic,is) * (ws_rel_vol(ic,is)**(-bclapp(ic,is)))
            ELSE
              mpot_act(ic,is) = mpot_sat(ic,is) * (ws_rel_vol(ic,is)**(-1._wp/mvgn(ic,is)) - 1._wp)**(1._wp/nvgn(ic,is))
            END IF
            mpot_act(ic,is) = MAX(-100._wp, mpot_act(ic,is) * ((1._wp + ck * ice_vol(ic,is))**2._wp))
          ELSE
            mpot_act(ic,is) = -100._wp
          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE get_soilhyd_cond_properties

#endif
END MODULE mo_hydro_process
