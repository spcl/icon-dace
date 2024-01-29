!> Contains the routines for the surface energy balance on LAKE lct_type.
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

MODULE mo_seb_lake
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_control,        ONLY: debug_on

  dsl4jsb_Use_processes SEB_, TURB_, RAD_, HYDRO_, A2L_
  dsl4jsb_Use_config(SEB_)

  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_energy_lake, update_surface_fluxes_lake

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_lake'

CONTAINS

  SUBROUTINE update_surface_energy_lake(tile, options)

    USE mo_phy_schemes,            ONLY: qsat_water, qsat_ice, surface_dry_static_energy
    USe mo_jsb_physical_constants, ONLY: cpd

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & t,                    &
      & s_star,               &
      & t_rad4,               &
      & t_eff4,               &
      & t_unfilt,             &
      & t_lwtr,               &
      & t_lice,               &
      & qsat_lwtr,            &
      & qsat_lice,            &
      & qsat_star,            &
      & s_lwtr,               &
      & s_lice,               &
      & heat_cap,             &
      & press_srf,            &
      & snow,                 &
      & latent_hflx_wtr,      &
      & sensible_hflx_wtr,    &
      & latent_hflx_ice,      &
      & sensible_hflx_ice,    &
      & hflx_form_ice,        &
      & hflx_cond_ice,        &
      & hflx_melt_lice,       &
      & hflx_melt_snow,       &
      & fract_lice,           &
      & depth_lice,           &
      & fract_snow_lice,      &
      & w_snow_lice,          &
      & w_snow,               &
      & rad_net_lwtr,         &
      & rad_net_lice

    ! Local variables
    REAL(wp) ::                     &
      & hcap_wtr(options%nc),       &
      & zero(options%nc),           &
      ! & lhflx_tmp(options%nc),      &
      & lhflx_tmp,      &
      & hflx_tmp(options%nc)
    REAL(wp), ALLOCATABLE :: &
      hcap_ice(:)

    LOGICAL ::                       &
      & l_ice(options%nc),           & ! Indicator for ice cover
      & l_ice_formation(options%nc), &
      & l_ice_melting(options%nc)
    LOGICAL :: config_l_ice_on_lakes, use_tmx

    INTEGER  :: iblk, ics, ice, nc, ic
    REAL(wp) :: dtime
    REAL(wp) ::         &
      & ml_depth,       & !< Depth of lake mixed layer
      & min_ice_depth,  & !< Minimum ice thickness to start ice formation
      & min_fract_lice, &
      & min_hcap_wtr

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_energy_lake'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_active(SEB_) .OR. .NOT. tile%contains_lake) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    use_tmx = model%config%use_tmx

    dsl4jsb_Get_config(SEB_)

    ml_depth      = dsl4jsb_Config(SEB_)%lake_mixed_layer_depth
    min_ice_depth = dsl4jsb_Config(SEB_)%lake_min_ice_depth
    config_l_ice_on_lakes = dsl4jsb_Config(SEB_)%l_ice_on_lakes

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(SEB_, t)              ! out
    IF (use_tmx) THEN
      dsl4jsb_Get_var2D_onChunk(SEB_, s_star)       ! out
      dsl4jsb_Get_var2D_onChunk(SEB_, t_rad4)       ! out
      dsl4jsb_Get_var2D_onChunk(SEB_, t_eff4)       ! out
      dsl4jsb_Get_var2D_onChunk(SEB_, qsat_star)    ! out
    END IF
    dsl4jsb_Get_var2D_onChunk(SEB_, t_unfilt)       ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lwtr)   ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, t_lwtr)         ! inout
    dsl4jsb_Get_var2D_onChunk(SEB_, qsat_lwtr)      ! out
    dsl4jsb_Get_var2D_onChunk(SEB_, s_lwtr)         ! out
    dsl4jsb_Get_var2D_onChunk(SEB_, heat_cap)       ! out

    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, w_snow)          ! out
    dsl4jsb_Get_var2D_onChunk(SEB_, latent_hflx_wtr)   ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, sensible_hflx_wtr) ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, fract_lice)        ! OUT

    IF (config_l_ice_on_lakes) THEN
      dsl4jsb_Get_var2D_onChunk(A2L_, snow)              ! in
      dsl4jsb_Get_var2D_onChunk(SEB_, t_lice)            ! inout
      dsl4jsb_Get_var2D_onChunk(SEB_, qsat_lice)         ! out
      dsl4jsb_Get_var2D_onChunk(SEB_, s_lice)            ! out
      dsl4jsb_Get_var2D_onChunk(SEB_, hflx_form_ice)     ! inout
      dsl4jsb_Get_var2D_onChunk(SEB_, hflx_cond_ice)
      dsl4jsb_Get_var2D_onChunk(SEB_, hflx_melt_lice)
      dsl4jsb_Get_var2D_onChunk(SEB_, hflx_melt_snow)

      dsl4jsb_Get_var2D_onChunk(SEB_, depth_lice)        ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_lice) ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, w_snow_lice)     ! inout

      dsl4jsb_Get_var2D_onChunk(SEB_, latent_hflx_ice)   ! in
      dsl4jsb_Get_var2D_onChunk(SEB_, sensible_hflx_ice) ! in

      dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lice)      ! in

      ALLOCATE(hcap_ice(nc))
      !$ACC ENTER DATA CREATE(hcap_ice)
    ELSE
      ALLOCATE(hflx_form_ice(nc))
      !$ACC ENTER DATA CREATE(hflx_form_ice)
    END IF

    !$ACC DATA CREATE(hcap_wtr, zero, l_ice, l_ice_formation, l_ice_melting, hflx_tmp)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      zero            (ic) = 0._wp
      hcap_wtr        (ic) = 0._wp
      l_ice_formation (ic) = .FALSE.
      l_ice_melting   (ic) = .FALSE.
      hflx_tmp(ic) = rad_net_lwtr(ic) + latent_hflx_wtr(ic) + sensible_hflx_wtr(ic) 
      IF (config_l_ice_on_lakes) THEN
        hcap_ice      (ic) = 0._wp
        l_ice         (ic) = fract_lice(ic) > 0.5_wp
      ELSE
        fract_lice    (ic) = 0._wp
        l_ice         (ic) = .FALSE.
        hflx_form_ice (ic) = 0._wp
      END IF
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    ! Calculate new water temperature and residual heat flux for ice formation over open water!
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      CALL calc_lake_wtr_temp( &
        & t_lwtr(ic),                                     & !< inout; Surface temperature over open water
        & hflx_tmp(ic),                                   & !< in;    Sum of all radiative and turbulent surface heat fluxes
        & hflx_form_ice(ic),                              & !< inout; Residual heat flux for ice formation (realized or unrealized)
        & hcap_wtr(ic),                                   & !< out;   Heat capacity of lake over open water
        & l_ice_formation(ic),                            & !< out;   Indicator for ice formation
        & l_ice_melting(ic),                              & !< inout; Indicator for ice melting = .FALSE.
        & l_ice(ic),                                      & !< in;    Indicator for ice cover
        & dtime,                                          & !< in;    Time step
        & ml_depth,                                       & !< in;    Mixed-layer depth
        & min_ice_depth)                                    !< in;    Minimum ice thickness
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    ! TODO How to put this on GPU? The following block is wrong!

      ! !$ACC KERNELS
    ! min_fract_lice = MINVAL(fract_lice)
    ! !$ACC END KERNELS

    ! !$ACC KERNELS
    ! min_hcap_wtr   = MINVAL(hcap_wtr)
    ! !$ACC END KERNELS

    ! IF (min_fract_lice <= 0.5_wp .AND. min_hcap_wtr < EPSILON(1._wp)) THEN
    !   CALL finish(TRIM(routine),'hcap_wtr can"t be zero')
    ! END IF

#ifndef _OPENACC
    IF (ANY(fract_lice(:) <= 0.5_wp .AND. hcap_wtr(:) < EPSILON(1._wp))) THEN
      CALL finish(TRIM(routine),'hcap_wtr can"t be zero')
    END IF
#endif

    IF (config_l_ice_on_lakes) THEN

      ! Calculate new ice thickness
      ! Tendency = freezing + conductive heat flux - ice melting + sublimation/deposition
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        CALL calc_lice_thickness(                           &
          & depth_lice(ic),                                 & !< inout; Ice thickness [m]
          & fract_lice(ic),                                 & !< inout; Indicator for ice cover
          & hflx_form_ice(ic),                              & !< inout; Heat flux for (complete) ice formation (l_ice_formation = T)
          & latent_hflx_ice(ic),                            & !< in;    Latent heat flux for ice
          & hflx_cond_ice(ic),                              & !< in;    Conductive heat flux for ice
          & hflx_melt_lice(ic),                             & !< inout; Heat flux from unrealized or complete ice melting
          & l_ice_formation(ic),                            & !< in;    Indicator for ice formation
          & l_ice_melting(ic),                              & !< inout; Indicator for ice melting
          & fract_snow_lice(ic),                            & !< in;    Fraction of snow on ice
          & min_ice_depth,                                  & !< in;    Mininum ice thickness
          & dtime                                           & !< in;    Time step
          & )
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
    
      ! Calculate new water temperature for the case of complete ice melting
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO ic=1,nc
          CALL calc_lake_wtr_temp( &
            & t_lwtr(ic),                                     & !< inout; Surface temperature over open water
            & zero(ic),                                       & !< in;    Sum of radiative and turbulent surface heat fluxes (not used)
            & hflx_melt_lice(ic),                             & !< inout; Residual heat flux from unrealized ice melting
            & hcap_wtr(ic),                                   & !< inout; Heat capacity of lake over open water (not changed)
            & l_ice_formation(ic),                            & !< inout; Indicator for ice formation (not changed)
            & l_ice_melting(ic),                              & !< in;    Indicator for ice melting = .TRUE.
            & l_ice(ic),                                      & !< in;    Indicator for ice cover
            & dtime,                                          & !< in;    Time step
            & ml_depth,                                       & !< in;    Mixed-layer depth
            & min_ice_depth)                                    !< in;    Minimum ice thickness
        END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      ! Calculate new ice temperature and heat fluxes from melting of snow and ice

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        lhflx_tmp = rad_net_lice(ic) + latent_hflx_ice(ic) + sensible_hflx_ice(ic)
        CALL calc_lake_ice_temp(                            &
          & t_lice(ic),                                     & ! inout
          & hcap_ice(ic),                                   & ! out
          & hflx_cond_ice(ic),                              & ! out
          & hflx_melt_snow(ic),                             & ! out
          & hflx_melt_lice(ic),                             & ! inout
          & lhflx_tmp,                                      & !< in;    Sum of all radiative and turbulent surface heat fluxes
          & latent_hflx_ice(ic),                            & !< in;
          & snow(ic),                                       & !< in; Snow fall
          & depth_lice(ic),                                 &
          & w_snow_lice(ic),                                &
          & fract_snow_lice(ic),                            &
          & min_ice_depth,                                  &
          & dtime)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

    END IF

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
    DO ic=1,nc
      qsat_lwtr(ic)  = qsat_water(t_lwtr(ic), press_srf(ic), use_convect_tables=.NOT. use_tmx)
      s_lwtr(ic)     = surface_dry_static_energy(t_lwtr(ic), qsat_lwtr(ic), cpd)
    END DO
    !$ACC END PARALLEL LOOP

    IF (config_l_ice_on_lakes) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO ic=1,nc
        qsat_lice(ic) = qsat_ice(t_lice(ic), press_srf(ic), use_convect_tables=.NOT. use_tmx)
        t(ic)         = (1._wp - fract_lice(ic)) * t_lwtr(ic) + fract_lice(ic) * t_lice(ic)
        t_unfilt(ic)  = t(ic)
        s_lice(ic)    = surface_dry_static_energy(t_lice(ic), qsat_lice(ic), cpd)
        IF (use_tmx) THEN
          s_star(ic) = (1._wp - fract_lice(ic)) * s_lwtr(ic)    + fract_lice(ic) * s_lice(ic)
          t_rad4(ic) = (1._wp - fract_lice(ic)) * t_lwtr(ic)**4 + fract_lice(ic) * t_lice(ic)**4
          t_eff4(ic) = t_rad4(ic)
          qsat_star(ic) = (1._wp - fract_lice(ic)) * qsat_lwtr(ic) + fract_lice(ic) * qsat_lice(ic)
        END IF
        heat_cap(ic)  = (1._wp - fract_lice(ic)) * hcap_wtr(ic) + fract_lice(ic) * hcap_ice(ic)
        w_snow(ic)    = fract_lice(ic) * w_snow_lice(ic)
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC EXIT DATA DELETE(hcap_ice)
      DEALLOCATE(hcap_ice)
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO ic=1,nc
        t(ic) = t_lwtr(ic)
        t_unfilt(ic) = t(ic)
        IF (use_tmx) THEN
          s_star(ic) = s_lwtr(ic)
          t_rad4(ic) = t_lwtr(ic)**4
          t_eff4(ic) = t_rad4(ic)
          qsat_star(ic) = qsat_lwtr(ic)
        END IF
        heat_cap(ic) = hcap_wtr(ic)
        w_snow(ic) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC EXIT DATA DELETE(hflx_form_ice)
      DEALLOCATE(hflx_form_ice)
    END IF

    !$ACC END DATA

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_energy_lake

  ! Calculate new surface water temperature
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_lake_wtr_temp(t_lwtr, hflx, hflx_resid, hcap, l_ice_formation, l_ice_melting, &
    &                                          l_ice, dtime, ml_depth, min_ice_depth)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, clw, rhoi, alf

    REAL(wp), INTENT(inout) :: t_lwtr            !< Surface temperature over open water
    REAL(wp), INTENT(in)    :: hflx             !< Sum of all radiative and turbulent surface heat fluxes
    REAL(wp), INTENT(inout) :: hflx_resid       !< Heat flux for ice formation or melting
                                                !! l_ice_melting = .TRUE.  : heat flux from complete ice melting (input)
                                                !!                           zero on output
                                                !!               = .FALSE. : heat flux from unrealized ice formation from
                                                !!                           previous time step (input),
                                                !!                           heat flux from realized or unrealized ice formation
                                                !!                           in this time step (output)
    REAL(wp), INTENT(inout) :: hcap             !< Heat capacity of lake over open water
    LOGICAL,  INTENT(inout) :: l_ice_formation  !< Indicator for ice formation
    LOGICAL,  INTENT(inout) :: l_ice_melting    !< Indicator for (complete) ice melting
    LOGICAL,  INTENT(in)    :: l_ice            !< Indicator for ice cover

    REAL(wp), INTENT(in)    :: dtime            !< Time step
    REAL(wp), INTENT(in)    :: ml_depth         !< Mixed-layer depth
    REAL(wp), INTENT(in)    :: min_ice_depth    !< Mininum ice thickness

    REAL(wp) :: &
      & t_star, &
      & t_freeze

    ! Lake tile is either not covered by ice or completely covered.
    IF (l_ice) RETURN

    IF (l_ice_melting) THEN
      t_lwtr = t_lwtr + dtime * hflx_resid / hcap
      hflx_resid = 0._wp
      l_ice_melting = .FALSE.
    ELSE
      ! Heat capacity of lake over open water
      ! Note: This is currently constant but in the future, different mixed-layer depths might be used
      hcap = rhoh2o * clw * ml_depth

      ! Water temperature below which a new slab of ice is formed
      t_freeze = tmelt - min_ice_depth * rhoi * alf / hcap

      ! Preliminary water temperature
      t_star = t_lwtr + dtime * (hflx + hflx_resid) / hcap

      IF (t_star >= tmelt) THEN              ! No ice formation
        t_lwtr = t_star
        hflx_resid = 0._wp
        l_ice_formation = .FALSE.
      ELSE                                   ! Realized or unrealized ice formation
        t_lwtr = tmelt
        hflx_resid = (t_star - tmelt) * hcap / dtime      ! < 0
        l_ice_formation = t_star <= t_freeze
        ! l_ice_formation = .TRUE.  : hflx_resid will be used to form a new slab of ice
        !                 = .FALSE. : cooling is not large enough to form a new slab of ice and
        !                             hflx_resid will be used as forcing flux in the next time
        !                             step to compute water temperature t_lwtr
      END IF

    END IF

  END SUBROUTINE calc_lake_wtr_temp

  ! Calculate new ice temperature, snow depth and heat fluxes
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_lake_ice_temp(t_lice, hcap_ice, hflx_cond_ice, hflx_melt_snow, hflx_melt_lice, &
    &                                          hflx, lhflx, snowfall, hi, hs, fs, hi_min, dtime)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, rhoi, rhoi, rhos, als, alf, ki, ks, ci, cs !, clw

    REAL(wp), INTENT(inout) :: t_lice
    REAL(wp), INTENT(out)   :: hcap_ice       !< Heat capacity of thin ice surface layer
    REAL(wp), INTENT(out)   :: hflx_cond_ice
    REAL(wp), INTENT(out)   :: hflx_melt_snow
    REAL(wp), INTENT(inout) :: hflx_melt_lice  !< in:  Residual heat flux from unrealized ice melting
                                              !! out: Heat flux from melting of ice
    REAL(wp), INTENT(in)    :: hflx           !< Sum of all surface heat fluxes (net radiation + sensible + latent)
    REAL(wp), INTENT(in)    :: lhflx          !< Latent surface heat flux for ice
    REAL(wp), INTENT(in)    :: snowfall       !< Snow fall                [kg m-2 s-1]
    REAL(wp), INTENT(in)    :: hi             !< Ice thickness [m]
    REAL(wp), INTENT(inout) :: hs             !< Snow depth [m water equivalent]
    REAL(wp), INTENT(in)    :: fs             !< Fraction of snow on ice
    REAL(wp), INTENT(in)    :: hi_min
    REAL(wp), INTENT(in)    :: dtime

    REAL(wp) :: hi_eff     !< Effective ice thickness including snow cover
    REAL(wp) :: hflx_melt  !< Excess heat available for melting snow and ice
    REAL(wp) :: delta_hs
    REAL(wp) :: delta_t_lice

    hcap_ice = rhoi * ci * hi_min + rhoh2o * cs * hs

    hi_eff = hi + (ki / ks) * (rhoh2o / rhos) * hs

    IF (hi < hi_min) THEN
      t_lice = tmelt
      hflx_cond_ice = 0._wp
      hflx_melt_snow = 0._wp
      hflx_melt_lice = 0._wp
      hs = 0._wp
      RETURN
    END IF

    ! Preliminary snow depth
    hs = hs + (dtime / rhoh2o) * ( snowfall           & ! Snowfall
      &                          + fs * lhflx / als   & ! Sublimation (< 0) resp. deposition (> 0)
      &                          )
    hs = MAX(hs, 0._wp)

    ! Preliminary ice temperature
    t_lice = ( hflx + hflx_melt_lice + hcap_ice * t_lice / dtime + ki * tmelt / hi_eff ) / ( hcap_ice / dtime + ki / hi_eff)

    hflx_melt_snow = 0._wp
    hflx_melt_lice  = 0._wp

    IF (t_lice > tmelt .AND. hs > 0._wp) THEN        ! Melting for snow first
      ! Excess heat which is available for melting of snow (first) and ice (as soon as snow is  melted away completely)
      hflx_melt = ( hcap_ice / dtime + ki /hi_eff ) * (t_lice - tmelt)
      ! Change in snow depth due to melting
      delta_hs = MIN(hs, dtime * hflx_melt / alf / rhoh2o)             ! positive
      ! Change in ice temperature resulting from snow melt
      hflx_melt_snow = delta_hs * rhoh2o *alf / dtime                  ! Heat that was used for melting the snow
      delta_t_lice = hflx_melt_snow / ( hcap_ice/dtime + ki/hi_eff)     ! positive
      t_lice = t_lice - delta_t_lice
      ! Update snow depth, eff. ice thickness and heat capacity of thin ice layer (including snow)
      hs = hs - delta_hs
      hi_eff = hi + (ki / ks) * (rhoh2o / rhos) * hs
      hcap_ice = rhoi * ci * hi_min + rhoh2o * cs * hs
    END IF
    IF (t_lice > tmelt) THEN                         ! Melting of ice after complete melting of snow
      hflx_melt_lice = ( hcap_ice/dtime + ki/hi_eff) * (t_lice - tmelt)
      t_lice = tmelt
    END IF

    ! Diagnose conductive heat flux for changing ice thickness in next time step
    hflx_cond_ice = - ki * (t_lice - tmelt) / hi_eff

  END SUBROUTINE calc_lake_ice_temp

  ! Calculation of ice depth and fraction
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_lice_thickness( &
    & hi,                                        &
    & fi,                                        &
    & hflx_form_ice,                             &
    & lhflx,                                     &
    & hflx_cond_ice,                             &
    & hflx_melt_lice,                            &
    & l_ice_formation,                           &
    & l_ice_melting,                             &
    & fs,                                        &
    & hi_min,                                    &
    & dtime                                      &
    & )

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: rhoi, alf, als

    REAL(wp), INTENT(inout) :: hi              !< Ice thickness [m]
    REAL(wp), INTENT(inout) :: fi              !< Ice fraction  []
    REAL(wp), INTENT(inout) :: hflx_form_ice   !< Residual heat flux from ice formation (complete or unrealized) [W m-2 = J m-2 s-1]
    REAL(wp), INTENT(in)    :: lhflx           !< Latent heat flux of ice  [W m-2 = J m-2 s-1]
    REAL(wp), INTENT(in)    :: hflx_cond_ice   !< Conductive heat flux of ice [W m-2 = J m-2 s-1]
    REAL(wp), INTENT(inout) :: hflx_melt_lice   !< Heat flux from melting of ice [W m-2 = J m-2 s-1]
                                               !! in:  flux from unrealized ice melting from previous time step
                                               !! out: flux from unrealized or complete ice melting during this time step
    LOGICAL,  INTENT(in)    :: l_ice_formation !< Indicator for ice formation
    LOGICAL,  INTENT(inout) :: l_ice_melting   !< Indicator for complete ice melting
    REAL(wp), INTENT(in)    :: fs              !< Fraction of snow on ice

    REAL(wp), INTENT(in)    :: hi_min          !< Minimum ice thickness [m]
    REAL(wp), INTENT(in)    :: dtime

    IF (l_ice_formation) THEN                              ! Formation of a new slab of ice

      hi = - hflx_form_ice * dtime / (rhoi * alf)
      fi = 1._wp
      hflx_form_ice = 0._wp

    ELSE IF (hi >= hi_min) THEN                            ! Ice already present

      hi = hi + (dtime / rhoi) * (   (hflx_cond_ice - hflx_melt_lice) / alf    &
        &                          + (1._wp - fs) * lhflx / als               &
        &                        )

      IF (hi >= hi_min) THEN
        hflx_melt_lice = 0._wp
        fi = 1._wp
        l_ice_melting = .FALSE.
      ELSE IF (hi <= 0._wp) THEN                              ! Complete melting
        hflx_melt_lice = - rhoi * alf * hi / dtime             ! >= 0
        hi = 0._wp
        fi = 0._wp
        l_ice_melting = .TRUE.
      ELSE                                                    ! Unrealized melting
        hflx_melt_lice = rhoi * alf * (hi_min - hi) / dtime    ! > 0
        hi = hi_min
        fi = 1._wp
        l_ice_melting = .FALSE.
      END IF

    END IF

  END SUBROUTINE calc_lice_thickness

  SUBROUTINE update_surface_fluxes_lake(tile, options)

    USE mo_phy_schemes,            ONLY: heat_transfer_coef ! q_effective, surface_dry_static_energy
    USE mo_jsb_physical_constants, ONLY: alv, als, cpd, cvd

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(A2L_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: sensible_hflx
    dsl4jsb_Real2D_onChunk :: latent_hflx
    dsl4jsb_Real2D_onChunk :: s_lwtr
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: t_acoef_wtr
    dsl4jsb_Real2D_onChunk :: t_bcoef_wtr
    dsl4jsb_Real2D_onChunk :: drag_wtr
    dsl4jsb_Real2D_onChunk :: latent_hflx_wtr
    dsl4jsb_Real2D_onChunk :: sensible_hflx_wtr
    dsl4jsb_Real2D_onChunk :: evapo_wtr
    dsl4jsb_Real2D_onChunk :: s_lice
    dsl4jsb_Real2D_onChunk :: t_acoef_ice
    dsl4jsb_Real2D_onChunk :: t_bcoef_ice
    dsl4jsb_Real2D_onChunk :: drag_ice
    dsl4jsb_Real2D_onChunk :: latent_hflx_ice
    dsl4jsb_Real2D_onChunk :: sensible_hflx_ice
    dsl4jsb_Real2D_onChunk :: evapo_ice
    dsl4jsb_Real2D_onChunk :: fract_lice
    dsl4jsb_Real2D_onChunk :: ch                 ! for new turb

    REAL(wp) ::                  &
      & s_air,     &  !< Dry static energy at lowest atmospheric level
      & heat_tcoef    !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc, ic
    REAL(wp) :: steplen
    LOGICAL  :: use_tmx

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_fluxes_lake'

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (.NOT. tile%Is_process_active(SEB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get model
    model => Get_model(tile%owner_model_id)
    use_tmx = model%config%use_tmx

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)

    IF (use_tmx) THEN
      dsl4jsb_Get_memory(TURB_)
      dsl4jsb_Get_var2D_onChunk(A2L_,   t_acoef)       ! IN
      dsl4jsb_Get_var2D_onChunk(A2L_,   t_bcoef)       ! IN
      dsl4jsb_Get_var2D_onChunk(TURB_, ch)               ! in
    ELSE
      dsl4jsb_Get_var2D_onChunk(A2L_,   t_acoef_wtr)       ! IN
      dsl4jsb_Get_var2D_onChunk(A2L_,   t_bcoef_wtr)       ! IN
      dsl4jsb_Get_var2D_onChunk(A2L_,   drag_wtr)        ! in
    END IF

    dsl4jsb_Get_var2D_onChunk(SEB_,   s_lwtr)            ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,   sensible_hflx)     ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   sensible_hflx_wtr) ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   latent_hflx)       ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   latent_hflx_wtr)   ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_wtr)         ! in

    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      IF (.NOT. use_tmx) THEN
        dsl4jsb_Get_var2D_onChunk(A2L_,   t_acoef_ice)       ! IN
        dsl4jsb_Get_var2D_onChunk(A2L_,   t_bcoef_ice)       ! IN
        dsl4jsb_Get_var2D_onChunk(A2L_,   drag_ice)      ! in
      END IF
      dsl4jsb_Get_var2D_onChunk(SEB_,   s_lice)            ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,   sensible_hflx_ice) ! out
      dsl4jsb_Get_var2D_onChunk(SEB_,   latent_hflx_ice)   ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_ice)         ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,   fract_lice)        ! in
    END IF

    ! ================================================================================================================================
    ! Surface fluxes for lake water

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(s_air, heat_tcoef)
    DO ic=1,nc

      IF (use_tmx) THEN
        s_air = t_bcoef(ic)  ! Old dry static energy at lowest atmospheric level
        heat_tcoef = ch(ic)  ! Transfer coefficient; TODO: distinguish between wtr and ice?
        sensible_hflx_wtr(ic) = heat_tcoef * (s_air - s_lwtr(ic)) * cvd / cpd ! Sensible energy flux
      ELSE
        s_air = t_acoef_wtr(ic) * s_lwtr(ic) + t_bcoef_wtr(ic)    ! New dry static energy at lowest atmospheric level by backsubstitution
        heat_tcoef = heat_transfer_coef(drag_wtr(ic), steplen)    ! Transfer coefficient
        sensible_hflx_wtr(ic) = heat_tcoef * (s_air - s_lwtr(ic)) ! Sensible heat flux
      END IF
      latent_hflx_wtr(ic) = alv * evapo_wtr(ic) ! Latent heat flux

    END DO
    !$ACC END PARALLEL

    ! ================================================================================================================================
    ! Surface fluxes for lake ice

    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG(STATIC: 1) VECTOR PRIVATE(s_air, heat_tcoef)
      DO ic=1,nc

        IF (use_tmx) THEN
          s_air = t_bcoef(ic) ! Old dry static energy at lowest atmospheric level
          heat_tcoef = ch(ic) ! Transfer coefficient; TODO: distinguish between wtr and ice?
          sensible_hflx_ice(ic) = heat_tcoef * (s_air - s_lice(ic)) * cvd / cpd ! Sensible energy flux
        ELSE
          s_air = t_acoef_ice(ic) * s_lice(ic) + t_bcoef_ice(ic)        ! New dry static energy at lowest atmospheric level by backsubstitution
          heat_tcoef = heat_transfer_coef(drag_ice(ic), steplen)        ! Transfer coefficient
          sensible_hflx_ice(ic) = heat_tcoef * (s_air - s_lice(ic)) ! Sensible heat flux
        END IF
        latent_hflx_ice(ic) = als * evapo_ice(ic) ! Latent heat flux

        ! ================================================================================================================================
        ! Surface fluxes for lake (wtr and ice aggregated)
        sensible_hflx(ic) = (1._wp - fract_lice(ic)) * sensible_hflx_wtr(ic) + fract_lice(ic) * sensible_hflx_ice(ic)
        latent_hflx  (ic) = (1._wp - fract_lice(ic)) * latent_hflx_wtr  (ic) + fract_lice(ic) * latent_hflx_ice  (ic)

      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO ic=1,nc
        sensible_hflx(ic) = sensible_hflx_wtr(ic)
        latent_hflx  (ic) = latent_hflx_wtr  (ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_fluxes_lake

#endif
END MODULE mo_seb_lake
