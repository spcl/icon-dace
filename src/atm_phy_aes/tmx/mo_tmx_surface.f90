!
! Interface for the surface component of the turbulent mixing package (tmx)
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!----------------------------
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!----------------------------

MODULE mo_tmx_surface_interface

  USE mo_kind, ONLY: wp, vp
  USE mo_exception, ONLY: finish
  USE mo_fortran_tools, ONLY: init
  USE mtime, ONLY: datetime
  USE mo_physical_constants, ONLY: rgrav, cpd, cvd, tmelt, Tf, stbo, rhos, alv, als, alf
  USE mo_turb_vdiff_params, ONLY: ckap
  USE mo_coupling_config,   ONLY: is_coupled_run
  USE mo_aes_phy_config,    ONLY: aes_phy_config  ! TODO: replace USE
  USE mo_tmx_field_class, ONLY: t_domain, isfc_oce, isfc_ice, isfc_lnd

  ! If land is present, JSBACH is currently the only surface scheme supported by AES physcis package
  ! USE mo_jsb_interface, ONLY: jsbach_interface
  USE mo_jsb_interface,     ONLY: jsbach_interface, jsbach_get_var

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_land, update_sea_ice, &
    & compute_sfc_sat_spec_humidity, compute_sfc_fluxes, compute_sfc_roughness, &
    & compute_lw_rad_net, compute_sw_rad_net, compute_albedo, &
    & compute_2m_temperature, compute_10m_wind

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_surface_interface'
  
CONTAINS

  SUBROUTINE update_land(jg, domain, datetime_old, dtime, cpd, &
    & dz, pres_srf, ptemp, pq, pres_air, rsfl, ssfl, &
    & rlds, rvds_dir, rnds_dir, rpds_dir, rvds_dif, rnds_dif, rpds_dif, &
    & cosmu0, wind, wind10m, rho, co2, &
    ! out
    & tsfc, tsfc_rad, tsfc_eff, q_snocpymlt, qsat, &
    & alb_vis_dir, &
    & alb_vis_dif, &
    & alb_nir_dir, &
    & alb_nir_dif, &
    & kh, km, kh_neutral, km_neutral)

    INTEGER, INTENT(in) :: &
      & jg
    TYPE(t_domain), INTENT(in), POINTER :: domain
    TYPE(datetime), INTENT(in), POINTER :: datetime_old ! date and time at beginning of this time step
    REAL(wp),       INTENT(in)          :: dtime
    REAL(wp),       INTENT(in)          :: cpd
    REAL(vp), INTENT(IN) :: &
      & dz        (:,:)     ! reference height in surface layer times 2
    REAL(wp), INTENT(IN) :: &
      & pres_srf  (:,:), &  ! surface pressure
      & ptemp     (:,:), &  ! temperature of lowest atmospheric level
      & pq        (:,:), &  ! humidity of lowest atmospheric level
      & pres_air  (:,:), &  ! pressure at lowest atmospheric level
      & rsfl      (:,:), &  ! surface rain flux, large-scale
      & ssfl      (:,:), &  ! surface snow flux, large-scale
      & rlds      (:,:), &
      & rvds_dir  (:,:), &
      & rnds_dir  (:,:), &
      & rpds_dir  (:,:), &
      & rvds_dif  (:,:), &
      & rnds_dif  (:,:), &
      & rpds_dif  (:,:), &
      & cosmu0    (:,:), &  ! cosine of zenith angle
      & wind      (:,:), &  ! wind speed at lowest level
      & wind10m   (:,:), &  ! wind speed at 10m
      & rho       (:,:), &  ! air density at surface
      & co2       (:,:)     ! CO2 at lowest atmospheric level

    REAL(wp), INTENT(out), OPTIONAL :: &
      & tsfc (:,:),       & ! new surface temperature
      & tsfc_rad (:,:),   & ! new surface radiative temperature
      & tsfc_eff (:,:),   & ! new surface effective temperature for rad heating
      & q_snocpymlt(:,:), & ! heat used to melt snow on canopy
      & qsat(:,:),        & ! saturated specific humidity at surface
      & alb_vis_dir(:,:), &
      & alb_vis_dif(:,:), &
      & alb_nir_dir(:,:), &
      & alb_nir_dif(:,:), &
      & kh         (:,:), & ! surface exchange coefficient (heat)
      & km         (:,:), & ! surface exchange coefficient (momentum)
      & kh_neutral (:,:), & ! neutral surface exchange coefficient (heat)
      & km_neutral (:,:)    ! neutral surface exchange coefficient (momentum)

    INTEGER :: jb, jc, jcs, jce
    REAL(wp), DIMENSION(domain%nproma) :: &
      & dz_srf, rain_tmp, snow_tmp, &
      & rvds, rnds, rpds, fract_par_diffuse, &
      & t_acoef, t_bcoef, q_acoef, q_bcoef

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_land'

!$OMP PARALLEL
    CALL init(tsfc)
    CALL init(qsat)
    IF (PRESENT(km)) THEN
      CALL init(tsfc_rad)
      CALL init(tsfc_eff)
      CALL init(q_snocpymlt)
      CALL init(alb_vis_dir)
      CALL init(alb_vis_dif)
      CALL init(alb_nir_dir)
      CALL init(alb_nir_dif)
      CALL init(kh)
      CALL init(km)
      CALL init(kh_neutral)
      CALL init(km_neutral)
    END IF
!$OMP END PARALLEL


    !$ACC DATA CREATE(dz_srf, rain_tmp, snow_tmp, rvds, rnds, rpds, fract_par_diffuse, t_acoef, t_bcoef, q_acoef, q_bcoef)

!$OMP PARALLEL DO PRIVATE(jb, jcs, jce, jc, dz_srf, rain_tmp, snow_tmp, rvds, rnds, rpds, fract_par_diffuse, &
!$OMP                     t_acoef, t_bcoef, q_acoef, q_bcoef) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c

      jcs = domain%i_startidx_c(jb)
      jce = domain%i_endidx_c(jb)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jc = jcs, jce
        dz_srf(jc) = dz(jc,jb) * 0.5_wp
        rain_tmp(jc) = rsfl(jc,jb)
        snow_tmp(jc) = ssfl(jc,jb)
        rvds(jc) = rvds_dir(jc,jb) + rvds_dif(jc,jb)
        rnds(jc) = rnds_dir(jc,jb) + rnds_dif(jc,jb)
        rpds(jc) = rpds_dir(jc,jb) + rpds_dif(jc,jb)
        IF (rpds(jc) > 0._wp) THEN
          fract_par_diffuse(jc) = rpds_dif(jc,jb) / rpds(jc)
        ELSE
          fract_par_diffuse(jc) = 0._wp
        END IF
        t_acoef(jc) = 0._wp
        t_bcoef(jc) = cpd * ptemp(jc,jb) ! dry static energy
        q_acoef(jc) = 0._wp
        q_bcoef(jc) = pq(jc,jb)
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC WAIT(1)

      IF (PRESENT(km)) THEN
        CALL jsbach_interface ( jg, jb, jcs, jce,                                         & ! in
          & datetime_old, dtime, dtime,                                                   & ! in
          & t_air             = ptemp(jcs:jce,jb),                                        & ! in
          & q_air             = pq(jcs:jce,jb),                                           & ! in
          & press_air         = pres_air(jcs:jce,jb),                                     & ! in
          & rain              = rain_tmp(jcs:jce),                                        & ! in
          & snow              = snow_tmp(jcs:jce),                                        & ! in
          & wind_air          = wind(jcs:jce,jb),                                         & ! in
          & wind_10m          = wind10m(jcs:jce,jb),                                      & ! in
          & lw_srf_down       = rlds(jcs:jce,jb),                                         & ! in
          & swvis_srf_down    = rvds(jcs:jce),                                            & ! in
          & swnir_srf_down    = rnds(jcs:jce),                                            & ! in
          & swpar_srf_down    = rpds(jcs:jce),                                            & ! in
          & fract_par_diffuse = fract_par_diffuse(jcs:jce),                               & ! in
          & dz_srf            = dz_srf(jcs:jce),                                          & ! in
          & press_srf         = pres_srf(jcs:jce,jb),                                     & ! in
          & rho_srf           = rho(jcs:jce,jb),                                          & ! in
          & t_acoef           = t_acoef(jcs:jce),                                         & ! in
          & t_bcoef           = t_bcoef(jcs:jce),                                         & ! in
          & q_acoef           = q_acoef(jcs:jce),                                         & ! in
          & q_bcoef           = q_bcoef(jcs:jce),                                         & ! in
          & cos_zenith_angle  = cosmu0(jcs:jce,jb),                                       & ! in
          & CO2_air           = co2(jcs:jce,jb),                                          & ! in
          & t_srf             = tsfc(jcs:jce,jb),                                         & ! out
          & t_srf_rad         = tsfc_rad(jcs:jce,jb),                                     & ! out
          & t_srf_eff         = tsfc_eff(jcs:jce,jb),                                     & ! out
          & q_snocpymlt       = q_snocpymlt(jcs:jce,jb),                                  & ! out
          & qsat_srf          = qsat(jcs:jce,jb),                                         & ! out
          & alb_vis_dir       = alb_vis_dir(jcs:jce,jb),                                  & ! out
          & alb_nir_dir       = alb_nir_dir(jcs:jce,jb),                                  & ! out
          & alb_vis_dif       = alb_vis_dif(jcs:jce,jb),                                  & ! out
          & alb_nir_dif       = alb_nir_dif(jcs:jce,jb),                                  & ! out
          & kh                = kh(jcs:jce,jb),                                           & ! out
          & km                = km(jcs:jce,jb),                                           & ! out
          & kh_neutral        = kh_neutral(jcs:jce,jb),                                   & ! out
          & km_neutral        = km_neutral(jcs:jce,jb)                                    & ! out
          ! & t_eff_srf         = ztsfc_lnd_eff(jcs:jce),                                   & ! out (T_s^eff) surface temp 
          !                                                                                     ! (effective, for longwave rad)
          ! & s_srf             = zcpt_lnd(jcs:jce),                                        & ! out (s_s^star, for vdiff scheme)
          ! & fact_q_air        = pcair(jcs:jce),                                           & ! out
          ! & fact_qsat_srf     = pcsat(jcs:jce),                                           & ! out
          ! & evapotrans        = zevap_lnd(jcs:jce),                                       & ! out
          ! & latent_hflx       = zlhflx_lnd(jcs:jce),                                      & ! out
          ! & sensible_hflx     = zshflx_lnd(jcs:jce),                                      & ! out
          ! & grnd_hflx         = zgrnd_hflx(jcs:jce, idx_lnd),                             & ! out
          ! & grnd_hcap         = zgrnd_hcap(jcs:jce, idx_lnd),                             & ! out
          ! & rough_h_srf       = z0h_lnd(jcs:jce),                                         & ! out
          ! & rough_m_srf       = z0m_tile(jcs:jce, idx_lnd),                               & ! out
          ! & q_snocpymlt       = q_snocpymlt(jcs:jce),                                     & ! out
          ! & co2_flux          = pco2_flux_tile(jcs:jce, idx_lnd)                          & ! out
        )
      ELSE
        CALL jsbach_interface ( jg, jb, jcs, jce,                                         & ! in
          & datetime_old, dtime, dtime,                                                   & ! in
          & t_air             = ptemp(jcs:jce,jb),                                        & ! in
          & q_air             = pq(jcs:jce,jb),                                           & ! in
          & press_air         = pres_air(jcs:jce,jb),                                     & ! in
          & rain              = rain_tmp(jcs:jce),                                        & ! in
          & snow              = snow_tmp(jcs:jce),                                        & ! in
          & wind_air          = wind(jcs:jce,jb),                                         & ! in
          & wind_10m          = wind10m(jcs:jce,jb),                                      & ! in
          & lw_srf_down       = rlds(jcs:jce,jb),                                         & ! in
          & swvis_srf_down    = rvds(jcs:jce),                                            & ! in
          & swnir_srf_down    = rnds(jcs:jce),                                            & ! in
          & swpar_srf_down    = rpds(jcs:jce),                                            & ! in
          & fract_par_diffuse = fract_par_diffuse(jcs:jce),                               & ! in
          & dz_srf            = dz_srf(jcs:jce),                                          & ! in
          & press_srf         = pres_srf(jcs:jce,jb),                                     & ! in
          & rho_srf           = rho(jcs:jce,jb),                                          & ! in
          & t_acoef           = t_acoef(jcs:jce),                                         & ! in
          & t_bcoef           = t_bcoef(jcs:jce),                                         & ! in
          & q_acoef           = q_acoef(jcs:jce),                                         & ! in
          & q_bcoef           = q_bcoef(jcs:jce),                                         & ! in
          & cos_zenith_angle  = cosmu0(jcs:jce,jb),                                       & ! in
          & CO2_air           = co2(jcs:jce,jb),                                          & ! in
          & t_srf             = tsfc(jcs:jce,jb),                                         & ! out
          & qsat_srf          = qsat(jcs:jce,jb)                                          & ! out
        )
      END IF
  
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

    !$ACC END DATA


  END SUBROUTINE update_land

  SUBROUTINE update_sea_ice(domain, dtime, cpd, &
    & old_tsfc, &
    & lwflx_net, swflx_net, lhflx, shflx, &
    & ssfl, ice_thickness, &
    & emissivity, &
    ! inout &
    & snow_thickness, &
    ! out &
    & new_tsfc, q_top, q_bot, &
    & albvisdir, albvisdif, albnirdir, albnirdif)

  USE mo_ice_interface, ONLY: ice_fast
  
    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), INTENT(in) :: dtime
    REAL(wp), INTENT(in) :: cpd
    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & old_tsfc,       &
      & lwflx_net,      &
      & swflx_net,      &
      & lhflx,          &
      & shflx,          &
      & ssfl,           &
      & ice_thickness,  &
      & emissivity

    REAL(wp), INTENT(inout), DIMENSION(:,:) :: &
      & snow_thickness

    REAL(wp), INTENT(out), DIMENSION(:,:) :: &
      & new_tsfc,  &  ! new surface temperature
      & q_top,     &
      & q_bot,     &
      & albvisdir, &
      & albvisdif, &
      & albnirdir, &
      & albnirdif

    INTEGER :: jb, jc, jcs, jce, kice
    REAL(wp), DIMENSION(domain%nproma) :: &
      & Tfw, nonsolar_flux, dnonsolar_flux_dt, &
      & T1, T2 !< Dummies, not used

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_sea_ice'

!$OMP PARALLEL
    CALL init(new_tsfc)
    CALL init(q_top)
    CALL init(q_bot)
    CALL init(albvisdir)
    CALL init(albvisdif)
    CALL init(albnirdir)
    CALL init(albnirdif)
    CALL init(T1)
    CALL init(T2)
!$OMP END PARALLEL


    kice = 1

    !$ACC DATA CREATE(Tfw, nonsolar_flux, dnonsolar_flux_dt, T1, T2)

!$OMP PARALLEL DO PRIVATE(jb, jcs, jce, jc, Tfw, nonsolar_flux, dnonsolar_flux_dt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      jcs = domain%i_startidx_c(jb)
      jce = domain%i_endidx_c(jb)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jc = jcs, jce
        Tfw(jc) = Tf
        nonsolar_flux(jc) = lwflx_net(jc,jb) + lhflx(jc,jb) + shflx(jc,jb)
        dnonsolar_flux_dt(jc) = -4._wp * emissivity(jc,jb) * stbo * old_tsfc(jc,jb)**3

        new_tsfc(jc,jb) = old_tsfc(jc,jb) - tmelt
      END DO
      !$ACC END PARALLEL LOOP

      CALL ice_fast(jcs, jce, domain%nproma, kice, dtime, &
        &   new_tsfc(:,jb),        & ! inout
        &   T1(:),                 & ! inout, dummy
        &   T2(:),                 & ! inout, dummy
        &   ice_thickness(:,jb),   & ! in
        &   snow_thickness(:,jb),  & ! in
        &   q_top(:,jb),           & ! out
        &   q_bot(:,jb),           & ! out
        &   swflx_net(:,jb),       & ! in
        &   nonsolar_flux(:),      & ! in
        &   dnonsolar_flux_dt(:),  & ! in
        &   Tfw(:),                & ! in
        &   albvisdir(:,jb),       & ! out
        &   albvisdif(:,jb),       & ! out
        &   albnirdir(:,jb),       & ! out
        &   albnirdif(:,jb),       & ! out
        &   lacc=.TRUE.)             ! in

      ! Update the thickness of snow on ice in atmosphere only simulation.
      ! In coupled experiments this is done by the ocean model in either
      ! ice_growth_zerolayer or ice_growth_winton.
      IF ( .NOT. is_coupled_run() ) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO jc = jcs, jce
          ! Snowfall on ice - no ice => no snow
          IF (ice_thickness(jc,jb) > 0._wp) THEN
            ! Snow only falls when it's below freezing
            IF (new_tsfc(jc,jb) < 0._wp) THEN
              snow_thickness(jc,jb) = snow_thickness(jc,jb) + ssfl(jc,jb) * dtime / rhos
            ENDIF
            ! Snow melt
            snow_thickness(jc,jb) = snow_thickness(jc,jb) - MIN( q_top(jc,jb) * dtime / (alf * rhos), snow_thickness(jc,jb) )
          ELSE
            snow_thickness(jc,jb) = 0._wp
          ENDIF
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jc=jcs,jce
        new_tsfc(jc,jb) = new_tsfc(jc,jb) + tmelt
      END DO
      !$ACC END PARALLEL LOOP

    END DO

!$OMP END PARALLEL DO

    !$ACC WAIT(1)

    !$ACC END DATA



  
  END SUBROUTINE update_sea_ice

  SUBROUTINE compute_lw_rad_net( &
    & domain,                    &
    & nvalid, indices,           &
    & emissivity,                &
    & rlds,                      &
    & tsfc,                      &
    & lwfl_net                   &
    & )

    USE mo_physical_constants,ONLY: stbo

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & emissivity ,     &
      & rlds,            &
      & tsfc
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: lwfl_net

    INTEGER  :: jb, jl, jls, js

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_lw_rad_net'

!$OMP PARALLEL
    CALL init(lwfl_net)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        lwfl_net(js,jb) = emissivity(js,jb) * (rlds(js,jb) - stbo * tsfc(js,jb)**4)
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_lw_rad_net
  !
  !=================================================================
  !
  SUBROUTINE compute_sw_rad_net( &
    & domain,                    &
    & nvalid, indices,           &
    & rvds_dir,                  &
    & rvds_dif,                  &
    & rnds_dir,                  &
    & rnds_dif,                  &
    & alb_vis_dir,               &
    & alb_vis_dif,               &
    & alb_nir_dir,               &
    & alb_nir_dif,               &
    & swfl_net                   &
    & )

    USE mo_physical_constants,ONLY: stbo

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & rvds_dir,            &
      & rvds_dif,            &
      & rnds_dir,            &
      & rnds_dif,            &
      & alb_vis_dir,         &
      & alb_vis_dif,         &
      & alb_nir_dir,         &
      & alb_nir_dif
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: swfl_net

    INTEGER  :: jb, jl, jls, js

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sw_rad_net'

!$OMP PARALLEL
    CALL init(swfl_net)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        swfl_net(js,jb) = &
          & rvds_dif(js,jb) * (1._wp - alb_vis_dif(js,jb)) + &
          & rvds_dir(js,jb) * (1._wp - alb_vis_dir(js,jb)) + &
          & rnds_dif(js,jb) * (1._wp - alb_nir_dif(js,jb)) + &
          & rnds_dir(js,jb) * (1._wp - alb_nir_dir(js,jb))
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_sw_rad_net
  !
  !=================================================================
  !
  SUBROUTINE compute_sfc_roughness( &
    & linit,                   &
    & domain,                  &
    & isfc,                    &
    & nvalid, indices,         &
    & rough_min, rough_oce, rough_ice, wind, km,  &
    & rough_h, rough_m         &
    & )

    USE mo_turb_vdiff_params, ONLY: cchar
  
    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    INTEGER,         INTENT(in) :: isfc
    !
    ! Input variables
    !
    LOGICAL :: linit
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), INTENT(in) :: rough_min, rough_oce, rough_ice
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      km ,     &
      wind
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: rough_h, rough_m

    INTEGER  :: jb, jl, jls, js

    REAL(wp), POINTER :: jsb_rough_m(:,:) => NULL(), jsb_rough_h(:,:) => NULL()
    REAL(wp) :: rough_tmp

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_roughness'

!$OMP PARALLEL
    CALL init(rough_h)
    CALL init(rough_m)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jl, jls, js, rough_tmp) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR PRIVATE(js, rough_tmp) ASYNC(1)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)

        IF (isfc == isfc_oce) THEN
          ! rough_m(js,jb) = MERGE(rough_oce, wind(js,jb) * km(js,jb) * cchar * rgrav, linit)
          ! rough_m(js,jb) = MAX(rough_min, rough_m(js,jb))
          rough_tmp = MERGE(rough_oce, wind(js,jb) * km(js,jb) * cchar * rgrav, linit)
          rough_tmp = MAX(rough_min, rough_tmp)
          rough_m(js,jb) = rough_tmp
          ! rough_m(js,jb) = 1.E-3_wp
          rough_h(js,jb) = rough_tmp
          ! rough_h(js,jb) = rough_m(js,jb)
        ELSE IF (isfc == isfc_ice) THEN
          ! Nothing to do
          rough_m(js,jb) = rough_ice
          rough_h(js,jb) = rough_ice
          ! rough_h(js,jb) = rough_m(js,jb)
        END IF

      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

    IF (isfc == isfc_lnd) THEN

      CALL jsbach_get_var('turb_rough_m', 1, ptr2d=jsb_rough_m)
      CALL jsbach_get_var('turb_rough_h', 1, ptr2d=jsb_rough_h)
!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c,domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR PRIVATE(js) ASYNC(1)
        DO jls = 1, nvalid(jb)
          js = indices(jls,jb)
          rough_m(js,jb) = jsb_rough_m(js,jb)
          rough_h(js,jb) = jsb_rough_h(js,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO
      NULLIFY(jsb_rough_m, jsb_rough_h)



    END IF

    !$ACC WAIT(1)

  END SUBROUTINE compute_sfc_roughness
  !
  !=================================================================
  !
  SUBROUTINE compute_sfc_sat_spec_humidity( &
    & linit, domain, isfc,     &
    & nvalid, indices,         &
    & ppsfc, ptsfc,            &
    & qsat                     &
    )

    USE mo_aes_thermo, ONLY: sat_pres_water, sat_pres_ice, specific_humidity

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    INTEGER,         INTENT(in) :: isfc
    !
    ! Input variables
    !
    LOGICAL, INTENT(in) :: linit
    INTEGER, INTENT(in) :: &
      & nvalid(:),         &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      ptsfc ,     &
      ppsfc
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: qsat

    INTEGER  :: jb, jls, js
    REAL(wp), POINTER :: jsb_qsat(:,:) => NULL()
    
    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_sat_spec_humidity'





!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c

      ! DO jl = 1, domain%nproma
      !   qsat(jl,jb) = 0._wp
      ! END DO

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        IF (isfc == isfc_oce) THEN
          qsat(js,jb) = specific_humidity(sat_pres_water(ptsfc(js,jb)),ppsfc(js,jb))
        ELSE IF (isfc == isfc_ice) THEN
          qsat(js,jb) = specific_humidity(sat_pres_ice(ptsfc(js,jb)),ppsfc(js,jb))
        ELSE IF (isfc == isfc_lnd .AND. linit) THEN
          qsat(js,jb) = specific_humidity(sat_pres_water(ptsfc(js,jb)),ppsfc(js,jb))
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO


    IF (isfc == isfc_lnd .AND. .NOT. linit) THEN
      CALL jsbach_get_var('seb_qsat_star', 1, ptr2d=jsb_qsat)
!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c,domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
        DO jls = 1, nvalid(jb)
          js = indices(jls,jb)
          qsat(js,jb) = jsb_qsat(js,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO
      NULLIFY(jsb_qsat)
    END IF


    !$ACC WAIT(1)

  END SUBROUTINE compute_sfc_sat_spec_humidity
  !
  !=================================================================
  !
  SUBROUTINE compute_sfc_fluxes( &
    & domain,                  &
    & isfc,                    &
    & nvalid, indices,         &
    & cpd, cvd,                &
    ! & ua, va, thetam1, qm1, wind, rho, qsat_sfc, theta_sfc, kh, km,  &
    & ua, va, ta, qm1, wind, u_sfc_oce, v_sfc_oce, rho, qsat_sfc, t_sfc, kh, km,  &
    & evapotrans, latent_hflx, sensible_hflx, ustress, vstress  &
    & )

    USE mo_turb_vdiff_params, ONLY: cchar
    USE mo_nh_testcases_nml,  ONLY: isrfc_type, shflx, lhflx
  
    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    INTEGER,         INTENT(in) :: isfc
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), INTENT(in) :: cpd, cvd
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & ua, &
      & va, &
      ! & thetam1, &
      & ta, &
      & qm1, &
      & kh ,     &
      & km,      &
      & rho,     &
      & qsat_sfc, &
      ! & theta_sfc, &
      & t_sfc, &
      & wind, &
      & u_sfc_oce, &
      & v_sfc_oce
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: evapotrans, latent_hflx, sensible_hflx, ustress, vstress

    INTEGER  :: jb, jl, jls, js

    ! REAL ::
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      &jsb_evapotrans_ptr => NULL(), jsb_latent_hflx_ptr => NULL(), jsb_sensible_hflx_ptr => NULL()

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_fluxes'

    ! CALL message(routine, 'Start')

!$OMP PARALLEL
      CALL init(evapotrans)
      CALL init(latent_hflx)
      CALL init(sensible_hflx)
      CALL init(ustress)
      CALL init(vstress)
!$OMP END PARALLEL

    IF (isrfc_type == 1) THEN
!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c,domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
        DO jls = 1, nvalid(jb)
          js = indices(jls,jb)
          latent_hflx(js,jb)   = -lhflx * alv * rho(js,jb)
          evapotrans(js,jb)    = -lhflx * rho(js,jb)
          sensible_hflx(js,jb) = -shflx * cvd * rho(js,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO
      !$ACC WAIT(1)
      RETURN
    END IF





    IF (isfc == isfc_lnd) THEN

      CALL jsbach_get_var('hydro_evapotrans',  1, ptr2d=jsb_evapotrans_ptr)
      CALL jsbach_get_var('seb_latent_hflx',   1, ptr2d=jsb_latent_hflx_ptr)
      CALL jsbach_get_var('seb_sensible_hflx', 1, ptr2d=jsb_sensible_hflx_ptr)

    END IF

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
      DO jls = 1, nvalid(jb)  
        js = indices(jls,jb)
        ! TODO: is the treatment of surface ocean current correct (cf. vdiff code)
        IF (isfc == isfc_oce) THEN
          evapotrans(js,jb) = rho(js,jb) * wind(js,jb) * kh(js,jb) * (qm1(js,jb) - qsat_sfc(js,jb))
          latent_hflx(js,jb) = alv * evapotrans(js,jb)
          sensible_hflx(js,jb) = cvd * rho(js,jb) * wind(js,jb) * kh(js,jb) * (ta(js,jb) - t_sfc(js,jb))
          ustress(js,jb) = rho(js,jb) * km(js,jb) * wind(js,jb) * (ua(js,jb) - u_sfc_oce(js,jb))
          vstress(js,jb) = rho(js,jb) * km(js,jb) * wind(js,jb) * (va(js,jb) - v_sfc_oce(js,jb))
        ELSE IF (isfc == isfc_ice) THEN
          evapotrans(js,jb) = rho(js,jb) * wind(js,jb) * kh(js,jb) * (qm1(js,jb) - qsat_sfc(js,jb))
          latent_hflx(js,jb) = als * evapotrans(js,jb)
          sensible_hflx(js,jb) = cvd * rho(js,jb) * wind(js,jb) * kh(js,jb) * (ta(js,jb) - t_sfc(js,jb))
          ustress(js,jb) = rho(js,jb) * km(js,jb) * wind(js,jb) * ua(js,jb)
          vstress(js,jb) = rho(js,jb) * km(js,jb) * wind(js,jb) * va(js,jb)
        ELSE IF (isfc == isfc_lnd) THEN
          evapotrans(js,jb) = jsb_evapotrans_ptr(js,jb)
          latent_hflx(js,jb) = jsb_latent_hflx_ptr(js,jb)
          sensible_hflx(js,jb) = jsb_sensible_hflx_ptr(js,jb)
          ustress(js,jb) = rho(js,jb) * km(js,jb) * wind(js,jb) * ua(js,jb)
          vstress(js,jb) = rho(js,jb) * km(js,jb) * wind(js,jb) * va(js,jb)
        END IF
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

    IF (isfc == isfc_lnd) THEN
      NULLIFY(jsb_evapotrans_ptr, jsb_latent_hflx_ptr, jsb_sensible_hflx_ptr)
    END IF

  END SUBROUTINE compute_sfc_fluxes
  !
  !=================================================================
  !
  SUBROUTINE compute_albedo(     &
    & domain,                    &
    & nvalid, indices,           &
    & rsds,                      &
    & rvds_dir,                  &
    & rvds_dif,                  &
    & rnds_dir,                  &
    & rnds_dif,                  &
    & alb_vis_dir,               &
    & alb_vis_dif,               &
    & alb_nir_dir,               &
    & alb_nir_dif,               &
    & albedo                     &
    & )

    USE mo_physical_constants,ONLY: stbo

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & rsds,                &
      & rvds_dir,            &
      & rvds_dif,            &
      & rnds_dir,            &
      & rnds_dif,            &
      & alb_vis_dir,         &
      & alb_vis_dif,         &
      & alb_nir_dir,         &
      & alb_nir_dif
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: albedo

    INTEGER  :: jb, jls, js
    REAL(wp) :: zalbvis, zalbnir, rvds, rnds

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_albedo'

!$OMP PARALLEL
    CALL init(albedo)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js, zalbvis, zalbnir, rvds, rnds) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)

        zalbvis = 0._wp
        rvds = rvds_dir(js,jb) + rvds_dif(js,jb)
        IF(rvds > 0._wp) THEN
          zalbvis = &
            & (alb_vis_dir(js,jb) * rvds_dir(js,jb) + alb_vis_dif(js,jb) * rvds_dif(js,jb)) &
            & / rvds
        END IF
        zalbnir = 0._wp
        rnds = rnds_dir(js,jb) + rnds_dif(js,jb)
        IF(rnds > 0._wp) THEN
          zalbnir = &
            & (alb_nir_dir(js,jb) * rnds_dir(js,jb) + alb_nir_dif(js,jb) * rnds_dif(js,jb)) &
            & / rnds
        END IF
        IF(rsds(js,jb) > 0._wp) THEN
          albedo(js,jb) = &
            & (zalbvis * rvds + zalbnir * rnds) &
            & / rsds(js,jb)
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_albedo
  !
  !=================================================================
  !
  SUBROUTINE compute_10m_wind( &
    & domain, isfc,               &
    & nvalid, indices, zf, zh,    &
    & ua, va, u_oce, v_oce,       &
    & moist_rich, km, km_neutral, &
    & u10m, v10m, wind10m)

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:),        &
      & isfc
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & zf, zh, &
      & ua, va, u_oce, v_oce, moist_rich, km, km_neutral
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & u10m, v10m, wind10m

    INTEGER :: jb, jls, js
    REAL(wp) :: zrat, zbm, zcbn, zcbs, zcbu, zmerge, zred

    ! to prevent floating-point arithmetic inconsistencies later in
    ! the interpolation to u 10m and 2m T/T_d: has been 0.01 before
    ! (Cray FP instead of IEEE 754 FP format)
    REAL(wp), PARAMETER :: zepsec = 0.028_wp

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_10m_wind'

!$OMP PARALLEL
    CALL init(u10m)
    CALL init(v10m)
    CALL init(wind10m)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js, zrat, zbm, zcbn, zcbs, zcbu, zmerge, zred) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js, zrat, zbm, zcbn, zcbs, zcbu, zmerge, zred)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        zbm    = 1._wp / MAX(zepsec, SQRT(km(js,jb)) / ckap)
        zrat   = 10._wp / (zf(js,jb) - zh(js,jb))
        zcbn   = LOG(1._wp + (EXP (km_neutral(js,jb)) - 1._wp) * zrat )
        zcbs   = -(km_neutral(js,jb) - zbm) * zrat
        zcbu   = -LOG(1._wp + (EXP (km_neutral(js,jb) - zbm) - 1._wp) * zrat)
        zmerge = MERGE(zcbs, zcbu, moist_rich(js,jb) > 0._wp)
        zred   = (zcbn + zmerge) / zbm
        u10m(js,jb) = zred * ua(js,jb)
        v10m(js,jb) = zred * va(js,jb)
        IF (isfc == isfc_oce) THEN
          wind10m(js,jb) = zred * SQRT((ua(js,jb) - u_oce(js,jb))**2 + (va(js,jb) - v_oce(js,jb))**2)
        ELSE
          wind10m(js,jb) = SQRT(u10m(js,jb)**2 + (v10m(js,jb)**2))
        END IF
      END DO
    !$ACC END PARALLEL LOOP
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_10m_wind
  !
  !=================================================================
  !
  SUBROUTINE compute_2m_temperature(              &
    & domain, isfc,                               &
    & nvalid, indices, zf, zh,                    &
    & tatm, tsfc,                                 &
    & moist_rich, kh, km, kh_neutral, km_neutral, &
    & t2m)

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:),        &
      & isfc
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & zf, zh, &
      & tatm, tsfc, moist_rich, kh, km, kh_neutral, km_neutral
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & t2m

    INTEGER :: jb, jls, js
    REAL(wp) :: zrat, zbm, zbh, zcbn, zcbs, zcbu, zmerge, zred

    ! to prevent floating-point arithmetic inconsistencies later in
    ! the interpolation to u 10m and 2m T/T_d: has been 0.01 before
    ! (Cray FP instead of IEEE 754 FP format)
    REAL(wp), PARAMETER :: zepsec = 0.028_wp

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_2m_temperature'

!$OMP PARALLEL
    CALL init(t2m)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js, zrat, zbm, zbh, zcbn, zcbs, zcbu, zmerge, zred) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js, zrat, zbm, zcbn, zcbs, zcbu, zmerge, zred)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        zrat = 2._wp / (zf(js,jb) - zh(js,jb))
        zbm  = 1._wp / MAX(zepsec, SQRT(km(js,jb)) / ckap)
        zbh  = 1._wp / MAX(zepsec, kh(js,jb) * zbm / ckap**2)
        IF (isfc == isfc_lnd) THEN
          zcbn   = LOG(1._wp + (EXP (kh_neutral(js,jb)) - 1._wp) * zrat )
          zcbs   = -(kh_neutral(js,jb) - zbh) * zrat
          zcbu   = -LOG(1._wp + (EXP (kh_neutral(js,jb) - zbh) - 1._wp) * zrat)
        ELSE
          zcbn   = LOG(1._wp + (EXP (km_neutral(js,jb)) - 1._wp) * zrat )
          zcbs   = -(km_neutral(js,jb) - zbh) * zrat
          zcbu   = -LOG(1._wp + (EXP (km_neutral(js,jb) - zbh) - 1._wp) * zrat)
        END IF
        zmerge = MERGE(zcbs, zcbu, moist_rich(js,jb) > 0._wp)
        zred   = (zcbn + zmerge) / zbh
        t2m(js,jb) = tsfc(js,jb) + zred * (tatm(js,jb) - tsfc(js,jb))
      END DO
    !$ACC END PARALLEL LOOP
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_2m_temperature
  !
  !=================================================================
  !
  SUBROUTINE compute_2m_dewpoint(              &
    & domain, isfc,                               &
    & nvalid, indices, zf, zh,                    &
    & tatm, tsfc,                                 &
    & moist_rich, kh, km, kh_neutral, km_neutral, &
    & t2m)

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:),        &
      & isfc
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & zf, zh, &
      & tatm, tsfc, moist_rich, kh, km, kh_neutral, km_neutral
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & t2m

    INTEGER :: jb, jls, js
    REAL(wp) :: zrat, zbm, zbh, zcbn, zcbs, zcbu, zmerge, zred

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_2m_dewpoint'

!$OMP PARALLEL DO PRIVATE(jb, jls, js, zrat, zbm, zbh, zcbn, zcbs, zcbu, zmerge, zred) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js, zrat, zbm, zcbn, zcbs, zcbu, zmerge, zred)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
      END DO
      !$ACC END PARALLEL LOOP
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_2m_dewpoint
  !
  !=================================================================
  !

END MODULE mo_tmx_surface_interface
