# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/tmx/mo_tmx_surface.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/tmx/mo_tmx_surface.f90"
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

! Interface for the surface component of the turbulent mixing package (tmx)

!----------------------------

# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/include/omp_definitions.inc" 1
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

# 16 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/tmx/mo_tmx_surface.f90" 2
!----------------------------

MODULE mo_tmx_surface_interface

  USE mo_kind, ONLY: wp, vp
  USE mo_exception, ONLY: finish
  USE mo_fortran_tools, ONLY: init
  USE mtime, ONLY: datetime
  USE mo_physical_constants, ONLY: grav, rgrav, tmelt, Tf, stbo, rhos, alf, cvv, clw, ci, vtmpc1, rd
  USE mo_aes_thermo, ONLY: &
    & lvc, lsc, &
    & sat_pres_water, sat_pres_ice, specific_humidity, dewpoint_temperature
  USE mo_turb_vdiff_params, ONLY: ckap
  USE mo_coupling_config,   ONLY: is_coupled_to_ocean
  USE mo_aes_phy_config,    ONLY: aes_phy_config  ! TODO: replace USE
  USE mo_tmx_field_class, ONLY: t_domain, isfc_oce, isfc_ice, isfc_lnd

  ! If land is present, JSBACH is currently the only surface scheme supported by AES physcis package





  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_land, update_sea_ice, &
    & compute_sfc_sat_spec_humidity, compute_sfc_fluxes, compute_sfc_roughness, &
    & compute_lw_rad_net, compute_sw_rad_net, compute_albedo, compute_energy_fluxes, &
    & compute_2m_temperature, compute_2m_humidity_and_dewpoint, compute_10m_wind

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_surface_interface'
  
CONTAINS

  SUBROUTINE update_land(jg, domain, datetime_old, dtime, cvd, &
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
    REAL(wp),       INTENT(in)          :: cvd
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
    CALL init(tsfc, lacc=.TRUE.)
    CALL init(qsat, lacc=.TRUE.)
    IF (PRESENT(km)) THEN
      CALL init(tsfc_rad, lacc=.TRUE.)
      CALL init(tsfc_eff, lacc=.TRUE.)
      CALL init(q_snocpymlt, lacc=.TRUE.)
      CALL init(alb_vis_dir, lacc=.TRUE.)
      CALL init(alb_vis_dif, lacc=.TRUE.)
      CALL init(alb_nir_dir, lacc=.TRUE.)
      CALL init(alb_nir_dif, lacc=.TRUE.)
      CALL init(kh, lacc=.TRUE.)
      CALL init(km, lacc=.TRUE.)
      CALL init(kh_neutral, lacc=.TRUE.)
      CALL init(km_neutral, lacc=.TRUE.)
    END IF
!$OMP END PARALLEL

# 254 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/tmx/mo_tmx_surface.f90"
    CALL finish(routine, "The JSBACH component is not activated")


  END SUBROUTINE update_land

  SUBROUTINE update_sea_ice(domain, dtime, &
    & old_tsfc, &
    & lwflx_net, swflx_net, lhflx, shflx, &
    & ssfl, ice_thickness, &
    & emissivity, &
    ! inout &
    & snow_thickness, &
    ! out &
    & new_tsfc, q_top, q_bot, &
    & albvisdir, albvisdif, albnirdir, albnirdif)




  
    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), INTENT(in) :: dtime
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

    !$ACC DATA CREATE(Tfw, nonsolar_flux, dnonsolar_flux_dt, T1, T2)

!$OMP PARALLEL
    CALL init(new_tsfc, lacc=.TRUE.)
    CALL init(q_top, lacc=.TRUE.)
    CALL init(q_bot, lacc=.TRUE.)
    CALL init(albvisdir, lacc=.TRUE.)
    CALL init(albvisdif, lacc=.TRUE.)
    CALL init(albnirdir, lacc=.TRUE.)
    CALL init(albnirdif, lacc=.TRUE.)
    CALL init(T1, lacc=.TRUE.)
    CALL init(T2, lacc=.TRUE.)
!$OMP END PARALLEL

# 389 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/tmx/mo_tmx_surface.f90"
    CALL finish(routine, "The ice process requires the ICON_OCEAN component")


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
    CALL init(lwfl_net, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
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
    CALL init(swfl_net, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
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
    CALL init(rough_h, lacc=.TRUE.)
    CALL init(rough_m, lacc=.TRUE.)
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
# 599 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/tmx/mo_tmx_surface.f90"
      CALL finish(routine, "The JSBACH component is not activated")

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


    IF (isfc == isfc_lnd) CALL finish(routine, "The JSBACH component is not activated")


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

# 681 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/tmx/mo_tmx_surface.f90"

    !$ACC WAIT(1)

  END SUBROUTINE compute_sfc_sat_spec_humidity
  !
  !=================================================================
  !
  SUBROUTINE compute_sfc_fluxes( &
    & domain,                  &
    & isfc,                    &
    & nvalid, indices,         &
    & cvd,                     &
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
    REAL(wp), INTENT(in) :: cvd
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
      CALL init(evapotrans, lacc=.TRUE.)
      CALL init(latent_hflx, lacc=.TRUE.)
      CALL init(sensible_hflx, lacc=.TRUE.)
      CALL init(ustress, lacc=.TRUE.)
      CALL init(vstress, lacc=.TRUE.)
!$OMP END PARALLEL

    IF (isrfc_type == 1) THEN
!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c,domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
        DO jls = 1, nvalid(jb)
          js = indices(jls,jb)
          latent_hflx(js,jb)   = -lhflx * rho(js,jb) * (lvc+(cvv-clw)*t_sfc(js,jb))
          evapotrans(js,jb)    = -lhflx * rho(js,jb)
          sensible_hflx(js,jb) = -shflx * cvd * rho(js,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO
      !$ACC WAIT(1)
      RETURN
    END IF


    IF (isfc == isfc_lnd) CALL finish(routine, "The JSBACH component is not activated")


    IF (isfc == isfc_lnd) THEN





    END IF

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
      DO jls = 1, nvalid(jb)  
        js = indices(jls,jb)
        ! TODO: is the treatment of surface ocean current correct (cf. vdiff code)
        IF (isfc == isfc_oce) THEN
          evapotrans(js,jb) = rho(js,jb) * wind(js,jb) * kh(js,jb) * (qm1(js,jb) - qsat_sfc(js,jb))
          latent_hflx(js,jb) = evapotrans(js,jb) * (lvc+(cvv-clw)*t_sfc(js,jb))
          sensible_hflx(js,jb) = cvd * rho(js,jb) * wind(js,jb) * kh(js,jb) * (ta(js,jb) - t_sfc(js,jb))
          ustress(js,jb) = rho(js,jb) * km(js,jb) * wind(js,jb) * (ua(js,jb) - u_sfc_oce(js,jb))
          vstress(js,jb) = rho(js,jb) * km(js,jb) * wind(js,jb) * (va(js,jb) - v_sfc_oce(js,jb))
        ELSE IF (isfc == isfc_ice) THEN
          evapotrans(js,jb) = rho(js,jb) * wind(js,jb) * kh(js,jb) * (qm1(js,jb) - qsat_sfc(js,jb))
          latent_hflx(js,jb) = evapotrans(js,jb) * (lsc+(cvv-ci)*t_sfc(js,jb))
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
  SUBROUTINE compute_energy_fluxes( &
    & domain,     &
    & cvv, cvd,   &
    & shfl,       &
    & evapotrans, &
    & ta,         &
    & rho,        &
    & ufts,       &
    & ufvs)

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    REAL(wp), INTENT(in) :: cvv, cvd
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & shfl,       &
      & evapotrans, &
      & ta,         &
      & rho
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & ufts, &
      & ufvs

    INTEGER :: ic, ib

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_energy_fluxes'

!$OMP PARALLEL
    CALL init(ufts, lacc=.TRUE.)
    CALL init(ufvs, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(ib, ic) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO ic = domain%i_startidx_c(ib), domain%i_endidx_c(ib)
        ufts(ic,ib) = shfl(ic,ib)
        ufvs(ic,ib) = ta(ic,ib) * evapotrans(ic,ib) * (cvv - cvd) 
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_energy_fluxes
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
    CALL init(albedo, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js, zalbvis, zalbnir, rvds, rnds) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
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
    CALL init(u10m, lacc=.TRUE.)
    CALL init(v10m, lacc=.TRUE.)
    CALL init(wind10m, lacc=.TRUE.)
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
    CALL init(t2m, lacc=.TRUE.)
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

    !$ACC WAIT(1)

  END SUBROUTINE compute_2m_temperature
  !
  !=================================================================
  !
  SUBROUTINE compute_2m_humidity_and_dewpoint( &
    & domain, nvalid, indices,    &
    & patm, psfc,                 &
    & tatm, t2m,                  &
    & qv_atm, qc_atm, qi_atm,     &
    & hus2m, dew2m)

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & patm, psfc,             & ! pressure in lowest model level and at surface
      & tatm, t2m,              & ! temperature in lowest model level and at surface
      & qv_atm, qc_atm, qi_atm    ! water vapor, cloud water and cloud ice in lowest model level
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & hus2m,                  & ! specific humidity at 2 meter
      & dew2m                     ! dewpoint temperature at 2 meter

    INTEGER :: jb, jls, js
    REAL(wp) :: &
      & qsat_atm, qsat_2m, & ! saturated humidity in lowest model level and at 2 meter
      & qrel_atm,          & ! relative humidity in lowest model level
      & pres2m               ! pressure at 2 meter

    REAL(wp), PARAMETER :: zephum = 0.05_wp ! epsilon for rel. humidity
    !$ACC DECLARE COPYIN(zephum)

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_2m_humidity'

!$OMP PARALLEL
    CALL init(hus2m, lacc=.TRUE.)
!$OMP END PARALLEL

    ! Note: Below it is assumed that relative humidity is constant with height. This should
    ! be revisited, see comments in https://gitlab.dkrz.de/icon/icon-mpim/-/merge_requests/341
    !
!$OMP PARALLEL DO PRIVATE(jb, jls, js, qsat_atm, qrel_atm, pres2m, qsat_2m) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
      !$ACC   PRIVATE(js, qsat_atm, qrel_atm, pres2m, qsat_2m)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        qsat_atm = specific_humidity(sat_pres_water(tatm(js,jb)),patm(js,jb))
        qrel_atm = MAX(zephum, qv_atm(js,jb) / qsat_atm)
        pres2m = psfc(js,jb) * &
          &  (1._wp - 2._wp * grav / ( rd * t2m(js,jb) * (1._wp + vtmpc1 * qv_atm(js,jb) - qc_atm(js,jb) - qi_atm(js,jb))))
        qsat_2m = specific_humidity(sat_pres_water(t2m(js,jb)),pres2m)
        hus2m(js,jb) = qrel_atm * qsat_2m
        dew2m(js,jb) = dewpoint_temperature(t2m(js,jb), hus2m(js,jb), pres2m)
      END DO
      !$ACC END PARALLEL LOOP
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_2m_humidity_and_dewpoint
  !
  !=================================================================
  !

END MODULE mo_tmx_surface_interface
