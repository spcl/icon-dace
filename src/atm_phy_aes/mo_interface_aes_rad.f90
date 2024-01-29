!
! Subroutine interface_aes_rad calls the radiative transfer scheme.
!
!
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
#include "omp_definitions.inc"
!----------------------------

MODULE mo_interface_aes_rad

  USE mo_kind,                 ONLY: wp
  USE mtime,                   ONLY: t_datetime => datetime, OPERATOR(<=), OPERATOR(>)

  USE mo_aes_phy_dims,         ONLY: aes_phy_dims
  USE mo_aes_phy_config,       ONLY: aes_phy_tc
  USE mo_aes_phy_memory,       ONLY: t_aes_phy_field, prm_field
  USE mo_snow_ice_reff,        ONLY: ice_reff_moss, snow_reff_funeedles, ice_reff_funeedles, snow_x, ice_x
  USE mo_run_config,           ONLY: iqs, iqi
  USE mo_aes_graupel,          ONLY: snow_number, snow_lambda, ice_number

#ifndef __NO_RTE_RRTMGP__
  USE mo_rte_rrtmgp_radiation, ONLY: rte_rrtmgp_radiation
#endif

  USE mo_timer,                ONLY: ltimer, timer_start, timer_stop, timer_rad

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_aes_rad

CONTAINS

  SUBROUTINE interface_aes_rad(jg, jb, jcs, jce)

    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    ! Pointers
    !
    TYPE(t_aes_phy_field), POINTER :: field

    ! Local variables
    !
    INTEGER  :: ntracer
    INTEGER  :: nlev
    INTEGER  :: nproma
    !
    TYPE(t_datetime), POINTER :: datetime
    LOGICAL  :: is_in_sd_ed_interval
    LOGICAL  :: is_active
    !
    INTEGER :: i1, i2, i3
    LOGICAL :: loland(aes_phy_dims(jg)%nproma)
    LOGICAL :: loglac(aes_phy_dims(jg)%nproma)
    !
    ! temp variable since non-contiguous slicing is not supported in OpenACC
    REAL(wp):: qtrc_phy(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev,aes_phy_dims(jg)%ntracer)

    IF (ltimer) CALL timer_start(timer_rad)

    ntracer = aes_phy_dims(jg)%ntracer
    nlev    = aes_phy_dims(jg)%nlev
    nproma  = aes_phy_dims(jg)%nproma

    datetime             => aes_phy_tc(jg)%datetime
    is_in_sd_ed_interval =  aes_phy_tc(jg)%is_in_sd_ed_interval_rad
    is_active            =  aes_phy_tc(jg)%is_active_rad

    ! associate pointers
    field => prm_field(jg)

    IF ( is_in_sd_ed_interval ) THEN
        !
#ifndef __NO_RTE_RRTMGP__
        IF ( is_active ) THEN
          !
          ! store ts_rad of this radiatiative transfer timestep in ts_rad_rt,
          ! so that it can be reused in radheat in the other timesteps
          !
          !DA TODO: remove derived % after merging with mc8
          !$ACC DATA PRESENT(field, field%ts_rad_rt, field%ts_rad, field%sftlf, field%sftgif) &
          !$ACC   PRESENT(field%qtrc_phy) & ! ACCWA (nvhpc on levante): to prevent illegal address during kernel execution
          !$ACC   CREATE(loland, loglac, qtrc_phy)
          !
          ! For some reason ACC kernels result in SEQ loop here
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO i2=1,SIZE(field%ts_rad_rt,2)
            DO i1=1,SIZE(field%ts_rad_rt,1)
              field%ts_rad_rt(i1,i2) = field%ts_rad(i1,i2)
            END DO
          END DO
          !$ACC END PARALLEL LOOP
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(3) ASYNC(1)
          DO i3=1,SIZE(qtrc_phy,3)
            DO i2=1,SIZE(qtrc_phy,2)
              DO i1=1,SIZE(qtrc_phy,1)
                qtrc_phy(i1,i2,i3) = field%qtrc_phy(i1,i2,jb,i3)
              END DO
            END DO
          END DO
          !$ACC END PARALLEL LOOP
          !
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
          DO i1=1,SIZE(loland,1)
            loland(i1) = field% sftlf (i1,jb) > 0._wp
            loglac(i1) = field% sftgif(i1,jb) > 0._wp
          END DO
          !$ACC END PARALLEL LOOP
          !
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1)
          DO i2=1,nlev
            DO i1=jcs,jce
               !calculate reff of snow here (further down, the variables are per m^2 instead of m^3)
               field%acsnc(i1,i2,jb) = snow_number(field%ta(i1,i2,jb), field%rho(i1,i2,jb), field%qtrc_phy(i1,i2,jb,iqs))
               field%acsnc(i1,i2,jb) = field%acsnc(i1,i2,jb) / snow_lambda(field%rho(i1,i2,jb), field%qtrc_phy(i1,i2,jb,iqs), field%acsnc(i1,i2,jb))
               field%x_snow(i1,i2,jb) = snow_x(field%ta(i1,i2,jb), field%rho(i1,i2,jb),field%qtrc_phy(i1,i2,jb,iqs))
               field%reff_snow(i1,i2,jb) = snow_reff_funeedles(field%ta(i1,i2,jb), field%rho(i1,i2,jb), field%qtrc_phy(i1,i2,jb,iqs))

               !calculate reff of ice here (further down, the variables are per m^2 instead of m^3)
               !diagnose cloud ice number concentration first
               field%acinc(i1,i2,jb) = ice_number(field%ta(i1,i2,jb), field%rho(i1,i2,jb)) * field%rho(i1,i2,jb)
               field%x_ice(i1,i2,jb) = ice_x(field%ta(i1,i2,jb), field%rho(i1,i2,jb), field%qtrc_phy(i1,i2,jb,iqi))
               field%reff_ice(i1,i2,jb) = ice_reff_moss(field%ta(i1,i2,jb), field%rho(i1,i2,jb), field%qtrc_phy(i1,i2,jb,iqi))
            END DO
          END DO
          !$ACC END PARALLEL LOOP
          !
          CALL rte_rrtmgp_radiation(                        &
              jg, jb, jcs, jce, nproma, nlev, ntracer,      &
              & loland         = loland                    ,&!< in  land-sea mask. (logical)
              & loglac         = loglac                    ,&!< in  glacier mask (logical)
              & this_datetime  = datetime                  ,&!< in  actual time step
              & pcos_mu0       = field%cosmu0_rt(:,jb)     ,&!< in  solar zenith angle
              & daylght_frc    = field%daylght_frc_rt(:,jb),&!in daylight fraction
              & alb_vis_dir    = field%albvisdir(:,jb)     ,&!< in  surface albedo for visible range, direct
              & alb_nir_dir    = field%albnirdir(:,jb)     ,&!< in  surface albedo for near IR range, direct
              & alb_vis_dif    = field%albvisdif(:,jb)     ,&!< in  surface albedo for visible range, diffuse
              & alb_nir_dif    = field%albnirdif(:,jb)     ,&!< in  surface albedo for near IR range, diffuse
              & emissivity     = field%emissivity(:,jb)    ,&!< in  surface longwave emissivity
              & tk_sfc         = field%ts_rad_rt(:,jb)     ,&!< in  grid box mean surface temperature
              & zf             = field%zf(:,:,jb)          ,&!< in  geometric height at full level      [m]
              & zh             = field%zh(:,:,jb)          ,&!< in  geometric height at half level      [m]
              & dz             = field%dz(:,:,jb)          ,&!< in  geometric height thickness of layer [m]
              & pp_hl          = field%phalf(:,:,jb)       ,&!< in  pressure at half levels at t-dt [Pa]
              & pp_fl          = field%pfull(:,:,jb)       ,&!< in  pressure at full levels at t-dt [Pa]
              & tk_fl          = field%ta(:,:,jb)          ,&!< in  tk_fl  = temperature at full level at t-dt
              & xm_air         = field%mair(:,:,jb)        ,&!< in  air mass in layer [kg/m2]
              & xq_trc         = qtrc_phy                  ,&!< in  tracer  mass fraction [kg/kg]
              & xv_ozn         = field%o3(:,:,jb)          ,&!< out ozone volume mixing ratio [mol/mol]
              !
              & reff_ice       = field%reff_ice(:,:,jb)    ,&!< in   effective radius of cloud ice [m]
              & tau_ice        = field% tau_ice(:,:,jb)    ,&!< inout cloud ice optical depth integrated over all bands
              & reff_snow      = field%reff_snow(:,:,jb)   ,&!< in   effective radius of snow [m]
              & tau_snow       = field% tau_snow(:,:,jb)   ,&!< inout snow optical depth integrated over all bands
              & cdnc           = field% acdnc(:,:,jb)      ,&!< in   cloud droplet number conc
              !
              & lw_dnw_clr     = field%rldcs_rt(:,:,jb)    ,&!< out  Clear-sky net longwave  at all levels
              & lw_upw_clr     = field%rlucs_rt(:,:,jb)    ,&!< out  Clear-sky net longwave  at all levels
              & sw_dnw_clr     = field%rsdcs_rt(:,:,jb)    ,&!< out  Clear-sky net shortwave at all levels
              & sw_upw_clr     = field%rsucs_rt(:,:,jb)    ,&!< out  Clear-sky net shortwave at all levels
              & lw_dnw         = field%rld_rt  (:,:,jb)    ,&!< out  All-sky net longwave  at all levels
              & lw_upw         = field%rlu_rt  (:,:,jb)    ,&!< out  All-sky net longwave  at all levels
              & sw_dnw         = field%rsd_rt  (:,:,jb)    ,&!< out  All-sky net longwave  at all levels
              & sw_upw         = field%rsu_rt  (:,:,jb)    ,&!< out  All-sky net longwave  at all levels
              !
              & vis_dn_dir_sfc = field%rvds_dir_rt(:,jb)   ,&!< out  all-sky downward direct visible radiation at surface
              & par_dn_dir_sfc = field%rpds_dir_rt(:,jb)   ,&!< out  all-sky downward direct PAR     radiation at surface
              & nir_dn_dir_sfc = field%rnds_dir_rt(:,jb)   ,&!< out  all-sky downward direct near-IR radiation at surface
              & vis_dn_dff_sfc = field%rvds_dif_rt(:,jb)   ,&!< out  all-sky downward diffuse visible radiation at surface
              & par_dn_dff_sfc = field%rpds_dif_rt(:,jb)   ,&!< out  all-sky downward diffuse PAR     radiation at surface
              & nir_dn_dff_sfc = field%rnds_dif_rt(:,jb)   ,&!< out  all-sky downward diffuse near-IR radiation at surface
              & vis_up_sfc     = field%rvus_rt    (:,jb)   ,&!< out  all-sky upward visible radiation at surface
              & par_up_sfc     = field%rpus_rt    (:,jb)   ,&!< out  all-sky upward PAR     radiation at surfac
              & nir_up_sfc     = field%rnus_rt    (:,jb)   ,&!< out  all-sky upward near-IR radiation at surface
              & aer_aod_533    = field%aer_aod_533 (:,:,jb),&!< out  aerosol optical density at 533 nm
              & aer_ssa_533    = field%aer_ssa_533 (:,:,jb),&!< out  single scattering albedo at 533 nm
              & aer_asy_533    = field%aer_asy_533 (:,:,jb),&!< out  asymmetrie factor at 533 nm
              & aer_aod_2325   = field%aer_aod_2325(:,:,jb),&!< out  aerosol optical density at 2325 nm
              & aer_ssa_2325   = field%aer_ssa_2325(:,:,jb),&!< out  single scattering albedo at 2325 nm
              & aer_asy_2325   = field%aer_asy_2325(:,:,jb),&!< out  asymmetrie factor at 2325 nm
              & aer_aod_9731   = field%aer_aod_9731(:,:,jb)&!< out  aerosol optical density at 9731
              & )
              !
          !$ACC WAIT(1)
          !$ACC END DATA
        !
        END IF
#endif
          !
       ELSE
          !
          ! LW
          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
          field%rldcs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net longwave  at all levels
          field%rlucs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net longwave  at all levels
          field%rld_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          field%rlu_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          !
          ! SW all
          field%rsdcs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net shortwave at all levels
          field%rsucs_rt(:,:,:)  = 0.0_wp !< out  Clear-sky net shortwave at all levels
          field%rsd_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          field%rsu_rt  (:,:,:)  = 0.0_wp !< out  All-sky net longwave  at all levels
          !
          ! SW vis, par and nir
          field%rvds_dir_rt(:,:) = 0.0_wp !< out  all-sky downward direct visible radiation at surface
          field%rpds_dir_rt(:,:) = 0.0_wp !< all-sky downward direct PAR     radiation at surface
          field%rnds_dir_rt(:,:) = 0.0_wp !< all-sky downward direct near-IR radiation at surface
          field%rvds_dif_rt(:,:) = 0.0_wp !< all-sky downward diffuse visible radiation at surface
          field%rpds_dif_rt(:,:) = 0.0_wp !< all-sky downward diffuse PAR     radiation at surface
          field%rnds_dif_rt(:,:) = 0.0_wp !< all-sky downward diffuse near-IR radiation at surface
          field%rvus_rt    (:,:) = 0.0_wp !< all-sky upward visible radiation at surface
          field%rpus_rt    (:,:) = 0.0_wp !< all-sky upward PAR     radiation at surfac
          field%rnus_rt    (:,:) = 0.0_wp !< all-sky upward near-IR radiation at surface
          !
          ! total cloud cover diagnostics
          field%aclcov(:,:)      = 0.0_wp !< out  total cloud cover
          ! cloud ice and snow optical depth integrated over all bands
          field%tau_ice(:,:,:) = 0.0_wp
          field%tau_snow(:,:,:) = 0.0_wp
          !$ACC END KERNELS
          !
       !
       END IF

     IF (ltimer) CALL timer_stop(timer_rad)

  END SUBROUTINE interface_aes_rad

END MODULE mo_interface_aes_rad
