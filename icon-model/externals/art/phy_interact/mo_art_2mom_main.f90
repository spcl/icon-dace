!
! mo_art_2mom_main
! This routine serves as wrapper for the process routines within
! mo_2mom_mcrph_main by A. Seifert
!
! NOTE/DISCLAIMER: The ART version of the 2MOM scheme does not grant
! exact reproducability of ICON 2MOM scheme results. There are differences
! due to clipping etc.
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

MODULE mo_art_2mom_main
! ICON
  USE mo_kind,                          ONLY: wp
! Microphysical Process Routines
  USE mo_2mom_mcrph_processes,          ONLY: cloud_freeze, vapor_dep_relaxation,             &
    &                                         ice_selfcollection, snow_selfcollection,        &
    &                                         graupel_selfcollection,                         &
    &                                         particle_particle_collection,                   &
    &                                         graupel_hail_conv_wet_gamlook,                  &
    &                                         ice_riming, snow_riming,rain_freeze_gamlook,    &
    &                                         ice_melting,snow_melting,graupel_melting,       &
    &                                         prepare_melting_lwf, particle_melting_lwf,      &
    &                                         hail_melting_simple, evaporation,               &
    &                                         particle_cloud_riming, particle_rain_riming,    &
    &                                         autoconversionKB, accretionKB,                  &
    &                                         rain_selfcollectionSB,autoconversionKK,         &
    &                                         accretionKK, autoconversionSB, accretionSB,     &
    &                                         rain_evaporation, ice_nucleation_homhet
  ! And some parameters declared in the process module that are required here. 
  USE mo_2mom_mcrph_processes,          ONLY: ice_typ, auto_typ
  USE mo_2mom_mcrph_main,               ONLY: qnc_const, rain_coeffs, snow_coeffs,         &
    &                                         sic_coeffs, gic_coeffs, gsc_coeffs,          &
    &                                         hic_coeffs, hsc_coeffs, hail_coeffs,         &
    &                                         cloud_coeffs, ice_coeffs, graupel_coeffs
  USE mo_2mom_mcrph_main,               ONLY: scr_coeffs, srr_coeffs, irr_coeffs,          &
    &                                         icr_coeffs, hrr_coeffs, grr_coeffs,          &
    &                                         hcr_coeffs, gcr_coeffs,                      &
    &                                         graupel_ltable1, graupel_ltable2,            &
    &                                         rain_ltable1, rain_ltable2, rain_ltable3,    &
    &                                         rain_nm1, rain_nm2, rain_nm3, rain_g1,       &
    &                                         rain_g2, graupel_nm1, graupel_nm2,           &
    &                                         graupel_g1, graupel_g2, rain_gfak
  USE mo_2mom_mcrph_types,              ONLY: particle, particle_frozen, particle_lwf, atmosphere, &
    &                                         ltabdminwgg, ltab_estick_ice, ltab_estick_snow, ltab_estick_parti
  USE mo_run_config,                    ONLY: ininact
! ART
  USE mo_art_nucleation_interface,      ONLY: art_nuc_warm_interface,art_nuc_cold_interface
  USE mo_art_config,                    ONLY: art_config
  
  IMPLICIT NONE
    
  PRIVATE
  
  PUBLIC :: art_clouds_2mom
  
  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clouds_2mom(ik_slice,dt,atmo,tke,tkvh,dz,cloud,rain,ice,snow,graupel,hail, &
  &                        p_trac,jg,jb)
!<
! SUBROUTINE art_clouds_2mom
! Main subroutine of the ART two-moment microphysics including aerosol-cloud-interactions.
! All individual processes are called in sequence and do their own time 
! integration for the mass and number densities, i.e., Marchuk-type operator 
! splitting. Temperature is only updated in the driver.
! Based on: Seifert and Beheng (2006) - A two-moment cloud microphysics parameterization 
!                                       for mixed-phase clouds. Part 1: Model description
! Part of Module: mo_art_2mom_main
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-10
! Modifications:
! 2015-12-09: Daniel Rieger, KIT
! - Adaptions to refactored 2MOM scheme
!>
  INTEGER,INTENT(in)          :: &
    &  ik_slice(4)                 !< Loop indices,ICON: (/jc_start,jc_end,jk_start,jk_end/)
  REAL(wp), INTENT(in)        :: &
    &  dt                          !< time step
  REAL(wp), INTENT(in)        :: &
    &  tke(:,:),                 & !< Tubulent kinetic energy
    &  tkvh(:,:),                & !< Turbulent diffusion coefficient for heat
    &  dz(:,:)                     !< Layer height
  TYPE(atmosphere),INTENT(inout) :: &
    &  atmo                        !< Atmospheric state
  CLASS(particle), INTENT(inout) :: &
    &  cloud, rain
  CLASS(particle_frozen), INTENT(inout) :: &
    &  ice, snow, graupel, hail
  REAL(wp), INTENT(inout), TARGET :: &
    &  p_trac(:,:,:)               !< Tracer fields
  INTEGER, INTENT (IN)        :: &
    &  jg, jb                      !< Domain and Block index
  ! Local variables
  INTEGER                     :: &
    &  k, i                        !< Loop indices
  REAL(wp)                    :: &
    &  dep_rate_ice(size(cloud%n,1),size(cloud%n,2)),  & !< Deposition rate for ice
    &  dep_rate_snow(size(cloud%n,1),size(cloud%n,2)), & !< Deposition rate for snow
    &  gmelting(size(cloud%n,1),size(cloud%n,2))         !< ambient atmospheric conditions for
                                                         !    melting in lwf scheme
  
  ! ----------------------------------
  ! --- Initializations
  ! ----------------------------------
  
  dep_rate_ice(:,:)  = 0.0_wp
  dep_rate_snow(:,:) = 0.0_wp
  
  ! ----------------------------------
  ! --- Calculate Microphysical Processes
  ! ----------------------------------
  
  ! Nucleation Warm
  SELECT CASE(art_config(jg)%iart_aci_warm)
    CASE(0)
      ! Constant cloud droplet number inwp_gscp=4
      cloud%n(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4)) = qnc_const
    CASE(1)
      ! Nucleation Warm ART
      CALL art_nuc_warm_interface(atmo%w,atmo%T,atmo%p,atmo%rho,tke,tkvh,dz,p_trac,      &
        &                         jg,jb,ik_slice(1),ik_slice(2),ik_slice(3),ik_slice(4), &
        &                         dt,cloud%x_min,art_config(jg)%lart_diag_out)
  END SELECT
  
  IF (art_config(jg)%iart_aci_warm > 0) THEN
    DO k=ik_slice(3),ik_slice(4)
      DO i=ik_slice(1),ik_slice(2)
        cloud%n(i,k) = MAX(cloud%n(i,k), cloud%q(i,k) / cloud%x_max)
        cloud%n(i,k) = MIN(cloud%n(i,k), cloud%q(i,k) / cloud%x_min)
      END DO
    END DO
  ENDIF
  
  ! Nucleation Cold
  SELECT CASE(art_config(jg)%iart_aci_cold)
    CASE(0)
      ! Nucleation originial 2MOM (KL06)
      CALL ice_nucleation_homhet(ik_slice,.false.,atmo,cloud,ice,n_inact = p_trac(:,:,ininact))
    CASE(1,2,3,4,5)
      ! Nucleation Cold ART
      CALL art_nuc_cold_interface(atmo%w,atmo%T,atmo%p,atmo%rho,tke,p_trac,ice%x_min,        &
        &                         ik_slice(1),ik_slice(2),ik_slice(3),ik_slice(4),jg,jb,     &
        &                         art_config(jg)%iart_aci_cold,art_config(jg)%lart_diag_out)
    CASE(6)
      ! BN09 with tracking of activated dust -> Use PDA13 (=6) het. nuc. scheme
      CALL art_nuc_cold_interface(atmo%w,atmo%T,atmo%p,atmo%rho,tke,p_trac,ice%x_min,        &
        &                         ik_slice(1),ik_slice(2),ik_slice(3),ik_slice(4),jg,jb,     &
        &                         art_config(jg)%iart_aci_cold,art_config(jg)%lart_diag_out, &
        &                         p_trac(:,:,ininact))
    CASE(7) ! TODO
!     ! KL06 with prognostic dust as input and relaxation of activated dust IN
!    CALL ice_nucleation_homhet(ik_slice,.true.,atmo,cloud,ice, &
!      &                        n_inact = p_trac(:,:,ininact) , &
!      &                        n_inpot = p_art_data(jg)%diag%ndust_tot(:,:,jb))
  END SELECT
  
  ! homogeneous freezing of cloud droplets
  CALL cloud_freeze(ik_slice, dt, cloud_coeffs, qnc_const, atmo, cloud, ice)
       
  DO k = ik_slice(3), ik_slice(4)
    DO i = ik_slice(1), ik_slice(2)
      ice%n(i,k) = MIN(ice%n(i,k), ice%q(i,k)/ice%x_min)
      ice%n(i,k) = MAX(ice%n(i,k), ice%q(i,k)/ice%x_max)
    END DO
  END DO
  
  ! depositional growth of all ice particles
  ! ( store deposition rate of ice and snow for conversion calculation in 
  !   ice_riming and snow_riming )
  CALL vapor_dep_relaxation(ik_slice,dt,ice_coeffs,snow_coeffs,graupel_coeffs,hail_coeffs,&
         &                    atmo,ice,snow,graupel,hail,dep_rate_ice,dep_rate_snow)
  
  ! ice-ice collisions
  CALL ice_selfcollection(ik_slice,dt,atmo,ice,snow,ice_coeffs,ltab_estick_ice)
  CALL snow_selfcollection(ik_slice,dt,atmo,snow,snow_coeffs,ltab_estick_snow)
  CALL particle_particle_collection(ik_slice, dt, atmo, ice, snow, sic_coeffs, ltab_estick_parti)
  
  CALL graupel_selfcollection(ik_slice, dt, atmo, graupel, graupel_coeffs)
  CALL particle_particle_collection(ik_slice, dt, atmo, ice, graupel, gic_coeffs, ltab_estick_parti)
  CALL particle_particle_collection(ik_slice, dt, atmo, snow, graupel, gsc_coeffs, ltab_estick_parti)
  
  IF (ice_typ > 1) THEN
    ! conversion of graupel to hail in wet growth regime
    CALL graupel_hail_conv_wet_gamlook(ik_slice, graupel_ltable1, graupel_ltable2,       &
                                       graupel_nm1, graupel_nm2, graupel_g1, graupel_g2, &
                                       ltabdminwgg, atmo, graupel, cloud, rain, ice, snow, hail)
    ! hail collisions
    CALL particle_particle_collection(ik_slice, dt, atmo, ice, hail, hic_coeffs, ltab_estick_parti)    ! Important?
    CALL particle_particle_collection(ik_slice, dt, atmo, snow, hail, hsc_coeffs, ltab_estick_parti)
  END IF
  
  ! riming of ice with cloud droplets and rain drops, and conversion to graupel
  CALL ice_riming(ik_slice, dt, icr_coeffs, irr_coeffs, atmo, ice, cloud,       &
    &             rain, graupel, dep_rate_ice)

  ! riming of snow with cloud droplets and rain drops, and conversion to graupel
  CALL snow_riming(ik_slice, dt, scr_coeffs, srr_coeffs, atmo, snow, cloud,     &
    &              rain, ice, graupel, dep_rate_snow)
  
  ! more riming
  IF (ice_typ > 1) THEN
     CALL particle_cloud_riming(ik_slice, dt, atmo, hail, hcr_coeffs, cloud, rain, ice)
     CALL particle_rain_riming(ik_slice, dt, atmo, hail, hrr_coeffs, rain, ice)
  END IF
  CALL particle_cloud_riming(ik_slice, dt, atmo, graupel, gcr_coeffs, cloud, rain, ice)
  CALL particle_rain_riming(ik_slice, dt, atmo, graupel, grr_coeffs, rain, ice)
  
  ! freezing of rain and conversion to ice/graupel/hail
  CALL rain_freeze_gamlook(ik_slice, dt, rain_ltable1, rain_ltable2, rain_ltable3, &
    &                      rain_nm1, rain_nm2, rain_nm3, rain_g1, rain_g2,         &
    &                      rain_coeffs,atmo,rain,ice,snow,graupel,hail)
  
  ! melting
  CALL ice_melting(ik_slice, atmo, ice, cloud, rain)
  CALL snow_melting(ik_slice,dt,snow_coeffs,atmo,snow,rain)
  ! melting of graupel and hail can be simple or LWF-based
  SELECT TYPE (graupel)
    TYPE is (particle_frozen) 
      CALL graupel_melting(ik_slice,dt,graupel_coeffs,atmo,graupel,rain)
    TYPE is (particle_lwf)
      CALL prepare_melting_lwf(ik_slice, atmo, gmelting)
      CALL particle_melting_lwf(ik_slice, dt, graupel, rain, gmelting)
    END SELECT
  SELECT TYPE (hail)
    TYPE is (particle_frozen) 
      CALL hail_melting_simple(ik_slice,dt,hail_coeffs,atmo,hail,rain)
    TYPE is (particle_lwf)
      CALL particle_melting_lwf(ik_slice, dt, hail, rain, gmelting)
  END SELECT
  
  ! evaporation from melting ice particles
  CALL evaporation(ik_slice, dt, atmo, snow, snow_coeffs)
  CALL evaporation(ik_slice, dt, atmo, graupel, graupel_coeffs)
  CALL evaporation(ik_slice, dt, atmo, hail, hail_coeffs)


  ! warm rain processes 
  ! (using something other than SB is somewhat inconsistent and not recommended)
  IF (auto_typ == 1) THEN
    CALL autoconversionKB(ik_slice, dt, cloud, rain)   ! Beheng (1994)
    CALL accretionKB(ik_slice, dt, cloud, rain)
    CALL rain_selfcollectionSB(ik_slice, dt, atmo, rain)
  ELSE IF (auto_typ == 2) THEN
    ! Khairoutdinov and Kogan (2000)
    ! (KK2000 originally assume a 25 micron size threshold)
    CALL autoconversionKK(ik_slice, dt, cloud, rain)
    CALL accretionKK(ik_slice, dt, cloud, rain)
    CALL rain_selfcollectionSB(ik_slice, dt, atmo, rain)
  ELSE IF (auto_typ == 3) THEN
    CALL autoconversionSB(ik_slice, dt, atmo, cloud_coeffs, cloud, rain)   ! Seifert and Beheng (2001)
    CALL accretionSB(ik_slice, dt, atmo, cloud, rain)
    CALL rain_selfcollectionSB(ik_slice, dt, atmo, rain)
  ENDIF
  
  ! evaporation of rain following Seifert (2008)
  CALL rain_evaporation(ik_slice, dt, rain_coeffs, rain_gfak, atmo, cloud, rain)
  
  ! ----------------------------------
  ! --- Size Limits for all Hydrometeors
  ! ----------------------------------
  IF (art_config(jg)%iart_aci_warm > 0) THEN
    DO k = ik_slice(3), ik_slice(4)
      DO i = ik_slice(1), ik_slice(2)
        cloud%n(i,k) = MIN(cloud%n(i,k), cloud%q(i,k)/cloud%x_min)
        cloud%n(i,k) = MAX(cloud%n(i,k), cloud%q(i,k)/cloud%x_max)
        ! Hard upper limit for cloud number conc.
        cloud%n(i,k) = MIN(cloud%n(i,k), 5000d6)
      ENDDO
    ENDDO
  ENDIF
  DO k = ik_slice(3), ik_slice(4)
    DO i = ik_slice(1), ik_slice(2)
      rain%n(i,k) = MIN(rain%n(i,k), rain%q(i,k)/rain%x_min)
      rain%n(i,k) = MAX(rain%n(i,k), rain%q(i,k)/rain%x_max)
    ENDDO
  ENDDO
  IF (ice_typ > 0) THEN
    DO k = ik_slice(3), ik_slice(4)
      DO i = ik_slice(1), ik_slice(2)
        ice%n(i,k)     = MIN(ice%n(i,k), ice%q(i,k)/ice%x_min)
        ice%n(i,k)     = MAX(ice%n(i,k), ice%q(i,k)/ice%x_max)
        snow%n(i,k)    = MIN(snow%n(i,k), snow%q(i,k)/snow%x_min)
        snow%n(i,k)    = MAX(snow%n(i,k), snow%q(i,k)/snow%x_max)
        graupel%n(i,k) = MIN(graupel%n(i,k), graupel%q(i,k)/graupel%x_min)
        graupel%n(i,k) = MAX(graupel%n(i,k), graupel%q(i,k)/graupel%x_max)
      ENDDO
    ENDDO
  ENDIF
  IF (ice_typ > 1) THEN
    DO k = ik_slice(3), ik_slice(4)
      DO i = ik_slice(1), ik_slice(2)
        hail%n(i,k) = MIN(hail%n(i,k), hail%q(i,k)/hail%x_min)
        hail%n(i,k) = MAX(hail%n(i,k), hail%q(i,k)/hail%x_max)
      ENDDO
    ENDDO
  ENDIF
  
END SUBROUTINE art_clouds_2mom
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_2mom_main
