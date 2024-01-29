!
! mo_art_2mom_driver
! This routine is an adapted version of the two_moment_mcrph by A. Seifert
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

MODULE mo_art_2mom_driver
! ICON
  USE mo_kind,                   ONLY: wp
  USE mo_run_config,             ONLY: iqr, iqi, iqs, iqg, iqh, iqc,         &
    &                                  iqnr,iqni,iqns,iqng,iqnh,iqnc,ininact
  USE mo_exception,              ONLY: finish, message
  USE mo_io_units,               ONLY: find_next_free_unit
  USE mo_physical_constants,     ONLY: alv,           & ! latent heat of vaporization
    &                                  als,           & ! latent heat of sublimation
    &                                  cpdr  => rcpd, & ! (spec. heat of dry air at const press)^-1
    &                                  cvdr  => rcvd    ! (spec. heat of dry air at const vol)^-1
  USE mo_2mom_mcrph_processes,   ONLY: sedi_icon_rain,sedi_icon_sphere, & ! sedimentation routines
    &                                  cfg_params
  USE mo_2mom_mcrph_config_default, ONLY: cfg_2mom_default
  USE mo_2mom_mcrph_main,        ONLY: init_2mom_scheme,                          &
    &                                  init_2mom_scheme_once,                     &
    &                                  snow_coeffs, ice_coeffs, graupel_coeffs,   &
    &                                  hail_coeffs, rain_coeffs
  USE mo_2mom_mcrph_util, ONLY:        init_dmin_wg_gr_ltab_equi,                 &
     &                                 dmin_wetgrowth_fit_check,                  &
     &                                 luse_dmin_wetgrowth_table,                 &
     &                                 lprintout_comp_table_fit
  USE mo_2mom_mcrph_types,       ONLY: ltabdminwgg,                               &
     &                                 atmosphere, particle, particle_frozen
! ART
  USE mo_art_clipping,           ONLY: art_clip_lt, art_clip_gt
  USE mo_art_2mom_main,          ONLY: art_clouds_2mom
  USE mo_art_2mom_prepare,       ONLY: art_prepare_2mom, art_post_2mom
  
  IMPLICIT NONE
    
  PRIVATE
  
  PUBLIC :: art_2mom_mcrph
  PUBLIC :: art_2mom_mcrph_init
  
  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_2mom_driver'
  INTEGER,          PARAMETER :: dbg_level = 25                   ! level for debug prints
  
  ! .. exponents for simple height dependency of terminal fall velocity
  REAL(wp), PARAMETER :: rho_vel    = 0.4e0_wp    !..exponent for density correction
  REAL(wp), PARAMETER :: rho_vel_c  = 0.2e0_wp    !..for cloud droplets
  REAL(wp), PARAMETER :: rho0       = 1.225_wp    !..Norm-Luftdichte
  
  INTEGER, PARAMETER  :: art_cloud_type = 2603
  
  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_2mom_mcrph(isize, ke, jg, jb, is, ie, ks, dt,   &
                          & dz, rho, pres, tke, p_trac, tk,    &
                          & w, prec_r, prec_i, prec_s,         &
                          & prec_g, prec_h, tkvh, l_cv)
!<
! SUBROUTINE art_2mom_mcrph
! Driver routine for the 2 moment microphysics
! Based on: Seifert and Beheng (2006) - A two-moment cloud microphysics parameterization 
!                                       for mixed-phase clouds. Part 1: Model description
! Part of Module: mo_art_2mom_driver
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-10
! Modifications:
! 2015-12-09: Daniel Rieger, KIT
! - Adaptions to refactored 2MOM scheme
!>
  ! Setup variables (Grid, timestep, looping)
  INTEGER, INTENT (in)            :: &
    &  isize, ke,                    & !< grid sizes
    &  jg,jb,                        & !< domain index (p_patch%id), block index
    &  is, ie, ks                      !< start/end indices
  REAL(wp), INTENT(in)            :: &
    &  dt                              !< time step
  ! Dynamical core variables
  REAL(wp), INTENT(in), TARGET    :: &
    &  dz(:,:),                      & !< Vertical layer thickness
    &  rho(:,:),                     & !< Density
    &  pres(:,:),                    & !< Pressure
    &  tke(:,:),                     & !< Tubulent kinetic energy
    &  w(:,:),                       & !< Vertical velocity
    &  tkvh(:,:)                       !< Turbulent diffusion coefficient for heat
  REAL(wp), INTENT(inout), TARGET :: &
    &  tk(:,:)                         !< Temperature
  ! Tracer fields
  REAL(wp), INTENT(inout), TARGET :: &
    &  p_trac(:,:,:)                   !< Tracer fields
  ! Precip rates, vertical profiles
  REAL(wp), INTENT (inout)        :: &
    &  prec_r(:),                    & !< Precipitation rate for rain
    &  prec_i(:),                    & !< Precipitation rate for ice
    &  prec_s(:),                    & !< Precipitation rate for snow
    &  prec_g(:),                    & !< Precipitation rate for graupel
    &  prec_h(:)                       !< Precipitation rate for hail
  ! Switches
  LOGICAL, INTENT (in)            :: &
    &  l_cv                            !< Use c_v (true) or c_p (false)
  ! Local Variables
  INTEGER                         :: &
    &  ii,kk,                        & !< loop indizes
    &  ntsedi                          !< for sedimentation sub stepping
  REAL(wp)                        :: &
    &  rdz(isize,ke),                & !< 1/dz
    &  rho_r(isize,ke),              & !< 1/rho
    &  q_liq_old(isize,ke),          & !< to store old values for latent heat calc (liquid water)
    &  q_vap_old(isize,ke),          & !< to store old values for latent heat calc (vapor)
    &  q_liq_new,q_vap_new,          & !< new values for latent heat calc (liquid,vapor)
    &  z_heat_cap_r,                 & !< reciprocal of cpdr or cvdr (depending on l_cv)
    &  convliq,convice,              & !< latent heat term for temperature equation
    &  hlp,                          & !< help variable
    &  tau_inact =  600.               !< relaxation time scale for activated IN number density
  REAL(wp), TARGET                :: &
    &  rhocorr(isize,ke),            & !< density dependency of particle fall speed
    &  rhocld(isize,ke)                !< density dependency of particle fall speed 
                                       !    for cloud droplets
  LOGICAL, PARAMETER              :: &
    &  clipping  = .TRUE.,           & !< not really necessary, just for cleanup
    &  explicit_solver = .TRUE.        !< explicit or semi-implicit solver
  TYPE(atmosphere)                :: &
    &  atmo                            !< Atmospheric state
  ! These are the fundamental hydrometeor particle variables for the two-moment scheme
  TYPE(particle), target          :: &
    &  cloud_hyd, rain_hyd
  TYPE(particle_frozen), target   :: &
    &  ice_frz, snow_frz,            &
    &  graupel_frz, hail_frz
!  TYPE(particle_lwf), target      :: &
!    &  graupel_lwf, hail_lwf
  ! Pointers to the derived types that are actually needed
  CLASS(particle), pointer        :: &
    &  cloud, rain
  CLASS(particle_frozen), pointer :: &
    &  ice, snow, graupel, hail
  INTEGER                         :: &
    &  ik_slice(4)                     !< Loop indices,ICON: (/jc_start,jc_end,jk_start,jk_end/)
  
  cloud => cloud_hyd
  rain  => rain_hyd
  ice   => ice_frz
  snow  => snow_frz
  
  ! DRIEG: THIS NEEDS TO BE ADAPTED AS SOON AS AXELS CHANGES ARE BACK IN ICON-DEV
!    IF (PRESENT(qgl)) THEN
!       graupel => graupel_lwf   ! gscp=4,5,6
!       hail => hail_lwf
!    ELSE
       graupel => graupel_frz   ! gscp=7
       hail => hail_frz
!    END IF
  
  ! Start/End indices
  ik_slice(1) = is
  ik_slice(2) = ie
  ik_slice(3) = ks
  ik_slice(4) = ke
  
  ! inverse of vertical layer thickness
  rdz(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4)) = &
    & 1._wp / dz(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4))
  
  IF (clipping) THEN
    ! Make sure there are no negative values in hydrometeor mass mixing ratios
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqr),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqi),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqs),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqg),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqh),0.0_wp)
  END IF
  
  IF (l_cv) THEN
    z_heat_cap_r = cvdr
  ELSE
    z_heat_cap_r = cpdr
  ENDIF
  
  DO kk = ik_slice(3), ik_slice(4)
    DO ii = ik_slice(1), ik_slice(2)
      ! ... 1/rho is used quite often
      rho_r(ii,kk) = 1.0 / rho(ii,kk)
      ! ... height dependency of terminal fall velocities
      hlp = log(rho(ii,kk)/rho0)
      rhocorr(ii,kk) = exp(-rho_vel*hlp)
      rhocld(ii,kk)  = exp(-rho_vel_c*hlp)
    END DO
  END DO
  
  
  !---------------------------------------------------------------------------------
  !--- set the particle types, but no calculations
  !---------------------------------------------------------------------------------
  
  CALL init_2mom_scheme(cloud,rain,ice,snow,graupel,hail)
  
  !---------------------------------------------------------------------------------
  !--- convert to densities and set pointerns to two-moment module
  !--- (pointers are used to avoid passing everything explicitly by argument and
  !--- to avoid local allocates within the OpenMP-loop, and keep everything on stack)
  !---------------------------------------------------------------------------------
  
  CALL art_prepare_2mom(atmo, cloud, rain, ice, snow, graupel, hail, &
    &                   rho, rhocorr, rhocld, pres, w, tk, p_trac, ik_slice)
!  
!  
!  
  IF (explicit_solver) THEN
    ! save old variables for latent heat calculation
    q_vap_old(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4)) =   &
      &  atmo%qv(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4))
    q_liq_old(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4)) =   &
      &   cloud%q(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4)) &
      & + rain%q(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4))

    ! this subroutine calculates all the microphysical sources and sinks
    CALL art_clouds_2mom(ik_slice,dt,atmo,tke,tkvh,dz,cloud,rain,ice,snow, &
      &                  graupel,hail,p_trac,jg,jb)
    
    WHERE(cloud%q(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4)) == 0.0_wp)  &
      &   cloud%n(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4)) = 0.0_wp
    
    ! latent heat term for temperature equation
    convice = z_heat_cap_r * als
    convliq = z_heat_cap_r * (alv-als)
    DO kk = ik_slice(3), ik_slice(4)
      DO ii = ik_slice(1), ik_slice(2)
        ! new variables
        q_vap_new = atmo%qv(ii,kk)
        q_liq_new = cloud%q(ii,kk) + rain%q(ii,kk)
        
        ! update temperature
        atmo%T(ii,kk) = atmo%T(ii,kk) - convice * rho_r(ii,kk) * (q_vap_new - q_vap_old(ii,kk))  &
             &                        + convliq * rho_r(ii,kk) * (q_liq_new - q_liq_old(ii,kk))
      ENDDO
    ENDDO
    
    ! if we solve explicitly, then sedimentation is done here after microphysics
    CALL art_hydrom_sedim_explicit()
    
  ELSE
    CALL finish ('art_2mom_mcrph', 'Semi-implicit solver not yet implemented into ART')
  ENDIF ! explicit_solver
  
  ! check_clouds routine for debug
  ! CALL check_clouds.....
  
  ! .. convert back and nullify two-moment pointers
  CALL art_post_2mom(atmo, cloud, rain, ice, snow, graupel, hail, &
    &                rho_r, p_trac, ik_slice)
!
  IF (clipping) THEN
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4), iqr),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4), iqi),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4), iqs),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4), iqg),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4), iqh),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqnr),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqni),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqns),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqng),0.0_wp)
    CALL art_clip_lt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqnh),0.0_wp)
  ENDIF
  
  !..relaxation of activated IN number density to zero
  DO kk = ik_slice(3), ik_slice(4)
    DO ii = ik_slice(1), ik_slice(2)
      IF(p_trac(ii,kk,iqi) == 0._wp) THEN
        p_trac(ii,kk,ininact) = p_trac(ii,kk,ininact) - p_trac(ii,kk,ininact)/tau_inact*dt
      END IF
    ENDDO
  ENDDO
  
  WHERE(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqc) < 1.0e-12_wp) &
    &   p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqnc) = 0.0_wp
  CALL art_clip_gt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqr),0.02_wp)
  CALL art_clip_gt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqi),0.02_wp)
  CALL art_clip_gt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqs),0.02_wp)
  CALL art_clip_gt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqg),0.02_wp)
  CALL art_clip_gt(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqh),0.02_wp)
  
  RETURN
  !
  ! end of driver routine, but many details are below in the contains-part of this subroutine
  !
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
  SUBROUTINE art_hydrom_sedim_explicit()
!<
! SUBROUTINE art_hydrom_sedim_explicit
! Explicit calculation of hydrometeor sedimentation
! Based on: 
! Part of Module: mo_art_2mom_driver
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-10
! Modifications:
! 2015-12-09: Daniel Rieger, KIT
! - Adaptions to refactored 2MOM scheme
!>
    REAL(wp) :: cmax
    INTEGER :: ii2
    REAL(wp) :: precrate3D(ik_slice(2),ik_slice(4))

    prec_r(:) = 0._wp
    prec_i(:) = 0._wp
    prec_s(:) = 0._wp
    prec_g(:) = 0._wp
    prec_h(:) = 0._wp
    cmax      = 0.0_wp
  
    CALL sedi_icon_rain(rain,rain_coeffs,p_trac(:,:,iqr),p_trac(:,:,iqnr),prec_r,  &
      &                 precrate3D(:,:), p_trac(:,:,iqc), &
      &                 rhocorr,rdz,dt,ik_slice(1),ik_slice(2),ik_slice(3),ik_slice(4),cmax)
    
    
    IF (art_cloud_type >= 1000) THEN
       
      IF (ANY(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqi)>0._wp)) THEN
        CALL sedi_icon_sphere (ice, ice_coeffs, p_trac(:,:,iqi),p_trac(:,:,iqni),prec_i, &
          &                    precrate3D, &
          &                    rhocorr,rdz,dt,ik_slice(1),ik_slice(2),ik_slice(3),ik_slice(4))
      ENDIF
      IF (ANY(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqs)>0._wp)) THEN
        CALL sedi_icon_sphere (snow,snow_coeffs,p_trac(:,:,iqs),p_trac(:,:,iqns),prec_s, &
          &                    precrate3D, &
          &                    rhocorr,rdz,dt,ik_slice(1),ik_slice(2),ik_slice(3),ik_slice(4))
      ENDIF
      
      IF (ANY(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqg)>0._wp)) THEN
        ntsedi = 1
        DO ii2=1,ntsedi
          CALL sedi_icon_sphere (graupel,graupel_coeffs,p_trac(:,:,iqg),p_trac(:,:,iqng), &
            &                    prec_g,precrate3D, rhocorr,rdz,dt/ntsedi,                &
            &                    ik_slice(1),ik_slice(2),ik_slice(3),ik_slice(4),cmax)
        ENDDO
      ENDIF
    ENDIF
    
    IF (art_cloud_type >= 2000) THEN
      prec_h(:) = 0.0_wp
      IF (ANY(p_trac(ik_slice(1):ik_slice(2),ik_slice(3):ik_slice(4),iqh)>0._wp)) THEN
        ntsedi = 1
        DO ii2=1,ntsedi
          CALL sedi_icon_sphere (hail,hail_coeffs,p_trac(:,:,iqh),p_trac(:,:,iqnh), &
            &                    prec_h,precrate3D,rhocorr,rdz,dt/ntsedi,           &
            &                    ik_slice(1),ik_slice(2),ik_slice(3),ik_slice(4),cmax)
        END DO
      ENDIF
    ENDIF

  END SUBROUTINE art_hydrom_sedim_explicit
!!
!!-------------------------------------------------------------------------
!!
END SUBROUTINE art_2mom_mcrph
!!
!!-------------------------------------------------------------------------
!!

SUBROUTINE art_2mom_mcrph_init(msg_level)
!<
! SUBROUTINE art_2mom_mcrph_init
! Initialization routine for the 2 moment microphysics
! Based on: Seifert and Beheng (2006) - A two-moment cloud microphysics parameterization 
!                                       for mixed-phase clouds. Part 1: Model description
! Part of Module: mo_art_2mom_driver
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-10
! Modifications:
! 2015-12-09: Daniel Rieger, KIT
! - Adaptions to refactored 2MOM scheme
!>
  INTEGER, INTENT(IN)           :: &
    &  msg_level                   !< message level
  ! Local variables
  TYPE(particle)                :: &
    &  cloud, rain
  TYPE(particle_frozen)         :: &
    &  ice, snow, graupel, hail
  INTEGER                       :: &
    &  unitnr                       !< input unit number
  
  ! Transfer the configuration parameters to the 2mom internal type instance:
  cfg_params = cfg_2mom_default

  ! .. set the particle types, and calculate some coefficients
  CALL init_2mom_scheme_once(cloud,rain,ice,snow,graupel,hail,art_cloud_type)

  IF (luse_dmin_wetgrowth_table .OR. lprintout_comp_table_fit) THEN
    unitnr = find_next_free_unit(10,30)
    CALL init_dmin_wg_gr_ltab_equi('dmin_wetgrowth_lookup', graupel, &
         unitnr, 61, ltabdminwgg, msg_level)
  END IF
  IF (.NOT. luse_dmin_wetgrowth_table) THEN
    ! check whether 4d-fit is consistent with graupel parameters
    IF (dmin_wetgrowth_fit_check(graupel)) THEN 
      CALL message (TRIM(routine), " Using 4d-fit for dmin_wetgrowth for "//TRIM(graupel%name))
    ELSE
      CALL finish(TRIM(routine),&
           & 'Error: luse_dmin_wetgrowth_table=.false., so 4D-fit should be used, '// &
           & 'but graupel parameters inconsistent with 4d-fit') 
    END IF
  END IF

  IF (msg_level>dbg_level) CALL message (TRIM(routine), " finished init_dmin_wetgrowth for "//TRIM(graupel%name))
  
END SUBROUTINE art_2mom_mcrph_init
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_2mom_driver
