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

! Two-moment mixed-phase bulk microphysics

!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"

MODULE mo_2mom_mcrph_driver

!------------------------------------------------------------------------------
!
! Description:
!
!   The subroutines in the module "gscp" calculate the rates of change of
!   temperature, cloud condensate and water vapor due to cloud microphysical
!   processes related to the formation of grid scale clouds and precipitation.
!   In the COSMO model the microphysical subroutines are either
!   called from "organize_gscp" or from "organize_physics" itself.
!
! - uses intrinsic gamma function. In gfortran you may want to compile with
!   the option -fall-intrinsics
!
!==============================================================================
!
! Declarations:
!
! Modules used:
!------------------------------------------------------------------------------
! Microphysical constants and variables
!------------------------------------------------------------------------------

USE mo_kind,                 ONLY: wp
USE mo_math_constants,       ONLY: pi
USE mo_physical_constants,   ONLY: &
    rhoh2o,           & ! density of liquid water
    alv,              & ! latent heat of vaporization
    als,              & ! latent heat of sublimation
    cpdr  => rcpd,    & ! (spec. heat of dry air at constant press)^-1
    cvdr  => rcvd,    & ! (spec. heat of dry air at const vol)^-1
    rho_ice => rhoice   ! density of pure ice

USE mo_satad,                ONLY: latent_heat_sublimation, latent_heat_melting

USE mo_exception,            ONLY: finish, message, message_text
USE mo_run_config,           ONLY: ldass_lhn

USE mo_timer,                ONLY:                                                              & 
                              timers_level, timer_start, timer_stop, timer_phys_2mom_dmin_init, &
                              timer_phys_2mom_prepost, timer_phys_2mom_proc, timer_phys_2mom_sedi


USE mo_reff_types,           ONLY: t_reff_calc

USE mo_2mom_mcrph_config,    ONLY: t_cfg_2mom

USE mo_2mom_mcrph_main,      ONLY:                                &
     &                        clouds_twomoment,                   &
     &                        atmosphere, particle, particle_frozen, particle_lwf, &
     &                        rain_coeffs, ice_coeffs, snow_coeffs, graupel_coeffs, hail_coeffs, &
     &                        ccn_coeffs, in_coeffs,                     &
     &                        init_2mom_scheme, init_2mom_scheme_once,   &
     &                        qnc_const

USE mo_2mom_mcrph_processes,  ONLY:                                &
     &                         sedi_vel_rain, sedi_vel_sphere,     &
     &                         sedi_icon_rain, sedi_icon_sphere, sedi_icon_sphere_lwf, &
     &                         particle_meanmass, sedi_vel_lwf,&
     &                         q_crit, cfg_params

USE mo_2mom_mcrph_config_default, ONLY: cfg_2mom_default

USE mo_2mom_mcrph_util, ONLY:                            &
     &                       init_dmin_wg_gr_ltab_equi,  &
     &                       dmin_wetgrowth_fit_check, luse_dmin_wetgrowth_table, lprintout_comp_table_fit

USE mo_2mom_mcrph_types, ONLY: ltabdminwgg, ltabdminwgh

USE mo_2mom_prepare, ONLY: prepare_twomoment, post_twomoment
USE mo_nwp_tuning_config,  ONLY: tune_sbmccn
USE mo_fortran_tools, ONLY: init
!==============================================================================

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_driver'
  INTEGER,          PARAMETER :: dbg_level = 25                   ! level for debug prints

  ! .. exponents for simple density of terminal fall velocity
  REAL(wp), PARAMETER :: rho_vel    = 0.4e0_wp    !..exponent for density correction
  REAL(wp), PARAMETER :: rho_vel_c  = 0.2e0_wp    !..for cloud droplets (in exact terms, this would be the ratio of dyn. visc. as function of T)
  REAL(wp), PARAMETER :: rho0       = 1.225_wp    !..reference air density

  INTEGER :: cloud_type, ccn_type

  LOGICAL :: lconstant_lh   ! Use constant latent heat (default .true.)

! UB: These settings should be converted into namelist parameters in the future!

!! Now in namelist phy_ctl!  INTEGER, PARAMETER :: i2mom_solver = 1  ! (0) explicit (1) semi-implicit solve
!!$  ! now this comes from cfg_params !  INTEGER, PARAMETER :: i2mom_solver = 1  ! (0) explicit (1) semi-implicit solve
  
  INTEGER, PARAMETER :: cloud_type_default_gscp4 = 2603, ccn_type_gscp4 = 7
  INTEGER, PARAMETER :: cloud_type_default_gscp5 = 2603, ccn_type_gscp5 = 8

  ! AS: For gscp=4 use 2103 with ccn_type = 1 (HDCP2 IN and CCN schemes)
  !     For gscp=5 use 2603 with ccn_type = 8 (PDA ice nucleation and Segal&Khain CCN activation)
  
  ! AS: Runs without hail, e.g, 1503 are buggy and give a segmentation fault.
  !     So far I was not able to identify the problem, needs more detailed debugging.
  
CONTAINS
  
  !==============================================================================
  !
  ! Two-moment mixed-phase bulk microphysics
  !
  ! original version by Axel Seifert, May 2003
  ! with modifications by Ulrich Blahak, August 2007

  !==============================================================================
  SUBROUTINE two_moment_mcrph(            &
                       isize,             & ! in: array size
                       ke,                & ! in: end level/array size
                       is,                & ! in: start index, optional
                       ie,                & ! in: end index, optional
                       ks,                & ! in: start index vertical , optional
                       dt,                & ! in: time step
                       dz,                & ! in: vertical layer thickness
                       hhl,               & ! in: height of half levels
                       rho,               & ! in: density
                       pres,              & ! in: pressure
                       tke,               & ! in: tke
                       qv,                & ! inout: specific humidity
                       qc, qnc,           & ! inout: cloud water
                       qr, qnr,           & ! inout: rain
                       qi, qni,           & ! inout: ice
                       qs, qns,           & ! inout: snow
                       qg, qng, qgl,      & ! inout: graupel
                       qh, qnh, qhl,      & ! inout: hail
                       nccn,              & ! inout: ccn
                       ninpot,            & ! inout: potential ice nuclei
                       ninact,            & ! inout: activated ice nuclei
                       tk,                & ! inout: temp
                       w,                 & ! inout: w
                       prec_r,            & ! inout: precip rate rain
                       prec_i,            & ! inout: precip rate ice
                       prec_s,            & ! inout: precip rate snow
                       prec_g,            & ! inout: precip rate graupel
                       prec_h,            & ! inout: precip rate hail
                       qrsflux,           & ! inout: 3D total precipitation rate
                       dtemp,             & ! inout: opt. temp increment
                       msg_level,         & ! in: msg_level
                       l_cv,              & ! in: switch for cv/cp
                       ithermo_water      ) ! in: thermodynamic option

    ! Declare variables in argument list
    
    INTEGER,            INTENT (IN)  :: isize, ke    ! grid sizes
    INTEGER,  OPTIONAL, INTENT (IN)  :: is, ie, ks   ! start/end indices

    REAL(wp), INTENT (IN)            :: dt           ! time step

    ! Dynamical core variables
    REAL(wp), DIMENSION(:,:), INTENT(IN), TARGET :: dz, rho, pres, w

    ! Optional Dynamical core variables
    REAL(wp), DIMENSION(:,:), INTENT(IN), POINTER :: tke

    REAL(wp), DIMENSION(:,:), INTENT(IN), TARGET :: hhl

    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET :: tk

    ! Microphysics variables
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) , TARGET :: &
         qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, ninact

    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               qgl, qhl

    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               nccn, ninpot

    ! Precip rates, vertical profiles
    REAL(wp), DIMENSION(:), INTENT (INOUT) :: &
         &               prec_r, prec_i, prec_s, prec_g, prec_h
    REAL(wp), DIMENSION(:,:), INTENT (INOUT) :: qrsflux
    
    REAL(wp), OPTIONAL, INTENT (INOUT)  :: dtemp(:,:)

    INTEGER,  INTENT (IN)             :: msg_level

    LOGICAL,  OPTIONAL,  INTENT (IN)  :: l_cv

    INTEGER,  OPTIONAL,  INTENT (IN)  :: ithermo_water

    ! ... Variables which are global in module_2mom_mcrph_main

    REAL(wp), TARGET, DIMENSION(isize,ke) ::        &
         &  rhocorr,       & ! density dependency of particle fall speed
         &  rhocld           ! density dependency of particle fall speed for cloud droplets

    REAL(wp) :: q_liq_old(isize,ke), q_vap_old(isize,ke)  ! to store old values for latent heat calc

    INTEGER  :: its,ite,kts,kte
    INTEGER  :: ii,kk
    INTEGER  :: ntsedi_rain, ntsedi_graupel, ntsedi_hail     ! for sedimentation sub stepping

    REAL(wp) :: q_liq_new,q_vap_new
    REAL(wp) :: zf,hlp,dtemp_loc
    REAL(wp) :: convliq,convice,led,lwe
    REAL(wp), PARAMETER :: tau_inact =  600.  ! relaxation time scale for activated IN number density
    REAL(wp), PARAMETER :: tau_inpot = 1800.  ! relaxation time scale for potential IN number density
    REAL(wp) :: in_bgrd            ! background profile of IN number density
    REAL(wp) :: z_heat_cap_r       ! reciprocal of cpdr or cvdr (depending on l_cv)
    REAL(wp) :: rdz(isize,ke), rho_r(isize,ke)

    LOGICAL :: lprogccn, lprogin, lprogmelt

    LOGICAL, PARAMETER :: debug     = .false.       !
    LOGICAL, PARAMETER :: clipping  = .true.        ! not really necessary, just for cleanup

    CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_driver'

    ! These structures include the pointers to the model arrays (which are automatic arrays
    ! of this driver subroutine). These structures live only for one time step and are
    ! different for the OpenMP threads. In contrast, the types like rain_coeffs, ice_coeffs,
    ! etc. that are declared in mo_2mom_mcrph_main live for the whole runtime and may include
    ! coefficients that are calculated once during initialization (and there is only one per
    ! mpi thread).
    TYPE(atmosphere)           :: atmo

    ! These are the fundamental hydrometeor particle variables for the two-moment scheme
    ! as they are used in the various options of the scheme
    TYPE(particle), target          :: cloud_hyd, rain_hyd
    TYPE(particle_frozen), target   :: ice_frz, snow_frz, graupel_frz, hail_frz
    TYPE(particle_lwf), target      :: graupel_lwf, hail_lwf

    ! Pointers to the derived types that are actually needed
    CLASS(particle), pointer        :: cloud, rain
    CLASS(particle_frozen), pointer :: ice, snow, graupel, hail

    INTEGER :: ik_slice(4)


    !$ACC DATA CREATE(q_liq_old, q_vap_old, rdz, rhocorr, rho_r, rhocld)

    lprogccn  = PRESENT(nccn)
    lprogin   = PRESENT(ninpot)
    lprogmelt = PRESENT(qgl)
    
#ifdef _OPENACC
    IF (lprogmelt) THEN
      CALL finish(routine, 'lprogmelt not available on GPU for two-moment microphysics')
    ENDIF
#endif

    IF (msg_level>5) CALL message (TRIM(routine), "called two_moment_mcrph")

    IF (PRESENT(ithermo_water)) THEN
       lconstant_lh = (ithermo_water == 0)
    ELSE  ! Default themodynamic is constant latent heat
       lconstant_lh = .true.
    END IF

    IF (lprogccn) THEN
       cloud_type = cloud_type_default_gscp5 + 10 * ccn_type
    ELSE
       cloud_type = cloud_type_default_gscp4 + 10 * ccn_type
    END IF

    cloud => cloud_hyd
    rain  => rain_hyd
    ice   => ice_frz
    snow  => snow_frz

    IF (lprogmelt) THEN
       graupel => graupel_lwf   ! with prognostic melting 
       hail => hail_lwf         ! of graupel and hail
    ELSE
       graupel => graupel_frz   ! simple melting
       hail => hail_frz
    END IF

    ! start/end indices
    IF (PRESENT(is)) THEN
      its = is
    ELSE
      its = 1
    END IF
    IF (PRESENT(ie)) THEN
      ite = ie
    ELSE
      ite = isize
    END IF
    IF (PRESENT(ks)) THEN
      kts = ks
    ELSE
      kts = 1
    END IF
    kte = ke

    IF (timers_level > 10) CALL timer_start(timer_phys_2mom_prepost) 

    ! inverse of vertical layer thickness
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO kk = kts,kte
      DO ii = its,ite
        rdz(ii,kk) = 1._wp / dz(ii,kk)

        IF (clipping) THEN
              IF(qr(ii,kk) < 0.0_wp) qr(ii,kk) = 0.0_wp
              IF(qi(ii,kk) < 0.0_wp) qi(ii,kk) = 0.0_wp
              IF(qs(ii,kk) < 0.0_wp) qs(ii,kk) = 0.0_wp
              IF(qg(ii,kk) < 0.0_wp) qg(ii,kk) = 0.0_wp
              IF(qh(ii,kk) < 0.0_wp) qh(ii,kk) = 0.0_wp
        END IF
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! Initialize qrsflux for LHN:
    IF (ldass_lhn) THEN
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO kk = kts,kte
        DO ii = its,ite
          qrsflux(ii,kk) = 0.0_wp
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    END IF

    IF (clipping) THEN
       IF (lprogmelt) THEN
          WHERE(qgl(its:ite,kts:kte) < 0.0_wp) qgl(its:ite,kts:kte) = 0.0_wp
          WHERE(qhl(its:ite,kts:kte) < 0.0_wp) qhl(its:ite,kts:kte) = 0.0_wp
       END IF
    END IF
    
    ! indices as used in two-moment scheme
    ik_slice(1) = its
    ik_slice(2) = ite
    ik_slice(3) = kts
    ik_slice(4) = kte

    IF (PRESENT(l_cv)) THEN
      IF (l_cv) THEN
        z_heat_cap_r = cvdr
      ELSE
        z_heat_cap_r = cpdr
      ENDIF
    ELSE
      z_heat_cap_r = cpdr
    ENDIF

    IF (msg_level>dbg_level) CALL message(TRIM(routine),'')

    IF (msg_level>dbg_level)THEN
       WRITE (message_text,'(1X,A,I4,3(A,L2))') &
            & "cloud_type = ",cloud_type,", lprogccn = ",lprogccn,", lprogin = ",lprogin,", lprogmelt = ",lprogmelt
       CALL message(TRIM(routine),TRIM(message_text))
    END IF

    IF (msg_level>dbg_level) CALL message(TRIM(routine), "prepare variables for 2mom")

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE(hlp) COLLAPSE(2)
    DO kk = kts, kte
       DO ii = its, ite

          ! ... 1/rho is used quite often
          rho_r(ii,kk) = 1.0 / rho(ii,kk)

          ! ... height dependency of terminal fall velocities
          hlp = LOG(MAX(rho(ii,kk),1e-6_wp)/rho0)
          rhocorr(ii,kk) = exp(-rho_vel*hlp)
          rhocld(ii,kk)  = exp(-rho_vel_c*hlp)

       END DO
    END DO
    !$ACC END PARALLEL

    ! .. set the particle types, but no calculations
    CALL init_2mom_scheme(cloud,rain,ice,snow,graupel,hail)
    !$ACC DATA COPYIN(cloud, rain, ice, snow, graupel, hail) &
    !$ACC   CREATE(atmo)

    ! .. convert to densities and set pointerns to two-moment module
    !    (pointers are used to avoid passing everything explicitly by argument and
    !     to avoid local allocates within the OpenMP-loop, and keep everything on stack)

    CALL prepare_twomoment(atmo, cloud, rain, ice, snow, graupel, hail, &
         rho, rhocorr, rhocld, pres, w, tk, hhl, tke, &
         nccn, ninpot, ninact, &
         qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, qgl, qhl, &
         lprogccn, lprogin, lprogmelt, its, ite, kts, kte)
    IF (timers_level > 10) CALL timer_stop(timer_phys_2mom_prepost)

    IF (msg_level>dbg_level) CALL message(TRIM(routine)," calling clouds_twomoment")

    IF (cfg_params%i2mom_solver.eq.0) THEN

#ifdef _OPENACC
      IF(PRESENT(dtemp)) CALL finish('clouds_twomoment:','dtemp not available on GPU')
#endif
      ! ... save old variables for latent heat calculation
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG COLLAPSE(2)
      DO kk = kts, kte
        DO ii = its, ite
          q_vap_old(ii,kk) = qv(ii,kk)

          IF (.NOT. lprogmelt) THEN
               q_liq_old(ii,kk) = qc(ii,kk) + qr(ii,kk)
          END IF

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      IF (lprogmelt) THEN
         q_liq_old(its:ite,kts:kte) = qc(its:ite,kts:kte) + qr(its:ite,kts:kte)  &
              &                     + qgl(its:ite,kts:kte) + qhl(its:ite,kts:kte)
      END IF

       IF (timers_level > 10) CALL timer_start(timer_phys_2mom_proc) 
       ! .. this subroutine calculates all the microphysical sources and sinks
       CALL clouds_twomoment(ik_slice, dt, lprogin, &
            atmo, cloud, rain, ice, snow, graupel, hail, ninact, nccn, ninpot)

       IF (timers_level > 10) CALL timer_stop(timer_phys_2mom_proc) 

       IF (lprogccn) THEN
#ifdef _OPENACC
        CALL finish('clouds_twomoment:','lprogccn not available on GPU')
#endif
        WHERE(qc(its:ite,kts:kte) == 0.0_wp) cloud%n(its:ite,kts:kte) = 0.0_wp
       END IF

       !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(led, lwe, convice, convliq, q_liq_new, dtemp_loc)
       DO kk = kts, kte
          DO ii = its, ite

            IF (lconstant_lh) THEN
              led = als
              lwe = (alv-als)
            ELSE
              led = latent_heat_sublimation(tk(ii,kk))
              lwe = latent_heat_melting(tk(ii,kk))
            END IF

            ! .. latent heat term for temperature equation
            convice = z_heat_cap_r * led
            convliq = z_heat_cap_r * lwe
            
            ! .. new variables
            q_vap_new = qv(ii,kk)
            if (lprogmelt) then
              q_liq_new = qr(ii,kk) + qc(ii,kk) + qgl(ii,kk) + qhl(ii,kk) 
            else
              q_liq_new = qr(ii,kk) + qc(ii,kk)
            end if
             
            ! .. update temperature
            dtemp_loc  = - convice * rho_r(ii,kk) * (q_vap_new - q_vap_old(ii,kk))  &
                 &       + convliq * rho_r(ii,kk) * (q_liq_new - q_liq_old(ii,kk))

            tk(ii,kk) = tk(ii,kk) + dtemp_loc

            IF(PRESENT(dtemp)) &
                 dtemp(ii,kk) = dtemp_loc

          ENDDO
       ENDDO
       !$ACC END PARALLEL

       IF (timers_level > 10) CALL timer_start(timer_phys_2mom_sedi) 

       IF (msg_level>dbg_level) CALL message(TRIM(routine)," calling sedimentation")

       ! .. if we solve explicitly, then sedimentation is done here after microphysics
       CALL sedimentation_explicit()
       IF (timers_level > 10) CALL timer_stop(timer_phys_2mom_sedi) 

    ELSE
      CALL clouds_twomoment_implicit ()
    END IF

    IF (timers_level > 10) CALL timer_start(timer_phys_2mom_prepost)    
    
    ! .. check for negative values
    IF (debug) CALL check_clouds()

    ! .. convert back and nullify two-moment pointers
    CALL post_twomoment(atmo, cloud, rain, ice, snow, graupel, hail, &
         rho_r, qnc, nccn, ninpot, ninact, &
         qv, qc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, qgl, qhl,  &
         lprogccn, lprogin, lprogmelt, its, ite, kts, kte)

    IF (clipping) THEN
      IF (lprogmelt) THEN
        WHERE(qgl(its:ite,kts:kte) < 0.0_wp) qgl(its:ite,kts:kte) = 0.0_wp
        WHERE(qhl(its:ite,kts:kte) < 0.0_wp) qhl(its:ite,kts:kte) = 0.0_wp
      END IF
    END IF

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO kk = kts,kte

      IF (clipping) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ii = its,ite
          IF ( qr(ii,kk) < 0.0_wp)  qr(ii,kk) = 0.0_wp
          IF ( qi(ii,kk) < 0.0_wp)  qi(ii,kk) = 0.0_wp
          IF ( qs(ii,kk) < 0.0_wp)  qs(ii,kk) = 0.0_wp
          IF ( qg(ii,kk) < 0.0_wp)  qg(ii,kk) = 0.0_wp
          IF ( qh(ii,kk) < 0.0_wp)  qh(ii,kk) = 0.0_wp
          IF (qnr(ii,kk) < 0.0_wp) qnr(ii,kk) = 0.0_wp
          IF (qni(ii,kk) < 0.0_wp) qni(ii,kk) = 0.0_wp
          IF (qns(ii,kk) < 0.0_wp) qns(ii,kk) = 0.0_wp
          IF (qng(ii,kk) < 0.0_wp) qng(ii,kk) = 0.0_wp
          IF (qnh(ii,kk) < 0.0_wp) qnh(ii,kk) = 0.0_wp
        ENDDO
      END IF

      IF (lprogccn) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zf)
        DO ii=its,ite
          zf = 0.5_wp*(hhl(ii,kk)+hhl(ii,kk+1))
          !..reset nccn for cloud-free grid points to background profile
          IF (qc(ii,kk) .LE. q_crit) THEN
            IF(zf > ccn_coeffs%z0) THEN
              nccn(ii,kk) = MAX(nccn(ii,kk),ccn_coeffs%Ncn0 &
                   * EXP((ccn_coeffs%z0 - zf)*(1._wp/ccn_coeffs%z1e)))
            ELSE
              nccn(ii,kk) = MAX(nccn(ii,kk),ccn_coeffs%Ncn0)
            END IF
          END IF
        END DO
      END IF

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ii=its,ite
        !..relaxation of activated IN number density to zero
        IF(qi(ii,kk) == 0) THEN
          ninact(ii,kk) = ninact(ii,kk) - ninact(ii,kk)*(1._wp/tau_inact)*dt
        END IF
      END DO

      IF (lprogin) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zf, in_bgrd)
        DO ii=its,ite
          zf = 0.5_wp*(hhl(ii,kk)+hhl(ii,kk+1))
          !..relaxation of potential IN number density to background profile
          IF(zf > in_coeffs%z0) THEN
            in_bgrd = in_coeffs%N0*EXP((in_coeffs%z0 - zf)*(1._wp/in_coeffs%z1e))
          ELSE
            in_bgrd = in_coeffs%N0
          END IF
          ninpot(ii,kk) = ninpot(ii,kk) - (ninpot(ii,kk)-in_bgrd)*(1._wp/tau_inpot)*dt
        END DO
      END IF

      IF (lprogccn) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ii = its,ite
          IF ( nccn(ii,kk) < 35e6_wp ) nccn(ii,kk) = 35e6_wp
        ENDDO
      END IF
  
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ii = its,ite
        IF(qc(ii,kk) < 1.0e-12_wp) qnc(ii,kk) = 0.0_wp
      ENDDO

    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA ! cloud, rain, ice, snow, graupel, hail, atmo
    !$ACC END DATA ! q_liq_old, q_vap_old, rdz, rhocorr, rho_r, rhocld

    IF (msg_level>dbg_level) CALL message(TRIM(routine), "two moment mcrph ends!")

    IF (timers_level > 10) CALL timer_stop(timer_phys_2mom_prepost) 

    RETURN
    !
    ! end of driver routine, but many details are below in the contains-part of this subroutine
    !
  CONTAINS

    SUBROUTINE clouds_twomoment_implicit()
      !
      ! semi-implicit solver for sedimentation including microphysics, the same
      ! approach is used in the COSMO microphysics, e.g, hydci_pp
      ! (see COSMO documentation for details)
      !
      
      ! a few 1d arrays, maybe we can reduce this later or we keep them ...
      real(wp), dimension(isize) :: &
           & qr_flux_now,qr_flux_new,qr_sum,vr_sedq_new,vr_sedq_now,qr_impl,xr_now, &
           & nr_flux_now,nr_flux_new,nr_sum,vr_sedn_new,vr_sedn_now,nr_impl,        &
           & qs_flux_now,qs_flux_new,qs_sum,vs_sedq_new,vs_sedq_now,qs_impl,xs_now, &
           & ns_flux_now,ns_flux_new,ns_sum,vs_sedn_new,vs_sedn_now,ns_impl,        &
           & qg_flux_now,qg_flux_new,qg_sum,vg_sedq_new,vg_sedq_now,qg_impl,xg_now, &
           & ng_flux_now,ng_flux_new,ng_sum,vg_sedn_new,vg_sedn_now,ng_impl,        &
           & qh_flux_now,qh_flux_new,qh_sum,vh_sedq_new,vh_sedq_now,qh_impl,xh_now, &
           & nh_flux_now,nh_flux_new,nh_sum,vh_sedn_new,vh_sedn_now,nh_impl,        &
           & qi_flux_now,qi_flux_new,qi_sum,vi_sedq_new,vi_sedq_now,qi_impl,xi_now, &
           & ni_flux_now,ni_flux_new,ni_sum,vi_sedn_new,vi_sedn_now,ni_impl         

      ! for lwf variables
      real(wp), dimension(isize) :: &
           & lh_flux_now,lh_flux_new,lh_sum,vh_sedl_new,vh_sedl_now,lh_impl, &
           & lg_flux_now,lg_flux_new,lg_sum,vg_sedl_new,vg_sedl_now,lg_impl

      REAL(wp), DIMENSION(isize,ke) :: rdzdt
      INTEGER :: i, ii, k, kk

      logical, parameter :: lmicro_impl = .true.  ! microphysics within semi-implicit sedimentation loop?

#ifdef _OPENACC
    IF (lprogmelt) THEN
      CALL finish(routine, 'lprogmelt not available on GPU for two-moment microphysics')
    ENDIF
#endif

      !$ACC DATA &
      !$ACC   CREATE(qr_flux_now, qr_flux_new, qr_sum, vr_sedq_new, vr_sedq_now, qr_impl, xr_now) &
      !$ACC   CREATE(nr_flux_now, nr_flux_new, nr_sum, vr_sedn_new, vr_sedn_now, nr_impl) &
      !$ACC   CREATE(qs_flux_now, qs_flux_new, qs_sum, vs_sedq_new, vs_sedq_now, qs_impl, xs_now) &
      !$ACC   CREATE(ns_flux_now, ns_flux_new, ns_sum, vs_sedn_new, vs_sedn_now, ns_impl) &
      !$ACC   CREATE(qg_flux_now, qg_flux_new, qg_sum, vg_sedq_new, vg_sedq_now, qg_impl, xg_now) &
      !$ACC   CREATE(ng_flux_now, ng_flux_new, ng_sum, vg_sedn_new, vg_sedn_now, ng_impl) &
      !$ACC   CREATE(qh_flux_now, qh_flux_new, qh_sum, vh_sedq_new, vh_sedq_now, qh_impl, xh_now) &
      !$ACC   CREATE(nh_flux_now, nh_flux_new, nh_sum, vh_sedn_new, vh_sedn_now, nh_impl) &
      !$ACC   CREATE(qi_flux_now, qi_flux_new, qi_sum, vi_sedq_new, vi_sedq_now, qi_impl, xi_now) &
      !$ACC   CREATE(ni_flux_now, ni_flux_new, ni_sum, vi_sedn_new, vi_sedn_now, ni_impl) &
      !$ACC   CREATE(lh_flux_now, lh_flux_new, lh_sum, vh_sedl_new, vh_sedl_now, lh_impl) &
      !$ACC   CREATE(lg_flux_now, lg_flux_new, lg_sum, vg_sedl_new, vg_sedl_now, lg_impl, rdzdt)

      if (.not.lmicro_impl) then
#ifdef _OPENACC
        CALL finish('clouds_twomoment_implicit', 'routine without lmicro_impl not available on GPU')
#endif        

        ! ... save old variables for latent heat calculation
        if (lprogmelt) then
          q_vap_old(its:ite,kts:kte) = qv(its:ite,kts:kte)
          q_liq_old(its:ite,kts:kte) = qc(its:ite,kts:kte) + qgl(its:ite,kts:kte) &
               &                     + qr(its:ite,kts:kte) + qhl(its:ite,kts:kte)
        else
          q_vap_old(its:ite,kts:kte) = qv(its:ite,kts:kte)
          q_liq_old(its:ite,kts:kte) = qc(its:ite,kts:kte) + qr(its:ite,kts:kte)
        end if

        ! .. this subroutine calculates all the microphysical sources and sinks
        CALL clouds_twomoment(ik_slice, dt, lprogin, atmo, cloud, rain, &
             ice, snow, graupel, hail, ninact, nccn, ninpot) 

        DO kk=kts,kte
          DO ii = its, ite
            ! .. latent heat term for temperature equation
            IF (lconstant_lh) THEN
              led = als
              lwe = (alv-als)
            ELSE
              led = latent_heat_sublimation(tk(ii,kk))
              lwe = latent_heat_melting(tk(ii,kk))
            END IF

            convice = z_heat_cap_r * led
            convliq = z_heat_cap_r * lwe
            
            q_vap_new = qv(ii,kk)
            if (lprogmelt) then
              q_liq_new = qr(ii,kk) + qc(ii,kk) + qgl(ii,kk) + qhl(ii,kk)
            else
              q_liq_new = qr(ii,kk) + qc(ii,kk)
            end if
            tk(ii,kk) = tk(ii,kk) - convice * rho_r(ii,kk) * (q_vap_new - q_vap_old(ii,kk))  &
                 &                + convliq * rho_r(ii,kk) * (q_liq_new - q_liq_old(ii,kk))
          ENDDO
        ENDDO

      end if

      ! clipping maybe not necessary
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k = kts,kte
        DO i = its,ite
          IF(qr(i,k) < 0.0_wp) qr(i,k) = 0.0_wp
          IF(qi(i,k) < 0.0_wp) qi(i,k) = 0.0_wp
          IF(qs(i,k) < 0.0_wp) qs(i,k) = 0.0_wp
          IF(qg(i,k) < 0.0_wp) qg(i,k) = 0.0_wp
          IF(qh(i,k) < 0.0_wp) qh(i,k) = 0.0_wp
          IF(qnr(i,k) < 0.0_wp) qnr(i,k) = 0.0_wp
          IF(qni(i,k) < 0.0_wp) qni(i,k) = 0.0_wp
          IF(qns(i,k) < 0.0_wp) qns(i,k) = 0.0_wp
          IF(qng(i,k) < 0.0_wp) qng(i,k) = 0.0_wp
          IF(qnh(i,k) < 0.0_wp) qnh(i,k) = 0.0_wp
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      if (lprogmelt) then
        WHERE(qgl(its:ite,kts:kte) < 0.0_wp) qgl(its:ite,kts:kte) = 0.0_wp
        WHERE(qhl(its:ite,kts:kte) < 0.0_wp) qhl(its:ite,kts:kte) = 0.0_wp
      end if

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k = kts,kte
        DO i = its,ite
            rdzdt(i,k) = 0.5_wp * rdz(i,k) * dt
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO i = its,ite
        qr_flux_now(i) = 0.0_wp
        nr_flux_now(i) = 0.0_wp
        qr_flux_new(i) = 0.0_wp
        nr_flux_new(i) = 0.0_wp

        qi_flux_now(i) = 0.0_wp
        ni_flux_now(i) = 0.0_wp
        qi_flux_new(i) = 0.0_wp
        ni_flux_new(i) = 0.0_wp

        qs_flux_now(i) = 0.0_wp
        ns_flux_now(i) = 0.0_wp
        qs_flux_new(i) = 0.0_wp
        ns_flux_new(i) = 0.0_wp

        qg_flux_now(i) = 0.0_wp
        ng_flux_now(i) = 0.0_wp
        qg_flux_new(i) = 0.0_wp
        ng_flux_new(i) = 0.0_wp

        qh_flux_now(i) = 0.0_wp
        nh_flux_now(i) = 0.0_wp
        qh_flux_new(i) = 0.0_wp
        nh_flux_new(i) = 0.0_wp
      ENDDO
      !$ACC END PARALLEL

      if (lprogmelt) then
        lg_flux_now(:) = 0.0_wp
        lg_flux_new(:) = 0.0_wp        
        lh_flux_now(:) = 0.0_wp
        lh_flux_new(:) = 0.0_wp        
      end if

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      do i=its,ite
        vr_sedn_new(i) = rain%vsedi_min
        vi_sedn_new(i) = ice%vsedi_min
        vs_sedn_new(i) = snow%vsedi_min
        vg_sedn_new(i) = graupel%vsedi_min
        vh_sedn_new(i) = hail%vsedi_min
        vr_sedq_new(i) = rain%vsedi_min
        vi_sedq_new(i) = ice%vsedi_min
        vs_sedq_new(i) = snow%vsedi_min
        vg_sedq_new(i) = graupel%vsedi_min
        vh_sedq_new(i) = hail%vsedi_min
      end do
      !$ACC END PARALLEL

      if (lprogmelt) then 
        do i=its,ite
          vg_sedl_new(i) = graupel%vsedi_min
          vh_sedl_new(i) = hail%vsedi_min
        end do
      end if

      ! here we simply assume that there is no cloud or precip in the uppermost level
      ! i.e. we start from kts+1 going down in physical space

!DIR$ IVDEP
      DO k=kts+1,kte

        !$ACC DATA PRESENT(rain, ice, snow, graupel, hail)
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        do i=its,ite
          xr_now(i) = particle_meanmass(rain, qr(i,k),qnr(i,k))
          xi_now(i) = particle_meanmass(ice, qi(i,k),qni(i,k))
          xs_now(i) = particle_meanmass(snow, qs(i,k),qns(i,k))
          xg_now(i) = particle_meanmass(graupel, qg(i,k),qng(i,k))
          xh_now(i) = particle_meanmass(hail, qh(i,k),qnh(i,k))
        end do
        !$ACC END PARALLEL
        !$ACC END DATA

        call sedi_vel_rain(rain,rain_coeffs,qr(:,k),xr_now,rhocorr(:,k),vr_sedn_now,vr_sedq_now,its,ite,qc(:,k),lacc=.TRUE.)
        call sedi_vel_sphere(ice,ice_coeffs,qi(:,k),xi_now,rhocorr(:,k),vi_sedn_now,vi_sedq_now,its,ite)
        call sedi_vel_sphere(snow,snow_coeffs,qs(:,k),xs_now,rhocorr(:,k),vs_sedn_now,vs_sedq_now,its,ite)
        if (lprogmelt) then
          call sedi_vel_lwf(graupel_lwf,graupel_coeffs,  &
               & qg(:,k),qgl(:,k),xg_now,rhocorr(:,k),vg_sedn_now,vg_sedq_now,vg_sedl_now,its,ite)
          call sedi_vel_lwf(hail_lwf,hail_coeffs,        &
               & qh(:,k),qhl(:,k),xh_now,rhocorr(:,k),vh_sedn_now,vh_sedq_now,vh_sedl_now,its,ite)
        else
          call sedi_vel_sphere(graupel,graupel_coeffs,qg(:,k),xg_now,rhocorr(:,k),vg_sedn_now,vg_sedq_now,its,ite)
          call sedi_vel_sphere(hail,hail_coeffs,qh(:,k),xh_now,rhocorr(:,k),vh_sedn_now,vh_sedq_now,its,ite)
        end if

        call implicit_core(qr(:,k), qr_sum,qr_impl,vr_sedq_new,vr_sedq_now,qr_flux_new,qr_flux_now,rdzdt(:,k),its,ite)
        call implicit_core(qnr(:,k),nr_sum,nr_impl,vr_sedn_new,vr_sedn_now,nr_flux_new,nr_flux_now,rdzdt(:,k),its,ite)
        call implicit_core(qi(:,k), qi_sum,qi_impl,vi_sedq_new,vi_sedq_now,qi_flux_new,qi_flux_now,rdzdt(:,k),its,ite)
        call implicit_core(qni(:,k),ni_sum,ni_impl,vi_sedn_new,vi_sedn_now,ni_flux_new,ni_flux_now,rdzdt(:,k),its,ite)
        call implicit_core(qs(:,k), qs_sum,qs_impl,vs_sedq_new,vs_sedq_now,qs_flux_new,qs_flux_now,rdzdt(:,k),its,ite)
        call implicit_core(qns(:,k),ns_sum,ns_impl,vs_sedn_new,vs_sedn_now,ns_flux_new,ns_flux_now,rdzdt(:,k),its,ite)
        call implicit_core(qg(:,k), qg_sum,qg_impl,vg_sedq_new,vg_sedq_now,qg_flux_new,qg_flux_now,rdzdt(:,k),its,ite)
        call implicit_core(qng(:,k),ng_sum,ng_impl,vg_sedn_new,vg_sedn_now,ng_flux_new,ng_flux_now,rdzdt(:,k),its,ite)
        call implicit_core(qh(:,k), qh_sum,qh_impl,vh_sedq_new,vh_sedq_now,qh_flux_new,qh_flux_now,rdzdt(:,k),its,ite)
        call implicit_core(qnh(:,k),nh_sum,nh_impl,vh_sedn_new,vh_sedn_now,nh_flux_new,nh_flux_now,rdzdt(:,k),its,ite)
        
        if (lprogmelt) then
          call implicit_core(qgl(:,k),lg_sum,lg_impl,vg_sedl_new,vg_sedl_now,lg_flux_new,lg_flux_now,rdzdt(:,k),its,ite)
          call implicit_core(qhl(:,k),lh_sum,lh_impl,vh_sedl_new,vh_sedl_now,lh_flux_new,lh_flux_now,rdzdt(:,k),its,ite)
        end if

        ! do microphysics on this k-level only (using the star-values)
        IF (lmicro_impl) THEN

          ! .. save old variables for latent heat calculation

          !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO ii = its, ite
            q_vap_old(ii,k) = qv(ii,k)
            if (lprogmelt) then
#ifndef _OPENACC
              q_liq_old(ii,k) = qr(ii,k) + qc(ii,k) + qgl(ii,k) + qhl(ii,k)
#endif
            else
              q_liq_old(ii,k) = qc(ii,k) + qr(ii,k)
            end if
          END DO
          !$ACC END PARALLEL

          ik_slice(3) = k
          ik_slice(4) = k
          CALL clouds_twomoment(ik_slice, dt, lprogin, &
               atmo, cloud, rain, ice, snow, graupel, hail, &
               ninact, nccn, ninpot)

          ! .. latent heat term for temperature equation
          !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(led, lwe, convice, convliq, q_liq_new, q_vap_new)
          DO ii = its, ite

            IF (lconstant_lh) THEN
              led = als
              lwe = (alv-als)
            ELSE
              led = latent_heat_sublimation(tk(ii,k))
              lwe = latent_heat_melting(tk(ii,k))
            END IF

            convice = z_heat_cap_r * led
            convliq = z_heat_cap_r * lwe

            q_vap_new  = qv(ii,k)
            if (lprogmelt) then
#ifndef _OPENACC
              q_liq_new = qr(ii,k) + qc(ii,k) + qgl(ii,k) + qhl(ii,k)
#endif
            else
              q_liq_new = qr(ii,k) + qc(ii,k)
            end if
            tk(ii,k)   = tk(ii,k) - convice * rho_r(ii,k) * (q_vap_new - q_vap_old(ii,k))  &
                 &                + convliq * rho_r(ii,k) * (q_liq_new - q_liq_old(ii,k))
          END DO
          !$ACC END PARALLEL

        END IF

        call implicit_time(qr(:,k), qr_sum,qr_impl,vr_sedq_new,vr_sedq_now,qr_flux_new,its,ite)
        call implicit_time(qnr(:,k),nr_sum,nr_impl,vr_sedn_new,vr_sedn_now,nr_flux_new,its,ite)
        call implicit_time(qi(:,k), qi_sum,qi_impl,vi_sedq_new,vi_sedq_now,qi_flux_new,its,ite)
        call implicit_time(qni(:,k),ni_sum,ni_impl,vi_sedn_new,vi_sedn_now,ni_flux_new,its,ite)
        call implicit_time(qs(:,k), qs_sum,qs_impl,vs_sedq_new,vs_sedq_now,qs_flux_new,its,ite)
        call implicit_time(qns(:,k),ns_sum,ns_impl,vs_sedn_new,vs_sedn_now,ns_flux_new,its,ite)
        call implicit_time(qg(:,k), qg_sum,qg_impl,vg_sedq_new,vg_sedq_now,qg_flux_new,its,ite)
        call implicit_time(qng(:,k),ng_sum,ng_impl,vg_sedn_new,vg_sedn_now,ng_flux_new,its,ite)
        call implicit_time(qh(:,k), qh_sum,qh_impl,vh_sedq_new,vh_sedq_now,qh_flux_new,its,ite)
        call implicit_time(qnh(:,k),nh_sum,nh_impl,vh_sedn_new,vh_sedn_now,nh_flux_new,its,ite)
        
        if (lprogmelt) then
          call implicit_time(qgl(:,k),lg_sum,lg_impl,vg_sedl_new,vg_sedl_now,lg_flux_new,its,ite)
          call implicit_time(qhl(:,k),lh_sum,lh_impl,vh_sedl_new,vh_sedl_now,lh_flux_new,its,ite)
        end if

        IF (ldass_lhn) THEN
          IF (lprogmelt) THEN
            qrsflux(:,k) = qr_flux_new + qi_flux_new + qs_flux_new + qg_flux_new + qh_flux_new + &
                           lg_flux_new + lh_flux_new
          ELSE
            !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO ii = its, ite
              qrsflux(ii,k) = qr_flux_new(ii) + qi_flux_new(ii) + qs_flux_new(ii) + qg_flux_new(ii) + qh_flux_new(ii)
            ENDDO
            !$ACC END PARALLEL
          END IF
        END IF

      END DO

      IF (lprogmelt) THEN
        ! implicit solver for LWF-scheme still has some issues
        prec_g(:) = MAX( qg_flux_new + lg_flux_new, 0.0_wp )
        prec_h(:) = MAX( qh_flux_new + lh_flux_new, 0.0_wp )
        prec_r(:) = qr_flux_new
        prec_i(:) = qi_flux_new
        prec_s(:) = qs_flux_new
      ELSE
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO i = its,ite
          prec_r(i) = qr_flux_new(i)
          prec_i(i) = qi_flux_new(i)
          prec_s(i) = qs_flux_new(i)
          prec_g(i) = qg_flux_new(i)
          prec_h(i) = qh_flux_new(i)
        ENDDO
        !$ACC END PARALLEL

      END IF

      IF (ldass_lhn) THEN
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
        !$ACC LOOP GANG VECTOR
        DO ii = its, ite
          qrsflux(ii,kte) = qr_flux_new(ii) + qi_flux_new(ii) + qs_flux_new(ii) + qg_flux_new(ii) + qh_flux_new(ii)
        ENDDO
        !$ACC END PARALLEL
      END IF

      !$ACC WAIT
      !$ACC END DATA ! DATA CREATE PRESENT

    END SUBROUTINE clouds_twomoment_implicit
   
   !
   ! sedimentation for explicit solver, i.e., sedimentation is done with an explicit
   ! flux-form semi-lagrangian scheme after the microphysics.
   !
   SUBROUTINE sedimentation_explicit()
    ! D.Rieger: the parameter lfullyexplicit needs to be set false, otherwise the nproma/mpi tests of buildbot are not passed
    LOGICAL, PARAMETER :: lfullyexplicit = .FALSE.
    REAL(wp) :: cmax, rdzmaxdt
    REAL(wp) :: prec3D_tmp(isize,ke)
    INTEGER :: ii, kk

#ifdef _OPENACC
    IF (lprogmelt) THEN
      CALL finish(routine, 'lprogmelt not available on GPU for two-moment microphysics')
    ENDIF
#endif

    !$ACC DATA CREATE(prec3D_tmp)

    cmax = 0.0_wp
    ! Use for sub-stepping of hydrometeors, lfullyexplicit needs to be set to TRUE
    IF (lfullyexplicit) THEN
      rdzmaxdt = maxval(rdz(its:ite,kts:kte)) * dt
      ntsedi_rain = ceiling(rain%vsedi_max*rdzmaxdt)
      ntsedi_graupel = ceiling(graupel%vsedi_max*rdzmaxdt)
      ntsedi_hail = ceiling(hail%vsedi_max*rdzmaxdt)
    ELSE
      ntsedi_rain = 1
      ntsedi_graupel = 1
      ntsedi_hail = 1
    ENDIF

     CALL init(prec_r, lacc=.TRUE., opt_acc_async=.TRUE.)
     CALL init(prec_i, lacc=.TRUE., opt_acc_async=.TRUE.)
     CALL init(prec_s, lacc=.TRUE., opt_acc_async=.TRUE.)
     CALL init(prec_g, lacc=.TRUE., opt_acc_async=.TRUE.)
     CALL init(prec_h, lacc=.TRUE., opt_acc_async=.TRUE.)

     ! The following IF ANY conditions are important only for performance on CPU and don't work with OpenACC
#ifndef _OPENACC
     IF (ANY(qr(its:ite,kts:kte)>0._wp)) THEN
#endif
      IF (ldass_lhn) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kk = kts,kte
          DO ii = its,ite
            prec3D_tmp(ii,kk) = 0.0_wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
      DO ii=1,ntsedi_rain
        CALL sedi_icon_rain(rain,rain_coeffs,qr,qnr,prec_r,prec3D_tmp,qc,rhocorr, &
          & rdz,dt/ntsedi_rain,its,ite,kts,kte,cmax,lacc=.TRUE.)
      END DO
      IF (ldass_lhn) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kk = kts,kte
          DO ii = its,ite
            qrsflux(ii,kk) = qrsflux(ii,kk) + prec3D_tmp(ii,kk)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
#ifndef _OPENACC
     END IF

     IF (ANY(qi(its:ite,kts:kte)>0._wp)) THEN
#endif
      IF (ldass_lhn) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kk = kts,kte
          DO ii = its,ite
            prec3D_tmp(ii,kk) = 0.0_wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
      CALL sedi_icon_sphere(ice,ice_coeffs,qi,qni,prec_i,prec3D_tmp,rhocorr,rdz,dt,its,ite,kts,kte,lacc=.TRUE.)
      IF (ldass_lhn) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kk = kts,kte
          DO ii = its,ite
            qrsflux(ii,kk) = qrsflux(ii,kk) + prec3D_tmp(ii,kk)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
#ifndef _OPENACC
     END IF

     IF (ANY(qs(its:ite,kts:kte)>0._wp)) THEN
#endif
      IF (ldass_lhn) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kk = kts,kte
          DO ii = its,ite
            prec3D_tmp(ii,kk) = 0.0_wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
      CALL sedi_icon_sphere(snow,snow_coeffs,qs,qns,prec_s,prec3D_tmp,rhocorr,rdz,dt,its,ite,kts,kte,lacc=.TRUE.)
      IF (ldass_lhn) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kk = kts,kte
          DO ii = its,ite
            qrsflux(ii,kk) = qrsflux(ii,kk) + prec3D_tmp(ii,kk)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
#ifndef _OPENACC
     END IF

     IF (ANY(qg(its:ite,kts:kte)>0._wp)) THEN
#endif
      IF (ldass_lhn) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kk = kts,kte
          DO ii = its,ite
            prec3D_tmp(ii,kk) = 0.0_wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
       IF (lprogmelt) THEN
         DO ii=1,ntsedi_graupel
           call sedi_icon_sphere_lwf(graupel_lwf,graupel_coeffs,qg,qng,qgl,&
                &                    prec_g,prec3D_tmp,rhocorr,rdz,dt/ntsedi_graupel,its,ite,kts,kte,cmax)
         END DO
       ELSE
         DO ii=1,ntsedi_graupel
           CALL sedi_icon_sphere(graupel,graupel_coeffs,qg,qng,prec_g,prec3D_tmp,rhocorr,rdz,dt/ntsedi_graupel, &
             & its,ite,kts,kte,cmax,lacc=.TRUE.)
         END DO
       END IF
      IF (ldass_lhn) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kk = kts,kte
          DO ii = its,ite
            qrsflux(ii,kk) = qrsflux(ii,kk) + prec3D_tmp(ii,kk)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
#ifndef _OPENACC
     END IF

     IF (ANY(qh(its:ite,kts:kte)>0._wp)) THEN
#endif
      IF (ldass_lhn) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kk = kts,kte
          DO ii = its,ite
            prec3D_tmp(ii,kk) = 0.0_wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
       IF (lprogmelt) THEN
         DO ii=1,ntsedi_hail
           call sedi_icon_sphere_lwf(hail_lwf,hail_coeffs,qh,qnh,qhl,&
                &                    prec_h,prec3D_tmp,rhocorr,rdz,dt/ntsedi_hail,its,ite,kts,kte,cmax)
         END DO
       ELSE
         DO ii=1,ntsedi_hail
           call sedi_icon_sphere(hail,hail_coeffs,qh,qnh,prec_h,prec3D_tmp,rhocorr,rdz,dt/ntsedi_hail, &
             & its,ite,kts,kte,cmax,lacc=.TRUE.)
         END DO
       END IF
      IF (ldass_lhn) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kk = kts,kte
          DO ii = its,ite
            qrsflux(ii,kk) = qrsflux(ii,kk) + prec3D_tmp(ii,kk)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
#ifndef _OPENACC
     END IF
#endif
     
     IF (msg_level > 100)THEN
       WRITE (message_text,'(1X,A,f8.2)') ' sedimentation_explicit  cmax = ',cmax
       CALL message(routine, message_text)
     END IF

    !$ACC WAIT
    !$ACC END DATA

   END SUBROUTINE sedimentation_explicit

    !
    ! check for negative values after microphysics
    !
    SUBROUTINE check_clouds()

      REAL(wp), PARAMETER :: meps = -1e-12

      IF (cloud_type.lt.2000) THEN
         IF (ANY(qh(its:ite,kts:kte)>0._wp)) THEN
            qh(its:ite,kts:kte)  = 0.0_wp
            WRITE (message_text,'(1X,A)') '  qh > 0, after cloud_twomoment for cloud_type < 2000'
            CALL message(routine,TRIM(message_text))
            CALL finish(TRIM(routine),'Error in two_moment_mcrph')
         END IF
         IF (ANY(qnh(its:ite,kts:kte)>0._wp)) THEN
            qnh(its:ite,kts:kte)  = 0.0_wp
            WRITE (message_text,'(1X,A)') '  qnh > 0, after cloud_twomoment for cloud_type < 2000'
            CALL message(routine,TRIM(message_text))
            CALL finish(TRIM(routine),'Error in two_moment_mcrph')
         END IF
      END IF
      IF (msg_level>dbg_level) CALL message(TRIM(routine), " test for negative values")
      IF (MINVAL(cloud%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, cloud%q < 0')
      ENDIF
      IF (MINVAL(rain%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, rain%q < 0')
      ENDIF
      IF (MINVAL(ice%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, ice%q < 0,')
      ENDIF
      IF (MINVAL(snow%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, snow%q < 0')
      ENDIF
      IF (MINVAL(graupel%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, graupel%q < 0')
      ENDIF
      IF (MINVAL(hail%q(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, hail%q < 0')
      ENDIF
      IF (MINVAL(cloud%n) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, cloud%n < 0')
      ENDIF
      IF (MINVAL(rain%n(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, rain%n < 0')
      ENDIF
      IF (MINVAL(ice%n(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, ice%n < 0')
      ENDIF
      IF (MINVAL(snow%n(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, snow%n < 0')
      ENDIF
      IF (MINVAL(graupel%n(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, graupel%n < 0')
      ENDIF
      IF (MINVAL(hail%n(its:ite,kts:kte)) < meps) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, hail%n < 0')
      ENDIF
    END subroutine check_clouds

  END SUBROUTINE two_moment_mcrph

  SUBROUTINE implicit_core(q_val,q_sum,q_impl,vsed_new,vsed_now,flux_new,flux_now,rdzdt,its,ite)

    REAL(wp), DIMENSION(:), INTENT(INOUT) :: &
         &    q_val,q_sum,q_impl,vsed_new,vsed_now,flux_new,flux_now,rdzdt
    INTEGER :: its,ite

    REAL(wp) :: q_star, flux_sum
    INTEGER  :: i

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE(q_star, flux_sum)
    DO i=its,ite

      ! new on r.h.s. is new value from level above
      vsed_new(i) = 0.5 * (vsed_now(i) + vsed_new(i))

      ! qflux_new, nflux_new are the updated flux values from the level above
      ! qflux_now, nflux_now are here the old (current time step) flux values from the level above
      ! In COSMO-Docu  {...} =  flux_(k-1),new + flux_(k-1),start
      flux_sum = flux_new(i) + flux_now(i)

      ! qflux_now, nflux_now are here overwritten with the current level
      flux_now(i) = min(vsed_now(i) * q_val(i),  flux_sum)    ! (loop dependency)
      flux_now(i) = max(flux_now(i),0.0_wp)                   ! maybe not necessary

      ! time integrated value without implicit weight
      q_sum(i)  = q_val(i)  + rdzdt(i) * (flux_sum - flux_now(i))      

      ! implicit weight
      q_impl(i) = 1.0_wp/(1.0_wp + vsed_new(i) * rdzdt(i))

      ! prepare for source term calculation
      q_star    = q_impl(i) * q_sum(i)       
      q_val(i)  = q_star                     ! source/sinks work on star-values
      q_sum(i)  = q_sum(i) - q_star           
    END DO
    !$ACC END PARALLEL
    
  END SUBROUTINE implicit_core
  
  SUBROUTINE implicit_time(q_val,q_sum,q_impl,vsed_new,vsed_now,flux_new,its,ite)

    INTEGER  :: its,ite
    REAL(wp), DIMENSION(:), INTENT(INOUT) :: &
         &    q_val,q_sum,q_impl,vsed_new,vsed_now,flux_new

    INTEGER  :: i
    
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO i=its,ite
      
      ! time integration
      q_val(i) =   MAX( 0.0_wp, q_impl(i)*(q_sum(i) + q_val(i)))    

      ! prepare for next level
      flux_new(i) = q_val(i) * vsed_new(i)     ! flux_(k),new
      vsed_new(i) = vsed_now(i)

    END DO
    !$ACC END PARALLEL
    
  END SUBROUTINE implicit_time

  !===========================================================================================

  SUBROUTINE two_moment_mcrph_init(igscp,N_cn0,z0_nccn,z1e_nccn,N_in0,z0_nin,z1e_nin,msg_level,cfg_2mom)

    INTEGER, INTENT(IN) :: igscp, msg_level

    REAL(wp), OPTIONAL, INTENT(OUT) ::             & ! for CCN and IN in case of gscp=5
         & N_cn0,z0_nccn,z1e_nccn,    &
         & N_in0,z0_nin,z1e_nin

    TYPE(particle)        :: cloud, rain
    TYPE(particle_frozen) :: ice, snow, graupel, hail
    TYPE(particle_lwf)    :: graupel_lwf, hail_lwf

    TYPE(t_cfg_2mom), OPTIONAL, INTENT(in) :: cfg_2mom

    INTEGER        :: unitnr

    ! Transfer the configuration parameters to the 2mom internal type instance:
    IF (PRESENT(cfg_2mom)) THEN
      cfg_params = cfg_2mom
    ELSE
      cfg_params = cfg_2mom_default
    END IF

    IF (msg_level>5) THEN
      CALL message (TRIM(routine), " Initialization of two-moment microphysics scheme") 
      WRITE(message_text,'(A,I5)')   "   inwp_gscp    = ",igscp ; CALL message(TRIM(routine),TRIM(message_text))
      WRITE(message_text,'(A,I5)')   "   i2mom_solver = ",cfg_params%i2mom_solver ; CALL message(TRIM(routine),TRIM(message_text))
      WRITE(message_text,'(A,L5)'  ) "   lconstant_lh = ",lconstant_lh ; CALL message(TRIM(routine),TRIM(message_text))
    END IF

    IF (PRESENT(N_cn0)) THEN
      IF (PRESENT(cfg_2mom)) THEN
        IF (cfg_2mom%ccn_type > 0) THEN
          ccn_type   = cfg_2mom%ccn_type
        ELSE 
          ccn_type   = ccn_type_gscp5
        END IF
      ELSE
        ccn_type   = ccn_type_gscp5
      END IF
      cloud_type = cloud_type_default_gscp5 + 10 * ccn_type
    ELSE
      IF (PRESENT(cfg_2mom)) THEN
        IF (cfg_2mom%ccn_type > 0) THEN
          ccn_type   = cfg_2mom%ccn_type
        ELSE 
          ccn_type   = ccn_type_gscp4
        END IF
      ELSE
        ccn_type   = ccn_type_gscp4
      END IF
      cloud_type = cloud_type_default_gscp4 + 10 * ccn_type
    END IF

    ! .. set the particle types, and calculate some coefficients
    IF (igscp == 7) THEN
       CALL init_2mom_scheme_once(cloud,rain,ice,snow,graupel_lwf,hail_lwf,cloud_type)
    ELSE
       CALL init_2mom_scheme_once(cloud,rain,ice,snow,graupel,hail,cloud_type)
    END IF

    IF (timers_level > 10) CALL timer_start(timer_phys_2mom_dmin_init)
    IF (luse_dmin_wetgrowth_table .OR. lprintout_comp_table_fit) THEN
      unitnr = 11
      IF (msg_level>5) CALL message (TRIM(routine), " Looking for dmin_wetgrowth table file for "//TRIM(graupel%name))
      CALL init_dmin_wg_gr_ltab_equi('dmin_wetgrowth_lookup', graupel, &
           unitnr, 61, ltabdminwgg, msg_level)
      IF (msg_level>5) CALL message (TRIM(routine), " Looking for dmin_wetgrowth table file for "//TRIM(hail%name))
      CALL init_dmin_wg_gr_ltab_equi('dmin_wetgrowth_lookup', hail, &
           unitnr, 61, ltabdminwgh, msg_level)
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
    
    IF (msg_level>dbg_level) CALL message (TRIM(routine), ' finished init_dmin_wetgrowth for '// &
         TRIM(graupel%name)//' and '//TRIM(hail%name))
    !..parameters for CCN and IN are set here. The 3D fields for prognostic CCN are then
    !  initialized in mo_nwp_phy_init.
    IF (timers_level > 10) CALL timer_stop(timer_phys_2mom_dmin_init)

    !..parameters for exponential decrease of N_ccn with height
    !  z0:  up to this height (m) constant unchanged value
    !  z1e: height interval at which N_ccn decreases by factor 1/e above z0_nccn
    
    ccn_coeffs%z0  = 4000.0_wp
    ccn_coeffs%z1e = 2000.0_wp

    ! min updraft speed for Segal&Khain activation
    ccn_coeffs%wcb_min = cfg_params%ccn_wcb_min

    ! characteristics of different kinds of CN
    ! (copied from COSMO 5.0 Segal & Khain nucleation subroutine)

    SELECT CASE(ccn_type)
    CASE(6)
      !... maritime case
      ccn_coeffs%Ncn0 = 100.0e6_wp   ! CN concentration at ground
      ccn_coeffs%Nmin =  35.0e6_wp   ! NOT relevant at the moment
      ccn_coeffs%lsigs = 0.4_wp      ! log(sigma_s)
      ccn_coeffs%R2    = 0.03_wp     ! in mum
      ccn_coeffs%etas  = 0.9_wp      ! soluble fraction
    CASE(7)
      !... intermediate case
      ccn_coeffs%Ncn0 = 250.0e6_wp
      ccn_coeffs%Nmin =  35.0e6_wp
      ccn_coeffs%lsigs = 0.4_wp
      ccn_coeffs%R2    = 0.03_wp       ! in mum
      ccn_coeffs%etas  = 0.8_wp        ! soluble fraction
    CASE(8)
      IF (tune_sbmccn < 1.0_wp) THEN
        !... maritime case
        ccn_coeffs%Ncn0 = 100.0e6_wp   ! CN concentration at ground
        ccn_coeffs%Nmin =  35.0e6_wp   ! NOT relevant at the moment
        ccn_coeffs%lsigs = 0.4_wp      ! log(sigma_s)
        ccn_coeffs%R2    = 0.03_wp     ! in mum
        ccn_coeffs%etas  = 0.9_wp      ! soluble fraction
      ELSE
        !... continental case
        ccn_coeffs%Ncn0 = 1700.0e6_wp
        ccn_coeffs%Nmin =   35.0e6_wp  ! NOT relevant at the moment
        ccn_coeffs%lsigs = 0.2_wp
        ccn_coeffs%R2    = 0.03_wp     ! in mum
        ccn_coeffs%etas  = 0.7_wp      ! soluble fraction
      END IF
    CASE(9)
      !... "polluted" continental
      ccn_coeffs%Ncn0 = 3200.0e6_wp
      ccn_coeffs%Nmin =   35.0e6_wp    ! NOT relevant at the moment
      ccn_coeffs%lsigs = 0.2_wp
      ccn_coeffs%R2    = 0.03_wp       ! in mum
      ccn_coeffs%etas  = 0.7_wp        ! soluble fraction
     CASE(1)
       !... dummy values
       ccn_coeffs%Ncn0  =  200.0e6_wp
       ccn_coeffs%Nmin  =   10.0e6_wp  ! NOT relevant at the moment
       ccn_coeffs%lsigs = 0.0_wp
       ccn_coeffs%R2    = 0.0_wp
       ccn_coeffs%etas  = 0.0_wp
    CASE DEFAULT
       CALL finish(TRIM(routine),'Error in two_moment_mcrph_init: Invalid value for ccn_type')
    END SELECT

    IF (cfg_params%ccn_Ncn0 > -900.0_wp) THEN
      ccn_coeffs%Ncn0 = cfg_params%ccn_Ncn0
    END IF

    IF (PRESENT(N_cn0)) THEN
      z0_nccn  = ccn_coeffs%z0
      z1e_nccn = ccn_coeffs%z1e
      N_cn0    = ccn_coeffs%Ncn0
    END IF
    
    WRITE(message_text,'(A)') "  CN properties:" ; CALL message(TRIM(routine),TRIM(message_text))
    WRITE(message_text,'(A,D10.3)') "    Ncn0 = ",ccn_coeffs%Ncn0 ; CALL message(TRIM(routine),TRIM(message_text))
    WRITE(message_text,'(A,D10.3)') "    z0   = ",ccn_coeffs%z0  ; CALL message(TRIM(routine),TRIM(message_text))
    WRITE(message_text,'(A,D10.3)') "    z1e  = ",ccn_coeffs%z1e ; CALL message(TRIM(routine),TRIM(message_text))

    IF (PRESENT(N_in0)) THEN

       in_coeffs%N0  = 200.0e6_wp ! this is currently just a scaling factor for the PDA scheme
       in_coeffs%z0  = 3000.0_wp
       in_coeffs%z1e = 1000.0_wp

       N_in0   = in_coeffs%N0
       z0_nin  = in_coeffs%z0
       z1e_nin = in_coeffs%z1e

       WRITE(message_text,'(A)') "  IN properties:" ; CALL message(TRIM(routine),TRIM(message_text))
       WRITE(message_text,'(A,D10.3)') "    Ncn0 = ",in_coeffs%N0  ; CALL message(TRIM(routine),TRIM(message_text))
       WRITE(message_text,'(A,D10.3)') "    z0   = ",in_coeffs%z0  ; CALL message(TRIM(routine),TRIM(message_text))
       WRITE(message_text,'(A,D10.3)') "    z1e  = ",in_coeffs%z1e ; CALL message(TRIM(routine),TRIM(message_text))
     END IF
     
    IF (msg_level>5) CALL message (TRIM(routine), " finished two_moment_mcrph_init successfully")
    !$ACC ENTER DATA COPYIN(ccn_coeffs, in_coeffs, cfg_params)
    !$ACC ENTER DATA COPYIN(ltabdminwgg, ltabdminwgh)
    !$ACC ENTER DATA COPYIN(ltabdminwgg%ltable, ltabdminwgg%x1, ltabdminwgg%x2, ltabdminwgg%x3, ltabdminwgg%x4) &
    !$ACC   COPYIN(ltabdminwgh%ltable, ltabdminwgh%x1, ltabdminwgh%x2, ltabdminwgh%x3, ltabdminwgh%x4)

  END SUBROUTINE two_moment_mcrph_init


  ! Subroutine that provides coefficients for the effective radius calculations
  ! consistent with two-moment microphysics
  SUBROUTINE two_mom_reff_coefficients( reff_calc ,return_fct)
    TYPE(t_reff_calc), INTENT(INOUT) ::  reff_calc                   ! Structure with options and coefficiencts
    LOGICAL          , INTENT(INOUT) ::  return_fct                  ! Return code of the subroutine

    ! These are the fundamental hydrometeor particle variables for the two-moment scheme
    TYPE(particle)        , TARGET   :: cloud_hyd, rain_hyd
    TYPE(particle_frozen) , TARGET   :: ice_frz, snow_frz, graupel_frz, hail_frz
    TYPE(particle_lwf)    , TARGET   :: graupel_lwf, hail_lwf

    ! Pointers to the derived types that are actually needed
    CLASS(particle)       , POINTER  :: cloud, rain
    CLASS(particle_frozen), POINTER  :: ice, snow, graupel, hail

    ! Parameters used in the paramaterization of reff (the same for all)
    CLASS(particle)       , POINTER  :: current_hyd
    REAL(wp)                         :: a_geo, b_geo, mu, nu
    REAL(wp)                         :: bf, bf2 
    LOGICAL                          :: monodisperse
        
    ! Check input return_fct
    IF (.NOT. return_fct) THEN
      WRITE (message_text,*) 'Reff: Function two_mom_provide_reff_coefficients entered with previous error'
      CALL message('',message_text)
      RETURN
    END IF

    ! We need to reinitiate the particles because they are deleted after every micro call
    cloud => cloud_hyd
    rain  => rain_hyd
    ice   => ice_frz
    snow  => snow_frz
    IF (reff_calc%microph_param == 7 ) THEN  ! Vivek Param. Frozen+Liquid
       graupel => graupel_lwf                ! gscp=7
       hail    => hail_lwf
    ELSE
       graupel => graupel_frz                ! gscp=4,5,6
       hail    => hail_frz
    END IF
    ! .. set the particle types, but no calculations
    CALL init_2mom_scheme(cloud,rain,ice,snow,graupel,hail)
   
    SELECT CASE ( reff_calc%hydrometeor )    ! Select Hydrometeor
    CASE (0)                                 ! Cloud water
      current_hyd => cloud
    CASE (1)  
      current_hyd => ice
    CASE (2)  
      current_hyd => rain
    CASE (3)  
      current_hyd => snow
    CASE (4)  
      current_hyd => graupel
    CASE (5)  
      current_hyd => hail
    END SELECT

    ! Extract properties of hydrometeor
    a_geo           = current_hyd%a_geo
    b_geo           = current_hyd%b_geo
    mu              = current_hyd%mu
    nu              = current_hyd%nu
    reff_calc%x_min = current_hyd%x_min
    reff_calc%x_max = current_hyd%x_max

    ! All DSD are polydisperse
    monodisperse    = .false.

    ! Overwrite monodisperse/polydisperse according to options
    SELECT CASE (reff_calc%dsd_type)
    CASE (1)
      monodisperse  = .true.
    CASE (2)
      monodisperse  = .false.
    END SELECT

    IF ( reff_calc%dsd_type == 2) THEN       ! Overwrite mu and nu coefficients
      mu            = reff_calc%mu
      nu            = reff_calc%nu
    END IF

    SELECT CASE ( reff_calc%reff_param )     ! Select Parameterization
    CASE(0)                                  ! Spheroids  Dge = c1 * x**[c2], which x = mean mass
      ! First calculate monodisperse
      reff_calc%reff_coeff(1)   = a_geo
      reff_calc%reff_coeff(2)   = b_geo

      ! Broadening for not monodisperse
      IF ( .NOT. monodisperse ) THEN 
        bf =  GAMMA( (3.0_wp * b_geo + nu + 1.0_wp)/ mu) / GAMMA( (2.0_wp * b_geo + nu + 1.0_wp)/ mu) * &
          & ( GAMMA( (nu + 1.0_wp)/ mu) / GAMMA( (nu + 2.0_wp)/ mu) )**b_geo

        reff_calc%reff_coeff(1) = reff_calc%reff_coeff(1)*bf        
      END IF      

    CASE (1)                                 ! Fu Random Hexagonal needles:  Dge = 1/(c1 * x**[c2] + c3 * x**[c4])
                                             ! Parameterization based on Fu, 1996; Fu et al., 1998; Fu ,2007
      ! First calculate monodisperse
      reff_calc%reff_coeff(1)   = SQRT( 3.0_wp *SQRT(3.0_wp) * rho_ice * a_geo / 8.0_wp )
      reff_calc%reff_coeff(2)   = (b_geo - 1.0_wp)/2.0_wp 
      reff_calc%reff_coeff(3)   = SQRT(3.0_wp)/4.0_wp/a_geo
      reff_calc%reff_coeff(4)   = -b_geo

      ! Broadening for not monodisperse. Generalized gamma distribution
      IF ( .NOT. monodisperse ) THEN 
        bf  =  GAMMA( ( b_geo + 2.0_wp * nu + 3.0_wp)/ mu/2.0_wp ) / GAMMA( (nu + 2.0_wp)/ mu) * &
           & ( GAMMA( (nu + 1.0_wp)/ mu) / GAMMA( (nu + 2.0_wp)/ mu) )**( (b_geo-1.0_wp)/2.0_wp)

        bf2 =  GAMMA( (-b_geo + nu + 2.0_wp)/ mu ) / GAMMA( (nu + 2.0_wp)/ mu) * &
           & ( GAMMA( (nu + 1.0_wp)/ mu) / GAMMA( (nu + 2.0_wp)/ mu) )**( -b_geo)

        reff_calc%reff_coeff(1) = reff_calc%reff_coeff(1)*bf
        reff_calc%reff_coeff(3) = reff_calc%reff_coeff(3)*bf2
      END IF

    END SELECT

  END SUBROUTINE two_mom_reff_coefficients

END MODULE mo_2mom_mcrph_driver
