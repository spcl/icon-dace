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

!------------------------------------------------------------------------------
!
! Description:
!
!   The coupling to radiation through the reff calculation and the
!   radar calculation make use of the two-moment formulations.
!   Hence, the number and mass moments are calculated from SBM and
!   are then passed to these subroutines. This introduces a spurious
!   dependency of SBM results on 2mom parameters (PSD, particle geometry).
!   A more consistent coupling that makes full use of the SBM particle
!   distributions and calculating reff directly from SBM is planned.

MODULE mo_sbm_driver

!==============================================================================
!
! Declarations:
!
! Modules used:
!------------------------------------------------------------------------------
! Microphysical constants and variables
!------------------------------------------------------------------------------

USE mo_kind,                 ONLY: wp

USE mo_exception,            ONLY: finish, message, message_text
USE mo_run_config,           ONLY: iqb_i, iqb_e

USE mo_2mom_mcrph_driver,    ONLY: two_moment_mcrph
USE mo_sbm_util,             ONLY: p_ff8i01,p_ff8i33

USE mo_sbm_main,             ONLY: warm_sbm
!==============================================================================

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*), PARAMETER :: routine = 'mo_sbm_driver'
  INTEGER,          PARAMETER :: dbg_level = 25                   ! level for debug prints

CONTAINS
  
  !==============================================================================
  !
  ! SBM warm phase microphysics
  !
  ! qnx in SBM is in units of 1/kg
  ! qx  in SBM is in units of kg/kg
  ! 
  !==============================================================================
  SUBROUTINE sbm(            &              ! used to be two_moment_mcrph
                       isize,             & ! in: array size
                       ke,                & ! in: end level/array size
                       is,                & ! in: start index, optional
                       ie,                & ! in: end index, optional
                       ks,                & ! in: start level ! needed for 2M
                       dt,                & ! in: time step
                       dz,                & ! in: vertical layer thickness
                       hhl,               & ! in: height of half levels
                       rho,               & ! in: density
                       pres,              & ! in: pressure
                       tke,               & ! in:  turbulent kinetic energy (on half levels, size nlev+1)
                       qv,                & ! inout: specific humidity (kg/kg atm_dyn_iconam/mo_nonhydro_state.f90)
                       qc, qnc,           & ! inout: cloud water (kg/kg, 1/kg atm_dyn_iconam/mo_nonhydro_state.f90)
                       qr, qnr,           & ! inout: rain (kg/kg, 1/kg atm_dyn_iconam/mo_nonhydro_state.f90)
                       qi, qni,           & ! inout: ice (kg/kg, 1/kg atm_dyn_iconam/mo_nonhydro_state.f90)
                       qs, qns,           & ! inout: snow (kg/kg, 1/kg atm_dyn_iconam/mo_nonhydro_state.f90)
                       qg, qng,           & ! inout: graupel (kg/kg, 1/kg atm_dyn_iconam/mo_nonhydro_state.f90)
                       qh, qnh,           & ! inout: hail (kg/kg, 1/kg atm_dyn_iconam/mo_nonhydro_state.f90)
!                      nccn,              & ! inout: ccn (1/kg atm_dyn_iconam/mo_nonhydro_state.f90)
!                      ninpot,            & ! inout: potential ice nuclei
                       ninact,            & ! inout: activated ice nuclei
                       tk,                & ! inout: temp 
                       w,                 & ! inout: w
                       prec_r,            & ! inout: precip rate rain
                       prec_i,            & ! inout: precip rate ice
                       prec_s,            & ! inout: precip rate snow
                       prec_g,            & ! inout: precip rate graupel
                       prec_h,            & ! inout: precip rate hail
                       qrsflux,           & ! inout: 3D total precipitation rate
!                      dtemp,             & ! inout: opt. temp increment
                       msg_level,         & ! in: msg_level
!                      l_cv,              & ! in: switch for cv/cp
                       ithermo_water,     & ! in: thermodynamic option - needed for 2M
                       qbin,              &
                       qv_before_satad,   & 
                       tk_before_satad,   & 
                       qv_old,            &
                       temp_old,          & 
!                      u,                 & ! in: u
!                      v,                 & ! in: v
                       exner,             & ! in: exner
!                      fr_land,           & ! in: fr_land
                       lsbm_warm_full)       !0-Piggy Backing with 2M, 1-full warm SBM

    ! Declare variables in argument list
    
    INTEGER,            INTENT (IN)  :: isize, ke    ! grid sizes
    INTEGER,  OPTIONAL, INTENT (IN)  :: is, ie, ks   ! start/end indices

    REAL(wp), INTENT (IN)            :: dt           ! time step

    ! Dynamical core variables
    REAL(wp), DIMENSION(:,:), INTENT(IN), TARGET :: dz, rho, pres, w, exner

    ! Optional Dynamical core variables
    REAL(wp), DIMENSION(:,:), INTENT(IN), POINTER :: tke

    REAL(wp), DIMENSION(:,:), INTENT(IN), TARGET :: hhl

    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET :: tk
!   REAL(wp), DIMENSION(:),   INTENT(IN), TARGET :: fr_land
    REAL(wp), DIMENSION(:,:), INTENT(IN), TARGET :: qv_before_satad, tk_before_satad, temp_old, qv_old !, theta_v_old, exner_old
!   REAL(KIND=wp), DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: extra_3d
    ! Microphysics variables
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) , TARGET :: &
         qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, ninact
    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) , TARGET :: &
         qbin

!   REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
!        &               qgl, qhl

!   REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
!        &               nccn, ninpot

    ! Precip rates, vertical profiles
    REAL(wp), DIMENSION(:), INTENT (INOUT) :: &
         &               prec_r, prec_i, prec_s, prec_g, prec_h
    REAL(wp), DIMENSION(:,:), INTENT (INOUT) :: qrsflux
    
!   REAL(wp), OPTIONAL, INTENT (INOUT)  :: dtemp(:,:)

    INTEGER,  INTENT (IN)             :: msg_level
    LOGICAL,  OPTIONAL, INTENT (IN)   :: lsbm_warm_full
    INTEGER,  OPTIONAL,  INTENT (IN)  :: ithermo_water

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::        &
         &  theta,         & ! potential temperature
         &  lh_rate,       &
         &  ce_rate,       &
         &  cldnucl_rate,  &
         &  qna_nucl,      &
         &  nccn2,         &
         &  theta_old,     &
         &  diag_satur_ba,diag_satur_aa,diag_satur_am,diag_supsat_out, &
         &  reff,reffc,reffr, &
         &  qv_sbm,qc_sbm,qr_sbm,qnc_sbm,qnr_sbm
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: prec_r_sbm

    INTEGER  :: its,ite,kts,kte
    INTEGER  :: ii,kk    !,bin,n_chem

!   CHARACTER(len=*), PARAMETER :: routine = 'mo_sbm_driver' !used to be: 'mo_2mom_mcrph_driver'

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
    
    kts = 1
    kte = ke
    
    ALLOCATE(theta_old(isize,ke))
    ALLOCATE(theta(isize,ke))
    ALLOCATE(nccn2(isize,ke))
    ALLOCATE(qna_nucl(isize,ke))
    ALLOCATE(lh_rate(isize,ke))
    ALLOCATE(ce_rate(isize,ke))
    ALLOCATE(cldnucl_rate(isize,ke))
    ALLOCATE(diag_satur_ba(isize,ke))
    ALLOCATE(diag_satur_aa(isize,ke))
    ALLOCATE(diag_satur_am(isize,ke))
    ALLOCATE(diag_supsat_out(isize,ke))
    ALLOCATE(reff(isize,ke))
    ALLOCATE(reffc(isize,ke))
    ALLOCATE(reffr(isize,ke))
    ALLOCATE(qv_sbm(isize,ke))
    ALLOCATE(qc_sbm(isize,ke))
    ALLOCATE(qr_sbm(isize,ke))
    ALLOCATE(qnc_sbm(isize,ke))
    ALLOCATE(qnr_sbm(isize,ke))
    ALLOCATE(prec_r_sbm(isize))

    DO ii = its, ite
      DO kk = kts, kte

        nccn2(ii,kk) = 0.0_wp   
        qna_nucl(ii,kk) = 0.0_wp 
        lh_rate(ii,kk) = 0.0_wp   
        ce_rate(ii,kk) = 0.0_wp    
        cldnucl_rate(ii,kk) = 0.0_wp

        qv_sbm(ii,kk)=qv_before_satad(ii,kk)
        theta(ii,kk) = tk_before_satad(ii,kk)/exner(ii,kk)

        qc_sbm(ii,kk) = 0.0_wp
        qr_sbm(ii,kk) = 0.0_wp
        qnc_sbm(ii,kk) = 0.0_wp
        qnr_sbm(ii,kk) = 0.0_wp
             
        diag_satur_ba(ii,kk)=0.0_wp
        diag_satur_aa(ii,kk)=0.0_wp
        diag_satur_am(ii,kk)=0.0_wp
        diag_supsat_out(ii,kk)=0.0_wp
      END DO
      prec_r_sbm(ii) = 0.0_wp
    END DO

    CALL WARM_SBM(dt=dt                &!in:    dt
                 ,dz8w=dz                 &!in:    vertical layer thickness
!                ,xland=fr_land           &!in:    land fraction 
                 ,rho_phy=rho             &!in:    density
                 ,p_phy=pres              &!in:    pressure
                 ,pi_phy=exner            &!in:    exner
                 ,w=w                     &!in:    velocities 
                 ,qv_old=qv_old           &
                 ,th_phy=theta            &!inout: theta. Check how to update prognostic theta_v
                 ,qv=qv_sbm               &
                 ,chem_new=qbin           &!inout: 99 mass bins
                 ,rainncv=prec_r_sbm      &!inout: 1 time step precipitation (mm/sec).    
                 ,qc=qc_sbm               &!inout: cloud water: input: 0 
                 ,qr=qr_sbm               &!inout: rain water:  input: 0
                 ,qnc=qnc_sbm             &!inout: cloud water concentration:input: 0
                 ,qnr=qnr_sbm             &!inout: rain water concentration: input: 0
                 ,qna=nccn2               &!inout: ccn concentration:   input: 0
                 ,qna_nucl=qna_nucl       &!inout: nucleated ccn concentration: input:0
                 ,lh_rate=lh_rate         &!inout: rate 1:      input: 0, output can go further to the model
                 ,ce_rate=ce_rate         &!inout: rate 2:      input: 0, output can go further to the model
                 ,cldnucl_rate=cldnucl_rate &!inout: rate 3:    input: 0, output can go further to the model
!                ,kde=kte &!in:    subdomain indeces
!                ,kme=kte &!in:    subdomain indeces
                 ,its=its,ite=ite, kts=kts,kte=kte &!in:    subdomain indeces
                 ,diag_satur_ba=diag_satur_ba          &
                 ,diag_satur_aa=diag_satur_aa          &
                 ,diag_satur_am=diag_satur_am          &
                 ,temp_old=temp_old                    &
                 ,temp_new=tk_before_satad             &
                 ,reff=reff    &
                 ,reffc=reffc  &
                 ,reffr=reffr  &
                 ,diag_supsat_out=diag_supsat_out) 

    IF (lsbm_warm_full) then ! if sbm only is used
      DO ii = its, ite
        DO kk = kts, kte
          qv(ii,kk)=qv_sbm(ii,kk)     !kg/kg
          qc(ii,kk)=qc_sbm(ii,kk)     !kg/kg
          qr(ii,kk)=qr_sbm(ii,kk)     !kg/kg
          qnc(ii,kk)=qnc_sbm(ii,kk)   !1/kg
          qnr(ii,kk)=qnr_sbm(ii,kk)   !1/kg
          tk(ii,kk)=theta(ii,kk)*exner(ii,kk)
          qi(ii,kk)=0.0_wp  ! inout: ice ok
          qni(ii,kk)=0.0_wp ! inout: ice ok
          qs(ii,kk)=0.0_wp  ! inout: snow ok
          qns(ii,kk)=0.0_wp ! inout: snow ok
          qg(ii,kk)=0.0_wp  ! inout: graupel ok
          qng(ii,kk)=0.0_wp ! inout: graupel ok
          qh(ii,kk)=0.0_wp  ! inout: hail ok
          qnh(ii,kk)=0.0_wp ! inout: hail ok
        END DO
        prec_r(ii)=prec_r_sbm(ii)
        prec_i(ii)=0.0_wp
        prec_s(ii)=0.0_wp
        prec_g(ii)=0.0_wp
        prec_h(ii)=0.0_wp
      END DO
    ELSE !Piggybacking (2mom-->dynamics, sbm-->output only)
      CALL two_moment_mcrph(                       &
                       isize  = isize, &!nproma,                &!in: array size
                       ke     = ke, &!nlev,                  &!in: end level/array size
                       is     = is, &!i_startidx,            &!in: start index
                       ie     = ie, &!i_endidx,              &!in: end index
                       ks     = ks, &!kstart_moist(jg),      &!in: start level
                       dt     = dt, &!tcall_gscp_jg ,        &!in: time step
                       dz     = dz, &!p_metrics%ddqz_z_full(:,:,jb),  &!in: vertical layer thickness
                       hhl    = hhl, &!p_metrics%z_ifc(:,:,jb),        &!in: height of half levels
                       rho    = rho, &!p_prog%rho(:,:,jb  )       ,    &!in:  density
                       pres   = pres, &!p_diag%pres(:,:,jb  )      ,    &!in:  pressure
                       tke    = tke, &!ptr_tke_loc, &!in:  turbulent kinetic energy (on half levels, size nlev+1)
                       qv     = qv, &!ptr_tracer (:,:,jb,iqv), &!inout:sp humidity
                       qc     = qc, &!ptr_tracer (:,:,jb,iqc), &!inout:cloud water
                       qnc    = qnc, &!ptr_tracer (:,:,jb,iqnc),&!inout: cloud droplet number
                       qr     = qr, &!ptr_tracer (:,:,jb,iqr), &!inout:rain
                       qnr    = qnr, &!ptr_tracer (:,:,jb,iqnr),&!inout:rain droplet number
                       qi     = qi, &!ptr_tracer (:,:,jb,iqi), &!inout: ice
                       qni    = qni, &!ptr_tracer (:,:,jb,iqni),&!inout: cloud ice number
                       qs     = qs, &!ptr_tracer (:,:,jb,iqs), &!inout: snow
                       qns    = qns, &!ptr_tracer (:,:,jb,iqns),&!inout: snow number
                       qg     = qg, &!ptr_tracer (:,:,jb,iqg), &!inout: graupel
                       qng    = qng, &!ptr_tracer (:,:,jb,iqng),&!inout: graupel number
                       qh     = qh, &!ptr_tracer (:,:,jb,iqh), &!inout: hail
                       qnh    = qnh, &!ptr_tracer (:,:,jb,iqnh),&!inout: hail number
                       ninact = ninact, &!ptr_tracer (:,:,jb,ininact), &!inout: IN number
                       tk     = tk, &!p_diag%temp(:,:,jb),            &!inout: temp 
                       w      = w, &!p_prog%w(:,:,jb),               &!inout: w
                       prec_r = prec_r, &!prm_diag%rain_gsp_rate (:,jb),  &!inout precp rate rain
                       prec_i = prec_i, &!prm_diag%ice_gsp_rate (:,jb),   &!inout precp rate ice
                       prec_s = prec_s, &!prm_diag%snow_gsp_rate (:,jb),  &!inout precp rate snow
                       prec_g = prec_g, &!prm_diag%graupel_gsp_rate (:,jb),&!inout precp rate graupel
                       prec_h = prec_h, &!prm_diag%hail_gsp_rate (:,jb),  &!inout precp rate hail
                       qrsflux= qrsflux, &!prm_diag%qrs_flux(:,:,jb),      & !inout: 3D precipitation flux for LHN
                       msg_level = msg_level,                   &
                       & l_cv=.TRUE.,                           &
                       & ithermo_water=ithermo_water) !atm_phy_nwp_config(jg)%ithermo_water ) !< in: latent heat choice

    END IF

  END SUBROUTINE sbm                 

END MODULE mo_sbm_driver ! used to be mo_2mom_mcrph_driver
