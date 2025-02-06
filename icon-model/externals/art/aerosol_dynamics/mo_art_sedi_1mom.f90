!
! mo_art_sedi_1mom
! This module provides the routine for calculating the sedimentation velocity
! of one-moment aerosol
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

MODULE mo_art_sedi_1mom
! ICON
  USE mo_kind,                          ONLY: wp,vp
  USE mo_fortran_tools,                 ONLY: set_acc_host_or_device, t_ptr_tracer
  USE mo_physical_constants,            ONLY: grav
! ART
  USE mo_art_data,                      ONLY: t_art_data
  USE mo_art_config,                    ONLY: t_art_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC   :: art_sedi_1mom

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_sedi_1mom(temp,pres,rho,rho_ic,wgtfac_c,wgtfacq_c,diameter_particle,rho_particle,turb_tracer, &
  &                      istart,iend,nlev,jb,itracer,art_config,p_art_data,p_mflx_contra_vsed,lacc)
!<
! SUBROUTINE art_sedi_1mom
! This subroutine calculates the sedimentation velocity of monodisperse tracer
! Based on: Rieger et al. (2015) - ICON-ART 1.0-a new online-coupled model system 
!                                  from the global to regional scale
! Part of Module: mo_art_sedi_1mom
! Author: Kristina Lundgren, KIT
! Initial Release: 2012-01-30
! Modifications:
! 2016-11-04: Daniel Rieger, KIT
! - Adapted to unified ART-physics structure
!>
  REAL(wp), INTENT(in)          :: &
    &  temp(:,:),                  & !< temperature (K)
    &  pres(:,:),                  & !< pressure (Pa)
    &  rho(:,:),                   & !< air densitry (kg m-3)
    &  rho_ic(:,:),                & !< air densitry interpolated to half levels (kg m-3)
    &  diameter_particle,          & !< particle diameter (m)
    &  rho_particle                  !< particle density (kg m-3)
  REAL(vp), INTENT(in)          :: &
    &  wgtfac_c(:,:),              & !< weighting factor for interpolation from full to half levels
    &  wgtfacq_c(:,:)                !< weighting factor for quadratic interpolation to surface
  TYPE(t_ptr_tracer),INTENT(in) :: &
    &  turb_tracer(:)                !< List of turbulent tracers
  INTEGER,INTENT(in)            :: &
    &  istart, iend,               & !< Start and end index of nproma loop
    &  nlev,                       & !< number of vertical levels
    &  jb,                         & !< block index
    &  itracer                       !< Index in tracer container
  TYPE(t_art_config)            :: &
    &  art_config                    !< ART configuration state
  TYPE(t_art_data),INTENT(inout):: &
    &  p_art_data                    !< ART data container
  REAL(wp),INTENT(INOUT)        :: &
    &  p_mflx_contra_vsed(:,:)       !< dim(nproma,nlevp1) inout is necessary as it is allocated as a pointer
  LOGICAL, OPTIONAL, INTENT(IN) :: &
    &  lacc
!Local variables
  REAL(wp),ALLOCATABLE          :: &
    &  vsed(:,:),                  & !< Sedimentation velocity at full levels
    &  vsed_ifc(:,:)                 !< Sedimentation velocity at half levels
  INTEGER                       :: &
    &  jk,jc,iturb,                & !< Loop indices
    &  nlevp1,                     & !< nlev+1
    &  jsp                           !< Index of tracer in container
  REAL(wp)                      :: &
    &  dver,                       & !< 
    &  rep,                        & !< particle Reynolds number
    &  cd,                         & !< drag coefficient after Fuchs(1964) and Friedlander (1977)
    &  nue,                        & !< kinematic viscosity (m2 s-1)
    &  vst,                        & !< STOKES SETTLING VELOCITY
    &  dp,                         & !< particle diameter in microns
    &  p0,                         & !< reference pressure for mean free path calculations
    &  t0,                         & !< reference temperature
    &  lambda0,                    & !< reference mean free path at p0,t0
    &  lambda,                     & !< mean free path
    &  kn,                         & !< knudsen number
    &  cc,                         & !< Cunningham slip correction factor
    &  rho_nlevp1                    !< store value of extrapolated rho at lowest half layer
  
  LOGICAL :: lzacc             ! OpenACC flag
  CALL set_acc_host_or_device(lzacc, lacc)

  nlevp1  = nlev + 1
  ALLOCATE(vsed(istart:iend,nlev),vsed_ifc(istart:iend,nlevp1))

  !$ACC DATA CREATE(vsed, vsed_ifc) PRESENT(pres, p_art_data, p_art_data%turb_fields, p_art_data%turb_fields%vdep) &
  !$ACC   PRESENT(p_mflx_contra_vsed, rho, rho_ic, temp, wgtfacq_c, wgtfac_c) IF(lzacc)

  dp      = diameter_particle * 1.0e6_wp  ! diameter in microns
  p0      = 101300.25_wp
  t0      = 293.15_wp
  lambda0 = 6.6E-8_wp                ! mean free path at t0, p0

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(cc, cd, dver, kn, lambda, nue, rep, vst)
  DO jk = 1, nlev
    DO jc = istart, iend
      dver=rho_particle/rho(jc,jk)
      nue = (1.458_wp*temp(jc,jk) * SQRT(temp(jc,jk)))     &
        & / ((temp(jc,jk)+110.4_wp)*rho(jc,jk)*(1.E06_wp))

      ! mean free path at grid point
      lambda = lambda0 * p0 * temp(jc,jk)  / (t0 * pres(jc,jk) )
      ! Knudsen number
      kn     = 2.0_wp * lambda / dp
      ! Cunningham slip correction (coefficients from seinfeld and pandis 2nd edition)
      cc     = 1.0_wp + kn * (1.257_wp + 0.4_wp*EXP( -1.1_wp / kn )) 
      vst    = (cc * grav*dver*dp*dp*(1.e-12_wp))/(18._wp*nue)

      ! Reynolds number
      rep = (vst*dp*(1.e-6_wp))/nue
      ! Drag coefficient
      cd = (24._wp/rep) * (1._wp+0.158_wp*(rep**(2._wp/3._wp)))

      vsed(jc,jk) = ((cc * 4._wp*rho_particle*grav*dp*(1.E-6_wp))/(3._wp*rho(jc,jk)*cd))
      vsed(jc,jk) = SQRT(vsed(jc,jk))
    ENDDO
  ENDDO
  !$ACC END PARALLEL

  !interpolate vsed to half levels and calculate mass flux here. Then use the PPM method for the tracer flux calculation.
  ! interpolate sedimentation velocity from full levels to half levels:
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP SEQ
  DO jk = 2, nlev  !at which level should we start??
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = istart,iend
      vsed_ifc(jc,jk)           = wgtfac_c(jc,jk)*vsed(jc,jk)           &
        &                       + (1._wp-wgtfac_c(jc,jk))*vsed(jc,jk-1)
      p_mflx_contra_vsed(jc,jk) = - vsed_ifc(jc,jk)*rho_ic(jc,jk)
    ENDDO  !jc
  ENDDO !jk

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(rho_nlevp1)
  DO jc = istart, iend
    vsed_ifc(jc,     1)           = vsed (jc,1)
    vsed_ifc(jc,nlevp1)           = vsed (jc,nlev)
    p_mflx_contra_vsed(jc,1)      = - vsed_ifc(jc,1)*rho_ic(jc,1)
    rho_nlevp1                    = wgtfacq_c(jc,1)*rho(jc,nlev)     &
      &                           + wgtfacq_c(jc,2)*rho(jc,nlev-1)   &
      &                           + wgtfacq_c(jc,3)*rho(jc,nlev-2)
    p_mflx_contra_vsed(jc,nlevp1) = - vsed_ifc(jc,nlevp1)*rho_nlevp1
  ENDDO
  !$ACC END PARALLEL

  ! drieg: quick solution, set vdep to vsed
  DO iturb = 1, art_config%nturb_tracer
    IF (turb_tracer(iturb)%idx_tracer == itracer) THEN
      ! Get the index of the turbulent tracer field
      jsp = iturb
    ENDIF
  ENDDO

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO jc = istart, iend
    p_art_data%turb_fields%vdep(jc,jb,jsp) = vsed(jc,nlev)
  END DO
  !$ACC END PARALLEL

  !$ACC WAIT
  !$ACC END DATA

  DEALLOCATE(vsed,vsed_ifc)


END SUBROUTINE art_sedi_1mom
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_sedi_1mom
