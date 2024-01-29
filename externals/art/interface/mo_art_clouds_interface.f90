!
! Provides interface to ART-routines dealing with aerosol-cloud-interactions
!
! This module provides an interface to a version of the Seifert and Beheng
! two-moment cloud microphysics scheme which incorporates prognostic aerosol
! as calculated by the ART routines.
! The interface is written in such a way, that ICON will compile and run
! properly, even if the ART-routines are not available at compile time.
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

MODULE mo_art_clouds_interface

  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_run_config,                    ONLY: lart
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_cldInt
  USE mo_2mom_mcrph_config,             ONLY: t_cfg_2mom
  USE mo_2mom_mcrph_processes,          ONLY: cfg_2mom_default, cfg_params
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_2mom_driver,               ONLY: art_2mom_mcrph,               &
                                          &   art_2mom_mcrph_init
  USE mo_art_prepare_aerosol,           ONLY: art_prepare_dust_KL06

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_clouds_interface_2mom
  PUBLIC  :: art_clouds_interface_2mom_init

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clouds_interface_2mom(isize, ke, jg, jb, is, ie, ks, dt, &
                                   & dz, rho, pres, tke, p_trac, tk,    &
                                   & w, prec_r, prec_i, prec_s,         &
                                   & prec_g, prec_h, tkvh, l_cv)
  !! Interface for ART: Aerosol-Cloud-Interactions
  !! @par Revision History
  !! Initial revision by Daniel Rieger, KIT (2014-11-10)
  ! Setup variables (Grid, timestep, looping)
  INTEGER,            INTENT (in) :: &
    &  isize, ke,                    & !< grid sizes
    &  jg, jb,                       & !< domain index (p_patch%id)
    &  is, ie, ks                      !< start/end indices
  REAL(wp), INTENT(in)            :: &
    &  dt                              !< time step
  ! Dynamical core variables
  REAL(wp), INTENT(in), TARGET    :: &
    &  dz(:,:),                      & !< Vertical layer thickness
    &  rho(:,:),                     & !< Density
    &  pres(:,:),                    & !< Pressure
    &  tke(:,:),                     & !< Turbulent kinetic energy
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
    
  IF (lart) THEN
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_cldInt)

    
    ! ----------------------------------
    ! --- Call of the coupled ART-twomoment microphysics
    ! ----------------------------------
    IF (art_config(jg)%iart_aci_cold == 6 .OR. art_config(jg)%iart_aci_cold == 7) THEN
      CALL art_prepare_dust_KL06(jg,jb,is,ie,ks,ke,rho,p_trac)
    ENDIF
    CALL art_2mom_mcrph(isize, ke, jg, jb, is, ie, ks, dt,           &
                        & dz, rho, pres, tke, p_trac(:,:,:), tk,    &
                        & w, prec_r, prec_i, prec_s,         &
                        & prec_g, prec_h, tkvh, l_cv)

    IF (timers_level > 3) CALL timer_stop(timer_art_cldInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  ELSE
    call finish('mo_art_clouds_interface:art_clouds_interface_2mom', &
         &      'Two moment micophysics with ART aerosol chosen (inwp_gscp=6), but lart=.FALSE.')
  ENDIF !lart

END SUBROUTINE art_clouds_interface_2mom
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clouds_interface_2mom_init(msg_level,cfg_2mom)
  !! Interface for ART: Aerosol-Cloud-Interactions Initialization
  !! @par Revision History
  !! Initial revision by Daniel Rieger, KIT (2014-11-10)
  INTEGER, INTENT(IN) :: &
    &  msg_level           !< message level

  TYPE(t_cfg_2mom), OPTIONAL, INTENT(in) :: cfg_2mom

  ! Transfer the configuration parameters to the 2mom internal type instance:
  IF (PRESENT(cfg_2mom)) THEN
    cfg_params = cfg_2mom
  ELSE
    cfg_params = cfg_2mom_default
  END IF
  
  IF (lart) THEN
    CALL art_2mom_mcrph_init(msg_level)
  ELSE
    call finish('mo_art_clouds_interface:art_clouds_interface_2mom_init', &
         &      'Two moment micophysics with ART aerosol chosen (inwp_gscp=6), but lart=.FALSE.')
  ENDIF !lart

END SUBROUTINE art_clouds_interface_2mom_init
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_clouds_interface
