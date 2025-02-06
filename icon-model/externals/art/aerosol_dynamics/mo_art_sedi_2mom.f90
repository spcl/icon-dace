!
! mo_art_sedi_2mom
! This module provides the calculation of the sedimentation
! velocities of log-normally distributed two-moment aerosol
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

MODULE mo_art_sedi_2mom
! ICON
  USE mo_kind,                          ONLY: wp, vp
  USE mo_physical_constants,            ONLY: grav ! gravitational acceleration [m s-2]
! ART
  USE mo_art_clipping,                  ONLY: art_clip_gt

IMPLICIT NONE
    
  PRIVATE
    
  PUBLIC :: art_calc_v_sed
  PUBLIC :: art_calc_sed_flx

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_v_sed(dyn_visc,rho_aero,diam_aero, exp_aero, knudsen_aero, &
                         & istart,iend,nlev,vsed0,vsed3)
!<
! SUBROUTINE art_calc_v_sed
! This subroutine calculates sedimentation velocities for 2-Moment aerosol
! Based on: N. Riemer - Numerische Simulationen zur Wirkung des Aerosols
!              (2002)   auf die troposphaerische Chemie und die Sichtweite
!                       Dissertation
!                       Fakultaet fuer Physik, Universitaet Karlsruhe
! Part of Module: mo_art_sedi_2mom
! Author: Daniel Rieger, KIT
! Initial Release: 2014-06-25
! Modifications:
! 2014-11-24: Daniel Rieger, KIT
! - Adapted to unified ART-physics structure
!>
  REAL(wp), INTENT(in)   :: &
    &  dyn_visc(:,:),       & !< Dynamic viscosity
    &  rho_aero(:,:),       & !< Aerosol density
    &  diam_aero(:,:),      & !< Aerosol median diameter (with respect to number conc.)
    &  exp_aero,            & !< Aerosol exponent
    &  knudsen_aero(:,:)      !< Aerosol knudsen number
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end of inner loop (nproma)
    &  nlev                   !< Number of vertical levels
  REAL(wp), INTENT(inout):: &
    &  vsed0(:,:),          & !< sedimentation velocity of zeroth moment
    &  vsed3(:,:)             !< sedimentation velocity of third moment
! Local variables
  REAL(wp)               :: &
    &  vsed_coeff             !< Eq. 3.69 factor for sedimentation velocity

  INTEGER ::                &
    &  jc, jk                 !< loop indizes
  
  !-----------------------------------------------------------------------------------------
  !--   Start Routine ----------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------
  
  vsed0(:,:) = 0.0_wp
  vsed3(:,:) = 0.0_wp
  
  DO jk=1,nlev
!NEC$ ivdep
    DO jc=istart,iend
      IF (diam_aero(jc,jk) > 0.0_wp) THEN
        ! ----------------------------------
        ! --- Calculate the sedimentation velocity
        ! ----------------------------------
  
        ! Eq. 3.69 factor
        vsed_coeff   = grav / ( 18.0_wp * dyn_visc(jc,jk) )    &
          &          * rho_aero(jc,jk) * (diam_aero(jc,jk))**2
        
        ! Eq. 3.69 for k=0
        vsed0(jc,jk) = vsed_coeff * (exp_aero**16.0_wp   &
          &          + 1.246_wp * knudsen_aero(jc,jk) * (exp_aero)**4.0_wp)
        ! Eq. 3.69 for k=3
        vsed3(jc,jk) = vsed_coeff * (exp_aero**64.0_wp   &
          &          + 1.246_wp * knudsen_aero(jc,jk) * (exp_aero)**28.0_wp)
        
        CALL art_clip_gt(vsed3(jc,jk),1.0_wp)
        CALL art_clip_gt(vsed0(jc,jk),1.0_wp)
      ELSE
        vsed0(jc,jk) = 0.0_wp
        vsed3(jc,jk) = 0.0_wp
      ENDIF
    ENDDO ! jc
  ENDDO ! jk
  
END SUBROUTINE art_calc_v_sed
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_sed_flx(vsed0,vsed3,wgtfac_c,wgtfacq_c, rho, rho_ic,      &
  &                         istart, iend, nlev, flx_ctra_sed0, flx_ctra_sed3)
!<
! SUBROUTINE art_calc_sed_flx
! This subroutine calculates mass fluxes out of the sedimentation
! velocities and handles interpolation for the bottom layer
! Based on: -
! Part of Module: mo_art_sedi_2mom
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-25
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp), INTENT(in)   :: &
    &  vsed0(:,:),          & !< Sedimentation velocity 0th moment
    &  vsed3(:,:),          & !< Sedimentation velocity 3rd moment
    &  rho(:,:),            & !< Density at full levels
    &  rho_ic(:,:)            !< Density at half levels
  REAL(vp), INTENT(in)   :: &
    &  wgtfac_c(:,:),       & !< weighting factor for interpolation from full to half levels
    &  wgtfacq_c(:,:)         !< weighting factor for quadratic interpolation to surface
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end of inner loop (nproma)
    &  nlev                   !< Number of vertical levels
  REAL(wp), INTENT(inout) :: &
    &  flx_ctra_sed0(:,:),  & !< flux due to sedimentation of zeroth moment
    &  flx_ctra_sed3(:,:)     !< flux due to sedimentation of third moment
! Local variables
  REAL(wp)               :: &
    &  rho_nlevp1             !< Artificial density at ground level (interpolated)
  INTEGER                :: &
    &  nlevp1,              & !< nlev + 1
    &  jc, jk                 !< Loop indices
  
  nlevp1 = nlev + 1
  
  ! ----------------------------------
  ! --- Interpolate on half levels and
  ! --- Calculate mass flux
  ! ----------------------------------    
  
  DO jk = 2, nlev ! all except top and bottom level
    DO jc = istart, iend
      ! First interpolate vsed to half levels and save it in flx_ctra_sed0/3
      flx_ctra_sed0(jc,jk) = wgtfac_c(jc,jk)*vsed0(jc,jk)      &
        &                  + (1._wp-wgtfac_c(jc,jk)) * vsed0(jc,jk-1)
      flx_ctra_sed3(jc,jk) = wgtfac_c(jc,jk)*vsed3(jc,jk)      &
        &                  + (1._wp-wgtfac_c(jc,jk)) * vsed3(jc,jk-1)
      ! Second calculate massflux
      flx_ctra_sed0(jc,jk)= - flx_ctra_sed0(jc,jk)*rho_ic(jc,jk)
      flx_ctra_sed3(jc,jk)= - flx_ctra_sed3(jc,jk)*rho_ic(jc,jk)
    ENDDO  ! jc
  ENDDO ! jk
  
  ! Calculate mass flux for top and bottom level
  DO jc =  istart, iend
    rho_nlevp1 = wgtfacq_c(jc,1)*rho(jc,nlev)     &
       &       + wgtfacq_c(jc,2)*rho(jc,nlev-1)   &
       &       + wgtfacq_c(jc,3)*rho(jc,nlev-2)
    flx_ctra_sed0(jc,1)     = - vsed0(jc,1)    * rho_ic(jc,1)
    flx_ctra_sed0(jc,nlevp1)= - vsed0(jc,nlev) * rho_nlevp1
    flx_ctra_sed3(jc,1)     = - vsed3(jc,1)    * rho_ic(jc,1)
    flx_ctra_sed3(jc,nlevp1)= - vsed3(jc,nlev) * rho_nlevp1
  ENDDO ! jc
  
END SUBROUTINE art_calc_sed_flx
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_sedi_2mom
