!
! mo_art_depo_2mom
! This module provides the calculation of the deposition
! velocities of log-normally distributed two-moment aerosol
! Based on N. Riemer - Numerische Simulationen zur Wirkung des Aerosols
!             (2002)   auf die troposphaerische Chemie und die Sichtweite
!                      Dissertation
!                      Fakultaet fuer Physik, Universitaet Karlsruhe
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

MODULE mo_art_depo_2mom
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_physical_constants,            ONLY: grav,ak ! gravitational acceleration [m s-2] / Boltzmann constant [J K-1]
  USE mo_fortran_tools,                 ONLY: t_ptr_tracer
! ART
  USE mo_art_aerosol_utilities,         ONLY: calc_aerodynamic_resistance,calc_ustar
  USE mo_art_clipping,                  ONLY: art_clip_gt

IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_calc_v_dep
  PUBLIC :: art_store_v_dep

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_v_dep(temp,temp_sfc,u,v,rho,tcm,tch,qv,qc,qr,theta_v,gz0,dz, &
  &                       dyn_visc,diam_aero,exp_aero, knudsen_aero,vsed0,vsed3, &
  &                       istart,iend,vdep0,vdep3)
!<
! SUBROUTINE art_calc_v_dep
! This subroutine calculates deposition velocities
! Based on: Riemer (2002)
! Part of Module: mo_art_depo_2mom
! Author: Daniel Rieger, KIT
! Initial Release: 2014-06-25
! Modifications:
! 2014-11-24: Daniel Rieger, KIT
! - Adapted to unified ART-physics structure
!>
  REAL(wp), INTENT(in)   :: &
    &  temp(:),             & !< Temperature in lowest layer [K]
    &  temp_sfc(:),         & !< Surface temperature [K]
    &  u(:),v(:),           & !< Horizontal wind in lowest layer [m s-1]
    &  rho(:),              & !< Air density [kg m-3]
    &  tcm(:), tch(:),      & !< Transfer coefficients for momentum and heat in lowest layer [--]
    &  qv(:),qc(:),qr(:),   & !< Mass mixing ratios for water vapor / cloud water / rain water [kg kg-1]
    &  theta_v(:),          & !< Virtual potential temperature [K]
    &  gz0(:),              & !< roughness length * g of the vertically not resolved canopy [m2 s-2]
    &  dz(:),               & !< Layer height [m]
    &  dyn_visc(:),         & !< Dynamic viscosity
    &  diam_aero(:),        & !< Aerosol median diameter (with respect to number conc.)
    &  exp_aero,            & !< Aerosol exponent
    &  knudsen_aero(:),     & !< Aerosol knudsen number
    &  vsed0(:),            & !< Sedimentation velocity of 0th moment in lowest layer [m s-1]
    &  vsed3(:)               !< Sedimentation velocity of 3rd moment in lowest layer [m s-1]
  INTEGER, INTENT(in)    :: &
    &  istart, iend           !< Start and end of inner loop (nproma)
  REAL(wp), INTENT(inout):: &
    &  vdep0(:), vdep3(:)     !< Deposition velocities for 0th / 3rd moment
! Local variables
  REAL(wp)               :: &
    &  diff_const,          & !< Eq. 3.68 factor
    &  diff_coeff0,         & !< Eq. 3.68 coefficient of Moment 0
    &  diff_coeff3,         & !< Eq. 3.68 coefficient of Moment 3
    &  schmidt0,schmidt3,   & !< Schmidt Number of Moment 0 / 3
    &  stokes0,stokes3,     & !< Stokes Number of Moment 0 / 3
    &  ustar,               & !< Friction velocity
    &  wstar,               & !< convective velocity (currently set constant)
    &  raerody,             & !< aerodynamical resistance
    &  nu,                  & !< nu = dynamic viscosity / density
    &  ustfac,              & !< Prefactor with ustar**2
    &  utscale,             & !< Prefactor
    &  twothird,            & !< 2.0/3.0
    &  rd0, rd3               !<

  INTEGER                :: &
    &  jc                     !< Loop index

  wstar    = 1.0_wp
  twothird = 2.0_wp / 3.0_wp

!NEC$ ivdep
  DO jc = istart, iend
    IF (diam_aero(jc) > 0._wp) THEN
      ! Eq. 3.68 factor
      diff_const   = ak * temp(jc) / (3.0_wp * pi * dyn_visc(jc))  / diam_aero(jc)
      ! Eq. 3.68 for k=0
      diff_coeff0 = diff_const * (exp_aero**(4.0_wp)   &
        &        + 1.246_wp * knudsen_aero(jc) * exp_aero**16.0_wp)
      ! Eq. 3.68 for k=3
      diff_coeff3 = diff_const * (exp_aero**(-20.0_wp) &
        &        + 1.246_wp * knudsen_aero(jc) * exp_aero**(-32.0_wp))

      ! Eq. 3.71 precalculate several factors
      CALL calc_ustar(ustar,tcm(jc),u(jc),v(jc))
      CALL calc_aerodynamic_resistance(ustar, theta_v(jc), gz0(jc), qv(jc), qc(jc),   &
        &                              qr(jc), temp(jc), temp_sfc(jc),dz(jc),tch(jc),tcm(jc),raerody)

      nu      = dyn_visc(jc) / rho(jc)
      ustfac  = ustar * ustar / (grav * nu)
      utscale = ustar + 0.24_wp * wstar * wstar / ustar

      ! ----------------------------------
      ! --- Moment 0
      ! ----------------------------------

      schmidt0 = nu / diff_coeff0 ! Schmidtzahl MOM0
      stokes0  = MAX(vsed0(jc) * ustfac , 0.01e0_wp)
      rd0 = 1.0_wp / (utscale * (schmidt0**(-twothird) + 10._wp**(-3._wp / stokes0))) ! Eq. 3.71 complete

      ! Eq. 3.70
! JF:       vdep0(jc) = vsed0(jc) + 1.0_wp / (raerody + rd0 + rd0 * raerody * vsed0(jc))  ! original version
! JF:       vdep0(jc) = 1.0_wp / (raerody + rd0 + rd0 * raerody * vsed0(jc))              ! version with parallel resistance term
      vdep0(jc) = 1.0_wp / (raerody + rd0)

      ! ----------------------------------
      ! --- Moment 3
      ! ----------------------------------

      schmidt3 = nu / diff_coeff3
      stokes3 = MAX( vsed3(jc) * ustfac, 0.01e0_wp)
      rd3 = 1.0_wp / (utscale * (schmidt3**(-twothird) + 10.0_wp**(-3.0_wp / stokes3))) ! Eq. 3.71 including impaction term
      ! Seinfeld and Pandis: For very large particles, impaction and interception lead
      !   to effective removal. The result is that 1/rb achieves a minimum in the range between 0.1
      !   and 1.0 mum diameter where none of the removal process is especially effective.

! JF:       vdep3(jc) = vsed3(jc) + 1.0_wp / (raerody + rd3 + rd3 * raerody * vsed3(jc))  ! original version
! JF:       vdep3(jc) = 1.0_wp / (raerody + rd3 + rd3 * raerody * vsed3(jc))              ! version with parallel resistance term
      vdep3(jc) = 1.0_wp / (raerody + rd3)

      CALL art_clip_gt(vdep0(jc),0.3_wp)
      CALL art_clip_gt(vdep3(jc),0.3_wp)
    ELSE
      vdep0(jc) = 0.0_wp
      vdep3(jc) = 0.0_wp
    ENDIF ! diam > 0
  ENDDO !jc

END SUBROUTINE art_calc_v_dep
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_store_v_dep(vdep0, vdep3, njsp, jsp, i_number_conc, nturb_tracer,jb, &
  &                        turb_tracer,istart,iend,vdep_turb)
!<
! SUBROUTINE art_store_v_dep
! This subroutine stores deposition velocities in the structure
! that will be passed to the turbulence routine
! Based on: Riemer (2002)
! Part of Module: mo_art_depo_2mom
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-25
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp), INTENT(in)   :: &
    &  vdep0(:),vdep3(:)      !< Deposition velocity of 0th and 3rd moment
  INTEGER, INTENT(in)    :: &
    &  njsp,                & !< Number of species contained in current mode
    &  jsp(:),              & !< Index array of mass species in current mode
    &  i_number_conc,       & !< Index of number species in current mode
    &  nturb_tracer           !< Number of tracers that underlie turbulence
  TYPE(t_ptr_tracer),INTENT(in) :: &
    &  turb_tracer(:,:)
  INTEGER,INTENT(in)     :: &
    &  jb,istart,iend
  REAL(wp), INTENT(inout):: &
    &  vdep_turb(:,:,:)       !< Deposition velocity storage array for turbulence routine
! Local variables
  INTEGER                :: &
    &  iturb,               & !< Counter for turb tracer array
    &  i,                   & !< Counter for species
    &  jc,                  & !< Counter inner loop
    &  jspturb                !< Index of species in turbulent tracer array

  ! ----------------------------------
  ! --- loop through the tracers contained in the mode ( i=0 number , i>0 mass mixing ratios)
  ! ----------------------------------
  DO i=0, njsp
    IF (i /= 0) THEN !< get index of mass mixing ratio of the species contained in this_mode
      DO iturb = 1, nturb_tracer
        IF (turb_tracer(jb,iturb)%idx_tracer == jsp(i)) THEN
          ! Get the index of the turbulent tracer field
          jspturb = iturb
        ENDIF
      ENDDO
      DO jc = istart, iend
        vdep_turb(jc,jb,jspturb) = vdep3(jc)
      ENDDO
    ELSE               !< get index of number mixing ratio of this_mode
      DO iturb = 1, nturb_tracer
        IF (turb_tracer(jb,iturb)%idx_tracer == i_number_conc) THEN
          ! Get the index of the turbulent tracer field
          jspturb = iturb
        ENDIF
      ENDDO
      DO jc = istart, iend
        vdep_turb(jc,jb,jspturb) = vdep0(jc)
      ENDDO
    ENDIF
  ENDDO

END SUBROUTINE art_store_v_dep
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_depo_2mom
