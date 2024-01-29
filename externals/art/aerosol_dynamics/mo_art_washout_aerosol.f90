!
! mo_art_washout_aerosol
! This module provides Washout for modal aerosol
!
! [J/K]     Boltzmann constant
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

MODULE mo_art_washout_aerosol
!ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_physical_constants,            ONLY: ak ! [J/K]     Boltzmann constant

IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_aerosol_washout

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_aerosol_washout(temp, p_trac0, rho_aero, diam_aero,sg_aero,       &
              &                qr,rho, dyn_visc,free_path, istart, iend, kstart, &
              &                nlev,lmass_specific, wash_rate_m0, wash_rate_m3,  &
              &                rrconv_2d, rrconv_3d, qnr, iart_aero_washout)
!<
! SUBROUTINE art_aerosol_washout
! Calculates washout for modal aerosol
! Based on: Rinke (2008) - Parametrisierung des Auswaschens von Aerosolpartikeln durch
!                          Niederschlag - Dissertation an der Fakultaet fuer Physik
!                          der Universitaet (TH) Karlsruhe
! Part of Module: mo_art_washout_aerosol
! Author: Daniel Rieger, KIT
! Initial Release: 2013-09-30
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - 2017-08-18: Andrea Steiner, DWD: new switch iart_aero_washout
! - 2017-12-22: Andrea Steiner, DWD: bugfixes

!>
  REAL(wp), INTENT(in) :: &
    &  temp(:,:),         & !< Temperature [K]
    &  p_trac0(:,:),      & !< Aerosol number concentration [# kg-1] or [# m-3] - Specified by lmass_specific
    &  rho_aero(:,:),     & !< Aerosol density [kg m-3]
    &  diam_aero(:,:),    & !< Aerosol median diameter [m]
    &  sg_aero,           & !< Aerosol distribution standard deviation
    &  qr(:,:),           & !< Rain water mass mixing ratio [kg kg-1]
    &  rho(:,:),          & !< Density [kg m-3]
    &  dyn_visc(:,:),     & !< Dynamic viscosity []
    &  free_path(:,:)       !< Mean free path []
  INTEGER, INTENT(in)  :: &
    &  istart, iend,      & !< Start and end indizes of inner loop
    &  kstart, nlev         !< Start and end indizes of vertical loop
  LOGICAL, INTENT(in)  :: &
    &  lmass_specific       !< .TRUE. if p_trac0 is [# kg-1], .FALSE. for [# m-3]
  REAL(wp),INTENT(inout) :: &
    &  wash_rate_m0(:,:), & !< Washout rate with respect to number concentration [# m-3 s-1]
    &  wash_rate_m3(:,:)    !< Washout rate with respect to mass concentration [m3 m-3 s-1]
  REAL(wp),INTENT(in),OPTIONAL :: &
    &  rrconv_2d(:),      & !< 2-Dimensional rain rate (only horizontal information) due to subgridscale convection
    &  rrconv_3d(:,:),    & !< 3-Dimensional rain rate (horizontal and vertical information) due to subgridscale convection
    &  qnr(:,:)             !< Rain water specific number [# kg-1] (rain drop number)
  INTEGER, INTENT(in),OPTIONAL :: &
    &  iart_aero_washout    !< Treatment of grid scale and convective precipitation:
                            !< 0:gscp+con; 1:gscp,con; 2:gscp,rcucov*con
! Local variables
  REAL(wp)              :: &
    &  n_dens_droplets,    & !< number density of droplets set according to liquid water content m-3
    &  numb_droplets,      & !< number of droplets (?)
    &  dnmbdt,             & !< change of number density with time by washout
    &  dnmbdtIMP,          & !< change of number density with time by impaction washout
    &  w_rain,             & !< vertical velocity of rain droplets
    &  gamma(8),           & !< parameters
    &  beta,               & !< parameter
    &  KN(11),             & !< parameters
    &  xeta(8),            & !< parameters
    &  gammaf(3),          & !< parameters
    &  rho_w,              & !< water density
    &  muw =1.792e-3

  REAL(wp),ALLOCATABLE :: &
    &  rr_conv(:,:),       & !< rain rate due to subgridscale convection (local)
    &  numb_aero(:,:),     & !< Aerosol number concentration [# m-3]
    &  lwp_kgm3(:,:)         !< liquid water content in kg m-3


  INTEGER               :: &
    &  i, k                  !< Loop indices

  parameter (rho_w=1000.0_wp)
  parameter (w_rain=10.0_wp)

  ALLOCATE(rr_conv(SIZE(p_trac0,1),SIZE(p_trac0,2)))
  ALLOCATE(numb_aero(SIZE(p_trac0,1),SIZE(p_trac0,2)))
  ALLOCATE(lwp_kgm3(SIZE(p_trac0,1),SIZE(p_trac0,2)))

  wash_rate_m0(:,:) = 0._wp
  wash_rate_m3(:,:) = 0._wp
  rr_conv(:,:)      = 0._wp

  DO k=kstart,nlev
    DO i=istart,iend

       IF (PRESENT(rrconv_2d)) rr_conv(i,k) = rrconv_2d(i)
       IF (PRESENT(rrconv_3d)) rr_conv(i,k) = rrconv_3d(i,k)

       IF (lmass_specific) THEN
         numb_aero(i,k) = p_trac0(i,k) * rho(i,k)
       ELSE
         numb_aero(i,k) = p_trac0(i,k)
       ENDIF

       ! Calculate liquid water content in kg m-3
       IF (PRESENT(iart_aero_washout)) THEN
         !------------------------------------
         !different ways to treat grid scale and convective precipitation:
         !iart_aero_washout=0: "gscp+con":
         !                     washout with gscp+con, only 1 call to washout routine
         !iart_aero_washout=1: "gscp,con":
         !                     washout with gscp, con separately, 2 calls to washout routine
         !iart_aero_washout=2: "gscp,con/rcucov":
         !                     washout with gscp, con separately, con scaled with rcucov,
         !                     2 calls to washout routine
         !NOTE: to indicate second call to washout routine, 100 is added to iart_aero_washout
         !      when art_aerosol_washout is called the second time
         !------------------------------------
         IF (iart_aero_washout == 0) THEN
           !first and only call to washout routine
           lwp_kgm3(i,k) = qr(i,k) * rho(i,k) + (rr_conv(i,k) / w_rain)
         ELSE IF (iart_aero_washout == 1 .OR. iart_aero_washout == 2) THEN
           !first call to washout routine with iart_aero_washout 1/2:
           !do washout only with grid scale precipitation
           lwp_kgm3(i,k) = qr(i,k) * rho(i,k)
         ELSE IF (iart_aero_washout == 101 .OR. iart_aero_washout == 102) THEN
           !second call to washout routine with iart_dust_waschot 1/2->101/102:
           !do washout only with convective precipitation
           lwp_kgm3(i,k) = (rr_conv(i,k) / w_rain)
         ENDIF

       ELSE !Default
         lwp_kgm3(i,k) = qr(i,k) * rho(i,k) + (rr_conv(i,k) / w_rain)
       ENDIF


       ! If aerosols and liq water are "present" --> do washout
       IF (numb_aero(i,k) >= 1.0_wp .AND. lwp_kgm3(i,k) > 1.0E-15_wp) THEN

         IF(PRESENT(qnr)) THEN
           n_dens_droplets = qnr(i,k) * rho(i,k)
         ELSE
           n_dens_droplets=1.E7_wp                            ! low precipitation
           IF (lwp_kgm3(i,k) > 0.001_wp) n_dens_droplets=5000.0_wp ! moderate precipitation
           IF (lwp_kgm3(i,k) > 0.005_wp) n_dens_droplets=500.0_wp  ! heavy precipitation
         ENDIF

         ! ----------------------------------
         ! --- Calculate the parameters for MARSHALL-PALMER Distribution
         ! ----------------------------------

         beta = (lwp_kgm3(i,k)/(pi*rho_w*n_dens_droplets))**(-1.0_wp/3.0_wp)
         numb_droplets = (n_dens_droplets**(4.0_wp/3.0_wp))*(lwp_kgm3(i,k)/(pi*rho_w))**(-1.0_wp/3.0_wp)


         gammaf(1) = (21.0_wp/16.0_wp) * (SQRT(2.0_wp)*pi) / 3.6256_wp   ! gammaf1=gamma(11/4)
         gammaf(2) = (5.0_wp/16.0_wp) * 3.6256_wp                        ! gammaf2=gamma(9/4)
         gammaf(3) = (231.0_wp/64.0_wp) * (SQRT(2.0_wp)*pi) / 3.6256_wp  ! gammaf3=gamma(15/4)

         xeta(1) = numb_droplets * (1.0_wp/beta**2.0_wp)
         xeta(2) = numb_droplets * (gammaf(1)/beta**(11.0_wp/4.0_wp))
         xeta(3) = numb_droplets * (gammaf(1)/beta**(11.0_wp/4.0_wp))
         xeta(4) = numb_droplets * (3.0_wp * SQRT(pi)/(4.0_wp*beta**(5.0_wp/2.0_wp)))
         xeta(5) = numb_droplets * (SQRT(pi)/(2.0_wp*beta**(3.0_wp/2.0_wp)))
         xeta(6) = numb_droplets * (gammaf(2)/beta**(9.0_wp/4.0_wp))
         xeta(7) = numb_droplets * ((15.0_wp*SQRT(pi))/(8.0_wp * beta**(7.0_wp/2.0_wp)))
         xeta(8) = numb_droplets * (gammaf(3)/(beta**(15.0_wp/4.0_wp)))

         gamma(1) = (ak * temp(i,k))/(6.0_wp*dyn_visc(i,k))
         gamma(2) = 0.4_wp * SQRT(130.0_wp) * pi/4.0_wp            &
            &   * (SQRT(2.0_wp*dyn_visc(i,k)/(rho(i,k)))) * (((ak*temp(i,k)*rho(i,k)) &
            &   /(3.0_wp*pi*dyn_visc(i,k)**2.0_wp))**(2.0_wp/3.0_wp))
         gamma(3) = 0.16_wp * SQRT(130.0_wp) * (pi/4.0_wp)         &
            &   * SQRT((2.0_wp*ak*temp(i,k))/(3.0_wp*pi*dyn_visc(i,k)))
         gamma(4) = 4.0_wp * (dyn_visc(i,k)/muw) * (pi/4.0_wp) * 130.0_wp
         gamma(5) = 4.0_wp * 130.0_wp * (pi/4.0_wp)
         gamma(6) = 4.0_wp * SQRT((rho(i,k))/(2.0_wp*dyn_visc(i,k))) * (pi/4.0_wp) &
            &   * (130.0_wp)**(3.0_wp/2.0_wp)
         ! Impaction
         gamma(7) = (pi/4.0_wp) * 130.0_wp * SQRT(rho_w/rho_aero(i,k))
         ! drieg: Bugfix: Added Minus
         gamma(8) = -0.9_wp * (pi/4.0_wp) * (((130.0_wp*9.0_wp*dyn_visc(i,k)*rho_w) &
            &   / (rho_aero(i,k)**(2.0_wp)))**(1.0_wp/2.0_wp))

         KN(1) = gamma(1) * xeta(1)
         KN(2) = gamma(1) * xeta(1) * (3.34_wp * free_path(i,k))
         KN(3) = gamma(2) * xeta(2) * ((3.34_wp * free_path(i,k))**(2.0_wp/3.0_wp))
         ! drieg: Bugfix: Added Minus in the Exponent (cf. Jung et al. 2003)
         KN(4) = gamma(2) * xeta(2) * (2.0_wp/3.0_wp) * ((3.34_wp * free_path(i,k))**(-1.0_wp/3.0_wp))
         KN(5) = gamma(3) * xeta(3) * ((3.34_wp * free_path(i,k))**(1.0_wp/2.0_wp))
         ! drieg: Bugfix: Added Minus in the Exponent (cf. Jung et al. 2003)
         KN(6) = gamma(3) * xeta(3) * (1.0_wp/2.0_wp) * ((3.34_wp * free_path(i,k))**(-1.0_wp/2.0_wp))
         KN(7) = gamma(4) * xeta(4)
         KN(8) = gamma(5) * xeta(5)
         KN(9) = gamma(6) * xeta(6)
         ! Impaction
         KN(10) = gamma(7) * xeta(7)
         KN(11) = gamma(8) * xeta(8)

         ! ----------------------------------
         ! --- Calculate update frequencies for washout
         ! ----------------------------------

         ! DIFFUSION and INTERCEPTION
         dnmbdt= ((KN(1) * numb_aero(i,k) * diam_aero(i,k)**(-1.0_wp)             &
             &           * (EXP((1.0_wp/2.0_wp)*(LOG(sg_aero))**(2.0_wp)))))      &
             & + ((KN(2) * numb_aero(i,k) * diam_aero(i,k)**(-2.0_wp)             &
             &           * (EXP((2.0_wp)*(LOG(sg_aero))**(2.0_wp)))))             &
             & + ((KN(3) * numb_aero(i,k) * diam_aero(i,k)**(-4.0_wp/3.0_wp)      &
             &           * (EXP((8.0_wp/9.0_wp)*(LOG(sg_aero))**(2.0_wp)))))      &
             & + ((KN(4) * numb_aero(i,k) * diam_aero(i,k)**(-1.0_wp/3.0_wp)      &
             &           * (EXP((1.0_wp/18.0_wp)*(LOG(sg_aero))**(2.0_wp)))))     &
             & + ((KN(5) * numb_aero(i,k) * diam_aero(i,k)**(-1.0_wp)             &
             &           * (EXP((1.0_wp/2.0_wp)*(LOG(sg_aero))**(2.0_wp)))))      &
             & +  (KN(6) * numb_aero(i,k))                                        &
             & + ((KN(7) * numb_aero(i,k) * diam_aero(i,k)                        &
             &           * (EXP((1.0_wp/2.0_wp)*(LOG(sg_aero))**(2.0_wp)))))      &
             & + ((KN(8) * numb_aero(i,k) * diam_aero(i,k)**(2.0_wp)              &
             &           * (EXP((2.0_wp)*(LOG(sg_aero))**(2.0_wp)))))             &
             & + ((KN(9) * numb_aero(i,k) * diam_aero(i,k)**(2.0_wp)              &
             &           * (EXP((2.0_wp)*(LOG(sg_aero))**(2.0_wp)))))

         ! IMPACTION
         dnmbdtIMP =  (KN(10)* numb_aero(i,k))                     &
             &     +  (KN(11)* numb_aero(i,k) * diam_aero(i,k)**(-1.0_wp)  &
             &               * (EXP((1.0_wp/2.0_wp)*(LOG(sg_aero))**(2.0_wp))))

         IF (dnmbdtIMP >= 0.0_wp) dnmbdt = dnmbdt + dnmbdtIMP

         wash_rate_m0(i,k) = (-1.0_wp)*dnmbdt
         wash_rate_m3(i,k) = (-1.0_wp)*dnmbdt * (diam_aero(i,k))**(3.0_wp) &
                &          * EXP(4.5_wp *(LOG(sg_aero))**(2.0_wp))
       ENDIF
    ENDDO
  ENDDO

  DEALLOCATE(rr_conv,numb_aero,lwp_kgm3)

END SUBROUTINE art_aerosol_washout
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_washout_aerosol
