!
! mo_art_psc_sts
! This module calculates the surface density and uptake coefficients of
! different heterogeneous reactions on STS polar stratospheric clouds
!
! Based on K. Carslaw - An analytic expression for the compositin of aqueous HNO3-H2SO4
!             (1995)    stratospheric aerosols including gas phase removal of HNO3
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

MODULE mo_art_psc_sts

! =====================================================================
!    FORTRAN77 CODE TO CALCULATE COMPOSITION OF AQUEOUS
!   HNO3/H2SO4/HCL/HBR/HOCL/HOBR STRATOSPHERIC AEROSOLS
!     (CARSLAW, LUO, PETER - GEOPHYS. RES. LETT., 1995)
!             CARSLAW@NIKE.MPCH-MAINZ.MPG.DE
!
!            DISTRIBUTION VERSION 2 (27 MAY 1996)
!      UPDATES FROM VERSION 1: HBR AND HOBR SOLUBILITIES
!
!         KEN CARSLAW
!         MAX-PLANCK-INSTITUT FUER CHEMIE
!         POSTFACH 3060
!         MAINZ 55020
!         GERMANY
!         TEL: (49) (0)6131 305 333
!         FAX: (49) (0)6131 305 328
!
! =====================================================================
! - hno3/h2so4 composition based in thermodynamic model of the system
!   hcl / hbr / hno3 / h2so4 / h2o
!   {carslaw et al, jpc, 1995 -
!    a thermodynamic model of the system hcl-hno3-h2so4-h2o,
!    including solubilities of hbr, from < 200 k to 328 k.
!    (also available in fortran77)}
! - hcl solubility parametrisation from luo el al., grl, 1995
! - hocl solubility parametrisation from huthwelker et al., jac, 1995
! *********************************************************************
! - the model is valid for 2e-5 mb < pw < 2e-3 mb
!   (pw is the water partial pressure)
! - the upper temperature limit is 240 k
!   (it extrapolates smoothly to higher temperatures, but compositions
!    are not reliable)
! - the lower temperature limit is 185 k, or tice - 3 k, which ever is
!   the higher
! - hno3 solubilities are calculated to a maximum temperature of 215 k
! - solubilities are calculated on a molality basis -
!   (moles of solute per kg of water)
! *********************************************************************
! - the solubilities of hcl, hocl, hbr and hobr are assumed not to
!   affect the uptake of hno3 or h2o.
!   this introduces only small errors at the very lowest temperatures,
!   but for a full calculation where interactions between all species
!   in solution are considered, use the model of carslaw et al.,
!   given above.
! =====================================================================

  ! ICON
  USE mo_kind,                    ONLY: wp
  USE mo_parallel_config,         ONLY: nproma
  USE mo_math_constants,          ONLY: pi
  USE mo_physical_constants,      ONLY: argas, amw, avo
  ! ART
  USE mo_art_psc_types,           ONLY: t_art_psc,  &
                                    &   mol_weight_H2SO4

  USE messy_mecca_kpp_global,     ONLY: ihs_N2O5_H2O , ihs_HOCl_HCl, ihs_ClNO3_HCl, &
                                    &   ihs_ClNO3_H2O, ihs_HOCl_HBr, ihs_HOBr_HCl , &
                                    &   ihs_HOBr_HBr , ihs_BrNO3_H2O
!==============================================================================



  IMPLICIT NONE

  PRIVATE


  PUBLIC :: art_psc_aerosol

!==============================================================================

CONTAINS

!************************************************************************

SUBROUTINE art_psc_aerosol (PSC, jb, kstart,kend, i_startidx, i_endidx, &
                  &         temp,pres,p_tracer_now)
!<
! SUBROUTINE art_psc_aerosol
! This subroutine calculates the surface density and uptake coefficients of
! different heterogeneous reactions on STS polar stratospheric clouds
! Based on Carslaw et al. (1995) see above
! - Eckstein (2013), diploma thesis, Simulation der polaren unteren
!   Stratosphaere mit den regionalen Chemie-Transport Modell COSMO-ART
! Part of Module: mo_art_psc_sts
! Author: Michael Weimer, KIT
! Initial Release: 2017-02-06
! Modifications:
!>

  IMPLICIT NONE
  
  TYPE(t_art_psc), INTENT(inout) :: &
     &  PSC                    !< structure for PSC meta data
  INTEGER, INTENT(in) ::            &
     &  jb, i_startidx, i_endidx,   &
     &   kstart,kend           !< loop indices
  REAL(wp), INTENT(in), DIMENSION(:,:,:)   ::  &
     &  temp, pres             !< air temperature (K), and air pressure (Pa)
  REAL(wp), INTENT(in) :: &
     &  p_tracer_now(:,:,:,:)  !< tracer number concentrations (# / cm3)
  ! local variables
  INTEGER ::    &
     &  jc, jk                 !< loop indices
  
  REAL(wp) ::              &
     &   Nconc2p,          &   !< conversion factor from number concentration (#/cm3) to partial
                               !  pressure (Pa)
     &   HNO3_max_Nconc,   &   !< maximum value for the HNO3 number concentration (#/cm3)
     &   mtot,             &   !< mass concentration in liquid aerosol (g / cm3)
     &   hhcl,             &   !< hhcl  = effective hcl henry's law constant 
                               !  in hno3-h2so4-h2o solution [mol / kg water / atm]
     &   pn,               &   !< pn    = hno3 vapour pressue [atm] 
                               !  (calculated after partitioning of hno3 to the aerosol)
     &   phbr,             &   !< (pbr) wohl phbr = hbr vapour pressue [atm] 
                               !  (calculated after partitioning of hbr to the aerosol)
     &   phcl,             &   !< (pcl) wohl phcl = hcl vapour pressue [atm] 
                               !  (calculated after partitioning of hcl to the aerosol)
     &   mnb,              &   !< mnb   = binary solution concentrations [mole / kg water]
     &   msb,              &   !< msb   = binary solution concentrations [mole / kg water]
     &   mhocl,            &   !< mhocl = concentration of hocl [mole / kg water]
     &   mhobr,            &   !< mhobr = concentration of hobr [mole / kg water]
     &   mcl,              &   !< mcl   = concentration of hcl [mole / kg water]
     &   mbr,              &   !< mbr   = concentration of hbr [mole / kg water]
     &   hsb,              &   !< hsb   = effective henry's law constant of hno3 
                               !          in h2so4-h20 solution [mol / kg water / atm]
     &   hnb,              &   !< hnb   = effective henry's law constant of hno3 
                               !          in hno3-h20 solution [mol / kg water / atm]
     &   vliq,             &   !< vliq  = specific aerosol volume (dimensionless)
     &   xsb, xnb,         &   !< xsb fraction of binary H2SO4/H2O and HNO3/H2O solution,
     &   nsms,             &   !< ns / ms
     &   mnb2, msb2,       &   !< mnb squared, msb squared
     &   sqrtws, sqrtwn,   &   !< sqrt of ws and wn
     &   ws2, wn2,         &   !< ws and wn squared
     &   k1hoclhcl,        &   !< first order rate constants for reaction in liquid [dm3/mol/s]
     &   k1hobrhcl,        &   !  ...
     &   k1hbrhobr,        &   !  ...
     &   k1hbrhocl,        &   !  ...
     &   ah2o,             &   !< water activity
     &   fhhbr,            &   !< effective henry constants converted to mol/dm3/atm
     &   fhhocl,           &   !< ...
     &   fhhobr,           &   !< ...
     &   p,                &   !< ratio of first-order loss rate coeff for ClONO2 + HCl 
                               !  to that of ClONO2 + H2O (see Hanson and Ravishankara, 1994)
     &   q, fq,            &   !< These come from Luo et al., 1995, 
                               !  but this publication does not exist (??)
     &   sqrtfunc,         &   !< auxillary variable for molality of H2SO4 in water
     &   rt,               &   !< atm2Pa / (R * T) with R = universal gas constant and T = temp.
     &   rmean,            &   !< effective radius of STS particles (cm)
     &   tn, tt,           &   !< tn: temperature (K), tt: partial pressure of H2SO4 in water (atm)
     &   tr, tr2, tr3,     &   !< fitting parameters for binary solutions (wrt temperature)
     &   pr, pr2, pr3,     &   !< fitting parameters for binary solutions (wrt pressure)
     &   a, b, c, cbar,    &   !< a,b,c: parameters for ternary solution;
                               !  cbar: parameter for uptake coefficient
     &   logvon760,        &   !< logarithms of numbers
     &   logvon10,         &   !  ...
     &   g0, gs, gcalc,    &   !< g0: uptake coefficient for ClONO2 + H2O without HCl present; 
                               !  gs: uptake coeff due to surface reactions;  
                               !  gcalc: uptake coeff due to bulk reactions
     &   ge, grxn,         &   !< ge: the overall uptake coeff for C1ONO2; 
                               !  grxn: auxillary variable to calculate 
                               !        uptake coefficient of BrONO2 and water
     &   logkorr, logt,    &   !< logkorr: auxillary variable to calculate tice;
                               !  logt: logarithm of temperature
     &   phi,              &   !< parameter for ternary solution
     &   dl,               &   !< diffusion constant of hocl in h2so4/hno3 solution (cm2 / s)
     &   dens,             &   !< the solution density [g/cm3]
     &   factor,           &   !< conversion factor of effective henry's law constants 
                               !  to mol/dm3/atm
     &   chclliq               !< liquid phase concentration [mol/dm3/atm]

  REAL (KIND=wp), PARAMETER :: &
     &  r=8.205e-5_wp,            & !< argas / atm2Pa
     &  ksur = 576.0_wp,          & !< surface reaction rate term for ClONO2
     &  alpha = 0.30_wp,          & !< constant for uptake coeff for ClONO2
     &  rho = 2.e3_wp               !< ratio of second-order rate coeff for ClONO2+HCl vs ClONO2+H2O


  ! fitting parameters
  REAL(wp), DIMENSION(10) :: &
    &  qnfit = (/14.5734_wp, 0.0615994_wp, -1.14895_wp, 0.691693_wp, -0.098863_wp,        &
    &              0.0051579_wp, 0.123472_wp, -0.115574_wp, 0.0110113_wp, 0.0097914_wp/), &
    &  qsfit = (/14.4700_wp, 0.0638795_wp, -3.29597_wp, 1.778224_wp, -0.223244_wp,        &
    &              0.0086486_wp, 0.536695_wp, -0.335164_wp, 0.0265153_wp, 0.0157550_wp/)
  REAL(wp), DIMENSION(7) :: &
    &  knfit = (/-39.136_wp, 6358.4_wp, 83.29_wp, -17650._wp, 198.53_wp, -11948._wp, -28.469_wp/), &
    &  ksfit = (/-21.661_wp, 2724.2_wp, 51.81_wp, -15732._wp, 47.004_wp,  -6969._wp, -4.6183_wp/)

  logvon760 = log(760.0_wp)
  logvon10  = log(10.0_wp)

  DO jk = kstart,kend
    DO jc = i_startidx,i_endidx
      PSC%liqsur(jc,jk,jb) = 0._wp
      PSC%cgaml(jc,jk,jb,:) = 0._wp
      PSC%parthno3(jc,jk)= 1.0_wp
      PSC%parthcl(jc,jk) = 1.0_wp
      PSC%parthbr(jc,jk) = 1.0_wp
    END DO
  END DO

  ! calculate partial pressure of relevant gaseous substances (in atm)
  DO jk = kstart,kend
    DO jc = i_startidx,i_endidx
      Nconc2p = argas * temp(jc,jk,jb) / avo * 1e6_wp

      ! 2.e-5 mb = 2.e-3 Pa and 2.e-3 mb = 2.e-1 Pa
      PSC%pw(jc,jk)     = min(max(PSC%H2O_Nconc_g(jc,jk) * Nconc2p,   &
                       &      2.e-3_wp),2.e-1_wp)/101325._wp


      HNO3_max_Nconc = 2.0e-8_wp * pres(jc,jk,jb) * avo / argas / temp(jc,jk,jb) * 1.e-6_wp
      IF (PSC%HNO3_Nconc_g(jc,jk) > HNO3_max_Nconc) THEN
        PSC%pn0(jc,jk)    =  HNO3_max_Nconc*Nconc2p/101325._wp
      ELSE
        PSC%pn0(jc,jk)    =  PSC%HNO3_Nconc_g(jc,jk)*Nconc2p/101325._wp
      END IF

      IF (PSC%iTRHCl /= 0) THEN
        PSC%phcl0(jc,jk)  =  p_tracer_now(jc,jk,jb,PSC%iTRHCl)*Nconc2p/101325._wp
      ELSE
        PSC%phcl0(jc,jk)  = 0.0_wp
      END IF

      IF (PSC%iTRHBr /= 0) THEN
        PSC%phbr0(jc,jk)  =  p_tracer_now(jc,jk,jb,PSC%iTRHBr)*Nconc2p/101325._wp
      ELSE
        PSC%phbr0(jc,jk)  = 0.0_wp
      END IF

      IF (PSC%iTRHOCl /= 0) THEN
        PSC%phocl0(jc,jk) =  p_tracer_now(jc,jk,jb,PSC%iTRHOCl)*Nconc2p/101325._wp
      ELSE
        PSC%phocl0(jc,jk)  = 0.0_wp
      END IF

      IF (PSC%iTRHOBr /= 0) THEN
        PSC%phobr0(jc,jk) =  p_tracer_now(jc,jk,jb,PSC%iTRHOBr)*Nconc2p/101325._wp
      ELSE
        PSC%phobr0(jc,jk)  = 0.0_wp
      END IF

      PSC%logpw(jc,jk)    = log(PSC%pw(jc,jk))
      logkorr   = (PSC%logpw(jc,jk)+logvon760)/logvon10
      PSC%tice(jc,jk)     = 2668.7_wp / ( 10.431_wp - logkorr)

      PSC%internal_temp(jc,jk) = temp(jc,jk,jb)

      IF ((temp(jc,jk,jb) < 185._wp)  .OR. (temp(jc,jk,jb) < PSC%tice(jc,jk)-3._wp)) THEN
        PSC%internal_temp(jc,jk) = MAX(185._wp,PSC%tice(jc,jk)-3._wp)
      END IF

      PSC%sqrtt(jc,jk)    = sqrt(PSC%internal_temp(jc,jk))
      ! H2SO4 is not fractionated (there is no sulphuric acid in liquid phase)
      ! so use the whole number concentration of H2SO4 for computation
      PSC%ns(jc,jk)   = p_tracer_now(jc,jk,jb,PSC%iTRH2SO4)*1.e6_wp / avo
    END DO
  END DO


  ! calculation of iexception based on Eckstein (2013);
  ! the parametrisation can explode under special circumstances included here
  DO jk = kstart,kend
    DO jc = i_startidx, i_endidx
      ! low numbers were one degree higher, high numbers one degree lower in new iexception
      IF (PSC%internal_temp(jc,jk) < 185._wp     &
             &   .OR. PSC%internal_temp(jc,jk) < PSC%tice(jc,jk)-3._wp   &
       &    .OR. PSC%internal_temp(jc,jk) > 240._wp .OR. PSC%pw(jc,jk)*1013.25_wp < 2.e-5_wp      &
       &    .OR. PSC%pw(jc,jk)*1013.25_wp > 2.e-3_wp                                              &
       &    .OR. ( PSC%internal_temp(jc,jk) >= 209._wp                                            &
       &    .AND. PSC%internal_temp(jc,jk) <= 215._wp                                             &
       &    .AND. PSC%internal_temp(jc,jk) >= (203.25_wp + 2.5e+05_wp*PSC%pn0(jc,jk) )            &
       &    .AND. PSC%internal_temp(jc,jk) <= (205.75_wp + 2.5e+05_wp*PSC%pn0(jc,jk) )            &
       &    .AND. PSC%internal_temp(jc,jk) >= (206.77_wp + 1.1e+08_wp*PSC%pw(jc,jk) ) )           &
       &    .OR.                                                                                  &
       &    ( PSC%internal_temp(jc,jk) <= 215._wp                                                 &
       &    .AND. PSC%internal_temp(jc,jk) >= (206.77_wp + 1.1e+08_wp*PSC%pw(jc,jk) )             &
       &    .AND. PSC%internal_temp(jc,jk) <= (209.27_wp + 1.1e+08_wp*PSC%pw(jc,jk) )             &
       &    .AND. PSC%internal_temp(jc,jk) <= (205.75_wp + 2.5e+05_wp*PSC%pn0(jc,jk) ) )          &
       &  ) THEN
        PSC%iexception(jc,jk) = 1
      ELSE 
        PSC%iexception(jc,jk) = 0
      END IF

      ! This is again just to see of there is an area, where the new iexception was used
      IF (                                                                                   &
       &    ( PSC%internal_temp(jc,jk) >= 209._wp .AND. PSC%internal_temp(jc,jk) <= 215._wp  &
       &    .AND. PSC%internal_temp(jc,jk) >= (203.25_wp + 2.5e+05_wp*PSC%pn0(jc,jk) )       &
       &    .AND. PSC%internal_temp(jc,jk) <= (205.75_wp + 2.5e+05_wp*PSC%pn0(jc,jk) )       &
       &    .AND. PSC%internal_temp(jc,jk) >= (206.77_wp + 1.1e+08_wp*PSC%pw(jc,jk) ) )      &
       &    .OR.                                                                             &
       &    ( PSC%internal_temp(jc,jk) <= 215._wp                                            &
       &    .AND. PSC%internal_temp(jc,jk) >= (206.77_wp + 1.1e+08_wp*PSC%pw(jc,jk) )        & 
       &    .AND. PSC%internal_temp(jc,jk) <= (209.27_wp + 1.1e+08_wp*PSC%pw(jc,jk) )        &
       &    .AND. PSC%internal_temp(jc,jk) <= (205.75_wp + 2.5e+05_wp*PSC%pn0(jc,jk) ) )     &
       &  ) THEN
        PSC%iexception(jc,jk) = 2
      END  IF

    END DO
  END DO

  DO jk = kstart,kend
    DO jc = i_startidx,i_endidx

      IF (PSC%iexception(jc,jk) == 0) THEN

        logt = log(PSC%internal_temp(jc,jk))
        tt = r * PSC%internal_temp(jc,jk) * PSC%ns(jc,jk)
        tr = 1.e4_wp / PSC%internal_temp(jc,jk) - 43.4782608_wp
        pr = PSC%logpw(jc,jk) + 18.4_wp
        rt = 1._wp / r / PSC%internal_temp(jc,jk)
        tr2 = tr**2
        tr3 = tr**3
        pr2 = pr**2
        pr3 = pr**3

        !==========
        !     the h2so4/h2o pure solution concentration
        !==========
        xsb = 1._wp / (2._wp * (ksfit(3) + ksfit(4)/PSC%internal_temp(jc,jk)))         &
           * ( - ksfit(1) - ksfit(2)/PSC%internal_temp(jc,jk)                    &
               - sqrt( (ksfit(1) + ksfit(2)/PSC%internal_temp(jc,jk))**2 - 4._wp*(ksfit(3)   &
                       +ksfit(4)/PSC%internal_temp(jc,jk))                    &
              * (ksfit(5) + ksfit(6)/PSC%internal_temp(jc,jk) + ksfit(7)*logt-PSC%logpw(jc,jk)) ) )
        ! EMAC: bH2SO4b in routine mz_psc_liq_bH2SO4b (molality of H2SO4 in
        ! water for binary solution / (mol/kg)  )
        msb = 55.51_wp * xsb / (1._wp - xsb)
        !     ---
        IF ((PSC%internal_temp(jc,jk) <= 215._wp) .AND. (PSC%pn0(jc,jk) > 0._wp)) THEN
        !==========
        !     the hno3/h2so4/h2o solution composition
        !==========
          hsb = qsfit(1)                                                   &
           &  + qsfit(2)*tr2                                               &
           &  + (qsfit(3) + qsfit(4)*tr + qsfit(5)*tr2 + qsfit(6)*tr3)*pr  &
           &  + (qsfit(7) + qsfit(8)*tr + qsfit(9)*tr2)*pr2                &
           &  + qsfit(10)*tr*pr3
          hsb = exp(hsb)
          xnb = 1._wp / (2._wp * (knfit(3) + knfit(4)/PSC%internal_temp(jc,jk)))   &
           &  * ( - knfit(1) - knfit(2)/PSC%internal_temp(jc,jk)                   &
           &      - sqrt( (knfit(1) + knfit(2)/PSC%internal_temp(jc,jk))**2        &
           &              -4._wp*(knfit(3)+knfit(4)/PSC%internal_temp(jc,jk))      &
           &    * (knfit(5) + knfit(6)/PSC%internal_temp(jc,jk)                    &
           &    + knfit(7)*logt-PSC%logpw(jc,jk))))
          ! EMAC: bHNO3b in routine mz_psc_liq_bHNO3b (molality of HNO3 in
          ! water for binary solution / (mol/kg) )
          mnb = 55.51_wp * xnb / (1._wp - xnb)
          hnb = qnfit(1)                                                   &
           &  + qnfit(2)*tr2                                               &
           &  + (qnfit(3) + qnfit(4)*tr + qnfit(5)*tr2 + qnfit(6)*tr3)*pr  &
           &  + (qnfit(7) + qnfit(8)*tr + qnfit(9)*tr2)*pr2                &
           &  + qnfit(10)*tr*pr3
          hnb = exp(hnb)
          mnb2 = mnb**2
          msb2 = msb**2

          a  = (tt*hnb*mnb2 - tt*hsb*mnb*msb - 2._wp*mnb2*msb                      &
          &  + mnb*msb2 + hnb*mnb*msb*PSC%pn0(jc,jk)  - hsb*msb2*PSC%pn0(jc,jk))   &
          &  / (mnb2 - mnb * msb)

          b  = msb                                                       &
          &  * (-2._wp*tt*hnb*mnb+tt*hsb*msb+mnb*msb-hnb*msb*PSC%pn0(jc,jk)) &
          &  /  (mnb - msb)

          c  = (tt*hnb*mnb*msb2)   &
          &  / (mnb - msb)

          sqrtfunc = -2._wp*a**3 +9._wp*a*b-27._wp*c
          phi = atan(sqrt(4._wp * (a**2 - 3._wp* b)**3   &
              &          - sqrtfunc**2_wp)                     &
              &    / sqrtfunc )

          IF (phi < 0._wp) phi = phi + pi

          ! EMAC: bH2SO4 in routine mz_psc_liq_bH2SO4 (molality of H2SO4 in
          ! water / (mol/kg)  )
          PSC%ms(jc,jk) = - (a + 2._wp * sqrt(a**2 - 3._wp * b)              &
           &        * cos((pi + phi) / 3._wp)) / 3._wp

          IF (PSC%ms(jc,jk) < 0.0_wp) THEN
            PSC%ms(jc,jk) = msb / 1.e5_wp 
          END IF
          ! EMAC: bHNO3 in routine mz_psc_liq_bHNO3 (molality of HNO3 in
          ! water for ternary solution / (mol/kg), is wrongly documented in EMAC )
          PSC%mn(jc,jk) = mnb * (1._wp - PSC%ms(jc,jk) / msb)
          ! EMAC: wnen in routine mz_psc_liq_wnen (auxiliary variable, 
          !       = m(H2O)/m(H2O) + m(H2SO4)/m(H2O) + m(HNO3)/m(H2O)
          PSC%wnen(jc,jk) = 1._wp + PSC%ms(jc,jk) * 0.098076_wp + PSC%mn(jc,jk) * 0.063012_wp
          ! EMAC: wH2SO4 in routine mz_psc_liq_wH2SO4 (mass fraction of
          ! H2SO4 in the liquid aerosol / (kg/kg)  )
          PSC%ws(jc,jk) = PSC%ms(jc,jk) * 0.098076_wp / PSC%wnen(jc,jk)
          ! EMAC: wHNO3 in routine mz_psc_liq_wHNO3 (mass fraction of HNO3
          ! in the liquid aerosol / (kg/kg)  )
          PSC%wn(jc,jk) = PSC%mn(jc,jk) * 0.063012_wp / PSC%wnen(jc,jk)
          pn = PSC%mn(jc,jk) / ((hnb*PSC%mn(jc,jk) + hsb*PSC%ms(jc,jk))        &
                        / (PSC%mn(jc,jk)+PSC%ms(jc,jk)))
          ! EMAC: partHNO3 in routine mz_psc_liq_partHNO3 (fraction of HNO3
          ! remaining in gas phase after partitioning into liquid (1 = no
          ! removal from gas phase to the aerosol)  )
          PSC%parthno3(jc,jk) = 1._wp - (PSC%pn0(jc,jk) - pn) / PSC%pn0(jc,jk)
          IF (PSC%ms(jc,jk) > msb .OR. PSC%parthno3(jc,jk) > 1._wp) then
            PSC%ms(jc,jk) = msb
            PSC%mn(jc,jk) = 0._wp
            PSC%wnen(jc,jk) = 1._wp + PSC%ms(jc,jk) * 0.098076_wp
            PSC%ws(jc,jk) = msb * 0.098076_wp / PSC%wnen(jc,jk)
            PSC%wn(jc,jk) = 0._wp
            PSC%parthno3(jc,jk) = 1._wp
          END IF
        ELSE
        !==========
        !     assume solution is pure h2so4/h2o
        !==========
           PSC%ms(jc,jk) = msb
           PSC%mn(jc,jk) = 0._wp
           PSC%wnen(jc,jk) = 1._wp + PSC%ms(jc,jk) * 0.098076_wp
           PSC%ws(jc,jk) = msb * 0.098076_wp / PSC%wnen(jc,jk)
           PSC%wn(jc,jk) = 0._wp
           PSC%parthno3(jc,jk) = 1._wp
        END IF

        nsms = PSC%ns(jc,jk) / PSC%ms(jc,jk)
        sqrtws = sqrt(PSC%ws(jc,jk))
        sqrtwn = sqrt(PSC%wn(jc,jk))
        ws2 = PSC%ws(jc,jk)**2
        wn2 = PSC%wn(jc,jk)**2
        !==========
        !     the solubility of hcl and hbr h*hcl / h*hbr [mol/kg/atm]
        !     (adapted from luo et al., grl, 1995)
        !
        !     calculated concentrations assume that hcl and hbr are trace
        !     components of the aerosol
        !==========
        ! EMAC: hHCl in routine mz_psc_liq_hHCl; effective HCl henrys law
        ! constant in HNO3-H2SO4-H2O solution / (mol/(kg*atm))
        hhcl = - (21._wp + 46.61_wp*PSC%wn(jc,jk) + 4.069_wp*PSC%ws(jc,jk)           &
          &      - 4.837_wp*sqrtwn + 2.186_wp*sqrtws                         &
          &      - 63._wp*wn2 - 40.17_wp*PSC%wn(jc,jk)*PSC%ws(jc,jk) - 1.571_wp*ws2) &
          &    - 1._wp/PSC%internal_temp(jc,jk)                                  &
          &    * (-7437._wp - 8327.8_wp*PSC%wn(jc,jk) + 1300.9_wp*PSC%ws(jc,jk)      &
          &             + 1087.2_wp*sqrtwn - 242.71_wp*sqrtws                &
          &             + 18749._wp*wn2 + 18500._wp*PSC%wn(jc,jk)*PSC%ws(jc,jk)      &
          &             + 5632._wp*ws2)                                      &
          &    - log(PSC%wn(jc,jk) + 0.61_wp*PSC%ws(jc,jk))                          &
          &    - log(36.461_wp / (1000._wp*PSC%wnen(jc,jk)))
        hhcl = exp(hhcl) * 1.013e3_wp
        mcl = (rt * PSC%phcl0(jc,jk)) / (nsms + rt / hhcl)
        ! EMAC: wHCl in routine mz_psc_liq_wHCl; mass fraction of HCl in the
        ! liquid aerosol / (kg/kg) (seems different here but is actually the
        ! same as in EMAC)
        PSC%wcl(jc,jk) = mcl * 0.036461_wp / PSC%wnen(jc,jk)
        phcl = mcl / hhcl
        IF (PSC%phcl0(jc,jk) > 0.0_wp) THEN
          PSC%parthcl(jc,jk) =  phcl / PSC%phcl0(jc,jk)
        ELSE
          PSC%parthcl(jc,jk) = 0.0_wp
        END IF
        !     ---
        ! EMAC: hHBr in routine mz_psc_liq_hHBr;  effective HBr henrys law
        ! constant in HNO3-H2SO4-H2O solution / (mol/(kg*atm))
        PSC%hhbr(jc,jk) = - (17.83_wp + 1.02_wp*PSC%wn(jc,jk) - 1.08_wp*PSC%ws(jc,jk)         &
          &             + 3.9_wp*sqrtwn + 4.38_wp*sqrtws                          &
          &             - 8.87_wp*wn2 - 17._wp*PSC%wn(jc,jk)*PSC%ws(jc,jk) + 3.73_wp*ws2) &
          &           - 1._wp/PSC%internal_temp(jc,jk)                                &
          &           * (-8220.5_wp - 362.76_wp*PSC%wn(jc,jk) + 658.93_wp*PSC%ws(jc,jk)   &
          &                    - 914._wp*sqrtwn - 955.3_wp*sqrtws                 &
          &                    + 9976.6_wp*wn2 + 19778.5_wp*PSC%wn(jc,jk)*PSC%ws(jc,jk)   &
          &                    + 7680._wp*ws2)                                    &
          &           - log(PSC%wn(jc,jk) + 0.41_wp*PSC%ws(jc,jk))                        &
          &           - log(80.918_wp / (1000._wp*PSC%wnen(jc,jk)))
        PSC%hhbr(jc,jk) = exp(PSC%hhbr(jc,jk)) * 1.013e3_wp
        mbr = (rt * PSC%phbr0(jc,jk)) / (nsms + rt / PSC%hhbr(jc,jk))
        phbr = mbr / PSC%hhbr(jc,jk)
        ! EMAC: partHBr in routine mz_psc_liq_partHBr;  fraction of HBr
        ! remaining in gas phase after partitioning into liquid (1 = no
        ! removal from gas phase to aerosols) 
        IF (PSC%phbr0(jc,jk) > 0.0_wp) THEN
          PSC%parthbr(jc,jk) = phbr / PSC%phbr0(jc,jk)
        ELSE
           PSC%parthbr(jc,jk) = 0.0_wp
        END IF
        !==========
        !     the solubility of hocl and hobr h*hocl / h*hobr [mol/kg/atm]
        !     (from huthwelker et al., jac, 1995)
        !
        !     as an approximation, huthwelker et al. assumed that h* depends
        !     upon total molality in hno3/h2so4 solution
        !
        !     solubility of hocl is low enough to ignore gas phase removal
        !
        !     solubility of hobr is not known for all stratospheric conditions.
        !     limited data (hanson and ravishankara 1995 at 210 K, 60 wt% h2so4
        !     indicate that hhobr = approx. 18 * hhocl.
        !     for hobr an effective henry's law constant  = 18. * hhocl is used.
        !==========
        ! EMAC: hHOCl in routine mz_psc_liq_hHOCl; effective HOCl henrys law
        ! constant in HNO3-H2SO4-H2O solution / (mol/(kg*atm))
        PSC%hhocl(jc,jk) = 6.4946_wp                                                  &
          &                - (-0.04107_wp + 54.56_wp/PSC%internal_temp(jc,jk))        &
          &                    * (PSC%ms(jc,jk) + PSC%mn(jc,jk))                      &
          &                - 5862._wp * (1._wp/298.15_wp - 1._wp/PSC%internal_temp(jc,jk))
        PSC%hhocl(jc,jk) = exp(PSC%hhocl(jc,jk))
        mhocl = (rt * PSC%phocl0(jc,jk)) / (nsms + rt / PSC%hhocl(jc,jk))
        ! EMAC: wHOCl in routine mz_psc_liq_wHOCl; mass fraction of HOCl in the
        ! liquid aerosol / (kg/kg) (seems different here but is actually the
        ! same as in EMAC)
        PSC%whocl(jc,jk) = mhocl * 0.05246_wp / PSC%wnen(jc,jk)
        !     ---
        ! EMAC: hHOBr in routine mz_psc_liq_hHOBr; effective HOBr henrys law
        ! constant in HNO3-H2SO4-H2O solution / (mol/(kg*atm))
        PSC%hhobr(jc,jk) = 18._wp * PSC%hhocl(jc,jk)
        mhobr = (rt * PSC%phobr0(jc,jk)) / (nsms + rt / PSC%hhobr(jc,jk))
        ! EMAC: wHOBr in routine mz_psc_liq_wHOBr; mass fraction of HBr in the
        ! liquid aerosol / (kg/kg) (seems different here but is actually the
        ! same as in EMAC)
        PSC%whobr(jc,jk) = mhobr * 0.09691_wp / PSC%wnen(jc,jk)

!     ---


         ! remove liquid phase from gaseous phase (for HNO3 and H2O)

         ! Formula from EMAC (mz_psc_liq_H2O)
         PSC%H2O_Nconc_l(jc,jk) = MAX(MIN((1._wp - PSC%wn(jc,jk) - PSC%ws(jc,jk))     &
                           &          * p_tracer_now(jc,jk,jb,PSC%iTRH2SO4)   &
                           &          / PSC%ws(jc,jk)  * mol_weight_H2SO4 / amw,  &
                           &        PSC%H2O_Nconc_g(jc,jk)),                  &
                           &      0._wp)
         PSC%H2O_Nconc_g(jc,jk) = PSC%H2O_Nconc_g(jc,jk) - PSC%H2O_Nconc_l(jc,jk)

         PSC%HNO3_Nconc_l(jc,jk,jb) = (1.0_wp - PSC%parthno3(jc,jk)) * PSC%HNO3_Nconc_g(jc,jk)
         PSC%HNO3_Nconc_g(jc,jk) = PSC%HNO3_Nconc_g(jc,jk) - PSC%HNO3_Nconc_l(jc,jk,jb)
         
         ! here removing the other substances from gas phase are missing. Or
         ! not? Oh, probably not, because cgaml are calculated herein using the
         ! parts of the substances and
         ! concentration change is computed by the heterogeneous reaction rates
      END IF

   END DO
 END DO

     

  !=========
  !     the liquid phase diffusion constant [cm2/s]
  !     dl = diffusion constant of hocl in h2so4/hno3 solution
  !          this is used also for hbr and hobr.
  !
  !     dens = the solution density [g/cm3]
  !=========
   CALL aerodendi(PSC%iexception,PSC% internal_temp,    &
           &     PSC%ws, PSC%wn, PSC%ms, PSC%mn, PSC%density, PSC%diff,i_startidx,i_endidx,kstart,kend)
!    ---
  DO jk = kstart,kend
    DO jc = i_startidx, i_endidx
!     ---
      IF (PSC%iexception(jc,jk) == 0) THEN
!     ---
        dl = PSC%diff(jc,jk)
        dens = PSC%density(jc,jk)
        !==========
        !     convert effective henry's law constants to mol/dm3/atm
        !==========
        factor = dens / PSC%wnen(jc,jk)
        fhhbr  = PSC%hhbr(jc,jk)  * factor
        fhhocl = PSC%hhocl(jc,jk) * factor
        fhhobr = PSC%hhobr(jc,jk) * factor
        tn = PSC%internal_temp(jc,jk)

        ! =========
        ! Algorithm to get volume density (cm**3 / cm**3), surface area density (cm**2 /
        ! cm**3) and mean (= effective) radius (cm)
        ! based on Grainger et al. (1995) and Hervig and Deshler (1998)
        ! ======= 
        ! Originally the volume density was calculated by
        ! vliq = PSC%ns(jc,jk) * 98.076e-6 / PSC%ws(jc,jk) / dens
        ! The factor of 1e-6 actually converts the value of vliq from cm^3 / m^3 to
        ! m^3 / m^3. (For a value in cm you have to divide it by 100 to get it in m.
        ! In three dimensions the factor is 1e-6, of course.

        ! mass concentration in liquid aerosol (g / cm3)
        mtot = (PSC%HNO3_Nconc_l(jc,jk,jb)   *  63.012_wp            &
         &      +  p_tracer_now(jc,jk,jb,PSC%iTRH2SO4) * 98.076_wp   &
         &      +  PSC%H2O_Nconc_l(jc,jk) * amw   ) &
         &    / avo

        ! in cm**3 / cm**3
        vliq = mtot / dens
    
        ! parameters derived from  Hervig and Deshler (1998)
        ! in cm**2 / cm**3
        PSC%liqsur(jc,jk,jb) = 6.068e-8_wp*(1.e+12_wp * vliq)**0.671_wp
        ! in cm
        rmean = 3._wp * vliq / PSC%liqsur(jc,jk,jb)
        PSC%radius_STS(jc,jk,jb) = rmean / 100._wp


        !=======================================================================
        ! ROUTINE TO CALCULATE UPTAKE COEFFICIENTS (GAMMA VALUES).
        ! GAMMA VALUES ARE INDICATED BY VARIABLES WITH PREFIX 'G', FOR
        ! EXAMPLE GHOCLHCL IS THE GAMMA VALUE OF HOCL DUE TO REACTION WITH
        ! HCL IN THE DROPLETS.
        ! FROM THE GAMMA VALUES, SECOND ORDER RATE CONSTANTS ARE CALCULATED.
        ! THESE HAVE THE PREFIX 'R' AND HAVE UNITS CM3 MOLECULE-1 S-1. FOR
        ! EXAMPLE, THE LOSS OF CLNO3 AND HCL DUE TO THE HETEROGENEOUS REACTION
        ! CLNO3+HCL -> CL2+HNO3 IS D(CLNO3)/DT (UNITS MOLECULE CM-3 S-1) =
        ! -RCLNO3HCL.[CLNO3].[HCL], WHERE [CLNO3] AND [HCL] ARE THE
        ! ****TOTAL**** AMOUNTS OF THESE SPECIES IN UNITS MOLECULE CM-3.
        !=======================================================================
        !=========
        !     liquid phase concentration [mol/dm3/atm]
        !=========
        chclliq = PSC%wcl(jc,jk) * dens / 0.036461_wp
        !=========
        !     first order rate constants for reaction in liquid [dm3/mol/s]
        !=========
        k1hoclhcl = 1.e5_wp * chclliq
        k1hbrhocl = 1.e6_wp * PSC%whocl(jc,jk)* dens / 0.05246_wp
        k1hbrhobr = 1.e7_wp * PSC%whobr(jc,jk) * dens / 0.09691_wp
        k1hobrhcl = 1.e5_wp * chclliq
        ! =====================================================================
        ! ========================== HOCL + HCL ===============================
        ! ======================= ON LIQUID AEROSOL ===========================
        ! =====================================================================
        q = rmean * sqrt(k1hoclhcl/dl)
        fq = (q + 0.312_wp*q**2) / (3._wp + q + 0.312_wp*q**2)
        cbar = 2008._wp * PSC%sqrtt(jc,jk)
        IF(fhhocl*k1hoclhcl>0.0_wp) THEN
          PSC%cgaml(jc,jk,jb,ihs_HOCl_HCl) = MAX(fq / (fq + cbar &
              &  / (4._wp * fhhocl * 0.082_wp * tn * sqrt(k1hoclhcl*dl))), 1.e-12_wp)
        ELSE
          PSC%cgaml(jc,jk,jb,ihs_HOCl_HCl) = 1.e-12_wp
        ENDIF
        ! ====================================================================
        ! =========================== HOBR + HCL =============================
        ! ======================= ON LIQUID AEROSOL ==========================
        ! ====================================================================
        q = rmean * sqrt(k1hobrhcl/dl)
        fq = (q + 0.312_wp*q**2) / (3._wp + q + 0.312_wp*q**2) 
        cbar = 1477._wp * PSC%sqrtt(jc,jk)
        IF(fhhobr*k1hobrhcl>0.0_wp) THEN
          PSC%cgaml(jc,jk,jb,ihs_HOBr_HCl) = MAX(fq / (fq + cbar     &
              &  / (4._wp * fhhobr * 0.082_wp * tn * sqrt(k1hobrhcl*dl))), 1.e-12_wp)
        ELSE
          PSC%cgaml(jc,jk,jb,ihs_HOBr_HCl) = 1.e-12_wp
        ENDIF
        ! =====================================================================
        ! ========================= HBR + HOCL ================================
        ! ======================= ON LIQUID AEROSOL ===========================
        ! =====================================================================
        q = rmean * sqrt(k1hbrhocl/dl)
        fq = (q + 0.312_wp*q**2) / (3._wp + q + 0.312_wp*q**2) 
        cbar = 1616._wp * PSC%sqrtt(jc,jk)
        IF(fhhbr*k1hbrhocl>0.0_wp) THEN
          PSC%cgaml(jc,jk,jb,ihs_HOCl_HBr) = MAX(fq / (fq + cbar   &
              &   / (4._wp * fhhbr * 0.082_wp * tn * sqrt(k1hbrhocl*dl))),1.e-12_wp)
        ELSE
          PSC%cgaml(jc,jk,jb,ihs_HOCl_HBr) = 1.e-12_wp
        ENDIF
        ! =====================================================================
        ! ========================= HBR + HOBR ================================
        ! ======================= ON LIQUID AEROSOL ===========================
        ! =====================================================================
        q = rmean * sqrt(k1hbrhobr/dl)
        fq = (q + 0.312_wp*q**2) / (3._wp + q + 0.312_wp*q**2) 
        cbar = 1616._wp * PSC%sqrtt(jc,jk)
        IF(fhhbr*k1hbrhobr>0.0_wp) THEN
          PSC%cgaml(jc,jk,jb,ihs_HOBr_HBr) = MAX(fq / (fq + cbar     &
              &  / (4._wp * fhhbr * 0.082_wp * tn * sqrt(k1hbrhobr*dl))),1.e-12_wp)
        ELSE
          PSC%cgaml(jc,jk,jb,ihs_HOBr_HBr) = 1.e-12_wp
        ENDIF

        ! =====================================================================
        ! ====================== CLONO2+HCL -> CL2+HNO3 =======================
        ! ====================== CLONO2+H2O -> HOCL+HNO3 ======================
        ! ========================= ON LIQUID AEROSOL =========================
        ! =====================================================================
        ! TAKEN DIRECTLY FROM HANSON AND RAVISHANKARA, J. PHYS. CHEM.,
        ! 98, 5728, 1994, EXCEPT FOR THE FUNCTION F - SEE FOLLOWING COMMENTS -
        ! AND THE HCL SOLUBILITY, WHICH IS CALCULATED ACCORDING TO LUO ET AL.
        !     ---
        ! THE FUNCTION F: THE FORM OF F USED BY HANSON AND RAVISHANKARA CAN
        ! 'EXPLODE' UNDER CERTAIN CONDITIONS. IT HAS BEEN REPLACED HERE BY
        ! A STABLE FUNCTION THAT IS ACCURATE WITHIN ABOUT 4%. THIS IS ALSO
        ! THE CASE FOR OTHER REACTIONS THAT FOLLOW.

        ! water activity
        ! remark: in EMAC 1013.25 is air pressure, but in Hanson and Ravishankara (1994), 
        !         1013.25 is used for stratospheric conditions
        ah2o = 1013.25_wp * PSC%pw(jc,jk)                                              &
               &  / (10_wp**(9.217_wp-2190.0_wp/(PSC%internal_temp(jc,jk)-12.70_wp)))
        g0=1.18e-4_wp+9.1e-3_wp*ah2o+0.50_wp*ah2o**2
        gs=ah2o*ksur*chclliq
        p=rho*chclliq/ah2o
        gcalc=g0*sqrt(1.0_wp+p)
        q = rmean * sqrt(ah2o) / 1.4e-6_wp 
        fq = (q + 0.312_wp*q**2) / (3._wp + q + 0.312_wp*q**2)
        ge=1.0_wp/(1.0_wp/(gs+fq*gcalc)+1.0_wp/alpha)
        PSC%cgaml(jc,jk,jb,ihs_ClNO3_HCl) = MAX(ge * (gs+fq*gcalc*p/(1._wp+p))   &
               &                         / (gs+fq*gcalc),1.e-12_wp)
        PSC%cgaml(jc,jk,jb,ihs_ClNO3_H2O) = MAX(ge - PSC%cgaml(jc,jk,jb,ihs_ClNO3_HCl),1.e-12_wp)
        ! =====================================================================
        ! ======================== BrONO2 + H2O ===============================
        ! ====================== ON LIQUID AEROSOL ============================
        ! =====================================================================
        ! FROM HANSON ET AL., JGR, 101, 9063-9069, 1996.
        ! =====================================================================
        grxn = 211._wp * ah2o**1.37_wp
        PSC%cgaml(jc,jk,jb,ihs_BrNO3_H2O) = MAX((0.84_wp * grxn) / (grxn + 0.84_wp),1.e-12_wp)

        PSC%cgaml(jc,jk,jb,ihs_N2O5_H2O) = 0.1_wp

!    ---
      ELSE
        PSC%radius_STS(jc,jk,jb) = 0._wp  !setting of radius_STS necessary for Restart
      END IF
!    ---
    END DO
  END DO

END SUBROUTINE art_psc_aerosol

!///////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////


SUBROUTINE aerodendi(iexception, blkta, &
                    ws, wn, ms, mn, density, diff, i_startidx, i_endidx,kstart, kend)
!<
! SUBROUTINE aerodendi
! This subroutine calculates diffusion constant of hocl in h2so4/hno3 solution
! and its solution density
! Based on:
! - Carslaw et al. (1995)
! - Luo et al. (1996)
! Part of Module: mo_art_psc_sts
! Author: Michael Weimer, KIT
! Initial Release: 2017-02-06
! Modifications:
!>
  IMPLICIT NONE
  INTEGER, INTENT(in) ::      &
     &  i_startidx, i_endidx, &  !< loop indices
     &  kstart, kend
  INTEGER, DIMENSION(:,:), INTENT(IN) ::  &
     &   iexception              !< exception handler
  REAL(wp), INTENT(in) :: &
    &  blkta(:,:)                !< temperature (K)                   
  REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) ::    &
    &  ws, wn, ms, mn            !< see in previous routine
  
  REAL (KIND=wp), DIMENSION(:,:), INTENT(OUT) ::   &
    &   density, diff            !< density (g / cm3)  and diffusion coefficient (cm2 / s) 
                                 !  of the ternary solution
  
  ! local variables
  INTEGER :: &
    &   jc,jk  !< loop indices
  
  ! fitting parameters
  REAL (KIND=wp) ::                               &
    &   c, yexp, a, xn, xb, t0, viss, visn, visc, &
    &   w, wh, v1, vs, vn, vmcal
  
  REAL (KIND=wp), DIMENSION(21), PARAMETER ::                                &
    &   xden = (/                                                            &
    &           2.393284e-02_wp,  -4.359335e-05_wp,      7.961181e-08_wp,    &
    &           -0.198716351_wp, 1.39564574e-03_wp,     -2.020633e-06_wp,    &
    &             0.51684706_wp,    -3.0539e-03_wp,      4.505475e-06_wp,    &
    &            -0.30119511_wp,   1.840408e-03_wp, -2.7221253742e-06_wp,    &
    &         -0.11331674116_wp,    8.47763e-04_wp,   -1.22336185e-06_wp,    &
    &              0.3455282_wp,    -2.2111e-03_wp,   3.503768245e-06_wp,    &
    &             -0.2315332_wp,    1.60074e-03_wp,    -2.5827835e-06_wp/)
  REAL (KIND=wp), DIMENSION(12), PARAMETER ::                                        &      
    &   xdiff = (/   623.8082_wp ,   5.221606_wp, -8.085769e-02_wp, 2.1769575e-4_wp, &
    &                154.3466_wp , -0.9521694_wp,-2.6749929e-03_wp, 1.984055e-4_wp,  &
    &                  0.02656_wp,  0.0019710_wp,  0.00023760e0_wp, 0._wp/)

  !=======================================================================
  !     density of ternary solution in g/cm3
  !
  !     ws, wn are wt fraction
  !     fitted to 0.05 < ws + wn < 0.70 wt fraction, but exprapolates well
  !
  !     temperature  > 185 K
  !=======================================================================

  ! EMAC: content of mz_psc_density, output dens in routine
  ! mz_psc_surface_liquid; mass density of liquid solution calculated
  ! according to Luo et al. (1996)

  ! This soubroutine is called much later in EMAC than here (--> possible
  ! difference??)
  DO jk = kstart,kend
    DO jc = i_startidx, i_endidx

      IF (iexception(jc,jk) == 0) THEN
        w = ws(jc,jk) + wn(jc,jk)
        wh = 1._wp - w
        ! These parametrisations are most probably based on the book by Soehnel and
        ! Novotny: Densities of aqueous solutions of inorganic
        ! substances, Elsevier, Amsterdam, 1985
        ! It cannot really be found in the book but the parameters are fitted
        ! by Luo et al. (1996) based on the values by Soehnel and Novotny.
        v1 = func(xden( 1), xden( 2), xden( 3), 0._wp, blkta(jc,jk))
        vs = func(xden( 4), xden( 5), xden( 6), 0._wp, blkta(jc,jk))      &
           + func(xden( 7), xden( 8), xden( 9), 0._wp, blkta(jc,jk))*w    &
           + func(xden(10), xden(11), xden(12), 0._wp, blkta(jc,jk))*w**2
        vn = func(xden(13), xden(14), xden(15), 0._wp, blkta(jc,jk))      &
           + func(xden(16), xden(17), xden(18), 0._wp, blkta(jc,jk))*w    &
           + func(xden(19), xden(20), xden(21), 0._wp, blkta(jc,jk))*w**2
        vmcal = wh*v1/18.016_wp + vs*ws(jc,jk)/98.08_wp + vn*wn(jc,jk)/63.016_wp
        density(jc,jk) = 0.001_wp / vmcal
      END IF
    END DO
  END DO

  !=======================================================================
  !     c   : wt % of h2so4
  !     visc: viscosity in si units (as mole ratio of h2so4 and hno3)
  !     t   : temperature in k
  !     diff: diffusivity of hocl in cm2 s-1
  !           calculated using the houghten cubic cell model
  !           with cell dimension 3.65 angstroms
  !           (see luo et al., grl, 1994)
  !           note this is in good agreement (+-10% with a composition
  !           dependent cell dimension as given in huthwelker et al., jac,
  !           1995)
  !=======================================================================

  ! EMAC: content of mz_psc_diff (also called much later there)
  DO jk = kstart,kend
    DO jc = i_startidx, i_endidx

      IF (iexception(jc,jk) == 0) THEN
        c = ws(jc,jk) * 100._wp
        yexp = -7.722133_wp + 3.773159e-2_wp * c
        a = exp(yexp)
        xn = 0.518604_wp
        xb = func(xdiff( 1), xdiff( 2), xdiff( 3), xdiff( 4), c)
        t0 = func(xdiff( 5), xdiff( 6), xdiff( 7), xdiff( 8), c)
        viss = funcviss(a, xn, xb, t0, blkta(jc,jk))
        a = func(xdiff( 9), xdiff(10), xdiff(11), xdiff(12), c)
        xn= -0.012750_wp*mn(jc,jk)
        xb = 735.7_wp
        t0 = 92.89_wp + 0.68480_wp*mn(jc,jk)
        visn = funcviss(a, xn, xb, t0, blkta(jc,jk))
        visc = (ms(jc,jk)*viss + mn(jc,jk)*visn) / (ms(jc,jk)+mn(jc,jk))
        diff(jc,jk) = 3.37e-14_wp * blkta(jc,jk) * density(jc,jk) / visc
      END IF
    END DO
  END DO
END SUBROUTINE aerodendi

FUNCTION func(k1,k2,k3,k4,c)
  IMPLICIT NONE
  REAL (KIND=wp), INTENT(IN) :: k1,k2,k3,k4,c
  REAL (KIND=wp) :: func

  func = k1 + k2*c + k3*c**2 + k4*c**3
END FUNCTION func

FUNCTION funcviss(k1,k2,k3,k4,c)
  IMPLICIT NONE
  REAL (KIND=wp), INTENT(IN) :: k1,k2,k3,k4,c
  REAL (KIND=wp) :: funcviss

 funcviss = 1.e-3_wp * k1 * c**k2 * exp( k3/(c-k4))
END FUNCTION funcviss

END MODULE mo_art_psc_sts
