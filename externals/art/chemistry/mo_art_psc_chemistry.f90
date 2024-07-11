!
! mo_art_psc_chemistry
! This module provides routines for calculation of the heterogeneous reaction
! rates on polar stratospheric clouds
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

MODULE mo_art_psc_chemistry

! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_parallel_config,               ONLY: nproma
  USE mo_run_config,                    ONLY: iqi
  USE mo_physical_constants,            ONLY: argas, avo, amw
  USE mo_math_constants,                ONLY: pi, sqrt2
! ART
  USE mo_art_psc_types,                 ONLY: t_art_psc
  USE messy_mecca_kpp_global,           ONLY: ihs_N2O5_H2O , ihs_HOCl_HCl, ihs_ClNO3_HCl, &
                                          &   ihs_ClNO3_H2O, ihs_N2O5_HCl, ihs_ClNO3_HBr, &
                                          &   ihs_BrNO3_HCl, ihs_HOCl_HBr, ihs_HOBr_HCl , &
                                          &   ihs_HOBr_HBr , ihs_BrNO3_H2O
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: art_psc_calc_k_het

  
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_psc_calc_k_het(PSC,p_tracer_now,temp,pres,jb,i_startidx,i_endidx)
!<
! SUBROUTINE art_psc_calc_k_het
! This subroutine calculates the heterogeneous reaction rates in the surface of PSCs
! Based on: 
! Part of Module: mo_art_psc_chemistry
! Author: Michael Weimer, KIT
! Initial Release: 2017-02-06
! Modifications: 
!>
  IMPLICIT NONE
  TYPE(t_art_psc), INTENT(inout) :: &
    &  PSC                            !< PSC meta data structure
  REAL(wp), INTENT(in) ::          &
    &  temp(:,:,:),                &  !< air temperature (K)
    &  pres(:,:,:),                &  !< air pressure (Pa)
    &  p_tracer_now(:,:,:,:)          !< tracer number concentration (#/cm3)
  INTEGER, INTENT(in) :: &
    &  jb,i_startidx,i_endidx !< loop indices
  !local variables
  INTEGER :: jc, jk     !< loop indices

  ! calculate general values that are necessary for all reaction rates
  DO jk = PSC%kstart, PSC%kend
    DO jc = i_startidx, i_endidx
      PSC%v_th(jc,jk) = SQRT(8._wp * temp(jc,jk,jb) * argas / pi * 1.e3_wp)
      PSC%mean_free_path(jc,jk) = 1._wp / (sqrt2 * pi * pres(jc,jk,jb)   &
             &                  / argas / temp(jc,jk,jb) * avo * 3.7e-10_wp**2)
    END DO
  END DO

  ! calculate radius and number density of ice
  CALL art_psc_calc_N_r_ice(PSC%dens_ice(:,:,jb),PSC%radius_ice(:,:,jb),  &
              &             p_tracer_now(:,:,jb,iqi),                     &
              &             PSC%kstart,PSC%kend,                          &
              &             i_startidx,i_endidx)


  !##################################################
  ! heterogeneous reactions with HCl
  !##################################################


  IF (PSC%iTRHCl /= 0) THEN
    ! number concentration of HCl
    DO jk = PSC%kstart, PSC%kend
      DO jc = i_startidx, i_endidx
        PSC%Nconc_higher(jc,jk) = p_tracer_now(jc,jk,jb,PSC%iTRHCl)
      END DO
    END DO

    ! ClONO2 + HCl --> HNO3 + Cl2
    IF (PSC%iTRClONO2 /= 0) THEN
      CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_ClNO3_HCl),PSC%kstart,PSC%kend,                 &
                &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                       &
                &            PSC%Nconc_higher,97.458_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_ClNO3_HCl), &
                &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                               &
                &            0.3_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                  &
                &            0.2_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
    END IF


    ! HOCl + HCl --> H2O + Cl2
    IF (PSC%iTRHOCl /= 0) THEN
      CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_HOCl_HCl),PSC%kstart,PSC%kend,                  &
                &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                       &
                &            PSC%Nconc_higher,52.460_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_HOCl_HCl),  &
                &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                               &
                &            0.2_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                  &
                &            0.1_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
    END IF

    ! HOBr + HCl --> H2O + Br2
    IF (PSC%iTRHOBr /= 0) THEN
      CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_HOBr_HCl),PSC%kstart,PSC%kend,                  &
                &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                       &
                &            PSC%Nconc_higher,96.911_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_HOBr_HCl),  &
                &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                               &
                &            0.3_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                  &
                &            0.1_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
    END IF

    ! N2O5 + HCl --> HNO3 + ClONO
    IF (PSC%iTRN2O5 /= 0) THEN
      CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_N2O5_HCl),PSC%kstart,PSC%kend,                  &
                &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                       &
                &            PSC%Nconc_higher,108.0104_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_N2O5_HCl),&
                &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                               &
                &            0.03_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                 &
                &            0.003_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
    END IF
    
    ! BrONO2 + HCl --> HNO3 + BrCl
    IF (PSC%iTRBrONO2 /= 0) THEN
      CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_BrNO3_HCl),PSC%kstart,PSC%kend,                 &
                &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                       &
                &            PSC%Nconc_higher,141.909_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_BrNO3_HCl),&
                &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                               &
                &            0.26_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                 &
                &            0.3_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
    END IF
    
  END IF  ! iTRHCl /= 0


  !###################################
  !  heterogeneous reactions with water
  !###################################

  ! number concentration of water vapour
  DO jk = PSC%kstart, PSC%kend
    DO jc = i_startidx, i_endidx
      PSC%Nconc_higher(jc,jk) = PSC%H2O_Nconc_g(jc,jk)
    END DO
  END DO

  ! ClONO2 + H2O --> HNO3 + HOCl
  IF (PSC%iTRClONO2 /= 0) THEN
    CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_ClNO3_H2O),PSC%kstart,PSC%kend,                 &
              &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                       &
              &            PSC%Nconc_higher,97.458_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_ClNO3_H2O), &
              &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                               &
              &            0.3_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                  &
              &            0.004_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
  END IF


  ! N2O5 + H2O --> 2 HNO3 
  IF (PSC%iTRN2O5 /= 0) THEN
    CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_N2O5_H2O),PSC%kstart,PSC%kend,                   &
              &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                        &
              &            PSC%Nconc_higher,108.0104_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_N2O5_H2O), &
              &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                                &
              &            0.027_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                 &
              &            0.0004_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
  END IF
  
  ! BrONO2 + H2O --> HNO3 + BrCl
  IF (PSC%iTRBrONO2 /= 0) THEN
    CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_BrNO3_H2O),PSC%kstart,PSC%kend,                  &
              &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                        &
              &            PSC%Nconc_higher,141.909_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_BrNO3_H2O), &
              &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                                &
              &            0.26_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                  &
              &            0.001_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
  END IF


  !###############################################
  ! heterogeneous reactions with HBr
  !###############################################

  IF (PSC%iTRHBr /= 0) THEN
    ! number concentration of HBr
    DO jk = PSC%kstart, PSC%kend
      DO jc = i_startidx, i_endidx
        PSC%Nconc_higher(jc,jk) = p_tracer_now(jc,jk,jb,PSC%iTRHBr)
      END DO
    END DO

    
    ! HOCl + HBr --> H2O + BrCl
    IF (PSC%iTRHOCl /= 0) THEN
      CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_HOCl_HBr),PSC%kstart,PSC%kend,                  &
                &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                       &
                &            PSC%Nconc_higher,52.460_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_HOCl_HBr),  &
                &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                               &
                &            0.3_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                  &
                &            0.3_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
    END IF


    ! HOBr + HBr --> H2O + Br2
    IF (PSC%iTRHOBr /= 0) THEN
      CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_HOBr_HBr),PSC%kstart,PSC%kend,                  &
                &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                       &
                &            PSC%Nconc_higher,96.911_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_HOBr_HBr),  &
                &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                               &
                &            0.1_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                  &
                &            0.1_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
    END IF

    
    ! ClONO2 + HBr --> HNO3 + BrCl
    IF (PSC%iTRClONO2 /= 0) THEN
      CALL art_psc_sum_k_het(PSC%k_het(:,:,jb,ihs_ClNO3_HBr), PSC%kstart,PSC%kend,                &
                &            i_startidx,i_endidx,PSC%khet_max,PSC%khet_min,                       &
                &            PSC%Nconc_higher,97.458_wp,PSC%v_th,PSC%cgaml(:,:,jb,ihs_ClNO3_HBr), &
                &            PSC%liqsur(:,:,jb),PSC%mean_free_path,                               &
                &            0.3_wp,PSC%radius_ice(:,:,jb),PSC%dens_ice(:,:,jb),                  &
                &            0.3_wp,PSC%radius_NAT(:,:,jb,:),PSC%dens_NAT(:,:,jb,:))
    END IF
    
  END IF  ! iTRHBr /= 0


END SUBROUTINE art_psc_calc_k_het
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_psc_sum_k_het(sum_khet,kstart,kend,i_startidx,i_endidx,khet_max,khet_min,     &
                &            Nconc1,M2,v_th,cgaml,liqsur,mean_free_path,                     &
                &            gamma_ice,radius_ice,dens_ice,                                  &
                &            gamma_NAT,radius_NAT,dens_NAT)
!<
! SUBROUTINE art_psc_sum_k_het
! This subroutine calculates sum of heterogeneous reaction rates on STS, NAT and ice PSCs
! Based on: 
! Part of Module: mo_art_psc_chemistry
! Author: Michael Weimer, KIT
! Initial Release: 2017-02-06
! Modifications: 
!>
  IMPLICIT NONE
  REAL(wp), INTENT(out) :: &
    &  sum_khet(:,:)          !< sum of the heterogeneous reaction rates (cm3 s-1 molec-1)

  INTEGER, INTENT(in)  ::  &
    &  kstart, kend,       &  !< loop indices
    &  i_startidx,         &
    &  i_endidx

  REAL(wp), INTENT(in) ::        &
    &  Nconc1(:,:),              & !< number concentration of the most abundant compound (cm-3)
    &  M2,                       & !< molar mass of the less abundant one (g mol-1)
    &  v_th(:,:),                & !< mean molecular speed without molar weight (m/s*sqrt(g/mol))
    &  cgaml(:,:),               & !< uptake coefficients in STS (-)
    &  liqsur(:,:),              & !< liquid surface area density (cm2 cm-3)
    &  mean_free_path(:,:),      & !< mean free path of the molecules in air (m)
    &  gamma_ice,                & !< reaction probability on ice (-)
    &  radius_ice(:,:),          & !< radius of ice particles (m)
    &  dens_ice(:,:),            & !< number density of ice particles (m-3)
    &  gamma_NAT,                & !< reaction probability on NAT (-)
    &  khet_max,                 & !< maximum and minimum value for
    &  khet_min,                 & !  heterogeneous reaction rate 
    &  radius_NAT(:,:,:),        & !< radius of NAT particles (m)
    &  dens_NAT(:,:,:)             !< number density of NAT particles (m-3)

  ! local variables
  INTEGER ::    &
    &  jc, jk     !< loop indices

  

  DO jk = kstart,kend
    DO jc = i_startidx,i_endidx
      ! heterogeneous rate on STS
      sum_khet(jc,jk) = art_psc_k_het_liq(v_th(jc,jk),cgaml(jc,jk),liqsur(jc,jk),  &
                           &              Nconc1(jc,jk),M2,khet_max,khet_min)

      ! heterogeneous rate on ice
      sum_khet(jc,jk) = sum_khet(jc,jk)  &
          &    + art_psc_k_het_ice(v_th(jc,jk),mean_free_path(jc,jk),gamma_ice, &
                    &              radius_ice(jc,jk), dens_ice(jc,jk),          &
                    &              Nconc1(jc,jk), M2, khet_max,khet_min)
      ! heterogeneous rate on NAT
      sum_khet(jc,jk) = sum_khet(jc,jk)  &
          &    + art_psc_k_het_NAT(v_th(jc,jk),mean_free_path(jc,jk),gamma_NAT, &
                    &              radius_NAT(jc,jk,:), dens_NAT(jc,jk,:),      &
                    &              Nconc1(jc,jk), M2, khet_max,khet_min)
    END DO
  END DO

END SUBROUTINE art_psc_sum_k_het
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION art_psc_k_het_liq(v_th,cgaml,liqsur,Nconc1,M2,khet_max,khet_min)
!<
! FUNCTION art_psc_k_het_liq
! This function calculates the heterogeneous reaction rate on liquid (STS)
! surface: k_het = 0.25 * v_th * gamma * liqsurdens / Nconc
! Based on: Hanson et al., JGR, 1996
! Part of Module: mo_art_psc_chemistry
! Author: Michael Weimer, KIT
! Initial Release: 2018-01-22
! Modifications: 
!>
  IMPLICIT NONE
  REAL(wp)  :: &
    &  art_psc_k_het_liq    !< heterogeneous reaction rate normalized with the
                            !  concentration of the most abundant species
                            !  (cm3 s-1 molec-1)
  REAL(wp), INTENT(in) :: &
    &  v_th,              & !< thermodynamic velocity of a particle in air 
                            !  without the molar weight (m / s * sqrt(g / mol))
    &  cgaml,             & !< uptake coefficient of the reaction (-)
    &  liqsur,            & !< liquid surface area density (cm2 cm-3)
    &  Nconc1,            & !< number concentration of the most abundant reactant (molec cm-3)
    &  M2,                & !< molar weights of the other reactant (g mol-1)
    &  khet_max, khet_min   !< boundaries for reaction rate to stay realistc 
                            !  (can explode if input number concentration is too low)
  ! local variables
  REAL(wp) :: &
    &  k_het_save

  k_het_save = 0.25_wp * v_th * cgaml * liqsur * 100._wp

  ! the "+ 1" ensures the denominator not to get zero
  art_psc_k_het_liq = k_het_save / SQRT(M2) / (Nconc1 + 1._wp)
  art_psc_k_het_liq = MAX(MIN(art_psc_k_het_liq,khet_max),khet_min)
END FUNCTION art_psc_k_het_liq
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION art_psc_k_het_ice(v_th,mean_free_path,gamma_react,radius,  &
                  &        Nconc_aero,Nconc1,M2,khet_max,khet_min)
!<
! FUNCTION art_psc_k_het_ice
! This function calculates the heterogeneous reaction rate on ice
! surface
! Based on: Rolf Mueller, PhD thesis, 1994 and references therein
! Part of Module: mo_art_psc_chemistry
! Author: Michael Weimer, KIT
! Initial Release: 2018-01-22
! Modifications: 
!>
  IMPLICIT NONE
  REAL(wp)  :: &
    &  art_psc_k_het_ice    !< heterogeneous reaction rate normalized with the
                            !  concentration of the most abundant species
                            !  (cm3 s-1 molec-1)
  REAL(wp), INTENT(in) :: &
    &  v_th,              & !< thermodynamic velocity of a particle in air 
                            !  without the molar weight (m / s * sqrt(g / mol))
    &  mean_free_path,    & !< mean free path in the air (m)
    &  gamma_react,       & !< reaction probability (-)
    &  radius,            & !< radius of the particles (m)
    &  Nconc_aero,        & !< number concentration of the particles (m-3)
    &  Nconc1,            & !< number concentration of the most abundant reactant (cm-3) 
                            !  ATTENTION ON UNITS!!!
    &  M2,                & !< molar weight of the other reactant (g mol-1)
    &  khet_max, khet_min   !< boundaries for reaction rate to stay realistc 
                            !  (can explode if input number concentration is too low)
  !local variables
  REAL(wp)   ::  &
    &  k_het_save

  k_het_save = gamma_react * pi * v_th  * radius**2 * Nconc_aero    &
      &        / (1._wp + (3._wp * gamma_react * radius)            &
      &                    / (4._wp * mean_free_path))

  ! the "+ 1" ensures the denominator not to get zero
  art_psc_k_het_ice = k_het_save / SQRT(M2) / (Nconc1 + 1._wp)
 
  art_psc_k_het_ice = MAX(MIN(art_psc_k_het_ice,khet_max),khet_min)
END FUNCTION art_psc_k_het_ice
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION art_psc_k_het_NAT(v_th,mean_free_path,gamma_react,radii,    &
                    &      NAT_densities,Nconc1,M2, khet_max, khet_min)
!<
! FUNCTION art_psc_k_het_NAT_kinpar
! This function calculates the heterogeneous reaction rate on NAT in case of the
! kinetic NAT parametrisation
! Based on: Rolf Mueller, PhD thesis, 1994 and references therein
! Part of Module: mo_art_psc_chemistry
! Author: Michael Weimer, KIT
! Initial Release: 2018-01-22
! Modifications: 
!>
  IMPLICIT NONE
  REAL(wp)  :: &
    &  art_psc_k_het_NAT  !< heterogeneous reaction rate normalized with the
                          !  concentration of the most abundant species
                          !  (cm3 s-1 molec-1)

  REAL(wp), INTENT(in) :: &
    &  v_th,              & !< thermodynamic velocity of a particle in air 
                            !  without the molar weight (m / s * sqrt(g / mol))
    &  mean_free_path,    & !< mean free path in the air (m)
    &  gamma_react,       & !< reaction probability (-)
    &  radii(:),          & !< radii of the NAT bins (m)
    &  NAT_densities(:),  & !< number concentration of the NAT bins (m-3)
    &  Nconc1,            & !< number concentration of the most abundant reactant (cm-3) 
                            !  !ATTENTION ON UNITS!!!
    &  M2,                & !< molar weight of the other reactant (g mol-1)
    &  khet_max, khet_min   !< boundaries for reaction rate to stay realistc 
                            !  (can explode if input number concentration is too low)
  ! local variables
  REAL(wp)   ::  &
    &  factor_khet, R

  factor_khet = gamma_react * pi * v_th

  R = SUM(radii(:)**2 * NAT_densities(:)              &
    &     / (1._wp + (3._wp * gamma_react * radii(:)) &
    &                 / (4._wp * mean_free_path)))

  ! the "+ 1" ensures the denominator not to get zero
  art_psc_k_het_NAT = factor_khet / SQRT(M2) / (Nconc1 + 1._wp) * R

  art_psc_k_het_NAT = MAX(MIN(art_psc_k_het_NAT,khet_max),khet_min)
END FUNCTION art_psc_k_het_NAT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_psc_calc_N_r_ice(N, r, qi, kstart,kend, i_startidx, i_endidx)
!<
! SUBROUTINE art_psc_calc_N_r_ice
! Calculate number concentration and radius of ice with the same formulas as in
! the one moment scheme of ICON in operational use
! Part of Module: mo_art_psc_chemistry
! Author: Michael Weimer, KIT
! Initial Release: 2017-02-02
! Modifications: 
!>
  REAL(wp), INTENT(out) :: &
    &  N(:,:),   &  !< particle number concentration (1 / m3)
    &  r(:,:)       !< particle radius (m)
  REAL(wp), INTENT(in) :: &
    &  qi(:,:)             !< ice number concentration (# / cm3)

  INTEGER, INTENT(in) :: &
    &  kstart,kend,      & !< loop indices
    &  i_startidx,       &
    &  i_endidx

  !local variables
  REAL(wp) ::  &
    &  m_i            !< particle mass (kg)
  INTEGER :: &
    &  jc, jk         !< loop indices

  DO jk = kstart, kend
    DO jc = i_startidx,i_endidx
      ! considering that the number concentration in the temperature range of ice
      ! formation in the stratosphere (.lt. -40 degC) is always the minimum value in
      ! the ICON formula, we use here just this minimum value of 250.e3 m-3
      N(jc,jk) = 250.e3_wp
      m_i = MAX(qi(jc,jk), 0._wp) * 1.e6_wp * amw * 1e-3_wp / avo / N(jc,jk)
      r(jc,jk) = (m_i / 130._wp)**(1._wp / 3._wp) / 2._wp
    END DO
  END DO

END SUBROUTINE art_psc_calc_N_r_ice
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
END MODULE mo_art_psc_chemistry
