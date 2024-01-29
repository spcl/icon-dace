!
! mo_art_psc_state
! This module provides routines for PSCs in ICON-ART
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

MODULE mo_art_psc_state

! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish, message, message_text, warning
  USE mo_physical_constants,            ONLY: argas, avo
  USE mo_math_constants,                ONLY: pi
! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_impl_constants,            ONLY: IART_QV, IART_QI, IART_QC
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_aerosol_utilities,         ONLY: art_air_properties
  USE mo_art_psc_types,                 ONLY: t_art_psc,       &
                                          &   Pa2torr,         &
                                          &   mol_weight_HNO3, &
                                          &   mol_weight_NAT,  &
                                          &   rho_NAT,         &
                                          &   pi4_3_rhoNAT    
  USE mo_art_psc_sts,                   ONLY: art_psc_aerosol
  USE mo_art_psc_sedimentation,         ONLY: art_psc_simpleUpwind_sedi,   &
                                          &   art_psc_KinPar_NAT_sedi_vel
  USE mo_art_psc_chemistry,             ONLY: art_psc_calc_k_het
  

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: art_psc_main
  PUBLIC :: art_psc_ice_formation
  ! These routines are set to public, rather for the standalone version, and
  ! not for the ICON-ART setup
  PUBLIC :: art_psc_calc_HNO3_gls
  PUBLIC :: art_psc_kin_NAT_param
  PUBLIC :: art_psc_NAT_formation
  PUBLIC :: art_calc_N_and_r_NAT
  PUBLIC :: art_psc_recombine_HNO3
  
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_psc_main(PSC,p_dtime,temp,pres,dz, cell_area,   &
                        jg, p_tracer_now)
!<
! SUBROUTINE art_psc_main
! Calculates everything of PSCs: formation, heterogeneous reaction rates, sedimentation etc.
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-09-14
! Modifications: 
! 2017-06-09: Michael Weimer, KIT
! - added initialisation of kinetic NAT parametrisation structure
!>
  IMPLICIT NONE
  TYPE(t_art_psc), INTENT(inout) :: &
    &  PSC                   !< structure with PSC meta information
  REAL(wp), INTENT(in) :: &
    &  p_dtime,           &  !< model time step (s)
    &  temp(:,:,:),       &  !< air temperature (K)
    &  pres(:,:,:),       &  !< air pressure (Pa)
    &  dz(:,:,:),         &  !< model layer height (m)
    &  cell_area(:,:)        !< area of grid cells (m2)
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_now(:,:,:,:) !< tracer number concentrations (# / cm3)
  INTEGER, INTENT(in) :: &
    &  jg                    !< patch id
  ! local variables
  INTEGER ::                  &
     &  i_startidx, i_endidx, &
     &  jb, n                !< loop indices
   TYPE(t_art_atmo), POINTER :: &
     &  art_atmo

   art_atmo => p_art_data(jg)%atmo


  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    CALL art_psc_calc_water_Nconc(PSC, p_art_data(jg)%chem%water_tracers,   &
                  &               IART_QV, IART_QC, IART_QI, jb)

    CALL art_psc_calc_HNO3_gls(PSC, p_tracer_now, jb, i_startidx,i_endidx)

    IF (PSC%NSB > 1) THEN
      CALL art_psc_kin_NAT_param(PSC%NSB,PSC,temp,pres,p_dtime,  &
                      &          p_tracer_now,       &
                      &          jb,i_startidx,i_endidx)


    ELSE
      CALL art_psc_NAT_formation(PSC,temp,               &
                    &            jb,i_startidx,i_endidx)

      CALL art_calc_N_and_r_NAT(PSC%radius_NAT(:,:,jb,1),              &
                        &       PSC%dens_NAT(:,:,jb,1),                &
                        &       PSC%HNO3_Nconc_s(:,:,jb),              &
                        &       PSC%total_max_density,                 &
                        &       PSC%r_min,                             &
                        &       PSC%kstart, PSC%kend, i_startidx, i_endidx)

    END IF


    CALL art_psc_aerosol(PSC,          &
                 &       jb,           &
                 &       PSC%kstart,   &
                 &       PSC%kend,     &
                 &       i_startidx,   &
                 &       i_endidx,     &
                 &       temp,         &
                 &       pres,         &
                 &       p_tracer_now)


    CALL art_psc_calc_k_het(PSC, p_tracer_now,           &
                   &        temp, pres,                  &
                   &        jb,i_startidx,i_endidx)

    ! calculate dynamic viscosity of the air
    CALL art_air_properties(pres(:,:,jb),temp(:,:,jb),                     &
                &           i_startidx,i_endidx,1,art_atmo%nlev,           &
                &           p_art_data(jg)%air_prop%art_free_path(:,:,jb), &
                &           p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb) )
!                &           i_endidx,1,art_atmo%nlev,jb,p_art_data(jg))

    IF (PSC%NSB > 1) THEN

      DO n = 1, PSC%NSB
        CALL art_psc_KinPar_NAT_sedi_vel(PSC%v_sed_NAT(i_startidx:i_endidx,:,jb,n), &
                          &     temp(:,:,jb),                                       &
                          &     pres(:,:,jb),                                       &
                          &     p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb),       &
                          &     PSC%radius_NAT(:,:,jb,n),                           &
                          &     PSC%dens_NAT(:,:,jb,n),                             &
                          &     i_startidx, i_endidx, PSC%kstart, art_atmo%nlev)

        CALL art_psc_simpleUpwind_sedi(p_tracer_now(:,:,jb,PSC%tracer_indices_bins(n)), &
                          &      PSC%v_sed_NAT(:,:,jb,n),                               &
                          &      dz,                                                    &
                          &      p_dtime, jb, PSC%kstart,                               &
                          &      PSC%kend, i_startidx, i_endidx,                        &
                          &      PSC%NAT_sedi_rel_diff(:,:,n),                          &
                          &      cell_area)

        PSC%v_sed_NAT_out(:,1:art_atmo%nlev,jb,n) = PSC%v_sed_NAT(:,1:art_atmo%nlev,jb,n)
      END DO
    ELSE
      CALL art_psc_KinPar_NAT_sedi_vel(PSC%v_sed_NAT(i_startidx:i_endidx,:,jb,1), &
                        &     temp(:,:,jb),                                       &
                        &     pres(:,:,jb),                                       &
                        &     p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb),       &
                        &     PSC%radius_NAT(:,:,jb,1),                           &
                        &     PSC%dens_NAT(:,:,jb,1),                             &
                        &     i_startidx, i_endidx, PSC%kstart, art_atmo%nlev)


        CALL art_psc_simpleUpwind_sedi(PSC%HNO3_Nconc_s(:,:,jb),     &
                          &      PSC%v_sed_NAT(:,:,jb,1),            &
                          &      dz,                                 &
                          &      p_dtime, jb, PSC%kstart,            &
                          &      PSC%kend, i_startidx, i_endidx,     &
                          &      PSC%NAT_sedi_rel_diff(:,:,1),       &
                          &      cell_area)


        PSC%v_sed_NAT_out(:,1:art_atmo%nlev,jb,1) = PSC%v_sed_NAT(:,1:art_atmo%nlev,jb,1)

    END IF


    CALL art_psc_recombine_HNO3(PSC,p_tracer_now,jb,art_atmo%nlev,i_startidx,i_endidx)
  END DO


END SUBROUTINE art_psc_main
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_psc_calc_water_Nconc(PSC, water_tracers, iqv, iqc, iqi, jb)
!<
! SUBROUTINE art_psc_calc_water_Nconc
! This subroutine transfers the water substances of ICON into the PSC structure
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-06-09
! Modifications:
!>
  IMPLICIT NONE
  TYPE(t_art_psc), INTENT(inout) ::  &
      &  PSC    !< PSC meta data structure
  REAL(wp), INTENT(in), TARGET :: &
      &  water_tracers(:,:,:,:)
  INTEGER, INTENT(in) :: &
      &  iqv, iqc, iqi  !< indices of gaseous, liquid and solid water in water_tracers
  INTEGER, INTENT(in) ::   &
      &  jb     !< block index

  PSC%H2O_Nconc_g => water_tracers(:,:,jb,iqv)

  PSC%H2O_Nconc_s => water_tracers(:,:,jb,iqi) 

  PSC%H2O_Nconc_l => water_tracers(:,:,jb,iqc) 
END SUBROUTINE art_psc_calc_water_Nconc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_psc_calc_HNO3_gls(PSC, p_tracer_now, jb,i_startidx,i_endidx)
!<
! SUBROUTINE art_psc_calc_HNO3_gls
! This subroutine transfers HNO3 in ICON-ART to the PSC module
! Output as number concentration (# / cm3)
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-06-30
! Modifications: 
!>
  IMPLICIT NONE
  TYPE(t_art_psc), INTENT(inout) ::  &
    &  PSC    !< PSC meta data structure
  REAL(wp), INTENT(in) :: &
    &  p_tracer_now(:,:,:,:)    !< number concentrationsof tracers (# / cm3)
  INTEGER, INTENT(in) ::   &
    &  jb,i_startidx,i_endidx   !< block index
  !local variables
  INTEGER :: jc, jk, n  !< loop indices

  IF (PSC%NSB < 2) THEN
    DO jk = PSC%kstart,PSC%kend
      DO jc = i_startidx, i_endidx
        PSC%HNO3_Nconc_g(jc,jk) = p_tracer_now(jc,jk,jb,PSC%iTRHNO3)
        PSC%HNO3_Nconc_l(jc,jk,jb) = 0._wp
        PSC%HNO3_Nconc_s(jc,jk,jb) = 0._wp
      END DO
    END DO
  ELSE
    DO jk = PSC%kstart,PSC%kend
      DO jc = i_startidx, i_endidx
        PSC%HNO3_Nconc_g(jc,jk) = p_tracer_now(jc,jk,jb,PSC%iTRHNO3)
        PSC%HNO3_Nconc_l(jc,jk,jb) = 0._wp
        PSC%HNO3_Nconc_s(jc,jk,jb) = 0._wp
      END DO
    END DO

    DO n = 1, PSC%NSB
      DO jk = PSC%kstart,PSC%kend
        DO jc = i_startidx, i_endidx
          PSC%HNO3_Nconc_s(jc,jk,jb) = PSC%HNO3_Nconc_s(jc,jk,jb) &
              &                  + p_tracer_now(jc,jk,jb,PSC%tracer_indices_bins(n))
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE art_psc_calc_HNO3_gls
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_psc_ice_formation(PSC,temp,pres,   &
               &                 jb,i_startidx,i_endidx)
!<
! SUBROUTINE art_psc_ice formation
! This subroutine calculates the ice formation in polar
! stratosphere
! Based on: Marti and Mauersberger, GRL, 1993
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-02-27
! Modifications: 
!>
  IMPLICIT NONE
  TYPE(t_art_psc), INTENT(inout) ::  &
    &  PSC             !< PSC meta data structure
  REAL(wp), INTENT(in) :: &
    &  temp(:,:,:),       & !< air temperature (K)
    &  pres(:,:,:)          !< air pressure (Pa)
  INTEGER, INTENT(in) ::              &
    &  jb, i_startidx,i_endidx    !< loop indices 

  !local variables
  INTEGER :: &
    &  jc, jk   !< loop indices
  REAL(wp), PARAMETER :: &
    &  A = -2663.5_wp,   & !< parameters of the parametrisation
    &  B = 12.537_wp       !< 
  REAL(wp) ::            &
    &  log_ice_pressure, & !< logarithm of saturation vapour pressure over ice
    &  partial_pres_H2O    !< partial pressure of total water (Pa)

  DO jk = PSC%kstart,PSC%kend
    DO jc = i_startidx,i_endidx
      log_ice_pressure = 0._wp

!      IF ((temp(jc,jk,jb) < 250._wp) .AND. (temp(jc,jk,jb) > 170._wp)) THEN
        log_ice_pressure = A / (temp(jc,jk,jb)) + B

        partial_pres_H2O = PSC%H2O_Nconc_g(jc,jk) * 1.e6_wp * argas * temp(jc,jk,jb) / avo

        ! zero is lower limit for ice content and gaseous water is upper limit
        PSC%ice_vmr_Marti(jc,jk,jb) = MIN(MAX((partial_pres_H2O - 10._wp ** log_ice_pressure)   &
                &              / pres(jc,jk,jb), 0._wp),                                        &
                &              PSC%H2O_Nconc_g(jc,jk) * 1.e6_wp * argas * temp(jc,jk,jb)        &
                &              / pres(jc,jk,jb) / avo)

        ! For comparing ICON microphysics with Marti and Mauersberger, this line
        ! is removed! It should be included again if you want to reintegrate it
        ! again within the PSC module
!      PSC%H2O_Nconc_g(jc,jk) = PSC%H2O_Nconc_g(jc,jk) - PSC%ice_vmr_Marti(jc,jk,jb) * pres * avo &
!      / argas / temp * 1.e-6_wp

!      END IF
    END DO
  END DO
END SUBROUTINE art_psc_ice_formation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_psc_NAT_formation(PSC,temp,       &
                   &             jb,i_startidx,i_endidx)
!<
! SUBROUTINE art_psc_NAT_formation
! This subroutine calculates the NAT formation in polar
! stratosphere (diagnostically)
! Based on: Hanson and Mauersberger, GRL, 1988
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-03-06
! Modifications: 
!>
  IMPLICIT NONE
  TYPE(t_art_psc), INTENT(inout) :: &
    &  PSC             !< PSC meta data structure
  REAL(wp), INTENT(in), DIMENSION(:,:,:) :: &
    &  temp            !< air temperature (K)
  INTEGER, INTENT(in) :: &
    &  jb,i_startidx, i_endidx    !< loop indices
  ! local variables
  REAL(wp) :: &
    &  partial_pres_HNO3_Pa  !< partial pressure of HNO3 in air (Pa)
  REAL(wp) :: &
    &  local_temp,        &  !< temperature but set to 180 K if temp is lower
    &  HNO3ice_coex_pres, &  !< saturation vapour pressure of HNO3 over NAT (Pa)
    &  HNO3_t                !< total HNO3 number concentration (# / cm3)
  INTEGER :: &
    & jc, jk                 !< loop indices


  DO jk = PSC%kstart,PSC%kend
    DO jc = i_startidx,i_endidx
      HNO3ice_coex_pres = 0._wp

      IF (temp(jc,jk,jb) < 200._wp) THEN
        local_temp = MAX(temp(jc,jk,jb),180._wp)

        HNO3_t = PSC%HNO3_Nconc_g(jc,jk)
        partial_pres_HNO3_Pa = HNO3_t  * argas * temp(jc,jk,jb)  / avo * 1.e6_wp

        HNO3ice_coex_pres = art_press_HNO3_over_NAT(local_temp,  &
                     &           PSC%H2O_Nconc_g(jc,jk))

        PSC%HNO3_Nconc_s(jc,jk,jb) =  MIN(MAX((partial_pres_HNO3_Pa - HNO3ice_coex_pres)  &
                                  &        / argas / temp(jc,jk,jb) * avo / 1.e6_wp, 0._wp),  &
                                  &      PSC%HNO3_Nconc_g(jc,jk))
        PSC%HNO3_Nconc_g(jc,jk) = PSC%HNO3_Nconc_g(jc,jk) - PSC%HNO3_Nconc_s(jc,jk,jb)

      END IF
    END DO
  END DO

END SUBROUTINE art_psc_NAT_formation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION art_press_HNO3_over_NAT(temp,H2O_gl)
!<
! FUNCTION art_press_HNO3_over_NAT
! This function calculates the saturation vapour pressure (in Pa) of HNO3 over
! NAT particles
! valid only in the temperature range 180 to 200 K
! Based on: Hanson and Mauersberger, GRL, 1988
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-06-14
! Modifications: 
!>
  IMPLICIT NONE
  REAL(wp), INTENT(in) ::   &
    &  temp,                &  !< air temperature (K)
    &  H2O_gl                  !< water number concentration including vapour and liquid (# / cm3)
  REAL(wp) ::   &
    &  art_press_HNO3_over_NAT !< output: saturation vapour pressure of HNO3 over NAT (Pa)
  !local variables
  REAL(wp) ::                    &
    &  m, b,                     & !< parameters of Hanson and Mauersberger
    &  partial_pres_water_torr,  & !< partial pressure of water (torr)
    &  log_HNO3ice_coex_pres       !< logarithm (base 10) of HNO3 saturation vapour pressure (torr)

  m = -2.7836_wp - 0.00088_wp * temp
  b = 38.9855_wp - 11397.0_wp / temp + 0.009179_wp * temp

  partial_pres_water_torr = H2O_gl * argas * temp  / avo * 1.e6_wp * Pa2torr

  log_HNO3ice_coex_pres = m * LOG10(partial_pres_water_torr) + b

  art_press_HNO3_over_NAT = 10._wp ** log_HNO3ice_coex_pres / Pa2torr
  
END FUNCTION art_press_HNO3_over_NAT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_psc_kin_NAT_param(NSB,PSC, temp, pres, dtime, p_tracer_now,   &
               &                 jb, i_startidx, i_endidx)
!<
! SUBROUTINE art_psc_kin_NAT_param
! This subroutine calculates NAT for given bins and exchanges between the bins
! (mass consistent)
! Based on: 
! - Carslaw et al., JGR, 2002
! - van den Broek et al., ACP, 2004
! - Kirner et al., ACP, 2011
! - EMAC subroutine K_PARA_NAT (used in msbm_physc) (of version 2.52)
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-06-14
! Modifications: 
!>
  IMPLICIT NONE
  
  INTEGER, INTENT(in) :: &
    &  NSB  !< number of size bins
  TYPE(t_art_psc), INTENT(inout) ::  &
    &  PSC  !< PSC meta data structure
  INTEGER, INTENT(in) ::   &
    &  jb, i_startidx, i_endidx  !< loop indices
  REAL(wp), INTENT(in) ::  &
    &  temp(:,:,:),        &  !< air temperature (K)
    &  pres(:,:,:)            !< air pressure (hPa)
  REAL(wp), INTENT(inout) ::   &
    &  p_tracer_now(:,:,:,:)     !< tracer number concentrations (#/cm3)
  REAL(wp), INTENT(in) :: &
    &  dtime                 !< advective model time step (s)
  !local variables
  REAL(wp) ::        &
    &  pHNO3,        & !< partial pressure of HNO3 (Pa)
    &  SVAPN,        & !< parameter based on vapour pressures of HNO3 and NAT to
                       !  decide whether NAT initially should be formed 
                       !  (yes if SVAPN > 1)
    &  m_nat,        & !< current mass of NAT (g)
    &  m_nat_old,    & !< temporarily saved mass of NAT (g)
    &  dens_diff,    & !< difference of current NAT number density
                       !  to limit number density of the bin (1 / cm3)
    &  new_mass_conc,& !< mass concentration of NAT after growth (g / cm3)
    &  nat_rad,      & !< radius of NAT (m)
    &  dens_nat_spheres!< number density of NAT  with current radius
                       !  (not NAT molecules) (mol / cm3)

  REAL(wp), DIMENSION(NSB) :: &
    &  mass_conc_diff  !< NAT mass concentration of the part that is larger than than the
                       !  limit of NAT number density (g / cm3)

  REAL(wp) ::            &
   &  rad_sq,            & !< radius squared (m2)
   &  diff_nat_Nconc,    & !< number concentrationo f difference between old 
                           !  and new particle mass (# / cm3)
   &  growth_factor,     & !< growth factor for increasing the radius temporarily 
                           !  (Carslaw et al., 2002) (m2 / s)
   &  total_nat_Nconc_save
  INTEGER ::  &
   &  jk, jc, n  !< loop indices

  DO jk = PSC%kstart,PSC%kend
    DO jc = i_startidx, i_endidx
   
      pHNO3 = PSC%HNO3_Nconc_g(jc,jk) * argas * temp(jc,jk,jb) / avo * 1e6_wp
      SVAPN = pHNO3 / art_press_HNO3_over_NAT(temp(jc,jk,jb) - PSC%NatFormThreshold, &
                                   &          PSC%H2O_Nconc_g(jc,jk))

      DO n = 1, NSB
        ! number density of NAT molecules in n'th bin (mol / cm3)
        dens_nat_spheres = p_tracer_now(jc,jk,jb,PSC%tracer_indices_bins(n)) / avo

        mass_conc_diff(n) = 0.0_wp ! initialization of mass to be transferred to next size bin
        dens_diff         = 0.0_wp ! initialization of xs particle no.density in size bin (n)
        nat_rad           = PSC%rbin_av(n)
        ! mass per particle [g] ; van den Broek et al. 1871-2
        m_nat_old         = pi4_3_rhoNAT * nat_rad**3.0_wp  

        IF (dens_nat_spheres < 0.0_wp) THEN
          dens_nat_spheres      = 0.0_wp
          PSC%dens_NAT(jc,jk,jb,n) = 0.0_wp
        ELSE
          ! calculates no.density of NAT molecules (# particles / cm3)
          ! (multiply number density of molecules in NAT by molar weight of NAT and divide 
          !  by particle mass of NAT particle with radius rbin_av)
          PSC%dens_NAT(jc,jk,jb,n) = dens_nat_spheres * mol_weight_NAT / m_nat_old 
        END IF

        ! Initialization of NAT particles at start of simulation
        IF ((n == 1) .AND. (PSC%dens_NAT(jc,jk,jb,n) < 1.5e-5_wp) .AND. (SVAPN >= 1._wp)) THEN
          PSC%dens_NAT(jc,jk,jb,n) = 1.5e-5_wp ! 1/cm**3
        END IF

        ! Growth factor can get negative if HNO3 partial pressure is below
        ! saturation vapour pressure wrt NAT. With this, NAT can also descrease
        ! in size. For temperatures greater than 220 K growth factor is set to
        ! nat_rad^2 / (2 * dtime) so that it is depleted completely in this grid
        ! box
        growth_factor = art_psc_NAT_growth_factor(nat_rad,               &
                                     &            pHNO3,                 &
                                     &            PSC%H2O_Nconc_g(jc,jk),&
                                     &            temp(jc,jk,jb),        &
                                     &            pres(jc,jk,jb),        &
                                     &            dtime)

        !--------------------------------
        ! In case that the size decrease is larger than the initial radius then the
        ! particle has completely evaporated
        !---------------------------------

        rad_sq = nat_rad**2.0_wp + 2.0_wp*growth_factor*dtime
        IF (rad_sq >= 0._wp) THEN
          nat_rad = sqrt(rad_sq)
        ELSE 
          nat_rad = 0.0_wp
        END IF
   
        ! determine new particle mass after growth
        ! calculates mass per particle [g] for new radius
        ! here the particle no. density remains the same, i.e. all particles
        ! experience identical radius change
        m_nat = pi4_3_rhoNAT * nat_rad**3.0_wp


        !----------------------------------------------------------------
        ! Each size bin is 'refilled' in the next section
        !----------------------------------------------------------------

        ! updates the number of particles (# particles / cm3)

        ! calculate mass concentration (of NAT molecules in all particles) after growth (g / cm3)
        new_mass_conc = PSC%dens_NAT(jc,jk,jb,n) * m_nat

        !----------------------------------------------------------------
        ! Transfer mass from the previous size bin if necessary
        !----------------------------------------------------------------

        IF (n > 1) THEN
          IF (mass_conc_diff(n-1) > 0.0_wp) THEN
             ! if mass has to be added from previous size bin: add it to the
             ! mass concentration and recalculate the number densities as above
             new_mass_conc = new_mass_conc + mass_conc_diff(n-1)
          END IF
        END IF


        ! same formula as above (but after growth and mass transfer)
        dens_nat_spheres = new_mass_conc / mol_weight_NAT
        ! number density of NAT particles after growth (dens_nat weighted with
        ! mass ratio before and after growth) (# particles / cm3)
        PSC%dens_NAT(jc,jk,jb,n) = new_mass_conc / m_nat_old


        !------------------------------------------------------------------------------------------
        ! Checks to see if the particle number density exceeds the threshold set
        ! for each size bin
        ! If so, the xs particles are essentially stored and added, as mass, to the
        ! next size bin up.
        ! No growth of stored particles occurs till they are added to the next bin.
        !------------------------------------------------------------------------------------------

        IF (PSC%dens_NAT(jc,jk,jb,n) > PSC%no_density_limit(n)) THEN
          ! how many particles over the limit
          dens_diff = PSC%dens_NAT(jc,jk,jb,n) - PSC%no_density_limit(n) 
   
          IF (dens_diff > 1.e-30_wp) THEN! filter to ensure -ve NAT concentrations don't occur
            ! calculate mass concentration which is above the bin limit number
            ! concentration and set number densities of NAT molecules and
            ! particles to the limit itself
            mass_conc_diff(n) = dens_diff * m_nat_old  ! [g / cm3]
            ! modify no. of NAT particles in bin (n)
            dens_nat_spheres = (PSC%no_density_limit(n) * m_nat_old) / mol_weight_NAT
            PSC%dens_NAT(jc,jk,jb,n) = PSC%no_density_limit(n)
          ELSE
            mass_conc_diff(n) = 0.0_wp
          END IF
        END IF

        p_tracer_now(jc,jk,jb,PSC%tracer_indices_bins(n)) = dens_nat_spheres * avo

        PSC%dens_NAT(jc,jk,jb,n) = PSC%dens_NAT(jc,jk,jb,n) * 1.e6_wp ! convert back to m**-3
      END DO ! loop over 9 bins

      !-------------------------------------------------------------------------------------------
      ! Warning given if an overflow of the last size bin occurs (during
      ! prolonged low temps)
      !-------------------------------------------------------------------------------------------

      IF (PSC%dens_NAT(jc,jk,jb,NSB) > PSC%no_density_limit(NSB) * 1.e6_wp) THEN
        WRITE(message_text,'(A, 3I5)') 'last BIN overflowing, jc, jk, jb ',  &
                &                     jc, jk, jb
        CALL warning('mo_art_psc_state:art_psc_kin_NAT_param',  &
               &     message_text)
      END IF
   
      !----------------------------------------------------------------------------------------
      ! Warning if the maximum particle radius becomes too big (i.e) greater than
      ! 25e-6 m
      !----------------------------------------------------------------------------------------
   
      IF (nat_rad > PSC%rbin_max(NSB)) then
        CALL message('mo_art_psc_state:art_psc_kin_NAT_param',  &
               &     'NAT particle is growing greater than maximum radius of last BIN.')
      END IF
   
      DO n = 1, PSC%NSB
        IF (PSC%dens_NAT(jc,jk,jb,n) > 5.0e-7_wp) THEN
           PSC%radius_NAT(jc,jk,jb,n)  = PSC%rbin_av(n)
        ELSE
           PSC%radius_NAT(jc,jk,jb,n) = 0.0_wp
        ENDIF
      END DO


      total_nat_Nconc_save = PSC%HNO3_Nconc_s(jc,jk,jb)
      PSC%HNO3_Nconc_s(jc,jk,jb) = 0._wp
      DO n = 1,NSB
        PSC%HNO3_Nconc_s(jc,jk,jb) = PSC%HNO3_Nconc_s(jc,jk,jb)        &
               &                   + p_tracer_now(jc,jk,jb,PSC%tracer_indices_bins(n))
      END DO

      diff_nat_Nconc = PSC%HNO3_Nconc_s(jc,jk,jb) - total_nat_Nconc_save

      PSC%HNO3_Nconc_g(jc,jk) = PSC%HNO3_Nconc_g(jc,jk) - diff_nat_Nconc

      IF (PSC%HNO3_Nconc_g(jc,jk) < 0._wp) THEN
!        WRITE(*,*) 'jc, jk', jc, jk
!        WRITE(*,*) 'tot_Nconc', PSC%HNO3_Nconc_g(jc,jk)
!        WRITE(*,*) 'diff_nat_Nconc', diff_nat_Nconc
        CALL finish('mo_art_psc_state:art_psc_kin_NAT_param',  &
               &    'HNO3_Nconc_g is lower than zero!')
      END IF

!      PSC%H2O_Nconc_g(jc,jk)    = PSC%H2O_Nconc_g(jc,jk) - 3._wp * diff_nat_Nconc

      IF (PSC%H2O_Nconc_g(jc,jk) < 0._wp) THEN
        CALL finish('mo_art_psc_state:art_psc_kin_NAT_param',  &
               &    'H2O_Nconc_g is lower than zero!')
      END IF
    END DO
  END DO
END SUBROUTINE art_psc_kin_NAT_param
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION art_psc_NAT_growth_factor(nat_rad,pHNO3,H2O_gl,temp,pres,dtime)
!<
! FUNCTION art_psc_NAT_growth_factor
! This function calculates the growth factor of NAT particles depending on
! diffusion coefficient and vapour pressures of HNO3 and temperature 
! Based on: 
! - Hall and Pruppacher, 1976, Journal of the Atmospheric Sciences (used for
!                                                                   diffusion coefficients)
! - Carslaw et al., JGR, 2002
! - EMAC function growth_factor (used in msbm_physc) (of version 2.52)
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-06-14
! Modifications: 
!>
  IMPLICIT NONE
  REAL(wp), INTENT(in) :: &
    &  nat_rad,           & !< radius of NAT particles (m)
    &  pHNO3,             & !< partial pressure of HNO3 (Pa)
    &  H2O_gl,            & !< number concentration of gaseous + liquid water (# / cm3)
    &  dtime,             & !< model time step (s)
    &  temp,              & !< air temperature (K)
    &  pres                 !< air pressure (Pa)
  REAL(wp) :: &
    &  art_psc_NAT_growth_factor !< growth factor (m2 / s)
  ! local variables
  REAL(wp) ::               &
    &  diff_H2O, diff_HNO3, & !< diffusion coefficients of H2O and HNO3 in air (m2 / s) 
    &  vHNO3,               & !< mean molecular speed of HNO3 in air (m / s)
    &  pHNO3_over_NAT,      & !< saturation vapour pressure of HNO3 over NAT (Pa)
    &  temp_corr,           & !< corrected temperature if not in range of art_press_HNO3_over_NAT 
                              !  (in K)
    &  D_star_HNO3            !< diffusion coefficient of HNO3 accouting for
                              !< "mass transfer noncontinuum effects for
                              !  particles with sizes similar to the mean free path (m2 / s)

  IF (temp > 220._wp) THEN
    art_psc_NAT_growth_factor = - nat_rad**2._wp / (2._wp * dtime)
  ELSE
    temp_corr = MAX(temp, 180._wp)

    ! from Hall and Pruppacher (1976)
    diff_H2O = 0.22_wp*((temp/273.15_wp)**1.94_wp)*(101325._wp/pres*100._wp) * 1.e-4_wp


    ! transfer diff_H2O to that of HNO3. The factor relies on Zhu et al., JAMS, 2015
    diff_HNO3 = diff_H2O * 0.466_wp

    
    ! mean molecular speed of HNO3 in the air
    ! (expected value of Maxwell-Boltzmann distribution)
    vHNO3 = sqrt((8.0_wp*argas*temp)/(pi*mol_weight_HNO3*1.0e-3_wp))
    ! should be only gasous water
    pHNO3_over_NAT = art_press_HNO3_over_NAT(temp_corr,H2O_gl)

    ! diff. coeff. of HNO3 in air, [m**2/s] accounting for mass transfer
    ! continuum effect (when size ~ mean free path)
    D_star_HNO3 = diff_HNO3/(1.0_wp + 4.0_wp * diff_HNO3 / (vhno3 * nat_rad))

    art_psc_NAT_growth_factor  &
      &  = D_star_HNO3*mol_weight_HNO3/(rho_NAT*temp*argas)*(pHNO3-pHNO3_over_NAT)
  END IF
END FUNCTION art_psc_NAT_growth_factor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_psc_recombine_HNO3(PSC,p_tracer_now,jb, nlev, i_startidx, i_endidx)
!<
! SUBROUTINE art_psc_recombine_HNO3
! Recombine the different phases of HNO3 to gaseous phase 
! and set the phases to zero which are computed diagnostically
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-07-03
! Modifications: 
!>
  IMPLICIT NONE
  TYPE(t_art_psc), INTENT(inout) :: &
    &  PSC                   !< PSC meta data structure
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_now(:,:,:,:) !< tracer number concentrations (# / cm3)
  INTEGER, INTENT(in) :: &
    &  jb, nlev,         &
    &  i_startidx, i_endidx  !< loop indices
  ! local variables
  INTEGER ::  &
    &  jc, jk   !< loop indices

  IF (PSC%NSB > 1) THEN
    ! In case of the kinetic NAT parametrisation everything concerning NAT PSCs
    ! is calculated within the NAT bin tracers themselves, so only the
    ! diagnostically calculated liquid part in STS comes back to the HNO3 tracer
    DO jk = PSC%kstart,nlev
      DO jc = i_startidx,i_endidx
        p_tracer_now(jc,jk,jb,PSC%iTRHNO3) = PSC%HNO3_Nconc_g(jc,jk) + PSC%HNO3_Nconc_l(jc,jk,jb)
        PSC%HNO3_Nconc_g(jc,jk) = p_tracer_now(jc,jk,jb,PSC%iTRHNO3)
        PSC%HNO3_Nconc_l(jc,jk,jb) = 0._wp
      END DO
    END DO
  ELSE
    ! in case of the diagnostic NAT parametrisation, the solid part (i.e. in
    ! NAT) is calculated diagnostically, too. However, sedimentation is
    ! calculated for the solid part, only. So, all parts have to be summed up to
    ! get the HNO3 tracer
    DO jk = PSC%kstart,nlev
      DO jc = i_startidx,i_endidx
        p_tracer_now(jc,jk,jb,PSC%iTRHNO3) = PSC%HNO3_Nconc_g(jc,jk)       &
                 &                           + PSC%HNO3_Nconc_l(jc,jk,jb)  &
                 &                           + PSC%HNO3_Nconc_s(jc,jk,jb)

        PSC%HNO3_Nconc_g(jc,jk) = p_tracer_now(jc,jk,jb,PSC%iTRHNO3)
      END DO
    END DO
  END IF
END SUBROUTINE art_psc_recombine_HNO3

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_calc_N_and_r_NAT(radius_NAT, dens_NAT,                  &
                   &            NAT_Nconc, total_max_density, r_min,   &
                   &            kstart, kend, i_startidx, i_endidx)
!<
! SUBROUTINE art_calc_massfrac_and_nc
! Calculate mass fraction of solid phase and number concentration of solid
! particles (NAT) from number concentration (#/cm3).
! Part of Module: mo_art_psc_state
! Author: Michael Weimer, KIT
! Initial Release: 2017-09-15
! Modifications: 
!>
  IMPLICIT NONE
  REAL(wp), INTENT(out) ::  &
    &  radius_NAT(:,:),     &  !< radius of NAT particles (m)
    &  dens_NAT(:,:)           !< number density of NAT particles (m-3)
  REAL(wp), INTENT(in) :: &
    &  NAT_Nconc(:,:),    &    !< number concentration of NAT (i.e. NAT phase HNO3) (#/cm3)
    &  total_max_density, &    !< maximum number density of NAT (cm-3)
    &  r_min                   !< minimum radius of NAT particles (m)
  INTEGER, INTENT(in) :: &
    &  kstart, kend, i_startidx, i_endidx  !< loop indices
  !local variables
  REAL(wp) ::              &
    &  N_solid_max,        &  !< maximum number density of NAT (m-3)
    &  molec_per_particle     !< number of molecules per particle
  INTEGER :: &
    &  jc,jk   !< loop indices

  N_solid_max = total_max_density * 1.e6_wp

  radius_NAT(:,:) = 0._wp
  dens_NAT(:,:) = 0._wp

  molec_per_particle = pi4_3_rhoNAT * r_min**3 / mol_weight_NAT * avo

  DO jk = kstart,kend
    DO jc = i_startidx, i_endidx
      dens_NAT(jc,jk) = NAT_Nconc(jc,jk) * 1.e6_wp / molec_per_particle

      IF (dens_NAT(jc,jk) >= N_solid_max) THEN

        radius_NAT(jc,jk) = (dens_NAT(jc,jk) / N_solid_max) ** (1._wp / 3._wp) &
                &         * r_min 

        dens_NAT(jc,jk) = N_solid_max

      ELSE
        radius_NAT(jc,jk) = r_min
      END IF

    END DO
  END DO

END SUBROUTINE art_calc_N_and_r_NAT

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
END MODULE mo_art_psc_state
