!
! mo_art_psc_sedimentation
! This module provides routines for PSC sedimentation
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

MODULE mo_art_psc_sedimentation

! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_physical_constants,            ONLY: grav, argas, avo
!ART
  USE mo_art_psc_types,                 ONLY: rho_NAT,       &
                                          &   mol_weight_NAT
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: art_psc_KinPar_NAT_sedi_vel
  PUBLIC :: art_psc_simpleUpwind_sedi

CONTAINS

!!
!!-------------------------------------------------------------------------
!!

SUBROUTINE art_psc_KinPar_NAT_sedi_vel(v_sed_NAT, temp, pres, dyn_visc,  &
                     &                 r_NAT, dens_NAT, &
                     &                 istart, iend, kstart, nlev)
!<
! SUBROUTINE art_psc_KinPar_NAT_sedi_vel
! Routine for calculation of sedimentation velocity of NAT particles
! Based on:
! - Carslaw et al., J. Geophys. Res., 2002
! Part of Module: mo_art_psc_sedimentation
! Author: Michael Weimer, KIT
! Initial Release: 2017-09-18
! Modifications:
!>
  IMPLICIT NONE
  REAL(wp), INTENT(in) ::  &
    &  temp(:,:),          &  !< air temperature (K)
    &  pres(:,:),          &  !< air pressure (Pa)
    &  dyn_visc(:,:),      &  !< dynamic viscosity of the air (kg m-1 s-1)
    &  r_NAT(:,:),         &  !< NAT radius (m)
    &  dens_NAT(:,:)          !< number density of NAT (m-3)
  INTEGER, INTENT(in) ::   &
    &  istart, iend,       &  !< loop indices
    &  kstart, nlev
  REAL(wp), DIMENSION(istart:iend,kstart-1:nlev+1), INTENT(out) :: &
    &  v_sed_NAT              !< NAT sedimentation velocity (m / s)
  !local variables
  INTEGER :: &
    & jc, jk                  !< loop indices
  REAL(wp), PARAMETER ::  &
    &  d_c = 4.2e-10_wp       !< collision diameter of HNO3 (m) (according to
                              !  Liu et al., J. Phys. Chem, 2007), in EMAC 3.0e-10 is used
  REAL(wp) ::               &
    &  mean_free_path_HNO3, & !< mean free path of HNO3 in air (m)
    &  Cc,                  & !< Cunningham factor (-)
    &  SedFactor              !< sedimentation factor (cf. Carslaw et al., 2002 and comment below)

  v_sed_NAT(:,kstart-1) = 0._wp
  v_sed_NAT(:,nlev+1)   = 0._wp

  DO jk = kstart,nlev
    DO jc = istart, iend
      IF (dens_NAT(jc,jk) > 0.1_wp) THEN
        mean_free_path_HNO3 = 1._wp / (SQRT(2._wp) * pi * pres(jc,jk)   &
                     &                  / argas / temp(jc,jk) * avo * d_c**2)

        Cc = 1._wp + mean_free_path_HNO3 / r_NAT(jc,jk)           &
              &  * (1.257_wp + 0.4_wp * EXP(-1.1_wp * r_NAT(jc,jk)&
              &                              / mean_free_path_HNO3))


        SedFactor = 2._wp * grav * (rho_NAT / 1000._wp) * Cc   &
           &        / 9._wp / dyn_visc(jc,jk)

        ! Attention: This differs from Carslaw et al. (2002) because in the
        ! FixedRad approach the radius of NAT particles is fixed and not
        ! time dependent
        v_sed_NAT(jc,jk) = SedFactor * r_NAT(jc,jk)**2.0_wp
      ELSE
        v_sed_NAT(jc,jk) = 0._wp
      END IF
    END DO
  END DO  

END SUBROUTINE art_psc_KinPar_NAT_sedi_vel

!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!

SUBROUTINE art_psc_simpleUpwind_sedi(HNO3_Nconc_s, v_sed_NAT, dz, p_dtime, jb,  &
                     &               kstart,kend, i_startidx, i_endidx,         &
                     &               NAT_sedi_rel_diff, cell_area)
!<
! SUBROUTINE art_psc_simpleUpwind_sedi
! Routine for performing sedimentation of a number concentration (#/cm3)
! Based on:
! - Buchholz (2005), Dissertation, MPI-C (but accounting for model levels of
!                                         non-equal distance)
! Part of Module: mo_art_psc_sedimentation
! Author: Michael Weimer, KIT
! Initial Release: 2017-09-18
! Modifications:
!>
  IMPLICIT NONE
  REAL(wp), INTENT(inout) :: &
    &  HNO3_Nconc_s(:,:)       !< HNO3 numer concentration within the NAT particles (#/cm3)

  REAL(wp), INTENT(in)    :: &
    &  v_sed_NAT(:,:),       & !< sedimentation velocity of NAT (m / s)
    &  dz(:,:,:),            & !< height of the model layers (m)
    &  p_dtime                 !< model time step (s)

  INTEGER, INTENT(in) :: &
    &  jb,               & !< loop indices
    &  kstart, kend,     &
    &  i_startidx,       &
    &  i_endidx

  REAL(wp), OPTIONAL, INTENT(inout) :: &
    &  NAT_sedi_rel_diff(:,:)   !< relative mass difference in the column before
                                !  to after sedimentation (-)

  REAL(wp), OPTIONAL, INTENT(in) :: &
    &  cell_area(:,:)           !< area of the model cells (m2)

  ! local variables
  INTEGER  :: &
    &  jc, jk            !< loop indices
  REAL(wp)  ::  &
    &  SedStep(i_startidx:i_endidx,kstart:kend),  &  !< height that a particle with the velocity 
                                                     !  moves downwards (m)
    &  change_Nconc(i_startidx:i_endidx,kstart:kend) !< change in number concentration in the 
                                                     !  respective layer (# / cm3)
  LOGICAL ::      &
    &  val_sed(i_startidx:i_endidx,kstart:kend)   !< array deciding if sedimentation happens or not
  REAL(wp) ::     &
    &  mass_NAT    !< mass after sedimentation

  ! calculate NAT mass of a column before sedimentation
  IF (PRESENT(NAT_sedi_rel_diff)) THEN
    DO jc = i_startidx, i_endidx

      NAT_sedi_rel_diff(jc,jb) = cell_area(jc,jb) * mol_weight_NAT             &
               &               * SUM(HNO3_Nconc_s(jc,:) * dz(jc,:,jb))
    END DO
  END IF


  ! initialisation
  DO jk = kstart,kend
    DO jc = i_startidx, i_endidx
      SedStep(jc,jk) = v_sed_NAT(jc,jk) * p_dtime
      val_sed(jc,jk) = (HNO3_Nconc_s(jc,jk) > 0._wp)
      change_Nconc(jc,jk) = 0._wp
    END DO
  END DO

  ! no sedimentation in uppermost and lowermost level
  val_sed(:,kend) = .FALSE.
  val_sed(:,kstart) = .FALSE.

  
  ! calculate the change in number concentration in the respective level (the molecules coming
  ! from above and which are removed from the current level)
  DO jk = MAX(kstart,2),kend
    DO jc = i_startidx, i_endidx
      IF (val_sed(jc,jk)) THEN
        change_Nconc(jc,jk) = - HNO3_Nconc_s(jc,jk)           &
                  &       * SedStep(jc,jk) / dz(jc,jk,jb)
      END IF

      IF (val_sed(jc,jk-1)) THEN
        change_Nconc(jc,jk) = change_Nconc(jc,jk) + HNO3_Nconc_s(jc,jk-1)  &
                   &        * SedStep(jc,jk-1) / dz(jc,jk,jb)
      END IF
      
    END DO
  END DO

  ! add the change in number concentraion to the NAT tracer
  DO jk = kstart,kend
    DO jc = i_startidx, i_endidx
      HNO3_Nconc_s(jc,jk) = HNO3_Nconc_s(jc,jk) + change_Nconc(jc,jk)
    END DO
  END DO

  ! calculate the relative difference of NAT mass
  IF (PRESENT(NAT_sedi_rel_diff)) THEN
    DO jc = i_startidx, i_endidx
      mass_NAT =  cell_area(jc,jb) * mol_weight_NAT  &
           &   * SUM(HNO3_Nconc_s(jc,:) * dz(jc,:,jb))

      IF (mass_NAT > 0._wp) THEN
        NAT_sedi_rel_diff(jc,jb) = (mass_NAT - NAT_sedi_rel_diff(jc,jb)) * 2._wp   &
                   &             / (mass_NAT + NAT_sedi_rel_diff(jc,jb))
      ELSE
        NAT_sedi_rel_diff(jc,jb) = 0._wp
      END IF
    END DO
  END IF

END SUBROUTINE art_psc_simpleUpwind_sedi
END MODULE mo_art_psc_sedimentation
