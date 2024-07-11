!
! mo_art_kinetic_constants
! This module provides the routines to calculate bi- and termolecular
! kinetic constants
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

MODULE mo_art_kinetic_constants
  USE mo_kind,                    ONLY: wp
  
  USE mo_art_data,                ONLY: p_art_data
  USE mo_art_atmo_data,           ONLY: t_art_atmo
  USE mo_art_wrapper_routines,    ONLY: art_get_indices_c
  USE mo_physical_constants,      ONLY: p0sl_bg


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_determine_bimolecular_kinetic_constant
  PUBLIC :: art_determine_termolecular_kinetic_constant
  PUBLIC :: art_CH4_CO_kinetic_constants
  PUBLIC :: art_get_CO_des_1d
  PUBLIC :: art_get_CH4_des_1d

CONTAINS

SUBROUTINE art_determine_bimolecular_kinetic_constant(k,temp,A,ER,add)
!<
! SUBROUTINE art_determine_bimolecular_kinetic_constant
! This subroutine calculates bimolecular kinetic constants
! Based on: Sander (2011), JPL Publications 10-6
! Part of Module: mo_art_kinetic_constants
! Author: Michael Weimer, KIT
! Initial Release: 2016-04-05
! Modifications:
!>
  IMPLICIT NONE
  REAL(wp), INTENT(in) :: &
    & temp    ! current temperature
  REAL(wp), INTENT(out) :: &
    & k       ! calculated kinetic constant
  REAL(wp), INTENT(in) :: &
    & A, ER   ! constants for bimolecular reaction as given in JPL database
              ! ("add" is nonzero only for reaction of acetone)
  REAL(wp), INTENT(in), OPTIONAL ::  &
    & add     ! "add" is nonzero only for reaction of acetone

  k = A * EXP(- ER / temp)

  IF (PRESENT(add)) THEN
    ! at least for reaction of OH + CH3COCH3 --> [products],
    ! a constant has to be added to the "normal" formular 
    ! of the kinetic constant
    k = k + add
  END IF

END SUBROUTINE art_determine_bimolecular_kinetic_constant
!!
!!-----------------------------------------------------------------------------
!!
SUBROUTINE  art_determine_termolecular_kinetic_constant(k,temp,k300_0,n,  &
                        &                               k300_inf,m,Nconc_M, ca)
!<
! SUBROUTINE art_determine_termolecular_kinetic_constant
! This subroutine calculates "pseudo"-termolecular kinetic constants
! Based on: Sander (2011), JPL Publications 10-6
! Part of Module: mo_art_kinetic_constants
! Author: Michael Weimer, KIT
! Initial Release: 2016-04-05
! Modifications:
!>
  IMPLICIT NONE
  REAL(wp), INTENT(in) :: &
    & temp    ! current temperature
  REAL(wp), INTENT(out) :: &
    & k       ! calculated kinetic constant
  REAL(wp), INTENT(in) :: &
    & k300_0, n, k300_inf, m  ! constants for termolecular reaction
                              ! (see JPL description for details)
  REAL(wp), INTENT(in) :: &
    & Nconc_M          ! number density of collision partner
                       ! in units of molecules / cm^3
  LOGICAL, INTENT(in) :: &
    & ca               ! formular for chemical activation? (see JPL description)

  ! local variables
  REAL(wp) :: &
    & k_0, k_inf, frac

   ! save some values which are needed more than one time
   k_0   = k300_0   * (temp / 300._wp)**(-n)
   k_inf = k300_inf * (temp / 300._wp)**(-m)

   IF (.NOT. ca) THEN
     k_0 = k_0 * Nconc_M
   ELSE
     k_inf = k_inf / Nconc_M
   END IF
   frac      = k_0 / k_inf

   ! calculate kinetic constant
   k = k_0 / (1._wp + frac) * 0.6_wp**((1._wp + (LOG10(frac))**(2._wp))**(-1._wp))

END SUBROUTINE art_determine_termolecular_kinetic_constant
!!
!!-----------------------------------------------------------------------------
!!

!!
!!-----------------------------------------------------------------------------
!!
SUBROUTINE art_CH4_CO_kinetic_constants(jg,k_CH4_OH,k_CO_OH,vmr2Nconc,temp)
!<
! SUBROUTINE 
! This subroutine calculates "pseudo"-termolecular kinetic constants
! Based on: Sander (2011), JPL Publications 10-6
! Part of Module: mo_art_kinetic_constants
! Author: Michael Weimer, KIT
! Initial Release: 2016-04-05
! Modifications:
!>
  IMPLICIT NONE
  INTEGER, INTENT(in)       :: &
    &  jg
  REAL(wp), INTENT(inout) :: &
    &  k_CH4_OH(:,:,:),      &
    &  k_CO_OH(:,:,:)
  REAL(wp), INTENT(in) :: &
    &  vmr2Nconc(:,:,:),  &
    &  temp(:,:,:)

  ! local variables
  INTEGER :: &
    &  jc,jk,jb, i_startidx, i_endidx
  REAL(wp) :: &
    &  N2_O2_Nconc, k_CO_OH1, k_CO_OH2
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo


  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    DO jk = 1,art_atmo%nlev
      DO jc = i_startidx,i_endidx
        N2_O2_Nconc = 0.9903 * vmr2Nconc(jc,jk,jb)
        CALL art_determine_termolecular_kinetic_constant(k_CO_OH1,        &
                                      &                  temp(jc,jk,jb),  &
                                      &                  5.9e-33_wp,      &
                                      &                  1.4_wp,          &
                                      &                  1.1e-12_wp,      &
                                      &                  -1.3_wp,         &
                                      &                  N2_O2_Nconc,     &
                                      &                  .FALSE.)
        CALL art_determine_termolecular_kinetic_constant(k_CO_OH2,        &
                                      &                  temp(jc,jk,jb),  &
                                      &                  1.5e-13_wp,      &
                                      &                  -0.6_wp,         &
                                      &                  2.1e9_wp,        &
                                      &                  -6.1_wp,         &
                                      &                  N2_O2_Nconc,     &
                                      &                  .TRUE.)

        k_CO_OH(jc,jk,jb)  = k_CO_OH1 + k_CO_OH2

        CALL art_determine_bimolecular_kinetic_constant  &
          & (k_CH4_OH(jc,jk,jb),temp(jc,jk,jb),2.45e-12_wp,1775._wp)
      END DO
    END DO
  END DO
  

END SUBROUTINE art_CH4_CO_kinetic_constants


SUBROUTINE art_get_CO_des_1d(des,pres)
!<
! SUBROUTINE 
! Based on:
! - KASIMA
! - Brasseur and Solomon?
! Part of Module: mo_art_kinetic_constants
! Author: Christian Stassen and  Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
!>
  IMPLICIT NONE

  REAL(wp), INTENT(out) :: des
  REAL(wp), INTENT(in) :: pres
  ! local variables
  REAL(wp) :: pres_alt

  pres_alt = -7000._wp * LOG(pres/p0sl_bg)
  des  = (                                                                              &
          &  2.0e-8_wp/(1._wp+EXP( (pres_alt-85.e3_wp)    /(0.4_wp*7000._wp)  ) )       &
          &        + 1.e-7_wp*EXP(-(pres_alt-70.e3_wp)**2._wp/(    7000._wp)**2._wp)    &
          &        + 2.e-6_wp*EXP(-(pres_alt-45.e3_wp)**2._wp/( 2._wp*7000._wp)**2._wp) &
          & )

END SUBROUTINE art_get_CO_des_1d


SUBROUTINE art_get_CH4_des_1d(des,pres)
!<
! SUBROUTINE 
! Based on:
! - KASIMA
! - Brasseur and Solomon?
! Part of Module: mo_art_kinetic_constants
! Author: Christian Stassen and  Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
!>
  IMPLICIT NONE

  REAL(wp), INTENT(out) :: des
  REAL(wp), INTENT(in) :: pres
  ! local variables
  REAL(wp) :: pres_alt

  pres_alt = -7000._wp * LOG(pres/p0sl_bg)
  des = (                                                                                        &
      &     1.20e-7_wp *          EXP(-(pres_alt - 44.e3_wp )**2._wp /(1.5_wp*7000._wp)**2._wp ) &
      &   + 2.60e-9_wp /( 1._wp + EXP(-(pres_alt - 20.e3_wp )        /7000._wp) )                &
      &   + 3.85e-6_wp /( 1._wp + EXP(-(pres_alt - 80.e3_wp )        /(0.5_wp*7000._wp)) )       &
      & )
END SUBROUTINE art_get_CH4_des_1d


END MODULE mo_art_kinetic_constants
