!
! mo_art_nucl_aero
! Calculation of sulfate nucleation
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

MODULE mo_art_nucl_aero
! ICON
  USE mo_kind,                          ONLY: wp
! ART

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_nucl_aero'

  PUBLIC :: art_nucl_kw

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_nucl_kw(relhum,temp,rho,ch2so4,istart,iend,kstart,kend,nucmass,totmasscond)
!<
! SUBROUTINE art_nucl_kw
! Calculation of sulfate nucleation based on Kerminen and Wexler method
! Based on: Riemer, 2002: Numerische Simulationen zur Wirkung des Aerosols auf die 
!                         troposphaerische Chemie und die Sichtweite
! Part of Module: mo_art_nucl_aero
! Author: Daniel Rieger, KIT
! Initial Release: 2017-05-09
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp), INTENT(in)           :: &
    &  relhum(:,:),                 & !< Relative humidity (%)
    &  temp(:,:),                   & !< Temperature (K)
    &  rho(:,:),                    & !< Air density (kg m-3)
    &  totmasscond(:,:)               !<Total condensated mass (mug kg-1)
!   &  ch2so4(:,:)                    !< Current H2SO4 concentration (mug kg-1)
  INTEGER, INTENT(in)            :: &
    &  istart,iend,                 & !< Start and end index of nproma loop
    &  kstart,kend                    !< Start and end index of vertical loop
  REAL(wp), INTENT(out)          :: &
    &  nucmass(:,:)                   !< Mass of SO4 that nucleates (mug kg-1)
  REAL(wp), INTENT(inout)          :: &
    &  ch2so4(:,:)                    !< Current H2SO4 concentration (mug kg-1)
!Local variables
  REAL(wp)                       :: &
    &  ckrit,                      & !< critical concentration (mug kg-1)
    &  kg2mug                        !< Factor kg/kg to amug/kg
  INTEGER                        :: &
    &  jc, jk                         !< Loop indices

  nucmass(:,:) = 0._wp
  kg2mug     =  1.E+09_wp

  DO jk=kstart,kend
!NEC$ ivdep
    DO jc=istart,iend
! First update of gas phase H2SO4 concentration
      ch2so4(jc,jk) = ch2so4(jc,jk) * kg2mug - totmasscond(jc,jk)
      ckrit = 0.16_wp * exp(0.1_wp * temp(jc,jk) - 3.5_wp * relhum(jc,jk) / 100._wp - 27.7_wp) / rho(jc,jk)
      IF( ch2so4(jc,jk) > ckrit) THEN
        nucmass(jc,jk) = ch2so4(jc,jk) - ckrit
        ch2so4(jc,jk) = ch2so4(jc,jk) - nucmass(jc,jk) 
      ENDIF
      ch2so4(jc,jk) = ch2so4(jc,jk) / kg2mug 
    ENDDO !jc
  ENDDO !jk

END SUBROUTINE art_nucl_kw
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_nucl_aero
