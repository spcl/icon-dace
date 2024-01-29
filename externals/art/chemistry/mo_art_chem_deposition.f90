!
! mo_art_chem_deposition
! This module provides the subroutines for chemical tracer deposition to the
! ground
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

MODULE mo_art_chem_deposition
  ! ICON
  USE mo_kind,                        ONLY: wp
  ! ART
  USE mo_art_data,                    ONLY: p_art_data
  USE mo_art_atmo_data,               ONLY: t_art_atmo
  USE mo_art_read_emissions,          ONLY: art_convert_emission_to_mmr
  USE mo_art_wrapper_routines,        ONLY: art_get_indices_c

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: art_CO2_deposition
CONTAINS

SUBROUTINE art_CO2_deposition(jg, CO2_mmr, p_dtime, art_atmo)
!<
! SUBROUTINE art_CO2_deposition
! This subroutine calculates the CO2 deposition to sea surface
! Based on:
! - Jacob, 1999
! Part of Module: mo_art_chem_deposition
! Author: Michael Weimer, KIT
! Initial Release: 2018-11-13
!>
  IMPLICIT NONE
  INTEGER, INTENT(in)  :: &
    &  jg                !< patch id
  REAL(wp), INTENT(inout) :: &
    &  CO2_mmr(:,:,:)    !< mass mixing ratio of CO2 (kg / kg)
  REAL(wp), INTENT(in) :: &
    &  p_dtime           !< model time step
  TYPE(t_art_atmo), INTENT(in) :: &
    &  art_atmo          !< ART atmo fields
  ! local variables
  REAL(wp), PARAMETER ::   &
    &  DEPOS_FACTOR = 0.5_wp    !< reduction factor for calculation of deposition 
                                !  from global mean emission of CO2 (should eventually be adapted)

  ! divide global average by oceanic area fraction on Earth's surface (is there a value in ICON?)
  ! The value 2.77e-9_wp roughly corresponds to 2012 in the emission datasets
  ! (end dates are 2008 and 2010...)
  ! Eventually, it should be combined with a actual ratio CO2(current_year) / CO2(2012) but this 
  ! is not stored in case of AES physics
  REAL(wp), PARAMETER :: &
    &  mean_emiss_kgm2s = 2.77e-9_wp / 0.707_wp          !< emission in kg / m2 / s
  INTEGER :: &
    & jc,jb, i_startidx, i_endidx

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    DO jc = i_startidx,i_endidx
      ! convert global mean emission to mass mixing ratio and subtract it 
      ! over the sea with a factor (50 % in Jacob, 1999)
      IF (.NOT. art_atmo%llsm(jc,jb)) THEN
        CALL art_convert_emission_to_mmr(p_art_data(jg)%chem%CO2_mmr_depos,   &
                       &                 mean_emiss_kgm2s,                    &
                       &                 art_atmo%temp(jc,:,jb),              &
                       &                 art_atmo%pres(jc,:,jb),              &
                       &                 art_atmo%dz(jc,:,jb),                &
                       &                 p_dtime,1,art_atmo%nlev)

        CO2_mmr(jc,art_atmo%nlev,jb) =  CO2_mmr(jc,art_atmo%nlev,jb)                         &
                       &                - DEPOS_FACTOR                                       &
                       &                  *  p_art_data(jg)%chem%CO2_mmr_depos(art_atmo%nlev)
      END IF
    END DO
  END DO
END SUBROUTINE art_CO2_deposition


END MODULE mo_art_chem_deposition
