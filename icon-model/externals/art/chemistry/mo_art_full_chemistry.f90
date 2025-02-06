!
! mo_art_full_chemistry
! This module provides the subroutines for chemical gas phase processes.
! This Module is only envoked if lart_mecca = .TRUE.
! Based on MECCA
! For reference see
! Sander, Rolf, et al. "Technical note: The new comprehensive
!        atmospheric chemistry module MECCA."
!        Atmospheric Chemistry and Physics 5.2 (2005): 445-450.
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

MODULE mo_art_full_chemistry
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: message
  USE mo_impl_constants,                ONLY: SUCCESS

  USE mtime,                            ONLY: datetime
!ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_chem_data,                 ONLY: t_art_chem
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE messy_mecca_kpp_parameters
  USE messy_cmn_photol_mem
  USE mo_mecicon,                       ONLY: mecicon_call
  USE mo_art_diagnostics,               ONLY: art_chem_calc_column
  
  

  IMPLICIT NONE

  PRIVATE

  PUBLIC   ::                        &
      &   art_loss_full_chemistry


CONTAINS

SUBROUTINE art_loss_full_chemistry(current_date,jg,p_dtime, p_tracer_now)
!<
! SUBROUTINE art_loss_full_chemistry
! This subroutine calculates tracer depletion for full gasphase chemistry
! (lart_mecca == .TRUE.)
! Based on: MECCA (?)
! Part of Module: mo_art_full_chemistry
! Author: Jennifer Schroeter, KIT
! Initial Release: 2015-??-??
! Modifications:
!>
  TYPE(datetime),POINTER,INTENT(in) :: &
    &  current_date                      !< mtime object with current date and time
  INTEGER, INTENT(in)               :: &
    &  jg                                !< patch on which computation is performed
  REAL(wp), INTENT(IN)              :: &
    &  p_dtime                           !< time step
  REAL(wp), INTENT(INOUT)           :: &
    &  p_tracer_now(:,:,:,:)             !< tracer mixing ratios at current time level (kg/kg)
  
!Local variables
  INTEGER                           :: &
    &  jc, jk, jb,                     & !< loop indizes
    &  i_startidx, i_endidx,           &
    &  ierror,                         &
    &  iTRO3_passive,                  &
    &  iTRNOy_passive
  TYPE(t_art_chem), POINTER :: &
    &  art_chem
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo
  INTEGER, POINTER :: &
    &  mapping_indices_kpp(:)


! ----------------------------------
! --- Initialisation
! ----------------------------------

  art_chem => p_art_data(jg)%chem
  art_atmo => p_art_data(jg)%atmo
  mapping_indices_kpp => art_chem%mecicon%utils%mapping_indices_kpp

  IF ((current_date%date%month == 5) .AND. (current_date%date%day == 1) &
   &  .AND. (current_date%time%hour > 0 ) .AND. (current_date%time%hour < 6 )) THEN
    CALL message('mo_art_full_chemistry','reinitializing passive O3 and NOy if they are present')

    CALL p_art_data(jg)%dict_tracer%get('TRO3_passive',iTRO3_passive, ierror)

    IF (ierror == SUCCESS) THEN
      p_tracer_now(:,:,:,iTRO3_passive) = p_tracer_now(:,:,:, mapping_indices_kpp(ind_O3))
    END IF

    CALL p_art_data(jg)%dict_tracer%get('TRNOy_passive',iTRNOy_passive, ierror)

    IF (ierror == SUCCESS) THEN
      p_tracer_now(:,:,:,iTRNOy_passive) = p_tracer_now(:,:,:, mapping_indices_kpp(ind_HNO3))  &
                                    &    + p_tracer_now(:,:,:, mapping_indices_kpp(ind_NO2))   &
                                    &    + p_tracer_now(:,:,:, mapping_indices_kpp(ind_NO))
    END IF
  END IF

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    DO jk = 1,art_atmo%nlev
      DO jc = i_startidx,i_endidx

        CALL mecicon_call(jc, jk, jb, jg,  p_dtime, art_atmo%pres(jc,jk,jb),&
                          art_atmo%temp(jc,jk,jb),p_tracer_now)

      ENDDO
    ENDDO
  ENDDO

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    CALL art_chem_calc_column(art_chem%mecicon%utils%o3_column(:,:,jb),            &
                  &            art_atmo%pres(:,:,jb),                              &
                  &           art_atmo%temp(:,:,jb), art_atmo%z_mc(:,:,jb),        &
                  &           p_tracer_now(:,:,jb, mapping_indices_kpp(ind_O3)),   &
                  &           i_startidx, i_endidx, art_atmo%nlev)

  ENDDO

END SUBROUTINE art_loss_full_chemistry
!!
!!-----------------------------------------------------------------------------
!!
END MODULE mo_art_full_chemistry
