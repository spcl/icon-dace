!
! mo_art_clean_up
! This module provides all clean up for ART structures
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

MODULE mo_art_clean_up
! ICON Routines
  USE mo_exception,                     ONLY: message
! ART Routines
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, delete_mode_list
  USE mo_art_collect_atmo_state,        ONLY: art_deallocate_atmo_state
  USE mo_art_deallocate_chemistry,      ONLY: art_deallocate_chemistry
  USE mo_art_prescribed_state,          ONLY: art_delete_prescr_list
  USE mo_art_diag_types,                ONLY: art_deallocate_diag
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_clean_up

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clean_up(n_dom)
!<
! SUBROUTINE art_clean_up
! This subroutine does the clean up for ART state
! Based on: -
! Part of Module: mo_art_clean_up
! Author: Daniel Rieger, KIT
! Initial Release: 2014-01-27
! Modifications:
! YYYY-MM-DD: <name>, <instituttion>
! - ...
!>
  INTEGER,INTENT(in) :: &
    &  n_dom              !< number of model domains
    
  INTEGER            :: &
    &  jg                 !< counter for the number of domains
    
  CALL message('','ART: Performing clean up')
  
  ! ----------------------------------
  ! --- 1.0 ART data structure p_art_data(n_dom)
  ! ----------------------------------
    
  IF(ALLOCATED(p_art_data)) THEN

      !$ACC WAIT
      DO jg=1,n_dom
        !$ACC EXIT DATA DELETE(p_art_data(jg)%air_prop%art_dyn_visc, p_art_data(jg)%air_prop%art_free_path) &
        !$ACC   DELETE(p_art_data(jg)%turb_fields%sv, p_art_data(jg)%turb_fields%vdep)
        DEALLOCATE(p_art_data(jg)%air_prop%art_free_path)
        DEALLOCATE(p_art_data(jg)%air_prop%art_dyn_visc)
        DEALLOCATE(p_art_data(jg)%turb_fields%sv)
        DEALLOCATE(p_art_data(jg)%turb_fields%vdep)
        CALL p_art_data(jg)%emiss%free()
        CALL art_delete_prescr_list(p_art_data(jg)%prescr_list)

        IF (art_config(jg)%lart_chem) THEN
          CALL art_deallocate_chemistry(jg, p_art_data(jg)%chem)
        END IF

        CALL art_deallocate_atmo_state(jg)

        CALL art_deallocate_diag(p_art_data(jg)%diag)

      ENDDO
        
      CALL message('','ART: Clean up: p_art_data')

      !$ACC WAIT
      !$ACC EXIT DATA DELETE(p_art_data)
      DEALLOCATE(p_art_data)

  END IF
    
  ! ----------------------------------
  ! --- 2.0 ART modes structure
  ! ----------------------------------
    
  DO jg = 1, n_dom
    ! Without aerosols, nothing will be associated
    IF(ASSOCIATED(p_mode_state(jg)%p_mode_list)) THEN
      CALL delete_mode_list(p_mode_state(jg)%p_mode_list)
    ENDIF
  ENDDO
    
  CALL message('','ART: Clean up: p_mode_state')
  DEALLOCATE(p_mode_state) ! This object is always allocated
  
  ! ----------------------------------
  ! --- 3.0 End of clean up
  ! ----------------------------------
    
  CALL message('','ART: Clean up for ART done')
    
END SUBROUTINE art_clean_up
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_clean_up
