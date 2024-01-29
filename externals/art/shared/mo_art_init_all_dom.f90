!
! mo_art_init_all_dom
! This module provides initialization routines for ICON-ART for all domains (n_dom)
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

MODULE mo_art_init_all_dom
! ICON
  USE mo_exception,                     ONLY: message,finish
! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_collect_atmo_state,        ONLY: art_collect_indices
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, mode_lists,    &
    &                                         n_mode_lists, max_mode_lists

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_init_all_dom

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_init_all_dom(n_dom)
!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------
!--- SUBROUTINE art_init_all_dom                 Author: Daniel Rieger
!--- First Release: 2013/05/16                   daniel.rieger@kit.edu
!--- This routine creates all structures depending on the number
!--- of domains (n_dom). CALL by: mo_art_init_interface
!--- Changelog:
!---            2013/05/16: Initial Release (Daniel Rieger)
!---------------------------------------------------------------------
  
  INTEGER,INTENT(in)  :: &
    &    n_dom             !< number of domains
    
  INTEGER             :: &
    &    istat,          & !< debug integer
    &    jg                !< patch ID
    
  ! ----------------------------------
  ! --- 1.0 Initialization
  ! ----------------------------------
  
  istat          = 0
  n_mode_lists   = 0       !< initialization of the counter of mode lists
  max_mode_lists = n_dom
  
  ! ----------------------------------
  ! --- 2.0 Create and allocate ART data structures
  ! ----------------------------------
  
  CALL message('','ART: Construction of data structure')
  
  ALLOCATE (p_art_data(n_dom),STAT=istat)
    IF (istat /= 0) THEN
      CALL finish('mo_art_init_all_dom:art_init_all_dom',          &
           &      'allocation of array p_art_data failed')
    ENDIF

  !$ACC ENTER DATA COPYIN(p_art_data)

  DO jg = 1,n_dom
    CALL art_collect_indices(jg)
  END DO

  ! ----------------------------------
  ! --- 3.0 Create mode linked list
  ! ----------------------------------
  
  CALL message('','ART: Construction of mode list')
  
  ALLOCATE(p_mode_state(n_dom),STAT=istat)
    IF (istat /= 0) THEN
      CALL finish('mo_art_init_all_dom:art_init_all_dom',          &
           &      'allocation of array p_mode_state failed')
    ELSE
      DO jg =1, n_dom
        p_mode_state(jg)%p_mode_list => NULL()
      ENDDO
    ENDIF
  
  ! ----------------------------------
  ! --- 3.0 Allocate memory buffer array for the mode list
  ! ----------------------------------
  
  CALL message('','ART: Allocating memory buffer array for mode lists')
  
  ALLOCATE(mode_lists(max_mode_lists),STAT=istat)
    IF (istat /= 0) THEN
      CALL finish('mo_art_init_all_dom:art_init_all_dom',          &
           &      'allocation of array mode_lists failed')
    ENDIF
  
END SUBROUTINE art_init_all_dom
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_init_all_dom
