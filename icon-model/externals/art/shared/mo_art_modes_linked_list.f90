!
! mo_art_modes_linked_list
! This module provides a linked list for modes
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

MODULE mo_art_modes_linked_list
! ICON
  USE mo_exception,                     ONLY: message,finish
! ART
  USE mo_art_impl_constants,            ONLY: IART_MODE_1MOM,                  &
                                          &   IART_MODE_RADIO, IART_MODE_VOLC, &
                                          &   IART_MODE_POLL, IART_MODE_2MOM,  &
                                          &   IART_EMISS2TRACER
  USE mo_art_modes,                     ONLY: t_mode_fields,                   &
                                          &   t_fields_2mom,t_fields_1mom,     &
                                          &   t_fields_pollen,t_fields_radio,  &
                                          &   t_fields_volc

  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_modes_linked_list'
  
  PUBLIC :: p_mode_state
  PUBLIC :: mode_lists
  PUBLIC :: max_mode_lists
  PUBLIC :: n_mode_lists
  PUBLIC :: new_mode_list
  PUBLIC :: append_mode
  PUBLIC :: delete_mode_list
  PUBLIC :: t_mode_state
  PUBLIC :: t_mode
  PUBLIC :: t_mode_list ! maybe temporary?
  
  !-----------------------------
  !-- Derived Types for modes
  !-----------------------------

  TYPE t_mode_list_intrinsic
    CHARACTER(LEN=12)    :: name         !< name of mode list
    INTEGER              :: nmodes       !< number of modes in list
    TYPE(t_mode),POINTER :: first_mode   !< reference to first
  END TYPE t_mode_list_intrinsic

  TYPE t_mode_list
    TYPE(t_mode_list_intrinsic), POINTER :: p
  END TYPE t_mode_list

  TYPE t_mode_state
    TYPE(t_mode_list),POINTER :: p_mode_list
  END TYPE t_mode_state

  TYPE t_mode
    CLASS(t_mode_fields),ALLOCATABLE   :: fields
    TYPE(t_mode), POINTER :: next_mode
  END TYPE t_mode


  !-----------------------------
  !-- Declarations
  !-----------------------------

  TYPE(t_mode_state),ALLOCATABLE :: & 
    &    p_mode_state(:)                 !< mode list dim[n_dom]
  
  TYPE(t_mode_list),ALLOCATABLE, TARGET, SAVE :: &
    &    mode_lists(:)                   !< memory buffer array

  INTEGER ::            &
    &   max_mode_lists, &                !< maximum number of mode lists
    &   n_mode_lists                     !< actual number of mode lists

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE new_mode_list(this_list)

  TYPE(t_mode_list), INTENT(inout) :: this_list
  ALLOCATE(this_list%p)
  
  !Set defaults:
  this_list%p%name       = ''
  this_list%p%nmodes     = 0
  this_list%p%first_mode => NULL()
END SUBROUTINE new_mode_list

SUBROUTINE delete_mode_list(this_list)

  TYPE(t_mode_list), INTENT(inout) :: this_list
  
  CALL message('','ART: Deconstruction of mode list: '//this_list%p%name)
  
  CALL delete_modes(this_list,this_list%p%first_mode)
  
  this_list%p%first_mode => NULL()

  IF (this_list%p%nmodes /= 0) THEN
    CALL finish ('delete_mode_list', 'List delete didnt work proper (nmodes)')
  ENDIF
  
END SUBROUTINE delete_mode_list
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE create_mode(this_list,current_mode,IART_MOMENT) !PRIVATE ROUTINE

  TYPE(t_mode_list),INTENT(inout) :: this_list
  TYPE(t_mode), POINTER           :: current_mode
  INTEGER,INTENT(in)              :: IART_MOMENT
  INTEGER :: istat

  ALLOCATE (current_mode, STAT=istat)
    IF (istat /= 0) THEN
      CALL finish('mo_art_modes_linked_list:create_mode',          &
           &      'allocation of mode failed')
    ENDIF
  
  this_list%p%nmodes = this_list%p%nmodes+1

  ! allocate fields depending on modal type
  ! and fill with initial values
  SELECT CASE (IART_MOMENT)
    CASE (IART_MODE_2MOM)
      ALLOCATE(t_fields_2mom   :: current_mode%fields)
    CASE (IART_MODE_1MOM)
      ALLOCATE(t_fields_1mom   :: current_mode%fields)
    CASE (IART_MODE_RADIO)
      ALLOCATE(t_fields_radio  :: current_mode%fields)
    CASE (IART_MODE_POLL)
      ALLOCATE(t_fields_pollen :: current_mode%fields)
    CASE (IART_MODE_VOLC)
      ALLOCATE(t_fields_volc   :: current_mode%fields)
    CASE (IART_EMISS2TRACER)
      !ALLOCATE(t_art_emiss2tracer   :: current_mode%fields)  ! in current design not possible
                                                              ! here, but will be done past this
    CASE DEFAULT
      CALL finish('mo_art_modes_linked_list:create_mode',          &
           &      'unknown modal type.')
  END SELECT

  current_mode%next_mode => NULL()
  
END SUBROUTINE create_mode
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE append_mode(this_list,new_mode,IART_MOMENT)
  
  TYPE(t_mode_list),INTENT(inout) :: this_list
  INTEGER,INTENT(in)           :: &
    &    IART_MOMENT
  TYPE(t_mode), POINTER        :: & ! INTENT(inout)
    &    new_mode
  ! Local variables
  TYPE(t_mode), POINTER        :: &
    &    current_mode

  ! insert as first mode if list is empty
  IF (.NOT. ASSOCIATED (this_list%p%first_mode)) THEN
    CALL create_mode (this_list, this_list%p%first_mode,IART_MOMENT)
    new_mode => this_list%p%first_mode

    RETURN
  ENDIF
    
  ! loop over list modes to find position
  current_mode => this_list%p%first_mode
  DO WHILE (ASSOCIATED(current_mode%next_mode)) 
    current_mode => current_mode%next_mode
  ENDDO
    
  ! insert mode between current_mode and current_mode%next_mode
  CALL create_mode (this_list, new_mode,IART_MOMENT)
  new_mode%next_mode     => current_mode%next_mode
  current_mode%next_mode => new_mode

END SUBROUTINE append_mode
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE delete_modes(this_list, this_list_mode)
  TYPE(t_mode_list),INTENT(inout) :: this_list
  TYPE(t_mode), POINTER           :: this_list_mode
    
  TYPE(t_mode), POINTER           :: this, next

  ! set next pointer to the first mode (which is intent(in) actually)
  next => this_list_mode
  ! disassociate the intent(in) pointer
  this_list_mode => NULL()
  
  DO
    IF (.NOT. ASSOCIATED(next)) EXIT
    this => next
    next => this%next_mode
    
    SELECT TYPE(fields=>this%fields)
      TYPE is (t_fields_2mom)
        CALL message('','ART: Deallocation of mode: '//fields%name)
        IF (ASSOCIATED(fields%third_moment))     DEALLOCATE(fields%third_moment)
        IF (ASSOCIATED(fields%mass))             DEALLOCATE(fields%mass)
        IF (ASSOCIATED(fields%density))          DEALLOCATE(fields%density)
!        IF (ASSOCIATED(fields%diameter)) DEALLOCATE(fields%diameter)! Is needed for diagnostics,
                                                                     ! ordering makes problems here
        IF (ASSOCIATED(fields%flx_contra_vsed0)) DEALLOCATE(fields%flx_contra_vsed0)
        IF (ASSOCIATED(fields%flx_contra_vsed3)) DEALLOCATE(fields%flx_contra_vsed3)
        IF (ASSOCIATED(fields%dmdt_condso4))     DEALLOCATE(fields%dmdt_condso4)
        IF (ASSOCIATED(fields%knudsen_nr))       DEALLOCATE(fields%knudsen_nr)
      CLASS is (t_fields_1mom)
        CALL message('','ART: Deallocation of mode: '//fields%name)
        SELECT TYPE(fields)
          CLASS is (t_fields_pollen)
            !$ACC EXIT DATA DELETE(fields%flx_contra_vsed)
          CLASS DEFAULT
            CALL message('','ART: Non-pollen fields are not allocated on GPU at the moment')
        END SELECT
        DEALLOCATE(fields%flx_contra_vsed)
      CLASS DEFAULT
        CALL message('','ART: Unknown mode field type')
    END SELECT
    ! One mode was deleted, so adjust the number of modes in container
    this_list%p%nmodes = this_list%p%nmodes-1
    
    DEALLOCATE (this)
  ENDDO
END SUBROUTINE delete_modes
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_modes_linked_list
