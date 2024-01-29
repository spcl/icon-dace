! Module provides a grouping structure for mtime events.
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

MODULE mo_event_manager

  USE mo_exception, ONLY: finish, message, message_text
  USE mtime

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: initEventManager 
  PUBLIC :: getModelReferenceDate
  PUBLIC :: addEventGroup
  PUBLIC :: getEventGroup 
  PUBLIC :: printEventGroup 
  PUBLIC :: getEventComponents

  TYPE event_group_list
    TYPE(eventgroup), POINTER :: group
  END TYPE event_group_list

  TYPE(event_group_list), ALLOCATABLE :: model_event_groups(:)

  INTEGER :: model_event_groups_list_member
  INTEGER :: model_event_groups_list_size = 16

  TYPE(datetime), POINTER :: model_reference_date => NULL()

  LOGICAL :: linitialized = .FALSE.

CONTAINS

  SUBROUTINE initEventManager(referenceDate)
    TYPE(datetime), POINTER, INTENT(in) :: referenceDate

    model_reference_date => newDatetime(referenceDate)

    ALLOCATE(model_event_groups(model_event_groups_list_size)) 
    model_event_groups_list_member = 0

    linitialized = .TRUE.

  END SUBROUTINE initEventManager

  FUNCTION getModelReferenceDate() RESULT(r)
    TYPE(datetime), POINTER :: r

    r => NULL()
    IF (linitialized) THEN
      r => model_reference_date
    ENDIF

  END FUNCTION getModelReferenceDate

  FUNCTION addEventGroup(group) RESULT(handle)
    INTEGER :: handle
    CHARACTER(len=*), INTENT(in) :: group
    TYPE(event_group_list), ALLOCATABLE :: tmp(:)
    CHARACTER(len=max_groupname_str_len) :: gstring
    INTEGER :: new_size

    IF (.NOT. linitialized) THEN
      CALL finish('', 'event manager not initialized.') 
    ENDIF

    IF (model_event_groups_list_member == model_event_groups_list_size) THEN
      CALL message('', 'reallocating event group list.') 
      new_size = 2*model_event_groups_list_size
      ALLOCATE(tmp(new_size))
      tmp(1:model_event_groups_list_size) = model_event_groups(:)
      CALL MOVE_ALLOC(tmp, model_event_groups)
      model_event_groups_list_size = new_size
      WRITE(message_text,'(a,i0)') 'new event group list size: ', &
           model_event_groups_list_size
      CALL message('', message_text)
    ENDIF

    model_event_groups_list_member = model_event_groups_list_member + 1

    model_event_groups(model_event_groups_list_member)%group => newEventGroup(TRIM(group))
    CALL getEventGroupName(model_event_groups(model_event_groups_list_member)%group, gstring)
    message_text = 'added event group: '//gstring
    CALL message('', message_text)

    handle = model_event_groups_list_member

  END FUNCTION addEventGroup

  FUNCTION getEventGroup(handle) RESULT(eventGroupListMember)
    TYPE(eventGroup), POINTER :: eventGroupListMember
    INTEGER, INTENT(in) :: handle
    IF (handle > model_event_groups_list_member) THEN
       eventGroupListMember => NULL()
     ELSE
       eventGroupListMember =>  model_event_groups(handle)%group
     ENDIF
  END FUNCTION getEventGroup

  SUBROUTINE printEventGroup(handle)
    INTEGER, INTENT(in) :: handle
    TYPE(eventGroup), POINTER :: currentEventGroup
    TYPE(event), POINTER :: currentEvent
    CHARACTER(len=max_eventname_str_len) :: estring
    CHARACTER(len=max_groupname_str_len) :: egstring
    INTEGER :: icount

    currentEventGroup => getEventGroup(handle)
    CALL getEventGroupName(currentEventGroup, egstring)

    currentEvent => getFirstEventFromEventGroup(model_event_groups(handle)%group)

    icount = 1
    WRITE(message_text,'(a,a)') 'Event list: ', TRIM(egstring)
    CALL message('', message_text)
    DO WHILE (ASSOCIATED(currentEvent))
      CALL eventToString(currentEvent, estring)
      WRITE(message_text,'(5x,i3,10x,a)') icount, TRIM(estring)
      CALL message('', message_text)
      currentEvent => getNextEventFromEventGroup(currentEvent)
      icount = icount+1
    ENDDO
  END SUBROUTINE printEventGroup

  SUBROUTINE getEventComponents(eventString, referenceDate, timeInterval, startDate, endDate)
    CHARACTER(len=max_repetition_str_len), INTENT(in) :: eventString 
    TYPE(datetime),  POINTER :: referenceDate
    TYPE(timedelta), POINTER :: timeInterval
    TYPE(datetime),  POINTER :: startDate
    TYPE(datetime),  POINTER :: endDate
    
    CHARACTER(len=max_repetition_str_len) :: r, s, e, d    
    LOGICAL :: lr, ls, le, ld    
    
    CALL splitRepetitionString(eventString, r, s, e, d, lr, ls, le, ld)
    
    IF (lr) THEN
      IF (getRepetitions(r) /= -1) THEN
        CALL message('', 'event setup should not have explicit repeat count.')
      ENDIF
    ENDIF
    
    IF (ls) THEN
      startDate => newDatetime(TRIM(s))
    ENDIF
    
    IF (le) THEN
      endDate => newDatetime(TRIM(e))
    ENDIF
    
    IF (ld) THEN
      timeInterval => newTimeDelta(TRIM(d))
    ELSE
      CALL finish('', 'time interval should be given.')
    ENDIF

  END SUBROUTINE getEventComponents

END MODULE mo_event_manager
