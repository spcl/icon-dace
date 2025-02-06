!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
MODULE mo_event_manager

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

    ALLOCATE (model_event_groups(model_event_groups_list_size))
    model_event_groups_list_member = 0

    linitialized = .TRUE.

  END SUBROUTINE initEventManager

  FUNCTION getModelReferenceDate() RESULT(r)
    TYPE(datetime), POINTER :: r

    r => NULL()
    IF (linitialized) THEN
      r => model_reference_date
    END IF

  END FUNCTION getModelReferenceDate

  FUNCTION addEventGroup(group) RESULT(handle)
    INTEGER :: handle
    CHARACTER(len=*), INTENT(in) :: group
    TYPE(event_group_list), ALLOCATABLE :: tmp(:)
    CHARACTER(len=max_groupname_str_len) :: gstring

    IF (.NOT. linitialized) THEN
      PRINT *, 'ERROR: event manager not initialized.'
      STOP
    END IF

    IF (model_event_groups_list_member == model_event_groups_list_size) THEN
      PRINT *, 'INFO: reallocating event group list.'
      ALLOCATE (tmp(model_event_groups_list_size))
      tmp(:) = model_event_groups(:)
      DEALLOCATE (model_event_groups)
      ALLOCATE (model_event_groups(2*model_event_groups_list_size))
      model_event_groups(:model_event_groups_list_size) = tmp(:)
      DEALLOCATE (tmp)
      model_event_groups_list_size = 2*model_event_groups_list_size
      PRINT *, 'INFO: new evcent group list size: ', model_event_groups_list_size
    END IF

    model_event_groups_list_member = model_event_groups_list_member + 1

    model_event_groups(model_event_groups_list_member)%group => newEventGroup(TRIM(group))
    CALL getEventGroupName(model_event_groups(model_event_groups_list_member)%group, gstring)
    PRINT *, 'INFO: Added event group: ', TRIM(gstring)

    handle = model_event_groups_list_member

  END FUNCTION addEventGroup

  FUNCTION getEventGroup(handle) RESULT(eventGroupListMember)
    TYPE(eventGroup), POINTER :: eventGroupListMember
    INTEGER, INTENT(in) :: handle
    IF (handle > model_event_groups_list_member) THEN
      eventGroupListMember => NULL()
    ELSE
      eventGroupListMember => model_event_groups(handle)%group
    END IF
  END FUNCTION getEventGroup

  SUBROUTINE printEventGroup(handle)
    INTEGER, INTENT(in) :: handle
    TYPE(eventGroup), POINTER :: currentEventGroup
    TYPE(event), POINTER :: currentEvent
    CHARACTER(len=max_eventname_str_len) :: estring
    CHARACTER(len=max_groupname_str_len) :: egstring

    currentEventGroup => getEventGroup(handle)
    CALL getEventGroupName(currentEventGroup, egstring)

    currentEvent => getFirstEventFromEventGroup(model_event_groups(handle)%group)

    PRINT *, 'Event list: ', TRIM(egstring)
    DO WHILE (ASSOCIATED(currentEvent))
      CALL eventToString(currentEvent, estring)
      PRINT *, '   event ', TRIM(estring)
      currentEvent => getNextEventFromEventGroup(currentEvent)
    END DO
  END SUBROUTINE printEventGroup

  SUBROUTINE getEventComponents(eventString, timeInterval, startDate, endDate)
    CHARACTER(len=max_repetition_str_len), INTENT(in) :: eventString
    TYPE(timedelta), POINTER :: timeInterval
    TYPE(datetime), POINTER :: startDate
    TYPE(datetime), POINTER :: endDate

    CHARACTER(len=max_repetition_str_len) :: r, s, e, d
    LOGICAL :: lr, ls, le, ld

    CALL splitRepetitionString(eventString, r, s, e, d, lr, ls, le, ld)

    IF (lr) THEN
      IF (getRepetitions(r) /= -1) THEN
        PRINT *, 'WARNING: event setup should not have explicit repeat count.'
      END IF
    END IF

    IF (ls) THEN
      startDate => newDatetime(TRIM(s))
    END IF

    IF (le) THEN
      endDate => newDatetime(TRIM(e))
    END IF

    IF (ld) THEN
      timeInterval => newTimeDelta(TRIM(d))
    ELSE
      PRINT *, 'ERROR: time interval should be given.'
      STOP
    END IF

  END SUBROUTINE getEventComponents

END MODULE mo_event_manager

