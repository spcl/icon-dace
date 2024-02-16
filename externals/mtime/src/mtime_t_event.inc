!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
  ! ================================================================================
  ! event section:
  ! ================================================================================

  ! generic assignment for constructors
  !
  FUNCTION t_event_assign_raw(name, referenceDate, firstdate, lastDate, interval, offset)
    TYPE(t_event) :: t_event_assign_raw
    CHARACTER(len=*), INTENT(in)           :: name
    CHARACTER(len=*), INTENT(in)           :: referenceDate
    CHARACTER(len=*), INTENT(in)           :: firstDate
    CHARACTER(len=*), INTENT(in)           :: lastDate
    CHARACTER(len=*), INTENT(in)           :: interval
    CHARACTER(len=*), INTENT(in), OPTIONAL :: offset

    CHARACTER(len=4)                       :: zeroOffset = "PT0S"

    t_event_assign_raw%eventName              = name
    t_event_assign_raw%eventReferenceDateTime = t_datetime(referenceDate)
    t_event_assign_raw%eventFirstDateTime     = t_datetime(firstDate)
    t_event_assign_raw%eventLastDateTime      = t_datetime(lastDate)
    t_event_assign_raw%eventInterval          = t_timedelta(interval)
    IF (PRESENT(offset)) THEN
      t_event_assign_raw%eventOffset = t_timedelta(offset)
    ELSE
      t_event_assign_raw%eventOffset = t_timedelta(zeroOffset)
    ENDIF

  END FUNCTION t_event_assign_raw

  FUNCTION t_event_assign_types(name, referenceDate, firstdate, lastDate, interval, offset)
    TYPE(t_event) :: t_event_assign_types
    CHARACTER(len=*), INTENT(in)            :: name
    TYPE(t_datetime), INTENT(in)            :: referenceDate
    TYPE(t_datetime), INTENT(in)            :: firstDate
    TYPE(t_datetime), INTENT(in)            :: lastDate
    TYPE(t_timedelta), INTENT(in)           :: interval
    TYPE(t_timedelta), INTENT(in), OPTIONAL :: offset

    CHARACTER(len=4)                       :: zeroOffset = "PT0S"

    t_event_assign_types%eventName              = name
    t_event_assign_types%eventReferenceDateTime = referenceDate
    t_event_assign_types%eventFirstDateTime     = firstDate
    t_event_assign_types%eventLastDateTime      = lastDate
    t_event_assign_types%eventInterval          = interval
    IF (PRESENT(offset)) THEN
      t_event_assign_types%eventOffset = offset
    ELSE
      t_event_assign_types%eventOffset = t_timedelta(zeroOffset)
    ENDIF

  END FUNCTION t_event_assign_types

  ! Iterate to next event in event group.
  !
  ! @return NULL() if no next event available.
  FUNCTION t_event_next_event(this)
    TYPE(t_event), POINTER :: t_event_next_event
    CLASS(t_event) :: this
    t_event_next_event => NULL()
    IF (ASSOCIATED(this%nextEventInGroup)) THEN
      t_event_next_event => this%nextEventInGroup
    ENDIF
  END FUNCTION t_event_next_event

  FUNCTION t_event_getId(this) RESULT(res)
    INTEGER(c_int64_t) :: res
    CLASS (t_event) :: this
    res = this%eventId
  END FUNCTION t_event_getId

  FUNCTION t_event_getName(this) RESULT(res)
    CHARACTER(len=max_event_str_len) :: res
    CLASS (t_event) :: this
    res = this%eventName
  END FUNCTION t_event_getName

  FUNCTION t_event_getFirstDatetime (this) RESULT(res)
    TYPE(t_datetime)        :: res
    CLASS(t_event)          :: this
    res = this%eventFirstDateTime
  END FUNCTION t_event_getFirstDatetime

  FUNCTION t_event_getInterval(this) RESULT(res)
    TYPE(t_timedelta)        :: res
    CLASS(t_event)           :: this
    res = this%eventInterval
  END FUNCTION t_event_getInterval

  FUNCTION t_event_getLastDatetime(this) RESULT(res)
    TYPE(t_datetime)        :: res
    CLASS(t_event)          :: this
    res = this%eventLastDateTime
  END FUNCTION t_event_getLastDatetime

  FUNCTION t_event_getNextOccurrenceDatetime(this, my_currentdatetime) RESULT(res)
    TYPE(t_datetime)             :: res
    CLASS(t_event)               :: this
    TYPE(t_datetime), INTENT(IN) :: my_currentdatetime
    res = this%triggerNextEventDateTime
  END FUNCTION t_event_getNextOccurrenceDatetime

  FUNCTION t_event_getPrevOccurrenceDatetime(this, my_currentdatetime) RESULT(res)
    TYPE(t_datetime)             :: res
    CLASS(t_event)               :: this
    TYPE(t_datetime), INTENT(IN) :: my_currentdatetime
    res = this%triggeredPreviousEventDateTime
  END FUNCTION t_event_getPrevOccurrenceDatetime


  FUNCTION t_event_is_active(this, my_datetime, plus_slack, minus_slack) result(ret)
    CLASS(t_event)              :: this
    TYPE(t_datetime)            :: my_datetime
    TYPE(t_timedelta), OPTIONAL :: plus_slack
    TYPE(t_timedelta), OPTIONAL :: minus_slack
    logical(c_bool) :: ret

    !FIXME
    STOP
    !ret = isCurrentEventActive(this, my_datetime%dt, plus_slack%td, minus_slack%td)

  END FUNCTION t_event_is_active

  ! ================================================================================
  ! event group section:
  ! ================================================================================

  FUNCTION t_eventGroup_constructor(name) RESULT(this_event_group)
    TYPE(t_eventGroup) :: this_event_group
    CHARACTER(len=*), INTENT(in) :: name
    event_group_id                        = event_group_id + 1
    this_event_group%event_group_id       = event_group_id
    this_event_group%event_group_name     = name
    this_event_group%first_event_in_group => NULL()
    this_event_group%last_event_in_group  => NULL()
  END FUNCTION t_eventGroup_constructor

  SUBROUTINE t_eventGroup_addToGroup(this, event_to_add)
    CLASS (t_eventGroup) :: this
    TYPE(t_event), TARGET :: event_to_add
    IF (.NOT. ASSOCIATED(this%last_event_in_group)) THEN
      this%first_event_in_group => event_to_add
      NULLIFY(this%first_event_in_group%nextEventInGroup)
    ELSE
      this%last_event_in_group%nextEventInGroup => event_to_add
      NULLIFY(event_to_add%nextEventInGroup)
    ENDIF
    this%last_event_in_group => event_to_add
  END SUBROUTINE t_eventGroup_addToGroup

  FUNCTION t_eventGroup_getGroupId(this) RESULT(group_id)
    INTEGER :: group_id
    CLASS(t_eventGroup) ::this
    group_id = this%event_group_id
  END FUNCTION t_eventGroup_getGroupId

  FUNCTION t_eventGroup_getGroupName(this) RESULT(name)
    CHARACTER(len=max_groupname_str_len) :: name
    CLASS(t_eventGroup) :: this
    name = this%event_group_name
  END FUNCTION t_eventGroup_getGroupName

  FUNCTION t_eventGroup_getFirstEvent(this) RESULT(event_ptr)
    TYPE(t_event), POINTER :: event_ptr
    CLASS(t_eventGroup) :: this
    event_ptr => this%first_event_in_group
  END FUNCTION t_eventGroup_getFirstEvent
