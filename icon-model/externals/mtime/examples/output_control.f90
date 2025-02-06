!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM output_control

  USE mtime

  IMPLICIT NONE

  TYPE(eventgroup), POINTER :: outputEventGroup => NULL()
  TYPE(event), POINTER :: outputEvent1 => NULL()
  TYPE(event), POINTER :: outputEvent2 => NULL()
  TYPE(event), POINTER :: outputEvent3 => NULL()
  TYPE(event), POINTER :: outputEvent4 => NULL()

  CHARACTER(len=max_groupname_str_len) :: egstring

  CHARACTER(len=*), PARAMETER :: reference_date = '2000-01-01T00:00:00'
  CHARACTER(len=*), PARAMETER :: start_date = '2010-01-01T00:00:00'
  CHARACTER(len=*), PARAMETER :: end_date = '2010-02-01T00:00:00'

  INTEGER :: error
  LOGICAL :: lret

  TYPE event_namelist_handling
    ! represented in ISO 9601:2004 R string
    CHARACTER(len=32) :: steps
    CHARACTER(len=32) :: start_date
    CHARACTER(len=32) :: end_date
    CHARACTER(len=32) :: duration
    ! additional model required numbers
    CHARACTER(len=32) :: reference_date
    CHARACTER(len=32) :: offset
    CHARACTER(len=32) :: lbound_handling
    CHARACTER(len=32) :: ubound_handling
  END TYPE event_namelist_handling

  CALL setCalendar(proleptic_gregorian)

  outputEventGroup => newEventGroup('output driver')
  CALL getEventGroupName(outputEventGroup, egstring)
  PRINT *, TRIM(egstring)

  outputEvent1 => newEvent('output1', reference_date, start_date, end_date, 'PT6H', errno=error)
  CALL checkError(error)
  lret = addEventToEventGroup(outputEvent1, outputEventGroup)

  outputEvent2 => newEvent('output2', reference_date, start_date, end_date, 'PT2H', errno=error)
  lret = addEventToEventGroup(outputEvent2, outputEventGroup)

  outputEvent3 => newEvent('output3', reference_date, start_date, end_date, 'PT4H')
  lret = addEventToEventGroup(outputEvent3, outputEventGroup)

  outputEvent4 => newEvent('output4', reference_date, start_date, end_date, 'PT10H')
  lret = addEventToEventGroup(outputEvent4, outputEventGroup)

  CALL printEventTriggerTimes(outputEvent1)
  CALL printEventTriggerTimes(outputEvent2)
  CALL printEventTriggerTimes(outputEvent3)
  CALL printEventTriggerTimes(outputEvent4)

CONTAINS

  SUBROUTINE printEventTriggerTimes(eventI)
    TYPE(event), POINTER :: eventI
    TYPE(datetime), POINTER :: cd => NULL()
    TYPE(datetime), POINTER :: sd => NULL()
    TYPE(datetime), POINTER :: ed => NULL()
    TYPE(timedelta), POINTER :: dt => NULL()
    CHARACTER(len=max_datetime_str_len) :: tstring
    CHARACTER(len=max_timedelta_str_len) :: dstring
    INTEGER :: i = 0
    sd => getEventFirstDateTime(eventI)
    CALL datetimeToString(sd, tstring)
    PRINT *, 'Start time    : ', TRIM(tstring)
    ed => getEventLastDateTime(eventI)
    CALL datetimeToString(ed, tstring)
    PRINT *, 'End time      : ', TRIM(tstring)
    dt => getEventInterval(eventI)
    CALL timedeltaToString(dt, dstring)
    PRINT *, 'Time step     : ', TRIM(dstring)
    cd => newDatetime(sd)
    CALL datetimeToString(cd, dstring)
    PRINT *, 'Current time  : ', TRIM(dstring)
    output_times_loop: DO
      ! open lower bound set
      cd = cd + dt
      ! open upper bound set
      IF (cd > ed) EXIT output_times_loop
      CALL datetimeToString(cd, dstring)
      PRINT '(a,i3,a,a)', 'Output step ', i, ' time: ', TRIM(dstring)
      ! closed lower bound set
      !cd = cd + dt
      i = i + 1
      ! closed upper bound set
      ! if (cd > ed) exit output_times_loop
    END DO output_times_loop
    PRINT *, '________________________________________________________'
    PRINT *, ''
    CALL deallocateDatetime(sd)
    CALL deallocateDatetime(ed)
    CALL deallocateDatetime(cd)
    CALL deallocateTimedelta(dt)
  END SUBROUTINE printEventTriggerTimes

  SUBROUTINE checkError(errno)
    INTEGER, INTENT(in) :: errno
    CHARACTER(len=max_mtime_error_str_len) :: estring
    IF (errno /= no_error) THEN
      CALL mtime_strerror(errno, estring)
      PRINT *, TRIM(estring)
      STOP 'ERROR: finish.'
    END IF
  END SUBROUTINE checkError

END PROGRAM output_control
