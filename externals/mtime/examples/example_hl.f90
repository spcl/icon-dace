!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM example

#ifndef __NVCOMPILER

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_int64_t
  USE mtime, ONLY: setCalendar, PROLEPTIC_GREGORIAN
  USE mtime_hl

  IMPLICIT NONE

  TYPE(t_datetime)  :: dt, dt2, dt3, dt4, dt5
  TYPE(t_timedelta) :: td, td1, td2
  TYPE(t_juliandelta) :: jd
  TYPE(t_julianday)   :: jday
  TYPE(divisionquotienttimespan) :: dqts

  CHARACTER(LEN=*), PARAMETER :: ERR_UNCAUGHT = "!!!!!!!!! ERROR WAS NOT CAUGHT !!!!!!!"
  INTEGER       :: test_number1, test_number2, test_result
  LOGICAL       :: lerror
  CHARACTER(c_char) :: c_sign

  WRITE (0, *) "example_hl : test example"

  WRITE (0, *) "  setting calendar."
  CALL setCalendar(PROLEPTIC_GREGORIAN)

  ! Instead of:  TYPE(datetime), POINTER :: dt
  !              dt => newDatetime("1970-01-01T00:00:00")
  !                [...]
  !              CALL deallocateDatetime(dt)

  ! TYPE(t_datetime)  :: dt
  WRITE (0, *) "  testing assignment of t_datetime."
  dt = t_datetime("1970-01-01T00:00:00")

  ! Instead of:  TYPE(timedelta), POINTER :: td
  !              td => newTimeDelta("PT1H1M1S")
  !                [...]
  !              CALL deallocateTimedelta(td)

  ! TYPE(t_timedelta) :: td
  WRITE (0, *) "  testing assignment of t_timedelta."
  td = t_timedelta("PT1H1M1S")

  ! Instead of:  CALL datetimeToString(dt, dstring1)
  !              CALL timedeltaToString(td, dstring2)
  !              WRITE (0,*) "compute ", dstring1, " + ", dstring2

  WRITE (0, *) "compute ", dt%toString(), " + ", td%toString()

  ! Adding time intervals is unchanged:

  dt = dt + td
  WRITE (0, *) "result: ", dt%toString()

  ! --- Further examples:

  ! subtraction of two dates
  dt2 = t_datetime("1970-01-01T00:00:00")
  td = dt - dt2
  WRITE (0, *) "subtraction of dates: time delta: ", td%toString()

  ! comparison of dates
  dt = t_datetime("1970-01-01T00:00:00")
  dt2 = t_datetime("1970-01-01T00:00:00")
  WRITE (0, *) dt%toString(), " == ", dt2%toString(), ": ", (dt == dt2)
  dt3 = t_datetime("1970-01-01T00:00:01")
  WRITE (0, *) dt%toString(), " == ", dt3%toString(), ": ", (dt == dt3)

  ! min / max
  dt4 = MIN(dt2, dt3)
  WRITE (0, *) "MIN(", dt2%toString(), ", ", dt3%toString(), "): ", dt4%toString()
  dt4 = MAX(dt2, dt3)
  WRITE (0, *) "MAX(", dt2%toString(), ", ", dt3%toString(), "): ", dt4%toString()

  ! interval assignment with milliseconds
  td = t_timedeltaFromMilliseconds(360000)
  WRITE (0, *) "interval assignment with milliseconds: ", td%toString()

  td = t_timedelta("PT1H")
  td = td*0.5_c_double
  WRITE (0, *) "PT1H * 0.5 = ", td%toString()
  td = t_timedelta("PT1H")
  td = 0.5_c_double*td
  WRITE (0, *) "PT1H * 0.5 = ", td%toString()

  td = t_timedelta("PT1H")
  td = td*2
  WRITE (0, *) "PT1H * 2 = ", td%toString()
  td = t_timedelta("PT1H")
  td = 2*td
  WRITE (0, *) "PT1H * 2 = ", td%toString()

  ! division
  td1 = t_timedelta("PT23H42M")
  dqts = td1%divideInSecondsBy(td)
  WRITE (0, *) td1%toString(), " / ", td%toString(), " = ", dqts%quotient, " plus stuff"

  ! toSeconds, toMilliSeconds
  WRITE (0, *) td%toString(), " is in seconds ", td%toSeconds(dt)
  WRITE (0, *) td%toString(), " is in milliseconds ", td%toMilliSeconds(dt)

  dt = t_datetime("1970-01-01T00:00:00")
  jday = dt%toJulianDay()
  dt2 = jday%toDateTime()
  PRINT *, dt%toString(), " = ", dt2%toString()

  ! register an error callback without stopping the application for
  ! our tests:
  CALL register_finish_mtime_procedure(error_callback)

  ! produce errors
  WRITE (0, *) 'error test: dt = t_datetime("1970--01-01T00:00:00")'
  lerror = .FALSE.
  dt = t_datetime("1970--01-01T00:00:00")
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  ! The following test cannot be handled since we removed
  ! the explicit "assign" implementation:
  !
  ! WRITE (0,*) 'error test: dt = dt5'
  ! lerror = .FALSE.
  ! dt = dt5
  ! IF (.NOT. lerror)  WRITE(0,*) ERR_UNCAUGHT

  WRITE (0, *) 'error test: dt5%toString()'
  lerror = .FALSE.
  WRITE (0, *) dt5%toString()
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  WRITE (0, *) 'error test: dt5%to_posix_string()'
  lerror = .FALSE.
  WRITE (0, *) dt5%toString("%s%d%LK")
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  WRITE (0, *) 'error test: dt = dt + td2'
  lerror = .FALSE.
  dt = dt + td2
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  WRITE (0, *) 'error test: dt = dt - td2'
  lerror = .FALSE.
  dt = dt - td2
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  WRITE (0, *) 'error test: dt = dt - td2'
  lerror = .FALSE.
  dt = dt - td2
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  WRITE (0, *) 'error test: td = t_timedelta(...)'
  lerror = .FALSE.
  td = t_timedelta("P1lK")
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  td = t_timedeltaFromMilliseconds(HUGE(INT(1)))
  WRITE (0, *) "t_timedeltaFromMilliseconds(HUGE(INT(1))) : HUGE(INT(1)) = ", HUGE(INT(1))
  WRITE (0, *) "td%toString() = ", td%toString()

  WRITE (0, *) 'error test: td = td * 0.000001D0'
  WRITE (0, *) 'td%toString() = ', td%toString()
  lerror = .FALSE.
  td = td*0.000001d0
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  ! FIXME: This does no longer work, probably it initially worked for the wrong
  ! reasons?
  WRITE (0, *) 'error test: td2 = td2 * 1'
  lerror = .FALSE.
  td2 = td2*1
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  ! FIXME: This does no longer work, probably it initially worked for the wrong
  ! reasons?
  WRITE (0, *) 'error test: td2%toString()'
  lerror = .FALSE.
  WRITE (0, *) 'td2%toString() = ', td2%toString()
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  ! --------------------------------------------------------------------------------
  ! Test cases for t_juliandelta

  c_sign = '+'
  jd = t_juliandelta(sign=c_sign, day=99_c_int64_t, ms=123_c_int64_t)
  WRITE (0, *) "error test: jd = t_juliandelta(sign='r', day='99', ms='123')"
  lerror = .FALSE.
  c_sign = 'r'
  jd = t_juliandelta(sign=c_sign, day=99_c_int64_t, ms=123_c_int64_t)
  IF (.NOT. lerror) WRITE (0, *) ERR_UNCAUGHT

  CALL event_tests
CONTAINS
  SUBROUTINE event_tests

    TYPE(t_eventgroup) :: outputEventGroup
    TYPE(t_event) :: outputEvent
    TYPE(t_event) :: checkpointEvent
    TYPE(t_event) :: restartEvent
    TYPE(t_event), POINTER :: currentEvent
    TYPE(t_datetime) :: dtt
    TYPE(t_timedelta) :: tdd
    CHARACTER(len=max_eventname_str_len) :: currentEventString
    LOGICAL :: lret
    CHARACTER(len=max_eventname_str_len) :: aa
    CHARACTER(len=max_groupname_str_len) :: bb
    CHARACTER(len=max_datetime_str_len)  :: current_date_string_tmp

    outputEventGroup = t_eventGroup('output driver')
    aa = outputEventGroup%getName()
    PRINT *, aa

    outputEvent = t_event('output', '2000-01-01T00:00:00', '2010-01-01T00:00:01', '2013-01-01T00:00:02', 'PT06H')
    CALL outputEventGroup%append(outputEvent)

!    dtt => getEventReferenceDateTime(outputEvent)
!    call datetimeToString(dtt, current_date_string_tmp)
!    print *, trim(current_date_string_tmp)

    dtt = outputEvent%getFirstDateTime()
    current_date_string_tmp = dtt%toString()
    PRINT *, TRIM(current_date_string_tmp)

    dtt = outputEvent%getLastDateTime()
    current_date_string_tmp = dtt%toString()
    PRINT *, TRIM(current_date_string_tmp)

    tdd = outputEvent%getInterval()
    current_date_string_tmp = tdd%toString()
    PRINT *, TRIM(current_date_string_tmp)

    checkpointEvent = t_event('checkpoint', '2010-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P05D')
    CALL outputEventGroup%append(checkpointEvent)

    restartEvent = t_event('restart', '2000-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P01M')
    CALL outputEventGroup%append(restartEvent)

    currentEvent => outputEventGroup%getFirstEvent()
    PRINT *, 'Event list: '
    DO WHILE (ASSOCIATED(currentEvent))
      currentEventString = currentEvent%getName()
      PRINT *, '   event: ', TRIM(currentEventString)
      ! AHHHHH
      currentEvent => currentEvent%nextEvent()
    END DO

    !FIXME
!    print *,'HELLO', restartEventgetEvent);
!    print *, 'GOOGLE', getEventisFirstInMonth(outputEvent)

!    type(t_datetime) :: current_date_test
!    current_date_test = t_datetime('2010-01-02T00:00:00')
!    tmp_date_test_1 = t_datetime('2000-01-01T01:00:00')
!    call getTriggeredPreviousEventAtDateTime(checkpointEvent, tmp_date_test_1)
!    call datetimeToString(tmp_date_test_1, current_date_string)
!    print *, current_date_string

    bb = outputEventGroup%getName()
    PRINT *, bb

!    call deallocateEventGroup(outputEventGroup)

  END SUBROUTINE event_tests

  SUBROUTINE error_callback(leading_text, message_text)
    CHARACTER(len=*), INTENT(in) :: leading_text
    CHARACTER(len=*), INTENT(in) :: message_text
    WRITE (0, *) TRIM(leading_text), ": ", TRIM(message_text)
    lerror = .TRUE.
  END SUBROUTINE error_callback

#endif

END PROGRAM example
