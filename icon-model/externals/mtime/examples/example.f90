!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM example

  USE mtime

  IMPLICIT NONE

  CHARACTER(len=max_calendar_str_len)  :: calendar_in_use
  CHARACTER(len=max_date_str_len)      :: test_date_string
  CHARACTER(len=max_time_str_len)      :: test_time_string
  CHARACTER(len=max_datetime_str_len)  :: start_date_string
  CHARACTER(len=max_datetime_str_len)  :: stop_date_string
  CHARACTER(len=max_timedelta_str_len) :: time_step_string
  CHARACTER(len=max_datetime_str_len)  :: current_date_string
  CHARACTER(len=32) :: fmtstr = '%x'

  TYPE(date), POINTER :: test_date
  TYPE(time), POINTER :: test_time

  TYPE(datetime), POINTER :: start_date
  TYPE(datetime), POINTER :: start_date_tmp
  TYPE(datetime), POINTER :: stop_date
  TYPE(datetime), POINTER :: end_date

  TYPE(timedelta), POINTER :: time_step

  TYPE(datetime), POINTER :: current_date

  TYPE(datetime), POINTER :: tmp_date_test_1

  ! setup of calendar

  CALL setCalendar(proleptic_gregorian) ! proleptic_gregorian/year_of_365_days/year_of_360_days.
  CALL calendarToString(calendar_in_use)
  PRINT *, 'string: >', TRIM(calendar_in_use), '< int: ', calendarType(), ' (expect 2)'

  ! test global variables

  PRINT '(1x,a,i0,a)', 'global var: >', no_of_sec_in_a_day, '< (expect 86400)'

  ! test date interfacce

  test_date => newdate(2020, 01, 30)
  !test_date => newdate('2012-01-15')
  IF (ASSOCIATED(test_date)) PRINT *, 'allocated test_date'
  CALL dateToString(test_date, test_date_string)
  PRINT *, TRIM(test_date_string)
  CALL deallocateDate(test_date)
  IF (.NOT. (ASSOCIATED(test_date))) PRINT *, 'deallocated test_date'

  ! test time interface

  !test_time => newTime('12:13:49.654')
  test_time => newTime(12, 13, 49, 904)
  IF (ASSOCIATED(test_time)) PRINT *, 'allocated test_time'
  CALL timeToString(test_time, test_time_string)
  PRINT *, TRIM(test_time_string)
  CALL deallocateTime(test_time)
  IF (.NOT. (ASSOCIATED(test_time))) PRINT *, 'deallocated test_time'

  ! test datetime and timeddelta interface

  !start_date_tmp => newDatetime('2012-09-01T02:10:00.000')
  start_date_tmp => newDatetime(2012, 09, 01, 02, 10, 0, 0)
  start_date => newDatetime(start_date_tmp)
  IF (ASSOCIATED(start_date)) PRINT *, 'allocated start_date'
  CALL datetimetoString(start_date, start_date_string)
  PRINT *, 'Model start time : ', TRIM(start_date_string)

  stop_date => newDateTime("2012-09-10T14:00:00.000"); 
  CALL datetimeToPosixString(stop_date, stop_date_string, fmtstr)
  PRINT *, 'Model stop time  : ', TRIM(stop_date_string)

  time_step => newTimedelta('PT12H')
  CALL timedeltaToString(time_step, time_step_string)
  PRINT *, 'Model time step  : ', TRIM(time_step_string)

  current_date => newDatetime('2012-09-01T02:10:00.000')
  ! copy operator - overload newDatetime using construct and copy!
  ! current_date = start_date
  CALL datetimeToString(current_date, current_date_string)
  PRINT *, 'Model time       : ', TRIM(current_date_string)

  time_integration: DO
    current_date = current_date + time_step
    CALL datetimeToString(current_date, current_date_string)
    PRINT *, 'Model time loop  : ', TRIM(current_date_string)
    IF (current_date >= stop_date) EXIT time_integration
  END DO time_integration

  CALL datetimeToString(current_date, current_date_string)
  PRINT *, 'Model stop time  : ', TRIM(current_date_string)

  ! cleanup of objects

  CALL deallocateDatetime(start_date)
  CALL deallocateDatetime(stop_date)
  CALL deallocateDatetime(current_date)
  CALL deallocateTimeDelta(time_step)

  CALL icon_tests

  CALL event_tests

  ! reset calendar

  CALL resetCalendar()

CONTAINS

  SUBROUTINE icon_tests

    TYPE(timedelta), POINTER             :: mtime_td => NULL()
    TYPE(datetime), POINTER             :: mtime_date => NULL()
    TYPE(datetime), POINTER               :: dt1 => NULL(), dt2 => NULL()
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN)  :: td_string
    CHARACTER(len=MAX_DATETIME_STR_LEN)   :: dstring
    TYPE(timedelta), POINTER              :: time_delta
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   ::  dt_string

    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td*0.3d0
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT1H1M1S * 0.3 = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td*0.5d0
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT1H1M1S * 0.5 = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td*1.5d0
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT1H1M1S * 1.5 = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td*2.0d0
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT1H1M1S * 2.0 = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT98765S")
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT98765S = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT0S")
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT0S = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT987654321S")
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT987654321S = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT7536.3S")
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT7536.3S = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT7536.0003S")
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT7536.0003S (expect <null>) = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    mtime_date => newDatetime("1970-01-01T00:00:00")
    mtime_td => newTimedelta("PT987654321S")
    mtime_date = mtime_date + mtime_td
    CALL datetimetostring(mtime_date, dstring)
    WRITE (0, *) "1970-01-01T00:00:00 + PT987654321S = ", TRIM(dstring)
    CALL deallocateDatetime(mtime_date)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT10000000S")
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT10000000S = ", TRIM(td_string)
    mtime_date => newDatetime("2014-06-01T00:00:00")
    CALL datetimetostring(mtime_date, dstring)
    WRITE (0, *) "2014-06-01T00:00:00 = ", TRIM(dstring)
    mtime_date = mtime_date + mtime_td
    CALL datetimetostring(mtime_date, dstring)
    WRITE (0, *) "2014-06-01T00:00:00 + PT10000000S = ", TRIM(dstring)
    WRITE (0, *) "Expect 10000000 s = ", getTotalSecondsTimedelta(mtime_td, mtime_date), " s"
    CALL deallocateDatetime(mtime_date)
    CALL deallocateTimedelta(mtime_td)

    dt1 => newDatetime("1979-03-01T00:00:00.000")
    dt2 => newDatetime("1979-01-01T01:00:00.000")
    mtime_td => newTimedelta("PT0S")
    mtime_td = getTimeDeltaFromDateTime(dt1, dt2)
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT01M27D23H = ", TRIM(td_string)
    CALL deallocateDatetime(dt1)
    CALL deallocateDatetime(dt2)
    CALL deallocateTimedelta(mtime_td)

    dt1 => newDatetime("1979-07-01T00:00:00.000")
    dt2 => newDatetime("1979-01-01T01:00:00.000")
    mtime_td => newTimedelta("PT0S")
    mtime_td = getTimeDeltaFromDateTime(dt1, dt2)
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT05M29D23H = ", TRIM(td_string)
    CALL deallocateDatetime(dt1)
    CALL deallocateDatetime(dt2)
    CALL deallocateTimedelta(mtime_td)

    dt1 => newDatetime("1979-12-01T00:00:00.000")
    dt2 => newDatetime("1979-01-01T01:00:00.000")
    mtime_td => newTimedelta("PT0S")
    mtime_td = getTimeDeltaFromDateTime(dt1, dt2)
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT10M29D23H = ", TRIM(td_string)
    CALL deallocateDatetime(dt1)
    CALL deallocateDatetime(dt2)
    CALL deallocateTimedelta(mtime_td)

    dt1 => newDatetime("1980-01-01T00:00:00.000")
    dt2 => newDatetime("1979-01-01T01:00:00.000")
    mtime_td => newTimedelta("PT0S")
    mtime_td = getTimeDeltaFromDateTime(dt1, dt2)
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "PT11M30D23H = ", TRIM(td_string)
    CALL deallocateDatetime(dt1)
    CALL deallocateDatetime(dt2)
    CALL deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("-PT1H")
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "-PT01H = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)
    mtime_td => newTimedelta('-', 0, 0, 0, 1, 0, 0, 0)
    CALL timedeltatostring(mtime_td, td_string)
    WRITE (0, *) "-PT01H = ", TRIM(td_string)
    CALL deallocateTimedelta(mtime_td)

    CALL setCalendar(PROLEPTIC_GREGORIAN)
    start_date => newDatetime("2017-07-01T00:00:00.000")
    end_date => newDatetime("2017-07-31T00:00:00.000")
    time_delta => newTimeDelta("P01D")
    time_delta = end_date - start_date
    CALL datetimeToString(start_date, dt_string)
    WRITE (0, *) "start_date = ", TRIM(dt_string)
    CALL datetimeToString(end_date, dt_string)
    WRITE (0, *) "end_date   = ", TRIM(dt_string)
    CALL timedeltaToString(time_delta, td_string)
    WRITE (0, *) "difference (P30D) = ", TRIM(td_string)
    CALL deallocateDatetime(start_date)
    CALL deallocateDatetime(end_date)
    CALL deallocateTimedelta(time_delta)

    start_date => newDatetime("2017-07-01T00:00:00.000")
    end_date => newDatetime("2017-08-01T00:00:00.000")
    time_delta => newTimeDelta("P01D")
    time_delta = end_date - start_date
    CALL datetimeToString(start_date, dt_string)
    WRITE (0, *) "start_date = ", TRIM(dt_string)
    CALL datetimeToString(end_date, dt_string)
    WRITE (0, *) "end_date   = ", TRIM(dt_string)
    CALL timedeltaToString(time_delta, td_string)
    WRITE (0, *) "difference (P01M) = ", TRIM(td_string)
    CALL deallocateDatetime(start_date)
    CALL deallocateDatetime(end_date)
    CALL deallocateTimedelta(time_delta)

  END SUBROUTINE icon_tests

  SUBROUTINE event_tests

    TYPE(eventgroup), POINTER :: outputEventGroup
    TYPE(event), POINTER :: outputEvent
    TYPE(event), POINTER :: checkpointEvent
    TYPE(event), POINTER :: restartEvent
    TYPE(event), POINTER :: currentEvent
    TYPE(datetime), POINTER :: dtt
    TYPE(timedelta), POINTER :: tdd
    CHARACTER(len=max_eventname_str_len) :: currentEventString
    LOGICAL :: lret
    CHARACTER(len=max_eventname_str_len) :: aa
    CHARACTER(len=max_groupname_str_len) :: bb
    CHARACTER(len=max_datetime_str_len)  :: current_date_string_tmp

    outputEventGroup => newEventGroup('output driver')
    CALL getEventGroupName(outputEventGroup, aa)
    PRINT *, aa

    outputEvent => newEvent('output', '2000-01-01T00:00:00', '2010-01-01T00:00:01', '2013-01-01T00:00:02', 'PT06H')
    lret = addEventToEventGroup(outputEvent, outputEventGroup)

    dtt => getEventReferenceDateTime(outputEvent)
    CALL datetimeToString(dtt, current_date_string_tmp)
    PRINT *, TRIM(current_date_string_tmp)

    dtt => getEventFirstDateTime(outputEvent)
    CALL datetimeToString(dtt, current_date_string_tmp)
    PRINT *, TRIM(current_date_string_tmp)

    dtt => getEventLastDateTime(outputEvent)
    CALL datetimeToString(dtt, current_date_string_tmp)
    PRINT *, TRIM(current_date_string_tmp)

    tdd => getEventInterval(outputEvent)
    CALL timedeltaToString(tdd, current_date_string_tmp)
    PRINT *, TRIM(current_date_string_tmp)

    checkpointEvent => newEvent('checkpoint', '2010-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P05D')
    lret = addEventToEventGroup(checkpointEvent, outputEventGroup)

    restartEvent => newEvent('restart', '2000-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P01M')
    lret = addEventToEventGroup(restartEvent, outputEventGroup)

    currentEvent => getFirstEventFromEventGroup(outputEventGroup)

    PRINT *, 'Event list: '
    DO WHILE (ASSOCIATED(currentEvent))
      CALL getEventName(currentEvent, currentEventString)
      PRINT *, '   event: ', TRIM(currentEventString)
      currentEvent => getNextEventFromEventGroup(currentEvent)
    END DO

    PRINT *, 'HELLO', getEventId(restartEvent); 
    PRINT *, 'GOOGLE', getEventisFirstInMonth(outputEvent)

    tmp_date_test_1 => newDatetime('2000-01-01T01:00:00')
    CALL getTriggeredPreviousEventAtDateTime(checkpointEvent, tmp_date_test_1)
    CALL datetimeToString(tmp_date_test_1, current_date_string)
    PRINT *, current_date_string

    CALL getEventGroupName(outputEventGroup, bb); 
    PRINT *, bb

    CALL deallocateEventGroup(outputEventGroup)

  END SUBROUTINE event_tests

END PROGRAM example
