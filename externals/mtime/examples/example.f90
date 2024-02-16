!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program example

  use mtime

  implicit none

  character(len=max_calendar_str_len)  :: calendar_in_use
  character(len=max_date_str_len)      :: test_date_string
  character(len=max_time_str_len)      :: test_time_string
  character(len=max_datetime_str_len)  :: start_date_string
  character(len=max_datetime_str_len)  :: stop_date_string
  character(len=max_timedelta_str_len) :: time_step_string
  character(len=max_datetime_str_len)  :: current_date_string
  character(len=32) :: fmtstr = '%x'

  type(date), pointer :: test_date 
  type(date), pointer :: test_date_tmp
  type(time), pointer :: test_time 

  type(datetime), pointer :: start_date 
  type(datetime), pointer :: start_date_tmp
  type(datetime), pointer :: start_date_test
  type(datetime), pointer :: stop_date
  type(datetime), pointer :: end_date

  type(timedelta), pointer :: time_step

  type(datetime), pointer :: current_date 
 
  type(datetime), pointer :: current_date_test
  type(datetime), pointer :: tmp_date_test_1
  type(datetime), pointer :: dt

  ! setup of calendar

  call setCalendar(proleptic_gregorian) ! proleptic_gregorian/year_of_365_days/year_of_360_days.
  call calendarToString(calendar_in_use)
  print *, 'string: >', TRIM(calendar_in_use), '< int: ', calendarType(), ' (expect 2)'

  ! test global variables

  print '(1x,a,i0,a)', 'global var: >', no_of_sec_in_a_day, '< (expect 86400)'
  
  ! test date interfacce

  test_date => newdate(2020,01,30)
  !test_date => newdate('2012-01-15')
  if (associated(test_date)) print *, 'allocated test_date'
  call dateToString(test_date, test_date_string)
  print *, trim(test_date_string)
  call deallocateDate(test_date)
  if (.not. (associated(test_date))) print *, 'deallocated test_date'

  ! test time interface

  !test_time => newTime('12:13:49.654')
  test_time => newTime(12,13,49,904)
  if (associated(test_time)) print *, 'allocated test_time'
  call timeToString(test_time, test_time_string)
  print *, trim(test_time_string)
  call deallocateTime(test_time)
  if (.not. (associated(test_time))) print *, 'deallocated test_time'


  ! test datetime and timeddelta interface

  !start_date_tmp => newDatetime('2012-09-01T02:10:00.000')
  start_date_tmp => newDatetime(2012,09,01,02,10,0,0)
  start_date => newDatetime(start_date_tmp)
  if (associated(start_date)) print *, 'allocated start_date'
  call datetimetoString(start_date, start_date_string)
  print *, 'Model start time : ', trim(start_date_string)


  stop_date => newDateTime("2012-09-10T14:00:00.000");
  call datetimeToPosixString(stop_date, stop_date_string, fmtstr)
  print *, 'Model stop time  : ', trim(stop_date_string)

  time_step => newTimedelta('PT12H')
  call timedeltaToString(time_step, time_step_string)
  print *, 'Model time step  : ', trim(time_step_string)

  current_date => newDatetime('2012-09-01T02:10:00.000')
  ! copy operator - overload newDatetime using construct and copy!
  ! current_date = start_date
  call datetimeToString(current_date, current_date_string)
  print *, 'Model time       : ', trim(current_date_string)

  time_integration: do 
    current_date = current_date + time_step
    call datetimeToString(current_date, current_date_string)
    print *, 'Model time loop  : ', trim(current_date_string)
    if (current_date >= stop_date) exit time_integration
  enddo time_integration
  
  call datetimeToString(current_date, current_date_string)
  print *, 'Model stop time  : ', trim(current_date_string)

  ! cleanup of objects

  call deallocateDatetime(start_date)
  call deallocateDatetime(stop_date)
  call deallocateDatetime(current_date)
  call deallocateTimeDelta(time_step)

  call icon_tests
  
  call event_tests
  
  ! reset calendar

  call resetCalendar()

contains

  subroutine icon_tests

    type(timedelta),  pointer             :: mtime_td => null()
    type(datetime),   pointer             :: mtime_date => null()
    type(datetime), pointer               :: dt1 => null(), dt2 => null()    
    character(len=MAX_TIMEDELTA_STR_LEN)  :: td_string
    character(len=MAX_DATETIME_STR_LEN)   :: dstring
    TYPE(timedelta), POINTER              :: time_delta
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   ::  dt_string
    
    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td * 0.3d0
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT1H1M1S * 0.3 = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)
    
    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td * 0.5d0
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT1H1M1S * 0.5 = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td * 1.5d0
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT1H1M1S * 1.5 = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT1H1M1S")
    mtime_td = mtime_td * 2.0d0
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT1H1M1S * 2.0 = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT98765S")
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT98765S = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT0S")
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT0S = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)
    
    mtime_td => newTimedelta("PT987654321S")
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT987654321S = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)
    
    mtime_td => newTimedelta("PT7536.3S")
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT7536.3S = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)
    
    mtime_td => newTimedelta("PT7536.0003S")
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT7536.0003S (expect <null>) = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)
    
    mtime_date => newDatetime("1970-01-01T00:00:00")
    mtime_td => newTimedelta("PT987654321S")
    mtime_date = mtime_date + mtime_td
    call datetimetostring(mtime_date, dstring)
    write (0,*) "1970-01-01T00:00:00 + PT987654321S = ", TRIM(dstring)
    call deallocateDatetime(mtime_date)
    call deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("PT10000000S")
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT10000000S = ", TRIM(td_string)
    mtime_date => newDatetime("2014-06-01T00:00:00")
    call datetimetostring(mtime_date, dstring)
    write (0,*) "2014-06-01T00:00:00 = ", TRIM(dstring)
    mtime_date = mtime_date + mtime_td
    call datetimetostring(mtime_date, dstring)
    write (0,*) "2014-06-01T00:00:00 + PT10000000S = ", TRIM(dstring)
    write (0,*) "Expect 10000000 s = ", getTotalSecondsTimedelta(mtime_td, mtime_date), " s"
    call deallocateDatetime(mtime_date)
    call deallocateTimedelta(mtime_td)

    dt1 => newDatetime("1979-03-01T00:00:00.000")
    dt2 => newDatetime("1979-01-01T01:00:00.000")
    mtime_td => newTimedelta("PT0S")
    mtime_td = getTimeDeltaFromDateTime(dt1, dt2)
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT01M27D23H = ", TRIM(td_string)
    call deallocateDatetime(dt1)
    call deallocateDatetime(dt2)    
    call deallocateTimedelta(mtime_td)

    dt1 => newDatetime("1979-07-01T00:00:00.000")
    dt2 => newDatetime("1979-01-01T01:00:00.000")
    mtime_td => newTimedelta("PT0S")
    mtime_td = getTimeDeltaFromDateTime(dt1, dt2)
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT05M29D23H = ", TRIM(td_string)
    call deallocateDatetime(dt1)
    call deallocateDatetime(dt2)    
    call deallocateTimedelta(mtime_td)

    dt1 => newDatetime("1979-12-01T00:00:00.000")
    dt2 => newDatetime("1979-01-01T01:00:00.000")
    mtime_td => newTimedelta("PT0S")
    mtime_td = getTimeDeltaFromDateTime(dt1, dt2)
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT10M29D23H = ", TRIM(td_string)
    call deallocateDatetime(dt1)
    call deallocateDatetime(dt2)    
    call deallocateTimedelta(mtime_td)

    dt1 => newDatetime("1980-01-01T00:00:00.000")
    dt2 => newDatetime("1979-01-01T01:00:00.000")
    mtime_td => newTimedelta("PT0S")
    mtime_td = getTimeDeltaFromDateTime(dt1, dt2)
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "PT11M30D23H = ", TRIM(td_string)
    call deallocateDatetime(dt1)
    call deallocateDatetime(dt2)    
    call deallocateTimedelta(mtime_td)

    mtime_td => newTimedelta("-PT1H")
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "-PT01H = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)
    mtime_td => newTimedelta('-',0,0,0,1,0,0,0)
    call timedeltatostring(mtime_td, td_string)
    write (0,*) "-PT01H = ", TRIM(td_string)
    call deallocateTimedelta(mtime_td)
    
    CALL setCalendar(PROLEPTIC_GREGORIAN)
    start_date => newDatetime("2017-07-01T00:00:00.000")
    end_date   => newDatetime("2017-07-31T00:00:00.000")
    time_delta => newTimeDelta("P01D")
    time_delta = end_date - start_date
    CALL datetimeToString(start_date, dt_string)
    WRITE (0,*) "start_date = ", TRIM(dt_string)
    CALL datetimeToString(end_date, dt_string)
    WRITE (0,*) "end_date   = ", TRIM(dt_string)
    CALL timedeltaToString(time_delta, td_string)
    WRITE (0,*) "difference (P30D) = ", TRIM(td_string)
    CALL deallocateDatetime(start_date)
    CALL deallocateDatetime(end_date)
    CALL deallocateTimedelta(time_delta)

    start_date => newDatetime("2017-07-01T00:00:00.000")
    end_date   => newDatetime("2017-08-01T00:00:00.000")
    time_delta => newTimeDelta("P01D")
    time_delta = end_date - start_date
    CALL datetimeToString(start_date, dt_string)
    WRITE (0,*) "start_date = ", TRIM(dt_string)
    CALL datetimeToString(end_date, dt_string)
    WRITE (0,*) "end_date   = ", TRIM(dt_string)
    CALL timedeltaToString(time_delta, td_string)
    WRITE (0,*) "difference (P01M) = ", TRIM(td_string)
    CALL deallocateDatetime(start_date)
    CALL deallocateDatetime(end_date)
    CALL deallocateTimedelta(time_delta)

  end subroutine icon_tests

  subroutine event_tests

    type(eventgroup), pointer :: outputEventGroup
    type(event), pointer :: outputEvent
    type(event), pointer :: checkpointEvent
    type(event), pointer :: restartEvent
    type(event), pointer :: currentEvent
    type(datetime), pointer :: dtt
    type(timedelta), pointer :: tdd
    character(len=max_eventname_str_len) :: currentEventString
    logical :: lret
    character(len=max_eventname_str_len) :: aa
    character(len=max_groupname_str_len) :: bb
    character(len=max_datetime_str_len)  :: current_date_string_tmp

    outputEventGroup => newEventGroup('output driver')
    call getEventGroupName(outputEventGroup, aa)
    print *, aa

    outputEvent => newEvent('output', '2000-01-01T00:00:00', '2010-01-01T00:00:01', '2013-01-01T00:00:02', 'PT06H')
    lret = addEventToEventGroup(outputEvent, outputEventGroup)

    dtt => getEventReferenceDateTime(outputEvent) 
    call datetimeToString(dtt, current_date_string_tmp)
    print *, trim(current_date_string_tmp)

    dtt => getEventFirstDateTime(outputEvent)
    call datetimeToString(dtt, current_date_string_tmp)
    print *, trim(current_date_string_tmp)

    dtt => getEventLastDateTime(outputEvent)
    call datetimeToString(dtt, current_date_string_tmp)
    print *, trim(current_date_string_tmp)

    tdd => getEventInterval(outputEvent)
    call timedeltaToString(tdd, current_date_string_tmp)
    print *, trim(current_date_string_tmp)

    checkpointEvent => newEvent('checkpoint', '2010-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P05D')
    lret = addEventToEventGroup(checkpointEvent, outputEventGroup)

    restartEvent => newEvent('restart', '2000-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P01M')
    lret = addEventToEventGroup(restartEvent, outputEventGroup)

    currentEvent => getFirstEventFromEventGroup(outputEventGroup)

    print *, 'Event list: '
    do while (associated(currentEvent))
        call getEventName(currentEvent, currentEventString)
        print *,'   event: ', trim(currentEventString)
        currentEvent => getNextEventFromEventGroup(currentEvent)
    enddo

    print *,'HELLO' ,getEventId(restartEvent);
   
    print *, 'GOOGLE', getEventisFirstInMonth(outputEvent) 

    !type(datetime), pointer :: current_date_test
    current_date_test => newDatetime('2010-01-02T00:00:00')
    tmp_date_test_1 => newDatetime('2000-01-01T01:00:00')
    call getTriggeredPreviousEventAtDateTime(checkpointEvent, tmp_date_test_1)
    call datetimeToString(tmp_date_test_1, current_date_string)
    print *, current_date_string


    call getEventGroupName(outputEventGroup, bb);
    print *, bb
    

    call deallocateEventGroup(outputEventGroup)

  end subroutine event_tests

end program example
