!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
module mtime_c_bindings
  !
  use, intrinsic :: iso_c_binding, only: c_bool, c_int, c_char, c_null_char, c_ptr, c_int64_t, &
       &                                 c_float, c_double, c_associated
  !
  use mtime_constants
  use mtime_error_handling
  !
  implicit none
  !
  public
  !
  ! Type, bind(c)
  !
  type, bind(c) :: event
    ! FIXME (some day in the future): This derived type should not specify both -
    ! the linked list element and the element itself.
    type(c_ptr) :: nextEventInGroup
    integer(c_int64_t) :: eventId
    type(c_ptr) :: eventName
    type(c_ptr) :: eventsLastEvaluationDateTime
    type(c_ptr) :: eventReferenceDatetime
    type(c_ptr) :: eventFirstDatetime
    type(c_ptr) :: eventLastDatetime
    type(c_ptr) :: eventInterval
    type(c_ptr) :: eventOffset
    logical(c_bool) :: neverTriggerEvent
    logical(c_bool) :: triggerCurrentEvent
    logical(c_bool) :: nextEventIsFirst
    logical(c_bool) :: lastEventWasFinal
    logical(c_bool) :: eventisFirstInDay
    logical(c_bool) :: eventisFirstInMonth
    logical(c_bool) :: eventisFirstInYear
    logical(c_bool) :: eventisLastInDay
    logical(c_bool) :: eventisLastInMonth
    logical(c_bool) :: eventisLastInYear
    type(c_ptr) :: triggerNextEventDateTime
    type(c_ptr) :: triggeredPreviousEventDateTime
  end type event
  !
  type, bind(c) :: eventgroup
    integer(c_int64_t) :: eventGroupId
    type(c_ptr) :: eventGroupName
    type(c_ptr) :: firstEventInGroup
  end type eventgroup
  !
  type, bind(c) :: julianday
    integer(c_int64_t) :: day  !< the actual Julian day
    integer(c_int64_t) :: ms   !< the milisecond on that particular day
  end type julianday
  !
  type, bind(c) :: juliandelta
    character(c_char)  :: sign
    integer(c_int64_t) :: day
    integer(c_int64_t) :: ms
  end type juliandelta
  !
  type, bind(c) :: date
    integer(c_int64_t) :: year
    integer(c_int) :: month
    integer(c_int) :: day
  end type date
  !
  type, bind(c) :: time
    integer(c_int) :: hour
    integer(c_int) :: minute
    integer(c_int) :: second
    integer(c_int) :: ms
  end type time
  !
  type, bind(c) :: datetime
    type(date) :: date
    type(time) :: time
  end type datetime
  !
  type, bind(c) :: timedelta
    integer(c_int) :: flag_std_form
    character(c_char) :: sign
    integer(c_int64_t) :: year
    integer(c_int) :: month
    integer(c_int) :: day
    integer(c_int) :: hour
    integer(c_int) :: minute
    integer(c_int) :: second
    integer(c_int) :: ms
  end type timedelta
  !
  type, bind(c) :: divisionquotienttimespan
    integer(c_int64_t) :: quotient;
    integer(c_int64_t) :: remainder_in_ms;
  end type divisionquotienttimespan
  !
  ! End Type, bind(c)
  !
  interface
    !
    subroutine setCalendar(ct) bind(c, name='initCalendar') !TESTED-OK
      import :: c_int
      integer(c_int), value :: ct
    end subroutine setCalendar
    !
    subroutine resetCalendar() bind(c, name='freeCalendar') !TESTED-OK
    end subroutine resetCalendar
    !
    function calendarType() bind(c, name='getCalendarType') !TESTED-OK
      import :: c_int
      integer(c_int) :: calendarType
    end function calendarType
    !
    function my_calendartostring(calendar) result(c_pointer) bind(c, name='calendarToString') !TESTED-OK
      import :: c_char, c_ptr
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: calendar
    end function my_calendartostring
    !
  end interface

  interface
    function my_newjuliandelta(sign, day, ms) result(c_pointer) bind(c, name='newJulianDelta')
      import :: c_int64_t, c_char, c_ptr
      type(c_ptr)               :: c_pointer
      character(c_char),  value :: sign
      integer(c_int64_t), value :: day
      integer(c_int64_t), value :: ms
    end function my_newjuliandelta
    !
    subroutine my_deallocatejuliandelta(jd) bind(c,name='deallocateJulianDelta')
      import :: c_ptr
      type(c_ptr), value :: jd
    end subroutine my_deallocatejuliandelta
  end interface

  interface
    function my_newjulianday(day, ms) result(c_pointer) bind(c, name='newJulianDay')
      import :: c_int64_t, c_ptr
      type(c_ptr) :: c_pointer
      integer(c_int64_t), value :: day
      integer(c_int64_t), value :: ms
    end function my_newjulianday
    !
    subroutine my_deallocatejulianday(jd) bind(c,name='deallocateJulianDay')
      import :: c_ptr
      type(c_ptr), value :: jd
    end subroutine my_deallocatejulianday
    !
    function my_replacejulianday(src, dest) result(ret_dest) bind(c, name='replaceJulianday')
      import :: c_ptr
      type(c_ptr) :: ret_dest
      type(c_ptr), value :: src
      type(c_ptr), value :: dest
    end function my_replacejulianday
    !
    function my_comparejulianday(op1, op2) result(ret) bind(c, name='compareJulianDay')
      import :: c_ptr, c_int
      integer(c_int) :: ret
      type(c_ptr), value :: op1
      type(c_ptr), value :: op2
    end function my_comparejulianday
    !
    function my_juliandaytostring(my_julianday, string) result(string_ptr) bind(c, name='juliandayToString')
      import :: c_ptr, c_char
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_julianday
      character(c_char), dimension(*) :: string
    end function my_juliandaytostring
    !
    function my_addjuliandeltatojulianday(my_julianday, my_juliandelta, ret_julianday) result(julianday_ptr) &
         &   bind(c, name='addJulianDeltaToJulianDay')
      import :: c_ptr
      type(c_ptr) :: julianday_ptr
      type(c_ptr), value :: my_julianday, my_juliandelta, ret_julianday
    end function my_addjuliandeltatojulianday
    !
    function my_substractjulianday(my_julianday1, my_julianday2, ret_juliandelta) result(juliandelta_ptr) &
         &   bind(c, name='substractJulianDay')
      import :: c_ptr
      type(c_ptr) :: juliandelta_ptr
      type(c_ptr), value :: my_julianday1, my_julianday2, ret_juliandelta
    end function my_substractjulianday
    !
  end interface
  !
  interface
    !
    function my_newdatefromstring(string) result(c_pointer) bind(c, name='newDate')
      import :: c_char, c_ptr
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: string
    end function my_newdatefromstring
    !
    function my_newrawdate(year, month, day) result(c_pointer) bind(c, name='newRawDate')
      import :: c_int64_t, c_int, c_ptr
      type(c_ptr) :: c_pointer
      integer(c_int64_t), value :: year
      integer(c_int), value :: month, day
    end function my_newrawdate
    !
    function my_constructandcopydate(d) result(c_pointer) bind(c,name='constructAndCopyDate')
      import :: c_ptr
      type(c_ptr) :: c_pointer
      type(c_ptr), value :: d
    end function my_constructandcopydate
    !
    subroutine my_deallocatedate(d) bind(c, name='deallocateDate')
      import :: c_ptr
      type(c_ptr), value :: d
    end subroutine my_deallocatedate
    !
    function my_replacedate(src, dest) result(ret_dest) bind(c, name='replaceDate')
      import :: c_ptr
      type(c_ptr) :: ret_dest
      type(c_ptr), value :: src , dest
    end function my_replacedate
    !
    function my_datetostring(my_date, string) result(string_ptr) bind(c, name='dateToString')
      import :: c_ptr, c_char
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_date
      character(c_char), dimension(*) :: string
    end function my_datetostring
    !
    function my_datetoposixstring(my_date, string, fmtstr) result(string_ptr) bind(c, name='dateToPosixString')
      import :: c_ptr, c_char
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_date
      character(c_char), dimension(*) :: string
      character(c_char), dimension(*) :: fmtstr
    end function my_datetoposixstring
    !
  end interface
  !
  interface
    !
    function my_newtimefromstring(string) result(c_pointer) bind(c, name='newTime')
      import :: c_char, c_ptr
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: string
    end function my_newtimefromstring
    !
    function my_newrawtime(hour, minute, second, ms) result(c_pointer) bind(c, name='newRawTime')
      import :: c_int, c_ptr
      type(c_ptr) :: c_pointer
      integer(c_int), value :: hour, minute, second, ms
    end function my_newrawtime
    !
    function my_constructandcopytime(t) result(c_pointer) bind(c,name='constructAndCopyTime')
      import :: c_ptr
      type(c_ptr) :: c_pointer
      type(c_ptr), value :: t
    end function my_constructandcopytime
    !
    subroutine my_deallocatetime(t) bind(c,name='deallocateTime')
      import :: c_ptr
      type(c_ptr), value :: t
    end subroutine my_deallocatetime
    !
    function my_replacetime(src, dest) result(ret_dest) bind(c, name='replaceTime')
      import :: c_ptr
      type(c_ptr) :: ret_dest
      type(c_ptr), value :: src , dest
    end function my_replacetime
    !
    function my_timetostring(my_time, string) result(string_ptr) bind(c, name='timeToString')
      import :: c_ptr, c_char
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_time
      character(c_char), dimension(*) :: string
    end function my_timetostring
    !
    function my_timetoposixstring(my_time, string, fmtstr) result(string_ptr) bind(c, name='timeToPosixString')
      import :: c_ptr, c_char
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_time
      character(c_char), dimension(*) :: string
      character(c_char), dimension(*) :: fmtstr
    end function my_timetoposixstring
    !
  end interface
  !
  interface
    !
    function my_newdatetime(string) result(c_pointer) bind(c, name='newDateTime')
      import :: c_ptr, c_char
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: string
    end function my_newdatetime
    !
    function my_newrawdatetime(year, month, day, hour, minute, second, ms) result(c_pointer) &
         &             bind(c, name='newRawDateTime')
      import :: c_int64_t, c_int, c_ptr
      type(c_ptr) :: c_pointer
      integer(c_int64_t), value :: year
      integer(c_int), value :: month, day, hour, minute, second, ms
    end function my_newrawdatetime
    !
    function my_constructandcopydatetime(dt) result(c_pointer) bind(c,name='constructAndCopyDateTime')
      import :: c_ptr
      type(c_ptr) :: c_pointer
      type(c_ptr), value :: dt
    end function my_constructandcopydatetime
    !
    subroutine my_deallocatedatetime(dt) bind(c,name='deallocateDateTime')
      import :: c_ptr
      type(c_ptr), value :: dt
    end subroutine my_deallocatedatetime
    !
    function my_comparedatetime(op1, op2) result(ret) bind(c, name='compareDatetime')
      import :: c_ptr, c_int
      integer(c_int) :: ret
      type(c_ptr), value :: op1
      type(c_ptr), value :: op2
    end function my_comparedatetime

    function my_replacedatetime(src, dest) result(ret_dest) bind(c, name='replaceDatetime')
      import :: c_ptr
      type(c_ptr) :: ret_dest
      type(c_ptr), value :: src
      type(c_ptr), value :: dest
    end function my_replacedatetime
    !
    function my_datetimetostring(my_time, string) result(string_ptr) bind(c, name='datetimeToString')
      import :: c_ptr, c_char
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_time
      character(c_char), dimension(*) :: string
    end function my_datetimetostring
    !
    function my_datetimetoposixstring(my_time, string, fmtstr) result(string_ptr) bind(c, name='datetimeToPosixString')
      import :: c_ptr, c_char
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_time
      character(c_char), dimension(*) :: string
      character(c_char), dimension(*) :: fmtstr
    end function my_datetimetoposixstring
    !
    function my_getnoofdaysinmonthdatetime(dt) bind(c, name='getNoOfDaysInMonthDateTime')
      import :: c_int, c_ptr
      integer(c_int) :: my_getnoofdaysinmonthdatetime
      type(c_ptr), value :: dt
    end function my_getnoofdaysinmonthdatetime
    !
    function my_getnoofdaysinyeardatetime(dt) bind(c, name='getNoOfDaysInYearDateTime')
      import :: c_int, c_ptr
      integer(c_int) :: my_getnoofdaysinyeardatetime
      type(c_ptr), value :: dt
    end function my_getnoofdaysinyeardatetime
    !
    function my_getdayofyearfromdatetime(dt) bind(c, name='getDayOfYearFromDateTime')
      import :: c_int, c_ptr
      integer(c_int) :: my_getdayofyearfromdatetime
      type(c_ptr), value :: dt
    end function my_getdayofyearfromdatetime
    !
    function my_getnoofsecondselapsedinmonthdatetime(dt) bind(c, name='getNoOfSecondsElapsedInMonthDateTime')
      import :: c_int64_t, c_ptr
      integer(c_int64_t) :: my_getnoofsecondselapsedinmonthdatetime
      type(c_ptr), value :: dt
    end function my_getnoofsecondselapsedinmonthdatetime
    !
    function my_getnoofsecondselapsedindaydatetime(dt) bind(c, name='getNoOfSecondsElapsedInDayDateTime')
      import :: c_int, c_ptr
      integer(c_int) :: my_getnoofsecondselapsedindaydatetime
      type(c_ptr), value :: dt
    end function my_getnoofsecondselapsedindaydatetime
    !
    function my_getjuliandayfromdatetime(dt, jd) bind(c, name='getJulianDayFromDateTime')
      import :: c_ptr
      type(c_ptr) :: my_getjuliandayfromdatetime
      type(c_ptr), value :: dt, jd
    end function my_getjuliandayfromdatetime
    !
    function my_getdatetimefromjulianday(jd, dt) bind(c, name='getDateTimeFromJulianDay')
      import :: c_ptr
      type(c_ptr) :: my_getdatetimefromjulianday
      type(c_ptr), value :: jd, dt
    end function my_getdatetimefromjulianday
    !
  end interface
  !
  interface
    !
    function my_newtimedeltafromstring(string) result(c_pointer) bind(c, name='newTimeDelta')
      import :: c_char, c_ptr
      type(c_ptr)                     :: c_pointer
      character(c_char), dimension(*) :: string
    end function my_newtimedeltafromstring
    !
    function my_newrawtimedelta(sign, year, month, day, hour, minute, second, ms) result(c_pointer) &
         &             bind(c, name='newRawTimeDelta')
      import :: c_int64_t, c_char, c_int, c_ptr
      type(c_ptr) :: c_pointer
      character(c_char), value :: sign
      integer(c_int64_t), value :: year
      integer(c_int), value :: month, day, hour, minute, second, ms
    end function my_newrawtimedelta
    !
    function my_constructandcopytimedelta(td) result(c_pointer) bind(c,name='constructAndCopyTimeDelta')
      import :: c_ptr
      type(c_ptr) :: c_pointer
      type(c_ptr), value :: td
    end function my_constructandcopytimedelta
    !
    subroutine my_deallocatetimedelta(dt) bind(c, name='deallocateTimeDelta')
      import :: c_ptr
      type(c_ptr), value :: dt
    end subroutine my_deallocatetimedelta
    !
    function my_comparetimedelta(op1, op2) result(ret) bind(c, name='compareTimeDelta')
      import :: c_ptr, c_int
      integer(c_int) :: ret
      type(c_ptr), value :: op1
      type(c_ptr), value :: op2
    end function my_comparetimedelta
    !
    function my_gettimedeltafromdate(my_date1,my_date2,timedelta_return) result(timedelta_ptr) &
         &    bind(c,name='getTimeDeltaFromDate')
      import :: c_ptr
      type(c_ptr) :: timedelta_ptr
      type(c_ptr), value :: my_date1, my_date2, timedelta_return
    end function my_gettimedeltafromdate
    !
    function my_gettimedeltafromdatetime(my_datetime1,my_datetime2,timedelta_return) result(timedelta_ptr) &
         &    bind(c,name='getTimeDeltaFromDateTime')
      import :: c_ptr
      type(c_ptr) :: timedelta_ptr
      type(c_ptr), value :: my_datetime1, my_datetime2, timedelta_return
    end function my_gettimedeltafromdatetime
    !
    function my_gettotalmillisecondstimedelta(my_timedelta, my_datetime) &
         &   bind(c, name='getTotalMilliSecondsTimeDelta')
      import :: c_ptr, c_int64_t
      integer(c_int64_t) :: my_gettotalmillisecondstimedelta
      type(c_ptr), value :: my_timedelta, my_datetime
    end function my_gettotalmillisecondstimedelta
    !
    function my_gettotalsecondstimedelta(my_timedelta, my_datetime) bind(c, name='getTotalSecondsTimeDelta')
      import :: c_ptr, c_int64_t
      integer(c_int64_t) :: my_gettotalsecondstimedelta
      type(c_ptr), value :: my_timedelta, my_datetime
    end function my_gettotalsecondstimedelta
    !
    function my_timedeltatostring(td, tostring) result(string_ptr) bind(c, name='timedeltaToString')
      import :: c_ptr, c_char
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: td
      character(c_char), dimension(*) :: tostring
    end function my_timedeltatostring
    !
    function my_addtimedeltatodatetime(my_datetime, my_timedelta, ret_datetime) result(datetime_ptr) &
         &   bind(c, name='addTimeDeltaToDateTime')
      import :: c_ptr
      type(c_ptr) :: datetime_ptr
      type(c_ptr), value :: my_datetime, my_timedelta, ret_datetime
    end function my_addtimedeltatodatetime
    !
    function my_addtimedeltatodate(my_date, my_timedelta, ret_date) result(date_ptr) &
         &   bind(c, name='addTimeDeltaToDate')
      import :: c_ptr
      type(c_ptr) :: date_ptr
      type(c_ptr), value :: my_date, my_timedelta, ret_date
    end function my_addtimedeltatodate
    !
    function my_elementwisescalarmultiplytimedelta(my_timedelta, lambda, scaled_timedelta) result(timedelta_return) &
         &   bind(c, name='elementwiseScalarMultiplyTimeDelta')
      import :: c_ptr, c_int64_t
      type(c_ptr) ::timedelta_return
      type(c_ptr), value :: my_timedelta, scaled_timedelta
      integer(c_int64_t), value :: lambda
    end function my_elementwisescalarmultiplytimedelta
    !
    function my_elementwisescalarmultiplytimedeltadp(my_timedelta, lambda, scaled_timedelta) result(timedelta_return) &
         &   bind(c, name='elementwiseScalarMultiplyTimeDeltaDP')
      import :: c_ptr, c_double
      type(c_ptr) ::timedelta_return
      type(c_ptr), value :: my_timedelta, scaled_timedelta
      real(c_double), value :: lambda
    end function my_elementwisescalarmultiplytimedeltadp
    !
    function my_elementwiseAddTimeDeltatoTimeDelta(my_timedelta1,my_timedelta2,added_timedelta) result(timedelta_return) &
         &   bind(c, name='elementwiseAddTimeDeltatoTimeDelta')
      import :: c_ptr
      type(c_ptr) :: timedelta_return
      type(c_ptr), value :: my_timedelta1,my_timedelta2,added_timedelta
    end function my_elementwiseAddTimeDeltatoTimeDelta
    !
    function my_modulotimedelta(a, p, quot) result(rem) bind(c, name='moduloTimedelta')
      import :: c_ptr, c_int64_t
      integer(c_int64_t) :: rem
      type(c_ptr), value :: a
      type(c_ptr), value :: p
      type(c_ptr), value :: quot
    end function my_modulotimedelta
    !
    function my_getptstringfromms(ms, tostring) result(string_ptr) bind(c, name='getPTStringFromMS')
      import :: c_ptr, c_int64_t, c_char
      type(c_ptr) :: string_ptr
      integer(c_int64_t), value :: ms
      character(c_char), dimension(*) :: tostring
    end function my_getptstringfromms
    !
    function my_getptstringfromsecondsint(s, tostring) result(string_ptr) bind(c, name='getPTStringFromSeconds')
      import :: c_ptr, c_int64_t, c_char
      type(c_ptr) :: string_ptr
      integer(c_int64_t), value :: s
      character(c_char), dimension(*) :: tostring
    end function my_getptstringfromsecondsint
    !
    function my_getptstringfromsecondsfloat(s, tostring) result(string_ptr) bind(c, name='getPTStringFromSecondsFloat')
      import :: c_ptr, c_float, c_char
      type(c_ptr) :: string_ptr
      real(c_float), value :: s
      character(c_char), dimension(*) :: tostring
    end function my_getptstringfromsecondsfloat
    !
    function my_getptstringfromsecondsdouble(s, tostring) result(string_ptr) bind(c, name='getPTStringFromSecondsDouble')
      import :: c_ptr, c_double, c_char
      type(c_ptr) :: string_ptr
      real(c_double), value :: s
      character(c_char), dimension(*) :: tostring
    end function my_getptstringfromsecondsdouble
    !
    function my_getptstringfromminutes(m, tostring) result(string_ptr) bind(c, name='getPTStringFromMinutes')
      import :: c_ptr, c_int64_t, c_char
      type(c_ptr) :: string_ptr
      integer(c_int64_t),value :: m
      character(c_char), dimension(*) :: tostring
    end function my_getptstringfromminutes
    !
    function my_getptstringfromhours(h, tostring) result(string_ptr) bind(c, name='getPTStringFromHours')
      import :: c_ptr, c_int64_t, c_char
      type(c_ptr) :: string_ptr
      integer(c_int64_t),value :: h
      character(c_char), dimension(*) :: tostring
    end function my_getptstringfromhours
    !
    function my_timedeltatojuliandelta(td,dt,jd) result(c_pointer) bind(c,name='timeDeltaToJulianDelta')
      import :: c_int64_t, c_ptr
      type(c_ptr) :: c_pointer
      type(c_ptr), value :: td
      type(c_ptr), value :: dt
      type(c_ptr), value :: jd
    end function my_timedeltatojuliandelta
    !
    function my_divideTimeDeltaInSeconds(dividend, divisor, quotient) result(ret_quotient) bind(c,name='divideTimeDeltaInSeconds')
      import :: c_ptr
      type(c_ptr) :: ret_quotient
      type(c_ptr), value :: dividend
      type(c_ptr), value :: divisor
      type(c_ptr), value :: quotient
    end function my_divideTimeDeltaInSeconds
    !
    function my_divideDatetimeDifferenceInSeconds(dt1, dt2, divisor, quotient) result(ret_quotient) &
         &                                                                     bind(c,name='divideDatetimeDifferenceInSeconds')
      import :: c_ptr
      type(c_ptr) :: ret_quotient
      type(c_ptr), value :: dt1
      type(c_ptr), value :: dt2
      type(c_ptr), value :: divisor
      type(c_ptr), value :: quotient
    end function my_divideDatetimeDifferenceInSeconds
    !
    function my_divideTwoDatetimeDiffsInSeconds(dt1_dividend, dt2_dividend, &
         &                                      dt1_divisor,  dt2_divisor, denominator, quotient) &
         &                                       result(ret_quotient) &
         &                                       bind(c,name='divideTwoDatetimeDiffsInSeconds')
      import :: c_ptr
      type(c_ptr) :: ret_quotient
      type(c_ptr), value :: dt1_dividend
      type(c_ptr), value :: dt2_dividend
      type(c_ptr), value :: dt1_divisor
      type(c_ptr), value :: dt2_divisor
      type(c_ptr), value :: denominator
      type(c_ptr), value :: quotient
    end function my_divideTwoDatetimeDiffsInSeconds
    !
  end interface
  !
  interface
    !
    function my_newevent(name, referenceDate, firstdate, lastDate, interval, offset) result(c_pointer) &
         &              bind(c, name='newEvent')
      import :: c_char, c_ptr
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: name
      character(c_char), dimension(*) :: referenceDate
      character(c_char), dimension(*) :: firstDate
      character(c_char), dimension(*) :: lastDate
      character(c_char), dimension(*) :: interval
      character(c_char), dimension(*) :: offset
    end function my_newevent
    !
    function my_neweventwithdatatypes(name, referenceDate, firstdate, lastDate, interval, offset) &
         &                           result(c_pointer) bind(c, name='newEventWithDataType')
      import :: c_char, c_ptr
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: name
      type(c_ptr), value :: referenceDate
      type(c_ptr), value :: firstDate
      type(c_ptr), value :: lastDate
      type(c_ptr), value :: interval
      type(c_ptr), value :: offset
    end function my_neweventwithdatatypes
    !
    subroutine my_deallocateevent(ev) bind(c,name='deallocateEvent')
      import :: c_ptr
      type(c_ptr), value :: ev
    end subroutine my_deallocateevent
    !
    function my_constructandcopyevent(my_event) result(event_copy) bind(c, name='constructAndCopyEvent')
      import :: c_ptr
      type(c_ptr) :: event_copy
      type(c_ptr), value :: my_event
    end function my_constructandcopyevent
    !
    function my_eventtostring(my_event, string) result(string_ptr) bind(c, name='eventToString')
      import :: c_ptr, c_char
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_event
      character(c_char), dimension(*) :: string
    end function my_eventtostring
    !
    function my_isCurrentEventActive(my_event, my_datetime, plus_slack, minus_slack) &
         &                      result(ret) bind(c, name='isCurrentEventActive')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      type(c_ptr), value :: my_datetime
      type(c_ptr), value :: plus_slack
      type(c_ptr), value :: minus_slack
      logical(c_bool) :: ret
    end function my_isCurrentEventActive
    !
    function my_iseventnextinnextday(my_event) result(ret) bind(c, name='iseventNextInNextDay')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      logical(c_bool) :: ret
    end function my_iseventnextinnextday
    !
    function my_iseventnextinnextmonth(my_event) result(ret) bind(c, name='iseventNextInNextMonth')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      logical(c_bool) :: ret
    end function my_iseventnextinnextmonth
    !
    function my_iseventnextinnextyear(my_event) result(ret) bind(c, name='iseventNextInNextYear')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      logical(c_bool) :: ret
    end function my_iseventnextinnextyear
    !
    function my_gettriggernexteventatdatetime(my_event,my_currentdatetime,my_datetime) result(c_pointer) &
         & bind(c, name='getTriggerNextEventAtDateTime')
      import :: c_ptr
      type(c_ptr) :: c_pointer
      type(c_ptr), value :: my_event
      type(c_ptr), value :: my_currentdatetime
      type(c_ptr), value :: my_datetime
    end function my_gettriggernexteventatdatetime
    !
    function my_gettriggeredpreviouseventatdatetime(my_event,my_datetime) result(c_pointer) &
         & bind(c, name='getTriggeredPreviousEventAtDateTime')
      import :: c_ptr
      type(c_ptr) :: c_pointer
      type(c_ptr), value :: my_event
      type(c_ptr), value :: my_datetime
    end function my_gettriggeredpreviouseventatdatetime
    !
    function my_geteventname(my_event, string) result(c_pointer) bind(c, name='getEventName')
      import :: c_ptr, c_char
      type(c_ptr) :: c_pointer
      type(c_ptr), value :: my_event
      character(c_char), dimension(*) :: string
    end function my_geteventname
    !
    function my_getnexteventisfirst(my_event) result(ret) bind(c, name='getNextEventIsFirst')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      logical(c_bool) :: ret
    end function my_getnexteventisfirst
    !
    function my_geteventisfirstinday(my_event) result(ret) bind(c, name='getEventisFirstInDay')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      logical(c_bool) :: ret
    end function my_geteventisfirstinday
    !
    function my_geteventisfirstinmonth(my_event) result(ret) bind(c, name='getEventisFirstInMonth')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      logical(c_bool) :: ret
    end function my_geteventisfirstinmonth
    !
    function my_geteventisfirstinyear(my_event) result(ret) bind(c, name='getEventisFirstInYear')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      logical(c_bool) :: ret
    end function my_geteventisfirstinyear
    !
    function my_geteventislastinday(my_event) result(ret) bind(c, name='getEventisLastInDay')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      logical(c_bool) :: ret
    end function my_geteventislastinday
    !
    function my_geteventislastinmonth(my_event) result(ret) bind(c, name='getEventisLastInMonth')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      logical(c_bool) :: ret
    end function my_geteventislastinmonth
    !
    function my_geteventislastinyear(my_event) result(ret) bind(c, name='getEventisLastInYear')
      import :: c_bool, c_ptr
      type(c_ptr), value :: my_event
      logical(c_bool) :: ret
    end function my_geteventislastinyear
    !
  end interface
  !
  interface
    !
    function my_neweventgroup(name) result(c_pointer) bind(c, name='newEventGroup')
      import :: c_char, c_ptr
      type(c_ptr) :: c_pointer
      character(c_char), dimension(*) :: name
    end function my_neweventgroup
    !
    subroutine my_deallocateeventgroup(evgrp) bind(c,name='deallocateEventGroup')
      import :: c_ptr
      type(c_ptr), value :: evgrp
    end subroutine my_deallocateeventgroup
    !
    function my_addeventtoeventgroup(my_event, my_eventgroup) result(ret) bind(c, name='addNewEventToEventGroup')
      import :: c_bool, c_ptr
      logical(c_bool) :: ret
      type(c_ptr), value :: my_event
      type(c_ptr), value :: my_eventgroup
    end function my_addeventtoeventgroup
    !
    function my_removeeventfromeventgroup(evname, evgrp) result(ret) bind(c, name='removeEventFromEventGroup')
      import :: c_bool, c_char, c_ptr
      logical(c_bool) :: ret
      character(c_char), dimension(*) :: evname
      type(c_ptr), value :: evgrp
    end function my_removeeventfromeventgroup
    !
    function my_geteventgroupname(my_eventgroup, string) result(string_ptr) bind(c, name='getEventGroupName')
      import :: c_ptr, c_char
      type(c_ptr) :: string_ptr
      type(c_ptr), value :: my_eventgroup
      character(c_char), dimension(*) :: string
    end function my_geteventgroupname
    !
  end interface
  !
  interface
    !
    function my_getRepetitions(repetitionString) bind(c, name='getRepetitions')
      import :: c_int, c_char
      integer(c_int) :: my_getRepetitions
      character(c_char), dimension(*) :: repetitionString
    end function my_getRepetitions
    !
    subroutine my_splitRepetitionString(recurringTimeInterval, repetitor, start, end, duration) &
         bind(c, name='splitRepetitionString')
      import :: c_char
      character(c_char), dimension(*) :: recurringTimeInterval
      character(c_char), dimension(*) :: repetitor
      character(c_char), dimension(*) :: start
      character(c_char), dimension(*) :: end
      character(c_char), dimension(*) :: duration
    end subroutine my_splitRepetitionString
    !
  end interface
  !


  interface handle_errno
    module procedure handle_errno_base
    module procedure handle_errno_cond
  end interface handle_errno


contains

  !___________________________________________________________________________
  ! auxiliary routine: handle error code.
  subroutine handle_errno_base(errno, routine_str, lineno)
    integer, intent(IN) :: errno
    integer, intent(IN) :: lineno
    character(LEN=*), intent(IN) :: routine_str
    character(len=max_mtime_error_str_len)     :: error_str
    if (errno /= no_error) then
      call mtime_strerror(errno, error_str)
      write (error_str,'(a,a,i0)') trim(error_str), " :: line ", lineno
      call finish_mtime(routine_str, error_str)
    end if
  end subroutine handle_errno_base


  !___________________________________________________________________________
  ! auxiliary routine: handle error code.
  subroutine handle_errno_cond(lcond, errno, routine_str, lineno)
    logical, intent(IN) :: lcond
    integer, intent(IN) :: errno
    integer, intent(IN) :: lineno

    character(LEN=*), intent(IN) :: routine_str
    character(len=max_mtime_error_str_len)     :: error_str
    if (lcond)  call handle_errno_base(errno, routine_str, lineno)
  end subroutine handle_errno_cond


end module mtime_c_bindings
