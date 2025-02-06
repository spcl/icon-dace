!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
MODULE mtime_c_bindings
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_bool, c_char, c_double, c_float, c_int, c_int64_t, c_ptr
  USE mtime_constants, ONLY: max_mtime_error_str_len
  USE mtime_error_handling, ONLY: finish_mtime, mtime_strerror, no_error
  !
  IMPLICIT NONE
  !
  PUBLIC
  !
  ! Type, bind(c)
  !
  TYPE, BIND(c) :: event
    ! FIXME (some day in the future): This derived type should not specify both -
    ! the linked list element and the element itself.
    TYPE(c_ptr) :: nextEventInGroup
    INTEGER(c_int64_t) :: eventId
    TYPE(c_ptr) :: eventName
    TYPE(c_ptr) :: eventsLastEvaluationDateTime
    TYPE(c_ptr) :: eventReferenceDatetime
    TYPE(c_ptr) :: eventFirstDatetime
    TYPE(c_ptr) :: eventLastDatetime
    TYPE(c_ptr) :: eventInterval
    TYPE(c_ptr) :: eventOffset
    LOGICAL(c_bool) :: neverTriggerEvent
    LOGICAL(c_bool) :: triggerCurrentEvent
    LOGICAL(c_bool) :: nextEventIsFirst
    LOGICAL(c_bool) :: lastEventWasFinal
    LOGICAL(c_bool) :: eventisFirstInDay
    LOGICAL(c_bool) :: eventisFirstInMonth
    LOGICAL(c_bool) :: eventisFirstInYear
    LOGICAL(c_bool) :: eventisLastInDay
    LOGICAL(c_bool) :: eventisLastInMonth
    LOGICAL(c_bool) :: eventisLastInYear
    TYPE(c_ptr) :: triggerNextEventDateTime
    TYPE(c_ptr) :: triggeredPreviousEventDateTime
  END TYPE event
  !
  TYPE, BIND(c) :: eventgroup
    INTEGER(c_int64_t) :: eventGroupId
    TYPE(c_ptr) :: eventGroupName
    TYPE(c_ptr) :: firstEventInGroup
  END TYPE eventgroup
  !
  TYPE, BIND(c) :: julianday
    INTEGER(c_int64_t) :: day  !< the actual Julian day
    INTEGER(c_int64_t) :: ms   !< the milisecond on that particular day
  END TYPE julianday
  !
  TYPE, BIND(c) :: juliandelta
    CHARACTER(c_char)  :: sign
    INTEGER(c_int64_t) :: day
    INTEGER(c_int64_t) :: ms
  END TYPE juliandelta
  !
  TYPE, BIND(c) :: date
    INTEGER(c_int64_t) :: year
    INTEGER(c_int) :: month
    INTEGER(c_int) :: day
  END TYPE date
  !
  TYPE, BIND(c) :: time
    INTEGER(c_int) :: hour
    INTEGER(c_int) :: minute
    INTEGER(c_int) :: second
    INTEGER(c_int) :: ms
  END TYPE time
  !
  TYPE, BIND(c) :: datetime
    TYPE(date) :: date
    TYPE(time) :: time
  END TYPE datetime
  !
  TYPE, BIND(c) :: timedelta
    INTEGER(c_int) :: flag_std_form
    CHARACTER(c_char) :: sign
    INTEGER(c_int64_t) :: year
    INTEGER(c_int) :: month
    INTEGER(c_int) :: day
    INTEGER(c_int) :: hour
    INTEGER(c_int) :: minute
    INTEGER(c_int) :: second
    INTEGER(c_int) :: ms
  END TYPE timedelta
  !
  TYPE, BIND(c) :: divisionquotienttimespan
    INTEGER(c_int64_t) :: quotient; 
    INTEGER(c_int64_t) :: remainder_in_ms; 
  END TYPE divisionquotienttimespan
  !
  ! End Type, bind(c)
  !
  INTERFACE
    !
    SUBROUTINE setCalendar(ct) BIND(c, name='initCalendar') !TESTED-OK
      IMPORT :: c_int
      INTEGER(c_int), VALUE :: ct
    END SUBROUTINE setCalendar
    !
    SUBROUTINE resetCalendar() BIND(c, name='freeCalendar') !TESTED-OK
    END SUBROUTINE resetCalendar
    !
    FUNCTION calendarType() BIND(c, name='getCalendarType') !TESTED-OK
      IMPORT :: c_int
      INTEGER(c_int) :: calendarType
    END FUNCTION calendarType
    !
    FUNCTION my_calendartostring(calendar) RESULT(c_pointer) BIND(c, name='calendarToString') !TESTED-OK
      IMPORT :: c_char, c_ptr
      TYPE(c_ptr) :: c_pointer
      CHARACTER(c_char), DIMENSION(*) :: calendar
    END FUNCTION my_calendartostring
    !
  END INTERFACE

  INTERFACE
    FUNCTION my_newjuliandelta(sign, day, ms) RESULT(c_pointer) BIND(c, name='newJulianDelta')
      IMPORT :: c_int64_t, c_char, c_ptr
      TYPE(c_ptr)               :: c_pointer
      CHARACTER(c_char), VALUE :: sign
      INTEGER(c_int64_t), VALUE :: day
      INTEGER(c_int64_t), VALUE :: ms
    END FUNCTION my_newjuliandelta
    !
    SUBROUTINE my_deallocatejuliandelta(jd) BIND(c, name='deallocateJulianDelta')
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE :: jd
    END SUBROUTINE my_deallocatejuliandelta
  END INTERFACE

  INTERFACE
    FUNCTION my_newjulianday(day, ms) RESULT(c_pointer) BIND(c, name='newJulianDay')
      IMPORT :: c_int64_t, c_ptr
      TYPE(c_ptr) :: c_pointer
      INTEGER(c_int64_t), VALUE :: day
      INTEGER(c_int64_t), VALUE :: ms
    END FUNCTION my_newjulianday
    !
    SUBROUTINE my_deallocatejulianday(jd) BIND(c, name='deallocateJulianDay')
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE :: jd
    END SUBROUTINE my_deallocatejulianday
    !
    FUNCTION my_replacejulianday(src, dest) RESULT(ret_dest) BIND(c, name='replaceJulianday')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: ret_dest
      TYPE(c_ptr), VALUE :: src
      TYPE(c_ptr), VALUE :: dest
    END FUNCTION my_replacejulianday
    !
    PURE FUNCTION my_comparejulianday(op1, op2) RESULT(ret) BIND(c, name='compareJulianDay')
      IMPORT :: c_ptr, c_int
      INTEGER(c_int) :: ret
      TYPE(c_ptr), VALUE, INTENT(in) :: op1
      TYPE(c_ptr), VALUE, INTENT(in) :: op2
    END FUNCTION my_comparejulianday
    !
    FUNCTION my_juliandaytostring(my_julianday, string) RESULT(string_ptr) BIND(c, name='juliandayToString')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      TYPE(c_ptr), VALUE :: my_julianday
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_juliandaytostring
    !
    FUNCTION my_addjuliandeltatojulianday(my_julianday, my_juliandelta, ret_julianday) RESULT(julianday_ptr) &
         &   BIND(c, name='addJulianDeltaToJulianDay')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: julianday_ptr
      TYPE(c_ptr), VALUE :: my_julianday, my_juliandelta, ret_julianday
    END FUNCTION my_addjuliandeltatojulianday
    !
    FUNCTION my_subtractjulianday(my_julianday1, my_julianday2, ret_juliandelta) RESULT(juliandelta_ptr) &
         &   BIND(c, name='subtractJulianDay')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: juliandelta_ptr
      TYPE(c_ptr), VALUE :: my_julianday1, my_julianday2, ret_juliandelta
    END FUNCTION my_subtractjulianday
    !
  END INTERFACE
  !
  INTERFACE
    !
    FUNCTION my_newdatefromstring(string) RESULT(c_pointer) BIND(c, name='newDate')
      IMPORT :: c_char, c_ptr
      TYPE(c_ptr) :: c_pointer
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_newdatefromstring
    !
    FUNCTION my_newrawdate(year, month, day) RESULT(c_pointer) BIND(c, name='newRawDate')
      IMPORT :: c_int64_t, c_int, c_ptr
      TYPE(c_ptr) :: c_pointer
      INTEGER(c_int64_t), VALUE :: year
      INTEGER(c_int), VALUE :: month, day
    END FUNCTION my_newrawdate
    !
    FUNCTION my_constructandcopydate(d) RESULT(c_pointer) BIND(c, name='constructAndCopyDate')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: c_pointer
      TYPE(c_ptr), VALUE :: d
    END FUNCTION my_constructandcopydate
    !
    SUBROUTINE my_deallocatedate(d) BIND(c, name='deallocateDate')
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE :: d
    END SUBROUTINE my_deallocatedate
    !
    FUNCTION my_replacedate(src, dest) RESULT(ret_dest) BIND(c, name='replaceDate')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: ret_dest
      TYPE(c_ptr), VALUE :: src, dest
    END FUNCTION my_replacedate
    !
    FUNCTION my_datetostring(my_date, string) RESULT(string_ptr) BIND(c, name='dateToString')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      TYPE(c_ptr), VALUE :: my_date
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_datetostring
    !
    FUNCTION my_datetoposixstring(my_date, string, fmtstr) RESULT(string_ptr) BIND(c, name='dateToPosixString')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      TYPE(c_ptr), VALUE :: my_date
      CHARACTER(c_char), DIMENSION(*) :: string
      CHARACTER(c_char), DIMENSION(*) :: fmtstr
    END FUNCTION my_datetoposixstring
    !
  END INTERFACE
  !
  INTERFACE
    !
    FUNCTION my_newtimefromstring(string) RESULT(c_pointer) BIND(c, name='newTime')
      IMPORT :: c_char, c_ptr
      TYPE(c_ptr) :: c_pointer
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_newtimefromstring
    !
    FUNCTION my_newrawtime(hour, minute, second, ms) RESULT(c_pointer) BIND(c, name='newRawTime')
      IMPORT :: c_int, c_ptr
      TYPE(c_ptr) :: c_pointer
      INTEGER(c_int), VALUE :: hour, minute, second, ms
    END FUNCTION my_newrawtime
    !
    FUNCTION my_constructandcopytime(t) RESULT(c_pointer) BIND(c, name='constructAndCopyTime')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: c_pointer
      TYPE(c_ptr), VALUE :: t
    END FUNCTION my_constructandcopytime
    !
    SUBROUTINE my_deallocatetime(t) BIND(c, name='deallocateTime')
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE :: t
    END SUBROUTINE my_deallocatetime
    !
    FUNCTION my_replacetime(src, dest) RESULT(ret_dest) BIND(c, name='replaceTime')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: ret_dest
      TYPE(c_ptr), VALUE :: src, dest
    END FUNCTION my_replacetime
    !
    FUNCTION my_timetostring(my_time, string) RESULT(string_ptr) BIND(c, name='timeToString')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      TYPE(c_ptr), VALUE :: my_time
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_timetostring
    !
    FUNCTION my_timetoposixstring(my_time, string, fmtstr) RESULT(string_ptr) BIND(c, name='timeToPosixString')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      TYPE(c_ptr), VALUE :: my_time
      CHARACTER(c_char), DIMENSION(*) :: string
      CHARACTER(c_char), DIMENSION(*) :: fmtstr
    END FUNCTION my_timetoposixstring
    !
  END INTERFACE
  !
  INTERFACE
    !
    FUNCTION my_newdatetime(string) RESULT(c_pointer) BIND(c, name='newDateTime')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: c_pointer
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_newdatetime
    !
    FUNCTION my_newrawdatetime(year, month, day, hour, minute, second, ms) RESULT(c_pointer) &
         &             BIND(c, name='newRawDateTime')
      IMPORT :: c_int64_t, c_int, c_ptr
      TYPE(c_ptr) :: c_pointer
      INTEGER(c_int64_t), VALUE :: year
      INTEGER(c_int), VALUE :: month, day, hour, minute, second, ms
    END FUNCTION my_newrawdatetime
    !
    FUNCTION my_constructandcopydatetime(dt) RESULT(c_pointer) BIND(c, name='constructAndCopyDateTime')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: c_pointer
      TYPE(c_ptr), VALUE :: dt
    END FUNCTION my_constructandcopydatetime
    !
    SUBROUTINE my_deallocatedatetime(dt) BIND(c, name='deallocateDateTime')
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE :: dt
    END SUBROUTINE my_deallocatedatetime
    !
    PURE FUNCTION my_comparedatetime(op1, op2) RESULT(ret) BIND(c, name='compareDatetime')
      IMPORT :: c_ptr, c_int
      INTEGER(c_int) :: ret
      TYPE(c_ptr), VALUE, INTENT(in) :: op1
      TYPE(c_ptr), VALUE, INTENT(in) :: op2
    END FUNCTION my_comparedatetime

    FUNCTION my_replacedatetime(src, dest) RESULT(ret_dest) BIND(c, name='replaceDatetime')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: ret_dest
      TYPE(c_ptr), VALUE :: src
      TYPE(c_ptr), VALUE :: dest
    END FUNCTION my_replacedatetime
    !
    FUNCTION my_datetimetostring(my_time, string) RESULT(string_ptr) BIND(c, name='datetimeToString')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      TYPE(c_ptr), VALUE :: my_time
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_datetimetostring
    !
    FUNCTION my_datetimetoposixstring(my_time, string, fmtstr) RESULT(string_ptr) BIND(c, name='datetimeToPosixString')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      TYPE(c_ptr), VALUE :: my_time
      CHARACTER(c_char), DIMENSION(*) :: string
      CHARACTER(c_char), DIMENSION(*) :: fmtstr
    END FUNCTION my_datetimetoposixstring
    !
    FUNCTION my_getnoofdaysinmonthdatetime(dt) BIND(c, name='getNoOfDaysInMonthDateTime')
      IMPORT :: c_int, c_ptr
      INTEGER(c_int) :: my_getnoofdaysinmonthdatetime
      TYPE(c_ptr), VALUE :: dt
    END FUNCTION my_getnoofdaysinmonthdatetime
    !
    FUNCTION my_getnoofdaysinyeardatetime(dt) BIND(c, name='getNoOfDaysInYearDateTime')
      IMPORT :: c_int, c_ptr
      INTEGER(c_int) :: my_getnoofdaysinyeardatetime
      TYPE(c_ptr), VALUE :: dt
    END FUNCTION my_getnoofdaysinyeardatetime
    !
    FUNCTION my_getdayofyearfromdatetime(dt) BIND(c, name='getDayOfYearFromDateTime')
      IMPORT :: c_int, c_ptr
      INTEGER(c_int) :: my_getdayofyearfromdatetime
      TYPE(c_ptr), VALUE :: dt
    END FUNCTION my_getdayofyearfromdatetime
    !
    FUNCTION my_getnoofsecondselapsedinmonthdatetime(dt) BIND(c, name='getNoOfSecondsElapsedInMonthDateTime')
      IMPORT :: c_int64_t, c_ptr
      INTEGER(c_int64_t) :: my_getnoofsecondselapsedinmonthdatetime
      TYPE(c_ptr), VALUE :: dt
    END FUNCTION my_getnoofsecondselapsedinmonthdatetime
    !
    FUNCTION my_getnoofsecondselapsedindaydatetime(dt) BIND(c, name='getNoOfSecondsElapsedInDayDateTime')
      IMPORT :: c_int, c_ptr
      INTEGER(c_int) :: my_getnoofsecondselapsedindaydatetime
      TYPE(c_ptr), VALUE :: dt
    END FUNCTION my_getnoofsecondselapsedindaydatetime
    !
    FUNCTION my_getjuliandayfromdatetime(dt, jd) BIND(c, name='getJulianDayFromDateTime')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: my_getjuliandayfromdatetime
      TYPE(c_ptr), VALUE :: dt, jd
    END FUNCTION my_getjuliandayfromdatetime
    !
    FUNCTION my_getdatetimefromjulianday(jd, dt) BIND(c, name='getDateTimeFromJulianDay')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: my_getdatetimefromjulianday
      TYPE(c_ptr), VALUE :: jd, dt
    END FUNCTION my_getdatetimefromjulianday
    !
  END INTERFACE
  !
  INTERFACE
    !
    FUNCTION my_newtimedeltafromstring(string) RESULT(c_pointer) BIND(c, name='newTimeDelta')
      IMPORT :: c_char, c_ptr
      TYPE(c_ptr)                     :: c_pointer
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_newtimedeltafromstring
    !
    FUNCTION my_newrawtimedelta(sign, year, month, day, hour, minute, second, ms) RESULT(c_pointer) &
         &             BIND(c, name='newRawTimeDelta')
      IMPORT :: c_int64_t, c_char, c_int, c_ptr
      TYPE(c_ptr) :: c_pointer
      CHARACTER(c_char), VALUE :: sign
      INTEGER(c_int64_t), VALUE :: year
      INTEGER(c_int), VALUE :: month, day, hour, minute, second, ms
    END FUNCTION my_newrawtimedelta
    !
    FUNCTION my_constructandcopytimedelta(td) RESULT(c_pointer) BIND(c, name='constructAndCopyTimeDelta')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: c_pointer
      TYPE(c_ptr), VALUE :: td
    END FUNCTION my_constructandcopytimedelta
    !
    SUBROUTINE my_deallocatetimedelta(dt) BIND(c, name='deallocateTimeDelta')
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE :: dt
    END SUBROUTINE my_deallocatetimedelta
    !
    PURE FUNCTION my_comparetimedelta(op1, op2) RESULT(ret) BIND(c, name='compareTimeDelta')
      IMPORT :: c_ptr, c_int
      INTEGER(c_int) :: ret
      TYPE(c_ptr), VALUE, INTENT(in) :: op1
      TYPE(c_ptr), VALUE, INTENT(in) :: op2
    END FUNCTION my_comparetimedelta
    !
    FUNCTION my_gettimedeltafromdate(my_date1, my_date2, timedelta_return) RESULT(timedelta_ptr) &
         &    BIND(c, name='getTimeDeltaFromDate')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: timedelta_ptr
      TYPE(c_ptr), VALUE :: my_date1, my_date2, timedelta_return
    END FUNCTION my_gettimedeltafromdate
    !
    FUNCTION my_gettimedeltafromdatetime(my_datetime1, my_datetime2, timedelta_return) RESULT(timedelta_ptr) &
         &    BIND(c, name='getTimeDeltaFromDateTime')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: timedelta_ptr
      TYPE(c_ptr), VALUE :: my_datetime1, my_datetime2, timedelta_return
    END FUNCTION my_gettimedeltafromdatetime
    !
    FUNCTION my_gettotalmillisecondstimedelta(my_timedelta, my_datetime) &
         &   BIND(c, name='getTotalMilliSecondsTimeDelta')
      IMPORT :: c_ptr, c_int64_t
      INTEGER(c_int64_t) :: my_gettotalmillisecondstimedelta
      TYPE(c_ptr), VALUE :: my_timedelta, my_datetime
    END FUNCTION my_gettotalmillisecondstimedelta
    !
    FUNCTION my_gettotalsecondstimedelta(my_timedelta, my_datetime) BIND(c, name='getTotalSecondsTimeDelta')
      IMPORT :: c_ptr, c_int64_t
      INTEGER(c_int64_t) :: my_gettotalsecondstimedelta
      TYPE(c_ptr), VALUE :: my_timedelta, my_datetime
    END FUNCTION my_gettotalsecondstimedelta
    !
    FUNCTION my_timedeltatostring(td, tostring) RESULT(string_ptr) BIND(c, name='timedeltaToString')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      TYPE(c_ptr), VALUE :: td
      CHARACTER(c_char), DIMENSION(*) :: tostring
    END FUNCTION my_timedeltatostring
    !
    FUNCTION my_addtimedeltatodatetime(my_datetime, my_timedelta, ret_datetime) RESULT(datetime_ptr) &
         &   BIND(c, name='addTimeDeltaToDateTime')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: datetime_ptr
      TYPE(c_ptr), VALUE :: my_datetime, my_timedelta, ret_datetime
    END FUNCTION my_addtimedeltatodatetime
    !
    FUNCTION my_addtimedeltatodate(my_date, my_timedelta, ret_date) RESULT(date_ptr) &
         &   BIND(c, name='addTimeDeltaToDate')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: date_ptr
      TYPE(c_ptr), VALUE :: my_date, my_timedelta, ret_date
    END FUNCTION my_addtimedeltatodate
    !
    FUNCTION my_elementwisescalarmultiplytimedelta(my_timedelta, lambda, scaled_timedelta) RESULT(timedelta_return) &
         &   BIND(c, name='elementwiseScalarMultiplyTimeDelta')
      IMPORT :: c_ptr, c_int64_t
      TYPE(c_ptr) ::timedelta_return
      TYPE(c_ptr), VALUE :: my_timedelta, scaled_timedelta
      INTEGER(c_int64_t), VALUE :: lambda
    END FUNCTION my_elementwisescalarmultiplytimedelta
    !
    FUNCTION my_elementwisescalarmultiplytimedeltadp(my_timedelta, lambda, scaled_timedelta) RESULT(timedelta_return) &
         &   BIND(c, name='elementwiseScalarMultiplyTimeDeltaDP')
      IMPORT :: c_ptr, c_double
      TYPE(c_ptr) ::timedelta_return
      TYPE(c_ptr), VALUE :: my_timedelta, scaled_timedelta
      REAL(c_double), VALUE :: lambda
    END FUNCTION my_elementwisescalarmultiplytimedeltadp
    !
    FUNCTION my_elementwiseAddTimeDeltatoTimeDelta(my_timedelta1, my_timedelta2, added_timedelta) RESULT(timedelta_return) &
         &   BIND(c, name='elementwiseAddTimeDeltatoTimeDelta')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: timedelta_return
      TYPE(c_ptr), VALUE :: my_timedelta1, my_timedelta2, added_timedelta
    END FUNCTION my_elementwiseAddTimeDeltatoTimeDelta
    !
    FUNCTION my_modulotimedelta(a, p, quot) RESULT(rem) BIND(c, name='moduloTimedelta')
      IMPORT :: c_ptr, c_int64_t
      INTEGER(c_int64_t) :: rem
      TYPE(c_ptr), VALUE :: a
      TYPE(c_ptr), VALUE :: p
      TYPE(c_ptr), VALUE :: quot
    END FUNCTION my_modulotimedelta
    !
    FUNCTION my_getptstringfromms(ms, tostring) RESULT(string_ptr) BIND(c, name='getPTStringFromMS')
      IMPORT :: c_ptr, c_int64_t, c_char
      TYPE(c_ptr) :: string_ptr
      INTEGER(c_int64_t), VALUE :: ms
      CHARACTER(c_char), DIMENSION(*) :: tostring
    END FUNCTION my_getptstringfromms
    !
    FUNCTION my_getptstringfromsecondsint(s, tostring) RESULT(string_ptr) BIND(c, name='getPTStringFromSeconds')
      IMPORT :: c_ptr, c_int64_t, c_char
      TYPE(c_ptr) :: string_ptr
      INTEGER(c_int64_t), VALUE :: s
      CHARACTER(c_char), DIMENSION(*) :: tostring
    END FUNCTION my_getptstringfromsecondsint
    !
    FUNCTION my_getptstringfromsecondsfloat(s, tostring) RESULT(string_ptr) BIND(c, name='getPTStringFromSecondsFloat')
      IMPORT :: c_ptr, c_float, c_char
      TYPE(c_ptr) :: string_ptr
      REAL(c_float), VALUE :: s
      CHARACTER(c_char), DIMENSION(*) :: tostring
    END FUNCTION my_getptstringfromsecondsfloat
    !
    FUNCTION my_getptstringfromsecondsdouble(s, tostring) RESULT(string_ptr) BIND(c, name='getPTStringFromSecondsDouble')
      IMPORT :: c_ptr, c_double, c_char
      TYPE(c_ptr) :: string_ptr
      REAL(c_double), VALUE :: s
      CHARACTER(c_char), DIMENSION(*) :: tostring
    END FUNCTION my_getptstringfromsecondsdouble
    !
    FUNCTION my_getptstringfromminutes(m, tostring) RESULT(string_ptr) BIND(c, name='getPTStringFromMinutes')
      IMPORT :: c_ptr, c_int64_t, c_char
      TYPE(c_ptr) :: string_ptr
      INTEGER(c_int64_t), VALUE :: m
      CHARACTER(c_char), DIMENSION(*) :: tostring
    END FUNCTION my_getptstringfromminutes
    !
    FUNCTION my_getptstringfromhours(h, tostring) RESULT(string_ptr) BIND(c, name='getPTStringFromHours')
      IMPORT :: c_ptr, c_int64_t, c_char
      TYPE(c_ptr) :: string_ptr
      INTEGER(c_int64_t), VALUE :: h
      CHARACTER(c_char), DIMENSION(*) :: tostring
    END FUNCTION my_getptstringfromhours
    !
    FUNCTION my_timedeltatojuliandelta(td, dt, jd) RESULT(c_pointer) BIND(c, name='timeDeltaToJulianDelta')
      IMPORT :: c_int64_t, c_ptr
      TYPE(c_ptr) :: c_pointer
      TYPE(c_ptr), VALUE :: td
      TYPE(c_ptr), VALUE :: dt
      TYPE(c_ptr), VALUE :: jd
    END FUNCTION my_timedeltatojuliandelta
    !
    FUNCTION my_divideTimeDeltaInSeconds(dividend, divisor, quotient) RESULT(ret_quotient) &
         &                                                            BIND(c, name='divideTimeDeltaInSeconds')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: ret_quotient
      TYPE(c_ptr), VALUE :: dividend
      TYPE(c_ptr), VALUE :: divisor
      TYPE(c_ptr), VALUE :: quotient
    END FUNCTION my_divideTimeDeltaInSeconds
    !
    FUNCTION my_divideDatetimeDifferenceInSeconds(dt1, dt2, divisor, quotient) RESULT(ret_quotient) &
         &                                                                     BIND(c, name='divideDatetimeDifferenceInSeconds')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: ret_quotient
      TYPE(c_ptr), VALUE :: dt1
      TYPE(c_ptr), VALUE :: dt2
      TYPE(c_ptr), VALUE :: divisor
      TYPE(c_ptr), VALUE :: quotient
    END FUNCTION my_divideDatetimeDifferenceInSeconds
    !
    FUNCTION my_divideTwoDatetimeDiffsInSeconds(dt1_dividend, dt2_dividend, &
         &                                      dt1_divisor, dt2_divisor, denominator, quotient) &
         &                                       RESULT(ret_quotient) &
         &                                       BIND(c, name='divideTwoDatetimeDiffsInSeconds')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: ret_quotient
      TYPE(c_ptr), VALUE :: dt1_dividend
      TYPE(c_ptr), VALUE :: dt2_dividend
      TYPE(c_ptr), VALUE :: dt1_divisor
      TYPE(c_ptr), VALUE :: dt2_divisor
      TYPE(c_ptr), VALUE :: denominator
      TYPE(c_ptr), VALUE :: quotient
    END FUNCTION my_divideTwoDatetimeDiffsInSeconds
    !
  END INTERFACE
  !
  INTERFACE
    !
    FUNCTION my_newevent(name, referenceDate, firstdate, lastDate, interval, offset) RESULT(c_pointer) &
         &              BIND(c, name='newEvent')
      IMPORT :: c_char, c_ptr
      TYPE(c_ptr) :: c_pointer
      CHARACTER(c_char), DIMENSION(*) :: name
      CHARACTER(c_char), DIMENSION(*) :: referenceDate
      CHARACTER(c_char), DIMENSION(*) :: firstDate
      CHARACTER(c_char), DIMENSION(*) :: lastDate
      CHARACTER(c_char), DIMENSION(*) :: interval
      CHARACTER(c_char), DIMENSION(*) :: offset
    END FUNCTION my_newevent
    !
    FUNCTION my_neweventwithdatatypes(name, referenceDate, firstdate, lastDate, interval, offset) &
         &                           RESULT(c_pointer) BIND(c, name='newEventWithDataType')
      IMPORT :: c_char, c_ptr
      TYPE(c_ptr) :: c_pointer
      CHARACTER(c_char), DIMENSION(*) :: name
      TYPE(c_ptr), VALUE :: referenceDate
      TYPE(c_ptr), VALUE :: firstDate
      TYPE(c_ptr), VALUE :: lastDate
      TYPE(c_ptr), VALUE :: interval
      TYPE(c_ptr), VALUE :: offset
    END FUNCTION my_neweventwithdatatypes
    !
    SUBROUTINE my_deallocateevent(ev) BIND(c, name='deallocateEvent')
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE :: ev
    END SUBROUTINE my_deallocateevent
    !
    FUNCTION my_constructandcopyevent(my_event) RESULT(event_copy) BIND(c, name='constructAndCopyEvent')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: event_copy
      TYPE(c_ptr), VALUE :: my_event
    END FUNCTION my_constructandcopyevent
    !
    FUNCTION my_eventtostring(my_event, string) RESULT(string_ptr) BIND(c, name='eventToString')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      TYPE(c_ptr), VALUE :: my_event
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_eventtostring
    !
    FUNCTION my_isCurrentEventActive(my_event, my_datetime, plus_slack, minus_slack) &
         &                      RESULT(ret) BIND(c, name='isCurrentEventActive')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      TYPE(c_ptr), VALUE :: my_datetime
      TYPE(c_ptr), VALUE :: plus_slack
      TYPE(c_ptr), VALUE :: minus_slack
      LOGICAL(c_bool) :: ret
    END FUNCTION my_isCurrentEventActive
    !
    FUNCTION my_iseventnextinnextday(my_event) RESULT(ret) BIND(c, name='iseventNextInNextDay')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      LOGICAL(c_bool) :: ret
    END FUNCTION my_iseventnextinnextday
    !
    FUNCTION my_iseventnextinnextmonth(my_event) RESULT(ret) BIND(c, name='iseventNextInNextMonth')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      LOGICAL(c_bool) :: ret
    END FUNCTION my_iseventnextinnextmonth
    !
    FUNCTION my_iseventnextinnextyear(my_event) RESULT(ret) BIND(c, name='iseventNextInNextYear')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      LOGICAL(c_bool) :: ret
    END FUNCTION my_iseventnextinnextyear
    !
    FUNCTION my_gettriggernexteventatdatetime(my_event, my_currentdatetime, my_datetime) RESULT(c_pointer) &
         & BIND(c, name='getTriggerNextEventAtDateTime')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: c_pointer
      TYPE(c_ptr), VALUE :: my_event
      TYPE(c_ptr), VALUE :: my_currentdatetime
      TYPE(c_ptr), VALUE :: my_datetime
    END FUNCTION my_gettriggernexteventatdatetime
    !
    FUNCTION my_gettriggeredpreviouseventatdatetime(my_event, my_datetime) RESULT(c_pointer) &
         & BIND(c, name='getTriggeredPreviousEventAtDateTime')
      IMPORT :: c_ptr
      TYPE(c_ptr) :: c_pointer
      TYPE(c_ptr), VALUE :: my_event
      TYPE(c_ptr), VALUE :: my_datetime
    END FUNCTION my_gettriggeredpreviouseventatdatetime
    !
    FUNCTION my_geteventname(my_event, string) RESULT(c_pointer) BIND(c, name='getEventName')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: c_pointer
      TYPE(c_ptr), VALUE :: my_event
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_geteventname
    !
    FUNCTION my_getnexteventisfirst(my_event) RESULT(ret) BIND(c, name='getNextEventIsFirst')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      LOGICAL(c_bool) :: ret
    END FUNCTION my_getnexteventisfirst
    !
    FUNCTION my_geteventisfirstinday(my_event) RESULT(ret) BIND(c, name='getEventisFirstInDay')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      LOGICAL(c_bool) :: ret
    END FUNCTION my_geteventisfirstinday
    !
    FUNCTION my_geteventisfirstinmonth(my_event) RESULT(ret) BIND(c, name='getEventisFirstInMonth')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      LOGICAL(c_bool) :: ret
    END FUNCTION my_geteventisfirstinmonth
    !
    FUNCTION my_geteventisfirstinyear(my_event) RESULT(ret) BIND(c, name='getEventisFirstInYear')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      LOGICAL(c_bool) :: ret
    END FUNCTION my_geteventisfirstinyear
    !
    FUNCTION my_geteventislastinday(my_event) RESULT(ret) BIND(c, name='getEventisLastInDay')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      LOGICAL(c_bool) :: ret
    END FUNCTION my_geteventislastinday
    !
    FUNCTION my_geteventislastinmonth(my_event) RESULT(ret) BIND(c, name='getEventisLastInMonth')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      LOGICAL(c_bool) :: ret
    END FUNCTION my_geteventislastinmonth
    !
    FUNCTION my_geteventislastinyear(my_event) RESULT(ret) BIND(c, name='getEventisLastInYear')
      IMPORT :: c_bool, c_ptr
      TYPE(c_ptr), VALUE :: my_event
      LOGICAL(c_bool) :: ret
    END FUNCTION my_geteventislastinyear
    !
  END INTERFACE
  !
  INTERFACE
    !
    FUNCTION my_neweventgroup(name) RESULT(c_pointer) BIND(c, name='newEventGroup')
      IMPORT :: c_char, c_ptr
      TYPE(c_ptr) :: c_pointer
      CHARACTER(c_char), DIMENSION(*) :: name
    END FUNCTION my_neweventgroup
    !
    SUBROUTINE my_deallocateeventgroup(evgrp) BIND(c, name='deallocateEventGroup')
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE :: evgrp
    END SUBROUTINE my_deallocateeventgroup
    !
    FUNCTION my_addeventtoeventgroup(my_event, my_eventgroup) RESULT(ret) BIND(c, name='addNewEventToEventGroup')
      IMPORT :: c_bool, c_ptr
      LOGICAL(c_bool) :: ret
      TYPE(c_ptr), VALUE :: my_event
      TYPE(c_ptr), VALUE :: my_eventgroup
    END FUNCTION my_addeventtoeventgroup
    !
    FUNCTION my_removeeventfromeventgroup(evname, evgrp) RESULT(ret) BIND(c, name='removeEventFromEventGroup')
      IMPORT :: c_bool, c_char, c_ptr
      LOGICAL(c_bool) :: ret
      CHARACTER(c_char), DIMENSION(*) :: evname
      TYPE(c_ptr), VALUE :: evgrp
    END FUNCTION my_removeeventfromeventgroup
    !
    FUNCTION my_geteventgroupname(my_eventgroup, string) RESULT(string_ptr) BIND(c, name='getEventGroupName')
      IMPORT :: c_ptr, c_char
      TYPE(c_ptr) :: string_ptr
      TYPE(c_ptr), VALUE :: my_eventgroup
      CHARACTER(c_char), DIMENSION(*) :: string
    END FUNCTION my_geteventgroupname
    !
  END INTERFACE
  !
  INTERFACE
    !
    FUNCTION my_getRepetitions(repetitionString) BIND(c, name='getRepetitions')
      IMPORT :: c_int, c_char
      INTEGER(c_int) :: my_getRepetitions
      CHARACTER(c_char), DIMENSION(*) :: repetitionString
    END FUNCTION my_getRepetitions
    !
    SUBROUTINE my_splitRepetitionString(recurringTimeInterval, repetitor, start, END, duration) &
      BIND(c, name='splitRepetitionString')
      IMPORT :: c_char
      CHARACTER(c_char), DIMENSION(*) :: recurringTimeInterval
      CHARACTER(c_char), DIMENSION(*) :: repetitor
      CHARACTER(c_char), DIMENSION(*) :: start
      CHARACTER(c_char), DIMENSION(*) :: END
      CHARACTER(c_char), DIMENSION(*) :: duration
    END SUBROUTINE my_splitRepetitionString
    !
  END INTERFACE
  !
  INTERFACE handle_errno
    MODULE PROCEDURE handle_errno_base
    MODULE PROCEDURE handle_errno_cond
  END INTERFACE handle_errno

CONTAINS

  !___________________________________________________________________________
  ! auxiliary routine: handle error code.
  SUBROUTINE handle_errno_base(errno, routine_str, lineno)
    INTEGER, INTENT(IN) :: errno
    INTEGER, INTENT(IN) :: lineno
    CHARACTER(LEN=*), INTENT(IN) :: routine_str
    CHARACTER(len=max_mtime_error_str_len)     :: error_str
    IF (errno /= no_error) THEN
      CALL mtime_strerror(errno, error_str)
      WRITE (error_str, '(a,a,i0)') TRIM(error_str), " :: line ", lineno
      CALL finish_mtime(routine_str, error_str)
    END IF
  END SUBROUTINE handle_errno_base

  !___________________________________________________________________________
  ! auxiliary routine: handle error code.
  SUBROUTINE handle_errno_cond(lcond, errno, routine_str, lineno)
    LOGICAL, INTENT(IN) :: lcond
    INTEGER, INTENT(IN) :: errno
    INTEGER, INTENT(IN) :: lineno

    CHARACTER(LEN=*), INTENT(IN) :: routine_str
    IF (lcond) CALL handle_errno_base(errno, routine_str, lineno)
  END SUBROUTINE handle_errno_cond

END MODULE mtime_c_bindings
