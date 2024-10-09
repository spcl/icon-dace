!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
!>
!! @brief Register a callback function and print mtime objects by this function
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________

MODULE mtime_print_by_callback
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  ! PUBLIC :: register_print_mtime_procedure
  PUBLIC :: register_finish_mtime_procedure
  PUBLIC :: print_mtime
  PUBLIC :: finish_mtime
  !
  ! ABSTRACT INTERFACE
  !   SUBROUTINE registered_print_mtime_procedure(leading_text, message_text)
  !     CHARACTER(len=*), INTENT(in) :: leading_text
  !     CHARACTER(len=*), INTENT(in) :: message_text
  !   END SUBROUTINE registered_print_mtime_procedure
  ! END INTERFACE
  !
  !> @cond DOXYGEN_IGNORE_THIS
  ! PROCEDURE(registered_print_mtime_procedure), POINTER :: print_message => NULL()
  !> @endcond DOXYGEN_IGNORE_THIS
  !
  INTERFACE print_mtime
    MODULE PROCEDURE finish_mtime_plain
!    module procedure print_mtime_datetime
!    module procedure print_mtime_timedelta
  END INTERFACE print_mtime
  !
  ABSTRACT INTERFACE
    SUBROUTINE registered_finish_mtime_procedure(leading_text, message_text)
      CHARACTER(len=*), INTENT(in) :: leading_text
      CHARACTER(len=*), INTENT(in) :: message_text
    END SUBROUTINE registered_finish_mtime_procedure
  END INTERFACE
  !
  !> @cond DOXYGEN_IGNORE_THIS
  PROCEDURE(registered_finish_mtime_procedure), POINTER :: finish_message => NULL()
  !> @endcond DOXYGEN_IGNORE_THIS
  !
  INTERFACE finish_mtime
    MODULE PROCEDURE finish_mtime_plain
  END INTERFACE finish_mtime
  !
CONTAINS
  !
!  !>
!  !! @brief Adds a callback function (subroutine) for printing a string
!  !!        representation of a datetime or timedelta. The callback function
!  !!        is required to have two string arguments, a leading text, and
!  !!        a message text.
!  !!
!  !! @param[in] message_procedure a subroutine with two character(len=*) arguments.
!  !!
!  SUBROUTINE register_print_mtime_procedure(message_procedure)
!    PROCEDURE(registered_print_mtime_procedure) :: message_procedure
!    print_message => message_procedure
!  END SUBROUTINE register_print_mtime_procedure
!  !>
!  !! @brief Print a datetime with associated text information by the provided callback function.
!  !!        Can be used via the generic print_mtime subroutine.
!  !!
!  !! @param[in] leading_text the leading information, eg. the caller
!  !! @param[in] message_text the message text provided
!  !! @param[in] this_datetime the mtime datetime to print
!  !!
!  subroutine print_mtime_datetime(leading_text, message_text, this_datetime)
!    character(len=*), intent(in) :: leading_text
!    character(len=*), intent(in) :: message_text
!    type(datetime), pointer :: this_datetime
!    !
!    character(len=max_datetime_str_len)  :: dstring
!    !
!    call datetimeToString(this_datetime, dstring)
!    call print_message(trim(leading_text), trim(message_text)//' '//trim(dstring))
!    !
!  end subroutine print_mtime_datetime
!  !>
!  !! @brief Print a timedelta with associated text information by the provided callback function.
!  !!        Can be used via the generic print_mtime subroutine.
!  !!
!  !! @param[in] leading_text the leading information, eg. the caller
!  !! @param[in] message_text the message text provided
!  !! @param[in] this_timedelta the mtime timedelta to print
!  !!        Can be used via the generic print_mtime subroutine.
!  !!
!  subroutine print_mtime_timedelta(leading_text, message_text, this_timedelta)
!    character(len=*), intent(in) :: leading_text
!    character(len=*), intent(in) :: message_text
!    type(timedelta), pointer :: this_timedelta
!    !
!    character(len=max_timedelta_str_len)  :: tdstring
!    !
!    call timedeltaToString(this_timedelta, tdstring)
!    call print_message(trim(leading_text), trim(message_text)//' '//trim(tdstring))
!    !
!  end subroutine print_mtime_timedelta
  !
  !>
  !! @brief Adds a callback function (subroutine) for printing a
  !!        string representation of a datetime or timedelta and
  !!        finish the program. The callback function is required to
  !!        have two string arguments, a leading text, and a message
  !!        text.
  !!
  !! @param[in] message_procedure a subroutine with two character(len=*) arguments.
  !!
  SUBROUTINE register_finish_mtime_procedure(finish_procedure)
    PROCEDURE(registered_finish_mtime_procedure) :: finish_procedure
    finish_message => finish_procedure
  END SUBROUTINE register_finish_mtime_procedure
  !>
  !! @brief Calling this procedure make mtime issuing an error message and finish.
  !!
  !! @param[in] leading_text the leading information, eg. the caller
  !! @param[in] message_text the message text provided
  !!
  SUBROUTINE finish_mtime_plain(leading_text, message_text)
    CHARACTER(len=*), INTENT(in) :: leading_text
    CHARACTER(len=*), INTENT(in) :: message_text
    IF (ASSOCIATED(finish_message)) THEN
      CALL finish_message(TRIM(leading_text), TRIM(message_text))
    ELSE
      WRITE (0, *) TRIM(leading_text), ": ", TRIM(message_text)
      STOP
    END IF
  END SUBROUTINE finish_mtime_plain
!  !>
!  !! @brief Print a datetime with associated text information by the provided callback function and finish program.
!  !!        Can be used via the generic finish_mtime subroutine.
!  !!
!  !! @param[in] leading_text the leading information, eg. the caller
!  !! @param[in] message_text the message text provided
!  !! @param[in] this_datetime the mtime datetime to print
!  !!
!  subroutine finish_mtime_datetime(leading_text, message_text, this_datetime)
!    character(len=*), intent(in) :: leading_text
!    character(len=*), intent(in) :: message_text
!    type(datetime), pointer :: this_datetime
!    !
!    character(len=max_datetime_str_len)  :: dstring
!    !
!    call datetimeToString(this_datetime, dstring)
!    call finish_message(trim(leading_text), trim(message_text)//' '//trim(dstring))
!    !
!  end subroutine finish_mtime_datetime
!  !>
!  !! @brief Print a timedelta with associated text information by the provided callback function and finish program.
!  !!        Can be used via the generic finish_mtime subroutine.
!  !!
!  !! @param[in] leading_text the leading information, eg. the caller
!  !! @param[in] message_text the message text provided
!  !! @param[in] this_timedelta the mtime timedelta to print
!  !!        Can be used via the generic print_mtime subroutine.
!  !!
!  subroutine finish_mtime_timedelta(leading_text, message_text, this_timedelta)
!    character(len=*), intent(in) :: leading_text
!    character(len=*), intent(in) :: message_text
!    type(timedelta), pointer :: this_timedelta
!    !
!    character(len=max_timedelta_str_len)  :: tdstring
!    !
!    call timedeltaToString(this_timedelta, tdstring)
!    call finish_message(trim(leading_text), trim(message_text)//' '//trim(tdstring))
!    !
!  end subroutine finish_mtime_timedelta
  !
END MODULE mtime_print_by_callback

MODULE mtime_error_handling

  USE mtime_print_by_callback

  IMPLICIT NONE

  PUBLIC

  !
  !> @cond DOXYGEN_IGNORE_THIS
  INTEGER, PARAMETER, PUBLIC :: no_error = 0

  INTEGER, PARAMETER, PUBLIC :: &
       &                calendar_calendartostring = 0*100 + 1

  INTEGER, PARAMETER, PUBLIC :: &
       &                general_arithmetic_error = 0*100 + 2

  INTEGER, PARAMETER, PUBLIC :: &
       &                julianday_newjulianday = 1*100 + 1, &
       &                julianday_juliandaytostring = 1*100 + 3

  INTEGER, PARAMETER, PUBLIC :: &
       &                date_newdatefromstring = 2*100 + 1, &
       &                date_newdatefromraw_yi8 = 2*100 + 2, &
       &                date_newdatefromraw = 2*100 + 3, &
       &                date_newdatefromconstructandcopy = 2*100 + 4, &
       &                date_replacedate = 2*100 + 6, &
       &                date_datetostring = 2*100 + 7, &
       &                date_datetoposixstring = 2*100 + 8

  INTEGER, PARAMETER, PUBLIC :: &
       &                time_newtimefromstring = 3*100 + 1, &
       &                time_newtimefromraw = 3*100 + 2, &
       &                time_newtimefromconstructandcopy = 3*100 + 3, &
       &                time_replacetime = 3*100 + 5, &
       &                time_timetostring = 3*100 + 6, &
       &                time_timetoposixstring = 3*100 + 7

  INTEGER, PARAMETER, PUBLIC :: &
       &                datetime_newdatetimefromstring = 4*100 + 1, &
       &                datetime_newdatetimefromraw_yi8 = 4*100 + 2, &
       &                datetime_newdatetimefromraw = 4*100 + 3, &
       &                datetime_newdatetimefromconstructandcopy = 4*100 + 4, &
       &                datetime_datetimetostring = 4*100 + 6, &
       &                datetime_datetimetoposixstring = 4*100 + 7, &
       &                datetime_getnoofdaysinmonthdatetime = 4*100 + 15, &
       &                datetime_getnoofdaysinyeardatetime = 4*100 + 16, &
       &                datetime_getdayofyearfromdatetime = 4*100 + 17, &
       &                datetime_getnoofsecondselapsedinmonthdatetime = 4*100 + 18, &
       &                datetime_getnoofsecondselapsedindaydatetime = 4*100 + 19, &
       &                datetime_getjuliandayfromdatetime = 4*100 + 20, &
       &                datetime_getdatetimefromjulianday = 4*100 + 21

  INTEGER, PARAMETER, PUBLIC :: &
       &                timedelta_newtimedeltafromstring = 5*100 + 1, &
       &                timedelta_newtimedeltafromraw = 5*100 + 2, &
       &                timedelta_newtimedeltafromraw_yi8 = 5*100 + 3, &
       &                timedelta_newtimedeltafromconstructandcopy = 5*100 + 4, &
       &                timedelta_gettotalmillisecondstimedelta = 5*100 + 8, &
       &                timedelta_gettotalsecondstimedelta = 5*100 + 9, &
       &                timedelta_timedeltatostring = 5*100 + 10, &
       &                timedelta_getptstring = 5*100 + 11

  INTEGER, PARAMETER, PUBLIC :: &
       &                events_newevent = 6*100 + 1, &
       &                events_eventtostring = 6*100 + 3, &
       &                events_gettriggernexteventatdatetime = 6*100 + 8, &
       &                events_gettriggeredpreviouseventatdatetime = 6*100 + 9, &
       &                events_geteventname = 6*100 + 11, &
       &                events_getFirstDatetime = 6*100 + 12, &
       &                events_geteventinterval = 6*100 + 13, &
       &                events_getlastdatetime = 6*100 + 14

  INTEGER, PARAMETER, PUBLIC :: &
       &                eventgroups_neweventgroup = 7*100 + 1, &
       &                eventgroups_geteventgroupname = 7*100 + 6
  !> @endcond DOXYGEN_IGNORE_THIS
  !

CONTAINS
  !>
  !! @brief returns error emssage associated to error number of mtime
  !!
  !! @param[in]   errno           error message number
  !!
  !! @param[out]  error_message   associated error message
  SUBROUTINE mtime_strerror(errno, error_message)
    INTEGER, INTENT(in) :: errno
    CHARACTER(len=*), INTENT(out) :: error_message

    SELECT CASE (errno)
    CASE (no_error)
      error_message = 'no error'

    CASE (general_arithmetic_error)
      error_message = 'error in arithmetic operation'

    CASE (calendar_calendartostring)
      error_message = 'could not retrieve the string in <calendartostring>'

    CASE (julianday_newjulianday)
      error_message = 'could not allocate julianday in <newjulianday>'
    CASE (julianday_juliandaytostring)
      error_message = 'could not retrieve the string in <juliandaytostring>'

    CASE (date_newdatefromstring)
      error_message = 'could not allocate date in <newdate>'
    CASE (date_newdatefromraw_yi8)
      error_message = 'could not allocate date in <newdate>'
    CASE (date_newdatefromraw)
      error_message = 'could not allocate date in <newdate>'
    CASE (date_newdatefromconstructandcopy)
      error_message = 'could not allocate date in <newdate>'
    CASE (date_replacedate)
      error_message = 'replace-date failed in <replacedate>'
    CASE (date_datetostring)
      error_message = 'could not retrieve the string in <datetostring>'
    CASE (date_datetoposixstring)
      error_message = 'could not retrieve the string in <datetoposixstring>'

    CASE (time_newtimefromstring)
      error_message = 'could not allocate time in <newtime>'
    CASE (time_newtimefromraw)
      error_message = 'could not allocate time in <newtime>'
    CASE (time_newtimefromconstructandcopy)
      error_message = 'could not allocate time in <newtime>'
    CASE (time_replacetime)
      error_message = 'replace-time failed in <replacetime>'
    CASE (time_timetostring)
      error_message = 'could not retrieve the string in <timetostring>'
    CASE (time_timetoposixstring)
      error_message = 'could not retrieve the string in <timetoposixstring>'

    CASE (datetime_newdatetimefromstring)
      error_message = 'could not allocate datetime in <newdatetime>'
    CASE (datetime_newdatetimefromraw_yi8)
      error_message = 'could not allocate datetime in <newdatetime>'
    CASE (datetime_newdatetimefromraw)
      error_message = 'could not allocate datetime in <newdatetime>'
    CASE (datetime_newdatetimefromconstructandcopy)
      error_message = 'could not allocate datetime in <newdatetime>'
    CASE (datetime_datetimetostring)
      error_message = 'could not retrieve the string in <timetostring>'
    CASE (datetime_datetimetoposixstring)
      error_message = 'could not retrieve the string in <timetoposixstring>'
    CASE (datetime_getnoofdaysinmonthdatetime)
      error_message = 'could not retrieve the no-of-days-in-month-datetime in <getnoofdaysinmonthdatetime>'
    CASE (datetime_getnoofdaysinyeardatetime)
      error_message = 'could not retrieve the no-of-days-in-year-datetime in <getnoofdaysinyeardatetime>'
    CASE (datetime_getdayofyearfromdatetime)
      error_message = 'could not calculate the day-of-year-from-date-time in <getdayofyearfromdatetime>'
    CASE (datetime_getnoofsecondselapsedinmonthdatetime)
      error_message = 'could not calculate the no-of-seconds-elapsed-in-month-datetime in <getnoofsecondselapsedinmonthdatetime>'
    CASE (datetime_getnoofsecondselapsedindaydatetime)
      error_message = 'could not calculate the no-of-seconds-elapsed-in-day-datetime in <getnoofsecondselapsedindaydatetime>'
    CASE (datetime_getjuliandayfromdatetime)
      error_message = 'could not retrieve the get-julian-day-from-datetime in <getjuliandayfromdatetime>'
    CASE (datetime_getdatetimefromjulianday)
      error_message = 'could not retrieve the get-datetime-from-julian-day in <getdatetimefromjulianday>'

    CASE (timedelta_newtimedeltafromstring)
      error_message = 'could not allocate timedelta in <newtimedelta>'
    CASE (timedelta_newtimedeltafromraw)
      error_message = 'could not allocate timedelta in <newtimedelta>'
    CASE (timedelta_newtimedeltafromraw_yi8)
      error_message = 'could not allocate timedelta in <newtimedelta>'
    CASE (timedelta_newtimedeltafromconstructandcopy)
      error_message = 'could not allocate timedelta in <newtimedelta>'
    CASE (timedelta_gettotalmillisecondstimedelta)
      error_message = 'error in calculating total-milli-seconds-timedelta in <gettotalmillisecondstimedelta>'
    CASE (timedelta_gettotalsecondstimedelta)
      error_message = 'error in calculating total-seconds-in-timedelta in <gettotalmillisecondstimedelta>'
    CASE (timedelta_timedeltatostring)
      error_message = 'could not retrieve the string in <timedeltatostring>'
    CASE (timedelta_getptstring)
      error_message = 'could not retrieve the string in <getptstring*>'

    CASE (events_newevent)
      error_message = 'could not allocate event in <newevent>'
    CASE (events_eventtostring)
      error_message = 'could not retrieve the string in <eventtostring>'
    CASE (events_gettriggernexteventatdatetime)
      error_message = 'failed to retrieve trigger-next-event-at-datetime in <gettriggernexteventatdatetime>'
    CASE (events_gettriggeredpreviouseventatdatetime)
      error_message = 'failed to retrieve triggered-previous-event-at-date-time in <gettriggeredpreviouseventatdatetime>'
    CASE (events_geteventname)
      error_message = 'failed to retrieve event-name in <geteventname>'
    CASE (events_getFirstDatetime)
      error_message = 'could not get the event first datetime'
    CASE (events_geteventinterval)
      error_message = 'could not get the event interval'
    CASE (events_getlastdatetime)
      error_message = 'could not retrieve the event last date'

    CASE (eventgroups_neweventgroup)
      error_message = 'could not allocate eventgroup in <neweventgroup>'
    CASE (eventgroups_geteventgroupname)
      error_message = 'could not retrieve the group-name in <geteventgroupname>'
    END SELECT
  END SUBROUTINE mtime_strerror
  !

END MODULE mtime_error_handling

