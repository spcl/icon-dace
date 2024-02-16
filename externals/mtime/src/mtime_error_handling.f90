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

module mtime_print_by_callback
  !
!  use mtime_datetime
!  use mtime_timedelta
  !
  implicit none
  !
  private
  !
  public :: register_print_mtime_procedure
  public :: register_finish_mtime_procedure
  public :: print_mtime
  public :: finish_mtime
  !
  abstract interface
    subroutine registered_print_mtime_procedure(leading_text, message_text)
      character(len=*), intent(in) :: leading_text
      character(len=*), intent(in) :: message_text
    end subroutine registered_print_mtime_procedure
  end interface
  !
  !> @cond DOXYGEN_IGNORE_THIS
  procedure(registered_print_mtime_procedure), pointer :: print_message => null()
  !> @endcond DOXYGEN_IGNORE_THIS
  !
  interface print_mtime
    module procedure finish_mtime_plain
!    module procedure print_mtime_datetime
!    module procedure print_mtime_timedelta
  end interface print_mtime
  !
  abstract interface
    subroutine registered_finish_mtime_procedure(leading_text, message_text)
      character(len=*), intent(in) :: leading_text
      character(len=*), intent(in) :: message_text
    end subroutine registered_finish_mtime_procedure
  end interface
  !
  !> @cond DOXYGEN_IGNORE_THIS
  procedure(registered_finish_mtime_procedure), pointer :: finish_message => null()
  !> @endcond DOXYGEN_IGNORE_THIS
  !
  interface finish_mtime
    module procedure finish_mtime_plain
  end interface finish_mtime
  !
contains
  !
  !>
  !! @brief Adds a callback function (subroutine) for printing a string
  !!        representation of a datetime or timedelta. The callback function
  !!        is required to have two string arguments, a leading text, and
  !!        a message text.
  !!
  !! @param[in] message_procedure a subroutine with two character(len=*) arguments.
  !!
  subroutine register_print_mtime_procedure(message_procedure)
    procedure(registered_print_mtime_procedure) :: message_procedure
    print_message => message_procedure
  end subroutine register_print_mtime_procedure
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
  subroutine register_finish_mtime_procedure(finish_procedure)
    procedure(registered_finish_mtime_procedure) :: finish_procedure
    finish_message => finish_procedure
  end subroutine register_finish_mtime_procedure
  !>
  !! @brief Calling this procedure make mtime issuing an error message and finish.
  !!
  !! @param[in] leading_text the leading information, eg. the caller
  !! @param[in] message_text the message text provided
  !!
  SUBROUTINE finish_mtime_plain(leading_text, message_text)
    character(len=*), intent(in) :: leading_text
    character(len=*), intent(in) :: message_text
    IF (ASSOCIATED(finish_message)) THEN
      CALL finish_message(TRIM(leading_text), TRIM(message_text))
    ELSE
      WRITE (0,*) TRIM(leading_text), ": ", TRIM(message_text)
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
end module mtime_print_by_callback


module mtime_error_handling

  use mtime_constants
  use mtime_print_by_callback
  
  implicit none
  public

  !
  !> @cond DOXYGEN_IGNORE_THIS
  INTEGER, PARAMETER, PUBLIC :: no_error                                        = 0,          &

       &                calendar_calendartostring                       = 0*100 +  1, &

       &                general_arithmetic_error                        = 0*100 +  2, &

       &                julianday_newjulianday                          = 1*100 +  1, &
       &                julianday_juliandaytostring                     = 1*100 +  3, &

       &                date_newdatefromstring                          = 2*100 +  1, &
       &                date_newdatefromraw_yi8                         = 2*100 +  2, &
       &                date_newdatefromraw                             = 2*100 +  3, &
       &                date_newdatefromconstructandcopy                = 2*100 +  4, &
       &                date_replacedate                                = 2*100 +  6, &
       &                date_datetostring                               = 2*100 +  7, &
       &                date_datetoposixstring                          = 2*100 +  8, &

       &                time_newtimefromstring                          = 3*100 +  1, &
       &                time_newtimefromraw                             = 3*100 +  2, &
       &                time_newtimefromconstructandcopy                = 3*100 +  3, &
       &                time_replacetime                                = 3*100 +  5, &
       &                time_timetostring                               = 3*100 +  6, &
       &                time_timetoposixstring                          = 3*100 +  7, &

       &                datetime_newdatetimefromstring                  = 4*100 +  1, &
       &                datetime_newdatetimefromraw_yi8                 = 4*100 +  2, &
       &                datetime_newdatetimefromraw                     = 4*100 +  3, &
       &                datetime_newdatetimefromconstructandcopy        = 4*100 +  4, &
       &                datetime_datetimetostring                       = 4*100 +  6, &
       &                datetime_datetimetoposixstring                  = 4*100 +  7, &
       &                datetime_getnoofdaysinmonthdatetime             = 4*100 + 15, &
       &                datetime_getnoofdaysinyeardatetime              = 4*100 + 16, &
       &                datetime_getdayofyearfromdatetime               = 4*100 + 17, &
       &                datetime_getnoofsecondselapsedinmonthdatetime   = 4*100 + 18, &
       &                datetime_getnoofsecondselapsedindaydatetime     = 4*100 + 19, &
       &                datetime_getjuliandayfromdatetime               = 4*100 + 20, &
       &                datetime_getdatetimefromjulianday               = 4*100 + 21, &

       &                timedelta_newtimedeltafromstring                = 5*100 +  1, &
       &                timedelta_newtimedeltafromraw                   = 5*100 +  2, &
       &                timedelta_newtimedeltafromraw_yi8               = 5*100 +  3, &
       &                timedelta_newtimedeltafromconstructandcopy      = 5*100 +  4, &
       &                timedelta_gettotalmillisecondstimedelta         = 5*100 +  8, &
       &                timedelta_gettotalsecondstimedelta              = 5*100 +  9, &
       &                timedelta_timedeltatostring                     = 5*100 + 10, &
       &                timedelta_getptstring                           = 5*100 + 11, &

       &                events_newevent                                 = 6*100 +  1, &
       &                events_eventtostring                            = 6*100 +  3, &
       &                events_gettriggernexteventatdatetime            = 6*100 +  8, &
       &                events_gettriggeredpreviouseventatdatetime      = 6*100 +  9, &
       &                events_geteventname                             = 6*100 + 11, &
       &                events_getFirstDatetime                         = 6*100 + 12, &
       &                events_geteventinterval                         = 6*100 + 13, &
       &                events_getlastdatetime                          = 6*100 + 14, &

       &                eventgroups_neweventgroup                       = 7*100 +  1, &
       &                eventgroups_geteventgroupname                   = 7*100 +  6
  !> @endcond DOXYGEN_IGNORE_THIS
  !

contains
  !>
  !! @brief returns error emssage associated to error number of mtime
  !!
  !! @param[in]   errno           error message number
  !!
  !! @param[out]  error_message   associated error message
  subroutine mtime_strerror(errno, error_message)
    integer, intent(in) :: errno
    character(len=*), intent(out) :: error_message

    select case (errno)
    case (no_error)
      error_message = 'no error'

    case (general_arithmetic_error)
      error_message = 'error in arithmetic operation'

    case (calendar_calendartostring)
      error_message = 'could not retrieve the string in <calendartostring>'

    case (julianday_newjulianday)
      error_message = 'could not allocate julianday in <newjulianday>'
    case (julianday_juliandaytostring)
      error_message = 'could not retrieve the string in <juliandaytostring>'

    case (date_newdatefromstring)
      error_message = 'could not allocate date in <newdate>'
    case (date_newdatefromraw_yi8)
      error_message = 'could not allocate date in <newdate>'
    case (date_newdatefromraw)
      error_message = 'could not allocate date in <newdate>'
    case (date_newdatefromconstructandcopy)
      error_message = 'could not allocate date in <newdate>'
    case (date_replacedate)
      error_message = 'replace-date failed in <replacedate>'
    case (date_datetostring)
      error_message = 'could not retrieve the string in <datetostring>'
    case (date_datetoposixstring)
      error_message = 'could not retrieve the string in <datetoposixstring>'

    case (time_newtimefromstring)
      error_message = 'could not allocate time in <newtime>'
    case (time_newtimefromraw)
      error_message = 'could not allocate time in <newtime>'
    case (time_newtimefromconstructandcopy)
      error_message = 'could not allocate time in <newtime>'
    case (time_replacetime)
      error_message = 'replace-time failed in <replacetime>'
    case (time_timetostring)
      error_message = 'could not retrieve the string in <timetostring>'
    case (time_timetoposixstring)
      error_message = 'could not retrieve the string in <timetoposixstring>'

    case (datetime_newdatetimefromstring)
      error_message = 'could not allocate datetime in <newdatetime>'
    case (datetime_newdatetimefromraw_yi8)
      error_message = 'could not allocate datetime in <newdatetime>'
    case (datetime_newdatetimefromraw)
      error_message = 'could not allocate datetime in <newdatetime>'
    case (datetime_newdatetimefromconstructandcopy)
      error_message = 'could not allocate datetime in <newdatetime>'
    case (datetime_datetimetostring)
      error_message = 'could not retrieve the string in <timetostring>'
    case (datetime_datetimetoposixstring)
      error_message = 'could not retrieve the string in <timetoposixstring>'
    case (datetime_getnoofdaysinmonthdatetime)
      error_message = 'could not retrieve the no-of-days-in-month-datetime in <getnoofdaysinmonthdatetime>'
    case (datetime_getnoofdaysinyeardatetime)
      error_message = 'could not retrieve the no-of-days-in-year-datetime in <getnoofdaysinyeardatetime>'
    case (datetime_getdayofyearfromdatetime)
      error_message = 'could not calculate the day-of-year-from-date-time in <getdayofyearfromdatetime>'
    case (datetime_getnoofsecondselapsedinmonthdatetime)
      error_message = 'could not calculate the no-of-seconds-elapsed-in-month-datetime in <getnoofsecondselapsedinmonthdatetime>'
    case (datetime_getnoofsecondselapsedindaydatetime)
      error_message = 'could not calculate the no-of-seconds-elapsed-in-day-datetime in <getnoofsecondselapsedindaydatetime>'
    case (datetime_getjuliandayfromdatetime)
      error_message = 'could not retrieve the get-julian-day-from-datetime in <getjuliandayfromdatetime>'
    case (datetime_getdatetimefromjulianday)
      error_message = 'could not retrieve the get-datetime-from-julian-day in <getdatetimefromjulianday>'

    case (timedelta_newtimedeltafromstring)
      error_message = 'could not allocate timedelta in <newtimedelta>'
    case (timedelta_newtimedeltafromraw)
      error_message = 'could not allocate timedelta in <newtimedelta>'
    case (timedelta_newtimedeltafromraw_yi8)
      error_message = 'could not allocate timedelta in <newtimedelta>'
    case (timedelta_newtimedeltafromconstructandcopy)
      error_message = 'could not allocate timedelta in <newtimedelta>'
    case (timedelta_gettotalmillisecondstimedelta)
      error_message = 'error in calculating total-milli-seconds-timedelta in <gettotalmillisecondstimedelta>'
    case (timedelta_gettotalsecondstimedelta)
      error_message = 'error in calculating total-seconds-in-timedelta in <gettotalmillisecondstimedelta>'
    case (timedelta_timedeltatostring)
      error_message = 'could not retrieve the string in <timedeltatostring>'
    case (timedelta_getptstring)
      error_message = 'could not retrieve the string in <getptstring*>'

    case (events_newevent)
      error_message = 'could not allocate event in <newevent>'
    case (events_eventtostring)
      error_message = 'could not retrieve the string in <eventtostring>'
    case (events_gettriggernexteventatdatetime)
      error_message = 'failed to retrieve trigger-next-event-at-datetime in <gettriggernexteventatdatetime>'
    case (events_gettriggeredpreviouseventatdatetime)
      error_message = 'failed to retrieve triggered-previous-event-at-date-time in <gettriggeredpreviouseventatdatetime>'
    case (events_geteventname)
      error_message = 'failed to retrieve event-name in <geteventname>'
    CASE(events_getFirstDatetime)
      error_message = 'could not get the event first datetime'
    CASE(events_geteventinterval)
      error_message = 'could not get the event interval'
    CASE(events_getlastdatetime)
      error_message = 'could not retrieve the event last date'

    case (eventgroups_neweventgroup)
      error_message = 'could not allocate eventgroup in <neweventgroup>'
    case (eventgroups_geteventgroupname)
      error_message = 'could not retrieve the group-name in <geteventgroupname>'
    end select
  end subroutine mtime_strerror
  !  

end module mtime_error_handling


