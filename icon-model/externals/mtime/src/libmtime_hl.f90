!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
!TODO: extract internal c binding routines (my_...) from libmtime.f90 into extra module
!      and derive old and new interface out of this

!> @file libmtime_hl.f90
!!
!! @brief Providing a high level interface for the Fortran language bindings of libmtime.
!!
!! The mtime library - at least in its first release - heavily uses
!! pointers for passing arguments to functions. This design causes
!! several inconveniences on the application side, e.g. the need for
!! explicit deallocation. Besides, the standard interface of the mtime
!! library provides numerous functions which are well suited for an
!! object-oriented implementation.
!!
!! This wrapper module attempts to provide the mtime functionality
!! with stack-based data structures. It does not refactor the libmtime
!! library itself but hides its allocate-deallocate code within
!! type-bound procedures.
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Florian Prill, DWD
!! @author  J.F. Engels, DKRZ
!!
!! TODOs:
!!  - Event: Wrappers for event functionality; direct use of c-bindings
!!  - Get rid of use of old mtime library
!!  - Expand examples_hl by event wrapper calls
!!  - Make December'18 branch compile, link and run basic test
!!  - Merge recent changes in ICON into our December'18 branch
!!  - Remove old mtime Fortran library
!!  - Change mtime C source, use stack variables
!!
!! @defgroup FortranHLBindings libmtime high-level Fortran language bindings
!! @{
!!
MODULE mtime_hl

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_associated, c_char, c_double, c_f_pointer, c_int, c_int32_t, &
    & c_int64_t, c_loc, c_null_char, c_ptr
  USE mtime_constants, ONLY: &
    & max_datetime_str_len, max_event_str_len, max_eventname_str_len, &
    & max_groupname_str_len, max_timedelta_str_len
  USE mtime_error_handling
  USE mtime_c_bindings, ONLY: &
    & date, datetime, divisionquotienttimespan, handle_errno, julianday, &
    & juliandelta, time, timedelta
  USE mtime_c_bindings, ONLY: &
    & my_addtimedeltatodatetime, my_comparedatetime, my_datetimetostring, &
    & my_datetoposixstring, my_deallocatedatetime, my_deallocatejuliandelta, &
    & my_deallocatetimedelta, my_elementwisescalarmultiplytimedelta, &
    & my_elementwisescalarmultiplytimedeltadp, my_getdayofyearfromdatetime, &
    & my_getjuliandayfromdatetime, my_getnoofdaysinmonthdatetime, &
    & my_getnoofdaysinyeardatetime, my_getnoofsecondselapsedindaydatetime, &
    & my_getnoofsecondselapsedinmonthdatetime, my_gettimedeltafromdate, &
    & my_newdatetime, my_newjuliandelta, my_newrawdatetime, &
    & my_newrawtimedelta, my_newtimedeltafromstring, my_timedeltatostring
  USE mtime, ONLY: &
    & OPERATOR(*), OPERATOR(/=), OPERATOR(<), OPERATOR(<=), OPERATOR(==), &
    & OPERATOR(>), OPERATOR(>=), deallocateDatetime, deallocateJulianDay, &
    & deallocateJulianDelta, deallocateTimedelta, &
    & divideDatetimeDifferenceInSeconds, divideTimeDeltaInSeconds, &
    & getDatetimeFromJulianDay, getPTStringFromMS, &
    & getTotalMilliSecondsTimeDelta, getTotalSecondsTimeDelta, newDatetime, &
    & newJulianDay, newJuliandelta, newTimedelta, timeDeltaToJulianDelta

  IMPLICIT NONE

#ifndef __NVCOMPILER

  PUBLIC :: t_datetime, t_timedelta, t_juliandelta, t_julianday, t_event, t_eventGroup
  PUBLIC :: t_timedeltaFromMilliseconds
  PUBLIC :: t_timedeltaFromSeconds
  PUBLIC :: min, max
  PUBLIC :: OPERATOR(*)

  !
  ! TODO: simply repeat the implementation of "divisionquotienttimespan" in
  !       order to disentangle the mtime_hl and the mtime Fortran modules.
  !
  !PUBLIC :: divisionquotienttimespan

  !> NOTE / TODO:
  !
  !  Why does this wrapper module *not* implement the assignment
  !  operator? - Lenghty answer:
  !
  !    The fundamental question is: Why do our "t_datetime",
  !    "t_timedelta", "t_xxx" types need assignment operators which
  !    apply "pointer magic"? Couldn't we simply copy the stack
  !    variables?
  !
  !    Here's why this is important: Since we are dealing with stack
  !    variables, these might be implicitly initialized with default
  !    values, causing problems with our "pointer magic". Imagine that
  !    in our application code we have a derived type as follows:
  !
  !      TYPE(mytype)
  !        [...]
  !        TYPE(t_datetime) :: mtime_date
  !      END TYPE
  !
  !    where "mtime_date" is an mtime datetime variable which is
  !    usually not used and therefore not explicitly initialized by
  !    the user. Now, when we have this operation:
  !
  !      TYPE(mytype) :: tmp_a, tmp_b
  !      tmp_a = tmp_b
  !
  !    then the application code crashes, because the assign operator
  !    of "t_datetime" attempts to create an mtime "datetime" pointer
  !    based on the uninitializd "tmp_b%dt".
  !

  !> Wrapper class for "mtime" data type "datetime".
  !
  TYPE t_datetime

    !private

    !> wrapped datatype of the standard mtime interface
    TYPE(datetime) :: dt = datetime(date(year=0_c_int64_t, month=0_c_int, day=0_c_int),&
         &                          time(hour=0_c_int, minute=0_c_int, second=0_c_int, ms=0_c_int))

  CONTAINS

    ! --- conversions

    PROCEDURE :: t_datetime_toString
    PROCEDURE :: t_datetime_to_posix_string
    GENERIC   :: toString => t_datetime_toString, t_datetime_to_posix_string

    PROCEDURE :: toJulianDay => t_datetime_toJulianDay

    ! --- inquire components

    PROCEDURE :: getDay => t_datetime_getDay

    ! --- derived quantities

    PROCEDURE :: daysInEntireMonth => t_datetime_daysInEntireMonth
    PROCEDURE :: daysInEntireYear => t_datetime_daysInEntireYear
    PROCEDURE :: elapsedDaysInYear => t_datetime_elapsedDaysInYear
    PROCEDURE :: elapsedSecondsInMonth => t_datetime_elapsedSecondsInMonth
    PROCEDURE :: elapsedSecondsInDay => t_datetime_elapsedSecondsInDay

    ! --- overloaded operators

    PROCEDURE :: add_timedelta => t_datetime_add_timedelta
    PROCEDURE :: sub_timedelta => t_datetime_sub_timedelta
    PROCEDURE :: sub_datetime => t_datetime_sub_datetime
    PROCEDURE :: equal_datetime => t_datetime_equal
    PROCEDURE :: not_equal_datetime => t_datetime_not_equal
    PROCEDURE :: less_than_datetime => t_datetime_less_than
    PROCEDURE :: greater_than_datetime => t_datetime_greater_than
    PROCEDURE :: less_or_equal_datetime => t_datetime_less_or_equal
    PROCEDURE :: greater_or_equal_datetime => t_datetime_greater_or_equal

    PROCEDURE :: get_c_pointer => t_datetime_get_c_pointer

    GENERIC   :: OPERATOR(+) => add_timedelta
    GENERIC   :: OPERATOR(-) => sub_timedelta
    GENERIC   :: OPERATOR(-) => sub_datetime
    GENERIC   :: OPERATOR(==) => equal_datetime
    GENERIC   :: OPERATOR(/=) => not_equal_datetime
    GENERIC   :: OPERATOR(<) => less_than_datetime
    GENERIC   :: OPERATOR(>) => greater_than_datetime
    GENERIC   :: OPERATOR(<=) => less_or_equal_datetime
    GENERIC   :: OPERATOR(>=) => greater_or_equal_datetime

  END TYPE t_datetime

  INTERFACE t_datetime
    MODULE PROCEDURE t_datetime_assign_string
    MODULE PROCEDURE t_datetime_assign_raw
  END INTERFACE t_datetime

  !> Wrapper class for "mtime" data type "timedelta".
  !
  TYPE t_timedelta

    !private

    !> wrapped datatype of the standard mtime interface
    TYPE(timedelta) :: td = timedelta(flag_std_form=0_c_int, sign='+', year=0_c_int64_t,   &
         &                                month=0_c_int, day=0_c_int, hour=0_c_int,         &
         &                                minute=0_c_int, second=0_c_int, ms=0_c_int)

  CONTAINS

    ! --- conversions

    PROCEDURE :: toString => t_timedelta_toString
    PROCEDURE :: toJulianDelta => t_timedelta_toJulianDelta

    ! --- inquire components

    ! todo: [...]

    ! --- derived quantities

    PROCEDURE :: toSeconds => t_timedelta_toSeconds

    ! t_timedelta_toMilliSeconds: todo: It would be convenient to have
    ! the reference date for this function as an optional argument;
    ! only in case of "non-well-definedness" an error should be
    ! thrown.
    PROCEDURE :: toMilliSeconds => t_timedelta_toMilliSeconds

    PROCEDURE :: divideInSecondsBy => t_timedelta_divideInSecondsBy

    ! --- overloaded operators

    ! note: the "+", "-" operators are not well-defined for timedelta
    ! objects!

    PROCEDURE :: equal_datetime => t_timedelta_equal
    PROCEDURE :: not_equal_datetime => t_timedelta_not_equal
    PROCEDURE :: less_than_datetime => t_timedelta_less_than
    PROCEDURE :: greater_than_datetime => t_timedelta_greater_than
    PROCEDURE :: less_or_equal_datetime => t_timedelta_less_than_or_equal
    PROCEDURE :: greater_or_equal_datetime => t_timedelta_greater_than_or_equal
    PROCEDURE :: scalar_multiply_long => t_timedelta_scalar_multiply_long
    PROCEDURE :: scalar_multiply_int => t_timedelta_scalar_multiply_int
    PROCEDURE :: scalar_multiply_real => t_timedelta_scalar_multiply_real

    PROCEDURE :: get_c_pointer => t_timedelta_get_c_pointer

    GENERIC   :: OPERATOR(==) => equal_datetime
    GENERIC   :: OPERATOR(/=) => not_equal_datetime
    GENERIC   :: OPERATOR(<) => less_than_datetime
    GENERIC   :: OPERATOR(>) => greater_than_datetime
    GENERIC   :: OPERATOR(<=) => less_or_equal_datetime
    GENERIC   :: OPERATOR(>=) => greater_or_equal_datetime
    GENERIC   :: OPERATOR(*) => scalar_multiply_long, scalar_multiply_int,         &
         &                                        scalar_multiply_real

  END TYPE t_timedelta

  INTERFACE t_timedelta
    MODULE PROCEDURE t_timedelta_assign_string
  END INTERFACE t_timedelta

  INTERFACE t_timedeltaFromMilliseconds
    MODULE PROCEDURE t_timedelta_assign_ms
    MODULE PROCEDURE t_timedelta_assign_ms_i8
  END INTERFACE t_timedeltaFromMilliseconds

  INTERFACE t_timedeltaFromSeconds
    MODULE PROCEDURE t_timedelta_assign_sec
    MODULE PROCEDURE t_timedelta_assign_sec_i8
  END INTERFACE t_timedeltaFromSeconds

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE t_timedelta_scalar_multiply_inv_long
    MODULE PROCEDURE t_timedelta_scalar_multiply_inv_int
    MODULE PROCEDURE t_timedelta_scalar_multiply_inv_real
  END INTERFACE OPERATOR(*)

  INTERFACE min
    MODULE PROCEDURE t_datetime_min
  END INTERFACE min

  INTERFACE max
    MODULE PROCEDURE t_datetime_max
  END INTERFACE max

  TYPE t_julianday

    !private

    !> wrapped datatype of the standard mtime interface
    TYPE(julianday) :: jd

  CONTAINS

    ! --- conversions

    PROCEDURE :: toDateTime => t_julianday_toDateTime

    ! --- inquire components

    PROCEDURE :: getDay => t_julianDay_getDay
    PROCEDURE :: getFractionOfDayInMS => t_julianday_getFractionOfDayInMS

  END TYPE t_julianday

  !> Wrapper class for "mtime" data type "juliandelta".
  !
  TYPE t_juliandelta

    !private

    !> wrapped datatype of the standard mtime interface
    TYPE(juliandelta) :: jd

  END TYPE t_juliandelta

  INTERFACE t_juliandelta
    MODULE PROCEDURE t_juliandelta_assign_raw
  END INTERFACE t_juliandelta

  TYPE t_event

    !private
    ! FIXME (some day in the future): This derived type should not specify both -
    ! the linked list element and the element itself.
    TYPE(t_event), POINTER               :: nextEventInGroup

    INTEGER(c_int64_t)                   :: eventId
    CHARACTER(len=max_eventname_str_len) :: eventName

    TYPE(t_datetime)                     :: eventsLastEvaluationDateTime
    TYPE(t_datetime)                     :: eventReferenceDateTime

    TYPE(t_datetime)                     :: eventFirstDateTime
    TYPE(t_datetime)                     :: eventLastDateTime

    TYPE(t_timedelta)                    :: eventInterval
    TYPE(t_timedelta)                    :: eventOffset

    LOGICAL                              :: neverTriggerEvent

    LOGICAL                              :: triggerCurrentEvent

    LOGICAL                              :: nextEventIsFirst
    LOGICAL                              :: lastEventWasFinal

    LOGICAL                              :: eventisFirstInDay
    LOGICAL                              :: eventisFirstInMonth
    LOGICAL                              :: eventisFirstInYear
    LOGICAL                              :: eventisLastInDay
    LOGICAL                              :: eventisLastInMonth
    LOGICAL                              :: eventisLastInYear

    TYPE(t_datetime)                     :: triggerNextEventDateTime
    TYPE(t_datetime)                     :: triggeredPreviousEventDateTime

  CONTAINS

    !> TODO: implement isAvtive ....
    ! PROCEDURE :: trigger                   => t_event_trigger
    PROCEDURE :: getFirstDatetime => t_event_getFirstDatetime
    PROCEDURE :: getInterval => t_event_getInterval
    PROCEDURE :: getLastDatetime => t_event_getLastDatetime
    PROCEDURE :: getNextOccurrenceDatetime => t_event_getNextOccurrenceDatetime
    PROCEDURE :: getPrevOccurrenceDatetime => t_event_getPrevOccurrenceDatetime

    PROCEDURE :: getName => t_event_getName
    PROCEDURE :: nextEvent => t_event_next_event

  END TYPE t_event

  INTERFACE t_event
    MODULE PROCEDURE t_event_assign_raw
    MODULE PROCEDURE t_event_assign_types
  END INTERFACE t_event

  TYPE t_eventGroup

    !private

    INTEGER(c_int64_t)                   :: event_group_id
    CHARACTER(len=max_groupname_str_len) :: event_group_name
    TYPE(t_event), POINTER               :: first_event_in_group
    TYPE(t_event), POINTER               :: last_event_in_group

  CONTAINS

    PROCEDURE :: append => t_eventGroup_addToGroup

    !> TODO: implement the removal of a event in a list
    !PROCEDURE :: remove        => t_eventGroup_removeFromGroup

    ! --- inquire components

    PROCEDURE :: getID => t_eventGroup_getGroupId
    PROCEDURE :: getName => t_eventGroup_getGroupName

    ! --- derived quantities

    PROCEDURE :: getFirstEvent => t_eventGroup_getFirstEvent

  END TYPE t_eventGroup

  INTERFACE t_eventGroup
    MODULE PROCEDURE t_eventGroup_constructor
  END INTERFACE t_eventGroup

  INTEGER :: event_group_id = 0
  INTEGER :: event_id = 0

CONTAINS

#include "mtime_t_datetime.inc"
#include "mtime_t_timedelta.inc"
#include "mtime_t_juliandelta.inc"
#include "mtime_t_event.inc"

#endif

END MODULE mtime_hl
!>
!! @}
