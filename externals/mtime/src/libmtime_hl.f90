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

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int32_t, c_int64_t, c_double, &
       &                                 c_null_char, c_ptr, c_loc, c_f_pointer
  USE mtime_c_bindings
  USE mtime_constants
  use mtime_error_handling
  USE mtime

  IMPLICIT NONE


  
  PUBLIC :: t_datetime, t_timedelta, t_juliandelta, t_julianday, t_event, t_eventGroup
  PUBLIC :: t_timedeltaFromMilliseconds
  PUBLIC :: t_timedeltaFromSeconds
  PUBLIC :: min, max
  PUBLIC :: OPERATOR(*)

  ! Re-export stuff
  public :: register_finish_mtime_procedure, finish_mtime  

  !
  ! TODO: simply repeat the implementation of "divisionquotienttimespan" in order to disentangle the mtime_hl and the mtime Fortran modules.
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
    GENERIC   :: toString                          => t_datetime_toString, t_datetime_to_posix_string

    PROCEDURE :: toJulianDay                       => t_datetime_toJulianDay

    ! --- inquire components

    PROCEDURE :: getDay                            => t_datetime_getDay

    ! --- derived quantities

    PROCEDURE :: daysInEntireMonth                 => t_datetime_daysInEntireMonth
    PROCEDURE :: daysInEntireYear                  => t_datetime_daysInEntireYear
    PROCEDURE :: elapsedDaysInYear                 => t_datetime_elapsedDaysInYear
    PROCEDURE :: elapsedSecondsInMonth             => t_datetime_elapsedSecondsInMonth
    PROCEDURE :: elapsedSecondsInDay               => t_datetime_elapsedSecondsInDay

    ! --- overloaded operators

    PROCEDURE :: add_timedelta                     => t_datetime_add_timedelta
    PROCEDURE :: sub_timedelta                     => t_datetime_sub_timedelta
    PROCEDURE :: sub_datetime                      => t_datetime_sub_datetime
    PROCEDURE :: equal_datetime                    => t_datetime_equal
    PROCEDURE :: not_equal_datetime                => t_datetime_not_equal
    PROCEDURE :: less_than_datetime                => t_datetime_less_than
    PROCEDURE :: greater_than_datetime             => t_datetime_greater_than
    PROCEDURE :: less_or_equal_datetime            => t_datetime_less_or_equal
    PROCEDURE :: greater_or_equal_datetime         => t_datetime_greater_or_equal

    PROCEDURE :: get_c_pointer                     => t_datetime_get_c_pointer

    GENERIC   :: OPERATOR(+)                       =>  add_timedelta
    GENERIC   :: OPERATOR(-)                       =>  sub_timedelta
    GENERIC   :: OPERATOR(-)                       =>  sub_datetime
    GENERIC   :: OPERATOR(==)                      =>  equal_datetime
    GENERIC   :: OPERATOR(/=)                      =>  not_equal_datetime
    GENERIC   :: OPERATOR(<)                       =>  less_than_datetime
    GENERIC   :: OPERATOR(>)                       =>  greater_than_datetime
    GENERIC   :: OPERATOR(<=)                      =>  less_or_equal_datetime
    GENERIC   :: OPERATOR(>=)                      =>  greater_or_equal_datetime

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
     TYPE(timedelta) :: td = timedelta( flag_std_form = 0_c_int, sign = '+', year = 0_c_int64_t,   &
          &                                month = 0_c_int, day = 0_c_int, hour = 0_c_int,         &
          &                                minute = 0_c_int, second = 0_c_int, ms = 0_c_int )


  CONTAINS

    ! --- conversions

    PROCEDURE :: toString                  => t_timedelta_toString
    PROCEDURE :: toJulianDelta             => t_timedelta_toJulianDelta


    ! --- inquire components

    ! todo: [...]


    ! --- derived quantities

    PROCEDURE :: toSeconds                 => t_timedelta_toSeconds

    ! t_timedelta_toMilliSeconds: todo: It would be convenient to have
    ! the reference date for this function as an optional argument;
    ! only in case of "non-well-definedness" an error should be
    ! thrown.
    PROCEDURE :: toMilliSeconds            => t_timedelta_toMilliSeconds

    PROCEDURE :: divideInSecondsBy         => t_timedelta_divideInSecondsBy


    ! --- overloaded operators

    ! note: the "+", "-" operators are not well-defined for timedelta
    ! objects!

    PROCEDURE :: equal_datetime            => t_timedelta_equal
    PROCEDURE :: not_equal_datetime        => t_timedelta_not_equal
    PROCEDURE :: less_than_datetime        => t_timedelta_less_than
    PROCEDURE :: greater_than_datetime     => t_timedelta_greater_than
    PROCEDURE :: less_or_equal_datetime    => t_timedelta_less_than_or_equal
    PROCEDURE :: greater_or_equal_datetime => t_timedelta_greater_than_or_equal
    PROCEDURE :: scalar_multiply_long      => t_timedelta_scalar_multiply_long
    PROCEDURE :: scalar_multiply_int       => t_timedelta_scalar_multiply_int
    PROCEDURE :: scalar_multiply_real      => t_timedelta_scalar_multiply_real

    PROCEDURE :: get_c_pointer             => t_timedelta_get_c_pointer

    GENERIC   :: OPERATOR(==)              =>  equal_datetime
    GENERIC   :: OPERATOR(/=)              =>  not_equal_datetime
    GENERIC   :: OPERATOR(<)               =>  less_than_datetime
    GENERIC   :: OPERATOR(>)               =>  greater_than_datetime
    GENERIC   :: OPERATOR(<=)              =>  less_or_equal_datetime
    GENERIC   :: OPERATOR(>=)              =>  greater_or_equal_datetime
    GENERIC   :: OPERATOR(*)               =>  scalar_multiply_long, scalar_multiply_int,         &
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

    PROCEDURE :: toDateTime           => t_julianday_toDateTime

    ! --- inquire components

    PROCEDURE :: getDay               => t_julianDay_getDay
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
    PROCEDURE :: getFirstDatetime          => t_event_getFirstDatetime
    PROCEDURE :: getInterval               => t_event_getInterval
    PROCEDURE :: getLastDatetime           => t_event_getLastDatetime
    PROCEDURE :: getNextOccurrenceDatetime => t_event_getNextOccurrenceDatetime
    PROCEDURE :: getPrevOccurrenceDatetime => t_event_getPrevOccurrenceDatetime

    PROCEDURE :: getName                   => t_event_getName
    PROCEDURE :: nextEvent                 => t_event_next_event
    PROCEDURE :: isActive                  => t_event_is_active

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

    PROCEDURE :: append        => t_eventGroup_addToGroup

    !> TODO: implement the removal of a event in a list
    !PROCEDURE :: remove        => t_eventGroup_removeFromGroup

    ! --- inquire components

    PROCEDURE :: getID         => t_eventGroup_getGroupId
    PROCEDURE :: getName       => t_eventGroup_getGroupName

    ! --- derived quantities

    PROCEDURE :: getFirstEvent => t_eventGroup_getFirstEvent

  END TYPE t_eventGroup

  INTERFACE t_eventGroup
    MODULE PROCEDURE t_eventGroup_constructor
  END INTERFACE t_eventGroup


  INTEGER :: event_group_id = 0
  INTEGER :: event_id       = 0

CONTAINS

!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
  ! ================================================================================
  ! datetime section:
  ! ================================================================================


  ! constructor for a datetime string
  !
  TYPE(t_datetime) FUNCTION t_datetime_assign_string(dt_string)
    CHARACTER(len=*), INTENT(in) :: dt_string
    TYPE(c_ptr)             :: c_pointer
    TYPE(datetime), POINTER :: dt_tmp
    INTEGER :: errno
    c_pointer = my_newdatetime(TRIM(ADJUSTL(dt_string))//c_null_char)
    CALL handle_errno(.NOT. c_ASSOCIATED(c_pointer), 4 *100 + 1, &
      & "externals/mtime/src/mtime_t_datetime.inc", &
      & 20)
    CALL c_f_pointer(c_pointer, dt_tmp)
    t_datetime_assign_string%dt = dt_tmp
    CALL my_deallocatedatetime(c_pointer)
  END FUNCTION t_datetime_assign_string

  ! constructor for a datetime
  !
  FUNCTION t_datetime_assign_raw(year, month, day, hour, minute, second, ms)  RESULT(res)
    TYPE(t_datetime) :: res

    INTEGER(c_int64_t), INTENT(in) :: year
    INTEGER(c_int),     INTENT(in) :: month, day, hour, minute, second, ms

    TYPE(datetime), POINTER :: dt_tmp
    TYPE(c_ptr)             :: c_pointer

    c_pointer = my_newrawdatetime(year, month, day, hour, minute, second, ms)
    CALL handle_errno(.NOT. c_ASSOCIATED(c_pointer), 4 *100 + 2, &
      & "externals/mtime/src/mtime_t_datetime.inc", &
      & 40)
    call c_f_pointer(c_pointer, dt_tmp)
    res%dt = dt_tmp
    CALL my_deallocatedatetime(c_pointer)
  END FUNCTION t_datetime_assign_raw


  ! Returns t_datetime objects day
  !
  FUNCTION t_datetime_getDay(this)
    INTEGER            :: t_datetime_getDay
    CLASS (t_datetime) :: this
    t_datetime_getDay = this%dt%date%day
  END FUNCTION t_datetime_getDay

  ! Convert t_datetime object to string.
  !
  FUNCTION t_datetime_toString(this) result(string)
    CHARACTER(len=max_datetime_str_len)  :: string
    CLASS (t_datetime)                   :: this
    type(c_ptr) :: c_pointer, dummy_ptr
    integer :: i

    string = ""
    
    c_pointer = this%get_c_pointer()
    CALL handle_errno(.not. c_associated(c_pointer), 0 * 100 + 2, &
      & "externals/mtime/src/mtime_t_datetime.inc", &
      & 68)
    dummy_ptr = my_datetimetostring(c_pointer, string)
    CALL handle_errno(.not. c_associated(dummy_ptr), 4 * 100 + 6, &
      & "externals/mtime/src/mtime_t_datetime.inc", &
      & 72)

    char_loop: do i = 1 , len(string)
       if (string(i:i) == c_null_char) exit char_loop
    end do char_loop
    string(i:len(string)) = ' '

    CALL my_deallocatedatetime(c_pointer)
  END FUNCTION t_datetime_toString

  ! Convert t_datetime object to string.
  !
  FUNCTION t_datetime_to_posix_string(this, format_string) result(string)
    CHARACTER(len=max_datetime_str_len)  :: string
    CHARACTER(len=*), INTENT(in)         :: format_string
    CLASS (t_datetime)                   :: this
    integer :: i
    type(c_ptr) :: c_pointer, dummy_ptr

    string = ""
    c_pointer = this%get_c_pointer()
    dummy_ptr = my_datetoposixstring(c_pointer, string, format_string)
    CALL handle_errno(.not. c_associated(dummy_ptr), 2 * 100 + 8, &
      & "externals/mtime/src/mtime_t_datetime.inc", &
      & 96)
    char_loop: do i = 1 , len(string)
      if (string(i:i) == c_null_char) exit char_loop
    end do char_loop
    string(i:len(string)) = ' '
    
    CALL my_deallocatedatetime(c_pointer)
  END FUNCTION t_datetime_to_posix_string

  FUNCTION t_datetime_toJulianDay(this) RESULT(jd)
    CLASS(t_datetime), INTENT(in) :: this
    TYPE(t_julianday), target :: jd
    type(c_ptr) :: c_pointer, dummy_ptr
    c_pointer = this%get_c_pointer()
    dummy_ptr = my_getjuliandayfromdatetime(c_pointer, c_loc(jd%jd))
    CALL handle_errno(.not. c_associated(dummy_ptr), 0 * 100 + 2, &
      & "externals/mtime/src/mtime_t_datetime.inc", &
      & 113)
    CALL my_deallocatedatetime(c_pointer)
  END FUNCTION t_datetime_toJulianDay

  ! Addition of time interval to datetime object.
  !
  FUNCTION t_datetime_add_timedelta(this, td) RESULT(dt_td_sum)
    TYPE (t_datetime)               :: dt_td_sum
    CLASS (t_datetime),  INTENT(in) :: this
    CLASS (t_timedelta), INTENT(in) :: td
    TYPE(datetime),  POINTER        :: dt_tmp

    type(c_ptr) :: c_pointer1, c_pointer2, dummy_ptr
    c_pointer1 = this%get_c_pointer()
    c_pointer2 = td%get_c_pointer()

    dummy_ptr = my_addtimedeltatodatetime(c_pointer1, c_pointer2, c_pointer1)
    CALL handle_errno(.not. c_associated(dummy_ptr), 0 * 100 + 2, &
      & "externals/mtime/src/mtime_t_datetime.inc", &
      & 132)
    call c_f_pointer(c_pointer1, dt_tmp)
    dt_td_sum%dt = dt_tmp
    CALL my_deallocatedatetime(c_pointer1)
    CALL my_deallocatedatetime(c_pointer2)    
  END FUNCTION t_datetime_add_timedelta

  ! Subtraction of time interval to datetime object.
  !
  FUNCTION t_datetime_sub_timedelta(this, td) RESULT(dt_td_sum)
    TYPE (t_datetime)               :: dt_td_sum
    CLASS (t_datetime),  INTENT(in) :: this
    type  (t_timedelta), INTENT(in) :: td
    TYPE(t_timedelta)              :: td_tmp
    TYPE(datetime), pointer        :: dt_tmp
    type(c_ptr) :: c_pointer1, c_pointer2, dummy_ptr

    td_tmp = td
    IF (td_tmp%td%sign == "+") THEN
      td_tmp%td%sign = "-"
    ELSE
      td_tmp%td%sign = "+"
    ENDIF
    
    c_pointer1 = this%get_c_pointer()
    c_pointer2 = td_tmp%get_c_pointer()

    dummy_ptr = my_addtimedeltatodatetime(c_pointer1, c_pointer2, c_pointer1)
    CALL handle_errno(.not. c_associated(dummy_ptr), 0 * 100 + 2, &
      & "externals/mtime/src/mtime_t_datetime.inc", &
      & 162)
    call c_f_pointer(c_pointer1, dt_tmp)
    dt_td_sum%dt = dt_tmp
    CALL my_deallocatedatetime(c_pointer1)
    CALL my_deallocatedatetime(c_pointer2)    
  END FUNCTION t_datetime_sub_timedelta

  ! Subtraction of two dates.
  !
  FUNCTION t_datetime_sub_datetime(this, dt) RESULT(dt_dt_diff)
    TYPE (t_timedelta), target :: dt_dt_diff
    CLASS (t_datetime),  INTENT(in), target :: this
    CLASS (t_datetime),  INTENT(in), target :: dt
    type(c_ptr) :: dummy_ptr
    dummy_ptr = my_gettimedeltafromdate(c_loc(this%dt),c_loc(dt%dt),c_loc(dt_dt_diff%td))
  END FUNCTION t_datetime_sub_datetime

  ! Overloaded operator: test for equivalence.
  !
  LOGICAL FUNCTION t_datetime_equal(this, dt) result(eq)
    CLASS (t_datetime),  INTENT(in), target :: this
    CLASS (t_datetime),  INTENT(in), target :: dt
    integer(c_int) :: ret
    ret = my_comparedatetime(c_loc(this%dt), c_loc(dt%dt))
    if (ret == 0) then
      eq = .true.
    else
      eq = .false.
    endif    
  END FUNCTION t_datetime_equal

  LOGICAL FUNCTION t_datetime_not_equal(this, dt)
    CLASS (t_datetime),  INTENT(in) :: this
    CLASS (t_datetime),  INTENT(in) :: dt
    t_datetime_not_equal = .not. (this%dt == dt%dt)
  END FUNCTION t_datetime_not_equal

  LOGICAL FUNCTION t_datetime_less_than(this, dt)  result(lt)
    CLASS (t_datetime),  INTENT(in), target :: this
    CLASS (t_datetime),  INTENT(in), target :: dt
    integer(c_int) :: ret
    ret = my_comparedatetime(c_loc(this%dt), c_loc(dt%dt))
    if (ret == -1) then
      lt = .true.
    else
      lt = .false.
    endif    
  END FUNCTION t_datetime_less_than

  LOGICAL FUNCTION t_datetime_greater_than(this, dt)  result(gt)
    CLASS (t_datetime),  INTENT(in), target :: this
    CLASS (t_datetime),  INTENT(in), target :: dt
    integer(c_int) :: ret
    ret = my_comparedatetime(c_loc(this%dt), c_loc(dt%dt))
    if (ret == 1) then
      gt = .true.
    else
      gt = .false.
    endif    
  END FUNCTION t_datetime_greater_than

  LOGICAL FUNCTION t_datetime_less_or_equal(this, dt)
    CLASS (t_datetime),  INTENT(in) :: this
    CLASS (t_datetime),  INTENT(in) :: dt
    t_datetime_less_or_equal = .not. (this > dt)
  END FUNCTION t_datetime_less_or_equal

  LOGICAL FUNCTION t_datetime_greater_or_equal(this, dt)
    CLASS (t_datetime),  INTENT(in) :: this
    CLASS (t_datetime),  INTENT(in) :: dt
    t_datetime_greater_or_equal = .not. (this < dt)
  END FUNCTION t_datetime_greater_or_equal

  FUNCTION t_datetime_daysInEntireMonth(this)
    CLASS (t_datetime), INTENT(in), target :: this
    INTEGER(c_int) :: t_datetime_daysInEntireMonth
    t_datetime_daysInEntireMonth = my_getnoofdaysinmonthdatetime(c_loc(this%dt))
    call handle_errno(t_datetime_daysInEntireMonth == 0, 4 * 100 + 15, &
      & "externals/mtime/src/mtime_t_datetime.inc", 240)
  END FUNCTION t_datetime_daysInEntireMonth

  FUNCTION t_datetime_daysInEntireYear(this)
    CLASS (t_datetime),  INTENT(in), target :: this
    INTEGER(c_int) :: t_datetime_daysInEntireYear
    t_datetime_daysInEntireYear = my_getnoofdaysinyeardatetime(c_loc(this%dt))
    CALL handle_errno(t_datetime_daysInEntireYear == 0, 4 * 100 + 16, &
      & "externals/mtime/src/mtime_t_datetime.inc", 248)
  END FUNCTION t_datetime_daysInEntireYear

  FUNCTION t_datetime_elapsedDaysInYear(this)
    CLASS (t_datetime),  INTENT(in), target :: this
    INTEGER(c_int) :: t_datetime_elapsedDaysInYear

    t_datetime_elapsedDaysInYear = my_getdayofyearfromdatetime(c_loc(this%dt))
    CALL handle_errno(t_datetime_elapsedDaysInYear == 0, 4 * 100 + 17, &
         &            "externals/mtime/src/mtime_t_datetime.inc", 257)
  END FUNCTION t_datetime_elapsedDaysInYear

  FUNCTION t_datetime_elapsedSecondsInMonth(this)
    CLASS (t_datetime),  INTENT(in), target :: this
    INTEGER(c_int64_t) :: t_datetime_elapsedSecondsInMonth
    t_datetime_elapsedSecondsInMonth = my_getnoofsecondselapsedinmonthdatetime(c_loc(this%dt))
    CALL handle_errno(t_datetime_elapsedSecondsInMonth == -1, 4 * 100 + 18, &
      & "externals/mtime/src/mtime_t_datetime.inc", 265)
  END FUNCTION t_datetime_elapsedSecondsInMonth

  FUNCTION t_datetime_elapsedSecondsInDay(this)
    CLASS (t_datetime),  INTENT(in), target :: this
    INTEGER(c_int64_t) :: t_datetime_elapsedSecondsInDay

    t_datetime_elapsedSecondsInDay = my_getnoofsecondselapsedindaydatetime(c_loc(this%dt))
    CALL handle_errno(t_datetime_elapsedSecondsInDay == -1, 4 * 100 + 19, &
     & "externals/mtime/src/mtime_t_datetime.inc", 274)
  END FUNCTION t_datetime_elapsedSecondsInDay


  FUNCTION t_datetime_get_c_pointer(this) RESULT(c_pointer)
    TYPE(c_ptr) :: c_pointer
    CLASS(t_datetime) :: this
    c_pointer = my_newrawdatetime(int(this%dt%date%year,c_int64_t), this%dt%date%month,     &
         &                        this%dt%date%day, this%dt%time%hour, this%dt%time%minute, &
         &                        this%dt%time%second, this%dt%time%ms)
    call handle_errno((.not. c_associated(c_pointer)), 4 * 100 + 3, &
      & "externals/mtime/src/mtime_t_datetime.inc", &
      & 286)
  END FUNCTION t_datetime_get_c_pointer


!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
  ! ================================================================================
  ! timedelta section:
  ! ================================================================================

! ToDo: Most of this file needs to be renovated for NOT USE mtime.

  ! constructor for timedelta string
  !
  TYPE(t_timedelta) FUNCTION t_timedelta_assign_string(td_string) !OK-tested
    CHARACTER(len=*), INTENT(in) :: td_string
    TYPE(timedelta), POINTER     :: td_tmp
    type(c_ptr) :: c_pointer

    c_pointer = my_newtimedeltafromstring(trim(adjustl(td_string))//c_null_char)
    if (.not. c_associated(c_pointer)) then
       call handle_errno(5 * 100 + 1, &
         & "externals/mtime/src/mtime_t_timedelta.inc", &
         & 22)
       t_timedelta_assign_string%td%sign = '?'
    else
       call c_f_pointer(c_pointer, td_tmp)
       t_timedelta_assign_string%td = td_tmp
       t_timedelta_assign_string%td%sign = td_tmp%sign
       call my_deallocatetimedelta(c_pointer)
    end if
  END FUNCTION t_timedelta_assign_string

  ! constructor for integer milliseconds (integer)
  !
  TYPE(t_timedelta) FUNCTION t_timedelta_assign_ms(td_ms)
    INTEGER, INTENT(in)                  :: td_ms
    TYPE(timedelta), POINTER             :: td_tmp
    CHARACTER(len=max_timedelta_str_len) :: td_string
    INTEGER                              :: errno
    CALL getptstringfromms(INT(td_ms,c_int64_t), td_string, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 42)
    td_tmp => newtimedelta(td_string, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 46)
    t_timedelta_assign_ms%td = td_tmp
    t_timedelta_assign_ms%td%sign = td_tmp%sign
    CALL deallocatetimedelta(td_tmp)
  END FUNCTION t_timedelta_assign_ms

  ! constructor for integer milliseconds (integer)
  !
  TYPE(t_timedelta) FUNCTION t_timedelta_assign_ms_i8(td_ms)
    INTEGER(c_int64_t), INTENT(in)              :: td_ms
    TYPE(timedelta), POINTER             :: td_tmp
    CHARACTER(len=max_timedelta_str_len) :: td_string
    INTEGER                              :: errno
    CALL getptstringfromms(INT(td_ms,c_int64_t), td_string, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 62)
    td_tmp => newtimedelta(td_string, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 66)
    t_timedelta_assign_ms_i8%td = td_tmp
    t_timedelta_assign_ms_i8%td%sign = td_tmp%sign
    CALL deallocatetimedelta(td_tmp)
  END FUNCTION t_timedelta_assign_ms_i8

  ! constructor for integer seconds (integer)
  !
  TYPE(t_timedelta) FUNCTION t_timedelta_assign_sec(td_sec)
    INTEGER, INTENT(in) :: td_sec
    t_timedelta_assign_sec = t_timedelta_assign_ms(td_sec*1000)
  END FUNCTION t_timedelta_assign_sec

  ! constructor for integer seconds (integer)
  !
  TYPE(t_timedelta) FUNCTION t_timedelta_assign_sec_i8(td_sec)
    INTEGER(c_int64_t), INTENT(in) :: td_sec
    t_timedelta_assign_sec_i8 = t_timedelta_assign_ms_i8(td_sec*1000)
  END FUNCTION t_timedelta_assign_sec_i8


  LOGICAL FUNCTION t_timedelta_equal(this, td)
    CLASS (t_timedelta),  INTENT(in) :: this
    CLASS (t_timedelta),  INTENT(in) :: td
    t_timedelta_equal = (this%td == td%td)
  END FUNCTION t_timedelta_equal

  LOGICAL FUNCTION t_timedelta_not_equal(this, td)
    CLASS (t_timedelta),  INTENT(in) :: this
    CLASS (t_timedelta),  INTENT(in) :: td
    t_timedelta_not_equal = (this%td /= td%td)
  END FUNCTION t_timedelta_not_equal

  LOGICAL FUNCTION t_timedelta_less_than(this, td)
    CLASS (t_timedelta),  INTENT(in) :: this
    CLASS (t_timedelta),  INTENT(in) :: td
    t_timedelta_less_than = (this%td < td%td)
  END FUNCTION t_timedelta_less_than

  LOGICAL FUNCTION t_timedelta_greater_than(this, td)
    CLASS (t_timedelta),  INTENT(in) :: this
    CLASS (t_timedelta),  INTENT(in) :: td
    t_timedelta_greater_than = (this%td > td%td)
  END FUNCTION t_timedelta_greater_than

  LOGICAL FUNCTION t_timedelta_less_than_or_equal(this, td)
    CLASS (t_timedelta),  INTENT(in) :: this
    CLASS (t_timedelta),  INTENT(in) :: td
    t_timedelta_less_than_or_equal = (this%td <= td%td)
  END FUNCTION t_timedelta_less_than_or_equal

  LOGICAL FUNCTION t_timedelta_greater_than_or_equal(this, td)
    CLASS (t_timedelta),  INTENT(in) :: this
    CLASS (t_timedelta),  INTENT(in) :: td
    t_timedelta_greater_than_or_equal = (this%td >= td%td)
  END FUNCTION t_timedelta_greater_than_or_equal


  FUNCTION t_timedelta_scalar_multiply_long (this, lambda) RESULT(scaled_td)
    TYPE(t_timedelta), TARGET :: scaled_td
    INTEGER(c_int64_t),       INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this
    TYPE(timedelta), POINTER             :: td_tmp, td_tmp2
    INTEGER :: errno
    td_tmp => newtimedelta(this%td, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 133)
    NULLIFY(td_tmp2)
    td_tmp2           = td_tmp * lambda
    IF (ASSOCIATED(td_tmp2)) THEN
      CALL handle_errno(general_arithmetic_error, &
        & "externals/mtime/src/mtime_t_timedelta.inc", &
        & 139)
      RETURN
    END IF
    scaled_td%td      = td_tmp2
    scaled_td%td%sign = td_tmp2%sign
    IF (ASSOCIATED(td_tmp))  CALL deallocatetimedelta(td_tmp)
    IF (ASSOCIATED(td_tmp2)) CALL deallocatetimedelta(td_tmp2)
  END FUNCTION t_timedelta_scalar_multiply_long

  FUNCTION t_timedelta_scalar_multiply_inv_long(lambda, this) RESULT(scaled_td)
    TYPE(t_timedelta), TARGET :: scaled_td
    INTEGER(c_int64_t),       INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this
    TYPE(timedelta), POINTER             :: td_tmp, td_tmp2
    INTEGER :: errno
    td_tmp => newtimedelta(this%td, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 157)
    NULLIFY(td_tmp2)
    td_tmp2 = td_tmp * lambda
    IF (ASSOCIATED(td_tmp2)) THEN
      CALL handle_errno(general_arithmetic_error, &
        & "externals/mtime/src/mtime_t_timedelta.inc", &
        & 163)
      RETURN
    END IF
    scaled_td%td      = td_tmp2
    scaled_td%td%sign = td_tmp2%sign
    IF (ASSOCIATED(td_tmp))  CALL deallocatetimedelta(td_tmp)
    IF (ASSOCIATED(td_tmp2)) CALL deallocatetimedelta(td_tmp2)
  END FUNCTION t_timedelta_scalar_multiply_inv_long

  FUNCTION t_timedelta_scalar_multiply_int (this, lambda) RESULT(scaled_td) !OK-tested
    TYPE(t_timedelta), TARGET :: scaled_td
    INTEGER(c_int32_t),       INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this
    TYPE(timedelta), POINTER             :: td_tmp
    TYPE(c_ptr) :: c_pointer, dummy_ptr, c_ptr_result

    dummy_ptr = my_elementwisescalarmultiplytimedelta(c_loc(this%td), int(lambda, c_int64_t), c_loc(scaled_td%td))
    IF (.NOT. c_associated(dummy_ptr)) THEN
      CALL handle_errno(0*100+2, &
        & "externals/mtime/src/mtime_t_timedelta.inc", &
        & 183)
      scaled_td%td%sign = '!'
    ENDIF

  END FUNCTION t_timedelta_scalar_multiply_int

  FUNCTION t_timedelta_scalar_multiply_inv_int(lambda, this) RESULT(scaled_td)
    TYPE(t_timedelta), TARGET :: scaled_td
    INTEGER(c_int32_t),       INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this

    scaled_td = t_timedelta_scalar_multiply_int(this, lambda)
  END FUNCTION t_timedelta_scalar_multiply_inv_int

  FUNCTION t_timedelta_scalar_multiply_real (this, lambda) RESULT(scaled_td) !OK-tested
    TYPE(t_timedelta),  TARGET             :: scaled_td
    REAL(c_double),             INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this
    TYPE(c_ptr)                            :: dummy_ptr

    dummy_ptr = my_elementwisescalarmultiplytimedeltadp(c_loc(this%td), lambda, c_loc(scaled_td%td))
    IF (.NOT. c_associated(dummy_ptr)) THEN
      CALL handle_errno(0*100+2, &
        & "externals/mtime/src/mtime_t_timedelta.inc", &
        & 207)
      scaled_td%td%sign = '!'
    ENDIF

  END FUNCTION t_timedelta_scalar_multiply_real

  FUNCTION t_timedelta_scalar_multiply_inv_real(lambda, this) RESULT(scaled_td)
    TYPE(t_timedelta), TARGET :: scaled_td
    REAL(c_double),       INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this

    scaled_td = t_timedelta_scalar_multiply_real(this, lambda)
  END FUNCTION t_timedelta_scalar_multiply_inv_real

  ! Convert t_timedelta object to string.
  !
  FUNCTION t_timedelta_toString(this) result(string)
    CHARACTER(len=max_timedelta_str_len)  :: string
    CLASS (t_timedelta)                   :: this
    type(c_ptr) :: c_pointer, dummy_ptr
    integer :: i

    string = ""

    c_pointer = this%get_c_pointer()
    CALL handle_errno(.not. c_associated(c_pointer), 0 * 100 + 2, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 234)
    dummy_ptr = my_timedeltatostring(c_pointer, string)
    CALL handle_errno(.not. c_associated(dummy_ptr), 4 * 100 + 6, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 238)

    char_loop: do i = 1 , len(string)
       if (string(i:i) == c_null_char) exit char_loop
    end do char_loop
    string(i:len(string)) = ' '

    CALL my_deallocatetimedelta(c_pointer) !deallocate timedelta
  END FUNCTION t_timedelta_toString

  FUNCTION t_timedelta_divideInSecondsBy (this, divisor, referenceDateTime) RESULT(quotient)
    CLASS(t_timedelta), INTENT(in)           :: this
    TYPE(t_timedelta),  INTENT(in)           :: divisor
    TYPE(t_datetime),   INTENT(IN), OPTIONAL :: referenceDateTime
    TYPE(divisionquotienttimespan)           :: quotient
    TYPE(timedelta),                POINTER  :: tmp_dividend
    TYPE(timedelta),                POINTER  :: tmp_divisor
    TYPE(datetime),                 POINTER  :: tmp_ref, tmp_dt
    TYPE(t_datetime)                         :: dt_tmp
    INTEGER                                  :: errno

    tmp_dividend => newtimedelta(this%td, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 262)
    tmp_divisor => newtimedelta(divisor%td, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 266)
    IF (PRESENT(referenceDateTime)) THEN
      tmp_ref => newDateTime(referenceDateTime%dt, errno)
      CALL handle_errno(errno, &
        & "externals/mtime/src/mtime_t_timedelta.inc", &
        & 271)
    ENDIF

    CALL divideTimeDeltaInSeconds(tmp_dividend, tmp_divisor, quotient, errno)

    IF (errno /= no_error) THEN
      IF (.NOT. PRESENT(referenceDateTime)) THEN
        CALL handle_errno(errno, &
          & "externals/mtime/src/mtime_t_timedelta.inc", &
          & 280)
      ELSE
        dt_tmp = referenceDateTime + this
        tmp_dt => newDateTime(dt_tmp%dt, errno)
        CALL handle_errno(errno, &
          & "externals/mtime/src/mtime_t_timedelta.inc", &
          & 286)
        CALL divideDatetimeDifferenceInSeconds(tmp_dt, tmp_ref, &
             &                                tmp_divisor, quotient, errno)
        CALL handle_errno(errno, &
          & "externals/mtime/src/mtime_t_timedelta.inc", &
          & 291)

        CALL deallocateDatetime(tmp_dt)
      ENDIF
    END IF

    CALL deallocatetimedelta(tmp_dividend)
    CALL deallocatetimedelta(tmp_divisor)
    IF (PRESENT(referenceDateTime)) THEN
      CALL deallocateDatetime(tmp_ref)
    ENDIF
  END FUNCTION t_timedelta_divideInSecondsBy

  FUNCTION t_timedelta_toSeconds (this, td) RESULT(seconds)
    CLASS(t_timedelta), INTENT(in) :: this
    TYPE(t_datetime),   INTENT(in) :: td
    INTEGER(c_int64_t)             :: seconds

    seconds = getTotalSecondsTimeDelta(this%td, td%dt)
  END FUNCTION t_timedelta_toSeconds

  FUNCTION t_timedelta_toMilliSeconds (this, td) RESULT(ms)
    CLASS(t_timedelta), INTENT(in) :: this
    TYPE(t_datetime),   INTENT(in) :: td
    INTEGER(c_int64_t)             :: ms

    ms = getTotalMilliSecondsTimeDelta(this%td, td%dt)
  END FUNCTION t_timedelta_toMilliSeconds

  FUNCTION t_timedelta_toJulianDelta(this, dt) RESULT(jd)
    TYPE(t_juliandelta) :: jd
    CLASS(t_timedelta), INTENT(in) :: this
    TYPE(t_datetime),   INTENT(in) :: dt
    TYPE(timedelta),    POINTER    :: td_tmp
    TYPE(datetime),     POINTER    :: dt_tmp
    TYPE(juliandelta),  POINTER    :: jd_tmp
    INTEGER                        :: errno
    td_tmp => newtimedelta(this%td, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 331)
    dt_tmp => newDatetime(dt%dt, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 335)
    jd_tmp => newJuliandelta('+', 0_c_int64_T, 0_c_int64_T, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 339)

    CALL timeDeltaToJulianDelta(td_tmp,dt_tmp,jd_tmp)
    jd%jd = jd_tmp

    CALL deallocatetimedelta(td_tmp)
    CALL deallocateDatetime(dt_tmp)
    CALL deallocateJuliandelta(jd_tmp)
  END FUNCTION t_timedelta_toJulianDelta

  FUNCTION t_julianDay_getDay (this) RESULT(d)
    CLASS(t_julianday), INTENT(in) :: this
    INTEGER(c_int64_t) :: d
    d = this%jd%day
  END FUNCTION t_julianDay_getDay

  FUNCTION t_julianday_getFractionOfDayInMS (this) RESULT(ms)
    CLASS(t_julianday), INTENT(in) :: this
    INTEGER(c_int64_t) :: ms
    ms = this%jd%ms
  END FUNCTION t_julianday_getFractionOfDayInMS

  FUNCTION t_julianday_toDateTime(this) RESULT(res)
    TYPE(t_datetime) :: res
    CLASS(t_julianday)           :: this
    TYPE(julianday),  POINTER    :: jd_tmp
    INTEGER                      :: errno
    TYPE(datetime),   POINTER    :: dt_tmp

    jd_tmp => newJulianDay(this%jd%day, this%jd%ms, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 371)
    dt_tmp => newDatetime("1970-01-01T00:00:00", errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 375)

    CALL  getDatetimeFromJulianDay(jd_tmp, dt_tmp, errno)
    CALL handle_errno(errno, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 380)

    res%dt = dt_tmp
    CALL deallocateJulianDay(jd_tmp)
    CALL deallocateDatetime(dt_tmp)
  END FUNCTION t_julianday_toDateTime

  FUNCTION t_datetime_min(a,b) RESULT(res)
    TYPE(t_datetime) :: a,b
    TYPE(t_datetime) :: res

    IF (a > b) THEN
      res = b
    ELSE
      res = a
    END IF
  END FUNCTION t_datetime_min

  FUNCTION t_datetime_max(a,b) RESULT(res)
    TYPE(t_datetime) :: a,b
    TYPE(t_datetime) :: res

    IF (a > b) THEN
      res = a
    ELSE
      res = b
    END IF
  END FUNCTION t_datetime_max

  FUNCTION t_timedelta_get_c_pointer(this) RESULT(c_pointer)
    TYPE(c_ptr) :: c_pointer
    CLASS(t_timedelta) :: this
    character(c_char) ::c_sign

    c_sign = this%td%sign(1:1)
    c_pointer = my_newrawtimedelta(c_sign, int(this%td%year,c_int64_t), this%td%month, this%td%day, &
         &                         this%td%hour, this%td%minute, this%td%second, this%td%ms)
    call handle_errno((.not. c_associated(c_pointer)), 5 * 100 + 2, &
      & "externals/mtime/src/mtime_t_timedelta.inc", &
      & 419)
  END FUNCTION t_timedelta_get_c_pointer
!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
  ! ================================================================================
  ! juliandelta section:
  ! ================================================================================

  ! generic assignment for constructors
  !
  FUNCTION t_juliandelta_assign_raw(sign,day, ms)
    TYPE(t_juliandelta) :: t_juliandelta_assign_raw
    CHARACTER(c_char),  INTENT(in) :: sign
    INTEGER(c_int64_t), INTENT(in) :: day
    INTEGER(c_int64_t), INTENT(in) :: ms
    type(c_ptr) :: c_pointer
    TYPE(juliandelta), pointer :: jd_tmp

    c_pointer = my_newjuliandelta(sign, day, ms)
    !print *,sign, c_pointer
    if (.not. c_associated(c_pointer)) then
       call handle_errno(1 * 100 + 1, "externals/mtime/src/mtime_t_juliandelta.inc", 22)
       t_juliandelta_assign_raw%jd%sign = 'L'
    else
       call c_f_pointer(c_pointer, jd_tmp)
       t_juliandelta_assign_raw%jd = jd_tmp
       t_juliandelta_assign_raw%jd%sign = jd_tmp%sign
       call my_deallocatejuliandelta(c_pointer)
    end if
  END FUNCTION t_juliandelta_assign_raw
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


END MODULE mtime_hl
!>
!! @}
