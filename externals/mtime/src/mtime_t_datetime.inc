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
      & __FILE__, &
      & __LINE__)
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
      & __FILE__, &
      & __LINE__)
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
      & __FILE__, &
      & __LINE__)
    dummy_ptr = my_datetimetostring(c_pointer, string)
    CALL handle_errno(.not. c_associated(dummy_ptr), 4 * 100 + 6, &
      & __FILE__, &
      & __LINE__)

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
      & __FILE__, &
      & __LINE__)
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
      & __FILE__, &
      & __LINE__)
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
      & __FILE__, &
      & __LINE__)
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
      & __FILE__, &
      & __LINE__)
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
      & __FILE__, __LINE__)
  END FUNCTION t_datetime_daysInEntireMonth

  FUNCTION t_datetime_daysInEntireYear(this)
    CLASS (t_datetime),  INTENT(in), target :: this
    INTEGER(c_int) :: t_datetime_daysInEntireYear
    t_datetime_daysInEntireYear = my_getnoofdaysinyeardatetime(c_loc(this%dt))
    CALL handle_errno(t_datetime_daysInEntireYear == 0, 4 * 100 + 16, &
      & __FILE__, __LINE__)
  END FUNCTION t_datetime_daysInEntireYear

  FUNCTION t_datetime_elapsedDaysInYear(this)
    CLASS (t_datetime),  INTENT(in), target :: this
    INTEGER(c_int) :: t_datetime_elapsedDaysInYear

    t_datetime_elapsedDaysInYear = my_getdayofyearfromdatetime(c_loc(this%dt))
    CALL handle_errno(t_datetime_elapsedDaysInYear == 0, 4 * 100 + 17, &
         &            __FILE__, __LINE__)
  END FUNCTION t_datetime_elapsedDaysInYear

  FUNCTION t_datetime_elapsedSecondsInMonth(this)
    CLASS (t_datetime),  INTENT(in), target :: this
    INTEGER(c_int64_t) :: t_datetime_elapsedSecondsInMonth
    t_datetime_elapsedSecondsInMonth = my_getnoofsecondselapsedinmonthdatetime(c_loc(this%dt))
    CALL handle_errno(t_datetime_elapsedSecondsInMonth == -1, 4 * 100 + 18, &
      & __FILE__, __LINE__)
  END FUNCTION t_datetime_elapsedSecondsInMonth

  FUNCTION t_datetime_elapsedSecondsInDay(this)
    CLASS (t_datetime),  INTENT(in), target :: this
    INTEGER(c_int64_t) :: t_datetime_elapsedSecondsInDay

    t_datetime_elapsedSecondsInDay = my_getnoofsecondselapsedindaydatetime(c_loc(this%dt))
    CALL handle_errno(t_datetime_elapsedSecondsInDay == -1, 4 * 100 + 19, &
     & __FILE__, __LINE__)
  END FUNCTION t_datetime_elapsedSecondsInDay


  FUNCTION t_datetime_get_c_pointer(this) RESULT(c_pointer)
    TYPE(c_ptr) :: c_pointer
    CLASS(t_datetime) :: this
    c_pointer = my_newrawdatetime(int(this%dt%date%year,c_int64_t), this%dt%date%month,     &
         &                        this%dt%date%day, this%dt%time%hour, this%dt%time%minute, &
         &                        this%dt%time%second, this%dt%time%ms)
    call handle_errno((.not. c_associated(c_pointer)), 4 * 100 + 3, &
      & __FILE__, &
      & __LINE__)
  END FUNCTION t_datetime_get_c_pointer


