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
    c_pointer = my_newdatetime(TRIM(ADJUSTL(dt_string))//c_null_char)
    CALL handle_errno(.NOT. C_ASSOCIATED(c_pointer), 4*100 + 1, &
      & __FILE__, &
      & __LINE__)
    CALL C_F_POINTER(c_pointer, dt_tmp)
    t_datetime_assign_string%dt = dt_tmp
    CALL my_deallocatedatetime(c_pointer)
  END FUNCTION t_datetime_assign_string

  ! constructor for a datetime
  !
  FUNCTION t_datetime_assign_raw(year, month, day, hour, minute, second, ms) RESULT(res)
    TYPE(t_datetime) :: res

    INTEGER(c_int64_t), INTENT(in) :: year
    INTEGER(c_int), INTENT(in) :: month, day, hour, minute, second, ms

    TYPE(datetime), POINTER :: dt_tmp
    TYPE(c_ptr)             :: c_pointer

    c_pointer = my_newrawdatetime(year, month, day, hour, minute, second, ms)
    CALL handle_errno(.NOT. C_ASSOCIATED(c_pointer), 4*100 + 2, &
      & __FILE__, &
      & __LINE__)
    CALL C_F_POINTER(c_pointer, dt_tmp)
    res%dt = dt_tmp
    CALL my_deallocatedatetime(c_pointer)
  END FUNCTION t_datetime_assign_raw

  ! Returns t_datetime objects day
  !
  FUNCTION t_datetime_getDay(this)
    INTEGER            :: t_datetime_getDay
    CLASS(t_datetime) :: this
    t_datetime_getDay = this%dt%date%day
  END FUNCTION t_datetime_getDay

  ! Convert t_datetime object to string.
  !
  FUNCTION t_datetime_toString(this) RESULT(string)
    CHARACTER(len=max_datetime_str_len)  :: string
    CLASS(t_datetime)                   :: this
    TYPE(c_ptr) :: c_pointer, dummy_ptr
    INTEGER :: i

    string = ""

    c_pointer = this%get_c_pointer()
    CALL handle_errno(.NOT. C_ASSOCIATED(c_pointer), 0*100 + 2, &
      & __FILE__, &
      & __LINE__)
    dummy_ptr = my_datetimetostring(c_pointer, string)
    CALL handle_errno(.NOT. C_ASSOCIATED(dummy_ptr), 4*100 + 6, &
      & __FILE__, &
      & __LINE__)

    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '

    CALL my_deallocatedatetime(c_pointer)
  END FUNCTION t_datetime_toString

  ! Convert t_datetime object to string.
  !
  FUNCTION t_datetime_to_posix_string(this, format_string) RESULT(string)
    CHARACTER(len=max_datetime_str_len)  :: string
    CHARACTER(len=*), INTENT(in)         :: format_string
    CLASS(t_datetime)                   :: this
    INTEGER :: i
    TYPE(c_ptr) :: c_pointer, dummy_ptr

    string = ""
    c_pointer = this%get_c_pointer()
    dummy_ptr = my_datetoposixstring(c_pointer, string, format_string)
    CALL handle_errno(.NOT. C_ASSOCIATED(dummy_ptr), 2*100 + 8, &
      & __FILE__, &
      & __LINE__)
    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '

    CALL my_deallocatedatetime(c_pointer)
  END FUNCTION t_datetime_to_posix_string

  FUNCTION t_datetime_toJulianDay(this) RESULT(jd)
    CLASS(t_datetime), INTENT(in) :: this
    TYPE(t_julianday), TARGET :: jd
    TYPE(c_ptr) :: c_pointer, dummy_ptr
    c_pointer = this%get_c_pointer()
    jd = t_julianday(julianday(0, 0))
    dummy_ptr = my_getjuliandayfromdatetime(c_pointer, C_LOC(jd%jd))
    CALL handle_errno(.NOT. C_ASSOCIATED(dummy_ptr), 0*100 + 2, &
      & __FILE__, &
      & __LINE__)
    CALL my_deallocatedatetime(c_pointer)
  END FUNCTION t_datetime_toJulianDay

  ! Addition of time interval to datetime object.
  !
  FUNCTION t_datetime_add_timedelta(this, td) RESULT(dt_td_sum)
    TYPE(t_datetime)               :: dt_td_sum
    CLASS(t_datetime), INTENT(in) :: this
    CLASS(t_timedelta), INTENT(in) :: td
    TYPE(datetime), POINTER        :: dt_tmp

    TYPE(c_ptr) :: c_pointer1, c_pointer2, dummy_ptr
    c_pointer1 = this%get_c_pointer()
    c_pointer2 = td%get_c_pointer()

    dummy_ptr = my_addtimedeltatodatetime(c_pointer1, c_pointer2, c_pointer1)
    CALL handle_errno(.NOT. C_ASSOCIATED(dummy_ptr), 0*100 + 2, &
      & __FILE__, &
      & __LINE__)
    CALL C_F_POINTER(c_pointer1, dt_tmp)
    dt_td_sum%dt = dt_tmp
    CALL my_deallocatedatetime(c_pointer1)
    CALL my_deallocatedatetime(c_pointer2)
  END FUNCTION t_datetime_add_timedelta

  ! Subtraction of time interval to datetime object.
  !
  FUNCTION t_datetime_sub_timedelta(this, td) RESULT(dt_td_sum)
    TYPE(t_datetime)               :: dt_td_sum
    CLASS(t_datetime), INTENT(in) :: this
    TYPE(t_timedelta), INTENT(in) :: td
    TYPE(t_timedelta)              :: td_tmp
    TYPE(datetime), POINTER        :: dt_tmp
    TYPE(c_ptr) :: c_pointer1, c_pointer2, dummy_ptr

    td_tmp = td
    IF (td_tmp%td%sign == "+") THEN
      td_tmp%td%sign = "-"
    ELSE
      td_tmp%td%sign = "+"
    END IF

    c_pointer1 = this%get_c_pointer()
    c_pointer2 = td_tmp%get_c_pointer()

    dummy_ptr = my_addtimedeltatodatetime(c_pointer1, c_pointer2, c_pointer1)
    CALL handle_errno(.NOT. C_ASSOCIATED(dummy_ptr), 0*100 + 2, &
      & __FILE__, &
      & __LINE__)
    CALL C_F_POINTER(c_pointer1, dt_tmp)
    dt_td_sum%dt = dt_tmp
    CALL my_deallocatedatetime(c_pointer1)
    CALL my_deallocatedatetime(c_pointer2)
  END FUNCTION t_datetime_sub_timedelta

  ! Subtraction of two dates.
  !
  FUNCTION t_datetime_sub_datetime(this, dt) RESULT(dt_dt_diff)
    TYPE(t_timedelta), TARGET :: dt_dt_diff
    CLASS(t_datetime), INTENT(in), TARGET :: this
    CLASS(t_datetime), INTENT(in), TARGET :: dt
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_gettimedeltafromdate(C_LOC(this%dt), C_LOC(dt%dt), C_LOC(dt_dt_diff%td))
  END FUNCTION t_datetime_sub_datetime

  ! Overloaded operator: test for equivalence.
  !
#ifndef MTIME_PURE_IF_C_LOC_IS_PURE
#  if defined(__NEC__) || (defined(NAGFOR) &&  __NAG_COMPILER_RELEASE <= 71)
!    NEC and older versions of NAG do not consider C_LOC as PURE
#    define MTIME_PURE_IF_C_LOC_IS_PURE
#  else
#    define MTIME_PURE_IF_C_LOC_IS_PURE PURE
#  endif
#endif
  MTIME_PURE_IF_C_LOC_IS_PURE LOGICAL FUNCTION t_datetime_equal(this, dt) RESULT(eq)
    CLASS(t_datetime), INTENT(in), TARGET :: this
    CLASS(t_datetime), INTENT(in), TARGET :: dt
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(C_LOC(this%dt), C_LOC(dt%dt))
    IF (ret == 0) THEN
      eq = .TRUE.
    ELSE
      eq = .FALSE.
    END IF
  END FUNCTION t_datetime_equal

  MTIME_PURE_IF_C_LOC_IS_PURE LOGICAL FUNCTION t_datetime_not_equal(this, dt)
    CLASS(t_datetime), INTENT(in) :: this
    CLASS(t_datetime), INTENT(in) :: dt
    t_datetime_not_equal = .NOT. (this%dt == dt%dt)
  END FUNCTION t_datetime_not_equal

  MTIME_PURE_IF_C_LOC_IS_PURE LOGICAL FUNCTION t_datetime_less_than(this, dt) RESULT(lt)
    CLASS(t_datetime), INTENT(in), TARGET :: this
    CLASS(t_datetime), INTENT(in), TARGET :: dt
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(C_LOC(this%dt), C_LOC(dt%dt))
    IF (ret == -1) THEN
      lt = .TRUE.
    ELSE
      lt = .FALSE.
    END IF
  END FUNCTION t_datetime_less_than

  MTIME_PURE_IF_C_LOC_IS_PURE LOGICAL FUNCTION t_datetime_greater_than(this, dt) RESULT(gt)
    CLASS(t_datetime), INTENT(in), TARGET :: this
    CLASS(t_datetime), INTENT(in), TARGET :: dt
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(C_LOC(this%dt), C_LOC(dt%dt))
    IF (ret == 1) THEN
      gt = .TRUE.
    ELSE
      gt = .FALSE.
    END IF
  END FUNCTION t_datetime_greater_than

  MTIME_PURE_IF_C_LOC_IS_PURE LOGICAL FUNCTION t_datetime_less_or_equal(this, dt)
    CLASS(t_datetime), INTENT(in) :: this
    CLASS(t_datetime), INTENT(in) :: dt
    t_datetime_less_or_equal = .NOT. (this > dt)
  END FUNCTION t_datetime_less_or_equal

  MTIME_PURE_IF_C_LOC_IS_PURE LOGICAL FUNCTION t_datetime_greater_or_equal(this, dt)
    CLASS(t_datetime), INTENT(in) :: this
    CLASS(t_datetime), INTENT(in) :: dt
    t_datetime_greater_or_equal = .NOT. (this < dt)
  END FUNCTION t_datetime_greater_or_equal

  FUNCTION t_datetime_daysInEntireMonth(this)
    CLASS(t_datetime), INTENT(in), TARGET :: this
    INTEGER(c_int) :: t_datetime_daysInEntireMonth
    t_datetime_daysInEntireMonth = my_getnoofdaysinmonthdatetime(C_LOC(this%dt))
    CALL handle_errno(t_datetime_daysInEntireMonth == 0, 4*100 + 15, &
      & __FILE__, __LINE__)
  END FUNCTION t_datetime_daysInEntireMonth

  FUNCTION t_datetime_daysInEntireYear(this)
    CLASS(t_datetime), INTENT(in), TARGET :: this
    INTEGER(c_int) :: t_datetime_daysInEntireYear
    t_datetime_daysInEntireYear = my_getnoofdaysinyeardatetime(C_LOC(this%dt))
    CALL handle_errno(t_datetime_daysInEntireYear == 0, 4*100 + 16, &
      & __FILE__, __LINE__)
  END FUNCTION t_datetime_daysInEntireYear

  FUNCTION t_datetime_elapsedDaysInYear(this)
    CLASS(t_datetime), INTENT(in), TARGET :: this
    INTEGER(c_int) :: t_datetime_elapsedDaysInYear

    t_datetime_elapsedDaysInYear = my_getdayofyearfromdatetime(C_LOC(this%dt))
    CALL handle_errno(t_datetime_elapsedDaysInYear == 0, 4*100 + 17, &
         &            __FILE__, __LINE__)
  END FUNCTION t_datetime_elapsedDaysInYear

  FUNCTION t_datetime_elapsedSecondsInMonth(this)
    CLASS(t_datetime), INTENT(in), TARGET :: this
    INTEGER(c_int64_t) :: t_datetime_elapsedSecondsInMonth
    t_datetime_elapsedSecondsInMonth = my_getnoofsecondselapsedinmonthdatetime(C_LOC(this%dt))
    CALL handle_errno(t_datetime_elapsedSecondsInMonth == -1, 4*100 + 18, &
      & __FILE__, __LINE__)
  END FUNCTION t_datetime_elapsedSecondsInMonth

  FUNCTION t_datetime_elapsedSecondsInDay(this)
    CLASS(t_datetime), INTENT(in), TARGET :: this
    INTEGER(c_int64_t) :: t_datetime_elapsedSecondsInDay

    t_datetime_elapsedSecondsInDay = my_getnoofsecondselapsedindaydatetime(C_LOC(this%dt))
    CALL handle_errno(t_datetime_elapsedSecondsInDay == -1, 4*100 + 19, &
     & __FILE__, __LINE__)
  END FUNCTION t_datetime_elapsedSecondsInDay

  FUNCTION t_datetime_get_c_pointer(this) RESULT(c_pointer)
    TYPE(c_ptr) :: c_pointer
    CLASS(t_datetime) :: this
    c_pointer = my_newrawdatetime(INT(this%dt%date%year, c_int64_t), this%dt%date%month,     &
         &                        this%dt%date%day, this%dt%time%hour, this%dt%time%minute, &
         &                        this%dt%time%second, this%dt%time%ms)
    CALL handle_errno((.NOT. C_ASSOCIATED(c_pointer)), 4*100 + 3, &
      & __FILE__, &
      & __LINE__)
  END FUNCTION t_datetime_get_c_pointer

