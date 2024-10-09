!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
!> @file libmtime.f90
!!
!! @brief Providing the Fortran language bindings for libmtime
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!! @defgroup FortranBindings libmtime Fortran language bindings
!! @{
!!
!___________________________________________________________________________________________________________
!>
!! @brief Provides the calendar to the users, abstracting the different calendars available.
!!
!! @details
!!
!! Three calendar types are provided:
!!
!!   - a proleptic Gregorian calendar
!!   - a calendar with 365 days per year without leap years
!!   - a calendar with 360 days per year and each month having 30 days
!!
!! To use a specific calendar a call to setCalendar with the
!! respective selector must be done. The implementation is based on a
!! singleton concept meaning that only one calendar can be active at a
!! time. To release a calendar a call to resetCalendar has to be done.
!!
!___________________________________________________________________________________________________________
MODULE mtime_calendar
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_associated, c_null_char, c_ptr
  USE mtime_constants, ONLY: max_calendar_str_len
  USE mtime_c_bindings, ONLY: my_calendartostring
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  ENUM, BIND(c)
    ENUMERATOR :: calendar_not_set = 0   ! calendar is not defined yet
    ENUMERATOR :: proleptic_gregorian = 1   ! proleptic Gregorian calendar
    ENUMERATOR :: year_of_365_days = 2   ! 365 day year without leap years
    ENUMERATOR :: year_of_360_days = 3   ! 360 day year with 30 day months
  END ENUM
  !
  PUBLIC :: calendar_not_set, proleptic_gregorian, year_of_365_days, year_of_360_days
  PUBLIC :: calendarToString
  !
CONTAINS
  !
#ifdef DOXYGEN_DOCUMENTATION_ONLY
  !>
  !! @brief Initialize a new calendar.
  !!
  !! setCalendar is done at the very begining to select one of the
  !! provided calendar libraries. It intializes the calendar to one of:
  !!
  !! - proleptic_gregorian
  !! - year_of_365_days
  !! - year_of_360_days
  !!
  !! The calendar type and hence it's behaviour (Calendar to Julian
  !! conversion and vice versa) is fixed for the lifetime of the
  !! selected calendar.  Attempts to change the calendar type on the
  !! fly is discouraged. The lib has built-in checks to reject
  !! change attempts at run time.  However, a calendar can be
  !! "re-initialized" after calling resetCalendar(), but this is not
  !! advised.
  !!
  !! MANTRA: Know what you are doing before you do it and do it
  !! right the first time.
  !!
#endif
  !>
  !! @brief convert the calendar identifier into a human readable string
  !!
  !! @param[out]       string      the calendar type verbose
  !! @param[out]       errno       optional, error message
  !!
  RECURSIVE SUBROUTINE calendarToString(string, errno) !TESTED-OK
    CHARACTER(len=max_calendar_str_len), INTENT(out) :: string
    INTEGER :: i
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_calendartostring(string)
    IF (.NOT. C_ASSOCIATED(c_pointer)) THEN
      string = '<null>'
    ELSE
      char_loop: DO i = 1, LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    END IF
    IF (PRESENT(errno)) errno = MERGE(0, 0*100 + 4, C_ASSOCIATED(c_pointer))
  END SUBROUTINE calendarToString
  !
END MODULE mtime_calendar

MODULE mtime_juliandelta
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_associated, c_char, c_f_pointer, c_int64_t, c_loc, c_ptr
  USE mtime_c_bindings, ONLY: &
    & juliandelta, my_deallocatejuliandelta, my_newjuliandelta
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: juliandelta
  PUBLIC :: newJulianDelta
  PUBLIC :: deallocateJulianDelta
  !
CONTAINS
  !
  RECURSIVE FUNCTION newJuliandelta(sign, day, ms, errno) RESULT(ret_juliandelta) !OK-TESTED.
    TYPE(juliandelta), POINTER :: ret_juliandelta
    CHARACTER(c_char), VALUE, INTENT(in) :: sign
    INTEGER(c_int64_t), INTENT(in) :: day
    INTEGER(c_int64_t), INTENT(in) :: ms
    TYPE(c_ptr)                    :: c_pointer
    INTEGER, OPTIONAL              :: errno

    c_pointer = my_newjuliandelta(sign, day, ms)
    IF (PRESENT(errno)) errno = MERGE(0, 1*100 + 1, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_juliandelta)
  END FUNCTION newJuliandelta
  !
  !>
  !! @brief destructor for a Julian delta
  RECURSIVE SUBROUTINE deallocateJuliandelta(my_juliandelta) !OK-TESTED.
    TYPE(juliandelta), POINTER :: my_juliandelta
    CALL my_deallocatejuliandelta(C_LOC(my_juliandelta))
    my_juliandelta => NULL()
  END SUBROUTINE deallocateJuliandelta
  !
END MODULE mtime_juliandelta
!>
!! @brief Julian Day Calendar and some operations supported on julian dates.
!!
!! @details
!___________________________________________________________________________________________________________
MODULE mtime_julianday
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_associated, c_f_pointer, c_int, c_int64_t, c_loc, c_null_char, c_ptr
  USE mtime_constants, ONLY: max_julianday_str_len
  USE mtime_c_bindings, ONLY: &
    & julianday, my_addjuliandeltatojulianday, my_comparejulianday, &
    & my_deallocatejulianday, my_juliandaytostring, my_newjulianday, &
    & my_replacejulianday, my_subtractjulianday
  USE mtime_juliandelta
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: julianday
  PUBLIC :: newJulianday
  PUBLIC :: deallocateJulianday
  PUBLIC :: juliandayToString
  PUBLIC :: ASSIGNMENT(=)
  PUBLIC :: OPERATOR(+)
  PUBLIC :: OPERATOR(-)
  PUBLIC :: OPERATOR(>)
  PUBLIC :: OPERATOR(<)
  PUBLIC :: OPERATOR(<=)
  PUBLIC :: OPERATOR(>=)
  PUBLIC :: OPERATOR(==)
  PUBLIC :: OPERATOR(/=)
  !
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE replacejulianday
  END INTERFACE ASSIGNMENT(=)
  !
  INTERFACE OPERATOR(+)
    MODULE PROCEDURE addjuliandeltatojulianday
    MODULE PROCEDURE addjuliandaytojuliandelta
  END INTERFACE OPERATOR(+)
  !
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE subtractjuliandayfromjulianday
  END INTERFACE OPERATOR(-)
  !
  INTERFACE OPERATOR(>)
    MODULE PROCEDURE julianday_gt
  END INTERFACE OPERATOR(>)
  !
  INTERFACE OPERATOR(<)
    MODULE PROCEDURE julianday_lt
  END INTERFACE OPERATOR(<)
  !
  INTERFACE OPERATOR(<=)
    MODULE PROCEDURE julianday_lt_or_eq
  END INTERFACE OPERATOR(<=)
  !
  INTERFACE OPERATOR(>=)
    MODULE PROCEDURE julianday_gt_or_eq
  END INTERFACE OPERATOR(>=)
  !
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE julianday_eq
  END INTERFACE OPERATOR(==)
  !
  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE julianday_ne
  END INTERFACE OPERATOR(/=)
  !
CONTAINS
  !
  !>
  !! @brief construct a new Julian date
  !!
  !! @param[in] day            the Julian day
  !! @param[in] ms             an integer denoting the actual milli seconds of a day
  !! @param[out]       errno       optional, error message
  !! @return    ret_julianday  a pointer of type(julianday)
  RECURSIVE FUNCTION newJulianday(day, ms, errno) RESULT(ret_julianday) !OK-TESTED.
    TYPE(julianday), POINTER :: ret_julianday
    INTEGER(c_int64_t), INTENT(in) :: day
    INTEGER(c_int64_t), INTENT(in) :: ms
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_newjulianday(day, ms)
    IF (PRESENT(errno)) errno = MERGE(0, 1*100 + 1, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_julianday)
  END FUNCTION newJulianday
  !
  !>
  !! @brief destructor for a Julian date
  !!
  !! @param     my_julianday   a pointer of type(julianday)
  RECURSIVE SUBROUTINE deallocateJulianday(my_julianday) !OK-TESTED.
    TYPE(julianday), POINTER :: my_julianday
    CALL my_deallocatejulianday(C_LOC(my_julianday))
    my_julianday => NULL()
  END SUBROUTINE deallocateJulianday
  !
  ! NOTE: Do not call the function using replacejulianday(.,.) directly; Use overloaded '=' instead as in "dest = src"
  RECURSIVE SUBROUTINE replacejulianday(dest, src) !OK-TESTED.
    TYPE(julianday), TARGET, INTENT(inout) :: dest
    TYPE(julianday), TARGET, INTENT(in) :: src
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_replacejulianday(C_LOC(src), C_LOC(dest))
  END SUBROUTINE replacejulianday
  !
  !>
  !! @brief get Julian day as a string.
  !!
  !! @param[in]  my_julianday   a pointer to type(julianday). The Julian day to be converted to a string
  !! @param[out] string         the Julian day verbose
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE juliandayToString(my_julianday, string, errno) !OK-TESTED.
    TYPE(julianday), POINTER :: my_julianday
    CHARACTER(len=max_julianday_str_len), INTENT(out) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_juliandaytostring(C_LOC(my_julianday), string)
    IF (PRESENT(errno)) errno = MERGE(0, 1*100 + 3, C_ASSOCIATED(dummy_ptr))
    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE juliandayToString
  !
  RECURSIVE FUNCTION addJuliandeltaToJulianday(op1, op2) RESULT(ret)
    TYPE(julianday), TARGET :: ret
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(juliandelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    ret = julianday(0, 0)
    dummy_ptr = my_addjuliandeltatojulianday(C_LOC(op1), C_LOC(op2), C_LOC(ret))
  END FUNCTION addJuliandeltaToJulianday
  !
  RECURSIVE FUNCTION addJuliandayToJuliandelta(op2, op1) RESULT(ret)
    TYPE(julianday), TARGET :: ret
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(juliandelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    ret = julianday(0, 0)
    dummy_ptr = my_addjuliandeltatojulianday(C_LOC(op1), C_LOC(op2), C_LOC(ret))
  END FUNCTION addJuliandayToJuliandelta
  !
  RECURSIVE FUNCTION subtractJuliandayFromJulianday(op1, op2) RESULT(ret)
    TYPE(juliandelta), TARGET :: ret
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    ret = juliandelta("+", 0, 0)
    dummy_ptr = my_subtractjulianday(C_LOC(op1), C_LOC(op2), C_LOC(ret))
  END FUNCTION subtractJuliandayFromJulianday
  !
#ifndef MTIME_PURE_IF_C_LOC_IS_PURE
#  if defined(__NEC__) || (defined(NAGFOR) &&  __NAG_COMPILER_RELEASE <= 71)
!    NEC and older versions of NAG do not consider C_LOC as PURE
#    define MTIME_PURE_IF_C_LOC_IS_PURE
#  else
#    define MTIME_PURE_IF_C_LOC_IS_PURE PURE
#  endif
#endif
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION julianday_gt(op1, op2) RESULT(gt)
    LOGICAL :: gt
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(C_LOC(op1), C_LOC(op2))
    IF (ret == 1) THEN
      gt = .TRUE.
    ELSE
      gt = .FALSE.
    END IF
  END FUNCTION julianday_gt
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION julianday_lt(op1, op2) RESULT(lt)
    LOGICAL :: lt
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(C_LOC(op1), C_LOC(op2))
    IF (ret == -1) THEN
      lt = .TRUE.
    ELSE
      lt = .FALSE.
    END IF
  END FUNCTION julianday_lt
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION julianday_lt_or_eq(op1, op2) RESULT(lt_or_eq)
    LOGICAL :: lt_or_eq
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(C_LOC(op1), C_LOC(op2))
    IF ((ret == 0) .OR. (ret == -1)) THEN
      lt_or_eq = .TRUE.
    ELSE
      lt_or_eq = .FALSE.
    END IF
  END FUNCTION julianday_lt_or_eq
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION julianday_gt_or_eq(op1, op2) RESULT(gt_or_eq)
    LOGICAL :: gt_or_eq
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(C_LOC(op1), C_LOC(op2))
    IF ((ret == 0) .OR. (ret == 1)) THEN
      gt_or_eq = .TRUE.
    ELSE
      gt_or_eq = .FALSE.
    END IF
  END FUNCTION julianday_gt_or_eq
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION julianday_eq(op1, op2) RESULT(eq)
    LOGICAL :: eq
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(C_LOC(op1), C_LOC(op2))
    IF (ret == 0) THEN
      eq = .TRUE.
    ELSE
      eq = .FALSE.
    END IF
  END FUNCTION julianday_eq
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION julianday_ne(op1, op2) RESULT(ne)
    LOGICAL :: ne
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(C_LOC(op1), C_LOC(op2))
    IF (ret /= 0) THEN
      ne = .TRUE.
    ELSE
      ne = .FALSE.
    END IF
  END FUNCTION julianday_ne
  !
END MODULE mtime_julianday
!>
!! @brief Date and some operations supported on Date.
!!
!! @details
!!
!___________________________________________________________________________________________________________
MODULE mtime_date
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_associated, c_f_pointer, c_int, c_int64_t, c_loc, c_null_char, c_ptr
  USE mtime_constants, ONLY: max_date_str_len
  USE mtime_c_bindings, ONLY: &
    & date, my_constructandcopydate, my_datetoposixstring, my_datetostring, &
    & my_deallocatedate, my_newdatefromstring, my_newrawdate, my_replacedate
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: date
  PUBLIC :: newdate
  PUBLIC :: deallocateDate
  PUBLIC :: replaceDate
  PUBLIC :: dateToString
  PUBLIC :: dateToPosixString
  !
  INTERFACE newdate
    MODULE PROCEDURE newdatefromstring
    MODULE PROCEDURE newdatefromraw
    MODULE PROCEDURE newdatefromraw_yi8
    MODULE PROCEDURE newdatefromconstructandcopy
  END INTERFACE newdate
  !
CONTAINS
  !
  !>
  !! @brief construct a new date
  !!
  !! @param[in] string         an ISO 8601 conforming date string
  !! @param[out]       errno       optional, error message
  !! @return    ret_date       a pointer of type(date)
  RECURSIVE FUNCTION newdatefromstring(string, errno) RESULT(ret_date)  !OK-TESTED.
    TYPE(date), POINTER :: ret_date
    CHARACTER(len=*), INTENT(in) :: string
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL :: errno
    c_pointer = my_newdatefromstring(TRIM(ADJUSTL(string))//c_null_char)
    IF (PRESENT(errno)) errno = MERGE(0, 2*100 + 1, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_date)
  END FUNCTION newdatefromstring
  !>
  !! @brief construct a new date from raw date
  !!
  !! @param[in] year           the year
  !! @param[in] month          the month
  !! @param[in] day            the day
  !! @param[out]       errno       optional, error message
  !! @return    ret_date       a pointer of type(date)
  RECURSIVE FUNCTION newdatefromraw_yi8(year, month, day, errno) RESULT(ret_date) !OK-TESTED.
    TYPE(date), POINTER :: ret_date
    INTEGER(c_int64_t), INTENT(in) :: year
    INTEGER(c_int), INTENT(in) :: month, day
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_newrawdate(year, month, day)
    IF (PRESENT(errno)) errno = MERGE(0, 2*100 + 2, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_date)
  END FUNCTION newdatefromraw_yi8
  !>
  !! @brief construct a new date from raw date
  !!
  !! @param[in] year           the year
  !! @param[in] month          the month
  !! @param[in] day            the day
  !! @param[out]       errno       optional, error message
  !! @return    ret_date       a pointer of type(date)
  RECURSIVE FUNCTION newdatefromraw(year, month, day, errno) RESULT(ret_date) !OK-TESTED.
    TYPE(date), POINTER :: ret_date
    INTEGER(c_int), INTENT(in) :: year, month, day
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_newrawdate(INT(year, c_int64_t), month, day)
    IF (PRESENT(errno)) errno = MERGE(0, 2*100 + 3, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_date)
  END FUNCTION newdatefromraw
  !>
  !! @brief construct a new date from an existing by construct and copy
  !!
  !! @param[in] src            a pointer of type(date)
  !! @param[out]       errno       optional, error message
  !! @return    ret_date       a pointer of type(date)
  RECURSIVE FUNCTION newdatefromconstructandcopy(src, errno) RESULT(dest) !OK-TESTED
    TYPE(date), POINTER :: dest
    TYPE(date), TARGET :: src
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_constructandcopydate(C_LOC(src))
    IF (PRESENT(errno)) errno = MERGE(0, 2*100 + 4, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, dest)
  END FUNCTION newdatefromconstructandcopy
  !>
  !! @brief destructor for a date
  !!
  !! @param[in] my_date        a pointer of type(date)
  RECURSIVE SUBROUTINE deallocateDate(my_date) !OK-TESTED.
    TYPE(date), POINTER :: my_date
    CALL my_deallocatedate(C_LOC(my_date))
    my_date => NULL()
  END SUBROUTINE deallocateDate
  !>
  !! @brief repace an existing date by a given one
  !!
  !! @param[in]  src            a pointer of type(date)
  !! @param[out] dest           a pointer of type(date)
  !! @param[out] errno          optional, error message
  RECURSIVE SUBROUTINE replaceDate(dest, src, errno)  !OK-TESTED.
    TYPE(date), TARGET, INTENT(inout) :: dest
    TYPE(date), TARGET, INTENT(in) :: src
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_replacedate(C_LOC(src), C_LOC(dest))
    IF (PRESENT(errno)) errno = MERGE(0, 2*100 + 6, C_ASSOCIATED(dummy_ptr))
  END SUBROUTINE replaceDate
  !>
  !! @brief Get Date as an extended string.
  !!
  !! DateToString returns a string in IS08601 compliant (and extended) format.
  !!
  !! @param[in]   my_date
  !!         A pointer to type date. The date to be converted to string.
  !!
  !! @param[out]  string
  !!         String where date is to be written.
  !!
  !! @param[out]  errno
  !!         Optional, error message
  RECURSIVE SUBROUTINE dateToString(my_date, string, errno) !OK-TESTED.
    TYPE(date), POINTER :: my_date
    CHARACTER(len=max_date_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_datetostring(C_LOC(my_date), string)
    IF (PRESENT(errno)) errno = MERGE(0, 2*100 + 7, C_ASSOCIATED(dummy_ptr))
    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE dateToString
  !>
  !! @brief Get date and return as a string.
  !!
  !! Only dates between and including 1582-10-15 TO 9999-12-31 supported.
  !!
  !! @param[in]  my_date
  !!         A pointer to type date. The date to be converted to string.
  !!
  !! @param[out]  string
  !!         String where date is to be written.
  !!
  !! @param[in]  fmtstr
  !!         Desired Format string. CRITICAL: Inappropriate fmt string will cause dump.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE dateToPosixString(my_date, string, fmtstr, errno) !OK-TESTED.
    TYPE(date), POINTER :: my_date
    CHARACTER(len=max_date_str_len) :: string
    CHARACTER(len=*) :: fmtstr
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_datetoposixstring(C_LOC(my_date), string, fmtstr)
    IF (PRESENT(errno)) errno = MERGE(0, 2*100 + 8, C_ASSOCIATED(dummy_ptr))
    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE dateToPosixString
  !
END MODULE mtime_date
!>
!! @brief Time and some operations supported on Time
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
MODULE mtime_time
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_associated, c_f_pointer, c_int, c_loc, c_null_char, c_ptr
  USE mtime_constants, ONLY: max_time_str_len
  USE mtime_c_bindings, ONLY: &
    & my_constructandcopytime, my_deallocatetime, my_newrawtime, &
    & my_newtimefromstring, my_replacetime, my_timetoposixstring, &
    & my_timetostring, time
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: time
  PUBLIC :: newtime
  PUBLIC :: deallocateTime
  PUBLIC :: replaceTime
  PUBLIC :: timeToString
  PUBLIC :: timeToPosixString
  !
  INTERFACE newtime
    MODULE PROCEDURE newtimefromstring
    MODULE PROCEDURE newtimefromraw
    MODULE PROCEDURE newtimefromconstructandcopy
  END INTERFACE newtime
  !
CONTAINS
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newtimefromstring(string, errno) RESULT(ret_time) !OK-TESTED.
    TYPE(time), POINTER :: ret_time
    CHARACTER(len=*), INTENT(in) :: string
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_newtimefromstring(TRIM(ADJUSTL(string))//c_null_char)
    IF (PRESENT(errno)) errno = MERGE(0, 3*100 + 1, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_time)
  END FUNCTION newtimefromstring
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newtimefromraw(hour, minute, second, ms, errno) RESULT(ret_time) !OK-TESTED.
    TYPE(time), POINTER :: ret_time
    INTEGER(c_int), INTENT(in) :: hour, minute, second, ms
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_newrawtime(hour, minute, second, ms)
    IF (PRESENT(errno)) errno = MERGE(0, 3*100 + 2, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_time)
  END FUNCTION newtimefromraw
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newtimefromconstructandcopy(src, errno) RESULT(dest) !OK-TESTED.
    TYPE(time), POINTER :: dest
    TYPE(time), TARGET :: src
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_constructandcopytime(C_LOC(src))
    IF (PRESENT(errno)) errno = MERGE(0, 3*100 + 3, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, dest)
  END FUNCTION newtimefromconstructandcopy
  !>
  !! @brief Destructor of Time.
  !!
  !! @param  my_time
  !!         A pointer to type time. my_time is deallocated.
  RECURSIVE SUBROUTINE deallocateTime(my_time) !OK-TESTED.
    TYPE(time), POINTER :: my_time
    CALL my_deallocatetime(C_LOC(my_time))
    my_time => NULL()
  END SUBROUTINE deallocateTime
  !>
  !! @brief COPY a time object.
  !!
  !! Routine replaceTime copies the contents of source Time into a Destination Time object.
  !!
  !! @param[in]  src
  !!         A pointer to type time. Copy "FROM" time object.
  !!
  !! @param[out]  dest
  !!      A pointer to type time. Copy "TO" time object.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE replaceTime(dest, src, errno) !OK-TESTED.
    TYPE(time), TARGET, INTENT(in) :: src
    TYPE(time), TARGET, INTENT(inout) :: dest
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_replacetime(C_LOC(src), C_LOC(dest))
    IF (PRESENT(errno)) errno = MERGE(0, 3*100 + 5, C_ASSOCIATED(dummy_ptr))
  END SUBROUTINE replaceTime
  !>
  !! @brief Get time as an extended string.
  !!
  !! timetoString returns a string in IS08601 compliant (and extended) format.
  !!
  !! @param[in]  my_time
  !!         A pointer to type time. The time to be converted to string.
  !!
  !! @param[out]  string
  !!         String where time is to be written.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE timeToString(my_time, string, errno) !OK-TESTED.
    TYPE(time), POINTER :: my_time
    CHARACTER(len=max_time_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_timetostring(C_LOC(my_time), string)
    IF (PRESENT(errno)) errno = MERGE(0, 3*100 + 6, C_ASSOCIATED(dummy_ptr))
    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE timeToString
  !>
  !! @brief Get time as a Posix formated string.
  !!
  !! @param[in]  my_time
  !!         A pointer to type time. The time to be converted to string.
  !!
  !! @param[out]  string
  !!         String where time is to be written.
  !!
  !! @param[in]  fmtstr
  !!         Desired Format string. CRITICAL: Inappropriate fmt string will cause dump.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE timeToPosixString(my_time, string, fmtstr, errno) !OK-TESTED.
    TYPE(time), POINTER :: my_time
    CHARACTER(len=max_time_str_len) :: string
    CHARACTER(len=32) :: fmtstr
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_timetoposixstring(C_LOC(my_time), string, fmtstr)
    IF (PRESENT(errno)) errno = MERGE(0, 3*100 + 7, C_ASSOCIATED(dummy_ptr))
    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE timeToPosixString
  !
END MODULE mtime_time
!>
!! @brief DateTime and some operations supported on DateTime.
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
MODULE mtime_datetime
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_associated, c_f_pointer, c_int, c_int64_t, c_loc, c_null_char, c_ptr
  USE mtime_constants, ONLY: max_datetime_str_len
  USE mtime_c_bindings, ONLY: &
    & datetime, my_comparedatetime, my_constructandcopydatetime, &
    & my_datetimetoposixstring, my_datetimetostring, my_deallocatedatetime, &
    & my_getdatetimefromjulianday, my_getdayofyearfromdatetime, &
    & my_getjuliandayfromdatetime, my_getnoofdaysinmonthdatetime, &
    & my_getnoofdaysinyeardatetime, my_getnoofsecondselapsedindaydatetime, &
    & my_getnoofsecondselapsedinmonthdatetime, my_newdatetime, &
    & my_newrawdatetime, my_replacedatetime
  USE mtime_julianday
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: datetime
  PUBLIC :: newDatetime
  PUBLIC :: deallocateDatetime
  PUBLIC :: datetimeToString
  PUBLIC :: datetimeToPosixString
  PUBLIC :: getNoOfDaysInMonthDateTime
  PUBLIC :: getNoOfDaysInYearDateTime
  PUBLIC :: getDayOfYearFromDateTime
  PUBLIC :: getNoOfSecondsElapsedInMonthDateTime
  PUBLIC :: getNoOfSecondsElapsedInDayDateTime
  PUBLIC :: getJulianDayFromDatetime
  PUBLIC :: getDatetimeFromJulianDay
  !
  PUBLIC :: min, max
  PUBLIC :: ASSIGNMENT(=)
  PUBLIC :: OPERATOR(>)
  PUBLIC :: OPERATOR(<)
  PUBLIC :: OPERATOR(<=)
  PUBLIC :: OPERATOR(>=)
  PUBLIC :: OPERATOR(==)
  PUBLIC :: OPERATOR(/=)
  !
  INTERFACE newDatetime
    MODULE PROCEDURE newdatetimefromstring
    MODULE PROCEDURE newdatetimefromraw
    MODULE PROCEDURE newdatetimefromraw_yi8
    MODULE PROCEDURE newdatetimefromconstructandcopy
  END INTERFACE newDatetime
  !
  INTERFACE min
    MODULE PROCEDURE datetime_min
  END INTERFACE min
  !
  INTERFACE max
    MODULE PROCEDURE datetime_max
  END INTERFACE max
  !
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE replacedatetime
  END INTERFACE ASSIGNMENT(=)
  !
  INTERFACE OPERATOR(>)
    MODULE PROCEDURE datetime_gt
  END INTERFACE OPERATOR(>)
  !
  INTERFACE OPERATOR(<)
    MODULE PROCEDURE datetime_lt
  END INTERFACE OPERATOR(<)
  !
  INTERFACE OPERATOR(<=)
    MODULE PROCEDURE datetime_lt_or_eq
  END INTERFACE OPERATOR(<=)
  !
  INTERFACE OPERATOR(>=)
    MODULE PROCEDURE datetime_gt_or_eq
  END INTERFACE OPERATOR(>=)
  !
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE datetime_eq
  END INTERFACE OPERATOR(==)
  !
  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE datetime_ne
  END INTERFACE OPERATOR(/=)
  !
CONTAINS
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newdatetimefromstring(string, errno) RESULT(ret_datetime) !OK-TESTED
    TYPE(datetime), POINTER :: ret_datetime
    CHARACTER(len=*), INTENT(in) :: string
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_newdatetime(TRIM(ADJUSTL(string))//c_null_char)
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 1, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_datetime)
  END FUNCTION newdatetimefromstring
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newdatetimefromraw_yi8(year, month, day, hour, minute, second, ms, errno) RESULT(ret_datetime) !OK-TESTED
    TYPE(datetime), POINTER :: ret_datetime
    INTEGER(c_int64_t), INTENT(in) :: year
    INTEGER(c_int), INTENT(in) :: month, day, hour, minute, second, ms
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_newrawdatetime(year, month, day, hour, minute, second, ms)
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 2, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_datetime)
  END FUNCTION newdatetimefromraw_yi8
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newdatetimefromraw(year, month, day, hour, minute, second, ms, errno) RESULT(ret_datetime) !OK-TESTED
    TYPE(datetime), POINTER :: ret_datetime
    INTEGER(c_int), INTENT(in) :: year, month, day, hour, minute, second, ms
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_newrawdatetime(INT(year, c_int64_t), month, day, hour, minute, second, ms)
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 3, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_datetime)
  END FUNCTION newdatetimefromraw
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newdatetimefromconstructandcopy(src, errno) RESULT(dest) !OK-TESTED.
    TYPE(datetime), POINTER :: dest
    TYPE(datetime), TARGET :: src
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_constructandcopydatetime(C_LOC(src))
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 4, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, dest)
  END FUNCTION newdatetimefromconstructandcopy

  ! @return Minimum of two date time. In case of equality we return @p a.
  !
  ! @note This function does not return a copy of one of the arguments
  ! but returns only the corresponding result point to avoid
  ! unnecessary deallocate calls.
  RECURSIVE FUNCTION datetime_min(a, b) RESULT(res)
    TYPE(datetime), POINTER :: res
    TYPE(datetime), POINTER :: a, b

    IF (a > b) THEN
      res => b
    ELSE
      res => a
    END IF
  END FUNCTION datetime_min

  ! @return Maximum of two date time. In case of equality we return @p a.
  !
  ! @note This function does not return a copy of one of the arguments
  ! but returns only the corresponding result point to avoid
  ! unnecessary deallocate calls.
  RECURSIVE FUNCTION datetime_max(a, b) RESULT(res)
    TYPE(datetime), POINTER :: res
    TYPE(datetime), POINTER :: a, b

    IF (a < b) THEN
      res => b
    ELSE
      res => a
    END IF
  END FUNCTION datetime_max

  !>
  !! @brief Destructor of DateTime.
  !!
  !! @param  my_datetime
  !!         A pointer to type datetime. my_datetime is deallocated.
  RECURSIVE SUBROUTINE deallocateDatetime(my_datetime) !OK-TESTED.
    TYPE(datetime), POINTER :: my_datetime
    CALL my_deallocatedatetime(C_LOC(my_datetime))
    NULLIFY (my_datetime)
  END SUBROUTINE deallocateDatetime
  !>
  !! @brief Get DateTime as a string.
  !!
  !! datetimeToString returns a string in IS08601 compliant (and extended) format.
  !!
  !! @param[in]  my_datetime
  !!         A pointer to struct _datetime. The datetime to be converted to string.
  !!
  !! @param[out]  string
  !!         String where datetime is to be written.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE datetimeToString(my_datetime, string, errno) !OK-TESTED
    TYPE(datetime), TARGET, INTENT(in) :: my_datetime
    CHARACTER(len=max_datetime_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_datetimetostring(C_LOC(my_datetime), string)
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 6, C_ASSOCIATED(dummy_ptr))
    IF (.NOT. C_ASSOCIATED(dummy_ptr)) THEN
      string = '<null>'
    ELSE
      char_loop: DO i = 1, LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    END IF
  END SUBROUTINE datetimeToString
  !>
  !! @brief Get DateTime in 'struct tm' format and return as a string.
  !!
  !! Only dates between and including 1582-10-15 TO 9999-12-31 supported.
  !!
  !! @param  my_datetime
  !!         An object of type datetime. The datetime to be converted to string.
  !!
  !! @param  string
  !!         String where datetime is to be written.
  !!
  !! @param  fmtstr
  !!         Desired Format string. CRITICAL: Inappropriate fmt string will cause dump.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE datetimeToPosixString(my_datetime, string, fmtstr, errno) !OK-TESTED.
    TYPE(datetime), TARGET, INTENT(in) :: my_datetime
    CHARACTER(len=max_datetime_str_len) :: string
    CHARACTER(len=*) :: fmtstr
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_datetimetoposixstring(C_LOC(my_datetime), string, fmtstr)
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 7, C_ASSOCIATED(dummy_ptr))
    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE datetimeToPosixString
  !
  ! NOTE: Do not call the function using replacedatetime(.,.) directly; Use overloaded '=' instead as in "dest = src"
  RECURSIVE SUBROUTINE replacedatetime(dest, src) !OK-TESTED.
    TYPE(datetime), TARGET, INTENT(inout) :: dest
    TYPE(datetime), TARGET, INTENT(in) :: src
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_replacedatetime(C_LOC(src), C_LOC(dest))
  END SUBROUTINE replacedatetime
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION datetime_gt(op1, op2) RESULT(gt) !OK-TESTED.
    LOGICAL :: gt
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(C_LOC(op1), C_LOC(op2))
    gt = (ret == 1_c_int)
  END FUNCTION datetime_gt
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION datetime_lt(op1, op2) RESULT(lt) !OK-TESTED.
    LOGICAL :: lt
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(C_LOC(op1), C_LOC(op2))
    lt = (ret == -1_c_int)
  END FUNCTION datetime_lt
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION datetime_lt_or_eq(op1, op2) RESULT(lt_or_eq) !OK-TESTED.
    LOGICAL :: lt_or_eq
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(C_LOC(op1), C_LOC(op2))
    lt_or_eq = (ret == 0_c_int) .OR. (ret == -1_c_int)
  END FUNCTION datetime_lt_or_eq
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION datetime_gt_or_eq(op1, op2) RESULT(gt_or_eq) !OK-TESTED
    LOGICAL :: gt_or_eq
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(C_LOC(op1), C_LOC(op2))
    gt_or_eq = (ret == 0_c_int) .OR. (ret == 1_c_int)
  END FUNCTION datetime_gt_or_eq
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION datetime_eq(op1, op2) RESULT(eq) !OK-TESTED.
    LOGICAL :: eq
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(C_LOC(op1), C_LOC(op2))
    eq = ret == 0_c_int
  END FUNCTION datetime_eq
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION datetime_ne(op1, op2) RESULT(ne) !OK-TESTED.
    LOGICAL :: ne
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(C_LOC(op1), C_LOC(op2))
    ne = ret /= 0_c_int
  END FUNCTION datetime_ne
  !>
  !! @brief Get nod (number of Days) in the month of DateTime.
  !!
  !! Routine getNoOfDaysInMonthDateTime returns number of days in the month of DateTime. This routine
  !! supports all calendar types.
  !!
  !! For eg. the number of days for 2001-10-15T00:00:00.000 will be 31 for Gregorian Calendar.
  !! Similarly, this value will be 30 for Calendar of type 360 day-Calendar.
  !!
  !! @param[in]  dt
  !!         A pointer to type datetime.
  !!
  !! @return nod
  !!         Integer value of nod. The value depends on the month and the calendar type. Zero indicates error.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION getNoOfDaysInMonthDateTime(dt, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    INTEGER(c_int) :: getNoOfDaysInMonthDateTime
    INTEGER, OPTIONAL:: errno
    getNoOfDaysInMonthDateTime = my_getnoofdaysinmonthdatetime(C_LOC(dt))
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 15, getNoOfDaysInMonthDateTime /= 0)
  END FUNCTION getNoOfDaysInMonthDateTime
  !>
  !! @brief Get number of days in the Year of DateTime.
  !!
  !! Routine getNoOfDaysInYearDateTime returns number of days in the Year of DateTime. This routine
  !! supports all calendar types.
  !!
  !! Number of days returned will depend on the calendar type and if applicable, leap v/s non leap year.
  !!
  !! @param[in]  dt
  !!         A pointer to type datetime.
  !!
  !! @return nod
  !!         Integer value of nod. The value depends on the year and the calendar type. Zero indicates error.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION getNoOfDaysInYearDateTime(dt, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    INTEGER(c_int) :: getNoOfDaysInYearDateTime
    INTEGER, OPTIONAL:: errno
    getNoOfDaysInYearDateTime = my_getnoofdaysinyeardatetime(C_LOC(dt))
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 16, getNoOfDaysInYearDateTime /= 0)
  END FUNCTION getNoOfDaysInYearDateTime
  !>
  !! @brief Get the 'day-of-year' value of a DateTime.
  !!
  !! Routine getDayOfYearFromDateTime returns Day of Year for the DateTime. This routine supports
  !! all Calendar types.
  !!
  !! For eg. the day of year value for 2001-10-15T00:00:00.000 will be 288 for Gregorian Calendar.
  !! Similarly, this value will be 285 for Calendar of type 360 day-Calendar.
  !!
  !! @param[in]  dt
  !!         A pointer to type datetime. Retrieve the 'day-of-year' from this DT object.
  !!
  !! @return doy
  !!         Integer value of doy. The value depends on the calendar type. Zero indicates error.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION getDayOfYearFromDateTime(dt, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    INTEGER(c_int) :: getDayOfYearFromDateTime
    INTEGER, OPTIONAL:: errno
    getDayOfYearFromDateTime = my_getdayofyearfromdatetime(C_LOC(dt))
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 17, getDayOfYearFromDateTime /= 0)
  END FUNCTION getDayOfYearFromDateTime
  !>
  !! @brief Get number of seconds elapsed in the month of DateTime.
  !!
  !! @param[in] dt
  !!         A pointer to type datetime.
  !!
  !! @return no_of_seconds
  !!         int(i8) value of no_of_seconds. -1 indicates error.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION getNoOfSecondsElapsedInMonthDateTime(dt, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    INTEGER(c_int64_t) :: getNoOfSecondsElapsedInMonthDateTime
    INTEGER, OPTIONAL:: errno
    getNoOfSecondsElapsedInMonthDateTime = my_getnoofsecondselapsedinmonthdatetime(C_LOC(dt))
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 18, getNoOfSecondsElapsedInMonthDateTime /= -1)
  END FUNCTION getNoOfSecondsElapsedInMonthDateTime
  !>
  !! @brief Get number of seconds elapsed in the day of DateTime.
  !!
  !! @param[in]  dt
  !!         A pointer to type datetime.
  !!
  !! @return no_of_seconds
  !!         int value of no_of_seconds. -1 indicates error.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION getNoOfSecondsElapsedInDayDateTime(dt, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    INTEGER(c_int) :: getNoOfSecondsElapsedInDayDateTime
    INTEGER, OPTIONAL:: errno
    getNoOfSecondsElapsedInDayDateTime = my_getnoofsecondselapsedindaydatetime(C_LOC(dt))
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 19, getNoOfSecondsElapsedInDayDateTime /= -1)
  END FUNCTION getNoOfSecondsElapsedInDayDateTime
  !>
  !! @brief Get the Julian Day from DateTime.
  !!
  !! The routine getJulianDayFromDateTime returns the equivalent Julian date to DateTime. Internally
  !! it calls translation routines based on Calendar type.
  !!
  !! @param[in]  dt
  !!         A pointer to type datetime. The DT's value is converted to julian day value.
  !!
  !! @param[out]  jd
  !!         A pointer to type julianday. JD where the converted value is stored.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE getJulianDayFromDatetime(dt, jd, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    TYPE(julianday), TARGET :: jd
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_getjuliandayfromdatetime(C_LOC(dt), C_LOC(jd))
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 20, C_ASSOCIATED(dummy_ptr))
  END SUBROUTINE getJulianDayFromDatetime
  !>
  !! @brief Get the DateTime from Julian Day.
  !!
  !! The routine getDateTimeFromJulianDay returns the equivalent DateTime to Julian date. Internally
  !! it calls translation routines based on Calendar type.
  !!
  !! @param[in]  jd
  !!         A pointer to type julianday. The JD's value is converted to julian day value.
  !!
  !! @param[out]  dt
  !!         A pointer to type datetime. The DT where the converted value is stored.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE getDatetimeFromJulianDay(jd, dt, errno)
    TYPE(julianday), TARGET, INTENT(in) :: jd
    TYPE(datetime), TARGET, INTENT(out) :: dt
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_getdatetimefromjulianday(C_LOC(jd), C_LOC(dt))
    IF (PRESENT(errno)) errno = MERGE(0, 4*100 + 21, C_ASSOCIATED(dummy_ptr))
  END SUBROUTINE getDatetimeFromJulianDay
  !
END MODULE mtime_datetime
!>
!! @brief TimeDelta and some operations supported on TimeDelta.
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
MODULE mtime_timedelta
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_associated, c_char, c_double, c_f_pointer, c_float, c_int, c_int32_t, &
    & c_int64_t, c_loc, c_null_char, c_ptr
  USE mtime_constants, ONLY: max_timedelta_str_len
  USE mtime_c_bindings, ONLY: &
    & divisionquotienttimespan, my_addtimedeltatodate, &
    & my_addtimedeltatodatetime, my_comparetimedelta, &
    & my_constructandcopytimedelta, my_deallocatetimedelta, &
    & my_dividedatetimedifferenceinseconds, my_dividetimedeltainseconds, &
    & my_dividetwodatetimediffsinseconds, &
    & my_elementwiseaddtimedeltatotimedelta, &
    & my_elementwisescalarmultiplytimedelta, &
    & my_elementwisescalarmultiplytimedeltadp, my_getptstringfromhours, &
    & my_getptstringfromminutes, my_getptstringfromms, &
    & my_getptstringfromsecondsdouble, my_getptstringfromsecondsfloat, &
    & my_getptstringfromsecondsint, my_gettimedeltafromdate, &
    & my_gettimedeltafromdatetime, my_gettimedeltafromdatetime, &
    & my_gettotalmillisecondstimedelta, my_gettotalsecondstimedelta, &
    & my_modulotimedelta, my_newrawtimedelta, my_newtimedeltafromstring, &
    & my_timedeltatojuliandelta, my_timedeltatostring, timedelta
  USE mtime_date
  USE mtime_datetime
  USE mtime_juliandelta
  USE mtime_time
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: timedelta
  PUBLIC :: newTimedelta
  PUBLIC :: deallocateTimedelta
  PUBLIC :: getTimeDeltaFromDate
  PUBLIC :: getTimeDeltaFromDateTime
  PUBLIC :: getTotalMilliSecondsTimeDelta
  PUBLIC :: getTotalSecondsTimeDelta
  PUBLIC :: timedeltaToString
  PUBLIC :: moduloTimedelta
  PUBLIC :: getPTStringFromMS
  PUBLIC :: getPTStringFromSeconds
  PUBLIC :: getPTStringFromMinutes
  PUBLIC :: getPTStringFromHours
  PUBLIC :: timeDeltaToJulianDelta
  PUBLIC :: divideTimeDeltaInSeconds
  PUBLIC :: divideTwoDatetimeDiffsInSeconds
  PUBLIC :: divideDatetimeDifferenceInSeconds
  PUBLIC :: OPERATOR(+)
  PUBLIC :: OPERATOR(-)
  PUBLIC :: OPERATOR(*)
  PUBLIC :: OPERATOR(>)
  PUBLIC :: OPERATOR(<)
  PUBLIC :: OPERATOR(<=)
  PUBLIC :: OPERATOR(>=)
  PUBLIC :: OPERATOR(==)
  PUBLIC :: OPERATOR(/=)
  !
  INTERFACE newTimedelta
    MODULE PROCEDURE newtimedeltafromstring
    MODULE PROCEDURE newtimedeltafromraw
    MODULE PROCEDURE newtimedeltafromraw_yi8
    MODULE PROCEDURE newtimedeltafromconstructandcopy
  END INTERFACE newTimedelta
  !
  INTERFACE OPERATOR(+)
    MODULE PROCEDURE addtimedeltatodatetime
    MODULE PROCEDURE adddatetimetotimedelta
    MODULE PROCEDURE addtimedeltatodate
    MODULE PROCEDURE adddatetotimedelta
    MODULE PROCEDURE elementwiseAddTimeDeltatoTimeDelta
  END INTERFACE OPERATOR(+)
  !
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE getTimeDeltaFromDate
    MODULE PROCEDURE getTimeDeltaFromDateTime
  END INTERFACE OPERATOR(-)
  !
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE elementwiseScalarMultiplyTimeDelta
    MODULE PROCEDURE elementwiseScalarMultiplyTimeDeltaInv
    MODULE PROCEDURE elementwiseScalarMultiplyTimeDelta_long
    MODULE PROCEDURE elementwiseScalarMultiplyTimeDeltaInv_long
    MODULE PROCEDURE elementwiseScalarMultiplyTimeDelta_real
  END INTERFACE OPERATOR(*)
  !
  INTERFACE OPERATOR(>)
    MODULE PROCEDURE timedelta_gt
  END INTERFACE OPERATOR(>)
  !
  INTERFACE OPERATOR(<)
    MODULE PROCEDURE timedelta_lt
  END INTERFACE OPERATOR(<)
  !
  INTERFACE OPERATOR(<=)
    MODULE PROCEDURE timedelta_lt_or_eq
  END INTERFACE OPERATOR(<=)
  !
  INTERFACE OPERATOR(>=)
    MODULE PROCEDURE timedelta_gt_or_eq
  END INTERFACE OPERATOR(>=)
  !
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE timedelta_eq
  END INTERFACE OPERATOR(==)
  !
  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE timedelta_ne
  END INTERFACE OPERATOR(/=)
  !
  INTERFACE getPTStringFromSeconds
    MODULE PROCEDURE getPTStringFromSecondsInt
    MODULE PROCEDURE getPTStringFromSecondsFloat
    MODULE PROCEDURE getPTStringFromSecondsDouble
  END INTERFACE getPTStringFromSeconds
  !
CONTAINS
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newtimedeltafromstring(string, errno) RESULT(ret_timedelta) !OK-TESTED.
    TYPE(timedelta), POINTER :: ret_timedelta
    CHARACTER(len=*), INTENT(in) :: string
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_newtimedeltafromstring(TRIM(ADJUSTL(string))//c_null_char)
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 1, C_ASSOCIATED(c_pointer))
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, ret_timedelta)
    ELSE
      NULLIFY (ret_timedelta)
    END IF
  END FUNCTION newtimedeltafromstring
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newtimedeltafromraw(sign, year, month, day, hour, minute, second, ms, errno) RESULT(ret_timedelta)
    TYPE(timedelta), POINTER :: ret_timedelta
    CHARACTER(len=*), INTENT(in) :: sign
    INTEGER(c_int), INTENT(in) :: year, month, day, hour, minute, second, ms
    TYPE(c_ptr) :: c_pointer
    CHARACTER(c_char) ::c_sign
    INTEGER, OPTIONAL:: errno
    c_sign = SIGN(1:1)
    c_pointer = my_newrawtimedelta(c_sign, INT(year, c_int64_t), month, day, hour, minute, second, ms)
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 2, C_ASSOCIATED(c_pointer))
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, ret_timedelta)
    ELSE
      NULLIFY (ret_timedelta)
    END IF
  END FUNCTION newtimedeltafromraw
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newtimedeltafromraw_yi8(sign, year, month, day, hour, minute, second, ms, errno) RESULT(ret_timedelta)
    TYPE(timedelta), POINTER :: ret_timedelta
    CHARACTER(len=*), INTENT(in) :: sign
    INTEGER(c_int64_t), INTENT(in) :: year
    INTEGER(c_int), INTENT(in) :: month, day, hour, minute, second, ms
    TYPE(c_ptr) :: c_pointer
    CHARACTER(c_char) ::c_sign
    INTEGER, OPTIONAL:: errno
    c_sign = SIGN(1:1)
    c_pointer = my_newrawtimedelta(c_sign, year, month, day, hour, minute, second, ms)
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 3, C_ASSOCIATED(c_pointer))
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, ret_timedelta)
    ELSE
      NULLIFY (ret_timedelta)
    END IF
  END FUNCTION newtimedeltafromraw_yi8
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newtimedeltafromconstructandcopy(src, errno) RESULT(dest) !OK-TESTED.
    TYPE(timedelta), POINTER :: dest
    TYPE(timedelta), TARGET :: src
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    c_pointer = my_constructandcopytimedelta(C_LOC(src))
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 4, C_ASSOCIATED(c_pointer))
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, dest)
    ELSE
      NULLIFY (dest)
    END IF
  END FUNCTION newtimedeltafromconstructandcopy
  !>
  !! @brief Destructor of TimeDelta.
  !!
  !! @param  my_timedelta
  !!         A pointer to typetimedelta. my_timedelta is deallocated.
  RECURSIVE SUBROUTINE deallocateTimedelta(my_timedelta) !OK-TESTED.
    TYPE(timedelta), POINTER :: my_timedelta
    CALL my_deallocatetimedelta(C_LOC(my_timedelta))
    NULLIFY (my_timedelta)
  END SUBROUTINE deallocateTimedelta
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION timedelta_gt(op1, op2) RESULT(gt)
    LOGICAL :: gt
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(C_LOC(op1), C_LOC(op2))
    gt = ret == 1_c_int
  END FUNCTION timedelta_gt
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION timedelta_lt(op1, op2) RESULT(lt)
    LOGICAL :: lt
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(C_LOC(op1), C_LOC(op2))
    lt = ret == -1_c_int
  END FUNCTION timedelta_lt
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION timedelta_lt_or_eq(op1, op2) RESULT(lt_or_eq)
    LOGICAL :: lt_or_eq
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(C_LOC(op1), C_LOC(op2))
    lt_or_eq = ret == 0_c_int .OR. ret == -1_c_int
  END FUNCTION timedelta_lt_or_eq
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION timedelta_gt_or_eq(op1, op2) RESULT(gt_or_eq)
    LOGICAL :: gt_or_eq
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(C_LOC(op1), C_LOC(op2))
    gt_or_eq = ret == 0_c_int .OR. ret == 1_c_int
  END FUNCTION timedelta_gt_or_eq
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION timedelta_eq(op1, op2) RESULT(eq)
    LOGICAL :: eq
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(C_LOC(op1), C_LOC(op2))
    eq = ret == 0_c_int
  END FUNCTION timedelta_eq
  !
  MTIME_PURE_IF_C_LOC_IS_PURE RECURSIVE FUNCTION timedelta_ne(op1, op2) RESULT(ne)
    LOGICAL :: ne
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(C_LOC(op1), C_LOC(op2))
    ne = ret /= 0_c_int
  END FUNCTION timedelta_ne
  !>
  !! @brief Get the TimeDelta between two Dates op1 and op2 as (op1-op22).
  !!
  !! Routine getTimeDeltaFromDate 'subtracts' two Dates and returns the TimeDelta between
  !! them. Internally, Dates are converted to DateTimes and then delta is calculated using
  !! getTimeDeltaFromDateTime().
  !!
  !! This routine  handles all supported Calendar types; i.e. the translation from Calendar date
  !! to Julian date and conversion from Julian Delta to normal TimeDetla is Calendar-type dependent.
  !! For eg. for Calendar type Gregorian, the TimeDelta between 2001-02-01 and 2001-01-01 will be 1 month.
  !! Similarly, for Calendar of type 360-Day-Calendar, the TimeDelta will be 1 month. It must be noted
  !! however, that the two dates differ by 31 and 30 days respectively.
  !!
  !! @param  op1
  !!         A pointer to type date.
  !!
  !! @param  op2
  !!         A pointer to type date.
  !!
  !! @return ret
  !!         A pointer to TimeDelta containing the result of subtraction.
  RECURSIVE FUNCTION getTimeDeltaFromDate(op1, op2) RESULT(ret) !OK-TESTED.
    TYPE(timedelta), TARGET :: ret
    TYPE(date), TARGET, INTENT(in)  :: op1, op2
    TYPE(c_ptr) :: dummy_ptr
    ret = timedelta(0, "+", 0, 0, 0, 0, 0, 0, 0)
    dummy_ptr = my_gettimedeltafromdate(C_LOC(op1), C_LOC(op2), C_LOC(ret))
  END FUNCTION getTimeDeltaFromDate
  !>
  !! @brief Get the TimeDelta between two DateTimes op1 and op2 as (op1-op2).
  !!
  !! Routine getTimeDeltaFromDateTime 'subtracts' two DateTime's and returns the TimeDelta between
  !! them. Each datetime is converted to an equivalent Julian Date. Subtraction is then performed
  !! on Julian axis. The "Julian delta" is finally converted back to normal calendar delta.
  !!
  !! This routine handles all supported Calendar types; i.e. the translation from Calendar date
  !! to Julian date and conversion from Julian Delta to normal TimeDetla is Calendar-type dependent.
  !! For eg. for Calendar type Gregorian, the TimeDelta between 2001-02-01T00:00:00.000 and
  !! 2001-01-01T00:00:00.000 will be 1 month. Similarly, for Calendar of type 360-Day-Calendar,
  !! the TimeDelta will be 1 month. It must be noted however, that the two dates differ by 31 and
  !! 30 days respectively.
  !!
  !! @param  op1
  !!         A pointer to type datetime.
  !!
  !! @param  op2
  !!         A pointer to type datetime.
  !!
  !! @return ret
  !!        A pointer to TimeDelta containing the result of subtraction.
  RECURSIVE FUNCTION getTimeDeltaFromDateTime(op1, op2) RESULT(ret) !OK-TESTED.
    TYPE(timedelta), TARGET :: ret
    TYPE(datetime), TARGET, INTENT(in) :: op1, op2
    TYPE(c_ptr) :: dummy_ptr
    ret = timedelta(0, "+", 0, 0, 0, 0, 0, 0, 0)
    dummy_ptr = my_gettimedeltafromdatetime(C_LOC(op1), C_LOC(op2), C_LOC(ret))
  END FUNCTION getTimeDeltaFromDateTime
  !>
  !! @brief Get total number of milliseconds in timedelta.
  !!
  !! Routine getTotalMilliSecondsTimeDelta returns the total number of milliseconds in TimeDelta.
  !! Notice that TimeDelta is not uniquely defined but depends on the definition of corresponding
  !! DateTime. TimeDelta is first converted to corresponding delta on the Julian axis. Julian delta
  !! is finally converted to the correct millisecond value.
  !!
  !! @param[in]  td
  !!         A pointer to type timedelta. Retrieve the number of milliseconds in this TD object.
  !!
  !! @param[in]  dt
  !!         A pointer to type datetime. Reference Datetime for the TD.
  !!
  !! @param[out] errno
  !!         Optional, error message
  !!
  !! @return totalmilliSeconds
  !!         Integer value of totalmilliSeconds. 0 indicates error.
  !!
  !! WARNING: TD 0 is error. If you know your TD is 0, ignore the error flag.
  RECURSIVE FUNCTION getTotalMilliSecondsTimeDelta(td, dt, errno)  !OK-TESTED.
    INTEGER(c_int64_t) :: getTotalMilliSecondsTimeDelta
    TYPE(timedelta), TARGET, INTENT(in):: td
    TYPE(datetime), TARGET, INTENT(in) :: dt
    INTEGER, OPTIONAL:: errno
    getTotalMilliSecondsTimeDelta = my_gettotalmillisecondstimedelta(C_LOC(td), C_LOC(dt))
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 8, getTotalMilliSecondsTimeDelta /= 0)
  END FUNCTION getTotalMilliSecondsTimeDelta
  !>
  !! @brief Get total number of seconds in timedelta.
  !!
  !! Routine getTotalSecondsTimeDelta returns the total number of seconds in TimeDelta. Notice that TimeDelta
  !! is not uniquely defined but depends on the definition of corresponding DateTime. Internally, number of seconds
  !! is calculated by calling the routine getTotalMilliSecondsTimeDelta() and then converting the millisecond value
  !! to seconds by dividing it by 1000.
  !!
  !! @param[in]  td
  !!         A pointer to struct _timedelta. Retrieve the number of seconds in this TD object.
  !!
  !! @param[in]  dt
  !!         A pointer to struct _datetime. Reference Datetime for the TD.
  !!
  !! @param[out] errno
  !!         Optional, error message
  !!
  !! @return totalSeconds
  !!         Integer value of totalSeconds. 0 indicates error.
  !!
  !! WARNING: TD 0 is error. If you know your TD is 0, ignore the error flag.
  RECURSIVE FUNCTION getTotalSecondsTimeDelta(td, dt, errno) !OK-TESTED.
    INTEGER(c_int64_t) :: getTotalSecondsTimeDelta
    TYPE(timedelta), TARGET, INTENT(in) :: td
    TYPE(datetime), TARGET, INTENT(in) :: dt
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    getTotalSecondsTimeDelta = my_gettotalsecondstimedelta(C_LOC(td), C_LOC(dt))
    ! error handling with "errno" not yet implemented.
  END FUNCTION getTotalSecondsTimeDelta
  !>
  !! @brief Get TimeDelta as an extended string.
  !!
  !! timedeltaToString returns a string in IS08601 compliant (and extended) format.
  !!
  !! @param[in] my_timedelta
  !!         A pointer to type timedelta. The timedelta to be converted to string.
  !!
  !! @param[out] string
  !!         String where timedelta is to be written.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE timedeltaToString(my_timedelta, string, errno) !OK-TESTED.
    TYPE(timedelta), TARGET :: my_timedelta
    CHARACTER(len=max_timedelta_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string = ''
    dummy_ptr = my_timedeltatostring(C_LOC(my_timedelta), string)
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 10, C_ASSOCIATED(dummy_ptr))
    IF (.NOT. C_ASSOCIATED(dummy_ptr)) THEN
      string = '<null>'
    ELSE
      char_loop: DO i = 1, LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    END IF
  END SUBROUTINE timedeltaToString
  !
  RECURSIVE FUNCTION addTimedeltaToDatetime(op1, op2) RESULT(ret) !OK-TESTED.
    TYPE(datetime), TARGET :: ret
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    ret = datetime(date(0, 0, 0), time(0, 0, 0, 0))
    dummy_ptr = my_addtimedeltatodatetime(C_LOC(op1), C_LOC(op2), C_LOC(ret))
  END FUNCTION addTimedeltaToDatetime
  !
  RECURSIVE FUNCTION addDatetimeToTimedelta(op2, op1) RESULT(ret) !OK-TESTED.
    TYPE(datetime), TARGET :: ret
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    ret = datetime(date(0, 0, 0), time(0, 0, 0, 0))
    dummy_ptr = my_addtimedeltatodatetime(C_LOC(op1), C_LOC(op2), C_LOC(ret))
  END FUNCTION addDatetimeToTimedelta
  !
  RECURSIVE FUNCTION addTimedeltaToDate(op1, op2) RESULT(ret) !OK-TESTED.
    TYPE(date), TARGET :: ret
    TYPE(date), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    ret = date(0, 0, 0)
    dummy_ptr = my_addtimedeltatodate(C_LOC(op1), C_LOC(op2), C_LOC(ret))
  END FUNCTION addTimedeltaToDate
  !
  RECURSIVE FUNCTION addDateToTimedelta(op2, op1) RESULT(ret) !OK-TESTED.
    TYPE(date), TARGET :: ret
    TYPE(date), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    ret = date(0, 0, 0)
    dummy_ptr = my_addtimedeltatodate(C_LOC(op1), C_LOC(op2), C_LOC(ret))
  END FUNCTION addDateToTimedelta
  !
  RECURSIVE FUNCTION elementwiseScalarMultiplyTimeDelta(base_td, ilambda) RESULT(scaled_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: scaled_td
    INTEGER(c_int32_t), INTENT(in) :: ilambda
    TYPE(timedelta), TARGET, INTENT(in) :: base_td
    INTEGER(c_int64_t) :: lambda
    TYPE(c_ptr) :: dummy_ptr
    lambda = INT(ilambda, c_int64_t)
    scaled_td = timedelta(0, "+", 0, 0, 0, 0, 0, 0, 0)
    dummy_ptr = my_elementwisescalarmultiplytimedelta(C_LOC(base_td), lambda, C_LOC(scaled_td))
  END FUNCTION elementwiseScalarMultiplyTimeDelta
  !
  RECURSIVE FUNCTION elementwiseScalarMultiplyTimeDeltaInv(ilambda, base_td) RESULT(scaled_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: scaled_td
    INTEGER(c_int32_t), INTENT(in) :: ilambda
    TYPE(timedelta), TARGET, INTENT(in) :: base_td
    INTEGER(c_int64_t) :: lambda
    TYPE(c_ptr) :: dummy_ptr
    lambda = INT(ilambda, c_int64_t)
    scaled_td = timedelta(0, "+", 0, 0, 0, 0, 0, 0, 0)
    dummy_ptr = my_elementwisescalarmultiplytimedelta(C_LOC(base_td), lambda, C_LOC(scaled_td))
  END FUNCTION elementwiseScalarMultiplyTimeDeltaInv
  !
  RECURSIVE FUNCTION elementwiseScalarMultiplyTimeDelta_long(base_td, lambda) RESULT(scaled_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: scaled_td
    INTEGER(c_int64_t), INTENT(in) :: lambda
    TYPE(timedelta), TARGET, INTENT(in) :: base_td
    TYPE(c_ptr) :: dummy_ptr
    scaled_td = timedelta(0, "+", 0, 0, 0, 0, 0, 0, 0)
    dummy_ptr = my_elementwisescalarmultiplytimedelta(C_LOC(base_td), lambda, C_LOC(scaled_td))
  END FUNCTION elementwiseScalarMultiplyTimeDelta_long
  !
  RECURSIVE FUNCTION elementwiseScalarMultiplyTimeDeltaInv_long(lambda, base_td) RESULT(scaled_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: scaled_td
    INTEGER(c_int64_t), INTENT(in) :: lambda
    TYPE(timedelta), TARGET, INTENT(in) :: base_td
    TYPE(c_ptr) :: dummy_ptr
    scaled_td = timedelta(0, "+", 0, 0, 0, 0, 0, 0, 0)
    dummy_ptr = my_elementwisescalarmultiplytimedelta(C_LOC(base_td), lambda, C_LOC(scaled_td))
  END FUNCTION elementwiseScalarMultiplyTimeDeltaInv_long
  !
  RECURSIVE FUNCTION elementwisescalarmultiplytimedelta_real(base_td, lambda) RESULT(scaled_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: scaled_td
    REAL(c_double), INTENT(in) :: lambda
    TYPE(timedelta), TARGET, INTENT(in) :: base_td
    TYPE(c_ptr) :: dummy_ptr
    scaled_td = timedelta(0, "+", 0, 0, 0, 0, 0, 0, 0)
    dummy_ptr = my_elementwisescalarmultiplytimedeltadp(C_LOC(base_td), lambda, C_LOC(scaled_td))
  END FUNCTION elementwisescalarmultiplytimedelta_real
  !
  RECURSIVE FUNCTION elementwiseAddTimeDeltatoTimeDelta(td1, td2) RESULT(added_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: added_td
    TYPE(timedelta), TARGET, INTENT(in) :: td1, td2
    TYPE(c_ptr) :: dummy_ptr
    added_td = timedelta(0, "+", 0, 0, 0, 0, 0, 0, 0)
    dummy_ptr = my_elementwiseaddtimedeltatotimedelta(C_LOC(td1), C_LOC(td2), C_LOC(added_td))
  END FUNCTION elementwiseAddTimeDeltatoTimeDelta
  !>
  !! @brief Returns modulo(a,p) and the quotient.
  !!
  !! @param[in]  a
  !!         A pointer to type timedelta.
  !!
  !! @param[in]  p
  !!         A pointer to type timedelta.
  !!
  !! @param[out]  quot
  !!         The quotient of a divided by p.
  !!
  !! @return rem
  !!       modulo(a, p)
  RECURSIVE FUNCTION moduloTimedelta(a, p, quot) RESULT(rem)
    TYPE(timedelta), TARGET, INTENT(in) :: a
    TYPE(timedelta), TARGET, INTENT(in) :: p
    INTEGER(c_int64_t), TARGET, INTENT(out) :: quot
    INTEGER(c_int64_t) :: rem
    rem = my_modulotimedelta(C_LOC(a), C_LOC(p), C_LOC(quot))
  END FUNCTION moduloTimedelta
  !>
  !! @brief Return a PT String corresponding to arbitrary number of milliseconds.
  !!
  !! getPTStringFromMS() translates ms values to ISO 8601 compliant timedelta string.
  !! Conversion of ms >= 86400000 and  ms <= -86400000 not supported.
  !!
  !! @param[in]  ms
  !!         An int64_t value to be translated.
  !!
  !! @param[out]  string
  !!         Translated string is written here.
  !!
  RECURSIVE SUBROUTINE getPTStringFromMS(ms, string, errno) !OK-TESTED.
    INTEGER(c_int64_t), INTENT(in) :: ms
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    dummy_ptr = my_getptstringfromms(ms, string)
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 11, C_ASSOCIATED(dummy_ptr))
    IF (.NOT. C_ASSOCIATED(dummy_ptr)) THEN
      string = '<null>'
    ELSE
      char_loop: DO i = 1, LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    END IF
  END SUBROUTINE getPTStringFromMS
  !
  RECURSIVE SUBROUTINE getPTStringFromSecondsInt(s, string, errno) !OK-TESTED.
    INTEGER(c_int64_t) :: s
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    dummy_ptr = my_getptstringfromsecondsint(s, string)
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 11, C_ASSOCIATED(dummy_ptr))
    IF (.NOT. C_ASSOCIATED(dummy_ptr)) THEN
      string = '<null>'
    ELSE
      char_loop: DO i = 1, LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    END IF
  END SUBROUTINE getPTStringFromSecondsInt
  !
  RECURSIVE SUBROUTINE getPTStringFromSecondsFloat(s, string, errno) !OK-TESTED.
    REAL(c_float) :: s
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    dummy_ptr = my_getptstringfromsecondsfloat(s, string)
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 11, C_ASSOCIATED(dummy_ptr))
    IF (.NOT. C_ASSOCIATED(dummy_ptr)) THEN
      string = '<null>'
    ELSE
      char_loop: DO i = 1, LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    END IF
  END SUBROUTINE getPTStringFromSecondsFloat
  !
  RECURSIVE SUBROUTINE getPTStringFromSecondsDouble(s, string, errno) !OK-TESTED.
    REAL(c_double) :: s
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    dummy_ptr = my_getptstringfromsecondsdouble(s, string)
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 11, C_ASSOCIATED(dummy_ptr))
    IF (.NOT. C_ASSOCIATED(dummy_ptr)) THEN
      string = '<null>'
    ELSE
      char_loop: DO i = 1, LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    END IF
  END SUBROUTINE getPTStringFromSecondsDouble
  !>
  !! @brief Return a PT String corresponding to arbitrary number of minutes.
  !!
  !! getPTStringFromMinutes() translates minutes values to ISO 8601 compliant timedelta string.
  !! Conversion of m >= 1440 and  m <= -1440 not supported.
  !!
  !! @param[in]  m
  !!         An int64_t value to be translated.
  !!
  !! @param[out] string
  !!         Translated string is written here.
  !!
  RECURSIVE SUBROUTINE getPTStringFromMinutes(m, string, errno) !OK-TESTED.
    INTEGER(c_int64_t) :: m
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    dummy_ptr = my_getptstringfromminutes(m, string)
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 11, C_ASSOCIATED(dummy_ptr))
    IF (.NOT. C_ASSOCIATED(dummy_ptr)) THEN
      string = '<null>'
    ELSE
      char_loop: DO i = 1, LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    END IF
  END SUBROUTINE getPTStringFromMinutes
  !>
  !! @brief Return a PT String corresponding to arbitrary number of Hours.
  !!
  !! getPTStringFromHours() translates hour values to ISO 8601 compliant timedelta string.
  !! Conversion of h >= 24 and  ms <= -24 not supported.
  !!
  !! @param[in]  h
  !!         An int64_t value to be translated.
  !!
  !! @param[out]  string
  !!         Translated string is written here.
  !!
  RECURSIVE SUBROUTINE getPTStringFromHours(h, string, errno) !OK-TESTED.
    INTEGER(c_int64_t)                   :: h
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    dummy_ptr = my_getptstringfromhours(h, string)
    IF (PRESENT(errno)) errno = MERGE(0, 5*100 + 11, C_ASSOCIATED(dummy_ptr))
    IF (.NOT. C_ASSOCIATED(dummy_ptr)) THEN
      string = '<null>'
    ELSE
      char_loop: DO i = 1, LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    END IF
  END SUBROUTINE getPTStringFromHours

  !>
  !! @brief Convert time delta to "Julian calendar delta".
  !!
  RECURSIVE SUBROUTINE timeDeltaToJulianDelta(td, dt, jd)
    TYPE(c_ptr) :: dummy_ptr
    TYPE(timedelta), TARGET, INTENT(in) :: td
    TYPE(datetime), TARGET, INTENT(in) :: dt
    TYPE(juliandelta), TARGET, INTENT(out) :: jd
    dummy_ptr = my_timedeltatojuliandelta(C_LOC(td), C_LOC(dt), C_LOC(jd))
  END SUBROUTINE timeDeltaToJulianDelta

  !>
  !! @brief division by seconds.
  !!
  !! @param[in]  dividend
  !!         A pointer to type timedelta
  !!
  !! @param[in]  divisor
  !!         A pointer to type timedelta
  !!
  !! @param[out]  quotient
  !!         A pointer to type divisionquotienttimespan
  !!
  RECURSIVE SUBROUTINE divideTimeDeltaInSeconds(dividend, divisor, quotient, errna)!OK-UNTESTED.
    TYPE(timedelta), TARGET, INTENT(in) :: dividend
    TYPE(timedelta), TARGET, INTENT(in) :: divisor
    TYPE(divisionquotienttimespan), TARGET, INTENT(out) :: quotient
    INTEGER, INTENT(out), OPTIONAL :: errna
    TYPE(c_ptr) :: dummy_ptr
    IF (PRESENT(errna)) errna = 0 ! FIXME: no_error
    dummy_ptr = my_dividetimedeltainseconds(C_LOC(dividend), C_LOC(divisor), C_LOC(quotient))
    IF (PRESENT(errna) .AND. .NOT. C_ASSOCIATED(dummy_ptr)) THEN
      errna = errna + 2  ! increment error number by 2, see below for an explanation.
    END IF
  END SUBROUTINE divideTimeDeltaInSeconds
  !>
  !! @brief division of two differences in datetimes.
  !!
  !! @param  dt1_dividend, dt2_dividend, dt1_divisor, dt2_divisor
  !!         Reference date (a pointer to struct _datetime).
  !!
  !! @param  intvlsec
  !!         Interval given in seconds.
  !!
  !! @return result of division. NULL indicates error.
  RECURSIVE SUBROUTINE divideTwoDatetimeDiffsInSeconds(dt1_dividend, dt2_dividend, &
      &                                      dt1_divisor, dt2_divisor,  &
      &                                      denominator, quotient)
    TYPE(datetime), TARGET, INTENT(in) :: dt1_dividend
    TYPE(datetime), TARGET, INTENT(in) :: dt2_dividend
    TYPE(datetime), TARGET, INTENT(in) :: dt1_divisor
    TYPE(datetime), TARGET, INTENT(in) :: dt2_divisor
    INTEGER(c_int64_t), TARGET, INTENT(out) :: denominator
    TYPE(divisionquotienttimespan), TARGET, INTENT(out) :: quotient
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_dividetwodatetimediffsinseconds(C_LOC(dt1_dividend), C_LOC(dt2_dividend),  &
        &                                                C_LOC(dt1_divisor), C_LOC(dt2_divisor),  &
        &                                                C_LOC(denominator), C_LOC(quotient))
  END SUBROUTINE divideTwoDatetimeDiffsInSeconds

  !>
  !! @brief division of an datetime interval by seconds.
  !!
  !! the datetime interval is calculated by dt1-dt2.
  !!
  !! @param[in]  dt1
  !!         A pointer to type datetime
  !!
  !! @param[in]  dt2
  !!         A pointer to type datetime
  !!
  !! @param[in]  divisor
  !!         A pointer to type timedelta
  !!
  !! @param[out]  quotient
  !!         A pointer to type divisionquotienttimespan
  !!
  RECURSIVE SUBROUTINE divideDatetimeDifferenceInSeconds(dt1, dt2, divisor, quotient, errna)
    TYPE(datetime), TARGET, INTENT(in) :: dt1
    TYPE(datetime), TARGET, INTENT(in) :: dt2
    TYPE(timedelta), TARGET, INTENT(in) :: divisor
    TYPE(divisionquotienttimespan), TARGET, INTENT(out) :: quotient
    INTEGER, INTENT(out), OPTIONAL :: errna
    TYPE(c_ptr) :: dummy_ptr
    IF (PRESENT(errna)) errna = 0 ! FIXME: no_error
    dummy_ptr = my_dividedatetimedifferenceinseconds(C_LOC(dt1), C_LOC(dt2), C_LOC(divisor), C_LOC(quotient))
    IF (PRESENT(errna) .AND. .NOT. C_ASSOCIATED(dummy_ptr)) THEN
      errna = errna + 2  ! increment error number by 2, see below for an explanation.
    END IF
  END SUBROUTINE divideDatetimeDifferenceInSeconds
  !
END MODULE mtime_timedelta
!>
!! @brief Definition of the basic event type and its methods.
!!
!! @details
!!
!___________________________________________________________________________________________________________
MODULE mtime_events
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_associated, c_bool, c_f_pointer, c_int64_t, c_loc, c_null_char, &
    & c_null_ptr, c_ptr
  USE mtime_constants, ONLY: max_eventname_str_len
  USE mtime_c_bindings, ONLY: &
    & event, my_constructandcopyevent, my_deallocateevent, my_eventtostring, &
    & my_geteventisfirstinday, my_geteventisfirstinmonth, &
    & my_geteventisfirstinyear, my_geteventislastinday, &
    & my_geteventislastinmonth, my_geteventislastinyear, my_geteventname, &
    & my_getnexteventisfirst, my_gettriggeredpreviouseventatdatetime, &
    & my_gettriggernexteventatdatetime, my_iscurrenteventactive, &
    & my_iseventnextinnextday, my_iseventnextinnextmonth, &
    & my_iseventnextinnextyear, my_newevent, my_neweventwithdatatypes
  USE mtime_datetime
  USE mtime_timedelta
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: event
  PUBLIC :: newEvent
  PUBLIC :: deallocateEvent
  PUBLIC :: eventToString
  PUBLIC :: isCurrentEventActive
  PUBLIC :: isEventNextInNextDay
  PUBLIC :: isEventNextInNextMonth
  PUBLIC :: isEventNextInNextYear
  PUBLIC :: getTriggerNextEventAtDateTime
  PUBLIC :: getTriggeredPreviousEventAtDateTime
  PUBLIC :: getEventId
  PUBLIC :: getEventName
  PUBLIC :: getEventReferenceDateTime
  PUBLIC :: getEventFirstDateTime
  PUBLIC :: getEventLastDateTime
  PUBLIC :: getEventInterval
  PUBLIC :: getNextEventIsFirst
  PUBLIC :: getEventisFirstInDay
  PUBLIC :: getEventisFirstInMonth
  PUBLIC :: getEventisFirstInYear
  PUBLIC :: getEventisLastInDay
  PUBLIC :: getEventisLastInMonth
  PUBLIC :: getEventisLastInYear
  !
  INTERFACE newEvent
    MODULE PROCEDURE newEventWithString
    MODULE PROCEDURE newEventWithDataTypes
    MODULE PROCEDURE constructAndCopyEvent
  END INTERFACE newEvent
  !
CONTAINS
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newEventWithString(name, referenceDate, firstdate, lastDate, interval, offset, errno) RESULT(ret_event) !OK-TESTED.
    TYPE(event), POINTER :: ret_event
    CHARACTER(len=*), INTENT(in)           :: name
    CHARACTER(len=*), INTENT(in)           :: referenceDate
    CHARACTER(len=*), INTENT(in)           :: firstDate
    CHARACTER(len=*), INTENT(in)           :: lastDate
    CHARACTER(len=*), INTENT(in)           :: interval
    CHARACTER(len=*), TARGET, INTENT(in), OPTIONAL :: offset
    CHARACTER(len=5), SAVE, TARGET :: zeroOffset = "PT00S" !Optional offset's proxy string.
    TYPE(c_ptr) :: c_pointer
    CHARACTER(len=:), POINTER :: offset_arg
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(offset)) THEN
      offset_arg => offset(VERIFY(offset, " "):LEN_TRIM(offset))
    ELSE
      offset_arg => zeroOffset
    END IF
    c_pointer = my_newevent(TRIM(ADJUSTL(name))//c_null_char,          &
         &                  TRIM(ADJUSTL(referenceDate))//c_null_char,      &
         &                  TRIM(ADJUSTL(firstDate))//c_null_char,          &
         &                  TRIM(ADJUSTL(lastDate))//c_null_char,           &
         &                  TRIM(ADJUSTL(interval))//c_null_char,           &
         &                  offset_arg//c_null_char)
    IF (PRESENT(errno)) errno = MERGE(0, 6*100 + 1, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_event)
  END FUNCTION newEventWithString
  !
  !! @param[out]       errno       optional, error message
  RECURSIVE FUNCTION newEventWithDataTypes(name, referenceDate, firstdate, lastDate, interval, offset, errno) RESULT(ret_event) !OK-TESTED.
    TYPE(event), POINTER :: ret_event
    CHARACTER(len=*), INTENT(in) :: name
    TYPE(datetime), POINTER, INTENT(in) :: referenceDate
    TYPE(datetime), POINTER, INTENT(in) :: firstDate
    TYPE(datetime), POINTER, INTENT(in) :: lastDate
    TYPE(timedelta), POINTER, INTENT(in) :: interval
    TYPE(timedelta), POINTER, OPTIONAL, INTENT(in) :: offset
    TYPE(c_ptr) :: c_pointer, offset_c
    INTEGER, OPTIONAL, INTENT(out) :: errno
    IF (PRESENT(offset)) THEN
      offset_c = C_LOC(offset)
    ELSE
      offset_c = c_null_ptr
    END IF
    c_pointer = my_neweventwithdatatypes(TRIM(ADJUSTL(name))//c_null_char, &
         &                  C_LOC(referenceDate), C_LOC(firstDate), C_LOC(lastDate), C_LOC(interval), offset_c)
    IF (PRESENT(errno)) errno = MERGE(0, 6*100 + 1, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_event)
  END FUNCTION newEventWithDataTypes
  !>
  !! @brief Destructor of Event.
  !!
  !! @param  my_event
  !!         A pointer to type event. my_event is deallocated.
  !!
  !! WARNING: If my_event was added to a group, this should never be called;
  !! use removeEventFromEventGroup instead.
  RECURSIVE SUBROUTINE deallocateEvent(my_event)
    TYPE(event), POINTER :: my_event
    CALL my_deallocateevent(C_LOC(my_event))
    my_event => NULL()
  END SUBROUTINE deallocateEvent
  !>
  RECURSIVE FUNCTION constructAndCopyEvent(my_event, errno) RESULT(ret_event) !OK-TESTED.
    TYPE(event), POINTER :: ret_event
    TYPE(event), TARGET, INTENT(in) :: my_event
    INTEGER, OPTIONAL, INTENT(out) :: errno
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_constructandcopyevent(C_LOC(my_event))
    IF (PRESENT(errno)) errno = MERGE(0, 6*100 + 1, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_event)
  END FUNCTION constructAndCopyEvent
  !>
  !! @brief Get Event as a string.
  !!
  !! @param[in]  my_event
  !!         A pointer to type event. The event to be converted to string.
  !!
  !! @param[out]  string
  !!         String where event is to be written.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE eventToString(my_event, string, errno) !TODO:C code still incomplete.
    TYPE(event), POINTER, INTENT(in) :: my_event
    CHARACTER(len=max_eventname_str_len), INTENT(out) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL, INTENT(out) :: errno
    dummy_ptr = my_eventtostring(C_LOC(my_event), string)
    IF (PRESENT(errno)) errno = MERGE(0, 6*100 + 3, C_ASSOCIATED(dummy_ptr))
    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE eventToString
  !>
  !! @brief Check if this event is active by comparing event's trigger time with current_dt.
  !!
  !!        The current_dt must lie in event's trigger time
  !!
  !!        (subject to optional specified slack: [Trigger_time -
  !!        minus_slack, Trigger_time + plus_slack]. Slacks can be
  !!        NULL. Always inclusive.)
  !!
  !!        The lib has no built-in clock but relies on
  !!        isCurrentEventActive(.) being called (polled) from the
  !!        application at fixed intervals.
  !!
  !! @param  my_event
  !!         A pointer to type event. This is the event being tested.
  !!
  !! @param  my_datetime
  !!         A variable of type datetime. This is the 'current' datetime of the system.
  !!
  !! @param  plus_slack
  !!         A pointer to type timedelta. Events are triggered between [actual_trigger_time, actual_trigger_time + plus_slack].
  !!         Sign MUST be '+'
  !!
  !! @param  minus_slack
  !!         A pointer to type timedelta. Events are triggered between [actual_trigger_time - minus_slack, actual_trigger_time].
  !!         Sign MUST be '+'
  !!
  !! @return ret
  !!        true/false indicating if the event is active.
  RECURSIVE FUNCTION isCurrentEventActive(my_event, my_datetime, plus_slack, minus_slack) RESULT(ret)
    TYPE(event), POINTER, INTENT(in) :: my_event
    TYPE(datetime), TARGET, INTENT(in) :: my_datetime
    TYPE(timedelta), POINTER, OPTIONAL, INTENT(in) :: plus_slack
    TYPE(timedelta), POINTER, OPTIONAL, INTENT(in) :: minus_slack
    TYPE(c_ptr) :: plus_slack_c, minus_slack_c
    LOGICAL(c_bool) :: ret
    IF (PRESENT(plus_slack)) THEN
      plus_slack_c = C_LOC(plus_slack)
    ELSE
      plus_slack_c = c_null_ptr
    END IF
    IF (PRESENT(minus_slack)) THEN
      minus_slack_c = C_LOC(minus_slack)
    ELSE
      minus_slack_c = c_null_ptr
    END IF
    ret = my_isCurrentEventActive(C_LOC(my_event), C_LOC(my_datetime), &
         &                        plus_slack_c, minus_slack_c)
  END FUNCTION isCurrentEventActive
  !>
  !! @brief Checks, if next event is on a new day
  !!
  !! @param[in] my_event
  !!        A pointer to a type event
  !!
  !! @returns ret
  !!        Logical: true if next event is on new day
  RECURSIVE FUNCTION isEventNextInNextDay(my_event) RESULT(ret)
    TYPE(event), POINTER, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_iseventnextinnextday(C_LOC(my_event))
  END FUNCTION isEventNextInNextDay
  !>
  !! @brief Checks, if next event is in a new month
  !!
  !! @param[in] my_event
  !!        A pointer to a type event
  !!
  !! @returns ret
  !!        Logical: true if next event is in a new month
  RECURSIVE FUNCTION iseventNextInNextMonth(my_event) RESULT(ret)
    TYPE(event), POINTER, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_iseventnextinnextmonth(C_LOC(my_event))
  END FUNCTION iseventNextInNextMonth
  !>
  !! @brief Checks, if next event is in a new year
  !!
  !! @param[in] my_event
  !!        A pointer to a type event
  !!
  !! @returns ret
  !!        Logical: true if next event is in a new year
  RECURSIVE FUNCTION iseventNextInNextYear(my_event) RESULT(ret)
    TYPE(event), POINTER, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_iseventnextinnextyear(C_LOC(my_event))
  END FUNCTION iseventNextInNextYear
  !>
  !! @brief Get the Datetime when 'this' event will be triggered next.
  !!
  !! WARNING: The value returned is with-respect-to current_dt and not
  !!          a true copy of triggerNextEventDateTime in the event
  !!          data structure.
  !!
  !! @param[in] my_event
  !!         A variable of type event. This is the event being queried.
  !!
  !! @param[in]  my_currentdatetime
  !!         A variable of type datetime. The next trigger datetime is copied here.
  !!
  !! @param[out] my_datetime
  !!         A variable of type datetime with next-trigger datetime.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE getTriggerNextEventAtDateTime(my_event, my_currentdatetime, my_datetime, errno)
    TYPE(event), TARGET, INTENT(in) :: my_event
    TYPE(datetime), TARGET, INTENT(in) :: my_currentdatetime
    TYPE(datetime), TARGET, INTENT(out) :: my_datetime
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL, INTENT(out) :: errno
    dummy_ptr = my_gettriggernexteventatdatetime(C_LOC(my_event), C_LOC(my_currentdatetime), C_LOC(my_datetime))
    IF (PRESENT(errno)) errno = MERGE(0, 6*100 + 8, C_ASSOCIATED(dummy_ptr))
  END SUBROUTINE getTriggerNextEventAtDateTime
  !>
  !! @brief Get the Datetime when 'this' event will be triggered last.
  !!
  !! NOTE: If the event was never tiggered, default value of
  !!       0-01-01T00:00:00.000 is returned.
  !!
  !! @param[in] my_event
  !!         A pointer to type event. This is the event being queried.
  !!
  !! @param[out] my_datetime
  !!         A variable of type datetime with last-trigger datetime.
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE getTriggeredPreviousEventAtDateTime(my_event, my_datetime, errno)
    TYPE(event), TARGET, INTENT(in) :: my_event
    TYPE(datetime), TARGET, INTENT(out) :: my_datetime
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_gettriggeredpreviouseventatdatetime(C_LOC(my_event), C_LOC(my_datetime))
    IF (PRESENT(errno)) errno = MERGE(0, 6*100 + 9, C_ASSOCIATED(dummy_ptr))
  END SUBROUTINE getTriggeredPreviousEventAtDateTime
  !>
  !! @brief get the event id
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret_evtid
  !!        the event id
  !!
  RECURSIVE FUNCTION getEventId(my_event) RESULT(ret_evtid) !OK-TESTED.
    INTEGER(c_int64_t) :: ret_evtid
    TYPE(event), POINTER :: my_event
    ret_evtid = my_event%eventId
  END FUNCTION getEventId
  !>
  !! @brief get the event name
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @param[out]       string      the name of the event
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE getEventName(my_event, string, errno) !OK-TESTED.
    TYPE(event), POINTER, INTENT(in) :: my_event
    CHARACTER(len=max_eventname_str_len), INTENT(out) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL :: errno
    dummy_ptr = my_geteventname(C_LOC(my_event), string)
    IF (PRESENT(errno)) errno = MERGE(0, 6*100 + 11, C_ASSOCIATED(dummy_ptr))
    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE getEventName
  !>
  !! @brief get the event reference date
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret_referenceDateTime
  !!        A pointer of type datetime. The event's reference date.
  !!
  RECURSIVE FUNCTION getEventReferenceDateTime(my_event) RESULT(ret_referenceDateTime) !OK-TESTED.
    TYPE(datetime), POINTER :: ret_referenceDateTime
    TYPE(event), POINTER, INTENT(in) :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%eventReferenceDatetime
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, ret_referenceDateTime)
    ELSE
      ret_referenceDateTime => NULL()
    END IF
  END FUNCTION getEventReferenceDateTime
  !>
  !! @brief get the event first date
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret_eventFirstDateTime
  !!        A pointer of type datetime. The event's first date.
  !!
  RECURSIVE FUNCTION getEventFirstDateTime(my_event) RESULT(ret_eventFirstDateTime) !OK-TESTED.
    TYPE(datetime), POINTER :: ret_eventFirstDateTime
    TYPE(event), POINTER, INTENT(in) :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%eventFirstDateTime
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, ret_eventFirstDateTime)
    ELSE
      ret_eventFirstDateTime => NULL()
    END IF
  END FUNCTION getEventFirstDateTime
  !>
  !! @brief get the event last date
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret_eventLastDateTime
  !!        A pointer of type datetime. The event's last date.
  !!
  RECURSIVE FUNCTION getEventLastDateTime(my_event) RESULT(ret_eventLastDateTime) !OK-TESTED.
    TYPE(datetime), POINTER :: ret_eventLastDateTime
    TYPE(event), POINTER, INTENT(in) :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%eventLastDateTime
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, ret_eventLastDateTime)
    ELSE
      ret_eventLastDateTime => NULL()
    END IF
  END FUNCTION getEventLastDateTime
  !>
  !! @brief get the event interval
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret_eventInterval
  !!        A pointer of type timedelta. The event's last date.
  !!
  RECURSIVE FUNCTION getEventInterval(my_event) RESULT(ret_eventInterval) !OK-TESTED.
    TYPE(timedelta), POINTER :: ret_eventInterval
    TYPE(event), POINTER, INTENT(in) :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%eventInterval
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, ret_eventInterval)
    ELSE
      ret_eventInterval => NULL()
    END IF
  END FUNCTION getEventInterval
  !>
  !! @brief Check if event is first
  !!
  !! @param[in] my_event
  !!        A reference of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is first
  !!
  RECURSIVE FUNCTION getNextEventIsFirst(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_getnexteventisfirst(C_LOC(my_event))
  END FUNCTION getNextEventIsFirst
  !>
  !! @brief Check if event is first in day
  !!
  !! @param[in] my_event
  !!        A reference of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is first in day
  !!
  RECURSIVE FUNCTION getEventisFirstInDay(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventisfirstinday(C_LOC(my_event))
  END FUNCTION getEventisFirstInDay
  !>
  !! @brief Check if event is first in month
  !!
  !! @param[in] my_event
  !!        A reference of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is first in month
  !!
  RECURSIVE FUNCTION getEventisFirstInMonth(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventisfirstinmonth(C_LOC(my_event))
  END FUNCTION getEventisFirstInMonth
  !>
  !! @brief Check if event is first in year
  !!
  !! @param[in] my_event
  !!        A reference of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is first in year
  !!
  RECURSIVE FUNCTION getEventisFirstInYear(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventisfirstinyear(C_LOC(my_event))
  END FUNCTION getEventisFirstInYear
  !>
  !! @brief Check if event is last in day
  !!
  !! @param[in] my_event
  !!        A reference of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is last in day
  !!

  RECURSIVE FUNCTION getEventisLastInDay(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventislastinday(C_LOC(my_event))
  END FUNCTION getEventisLastInDay
  !>
  !! @brief Check if event is last in month
  !!
  !! @param[in] my_event
  !!        A reference of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is last in month
  !!
  RECURSIVE FUNCTION getEventisLastInMonth(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventislastinmonth(C_LOC(my_event))
  END FUNCTION getEventisLastInMonth
  !>
  !! @brief Check if event is last in year
  !!
  !! @param[in] my_event
  !!        A reference of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is last in year
  !!
  RECURSIVE FUNCTION getEventisLastInYear(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventislastinyear(C_LOC(my_event))
  END FUNCTION getEventisLastInYear
  !
END MODULE mtime_events
!>
!! @brief Event-groups which contains a list of events.
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
MODULE mtime_eventgroups
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: &
    & c_associated, c_f_pointer, c_int64_t, c_loc, c_null_char, c_ptr
  USE mtime_constants, ONLY: max_groupname_str_len
  USE mtime_c_bindings, ONLY: &
    & eventgroup, my_addeventtoeventgroup, my_deallocateeventgroup, &
    & my_geteventgroupname, my_neweventgroup, my_removeeventfromeventgroup
  USE mtime_events
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: eventgroup
  PUBLIC :: newEventGroup
  PUBLIC :: deallocateEventGroup
  PUBLIC :: addEventToEventGroup
  PUBLIC :: removeEventFromEventGroup
  PUBLIC :: getEventGroupId
  PUBLIC :: getEventGroupName
  PUBLIC :: getFirstEventFromEventGroup
  PUBLIC :: getNextEventFromEventGroup
  !
CONTAINS
  !>
  !! @brief Construct new event-Group using a string.
  !!
  !! @param[in]  name
  !!         This string contains the name of event group.
  !!
  !! @param[out]       errno       optional, error message
  !!
  !! @return ret_eventgroup
  !!         A pointer to an initialized event-Group.
  !!
  RECURSIVE FUNCTION newEventGroup(name, errno) RESULT(ret_eventgroup) !OK-TESTED.
    TYPE(eventgroup), POINTER :: ret_eventgroup
    CHARACTER(len=*), INTENT(in) :: name
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL, INTENT(out):: errno
    c_pointer = my_neweventgroup(TRIM(ADJUSTL(name))//c_null_char)
    IF (PRESENT(errno)) errno = MERGE(0, 7*100 + 1, C_ASSOCIATED(c_pointer))
    CALL C_F_POINTER(c_pointer, ret_eventgroup)
  END FUNCTION newEventGroup
  !>
  !! @brief Destructor of EventGroup.
  !!
  !! @param[in] my_eventgroup
  !!         A pointer to type eventGroup. my_eventgroup is deallocated.
  !!
  RECURSIVE SUBROUTINE deallocateEventGroup(my_eventgroup) !OK-TESTED.
    TYPE(eventgroup), POINTER :: my_eventgroup
    CALL my_deallocateeventgroup(C_LOC(my_eventgroup))
    my_eventgroup => NULL()
  END SUBROUTINE deallocateEventGroup
  !>
  !! @brief Add new event to an eventgroup.
  !!
  !! @param  my_event
  !!         A reference to type event. The event to be added.
  !!
  !! @param  my_eventgroup
  !!         A reference to type eventgroup. The eventgroup where the event is added.
  !!
  !! @return ret
  !!         true/false indicating success or failure of addition.
  RECURSIVE FUNCTION addEventToEventGroup(my_event, my_eventgroup) RESULT(ret) !OK-TESTED.
    LOGICAL :: ret
    TYPE(event), TARGET, INTENT(in) :: my_event
    TYPE(eventgroup), TARGET, INTENT(inout) :: my_eventgroup
    ret = my_addeventtoeventgroup(C_LOC(my_event), C_LOC(my_eventgroup))
  END FUNCTION addEventToEventGroup
  !>
  !! @brief Remove event from eventgroup. CRITICAL: Also, deallocate the event.
  !!
  !! @param  my_name
  !!         The name of event to be removed.
  !!
  !! @param  my_eventgroup
  !!         A pointer to  type eventgroup. The eventgroup to which this event belongs.
  !!
  !! @return ret
  !!         true/false indicating success or failure of removal.
  RECURSIVE FUNCTION removeEventfromEventGroup(my_name, my_eventgroup) RESULT(ret) !OK-TESTED.
    LOGICAL :: ret
    CHARACTER(len=*), INTENT(in) :: my_name
    TYPE(eventgroup), POINTER, INTENT(inout) :: my_eventgroup
    ret = my_removeeventfromeventgroup(TRIM(ADJUSTL(my_name))//c_null_char, C_LOC(my_eventgroup))
  END FUNCTION removeEventFromEventGroup
  !>
  !! @brief Get event group id
  !!
  !! @param[in]  my_eventgroup
  !!         A pointer to  type eventgroup.
  !!
  !! @return ret_grpid
  !!         The event group id
  RECURSIVE FUNCTION getEventGroupId(my_eventgroup) RESULT(ret_grpid) !OK-TESTED.
    INTEGER(c_int64_t) :: ret_grpid
    TYPE(eventgroup), POINTER, INTENT(in) :: my_eventgroup
    ret_grpid = my_eventgroup%eventGroupId
  END FUNCTION getEventGroupId
  !>
  !! @brief get the event group name
  !!
  !! @param[in] my_eventgroup
  !!        A pointer of type event.
  !!
  !! @param[out]       string      the name of the event group
  !!
  !! @param[out]       errno       optional, error message
  RECURSIVE SUBROUTINE getEventGroupName(my_eventgroup, string, errno)  !TESTED-OK.
    TYPE(eventgroup), POINTER, INTENT(in) :: my_eventgroup
    CHARACTER(len=max_groupname_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    dummy_ptr = my_geteventgroupname(C_LOC(my_eventgroup), string)
    IF (PRESENT(errno)) errno = MERGE(0, 7*100 + 6, C_ASSOCIATED(dummy_ptr))
    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE getEventGroupName
  !>
  !! @brief get the first event in event group
  !!
  !! @param[in] my_eventgroup
  !!        A pointer of type eventgroup.
  !!
  !! @returns ret_event
  !!        A pointer of type event. The first event in eventgroup
  RECURSIVE FUNCTION getFirstEventFromEventGroup(my_eventgroup) RESULT(ret_event) !OK-TESTED.
    TYPE(event), POINTER :: ret_event
    TYPE(eventgroup), POINTER, INTENT(in) :: my_eventgroup
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_eventgroup%firstEventInGroup
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, ret_event)
    ELSE
      ret_event => NULL()
    END IF
  END FUNCTION getFirstEventFromEventGroup
  !>
  !! @brief get the next event in an event group an event belongs to
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret_event
  !!        A pointer of type event. The next event in eventgroup
  RECURSIVE FUNCTION getNextEventFromEventGroup(my_event) RESULT(ret_event) !OK-TESTED.
    TYPE(event), POINTER :: ret_event
    TYPE(event), POINTER, INTENT(in) :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%nextEventInGroup
    IF (C_ASSOCIATED(c_pointer)) THEN
      CALL C_F_POINTER(c_pointer, ret_event)
    ELSE
      ret_event => NULL()
    END IF
  END FUNCTION getNextEventFromEventGroup
  !
END MODULE mtime_eventgroups
!>
!! @brief Support for handling ISO 8601:2004 repeated time interval strings.
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
MODULE mtime_utilities
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_null_char
  USE mtime_c_bindings, ONLY: my_getrepetitions, my_splitrepetitionstring
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: splitRepetitionString
  PUBLIC :: getRepetitions
  PUBLIC :: max_repetition_str_len
  !
  !> provides a string length for the maximum handable repetition string length input
  INTEGER, PARAMETER :: max_repetition_str_len = 132
  !
CONTAINS
  !>
  !! @brief Extract number of repetitions from repetition string part.
  !!
  !! @param[in] repetitionString
  !!         A repetition string starting with 'R'.
  !!         A string literal can be accepted.
  !!
  !! @return r
  !!         An int representing the number of repetitions.
  !!
  RECURSIVE FUNCTION getRepetitions(repetitionString) RESULT(r)
    INTEGER :: r
    CHARACTER(len=*), INTENT(in) :: repetitionString
    r = my_getRepetitions(TRIM(ADJUSTL(repetitionString))//c_null_char)
  END FUNCTION getRepetitions
  !>
  !! @brief Split ISO 8601:2004 repeated time interval strings into base components
  !!
  !! @param[in] recurringTimeInterval
  !!         The string should contain an  ISO 8601:2004 repeated time
  !!         interval string.
  !!         A string literal can be accepted.
  !! @param[out] repetitor
  !!         Contains the repetitor part of the input string.
  !! @param[out] start
  !!         Contains the start date part of the input string.
  !! @param[out] end
  !!         Contains the end date part of the input string.
  !! @param[out] duration
  !!         Contains the duration part of the input string.
  !! @param[out] lrepetitor
  !!         Logical: true, if repetion is available
  !! @param[out] lstart
  !!         Logical: true, if start is available
  !! @param[out] lend
  !!         Logical: true, if end is available
  !! @param[out] lduration
  !!         Logical: true, if duration is available
  !!
  RECURSIVE SUBROUTINE splitRepetitionString(recurringTimeInterval, &
                                           & repetitor, start, END, duration, &
                                           & lrepetitor, lstart, lend, lduration)
    CHARACTER(len=*), INTENT(in) :: recurringTimeInterval
    CHARACTER(len=*), INTENT(out) :: repetitor
    CHARACTER(len=*), INTENT(out) :: start
    CHARACTER(len=*), INTENT(out) :: END
    CHARACTER(len=*), INTENT(out) :: duration

    LOGICAL, INTENT(out) :: lrepetitor
    LOGICAL, INTENT(out) :: lstart
    LOGICAL, INTENT(out) :: lend
    LOGICAL, INTENT(out) :: lduration

    INTEGER :: i

    CALL my_splitRepetitionString(TRIM(ADJUSTL(recurringTimeInterval))//c_null_char, repetitor, start, END, duration)

    lrepetitor = repetitor(1:1) /= c_null_char
    char_loop1: DO i = 1, LEN(repetitor)
      IF (repetitor(i:i) == c_null_char) EXIT char_loop1
    END DO char_loop1
    repetitor(i:LEN(repetitor)) = ' '

    lstart = start(1:1) /= c_null_char
    char_loop2: DO i = 1, LEN(start)
      IF (start(i:i) == c_null_char) EXIT char_loop2
    END DO char_loop2
    start(i:LEN(start)) = ' '

    lend = END(1:1) /= c_null_char
    char_loop3: DO i = 1, LEN(END)
      IF (END(i:i) == c_null_char) EXIT char_loop3
    END DO char_loop3
    END(i:LEN(END)) = ' '

    lduration = duration(1:1) /= c_null_char
    char_loop4: DO i = 1, LEN(duration)
      IF (duration(i:i) == c_null_char) EXIT char_loop4
    END DO char_loop4
    duration(i:LEN(duration)) = ' '

  END SUBROUTINE splitRepetitionString
!
END MODULE mtime_utilities
!>
!! @}
!>
!! @brief mtime is a compound module making all library components accessible via one module.
!!
!! @details
!!
!! @author  Luis Kornblueh, Max Planck Institute for Meteorology
!! @author  Rahul Sinha, Max Planck Institute for Meteorology
!!
!___________________________________________________________________________________________________________
MODULE mtime
  !
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  !
  USE mtime_constants
  USE mtime_error_handling
  USE mtime_c_bindings, ONLY: &
    & calendarType, divisionquotienttimespan, resetCalendar, setCalendar
  USE mtime_calendar
  USE mtime_juliandelta
  USE mtime_julianday
  USE mtime_date
  USE mtime_time
  USE mtime_datetime
  USE mtime_timedelta
  USE mtime_events
  USE mtime_eventgroups
  USE mtime_utilities
  !
  IMPLICIT NONE
  !
  PUBLIC
  !
#ifdef DOXYGEN_DOCUMENTATION_ONLY
  INTEGER, PARAMETER :: no_of_sec_in_a_day = 86400  !!< number of seconds per day, defined in C
  INTEGER, PARAMETER :: no_of_sec_in_a_hour = 3600  !!< number of seconds per hour, defined in C
  INTEGER, PARAMETER :: no_of_sec_in_a_minute = 60  !!< number of seconds per minute, defined in C
  !
  INTEGER, PARAMETER :: no_of_ms_in_a_day = 86400000  !!< number of milli-seconds per day, defined in C
  INTEGER, PARAMETER :: no_of_ms_in_half_day = 43200000  !!< number of milli-seconds per 12 hours, defined in C
  INTEGER, PARAMETER :: no_of_ms_in_a_hour = 3600000  !!< number of milli-seconds per hour, defined in C
  INTEGER, PARAMETER :: no_of_ms_in_a_minute = 60000  !!< number of milli-seconds per minute, defined in C
  INTEGER, PARAMETER :: no_of_ms_in_a_second = 1000  !!< number of milli-seconds per second, defined in C
#endif
  !
  !> @cond DOXYGEN_IGNORE_THIS
  INTEGER(c_int), BIND(c, name='NO_OF_SEC_IN_A_DAY') :: no_of_sec_in_a_day
  INTEGER(c_int), BIND(c, name='NO_OF_SEC_IN_A_HOUR') :: no_of_sec_in_a_hour
  INTEGER(c_int), BIND(c, name='NO_OF_SEC_IN_A_MINUTE') :: no_of_sec_in_a_minute
  !
  INTEGER(c_int), BIND(c, name='NO_OF_MS_IN_A_DAY') :: no_of_ms_in_a_day
  INTEGER(c_int), BIND(c, name='NO_OF_MS_IN_HALF_DAY') :: no_of_ms_in_half_day
  INTEGER(c_int), BIND(c, name='NO_OF_MS_IN_A_HOUR') :: no_of_ms_in_a_hour
  INTEGER(c_int), BIND(c, name='NO_OF_MS_IN_A_MINUTE') :: no_of_ms_in_a_minute
  INTEGER(c_int), BIND(c, name='NO_OF_MS_IN_A_SECOND') :: no_of_ms_in_a_second
  !> @endcond DOXYGEN_IGNORE_THIS

END MODULE mtime
