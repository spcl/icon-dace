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
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_char, c_null_char, c_ptr, c_associated, c_int64_t
  USE mtime_c_bindings
  USE mtime_error_handling
  USE mtime_constants
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  ENUM, BIND(c)
    ENUMERATOR :: calendar_not_set    = 0   ! calendar is not defined yet
    ENUMERATOR :: proleptic_gregorian = 1   ! proleptic Gregorian calendar
    ENUMERATOR :: year_of_365_days    = 2   ! 365 day year without leap years
    ENUMERATOR :: year_of_360_days    = 3   ! 360 day year with 30 day months
  END ENUM
  !
  PUBLIC :: calendar_not_set, proleptic_gregorian, year_of_365_days, year_of_360_days
  PUBLIC :: setCalendar
  PUBLIC :: resetCalendar
  PUBLIC :: calendarType
  PUBLIC :: calendarToString
  !
  !
CONTAINS
  !
  !>
  !! @brief convert the calendar identifier into a human readable string
  !!
  !! @param[out]       string      the calendar type verbose
  !! @param[out]       errno       optional, error message
  !!
  SUBROUTINE calendarToString(string, errno) !TESTED-OK
    CHARACTER(len=max_calendar_str_len), INTENT(out) :: string
    INTEGER :: i
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_calendartostring(string)
    IF (.NOT. c_ASSOCIATED(c_pointer)) THEN
      IF (PRESENT(errno)) errno = 0 * 100 + 4
      string = '<null>'
    ELSE
      char_loop: DO i = 1 , LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    ENDIF
  END SUBROUTINE calendarToString
  !
END MODULE mtime_calendar


MODULE mtime_juliandelta
  !
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t, c_char, c_ptr, c_loc, &
       &                                 c_associated, c_f_pointer
  USE mtime_c_bindings
  USE mtime_error_handling
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
  FUNCTION newJuliandelta(sign,day, ms, errno) RESULT(ret_juliandelta) !OK-TESTED.
    TYPE(juliandelta), POINTER :: ret_juliandelta
    CHARACTER(c_char),  VALUE, INTENT(in) :: sign
    INTEGER(c_int64_t), INTENT(in) :: day
    INTEGER(c_int64_t), INTENT(in) :: ms
    TYPE(c_ptr)                    :: c_pointer
    INTEGER, OPTIONAL              :: errno

    IF (PRESENT(errno)) errno = 0
    c_pointer = my_newjuliandelta(sign,day, ms)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 1 * 100 + 1
    CALL c_f_POINTER(c_pointer, ret_juliandelta)
  END FUNCTION newJuliandelta
  !
  !>
  !! @brief destructor for a Julian delta
  SUBROUTINE deallocateJuliandelta(my_juliandelta) !OK-TESTED.
    TYPE(juliandelta), POINTER :: my_juliandelta
    CALL my_deallocatejuliandelta(c_LOC(my_juliandelta))
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
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t, c_char, c_null_char, c_ptr, &
       &                                 c_loc, c_f_pointer, c_associated
  USE mtime_c_bindings
  USE mtime_error_handling
  USE mtime_constants
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
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
  INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE replacejulianday
  END INTERFACE ASSIGNMENT (=)
  !
  INTERFACE OPERATOR (+)
    MODULE PROCEDURE addjuliandeltatojulianday
    MODULE PROCEDURE addjuliandaytojuliandelta
  END INTERFACE OPERATOR (+)
  !
  INTERFACE OPERATOR (-)
    MODULE PROCEDURE substractjuliandayfromjulianday
  END INTERFACE OPERATOR (-)
  !
  INTERFACE OPERATOR (>)
    MODULE PROCEDURE julianday_gt
  END INTERFACE OPERATOR (>)
  !
  INTERFACE OPERATOR (<)
    MODULE PROCEDURE julianday_lt
  END INTERFACE OPERATOR (<)
  !
  INTERFACE OPERATOR (<=)
    MODULE PROCEDURE julianday_lt_or_eq
  END INTERFACE OPERATOR (<=)
  !
  INTERFACE OPERATOR (>=)
    MODULE PROCEDURE julianday_gt_or_eq
  END INTERFACE OPERATOR (>=)
  !
  INTERFACE OPERATOR (==)
    MODULE PROCEDURE julianday_eq
  END INTERFACE OPERATOR (==)
  !
  INTERFACE OPERATOR (/=)
    MODULE PROCEDURE julianday_ne
  END INTERFACE OPERATOR (/=)
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
  FUNCTION newJulianday(day, ms, errno) RESULT(ret_julianday) !OK-TESTED.
    TYPE(julianday), POINTER :: ret_julianday
    INTEGER(c_int64_t), INTENT(in) :: day
    INTEGER(c_int64_t), INTENT(in) :: ms
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_newjulianday(day, ms)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 1 * 100 + 1
    CALL c_f_POINTER(c_pointer, ret_julianday)
  END FUNCTION newJulianday
  !
  !>
  !! @brief destructor for a Julian date
  !!
  !! @param     my_julianday   a pointer of type(julianday)
  SUBROUTINE deallocateJulianday(my_julianday) !OK-TESTED.
    TYPE(julianday), POINTER :: my_julianday
    CALL my_deallocatejulianday(c_LOC(my_julianday))
    my_julianday => NULL()
  END SUBROUTINE deallocateJulianday
  !
  ! NOTE: Do not call the function using replacejulianday(.,.) directly; Use overloaded '=' instead as in "dest = src"
  SUBROUTINE replacejulianday(dest, src) !OK-TESTED.
    TYPE(julianday), TARGET, INTENT(inout) :: dest
    TYPE(julianday), TARGET, INTENT(in) :: src
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_replacejulianday(c_LOC(src), c_LOC(dest))
  END SUBROUTINE replacejulianday
  !
  !>
  !! @brief get Julian day as a string.
  !!
  !! @param[in]  my_julianday   a pointer to type(julianday). The Julian day to be converted to a string
  !! @param[out] string         the Julian day verbose
  !! @param[out]       errno       optional, error message
  SUBROUTINE juliandayToString(my_julianday, string, errno) !OK-TESTED.
    TYPE(julianday), POINTER :: my_julianday
    CHARACTER(len=max_julianday_str_len), INTENT(out) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_juliandaytostring(c_LOC(my_julianday), string)
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 1 * 100 + 3
    char_loop: DO i = 1 , LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE juliandayToString
  !
  FUNCTION addJuliandeltaToJulianday(op1, op2) RESULT(ret)
    TYPE(julianday), TARGET :: ret
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(juliandelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_addjuliandeltatojulianday(c_LOC(op1), c_LOC(op2), c_LOC(ret))
  END FUNCTION addJuliandeltaToJulianday
  !
  FUNCTION addJuliandayToJuliandelta(op2, op1) RESULT(ret)
    TYPE(julianday), TARGET :: ret
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(juliandelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_addjuliandeltatojulianday(c_LOC(op1), c_LOC(op2), c_LOC(ret))
  END FUNCTION addJuliandayToJuliandelta
  !
  FUNCTION substractJuliandayFromJulianday(op1, op2) RESULT(ret)
    TYPE(juliandelta), TARGET :: ret
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_substractjulianday(c_LOC(op1), c_LOC(op2), c_LOC(ret))
  END FUNCTION substractJuliandayFromJulianday
  !
  FUNCTION julianday_gt(op1, op2) RESULT(gt)
    LOGICAL :: gt
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(c_LOC(op1), c_LOC(op2))
    IF (ret == 1) THEN
      gt = .TRUE.
    ELSE
      gt = .FALSE.
    ENDIF
  END FUNCTION julianday_gt
  !
  FUNCTION julianday_lt(op1, op2) RESULT(lt)
    LOGICAL :: lt
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(c_LOC(op1), c_LOC(op2))
    IF (ret == -1) THEN
      lt = .TRUE.
    ELSE
      lt = .FALSE.
    ENDIF
  END FUNCTION julianday_lt
  !
  FUNCTION julianday_lt_or_eq(op1, op2) RESULT(lt_or_eq)
    LOGICAL :: lt_or_eq
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(c_LOC(op1), c_LOC(op2))
    IF ((ret == 0) .OR. (ret == -1)) THEN
      lt_or_eq = .TRUE.
    ELSE
      lt_or_eq = .FALSE.
    ENDIF
  END FUNCTION julianday_lt_or_eq
  !
  FUNCTION julianday_gt_or_eq(op1, op2) RESULT(gt_or_eq)
    LOGICAL :: gt_or_eq
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(c_LOC(op1), c_LOC(op2))
    IF ((ret == 0) .OR. (ret == 1)) THEN
      gt_or_eq = .TRUE.
    ELSE
      gt_or_eq = .FALSE.
    ENDIF
  END FUNCTION julianday_gt_or_eq
  !
  FUNCTION julianday_eq(op1, op2) RESULT(eq)
    LOGICAL :: eq
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(c_LOC(op1), c_LOC(op2))
    IF (ret == 0) THEN
      eq = .TRUE.
    ELSE
      eq = .FALSE.
    ENDIF
  END FUNCTION julianday_eq
  !
  FUNCTION julianday_ne(op1, op2) RESULT(ne)
    LOGICAL :: ne
    TYPE(julianday), TARGET, INTENT(in) :: op1
    TYPE(julianday), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparejulianday(c_LOC(op1), c_LOC(op2))
    IF (ret /= 0) THEN
      ne = .TRUE.
    ELSE
      ne = .FALSE.
    ENDIF
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
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_int64_t, c_char, c_ptr, c_null_char, &
                                       & c_loc, c_f_pointer, c_associated
  USE mtime_c_bindings
  USE mtime_constants
  USE mtime_error_handling
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: newdate
  PUBLIC :: deallocateDate
  PUBLIC :: replaceDate
  PUBLIC :: dateToString
  PUBLIC :: dateToPosixString
  !
  !
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
  FUNCTION newdatefromstring(string,errno) RESULT(ret_date)  !OK-TESTED.
    TYPE(date), POINTER :: ret_date
    CHARACTER(len=*), INTENT(in) :: string
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL :: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_newdatefromstring(TRIM(ADJUSTL(string))//c_null_char)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 2 * 100 + 1
    CALL c_f_POINTER(c_pointer, ret_date)
  END FUNCTION newdatefromstring
  !>
  !! @brief construct a new date from raw date
  !!
  !! @param[in] year           the year
  !! @param[in] month          the month
  !! @param[in] day            the day
  !! @param[out]       errno       optional, error message
  !! @return    ret_date       a pointer of type(date)
  FUNCTION newdatefromraw_yi8(year, month, day, errno) RESULT(ret_date) !OK-TESTED.
    TYPE(date), POINTER :: ret_date
    INTEGER(c_int64_t), INTENT(in) :: year
    INTEGER(c_int), INTENT(in) :: month, day
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_newrawdate(year, month, day)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 2 * 100 + 2
    CALL c_f_POINTER(c_pointer, ret_date)
  END FUNCTION newdatefromraw_yi8
  !>
  !! @brief construct a new date from raw date
  !!
  !! @param[in] year           the year
  !! @param[in] month          the month
  !! @param[in] day            the day
  !! @param[out]       errno       optional, error message
  !! @return    ret_date       a pointer of type(date)
  FUNCTION newdatefromraw(year, month, day, errno) RESULT(ret_date) !OK-TESTED.
    TYPE(date), POINTER :: ret_date
    INTEGER(c_int), INTENT(in) :: year, month, day
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_newrawdate(INT(year,c_int64_t), month, day)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 2 * 100 + 3
    CALL c_f_POINTER(c_pointer, ret_date)
  END FUNCTION newdatefromraw
  !>
  !! @brief construct a new date from an existing by construct and copy
  !!
  !! @param[in] src            a pointer of type(date)
  !! @param[out]       errno       optional, error message
  !! @return    ret_date       a pointer of type(date)
  FUNCTION newdatefromconstructandcopy(src, errno) RESULT(dest) !OK-TESTED
    TYPE(date), POINTER :: dest
    TYPE(date), TARGET :: src
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_constructandcopydate(c_LOC(src))
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 2 * 100 + 4
    CALL c_f_POINTER(c_pointer, dest)
  END FUNCTION newdatefromconstructandcopy
  !>
  !! @brief destructor for a date
  !!
  !! @param[in] my_date        a pointer of type(date)
  SUBROUTINE deallocateDate(my_date) !OK-TESTED.
    TYPE(date), POINTER :: my_date
    CALL my_deallocatedate(c_LOC(my_date))
    my_date => NULL()
  END SUBROUTINE deallocateDate
  !>
  !! @brief repace an existing date by a given one
  !!
  !! @param[in]  src            a pointer of type(date)
  !! @param[out] dest           a pointer of type(date)
  !! @param[out] errno          optional, error message
  SUBROUTINE replaceDate(dest, src, errno)  !OK-TESTED.
    TYPE(date), TARGET, INTENT(inout) :: dest
    TYPE(date), TARGET, INTENT(in) :: src
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_replacedate(c_LOC(src), c_LOC(dest))
    IF ((.NOT. c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 2 * 100 + 6
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
  SUBROUTINE dateToString(my_date, string, errno) !OK-TESTED.
    TYPE(date), POINTER :: my_date
    CHARACTER(len=max_date_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_datetostring(c_LOC(my_date), string)
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 2 * 100 + 7
    char_loop: DO i = 1 , LEN(string)
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
  SUBROUTINE dateToPosixString(my_date, string, fmtstr, errno) !OK-TESTED.
    TYPE(date), POINTER :: my_date
    CHARACTER(len=max_date_str_len) :: string
    CHARACTER(len=*) :: fmtstr
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_datetoposixstring(c_LOC(my_date), string, fmtstr)
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 2 * 100 + 8
    char_loop: DO i = 1 , LEN(string)
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
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_char, c_ptr, c_null_char, &
       &                                 c_loc, c_f_pointer, c_associated
  USE mtime_constants
  USE mtime_c_bindings
  USE mtime_error_handling
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: newtime
  PUBLIC :: deallocateTime
  PUBLIC :: replaceTime
  PUBLIC :: timeToString
  PUBLIC :: timeToPosixString
  !
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
  FUNCTION newtimefromstring(string, errno) RESULT(ret_time) !OK-TESTED.
    TYPE(time), POINTER :: ret_time
    CHARACTER(len=*), INTENT(in) :: string
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_newtimefromstring(TRIM(ADJUSTL(string))//c_null_char)
    IF ((.NOT.c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 3 * 100 + 1
    CALL c_f_POINTER(c_pointer, ret_time)
  END FUNCTION newtimefromstring
  !
  !! @param[out]       errno       optional, error message
  FUNCTION newtimefromraw(hour, minute, second, ms, errno) RESULT(ret_time) !OK-TESTED.
    TYPE(time), POINTER :: ret_time
    INTEGER(c_int), INTENT(in) :: hour, minute, second, ms
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_newrawtime(hour, minute, second, ms)
    IF ((.NOT.c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 3 * 100 + 2
    CALL c_f_POINTER(c_pointer, ret_time)
  END FUNCTION newtimefromraw
  !
  !! @param[out]       errno       optional, error message
  FUNCTION newtimefromconstructandcopy(src, errno) RESULT(dest) !OK-TESTED.
    TYPE(time), POINTER :: dest
    TYPE(time), TARGET :: src
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_constructandcopytime(c_LOC(src))
    IF ((.NOT.c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 3 * 100 + 3
    CALL c_f_POINTER(c_pointer, dest)
  END FUNCTION newtimefromconstructandcopy
  !>
  !! @brief Destructor of Time.
  !!
  !! @param  my_time
  !!         A pointer to type time. my_time is deallocated.
  SUBROUTINE deallocateTime(my_time) !OK-TESTED.
    TYPE(time), POINTER :: my_time
    CALL my_deallocatetime(c_LOC(my_time))
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
  SUBROUTINE replaceTime(dest,src,errno) !OK-TESTED.
    TYPE(time), TARGET, INTENT(in) :: src
    TYPE(time), TARGET, INTENT(inout) :: dest
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_replacetime(c_LOC(src), c_LOC(dest))
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 3 * 100 + 5
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
  SUBROUTINE timeToString(my_time, string, errno) !OK-TESTED.
    TYPE(time), POINTER :: my_time
    CHARACTER(len=max_time_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_timetostring(c_LOC(my_time), string)
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 3 * 100 + 6
    char_loop: DO i = 1 , LEN(string)
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
  SUBROUTINE timeToPosixString(my_time, string, fmtstr, errno) !OK-TESTED.
    TYPE(time), POINTER :: my_time
    CHARACTER(len=max_time_str_len) :: string
    CHARACTER(len=32) :: fmtstr
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_timetoposixstring(c_LOC(my_time), string, fmtstr)
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 3 * 100 + 7
    char_loop: DO i = 1 , LEN(string)
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
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_int64_t, c_char, c_null_char, c_ptr, &
       &                                 c_loc, c_f_pointer, c_associated
  USE mtime_constants
  USE mtime_c_bindings
  USE mtime_error_handling
  !
  USE mtime_julianday
  USE mtime_date
  USE mtime_time
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
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
  INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE replacedatetime
  END INTERFACE ASSIGNMENT (=)
  !
  INTERFACE OPERATOR (>)
    MODULE PROCEDURE datetime_gt
  END INTERFACE OPERATOR (>)
  !
  INTERFACE OPERATOR (<)
    MODULE PROCEDURE datetime_lt
  END INTERFACE OPERATOR (<)
  !
  INTERFACE OPERATOR (<=)
    MODULE PROCEDURE datetime_lt_or_eq
  END INTERFACE OPERATOR (<=)
  !
  INTERFACE OPERATOR (>=)
    MODULE PROCEDURE datetime_gt_or_eq
  END INTERFACE OPERATOR (>=)
  !
  INTERFACE OPERATOR (==)
    MODULE PROCEDURE datetime_eq
  END INTERFACE OPERATOR (==)
  !
  INTERFACE OPERATOR (/=)
    MODULE PROCEDURE datetime_ne
  END INTERFACE OPERATOR (/=)
  !
CONTAINS
  !
  !! @param[out]       errno       optional, error message
  FUNCTION newdatetimefromstring(string, errno) RESULT(ret_datetime) !OK-TESTED
    TYPE(datetime), POINTER :: ret_datetime
    CHARACTER(len=*), INTENT(in) :: string
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_newdatetime(TRIM(ADJUSTL(string))//c_null_char)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 4 * 100 + 1
    CALL c_f_POINTER(c_pointer, ret_datetime)
  END FUNCTION newdatetimefromstring
  !
  !! @param[out]       errno       optional, error message
  FUNCTION newdatetimefromraw_yi8(year, month, day, hour, minute, second, ms, errno) RESULT(ret_datetime) !OK-TESTED
    TYPE(datetime), POINTER :: ret_datetime
    INTEGER(c_int64_t), INTENT(in) :: year
    INTEGER(c_int), INTENT(in) :: month, day, hour, minute, second, ms
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_newrawdatetime(year, month, day, hour, minute, second, ms)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 4 * 100 + 2
    CALL c_f_POINTER(c_pointer, ret_datetime)
  END FUNCTION newdatetimefromraw_yi8
  !
  !! @param[out]       errno       optional, error message
  FUNCTION newdatetimefromraw(year, month, day, hour, minute, second, ms, errno) RESULT(ret_datetime) !OK-TESTED
    TYPE(datetime), POINTER :: ret_datetime
    INTEGER(c_int), INTENT(in) :: year, month, day, hour, minute, second, ms
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_newrawdatetime(INT(year,c_int64_t), month, day, hour, minute, second, ms)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 4 * 100 + 3
    CALL c_f_POINTER(c_pointer, ret_datetime)
  END FUNCTION newdatetimefromraw
  !
  !! @param[out]       errno       optional, error message
  FUNCTION newdatetimefromconstructandcopy(src, errno) RESULT(dest) !OK-TESTED.
    TYPE(datetime), POINTER :: dest
    TYPE(datetime), TARGET :: src
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_constructandcopydatetime(c_LOC(src))
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 4 * 100 + 4
    CALL c_f_POINTER(c_pointer, dest)
  END FUNCTION newdatetimefromconstructandcopy


  ! @return Minimum of two date time. In case of equality we return @p a.
  !
  ! @note This function does not return a copy of one of the arguments
  ! but returns only the corresponding result point to avoid
  ! unnecessary deallocate calls.
  FUNCTION datetime_min(a,b) RESULT(res)
    TYPE(datetime), POINTER :: res
    TYPE(datetime), POINTER :: a,b

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
  FUNCTION datetime_max(a,b) RESULT(res)
    TYPE(datetime), POINTER :: res
    TYPE(datetime), POINTER :: a,b

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
  SUBROUTINE deallocateDatetime(my_datetime) !OK-TESTED.
    TYPE(datetime), POINTER :: my_datetime
    CALL my_deallocatedatetime(c_LOC(my_datetime))
    NULLIFY(my_datetime)
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
  SUBROUTINE datetimeToString(my_datetime, string, errno) !OK-TESTED
    TYPE(datetime), TARGET, INTENT(in) :: my_datetime
    CHARACTER(len=max_datetime_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_datetimetostring(c_LOC(my_datetime), string)
    IF (.NOT. c_ASSOCIATED(dummy_ptr)) THEN
      IF (PRESENT(errno)) errno = 4 * 100 + 6
      string ='<null>'
    ELSE
      char_loop: DO i = 1 , LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    ENDIF
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
  SUBROUTINE datetimeToPosixString(my_datetime, string, fmtstr, errno) !OK-TESTED.
    TYPE(datetime), TARGET, INTENT(in) :: my_datetime
    CHARACTER(len=max_datetime_str_len) :: string
    CHARACTER(len=*) :: fmtstr
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_datetimetoposixstring(c_LOC(my_datetime), string, fmtstr)
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 4 * 100 + 7
    char_loop: DO i = 1 , LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '
  END SUBROUTINE datetimeToPosixString
  !
  ! NOTE: Do not call the function using replacedatetime(.,.) directly; Use overloaded '=' instead as in "dest = src"
  SUBROUTINE replacedatetime(dest, src) !OK-TESTED.
    TYPE(datetime), TARGET, INTENT(inout) :: dest
    TYPE(datetime), TARGET, INTENT(in) :: src
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_replacedatetime(c_LOC(src), c_LOC(dest))
  END SUBROUTINE replacedatetime
  !
  FUNCTION datetime_gt(op1, op2) RESULT(gt) !OK-TESTED.
    LOGICAL :: gt
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(c_LOC(op1), c_LOC(op2))
    IF (ret == 1) THEN
      gt = .TRUE.
    ELSE
      gt = .FALSE.
    ENDIF
  END FUNCTION datetime_gt
  !
  FUNCTION datetime_lt(op1, op2) RESULT(lt) !OK-TESTED.
    LOGICAL :: lt
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(c_LOC(op1), c_LOC(op2))
    IF (ret == -1) THEN
      lt = .TRUE.
    ELSE
      lt = .FALSE.
    ENDIF
  END FUNCTION datetime_lt
  !
  FUNCTION datetime_lt_or_eq(op1, op2) RESULT(lt_or_eq) !OK-TESTED.
    LOGICAL :: lt_or_eq
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(c_LOC(op1), c_LOC(op2))
    IF ((ret == 0) .OR. (ret == -1)) THEN
      lt_or_eq = .TRUE.
    ELSE
      lt_or_eq = .FALSE.
    ENDIF
  END FUNCTION datetime_lt_or_eq
  !
  FUNCTION datetime_gt_or_eq(op1, op2) RESULT(gt_or_eq) !OK-TESTED
    LOGICAL :: gt_or_eq
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(c_LOC(op1), c_LOC(op2))
    IF ((ret == 0) .OR. (ret == 1)) THEN
      gt_or_eq = .TRUE.
    ELSE
      gt_or_eq = .FALSE.
    ENDIF
  END FUNCTION datetime_gt_or_eq
  !
  FUNCTION datetime_eq(op1, op2) RESULT(eq) !OK-TESTED.
    LOGICAL :: eq
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(c_LOC(op1), c_LOC(op2))
    IF (ret == 0) THEN
      eq = .TRUE.
    ELSE
      eq = .FALSE.
    ENDIF
  END FUNCTION datetime_eq
  !
  FUNCTION datetime_ne(op1, op2) RESULT(ne) !OK-TESTED.
    LOGICAL :: ne
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(datetime), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparedatetime(c_LOC(op1), c_LOC(op2))
    IF (ret /= 0) THEN
      ne = .TRUE.
    ELSE
      ne = .FALSE.
    ENDIF
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
  FUNCTION getNoOfDaysInMonthDateTime(dt, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    INTEGER(c_int) :: getNoOfDaysInMonthDateTime
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    getNoOfDaysInMonthDateTime = my_getnoofdaysinmonthdatetime(c_LOC(dt))
    IF (getNoOfDaysInMonthDateTime == 0 .AND. PRESENT(errno)) errno = 4 * 100 + 15
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
  FUNCTION getNoOfDaysInYearDateTime(dt, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    INTEGER(c_int) :: getNoOfDaysInYearDateTime
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    getNoOfDaysInYearDateTime = my_getnoofdaysinyeardatetime(c_LOC(dt))
    IF (getNoOfDaysInYearDateTime == 0 .AND. PRESENT(errno)) errno = 4 * 100 + 16
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
  FUNCTION getDayOfYearFromDateTime(dt, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    INTEGER(c_int) :: getDayOfYearFromDateTime
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    getDayOfYearFromDateTime = my_getdayofyearfromdatetime(c_LOC(dt))
    IF (getDayOfYearFromDateTime == 0 .AND. PRESENT(errno)) errno = 4 * 100 + 17
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
  FUNCTION getNoOfSecondsElapsedInMonthDateTime(dt, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    INTEGER(c_int64_t) :: getNoOfSecondsElapsedInMonthDateTime
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    getNoOfSecondsElapsedInMonthDateTime = my_getnoofsecondselapsedinmonthdatetime(c_LOC(dt))
    IF (getNoOfSecondsElapsedInMonthDateTime == -1 .AND. PRESENT(errno)) errno = 4 * 100 + 18
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
  FUNCTION getNoOfSecondsElapsedInDayDateTime(dt, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    INTEGER(c_int) :: getNoOfSecondsElapsedInDayDateTime
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    getNoOfSecondsElapsedInDayDateTime = my_getnoofsecondselapsedindaydatetime(c_LOC(dt))
    IF (getNoOfSecondsElapsedInDayDateTime == -1 .AND. PRESENT(errno)) errno = 4 * 100 + 19
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
  SUBROUTINE getJulianDayFromDatetime(dt, jd, errno) !OK-TESTED.
    TYPE(datetime), TARGET :: dt
    TYPE(julianday), TARGET :: jd
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_getjuliandayfromdatetime(c_LOC(dt), c_LOC(jd))
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 4 * 100 + 20
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
  SUBROUTINE getDatetimeFromJulianDay(jd, dt, errno)
    TYPE(julianday), TARGET, INTENT(in) :: jd
    TYPE(datetime), TARGET, INTENT(out) :: dt
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_getdatetimefromjulianday(c_LOC(jd), c_LOC(dt))
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 4 * 100 + 21
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
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_int32_t, c_int64_t, c_float, c_double, c_char, &
       &                                 c_ptr, c_null_char, c_loc, c_f_pointer, c_associated
  USE mtime_c_bindings
  USE mtime_constants
  USE mtime_error_handling
  !
  USE mtime_datetime
  USE mtime_date
  USE mtime_juliandelta
  !use mtime_time
  !
  IMPLICIT NONE
  !
!  private
  !
  PUBLIC :: divisionquotienttimespan
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
  !
  !
  !
  INTERFACE newTimedelta
    MODULE PROCEDURE newtimedeltafromstring
    MODULE PROCEDURE newtimedeltafromraw
    MODULE PROCEDURE newtimedeltafromraw_yi8
    MODULE PROCEDURE newtimedeltafromconstructandcopy
  END INTERFACE newTimedelta
  !
  INTERFACE OPERATOR (+)
    MODULE PROCEDURE addtimedeltatodatetime
    MODULE PROCEDURE adddatetimetotimedelta
    MODULE PROCEDURE addtimedeltatodate
    MODULE PROCEDURE adddatetotimedelta
    MODULE PROCEDURE elementwiseAddTimeDeltatoTimeDelta
  END INTERFACE OPERATOR (+)
  !
  INTERFACE OPERATOR (-)
    MODULE PROCEDURE getTimeDeltaFromDate
    MODULE PROCEDURE getTimeDeltaFromDateTime
  END INTERFACE OPERATOR (-)
  !
  INTERFACE OPERATOR (*)
    MODULE PROCEDURE elementwiseScalarMultiplyTimeDelta
    MODULE PROCEDURE elementwiseScalarMultiplyTimeDeltaInv
    MODULE PROCEDURE elementwiseScalarMultiplyTimeDelta_long
    MODULE PROCEDURE elementwiseScalarMultiplyTimeDeltaInv_long
    MODULE PROCEDURE elementwiseScalarMultiplyTimeDelta_real
  END INTERFACE OPERATOR (*)
  !
  INTERFACE OPERATOR (>)
  MODULE PROCEDURE timedelta_gt
  END INTERFACE OPERATOR (>)
  !
  INTERFACE OPERATOR (<)
  MODULE PROCEDURE timedelta_lt
  END INTERFACE OPERATOR (<)
  !
  INTERFACE OPERATOR (<=)
  MODULE PROCEDURE timedelta_lt_or_eq
  END INTERFACE OPERATOR (<=)
  !
  INTERFACE OPERATOR (>=)
  MODULE PROCEDURE timedelta_gt_or_eq
  END INTERFACE OPERATOR (>=)
  !
  INTERFACE OPERATOR (==)
  MODULE PROCEDURE timedelta_eq
  END INTERFACE OPERATOR (==)
  !
  INTERFACE OPERATOR (/=)
  MODULE PROCEDURE timedelta_ne
  END INTERFACE OPERATOR (/=)
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
  FUNCTION newtimedeltafromstring(string, errno) RESULT(ret_timedelta) !OK-TESTED.
    TYPE(timedelta), POINTER :: ret_timedelta
    CHARACTER(len=*), INTENT(in) :: string
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    ret_timedelta => NULL()
    c_pointer = my_newtimedeltafromstring(TRIM(ADJUSTL(string))//c_null_char)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) THEN
      errno = 5 * 100 + 1
    ELSE
      CALL c_f_POINTER(c_pointer, ret_timedelta)
    ENDIF
  END FUNCTION newtimedeltafromstring
  !
  !! @param[out]       errno       optional, error message
  FUNCTION newtimedeltafromraw(sign, year, month, day, hour, minute, second, ms, errno) RESULT(ret_timedelta)
    TYPE(timedelta), POINTER :: ret_timedelta
    CHARACTER(len=*), INTENT(in) :: sign
    INTEGER(c_int), INTENT(in) :: year, month, day, hour, minute, second, ms
    TYPE(c_ptr) :: c_pointer
    CHARACTER(c_char) ::c_sign
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_sign = SIGN(1:1)
    c_pointer = my_newrawtimedelta(c_sign, INT(year,c_int64_t), month, day, hour, minute, second, ms)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 5 * 100 + 2
    CALL c_f_POINTER(c_pointer, ret_timedelta)
  END FUNCTION newtimedeltafromraw
  !
  !! @param[out]       errno       optional, error message
  FUNCTION newtimedeltafromraw_yi8(sign, year, month, day, hour, minute, second, ms, errno) RESULT(ret_timedelta)
    TYPE(timedelta), POINTER :: ret_timedelta
    CHARACTER(len=*), INTENT(in) :: sign
    INTEGER(c_int64_t), INTENT(in) :: year
    INTEGER(c_int), INTENT(in) :: month, day, hour, minute, second, ms
    TYPE(c_ptr) :: c_pointer
    CHARACTER(c_char) ::c_sign
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_sign = SIGN(1:1)
    c_pointer = my_newrawtimedelta(c_sign, year, month, day, hour, minute, second, ms)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 5 * 100 + 3
    CALL c_f_POINTER(c_pointer, ret_timedelta)
  END FUNCTION newtimedeltafromraw_yi8
  !
  !! @param[out]       errno       optional, error message
  FUNCTION newtimedeltafromconstructandcopy(src, errno) RESULT(dest) !OK-TESTED.
    TYPE(timedelta), POINTER :: dest
    TYPE(timedelta), TARGET :: src
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_constructandcopytimedelta(c_LOC(src))
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 5 * 100 + 4
    CALL c_f_POINTER(c_pointer, dest)
  END FUNCTION newtimedeltafromconstructandcopy
  !>
  !! @brief Destructor of TimeDelta.
  !!
  !! @param  my_timedelta
  !!         A pointer to typetimedelta. my_timedelta is deallocated.
  SUBROUTINE deallocateTimedelta(my_timedelta) !OK-TESTED.
    TYPE(timedelta), POINTER :: my_timedelta
    CALL my_deallocatetimedelta(c_LOC(my_timedelta))
    NULLIFY(my_timedelta)
  END SUBROUTINE deallocateTimedelta
  !
  FUNCTION timedelta_gt(op1, op2) RESULT(gt)
    LOGICAL :: gt
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(c_LOC(op1), c_LOC(op2))
    IF (ret == 1) THEN
      gt = .TRUE.
    ELSE
      gt = .FALSE.
    ENDIF
  END FUNCTION timedelta_gt
  !
  FUNCTION timedelta_lt(op1, op2) RESULT(lt)
    LOGICAL :: lt
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(c_LOC(op1), c_LOC(op2))
    IF (ret == -1) THEN
      lt = .TRUE.
    ELSE
      lt = .FALSE.
    ENDIF
  END FUNCTION timedelta_lt
  !
  FUNCTION timedelta_lt_or_eq(op1, op2) RESULT(lt_or_eq)
    LOGICAL :: lt_or_eq
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(c_LOC(op1), c_LOC(op2))
    IF ((ret == 0) .OR. (ret == -1)) THEN
      lt_or_eq = .TRUE.
    ELSE
      lt_or_eq = .FALSE.
    ENDIF
  END FUNCTION timedelta_lt_or_eq
  !
  FUNCTION timedelta_gt_or_eq(op1, op2) RESULT(gt_or_eq)
    LOGICAL :: gt_or_eq
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(c_LOC(op1), c_LOC(op2))
    IF ((ret == 0) .OR. (ret == 1)) THEN
      gt_or_eq = .TRUE.
    ELSE
      gt_or_eq = .FALSE.
    ENDIF
  END FUNCTION timedelta_gt_or_eq
  !
  FUNCTION timedelta_eq(op1, op2) RESULT(eq)
    LOGICAL :: eq
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(c_LOC(op1), c_LOC(op2))
    IF (ret == 0) THEN
      eq = .TRUE.
    ELSE
      eq = .FALSE.
    ENDIF
  END FUNCTION timedelta_eq
  !
  FUNCTION timedelta_ne(op1, op2) RESULT(ne)
    LOGICAL :: ne
    TYPE(timedelta), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    INTEGER(c_int) :: ret
    ret = my_comparetimedelta(c_LOC(op1), c_LOC(op2))
    IF (ret /= 0) THEN
      ne = .TRUE.
    ELSE
      ne = .FALSE.
    ENDIF
  END FUNCTION timedelta_ne
  !>
  !! @brief Get the TimeDelta between two Dates op1 and op2 as (op1-op22).
  !!
  !! Routine getTimeDeltaFromDate 'substracts' two Dates and returns the TimeDelta between
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
  !!         A pointer to TimeDelta containing the result of substraction.
  FUNCTION getTimeDeltaFromDate(op1,op2) RESULT(ret) !OK-TESTED.
    TYPE(timedelta), TARGET :: ret
    TYPE(date), TARGET, INTENT(in)  :: op1, op2
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_gettimedeltafromdate(c_LOC(op1),c_LOC(op2),c_LOC(ret))
  END FUNCTION getTimeDeltaFromDate
  !>
  !! @brief Get the TimeDelta between two DateTimes op1 and op2 as (op1-op2).
  !!
  !! Routine getTimeDeltaFromDateTime 'substracts' two DateTime's and returns the TimeDelta between
  !! them. Each datetime is converted to an equivalent Julian Date. Substraction is then performed
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
  !!        A pointer to TimeDelta containing the result of substraction.
  FUNCTION getTimeDeltaFromDateTime(op1,op2) RESULT(ret) !OK-TESTED.
    TYPE(timedelta), TARGET :: ret
    TYPE(datetime), TARGET, INTENT(in) :: op1, op2
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_gettimedeltafromdatetime(c_LOC(op1),c_LOC(op2),c_LOC(ret))
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
  FUNCTION getTotalMilliSecondsTimeDelta(td, dt, errno)  !OK-TESTED.
    INTEGER(c_int64_t) :: getTotalMilliSecondsTimeDelta
    TYPE(timedelta), TARGET, INTENT(in):: td
    TYPE(datetime), TARGET, INTENT(in) :: dt
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    getTotalMilliSecondsTimeDelta = my_gettotalmillisecondstimedelta(c_LOC(td), c_LOC(dt))
    IF ((getTotalMilliSecondsTimeDelta == 0) .AND. PRESENT(errno)) errno = 5 * 100 + 8
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
  FUNCTION getTotalSecondsTimeDelta(td, dt, errno) !OK-TESTED.
    INTEGER(c_int64_t) :: getTotalSecondsTimeDelta
    TYPE(timedelta), TARGET, INTENT(in) :: td
    TYPE(datetime), TARGET, INTENT(in) :: dt
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    getTotalSecondsTimeDelta = my_gettotalsecondstimedelta(c_LOC(td), c_LOC(dt))
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
  SUBROUTINE timedeltaToString(my_timedelta, string, errno) !OK-TESTED.
    TYPE(timedelta), TARGET :: my_timedelta
    CHARACTER(len=max_timedelta_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    string(:) = ''
    dummy_ptr = my_timedeltatostring(c_LOC(my_timedelta), string)
    IF (.NOT. c_ASSOCIATED(dummy_ptr)) THEN
      IF (PRESENT(errno)) errno = 5 * 100 + 10
      string ='<null>'
    ELSE
      char_loop: DO i = 1 , LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    ENDIF
  END SUBROUTINE timedeltaToString
  !
  FUNCTION addTimedeltaToDatetime(op1, op2) RESULT(ret) !OK-TESTED.
    TYPE(datetime), TARGET :: ret
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_addtimedeltatodatetime(c_LOC(op1), c_LOC(op2), c_LOC(ret))
  END FUNCTION addTimedeltaToDatetime
  !
  FUNCTION addDatetimeToTimedelta(op2, op1) RESULT(ret) !OK-TESTED.
    TYPE(datetime), TARGET :: ret
    TYPE(datetime), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_addtimedeltatodatetime(c_LOC(op1), c_LOC(op2), c_LOC(ret))
  END FUNCTION addDatetimeToTimedelta
  !
  FUNCTION addTimedeltaToDate(op1, op2) RESULT(ret) !OK-TESTED.
    TYPE(date), TARGET :: ret
    TYPE(date), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_addtimedeltatodate(c_LOC(op1), c_LOC(op2), c_LOC(ret))
  END FUNCTION addTimedeltaToDate
  !
  FUNCTION addDateToTimedelta(op2, op1) RESULT(ret) !OK-TESTED.
    TYPE(date), TARGET :: ret
    TYPE(date), TARGET, INTENT(in) :: op1
    TYPE(timedelta), TARGET, INTENT(in) :: op2
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_addtimedeltatodate(c_LOC(op1), c_LOC(op2), c_LOC(ret))
  END FUNCTION addDateToTimedelta
  !
  FUNCTION elementwiseScalarMultiplyTimeDelta(base_td, ilambda) RESULT(scaled_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: scaled_td
    INTEGER(c_int32_t), INTENT(in) :: ilambda
    TYPE(timedelta), TARGET, INTENT(in) :: base_td
    INTEGER(c_int64_t) :: lambda
    TYPE(c_ptr) :: dummy_ptr
    lambda = INT(ilambda, c_int64_t)
    dummy_ptr = my_elementwisescalarmultiplytimedelta(c_LOC(base_td), lambda, c_LOC(scaled_td))
  END FUNCTION elementwiseScalarMultiplyTimeDelta
  !
  FUNCTION elementwiseScalarMultiplyTimeDeltaInv(ilambda, base_td) RESULT(scaled_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: scaled_td
    INTEGER(c_int32_t), INTENT(in) :: ilambda
    TYPE(timedelta), TARGET, INTENT(in) :: base_td
    INTEGER(c_int64_t) :: lambda
    TYPE(c_ptr) :: dummy_ptr
    lambda = INT(ilambda, c_int64_t)
    dummy_ptr = my_elementwisescalarmultiplytimedelta(c_LOC(base_td), lambda, c_LOC(scaled_td))
  END FUNCTION elementwiseScalarMultiplyTimeDeltaInv
  !
  FUNCTION elementwiseScalarMultiplyTimeDelta_long(base_td, lambda) RESULT(scaled_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: scaled_td
    INTEGER(c_int64_t), INTENT(in) :: lambda
    TYPE(timedelta), TARGET, INTENT(in) :: base_td
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_elementwisescalarmultiplytimedelta(c_LOC(base_td), lambda, c_LOC(scaled_td))
  END FUNCTION elementwiseScalarMultiplyTimeDelta_long
  !
  FUNCTION elementwiseScalarMultiplyTimeDeltaInv_long(lambda, base_td) RESULT(scaled_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: scaled_td
    INTEGER(c_int64_t), INTENT(in) :: lambda
    TYPE(timedelta), TARGET, INTENT(in) :: base_td
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_elementwisescalarmultiplytimedelta(c_LOC(base_td), lambda, c_LOC(scaled_td))
  END FUNCTION elementwiseScalarMultiplyTimeDeltaInv_long
  !
  FUNCTION elementwisescalarmultiplytimedelta_real(base_td, lambda) RESULT(scaled_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: scaled_td
    REAL(c_double), INTENT(in) :: lambda
    TYPE(timedelta), TARGET, INTENT(in) :: base_td
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_elementwisescalarmultiplytimedeltadp(c_LOC(base_td), lambda, c_LOC(scaled_td))
  END FUNCTION elementwisescalarmultiplytimedelta_real
  !
  FUNCTION elementwiseAddTimeDeltatoTimeDelta(td1, td2) RESULT(added_td) !OK-TESTED.
    TYPE(timedelta), TARGET :: added_td
    TYPE(timedelta), TARGET, INTENT(in) :: td1, td2
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_elementwiseaddtimedeltatotimedelta(c_LOC(td1), c_LOC(td2), c_LOC(added_td))
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
  FUNCTION moduloTimedelta(a, p, quot) RESULT(rem)
    TYPE(timedelta), TARGET, INTENT(in) :: a
    TYPE(timedelta), TARGET, INTENT(in) :: p
    INTEGER(c_int64_t), TARGET, INTENT(out) :: quot
    INTEGER(c_int64_t) :: rem
    rem = my_modulotimedelta(c_LOC(a), c_LOC(p), c_LOC(quot))
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
  SUBROUTINE getPTStringFromMS(ms, string, errno) !OK-TESTED.
    INTEGER(c_int64_t), INTENT(in) :: ms
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    string(:) = ''
    dummy_ptr = my_getptstringfromms(ms, string)
    IF (.NOT. c_ASSOCIATED(dummy_ptr)) THEN
      IF (PRESENT(errno)) errno = 5 * 100 + 11
      string = '<null>'
    ELSE
      char_loop: DO i = 1 , LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    ENDIF
  END SUBROUTINE getPTStringFromMS
  !
  SUBROUTINE getPTStringFromSecondsInt(s, string, errno) !OK-TESTED.
    INTEGER(c_int64_t) :: s
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_getptstringfromsecondsint(s, string)
    IF (.NOT. c_ASSOCIATED(dummy_ptr)) THEN
      IF (PRESENT(errno)) errno = 5 * 100 + 11
      string = '<null>'
    ELSE
      char_loop: DO i = 1 , LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    ENDIF
  END SUBROUTINE getPTStringFromSecondsInt
  !
  SUBROUTINE getPTStringFromSecondsFloat(s, string, errno) !OK-TESTED.
    REAL(c_float) :: s
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_getptstringfromsecondsfloat(s, string)
    IF (.NOT. c_ASSOCIATED(dummy_ptr)) THEN
      IF (PRESENT(errno)) errno = 5 * 100 + 11
      string = '<null>'
    ELSE
      char_loop: DO i = 1 , LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    ENDIF
  END SUBROUTINE getPTStringFromSecondsFloat
  !
  SUBROUTINE getPTStringFromSecondsDouble(s, string, errno) !OK-TESTED.
    REAL(c_double) :: s
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_getptstringfromsecondsdouble(s, string)
    IF (.NOT. c_ASSOCIATED(dummy_ptr)) THEN
      IF (PRESENT(errno)) errno = 5 * 100 + 11
      string = '<null>'
    ELSE
      char_loop: DO i = 1 , LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    ENDIF
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
  SUBROUTINE getPTStringFromMinutes(m, string, errno) !OK-TESTED.
    INTEGER(c_int64_t) :: m
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_getptstringfromminutes(m, string)
    IF (.NOT. c_ASSOCIATED(dummy_ptr)) THEN
      IF (PRESENT(errno)) errno = 5 * 100 + 11
      string = '<null>'
    ELSE
      char_loop: DO i = 1 , LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    ENDIF
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
  SUBROUTINE getPTStringFromHours(h, string, errno) !OK-TESTED.
    INTEGER(c_int64_t)                   :: h
    CHARACTER(len=*) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    string(:) = ''
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_getptstringfromhours(h, string)
    IF (.NOT. c_ASSOCIATED(dummy_ptr)) THEN
      IF (PRESENT(errno)) errno = 5 * 100 + 11
      string = '<null>'
    ELSE
      char_loop: DO i = 1 , LEN(string)
        IF (string(i:i) == c_null_char) EXIT char_loop
      END DO char_loop
      string(i:LEN(string)) = ' '
    ENDIF
  END SUBROUTINE getPTStringFromHours

  !>
  !! @brief Convert time delta to "Julian calendar delta".
  !!
  SUBROUTINE timeDeltaToJulianDelta(td,dt,jd)
    TYPE(c_ptr) :: dummy_ptr
    TYPE(timedelta), TARGET, INTENT(in) :: td
    TYPE(datetime),  TARGET, INTENT(in) :: dt
    TYPE(juliandelta),  TARGET, INTENT(out) :: jd
    dummy_ptr = my_timedeltatojuliandelta(c_LOC(td),c_LOC(dt),c_LOC(jd))
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
  SUBROUTINE divideTimeDeltaInSeconds(dividend, divisor, quotient, errna)!OK-UNTESTED.
    TYPE(timedelta), TARGET, INTENT(in) :: dividend
    TYPE(timedelta), TARGET, INTENT(in) :: divisor
    TYPE(divisionquotienttimespan), TARGET, INTENT(out) :: quotient
    INTEGER, INTENT(out), OPTIONAL :: errna
    TYPE(c_ptr) :: dummy_ptr
    IF (PRESENT(errna)) errna = 0 ! FIXME: no_error
    dummy_ptr = my_dividetimedeltainseconds(c_LOC(dividend), c_LOC(divisor), c_LOC(quotient))
    IF (PRESENT(errna) .AND. .NOT. c_ASSOCIATED(dummy_ptr)) THEN
      errna = errna + 2  ! increment error number by 2, see below for an explanation.
    ENDIF
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
  SUBROUTINE divideTwoDatetimeDiffsInSeconds(dt1_dividend,dt2_dividend, &
      &                                      dt1_divisor, dt2_divisor,  &
      &                                      denominator, quotient)
    TYPE(datetime), TARGET, INTENT(in) :: dt1_dividend
    TYPE(datetime), TARGET, INTENT(in) :: dt2_dividend
    TYPE(datetime), TARGET, INTENT(in) :: dt1_divisor
    TYPE(datetime), TARGET, INTENT(in) :: dt2_divisor
    INTEGER(c_int64_t), TARGET, INTENT(out) :: denominator
    TYPE(divisionquotienttimespan), TARGET, INTENT(out) :: quotient
    TYPE(c_ptr) :: dummy_ptr
    dummy_ptr = my_dividetwodatetimediffsinseconds(c_LOC(dt1_dividend), c_LOC(dt2_dividend),  &
        &                                                c_LOC(dt1_divisor), c_LOC(dt2_divisor),  &
        &                                                c_LOC(denominator), c_LOC(quotient))
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
  SUBROUTINE divideDatetimeDifferenceInSeconds(dt1, dt2, divisor, quotient, errna)
    TYPE(datetime), TARGET, INTENT(in) :: dt1
    TYPE(datetime), TARGET, INTENT(in) :: dt2
    TYPE(timedelta), TARGET, INTENT(in) :: divisor
    TYPE(divisionquotienttimespan), TARGET, INTENT(out) :: quotient
    INTEGER, INTENT(out), OPTIONAL :: errna
    TYPE(c_ptr) :: dummy_ptr
    IF (PRESENT(errna)) errna = 0 ! FIXME: no_error
    dummy_ptr = my_dividedatetimedifferenceinseconds(c_LOC(dt1), c_LOC(dt2), c_LOC(divisor), c_LOC(quotient))
    IF (PRESENT(errna) .AND. .NOT. c_ASSOCIATED(dummy_ptr)) THEN
      errna = errna + 2  ! increment error number by 2, see below for an explanation.
    ENDIF
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
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t, c_char, c_null_char, c_bool, c_ptr, &
       &                                 c_null_ptr, c_loc, c_f_pointer, c_associated
  USE mtime_c_bindings
  USE mtime_constants
  USE mtime_error_handling
  !
  USE mtime_datetime
  USE mtime_timedelta
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
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
  !
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
  FUNCTION newEventWithString(name, referenceDate, firstdate, lastDate, interval, offset, errno) RESULT(ret_event) !OK-TESTED.
    TYPE(event), POINTER :: ret_event
    CHARACTER(len=*), INTENT(in)           :: name
    CHARACTER(len=*), INTENT(in)           :: referenceDate
    CHARACTER(len=*), INTENT(in)           :: firstDate
    CHARACTER(len=*), INTENT(in)           :: lastDate
    CHARACTER(len=*), INTENT(in)           :: interval
    CHARACTER(len=*), INTENT(in), OPTIONAL :: offset
    CHARACTER(len=32)                      :: zeroOffset = "PT00S" !Optional offset's proxy string.
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    IF(PRESENT(offset)) THEN !offset is not null.
         c_pointer = my_newevent(TRIM(ADJUSTL(name))//c_null_char,          &
         &                  TRIM(ADJUSTL(referenceDate))//c_null_char,      &
         &                  TRIM(ADJUSTL(firstDate))//c_null_char,          &
         &                  TRIM(ADJUSTL(lastDate))//c_null_char,           &
         &                  TRIM(ADJUSTL(interval))//c_null_char,           &
         &                  TRIM(ADJUSTL(offset))//c_null_char)
    ELSE                     !offset is null. Pass a 0 offset.
         c_pointer = my_newevent(TRIM(ADJUSTL(name))//c_null_char,    &
         &                  TRIM(ADJUSTL(referenceDate))//c_null_char,      &
         &                  TRIM(ADJUSTL(firstDate))//c_null_char,          &
         &                  TRIM(ADJUSTL(lastDate))//c_null_char,           &
         &                  TRIM(ADJUSTL(interval))//c_null_char,           &
         &                  TRIM(ADJUSTL(zeroOffset))//c_null_char)
    END IF
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 6 * 100 + 1
    CALL c_f_POINTER(c_pointer, ret_event)
  END FUNCTION newEventWithString
  !
  !! @param[out]       errno       optional, error message
  FUNCTION newEventWithDataTypes(name, referenceDate, firstdate, lastDate, interval, offset, errno) RESULT(ret_event) !OK-TESTED.
    TYPE(event), POINTER :: ret_event
    CHARACTER(len=*), INTENT(in) :: name
    TYPE(datetime), POINTER :: referenceDate
    TYPE(datetime), POINTER :: firstDate
    TYPE(datetime), POINTER :: lastDate
    TYPE(timedelta), POINTER :: interval
    TYPE(timedelta), POINTER, OPTIONAL :: offset
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    IF (PRESENT(offset)) THEN
    c_pointer = my_neweventwithdatatypes(TRIM(ADJUSTL(name))//c_null_char, &
         &                  c_LOC(referenceDate), c_LOC(firstDate), c_LOC(lastDate), c_LOC(interval), c_LOC(offset))
    ELSE
    c_pointer = my_neweventwithdatatypes(TRIM(ADJUSTL(name))//c_null_char, &
         &                  c_LOC(referenceDate), c_LOC(firstDate), c_LOC(lastDate), c_LOC(interval), c_null_ptr)
    END IF
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 6 * 100 + 1
    CALL c_f_POINTER(c_pointer, ret_event)
  END FUNCTION newEventWithDataTypes
  !>
  !! @brief Destructor of Event.
  !!
  !! @param  my_event
  !!         A pointer to type event. my_event is deallocated.
  !!
  !! WARNING: If my_event was added to a group, this should never be called;
  !! use removeEventFromEventGroup instead.
  SUBROUTINE deallocateEvent(my_event)
    TYPE(event), POINTER :: my_event
    CALL my_deallocateevent(c_LOC(my_event))
    my_event => NULL()
  END SUBROUTINE deallocateEvent
  !>
  FUNCTION constructAndCopyEvent(my_event, errno) RESULT(ret_event) !OK-TESTED.
    TYPE(event), POINTER :: ret_event
    TYPE(event), TARGET  :: my_event
    INTEGER, OPTIONAL:: errno
    TYPE(c_ptr) :: c_pointer
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_constructandcopyevent(c_LOC(my_event))
    IF ((.NOT. C_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 6 * 100 + 1
    CALL c_f_POINTER(c_pointer, ret_event)
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
  SUBROUTINE eventToString(my_event, string, errno) !TODO:C code still incomplete.
    TYPE(event), POINTER :: my_event
    CHARACTER(len=max_eventname_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_eventtostring(c_LOC(my_event), string)
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 6 * 100 + 3
    char_loop: DO i = 1 , LEN(string)
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
  !!         A pointer to type datetime. This is the 'current' datetime of the system.
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
  FUNCTION isCurrentEventActive(my_event, my_datetime, plus_slack, minus_slack) RESULT(ret)
    TYPE(event), POINTER :: my_event
    TYPE(datetime), TARGET :: my_datetime
    TYPE(timedelta), POINTER, OPTIONAL :: plus_slack
    TYPE(timedelta), POINTER, OPTIONAL :: minus_slack
    LOGICAL(c_bool) :: ret
    IF (PRESENT(plus_slack) .AND. PRESENT(minus_slack)) THEN
      ret = my_isCurrentEventActive(c_LOC(my_event), c_LOC(my_datetime), &
            &                   c_LOC(plus_slack),  c_LOC(minus_slack))

    ELSE IF (PRESENT(plus_slack) .AND. .NOT.PRESENT(minus_slack)) THEN
      ret = my_isCurrentEventActive(c_LOC(my_event), c_LOC(my_datetime), &
            &                           c_LOC(plus_slack),  c_null_ptr)

    ELSE IF (.NOT.PRESENT(plus_slack) .AND. PRESENT(minus_slack)) THEN
      ret = my_isCurrentEventActive(c_LOC(my_event), c_LOC(my_datetime), &
            &                          c_null_ptr,  c_LOC(minus_slack))

    ELSE
      ret = my_isCurrentEventActive(c_LOC(my_event), c_LOC(my_datetime), &
            &                                   c_null_ptr, c_null_ptr)

    END IF
  END FUNCTION isCurrentEventActive
  !>
  !! @brief Checks, if next event is on a new day
  !!
  !! @param[in] my_event
  !!        A pointer to a type event
  !!
  !! @returns ret
  !!        Logical: true if next event is on new day
  FUNCTION isEventNextInNextDay(my_event) RESULT(ret)
    TYPE(event), POINTER :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_iseventnextinnextday(c_LOC(my_event))
  END FUNCTION isEventNextInNextDay
  !>
  !! @brief Checks, if next event is in a new month
  !!
  !! @param[in] my_event
  !!        A pointer to a type event
  !!
  !! @returns ret
  !!        Logical: true if next event is in a new month
  FUNCTION iseventNextInNextMonth(my_event) RESULT(ret)
    TYPE(event), POINTER :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_iseventnextinnextmonth(c_LOC(my_event))
  END FUNCTION iseventNextInNextMonth
  !>
  !! @brief Checks, if next event is in a new year
  !!
  !! @param[in] my_event
  !!        A pointer to a type event
  !!
  !! @returns ret
  !!        Logical: true if next event is in a new year
  FUNCTION iseventNextInNextYear(my_event) RESULT(ret)
    TYPE(event), POINTER :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_iseventnextinnextyear(c_LOC(my_event))
  END FUNCTION iseventNextInNextYear
  !>
  !! @brief Get the Datetime when 'this' event will be triggered next.
  !!
  !! WARNING: The value returned is with-respect-to current_dt and not
  !!          a true copy of triggerNextEventDateTime in the event
  !!          data structure.
  !!
  !! @param[in] my_event
  !!         A pointer to type event. This is the event being queried.
  !!
  !! @param[in]  my_currentdatetime
  !!         A pointer to type datetime. The next trigger datetime is copied here.
  !!
  !! @param[out] my_datetime
  !!         A pointer to type datetime with next-trigger datetime.
  !!
  !! @param[out]       errno       optional, error message
  SUBROUTINE getTriggerNextEventAtDateTime(my_event,my_currentdatetime,my_datetime, errno)
    TYPE(event), TARGET, INTENT(in) :: my_event
    TYPE(datetime), TARGET, INTENT(in) :: my_currentdatetime
    TYPE(datetime), TARGET :: my_datetime
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_gettriggernexteventatdatetime(c_LOC(my_event),c_LOC(my_currentdatetime),c_LOC(my_datetime))
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 6 * 100 + 8
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
  !!         A pointer to type datetime with last-trigger datetime.
  !!
  !! @param[out]       errno       optional, error message
  SUBROUTINE getTriggeredPreviousEventAtDateTime(my_event,my_datetime,errno)
    TYPE(event), TARGET, INTENT(in) :: my_event
    TYPE(datetime), TARGET :: my_datetime
    TYPE(c_ptr) :: dummy_ptr
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_gettriggeredpreviouseventatdatetime(c_LOC(my_event),c_LOC(my_datetime))
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 6 * 100 + 9
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
  FUNCTION getEventId(my_event) RESULT(ret_evtid) !OK-TESTED.
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
  SUBROUTINE getEventName(my_event, string, errno) !OK-TESTED.
    TYPE(event), POINTER :: my_event
    CHARACTER(len=max_eventname_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL :: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_geteventname(c_LOC(my_event), string)
    IF ((.NOT.c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 6 * 100 + 11
    char_loop: DO i = 1 , LEN(string)
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
  FUNCTION getEventReferenceDateTime(my_event) RESULT(ret_referenceDateTime) !OK-TESTED.
    TYPE(datetime), POINTER :: ret_referenceDateTime
    TYPE(event), POINTER :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%eventReferenceDatetime
    IF (c_ASSOCIATED(c_pointer)) THEN
      CALL c_f_POINTER(c_pointer, ret_referenceDateTime)
    ELSE
      ret_referenceDateTime => NULL()
    ENDIF
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
  FUNCTION getEventFirstDateTime(my_event) RESULT(ret_eventFirstDateTime) !OK-TESTED.
    TYPE(datetime), POINTER :: ret_eventFirstDateTime
    TYPE(event), POINTER :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%eventFirstDateTime
    IF (c_ASSOCIATED(c_pointer)) THEN
      CALL c_f_POINTER(c_pointer, ret_eventFirstDateTime)
    ELSE
      ret_eventFirstDateTime => NULL()
    ENDIF
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
  FUNCTION getEventLastDateTime(my_event) RESULT(ret_eventLastDateTime) !OK-TESTED.
    TYPE(datetime), POINTER :: ret_eventLastDateTime
    TYPE(event), POINTER :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%eventLastDateTime
    IF (c_ASSOCIATED(c_pointer)) THEN
      CALL c_f_POINTER(c_pointer, ret_eventLastDateTime)
    ELSE
      ret_eventLastDateTime => NULL()
    ENDIF
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
  FUNCTION getEventInterval(my_event) RESULT(ret_eventInterval) !OK-TESTED.
    TYPE(timedelta), POINTER :: ret_eventInterval
    TYPE(event), POINTER :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%eventInterval
    IF (c_ASSOCIATED(c_pointer)) THEN
      CALL c_f_POINTER(c_pointer, ret_eventInterval)
    ELSE
      ret_eventInterval => NULL()
    ENDIF
  END FUNCTION getEventInterval
  !>
  !! @brief Check if event is first
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is first
  !!
  FUNCTION getNextEventIsFirst(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_getnexteventisfirst(c_LOC(my_event))
  END FUNCTION getNextEventIsFirst
  !>
  !! @brief Check if event is first in day
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is first in day
  !!
  FUNCTION getEventisFirstInDay(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventisfirstinday(c_LOC(my_event))
  END FUNCTION getEventisFirstInDay
  !>
  !! @brief Check if event is first in month
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is first in month
  !!
  FUNCTION getEventisFirstInMonth(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventisfirstinmonth(c_LOC(my_event))
  END FUNCTION getEventisFirstInMonth
  !>
  !! @brief Check if event is first in year
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is first in year
  !!
  FUNCTION getEventisFirstInYear(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventisfirstinyear(c_LOC(my_event))
  END FUNCTION getEventisFirstInYear
  !>
  !! @brief Check if event is last in day
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is last in day
  !!

  FUNCTION getEventisLastInDay(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventislastinday(c_LOC(my_event))
  END FUNCTION getEventisLastInDay
  !>
  !! @brief Check if event is last in month
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is last in month
  !!
  FUNCTION getEventisLastInMonth(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventislastinmonth(c_LOC(my_event))
  END FUNCTION getEventisLastInMonth
  !>
  !! @brief Check if event is last in year
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret
  !!        Logical: true if event is last in year
  !!
  FUNCTION getEventisLastInYear(my_event) RESULT(ret)
    TYPE(event), TARGET, INTENT(in) :: my_event
    LOGICAL(c_bool) :: ret
    ret = my_geteventislastinyear(c_LOC(my_event))
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
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t, c_ptr, c_char, c_null_char, c_bool, &
       &                                 c_loc, c_f_pointer, c_associated
  USE mtime_c_bindings
  USE mtime_error_handling
  USE mtime_constants
  !
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
  !
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
  FUNCTION newEventGroup(name, errno) RESULT(ret_eventgroup) !OK-TESTED.
    TYPE(eventgroup), POINTER :: ret_eventgroup
    CHARACTER(len=*), INTENT(in) :: name
    TYPE(c_ptr) :: c_pointer
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    c_pointer = my_neweventgroup(TRIM(ADJUSTL(name))//c_null_char)
    IF ((.NOT. c_ASSOCIATED(c_pointer)) .AND. PRESENT(errno)) errno = 7 * 100 + 1
    CALL c_f_POINTER(c_pointer, ret_eventgroup)
  END FUNCTION newEventGroup
  !>
  !! @brief Destructor of EventGroup.
  !!
  !! @param[in] my_eventgroup
  !!         A pointer to type eventGroup. my_eventgroup is deallocated.
  !!
  SUBROUTINE deallocateEventGroup(my_eventgroup) !OK-TESTED.
    TYPE(eventgroup), POINTER :: my_eventgroup
    CALL my_deallocateeventgroup(c_LOC(my_eventgroup))
    my_eventgroup => NULL()
  END SUBROUTINE deallocateEventGroup
  !>
  !! @brief Add new event to an eventgroup.
  !!
  !! @param  my_event
  !!         A pointer to type event. The event to be added.
  !!
  !! @param  my_eventgroup
  !!         A pointer to type eventgroup. The eventgroup where the event is added.
  !!
  !! @return ret
  !!         true/false indicating success or failure of addition.
  FUNCTION addEventToEventGroup(my_event, my_eventgroup) RESULT(ret) !OK-TESTED.
    LOGICAL :: ret
    TYPE(event), POINTER :: my_event
    TYPE(eventgroup), POINTER :: my_eventgroup
    ret = my_addeventtoeventgroup(c_LOC(my_event), c_LOC(my_eventgroup))
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
  FUNCTION removeEventfromEventGroup(my_name, my_eventgroup) RESULT(ret) !OK-TESTED.
    LOGICAL :: ret
    CHARACTER(len=*), INTENT(in) :: my_name
    TYPE(eventgroup), POINTER :: my_eventgroup
    ret = my_removeeventfromeventgroup(TRIM(ADJUSTL(my_name))//c_null_char, c_LOC(my_eventgroup))
  END FUNCTION removeEventFromEventGroup
  !>
  !! @brief Get event group id
  !!
  !! @param[in]  my_eventgroup
  !!         A pointer to  type eventgroup.
  !!
  !! @return ret_grpid
  !!         The event group id
  FUNCTION getEventGroupId(my_eventgroup) RESULT(ret_grpid) !OK-TESTED.
    INTEGER(c_int64_t) :: ret_grpid
    TYPE(eventgroup), POINTER :: my_eventgroup
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
  SUBROUTINE getEventGroupName(my_eventgroup, string, errno)  !TESTED-OK.
    TYPE(eventgroup), POINTER :: my_eventgroup
    CHARACTER(len=max_groupname_str_len) :: string
    TYPE(c_ptr) :: dummy_ptr
    INTEGER :: i
    INTEGER, OPTIONAL:: errno
    IF (PRESENT(errno)) errno = 0
    dummy_ptr = my_geteventgroupname(c_LOC(my_eventgroup), string)
    IF ((.NOT. c_ASSOCIATED(dummy_ptr)) .AND. PRESENT(errno)) errno = 7 * 100 + 6
    char_loop: DO i = 1 , LEN(string)
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
  FUNCTION getFirstEventFromEventGroup(my_eventgroup) RESULT(ret_event) !OK-TESTED.
    TYPE(event), POINTER :: ret_event
    TYPE(eventgroup), POINTER :: my_eventgroup
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_eventgroup%firstEventInGroup
    IF (c_ASSOCIATED(c_pointer)) THEN
      CALL c_f_POINTER(c_pointer, ret_event)
    ELSE
      ret_event => NULL()
    ENDIF
  END FUNCTION getFirstEventFromEventGroup
  !>
  !! @brief get the next event in an event group an event belongs to
  !!
  !! @param[in] my_event
  !!        A pointer of type event.
  !!
  !! @returns ret_event
  !!        A pointer of type event. The next event in eventgroup
  FUNCTION getNextEventFromEventGroup(my_event) RESULT(ret_event) !OK-TESTED.
    TYPE(event), POINTER :: ret_event
    TYPE(event), POINTER :: my_event
    TYPE(c_ptr) :: c_pointer
    c_pointer = my_event%nextEventInGroup
    IF (c_ASSOCIATED(c_pointer)) THEN
      CALL c_f_POINTER(c_pointer, ret_event)
    ELSE
      ret_event => NULL()
    ENDIF
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
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_char, c_null_char
  USE mtime_c_bindings
  USE mtime_error_handling
  USE mtime_constants
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
  FUNCTION getRepetitions(repetitionString) RESULT(r)
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
  SUBROUTINE splitRepetitionString(recurringTimeInterval, repetitor, start, END, duration, lrepetitor, lstart, lend, lduration)
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

    lrepetitor = .FALSE.
    lstart = .FALSE.
    lend = .FALSE.
    lduration = .FALSE.

    CALL my_splitRepetitionString(TRIM(ADJUSTL(recurringTimeInterval))//c_null_char, repetitor, start, END, duration)

    IF (repetitor(1:1) /= c_null_char) lrepetitor = .TRUE.
    char_loop1: DO i = 1 , LEN(repetitor)
      IF (repetitor(i:i) == c_null_char) EXIT char_loop1
    END DO char_loop1
    repetitor(i:LEN(repetitor)) = ' '

    IF (start(1:1) /= c_null_char) lstart = .TRUE.
    char_loop2: DO i = 1 , LEN(start)
      IF (start(i:i) == c_null_char) EXIT char_loop2
    END DO char_loop2
    start(i:LEN(start)) = ' '

    IF (END(1:1) /= c_null_char) lend = .TRUE.
    char_loop3: DO i = 1 , LEN(END)
      IF (END(i:i) == c_null_char) EXIT char_loop3
    END DO char_loop3
    END(i:LEN(END)) = ' '

    IF (duration(1:1) /= c_null_char) lduration = .TRUE.
    char_loop4: DO i = 1 , LEN(duration)
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
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_int64_t, c_bool, c_ptr, c_char
  !
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
  USE mtime_print_by_callback
  USE mtime_error_handling
  USE mtime_c_bindings
  !
  IMPLICIT NONE
  !
  PUBLIC

  !
  !
  !> @cond DOXYGEN_IGNORE_THIS
  INTEGER(c_int), BIND(c,name='NO_OF_SEC_IN_A_DAY') :: no_of_sec_in_a_day
  INTEGER(c_int), BIND(c,name='NO_OF_SEC_IN_A_HOUR') :: no_of_sec_in_a_hour
  INTEGER(c_int), BIND(c,name='NO_OF_SEC_IN_A_MINUTE') :: no_of_sec_in_a_minute
  !
  INTEGER(c_int), BIND(c,name='NO_OF_MS_IN_A_DAY') :: no_of_ms_in_a_day
  INTEGER(c_int), BIND(c,name='NO_OF_MS_IN_HALF_DAY') :: no_of_ms_in_half_day
  INTEGER(c_int), BIND(c,name='NO_OF_MS_IN_A_HOUR') :: no_of_ms_in_a_hour
  INTEGER(c_int), BIND(c,name='NO_OF_MS_IN_A_MINUTE') :: no_of_ms_in_a_minute
  INTEGER(c_int), BIND(c,name='NO_OF_MS_IN_A_SECOND') :: no_of_ms_in_a_second
  !> @endcond DOXYGEN_IGNORE_THIS


END MODULE mtime
