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
    TYPE(c_ptr) :: c_pointer

    c_pointer = my_newtimedeltafromstring(TRIM(ADJUSTL(td_string))//c_null_char)
    IF (.NOT. C_ASSOCIATED(c_pointer)) THEN
      CALL handle_errno(5*100 + 1, &
        & __FILE__, &
        & __LINE__)
      t_timedelta_assign_string%td%sign = '?'
    ELSE
      CALL C_F_POINTER(c_pointer, td_tmp)
      t_timedelta_assign_string%td = td_tmp
      t_timedelta_assign_string%td%sign = td_tmp%sign
      CALL my_deallocatetimedelta(c_pointer)
    END IF
  END FUNCTION t_timedelta_assign_string

  ! constructor for integer milliseconds (integer)
  !
  TYPE(t_timedelta) FUNCTION t_timedelta_assign_ms(td_ms)
    INTEGER, INTENT(in)                  :: td_ms
    TYPE(timedelta), POINTER             :: td_tmp
    CHARACTER(len=max_timedelta_str_len) :: td_string
    INTEGER                              :: errno
    CALL getPTStringFromMS(INT(td_ms, c_int64_t), td_string, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    td_tmp => newTimedelta(td_string, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    t_timedelta_assign_ms%td = td_tmp
    t_timedelta_assign_ms%td%sign = td_tmp%sign
    CALL deallocateTimedelta(td_tmp)
  END FUNCTION t_timedelta_assign_ms

  ! constructor for integer milliseconds (integer)
  !
  TYPE(t_timedelta) FUNCTION t_timedelta_assign_ms_i8(td_ms)
    INTEGER(c_int64_t), INTENT(in)              :: td_ms
    TYPE(timedelta), POINTER             :: td_tmp
    CHARACTER(len=max_timedelta_str_len) :: td_string
    INTEGER                              :: errno
    CALL getPTStringFromMS(INT(td_ms, c_int64_t), td_string, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    td_tmp => newTimedelta(td_string, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    t_timedelta_assign_ms_i8%td = td_tmp
    t_timedelta_assign_ms_i8%td%sign = td_tmp%sign
    CALL deallocateTimedelta(td_tmp)
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
    CLASS(t_timedelta), INTENT(in) :: this
    CLASS(t_timedelta), INTENT(in) :: td
    t_timedelta_equal = (this%td == td%td)
  END FUNCTION t_timedelta_equal

  LOGICAL FUNCTION t_timedelta_not_equal(this, td)
    CLASS(t_timedelta), INTENT(in) :: this
    CLASS(t_timedelta), INTENT(in) :: td
    t_timedelta_not_equal = (this%td /= td%td)
  END FUNCTION t_timedelta_not_equal

  LOGICAL FUNCTION t_timedelta_less_than(this, td)
    CLASS(t_timedelta), INTENT(in) :: this
    CLASS(t_timedelta), INTENT(in) :: td
    t_timedelta_less_than = (this%td < td%td)
  END FUNCTION t_timedelta_less_than

  LOGICAL FUNCTION t_timedelta_greater_than(this, td)
    CLASS(t_timedelta), INTENT(in) :: this
    CLASS(t_timedelta), INTENT(in) :: td
    t_timedelta_greater_than = (this%td > td%td)
  END FUNCTION t_timedelta_greater_than

  LOGICAL FUNCTION t_timedelta_less_than_or_equal(this, td)
    CLASS(t_timedelta), INTENT(in) :: this
    CLASS(t_timedelta), INTENT(in) :: td
    t_timedelta_less_than_or_equal = (this%td <= td%td)
  END FUNCTION t_timedelta_less_than_or_equal

  LOGICAL FUNCTION t_timedelta_greater_than_or_equal(this, td)
    CLASS(t_timedelta), INTENT(in) :: this
    CLASS(t_timedelta), INTENT(in) :: td
    t_timedelta_greater_than_or_equal = (this%td >= td%td)
  END FUNCTION t_timedelta_greater_than_or_equal

  FUNCTION t_timedelta_scalar_multiply_long(this, lambda) RESULT(scaled_td)
    TYPE(t_timedelta), TARGET :: scaled_td
    INTEGER(c_int64_t), INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this
    TYPE(timedelta), POINTER             :: td_tmp, td_tmp2
    INTEGER :: errno
    td_tmp => newTimedelta(this%td, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    NULLIFY (td_tmp2)
    td_tmp2 = td_tmp*lambda
    IF (ASSOCIATED(td_tmp2)) THEN
      CALL handle_errno(general_arithmetic_error, &
        & __FILE__, &
        & __LINE__)
      RETURN
    END IF
    scaled_td%td = td_tmp2
    scaled_td%td%sign = td_tmp2%sign
    IF (ASSOCIATED(td_tmp)) CALL deallocateTimedelta(td_tmp)
    IF (ASSOCIATED(td_tmp2)) CALL deallocateTimedelta(td_tmp2)
  END FUNCTION t_timedelta_scalar_multiply_long

  FUNCTION t_timedelta_scalar_multiply_inv_long(lambda, this) RESULT(scaled_td)
    TYPE(t_timedelta), TARGET :: scaled_td
    INTEGER(c_int64_t), INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this
    TYPE(timedelta), POINTER             :: td_tmp, td_tmp2
    INTEGER :: errno
    td_tmp => newTimedelta(this%td, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    NULLIFY (td_tmp2)
    td_tmp2 = td_tmp*lambda
    IF (ASSOCIATED(td_tmp2)) THEN
      CALL handle_errno(general_arithmetic_error, &
        & __FILE__, &
        & __LINE__)
      RETURN
    END IF
    scaled_td%td = td_tmp2
    scaled_td%td%sign = td_tmp2%sign
    IF (ASSOCIATED(td_tmp)) CALL deallocateTimedelta(td_tmp)
    IF (ASSOCIATED(td_tmp2)) CALL deallocateTimedelta(td_tmp2)
  END FUNCTION t_timedelta_scalar_multiply_inv_long

  FUNCTION t_timedelta_scalar_multiply_int(this, lambda) RESULT(scaled_td) !OK-tested
    TYPE(t_timedelta), TARGET :: scaled_td
    INTEGER(c_int32_t), INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this
    TYPE(c_ptr) :: dummy_ptr

    dummy_ptr = my_elementwisescalarmultiplytimedelta(C_LOC(this%td), INT(lambda, c_int64_t), C_LOC(scaled_td%td))
    IF (.NOT. C_ASSOCIATED(dummy_ptr)) THEN
      CALL handle_errno(0*100 + 2, &
        & __FILE__, &
        & __LINE__)
      scaled_td%td%sign = '!'
    END IF

  END FUNCTION t_timedelta_scalar_multiply_int

  FUNCTION t_timedelta_scalar_multiply_inv_int(lambda, this) RESULT(scaled_td)
    TYPE(t_timedelta), TARGET :: scaled_td
    INTEGER(c_int32_t), INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this

    scaled_td = t_timedelta_scalar_multiply_int(this, lambda)
  END FUNCTION t_timedelta_scalar_multiply_inv_int

  FUNCTION t_timedelta_scalar_multiply_real(this, lambda) RESULT(scaled_td) !OK-tested
    TYPE(t_timedelta), TARGET             :: scaled_td
    REAL(c_double), INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this
    TYPE(c_ptr)                            :: dummy_ptr

    dummy_ptr = my_elementwisescalarmultiplytimedeltadp(C_LOC(this%td), lambda, C_LOC(scaled_td%td))
    IF (.NOT. C_ASSOCIATED(dummy_ptr)) THEN
      CALL handle_errno(0*100 + 2, &
        & __FILE__, &
        & __LINE__)
      scaled_td%td%sign = '!'
    END IF

  END FUNCTION t_timedelta_scalar_multiply_real

  FUNCTION t_timedelta_scalar_multiply_inv_real(lambda, this) RESULT(scaled_td)
    TYPE(t_timedelta), TARGET :: scaled_td
    REAL(c_double), INTENT(in) :: lambda
    CLASS(t_timedelta), TARGET, INTENT(in) :: this

    scaled_td = t_timedelta_scalar_multiply_real(this, lambda)
  END FUNCTION t_timedelta_scalar_multiply_inv_real

  ! Convert t_timedelta object to string.
  !
  FUNCTION t_timedelta_toString(this) RESULT(string)
    CHARACTER(len=max_timedelta_str_len)  :: string
    CLASS(t_timedelta)                   :: this
    TYPE(c_ptr) :: c_pointer, dummy_ptr
    INTEGER :: i

    string = ""

    c_pointer = this%get_c_pointer()
    CALL handle_errno(.NOT. C_ASSOCIATED(c_pointer), 0*100 + 2, &
      & __FILE__, &
      & __LINE__)
    dummy_ptr = my_timedeltatostring(c_pointer, string)
    CALL handle_errno(.NOT. C_ASSOCIATED(dummy_ptr), 4*100 + 6, &
      & __FILE__, &
      & __LINE__)

    char_loop: DO i = 1, LEN(string)
      IF (string(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    string(i:LEN(string)) = ' '

    CALL my_deallocatetimedelta(c_pointer) !deallocate timedelta
  END FUNCTION t_timedelta_toString

  FUNCTION t_timedelta_divideInSecondsBy(this, divisor, referenceDateTime) RESULT(quotient)
    CLASS(t_timedelta), INTENT(in)           :: this
    TYPE(t_timedelta), INTENT(in)           :: divisor
    TYPE(t_datetime), INTENT(IN), OPTIONAL :: referenceDateTime
    TYPE(divisionquotienttimespan)           :: quotient
    TYPE(timedelta), POINTER  :: tmp_dividend
    TYPE(timedelta), POINTER  :: tmp_divisor
    TYPE(datetime), POINTER  :: tmp_ref, tmp_dt
    TYPE(t_datetime)                         :: dt_tmp
    INTEGER                                  :: errno

    tmp_dividend => newTimedelta(this%td, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    tmp_divisor => newTimedelta(divisor%td, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    IF (PRESENT(referenceDateTime)) THEN
      tmp_ref => newDateTime(referenceDateTime%dt, errno)
      CALL handle_errno(errno, &
        & __FILE__, &
        & __LINE__)
    END IF

    CALL divideTimeDeltaInSeconds(tmp_dividend, tmp_divisor, quotient, errno)

    IF (errno /= no_error) THEN
      IF (.NOT. PRESENT(referenceDateTime)) THEN
        CALL handle_errno(errno, &
          & __FILE__, &
          & __LINE__)
      ELSE
        dt_tmp = referenceDateTime + this
        tmp_dt => newDateTime(dt_tmp%dt, errno)
        CALL handle_errno(errno, &
          & __FILE__, &
          & __LINE__)
        CALL divideDatetimeDifferenceInSeconds(tmp_dt, tmp_ref, &
             &                                tmp_divisor, quotient, errno)
        CALL handle_errno(errno, &
          & __FILE__, &
          & __LINE__)

        CALL deallocateDatetime(tmp_dt)
      END IF
    END IF

    CALL deallocateTimedelta(tmp_dividend)
    CALL deallocateTimedelta(tmp_divisor)
    IF (PRESENT(referenceDateTime)) THEN
      CALL deallocateDatetime(tmp_ref)
    END IF
  END FUNCTION t_timedelta_divideInSecondsBy

  FUNCTION t_timedelta_toSeconds(this, td) RESULT(seconds)
    CLASS(t_timedelta), INTENT(in) :: this
    TYPE(t_datetime), INTENT(in) :: td
    INTEGER(c_int64_t)             :: seconds

    seconds = getTotalSecondsTimeDelta(this%td, td%dt)
  END FUNCTION t_timedelta_toSeconds

  FUNCTION t_timedelta_toMilliSeconds(this, td) RESULT(ms)
    CLASS(t_timedelta), INTENT(in) :: this
    TYPE(t_datetime), INTENT(in) :: td
    INTEGER(c_int64_t)             :: ms

    ms = getTotalMilliSecondsTimeDelta(this%td, td%dt)
  END FUNCTION t_timedelta_toMilliSeconds

  FUNCTION t_timedelta_toJulianDelta(this, dt) RESULT(jd)
    TYPE(t_juliandelta) :: jd
    CLASS(t_timedelta), INTENT(in) :: this
    TYPE(t_datetime), INTENT(in) :: dt
    TYPE(timedelta), POINTER    :: td_tmp
    TYPE(datetime), POINTER    :: dt_tmp
    TYPE(juliandelta), POINTER    :: jd_tmp
    INTEGER                        :: errno
    td_tmp => newTimedelta(this%td, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    dt_tmp => newDatetime(dt%dt, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    jd_tmp => newJuliandelta('+', 0_c_int64_t, 0_c_int64_t, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)

    CALL timeDeltaToJulianDelta(td_tmp, dt_tmp, jd_tmp)
    jd%jd = jd_tmp

    CALL deallocateTimedelta(td_tmp)
    CALL deallocateDatetime(dt_tmp)
    CALL deallocateJuliandelta(jd_tmp)
  END FUNCTION t_timedelta_toJulianDelta

  FUNCTION t_julianDay_getDay(this) RESULT(d)
    CLASS(t_julianday), INTENT(in) :: this
    INTEGER(c_int64_t) :: d
    d = this%jd%day
  END FUNCTION t_julianDay_getDay

  FUNCTION t_julianday_getFractionOfDayInMS(this) RESULT(ms)
    CLASS(t_julianday), INTENT(in) :: this
    INTEGER(c_int64_t) :: ms
    ms = this%jd%ms
  END FUNCTION t_julianday_getFractionOfDayInMS

  FUNCTION t_julianday_toDateTime(this) RESULT(res)
    TYPE(t_datetime) :: res
    CLASS(t_julianday)           :: this
    TYPE(julianday), POINTER    :: jd_tmp
    INTEGER                      :: errno
    TYPE(datetime), POINTER    :: dt_tmp

    jd_tmp => newJulianDay(this%jd%day, this%jd%ms, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)
    dt_tmp => newDatetime("1970-01-01T00:00:00", errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)

    CALL getDatetimeFromJulianDay(jd_tmp, dt_tmp, errno)
    CALL handle_errno(errno, &
      & __FILE__, &
      & __LINE__)

    res%dt = dt_tmp
    CALL deallocateJulianDay(jd_tmp)
    CALL deallocateDatetime(dt_tmp)
  END FUNCTION t_julianday_toDateTime

  FUNCTION t_datetime_min(a, b) RESULT(res)
    TYPE(t_datetime) :: a, b
    TYPE(t_datetime) :: res

    IF (a > b) THEN
      res = b
    ELSE
      res = a
    END IF
  END FUNCTION t_datetime_min

  FUNCTION t_datetime_max(a, b) RESULT(res)
    TYPE(t_datetime) :: a, b
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
    CHARACTER(c_char) ::c_sign

    c_sign = this%td%SIGN(1:1)
    c_pointer = my_newrawtimedelta(c_sign, INT(this%td%year, c_int64_t), this%td%month, this%td%day, &
         &                         this%td%hour, this%td%minute, this%td%second, this%td%ms)
    CALL handle_errno((.NOT. C_ASSOCIATED(c_pointer)), 5*100 + 2, &
      & __FILE__, &
      & __LINE__)
  END FUNCTION t_timedelta_get_c_pointer
