!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM duration

  USE mtime
  USE mo_kind, ONLY: i8

  IMPLICIT NONE

  CALL setCalendar(proleptic_gregorian)

  CALL testTimedelta('PT12H', 'PT12H')
  CALL testTimedelta('PT06H', 'PT06H')
  CALL testTimedelta('PT6H', 'PT06H')
  CALL testTimedelta('PT2M', 'PT02M')
  CALL testTimedelta('PT01H', 'PT01H')
  CALL testTimedelta('P01M', 'P01M')
  CALL testTimedelta('P2M', 'P02M')
  CALL testTimedelta('PT10S', 'PT10.000S')
  CALL testTimedelta('PT1S', 'PT01.000S')
  CALL testTimedelta('PT100S', 'PT100S.000S')
  CALL testTimedelta('PT100.5S', 'PT100S.500S')
  CALL testTimedelta('PT934565S', 'PT934565.000S')
  CALL testTimedelta('P02DT06H', 'P02DT06H')
  CALL testTimedelta('P2DT6H', 'P02DT06H')
  CALL testTimedelta('P200D', 'P200D')

  CALL testTimedeltaComp('PT24M', 'PT06H00M00S', '<', .TRUE.)
  CALL testTimedeltaComp('PT43200S', 'PT18H00M00S', '<', .TRUE.)
  CALL testTimedeltaComp('PT345600S', 'P3DT18H12M01S', '>', .TRUE.)

  CALL testDatetimeDiffComp('2003-12-31T12:00:00', '2003-12-31T00:00:00', &
               &            '2004-01-01T00:00:00', '2003-12-31T00:00:00', 0.5)
  CALL testDatetimeDiffComp('2004-02-15T12:00:00', '2004-02-01T00:00:00', &
               &            '2004-03-01T00:00:00', '2004-02-01T00:00:00', 0.5)
CONTAINS

  SUBROUTINE testTimedeltaComp(td_A, td_B, op, expected)
    CHARACTER(len=*), INTENT(in) :: td_a
    CHARACTER(len=*), INTENT(in) :: td_b
    CHARACTER(len=*), INTENT(in) :: op
    LOGICAL, INTENT(in) :: expected

    TYPE(timedelta), POINTER   :: mtime_A, mtime_B

    mtime_A => newTimedelta(td_a)
    mtime_B => newTimedelta(td_b)

    SELECT CASE (op)
    CASE ('<')
      PRINT *, td_a, ' < ', td_b, ' ? ', mtime_A < mtime_B, ' expected ', expected
    CASE ('>')
      PRINT *, td_a, ' > ', td_b, ' ? ', mtime_A > mtime_B, ' expected ', expected
    CASE default
      PRINT *, 'Unknown comparison operator.'
    END SELECT

  END SUBROUTINE testTimedeltaComp

  SUBROUTINE testDatetimeDiffComp(dt1_dividend, dt2_dividend, &
      &                           dt1_divisor, dt2_divisor, expected)
    CHARACTER(len=*), INTENT(in) :: dt1_dividend
    CHARACTER(len=*), INTENT(in) :: dt2_dividend
    CHARACTER(len=*), INTENT(in) :: dt1_divisor
    CHARACTER(len=*), INTENT(in) :: dt2_divisor
    REAL, INTENT(in) :: expected

    INTEGER  :: ierr, ierrs
    INTEGER(i8)  :: denominator = 0

    TYPE(divisionquotienttimespan) :: quotient
    TYPE(datetime), POINTER   :: dt1_num => NULL()
    TYPE(datetime), POINTER   :: dt2_num => NULL()
    TYPE(datetime), POINTER   :: dt1_denom => NULL()
    TYPE(datetime), POINTER   :: dt2_denom => NULL()
    TYPE(timedelta), POINTER   :: td => NULL()

    ierrs = 0
    dt1_num => newDatetime(dt1_dividend, errno=ierr)
    ierrs = ierrs + ierr
    dt2_num => newDatetime(dt2_dividend, errno=ierr)
    ierrs = ierrs + ierr
    dt1_denom => newDatetime(dt1_divisor, errno=ierr)
    ierrs = ierrs + ierr
    dt2_denom => newDatetime(dt2_divisor, errno=ierr)
    ierrs = ierrs + ierr

    CALL divideTwoDatetimeDiffsInSeconds(dt1_num, dt2_num, dt1_denom, dt2_denom, denominator, quotient)
    PRINT *, '(', TRIM(dt1_dividend), ' - ', TRIM(dt2_dividend), ') / ',  &
      &      '(', TRIM(dt1_divisor), ' - ', TRIM(dt2_divisor), '):  ', REAL(quotient%remainder_in_ms)/1000./REAL(denominator)
    PRINT *, 'expected: ', expected

    CALL deallocateTimeDelta(td)

  END SUBROUTINE testDatetimeDiffComp

  SUBROUTINE testTimedelta(interval, expected)
    CHARACTER(len=*), INTENT(in) :: interval, expected
    CHARACTER(len=max_timedelta_str_len) :: dstring
    CHARACTER(len=max_mtime_error_str_len) :: estring
    TYPE(timedelta), POINTER :: d
    INTEGER :: error
    d => newTimedelta(TRIM(interval))
    error = 0
    CALL timedeltaToString(d, dstring, error)
    IF (error /= no_error) THEN
      CALL mtime_strerror(error, estring)
      PRINT *, 'ERROR: ', TRIM(estring), ' input: ', TRIM(interval)
    ELSE
      PRINT *, 'timedelta: input ', TRIM(interval), ' expected ', TRIM(expected), ': ', TRIM(dstring)
    END IF
    CALL deallocateTimedelta(d)
  END SUBROUTINE testTimedelta

END PROGRAM duration
