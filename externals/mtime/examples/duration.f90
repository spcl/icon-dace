!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program duration
    
  use mtime
  USE mtime_error_handling
  use mo_kind , only: i8

  implicit none

  call setCalendar(proleptic_gregorian)

  call testTimedelta('PT12H', 'PT12H')
  call testTimedelta('PT06H', 'PT06H')
  call testTimedelta('PT6H', 'PT06H')
  call testTimedelta('PT2M', 'PT02M')
  call testTimedelta('PT01H', 'PT01H')
  call testTimedelta('P01M', 'P01M')
  call testTimedelta('P2M', 'P02M')
  call testTimedelta('PT10S', 'PT10.000S')
  call testTimedelta('PT1S', 'PT01.000S')
  call testTimedelta('PT100S', 'PT100S.000S')
  call testTimedelta('PT100.5S', 'PT100S.500S')
  call testTimedelta('PT934565S', 'PT934565.000S')
  call testTimedelta('P02DT06H', 'P02DT06H')
  call testTimedelta('P2DT6H', 'P02DT06H')
  call testTimedelta('P200D', 'P200D')

  call testTimedeltaComp('PT24M', 'PT06H00M00S', '<', .TRUE.)
  call testTimedeltaComp('PT43200S', 'PT18H00M00S', '<', .TRUE.)
  call testTimedeltaComp('PT345600S', 'P3DT18H12M01S', '>', .TRUE.) 

  call testDatetimeDiffComp('2003-12-31T12:00:00','2003-12-31T00:00:00', &
               &            '2004-01-01T00:00:00','2003-12-31T00:00:00',0.5)
  call testDatetimeDiffComp('2004-02-15T12:00:00','2004-02-01T00:00:00', &
               &            '2004-03-01T00:00:00','2004-02-01T00:00:00',0.5)
contains

  subroutine testTimedeltaComp(td_A, td_B, op, expected)
    character(len=*), intent(in) :: td_a
    character(len=*), intent(in) :: td_b
    character(len=*), intent(in) :: op
    logical, intent(in) :: expected

    TYPE(timedelta),     POINTER   :: mtime_A, mtime_B

    mtime_A => newTimedelta(td_a)
    mtime_B => newTimedelta(td_b)

    select case (op)
    case ('<')
      print *, td_a, ' < ', td_b, ' ? ', mtime_A < mtime_B, ' expected ', expected      
    case ('>')
      print *, td_a, ' > ', td_b, ' ? ', mtime_A > mtime_B, ' expected ', expected      
    case default
      print *, 'Unknown comparison operator.'
    end select
    
  end subroutine testTimedeltaComp
  
  subroutine testDatetimeDiffComp(dt1_dividend, dt2_dividend, &
      &                           dt1_divisor, dt2_divisor, expected)
    character(len=*), intent(in) :: dt1_dividend
    character(len=*), intent(in) :: dt2_dividend
    character(len=*), intent(in) :: dt1_divisor
    character(len=*), intent(in) :: dt2_divisor
    REAL, intent(in) :: expected

    integer  :: ierr, ierrs
    integer(i8)  :: denominator = 0

    TYPE(divisionquotienttimespan) :: quotient
    TYPE(datetime),     POINTER   :: dt1_num => NULL()
    TYPE(datetime),     POINTER   :: dt2_num => NULL()
    TYPE(datetime),     POINTER   :: dt1_denom => NULL()
    TYPE(datetime),     POINTER   :: dt2_denom => NULL()
    TYPE(timedelta),    POINTER   :: td => NULL()


    ierrs = 0
    dt1_num => newDatetime(dt1_dividend, errno = ierr)
    ierrs = ierrs + ierr
    dt2_num => newDatetime(dt2_dividend, errno = ierr)
    ierrs = ierrs + ierr
    dt1_denom => newDatetime(dt1_divisor, errno = ierr)
    ierrs = ierrs + ierr
    dt2_denom => newDatetime(dt2_divisor, errno = ierr)
    ierrs = ierrs + ierr

    CALL divideTwoDatetimeDiffsInSeconds(dt1_num,dt2_num,dt1_denom,dt2_denom,denominator,quotient)
    print *, '(', TRIM(dt1_dividend), ' - ', TRIM(dt2_dividend), ') / ',  &
      &      '(', TRIM(dt1_divisor), ' - ', TRIM(dt2_divisor), '):  ', REAL(quotient%remainder_in_ms) / 1000. / REAL(denominator)
    print *, 'expected: ', expected

    CALL deallocateTimeDelta(td)

  end subroutine testDatetimeDiffComp
  
  subroutine testTimedelta(interval, expected)
    character(len=*), intent(in) :: interval, expected
    character(len=max_timedelta_str_len) :: dstring
    character(len=max_mtime_error_str_len) :: estring
    type(timedelta), pointer :: d
    integer :: error 
    d => newTimedelta(trim(interval))
    error = 0
    call timedeltaToString(d, dstring, error)
    if (error /= no_error) then
      call mtime_strerror(error, estring)
      print *, 'ERROR: ', trim(estring), ' input: ', trim(interval)
    else
      print *, 'timedelta: input ', trim(interval), ' expected ', trim(expected), ': ', trim(dstring)
    endif
    call deallocateTimedelta(d)
  end subroutine testTimedelta

end program duration
