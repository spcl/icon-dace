!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM tas

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_long

  USE mtime

  IMPLICIT NONE

  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)

  CALL setCalendar(year_of_365_days)

  moduloTest: BLOCK

    TYPE(timedelta), POINTER :: td1 => NULL()
    TYPE(timedelta), POINTER :: td2 => NULL()
    INTEGER(i8) :: rem, quot

    td1 => newTimedelta('PT2H')
    td2 => newTimedelta('PT10M')

    rem = moduloTimedelta(td1, td2, quot)

    PRINT *, rem, quot

  END BLOCK moduloTest

  addSeconds: BLOCK

    TYPE(datetime), POINTER :: dti => NULL()
    TYPE(datetime), POINTER :: dto => NULL()
    INTEGER :: secs
    CHARACTER(len=max_datetime_str_len)  :: date_string

    dti => newdatetime('2014-03-04T01:40:00')

    CALL datetimeToString(dti, date_string)
    PRINT *, date_string

    secs = 600

    dto => datetimeaddseconds(dti, secs)

    CALL datetimeToString(dto, date_string)
    PRINT *, date_string

  END BLOCK addSeconds

  divideBySeconds: BLOCK

    TYPE(datetime), POINTER :: dt1 => NULL()
    TYPE(datetime), POINTER :: dt2 => NULL()

    CHARACTER(len=max_timedelta_str_len) :: ctd
    TYPE(timedelta), POINTER :: tdividend => NULL()
    TYPE(timedelta), POINTER :: tdivisor => NULL()

    TYPE(divisionquotienttimespan) :: tq

    dt1 => newdatetime('2014-03-04T00:00:00')
    dt2 => newdatetime('2014-03-07T00:00:00')

    tdividend => newtimedelta('PT0S')

    tdividend = dt2 - dt1

    CALL timedeltatostring(tdividend, ctd)
    PRINT *, ctd

    tdivisor => newtimedelta('PT600S')

    CALL timedeltatostring(tdivisor, ctd)
    PRINT *, ctd

    CALL divideTimeDeltaInSeconds(tdividend, tdivisor, tq)

    PRINT *, tq

  END BLOCK divideBySeconds

  compute_step_test: BLOCK

    TYPE(datetime), POINTER  :: start => NULL()
    TYPE(datetime), POINTER  :: END => NULL()
    TYPE(datetime), POINTER  :: current => NULL()

    TYPE(timedelta), POINTER :: delta => NULL()

    REAL :: dtime = 150.0

    INTEGER :: iadv_rcf = 4
    INTEGER :: step_offset = 0

    INTEGER :: step
    CHARACTER(len=max_datetime_str_len) :: exact_date

    CALL setCalendar(year_of_365_days)

    start => newDatetime('2014-01-01T00:00:00')
    END => newDatetime('2014-01-10T00:00:00')
    current => newDatetime('2014-01-03T12:40:00')

    delta => newtimedelta('PT600S')

    CALL compute_step(current, start, END, dtime, iadv_rcf, delta, step_offset, step, exact_date)

    PRINT *, 'Step: ', step
    PRINT *, 'Date: ', exact_date

  END BLOCK compute_step_test

CONTAINS

  FUNCTION datetimeaddseconds(refdt, intvlsec) RESULT(ret_datetime)
    TYPE(datetime), POINTER :: ret_datetime

    TYPE(datetime), POINTER :: refdt
    INTEGER, INTENT(in) :: intvlsec

    CHARACTER(len=max_timedelta_str_len) :: csec
    TYPE(timedelta), POINTER :: vlsec => NULL()

    CALL getptstringfromseconds(INT(intvlsec, c_long), csec)
    vlsec => newtimedelta(csec)

    ret_datetime => newDatetime("0000-01-01T00:00:00.000"); 
    ret_datetime = refdt + vlsec

    CALL deallocatetimedelta(vlsec)

  END FUNCTION datetimeaddseconds

  SUBROUTINE compute_step(mtime_current, mtime_begin, mtime_end, dtime, iadv_rcf, delta, step_offset, step, exact_date)

    TYPE(datetime), POINTER                         :: mtime_current       !< input date to translated into step
    TYPE(datetime), POINTER                         :: mtime_begin         !< begin of run (note: restart cases!)
    TYPE(datetime), POINTER                         :: mtime_end           !< end of run

    REAL, INTENT(in)  :: dtime               !< [s] length of a time step
    INTEGER, INTENT(in)  :: iadv_rcf            !< advection step: frequency ratio
    TYPE(timedelta), POINTER                         :: delta
    INTEGER, INTENT(in)  :: step_offset

    INTEGER, INTENT(out) :: step                !< result: corresponding simulations step
    CHARACTER(len=max_datetime_str_len), INTENT(out) :: exact_date          !< result: corresponding simulation date

    ! local variables

    INTEGER                              :: i
    INTEGER                              :: intvlsec

    TYPE(datetime), POINTER              :: mtime_step

    CHARACTER(len=max_datetime_str_len)  :: dt_string
    CHARACTER(len=max_timedelta_str_len) :: td_String

    TYPE(timedelta), POINTER             :: tddiff => NULL()

    TYPE(divisionquotienttimespan)      :: tq

    TYPE(timedelta), POINTER             :: vlsec => NULL()

    ! first, we compute the dynamic time step which is *smaller* than
    ! the desired date "mtime_date1"

    intvlsec = INT(dtime)
    CALL getptstringfromseconds(INT(intvlsec, c_long), td_string)
    vlsec => newtimedelta(td_string)

    tddiff => newtimedelta('PT0S')
    tddiff = mtime_current - mtime_begin

    CALL dividetimedeltainseconds(tddiff, vlsec, tq)
    step = INT(tq%quotient)

    mtime_step => newDatetime("0000-01-01T00:00:00.000"); 
    IF (step >= 0) THEN

      mtime_step = mtime_begin + step*vlsec

      ! starting from this step, we make (at most iadv_rcf) steps
      ! until we are *greater* than the desired date "mtime_date1" and
      ! we have reached an advection time step
      loop: DO i = 1, iadv_rcf
        !        if (ldebug) then
        !        dt_string = ''
        CALL datetimetostring(mtime_step, dt_string)
        WRITE (0, *) 'mtime_step = ', TRIM(dt_string)
        CALL datetimetostring(mtime_current, dt_string)
        WRITE (0, *) 'mtime_current = ', TRIM(dt_string)
        WRITE (0, *) ''
        !        end if

        IF ((mtime_step >= mtime_current) .AND. (MOD(step, iadv_rcf) == 0) .OR. (mtime_step == mtime_end)) THEN
          EXIT loop
        END IF

        mtime_step = mtime_step + delta
        step = step + 1

      END DO loop

      CALL datetimetostring(mtime_step, exact_date)
      CALL deallocatedatetime(mtime_step)

    END IF

    ! then we add the offset "jstep0" (nonzero for restart cases):

    step = step + step_offset

  END SUBROUTINE compute_step

END PROGRAM
