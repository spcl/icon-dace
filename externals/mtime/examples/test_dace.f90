!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM test_mtime

  USE mtime
  USE mo_event_manager

  IMPLICIT NONE

  CHARACTER(len=max_timedelta_str_len)   :: td_string
  CHARACTER(len=max_datetime_str_len)    :: dt_string

  CALL setcalendar(proleptic_gregorian)
  CALL test1()
  CALL test2()

CONTAINS

  SUBROUTINE test1()

    TYPE(timedelta), POINTER :: dt, t0, t1, t2, t3

    WRITE (0, *) "test1:"

    dt => newtimedelta("pt24.000s")

    t0 => newtimedelta("pt3600.000s")

    t1 => newtimedelta("pt3600.000s")
    t2 => newtimedelta("pt60m")
    t3 => newtimedelta("pt1h")

    t1 = t1 + dt
    t2 = t2 + dt
    t3 = t3 + dt

    CALL timedeltatostring(t0, td_string)
    WRITE (0, *) "  t0 : ", td_string, "reference"
    CALL timedeltatostring(dt, td_string)
    WRITE (0, *) "  dt : ", td_string, "modification"
    CALL timedeltatostring(t1, td_string)
    WRITE (0, *) "  t1 : ", td_string, "???"
    CALL timedeltatostring(t2, td_string)
    WRITE (0, *) "  t2 : ", td_string, "ok"
    CALL timedeltatostring(t3, td_string)
    WRITE (0, *) "  t3 : ", td_string, "ok"
  END SUBROUTINE test1
  !-
  SUBROUTINE test2()

    TYPE(event), POINTER :: mec_event => NULL()

    TYPE(datetime), POINTER :: mec_refdate => NULL()

    TYPE(eventgroup), POINTER :: mec_eventgroup => NULL()

    TYPE(datetime), POINTER :: mec_startdate => NULL()
    TYPE(datetime), POINTER :: mec_enddate => NULL()

    TYPE(timedelta), POINTER :: mec_start => NULL()
    TYPE(timedelta), POINTER :: mec_stop => NULL()
    TYPE(timedelta), POINTER :: mec_interval => NULL()

    TYPE(timedelta), POINTER :: time_step => NULL()

    INTEGER                   :: mec_events
    INTEGER                   :: ierr
    LOGICAL                   :: lret

    TYPE(datetime), POINTER :: mtime

    INTEGER                   :: i

    ! offset for start from reference date
    mec_start => newtimedelta("pt0s")
    ! offset for stop from reference date
    mec_stop => newtimedelta("pt3600s")
    mec_interval => newtimedelta("pt300s")

    time_step => newtimedelta("pt24s")

    mec_refdate => newdatetime("2016-05-29t00:00:00.000")
    mec_startdate => newdatetime(mec_refdate)
    mec_enddate => newdatetime(mec_refdate)
    mec_startdate = mec_startdate + mec_start
    mec_enddate = mec_enddate + mec_stop

    WRITE (0, *)
    WRITE (0, *) "test2:"
    CALL timedeltatostring(time_step, td_string)
    WRITE (0, *) "model time step   : ", td_string
    CALL timedeltatostring(mec_start, td_string)
    WRITE (0, *) "mec start time    : ", td_string
    CALL timedeltatostring(mec_stop, td_string)
    WRITE (0, *) "mec stop time     : ", td_string
    CALL timedeltatostring(mec_interval, td_string)
    WRITE (0, *) "mec interval      : ", td_string

    WRITE (0, *)

    CALL datetimetostring(mec_refdate, dt_string)
    WRITE (0, *) "mec reference date: ", dt_string
    CALL datetimetostring(mec_startdate, dt_string)
    WRITE (0, *) "mec start date    : ", dt_string
    CALL datetimetostring(mec_enddate, dt_string)
    WRITE (0, *) "mec end date      : ", dt_string

    WRITE (0, *)

    WRITE (0, *) "checking event management"

    CALL initeventmanager(mec_refdate)

    mec_events = addeventgroup('meceventgroup')
    mec_eventgroup => geteventgroup(mec_events)
    mec_enddate = mec_enddate + time_step
    mec_event => newevent('mec', mec_refdate, mec_startdate, mec_enddate, &
                          mec_interval, errno=ierr)
    lret = addeventtoeventgroup(mec_event, mec_eventgroup)
    WRITE (0, *) "addeventtoeventgroup returns:", lret
    CALL printeventgroup(mec_events)

    mtime => newdatetime(mec_startdate)
    i = 0
    DO
      IF (iscurrenteventactive(mec_event, mtime, plus_slack=time_step)) THEN
        i = i + 1
        CALL datetimetostring(mtime, dt_string)
        WRITE (0, *) "mec will be called on: ", TRIM(dt_string)
      END IF

      IF (mtime >= mec_enddate) THEN
        EXIT
      END IF
      mtime = mtime + time_step
    END DO

    WRITE (0, *) "check_dace_timer: total mec calls:", i, "(expected: 13)"

  END SUBROUTINE test2

END PROGRAM test_mtime
