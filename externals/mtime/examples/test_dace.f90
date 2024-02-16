!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program test_mtime

  use mtime
  use mo_event_manager

  implicit none
  
  character(len=max_timedelta_str_len)   :: td_string
  character(len=max_datetime_str_len)    :: dt_string

  call setcalendar(proleptic_gregorian)
  call test1 ()
  call test2 ()
  
contains

  subroutine test1 ()

    type(timedelta),  pointer :: dt, t0, t1, t2, t3

    write (0,*) "test1:"

    dt => newtimedelta("pt24.000s")

    t0 => newtimedelta("pt3600.000s")

    t1 => newtimedelta("pt3600.000s")
    t2 => newtimedelta("pt60m")
    t3 => newtimedelta("pt1h")

    t1 =  t1 + dt
    t2 =  t2 + dt
    t3 =  t3 + dt

    call timedeltatostring (t0, td_string)
    write(0,*) "  t0 : ", td_string, "reference"
    call timedeltatostring (dt, td_string)
    write(0,*) "  dt : ", td_string, "modification"
    call timedeltatostring (t1, td_string)
    write(0,*) "  t1 : ", td_string, "???"
    call timedeltatostring (t2, td_string)
    write(0,*) "  t2 : ", td_string, "ok"
    call timedeltatostring (t3, td_string)
    write(0,*) "  t3 : ", td_string, "ok"
  end subroutine test1
  !-
  subroutine test2 ()

    type(event),      pointer :: mec_event      => null()

    type(datetime),   pointer :: mec_refdate    => null()

    type(eventgroup), pointer :: mec_eventgroup => null()

    type(datetime),   pointer :: mec_startdate  => null()
    type(datetime),   pointer :: mec_enddate    => null()

    type(timedelta),  pointer :: mec_start      => null()
    type(timedelta),  pointer :: mec_stop       => null()
    type(timedelta),  pointer :: mec_interval   => null()

    type(timedelta),  pointer :: time_step          => null()
    
    integer                   :: mec_events
    integer                   :: ierr
    logical                   :: lret

    type(datetime),   pointer :: mtime

    integer                   :: i

    ! offset for start from reference date
    mec_start      => newtimedelta("pt0s")
    ! offset for stop from reference date
    mec_stop       => newtimedelta("pt3600s")
    mec_interval   => newtimedelta("pt300s")

    time_step      => newtimedelta("pt24s")

    mec_refdate    => newdatetime ("2016-05-29t00:00:00.000")
    mec_startdate  => newdatetime (mec_refdate)
    mec_enddate    => newdatetime (mec_refdate)
    mec_startdate  =  mec_startdate + mec_start
    mec_enddate    =  mec_enddate   + mec_stop

    write (0,*)
    write (0,*) "test2:"
    call timedeltatostring (time_step,    td_string)
    write (0,*) "model time step   : ",   td_string
    call timedeltatostring (mec_start,    td_string)
    write (0,*) "mec start time    : ",   td_string
    call timedeltatostring (mec_stop,     td_string)
    write (0,*) "mec stop time     : ",   td_string
    call timedeltatostring (mec_interval, td_string)
    write (0,*) "mec interval      : ",   td_string

    write (0,*)

    call datetimetostring (mec_refdate,   dt_string)
    write (0,*) "mec reference date: ",   dt_string
    call datetimetostring (mec_startdate, dt_string)
    write (0,*) "mec start date    : ",   dt_string
    call datetimetostring (mec_enddate,   dt_string)
    write (0,*) "mec end date      : ",   dt_string

    write (0,*)

    write (0,*) "checking event management"

    call initeventmanager (mec_refdate)

    mec_events     =  addeventgroup('meceventgroup')
    mec_eventgroup => geteventgroup(mec_events)
    mec_enddate    =  mec_enddate   + time_step
    mec_event      => newevent('mec', mec_refdate, mec_startdate, mec_enddate, &
                                      mec_interval, errno=ierr)
    lret = addeventtoeventgroup(mec_event, mec_eventgroup)
    write (0,*) "addeventtoeventgroup returns:", lret
    call printeventgroup (mec_events)

    mtime => newdatetime (mec_startdate)
    i = 0
    do
      if (iscurrenteventactive (mec_event, mtime, plus_slack=time_step)) then
        i = i + 1
        call datetimetostring (mtime, dt_string)
        write (0,*) "mec will be called on: ", trim (dt_string)
      end if
      
      if (mtime >= mec_enddate) then
        exit
      end if
      mtime = mtime + time_step
    end do
    
    write(0,*) "check_dace_timer: total mec calls:", i, "(expected: 13)"
    
  end subroutine test2

end program test_mtime
