!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program tas

  use, intrinsic :: iso_c_binding, only: c_long

  use mtime

  implicit none
  
  integer, parameter :: i8 = selected_int_kind(14)
  
  call setCalendar(year_of_365_days)

  moduloTest: block

    type(timedelta), pointer :: td1 => null()
    type(timedelta), pointer :: td2 => null()
    integer(i8) :: rem, quot

    td1 => newTimedelta('PT2H')
    td2 => newTimedelta('PT10M')
    
    rem = moduloTimedelta(td1, td2, quot)

    print *, rem, quot
    
  end block moduloTest
  
  addSeconds: block 

    type(datetime), pointer :: dti => null()
    type(datetime), pointer :: dto => null()
    integer :: secs
    character(len=max_datetime_str_len)  :: date_string

    dti => newdatetime('2014-03-04T01:40:00')

    call datetimeToString(dti, date_string)
    print *, date_string

    secs = 600

    dto => datetimeaddseconds(dti, secs) 

    call datetimeToString(dto, date_string)
    print *, date_string

  end block addSeconds

  divideBySeconds: block

    type(datetime), pointer :: dt1 => null()
    type(datetime), pointer :: dt2 => null()

    character(len=max_timedelta_str_len) :: ctd
    type(timedelta), pointer :: tdividend => null()
    type(timedelta), pointer :: tdivisor => null()

    type(divisionquotienttimespan) :: tq 

    dt1 => newdatetime('2014-03-04T00:00:00')
    dt2 => newdatetime('2014-03-07T00:00:00')

    tdividend => newtimedelta('PT0S')

    tdividend = dt2 - dt1

    call timedeltatostring(tdividend, ctd)
    print *, ctd

    tdivisor => newtimedelta('PT600S')

    call timedeltatostring(tdivisor, ctd)
    print *, ctd

    call divideTimeDeltaInSeconds(tdividend, tdivisor, tq)

    print *, tq

  end block divideBySeconds

  compute_step_test: block

    type(datetime), pointer  :: start   => null()
    type(datetime), pointer  :: end     => null()
    type(datetime), pointer  :: current => null()

    type(timedelta), pointer :: delta   => null()

    real :: dtime = 150.0

    integer :: iadv_rcf    = 4
    integer :: step_offset = 0

    integer :: step
    character(len=max_datetime_str_len) :: exact_date

    call setCalendar(year_of_365_days)    

    start   => newDatetime('2014-01-01T00:00:00')
    end => newDatetime('2014-01-10T00:00:00')
    current => newDatetime('2014-01-03T12:40:00')

    delta => newtimedelta('PT600S')

    call compute_step(current, start, end, dtime, iadv_rcf, delta, step_offset, step, exact_date)

    print *, 'Step: ', step
    print *, 'Date: ', exact_date
    
  end block compute_step_test

contains

  function datetimeaddseconds(refdt, intvlsec) result(ret_datetime)
    type(datetime), pointer :: ret_datetime
    
    type(datetime), pointer :: refdt
    integer, intent(in) :: intvlsec
    
    character(len=max_timedelta_str_len) :: csec    
    type(timedelta), pointer :: vlsec => null()
    
    call getptstringfromseconds(int(intvlsec,c_long), csec)
    vlsec => newtimedelta(csec)
    
    ret_datetime => newDatetime("0000-01-01T00:00:00.000"); 
    
    ret_datetime = refdt + vlsec

    call deallocatetimedelta(vlsec)
    
  end function datetimeaddseconds
  
  subroutine compute_step(mtime_current, mtime_begin, mtime_end, dtime, iadv_rcf, delta, step_offset, step, exact_date)
    
    type(datetime),  pointer                         :: mtime_current       !< input date to translated into step
    type(datetime),  pointer                         :: mtime_begin         !< begin of run (note: restart cases!)
    type(datetime),  pointer                         :: mtime_end           !< end of run
    
    real,                                intent(in)  :: dtime               !< [s] length of a time step
    integer,                             intent(in)  :: iadv_rcf            !< advection step: frequency ratio
    type(timedelta), pointer                         :: delta
    integer,                             intent(in)  :: step_offset
    
    integer,                             intent(out) :: step                !< result: corresponding simulations step
    character(len=max_datetime_str_len), intent(out) :: exact_date          !< result: corresponding simulation date
    
    ! local variables

    integer                              :: i
    integer                              :: intvlsec
    
    type(datetime), pointer              :: mtime_step
    
    character(len=max_datetime_str_len)  :: dt_string
    character(len=max_timedelta_str_len) :: td_String    
    
    type(timedelta), pointer             :: tddiff => null()
    
    type(divisionquotienttimespan)      :: tq 
    
    type(timedelta), pointer             :: vlsec => null()
    
    ! first, we compute the dynamic time step which is *smaller* than
    ! the desired date "mtime_date1"
    
    intvlsec = int(dtime)
    call getptstringfromseconds(int(intvlsec,c_long), td_string)
    vlsec => newtimedelta(td_string)
    
    tddiff => newtimedelta('PT0S')
    tddiff = mtime_current - mtime_begin
    
    call dividetimedeltainseconds(tddiff, vlsec, tq)
    step = tq%quotient
    
    mtime_step => newDatetime("0000-01-01T00:00:00.000"); 
    
    if (step >= 0) then
    
      mtime_step = mtime_begin + step * vlsec

      ! starting from this step, we make (at most iadv_rcf) steps
      ! until we are *greater* than the desired date "mtime_date1" and
      ! we have reached an advection time step
      loop : do i = 1, iadv_rcf
        !        if (ldebug) then
        !        dt_string = ''
        call datetimetostring(mtime_step, dt_string)
        write (0,*) 'mtime_step = ', trim(dt_string)
        call datetimetostring(mtime_current, dt_string)
        write (0,*) 'mtime_current = ', trim(dt_string)
        write (0,*) ''
        !        end if
        
        if ((mtime_step >= mtime_current) .and. (mod(step, iadv_rcf) == 0) .or. (mtime_step == mtime_end)) then
          exit loop
        end if
        
        mtime_step = mtime_step + delta
        step = step + 1
        
      end do loop
      
      call datetimetostring(mtime_step, exact_date)
      call deallocatedatetime(mtime_step)

    end if
    
    ! then we add the offset "jstep0" (nonzero for restart cases):
    
    step        = step + step_offset
    
  end subroutine compute_step
  
end program
