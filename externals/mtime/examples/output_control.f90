!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program output_control

  use mtime

  implicit none

  type(eventgroup), pointer :: outputEventGroup => null()
  type(event), pointer :: outputEvent1 => null()
  type(event), pointer :: outputEvent2 => null()
  type(event), pointer :: outputEvent3 => null()
  type(event), pointer :: outputEvent4 => null()

  character(len=max_groupname_str_len) :: egstring

  character(len=*), parameter :: reference_date = '2000-01-01T00:00:00'
  character(len=*), parameter :: start_date = '2010-01-01T00:00:00'
  character(len=*), parameter :: end_date = '2010-02-01T00:00:00'

  integer :: error
  logical :: lret

  type event_namelist_handling
    ! represented in ISO 9601:2004 R string 
    character(len=32) :: steps
    character(len=32) :: start_date
    character(len=32) :: end_date
    character(len=32) :: duration
    ! additional model required numbers
    character(len=32) :: reference_date
    character(len=32) :: offset
    character(len=32) :: lbound_handling
    character(len=32) :: ubound_handling
  end type event_namelist_handling

  call setCalendar(proleptic_gregorian)

  outputEventGroup => newEventGroup('output driver')
  call getEventGroupName(outputEventGroup, egstring)
  print *, trim(egstring)

  outputEvent1 => newEvent('output1', reference_date, start_date, end_date, 'PT6H',errno=error)
  call checkError(error)
  lret = addEventToEventGroup(outputEvent1, outputEventGroup)

  outputEvent2 => newEvent('output2', reference_date, start_date, end_date, 'PT2H',errno=error)
  lret = addEventToEventGroup(outputEvent2, outputEventGroup)

  outputEvent3 => newEvent('output3', reference_date, start_date, end_date, 'PT4H')
  lret = addEventToEventGroup(outputEvent3, outputEventGroup)

  outputEvent4 => newEvent('output4', reference_date, start_date, end_date, 'PT10H')
  lret = addEventToEventGroup(outputEvent4, outputEventGroup)

  call printEventTriggerTimes(outputEvent1)
  call printEventTriggerTimes(outputEvent2)
  call printEventTriggerTimes(outputEvent3)
  call printEventTriggerTimes(outputEvent4)

contains

  subroutine printEventTriggerTimes(eventI)
    type(event), pointer :: eventI
    type(datetime), pointer :: cd => null()
    type(datetime), pointer :: sd => null()
    type(datetime), pointer :: ed => null()
    type(timedelta), pointer :: dt => null()
    character(len=max_datetime_str_len) :: tstring
    character(len=max_timedelta_str_len) :: dstring
    integer :: i = 0
    sd => getEventFirstDateTime(eventI)
    call datetimeToString(sd, tstring)
    print *, 'Start time    : ', trim(tstring)
    ed => getEventLastDateTime(eventI)
    call datetimeToString(ed, tstring)
    print *, 'End time      : ', trim(tstring)
    dt => getEventInterval(eventI)
    call timedeltaToString(dt, dstring)
    print *, 'Time step     : ', trim(dstring)
    cd => newDatetime(sd)
    call datetimeToString(cd, dstring)
    print *, 'Current time  : ', trim(dstring)
    output_times_loop: do 
      ! open lower bound set
      cd = cd + dt
      ! open upper bound set
      if (cd > ed) exit output_times_loop
      call datetimeToString(cd, dstring)
      print '(a,i3,a,a)', 'Output step ', i, ' time: ', trim(dstring)
      ! closed lower bound set
      !cd = cd + dt
      i = i + 1
      ! closed upper bound set
      ! if (cd > ed) exit output_times_loop
    enddo output_times_loop
    print *,'________________________________________________________'
    print *,''
    call deallocateDatetime(sd)
    call deallocateDatetime(ed)
    call deallocateDatetime(cd)
    call deallocateTimedelta(dt)
  end subroutine printEventTriggerTimes

  subroutine checkError(errno)
    integer, intent(in) :: errno
    character(len=max_mtime_error_str_len) :: estring
    if (errno /= no_error) then
      call mtime_strerror(errno, estring)
      print *, trim(estring)
      stop 'ERROR: finish.'
    endif
  end subroutine checkError


end program output_control
