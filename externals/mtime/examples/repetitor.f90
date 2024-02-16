!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program repetitor
  
  use mtime

  implicit none

  character(len=64) :: case1 = 'R2/20130304T000000/20130504T030405'
  character(len=64) :: case2 = 'R3/8000101T230516/P10D'
  character(len=64) :: case3 = 'R/PT2M/-34560101T040404'
  character(len=64) :: case4 = 'R4/P1Y'
  character(len=64) :: case5 = 'R012/2013-03-04T00:00:00.546/2013-05-04T03:04:05'
  character(len=64) :: case6 = 'R3/800-01-01T23:05:16/P10D'
  character(len=64) :: case7 = 'R/PT2M/-347856-01-01T04:04:04'
  character(len=64) :: case8 = 'R4/P01Y'

  call setCalendar(proleptic_gregorian);

  call testCase(case1)
  call testCase(case2)
  call testCase(case3)
  call testCase(case4)
  call testCase(case5)
  call testCase(case6)
  call testCase(case7)
  call testCase(case8)

contains 
  
  subroutine testCase(caseX)
    character(len=*), intent(in) :: caseX
    
    character(len=MAX_REPETITION_STR_LEN) :: repetitor
    character(len=MAX_REPETITION_STR_LEN) :: start    
    character(len=MAX_REPETITION_STR_LEN) :: end      
    character(len=MAX_REPETITION_STR_LEN) :: duration 
    
    type(timedelta), pointer :: d
    type(datetime), pointer :: s, e
    
    character(len=MAX_DATETIME_STR_LEN) :: dstring
    character(len=MAX_TIMEDELTA_STR_LEN) :: tdstring
    
    logical :: lrepetitor, lstart, lend, lduration
    
    print *, 'Input ', caseX
    
    call splitRepetitionString(caseX, repetitor, start, end, duration, &
         &                            lrepetitor, lstart, lend, lduration);
    
    if (lrepetitor) then
      print *, 'Repetitor: ', trim(repetitor)
      print *, 'Repetitions ', getRepetitions(repetitor)
    endif
    if (lstart) then
      print *,'Start: ', trim(start)
      s => newDateTime(start)
      call datetimeToString(s, dstring)
      print *, 'Start rev: ', trim(dstring)
      call deallocateDateTime(s)
    endif
    if (lend) then
      print *, 'End: ', trim(end)
      e => newDateTime(end)
      call datetimeToString(e, dstring)
      print *, 'End rev: ', trim(dstring)
      call deallocateDateTime(e)
    endif
    if (lduration) then
      print *, 'Duration: ', trim(duration)
      d => newTimeDelta(duration)
      call timedeltaToString(d, tdstring)
      print *, 'Duration rev: ', trim(tdstring)
      call deallocateTimeDelta(d)
    endif

    print *, ''

  end subroutine testCase

end program repetitor
