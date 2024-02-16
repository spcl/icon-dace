!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program simulate_iau

  use mtime
  use mtime_hl

  implicit none

  type(t_datetime)  :: start_date, stop_date, current_date, previous_date
  type(t_timedelta) :: time_step, iau_time_shift

  logical, parameter :: iterate_iau = .TRUE.
  integer :: iau_iter

  real :: dt_shift, dtime

  integer :: jstep, jstep0, jstep_shift

  write (0,*) "Start execution, set calendar ..."

  call setCalendar(proleptic_gregorian)

  write (0,*) "Assign values ..."

  start_date = t_datetime("2016-01-01T00:00:00")
  stop_date  = t_datetime("2016-01-02T00:00:00")

  time_step = t_timedelta("PT15M")
  dtime = 900.0

  iau_time_shift = t_timedelta("-PT1H30M")
  dt_shift = -5400.0

  write (0,*) "Prepare time loop ..."

  current_date = start_date
  current_date = current_date + iau_time_shift

  write (0,*) '           start date ', start_date%toString()
  write (0,*) '   shifted start date ', current_date%toString()

  iau_iter = MERGE(1, 0, iterate_iau)

  if (iterate_iau) then
    jstep_shift = nint(dt_shift/dtime)
  else
    jstep_shift = 0
  endif

  previous_date = current_date

  jstep0 = 0
  jstep = (jstep0 + 1) + jstep_shift

  write (0,*) "Start time loop ..."

  time_loop: do

    current_date = current_date + time_step

    write (0,*) "   Time loop ", current_date%toString(), jstep

    if ( (current_date%getDay() /= previous_date%getDay()) .and. .not. (jstep == 0 .and. iau_iter == 1) ) then
      previous_date = current_date
    endif

    write (0,*) '   --- integrate nh input: ', current_date%toString(), jstep-jstep_shift, iau_iter

    if (current_date >= stop_date) then
      exit time_loop
    end if

    if (jstep == 0 .and. iau_iter == 1) then
      iau_iter = 2
      jstep = (jstep0 + 1) + jstep_shift
      current_date = current_date + iau_time_shift
    else
      jstep = jstep + 1
    endif

  enddo time_loop

end program simulate_iau
