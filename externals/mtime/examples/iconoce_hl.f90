!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM event_testing

#ifndef __NVCOMPILER

  USE mtime, ONLY: setcalendar, proleptic_gregorian
  USE mtime_hl

  IMPLICIT NONE

  TYPE(t_eventGroup) :: restart_events

  TYPE(t_event) :: restart
  TYPE(t_event) :: checkpoint

  TYPE(t_datetime)  :: exp_start_date, exp_end_date, exp_ref_date
  TYPE(t_timedelta) :: restart_interval, checkpoint_interval

  CALL setcalendar(proleptic_gregorian)

  exp_start_date = t_datetime("2016-08-01T00:00:00")
  exp_end_date = t_datetime("2016-09-10T00:00:00")

  exp_ref_date = exp_start_date

  restart_interval = t_timedelta("P1D")
  checkpoint_interval = t_timedelta("PT1H")

  restart_events = t_eventGroup("restart events")

  restart = t_event("restart", exp_ref_date, exp_start_date, exp_end_date, restart_interval)
  CALL restart_events%append(restart)
  checkpoint = t_event("checkpoint", exp_ref_date, exp_start_date, exp_end_date, checkpoint_interval)
  CALL restart_events%append(checkpoint)

#endif

END PROGRAM event_testing

! module mo_event_manager

!   use mtime_hl
!   use mtime, only: max_groupname_str_len, &
!        &           max_eventname_str_len, &
!        &           max_repetition_str_len

!   implicit none

!   private

!   public :: initEventManager
!   public :: getModelReferenceDate
!   public :: addEventGroup
!   public :: getEventGroup
!   public :: printEventGroup
!   public :: getEventComponents

!   type event_group_list
!     type(t_event_group), pointer :: group
!   end type event_group_list

!   type(event_group_list), allocatable :: model_event_groups(:)

!   integer :: model_event_groups_list_member
!   integer :: model_event_groups_list_size = 16

!   type(t_datetime) :: model_reference_date

!   logical :: linitialized = .false.

! contains

!   subroutine initEventManager(referenceDate)
!     type(t_datetime), intent(in) :: referenceDate

!     model_reference_date = t_datetime(referenceDate)

!     allocate(model_event_groups(model_event_groups_list_size))
!     model_event_groups_list_member = 0

!     linitialized = .true.

!   end subroutine initEventManager

!   function getModelReferenceDate() result(r)
!     type(t_datetime) :: r

!     if (linitialized) then
!       r = model_reference_date
!     endif

!   end function getModelReferenceDate

!   function addEventGroup(group) result(handle)
!     integer :: handle
!     character(len=*), intent(in) :: group
!     type(event_group_list), allocatable :: tmp(:)
!     character(len=max_groupname_str_len) :: gstring

!     if (.not. linitialized) then
!       print *, 'ERROR: event manager not initialized.'
!       stop
!     endif

!     if (model_event_groups_list_member == model_event_groups_list_size) then
!       print *, 'INFO: reallocating event group list.'
!       allocate(tmp(model_event_groups_list_size))
!       tmp(:) = model_event_groups(:)
!       deallocate(model_event_groups)
!       allocate(model_event_groups(2*model_event_groups_list_size))
!       model_event_groups(:model_event_groups_list_size) = tmp(:)
!       deallocate(tmp)
!       model_event_groups_list_size = 2*model_event_groups_list_size
!       print *, 'INFO: new evcent group list size: ', model_event_groups_list_size
!     endif

!     model_event_groups_list_member = model_event_groups_list_member + 1

!     model_event_groups(model_event_groups_list_member)%group = t_event_group(group)
!     call getEventGroupName(model_event_groups(model_event_groups_list_member)%group, gstring)
!     print *, 'INFO: Added event group: ', trim(gstring)

!     handle = model_event_groups_list_member

!   end function addEventGroup

!   function getEventGroup(handle) result(eventGroupListMember)
!     type(t_event_group) :: eventGroupListMember
!     integer, intent(in) :: handle
!     if (handle <= model_event_groups_list_member) then
!        eventGroupListMember =  model_event_groups(handle)%group
!      endif
!   end function getEventGroup

!   subroutine printEventGroup(handle)
!     integer, intent(in) :: handle
!     type(t_event_group) :: currentEventGroup
!     type(t_event) :: currentEvent
!     character(len=max_eventname_str_len) :: estring
!     character(len=max_groupname_str_len) :: egstring

!     currentEventGroup => getEventGroup(handle)
!     call getEventGroupName(currentEventGroup, egstring)

!     currentEvent => getFirstEventFromEventGroup(model_event_groups(handle)%group)

!     print *, 'Event list: ', trim(egstring)
!     do while (associated(currentEvent))
!       call eventToString(currentEvent, estring)
!       print *,'   event ', trim(estring)
!       currentEvent = currentEvent%getNextEvent()
!     enddo
!   end subroutine printEventGroup

!   subroutine getEventComponents(eventString, referenceDate, timeInterval, startDate, endDate)
!     character(len=max_repetition_str_len), intent(in) :: eventString
!     type(t_datetime) :: referenceDate
!     type(t_timedelta :: timeInterval
!     type(t_datetime) :: startDate
!     type(t_datetime) :: endDate

!     character(len=max_repetition_str_len) :: r, s, e, d
!     logical :: lr, ls, le, ld

!     call splitRepetitionString(eventString, r, s, e, d, lr, ls, le, ld)

!     if (lr) then
!       if (getRepetitions(r) /= -1) then
!         print *, 'WARNING: event setup should not have explicit repeat count.'
!       endif
!     endif

!     if (ls) then
!       startDate = t_datetime(s)
!     endif

!     if (le) then
!       endDate = t_datetime(e)
!     endif

!     if (ld) then
!       timeInterval = t_timeDelta(d)
!     else
!       print *, 'ERROR: time interval should be given.'
!       stop
!     endif

!   end subroutine getEventComponents

! end module mo_event_manager
! !__________________________________________________________________________________________________
! !
! program iconoce

!   use mtime_hl
!   use mo_event_manager

!   implicit none

!   type(t_datetime) :: experiment_reference_date

!   type(t_datetime) :: checkpoint_reference_date

!   type(t_datetime) :: experiment_start_date
!   type(t_datetime) :: experiment_end_date

!   type(t_datetime) :: start_date
!   type(t_datetime) :: stop_date

!   type(t_datetime) :: next_checkpoint_date
!   type(t_datetime) :: next_restart_date

!   type(t_timedelta) :: dynamic_time_step

!   type(t_timedelta) :: checkpoint_time_step
!   type(t_timedelta) :: restart_time_step

!   type(t_timedelta) :: coupling_time_step
!   type(t_timedelta) :: model_time_step

!   type(t_datetime) :: current_date

!   type(t_event_group) :: outputEventGroup

!   type(t_event) :: checkpointEvent
!   type(t_event) :: restartEvent

!   character(len=max_calendar_str_len)  :: calendar_in_use
!   character(len=max_datetime_str_len)  :: dstring
!   character(len=max_timedelta_str_len) :: tdstring

!   character(len=max_mtime_error_str_len) :: errstring
!   ! namelist variables

!   character(len=max_calendar_str_len) :: calendar

!   character(len=max_datetime_str_len) :: experimentReferenceDate = ''
!   character(len=max_datetime_str_len) :: experimentStartDate = ''
!   character(len=max_datetime_str_len) :: experimentEndDate

!   character(len=max_datetime_str_len) :: startDate

!   character(len=max_timedelta_str_len) :: modelTimeStep

!   character(len=max_repetition_str_len) :: checkpointTimeInterval
!   character(len=max_repetition_str_len) :: restartTimeInterval

!   character(len=max_repetition_str_len) :: couplingTimeInterval

!   character(len=132) :: error_message

!   integer :: iunit, icalendar, ierror
!   integer :: outputEvents
!   logical :: lret, isRestart, isRestartTimeRelative

!   character(len=max_groupname_str_len) :: egstring

!   type(t_datetime) :: checkpointRefDate
!   type(t_datetime) :: checkpointStartDate
!   type(t_datetime) :: checkpointEndDate
!   type(t_timedelta) :: checkpointInterval

!   type(t_datetime) :: restartRefDate
!   type(t_datetime) :: restartStartDate
!   type(t_datetime) :: restartEndDate
!   type(t_timedelta) :: restartInterval

!   !________________________________________________________________________________________________
!   !

!   namelist /timeControl/ &
!        &    calendar, &
!        &    experimentReferenceDate, &
!        &    experimentStartDate, &
!        &    experimentEndDate, &
!        &    modelTimeStep, &
!        &    checkpointTimeInterval, &
!        &    restartTimeInterval, &
!        &    couplingTimeInterval, &
!        &    isRestart, &
!        &    isRestartTimeRelative

!   open (file='examples/iconoce.nml', newunit=iunit, iostat=ierror)
!   if (ierror /= 0) then
!     print *, 'ERROR: could not open namelist file.'
!     stop
!   else
!     read (unit=iunit, nml=timeControl, iostat=ierror, iomsg=error_message)
!     if (ierror /= 0) then
!       print *, 'ERROR: could not read namelist file.'
!       print *, '       ', trim(error_message)
!       stop
!     endif
!     close (unit=iunit)
!   endif

!   !________________________________________________________________________________________________
!   !

!   select case (toLower(calendar))
!   case ('proleptic gregorian')
!     icalendar  = proleptic_gregorian
!   case ('365 day year')
!     icalendar = year_of_365_days
!   case ('360 day year')
!     icalendar = year_of_360_days
!   case default
!     icalendar = calendar_not_set
!     print *, 'ERROR: calendar ', trim(calendar), ' not available/unknown.'
!     stop
!   end select

!   call setCalendar(icalendar)
!   call calendarToString(calendar_in_use)
!   print *, 'Calendar: ', trim(calendar_in_use)

!   print *, ''

!   !________________________________________________________________________________________________
!   !

!   if (experimentReferenceDate /= '') then
!     experiment_reference_date => newDatetime(experimentReferenceDate)
!   endif

!   if (isRestart) then
!     call readRestart(startDate)
!     start_date => newDatetime(startDate)
!   else
!     start_date => newDatetime(experimentStartDate)
!   endif

!   experiment_start_date => newDatetime(experimentStartDate)

!   if (isRestartTimeRelative) then
!     checkpoint_reference_date => newDatetime(start_date)
!   else
!     checkpoint_reference_date => newDatetime(experiment_reference_date)
!   endif

!   if (associated(experiment_reference_date)) then
!     call initEventManager(experiment_reference_date)
!   else
!     call initEventManager(experiment_start_date)
!   endif

!   experiment_reference_date => getModelReferenceDate()

!   print *, 'Experiment reference date: ', experiment_reference_date%toString()
!   print *, 'Experiment start date    : ', experiment_start_date%toString()

!   experiment_end_date = t_datetime(experimentEndDate)
!   print *, 'Experiment end date      : ', experiment_end_date%toString()

!   print *, ''

!   !________________________________________________________________________________________________
!   !
!   ! event_group_setup: block

!   outputEvents =  addEventGroup('outputEventGroup')
!   outputEventGroup => getEventGroup(outputEvents)
!   print *, 'output event group handler: ', outputEvents
!   call getEventGroupName(outputEventGroup, egstring)
!   print *, 'output event group name   : ', trim(egstring)
!   print *, ''

!   ! end block event_group_setup
!   !________________________________________________________________________________________________
!   !
!   ! checkpoint_restart_time_intervals: block

!   checkpointRefDate   => checkpoint_reference_date
!   checkpointStartDate => experiment_start_date
!   checkpointEndDate   => experiment_end_date
!   call getEventComponents(checkpointTimeInterval, checkpointRefDate, &
!        &                  checkpointInterval, checkpointStartDate, checkpointEndDate)
!   checkpointEvent => newEvent('checkpoint', checkpointRefDate, &
!        &                      checkpointStartDate, checkpointEndDate, checkpointInterval, &
!        &                      errno=ierror)
!   if (ierror /= no_Error) then
!     CALL mtime_strerror(ierror, errstring)
!     print *, 'ERROR: ', trim(errstring)
!     stop
!   endif
!   lret = addEventToEventGroup(checkpointEvent, outputEventGroup)

!   restartRefDate   => checkpoint_reference_date
!   restartStartDate => experiment_start_date
!   restartEndDate   => experiment_end_date
!   call getEventComponents(restartTimeInterval, restartRefDate, &
!        &                  restartInterval, restartStartDate, restartEndDate)
!   restartEvent => newEvent('restart', restartRefDate, &
!        &                   restartStartDate, restartEndDate, restartInterval, &
!        &                   errno=ierror)
!   if (ierror /= no_Error) then
!     CALL mtime_strerror(ierror, errstring)
!     print *, 'ERROR: ', trim(errstring)
!     stop
!   endif
!   lret = addEventToEventGroup(restartEvent, outputEventGroup)

!   ! end block checkpoint_restart_time_intervals

!   call printEventGroup(outputEvents)

!   !________________________________________________________________________________________________
!   !

!   print *, ''
!   model_time_step => newTimedelta(modelTimeStep)
!   call timedeltaToString(model_time_step, tdstring)
!   print *, 'Dynamics (basic model) time step: ', trim(tdstring)
!   print *, ''

!   !________________________________________________________________________________________________
!   !

!   current_date => newDatetime(start_date)
!   stop_date => newDatetime(start_date)
!   stop_date = stop_date + getEventInterval(restartEvent)

!   !________________________________________________________________________________________________
!   !
!   ! check_time_interval_consistency: block
!   !..............................................................................................
!   ! 1. check, if restart is in the experiments time interval
!   !
!   if (stop_date > experiment_end_date) then
!     print *, 'WARNING: run would not create a restart file.'
!     print *, '         Reset the stop_date to experiment end_date.'
!     stop_date => experiment_end_date
!   endif
!   !..............................................................................................
!   ! 2. check, if checkpointing is
!   !
!   next_checkpoint_date => newDatetime(start_date)
!   next_checkpoint_date = next_checkpoint_date + getEventInterval(checkpointEvent)
!   call datetimeToString(next_checkpoint_date, dstring)
!   print *, 'First checkpoint date: ', trim(dstring)
!   !..............................................................................................
!   ! 3. check, if restarting is
!   !
!   next_restart_date => newDatetime(start_date)
!   next_restart_date = next_restart_date + getEventInterval(restartEvent)
!   call datetimeToString(next_restart_date, dstring)
!   print *, 'First restart date: ', trim(dstring)
!   !..............................................................................................
!   !end block check_time_interval_consistency
!   !________________________________________________________________________________________________
!   !

!   call datetimeToString(current_date, dstring)
!   print *, 'Model date starting the time integration loop: ', trim(dstring)

!   time_integration: do
!     !............................................................................................
!     ! print date and time
!     call datetimeToString(current_date, dstring)
! !    print *, 'Model time loop  : ', trim(dstring)
!     !............................................................................................
!     ! initiate restart
!     if ((isCurrentEventActive(restartEvent, current_date) .and. start_date /= current_date) &
!          .or. current_date == experiment_end_date) then
!       print *, 'INFO: write restart.'
!       call writeRestart(current_date)
!       exit time_integration
!     endif
!     !............................................................................................
!     ! initiate checkpoint, we do not checkpoint/restart
!     if (isCurrentEventActive(checkpointEvent, current_date) .and. start_date /= current_date) then
!       print *, 'INFO: write checkpoint.'
!       call writeRestart(current_date)
!     endif
!     !............................................................................................
!     ! calculate next date and time
!     current_date = current_date + model_time_step
!     !............................................................................................
!     ! if new date and time is larger than end of run exit time integration: should never hit
!     if (current_date > stop_date) exit time_integration
!   enddo time_integration

!   call datetimeToString(current_date, dstring)
!   print *, 'Model date leaving the time integration loop : ', trim(dstring)

!   !________________________________________________________________________________________________
!   !

! contains

!   !________________________________________________________________________________________________
!   !
!   subroutine writeRestart(currentDate)
!     type(t_datetime), pointer :: currentDate
!     character(len=max_datetime_str_len)  :: dstring
!     character(len=max_datetime_str_len+16)  :: filename

!     integer :: iunit, ierror

!     call datetimeToString(currentDate, dstring)
!     print *, '      Writing restart/checkpoint file for ', trim(dstring)

!     write (filename,'(a,a,a)') 'restart_oce_', trim(dstring), '.dat'

!     open (file='examples/'//trim(filename), newunit=iunit, iostat=ierror)
!     if (ierror /= 0) then
!       print *, 'ERROR: could not open restart file for writing.'
!       stop
!     else
!       write (unit=iunit, iostat=ierror, iomsg=error_message, fmt='(a,a)') &
!            & 'restart: ', trim(dstring)
!       if (ierror /= 0) then
!         print *, 'ERROR: could not write restart/checkpoint file.'
!         print *, '       ', trim(error_message)
!         stop
!       else
!         close (unit=iunit)
!         !..........................................................................................
!         !
!         ierror = 0
!         error_message = ''
!         call execute_command_line('rm -f examples/restart_oce.dat', &
!              &                    cmdstat=ierror, cmdmsg=error_message)
!         if (ierror /= 0) then
!           print *, 'ERROR: could not remove previous soft link restart/checkpoint file.'
!           print *, '       ', trim(error_message)
!           stop
!         endif
!         !..........................................................................................
!         !
!         ierror = 0
!         error_message = ''
!         call execute_command_line('ln -s '//trim(filename)//' examples/restart_oce.dat', &
!              &                     cmdstat=ierror, cmdmsg=error_message)
!         if (ierror /= 0) then
!           print *, 'ERROR: could not soft link restart/checkpoint file.'
!           print *, '       ', trim(error_message)
!           stop
!         endif
!       endif
!     endif

!   end subroutine writeRestart
!   !________________________________________________________________________________________________
!   !
!   subroutine readRestart(currentDate)
!     character(len=max_datetime_str_len), intent(out) :: currentDate
!     character(len=132) :: line

!     integer :: iunit, ierror

!     open (file='examples/restart_oce.dat', status='old', newunit=iunit, iostat=ierror)
!     if (ierror /= 0) then
!       print *, 'ERROR: could not open restart file for reading'
!       print *, '       check isRestart in namelist.'
!       stop
!     else
!       read (unit=iunit, iostat=ierror, iomsg=error_message, fmt='(a)') line
!       if (ierror /= 0) then
!         print *, 'ERROR: could not read restart/checkpoint file.'
!         print *, '       ', trim(error_message)
!         stop
!       endif
!       currentDate=line(10:33)
!       close (unit=iunit)
!     endif

!     print *, 'Read restart/checkpoint file for ', trim(currentDate)

!   end subroutine readRestart
!   !________________________________________________________________________________________________
!   !

!   pure function toLower (str) result (string)

!     character(*), intent(in) :: str
!     character(len(str))      :: string

!     integer :: ic, i

!     character(len=26), parameter :: capitel = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!     character(len=26), parameter :: lower   = 'abcdefghijklmnopqrstuvwxyz'

!     string = str
!     do i = 1, LEN_TRIM(str)
!       ic = INDEX(capitel, str(i:i))
!       if (ic > 0) string(i:i) = lower(ic:ic)
!     end do

!   end function toLower

!   !________________________________________________________________________________________________
!   !

! end program iconoce
