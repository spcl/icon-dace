!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program iconoce

  use mtime
  use mo_event_manager

  implicit none

  type(datetime), pointer :: experiment_reference_date => null()

  type(datetime), pointer :: checkpoint_reference_date => null()  

  type(datetime), pointer :: experiment_start_date => null()
  type(datetime), pointer :: experiment_end_date => null()

  type(datetime), pointer :: start_date => null() 
  type(datetime), pointer :: stop_date => null() 

  type(datetime), pointer :: next_checkpoint_date => null() 
  type(datetime), pointer :: next_restart_date => null() 

  type(timedelta), pointer :: dynamic_time_step => null()

  type(timedelta), pointer :: checkpoint_time_step => null()
  type(timedelta), pointer :: restart_time_step => null()

  type(timedelta), pointer :: coupling_time_step => null()
  type(timedelta), pointer :: model_time_step => null()

  type(datetime), pointer :: current_date => null() 

  type(eventgroup), pointer :: outputEventGroup => null()

  type(event), pointer :: checkpointEvent => null()
  type(event), pointer :: restartEvent => null()


  character(len=max_calendar_str_len)  :: calendar_in_use
  character(len=max_datetime_str_len)  :: dstring
  character(len=max_timedelta_str_len) :: tdstring

  character(len=max_mtime_error_str_len  ) :: errstring
  ! namelist variables

  character(len=max_calendar_str_len) :: calendar

  character(len=max_datetime_str_len) :: experimentReferenceDate = ''   
  character(len=max_datetime_str_len) :: experimentStartDate = ''  
  character(len=max_datetime_str_len) :: experimentEndDate

  character(len=max_datetime_str_len) :: startDate     
                                                                 
  character(len=max_timedelta_str_len) :: modelTimeStep         
                                                                 
  character(len=max_repetition_str_len) :: checkpointTimeInterval
  character(len=max_repetition_str_len) :: restartTimeInterval   

  character(len=max_repetition_str_len) :: couplingTimeInterval

  character(len=132) :: error_message

  integer :: iunit, icalendar, ierror
  integer :: outputEvents
  logical :: lret, isRestart, isRestartTimeRelative

  character(len=max_groupname_str_len) :: egstring

  type(datetime),  pointer :: checkpointRefDate => null()
  type(datetime),  pointer :: checkpointStartDate => null()
  type(datetime),  pointer :: checkpointEndDate => null()
  type(timedelta), pointer :: checkpointInterval => null()
  
  type(datetime),  pointer :: restartRefDate => null()
  type(datetime),  pointer :: restartStartDate => null()
  type(datetime),  pointer :: restartEndDate => null()
  type(timedelta), pointer :: restartInterval => null()
  
  !________________________________________________________________________________________________
  !

  namelist /timeControl/ &
       &    calendar, & 
       &    experimentReferenceDate, &
       &    experimentStartDate, &   
       &    experimentEndDate, &        
       &    modelTimeStep, &            
       &    checkpointTimeInterval, &   
       &    restartTimeInterval, &      
       &    couplingTimeInterval, &
       &    isRestart, &
       &    isRestartTimeRelative

  open (file='examples/iconoce.nml', newunit=iunit, iostat=ierror)
  if (ierror /= 0) then
    print *, 'ERROR: could not open namelist file.' 
    stop 
  else
    read (unit=iunit, nml=timeControl, iostat=ierror, iomsg=error_message)
    if (ierror /= 0) then
      print *, 'ERROR: could not read namelist file.'
      print *, '       ', trim(error_message)  
      stop 
    endif
    close (unit=iunit)
  endif

  !________________________________________________________________________________________________
  !

  select case (toLower(calendar))
  case ('proleptic gregorian')
    icalendar  = proleptic_gregorian
  case ('365 day year')  
    icalendar = year_of_365_days
  case ('360 day year')  
    icalendar = year_of_360_days
  case default
    icalendar = calendar_not_set
    print *, 'ERROR: calendar ', trim(calendar), ' not available/unknown.' 
    stop 
  end select

  call setCalendar(icalendar)
  call calendarToString(calendar_in_use)
  print *, 'Calendar: ', trim(calendar_in_use)
  
  print *, ''

  !________________________________________________________________________________________________
  !

  if (experimentReferenceDate /= '') then
    experiment_reference_date => newDatetime(experimentReferenceDate)
  endif

  if (isRestart) then
    call readRestart(startDate)
    start_date => newDatetime(startDate)
  else
    start_date => newDatetime(experimentStartDate)    
  endif

  experiment_start_date => newDatetime(experimentStartDate)

  if (isRestartTimeRelative) then
    checkpoint_reference_date => newDatetime(start_date)    
  else
    checkpoint_reference_date => newDatetime(experiment_reference_date)    
  endif
  
  if (associated(experiment_reference_date)) then
    call initEventManager(experiment_reference_date)
  else
    call initEventManager(experiment_start_date)
  endif

  experiment_reference_date => getModelReferenceDate()

  call datetimeToString(experiment_reference_date, dstring)
  print *, 'Experiment reference date: ', dstring

  call datetimeToString(experiment_start_date, dstring)
  print *, 'Experiment start date    : ', dstring

  experiment_end_date => newDatetime(experimentEndDate)
  call datetimeToString(experiment_end_date, dstring)
  print *, 'Experiment end date      : ', dstring
  
  print *, ''

  !________________________________________________________________________________________________
  !
  ! event_group_setup: block

  outputEvents =  addEventGroup('outputEventGroup')
  outputEventGroup => getEventGroup(outputEvents)
  print *, 'output event group handler: ', outputEvents
  call getEventGroupName(outputEventGroup, egstring)
  print *, 'output event group name   : ', trim(egstring)    
  print *, ''

  ! end block event_group_setup
  !________________________________________________________________________________________________
  !
  ! checkpoint_restart_time_intervals: block

  checkpointRefDate   => checkpoint_reference_date
  checkpointStartDate => experiment_start_date
  checkpointEndDate   => experiment_end_date
  call getEventComponents(checkpointTimeInterval, checkpointRefDate, &
       &                  checkpointInterval, checkpointStartDate, checkpointEndDate)
  checkpointEvent => newEvent('checkpoint', checkpointRefDate, &
       &                      checkpointStartDate, checkpointEndDate, checkpointInterval, &
       &                      errno=ierror)
  if (ierror /= no_Error) then
    CALL mtime_strerror(ierror, errstring)
    print *, 'ERROR: ', trim(errstring)
    stop
  endif
  lret = addEventToEventGroup(checkpointEvent, outputEventGroup)
  
  restartRefDate   => checkpoint_reference_date
  restartStartDate => experiment_start_date
  restartEndDate   => experiment_end_date
  call getEventComponents(restartTimeInterval, restartRefDate, &
       &                  restartInterval, restartStartDate, restartEndDate)
  restartEvent => newEvent('restart', restartRefDate, &
       &                   restartStartDate, restartEndDate, restartInterval, &
       &                   errno=ierror)
  if (ierror /= no_Error) then
    CALL mtime_strerror(ierror, errstring)
    print *, 'ERROR: ', trim(errstring)
    stop
  endif
  lret = addEventToEventGroup(restartEvent, outputEventGroup)

  ! end block checkpoint_restart_time_intervals
  
  call printEventGroup(outputEvents)
  
  !________________________________________________________________________________________________
  !

  print *, ''
  model_time_step => newTimedelta(modelTimeStep)
  call timedeltaToString(model_time_step, tdstring)
  print *, 'Dynamics (basic model) time step: ', trim(tdstring)
  print *, ''

  !________________________________________________________________________________________________
  ! 

  current_date => newDatetime(start_date)
  stop_date => newDatetime(start_date) 
  stop_date = stop_date + getEventInterval(restartEvent)

  !________________________________________________________________________________________________
  !
  ! check_time_interval_consistency: block
  !..............................................................................................
  ! 1. check, if restart is in the experiments time interval
  !
  if (stop_date > experiment_end_date) then
    print *, 'WARNING: run would not create a restart file.'
    print *, '         Reset the stop_date to experiment end_date.'
    stop_date => experiment_end_date
  endif
  !..............................................................................................
  ! 2. check, if checkpointing is
  !
  next_checkpoint_date => newDatetime(start_date)
  next_checkpoint_date = next_checkpoint_date + getEventInterval(checkpointEvent)
  call datetimeToString(next_checkpoint_date, dstring)
  print *, 'First checkpoint date: ', trim(dstring)
  !..............................................................................................
  ! 3. check, if restarting is
  !
  next_restart_date => newDatetime(start_date)
  next_restart_date = next_restart_date + getEventInterval(restartEvent)
  call datetimeToString(next_restart_date, dstring)
  print *, 'First restart date: ', trim(dstring)
  !..............................................................................................
  !end block check_time_interval_consistency
  !________________________________________________________________________________________________
  !
  
  call datetimeToString(current_date, dstring)
  print *, 'Model date starting the time integration loop: ', trim(dstring)
  
  time_integration: do 
    !............................................................................................
    ! print date and time
    call datetimeToString(current_date, dstring)
!    print *, 'Model time loop  : ', trim(dstring)
    !............................................................................................
    ! initiate restart
    if ((isCurrentEventActive(restartEvent, current_date) .and. start_date /= current_date) &
         .or. current_date == experiment_end_date) then
      print *, 'INFO: write restart.'
      call writeRestart(current_date)
      exit time_integration
    endif
    !............................................................................................
    ! initiate checkpoint, we do not checkpoint/restart
    if (isCurrentEventActive(checkpointEvent, current_date) .and. start_date /= current_date) then
      print *, 'INFO: write checkpoint.'
      call writeRestart(current_date)
    endif
    !............................................................................................
    ! calculate next date and time
    current_date = current_date + model_time_step
    !............................................................................................
    ! if new date and time is larger than end of run exit time integration: should never hit
    if (current_date > stop_date) exit time_integration
  enddo time_integration
  
  call datetimeToString(current_date, dstring)
  print *, 'Model date leaving the time integration loop : ', trim(dstring)
    
  !________________________________________________________________________________________________
  !

contains

  !________________________________________________________________________________________________
  !
  subroutine writeRestart(currentDate)
    type(datetime), pointer :: currentDate
    character(len=max_datetime_str_len)  :: dstring
    character(len=max_datetime_str_len+16)  :: filename

    integer :: iunit, ierror

    call datetimeToString(currentDate, dstring)
    print *, '      Writing restart/checkpoint file for ', trim(dstring)

    write (filename,'(a,a,a)') 'restart_oce_', trim(dstring), '.dat'

    open (file='examples/'//trim(filename), newunit=iunit, iostat=ierror)
    if (ierror /= 0) then
      print *, 'ERROR: could not open restart file for writing.' 
      stop 
    else
      write (unit=iunit, iostat=ierror, iomsg=error_message, fmt='(a,a)') &
           & 'restart: ', trim(dstring)
      if (ierror /= 0) then
        print *, 'ERROR: could not write restart/checkpoint file.'
        print *, '       ', trim(error_message)  
        stop 
      else
        close (unit=iunit)
        !..........................................................................................
        !
        ierror = 0
        error_message = ''
        call execute_command_line('rm -f examples/restart_oce.dat', &
             &                    cmdstat=ierror, cmdmsg=error_message)
        if (ierror /= 0) then
          print *, 'ERROR: could not remove previous soft link restart/checkpoint file.'
          print *, '       ', trim(error_message)  
          stop
        endif
        !..........................................................................................
        !
        ierror = 0
        error_message = ''
        call execute_command_line('ln -s '//trim(filename)//' examples/restart_oce.dat', &
             &                     cmdstat=ierror, cmdmsg=error_message)
        if (ierror /= 0) then
          print *, 'ERROR: could not soft link restart/checkpoint file.'
          print *, '       ', trim(error_message)  
          stop
        endif
      endif
    endif

  end subroutine writeRestart
  !________________________________________________________________________________________________
  !
  subroutine readRestart(currentDate)
    character(len=max_datetime_str_len), intent(out) :: currentDate
    character(len=132) :: line
    
    integer :: iunit, ierror
    
    open (file='examples/restart_oce.dat', status='old', newunit=iunit, iostat=ierror)
    if (ierror /= 0) then
      print *, 'ERROR: could not open restart file for reading'
      print *, '       check isRestart in namelist.' 
      stop 
    else
      read (unit=iunit, iostat=ierror, iomsg=error_message, fmt='(a)') line
      if (ierror /= 0) then
        print *, 'ERROR: could not read restart/checkpoint file.'
        print *, '       ', trim(error_message)  
        stop 
      endif
      currentDate=line(10:33)
      close (unit=iunit)
    endif

    print *, 'Read restart/checkpoint file for ', trim(currentDate)

  end subroutine readRestart
  !________________________________________________________________________________________________
  !
  
  pure function toLower (str) result (string)
    
    character(*), intent(in) :: str
    character(len(str))      :: string
    
    integer :: ic, i
    
    character(len=26), parameter :: capitel = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=26), parameter :: lower   = 'abcdefghijklmnopqrstuvwxyz'
    
    string = str
    do i = 1, LEN_TRIM(str)
      ic = INDEX(capitel, str(i:i))
      if (ic > 0) string(i:i) = lower(ic:ic)
    end do
    
  end function toLower

  !________________________________________________________________________________________________
  !

end program iconoce
