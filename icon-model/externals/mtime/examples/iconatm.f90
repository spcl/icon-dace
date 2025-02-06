!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM iconatm

  USE mtime
  USE mo_event_manager

  IMPLICIT NONE

  TYPE(datetime), POINTER :: experiment_reference_date => NULL()

  TYPE(datetime), POINTER :: experiment_start_date => NULL()
  TYPE(datetime), POINTER :: experiment_end_date => NULL()

  TYPE(datetime), POINTER :: start_date => NULL()
  TYPE(datetime), POINTER :: stop_date => NULL()

  TYPE(datetime), POINTER :: next_checkpoint_date => NULL()
  TYPE(datetime), POINTER :: next_restart_date => NULL()

  TYPE(timedelta), POINTER :: model_time_step => NULL()

  TYPE(datetime), POINTER :: current_date => NULL()

  TYPE(eventgroup), POINTER :: outputEventGroup => NULL()
  TYPE(eventgroup), POINTER :: physicsEventGroup => NULL()

  TYPE(event), POINTER :: checkpointEvent => NULL()
  TYPE(event), POINTER :: restartEvent => NULL()

  TYPE(event), POINTER :: advectionEvent => NULL()

  CHARACTER(len=max_calendar_str_len)  :: calendar_in_use
  CHARACTER(len=max_datetime_str_len)  :: dstring
  CHARACTER(len=max_timedelta_str_len) :: tdstring

  ! namelist variables

  CHARACTER(len=max_calendar_str_len) :: calendar

  CHARACTER(len=max_datetime_str_len) :: experimentReferenceDate = ''
  CHARACTER(len=max_datetime_str_len) :: experimentStartDate
  CHARACTER(len=max_datetime_str_len) :: experimentEndDate

  CHARACTER(len=max_timedelta_str_len) :: modelTimeStep

  CHARACTER(len=max_repetition_str_len) :: advectionTimeInterval
  CHARACTER(len=max_repetition_str_len) :: radiationTimeInterval
  CHARACTER(len=max_repetition_str_len) :: convectionTimeInterval
  CHARACTER(len=max_repetition_str_len) :: cloudTimeInterval
  CHARACTER(len=max_repetition_str_len) :: ssoTimeInterval
  CHARACTER(len=max_repetition_str_len) :: gwdragTimeInterval
  CHARACTER(len=max_repetition_str_len) :: checkpointTimeInterval
  CHARACTER(len=max_repetition_str_len) :: restartTimeInterval
  CHARACTER(len=max_repetition_str_len) :: couplingTimeInterval

  CHARACTER(len=132) :: error_message

  INTEGER :: iunit, icalendar, ierror
  INTEGER :: outputEvents, physicsEvents
  LOGICAL :: lret

  CHARACTER(len=max_groupname_str_len) :: egstring

  TYPE(datetime), POINTER :: checkpointRefDate => NULL()
  TYPE(datetime), POINTER :: checkpointStartDate => NULL()
  TYPE(datetime), POINTER :: checkpointEndDate => NULL()
  TYPE(timedelta), POINTER :: checkpointInterval => NULL()

  TYPE(datetime), POINTER :: restartRefDate => NULL()
  TYPE(datetime), POINTER :: restartStartDate => NULL()
  TYPE(datetime), POINTER :: restartEndDate => NULL()
  TYPE(timedelta), POINTER :: restartInterval => NULL()

  TYPE(datetime), POINTER :: advectionRefDate => NULL()
  TYPE(datetime), POINTER :: advectionStartDate => NULL()
  TYPE(datetime), POINTER :: advectionEndDate => NULL()
  TYPE(timedelta), POINTER :: advectionInterval => NULL()

  !________________________________________________________________________________________________
  !

  NAMELIST /timeControl/ &
       &    calendar, &
       &    experimentReferenceDate, &
       &    experimentStartDate, &
       &    experimentEndDate, &
       &    modelTimeStep, &
       &    advectionTimeInterval, &
       &    radiationTimeInterval, &
       &    convectionTimeInterval, &
       &    cloudTimeInterval, &
       &    ssoTimeInterval, &
       &    gwdragTimeInterval, &
       &    checkpointTimeInterval, &
       &    restartTimeInterval, &
       &    couplingTimeInterval

  OPEN (file='iconatm.nml', newunit=iunit, iostat=ierror)
  IF (ierror /= 0) THEN
    PRINT *, 'ERROR: could not open namelist file.'
    STOP
  ELSE
    READ (unit=iunit, nml=timeControl, iostat=ierror, iomsg=error_message)
    IF (ierror /= 0) THEN
      PRINT *, 'ERROR: could not read namelist file.'
      PRINT *, '       ', TRIM(error_message)
      STOP
    END IF
    CLOSE (unit=iunit)
  END IF

  !________________________________________________________________________________________________
  !

  SELECT CASE (toLower(calendar))
  CASE ('proleptic gregorian')
    icalendar = proleptic_gregorian
  CASE ('365 day year')
    icalendar = year_of_365_days
  CASE ('360 day year')
    icalendar = year_of_360_days
  CASE default
    icalendar = calendar_not_set
    PRINT *, 'ERROR: calendar ', TRIM(calendar), ' not available/unknown.'
    STOP
  END SELECT

  CALL setCalendar(icalendar)
  CALL calendarToString(calendar_in_use)
  PRINT *, 'Calendar: ', TRIM(calendar_in_use)

  PRINT *, ''

  !________________________________________________________________________________________________
  !

  experiment_start_date => newDatetime(experimentStartDate)

  IF (experimentReferenceDate /= '') THEN
    experiment_reference_date => newDatetime(experimentReferenceDate)
  END IF

  IF (ASSOCIATED(experiment_reference_date)) THEN
    CALL initEventManager(experiment_reference_date)
  ELSE
    CALL initEventManager(experiment_start_date)
  END IF

  experiment_reference_date => getModelReferenceDate()

  CALL datetimeToString(experiment_reference_date, dstring)
  PRINT *, 'Experiment reference date: ', dstring

  CALL datetimeToString(experiment_start_date, dstring)
  PRINT *, 'Experiment start date    : ', dstring

  experiment_end_date => newDatetime(experimentEndDate)
  CALL datetimeToString(experiment_end_date, dstring)
  PRINT *, 'Experiment end date      : ', dstring

  PRINT *, ''

  !________________________________________________________________________________________________
  !
  ! event_group_setup: block

  outputEvents = addEventGroup('outputEventGroup')
  outputEventGroup => getEventGroup(outputEvents)
  PRINT *, 'output event group handler: ', outputEvents
  CALL getEventGroupName(outputEventGroup, egstring)
  PRINT *, 'output event group name   : ', TRIM(egstring)
  PRINT *, ''

  physicsEvents = addEventGroup('physicsEventGroup')
  physicsEventGroup => getEventGroup(physicsEvents)
  PRINT *, 'physics event group handler: ', physicsEvents
  CALL getEventGroupName(physicsEventGroup, egstring)
  PRINT *, 'physics event group name   : ', TRIM(egstring)
  PRINT *, ''

  ! end block event_group_setup
  !________________________________________________________________________________________________
  !
  ! checkpoint_restart_time_intervals: block

  checkpointRefDate => experiment_reference_date
  checkpointStartDate => experiment_start_date
  checkpointEndDate => experiment_end_date
  CALL getEventComponents(checkpointTimeInterval, checkpointInterval, checkpointStartDate, checkpointEndDate)
  checkpointEvent => newEvent('checkpoint', checkpointRefDate, checkpointStartDate, checkpointEndDate, &
  & checkpointInterval, errno=ierror)
  IF (ierror /= no_Error) THEN
    PRINT *, 'ERROR: ', ierror
    STOP
  END IF
  lret = addEventToEventGroup(checkpointEvent, outputEventGroup)

  restartRefDate => experiment_reference_date
  restartStartDate => experiment_start_date
  restartEndDate => experiment_end_date
  CALL getEventComponents(restartTimeInterval, restartInterval, restartStartDate, restartEndDate)
  restartEvent => newEvent('restart', restartRefDate, restartStartDate, restartEndDate, restartInterval,  &
  & errno=ierror)
  IF (ierror /= no_Error) THEN
    PRINT *, 'ERROR: ', ierror
    STOP
  END IF
  lret = addEventToEventGroup(restartEvent, outputEventGroup)

  ! end block checkpoint_restart_time_intervals

  CALL printEventGroup(outputEvents)

  !________________________________________________________________________________________________
  !
  ! physics_time_intervals: block

  advectionRefDate => experiment_reference_date
  advectionStartDate => experiment_start_date
  advectionEndDate => experiment_end_date

  CALL getEventComponents(advectionTimeInterval, advectionInterval, advectionStartDate, advectionEndDate)
  advectionEvent => newEvent('advection', advectionRefDate, advectionStartDate, advectionEndDate, advectionInterval, &
  & errno=ierror)
  IF (ierror /= no_Error) THEN
    PRINT *, 'ERROR: ', ierror
    STOP
  END IF
  lret = addEventToEventGroup(advectionEvent, physicsEventGroup)

  CALL printEventGroup(physicsEvents)

  ! end block physics_time_intervals
  !________________________________________________________________________________________________
  !

  PRINT *, ''
  model_time_step => newTimedelta(modelTimeStep)
  CALL timedeltaToString(model_time_step, tdstring)
  PRINT *, 'Dynamics (basic model) time step: ', TRIM(tdstring)
  PRINT *, ''

  !________________________________________________________________________________________________
  !

  start_date => newDatetime(experiment_start_date)
  current_date => newDatetime(start_date)
  stop_date => newDatetime(start_date)
  stop_date = stop_date + getEventInterval(restartEvent)

  !________________________________________________________________________________________________
  !
  ! check_time_interval_consistency: block

  !..............................................................................................
  ! 1. check, if restart is in the experiments time interval
  !
  IF (stop_date > experiment_end_date) THEN
    PRINT *, 'WARNING: run would not create a restart file.'
    PRINT *, '         Reset the stop_date to experiment end_date.'
    stop_date => experiment_end_date
  END IF
  !..............................................................................................
  ! 2. check, if checkpointing is
  !
  next_checkpoint_date => newDatetime(start_date)
  next_checkpoint_date = next_checkpoint_date + getEventInterval(checkpointEvent)
  CALL datetimeToString(next_checkpoint_date, dstring)
  PRINT *, 'First checkpoint date: ', TRIM(dstring)
  !..............................................................................................
  ! 3. check, if restarting is
  !
  next_restart_date => newDatetime(start_date)
  next_restart_date = next_restart_date + getEventInterval(restartEvent)
  CALL datetimeToString(next_restart_date, dstring)
  PRINT *, 'First restart date: ', TRIM(dstring)
  !..............................................................................................

  ! end block check_time_interval_consistency
  !________________________________________________________________________________________________
  !

  CALL datetimeToString(current_date, dstring)
  PRINT *, 'Model date starting the time integration loop: ', TRIM(dstring)

  time_integration: DO
    !............................................................................................
    ! print date and time
    CALL datetimeToString(current_date, dstring)
    PRINT *, 'Model time loop  : ', TRIM(dstring)
    !............................................................................................
    ! need to run advection
    IF (isCurrentEventActive(advectionEvent, current_date)) THEN
      PRINT *, 'Calculate advection: ', TRIM(dstring)
    END IF
    !............................................................................................
    ! initiate restart
    IF ((isCurrentEventActive(restartEvent, current_date) .AND. start_date /= current_date) &
        .OR. current_date == experiment_end_date) THEN
      CALL writeRestart(current_date)
      PRINT *, 'INFO: write restart.'
      EXIT time_integration
    END IF
    !............................................................................................
    ! initiate checkpoint, we do not checkpoint/restart
    IF (isCurrentEventActive(checkpointEvent, current_date) .AND. start_date /= current_date) THEN
      CALL writeRestart(current_date)
      PRINT *, 'INFO: write checkpoint.'
    END IF
    !............................................................................................
    ! calculate next date and time
    current_date = current_date + model_time_step
    !............................................................................................
    ! if new date and time is larger than end of run exit time integration: should never hit
    IF (current_date > stop_date) EXIT time_integration
  END DO time_integration

  CALL datetimeToString(current_date, dstring)
  PRINT *, 'Model date leaving the time integration loop : ', TRIM(dstring)

  !________________________________________________________________________________________________
  !

CONTAINS

  !________________________________________________________________________________________________
  !
  SUBROUTINE writeRestart(currentDate)
    TYPE(datetime), POINTER :: currentDate
    CHARACTER(len=max_datetime_str_len)  :: dstring
    CHARACTER(len=max_datetime_str_len + 12)  :: filename

    INTEGER :: iunit, ierror

    CALL datetimeToString(currentDate, dstring)
    PRINT *, 'Write restart/ceckpoint file for ', TRIM(dstring)

    WRITE (filename, '(a,a,a)') 'restart_', TRIM(dstring), '.dat'

    OPEN (file=TRIM(filename), newunit=iunit, iostat=ierror)
    IF (ierror /= 0) THEN
      PRINT *, 'ERROR: could not open namelist file.'
      STOP
    ELSE
      WRITE (unit=iunit, iostat=ierror, iomsg=error_message, fmt='(a,a)') &
           &  'restart: ', TRIM(dstring)
      IF (ierror /= 0) THEN
        PRINT *, 'ERROR: could not write restart/checkpoint file.'
        PRINT *, '       ', TRIM(error_message)
        STOP
      END IF
      CLOSE (unit=iunit)
    END IF

  END SUBROUTINE writeRestart
  !________________________________________________________________________________________________
  !

  PURE FUNCTION toLower(str) RESULT(string)

    CHARACTER(*), INTENT(in) :: str
    CHARACTER(LEN(str))      :: string

    INTEGER :: ic, i

    CHARACTER(len=26), PARAMETER :: capitel = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(len=26), PARAMETER :: lower = 'abcdefghijklmnopqrstuvwxyz'

    string = str
    DO i = 1, LEN_TRIM(str)
      ic = INDEX(capitel, str(i:i))
      IF (ic > 0) string(i:i) = lower(ic:ic)
    END DO

  END FUNCTION toLower

  !________________________________________________________________________________________________
  !

!   subroutine event_tests

!     type(eventgroup), pointer :: outputEventGroup
!     type(event), pointer :: outputEvent
!     type(event), pointer :: currentEvent
!     type(datetime), pointer :: dtt
!     type(timedelta), pointer :: tdd
!     character(len=max_eventname_str_len) :: currentEventString
!     logical :: lret
!     character(len=max_eventname_str_len) :: aa
!     character(len=max_groupname_str_len) :: bb
!     character(len=max_datetime_str_len)  :: current_date_string_tmp

!     outputEvent => newEvent('output', '2000-01-01T00:00:00', '2010-01-01T00:00:01', '2013-01-01T00:00:02', 'PT06H')
!     lret = addEventToEventGroup(outputEvent, outputEventGroup)

!     dtt => getEventReferenceDateTime(outputEvent)
!     call datetimeToString(dtt, current_date_string_tmp)
!     print *, trim(current_date_string_tmp)

!     dtt => getEventFirstDateTime(outputEvent)
!     call datetimeToString(dtt, current_date_string_tmp)
!     print *, trim(current_date_string_tmp)

!     dtt => getEventLastDateTime(outputEvent)
!     call datetimeToString(dtt, current_date_string_tmp)
!     print *, trim(current_date_string_tmp)

!     tdd => getEventInterval(outputEvent)
!     call timedeltaToString(tdd, current_date_string_tmp)
!     print *, trim(current_date_string_tmp)

!     checkpointEvent => newEvent('checkpoint', '2010-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P05D')
!     lret = addEventToEventGroup(checkpointEvent, outputEventGroup)

!     restartEvent => newEvent('restart', '2000-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P01M')
!     lret = addEventToEventGroup(restartEvent, outputEventGroup)

!     currentEvent => getFirstEventFromEventGroup(outputEventGroup)

!     print *, 'Event list: '
!     do while (associated(currentEvent))
!         call getEventName(currentEvent, currentEventString)
!         print *,'   event: ', trim(currentEventString)
!         currentEvent => getNextEventFromEventGroup(currentEvent)
!     enddo

!     print *,'HELLO' ,getEventId(restartEvent);

!     print *, 'GOOGLE', getEventisFirstInMonth(outputEvent)

!     !type(datetime), pointer :: current_date_test
!     current_date_test => newDatetime('2010-01-02T00:00:00')
!     tmp_date_test_1 => newDatetime('2000-01-01T01:00:00')
!     call getTriggeredPreviousEventAtDateTime(checkpointEvent, tmp_date_test_1)
!     call datetimeToString(tmp_date_test_1, current_date_string)
!     print *, current_date_string

!     call getEventGroupName(outputEventGroup, bb);
!     print *, bb

!     call deallocateEventGroup(outputEventGroup)

!   end subroutine event_tests

END PROGRAM iconatm
