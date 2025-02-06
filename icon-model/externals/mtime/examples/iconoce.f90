!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM iconoce

  USE mtime
  USE mo_event_manager

  IMPLICIT NONE

  TYPE(datetime), POINTER :: experiment_reference_date => NULL()

  TYPE(datetime), POINTER :: checkpoint_reference_date => NULL()

  TYPE(datetime), POINTER :: experiment_start_date => NULL()
  TYPE(datetime), POINTER :: experiment_end_date => NULL()

  TYPE(datetime), POINTER :: start_date => NULL()
  TYPE(datetime), POINTER :: stop_date => NULL()

  TYPE(datetime), POINTER :: next_checkpoint_date => NULL()
  TYPE(datetime), POINTER :: next_restart_date => NULL()

  TYPE(timedelta), POINTER :: model_time_step => NULL()

  TYPE(datetime), POINTER :: current_date => NULL()

  TYPE(eventgroup), POINTER :: outputEventGroup => NULL()

  TYPE(event), POINTER :: checkpointEvent => NULL()
  TYPE(event), POINTER :: restartEvent => NULL()

  CHARACTER(len=max_calendar_str_len)  :: calendar_in_use
  CHARACTER(len=max_datetime_str_len)  :: dstring
  CHARACTER(len=max_timedelta_str_len) :: tdstring

  CHARACTER(len=max_mtime_error_str_len) :: errstring
  ! namelist variables

  CHARACTER(len=max_calendar_str_len) :: calendar

  CHARACTER(len=max_datetime_str_len) :: experimentReferenceDate = ''
  CHARACTER(len=max_datetime_str_len) :: experimentStartDate = ''
  CHARACTER(len=max_datetime_str_len) :: experimentEndDate

  CHARACTER(len=max_datetime_str_len) :: startDate

  CHARACTER(len=max_timedelta_str_len) :: modelTimeStep

  CHARACTER(len=max_repetition_str_len) :: checkpointTimeInterval
  CHARACTER(len=max_repetition_str_len) :: restartTimeInterval

  CHARACTER(len=max_repetition_str_len) :: couplingTimeInterval

  CHARACTER(len=132) :: error_message

  INTEGER :: iunit, icalendar, ierror
  INTEGER :: outputEvents
  LOGICAL :: lret, isRestart, isRestartTimeRelative

  CHARACTER(len=max_groupname_str_len) :: egstring

  TYPE(datetime), POINTER :: checkpointRefDate => NULL()
  TYPE(datetime), POINTER :: checkpointStartDate => NULL()
  TYPE(datetime), POINTER :: checkpointEndDate => NULL()
  TYPE(timedelta), POINTER :: checkpointInterval => NULL()

  TYPE(datetime), POINTER :: restartRefDate => NULL()
  TYPE(datetime), POINTER :: restartStartDate => NULL()
  TYPE(datetime), POINTER :: restartEndDate => NULL()
  TYPE(timedelta), POINTER :: restartInterval => NULL()

  !________________________________________________________________________________________________
  !

  NAMELIST /timeControl/ &
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

  OPEN (file='iconoce.nml', newunit=iunit, iostat=ierror)
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

  IF (experimentReferenceDate /= '') THEN
    experiment_reference_date => newDatetime(experimentReferenceDate)
  END IF

  IF (isRestart) THEN
    CALL readRestart(startDate)
    start_date => newDatetime(startDate)
  ELSE
    start_date => newDatetime(experimentStartDate)
  END IF

  experiment_start_date => newDatetime(experimentStartDate)

  IF (isRestartTimeRelative) THEN
    checkpoint_reference_date => newDatetime(start_date)
  ELSE
    checkpoint_reference_date => newDatetime(experiment_reference_date)
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

  ! end block event_group_setup
  !________________________________________________________________________________________________
  !
  ! checkpoint_restart_time_intervals: block

  checkpointRefDate => checkpoint_reference_date
  checkpointStartDate => experiment_start_date
  checkpointEndDate => experiment_end_date
  CALL getEventComponents(checkpointTimeInterval, &
       &                  checkpointInterval, checkpointStartDate, checkpointEndDate)
  checkpointEvent => newEvent('checkpoint', checkpointRefDate, &
       &                      checkpointStartDate, checkpointEndDate, checkpointInterval, &
       &                      errno=ierror)
  IF (ierror /= no_Error) THEN
    CALL mtime_strerror(ierror, errstring)
    PRINT *, 'ERROR: ', TRIM(errstring)
    STOP
  END IF
  lret = addEventToEventGroup(checkpointEvent, outputEventGroup)

  restartRefDate => checkpoint_reference_date
  restartStartDate => experiment_start_date
  restartEndDate => experiment_end_date
  CALL getEventComponents(restartTimeInterval, &
       &                  restartInterval, restartStartDate, restartEndDate)
  restartEvent => newEvent('restart', restartRefDate, &
       &                   restartStartDate, restartEndDate, restartInterval, &
       &                   errno=ierror)
  IF (ierror /= no_Error) THEN
    CALL mtime_strerror(ierror, errstring)
    PRINT *, 'ERROR: ', TRIM(errstring)
    STOP
  END IF
  lret = addEventToEventGroup(restartEvent, outputEventGroup)

  ! end block checkpoint_restart_time_intervals

  CALL printEventGroup(outputEvents)

  !________________________________________________________________________________________________
  !

  PRINT *, ''
  model_time_step => newTimedelta(modelTimeStep)
  CALL timedeltaToString(model_time_step, tdstring)
  PRINT *, 'Dynamics (basic model) time step: ', TRIM(tdstring)
  PRINT *, ''

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
  !end block check_time_interval_consistency
  !________________________________________________________________________________________________
  !

  CALL datetimeToString(current_date, dstring)
  PRINT *, 'Model date starting the time integration loop: ', TRIM(dstring)

  time_integration: DO
    !............................................................................................
    ! print date and time
    CALL datetimeToString(current_date, dstring)
!    print *, 'Model time loop  : ', trim(dstring)
    !............................................................................................
    ! initiate restart
    IF ((isCurrentEventActive(restartEvent, current_date) .AND. start_date /= current_date) &
        .OR. current_date == experiment_end_date) THEN
      PRINT *, 'INFO: write restart.'
      CALL writeRestart(current_date)
      EXIT time_integration
    END IF
    !............................................................................................
    ! initiate checkpoint, we do not checkpoint/restart
    IF (isCurrentEventActive(checkpointEvent, current_date) .AND. start_date /= current_date) THEN
      PRINT *, 'INFO: write checkpoint.'
      CALL writeRestart(current_date)
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
    CHARACTER(len=max_datetime_str_len + 16)  :: filename

    INTEGER :: iunit, ierror

    CALL datetimeToString(currentDate, dstring)
    PRINT *, '      Writing restart/checkpoint file for ', TRIM(dstring)

    WRITE (filename, '(a,a,a)') 'restart_oce_', TRIM(dstring), '.dat'

    OPEN (file=TRIM(filename), newunit=iunit, iostat=ierror)
    IF (ierror /= 0) THEN
      PRINT *, 'ERROR: could not open restart file for writing.'
      STOP
    ELSE
      WRITE (unit=iunit, iostat=ierror, iomsg=error_message, fmt='(a,a)') &
           & 'restart: ', TRIM(dstring)
      IF (ierror /= 0) THEN
        PRINT *, 'ERROR: could not write restart/checkpoint file.'
        PRINT *, '       ', TRIM(error_message)
        STOP
      ELSE
        CLOSE (unit=iunit)
        !..........................................................................................
        !
        ierror = 0
        error_message = ''
        CALL execute_command_line('rm -f restart_oce.dat', &
             &                    cmdstat=ierror, cmdmsg=error_message)
        IF (ierror /= 0) THEN
          PRINT *, 'ERROR: could not remove previous soft link restart/checkpoint file.'
          PRINT *, '       ', TRIM(error_message)
          STOP
        END IF
        !..........................................................................................
        !
        ierror = 0
        error_message = ''
        CALL execute_command_line('ln -s '//TRIM(filename)//' restart_oce.dat', &
             &                     cmdstat=ierror, cmdmsg=error_message)
        IF (ierror /= 0) THEN
          PRINT *, 'ERROR: could not soft link restart/checkpoint file.'
          PRINT *, '       ', TRIM(error_message)
          STOP
        END IF
      END IF
    END IF

  END SUBROUTINE writeRestart
  !________________________________________________________________________________________________
  !
  SUBROUTINE readRestart(currentDate)
    CHARACTER(len=max_datetime_str_len), INTENT(out) :: currentDate
    CHARACTER(len=132) :: line

    INTEGER :: iunit, ierror

    OPEN (file='restart_oce.dat', status='old', newunit=iunit, iostat=ierror)
    IF (ierror /= 0) THEN
      PRINT *, 'ERROR: could not open restart file for reading'
      PRINT *, '       check isRestart in namelist.'
      STOP
    ELSE
      READ (unit=iunit, iostat=ierror, iomsg=error_message, fmt='(a)') line
      IF (ierror /= 0) THEN
        PRINT *, 'ERROR: could not read restart/checkpoint file.'
        PRINT *, '       ', TRIM(error_message)
        STOP
      END IF
      currentDate = line(10:33)
      CLOSE (unit=iunit)
    END IF

    PRINT *, 'Read restart/checkpoint file for ', TRIM(currentDate)

  END SUBROUTINE readRestart
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

END PROGRAM iconoce
