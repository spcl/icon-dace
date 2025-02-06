! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Creation and destruction of mtime events for the wave model

MODULE mo_wave_events

  USE mo_exception,                ONLY: message, message_text, finish
  USE mtime,                       ONLY: datetime, timedelta, newTimedelta,              &
    &                                    deallocateTimedelta, no_Error, mtime_strerror,  &
    &                                    datetimeToString, timedeltaToString,            &
    &                                    MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN,    &
    &                                    MAX_MTIME_ERROR_STR_LEN, MAX_EVENTNAME_STR_LEN, &
    &                                    event, eventGroup, newEvent, addEventToEventGroup
  USE mo_event_manager,            ONLY: initEventManager, addEventGroup, getEventGroup, printEventGroup
  USE mo_time_config,              ONLY: t_time_config

  IMPLICIT NONE

  PRIVATE

  ! subroutine
  PUBLIC :: create_wave_events

  ! events
  PUBLIC :: waveDummyEvent
  PUBLIC :: waveCheckpointEvent
  PUBLIC :: waveRestartEvent
  PUBLIC :: waveEventGroup

  TYPE(event),      POINTER :: waveDummyEvent      => NULL()
  TYPE(event),      POINTER :: waveCheckpointEvent => NULL()
  TYPE(event),      POINTER :: waveRestartEvent    => NULL()
  TYPE(eventGroup), POINTER :: waveEventGroup      => NULL()

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_events'

CONTAINS

  !>
  !! Create mtime events for wave model
  !!
  !! This routine creates mtime events for the wave model and
  !! puts them into suitable event groups.
  !!
  SUBROUTINE create_wave_events (time_config)

    TYPE(t_time_config), INTENT(IN) :: time_config  !< information for time control

    TYPE(datetime), POINTER         :: eventStartDate    => NULL(), &
      &                                eventEndDate      => NULL(), &
      &                                eventRefDate      => NULL(), &
      &                                checkpointRefDate => NULL(), &
      &                                restartRefDate    => NULL()
    TYPE(timedelta), POINTER        :: eventInterval     => NULL()
    INTEGER                         :: waveEventsGroupID   !< might be changed to PUBLIC module variable,
                                                           !< if needed outside of this module
    INTEGER                             :: ierr
    LOGICAL                             :: lret
    CHARACTER(len=MAX_EVENTNAME_STR_LEN):: eventName

    !
    ! Create an event manager, ie. a collection of different events
    !
    CALL initEventManager(time_config%tc_exp_refdate)

    ! create an event group for the wave model
    waveEventsGroupID = addEventGroup('waveEventGroup')
    waveEventGroup    => getEventGroup(waveEventsGroupID)

    !
    ! create an example dummy event
    !
    eventRefDate   => time_config%tc_exp_startdate
    eventStartDate => time_config%tc_exp_startdate
    eventEndDate   => time_config%tc_exp_stopdate
    eventInterval  => newTimedelta("PT01h")
    eventName      =  'waveDummy'

    waveDummyEvent => newEvent(TRIM(eventName), eventRefDate, eventStartDate, &
      &                           eventEndDate, eventInterval, errno=ierr)
    IF (ierr /= no_Error) THEN
      CALL event_write_dbg_info(eventName, eventRefDate, eventStartDate, eventEndDate, eventInterval)
    ENDIF

    !
    lret = addEventToEventGroup(waveDummyEvent, waveEventGroup)
    CALL deallocateTimedelta(eventInterval)
    NULLIFY(eventRefDate,eventStartDate,eventEndDate)


    ! for debugging purposes the reference (anchor) date for checkpoint
    ! and restart may be switched to be relative to current jobs start
    ! date instead of the experiments start date.
    IF (time_config%is_relative_time) THEN
      checkpointRefDate => time_config%tc_startdate
      restartRefDate    => time_config%tc_startdate
    ELSE
      checkpointRefDate => time_config%tc_exp_startdate
      restartRefDate    => time_config%tc_exp_startdate
    ENDIF

    !
    ! create checkpoint event
    !
    eventStartDate => time_config%tc_exp_startdate
    eventEndDate   => time_config%tc_exp_stopdate
    eventInterval  => time_config%tc_dt_checkpoint
    eventName      =  'waveCheckpoint'

    waveCheckpointEvent => newEvent(TRIM(eventName), checkpointRefDate, eventStartDate, &
      &                         eventEndDate, eventInterval, errno=ierr)
    IF (ierr /= no_Error) THEN
      CALL event_write_dbg_info(eventName, checkpointRefDate, eventStartDate, eventEndDate, eventInterval)
    ENDIF
    !
    lret = addEventToEventGroup(waveCheckpointEvent, waveEventGroup)
    NULLIFY(checkpointRefDate,eventStartDate,eventEndDate,eventInterval)

    !
    ! create restart event
    !
    eventStartDate => time_config%tc_exp_startdate
    eventEndDate   => time_config%tc_exp_stopdate
    eventInterval  => time_config%tc_dt_restart
    eventName      =  'waveRestart'

    waveRestartEvent => newEvent(TRIM(eventName), restartRefDate, eventStartDate, &
      &                      eventEndDate, eventInterval, errno=ierr)
    IF (ierr /= no_Error) THEN
      CALL event_write_dbg_info(eventName, restartRefDate, eventStartDate, eventEndDate, eventInterval)
    ENDIF
    !
    lret = addEventToEventGroup(waveRestartEvent, waveEventGroup)
    NULLIFY(restartRefDate,eventStartDate,eventEndDate,eventInterval)


    CALL printEventGroup(waveEventsGroupID)

  END SUBROUTINE create_wave_events


  ! Write mtime event debugging information
  !
  ! Writes mtime event debugging information to stdout.
  !
  SUBROUTINE event_write_dbg_info (eventName, eventRefDate, eventStartDate, eventEndDate, eventInterval)

    CHARACTER(len=*), PARAMETER :: routine = modname//':event_write_dbg_info'

    CHARACTER(len=MAX_EVENTNAME_STR_LEN), INTENT(IN):: eventName
    TYPE(datetime),  INTENT(IN):: eventRefDate
    TYPE(datetime),  INTENT(IN):: eventStartDate
    TYPE(datetime),  INTENT(IN):: eventEndDate
    TYPE(timedelta), INTENT(IN):: eventInterval

    ! local
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN)   :: td_string
    CHARACTER(len=MAX_DATETIME_STR_LEN)    :: dt_string
    CHARACTER(len=MAX_MTIME_ERROR_STR_LEN) :: errstring
    INTEGER:: ierr

    ! give an elaborate error message:
    CALL datetimeToString(eventRefDate,    dt_string)
    WRITE(message_text,'(a,a)') "event reference date: ",  dt_string
    CALL message(routine, message_text)
    !
    CALL datetimeToString(eventStartDate,  dt_string)
    WRITE(message_text,'(a,a)') "event start date: ",  dt_string
    CALL message(routine, message_text)
    !
    CALL datetimeToString(eventEndDate,    dt_string)
    WRITE(message_text,'(a,a)') "event end date: ",  dt_string
    CALL message(routine, message_text)
    !
    CALL timedeltaToString(eventInterval,  td_string)
    WRITE(message_text,'(a,a)') "event interval: ",  td_string
    CALL message(routine, message_text)
    !
    CALL mtime_strerror(ierr, errstring)
    WRITE(message_text,'(a,a)') TRIM(eventName),": ",errstring
    CALL finish(routine, message_text)

  END SUBROUTINE event_write_dbg_info

END MODULE mo_wave_events

