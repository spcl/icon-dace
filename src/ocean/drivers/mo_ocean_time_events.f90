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

! Contains the ocean time keeping and events

MODULE mo_ocean_time_events
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mtime,                     ONLY: datetime, datetimeToString, deallocateDatetime,              &
       &                               timedelta, newTimedelta, deallocateTimedelta,                &
       &                               MAX_DATETIME_STR_LEN, newDatetime,                           &
       &                               MAX_MTIME_ERROR_STR_LEN, no_error, mtime_strerror,           &
       &                               OPERATOR(-), OPERATOR(+), OPERATOR(>), OPERATOR(*),          &
       &                               ASSIGNMENT(=), OPERATOR(==), OPERATOR(>=), OPERATOR(/=),     &
       &                               event, eventGroup, newEvent,                                 &
       &                               addEventToEventGroup, isCurrentEventActive
  USE mo_event_manager,          ONLY: initEventManager, addEventGroup, getEventGroup, printEventGroup  
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_master_config,          ONLY: isRestart
  USE mo_time_config,            ONLY: time_config
  USE mo_run_config,             ONLY: nsteps, config_dtime => dtime
  USE mo_io_config,              ONLY: write_last_restart

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_ocean_time_events
  PUBLIC :: ocean_time_nextStep
  PUBLIC :: isCheckpoint
  PUBLIC :: isEndOfThisRun
  PUBLIC :: getCurrentDate_to_String
  PUBLIC :: isStartdate
  PUBLIC :: newNullDatetime
  
  PUBLIC :: get_OceanCurrentTime_Pointer 


  CHARACTER(LEN=20)  :: str_module = 'mo_ocean_time_events'  ! Output of module for 1 line debug
  !-------------------------------------------------------------------------

  TYPE(eventGroup), POINTER           :: checkpointEventGroup => NULL()

  TYPE(timedelta), POINTER            :: ocean_time_step => NULL()

  TYPE(datetime), POINTER             :: ocean_current_time     => NULL()
  TYPE(datetime), POINTER             :: ocean_previous_time    => NULL()
  TYPE(datetime), POINTER             :: eventRefDate      => NULL(), &
        &                                eventStartDate    => NULL(), &
        &                                eventEndDate      => NULL()
  TYPE(datetime), POINTER             :: checkpointRefDate => NULL(), &
        &                                restartRefDate    => NULL()

  TYPE(datetime), POINTER             :: return_current_time => NULL()

  TYPE(event), POINTER                :: checkpointEvent => NULL()
  TYPE(event), POINTER                :: restartEvent    => NULL()
  
  INTEGER                             :: checkpointEvents
    
  REAL(wp):: ocean_dtime           !< [s] length of a time step

CONTAINS

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  ! set events, group and the events
  SUBROUTINE init_ocean_time_events()

    CHARACTER(len=*), PARAMETER :: method_name = "init_ocean_time_events"
    
    INTEGER :: ierr
    LOGICAL :: return_status
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: dstring
    CHARACTER(len=MAX_MTIME_ERROR_STR_LEN) :: errstring
    TYPE(timedelta), POINTER               :: eventInterval   => NULL()

    ocean_current_time  => time_config%tc_current_date
    ocean_dtime         = config_dtime
    ocean_previous_time => newDatetime(ocean_current_time) ! - ocean_dtime
    eventRefDate        => time_config%tc_exp_refdate
    eventStartDate      => time_config%tc_exp_startdate
    eventEndDate        => time_config%tc_exp_stopdate
    ocean_time_step     => time_config%tc_dt_model

    ! the time varibales for returning 
    return_current_time => newNullDatetime()


    ! for debugging purposes the referenece (anchor) date for checkpoint
    ! and restart may be switched to be relative to current jobs start
    ! date instead of the experiments start date.
    
    IF (time_config%is_relative_time) THEN
      checkpointRefDate => time_config%tc_startdate
      restartRefDate    => time_config%tc_startdate
    ELSE
      checkpointRefDate => time_config%tc_exp_startdate
      restartRefDate    => time_config%tc_exp_startdate
    ENDIF
    
    ! create an event manager, ie. a collection of different events
    CALL initEventManager(time_config%tc_exp_refdate)

    ! --- create an event group for checkpointing and restart
    checkpointEvents =  addEventGroup('checkpointEventGroup')
    checkpointEventGroup => getEventGroup(checkpointEvents)
    
    ! --- --- create checkpointing event
    eventInterval  => time_config%tc_dt_checkpoint
    checkpointEvent => newEvent('checkpoint', checkpointRefDate, eventStartDate, eventEndDate, eventInterval, errno=ierr)
    IF (ierr /= no_Error) THEN
       CALL mtime_strerror(ierr, errstring)
       CALL finish(method_name, errstring)
    ENDIF
    return_status = addEventToEventGroup(checkpointEvent, checkpointEventGroup)

    ! --- --- create restart event, ie. checkpoint + model stop
    eventInterval  => time_config%tc_dt_restart
    restartEvent => newEvent('restart', restartRefDate, eventStartDate, eventEndDate, eventInterval, errno=ierr)
    IF (ierr /= no_Error) THEN
       CALL mtime_strerror(ierr, errstring)
       CALL finish('perform_ho_timeloop', errstring)
    ENDIF
    return_status = addEventToEventGroup(restartEvent, checkpointEventGroup)

    ! print some info
    CALL message('','')
    CALL datetimeToString(ocean_current_time, dstring)
    WRITE(message_text,'(a,a)') 'Start date of this run: ', dstring
    CALL message('',message_text)
    CALL datetimeToString(time_config%tc_stopdate, dstring)
    WRITE(message_text,'(a,a)') 'Stop date of this run:  ', dstring
    CALL message('',message_text)
    CALL message('','')

   CALL printEventGroup(checkpointEvents)
  
  END SUBROUTINE init_ocean_time_events
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  FUNCTION newNullDatetime()
    TYPE(datetime), POINTER:: newNullDatetime 
    
    newNullDatetime => newDatetime('0001-01-01T00:00:00')
  END FUNCTION newNullDatetime
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  FUNCTION ocean_time_nextStep() 
    TYPE(datetime), POINTER:: ocean_time_nextStep 
     
    ocean_previous_time = ocean_current_time
    ocean_current_time = ocean_current_time + ocean_time_step
    return_current_time = ocean_current_time
    ocean_time_nextStep => return_current_time

  END FUNCTION ocean_time_nextStep
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  LOGICAL FUNCTION isStartdate()
    isStartdate    = (time_config%tc_startdate == ocean_current_time)
  END FUNCTION isStartdate
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  LOGICAL FUNCTION isCheckpoint()

    LOGICAL :: isRestart,isThisCheckpoint,doWriteRestart
    
    isCheckpoint = .false.

    isRestart      = isCurrentEventActive(restartEvent, ocean_current_time)
    isThisCheckpoint = isCurrentEventActive(checkpointEvent, ocean_current_time)
    doWriteRestart = time_config%tc_write_restart
        
    isCheckpoint = (isRestart .OR. isThisCheckpoint) .AND. .NOT. isStartdate() .AND. doWriteRestart

    isCheckpoint = isCheckpoint .OR. (write_last_restart .AND. isEndOfThisRun())

!        isExpStopdate  = (time_config%tc_exp_stopdate == ocean_current_time)
!         IF ( &
!              !  if normal checkpoint or restart cycle has been reached, i.e. checkpoint+model stop
!              &         (isRestart .OR. isThisCheckpoint)                     &
!              &  .AND.                                                        &
!              !  and the current date differs from the start date
!              &        .NOT. isStartdate                                    &
!              &  .AND.                                                        &
!              !  and end of run has not been reached or restart writing has been disabled
!              &        (.NOT.isExpStopdate .OR.doWriteRestart)          &
!              & ) THEN
!           isCheckpoint = .TRUE.
!         END IF
 
  END FUNCTION isCheckpoint
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  LOGICAL FUNCTION isEndOfThisRun()
    isEndOfThisRun = ocean_current_time >= time_config%tc_stopdate
  END FUNCTION isEndOfThisRun
  !-------------------------------------------------------------------------
 
  !-------------------------------------------------------------------------
  FUNCTION getCurrentDate_to_String() 
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: getCurrentDate_to_String
   
    CALL datetimeToString(ocean_current_time, getCurrentDate_to_String)

  END FUNCTION getCurrentDate_to_String
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  ! this will be removed, do not use !
  FUNCTION get_OceanCurrentTime_Pointer()
    TYPE(datetime), POINTER:: get_OceanCurrentTime_Pointer

    return_current_time = ocean_current_time
    get_OceanCurrentTime_Pointer => return_current_time

  END FUNCTION get_OceanCurrentTime_Pointer
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  FUNCTION newOceanEvent(description, dtime_string)
    TYPE(event), POINTER :: newOceanEvent
    CHARACTER(LEN=*) :: description
    CHARACTER(LEN=*) :: dtime_string

    TYPE(timedelta), POINTER ::  dtime
    INTEGER :: ierr

    dtime => newTimedelta(dtime_string)
    newOceanEvent => newEvent(description, eventStartDate, eventStartDate, eventEndDate, dtime, errno=ierr)

  END FUNCTION newOceanEvent
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! this is for "preprocessing" tasks
  LOGICAL FUNCTION isOceanPreEventActive(oceanEvent)
    TYPE(event), POINTER :: oceanEvent

    CHARACTER(LEN=32)               :: datestring
   ! the current time is advanced in the beginning of the loop, so for forcing we need the previous one
    isOceanPreEventActive = isCurrentEventActive(oceanEvent, ocean_previous_time) 
!     CALL datetimeToString(ocean_previous_time, datestring)
!     write(0,*) "Ocean Event at", datestring, " is ", isOceanPreEventActive

  END FUNCTION isOceanPreEventActive
  !-------------------------------------------------------------------------

END MODULE mo_ocean_time_events
