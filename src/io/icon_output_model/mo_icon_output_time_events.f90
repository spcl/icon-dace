! Contains the icon_output time keeping and events
!
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

MODULE mo_icon_output_time_events
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

  PUBLIC :: init_icon_output_time_events
  PUBLIC :: icon_output_time_nextStep
  PUBLIC :: getCurrentDate_to_String
  PUBLIC :: isStartdate
  PUBLIC :: newNullDatetime
  PUBLIC :: isEndOfThisRun
  
  PUBLIC :: get_icon_outputCurrentTime_Pointer 

!   PUBLIC :: set_icon_outputCurrentTime
  PUBLIC :: get_icon_outputCurrentTime
  
  CHARACTER(LEN=20)  :: str_module = 'mo_icon_output_time_events'  ! Output of module for 1 line debug
  !-------------------------------------------------------------------------

  TYPE(eventGroup), POINTER           :: checkpointEventGroup => NULL()

  TYPE(timedelta), POINTER            :: icon_output_time_step => NULL()

  TYPE(datetime), POINTER             :: icon_output_current_time     => NULL()
  TYPE(datetime), POINTER             :: icon_output_previous_time    => NULL()
  TYPE(datetime), POINTER             :: eventRefDate      => NULL(), &
        &                                eventStartDate    => NULL(), &
        &                                eventEndDate      => NULL()
  TYPE(datetime), POINTER             :: checkpointRefDate => NULL(), &
        &                                restartRefDate    => NULL()

  TYPE(datetime), POINTER             :: return_current_time => NULL()

  TYPE(event), POINTER                :: checkpointEvent => NULL()
  TYPE(event), POINTER                :: restartEvent    => NULL()
  
  INTEGER                             :: checkpointEvents
    
  REAL(wp):: icon_output_dtime           !< [s] length of a time step

CONTAINS

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  ! set events, group and the events
  SUBROUTINE init_icon_output_time_events()

    CHARACTER(len=*), PARAMETER :: method_name = "init_icon_output_time_events"
    
    INTEGER :: ierr
    LOGICAL :: return_status
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: dstring
    CHARACTER(len=MAX_MTIME_ERROR_STR_LEN) :: errstring
    TYPE(timedelta), POINTER               :: eventInterval   => NULL()

    icon_output_current_time  => time_config%tc_current_date
    icon_output_dtime         = config_dtime
    icon_output_previous_time => newDatetime(icon_output_current_time) ! - icon_output_dtime
    eventRefDate        => time_config%tc_exp_refdate
    eventStartDate      => time_config%tc_exp_startdate
    eventEndDate        => time_config%tc_exp_stopdate
    icon_output_time_step     => time_config%tc_dt_model

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
    CALL datetimeToString(icon_output_current_time, dstring)
    WRITE(message_text,'(a,a)') 'Start date of this run: ', dstring
    CALL message('',message_text)
    CALL datetimeToString(time_config%tc_stopdate, dstring)
    WRITE(message_text,'(a,a)') 'Stop date of this run:  ', dstring
    CALL message('',message_text)
    CALL message('','')

   CALL printEventGroup(checkpointEvents)
  
  END SUBROUTINE init_icon_output_time_events
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  FUNCTION newNullDatetime()
    TYPE(datetime), POINTER:: newNullDatetime 
    
    newNullDatetime => newDatetime('0001-01-01T00:00:00')
  END FUNCTION newNullDatetime
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  FUNCTION icon_output_time_nextStep() 
    TYPE(datetime), POINTER:: icon_output_time_nextStep 
     
    icon_output_previous_time = icon_output_current_time
    icon_output_current_time = icon_output_current_time + icon_output_time_step
    return_current_time = icon_output_current_time
    icon_output_time_nextStep => return_current_time

  END FUNCTION icon_output_time_nextStep
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  LOGICAL FUNCTION isStartdate()
    isStartdate    = (time_config%tc_startdate == icon_output_current_time)
  END FUNCTION isStartdate
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  LOGICAL FUNCTION isEndOfThisRun()
    isEndOfThisRun = icon_output_current_time >= time_config%tc_stopdate
  END FUNCTION isEndOfThisRun
  !-------------------------------------------------------------------------
 
  !-------------------------------------------------------------------------
  FUNCTION getCurrentDate_to_String() 
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: getCurrentDate_to_String
   
    CALL datetimeToString(icon_output_current_time, getCurrentDate_to_String)

  END FUNCTION getCurrentDate_to_String
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  ! this will be removed, do not use !
  FUNCTION get_icon_outputCurrentTime_Pointer()
    TYPE(datetime), POINTER:: get_icon_outputCurrentTime_Pointer

    return_current_time = icon_output_current_time
    get_icon_outputCurrentTime_Pointer => return_current_time

  END FUNCTION get_icon_outputCurrentTime_Pointer
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE set_icon_outputCurrentTime( current_time )
    TYPE(datetime), INTENT(IN) :: current_time

    icon_output_current_time = current_time

  END SUBROUTINE set_icon_outputCurrentTime
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  FUNCTION get_icon_outputCurrentTime()
    TYPE(datetime)         :: get_icon_outputCurrentTime

    get_icon_outputCurrentTime = icon_output_current_time

  END FUNCTION get_icon_outputCurrentTime
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  FUNCTION new_icon_output_event(description, dtime_string)
    TYPE(event), POINTER :: new_icon_output_event
    CHARACTER(LEN=*) :: description
    CHARACTER(LEN=*) :: dtime_string

    TYPE(timedelta), POINTER ::  dtime
    INTEGER :: ierr

    dtime => newTimedelta(dtime_string)
    new_icon_output_event => newEvent(description, eventStartDate, eventStartDate, eventEndDate, dtime, errno=ierr)

  END FUNCTION new_icon_output_event
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! this is for "preprocessing" tasks
  LOGICAL FUNCTION isicon_outputPreEventActive(icon_outputEvent)
    TYPE(event), POINTER :: icon_outputEvent

    CHARACTER(LEN=32)               :: datestring
   ! the current time is advanced in the beginning of the loop, so for forcing we need the previous one
    isicon_outputPreEventActive = isCurrentEventActive(icon_outputEvent, icon_output_previous_time) 
!     CALL datetimeToString(icon_output_previous_time, datestring)
!     write(0,*) "icon_output Event at", datestring, " is ", isicon_outputPreEventActive

  END FUNCTION isicon_outputPreEventActive
  !-------------------------------------------------------------------------

END MODULE mo_icon_output_time_events
