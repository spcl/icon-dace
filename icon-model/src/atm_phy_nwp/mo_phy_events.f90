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

! Time management for parameterized physical processes
!
! Setup mtime events and event groups for physical processes

MODULE mo_phy_events

  USE mo_kind,                     ONLY: wp
  USE mo_impl_constants,           ONLY: max_char_length, SUCCESS
  USE mo_exception,                ONLY: finish, message, message_text
  USE mtime,                       ONLY: datetime, newDatetime, timedelta, &
    &                                    datetimeToString, timedeltaToString, &
    &                                    event, newEvent, isCurrentEventActive,      &
    &                                    MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN, &
    &                                    MAX_EVENTNAME_STR_LEN, &
    &                                    OPERATOR(+), OPERATOR(-), OPERATOR(==), OPERATOR(<=), &
    &                                    getEventName, &
    &                                    getEventLastDateTime, getEventFirstDateTime,       &
    &                                    getEventInterval, getTotalSecondsTimedelta,        &
    &                                    getTriggerNextEventAtDateTime, deallocateEvent,    &
    &                                    deallocateDatetime
!!$                                   getTriggeredPreviousEventAtDateTime
  USE mo_util_mtime,               ONLY: mtime_timedelta_from_fseconds
  USE mo_util_table,               ONLY: t_table, initialize_table, add_table_column, &
    &                                    set_table_entry, print_table, finalize_table
  USE mo_run_config,               ONLY: msg_level
  USE mo_parallel_config,          ONLY: proc0_offloading
  USE mo_mpi,                      ONLY: my_process_is_stdio, p_bcast, p_comm_work, p_io
  USE mo_restart_nml_and_att,      ONLY: getAttributesForRestarting
  USE mo_key_value_store,          ONLY: t_key_value_store

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_nh_events'


  ! types
  PUBLIC :: t_phyProcFast
  PUBLIC :: t_phyProcSlow
  PUBLIC :: t_phyProcGroup

  ! subroutine
  PUBLIC :: mtime_ctrl_physics


  ! base type for individual physical process
  TYPE t_phyProcBase
    INTEGER                         :: id            !< process handle
    CHARACTER(len=MAX_CHAR_LENGTH)  :: name          !< name of physical process
    LOGICAL                         :: is_enabled    !< is this process considered at all
    TYPE(datetime)                  :: startDate     !< process startDateTime
    TYPE(datetime)                  :: endDate       !< process endDateTime
    TYPE(timedelta)                 :: dt            !< calling interval
    LOGICAL                         :: reqInit       !< requires initialization call (TRUE/FALSE)
    LOGICAL                         :: inclStart     !< include startDate
                                                     !< events are triggered between 
                                                     !< [startDate,endDate], if .TRUE.
                                                     !< ]startDate,endDate], if .FALSE.
    TYPE(timedelta)                 :: plusSlack     !< Events are triggered between 
                                                     ! [actual_trigger_time, actual_trigger_time + plus_slack]
    !
    ! mtime events
    TYPE(event), POINTER            :: ev_ptr        !< mtime event pointer
    TYPE(datetime)                  :: lastActive    !< last event date/time
  CONTAINS
    !
    ! initialize new process
    PROCEDURE  :: initialize => phyProcBase_initialize
    !
    ! reinitialize process
    PROCEDURE  :: reinitEvent => phyProcBase_reinitEvent
    !
    ! query, whether an initialization call should be triggered
    PROCEDURE  :: doInit => phyProcBase_doInit
    !
    ! query, whether the process should be triggered
    PROCEDURE  :: isActive => phyProcBase_isActive
    !
    ! get last trigger date
    PROCEDURE  :: getLastActive => phyProcBase_getLastActive
    !
    ! get last trigger date in PTString-Format
    PROCEDURE  :: getLastActivePTString => phyProcBase_getLastActivePTString
    !
    ! get time elapsed since last trigger event
    PROCEDURE  :: getElapsedTime => phyProcBase_getElapsedTime 
    !
    ! get time elapsed since last trigger event in PTString-Format
    PROCEDURE  :: getElapsedTimePTString => phyProcBase_getElapsedTimePTString
    !
    ! get next trigger date
    PROCEDURE  :: getNextActive => phyProcBase_getNextActive
    !
    ! is next trigger date within selected time range
    PROCEDURE  :: isNextTriggerTimeInRange_nobcast => phyProcBase_isNextTriggerTimeInRange_nobcast
    ! special broadcasting variant
    PROCEDURE  :: isNextTriggerTimeInRange_bcast   => phyProcBase_isNextTriggerTimeInRange_bcast
    !
    GENERIC    :: isNextTriggerTimeInRange         => isNextTriggerTimeInRange_nobcast, &
      &                                               isNextTriggerTimeInRange_bcast
    !
    ! finalization routine
    PROCEDURE  :: final    => phyProcBase_final
  END TYPE t_phyProcBase


  ! fast physics process
  TYPE, EXTENDS(t_phyProcBase) :: t_phyProcFast
  END TYPE t_phyProcFast

  ! slow physics process
  TYPE, EXTENDS(t_phyProcBase) :: t_phyProcSlow
  END TYPE t_phyProcSlow

  ! base type for creating an array/group of physical processes
  TYPE t_phyProcArr
    CLASS(t_phyProcBase), POINTER :: p
  END TYPE t_phyProcArr


  ! physical process group
  TYPE t_phyProcGroup
    CHARACTER(len=MAX_CHAR_LENGTH)    :: grpName              !< group name 
    INTEGER                           :: pid                  !< patch ID
    INTEGER                           :: ncontained           !< number of group members
    TYPE(t_phyProcArr), ALLOCATABLE   :: proc(:)              !< group of physical processes
  CONTAINS
    !
    ! construct group of physical processes
    PROCEDURE  :: construct => phyProcGroup_construct
    !
    ! add physical process to group
    PROCEDURE  :: addToGroup => phyProcGroup_addToGroup
    !
    ! reinitialize group's mtime events
    ! e.g. necessary when using IAU iteration.
    PROCEDURE  :: reinitEvents => phyProcGroup_reinitEvents
    !
    ! print basic group setup
    PROCEDURE  :: printSetup => phyProcGroup_printSetup
    !
    ! print group's current status
    PROCEDURE  :: printStatus => phyProcGroup_printStatus
    !
    ! for restart purposes
    ! translate object components into a format that can be stored in a file
    PROCEDURE  :: serialize => phyProcGroup_serialize
    !
    ! for restart purposes
    ! De-serialize object components coming from the restart file
    PROCEDURE  :: deserialize => phyProcGroup_deserialize
    !
    ! finalization routine
    PROCEDURE  :: finalize => phyProcGroup_finalize
  END TYPE t_phyProcGroup

CONTAINS




  !>
  !! Initializes a new physical process
  !!
  !! Initializes a new physical process.
  !! Creates mtime event.
  !!
  SUBROUTINE phyProcBase_initialize (phyProc, name, id, is_enabled, startDate, endDate, dt, &
    &                                plusSlack, optReqInit, optInclStart)
    CLASS(t_phyProcBase)          , INTENT(INOUT) :: phyProc     !< passed-object dummy argument
    CHARACTER(len=*)              , INTENT(IN)    :: name        !< name of physical process
    INTEGER                       , INTENT(IN)    :: id          !< process handle
    LOGICAL                       , INTENT(IN)    :: is_enabled  !< is this process considered at all
    TYPE(datetime) , TARGET       , INTENT(IN)    :: startDate   !< process startDateTime
    TYPE(datetime) , TARGET       , INTENT(IN)    :: endDate     !< process endDateTime
    TYPE(timedelta), TARGET       , INTENT(IN)    :: dt          !< calling interval
    TYPE(timedelta)               , INTENT(IN)    :: plusSlack   !< Events are triggered between
                                                                 ! [actual_trigger_time, actual_trigger_time + plus_slack]
    LOGICAL        , OPTIONAL     , INTENT(IN)    :: optReqInit  !< requires initialization call (TRUE/FALSE)
    LOGICAL        , OPTIONAL     , INTENT(IN)    :: optInclStart!< if .TRUE., startDate is included

    !
    ! local
    TYPE(datetime) , POINTER   :: eventRefDate     => NULL(), &
      &                           eventStartDate   => NULL(), &
      &                           eventEndDate     => NULL()
    TYPE(timedelta), POINTER   :: eventInterval    => NULL()
    TYPE(datetime) , POINTER   :: lastActive_ptr   => NULL()

    INTEGER :: ierr                                    ! error flag
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":phyProcBase_initialize"
  !-----------------------------------------------------------------

    ! fill in user input
    phyProc%id         = id
    phyProc%name       = TRIM(name)
    phyProc%is_enabled = is_enabled
    phyProc%startDate  = startDate
    phyProc%endDate    = endDate
    phyProc%dt         = dt
    IF (PRESENT(optReqInit)) THEN
      phyProc%reqInit = optReqInit
    ELSE
      phyProc%reqInit = .FALSE.
    ENDIF
    IF (PRESENT(optInclStart)) THEN
      phyProc%inclStart = optInclStart
    ELSE
      phyProc%inclStart = .TRUE.
    ENDIF
    ! cross check: if the process at hand is disabled, 
    ! - the endDate is set equal to the startDate
    ! - the startDate is excluded, in order to have an empty set.
    IF (.NOT. phyProc%is_enabled) THEN
      phyProc%endDate   = phyProc%startDate
      phyProc%inclStart = .FALSE.
    ENDIF
    !
! GNU: Problem with derived type optional parameters in type-bound procedures?
! Therefore, I made it a non-optional variable
    phyProc%plusSlack = plusSlack


    ! create new event
    eventStartDate => startDate
    eventEndDate   => endDate
    eventRefDate   => startDate  ! should this be set to the model start date
    eventInterval  => dt
    !
    phyProc%ev_ptr => newEvent(phyProc%name, eventRefDate, eventStartDate, &
      &                        eventEndDate, eventInterval)

    ! dummy initialization
    lastActive_ptr => newDatetime("1111-01-01T00:00:00.000")
    phyProc%lastActive = lastActive_ptr

    IF (ASSOCIATED(lastActive_ptr)) THEN
      CALL deallocateDatetime(lastActive_ptr)
    END IF

  END SUBROUTINE phyProcBase_initialize


  !>
  !! Re-initialize mtime event
  !!
  !! Re-initialize mtime event for the physical process at hand.
  !! The event is re-created based on the already available information
  !! phyProc%startDate
  !! phyProc%endDate
  !! phyProc%dt
  !! This re-initialization is, e.g., necessary if the IAU phase is iterated.
  !!
  SUBROUTINE phyProcBase_reinitEvent (phyProc)
    CLASS(t_phyProcBase), TARGET , INTENT(INOUT) :: phyProc  !< passed-object dummy argument

    !
    ! local
    TYPE(datetime) , POINTER   :: eventRefDate     => NULL(), &
      &                           eventStartDate   => NULL(), &
      &                           eventEndDate     => NULL()
    TYPE(timedelta), POINTER   :: eventInterval    => NULL()
    TYPE(datetime) , POINTER   :: lastActive_ptr   => NULL()

    INTEGER :: ierr                                    ! error flag
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":phyProcBase_initialize"
  !-----------------------------------------------------------------

    IF (ASSOCIATED(phyProc%ev_ptr)) THEN
      CALL deallocateEvent(phyProc%ev_ptr)
    ENDIF

    ! create new event
    eventStartDate => phyProc%startDate
    eventEndDate   => phyProc%endDate
    eventRefDate   => phyProc%startDate  ! should this be set to the model start date
    eventInterval  => phyProc%dt
    !
    phyProc%ev_ptr => newEvent(phyProc%name, eventRefDate, eventStartDate, &
      &                        eventEndDate, eventInterval, errno=ierr)

    ! dummy initialization
    lastActive_ptr => newDatetime("1111-01-01T00:00:00.000")
    phyProc%lastActive = lastActive_ptr
    CALL deallocateDatetime(lastActive_ptr)

  END SUBROUTINE phyProcBase_reinitEvent


  !>
  !! Checks, whether an initialization should be triggered
  !!
  !! Checks, whether an initialization should be triggered for 
  !! the physical process at hand.
  !!
  !!
  LOGICAL FUNCTION phyProcBase_doInit (phyProc, mtime_current)

    CLASS(t_phyProcBase)   , INTENT(INOUT) :: phyProc        !< passed-object dummy argument
    TYPE(datetime)         , INTENT(IN)    :: mtime_current  !< current_datetime

  !-----------------------------------------------------------------


    IF (phyProc%is_enabled) THEN
      IF (phyProc%reqInit) THEN
        phyProcBase_doInit = .TRUE.
        phyProc%lastActive = mtime_current
      ELSE
        phyProcBase_doInit = .FALSE.
      ENDIF
    ELSE
      phyProcBase_doInit = .FALSE.
    ENDIF

  END FUNCTION phyProcBase_doInit



  !>
  !! Checks, whether the process should be triggered
  !!
  !! Checks, whether the process should be triggered.
  !!
  !!
  LOGICAL FUNCTION phyProcBase_isActive (phyProc, mtime_current, lasync)

    CLASS(t_phyProcBase)   , INTENT(INOUT) :: phyProc        !< passed-object dummy argument
    TYPE(datetime)         , INTENT(IN)    :: mtime_current  !< current_datetime
    LOGICAL, INTENT(IN), OPTIONAL          :: lasync         !< if present, broadcast is done in calling routine

    ! local
    TYPE(timedelta), TARGET  :: plusSlack
    TYPE(timedelta), POINTER :: plusSlack_ptr    => NULL()
    LOGICAL                  :: is_active
  !-----------------------------------------------------------------

    plusSlack     =  phyProc%plusSlack
    plusSlack_ptr => plusSlack
    ! In NEC hybrid mode, execute isCurrentEventActive only on PE0 and broadcast result afterwards
    ! This serves as a workaround for slow mtime on the NEC VEs
    IF (.NOT. proc0_offloading .OR. my_process_is_stdio()) THEN
      is_active = isCurrentEventActive(phyProc%ev_ptr, mtime_current, plus_slack=plusSlack_ptr)
    ENDIF
    IF (.NOT. PRESENT(lasync) .AND. proc0_offloading) CALL p_bcast(is_active, p_io, p_comm_work)

    phyProcBase_isActive = is_active

    ! Mtime currently does not support choosing between open or closed intervals.
    ! Therefore, the result of isCurrentEventActive is overwritten 
    ! by phyProc%inclStart, if this is the startDate.
    IF (mtime_current == phyProc%startDate) THEN
      phyProcBase_isActive = phyProc%inclStart
    ENDIF

    ! store last trigger date
    IF (phyProcBase_isActive)  phyProc%lastActive = mtime_current

  END FUNCTION phyProcBase_isActive


!!$  !>
!!$  !! Get last trigger date
!!$  !!
!!$  !! Get last trigger date
!!$  !! Currently not restart-safe!
!!$  !!
!!$  TYPE(datetime) FUNCTION phyProcBase_getLastActive (phyProc)
!!$
!!$    CLASS(t_phyProcBase), INTENT(IN) :: phyProc       !< passed-object dummy argument
!!$
!!$    ! local
!!$    TYPE(datetime) :: lastActive
!!$    INTEGER        :: ierr
!!$  !-----------------------------------------------------------------
!!$
!!$    CALL getTriggeredPreviousEventAtDateTime(phyProc%ev_ptr, lastActive, ierr)
!!$
!!$    phyProcBase_getLastActive = lastActive
!!$
!!$  END FUNCTION phyProcBase_getLastActive


  !>
  !! Get last trigger date
  !!
  !! Get last trigger date
  !! Trivial function so far. Last trigger date is part of  
  !! t_phyProcBase. Actually it would be cleaner to use the 
  !! mtime query function getTriggeredPreviousEventAtDateTime.
  !! However, this function is not restart-safe. The information 
  !! about the last trigger date is lost after restart. 
  !!
  TYPE(datetime) FUNCTION phyProcBase_getLastActive (phyProc)

    CLASS(t_phyProcBase), INTENT(IN) :: phyProc    !< passed-object dummy argument

  !-----------------------------------------------------------------

    phyProcBase_getLastActive = phyProc%lastActive

  END FUNCTION phyProcBase_getLastActive


  !>
  !! Get last trigger date in PTString-Format
  !!
  !! Get last trigger date in PTString-Format
  !!
  CHARACTER(len=MAX_DATETIME_STR_LEN) FUNCTION phyProcBase_getLastActivePTString (phyProc)

    CLASS(t_phyProcBase)   , INTENT(IN) :: phyProc    !< passed-object dummy argument

    ! local
    TYPE(datetime) :: lastActive
    INTEGER :: ierr
  !-----------------------------------------------------------------

    lastActive = phyProc%getLastActive()

    CALL datetimeToString(lastActive, phyProcBase_getLastActivePTString, ierr)

  END FUNCTION phyProcBase_getLastActivePTString



  !>
  !! Compute time elapsed since last trigger event
  !!
  !! Compute time elapsed since last trigger event.
  !!
  !!
  TYPE(timedelta) FUNCTION phyProcBase_getElapsedTime (phyProc, mtime_current)

    CLASS(t_phyProcBase)   , INTENT(IN) :: phyProc        !< passed-object dummy argument
    TYPE(datetime)         , INTENT(IN) :: mtime_current  !< current_datetime

  !-----------------------------------------------------------------

    phyProcBase_getElapsedTime = mtime_current - phyProc%getLastActive()

  END FUNCTION phyProcBase_getElapsedTime


  !>
  !! Get time elapsed since last trigger event in PTString-format
  !!
  !! Get time elapsed since last trigger event in PTString-format
  !!
  CHARACTER(len=MAX_TIMEDELTA_STR_LEN) FUNCTION phyProcBase_getElapsedTimePTString (phyProc, mtime_current)

    CLASS(t_phyProcBase)   , INTENT(IN) :: phyProc        !< passed-object dummy argument
    TYPE(datetime)         , INTENT(IN) :: mtime_current  !< current_datetime

    ! local
    TYPE(timedelta), TARGET  :: elapsedTime
    INTEGER :: ierr
  !-----------------------------------------------------------------

    elapsedTime = phyProc%getElapsedTime(mtime_current)

    CALL timedeltaToString(elapsedTime, phyProcBase_getElapsedTimePTString, ierr)

  END FUNCTION phyProcBase_getElapsedTimePTString



  !>
  !! Get next trigger date.
  !!
  !! Get next trigger date.
  !!
  TYPE(datetime) FUNCTION phyProcBase_getNextActive (phyProc, mtime_current)

    CLASS(t_phyProcBase)   , INTENT(IN) :: phyProc        !< passed-object dummy argument
    TYPE(datetime), POINTER, INTENT(IN) :: mtime_current  !< current_datetime

    ! local
    TYPE(datetime) :: nextActive
    INTEGER        :: ierr
  !-----------------------------------------------------------------


    CALL getTriggerNextEventAtDateTime(phyProc%ev_ptr,mtime_current,nextActive,ierr)

    phyProcBase_getNextActive = nextActive

  END FUNCTION phyProcBase_getNextActive


  !>
  !! Is the next trigger date within selected time range.
  !!
  !! Is the next trigger date within selected time range.
  !! The time range is given by [mtime_current, mtime_current+slack].
  !!
  LOGICAL FUNCTION phyProcBase_isNextTriggerTimeInRange_nobcast (phyProc, mtime_current, slack)

    CLASS(t_phyProcBase)   , INTENT(IN) :: phyProc        !< passed-object dummy argument
    TYPE(datetime)         , INTENT(IN) :: mtime_current  !< current_datetime
    TYPE(timedelta)        , INTENT(IN) :: slack          !< time frame

    ! local
    TYPE(datetime) :: nextActive
    TYPE(datetime) :: timeFrameEnd
    INTEGER        :: ierr
    LOGICAL        :: lret
  !-----------------------------------------------------------------

    IF (.NOT. proc0_offloading .OR. my_process_is_stdio()) THEN
      CALL getTriggerNextEventAtDateTime(phyProc%ev_ptr, mtime_current, nextActive, ierr)
      timeFrameEnd = mtime_current + slack
      lret = nextActive <= timeFrameEnd
    ENDIF

    IF (proc0_offloading) CALL p_bcast(lret, p_io, p_comm_work)

    phyProcBase_isNextTriggerTimeInRange_nobcast = lret

  END FUNCTION phyProcBase_isNextTriggerTimeInRange_nobcast


  !>
  !! Is the next trigger date within selected time range.
  !!
  !! Is the next trigger date within selected time range.
  !! The time range is given by [mtime_current, mtime_current+slack].
  !!
  !! This version broadcasts mtime_current and nextActive from 
  !! p_source. It is particularly suited for the case of nonzero 
  !! patch weights.
  !!
  LOGICAL FUNCTION phyProcBase_isNextTriggerTimeInRange_bcast (phyProc, mtime_current, slack, p_source, comm)

    CLASS(t_phyProcBase)   , INTENT(IN)    :: phyProc        !< passed-object dummy argument
    TYPE(datetime)         , INTENT(INOUT) :: mtime_current  !< current_datetime
    TYPE(timedelta)        , INTENT(IN)    :: slack          !< time frame
    INTEGER                , INTENT(IN)    :: p_source
    INTEGER                , INTENT(IN)    :: comm

    ! local
    TYPE(datetime) :: nextActive
    TYPE(datetime) :: timeFrameEnd
    INTEGER        :: ierr
  !-----------------------------------------------------------------

    CALL p_bcast(mtime_current, p_source, comm)
    !
    CALL getTriggerNextEventAtDateTime(phyProc%ev_ptr, mtime_current, nextActive, ierr)
    !
    CALL p_bcast(nextActive, p_source, comm)

    timeFrameEnd = mtime_current + slack
    phyProcBase_isNextTriggerTimeInRange_bcast = (nextActive <= timeFrameEnd)

  END FUNCTION phyProcBase_isNextTriggerTimeInRange_bcast



  !>
  !! Finalize physical process
  !!
  SUBROUTINE phyProcBase_final    (phyProc)
    CLASS(t_phyProcBase), INTENT(INOUT) :: phyProc    !< passed-object dummy argument

    ! local
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":phyProcBase_final"
  !-----------------------------------------------------------------

    IF (ASSOCIATED(phyProc%ev_ptr)) THEN
      CALL deallocateEvent(phyProc%ev_ptr)
    ENDIF

  END SUBROUTINE phyProcBase_final




  !>
  !! Constructor for variable of type t_phyProcGroup
  !!
  !! Constructor for variable of type t_phyProcGroup
  !!
  SUBROUTINE phyProcGroup_construct (phyProcGrp, pid, grpName, grpSize)
    CLASS(t_phyProcGroup)         , INTENT(INOUT)   :: phyProcGrp !< passed-object dummy argument
    INTEGER                       , INTENT(IN)      :: pid        !< patch ID
    CHARACTER(len=*)              , INTENT(IN)      :: grpName    !< physical process group name
    INTEGER                       , INTENT(IN)      :: grpSize    !< size of group to be constructed

    ! local
    INTEGER :: error                                   ! error flag
    INTEGER :: i                                       ! loop index
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":phyProcGroup_construct"
  !-----------------------------------------------------------------
 

    ! store user input
    phyProcGrp%pid        = pid
    phyProcGrp%grpName    = TRIM(grpName)
    phyProcGrp%ncontained = 0


    ! allocate group according to the given size
    ALLOCATE(phyProcGrp%proc(grpSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    ! disassociate
    DO i=1, UBOUND(phyProcGrp%proc,1)
      phyProcGrp%proc(i)%p => NULL()
    ENDDO

  END SUBROUTINE phyProcGroup_construct



  !>
  !! Add process to group
  !!
  !! Add process to group.
  !!
  SUBROUTINE phyProcGroup_addToGroup (phyProcGrp, phyProc)
    CLASS(t_phyProcGroup)       , INTENT(INOUT)   :: phyProcGrp !< passed-object dummy argument
    CLASS(t_phyProcBase), TARGET, INTENT(IN)      :: phyProc    !< physical process

    ! local
    INTEGER :: id                                      ! process ID
    INTEGER :: error                                   ! error flag
    TYPE(t_phyProcArr), ALLOCATABLE :: tmp_proc(:)     ! temporary array
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":phyProcGroup_addToGroup"
  !-----------------------------------------------------------------
 
    id = phyProc%id

    phyProcGrp%ncontained = phyProcGrp%ncontained + 1

    ! expand process array, if necessary.
    IF (phyProcGrp%ncontained > UBOUND(phyProcGrp%proc,1)) THEN
      ALLOCATE(tmp_proc(phyProcGrp%ncontained), STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
      
      tmp_proc(1:UBOUND(phyProcGrp%proc,1)) = phyProcGrp%proc
      CALL move_alloc(from=tmp_proc, to=phyProcGrp%proc)
      WRITE(message_text,'(a,a,a,i3)') 'Increased size of process group ',TRIM(phyProcGrp%grpName), &
        &                               ' to ', SIZE(phyProcGrp%proc,1)
      CALL message('', message_text)
    ENDIF
    ! Sanity checks
    IF (id > UBOUND(phyProcGrp%proc,1)) THEN
      WRITE(message_text,'(a,i3,a,i3)') 'requested storage index ',id,                &
        &                               ' exceeds upper bound of storage container ', &
        &                               UBOUND(phyProcGrp%proc,1)
      CALL message('', TRIM(message_text))
    ENDIF

    ! add to group
    phyProcGrp%proc(id)%p => phyProc

  END SUBROUTINE phyProcGroup_addToGroup



  !>
  !! All group's mtime events are re-initialized
  !!
  !! The mtime-events of all group members are re-initialized 
  !! by calling the member-specific routine reinitEvent.
  !!
  SUBROUTINE phyProcGroup_reinitEvents (phyProcGrp)
    CLASS(t_phyProcGroup), INTENT(INOUT)     :: phyProcGrp  !< passed-object dummy argument

    ! local
    INTEGER :: iproc                                   ! process ID
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":phyProcGroup_addToGroup"
  !-----------------------------------------------------------------

    IF (msg_level >= 12) THEN
      WRITE(message_text,'(a,a,a,i2)') 'Reinitialize events for ', TRIM(phyProcGrp%grpName), &
        &                              ' on patch ', phyProcGrp%pid
      CALL message('', TRIM(message_text))
    ENDIF

    DO iproc = 1, UBOUND(phyProcGrp%proc,1)
      IF (.NOT. ASSOCIATED(phyProcGrp%proc(iproc)%p)) CYCLE
      !
      ! re-initialize event for the process at hand
      CALL phyProcGrp%proc(iproc)%p%reinitEvent()
    ENDDO

  END SUBROUTINE phyProcGroup_reinitEvents



  !>
  !! Translate object into a format that can be stored in a file
  !!
  !! Translate object into a format that can be stored in a file.
  !! Currently this is only done for those components that must be 
  !! stored for restart.
  !! The following components are serialized:
  !! phyProcGrp%proc(:)%p%lastActive
  !! It is transformed into 'seconds since last trigger date'. 
  !! From this, phyProcGrp%proc(:)%p%lastActive can be restored after restart.
  !!
  SUBROUTINE phyProcGroup_serialize (phyProcGrp, mtime_current, elapsedTime)
    CLASS(t_phyProcGroup)         , INTENT(IN)   :: phyProcGrp     !< passed-object dummy argument
    TYPE(datetime), POINTER       , INTENT(IN)   :: mtime_current  !< current_datetime
    REAL(wp), ALLOCATABLE         , INTENT(INOUT):: elapsedTime(:)

    ! local
    INTEGER         :: iproc  ! loop counter
    INTEGER         :: ierr
    TYPE(timedelta) :: td ! elapsed time since last trigger date
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":phyProcGroup_serialize"

  !-----------------------------------------------------------------
 
    IF (.NOT.ALLOCATED(elapsedTime)) THEN
      ALLOCATE(elapsedTime(UBOUND(phyProcGrp%proc,1)), STAT=ierr)
      IF(ierr /= SUCCESS) CALL finish(routine, "memory allocation failure")
    ENDIF

    ! initialize
    elapsedTime(:) = -999._wp

    ! get elapsed time since last trigger date in s
    IF (.NOT. proc0_offloading .OR. my_process_is_stdio()) THEN
      DO iproc = 1, UBOUND(phyProcGrp%proc,1)
        IF (.NOT. ASSOCIATED(phyProcGrp%proc(iproc)%p)) CYCLE
        IF (.NOT. phyProcGrp%proc(iproc)%p%is_enabled) CYCLE
        !
        td = phyProcGrp%proc(iproc)%p%getElapsedTime(mtime_current)
        elapsedTime(iproc) = REAL(getTotalSecondsTimedelta(td,mtime_current),wp)
      ENDDO
    ENDIF
    IF (proc0_offloading) CALL p_bcast(elapsedTime, p_io, p_comm_work)

  END SUBROUTINE phyProcGroup_serialize



  !>
  !! De-serialize object components coming from the restart file
  !!
  !! De-serialize object components coming from the restart file.
  !! So far, the following components are de-serialized:
  !! elapsedTime
  !! From this, phyProcGrp%proc(:)%p%lastActive is restored.
  !!
  SUBROUTINE phyProcGroup_deserialize (phyProcGrp, mtime_current, optAttnamePrefix)
    CLASS(t_phyProcGroup)         , INTENT(INOUT):: phyProcGrp     !< passed-object dummy argument
    !> current_datetime
    TYPE(datetime)                , INTENT(IN)   :: mtime_current
    CHARACTER(LEN=*), OPTIONAL    , INTENT(IN)   :: optAttnamePrefix !< optional prefix in attribute name

    ! local
    INTEGER                               :: iproc               ! loop conter
    REAL(wp)                              :: elapsedTime
    TYPE(timedelta)                       :: mtime_elapsedTime   ! elapsed time in mtime format
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CHARACTER(len=MAX_CHAR_LENGTH)        :: attname             ! attribute name
    CHARACTER(LEN=*), PARAMETER           :: routine = modname//":phyProcGroup_deserialize"
  !-----------------------------------------------------------------
 
    CALL getAttributesForRestarting(restartAttributes)

    IF (restartAttributes%is_init) THEN
      DO iproc=1,UBOUND(phyProcGrp%proc,1)
        IF (.NOT. ASSOCIATED(phyProcGrp%proc(iproc)%p)) CYCLE
        IF (.NOT. phyProcGrp%proc(iproc)%p%is_enabled) CYCLE

        IF (.NOT. PRESENT(optAttnamePrefix)) THEN
          WRITE(attname,'(a,i2.2,a,i2.2)') 't_elapsed_phy_DOM',phyProcGrp%pid,'_PHY',iproc
        ELSE
          WRITE(attname,'(a,i2.2,a,i2.2)') TRIM(optAttnamePrefix),phyProcGrp%pid,'_PHY',iproc
        ENDIF
        CALL restartAttributes%get(attname, elapsedTime)

        ! Note that elapsedTime is multiplied by -1, since we only have a 
        ! '+' operator available.
        CALL mtime_timedelta_from_fseconds(-elapsedTime, mtime_current, &
          &                                mtime_elapsedTime)
        phyProcGrp%proc(iproc)%p%lastActive = mtime_current + mtime_elapsedTime
        !
      ENDDO
    ENDIF  ! associated

  END SUBROUTINE phyProcGroup_deserialize


  !>
  !! Print setup for all group members
  !!
  !! Print setup for all group members.
  !! I.e. show details about the mtime-events.
  !!
  SUBROUTINE phyProcGroup_printSetup (phyProcGrp)
    !
    CLASS(t_phyProcGroup), INTENT(INOUT) :: phyProcGrp   !< passed-object dummy argument

    ! local variables
    INTEGER         :: iev             ! loop counter
    TYPE(t_table)   :: table
    INTEGER         :: irow            ! row to fill
    !
    CHARACTER(LEN = *), PARAMETER :: eventNameCol = "eventName", &
      &                              startDateCol = "startDate", &
      &                              endDateCol   = "endDate",   &
      &                              intvlCol     = "intvl",     &
      &                              enabledCol   = "enabled"

    CHARACTER(len=MAX_EVENTNAME_STR_LEN) :: eventName
    CHARACTER(len=MAX_DATETIME_STR_LEN)  :: startDate_str, endDate_str
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: intvl_str
    CHARACTER(len=2)                     :: inclStart_str, inclEnd_str  ! open/closed intervals

    INTEGER :: ierr
#ifdef __INTEL_COMPILER
    INTEGER :: nrows ! HB : workaround for ifort16.0 weirdness aka bug
#endif
    TYPE(datetime), POINTER   :: startDate_ptr, endDate_ptr
    TYPE(timedelta), POINTER  :: intvl_ptr
    CHARACTER(LEN=3)          :: enabled_str   ! is the process enabled at all (TRUE/FALSE)
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":events_print_status"
    !--------------------------------------------------------------------------

    ! will only be executed by stdio process
    IF(.NOT. my_process_is_stdio()) RETURN

    ! could this be transformed into a table header? 
    WRITE(message_text,'(a,a,a,i2)') 'Event-setup for ',TRIM(phyProcGrp%grpName), &
      &                              ' on patch ', phyProcGrp%pid
    CALL message('', TRIM(message_text))


    ! table-based output
    CALL initialize_table(table)
    ! the latter is no longer mandatory
#ifdef __INTEL_COMPILER
    nrows = 0
    DO iev=1,UBOUND(phyProcGrp%proc,1)   
      IF (ASSOCIATED(phyProcGrp%proc(iev)%p)) &
        & nrows = nrows + 1
    END DO
    CALL add_table_column(table, eventNameCol, nrows)
    CALL add_table_column(table, enabledCol, nrows)
    CALL add_table_column(table, startDateCol, nrows)
    CALL add_table_column(table, endDateCol, nrows)
    CALL add_table_column(table, intvlCol, nrows)
#else
    CALL add_table_column(table, eventNameCol)
    CALL add_table_column(table, enabledCol)
    CALL add_table_column(table, startDateCol)
    CALL add_table_column(table, endDateCol)
    CALL add_table_column(table, intvlCol)
#endif


    ! initialize counter
    irow = 0

    DO iev=1,UBOUND(phyProcGrp%proc,1)
      ! print NWP-phy event information
      !
      IF (.NOT. ASSOCIATED(phyProcGrp%proc(iev)%p)) CYCLE

      irow = irow + 1
      !
      CALL getEventName(phyProcGrp%proc(iev)%p%ev_ptr, eventName, ierr) 
      CALL set_table_entry(table,irow,eventNameCol, ADJUSTL(TRIM(eventName)))
      !
      ! print whether the process is enabled
      enabled_str = MERGE('T','F', phyProcGrp%proc(iev)%p%is_enabled)
      CALL set_table_entry(table,irow,enabledCol, ADJUSTL(TRIM(enabled_str)))
      !
      ! print event start date
      startDate_ptr => getEventFirstDateTime(phyProcGrp%proc(iev)%p%ev_ptr)
      CALL dateTimeToString(startDate_ptr, startDate_str, ierr)
      inclStart_str = MERGE("[","]",phyProcGrp%proc(iev)%p%inclStart)
      CALL set_table_entry(table,irow,startDateCol, TRIM(inclStart_str)//ADJUSTL(TRIM(startDate_str)))
      !
      ! print event end date
      endDate_ptr   => getEventLastDateTime(phyProcGrp%proc(iev)%p%ev_ptr)
      CALL dateTimeToString(endDate_ptr, endDate_str, ierr)
      inclEnd_str = "]"
      CALL set_table_entry(table,irow,endDateCol, ADJUSTL(TRIM(endDate_str))//TRIM(inclEnd_str))
      !
      ! print event interval
      intvl_ptr => getEventInterval(phyProcGrp%proc(iev)%p%ev_ptr)
      CALL timedeltaToString(intvl_ptr, intvl_str, ierr)
      CALL set_table_entry(table,irow,intvlCol, ADJUSTL(TRIM(intvl_str)))

    ENDDO

    CALL print_table(table, opt_delimiter=' | ')
    CALL finalize_table(table)

    WRITE (0,*) " " ! newline
  END SUBROUTINE phyProcGroup_printSetup


  !>
  !! Print status of all group members
  !!
  !! Print status of all group members.
  !! I.e. details about the mtime-events are printed.
  !!
  SUBROUTINE phyProcGroup_printStatus (phyProcGrp, mtime_current)
    !
    CLASS(t_phyProcGroup)  , INTENT(INOUT) :: phyProcGrp    !< passed-object dummy argument
    TYPE(datetime)         , INTENT(IN)    :: mtime_current !< current_datetime

    ! local variables
    INTEGER         :: iev             ! loop counter
    TYPE(t_table)   :: table
    INTEGER         :: irow            ! row to fill
    !
    CHARACTER(LEN = *), PARAMETER :: eventNameCol   = "eventName",   & 
      &                              triggerNowCol  = "trigger now", &
      &                              lastActiveCol  = "lastActive",  &
      &                              elapsedTimeCol = "elapsedTime"

    !
    CHARACTER(len=MAX_EVENTNAME_STR_LEN) :: eventName
    CHARACTER(len=MAX_DATETIME_STR_LEN)  :: lastActive_str
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: elapsed_str

    INTEGER :: ierr
    CHARACTER(LEN=3)          :: trigger_str   ! triggered this timestep (TRUE/FALSE)
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":events_print_status"
    !--------------------------------------------------------------------------

    ! will only be executed by stdio process
    IF(.NOT. my_process_is_stdio()) RETURN

    ! could this be transformed into a table header? 
    WRITE(message_text,'(a,a,i2)') TRIM(phyProcGrp%grpName),' events debug output for patch ', &
      &                          phyProcGrp%pid
    CALL message('', TRIM(message_text))


    ! table-based output
    CALL initialize_table(table)
    ! the latter is no longer mandatory
    CALL add_table_column(table, eventNameCol)
    CALL add_table_column(table, triggerNowCol)
    CALL add_table_column(table, lastActiveCol)
    CALL add_table_column(table, elapsedTimeCol)


    ! initialize counter
    irow = 0

    DO iev=1,UBOUND(phyProcGrp%proc,1)
      ! print NWP-phy event information
      !
      IF (.NOT. ASSOCIATED(phyProcGrp%proc(iev)%p)) CYCLE

      irow = irow + 1
      !
      CALL getEventName(phyProcGrp%proc(iev)%p%ev_ptr, eventName, ierr) 
      CALL set_table_entry(table,irow,eventNameCol, ADJUSTL(TRIM(eventName)))
      !
      ! is the process triggered at the current time step
      trigger_str = MERGE('T','F',phyProcGrp%proc(iev)%p%getLastActive() == mtime_current)
      CALL set_table_entry(table,irow,triggerNowCol, ADJUSTL(TRIM(trigger_str)))
      !
      ! print date of last triggering-event
      lastActive_str = phyProcGrp%proc(iev)%p%getLastActivePTString()
      CALL set_table_entry(table,irow,lastActiveCol, ADJUSTL(TRIM(lastActive_str)))
      !
      ! print elapsed time since last triggering-event
      elapsed_str = phyProcGrp%proc(iev)%p%getElapsedTimePTString(mtime_current)
      CALL set_table_entry(table,irow,elapsedTimeCol, TRIM(elapsed_str))

    ENDDO

    CALL print_table(table, opt_delimiter=' | ')
    CALL finalize_table(table)

    WRITE (0,*) " " ! newline
  END SUBROUTINE phyProcGroup_printStatus


  !>
  !! Finalization for variable of type t_phyProcGroup
  !!
  !! Finalization for variable of type t_phyProcGroup
  !!
  SUBROUTINE phyProcGroup_finalize (phyProcGrp)
    CLASS(t_phyProcGroup) , INTENT(INOUT) :: phyProcGrp  !< passed-object dummy argument

    ! local
    INTEGER :: iproc                                 ! loop counter
    INTEGER :: ierrstat                              ! error flag
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":phyProcGroup_finalize"
  !-----------------------------------------------------------------
 
    DO iproc=1,UBOUND(phyProcGrp%proc,1)
      IF (.NOT. ASSOCIATED(phyProcGrp%proc(iproc)%p)) CYCLE
      CALL phyProcGrp%proc(iproc)%p%final()
    ENDDO

    IF (ALLOCATED(phyProcGrp%proc))  DEALLOCATE(phyProcGrp%proc, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE phyProcGroup_finalize


  !-------------------------------------------------------------------------
  !>
  !! Physics time control
  !!
  !! Physics time control. This function returns a 1D array
  !! of size SIZE(phyProcs%proc,1) and type LOGICAL with one entry per  
  !! physical process. 
  !! If lcall_phy(iphys)=.TRUE., the physical process at hand must be 
  !! called this time step. If lcall_phy(iphys)=.FALSE., it must not be called.
  !! The decision making is based upon mtime events, which are initialized in 
  !! mo_atm_phy_nwp_config:configure_atm_phy_nwp.
  !!
  SUBROUTINE mtime_ctrl_physics ( phyProcs, mtime_current, isInit, lcall_phy )

    TYPE(t_phyProcGroup)    , INTENT(INOUT):: phyProcs       !< physics group
    TYPE(datetime)          , INTENT(IN)   :: mtime_current  !< current_datetime
    LOGICAL                 , INTENT(IN)   :: isInit         !< if TRUE, special settings 
                                                             !  for lcall_phy which are adjusted 
                                                             !  for the physics initialization phase
    LOGICAL                , INTENT(INOUT) :: lcall_phy(:)   !< trigger information

    ! local
    INTEGER :: iproc      ! loop counter
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":mtime_ctrl_physics"
    LOGICAL :: lasync = .TRUE.
  !-----------------------------------------------------------------

    ! intialization
    lcall_phy(:) =.FALSE.

    IF (isInit) THEN
      DO iproc = 1, UBOUND(phyProcs%proc,1)
        IF (.NOT. ASSOCIATED(phyProcs%proc(iproc)%p)) CYCLE
        lcall_phy(iproc) = phyProcs%proc(iproc)%p%doInit(mtime_current)
      ENDDO
    ELSE
      DO iproc = 1, UBOUND(phyProcs%proc,1)
        IF (.NOT. ASSOCIATED(phyProcs%proc(iproc)%p)) CYCLE
        lcall_phy(iproc) = phyProcs%proc(iproc)%p%isActive(mtime_current,lasync)
      ENDDO
      IF (proc0_offloading) CALL p_bcast(lcall_phy, p_io, p_comm_work)
    ENDIF

    !$ACC UPDATE DEVICE(lcall_phy) ASYNC(1)

    ! debug output
    IF (msg_level >= 13) THEN
      CALL phyProcs%printStatus(mtime_current)
    ENDIF

  END SUBROUTINE mtime_ctrl_physics


END MODULE mo_phy_events

