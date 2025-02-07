! mo_time_manager.f90 - Conduct time management operations
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_time_manager
!BOP
  ! !MODULE: mo_time_manager
  ! !DESCRIPTION:
  !   time managment operations

  ! !REVISION HISTORY:
  ! Version: E5/R1.07+ 04-February-2003
  !
  ! Authors:
  !   I. Kirchner,  MPI, May 2000
  !   L. Kornblueh, MPI, September 2000
  !   I. Kirchner,  MPI, November 2000, revision
  !   I. Kirchner,  MPI, March 2001, revision
  !   I. Kirchner,  MPI, August 2002, add feature lock/unlock manager and rewind manager
  !   I. Kirchner,  FUB, February 2003, revision/code review
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.
  !
  ! !USES:
  USE mo_kind,            ONLY: dp
  USE mo_exception,       ONLY: finish, message
  USE mo_time_conversion, ONLY: TC_convert, TC_get, OPERATOR(==), OPERATOR(<), &
                                time_days, time_intern, time_native, &
                                add_date, print_date, TC_PRN_NATIVE
  USE mo_time_base,       ONLY: idaylen
!EOP
!BOC
!BOX
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: STR_LEN_A = 100
  INTEGER, PARAMETER :: STR_LEN_C = 256

  CHARACTER(len=STR_LEN_C) :: message_text
!EOX
  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: time_manager
  TYPE time_manager
    PRIVATE                                               ! no external access possible
    LOGICAL              :: init         = .FALSE.        ! manager state for access control
    LOGICAL              :: freeze       = .FALSE.        ! lock/unlock the manager
    CHARACTER(STR_LEN_A) :: label        = ''             ! short description of the manager
    TYPE (time_days)     :: start_date                    ! initial date at pos=0
    INTEGER              :: pos          = 0              ! present step of the manager
    REAL(dp)             :: delta_time   = 0._dp          ! distance between two points in seconds
  END TYPE time_manager

  ! *************** INTERFACE SUBROUTINES AND FUNCTIONS
  !
  ! input: label, startdate, increment
  !                different formats
  !        only a date will reset the start date
  !        only a increment reset the increment and the start date
  !                           
  PUBLIC :: manager_init    !   set/reset manager state
  !   PARAMETERS(
  !     manager [time_manager] (inp, out) the time axis manager
  !     label   [character]    (inp)      name of the manager
  !     date    [time_days|time_intern|time_native] 
  !                            (inp)      start date
  !     second  [integer]      (inp)      distance between two points
  !     zsecond [real]         (inp)      distance between two points
  !     offset  [integer]      (inp)      point offset moving position
  !     lfreeze [logical]      (inp)      lock/unlock feature
  !   )
  !
  !   separation of functions using the type and number of parameter
  !
  INTERFACE manager_init
    !   manager, label, date, second  -> initialisation
    MODULE PROCEDURE TE_manager_init_days
    MODULE PROCEDURE TE_manager_init_native
    MODULE PROCEDURE TE_manager_init_intern
    !   manager, date                 -> reset start date, adjust position
    MODULE PROCEDURE TE_manager_reinit_days
    MODULE PROCEDURE TE_manager_reinit_native
    MODULE PROCEDURE TE_manager_reinit_intern
    !   manager, zsecond              -> reset increment, adjust start date
    MODULE PROCEDURE TE_manager_reinit_incr
    !   manager, offset               -> set new position using offset
    MODULE PROCEDURE TE_manager_reinit_step
    !   manager, lfreeze              -> lock/unlock time manager
    MODULE PROCEDURE TE_manager_freeze
  END INTERFACE

  PUBLIC :: rewind_manager
  PUBLIC :: manager_print
  PUBLIC :: manager_state
  !
  ! manager_state [subroutine, interface]
  !   get informations about the manager state
  !   (
  !   manager  [time_manager] (inp)      the time axis manager
  !   date     [time_days]    (inp, out) calendar date
  !   position [integer]      (out)      position on the manager
  !   lequal   [logical]      (out)      date fit to manager position
  !   )
  !   manager, date, position         -> get date for position at manager
  !   manager, date                   -> get present date of manager
  !   manager, date, position, lequal -> get nearest position for date
  !   manager, dtime                  -> get delta time
  !   manager, pos                    -> get time step
  !
  INTERFACE manager_state
    MODULE PROCEDURE TE_manager_any_date
    MODULE PROCEDURE TE_manager_current_date
    MODULE PROCEDURE TE_manager_step_at_date
    MODULE PROCEDURE TE_manager_dtime
    MODULE PROCEDURE TE_manager_pos
  END INTERFACE
!EOC
CONTAINS

  !---------------------------------------------------------------------
  ! manager manipulations
  !
!BOP
  ! !IROUTINE: TE_manager_init_days
  ! !INTERFACE:
  SUBROUTINE TE_manager_init_days (manager, name, date, incr)
    ! !DESCRIPTION:
    ! initialisation - reinitialisation of a time manager

    ! !INPUT/OUTPUT PARAMETERS:
    TYPE (time_manager), INTENT(inout) :: manager
    CHARACTER(len=*),    INTENT(in)    :: name
    TYPE (time_days),    INTENT(in)    :: date
    INTEGER,             INTENT(in)    :: incr
!EOP
    INTEGER    :: is

    IF (manager%freeze) THEN
      WRITE(message_text,*) 'Time manager <',&
           TRIM(manager%label),'> is frozen, skip initialization.'

    ELSE IF (manager%init) THEN
      WRITE(message_text,*) 'Time manager <',&
           TRIM(manager%label), '> was initialzed before.'

    ELSE
      DO is=1,STR_LEN_A
        manager%label(is:is) = ''
      END DO
      is = MIN(LEN(TRIM(name)),STR_LEN_A)
      manager%label(1:is)  = name(1:is)
      manager%start_date   = date
      manager%pos          = 0
      manager%delta_time   = REAL(incr,dp)
      manager%init         = .TRUE.
      WRITE(message_text,*) &
           'Basic initialization of time manager <',&
           TRIM(manager%label),'> done.'

    END IF
    CALL message('',message_text)

  END SUBROUTINE TE_manager_init_days

  SUBROUTINE TE_manager_init_native (manager, name, date, incr)
    TYPE (time_manager), INTENT(inout) :: manager
    CHARACTER(len=*),    INTENT(in)    :: name
    TYPE (time_native),  INTENT(in)    :: date
    INTEGER,             INTENT(in)    :: incr

    TYPE (time_days)                   :: my_date

    CALL TC_convert (date,                    my_date)
    CALL TE_manager_init_days (manager, name, my_date, incr)

  END SUBROUTINE TE_manager_init_native


  SUBROUTINE TE_manager_init_intern (manager, name, date, incr)
    TYPE (time_manager), INTENT(inout) :: manager
    CHARACTER(len=*), INTENT(in)       :: name
    TYPE (time_intern), INTENT(in)     :: date
    INTEGER, INTENT(in)                :: incr

    TYPE (time_days)                   :: my_date

    CALL TC_convert (date,                    my_date)
    CALL TE_manager_init_days (manager, name, my_date, incr)

  END SUBROUTINE TE_manager_init_intern



  !------------------------------------------------------------
  ! reset the manager start date and correct the present time step
  ! the manager position in real time should not changed for consistency

  ! reset start date

  SUBROUTINE TE_manager_reinit_days (manager, new_date)
    TYPE (time_manager), INTENT(inout) :: manager
    TYPE (time_days),    INTENT(in)    :: new_date

    TYPE (time_days) :: current_date
    REAL(dp)         :: zsecs
    INTEGER          :: iday1, isec1, iday2, isec2, new_step

    IF (.NOT. manager%init) THEN
      CALL finish('mo_time_manager:TE_manager_reinit_days',&
           'Time manager was not initialized, no reinit possible.')

    ELSE IF (manager%freeze) THEN
      WRITE(message_text,*) 'Time manager <',&
           TRIM(manager%label),'> is frozen, reinit skiped.'
      CALL message('',message_text)

    ELSE

      ! get current date of manager
      CALL TE_manager_current_date (manager, current_date)

      ! get time since new start date
      CALL TC_get (current_date, iday1, isec1)
      CALL TC_get (new_date,     iday2, isec2)
      zsecs = REAL(idaylen,dp)*(iday1-iday2) + (isec1-isec2)

      ! get new manager step
      new_step = INT((zsecs+0.0001_dp)/manager%delta_time)

      IF (new_step < 0) &
           CALL message('Warning','New time manager position < 0')

      manager %pos        = new_step
      manager %start_date = new_date          ! correct the start date

    END IF

  END SUBROUTINE TE_manager_reinit_days


  SUBROUTINE TE_manager_reinit_native (manager, date)
    TYPE (time_manager), INTENT(inout) :: manager
    TYPE (time_native),  INTENT(in)    :: date

    TYPE (time_days)                   :: my_date

    CALL TC_convert (date,                my_date)
    CALL TE_manager_reinit_days (manager, my_date)

  END SUBROUTINE TE_manager_reinit_native


  SUBROUTINE TE_manager_reinit_intern (manager, date)
    TYPE (time_manager), INTENT(inout) :: manager
    TYPE (time_intern),  INTENT(in)    :: date

    TYPE (time_days)                   :: my_date

    CALL TC_convert (date,                my_date)
    CALL TE_manager_reinit_days (manager, my_date)

  END SUBROUTINE TE_manager_reinit_intern



  !------------------------------------------------------------
  ! reset increment between two steps

  SUBROUTINE TE_manager_reinit_incr (manager, incr)
    TYPE (time_manager), INTENT(inout) :: manager
    REAL(dp),            INTENT(in)    :: incr

    REAL(dp)         :: zsecs
    TYPE (time_days) :: current_date
    INTEGER          :: idays, isecs

    ! find a new start date
    ! the position of the manager will be conserved

    IF (manager%freeze) THEN
      WRITE(message_text,*) 'Time manager <',&
           TRIM(manager%label),'> is frozen, increment not changed.'
      CALL message('',message_text)

    ELSE

      ! distance from present point backward with new incr
      zsecs = manager%pos * incr
      idays = INT((zsecs+0.0001_dp)/REAL(idaylen,dp))
      isecs = INT(zsecs - REAL(idays,dp)*REAL(idaylen,dp)+0.0001_dp)
      idays = -idays
      isecs = -isecs

      ! get the present date
      CALL TE_manager_current_date (manager, current_date)

      ! calculates the new start date
      CALL add_date (idays, isecs,           current_date)

      manager %start_date = current_date
      manager %delta_time = REAL(incr,dp)     ! set the new increment

    END IF

  END SUBROUTINE TE_manager_reinit_incr

  !------------------------------------------------------------
  ! reset position

  SUBROUTINE TE_manager_reinit_step (manager, offset)
    TYPE (time_manager), INTENT(inout) :: manager
    INTEGER,             INTENT(in)    :: offset

    INTEGER :: new_step

    new_step = manager%pos + offset
    IF (new_step < 0) &
         CALL finish('mo_time_manager:TE_manager_reinit_step',&
         'negative steps invalid')

    manager %pos = new_step

  END SUBROUTINE TE_manager_reinit_step


  SUBROUTINE TE_manager_freeze (manager, lfreeze)
    TYPE (time_manager), INTENT(inout) :: manager
    LOGICAL,             INTENT(in)    :: lfreeze

    IF (.NOT. manager%init) THEN
      WRITE(message_text,*) &
           'Can not freeze state of time manager <',&
           TRIM(manager%label),&
           '>, missing initialization before.'
      CALL finish('mo_time_manager:TE_manager_freeze',message_text)

    ELSE IF (manager%freeze) THEN
      WRITE(message_text,*) 'Time manager <',&
           TRIM(manager%label),'> was frozen before.'
      CALL message('Warning',message_text)

    ELSE IF(lfreeze) THEN
      manager %freeze = .TRUE.
      WRITE(message_text,*) &
           'Initial state of time manager <',&
           TRIM(manager%label),'> is locked now!'
      CALL message('',message_text)

    END IF

    IF (.NOT.lfreeze) THEN
      manager %freeze = .FALSE.
      WRITE(message_text,*) 'The time manager <',&
           TRIM(manager%label),'> is unlocked now!'
      CALL message('',message_text)

    END IF

  END SUBROUTINE TE_manager_freeze

  !------------------------------------------------------------
  ! print out manager settings

  SUBROUTINE manager_print (manager, short)
    TYPE (time_manager), INTENT(in) :: manager
    LOGICAL, OPTIONAL,   INTENT(in) :: short

    TYPE (time_days)   :: current_date
    CHARACTER(len=256) :: my_message

    LOGICAL :: my_short

    my_short = .FALSE. ; IF (PRESENT(short)) my_short = short

    IF (manager%init) THEN

      IF (my_short) THEN         ! short message
        WRITE(message_text,*) &
             'State of >>',TRIM(manager%label),'<<'
        CALL message('',message_text)
        WRITE(message_text,*) &
             'Step: ',manager%pos,' dtime [s]: ',manager%delta_time
        CALL message('',message_text)

      ELSE                       ! long message
        WRITE(message_text,*) &
             'State of manager >>',TRIM(manager%label),'<<'
        CALL message('',message_text)
        WRITE(message_text,*) &
             ' Step counter : ', manager%pos
        CALL message('',message_text)
        WRITE(message_text,*) &
             ' Time step [s]: ', manager%delta_time
        CALL message('',message_text)

      END IF

      CALL print_date (manager%start_date,TC_PRN_NATIVE,mess=my_message)
      WRITE(message_text,*) 'Initial date: ',TRIM(my_message)
      CALL message('',message_text)

      CALL TE_manager_current_date (manager, current_date)
      CALL print_date (current_date,TC_PRN_NATIVE,mess=my_message)
      WRITE(message_text,*) 'Current date: ',TRIM(my_message)
      CALL message('',message_text)

    ELSE

      WRITE(message_text,*) 'Time manager >>',TRIM(manager%label), &
           '<< was not initialized, can not print settings.'
      CALL message('Warning',message_text)

    END IF
    CALL message('',' ')
    
  END SUBROUTINE manager_print


  !------------------------------------------------------------
  ! get date from manager

  SUBROUTINE TE_manager_any_date (manager, date, point)
    ! calculates date at a given timestep of the manager
    TYPE (time_manager), INTENT(in)  :: manager
    TYPE (time_days),    INTENT(out) :: date
    INTEGER,             INTENT(in)  :: point

    REAL(dp)     :: zsecs
    INTEGER      :: idays, isecs

    date  = manager %start_date      ! load start date

    ! calculates new offset
    IF (point < 0) THEN
      zsecs = REAL(point,dp)*manager%delta_time - 0.0001_dp
    ELSE
      zsecs = REAL(point,dp)*manager%delta_time + 0.0001_dp
    END IF
    idays = INT(zsecs/REAL(idaylen,dp))

    IF (zsecs < 0.0_dp) THEN
      isecs = INT(zsecs - REAL(idays,dp)*REAL(idaylen,dp) - 0.0001_dp)
    ELSE
      isecs = INT(zsecs - REAL(idays,dp)*REAL(idaylen,dp) + 0.0001_dp)
    END IF

    ! final calculation of date at point
    CALL add_date (idays, isecs, date)

  END SUBROUTINE TE_manager_any_date


  SUBROUTINE TE_manager_current_date (manager, date)
    TYPE (time_manager), INTENT(in)  :: manager
    TYPE (time_days),    INTENT(out) :: date

    CALL TE_manager_any_date (manager, date, manager%pos)

  END SUBROUTINE TE_manager_current_date


  !------------------------------------------------------------
  ! get nearest manager step for special date

  SUBROUTINE TE_manager_step_at_date (manager, current_date, steps, lfit)
    ! get the time step (pointer) for a known date/time
    TYPE (time_manager), INTENT(in)  :: manager
    TYPE (time_days),    INTENT(in)  :: current_date
    INTEGER,             INTENT(out) :: steps
    LOGICAL,             INTENT(out) :: lfit

    TYPE (time_days) :: check_date
    INTEGER          :: iday1, isec1, iday2, isec2, iday, isecond
    REAL(dp)         :: zsecs

    ! calculate time difference in seconds
    CALL TC_get (current_date,       iday1, isec1)
    CALL TC_get (manager%start_date, iday2, isec2)

    ! calculates the difference between the dates in seconds
    zsecs = REAL(iday1-iday2,dp) * idaylen +  REAL(isec1-isec2,dp)
    IF (zsecs < 0.0_dp) THEN
      zsecs = zsecs-0.001_dp
    ELSE
      zsecs = zsecs+0.001_dp
    END IF
    steps = INT(zsecs / manager%delta_time)

    ! recalculate the present day again using the number of steps found
    zsecs = REAL(manager%delta_time*steps,dp)
    IF (zsecs < 0) THEN
      zsecs   = zsecs - 0.001_dp
    ELSE
      zsecs   = zsecs + 0.001_dp
    END IF
    iday    = INT(zsecs/idaylen)
    isecond = MOD(INT(zsecs),idaylen)

    check_date = manager %start_date
    CALL add_date (iday, isecond, check_date)

    IF (check_date == current_date) THEN
      lfit = .TRUE.
    ELSE
      lfit = .FALSE.
      CALL message('TE_manager_step_at_date',&
           'date does not fit into time grid')
    END IF

  END SUBROUTINE TE_manager_step_at_date


!BOP
  ! !IROUTINE: TE_manager_dtime
  ! !INTERFACE:
  SUBROUTINE TE_manager_dtime (manager,delta_time)

    ! !DESCRIPTION:
    ! return the manager increment (time step)
    ! !INPUT/OUTPUT PARAMETERS:
    TYPE (time_manager), INTENT(in)  :: manager
    REAL(dp),            INTENT(out) :: delta_time
!EOP
    IF (manager%init) THEN
      delta_time = manager%delta_time
    ELSE
      CALL manager_print(manager)
      CALL finish('mo_time_manager:TE_manager_dtime','manager not initialized')
    END IF

  END SUBROUTINE TE_manager_dtime

!BOP
  ! !IROUTINE: TE_manager_pos
  ! !INTERFACE:
  SUBROUTINE TE_manager_pos (manager,position)

    ! !DESCRIPTION:
    ! return the manager position
    ! !INPUT/OUTPUT PARAMETERS:
    TYPE (time_manager), INTENT(in)  :: manager
    INTEGER            , INTENT(out) :: position
!EOP
    IF (manager%init) THEN
      position = manager %pos
    ELSE
      CALL manager_print(manager)
      CALL finish('mo_time_manager:TE_manager_dtime','manager not initialized')
    END IF

  END SUBROUTINE TE_manager_pos

!BOP
  ! !IROUTINE: rewind_manager
  ! !INTERFACE:

  SUBROUTINE rewind_manager (manager, new_date, fix_start_date)

    ! !DESCRIPTION:
    ! reset the present date of the time manager
    ! !INPUT/OUTPUT PARAMETERS:
    TYPE (time_manager), INTENT(inout) :: manager         ! the time manager
    TYPE (time_days),    INTENT(in)    :: new_date        ! set position to that date
    LOGICAL,             INTENT(in)    :: fix_start_date  ! TRUE will change the time step
!EOP
    INTEGER          :: isecs, idays, new_step, old_step, offset
    TYPE (time_days) :: wrk_day
    LOGICAL          :: lfit

    IF (manager%freeze) THEN

      WRITE(message_text,*) 'Operation not allowed on time manager',&
           TRIM(manager%label),': unfreeze first!'
      CALL finish('mo_time_manager:rewind_manager',message_text)

    ELSE IF (fix_start_date) THEN
      ! *********** reset manager without changing the start date

      IF (new_date < manager%start_date) &
           CALL finish('mo_time_manager:rewind_manager','new date before start date')

      CALL manager_state(manager,old_step)                            ! get present position
      CALL TE_manager_step_at_date(manager, new_date, new_step, lfit) ! get new position

      IF (.NOT. lfit) &
           CALL finish('mo_time_manager:rewind_manager',&
           'new date do not fit to a time step')

      offset = new_step - old_step
      CALL manager_init(manager,offset)                               ! set to new position

    ELSE
      ! ************ reset the manager with changing the start date

      CALL manager_state (manager, wrk_day)          ! get current_date of manager
      CALL TC_get                 (wrk_day, idays, isecs)
      wrk_day = new_date
      CALL add_date (-idays, -isecs, wrk_day)        ! wrk_day = new_date - current_date
      CALL TC_get                   (wrk_day, idays, isecs)
      CALL add_date (idays,isecs,manager%start_date) ! correct start date of manager

      ! the time step is not changed

    END IF

  END SUBROUTINE rewind_manager

END MODULE mo_time_manager
