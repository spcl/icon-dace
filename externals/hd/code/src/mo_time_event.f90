! mo_time_event.f90 - Set, control and manipulate time events
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_time_event

  !+
  !
  ! mo_time_event [module]
  !
  !   implementation of EVENTs, their control and manipulations with
  !
  ! Version: E5/R1.07+ 03-February-2003
  !
  ! Authors:
  !   I.Kirchner,  MPI Hamburg, Nov-99
  !   I.Kirchner,  MPI, May 2000, revison
  !   L.Kornblueh, MPI, August 2000, minor updates
  !   I.Kirchner,  MPI, October 2000, update
  !   I.Kirchner,  MPI, March/September 2001, revision
  !   I.Kirchner,  FUB, Feb/2003, revision/code review
  !
  ! external modules
  !    mo_kind
  !    mo_exception
  !    mo_time_conversion
  !    mo_time_base
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.
  !
  !-
  USE mo_kind,            ONLY: dp
  USE mo_exception,       ONLY: finish, message
  USE mo_time_conversion, ONLY: TC_convert, TC_set, TC_get, &
                                add_date, print_date, &
                                time_days, time_native, time_intern, &
                                OPERATOR(==), OPERATOR(<), OPERATOR(>)
  USE mo_time_base,       ONLY: IDAYLEN, get_calendar_type, JULIAN, CYL360, &
                                Get_JulianMonLen, Get_Ly360MonLen

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: STR_LEN_A = 20
  INTEGER, PARAMETER :: STR_LEN_B = 40
  INTEGER, PARAMETER :: STR_LEN_C = 256

  CHARACTER(len=STR_LEN_C) :: m_text
  CHARACTER(len=STR_LEN_C) :: date_text

  !+
  ! **************** parameters ------------------------------------------------
  !
  ! TIME_INC_*     predefined counter units
  !
  CHARACTER(len=*), PUBLIC, PARAMETER :: &
       TIME_INC_SECONDS = 'seconds'  ,&!
       TIME_INC_MINUTES = 'minutes'  ,&!
       TIME_INC_HOURS   = 'hours'    ,&!
       TIME_INC_DAYS    = 'days'     ,&!
       TIME_INC_MONTHS  = 'months'   ,&!
       TIME_INC_YEARS   = 'years'      !

  ! TRIG_*  type of trigger adjustment of an event
  !
  CHARACTER(len=*), PUBLIC, PARAMETER :: &
       TRIG_FIRST  = 'first'  ,&! trigger in first step of counter unit
       TRIG_LAST   = 'last'   ,&! trigger in last step of counter unit
       TRIG_EXACT  = 'exact'  ,&! trigger without adjustment in side the unit
       TRIG_NONE   = 'off'     ! event trigger non active

  ! **************** structures ------------------------------------------------
  !
  TYPE, PUBLIC :: io_time_event            ! external given event properties
    INTEGER                  :: counter    = 0                ! No. of steps in given unit
    CHARACTER(len=STR_LEN_A) :: unit       = TIME_INC_SECONDS ! counter unit type
    CHARACTER(len=STR_LEN_A) :: adjustment = TRIG_EXACT       ! adjustment in side the unit
    INTEGER                  :: offset     = 0            ! offset to initial date in seconds
  END TYPE io_time_event

  TYPE, PUBLIC :: time_event           !   hold all event relevant informations
    PRIVATE
    LOGICAL                  :: init   = .FALSE.    ! event state for access control
    CHARACTER(len=STR_LEN_B) :: label  = ''         ! short description of the event
    INTEGER                  :: count  = 1          ! increment between the action
    CHARACTER(len=STR_LEN_A) :: unit   = TIME_INC_SECONDS  ! name of the basic units
    INTEGER                  :: offset = 0          ! offset (seconds)
    CHARACTER(len=STR_LEN_A) :: adjust = TRIG_EXACT ! type of triggering
    REAL(dp)                 :: delta      = 2.0_dp ! adjustment interval in seconds
    REAL(dp)                 :: half_delta = 1.0_dp
    LOGICAL                  :: active = .FALSE.    ! event active true=on false=off

    TYPE (time_days)         ::  initial_date     ! initial data for event
    TYPE (time_days)         ::    cycle_date     ! without offset
    TYPE (time_days)         :: previous_trigger  ! previous trigger date
    TYPE (time_days)         ::  current_trigger  ! current trigger date
    TYPE (time_days)         ::     next_trigger  ! next trigger date

    LOGICAL                  :: flag1 = .FALSE.   ! flags for special purpose
    LOGICAL                  :: flag2 = .FALSE.
    LOGICAL                  :: flag3 = .FALSE.
    LOGICAL                  :: flag4 = .FALSE.
    LOGICAL                  :: flag5 = .FALSE.

  END TYPE time_event

  ! *************** INTERFACE SUBROUTINES AND FUNCTIONS ------------------------
  !

  PUBLIC :: event_init        !   initialise parts of an event
  INTERFACE event_init
    MODULE PROCEDURE TE_event_init
                              ! (IO:event,I:name,I:count,I:unit,I:adj[,I:delta])
    MODULE PROCEDURE TE_event_set_day        ! (IO:event,I:time_days)
    MODULE PROCEDURE TE_event_set_day_intern ! (IO:event,I:time_intern)
    MODULE PROCEDURE TE_event_set_day_native ! (IO:event,I:time_native)
    MODULE PROCEDURE TE_event_toggle         ! (IO:event,I:time_days[,I:lshift])
  END INTERFACE

  PUBLIC :: event_reinit     ! reset initial trigger dates
  INTERFACE event_reinit
    MODULE PROCEDURE TE_event_reset_day        ! (IO:event,I:time_days)
    MODULE PROCEDURE TE_event_reset_day_intern ! (IO:event,I:time_intern)
    MODULE PROCEDURE TE_event_reset_day_native ! (IO:event,I:time_native)
  END INTERFACE

  PUBLIC :: event_print                     ! (I:event[,I:format]) print event
  INTERFACE event_print
    MODULE PROCEDURE TE_print_event
  END INTERFACE

  PUBLIC :: TE_print_event_name             ! (I:event) returns event name 

  PUBLIC :: event_state                     ! get the state of event
  INTERFACE event_state
                                            !  get interval
    MODULE PROCEDURE I_event_steps          ! int<-(I:event,I:delta) in steps
    MODULE PROCEDURE I_event_seconds        ! int<-(I:event)  in seconds
    MODULE PROCEDURE I_event_seconds_next   ! int<-(I:event,I:true) future (s)
                                            !  evaluate next trigger, check date
    MODULE PROCEDURE L_event_trigger_days   ! log<-(I:event,I:time_days)
    MODULE PROCEDURE L_event_trigger_intern ! log<-(I:event,I:time_intern)
    MODULE PROCEDURE L_event_trigger_native ! log<-(I:event,I:time_native)
  END INTERFACE

  PUBLIC  :: event_eval         ! int<-(I:event,I:delta) fit delta into interval
  PUBLIC  :: event_next_date    ! (I:event,O:next_trigger)
  PUBLIC  :: event_current_date ! (I:event,O:current_trigger)

  PRIVATE :: TE_get_next_trigger ! (I:event) find next trigger

  !-
CONTAINS

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE TE_event_init (event, name, count, unit, adjust, delta, offset) !**

    ! initialize elements of an event structure
    !-
    TYPE (time_event), INTENT(inout) :: event   ! new event
    CHARACTER(len=*),  INTENT(in)    :: name    ! name of event
    INTEGER,           INTENT(inout) :: count   ! no. of units between events
    CHARACTER(len=*),  INTENT(inout) :: unit    ! unit
    CHARACTER(len=*),  INTENT(in)    :: adjust  ! type of adjustment
    REAL(dp), OPTIONAL,INTENT(in)    :: delta   ! adjustment interval [seconds]
    INTEGER, OPTIONAL, INTENT(in)    :: offset  ! offset relativ to initial date

    INTEGER        :: is, isec1, isec2
    REAL(dp)       :: smallest

    IF (event%init) THEN
      CALL message('','Event was initialized, nothing changed')

    ELSE
      DO is=1,STR_LEN_B
        event%label(is:is) = ''
      END DO
      is = MIN(LEN(TRIM(name)),STR_LEN_B)
      event%label(1:is) = name(1:is)

      ! **** evaluate offset ---------------------------------------------------

      event%offset = 0
      IF (PRESENT(offset)) event%offset = offset
      IF (event%offset /= 0) THEN

        IF (adjust /= TRIG_EXACT) THEN

          SELECT CASE(unit)
          CASE(TIME_INC_MONTHS,TIME_INC_YEARS)
            CALL finish('mo_time_event:TE_event_init','combination of offset/unit invalid')
          END SELECT

          ! correct unit, depenend on the offset
          IF (MOD(event%offset,24*3600) == 0) THEN ! it's fine

          ELSE IF (MOD(event%offset,3600) == 0) THEN
            SELECT CASE(unit)
            CASE(TIME_INC_DAYS);  unit = TIME_INC_HOURS;   count = count*24
            END SELECT

          ELSE IF (MOD(event%offset,60) == 0) THEN
            SELECT CASE(unit)
            CASE(TIME_INC_DAYS);  unit = TIME_INC_MINUTES; count = count*24*60
            CASE(TIME_INC_HOURS); unit = TIME_INC_MINUTES; count = count*60
            END SELECT
            
          ELSE 
            SELECT CASE(unit)
            CASE(TIME_INC_DAYS);    count = count*24*3600
            CASE(TIME_INC_HOURS);   count = count*3600
            CASE(TIME_INC_MINUTES); count = count*60
            END SELECT
            unit  = TIME_INC_SECONDS

          END IF

        END IF

      END IF
      event%count = count

      ! **** define the smallest delta needed for adjustment -------------------

      SELECT CASE(unit)
      CASE(TIME_INC_SECONDS);  smallest = 1.5_dp
      CASE(TIME_INC_MINUTES);  smallest = 60.0_dp
      CASE(TIME_INC_HOURS);    smallest = 3600.0_dp
      CASE(TIME_INC_DAYS);     smallest = 24.0_dp*3600.0_dp
      CASE(TIME_INC_MONTHS);   smallest = 28.0_dp*24.0_dp*3600.0_dp
      CASE(TIME_INC_YEARS);    smallest = 360.0_dp*28.0_dp*24.0_dp*3600.0_dp
      CASE default
        m_text = 'Counter unit unknown ::' // TRIM(unit)
        CALL finish('mo_time_event:TE_event_init',m_text)
      END SELECT
      smallest   = smallest * REAL(count,dp) - 1.0_dp
      event%unit = unit

      ! *** evaluate epsilon interval and adjustment ---------------------------

      SELECT CASE(adjust)
      CASE (TRIG_FIRST, TRIG_LAST,TRIG_EXACT,TRIG_NONE)
        IF (PRESENT(delta)) THEN
          IF (0 < delta .AND. delta <= smallest) THEN
            event%delta = delta
          ELSE
            event%delta = smallest
            WRITE(m_text,*) 'Preset adjustment interval ::',smallest
            CALL message('',m_text)
          END IF

        ELSE
          event%delta = smallest
          WRITE(m_text,*) 'Preset adjustment interval ::',smallest
          CALL message('',m_text)

        END IF

      CASE default
        m_text = 'Event adjustment not defined ::' // TRIM(adjust)
        CALL finish('mo_time_event:TE_event_init',m_text)

      END SELECT
      
      SELECT CASE(adjust)
      CASE (TRIG_NONE)
        event %adjust = TRIG_EXACT
        event %active = .FALSE.
      CASE default
        event %adjust = TRIM(adjust)
        event %active = .TRUE.
      END SELECT

      ! **** calculates half delta in full seconds -----------------------------

      IF (event%delta < 1.0_dp) &
           CALL finish('mo_time_event:TC_event_init','adjustment interval too small')

      isec1 = INT(0.5_dp * event%delta)
      isec2 = INT(event%delta) - isec1
      event%half_delta = REAL(isec1,dp)
      IF (isec1 < isec2 ) event%half_delta = REAL(isec2,dp)

      ! **** set initial date to calendar start --------------------------------

      CALL TC_set (0,0,event %initial_date)
      CALL TC_set (0,0,event %cycle_date)
      CALL TC_set (0,0,event %current_trigger)
      CALL TC_set (0,0,event %previous_trigger)
      CALL TC_set (0,0,event %next_trigger)

      event %init    = .TRUE.
    END IF

  END SUBROUTINE TE_event_init


  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE TE_event_toggle (event, adjust)  !*********************************

    ! reset the adjustment type
    ! used for activate/deactive an event
    !-
    TYPE (time_event), INTENT(inout) :: event   ! new event
    CHARACTER(len=*),  INTENT(in)    :: adjust  ! type of adjustment

    SELECT CASE(adjust)
    CASE (TRIG_FIRST, TRIG_LAST, TRIG_EXACT)
      m_text = 'Reset event adjustment for ... ' // TRIM(event%label)
      event%adjust = TRIM(adjust)
      event%active = .TRUE.

    CASE (TRIG_NONE)
      m_text = 'Deactivate event ... ' // TRIM(event%label)
      event%active = .FALSE.

    CASE default
       m_text = 'Event adjustment not defined <' // TRIM(adjust) // '>'
       CALL finish('mo_time_event:TE_event_toggle',m_text)

    END SELECT
    CALL message('',m_text)


  END SUBROUTINE TE_event_toggle


  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE TE_event_set_day (event, date, lshift) !***************************

    ! set first trigger date
    !-
    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_days),  INTENT(in)    :: date
    LOGICAL, OPTIONAL, INTENT(in)    :: lshift

    LOGICAL :: llshift

    llshift = .FALSE. ; IF (PRESENT(lshift)) llshift = lshift

    IF (.NOT. event%init) THEN
      m_text = 'Initialize event <' // TRIM(event%label) // '>first'
      CALL finish('mo_time_event:mo_time_event:TE_event_set_day',m_text)
      
    ELSE

      IF (llshift) THEN
        event% previous_trigger = event% current_trigger
        
      ELSE
        event%  initial_date    = date
        event%    cycle_date    = date
        event% previous_trigger = date
      END IF
      event%   current_trigger  = date

      CALL TE_get_next_trigger(event)

    END IF

  END SUBROUTINE TE_event_set_day

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE TE_event_set_day_intern (event, date, lshift)  !*******************
    !-
    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_intern), INTENT(in)   :: date
    LOGICAL, OPTIONAL                :: lshift

    TYPE (time_days) :: my_date

    CALL TC_convert(date, my_date)
    IF (PRESENT(lshift)) THEN
      CALL TE_event_set_day(event, my_date, lshift)
    ELSE
      CALL TE_event_set_day(event, my_date)
    END IF

  END SUBROUTINE TE_event_set_day_intern


  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE TE_event_set_day_native (event, date, lshift) !********************
    !-
    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_native), INTENT(in)   :: date
    LOGICAL, OPTIONAL                :: lshift

    TYPE (time_days) :: my_date

    CALL TC_convert(date, my_date)
    IF (PRESENT(lshift)) THEN
      CALL TE_event_set_day(event, my_date, lshift)
    ELSE
      CALL TE_event_set_day(event, my_date)
    END IF

  END SUBROUTINE TE_event_set_day_native


  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE TE_event_reset_day (event, date)  !********************************

    ! reset first trigger date without evaluation of next trigger
    !-
    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_days), INTENT(in)     :: date
 
    IF (event%next_trigger < date) THEN
      m_text = 'Next trigger date of <' // TRIM(event%label) // &
           '> in the past. Reset not possible.'
      CALL finish('mo_time_event:mo_time_event:TE_event_reset_day',m_text)
    END IF
 
    event% previous_trigger = date
    event%  current_trigger = date
 
  END SUBROUTINE TE_event_reset_day

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE TE_event_reset_day_intern (event, date)  !*************************
    !-
    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_intern), INTENT(in)   :: date
 
    TYPE (time_days) :: my_date
 
    CALL TC_convert (date,          my_date)
    CALL TE_event_reset_day (event, my_date)
 
  END SUBROUTINE TE_event_reset_day_intern
 
  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE TE_event_reset_day_native (event, date) !**************************
    !-
    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_native), INTENT(in)   :: date
 
    TYPE (time_days) :: my_date
 
    CALL TC_convert (date,          my_date)
    CALL TE_event_reset_day (event, my_date)
 
  END SUBROUTINE TE_event_reset_day_native
 
  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE TE_print_event (event, short)  !***********************************

    ! print event contents
    !-
    TYPE (time_event), INTENT(in) :: event
    LOGICAL, OPTIONAL, INTENT(in) :: short

    TYPE (time_native) :: my_date

    IF (event%init) THEN

      IF (PRESENT(short)) THEN
        IF(short) THEN
          WRITE(m_text,*) 'Event <',TRIM(event%label),&
              '> : interval ',event%count,' ',TRIM(event%unit),&
              ' : adjustment ',TRIM(event%adjust),' : offset[sec] ',event%offset

        ELSE
          m_text = 'trigger event >>' // TRIM(event%label) // '<<'

        END IF
        CALL message('',m_text)

      ELSE
        IF (.NOT. event%active) THEN
          WRITE(m_text,*) 'Event <',TRIM(event%label),'> ... not active'

        ELSE
          CALL message('',' ')
          WRITE(m_text,*) &
               'State of event <',TRIM(event%label),'> ... initialized'
          CALL message('',m_text)

          SELECT CASE(event%unit)
          CASE(TIME_INC_SECONDS)
            WRITE(m_text,*) ' trigger each ',event%count,' seconds'

          CASE(TIME_INC_MINUTES)
            WRITE(m_text,*) ' trigger each ',event%count,' minutes'

          CASE(TIME_INC_HOURS)
            WRITE(m_text,*) ' trigger each ',event%count,' hours'

          CASE(TIME_INC_DAYS)
            WRITE(m_text,*) ' trigger each ',event%count,' days'

          CASE(TIME_INC_MONTHS)
            WRITE(m_text,*) ' trigger each ',event%count,' months'

          CASE(TIME_INC_YEARS)
            WRITE(m_text,*) ' trigger each ',event%count,' years'

          CASE default
            WRITE(m_text,*) 'Counter unit unknown ::',event%unit
            CALL finish('mo_time_event:TE_print_event',m_text)
          END SELECT
          CALL message('',m_text)
             
          CALL TC_convert(event% previous_trigger,my_date)
          CALL print_date(my_date,mess=date_text)
          m_text = '  last event trigger date is ...    ' // TRIM(date_text)
          CALL message('',m_text)

          CALL TC_convert(event% current_trigger,my_date)
          CALL print_date(my_date,mess=date_text)
          IF (event% current_trigger == event% previous_trigger) THEN
            m_text = '  initial event date is ...         ' // TRIM(date_text)
            CALL message('',m_text)
          ELSE
            WRITE(m_text,*) '  time between two triggers: ',&
                 I_event_seconds(event),' seconds'
            CALL message('',m_text)
            m_text = '  present trigger date is ...       ' // TRIM(date_text)
            CALL message('',m_text)
          END IF

          CALL TC_convert(event% next_trigger,my_date)
          CALL print_date(my_date,mess=date_text)
          WRITE(m_text,*)    '  next trigger date is ...         ',&
               TRIM(date_text),' (offset of ',event%offset,' seconds included)'
          CALL message('',m_text)

          SELECT CASE(event%adjust)
          CASE(TRIG_FIRST)
            WRITE(m_text,*) '  adjustment at beginning of unit,',&
                 ' interval: ( trigger : trigger + ', NINT(event%delta),'s )'
                
          CASE(TRIG_LAST)
            WRITE(m_text,*) '  adjustment at end of unit,',&
                 ' interval: ( trigger - ', NINT(event%delta),'s : trigger )'

          CASE(TRIG_EXACT)
            WRITE(m_text,*) '  no adjustment,',&
                 ' interval: ( trigger - ',INT(event%half_delta),&
                 's : trigger + ',INT(event%half_delta)-1,'s )'

          END SELECT

        END IF
        CALL message('',m_text)
        CALL message('',' ')

      END IF

    ELSE
      CALL message('Warning','no printout of uninitialized event')

    END IF


  END SUBROUTINE TE_print_event

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION TE_print_event_name (ev) RESULT (str)  !*****************************

    ! return event name
    !-
    TYPE(time_event), INTENT(in) :: ev
    CHARACTER(len=STR_LEN_B)     :: str

    str = TRIM(ev%label)

  END FUNCTION TE_print_event_name


  ! ----------------------------------------------------------------------------
  !+
  FUNCTION I_event_steps (event, delta) RESULT (ix)  !**************************

    ! get out interval between two event trigger points in steps
    !-
    TYPE (time_event), INTENT(in) :: event
    REAL(dp),          INTENT(in) :: delta   ! seconds of special unit
    INTEGER                       :: ix

    INTEGER :: isecs

    isecs = I_event_seconds(event)
    ix = INT(REAL(isecs+0.0001_dp,dp)/delta) ! convert into steps

  END FUNCTION I_event_steps

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION I_event_seconds (event) RESULT (ix) !********************************

    ! get out interval between two event trigges in seconds
    !-
    TYPE (time_event), INTENT(in) :: event
    INTEGER                       :: ix

    INTEGER           :: iday, isec
    TYPE (time_days)  :: my_date

    ! get difference between last and present trigger date
    CALL TC_get (event%previous_trigger, iday, isec)
    iday = - iday
    isec = - isec
    my_date = event%current_trigger
    CALL add_date (iday, isec, my_date)

    ! check result
    CALL TC_get(my_date, iday, isec)
    IF (iday < 0 .OR. isec < 0) &
         CALL finish('mo_time_event:I_event_seconds','Event interval < 0')

    ! convert into seconds
    ix = iday*idaylen+isec

    ! at initial time the interval is zero
    ! this (may be) is important for accumulation

  END FUNCTION I_event_seconds

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION I_event_seconds_next (event, lnext) RESULT (ix)  !*******************

    ! return distance between present and next trigger
    !-
    TYPE (time_event), INTENT(in) :: event
    LOGICAL,           INTENT(in) :: lnext
    INTEGER                       :: ix

    INTEGER           :: iday, isec
    TYPE (time_days)  :: my_date

    ix = 0
    IF (lnext) THEN

      CALL get_difference (event%current_trigger, event%next_trigger, my_date)

      CALL TC_get(my_date, iday, isec)
      IF (iday < 0 .OR. isec < 0) &
           CALL finish('mo_time_event:I_event_seconds_next','Event interval < 0')

      ! convert into seconds
      ix = INT(iday*idaylen + isec + 0.0001_dp)

    END IF

  END FUNCTION I_event_seconds_next

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE get_difference (date1, date2, date3)
    TYPE (time_days), INTENT(in)  :: date1, date2
    TYPE (time_days), INTENT(out) :: date3
    !-
    INTEGER           :: iday, isec

    CALL TC_get (date1, iday,  isec)
                                     date3 = date2
    CALL add_date     (-iday, -isec, date3)

  END SUBROUTINE get_difference

  ! ----------------------------------------------------------------------------
  !+
  LOGICAL FUNCTION L_event_trigger_intern (event, date)  !**********************

    ! ask an event for a trigger
    !-
    TYPE (time_event),  INTENT(inout) :: event
    TYPE (time_intern), INTENT(in)    :: date

    TYPE (time_days) :: my_date

    CALL TC_convert (date,                                my_date)
    L_event_trigger_intern = L_event_trigger_days (event, my_date)

  END FUNCTION L_event_trigger_intern

  ! ----------------------------------------------------------------------------
  !+
  LOGICAL FUNCTION L_event_trigger_native (event, date)  !**********************
    !-
    TYPE (time_event),  INTENT(inout) :: event
    TYPE (time_native), INTENT(in)    :: date

    TYPE (time_days) :: my_date

    CALL TC_convert (date,                                my_date)
    L_event_trigger_native = L_event_trigger_days (event, my_date)

  END FUNCTION L_event_trigger_native


  ! ----------------------------------------------------------------------------
  !+
  LOGICAL FUNCTION L_event_trigger_days (event, date)  !************************
    !-
    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_days),  INTENT(in)    :: date

    REAL(dp)         :: zsec  ! length of event increment in seconds
    INTEGER          :: iday1, isec1, iday2, isec2
    TYPE (time_days) :: adj_start, adj_stop

    L_event_trigger_days = .FALSE.

    IF (.NOT. event%active) RETURN

    ! *** check if the date will fit into the trigger interval -----------------

    zsec = event%delta
    IF (event%adjust == TRIG_EXACT) zsec = event%half_delta
    iday1 = INT((zsec+0.0001_dp)/REAL(idaylen,dp))
    isec1 = INT(zsec - REAL(iday1,dp)*REAL(idaylen,dp) + 0.0001_dp)

    IF (event%adjust == TRIG_EXACT) THEN
       zsec  = event%half_delta-1.0_dp
       iday2 = INT((zsec+0.0001_dp)/REAL(idaylen),dp)
       isec2 = INT(zsec - REAL(iday2,dp)*REAL(idaylen,dp) + 0.0001_dp)
    END IF

    ! *** detect first or last second in given unit ----------------------------

    adj_start = event% next_trigger
    adj_stop  = event% next_trigger

    SELECT CASE(event%adjust)
    CASE(TRIG_FIRST);   CALL add_date( iday1, isec1, adj_stop)
    CASE(TRIG_LAST);    CALL add_date(-iday1,-isec1, adj_start)
    CASE(TRIG_EXACT)
      CALL add_date( iday1, isec1, adj_stop)
      CALL add_date(-iday2,-isec2, adj_start)
    END SELECT

    ! *** check the position of the present date -------------------------------

    IF (adj_stop < date) THEN
      m_text = 'Warning: event trigger not longer in future <'&
           // TRIM(event%label) // '>'
      CALL message('L_event_trigger_days',m_text)

    ELSE                                    ! adj_start <= my_date <= adj_stop
      L_event_trigger_days = &
            (adj_start==date) .OR. (date==adj_stop) .OR. &
           ((adj_start <date).AND. (date <adj_stop))

      IF (L_event_trigger_days) CALL TE_event_set_day (event, date, .TRUE.)
    END IF

  END FUNCTION L_event_trigger_days


  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE TE_get_next_trigger (event)  !*************************************

    ! calculates next trigger from present trigger
    !-
    TYPE (time_event), INTENT(inout) :: event

    TYPE (time_native)               :: my_date
    LOGICAL                          :: first_call

    ! calculates the end of the next time interval
    !
    CALL TC_convert (event%cycle_date,           my_date)
    CALL get_next_date (event%unit, event%count, my_date)

    ! adjust the interval date
    !
    first_call = (event%previous_trigger == event%current_trigger)
    CALL adjust_date(first_call, event%adjust, event%unit, event%delta, my_date)

    ! store the adjusted interval end
    !
    CALL TC_convert (my_date, event%cycle_date)

    ! perform the offset calculation
    !
    CALL add_date (0, event%offset, my_date)
    CALL TC_convert                (my_date, event%next_trigger)

    IF (.NOT. event%current_trigger < event%next_trigger) THEN
      m_text = 'Event <' // TRIM(event%label) // '> [current date > next date]'
      CALL message('',m_text)
      CALL print_date (event% current_trigger,'native')
      CALL print_date (event%    next_trigger,'native')
      CALL finish('mo_time_event:TE_get_next_trigger',&
           'can not evaluate next trigger right')
    END IF

  CONTAINS

    SUBROUTINE get_next_date (unit, count, my_date)
      CHARACTER (len=*),  INTENT(in)    :: unit
      INTEGER,            INTENT(in)    :: count
      TYPE (time_native), INTENT(inout) :: my_date 
      !-
      INTEGER        :: yr, mo, dy, hr, mn, se, iday, isec
      REAL (kind=dp) :: zsec

      CALL TC_get(my_date, yr, mo, dy, hr, mn, se)

      SELECT CASE(unit)
      CASE(TIME_INC_SECONDS, TIME_INC_MINUTES, TIME_INC_HOURS, TIME_INC_DAYS)

        ! convert seconds into days and seconds

        SELECT CASE(unit)       ! calculate increment 
        CASE(TIME_INC_SECONDS);  zsec =           1._dp
        CASE(TIME_INC_MINUTES);  zsec =          60._dp
        CASE(TIME_INC_HOURS);    zsec =        3600._dp
        CASE(TIME_INC_DAYS);     zsec = 24._dp*3600._dp
        END SELECT

        zsec = zsec * REAL(count,dp)
        iday = INT((zsec+0.0001_dp)/REAL(IDAYLEN,dp))
        isec = INT(zsec - REAL(iday*IDAYLEN,dp) + 0.0001_dp)

        ! add  day and time to last adjustment trigger date
        CALL add_date(iday, isec, my_date)

      CASE(TIME_INC_MONTHS)
        mo = mo + count
        yr = yr + INT((mo-0.5_dp)/12)
        mo = MOD(mo,12); IF (mo == 0) mo = 12
        SELECT CASE (get_calendar_type())
        CASE (JULIAN)
          dy = MIN(dy,Get_JulianMonLen(yr,mo))
        CASE (CYL360)
          dy = MIN(dy,Get_Ly360MonLen())
        END SELECT
        CALL TC_set(yr, mo, dy, hr, mn, se, my_date)

      CASE(TIME_INC_YEARS)
        yr = yr + count
        SELECT CASE (get_calendar_type())
        CASE (JULIAN)
          dy = MIN(dy,Get_JulianMonLen(yr,mo))
        CASE (CYL360)
          dy = MIN(dy,Get_Ly360MonLen())
        END SELECT
        CALL TC_set(yr, mo, dy, hr, mn, se, my_date)

      END SELECT

    END SUBROUTINE get_next_date

    !+
    SUBROUTINE adjust_date (first_call, adjust, unit, delta, my_date)
      LOGICAL,            INTENT(in)    :: first_call
      CHARACTER (len=*),  INTENT(in)    :: adjust, unit
      REAL (kind=dp),     INTENT(in)    :: delta
      TYPE (time_native), INTENT(inout) :: my_date 
      !-
      INTEGER :: yr, mo, dy, hr, mn, se, iday, isec
      
      SELECT CASE(adjust)
      CASE(TRIG_FIRST)

        CALL TC_get (my_date, yr, mo, dy, hr, mn, se)

        ! adjust next trigger date to first delta environment in unit
        ! interval = <first second, first second + delta>

        SELECT CASE(unit)
        CASE(TIME_INC_MINUTES);                               se = 0
        CASE(TIME_INC_HOURS);                         mn = 0; se = 0
        CASE(TIME_INC_DAYS);                  hr = 0; mn = 0; se = 0
        CASE(TIME_INC_MONTHS);        dy = 1; hr = 0; mn = 0; se = 0
        CASE(TIME_INC_YEARS); mo = 1; dy = 1; hr = 0; mn = 0; se = 0
        END SELECT
        CALL TC_set          (yr, mo, dy, hr, mn, se, my_date)

      CASE(TRIG_LAST)

        IF (first_call) THEN ! skip back (one fraction - delta_adjustment)

          iday = 0
          isec = INT(delta)
          CALL add_date (iday, isec, my_date)

          CALL TC_get (my_date, yr, mo, dy, hr, mn, se)

          SELECT CASE(unit)
          CASE(TIME_INC_MINUTES,TIME_INC_HOURS,TIME_INC_DAYS)

            SELECT CASE(unit)
            CASE(TIME_INC_MINUTES); isec =   -60
            CASE(TIME_INC_HOURS);   isec = -3600
            CASE(TIME_INC_DAYS);    iday =    -1; isec = 0
            END SELECT
            CALL add_date (iday, isec, my_date)
            CALL TC_get               (my_date, yr, mo, dy, hr, mn, se)

          CASE(TIME_INC_MONTHS)           ! skip back one month
            IF (mo == 1) THEN
              mo = 12; yr = yr - 1
            ELSE
              mo = mo - 1
            END IF

          CASE(TIME_INC_YEARS)
            yr = yr - 1

          END SELECT

        ELSE
          CALL TC_get (my_date, yr, mo, dy, hr, mn, se)

        END IF

        ! adjust next trigger date to last environment in unit
        ! interval = <last second - delta, last second>
        SELECT CASE(unit)
        CASE(TIME_INC_MINUTES);                                    se = 59
        CASE(TIME_INC_HOURS);                             mn = 59; se = 59
        CASE(TIME_INC_DAYS);                     hr = 23; mn = 59; se = 59
        CASE(TIME_INC_MONTHS);                   hr = 23; mn = 59; se = 59
        SELECT CASE (get_calendar_type())
        CASE (JULIAN)
          dy = Get_JulianMonLen(yr,mo)  ! reset to the end of the month
        CASE (CYL360)
          dy = Get_Ly360MonLen()        ! reset to the end of the month
        END SELECT

        CASE(TIME_INC_YEARS);  mo = 12; dy = 31; hr = 23; mn = 59; se = 59
        END SELECT

        CALL TC_set            (yr, mo, dy, hr, mn, se, my_date)

      END SELECT

    END SUBROUTINE adjust_date

  END SUBROUTINE TE_get_next_trigger

  
  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE event_next_date (event, date) !************************************
    !-
    TYPE (time_event), INTENT (in)  :: event
    TYPE (time_days),  INTENT (out) :: date

    date = event %next_trigger

  END SUBROUTINE event_next_date

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE event_current_date (event, date) !*********************************
    !-
    TYPE (time_event), INTENT (in)  :: event
    TYPE (time_days),  INTENT (out) :: date

    date = event %current_trigger

  END SUBROUTINE event_current_date


  ! ----------------------------------------------------------------------------
  !+
  FUNCTION event_eval (event, delta) RESULT (slen) !****************************

    ! check the multiple of delta in event interval
    !-
    TYPE(time_event),   INTENT(in) :: event
    REAL(dp), OPTIONAL, INTENT(in) :: delta
    INTEGER                        :: slen

    INTEGER            :: fit

    slen = -1

    IF (.NOT. event%init) THEN
      CALL message('','event not initialised, no unit given.')

    ELSE
      SELECT CASE(event%unit)
      CASE(TIME_INC_SECONDS);   slen = event%count
      CASE(TIME_INC_MINUTES);   slen = event%count * 60
      CASE(TIME_INC_HOURS);     slen = event%count * 60 * 60
      CASE(TIME_INC_DAYS);      slen = event%count * 60 * 60 * 24

      CASE(TIME_INC_MONTHS, TIME_INC_YEARS)
        CALL message('','Exact value undefined.')
        CALL message('','Event trigger interval for months or years may varied.') 

      CASE default
        m_text = 'Counter unit unknown ::' // TRIM(event%unit)
        CALL finish('mo_time_event:event_eval',m_text)
      END SELECT

      IF (slen > 0 .AND. PRESENT(delta)) THEN
        fit = MOD(slen,INT(delta))
        IF (fit /= 0) THEN
          WRITE(m_text,*) 'Event <',TRIM(event%label),&
               '> interval not a multiple of ',INT(delta)
          CALL message('Warning',m_text)
          slen = -fit
        END IF
      END IF

    END IF

  END FUNCTION event_eval

END MODULE mo_time_event
!
! ------------------------------------------------------------------------------
