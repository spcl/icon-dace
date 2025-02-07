! mo_time_control.f90 - Time controller: Date and Time Control Interface Module
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!
! This file has been modified for the use in HD-model. 
!_________________________________________

MODULE mo_time_control

  ! ----------------------------------------------------------------------------
  !+
  !
  ! ECHAM Date and Time Control Interface Module
  ! --------------------------------------------
  !
  ! * interface between the time control modules and ECHAM
  !
  ! * definitions of global date/time constants
  !
  ! * the date/time functions are controlled by three structures
  !
  !  TIME_EVENT ..... time dependend handling of operations
  !
  !  TIME_MANAGER ... time axis informations of the model
  !                   like the start date, present position (= the time step)
  !                   time interval between two steps
  !
  !  TIME_DAYS ...... specific dates of the model history
  !                   like start date, stop date, ...
  !
  ! Authors:
  !
  ! I. Kirchner, MPI, April 2000
  ! I. Kirchner, MPI, October/December 2000
  ! I. Kirchner, MPI, March 2001, revision
  ! S. Legutke,  MPI, Aug   2001, separate events for model coupling read/write
  ! I. Kirchner, MPI, Sep/Oct 2001, revision
  ! I. Kirchner, MPI, Aug 2002, add change_present_date
  ! I. Kirchner, FUB, February 2003, revision/code review
  ! L. Kornblueh MPI, April 2003, additional features for debugging and testing
  ! U. Schulzweida, MPI, March 2007, added weights for daily interpolation
  ! S. Rast, MPI, April 2010, added weights for interpolation with respect to 
  !                           radiation time step
  ! D. Klocke, MPI, Nov 2010, changed time step lenght for first time step in
  !                           NWP mode and reset switch after initial time step
  ! V. Gayler, MPI, May 2011, added function get_interval_seconds_next
  ! S. Hagemann, Hereon, Aug. 2024, lrad set to .False. as no radiation in HD model
  !
  ! external modules
  !
  ! ECHAM specific modules
  !   mo_kind
  !   mo_control
  !   mo_machine
  !   mo_exception
  !   mo_mpi
  !   mo_constants
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.
  !
  !-
  USE mo_kind,            ONLY: dp
  USE mo_machine,         ONLY: prec
  USE mo_control,         ONLY: lcouple, lhd, nn, lnwp, nlev, ldiagamip
  USE mo_exception,       ONLY: finish, message, message_text
  USE mo_mpi,             ONLY: p_bcast, p_all_comm
  USE mo_constants,       ONLY: api
  !+
  !
  ! date/time handling modules
  !   mo_time_conversion
  !   mo_time_manager
  !   mo_time_event
  !   mo_time_base
  !
  !-
  USE mo_time_conversion, ONLY: time_days, time_native, time_intern, month_len,&
                                TC_get, TC_set, TC_convert,                    &
                                OPERATOR(==), OPERATOR(<), OPERATOR(>),        &
                                TC_PRN_NATIVE, print_date, add_date, year_len, &
                                day_len, IMerge_HMS2Sec

  USE mo_time_manager,    ONLY: time_manager, manager_state, manager_init,     &
                                manager_print, rewind_manager

  USE mo_time_event,      ONLY: time_event, event_state, io_time_event,        &
                                event_init, event_reinit, event_print,         &
                                event_next_date, event_current_date,           &
                                event_eval, TE_print_event_name,               &
                                TRIG_LAST, TRIG_FIRST, TRIG_EXACT, TRIG_NONE,  &
                                TIME_INC_SECONDS, TIME_INC_MINUTES,            &
                                TIME_INC_HOURS , TIME_INC_DAYS,                &
                                TIME_INC_MONTHS, TIME_INC_YEARS

  USE mo_time_base,       ONLY: IDAYLEN, get_calendar_type, JULIAN, CYL360,    &
                                julian_date, Set_JulianDay, Get_JulianYearDay, &
                                ly360_date,  Set_LY360Day,  Get_Ly360YearDay,  &
                                sec2frac


  IMPLICIT NONE

  PUBLIC

  ! ----------------------------------------------------------------------------
  ! internal used variables

  INTEGER, PARAMETER, PRIVATE :: &
       STR_LEN_A = 20      ,&! predefined string len A
       STR_LEN_B = 256       ! predefined string len B

  ! the adjustment of events for RERUN is dependent on the trigger step
  ! the trigger step can be the present date or the next date
  !
  CHARACTER(len=*), PARAMETER, PRIVATE :: &
       EV_TLEV_PRES    = 'present'        ,&! check event with present date
       EV_TLEV_NEXT    = 'next'           ,&! check event with next date
       TIME_INC_STEPS  = 'steps'          ,&! special event interval unit
       TIME_INC_ALWAYS = 'always'           ! special event always used  

  ! define here the default format for date print outs
  ! possible values are TC_PRN_DAYS, TC_PRN_INTERN, TC_PRN_NATIVE
  !
  CHARACTER(len=*), PARAMETER, PRIVATE :: PRN_DATE_FORMAT = TC_PRN_NATIVE 

  CHARACTER(len=STR_LEN_B), PRIVATE    :: m_text

  INTEGER, PARAMETER ::    & ! define label types formatting date/time in filenames
       FLT_YM         = 2, &
       FLT_YMD        = 3, &
       FLT_YMDH       = 4, &
       FLT_YMDHM      = 5, &
       FLT_YMDHMS     = 6, &
       FLT_ISO_YMD    = 7, &
       FLT_ISO_YMDHMS = 8, &
       FLT_NWP_YMD    = 9

  ! ----------------------------------------------------------------------------
  !+
  !
  ! ----------------------------------------------------------------------------
  ! ***** definition of DATE/TIME variables/constants

  INTEGER, PARAMETER :: NDAYLEN  = IDAYLEN     ! transfer the day length
  LOGICAL            :: ldebugev = .FALSE.     ! .T. print more infos

  REAL(dp)           :: delta_time    = 0.0_dp ! distance of adjacent times
  REAL(dp)           :: time_step_len = 0.0_dp ! forecast time step,
                                               ! at beginning equal delta_time
                                               ! else 2*delta_time

  ! dt_initial  corresponds to the initial condition date
  ! dt_start    defines the start of an experiment
  !
  TYPE(time_days),SAVE :: initial_date       ! should set at initial time from file
  !
  INTEGER, PARAMETER   :: INIT_STEP   = 0    ! initial time step
  INTEGER              :: dt_start(6) = 0    ! (runctl) start date of experiment
                                             ! meaning (yr, mo, dy, hr, mi, se)
  TYPE(time_days),SAVE :: start_date         ! transformed start date
  LOGICAL              :: lstart    = .TRUE. ! .TRUE. for the first time step

  LOGICAL            :: lfirst_day = .TRUE.  ! .TRUE. during the first day
  LOGICAL            :: l2nd_day   = .TRUE.  ! .TRUE. during the first+second day

  INTEGER              :: dt_resume(6) = 0   ! user defined restart time

  INTEGER              :: dt_stop(6) = 0     ! (runctl) stop experiment here
                                             ! meaning (yr, mo, dy, hr, mi, se)
  TYPE(time_days),SAVE :: stop_date          ! transformed stop date
  LOGICAL              :: lbreak   = .FALSE. ! .TRUE. at end of one time segment
  LOGICAL              :: lstop    = .FALSE. ! .TRUE. during the last time step
  LOGICAL              :: lresume = .FALSE.  ! .TRUE. during rerun step
  TYPE(time_days),SAVE ::  resume_date       ! transformed rerun date

  TYPE(time_days),SAVE ::  previous_date     ! date at (time - delta_time)
  TYPE(time_days),SAVE ::   current_date     ! date at (time)
  TYPE(time_days),SAVE ::      next_date     ! date at (time + delta_time)

  TYPE (time_days),SAVE :: radiation_date    ! date corresponding to rad_calc
  TYPE (time_days),SAVE :: prev_radiation_date ! date corresponding to previous
                                             ! radiation calc (even for first time step
  LOGICAL           :: l_orbvsop87 = .TRUE.  ! .TRUE. : orbit routine from vsop87
                                             ! .FALSE.: orbit routine pcmdi (AMIP)
                                             ! negative number force default value
  INTEGER          :: no_days      = -1      ! stop at (init/restart_time + no_days)
  INTEGER          :: no_steps     = -1      ! stop at (init_time + no_steps)

  INTEGER          :: no_cycles        =  1         ! maximal number of rerun intervals
  INTEGER, PRIVATE :: nmcount          =  1         ! rerun loop counter
  LOGICAL, PRIVATE :: l_need_trigfiles = .FALSE.    ! needed t get restart cycles corrrect
  LOGICAL          :: lfirst_cycle     = .TRUE.     ! .TRUE. during the first rerun cycle
  LOGICAL          :: lstop_rerun      = .FALSE.    ! couple lstop and l_putrerun


  ! ----------------------------------------------------------------------------
  ! ***** time axis and event definitions

  TYPE(time_manager),SAVE :: echam_time  !  (internal) echam time axis

  ! external modifcation of events is handled with io_time_events
  ! and transformed to the internal structure time_event
  ! at the beginning of the time loop in STEPON all
  ! events are checked and a corresponding logical is set

  ! ----------------------------------------------------------------------------
  ! *** RUNCTL events

  INTEGER, PARAMETER     :: NO_PUTDATA = 30  ! number of independent io events
  INTEGER                :: idx_putdata = 1  ! last allocated putdata event
  LOGICAL                ::  l_putdata(NO_PUTDATA) = .FALSE.
  TYPE(time_event),SAVE  :: ev_putdata(NO_PUTDATA)
  TYPE(io_time_event),SAVE ::    putdata    &! post-processing frequency
       = &
       io_time_event(12,TIME_INC_HOURS,TRIG_FIRST,0)

  INTEGER               :: timelabel_type = FLT_ISO_YMD
  LOGICAL               ::  l_trigfiles = .FALSE.
  TYPE(time_event),SAVE :: ev_trigfiles
  TYPE(io_time_event),SAVE ::    trigfiles  &! generate new file names
       = &
       io_time_event(1,TIME_INC_MONTHS,TRIG_FIRST,0)

  LOGICAL               ::  l_putrerun  = .FALSE.
  TYPE(time_event),SAVE :: ev_putrerun
  TYPE(io_time_event),SAVE ::    putrerun   &! rerun write-up interval
       = &
       io_time_event(1,TIME_INC_MONTHS,TRIG_LAST,0)

  ! ----------------------------------------------------------------------------
  ! atmosphere ocean coupling 

  LOGICAL               ::  l_getocean = .FALSE.
  TYPE(time_event),SAVE :: ev_getocean
  TYPE(io_time_event),SAVE ::    getocean   &! getting data from the ocean
       = &
       io_time_event (1,TIME_INC_DAYS,TRIG_NONE,0)

  LOGICAL               ::  l_putocean = .FALSE.
  TYPE(time_event),SAVE :: ev_putocean
  TYPE(io_time_event),SAVE ::    putocean   &! transfer data to the ocean
       = &
       io_time_event (1,TIME_INC_DAYS,TRIG_NONE,0)

  ! coupling with hydrological discharge model

  LOGICAL               ::  l_puthd = .FALSE. ! putting data to the ocean
  TYPE(time_event),SAVE    :: ev_puthd
  TYPE(io_time_event),SAVE ::    puthd       &! transfer data to HD model
       = &
       io_time_event(1,TIME_INC_DAYS,TRIG_NONE,0)

  LOGICAL               ::  l_gethd = .FALSE. ! getting data from the ocean
  TYPE(time_event),SAVE :: ev_gethd
  TYPE(io_time_event),SAVE ::    gethd       &! get data from HD model
       = &
       io_time_event(1,TIME_INC_DAYS,TRIG_NONE,0)

  ! ----------------------------------------------------------------------------
  ! subjob handling

  LOGICAL               :: lsub     = .FALSE. ! .TRUE. submit subjobs
  INTEGER               :: nsub     = 0       ! user defined number of subjobs
  INTEGER, PARAMETER    :: NSUB_MAX = 9       ! maximal number of subjobs
  LOGICAL               :: subflag(NSUB_MAX) = .FALSE. ! bind subjobs to output streams

  LOGICAL               ::  l_trigjob(NSUB_MAX) = .FALSE.
  TYPE(time_event),SAVE :: ev_trigjob(NSUB_MAX)
  TYPE(io_time_event),SAVE ::    trigjob(NSUB_MAX) &! subjob submit-interval
       = &
       io_time_event (1,TIME_INC_MONTHS,TRIG_NONE,0)

  ! ----------------------------------------------------------------------------
  ! *** AMIP2 diagnostic events

  LOGICAL               ::  l_diagamip = .FALSE.
  TYPE(time_event),SAVE :: ev_diagamip
  TYPE(io_time_event),SAVE ::    diagamip         &! event at first step of a day
       = &
       io_time_event(1,TIME_INC_DAYS,TRIG_NONE,0)

  ! ----------------------------------------------------------------------------
  ! *** DYNCTL events

  LOGICAL               ::  l_diagdyn = .FALSE.
  TYPE(time_event),SAVE :: ev_diagdyn
  TYPE(io_time_event),SAVE ::    diagdyn           &! freq. of dynamical diag.
       = &
       io_time_event(5,TIME_INC_DAYS,TRIG_NONE,0)

  LOGICAL               ::  l_diagvert = .FALSE.
  TYPE(time_event),SAVE :: ev_diagvert
  TYPE(io_time_event),SAVE ::    diagvert          &! freq. of vertical dyn. diag.
       = &
       io_time_event(5,TIME_INC_DAYS,TRIG_NONE,0)

  ! ----------------------------------------------------------------------------
  ! *** RADCTL events

  LOGICAL               ::  l_trigrad = .TRUE.
  TYPE(time_event),SAVE :: ev_trigrad
  TYPE(io_time_event),SAVE ::    trigrad           &! freq. of full radiation
       = &
       io_time_event (2,TIME_INC_HOURS,TRIG_FIRST,0)

  LOGICAL               ::  l_trigradm1= .TRUE.     ! time step before full radiation
  TYPE(time_event),SAVE :: ev_trigradm1
  TYPE(io_time_event),SAVE ::    trigradm1         &! freq. of full radiation
       = &
       io_time_event (2,TIME_INC_HOURS,TRIG_FIRST,0)

  ! ----------------------------------------------------------------------------
  ! ***** definitions for nudging

  INTEGER         :: dt_nudg_start(6) = 0  ! first nudging time
  TYPE(time_days),SAVE ::    nudg_start         ! transformed date
  INTEGER         :: dt_nudg_stop(6)  = 0  ! last nudging time
  TYPE(time_days),SAVE ::    nudg_stop          ! transformed date

  TYPE(time_days), SAVE :: &
       ndg_date0, ndg_date1, ndg_date2, ndg_date3   ! nudging time window
  TYPE(time_days), SAVE :: ndg_inp_date                   ! nudging input date
  TYPE(time_days), SAVE :: nudg_sst_time                  ! sst record date


  ! ----------------------------------------------------------------------------
  ! ***** Functions/Subroutines
  !-

  PUBLIC  :: p_bcast_event    ! (IO:io_time_event,I:isrc[,I:icom]) broadcast event
  PUBLIC  :: print_events     ! prints namelist changeable events with ldebugev enabled
  PUBLIC  :: ec_manager_init  ! (I:deltatime,I:timestep) initialize time manager

  PUBLIC :: init_manager       ! () reinitialize time manager
  PUBLIC :: change_current_date! (I:date) change the manager present date 
  PUBLIC :: init_times         ! () initialize all dates
  PUBLIC :: init_events        ! () initialize all events
  PUBLIC :: get_new_ev_putdata ! (I:io_event) allocate new putdata event

  PUBLIC :: init_nudgingtime_a ! () set nudging period
  PUBLIC :: init_nudgingtime_b ! () adjust time manager for nudging
  PUBLIC :: init_nudgingtime_c ! (I:linear) calculate nudging restart date
  PUBLIC :: init_nudgingtime_d ! () foreward calculation for next open

  PUBLIC :: get_step_from_header
  PUBLIC :: nudging_date_fit
  PUBLIC :: skip_nudging_read

  PUBLIC :: time_set     ! evaluation of events at beginning of the time loop
  PUBLIC :: time_reset   ! evaluation of events at the end of the time loop

  PUBLIC :: echam_ev_print        ! (I:event[,I:format]) print out event
  PUBLIC :: get_interval_steps    ! (I:event) returns steps
  PUBLIC :: get_interval_seconds  ! (I:event) returns seconds
  PUBLIC :: get_interval_seconds_next ! (I:event) returns seconds till next event

  INTERFACE write_date            ! print out date formatted
    MODULE PROCEDURE write_date_days    ! (I:time_days[,I:text])
    MODULE PROCEDURE write_date_intern  ! (I:time_intern[,I:text])
    MODULE PROCEDURE write_date_native  ! (I:time_native[,I:text])
  END INTERFACE

  PUBLIC :: set_delta_time   ! (I:trunc) preset delta time

  INTERFACE set_start_date   ! set start date from echam time manager
    MODULE PROCEDURE set_start_date_days    ! (I:time_days)
    MODULE PROCEDURE set_start_date_intern  ! (I:time_intern)
    MODULE PROCEDURE set_start_date_native  ! (I:time_native)
  END INTERFACE

  INTERFACE get_start_date   ! get start date from echam time manager
    MODULE PROCEDURE get_start_date_days   ! (O:time_days)
    MODULE PROCEDURE get_start_date_intern ! (O:time_intern)
    MODULE PROCEDURE get_start_date_native ! (O:time_native)
  END INTERFACE

  PUBLIC  :: get_date_components ! (I:time_days,O:yr,mo,dy,hr,mn,se)
  PUBLIC  :: get_forecast_hours  ! used in forecast mode
  PRIVATE :: get_delta_time               ! get delta time from echam time manager
  PUBLIC  :: get_time_step       ! provide present time step

  PUBLIC :: get_month_len   ! (I:year,I:month) provide length of month
  PUBLIC :: get_year_len    ! ([I:year]) provide length of the year
  PUBLIC :: get_year_day    ! (I:time_days) calculates the day in the year

  INTERFACE str_date  ! convert time information into a string
    MODULE PROCEDURE str_date_days   ! (I:format,I:time_days)
    MODULE PROCEDURE str_date_intern ! (I:format,I:time_intern)
    MODULE PROCEDURE str_date_native ! (I:format,I:time_native)
  END INTERFACE

  PUBLIC :: inp_convert_date  ! (I:yymmdd,I:hhmmss,O:time_days)
  PUBLIC :: out_convert_date  ! (I:time_days,O:yymmdd,I:hhmmss)

  PUBLIC :: get_clock ! provide the daytime fraction for orbit
  PUBLIC :: get_cycle ! provide the current cycle   

  PUBLIC :: day_difference  ! Calculate the difference in days for two ISO dates

CONTAINS
  !+
  ! ------------------------------------------------------------------------------
  !
  ! low level initialization, predefinition of structures
  !

  ! ----------------------------------------------------------------------------
  !
  SUBROUTINE p_bcast_event (event, p_source, comm)  !***************************

    ! distribute events to all nodes
    !-
    TYPE(io_time_event), INTENT(inout) :: event
    INTEGER,             INTENT(in)    :: p_source
    INTEGER, OPTIONAL,   INTENT(in)    :: comm

    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL p_bcast(event%counter,    p_source, p_comm)
    CALL p_bcast(event%unit,       p_source, p_comm)
    CALL p_bcast(event%adjustment, p_source, p_comm)
    CALL p_bcast(event%offset,     p_source, p_comm)

  END SUBROUTINE p_bcast_event

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE print_events               !***** called in CONTROL ***********

    ! print all ECHAM event structures
    !-
    INTEGER :: i

    ! **** external controled events -------------------------------------------

    IF (ldebugev) THEN

      CALL message('TIMECONTROL',' ********* default settings ***********')

      CALL message('',''); CALL message('','&RUNCTL')
      WRITE(m_text,*) 'LRESUME=    ',lresume,   ','; CALL message('',m_text)
      WRITE(m_text,*) 'LDEBUGEV=   ',ldebugev,  ','; CALL message('',m_text)
      WRITE(m_text,*) 'L_ORBVSOP87= ',l_orbvsop87,','; CALL message('',m_text)
      WRITE(m_text,*) 'DELTA_TIME= ',delta_time,','; CALL message('',m_text)
      WRITE(m_text,*) 'DT_START=   ',dt_start,  ','; CALL message('',m_text)
      WRITE(m_text,*) 'DT_RESUME=  ',dt_resume, ','; CALL message('',m_text)
      WRITE(m_text,*) 'DT_STOP=    ',dt_stop,   ','; CALL message('',m_text)
      WRITE(m_text,*) 'TRIGFILES=  ',trigfiles, ','; CALL message('',m_text)
      WRITE(m_text,*) 'PUTDATA=    ',putdata,   ','; CALL message('',m_text)
      WRITE(m_text,*) 'PUTRERUN=   ',putrerun,  ','; CALL message('',m_text)
      WRITE(m_text,*) 'PUTOCEAN=   ',putocean,  ','; CALL message('',m_text)
      WRITE(m_text,*) 'GETOCEAN=   ',getocean,  ','; CALL message('',m_text)
      WRITE(m_text,*) 'PUTHD=      ',puthd,     ','; CALL message('',m_text)
      WRITE(m_text,*) 'GETHD=      ',gethd,     ','; CALL message('',m_text)
      WRITE(m_text,*) 'NO_CYCLES=  ',no_cycles, ','; CALL message('',m_text)
      WRITE(m_text,*) 'NO_DAYS=    ',no_days,   ','; CALL message('',m_text)
      WRITE(m_text,*) 'NO_STEPS=   ',no_steps,  ','; CALL message('',m_text)
      WRITE(m_text,*) 'LSUB=       ',lsub,      ','; CALL message('',m_text)
      WRITE(m_text,*) 'NSUB=       ',nsub,      ','; CALL message('',m_text)
      WRITE(m_text,*) 'SUBFLAG=    ',subflag,   ','; CALL message('',m_text)
      WRITE(m_text,*) 'LDIAGAMIP=  ',ldiagamip;      CALL message('',m_text)

      DO i=1,NSUB_MAX
        WRITE(m_text,*)&
             'TRIGJOB(',i,')= ',TRIGJOB(i),','; CALL message('',m_text)
      END DO

      CALL message('',''); CALL message('','&DYNCTL')
      WRITE(m_text,*) 'DIAGDYN=    ',diagdyn,   ','; CALL message('',m_text)
      WRITE(m_text,*) 'DIAGVERT=   ',diagvert,  ','; CALL message('',m_text)

      CALL message('',''); CALL message('','&RADCTL')
      WRITE(m_text,*) 'TRIGRAD=    ',trigrad,   ','; CALL message('',m_text)

      CALL message('','')
      CALL message('TIMECONTROL',' **************************************')

    END IF

  END SUBROUTINE print_events


  !+
  ! ------------------------------------------------------------------------------
  !
  ! initialization during start/restart of the model
  !

  ! ----------------------------------------------------------------------------
  !
  SUBROUTINE ec_manager_init ( timestep, step )      !***** called in IO_init **

    ! echam time manager initialization, performed in IO_init
    !-
    INTEGER :: timestep   ! delta time in seconds
    INTEGER :: step       ! model time step

    CALL manager_init(echam_time,'echam time manager',start_date,timestep)

    IF (lresume) THEN
      CALL write_date(resume_date,'Experiment resumed at: ')

      IF (lfirst_cycle) THEN
        CALL manager_init(echam_time, step+1)
        CALL message('','Set manager to resumed time step + 1.')
      END IF

    ELSE
      CALL manager_state(echam_time,resume_date,step)

    END IF

  END SUBROUTINE ec_manager_init


  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE init_manager                !****** called in INITIALIZE **********

    ! time manager is initialized in IO_init, now only corrections
    !-

    ! init_manager is only called during the first rerun cycle

    REAL(dp)          :: zdtold
    INTEGER           :: istep
    TYPE(time_native) :: date_nat
    TYPE(time_days)   :: date_day
!!    INTEGER           :: days, seconds

    ! *** check and initialize delta time --------------------------------------

    delta_time = MAX(delta_time,0.0_dp)    ! time step must be positive
    IF (delta_time < prec) delta_time = set_delta_time(nn,nlev)

    CALL manager_state(echam_time,zdtold)

    IF (ABS(zdtold-delta_time) > 0.0_dp) THEN

      IF (lresume) THEN
        CALL message('','Attention: length of timestep has been changed!')
        WRITE (m_text,*) ' Old timestep was ', zdtold, ' seconds'
        CALL message('',m_text)
        WRITE (m_text,*) ' New timestep is  ', delta_time, ' seconds'
        CALL message('',m_text)
      END IF
      CALL manager_init(echam_time,delta_time)

    END IF

    IF (lresume) THEN

      lstart = .FALSE.

      ! *** check if resume date is changed ------------------------------------

      IF (SUM(dt_resume(:)) /= 0) THEN

        CALL TC_set(&
             dt_resume(1), dt_resume(2), dt_resume(3), &
             dt_resume(4), dt_resume(5), dt_resume(6), date_nat)
        CALL TC_convert(date_nat, date_day)
        CALL write_date(date_day,'Reset resume date from namelist: ')

        ! reset the present date with changing the start date
        ! the time step is conserved
        CALL rewind_manager(echam_time, date_day, .FALSE.)

      END IF

    END IF


    ! *** check if start date is changed ---------------------------------------

    IF (lnwp) THEN
      IF (.NOT. lresume) THEN
        CALL finish('mo_time_control:init_manager','NWP needs to start from restart file.')
      END IF
      IF (SUM(dt_start(:)) /= 0) THEN      
        CALL message('mo_time_control:init_manager','NWP ignores setting of DT_START in runctl.')
      END IF
      start_date = resume_date
      CALL write_date(start_date,'Forecast starts at ')
      CALL set_start_date(start_date)        ! reset the time manager start date
    ELSE

      IF (SUM(dt_start(:)) /= 0) THEN
        CALL TC_set(&
             dt_start(1), dt_start(2), dt_start(3), &
             dt_start(4), dt_start(5), dt_start(6), date_nat)
        CALL TC_convert(date_nat, start_date)

        CALL write_date(start_date,'Start date replaced by namelist start date: ')

        CALL set_start_date(start_date)        ! reset the time manager start date

        IF (.NOT. lresume) THEN

          !       set the time manager to the init_step now
          !       reset time step at beginning of experiment
          !       count backward

          istep = INIT_STEP - get_time_step()
          CALL manager_init(echam_time,istep)
        ELSE

          CALL manager_state(echam_time, date_day)
          IF (date_day < start_date) THEN
            CALL write_date(date_day,'Current date ...')
            CALL write_date(start_date,'Start date ...')
            CALL finish('mo_time_control:init_manager','Start date in future')

          ELSE IF (date_day == start_date) THEN
            lresume = .FALSE.
            CALL message('init_manager','Set start date to resume date, force initial run')
          END IF

        END IF

      END IF

    END IF

    ! ***  preliminary initializion of  previous date and next date ------------
    !      may be changed during nudging

    istep = get_time_step()
    CALL manager_state(echam_time,previous_date,istep-1)
    CALL manager_state(echam_time,next_date,    istep+1)

    CALL message('','Time step and start date evaluation done.') 

  END SUBROUTINE init_manager

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE change_current_date(new_date)
    ! rewind the time manager to a specific date without changing the start time

    TYPE(time_days), INTENT(in) :: new_date

    LOGICAL :: old_lresume
    INTEGER :: istep
 
    CALL manager_init(echam_time, .FALSE.)   !  unlock time manager
    CALL rewind_manager(echam_time, new_date, .TRUE.) 

    CALL write_date(new_date,'Reset resume date internal: ')
 
    CALL manager_init(echam_time, .FALSE.)  !  lock time manager
    CALL manager_print(echam_time)

    ! correct dates
    istep = get_time_step()
    CALL manager_state(echam_time, current_date)
    CALL manager_state(echam_time,    next_date,istep+1)

    ! reset events
    old_lresume = lresume
    lresume = .TRUE.
    CALL init_events
    lresume = old_lresume

    lstart = (istep == INIT_STEP) 
    IF (lstart) THEN
      time_step_len = delta_time
      lresume = .FALSE.
    ELSE
      time_step_len = 2._dp*delta_time
    END IF
 
  END SUBROUTINE change_current_date

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE init_times   !**** called in INITIALIZE ***************************

    ! initialize dates
    !-
    INTEGER           :: istep, incr
    TYPE(time_native) :: date_nat

    ! *** no changes of time manager possible from here ------------------------

    IF (lfirst_cycle) CALL manager_init(echam_time,.TRUE.)
    CALL manager_print(echam_time)

    ! *** define time stepping -------------------------------------------------

    istep  = get_time_step()
    
!    write(0,*) '+++++HD: mo_time_control: sub.init_times: INIT_STEP=',INIT_STEP
    
    lstart = (istep == INIT_STEP) 
    IF (lstart) THEN
      time_step_len = delta_time
    ELSE
      time_step_len = 2._dp*delta_time
    END IF

    ! *** start date of manager can be changed ---------------------------------
    !     final definition of start date here

    CALL manager_state(echam_time,current_date)
    CALL get_start_date(start_date)

    ! *** stop of experiment evaluated only during first rerun cycle -----------
    !     evaluation in the following order (highest priority left)
    !     NO_STEPS or NO_DAYS or DT_STOP or default

    IF (lfirst_cycle) THEN

      IF (SUM(dt_stop(:)) /= 0) THEN                    ! evaluate DT_STOP -----
        CALL TC_set(&
             dt_stop(1), dt_stop(2), dt_stop(3), &
             dt_stop(4), dt_stop(5), dt_stop(6), date_nat)
        CALL TC_convert(date_nat, stop_date)
        CALL message('','Using DT_STOP for model stop.')

      ELSE IF (no_steps < 0 .AND. no_days < 0) THEN     ! set default ----------
        no_steps = 10
        WRITE(m_text,*) & 
             'no predefined NO_DAYS/NO_STEPS: set NO_STEPS = ',no_steps
        CALL message('',m_text)

      ELSE IF (no_steps == 0 .OR. no_days ==0) THEN     ! not allowed ----------
        CALL finish('mo_time_control:init_times','NO_DAYS or NO_STEPS equal zero')

      END IF

      IF (lnwp) no_days = MIN(no_days,10)  ! stop after maximal 10 days

      IF (no_days > 0) THEN                             ! evaluate NO_DAYS -----
        stop_date = current_date
        CALL add_date(no_days, 0 , stop_date)
        CALL message('','Using NO_DAYS for model stop.')
      END IF

      IF (no_steps > 0) THEN                            ! evaluate NO_STEPS ----
        stop_date = start_date
        incr = no_steps * INT(delta_time)
        CALL add_date(0, incr, stop_date)
        CALL message('','Using NO_STEPS for model stop.')
      END IF

    ELSE

      CALL message('','No evaluation of model stop date.') 

    END IF

    ! *** check date order -----------------------------------------------------

    IF (.NOT. (start_date < stop_date)) THEN
      CALL write_date(start_date,' Start date: ')
      CALL write_date(stop_date, ' Stop date : ')
      CALL finish('mo_time_control:init_times','Start date larger/equal than stop date ....')
 
    ELSE IF (stop_date < current_date) THEN
      CALL write_date(current_date,' Current date: ')
      CALL write_date(stop_date,   ' Stop date   : ')
      CALL finish('mo_time_control:init_times','Current date larger than stop date ....')
 
    END IF
 
    CALL write_date(stop_date,'Stop experiment at: ')
 
    ! *** initialize previous date and next date -------------------------------
 
    CALL manager_state(echam_time,previous_date,istep-1)
    CALL write_date              (previous_date,'Previous date: ')

    CALL manager_state(echam_time,next_date,    istep+1)
    CALL write_date              (next_date,    'Next date    : ')

    CALL message ('','END SUBROUTINE init_times')
 
  END SUBROUTINE init_times


  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE init_events     !********* called in INITIALIZE *******************

    ! initialize all events
    !-
    CHARACTER(len=40) :: sub_txt
    INTEGER           :: i

    ! *** check weather prediction mode ----------------------------------------

    IF (lnwp) THEN
      IF (putdata%unit /= TIME_INC_HOURS) THEN
        CALL message('','Only HOURS as interval in putdata in LNWP mode are valid.')
      END IF
      trigfiles = putdata
    END IF

    ! *** define the time label for the output file name generation ------------

    SELECT CASE(trigfiles%unit)
    CASE(TIME_INC_MONTHS);  timelabel_type = FLT_YM
    CASE(TIME_INC_DAYS);    timelabel_type = FLT_YMD
    CASE(TIME_INC_HOURS);   timelabel_type = FLT_YMDH
    CASE(TIME_INC_MINUTES); timelabel_type = FLT_YMDHM
    CASE(TIME_INC_SECONDS); timelabel_type = FLT_YMDHMS
    CASE default
      timelabel_type = FLT_ISO_YMD
      CALL message('','Use default time label in file names (YYYYMM_DD).')
    END SELECT

    ! *** History and postprocessing events ------------------------------------

    CALL echam_ev_init(ev_putrerun, putrerun, 'rerun interval',    EV_TLEV_NEXT)
    CALL echam_ev_init(ev_trigfiles,trigfiles,'rebuild file names',EV_TLEV_NEXT)
    CALL echam_ev_init(ev_putdata(1),  putdata,  'output interval',EV_TLEV_NEXT)

    ! *** Coupling events ------------------------------------------------------

    IF ( lcouple .AND. (&
         (TRIM(putocean%adjustment) == TRIM(TRIG_NONE)) .OR.  &
         (TRIM(getocean%adjustment) == TRIM(TRIG_NONE)) &
         ))     &
         CALL finish('mo_time_control:init_events','definition of OCEAN coupling inconsistent')
 
    CALL echam_ev_init(ev_getocean,getocean,'couple get-from-ocean',EV_TLEV_PRES)
    CALL echam_ev_init(ev_putocean,putocean,'couple put-to-ocean',  EV_TLEV_NEXT)
 
    IF (lcouple .AND. lhd) THEN
      puthd = putocean
      gethd = getocean
      CALL message('',&
           'Synchronization of ocean and hydrological discharge model done.')
    END IF

    IF (lhd .AND. ( &
         (TRIM(gethd%adjustment) == TRIM(TRIG_NONE)) .OR.  &
         (TRIM(puthd%adjustment) == TRIM(TRIG_NONE)) &
         ) ) &
         CALL finish('mo_time_control:init_events','definition of HD events inconsistent')
    CALL echam_ev_init(ev_gethd,gethd,'couple get-from-hd',EV_TLEV_PRES)
    CALL echam_ev_init(ev_puthd,puthd,'couple put-to-hd',  EV_TLEV_NEXT)

    ! *** subjob control -------------------------------------------------------

    nsub = MIN(MAX(nsub,0),NSUB_MAX)
    IF (nsub == 0) lsub = .FALSE.
    DO i=1,nsub
      WRITE(sub_txt,*) 'subjob no ',i,' interval'
      CALL echam_ev_init(ev_trigjob(i),trigjob(i),sub_txt,EV_TLEV_NEXT)
    END DO
 
    ! *** radiation ------------------------------------------------------------
    CALL echam_ev_init(ev_trigrad,trigrad,'radiation computation', EV_TLEV_PRES)
    trigradm1 = trigrad
    CALL echam_ev_init(ev_trigradm1,trigradm1,'albedo computation', EV_TLEV_NEXT)

    ! *** diagnostics ----------------------------------------------------------
    IF (ldiagamip) diagamip%adjustment = TRIG_FIRST
    CALL echam_ev_init(ev_diagamip,diagamip, 'AMIP2 table 4 diagnostics',EV_TLEV_PRES)
    CALL echam_ev_init(ev_diagdyn,diagdyn, 'dynamical diagnostics',EV_TLEV_PRES)
    CALL echam_ev_init(ev_diagvert,diagvert,'vertical dyn-diag',   EV_TLEV_PRES)

    CALL message ('',' ')

  END SUBROUTINE init_events

  SUBROUTINE echam_ev_init (event, io_ev, ev_name, eval_date)  !**************

    ! the initial state of an event is calculated
    ! it is dependent on the initial date of a run
    ! at rerun the event triggers will be recalculated from beginning

      TYPE(time_event),    INTENT(inout) :: event     ! will be initialized
      TYPE(io_time_event), INTENT(inout) :: io_ev     ! pass external settings
      CHARACTER(len=*)                   :: ev_name   ! submit a name
      CHARACTER(len=*), OPTIONAL         :: eval_date ! the evaluation date
                                                      ! (next or present)
 
      TYPE(time_days)          :: my_date
      REAL(dp)                 :: zdtime, zztime
      CHARACTER(len=STR_LEN_A) :: newunit
      INTEGER                  :: newcount, idtime

      zztime = get_delta_time()
      zdtime = zztime - 1.0_dp
      idtime = NINT(zztime)

      ! *** convert from steps into other units --------------------------------

      IF (io_ev%unit == TIME_INC_STEPS .OR. io_ev%unit == TIME_INC_ALWAYS) THEN

        IF (io_ev%unit == TIME_INC_ALWAYS) THEN 
          ! for debugging and testing only
          IF (no_cycles == 0) THEN
            io_ev%counter = 1
          ELSE
            IF (ev_name(1:14) == 'rerun interval') THEN
              io_ev%counter = 0
              io_ev%adjustment = TRIG_NONE
              lstop_rerun = .TRUE.
            ELSE
              io_ev%counter = 1
            END IF
          END IF
        END IF

        CALL convert_steps2inter(io_ev%counter,zztime,newunit,newcount)
        io_ev%unit    = TRIM(newunit)
        io_ev%counter = newcount

        SELECT CASE(io_ev%unit)
        CASE(TIME_INC_SECONDS)
          CALL message('','Convert event interval from STEPS into SECONDS')

        CASE(TIME_INC_MINUTES)
          CALL message('','Convert event interval from STEPS into MINUTES')

        END SELECT

      END IF

      ! *** initialize event with external settings ----------------------------

      CALL event_init(event, ev_name, &
           io_ev%counter, io_ev%unit, io_ev%adjustment, zdtime, io_ev%offset)

      ! *** check time stepping against event intervals ------------------------

      SELECT CASE(io_ev%adjustment)
      CASE(TRIG_FIRST,TRIG_LAST)

        SELECT CASE(io_ev%unit)
        CASE(TIME_INC_SECONDS, TIME_INC_MINUTES, TIME_INC_HOURS, TIME_INC_DAYS)
          IF (event_eval(event,zztime) < 0) THEN
            CALL echam_ev_print(event)
            CALL finish('mo_time_control:echam_ev_init',&
                 'Event counter mismatch with time stepping.')
          END IF

        CASE(TIME_INC_MONTHS, TIME_INC_YEARS)
          CALL message('',&
               'No time stepping mismatch check defined for MONTHS and YEARS.')

        CASE default
          CALL finish('mo_time_control:echam_ev_init','Event unit not defined.')

        END SELECT
      END SELECT

      ! *** define initial date and trigger dates ------------------------------

      CALL get_start_date(my_date)
      ! for debugging and testing only
      IF (no_cycles == 0) my_date = resume_date
      IF (eval_date == EV_TLEV_NEXT) THEN
        ! next_date event never triggered at the second time step

        CALL add_date(0,idtime,my_date)
        CALL event_init(event,my_date)
        CALL get_start_date(my_date)
        CALL event_reinit(event,my_date)

      ELSE         
        CALL event_init(event,my_date)

      END IF

      ! *** find next trigger date starting at initial date --------------------

      IF (lresume) THEN

        !================= preliminar
        ! the finding of the next possible trigger can take a lot of
        ! time if the present date is very fare from the start date
        !
        ! revision needed with a new I/O concept:
        !     event elements must be available in a rerun file
        
        DO
          CALL event_next_date(event,my_date)

          IF (PRESENT(eval_date)) THEN
            SELECT CASE(eval_date)
            CASE(EV_TLEV_PRES)          ! check with current date
              IF (.NOT.(current_date > my_date)) EXIT
              
            CASE(EV_TLEV_NEXT)          ! check with next date
              IF (.NOT.   (next_date > my_date)) EXIT
              
            END SELECT

          ELSE
            IF (resume_date < my_date) EXIT

          END IF

          ! next date smaller as current date, rotate all dates
          CALL event_init(event,my_date,.TRUE.)

        END DO

      END IF

      ! *** print out event settings -------------------------------------------
      
      CALL echam_ev_print(event,ldebugev)

    END SUBROUTINE echam_ev_init

    SUBROUTINE convert_steps2inter (steps, dt, unit, count)  !****************
        
      ! convert steps into normal event units

      INTEGER,          INTENT(in)  :: steps  ! interval in steps
      REAL(dp),         INTENT(in)  :: dt     ! length of on step in seconds
      CHARACTER(len=*), INTENT(out) :: unit   ! new unit
      INTEGER,          INTENT(out) :: count  ! new counter

      INTEGER :: isdt, imdt, ihdt, iddt, isteps
      REAL(dp):: rdt

      INTEGER, PARAMETER :: imax = 2000000000

      isteps = steps
      isdt   = NINT(dt)
      rdt    = REAL(isdt,dp)

      IF ( (dt-rdt) > 0.0_dp) THEN ! delta time has fractional seconds
        WRITE(m_text,*) 'Delta time can not converted to integer ',rdt,dt
        CALL finish('mo_time_control:convert_steps2inter',m_text)
      END IF

      IF (isteps > imax/NINT(dt) ) THEN
        WRITE(m_text,*) 'step number too large: ',&
             ' maximum is ',imax/NINT(dt),' < ',isteps
        CALL finish('mo_time_control:convert_steps2inter',m_text)
      END IF

      isdt = NINT(dt*isteps)
      imdt = isdt /        60
      ihdt = isdt /    (60*60)
      iddt = isdt / (24*60*60)
      
      unit  = ''
      count = 0

      IF (iddt * 24*3600 == isdt) THEN      ! multiple of seconds per day
        unit  = TIME_INC_DAYS;    count = iddt

      ELSE IF (ihdt * 3600 == isdt) THEN    ! multiple of seconds per hour
        unit  = TIME_INC_HOURS;   count = ihdt

      ELSE IF (imdt * 60 == isdt) THEN      ! multiple of seconds per minute
        unit  = TIME_INC_MINUTES; count = imdt

      END IF

    END SUBROUTINE convert_steps2inter

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION get_new_ev_putdata(io_ev) RESULT(index)

    ! allocate new putdata events
    !+
    TYPE(io_time_event), INTENT(inout) :: io_ev
    INTEGER                            :: index
    CHARACTER(len=40)                  :: event_name = ' '

    index = 1    ! predefined

    IF (idx_putdata == NO_PUTDATA) &
         CALL finish('mo_time_control:get_new_ev_putdata','putdata events field full')

    idx_putdata = idx_putdata + 1
    WRITE(event_name,'(a,i3,a)') 'stream ',idx_putdata,' event'

    CALL echam_ev_init(ev_putdata(idx_putdata), io_ev, TRIM(event_name), EV_TLEV_NEXT)

    index = idx_putdata

  END FUNCTION get_new_ev_putdata


  !+
  ! ------------------------------------------------------------------------------
  !
  ! nudging interface procedures
  !

  ! ----------------------------------------------------------------------------
  !
  FUNCTION init_nudgingtime_a() RESULT (lnudgstop) !*** called in NUDGING_INIT *

    ! evaluate nudging period
    !-

    LOGICAL           :: lnudgstop
    TYPE(time_native) :: date_nat

    ! *** define nudging start date --------------------------------------------

    IF (SUM(dt_nudg_start(:)) /= 0) THEN
      CALL TC_set(&
           dt_nudg_start(1),dt_nudg_start(2),&
           dt_nudg_start(3),dt_nudg_start(4), &
           dt_nudg_start(5),dt_nudg_start(6), date_nat)
      CALL TC_convert(date_nat, nudg_start)
      CALL write_date(nudg_start,' Reset nudging start date from namelist: ')

    ELSE
      nudg_start = start_date
      CALL write_date(nudg_start,' Start with nudging at: ')

    END IF

    ! *** find nudging stop date -----------------------------------------------

    IF (SUM(dt_nudg_stop(:)) /= 0) THEN
      CALL TC_set(&
           dt_nudg_stop(1),dt_nudg_stop(2),&
           dt_nudg_stop(3),dt_nudg_stop(4), &
           dt_nudg_stop(5),dt_nudg_stop(6), date_nat)
      CALL TC_convert(date_nat, nudg_stop)
      CALL write_date(nudg_stop,' Reset nudging stop date from namelist: ')
      lnudgstop = .TRUE.

    ELSE
      CALL message('',' No external nudging stop time defined.')
      lnudgstop = .FALSE.

    END IF

    IF (lnudgstop .AND. (nudg_stop < nudg_start)) &
         CALL finish('mo_time_control:init_nudgingtime_a',&
         'Nudging stop time before nudging start time.')

    ! *** initialize nudging window dates --------------------------------------

    IF (lfirst_cycle) THEN
      CALL TC_set(0,0,ndg_date0)
      CALL TC_set(0,0,ndg_date1)
      CALL TC_set(0,0,ndg_date2)
      CALL TC_set(0,0,ndg_date3)
    END IF

  END FUNCTION init_nudgingtime_a

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE init_nudgingtime_b  !*** called in NUDGING_INIT *******************

    ! adjust time manager for nudging
    !-

    INTEGER :: istep
    LOGICAL :: lfit

    CALL write_date(ndg_inp_date,' Start reading nudging data at : ')
    CALL write_date(nudg_start,  ' Start nudging at              : ')

    ! *** check first possible nudging date with model time --------------------

    IF (nudg_start < start_date) THEN
      CALL message('',&
           ' Initial date of model newer - no adjustement performed.')

    ELSE IF (nudg_start < next_date) THEN
      CALL message('',&
           ' Next progostic date can be nudged - no adjustment performed.')

    ELSE
      CALL message('',' Adjust nudging using the next prognostic date.')

      ! *** which time step corresponds to the reference date ? ----------------

      CALL manager_state(echam_time,ndg_inp_date,istep,lfit)

      IF (.NOT. lfit) THEN
        CALL message('',&
             'Nudging reference date mismatch with time stepping - correction done.')
        istep = istep + 1
      END IF

      ! *** adjustment of the time manager -------------------------------------

      IF (lstart) THEN                ! reset the initial date
        start_date = ndg_inp_date
        CALL set_start_date(start_date)         
        istep = INIT_STEP - get_time_step()

      ELSE
        istep = istep - get_time_step() + 1

      END IF

      CALL manager_init(echam_time,istep)
      CALL message('',' Reset echam time manager position ...')

    END IF

  END SUBROUTINE init_nudgingtime_b

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE init_nudgingtime_c (llin) !*** called in NUDGING_INIT *************

    ! define the nudging restart date due to backward calculation
    !-
    LOGICAL :: llin        ! .TRUE. for linear interpolation

    INTEGER :: iday, isec
    LOGICAL :: lback, lset

    lback = .FALSE.
    lset  = .FALSE.

    ! *** control time interpolation type --------------------------------------

    IF (llin) THEN        ! linear time interpolation

      IF (next_date < ndg_date0) THEN
        lback = .TRUE.
      ELSE
        ndg_inp_date = ndg_date0
        lset = .TRUE.
      END IF

    ELSE                  ! non-linear time interpolation

      IF (next_date < ndg_date1) THEN
        lback = .TRUE.
      ELSE
        ndg_inp_date = ndg_date0
        lset = .TRUE.
      END IF

    END IF

    ! *** perform backward calculation -----------------------------------------

    IF (lback) THEN

      CALL TC_get  (ndg_date0,iday, isec)
      CALL add_date         (-iday,-isec,ndg_date1)  ! ndg_date1 - ndg_date0

      CALL TC_get (ndg_date1,iday,isec)

!!! debugging
      WRITE(m_text,*) 'delta time for backward calculation: ',iday,isec
      CALL message('init_nudgingtime_c',m_text)
!!! debugging

      DO
        CALL add_date(-iday,-isec,ndg_date0)   ! ndg_date0 - delta
!!! debugging
        CALL message('init_nudgingtime_c','skip back one delta time')
!!! debugging

        IF (.NOT. (ndg_date0 > next_date)) THEN
          IF (.NOT. llin) CALL add_date(-iday,-isec,ndg_date0)
          ndg_inp_date = ndg_date0
          lset = .TRUE.
          EXIT
        END IF
      END DO

    END IF

    IF (.NOT. lset) &
         CALL finish('mo_time_control:init_nudgingtime_c','Evaluation of nudging date failed.')

  END SUBROUTINE init_nudgingtime_c

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE init_nudgingtime_d  !*** called in NUDGING_IO *********************

    ! evaluate next nudging time for open
    ! until the end of one nudging data block is reached
    !-
    INTEGER :: iday, isec

    ! *** estimate time step interval

    CALL TC_get  (ndg_date1,iday, isec)
    ndg_inp_date = ndg_date2
    CALL add_date         (-iday,-isec,ndg_inp_date)   ! ndg_date2 - ndg_date1

    CALL TC_get (ndg_inp_date,iday,isec)
    ndg_inp_date = ndg_date2
    CALL add_date            (iday,isec,ndg_inp_date)  ! ndg_date2 + delta

    CALL write_date(ndg_inp_date,'Search next nudging data set at: ')

  END SUBROUTINE init_nudgingtime_d


  ! ----------------------------------------------------------------------------
  !+
  FUNCTION get_step_from_header (ymd, hms) RESULT (istep) !*** see NUDGING_IO **

    ! calculates the time step 
    ! corresponding to the nudging date in data record header
    !-
    INTEGER :: ymd, hms, istep
    LOGICAL :: lfit

    CALL inp_convert_date(ymd, hms, ndg_date3)
    CALL manager_state(echam_time, ndg_date3, istep, lfit)
    IF (.NOT.lfit) &
         CALL message('get_step_from_header',&
         'WARNING: Nudging date DO not fit to time stepping')

  END FUNCTION get_step_from_header

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION nudging_date_fit() RESULT (yes)   !*** called in NUDGING_IO *********

    ! check the nudging step after reading of new data
    !-
    LOGICAL :: yes

    yes = .FALSE.
    IF( next_date < ndg_date1) THEN      ! missing data at the beginning
      CALL write_date(ndg_date1,'Date in first DATA record: ')
      CALL finish('mo_time_control:nudging_date_fit','missing first nudging date')

    ELSE IF (                           &! window fit the calculation time step
         (ndg_date1 < next_date .OR. ndg_date1 == next_date) .AND. &
         (next_date < ndg_date2) ) THEN
      yes = .TRUE.

    END IF

  END FUNCTION nudging_date_fit

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION skip_nudging_read() RESULT (yes)   !*** called in NUDGING ***********

    ! check the time step at starting point of reading
    !-
    LOGICAL :: yes

    yes = &
         (ndg_date1 < next_date .OR. ndg_date1 == next_date) .AND. &
         (next_date < ndg_date2)

  END FUNCTION skip_nudging_read

  !+
  ! ------------------------------------------------------------------------------
  !
  ! event control inside the time loop (-> STEPON)
  !

  ! ----------------------------------------------------------------------------
  !
  SUBROUTINE time_set    !********** called in STEPON **************************

    ! evaluate events and date/time control elements at the
    ! beginning of the time loop in STEPON
    !-
    LOGICAL :: lrad = .FALSE.
    INTEGER         :: istep, isj, iput
    TYPE(time_days) :: my_day

    ! *** find current time step and current/next date -------------------------

    CALL manager_state(echam_time, istep)
    CALL manager_state(echam_time, current_date)
    CALL manager_state(echam_time,    next_date,istep+1)

    ! *** evaluate events ------------------------------------------------------
!ik_bugfix_R1.02c
    DO iput=1,idx_putdata
      l_putdata(iput)   = event_state(ev_putdata(iput),   next_date)
    END DO

!ik_bugfix_R1.02a
    l_trigfiles = event_state(ev_trigfiles, next_date) .OR. lfirst_cycle .OR. l_need_trigfiles

    IF(l_need_trigfiles) THEN
      l_need_trigfiles = .FALSE.
    ENDIF

    l_putrerun  = event_state(ev_putrerun,  next_date)

    l_getocean = lcouple .AND.(event_state(ev_getocean, current_date).OR.lstart)
    l_putocean = lcouple .AND. event_state(ev_putocean,    next_date)

    l_gethd = lhd .AND.(event_state(ev_gethd, current_date).OR.lstart)
    l_puthd = lhd .AND. event_state(ev_puthd,    next_date)

    DO isj = 1,nsub
      l_trigjob(isj) = event_state(ev_trigjob(isj), next_date)
    ENDDO

    l_diagamip  = event_state(ev_diagamip,  next_date)

    l_diagdyn  = event_state(ev_diagdyn,  current_date) .OR. lstart
    l_diagvert = event_state(ev_diagvert, current_date) .OR. lstart

    l_trigrad = event_state(ev_trigrad, current_date) .OR. lstart
    l_trigradm1 = event_state(ev_trigradm1, next_date) .OR. lstart

    IF(.NOT. lrad) THEN
      l_trigrad=.FALSE.
    ENDIF

    ! *** check end of first/second day ----------------------------------------

    IF (lfirst_day) THEN
      CALL get_start_date(my_day)
      CALL add_date(1,0,my_day)
!      lfirst_day = next_date < my_day     ! new version
      lfirst_day = previous_date < my_day  ! r0.10 solution
    END IF

    IF (l2nd_day) THEN
      CALL get_start_date(my_day)
      CALL add_date(2,0,my_day)
!      l2nd_day = next_date < my_day     ! new version
      l2nd_day = previous_date < my_day  ! r0.10 solution
    END IF

    ! *** evaluate the stop of the model ---------------------------------------
    ! HAG - original programming was that the first time step of the next day after the run period was run
    !       However, in HD the model stops in the end of run period
!!    lstop = (stop_date<next_date .OR. stop_date==next_date)
    lstop = (stop_date<current_date .OR. stop_date==current_date)

    IF (l_putrerun) THEN
      IF (no_cycles >= 1) THEN
        lbreak  = (nmcount >= no_cycles)  ! break rerun cycles
        nmcount = nmcount + 1             ! count rerun intervals
        l_need_trigfiles = .TRUE.
      ELSE
        lbreak  = .FALSE.                 ! break rerun cycles
      END IF
    ELSE
      lbreak = .FALSE.
      ! in 'model debug mode (no_cycles = 0)', put always rerun file before
      ! leaving integration loop in stepon
      IF (no_cycles == 0) THEN
        IF (lstop) &
           l_putrerun = .TRUE. 
      ELSE
        IF (lstop) THEN
          IF (lstop_rerun) THEN
            l_putrerun = .TRUE. 
          ELSE
            CALL message('Warning','stop model without rerun file generation!')
          END IF
        END IF
      END IF
    END IF

    ! *** print settings -------------------------------------------------------

    IF (ldebugev) THEN
      CALL message('time_set','----------------------------------------')
      CALL write_date(current_date,'Current date: ')
      CALL message('','Events triggered with current date ...')
      IF (l_diagamip)  CALL echam_ev_print(ev_diagamip)
      IF (l_diagdyn)   CALL echam_ev_print(ev_diagdyn)
      IF (l_diagvert)  CALL echam_ev_print(ev_diagvert)
      IF (l_trigrad)   CALL echam_ev_print(ev_trigrad)
      IF (l_getocean)  CALL echam_ev_print(ev_getocean)
      IF (l_gethd)     CALL echam_ev_print(ev_gethd)


      CALL write_date(next_date,   'Next date   : ')
      CALL message('','Events triggered with next date ...')
      DO iput=1,idx_putdata
        IF (l_putdata(iput))   CALL echam_ev_print(ev_putdata(iput))
      END DO
      IF (l_trigfiles) CALL echam_ev_print(ev_trigfiles)
      IF (l_putrerun)  CALL echam_ev_print(ev_putrerun)
      DO isj = 1,nsub
        IF (l_trigjob(isj)) CALL echam_ev_print(ev_trigjob(isj))
      END DO
      IF (l_trigradm1) CALL echam_ev_print(ev_trigradm1)
      IF (l_diagamip)  CALL echam_ev_print(ev_diagamip)
      IF (l_putocean)  CALL echam_ev_print(ev_putocean)
      IF (l_puthd)     CALL echam_ev_print(ev_puthd)

      IF(lfirst_day) CALL message('','time step during the first day')
      IF(l2nd_day)   CALL message('','time step during the first two days')

    END IF

    IF (lstop) THEN
      CALL write_date(next_date,'Stop model, last prognostic date/time is: ')
    ELSE IF (lbreak) THEN
      CALL write_date(next_date,'Interrupt model, last prognostic date/time is: ')
    END IF


    ! *** calculate integration interval ---------------------------------------

    delta_time = get_delta_time()
    IF (lstart) THEN
      time_step_len = delta_time
    ELSE IF (lnwp) THEN
      time_step_len = delta_time
    ELSE
      time_step_len = 2.0_dp*delta_time
    END IF

    ! *** prepare next radiation date ------------------------------------------

    IF (l_trigrad) THEN
      CALL radiation_time
      IF (lrad) CALL write_date(radiation_date,'Radiation calculated for : ')
!LK      IF (lrad) CALL write_date(prev_radiation_date,'Previous radiation date:')
    END IF

    ! *** set weighting factors for time interpolation of sst and sea ice ------
    CALL time_weights

    IF (ldebugev) CALL message('',' ')

  CONTAINS
    ! --------------------------------------------------------------------------
    !+
    SUBROUTINE time_weights  !**************************************************
      
      ! calculates weighting factores for clsst2 and ozone
      !-
      
      USE mo_interpo, ONLY : wgt1, wgt2, nmw1, nmw2, nmw1cl, nmw2cl, &
                             wgtd1, wgtd2, ndw1, ndw2
      USE mo_interpo, ONLY : wgt1_m, wgt2_m, nmw1_m, nmw2_m

      TYPE (time_native) :: date_monm1, date_mon, date_monp1
      INTEGER   :: yr, mo, dy, hr, mn, se
      INTEGER   :: days, seconds, isec
      INTEGER   :: imp1, imm1, imp1cl, imm1cl, imlenm1, imlen, imlenp1
      REAL (dp) :: zsec, zdayl
      REAL (dp) :: zmohlf, zmohlfp1, zmohlfm1
      REAL (dp) :: zdh, zdhp1, zdhm1
      
      ! ***  set calendar related parameters -----------------------------------
      
      CALL TC_convert(next_date,date_mon)
      CALL TC_get (date_mon, yr, mo, dy, hr, mn, se)
      
      ! month index for AMIP data  (0..13)
      imp1 = mo+1
      imm1 = mo-1
      
      ! month index for cyclic climatological data (1..12)
      imp1cl = mo+1
      imm1cl = mo-1
      IF (imp1cl > 12) imp1cl= 1
      IF (imm1cl <  1) imm1cl=12
      
      ! *** determine length of months and position within current month -------
      
      CALL TC_set(yr, imm1cl, 1, 0, 0, 0, date_monm1)
      CALL TC_set(yr, imp1cl, 1, 0, 0, 0, date_monp1)
      imlenm1 = month_len(date_monm1)
      imlen   = month_len(date_mon)
      imlenp1 = month_len(date_monp1)
      
      zdayl    = REAL(NDAYLEN,dp)
      zmohlfm1 = imlenm1*zdayl*0.5_dp
      zmohlf   = imlen  *zdayl*0.5_dp
      zmohlfp1 = imlenp1*zdayl*0.5_dp
      
      ! *** weighting factors for first/second half of month -------------------
      
      nmw1   = mo
      nmw1cl = mo
      
      ! seconds in the present month
      CALL TC_get (next_date, days, seconds)
      isec = (dy-1)*NDAYLEN + seconds
      zsec = REAL(isec,dp)
      
      IF(zsec <= zmohlf) THEN                     ! first part of month
        wgt1   = (zmohlfm1+zsec)/(zmohlfm1+zmohlf)
        wgt2   = 1._dp-wgt1
        nmw2   = imm1
        nmw2cl = imm1cl
      ELSE                                        ! second part of month
        wgt2   = (zsec-zmohlf)/(zmohlf+zmohlfp1)
        wgt1   = 1._dp-wgt2
        nmw2   = imp1
        nmw2cl = imp1cl
      ENDIF
      
      ! *** weighting factors for first/second half of day -------------------
      
      ndw1   = 2
      
      zsec = REAL(seconds,dp)
      zdh   = 12._dp*3600._dp
      zdhm1 = zdh
      zdhp1 = zdh
      IF( zsec <= zdh ) THEN                     ! first part of day
        wgtd1  = (zdhm1+zsec)/(zdhm1+zdh)
        wgtd2  = 1._dp-wgtd1
        ndw2   = 1
      ELSE                                       ! second part of day
        wgtd2  = (zsec-zdh)/(zdh+zdhp1)
        wgtd1  = 1._dp-wgtd2
        ndw2   = 3
      ENDIF
      
      !----------------------------------------------
      ! *** weighting factors for radiation time step
      !----------------------------------------------
      
      IF (l_trigrad) THEN
        
        CALL TC_convert(radiation_date,date_mon)
        CALL TC_get (date_mon, yr, mo, dy, hr, mn, se)
        
        ! month index for AMIP data  (0..13)
        imp1 = mo+1
        imm1 = mo-1
        
        ! month index for cyclic climatological data (1..12)
        imp1cl = mo+1
        imm1cl = mo-1
        IF (imp1cl > 12) imp1cl= 1
        IF (imm1cl <  1) imm1cl=12
        
        ! *** determine length of months and position within current month -------
        
        CALL TC_set(yr, imm1cl, 1, 0, 0, 0, date_monm1)
        CALL TC_set(yr, imp1cl, 1, 0, 0, 0, date_monp1)
        imlenm1 = month_len(date_monm1)
        imlen   = month_len(date_mon)
        imlenp1 = month_len(date_monp1)
        
        zdayl    = REAL(NDAYLEN,dp)
        zmohlfm1 = imlenm1*zdayl*0.5_dp
        zmohlf   = imlen  *zdayl*0.5_dp
        zmohlfp1 = imlenp1*zdayl*0.5_dp
        
        ! *** weighting factors for first/second half of month -------------------
        
        nmw1_m   = mo
        
        ! seconds in the present month
        CALL TC_get (radiation_date, days, seconds)
        isec = (dy-1)*NDAYLEN + seconds
        zsec = REAL(isec,dp)
      
        IF(zsec <= zmohlf) THEN                     ! first part of month
          wgt1_m   = (zmohlfm1+zsec)/(zmohlfm1+zmohlf)
          wgt2_m   = 1._dp-wgt1_m
          nmw2_m   = imm1
        ELSE                                        ! second part of month
          wgt2_m   = (zsec-zmohlf)/(zmohlf+zmohlfp1)
          wgt1_m   = 1._dp-wgt2_m
          nmw2_m   = imp1
        ENDIF
        
      END IF
      
    END SUBROUTINE time_weights

  END SUBROUTINE time_set

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE time_reset   !********* called in STEPON **************************

    ! evaluate actions at the end of the time step loop in STEPON
    !-
    INTEGER :: istep

    ! *** reset time manager and dates -----------------------------------------

    CALL manager_state(echam_time,previous_date)
    CALL manager_init (echam_time,1)         ! increment time step

    ! *** redefine time window -------------------------------------------------

    CALL manager_state(echam_time, istep)
    CALL manager_state(echam_time, current_date)
    CALL manager_state(echam_time, next_date,istep+1)

    ! *** reset switches -------------------------------------------------------

    lstart       = .FALSE.
    lresume      = .FALSE.
    lfirst_cycle = .FALSE.
    lnwp         = .FALSE.

  END SUBROUTINE time_reset


  !+
  ! ------------------------------------------------------------------------------
  !
  ! utility procedures
  !

  ! ----------------------------------------------------------------------------
  !
  SUBROUTINE echam_ev_print (event, lformat)  !*********************************

    ! print out event contents
    !-
    TYPE(time_event), INTENT(in) :: event    ! print it's contents
    LOGICAL, OPTIONAL            :: lformat  ! control the format

    IF (PRESENT(lformat)) THEN
      IF (lformat) THEN
        CALL event_print(event)           ! extensive information
      ELSE
        CALL event_print(event,.TRUE.)    ! one row enhanced
      END IF

    ELSE
      CALL event_print(event,.FALSE.)     ! one row very short

    END IF

  END SUBROUTINE echam_ev_print


  ! ----------------------------------------------------------------------------
  !+
  FUNCTION get_interval_steps (ev) RESULT (steps)  !****************************

    ! returns number of steps between last and present event trigger
    ! the trigger interval is given in echam steps
    !-
    TYPE(time_event) :: ev    ! use triggers from me
    INTEGER          :: steps

    steps = event_state(ev,delta_time)

    IF (ldebugev) THEN
      WRITE(m_text,*) 'return ',steps,&
           ' steps from event <',TRIM(TE_print_event_name(ev)),'>'
      CALL message('',m_text)
    END IF

  END FUNCTION get_interval_steps


  ! ----------------------------------------------------------------------------
  !+
  FUNCTION get_interval_seconds (ev) RESULT (sec)  !****************************

    ! returns number of seconds between last and present event trigger
    ! the trigger interval is given in seconds
    !-

    TYPE(time_event) :: ev    ! use triggers from me
    INTEGER          :: sec

    sec = event_state(ev)

#ifdef DEBUG
    IF (ldebugev) THEN
      WRITE(m_text,*) 'return ',sec,&
           ' seconds from event <',TRIM(TE_print_event_name(ev)),'>'
      CALL message('',m_text)
    END IF
#endif

  END FUNCTION get_interval_seconds

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION get_interval_seconds_next (ev) RESULT (sec)  !***********************

    ! returns number of seconds between present and next event trigger
    ! the trigger interval is given in seconds
    !-

    TYPE(time_event) :: ev    ! use triggers from me
    INTEGER          :: sec
    LOGICAL          :: l_next

    l_next = .TRUE.
    sec = event_state(ev, l_next)

  END FUNCTION get_interval_seconds_next

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE write_date_days (dat_day, text)  !*********************************

    ! write date in constant format to standard output
    ! input date can be different declared
    !-
    TYPE (time_days), INTENT(in)           :: dat_day
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: text
    CHARACTER(len=STR_LEN_B)               :: date_mess1, date_mess2

    CALL print_date(dat_day,PRN_DATE_FORMAT,mess=date_mess1)
    IF (PRESENT(text)) THEN
      !
      ! problem with Linux, cannot copy string to itself
      !
      date_mess2 = TRIM(text) // ' ' // TRIM(date_mess1)
      date_mess1 = TRIM(date_mess2)
    END IF
    CALL message('',date_mess1)

  END SUBROUTINE write_date_days

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE write_date_intern (dat_int, text)   !******************************
    !-
    TYPE (time_intern), INTENT(in)          :: dat_int
    CHARACTER(len=*),  OPTIONAL, INTENT(in) :: text
    CHARACTER(len=STR_LEN_B)                :: date_mess1, date_mess2

    CALL print_date(dat_int,PRN_DATE_FORMAT,mess=date_mess1)
    IF (PRESENT(text)) THEN
      date_mess2 = TRIM(text) // ' ' // TRIM(date_mess1)
      date_mess1 = TRIM(date_mess2)
    END IF
    CALL message('',date_mess1)

  END SUBROUTINE write_date_intern

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE write_date_native (dat_nat, text)   !******************************
    !-
    TYPE (time_native), INTENT(in)         :: dat_nat
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: text
    CHARACTER(len=STR_LEN_B)               :: date_mess1, date_mess2

    CALL print_date(dat_nat,PRN_DATE_FORMAT,mess=date_mess1)
    IF (PRESENT(text)) THEN
      date_mess2 = TRIM(text) // ' ' // TRIM(date_mess1)
      date_mess1 = TRIM(date_mess2)
    END IF
    CALL message('',date_mess1)

  END SUBROUTINE write_date_native


  ! ----------------------------------------------------------------------------
  !+
  REAL(dp) FUNCTION set_delta_time (truncation,nlevel)

    ! preset time stepping (in seconds) dependent on the truncation
    !-
    INTEGER, INTENT(in) :: truncation, nlevel

    IF (nlevel == 19 .OR. nlevel == 11) THEN
      IF      (truncation == 31) THEN; set_delta_time = 2400.0_dp
      ELSE IF (truncation == 42) THEN; set_delta_time = 1800.0_dp
      ELSE
        CALL finish ('set_delta_time', 'Truncation not supported.')
      END IF
    ELSE IF (nlevel == 31) THEN
      IF      (truncation == 31) THEN; set_delta_time = 1800.0_dp
      ELSEIF  (truncation == 63) THEN; set_delta_time =  720.0_dp
      ELSE
        CALL finish ('set_delta_time', 'Truncation not supported.')
      END IF
    ELSE IF (nlevel == 39) THEN
      ! Middle atmosphere version
      IF      (truncation == 31) THEN; set_delta_time =  900.0_dp
      ELSE
       CALL finish ('set_delta_time', 'Truncation not supported.')
      END IF
    ELSE IF (nlevel == 47) THEN
      ! Middle atmosphere version
      IF      (truncation == 31) THEN; set_delta_time =  900.0_dp
      ELSE IF (truncation == 63) THEN; set_delta_time =  600.0_dp
      ELSE IF (truncation ==127) THEN; set_delta_time =  300.0_dp
      ELSE
        CALL finish ('set_delta_time', 'Truncation not supported.')
      END IF
    ELSE IF (nlevel == 95) THEN  
      ! Middle atmosphere version
      IF      (truncation == 63) THEN; set_delta_time =  450.0_dp
      ELSE IF (truncation ==127) THEN; set_delta_time =  240.0_dp
      ELSE IF (truncation ==255) THEN; set_delta_time =  120.0_dp
      ELSE
        CALL finish ('set_delta_time', 'Truncation not supported.')
      END IF
    ELSE
      CALL finish ('set_delta_time', 'Truncation not supported.')
    END IF

  END FUNCTION set_delta_time


  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE set_start_date_days (date) !***************************************
    !-
    TYPE(time_days), INTENT(in) :: date

    CALL manager_init(echam_time,date)

  END SUBROUTINE set_start_date_days


  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE get_start_date_days (date)  !**************************************
    !-
    TYPE(time_days),INTENT(out) :: date

    CALL manager_state(echam_time,date,init_step)

  END SUBROUTINE get_start_date_days

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE set_start_date_intern (my_date)  !*********************************
    !-
    TYPE(time_intern), INTENT(in) :: my_date
    TYPE(time_days)               :: date

    CALL TC_convert(my_date, date)
    CALL set_start_date_days(date)

  END SUBROUTINE set_start_date_intern

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE get_start_date_intern (my_date)  !*********************************
    !-
    TYPE(time_intern), INTENT(out) :: my_date
    TYPE(time_days)                :: date

    CALL get_start_date_days(date)
    CALL TC_convert(date, my_date)

  END SUBROUTINE get_start_date_intern

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE set_start_date_native (my_date)  !*********************************
    !-
    TYPE(time_native), INTENT(in) :: my_date
    TYPE(time_days)               :: date

    CALL TC_convert(my_date, date)
    CALL set_start_date_days(date)

  END SUBROUTINE set_start_date_native

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE get_start_date_native (my_date)   !********************************
    !-
    TYPE(time_native), INTENT(out) :: my_date
    TYPE(time_days)                :: date

    CALL get_start_date_days(date)
    CALL TC_convert(date, my_date)

  END SUBROUTINE get_start_date_native


  SUBROUTINE get_date_components(time, year, month, day, hour, minute, second)
    TYPE (time_days)  ,INTENT(in)   :: time
    INTEGER ,OPTIONAL ,INTENT(out)  :: year, month, day, hour, minute, second

    TYPE (time_native) :: my_time
    INTEGER            :: yr,mo,dy,hr,mn,se

    CALL TC_convert(time,my_time)
    CALL TC_get(my_time,yr,mo,dy,hr,mn,se)

    IF (PRESENT(year  )) year   = yr
    IF (PRESENT(month )) month  = mo
    IF (PRESENT(day   )) day    = dy
    IF (PRESENT(hour  )) hour   = hr
    IF (PRESENT(minute)) minute = mn
    IF (PRESENT(second)) second = se

  END SUBROUTINE get_date_components

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION get_forecast_hours() RESULT (hours)  !*******************************

    !  calculates the time between start and next date in hours
    !-
    INTEGER         :: hours

    TYPE(time_days) :: mydate
    INTEGER         :: iday, isec, isteps
    INTEGER         :: current_day, current_sec
    INTEGER         :: start_day, start_sec

    hours = 0

    CALL event_current_date   (ev_putdata(1), mydate)
    CALL write_date(mydate,'1) ........... ')

    CALL TC_get   (start_date, iday, isec)
    CALL write_date(start_date,'2) ........... ')

    CALL TC_get   (mydate, current_day, current_sec)
    CALL TC_get   (start_date, start_day, start_sec)

    CALL add_date            (-iday,-isec, mydate)
    CALL TC_get                           (mydate, iday, isec)

    isteps = INT((iday*3600+isec)/delta_time)
    isec   = isteps * INT(delta_time)

    write (0,*) iday, isec

    ! convert into hours

    IF (MOD(isec,3600) > 0) &
         CALL finish('mo_time_control:get_forecast_hours',&
         'Forecast time not multiple of hours.')

    hours = INT(isec/3600)

    IF (hours > 744) THEN
      CALL message('','NWP mode makes no sense for this time range.')
      CALL message('','Please change to the climate mode.')
      CALL finish('mo_time_control:get_forecast_hours','Run terminated.')
    END IF

  END FUNCTION get_forecast_hours


  ! ----------------------------------------------------------------------------
  !+
  REAL(dp) FUNCTION get_delta_time()  !*****************************************

    ! generate the time interval from echam_time manager
    !-
    REAL(dp) :: atime

    CALL manager_state(echam_time,atime)
    get_delta_time = atime

  END FUNCTION get_delta_time


  ! ----------------------------------------------------------------------------
  !+
  FUNCTION get_time_step() RESULT (istep)  !************************************

    ! provide time step from echam time manager
    !-
    INTEGER :: istep

    CALL manager_state(echam_time,istep)

  END FUNCTION get_time_step

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION get_month_len (yr, mo) RESULT (mlen)  !******************************

    ! returns the number of days of a month
    !-
    INTEGER            :: yr, mo, mlen
    TYPE (time_native) :: date

    CALL TC_set(yr, mo, 1, 0, 0, 0, date)
    mlen = month_len(date)

  END FUNCTION get_month_len


  ! ----------------------------------------------------------------------------
  !+
  FUNCTION get_year_len (year) RESULT (yearl)  !********************************

    ! returns the number of days in the year
    ! without parameter the length of a mean year is given
    !-
    REAL(dp)          :: yearl
    INTEGER, OPTIONAL :: year

    TYPE(time_native) :: my_date
    TYPE(time_days)   :: my_day

    IF (PRESENT(year)) THEN      ! use the absolut number of days
      CALL TC_set(year,1,1,0,0,0,my_date)
      CALL TC_convert(my_date, my_day)
      yearl = REAL(year_len(my_day),dp)

    ELSE                         ! get the yearl length with fractional day
      yearl = year_len()

    END IF

  END FUNCTION get_year_len


  ! ----------------------------------------------------------------------------
  !+
  FUNCTION get_year_day (date) RESULT (dayno)   !*******************************

    ! returns the day of the year for a given date
    ! the seconds of the day are contained fractional
    !-
    TYPE(time_days) :: date    ! evaluate for this date
    REAL(dp)        :: dayno

    INTEGER   :: yr, mo, dy, hr, mn, se
    INTEGER   :: jyr, jyday, jsec, isec
    TYPE(julian_date) :: juldate
    TYPE(ly360_date) :: ly360date

    CALL get_date_components(date,yr,mo,dy,hr,mn,se)
    isec = IMerge_HMS2Sec(hr, mn, se)

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      CALL Set_JulianDay(yr,mo,dy,isec,juldate)
      CALL Get_JulianYearDay(juldate, jyr,jyday,jsec)
    CASE (CYL360)
      CALL Set_Ly360Day(yr,mo,dy,isec,ly360date)
      CALL Get_Ly360YearDay(ly360date, jyr,jyday,jsec)
    END SELECT

    dayno = REAL(jyday,dp) + sec2frac(jsec)

  END FUNCTION get_year_day


  ! ----------------------------------------------------------------------------
  !+
  FUNCTION str_date_intern (no, date) RESULT (label)  !*************************

    ! generate strings with time information
    !-
    INTEGER           :: no   ! label type
    TYPE(time_intern) :: date ! convert this date
    TYPE(time_native) :: my_date
    CHARACTER(len=19) :: label

    CALL TC_convert(date, my_date)
    label = str_date_native(no, my_date)

  END FUNCTION str_date_intern

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION str_date_days (no, date) RESULT (label)  !***************************
    !-
    INTEGER           :: no   ! label type
    TYPE(time_days)   :: date ! convert this date
    TYPE(time_native) :: my_date
    CHARACTER(len=19) :: label

    CALL TC_convert(date, my_date)
    label = str_date_native(no, my_date)

  END FUNCTION str_date_days

  ! ----------------------------------------------------------------------------
  !+
  FUNCTION str_date_native (no, date) RESULT (label)  !*************************
    !-
    INTEGER           :: no   ! label type
    TYPE(time_native) :: date ! convert this date

    CHARACTER(len=23) :: label, label_string

    INTEGER :: yr, mo, dy, hr, mn, se
    INTEGER :: my_no

    label_string = ''

    my_no = MAX(MIN(no,FLT_NWP_YMD),FLT_YMD)  ! as minimum return year+month+day

    CALL TC_get(date,yr,mo,dy,hr,mn,se)

    SELECT CASE (my_no)
    CASE (FLT_YM)
      WRITE(label_string,'(i8.4,i2.2)')         yr,mo
    CASE (FLT_YMD)
      WRITE(label_string,'(i8.4,i2.2,a1,i2.2)') yr,mo,'.',dy
    CASE (FLT_YMDH)
      WRITE(label_string,'(i8.4,i2.2,2(a1,i2.2))') &
                                               yr,mo,'.',dy,'_',hr
    CASE (FLT_YMDHM)
      WRITE(label_string,'(i8.4,i2.2,3(a1,i2.2))') &
                                               yr,mo,'.',dy,'_',hr,':',mn
    CASE (FLT_YMDHMS)
      WRITE(label_string,'(i8.4,i2.2,4(a1,i2.2))') &
                                               yr,mo,'.',dy,'_',hr,':',mn,':',se
    CASE (FLT_ISO_YMD)
      WRITE(label_string,'(i8.4,a1,i2.2,a1,i2.2)') yr,'-',mo,'-',dy
    CASE (FLT_ISO_YMDHMS)
      WRITE(label_string,'(i8.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') yr,'-',mo,'-',dy, '_',hr,':',mn,':',se
    CASE (FLT_NWP_YMD)
      WRITE(label_string,'(i8.4,i2.2,i2.2)') yr,mo,dy
    CASE default
      WRITE(label_string,'(i8.4,a1,i2.2,a1,i2.2)') yr,'-',mo,'-',dy
    END SELECT

    label=ADJUSTL(label_string)

  END FUNCTION str_date_native


  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE inp_convert_date (d1, d2, date)  !*********************************

    ! convert date and time from initial/rerun files to date
    !-
    INTEGER, INTENT(in)          :: d1, d2   ! YYMMDD, HHMMSS
    TYPE(time_days), INTENT(out) :: date     ! transformed date
    TYPE(time_intern)            :: io_date

    CALL TC_set (d1, d2, io_date)
    CALL TC_convert(io_date, date)

  END SUBROUTINE inp_convert_date

  ! ----------------------------------------------------------------------------
  !+
  SUBROUTINE out_convert_date (date, d1, d2)  !*********************************

    ! convert date and time from date to initial/rerun file date
    !-
    TYPE(time_days),INTENT(in) :: date   ! inserted date
    INTEGER, INTENT(out)       :: d1, d2 ! transformed YYMMDD, HHMMSS
    TYPE(time_intern)          :: io_date

    CALL TC_convert(date,io_date)
    CALL TC_get(io_date,d1,d2)

  END SUBROUTINE out_convert_date

  ! ----------------------------------------------------------------------------
  !+
  !*********************************
  SUBROUTINE day_difference(date_start, date_end, year_start, month_start, day_start, ndays) 
  !*********************************

    ! Convert the ASCII ISO dates date_start & date_end into the respective ymd values and 
    ! calculate the difference in days including date_end
    !
    ! Stefan Hagemann - Helmholtz-Zentrum Hereon - August 2023
    !
    CHARACTER(LEN=10), INTENT(IN)  :: date_start  ! start date, format YYYYMMDD or YYYY-MM-DD
    CHARACTER(LEN=10), INTENT(IN)  :: date_end    ! end date, format YYYYMMDD or YYYY-MM-DD
    INTEGER, INTENT(OUT)           :: year_start  ! start year
    INTEGER, INTENT(OUT)           :: month_start ! start month
    INTEGER, INTENT(OUT)           :: day_start   ! start day
    INTEGER, INTENT(OUT)           :: ndays       ! difference in days including date_end

    INTEGER           :: ymd1, ymd2, day1, day2, sec1, sec2, ios
    INTEGER           :: yyyy, mm, dd
    INTEGER           :: zero = 0
    TYPE(time_days)   :: date1
    TYPE(time_days)   :: date2
    TYPE(time_native) :: date_nat

    !
    ! *** Convert ASCII dates to ymd numbers
    IF (LEN_TRIM(date_start).EQ.8) THEN
      READ(date_start,'(I8)', iostat=ios) ymd1
      IF (ios.NE.0) THEN
        WRITE (message_text,*) 'Error while reading date_start --> checkISO-format: either YYYYMMDD or YYYY-MM-DD'
        CALL finish ('config_hd', message_text)
      ENDIF
      CALL inp_convert_date (ymd1, zero, date1)
      CALL TC_convert(date1, date_nat)
      CALL TC_get    (date_nat, year_start, month_start, day_start)
    ELSE IF (LEN_TRIM(date_start).EQ.10) THEN
      READ(date_start,'(I4,1X,I2,1X,I2)', iostat=ios) year_start, month_start, day_start
      IF (ios.NE.0) THEN
        WRITE (message_text,*) 'Error while reading date_start --> checkISO-format: either YYYY-MM-DD or YYYYMMDD'
        CALL finish ('config_hd', message_text)
      ENDIF
      ymd1 = 10000 * year_start + 100 * month_start + day_start
      CALL inp_convert_date (ymd1, zero, date1)
    ELSE 
       WRITE (message_text,*) 'date_start has no ISO-format, neither YYYYMMDD nor YYYY-MM-DD -> Error'
       CALL finish ('config_hd', message_text)
    ENDIF	
!
    IF (LEN_TRIM(date_end).EQ.8) THEN
      READ(date_end,'(I8)', iostat=ios) ymd2
      IF (ios.NE.0) THEN
        WRITE (message_text,*) 'Error while reading date_end --> checkISO-format: either YYYYMMDD or YYYY-MM-DD'
        CALL finish ('config_hd', message_text)
      ENDIF
    ELSE IF (LEN_TRIM(date_end).EQ.10) THEN
      READ(date_end,'(I4,1X,I2,1X,I2)', iostat=ios) yyyy, mm, dd
      IF (ios.NE.0) THEN
        WRITE (message_text,*) 'Error while reading date_end --> checkISO-format: either YYYY-MM-DD or YYYYMMDD'
        CALL finish ('config_hd', message_text)
      ENDIF
      ymd2 = 10000 * yyyy + 100 * mm + dd
    ELSE 
       WRITE (message_text,*) 'date_end has no ISO-format, neither YYYYMMDD nor YYYY-MM-DD -> Error'
       CALL finish ('config_hd', message_text)
    ENDIF	
!
!   *** Calculate difference in days
    CALL inp_convert_date (ymd2, zero, date2)

    ! Calculate days from ymd1 to ymd2 incl. ymd2 
    CALL TC_get(date1, day1, sec1)     
    CALL TC_get(date2, day2, sec2)     
    ndays = day2 - day1 + 1 + (sec2-sec1)/86400._dp

  END SUBROUTINE day_difference

  !+
  ! ------------------------------------------------------------------------------
  !
  ! functions needed for orbit
  !
  !---------------------------------------------------------------------------
  !>
  !! GET_ORBIT_TIMES: Returns date for orbit model and day fraction
  !!
  !! @par Description
  !! The routine originates from a routine in ECHAM5 called prerad.  The
  !! purpose is to return the day fraction and what we call the orbit_day. The
  !! orbit_day must be cast in a format appropriate for the orbit model being
  !! used, and consistent with control flags allowing for things like a 360 d
  !! year, or perpetual month runs, and or a perpetual year run.  Currently
  !! two orbital models are being used, one (vsop87) requires the Julian date,
  !! the other requires days since the vernal equinox measured in radians.
  !! The procedure calculates the diurnal fraction and orbit day for either
  !! the current time or the radiation time, the latter is the current time
  !! plus half of the interval to the next radiation step and is used when
  !! radiation is not updated at every timestep. So doing ensures consistency
  !! with calculated values and those later used in "radheat"
  !!
  !! @par Revsision History
  !! Abstracted from prerad by B Stevens (2009-08)
  !!
  SUBROUTINE get_orbit_times(lrad_date, lyr_perp, nmonth, yr_perp            &
       , time_of_day, orbit_date)

    LOGICAL, INTENT (IN)    :: lrad_date, lyr_perp
    INTEGER, INTENT (IN)    :: nmonth, yr_perp
    REAL (dp), INTENT (OUT) :: time_of_day, orbit_date

    TYPE(julian_date) :: date_now, date_pal
    TYPE(ly360_date)  :: idate_format
    TYPE(time_days)   :: valid_date

    INTEGER  :: iyr, imo, idy, isec, jsec
    REAL(dp) :: zdy, zdy_mar0, zscr

    if (lrad_date) then
      valid_date = radiation_date
    else
      valid_date = current_date
    end if

    CALL get_date_components(valid_date, year=iyr, month=imo, day=idy)
    CALL TC_get(valid_date, second=jsec)
    time_of_day = (REAL(jsec, dp)/day_len())*2.0_dp*api
    !
    ! Calculate orbital model input for a real orbit, with the possibility
    ! of a perpetual year, as determined by (lyr_perp, yr_perp)
    ! --------------------------------
    IF (l_orbvsop87) THEN
      if (lyr_perp) iyr = yr_perp  ! use the specified yr for perptual yr case
      CALL Set_JulianDay(iyr, imo, idy, jsec, date_now, lperpetual_year=lyr_perp)
      orbit_date = date_now%day + date_now%fraction
      !
      ! Calculate orbital model imput for an idealized orbit with exception
      ! handling. Exceptions include:  perpetual month experiments (where the
      ! orbital parameters are fixed on the middle point of the month), and an
      ! artificial 360 day calendar.  For this orbital model the day must be
      ! converted to days from vernal equinox, and to conform with CMIP specs.
      ! orbital positions are based on the days elapsed since 1900-01-01.
      ! --------------------------------
    ELSE
      SELECT CASE (get_calendar_type())
      CASE (JULIAN)
        CALL Set_JulianDay(1900, 1, 1, 0, date_pal)
        IF (nmonth /= 0) then
          idy  = get_month_len(1987,nmonth)
          isec = INT(MOD(idy,2)*IDAYLEN*0.5_dp)
          idy  = idy/2 + 1
          CALL Set_JulianDay(1987, nmonth, idy, isec, date_now)
        ELSE
          CALL Set_JulianDay(iyr, imo, idy, jsec, date_now)
        END IF
      CASE (CYL360)
        CALL Set_Ly360Day(1900, 1,   1,     0, idate_format)
        date_pal%day      = REAL(idate_format%day,dp)
        date_pal%fraction = idate_format%fraction
        CALL Set_Ly360Day(iyr, imo, idy, jsec, idate_format)
        date_now%day      = REAL(idate_format%day,dp)
        date_now%fraction = idate_format%fraction
      END SELECT
      zdy = (date_now%day+date_now%fraction)-(date_pal%day+date_pal%fraction)
      !
      ! Here is we convert to days since vernal equinox
      ! --------------------------------
      zdy_mar0   = 78.41_dp - 0.0078_dp*(1900-1987) + 0.25_dp*MOD(1900,4)
      zscr       = zdy + get_year_len() - zdy_mar0
      orbit_date = MOD(zscr/get_year_len(),1.0_dp)*2.0_dp*api
    END IF

  END SUBROUTINE get_orbit_times

  SUBROUTINE radiation_time

    INTEGER       :: iradhlen, iradlen
      radiation_date = current_date
!++jsr+mag assure that radiation is calculated in the middle of two time steps
!  -0.5*delta_time was added
      iradhlen = INT(0.5_dp*(event_state(ev_trigrad,.TRUE.)-delta_time))
!--jsr-mag
      CALL add_date(0,iradhlen,radiation_date)
      prev_radiation_date=radiation_date
      iradlen=-event_state(ev_trigrad,.TRUE.)
      CALL add_date(0,iradlen,prev_radiation_date)
    END SUBROUTINE radiation_time

  ! ----------------------------------------------------------------------------
  !
  FUNCTION get_clock (date) RESULT (zclock)  !**********************************

    ! calculates the daytime from present time or radiation time
    !-
    TYPE(time_days)  :: date    ! get clock from this date
    REAL(dp)         :: zclock

    INTEGER :: jday, jsec

    CALL TC_get (date, jday, jsec)
    zclock = REAL(jsec, dp)/day_len()*2*api

  END FUNCTION get_clock

  FUNCTION get_cycle() RESULT (ncycle)
    INTEGER :: ncycle
    ncycle = nmcount
  END FUNCTION get_cycle

END MODULE mo_time_control
!
! ------------------------------------------------------------------------------
