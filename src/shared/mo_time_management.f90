! This module contains routines to set start date, end date etc.
!
! These routines are necessary, since (for reasons of backward
! compatibility) different naemlist settings may be used to specify
! this date information.
!
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

MODULE mo_time_management

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int32_t
  USE mo_kind,                     ONLY: wp, i8
  USE mo_master_config,            ONLY: lrestart_write_last
  USE mo_parallel_config,          ONLY: num_restart_procs
  USE mo_util_string,              ONLY: tolower, int2string
  USE mtime,                       ONLY: MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN,      &
    &                                    MAX_CALENDAR_STR_LEN,                             &
    &                                    OPERATOR(>),OPERATOR(/=), OPERATOR(-),            &
    &                                    datetime, newDatetime, deallocateDatetime,        &
    &                                    timedelta, newTimedelta, deallocateTimedelta,     &
    &                                    getPTStringFromMS, getPTStringFromSeconds, min,   &
    &                                    deallocateTimedelta, OPERATOR(+), OPERATOR(==),   &
    &                                    timedeltatostring, datetimetostring,              &
    &                                    OPERATOR(*), OPERATOR(<), OPERATOR(<=),           &
    &                                    calendarToString, divideTimeDeltaInSeconds,       &
    &                                    divisionquotienttimespan,                         &
    &                                    mtime_proleptic_gregorian => proleptic_gregorian, &
    &                                    mtime_year_of_365_days => year_of_365_days,       &
    &                                    mtime_year_of_360_days => year_of_360_days,       &
    &                                    getTotalMilliSecondsTimeDelta
  USE mo_time_config,              ONLY: dt_restart, time_config,                          &
    &                                    ini_datetime_string, is_relative_time,            &
    &                                    time_nml_icalendar => icalendar,                  &
    &                                    restart_calendar, restart_ini_datetime_string,    &
    &                                    set_calendar, set_is_relative_time,               &
    &                                    set_tc_dt_model, calendar_index2string
  USE mo_run_config,               ONLY: dtime, mtime_modelTimeStep => modelTimeStep
  USE mo_master_control,           ONLY: atmo_process, get_my_process_type
  USE mo_impl_constants,           ONLY: max_dom, inh_atmosphere,                          &
    &                                    dtime_proleptic_gregorian => proleptic_gregorian, &
    &                                    dtime_cly360              => cly360,              &
    &                                    dtime_julian_gregorian    => julian_gregorian
  USE mo_dynamics_config,          ONLY: iequations
  USE mo_exception,                ONLY: message, message_text, finish
  USE mo_grid_config,              ONLY: patch_weight, grid_rescale_factor, n_dom,         &
    &                                    start_time, lrescale_timestep
  USE mo_io_config,                ONLY: dt_checkpoint
  USE mo_master_config,            ONLY: experimentReferenceDate,                          &
    &                                    experimentStartDate,                              &
    &                                    checkpointTimeIntval, restartTimeIntval,          &
    &                                    experimentStopDate, isRestart,                    &
    &                                    master_nml_calendar => calendar_str
  USE mo_time_config,              ONLY: set_tc_exp_refdate, set_tc_exp_startdate,         &
    &                                    set_tc_exp_stopdate, set_tc_startdate,            &
    &                                    set_tc_stopdate,                                  &
    &                                    set_tc_dt_restart, set_tc_dt_checkpoint,          &
    &                                    set_tc_current_date, set_tc_write_restart,        &
    &                                    end_datetime_string
  USE mo_restart_nml_and_att,      ONLY: getAttributesForRestarting
  USE mo_key_value_store,          ONLY: t_key_value_store

#ifndef __NO_ICON_ATMO__
  USE mo_nonhydrostatic_config,    ONLY: divdamp_order
  USE mo_atm_phy_nwp_config,       ONLY: atm_phy_nwp_config
  USE mo_initicon_config,          ONLY: timeshift 
#endif


  IMPLICIT NONE

  PUBLIC :: compute_timestep_settings
  PUBLIC :: compute_restart_settings
  PUBLIC :: compute_date_settings  

  PRIVATE

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_time_management'

  LOGICAL :: ldebug = .FALSE.
  
CONTAINS 

  !---------------------------------------------------------------------------------------
  !> Set time step for this run.
  !
  !  This routine is necessary, since (for reasons of backward
  !  compatibility) different naemlist settings may be used to specify
  !  the date information.
  !
  !  Initial revision:  10/2015 : F. Prill, DWD
  !
  SUBROUTINE compute_timestep_settings()
    ! local variables
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::compute_timestep_settings'

    TYPE(timedelta), POINTER             ::  dtime1, dtime2, zero_dt
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) ::  dtime_str2, dtime_string
    INTEGER(i8)                          ::  dtime_ms
    REAL(wp)                             ::  dtime_real
    INTEGER                              ::  jg
    TYPE(datetime), POINTER              ::  reference_dt

    ! --------------------------------------------------------------
    ! PART I: Collect time step as ISO8601 string
    ! --------------------------------------------------------------

    ! The time step is determined either by the "dtime" parameter or
    ! by the "modelTimeStep" parameter from the namelist "run_nml".
    !
    dtime_string = mtime_modelTimeStep
    
    dtime_real = dtime
    IF (dtime_real > 0._wp) THEN
      dtime_ms   = NINT(dtime_real*1000, i8)
      CALL getPTStringFromMS(dtime_ms, dtime_str2)
      
      IF (dtime_string == "") THEN
        dtime_string = dtime_str2
      ELSE
        ! Obviously, both namelist parameters for the time step have
        ! been used. We need to test for equality.
        dtime1 => newTimedelta(dtime_string)
        dtime2 => newTimedelta(dtime_str2)
        IF (dtime1 /= dtime2) THEN
          CALL finish(routine, "Two inconsistent time steps have been specified!")
        END IF
        CALL deallocateTimedelta(dtime1)
        CALL deallocateTimedelta(dtime2)
      END IF
    END IF
    dtime1 => newTimedelta(dtime_string)

    ! Furthermore, the time step may be rescaled by the
    ! "grid_rescale_factor".
    !
    IF (lrescale_timestep .AND. grid_rescale_factor /= 1.0_wp) THEN
      dtime1 = dtime1 * grid_rescale_factor
      CALL timedeltaToString(dtime1, dtime_string)
      IF (dtime_real > 0._wp)  dtime_real = dtime_real * grid_rescale_factor

#ifndef __NO_ICON_ATMO__
      IF (get_my_process_type() == atmo_process) THEN
        DO jg=1,max_dom
          atm_phy_nwp_config(jg)%dt_conv = &
            atm_phy_nwp_config(jg)%dt_conv * grid_rescale_factor
          atm_phy_nwp_config(jg)%dt_rad  = &
            atm_phy_nwp_config(jg)%dt_rad  * grid_rescale_factor
          atm_phy_nwp_config(jg)%dt_sso  = &
            atm_phy_nwp_config(jg)%dt_sso  * grid_rescale_factor
          atm_phy_nwp_config(jg)%dt_gwd  = &
            atm_phy_nwp_config(jg)%dt_gwd  * grid_rescale_factor
        END DO
      END IF
#endif

    ELSE
      CALL timedeltaToString(dtime1, dtime_string)
    END IF

    ! consistency check
    zero_dt => newTimedelta("PT0S") ! mtime object for zero
    IF (dtime1 <= zero_dt) CALL finish(routine,'"dtime" must be positive')

    CALL deallocateTimedelta( zero_dt )

    ! --------------------------------------------------------------
    ! PART II: Convert ISO8601 string into "mtime" and old REAL
    ! --------------------------------------------------------------

    CALL set_tc_dt_model(dtime_string) ! dyn. time step  on the global grid
    IF (dtime_real > 0._wp) THEN
      ! In case that we came from the REAL-valued namelist setting of
      ! the time step we try to avoid rounding errors in floating
      ! point arithmetic
      dtime = dtime_real
    ELSE
      ! For conversion to milliseconds, we need an anchor date,
      ! although we know that "dtime" is too small for this to be
      ! relevant.
      reference_dt => newDatetime("1980-06-01T00:00:00.000")
      dtime = REAL(getTotalMilliSecondsTimeDelta(dtime1, reference_dt),wp)/1000._wp
      CALL deallocateDatetime(reference_dt)
    END IF

    CALL deallocateTimedelta( dtime1  )
    
    ! --------------------------------------------------------------
    ! PART III: Print time step
    ! --------------------------------------------------------------

    CALL message('','')
    WRITE(message_text,'(a,a)') 'Model time step          : ', TRIM(dtime_string)
    CALL message('',message_text)
    CALL message('','')

  END SUBROUTINE compute_timestep_settings

  !---------------------------------------------------------------------------------------
  !> Set restart and checkpoint interval for this run.
  !
  !  This routine is necessary, since (for reasons of backward
  !  compatibility) different naemlist settings may be used to specify
  !  the date information.
  !
  !  Initial revision:  10/2015 : F. Prill, DWD
  !
  SUBROUTINE compute_restart_settings()
    ! local variables:
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::compute_restart_settings'
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: checkpt_intvl_string, checkpt_intvl2,       &
                                            restart_intvl_string, dtime_string
    TYPE(timedelta), POINTER             :: mtime_0h, mtime_2_5h, mtime_dt_checkpoint,  &
      &                                     mtime_dt_restart, mtime_dom_start,          &
      &                                     tmp_td1, tmp_td2
    TYPE(datetime), POINTER              :: reference_dt
    INTEGER                              :: jg
    
    ! --------------------------------------------------------------
    ! PART I: Collect the restart and checkpoint intervals 
    !         as ISO8601 strings
    ! --------------------------------------------------------------

    ! --- --- RESTART INTERVAL:
    !
    !         This time interval specifies when the run is supposed to
    !         save its state to a file and stop, later to be resumed.
    !
    !         The restart interval is set by the namelist parameter
    !         "restartTimeIntval" (in the namelist
    !         "master_time_control_nml") - the new way -, and
    !         "dt_restart" (in "time_nml") - the old way.
    !
    !         TODO: The restart interval needs to be multiple of the model
    !               time steps.
    !
    CALL set_tc_write_restart(.TRUE.)          
    IF (.NOT. lrestart_write_last) THEN
      restartTimeIntval = 'PT0.000S'
      dt_restart = 0.0_wp
      CALL set_tc_write_restart(.FALSE.)      
      restart_intvl_string = restartTimeIntval
    ELSE
      IF (restartTimeIntval /= "") THEN
        IF (dt_restart > 0.0_wp) THEN
          ! Comparison of intervals cannot be done with strings as there are multiple options
          ! for expressing an interval
          tmp_td1 => newTimedelta(restartTimeIntval)
          CALL getPTStringFromSeconds(dt_restart, restart_intvl_string)
          tmp_td2 => newTimedelta(restart_intvl_string)
          IF (tmp_td1 == tmp_td2) THEN
            restart_intvl_string = restartTimeIntval
          ELSE
            ! if both are set but inconsistent, finish with error message
            CALL finish(routine, "Inconsistent setting of restart interval: " &
                 &               //TRIM(restart_intvl_string)//"/"//TRIM(restartTimeIntval))
          ENDIF
          CALL deallocateTimedelta(tmp_td1)
          CALL deallocateTimedelta(tmp_td2)
        ELSE
          ! use restartTimeIntval
          restart_intvl_string = restartTimeIntval
        END IF
      ELSE
        IF (dt_restart > 0.0_wp) THEN
          ! use dt_restart
          CALL getPTStringFromSeconds(dt_restart, restart_intvl_string)
        ELSE
          restart_intvl_string = 'PT0.000S'
          CALL set_tc_write_restart(.FALSE.)
        END IF
      END IF
    END IF

    mtime_dt_restart    => newTimedelta(restart_intvl_string)
    
    ! --- --- CHECKPOINT INTERVAL:
    !
    !         This time interval specifies when the run is supposed to
    !         save its state to a file (but not to stop afterwards).
    !
    !         The checkpoint interval is set by the namelist parameter
    !         "checkpointTimeIntval" (in the namelist
    !         "master_time_control_nml") - the new way -, and
    !         "dt_checkpoint" (in "io_nml") - the old way.
    !
    !         TODO: The checkpoint interval needs to be multiple of the
    !               model time steps.
    !
    checkpt_intvl_string = checkpointTimeIntval
    IF (dt_checkpoint > 0._wp) THEN
      WRITE (checkpt_intvl2, '(a,i0,a)') "PT", INT(dt_checkpoint), "S"
      IF (checkpt_intvl_string == "") THEN
        checkpt_intvl_string = checkpt_intvl2
      ELSE
        tmp_td1 => newTimedelta(checkpt_intvl_string)
        tmp_td2 => newTimedelta(checkpt_intvl2)        
        IF (.NOT. (tmp_td1 < tmp_td2) .AND. .NOT. (tmp_td2 < tmp_td1)) THEN
          checkpt_intvl_string = checkpt_intvl2
        ELSE
          CALL finish(routine, "Inconsistent setting of checkpoint interval: " &
               &               //TRIM(checkpt_intvl_string)//"/"//TRIM(checkpt_intvl2))
        END IF
        CALL deallocateTimedelta(tmp_td1)
        CALL deallocateTimedelta(tmp_td2)
      END IF
    END IF
    ! if "checkpt_intvl_string" still unspecified: set default to no checkpoint
    IF (checkpt_intvl_string == "") THEN
      checkpt_intvl_string = 'PT0.000S'
    END IF

    mtime_dt_checkpoint => newTimedelta(checkpt_intvl_string)

    ! --------------------------------------------------------------
    ! PART II: Convert ISO8601 string into "mtime" and old REAL
    ! --------------------------------------------------------------

    CALL set_tc_dt_restart     ( restart_intvl_string )
    CALL set_tc_dt_checkpoint  ( checkpt_intvl_string )

    ! For conversion to (milli-)seconds, we need an anchor date:
    !
    reference_dt => newDatetime("1980-06-01T00:00:00.000")
    dt_restart    = REAL(getTotalMilliSecondsTimeDelta(mtime_dt_restart,    reference_dt),wp)/1000._wp
    dt_checkpoint = REAL(getTotalMilliSecondsTimeDelta(mtime_dt_checkpoint, reference_dt),wp)/1000._wp
    CALL deallocateDatetime(reference_dt)

    ! --------------------------------------------------------------
    ! PART III: Print restart and checkpoint intervals
    ! --------------------------------------------------------------

    CALL message('','')    
    WRITE(message_text,'(a,a)') 'Checkpoint interval      : ', checkpt_intvl_string
    CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'Restart interval         : ', restart_intvl_string
    CALL message('',message_text)
    CALL message('','')

    ! --------------------------------------------------------------
    ! PART IV:  Consistency checks 
    ! --------------------------------------------------------------
    !
    ! When increased sound-wave and gravity-wave damping is chosen
    ! during the spinup phase (i.e. divdamp_order = 24),
    ! checkpointing/restarting is not allowed earlier than three hours
    ! into the integration because the results would not be
    ! bit-identical in this case
    !

#ifndef __NO_ICON_ATMO__
    mtime_0h => newTimedelta("PT0S")
    IF (mtime_dt_checkpoint /= mtime_0h) THEN
      mtime_2_5h => newTimedelta("PT02H30M")
      IF ((iequations == inh_atmosphere) .AND. &
        & (divdamp_order == 24) .AND. .NOT. isRestart() .AND. &
        & (mtime_dt_checkpoint < mtime_2_5h)) THEN
        WRITE(message_text,'(a)') &
             &  'dt_checkpoint < 2.5 hours not allowed in combination with divdamp_order = 24'
        CALL finish(routine, message_text)
      ENDIF
      CALL deallocateTimedelta(mtime_2_5h)
    ENDIF
    CALL deallocateTimedelta(mtime_0h)    
#endif

    ! Writing a checkpoint file exactly at the start time of a nest is
    ! not allowed:
    !
    DO jg =1,n_dom
      IF (INT(start_time(jg)) <= 0)  CYCLE
      WRITE (dtime_string, '(a,i0,a)') "PT", INT(start_time(jg)), "S"
      mtime_dom_start => newTimedelta(dtime_string)
      IF (mtime_dom_start == mtime_dt_checkpoint) THEN
        WRITE(message_text,'(a)') &
          &  'writing a checkpoint file exactly at the start time of a nest is not allowed'
        CALL finish(routine, message_text)
      END IF
      CALL deallocateTimedelta(mtime_dom_start)
    END DO

    CALL deallocateTimedelta(mtime_dt_checkpoint)
    CALL deallocateTimedelta(mtime_dt_restart)

  END SUBROUTINE compute_restart_settings


  !---------------------------------------------------------------------------------------
  !> Set start date, end date etc. for this run.
  !
  !  This routine is necessary, since (for reasons of backward
  !  compatibility) different naemlist settings may be used to specify
  !  this date information.
  !
  !  This subroutine checks for contradictory namelist settings and
  !  finally creates 
  !   - data objects of type "t_datetime" 
  !   - data objects from  the mtime library.
  !   - the "nsteps" time loop count
  !
  !  As input this subroutine requires final settings dtime and
  !  dt_restart (these may have been modified in the crosscheck
  !  routine in "mo_nml_crosscheck").
  !
  !  Initial revision:  10/2015 : F. Prill, DWD
  !
  SUBROUTINE compute_date_settings(model_string, dt_restart, nsteps)
    CHARACTER(LEN=*),    INTENT(IN)    :: model_string       !< e.g.  "atm"
    REAL(wp),            INTENT(IN)    :: dt_restart         !< Length of restart cycle in seconds
    INTEGER,             INTENT(INOUT) :: nsteps
    ! local variables
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::compute_date_settings'

    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   ::  ini_datetime1, end_datetime1, dstring
    CHARACTER(LEN=32)                     ::  ini_datetime2, end_datetime2
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   ::   &
      &                                       stop_datetime_string,          & !< run stop date
      &                                       exp_start_datetime_string,     & !< experiment start date
      &                                       exp_stop_datetime_string,      & !< experiment stop date
      &                                       exp_ref_datetime_string,       & !< experiment reference date
      &                                       cur_datetime_string
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  ::  td_string
    TYPE(datetime),  POINTER              ::  mtime_start, mtime_stop,       &
      &                                       mtime_exp_stop,                &
      &                                       mtime_restart_stop,            &
      &                                       mtime_nsteps_stop,             &
      &                                       tmp_dt1, tmp_dt2 
    TYPE(timedelta), POINTER              ::  mtime_dt_restart, mtime_dtime, &
      &                                       mtime_td
    TYPE(divisionquotienttimespan)        ::  mtime_quotient
    INTEGER                               ::  mtime_calendar, dtime_calendar,&
      &                                       errno, tlen1, tlen2
    CHARACTER(len=MAX_CALENDAR_STR_LEN)   ::  calendar1, calendar2, calendar
    TYPE(t_key_value_store), POINTER ::  restartAttributes
    CHARACTER(LEN=:), ALLOCATABLE :: start_datetime_string !< run start date
#ifndef __NO_ICON_ATMO__
    REAL(wp)                              :: zdt_shift            ! rounded dt_shift
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: dt_shift_string
#endif

    ! --------------------------------------------------------------
    ! PART I: Collect all the dates as ISO8601 strings
    ! --------------------------------------------------------------

    ! --- Prelimary: set calendars in old and new date module
    !
    !     The calendar may be set by the namelist parameter "calendar"
    !     in "time_nml" or the namelist parameter "calendar" in
    !     "master_time_control_nml".

    ! Convert the calendar setting (which is an integer value for this
    ! namelist) into a string. The naming scheme is then compatible
    ! with concurrent namelist settings of the calendar (mtime):
    calendar1 = calendar_index2string(time_nml_icalendar)
    calendar2 = master_nml_calendar
    tlen1 = LEN_TRIM(calendar1)
    tlen2 = LEN_TRIM(calendar2)
    IF (tlen2 /= 0) THEN
      calendar = calendar2
      IF (tlen1 /= 0) THEN
        ! both settings were used; we need to test for equality
        IF (tolower(calendar1) /= tolower(calendar2))  &
             &  CALL finish(routine, "Inconsistent setting of calendar")
      END IF
    ELSE IF (tlen1 /= 0) THEN
      calendar = calendar1
    ELSE
      calendar = ''
    END IF
    SELECT CASE (toLower(calendar))
    CASE ('julian gregorian')
      dtime_calendar  = dtime_julian_gregorian
      mtime_calendar  = -1
      CALL finish(routine, "Julian-Gregorian calendar unsupported by mtime library!")
    CASE ('proleptic gregorian')
      dtime_calendar  = dtime_proleptic_gregorian
      mtime_calendar  = mtime_proleptic_gregorian
    CASE ('365 day year')  
      dtime_calendar  = -1
      mtime_calendar  = mtime_year_of_365_days
      !      CALL finish(routine, "Year-of-365-days calendar unsupported by datetime module!")
    CASE ('360 day year')  
      dtime_calendar  = dtime_cly360
      mtime_calendar  = mtime_year_of_360_days
    CASE default
      dtime_calendar       = dtime_proleptic_gregorian
      mtime_calendar       = mtime_proleptic_gregorian
      time_nml_icalendar   = dtime_proleptic_gregorian
      CALL message('','No calendar selected! Use default proleptic Gregorian.')
    END SELECT

    ! --- Now, we merge ISO time stamp strings from the namelists
    !     "time_nml" and "master_time_control_nml". If both namelist
    !     specify the same setting, then we compare the dates and
    !     throw an error in case of inconsistency.  We take also into
    !     account the namelist settings for "nsteps" and "dt_restart".
    !
    !    ***************************  EXPERIMENT RUN   *************
    !    *                                                         *
    !    *                    ***  JOB RUN   ***                   *
    !    *                    *                *                   *
    !  --|--------------------[----------------]-------------------|---------------> (time axis)
    !   exp_start           start            stop                exp_stop
    !                         |----------------|
    !                          restart interval
    !
    ! --- --- EXPERIMENT START DATE:
    !
    !         This is the start date for the whole experiment, which
    !         may be much earlier than the start of the current run.
    !
    !         The experiment start date is set by the namelist
    !         parameter "experimentStartDate" (in namelist
    !         "master_time_control_nml") - the new way -, and
    !         "ini_datetime_string" (in "time_nml") - the old way.
    !
    !         Note that if the namelist parameter
    !         "experimentStartDate" (in "master_time_control_nml") is
    !         unspecified, but the namelist parameter
    !         "experimentReferenceDate" (see also below) is set, then
    !         set the experiment start date to the latter.
    !
    exp_start_datetime_string = ''
    ini_datetime1 = experimentStartDate
    ini_datetime2 = ini_datetime_string
    IF (ini_datetime1 /= "")  exp_start_datetime_string = ini_datetime1
    IF (ini_datetime2 /= "")  exp_start_datetime_string = ini_datetime2
    IF (ini_datetime1 /= "" .AND. LEN_TRIM(ini_datetime2) > 0) THEN
      ! both settings were used; we need to test for equality
      tmp_dt1 => newDatetime(ini_datetime1)
      tmp_dt2 => newDatetime(ini_datetime2)      
      IF (.NOT. (tmp_dt1 == tmp_dt2)) THEN
        CALL finish(routine, "Inconsistent setting of experiment start date: " &
             &               //TRIM(ini_datetime1)//"/"//TRIM(ini_datetime2))
      END IF
      CALL deallocateDatetime(tmp_dt1)
      CALL deallocateDatetime(tmp_dt2)      
    END IF
    IF ( (ini_datetime1 == "")        .AND.  &
      &  (ini_datetime2 == "")        .AND.  &
      &  (experimentReferenceDate /= "")) THEN
      ! set the experiment start date to the reference date
      exp_start_datetime_string = experimentReferenceDate
    END IF
    ! throw an error, if no start date has been specified at all
    IF (exp_start_datetime_string == "") THEN
      CALL finish(routine, "No experiment start date has been set!")
    END IF


    ! --- --- EXPERIMENT STOP DATE:
    !
    !         This is when the whole experiment ends; note that the
    !         current run may stop before this end date has been
    !         reached (and be restarted later).
    !
    !         The experiment stop date is set by the namelist
    !         parameter "experimentStopDate" (in namelist
    !         "master_time_control_nml") - the new way -, and
    !         "end_datetime_string" (in "time_nml") - the old way.
    !
    exp_stop_datetime_string = ''
    end_datetime1 = experimentStopDate
    end_datetime2 = end_datetime_string
    IF (end_datetime1 /= "")  exp_stop_datetime_string = end_datetime1
    IF (end_datetime2 /= "")  exp_stop_datetime_string = end_datetime2
    IF (end_datetime1 /= "" .AND. end_datetime2 /= "") THEN
      ! both settings were used; we need to test for equality
      tmp_dt1 => newDatetime(end_datetime1)
      tmp_dt2 => newDatetime(end_datetime2)      
      IF (.NOT. (tmp_dt1 == tmp_dt2)) THEN
        CALL finish(routine, "Inconsistent setting of experiment stop date: " &
             &               //TRIM(end_datetime1)//"/"//TRIM(end_datetime2))
      END IF
      CALL deallocateDatetime(tmp_dt1)
      CALL deallocateDatetime(tmp_dt2)      
    END IF
    ! throw an error, if no start date has been specified at all
    !
    ! Note: this error check is disabled since some scripts do not
    !       specify the experiment stop date.
    !
    !  IF (TRIM(exp_stop_datetime_string) == "") THEN
    !    CALL finish(routine, "No experiment stop date has been set!")
    !  END IF

    ! --- --- REFERENCE DATE:
    !
    !         This specifies the reference date for the calendar in
    !         use (some kind of "anchor date" on the time line).
    !
    !         If the namelist parameter "experimentReferenceDate (in
    !         "master_time_control_nml") is unspecified, then the
    !         reference date is set to the experiment start date.
    !
    IF (experimentReferenceDate == '') THEN
      exp_ref_datetime_string = exp_start_datetime_string
    ELSE
      exp_ref_datetime_string = experimentReferenceDate
    END IF

    ! --- --- START DATE:
    !
    !         This is the start date of the current run.
    !
    !         a) in case of restart:  start date := restart attribute "tc_startdate"
    !         b) no-restart case:     start date := experiment start date
    !
    !         Special case: In a restart situation, if the calendar or
    !         initial date/time is different from those in the restart
    !         file, we regard this integration as a new one with its
    !         own calendar.  Model time at which the previous run
    !         stopped is thus not relevant.  Simulation will start
    !         from the user-specified initial date/time, which is also
    !         the current model date/time.
    !
    IF (isRestart()) THEN

      IF ((restart_calendar /= dtime_calendar) .AND. &
        & (restart_calendar /= -1)) THEN
        
        CALL message('','Restart calendar is not matching, fallback to experiment start date: '//&
          &int2string(restart_calendar)//' /= '//int2string(dtime_calendar))
        start_datetime_string = exp_start_datetime_string

      ELSE IF (restart_ini_datetime_string /= ini_datetime_string) THEN

        CALL message('','Restart ini date is not matching, fallback to experiment start date: '//&
        &TRIM(restart_ini_datetime_string)//' /= '//TRIM(ini_datetime_string))

      ELSE
        CALL message('','Read restart file meta data ...')
        CALL getAttributesForRestarting(restartAttributes)
        IF (restartAttributes%is_init) THEN
          CALL restartAttributes%get('tc_startdate', start_datetime_string)
        ELSE
          CALL finish(routine, "Could not retrieve tc_startdate from restart file!")
        ENDIF
      END IF

    ELSE
      CALL message('','This is not a RESTART run ...')
      start_datetime_string = exp_start_datetime_string
    ENDIF

    ! --- --- CURRENT DATE:
    !
    !         This is the current model date, that is stepped forward
    !         during the integration loop. We begin by setting it to
    !         the start date computed above.
    !
    cur_datetime_string = start_datetime_string


    IF (model_string == 'atm') THEN
#ifndef __NO_ICON_ATMO__
      !
      ! timeshift-operations for CURRENT DATE
      !
      ! A timeshift will be used to shift the current model date, and thus 
      ! the actual start date by timeshift%dt_shift backwards in time. 
      ! This is required for the Incremental Analysis Update (IAU) procedure 
      ! which is used during the initialization phase of the atmospheric 
      ! model component, in order to filter spurious noise.
      !
      !
      ! Round dt_shift to the nearest integer multiple of the advection time step
      !
      IF (timeshift%dt_shift < 0._wp) THEN
        zdt_shift = REAL(NINT(timeshift%dt_shift/dtime),wp)*dtime
        IF (ABS((timeshift%dt_shift-zdt_shift)/zdt_shift) > 1.e-10_wp) THEN
          WRITE(message_text,'(a,f10.3,a)') '*** WARNING: dt_shift adjusted to ', zdt_shift, &
            &                               ' s in order to be an integer multiple of dtime ***'
          CALL message('',message_text)
        ENDIF
        timeshift%dt_shift = zdt_shift
      END IF
      !
      ! transform timeshift to mtime-format
      !
      CALL getPTStringFromSeconds(timeshift%dt_shift, dt_shift_string)
      timeshift%mtime_shift => newTimedelta(dt_shift_string)
      WRITE(message_text,'(a,a)') 'IAU time shift: ', TRIM(dt_shift_string)
      !
      CALL getPTStringFromSeconds(ABS(timeshift%dt_shift), dt_shift_string)
      timeshift%mtime_absshift => newTimedelta(dt_shift_string)
      CALL message('',message_text)
#endif
    ENDIF





    ! --- --- STOP DATE:
    !
    !         This is the date when the current run stops time stepping.
    !
    !         The simulation stops at the next (i.e. earliest) of the
    !         following dates:
    !         a) experiment stop date
    !         b) start date + dt_restart
    !         c) start date + nsteps*dtime
    !
    NULLIFY(mtime_stop)
    NULLIFY(mtime_exp_stop)
    NULLIFY(mtime_restart_stop)
    NULLIFY(mtime_nsteps_stop)
    
    ! dtime is always available, maybe the default only
    CALL getPTStringFromMS(INT(dtime*1000,i8), td_string)
    mtime_dtime => newTimeDelta(td_string)
    IF (.NOT. ASSOCIATED(mtime_dtime))  CALL finish(routine, "Error in conversion of dtime to mtime!")

    mtime_start        => newDatetime(start_datetime_string, errno)
    IF (errno /= 0)  CALL finish(routine, "Error in conversion of start date: "//start_datetime_string)

    IF (exp_stop_datetime_string /= "") THEN
      mtime_exp_stop     => newDatetime(exp_stop_datetime_string, errno)
      IF (errno /= 0)  CALL finish(routine, "Error in conversion of exp stop date: "//exp_stop_datetime_string)
    END IF

    mtime_dt_restart   => newTimedelta(time_config%tc_dt_restart,errno)
    IF (errno /= 0)  CALL finish(routine, "Error in conversion of dt_restart!")

    mtime_restart_stop => newDatetime(mtime_start, errno)
    IF (errno /= 0)  CALL finish(routine, "Error in initialization of restart date")
    mtime_restart_stop =  mtime_restart_stop + mtime_dt_restart

    ! in case the restart is not given, mtime_restart_stop is deallocated

    IF (mtime_restart_stop == mtime_start) CALL deallocateDatetime(mtime_restart_stop)
    
    IF (nsteps >= 0) THEN   
      mtime_nsteps_stop  => newDatetime(mtime_start, errno)
      IF (errno /= 0)  CALL finish(routine, "Error in initialization of nsteps stop date")
      mtime_nsteps_stop = mtime_nsteps_stop + mtime_dtime * INT(nsteps,c_int32_t)
    ELSE
      IF (exp_stop_datetime_string /= "") THEN
        mtime_nsteps_stop  => newDatetime(mtime_exp_stop, errno)
        IF (errno /= 0)  CALL finish(routine, "Error in initialization of nsteps  stop date")
      END IF
    END IF

    IF (end_datetime_string /= "") THEN
      mtime_stop => newDatetime(end_datetime_string)
    ELSE
    
      !  in case stop date not given, we could simply set      !
      !   mtime_stop => newDatetime(MIN(MIN(mtime_exp_stop, mtime_restart_stop), mtime_nsteps_stop))
      !
      ! but we need to check cases where one or two of these dates have
      ! not been specified by the user...
      
      IF (ASSOCIATED(mtime_nsteps_stop)) THEN

        IF (ASSOCIATED(mtime_exp_stop)) THEN
          IF (ASSOCIATED(mtime_restart_stop)) THEN
            mtime_stop => newDatetime(MIN(MIN(mtime_exp_stop, mtime_restart_stop), mtime_nsteps_stop))
          ELSE
            mtime_stop => newDatetime(MIN(mtime_exp_stop, mtime_nsteps_stop))
          END IF
        ELSE
          IF (ASSOCIATED(mtime_restart_stop)) THEN
            mtime_stop => newDatetime(MIN(mtime_nsteps_stop, mtime_restart_stop))
          ELSE
            mtime_stop => newDatetime(mtime_nsteps_stop)
          END IF
        END IF

      ELSE

        IF (ASSOCIATED(mtime_exp_stop)) THEN
          IF (ASSOCIATED(mtime_restart_stop)) THEN
            mtime_stop => newDatetime(MIN(mtime_exp_stop, mtime_restart_stop))
          ELSE
            mtime_stop => newDatetime(mtime_exp_stop)
          END IF
        ELSE
          ! if neither "exp_stop_date" nor "nsteps" are given: throw
          ! an error regardless of the state of "restart_stop"
          CALL finish(routine, "Error in initialization of stop date")
        END IF

      END IF

    END IF

    CALL datetimeToString(mtime_stop, stop_datetime_string)

    ! if it has not been specified by the user, set experiment stop
    ! date to stop date:
    IF (.NOT. ASSOCIATED(mtime_exp_stop)) THEN
      mtime_exp_stop => newDatetime(stop_datetime_string)
    END IF

    CALL datetimeToString(mtime_exp_stop, exp_stop_datetime_string)
    
    ! consistency checks:
    !
    IF (mtime_stop < mtime_start) THEN
      CALL finish(routine, 'The end date and time must not be '// &
        &                  'before the current date and time')
    END IF

    CALL deallocateDatetime(mtime_start)
    CALL deallocateDatetime(mtime_stop)    
    CALL deallocateTimedelta(mtime_dtime)
    
    ! If a restart event occurs, check for unsupported combinations of
    ! namelist settings:
    IF (ASSOCIATED(mtime_nsteps_stop) .AND. ASSOCIATED(mtime_restart_stop)) THEN
      IF (.NOT. (mtime_nsteps_stop < mtime_restart_stop)) THEN
        ! processor splitting cannot be combined with synchronous restart:
        IF ((num_restart_procs == 0) .AND. ANY(patch_weight(1:) > 0._wp)) THEN
          CALL finish(routine, "Processor splitting cannot be combined with synchronous restart!")
        END IF
      END IF
    END IF

    IF (ASSOCIATED(mtime_exp_stop))  CALL deallocateDatetime(mtime_exp_stop)
    IF (ASSOCIATED(mtime_restart_stop))  CALL deallocateDatetime(mtime_restart_stop)
    IF (ASSOCIATED(mtime_nsteps_stop))  CALL deallocateDatetime(mtime_nsteps_stop)
    IF (INT(dt_restart,i8) > 0)  CALL deallocateTimedelta(mtime_dt_restart)

    ! --- --- NSTEPS
    !
    !         This specifies the number of steps in the time loop
    !         (INTEGER).
    !
    !         If the namelist parameter "nsteps" has not been
    !         specified in the namelist "run_nml", then this value is
    !         computed from the stop date and the start date. On the
    !         other hand, if the namelist parameter "nsteps" is
    !         specified in the namelist "run_nml", then this defines
    !         the end date (see above).
    !
    IF (nsteps < 0) THEN   
      ! User did not specified a value, we need to compute "nsteps" as
      ! (stop date - start date)/dtime:
      mtime_start => newDatetime(start_datetime_string)
      mtime_stop  => newDatetime(stop_datetime_string)
      mtime_td    => newTimedelta("PT0S")
      mtime_td    =  mtime_stop - mtime_start
      IF (mtime_td%year == 0 .AND. mtime_td%month == 0) THEN
        CALL getPTStringFromMS(INT(dtime*1000._wp,i8), td_string)
        mtime_dtime => newTimedelta(td_string)
        CALL divideTimeDeltaInSeconds(mtime_td, mtime_dtime, mtime_quotient)
        nsteps = INT(mtime_quotient%quotient)
        CALL deallocateTimedelta(mtime_dtime)
      ELSE
        CALL message('', 'Warning - cannot calculate nsteps unambiguous, set 0!')
        nsteps = 0
      END IF
      CALL deallocateTimedelta(mtime_td)
      CALL deallocateDatetime(mtime_start)
      CALL deallocateDatetime(mtime_stop)
    END IF

    ! --------------------------------------------------------------
    ! PART II: Convert ISO8601 strings into "mtime" and "t_datetime" 
    ! --------------------------------------------------------------

    ! --- Second, create date-time objects from the mtime library

    CALL set_tc_startdate    ( start_datetime_string     )
    CALL set_tc_stopdate     ( stop_datetime_string      )
    CALL set_tc_exp_startdate( exp_start_datetime_string )
    CALL set_tc_exp_stopdate ( exp_stop_datetime_string  )
    CALL set_tc_exp_refdate  ( exp_ref_datetime_string   )
    CALL set_tc_current_date ( cur_datetime_string       )

    IF (model_string == 'atm') THEN
#ifndef __NO_ICON_ATMO__
      ! add IAU time shift to current date
      IF (.NOT. isRestart()) THEN
        IF (timeshift%dt_shift < 0._wp) THEN
          time_config%tc_current_date = time_config%tc_current_date + timeshift%mtime_shift
        ENDIF
      ENDIF
#endif
    ENDIF


    ! --- Finally, store the same information in a "t_datetime" data
    !     structure
    CALL set_calendar        (dtime_calendar)
    CALL set_is_relative_time(is_relative_time)

    ! --------------------------------------------------------------
    ! PART III: Print all date and time components
    ! --------------------------------------------------------------

    IF (ldebug) THEN

      CALL message('DEBUG','')
      CALL message('DEBUG','Calendar: '//TRIM(master_nml_calendar))      
      CALL message('DEBUG','Calendar: '//TRIM(calendar_index2string(time_nml_icalendar))//' (deprecated interface)')
      CALL message('DEBUG','')

      WRITE(message_text,'(a,a)') 'Model time step         : ', TRIM(mtime_modelTimeStep)
      CALL message('DEBUG',message_text)
      WRITE(message_text,'(a,g0,a)') 'Model time step         : ', dtime, ' (deprecated interface)' 
      CALL message('DEBUG',message_text)
      WRITE(message_text,'(a,a)') 'Checkpoint time interval: ', TRIM(checkpointTimeIntval)
      CALL message('DEBUG',message_text)
      WRITE(message_text,'(a,g0,a)') 'Checkpoint time interval: ', dt_checkpoint, ' (deprecated interface)' 
      CALL message('DEBUG',message_text)
      WRITE(message_text,'(a,a)') 'Restart time interval   : ', TRIM(restartTimeIntval)
      CALL message('DEBUG',message_text)
      WRITE(message_text,'(a,g0,a)') 'Restart time interval   : ', dt_restart, ' (deprecated interface)' 
      CALL message('DEBUG',message_text)
      CALL message('DEBUG','')
      WRITE(message_text,'(a,a)') 'Experiment reference date: ', TRIM(experimentReferenceDate)
      CALL message('DEBUG',message_text)
      WRITE(message_text,'(a,a)') 'Experiment start date    : ', TRIM(experimentStartDate)
      CALL message('DEBUG',message_text)
      WRITE(message_text,'(a,a)') 'Experiment stop date     : ', TRIM(experimentStopDate)
      CALL message('DEBUG',message_text)
      CALL message('DEBUG','')
      
      CALL message('',message_text)
      IF (isRestart()) THEN
        WRITE(message_text,'(a,a,a)') 'Start date      : ', TRIM(start_datetime_string), ' (deduced, restart run)'
      ELSE
        WRITE(message_text,'(a,a,a)') 'Start date      : ', TRIM(start_datetime_string), ' (deduced)'
      END IF
      CALL message('DEBUG',message_text)
      WRITE(message_text,'(a,a)')     'Stop date       : ', TRIM(end_datetime_string)
      CALL message('DEBUG',message_text)
      WRITE(message_text,'(a,i0,a)')  'Stop date, steps: ', nsteps, ' (deprecated interface)'  
      CALL message('DEBUG',message_text)
      CALL message('DEBUG','')
    END IF
    
    CALL message('','')
    CALL calendarToString(dstring)
    CALL message('','Calendar: '//TRIM(dstring))
    CALL message('','')

    WRITE(message_text,'(a,a)')     'Experiment reference date: ', TRIM(exp_ref_datetime_string)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)')     'Experiment start date    : ', TRIM(exp_start_datetime_string)
    CALL message('',message_text)
    WRITE(message_text,'(a,a)')     'Experiment stop date     : ', TRIM(exp_stop_datetime_string)
    CALL message('',message_text)
    CALL message('','')

    CALL message('',message_text)
    IF (isRestart()) THEN
      WRITE(message_text,'(a,a,a)') 'Start date    : ', TRIM(start_datetime_string), ' (restart run)'
    ELSE
      WRITE(message_text,'(a,a)')   'Start date    : ', TRIM(start_datetime_string)
    END IF
    CALL message('',message_text)
    WRITE(message_text,'(a,a)')     'Stop date     : ', TRIM(stop_datetime_string)
    CALL message('',message_text)
    CALL message('','')
    
  END SUBROUTINE compute_date_settings
 
END MODULE mo_time_management
