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

MODULE mo_time_config

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish
  USE mtime,                    ONLY: datetime, timedelta, newDatetime, newTimedelta, &
    &                                 deallocateDatetime, MAX_CALENDAR_STR_LEN,       &
    &                                 OPERATOR(*)
  USE mo_impl_constants,        ONLY: proleptic_gregorian, julian_gregorian, cly360,  &
    &                                 SUCCESS
  USE mo_util_string,           ONLY: tolower
  USE mo_model_domain,          ONLY: t_patch
 
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calendar_index2string
  PUBLIC :: ini_datetime_string, end_datetime_string, icalendar
  PUBLIC :: restart_ini_datetime_string, restart_end_datetime_string, restart_calendar
  PUBLIC :: dt_restart, is_relative_time
  PUBLIC :: t_time_config, time_config
  PUBLIC :: set_tc_exp_refdate
  PUBLIC :: set_tc_exp_startdate, set_tc_exp_stopdate
  PUBLIC :: set_tc_startdate, set_tc_stopdate
  PUBLIC :: set_tc_dt_checkpoint
  PUBLIC :: set_tc_dt_restart
  PUBLIC :: set_tc_current_date
  PUBLIC :: set_is_relative_time
  PUBLIC :: set_calendar
  PUBLIC :: set_tc_dt_model
  PUBLIC :: set_tc_write_restart
  PUBLIC :: get_dynamics_timestep
  
  !> namelist parameters (as raw character strings):
  !
  ! these are the namelist settings originating from the restart file:
  CHARACTER(len=32)                   :: restart_ini_datetime_string
  CHARACTER(len=32)                   :: restart_end_datetime_string
  INTEGER                             :: restart_calendar
  !
  !  these are namelist setting which may originate from the restart
  !  file, but with user modifications in the current run:
  INTEGER                             :: icalendar
  CHARACTER(len=32)                   :: ini_datetime_string
  CHARACTER(len=32)                   :: end_datetime_string
  REAL(wp)                            :: dt_restart          !< Length of restart cycle in seconds
  LOGICAL                             :: is_relative_time

  !>
  !! Derived type containing information for time control. 
  !!
  TYPE t_time_config

    ! from namelist 

    REAL(wp)         :: dt_restart         !< Length of restart cycle in seconds
    INTEGER          :: calendar           !< calendar type

    ! not directly from namelist  

    !> LOGICAL is_relative_time: .TRUE., if time loop shall start with
    !> step 0 regardless whether we are in a standard run or in a
    !> restarted run (which means re-initialized run):
    LOGICAL          ::  is_relative_time

    ! whole experiment and single run time information
    ! -----------------------------------------------------------------
    !
    ! experiment
    
    TYPE(datetime),  POINTER :: tc_exp_refdate   => NULL()

    TYPE(datetime),  POINTER :: tc_exp_startdate => NULL()
    TYPE(datetime),  POINTER :: tc_exp_stopdate  => NULL()

    ! single run 

    TYPE(datetime),  POINTER :: tc_startdate     => NULL()
    TYPE(datetime),  POINTER :: tc_stopdate      => NULL()

    ! current model date

    TYPE(datetime),  POINTER :: tc_current_date  => NULL()

    ! check point and restart time interval (needs to be equal for all
    ! component models)

    TYPE(timedelta), POINTER :: tc_dt_checkpoint => NULL()
    TYPE(timedelta), POINTER :: tc_dt_restart    => NULL()

    ! in case no restart time interval is given, no restart
    ! should be written - assumption is that this is the default

    LOGICAL :: tc_write_restart = .TRUE.

    ! well, the model's timestep
    
    TYPE(timedelta), POINTER     :: tc_dt_model => NULL() ! dynamics time step  on the global grid in mtime format

  CONTAINS
    ! create a copy of t_time_config object
    PROCEDURE  :: copy      => time_config_copy
    ! destruct t_time_config object
    PROCEDURE  :: destruct  => time_config_destruct
  END TYPE t_time_config
  !>
  !!
  !! The actual variable
  !!
  TYPE(t_time_config), TARGET, SAVE :: time_config

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_time_config"

CONTAINS

  !> Create manual deep copy of a t_time_config object
  !
  SUBROUTINE time_config_copy(me, tc_new)
    CLASS(t_time_config)               :: me
    TYPE(t_time_config), INTENT(INOUT) :: tc_new

    INTEGER :: ist
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::time_config_copy'

    ALLOCATE(tc_new%tc_exp_refdate, tc_new%tc_exp_startdate, tc_new%tc_exp_stopdate, &
      &      tc_new%tc_startdate, tc_new%tc_stopdate, tc_new%tc_current_date,        &
      &      tc_new%tc_dt_checkpoint, tc_new%tc_dt_restart, tc_new%tc_dt_model,      &
      &      stat=ist)
    IF(ist /= SUCCESS) CALL finish(routine, "memory allocation failure")

    tc_new%dt_restart       = me%dt_restart
    tc_new%calendar         = me%calendar
    tc_new%is_relative_time = me%is_relative_time
    tc_new%tc_exp_refdate   = me%tc_exp_refdate
    tc_new%tc_exp_startdate = me%tc_exp_startdate
    tc_new%tc_exp_stopdate  = me%tc_exp_stopdate
    tc_new%tc_startdate     = me%tc_startdate
    tc_new%tc_stopdate      = me%tc_stopdate
    tc_new%tc_current_date  = me%tc_current_date
    tc_new%tc_dt_checkpoint = me%tc_dt_checkpoint
    tc_new%tc_dt_restart    = me%tc_dt_restart
    tc_new%tc_write_restart = me%tc_write_restart
    tc_new%tc_dt_model      = me%tc_dt_model

  END SUBROUTINE time_config_copy


  !> Destructs an object of type t_time_config
  !
  SUBROUTINE time_config_destruct(me)
    CLASS(t_time_config) :: me

    INTEGER :: ist
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::time_config_copy'

    IF (ASSOCIATED(me%tc_exp_refdate)) DEALLOCATE(me%tc_exp_refdate, stat=ist)
    IF(ist /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    NULLIFY(me%tc_exp_refdate)
    !
    IF (ASSOCIATED(me%tc_exp_startdate)) DEALLOCATE(me%tc_exp_startdate, stat=ist)
    IF(ist /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    NULLIFY(me%tc_exp_startdate)
    !
    IF (ASSOCIATED(me%tc_exp_stopdate)) DEALLOCATE(me%tc_exp_stopdate, stat=ist)
    IF(ist /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    NULLIFY(me%tc_exp_stopdate)
    !
    IF (ASSOCIATED(me%tc_startdate)) DEALLOCATE(me%tc_startdate, stat=ist)
    IF(ist /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    NULLIFY(me%tc_startdate)
    !
    IF (ASSOCIATED(me%tc_stopdate)) DEALLOCATE(me%tc_stopdate, stat=ist)
    IF(ist /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    NULLIFY(me%tc_stopdate)
    !
    IF (ASSOCIATED(me%tc_current_date)) DEALLOCATE(me%tc_current_date, stat=ist)
    IF(ist /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    NULLIFY(me%tc_current_date)
    !
    IF (ASSOCIATED(me%tc_dt_checkpoint)) DEALLOCATE(me%tc_dt_checkpoint, stat=ist)
    IF(ist /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    NULLIFY(me%tc_dt_checkpoint)
    !
    IF (ASSOCIATED(me%tc_dt_restart)) DEALLOCATE(me%tc_dt_restart, stat=ist)
    IF(ist /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    NULLIFY(me%tc_dt_restart)
    !
    IF (ASSOCIATED(me%tc_dt_model)) DEALLOCATE(me%tc_dt_model, stat=ist)
    IF(ist /= SUCCESS) CALL finish(routine, "memory deallocation failure")
    NULLIFY(me%tc_dt_model)

  END SUBROUTINE time_config_destruct


  !> Convert the calendar setting (which is an integer value for this
  !  namelist) into a string. The naming scheme is then compatible
  !  with concurrent namelist settings of the calendar (mtime).
  !
  FUNCTION calendar_index2string(icalendar) RESULT(ret)
    CHARACTER(LEN=MAX_CALENDAR_STR_LEN) :: ret
    INTEGER, INTENT(IN) :: icalendar
    ! local variables
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::calendar_index2string'

    ret = ""
    SELECT CASE(icalendar)
    CASE(julian_gregorian)
      ret = 'julian gregorian'
    CASE(proleptic_gregorian)
      ret = 'proleptic gregorian'
    CASE(cly360)
      ret = '360 day year'
    END SELECT
  END FUNCTION calendar_index2string

  !> Convert the calendar setting (which is an integer value for this
  !  namelist) into a string. The naming scheme is then compatible
  !  with concurrent namelist settings of the calendar (mtime).
  !
  FUNCTION calendar_string2index(cal_str) RESULT(ret)
    INTEGER :: ret
    CHARACTER(LEN=*), INTENT(IN) :: cal_str
    ! local variables
    CHARACTER(len=*), PARAMETER ::  routine = modname//'::calendar_string2index'

    ret = -1
    IF (TRIM(tolower(cal_str)) == 'julian gregorian') THEN
      ret = julian_gregorian
    ELSE IF (TRIM(tolower(cal_str)) == 'proleptic gregorian') THEN
      ret = proleptic_gregorian
    ELSE IF (TRIM(tolower(cal_str)) == '360 day year') THEN
      ret = cly360
    END IF
  END FUNCTION calendar_string2index



  SUBROUTINE set_tc_exp_refdate(experimentReferenceDate)   
    CHARACTER(len=*), INTENT(in) :: experimentReferenceDate   
    time_config%tc_exp_refdate => newDatetime(experimentReferenceDate)
  END SUBROUTINE set_tc_exp_refdate

  SUBROUTINE set_tc_exp_startdate(experimentStartDate)   
    CHARACTER(len=*), INTENT(in) :: experimentStartDate   
    time_config%tc_exp_startdate => newDatetime(experimentStartDate)   
  END SUBROUTINE set_tc_exp_startdate

  SUBROUTINE set_tc_exp_stopdate(experimentStopDate)
    CHARACTER(len=*), INTENT(in) :: experimentStopDate
    time_config%tc_exp_stopdate => newDatetime(experimentStopDate)
  END SUBROUTINE set_tc_exp_stopdate

  SUBROUTINE set_tc_startdate(startdate)
    CHARACTER(len=*), INTENT(in) :: startdate
    time_config%tc_startdate => newDatetime(startdate)
  END SUBROUTINE set_tc_startdate

  SUBROUTINE set_tc_stopdate(stopdate)
    CHARACTER(len=*), INTENT(in) :: stopdate
    time_config%tc_stopdate => newDatetime(stopdate)
  END SUBROUTINE set_tc_stopdate

  SUBROUTINE set_tc_dt_checkpoint(checkpointTimeIntval)
    CHARACTER(len=*), INTENT(in) :: checkpointTimeIntval
    time_config%tc_dt_checkpoint => newTimedelta(checkpointTimeIntval)
  END SUBROUTINE set_tc_dt_checkpoint
  
  SUBROUTINE set_tc_dt_restart(restartTimeIntval)
    CHARACTER(len=*), INTENT(in) :: restartTimeIntval   
    time_config%tc_dt_restart => newTimedelta(restartTimeIntval)
  END SUBROUTINE set_tc_dt_restart

  SUBROUTINE set_tc_current_date(current_date)
    CHARACTER(len=*), INTENT(in) :: current_date
    IF (ASSOCIATED(time_config%tc_current_date)) THEN
      CALL deallocateDatetime(time_config%tc_current_date)
    END IF
    time_config%tc_current_date => newDatetime(current_date)
  END SUBROUTINE set_tc_current_date

  SUBROUTINE set_calendar(icalendar)
    INTEGER, INTENT(IN) :: icalendar
    time_config%calendar = icalendar
  END SUBROUTINE set_calendar

  SUBROUTINE set_is_relative_time(lvalue)
    LOGICAL, INTENT(IN) :: lvalue
    time_config%is_relative_time = lvalue
  END SUBROUTINE set_is_relative_time

  SUBROUTINE set_tc_dt_model(modelTimeStep)
    CHARACTER(len=*), INTENT(in) :: modelTimeStep
    time_config%tc_dt_model => newTimedelta(modelTimeStep)
  END SUBROUTINE set_tc_dt_model

  SUBROUTINE set_tc_write_restart(writeRestart)
    LOGICAL, INTENT(in) :: writeRestart
    time_config%tc_write_restart = writeRestart
  END SUBROUTINE set_tc_write_restart
  
  !> Compute the effective dynamics timestep for the given patch.
  !! The timestep is halved for every level of nesting compared to the "global" timestep.
  TYPE(timedelta) FUNCTION get_dynamics_timestep(patch) RESULT(dt_dyn)
    TYPE(t_patch), INTENT(IN) :: patch !< Patch whose timestep will be returned.

    INTEGER :: k

    dt_dyn = time_config%tc_dt_model
    DO k = 1, patch%nest_level
      dt_dyn = dt_dyn * 0.5_wp
    END DO

  END FUNCTION get_dynamics_timestep

END MODULE mo_time_config

