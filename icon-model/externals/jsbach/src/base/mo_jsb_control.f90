!> Contains methods for JSBACH master control.
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_control
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_io_units,  ONLY: filename_max
  USE mo_exception, ONLY: finish
  USE mo_util,      ONLY: int2string
  USE mo_timer,     ONLY: new_timer

  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
  USE mo_jsb_time_iface,     ONLY: l_timer_host !< timers on in host model?
  USE mo_jsb_io_iface,       ONLY: ldebugio

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_jsb_master_control
  PUBLIC :: jsbach_runs_standalone, jsbach_is_restarted, force_from_observations
  PUBLIC :: debug_on, timer_on, l_timer_host, get_debug_memory_level
  PUBLIC :: model_base_dir
  ! t_jsb_model_nml needs to be public for the NEC compiler
  PUBLIC :: t_jsb_model_nml, jsb_models_nml, get_no_of_models
  PUBLIC :: timer_jsbach, timer_aggregate, timer_integrate, timer_process_task, timer_process_msg, timer_forcing
  PUBLIC :: max_no_of_models

  ! Needed for inlining timer_on and debug_on on NEC.
  PUBLIC :: debug_level, timer_level
  PROTECTED :: debug_level, timer_level

  LOGICAL,           SAVE :: is_standalone
  LOGICAL,           SAVE :: restart_jsbach
  CHARACTER(len=99), SAVE :: model_base_dir
  INTEGER,           SAVE :: debug_level
  INTEGER,           SAVE :: debug_io
  INTEGER,           SAVE :: debug_memory_level
  INTEGER,           SAVE :: timer_level
  LOGICAL,           SAVE :: l_force_from_obs

  TYPE t_jsb_model_nml
    INTEGER            :: model_id = 0
    CHARACTER(len=30)  :: model_name
    CHARACTER(len=SHORT_NAME_LEN) :: model_shortname
    CHARACTER(len=132) :: model_description
    CHARACTER(len=filename_max) :: model_namelist_filename
  END TYPE t_jsb_model_nml

  INTEGER, PARAMETER    :: max_no_of_models = 10
  INTEGER :: no_of_models = 0
  TYPE(t_jsb_model_nml), SAVE :: jsb_models_nml(max_no_of_models)

  INTEGER, SAVE, DIMENSION(max_no_of_models) :: &
    & timer_jsbach, timer_aggregate, timer_integrate, timer_process_task, timer_process_msg, timer_forcing

  CHARACTER(len=filename_max) :: master_namelist_filename = 'jsbach_master.namelist'

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_control'

CONTAINS

  SUBROUTINE init_jsb_master_control(namelist_filename)

    USE mo_jsb_io,      ONLY: init_jsb_io
    USE mo_jsb_version, ONLY: jsbach_init_version

    CHARACTER(len=*), OPTIONAL, INTENT(in) :: namelist_filename

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_jsb_master_control'

    IF (PRESENT(namelist_filename)) master_namelist_filename = TRIM(namelist_filename)

    CALL read_jsb_control_namelist(TRIM(master_namelist_filename)) ! This sets *no_of_models*

    CALL jsbach_init_version(jsbach_runs_standalone())

    CALL init_jsb_io(TRIM(master_namelist_filename))

    IF (timer_on() .OR. l_timer_host) THEN
      IF (no_of_models > 1) THEN
        DO i=1,no_of_models
          timer_jsbach(i)  = new_timer('jsbach_'//int2string(i))
          timer_forcing(i) = new_timer('jsb:forcing_'//int2string(i))
        END DO
      ELSE
        timer_jsbach(1)  = new_timer('jsbach')
        timer_forcing(1) = new_timer('jsb:forcing')
      END IF
    END IF
    IF (timer_on()) THEN
      IF (no_of_models > 1) THEN
        DO i=1,no_of_models
          timer_aggregate(i)    = new_timer('jsb:aggregate_'//int2string(i))
          timer_integrate(i)    = new_timer('jsb:integrate_'//int2string(i))
          timer_process_task(i) = new_timer('jsb:process_task_'//int2string(i))
          timer_process_msg(i)  = new_timer('jsb:process_msg_'//int2string(i))
        END DO
      ELSE
        timer_aggregate(1)    = new_timer('jsb:aggregate')
        timer_integrate(1)    = new_timer('jsb:integrate')
        timer_process_task(1) = new_timer('jsb:process_task')
        timer_process_msg(1)  = new_timer('jsb:process_msg')
      END IF
    END IF

  END SUBROUTINE init_jsb_master_control

  !>
  !! Reads JSBACH control parameters.
  !!
  !! Reads namelist 'jsb_control_nml' and namelists 'jsb_model_nml' from namelist file.
  !!
  SUBROUTINE read_jsb_control_namelist(namelist_filename)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CHARACTER(len=*), INTENT(in) :: namelist_filename

    INTEGER            :: model_id
    CHARACTER(len=30)  :: model_name
    CHARACTER(len=SHORT_NAME_LEN) :: model_shortname
    CHARACTER(len=132) :: model_description
    CHARACTER(len=filename_max) :: model_namelist_filename

    INTEGER :: nml_handler, nml_unit, istat, i
    LOGICAL :: rwnd

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_jsb_control_namelist'

    NAMELIST /jsb_control_nml/ &
      is_standalone,           &
      restart_jsbach,          &
      model_base_dir,          &
      timer_level,             &
      debug_level,             &
      debug_io,                &
      debug_memory_level,      &
      l_force_from_obs

    NAMELIST /jsb_model_nml/   &
      model_id,                &
      model_name,              &
      model_shortname,         &
      model_description,       &
      model_namelist_filename

    ! Namelist defaults
    is_standalone      = .TRUE.
    restart_jsbach     = .FALSE.
    model_base_dir     = '.'
    timer_level        = 0
    debug_level        = 0
    debug_io           = 0
    debug_memory_level = 0
    l_force_from_obs   = .TRUE.

    nml_handler = open_nml(TRIM(namelist_filename))

    nml_unit = position_nml('jsb_control_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_control_nml)

    ! Get number of models first
    rwnd = .TRUE.
    no_of_models = 0
    DO
      nml_unit = position_nml('jsb_model_nml', nml_handler, lrewind=rwnd, STATUS=istat)
      IF (istat /= POSITIONED) EXIT
      rwnd=.FALSE.
      no_of_models = no_of_models + 1
      READ(nml_unit, jsb_model_nml)
    END DO

    IF (no_of_models == 0 ) THEN
      CALL finish(routine, 'Namelist jsb_model_nml missing in '//TRIM(namelist_filename)// &
        & '. One jsb_model_nml section is needed for each model incidence.')
    END IF
    IF (no_of_models > max_no_of_models) THEN
      CALL finish(routine, 'no_of_models > max_no_of_models')
    ENDIF

    IF (is_standalone .AND. .NOT. l_force_from_obs) THEN
      CALL finish(routine, 'l_force_from_obs needs to be true - also with forcing data generated from '// &
        & 'model output. Forcing from atmosphere interface data currently not implemented.')
    END IF

    ldebugio = debug_io > 0

    ! Now read namelists again and set parameters
    ! Note: no_of_models needs to be known here (therfore extra loop above) in order to test for a valid model_id below
    rwnd = .TRUE.
    DO i=1,no_of_models
      nml_unit = position_nml('jsb_model_nml', nml_handler, lrewind=rwnd, STATUS=istat)
      IF (istat /= POSITIONED) THEN
        CALL finish(TRIM(routine), 'jsb_model_nml not found in namelist')
      END IF
      rwnd=.FALSE.

      model_id = i
      model_name = ''
      model_shortname = ''
      model_description = ''
      model_namelist_filename = ''

      READ(nml_unit, jsb_model_nml)

      IF (jsb_models_nml(model_id)%model_id > 0) THEN
        CALL finish(TRIM(routine), 'Model id already taken.')
      END IF

      IF (model_id > no_of_models) THEN
        CALL finish(TRIM(routine), 'Model id larger than number of models.')
      END IF

      jsb_models_nml(i)%model_id                = model_id
      jsb_models_nml(i)%model_name              = model_name
      jsb_models_nml(i)%model_shortname         = model_shortname
      jsb_models_nml(i)%model_description       = model_description
      jsb_models_nml(i)%model_namelist_filename = model_namelist_filename
    END DO

    CALL close_nml(nml_handler)

  END SUBROUTINE read_jsb_control_namelist

  LOGICAL FUNCTION jsbach_runs_standalone()

    jsbach_runs_standalone = is_standalone

  END FUNCTION jsbach_runs_standalone

  LOGICAL FUNCTION force_from_observations()

    force_from_observations = l_force_from_obs

  END FUNCTION force_from_observations

  LOGICAL FUNCTION jsbach_is_restarted()

    ! TBD: get this from parent model if not standalone
    jsbach_is_restarted = restart_jsbach

  END FUNCTION jsbach_is_restarted

  INTEGER FUNCTION get_no_of_models()

    get_no_of_models = no_of_models

  END FUNCTION get_no_of_models

  LOGICAL FUNCTION debug_on(level)

    CHARACTER(len=*), OPTIONAL, INTENT(in) :: level

    CHARACTER(len=*), PARAMETER :: routine = modname//':debug_on'

    IF (PRESENT(level)) THEN
      SELECT CASE (TRIM(level))
      CASE ('basic')
        debug_on = debug_level >= 1
      CASE ('detail')
        debug_on = debug_level >= 2
      CASE ('hsm')
        debug_on = debug_level >= 3
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Unknown debug level')
      END SELECT
    ELSE
      debug_on = debug_level > 0
    END IF

  END FUNCTION debug_on

  INTEGER FUNCTION get_debug_memory_level()

    get_debug_memory_level = debug_memory_level

  END FUNCTION get_debug_memory_level

  LOGICAL FUNCTION timer_on(level)

    CHARACTER(len=*), OPTIONAL, INTENT(in) :: level

    CHARACTER(len=*), PARAMETER :: routine = modname//':timer_on'

    IF (PRESENT(level)) THEN
      SELECT CASE (TRIM(level))
      CASE ('basic')
        timer_on = timer_level >= 1
      CASE ('detail')
        timer_on = timer_level >= 2
      CASE ('all')
        timer_on = timer_level >= 3
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Unknown timer level')
      END SELECT
    ELSE
      timer_on = timer_level > 0
    END IF

  END FUNCTION timer_on

#endif
END MODULE mo_jsb_control
