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

MODULE mo_wave
  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: message, finish
  USE mo_model_domain,          ONLY: p_patch
  USE mo_master_config,         ONLY: isRestart
  USE mo_master_control,        ONLY: get_my_process_name
  USE mo_grid_config,           ONLY: n_dom, start_time, end_time
  USE mo_wave_state,            ONLY: construct_wave_state, destruct_wave_state
  USE mo_wave_forcing_state,    ONLY: construct_wave_forcing_state, destruct_wave_forcing_state
  USE mo_time_config,           ONLY: time_config
  USE mo_util_mtime,            ONLY: getElapsedSimTimeInSeconds
  USE mo_output_event_types,    ONLY: t_sim_step_info
  USE mo_name_list_output_init, ONLY: init_name_list_output, parse_variable_groups, &
    &                                 output_file, create_vertical_axes
  USE mo_var_list_register_utils, ONLY: vlr_print_groups
  USE mo_run_config,            ONLY: output_mode, msg_level
  USE mo_io_config,             ONLY: configure_io
  USE mo_mpi,                   ONLY: my_process_is_stdio
  USE mo_wave_stepping,         ONLY: perform_wave_stepping
  USE mo_wave_events,           ONLY: create_wave_events
  USE mo_opt_diagnostics,       ONLY: construct_opt_diag, destruct_opt_diag
  USE mo_pp_scheduler,          ONLY: pp_scheduler_init, pp_scheduler_finalize
  ! restart
  USE mo_restart,               ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
  USE mo_load_restart,          ONLY: read_restart_files
  USE mo_restart_nml_and_att,   ONLY: getAttributesForRestarting
  USE mo_key_value_store,       ONLY: t_key_value_store
  USE mo_timer,                 ONLY: timers_level, timer_start, timer_stop, timer_read_restart

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wave


CONTAINS

  SUBROUTINE wave

    CHARACTER(*), PARAMETER :: routine = "mo_wave:wave"
    CLASS(t_RestartDescriptor), POINTER  :: restartDescriptor

#ifdef _OPENACC
    CALL finish(routine, "wave_process not ported to GPU yet")
#endif
    CALL construct_wave()

    restartDescriptor => createRestartDescriptor(TRIM(get_my_process_name()))

    CALL perform_wave_stepping(time_config       = time_config, &
      &                        restartDescriptor = restartDescriptor)

    CALL deleteRestartDescriptor(restartDescriptor)

    CALL destruct_wave()

    CALL message(TRIM(routine),'finished')

  END SUBROUTINE wave


  SUBROUTINE construct_wave()

    CHARACTER(*), PARAMETER :: routine = "mo_wave:construct_wave"

    TYPE(t_sim_step_info) :: sim_step_info
    TYPE(t_key_value_store), POINTER :: restartAttributes
    INTEGER :: jg
    REAL(wp):: sim_time


    ! calculate elapsed simulation time in seconds
    sim_time = getElapsedSimTimeInSeconds(time_config%tc_current_date)

    DO jg=1, n_dom
      p_patch(jg)%ldom_active &
           =        (jg <= 1 .OR. start_time(jg) <= sim_time) &
           &  .AND. end_time(jg) > sim_time
    END DO

    CALL construct_wave_state(p_patch(1:),n_timelevels=2)

    CALL construct_wave_forcing_state(p_patch(1:))

    ! Add optional diagnostic variable list (which might remain empty)
    ! primarily for lat-lon output
    CALL construct_opt_diag(p_patch(1:), .FALSE.)


    ! create wave events (restart, checkpointing)
    CALL create_wave_events(time_config)

    !------------------------------------------------------------------
    ! Prepare initial conditions for time integration.
    !------------------------------------------------------------------
    !
    IF (isRestart()) THEN
      !
      ! This is a resumed integration. Read model state from restart file(s).
      !
      IF (timers_level > 4) CALL timer_start(timer_read_restart)
      !
      DO jg = 1,n_dom
        IF (p_patch(jg)%ldom_active) THEN
          CALL read_restart_files( p_patch(jg), n_dom)
        END IF
      END DO
      !
      CALL message(routine,'normal exit from read_restart_files')
      !
      IF (timers_level > 4) CALL timer_stop(timer_read_restart)
    ELSE
      !
      ! Initialize with real atmospheric data
      !
      ! TO BE IMPLEMENTED
    ENDIF

    !------------------------------------------------------------------
    ! Prepare output file
    !------------------------------------------------------------------
    CALL configure_io()   ! set n_chkpt and n_diag, which control
                          ! writing of restart files and tot_int diagnostics.

    IF (output_mode%l_nml) THEN
       CALL parse_variable_groups()
    END IF

    ! setup of post-processing job queue, e.g. setup of optional
    ! diagnostic quantities like pz-level interpolation ot lat-lon output
    CALL pp_scheduler_init(l_init_prm_diag=.FALSE.)

    ! If async IO is in effect, init_name_list_output is a collective call
    ! with the IO procs and effectively starts async IO
    IF (output_mode%l_nml) THEN
       ! compute sim_start, sim_end
       sim_step_info%sim_start = time_config%tc_exp_startdate
       sim_step_info%sim_end = time_config%tc_exp_stopdate
       sim_step_info%run_start = time_config%tc_startdate
       sim_step_info%restart_time = time_config%tc_stopdate

       sim_step_info%dtime  = time_config%get_model_timestep_sec(p_patch(1)%nest_level)
       sim_step_info%jstep0 = 0

       CALL getAttributesForRestarting(restartAttributes)
       ! get start counter for time loop from restart file:
       IF (restartAttributes%is_init) THEN
         CALL restartAttributes%get("jstep", sim_step_info%jstep0)
       ENDIF

       CALL init_name_list_output(sim_step_info)

       CALL create_vertical_axes(output_file)
    END IF

    !-------------------------------------------------------!
    !  (Optional) detailed print-out of some variable info  !
    !-------------------------------------------------------!
    ! variable group information
    IF (my_process_is_stdio() .AND. (msg_level >= 15)) THEN
       CALL vlr_print_groups(idom=1,                            &
            &                opt_latex_fmt           = .FALSE., &
            &                opt_reduce_trailing_num = .TRUE.,  &
            &                opt_skip_trivial        = .TRUE.)
    END IF


    CALL message(TRIM(routine),'finished')

  END SUBROUTINE construct_wave


  SUBROUTINE destruct_wave()

    CHARACTER(*), PARAMETER :: routine = "destruct_wave"

    ! Destruction of post-processing job queue
    CALL pp_scheduler_finalize()

    ! Delete optional diagnostics
    CALL destruct_opt_diag()

    CALL destruct_wave_state()

    CALL destruct_wave_forcing_state()

    CALL message(TRIM(routine),'finished')

  END SUBROUTINE destruct_wave


END MODULE mo_wave
