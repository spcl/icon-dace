!> ICON driver for the JSBACH standalone model
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
!>#### Main routine for running ICON-Land standalone
!>
!> `jsbach_model()` is called from the ICON program `src/drivers/icon:icon`
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsbach_model

#ifdef __ICON__
#ifndef __NO_JSBACH__

  USE mo_exception, ONLY: message, message_text, finish
  USE mo_kind,      ONLY: wp

  USE mtime,        ONLY: datetime, newDatetime, datetimeToString, timedelta, newTimedelta, getPTStringFromSeconds, &
    &                     mtime_strerror, timedeltaToString, deallocateDatetime, deallocateTimedelta, &
    &                     OPERATOR(+), OPERATOR(==), OPERATOR(<), OPERATOR(>), OPERATOR(>=), &
    &                     event, eventGroup, newEvent, addEventToEventGroup, &
    &                     NO_ERROR, MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN, MAX_MTIME_ERROR_STR_LEN
  USE mo_timer,     ONLY: init_timer, timer_start, timer_stop, timers_level,    &
    &                     timer_model_init, print_timer,                        &
    &                     timer_domain_decomp, timer_compute_coeffs,            &
    &                     timer_read_restart, timer_total
  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_class,       ONLY: get_model

  USE mo_mpi,       ONLY: i_am_accel_node

  IMPLICIT NONE

  PRIVATE

#ifdef HAVE_CDI_PIO
    INCLUDE 'cdipio.inc'
#endif

  PUBLIC :: jsbach_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsbach_model'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> main routine for running ICON-Land standalone using ICON infrastructure
  !>
  !>   'jsbach_model()' is called from the ICON program `src/drivers/icon:icon` (ICON-MPIM repo)
  !>
  SUBROUTINE jsbach_model(jsb_namelist_filename, shr_namelist_filename)

    USE mo_master_control,          ONLY: jsbach_process, get_my_process_name
    USE mo_master_config,           ONLY: isRestart

    USE mo_load_restart,            ONLY: read_restart_header, read_restart_files
    USE mo_restart_nml_and_att,     ONLY: getAttributesForRestarting
    USE mo_key_value_store,         ONLY: t_key_value_store

    USE mo_read_namelists_iface,    ONLY: read_infrastructure_namelists_for_jsbach
    USE mo_namelist,                ONLY: open_nml_output, close_nml_output
    USE mo_nml_annotate,            ONLY: log_nml_settings
    USE mo_name_list_output_init,   ONLY: read_name_list_output_namelists, init_name_list_output, parse_variable_groups, &
      &                                   output_file, create_vertical_axes
    USE mo_level_selection,         ONLY: create_mipz_level_selections

    USE mo_time_management,          ONLY: compute_timestep_settings,                        &
      &                                    compute_restart_settings,                         &
      &                                    compute_date_settings
    USE mo_time_config,              ONLY: time_config, dt_restart
    USE mo_run_config,               ONLY: configure_run, dtime, nsteps, ltimer, output_mode, &
      &                                    nshift, num_lev, msg_level, grid_generatingCenter, grid_generatingSubcenter
    USE mo_event_manager,            ONLY: initEventManager, addEventGroup, getEventGroup, printEventGroup

    USE mo_io_config,               ONLY: restartWritingParameters, configure_io

    USE mo_memory_log,              ONLY: memory_log_add

    USE mo_mpi,                     ONLY: stop_mpi, my_process_is_io, my_process_is_work, my_process_is_stdio, &
      &                                   set_mpi_work_communicators, process_mpi_io_size
    USE mo_parallel_config,         ONLY: p_test_run, num_test_pe, l_test_openmp,                              &
      &                                   num_io_procs, pio_type,                     &
      &                                   num_prefetch_proc, nproma, proc0_offloading, check_parallel_configuration

    USE mo_zaxis_type,              ONLY: zaxisTypeList, t_zaxisTypeList
    USE mo_jsb_base,                ONLY: jsbach_setup_models, jsbach_setup_tiles
    USE mo_jsb_time,                ONLY: timesteps_per_day
    USE mo_jsb_model_init,          ONLY: jsbach_setup_grid, jsbach_init, jsbach_init_after_restart
    USE mo_jsb_model_final,         ONLY: jsbach_finalize
    USE mo_jsb_version,             ONLY: jsbach_label_run
    USE mo_jsb_control,             ONLY: get_no_of_models, force_from_observations
    USE mo_jsb4_forcing,            ONLY: setup_forcing, init_forcing, finalize_external_forcing

    USE mo_build_decomposition,     ONLY: build_decomposition

    ! I/O
    USE mo_restart,                 ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor, detachRestartProcs
    USE mo_name_list_output,        ONLY: name_list_io_main_proc, write_name_list_output, close_name_list_output, &
      &                                   istime4name_list_output
    USE mo_output_event_handler,    ONLY: get_current_jfile
    USE mo_icon_output_tools,       ONLY: init_io_processes
#ifdef HAVE_CDI_PIO
    USE mo_mpi,                     ONLY: mpi_comm_null, p_comm_work_io
    USE mo_name_list_output_init,   ONLY: init_cdipio_cb
    USE mo_name_list_output,        ONLY: write_ready_files_cdipio
    USE yaxt,                       ONLY: xt_initialize, xt_initialized
    USE mo_cdi,                     ONLY: namespacegetactive, namespaceSetActive
    USE mo_cdi_pio_interface,       ONLY: nml_io_cdi_pio_namespace, &
      &                                   cdi_base_namespace, &
      &                                   nml_io_cdi_pio_client_comm, &
      &                                   nml_io_cdi_pio_conf_handle
#endif
    USE mo_name_list_output_config, ONLY: use_async_name_list_io
    USE mo_output_event_types,      ONLY: t_sim_step_info
    USE mo_impl_constants,          ONLY: SUCCESS, pio_type_async, pio_type_cdipio
    USE mo_memory_log,              ONLY: memory_log_terminate
    USE mo_util_mtime,              ONLY: getElapsedSimTimeInSeconds, is_event_active

    USE mo_grid_config,             ONLY: n_dom, n_dom_start
    USE mo_model_domain,            ONLY: p_patch
    USE mo_model_domain,            ONLY: p_patch_local_parent
    USE mo_icon_comm_interface,     ONLY: construct_icon_communication, destruct_icon_communication
    USE mo_complete_subdivision,    ONLY: setup_phys_patches
    USE mo_gribout_config,          ONLY: configure_gribout

    USE mo_vertical_coord_table,    ONLY: vct_a, vct_b, vct, allocate_vct_atmo
    ! Horizontal interpolation
    USE mo_intp_lonlat_types,       ONLY: lonlat_grids
    USE mo_intp_data_strc,          ONLY: p_int_state, p_int_state_local_parent
    USE mo_interpol_config,         ONLY: configure_interpolation
    USE mo_intp_state,              ONLY: construct_2d_interpol_state,                          &
      &                                   destruct_2d_interpol_state !, transfer_interpol_state
    USE mo_intp_lonlat,             ONLY: compute_lonlat_intp_coeffs
    USE mo_gridref_config,          ONLY: configure_gridref

    USE mo_alloc_patches,           ONLY: destruct_patches, destruct_comm_patterns
    USE mo_action,                  ONLY: ACTION_RESET, reset_act
    USE mo_var_list_register_utils, ONLY: vlr_print_groups, vlr_print_vls

    USE mo_jsb_convect_tables_iface,ONLY: init_convect_tables

    USE mo_derived_variable_handling, ONLY: update_statistics
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(LEN=*), INTENT(in) :: jsb_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER   :: model

    TYPE(t_key_value_store),      POINTER :: restartAttributes
    CLASS(t_RestartDescriptor),   POINTER :: restartDescriptor
    INTEGER :: jstep0, jstep, dedicatedRestartProcs
    INTEGER :: no_models, jg
    INTEGER :: i
    CHARACTER(len=MAX_MTIME_ERROR_STR_LEN) :: errstring
    INTEGER                                :: error_status

#ifdef HAVE_CDI_PIO
    INTEGER :: prev_cdi_namespace
#endif

    ! Time and event management
    TYPE(datetime), POINTER       :: mtime_current
    TYPE(datetime), POINTER       :: mtime_old
    REAL(wp) :: sim_time          ! elapsed simulation time
    TYPE(datetime),  POINTER :: eventStartDate    => NULL(), &
      &                         eventEndDate      => NULL()
    TYPE(datetime),  POINTER :: checkpointRefDate => NULL(), &
      &                         restartRefDate    => NULL()
    TYPE(event), POINTER                   :: checkpointEvent   => NULL()
    TYPE(event), POINTER                   :: restartEvent      => NULL()
    INTEGER                                :: checkpointEvents
    TYPE(eventGroup), POINTER              :: checkpointEventGroup => NULL()
    TYPE(timedelta), POINTER               :: eventInterval     => NULL()
    TYPE(timedelta), POINTER               :: model_time_step   => NULL()
    TYPE(t_sim_step_info) :: sim_step_info
    LOGICAL                                :: lret, lwrite_checkpoint
    LOGICAL                                :: l_isStartdate, l_isRestart, l_isExpStopdate, l_isCheckpoint, l_doWriteRestart
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: dt_string, dstring
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)   :: td_string
    LOGICAL                                :: lprint_timestep, l_nml_output

    INTEGER, ALLOCATABLE :: output_jfile(:) ! number of output files

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_model'
    ! ----------------------------------------------------------------------------------------------------- !

    CALL message('======== ', 'Starting the ICON-Land standalone model ...')
    CALL message('', '')
    CALL message('Namelist files', TRIM(jsb_namelist_filename) // ', ' // TRIM(shr_namelist_filename))
    CALL message('', '')

    IF (n_dom > 1) CALL finish(routine, 'Only one atmospheric domain allowed')

    ! ----------------------------------------------------------------------------------------------------- !
    !> 0.1 Initialize global registry of lon-lat grids
    !>
    CALL lonlat_grids%init()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 0.2 If this is a resumed or warm-start run...
    !>
    restartAttributes => NULL()
    IF (isRestart()) THEN
      CALL message('','Read restart file meta data ...')
      ! Reads attributes and namelists for all available domains from restart file.
      ! These namelist values will overwrite the model default, and will later be overwritten if the user has
      ! specified something different for this integration.
      CALL read_restart_header(TRIM(get_my_process_name()))
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.0 read namelists
    !>

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.1 read namelists required for the ICON infrastructure and fill the
    !>     corresponding sections of the ICON configuration states
    !>
    ! Create a new file in which all the namelist variables and their
    ! actual values used in the model run will be stored.
    IF (my_process_is_stdio()) CALL open_nml_output(TRIM(jsb_namelist_filename)//'_output')

    CALL read_infrastructure_namelists_for_jsbach(TRIM(jsb_namelist_filename), TRIM(shr_namelist_filename))

    ! The following can't be called from the previous subroutine since it's causing a circular dependency.
    CALL read_name_list_output_namelists (TRIM(jsb_namelist_filename))

    ! Close the file in which all the namelist variables and their
    ! actual values were stored and write an annotated table of all namelist settings to a text file
    IF (my_process_is_stdio()) THEN
      CALL close_nml_output()
      CALL log_nml_settings(TRIM(jsb_namelist_filename)//".log")
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.2 Cross-check namelist setups
    !>
    CALL check_parallel_configuration()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 2.0 Set times and events
    !>
    CALL compute_timestep_settings()
    CALL compute_restart_settings()
    CALL compute_date_settings("jsbach", dt_restart, nsteps)
    CALL initEventManager(time_config%tc_exp_refdate)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 3.0 Call configure_run to finish filling the run_config state.
    !>
    !>    This needs to be done very early (but anyway after atm_crosscheck)
    !>    because some components of the state, e.g., num_lev, may be
    !>    modified in this subroutine which affects the following CALLs.
    !>
    CALL configure_run( )

    ! ----------------------------------------------------------------------------------------------------- !
    !> 4.0 intialization
    !>

    ! ----------------------------------------------------------------------------------------------------- !
    !> 4.1 Initialize the mpi work groups
    !>
    num_prefetch_proc = 0
    CALL restartWritingParameters(opt_dedicatedProcCount = dedicatedRestartProcs)
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, &
         &                          num_io_procs, dedicatedRestartProcs, &
         &                          jsbach_process, num_prefetch_proc, num_test_pe,      &
         &                          pio_type, num_dio_procs=0)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 4.2 Initialize various timers
    !>
    IF (ltimer) CALL init_timer
    IF (timers_level > 1) CALL timer_start(timer_model_init)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 4.3 initialize dynamic list of vertical axes
    !>
    zaxisTypeList = t_zaxisTypeList()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 5.0 Setup the model in ICON-Land: read namelists, configure models for each domain
    !>
    !>   @NOTE This has to be after (!) the ICON zaxes have been created in the above line but
    !>         before (!) the restart PEs are detached a few lines below since ICON-Land
    !>         adds its zaxes to zaxisTypeList
    !>
    CALL jsbach_setup_models(shr_namelist_filename)
    no_models = get_no_of_models()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 6.0 I/O initialization
    !>
    ! This won't RETURN on dedicated restart PEs, starting their main loop instead.
    CALL detachRestartProcs(ltimer)
    CALL init_io_processes()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 7.0 Import patches, perform domain decomposition
    !>
    !>   define the horizontal and vertical grids since they are aready
    !>   needed for some derived control parameters. This includes
    !>   - patch import
    !>   - domain decompistion
    !>   - vertical coordinates
    !>
    IF (timers_level > 4) CALL timer_start(timer_domain_decomp)
    CALL build_decomposition(num_lev, nshift, is_ocean_decomposition = .false.)
    IF (timers_level > 4) CALL timer_stop(timer_domain_decomp)

    p_patch(:)%ldom_active = .TRUE.

    ! ----------------------------------------------------------------------------------------------------- !
    !> 8.0 Construct interpolation state, compute interpolation coefficients.
    !>
    IF (timers_level > 4) CALL timer_start(timer_compute_coeffs)
    CALL configure_interpolation( n_dom, p_patch(1:)%level, &
                                  p_patch(1:)%geometry_info )

    CALL configure_gridref(n_dom, p_patch(1:)%geometry_info%mean_characteristic_length)

    ! Allocate array for interpolation state
    ALLOCATE( p_int_state(n_dom_start:n_dom), &
            & STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(routine, 'allocation for ptr_int_state failed')
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 9.0 Initialize icon_comm_lib
    !>
    CALL construct_icon_communication(p_patch, n_dom)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 10.0 Setup the information for the physical patches
    !>
    CALL setup_phys_patches

    ! Copy needed parts of the patches to GPU
    DO jg=1,no_models
      !$ACC ENTER DATA &
      !$ACC   COPYIN(p_patch(jg)%cells%decomp_info%owner_mask)
    END DO

    ! Constructing data for lon-lat interpolation
    ! CALL compute_lonlat_intp_coeffs(p_patch(1:), p_int_state(1:)) ! inout, inout

    CALL allocate_vct_atmo(p_patch(1)%nlevp1)

    IF (timers_level > 4) CALL timer_stop(timer_compute_coeffs)

    CALL init_convect_tables ! necessary for FUNCTION qsat

    ! ----------------------------------------------------------------------------------------------------- !
    !> 11.0 setup and init external forcing
    !>
    IF (force_from_observations()) CALL setup_forcing(no_models)  ! @TODO this might always be TRUE with icon-land standalone
    DO jg=1,no_models
      model => Get_model(jg)
      CALL jsbach_setup_grid(jg, p_patch(jg), type='icon')
      CALL jsbach_setup_tiles(jg)
      CALL jsbach_init(jg)
      IF (force_from_observations()) CALL init_forcing(jg, model%namelist_filename, dtime, dtime)
    END DO
    CALL jsbach_label_run(.TRUE.)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 12.0 restart a simulation
    !>
    IF (isRestart()) THEN
      !
      ! This is a resumed integration. Read model state from restart file(s).
      !
      IF (timers_level > 4) CALL timer_start(timer_read_restart)
      DO jg = 1,n_dom
        IF (p_patch(jg)%ldom_active) THEN
          CALL read_restart_files( p_patch(jg), n_dom)
        END IF
      END DO
      CALL message(TRIM(routine),'normal exit from read_restart_files')

      IF (timers_level > 4) CALL timer_stop(timer_read_restart)

      DO jg=1,no_models
        CALL jsbach_init_after_restart(jg)
      END DO

    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 13.0 Prepare output
    !>
    CALL configure_gribout(grid_generatingCenter, grid_generatingSubcenter, n_dom)
    CALL configure_io()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 13.1 Map the variable groups given in the output namelist onto the corresponding variable subsets:
    !>
    !>   @NOTE all add_vars must be finished before calling this routine.
    !>
    IF (output_mode%l_nml) THEN
      CALL parse_variable_groups()
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 13.2 If async IO is in effect, init_name_list_output is a collective call
    !>      with the IO procs and effectively starts async IO
    !>
    IF (output_mode%l_nml) THEN
      ! compute sim_start, sim_end
      sim_step_info%sim_start    = time_config%tc_exp_startdate
      sim_step_info%sim_end      = time_config%tc_exp_stopdate
      sim_step_info%run_start    = time_config%tc_startdate
      sim_step_info%restart_time = time_config%tc_stopdate

      sim_step_info%dtime  = dtime
      sim_step_info%jstep0 = 0

      CALL getAttributesForRestarting(restartAttributes)
      ! get start counter for time loop from restart file:
      IF (restartAttributes%is_init) &
        & CALL restartAttributes%get("jstep", sim_step_info%jstep0)
      CALL init_name_list_output(sim_step_info)
      CALL create_mipz_level_selections(output_file)
      CALL create_vertical_axes(output_file)

    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 13.3 Initialize reset-Action, i.e. assign variables to action object
    !>
    CALL reset_act%initialize(ACTION_RESET)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 13.4 Optional: detailed print-out of some variable info  !
    !>
    IF (my_process_is_stdio() .AND. (msg_level >= 15)) THEN
      DO jg=1,n_dom
        CALL vlr_print_groups(idom=jg)
      END DO
    END IF

    IF (timers_level > 1) CALL timer_stop(timer_model_init)

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> End of model initialization
    !>
    ! ----------------------------------------------------------------------------------------------------- !

#if defined( _OPENACC )
    i_am_accel_node = my_process_is_work()
#endif

    ! ----------------------------------------------------------------------------------------------------- !
    !> 14.0 get and write out some of the initial values
    !>
    mtime_current => time_config%tc_current_date
    IF (.NOT. isRestart() .AND. (mtime_current >= time_config%tc_exp_startdate)) THEN
      ! ..
      CALL update_statistics
      ! --
      IF (output_mode%l_nml) THEN
          CALL write_name_list_output(jstep=0, lacc=i_am_accel_node)
      END IF

    END IF ! not isRestart()

    IF (ltimer) CALL timer_start(timer_total)

    ! calculate elapsed simulation time in seconds
    sim_time = getElapsedSimTimeInSeconds(mtime_current)

    ! allocate temporary variable for restarting purposes
    IF (output_mode%l_nml) THEN
      ALLOCATE(output_jfile(SIZE(output_file)), STAT=error_status)
      IF (error_status /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed!')
    ENDIF

    mtime_old     => newDatetime(mtime_current) ! in

    restartDescriptor => createRestartDescriptor("jsbach")
    jstep0 = 0
    CALL getAttributesForRestarting(restartAttributes)
    IF (isRestart()) THEN
      ! get start counter for time loop from restart file:
      CALL restartAttributes%get("jstep", jstep0)
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 14.1 for debug purposes print var lists: for msg_level >= 13 short and for >= 20 long format
    !>
    IF  (msg_level >= 13) CALL vlr_print_vls(lshort=(msg_level < 20))

    ! ----------------------------------------------------------------------------------------------------- !
    !> 15.0 set events, group and the events
    !>
    CALL message('','')
    eventStartDate => time_config%tc_exp_startdate
    eventEndDate   => time_config%tc_exp_stopdate
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
    ! --- create an event group for checkpointing and restart
    checkpointEvents =  addEventGroup('checkpointEventGroup')
    checkpointEventGroup => getEventGroup(checkpointEvents)
    ! --- --- create checkpointing event
    eventInterval  => time_config%tc_dt_checkpoint
    checkpointEvent => newEvent('checkpoint', checkpointRefDate, eventStartDate, eventEndDate, eventInterval, errno=error_status)
    IF (error_status /= NO_ERROR) THEN
      ! give an elaborate error message:
      CALL datetimeToString(checkpointRefDate, dt_string)
      WRITE (0,*) "event reference date: ",    dt_string
      CALL datetimeToString(eventStartDate,    dt_string)
      WRITE (0,*) "event start date    : ",    dt_string
      CALL datetimeToString(eventEndDate,      dt_string)
      WRITE (0,*) "event end date      : ",    dt_string
      CALL timedeltaToString(eventInterval,    td_string)
      WRITE (0,*) "event interval      : ",    td_string
      CALL mtime_strerror(error_status, errstring)
      CALL finish('perform_nh_timeloop', "event 'checkpoint': "//errstring)
    ENDIF
    lret = addEventToEventGroup(checkpointEvent, checkpointEventGroup)
    ! --- --- create restart event, ie. checkpoint + model stop
    eventInterval  => time_config%tc_dt_restart
    restartEvent => newEvent('restart', restartRefDate, eventStartDate, eventEndDate, eventInterval, errno=error_status)
    IF (error_status /= NO_ERROR) THEN
      ! give an elaborate error message:
      CALL datetimeToString(restartRefDate, dt_string)
      WRITE (0,*) "event reference date: ", dt_string
      CALL datetimeToString(eventStartDate, dt_string)
      WRITE (0,*) "event start date    : ", dt_string
      CALL datetimeToString(eventEndDate,   dt_string)
      WRITE (0,*) "event end date      : ", dt_string
      CALL timedeltaToString(eventInterval, td_string)
      WRITE (0,*) "event interval      : ", td_string
      CALL mtime_strerror(error_status, errstring)
      CALL finish('perform_nh_timeloop', "event 'restart': "//errstring)
    ENDIF
    lret = addEventToEventGroup(restartEvent, checkpointEventGroup)
    CALL printEventGroup(checkpointEvents)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 16.0 set time loop properties
    !>
    model_time_step => time_config%tc_dt_model

    CALL message('','')
    CALL datetimeToString(mtime_current, dstring)
    WRITE(message_text,'(a,a)') 'Start date of this run: ', dstring
    CALL message('',message_text)
    CALL datetimeToString(time_config%tc_stopdate, dstring)
    WRITE(message_text,'(a,a)') 'Stop date of this run:  ', dstring
    CALL message('',message_text)
    CALL message('','')

    jstep = jstep0+1

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 17.0 time loop for running ICON-Land
    !>
    TIME_LOOP: DO
      ! optional memory loggin
      CALL memory_log_add

      !TODO Compute mtime_forcing
      mtime_current = mtime_current + model_time_step

      ! store state of output files for restarting purposes
      IF (output_mode%l_nml .AND. jstep>=0 ) THEN
        DO i=1,SIZE(output_file)
          output_jfile(i) = get_current_jfile(output_file(i)%out_event)
        END DO
      ENDIF

      IF (msg_level > 2) THEN
        lprint_timestep = .TRUE.
      ELSE
        lprint_timestep = MOD(jstep,timesteps_per_day(dtime,1)) == 0
      ENDIF
      ! always print the first and the last time step
      lprint_timestep = lprint_timestep .OR. (jstep == jstep0+1) .OR. (jstep == jstep0+nsteps)

      IF (lprint_timestep) THEN
        WRITE(message_text,'(a,i8,a,i0,a,5(i2.2,a),i3.3)') &
              &             'Time step: ', jstep, ' model time ',                                &
              &             mtime_current%date%year,   '-', mtime_current%date%month,    '-',    &
              &             mtime_current%date%day,    ' ', mtime_current%time%hour,     ':',    &
              &             mtime_current%time%minute, ':', mtime_current%time%second,   '.',    &
              &             mtime_current%time%ms
        CALL message('',message_text)
      ENDIF

      ! set namelist output flag
      l_nml_output = output_mode%l_nml .AND. jstep >= 0 .AND. istime4name_list_output(jstep)

      ! --------------------------------------------------------------------------------------------------- !
      !> 17.1 run ICON-Land for one timestep
      !>
      DO jg = 1,no_models
        CALL run_one_timestep( &
          & jg, &                   ! model
          & p_patch(jg)%nblks_c, &  ! ..
          & dtime, &                ! length of timestep in seconds
          & mtime_current)          ! ..
      END DO

      ! update accumulated values
      CALL update_statistics

      IF (l_nml_output) THEN
        CALL write_name_list_output(jstep, lacc=i_am_accel_node)
      ENDIF

      ! re-initialize MAX/MIN fields with 'resetval'
      !   @NOTE must be done AFTER output
      CALL reset_act%execute(slack=dtime, mtime_date=mtime_current)

      ! --------------------------------------------------------------------------------------------------- !
      !> 17.2 Write restart file
      !>
      !>  check whether time has come for writing restart file
      !>
      ! default is to assume we do not write a checkpoint/restart file
      lwrite_checkpoint = .FALSE.
      ! if the model is not supposed to write output, do not write checkpoints
      IF (.NOT. output_mode%l_none ) THEN
        ! to clarify the decision tree we use shorter and more expressive names:
        l_isStartdate    = (time_config%tc_startdate == mtime_current)
        l_isExpStopdate  = (time_config%tc_exp_stopdate == mtime_current)
        l_isRestart      = is_event_active(restartEvent, mtime_current, proc0_offloading)
        l_isCheckpoint   = is_event_active(checkpointEvent, mtime_current, proc0_offloading)
        l_doWriteRestart = time_config%tc_write_restart
        ! ...
        IF ( &
              !  if normal checkpoint or restart cycle has been reached, i.e. checkpoint+model stop
              &         (l_isRestart .OR. l_isCheckpoint)                     &
              &  .AND.                                                        &
              !  and the current date differs from the start date
              &        .NOT. l_isStartdate                                    &
              &  .AND.                                                        &
              !  and end of run has not been reached or restart writing has been disabled
              &        (.NOT. l_isExpStopdate .OR. l_doWriteRestart)          &
              & ) THEN
          lwrite_checkpoint = .TRUE.
        END IF
      END IF
      ! ..
      IF (lwrite_checkpoint) THEN
        DO jg=1,n_dom
          CALL restartDescriptor%updatePatch(p_patch(jg), opt_ndom=n_dom)
        END DO
        ! finally write restart file
        CALL restartDescriptor%writeRestart(mtime_current, jstep, opt_output_jfile = output_jfile)
      END IF  ! lwrite_checkpoint
      ! ...
      IF (mtime_current >= time_config%tc_stopdate) THEN
        ! leave time loop
        EXIT TIME_LOOP
      END IF

      ! increase timestep by one
      jstep = jstep + 1
      ! step_in_forcing_file = step_in_forcing_file + 1
      sim_time = getElapsedSimTimeInSeconds(mtime_current)
    END DO TIME_LOOP

#if defined( _OPENACC )
    i_am_accel_node = .FALSE.
#endif

    ! ----------------------------------------------------------------------------------------------------- !
    !> 18.0  ... @TODO describe
    !>
    CALL deleteRestartDescriptor(restartDescriptor)

    ! ..
    IF (ltimer) CALL timer_stop(timer_total)

    ! ..
    IF (output_mode%l_nml) THEN
      DEALLOCATE(output_jfile, STAT=error_status)
      IF (error_status /= SUCCESS)  CALL finish (routine, 'DEALLOCATE failed!')
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 19.0 deallocation of structure and variables
    !>
    ! @TODO CALL deallocateDatetime(mtime_forcing)
    CALL deallocateDatetime(mtime_old)

    ! Deallocate global registry for lon-lat grids
    CALL lonlat_grids%finalize()

    ! Destruct communication patterns
    CALL destruct_comm_patterns( p_patch, p_patch_local_parent )

    ! Deallocate grid patches
    DO jg=1,no_models
      !$ACC EXIT DATA &
      !$ACC   DELETE(p_patch(jg)%cells%decomp_info%owner_mask)
    END DO
    CALL destruct_patches( p_patch )

    ! ..
    DEALLOCATE( p_patch, STAT=error_status )
    IF (error_status/=SUCCESS) THEN
      CALL finish(routine, 'deallocate for patch array failed')
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 20.0 close memory logging files and finalize model
    !>
    CALL memory_log_terminate
    CALL destruct_icon_communication()
    CALL finalize_external_forcing(no_models)
    CALL jsbach_finalize()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 21.0 Delete variable lists
    !>
    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
    ENDIF
#ifdef HAVE_CDI_PIO
    IF (pio_type == pio_type_cdipio) THEN
      prev_cdi_namespace = namespaceGetActive()
      CALL namespaceSetActive(nml_io_cdi_pio_namespace)
      CALL pioFinalize
      CALL namespaceSetActive(prev_cdi_namespace)
    END IF
#endif

    ! ----------------------------------------------------------------------------------------------------- !
    !> 22.0 finish message and print performance timers
    !>
    CALL message(routine, 'clean-up finished')
    IF (ltimer) CALL print_timer

  END SUBROUTINE jsbach_model

  ! ======================================================================================================= !
  !>
  !> run ICON-Land for one timestep
  !>
  SUBROUTINE run_one_timestep(model_id, nblks, dtime, datetime_new)

    USE mo_jsb_control,               ONLY: jsbach_is_restarted, force_from_observations, &
      &                                     debug_on, timer_on, l_timer_host, timer_jsbach, timer_forcing
    USE mo_jsb_time,                  ONLY: is_time_experiment_start, is_time_restart, get_secs_of_day, get_year, is_newyear
    USE mo_jsb_model_class,           ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
    USE mo_jsb_class,                 ONLY: get_model
    USE mo_jsb_grid,                  ONLY: Get_grid
    USE mo_jsb_grid_class,            ONLY: t_jsb_grid
    USE mo_jsb_tile_class,            ONLY: t_jsb_tile_abstract
    USE mo_jsb_interface,             ONLY: jsbach_start_timestep, jsbach_finish_timestep
    USE mo_jsb4_forcing,              ONLY: forcing_options, get_standalone_driver
    USE mo_orbit_solar,               ONLY: compute_cos_zenith_angle
    USE mo_jsb_physical_constants,    ONLY: molarMassDryAir, molarMassCO2
    USE mo_jsb_parallel,              ONLY: Get_omp_thread
    USE mo_jsb_subset,                ONLY: ON_CHUNK

#ifndef __NO_QUINCY__
    USE mo_iq_atm2land_process,       ONLY: update_local_time_and_daytime_counter, update_slow_sb_pool_accelerator_bookkeeping
    dsl4jsb_Use_processes A2L_, HYDRO_, SPQ_
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_config(SPQ_)
#else
    dsl4jsb_Use_processes A2L_, HYDRO_
#endif
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                INTENT(in)            :: model_id
    INTEGER,                INTENT(in)            :: nblks
    REAL(wp)              , INTENT(in)            :: dtime          !< time step
    TYPE(datetime)        , POINTER               :: datetime_new   !< date and time at the end of this time step
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=max_timedelta_str_len) :: dtime_string     !< time delta as string
    TYPE(timedelta)         , POINTER    :: dtime_mtime      !< time delta as mtime variable
    TYPE(datetime)          , POINTER    :: datetime_old     !< date and time at the beginning of this time step
    TYPE(datetime)          , POINTER    :: datetime_old_old !< date and time one time step before datetime_old
    TYPE(datetime)          , POINTER    :: datetime_forcing

    INTEGER  :: nproma
    INTEGER  :: iblk                           !< index of block loop
    INTEGER  :: startblk, endblk               !< start and end indices of block   loop
    INTEGER  :: startidx(nblks), endidx(nblks) !< start and end indices of columns loop

    TYPE(t_jsb_model), POINTER           :: model
    TYPE(t_jsb_grid),  POINTER           :: grid       ! Horizontal grid
    CLASS(t_jsb_tile_abstract), POINTER  :: tile

    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Real2D_onDomain :: evapotrans
    dsl4jsb_Real2D_onDomain :: evapopot
    dsl4jsb_Real2D_onDomain :: co2_air
    dsl4jsb_Real2D_onDomain :: co2_air_mol
#ifndef __NO_QUINCY__
    ! quincy
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_config(SPQ_)

    dsl4jsb_Real2D_onDomain :: swpar_srf_down
    dsl4jsb_Real2D_onDomain :: daytime_counter
    dsl4jsb_Real2D_onDomain :: local_time_day_seconds
    dsl4jsb_Real2D_onDomain :: slow_sb_pool_accelerator_execution_counter
    dsl4jsb_Real2D_onDomain :: slow_sb_pool_accelerator_execute

#endif

    REAL(wp), POINTER     :: lon(:,:), lat(:,:), coslat(:,:), sinlat(:,:)
    REAL(wp), ALLOCATABLE :: evapo_act2pot(:,:)

    LOGICAL :: use_external_forcing               ! use forcing from external forcing (e.g., observation based data)
    LOGICAL :: l_first_forcing_file    = .TRUE.
    LOGICAL :: l_shift_one_timestep    = .TRUE.
    LOGICAL :: l_first_step_in_restart = .FALSE.

    INTEGER :: nc, ics, ice, ic
    INTEGER :: global_seconds_day, current_year

    LOGICAL :: save_i_am_accel_node

    ! IQ spin-up accelerator model configurations
    INTEGER :: accelerator_max_executions !< bookkeeping configurations for the iq spin-up accelerator: max number of executions
    INTEGER :: accelerator_frequency      !< bookkeeping configurations for the iq spin-up accelerator: frequency of executions
    INTEGER :: accelerator_start_year     !< bookkeeping configurations for the iq spin-up accelerator: start year of executions
    ! IQ: further configuration
    LOGICAL :: flag_snow = .TRUE. !< quincy boolean indicating if simulations should include snow (else read snow is added to rain)

    CHARACTER(len=*), PARAMETER :: routine = modname//':run_one_timestep'
    ! ----------------------------------------------------------------------------------------------------- !

    ! use forcing from external forcing, such as, observation based data (e.g., CRUJRA)
    !                                             or output from another model simulation
    ! see namelist 'jsb_control_nml' and 'mo_jsb_control:l_force_from_obs'
    use_external_forcing = force_from_observations()

    IF (jsbach_is_restarted()) THEN
      l_first_forcing_file = .FALSE.
      l_shift_one_timestep = .FALSE.
    ELSE
      l_shift_one_timestep = .TRUE.
    END IF

    CALL getPTStringFromSeconds(-dtime, dtime_string)
    dtime_mtime      => newTimedelta(dtime_string)
    datetime_old     => newDatetime(datetime_new)
    datetime_old     =  datetime_new + dtime_mtime
    datetime_old_old => newDatetime(datetime_old)
    datetime_old_old =  datetime_old + dtime_mtime
    CALL deallocateTimedelta(dtime_mtime)

    IF (is_time_restart(datetime_old)) THEN
      l_first_step_in_restart = .TRUE.
    END IF

    datetime_forcing => newDatetime(datetime_old)
    IF (l_first_forcing_file) THEN
      IF (forcing_options(model_id)%forcing_steps_per_day == 1  .AND. .NOT. use_external_forcing ) THEN
        ! Go to next time step of next file
        CALL getPTStringFromSeconds(forcing_options(model_id)%forcing_synchron_factor * dtime, dtime_string)
        dtime_mtime  => newTimedelta(dtime_string)
        datetime_forcing  =  datetime_forcing + dtime_mtime
        CALL deallocateTimedelta(dtime_mtime)
      END IF
    ELSE
      IF (.NOT. l_first_step_in_restart  .AND. .NOT. use_external_forcing) THEN
        datetime_forcing = datetime_new
      END IF
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> start timestep
    !>
    CALL jsbach_start_timestep(model_id, datetime_old, dtime)

    model => get_model(model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()

    CALL model%Get_top_tile(tile)

    startblk = grid%Get_blk_start()
    endblk   = grid%Get_blk_end()
    DO iblk = startblk, endblk
      startidx(iblk) = grid%Get_col_start(iblk)
      endidx  (iblk) = grid%Get_col_end(iblk)
    END DO
    !$ACC DATA COPYIN(startidx, endidx) ASYNC(1)

    IF (use_external_forcing) THEN

      dsl4jsb_Get_memory(A2L_)

      lon => grid%lon
      lat => grid%lat
      coslat => grid%coslat
      sinlat => grid%sinlat

      dsl4jsb_Get_var2d_onDomain(A2L_,   co2_air_mol)
      dsl4jsb_Get_var2d_onDomain(A2L_,   co2_air)

      CALL compute_cos_zenith_angle(datetime_old, grid%patch, dsl4jsb_var2D_onDomain(A2L_, cos_zenith_angle))
      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(A2L_, cos_zenith_angle)) ASYNC(1)
      ! Note: inquire_declination should only be called each time step after compute_cos_zenith_angle

      ! TODO: if really required calculate daily sums of evapotrans and evapopot in jsbach and get them here (see comment below)
      ! -> would only be required if no Qair forcing is given and qair needs to be calculated from temperature!
      ALLOCATE(evapo_act2pot(nproma, nblks))
      !$ACC ENTER DATA CREATE(evapo_act2pot)
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        dsl4jsb_Get_memory(HYDRO_)
        dsl4jsb_Get_var2d_onDomain(HYDRO_, evapotrans)
        dsl4jsb_Get_var2d_onDomain(HYDRO_, evapopot)

        IF (is_time_experiment_start(datetime_old)) THEN

          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO iblk = 1, nblks
            DO ic = 1, nproma
              evapo_act2pot(ic, iblk) = 1._wp
            ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

        ELSE

          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO iblk = 1, nblks
            DO ic = 1, nproma
              evapo_act2pot(ic, iblk) =   MIN(evapotrans(ic, iblk), -EPSILON(1._wp)) &
                &                       / MIN(evapopot  (ic, iblk), -EPSILON(1._wp))
            ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

        END IF
      CASE (MODEL_QUINCY)
        ! Note: if neither qair nor rh are given as forcing (i.e. forcing_options(model_id)%type_of_qair_forcing = NONE,
        !       which so far never is the case), this would need to be calculated
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO iblk = 1, nblks
          DO ic = 1, nproma
            evapo_act2pot(ic, iblk) = 1._wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      END SELECT

      IF (timer_on() .OR. l_timer_host) CALL timer_start(timer_forcing(model%id))

      ! --------------------------------------------------------------------------------------------------- !
      !> get variables (that are in the coupled simulation provided via the model interface by the host model)
      !>   including dimensional data and data from external forcing, i.e., forcing input files
      !>   see: 'mo_jsb_interface:interface_full'
      !>

      ! We need to turn off i_am_accel_node for the forcing module to do host-based file reading
      ! This does not affect the OpenACC clauses themselves in the forcing module!
      save_i_am_accel_node = i_am_accel_node
      i_am_accel_node = .FALSE.

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        CALL get_standalone_driver( &
          ! in
          & model_id,                                                                  &
          & model%config%model_scheme,                                                 &
          & flag_snow,                                                                 &
          & startblk, endblk, startidx, endidx, datetime_old, datetime_new,            &
          & dsl4jsb_var2D_onDomain(HYDRO_, elevation),                                 &
          & sinlat, coslat, lon, evapo_act2pot,                                        &
          & dsl4jsb_var2D_onDomain(A2L_,  cos_zenith_angle),                           &
          ! out
          & dsl4jsb_var2D_onDomain(A2L_, CO2_air),                                     &
          & dsl4jsb_var2D_onDomain(A2L_, t_air),                                       &
          & dsl4jsb_var2D_onDomain(A2L_, q_air),                                       &
          & dsl4jsb_var2D_onDomain(A2L_, rain),                                        &
          & dsl4jsb_var2D_onDomain(A2L_, snow),                                        &
          & dsl4jsb_var2D_onDomain(A2L_, wind_air),                                    &
          & dsl4jsb_var2D_onDomain(A2L_, wind_10m),                                    &
          & dsl4jsb_var2D_onDomain(A2L_, lw_srf_down),                                 &
          & dsl4jsb_var2D_onDomain(A2L_, swvis_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, swnir_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, swpar_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, fract_par_diffuse),                           &
          & dsl4jsb_var2D_onDomain(A2L_, press_srf)                                    &
          & )
#ifndef __NO_QUINCY__
      CASE (MODEL_QUINCY)
        dsl4jsb_Get_memory(SPQ_)
        dsl4jsb_Get_config(SPQ_)
        flag_snow = dsl4jsb_Config(SPQ_)%flag_snow

        CALL get_standalone_driver( &
          ! in
          & model_id,                                                                  &
          & model%config%model_scheme,                                                 &
          & flag_snow,                                                                 &
          & startblk, endblk, startidx, endidx, datetime_old, datetime_new,            &
          & dsl4jsb_var2D_onDomain(SPQ_, elevation),                                   &
          & sinlat, coslat, lon, evapo_act2pot,                                        &
          & dsl4jsb_var2D_onDomain(A2L_,  cos_zenith_angle),                           &
          ! out
          & dsl4jsb_var2D_onDomain(A2L_, CO2_air),                                     &
          & dsl4jsb_var2D_onDomain(A2L_, t_air),                                       &
          & dsl4jsb_var2D_onDomain(A2L_, q_air),                                       &
          & dsl4jsb_var2D_onDomain(A2L_, rain),                                        &
          & dsl4jsb_var2D_onDomain(A2L_, snow),                                        &
          & dsl4jsb_var2D_onDomain(A2L_, wind_air),                                    &
          & dsl4jsb_var2D_onDomain(A2L_, wind_10m),                                    &
          & dsl4jsb_var2D_onDomain(A2L_, lw_srf_down),                                 &
          & dsl4jsb_var2D_onDomain(A2L_, swvis_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, swnir_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, swpar_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, fract_par_diffuse),                           &
          & dsl4jsb_var2D_onDomain(A2L_, press_srf),                                   &
          & dsl4jsb_var2D_onDomain(A2L_, nhx_deposition),                              &
          & dsl4jsb_var2D_onDomain(A2L_, noy_deposition),                              &
          & dsl4jsb_var2D_onDomain(A2L_, nhx_n15_deposition),                          &
          & dsl4jsb_var2D_onDomain(A2L_, noy_n15_deposition),                          &
          & dsl4jsb_var2D_onDomain(A2L_, p_deposition)                                 &
          & )
#endif
      END SELECT

      ! Turn it on again for the rest of the time step
      i_am_accel_node = save_i_am_accel_node

      IF (timer_on() .OR. l_timer_host) CALL timer_stop(timer_forcing(model%id))

      ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO iblk = 1, nblks
        DO ic = 1, nproma
          co2_air_mol(ic, iblk) = co2_air(ic, iblk) * molarMassDryAir / molarMassCO2
        ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP

#ifndef __NO_QUINCY__
      ! --------------------------------------------------------------------------------------------------- !
      !> code specific to the QUINCY model
      !>
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)
#ifdef _OPENACC
        CALL finish(routine, 'QUINCY standalone not ported to GPU, yet. Stop.')
#endif
        ! convert CO2 from "molar ratio (volume)" to "co2 mixing ratio ppmv"
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO iblk = 1, nblks
          DO ic = 1, nproma
            dsl4jsb_var_ptr(A2L_, CO2_mixing_ratio)(ic, iblk) = dsl4jsb_var_ptr(A2L_, CO2_air_mol)(ic, iblk) * 1000000._wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

        ! Update local time and daytime counter
        dsl4jsb_Get_var2d_onDomain(A2L_, daytime_counter)
        dsl4jsb_Get_var2d_onDomain(A2L_, local_time_day_seconds)
        dsl4jsb_Get_var2d_onDomain(A2L_, swpar_srf_down)
        global_seconds_day = get_secs_of_day(datetime_new)
        CALL update_local_time_and_daytime_counter( &
          &     global_seconds_day, dtime, lon, swpar_srf_down, daytime_counter, local_time_day_seconds)

        IF (model%config%flag_slow_sb_pool_spinup_accelerator) THEN
          dsl4jsb_Get_var2d_onDomain(A2L_, slow_sb_pool_accelerator_execute)

          ! in most years the spin-up will not be accelerated
          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
          DO iblk = 1, nblks
            DO ic = 1, nproma
              slow_sb_pool_accelerator_execute(ic, iblk) = 0.0_wp
            ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP

          IF(is_newyear(datetime_new, dtime)) THEN
            current_year = get_year(datetime_new)
            dsl4jsb_Get_var2d_onDomain(A2L_, slow_sb_pool_accelerator_execution_counter)
            accelerator_max_executions = model%config%slow_sb_pool_spinup_accelerator_max_executions
            accelerator_frequency = model%config%slow_sb_pool_spinup_accelerator_frequency
            accelerator_start_year = model%config%slow_sb_pool_spinup_accelerator_start_year
            !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2)
            DO iblk = 1, nblks
              DO ic = 1, nproma
                CALL update_slow_sb_pool_accelerator_bookkeeping( dtime, current_year, &
                  & accelerator_max_executions, accelerator_frequency, accelerator_start_year, &
                  & slow_sb_pool_accelerator_execution_counter(ic, iblk), slow_sb_pool_accelerator_execute(ic, iblk))
              ENDDO
            ENDDO
            !$ACC END PARALLEL LOOP
          END IF
        END IF
      END SELECT
#endif

      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(evapo_act2pot)
      DEALLOCATE(evapo_act2pot)

    END IF ! IF (use_external_forcing) THEN

!$OMP PARALLEL DO PRIVATE(iblk,ics,ice,nc) IF(nblks > 1)
    ! ----------------------------------------------------------------------------------------------------- !
    !> loop over the chunks of the block
    !>
    DO iblk = startblk, endblk
      ! index of gridcell at start of the chunk
      ics = startidx(iblk)
      ! index of gridcell at end of the chunk
      ice = endidx  (iblk)

      IF (ics > ice) CYCLE

      IF (timer_on() .OR. l_timer_host) CALL timer_start(timer_jsbach(model%id))

      ! calc size of chunk, i.e., number of gridcells
      nc = ice - ics + 1

      IF (debug_on('basic') .AND. iblk == 1) CALL message( TRIM(routine), 'Updating '//TRIM(model%name))

      ! Update options
      CALL model%Set_options(iblk=iblk, ics=ics, ice=ice, nc=nc, current_datetime=datetime_old, &
        &                    dtime=dtime, steplen=dtime, alpha=1._wp)
      CALL model%Set_subset(type=ON_CHUNK, iblk=iblk, ics=ics, ice=ice)
      CALL model%Associate_var_pointers(ics, ice, iblk, iblk)

      ! Run tasks in queue
      CALL model%Run_tasks(debug_on('hsm') .AND. iblk == 1)

      ! @TODO cleanup / remove ?
      !!$    IF (iblk == 1)                                                              &
      !!$      & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//          &
      !!$      &                             ' on PE '//int2string(p_pe)//' completed.', &
      !!$      &              all_print=.TRUE.)

      IF (debug_on('basic') .AND. iblk == 1) &
        & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//' completed.')

      IF (timer_on() .OR. l_timer_host) CALL timer_stop(timer_jsbach(model%id))

    END DO
!$OMP END PARALLEL DO

    ! ----------------------------------------------------------------------------------------------------- !
    !> finish timestep
    !>
    CALL jsbach_finish_timestep(model_id, datetime_old, dtime)

    CALL deallocateDatetime(datetime_old)
    CALL deallocateDatetime(datetime_old_old)
    CALL deallocateDatetime(datetime_forcing)

    !$ACC END DATA

  END SUBROUTINE run_one_timestep

#endif
#endif
END MODULE mo_jsbach_model
