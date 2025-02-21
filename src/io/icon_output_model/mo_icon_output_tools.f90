! @page pagecontrolmodelf901 Main program for the ICON icon_output model
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

MODULE mo_icon_output_tools

  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,     ONLY: pio_type
  USE mo_mpi,                 ONLY: process_mpi_io_size, stop_mpi, my_process_is_io
  USE mo_impl_constants,      ONLY: pio_type_async, pio_type_cdipio
  USE mo_timer,               ONLY: timer_stop, timer_model_init, timers_level
  USE mo_name_list_output,    ONLY: name_list_io_main_proc
  USE mo_name_list_output_init, ONLY: init_name_list_output, parse_variable_groups, &
    &                                 create_vertical_axes, output_file
  USE mo_name_list_output_config,  ONLY: use_async_name_list_io
  USE mo_level_selection, ONLY: create_mipz_level_selections
  USE mo_run_config,          ONLY: output_mode, dtime
  USE mo_restart_nml_and_att,  ONLY: getAttributesForRestarting
  USE mo_output_event_types,   ONLY: t_sim_step_info
  USE mo_key_value_store,     ONLY: t_key_value_store
  USE mo_time_config,         ONLY: time_config

  !-------------------------------------------------------------
 
  IMPLICIT NONE

  PRIVATE

    PUBLIC :: init_io_processes
    PUBLIC :: prepare_output

  CONTAINS
  !-------------------------------------------------------------------------
  SUBROUTINE init_io_processes()
    TYPE(t_sim_step_info)   :: sim_step_info
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CHARACTER(*), PARAMETER :: routine = 'mo_icon_output_tools:init_io_processes'

    IF (process_mpi_io_size > 0 .AND. pio_type == pio_type_async) THEN
      ! Decide whether async vlist or name_list IO is to be used,
      ! only one of both may be enabled!
      IF (output_mode%l_nml) THEN
        ! -----------------------------------------
        ! asynchronous I/O
        ! -----------------------------------------
        !
        use_async_name_list_io = .TRUE.
        CALL message(routine,'asynchronous namelist I/O scheme is enabled.')
        ! consistency check
        IF (my_process_is_io()) THEN
          ! Stop timer which is already started but would not be stopped
          ! since xxx_io_main_proc never returns
          IF (timers_level > 1) CALL timer_stop(timer_model_init)

          ! compute sim_start, sim_end
          sim_step_info%sim_start = time_config%tc_exp_startdate
          sim_step_info%sim_end = time_config%tc_exp_stopdate
          sim_step_info%run_start = time_config%tc_startdate
          sim_step_info%restart_time = time_config%tc_stopdate
          sim_step_info%dtime = dtime

          CALL getAttributesForRestarting(restartAttributes)
          IF (restartAttributes%is_init) THEN

            ! get start counter for time loop from restart file:
            CALL restartAttributes%get("jstep", sim_step_info%jstep0)
          ELSE
            sim_step_info%jstep0 = 0
          END IF
          ! If we belong to the I/O PEs just call xxx_io_main_proc before
          ! reading patches.  This routine will never return
          CALL name_list_io_main_proc(sim_step_info)
        END IF
      ELSE IF (my_process_is_io()) THEN
        ! Shut down MPI
        CALL stop_mpi
        STOP
      ENDIF
    ELSE IF (process_mpi_io_size > 0 .AND. pio_type == pio_type_cdipio) THEN

      CALL finish(routine, 'CDI-PIO requested but unavailable')
    ELSE
      ! -----------------------------------------
      ! non-asynchronous I/O (performed by PE #0)
      ! -----------------------------------------
      !
      IF (output_mode%l_nml) THEN
        CALL message(routine, 'synchronous namelist I/O scheme is enabled.')
      ENDIF
      IF (my_process_is_io()) THEN
        ! Shut down MPI
        CALL stop_mpi
        STOP
      ENDIF
    END IF
  END SUBROUTINE init_io_processes
  !-------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE prepare_output()
    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:prepare_output"
    TYPE(t_sim_step_info)               :: sim_step_info
    TYPE(t_key_value_store), POINTER :: restartAttributes

    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------

    IF (output_mode%l_nml) THEN
!       WRITE(0,*)'process_mpi_io_size:',process_mpi_io_size
!       IF (process_mpi_io_size > 0) use_async_name_list_io = .TRUE.
      CALL parse_variable_groups()
      ! compute sim_start, sim_end
      sim_step_info%sim_start = time_config%tc_exp_startdate
      sim_step_info%sim_end = time_config%tc_exp_stopdate
      sim_step_info%run_start = time_config%tc_startdate
      sim_step_info%restart_time = time_config%tc_stopdate

      sim_step_info%dtime      = dtime
      sim_step_info%jstep0 = 0

      CALL getAttributesForRestarting(restartAttributes)
      ! get start counter for time loop from restart file:
      IF (restartAttributes%is_init) &
        & CALL restartAttributes%get("jstep", sim_step_info%jstep0)
      CALL init_name_list_output(sim_step_info, opt_lprintlist=.TRUE.,opt_l_is_ocean=.TRUE.)
      CALL create_mipz_level_selections(output_file)
      CALL create_vertical_axes(output_file)
    ENDIF
  
  END SUBROUTINE prepare_output
  !--------------------------------------------------------------------------

END MODULE mo_icon_output_tools

