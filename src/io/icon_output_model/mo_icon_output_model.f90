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
#include "iconfor_dsl_definitions.inc"
!=============================================================================================
MODULE mo_icon_output_model

  USE mo_exception,           ONLY: message, message_text
  USE mo_master_control,      ONLY: get_my_process_name, get_my_process_type
  USE mo_master_config,       ONLY: isRestart
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs, &
       &                            pio_type, num_test_pe, num_prefetch_proc
  USE mo_mpi,                 ONLY: set_mpi_work_communicators
#ifdef HAVE_CDI_PIO
  USE mo_impl_constants,      ONLY: pio_type_cdipio
  USE mo_cdi_pio_interface,   ONLY: nml_io_cdi_pio_namespace
  USE mo_cdi,                 ONLY: namespaceGetActive, namespaceSetActive
#endif
  USE mo_timer,               ONLY: init_timer, timer_start, timer_stop, print_timer, &
       &                            timer_model_init
  USE mo_memory_log,          ONLY: memory_log_terminate
  USE mtime,                  ONLY: datetimeToString, datetime
  USE mo_name_list_output,    ONLY: close_name_list_output, write_name_list_output
  USE mo_zaxis_type,          ONLY: zaxisTypeList, t_zaxisTypeList

  !  USE mo_advection_config,    ONLY: configure_advection
  USE mo_run_config,          ONLY: output_mode !,ldynamics, ltransport
  USE mo_gribout_config,      ONLY: configure_gribout
  ! Control parameters: run control, dynamics, i/o
  USE mo_run_config,          ONLY: &
    & ltimer,                 & !    :
    & num_lev,                &
    & nshift,                 &
    & grid_generatingcenter,  & ! grid generating center
    & grid_generatingsubcenter  ! grid generating subcenter
  ! Horizontal grid
  USE mo_model_domain,        ONLY: p_patch_local_parent
  USE mo_grid_config,         ONLY: n_dom
!  USE mo_grid_config,         ONLY: use_dummy_cell_closure
  USE mo_build_decomposition, ONLY: build_decomposition
  USE mo_complete_subdivision,ONLY: setup_phys_patches
  USE mo_util_dbg_prnt,       ONLY: init_dbg_index
  USE mo_alloc_patches,        ONLY: destruct_comm_patterns
  USE mo_icon_output_read_namelists, ONLY: read_icon_output_namelists
  USE mo_load_restart,         ONLY: read_restart_header
  USE mo_icon_comm_interface,  ONLY: construct_icon_communication, destruct_icon_communication
  USE mo_io_config,            ONLY: restartWritingParameters
!  USE mo_io_config,            ONLY: write_initial_state
  USE mo_icon_output_time_events,   ONLY: init_icon_output_time_events, isEndOfThisRun, &
    & icon_output_time_nextStep, newNullDatetime
  !-------------------------------------------------------------
  USE mo_icon_output_tools,    ONLY: init_io_processes, prepare_output
  ! For the coupling
  USE mo_coupling_config,           ONLY: is_coupled_run
  USE mo_icon_output_coupling,      ONLY: construct_icon_output_coupling, destruct_icon_output_coupling

  !-------------------------------------------------------------
  USE mo_icon_output_variables, ONLY: construct_icon_output_variables, destruct_icon_output_variables, &
       & patch_3d

#ifndef __NO_ICON_COMIN__
  USE mo_mpi,               ONLY: p_comm_comin
  USE comin_host_interface, ONLY: mpi_handshake_dummy
#endif


  IMPLICIT NONE

  PRIVATE
#ifdef HAVE_CDI_PIO
  INCLUDE 'cdipio.inc'
#endif

    PUBLIC :: icon_output_driver

  CONTAINS


  !--------------------------------------------------------------------------
  !>
  SUBROUTINE icon_output_driver(icon_output_namelist_filename,shr_namelist_filename)
    CHARACTER(LEN=*), INTENT(in) :: icon_output_namelist_filename,shr_namelist_filename

    INTEGER :: step
    CHARACTER(*), PARAMETER :: method_name = "mo_icon_output_model:icon_output_driver"

    CALL message("Hi, this the ICON Output Model. My name is", TRIM(get_my_process_name()))
    !-------------------------------------------------------------------
    IF (isRestart()) THEN
      CALL read_restart_header(TRIM(get_my_process_name()) )
    END IF

    !-------------------------------------------------------------------
    ! initialize dynamic list of vertical axes
    !-------------------------------------------------------------------
    zaxisTypeList = t_zaxisTypeList()

    !-------------------------------------------------------------------
    CALL construct_icon_output_model(icon_output_namelist_filename,shr_namelist_filename)    
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
!     IF (isRestart() .OR. initialize_fromRestart) THEN
!       icon_output_initFromRestart_OVERRIDE = initialize_fromRestart
!       ! This is an resumed integration. Read model state from restart file(s).
!       CALL read_restart_files( patch_3d%p_patch_2d(1) )
!       CALL message(TRIM(method_name),'normal exit from read_restart_files')
!     END IF ! isRestart()
    !------------------------------------------------------------------


    CALL prepare_output()

    !------------------------------------------------------------------
    ! write initial state
    !------------------------------------------------------------------    
!     IF (output_mode%l_nml .and. .true.) THEN
!     IF (output_mode%l_nml .AND. write_initial_state) THEN
!       CALL write_initial_icon_output_timestep(patch_3d,icon_output_state(1),v_oce_sfc,v_sea_ice, operators_coefficients)
!     ENDIF
    !------------------------------------------------------------------
    step=0
    CALL timestep_icon_output(output_step=step)
    !------------------------------------------------------------------

    CALL print_timer()

    !------------------------------------------------------------------
    !  cleaning up process

    CALL destruct_icon_output_model()

  END SUBROUTINE icon_output_driver
  !--------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------
  SUBROUTINE timestep_icon_output(output_step)
    INTEGER :: output_step
    
    TYPE(datetime), POINTER  :: current_time => NULL()
    CHARACTER(LEN=32)        :: datestring
    CHARACTER(*), PARAMETER  :: method_name = "timestep_icon_output"
  
    CALL message(method_name, "starts...")
    current_time => newNullDatetime()
  
    DO WHILE(.not. isEndOfThisRun())
      output_step = output_step + 1
      ! update model date and time mtime based
      current_time => icon_output_time_nextStep()

      CALL datetimeToString(current_time, datestring)
      WRITE(message_text,'(2a,i10,2a)') TRIM(get_my_process_name()), ' begin of timestep =',output_step,'  datetime:  ', datestring
      CALL message (TRIM(method_name), message_text)
      
      IF (output_mode%l_nml) CALL write_name_list_output(output_step)
 
    ENDDO
    
    
    CALL message(method_name, "ended")
  
  END SUBROUTINE timestep_icon_output
  !--------------------------------------------------------------------------
 
  !--------------------------------------------------------------------------
  !>
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_icon_output_model()
    CHARACTER(*), PARAMETER :: method_name = "mo_icon_output_model:destruct_icon_output_model"
#ifdef HAVE_CDI_PIO
    INTEGER :: prev_cdi_namespace
#endif

    !  cleaning up process
    CALL message(TRIM(method_name),'start to clean up')
    CALL destruct_comm_patterns( patch_3d%p_patch_2d, p_patch_local_parent )
    CALL destruct_icon_output_variables()
    IF (output_mode%l_nml) CALL close_name_list_output
#ifdef HAVE_CDI_PIO
    IF (pio_type == pio_type_cdipio) THEN
      prev_cdi_namespace = namespaceGetActive()
      CALL namespaceSetActive(nml_io_cdi_pio_namespace)
      CALL pioFinalize
      CALL namespaceSetActive(prev_cdi_namespace)
    END IF
#endif
    CALL destruct_icon_communication()
    IF ( is_coupled_run() ) THEN
      CALL destruct_icon_output_coupling ()
    END IF

    ! close memory logging files
    CALL memory_log_terminate
    CALL message(TRIM(method_name),'clean-up finished')
  END SUBROUTINE destruct_icon_output_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !!
  !! It does not include the restart processes, these are called from the calling method_name icon_output_driver
  !!
!<Optimize:inUse>
  SUBROUTINE construct_icon_output_model(icon_output_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: icon_output_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_icon_output_model:construct_icon_output_model"
    INTEGER :: dedicatedRestartProcs, comp_id
!    INTEGER :: ist
    !-------------------------------------------------------------------

    CALL message(method_name, "starts...")
    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------
    CALL read_icon_output_namelists(icon_output_namelist_filename,shr_namelist_filename)

    !-------------------------------------------------------------------
    CALL init_icon_output_time_events()

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL restartWritingParameters(opt_dedicatedProcCount = dedicatedRestartProcs)
!orig
!    write(0,*)'construct_icon_output_model:pio_type=',pio_type
!    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs, &
!      &                             dedicatedRestartProcs, num_test_pe, pio_type)
!orig
!pa
!pa    
!pa    write(0,*)'construct_icon_output_model:pio_type=',pio_type
!pa    write(0,*)'construct_icon_output_model:restartProcs=',dedicatedRestartProcs
    comp_id = get_my_process_type()
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, &
         &                          num_io_procs, dedicatedRestartProcs, &
         &                          comp_id,num_prefetch_proc, num_test_pe,      &
         &                          pio_type)
    !pa

#ifndef __NO_ICON_COMIN__
    ! we dont participate at comin (yet) but we need to be friendly and shake hands
    CALL mpi_handshake_dummy(p_comm_comin)
#endif

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    CALL init_timer

    IF (ltimer) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! 4. Setup IO procs
    !-------------------------------------------------------------------
    ! If we belong to the I/O PEs just call xxx_io_main_proc before
    ! reading patches.  This routine will never return
    CALL init_io_processes()

    ! 4. Import patches
    !-------------------------------------------------------------------
    ! FixMe: is_ocean_decomposition =.TRUE. should be parametrized
    CALL build_decomposition(num_lev,nshift, is_ocean_decomposition =.TRUE., &
      & patch_3d=patch_3d)
    CALL construct_icon_communication(patch_3d%p_patch_2d(:), n_dom=1)
    ! Setup the information for the physical patches
    CALL setup_phys_patches

    ! we need the nnow info
!     CALL configure_dynamics ( n_dom, ldynamics, ltransport )

    CALL construct_icon_output_variables()

    IF ( is_coupled_run() ) THEN
      CALL construct_icon_output_coupling()
    END IF

    !------------------------------------------------------------------
    ! step 5b: allocate state variables
    !---------------------------------------------------------------------
    ! 9. Horizontal and vertical grid(s) are now defined.
    !    Assign values to derived variables in the configuration states

    CALL configure_gribout(grid_generatingcenter, grid_generatingsubcenter, n_dom)


    CALL init_dbg_index(patch_3d%p_patch_2d(1))!(patch_2D(1))
    !---------------------------------------------------------------------    
!     IF (use_dummy_cell_closure) CALL create_dummy_cell_closure(patch_3d)

    ! initialize icon_output indices for debug output (including 3-dim lsm)
!     CALL init_oce_index(patch_3d%p_patch_2d,patch_3d, icon_output_state, ext_data )


    IF (ltimer) CALL timer_stop(timer_model_init)
    CALL message(method_name, "ended")

  END SUBROUTINE construct_icon_output_model
  !--------------------------------------------------------------------------


END MODULE mo_icon_output_model

