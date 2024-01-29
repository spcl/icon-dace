! @page pagecontrolmodelf901 Main program for the ICON hamocc model
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

MODULE mo_hamocc_model

  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_master_config,       ONLY: isRestart
  USE mo_master_control,      ONLY: hamocc_process, get_my_process_name
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs, &
       &                            pio_type, num_test_pe
  USE mo_mpi,                 ONLY: set_mpi_work_communicators, process_mpi_io_size, &
       &                            stop_mpi, my_process_is_io, my_process_is_mpi_test,   &
       &                            process_mpi_io_size
  USE mo_timer,               ONLY: init_timer, timer_start, timer_stop, print_timer, &
       &  timer_model_init, timer_total
  USE mo_memory_log,              ONLY: memory_log_terminate
  USE mtime,                  ONLY: MAX_DATETIME_STR_LEN, datetimeToString
  USE mo_name_list_output_init, ONLY: init_name_list_output, parse_variable_groups, &
    &                                 output_file
  USE mo_name_list_output,    ONLY: close_name_list_output
  USE mo_name_list_output_config,  ONLY: use_async_name_list_io
  USE mo_level_selection, ONLY: create_mipz_level_selections
  USE mo_dynamics_config,     ONLY: configure_dynamics
  USE mo_zaxis_type,          ONLY: zaxisTypeList, t_zaxisTypeList

  !  USE mo_advection_config,    ONLY: configure_advection
 USE mo_gribout_config,      ONLY: configure_gribout

  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_run_config,          ONLY: &
    & test_mode,              &
    & ldynamics,              &
    & ltransport,             &
    & dtime,                  & !    :
    & ltimer,                 & !    :
    & num_lev,                &
    & nshift,                 &
    & grid_generatingcenter,  & 
    & grid_generatingsubcenter, &
    & configure_run, output_mode 

  USE mo_ocean_nml_crosscheck,   ONLY: ocean_crosscheck
  USE mo_ocean_nml,              ONLY: i_sea_ice, no_tracer, use_omip_forcing, lhamocc, &
    & initialize_fromRestart, n_zlev

  USE mo_model_domain,        ONLY: t_patch_3d, p_patch_local_parent

  ! Horizontal grid
  !
  USE mo_grid_config,         ONLY: n_dom, use_dummy_cell_closure

  USE mo_build_decomposition, ONLY: build_decomposition
  USE mo_complete_subdivision,ONLY: setup_phys_patches

  USE mo_ocean_ext_data,      ONLY: ext_data, construct_ocean_ext_data, destruct_ocean_ext_data
  USE mo_bgc_bcond,           ONLY: construct_bgc_ext_data, destruct_bgc_ext_data, ext_data_bgc
  USE mo_ocean_types,           ONLY: t_hydro_ocean_state, &
    & t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_ocean_state,           ONLY:  v_base, &
    & construct_hydro_ocean_base, &! destruct_hydro_ocean_base, &
    & construct_patch_3d, destruct_patch_3d, &
    & construct_ocean_var_lists
    
  USE mo_ocean_initialization, ONLY: init_ho_base, &
    & init_ho_basins, init_coriolis_oce, init_patch_3d,   &
    & init_patch_3d
  USE mo_hamocc_output,        ONLY: construct_hamocc_var_lists, construct_hamocc_state, &
    & destruct_hamocc_state      
  USE mo_ocean_initial_conditions,  ONLY:  init_ocean_bathymetry
  USE mo_ocean_check_tools,     ONLY: init_oce_index
  USE mo_util_dbg_prnt,       ONLY: init_dbg_index, dbg_print
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_operator_ocean_coeff_3d,ONLY: construct_operators_coefficients, &
    & destruct_operators_coefficients

   USE mo_impl_constants,      ONLY: success
  
  USE mo_alloc_patches,        ONLY: destruct_patches, destruct_comm_patterns
  USE mo_ocean_read_namelists, ONLY: read_ocean_namelists
  USE mo_load_restart,         ONLY: read_restart_header, read_restart_files
  USE mo_key_value_store,      ONLY: t_key_value_store
  USE mo_restart_nml_and_att,  ONLY: getAttributesForRestarting, ocean_initFromRestart_OVERRIDE
  USE mo_ocean_patch_setup,    ONLY: complete_ocean_patch
  USE mo_icon_comm_interface,  ONLY: construct_icon_communication, destruct_icon_communication
  USE mo_output_event_types,   ONLY: t_sim_step_info
  USE mo_grid_tools,           ONLY: create_dummy_cell_closure
  USE mo_io_config,            ONLY: restartWritingParameters
  USE mo_bgc_icon_comm,        ONLY: hamocc_state
  USE mo_ocean_time_events,    ONLY: init_ocean_time_events, getCurrentDate_to_String,          &
    & ocean_time_nextStep, isCheckpoint, isEndOfThisRun, newNullDatetime,  &
    & get_OceanCurrentTime_Pointer
  USE mtime,                   ONLY: datetime, datetimeToString, datetimeToPosixString, MAX_DATETIME_STR_LEN, deallocateDatetime

  USE mo_derived_variable_handling, ONLY: update_statistics
  USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
  USE mo_restart,                ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
 
  !-------------------------------------------------------------
  USE mo_construct_icon_hamocc, ONLY: construct_icon_hamocc, destruct_icon_hamocc, init_icon_hamocc
  USE mo_ocean_hamocc_couple_state, ONLY:   t_ocean_to_hamocc_state, t_hamocc_to_ocean_state, &
    & t_ocean_transport_state, t_hamocc_ocean_state, &
    & hamocc_ocean_state, hamocc_ocean_state_list, &
    & construct_hamocc_ocean_state, destruct_hamocc_ocean_state
  USE mo_hamocc_ocean_physics,   ONLY: tracer_biochemistry_transport

  USE mo_ocean_hamocc_interface, ONLY: hamocc_to_ocean_init, hamocc_to_ocean_end, hamocc_to_ocean_interface
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_ocean_math_operators,   ONLY: update_height_hamocc
#ifdef YAC_coupling
  USE mo_io_coupling_frame,      ONLY: construct_io_coupling, destruct_io_coupling
#endif
  USE mo_icon_output_tools,      ONLY: init_io_processes, prepare_output
#ifndef __NO_ICON_COMIN__
  USE mo_mpi,               ONLY: p_comm_comin
  USE comin_host_interface, ONLY: mpi_handshake_dummy
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hamocc_model
  PUBLIC :: patch_3d, ocean_operators_coefficients

  TYPE(t_patch_3d), POINTER                       :: patch_3d
  TYPE(t_operator_coeff), TArGET                  :: ocean_operators_coefficients ! needed for running the transport
  TYPE(t_solverCoeff_singlePrecision), TARGET     :: solverCoefficients_sp
  CLASS(t_RestartDescriptor), POINTER :: restartDescriptor

  CONTAINS


  !--------------------------------------------------------------------------
  !>
  SUBROUTINE hamocc_model(hamocc_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: hamocc_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_hamocc_model:hamocc_model"

    !-------------------------------------------------------------------
    IF (isRestart()) THEN
      CALL read_restart_header("hamocc")
    END IF

    !-------------------------------------------------------------------
    ! initialize dynamic list of vertical axes
    !-------------------------------------------------------------------

    zaxisTypeList = t_zaxisTypeList()

    !-------------------------------------------------------------------
    CALL construct_hamocc_model(hamocc_namelist_filename,shr_namelist_filename)
 
    CALL hamocc_to_ocean_init()
    !-------------------------------------------------------------------
    IF (isRestart()) THEN ! .OR. initialize_fromRestart) THEN
      ocean_initFromRestart_OVERRIDE = initialize_fromRestart
      ! This is an resumed integration. Read model state from restart file(s).
      CALL read_restart_files( patch_3d%p_patch_2d(1) )
      CALL message(method_name,'normal exit from read_restart_files')
      !ELSE
      !  Prepare the initial conditions:
      !  forcing is part of the restart file
    END IF ! isRestart()
    
          

    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------
    CALL prepare_output()

    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------
    CALL hamocc_timestep()
    
    
    !------------------------------------------------------------------
    ! This is the end...
    !------------------------------------------------------------------   
    CALL print_timer()
    CALL hamocc_to_ocean_end()
    
    CALL destruct_hamocc_model()
     
  END SUBROUTINE hamocc_model
  !--------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------
  SUBROUTINE hamocc_timestep()
 
    TYPE(datetime), POINTER        :: current_time     => NULL()
    CHARACTER(LEN=32)              :: datestring
    INTEGER :: jstep,  jstep0
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CHARACTER(*), PARAMETER :: method_name = "mo_hamocc_model:hamocc_timestep"
    
    current_time => newNullDatetime()
    jstep0 = 0
    CALL getAttributesForRestarting(restartAttributes)
    IF (restartAttributes%is_init) THEN
      ! get start counter for time loop from restart file:
      CALL restartAttributes%get("jstep", jstep0)
    END IF
    
    jstep = jstep0
    
    
    CALL timer_start(timer_total)
    CALL output_hamocc(jstep, get_OceanCurrentTime_Pointer())
    
    ! timestep ...
    !-------------------------------------------------------------------------
    TIME_LOOP: DO
      jstep = jstep + 1    
      current_time = ocean_time_nextStep()

      CALL datetimeToString(current_time, datestring)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
      CALL message (method_name, message_text)

      CALL hamocc_to_ocean_interface()
      CALL update_height_hamocc( patch_3D, ocean_operators_coefficients, hamocc_ocean_state%ocean_to_hamocc_state%h_old_withIce)
      CALL tracer_biochemistry_transport(hamocc_ocean_state, ocean_operators_coefficients, current_time)

      !-------------------------------------------------------------------------
      CALL output_hamocc(jstep, current_time)
      
      CALL update_time_indices
      !-------------------------------------------------------------------------
      IF (isCheckpoint()) THEN
        IF (.NOT. output_mode%l_none ) THEN
          call message(method_name, "writeRestart")
          !
          ! For multifile restart (restart_write_mode = "joint procs multifile")
          ! the domain flag must be set to .TRUE. in order to activate the domain,
          ! even though we have one currently in the ocean. Without this the
          ! processes won't write out their data into a the patch restart files.
          patch_3d%p_patch_2d(1)%ldom_active = .TRUE.
          !
          CALL restartDescriptor%updatePatch(patch_3d%p_patch_2d(1), &
                                            &opt_nice_class=1, &
                                            &opt_ocean_zlevels=n_zlev, &
                                            &opt_ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:), &
                                            &opt_ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
          CALL restartDescriptor%writeRestart(current_time, jstep)
        END IF
      END IF
      !-------------------------------------------------------------------------            
      IF (isEndOfThisRun()) THEN
        ! leave time loop
        EXIT TIME_LOOP
      END IF

    ENDDO TIME_LOOP
    !-------------------------------------------------------------------------
    CALL timer_stop(timer_total)
    
    CALL deallocateDatetime(current_time)   


  END SUBROUTINE hamocc_timestep
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE output_hamocc(jstep, this_datetime)
    INTEGER, INTENT(in)      :: jstep
    TYPE(datetime), POINTER  :: this_datetime
  
    CHARACTER(len=32) :: fmtstr
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: datestring
 
    !------------------------------------------------------------------
   !CALL calc_slow_oce_diagnostics( patch_3d       , &
   !  &                             ocean_state(jg), &
   !  &                             surface_fluxes      , &
   !  &                             sea_ice          , &
   !  &                             jstep-jstep0   , &
   !  &                             this_datetime) ! , &
          ! &                             oce_ts)
     !------------------------------------------------------------------
    CALL update_statistics
    !------------------------------------------------------------------
    IF (istime4name_list_output(jstep) ) THEN
      fmtstr = '%Y-%m-%d %H:%M:%S'
      call datetimeToPosixString(this_datetime, datestring, fmtstr)
      WRITE(message_text,'(a,a)') 'Write output at:', TRIM(datestring)
      CALL message ("output_hamocc",message_text)
      
      IF (output_mode%l_nml) CALL write_name_list_output(jstep)
    ENDIF
    !------------------------------------------------------------------
  END SUBROUTINE output_hamocc
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE update_time_indices()

  INTEGER             :: n_temp
    ! Step 7: Swap time indices before output
    !         half time levels of semi-implicit Adams-Bashforth timestepping are
    !         stored in auxiliary arrays g_n and g_nimd of p_diag%aux
    n_temp    = nold(1)
    nold(1)  = nnew(1)
    nnew(1)  = n_temp
  END SUBROUTINE update_time_indices
  !-------------------------------------------------------------------------
 
  !--------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE destruct_hamocc_model()

    CHARACTER(*), PARAMETER :: method_name = "destruct_hamocc_model"

    INTEGER :: error_status

    !------------------------------------------------------------------
    !  cleaning up process
    !------------------------------------------------------------------
    CALL message(method_name,'start to clean up')
    
    CALL restartDescriptor%destruct()

    CALL destruct_hamocc_ocean_state()
    CALL destruct_icon_hamocc()

    ! Destruct communication patterns
    CALL destruct_comm_patterns( patch_3d%p_patch_2d, p_patch_local_parent )

    !The 3D-ocean version of previous calls
    CALL destruct_patches( patch_3d%p_patch_2d )
    CALL destruct_patches( p_patch_local_parent )
    NULLIFY( patch_3d%p_patch_2d )
    CALL destruct_patch_3d( patch_3d )


    ! Delete variable lists

    IF (output_mode%l_nml) CALL close_name_list_output
   
#ifdef YAC_coupling
     CALL destruct_io_coupling ( get_my_process_name() )
#endif

    CALL destruct_icon_communication()

    CALL destruct_operators_coefficients(ocean_operators_coefficients, solverCoefficients_sp)

    ! close memory logging files
    CALL memory_log_terminate

    CALL message(method_name,'clean-up finished')

  END SUBROUTINE destruct_hamocc_model
  !--------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------
  !>
  !!
!<Optimize:inUse>
  SUBROUTINE construct_hamocc_model(hamocc_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: hamocc_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_hamocc_model:construct_hamocc_model"
    INTEGER :: ist, error_status, dedicatedRestartProcs
    !-------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------
    ! for the moment we just read the hamocc nemalis
    CALL read_ocean_namelists(hamocc_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------
    CALL ocean_crosscheck()

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run

    !-------------------------------------------------------------------
    CALL init_ocean_time_events()


    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL restartWritingParameters(opt_dedicatedProcCount = dedicatedRestartProcs)
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs, &
      &  dedicatedRestartProcs, my_comp_id=hamocc_process, num_test_pe=num_test_pe, pio_type=pio_type)

#ifndef __NO_ICON_COMIN__
    ! we dont participate at comin (yet) but we need to be friendly and shake hands
    CALL mpi_handshake_dummy(p_comm_comin)
#endif

#ifdef YAC_coupling
      ! The initialisation of YAC needs to be called by all (!) MPI processes
      ! in MPI_COMM_WORLD.
      ! construct_io_coupling needs to be called before init_name_list_output
      ! due to calling sequence in subroutine atmo_model for other atmosphere
      ! processes
      CALL construct_io_coupling ( get_my_process_name()  )
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
    CALL build_decomposition(num_lev,nshift, is_ocean_decomposition =.TRUE., &
      & patch_3d=patch_3d)
    CALL construct_icon_communication(patch_3d%p_patch_2d(:), n_dom=1)
    CALL complete_ocean_patch(patch_3d%p_patch_2d(1))
    ! Setup the information for the physical patches
    CALL setup_phys_patches

    ! we need the nnow info
    CALL configure_dynamics ( n_dom, ldynamics, ltransport )

    !------------------------------------------------------------------
    CALL construct_ocean_var_lists(patch_3d%p_patch_2d(1))
    !------------------------------------------------------------------
    ! step 5b: allocate state variables
    !------------------------------------------------------------------
!     ALLOCATE (ocean_state(n_dom), stat=ist)
!     IF (ist /= success) THEN
!       CALL finish(method_name,'allocation for ocean_state failed')
!     ENDIF

    !if(lhamocc)then
    !ALLOCATE (hamocc_state(n_dom), stat=ist)
    !IF (ist /= success) THEN
    !  CALL finish(method_name,'allocation for hamocc_state failed')
    !ENDIF
    !ENDIF
    !---------------------------------------------------------------------
    ! 9. Horizontal and vertical grid(s) are now defined.
    !    Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------

    CALL configure_gribout(grid_generatingcenter, grid_generatingsubcenter, n_dom)

    !------------------------------------------------------------------
    ! 10. Create and optionally read external data fields
    !------------------------------------------------------------------
     ALLOCATE (ext_data(1), stat=error_status)
     IF (error_status /= success) THEN
       CALL finish(method_name,'allocation for ext_data failed')
    ENDIF
    
! 
!     ! allocate memory for oceanic external data and
!     ! optionally read those data from netCDF file.
    CALL construct_ocean_ext_data(patch_3d%p_patch_2d(1:), ext_data)
    ! initial analytic bathymetry via namelist
    CALL init_ocean_bathymetry(patch_3d=patch_3d,  &
      & cells_bathymetry=ext_data(1)%oce%bathymetry_c(:,:))    


    CALL construct_patch_3d(patch_3d)

    CALL construct_hydro_ocean_base(patch_3d%p_patch_2d(1), v_base)
    CALL init_ho_base (patch_3d%p_patch_2d(1), ext_data(1), v_base)
    CALL init_ho_basins(patch_3d%p_patch_2d(1), v_base) ! This initializes the wet_c,..., for all cells ! unbelievable !
    CALL init_coriolis_oce(patch_3d%p_patch_2d(1) )
    CALL init_patch_3d    (patch_3d, ext_data(1), v_base)
    CALL construct_operators_coefficients     ( patch_3d, ocean_operators_coefficients, solverCoefficients_sp)

    IF (use_dummy_cell_closure) CALL create_dummy_cell_closure(patch_3d)
    !------------------------------------------------------------------
    
    CALL construct_icon_hamocc(patch_3d, ext_data(1))
    CALL construct_hamocc_ocean_state(patch_3d)   
    
    !---------------------------------------------------------------------
    restartDescriptor => createRestartDescriptor("hamocc")
    !---------------------------------------------------------------------

    IF (ltimer) CALL timer_stop(timer_model_init)

  END SUBROUTINE construct_hamocc_model
  !--------------------------------------------------------------------------

END MODULE mo_hamocc_model

