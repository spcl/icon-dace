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

! @brief Main program for the ICON atmospheric model

MODULE mo_atmo_model

  ! basic modules
  USE mo_exception,               ONLY: message, finish
  USE mo_mpi,                     ONLY: set_mpi_work_communicators,       &
    &                                   my_process_is_pref, process_mpi_pref_size,   &
    &                                   my_process_is_work,      &
    &                                   my_process_is_mpi_test
  USE mo_timer,                   ONLY: init_timer, timer_start, timer_stop,                  &
    &                                   timers_level, timer_model_init, timer_coupling,       &
    &                                   timer_domain_decomp, timer_compute_coeffs,            &
    &                                   timer_ext_data, print_timer
#ifdef HAVE_RADARFWO
  USE mo_mpi, ONLY: my_process_is_mpi_test, my_process_is_radario,        &
                    get_my_mpi_all_id, process_mpi_radario_size,                              &
                    MPI_COMM_NULL, MPI_UNDEFINED, p_comm_work, get_my_mpi_work_id,            &
                    process_mpi_all_comm, num_work_procs, process_mpi_all_workroot_id,        &
                    process_mpi_all_radarioroot_id, process_mpi_all_size
  USE mo_emvorado_init,           ONLY: prep_emvorado_domains, init_emvorado_mpi
  USE mo_emvorado_interface,      ONLY: radar_mpi_barrier
#ifndef NOMPI
  USE mo_emvorado_interface,      ONLY: exchg_with_detached_emvorado_io,                      &
       &                                detach_emvorado_io
#endif
#endif
  USE mo_parallel_config,         ONLY: p_test_run, num_test_pe, l_test_openmp, num_io_procs, &
    &                                   proc0_shift, num_prefetch_proc, pio_type, num_io_procs_radar, &
    &                                   ignore_nproma_use_nblocks_c, ignore_nproma_use_nblocks_e,     &
    &                                   nproma, update_nproma_for_io_procs
  USE mo_master_config,           ONLY: isRestart
  USE mo_memory_log,              ONLY: memory_log_terminate
#ifndef NOMPI
#if defined(__GET_MAXRSS__)
  USE mo_mpi,                     ONLY: get_my_mpi_all_id
  USE mo_util_sysinfo,            ONLY: util_get_maxrss
#endif
#endif
  USE mo_coupling_config,         ONLY: is_coupled_run
  USE mo_coupling_utils,          ONLY: cpl_construct, cpl_destruct
  USE mo_impl_constants,          ONLY: SUCCESS, inwp, LSS_JSBACH, &
    &                                   min_rlcell_int, min_rlcell, max_rlcell, &
    &                                   min_rlvert_int, min_rlvert, max_rlvert, &
    &                                   min_rledge_int, min_rledge, max_rledge
  USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_zaxis_type,              ONLY: zaxisTypeList, t_zaxisTypeList
  USE mo_load_restart,            ONLY: read_restart_header

  ! namelist handling; control parameters: run control, dynamics
  USE mo_read_namelists,          ONLY: read_atmo_namelists
  USE mo_nml_crosscheck,          ONLY: atm_crosscheck
  USE mo_initicon_config,         ONLY: configure_initicon
  USE mo_io_config,               ONLY: restartWritingParameters
  USE mo_lnd_nwp_config,          ONLY: configure_lnd_nwp
  USE mo_dynamics_config,         ONLY: configure_dynamics, lmoist_thdyn
  USE mo_run_config,              ONLY: configure_run,                                        &
    &                                   ltimer, ltestcase,                                    &
    &                                   ldynamics, ltransport,                                &
    &                                   nshift,                                               &
    &                                   num_lev,                                              &
    &                                   msg_level,                                            &
    &                                   grid_generatingCenter,                                & ! grid generating center
    &                                   grid_generatingSubcenter,                             & ! grid generating subcenter
    &                                   iforcing, luse_radarfwo,                              &
    &                                   iqc, iqt, iqv, iqi, iqs, iqr, ltimer,                 &
    &                                   iqni, iqg, iqm_max, iqtke, iqh, iqnr, iqns, iqng,     &
    &                                   iqnh, iqnc, iqgl, iqhl, inccn, ininact, ininpot,      &
    &                                   lart, nqtendphy, ntracer,                             &
    &                                   iqbin, iqb_i, iqb_e, iqb_s
  USE mo_gribout_config,          ONLY: configure_gribout
  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config
#ifndef __NO_JSBACH__
  USE mo_aes_phy_config,          ONLY: aes_phy_config
  USE mo_master_control,          ONLY: master_namelist_filename
  USE mo_jsb_base,                ONLY: jsbach_setup => jsbach_setup_models, jsbach_setup_tiles
  USE mo_jsb_model_init,          ONLY: jsbach_setup_grid
  USE mo_jsb_model_final,         ONLY: jsbach_finalize
#endif
  USE mo_master_control,          ONLY: atmo_process, get_my_process_name

  ! time stepping
  USE mo_atmo_nonhydrostatic,     ONLY: atmo_nonhydrostatic

  USE mo_nh_testcases,            ONLY: init_nh_testtopo

  USE mo_alloc_patches,           ONLY: destruct_patches, destruct_comm_patterns

  ! advection
  USE mo_advection_config,        ONLY: advection_config
  USE mo_advection_utils,         ONLY: init_tracer_settings

  ! horizontal grid, domain decomposition, memory
  USE mo_grid_config,             ONLY: n_dom, n_dom_start, n_phys_dom, l_scm_mode
  USE mo_model_domain,            ONLY: p_patch, p_patch_local_parent
  USE mo_build_decomposition,     ONLY: build_decomposition
  USE mo_complete_subdivision,    ONLY: setup_phys_patches
  USE mo_icon_comm_interface,     ONLY: construct_icon_communication,                         &
    &                                   destruct_icon_communication
  ! Vertical grid
  USE mo_vertical_coord_table,    ONLY: vct_a, vct_b, vct, allocate_vct_atmo
  USE mo_init_vgrid,              ONLY: nflatlev
  USE mo_util_vgrid,              ONLY: construct_vertical_grid

  ! external data, physics
  USE mo_ext_data_state,          ONLY: ext_data, construct_ext_data, destruct_ext_data
  USE mo_ext_data_init,           ONLY: init_ext_data

  USE mo_diffusion_config,        ONLY: configure_diffusion

  ! horizontal interpolation
  USE mo_interpol_config,         ONLY: configure_interpolation
  USE mo_gridref_config,          ONLY: configure_gridref
  USE mo_intp_state,              ONLY: construct_2d_interpol_state,                          &
    &                                   destruct_2d_interpol_state, transfer_interpol_state
  USE mo_grf_intp_state,          ONLY: construct_2d_gridref_state,                           &
    &                                   destruct_2d_gridref_state, transfer_grf_state,        &
    &                                   create_grf_index_lists, destruct_interpol_patterns
  USE mo_intp_data_strc,          ONLY: p_int_state, p_int_state_local_parent
  USE mo_intp_lonlat_types,       ONLY: lonlat_grids
  USE mo_grf_intp_data_strc,      ONLY: p_grf_state, p_grf_state_local_parent
  USE mo_intp_lonlat,             ONLY: compute_lonlat_intp_coeffs

  ! I/O
  USE mo_restart,                 ONLY: detachRestartProcs
  USE mo_icon_output_tools,       ONLY: init_io_processes
  USE mtime,                      ONLY: OPERATOR(<), OPERATOR(+)
#ifndef NOMPI
  ! Prefetching
  USE mo_async_latbc,             ONLY: prefetch_main_proc
#endif
  ! ART
  USE mo_art_config,              ONLY: ctracer_art
#ifdef __ICON_ART
  USE mo_art_config,              ONLY: art_config
  USE mo_art_init_interface,      ONLY: art_init_interface, &
    &                                   art_calc_ntracer_and_names
#endif
  USE mo_sync,                    ONLY: global_max
  USE mo_check_ext_constants,     ONLY: check_ext_constants

#ifndef __NO_ICON_COMIN__
  USE comin_host_interface,       ONLY: comin_parallel_mpi_handshake,     &
    &                                   mpi_handshake_dummy
  USE mo_comin_config,            ONLY: comin_config
  USE mo_comin_adapter,           ONLY: icon_expose_descrdata
  USE mo_time_config,             ONLY: time_config
  USE mo_run_config,              ONLY: number_of_grid_used
  USE mo_util_vgrid_types,        ONLY: vgrid_buffer
  USE mo_run_config,              ONLY: dtime
  USE mo_grid_config,             ONLY: start_time, end_time
  USE mo_mpi,                     ONLY: p_comm_comin
  USE mo_impl_constants,          ONLY: max_dom
#endif

  !-------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE
#ifdef HAVE_CDI_PIO
  INCLUDE 'cdipio.inc'
#endif

  PUBLIC :: atmo_model, construct_atmo_model, destruct_atmo_model

CONTAINS


  !-------------------------------------------------------------------
  !>
  SUBROUTINE atmo_model(atm_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: atm_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:atmo_model"

#ifndef NOMPI
#if defined(__SX__)
    INTEGER  :: maxrss
#endif
#endif

    !---------------------------------------------------------------------
    ! construct the atmo model
    CALL construct_atmo_model(atm_namelist_filename,shr_namelist_filename)


    !---------------------------------------------------------------------
    ! constructs the nonhydrostatic atmospheric model
    ! constructs the coupler
    ! performs integration
    ! destructs the coupler
    ! destructs the nonhydrostatic atmospheric model
    !---------------------------------------------------------------------
    CALL atmo_nonhydrostatic ()


    ! print performance timers:
    IF (ltimer) CALL print_timer

#ifdef HAVE_RADARFWO
#ifndef __SCT__
    IF (process_mpi_radario_size > 0) THEN
      IF (ltimer) THEN
#ifndef NOMPI
        ! To synchronize the EMVORADO async. IO timing output, so that this
        !  timing printout appears after the "normal" timing printout for the workers
        CALL radar_mpi_barrier ()
#endif
      END IF
    END IF
#endif
#endif

    !---------------------------------------------------------------------
    ! 13. Integration finished. Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_model ()

    !---------------------------------------------------------------------
    ! (optional:) write resident set size from OS
#ifndef NOMPI
#if defined(__GET_MAXRSS__)
    IF (msg_level >= 16) THEN
      CALL util_get_maxrss(maxrss)
      PRINT  *, "PE #", get_my_mpi_all_id(), &
        &    ": MAXRSS (MiB) = ", maxrss
    END IF
#endif
#endif

  END SUBROUTINE atmo_model
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  !>
  SUBROUTINE construct_atmo_model(atm_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: atm_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename
    ! local variables
    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:construct_atmo_model"
    INTEGER                 :: jg, jgp, error_status, dedicatedRestartProcs
    CHARACTER(len=1000)     :: message_text = ''
    INTEGER                 :: icomm_cart, my_cart_id, nproma_max, iart_ntracer

    ! initialize global registry of lon-lat grids
    CALL lonlat_grids%init()

    !---------------------------------------------------------------------
    ! 0. If this is a resumed or warm-start run...
    !---------------------------------------------------------------------

    IF (isRestart()) THEN
      CALL message('','Read restart file meta data ...')
      CALL read_restart_header(get_my_process_name())
    ENDIF

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL read_atmo_namelists(atm_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------

    CALL atm_crosscheck

#ifdef MESSY
    CALL messy_setup
#endif

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some components of the state, e.g., num_lev, may be
    !    modified in this subroutine which affects the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run( )


    ! complete initicon config-state
    CALL configure_initicon()

    ! check whether the external impl_constants are still the same as their original values
    CALL check_ext_constants()
    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL restartWritingParameters(opt_dedicatedProcCount = dedicatedRestartProcs)

#ifdef HAVE_RADARFWO
    IF (iforcing == inwp .AND. ANY(luse_radarfwo(1:n_dom))) THEN
      CALL prep_emvorado_domains (n_dom, luse_radarfwo(1:n_dom))
    ENDIF
#endif

    CALL set_mpi_work_communicators(p_test_run, l_test_openmp,                    &
         &                          num_io_procs, dedicatedRestartProcs,          &
         &                          atmo_process, num_prefetch_proc, num_test_pe, &
         &                          pio_type,                                     &
         &                          num_io_procs_radar=num_io_procs_radar,        &
         &                          radar_flag_doms_model=luse_radarfwo(1:n_dom), &
         &                          num_dio_procs=proc0_shift)

#ifndef __NO_ICON_COMIN__
    ! Non-work PEs dont participate in the plugin comms
    IF (my_process_is_work()) THEN
      CALL comin_parallel_mpi_handshake(p_comm_comin, &
           & comin_config%plugin_list(1:comin_config%nplugins)%comm, TRIM(get_my_process_name()))
    ELSE
      CALL mpi_handshake_dummy(p_comm_comin)
    ENDIF
#endif


#ifdef HAVE_RADARFWO
!! EMVORADO MPI initialization is also needed in case of output of grid point reflectivities using
!! Mie- or Tmatrix-scattering from EMVORADO.
!! Because in case of luse_radarfwo(:) = .FALSE. this initialization does only very few things, we call it
!! in any case for the NWP ICON:
    IF (iforcing == inwp) THEN
      message_text(:) = ' '
#ifdef NOMPI
      CALL init_emvorado_mpi ( luse_radarfwo(1:n_dom), & ! INPUT
                               MPI_COMM_NULL, 0, 1,    & ! INPUT
                               MPI_COMM_NULL, 0, 1,    & ! INPUT
                               .TRUE.,                 & ! INPUT
                               0, 0, 0,                & ! INPUT
                               error_status, message_text )
#else
      IF (my_process_is_work()) THEN
        icomm_cart = p_comm_work
        my_cart_id = get_my_mpi_work_id()
      ELSE
        icomm_cart = MPI_COMM_NULL
        my_cart_id = MPI_UNDEFINED
      END IF
      CALL init_emvorado_mpi ( luse_radarfwo(1:n_dom),                        & ! INPUT
           process_mpi_all_comm, get_my_mpi_all_id(), process_mpi_all_size,   & ! INPUT
           icomm_cart, my_cart_id, num_work_procs,                            & ! INPUT
           my_process_is_work(),                                              & ! INPUT
           process_mpi_radario_size, process_mpi_all_workroot_id,             & ! INPUT
           process_mpi_all_radarioroot_id,                                    & ! INPUT
           error_status, message_text)
#endif
      IF (error_status /= 0) THEN
        CALL finish ('init_emvorado_mpi', TRIM(message_text))
      END IF
    ENDIF
#endif

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer
    IF (timers_level > 1) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! 3.3 construct basic coupler
    !-------------------------------------------------------------------

    IF (is_coupled_run()) THEN
      IF (ltimer) CALL timer_start(timer_coupling)
      CALL cpl_construct()
      IF (ltimer) CALL timer_stop(timer_coupling)
    END IF

    !-------------------------------------------------------------------
    ! initialize dynamic list of vertical axes
    !-------------------------------------------------------------------

    zaxisTypeList = t_zaxisTypeList()

#ifndef __NO_JSBACH__
    ! Setup JSBACH: read namelists, configure models for each domain
    ! This has to be after (!) the ICON zaxes have been created in the above line but
    ! before (!) the restart PEs are detached a few lines below since JSBACH
    ! adds its zaxes to zaxisTypeList
    IF (ANY(aes_phy_config(:)%ljsb) &
        & .OR. ANY(atm_phy_nwp_config(1:n_dom)%inwp_surface == LSS_JSBACH)) THEN
      ! Do basic initialization of JSBACH
      CALL jsbach_setup(master_namelist_filename)
    END IF
#endif

    !------------------
    ! Next, define the horizontal and vertical grids since they are already
    ! needed for some derived control parameters. This includes
    ! - patch import
    ! - domain decompistion
    ! - vertical coordinates
    !-------------------------------------------------------------------
    ! 4. Import patches, perform domain decomposition
    !-------------------------------------------------------------------

    IF (timers_level > 4) CALL timer_start(timer_domain_decomp)
    ! Only do the decomposition for relevant processes
    IF (my_process_is_work() .OR. my_process_is_mpi_test()) THEN
      CALL build_decomposition(num_lev, nshift, is_ocean_decomposition = .FALSE.)
    ENDIF

    IF (ignore_nproma_use_nblocks_c .OR. ignore_nproma_use_nblocks_e) THEN
      nproma_max = global_max(nproma)
      CALL update_nproma_for_io_procs(nproma_max)
    ENDIF

    IF (timers_level > 4) CALL timer_stop(timer_domain_decomp)

    !-------------------------------------------------------------------
    ! 5. I/O initialization
    !-------------------------------------------------------------------

    ! This won't RETURN on dedicated restart PEs, starting their main loop instead.
    CALL detachRestartProcs(timers_level > 1)

#ifndef NOMPI
    ! If we belong to the prefetching PEs just call prefetch_main_proc before reading patches.
    ! This routine will never return
    IF (process_mpi_pref_size > 0) THEN
      num_prefetch_proc = 1
      CALL message(routine, 'asynchronous input prefetching is enabled.')
      IF (my_process_is_pref()) CALL prefetch_main_proc
    ENDIF
#endif


    CALL init_io_processes()

#ifdef HAVE_RADARFWO
#ifndef NOMPI
    IF ( .NOT. my_process_is_mpi_test() .AND. iforcing == inwp .AND. &
         ANY(luse_radarfwo(1:n_dom)) .AND. num_io_procs_radar > 0   ) THEN

      ! -------------------------------------------------------------------------
      ! Radar forward operator EMVORADO
      ! -------------------------------------------------------------------------
      !
      ! - Initialize data here which must be available on
      !   all radar PEs (workers + radario).
      !
      ! - also detach the separate radario-procs after this initialisation here.
      !
      ! -------------------------------------------------------------------------

      IF (my_process_is_radario()) THEN

        ! We have to call configure_gribout(...) also on the separate radar-IO-procs
        !  due to composite-output. For this, grid_generatingCenter and grid_generatingSubcenter
        !  have to be received on the radar-IO procs from the workers first. Both things are done here.
        ! The counterpart on the worker nodes is called after model initialization below,
        !  because only then the grid_generatingCenter and grid_generatingSubcenter are known:
        CALL exchg_with_detached_emvorado_io (grid_generatingCenter, grid_generatingSubcenter)

        CALL configure_gribout(grid_generatingCenter, grid_generatingSubcenter, n_dom)

        ! This is a pure radar IO PE, so the below "detach_emvorado_io" will never return.
        ! So we STOP the already started timer "timer_model_init", because it is useless on this PE:
        IF (timers_level > 1) CALL timer_stop(timer_model_init)

        CALL detach_emvorado_io (n_dom, luse_radarfwo(1:n_dom))

      END IF
    END IF
#endif
#endif

    !--------------------------------------------------------------------------------
    ! 6. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------

    IF (timers_level > 4) CALL timer_start(timer_compute_coeffs)
    CALL configure_interpolation( n_dom, p_patch(1:)%level, &
                                  p_patch(1:)%geometry_info )

    CALL configure_gridref(n_dom, p_patch(1:)%geometry_info%mean_characteristic_length)


    ! Allocate array for interpolation state

    ALLOCATE( p_int_state(n_dom_start:n_dom), &
            & p_grf_state(n_dom_start:n_dom), &
            & STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(routine, 'allocation for ptr_int_state failed')
    ENDIF

    ALLOCATE( p_int_state_local_parent(n_dom_start+1:n_dom), &
         &    p_grf_state_local_parent(n_dom_start+1:n_dom), &
         &    STAT=error_status)
    IF (error_status /= SUCCESS) &
         CALL finish(routine, 'allocation for local parents failed')

    ! Construct interpolation state
    ! Please note that for parallel runs the divided state is constructed here
    CALL construct_2d_interpol_state(p_patch, p_int_state)

    ! Transfer interpolation state to local parent
    DO jg = n_dom_start+1, n_dom
      jgp = p_patch(jg)%parent_id
      CALL transfer_interpol_state(p_patch(jgp),p_patch_local_parent(jg), &
           &  p_int_state(jgp), p_int_state_local_parent(jg))
    ENDDO

    !-----------------------------------------------------------------------------
    ! 7. Construct grid refinment state, compute coefficients
    !-----------------------------------------------------------------------------
    ! For the NH model, the initialization routines called from
    ! construct_2d_gridref_state require the metric terms to be present

    IF (n_dom_start==0 .OR. n_dom > 1) THEN

      ! Construct gridref state
      ! For non-parallel runs (either 1 PE only or on the test PE) this is done as usual,
      ! for parallel runs, the main part of the gridref state is constructed on the
      ! local parent with the following call
      CALL construct_2d_gridref_state (p_patch(n_dom_start:),                  &
        &                              p_patch_local_parent(n_dom_start+1:),   &
        &                              p_grf_state(n_dom_start:),              &
        &                              p_grf_state_local_parent(n_dom_start+1:))

      ! Transfer gridref state from local parent to p_grf_state
      DO jg = n_dom_start+1, n_dom
        jgp = p_patch(jg)%parent_id
        CALL transfer_grf_state(p_patch(jgp), p_patch_local_parent(jg),         &
          &                     p_grf_state(jgp), p_grf_state_local_parent(jg), &
          &                     p_patch(jg)%parent_child_index)
      ENDDO
    ENDIF


    !-------------------------------------------------------------------
    ! Initialize icon_comm_lib
    !-------------------------------------------------------------------
    !    IF (use_icon_comm) THEN
    CALL construct_icon_communication(p_patch, n_dom)
    !    ENDIF


    !--------------------------------------------
    ! Setup the information for the physical patches
    CALL setup_phys_patches

    !-------------------------------------------------------------------
    ! 8. Constructing data for lon-lat interpolation
    !-------------------------------------------------------------------

    CALL compute_lonlat_intp_coeffs(p_patch(1:), p_int_state(1:))

    IF (n_dom_start==0 .OR. n_dom > 1) THEN
      CALL create_grf_index_lists(p_patch, p_grf_state, p_int_state)
    ENDIF
    IF (timers_level > 4) CALL timer_stop(timer_compute_coeffs)


    !---------------------------------------------------------------------
    ! Prepare dynamics and land
    !---------------------------------------------------------------------

    CALL configure_dynamics ( n_dom, ldynamics, ltransport )

    IF (iforcing == inwp) THEN ! set dimensions of tile-based variables
      CALL configure_lnd_nwp ( &
          & ANY(atm_phy_nwp_config(1:n_dom)%inwp_surface == LSS_JSBACH) &
        )
    ENDIF

    !------------------------------------------------------------------
    ! Create and optionally read external data fields
    !------------------------------------------------------------------

    ! allocate memory for atmospheric/oceanic external data and
    ! optionally read those data from netCDF file.
    IF (timers_level > 4) CALL timer_start(timer_ext_data)
    CALL construct_ext_data (p_patch(1:), ext_data)
    CALL init_ext_data (p_patch(1:), p_int_state(1:), ext_data)
    IF (timers_level > 4) CALL timer_stop(timer_ext_data)

    !---------------------------------------------------------------------
    ! Import vertical grid/ define vertical coordinate
    !---------------------------------------------------------------------

    CALL allocate_vct_atmo(p_patch(1)%nlevp1)
    IF (ltestcase .AND. (.NOT. l_scm_mode)) THEN
      CALL init_nh_testtopo(p_patch(1:), ext_data)   ! set analytic topography
      ! for single column model (SCM) the topography is read in ext_data_init from SCM input file
    ENDIF
    CALL construct_vertical_grid(p_patch(1:), p_int_state(1:), ext_data, &
      &                          vct_a, vct_b, vct, nflatlev)


    !---------------------------------------------------------------------
    ! Horizontal and vertical grid(s) are now defined.
    ! Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------

    CALL configure_diffusion( n_dom, p_patch(1:)%parent_id )

    CALL configure_gribout(grid_generatingCenter, grid_generatingSubcenter, n_phys_dom)



    !-------------------------------------------------------------------------
    ! EMVORADO: The worker part of the communication with radar IO PEs
    ! To send informations on grid_generatingCenter, grid_generatingSubcenter
    !-------------------------------------------------------------------------

#ifdef HAVE_RADARFWO
#ifndef NOMPI
    IF ( .NOT. my_process_is_mpi_test() .AND. iforcing == inwp .AND. &
         ANY(luse_radarfwo(1:n_dom)) .AND. num_io_procs_radar > 0   ) THEN

      ! Only workers reach this point here. Send
      !  grid_generatingCenter and grid_generatingSubcenter to the radar I/O PEs:
      CALL exchg_with_detached_emvorado_io(grid_generatingCenter, grid_generatingSubcenter)

    END IF
#endif
#endif


#ifndef __NO_JSBACH__
    ! Setup horizontal grids and tiles for JSBACH
    DO jg=1,n_dom
      IF (aes_phy_config(jg)%ljsb .OR. atm_phy_nwp_config(jg)%inwp_surface == LSS_JSBACH) THEN
        CALL jsbach_setup_grid( jg, p_patch(jg), type='icon') !< in
        CALL jsbach_setup_tiles(jg)
      END IF
    END DO
#endif

#ifdef MESSY
    CALL messy_initialize(n_dom)
    CALL messy_new_tracer
#endif

    iart_ntracer = 0
#ifdef __ICON_ART
    IF (lart) THEN
      ! determine number of ART-tracers (by reading given XML-Files)
      ! * art_config(1)%iart_ntracer
      CALL art_calc_ntracer_and_names()
      iart_ntracer = art_config(1)%iart_ntracer
    END IF
#endif

#ifndef __NO_ICON_COMIN__
    ! expose descriptive data structures
    CALL icon_expose_descrdata( &
      &  msg_level = msg_level, &
      &  n_dom = n_dom, &
      &  max_dom = max_dom, &
      &  nproma = nproma, &
      &  min_rlcell_int = min_rlcell_int, &
      &  min_rlcell = min_rlcell, &
      &  max_rlcell = max_rlcell, &
      &  min_rlvert_int = min_rlvert_int, &
      &  min_rlvert = min_rlvert, &
      &  max_rlvert = max_rlvert, &
      &  min_rledge_int = min_rledge_int, &
      &  min_rledge = min_rledge, &
      &  max_rledge = max_rledge, &
      &  grf_bdywidth_c = grf_bdywidth_c, &
      &  grf_bdywidth_e = grf_bdywidth_e, &
      &  lrestart = isRestart(), &
      &  vct_a = vct_a, &
      &  p_patch = p_patch(1:), &
      &  number_of_grid_used = number_of_grid_used, &
      &  vgrid_buffer = vgrid_buffer, &
      &  start_time = start_time, &
      &  end_time = end_time, &
      &  time_config = time_config, &
      &  dtime = dtime )
#endif

    CALL init_tracer_settings(iforcing, n_dom, ltransport,                 &
      &                       atm_phy_nwp_config(:)%inwp_turb,             &
      &                       atm_phy_nwp_config(:)%inwp_gscp,             &
      &                       atm_phy_nwp_config(:)%inwp_convection,       &
      &                       lart, iart_ntracer, ctracer_art,             &
      &                       advection_config,                            &
      &                       iqv, iqc, iqi, iqr, iqs, iqt, iqg, iqni,     &
      &                       iqh, iqnr, iqns, iqng, iqnh, iqnc,           &
      &                       iqgl, iqhl, inccn, ininact, ininpot,         &
      &                       iqtke, iqm_max, ntracer, nqtendphy,          &
      &                       atm_phy_nwp_config(:)%nclass_gscp,           &
      &                       iqbin, iqb_i, iqb_e, iqb_s)

    IF (lmoist_thdyn .AND. .NOT.(iqv>0)) THEN
      CALL finish( routine, 'Trying to run moist thermodynamics without moisture')
    ENDIF

#ifdef __ICON_ART
    !------------------------------------------------------------------
    ! 11. Create ART data fields
    !------------------------------------------------------------------

    CALL art_init_interface(n_dom,'construct',ntracer,advection_config(1)%npassive_tracer)
#endif

    !------------------------------------------------------------------

    IF (timers_level > 1) CALL timer_stop(timer_model_init)

  END SUBROUTINE construct_atmo_model

  !-------------------------------------------------------------------
  !>
  SUBROUTINE destruct_atmo_model()

    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:destruct_atmo_model"

    INTEGER :: error_status

    ! Destruct external data state
    CALL destruct_ext_data()

    ! destruct interpolation patterns generate in create_grf_index_lists
    IF (n_dom_start==0 .OR. n_dom > 1) THEN
      CALL destruct_interpol_patterns(p_patch)
    END IF

    ! Deconstruct grid refinement state
    IF ( n_dom_start==0 .OR. n_dom > 1) THEN
      CALL destruct_2d_gridref_state( p_patch, p_patch_local_parent, p_grf_state, &
        &                             p_grf_state_local_parent )

      ! Note that p_grf_state_local_parent is not allocated for n_dom=1 without radiation grid
      DEALLOCATE (p_grf_state_local_parent, STAT=error_status)
      IF (error_status /= SUCCESS) THEN
        CALL finish(routine, 'deallocation for p_grf_state_local_parent failed')
      ENDIF
    ENDIF

    DEALLOCATE (p_grf_state, STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(routine, 'deallocation for ptr_grf_state failed')
    ENDIF

    ! Deallocate interpolation fields
    CALL destruct_2d_interpol_state( p_int_state )
    IF (msg_level>5) CALL message(routine,'destruct_2d_interpol_state is done')

    DEALLOCATE (p_int_state, STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(routine, 'deallocation for ptr_int_state failed')
    ENDIF

    ! Deallocate global registry for lon-lat grids
    CALL lonlat_grids%finalize()

    ! Destruct communication patterns
    CALL destruct_comm_patterns( p_patch, p_patch_local_parent )

    ! Deallocate grid patches
    CALL destruct_patches( p_patch )
    CALL destruct_patches( p_patch_local_parent )
    IF (msg_level>5) CALL message(routine, 'destruct_patches is done')

    DEALLOCATE( p_patch, STAT=error_status )
    IF (error_status/=SUCCESS) THEN
      CALL finish(routine, 'deallocate for patch array failed')
    ENDIF

    ! destruct basic coupler
    IF (is_coupled_run()) THEN
      CALL cpl_destruct()
    END IF

    ! close memory logging files
    CALL memory_log_terminate

!    IF (use_icon_comm) THEN
      CALL destruct_icon_communication()
!    ENDIF

#ifdef __ICON_ART
    ! Destruct ART data fields
    CALL art_init_interface(n_dom,'destruct',ntracer,advection_config(1)%npassive_tracer)
#endif

#ifndef __NO_JSBACH__
    CALL jsbach_finalize()
#endif
    CALL message(routine, 'clean-up finished')

  END SUBROUTINE destruct_atmo_model
  !-------------------------------------------------------------------

END MODULE mo_atmo_model
