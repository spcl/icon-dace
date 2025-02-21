! @brief Main program for the ICON atmospheric model
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

MODULE mo_atmo_model

  ! basic modules
  USE mo_exception,               ONLY: message, finish
  USE mo_mpi,                     ONLY: set_mpi_work_communicators,       &
    &                                   my_process_is_pref, process_mpi_pref_size,   &
    &                                   my_process_is_work,      &
    &                                   my_process_is_mpi_test
  USE mo_timer,                   ONLY: init_timer, timer_start, timer_stop,                  &
    &                                   timers_level, timer_model_init,                       &
    &                                   timer_domain_decomp, timer_compute_coeffs,            &
    &                                   timer_ext_data, print_timer
  USE mo_parallel_config,         ONLY: p_test_run, num_test_pe, l_test_openmp, num_io_procs, &
    &                                   proc0_shift, num_prefetch_proc, pio_type, num_io_procs_radar, &
    &                                   ignore_nproma_use_nblocks_c, nproma, update_nproma_for_io_procs
  USE mo_master_config,           ONLY: isRestart
  USE mo_memory_log,              ONLY: memory_log_terminate
  USE mo_impl_constants,          ONLY: SUCCESS, inh_atmosphere, inwp, LSS_JSBACH, LSS_TERRA,  &
       &                                   min_rlcell_int, min_rlcell
  USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_zaxis_type,              ONLY: zaxisTypeList, t_zaxisTypeList
  USE mo_load_restart,            ONLY: read_restart_header

  ! namelist handling; control parameters: run control, dynamics
  USE mo_read_namelists,          ONLY: read_atmo_namelists
  USE mo_nml_crosscheck,          ONLY: atm_crosscheck
  USE mo_initicon_config,         ONLY: configure_initicon
  USE mo_io_config,               ONLY: restartWritingParameters
  USE mo_lnd_nwp_config,          ONLY: configure_lnd_nwp
  USE mo_dynamics_config,         ONLY: configure_dynamics, iequations, lmoist_thdyn
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
  USE mo_aes_phy_config,          ONLY: aes_phy_config
  USE mo_master_control,          ONLY: master_namelist_filename
  USE mo_jsb_base,                ONLY: jsbach_setup => jsbach_setup_models, jsbach_setup_tiles
  USE mo_jsb_model_init,          ONLY: jsbach_setup_grid
  USE mo_jsb_model_final,         ONLY: jsbach_finalize
  USE mo_master_control,          ONLY: atmo_process

  ! time stepping
  USE mo_atmo_nonhydrostatic,     ONLY: atmo_nonhydrostatic, construct_atmo_nonhydrostatic

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
  USE mo_ext_data_state,          ONLY: ext_data, destruct_ext_data
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

  ! coupling

  ! I/O
  USE mo_restart,                 ONLY: detachRestartProcs
  USE mo_icon_output_tools,       ONLY: init_io_processes
  USE mtime,                      ONLY: OPERATOR(<), OPERATOR(+)
  USE mo_async_latbc_types,       ONLY: t_latbc_data
  ! ART
  USE mo_art_config,              ONLY: ctracer_art
  USE mo_sync,                    ONLY: global_max


  !-------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: atmo_model, construct_atmo_model, destruct_atmo_model

CONTAINS


  !-------------------------------------------------------------------
  !>
  SUBROUTINE atmo_model(atm_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: atm_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:atmo_model"

    TYPE(t_latbc_data) :: latbc !< data structure for async latbc prefetching


    !---------------------------------------------------------------------
    ! construct the atmo model
    CALL construct_atmo_model(atm_namelist_filename,shr_namelist_filename)

    SELECT CASE(iequations)

    CASE(inh_atmosphere)
      CALL construct_atmo_nonhydrostatic(latbc)

    CASE DEFAULT
      CALL finish(routine, 'unknown choice for iequations.')
    END SELECT

    !---------------------------------------------------------------------
    ! construct the coupler
    !


    !---------------------------------------------------------------------
    ! 12. The hydrostatic model has been deleted. Only the non-hydrostatic
    !     model is available.
    !---------------------------------------------------------------------
    SELECT CASE(iequations)
    CASE(inh_atmosphere)
      CALL atmo_nonhydrostatic(latbc)

    CASE DEFAULT
      CALL finish(routine, 'unknown choice for iequations.')
    END SELECT

    ! print performance timers:
    IF (ltimer) CALL print_timer


    !---------------------------------------------------------------------
    ! 13. Integration finished. Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_model ()

    !---------------------------------------------------------------------
    ! (optional:) write resident set size from OS

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
      CALL read_restart_header("atm")
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


    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some components of the state, e.g., num_lev, may be
    !    modified in this subroutine which affects the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run( )


    ! complete initicon config-state
    CALL configure_initicon()


    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL restartWritingParameters(opt_dedicatedProcCount = dedicatedRestartProcs)


    CALL set_mpi_work_communicators(p_test_run, l_test_openmp,                    &
         &                          num_io_procs, dedicatedRestartProcs,          &
         &                          atmo_process, num_prefetch_proc, num_test_pe, &
         &                          pio_type,                                     &
         &                          num_io_procs_radar=num_io_procs_radar,        &
         &                          radar_flag_doms_model=luse_radarfwo(1:n_dom), &
         &                          num_dio_procs=proc0_shift)




    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer
    IF (timers_level > 1) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! initialize dynamic list of vertical axes
    !-------------------------------------------------------------------

    zaxisTypeList = t_zaxisTypeList()

    ! Setup JSBACH: read namelists, configure models for each domain
    ! This has to be after (!) the ICON zaxes have been created in the above line but
    ! before (!) the restart PEs are detached a few lines below since JSBACH
    ! adds its zaxes to zaxisTypeList
    IF (ANY(aes_phy_config(:)%ljsb) &
        & .OR. ANY(atm_phy_nwp_config(1:n_dom)%inwp_surface == LSS_JSBACH)) THEN
      ! Do basic initialization of JSBACH
      CALL jsbach_setup(master_namelist_filename)
    END IF

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

    IF (ignore_nproma_use_nblocks_c) THEN
      nproma_max = global_max(nproma)
      CALL update_nproma_for_io_procs(nproma_max)
    ENDIF

    IF (timers_level > 4) CALL timer_stop(timer_domain_decomp)

    !-------------------------------------------------------------------
    ! 5. I/O initialization
    !-------------------------------------------------------------------

    ! This won't RETURN on dedicated restart PEs, starting their main loop instead.
    CALL detachRestartProcs(timers_level > 1)



    CALL init_io_processes()


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
    ALLOCATE (ext_data(n_dom), STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(routine, 'allocation for ext_data failed')
    ENDIF

    ! allocate memory for atmospheric/oceanic external data and
    ! optionally read those data from netCDF file.
    IF (timers_level > 4) CALL timer_start(timer_ext_data)
    CALL init_ext_data (p_patch(1:), p_int_state(1:), ext_data)
    IF (timers_level > 4) CALL timer_stop(timer_ext_data)

    !---------------------------------------------------------------------
    ! Import vertical grid/ define vertical coordinate
    !---------------------------------------------------------------------

    CALL allocate_vct_atmo(p_patch(1)%nlevp1)
    IF (iequations == inh_atmosphere .AND. ltestcase .AND. (.NOT. l_scm_mode)) THEN
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



    ! Setup horizontal grids and tiles for JSBACH
    DO jg=1,n_dom
      IF (aes_phy_config(jg)%ljsb .OR. atm_phy_nwp_config(jg)%inwp_surface == LSS_JSBACH) THEN
        CALL jsbach_setup_grid( jg, p_patch(jg), type='icon') !< in
        CALL jsbach_setup_tiles(jg)
      END IF
    END DO


    iart_ntracer = 0


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


    !------------------------------------------------------------------

    IF (timers_level > 1) CALL timer_stop(timer_model_init)

  END SUBROUTINE construct_atmo_model

  !-------------------------------------------------------------------
  !>
  SUBROUTINE destruct_atmo_model()

    CHARACTER(*), PARAMETER :: routine = "mo_atmo_model:destruct_atmo_model"

    INTEGER :: error_status
    ! Destruct external data state

    CALL destruct_ext_data
    IF (msg_level > 5) CALL message(routine, 'destruct_ext_data is done')

    ! deallocate ext_data array
    DEALLOCATE(ext_data, stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(routine, 'deallocation of ext_data')
    ENDIF


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

    ! close memory logging files
    CALL memory_log_terminate

!    IF (use_icon_comm) THEN
      CALL destruct_icon_communication()
!    ENDIF


    CALL jsbach_finalize()
    CALL message(routine, 'clean-up finished')

  END SUBROUTINE destruct_atmo_model
  !-------------------------------------------------------------------

END MODULE mo_atmo_model
