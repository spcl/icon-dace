! Testbed for modifications requiresd to enable zstar
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

!----------------------------
#include "icon_definitions.inc"
!----------------------------
MODULE mo_ocean_testbed_zstar
  !-------------------------------------------------------------------------
 
  USE mo_kind,                      ONLY: wp
  USE mo_parallel_config,           ONLY: nproma
  USE mo_sync,                      ONLY: sync_e, sync_c, sync_c1, sync_patch_array, &
    & sync_patch_array_mult, global_sum_array

  USE mo_impl_constants,            ONLY: sea_boundary, max_char_length, min_dolic
  USE mo_dbg_nml,                   ONLY: idbg_mxmn
  USE mo_ocean_nml, ONLY: n_zlev, solver_tolerance,&
    & l_with_vert_tracer_advection, OceanReferenceDensity, OceanReferenceDensity_inv, &
    & para_surfRelax_Temp, para_surfRelax_Salt, l_relaxsal_ice, i_sea_ice, &
    & type_surfRelax_Temp, type_surfRelax_Salt, Analytical_Forcing, &
    & OMIP_FluxFromFile, lhamocc, lfb_bgc_oce, lswr_jerlov, &
    & limit_elevation, &
    & ab_const, ab_beta, ab_gam, iswm_oce, iforc_oce, &
    & no_tracer, l_rigid_lid, l_edge_based,               &
    & use_absolute_solver_tolerance, solver_max_restart_iterations, &
    & solver_max_iter_per_restart, dhdtw_abort, select_transfer, &
    & select_solver, select_gmres, select_gmres_r, select_mres, &
    & select_gmres_mp_r, select_cg, select_cgj, select_bcgs, &
    & select_legacy_gmres, use_continuity_correction, select_cg_mp, &
    & solver_max_iter_per_restart_sp, solver_tolerance_sp, No_Forcing, &
    & MASS_MATRIX_INVERSION_TYPE,            &
    & MASS_MATRIX_INVERSION_ADVECTION, solver_tolerance_comp, &
    & MASS_MATRIX_INVERSION_ALLTERMS, &
    & PPscheme_type, PPscheme_ICON_Edge_vnPredict_type, &
    & solver_FirstGuess, MassMatrix_solver_tolerance,     &
    & createSolverMatrix, l_solver_compare, solver_comp_nsteps, &
    & Cartesian_Mixing, GMRedi_configuration, atm_pressure_included_in_ocedyn, &
    & cfl_check, cfl_threshold, cfl_stop_on_violation,   &
    & cfl_write, Coupled_FluxFromAtmo, use_tides
  USE mo_run_config,                ONLY: dtime, debug_check_level, nsteps, output_mode
  USE mo_timer, ONLY: timer_start, timer_stop, timers_level, timer_extra1, &
    & timer_extra2, timer_extra3, timer_extra4, timer_ab_expl, timer_ab_rhs4sfc, &
    & timer_total, timer_solve_ab, timer_upd_flx, timer_vert_veloc, &
    & timer_extra20, timer_extra21, timer_normal_veloc, timer_scalar_prod_veloc

  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_physical_constants,        ONLY: grav, clw, rho_ref, Tf
  USE mo_ocean_initialization,      ONLY: is_initial_timestep
  USE mo_ocean_types, ONLY: t_hydro_ocean_state
  USE mo_ocean_time_events,   ONLY: ocean_time_nextStep, isCheckpoint, isEndOfThisRun, newNullDatetime
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,            ONLY: t_external_data
  USE mo_exception,                 ONLY: message, finish, warning, message_text
  USE mo_util_dbg_prnt,             ONLY: dbg_print, debug_print_MaxMinMean, debug_printValue
  USE mo_ocean_boundcond,           ONLY: VelocityBottomBoundaryCondition_onBlock, top_bound_cond_horz_veloc
  USE mo_ocean_thermodyn,           ONLY: calculate_density, calc_internal_press_grad, calc_potential_density
  USE mo_ocean_physics_types,       ONLY: t_ho_params
  USE mo_ocean_pp_scheme,           ONLY: ICON_PP_Edge_vnPredict_scheme
  USE mo_ocean_surface_types,       ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_scalar_product, ONLY: map_edges2edges_viacell_3d_const_z, map_edges2edges_viacell_2D_per_level, &
    & calc_scalar_product_veloc_3d, map_edges2edges_sc_zstar, map_edges2edges_3d_zstar
  USE mo_ocean_math_operators, ONLY: div_oce_3D_onTriangles_onBlock, smooth_onCells, div_oce_3d, &
    & grad_fd_norm_oce_2d_onblock, grad_fd_norm_oce_2d_3d, div_oce_3D_general_onBlock, &
    & div_oce_3D_onTriangles_onBlock, update_height_depdendent_variables, &
    & check_cfl_horizontal, check_cfl_vertical

  USE mo_ocean_velocity_advection, ONLY: veloc_adv_horz_mimetic, veloc_adv_vert_mimetic
  USE mo_ocean_velocity_diffusion, ONLY: velocity_diffusion, velocity_diffusion_vertical_implicit_onBlock
  USE mo_ocean_types, ONLY: t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_grid_subset, ONLY: t_subset_range, get_index_range
  USE mo_grid_config, ONLY: n_dom
  USE mo_mpi, ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_statistics, ONLY: global_minmaxmean, print_value_location
  USE mo_ocean_solve, ONLY: t_ocean_solve
  USE mo_ocean_solve_lhs_type, ONLY: t_lhs_agen
  USE mo_ocean_solve_transfer, ONLY: t_transfer
  USE mo_ocean_solve_trivial_transfer, ONLY: t_trivial_transfer
  USE mo_ocean_solve_subset_transfer, ONLY: t_subset_transfer
  USE mo_ocean_solve_aux, ONLY: t_ocean_solve_parm, solve_gmres, solve_cg, solve_mres, &
   & solve_precon_none, solve_precon_jac, solve_bcgs, solve_legacy_gmres, &
   & solve_trans_scatter, solve_trans_compact, solve_cell, solve_edge, solve_invalid
  USE mo_primal_flip_flop_lhs, ONLY: t_primal_flip_flop_lhs
  USE mo_surface_height_lhs, ONLY: t_surface_height_lhs
  USE mo_surface_height_lhs_zstar, ONLY: t_surface_height_lhs_zstar
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_sea_ice_types,          ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_name_list_output_init,  ONLY: isRegistered
  USE mo_hamocc_types,           ONLY: t_hamocc_state
  USE mtime,                     ONLY: datetime, timedelta, newTimedelta, getNoOfSecondsElapsedInDayDateTime, &
    & datetimeToString, deallocateTimedelta
  USE mo_ocean_tracer_transport_types,  ONLY: t_ocean_transport_state, t_ocean_tracer, t_tracer_collection
  USE mo_math_constants,         ONLY: pi, pi_2, rad2deg, deg2rad, dbl_eps
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff, no_primal_edges
  USE mo_ocean_tracer_transport_vert, ONLY: advect_flux_vertical 
  USE mo_memory_log,             ONLY: memory_log_add 
  USE mo_ocean_surface_refactor, ONLY: update_ocean_surface_refactor
  
  USE mo_ocean_physics,         ONLY: update_ho_params
  USE mo_ocean_ab_timestepping_mimetic,  ONLY: calculate_explicit_term_ab, fill_rhs4surface_eq_ab, &
    & solve_free_sfc_ab_mimetic
  USE mo_hydro_ocean_run, ONLY: update_time_g_n, update_time_indices
  USE mo_ocean_output, ONLY: output_ocean
  USE mo_ocean_ab_timestepping,  ONLY: calc_vert_velocity, calc_normal_velocity_ab
  USE mo_ocean_tracer,           ONLY: advect_ocean_tracers
  USE mo_ocean_diagnostics,      ONLY: diag_heat_salt_tendency
  USE mo_ocean_bulk_forcing,  ONLY: apply_surface_relaxation, update_ocean_surface_stress
  USE mo_swr_absorption,     ONLY: dynamic_swr_absorption
  USE mo_ocean_diagnostics,      ONLY: calc_fast_oce_diagnostics, calc_psi
  USE mo_derived_variable_handling, ONLY: update_statistics
  USE mo_master_config,          ONLY: isRestart
  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_restart,                ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
  USE mo_ice_fem_interface,      ONLY: ice_fem_init_vel_restart, ice_fem_update_vel_restart
  USE mo_ocean_ab_timestepping_zstar,  ONLY: calc_vert_velocity_bottomup_zstar, calc_normal_velocity_ab_zstar, &
    & solve_free_surface_eq_zstar, update_zstar_variables
  USE mo_ocean_tracer_zstar, ONLY:advect_individual_tracers_zstar, advect_ocean_tracers_zstar
  USE mo_swr_absorption, ONLY: subsurface_swr_absorption_zstar
  USE mo_ocean_tracer_dev, ONLY: advect_ocean_tracers_GMRedi_zstar
  USE mo_ocean_nudging,          ONLY: nudge_ocean_tracers
!  USE mo_ocean_atmo_coupling,   ONLY: couple_ocean_toatmo_fluxes
  !-------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ocean_test_zstar_advection
  PUBLIC :: test_stepping_zstar 
  PUBLIC :: test_stepping_z
 
  INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug

  !-------------------------------------------------------------------------
CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! FIXME: Notes and overview
  !! Routines for testing advection with zstar are below
  !! Routines have been copied from other files and modified for zstar
  !! 1. Ideally, coefficients for the modified discretization should be calculated at one place
  !! and passed as arguments to repeat calculations
  !! 2. One routine can be used to calculate both low and high order flux for speedup
  !! 3. Variable to calculate depth needs to be clarified
  !! 4. This can be converted to using generalized vertical co-ordinates by using only 
  !! the coefficient of dz as variables to be modified
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !> Setup a test case that advects tracers for testing with zstar
  !  Should start by using the low order horizontal advection
  SUBROUTINE ocean_test_zstar_advection( patch_3d, ocean_state, &
    & this_datetime, ocean_surface, physics_parameters,             &
    & ocean_ice,operators_coefficients)
    
    TYPE(t_patch_3d), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: ocean_surface
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(in)          :: operators_coefficients
    
    ! local variables
    TYPE (t_hamocc_state)        :: hamocc_State
    INTEGER :: jstep, jg
    !LOGICAL                         :: l_outputtime
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0 ! start counter for time loop
    INTEGER :: i
    INTEGER :: tracer_index 
    TYPE(timedelta), POINTER :: model_time_step => NULL()
    
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection
    TYPE(t_ocean_transport_state)                    :: transport_state
    
    INTEGER  :: jb, jc, je, level 
    REAL(wp) :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: delta_t, delta_z,delta_z_new, delta_z1,delta_z_new1
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: flux_horz(nproma,n_zlev, patch_3d%p_patch_2D(1)%nblks_e)
    REAL(wp) :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: top_bc(nproma)
    INTEGER  :: start_index, end_index
    INTEGER  :: start_cell_index, end_cell_index
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_low (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_high(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z2(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER  :: bt_level 
    REAL(wp) :: H_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp) :: stretch_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: st1, st2 
    REAL(wp) :: eta_0(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: eta_1(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER, DIMENSION(:,:,:), POINTER :: idx, blk
    INTEGER  :: id1, id2, bl1, bl2 

    TYPE(t_ocean_tracer), POINTER :: new_tracer
    TYPE(t_ocean_tracer), POINTER :: old_tracer

    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    
    REAL(wp) :: temp(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & method_name = 'mo_ocean_testbed_modules:ocean_test_zstar_advection'
    !------------------------------------------------------------------
    
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2d%edges%in_domain
    idx             => patch_3D%p_patch_2D(1)%edges%cell_idx
    blk             => patch_3D%p_patch_2D(1)%edges%cell_blk
 
    CALL datetimeToString(this_datetime, datestring)

    ! IF (ltimer) CALL timer_start(timer_total)
    CALL timer_start(timer_total)
    
    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(method_name), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    !! sea surface height type 201 set explicitly since we don't want
    !! the grid to change
    eta = 0.  
    ! Initialize eta for zstar
    ! #slo#: simple elevation between 30W and 30E (pi/3.)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        IF ( patch_3d%lsm_c(jc, 1, jb) <= sea_boundary ) THEN
           eta(jc, jb) = 10.0_wp * &
            & SIN(patch_2d%cells%center(jc, jb)%lon * 6.0_wp) &
            & * COS(patch_2d%cells%center(jc, jb)%lat * 3.0_wp)
        ENDIF
      END DO
    END DO
    
    !------------------------------------------------------------------
    stretch_c = 1.0_wp
    stretch_e = 1.0_wp
    !------------------------------------------------------------------
 
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, bt_lev) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        
        bt_level = patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      
 
    !------------------------------------------------------------------
        !! Initialize only as a placeholder to call subroutine
        eta_0(jc, jb)      = eta(jc, jb)
        eta_1(jc, jb)      = eta(jc, jb)
    !------------------------------------------------------------------
        eta_c(jc, jb)      = eta(jc, jb)
        H_c  (jc, jb)      = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_level + 1, jb)
        if ( patch_3D%lsm_c(jc, 1, jb) <= sea_boundary ) THEN
          stretch_c(jc, jb)  = (H_c(jc, jb) + eta_c(jc, jb))/H_c(jc, jb) 
        else 
          stretch_c(jc, jb)  = 1.0_wp
        ENDIF
        stretch_c_new(jc, jb) = stretch_c(jc, jb)

      END DO
    END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_MASTER
    CALL sync_patch_array(sync_c, patch_2D, stretch_c)
    CALL sync_patch_array(sync_c, patch_2D, stretch_c_new)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER


!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, je, jk, id1, id2, bl1, bl2, st1, st2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      DO je = start_index, end_index
        id1 = idx(je, jb, 1)
        id2 = idx(je, jb, 2)
        bl1 = blk(je, jb, 1)
        bl2 = blk(je, jb, 2)
 
        st1 = stretch_c(id1, bl1) 
        st2 = stretch_c(id2, bl2) 

        !! FIXME: There seem to be edge cases where this does not work
        IF(patch_3D%lsm_e(je, 1, jb) <= sea_boundary)THEN
          stretch_e(je, jb) = 0.5_wp*(st1 + st2)
        ELSE
          stretch_e(je, jb) = 1.0_wp
        ENDIF


      ENDDO
    END DO
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_MASTER
    CALL sync_patch_array(sync_e, patch_2D, stretch_e)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER




!    !---------------------------------------------------------------------
!    !-FIXME: test divergence of constant fn
!    !---------------------------------------------------------------------
!  
!    ! calc_vert_vel uses vn_time_weighter instead of vn
!    ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn
!    
!    !! Update mass_flux and w 
!    CALL calc_vert_velocity_bottomup_zstar( patch_3d, ocean_state(jg), operators_coefficients, &
!      & stretch_c, stretch_e, eta_0, eta_1)
!
!    ! fill transport_state
!    transport_state%patch_3d    => patch_3d
!    transport_state%h_old       => ocean_state(jg)%p_prog(nold(1))%h
!    transport_state%h_new       => ocean_state(jg)%p_prog(nnew(1))%h
!    transport_state%vn          => ocean_state(jg)%p_prog(nold(1))%vn
!    transport_state%mass_flux_e => ocean_state(jg)%p_diag%mass_flx_e
!    transport_state%w           => ocean_state(jg)%p_diag%w
!
!    temp = 1.0_wp
!
!    CALL advect_flux_vertical( patch_3d,&
!        & temp, &
!        & transport_state,                           &
!        & operators_coefficients,                     &
!        & div_adv_flux_vert)
!
!    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
!      DO jc = start_cell_index, end_cell_index
!        div_adv_flux_vert(jc, :, jb) = stretch_c(jc, jb)*div_adv_flux_vert(jc, :, jb)
!      END DO
!    END DO
!
!    !---------------------------------------------------------------------
!    !-Horizontal  advection
!    !---------------------------------------------------------------------
!    CALL upwind_zstar_hflux_oce( patch_3d,  &
!      & temp, &
!      & transport_state%mass_flux_e,         &
!      & z_adv_flux_h)                         
! 
!    !Calculate divergence of advective fluxes
!    CALL div_oce_3d( z_adv_flux_h, patch_3D, operators_coefficients%div_coeff, &
!      & div_adv_flux_horz, subset_range=cells_in_domain )
!
!    call map_edges2edges_sc_zstar( patch_3d, transport_state%vn, temp, &
!      & operators_coefficients, stretch_e, z_adv_flux_h)
!!
!    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
!      DO jc = start_cell_index, end_cell_index
!        DO level = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
!
!          if ( (jb == 4) .AND. (jc == 10) ) THEN
!!            write(*, *) level, div_adv_flux_horz(jc,level,jb) , temp(jc,level,jb)
!!              & transport_state%w(jb, level + 1, jc) - transport_state%w(jb, level, jc), &
!!              & patch_3d%p_patch_1d(1)%prism_thick_c(jb, level, jc)
!          ENDIF
!    
!        ENDDO
!      ENDDO
!    ENDDO
! 
!    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
!      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
!      DO je = start_index, end_index
!        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
!           if ( (jb == 4) .AND. (je == 10) ) THEN
!!            write(*, *) level, transport_state%w(je,level,jb) 
!          ENDIF
!    
!        ENDDO
!      ENDDO
!    ENDDO
!
!
!    !---------------------------------------------------------------------
!    !-FIXME: test end 
!    !---------------------------------------------------------------------
 
    jstep0 = 0
    DO jstep = (jstep0+1), (jstep0+nsteps)

        CALL datetimeToString(this_datetime, datestring)
        WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
        CALL message (TRIM(method_name), message_text)
 
        old_tracer_collection => ocean_state(jg)%p_prog(nold(1))%tracer_collection
        new_tracer_collection => ocean_state(jg)%p_prog(nnew(1))%tracer_collection

        IF (no_tracer>=1) THEN

          ! calc_vert_vel uses vn_time_weighter instead of vn
          ocean_state(jg)%p_diag%vn_time_weighted = ocean_state(jg)%p_prog(nold(1))%vn

          !! Update w and mass_flx_e for tracer advection
          !! FIXME: Needs to be tested with correct vertical velocity
          CALL calc_vert_velocity_bottomup_zstar( patch_3d, ocean_state(jg), operators_coefficients, &
            & stretch_c, stretch_e)

          ! fill transport_state
          transport_state%patch_3d    => patch_3d
          transport_state%h_old       => ocean_state(jg)%p_prog(nold(1))%h
          transport_state%h_new       => ocean_state(jg)%p_prog(nnew(1))%h
          transport_state%vn          => ocean_state(jg)%p_prog(nold(1))%vn
          transport_state%w           => ocean_state(jg)%p_diag%w
          transport_state%mass_flux_e => ocean_state(jg)%p_diag%mass_flx_e

          ! fill diffusion coefficients
          old_tracer_collection%tracer(1)%hor_diffusion_coeff => physics_parameters%TracerDiffusion_coeff(:,:,:,1)
          old_tracer_collection%tracer(1)%ver_diffusion_coeff => physics_parameters%a_tracer_v(:,:,:,1)
          DO i = 2, old_tracer_collection%no_of_tracers
            old_tracer_collection%tracer(i)%hor_diffusion_coeff => physics_parameters%TracerDiffusion_coeff(:,:,:,2)
            old_tracer_collection%tracer(i)%ver_diffusion_coeff => physics_parameters%a_tracer_v(:,:,:,2)
          ENDDO
        ENDIF
        !------------------------------------------------------------------------

        DO tracer_index = 1, old_tracer_collection%no_of_tracers
          
          old_tracer => old_tracer_collection%tracer(tracer_index)
          new_tracer => new_tracer_collection%tracer(tracer_index)
          IF ( old_tracer%is_advected) THEN
           
            call advect_individual_tracers_zstar( patch_3d, transport_state, &
              & operators_coefficients, stretch_e, stretch_c, stretch_c_new, old_tracer, new_tracer)

          ENDIF
        END DO

        ! One integration cycle finished on the lowest grid level (coarsest
        ! resolution). Set model time.
        model_time_step => newTimedelta('+', 0, 0, 0, 0, 0, NINT(dtime), 0)
!        this_datetime = this_datetime + model_time_step
        CALL deallocateTimedelta(model_time_step) 
          
        CALL output_ocean( patch_3d, &
          & ocean_state,             &
          & this_datetime,                &
          & ocean_surface,          &
          & ocean_ice,               &
          & jstep, jstep0)

        ! Shift time indices for the next loop
        ! this HAS to ge into the restart files, because the start with the following loop
        CALL update_time_indices(jg)
        ! update intermediate timestepping variables for the tracers
        ! velocity
        ocean_state(jg)%p_aux%g_nm1 = ocean_state(jg)%p_aux%g_n
        ocean_state(jg)%p_aux%g_n   = 0.0_wp
      
        CALL update_time_g_n(ocean_state(jg))

    END DO
    
    CALL timer_stop(timer_total)
    
  END SUBROUTINE ocean_test_zstar_advection
  !-------------------------------------------------------------------------

 

 
   

  !-------------------------------------------------------------------------
  SUBROUTINE tracer_transport_zstar(patch_3d, ocean_state, p_as, sea_ice, &
      & p_oce_sfc, p_phys_param, operators_coefficients, current_time, &
      & stretch_c, stretch_e, stretch_c_new)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: p_phys_param
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(datetime), POINTER, INTENT(in)              :: current_time
    REAL(wp), INTENT(IN) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    REAL(wp), INTENT(IN) :: stretch_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
 
    TYPE(t_ocean_transport_state)                    :: transport_state
    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection

    INTEGER :: i

    !------------------------------------------------------------------------
    !Tracer transport
 
    old_tracer_collection => ocean_state%p_prog(nold(1))%tracer_collection
    new_tracer_collection => ocean_state%p_prog(nnew(1))%tracer_collection

 
    !------------------------------------------------------------------------
    IF (no_tracer>=1) THEN
  
      ! fill transport_state
      transport_state%patch_3d    => patch_3d
      transport_state%h_old       => ocean_state%p_prog(nold(1))%h
      transport_state%h_new       => ocean_state%p_prog(nnew(1))%h
      transport_state%w           => ocean_state%p_diag%w  ! w_time_weighted
      transport_state%mass_flux_e => ocean_state%p_diag%mass_flx_e
      transport_state%vn          => ocean_state%p_diag%vn_time_weighted
      ! fill boundary conditions
      old_tracer_collection%tracer(1)%top_bc => p_oce_sfc%TopBC_Temp_vdiff
      IF (no_tracer > 1) &
        old_tracer_collection%tracer(2)%top_bc => p_oce_sfc%TopBC_Salt_vdiff
  
      ! fill diffusion coefficients
      DO i = 1, old_tracer_collection%no_of_tracers
          old_tracer_collection%tracer(i)%hor_diffusion_coeff => p_phys_param%TracerDiffusion_coeff(:,:,:,i)
          old_tracer_collection%tracer(i)%ver_diffusion_coeff => p_phys_param%a_tracer_v(:,:,:,i)
      ENDDO
      
    ENDIF
    !------------------------------------------------------------------------
  
    !------------------------------------------------------------------------
    ! FIXME zstar: GM diffusion not implemented
    ! transport tracers and diffuse them
    IF (no_tracer>=1) THEN

      IF (GMRedi_configuration==Cartesian_Mixing) THEN
        !! Note that zstar has no horizontal diffusion
        CALL advect_ocean_tracers_zstar(old_tracer_collection, new_tracer_collection, &
          & transport_state, operators_coefficients, stretch_e, stretch_c, stretch_c_new)
      ELSE
        CALL  advect_ocean_tracers_GMRedi_zstar(old_tracer_collection, new_tracer_collection, &
          &  ocean_state, transport_state, p_phys_param, operators_coefficients, &
          &  stretch_c, stretch_e, stretch_c_new)
      ENDIF
 
    ENDIF
    !------------------------------------------------------------------------


    !! FIXME zstar: hamocc interface not implemented
!    CALL ocean_to_hamocc_interface(ocean_state, transport_state, &
!      & p_oce_sfc, p_as, sea_ice, p_phys_param, operators_coefficients, current_time)

  END SUBROUTINE tracer_transport_zstar
  !-------------------------------------------------------------------------




  
  SUBROUTINE test_stepping_zstar( patch_3d, ocean_state, p_ext_data,  &
    & this_datetime, p_oce_sfc, p_phys_param, &
    & p_as, p_atm_f, sea_ice, &
    & hamocc_state,operators_coefficients,solvercoeff_sp)
    
    TYPE(t_patch_3d ), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: p_phys_param 
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: p_atm_f
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_hamocc_state),     INTENT(inout)          :: hamocc_state
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(t_solvercoeff_singleprecision), INTENT(inout) :: solvercoeff_sp

    !! FIXME: after latest master merge huge changes in solver
    !! NEeds to be fixed
    
!    ! local variables
!    INTEGER :: jstep, jg
!    INTEGER :: jb, jc, je, bt_lev 
!    INTEGER :: start_index, end_index 
!    INTEGER :: i 
!    CHARACTER(LEN=32)               :: datestring
!    TYPE(t_patch), POINTER :: patch_2d
!    INTEGER :: jstep0 ! start counter for time loop
!    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
!    CHARACTER(LEN=max_char_length), PARAMETER :: &
!      & routine = 'mo_ocean_testbed_modules:test_zstar_core'
!    
!    REAL(wp) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! Surface height at cell
!    REAL(wp) :: H_c  (nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! Column depth at cell 
!    REAL(wp) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! stretch factor 
!    REAL(wp) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e)           !! 
!    REAL(wp) :: stretch_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! stretch factor 
!    REAL(wp) :: st1, st2 
!    REAL(wp) :: eta_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! Surface height after time step 
!    
!    REAL(wp) :: eta_p(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! Surface height at cell
!
!    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection
!    TYPE(t_ocean_transport_state)                    :: transport_state
! 
!    TYPE(datetime), POINTER             :: current_time     => NULL()
!    TYPE(t_ocean_solve), POINTER :: solve, solve_comp
!
!    CLASS(t_destructible), POINTER :: free_sfc_solver => NULL()
!    CLASS(t_destructible), POINTER :: free_sfc_solver_comp => NULL()
!    CLASS(t_destructible), POINTER :: free_sfc_solver_lhs => NULL()
!    CLASS(t_destructible), POINTER :: free_sfc_solver_trans => NULL()
!    CLASS(t_destructible), POINTER :: free_sfc_solver_comp_trans => NULL()
!    TYPE(t_surface_height_lhs_zstar), POINTER :: lhs_sh => NULL()
!!
!    TYPE(t_subset_range), POINTER :: owned_cells, owned_edges
!    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
!    INTEGER, DIMENSION(:,:,:), POINTER :: idx, blk
!    INTEGER  :: id1, id2, bl1, bl2 
!
!    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
!    CLASS(t_RestartDescriptor), POINTER :: restartDescriptor
! 
!    REAL(wp), PARAMETER ::  min_top_height = 0.05_wp !we have to have at least 5cm water on topLevel of sea cells)
!    INTEGER :: n_it, n_it_sp, ret_status
!    INTEGER  :: rn
!    INTEGER :: it, level 
!    REAL(wp) :: verticalMeanFlux(n_zlev+1)
!    REAL(wp) :: mean_height, old_mean_height
!    INTEGER :: return_status
!
!    CHARACTER(LEN=12)  :: str_module = 'zstar_dyn'  ! Output of module for 1 line debug)
!
!    TYPE(t_subset_range), POINTER :: edges_indomain
! 
!    !------------------------------------------------------------------
!    ! no grid refinement allowed here so far
!    !------------------------------------------------------------------
!    IF (n_dom > 1 ) THEN
!      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
!    END IF
!    jg = n_dom
!
!    
!    jstep0 = 0
!
!    !! FIXME zstar: Restart files will not have eta_c information
!    restartAttributes => getAttributesForRestarting()
!    IF (ASSOCIATED(restartAttributes)) THEN
!      ! get start counter for time loop from restart file:
!      jstep0 = restartAttributes%getInteger("jstep")
!    END IF
!    IF (isRestart() .AND. mod(nold(jg),2) /=1 ) THEN
!      ! swap the g_n and g_nm1
!      CALL update_time_g_n(ocean_state(jg))
!    ENDIF
!
!    restartDescriptor => createRestartDescriptor("oce")
!
!    !------------------------------------------------------------------
!    patch_2d        => patch_3d%p_patch_2d(jg)
! 
!    !------------------------------------------------------------------
! 
!    IF ( .NOT. isRestart()  ) THEN
!        ocean_state(jg)%p_prog(nold(1))%stretch_c = 1.0_wp
!    ENDIF
!   
!    !------------------------------------------------------------------
!    jstep  = jstep0
!
!    ! local time var to be passed along, so the global is kept safe 
!    current_time => newNullDatetime()
!    !------------------------------------------------------------------
!
!    !! Start time stepping
!    DO 
!      ! optional memory loggin
!      CALL memory_log_add
!
!      jstep = jstep + 1
!      ! update model date and time mtime based
!      current_time = ocean_time_nextStep()
!  
!      CALL datetimeToString(current_time, datestring)
!      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =', jstep , '  datetime:  ', datestring
!      CALL message (TRIM(routine), message_text)
!     
!      !! Get kinetic energy
!      start_timer(timer_scalar_prod_veloc,2)
!      CALL calc_scalar_product_veloc_3d( patch_3d,  &
!        & ocean_state(jg)%p_prog(nold(1))%vn,         &
!        & ocean_state(jg)%p_diag,                     &
!        & operators_coefficients)
!       stop_timer(timer_scalar_prod_veloc,2)
!     
!      !! Updates velocity, tracer boundary condition
!      !! Changes height based on ice etc
!      !! Ice eqn sends back heat fluxes and volume fluxes
!      !! Tracer relaxation and surface flux boundary conditions
!      start_timer(timer_upd_flx,3)
!      CALL update_ocean_surface_refactor( patch_3d, ocean_state(jg), p_as, sea_ice, p_atm_f, p_oce_sfc, &
!           & current_time, operators_coefficients, ocean_state(jg)%p_prog(nold(1))%eta_c, &
!           & ocean_state(jg)%p_prog(nold(1))%stretch_c)
!      stop_timer(timer_upd_flx,3)
!
!      !------------------------------------------------------------------
! 
!      !! Update stretch variables 
!      CALL update_zstar_variables( patch_3d, ocean_state(jg), operators_coefficients, &
!        & ocean_state(jg)%p_prog(nold(1))%eta_c, &
!        & ocean_state(jg)%p_prog(nold(1))%stretch_c, stretch_e)
!
!      !---------------------------------------------------------------------
!      ! compute tidal potential
!      ! FIXME zstar: Not modified for zstar
!      IF (use_tides) THEN
!        CALL tide(patch_3d,current_time,ocean_state(jg)%p_aux%bc_tides_potential)
!      ENDIF
!      ! total top potential
!      IF (atm_pressure_included_in_ocedyn) THEN
!        ocean_state(jg)%p_aux%bc_total_top_potential = ocean_state(jg)%p_aux%bc_tides_potential  &
!          & + p_as%pao * OceanReferenceDensity_inv
!      ELSE
!        ocean_state(jg)%p_aux%bc_total_top_potential = ocean_state(jg)%p_aux%bc_tides_potential
!      ENDIF 
!      !---------DEBUG DIAGNOSTICS-------------------------------------------
!      idt_src=3  ! output print level (1-5, fix)
!      CALL dbg_print('on entry: h-old'           ,ocean_state(jg)%p_prog(nold(1))%eta_c ,str_module,idt_src, &
!        & patch_2d%cells%owned )
!      CALL dbg_print('on entry: h-new'           ,ocean_state(jg)%p_prog(nnew(1))%eta_c ,str_module,idt_src, &
!        & patch_2d%cells%owned )
!      CALL dbg_print('HydOce: ScaProdVel kin'    ,ocean_state(jg)%p_diag%kin        ,str_module,idt_src, &
!        & patch_2d%cells%owned )
!      CALL dbg_print('HydOce: ScaProdVel ptp_vn' ,ocean_state(jg)%p_diag%ptp_vn     ,str_module,idt_src, &
!        & patch_2d%edges%owned )
!      CALL dbg_print('HydOce: fu10'              ,p_as%fu10                         ,str_module,idt_src, &
!        & in_subset=patch_2d%cells%owned)
!      CALL dbg_print('HydOce: concsum'           ,sea_ice%concsum                   ,str_module,idt_src, &
!        & in_subset=patch_2d%cells%owned)
!      !------------------------------------------------------------------------
!
!      !! For updating ocean physics parameterization
!      CALL update_ho_params(patch_3d, ocean_state(jg), p_as%fu10, sea_ice%concsum, p_phys_param, &
!        & operators_coefficients, p_atm_f, p_oce_sfc)
!
!      !------------------------------------------------------------------------
!      ! solve for new free surface
!      !------------------------------------------------------------------------
!      start_timer(timer_solve_ab,1)
!      CALL solve_free_surface_eq_zstar( patch_3d, ocean_state, p_ext_data,  &
!        & p_oce_sfc , p_as, p_phys_param, operators_coefficients, solvercoeff_sp, &
!        & jstep, ocean_state(jg)%p_prog(nold(1))%eta_c, ocean_state(jg)%p_prog(nold(1))%stretch_c, &
!        & stretch_e, ocean_state(jg)%p_prog(nnew(1))%eta_c, ocean_state(jg)%p_prog(nnew(1))%stretch_c)
!
!      stop_timer(timer_solve_ab,1)
!
!      !------------------------------------------------------------------------
!      ! Step 4: calculate final normal velocity from predicted horizontal
!      ! velocity vn_pred and updated surface height
!      CALL calc_normal_velocity_ab_zstar(patch_3d, ocean_state(jg), operators_coefficients, &
!        & ocean_state(jg)%p_prog(nnew(1))%eta_c)
!     
!      !------------------------------------------------------------------------
!      ! Step 5: calculate vertical velocity and mass_flx_e from continuity equation under
!      ! incompressiblity condition in the non-shallow-water case
!      CALL calc_vert_velocity_bottomup_zstar( patch_3d, ocean_state(jg),operators_coefficients, & 
!        & ocean_state(jg)%p_prog(nnew(1))%stretch_c, stretch_e)
!      !------------------------------------------------------------------------
!   
!      CALL tracer_transport_zstar(patch_3d, ocean_state(jg), p_as, sea_ice, &
!        & p_oce_sfc, p_phys_param, operators_coefficients, current_time, &
!        & ocean_state(jg)%p_prog(nold(1))%stretch_c, stretch_e, ocean_state(jg)%p_prog(nnew(1))%stretch_c)
!      !------------------------------------------------------------------------
!     
!      !! Store in temporary variables to assign to nold
!      eta_c_new     = ocean_state(jg)%p_prog(nnew(1))%eta_c
!      stretch_c_new = ocean_state(jg)%p_prog(nnew(1))%stretch_c
!      !------------------------------------------------------------------------
! 
!      !! FIXME zstar: Diagnostics does not use zstar
!      CALL calc_fast_oce_diagnostics( patch_2d, &
!          & patch_3d, &
!          & ocean_state(1), &
!          & patch_3d%p_patch_1d(1)%dolic_c, &
!          & patch_3d%p_patch_1d(1)%prism_thick_c, &
!          & patch_3d%p_patch_1d(1)%zlev_m, &
!          & ocean_state(jg)%p_diag, &
!          & ocean_state(jg)%p_prog(nnew(1))%eta_c, &
!          & ocean_state(jg)%p_prog(nnew(1))%vn, &
!          & ocean_state(jg)%p_prog(nnew(1))%tracer, &
!          & p_atm_f, &
!          & p_oce_sfc, &
!          & sea_ice) 
!
!      !------------------------------------------------------------------------
!        
!      CALL update_statistics
!  
!      CALL output_ocean( patch_3d, &
!        & ocean_state,             &
!        & current_time,                &
!        & p_oce_sfc,          &
!        & sea_ice,               &
!        & jstep, jstep0)
!        
!      CALL reset_statistics
!
!      ! Shift time indices for the next loop
!      ! this HAS to ge into the restart files, because the start with the following loop
!      CALL update_time_indices(jg)
!  
!      ! update intermediate timestepping variables for the tracers
!      CALL update_time_g_n(ocean_state(jg))
!
!      !!ICON_OMP PARALLEL WORKSHARE
!      ocean_state(jg)%p_prog(nold(1))%eta_c(:, :)     = eta_c_new(:, :) 
!      ocean_state(jg)%p_prog(nold(1))%stretch_c(:, :) = stretch_c_new(:, :) 
!      !!ICON_OMP END PARALLEL WORKSHARE
! 
!      ! check whether time has come for writing restart file
!      IF (isCheckpoint()) THEN
!        IF (.NOT. output_mode%l_none ) THEN
!          !
!          ! For multifile restart (restart_write_mode = "joint procs multifile")
!          ! the domain flag must be set to .TRUE. in order to activate the domain,
!          ! even though we have one currently in the ocean. Without this the
!          ! processes won't write out their data into a the patch restart files.
!          !
!          patch_2d%ldom_active = .TRUE.
!          !
!          IF (i_ice_dyn == 1) CALL ice_fem_update_vel_restart(patch_2d, sea_ice) ! write FEM vel to restart or checkpoint file
!          !! FIXME zstar: Restar files will not have eta_c information
!          CALL restartDescriptor%updatePatch(patch_2d, &
!                                            &opt_nice_class=1, &
!                                            &opt_ocean_zlevels=n_zlev, &
!                                            &opt_ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:), &
!                                            &opt_ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
!          CALL restartDescriptor%writeRestart(current_time, jstep)
!        END IF
!      END IF
!
!
!      IF (isEndOfThisRun()) THEN
!        ! leave time loop
!        RETURN
!      END IF
!
!    END DO



  END SUBROUTINE test_stepping_zstar


  !-------------------------------------------------------------------------
  !>
  !! Test time stepping with z levels
  SUBROUTINE test_stepping_z( patch_3d, ocean_state, p_ext_data,  &
    & this_datetime, p_oce_sfc, p_phys_param, &
    & p_as, p_atm_f, sea_ice, &
    & hamocc_state,operators_coefficients,solvercoeff_sp)

    TYPE(t_patch_3d ), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: p_ext_data(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: p_phys_param 
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: p_atm_f
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_hamocc_state),     INTENT(inout)          :: hamocc_state
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(t_solvercoeff_singleprecision), INTENT(inout) :: solvercoeff_sp
    

    ! local variables
    INTEGER :: jstep, jg
    INTEGER :: i 
    CHARACTER(LEN=32)               :: datestring
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: jstep0 ! start counter for time loop

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_testbed_modules:test_core'

    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection
    TYPE(t_ocean_transport_state)                    :: transport_state

    TYPE(datetime), POINTER             :: current_time     => NULL()
    !! FIXME: after latest master merge huge changes in solver
    !! needs to be fixed
!    TYPE(t_ocean_solve), POINTER :: solve, solve_comp
!    CLASS(t_destructible), POINTER :: free_sfc_solver => NULL()
!    CLASS(t_destructible), POINTER :: free_sfc_solver_comp => NULL()
!    CLASS(t_destructible), POINTER :: free_sfc_solver_lhs => NULL()
!    CLASS(t_destructible), POINTER :: free_sfc_solver_trans => NULL()
!    CLASS(t_destructible), POINTER :: free_sfc_solver_comp_trans => NULL()
!    TYPE(t_surface_height_lhs_zstar), POINTER :: lhs_sh => NULL()
!
!    REAL(wp) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e)           !! 
!
!    TYPE(t_subset_range), POINTER :: owned_cells, owned_edges
!    REAL(wp), PARAMETER ::  min_top_height = 0.05_wp !we have to have at least 5cm water on topLevel of sea cells)
!    INTEGER :: n_it, n_it_sp, ret_status
!    INTEGER :: rho_switch 
!    REAL(wp) :: rn, minmaxmean(3)
!    INTEGER  :: return_status 
!    CHARACTER(LEN=12)  :: str_module = 'core_test'  ! Output of module for 1 line debug)
!
!    !------------------------------------------------------------------
!    patch_2D        => patch_3d%p_patch_2d(1)
!    rho_switch      = 1 ! 0: default 1: linear 
!
!    !------------------------------------------------------------------
!    ! no grid refinement allowed here so far
!    !------------------------------------------------------------------
!
!    IF (n_dom > 1 ) THEN
!      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
!    END IF
!
!    jg = n_dom
!    patch_2d => patch_3d%p_patch_2d(jg)
!
!    !------------------------------------------------------------------
!
!    jstep0 = 0
!    jstep  = jstep0 
!    ! local time var to be passed along, so the global is kept safe 
!    current_time => newNullDatetime()
!    !! Start time stepping
!
!    DO
!
!      ! optional memory loggin
!      CALL memory_log_add
!
!      jstep = jstep + 1
!
!      ! update model date and time mtime based
!      current_time = ocean_time_nextStep()
!
!      CALL datetimeToString(current_time, datestring)
!      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep =',jstep,'  datetime:  ', datestring
!      CALL message (TRIM(routine), message_text)
!
!      CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, solvercoeff_sp)
!
!      CALL calc_scalar_product_veloc_3d( patch_3d,  &
!        & ocean_state(jg)%p_prog(nold(1))%vn,         &
!        & ocean_state(jg)%p_diag,                     &
!        & operators_coefficients)
!
!      !! Updates velocity, tracer boundary condition
!      !! Changes height based on ice etc
!      CALL update_ocean_surface_refactor( patch_3d, ocean_state(jg), p_as, sea_ice, p_atm_f, p_oce_sfc, &
!           & current_time, operators_coefficients)
!  
!      CALL update_height_depdendent_variables( patch_3d, ocean_state(jg), p_ext_data(jg), operators_coefficients, solvercoeff_sp)
!
!      !---------------------------------------------------------------------
!
!      CALL update_ho_params(patch_3d, ocean_state(jg), p_as%fu10, sea_ice%concsum, p_phys_param, &
!        & operators_coefficients, p_atm_f, p_oce_sfc)
!
!      !---------------------------------------------------------------------
!      CALL solve_free_sfc_ab_mimetic (patch_3d, ocean_state(jg), p_ext_data(jg), &
!        & p_as, p_oce_sfc, p_phys_param, jstep, operators_coefficients, solvercoeff_sp, return_status)
!      !------------------------------------------------------------------------
!      ! Step 4: calculate final normal velocity from predicted horizontal
!      ! velocity vn_pred and updated surface height
!      CALL calc_normal_velocity_ab(patch_3d, ocean_state(jg),&
!        & operators_coefficients, solvercoeff_sp,  p_ext_data(jg), p_phys_param)
!      !------------------------------------------------------------------------
!      ! Step 5: calculate vertical velocity and mass_flx_e from continuity equation under
!      ! incompressiblity condition in the non-shallow-water case
!      CALL calc_vert_velocity( patch_3d, ocean_state(jg),operators_coefficients)
!      !------------------------------------------------------------------------
!      !------------------------------------------------------------------------
!      !Tracer transport
!      old_tracer_collection => ocean_state(jg)%p_prog(nold(1))%tracer_collection
!      new_tracer_collection => ocean_state(jg)%p_prog(nnew(1))%tracer_collection
!      !------------------------------------------------------------------------
!
!      IF (no_tracer>=1) THEN
!
!        ! fill transport_state
!        transport_state%patch_3d    => patch_3d
!        transport_state%h_old       => ocean_state(jg)%p_prog(nold(1))%h
!        transport_state%h_new       => ocean_state(jg)%p_prog(nnew(1))%h
!        transport_state%w           => ocean_state(jg)%p_diag%w  ! w_time_weighted
!        transport_state%mass_flux_e => ocean_state(jg)%p_diag%mass_flx_e
!        transport_state%vn          => ocean_state(jg)%p_diag%vn_time_weighted
!
!        ! fill boundary conditions
!        old_tracer_collection%tracer(1)%top_bc => p_oce_sfc%TopBC_Temp_vdiff
!
!        IF (no_tracer > 1) &
!          old_tracer_collection%tracer(2)%top_bc => p_oce_sfc%TopBC_Salt_vdiff
!
!        ! fill diffusion coefficients
!        DO i = 1, old_tracer_collection%no_of_tracers
!            old_tracer_collection%tracer(i)%hor_diffusion_coeff => p_phys_param%TracerDiffusion_coeff(:,:,:,i)
!            old_tracer_collection%tracer(i)%ver_diffusion_coeff => p_phys_param%a_tracer_v(:,:,:,i)
!        ENDDO
!      ENDIF
!
!      !------------------------------------------------------------------------
!      !------------------------------------------------------------------------
!      ! transport tracers and diffuse them
!
!      IF (no_tracer>=1) THEN
!          CALL advect_ocean_tracers(old_tracer_collection, new_tracer_collection, transport_state, &
!            operators_coefficients)
!      ENDIF
!
!      !------------------------------------------------------------------------
!      !------------------------------------------------------------------------
!
!      CALL output_ocean( patch_3d, &
!        & ocean_state,             &
!        & current_time,                &
!        & p_oce_sfc,          &
!        & sea_ice,               &
!        & jstep, jstep0)
!
!      ! Shift time indices for the next loop
!      ! this HAS to ge into the restart files, because the start with the following loop
!      CALL update_time_indices(jg)
!
!      ! update intermediate timestepping variables for the tracers
!      CALL update_time_g_n(ocean_state(jg))
!
!      IF (isEndOfThisRun()) THEN
!        ! leave time loop
!
!        RETURN
!      END IF
!        
!    END DO

  END SUBROUTINE test_stepping_z

END MODULE mo_ocean_testbed_zstar
