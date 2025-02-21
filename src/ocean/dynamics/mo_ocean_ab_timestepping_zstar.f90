! Contains the implementation of the semi-implicit Adams-Bashforth timestepping
! for the ICON ocean model using the z* vertical co-ordinate.
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

MODULE mo_ocean_ab_timestepping_zstar
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
    & press_grad_type
  USE mo_run_config,                ONLY: dtime, debug_check_level, nsteps, output_mode
  USE mo_timer, ONLY: timer_start, timer_stop, timers_level, timer_extra1, &
    & timer_extra2, timer_extra3, timer_extra4, timer_ab_expl, timer_ab_rhs4sfc, timer_total

  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_physical_constants,        ONLY: grav, clw, rho_ref, Tf
  USE mo_ocean_initialization,      ONLY: is_initial_timestep
  USE mo_ocean_types, ONLY: t_hydro_ocean_state
  USE mo_ocean_time_events,   ONLY: ocean_time_nextStep, isCheckpoint, isEndOfThisRun, newNullDatetime
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,            ONLY: t_external_data
  USE mo_exception,                 ONLY: message, finish, warning, message_text
  USE mo_util_dbg_prnt,             ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_boundcond,           ONLY: VelocityBottomBoundaryCondition_onBlock, top_bound_cond_horz_veloc
  USE mo_ocean_thermodyn,           ONLY: calculate_density, &
    & calc_internal_press_grad_zstar, calc_internal_press_grad_zstar_chain
  USE mo_ocean_physics_types,       ONLY: t_ho_params
  USE mo_ocean_pp_scheme,           ONLY: ICON_PP_Edge_vnPredict_scheme
  USE mo_ocean_surface_types,       ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_scalar_product, ONLY: map_edges2edges_viacell_3d_const_z, map_edges2edges_viacell_2D_per_level, &
    & calc_scalar_product_veloc_3d, map_edges2edges_sc_zstar, map_edges2edges_3d_zstar, &
    & map_scalar_center2prismtop, map_edges2cell_3d
  USE mo_ocean_math_operators, ONLY: div_oce_3D_onTriangles_onBlock, smooth_onCells, div_oce_3d, &
    & grad_fd_norm_oce_2d_onblock, grad_fd_norm_oce_2d_3d, div_oce_3D_general_onBlock, &
    & div_oce_3D_onTriangles_onBlock, update_height_depdendent_variables
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
  USE mo_ocean_surface_refactor, ONLY: update_ocean_surface_refactor, update_atmos_fluxes
  USE mo_ocean_physics,         ONLY: update_ho_params
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
  USE mo_ocean_math_utils,       ONLY: solve_tridiag_block
  USE mo_fortran_tools,          ONLY: set_acc_host_or_device

  USE mo_ocean_physics_types,ONLY: v_params

  use ieee_arithmetic ! For ieee_is_nan
 
  !-------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_vert_velocity_bottomup_zstar
  PUBLIC :: calc_normal_velocity_ab_zstar
  PUBLIC :: solve_free_surface_eq_zstar
  PUBLIC :: update_zstar_variables

  INTERFACE velocity_diffusion_vertical_implicit_zstar



    MODULE PROCEDURE velocity_diffusion_vertical_implicit_zstar_scalar

  END INTERFACE

  LOGICAL :: eliminate_upper_diag = .true.

! solver object (free ocean surface)
  TYPE(t_ocean_solve) :: free_sfc_solver
! solver object (free ocean surface)
  TYPE(t_ocean_solve) :: free_sfc_solver_comp
! left-hand-side object (free ocean surface)
! communication infrastructure object (free ocean surface)
  TYPE(t_trivial_transfer), TARGET :: free_sfc_solver_trans_triv
  TYPE(t_subset_transfer), TARGET :: free_sfc_solver_trans_sub
    
  TYPE(t_surface_height_lhs_zstar), TARGET :: lhs_zstar
!
 
  !-------------------------------------------------------------------------
CONTAINS

 

  !! Update stretching variables based on new surface height
  SUBROUTINE update_zstar_variables( patch_3d, ocean_state, operators_coefficients, &
      & eta_c, stretch_c, stretch_e, lacc)
    
    TYPE(t_patch_3d ), POINTER, INTENT(in)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT) :: ocean_state
    TYPE(t_operator_coeff), INTENT(in)              :: operators_coefficients
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(INOUT) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp), INTENT(INOUT) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    
    TYPE(t_patch), POINTER :: patch_2d
    
    REAL(wp) :: H_c  (nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! Column depth at cell 

    REAL(wp) :: z_depth  (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: temp     (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp) :: flux_vz  (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp) :: flux_v   (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp) :: div_vz   (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: div_v    (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER  :: jb, jc, je, bt_lev, level, jk, start_level 
    INTEGER  :: start_index, end_index 

    REAL(wp) :: w_temp(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp) :: w_edg(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp) :: w_deriv(nproma, n_zlev + 1, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 

    INTEGER, DIMENSION(:,:,:), POINTER :: idx, blk
    TYPE(t_subset_range), POINTER :: all_cells, all_edges 
    INTEGER  :: id1, id2, bl1, bl2 
    INTEGER  :: edge_1_index, edge_1_block, edge_2_index, edge_2_block, edge_3_index, edge_3_block
    REAL(wp) :: st1, st2, st3 
    REAL(wp) :: dz_dt 
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    idx             => patch_3D%p_patch_2D(1)%edges%cell_idx
    blk             => patch_3D%p_patch_2D(1)%edges%cell_blk
    
    all_cells => patch_2d%cells%ALL
    all_edges => patch_2d%edges%ALL
    !------------------------------------------------------------------

    !$ACC DATA CREATE(H_c, z_depth, w_temp, w_edg, w_deriv) IF(lzacc)

!ICON_OMP_MASTER
      CALL sync_patch_array(sync_c, patch_2D, eta_c    )
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    z_depth(:,:,:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, bt_lev) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        
        bt_lev = patch_3d%p_patch_1d(1)%dolic_c(jc, jb)      
 
        !! Get physical depth of levels for later use
        !! FIXME TODO: Something is wrong here
        IF(patch_3D%lsm_c(jc, 1, jb) <= sea_boundary)THEN
          z_depth(jc, 1:bt_lev, jb) =  -1.0_wp*patch_3d%p_patch_1d(1)%depth_CellMiddle(jc, 1:bt_lev, jb) &
            & *stretch_c(jc, jb) + eta_c(jc, jb)  
        ENDIF
 
        H_c  (jc, jb)      = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, jb)
        if ( patch_3D%lsm_c(jc, 1, jb) <= sea_boundary ) THEN
          stretch_c(jc, jb)  = (H_c(jc, jb) + eta_c(jc, jb))/H_c(jc, jb) 
        else 
          stretch_c(jc, jb)  = 1.0_wp
        ENDIF
          
      END DO
      !$ACC END PARALLEL LOOP
    END DO ! blockNo
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_MASTER
      CALL sync_patch_array(sync_c, patch_2D, stretch_c)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER

!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, je, &
!ICON_OMP  id1, id2, bl1, bl2, st1, st2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO je = start_index, end_index
        
        !Get indices of two adjacent triangles
        id1 = idx(je, jb, 1)
        bl1 = blk(je, jb, 1)
        id2 = idx(je, jb, 2)
        bl2 = blk(je, jb, 2)

        !! Get coefficient for scalar by getting magnitude of vector coefficient
        st1 = DOT_PRODUCT(operators_coefficients%edge2cell_coeff_cc_t(je, 1, jb, 1)%x, &
          & operators_coefficients%edge2cell_coeff_cc_t(je, 1, jb ,1)%x)  
        st1 = SQRT(st1)
        st2 = DOT_PRODUCT(operators_coefficients%edge2cell_coeff_cc_t(je, 1, jb, 2)%x, &
          & operators_coefficients%edge2cell_coeff_cc_t(je, 1, jb, 2)%x)
        st2 = SQRT(st2)

        IF(patch_3D%lsm_e(je, 1, jb) <= sea_boundary)THEN
          !! Bringing in averaging inline with how dz_e is calculated
          stretch_e(je, jb) = 0.5*stretch_c(id1, bl1) + 0.5*stretch_c(id2, bl2)
        ELSE
          stretch_e(je, jb) = 1.0_wp
        ENDIF
        
      END DO
      !$ACC END PARALLEL LOOP
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO    


!ICON_OMP_MASTER
      CALL sync_patch_array(sync_e, patch_2D, stretch_e)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER


  !-------------------------------------------------------------------------
  !! Transform w* to w
  !! w* = (w - v.grad z)*(dz/dz*)^(-1)
  !-------------------------------------------------------------------------

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    w_edg(:,:,:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)

!ICON_OMP_DO PRIVATE(start_index, end_index, je, id1, bl1, id2, bl2, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_index, end_index)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO je = start_index, end_index

        id1=patch_2D%edges%cell_idx(je,jb,1)
        bl1=patch_2D%edges%cell_blk(je,jb,1)
        id2=patch_2D%edges%cell_idx(je,jb,2)
        bl2=patch_2D%edges%cell_blk(je,jb,2)

        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)

          w_edg(je,jk,jb)=&
            & (z_depth(id2,jk,bl2)-z_depth(id1,jk,bl1))*operators_coefficients%grad_coeff(je,jk,jb)* &
            & ocean_state%p_prog(nold(1))%vn(je, jk, jb)
        END DO
          
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    w_temp(:,:,:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, st1, st2, st3, &
!ICON_OMP edge_1_index, edge_1_block, edge_2_index, edge_2_block, edge_3_index, edge_3_block,  &
!ICON_OMP level) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        IF(patch_3D%lsm_c(jc, 1, jb) <= sea_boundary)THEN

          edge_1_index = patch_2d%cells%edge_idx(jc,jb,1)
          edge_1_block = patch_2d%cells%edge_blk(jc,jb,1)
          edge_2_index = patch_2d%cells%edge_idx(jc,jb,2)
          edge_2_block = patch_2d%cells%edge_blk(jc,jb,2)
          edge_3_index = patch_2d%cells%edge_idx(jc,jb,3)
          edge_3_block = patch_2d%cells%edge_blk(jc,jb,3)
  
          st1 = DOT_PRODUCT(operators_coefficients%edge2cell_coeff_cc(jc, 1, jb, 1)%x, &
            & operators_coefficients%edge2cell_coeff_cc(jc, 1, jb ,1)%x)  
          st2 = DOT_PRODUCT(operators_coefficients%edge2cell_coeff_cc(jc, 1, jb, 2)%x, &
            & operators_coefficients%edge2cell_coeff_cc(jc, 1, jb ,2)%x)  
          st3 = DOT_PRODUCT(operators_coefficients%edge2cell_coeff_cc(jc, 1, jb, 3)%x, &
            & operators_coefficients%edge2cell_coeff_cc(jc, 1, jb ,3)%x)  
  
          st1 = st1/( st1 + st2 + st3 )
          st2 = st2/( st1 + st2 + st3 )
          st3 = st3/( st1 + st2 + st3 )
  
          DO level = 1, MIN(patch_3D%p_patch_1D(1)%dolic_c(jc,jb), n_zlev)
            
            !! dz/dt = d(eta)/dt * (1 + z*/H )
            !! Remember z* = 0:-H, hence the negative sign
            dz_dt = ( ( ocean_state%p_prog(nnew(1))%eta_c(jc, jb) - &
            & ocean_state%p_prog(nold(1))%eta_c(jc, jb) )/dtime ) * &
            & ( 1.0_wp - patch_3d%p_patch_1d(1)%depth_CellMiddle(jc, level, jb)/H_c(jc, jb))
  
            w_temp(jc, level, jb) =  ( w_edg(edge_1_index, level, edge_1_block)*st1 + &
              & w_edg(edge_2_index, level, edge_2_block)*st2 +                        &
              & w_edg(edge_3_index, level, edge_3_block)*st3 ) + dz_dt                         
          END DO

        END IF
      END DO
      !$ACC END PARALLEL LOOP

    END DO ! blockNo = all_cells%start_block, all_cells%end_block
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    w_deriv(:,:,:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)

!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, level)
!ICON_OMP_DEFAULT_SCHEDULE     
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        DO jk = 1, patch_3D%p_patch_1d(1)%dolic_c(jc, jb)-1
          w_deriv(jc ,jk + 1, jb) &
          & = 0.5_wp*( w_temp(jc, jk, jb)    &
          & +          w_temp(jc, jk + 1, jb))              
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    !! FIXME: The below usage gives an error
    !! Maybe it is because of all_cells vs cells_in_domain?
!    CALL map_scalar_center2prismtop(patch_3d, w_temp, operators_coefficients, w_deriv)

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    ocean_state%p_diag%w_deriv(:,:,:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        IF(patch_3D%lsm_c(jc, 1, jb) <= sea_boundary)THEN
          DO jk = 1, MIN(patch_3D%p_patch_1D(1)%dolic_c(jc,jb), n_zlev)
            ocean_state%p_diag%w_deriv(jc, jk, jb) =  ocean_state%p_diag%w(jc, jk, jb) &
              & *stretch_c(jc, jb) &
              & + w_deriv(jc, jk, jb)
          END DO
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO ! blockNo
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO
    

!    write(*, *) '1', maxval(ocean_state%p_diag%w), maxval(ocean_state%p_diag%w_deriv)
!    write(*, *) '2', minval(ocean_state%p_diag%w), minval(ocean_state%p_diag%w_deriv)
!    write(*, *) '3', SUM(ocean_state%p_diag%w), SUM(ocean_state%p_diag%w_deriv)
!
!    CALL dbg_print('on entry: wstar', ocean_state%p_diag%w,'check', 1, in_subset=cells_in_domain)
!    CALL dbg_print('on entry: w    ', w_edg,'check', 1, in_subset=edges_in_domain)

    !$ACC END DATA
  END SUBROUTINE update_zstar_variables
 

   
  !! Init variables related to surface height elliptic solver
  SUBROUTINE init_free_sfc(patch_3d, ocean_state, op_coeffs, solverCoeff_sp, str_e )
    TYPE(t_patch_3d ),POINTER, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT) :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN), TARGET :: op_coeffs
    TYPE(t_solverCoeff_singlePrecision), INTENT(in), TARGET :: solverCoeff_sp
    REAL(wp), INTENT(IN), CONTIGUOUS :: str_e(:,:)
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_trivial_transfer), POINTER :: trans_triv
    TYPE(t_subset_transfer), POINTER :: trans_subs
    TYPE(t_ocean_solve), POINTER :: solve, solve_comp
    CLASS(t_lhs_agen), POINTER :: lhs
    CLASS(t_transfer), POINTER :: trans
    TYPE(t_ocean_solve_parm) :: par, par_sp
    INTEGER :: trans_mode, sol_type
    CHARACTER(len=*), PARAMETER :: method_name='mo_ocean_ab_timestepping_mimetic:init_free_sfc_ab_mimetic'

    IF (free_sfc_solver%is_init) RETURN
    patch_2D => patch_3d%p_patch_2d(1)
    
    CALL lhs_zstar%construct(patch_3d, ocean_state%p_diag%thick_e, &
      & op_coeffs, solverCoeff_sp, str_e)

!prepare init of solver
    CALL par%init(solve_precon_none, 1, 800, patch_2d%cells%in_domain%end_block, &
      & patch_2D%alloc_cell_blocks, nproma, patch_2d%cells%in_domain%end_index, &
      & solver_tolerance, use_absolute_solver_tolerance)
    par_sp%nidx = -1 ! indicates not using sp-solver
    sol_type = solve_gmres
! decide which solver type to use
    SELECT CASE(select_solver)
    CASE(select_gmres) ! GMRES
    CASE(select_gmres_r) ! GMRES-restart
      par%m = solver_max_iter_per_restart
      par%nr = solver_max_restart_iterations
    CASE(select_gmres_mp_r) ! GMRES-restart(sp+wp)
      par%nr = solver_max_restart_iterations
      par_sp = par
      par_sp%m = solver_max_iter_per_restart_sp
      par_sp%tol = REAL(solver_tolerance_sp, wp)
      par%m = solver_max_iter_per_restart
    CASE(select_cg) ! CG (Fletcher-Reeves)
      sol_type = solve_cg
    CASE(select_cg_mp) ! CG (Fletcher-Reeves, sp+wp)
      sol_type = solve_cg
      par_sp = par
      par_sp%tol = REAL(solver_tolerance_sp, wp)
    CASE(select_cgj) ! CG + Jacobi preconditioner
      sol_type = solve_cg
      par%pt = solve_precon_jac
    CASE(select_bcgs) ! BiCG-Stab
      sol_type = solve_bcgs
    CASE(select_legacy_gmres)
      sol_type = solve_legacy_gmres
      par%m = solver_max_iter_per_restart
      par_sp%nr = solver_max_restart_iterations
    CASE(select_mres)
      sol_type = solve_mres
      par%nr = (par%m + 18) / 19
      par%m = 19
    CASE DEFAULT
      CALL finish(method_name, "Unknown solver")
    END SELECT ! solver
! init lhs object
! allocate and init communication infrastructure object 
    SELECT CASE(select_transfer)
    CASE(0) ! all ocean workers are involved in solving (input is just copied to internal arrays)
      CALL free_sfc_solver_trans_triv%construct(solve_cell, patch_2D) ! solve only on a subset of workers
      CALL free_sfc_solver%construct(sol_type, par, par_sp, lhs_zstar, free_sfc_solver_trans_triv)
    CASE DEFAULT ! solve only on a subset of workers
      trans_mode = MERGE(solve_trans_compact, solve_trans_scatter, select_transfer .GT. 0)
      CALL free_sfc_solver_trans_sub%construct(solve_cell, patch_2D, ABS(select_transfer), trans_mode)
      CALL free_sfc_solver%construct(sol_type, par, par_sp, lhs_zstar, free_sfc_solver_trans_sub)
    END SELECT

  END SUBROUTINE init_free_sfc


  !-------------------------------------------------------------------------
  !>
  !! Computation of new velocity in Adams-Bashforth timestepping.
  !!
  !! Adapted to zstar
  !!
  SUBROUTINE calc_normal_velocity_ab_zstar(patch_3d,ocean_state, op_coeffs, eta, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff),INTENT(in)    :: op_coeffs
    REAL(wp), INTENT(in)                 :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL        :: lacc
    !  local variables
    INTEGER  :: start_edge_index, end_edge_index, je, jk, blockNo
    REAL(wp) :: gdt_x_ab_beta, one_minus_ab_gam
    REAL(wp) :: z_grad_h(nproma,patch_3d%p_patch_2d(1)%nblks_e), z_grad_h_block(nproma)
    LOGICAL :: lzacc
    TYPE(t_subset_range), POINTER :: edges_in_domain, owned_edges
    TYPE(t_patch), POINTER :: patch
    REAL(wp), POINTER :: vn_old(:,:,:), vn_new(:,:,:), vn_time_weighted(:,:,:), vn_pred(:,:,:), grad_coeff(:,:,:)
    !----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    !-----------------------------------------------------------------------
    patch           => patch_3d%p_patch_2d(1)
    edges_in_domain => patch%edges%in_domain
    owned_edges     => patch%edges%owned
    one_minus_ab_gam = 1.0_wp - ab_gam
    gdt_x_ab_beta =  grav * dtime * ab_beta

    vn_old => ocean_state%p_prog(nold(1))%vn
    vn_new => ocean_state%p_prog(nnew(1))%vn
    vn_time_weighted => ocean_state%p_diag%vn_time_weighted
    vn_pred => ocean_state%p_diag%vn_pred
    grad_coeff => op_coeffs%grad_coeff

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(z_grad_h_block) IF(lzacc)

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_edge_index,end_edge_index, je, jk, z_grad_h_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      ! Step 1) Compute normal derivative of new surface height
      CALL grad_fd_norm_oce_2d_onBlock(eta, patch, grad_coeff(:,1, blockNo), &
        & z_grad_h_block(:), start_edge_index, end_edge_index, blockNo, lacc=lzacc)

      ! Step 2) Calculate the new velocity from the predicted one and the new surface height

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO je = start_edge_index, end_edge_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          vn_new(je,jk,blockNo) = (vn_pred(je,jk,blockNo) &
            & - gdt_x_ab_beta * z_grad_h_block(je))

          vn_time_weighted(je,jk,blockNo) = ab_gam * vn_new(je,jk,blockNo) &
            & + one_minus_ab_gam * vn_old(je,jk,blockNo)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
    END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO

    CALL sync_patch_array_mult(sync_e, patch, 2, vn_new, vn_time_weighted)

    !$ACC END DATA
  END SUBROUTINE calc_normal_velocity_ab_zstar
  !-------------------------------------------------------------------------
 
  
  !-------------------------------------------------------------------------
  !>
  !! Computation of new vertical velocity using continuity equation
  !! Calculate diagnostic vertical velocity from horizontal velocity using the
  !! incommpressibility condition in the continuity equation.
  !! vertical velocity is integrated from bottom to topLevel
  !! vertical velocity is negative for positive divergence
  !! of horizontal velocity
  !!
  !! Adapted to zstar
  !!
  SUBROUTINE calc_vert_velocity_bottomup_zstar( patch_3d, ocean_state, op_coeffs, stretch_c, stretch_e, lacc)
    TYPE(t_patch_3d), TARGET :: patch_3d       ! patch on which computation is performed
    TYPE(t_hydro_ocean_state) :: ocean_state
    TYPE(t_operator_coeff), INTENT(in) :: op_coeffs
    REAL(wp), INTENT(IN)               :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp), INTENT(IN)               :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    LOGICAL, INTENT(in), OPTIONAL      :: lacc
    ! Local variables
    INTEGER :: jc, jk, blockNo, start_index, end_index
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells, cells_owned
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), POINTER :: vertical_velocity(:,:,:)
    REAL(wp) :: div_m_c(nproma, n_zlev)
    REAL(wp) :: deta_dt, H_c 
    INTEGER  :: bt_lev
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)
 
    !-----------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain  => patch_2D%cells%in_domain
    cells_owned      => patch_2D%cells%owned
    all_cells        => patch_2D%cells%all
    edges_in_domain  => patch_2D%edges%in_domain
    vertical_velocity=> ocean_state%p_diag%w
    ! due to nag -nan compiler-option:
    !------------------------------------------------------------------
    ! Step 1) Calculate divergence of horizontal velocity at all levels
    !------------------------------------------------------------------

    !$ACC DATA CREATE(div_m_c) IF(lzacc)

    !-------------------------------------------------------------------------------
    CALL map_edges2edges_3d_zstar( patch_3d, ocean_state%p_diag%vn_time_weighted, op_coeffs, &
      & stretch_e, ocean_state%p_diag%mass_flx_e, lacc=lzacc)
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      CALL div_oce_3D_onTriangles_onBlock(ocean_state%p_diag%mass_flx_e, patch_3D, op_coeffs%div_coeff, &
        & ocean_state%p_diag%div_mass_flx_c(:,:,blockNo), blockNo=blockNo, start_index=start_index, &
        & end_index=end_index, start_level=1, end_level=n_zlev, lacc=lzacc)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      div_m_c(:, :) = ocean_state%p_diag%div_mass_flx_c(:,:,blockNo)
      !$ACC END KERNELS
      !$ACC WAIT(1)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        !use bottom boundary condition for vertical velocity at bottom of prism
        ! this should be awlays zero
        
        bt_lev = patch_3d%p_patch_1d(1)%dolic_c(jc, blockNo)      
        H_c    = patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, blockNo)

        deta_dt = -SUM(div_m_c(jc, 1:patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)))
        vertical_velocity(jc, patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) + 1, blockNo) = 0.0_wp
        !$ACC LOOP SEQ
        DO jk = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), 1, -1
          vertical_velocity(jc,jk,blockNo) = vertical_velocity(jc,jk+1,blockNo) - &
            & ( ocean_state%p_diag%div_mass_flx_c(jc,jk,blockNo) + &
            & (1.0_wp/H_c)*                       &
            & deta_dt*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc, jk, blockNo) &
            & )/stretch_c(jc, blockNo)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
    END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO
    
    CALL sync_patch_array(sync_c,patch_2D,vertical_velocity)

    !$ACC END DATA
  END SUBROUTINE calc_vert_velocity_bottomup_zstar
  !-------------------------------------------------------------------------
 
  
  !-------------------------------------------------------------------------
  !>
  !! Computation of velocity predictor in Adams-Bashforth timestepping.
  !!
  SUBROUTINE calculate_explicit_term_zstar( patch_3d, ocean_state, p_phys_param,&
    & is_first_timestep, op_coeffs, p_as, stretch_c, eta_c, stretch_e)
    TYPE(t_patch_3d ), POINTER, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE (t_ho_params)                   :: p_phys_param
    LOGICAL,INTENT(in)                   :: is_first_timestep
    TYPE(t_operator_coeff), INTENT(IN), TARGET :: op_coeffs
    TYPE(t_atmos_for_ocean), INTENT(inout) :: p_as
    REAL(wp), INTENT(IN) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! stretch factor 
    REAL(wp), INTENT(IN) :: eta_c    (nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! stretch factor 
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    TYPE(t_subset_range), POINTER :: owned_edges, owned_cells
    
    owned_edges     => patch_3d%p_patch_2d(n_dom)%edges%owned
    owned_cells     => patch_3d%p_patch_2d(n_dom)%cells%owned

    ! STEP 1: horizontal advection
    CALL veloc_adv_horz_mimetic( patch_3d,         &
      & ocean_state%p_prog(nold(1))%vn,    &
      & ocean_state%p_prog(nold(1))%vn,    &
      & ocean_state%p_diag,                &
      & ocean_state%p_diag%veloc_adv_horz, &
      & op_coeffs)
    
    ! STEP 2: compute 3D contributions: gradient of hydrostatic pressure and vertical velocity advection
      
!    ! calculate density from EOS using temperature and salinity at timelevel n
!    CALL calculate_density_zstar( patch_3d,                              &
!     & ocean_state%p_prog(nold(1))%tracer(:,:,:,1:no_tracer), eta_c, stretch_c, &
!     & ocean_state%p_diag%rho(:,:,:) )
    !! z ignores any surface variations for density 
    !! To simplify the implementation for zstar and to keep
    !! well-balancedness we use the same subroutine for zstar
    CALL calculate_density( patch_3d,                         &
     & ocean_state%p_prog(nold(1))%tracer(:,:,:,1:no_tracer),&
     & ocean_state%p_diag%rho(:,:,:) )

   IF ( press_grad_type .EQ. 0 ) THEN
    CALL calc_internal_press_grad_zstar( patch_3d,&
       &                          ocean_state%p_diag%rho,&
       &                          ocean_state%p_diag%press_hyd,& 
       &                          ocean_state%p_aux%bc_total_top_potential, &
       &                          op_coeffs%grad_coeff,  &
       &                          stretch_c,             &
       &                          ocean_state%p_diag%press_grad)     
   ELSE IF ( press_grad_type .EQ. 1 ) THEN
    CALL calc_internal_press_grad_zstar_chain( patch_3d,&
       &                          ocean_state%p_diag%rho,&
       &                          ocean_state%p_diag%press_hyd,& 
       &                          ocean_state%p_aux%bc_total_top_potential, &
       &                          op_coeffs%grad_coeff,  &
       &                          stretch_c,             &
       &                          ocean_state%p_diag%press_grad)     
   END IF
    ! calculate vertical velocity advection
    !! All derivatives are calculated from level = 2
    !! Level=1 derivatives seem to be assumed to be 0
    !! With this assumption no changes would be required
    CALL veloc_adv_vert_mimetic(          &
      & patch_3d,                         &
      & ocean_state%p_diag,op_coeffs,     &
      & ocean_state%p_diag%veloc_adv_vert )

    ! STEP 3: compute harmonic or biharmoic laplacian diffusion of velocity.
    !         This term is discretized explicitly. Order and form of the laplacian
    !         are determined in mo_oce_diffusion according to namelist settings
    !! FIXME: Has not been ported so is only an approximation
    CALL velocity_diffusion(patch_3d,              &
      & ocean_state%p_prog(nold(1))%vn, &
      & p_phys_param,            &
      & ocean_state%p_diag,op_coeffs,  &
      & ocean_state%p_diag%laplacian_horz)
    
    CALL explicit_vn_pred_zstar( patch_3d, ocean_state, op_coeffs, p_phys_param, is_first_timestep, &
      & eta_c, stretch_e)
    
  END SUBROUTINE calculate_explicit_term_zstar

  !-------------------------------------------------------------------------
  SUBROUTINE calculate_explicit_term_g_n_onBlock_zstar( patch_3d, ocean_state, is_first_timestep, &
    & start_edge_index, end_edge_index, blockNo)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    LOGICAL,INTENT(in)                   :: is_first_timestep
    INTEGER,INTENT(in)                   :: start_edge_index, end_edge_index, blockNo
    INTEGER :: je, jk
    !-----------------------------------------------------------------------
    
    DO je = start_edge_index, end_edge_index
      DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        ocean_state%p_aux%g_n(je, jk, blockNo) = &
          & - ocean_state%p_diag%press_grad    (je, jk, blockNo)  &
          & - ocean_state%p_diag%grad          (je, jk, blockNo)  &            
          & - ocean_state%p_diag%veloc_adv_horz(je, jk, blockNo)  &
          & - ocean_state%p_diag%veloc_adv_vert(je, jk, blockNo)  &
          & + ocean_state%p_diag%laplacian_horz(je, jk, blockNo)  
      END DO
    END DO
    
    IF(is_first_timestep)THEN
      ocean_state%p_aux%g_nimd(1:nproma,1:n_zlev, blockNo) = &
        & ocean_state%p_aux%g_n(1:nproma,1:n_zlev,blockNo)
    ELSE
      DO je = start_edge_index, end_edge_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          ocean_state%p_aux%g_nimd(je, jk,blockNo)                          &
            & = (1.5_wp+ab_const) * ocean_state%p_aux%g_n(je, jk,blockNo)   &
            & - (0.5_wp+ab_const) * ocean_state%p_aux%g_nm1(je, jk,blockNo)
        END DO
      END DO
    ENDIF

  END SUBROUTINE calculate_explicit_term_g_n_onBlock_zstar
 

  !-------------------------------------------------------------------------
  SUBROUTINE calculate_explicit_vn_pred_3D_onBlock_zstar( patch_3d, ocean_state, z_gradh_e, &
    & start_edge_index, end_edge_index, blockNo, stretch_e)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    REAL(wp) :: z_gradh_e(nproma)
    INTEGER, INTENT(in) :: start_edge_index, end_edge_index, blockNo
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    INTEGER :: je, jk, bottom_level
    !-----------------------------------------------------------------------
    DO je = start_edge_index, end_edge_index
      DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        ocean_state%p_diag%vn_pred(je,jk,blockNo) = ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)  &
          & + dtime*(ocean_state%p_aux%g_nimd(je,jk,blockNo) &
          & - z_gradh_e(je))
      END DO
    END DO
        
    CALL VelocityBottomBoundaryCondition_onBlock(patch_3d, &
      & blockNo,start_edge_index, end_edge_index, &
      & ocean_state%p_prog(nold(1))%vn(:,:,blockNo), &
      & ocean_state%p_diag%vn_pred(:,:,blockNo),     &
      & ocean_state%p_aux%bc_bot_vn(:,blockNo))

    !IF surface forcing applied as topLevel boundary condition to vertical diffusion
    !The surface forcing is applied as volume forcing at rhs,
    !i.e. if it part of explicit term in momentum and tracer eqs.
    !in this case, topLevel boundary ondition of vertical Laplacians are homogeneous.
    !Below is the code that adds surface forcing to explicit term of momentum eq.
    ! The top BC is equivalent to using it in the vertical diffusion implicit solver
    DO je = start_edge_index, end_edge_index
      IF(patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)>=min_dolic) THEN
        ocean_state%p_diag%vn_pred(je,1,blockNo) =  ocean_state%p_diag%vn_pred(je,1,blockNo)      &
          & + dtime*ocean_state%p_aux%bc_top_vn(je,blockNo)                                  &
          & /(stretch_e(je, blockNo) * patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,1,blockNo) )
        bottom_level = patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        !! Using a stable version for boundary modification
        !! If the velocity becomes high
        IF ( abs( ocean_state%p_diag%vn_pred(je,bottom_level,blockNo) ) > 5. ) THEN
          ocean_state%p_diag%vn_pred(je,bottom_level,blockNo)           &
              & = ocean_state%p_diag%vn_pred(je,bottom_level,blockNo)/    &
              & ( 1.0_wp + dtime* v_params%bottom_drag_coeff*             &
              & abs(ocean_state%p_diag%vn_pred(je,bottom_level,blockNo))  &
              & /( stretch_e(je, blockNo) * patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,bottom_level,blockNo) ) )
        ELSE
          ocean_state%p_diag%vn_pred(je,bottom_level,blockNo)                                  &
            & = ocean_state%p_diag%vn_pred(je,bottom_level,blockNo)                            &
            & - dtime*ocean_state%p_aux%bc_bot_vn(je,blockNo)                                  &
            & /( stretch_e(je, blockNo) * patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,bottom_level,blockNo) )
        END IF
      ENDIF
    END DO
  END SUBROUTINE calculate_explicit_vn_pred_3D_onBlock_zstar

  
  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for horizontal velocity fields
  !!by inverting a scalar field..
  !------------------------------------------------------------------------
  SUBROUTINE velocity_diffusion_vertical_implicit_zstar_vector( &
    & patch_3d,                            &
    & velocity,                            &
    & a_v, stretch_e,                      &
    & operators_coefficients,              &
    & start_index, end_index, edge_block)

    TYPE(t_patch_3d), INTENT(IN) :: patch_3d
    REAL(wp), INTENT(INOUT) :: velocity(:,:)   ! on edges, (nproma, levels)
    REAL(wp), INTENT(IN) :: a_v(:,:)      ! on edges, (nproma, levels)
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor
    TYPE(t_operator_coeff), INTENT(IN) :: operators_coefficients
    INTEGER , INTENT(IN):: start_index, end_index, edge_block
    !
    REAL(wp) :: inv_prism_thickness, inv_prisms_center_distance, inv_prisms_center_distance_next
    REAL(wp) :: a(nproma, 1:n_zlev), b(nproma, 1:n_zlev), c(nproma, 1:n_zlev), d(nproma, 1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    REAL(wp) :: inv_str_e

    INTEGER, POINTER :: dolic_e(:)
    REAL(wp), POINTER :: inv_prism_thick_e(:,:)
    REAL(wp), POINTER :: inv_prism_center_dist_e(:,:)

    INTEGER :: max_bottom_level
    INTEGER :: bottom_level
    INTEGER :: edge_index, level

    !-----------------------------------------------------------------------

    dolic_e => patch_3d%p_patch_1d(1)%dolic_e(:,edge_block)
    inv_prism_thick_e => patch_3d%p_patch_1d(1)%inv_prism_thick_e(:,:,edge_block)
    inv_prism_center_dist_e => patch_3d%p_patch_1d(1)%inv_prism_center_dist_e(:,:,edge_block)

    max_bottom_level = MAXVAL(dolic_e(start_index:end_index))

    ! Nothing to diffuse?
    IF (max_bottom_level < 2) RETURN

    DO edge_index = start_index, end_index
      bottom_level = dolic_e(edge_index)

      IF (bottom_level < 1) THEN
        CYCLE
      ELSE IF (bottom_level == 1) THEN
        a(edge_index, 1) = 0.0_wp
        b(edge_index, 1) = 1.0_wp
        c(edge_index, 1) = 0.0_wp
        d(edge_index, 1) = velocity(edge_index, 1)
      ELSE
        inv_str_e = 1._wp/stretch_e(edge_index, edge_block)

        inv_prism_thickness = inv_str_e * inv_prism_thick_e(edge_index, 1)
        inv_prisms_center_distance = inv_str_e * inv_prism_center_dist_e(edge_index, 2)

        a(edge_index, 1) = 0.0_wp
        c(edge_index, 1) = -a_v(edge_index, 2) * inv_prism_thickness * inv_prisms_center_distance * dtime
        b(edge_index, 1) = 1.0_wp - c(edge_index, 1)
        d(edge_index, 1) = velocity(edge_index, 1)

        inv_prism_thickness = inv_str_e * inv_prism_thick_e(edge_index, bottom_level)
        inv_prisms_center_distance = inv_str_e * inv_prism_center_dist_e(edge_index, bottom_level)

        a(edge_index, bottom_level) = -a_v(edge_index,bottom_level) * inv_prism_thickness * inv_prisms_center_distance * dtime
        c(edge_index, bottom_level) = 0.0_wp
        b(edge_index, bottom_level) = 1.0_wp - a(edge_index, bottom_level)
        d(edge_index, bottom_level) = velocity(edge_index, bottom_level)
      END IF
    END DO

    DO level=2, max_bottom_level-1
      DO edge_index = start_index, end_index
        IF (level > dolic_e(edge_index) - 1) CYCLE

        inv_str_e = 1._wp/stretch_e(edge_index, edge_block)

        inv_prism_thickness = inv_str_e * inv_prism_thick_e(edge_index, level)
        inv_prisms_center_distance = inv_str_e * inv_prism_center_dist_e(edge_index, level)
        inv_prisms_center_distance_next = inv_str_e * inv_prism_center_dist_e(edge_index, level+1)

        a(edge_index, level) = - a_v(edge_index,level) * inv_prism_thickness * inv_prisms_center_distance * dtime
        c(edge_index, level) = - a_v(edge_index,level+1) * inv_prism_thickness * inv_prisms_center_distance_next * dtime
        b(edge_index, level) = 1.0_wp - a(edge_index, level) - c(edge_index, level)
        d(edge_index, level) = velocity(edge_index, level)
      END DO
    END DO

    CALL solve_tridiag_block ( &
        & a=a(start_index:end_index,:), &
        & b=b(start_index:end_index,:), &
        & c=c(start_index:end_index,:), &
        & d=d(start_index:end_index,:), &
        & x=velocity(start_index:end_index,:), &
        & n=dolic_e(start_index:end_index), &
        & eliminate_upper=eliminate_upper_diag &
      )

  END SUBROUTINE velocity_diffusion_vertical_implicit_zstar_vector
  !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for horizontal velocity fields
  !!by inverting a scalar field..
  !------------------------------------------------------------------------
  SUBROUTINE velocity_diffusion_vertical_implicit_zstar_scalar( &
    & patch_3d,                            &
    & velocity,                            &
    & a_v, stretch_e,                      &
    & operators_coefficients,              &
    & start_index, end_index, edge_block)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)              :: velocity(:,:)   ! on edges, (nproma, levels)
    REAL(wp), INTENT(inout)              :: a_v(:,:)      ! on edges, (nproma, levels)
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    TYPE(t_operator_coeff),INTENT(IN) ,TARGET :: operators_coefficients
    INTEGER , INTENT(in):: start_index, end_index, edge_block
    !
    REAL(wp) :: dt_inv
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev), diagonal_product
    REAL(wp) :: column_velocity(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    REAL(wp) :: inv_str_e 

    INTEGER :: bottom_level
    INTEGER :: edge_index, level

    !-----------------------------------------------------------------------
    dt_inv = 1.0_wp/dtime

    DO edge_index = start_index, end_index
      bottom_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)

      inv_str_e    = 1._wp/stretch_e(edge_index, edge_block)

      IF (bottom_level < 2 ) CYCLE ! nothing to diffuse

      ! Note : the inv_prism_thick_e, inv_prism_center_dist_e should be updated in calculate_thickness
      DO level=1, bottom_level
        inv_prism_thickness(level)        = inv_str_e*patch_3d%p_patch_1d(1)%inv_prism_thick_e(edge_index,level,edge_block)
        inv_prisms_center_distance(level) = inv_str_e*patch_3d%p_patch_1d(1)%inv_prism_center_dist_e(edge_index,level,edge_block)

        column_velocity(level) = velocity(edge_index,level)
      ENDDO


      !------------------------------------
      ! Fill triangular matrix
      ! b is diagonal, a is the upper diagonal, c is the lower
      !   top level
      a(1) = 0.0_wp
      c(1) = -a_v(edge_index,2) * inv_prism_thickness(1) &
        &* inv_prisms_center_distance(2)*dtime
      b(1) = 1.0_wp - c(1)
      DO level = 2, bottom_level-1
        a(level) = - a_v(edge_index,level)   * inv_prism_thickness(level) &
          &* inv_prisms_center_distance(level)*dtime
        c(level) = - a_v(edge_index,level+1) * inv_prism_thickness(level) &
          &* inv_prisms_center_distance(level+1)*dtime
        b(level) = 1.0_wp - a(level) - c(level)
      END DO
      ! bottom
      a(bottom_level) = -a_v(edge_index,bottom_level) * inv_prism_thickness(bottom_level)* &
        & inv_prisms_center_distance(bottom_level)*dtime
      b(bottom_level) = 1.0_wp - a(bottom_level)
      c(bottom_level) = 0.0_wp

      !! TODO: Check whether this makes sense for z*
      IF (eliminate_upper_diag) THEN
        ! solve the tridiagonal matrix by eliminating c (the upper diagonal) 
        DO level = bottom_level-1, 1, -1
            fact(level)=c(level)/b(level+1)
            b(level)=b(level)-a(level+1)*fact(level)
            c(level) = 0.0_wp
            column_velocity(level) = column_velocity(level) - fact(level)*column_velocity(level+1)
        ENDDO
  
        velocity(edge_index,1) = column_velocity(1)/b(1)
        DO level = 2, bottom_level
            velocity(edge_index,level) = (column_velocity(level) - &
              a(level)*  velocity(edge_index,level-1)) / b(level)    
        ENDDO
  
      ELSE
          ! solve the tridiagonal matrix by eliminating a (the lower diagonal) 
          DO level=2, bottom_level
            fact(level)=a(level)/b(level-1)
            b(level)=b(level)-c(level-1)*fact(level)
            a(level) = 0.0_wp
            column_velocity(level) = column_velocity(level) - fact(level)*column_velocity(level-1)
          ENDDO
          velocity(edge_index,bottom_level) = column_velocity(bottom_level)/b(bottom_level)
          DO level=bottom_level-1,1,-1
             velocity(edge_index,level) = (column_velocity(level) - &
              c(level) * velocity(edge_index,level+1)) / b(level)    
          ENDDO                 
      
      ENDIF



    END DO ! edge_index = start_index, end_index

  END SUBROUTINE velocity_diffusion_vertical_implicit_zstar_scalar
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE explicit_vn_pred_zstar( patch_3d, ocean_state, op_coeffs, p_phys_param, &
      & is_first_timestep, eta, stretch_e)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN) :: op_coeffs
    TYPE (t_ho_params) :: p_phys_param
    LOGICAL, INTENT(in)  :: is_first_timestep
    REAL(wp), INTENT(IN) :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    REAL(wp) :: z_gradh_e(nproma)
    TYPE(t_subset_range), POINTER :: edges_in_domain
    INTEGER :: start_edge_index, end_edge_index, blockNo
    TYPE(t_patch), POINTER :: patch_2D

    patch_2D        => patch_3d%p_patch_2d(n_dom)
    edges_in_domain => patch_3d%p_patch_2d(n_dom)%edges%in_domain
    ! STEP 4: calculate weighted gradient of surface height at previous timestep
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, z_gradh_e) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      z_gradh_e(:)  = 0.0_wp
      CALL grad_fd_norm_oce_2d_onBlock( eta, patch_2D, &
        & op_coeffs%grad_coeff(:,1,blockNo), z_gradh_e(:), start_edge_index, end_edge_index, blockNo)
      z_gradh_e(start_edge_index:end_edge_index) = &
         & (1.0_wp-ab_beta) * grav * z_gradh_e(start_edge_index:end_edge_index)
      ! STEP 5:
      CALL calculate_explicit_term_g_n_onBlock_zstar( patch_3d, ocean_state, is_first_timestep, &
        & start_edge_index, end_edge_index, blockNo)
        
      CALL calculate_explicit_vn_pred_3D_onBlock_zstar( patch_3d, ocean_state, z_gradh_e(:),    &
      & start_edge_index, end_edge_index, blockNo, stretch_e)
      ! calculate vertical friction, ie p_phys_param%a_veloc_v
      ! FIXME zstar: This is not modified for zstar
      ! requires density gradient, where is that calculated 
      IF (PPscheme_type == PPscheme_ICON_Edge_vnPredict_type) &
        CALL ICON_PP_Edge_vnPredict_scheme(patch_3d, blockNo, start_edge_index, end_edge_index, &
          & ocean_state, ocean_state%p_diag%vn_pred(:,:,blockNo))
      !In 3D case implicit vertical velocity diffusion is chosen
!      eliminate_upper_diag = .not. eliminate_upper_diag ! switch the methods 
      CALL velocity_diffusion_vertical_implicit_zstar(patch_3d, &
        & ocean_state%p_diag%vn_pred(:,:,blockNo), p_phys_param%a_veloc_v(:,:,blockNo), stretch_e, &
        & op_coeffs, start_edge_index, end_edge_index, blockNo)

    END DO
!ICON_OMP_END_PARALLEL_DO
  END SUBROUTINE explicit_vn_pred_zstar
  
  
  !-------------------------------------------------------------------------
  !>
  !!  Calculation of right-hand side of elliptic surface equation.
  !!  This is used in semi implicit timelevel stepping.
  !!
  SUBROUTINE fill_rhs4surface_eq_zstar( patch_3d, ocean_state, p_oce_sfc, &
      & op_coeffs, stretch_e, eta)
    ! Patch on which computation is performed
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_ocean_surface)                :: p_oce_sfc
    TYPE(t_operator_coeff), INTENT(IN) :: op_coeffs
    REAL(wp), INTENT(IN) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
    REAL(wp), INTENT(IN) :: eta(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: jc, blockNo, jk, je
    REAL(wp) ::inv_gdt2    !, delta_z
    REAL(wp) :: z_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: div_z_depth_int_c(nproma)
    REAL(wp) :: div_z_c(nproma,n_zlev)
    REAL(wp) :: z_vn_ab(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, owned_edges
    TYPE(t_patch), POINTER :: patch_2D

    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_3d%p_patch_2d(1)%cells%in_domain
    edges_in_domain => patch_3d%p_patch_2d(1)%edges%in_domain
    owned_edges     => patch_3d%p_patch_2d(1)%edges%owned
    inv_gdt2 = 1.0_wp / (grav*dtime*dtime)
    z_vn_ab(:,:,:edges_in_domain%start_block-1) = 0._wp

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      z_vn_ab(:,:,blockNo)  = 0.0_wp
      DO je = start_edge_index, end_edge_index
!DIR$ SIMD
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          z_vn_ab(je,jk,blockNo)=ab_gam * ocean_state%p_diag%vn_pred(je,jk,blockNo) &
            & + (1.0_wp -ab_gam) * ocean_state%p_prog(nold(1))%vn(je,jk,blockNo)
        END DO
      ENDDO
    END DO
!ICON_OMP_END_PARALLEL_DO

    z_vn_ab(:,:,edges_in_domain%end_block+1:) = 0._wp
    !
    ! calculate depth-integrated velocity z_e
    !  - edge-based and cell-based
    !  - 3d and 2d (surface)
    !-------------------------------------------------------------------------------
      
    CALL sync_patch_array(sync_e, patch_2D, z_vn_ab)
    CALL map_edges2edges_3d_zstar( patch_3d, z_vn_ab, op_coeffs, stretch_e, z_e )
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, jc, jk, div_z_depth_int_c, div_z_c) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      CALL div_oce_3D_onTriangles_onBlock( z_e, patch_3D, op_coeffs%div_coeff, div_z_c, &
        & blockNo=blockNo, start_index=start_cell_index, end_index=end_cell_index,      &
        & start_level=1, end_level=n_zlev)
      ! integrate div on columns
      div_z_depth_int_c(:) = 0.0_wp
      DO jc = start_cell_index, end_cell_index
        div_z_depth_int_c(jc) = SUM(div_z_c(jc, 1:patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)))
      END DO
      ocean_state%p_aux%p_rhs_sfc_eq(:,blockNo) = 0.0_wp
      DO jc = start_cell_index, end_cell_index
        IF (patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) > 0) THEN
          ocean_state%p_aux%p_rhs_sfc_eq(jc,blockNo) = ( ( eta(jc,blockNo) &
            & - dtime * div_z_depth_int_c(jc)) * inv_gdt2)
        ENDIF
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO
  
  END SUBROUTINE fill_rhs4surface_eq_zstar

 
   

  SUBROUTINE solve_free_surface_eq_zstar( patch_3d, ocean_state, p_ext_data,  &
    & p_oce_sfc , p_as, p_phys_param, operators_coefficients, solvercoeff_sp, &
    & timestep, eta_c, stretch_c, stretch_e, eta_c_new, stretch_c_new)

      TYPE(t_patch_3d ), POINTER, INTENT(in)             :: patch_3d
      TYPE(t_hydro_ocean_state), TARGET, INTENT(inout)   :: ocean_state(n_dom)
      TYPE(t_external_data), TARGET, INTENT(in)          :: p_ext_data(n_dom)
      TYPE(t_ocean_surface)                              :: p_oce_sfc
      TYPE(t_atmos_for_ocean),  INTENT(inout)            :: p_as
      TYPE(t_ho_params)                                  :: p_phys_param 
      !! NOTE: TARGET had to be added below because of a dangling pointer error
      !! raised by NAG compiler. Possibly when sending to init_free_sfc
      !! since it was not a pointer, the lifetime of the variable was limited
      !! to this subroutine only hence leading to dangling pointer
      TYPE(t_operator_coeff), TARGET, INTENT(inout)      :: operators_coefficients
      TYPE(t_solvercoeff_singleprecision), INTENT(inout) :: solvercoeff_sp
   
      INTEGER , INTENT(IN   ) :: timestep 
      REAL(wp), INTENT(IN   ) :: eta_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht 
      REAL(wp), INTENT(IN   ) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
      REAL(wp), INTENT(INOUT) :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor 
      REAL(wp), INTENT(INOUT) :: eta_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
      REAL(wp), INTENT(INOUT) :: stretch_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
      TYPE(t_patch), POINTER :: patch_2d
    
      TYPE(t_subset_range), POINTER :: owned_cells, owned_edges
    
      INTEGER :: n_it, n_it_sp, ret_status
    
      CHARACTER(LEN=max_char_length) :: string

      CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_timestepping_zstar:solve_free_sfc'
    
      CHARACTER(LEN=12)  :: str_module = 'zstar_surf'  ! Output of module for 1 line debug)
    
      REAL(wp) :: rn, minmaxmean(3)
      
      REAL(wp) :: temp_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) 

      !------------------------------------------------------------------------
      patch_2d        => patch_3d%p_patch_2d(1)
      owned_cells => patch_2D%cells%owned
      owned_edges => patch_2D%edges%owned

      !------------------------------------------------------------------------
      ! solve for new free surface
      !------------------------------------------------------------------------
      
      !! The RHS can be filled explicitly here, however, the LHS requires
      !! multiplying the Beta*Gamma term with (H+eta)/H
      !! The height goes in using map_edges2edges_viacell_2D 
      !! init sfc solver related objects, if necessary
      IF (.NOT.free_sfc_solver%is_init) &
        CALL init_free_sfc(patch_3d, ocean_state(1), operators_coefficients, solverCoeff_sp, stretch_e)  
      
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      CALL dbg_print('test    : h var'   ,ocean_state(1)%p_prog(nold(1))%h ,str_module, 2, in_subset=owned_cells)
      CALL dbg_print('on entry: h-new'   ,eta_c_new                         ,str_module, 2, in_subset=owned_cells)
      CALL dbg_print('on entry: vn-new'  ,ocean_state(1)%p_prog(nnew(1))%vn,str_module, 2, in_subset=owned_edges)

      ! Apply windstress
      CALL top_bound_cond_horz_veloc(patch_3d, ocean_state(1), operators_coefficients, p_oce_sfc)
      
      CALL calculate_explicit_term_zstar(patch_3d, ocean_state(1), p_phys_param, &
        & is_initial_timestep(timestep), operators_coefficients, p_as, stretch_c, eta_c, stretch_e)
      
      ! Calculate RHS of surface equation
      CALL fill_rhs4surface_eq_zstar(patch_3d, ocean_state(1), p_oce_sfc, &
        & operators_coefficients, stretch_e, eta_c)

      !! Update stretching co-efficient for LHS
      CALL lhs_zstar%update(stretch_e)

      ! Solve surface equation with solver

      !!ICON_OMP PARALLEL WORKSHARE
      free_sfc_solver%x_loc_wp(:,:) = eta_c(:, :)
      !!ICON_OMP END PARALLEL WORKSHARE

      free_sfc_solver%b_loc_wp => ocean_state(1)%p_aux%p_rhs_sfc_eq
  
      ! call solver
      CALL free_sfc_solver%solve(n_it, n_it_sp)
      rn = MERGE(free_sfc_solver%res_loc_wp(1), 0._wp, n_it .NE. 0)
      
      ! output of sum of iterations every timestep
      IF (idbg_mxmn >= 0) THEN
        IF (n_it_sp .NE. -2) THEN
          WRITE(string,'(2(a,i4),2(a,e28.20),a)') &
            & 'SUM of ocean_solve iteration(sp,wp) = (', &
            & n_it_sp - 1, ', ', n_it - 1, ') , residual = (', &
            & free_sfc_solver%res_loc_wp(1), ', ', rn, ')'
        ELSE
          WRITE(string,'(a,i4,a,e28.20)') &
            &'SUM of ocean_solve iteration =', n_it - 1, ', residual =', rn
        END IF
        CALL message('ocean_solve('//TRIM(free_sfc_solver%sol_type_name)//'): surface height',TRIM(string))
      ENDIF
 
      IF (rn > solver_tolerance) THEN
        ret_status = 2
        CALL warning(routine, "NOT YET CONVERGED !!")
        RETURN
      ENDIF
      !!ICON_OMP PARALLEL WORKSHARE
      eta_c_new(:,:) = free_sfc_solver%x_loc_wp(:,:)
      !!ICON_OMP END PARALLEL WORKSHARE
 
      minmaxmean(:) = global_minmaxmean(values=eta_c_new, in_subset=owned_cells)
      IF ( abs(minmaxmean(1)) >  300 ) CALL finish("Surface height too large!!")


      IF ( ieee_is_nan(minmaxmean(1)) ) CALL finish("Surface height too large!!")


      IF (createSolverMatrix) &
        CALL free_sfc_solver%dump_matrix(timestep)
  
      !-------- end of solver ---------------
      !---------------------------------------------------------------------
 
      !! stretch_e should not be updated so we are replacing with dummy variable
      !! If updated it will calculate the wrong vn_time_weighted
      CALL update_zstar_variables( patch_3d, ocean_state(1), operators_coefficients, &
        & eta_c_new, stretch_c_new, temp_e)

      !------------------------------------------------------------------------
      ! end solve for free surface and update stretching
      !------------------------------------------------------------------------
     
  END SUBROUTINE solve_free_surface_eq_zstar


  
END MODULE mo_ocean_ab_timestepping_zstar
