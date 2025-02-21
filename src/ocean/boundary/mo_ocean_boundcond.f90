! Contains the implementation of the top and bottom ocean boundary conditions
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
MODULE mo_ocean_boundcond
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: max_char_length, sea_boundary, sea, min_rlcell, min_dolic
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ocean_nml,          ONLY: iswm_oce, i_bc_veloc_top, i_bc_veloc_bot, forcing_smooth_steps, &
    & forcing_windstress_u_type, OceanReferenceDensity, forcing_windStress_weight
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_run_config,         ONLY: dtime
  USE mo_exception,          ONLY: message, finish
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_ocean_types,        ONLY: t_hydro_ocean_state
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_scalar_product,     ONLY: map_cell2edges_3D
  USE mo_ocean_surface_types,ONLY: t_ocean_surface
  USE mo_ocean_physics_types,ONLY: t_ho_params, v_params
  USE mo_math_types,         ONLY: t_cartesian_coordinates
  USE mo_math_utilities,     ONLY: gvec2cvec
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range
  USE mo_sync,               ONLY: SYNC_E, sync_patch_array
  USE mo_master_config,      ONLY: isRestart
  USE mo_fortran_tools,      ONLY: set_acc_host_or_device

  IMPLICIT NONE
  
  PRIVATE
   
  
  PUBLIC :: bot_bound_cond_horz_veloc
  PUBLIC :: top_bound_cond_horz_veloc
  PUBLIC :: VelocityBottomBoundaryCondition_onBlock
  !PUBLIC :: bot_bound_cond_vert_veloc
  !PUBLIC :: top_bound_cond_vert_veloc
  
  PUBLIC :: top_bound_cond_tracer
  PUBLIC :: bot_bound_cond_tracer

  INTEGER, PARAMETER :: top=1
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_scalar_product'
  CHARACTER(len=12)  :: str_module = 'oceBoundCond'  ! Output of module for 1 line debug
  INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug
  INTEGER            :: current_step = 0

CONTAINS

  SUBROUTINE top_bound_cond_horz_veloc( patch_3D, ocean_state, p_op_coeff, p_oce_sfc, lacc )
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN)       :: patch_3D
    TYPE(t_hydro_ocean_state), INTENT(inout)   :: ocean_state            ! ocean state variable
    TYPE(t_operator_coeff), INTENT(IN)         :: p_op_coeff
    TYPE(t_ocean_surface)                      :: p_oce_sfc       ! external data
    LOGICAL, INTENT(in), OPTIONAL              :: lacc

    LOGICAL :: lzacc
    CHARACTER(len=*), PARAMETER :: routine = modname//':top_bound_cond_horz_veloc'

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (forcing_windstress_u_type > 100 .OR. forcing_windstress_u_type == 0) THEN

      ! analytic wind
      CALL top_bound_cond_horz_veloc_onEdges( patch_3D, ocean_state, p_op_coeff)
    ELSE
      ! OMIP wind
      CALL top_bound_cond_horz_veloc_fromCells( patch_3D, ocean_state, p_op_coeff, p_oce_sfc, lacc=lzacc)
    ENDIF

  END SUBROUTINE top_bound_cond_horz_veloc
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes top boundary condition for horizontal velocity. This information
  !! is required for cell velocity, i.e. we have to prescribe values for
  !! du_c/dz and dv_c/dz.
  !!
  !! The forcing fluxes are provided by mo_ho_forcing
  !!
  !<Optimize:inUse>
  SUBROUTINE top_bound_cond_horz_veloc_onEdges( patch_3D, ocean_state, p_op_coeff) 
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3D
    TYPE(t_hydro_ocean_state), INTENT(inout)   :: ocean_state            ! ocean state variable
    TYPE(t_operator_coeff), INTENT(IN)         :: p_op_coeff
    
    !Local variables
    INTEGER :: je, jb
    INTEGER :: start_index, end_index
    REAL(wp):: z_scale(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: smooth_coeff, u_diff, v_diff
    TYPE(t_subset_range), POINTER :: all_edges   
    TYPE(t_patch), POINTER        :: patch_2D
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !& routine = ('mo_ocean_boundcond:top_bound_cond_veloc')
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_edges  => patch_2D%edges%all
    !-----------------------------------------------------------------------
    IF (isRestart() .OR. i_bc_veloc_top /= 4) THEN
      smooth_coeff = 1.0_wp
    ELSE
      smooth_coeff = MIN(REAL(current_step, wp) / REAL(forcing_smooth_steps, wp), 1.0_wp)
      current_step = current_step + 1
    ENDIF
    
    ! Modification of surface wind forcing according to surface boundary condition
!ICON_OMP_PARALLEL
    IF(iswm_oce == 1)THEN
!ICON_OMP_DO PRIVATE(start_index, end_index, je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
          DO je = start_index, end_index
            z_scale(je,jb) = 1.0_wp / (OceanReferenceDensity*ocean_state%p_diag%thick_e(je,jb))
        ENDDO
      ENDDO
!ICON_OMP_END_DO 
    ELSEIF(iswm_oce /= 1)THEN
!ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_edges%start_block, all_edges%end_block
        z_scale(:,jb) = 1.0_wp / OceanReferenceDensity
      ENDDO
!ICON_OMP_END_DO
    ENDIF

    SELECT CASE (i_bc_veloc_top)

   CASE (0)
     ! this is the same as forcing_windstress_u/v_type = 0
!ICON_OMP_DO PRIVATE(start_index, end_index, je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        DO je = start_index, end_index
          ocean_state%p_aux%bc_top_vn(je,jb) = 0.0_wp
        END DO
      END DO
!ICON_OMP_END_DO


    CASE (1,4)
!ICON_OMP_DO PRIVATE(start_index, end_index, je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        DO je = start_index, end_index
!           stress_coeff = smooth_coeff * z_scale(je,jb)
          ocean_state%p_aux%bc_top_vn(je,jb) = &
            & ocean_state%p_aux%bc_top_WindStress(je,jb) * smooth_coeff * z_scale(je,jb)  * forcing_windStress_weight
        END DO
      END DO
!ICON_OMP_END_DO

    CASE default
      CALL finish("top_bound_cond_horz_veloc", "unknown i_bc_veloc_top")

    END SELECT
!ICON_OMP_END_PARALLEL
!     CALL map_cell2edges_3D( patch_3D, ocean_state%p_aux%bc_top_veloc_cc,ocean_state%p_aux%bc_top_vn,p_op_coeff,level=1)
    ! CALL sync_patch_array(SYNC_E, patch_3D%p_patch_2D(1), ocean_state%p_aux%bc_top_vn)

    !---------Debug Diagnostics-------------------------------------------
!     idt_src=2  ! output print level (1-5, fix)
!     CALL dbg_print('top bound.cond. u' ,ocean_state%p_aux%bc_top_u ,str_module,idt_src, in_subset=patch_2D%cells%owned)
!     CALL dbg_print('top bound.cond. v' ,ocean_state%p_aux%bc_top_v ,str_module,idt_src, in_subset=patch_2D%cells%owned)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('top bound.cond. vn',ocean_state%p_aux%bc_top_vn,str_module,idt_src, in_subset=patch_2D%edges%owned)
    !---------------------------------------------------------------------
    
  END SUBROUTINE top_bound_cond_horz_veloc_onEdges
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE top_bound_cond_horz_veloc_fromCells( patch_3D, ocean_state, p_op_coeff, p_oce_sfc, lacc)  !  , &
 !  & top_bc_u_c, top_bc_v_c, top_bc_u_cc )
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN)       :: patch_3D
    TYPE(t_hydro_ocean_state), INTENT(inout)   :: ocean_state            ! ocean state variable
    TYPE(t_operator_coeff), INTENT(IN)         :: p_op_coeff
    TYPE(t_ocean_surface)                      :: p_oce_sfc       ! external data
 !  REAL(wp)                                   :: top_bc_u_c(:,:) ! Top boundary condition
 !  REAL(wp)                                   :: top_bc_v_c(:,:) ! dim: (nproma,alloc_cell_blocks)
 !  TYPE(t_cartesian_coordinates), INTENT(inout) :: top_bc_u_cc(:,:)
    LOGICAL, INTENT(in), OPTIONAL              :: lacc

    !Local variables
    INTEGER :: jc, jb
    INTEGER :: start_index, end_index
    REAL(wp):: z_scale(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: smooth_coeff, u_diff, v_diff, stress_coeff
    LOGICAL :: lzacc
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: patch_2D
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !& routine = ('mo_ocean_boundcond:top_bound_cond_veloc')
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_cells => patch_2D%cells%all
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(z_scale) IF(lzacc)

    IF (isRestart()) THEN
      smooth_coeff = 1.0_wp
    ELSE
      smooth_coeff = MIN(REAL(current_step, wp) / REAL(forcing_smooth_steps, wp), 1.0_wp)
      current_step = current_step + 1
    ENDIF

    ! Modification of surface wind forcing according to surface boundary condition
!ICON_OMP_PARALLEL
    IF(iswm_oce == 1)THEN
      !z_scale(:,:) = v_base%del_zlev_m(1)*OceanReferenceDensity
      !z_scale(:,:) = OceanReferenceDensity*patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,1,:)
!ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        z_scale(:,jb) = 1.0_wp / (OceanReferenceDensity*ocean_state%p_diag%thick_c(:,jb))
        !$ACC END KERNELS
      END DO
      !$ACC WAIT(1)
!ICON_OMP_END_DO
    ELSEIF(iswm_oce /= 1)THEN
!ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        z_scale(:,jb) = 1.0_wp / OceanReferenceDensity
        !$ACC END KERNELS
      END DO
      !$ACC WAIT(1)
!ICON_OMP_END_DO
    ENDIF

    ! set to zero (NAG)
!    top_bc_u_c (:,:)      = 0.0_wp
!    top_bc_v_c (:,:)      = 0.0_wp
!    top_bc_u_cc(:,:)%x(1) = 0.0_wp
!    top_bc_u_cc(:,:)%x(2) = 0.0_wp
!    top_bc_u_cc(:,:)%x(3) = 0.0_wp

    SELECT CASE (i_bc_veloc_top)

    CASE (0)

   !  ! CALL message (TRIM(routine),'ZERO top velocity boundary conditions chosen')
!ICON_OMP_DO PRIVATE(start_index, end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block

       !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
       ocean_state%p_aux%bc_top_u(:,jb)          =0.0_wp
       ocean_state%p_aux%bc_top_v(:,jb)          =0.0_wp
       !$ACC END KERNELS

        CALL get_index_range(all_cells, jb, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = start_index, end_index
!          ocean_state%p_aux%bc_top_u(jc,jb)          =0.0_wp
!          ocean_state%p_aux%bc_top_v(jc,jb)          =0.0_wp
          ocean_state%p_aux%bc_top_veloc_cc(jc,jb)%x =0.0_wp
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)
!ICON_OMP_END_DO

    CASE (1) ! Forced by wind stress stored in p_oce_sfc

      ! CALL message (TRIM(routine),'(1) top velocity boundary condition: use surface wind stress')
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, stress_coeff) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = start_index, end_index
          IF (patch_3d%p_patch_1d(1)%dolic_c(jc, jb) > 0) THEN
            stress_coeff = z_scale(jc,jb)
            ocean_state%p_aux%bc_top_u(jc,jb)          = p_oce_sfc%TopBC_WindStress_u(jc,jb)    * stress_coeff
            ocean_state%p_aux%bc_top_v(jc,jb)          = p_oce_sfc%TopBC_WindStress_v(jc,jb)    * stress_coeff
            ocean_state%p_aux%bc_top_veloc_cc(jc,jb)%x = p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x * stress_coeff
          ENDIF
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)
!ICON_OMP_END_DO

 !  CASE (2) ! Forced by difference between wind velocity stored in p_oce_sfc and ocean velocity at top layer
 !    ! TODO: topBoundCond_windStress_u is not a velocity, but boundary condition for diffusion
 !    !       to subtract velocity from wind stress is unphysical - to be checked
 !    !  option disabled (#slo#, 2014-04)

 !    ! CALL message (TRIM(routine),'(2) top velocity boundary condition: use forc-u minus U(1) ')
 !    DO jb = all_cells%start_block, all_cells%end_block
 !      CALL get_index_range(all_cells, jb, start_index, end_index)
 !      DO jc = start_index, end_index
 !        IF(patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
 !        top_bc_u_c(jc,jb)    = ( p_oce_sfc%TopBC_WindStress_u(jc,jb)   &
 !          & - ocean_state%p_diag%u(jc,1,jb) ) / z_scale(jc,jb)
 !        top_bc_v_c(jc,jb)    = ( p_oce_sfc%TopBC_WindStress_v(jc,jb)   &
 !          & - ocean_state%p_diag%v(jc,1,jb) ) / z_scale(jc,jb)
 !        top_bc_u_cc(jc,jb)%x = ( p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x &
 !          & - ocean_state%p_diag%p_vn(jc,1,jb)%x)/z_scale(jc,jb)
 !        ENDIF
 !      END DO
 !    END DO

 !  CASE (3) ! Forced by difference between wind velocity stored in p_oce_sfc and ocean velocity at top layer as in 2
 !           ! but gradually increase the for forcing_smooth_steps
 !    ! TODO: topBoundCond_windStress_u is not a velocity, but boundary condition for diffusion
 !    !       to subtract velocity from wind stress is unphysical - to be checked
 !    !  option disabled (#slo#, 2014-04)

 !    IF (isRestart()) THEN
 !      smooth_coeff = 1.0_wp
 !    ELSE
 !      smooth_coeff = MIN(REAL(current_step, wp) / REAL(forcing_smooth_steps, wp), 1.0_wp)
 !      current_step = current_step + 1
 !    ENDIF

 !    ! CALL message (TRIM(routine),'(2) top velocity boundary condition: use forc-u minus U(1) ')
 !    DO jb = all_cells%start_block, all_cells%end_block
 !      CALL get_index_range(all_cells, jb, start_index, end_index)
 !      DO jc = start_index, end_index
 !        IF(patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
 !          u_diff = p_oce_sfc%TopBC_WindStress_u(jc,jb) - ocean_state%p_diag%u(jc,1,jb)
 !          v_diff = p_oce_sfc%TopBC_WindStress_v(jc,jb) - ocean_state%p_diag%v(jc,1,jb)

 !          stress_coeff = smooth_coeff / z_scale(jc,jb)

 !          top_bc_u_c(jc,jb)    = u_diff * stress_coeff
 !          top_bc_v_c(jc,jb)    = v_diff * stress_coeff
 !          top_bc_u_cc(jc,jb)%x = ( p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x &
 !            & - ocean_state%p_diag%p_vn(jc,1,jb)%x) * stress_coeff
 !        ENDIF
 !      END DO
 !    END DO

    CASE (4) ! gradually increase the for forcing_smooth_steps

    ! CALL message (TRIM(routine),'(2) top velocity boundary condition: use forc-u minus U(1) ')
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, stress_coeff) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = start_index, end_index
          IF(patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN

            stress_coeff = smooth_coeff * z_scale(jc,jb) * forcing_windStress_weight

            ocean_state%p_aux%bc_top_u(jc,jb)          = p_oce_sfc%TopBC_WindStress_u(jc,jb)    * stress_coeff
            ocean_state%p_aux%bc_top_v(jc,jb)          = p_oce_sfc%TopBC_WindStress_v(jc,jb)    * stress_coeff
            ocean_state%p_aux%bc_top_veloc_cc(jc,jb)%x = p_oce_sfc%TopBC_WindStress_cc(jc,jb)%x * stress_coeff

         ENDIF
       END DO
       !$ACC END PARALLEL LOOP
     END DO
     !$ACC WAIT(1)
!ICON_OMP_END_DO

    CASE default
      CALL finish("top_bound_cond_horz_veloc", "unknown i_bc_veloc_top")

    END SELECT
!ICON_OMP_END_PARALLEL

    CALL map_cell2edges_3D( patch_3D, ocean_state%p_aux%bc_top_veloc_cc,ocean_state%p_aux%bc_top_vn, &
      p_op_coeff,level=1,lacc=lzacc)
    ! CALL sync_patch_array(SYNC_E, patch_3D%p_patch_2D(1), ocean_state%p_aux%bc_top_vn)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('top bound.cond. u' ,ocean_state%p_aux%bc_top_u ,str_module,idt_src, in_subset=patch_2D%cells%owned)
    CALL dbg_print('top bound.cond. v' ,ocean_state%p_aux%bc_top_v ,str_module,idt_src, in_subset=patch_2D%cells%owned)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('top bound.cond. vn',ocean_state%p_aux%bc_top_vn,str_module,idt_src, in_subset=patch_2D%edges%owned)
    !---------------------------------------------------------------------

    !$ACC END DATA
  END SUBROUTINE top_bound_cond_horz_veloc_fromCells
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE VelocityBottomBoundaryCondition_onBlock(patch_3d, &
    & blockNo, start_edge_index, end_edge_index, &
    & vn_old, vn_pred, &
    & bc_bot_vn, lacc)
      
    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: patch_3D
    INTEGER, INTENT(in)                      :: blockNo, start_edge_index, end_edge_index
    REAL(wp)                                 :: vn_old(:,:), vn_pred(:,:)
    REAL(wp)                                 :: bc_bot_vn(:)
    LOGICAL, INTENT(in), OPTIONAL            :: lacc
    ! Local variables
    INTEGER :: bottom_level, je
    REAL(wp) :: norm, vn_max, vn
    LOGICAL :: lzacc
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_ocean_boundcond:VelocityBottomBoundaryCondition_onBlock')
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(v_params) IF(lzacc)

    SELECT CASE (i_bc_veloc_bot)

    CASE(0)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO je = start_edge_index, end_edge_index
        bc_bot_vn(je) = 0.0_wp
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

    CASE(1)!Bottom friction

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO je = start_edge_index, end_edge_index
        bottom_level =  patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)
        IF ( bottom_level > 0 ) THEN  ! wet points only
          vn_max = MAX(ABS(vn_old(je,bottom_level)), ABS(vn_pred(je,bottom_level)), &
            & ABS(vn_old(je,bottom_level) - vn_pred(je,bottom_level)))
          norm  = SQRT(vn_max * vn_max)
          bc_bot_vn(je) = v_params%bottom_drag_coeff * &
            & norm * vn_pred(je,bottom_level)
        ENDIF
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

    CASE(2)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO je = start_edge_index, end_edge_index
        bottom_level =  patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)
        IF ( bottom_level > 0 ) THEN  ! wet points only	
          !vn    = vn_pred(je,bottom_level)!
          vn    = vn_old(je,bottom_level)
          norm  = SQRT(vn * vn)
          bc_bot_vn(je) = v_params%bottom_drag_coeff * &
            & norm * vn_old(je,bottom_level)
        END IF
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

    CASE(3) !Bottom friction and topographic slope
      CALL message (TRIM(routine), &
        & 'TOPOGRAPHY_SLOPE bottom velocity boundary conditions not implemented yet')
      CALL finish (TRIM(routine), 'TOPOGRAPHY_SLOPE bottom velocity boundary conditions not implemented yet')
    CASE default
      CALL message (TRIM(routine),'choosen wrong bottom velocity boundary conditions')
    END SELECT

    !$ACC END DATA

  END SUBROUTINE VelocityBottomBoundaryCondition_onBlock
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes bottom boundary condition for horizontal velocity.
  !!
  !! Computes bottom boundary condition for horizontal velocity. This information
  !! is required for the cell velocity vector, i.e. we have to prescribe values for
  !! du_c/dz and dv_c/dz.
  !!
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE bot_bound_cond_horz_veloc( patch_3D, ocean_state, physics_parameters, p_op_coeff)
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3D
    TYPE(t_hydro_ocean_state), INTENT(inout) :: ocean_state            ! ocean state variable
    TYPE(t_ho_params), INTENT(in)            :: physics_parameters    ! physical parameters
    TYPE(t_operator_coeff), INTENT(IN)       :: p_op_coeff
    !REAL(wp), INTENT(in)                     :: div_coeff(:,:,:,:)
    
    ! Local variables
    INTEGER :: start_index, end_index, bottom_level, jb, je
    REAL(wp) :: norm, vn
    TYPE(t_subset_range), POINTER :: all_edges
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_ocean_boundcond:bot_bound_cond_veloc')
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    all_edges      => patch_2D%edges%all

!ICON_OMP_PARALLEL
    SELECT CASE (i_bc_veloc_bot)

    CASE(0)
!ICON_OMP_DO PRIVATE(start_index, end_index, je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        DO je = start_index, end_index
          ocean_state%p_aux%bc_bot_vn(je,jb) = 0.0_wp
        END DO
      END DO
!ICON_OMP_END_DO

!       DO jb = all_cells%start_block, all_cells%end_block
!         CALL get_index_range(all_cells, jb, start_index, end_index)
!         DO jc = start_index, end_index
!           ocean_state%p_aux%bc_bot_veloc_cc(jc,jb)%x = 0.0_wp
!         END DO
!       END DO

    CASE(1)!Bottom friction
 
!ICON_OMP_DO PRIVATE(start_index, end_index, je, bottom_level, vn, norm) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        DO je = start_index, end_index
          bottom_level =  patch_3D%p_patch_1D(1)%dolic_e(je,jb)
          IF ( bottom_level > 0 ) THEN  ! wet points only
            vn = ocean_state%p_prog(nold(1))%vn(je,bottom_level,jb)
            norm  = SQRT(vn * vn)
            ocean_state%p_aux%bc_bot_vn(je,jb) = physics_parameters%bottom_drag_coeff * &
              & norm * vn
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO

!       DO jb = all_cells%start_block, all_cells%end_block
!         CALL get_index_range(all_cells, jb, start_index, end_index)
!         DO jc = start_index, end_index
!           bottom_level =  patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
!           IF ( bottom_level > 0 ) THEN  ! wet points only
!             z_norm  = SQRT(2.0_wp * ocean_state%p_diag%kin(jc,bottom_level,jb))
!             ocean_state%p_aux%bc_bot_veloc_cc(jc,jb)%x = &
!               & physics_parameters%bottom_drag_coeff * z_norm * ocean_state%p_diag%p_vn(jc,bottom_level,jb)%x
! 
!           END IF
!         END DO
!       END DO

    CASE(2) !Bottom friction and topographic slope
      CALL message (TRIM(routine), &
        & 'TOPOGRAPHY_SLOPE bottom velocity boundary conditions not implemented yet')
      CALL finish (TRIM(routine), 'TOPOGRAPHY_SLOPE bottom velocity boundary conditions not implemented yet') 
! !       DO jb = all_cells%start_block, all_cells%end_block
! !         CALL get_index_range(all_cells, jb, start_index, end_index)
! !         DO jc = start_index, end_index
! !           
! !           !bottom_level = v_base%dolic_c(jc,jb)
! !           bottom_level = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
! !           IF ( bottom_level >= min_dolic ) THEN  ! wet points only
! !             
! !             z_norm  = SQRT(2.0_wp*ocean_state%p_diag%kin(jc,bottom_level,jb))
! !             
! !             ocean_state%p_aux%bc_bot_veloc_cc(jc,jb)%x =&
! !               & physics_parameters%bottom_drag_coeff*z_norm*ocean_state%p_diag%p_vn(jc,bottom_level,jb)%x
! !             
! !             !Only for RBF relevant: there should be an if, PK
! !             ocean_state%p_aux%bc_bot_u(jc,jb)=&
! !               & physics_parameters%bottom_drag_coeff*z_norm*ocean_state%p_diag%u(jc,bottom_level,jb)
! !             ocean_state%p_aux%bc_bot_v(jc,jb)=&
! !               & physics_parameters%bottom_drag_coeff*z_norm*ocean_state%p_diag%v(jc,bottom_level,jb)
! !             
! !           END IF
! !         END DO
! !       END DO
! !       
! !       !z_depth(:,1,:)=ocean_state%p_diag%thick_e
! !       z_depth(:,1,:)=patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_e(:,1,:)
! !       CALL div_oce_3d( z_depth, patch_2D, p_op_coeff%div_coeff, z_div_depth, opt_slev=1,opt_elev=1 )
! ! 
! !       
! !       ! LL: the whole loop seems to be doing nothing,
! !       !     no sync
! !       DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !         CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
! !         DO je = i_startidx_e, i_endidx_e
! !           
! !           !Get indices of two adjacent triangles
! !           il_c1 = patch_2D%edges%cell_idx(je,jb,1)
! !           ib_c1 = patch_2D%edges%cell_blk(je,jb,1)
! !           il_c2 = patch_2D%edges%cell_idx(je,jb,2)
! !           ib_c2 = patch_2D%edges%cell_blk(je,jb,2)
! !           
! !           !z_dolic_c1 = v_base%dolic_c(il_c1,ib_c1)
! !           !z_dolic_c2 = v_base%dolic_c(il_c2,ib_c2)
! ! 
! !           z_dolic_c1 = patch_3D%p_patch_1D(1)%dolic_c(il_c1,ib_c1)
! !           z_dolic_c2 = patch_3D%p_patch_1D(1)%dolic_c(il_c2,ib_c2)
! !           ! LL: this is not used
! !           IF(z_dolic_c1 >= min_dolic .AND. z_dolic_c2 >= min_dolic) THEN
! !             z_grad_u(je,1,jb)%x = (ocean_state%p_diag%p_vn(il_c2,z_dolic_c2,ib_c2)%x &
! !               & - ocean_state%p_diag%p_vn(il_c1,z_dolic_c1,ib_c1)%x) &
! !               & / patch_2D%edges%dual_edge_length(je,jb)
! !           ELSE
! !             z_grad_u(je,1,jb)%x = 0.0_wp
! !           ENDIF
! ! 
! !           ! LL: this is not used
! !           z_e(je,1,jb) = &
! !             & DOT_PRODUCT(patch_2D%edges%primal_cart_normal(je,jb)%x,z_grad_u(je,1,jb)%x)
! !           
! !         END DO
! !       END DO
! !       
! !       CALL map_cell2edges_2d(patch_3D, ocean_state%p_aux%bc_bot_veloc_cc, ocean_state%p_aux%bc_bot_vn,p_op_coeff)
! !       CALL sync_patch_array(SYNC_E, patch_2D, ocean_state%p_aux%bc_bot_v)
! !       
! !       !ocean_state%p_aux%bc_bot_vn(:,:) = ocean_state%p_aux%bc_bot_vn(:,:) - z_e(:,1,:)
! !       
    CASE default
      CALL message (TRIM(routine),'choosen wrong bottom velocity boundary conditions')
    END SELECT
!ICON_OMP_END_PARALLEL

!     CALL map_cell2edges_3D( patch_3D, ocean_state%p_aux%bc_bot_veloc_cc,ocean_state%p_aux%bc_bot_vn,p_op_coeff,level=1)
    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('bot bound.cond. vn'          ,ocean_state%p_aux%bc_bot_vn     ,str_module,idt_src, &
      in_subset=patch_2D%edges%owned)
    !---------------------------------------------------------------------
    
  END SUBROUTINE bot_bound_cond_horz_veloc
  !-------------------------------------------------------------------------
  
! !   !-------------------------------------------------------------------------
! !   !>
! !   !! Computes bottom boundary condition for vertical velocity.
! !   !! sbr calulates  Pu dot P (nabla H), this corresponds to
! !   !! continuous top boundary conditiopn u dot nabla H
! !   !!
! !   !!  mpi parallelized LL
! !   !!
! !   SUBROUTINE bot_bound_cond_vert_veloc( patch_2D, patch_3D, ocean_state, bot_bc_w )
! !     !
! !     TYPE(t_patch), TARGET, INTENT(in) :: patch_2D     !  patch on which computation is performed
! !     TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3D
! !     !
! !     ! Normal verlocity at edges
! !     TYPE(t_hydro_ocean_state), TARGET :: ocean_state
! !     !DR TYPE(external_data), INTENT(in) :: p_ext_data  !< external data
! !     !
! !     ! Bottom boundary condition at cells
! !     REAL(wp), INTENT(inout)           :: bot_bc_w(:,:) ! dim: (nproma,alloc_cell_blocks)
! !     !
! !     ! Local variables
! !     INTEGER :: jb, jc, je, i_dolic
! !     INTEGER :: start_index, end_index
! !     INTEGER :: i_startidx_e, i_endidx_e
! !     REAL(wp) :: z_grad_h(nproma,1,patch_2D%nblks_e)
! !     INTEGER, DIMENSION(:,:,:),POINTER :: iidx, iblk
! !     INTEGER, DIMENSION(:,:),  POINTER :: p_dolic
! !     REAL(wp), DIMENSION(:,:,:),   POINTER :: p_bathy
! !     TYPE(t_cartesian_coordinates) :: z_grad_h_cc(nproma,1,patch_2D%alloc_cell_blocks)
! ! 
! !     TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells    
! !     !-----------------------------------------------------------------------
! !     edges_in_domain => patch_2D%edges%in_domain
! !     all_cells => patch_2D%cells%all
! !         
! !     bot_bc_w(:,:) = 0.0_wp
! !     z_grad_h_cc(nproma,1,patch_2D%alloc_cell_blocks)%x = 0.0_wp
! !     
! !     iidx      => patch_2D%edges%cell_idx
! !     iblk      => patch_2D%edges%cell_blk
! !     !p_bathy   => v_base%zlev_m
! !     !p_dolic   => v_base%dolic_c
! !   
! !     p_bathy   => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c
! !     p_dolic   => patch_3D%p_patch_1D(1)%dolic_c
! ! 
! !     !----------------------------------------
! !     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
! !       DO je = i_startidx_e, i_endidx_e
! !       
! !         !i_dolic = v_base%dolic_e(je,jb)
! !         i_dolic =patch_3D%p_patch_1D(1)%dolic_e(je,jb)
! !         !IF ( v_base%lsm_e(je,i_dolic,jb) <= sea ) THEN
! !         IF(patch_3D%lsm_e(je,i_dolic,jb) <= sea_boundary)THEN
! !           z_grad_h(je,1,jb) =  &
! !             & ( p_bathy(iidx(je,jb,2),i_dolic,iblk(je,jb,2)) &
! !             & -  p_bathy(iidx(je,jb,1),i_dolic,iblk(je,jb,1))) &
! !             & * patch_2D%edges%inv_dual_edge_length(je,jb)
! !         ELSE
! !           z_grad_h(je,1,jb) =  0.0_wp
! !         ENDIF
! !         
! !       ENDDO
! !     END DO
! !     CALL sync_patch_array(SYNC_E, patch_2D, z_grad_h(:,:,:))
! !     !----------------------------------------
! !     
! !     !----------------------------------------
! ! !     CALL map_edges2cell( patch_2D, &
! ! !       & z_grad_h,&
! ! !       & z_grad_h_cc,&
! ! !       & opt_slev=1, opt_elev=1)
! ! !     CALL map_edges2cell_3D( patch_2D, &
! ! !       & z_grad_h,&
! ! !       & z_grad_h_cc,&
! ! !       & opt_slev=1, opt_elev=1)
! !     !----------------------------------------
! !     
! !     DO jb = all_cells%start_block, all_cells%end_block
! !       CALL get_index_range(all_cells, jb, start_index, end_index)
! !       DO jc = start_index, end_index
! !         !calulate  Pu dot P (nabla H), this corresponds to continuous top boundary condition u dot nabla H
! !         bot_bc_w(jc,jb) = -DOT_PRODUCT(z_grad_h_cc(jc,1,jb)%x,&
! !           & ocean_state%p_diag%p_vn(jc,1,jb)%x)
! !       END DO
! !     END DO
! !     
! !   END SUBROUTINE bot_bound_cond_vert_veloc
! !   !-------------------------------------------------------------------------
! !   
! !   !-------------------------------------------------------------------------
! !   !>
! !   !! Computes top boundary condition for vertical velocity.
! !   !! sbr calulates (h^(n+1)-h^n)/dt + Pu dot P (nabla h), this corresponds to
! !   !! continuous top boundary condition d_t h +u dot nabla h
! !   !!
! !   !!   no-mpi parallelized
! !   !!
! !   SUBROUTINE top_bound_cond_vert_veloc( patch_2D, ocean_state, top_bc_w, timestep)!, p_int )
! !     !
! !     TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
! !     TYPE(t_hydro_ocean_state), TARGET :: ocean_state
! !     REAL(wp),POINTER                  :: grad_coeff(:,:,:)
! !     REAL(wp), INTENT(inout)           :: top_bc_w(nproma,patch_2D%alloc_cell_blocks)
! !     INTEGER                           :: timestep
! !     !TYPE(t_int_state),TARGET,INTENT(in), OPTIONAL :: p_int
! !     
! !     ! Local variables
! !     INTEGER :: jb, jc
! !     INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! !     INTEGER :: rl_start, rl_end
! !     REAL(wp) :: z_grad_h(nproma,1,patch_2D%nblks_e)
! !     REAL(wp) :: z_u_times_gradh_c
! !     TYPE(t_cartesian_coordinates) :: z_grad_h_cc_vec(1:nproma,1,1:patch_2D%alloc_cell_blocks)
! !     REAL(wp) :: grad_h_u(1:nproma,1,1:patch_2D%alloc_cell_blocks)
! !     REAL(wp) :: grad_h_v(1:nproma,1,1:patch_2D%alloc_cell_blocks)
! !     ! CHARACTER(len=max_char_length), PARAMETER :: &
! !     !          & routine = ('mo_ocean_boundcond:bot_bound_cond_veloc')
! !     !-----------------------------------------------------------------------
! !     rl_start = 1
! !     rl_end = min_rlcell
! !     i_startblk = patch_2D%cells%start_blk(rl_start,1)
! !     i_endblk   = patch_2D%cells%end_blk(rl_end,1)
! !     
! !     top_bc_w(:,:) = 0.0_wp
! !  !  z_grad_h_cc_vec(nproma,1,patch_2D%alloc_cell_blocks)%x(:) = 0.0_wp
! !     
! !     DO jb = i_startblk, i_endblk
! !       CALL get_indices_c(patch_2D, jb, i_startblk, i_endblk,&
! !         & i_startidx, i_endidx, rl_start, rl_end)
! !       DO jc = i_startidx, i_endidx
! !         z_grad_h_cc_vec(nproma,1,patch_2D%alloc_cell_blocks)%x(:) = 0.0_wp
! !       END DO
! !     END DO
! !     
! !     !calculate normal derivative of new height field
! !     !CALL grad_fd_norm_oce_2d(ocean_state%p_prog(nnew(1))%h, &
! !     !  & patch_2D,                 &
! !     !  & z_grad_h(:,1,:))
! !      CALL grad_fd_norm_oce_2d_3d( ocean_state%p_prog(nnew(1))%h,&
! !                                 & patch_2D,               &
! !                                 & grad_coeff(:,1,:),     &
! !                                 & z_grad_h(:,1,:))
! ! ! CALL grad_fd_norm_oce_2d_3D( ocean_state%p_prog(nold(1))%h,     &
! ! !          &                  p_patch_horz,                  &
! ! !          &                  p_op_coeff%grad_coeff(:,1,:),  &
! ! !          &                  z_gradh_e(:,1,:))
! !     CALL sync_patch_array(SYNC_E, patch_2D, z_grad_h(:,1,:))
! !     
! !     IF(discretization_scheme==1)THEN
! ! !       CALL map_edges2cell( patch_2D,        &
! ! !         & z_grad_h,       &
! ! !         & z_grad_h_cc_vec,&
! ! !       !                         & ocean_state%p_diag%h_e,&
! ! !         & opt_slev=1,opt_elev=1 )
! !       
! ! ! !     ELSEIF(discretization_scheme==2)THEN
! ! ! !       
! ! ! !       CALL rbf_vec_interpol_cell( z_grad_h,&
! ! ! !         & patch_2D,    &
! ! ! !         & p_int,      &
! ! ! !         & grad_h_u,   &
! ! ! !         & grad_h_v,   &
! ! ! !         & opt_slev=1, opt_elev=1)
! ! ! !       
! ! ! !       DO jb = i_startblk, i_endblk
! ! ! !         CALL get_indices_c(patch_2D, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
! ! ! !           & rl_start, rl_end)
! ! ! !         DO jc = i_startidx, i_endidx
! ! ! !           CALL gvec2cvec(grad_h_u(jc,1,jb),        &
! ! ! !             & grad_h_v(jc,1,jb),        &
! ! ! !             & patch_2D%cells%center(jc,jb)%lon,&
! ! ! !             & patch_2D%cells%center(jc,jb)%lat,&
! ! ! !             & z_grad_h_cc_vec(jc,1,jb)%x(1),&
! ! ! !             & z_grad_h_cc_vec(jc,1,jb)%x(2),&
! ! ! !             & z_grad_h_cc_vec(jc,1,jb)%x(3) )
! ! ! !           ! if(jb==900)then
! ! ! !           ! write(*,*)'top w',grad_h_u(jc,1,jb),grad_h_v(jc,1,jb),&
! ! ! !           ! &z_grad_h_cc_vec(jc,1,jb)%x
! ! ! !           ! endif
! ! ! !         END DO
! ! ! !       END DO
! !     ENDIF
! !     
! !     !CALL message (TRIM(routine),'ZERO bottom velocity boundary conditions chosen')
! !     IF(timestep>1)THEN
! !       DO jb = i_startblk, i_endblk
! !         CALL get_indices_c(patch_2D, jb, i_startblk, i_endblk,&
! !           & i_startidx, i_endidx, rl_start, rl_end)
! !         DO jc = i_startidx, i_endidx
! !           !calulate  Pu dot P (nabla h), this corresponds to continuous top boundary condition u dot nabla h
! !           !z_h_u(jc,1,jb)*z_u_v(jc,1,jb) + z_h_v(jc,1,jb)*z_v_v(jc,1,jb)
! !           z_u_times_gradh_c = DOT_PRODUCT(z_grad_h_cc_vec(jc,1,jb)%x,ocean_state%p_diag%p_vn(jc,1,jb)%x)
! !           
! !           top_bc_w(jc,jb) = (ocean_state%p_prog(nnew(1))%h(jc,jb) - ocean_state%p_prog(nold(1))%h(jc,jb))/dtime&
! !             & + z_u_times_gradh_c
! !           !write(*,*)'top bc W:',jc,jb,ocean_state%p_diag%p_vn(jc,1,jb)%x
! !           !ocean_state%p_prog(nnew(1))%h(jc,jb), ocean_state%p_prog(nold(1))%h(jc,jb),&
! !           !          & z_u_times_gradh_c p_diag%p_vn(jc,1,jb)%x
! !         END DO
! !       END DO
! !     ELSE
! !       DO jb = i_startblk, i_endblk
! !         CALL get_indices_c(patch_2D, jb, i_startblk, i_endblk,&
! !           & i_startidx, i_endidx, rl_start, rl_end)
! !         DO jc = i_startidx, i_endidx
! !           !calulate  Pu dot P (nabla h), this corresponds to continuous top boundary condition u dot nabla h
! !           !z_h_u(jc,1,jb)*z_u_v(jc,1,jb) + z_h_v(jc,1,jb)*z_v_v(jc,1,jb)
! !           
! !           top_bc_w(jc,jb) = DOT_PRODUCT(z_grad_h_cc_vec(jc,1,jb)%x,ocean_state%p_diag%p_vn(jc,1,jb)%x)
! !           
! !           !write(*,*)'top bc W:',jc,jb,top_bc_w(jc,jb)!ocean_state%p_diag%p_vn(jc,1,jb)%x
! !           !ocean_state%p_prog(nnew(1))%h(jc,jb), ocean_state%p_prog(nold(1))%h(jc,jb),&
! !           !          & z_u_times_gradh_c p_diag%p_vn(jc,1,jb)%x
! !         END DO
! !       END DO
! !       
! !     ENDIF
! !   END SUBROUTINE top_bound_cond_vert_veloc
! !   !-------------------------------------------------------------------------

!-------------------------------------------------------------------------
  !>
  !! Computes top boundary condition for tracer specified by tracer_id.
  !! d C/dz
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE top_bound_cond_tracer( patch_2D, pstate_oce, tracer_id, p_oce_sfc, top_bc_tracer)
    
    TYPE(t_patch)    , TARGET, INTENT(in) :: patch_2D             ! patch on which computation is performed
    TYPE(t_hydro_ocean_state), INTENT(in) :: pstate_oce          ! ocean state variable
    INTEGER, INTENT(in)                   :: tracer_id
    TYPE(t_ocean_surface), INTENT(in)     :: p_oce_sfc
    REAL(wp), INTENT(inout)               :: top_bc_tracer(:,:,:) !Top boundary condition at cells for all tracers
    !
    !Local variables
    INTEGER :: jc, jb
    INTEGER :: start_index, end_index

    TYPE(t_subset_range), POINTER :: all_cells
    
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_ocean_boundcond:top_bound_cond_tracer')
    !-----------------------------------------------------------------------
    all_cells => patch_2D%cells%all
    
    IF (tracer_id == 1) THEN
!ICON_OMP_PARALLEL_DO  PRIVATE(start_index, end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
          top_bc_tracer(jc,jb, tracer_id) = p_oce_sfc%TopBC_Temp_vdiff(jc,jb)
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

     CALL dbg_print('top bound.cond.temp' ,top_bc_tracer(:,:,tracer_id), str_module, 2, in_subset=patch_2D%cells%owned)

    ELSE IF (tracer_id == 2) THEN
!ICON_OMP_PARALLEL_DO  PRIVATE(start_index, end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
          top_bc_tracer(jc,jb, tracer_id) = p_oce_sfc%TopBC_Salt_vdiff(jc,jb)
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO
    ELSE
    
    !default is here a homogeneous boundary condition: thisis subject to change.
!ICON_OMP_PARALLEL_DO  PRIVATE(start_index, end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
          top_bc_tracer(jc,jb, tracer_id) = 0.0_wp!p_oce_sfc%TopBC_Temp_vdiff(jc,jb)
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO    
      !CALL finish("top_bound_cond_tracer", "unknown boundary condition for tracer_id>2")
    END IF

    !---------Debug Diagnostics-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('top bound.cond.tracer' ,top_bc_tracer(:,:,tracer_id), str_module, idt_src, in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------
    
  END SUBROUTINE top_bound_cond_tracer
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes bottom boundary condition for tracer specified by tracer_id.
  !!
  !!
  !!  mpi parallelized LL (no sync required)
  !!
  SUBROUTINE bot_bound_cond_tracer( patch_2D, pstate_oce, tracer_id, bot_bc_tracer)
    
    TYPE(t_patch)    , TARGET, INTENT(in) :: patch_2D              ! patch on which computation is performed
    TYPE(t_hydro_ocean_state), INTENT(in) :: pstate_oce           ! ocean state variable
    INTEGER, INTENT(in)                   :: tracer_id
    REAL(wp), INTENT(inout)                 :: bot_bc_tracer(:,:,:) !Bottom boundary condition at cells for all tracers
    
    !Local variables
    INTEGER :: jc, jb
    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    all_cells => patch_2D%cells%all
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      DO jc = start_index, end_index
        bot_bc_tracer(jc,jb, tracer_id) = 0.0_wp
      END DO
    END DO
  END SUBROUTINE bot_bound_cond_tracer
  !-------------------------------------------------------------------------
  
END MODULE mo_ocean_boundcond
