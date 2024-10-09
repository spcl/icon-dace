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

! Contains the implementation of the mathematical operators for the ocean.

MODULE mo_ocean_testbed_coriolis
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_run_config,         ONLY: test_mode
  USE mo_impl_constants,     ONLY: boundary, sea, sea_boundary !,sea,land, land_boundary, sea, max_char_length, &
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce
  USE mo_dynamics_config,    ONLY: nold
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  !USE mo_exception,          ONLY: finish, message
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_div, timer_grad
  USE mo_math_types,         ONLY: t_cartesian_coordinates
  USE mo_math_utilities,      ONLY: vector_product
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_grid_config,         ONLY: n_dom
  USE mo_ocean_math_operators

  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_physics_types, ONLY: t_ho_params
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff, no_primal_edges, no_dual_edges


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=12)           :: str_module    = 'test_coriolis  '  
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug


  PUBLIC :: test_nonlinear_coriolis_3d

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Note:at it only works in sequential mode !
  SUBROUTINE test_nonlinear_coriolis_3d(patch_3D, ocean_state, operators_coefficients)
    TYPE(t_patch_3d ),TARGET,INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(t_operator_coeff),TARGET, INTENT(in)   :: operators_coefficients

    REAL(wp) :: vort_flux(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: vort_v   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), POINTER :: vn(:,:,:)
    TYPE(t_cartesian_coordinates), POINTER  :: vn_dual(:,:,:)
    TYPE(t_patch), POINTER :: patch_2d

    patch_2d   => patch_3d%p_patch_2d(1)
    vn      => ocean_state%p_prog(nold(1))%vn
    vn_dual => ocean_state%p_diag%p_vn_dual
    
    vn = 1.0_wp
    vort_v  = 0.0_wp
    patch_2d%verts%f_v = 1.0_wp    
    
    ! vn_dual = 0.0_wp
    ! CALL rot_vertex_ocean_3d( patch_3d, vn, p_vn_dual, operators_coefficients, vort_v)
    
    write(0,*) "-------------- original coriolis budget -----------------------"
    CALL testbed_nonlinear_coriolis_3d(patch_3d, vn, vort_v, &
       & operators_coefficients, vort_flux)
    CALL diagnose_coriolis_energy(patch_3d, vn, vort_flux)

    write(0,*) "-------------- fast coriolis budget -----------------------"
    CALL testbed_nonlinear_coriolis_3d_fast(patch_3d, vn, vort_v, &
       & operators_coefficients, vort_flux)
    CALL diagnose_coriolis_energy(patch_3d, vn, vort_flux)

  END SUBROUTINE test_nonlinear_coriolis_3d
  !-------------------------------------------------------------------------
  !>
  !! Optimized version of the nonlinear_coriolis_3d
  !! Note: intel -o3 will produce very different results when using this version,
  !!  the model exhibits great sensitivity using the AtlanticBoxACC setup
  !<Optimize:inUse>
  SUBROUTINE testbed_nonlinear_coriolis_3d_fast(patch_3d, vn, vort_v, &
    & operators_coefficients, vort_flux)
    TYPE(t_patch_3d ),TARGET,INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    ! TYPE(t_cartesian_coordinates), INTENT(inout)  :: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout)                    :: vort_v   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in)          :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: startLevel! , endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, vertex_edge
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk
    REAL(wp) :: this_vort_flux(n_zlev, 2) ! for each of the two vertices

    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    startLevel    = 1
    ! endLevel    = n_zlev


!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,level,je,start_edge_index,end_edge_index, this_vort_flux, &
!ICON_OMP  vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk,  &
!ICON_OMP vertex_edge, ictr) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      ! vort_flux(:,:,blockNo) = 0.0_wp

      edge_idx_loop: DO je =  start_edge_index, end_edge_index

        this_vort_flux(:,:) = 0.0_wp

        vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
        vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
        vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
        vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)

        ! vertex 1
        ictr = 0
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex1_idx,vertex1_blk)!no_dual_cell_edges

          ictr =ictr+1

          DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

            this_vort_flux(level, 1) =  this_vort_flux(level, 1) + &
              & vn( patch_2d%verts%edge_idx(vertex1_idx,vertex1_blk,vertex_edge), level, &
              &      patch_2d%verts%edge_blk(vertex1_idx,vertex1_blk,vertex_edge))  * &
              &  operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)

          ENDDO

        END DO ! edges of this vertex

        ! vertex 2
        ictr = no_dual_edges
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex2_idx,vertex2_blk)!no_dual_cell_edges

          ictr =ictr+1

          DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

            this_vort_flux(level, 2) =  this_vort_flux(level, 2) + &
              & vn( patch_2d%verts%edge_idx(vertex2_idx,vertex2_blk,vertex_edge), level, &
              &      patch_2d%verts%edge_blk(vertex2_idx,vertex2_blk,vertex_edge))  * &
              &  operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)

          ENDDO

        END DO ! edges of this vertex

        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          vort_flux(je,level,blockNo) =  &
            & this_vort_flux(level,1) * &
            & (vort_v(vertex1_idx, level, vertex1_blk) + patch_2d%verts%f_v(vertex1_idx, vertex1_blk))  &
            & + this_vort_flux(level,2) * &
            & (vort_v(vertex2_idx, level, vertex2_blk) + patch_2d%verts%f_v(vertex2_idx, vertex2_blk))

        ENDDO

      END DO edge_idx_loop

    END DO ! blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE testbed_nonlinear_coriolis_3d_fast
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Note: it only works in sequential mode !
  SUBROUTINE testbed_nonlinear_coriolis_3d(patch_3d, vn, vort_v, &
    & operators_coefficients, vort_flux)
    TYPE(t_patch_3d ),TARGET,INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    ! TYPE(t_cartesian_coordinates), INTENT(inout)  :: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout)                    :: vort_v   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in)          :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: startLevel, endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    REAL(wp) :: vort_global, thick_edge, thick_vert

    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    startLevel    = 1
    endLevel    = n_zlev


    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      level_loop: DO level = startLevel, endLevel

        edge_idx_loop: DO je =  start_edge_index, end_edge_index

          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN

            vort_flux(je,level,blockNo) = 0.0_wp

            DO neighbor=1,2
              IF(neighbor==1) ictr = 0
              IF(neighbor==2) ictr = no_dual_edges

              il_v = patch_2d%edges%vertex_idx(je,blockNo,neighbor)
              ib_v = patch_2d%edges%vertex_blk(je,blockNo,neighbor)

              vort_global = (vort_v(il_v,level,ib_v) + patch_2d%verts%f_v(il_v,ib_v))

              thick_vert=0.0_wp
              DO vertex_edge=1, patch_2d%verts%num_edges(il_v,ib_v)
                il_e = patch_2d%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = patch_2d%verts%edge_blk(il_v,ib_v,vertex_edge)

                thick_edge = patch_3D%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
                thick_vert = thick_vert+thick_edge/patch_2d%verts%num_edges(il_v,ib_v)
              END DO

              DO vertex_edge=1, patch_2d%verts%num_edges(il_v,ib_v)!no_dual_cell_edges

                ictr =ictr+1

                il_e = patch_2d%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = patch_2d%verts%edge_blk(il_v,ib_v,vertex_edge)

                thick_edge = patch_3D%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
                vort_flux(je,level,blockNo) =  vort_flux(je,level,blockNo)+vn(il_e,level,ib_e)*vort_global&
                  & *operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)*(thick_edge/thick_vert)
              END DO
            END DO
            
          ELSE
            vort_flux(je,level,blockNo)= 0.0_wp
          ENDIF
          
        END DO edge_idx_loop
      END DO level_loop
    END DO ! blockNo = edges_inDomain%start_blockold, edges_inDomain%end_block
  END SUBROUTINE testbed_nonlinear_coriolis_3d
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Note: it only works in sequential mode !
  SUBROUTINE diagnose_coriolis_energy(patch_3d, vn, vort_flux)
    TYPE(t_patch_3d ),TARGET,INTENT(in)        :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    
    REAL(wp) :: vort_flux_budget_perlevel(n_zlev), volume(n_zlev)
    REAL(wp) :: vort_flux_budget_perlevel_pos(n_zlev),vort_flux_budget_perlevel_neg(n_zlev)
    REAL(wp) :: in_product
    
    INTEGER :: startLevel, endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    startLevel    = 1
    endLevel    = n_zlev
    !-----------------------------------------------------------------------
    !Diagnostic of energetic neutrality
    vort_flux_budget_perlevel(1:n_zlev)=0.0_wp
    vort_flux_budget_perlevel_pos(1:n_zlev)=0.0_wp
    vort_flux_budget_perlevel_neg(1:n_zlev)=0.0_wp
    volume(1:n_zlev)=0.0_wp

    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      DO je =  start_edge_index, end_edge_index
        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          in_product = vort_flux(je,level,blockNo) * vn(je,level,blockNo)
          
          vort_flux_budget_perlevel(level)=vort_flux_budget_perlevel(level)&
            & +in_product &
            & *patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
            & *patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)

          IF(in_product >0.0_wp) THEN
            vort_flux_budget_perlevel_pos(1:n_zlev)= vort_flux_budget_perlevel_pos(1:n_zlev)&
              & + in_product &
              &   * patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
              &   * patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
          ELSEIF(in_product <0.0_wp) THEN
            vort_flux_budget_perlevel_neg(1:n_zlev)= vort_flux_budget_perlevel_neg(1:n_zlev)&
              & + in_product &
              & * patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
              & * patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
          ENDIF
          
          volume(level)=volume(level)+patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
            & *patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
          
       END DO
      END DO
    END DO ! blockNo = edges_inDomain%start_blockold, edges_inDomain%end_block

    DO level = startLevel, endLevel
      write(0,*)'global vorticity budget', level, &
        & vort_flux_budget_perlevel(level)/volume(level)
      write(0,*)'global vorticity budget pos:neg', level,&
        & vort_flux_budget_perlevel_pos(level), vort_flux_budget_perlevel_neg(level), &
        & (vort_flux_budget_perlevel_pos(level) + vort_flux_budget_perlevel_neg(level))/volume(level)
    END DO
  END SUBROUTINE diagnose_coriolis_energy
  !-------------------------------------------------------------------------


END MODULE mo_ocean_testbed_coriolis

