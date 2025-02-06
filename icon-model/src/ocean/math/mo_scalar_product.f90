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

! The module <i>mo_scalar_product</i>
! implements discrete scalar products which depend on the grid geometry only
! are used to formulate the primitive equations in weak form.
! The coefficients used for projections from triangular edges to cell centers
! and vice versa are contained in the ocean part of the model domain and
! are calculated in <i>mo_ocean_topo</i>.

!----------------------------
#include "omp_definitions.inc"
#include "icon_definitions.inc"
!----------------------------
MODULE mo_scalar_product
  !-------------------------------------------------------------------------
  !
  USE mo_kind,               ONLY: wp, sp
  USE mo_exception,          ONLY: finish
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: sea_boundary, sea, min_dolic
  USE mo_model_domain,       ONLY: t_patch, t_patch_3d
  USE mo_ocean_types,        ONLY: t_hydro_ocean_diag, t_solvercoeff_singleprecision
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce, fast_performance_level,l_ANTICIPATED_VORTICITY,&
    & KineticEnergy_type, KineticEnergy_onPrimalGrid, KineticEnergy_onDualGrid, NoKineticEnergy
  USE mo_math_types,         ONLY: t_cartesian_coordinates
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff, no_primal_edges, no_dual_edges, &
    & Get3DVectorTo2DLocal_array3D
  USE mo_ocean_math_operators,ONLY: rot_vertex_ocean_3d, map_edges2vert_3d
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_e, sync_v,sync_patch_array,sync_c, sync_patch_array_mult!,  & sync_idx, global_max
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_timer,               ONLY: timers_level, timer_start, timer_stop, timer_extra25, timer_extra26, timer_extra27, &
    & timer_extra28, timer_extra29
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_scalar_product'

  PUBLIC :: calc_scalar_product_veloc_3d
  PUBLIC :: nonlinear_coriolis_3d
  PUBLIC :: map_edges2vert_3d
  PUBLIC :: map_edges2cell_3d
  PUBLIC :: map_cell2edges_3d
  !PRIVATE :: map_edges2edges_viacell_3d
  PUBLIC :: map_edges2edges_viavert_3d
  PUBLIC :: map_edges2edges_viacell_3d_const_z
  PUBLIC :: map_edges2edges_3d_zstar
  PUBLIC :: map_edges2edges_sc_zstar
  PUBLIC :: map_edges2edges_viacell_2d_constZ_onTriangles_sp
  PUBLIC :: map_edges2edges_viacell_2D_zstar
  PUBLIC :: map_vec_prismtop2center_on_block
  PUBLIC :: map_scalar_prismtop2center
  PUBLIC :: map_scalar_center2prismtop
  PUBLIC :: map_edges2edges_viacell_2D_per_level
  PUBLIC :: map_edges2cell_with_height_3d
  PUBLIC :: map_scalar_prismtop2center_onBlock
  PUBLIC :: map_vector_center2prismtop_onBlock


  PUBLIC :: map_scalar_center2prismtop_GM
  PUBLIC :: map_scalar_prismtop2center_GM
  PUBLIC :: map_vec_prismtop2center_on_block_GM

!   INTERFACE map_edges2edges_viacell_3d
!     MODULE PROCEDURE map_edges2edges_viacell_3d_1lev
!     MODULE PROCEDURE map_edges2edges_viacell_3d_mlev
!   END INTERFACE

  INTERFACE nonlinear_coriolis_3d_fast
#if defined(__LVECTOR__) && !defined(__LVEC_BITID__)
    MODULE PROCEDURE nonlinear_coriolis_3d_fast_vector
#else
    MODULE PROCEDURE nonlinear_coriolis_3d_fast_scalar
#endif
  END INTERFACE

  INTERFACE map_edges2edges_viacell_3d_const_z
    MODULE PROCEDURE map_edges2edges_viacell_2D_constZ
!     MODULE PROCEDURE map_edges2edges_viacell_2D_constZs
    MODULE PROCEDURE map_edges2edges_viacell_3d_mlev_const_z
    MODULE PROCEDURE map_edges2edges_viacell_3d_mlev_constZs
  END INTERFACE


  INTERFACE map_cell2edges_3d
    MODULE PROCEDURE map_cell2edges_3d_1level
    MODULE PROCEDURE map_cell2edges_3d_mlevels
  END INTERFACE

  INTERFACE map_edges2cell_3d
    !MODULE PROCEDURE map_edges2cell_with_height_3d
    MODULE PROCEDURE map_edges2cell_no_height_3d
  END INTERFACE
  CHARACTER(LEN=16)          :: module_name    = 'mo_scalarprod'  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug
CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!
  !<Optimize:inUse>
  SUBROUTINE calc_scalar_product_veloc_3d( patch_3d, vn_e,&
    & p_diag, operators_coefficients, lacc )

    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)      :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_hydro_ocean_diag)  :: p_diag
    TYPE(t_operator_coeff), INTENT(in)    :: operators_coefficients
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: cell_index, blockNo, level
    INTEGER :: edge_idx, edge_blk, n
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    REAL(wp) :: w1, w2, w3, w4, w
    TYPE(t_subset_range), POINTER :: all_cells ! , cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    CHARACTER(*), PARAMETER :: method_name="calc_scalar_product_veloc_3d"
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    all_cells       => patch_2d%cells%all

    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
!     cells_in_domain => patch_2d%cells%in_domain
    !-----------------------------------------------------------------------
    startLevel = 1
    endLevel = n_zlev

    start_detail_timer(timer_extra25,5)
    CALL map_edges2vert_3d(patch_3d%p_patch_2d(1), vn_e, operators_coefficients%edge2vert_coeff_cc, &
      & p_diag%p_vn_dual, lacc = lzacc)
    stop_detail_timer(timer_extra25,5)

    start_detail_timer(timer_extra26,5)
    !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
    CALL map_edges2cell_3d(patch_3d, vn_e, operators_coefficients, p_diag%p_vn, lacc = lzacc) !, subset_range=cells_in_domain)
    stop_detail_timer(timer_extra26,5)

    start_detail_timer(timer_extra27,5)
    CALL map_edges2edges_viacell_3d_const_z( patch_3D,       &
                                    & vn_e,                  &
                                    & operators_coefficients,&
                                    & p_diag%ptp_vn,         &
                                    & lacc = lzacc )
    stop_detail_timer(timer_extra27,5)

    CALL sync_patch_array(sync_e, patch_2d, p_diag%ptp_vn)

    start_detail_timer(timer_extra28,5)
    !--------------------------------------------------------------
    !calculate kinetic energy
    ! First option is our default configuration
    SELECT CASE(KineticEnergy_type)
    CASE (NoKineticEnergy)
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO cell_index =  start_cell_index, end_cell_index
          DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)

            p_diag%kin(cell_index,level,blockNo) = 0.0_wp

          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
      !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    CASE (KineticEnergy_onDualGrid)
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)

!         p_diag%kin(:,:,blockNo)=0.0_wp
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO cell_index =  start_cell_index, end_cell_index
          DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)

            p_diag%kin(cell_index,level,blockNo) = 0.5_wp * &
              & DOT_PRODUCT(p_diag%p_vn(cell_index,level,blockNo)%x, p_diag%p_vn(cell_index,level,blockNo)%x)

          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
      !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    CASE (KineticEnergy_onPrimalGrid)
      IF (no_primal_edges /= 3) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level, n, &
!ICON_OMP edge_idx, edge_blk, w1, w2, w3, w4, w) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)

          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          p_diag%kin(:,:,blockNo)=0.0_wp
          !$ACC END KERNELS

          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          DO cell_index =  start_cell_index, end_cell_index
            w = 0.0_wp
            DO n=1,patch_2d%cells%max_connectivity
              edge_idx = edge_of_cell_idx(cell_index,blockNo,n)
              edge_blk = edge_of_cell_blk(cell_index,blockNo,n)
              IF (edge_idx > 0) THEN
                w1 = patch_2d%edges%area_edge(edge_idx, edge_blk)
                w  = w + w1
                DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
                  p_diag%kin(cell_index,level,blockNo) = p_diag%kin(cell_index,level,blockNo) + &
                    & vn_e(edge_idx, level, edge_blk)**2 * w1
                ENDDO
              ENDIF
            ENDDO

            DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
              p_diag%kin(cell_index,level,blockNo) = & ! this is doubled for some metaphysocal reason 0.5_wp * &
                (p_diag%kin(cell_index,level,blockNo) / w)! -   &
               ! ((p_diag%w(cell_index,level,blockNo)+p_diag%w(cell_index,level+1,blockNo))*0.5_wp)**2)
            END DO

          END DO
          !$ACC END PARALLEL LOOP
        END DO ! block
        !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO
!         CALL finish("calc_scalar_product_veloc_3d", "NONLINEAR_CORIOLIS_PRIMAL_GRID works only for triangles")
      ELSE
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level, w1, w2, w3, w4, w) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)

  !         p_diag%kin(:,:,blockNo)=0.0_wp

          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          DO cell_index =  start_cell_index, end_cell_index
            w1 =  patch_2d%edges%area_edge      &
              &      (edge_of_cell_idx(cell_index,blockNo,1),edge_of_cell_blk(cell_index,blockNo,1))
            w2 =  patch_2d%edges%area_edge      &
              &      (edge_of_cell_idx(cell_index,blockNo,2),edge_of_cell_blk(cell_index,blockNo,2))
            w3 =  patch_2d%edges%area_edge      &
              &      (edge_of_cell_idx(cell_index,blockNo,3),edge_of_cell_blk(cell_index,blockNo,3))
  !           w4 = 2.0_wp * patch_2d%cells%area(cell_index,blockNo)
            w = w1 + w2 + w3
            DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)

              p_diag%kin(cell_index,level,blockNo)= (&!p_diag%kin(cell_index,level,blockNo)&
                & vn_e(edge_of_cell_idx(cell_index,blockNo,1),level,edge_of_cell_blk(cell_index,blockNo,1))**2 &
                !*p_diag%ptp_vn                         (edge_of_cell_idx(cell_index,blockNo,1),level,edge_of_cell_blk(cell_index,blockNo,1))&
                & *  w1 &
                &+&
                & vn_e(edge_of_cell_idx(cell_index,blockNo,2),level,edge_of_cell_blk(cell_index,blockNo,2))**2 &
                !&*p_diag%ptp_vn                         (edge_of_cell_idx(cell_index,blockNo,2),level,edge_of_cell_blk(cell_index,blockNo,2))&
                & * w2 &
                &+&
                & vn_e(edge_of_cell_idx(cell_index,blockNo,3),level,edge_of_cell_blk(cell_index,blockNo,3))**2 &
                ! &*p_diag%ptp_vn                         (edge_of_cell_idx(cell_index,blockNo,3),level,edge_of_cell_blk(cell_index,blockNo,3))&
                & * w3 &
                & ) / w

!               p_diag%kin(cell_index,level,blockNo) = & ! this is doubled for some metaphysocal reason 0.5_wp * &
!                 & (p_diag%kin(cell_index,level,blockNo) - &
!                 & ((p_diag%w(cell_index,level,blockNo)+p_diag%w(cell_index,level+1,blockNo))*0.5_wp)**2)

            END DO
          END DO
          !$ACC END PARALLEL LOOP
        END DO ! block
        !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO
      ENDIF

    CASE default
      CALL finish(method_name,"unknown KineticEnergy_type")

    END SELECT

    stop_detail_timer(timer_extra28,5)
!       CALL sync_patch_array(sync_c, patch_2d,p_diag%kin)

      !convert cartesian velocity vector p_diag%p_vn(cell_index,level,blockNo)%x to geographical coordinate system
      !for output, sea-ice and coupling
    start_detail_timer(timer_extra29,5)
    CALL Get3DVectorTo2DLocal_array3D(vector=p_diag%p_vn, &
      & position_local=patch_2d%cells%center, &
      & levels=patch_3d%p_patch_1d(1)%dolic_c, &
      & subset=all_cells,                      &
      & geometry_info=patch_2d%geometry_info, &
      & x=p_diag%u, y=p_diag%v, lacc=lzacc)
    stop_detail_timer(timer_extra29,5)

!       DO cell_index =  start_cell_index, end_cell_index
!         DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
!           CALL Get3DVectorToPlanarLocal(p_diag%p_vn(cell_index,level,blockNo), &
!             & patch_2d%cells%center(cell_index,blockNo), patch_2d%geometry_info, &
!             & p_diag%u(cell_index,level,blockNo), p_diag%v(cell_index,level,blockNo) )
!         END DO
!       END DO

    !--------------------------------------------------------------

    CALL dbg_print('veloc_3d: vn_e',                   vn_e,      module_name,3, &
          patch_2D%edges%owned )
    CALL dbg_print('veloc_3d: kin energy',             p_diag%kin,module_name,3, &
          patch_2D%cells%owned )
!     CALL dbg_print('veloc_3d: w', p_diag%w,module_name,1, &
!           patch_2D%cells%owned )

   CALL dbg_print('veloc_3d: East-West  :u',p_diag%u,  module_name,1, &
         patch_2D%cells%owned )
   CALL dbg_print('veloc_3d: North-South :v',p_diag%v, module_name,1, &
         patch_2D%cells%owned )


  END SUBROUTINE calc_scalar_product_veloc_3d
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Optimized version of the nonlinear_coriolis_3d
  !! Note: intel -o3 will produce very different results when using this version,
  !!  the model exhibits great sensitivity using the AtlanticBoxACC setup
  !! Note: vn should be zero on land edges, otherwise it will give wrong result.
  !!
  !!Option for Anticipated Vorticity added. Code for this option is mostly identical to
  !! code without this option, but optimization suggests a bit of code doubling.
  !<Optimize:inUse>
  SUBROUTINE nonlinear_coriolis_3d_fast_scalar(patch_3d, vn, p_vn_dual, vort_v, &
    & operators_coefficients, vort_flux, lacc)
    TYPE(t_patch_3d ), POINTER, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(inout):: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout)                    :: vort_v   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in), TARGET :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL, INTENT(IN), OPTIONAL              :: lacc

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: startLevel! , endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, vertex_edge
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk
    INTEGER ::  edgeOfVertex_index, edgeOfVertex_block
    LOGICAL :: lzacc
    REAL(wp):: this_vort_flux(n_zlev, 2) ! for each of the two vertices
    REAL(wp):: thick_edge(n_zlev,2), thick_vert(n_zlev,2)
    REAL(wp):: numOfEdges(n_zlev,2)
!     REAL(wp):: this_edge_length

    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp) :: vort_flux_old(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    startLevel    = 1
    ! endLevel    = n_zlev

    CALL set_acc_host_or_device(lzacc, lacc)

    CALL rot_vertex_ocean_3d( patch_3d, vn, p_vn_dual, operators_coefficients, vort_v, lacc=lzacc)
    ! sync not needed here, but used for example for the Leith
    CALL sync_patch_array(SYNC_V, patch_2D, vort_v)

    IF(.NOT.l_ANTICIPATED_VORTICITY)THEN
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,level,je,start_edge_index,end_edge_index, this_vort_flux, &
!ICON_OMP  vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk,  &
!ICON_OMP vertex_edge, ictr, thick_vert, numOfEdges, edgeOfVertex_index, edgeOfVertex_block, &
!ICON_OMP thick_edge) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      ! vort_flux(:,:,blockNo) = 0.0_wp

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(numOfEdges, this_vort_flux) ASYNC(1) IF(lzacc)
      DO je =  start_edge_index, end_edge_index

        this_vort_flux(:,:) = 0.0_wp

        vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
        vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
        vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
        vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)

!         this_edge_length  = patch_2d%edges%dual_edge_length(je,blockNo)
        ! vertex 1
        ictr = 0

!         thick_vert(1:n_zlev,1)=0.0_wp
!         thick_vert(1:n_zlev,2)=0.0_wp
        numOfEdges(1:n_zlev,1)=0.0_wp
        numOfEdges(1:n_zlev,2)=0.0_wp

        !$ACC LOOP SEQ
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex1_idx,vertex1_blk)

          ictr =ictr+1
          edgeOfVertex_index = patch_2d%verts%edge_idx(vertex1_idx,vertex1_blk,vertex_edge)
          edgeOfVertex_block = patch_2d%verts%edge_blk(vertex1_idx,vertex1_blk,vertex_edge)

          DO level = startLevel, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,blockNo), &
            &                        patch_3d%p_patch_1d(1)%dolic_e(edgeOfVertex_index, edgeOfVertex_block))

            numOfEdges(level,1) = numOfEdges(level,1) + 1.0_wp
!             thick_edge(level,1) = patch_3D%p_patch_1d(1)%prism_thick_e(                  &
!               & edgeOfVertex_index, level, edgeOfVertex_block)

            ! This should be calculated in the calculate_thickness routine
!             thick_vert(level,1) = thick_vert(level,1)+thick_edge(level,1)

            this_vort_flux(level, 1) =  this_vort_flux(level, 1) +                        &
              & vn( edgeOfVertex_index, level, edgeOfVertex_block)                        &
              & * operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr) !  &
!               & * thick_edge(level,1) !* &
!               & patch_2d%edges%inv_dual_edge_length(edgeOfVertex_index, edgeOfVertex_block) &
!               & * this_edge_length

          ENDDO

        END DO ! edges of this vertex

        ! vertex 2
        ictr = no_dual_edges
        !$ACC LOOP SEQ
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex2_idx,vertex2_blk)!no_dual_cell_edges

          ictr =ictr+1
          edgeOfVertex_index = patch_2d%verts%edge_idx(vertex2_idx,vertex2_blk,vertex_edge)
          edgeOfVertex_block = patch_2d%verts%edge_blk(vertex2_idx,vertex2_blk,vertex_edge)

          DO level = startLevel, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,blockNo), &
            &                        patch_3d%p_patch_1d(1)%dolic_e(edgeOfVertex_index, edgeOfVertex_block))

            numOfEdges(level,2) = numOfEdges(level,2) + 1.0_wp
!             thick_edge(level,2) = patch_3D%p_patch_1d(1)%prism_thick_e(                  &
!               & edgeOfVertex_index, level, edgeOfVertex_block)


            ! This should be calculated in the calculate_thickness routine
!             thick_vert(level,2) = thick_vert(level,2)+thick_edge(level,2)

            this_vort_flux(level, 2) =  this_vort_flux(level, 2) +                        &
              & vn( edgeOfVertex_index, level, edgeOfVertex_block)                        &
              & * operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr) !  &
!               & * thick_edge(level,2) !* &
!               & patch_2d%edges%inv_dual_edge_length(edgeOfVertex_index, edgeOfVertex_block) &
!               & * this_edge_length

          ENDDO
        END DO ! edges of this vertex

        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          vort_flux(je,level,blockNo) =  &
            & this_vort_flux(level,1) &
            ! divide by the vertex thickness
!             &  * numOfEdges(level,1) / thick_vert(level,1)  &
            &  * (vort_v(vertex1_idx, level, vertex1_blk) + patch_2d%verts%f_v(vertex1_idx, vertex1_blk))  &
            & + this_vort_flux(level,2)  &
            ! divide by the vertex thickness
!             &  * numOfEdges(level,2) / thick_vert(level,2)  &
            &  * (vort_v(vertex2_idx, level, vertex2_blk) + patch_2d%verts%f_v(vertex2_idx, vertex2_blk))

        ENDDO

      END DO
      !$ACC END PARALLEL LOOP

    END DO ! blockNo = edges_inDomain%start_block, edges_inDomain%end_block
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

  ELSEIF(l_ANTICIPATED_VORTICITY)THEN

    !This is only for debugging needed
    vort_flux_old(:,:,:)=0.0_wp

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,level,je,start_edge_index,end_edge_index, this_vort_flux, &
!ICON_OMP  vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk,  &
!ICON_OMP vertex_edge, ictr, thick_vert, numOfEdges, edgeOfVertex_index, edgeOfVertex_block, &
!ICON_OMP thick_edge) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      ! vort_flux(:,:,blockNo) = 0.0_wp

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) &
      !$ACC   PRIVATE(numOfEdges, this_vort_flux, thick_edge, thick_vert) ASYNC(1) IF(lzacc)
      DO je =  start_edge_index, end_edge_index

        this_vort_flux(:,:) = 0.0_wp

        vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
        vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
        vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
        vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)

        ! vertex 1
        ictr = 0

        thick_vert(1:n_zlev,1)=0.0_wp
        thick_vert(1:n_zlev,2)=0.0_wp
        numOfEdges(1:n_zlev,1)=0.0_wp
        numOfEdges(1:n_zlev,2)=0.0_wp

        !$ACC LOOP SEQ
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex1_idx,vertex1_blk)

          ictr =ictr+1
          edgeOfVertex_index = patch_2d%verts%edge_idx(vertex1_idx,vertex1_blk,vertex_edge)
          edgeOfVertex_block = patch_2d%verts%edge_blk(vertex1_idx,vertex1_blk,vertex_edge)

          DO level = startLevel, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,blockNo), &
            &                        patch_3d%p_patch_1d(1)%dolic_e(edgeOfVertex_index, edgeOfVertex_block))

            numOfEdges(level,1) = numOfEdges(level,1) + 1.0_wp
            thick_edge(level,1) = patch_3D%p_patch_1d(1)%prism_thick_e(                  &
              & edgeOfVertex_index, level, edgeOfVertex_block)

            ! This should be calculated in the calculate_thickness routine
            thick_vert(level,1) = thick_vert(level,1)+thick_edge(level,1)

            this_vort_flux(level, 1) =  this_vort_flux(level, 1) +                        &
              & vn( edgeOfVertex_index, level, edgeOfVertex_block)                        &
              & * operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)   &
              & * thick_edge(level,1)

          ENDDO

        END DO ! edges of this vertex

        ! vertex 2
        ictr = no_dual_edges
        !$ACC LOOP SEQ
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex2_idx,vertex2_blk)!no_dual_cell_edges

          ictr =ictr+1
          edgeOfVertex_index = patch_2d%verts%edge_idx(vertex2_idx,vertex2_blk,vertex_edge)
          edgeOfVertex_block = patch_2d%verts%edge_blk(vertex2_idx,vertex2_blk,vertex_edge)

          DO level = startLevel, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,blockNo), &
            &                        patch_3d%p_patch_1d(1)%dolic_e(edgeOfVertex_index, edgeOfVertex_block))

            numOfEdges(level,2) = numOfEdges(level,2) + 1.0_wp
            thick_edge(level,2) = patch_3D%p_patch_1d(1)%prism_thick_e(                  &
              & edgeOfVertex_index, level, edgeOfVertex_block)


            ! This should be calculated in the calculate_thickness routine
            thick_vert(level,2) = thick_vert(level,2)+thick_edge(level,2)

            this_vort_flux(level, 2) =  this_vort_flux(level, 2) +                        &
              & vn( edgeOfVertex_index, level, edgeOfVertex_block)                        &
              & * operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)   &
              & * thick_edge(level,2)

          ENDDO
        END DO ! edges of this vertex

        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          vort_flux(je,level,blockNo) =  &
            & this_vort_flux(level,1) &
            ! divide by the vertex thickness
            &  * numOfEdges(level,1) / thick_vert(level,1)  &
            &  * (vort_v(vertex1_idx, level, vertex1_blk) + patch_2d%verts%f_v(vertex1_idx, vertex1_blk))  &
            & + this_vort_flux(level,2)  &
            ! divide by the vertex thickness
            &  * numOfEdges(level,2) / thick_vert(level,2)  &
            &  * (vort_v(vertex2_idx, level, vertex2_blk) + patch_2d%verts%f_v(vertex2_idx, vertex2_blk))

          vort_flux_old(je,level,blockNo)=vort_flux(je,level,blockNo)

          vort_flux(je,level,blockNo) =  vort_flux(je,level,blockNo)-&
               &   (this_vort_flux(level,1) * numOfEdges(level,1) / thick_vert(level,1)  &
               & + this_vort_flux(level,2)  * numOfEdges(level,2) / thick_vert(level,2)) &
               & *0.5_wp*(vort_v(vertex2_idx, level, vertex2_blk)-vort_v(vertex1_idx, level, vertex1_blk))&
                &/patch_2d%edges%primal_edge_length(je,blockNo)

!write(123,*)'details',0.5_wp*(vort_v(vertex2_idx, level, vertex2_blk)-vort_v(vertex1_idx, level, vertex1_blk))&
!&/patch_2d%edges%primal_edge_length(je,blockNo),&
!&this_vort_flux(level,1) * numOfEdges(level,1) / thick_vert(level,1)  &
!& + this_vort_flux(level,2)  * numOfEdges(level,2) / thick_vert(level,2),&
!&(this_vort_flux(level,1) * numOfEdges(level,1) / thick_vert(level,1)  &
!              & + this_vort_flux(level,2)  * numOfEdges(level,2) / thick_vert(level,2)) &
!              & *0.5_wp*(vort_v(vertex2_idx, level, vertex2_blk)-vort_v(vertex1_idx, level, vertex1_blk))&
!              &/patch_2d%edges%primal_edge_length(je,blockNo)
        ENDDO

      END DO
      !$ACC END PARALLEL LOOP

    END DO ! blockNo = edges_inDomain%start_block, edges_inDomain%end_block
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

  ENDIF


  END SUBROUTINE nonlinear_coriolis_3d_fast_scalar
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Optimized version of the nonlinear_coriolis_3d
  !! Note: intel -o3 will produce very different results when using this version,
  !!  the model exhibits great sensitivity using the AtlanticBoxACC setup
  !! Note: vn should be zero on land edges, otherwise it will give wrong result.
  !!
  !!Option for Anticipated Vorticity added. Code for this option is mostly identical to
  !! code without this option, but optimization suggests a bit of code doubling.
  !<Optimize:inUse>
  SUBROUTINE nonlinear_coriolis_3d_fast_vector (patch_3d, vn, p_vn_dual, vort_v, &
    & operators_coefficients, vort_flux, lacc)
    TYPE(t_patch_3d), POINTER, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout) :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(inout) :: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout) :: vort_v(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in), TARGET :: operators_coefficients
    REAL(wp), INTENT(inout) :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: startLevel! , endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: jv
    INTEGER :: start_index, end_index
    INTEGER :: vertex_edge
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk
    INTEGER :: edgeOfVertex_index, edgeOfVertex_block
    LOGICAL :: lzacc

    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_subset_range), POINTER :: verts_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER, POINTER :: dolic_e(:,:)

    INTEGER, PARAMETER :: MAX_DUAL_EDGES = 6

    REAL(wp) :: vert_vort(MAX_DUAL_EDGES, nproma, patch_3d%p_patch_2D(1)%nblks_v)

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    verts_in_domain => patch_2d%verts%in_domain
    dolic_e => patch_3d%p_patch_1d(1)%dolic_e
    !-----------------------------------------------------------------------
    startLevel    = 1
    ! endLevel    = n_zlev

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(vert_vort) IF(lzacc)

    CALL rot_vertex_ocean_3d(patch_3d, vn, p_vn_dual, operators_coefficients, vort_v, lacc=lzacc)
    ! sync not needed here, but used for example for the Leith
    CALL sync_patch_array(SYNC_V, patch_2D, vort_v)

    IF (.NOT. l_ANTICIPATED_VORTICITY) THEN
      DO level = startLevel, n_zlev
        DO blockNo = verts_in_domain%start_block, verts_in_domain%end_block
          CALL get_index_range (verts_in_domain, blockNo, start_index, end_index)

          !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          DO jv = start_index, end_index
            !NEC$ unroll_complete
            DO vertex_edge = 1, MAX_DUAL_EDGES
              edgeOfVertex_index = patch_2d%verts%edge_idx(jv,blockNo,vertex_edge)
              edgeOfVertex_block = patch_2d%verts%edge_blk(jv,blockNo,vertex_edge)

              IF (level <= dolic_e(edgeOfVertex_index, edgeOfVertex_block) .AND. &
                  & vertex_edge <= patch_2d%verts%num_edges(jv, blockNo)) THEN
                vert_vort(vertex_edge, jv, blockNo) = &
                    & vn(edgeOfVertex_index, level, edgeOfVertex_block) &
                    & * (vort_v(jv, level, blockNo) + patch_2d%verts%f_v(jv, blockNo))
              ELSE
                vert_vort(vertex_edge, jv, blockNo) = 0._wp
              END IF
            END DO
          END DO
          !$ACC END PARALLEL LOOP
        END DO

        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)

          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          vort_flux(start_index:end_index,level,blockNo) = 0.0_wp
          !$ACC END KERNELS

          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          DO je = start_index, end_index
            IF (level > dolic_e(je,blockNo)) CYCLE

            vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
            vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
            vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
            vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)

            !NEC$ unroll_complete
            DO vertex_edge = 1, MAX_DUAL_EDGES
              ! Trust that indices > patch_2d%verts%num_edges(vertex1_idx,vertex1_blk) are zero in edge2edge
              vort_flux(je, level, blockNo) = vort_flux(je, level, blockNo) + &
                  & vert_vort(vertex_edge, vertex1_idx, vertex1_blk) &
                  & * operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,vertex_edge) &
                  & + vert_vort(vertex_edge, vertex2_idx, vertex2_blk) &
                  & * operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,no_dual_edges + vertex_edge)
            END DO
          END DO
          !$ACC END PARALLEL LOOP
        END DO
      END DO
      !$ACC WAIT(1)

    ELSEIF(l_ANTICIPATED_VORTICITY)THEN
      CALL finish('nonlinear_coriolis_3d_fast_vector', 'l_ANTICIPATED_VORTICITY=.TRUE. not vectorized.')
    ENDIF

    !$ACC END DATA
  END SUBROUTINE nonlinear_coriolis_3d_fast_vector
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Note: vn must habve been synced before this routine
  !! the resulting vort_v is synced,
  !! vort_flux id calculated on edges in_domain
!<Optimize:inUse>
  SUBROUTINE nonlinear_coriolis_3d(patch_3d, vn, p_vn_dual, vort_v, &
    & operators_coefficients, vort_flux, lacc)
    TYPE(t_patch_3d ),POINTER,INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(inout)  :: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout)                    :: vort_v   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),TARGET,INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL, INTENT(IN), OPTIONAL              :: lacc

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: startLevel, endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    LOGICAL :: lzacc
    REAL(wp) :: vort_global, thick_edge, thick_vert
    REAL(wp) :: vort_flux_budget_perlevel(n_zlev), volume(n_zlev)
    REAL(wp) :: vort_flux_budget_perlevel_pos(n_zlev),vort_flux_budget_perlevel_neg(n_zlev)
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    CHARACTER(len=*), PARAMETER :: routine = modname//':nonlinear_coriolis_3d'

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
      IF (fast_performance_level > 10) THEN
        CALL nonlinear_coriolis_3d_fast(patch_3d, vn, p_vn_dual, vort_v, &
                  & operators_coefficients, vort_flux, lacc=lzacc)
        RETURN
      ENDIF

#ifdef _OPENACC
    IF (lzacc) CALL finish(routine, 'OpenACC version for fast_performance_level <= 10 currently not tested/validated')
#endif

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    startLevel    = 1
    endLevel    = n_zlev

    !vn=0.0_wp
    !vn(3,1:n_zlev,5)=1.0E14

    CALL rot_vertex_ocean_3d( patch_3d, vn, p_vn_dual, operators_coefficients, vort_v)
    ! sync not needed here, but used for example for the Leith
    CALL sync_patch_array(SYNC_V, patch_2D, vort_v)


    ! !$OMP PARALLEL
    ! !$OMP DO PRIVATE(blockNo,level,je,start_edge_index,end_edge_index)
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      level_loop: DO level = startLevel, endLevel

        edge_idx_loop: DO je =  start_edge_index, end_edge_index

          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN

            vort_flux(je,level,blockNo) = 0.0_wp

            !$ACC LOOP SEQ
            DO neighbor=1,2
              IF(neighbor==1) ictr = 0
              IF(neighbor==2) ictr = no_dual_edges

              il_v = patch_2d%edges%vertex_idx(je,blockNo,neighbor)
              ib_v = patch_2d%edges%vertex_blk(je,blockNo,neighbor)

              vort_global = (vort_v(il_v,level,ib_v) + patch_2d%verts%f_v(il_v,ib_v))
              !$ACC LOOP SEQ
              DO vertex_edge=1, patch_2d%verts%num_edges(il_v,ib_v)

                ictr =ictr+1

                il_e = patch_2d%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = patch_2d%verts%edge_blk(il_v,ib_v,vertex_edge)

!                 thick_edge = patch_3D%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
                vort_flux(je,level,blockNo) =  vort_flux(je,level,blockNo)+vn(il_e,level,ib_e)*vort_global&
                  & *operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)!*(thick_edge/thick_vert)
              END DO
            END DO
          ELSE
            vort_flux(je,level,blockNo)= 0.0_wp
          ENDIF
        END DO edge_idx_loop
      END DO level_loop
      !$ACC END PARALLEL LOOP
    END DO ! blockNo = edges_inDomain%start_blockold, edges_inDomain%end_block
    !$ACC WAIT(1)
    ! !$OMP END DO NOWAIT
    ! !$OMP END PARALLEL

!     !Diagnostic of energetic neutrality
!     vort_flux_budget_perlevel(1:n_zlev)=0.0_wp
!     vort_flux_budget_perlevel_pos(1:n_zlev)=0.0_wp
!     vort_flux_budget_perlevel_neg(1:n_zlev)=0.0_wp
!     ! !$OMP PARALLEL
!     ! !$OMP DO PRIVATE(blockNo,level,je,start_edge_index,end_edge_index)
!     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!       CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
!       DO level = startLevel, endLevel
!       volume(level)=0.0_wp
!
!         DO je =  start_edge_index, end_edge_index
!
!           IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN
!           vort_flux_budget_perlevel(level)=vort_flux_budget_perlevel(level)&
!           &+vort_flux(je,level,blockNo)*vn(je,level,blockNo)&
!           &*patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)!&
!           !&*patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
!
!           IF(vort_flux(je,level,blockNo)>0.0_wp)THEN
!            vort_flux_budget_perlevel_pos(1:n_zlev)= vort_flux_budget_perlevel_pos(1:n_zlev)&
!            &+vort_flux(je,level,blockNo)*vn(je,level,blockNo)&
!           &*patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
!           &*patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
!           ELSEIF(vort_flux(je,level,blockNo)<0.0_wp)THEN
!            vort_flux_budget_perlevel_neg(1:n_zlev)= vort_flux_budget_perlevel_neg(1:n_zlev)&
!            &+vort_flux(je,level,blockNo)*vn(je,level,blockNo)&
!           &*patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
!           &*patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
!
!           ENDIF
!           volume(level)=volume(level)+patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
!           &*patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
!           ENDIF
!        END DO
!       END DO
!     END DO ! blockNo = edges_inDomain%start_blockold, edges_inDomain%end_block
!     ! !$OMP END DO NOWAIT
!     ! !$OMP END PARALLEL
!
!     DO level = startLevel, endLevel
!     write(0,*)'global vorticity budget', level,&
!     &vort_flux_budget_perlevel(level)/volume(level),&
!     &vort_flux_budget_perlevel_pos(level)/volume(level),&
!     &vort_flux_budget_perlevel_neg(level)/volume(level)
!     write(0,*)'global vorticity budget pos+neg', level,&
!     &vort_flux_budget_perlevel_pos(level)/volume(level)+vort_flux_budget_perlevel_neg(level)/volume(level)
!     END DO
!
!     stop
  END SUBROUTINE nonlinear_coriolis_3d
  !-------------------------------------------------------------------------


  !-----------------------------------------------------------------------
  ! map_edges2cell_without_height: This sbr is used for kinetic energy calculation
  !>
  !!
  !<Optimize:inUse>
  SUBROUTINE map_edges2cell_no_height_3d( patch_3d, vn_e, operators_coefficients, p_vn_c, opt_startLevel, &
    & opt_endLevel, subset_range, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    TYPE(t_cartesian_coordinates)              :: p_vn_c(:,:,:)  ! output (nproma,n_zlev,alloc_cell_blocks)
    ! intent(inout) for nag compiler
    INTEGER, INTENT(in), OPTIONAL :: opt_startLevel,  opt_endLevel      ! optional vertical start level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: il_e, ib_e
    INTEGER :: cell_index, blockNo, level, ie
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc
    CHARACTER(len=*), PARAMETER :: routine = modname//':map_edges2cell_no_height_3d'
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (no_primal_edges == 3) THEN
      CALL map_edges2cell_no_height_3d_onTriangles( patch_3d, vn_e, operators_coefficients, p_vn_c, opt_startLevel, &
        & opt_endLevel, subset_range, lacc = lzacc )
      RETURN
    ENDIF

#ifdef _OPENACC
    IF (lzacc) CALL finish(routine, 'OpenACC version for no_primal_edges /= 3 currently not tested/validated')
#endif

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2d%cells%ALL
    ENDIF
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_startLevel) ) THEN
      startLevel = opt_startLevel
    ELSE
      startLevel = 1
    END IF
    IF ( PRESENT(opt_endLevel) ) THEN
      endLevel = opt_endLevel
    ELSE
      endLevel = n_zlev
    END IF

    !Calculation of Pv in cartesian coordinates
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      cell_idx_loop: DO cell_index =  start_cell_index, end_cell_index

        endLevel=MIN(endLevel,patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo))
        !$ACC LOOP SEQ
        level_loop: DO level = startLevel, endLevel
           !calculate velocity reconstruction at cell center
          p_vn_c(cell_index,level,blockNo)%x = 0.0_wp

          !$ACC LOOP SEQ
          DO ie=1, no_primal_edges
            il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
            ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)
            IF (il_e > 0) THEN
              p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x&
                & + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,ie)%x&
                & * vn_e(il_e,level,ib_e)*patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
            ENDIF
          END DO

          p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x &
              & / (operators_coefficients%fixed_vol_norm(cell_index,level,blockNo) &
              &   *patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo))

!           IF (abs(p_vn_c(cell_index,level,blockNo)%x(3)) > 0.5) THEN
!             write(0,*) "---------------------------"
!             write(0,*) blockNo, cell_index, " velocity at cell:", p_vn_c(cell_index,level,blockNo)%x
!             p_vn_c(cell_index,level,blockNo)%x = 0.0_wp
!             DO ie=1, no_primal_edges
!               il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
!               ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)
!               IF (il_e > 0) THEN
!                 write(0,*) ie, ":"
!                 write(0,*) " ", patch_2d%edges%cartesian_center(il_e,ib_e)%x
!                 write(0,*) " ", patch_2d%verts%cartesian( &
!                   patch_2d%edges%vertex_idx(il_e,ib_e,1), patch_2d%edges%vertex_blk(il_e,ib_e,1)) %x
!                 write(0,*) " ", patch_2d%verts%cartesian( &
!                   patch_2d%edges%vertex_idx(il_e,ib_e,2), patch_2d%edges%vertex_blk(il_e,ib_e,2)) %x
!                 write(0,*) " ", (patch_2d%verts%cartesian( &
!                   patch_2d%edges%vertex_idx(il_e,ib_e,2), patch_2d%edges%vertex_blk(il_e,ib_e,2))%x + &
!                   patch_2d%verts%cartesian( &
!                   patch_2d%edges%vertex_idx(il_e,ib_e,2), patch_2d%edges%vertex_blk(il_e,ib_e,2))%x) * 0.5_wp
!
!                 write(0,*) " "
!                 write(0,*) " ", patch_2d%cells%cartesian_center(cell_index,blockNo)%x
!                 write(0,*) " ", patch_2d%edges%cartesian_center(il_e,ib_e)%x - &
!                   & patch_2d%cells%cartesian_center(cell_index,blockNo)%x
!
!                 write(0,*) " "
!                 write(0,*) "   ", vn_e(il_e,level,ib_e)
!
!                 p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x&
!                   & + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,ie)%x&
!                   & * vn_e(il_e,level,ib_e)*patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
!                 write(0,*) "   velocity:", p_vn_c(cell_index,level,blockNo)%x
!               ENDIF
!             END DO
!           ENDIF

        END DO level_loop
      END DO  cell_idx_loop
      !$ACC END PARALLEL LOOP

    END DO ! blockNo = all_cells%start_block, all_cells%end_block
    !$ACC WAIT(1)

!     CALL dbg_print('x(1)', p_vn_c(:,:,:)%x(1), &
!       & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     CALL dbg_print('x(2)', p_vn_c(:,:,:)%x(2), &
!       & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     CALL dbg_print('x(3)', p_vn_c(:,:,:)%x(3), &
!       & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)

  END SUBROUTINE map_edges2cell_no_height_3d
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !: This sbr is used for kinetic energy calculation
  SUBROUTINE map_edges2cell_no_height_3d_onTriangles( patch_3d, vn_e, operators_coefficients, p_vn_c, &
    & opt_startLevel, opt_endLevel, subset_range, lacc )

    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    TYPE(t_cartesian_coordinates)              :: p_vn_c(:,:,:)  ! output (nproma,n_zlev,alloc_cell_blocks)
    ! intent(inout) for nag compiler
    INTEGER, INTENT(in), OPTIONAL              :: opt_startLevel,  opt_endLevel      ! optional vertical start level
    TYPE(t_subset_range), TARGET,  OPTIONAL    :: subset_range
    LOGICAL, INTENT(IN), OPTIONAL              :: lacc
    !Local variables
    INTEGER :: startLevel, endLevel, max_level
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: edge_1_index, edge_1_block, edge_2_index, edge_2_block, edge_3_index, edge_3_block
    INTEGER :: cell_index, blockNo, level
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER, DIMENSION(:,:), POINTER           :: dolic_c
    LOGICAL :: lzacc
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    dolic_c    => patch_3d%p_patch_1d(1)%dolic_c
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2d%cells%ALL
    ENDIF
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_startLevel) ) THEN
      startLevel = opt_startLevel
    ELSE
      startLevel = 1
    END IF
    IF ( PRESENT(opt_endLevel) ) THEN
      endLevel = opt_endLevel
    ELSE
      endLevel = n_zlev
    END IF

    CALL set_acc_host_or_device(lzacc, lacc)

    !Calculation of Pv in cartesian coordinates
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, &
!ICON_OMP edge_1_index, edge_1_block, edge_2_index, edge_2_block, edge_3_index, edge_3_block,  &
!ICON_OMP level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      max_level = MIN(MAXVAL(dolic_c(start_cell_index:end_cell_index,blockNo)),endLevel)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
#if defined(__LVECTOR__) && !defined(__LVEC_BITID__)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC   PRIVATE(edge_1_index, edge_1_block, edge_2_index, edge_2_block, edge_3_index, edge_3_block)
      DO level = startLevel, max_level
        DO cell_index = start_cell_index, end_cell_index
          IF (dolic_c(cell_index,blockNo) < level) CYCLE
#else
      !$ACC LOOP GANG VECTOR
      DO cell_index = start_cell_index, end_cell_index
#endif

        edge_1_index = patch_2d%cells%edge_idx(cell_index,blockNo,1)
        edge_1_block = patch_2d%cells%edge_blk(cell_index,blockNo,1)
        edge_2_index = patch_2d%cells%edge_idx(cell_index,blockNo,2)
        edge_2_block = patch_2d%cells%edge_blk(cell_index,blockNo,2)
        edge_3_index = patch_2d%cells%edge_idx(cell_index,blockNo,3)
        edge_3_block = patch_2d%cells%edge_blk(cell_index,blockNo,3)

#if defined(__LVECTOR__) && !defined(__LVEC_BITID__)
#else
        !$ACC LOOP SEQ
        DO level = startLevel, MIN(dolic_c(cell_index,blockNo), endLevel)
#endif

          p_vn_c(cell_index,level,blockNo)%x =                                            &
            & (  operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,1)%x  &
            &      * vn_e(edge_1_index,level,edge_1_block)                                &
            &*patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_index,level,edge_1_block)       &
            &  + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,2)%x  &
            &      * vn_e(edge_2_index,level,edge_2_block)                                &
            &*patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_index,level,edge_2_block)       &
            &  + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,3)%x  &
            &       * vn_e(edge_3_index,level,edge_3_block)                               &
            &*patch_3d%p_patch_1d(1)%prism_thick_e(edge_3_index,level,edge_3_block))      &
            & / (operators_coefficients%fixed_vol_norm(cell_index,level,blockNo)          &
            &    * patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo))
        END DO
      END DO
      !$ACC END PARALLEL

    END DO ! blockNo = all_cells%start_block, all_cells%end_block
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2cell_no_height_3d_onTriangles
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! map_edges2cell_with_height: This sbr is used by the GMRedi
  !>
  !!
  !<Optimize:inUse>
  SUBROUTINE map_edges2cell_with_height_3d( patch_3d, vn_e, operators_coefficients, p_vn_c, opt_startLevel, &
    & opt_endLevel, subset_range)

    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    TYPE(t_cartesian_coordinates)              :: p_vn_c(:,:,:)  ! output (nproma,n_zlev,alloc_cell_blocks)
    ! intent(inout) for nag compiler
    INTEGER, INTENT(in), OPTIONAL :: opt_startLevel,  opt_endLevel      ! optional vertical start level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: il_e, ib_e
    INTEGER :: cell_index, blockNo, level, ie
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
     IF (patch_3D%p_patch_2D(1)%cells%max_connectivity == 3.and. fast_performance_level > 10) THEN
       CALL map_edges2cell_with_height_3d_onTriangles( patch_3d, vn_e, operators_coefficients, p_vn_c, opt_startLevel, &
         & opt_endLevel, subset_range)
       RETURN
     ENDIF
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2d%cells%ALL
    ENDIF
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_startLevel) ) THEN
      startLevel = opt_startLevel
    ELSE
      startLevel = 1
    END IF
    IF ( PRESENT(opt_endLevel) ) THEN
      endLevel = opt_endLevel
    ELSE
      endLevel = n_zlev
    END IF

    !Calculation of Pv in cartesian coordinates
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)

      cell_idx_loop: DO cell_index =  start_cell_index, end_cell_index

        endLevel=MIN(endLevel, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo))
        level_loop: DO level = startLevel, endLevel
           !calculate velocity reconstruction at cell center
          p_vn_c(cell_index,level,blockNo)%x = 0.0_wp

          DO ie=1, patch_3D%p_patch_2D(1)%cells%max_connectivity!no_primal_edges
            il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
            ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)

            IF (il_e > 0) THEN
            p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x&
              & + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,ie)%x&
              & * vn_e(il_e,level,ib_e)*(patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e))
            ENDIF

          END DO

          p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x &
              & / operators_coefficients%fixed_vol_norm(cell_index,level,blockNo)



        END DO level_loop
      END DO  cell_idx_loop

    END DO ! blockNo = all_cells%start_block, all_cells%end_block
! CALL dbg_print('thick_e', patch_3d%p_patch_1d(1)%prism_thick_e, &
!       & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%edges%owned)
!     CALL dbg_print('x(1)', p_vn_c(:,:,:)%x(1), &
!       & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     CALL dbg_print('x(2)', p_vn_c(:,:,:)%x(2), &
!       & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     CALL dbg_print('x(3)', p_vn_c(:,:,:)%x(3), &
!       & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)

  END SUBROUTINE map_edges2cell_with_height_3d
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !: This sbr is used by GMRedi
  SUBROUTINE map_edges2cell_with_height_3d_onTriangles( patch_3d, vn_e, operators_coefficients, p_vn_c, &
    & opt_startLevel, opt_endLevel, subset_range)

    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    TYPE(t_cartesian_coordinates)              :: p_vn_c(:,:,:)  ! output (nproma,n_zlev,alloc_cell_blocks)
    ! intent(inout) for nag compiler
    INTEGER, INTENT(in), OPTIONAL :: opt_startLevel,  opt_endLevel      ! optional vertical start level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: edge_1_index, edge_1_block, edge_2_index, edge_2_block, edge_3_index, edge_3_block
    INTEGER :: cell_index, blockNo, level
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2d%cells%ALL
    ENDIF
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_startLevel) ) THEN
      startLevel = opt_startLevel
    ELSE
      startLevel = 1
    END IF
    IF ( PRESENT(opt_endLevel) ) THEN
      endLevel = opt_endLevel
    ELSE
      endLevel = n_zlev
    END IF

    !Calculation of Pv in cartesian coordinates
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, &
!ICON_OMP edge_1_index, edge_1_block, edge_2_index, edge_2_block, edge_3_index, edge_3_block,  &
!ICON_OMP level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      DO cell_index = start_cell_index, end_cell_index

        edge_1_index = patch_2d%cells%edge_idx(cell_index,blockNo,1)
        edge_1_block = patch_2d%cells%edge_blk(cell_index,blockNo,1)
        edge_2_index = patch_2d%cells%edge_idx(cell_index,blockNo,2)
        edge_2_block = patch_2d%cells%edge_blk(cell_index,blockNo,2)
        edge_3_index = patch_2d%cells%edge_idx(cell_index,blockNo,3)
        edge_3_block = patch_2d%cells%edge_blk(cell_index,blockNo,3)


        DO level = startLevel, MIN(patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo), endLevel)
          p_vn_c(cell_index,level,blockNo)%x =                                            &
            & (  operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,1)%x  &
            &      * vn_e(edge_1_index,level,edge_1_block)                                &
            &*patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_index,level,edge_1_block)       &
            &  + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,2)%x  &
            &      * vn_e(edge_2_index,level,edge_2_block)                                &
            &*patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_index,level,edge_2_block)       &
            &  + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,3)%x  &
            &       * vn_e(edge_3_index,level,edge_3_block)                               &
            &*patch_3d%p_patch_1d(1)%prism_thick_e(edge_3_index,level,edge_3_block))      &
            & / operators_coefficients%fixed_vol_norm(cell_index,level,blockNo)

        END DO

      END DO

    END DO ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2cell_with_height_3d_onTriangles
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  !This sbr is used in continuity equation for the calculation of vertical velocity, and
  !in setting up he right hand side for the free-surface equation.
  !With an additional (tracer)-scalar it is used optionally in tracer transport.
  !It has a 2d-version which is used in left-hand-side of free-surface equation.
  SUBROUTINE map_edges2edges_viacell_3d_mlev_const_z( patch_3d, vn_e, operators_coefficients, out_vn_e, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL, INTENT(in), OPTIONAL        :: lacc
    !Local variables
    INTEGER :: startLevel, endLevel, start_block, end_block
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, blockNo, level, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc
    ! Pointers needed for GPU/OpenACC
    INTEGER , POINTER :: lsm_e(:,:,:)
    INTEGER , POINTER :: cell_idx(:,:,:)
    INTEGER , POINTER :: cell_blk(:,:,:)
    INTEGER , POINTER :: edge_idx(:,:,:)
    INTEGER , POINTER :: edge_blk(:,:,:)
    REAL(wp), POINTER :: prism_thick_e(:,:,:)
    REAL(wp), POINTER :: coeffs(:,:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':map_edges2edges_viacell_3d_mlev_const_z'
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF ( patch_2d%cells%max_connectivity == 3 .and. fast_performance_level > 10) THEN
      !$ACC DATA COPYIN(patch_3d%p_patch_2d(1)%nblks_e) IF(lzacc)
      CALL map_edges2edges_viacell_3d_mlev_constZ_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e, lacc=lzacc)
      !$ACC END DATA
      RETURN
    ENDIF
    !-----------------------------------------------------------------------
#ifdef _OPENACC
    IF (lzacc) CALL finish(routine, 'OpenACC version for max_connectivity /=3 or perf_level currently not tested/validated')
#endif

    edges_inDomain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    start_block = edges_inDomain%start_block
    end_block = edges_inDomain%end_block

    lsm_e    => patch_3d%lsm_e
    cell_idx => patch_2d%edges%cell_idx
    cell_blk => patch_2d%edges%cell_blk
    edge_idx => patch_2d%cells%edge_idx
    edge_blk => patch_2d%cells%edge_blk
    prism_thick_e => patch_3d%p_patch_1d(1)%prism_thick_e
    coeffs => operators_coefficients%edge2edge_viacell_coeff

    !$ACC DATA COPYIN(cell_blk, cell_idx, coeffs, edge_blk, edge_idx, lsm_e, prism_thick_e, vn_e) &
    !$ACC   COPY(out_vn_e) IF(lzacc)

    DO blockNo = start_block, end_block
      CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      out_vn_e(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      level_loop_e: DO level = startLevel, endLevel
        edge_idx_loop: DO je =  start_edge_index, end_edge_index
          IF (lsm_e(je,level,blockNo) == sea) THEN
            out_vn_e(je,level,blockNo)= 0.0_wp

            ictr = 0
            il_c = cell_idx(je,blockNo,1)
            ib_c = cell_blk(je,blockNo,1)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,level,ib_c)

            !$ACC LOOP SEQ
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = edge_idx(il_c,ib_c,ie)
              ib_e = edge_blk(il_c,ib_c,ie)
              IF (il_e > 0) THEN
                thick_edge = prism_thick_e(il_e,level,ib_e)

                out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
                  & +vn_e(il_e,level,ib_e)*coeffs(je,level,blockNo,ictr)&
                  & *thick_edge
              ENDIF

            END DO

            ictr = no_primal_edges
            il_c = cell_idx(je,blockNo,2)
            ib_c = cell_blk(je,blockNo,2)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,level,ib_c)

            !$ACC LOOP SEQ
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = edge_idx(il_c,ib_c,ie)
              ib_e = edge_blk(il_c,ib_c,ie)
              IF (il_e > 0) THEN
                thick_edge = prism_thick_e(il_e,level,ib_e)

                out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
                  & +vn_e(il_e,level,ib_e)*coeffs(je,level,blockNo,ictr)&
                  & *thick_edge
              ENDIF

            END DO
          ENDIF
        END DO edge_idx_loop
      END DO level_loop_e
      !$ACC END PARALLEL
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE map_edges2edges_viacell_3d_mlev_const_z
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! the map_edges2edges_viacell_3d_mlev_constZ optimized for triangles
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_3d_mlev_constZ_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, blockNo, level, start_block, end_block

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
    INTEGER :: edge_11_block, edge_12_block, edge_13_block
    INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
    INTEGER :: edge_21_block, edge_22_block, edge_23_block

    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc
    ! Pointers needed for GPU/OpenACC
    INTEGER , POINTER :: cell_idx(:,:,:)
    INTEGER , POINTER :: cell_blk(:,:,:)
    INTEGER , POINTER :: edge_idx(:,:,:)
    INTEGER , POINTER :: edge_blk(:,:,:)
    REAL(wp), POINTER :: prism_thick_e(:,:,:)
    INTEGER, POINTER :: dolic_e(:,:)
    REAL(wp), POINTER :: coeffs(:,:,:,:)
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')

    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    start_block = edges_in_domain%start_block
    end_block = edges_in_domain%end_block

    cell_idx => patch_2d%edges%cell_idx
    cell_blk => patch_2d%edges%cell_blk
    edge_idx => patch_2d%cells%edge_idx
    edge_blk => patch_2d%cells%edge_blk
    prism_thick_e => patch_3d%p_patch_1d(1)%prism_thick_e
    dolic_e => patch_3d%p_patch_1d(1)%dolic_e
    coeffs => operators_coefficients%edge2edge_viacell_coeff
    !-----------------------------------------------------------------------
    ! All derived variables used to specify the size of OpenACC arrays in the declaration have to be present
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, je, cell_1_index, cell_1_block, &
!ICON_OMP   cell_2_index, cell_2_block, edge_11_index, edge_12_index, edge_13_index, &
!ICON_OMP  edge_11_block, edge_12_block, edge_13_block, edge_21_index, edge_22_index, &
!ICON_OMP  edge_23_index, edge_21_block, edge_22_block, edge_23_block, level)  ICON_OMP_DEFAULT_SCHEDULE

    DO blockNo = start_block, end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      out_vn_e(:, :, blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO je =  start_edge_index, end_edge_index

        IF (dolic_e(je,blockNo) < 1) CYCLE

        cell_1_index = cell_idx(je,blockNo,1)
        cell_1_block = cell_blk(je,blockNo,1)
        cell_2_index = cell_idx(je,blockNo,2)
        cell_2_block = cell_blk(je,blockNo,2)

        edge_11_index = edge_idx(cell_1_index, cell_1_block, 1)
        edge_12_index = edge_idx(cell_1_index, cell_1_block, 2)
        edge_13_index = edge_idx(cell_1_index, cell_1_block, 3)
        edge_11_block = edge_blk(cell_1_index, cell_1_block, 1)
        edge_12_block = edge_blk(cell_1_index, cell_1_block, 2)
        edge_13_block = edge_blk(cell_1_index, cell_1_block, 3)

        edge_21_index = edge_idx(cell_2_index, cell_2_block, 1)
        edge_22_index = edge_idx(cell_2_index, cell_2_block, 2)
        edge_23_index = edge_idx(cell_2_index, cell_2_block, 3)
        edge_21_block = edge_blk(cell_2_index, cell_2_block, 1)
        edge_22_block = edge_blk(cell_2_index, cell_2_block, 2)
        edge_23_block = edge_blk(cell_2_index, cell_2_block, 3)

        ! levels
        DO level = startLevel, dolic_e(je,blockNo)

          out_vn_e(je, level, blockNo) =  &
            & (  vn_e(edge_11_index, level, edge_11_block) * coeffs(je, level, blockNo, 1)               &
            & * prism_thick_e(edge_11_index, level, edge_11_block)  +  &
            & vn_e(edge_12_index, level, edge_12_block) * coeffs(je, level, blockNo, 2)               &
            & * prism_thick_e(edge_12_index, level, edge_12_block)  +  &
            & vn_e(edge_13_index, level, edge_13_block) * coeffs(je, level, blockNo, 3)               &
            & * prism_thick_e(edge_13_index, level, edge_13_block)     &
            & )                               &
            & + &
            & (  vn_e(edge_21_index, level, edge_21_block) * coeffs(je, level, blockNo, 4)               &
            & * prism_thick_e(edge_21_index, level, edge_21_block)  +  &
            & vn_e(edge_22_index, level, edge_22_block) * coeffs(je, level, blockNo, 5)               &
            & * prism_thick_e(edge_22_index, level, edge_22_block)  +  &
            & vn_e(edge_23_index, level, edge_23_block) * coeffs(je, level, blockNo, 6)               &
            & * prism_thick_e(edge_23_index, level, edge_23_block)     &
            & )

        END DO ! levels

      END DO
      !$ACC END PARALLEL

    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)

!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

  END SUBROUTINE map_edges2edges_viacell_3d_mlev_constZ_onTriangles
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! the map_edges2edges_viacell_3d_mlev_constZ optimized for triangles
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_3d_zstar( patch_3d, vn_e, operators_coefficients, stretch_e, out_vn_e, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                 :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL, INTENT(in), OPTIONAL        :: lacc
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, blockNo, level

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
    INTEGER :: edge_11_block, edge_12_block, edge_13_block
    INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
    INTEGER :: edge_21_block, edge_22_block, edge_23_block

    REAL(wp), POINTER :: coeffs(:,:,:,:)
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    coeffs => operators_coefficients%edge2edge_viacell_coeff
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, je, cell_1_index, cell_1_block, &
!ICON_OMP   cell_2_index, cell_2_block, edge_11_index, edge_12_index, edge_13_index, &
!ICON_OMP  edge_11_block, edge_12_block, edge_13_block, edge_21_index, edge_22_index, &
!ICON_OMP  edge_23_index, edge_21_block, edge_22_block, edge_23_block, level)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      out_vn_e(:, :, blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO je =  start_edge_index, end_edge_index

        IF (patch_3d%p_patch_1d(1)%dolic_e(je,blockNo) < 1) CYCLE

        cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
        cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)

        edge_11_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
        edge_12_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
        edge_13_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
        edge_11_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
        edge_12_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
        edge_13_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)

        edge_21_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
        edge_22_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
        edge_23_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
        edge_21_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
        edge_22_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
        edge_23_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)

        ! levels
        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          out_vn_e(je, level, blockNo) =  &
            & (  vn_e(edge_11_index, level, edge_11_block) * coeffs(je, level, blockNo, 1)    &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_11_index, level, edge_11_block)     &
            & * stretch_e(edge_11_index, edge_11_block)  +                                    &
            & vn_e(edge_12_index, level, edge_12_block) * coeffs(je, level, blockNo, 2)       &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_12_index, level, edge_12_block)     &
            & * stretch_e(edge_12_index, edge_12_block)  +                                    &
            & vn_e(edge_13_index, level, edge_13_block) * coeffs(je, level, blockNo, 3)       &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_13_index, level, edge_13_block)     &
            & * stretch_e(edge_13_index, edge_13_block)                                       &
            & )                                                                               &
            & + &
            & (  vn_e(edge_21_index, level, edge_21_block) * coeffs(je, level, blockNo, 4)    &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_21_index, level, edge_21_block)     &
            & * stretch_e(edge_21_index, edge_21_block)  +                                    &
            & vn_e(edge_22_index, level, edge_22_block) * coeffs(je, level, blockNo, 5)       &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_22_index, level, edge_22_block)     &
            & * stretch_e(edge_22_index, edge_22_block)  +                                    &
            & vn_e(edge_23_index, level, edge_23_block) * coeffs(je, level, blockNo, 6)       &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_23_index, level, edge_23_block)     &
            & * stretch_e(edge_23_index, edge_23_block)                                       &
            & )

        END DO ! levels

      END DO
      !$ACC END PARALLEL

    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
  END SUBROUTINE map_edges2edges_3d_zstar
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  !Same sbr as above, now with an additional scalar (e.g. tracer). It is used optionally in tracer transport.
  SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs( patch_3d, vn_e, operators_coefficients, out_vn_e, scalar, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)       :: scalar(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    !Local variables
    INTEGER :: startLevel, endLevel, start_block, end_block
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, blockNo, level, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc
    ! Pointers needed for GPU/OpenACC
    INTEGER, POINTER :: lsm_e(:,:,:)
    INTEGER , POINTER :: cell_idx(:,:,:)
    INTEGER , POINTER :: cell_blk(:,:,:)
    INTEGER , POINTER :: edge_idx(:,:,:)
    INTEGER , POINTER :: edge_blk(:,:,:)
    REAL(wp), POINTER :: prism_thick_e(:,:,:)
    REAL(wp), POINTER :: edge2edge_viacell_coeff(:,:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':map_edges2edges_viacell_3d_mlev_constZs'

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF ( patch_2d%cells%max_connectivity == 3) THEN
      IF (fast_performance_level > 10 ) THEN
        CALL map_edges2edges_viacell_3d_mlev_constZs_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e, scalar, lacc=lzacc)
        RETURN
      ENDIF
    ENDIF

#ifdef _OPENACC
    IF (lzacc) CALL finish(routine, 'OpenACC version for max_connectivity /= 3 currently not tested/validated')
#endif

    !-----------------------------------------------------------------------
    edges_inDomain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    start_block = edges_inDomain%start_block
    end_block = edges_inDomain%end_block

    lsm_e => patch_3d%lsm_e
    cell_idx => patch_2d%edges%cell_idx
    cell_blk => patch_2d%edges%cell_blk
    edge_idx => patch_2d%cells%edge_idx
    edge_blk => patch_2d%cells%edge_blk
    prism_thick_e => patch_3d%p_patch_1d(1)%prism_thick_e
    edge2edge_viacell_coeff => operators_coefficients%edge2edge_viacell_coeff

    !$ACC DATA COPYIN(cell_idx, cell_blk, edge2edge_viacell_coeff, edge_idx, edge_blk, lsm_e, prism_thick_e) &
    !$ACC   COPYIN(scalar, vn_e) &
    !$ACC   COPY(out_vn_e) IF(lzacc)

    DO blockNo = start_block, end_block
      CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      out_vn_e(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      level_loop_e2: DO level = startLevel, endLevel
        edge_idx_loop2: DO je =  start_edge_index, end_edge_index
          IF (lsm_e(je,level,blockNo) == sea) THEN
            out_vn_e(je,level,blockNo)= 0.0_wp

            ictr = 0
            il_c        = cell_idx(je,blockNo,1)
            ib_c        = cell_blk(je,blockNo,1)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,level,ib_c)
            scalar_cell = scalar(il_c,level,ib_c)

            !$ACC LOOP SEQ
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = edge_idx(il_c,ib_c,ie)
              ib_e = edge_blk(il_c,ib_c,ie)
              thick_edge = prism_thick_e(il_e,level,ib_e)

              out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
                & +vn_e(il_e,level,ib_e)*edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                & *thick_edge*scalar_cell

            END DO

            ictr        = no_primal_edges
            il_c        = cell_idx(je,blockNo,2)
            ib_c        = cell_blk(je,blockNo,2)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,level,ib_c)
            scalar_cell = scalar(il_c,level,ib_c)

            !$ACC LOOP SEQ
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = edge_idx(il_c,ib_c,ie)
              ib_e = edge_blk(il_c,ib_c,ie)
              thick_edge = prism_thick_e(il_e,level,ib_e)

              out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
                & +vn_e(il_e,level,ib_e)*edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                & *thick_edge*scalar_cell

            END DO
          ENDIF
        END DO edge_idx_loop2
      END DO level_loop_e2
      !$ACC END PARALLEL
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! the map_edges2edges_viacell_3d_mlev_constZs optimized for triangles
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e,scalar, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)                 :: scalar(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL        :: lacc
    !Local variables
    INTEGER :: startLevel, endLevel, start_block, end_block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp), POINTER :: coeffs(:,:,:,:)
    LOGICAL :: lzacc

    ! omp private
    ! defined in a struct for omp clearness
    TYPE omp_local_private

      INTEGER :: start_edge_index, end_edge_index
      INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
      INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
      INTEGER :: edge_11_block, edge_12_block, edge_13_block
      INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
      INTEGER :: edge_21_block, edge_22_block, edge_23_block

    END TYPE omp_local_private

!     TYPE(omp_local_private) :: omp_this
    INTEGER :: je, blockNo, level, max_dolic_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
    INTEGER :: edge_11_block, edge_12_block, edge_13_block
    INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
    INTEGER :: edge_21_block, edge_22_block, edge_23_block
    ! Pointers needed for GPU/OpenACC
    INTEGER , POINTER :: cell_idx(:,:,:)
    INTEGER , POINTER :: cell_blk(:,:,:)
    INTEGER , POINTER :: edge_idx(:,:,:)
    INTEGER , POINTER :: edge_blk(:,:,:)
    REAL(wp), POINTER :: prism_thick_e(:,:,:)
    INTEGER, POINTER :: dolic_e(:,:)
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')

    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    start_block = edges_in_domain%start_block
    end_block = edges_in_domain%end_block
    coeffs => operators_coefficients%edge2edge_viacell_coeff

    cell_idx => patch_2d%edges%cell_idx
    cell_blk => patch_2d%edges%cell_blk
    edge_idx => patch_2d%cells%edge_idx
    edge_blk => patch_2d%cells%edge_blk
    prism_thick_e => patch_3d%p_patch_1d(1)%prism_thick_e
    dolic_e => patch_3d%p_patch_1d(1)%dolic_e
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE( je, level, start_edge_index, end_edge_index, &
!ICON_OMP cell_1_index, cell_2_index, cell_1_block, cell_2_block, &
!ICON_OMP edge_11_index, edge_12_index, edge_13_index, &
!ICON_OMP edge_11_block, edge_12_block, edge_13_block, &
!ICON_OMP edge_21_index, edge_22_index, edge_23_index, &
!ICON_OMP edge_21_block, edge_22_block, edge_23_block) ICON_OMP_DEFAULT_SCHEDULE

    DO blockNo = start_block, end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

! Note that this loop structure is currenlty suboptimal on Aurora

#ifdef __LVECTOR__
      max_dolic_e = -1
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) REDUCTION(MAX: max_dolic_e) IF(lzacc)
      DO je = start_edge_index, end_edge_index
        max_dolic_e = MAX(max_dolic_e, dolic_e(je,blockNo))
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO level = startLevel, max_dolic_e
        DO je = start_edge_index, end_edge_index
          IF (dolic_e(je,blockNo) < level) CYCLE
#else           
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO je =  start_edge_index, end_edge_index
#endif
        cell_1_index = cell_idx(je,blockNo,1)
        cell_1_block = cell_blk(je,blockNo,1)
        cell_2_index = cell_idx(je,blockNo,2)
        cell_2_block = cell_blk(je,blockNo,2)

        edge_11_index = edge_idx(cell_1_index, cell_1_block, 1)
        edge_12_index = edge_idx(cell_1_index, cell_1_block, 2)
        edge_13_index = edge_idx(cell_1_index, cell_1_block, 3)
        edge_11_block = edge_blk(cell_1_index, cell_1_block, 1)
        edge_12_block = edge_blk(cell_1_index, cell_1_block, 2)
        edge_13_block = edge_blk(cell_1_index, cell_1_block, 3)

        edge_21_index = edge_idx(cell_2_index, cell_2_block, 1)
        edge_22_index = edge_idx(cell_2_index, cell_2_block, 2)
        edge_23_index = edge_idx(cell_2_index, cell_2_block, 3)
        edge_21_block = edge_blk(cell_2_index, cell_2_block, 1)
        edge_22_block = edge_blk(cell_2_index, cell_2_block, 2)
        edge_23_block = edge_blk(cell_2_index, cell_2_block, 3)

#ifndef __LVECTOR__
        ! levels
        DO level = startLevel, dolic_e(je,blockNo)
#endif
          out_vn_e(je, level, blockNo) =  &
            & (  vn_e(edge_11_index, level, edge_11_block) * coeffs(je, level, blockNo, 1)      &
            &    * prism_thick_e(edge_11_index, level, edge_11_block)    &
            &  + vn_e(edge_12_index, level, edge_12_block) * coeffs(je, level, blockNo, 2)      &
            &    * prism_thick_e(edge_12_index, level, edge_12_block)    &
            &  + vn_e(edge_13_index, level, edge_13_block) * coeffs(je, level, blockNo, 3)      &
            &    * prism_thick_e(edge_13_index, level, edge_13_block)    &
            & ) * scalar(cell_1_index, level, cell_1_block)                                     &
            & + &
            & (  vn_e(edge_21_index, level, edge_21_block) * coeffs(je, level, blockNo, 4)           &
            &   * prism_thick_e(edge_21_index, level, edge_21_block)  +  &
            &  vn_e(edge_22_index, level, edge_22_block) * coeffs(je, level, blockNo, 5)             &
            &  * prism_thick_e(edge_22_index, level, edge_22_block)  +   &
            &  vn_e(edge_23_index, level, edge_23_block) * coeffs(je, level, blockNo, 6)             &
            & * prism_thick_e(edge_23_index, level, edge_23_block)       &
            & ) * scalar(cell_2_index, level, cell_2_block)

        END DO !levels

      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs_onTriangles
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! the map_edges2edges_viacell_3d_mlev_constZ optimized for triangles
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_sc_zstar( patch_3d, vn_e, scalar, operators_coefficients, stretch_e, out_vn_e, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)                 :: scalar(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                 :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL, INTENT(in), OPTIONAL        :: lacc
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, blockNo, level
    LOGICAL :: lzacc

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
    INTEGER :: edge_11_block, edge_12_block, edge_13_block
    INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
    INTEGER :: edge_21_block, edge_22_block, edge_23_block

    REAL(wp), POINTER :: coeffs(:,:,:,:)
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    coeffs => operators_coefficients%edge2edge_viacell_coeff
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(patch_3d%p_patch_2d(1)%nblks_e, patch_3d%p_patch_2d(1)%alloc_cell_blocks) IF(lzacc)
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, je, cell_1_index, cell_1_block, &
!ICON_OMP   cell_2_index, cell_2_block, edge_11_index, edge_12_index, edge_13_index, &
!ICON_OMP  edge_11_block, edge_12_block, edge_13_block, edge_21_index, edge_22_index, &
!ICON_OMP  edge_23_index, edge_21_block, edge_22_block, edge_23_block, level)  ICON_OMP_DEFAULT_SCHEDULE

    !$ACC DATA COPYIN(patch_2d%edges%cell_blk, patch_2d%edges%cell_idx, coeffs) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%dolic_e, patch_2d%cells%edge_blk, patch_2d%cells%edge_idx) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%prism_thick_e, vn_e, scalar, stretch_e) &
    !$ACC   COPY(out_vn_e) IF(lzacc)
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      out_vn_e(:, :, blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO je =  start_edge_index, end_edge_index

        IF (patch_3d%p_patch_1d(1)%dolic_e(je,blockNo) < 1) CYCLE

        cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
        cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)

        edge_11_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
        edge_12_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
        edge_13_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
        edge_11_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
        edge_12_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
        edge_13_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)

        edge_21_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
        edge_22_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
        edge_23_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
        edge_21_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
        edge_22_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
        edge_23_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)

        ! levels
        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          out_vn_e(je, level, blockNo) =  &
            & (  vn_e(edge_11_index, level, edge_11_block) * coeffs(je, level, blockNo, 1)    &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_11_index, level, edge_11_block)     &
            & * stretch_e(edge_11_index, edge_11_block)  +                                    &
            & vn_e(edge_12_index, level, edge_12_block) * coeffs(je, level, blockNo, 2)       &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_12_index, level, edge_12_block)     &
            & * stretch_e(edge_12_index, edge_12_block)  +                                    &
            & vn_e(edge_13_index, level, edge_13_block) * coeffs(je, level, blockNo, 3)       &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_13_index, level, edge_13_block)     &
            & * stretch_e(edge_13_index, edge_13_block)                                       &
            & )* scalar(cell_1_index, level, cell_1_block)                                    &
            & + &
            & (  vn_e(edge_21_index, level, edge_21_block) * coeffs(je, level, blockNo, 4)    &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_21_index, level, edge_21_block)     &
            & * stretch_e(edge_21_index, edge_21_block)  +                                    &
            & vn_e(edge_22_index, level, edge_22_block) * coeffs(je, level, blockNo, 5)       &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_22_index, level, edge_22_block)     &
            & * stretch_e(edge_22_index, edge_22_block)  +                                    &
            & vn_e(edge_23_index, level, edge_23_block) * coeffs(je, level, blockNo, 6)       &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_23_index, level, edge_23_block)     &
            & * stretch_e(edge_23_index, edge_23_block)                                       &
            & ) * scalar(cell_2_index, level, cell_2_block)

        END DO ! levels

      END DO
      !$ACC END PARALLEL

    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)
    !$ACC END DATA

!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

  !$ACC END DATA
  END SUBROUTINE map_edges2edges_sc_zstar
  !-----------------------------------------------------------------------------




  !-----------------------------------------------------------------------------
!<Optimize:inUse>
  !This sbr is used in left-hand-side of free-surface equation for certain solver configurations.
  SUBROUTINE map_edges2edges_viacell_2D_constZ( patch_3d, vn_e, operators_coefficients, out_vn_e, lacc)!&   subset_range)

    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    !Local variables
    INTEGER :: startLevel, endLevel, start_block, end_block
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, blockNo, level, ie
    !REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell!, thick_frac
    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER :: patch_2d
    ! Pointers needed for GPU/OpenACC
    INTEGER, POINTER :: lsm_e(:,:,:)
    INTEGER , POINTER :: cell_idx(:,:,:)
    INTEGER , POINTER :: cell_blk(:,:,:)
    INTEGER , POINTER :: edge_idx(:,:,:)
    INTEGER , POINTER :: edge_blk(:,:,:)
    REAL(wp), POINTER :: prism_thick_e(:,:,:)
    INTEGER, POINTER :: dolic_e(:,:)
    REAL(wp), POINTER :: coeffs(:,:,:,:)
    LOGICAL :: lzacc

    CHARACTER(len=*), PARAMETER :: routine = modname//':map_edges2edges_viacell_2D_constZ'

    patch_2d   => patch_3d%p_patch_2d(1)

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    IF ( patch_2d%cells%max_connectivity == 3 .and. fast_performance_level > 10 ) THEN
      CALL map_edges2edges_viacell_2D_constZ_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e, lacc=lzacc)
      RETURN
    ENDIF
    !-----------------------------------------------------------------------

#ifdef _OPENACC
    IF (lzacc) CALL finish(routine, 'OpenACC version currently not tested/validated')
#endif

    edges_inDomain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = startLevel!n_zlev
    start_block = edges_inDomain%start_block
    end_block = edges_inDomain%end_block

    lsm_e => patch_3d%lsm_e
    cell_idx => patch_2d%edges%cell_idx
    cell_blk => patch_2d%edges%cell_blk
    edge_idx => patch_2d%cells%edge_idx
    edge_blk => patch_2d%cells%edge_blk
    prism_thick_e => patch_3d%p_patch_1d(1)%prism_thick_e
    dolic_e => patch_3d%p_patch_1d(1)%dolic_e
    coeffs => operators_coefficients%edge2edge_viacell_coeff

!$OMP PARALLEL DO PRIVATE(start_edge_index,end_edge_index,je,endLevel,level,&
!$OMP & ictr,il_c,ib_c,ie,il_e,ib_e,thick_edge) SCHEDULE(GUIDED)

    !$ACC DATA COPYIN(cell_blk, cell_idx, coeffs, dolic_e, edge_blk, edge_idx, lsm_e, prism_thick_e, vn_e) &
    !$ACC   COPY(out_vn_e) IF(lzacc)

    DO blockNo = start_block, end_block
      CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      out_vn_e(:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      edge_idx_loop: DO je =  start_edge_index, end_edge_index
        endLevel = dolic_e(je,blockNo)
        level_loop_e: DO level = startLevel, endLevel

          IF (lsm_e(je,level,blockNo) == sea) THEN

            ictr = 0
            il_c = cell_idx(je,blockNo,1)
            ib_c = cell_blk(je,blockNo,1)
            !thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)

            !$ACC LOOP SEQ
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = edge_idx(il_c,ib_c,ie)
              ib_e = edge_blk(il_c,ib_c,ie)
              IF (il_e > 0) THEN
                thick_edge = prism_thick_e(il_e,level,ib_e)

                out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
                  & +vn_e(il_e,ib_e)*coeffs(je,level,blockNo,ictr)&
                  & *thick_edge
              ENDIF

            END DO

            ictr = no_primal_edges
            il_c = cell_idx(je,blockNo,2)
            ib_c = cell_blk(je,blockNo,2)
            !thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)

            !$ACC LOOP SEQ
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = edge_idx(il_c,ib_c,ie)
              ib_e = edge_blk(il_c,ib_c,ie)
              IF (il_e > 0) THEN
                thick_edge = prism_thick_e(il_e,level,ib_e)

                out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
                  & +vn_e(il_e,ib_e)*coeffs(je,level,blockNo,ictr)&
                  & *thick_edge
              ENDIF

            END DO
          ENDIF
        END DO  level_loop_e
      END DO edge_idx_loop
      !$ACC END PARALLEL
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE map_edges2edges_viacell_2D_constZ
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_2D_constZ_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e, lacc)!&   subset_range)

    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL, INTENT(in), OPTIONAL              :: lacc

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, blockNo, start_edge_index, end_edge_index, level
    INTEGER :: start_block, end_block
    LOGICAL :: lzacc

    TYPE(t_subset_range), POINTER :: edges_indomain
    TYPE(t_patch), POINTER :: patch_2d
    ! Pointers needed for GPU/OpenACC
    INTEGER , POINTER :: cell_idx(:,:,:)
    INTEGER , POINTER :: cell_blk(:,:,:)
    INTEGER , POINTER :: edge_idx(:,:,:)
    INTEGER , POINTER :: edge_blk(:,:,:)
    INTEGER, POINTER :: dolic_e(:,:)
    REAL(wp), POINTER :: coeffs(:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':map_edges2edges_viacell_2D_constZ_onTriangles'

    !-----------------------------------------------------------------------
!     IF (no_primal_edges /= 3) &
!       & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_indomain    => patch_2d%edges%in_domain
    coeffs        => operators_coefficients%edge2edge_viacell_coeff_all

    cell_idx => patch_2d%edges%cell_idx
    cell_blk => patch_2d%edges%cell_blk
    edge_idx => patch_2d%cells%edge_idx
    edge_blk => patch_2d%cells%edge_blk
    dolic_e => patch_3d%p_patch_1d(1)%dolic_e

    start_block = edges_indomain%start_block
    end_block = edges_indomain%end_block

    CALL set_acc_host_or_device(lzacc, lacc)

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, level, cell_1_index, cell_1_block, &
!ICON_OMP  cell_2_index, cell_2_block, edge_1_1_index, edge_1_2_index, edge_1_3_index, &
!ICON_OMP  edge_1_1_block, edge_1_2_block, edge_1_3_block, edge_2_1_index, edge_2_2_index, &
!ICON_OMP  edge_2_3_index, edge_2_1_block, edge_2_2_block, edge_2_3_block)  ICON_OMP_DEFAULT_SCHEDULE

    !$ACC DATA COPYIN(cell_blk, cell_idx, coeffs, dolic_e, edge_blk, edge_idx, vn_e) &
    !$ACC   COPY(out_vn_e) IF(lzacc)

    DO blockNo = start_block, end_block
      CALL get_index_range(edges_indomain, blockNo, start_edge_index, end_edge_index)
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO je = start_edge_index, end_edge_index

        out_vn_e(je,blockNo) = 0.0_wp

        DO level = 1, MIN(1,dolic_e(je,blockNo))
!         IF (patch_3d%lsm_e(je,1,blockNo) == sea) THEN ! this if should be removed

          ! get the two cells of the edge
          cell_1_index = cell_idx(je,blockNo,1)
          cell_1_block = cell_blk(je,blockNo,1)
          cell_2_index = cell_idx(je,blockNo,2)
          cell_2_block = cell_blk(je,blockNo,2)

          ! get the six edges of the two cells
          edge_1_1_index = edge_idx(cell_1_index, cell_1_block, 1)
          edge_1_2_index = edge_idx(cell_1_index, cell_1_block, 2)
          edge_1_3_index = edge_idx(cell_1_index, cell_1_block, 3)
          edge_2_1_index = edge_idx(cell_2_index, cell_2_block, 1)
          edge_2_2_index = edge_idx(cell_2_index, cell_2_block, 2)
          edge_2_3_index = edge_idx(cell_2_index, cell_2_block, 3)
          edge_1_1_block = edge_blk(cell_1_index, cell_1_block, 1)
          edge_1_2_block = edge_blk(cell_1_index, cell_1_block, 2)
          edge_1_3_block = edge_blk(cell_1_index, cell_1_block, 3)
          edge_2_1_block = edge_blk(cell_2_index, cell_2_block, 1)
          edge_2_2_block = edge_blk(cell_2_index, cell_2_block, 2)
          edge_2_3_block = edge_blk(cell_2_index, cell_2_block, 3)

          out_vn_e(je,blockNo) = &
            & vn_e(edge_1_1_index, edge_1_1_block) * coeffs(1, je, blockNo) + &
            & vn_e(edge_1_2_index, edge_1_2_block) * coeffs(2, je, blockNo) + &
            & vn_e(edge_1_3_index, edge_1_3_block) * coeffs(3, je, blockNo) + &
            & vn_e(edge_2_1_index, edge_2_1_block) * coeffs(4, je, blockNo) + &
            & vn_e(edge_2_2_index, edge_2_2_block) * coeffs(5, je, blockNo) + &
            & vn_e(edge_2_3_index, edge_2_3_block) * coeffs(6, je, blockNo)

        END DO
      END DO
      !$ACC END PARALLEL
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)
    !$ACC END DATA

!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2edges_viacell_2D_constZ_onTriangles
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_2d_constZ_onTriangles_sp( patch_3d, in_vn_e, operators_coefficients, out_vn_e )!&   subset_range)

    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(sp), INTENT(in)                       :: in_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_solvercoeff_singleprecision), INTENT(in)  :: operators_coefficients
    REAL(sp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, blockNo, start_edge_index, end_edge_index

    REAL(sp), POINTER :: all_coeffs(:,:,:)

    TYPE(t_subset_range), POINTER :: edges_indomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
!     IF (no_primal_edges /= 3) &
!       & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_indomain    => patch_2d%edges%in_domain
    all_coeffs        => operators_coefficients%edge2edge_viacell_coeff_all

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, cell_1_index, cell_1_block, &
!ICON_OMP  cell_2_index, cell_2_block, edge_1_1_index, edge_1_2_index, edge_1_3_index, &
!ICON_OMP  edge_1_1_block, edge_1_2_block, edge_1_3_block, edge_2_1_index, edge_2_2_index, &
!ICON_OMP  edge_2_3_index, edge_2_1_block, edge_2_2_block, edge_2_3_block)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_indomain%start_block, edges_indomain%end_block
      CALL get_index_range(edges_indomain, blockNo, start_edge_index, end_edge_index)

      DO je = start_edge_index, end_edge_index

        out_vn_e(je,blockNo) = 0.0_wp

        ! get the two cells of the edge
        cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
        cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)

        ! get the six edges of the two cells
        edge_1_1_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
        edge_1_2_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
        edge_1_3_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
        edge_2_1_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
        edge_2_2_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
        edge_2_3_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
        edge_1_1_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
        edge_1_2_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
        edge_1_3_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
        edge_2_1_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
        edge_2_2_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
        edge_2_3_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)

        out_vn_e(je,blockNo) = &
          & in_vn_e(edge_1_1_index, edge_1_1_block) * all_coeffs(1, je, blockNo) + &
          & in_vn_e(edge_1_2_index, edge_1_2_block) * all_coeffs(2, je, blockNo) + &
          & in_vn_e(edge_1_3_index, edge_1_3_block) * all_coeffs(3, je, blockNo) + &
          & in_vn_e(edge_2_1_index, edge_2_1_block) * all_coeffs(4, je, blockNo) + &
          & in_vn_e(edge_2_2_index, edge_2_2_block) * all_coeffs(5, je, blockNo) + &
          & in_vn_e(edge_2_3_index, edge_2_3_block) * all_coeffs(6, je, blockNo)

      END DO
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2edges_viacell_2d_constZ_onTriangles_sp
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !-map_edges2edges modified to use with zstar
  SUBROUTINE map_edges2edges_viacell_2D_zstar( patch_3d, in_vn_e, operators_coefficients, stretch_e, out_vn_e, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                       :: in_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(in)                       :: stretch_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    LOGICAL, INTENT(in), OPTIONAL              :: lacc

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, blockNo, start_edge_index, end_edge_index
    LOGICAL :: lzacc
    
    REAL(wp), POINTER :: all_coeffs(:,:,:)

    TYPE(t_subset_range), POINTER :: edges_indomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
!     IF (no_primal_edges /= 3) &
!       & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    !-----------------------------------------------------------------------
    patch_2d          => patch_3d%p_patch_2d(1)
    edges_indomain    => patch_2d%edges%in_domain
    all_coeffs        => operators_coefficients%edge2edge_viacell_coeff_all

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(patch_3d%p_patch_2d(1)%nblks_e) &
    !$ACC   COPYIN(patch_2d%cells%edge_idx, patch_2d%cells%edge_blk, patch_3d%p_patch_1d(1)%dolic_e) &
    !$ACC   COPYIN(patch_2d%edges%cell_idx, patch_2d%edges%cell_blk, all_coeffs, in_vn_e, stretch_e) &
    !$ACC   COPY(out_vn_e) IF(lzacc)

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, cell_1_index, cell_1_block, &
!ICON_OMP  cell_2_index, cell_2_block, edge_1_1_index, edge_1_2_index, edge_1_3_index, &
!ICON_OMP  edge_1_1_block, edge_1_2_block, edge_1_3_block, edge_2_1_index, edge_2_2_index, &
!ICON_OMP  edge_2_3_index, edge_2_1_block, edge_2_2_block, edge_2_3_block)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_indomain%start_block, edges_indomain%end_block
      CALL get_index_range(edges_indomain, blockNo, start_edge_index, end_edge_index)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO je = start_edge_index, end_edge_index

        IF (patch_3d%p_patch_1d(1)%dolic_e(je,blockNo) < 1) CYCLE

        out_vn_e(je,blockNo) = 0.0_wp

        ! get the two cells of the edge
        cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
        cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)

        ! get the six edges of the two cells
        edge_1_1_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
        edge_1_2_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
        edge_1_3_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
        edge_2_1_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
        edge_2_2_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
        edge_2_3_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
        edge_1_1_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
        edge_1_2_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
        edge_1_3_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
        edge_2_1_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
        edge_2_2_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
        edge_2_3_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)
!          if ( stretch_e(edge_1_1_index, edge_1_1_block ) .NE. 1.0_wp   ) THEN
!            write(*, *) edge_1_1_index, edge_1_1_block, stretch_e(edge_1_1_index, edge_1_1_block)
!          ENDIF
        out_vn_e(je,blockNo) = &
          & in_vn_e(edge_1_1_index, edge_1_1_block) * all_coeffs(1, je, blockNo) &
          & * stretch_e(edge_1_1_index, edge_1_1_block) + &
          & in_vn_e(edge_1_2_index, edge_1_2_block) * all_coeffs(2, je, blockNo) &
          & * stretch_e(edge_1_2_index, edge_1_2_block) + &
          & in_vn_e(edge_1_3_index, edge_1_3_block) * all_coeffs(3, je, blockNo) &
          & * stretch_e(edge_1_3_index, edge_1_3_block) + &
          & in_vn_e(edge_2_1_index, edge_2_1_block) * all_coeffs(4, je, blockNo) &
          & * stretch_e(edge_2_1_index, edge_2_1_block) + &
          & in_vn_e(edge_2_2_index, edge_2_2_block) * all_coeffs(5, je, blockNo) &
          & * stretch_e(edge_2_2_index, edge_2_2_block) + &
          & in_vn_e(edge_2_3_index, edge_2_3_block) * all_coeffs(6, je, blockNo) &
          & * stretch_e(edge_2_3_index, edge_2_3_block)

      END DO
      !$ACC END PARALLEL LOOP
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    !$ACC END DATA
  END SUBROUTINE map_edges2edges_viacell_2D_zstar
  !-----------------------------------------------------------------------------




  !-----------------------------------------------------------------------------
!<Optimize:inUse>
  !This sbr is used in left-hand-side of free-surface equation for certain solver configurations.
  SUBROUTINE map_edges2edges_viacell_2D_per_level( patch_3d, vn_e, operators_coefficients, out_vn_e, level)!&   subset_range)

    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER, INTENT(IN)                        :: level
    !Local variables
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, blockNo, ie
    REAL(wp) :: thick_edge, thick_cell
    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER :: patch_2d

    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
!     IF ( patch_2d%cells%max_connectivity == 3 .and. fast_performance_level > 10 ) THEN
!       CALL map_edges2edges_viacell_2D_per_level_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e, level)!&   subset_range)
!       RETURN
!     ENDIF
    !-----------------------------------------------------------------------

    edges_inDomain => patch_2d%edges%in_domain

    DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
      CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)

      out_vn_e(:,blockNo) = 0.0_wp

      edge_idx_loop: DO je =  start_edge_index, end_edge_index

          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN

            ictr = 0
            il_c = patch_2d%edges%cell_idx(je,blockNo,1)
            ib_c = patch_2d%edges%cell_blk(je,blockNo,1)
            thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)


            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
              !IF(il_e==je.and.ib_e==blockNo)THEN
              out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
                & +vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                & *(thick_edge/thick_cell)
              !ENDIF
            END DO

            ictr = no_primal_edges
            il_c = patch_2d%edges%cell_idx(je,blockNo,2)
            ib_c = patch_2d%edges%cell_blk(je,blockNo,2)
            thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
              !IF(il_e==je.and.ib_e==blockNo)THEN
              out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
                & +vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                & *(thick_edge/thick_cell)
              !ENDIF
            END DO
          ENDIF
      END DO edge_idx_loop
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block

  END SUBROUTINE map_edges2edges_viacell_2D_per_level
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_2D_per_level_onTriangles( patch_3d, in_vn_e, operators_coefficients, out_vn_e, level)!&   subset_range)

    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                       :: in_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER, INTENT(IN)                        :: level

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, blockNo, start_edge_index, end_edge_index

    REAL(wp), POINTER :: all_coeffs(:,:,:)
    TYPE(t_subset_range), POINTER :: edges_indomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_indomain    => patch_2d%edges%in_domain
    all_coeffs        => operators_coefficients%edge2edge_viacell_coeff_all

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, cell_1_index, cell_1_block, &
!ICON_OMP   cell_2_index, cell_2_block, edge_1_1_index, edge_1_2_index, edge_1_3_index, &
!ICON_OMP  edge_1_1_block, edge_1_2_block, edge_1_3_block, edge_2_1_index, edge_2_2_index, &
!ICON_OMP  edge_2_3_index, edge_2_1_block, edge_2_2_block, edge_2_3_block)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_indomain%start_block, edges_indomain%end_block
      CALL get_index_range(edges_indomain, blockNo, start_edge_index, end_edge_index)

      DO je = start_edge_index, end_edge_index

        out_vn_e(je,blockNo) = 0.0_wp

        IF (patch_3d%lsm_e(je,1,blockNo) == sea) THEN ! this if should be removed

          ! get the two cells of the edge
          cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
          cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
          cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
          cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)

          ! get the six edges of the two cells
          edge_1_1_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
          edge_1_2_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
          edge_1_3_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
          edge_2_1_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
          edge_2_2_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
          edge_2_3_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
          edge_1_1_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
          edge_1_2_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
          edge_1_3_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
          edge_2_1_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
          edge_2_2_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
          edge_2_3_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)
          !not that all_coeffs already includes the edge thickness
          out_vn_e(je,blockNo) = &
            & in_vn_e(edge_1_1_index, edge_1_1_block) * all_coeffs(1, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_1_index,level, edge_1_1_block)&
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block)     &
            &+ &
            & in_vn_e(edge_1_2_index, edge_1_2_block) * all_coeffs(2, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_2_index,level, edge_1_2_block)&
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block)     &
            &+ &
            & in_vn_e(edge_1_3_index, edge_1_3_block) * all_coeffs(3, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_3_index,level, edge_1_3_block)&
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block)     &
            &+ &
            & in_vn_e(edge_2_1_index, edge_2_1_block) * all_coeffs(4, je, blockNo)       &
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block)     &
            &+ &
            & in_vn_e(edge_2_2_index, edge_2_2_block) * all_coeffs(5, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_2_index,level, edge_2_2_block)&
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block)     &
            &+ &
            & in_vn_e(edge_2_3_index, edge_2_3_block) * all_coeffs(6, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_3_index,level, edge_2_3_block)&
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block)

        ENDIF
      END DO
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2edges_viacell_2D_per_level_onTriangles
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_2D_per_level_onTriangles_sp( patch_3d, in_vn_e, operators_coefficients, out_vn_e, level)!&   subset_range)

    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(sp), INTENT(in)                       :: in_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_solvercoeff_singleprecision), INTENT(in)  :: operators_coefficients
    REAL(sp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER, INTENT(IN)                        :: level

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, blockNo, start_edge_index, end_edge_index

    REAL(sp), POINTER :: all_coeffs(:,:,:)
    !REAL(wp) :: thick_edge, thick_cell
    TYPE(t_subset_range), POINTER :: edges_indomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_indomain    => patch_2d%edges%in_domain
    all_coeffs        => operators_coefficients%edge2edge_viacell_coeff_all

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, cell_1_index, cell_1_block, &
!ICON_OMP   cell_2_index, cell_2_block, edge_1_1_index, edge_1_2_index, edge_1_3_index, &
!ICON_OMP  edge_1_1_block, edge_1_2_block, edge_1_3_block, edge_2_1_index, edge_2_2_index, &
!ICON_OMP  edge_2_3_index, edge_2_1_block, edge_2_2_block, edge_2_3_block)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_indomain%start_block, edges_indomain%end_block
      CALL get_index_range(edges_indomain, blockNo, start_edge_index, end_edge_index)

      DO je = start_edge_index, end_edge_index

        out_vn_e(je,blockNo) = 0.0_wp

        IF (patch_3d%lsm_e(je,1,blockNo) == sea) THEN ! this if should be removed

          ! get the two cells of the edge
          cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
          cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
          cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
          cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)

          ! get the six edges of the two cells
          edge_1_1_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
          edge_1_2_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
          edge_1_3_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
          edge_2_1_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
          edge_2_2_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
          edge_2_3_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
          edge_1_1_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
          edge_1_2_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
          edge_1_3_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
          edge_2_1_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
          edge_2_2_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
          edge_2_3_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)
          !not that all_coeffs already includes the edge thickness
          out_vn_e(je,blockNo) = &
            & in_vn_e(edge_1_1_index, edge_1_1_block) * all_coeffs(1, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_1_index,level, edge_1_1_block)&
            &* (patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block))     &
            &+ &
            & in_vn_e(edge_1_2_index, edge_1_2_block) * all_coeffs(2, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_2_index,level, edge_1_2_block)&
            & * (patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block))     &
            &+ &
            & in_vn_e(edge_1_3_index, edge_1_3_block) * all_coeffs(3, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_3_index,level, edge_1_3_block)&
            & * (patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block))     &
            &+ &
            & in_vn_e(edge_2_1_index, edge_2_1_block) * all_coeffs(4, je, blockNo)       &
            & * (patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block))     &
            &+ &
            & in_vn_e(edge_2_2_index, edge_2_2_block) * all_coeffs(5, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_2_index,level, edge_2_2_block)&
            &/(patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block))     &
            &+ &
            & in_vn_e(edge_2_3_index, edge_2_3_block) * all_coeffs(6, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_3_index,level, edge_2_3_block)&
            & * (patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block))

        ENDIF
      END DO
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2edges_viacell_2D_per_level_onTriangles_sp
  !-----------------------------------------------------------------------------



  !-----------------------------------------------------------------------------
  SUBROUTINE map_edges2edges_viavert_3d(patch_3d, vn, p_vn_dual,operators_coefficients, vort_flux)
    TYPE(t_patch_3d ),TARGET,INTENT(in)        :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(inout)  :: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in)          :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    !Local variables
    INTEGER :: startLevel, endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_inDomain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    startLevel    = 1
    endLevel    = n_zlev

    DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
      CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)

      level_loop: DO level = startLevel, endLevel

        vort_flux(:,:,blockNo) = 0.0_wp
        edge_idx_loop: DO je =  start_edge_index, end_edge_index

          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN

            DO neighbor=1,2
              IF(neighbor==1) ictr = 0
              IF(neighbor==2) ictr = no_dual_edges

              il_v = patch_2d%edges%vertex_idx(je,blockNo,neighbor)
              ib_v = patch_2d%edges%vertex_blk(je,blockNo,neighbor)

              DO vertex_edge=1, patch_2d%verts%num_edges(il_v,ib_v)!no_dual_cell_edges

                ictr =ictr+1

                il_e = patch_2d%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = patch_2d%verts%edge_blk(il_v,ib_v,vertex_edge)

                vort_flux(je,level,blockNo) =  vort_flux(je,level,blockNo)+vn(il_e,level,ib_e)&
                  & *operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)
              END DO
            END DO

          ELSE
            vort_flux(je,level,blockNo)= 0.0_wp
          ENDIF ! (v_base%lsm_e(je,level,blockNo) <= sea_boundary)
        END DO edge_idx_loop
      END DO level_loop
    END DO ! blockNo = edges_inDomain%start_block, edges_inDomain%end_block

  END SUBROUTINE map_edges2edges_viavert_3d
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !<Optimize:inUse>
  SUBROUTINE map_cell2edges_3d_mlevels( patch_3d, p_vn_c, ptp_vn, operators_coefficients,&
    & opt_startLevel, opt_endLevel, subset_range, lacc )

    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:,:)    ! input vector (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)                      :: ptp_vn(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)                     :: operators_coefficients
    INTEGER, INTENT(in), OPTIONAL :: opt_startLevel        ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_endLevel        ! optional vertical end level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, blockNo, level
    INTEGER :: cell_1_index,cell_1_block, cell_2_index,cell_2_block
    LOGICAL :: lzacc
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    edges_in_domain   => patch_2d%edges%in_domain


    ! check optional arguments
    IF ( PRESENT(opt_startLevel) ) THEN
      startLevel = opt_startLevel
    ELSE
      startLevel = 1
    END IF
    IF ( PRESENT(opt_endLevel) ) THEN
      endLevel = opt_endLevel
    ELSE
      endLevel = n_zlev
    END IF

    CALL set_acc_host_or_device(lzacc, lacc)

    ! calculation of transposed P^TPv from Pv (incart coord)
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index,je, level,  &
!ICON_OMP cell_1_index, cell_1_block, cell_2_index, cell_2_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      ptp_vn(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO je = start_edge_index, end_edge_index
        !Get indices of two adjacent triangles
        cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
        cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)
        DO level = 1, MIN(endlevel,patch_3d%p_patch_1d(1)%dolic_e(je,blockNo))

          ptp_vn(je,level,blockNo) =&
            & DOT_PRODUCT(p_vn_c(cell_1_index,level,cell_1_block)%x,&
            & operators_coefficients%edge2cell_coeff_cc_t(je,level,blockNo,1)%x)&
            & +DOT_PRODUCT(p_vn_c(cell_2_index,level,cell_2_block)%x,&
            & operators_coefficients%edge2cell_coeff_cc_t(je,level,blockNo,2)%x)

        END DO
      END DO
      !$ACC END PARALLEL LOOP

    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    ! sync the result if necessary
    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2d, ptp_vn)
    ENDIF

  END SUBROUTINE map_cell2edges_3d_mlevels
  !-----------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !<Optimize:inUse:done>
  SUBROUTINE map_cell2edges_3d_1level( patch_3d, p_vn_c, ptp_vn,operators_coefficients, level, subset_range, lacc )

    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:)    ! input vector (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)                      :: ptp_vn(:,:)    ! output vector (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)                     :: operators_coefficients
    INTEGER, INTENT(in) :: level          ! vertical level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    !Local variables
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, blockNo
    INTEGER :: cell_1_index,cell_1_block, cell_2_index,cell_2_block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    edges_in_domain   => patch_2d%edges%in_domain

    ! calculation of transposed P^TPv from Pv (incart coord)
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, &
!ICON_OMP  cell_1_index,cell_1_block, cell_2_index,cell_2_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      ptp_vn(:,blockNo) = 0.0_wp
      !$ACC END KERNELS
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO je =  start_edge_index, end_edge_index

        IF (patch_3d%p_patch_1d(1)%dolic_e(je, blockNo) > 0) THEN
          !Get indices of two adjacent triangles
          cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
          cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
          cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
          cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)
          ptp_vn(je,blockNo) = &
            & DOT_PRODUCT(p_vn_c(cell_1_index,cell_1_block)%x,               &
            & operators_coefficients%edge2cell_coeff_cc_t(je,level,blockNo,1)%x)  &
            & +DOT_PRODUCT(p_vn_c(cell_2_index,cell_2_block)%x,              &
            & operators_coefficients%edge2cell_coeff_cc_t(je,level,blockNo,2)%x)

        ENDIF

      END DO
      !$ACC END PARALLEL LOOP
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    ! sync the result if necessary
    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2d, ptp_vn)
    ENDIF

  END SUBROUTINE map_cell2edges_3d_1level
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE maps for a fluid column a scalar value from the top/bottom of a 3D prism to the central level of the prism.
  !!
  !!
  !pgi$r opt 1
  SUBROUTINE map_vec_prismtop2center_on_block(patch_3d, vec_top, vec_center, &
    & blockNo, start_cell_index, end_cell_index, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(in)            :: patch_3d
    TYPE(t_cartesian_coordinates), INTENT(in)       :: vec_top(:,:) ! (nproma, n_zlev+1)
    TYPE(t_cartesian_coordinates), INTENT(inout)    :: vec_center(:,:) ! (nproma, n_zlev) ! out
    INTEGER, INTENT(in)                             :: blockNo, start_cell_index, end_cell_index
    LOGICAL, INTENT(IN), OPTIONAL                   :: lacc

    !Local variables
    INTEGER :: level, jc!,jb
    INTEGER :: start_level, end_level
    LOGICAL :: lzacc
    REAL(wp), POINTER :: prism_center_distance(:,:), prism_thick(:,:)
    !-------------------------------------------------------------------------------
    start_level = 1
!     vec_center(1:nproma,1:n_zlev,blockNo)%x(1)=0.0_wp
!     vec_center(1:nproma,1:n_zlev,blockNo)%x(2)=0.0_wp
!     vec_center(1:nproma,1:n_zlev,blockNo)%x(3)=0.0_wp
   !-------------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! these do not include the height
    prism_center_distance => patch_3D%p_patch_1D(1)%constantPrismCenters_Zdistance  (:,:,blockNo)
    prism_thick           => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,blockNo)

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO jc = start_cell_index, end_cell_index
      end_level  = patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)
!       IF ( end_level >=min_dolic ) THEN
        !$ACC LOOP SEQ
        DO level = start_level, end_level
          vec_center(jc,level)%x &
          & = (prism_center_distance(jc,level)   * vec_top(jc,level)%x    &
          & +  prism_center_distance(jc,level+1) * vec_top(jc,level+1)%x) &
          & / (2.0_wp*prism_thick(jc,level))

        END DO

    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
!     CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vec_center(:,:,:)%x(1))
!     CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vec_center(:,:,:)%x(2))
!     CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vec_center(:,:,:)%x(3))
  END SUBROUTINE map_vec_prismtop2center_on_block
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE maps for a fluid column a scalar value from the top/bottom of a 3D prism to the central level of the prism.
  !!
  !!
  SUBROUTINE map_scalar_prismtop2center(patch_3d, scalar_top, p_op_coeff,scalar_center)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    REAL(wp)                                         :: scalar_top(nproma, n_zlev+1,patch_3D%p_patch_2d(1)%nblks_c)
    TYPE(t_operator_coeff),INTENT(in)                :: p_op_coeff
    REAL(wp)                                         :: scalar_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_c)

    !Local variables
    INTEGER :: level, blockNo, jc!,jb
    INTEGER :: start_cell_index, end_cell_index!, cell_index
    !INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level!, level
    TYPE(t_subset_range), POINTER :: cells_in_domain!, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), POINTER ::  prism_center_distance(:,:),prism_thick(:,:)
    ! INTEGER :: dolic
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    !edges_in_domain => patch_2D%edges%in_domain
    start_level = 1

   !-------------------------------------------------------------------------------
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, level,&
!ICON_OMP prism_center_distance,prism_thick, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

      ! these do not include the height
      ! if they are used for the GM-Redi then we should re-consider if height is needed
      prism_center_distance => patch_3D%p_patch_1D(1)%constantPrismCenters_Zdistance  (:,:,blockNo)
      prism_thick => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,blockNo)

      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)-1
          scalar_center(jc,level,blockNo) &
          & = (prism_center_distance(jc,level)   * scalar_top(jc,level,blockNo)    &
          & +  prism_center_distance(jc,level+1) * scalar_top(jc,level+1,blockNo)) &
          & / (2.0_wp*prism_thick(jc,level))
!write(123,*)'details',level,prism_center_distance(jc,level),&
!& scalar_top(jc,level,blockNo),&
!          &prism_center_distance(jc,level+1) , scalar_top(jc,level+1,blockNo),prism_thick(jc,level)
        END DO

      END DO
    END DO

!ICON_OMP_END_PARALLEL_DO
   !CALL sync_patch_array(sync_c, patch_2D, scalar_center)
  END SUBROUTINE map_scalar_prismtop2center
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE map_scalar_prismtop2center_onBlock(patch_3d, scalar_top, scalar_center, &
    & blockNo, start_cell_index, end_cell_index)
    TYPE(t_patch_3d ),TARGET                         :: patch_3d
    REAL(wp)                                         :: scalar_top(:,:) ! (nproma, n_zlev+1)
    REAL(wp)                                         :: scalar_center(:,:) !nproma, n_zlev)
    INTEGER, INTENT(in) :: blockNo, start_cell_index, end_cell_index

    !Local variables
    INTEGER :: level, jc!,jb
    REAL(wp), POINTER ::  prism_center_distance(:,:),prism_thick(:,:)

    !-------------------------------------------------------------------------------

    ! these do not include the height
    ! if they are used for the GM-Redi then we should re-consider if height is needed
    prism_center_distance => patch_3D%p_patch_1D(1)%constantPrismCenters_Zdistance  (:,:,blockNo)
    prism_thick => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,blockNo)

    DO jc = start_cell_index, end_cell_index
      DO level = 1, patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)
        scalar_center(jc,level) &
          & = (prism_center_distance(jc,level)   * scalar_top(jc,level)    &
          & +  prism_center_distance(jc,level+1) * scalar_top(jc,level+1)) &
          & / (2.0_wp*prism_thick(jc,level))

      END DO
    END DO

  END SUBROUTINE map_scalar_prismtop2center_onBlock
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE maps within a fluid column a scalar value from the the central level of the prism to top/bottom of a 3D prism.
  !!
  !!
  SUBROUTINE map_scalar_center2prismtop(patch_3d, scalar_center, p_op_coeff,scalar_top)
    TYPE(t_patch_3d ),TARGET, INTENT(in   )          :: patch_3d
    REAL(wp), INTENT(in)                             :: scalar_center(:,:,:)
    TYPE(t_operator_coeff),INTENT(in)                :: p_op_coeff
    REAL(wp), INTENT(inout)                          :: scalar_top(:,:,:)

    !Local variables
    INTEGER :: blockNo
    INTEGER :: start_cell_index, end_cell_index, cell_index
    !INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level
    TYPE(t_subset_range), POINTER :: cells_in_domain!, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !REAL(wp), POINTER ::  prism_center_distance(:,:),prism_thick(:,:)
!     INTEGER :: dolic
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    !edges_in_domain => patch_2D%edges%in_domain
    start_level = 1

   !-------------------------------------------------------------------------------
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, level)
!ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

      scalar_top(:,:,blockNo) = 0.0_wp
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      DO cell_index = start_cell_index, end_cell_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(cell_index,blockNo)-1
          scalar_top(cell_index,level+1,blockNo) &
          & = 0.5_wp*( scalar_center(cell_index,level,blockNo)    &
          & +          scalar_center(cell_index,level+1,blockNo))
        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_scalar_center2prismtop
  !-------------------------------------------------------------------------





  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE maps for a fluid column a scalar value from the top/bottom of a 3D prism to the central level of the prism.
  !!
  !!
  SUBROUTINE map_vec_prismtop2center_on_block_GM(patch_3d, vec_top, vec_center, &
    & blockNo, start_cell_index, end_cell_index)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_cartesian_coordinates), INTENT(in)        :: vec_top(nproma, n_zlev)
    TYPE(t_cartesian_coordinates), INTENT(inout)     :: vec_center(nproma, n_zlev) ! out
!     TYPE(t_operator_coeff), INTENT(in)               :: p_op_coeff
    INTEGER, INTENT(in)                              :: blockNo, start_cell_index, end_cell_index

    !Local variables
    INTEGER :: level, jc!,jb
    INTEGER :: start_level, end_level
    REAL(wp), POINTER :: prism_center_distance(:,:), prism_thick(:,:)
    !-------------------------------------------------------------------------------
    start_level = 1
!     vec_center(1:nproma,1:n_zlev,blockNo)%x(1)=0.0_wp
!     vec_center(1:nproma,1:n_zlev,blockNo)%x(2)=0.0_wp
!     vec_center(1:nproma,1:n_zlev,blockNo)%x(3)=0.0_wp
   !-------------------------------------------------------------------------------
    ! these do not include the height
    !prism_center_distance => patch_3D%p_patch_1D(1)%constantPrismCenters_Zdistance  (:,:,blockNo)
    !prism_thick           => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,blockNo)

    DO jc = start_cell_index, end_cell_index
      end_level  = patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)

        DO level = start_level, end_level-1
          vec_center(jc,level)%x &
          & = 0.5_wp*( vec_top(jc,level)%x + vec_top(jc,level+1)%x)
        END DO

    END DO
!     CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vec_center(:,:,:)%x(1))
!     CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vec_center(:,:,:)%x(2))
!     CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vec_center(:,:,:)%x(3))
  END SUBROUTINE map_vec_prismtop2center_on_block_GM
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE maps for a fluid column a scalar value from the top/bottom of a 3D prism to the central level of the prism.
  !!
  !!
  SUBROUTINE map_scalar_prismtop2center_GM(patch_3d, scalar_top, p_op_coeff, scalar_center)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    REAL(wp)                                         :: scalar_top(nproma, n_zlev+1,patch_3D%p_patch_2d(1)%nblks_c)
    TYPE(t_operator_coeff),INTENT(in)                :: p_op_coeff
    REAL(wp)                                         :: scalar_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_c)

    !Local variables
    INTEGER :: level, blockNo, jc!,jb
    INTEGER :: start_cell_index, end_cell_index!, cell_index
    !INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level!, level
    TYPE(t_subset_range), POINTER :: cells_in_domain!, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !REAL(wp), POINTER ::  prism_center_distance(:,:),prism_thick(:,:)
    ! INTEGER :: dolic
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    !edges_in_domain => patch_2D%edges%in_domain
    start_level = 1

   !-------------------------------------------------------------------------------
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index)
!ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

      ! these do not include the height
      ! if they are used for the GM-Redi then we should re-consider if height is needed
      !prism_center_distance => patch_3D%p_patch_1D(1)%constantPrismCenters_Zdistance  (:,:,blockNo)
      !prism_thick => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,blockNo)

      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)-1
          scalar_center(jc,level,blockNo) &
          & = 0.5_wp*( scalar_top(jc,level,blockNo) + scalar_top(jc,level+1,blockNo))

        END DO

      END DO
    END DO

!ICON_OMP_END_PARALLEL_DO
   !CALL sync_patch_array(sync_c, patch_2D, scalar_center)
  END SUBROUTINE map_scalar_prismtop2center_GM
  !------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE maps within a fluid column a scalar value from the the central level of the prism to top/bottom of a 3D prism.
  !!
  !!
  SUBROUTINE map_scalar_center2prismtop_GM(patch_3d, scalar_center, p_op_coeff, scalar_top)
    TYPE(t_patch_3d ),TARGET, INTENT(in   )          :: patch_3d
    REAL(wp), INTENT(inout)                          :: scalar_center(:,:,:)
    TYPE(t_operator_coeff),INTENT(in)                :: p_op_coeff
    REAL(wp), INTENT(inout)                          :: scalar_top(:,:,:)

    !Local variables
    INTEGER :: blockNo
    INTEGER :: start_cell_index, end_cell_index, cell_index
    !INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level
    TYPE(t_subset_range), POINTER :: cells_in_domain!, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), POINTER ::  prism_center_distance(:,:),prism_thick(:,:)
!     INTEGER :: dolic
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    !edges_in_domain => patch_2D%edges%in_domain
    start_level = 1

   !-------------------------------------------------------------------------------
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, level,&
!ICON_OMP prism_center_distance,prism_thick) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      prism_center_distance => patch_3D%p_patch_1D(1)%constantPrismCenters_Zdistance(:,:,blockNo)
      prism_thick => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,blockNo)

      scalar_top(:,:,blockNo) = 0.0_wp
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      DO cell_index = start_cell_index, end_cell_index
!        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(cell_index,blockNo)-1
!        scalar_center(cell_index,level,blockNo)=level
!        END DO

        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(cell_index,blockNo)-1
          scalar_top(cell_index,level+1,blockNo) &
          & = ( prism_thick(cell_index,level)  * scalar_center(cell_index,level,blockNo)      &
          & +   prism_thick(cell_index,level+1)* scalar_center(cell_index,level+1,blockNo))&
          & /(2.0_wp*prism_center_distance(cell_index,level))
!write(123,*)'details',level, scalar_top(cell_index,level+1,blockNo),&
!&scalar_center(cell_index,level,blockNo),prism_thick(cell_index,level),&
!&scalar_center(cell_index,level+1,blockNo),prism_thick(cell_index,level+1),&
!!&(prism_thick(cell_index,level)+prism_thick(cell_index,level+1)),&
!2.0_wp*prism_center_distance(cell_index,level)
        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO
!stop
  END SUBROUTINE map_scalar_center2prismtop_GM
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE map_vector_center2prismtop_onBlock(patch_3d, vector_center, vector_top, &
    & blockNo, start_cell_index, end_cell_index)
    TYPE(t_patch_3d ),TARGET          :: patch_3d
    TYPE(t_cartesian_coordinates)     :: vector_center(:,:) ! (nproma, n_zlev,blockNo) ! out
    TYPE(t_cartesian_coordinates)     :: vector_top(:,:)
    INTEGER, INTENT(in) :: blockNo, start_cell_index, end_cell_index

    !Local variables
    INTEGER :: cell_index, level
    REAL(wp), POINTER ::  inv_prism_center_distance(:,:),prism_thick(:,:)

    !-------------------------------------------------------------------------------
    inv_prism_center_distance => patch_3D%p_patch_1D(1)%constantPrismCenters_invZdistance  (:,:,blockNo)
    prism_thick => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,blockNo)
    !-------------------------------------------------------------------------------
    vector_top(:,:)%x(1) = 0.0_wp
    vector_top(:,:)%x(2) = 0.0_wp
    vector_top(:,:)%x(3) = 0.0_wp
    DO cell_index = start_cell_index, end_cell_index
      ! the top moves only up-down
!       vector_top(cell_index,1)%x = vector_center(cell_index,1)%x !* 0.5_wp
      vector_top(cell_index,1)%x = 0.0_wp

      DO level = 2, patch_3D%p_patch_1d(1)%dolic_c(cell_index,blockNo)
        vector_top(cell_index,level)%x = &
        & (  vector_center(cell_index,level-1)%x * prism_thick(cell_index,level-1)  &
        &  + vector_center(cell_index,level)%x   * prism_thick(cell_index,level))   &
        &  * 2.0_wp * inv_prism_center_distance(cell_index,level)
      END DO
    END DO

  END SUBROUTINE map_vector_center2prismtop_onBlock
  !-------------------------------------------------------------------------


  ! !   !-----------------------------------------------------------------------------
!   !-----------------------------------------------------------------------------
!   SUBROUTINE map_edges2edges_viacell_3d_mlev( patch_3d, vn_e, operators_coefficients, out_vn_e,scalar, &
!     & opt_startLevel, opt_endLevel, subset_range)
!
!     TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
!     REAL(wp), INTENT(in)                       :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
!     TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
!     REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
!     REAL(wp), INTENT(in), OPTIONAL :: scalar(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     INTEGER, INTENT(in), OPTIONAL :: opt_startLevel       ! optional vertical start level
!     INTEGER, INTENT(in), OPTIONAL :: opt_endLevel       ! optional vertical end level
!     TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
!     !Local variables
!     INTEGER :: startLevel, endLevel
!     INTEGER :: start_edge_index, end_edge_index
!     INTEGER :: il_e, ib_e, il_c, ib_c, ictr!, neighbor
!     INTEGER :: je, blockNo, level, ie
!     REAL(wp) :: scalar_cell
!     REAL(wp) :: thick_edge, thick_cell, thick_frac
!     TYPE(t_subset_range), POINTER :: edges_inDomain
!     TYPE(t_patch), POINTER :: patch_2d
!     !-----------------------------------------------------------------------
!     patch_2d   => patch_3d%p_patch_2d(1)
!     !-----------------------------------------------------------------------
!     IF ( PRESENT(subset_range) ) THEN
!       edges_inDomain => subset_range
!     ELSE
!       edges_inDomain => patch_2d%edges%in_domain
!     ENDIF
!     !-----------------------------------------------------------------------
!     IF ( PRESENT(opt_startLevel) ) THEN
!       startLevel = opt_startLevel
!     ELSE
!       startLevel = 1
!     END IF
!     IF ( PRESENT(opt_endLevel) ) THEN
!       endLevel = opt_endLevel
!     ELSE
!       endLevel = n_zlev
!     END IF
!
!     IF(.NOT.PRESENT(scalar))THEN
!
!       DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!         CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
!         out_vn_e(:,:,blockNo)= 0.0_wp
!
!         level_loop_e: DO level = startLevel, endLevel
!           edge_idx_loop: DO je =  start_edge_index, end_edge_index
!             IF (patch_3d%lsm_e(je,level,blockNo) <= sea_boundary) THEN
!               ictr            = 0
!               !IF(patch_3D%lsm_e(je,level,blockNo) == sea)THEN
!               il_c = patch_2d%edges%cell_idx(je,blockNo,1)
!               ib_c = patch_2d%edges%cell_blk(je,blockNo,1)
!               thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
!
!               DO ie=1, no_primal_edges
!                 ictr =ictr+1
!                 il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!                 ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!
!                 thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
!                 thick_frac = thick_edge/thick_cell
!                 out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
!                   & +vn_e(il_e,level,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr))*thick_frac
!
!               END DO
!
!               ictr = no_primal_edges
!               il_c = patch_2d%edges%cell_idx(je,blockNo,2)
!               ib_c = patch_2d%edges%cell_blk(je,blockNo,2)
!               thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
!
!               DO ie=1, no_primal_edges
!                 ictr =ictr+1
!                 il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!                 ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!
!                 thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
!                 thick_frac = thick_edge/thick_cell
!                 out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
!                   & +vn_e(il_e,level,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr))*thick_frac
!
!               END DO
!             ENDIF
!           END DO edge_idx_loop
!         END DO level_loop_e
!       END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!
!     ELSEIF(PRESENT(scalar))THEN
!
!       DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!         CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
!
!         level_loop_e2: DO level = startLevel, endLevel
!           edge_idx_loop2: DO je =  start_edge_index, end_edge_index
!             IF (patch_3d%lsm_e(je,level,blockNo) <= sea_boundary) THEN
!               ictr            = 0
!               out_vn_e(je,level,blockNo)= 0.0_wp
!               il_c        = patch_2d%edges%cell_idx(je,blockNo,1)
!               ib_c        = patch_2d%edges%cell_blk(je,blockNo,1)
!               scalar_cell = scalar(il_c,level,ib_c)
!               thick_cell  = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
!
!               DO ie=1, no_primal_edges
!                 ictr =ictr+1
!                 il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!                 ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!                 thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
!                 thick_frac = thick_edge/thick_cell
!
!                 out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
!                   & +vn_e(il_e,level,ib_e)*scalar_cell   &
!                   & *(operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr))*thick_frac
!               END DO
!
!
!               ictr        = no_primal_edges
!               il_c        = patch_2d%edges%cell_idx(je,blockNo,2)
!               ib_c        = patch_2d%edges%cell_blk(je,blockNo,2)
!               scalar_cell = scalar(il_c,level,ib_c)
!               thick_cell  = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
!
!               DO ie=1, no_primal_edges
!                 ictr =ictr+1
!                 il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!                 ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!                 thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
!                 thick_frac = thick_edge/thick_cell
!
!                 out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
!                   & +vn_e(il_e,level,ib_e)*scalar_cell   &
!                   & *(operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr))*thick_frac
!
!               END DO
!
!             ENDIF
!           END DO edge_idx_loop2
!         END DO level_loop_e2
!       END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!     ENDIF
!
!   END SUBROUTINE map_edges2edges_viacell_3d_mlev
!   !-----------------------------------------------------------------------------
!
!   !-----------------------------------------------------------------------------
!   SUBROUTINE map_edges2edges_viacell_3d_1lev( patch_3d, vn_e, operators_coefficients, out_vn_e,scalar,scalar_e, level, &
!     & subset_range)
!
!     TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
!     REAL(wp), INTENT(in)                       :: vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
!     TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
!     REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
!     REAL(wp), INTENT(in), OPTIONAL :: scalar(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     REAL(wp), INTENT(in), OPTIONAL :: scalar_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
!     INTEGER, INTENT(in), OPTIONAL :: level       ! optional vertical start level
!     TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
!     !Local variables
!     INTEGER :: startLevel, endLevel
!     INTEGER :: start_edge_index, end_edge_index
!     INTEGER :: il_e, ib_e, il_c, ib_c, ictr
!     INTEGER :: je, blockNo, ie, lev
!     REAL(wp) :: scalar_cell
!     TYPE(t_subset_range), POINTER :: edges_inDomain
!     TYPE(t_patch), POINTER :: patch_2d
!     !-----------------------------------------------------------------------
!     patch_2d   => patch_3d%p_patch_2d(1)
!     !-----------------------------------------------------------------------
!     IF ( PRESENT(subset_range) ) THEN
!       edges_inDomain => subset_range
!     ELSE
!       edges_inDomain => patch_2d%edges%in_domain
!     ENDIF
!     !--------------------------
!     IF ( PRESENT(level) ) THEN
!       lev = level
!     ELSE
!       lev = 1
!     END IF
!
!     IF(.NOT.PRESENT(scalar))THEN
!
!       DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!         CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
!
!         out_vn_e(:,blockNo) = 0.0_wp
!         edge_idx_loop: DO je =  start_edge_index, end_edge_index
!           IF (patch_3d%lsm_e(je,lev,blockNo) <= sea_boundary) THEN
!             ictr          = 0
!             il_c = patch_2d%edges%cell_idx(je,blockNo,1)
!             ib_c = patch_2d%edges%cell_blk(je,blockNo,1)
!
!             !thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)
!             DO ie=1, no_primal_edges
!               ictr =ictr+1
!               il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!               ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!               !thick_edge=patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
!               !thick_frac=thick_edge/thick_cell
!               out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
!                 & +vn_e(il_e,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,lev,blockNo,ictr))!*thick_frac
!
!             END DO
!
!             ictr = no_primal_edges
!             il_c = patch_2d%edges%cell_idx(je,blockNo,2)
!             ib_c = patch_2d%edges%cell_blk(je,blockNo,2)
!
!             !thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)
!             DO ie=1, no_primal_edges
!               ictr =ictr+1
!               il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!               ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!               !thick_edge=patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
!               !thick_frac=thick_edge/thick_cell
!               out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
!                 & +vn_e(il_e,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,lev,blockNo,ictr))!*thick_frac
!
!             END DO
!
!           ENDIF
!         END DO edge_idx_loop
!       END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!
!     ELSEIF(PRESENT(scalar))THEN
!
!       DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!         CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
!
!         edge_idx_loop2: DO je =  start_edge_index, end_edge_index
!           IF (patch_3d%lsm_e(je,lev,blockNo) <= sea_boundary) THEN
!             ictr = 0
!             out_vn_e(je,blockNo) = 0.0_wp
!             il_c        = patch_2d%edges%cell_idx(je,blockNo,1)
!             ib_c        = patch_2d%edges%cell_blk(je,blockNo,1)
!             scalar_cell = scalar(il_c,ib_c)
!
!             DO ie=1, no_primal_edges
!               ictr =ictr+1
!               il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!               ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!
!               out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
!                 & +vn_e(il_e,ib_e)*scalar_cell &
!                 & *(operators_coefficients%edge2edge_viacell_coeff(je,lev,blockNo,ictr))!*thick_frac
!
!             END DO
!             ictr        = no_primal_edges
!             il_c        = patch_2d%edges%cell_idx(je,blockNo,2)
!             ib_c        = patch_2d%edges%cell_blk(je,blockNo,2)
!             scalar_cell = scalar(il_c,ib_c)
!             !thick_cell  = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)
!
!             DO ie=1, no_primal_edges
!               ictr =ictr+1
!               il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!               ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!
!               !thick_edge=patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
!               !thick_frac=thick_edge/thick_cell
!               out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
!                 & +vn_e(il_e,ib_e)*scalar_cell&
!                 & *(operators_coefficients%edge2edge_viacell_coeff(je,lev,blockNo,ictr))!*thick_frac
!
!             END DO
!           ENDIF
!         END DO edge_idx_loop2
!       END DO
!     ENDIF
!   END SUBROUTINE map_edges2edges_viacell_3d_1lev
!   !-----------------------------------------------------------------------------
!   !-----------------------------------------------------------------------------
!   !>
!   !!
!   SUBROUTINE map_edges2cell_with_height_3d( patch_3d, vn_e, operators_coefficients, p_vn_c, h_e,&
!     & opt_startLevel, opt_endLevel, subset_range)
!
!     TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
!     REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
!     ! 3D case: h_e is surface endLevelation at edges
!     TYPE(t_cartesian_coordinates),INTENT(inout):: p_vn_c(:,:,:)  ! outputput (nproma,n_zlev,alloc_cell_blocks)
!     REAL(wp), INTENT(in)                       :: h_e(:,:)       ! SW-case: h_e is thickness at edges
!     TYPE(t_operator_coeff)                     :: operators_coefficients
!     INTEGER, INTENT(in), OPTIONAL :: opt_startLevel       ! optional vertical start level
!     INTEGER, INTENT(in), OPTIONAL :: opt_endLevel       ! optional vertical end level
!     TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
!
!     !Local variables
!     !INTEGER, PARAMETER :: no_primal_edges = 3
!     INTEGER :: startLevel, endLevel
!     INTEGER :: start_cell_index, end_cell_index
!     INTEGER :: il_e, ib_e
!     INTEGER :: cell_index, blockNo, level, ie!,je
!     REAL(wp) :: z_weight
!     REAL(wp) :: z_thick_e
!
!     TYPE(t_subset_range), POINTER :: all_cells
!     !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
!     !  & routine = ('mo_scalar_product:primal_map_e2c')
!     TYPE(t_patch), POINTER :: patch_2d
!     !-----------------------------------------------------------------------
!     patch_2d   => patch_3d%p_patch_2d(1)
!     !-----------------------------------------------------------------------
!     IF ( PRESENT(subset_range) ) THEN
!       all_cells => subset_range
!     ELSE
!       all_cells => patch_2d%cells%ALL
!     ENDIF
!     !-----------------------------------------------------------------------
!     IF ( PRESENT(opt_startLevel) ) THEN
!       startLevel = opt_startLevel
!     ELSE
!       startLevel = 1
!     END IF
!     IF ( PRESENT(opt_endLevel) ) THEN
!       endLevel = opt_endLevel
!     ELSE
!       endLevel = n_zlev
!     END IF
!
!     IF ( iswm_oce == 1 ) THEN
!
!       !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
!       DO blockNo = all_cells%start_block, all_cells%end_block
!         CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
!         level_loop_swm: DO level = startLevel, endLevel
!           cell_idx_loop_swm: DO cell_index =  start_cell_index, end_cell_index
!             !calculate velocity reconstruction at cell center
!             z_weight           = 0.0_wp
!             p_vn_c(cell_index,level,blockNo)%x = 0.0_wp
!             DO ie=1, no_primal_edges
!
!               il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
!               ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)
!
!               z_thick_e =patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(il_e,level,ib_e)&
!                 & + h_e(il_e,ib_e)
!               z_weight = z_weight + operators_coefficients%variable_vol_norm(cell_index,level,blockNo,ie) * z_thick_e
!
!               p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x&
!                 & + operators_coefficients%edge2cell_coeff_cc_dyn(cell_index,level,blockNo,ie)%x&
!                 & * vn_e(il_e,level,ib_e)* z_thick_e
!             END DO
!             IF( z_weight/=0.0_wp)THEN
!               p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x / z_weight
!             ELSE
!               p_vn_c(cell_index,level,blockNo)%x=0.0_wp
!             ENDIF
!           END DO cell_idx_loop_swm
!         END DO level_loop_swm
!       END DO ! blockNo = all_cells%start_block, all_cells%end_block
!
!     ELSEIF( iswm_oce /= 1 ) THEN
!
!       !Step 1: Calculation of Pv in cartesian coordinates
!       DO blockNo = all_cells%start_block, all_cells%end_block
!         CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
!
!         !We are dealing with the surface layer first
!         cell_idx_loop_top: DO cell_index =  start_cell_index, end_cell_index
!           z_weight             = 0.0_wp
!           p_vn_c(cell_index,startLevel,blockNo)%x = 0.0_wp
!           DO ie=1, no_primal_edges
!
!             il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
!             ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)
!
!             z_thick_e = patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(il_e,startLevel,ib_e)&
!               & + h_e(il_e,ib_e)
!             z_weight = z_weight + operators_coefficients%variable_vol_norm(cell_index,startLevel,blockNo,ie) * z_thick_e
!
!             p_vn_c(cell_index,startLevel,blockNo)%x = p_vn_c(cell_index,startLevel,blockNo)%x&
!               & + operators_coefficients%edge2cell_coeff_cc_dyn(cell_index,1,blockNo,ie)%x&
!               & * vn_e(il_e,startLevel,ib_e) * z_thick_e
!
!           END DO
!
!           IF(z_weight/=0.0_wp)THEN
!             p_vn_c(cell_index,startLevel,blockNo)%x = p_vn_c(cell_index,startLevel,blockNo)%x / z_weight
!           ELSE
!             p_vn_c(cell_index,startLevel,blockNo)%x=0.0_wp
!           ENDIF
!         END DO cell_idx_loop_top
!
!         !Now we calculate at the levels below the surface
!         level_loop: DO level = startLevel+1, endLevel
!           cell_idx_loop: DO cell_index =  start_cell_index, end_cell_index
!             p_vn_c(cell_index,level,blockNo)%x = 0.0_wp
!             !z_weight = 0.0_wp
!             DO ie=1, no_primal_edges
!
!               il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
!               ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)
!               p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x&
!                 & + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,ie)%x&
!                 & * vn_e(il_e,level,ib_e)
!             END DO
!           END DO cell_idx_loop
!         END DO level_loop
!       END DO ! blockNo = all_cells%start_block, all_cells%end_block
!     ENDIF
!
!     ! LL no sync required
!
!   END SUBROUTINE map_edges2cell_with_height_3d
!   !----------------------------------------------------------------

  ! !   !-----------------------------------------------------------------------------
  ! !   !>
  ! !   !! Discrete mapping of cell-based vectors to edges on the primal grid.
  ! !   !!
  ! !   !!
  ! !   !!  mpi parallelized LL, result not synced
  ! !   SUBROUTINE map_cell2edges_2d( patch_3D, p_vn_c, ptp_vn, operators_coefficients)
  ! !
  ! !     TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3D
  ! !     TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:)    ! input vector (nproma,n_zlev,alloc_cell_blocks)
  ! !     REAL(wp), INTENT(inout)                    :: ptp_vn(:,:)    ! output vector (nproma,n_zlev,nblks_e)
  ! !     TYPE(t_operator_coeff)                     :: operators_coefficients
  ! !
  ! !     !Local variables
  ! !     INTEGER :: start_edge_index, end_edge_index
  ! !     INTEGER :: cell_1_index, cell_1_block, cell_2_index, cell_2_block
  ! !     INTEGER :: je, blockNo
  ! !     TYPE(t_subset_range), POINTER :: edges_in_domain
  ! !     TYPE(t_patch), POINTER        :: patch_2D
  ! !     !-----------------------------------------------------------------------
  ! !     patch_2D   => patch_3D%p_patch_2D(1)
  ! !     !-----------------------------------------------------------------------
  ! !     !CALL message (TRIM(routine), 'start')
  ! !     edges_in_domain   => patch_2D%edges%in_domain
  ! !
  ! !     ! calculation of transposed P^TPv from Pv (incart coord)
  ! !     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
  ! !       CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
  ! !
  ! !       edge_idx_loop: DO je =  start_edge_index, end_edge_index
  ! !
  ! !         IF(patch_3D%lsm_e(je,1,blockNo) <= sea_boundary)THEN
  ! !           !Get indices of two adjacent triangles
  ! !           cell_1_index = patch_2D%edges%cell_idx(je,blockNo,1)
  ! !           cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
  ! !           cell_2_index = patch_2D%edges%cell_idx(je,blockNo,2)
  ! !           cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)
  ! !
  ! !           ptp_vn(je,blockNo) = &
  ! !             &  DOT_PRODUCT(p_vn_c(cell_1_index,cell_1_block)%x, operators_coefficients%edge2cell_coeff_cc_t(je,1,blockNo,1)%x)&
  ! !             & +DOT_PRODUCT(p_vn_c(cell_2_index,cell_2_block)%x, operators_coefficients%edge2cell_coeff_cc_t(je,1,blockNo,2)%x)
  ! !         ELSE
  ! !           ptp_vn(je,blockNo) = 0.0_wp
  ! !         ENDIF
  ! !
  ! !       END DO edge_idx_loop
  ! !     END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
  ! !
  ! ! !     CALL sync_patch_array(SYNC_E, patch_2D, ptp_vn)
  ! !
  ! !   END SUBROUTINE map_cell2edges_2d
  !-----------------------------------------------------------------------------



END MODULE mo_scalar_product
