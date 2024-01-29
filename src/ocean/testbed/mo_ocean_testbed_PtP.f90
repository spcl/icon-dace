! Contains the implementation of the mathematical operators for the ocean.
!
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
#include "omp_definitions.inc"
#include "iconfor_dsl_definitions.inc"

MODULE mo_ocean_testbed_PtP
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_run_config,         ONLY: test_mode
  USE mo_math_constants
  USE mo_physical_constants
  USE mo_impl_constants,     ONLY: boundary, sea, sea_boundary !,sea,land, land_boundary, sea, max_char_length, &
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce
  USE mo_dynamics_config,    ONLY: nold, nnew
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  !USE mo_exception,          ONLY: finish, message
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_div, timer_grad
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_math_utilities,      ONLY: vector_product
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_grid_config,         ONLY: n_dom
  USE mo_ocean_math_operators
  USE mo_scalar_product,      ONLY:  map_edges2edges_viacell_3d_const_z, map_edges2cell_3d

  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_physics_types, ONLY: t_ho_params
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff, no_primal_edges, no_dual_edges
  USE mo_ocean_math_operators,ONLY: div_oce_3d
  USE mo_statistics

  IMPLICIT NONE

  PRIVATE
  REAL(wp), POINTER                    :: vn(:,:,:)
  REAL(wp), POINTER                    :: PtPvn(:,:,:)
  REAL(wp), POINTER                    :: u(:,:,:), v(:,:,:)

  TYPE(t_patch_3D ),POINTER           :: patch_3d
  TYPE(t_patch), POINTER              :: patch_2d
  TYPE(t_hydro_ocean_state), POINTER  :: ocean_state
  TYPE(t_operator_coeff),POINTER      :: operators_coefficients

  CHARACTER(len=12)           :: str_module    = 'test_PtP'      ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

    onCells_Type(t_cartesian_coordinates) :: &
      & Pvn              ! reconstructed velocity at cell center in cartesian coordinates

  PUBLIC :: test_PtP

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Note: not tested with MPI !
  SUBROUTINE test_PtP(in_patch_3D, in_ocean_state, in_operators_coefficients, test_mode)
    TYPE(t_patch_3d ),TARGET,INTENT(in) :: in_patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: in_ocean_state
    TYPE(t_operator_coeff),TARGET, INTENT(in)   :: in_operators_coefficients
    INTEGER :: test_mode

    ! initialize the div test memory
    patch_3d   => in_patch_3D
    patch_2d   => patch_3d%p_patch_2d(1)
    operators_coefficients => in_operators_coefficients

    ocean_state  => in_ocean_state
    vn           => ocean_state%p_diag%mass_flx_e
    PtPvn        => ocean_state%p_diag%ptp_vn
    Pvn          => ocean_state%p_diag%p_vn 
    u            => ocean_state%p_diag%u
    v            => ocean_state%p_diag%v
     
    SELECT CASE (test_mode) ! 100 - 999

      CASE (114)
        CALL test_PtP_onTorus_divTest()
       
      CASE (115)
        CALL test_PtP_onIcos_divTest()
       
      CASE DEFAULT
        CALL finish("test_PtP", "Unknown test_mode")

    END SELECT

  END SUBROUTINE test_PtP
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_PtP_onTorus_divTest()

    write(0,*) "-------------- test PtP torus div test 1 -----------------------"
    CALL fill_vn_onTorus_divTest_1() 
    CALL diagnose_PtP_results()

    write(0,*) "-------------- test PtP torus div test 2 -----------------------"
    CALL fill_vn_onTorus_divTest_2() 
    CALL diagnose_PtP_results()


    write(0,*) "-------------- test PtP torus flow test 1 -----------------------"
    CALL fill_vn_onTorus_flowTest_1()
    CALL diagnose_PtP_results()

    write(0,*) "-------------- test PtP torus flow test 2 -----------------------"
    CALL fill_vn_onTorus_flowTest_2()
    CALL diagnose_PtP_results()

  END SUBROUTINE test_PtP_onTorus_divTest
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_PtP_onIcos_divTest()

    write(0,*) "-------------- test PtP Icos div test 1 -----------------------"
    CALL fill_vn_onIcos_divTest_1() 
    CALL diagnose_PtP_results()

 
  END SUBROUTINE test_PtP_onIcos_divTest
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE diagnose_PtP_results()

    REAL(wp), POINTER                    :: Pvn_x(:,:,:)
    !-----------------------------------------------------------------------
    CALL map_edges2edges_viacell_3d_const_z( patch_3d, vn(:,:,:), operators_coefficients, PtPvn(:,:,:) )
    CALL map_edges2cell_3d(patch_3d, vn, operators_coefficients, Pvn) !, subset_range=cells_in_domain)

    write(0,*) "====================================================================="
!     CALL print_PtP_coefficients
    
    write(0,*) "====================================================================="
!     CALL dbg_print('thick_c',patch_3d%p_patch_1d(1)%prism_thick_c,str_module, 1, in_subset=patch_2d%cells%owned)
    CALL dbg_print('thick_e',patch_3d%p_patch_1d(1)%prism_thick_e,str_module, 1, in_subset=patch_2d%edges%owned)
    CALL dbg_print('vn',vn,str_module,1, in_subset=patch_2d%edges%owned)
    CALL dbg_print('PtPvn',PtPvn,str_module,1, in_subset=patch_2d%edges%owned)
    write(0,*) "innder product:", inner_product(vn, vn), inner_product(vn, PtPvn)

!     u = Pvn(1,1,1)%x(1)
!     v = Pvn(1,1,1)%x(2)
!     CALL dbg_print('Pvn 1', Pvn(:,:,:)%x(1),str_module,1, in_subset=patch_2d%cells%owned)
!     CALL dbg_print('Pvn 2', Pvn(:,:,:)%x(2),str_module,1, in_subset=patch_2d%cells%owned)

    write(0,*) "====================================================================="
    write(0,*) 
!     write(0,*) "============== vn ===================================="
!     CALL print_vn(vn)
!     write(0,*) "============== PtPvn ===================================="
!     CALL print_vn(PtPvn)
  
  END SUBROUTINE diagnose_PtP_results
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  REAL(wp) FUNCTION inner_product(vn, PtPvn)
    REAL(wp), INTENT(in) :: vn(:,:,:), PtPvn(:,:,:)

    inner_product = subset_sum(vn(:,1,:), PtPvn(:,1,:), patch_2d%edges%owned)

  END FUNCTION inner_product
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE print_vn(vn)
  
    REAL(wp), INTENT(in) :: vn(:,:,:)
    !-----------------------------------------------------------------------
 
    write(0,*) vn(:,1,:)

  END SUBROUTINE print_vn
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE fill_vn_onIcos_divTest_1()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: v1, v2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, v1, v2) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        v1 = patch_2d%edges%primal_normal(j,block)%v1
        v2 = patch_2d%edges%primal_normal(j,block)%v2

        IF (abs(v1) < 1.0e-16_wp) THEN 
          ! horizontal edge, or right_diagonal edge
          vn(j,1,block) = SIGN(1.0_wp, v2)
        ELSE
          vn(j,1,block) = SIGN(1.0_wp, v1)
        ENDIF

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_onIcos_divTest_1
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE fill_vn_onTorus_flowTest_1()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: v1, v2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, v1, v2) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        v1 = patch_2d%edges%primal_normal(j,block)%v1
        v2 = patch_2d%edges%primal_normal(j,block)%v2

!         IF (v1 == 0.0_wp .or. v2 < 0.0_wp) THEN 
          ! horizontal edge, or right_diagonal edge
          vn(j,1,block) = 1.0_wp
!         ELSE
!           vn(j,1,block) = -1.0_wp
!         ENDIF

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_onTorus_flowTest_1
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE fill_vn_onTorus_flowTest_2()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: v1, v2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, v1, v2) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        v1 = patch_2d%edges%primal_normal(j,block)%v1
        v2 = patch_2d%edges%primal_normal(j,block)%v2

        IF (v1 == 0.0_wp) THEN 
          ! horizontal edge
          vn(j,1,block) = 0.0_wp
        ELSE
          vn(j,1,block) = 1.0_wp
        ENDIF

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_onTorus_flowTest_2
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE fill_vn_onTorus_divTest_1()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: v1, v2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, v1, v2) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        v1 = patch_2d%edges%primal_normal(j,block)%v1
        v2 = patch_2d%edges%primal_normal(j,block)%v2

        IF (v1 == 0.0_wp .or. v2 < 0.0_wp) THEN 
          ! horizontal edge, or right_diagonal edge
          vn(j,1,block) = 1.0_wp
        ELSE
          vn(j,1,block) = -1.0_wp
        ENDIF

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_onTorus_divTest_1
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE fill_vn_onTorus_divTest_2()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: v1, v2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, v1, v2) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        v1 = patch_2d%edges%primal_normal(j,block)%v1
        v2 = patch_2d%edges%primal_normal(j,block)%v2

        IF (v1 == 0.0_wp) THEN 
          ! horizontal edge, or right_diagonal edge
          vn(j,1,block) = 1.1_wp
        ELSEIF (v2 < 0.0_wp) THEN 
          ! horizontal edge, or right_diagonal edge
          vn(j,1,block) = 1.0_wp
        ELSE
          vn(j,1,block) = -1.0_wp
        ENDIF

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_onTorus_divTest_2
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE print_PtP_coefficients()

    TYPE(t_subset_range), POINTER :: all_edges
    INTEGER :: block, j, start_index,  end_index

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
    INTEGER :: edge_11_block, edge_12_block, edge_13_block
    INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
    INTEGER :: edge_21_block, edge_22_block, edge_23_block
    REAL(wp) :: abs_sum
   !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    !-----------------------------------------------------------------------

    write(0,*) "-------------- PtP coefficients -----------------------"


    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        
        IF (patch_3d%p_patch_1d(1)%dolic_e(j,block) <= 0) CYCLE 

        cell_1_index = patch_2d%edges%cell_idx(j,block,1)
        cell_1_block = patch_2d%edges%cell_blk(j,block,1)
        cell_2_index = patch_2d%edges%cell_idx(j,block,2)
        cell_2_block = patch_2d%edges%cell_blk(j,block,2)

        IF (patch_3d%lsm_c(cell_1_index,1,cell_1_block) /= sea) CYCLE
        IF (patch_3d%lsm_c(cell_2_index,1,cell_2_block) /= sea) CYCLE

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
  
        abs_sum = &
            &   ABS(operators_coefficients%edge2edge_viacell_coeff(j,1,block,1)) & 
            & + ABS(operators_coefficients%edge2edge_viacell_coeff(j,1,block,2)) & 
            & + ABS(operators_coefficients%edge2edge_viacell_coeff(j,1,block,3)) & 
            & + ABS(operators_coefficients%edge2edge_viacell_coeff(j,1,block,4)) & 
            & + ABS(operators_coefficients%edge2edge_viacell_coeff(j,1,block,5)) & 
            & + ABS(operators_coefficients%edge2edge_viacell_coeff(j,1,block,6))

          WRITE(0,*) "> ", &
            & operators_coefficients%edge2edge_viacell_coeff(j,1,block,1), & 
            & operators_coefficients%edge2edge_viacell_coeff(j,1,block,2), & 
            & operators_coefficients%edge2edge_viacell_coeff(j,1,block,3), & 
            & operators_coefficients%edge2edge_viacell_coeff(j,1,block,4), & 
            & operators_coefficients%edge2edge_viacell_coeff(j,1,block,5), & 
            & operators_coefficients%edge2edge_viacell_coeff(j,1,block,6), &
            & " = ", abs_sum

      ENDDO
    ENDDO

    write(0,*) "-------------------------------------------------------"

  END SUBROUTINE print_PtP_coefficients
  !-------------------------------------------------------------------------

END MODULE mo_ocean_testbed_PtP

