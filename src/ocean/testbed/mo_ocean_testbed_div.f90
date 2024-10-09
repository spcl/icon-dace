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

#include "omp_definitions.inc"
#include "iconfor_dsl_definitions.inc"

MODULE mo_ocean_testbed_div
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_run_config,         ONLY: test_mode
  USE mo_math_constants,     ONLY: pi
  USE mo_physical_constants, ONLY: earth_radius
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
  USE mo_scalar_product,      ONLY:  map_edges2edges_viacell_3d_const_z


  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_physics_types, ONLY: t_ho_params
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff, no_primal_edges, no_dual_edges
  USE mo_ocean_math_operators,ONLY: div_oce_3d, rot_vertex_ocean_3D
  USE mo_statistics
  
  USE mo_name_list_output,       ONLY: write_name_list_output
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode

  IMPLICIT NONE

  PRIVATE
  REAL(wp), POINTER                    :: vn(:,:,:), u(:,:,:), v(:,:,:)
  REAL(wp), POINTER                    :: div_model(:,:,:)
  REAL(wp), POINTER                    :: div_analytic(:,:,:)
  REAL(wp), POINTER                    :: div_diff(:,:,:)
  REAL(wp), POINTER                    :: PtPvn(:,:,:)
  REAL(wp), POINTER                    :: u_vert(:,:), v_vert(:,:), div_vert(:,:)
  REAL(wp), POINTER                    :: divPtP(:,:,:), divPtP_diff(:,:,:)
  onVertices :: vort                   ! vorticity at triangle vertices. Unit [1/s]
  onVertices_Type(t_cartesian_coordinates) :: vn_dual


  REAL(wp) :: minmaxmean(3), L2Diff, L2DivAn, LInfDiff, LInfDivAn

  TYPE(t_patch_3D ),POINTER           :: patch_3d
  TYPE(t_patch), POINTER              :: patch_2d
  TYPE(t_hydro_ocean_state), POINTER  :: ocean_state
  TYPE(t_operator_coeff),POINTER      :: operators_coefficients

  CHARACTER(len=12)           :: str_module    = 'test_div'      ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug


  PUBLIC :: test_div

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Note: not tested with MPI !
  SUBROUTINE test_div(in_patch_3D, in_ocean_state, in_operators_coefficients, test_mode)
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
    div_analytic => ocean_state%p_diag%div_mass_flx_c

    div_model    => ocean_state%p_diag%div_model   
    div_diff     => ocean_state%p_diag%div_diff
    PtPvn        => ocean_state%p_diag%ptp_vn
    divPtP       => ocean_state%p_diag%divPtP
    divPtP_diff  => ocean_state%p_diag%divPtP_diff
    vort         => ocean_state%p_diag%vort
    vn_dual      => ocean_state%p_diag%p_vn_dual
    
    ALLOCATE(u_vert(nproma,patch_2D%nblks_v), v_vert(nproma,patch_2D%nblks_v), div_vert(nproma,patch_2D%nblks_v))


    SELECT CASE (test_mode) ! 100 - 999

      CASE (103)
        CALL test_div_accuracy_onSphere_HeikesRandall()
       
      CASE (104)
        CALL test_div_accuracy_onPlane_HeikesRandall()
            
      CASE (105)
        CALL test_div_accuracy_onPlane_Hui()

      CASE (106)
        CALL test_div_accuracy_onPlane_Basic()

      CASE (107)
        CALL test_div_accuracy_onPlane_Basic_Accurate()

      CASE (109)
        CALL test_vn_accuracy_onPlaneQuads_Hui()

      CASE (110)
        CALL test_div_onTorus()

      CASE DEFAULT
        CALL finish("test_div", "Unknown test_mode")

    END SELECT

    DEALLOCATE(u_vert, v_vert, div_vert)

  END SUBROUTINE test_div
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_div_accuracy_onPlane_Basic()

    write(0,*) "-------------- test div accuracy plane Basic -----------------------"

!     CALL fill_vn_divAnalytic_plane_basic2()
!     CALL fill_vn_divAnalytic_plane_basic2_onVertices()
    CALL fill_vn_divAnalytic_plane_basic2_3order()
 
    CALL diagnose_div_accuracy()

  END SUBROUTINE test_div_accuracy_onPlane_Basic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_div_accuracy_onPlane_Basic_Accurate()

    write(0,*) "-------------- test div accuracy plane Basic Accurate -----------------------"

    CALL fill_vn_divAnalytic_plane_basic2_accurate()
    vn => PtPvn
    CALL fill_vn_divAnalytic_plane_basic2()
    return

    CALL diagnose_div_accuracy()

  END SUBROUTINE test_div_accuracy_onPlane_Basic_Accurate
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_div_accuracy_onSphere_HeikesRandall()

    write(0,*) "-------------- test div accuracy onSphere HeikesRandall -----------------------"

    CALL fill_vn_divAnalytic_onSphere_HeikesRandall
 
    CALL diagnose_div_accuracy()

  END SUBROUTINE test_div_accuracy_onSphere_HeikesRandall
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_div_accuracy_onPlane_HeikesRandall()

    write(0,*) "-------------- test div accuracy plane HeikesRandall -----------------------"

    CALL fill_vn_divAnalytic_plane_HeikesRandall()
 
    CALL diagnose_div_accuracy()

  END SUBROUTINE test_div_accuracy_onPlane_HeikesRandall
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_div_accuracy_onPlane_Hui()


    write(0,*) "-------------- test div accuracy plane Hui -----------------------"

    CALL fill_vn_divAnalytic_plane_Hui()
    CALL fill_vn_divAnalytic_plane_Hui_onVertices
    CALL interpolate_vn_fromVerticesEdges

!     CALL fill_vn_divAnalytic_fromVertices

!     CALL fill_vnAveraged_planeQuads_Hui
!     CALL fill_divAnalyticAverage_planeQuads_Hui
 
    CALL diagnose_div_accuracy()

  END SUBROUTINE test_div_accuracy_onPlane_Hui
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  ! 
  SUBROUTINE test_div_onTorus()
  
    CALL fill_div_onTorus()
 
    CALL diagnose_div_accuracy()
     
    IF (output_mode%l_nml) CALL write_name_list_output(0)

  END SUBROUTINE test_div_onTorus
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE diagnose_div_accuracy()

    TYPE(t_subset_range), POINTER :: owned_cells
    !-----------------------------------------------------------------------
    owned_cells => patch_2d%cells%owned
    !-----------------------------------------------------------------------
    CALL map_edges2edges_viacell_3d_const_z( patch_3d, vn(:,:,:), operators_coefficients, PtPvn(:,:,:) )

    CALL div_oce_3d( vn, patch_3D, operators_coefficients%div_coeff, div_model)
    CALL div_oce_3d( PtPvn, patch_3D, operators_coefficients%div_coeff, divPtP)
    CALL rot_vertex_ocean_3D( patch_3D, vn, vn_dual, operators_coefficients, vort)
    
!     CALL dbg_print('prism_thick_e',patch_3d%p_patch_1d(1)%prism_thick_e, &
!         & str_module,1, in_subset=patch_2d%edges%owned)
!     CALL dbg_print('h',ocean_state%p_prog(nold(1))%h,str_module,1, in_subset=patch_2d%cells%owned)
! 
!     CALL dbg_print('thick_c',ocean_state%p_diag%thick_c,str_module,1, in_subset=patch_2d%cells%owned)

    div_diff    = div_model - div_analytic
    divPtP_diff = divPtP    - div_analytic

    write(0,*) "=================================="
    CALL dbg_print('rot_coeff',operators_coefficients%rot_coeff(:,:,:,1), str_module,1, in_subset=patch_2d%verts%owned)
    CALL dbg_print('rot_coeff',operators_coefficients%rot_coeff(:,:,:,2), str_module,1, in_subset=patch_2d%verts%owned)
    CALL dbg_print('rot_coeff',operators_coefficients%rot_coeff(:,:,:,3), str_module,1, in_subset=patch_2d%verts%owned)
    CALL dbg_print('rot_coeff',operators_coefficients%rot_coeff(:,:,:,4), str_module,1, in_subset=patch_2d%verts%owned)
    CALL dbg_print('rot_coeff',operators_coefficients%rot_coeff(:,:,:,5), str_module,1, in_subset=patch_2d%verts%owned)
    CALL dbg_print('rot_coeff',operators_coefficients%rot_coeff(:,:,:,6), str_module,1, in_subset=patch_2d%verts%owned)
     
    write(0,*) "=================================="
    CALL dbg_print('div_coeff',operators_coefficients%div_coeff(:,:,:,1), str_module,1, in_subset=patch_2d%cells%owned)
    CALL dbg_print('div_coeff',operators_coefficients%div_coeff(:,:,:,2), str_module,1, in_subset=patch_2d%cells%owned)
    CALL dbg_print('div_coeff',operators_coefficients%div_coeff(:,:,:,3), str_module,1, in_subset=patch_2d%cells%owned)

    write(0,*) "=================================="
    CALL dbg_print('vn',vn,str_module,1, in_subset=patch_2d%edges%owned)
    CALL dbg_print('PtPvn',PtPvn,str_module,1, in_subset=patch_2d%edges%owned)
    CALL dbg_print('div_analytic',div_analytic,str_module,1, in_subset=patch_2d%cells%owned)
    CALL dbg_print('div_model',div_model,str_module,1, in_subset=patch_2d%cells%owned)
    CALL dbg_print('div_diff',div_diff,str_module,1, in_subset=patch_2d%cells%owned)
    CALL dbg_print('divPtP',divPtP,str_module,1, in_subset=patch_2d%cells%owned)
    CALL dbg_print('divPtP_diff',divPtP_diff,str_module,1, in_subset=patch_2d%cells%owned)
   
    L2Diff  = L2Norm(div_diff(:,1,:), owned_cells)
    L2DivAn = L2Norm(div_analytic(:,1,:), owned_cells)
    L2Diff = L2Diff / L2DivAn

    LInfDiff  = LInfNorm(div_diff(:,1,:), owned_cells)
    LInfDivAn = LInfNorm(div_analytic(:,1,:), owned_cells)
    LInfDiff  = LInfDiff / LInfDivAn

    write(0,*) "=================================="
    write(0,*) " Resolution: ", patch_2d%geometry_info%mean_characteristic_length
    write(0,*) "=================================="
    write(0,*)  "    Without PtP"
    write(0,*) "L2 error:",   L2Diff
    write(0,*) "LInf error:", LInfDiff
 
    L2Diff  = L2Norm(divPtP_diff(:,1,:), owned_cells)
    L2Diff = L2Diff / L2DivAn

    LInfDiff  = LInfNorm(divPtP_diff(:,1,:), owned_cells)
    LInfDiff  = LInfDiff / LInfDivAn

    write(0,*) "=================================="
    write(0,*)  "      With PtP"
    write(0,*) "L2 error:",   L2Diff
    write(0,*) "LInf error:", LInfDiff
    write(0,*) "=================================="

  END SUBROUTINE diagnose_div_accuracy
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  ! u = - cos^3(lat) * sin^2(lon)
  ! v = -4 * cos^3(lat) * sin(lat) * sin(lon) * cos(lon)
  ! div = 1/(r*cos(lat) * d(u)/d(lon) + 1/(r*cos(lat) * d(v*cos(lat))/d(lat) = 1/r *
  !  [(-cos^2(lat) * sin(2*lon)) 
  !   -(2*sin(2*lon) * cos^2(lat) * (cos^2(lat) - 4*sin^2(lat) )) =
  !    
  !  - sin(2*lon)/ r *
  ! (cos^2(lat) + (2 * cos^2(lat) *  (cos^2(lat) - 4*sin^2(lat) ))) = 
  ! - sin(2*lon) * cos^2(lat) / r * 
  ! (1 + 2*(cos^2(lat) - 4*sin^2(lat) ) )
  ! 
  SUBROUTINE fill_vn_divAnalytic_onSphere_HeikesRandall()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, u, v
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%edges%center(j,block)%lon
        lat = patch_2d%edges%center(j,block)%lat
        
        u = -cos(lat)**3 * sin(lon)**2
        v = -2.0_wp * cos(lat)**3 * sin(lat) * sin(2.0_wp*lon)
        vn(j,1,block) =  u * patch_2d%edges%primal_normal(j,block)%v1 &
                       & + v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
        
        div_analytic (j,1,block) = &
          ! - sin(2*lon) * cos^2(lat) / r * 
          ! (1 + 2*(cos^2(lat) - 4*sin^2(lat) ) )
          & - (sin(2.0_wp*lon) * cos(lat)**2 / earth_radius) *      &
          & (1.0_wp + 2.0_wp *(cos(lat)**2 - 4.0_wp*sin(lat)**2 ) ) 

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_onSphere_HeikesRandall
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! u = - cos^3(lat) * sin^2(lon)
  ! v = -4 * cos^3(lat) * sin(lat) * sin(lon) * cos(lon)
  ! div = [-cos^3(lat) * sin(2*lon)] +
  !       [2*sin(2*lon)*cos^2(lat) * (3*sin^2(lat)-cos^2(lat))]
  ! 
  SUBROUTINE fill_vn_divAnalytic_plane_HeikesRandall()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, u, v
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%edges%center(j,block)%lon
        lat = patch_2d%edges%center(j,block)%lat
        
        u = -cos(lat)**3 * sin(lon)**2
        v = -2.0_wp * cos(lat)**3 * sin(lat) * sin(2.0_wp*lon)
        vn(j,1,block) =  u * patch_2d%edges%primal_normal(j,block)%v1 &
                       & + v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
        
        div_analytic (j,1,block) = &
          !  [-cos^3(lat) * sin(2*lon)] +
          !  [2*sin(2*lon)*cos^2(lat) * (3*sin^2(lat)-cos^2(lat))]
          &( sin(2.0_wp * lon) * cos(lat)**2 ) *      &
          & (-cos(lat) + (6.0_wp * sin(lat)**2 - 2.0_wp * cos(lat)**2 ) ) / earth_radius

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_plane_HeikesRandall
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! u = 1/4 * sqrt(105/(2*pi)) cos(2*lon) cos^2(lat) sin(lat)
  ! v = -1/2 * sqrt(15/(2*pi) * cos(lon) * cos(lat) * sin(lat)
  ! div = -1/(2*sqr(2*pi)) (sqrt(105) * sin(2*lon) * cos^2(lat) * sin(lat) + sqrt(15) * cos(lon) cos(2*lon))
  ! 
  SUBROUTINE fill_vn_divAnalytic_plane_Hui()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, u, v
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%edges%center(j,block)%lon
        lat = patch_2d%edges%center(j,block)%lat
        
        u = 0.25_wp * sqrt_105 / sqrt_pi_x_2 * cos(2.0_wp*lon) * cos(lat)**2 * sin(lat)
        v = - 0.5_wp * sqrt_15 / sqrt_pi_x_2 * cos(lon) * cos(lat) * sin(lat) 
        vn(j,1,block) =  u * patch_2d%edges%primal_normal(j,block)%v1 &
                       & + v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
        
        div_analytic (j,1,block) = &
           & - 1.0_wp/(2._wp*sqrt_pi_x_2 * earth_radius) * (sqrt_105 * sin(2.0_wp * lon) * cos(lat)**2 * sin(lat) + &
             sqrt_15 * cos(lon) * cos(2.0_wp*lat))     
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_plane_Hui
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! 
  SUBROUTINE fill_div_onTorus()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, x, y, u, v
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, x, y, u, v, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        x = patch_2d%edges%cartesian_center(j,block)%x(1)
        y = patch_2d%edges%cartesian_center(j,block)%x(2)
 
        lon = patch_2d%edges%center(j,block)%lon
        lat = patch_2d%edges%center(j,block)%lat
        
!   1st case
        !         u = sin(x*pi)
!         v = cos(y*pi) ! note we need at least to be continues at +-pi/2
        
!   2nd case
!         u = sin(x*pi)*cos(y*pi)
!         v = cos(y*pi)**2
         
         
         
!   3 case (from Hui but with x,y)
        u = 0.25_wp * sqrt_105 / sqrt_pi_x_2 * cos(2.0_wp*x*pi) * cos(y*pi)**2 * sin(y*pi)
        v = - 0.5_wp * sqrt_15 / sqrt_pi_x_2 * cos(x*pi) * cos(y*pi) * sin(y*pi) 

!   4 case (from Hui with lon/lat)
!         u = 0.25_wp * sqrt_105 / sqrt_pi_x_2 * cos(2.0_wp*lon) * cos(lat)**2 * sin(lat)
!         v = - 0.5_wp * sqrt_15 / sqrt_pi_x_2 * cos(lon) * cos(lat) * sin(lat) 

        vn(j,1,block) =  u * patch_2d%edges%primal_normal(j,block)%v1 &
                     & + v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, x, y, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        x = patch_2d%cells%cartesian_center(j,block)%x(1)
        y = patch_2d%cells%cartesian_center(j,block)%x(2)
 
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
        
!   1st case
!         div_analytic (j,1,block) = (cos(x*pi) - sin(y*pi)) * pi 

!   2nd case
!         div_analytic (j,1,block) = (cos(x*pi) * cos(y*pi) - 2.0_wp * cos(y*pi) * sin(y*pi)) * pi 


!   3 case (from Hui but with x,y)
         div_analytic (j,1,block) = &
           & - 1.0_wp/(2._wp*sqrt_pi_x_2) * (sqrt_105 * sin(2.0_wp * x*pi) * cos(y*pi)**2 * sin(y*pi) + &
             sqrt_15 * cos(x*pi) * cos(2.0_wp*y*pi)) * pi  
             
 !   4 case (from Hui with lon/lat)
!         div_analytic (j,1,block) = &
!            & -1.0_wp/(2._wp*sqrt_pi_x_2) * (sqrt_105 * sin(2.0_wp * lon) * cos(lat)**2 * sin(lat) + &
!              sqrt_15 * cos(lon) * cos(2.0_wp*lat)) * (2.0_wp * pi / patch_2D%geometry_info%domain_length)  
            
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO
     

  END SUBROUTINE fill_div_onTorus
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  ! Su dy = 1/12 * sqrt(105/(2*pi)) cos(2*lon_0) (cos^3(lat_0) - cos^3(lat_1))
  ! Sv dx = 1/4 * sqrt(15/(2*pi) * sin(2*lat_0) * (sin(lon_0) - sin(lon_1))
  ! 
  SUBROUTINE fill_vnAveraged_planeQuads_Hui()
  
    TYPE(t_subset_range), POINTER :: all_edges
    REAL(wp) :: lon_0, lon_1, lat_0, lat_1, swap
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index
    INTEGER :: vertex_1_idx,vertex_1_blk,vertex_2_idx,vertex_2_blk
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon_0, lon_1, lat_0, lat_1, swap, &
!ICON_OMP vertex_1_idx,vertex_1_blk,vertex_2_idx,vertex_2_blk) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        vertex_1_idx = patch_2d%edges%vertex_idx(j,block,1)
        vertex_1_blk = patch_2d%edges%vertex_blk(j,block,1)
        vertex_2_idx = patch_2d%edges%vertex_idx(j,block,2)
        vertex_2_blk = patch_2d%edges%vertex_blk(j,block,2)
      
        lat_0 = MIN(patch_2d%verts%vertex(vertex_1_idx,vertex_1_blk)%lat, patch_2d%verts%vertex(vertex_2_idx,vertex_2_blk)%lat)
        lat_1 = MAX(patch_2d%verts%vertex(vertex_1_idx,vertex_1_blk)%lat, patch_2d%verts%vertex(vertex_2_idx,vertex_2_blk)%lat)
        lon_0 = MIN(patch_2d%verts%vertex(vertex_1_idx,vertex_1_blk)%lon, patch_2d%verts%vertex(vertex_2_idx,vertex_2_blk)%lon)
        lon_1 = MAX(patch_2d%verts%vertex(vertex_1_idx,vertex_1_blk)%lon, patch_2d%verts%vertex(vertex_2_idx,vertex_2_blk)%lon)       

        if (lon_1 - lon_0 > pi) then
          swap=lon_1
          lon_1=lon_0
          lon_0=swap
        endif
!         u = 0.25_wp * sqrt_105 / sqrt_pi_x_2 * cos(2.0_wp*lon) * cos(lat)**2 * sin(lat)
!         v = - 0.5_wp * sqrt_15 / sqrt_pi_x_2 * cos(lon) * cos(lat) * sin(lat) 

        IF (lat_0 == lat_1) THEN
          ! horizontal edges
          ! Sv dx = 1/4 * sqrt(15/(2*pi) * sin(2*lat_0) * (sin(lon_0) - sin(lon_1))
          vn(j,1,block) =  0.25_wp * sqrt_15/sqrt_pi_x_2 * sin(2.0_wp*lat_0) * (sin(lon_0) - sin(lon_1)) &
            & * patch_2d%edges%primal_normal(j,block)%v2
        ELSE
          ! vertical edge
          ! Su dy = 1/12 * sqrt(105/(2*pi)) cos(2*lon_0) (cos^3(lat_0) - cos^3(lat_1))
          vn(j,1,block) =  1.0_wp/12.0_wp * sqrt_105/sqrt_pi_x_2 * cos(2.0_wp*lon_0) * (cos(lat_0)**3 - cos(lat_1)**3) &
            & * patch_2d%edges%primal_normal(j,block)%v1
        ENDIF
        vn(j,1,block) = vn(j,1,block)  * earth_radius / patch_2d%edges%primal_edge_length(j,block)

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vnAveraged_planeQuads_Hui
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! u = 1/4 * sqrt(105/(2*pi)) cos(2*lon) cos^2(lat) sin(lat)
  ! v = -1/2 * sqrt(15/(2*pi) * cos(lon) * cos(lat) * sin(lat)
  !   SS(du/dx)dxdy = -1/12 sqrt(105/2*pi) * (cos(2*lon_1) - cos(2*lon_0)) * (c0s^3(lat_1) - cos^3(lat_0))
  !   SS(dv/dy)dydx = -1/4 * sqrt(15/2*pi) * (sin(2*lat_1) - sin(2*lat_0)) * (sin(lon_1) - sin(lon_0)) 
  ! 
  SUBROUTINE fill_divAnalyticAverage_planeQuads_Hui()
  
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lat_0,lat_1,lon_0,lon_1,swap
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index, neigbor, vertex_idx, vertex_blk
    !-----------------------------------------------------------------------
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)


!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, neigbor, vertex_idx, vertex_blk, &
!ICON_OMP lat_0,lat_1,lon_0,lon_1,swap) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lat_0 =  999.0_wp
        lat_1 = -999.0_wp
        lon_0 =  999.0_wp
        lon_1 = -999.0_wp
        DO neigbor=1, 4 ! only works for quads    
          vertex_idx = patch_2d%cells%vertex_idx(j,block,neigbor)
          vertex_blk = patch_2d%cells%vertex_blk(j,block,neigbor)
      
          lat_0 = MIN(lat_0,patch_2d%verts%vertex(vertex_idx,vertex_blk)%lat)
          lat_1 = MAX(lat_1,patch_2d%verts%vertex(vertex_idx,vertex_blk)%lat)
          lon_0 = MIN(lon_0,patch_2d%verts%vertex(vertex_idx,vertex_blk)%lon+pi)
          lon_1 = MAX(lon_1,patch_2d%verts%vertex(vertex_idx,vertex_blk)%lon+pi)
          
        END DO
        if (lon_1 - lon_0 > pi) then
          swap=lon_1
          lon_1=lon_0
          lon_0=swap
        endif

        div_analytic (j,1,block) = &
  !   SS(du/dx)dxdy = -1/12 sqrt(105/2*pi) * (cos(2*lon_1) - cos(2*lon_0)) * (c0s^3(lat_1) - cos^3(lat_0))
          & -1.0_wp/12.0_wp * sqrt_105/sqrt_pi_x_2 * (cos(2.0_wp*lon_1) - cos(2.0_wp*lon_0)) * (cos(lat_1)**3 - cos(lat_0)**3) + &
  !   SS(dv/dy)dydx = -1/4 * sqrt(15/2*pi) * (sin(2*lat_1) - sin(2*lat_0)) * (sin(lon_1) - sin(lon_0)) 
          & 0.25_wp * sqrt_15/sqrt_pi_x_2 * (sin(2.0_wp*lat_1) - sin(2.0_wp*lat_0)) * (sin(lon_1)-sin(lon_0))

        div_analytic(j,1,block) = div_analytic(j,1,block) * earth_radius / patch_2d%cells%area(j,block) 
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_divAnalyticAverage_planeQuads_Hui
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE fill_vn_divAnalytic_fromVertices()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: block, j, start_index,  end_index
    INTEGER :: vertex_1_idx, vertex_1_blk, vertex_2_idx, vertex_2_blk
    INTEGER :: num_verts, neigbor
    REAL(wp) :: u, v
    INTEGER, POINTER :: vertexOfEdge_idx(:,:,:)
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    vertexOfEdge_idx => patch_2d%edges%vertex_idx
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, u, v, vertex_1_idx, &
!ICON_OMP vertex_1_blk, vertex_2_idx, vertex_2_blk) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        vertex_1_idx = patch_2d%edges%vertex_idx(j,block,1)
        vertex_1_blk = patch_2d%edges%vertex_blk(j,block,1)
        vertex_2_idx = patch_2d%edges%vertex_idx(j,block,2)
        vertex_2_blk = patch_2d%edges%vertex_blk(j,block,2)

        u = (u_vert(vertex_1_idx,vertex_1_blk) + u_vert(vertex_2_idx,vertex_2_blk)) * 0.5_wp
        v = (v_vert(vertex_1_idx,vertex_1_blk) + v_vert(vertex_2_idx,vertex_2_blk)) * 0.5_wp
        vn(j,1,block) =  u * patch_2d%edges%primal_normal(j,block)%v1 &
                       & + v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, vertex_1_idx, vertex_1_blk, neigbor, num_verts) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      div_analytic (:,1,block) = 0.0_wp
      DO j =  start_index, end_index
        num_verts = patch_2D%cells%num_edges(j,block)
        DO neigbor=1, num_verts     
          vertex_1_idx = patch_2d%cells%vertex_idx(j,block,neigbor)
          vertex_1_blk = patch_2d%cells%vertex_blk(j,block,neigbor)

          div_analytic (j,1,block) = div_analytic (j,1,block) + div_vert(vertex_1_idx,vertex_1_blk)
        END DO
        div_analytic (j,1,block) = div_analytic (j,1,block) / REAL(num_verts,wp)

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_fromVertices
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE interpolate_vn_fromVerticesEdges()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: block, j, start_index,  end_index
    INTEGER :: vertex_1_idx, vertex_1_blk, vertex_2_idx, vertex_2_blk
    INTEGER :: num_verts, neigbor
    REAL(wp) :: u, v
    INTEGER, POINTER :: vertexOfEdge_idx(:,:,:)
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    vertexOfEdge_idx => patch_2d%edges%vertex_idx
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, u, v, vertex_1_idx, &
!ICON_OMP vertex_1_blk, vertex_2_idx, vertex_2_blk) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        vertex_1_idx = patch_2d%edges%vertex_idx(j,block,1)
        vertex_1_blk = patch_2d%edges%vertex_blk(j,block,1)
        vertex_2_idx = patch_2d%edges%vertex_idx(j,block,2)
        vertex_2_blk = patch_2d%edges%vertex_blk(j,block,2)

        u = (u_vert(vertex_1_idx,vertex_1_blk) + u_vert(vertex_2_idx,vertex_2_blk)) * 0.5_wp
        v = (v_vert(vertex_1_idx,vertex_1_blk) + v_vert(vertex_2_idx,vertex_2_blk)) * 0.5_wp
        
        vn(j,1,block) =  0.75_wp * vn(j,1,block)  + &
          & 0.25_wp * (  u * patch_2d%edges%primal_normal(j,block)%v1 &
          &            + v * patch_2d%edges%primal_normal(j,block)%v2)

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE interpolate_vn_fromVerticesEdges
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  ! u = 1/4 * sqrt(105/(2*pi)) cos(2*lon) cos^2(lat) sin(lat)
  ! v = -1/2 * sqrt(15/(2*pi) * cos(lon) * cos(lat) * sin(lat)
  ! div = -1/(2*sqr(2*pi)) (sqrt(105) * sin(2*lon) * cos^2(lat) * sin(lat) + sqrt(15) * cos(lon) cos(2*lon))
  ! 
  SUBROUTINE fill_vn_divAnalytic_plane_Hui_onVertices()
  
    TYPE(t_subset_range), POINTER :: all_verts
    REAL(wp) :: lon, lat, u, v
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_verts => patch_2d%verts%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%verts%vertex(j,block)%lon
        lat = patch_2d%verts%vertex(j,block)%lat
        
        u_vert(j,block) = 0.25_wp * sqrt_105 / sqrt_pi_x_2 * cos(2.0_wp*lon) * cos(lat)**2 * sin(lat)
        v_vert(j,block) = - 0.5_wp * sqrt_15 / sqrt_pi_x_2 * cos(lon) * cos(lat) * sin(lat) 
        div_vert (j,block) = &
           & - 1.0_wp/(2._wp*sqrt_pi_x_2 * earth_radius) * (sqrt_105 * sin(2.0_wp * lon) * cos(lat)**2 * sin(lat) + &
             sqrt_15 * cos(lon) * cos(2.0_wp*lat))     
 
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_plane_Hui_onVertices
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! u = sin(lon)
  ! v = 0
  ! div = cos(lon)
  ! 
  SUBROUTINE fill_vn_divAnalytic_plane_basic1()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, u, v
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%edges%center(j,block)%lon
        lat = patch_2d%edges%center(j,block)%lat
        
        u = sin(lon)
        v = 0.0_wp 
        vn(j,1,block) =  u * patch_2d%edges%primal_normal(j,block)%v1 &
                       & + v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
!         lon = patch_2d%cells%cartesian_center(j,block)%x(1)
!         lat = patch_2d%cells%cartesian_center(j,block)%x(2)
        
        div_analytic (j,1,block) = cos(lon) / earth_radius
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_plane_basic1
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! u = 1/4 * sqrt(105/(2*pi)) cos(2*lon) cos^2(lat) sin(lat)
  ! v = -1/2 * sqrt(15/(2*pi) * cos(lon) * cos(lat) * sin(lat)
  ! div = -1/(2*sqr(2*pi)) (sqrt(105) * sin(2*lon) * cos^2(lat) * sin(lat) + sqrt(15) * cos(lon) cos(2*lon))
  ! 
  SUBROUTINE fill_vn_divAnalytic_plane_basic2_onVertices()
  
    TYPE(t_subset_range), POINTER :: all_verts
    REAL(wp) :: lon, lat, u, v
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_verts => patch_2d%verts%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%verts%vertex(j,block)%lon
        lat = patch_2d%verts%vertex(j,block)%lat
        
        u = 0.0_wp 
        v = cos(lat)  
        u_vert(j,block) = 0.0_wp
        v_vert(j,block) = cos(lat) 
        div_vert (j,block) = -sin(lat) / earth_radius
 
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

    CALL fill_vn_divAnalytic_fromVertices

  END SUBROUTINE fill_vn_divAnalytic_plane_basic2_onVertices
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! u = 0
  ! v = cos(lat)
  ! div = -sin(lat)
  ! 
  SUBROUTINE fill_vn_divAnalytic_plane_basic2()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, u, v
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%edges%center(j,block)%lon
        lat = patch_2d%edges%center(j,block)%lat
        
        u = 0.0_wp 
        v = cos(lat)  
        vn(j,1,block) =  u * patch_2d%edges%primal_normal(j,block)%v1 &
                       & + v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
!         lon = patch_2d%cells%cartesian_center(j,block)%x(1)
!         lat = patch_2d%cells%cartesian_center(j,block)%x(2)
        
        div_analytic (j,1,block) = -sin(lat) / earth_radius
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_plane_basic2
  !-------------------------------------------------------------------------

   !-------------------------------------------------------------------------
  ! u = 0
  ! v = cos(lat)
  ! integrated: Sv =  = sin(lat)
  ! div = -sin(lat)
  ! 
  SUBROUTINE fill_vn_divAnalytic_plane_basic2_accurate()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, lat_1, lat_2, v
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index
    INTEGER :: vertex_1_idx, vertex_1_blk, vertex_2_idx, vertex_2_blk
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lat, lat_1, lat_2, v, &
!ICON_OMP vertex_1_idx, vertex_1_blk, vertex_2_idx, vertex_2_blk) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        vertex_1_idx = patch_2d%edges%vertex_idx(j,block,1)
        vertex_1_blk = patch_2d%edges%vertex_blk(j,block,1)
        vertex_2_idx = patch_2d%edges%vertex_idx(j,block,2)
        vertex_2_blk = patch_2d%edges%vertex_blk(j,block,2)
      
        lat = patch_2d%edges%center(j,block)%lat
        lat_1 = patch_2d%verts%vertex(vertex_1_idx,vertex_1_blk)%lat
        lat_2 = patch_2d%verts%vertex(vertex_2_idx,vertex_2_blk)%lat
!         lon_1 = patch_2d%verts%vertex(vertex_1_idx,vertex_1_blk)%lon
!         lon_2 = patch_2d%verts%vertex(vertex_2_idx,vertex_2_blk)%lat
     
        IF (lat_1 == lat_2) THEN
          v = cos(lat_1)
        ELSE
          v = (sin(lat_2) - sin(lat_1)) * patch_2d%edges%tangent_orientation(j,block) &
            & * ABS(patch_2d%edges%dual_normal(j,block)%v2 / patch_2d%edges%dual_normal(j,block)%v1) &
            & * earth_radius / patch_2d%edges%primal_edge_length(j,block)
        ENDIF

        vn(j,1,block) = v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
!         lon = patch_2d%cells%cartesian_center(j,block)%x(1)
!         lat = patch_2d%cells%cartesian_center(j,block)%x(2)
        
        div_analytic (j,1,block) = -sin(lat) / earth_radius
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_plane_basic2_accurate
  !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
  ! u = 0
  ! v = cos(lat)
  ! div = -sin(lat)
  ! 
  ! second order term= d^2v/(dx)^2 * (n_1)^2 + 2* d^2v/d(xy) * (n_1 * n_2) + d^2v/(dy)^2 * n_2^2
  !  d^2v/(dx)^2 = =
  !  d^2v/d(xy)  = 0
  !  d^2v/(dy)^2 = -cos(lat)
  !
  ! Not tested !
  SUBROUTINE fill_vn_divAnalytic_plane_basic2_3order()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, u, v,n2
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2, secondOrderV
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v, n2, secondOrderV) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%edges%center(j,block)%lon
        lat = patch_2d%edges%center(j,block)%lat
        n2 = patch_2d%edges%dual_normal(j,block)%v2

        u = 0.0_wp 
        v = cos(lat)  
        ! 1/2 * v'' * 1/2 * 1/3 * len^3 / len^2
        secondOrderV = -cos(lat) * &
          & (n2 * patch_2d%edges%primal_edge_length(j,block) / earth_radius)**2 &
          & / 6.0_wp

        vn(j,1,block) = (v + secondOrderV) * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
!         lon = patch_2d%cells%cartesian_center(j,block)%x(1)
!         lat = patch_2d%cells%cartesian_center(j,block)%x(2)
        
        div_analytic (j,1,block) = -sin(lat) / earth_radius
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_plane_basic2_3order
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! u = 0
  ! v = 1
  ! div = 0
  ! 
  SUBROUTINE fill_vn_divAnalytic_plane_basic3()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, u, v
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%edges%center(j,block)%lon
        lat = patch_2d%edges%center(j,block)%lat
        
        u = 0.0_wp 
        v = 1.0_wp
        vn(j,1,block) =  u * patch_2d%edges%primal_normal(j,block)%v1 &
                       & + v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
!         lon = patch_2d%cells%cartesian_center(j,block)%x(1)
!         lat = patch_2d%cells%cartesian_center(j,block)%x(2)
        
        div_analytic (j,1,block) = 0.0_wp
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_plane_basic3
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! u = 1
  ! v = 0
  ! div = 0
  ! 
  SUBROUTINE fill_vn_divAnalytic_plane_basic4()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, u, v
    REAL(wp) :: sqrt_105, sqrt_15, sqrt_pi_x_2
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------
    sqrt_105 = SQRT(105.0_wp)
    sqrt_15  = SQRT(15.0_wp)
    sqrt_pi_x_2 = SQRT(pi*2.0_wp)

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%edges%center(j,block)%lon
        lat = patch_2d%edges%center(j,block)%lat
        
        u = 1.0_wp 
        v = 0.0_wp
        vn(j,1,block) =  u * patch_2d%edges%primal_normal(j,block)%v1 &
                       & + v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
        
        div_analytic (j,1,block) = 0.0_wp
      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_plane_basic4
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  ! u = 0
  ! v = -1 (=1 in spherical)
  ! div = 1/(r*cos(lat)) * sin(lat)
  SUBROUTINE fill_vn_divAnalytic_onSphere_test1()
  
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: lon, lat, u, v
    INTEGER :: block, j, start_index,  end_index
    !-----------------------------------------------------------------------
    all_edges => patch_2d%edges%all
    all_cells => patch_2d%cells%all
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat, u, v) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%edges%center(j,block)%lon
        lat = patch_2d%edges%center(j,block)%lat
        
        u = 0.0_wp 
        v = -1.0_wp
        vn(j,1,block) =  u * patch_2d%edges%primal_normal(j,block)%v1 &
                       & + v * patch_2d%edges%primal_normal(j,block)%v2

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(block,j,start_index,end_index, lon, lat) ICON_OMP_DEFAULT_SCHEDULE
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_index, end_index)
      DO j =  start_index, end_index
        lon = patch_2d%cells%center(j,block)%lon
        lat = patch_2d%cells%center(j,block)%lat
        
        div_analytic (j,1,block) =  sin(lat) / (earth_radius * cos(lat))  

      END DO
    END DO 
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE fill_vn_divAnalytic_onSphere_test1
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE test_vn_accuracy_onPlaneQuads_Hui()
    REAL(wp), POINTER :: vn_analytic(:,:,:), vn_point(:,:,:), vn_diff(:,:,:)
    TYPE(t_subset_range), POINTER :: owned_edges

    owned_edges => patch_2d%edges%owned
    !-----------------------------------------------------------------------
    vn_point     => ocean_state%p_diag%mass_flx_e
    vn_analytic  => ocean_state%p_diag%vn_pred
    vn_diff      => ocean_state%p_diag%vn_pred_ptp

    write(0,*) "-------------- test vn accuracy plane -----------------------"

    vn => vn_analytic
    CALL fill_vnAveraged_planeQuads_Hui

    vn => vn_point
    CALL fill_vn_divAnalytic_plane_Hui()
    
    vn_diff = vn_point - vn_analytic
    !-----------------------------------------------------------------------


    CALL dbg_print('vn_point',vn_point,str_module,1, in_subset=owned_edges)
    CALL dbg_print('vn_analytic',vn_analytic,str_module,1, in_subset=owned_edges)
    CALL dbg_print('vn_diff',vn_diff,str_module,1, in_subset=owned_edges)
    
    L2Diff  = L2Norm(vn_diff(:,1,:), owned_edges)
    L2DivAn = L2Norm(vn_analytic(:,1,:), owned_edges)
    L2Diff = L2Diff / L2DivAn

    LInfDiff  = LInfNorm(vn_diff(:,1,:), owned_edges)
    LInfDivAn = LInfNorm(vn_analytic(:,1,:), owned_edges)
    LInfDiff  = LInfDiff / LInfDivAn

    write(0,*) "=================================="
    write(0,*) "L2 error:",   L2Diff
    write(0,*) "LInf error:", LInfDiff
    write(0,*) "=================================="

  END SUBROUTINE test_vn_accuracy_onPlaneQuads_Hui
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

END MODULE mo_ocean_testbed_div

