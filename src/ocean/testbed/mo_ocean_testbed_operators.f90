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

MODULE mo_ocean_testbed_operators
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
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_div, timer_grad
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_math_utilities,      ONLY: vector_product
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_grid_config,         ONLY: n_dom
  USE mo_ocean_math_operators

  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_physics_types, ONLY: t_ho_params
  USE mo_sea_ice_types,       ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_ocean_surface_types, ONLY: t_ocean_surface
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff, no_primal_edges, no_dual_edges
  USE mo_ocean_testbed_coriolis, ONLY: test_nonlinear_coriolis_3d
  USE mo_ocean_testbed_div,   ONLY: test_div
  USE mo_ocean_testbed_PtP,   ONLY: test_PtP


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=12)           :: str_module    = 'oceMathOps  '  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug


  PUBLIC :: ocean_test_operators

CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_operators( namelist_filename, shr_namelist_filename, &
    & patch_3d, ocean_state, external_data,                     &
    & surface_fluxes, physics_parameters,             &
    & oceans_atmosphere_fluxes, ocean_ice, operators_coefficients)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data(n_dom)
    TYPE(t_ocean_surface)                            :: surface_fluxes
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: oceans_atmosphere_fluxes
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),TARGET,   INTENT(inout)   :: operators_coefficients

    CHARACTER(LEN=*), PARAMETER ::  method_name = "ocean_test_operators"

    SELECT CASE (test_mode) ! 100 - 999

      CASE (101)
        CALL test_nonlinear_coriolis_3d(patch_3D, ocean_state(1), operators_coefficients)
      CASE (102)
        CALL operator_test_old(patch_3D, ocean_state(1), operators_coefficients)
      CASE (103:113)
        CALL test_div(patch_3D, ocean_state(1), operators_coefficients, test_mode)
      CASE (114:115)
        CALL test_PtP(patch_3D, ocean_state(1), operators_coefficients, test_mode)
      
      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT

  END SUBROUTINE ocean_test_operators
  !-------------------------------------------------------------------------


  !-----------------------------------------------------------------------
  SUBROUTINE operator_test_old( patch_3d, ocean_state, operators_coefficients)!, vn_e, trac_c)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff),TARGET, INTENT(in)   :: operators_coefficients
    
    !Local variables
    INTEGER :: start_index, end_index
    INTEGER :: jc, je, jv, jk, jb,edge_index,edge_block
    REAL(wp) :: curl_integral(1:n_zlev), div_integral(1:n_zlev), lhs(1:n_zlev),rhs(1:n_zlev)
    REAL(wp) :: grad(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: div (nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp) :: curl(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp) :: curlgrad(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_v)    
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain, verts_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_cartesian_coordinates) :: vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)    
    REAL(wp), POINTER :: vn_e(:,:,:)!(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), POINTER :: trac_c(:,:,:)!(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    

    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    cells_in_domain => patch_2D%cells%in_domain
    verts_in_domain => patch_2D%verts%in_domain
    !-------------------------------------------------------------------------------
    vn_e  => ocean_state%p_prog(1)%vn
    trac_c=>ocean_state%p_prog(1)%tracer(:,:,:,1) 	

    grad    (1:nproma,1:n_zlev,1:patch_2D%nblks_e)=0.0_wp
    div     (1:nproma,1:n_zlev,1:patch_2D%nblks_c)=0.0_wp
    curl    (1:nproma,1:n_zlev,1:patch_2D%nblks_v)=0.0_wp
    curlgrad(1:nproma,1:n_zlev,1:patch_2D%nblks_v)=0.0_wp
    
    
    !vn_e=0.0_wp
    vn_e(1:5,1,1)=10.0_wp
    trac_c=35.0_wp
    CALL map_edges2vert_3d(patch_2D, vn_e, operators_coefficients%edge2vert_coeff_cc, vn_dual)
    
    !calculate gradient of curl and curl
    CALL grad_fd_norm_oce_3d( trac_c, patch_3D, operators_coefficients%grad_coeff, grad)
    CALL rot_vertex_ocean_3d( patch_3D, grad, vn_dual, operators_coefficients, curlgrad)
    CALL rot_vertex_ocean_3d( patch_3D, vn_e, vn_dual, operators_coefficients, curl)
    
    !domain integrated curl
    curl_integral=0.0_wp
    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, start_index, end_index)
      DO jk = 1, n_zlev
        DO jv = start_index, end_index
        curl_integral(jk) = curl_integral(jk) + curl(jv,jk,jb)*patch_2D%verts%dual_area(jv,jb)

        END DO
      END DO
    END DO
    write(0,*)'OPERATOR-TEST:curl-of-gradient:',maxval(curlgrad),minval(curlgrad)
    Do jk=1,n_zlev
    write(0,*)'OPERATOR-TEST:curl-integral:',jk,curl_integral(jk)
    END DO
    !calculate domain integrated div
    CALL div_oce_3d( vn_e, patch_3D,operators_coefficients%div_coeff, div)
    div_integral=0.0_wp
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
          DO je=1,3
          edge_index = patch_2D%cells%edge_idx(jc,jb, je)
          edge_block = patch_2D%cells%edge_blk(jc,jb,je)
          
             div_integral(jk)= div_integral(jk)+(vn_e(edge_index,jk, edge_block)* &
             &patch_2D%edges%primal_edge_length(edge_index, edge_block)/patch_2D%cells%area(jc,jb)       * &
            & patch_2D%cells%edge_orientation(jc,jb, je))* patch_2D%cells%area(jc,jb)
IF(jk==1)&
&write(123,*)'details:',jk,jc,jb,je,edge_index,edge_block,vn_e(edge_index,jk, edge_block)* &
!             &patch_2D%edges%primal_edge_length(edge_index, edge_block)        * &
            & patch_2D%cells%edge_orientation(jc,jb, je)             
          END DO
          div_integral(jk)= div_integral(jk)
IF(jk==1)&          
&write(123,*)'details2',jk,jc,jb,div_integral(jk)         
!          div_integral(jk) = div_integral(jk) + div(jc,jk,jb)*patch_2D%cells%area(jc,jb)                            
        END DO
      END DO
    END DO
    write(123,*)'-----------------------'
    DO jk=1,n_zlev
    write(0,*)'OPERATOR-TEST: div-integral:',jk,div_integral(jk)
    END DO

    !product rule for divergence
    CALL div_oce_3d( vn_e, patch_3D,operators_coefficients%div_coeff, div)
    CALL grad_fd_norm_oce_3d( trac_c, patch_3D, operators_coefficients%grad_coeff, grad)

    lhs=0.0_wp
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jk = 1, n_zlev
        DO jc = start_index, end_index    
          lhs(jk)=lhs(jk)+trac_c(jc,jk,jb)*div(jc,jk,jb)*patch_2D%cells%area(jc,jb)        
        END DO
      END DO
    END DO
    rhs=0.0_wp
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_index, end_index)
      DO jk = 1, n_zlev
        DO je = start_index, end_index
        rhs(jk)=rhs(jk)+grad(je,jk,jb)*vn_e(je,jk,jb)*patch_2D%edges%area_edge(je,jb)
        !rhs=rhs+grad(je,jk,jb)*vn_e(je,jk,jb)&
        !&*patch_2D%edges%primal_edge_length(je,jb)*patch_2D%edges%dual_edge_length(je,jb)
        END DO
      END DO
    END DO
    Do jk=1,n_zlev
      write(0,*)'OPERATOR-TEST: integration-by-parts:',jk,lhs(jk),rhs(jk),lhs(jk)+rhs(jk)          
      write(0,*)'OPERATOR-TEST: integration-by-parts2:',jk,&
      &maxval(grad(:,jk,:)),minval(grad(:,jk,:)),maxval(vn_e(:,jk,:)),minval(vn_e(:,jk,:))          
      
    END DO  
  END SUBROUTINE operator_test_old
  !-------------------------------------------------------------------------------

END MODULE mo_ocean_testbed_operators

