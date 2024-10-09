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
! #include "omp_definitions.inc"
!----------------------------

MODULE mo_ice_new_dynamics
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_run_config,          ONLY: dtime, ltimer,nsteps
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_momentum, timer_ice_interp, timer_ice_advection

  USE mo_exception,           ONLY: message, message_text, finish

! USE mo_grid_config,         ONLY: l_limited_area, n_dom   ! for now sea-ice works on global domain-only
  USE mo_parallel_config,     ONLY: nproma
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_patch_array_mult
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_impl_constants,      ONLY: sea_boundary,sea, land, boundary

  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_dynamics_config,     ONLY: nold
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_nml,           ONLY: atm_pressure_included_in_icedyn, ssh_in_icedyn_type, &
    &                              vert_cor_type,&
    &                               ice_free_drift_only,ice_laplace_dynamics, &
    &                               ice_stabilization
  USE mo_ocean_surface_types, ONLY: t_atmos_for_ocean, t_ocean_surface
  USE mo_physical_constants,  ONLY: grav, rho_ref, sfc_press_pascal, rhoi, rhos
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_sea_ice_nml,         ONLY: i_ice_advec, pstar, n_ice_iter, delta_min, Cd_io, Cd_ia
  USE mo_ice_fem_advection,   ONLY: fct_ice_solve, ice_TG_rhs
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_math_utilities,      ONLY: gvec2cvec, cvec2gvec
  USE mo_scalar_product,      ONLY: map_cell2edges_3D, map_edges2cell_3D,map_edges2edges_viavert_3d
  USE mo_icon_interpolation_scalar, ONLY: verts2cells_scalar
  USE mo_ice_fem_interpolation, ONLY: map_edges2verts, map_verts2edges,                 &
                                      gvec2cvec_c_2d, cvec2gvec_c_2d,                   &
                                      rotate_cvec_v, gvec2cvec_v_fem, cvec2gvec_v_fem,  &
                                      cells2verts_scalar_seaice

  USE mo_ocean_math_operators,  ONLY: div_oce_3D

  USE mo_ocean_nml,          ONLY: n_zlev


  USE mo_name_list_output,       ONLY: write_name_list_output
  
  IMPLICIT NONE

  PUBLIC  :: ice_new_dynamics
  PUBLIC  :: interface_boundary_edge_marker
  PUBLIC  :: interface_boundary_cell_marker



  CHARACTER(len=12)           :: str_module    = 'NewIceDyn'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1         ! Level of detail for 1 line debug

CONTAINS

!------------------------------------------------------------------------
!
!>  Wrapper for the call to the AWI FEM ice model
!!  We first remap the neccesary inputs, then call the momentum solver (EVPdynamics)
!!  and map the resulting velocity onto edges and cell centres.
!!

  SUBROUTINE ice_new_dynamics( p_patch_3D, p_ice, p_os, p_as, atmos_fluxes, p_op_coeff, p_oce_sfc)

  
    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)     :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff
    TYPE(t_atmos_for_ocean),  INTENT(INOUT)  :: p_as
    TYPE(t_subset_range),     POINTER        :: all_cells
    TYPE(t_subset_range), POINTER            :: all_edges
    TYPE(t_subset_range), POINTER            :: owned_cells
    TYPE(t_ocean_surface),    INTENT(INOUT)  :: p_oce_sfc
    TYPE(t_subset_range), POINTER :: edges_in_domain


    ! Local variables
    TYPE(t_patch), POINTER :: patch_2D
   ! REAL(wp),      POINTER :: ssh(:,:) ! sea surface height (input only)         [m]
    !REAL(wp), ALLOCATABLE  :: ssh_reduced(:,:) ! reduced sea surface height to take slp coupling into account     [m]
    REAL(wp):: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)::  s11(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s22(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s21(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s12(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s31(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &s32(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &sigma_I(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &sigma_II(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &zeta_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &zeta_stabi(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &cell_area_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &x1_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &x2_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks),&
         &x3_c(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)



    REAL(wp) ::  boundary_edge_marker(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  mass(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  atm_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  atm_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  u_ocean_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  u_ocean_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  Au_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  Au_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  h_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  A_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  s_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  S_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  S_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  S_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  S_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  ice_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  ice_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  zeta_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)


    TYPE(t_cartesian_coordinates) :: p_tau_n_c(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: ocean_c(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
 REAL(wp)                      :: tau_n(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
     Real(wp)                      :: ocean_n(nproma,p_patch_3D%p_patch_2D(1)%nblks_e) 
     TYPE(t_cartesian_coordinates) :: p_tau_n_dual(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
     TYPE(t_cartesian_coordinates) :: p_ocean_n_dual(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)

   INTEGER ::time_iter,outer_iter, time, edge_block_i,start_index,end_index,edge_index_i, cell_block, cell_index, neigbor,&
        &edge_block_1,edge_block_2,edge_block_3,edge_index_1,edge_index_2,edge_index_3, inner_iter&
       &,cell_index1,cell_index2,cell_block1,cell_block2,doy,vert_index,vert_block,vert_index_1,vert_index_2,vert_block_1,vert_block_2 !, vtkeverystep, vtk_iter
   



   REAL(wp) :: a,x,y,z,x_0,y_0, nx,ny,tx,ty,Cda,Cdw, alpha_evp, beta_evp,nix,niy,niz,tix,tiy,tiz&
        &,T,e,Delta, nx_cc,ny_cc, tx_cc,ty_cc, e1,e2,e3,s,R,nz,tz, nix1,nix2,nix3,niy2,niy3,niy1,O1,O2,O3,x123,grad_x,grad_y
   REAL(wp) ::  Oi, Oj, ei, e11, e22, e21, e12, e31,e32,P, zeta,eta, x1,x2,x3, y1,y2,y3, tix1,tix2,tix3,tiy2,tiy3,tiy1,ux,uy,uz,&
        &sc_pn, l_sn, nix_1,niy_1,niz_1, tix_1,tiy_1,tiz_1,nix_l, niy_l, tix_l,tiy_l,U,V,W, U_ana, V_ana,W_ana, ice_z, u_rel,v_rel,w_rel,ocean_z,&
        &atop_n,atop_t,abot_n,abot_t,btop_n,btop_t,bbot_n,bbot_t,aaa,bbb,x_shift,y_shift,L,diff_n,diff_t,vdw,tp,U_xx,U_yy,atm_u,atm_v,mass_ice
   REAL(wp) :: alpha,vmax, ws,wx,wy,mx,my,rw, C_imp, vda, weight      
       
   CHARACTER(*), PARAMETER :: &
          method_name = 'mo_ice_new_dynamics:ice_new_dynamics'
   !--------------------------------------------------------------------------------------------------
   
    patch_2D => p_patch_3D%p_patch_2D(1)
    all_cells     => p_patch_3d%p_patch_2d(1)%cells%all
    all_edges =>     p_patch_3d%p_patch_2d(1)%edges%all
    owned_cells => patch_2D%cells%owned
    edges_in_domain   => patch_2D%edges%in_domain


    ocean_c(:,:)%x(1)=0.0_wp
    ocean_c(:,:)%x(2)=0.0_wp
    ocean_c(:,:)%x(3)=0.0_wp
    p_tau_n_c(:,:)%x(1)=0.0_wp
    p_tau_n_c(:,:)%x(2)=0.0_wp
    p_tau_n_c(:,:)%x(3)=0.0_wp
    boundary_cell_marker(:,:,:)=0.0_wp
    boundary_edge_marker(:,:)=0.0_wp
    mass(:,:)=0.0_wp
    atm_n(:,:)=0.0_wp
    atm_t(:,:)=0.0_wp
    u_ocean_n(:,:)=0.0_wp
    u_ocean_t(:,:)=0.0_wp
    ice_x(:,:)=0.0_wp
    ice_y(:,:)=0.0_wp
    Au_n(:,:)=0.0_wp
    Au_t(:,:)=0.0_wp
    A_e(:,:)=0.0_wp
    h_e(:,:)=0.0_wp
    s_e(:,:)=0.0_wp
    S_x(:,:)=0.0_wp
    S_y(:,:)=0.0_wp
    S_t(:,:)=0.0_wp
    S_n(:,:)=0.0_wp
    zeta_e(:,:)=0.0_wp



    Cdw=rho_ref*Cd_io  ! reference_density * Ice-ocean drag coefficient
    Cda=1.3_wp*Cd_ia   ! air density * Ice-atmosphere drag coefficient



    alpha_evp=1500 !1500

    beta_evp=alpha_evp
    
    doy=0


!     CALL dbg_print('start 0 vn', p_ice%vn_e ,str_module, 2, in_subset=patch_2D%edges%owned)
    
    DO cell_block = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, cell_block, start_index, end_index)
!       DO cell_index = start_index, end_index
          s11(:,:,cell_block)=0.0_wp
          s12(:,:,cell_block)=0.0_wp
          s22(:,:,cell_block)=0.0_wp
          s21(:,:,cell_block)=0.0_wp
          s32(:,:,cell_block)=0.0_wp
          s31(:,:,cell_block)=0.0_wp
          sigma_I(:,:,cell_block)=0.0_wp
          sigma_II(:,:,cell_block)=0.0_wp
          zeta_c(:,:,cell_block)=0.0_wp
          zeta_stabi(:,:,cell_block)=0.0_wp
          p_ice%Delta(:,cell_block)=0.0
          cell_area_c(:,:,cell_block)=0.0
          x1_c(:,:,cell_block)=0.0
          x2_c(:,:,cell_block)=0.0
          x3_c(:,:,cell_block)=0.0
!       ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block

!Boundary_cell_marker and boundary_edge_marker are used to marke the cells and edge where sea ice is present. This changes in every time step  
    CALL interface_boundary_cell_marker(boundary_cell_marker, p_patch_3D, p_ice)
    CALL interface_boundary_edge_marker(boundary_edge_marker,boundary_cell_marker, p_patch_3D, p_ice)
    CALL cell_area(cell_area_c, p_patch_3D)  ! why in the timeloop ?this can be moved to init file, this is needed for the integration of the stress tensor
    CALL init_mass_matrix(mass,cell_area_c,p_patch_3D) !this can be moved to init file, This is needed for the mass lumping in the momentum equation 

!     CALL dbg_print('start 0.1 vn', p_ice%vn_e ,str_module, 2, in_subset=patch_2D%edges%owned)

    !Initialize wind/ocean

    !The loop calculates the midpoint of the flat traingle. This is needed to project the surface normal and tangetial vectors into the tangential plane that contains the midpoint x_c=(x1_c,x2_c,x3_c). This is needed to calculate the stress tensor
       DO cell_block = all_cells%start_block, owned_cells%end_block
          CALL get_index_range(all_cells, cell_block, start_index, end_index)
          DO cell_index = start_index, end_index
             
             DO neigbor=1,3 !no_primal_edges                                                                                                                                                                           
                edge_index_i =patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
                edge_block_i = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)
                
                x = patch_2D%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)
                y = patch_2D%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)
                z = patch_2D%edges%cartesian_center(edge_index_i,edge_block_i)%x(3)
                
                x1_c(cell_index,1,cell_block)=x1_c(cell_index,1,cell_block) + 1.0_wp/3.0_wp * x
                x2_c(cell_index,1,cell_block)=x2_c(cell_index,1,cell_block) + 1.0_wp/3.0_wp * y
                x3_c(cell_index,1,cell_block)=x3_c(cell_index,1,cell_block) + 1.0_wp/3.0_wp * z
            
                
            ENDDO

            IF ( boundary_cell_marker(cell_index,1,cell_block)>1.0_wp) THEN
              WRITE(message_text,'(a,i10)') 'assert message cellmarker',MAXVAL(boundary_cell_marker)
              CALL finish(method_name,  message_text)
            ENDIF

             x123 = sqrt(x1_c(cell_index,1,cell_block)*x1_c(cell_index,1,cell_block) + &
                  &x2_c(cell_index,1,cell_block)*x2_c(cell_index,1,cell_block) + &
                  &x3_c(cell_index,1,cell_block)*x3_c(cell_index,1,cell_block))
          
             x1_c(cell_index,1,cell_block)=x1_c(cell_index,1,cell_block)/x123
             x2_c(cell_index,1,cell_block)=x2_c(cell_index,1,cell_block)/x123
             x3_c(cell_index,1,cell_block)=x3_c(cell_index,1,cell_block)/x123
          ENDDO
       ENDDO
       
      CALL sync_patch_array_mult(sync_c, patch_2D, 3, x1_c, x2_c,x3_c) 
       
!        CALL sync_patch_array(SYNC_C, patch_2D, x1_c)
!        CALL sync_patch_array(SYNC_C, patch_2D, x2_c)
!        CALL sync_patch_array(SYNC_C, patch_2D, x3_c)

    !**************************************************************
    ! (1) Convert lat-lon wind ocean  stress to cartesian coordinates
    !**************************************************************
       CALL  gvec2cvec_c_2d(p_patch_3D, atmos_fluxes%stress_x, atmos_fluxes%stress_y, p_tau_n_c)
!       CALL dbg_print('start 0.2 vn', p_ice%vn_e ,str_module, 2, in_subset=patch_2D%edges%owned)
                      
!        CALL sync_patch_array(SYNC_C, patch_2D, p_tau_n_c(:,:)%x(1))
!        CALL sync_patch_array(SYNC_C, patch_2D, p_tau_n_c(:,:)%x(2))
!        CALL sync_patch_array(SYNC_C, patch_2D, p_tau_n_c(:,:)%x(3))

       CALL gvec2cvec_c_2d(p_patch_3D, p_os%p_diag%u(:,1,:), p_os%p_diag%v(:,1,:), ocean_c)
 
!        CALL dbg_print('start 0.21 vn', p_ice%vn_e ,str_module, 2, in_subset=patch_2D%edges%owned)
 
!       CALL sync_patch_array(SYNC_C, patch_2D, ocean_c(:,:)%x(1))
!       CALL sync_patch_array(SYNC_C, patch_2D, ocean_c(:,:)%x(2))
!       CALL sync_patch_array(SYNC_C, patch_2D, ocean_c(:,:)%x(3))
       
!      CALL dbg_print('start 0.3 vn', p_ice%vn_e ,str_module, 2, in_subset=patch_2D%edges%owned)
       
    !**************************************************************
    ! (2) Interpolate 3D wind stress and ocean  from cell centers to edges
    !calculating normal and tangent component of wind/ocean
    !**************************************************************

#ifdef NAGFOR
       tau_n = 0.0_wp
       ocean_n=0.0_wp
#endif
        CALL map_cell2edges_3D(p_patch_3D, p_tau_n_c, tau_n ,p_op_coeff, 1)
        CALL sync_patch_array(SYNC_E, patch_2D, tau_n)
  
        !**************************************************************                                                                                                               
        ! (3) Interpolate 3D wind stress from normal value to vertices 3D                                                                                                                       
        !**************************************************************                                                                                                               
#ifdef NAGFOR
    p_tau_n_dual(:,:)%x(1) = 0.0_wp
    p_tau_n_dual(:,:)%x(2) = 0.0_wp
    p_tau_n_dual(:,:)%x(3) = 0.0_wp

    p_ocean_n_dual(:,:)%x(1) = 0.0_wp
    p_ocean_n_dual(:,:)%x(2) = 0.0_wp
    p_ocean_n_dual(:,:)%x(3) = 0.0_wp
#endif


    DO edge_block_i = all_edges%start_block, all_edges%end_block
       CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
       DO edge_index_i =  start_index, end_index
          ocean_n(edge_index_i,edge_block_i)=p_os%p_prog(nold(1))%vn(edge_index_i,1,edge_block_i)
       ENDDO
    ENDDO
       
!     CALL sync_patch_array(SYNC_E, patch_2D, ocean_n)
!     CALL dbg_print('start 1 vn', p_ice%vn_e ,str_module, 2, in_subset=patch_2D%edges%owned)
           
    CALL map_edges2verts(patch_2D, tau_n, p_op_coeff%edge2vert_coeff_cc, p_tau_n_dual)
    CALL map_edges2verts(patch_2D, ocean_n, p_op_coeff%edge2vert_coeff_cc, p_ocean_n_dual)

!     CALL sync_patch_array(SYNC_V, patch_2D, p_tau_n_dual(:,:)%x(1))
!     CALL sync_patch_array(SYNC_V, patch_2D, p_tau_n_dual(:,:)%x(2))
!     CALL sync_patch_array(SYNC_V, patch_2D, p_tau_n_dual(:,:)%x(3))
         
    !**************************************************************
    ! (4.1) Interpolate 3D wind stress from vertices to 3D vector on edges to calculate normal and tangential componantes.
    !The step via vertices is nesessarty to get a smooth representation of the wind field.
    !(4.2) The same is done for the ocean                                                        
    !**************************************************************   
    DO edge_block_i = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, edge_block_i, start_index, end_index)
      DO edge_index_i =  start_index, end_index

        IF ( boundary_edge_marker(edge_index_i,edge_block_i)>1.0_wp) THEN
          WRITE(message_text,'(a,i10)') 'assert message edgemarker',MAXVAL(boundary_edge_marker)
          CALL finish(method_name,  message_text)
        ENDIF

        ice_x(edge_index_i,edge_block_i)=p_ice%vn_e(edge_index_i,edge_block_i) &
            *boundary_edge_marker(edge_index_i,edge_block_i)
        ice_y(edge_index_i,edge_block_i)=p_ice%vt_e(edge_index_i,edge_block_i) &
            *boundary_edge_marker(edge_index_i,edge_block_i)

        nix=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)
        niy=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)
        niz=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)

        tix=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)
        tiy=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)
        tiz=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)
        
        vert_index_1= patch_2D%edges%vertex_idx(edge_index_i,edge_block_i,1)
        vert_index_2=patch_2D%edges%vertex_idx(edge_index_i,edge_block_i,2)
        vert_block_1=patch_2D%edges%vertex_blk(edge_index_i,edge_block_i,1)
        vert_block_2=patch_2D%edges%vertex_blk(edge_index_i,edge_block_i,2)

        atm_n(edge_index_i,edge_block_i)=0.0_wp
        atm_t(edge_index_i,edge_block_i)=0.0_wp
        u_ocean_n(edge_index_i,edge_block_i)=0.0_wp
        u_ocean_t(edge_index_i,edge_block_i)=0.0_wp
        
        IF (p_patch_3d%lsm_e(edge_index_i,1,edge_block_i) <= sea_boundary) THEN
            atm_n(edge_index_i,edge_block_i)=&
                  &( p_tau_n_dual(vert_index_1,vert_block_1)%x(1)+ p_tau_n_dual(vert_index_2,vert_block_2)%x(1))*0.5_wp*nix+&
                  &(p_tau_n_dual(vert_index_1,vert_block_1)%x(2)+ p_tau_n_dual(vert_index_2,vert_block_2)%x(2))*0.5_wp*niy+&
                  &(p_tau_n_dual(vert_index_1,vert_block_1)%x(3)+ p_tau_n_dual(vert_index_2,vert_block_2)%x(3))*0.5_wp*niz
      
            
              atm_t(edge_index_i,edge_block_i)=&
                  &(p_tau_n_dual(vert_index_1,vert_block_1)%x(1)+ p_tau_n_dual(vert_index_2,vert_block_2)%x(1))*0.5_wp*tix+&
                  &(p_tau_n_dual(vert_index_1,vert_block_1)%x(2)+ p_tau_n_dual(vert_index_2,vert_block_2)%x(2))*0.5_wp*tiy+&
                  &(p_tau_n_dual(vert_index_1,vert_block_1)%x(3)+ p_tau_n_dual(vert_index_2,vert_block_2)%x(3))*0.5_wp*tiz
                  

              u_ocean_n(edge_index_i,edge_block_i)=&
                  &(p_ocean_n_dual(vert_index_1,vert_block_1)%x(1)+ p_ocean_n_dual(vert_index_2,vert_block_2)%x(1))*0.5_wp*nix+&
                  &(p_ocean_n_dual(vert_index_1,vert_block_1)%x(2)+ p_ocean_n_dual(vert_index_2,vert_block_2)%x(2))*0.5_wp*niy+&
                  &(p_ocean_n_dual(vert_index_1,vert_block_1)%x(3)+ p_ocean_n_dual(vert_index_2,vert_block_2)%x(3))*0.5_wp*niz


              u_ocean_t(edge_index_i,edge_block_i)=&
                  &(p_ocean_n_dual(vert_index_1,vert_block_1)%x(1)+ p_ocean_n_dual(vert_index_2,vert_block_2)%x(1))*0.5_wp*tix+&
                  &(p_ocean_n_dual(vert_index_1,vert_block_1)%x(2)+ p_ocean_n_dual(vert_index_2,vert_block_2)%x(2))*0.5_wp*tiy+&
                  &(p_ocean_n_dual(vert_index_1,vert_block_1)%x(3)+ p_ocean_n_dual(vert_index_2,vert_block_2)%x(3))*0.5_wp*tiz
        

        ENDIF

      ENDDO
    ENDDO
       
!     CALL sync_patch_array_mult(sync_e, patch_2D, 4, atm_n, atm_t, u_ocean_n,  u_ocean_t)
          
    CALL sync_patch_array(SYNC_E, patch_2D, atm_n)
    CALL sync_patch_array(SYNC_E, patch_2D, atm_t)
    CALL sync_patch_array(SYNC_E, patch_2D, u_ocean_n)
    CALL sync_patch_array(SYNC_E, patch_2D, u_ocean_t)
 
!        CALL dbg_print('start 2 vn', p_ice%vn_e ,str_module, 2, in_subset=patch_2D%edges%owned)
   
!       WRITE(message_text,'(a,3i10)') 'ocean/atm/ice_old',MAXVAL(u_ocean_n) , MAXVAL(atm_n), MAXVAL(ice_x)
!        write(0,*) "ocean/atm/ice_old:", p_ice%ice_iter, maxval(u_ocean_n) , maxval(atm_n), maxval(ice_x)

    
    DO outer_iter=1, n_ice_iter
        
!        write(0,*) "outer_iter=", outer_iter
       !This loop calculates the stresses for the ice 
       DO cell_block = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, cell_block, start_index, end_index)
          DO cell_index = start_index, end_index

            IF (p_patch_3d%lsm_c(cell_index,1,cell_block) <= sea_boundary) THEN

               IF (boundary_cell_marker(cell_index,1,cell_block)>0.0) THEN
                  !here the strain rates are computed
                CALL compute_2Dvector_grad(cell_index,cell_block,e11,e12,e21,e22,cell_area_c(cell_index,1,cell_block),&
                     &x1_c(cell_index,1,cell_block),x2_c(cell_index,1,cell_block),x3_c(cell_index,1,cell_block),p_patch_3D, p_ice)

                p_ice%Delta(cell_index,cell_block)=dsqrt(delta_min*delta_min+ e12*e12 + 1.25_wp*(e11*e11+e22*e22)+1.5_wp*e11*e22)

                P=p_ice%hi(cell_index,1,cell_block)*Pstar*EXP(-20.0_wp*(1.0_wp-p_ice%conc(cell_index,1,cell_block)))

                zeta_c(cell_index,1,cell_block)=P/(2*max(p_ice%Delta(cell_index,cell_block),0.000000002_wp))!*boundary_cell_marker(cell_index,1,cell_block)

                zeta_stabi(cell_index,1,cell_block)=P/(2* max(p_ice%Delta(cell_index,cell_block),0.000000002_wp*10.0_wp))!*&
                    ! &boundary_cell_marker(cell_index,1,cell_block)

                  !this updates the stresses   
                Call mEVP(alpha_evp,e11,e12,e21,e22,s11(cell_index,1,cell_block),s12(cell_index,1,cell_block),&
                     &s21(cell_index,1,cell_block),s22(cell_index,1,cell_block),sigma_I(cell_index,1,cell_block),&
                     &sigma_II(cell_index,1,cell_block),zeta_c(cell_index,1,cell_block),P)
             endif
           ENDIF

          ENDDO ! cell_index = start_index, end_index

       ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
       
!         CALL sync_patch_array(SYNC_C, patch_2D, s11(:,:,:))
!         CALL sync_patch_array(SYNC_C, patch_2D, s12(:,:,:))
!         CALL sync_patch_array(SYNC_C, patch_2D, s21(:,:,:))
!         CALL sync_patch_array(SYNC_C, patch_2D, s22(:,:,:))
!         CALL sync_patch_array(SYNC_C, patch_2D, sigma_I(:,:,:))
!         CALL sync_patch_array(SYNC_C, patch_2D, sigma_II(:,:,:))
!         CALL sync_patch_array(SYNC_C, patch_2D, zeta_c(:,:,:))
!         CALL sync_patch_array(SYNC_C, patch_2D, zeta_stabi(:,:,:))
!         CALL sync_patch_array(SYNC_C, patch_2D, p_ice%Delta)
        
       !This loop averages the cell values to the edges.This is needed for the momentum equation
       DO edge_block_i = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, edge_block_i, start_index, end_index)
          DO edge_index_i =  start_index, end_index
             
             zeta_e(edge_index_i,edge_block_i)=0.0_wp
             h_e(edge_index_i,edge_block_i)=0.0_wp
             A_e(edge_index_i,edge_block_i)=0.0_wp
             s_e(edge_index_i,edge_block_i)=0.0_wp
             
           IF (p_patch_3d%lsm_e(edge_index_i,1,edge_block_i) <= sea_boundary) THEN
             cell_index1 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,1)
             cell_block1 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,1)
             cell_index2 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,2)
             cell_block2 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,2)

              if (boundary_cell_marker(cell_index1,1,cell_block1)+boundary_cell_marker(cell_index2,1,cell_block2)>1.0_wp) then
                weight=0.5_wp
             else
                weight=1.0_wp
             endif


             
             zeta_e(edge_index_i,edge_block_i)=weight*boundary_edge_marker(edge_index_i,edge_block_i)*(zeta_stabi(cell_index1,1,cell_block1)+zeta_stabi(cell_index2,1,cell_block2))
             
             h_e(edge_index_i,edge_block_i)=weight*boundary_edge_marker(edge_index_i,edge_block_i)*(p_ice%hi(cell_index1,1,cell_block1)+&
                  &p_ice%hi(cell_index2,1,cell_block2))
             
             
             A_e(edge_index_i,edge_block_i)=weight*boundary_edge_marker(edge_index_i,edge_block_i)*(p_ice%conc(cell_index1,1,cell_block1)+p_ice%conc(cell_index2,1,cell_block2))


             s_e(edge_index_i,edge_block_i)=weight*boundary_edge_marker(edge_index_i,edge_block_i)*( p_ice%hs(cell_index1,1,cell_block1)+&
            
                  &p_ice%hs(cell_index2,1,cell_block2))! *p_ice%conc(cell_index2,1,cell_block2))

             
          ENDIF

          ENDDO
       ENDDO
       
      IF (ice_free_drift_only) THEN
  
        Au_n=0.0_wp
        Au_t=0.0_wp

        S_n=0.0_wp
        S_t=0.0_wp
      ELSE
      CALL sync_patch_array(SYNC_E, patch_2D, zeta_e)
      CALL sync_patch_array(SYNC_E, patch_2D, h_e)
      CALL sync_patch_array(SYNC_E, patch_2D, A_e)
      CALL sync_patch_array(SYNC_E, patch_2D, s_e)
!       IF ( outer_iter==2) &
!        CALL finish('sea-ice dynamics','Test sync ok')
       
       Au_n=0.0_wp
       Au_t=0.0_wp

       !This function computes div(sigma) 
         CALL compute_sigma(Au_n,Au_t,boundary_cell_marker,h_e,A_e,s_e,s11,s12,s21,s22,cell_area_c,x1_c,x2_c,x3_c,p_patch_3D,p_ice,p_os)
  
       CALL sync_patch_array(SYNC_E, patch_2D, Au_n)
       CALL sync_patch_array(SYNC_E, patch_2D, Au_t)
!        CALL dbg_print(' Au_n '  , Au_n , str_module, 2, in_subset=patch_2D%cells%owned)
!        CALL dbg_print(' Au_t '  , Au_t , str_module, 2, in_subset=patch_2D%cells%owned)
       
!        IF ( outer_iter==2) &
!         CALL finish('sea-ice dynamics','Test sync ok')
 
       S_x=0.0_wp
       S_y=0.0_wp
         S_n=0.0_wp
         S_t=0.0_wp
         IF ( ice_stabilization ) THEN
!        CALL sync_patch_array(SYNC_E, patch_2D, p_ice%vn_e)
!        CALL sync_patch_array(SYNC_E, patch_2D, p_ice%vt_e)
!        CALL sync_patch_array(SYNC_C, patch_2D,x1_c)
!        CALL sync_patch_array(SYNC_C, patch_2D,x2_c)
!        CALL sync_patch_array(SYNC_C, patch_2D,x3_c)
!        CALL sync_patch_array(SYNC_C, patch_2D,boundary_cell_marker)
       
!        IF ( outer_iter==2) &
!          CALL finish('sea-ice dynamics','Test sync ok')
        
       Call  Stabilization_sum(S_x,S_y,boundary_cell_marker,x1_c,x2_c,x3_c,p_ice,p_patch_3D)
      
       CALL sync_patch_array(SYNC_E, patch_2D, S_x)
       CALL sync_patch_array(SYNC_E, patch_2D, S_y)
!        IF ( outer_iter==2) &
!          CALL finish('sea-ice dynamics','Test sync ok')

       
       DO edge_block_i = all_edges%start_block, all_edges%end_block
         CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
         DO edge_index_i =  start_index, end_index
           S_x(edge_index_i,edge_block_i)=S_x(edge_index_i,edge_block_i)*boundary_edge_marker(edge_index_i,edge_block_i)
           S_y(edge_index_i,edge_block_i)=S_y(edge_index_i,edge_block_i)*boundary_edge_marker(edge_index_i,edge_block_i)
         ENDDO
       ENDDO

       S_n=0.0_wp
       S_t=0.0_wp
       CALL  Stabilization(S_x,S_y,S_n,S_t,boundary_cell_marker,x1_c,x2_c,x3_c,zeta_e,p_patch_3D)
       
       CALL sync_patch_array(SYNC_E, patch_2D, S_n)
       CALL sync_patch_array(SYNC_E, patch_2D, S_t)
!        IF ( outer_iter==2) &
!          CALL finish('sea-ice dynamics','Test sync ok')
         ENDIF ! ice_stabilization

       ENDIF ! not ice_free_drift_only

       !solve mEVP
        DO edge_block_i = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
          DO edge_index_i =  start_index, end_index
             
             
             if(p_patch_3D%lsm_e(edge_index_i,1,edge_block_i) <= sea_boundary ) THEN
             if(A_e(edge_index_i,edge_block_i)>0.01_wp) then
                
             x = patch_2D%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)
             y = patch_2D%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)
             z = patch_2D%edges%cartesian_center(edge_index_i,edge_block_i)%x(3)

             nix=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)
             niy=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)
             niz=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)

             tix=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)
             tiy=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)
             tiz=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)

                         
             
             !vwd update ocean stress
             diff_n= u_ocean_n(edge_index_i,edge_block_i)-p_ice%vn_e(edge_index_i,edge_block_i)
             diff_t= u_ocean_t(edge_index_i,edge_block_i)-p_ice%vt_e(edge_index_i,edge_block_i)
             
             vdw=sqrt(diff_n*diff_n+diff_t*diff_t)*boundary_edge_marker(edge_index_i,edge_block_i)

             mass_ice=MAX(rhoi*h_e(edge_index_i,edge_block_i)+rhos*s_e(edge_index_i,edge_block_i),9.0_wp)
             C_imp=1.0_wp/(beta_evp+1.0_wp+vdw*Cdw*dtime/mass_ice)

               p_ice%vn_e(edge_index_i,edge_block_i)=C_imp*ice_x(edge_index_i,edge_block_i)+&
                  &C_imp*beta_evp*p_ice%vn_e(edge_index_i,edge_block_i)+&
                  &dtime/mass_ice*C_imp*&
                  & (A_e(edge_index_i,edge_block_i)*atm_n(edge_index_i,edge_block_i)+&
                  &A_e(edge_index_i,edge_block_i)*vdw*Cdw*u_ocean_n(edge_index_i,edge_block_i)&
                  &-Au_n(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i))&
                  &+C_imp*dtime/mass_ice*S_n(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i)&
                  &+patch_2D%edges%f_e(edge_index_i,edge_block_i)/dtime*(p_ice%vt_e(edge_index_i,edge_block_i)&
                  &-MERGE(u_ocean_t(edge_index_i,edge_block_i), 0.0_wp, ssh_in_icedyn_type == 2))

               p_ice%vt_e(edge_index_i,edge_block_i)= C_imp*ice_y(edge_index_i,edge_block_i)+&
                  &C_imp*beta_evp*p_ice%vt_e(edge_index_i,edge_block_i)+&
                  &C_imp*dtime/mass_ice*&
                  &(A_e(edge_index_i,edge_block_i)*atm_t(edge_index_i,edge_block_i)+&
                  & A_e(edge_index_i,edge_block_i)*u_ocean_t(edge_index_i,edge_block_i)&
                  &-Au_t(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i))&
                  &+C_imp*dtime/mass_ice*S_t(edge_index_i,edge_block_i)/mass(edge_index_i,edge_block_i)&
                  &-patch_2D%edges%f_e(edge_index_i,edge_block_i)/dtime*(p_ice%vn_e(edge_index_i,edge_block_i)&
                  &-MERGE(u_ocean_n(edge_index_i,edge_block_i),0.0_wp,ssh_in_icedyn_type == 2))


                  IF ( h_e(edge_index_i,edge_block_i)==0.0_wp) THEN
                    WRITE(message_text,'(a,i10)') 'assert message h_e',h_e(edge_index_i,edge_block_i)
                    CALL finish(method_name,  message_text)
                  ENDIF

               endif
            endif

            if(A_e(edge_index_i,edge_block_i)<=0.01_wp) then
               p_ice%vn_e(edge_index_i,edge_block_i)=0.0_wp
               p_ice%vt_e(edge_index_i,edge_block_i)=0.0_wp
            endif

          if(p_patch_3D%lsm_e(edge_index_i,1,edge_block_i) >sea_boundary ) THEN
             p_ice%vn_e(edge_index_i,edge_block_i)=0.0_wp
             p_ice%vt_e(edge_index_i,edge_block_i)=0.0_wp
          endif
          ENDDO
       ENDDO

       CALL sync_patch_array(SYNC_E, patch_2D, p_ice%vn_e)
       CALL sync_patch_array(SYNC_E, patch_2D, p_ice%vt_e)

!        IF ( outer_iter==2) &
!          CALL finish('sea-ice dynamics','Test sync ok')
          
    ENDDO !outer iter
    
!     CALL dbg_print('end outloop vn', p_ice%vn_e ,str_module, 2, in_subset=patch_2D%edges%owned)
  
!     CALL finish('sea-ice dynamics','Test sync ok')
!     CALL sync_patch_array(SYNC_E, patch_2D, p_ice%vn_e)
!     CALL sync_patch_array(SYNC_E, patch_2D, p_ice%vt_e)
!     CALL finish('sea-ice dynamics','Test sync ok')
   
  ! write(0,*) "ice_old/ice_up:",  maxval(ice_x), maxval( p_ice%vn_e)
 !write(0,*) "huhu endout"

      !visualisierung
      
      DO cell_block = owned_cells%start_block, owned_cells%end_block
         CALL get_index_range(owned_cells, cell_block, start_index, end_index)
         DO cell_index = start_index, end_index
            DO neigbor=1,3 !no_primal_edges
               edge_index_i = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
               edge_block_i = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)
                              
               x = patch_2D%edges%cartesian_center(edge_index_i,edge_block_i)%x(1)
               y = patch_2D%edges%cartesian_center(edge_index_i,edge_block_i)%x(2)
               z = patch_2D%edges%cartesian_center(edge_index_i,edge_block_i)%x(3)

               nix=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)
               niy=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)
               niz=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)
               
               tix=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)
               tiy=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)
               tiz=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)

               ice_x(edge_index_i,edge_block_i)=boundary_edge_marker(edge_index_i,edge_block_i)*&
                    &(u_ocean_t(edge_index_i,edge_block_i)*tix+u_ocean_n(edge_index_i,edge_block_i)*nix)   
               ice_y(edge_index_i,edge_block_i)=boundary_edge_marker(edge_index_i,edge_block_i)*&
                    & (p_ice%vn_e(edge_index_i,edge_block_i)*nix+p_ice%vt_e(edge_index_i,edge_block_i)*tix)
               
                
             ENDDO
          ENDDO
       ENDDO


!     CALL dbg_print('end vn', p_ice%vn_e ,str_module, 2, in_subset=patch_2D%edges%owned)


  END SUBROUTINE ice_new_dynamics





  SUBROUTINE mEVP( alpha_evp,e11,e12,e21,e22,s11,s12,s21,s22,sigma_I,sigma_II,zeta,P)

    REAL(wp), INTENT(in) :: e11
    REAL(wp), INTENT(in) :: e12
    REAL(wp), INTENT(in) :: e21
    REAL(wp), INTENT(in) :: e22
    REAL(wp), INTENT(in) :: alpha_evp

    REAL(wp), INTENT(inout) :: s12
    REAL(wp), INTENT(inout) :: s21
    REAL(wp), INTENT(out) :: s11
    REAL(wp), INTENT(out) :: s22
    REAL(wp), INTENT(inout) :: sigma_I
    REAL(wp), INTENT(inout) :: sigma_II
    REAL(wp), INTENT(in) :: zeta
    REAL(wp), INTENT(in) :: P
    REAL(wp)  :: B_evp,C_evp

    B_evp=(alpha_evp-1.0_wp)/alpha_evp! implizit alpha_evp/(1.0_wp+alpha_evp)
    C_evp=1.0_wp/alpha_evp !1.0_wp/(1.0_wp+alpha_evp)

    sigma_I= B_evp*sigma_I+C_evp*&
         &(2.0_wp*zeta*(e11+e22)-P)

    sigma_II=B_evp*sigma_II+C_evp*&
         &0.5_wp*zeta*(e11-e22)

    s12=B_evp*s12+C_evp*0.5_wp*zeta*e12
    s21=B_evp*s21+C_evp*0.5_wp*zeta*e21

    s11=0.5_wp*(sigma_I+sigma_II)
    s22=0.5_wp*(sigma_I-sigma_II)



  END SUBROUTINE mEVP



   SUBROUTINE interface_boundary_cell_marker(boundary_cell_marker, p_patch_3D, p_ice)
     TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
     TYPE(t_patch),  POINTER :: patch_2D
     TYPE(t_subset_range), POINTER :: all_cells
     TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
   !--------------------------------------------------------------------

    REAL(wp),TARGET,INTENT(inout) :: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

!    Real(wp):: x,y,z
    INTEGER :: cell_block,start_index,end_index,cell_index
    !-------------------------------------------------------------------
    patch_2D         => p_patch_3D%p_patch_2D(1)
    all_cells        =>patch_2D%cells%all

!    CALL sync_patch_array(SYNC_C, patch_2D, p_ice%conc )

    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index =  start_index, end_index
          boundary_cell_marker(cell_index,1,cell_block)=0.0_wp
          if(p_patch_3D%lsm_c(cell_index,1,cell_block) <= sea_boundary ) THEN  
             if ( p_ice%conc(cell_index,1,cell_block)>0.0_wp) then
                boundary_cell_marker(cell_index,1,cell_block)=1.0_wp
             endif
          endif
          
       END DO
    END DO

    CALL sync_patch_array(SYNC_C, patch_2D, boundary_cell_marker(:,1,:))


  END SUBROUTINE interface_boundary_cell_marker


  SUBROUTINE interface_boundary_edge_marker(boundary_edge_marker,boundary_cell_marker, p_patch_3D, p_ice)
    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D


    TYPE(t_subset_range), POINTER            :: edges_in_domain
    TYPE(t_patch),  POINTER                  :: patch_2D
    TYPE(t_sea_ice), INTENT(INOUT)           :: p_ice
    !-----------------------------------------------------------------------
    REAL(wp),TARGET,INTENT(in) :: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp),TARGET, INTENT(inout) :: boundary_edge_marker(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)

    INTEGER :: cell_block1,cell_block2,start_index,end_index,cell_index1,cell_index2, edge_index_i,&
         &edge_block_i,cell_index_1,cell_index_2, doy
    !-----------------------------------------------------------------------
    patch_2D         => p_patch_3D%p_patch_2D(1)
!    all_edges => patch_2D%edges%all
!    owned_edges => patch_2D%edges%owned
    edges_in_domain   => patch_2D%edges%in_domain


    DO edge_block_i = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, edge_block_i, start_index, end_index)
      DO edge_index_i =  start_index, end_index

        boundary_edge_marker(edge_index_i,edge_block_i)=0.0_wp

        IF (p_patch_3d%lsm_e(edge_index_i,1,edge_block_i) <= sea_boundary) THEN
          cell_index1 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,1)
          cell_block1 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,1)
          cell_index2 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,2)
          cell_block2 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,2)

          doy=boundary_cell_marker(cell_index1,1,cell_block1)+boundary_cell_marker(cell_index2,1,cell_block2)

          IF( (doy > 0) ) THEN
            boundary_edge_marker(edge_index_i,edge_block_i)=1.0_wp
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    CALL sync_patch_array(SYNC_E, patch_2D, boundary_edge_marker)


  END SUBROUTINE interface_boundary_edge_marker

  SUBROUTINE cell_area(cell_area_c,p_patch_3D)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    TYPE(t_patch),  POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: owned_cells

    REAL(wp),TARGET,INTENT(inout) :: cell_area_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    INTEGER :: cell_block,start_index,end_index,cell_index,&
         &edge_index_1,edge_block_1,edge_index_2,edge_block_2,edge_index_3,edge_block_3
    REAL(wp) :: s,e1,e2,e3,n1x,n2x,n3x,n1y,n2y,n3y,n1z,n2z,n3z

    patch_2D         => p_patch_3D%p_patch_2D(1)
!    all_cells => patch_2D%cells%all
    owned_cells => patch_2D%cells%owned

    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        edge_index_1 = patch_2D%cells%edge_idx(cell_index, cell_block, 1)
        edge_block_1 = patch_2D%cells%edge_blk(cell_index, cell_block, 1)

        edge_index_2 = patch_2D%cells%edge_idx(cell_index, cell_block, 2)
        edge_block_2 = patch_2D%cells%edge_blk(cell_index, cell_block, 2)

        edge_index_3 = patch_2D%cells%edge_idx(cell_index, cell_block, 3)
        edge_block_3 = patch_2D%cells%edge_blk(cell_index, cell_block, 3)

        e1=patch_2D%edges%primal_edge_length(edge_index_1,edge_block_1)
        e2=patch_2D%edges%primal_edge_length(edge_index_2,edge_block_2)
        e3=patch_2D%edges%primal_edge_length(edge_index_3,edge_block_3)


        s=(e1+e2+e3)*0.5_wp

        cell_area_c(cell_index,1,cell_block)=SQRT(s*(s-e1)*(s-e2)*(s-e3))
        !patch_2D%cells%area(cell_index,cell_block)

      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block

    CALL sync_patch_array(SYNC_C, patch_2D, cell_area_c)

  END SUBROUTINE cell_area


  SUBROUTINE init_mass_matrix(mass,cell_area_c,p_patch_3D)

     TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
     TYPE(t_patch),  POINTER :: patch_2D
     TYPE(t_subset_range), POINTER :: all_cells
     TYPE(t_subset_range), POINTER :: owned_cells

     REAL(wp),TARGET,INTENT(in) :: cell_area_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
     INTEGER :: cell_block,start_index,end_index,cell_index, edge_index_i,edge_block_i, neigbor
     REAL(wp), TARGET, INTENT(inout) :: mass(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)

     patch_2D         => p_patch_3D%p_patch_2D(1)
     all_cells => patch_2D%cells%all
     owned_cells => patch_2D%cells%owned

     CALL sync_patch_array(SYNC_E, patch_2D, mass)

     DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index

         DO neigbor=1,3!patch_2D%cells%num_edges(cell_index,cell_block)!no_primal_edges
           edge_index_i = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
           edge_block_i = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

           mass(edge_index_i,edge_block_i)=mass(edge_index_i,edge_block_i)&
                  &+cell_area_c(cell_index,1,cell_block)/3.0_wp
         !    write(0,*) 'mass', mass(edge_index_i,edge_block_i),cell_area_c(cell_index,1,cell_block)
         ENDDO !neigbor=1,patch_2D%num_edges i
        !  write(0,*) 'mass', mass(edge_index_i,edge_block_i),cell_area_c(cell_index,1,cell_block)

       ENDDO ! cell_index = start_index, end_index
     ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block

     CALL sync_patch_array(SYNC_E, patch_2D, mass)

  END SUBROUTINE init_mass_matrix




  SUBROUTINE compute_2Dvector_grad(cell_index,cell_block,uxx,uxy,uyx,uyy,cell_area,x1_c,x2_c,x3_c,p_patch_3D, p_ice)
    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_patch),  POINTER                  :: patch_2D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_subset_range), POINTER            :: all_edges
    TYPE(t_subset_range),     POINTER        :: all_cells
    TYPE(t_subset_range), POINTER            :: owned_cells

    INTEGER, INTENT(in) :: cell_index
    INTEGER, INTENT(in) :: cell_block
    REAL(wp), INTENT(in) :: cell_area
    REAL(wp), INTENT(in) :: x1_c
    REAL(wp), INTENT(in) :: x2_c
    REAL(wp), INTENT(in) :: x3_c
    REAL(wp), INTENT(out) :: uxx,uyy, uxy,uyx



    INTEGER :: edge_block_1,edge_index_1,&
         &edge_index_i, edge_block_i, neigbor,edge_index_j, edge_block_j, neigbor_j

    REAL(wp) :: ei,h_ei,Pi,Oi,l_sn,sc_pn,l_st,tix_1,tiy_1,tiz_1,tjx_1,tjy_1,tjz_1,correct_n,correct_t,&
         &nix_l,niy_l,tix_l, tiy_l,nix,niz,niy,tix,tiy,tiz,O1,nix_1,niy_1,niz_1,njx_1,njy_1,njz_1,&
         & U,V,Z,uxx_l,uxy_l,uyy_l,uyx_l,uzx_l,uzy_l, tix_k,tiy_k,tiz_k,nix_k,niy_k,niz_k, Ux,Vx,Uy,Vy,Uz,Vz 

    

    patch_2D   => p_patch_3D%p_patch_2D(1)

    all_cells => patch_2D%cells%all
    owned_cells => patch_2D%cells%owned
    all_edges => patch_2D%edges%all

    uxx_l=0.0_wp
    uxy_l=0.0_wp
    uyx_l=0.0_wp
    uyy_l=0.0_wp
    uzx_l=0.0_wp
    uzy_l=0.0_wp

    DO neigbor=1,3 !no_primal_edges                                                                                                                                                                               

       edge_index_i =patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
       edge_block_i =patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

       ei=patch_2D%edges%primal_edge_length(edge_index_i,edge_block_i)

       h_ei = 2.0_wp*cell_area/ei
       Pi=2.0_wp/h_ei

       Oi=patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)

       !       call  cell_local_N_T(edge_index_i,edge_block_i,neigbor,Oi,nix_l,niy_l,tix_l,tiy_l, x1_c,x2_c,x3_c,p_patch_3D)
       
       !Alt-----------------
       nix=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)
       niy=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)
       niz=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)
                   
       tix=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)
       tiy=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)
       tiz=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)

       !       U=p_ice%vn_e(edge_index_i,edge_block_i)*nix+&                                                                                 
       !            &p_ice%vt_e(edge_index_i,edge_block_i)*tix  
        
       !      V=p_ice%vn_e(edge_index_i,edge_block_i)*niy+&
       !           &p_ice%vt_e(edge_index_i,edge_block_i)*tiy
       
        !     Z=p_ice%vn_e(edge_index_i,edge_block_i)*niz+&
        !        &p_ice%vt_e(edge_index_i,edge_block_i)*tiz    
       nix=nix*Oi
       niy=niy*Oi
       niz=niz*Oi
       
       tix=tix*Oi
       tiy=tiy*Oi
       tiz=tiz*Oi
                   
       if(neigbor==1) then
          O1=Oi

          sc_pn=nix*x1_c + niy*x2_c + niz*x3_c
          nix_1=nix-sc_pn*x1_c
          niy_1=niy-sc_pn*x2_c
          niz_1=niz-sc_pn*x3_c
          l_sn=sqrt(nix_1*nix_1 + niy_1*niy_1 + niz_1*niz_1)
          nix_1=nix_1/l_sn
          niy_1=niy_1/l_sn
          niz_1=niz_1/l_sn
          
          sc_pn=tix*x1_c + tiy*x2_c + tiz*x3_c
          tix_1=tix-sc_pn*x1_c
          tiy_1=tiy-sc_pn*x2_c
          tiz_1=tiz-sc_pn*x3_c
          l_sn=sqrt(tix_1*tix_1 + tiy_1*tiy_1 + tiz_1*tiz_1)
          tix_1=tix_1/l_sn
          tiy_1=tiy_1/l_sn
          tiz_1=tiz_1/l_sn
       endif
       
       !Projection in Ebene     
       sc_pn=nix*x1_c + niy*x2_c + niz*x3_c
       nix=nix-sc_pn*x1_c
       niy=niy-sc_pn*x2_c
       niz=niz-sc_pn*x3_c
       l_sn=sqrt(nix*nix + niy*niy + niz*niz)
       nix=nix/l_sn
       niy=niy/l_sn
       niz=niz/l_sn
       
       sc_pn=tix*x1_c + tiy*x2_c + tiz*x3_c
       tix=tix-sc_pn*x1_c
       tiy=tiy-sc_pn*x2_c
       tiz=tiz-sc_pn*x3_c          
       l_sn=sqrt(tix*tix + tiy*tiy + tiz*tiz)
       tix=tix/l_sn
       tiy=tiy/l_sn
       tiz=tiz/l_sn

       !Projektion ende 
       nix_l =   tix_1*nix + tiy_1*niy + tiz_1*niz
       niy_l = - nix_1*nix - niy_1*niy - niz_1*niz
        
       tix_l =   tix_1*tix + tiy_1*tiy + tiz_1*tiz
       tiy_l = - nix_1*tix - niy_1*tiy - niz_1*tiz
       !ENDE ALT 
       !      !neu mit nix_l+ niy_l
       U=(p_ice%vn_e(edge_index_i,edge_block_i)*nix_l + p_ice%vt_e(edge_index_i,edge_block_i)*tix_l)*Oi          
       V=(p_ice%vn_e(edge_index_i,edge_block_i)*niy_l + p_ice%vt_e(edge_index_i,edge_block_i)*tiy_l)*Oi

       IF (ice_laplace_dynamics) THEN
       !Komponentenweise 2d Laplace
         uxx_l=uxx_l+1.0_wp*nix_l*Pi*U
         uxy_l=uxy_l+1.0_wp*niy_l*Pi*U
       
         uyx_l=uyx_l+1.0_wp*nix_l*Pi*V
         uyy_l=uyy_l+1.0_wp*niy_l*Pi*V
       ELSE
        
       ! strain rate tensor: 0.5*(nabal u +nabla u^T)
       uxx_l=uxx_l+1.0_wp*nix_l*Pi*U
       uxy_l=uxy_l+0.5_wp*niy_l*Pi*U+0.5_wp*nix_l*Pi*V
       uyx_l=uxy_l
       uyy_l=uyy_l+1.0_wp*niy_l*Pi*V

       ENDIF
       !       uzx_l=uzx_l+1.0_wp*nix_l*Pi*Z
       !       uzy_l=uzy_l+1.0_wp*niy_l*Pi*Z
       
    ENDDO !neigbor=1,patch_2D%num_edges i
    !%% Visulaisierung__________________________________________
    !trafo kartesisch Visu
    ! neu 
    ! ! Ux    Vx  0
    ! ! Uy    Vy  0
    ! ! Uz    Vz  0

    Ux=tix_1*uxx_l-nix_1*uyx_l
    Vx=tix_1*uxy_l-nix_1*uyy_l
    Uy=tiy_1*uxx_l-niy_1*uyx_l
    Vy=tiy_1*uxy_l-niy_1*uyy_l

    !Uz=tiz_1*uxx_l-niz_1*uyx_l
    !Vz=tiz_1*uxy_l-niz_1*uyy_l

    !uxx=uxx_l * tix_1 - uxy_l * nix_1
    !uyx=uyx_l * tiy_1 - uyy_l * niy_1       
    !neu:
    uxx=Ux * tix_1 - Vx * nix_1   
    ! uyx=Uy * tiy_1 - Vy * niy_1
    uxy=Ux * tiy_1 - Vx * niy_1
    !neu:
    !uzx=Uz * tiz_1 - Vz * niz_1    
    ! uzx=uzx_l * tiz_1 - uzy_l * niz_1
    
!    p_ice%hi(cell_index,1,cell_block)  = uxx!nix_1*sc_pn+tix_1*l_sn
!    p_ice%conc(cell_index,1,cell_block)= uxy!uyx_l * tiy_1 - uyy_l * niy_1 
    !return value 3d Laplace
    
    !     uxx=Ux!uxx_l
    !     uxy=Vx!uxy_l
    !     uyx=Uy!uyx_l
    !     uyy=Vy!uyy_l
    !     uzx=Uz!uzx_l
    !     uzy=Vz!uzy_l
    !________________________________________________visu
    !return value 2d Laplace
    uxx=uxx_l
    uxy=uxy_l
    uyx=uyx_l
    uyy=uyy_l
    !uzx=uzx_l                                                                                                                                                          
    !uzy=uzy_l

  
  END SUBROUTINE compute_2Dvector_grad




  SUBROUTINE compute_sigma(Au_n,Au_t,boundary_cell_marker,h_e,A_e,s_e,s11,s12,s21,s22,cell_area_c,x1_c,x2_c,x3_c,p_patch_3D,p_ice,p_os)
    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_patch),  POINTER                  :: patch_2D
    TYPE(t_sea_ice),          INTENT(IN)     :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    TYPE(t_subset_range), POINTER            :: all_edges
    TYPE(t_subset_range),     POINTER        :: all_cells
    TYPE(t_subset_range), POINTER            :: owned_cells

    REAL(wp),TARGET,INTENT(in) :: x1_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: x2_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: x3_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: s11(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: s12(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: s22(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),TARGET,INTENT(in) :: s21(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
 
    REAL(wp), TARGET, INTENT(in) ::  h_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  A_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) ::  s_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)

    REAL(wp), TARGET, INTENT(out) ::  Au_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(out) ::  Au_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)

    REAL(wp),TARGET,INTENT(in) :: cell_area_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    INTEGER :: edge_block_1,edge_index_1,cell_index,cell_block,start_index,end_index,&
         &edge_index_i, edge_block_i, neigbor,edge_index_j, edge_block_j, neigbor_j, cell_index1,cell_index2,cell_block1,cell_block2

    REAL(wp) :: ei,h_ei,Pi,Oi,l_sn,sc_pn,tix_1,tiy_1,tiz_1,tjx_1,tjy_1,tjz_1,&
         &nix_l,niy_l,nix,niz,niy,tix,tiy,tiz,tix_l, tiy_l,O1,nix_1,niy_1,niz_1,njx_1,njy_1,njz_1,&
         &cell_area,x1,x2,x3,ssh,marker,mass_ice


    Au_n(:,:)=0.0_wp
    Au_t(:,:)=0.0_wp

    patch_2D   => p_patch_3D%p_patch_2D(1)

    all_cells => patch_2D%cells%all
    owned_cells => patch_2D%cells%owned
    all_edges => patch_2D%edges%all

    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index

         IF ((boundary_cell_marker(cell_index,1,cell_block).GT.0.0_wp) .AND. &
              (p_patch_3d%lsm_c(cell_index,1,cell_block) <= sea_boundary)) THEN

         cell_area = cell_area_c(cell_index,1,cell_block)
         x1=x1_c(cell_index,1,cell_block)
         x2=x2_c(cell_index,1,cell_block)
         x3=x3_c(cell_index,1,cell_block)

         IF (ssh_in_icedyn_type == 1) THEN
           IF (vert_cor_type == 0) THEN
             ssh=p_os%p_prog(nold(1))%h(cell_index,cell_block)
           ELSEIF (vert_cor_type == 1) THEN
             ssh=p_os%p_prog(nold(1))%eta_c(cell_index,cell_block) + p_ice%draftave(cell_index,cell_block)
           ENDIF
         ELSEIF (ssh_in_icedyn_type ==0) THEN
           ssh = 0.0_wp
         ENDIF
         DO neigbor = 1,3 !no_primal_edges
           edge_index_i =patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
           edge_block_i =patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

           ei=patch_2D%edges%primal_edge_length(edge_index_i,edge_block_i)

           h_ei = 2.0_wp*cell_area/ei
           Pi=2.0_wp/h_ei   !fixme use a better name than pi

           Oi=patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)

           nix=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
           niy=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
           niz=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi

           tix=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
           tiy=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
           tiz=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi

           IF(neigbor==1) THEN
             O1=Oi
             sc_pn=nix*x1 + niy*x2 + niz*x3
             nix_1=nix-sc_pn*x1
             niy_1=niy-sc_pn*x2
             niz_1=niz-sc_pn*x3
             l_sn=SQRT(nix_1*nix_1 + niy_1*niy_1 + niz_1*niz_1)
             nix_1=nix_1/l_sn
             niy_1=niy_1/l_sn
             niz_1=niz_1/l_sn

             sc_pn=tix*x1 + tiy*x2 + tiz*x3
             tix_1=tix-sc_pn*x1
             tiy_1=tiy-sc_pn*x2
             tiz_1=tiz-sc_pn*x3
             l_sn=SQRT(tix_1*tix_1 + tiy_1*tiy_1 + tiz_1*tiz_1)
             tix_1=tix_1/l_sn
             tiy_1=tiy_1/l_sn
             tiz_1=tiz_1/l_sn
           ENDIF

           sc_pn=nix*x1 + niy*x2 + niz*x3
           nix=nix-sc_pn*x1
           niy=niy-sc_pn*x2
           niz=niz-sc_pn*x3
           l_sn=SQRT(nix*nix + niy*niy + niz*niz)
           nix=nix/l_sn
           niy=niy/l_sn
           niz=niz/l_sn

           sc_pn=tix*x1 + tiy*x2 + tiz*x3
           tix=tix-sc_pn*x1
           tiy=tiy-sc_pn*x2
           tiz=tiz-sc_pn*x3
           l_sn=SQRT(tix*tix + tiy*tiy + tiz*tiz)
           tix=tix/l_sn
           tiy=tiy/l_sn
           tiz=tiz/l_sn

           nix_l =   tix_1*nix + tiy_1*niy + tiz_1*niz
           niy_l = - nix_1*nix - niy_1*niy - niz_1*niz

           tix_l =   tix_1*tix + tiy_1*tiy + tiz_1*tiz
           tiy_l = - nix_1*tix - niy_1*tiy - niz_1*tiz
    !need to think about this if       
    !       IF (p_patch_3d%lsm_e(edge_index_i,1,edge_block_i) <= sea_boundary) THEN
           cell_index1 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,1)
           cell_block1 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,1)
           cell_index2 = patch_2D%edges%cell_idx(edge_index_i,edge_block_i,2)
           cell_block2 = patch_2D%edges%cell_blk(edge_index_i,edge_block_i,2)

           marker=boundary_cell_marker(cell_index1,1,cell_block1)+boundary_cell_marker(cell_index2,1,cell_block2)
           mass_ice=0.0_wp
           IF(marker>1)THEN
              mass_ice=max(rhoi*h_e(edge_index_i,edge_block_i)+s_e(edge_index_i,edge_block_i),9.0_wp)*9.81_wp
           ENDIF
           If(A_e(edge_index_i,edge_block_i)<0.01)THEN
              mass_ice=0.0_wp
           ENDIF
             Au_n(edge_index_i,edge_block_i)=Au_n(edge_index_i,edge_block_i)+&
                  &Pi*cell_area*Oi*(&
                  & s11(cell_index,1,cell_block)*nix_l*nix_l+s12(cell_index,1,cell_block)*nix_l*niy_l&
                  &+s21(cell_index,1,cell_block)*niy_l*nix_l+s22(cell_index,1,cell_block)*niy_l*niy_l)&
                  !this is the sea surface evelation
                  &+MERGE(0.0_wp,ssh*mass_ice*ei*Oi,ssh_in_icedyn_type == 2)

             Au_t(edge_index_i,edge_block_i)=Au_t(edge_index_i,edge_block_i)+&
                  &Pi*cell_area*Oi*(&
                  & s11(cell_index,1,cell_block)*tix_l*nix_l+s12(cell_index,1,cell_block)*tix_l*niy_l&
                  &+s21(cell_index,1,cell_block)*tiy_l*nix_l+s22(cell_index,1,cell_block)*tiy_l*niy_l)
     !      ENDIF
         ENDDO !neigbor=1,patch_2D%num_edges i
       ENDIF
     ENDDO
   ENDDO


  END SUBROUTINE compute_sigma


  SUBROUTINE Stabilization_sum(S_x,S_y,boundary_cell_marker,x1_c,x2_c,x3_c,p_ice,p_patch_3D)
   TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
   TYPE(t_sea_ice),          INTENT(IN)  :: p_ice
   TYPE(t_patch),  POINTER                  :: patch_2D
   
   REAL(wp), TARGET, INTENT(inout) ::  S_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
   REAL(wp), TARGET, INTENT(inout) ::  S_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
   REAL(wp),TARGET,INTENT(in) :: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2d(1)%alloc_cell_blocks)
   REAL(wp),TARGET, INTENT(in) :: x1_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) 
   REAL(wp),TARGET, INTENT(in) :: x2_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) 
   REAL(wp),TARGET, INTENT(in) :: x3_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
   
   REAL(wp) ::nx1,ny1,nx2,ny2,nx3,ny3,tx1,ty1,tx2,ty2,tx3,ty3
   
   INTEGER :: cell_block, start_index,  end_index, cell_index,&
        &edge_index_1,edge_index_2,edge_index_3,edge_block_1,edge_block_2,edge_block_3
   
   INTEGER :: edge_index_i, edge_block_i, neigbor
   
   REAL(wp) :: Oi,l_sn,sc_pn,tix_1,tiy_1,tiz_1,&
        &nix,niz,niy,tix,tiy,tiz,O1,O2,O3,nix_1,niy_1,niz_1,&
        &cell_area,x1,x2,x3,  Sx1,Sx2,Sx3,Sy1,Sy2,Sy3, a1,a2,a3,z1,z2,z3
   REAL(wp) :: v2_x,v2_y,v2_z,v3_x,v3_y, v3_z, p1,p2,p3,sc_p
   
   TYPE(t_subset_range), POINTER :: owned_cells
   TYPE(t_subset_range), POINTER :: all_edges
   TYPE(t_subset_range), POINTER :: all_cells
    
   patch_2D   => p_patch_3D%p_patch_2D(1)
   owned_cells =>patch_2D%cells%owned
   all_cells =>patch_2D%cells%all
   all_edges => patch_2D%edges%all
   
   DO cell_block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

         if(p_patch_3D%lsm_c(cell_index,1,cell_block) <= sea_boundary ) THEN    
         if(boundary_cell_marker(cell_index,1,cell_block)>0.0)then
         x1=x1_c(cell_index,1,cell_block)
         x2=x2_c(cell_index,1,cell_block)
         x3=x3_c(cell_index,1,cell_block)
         
         DO neigbor=1,3 !no_primal_edges
            edge_index_i =patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
            edge_block_i =patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

            Oi=patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)

            nix=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
            niy=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
            niz=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi
            
            tix=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
            tiy=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
            tiz=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi
            
            if(neigbor==1) then
               O1=Oi
               sc_pn=nix*x1 + niy*x2 + niz*x3
               nix_1=nix-sc_pn*x1
               niy_1=niy-sc_pn*x2
               niz_1=niz-sc_pn*x3
               l_sn=sqrt(nix_1*nix_1 + niy_1*niy_1 + niz_1*niz_1)
               nix_1=nix_1/l_sn
               niy_1=niy_1/l_sn
               niz_1=niz_1/l_sn
               
               sc_pn=tix*x1 + tiy*x2 + tiz*x3
               tix_1=tix-sc_pn*x1
               tiy_1=tiy-sc_pn*x2
               tiz_1=tiz-sc_pn*x3
               l_sn=sqrt(tix_1*tix_1 + tiy_1*tiy_1 + tiz_1*tiz_1)
               tix_1=tix_1/l_sn
               tiy_1=tiy_1/l_sn
               tiz_1=tiz_1/l_sn
            endif

            sc_pn=nix*x1 + niy*x2 + niz*x3
            nix=nix-sc_pn*x1
            niy=niy-sc_pn*x2
            niz=niz-sc_pn*x3
            l_sn=sqrt(nix*nix + niy*niy + niz*niz)
            nix=nix/l_sn
            niy=niy/l_sn
            niz=niz/l_sn
            
            sc_pn=tix*x1 + tiy*x2 + tiz*x3
            tix=tix-sc_pn*x1
            tiy=tiy-sc_pn*x2
            tiz=tiz-sc_pn*x3
            l_sn=sqrt(tix*tix + tiy*tiy + tiz*tiz)
            tix=tix/l_sn
            tiy=tiy/l_sn
            tiz=tiz/l_sn
            
            if(neigbor==1) then
               nx1 =(   tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
               ny1 =( - nix_1*nix - niy_1*niy - niz_1*niz)*Oi
               
               tx1 = (  tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
               ty1 = (- nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
            endif
            if(neigbor==2) then
               O2=Oi
               nx2 =  ( tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
               ny2 = (- nix_1*nix - niy_1*niy - niz_1*niz)*Oi
               
               tx2 =  ( tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
               ty2 = (- nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
            endif
            if(neigbor==3) then
               O3=Oi
               nx3 =  ( tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
               ny3 = (- nix_1*nix - niy_1*niy - niz_1*niz)*Oi
               
               tx3 =  ( tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
               ty3 = (- nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
            endif
         ENDDO !neigbor=1,patch_2D%num_edges i     
         
         edge_index_1 = patch_2D%cells%edge_idx(cell_index, cell_block, 1)
         edge_block_1 = patch_2D%cells%edge_blk(cell_index, cell_block, 1)

         edge_index_2 = patch_2D%cells%edge_idx(cell_index, cell_block, 2)
         edge_block_2 = patch_2D%cells%edge_blk(cell_index, cell_block, 2)

         edge_index_3 = patch_2D%cells%edge_idx(cell_index, cell_block, 3)
         edge_block_3 = patch_2D%cells%edge_blk(cell_index, cell_block, 3)

         ! nx1=patch_2D%edges%primal_normal(edge_index_1,edge_block_1)%v1
         ! ny1=patch_2D%edges%primal_normal(edge_index_1,edge_block_1)%v2
         ! tx1=patch_2D%edges%dual_normal(edge_index_1,edge_block_1)%v1
         ! ty1=patch_2D%edges%dual_normal(edge_index_1,edge_block_1)%v2
         
         !nx2=patch_2D%edges%primal_normal(edge_index_2,edge_block_2)%v1
         !ny2=patch_2D%edges%primal_normal(edge_index_2,edge_block_2)%v2
         !tx2=patch_2D%edges%dual_normal(edge_index_2,edge_block_2)%v1
         !ty2=patch_2D%edges%dual_normal(edge_index_2,edge_block_2)%v2
         
         
         ! nx3=patch_2D%edges%primal_normal(edge_index_3,edge_block_3)%v1
         ! ny3=patch_2D%edges%primal_normal(edge_index_3,edge_block_3)%v2
         ! tx3=patch_2D%edges%dual_normal(edge_index_3,edge_block_3)%v1
         ! ty3=patch_2D%edges%dual_normal(edge_index_3,edge_block_3)%v2


         ! [ u ] auf einer Zelle fuer alle 3 Kanten kartesisch (lokal)
         z1= patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(3)
         z2= patch_2D%edges%cartesian_center(edge_index_2,edge_block_2)%x(3)
         z3= patch_2D%edges%cartesian_center(edge_index_3,edge_block_3)%x(3)

         a1=1.0_wp
         a2=1.0_wp
         a3=1.0_wp

         v2_x = patch_2D%edges%cartesian_center(edge_index_2,edge_block_2)%x(1)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(1)
         v2_y = patch_2D%edges%cartesian_center(edge_index_2,edge_block_2)%x(2)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(2)
         v2_z = patch_2D%edges%cartesian_center(edge_index_2,edge_block_2)%x(3)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(3)

         v3_x = patch_2D%edges%cartesian_center(edge_index_3,edge_block_3)%x(1)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(1)
         v3_y = patch_2D%edges%cartesian_center(edge_index_3,edge_block_3)%x(2)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(2)
         v3_z = patch_2D%edges%cartesian_center(edge_index_3,edge_block_3)%x(3)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(3)

         p1=v2_y*v3_z-v2_z*v3_y
         p2=v2_z*v3_x-v2_x*v3_z
         p3=v2_x*v3_y-v2_y*v3_x

         sc_p= p1* x1_c(cell_index,1,cell_block)+p2*x2_c(cell_index,1,cell_block)+p3*x3_c(cell_index,1,cell_block)

         if(sc_p < 0.0) then
           a1=-1.0_wp
           a2=-1.0_wp
           a3=-1.0_wp
         endif 


         Sx1 = p_ice%vn_e(edge_index_2,edge_block_2)*(nx2)-p_ice%vn_e(edge_index_3,edge_block_3)*(nx3)+&
              &p_ice%vt_e(edge_index_2,edge_block_2)*(tx2)-p_ice%vt_e(edge_index_3,edge_block_3)*(tx3)
         
         Sx2 = p_ice%vn_e(edge_index_3,edge_block_3)*(nx3)-p_ice%vn_e(edge_index_1,edge_block_1)*(nx1)+&
              &p_ice%vt_e(edge_index_3,edge_block_3)*(tx3)-p_ice%vt_e(edge_index_1,edge_block_1)*(tx1)
         
         Sx3 = p_ice%vn_e(edge_index_1,edge_block_1)*(nx1)-p_ice%vn_e(edge_index_2,edge_block_2)*(nx2)+&
              &p_ice%vt_e(edge_index_1,edge_block_1)*(tx1)-p_ice%vt_e(edge_index_2,edge_block_2)*(tx2) 
         
         Sy1 = p_ice%vn_e(edge_index_2,edge_block_2)*(ny2)-p_ice%vn_e(edge_index_3,edge_block_3)*(ny3)+&
              &p_ice%vt_e(edge_index_2,edge_block_2)*(ty2)-p_ice%vt_e(edge_index_3,edge_block_3)*(ty3)
         
         Sy2 = p_ice%vn_e(edge_index_3,edge_block_3)*(ny3)-p_ice%vn_e(edge_index_1,edge_block_1)*(ny1)+&
              &p_ice%vt_e(edge_index_3,edge_block_3)*(ty3)-p_ice%vt_e(edge_index_1,edge_block_1)*(ty1)
         
         Sy3 = p_ice%vn_e(edge_index_1,edge_block_1)*(ny1)-p_ice%vn_e(edge_index_2,edge_block_2)*(ny2)+&
              &p_ice%vt_e(edge_index_1,edge_block_1)*(ty1)-p_ice%vt_e(edge_index_2,edge_block_2)*(ty2)                            

         ! normal- und tangentialkompoenente aufaddieren auf die Kanten, evtl. mit Oi multiplizieren?
         S_x(edge_index_1,edge_block_1)=S_x(edge_index_1,edge_block_1)+ (Sx1 * nx1 + Sy1 * ny1)*a1
         S_x(edge_index_2,edge_block_2)=S_x(edge_index_2,edge_block_2)+ (Sx2 * nx2 + Sy2 * ny2)*a2
         S_x(edge_index_3,edge_block_3)=S_x(edge_index_3,edge_block_3)+ (Sx3 * nx3 + Sy3 * ny3)*a3
         
         S_y(edge_index_1,edge_block_1)=S_y(edge_index_1,edge_block_1)+ (Sx1 * tx1 + Sy1 * ty1)*a1
         S_y(edge_index_2,edge_block_2)=S_y(edge_index_2,edge_block_2)+ (Sx2 * tx2 + Sy2 * ty2)*a2
         S_y(edge_index_3,edge_block_3)=S_y(edge_index_3,edge_block_3)+ (Sx3 * tx3 + Sy3 * ty3)*a3
         

         
         
!          S_x(edge_index_1,edge_block_1)=S_x(edge_index_1,edge_block_1)+&
!               &(p_ice%vn_e(edge_index_2,edge_block_2)*(nx2)-p_ice%vn_e(edge_index_3,edge_block_3)*(nx3)+&
!               & p_ice%vt_e(edge_index_2,edge_block_2)*(tx2)-p_ice%vt_e(edge_index_3,edge_block_3)*(tx3))
!          
!          S_x(edge_index_2,edge_block_2)=S_x(edge_index_2,edge_block_2)+&
!               &p_ice%vn_e(edge_index_3,edge_block_3)*(nx3)-p_ice%vn_e(edge_index_1,edge_block_1)*(nx1)+&
!               &p_ice%vt_e(edge_index_3,edge_block_3)*(tx3)-p_ice%vt_e(edge_index_1,edge_block_1)*(tx1)
!          
!          S_x(edge_index_3,edge_block_3)=S_x(edge_index_3,edge_block_3)+&
!               &p_ice%vn_e(edge_index_1,edge_block_1)*(nx1)-p_ice%vn_e(edge_index_2,edge_block_2)*(nx2)+&
!               &p_ice%vt_e(edge_index_1,edge_block_1)*(tx1)-p_ice%vt_e(edge_index_2,edge_block_2)*(tx2) 
!          
!          S_y(edge_index_1,edge_block_1)=S_y(edge_index_1,edge_block_1)+&
!               &p_ice%vn_e(edge_index_2,edge_block_2)*(ny2)-p_ice%vn_e(edge_index_3,edge_block_3)*(ny3)+&
!               &p_ice%vt_e(edge_index_2,edge_block_2)*(ty2)-p_ice%vt_e(edge_index_3,edge_block_3)*(ty3)
!          
!          S_y(edge_index_2,edge_block_2)=S_y(edge_index_2,edge_block_2)+&
!               &p_ice%vn_e(edge_index_3,edge_block_3)*(ny3)-p_ice%vn_e(edge_index_1,edge_block_1)*(ny1)+&
!               &p_ice%vt_e(edge_index_3,edge_block_3)*(ty3)-p_ice%vt_e(edge_index_1,edge_block_1)*(ty1)
!          
!          S_y(edge_index_3,edge_block_3)=S_y(edge_index_3,edge_block_3)+&
!               &p_ice%vn_e(edge_index_1,edge_block_1)*(ny1)-p_ice%vn_e(edge_index_2,edge_block_2)*(ny2)+&
!               &p_ice%vt_e(edge_index_1,edge_block_1)*(ty1)-p_ice%vt_e(edge_index_2,edge_block_2)*(ty2)                            

      endif
   endif
   ENDDO ! cell_index = start_index, end_index
ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block


 
  END SUBROUTINE Stabilization_sum


  SUBROUTINE Stabilization(S_x,S_y,S_n,S_t,boundary_cell_marker,x1_c,x2_c,x3_c,zeta_e,p_patch_3D)

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_patch),  POINTER                  :: patch_2D

    REAL(wp), TARGET, INTENT(out) :: S_n(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(out) :: S_t(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: S_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: S_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), TARGET, INTENT(in) :: boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), TARGET, INTENT(in) :: x1_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), TARGET, INTENT(in) :: x2_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), TARGET, INTENT(in) :: x3_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), TARGET, INTENT(in) :: zeta_e(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)




    REAL(wp) ::nx1,ny1,nx2,ny2,nx3,ny3,tx1,ty1,tx2,ty2,tx3,ty3,ei,alpha,O2,O3

    INTEGER :: cell_block, start_index,  end_index, cell_index,&
         &edge_index_1,edge_index_2,edge_index_3,edge_block_1,edge_block_2,edge_block_3
    INTEGER :: edge_index_i, edge_block_i, neigbor

    REAL(wp) :: Oi,l_sn,sc_pn,tix_1,tiy_1,tiz_1,&
         &nix,niz,niy,tix,tiy,tiz,O1,nix_1,niy_1,niz_1,&
         &cell_area,x1,x2,x3, Sx1,Sx2,Sx3,Sy1,Sy2,Sy3, a1,a2,a3,z1,z2,z3, zeta_1,zeta_2,zeta_3
    REAL(wp) :: v2_x, v2_y,v2_z,v3_x,v3_y,v3_z,p1,p2,p3,sc_p


    TYPE(t_subset_range), POINTER :: owned_cells
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: all_edges

    S_n(:,:)=0.0_wp
    S_t(:,:)=0.0_wp


    patch_2D   => p_patch_3D%p_patch_2D(1)
    owned_cells => patch_2D%cells%owned
    all_cells => patch_2D%cells%all
    all_edges => patch_2D%edges%all

    DO cell_block = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, cell_block, start_index, end_index)
       DO cell_index = start_index, end_index

         !IF ((boundary_cell_marker(cell_index,1,cell_block).GT.0.0_wp) .AND. &
         !     (p_patch_3d%lsm_c(cell_index,1,cell_block) <= sea_boundary)) THEN
          if(p_patch_3D%lsm_c(cell_index,1,cell_block) <= sea_boundary ) THEN
          if(boundary_cell_marker(cell_index,1,cell_block)>0.0) then 
          
           x1=x1_c(cell_index,1,cell_block)
           x2=x2_c(cell_index,1,cell_block)
           x3=x3_c(cell_index,1,cell_block)

           DO neigbor=1,3 !no_primal_edges

             edge_index_i =patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
             edge_block_i =patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

             Oi=patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)

             nix=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
             niy=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
             niz=patch_2D%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi

             tix=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)*Oi
             tiy=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)*Oi
             tiz=patch_2D%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)*Oi

             IF(neigbor==1) THEN
               O1=Oi
               sc_pn=nix*x1 + niy*x2 + niz*x3

               nix_1=nix-sc_pn*x1
               niy_1=niy-sc_pn*x2
               niz_1=niz-sc_pn*x3
               l_sn=SQRT(nix_1*nix_1 + niy_1*niy_1 + niz_1*niz_1)
               nix_1=nix_1/l_sn
               niy_1=niy_1/l_sn
               niz_1=niz_1/l_sn

               sc_pn=tix*x1 + tiy*x2 + tiz*x3
               tix_1=tix-sc_pn*x1
               tiy_1=tiy-sc_pn*x2
               tiz_1=tiz-sc_pn*x3
               l_sn=SQRT(tix_1*tix_1 + tiy_1*tiy_1 + tiz_1*tiz_1)
               tix_1=tix_1/l_sn
               tiy_1=tiy_1/l_sn
               tiz_1=tiz_1/l_sn
             ENDIF

             sc_pn=nix*x1 + niy*x2 + niz*x3
             nix=nix-sc_pn*x1
             niy=niy-sc_pn*x2
             niz=niz-sc_pn*x3
             l_sn=SQRT(nix*nix + niy*niy + niz*niz)
             nix=nix/l_sn
             niy=niy/l_sn
             niz=niz/l_sn

             sc_pn=tix*x1 + tiy*x2 + tiz*x3
             tix=tix-sc_pn*x1
             tiy=tiy-sc_pn*x2
             tiz=tiz-sc_pn*x3
             l_sn=SQRT(tix*tix + tiy*tiy + tiz*tiz)
             tix=tix/l_sn
             tiy=tiy/l_sn
             tiz=tiz/l_sn

             IF(neigbor==1) THEN
               nx1 = (   tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
               ny1 = ( - nix_1*nix - niy_1*niy - niz_1*niz)*Oi

               tx1 =   ( tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
               ty1 = ( - nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
             ENDIF
             IF(neigbor==2) THEN
               O2=Oi
               nx2 =   ( tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
               ny2 = ( - nix_1*nix - niy_1*niy - niz_1*niz)*Oi

               tx2 =   ( tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
               ty2 = ( - nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
             ENDIF
             IF(neigbor==3) THEN
               O3=Oi
               nx3 =   ( tix_1*nix + tiy_1*niy + tiz_1*niz)*Oi
               ny3 = ( - nix_1*nix - niy_1*niy - niz_1*niz)*Oi

               tx3 =   ( tix_1*tix + tiy_1*tiy + tiz_1*tiz)*Oi
               ty3 = ( - nix_1*tix - niy_1*tiy - niz_1*tiz)*Oi
             ENDIF
          ENDDO !neigbor=1,patch_2D%num_edges i
          !                     alpha=-zeta_c(cell_index, 1,cell_block)*1.e-13_wp/3.0_wp!*16.0_wp

          edge_index_1 = patch_2D%cells%edge_idx(cell_index, cell_block, 1)
          edge_block_1 = patch_2D%cells%edge_blk(cell_index, cell_block, 1)

          edge_index_2 = patch_2D%cells%edge_idx(cell_index, cell_block, 2)
          edge_block_2 = patch_2D%cells%edge_blk(cell_index, cell_block, 2)

          edge_index_3 = patch_2D%cells%edge_idx(cell_index, cell_block, 3)
          edge_block_3 = patch_2D%cells%edge_blk(cell_index, cell_block, 3)

          ei=patch_2D%edges%primal_edge_length(edge_index_1,edge_block_1)

          alpha=-1.0_wp/ei !*0.1


          !!! Sprung ueber die Kanten im lokalen kart. System
          Sx1 = S_x(edge_index_1,edge_block_1)*nx1 + S_y(edge_index_1,edge_block_1)*tx1
          Sx2 = S_x(edge_index_2,edge_block_2)*nx2 + S_y(edge_index_2,edge_block_2)*tx2
          Sx3 = S_x(edge_index_3,edge_block_3)*nx3 + S_y(edge_index_3,edge_block_3)*tx3

          Sy1 = S_x(edge_index_1,edge_block_1)*ny1 + S_y(edge_index_1,edge_block_1)*ty1
          Sy2 = S_x(edge_index_2,edge_block_2)*ny2 + S_y(edge_index_2,edge_block_2)*ty2
          Sy3 = S_x(edge_index_3,edge_block_3)*ny3 + S_y(edge_index_3,edge_block_3)*ty3

          !! Mit testfunktion
          z1= patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(3)
          z2= patch_2D%edges%cartesian_center(edge_index_2,edge_block_2)%x(3)
          z3= patch_2D%edges%cartesian_center(edge_index_3,edge_block_3)%x(3)

          !!
          zeta_1=zeta_e(edge_index_1,edge_block_1)
          zeta_2=zeta_e(edge_index_2,edge_block_2)
          zeta_3=zeta_e(edge_index_3,edge_block_3)

!          write(0,*) zeta_1,zeta_2,zeta_3

          a1=1.0_wp
          a2=1.0_wp
          a3=1.0_wp

          v2_x = patch_2D%edges%cartesian_center(edge_index_2,edge_block_2)%x(1)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(1)
          v2_y = patch_2D%edges%cartesian_center(edge_index_2,edge_block_2)%x(2)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(2)
          v2_z = patch_2D%edges%cartesian_center(edge_index_2,edge_block_2)%x(3)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(3)

          v3_x = patch_2D%edges%cartesian_center(edge_index_3,edge_block_3)%x(1)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(1)
          v3_y = patch_2D%edges%cartesian_center(edge_index_3,edge_block_3)%x(2)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(2)
          v3_z = patch_2D%edges%cartesian_center(edge_index_3,edge_block_3)%x(3)-patch_2D%edges%cartesian_center(edge_index_1,edge_block_1)%x(3)

          p1=v2_y*v3_z-v2_z*v3_y
          p2=v2_z*v3_x-v2_x*v3_z
          p3=v2_x*v3_y-v2_y*v3_x

          sc_p= p1* x1_c(cell_index,1,cell_block)+p2*x2_c(cell_index,1,cell_block)+p3*x3_c(cell_index,1,cell_block)


          if(sc_p < 0.0) then
            a1=-1.0_wp
            a2=-1.0_wp
            a3=-1.0_wp
          endif 




          S_n(edge_index_1,edge_block_1)=S_n(edge_index_1,edge_block_1)+ei/3.0_wp*alpha*&
               &a1*((Sx3*nx1+Sy3*ny1)*zeta_3 -(Sx2*nx1+Sy2*ny1)*zeta_2)

          S_t(edge_index_1,edge_block_1)=S_t(edge_index_1,edge_block_1)+ei/3.0_wp*alpha*&
               &a1*((Sx3*tx1+Sy3*ty1)*zeta_3 -(Sx2*tx1+Sy2*ty1)*zeta_2)

          S_n(edge_index_2,edge_block_2)=S_n(edge_index_2,edge_block_2)+ei/3.0_wp*alpha*&
               &a2*((Sx1*nx2+Sy1*ny2)*zeta_1 -(Sx3*nx2+Sy3*ny2)*zeta_3)

          S_t(edge_index_2,edge_block_2)=S_t(edge_index_2,edge_block_2)+ei/3.0_wp*alpha*&
               &a2*((Sx1*tx2+Sy1*ty2)*zeta_1 - (Sx3*tx2+Sy3*ty2)*zeta_3)

          S_n(edge_index_3,edge_block_3)=S_n(edge_index_3,edge_block_3)+ei/3.0_wp*alpha*&
               &a3*((Sx2*nx3+Sy2*ny3)*zeta_2 - (Sx1*nx3+Sy1*ny3)*zeta_1)

          S_t(edge_index_3,edge_block_3)=S_t(edge_index_3,edge_block_3)+ei/3.0_wp*alpha*&
               &a3*((Sx2*tx3+Sy2*ty3)*zeta_2 - (Sx1*tx3+Sy1*ty3)*zeta_1)


           !S_n(edge_index_1,1,edge_block_1)=S_n(edge_index_1,1,edge_block_1)+ei/3.0_wp*alpha*&
           !    &(zeta_e(edge_index_3,1,edge_block_3)*S_x(edge_index_3,1,edge_block_3)*nx1+&
           !    &zeta_e(edge_index_3,1,edge_block_3)*S_y(edge_index_3,1,edge_block_3)*ny1&
           !    &-zeta_e(edge_index_2,1,edge_block_2)*S_x(edge_index_2,1,edge_block_2)*nx1&
           !    &-zeta_e(edge_index_2,1,edge_block_2)*S_y(edge_index_2,1,edge_block_2)*ny1)

!             !alt aus 2d
!             S_n(edge_index_1,edge_block_1)=S_n(edge_index_1,edge_block_1)+ei/3.0_wp*alpha*&
!                  &(S_x(edge_index_3,edge_block_3)*nx1+S_y(edge_index_3,edge_block_3)*ny1&
!                  &-S_x(edge_index_2,edge_block_2)*nx1-S_y(edge_index_2,edge_block_2)*ny1)
!
!             S_t(edge_index_1,edge_block_1)=S_t(edge_index_1,edge_block_1)+ei/3.0_wp*alpha*&
!                  &(S_x(edge_index_3,edge_block_3)*tx1+S_y(edge_index_3,edge_block_3)*ty1&
!                  &-S_x(edge_index_2,edge_block_2)*tx1-S_y(edge_index_2,edge_block_2)*ty1)
!
!             S_n(edge_index_2,edge_block_2)=S_n(edge_index_2,edge_block_2)+ei/3.0_wp*alpha*&
!                  &(S_x(edge_index_1,edge_block_1)*nx2+S_y(edge_index_1,edge_block_1)*ny2&
!                  &-S_x(edge_index_3,edge_block_3)*nx2-S_y(edge_index_3,edge_block_3)*ny2)
!
!             S_t(edge_index_2,edge_block_2)=S_t(edge_index_2,edge_block_2)+ei/3.0_wp*alpha*&
!                  &(S_x(edge_index_1,edge_block_1)*tx2+S_y(edge_index_1,edge_block_1)*ty2&
!                  &-S_x(edge_index_3,edge_block_3)*tx2-S_y(edge_index_3,edge_block_3)*ty2)
!
!             S_n(edge_index_3,edge_block_3)=S_n(edge_index_3,edge_block_3)+ei/3.0_wp*alpha*&
!                  &(S_x(edge_index_2,edge_block_2)*nx3+S_y(edge_index_2,edge_block_2)*ny3&
!                  &-S_x(edge_index_1,edge_block_1)*nx3-S_y(edge_index_1,edge_block_1)*ny3)
!
!             S_t(edge_index_3,edge_block_3)=S_t(edge_index_3,edge_block_3)+ei/3.0_wp*alpha*&
!                  &(S_x(edge_index_2,edge_block_2)*tx3+S_y(edge_index_2,edge_block_2)*ty3&
!                  &-S_x(edge_index_1,edge_block_1)*tx3-S_y(edge_index_1,edge_block_1)*ty3)
       ENDIF
    ENDIF
 ENDDO ! cell_index = start_index, end_index

ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block

END SUBROUTINE Stabilization


  

END MODULE mo_ice_new_dynamics
