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

! Contains the implementation of velocity and tracer diffusion for the ICON ocean model.

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_velocity_diffusion
  
  USE mo_kind,                ONLY: wp
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_impl_constants,      ONLY: boundary, sea_boundary, min_dolic ! ,max_char_length
  USE mo_parallel_config,     ONLY: nproma
  USE mo_ocean_nml,           ONLY: n_zlev, iswm_oce, VelocityDiffusion_order, laplacian_form, &
    &  HarmonicVort_weight, HarmonicDiv_weight, BiharmonicVort_weight, BiharmonicDiv_weight

  USE mo_run_config,          ONLY: dtime
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, t_hydro_ocean_aux
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_ocean_physics_types, ONLY: t_ho_params
  USE mo_scalar_product,      ONLY: map_cell2edges_3d, map_edges2edges_viacell_3d_const_z
  USE mo_ocean_math_operators,ONLY: div_oce_3d, rot_vertex_ocean_3d,&
    & map_edges2vert_3d, grad_fd_norm_oce_3D, grad_vector, div_vector_onTriangle
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, sync_patch_array_mult
  USE mo_exception,           ONLY: finish !, message_text, message
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE
  
  PRIVATE
  
  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(LEN=12)           :: str_module    = 'oceDiffusion'  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  !
  ! PUBLIC INTERFACE
  !
  INTEGER, PARAMETER :: top=1
  PUBLIC :: velocity_diffusion
  PUBLIC :: velocity_diffusion_vertical_implicit
  PUBLIC :: velocity_diffusion_vertical_implicit_onBlock
  PUBLIC :: veloc_diff_harmonic_div_grad

  LOGICAL :: eliminate_upper_diag = .true.
  
CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates horizontal diffusion of edge velocity via laplacian diffusion
  !!    implemented as P^T div( K_H grad P v).
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE velocity_diffusion( patch_3D, vn_in, physics_parameters, p_diag,operators_coeff, &
    & laplacian_vn_out, lacc )
    
    TYPE(t_patch_3d ),TARGET :: patch_3D ! INTENT(in)
    REAL(wp)                 :: vn_in(:,:,:)! INTENT(in)
    TYPE(t_ho_params)        :: physics_parameters !mixing parameters INTENT(in)
    TYPE(t_hydro_ocean_diag) :: p_diag! INTENT(in)
    TYPE(t_operator_coeff),INTENT(IN)   :: operators_coeff! INTENT(in)
    REAL(wp)                 :: laplacian_vn_out(:,:,:)! INTENT(out)
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    
    !Local variables
    REAL(wp) :: z_lapl(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    !INTEGER  :: level
    LOGICAL :: lzacc
    CHARACTER(*), PARAMETER :: method_name = "velocity_diffusion"
    !-------------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')

    CALL set_acc_host_or_device(lzacc, lacc)

    IF(VelocityDiffusion_order==1)THEN
      
      !divgrad laplacian is chosen
      IF(laplacian_form==2)THEN
#ifdef _OPENACC
        IF (lzacc) CALL finish(method_name, 'OpenACC version currently for laplacian_form==2 not implemented')
#endif
        CALL finish(method_name, "form of harmonic Laplacian not recommended")
        CALL veloc_diff_harmonic_div_grad( patch_3D,      &
          & physics_parameters%HarmonicViscosity_coeff,   &
          & p_diag,    &
          & operators_coeff,&
          & laplacian_vn_out)
        
      ELSEIF(laplacian_form==1)THEN
        ! inUse
        CALL veloc_diff_harmonic_curl_curl(     &
          & patch_3D=patch_3D,                  &
          & u_vec_e=vn_in,                      &
          & vort=p_diag%vort,                   &
          & div_coeff=operators_coeff%div_coeff,&
          & HarmonicDiffusion=laplacian_vn_out, &
          & k_h=physics_parameters%HarmonicViscosity_coeff, &
          & lacc=lzacc)
       
      CALL dbg_print('laplacian_vn_out:', laplacian_vn_out,str_module,4, &
        & in_subset=patch_3D%p_patch_2D(1)%edges%owned)
  
      ENDIF
      
    ELSEIF(VelocityDiffusion_order==2 .or. VelocityDiffusion_order==21)THEN
      IF(laplacian_form==2)THEN
#ifdef _OPENACC
        IF (lzacc) CALL finish(method_name, 'OpenACC version currently for laplacian_form==2 not implemented')
#endif
        !CALL finish("mo_ocean_velocity_diffusion:velocity_diffusion", "form of biharmonic Laplacian not recommended")
        CALL veloc_diff_biharmonic_div_grad( patch_3D,   &
          & physics_parameters,      &
          & p_diag,       &
          & operators_coeff,   &
          & laplacian_vn_out)
          
          
      ELSEIF(laplacian_form==1)THEN
        
        CALL veloc_diff_biharmonic_curl_curl( &
          & patch_3D,            &
          & physics_parameters,  &
          & vn_in,               &
          & p_diag%vort,         &
          & operators_coeff,     &
          & laplacian_vn_out,    &
          & lacc=lzacc)
      ENDIF
    ELSEIF(VelocityDiffusion_order==0)THEN
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      laplacian_vn_out = 0.0_wp
      !$ACC END KERNELS
      !$ACC WAIT(1)
    ELSE
      CALL finish(method_name, "unknown VelocityDiffusion_order")
    ENDIF
    
  END SUBROUTINE velocity_diffusion
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates horizontal diffusion of edge velocity via laplacian diffusion
  !!    implemented as P^T div( K_H grad P v).
  !!
  !!
  SUBROUTINE veloc_diff_harmonic_div_grad( patch_3D, grad_coeff, p_diag,&
    & operators_coeff, laplacian_vn_out)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp)                               :: grad_coeff(:,:,:) ! grad_coeff contains the mixing parameters INTENT(in)
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(in) :: operators_coeff
    REAL(wp), INTENT(inout)           :: laplacian_vn_out(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    
    !Local variables
    INTEGER :: start_level, end_level
    INTEGER :: level, blockNo, edge_index,cell_index
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: start_index, end_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: idx_cartesian
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_cartesian_coordinates) :: z_grad_u    (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: z_div_grad_u(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells, all_edges, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-------------------------------------------------------------------------------
    patch_2D         => patch_3D%p_patch_2d(1)
    all_cells       => patch_2D%cells%ALL
    all_edges       => patch_2D%edges%ALL
    edges_in_domain => patch_2D%edges%in_domain
    !-------------------------------------------------------------------------------
    start_level = 1
    end_level = n_zlev
    

    !-------------------------------------------------------------------------------
    ! Note that HarmonicViscosity_coef is divided by dual_edge_length
    CALL grad_vector(cellVector=p_diag%p_vn, patch_3D=patch_3d, &
      grad_coeff=grad_coeff, gradVector=z_grad_u)

    CALL div_vector_onTriangle(patch_3d=patch_3D, edgeVector=z_grad_u, &
      & divVector=z_div_grad_u, div_coeff=operators_coeff%div_coeff)


!     laplacian_vn_out(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e) = 0.0_wp
!     
!     ! loop over cells in local domain + halo
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!       DO level = start_level, end_level
!         DO cell_index = start_index, end_index
!           z_div_grad_u(cell_index,level,blockNo)%x =  0.0_wp
!         END DO
!       END DO
!     END DO
!     
!     ! loop over edges in local domain + halo
!     DO blockNo = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
!       DO level = start_level, end_level
!         DO edge_index = start_edge_index, end_edge_index
!           z_grad_u(edge_index,level,blockNo)%x = 0.0_wp
!         ENDDO
!       END DO
!     END DO
!     
!     !-------------------------------------------------------------------------------------------------------
!     !Step 1: Calculate gradient of cell velocity vector.
!     !Result is a gradient vector, located at edges
!     !Step 2: Multiply each component of gradient vector with mixing coefficients
!     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!       CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
!       
!       DO level = start_level, end_level
!         DO edge_index = start_edge_index, end_edge_index
!           !IF ( v_base%lsm_e(edge_index,level,blockNo) <= sea_boundary ) THEN
!           IF (patch_3D%lsm_e(edge_index,level,blockNo) <= sea_boundary) THEN
!             !Get indices of two adjacent triangles
!             il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
!             ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
!             il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
!             ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)
!             
!             z_grad_u(edge_index,level,blockNo)%x = &
!               & physics_parameters%HarmonicViscosity_coeff(edge_index,level,blockNo)*  &
!               & ( p_diag%p_vn(il_c2,level,ib_c2)%x &
!               & - p_diag%p_vn(il_c1,level,ib_c1)%x)&
!               & / patch_2D%edges%dual_edge_length(edge_index,blockNo)
!             
!           ENDIF
!         ENDDO
!       END DO
!     END DO
!     DO idx_cartesian = 1,3
!       CALL sync_patch_array(sync_e, patch_2D,z_grad_u(:,:,:)%x(idx_cartesian) )
!     END DO
!     
!     !Step 2: Apply divergence to each component of mixing times gradient vector
!     iidx => patch_2D%cells%edge_idx
!     iblk => patch_2D%cells%edge_blk
!     
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!       
! #ifdef __SX__
! !CDIR UNROLL=6
! #endif
!       DO level = start_level, end_level
!         DO cell_index = start_index, end_index
!           
!           IF (patch_3D%lsm_c(cell_index,level,blockNo) >= boundary) THEN
!             z_div_grad_u(cell_index,level,blockNo)%x = 0.0_wp
!           ELSE
!             z_div_grad_u(cell_index,level,blockNo)%x =  &
!               & z_grad_u(iidx(cell_index,blockNo,1),level,iblk(cell_index,blockNo,1))%x&
!               & * operators_coeff%div_coeff(cell_index,level,blockNo,1)+&
!               & z_grad_u(iidx(cell_index,blockNo,2),level,iblk(cell_index,blockNo,2))%x&
!               & * operators_coeff%div_coeff(cell_index,level,blockNo,2)+&
!               & z_grad_u(iidx(cell_index,blockNo,3),level,iblk(cell_index,blockNo,3))%x&
!               & * operators_coeff%div_coeff(cell_index,level,blockNo,3)
!             
!           ENDIF
!         END DO
!       END DO
!     END DO
!     DO idx_cartesian = 1,3
!       CALL sync_patch_array(sync_c, patch_2D,z_div_grad_u(:,:,:)%x(idx_cartesian) )
!     END DO
    CALL sync_patch_array_mult(sync_c, patch_2D, 3, &
      z_div_grad_u(:,:,:)%x(1), z_div_grad_u(:,:,:)%x(2),z_div_grad_u(:,:,:)%x(3))
    
    !Step 3: Map divergence back to edges
    CALL map_cell2edges_3d( patch_3D, z_div_grad_u, laplacian_vn_out, operators_coeff)
    
  END SUBROUTINE veloc_diff_harmonic_div_grad
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates horizontal diffusion of edge velocity via bilaplacian diffusion
  !!    implemented as P^T (div K grad(divgrad P v)). Due to the position of the mixing matrix
  !!    which is following MPI-OM, the operator can not be written as simple iteration of
  !!    the laplacian in divgrad form.
  !!
  !!
  SUBROUTINE veloc_diff_biharmonic_div_grad0( patch_3D, physics_parameters, p_diag,&
    & operators_coeff, laplacian_vn_out)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    !REAL(wp), INTENT(in)              :: vn_in(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_ho_params), INTENT(in)     :: physics_parameters !mixing parameters
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(in) :: operators_coeff
    REAL(wp), INTENT(inout)           :: laplacian_vn_out(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    
    !Local variables
    INTEGER :: start_level, end_level
    INTEGER :: level, blockNo, edge_index,cell_index
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: start_index, end_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: idx_cartesian
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_cartesian_coordinates) :: z_grad_u          (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: z_div_grad_u      (nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, cells_oneEdgeInDomain
    TYPE(t_subset_range), POINTER :: all_edges, edges_in_domain, edges_gradIsCalculable
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_ocean_velocity_diffusion:velocity_diffusion_horz')
    !-------------------------------------------------------------------------------
    patch_2D               => patch_3D%p_patch_2d(1)
    all_cells              => patch_2D%cells%ALL
    cells_in_domain        => patch_2D%cells%in_domain
    cells_oneEdgeInDomain  => patch_2D%cells%one_edge_in_domain
    all_edges              => patch_2D%edges%ALL
    edges_in_domain        => patch_2D%edges%in_domain
    edges_gradIsCalculable => patch_2D%edges%gradIsCalculable
    !-------------------------------------------------------------------------------
    start_level = 1
    end_level = n_zlev
    
#ifdef NAGFOR
     z_div_grad_u(:,:,:)%x(1) = 0.0_wp
     z_div_grad_u(:,:,:)%x(2) = 0.0_wp
     z_div_grad_u(:,:,:)%x(3) = 0.0_wp
#endif
!     laplacian_vn_out  (1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e) = 0.0_wp
    
    ! loop over cells in local domain + halo
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!       DO level = start_level, end_level
!         DO cell_index = start_index, end_index
!           z_div_grad_u(cell_index,level,blockNo)%x =  0.0_wp
!         END DO
!       END DO
!     END DO
    
!ICON_OMP_PARALLEL PRIVATE(iidx, iblk )
    !-------------------------------------------------------------------------------------------------------
    !Step 1: Calculate gradient of cell velocity vector.
    !Result is a gradient vector, located at edges
!ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2 ) ICON_OMP_DEFAULT_SCHEDULE
    ! DO blockNo = edges_gradIsCalculable%start_block, edges_gradIsCalculable%end_block
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)

      DO edge_index = start_edge_index, end_edge_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)

          !Get indices of two adjacent triangles
          il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
          ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
          il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
          ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)

          z_grad_u(edge_index,level,blockNo)%x = & 
            & (p_diag%p_vn(il_c2,level,ib_c2)%x - p_diag%p_vn(il_c1,level,ib_c1)%x)&
            & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
        ENDDO
        ! zero the land levels
        DO level = patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)+1, end_level
          z_grad_u(edge_index,level,blockNo)%x = 0.0_wp
        ENDDO
      END DO
    END DO
!ICON_OMP_END_DO
    
!     DO idx_cartesian = 1,3
!       CALL sync_patch_array(sync_e, patch_2D,z_grad_u(:,:,:)%x(idx_cartesian) )
!     END DO
    
    !Step 2: Apply divergence to each component of gradient vector
    iidx => patch_2D%cells%edge_idx
    iblk => patch_2D%cells%edge_blk
    
!ICON_OMP_DO PRIVATE(start_index, end_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_oneEdgeInDomain%start_block, cells_oneEdgeInDomain%end_block
    ! DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(cells_oneEdgeInDomain, blockNo, start_index, end_index)
      DO cell_index = start_index, end_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(cell_index, blockNo)
        
          z_div_grad_u(cell_index,level,blockNo)%x =  &
            & z_grad_u(iidx(cell_index,blockNo,1),level,iblk(cell_index,blockNo,1))%x &
            & * operators_coeff%div_coeff(cell_index,level,blockNo,1) + &
            & z_grad_u(iidx(cell_index,blockNo,2),level,iblk(cell_index,blockNo,2))%x &
            & * operators_coeff%div_coeff(cell_index,level,blockNo,2) + &
            & z_grad_u(iidx(cell_index,blockNo,3),level,iblk(cell_index,blockNo,3))%x &
            & * operators_coeff%div_coeff(cell_index,level,blockNo,3)
            
        END DO

        ! this is only needed for the sync below, when running in parallel test mode !
        DO level = patch_3D%p_patch_1d(1)%dolic_c(cell_index, blockNo)+1, end_level
          z_div_grad_u(cell_index,level,blockNo)%x = 0.0_wp
        ENDDO
        
      END DO
    END DO
!ICON_OMP_END_DO
    
!     DO idx_cartesian = 1,3
!       CALL sync_patch_array(sync_c, patch_2D,z_div_grad_u(:,:,:)%x(idx_cartesian) )
!     END DO
    
    !Step 4: Repeat the application of div and grad and take the mixing coefficients into account
    !First the grad of previous result
    !now times the mixiing/friction coefficient
!ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2 ) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO edge_index = start_edge_index, end_edge_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
          
            !Get indices of two adjacent triangles
            il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
            ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
            il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
            ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)
            
            z_grad_u(edge_index,level,blockNo)%x =  &
              & ( z_div_grad_u(il_c2,level,ib_c2)%x      &
              &   - z_div_grad_u(il_c1,level,ib_c1)%x)   &
              & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
              
        ENDDO
      END DO
    END DO
!ICON_OMP_END_DO
!     DO idx_cartesian = 1,3
!       CALL sync_patch_array(sync_e, patch_2D,z_grad_u(:,:,:)%x(idx_cartesian) )
!     END DO
    
    !Step 5: Apply divergence to each component of gradient vector
    iidx => patch_2D%cells%edge_idx
    iblk => patch_2D%cells%edge_blk
    
!ICON_OMP_DO PRIVATE(start_index, end_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      DO cell_index = start_index, end_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(cell_index, blockNo)
          
          z_div_grad_u(cell_index,level,blockNo)%x =  &
            & -physics_parameters%BiharmonicViscosity_coeff(edge_index,level,blockNo) *( & ! take the negative div in order to avoid the negation of the laplacian
            & z_grad_u(iidx(cell_index,blockNo,1),level,iblk(cell_index,blockNo,1))%x &
            & * operators_coeff%div_coeff(cell_index,level,blockNo,1)+ &
            & z_grad_u(iidx(cell_index,blockNo,2),level,iblk(cell_index,blockNo,2))%x &
            & * operators_coeff%div_coeff(cell_index,level,blockNo,2)+ &
            & z_grad_u(iidx(cell_index,blockNo,3),level,iblk(cell_index,blockNo,3))%x &
            & * operators_coeff%div_coeff(cell_index,level,blockNo,3))
              
        END DO
      END DO
    END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL


    CALL sync_patch_array_mult(sync_c, patch_2D, 3,  &
     & z_div_grad_u(:,:,:)%x(1),  &
     & z_div_grad_u(:,:,:)%x(2),  &
     & z_div_grad_u(:,:,:)%x(3)   )
    
    
    !Step 6: Map divergence back to edges
    CALL map_cell2edges_3d( patch_3D, z_div_grad_u,laplacian_vn_out,operators_coeff)! requires cells_oneEdgeInDomain

!     ! this is not needed since we take the negative div
!     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!       CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
!       DO level = start_level, end_level
!         DO edge_index = start_edge_index, end_edge_index
!           
!           laplacian_vn_out(edge_index,level,blockNo)&
!           &=physics_parameters%k_veloc_h(edge_index,level,blockNo)*laplacian_vn_out(edge_index,level,blockNo)
!  
!         ENDDO
!       END DO
!     END DO
    
    
!        DO level=1,n_zlev
!         write(*,*)'Biharmonic divgrad',level,maxval(laplacian_vn_out(:,level,:)),&
!         &minval(laplacian_vn_out(:,level,:))
!        END DO
     CALL sync_patch_array(SYNC_E, patch_2D, laplacian_vn_out)
  END SUBROUTINE veloc_diff_biharmonic_div_grad0
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates horizontal diffusion of edge velocity via bilaplacian diffusion
  !!    implemented as P^T (div K grad(divgrad P v)). Due to the position of the mixing matrix
  !!    which is following MPI-OM, the operator can not be written as simple iteration of
  !!    the laplacian in divgrad form.
  !!
  !!
  SUBROUTINE veloc_diff_biharmonic_div_grad( patch_3D, physics_parameters, p_diag,&
    & operators_coeff, laplacian_vn_out)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    !REAL(wp), INTENT(in)              :: vn_in(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_ho_params), INTENT(in)     :: physics_parameters !mixing parameters
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(in) :: operators_coeff
    REAL(wp), INTENT(inout)           :: laplacian_vn_out(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    
    !Local variables
    INTEGER :: start_level, end_level
    INTEGER :: level, blockNo, edge_index,cell_index
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: start_index, end_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: idx_cartesian
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    REAL(wp):: z_grad_u_normal(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp):: z_grad_u_normal_ptp(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)   
    REAL(wp):: grad_div_e(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e) 
    REAL(wp):: div_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)        
    TYPE(t_cartesian_coordinates) :: z_grad_u          (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    !TYPE(t_cartesian_coordinates) :: z_div_grad_u      (nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, cells_oneEdgeInDomain
    TYPE(t_subset_range), POINTER :: all_edges, edges_in_domain, edges_gradIsCalculable
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_ocean_velocity_diffusion:velocity_diffusion_horz')
    !-------------------------------------------------------------------------------
    patch_2D               => patch_3D%p_patch_2d(1)
    all_cells              => patch_2D%cells%ALL
    cells_in_domain        => patch_2D%cells%in_domain
    cells_oneEdgeInDomain  => patch_2D%cells%one_edge_in_domain
    all_edges              => patch_2D%edges%ALL
    edges_in_domain        => patch_2D%edges%in_domain
    edges_gradIsCalculable => patch_2D%edges%gradIsCalculable
    !-------------------------------------------------------------------------------
    start_level = 1
    end_level = n_zlev
    
    z_grad_u_normal    (1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)          =0.0_wp
    z_grad_u_normal_ptp(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)          =0.0_wp    
    grad_div_e         (1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)          =0.0_wp
    div_c              (1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp 
    
!#ifdef NAGFOR
!     z_div_grad_u(:,:,:)%x(1) = 0.0_wp
!     z_div_grad_u(:,:,:)%x(2) = 0.0_wp
!     z_div_grad_u(:,:,:)%x(3) = 0.0_wp
!#endif
!     laplacian_vn_out  (1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e) = 0.0_wp
    
    
    !-------------------------------------------------------------------------------------------------------
    !Step 1: Calculate gradient of cell velocity vector.
    !Result is a gradient vector, located at edges
!ICON_OMP_DO PRIVATE(start_edge_index,end_edge_index,edge_index,level,il_c1,ib_c1,il_c2,ib_c2 ) ICON_OMP_DEFAULT_SCHEDULE
    ! DO blockNo = edges_gradIsCalculable%start_block, edges_gradIsCalculable%end_block
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)

      DO edge_index = start_edge_index, end_edge_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)

          !Get indices of two adjacent triangles
          il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
          ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
          il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
          ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)

          z_grad_u(edge_index,level,blockNo)%x = & !physics_parameters%k_veloc_h(edge_index,level,blockNo)*&
            & (p_diag%p_vn(il_c2,level,ib_c2)%x - p_diag%p_vn(il_c1,level,ib_c1)%x)&
            & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
            
            
            z_grad_u_normal(edge_index,level,blockNo)&
            &=DOT_PRODUCT(z_grad_u(edge_index,level,blockNo)%x, patch_2D%edges%primal_cart_normal(edge_index,blockNo)%x)
            
        ENDDO
        ! zero the land levels
        !DO level = patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)+1, end_level
        !  z_grad_u(edge_index,level,blockNo)%x = 0.0_wp
        !ENDDO
      END DO
    END DO
!ICON_OMP_END_DO
    CALL sync_patch_array(SYNC_E, patch_2D, z_grad_u_normal)
    
    !CALL map_edges2edges_viacell_3d_const_z( patch_3d, z_grad_u_normal, operators_coeff, &
    !    & z_grad_u_normal_ptp)   
    !CALL div_oce_3D( z_grad_u_normal_ptp, patch_3D, operators_coeff%div_coeff, div_c)
    CALL div_oce_3D( z_grad_u_normal, patch_3D, operators_coeff%div_coeff, div_c)
    CALL grad_fd_norm_oce_3D( div_c, patch_3D, operators_coeff%grad_coeff, grad_div_e)

! !     !Step 4: Repeat the application of div and grad and take the mixing coefficients into account
! !     !First the grad of previous result
! !     !now times the mixiing/friction coefficient
! ! !ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, edge_index, level, &
! ! !ICON_OMP il_c1, ib_c1, il_c2, ib_c2 ) ICON_OMP_DEFAULT_SCHEDULE
! !     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
! !       DO edge_index = start_edge_index, end_edge_index
! !         DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
! !           
! !             grad_div_e(edge_index,level,blockNo) &
! !             &= sqrt(physics_parameters%k_veloc_h(edge_index,level,blockNo)) * &
! !             & grad_div_e(edge_index,level,blockNo)
! !               
! !         ENDDO
! !       END DO
! !     END DO
! ! !ICON_OMP_END_DO
   CALL sync_patch_array(SYNC_E, patch_2D, grad_div_e)

    CALL div_oce_3D( grad_div_e, patch_3D, operators_coeff%div_coeff, div_c)

!ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, edge_index, level, il_c1, ib_c1, il_c2, ib_c2 ) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO edge_index = start_edge_index, end_edge_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
          
          il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
          ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
          il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
          ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)
          
          laplacian_vn_out(edge_index,level,blockNo) &
            &= -0.5_wp*physics_parameters%BiharmonicViscosity_coeff(edge_index,level,blockNo) &
            & * (div_c(il_c1,level,ib_c1)+div_c(il_c2,level,ib_c2))
              
        END DO
      END DO
    END DO
!ICON_OMP_END_DO
    !!Step 6: Map divergence back to edges
    !CALL map_cell2edges_3d( patch_3D, div_c,laplacian_vn_out,operators_coeff)! requires cells_oneEdgeInDomain
    
        DO level=1,n_zlev
         write(*,*)'Biharmonic divgrad',level,maxval(laplacian_vn_out(:,level,:)),&
         &minval(laplacian_vn_out(:,level,:))
        END DO
     CALL sync_patch_array(SYNC_E, patch_2D, laplacian_vn_out)
  END SUBROUTINE veloc_diff_biharmonic_div_grad
  !-------------------------------------------------------------------------
   
  
  !-------------------------------------------------------------------------
  !>
  !!  Computes  laplacian of a vector field in curl curl form.
  !!
  !! input:  lives on edges (velocity points)
  !! output: lives on edges
  !!
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
  SUBROUTINE veloc_diff_harmonic_curl_curl( patch_3D, u_vec_e, vort, div_coeff, &
    & nabla2_vec_e, HarmonicDiffusion, k_h, lacc )
    !
    TYPE(t_patch_3d ),TARGET      :: patch_3D  ! INTENT(in)
    REAL(wp)                      :: u_vec_e(:,:,:)  ! INTENT(in)
    REAL(wp)                      :: vort   (:,:,:) ! INTENT(in)
    REAL(wp)                      :: div_coeff(:,:,:,:) ! INTENT(in)
!     TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_dual    (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_v)
    REAL(wp),OPTIONAL             :: nabla2_vec_e(:,:,:) !  ! INTENT(out)
    REAL(wp),OPTIONAL             :: HarmonicDiffusion(:,:,:) !  ! INTENT(out)
    REAL(wp),OPTIONAL             :: k_h(:,:,:) ! INTENT(in)
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    !
    !Local variables
    INTEGER :: start_level, end_level     ! vertical start and end level
    INTEGER :: edge_index, level, blockNo
    INTEGER :: start_index, end_index
    LOGICAL :: lzacc
    REAL(wp) :: z_div_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)!, &
    REAL(wp) :: nabla2(nproma,n_zlev)
    !REAL(wp) ::  z_vn_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    !REAL(wp) ::  z_rot_v(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk, ividx, ivblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    ! note that this will go through the lateral boundaries
    patch_2D        => patch_3D%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-----------------------------------------------------------------------
    IF (present(HarmonicDiffusion) .and. .not. present(k_h)) THEN
      CALL finish('veloc_diff_harmonic_curl_curl','present(HarmonicDiffusion) .and. .not. present(k_h)')
    ENDIF
    IF (.not. present(HarmonicDiffusion) .and. present(k_h)) THEN
      CALL finish('veloc_diff_harmonic_curl_curl','.not. present(HarmonicDiffusion) .and. present(k_h)')
    ENDIF
    !-----------------------------------------------------------------------

    start_level = 1    
    icidx => patch_2D%edges%cell_idx
    icblk => patch_2D%edges%cell_blk
    ividx => patch_2D%edges%vertex_idx
    ivblk => patch_2D%edges%vertex_blk

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(z_div_c, nabla2) IF(lzacc)

    ! compute divergence of vector field
    ! z_div_c(:,:,patch_2D%alloc_cell_blocks) = 0.0_wp
    
#ifdef NAGFOR
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    z_div_c(:,:,:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)
#endif

    ! vn is synced on all edges
    CALL div_oce_3d( u_vec_e, patch_3D, div_coeff, z_div_c, subset_range=patch_2D%cells%all, lacc=lzacc)
!     CALL sync_patch_array(sync_c,patch_2D,z_div_c)
    
    ! compute rotation of vector field for the ocean
    !CALL rot_vertex_ocean_3D( patch_2D, u_vec_e, p_vn_dual, operators_coeff, z_rot_v)!
    !CALL sync_patch_array(SYNC_V,patch_2D,z_rot_v)
    !z_rot_v=vort
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, edge_index, level, nabla2) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      nabla2(:,:) = 0.0_wp
      !$ACC END KERNELS
      !$ACC WAIT(1)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO edge_index = start_index, end_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)

          !              write(0, *) "==============================="
          !              write(0, *) "0",  edge_index,level,blockNo
          !              write(0, *) "1",  p_patch_3D%wet_e(edge_index,level,blockNo)
          !              write(0,*)  "2",  k_h(edge_index,level,blockNo)
          !              write(0,*)  "3",  patch_2D%edges%tangent_orientation(edge_index,blockNo)
          !              write(0,*)  "4",  vort(ividx(edge_index,blockNo,2),level,ivblk(edge_index,blockNo,2))
          !              write(0,*)  "5",  vort(ividx(edge_index,blockNo,1),level,ivblk(edge_index,blockNo,1))
          !              write(0,*)  "6",  patch_2D%edges%inv_primal_edge_length(edge_index,blockNo)
          !              write(0,*)  "7",  z_div_c(icidx(edge_index,blockNo,2),level,icblk(edge_index,blockNo,2))
          !              write(0,*)  "8",  z_div_c(icidx(edge_index,blockNo,1),level,icblk(edge_index,blockNo,1))
          !              write(0,*)  "9",  patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
          !IF(v_base%lsm_e(edge_index,level,blockNo) < land_boundary)THEN

          nabla2(edge_index,level) = &   ! patch_3D%wet_e(edge_index,level,blockNo)* &
            & patch_2D%edges%tangent_orientation(edge_index,blockNo) *                      &
            & ( vort(ividx(edge_index,blockNo,2),level,ivblk(edge_index,blockNo,2))         &
            & - vort(ividx(edge_index,blockNo,1),level,ivblk(edge_index,blockNo,1)) )       &
            & * patch_2D%edges%inv_primal_edge_length(edge_index,blockNo) * HarmonicVort_weight   &
            & +                                                &
            & ( z_div_c(icidx(edge_index,blockNo,2),level,icblk(edge_index,blockNo,2))      &
            & - z_div_c(icidx(edge_index,blockNo,1),level,icblk(edge_index,blockNo,1)) )    &
            & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo) * HarmonicDiv_weight

        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      IF (present(nabla2_vec_e)) THEN
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        nabla2_vec_e(:,:,blockNo) = nabla2(:,:)
        !$ACC END KERNELS
        !$ACC WAIT(1)
      END IF
      IF (present(HarmonicDiffusion)) THEN
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        HarmonicDiffusion(:,:,blockNo) = nabla2(:,:) * k_h(:,:,blockNo)
        !$ACC END KERNELS
        !$ACC WAIT(1)
      END IF
    END DO
!ICON_OMP_END_PARALLEL_DO

   !$ACC END DATA
  END SUBROUTINE veloc_diff_harmonic_curl_curl
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!  Computes  Bilaplacian of a vector field.The placement of the mixing tensor
  !! follows the practice in NEMO. There is no corresponding operator in MPI-OM
  !!
  !! input:  lives on edges (velocity points)
  !! output: lives on edges
  !!
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
  SUBROUTINE veloc_diff_biharmonic_curl_curl( patch_3D,physics_parameters,u_vec_e,vort, operators_coeff,&
    & nabla4_vec_e, lacc)

    TYPE(t_patch_3d ),TARGET      :: patch_3D ! INTENT(in)
    TYPE(t_ho_params)             :: physics_parameters !mixing parameters  INTENT(in)
    REAL(wp)                      :: u_vec_e(:,:,:) ! INTENT(in)
    REAL(wp)                      :: vort   (:,:,:) ! INTENT(in)
    TYPE(t_operator_coeff),INTENT(IN) :: operators_coeff ! INTENT(in)
    REAL(wp)                      :: nabla4_vec_e(:,:,:) ! INTENT(out)
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    
    !Local variables
    REAL(wp), POINTER                      :: k_h(:,:,:)
    INTEGER :: start_level, end_level     ! vertical start and end level
    INTEGER :: edge_index, level, blockNo
    INTEGER :: start_index, end_index
    REAL(wp) ::  z_div_c   (nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) ::  z_rot_v   (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_v)
    REAL(wp) ::  z_nabla2_e(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
!     REAL(wp) ::  HarmonicDiffusion(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: h_e        (nproma,patch_3D%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates)  :: p_nabla2_dual(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_v)
#ifdef _OPENACC
    REAL(wp), DIMENSION(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v) ::p_nabla2_dual_x, p_nabla2_dual_y, p_nabla2_dual_z
#endif _OPENACC
    INTEGER,  DIMENSION(:,:,:), POINTER :: icidx, icblk, ividx, ivblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(z_div_c, z_rot_v, z_nabla2_e, p_nabla2_dual) &
    !$ACC   CREATE(p_nabla2_dual_x, p_nabla2_dual_y, p_nabla2_dual_z) IF(lzacc)

    !-----------------------------------------------------------------------
    ! note that this will go through the lateral boundaries
    patch_2D         => patch_3D%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    k_h => physics_parameters%BiharmonicViscosity_coeff
    !-----------------------------------------------------------------------
    
    start_level = 1
    end_level = n_zlev
    
    icidx => patch_2D%edges%cell_idx
    icblk => patch_2D%edges%cell_blk
    ividx => patch_2D%edges%vertex_idx
    ivblk => patch_2D%edges%vertex_blk
    
#ifdef NAGFOR
    ! this is only for sync with nag
    !$ACC KERNELS DEFAUL(PRESENT) ASYNC(1) IF(lzacc)
    z_nabla2_e(:,:,:) = 0.0_wp
    z_div_c(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks) =0.0_wp
    p_nabla2_dual(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v)%x(1)=0.0_wp
    p_nabla2_dual(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v)%x(2)=0.0_wp
    p_nabla2_dual(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v)%x(3)=0.0_wp
#ifdef _OPENACC
    p_nabla2_dual_x(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v)=0.0_wp
    p_nabla2_dual_y(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v)=0.0_wp
    p_nabla2_dual_z(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v)=0.0_wp
#endif _OPENACC
    !$ACC END KERNELS
    !$ACC WAIT(1)
#endif
!     z_rot_v(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v) =0.0_wp

!     IF (VelocityDiffusion_order==21) THEN
!       CALL veloc_diff_harmonic_curl_curl(     &
!         & patch_3D=patch_3D,                  &
!         & u_vec_e=u_vec_e,                    &
!         & vort=vort,                          &
!         & div_coeff=operators_coeff%div_coeff,&
!         & HarmonicDiffusion=HarmonicDiffusion,&
!         & nabla2_vec_e=z_nabla2_e,            &
!         & k_h=physics_parameters%HarmonicViscosity_coeff)
!     ELSE
      CALL veloc_diff_harmonic_curl_curl(     &
        & patch_3D=patch_3D,                  &
        & u_vec_e=u_vec_e,                    &
        & vort=vort,                          &
        & div_coeff=operators_coeff%div_coeff,&
        & nabla2_vec_e=z_nabla2_e,            &
        & lacc=lzacc)
!     ENDIF
  
    CALL sync_patch_array(sync_e,patch_2D,z_nabla2_e)
      
    ! compute divergence of vector field
    !     CALL div_oce_3d( u_vec_e, patch_2D, operators_coeff%div_coeff, z_div_c)
    !     ! DO level = start_level, end_level
    !     ! write(*,*)'vort1:',level,maxval(vort(:,level,:)),minval(vort(:,level,:))
    !     ! END DO
    !     ! compute rotation of vector field for the ocean
    !     !CALL rot_vertex_ocean_3D( patch_2D, u_vec_e, p_vn_dual, operators_coeff, z_rot_v)!
    !     !z_rot_v=vort
    !     !
    !     !  loop through all patch edges (and blocks)
    !     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !       CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
    !       DO edge_index = start_index, end_index
    !         DO level = start_level, end_level
    !           !DO edge_index = start_index, end_index
    !           !IF(v_base%lsm_e(edge_index,level,blockNo) < land_boundary)THEN
    !           z_nabla2_e(edge_index,level,blockNo) =  &
    !             & v_base%wet_e(edge_index,level,blockNo)*     &
    !             & (patch_2D%edges%tangent_orientation(edge_index,blockNo) *  &
    !             & ( vort(ividx(edge_index,blockNo,2),level,ivblk(edge_index,blockNo,2))  &
    !             & - vort(ividx(edge_index,blockNo,1),level,ivblk(edge_index,blockNo,1)) )  &
    !             & * patch_2D%edges%inv_primal_edge_length(edge_index,blockNo))  &
    !             & +v_base%wet_e(edge_index,level,blockNo)*&
    !             & (( z_div_c(icidx(edge_index,blockNo,2),level,icblk(edge_index,blockNo,2))    &
    !             & - z_div_c(icidx(edge_index,blockNo,1),level,icblk(edge_index,blockNo,1)) )  &
    !             & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo))
    !            z_nabla2_e(edge_index,level,blockNo)=sqrt(k_h(edge_index,level,blockNo))*z_nabla2_e(edge_index,level,blockNo)
    !         END DO
    !       END DO
    !     END DO
    
    ! compute divergence of vector field
    CALL div_oce_3d( z_nabla2_e, patch_3D, operators_coeff%div_coeff, z_div_c, &
      & subset_range=patch_2D%cells%all, lacc=lzacc)
!     CALL sync_patch_array(sync_c,patch_2D,z_div_c)
    
    ! compute rotation of vector field for the ocean
    CALL map_edges2vert_3d( patch_2D, &
      & z_nabla2_e,&
      & operators_coeff%edge2vert_coeff_cc,&
      & p_nabla2_dual, lacc=lzacc)

#ifdef _OPENACC
    ! FIXME 2024-09 DKRZ-dzo: The call to sync_patch_array fails on GPU when the entries
    !                         p_nabla2_dual(:,:,:)%x(1) to x(3) are passed directly.
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    p_nabla2_dual_x = p_nabla2_dual(:,:,:)%x(1)
    p_nabla2_dual_y = p_nabla2_dual(:,:,:)%x(2)
    p_nabla2_dual_z = p_nabla2_dual(:,:,:)%x(3)
    !$ACC END KERNELS
    !$ACC WAIT(1)

    CALL sync_patch_array_mult(sync_v, patch_2D, 3, &
      & p_nabla2_dual_x, p_nabla2_dual_y, p_nabla2_dual_z)

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    p_nabla2_dual(:,:,:)%x(1) = p_nabla2_dual_x
    p_nabla2_dual(:,:,:)%x(2) = p_nabla2_dual_y
    p_nabla2_dual(:,:,:)%x(3) = p_nabla2_dual_z
    !$ACC END KERNELS
    !$ACC WAIT(1)
#else
    CALL sync_patch_array_mult(sync_v, patch_2D, 3, &
      & p_nabla2_dual(:,:,:)%x(1), p_nabla2_dual(:,:,:)%x(2), p_nabla2_dual(:,:,:)%x(3))
#endif

    CALL rot_vertex_ocean_3d( patch_3D, z_nabla2_e, p_nabla2_dual, operators_coeff, z_rot_v, lacc=lzacc)
    
    !combine divergence and vorticity
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index,end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO edge_index = start_index, end_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
          
          nabla4_vec_e(edge_index,level,blockNo) =  &    ! patch_3D%wet_e(edge_index,level,blockNo) *&
            &  - k_h(edge_index,level,blockNo) * ( &
            & (patch_2D%edges%tangent_orientation(edge_index,blockNo) *  &
            & ( z_rot_v(ividx(edge_index,blockNo,2),level,ivblk(edge_index,blockNo,2))  &
            & - z_rot_v(ividx(edge_index,blockNo,1),level,ivblk(edge_index,blockNo,1)) )  &
            & * patch_2D%edges%inv_primal_edge_length(edge_index,blockNo)) * BiharmonicVort_weight  &
            & + & !patch_3D%wet_e(edge_index,level,blockNo)  *                    &
            & (( z_div_c(icidx(edge_index,blockNo,2),level,icblk(edge_index,blockNo,2))    &
            & - z_div_c(icidx(edge_index,blockNo,1),level,icblk(edge_index,blockNo,1)) )  &
            & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo) * BiharmonicDiv_weight ))

        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)
!ICON_OMP_END_DO

    IF (VelocityDiffusion_order==21) THEN
!ICON_OMP_DO PRIVATE(start_index,end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO edge_index = start_index, end_index
          DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)

          nabla4_vec_e(edge_index,level,blockNo) =      &
            & nabla4_vec_e(edge_index,level,blockNo) + &
            & z_nabla2_e(edge_index,level,blockNo) *    &
            & physics_parameters%HarmonicViscosity_coeff(edge_index,level,blockNo)

          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)
!ICON_OMP_END_DO NOWAIT
    ENDIF
!ICON_OMP_END_PARALLEL

    !$ACC END DATA
  END SUBROUTINE veloc_diff_biharmonic_curl_curl
  
  !-------------------------------------------------------------------------
  !
  !!Subroutine implements implicit vertical diffusion for horizontal velocity fields
  !!by inverting a scalar field..
  !>
  !!
  !------------------------------------------------------------------------
  SUBROUTINE velocity_diffusion_vertical_implicit( patch_3d,           &
    & velocity, a_v,                                                   &
    & operators_coefficients ) !,  &
    ! & diff_column)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)              :: velocity(:,:,:)   ! on edges
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)
    TYPE(t_operator_coeff),INTENT(IN) ,TARGET :: operators_coefficients
    !
    INTEGER :: start_index, end_index, edge_block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    
    !-----------------------------------------------------------------------
    edges_in_domain       =>  patch_3d%p_patch_2d(1)%edges%in_domain
    !-----------------------------------------------------------------------
   
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO edge_block = edges_in_domain%start_block, edges_in_domain%end_block    
      CALL get_index_range(edges_in_domain, edge_block, start_index, end_index)
      
      CALL velocity_diffusion_vertical_implicit_onBlock(  &
        & patch_3d,                                       &
        & velocity(:,:,edge_block),                       &
        & a_v(:,:,edge_block),                            &
        & operators_coefficients,                         &
        & start_index, end_index, edge_block)
        
    END DO
!ICON_OMP_END_PARALLEL_DO 

    eliminate_upper_diag = .not. eliminate_upper_diag ! switch the methods
    
  END SUBROUTINE velocity_diffusion_vertical_implicit
  !------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !
  !!Subroutine implements implicit vertical diffusion for horizontal velocity fields
  !!by inverting a scalar field..
  !>
  !! sbr identical to previous one, except for homogeneous boundary conditions
  !!
  !------------------------------------------------------------------------
!<Optimize:inUse>
#ifdef __LVECTOR__
  SUBROUTINE velocity_diffusion_vertical_implicit_onBlock( &
    & patch_3d,                            &
    & velocity,                            &
    & a_v,                                 &
    & operators_coefficients,              &
    & start_index, end_index, edge_block, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)              :: velocity(:,:)   ! on edges, (nproma, levels)
    REAL(wp), INTENT(inout)              :: a_v(:,:)      ! on edges, (nproma, levels)
    TYPE(t_operator_coeff),INTENT(IN) ,TARGET :: operators_coefficients
    INTEGER , INTENT(in):: start_index, end_index, edge_block
    LOGICAL, INTENT(in), OPTIONAL        :: lacc
    !
!     REAL(wp) :: dt_inv
    REAL(wp) :: inv_prism_thickness(nproma,1:n_zlev), inv_prisms_center_distance(nproma,1:n_zlev)
    REAL(wp) :: a(nproma,1:n_zlev), b(nproma,1:n_zlev), c(nproma,1:n_zlev)
    REAL(wp) :: column_velocity(nproma,1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)

    INTEGER :: bottom_level(nproma)
    INTEGER :: edge_index, level, max_end_level
    TYPE(t_patch), POINTER :: patch_2d
    CHARACTER(LEN=*), PARAMETER :: routine ='velocity_diffusion_vertical_implicit_onBlock (lvector)'
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    max_end_level = MAXVAL(patch_3d%p_patch_1d(1)%dolic_e(start_index:end_index,edge_block))

    !$ACC DATA CREATE(inv_prism_thickness, inv_prisms_center_distance, a, b, c) &
    !$ACC   CREATE(column_velocity, fact, bottom_level) IF(lzacc)

    !-----------------------------------------------------------------------
!     dt_inv = 1.0_wp/dtime
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO edge_index = start_index, end_index
      bottom_level(edge_index) = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
    ENDDO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
    ! Note : the inv_prism_thick_e, inv_prism_center_dist_e should be updated in calculate_thickness
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO level=1, max_end_level
      DO edge_index = start_index, end_index
        IF (bottom_level(edge_index) < 2 .OR. level > bottom_level(edge_index)) CYCLE ! nothing to diffuse
       
        inv_prism_thickness(edge_index,level)        = patch_3d%p_patch_1d(1)%inv_prism_thick_e(edge_index,level,edge_block)
        inv_prisms_center_distance(edge_index,level) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_e(edge_index,level,edge_block)
        
        column_velocity(edge_index,level) = velocity(edge_index,level)
               
      END DO ! edge_index = start_index, end_index
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    !------------------------------------
    ! Fill triangular matrix
    ! b is diagonal, a is the upper diagonal, c is the lower
    !   top level
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO edge_index = start_index, end_index
      a(edge_index,1) = 0.0_wp
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
    
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO edge_index = start_index, end_index
      IF (bottom_level(edge_index) < 2) CYCLE ! nothing to diffuse
            
      c(edge_index,1) = -a_v(edge_index,2) * & 
          & inv_prism_thickness(edge_index,1) * inv_prisms_center_distance(edge_index,2)*dtime
      b(edge_index,1) = 1.0_wp - c(edge_index,1)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
    
!     c(start_index:end_index,1) = -a_v(start_index:end_index,2) * & 
!       & inv_prism_thickness(start_index:end_index,1) * inv_prisms_center_distance(start_index:end_index,2)*dtime
!     b(start_index:end_index,1) = 1.0_wp - c(start_index:end_index,1)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO level = 2, max_end_level-1
      DO edge_index = start_index, end_index
!         bottom_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
        IF (bottom_level(edge_index) < 2 .OR. level > bottom_level(edge_index)-1) CYCLE ! nothing to diffuse
        a(edge_index,level) = - a_v(edge_index,level)   * inv_prism_thickness(edge_index,level) * inv_prisms_center_distance(edge_index,level)*dtime
        c(edge_index,level) = - a_v(edge_index,level+1) * inv_prism_thickness(edge_index,level) * inv_prisms_center_distance(edge_index,level+1)*dtime
        b(edge_index,level) = 1.0_wp - a(edge_index,level) - c(edge_index,level)
      END DO ! edge_index = start_index, end_index
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
      ! bottom
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO edge_index = start_index, end_index
!       bottom_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
      IF (bottom_level(edge_index) < 2) CYCLE
      a(edge_index,bottom_level(edge_index)) = -a_v(edge_index,bottom_level(edge_index)) *  &
        & inv_prism_thickness(edge_index,bottom_level(edge_index)) * &
        & inv_prisms_center_distance(edge_index,bottom_level(edge_index)) * dtime
      b(edge_index,bottom_level(edge_index)) = 1.0_wp - a(edge_index,bottom_level(edge_index))
      c(edge_index,bottom_level(edge_index)) = 0.0_wp
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    IF (eliminate_upper_diag) THEN
        ! solve the tridiagonal matrix by eliminating c (the upper diagonal) 
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO level = max_end_level-1, 1, -1
        !$ACC LOOP GANG VECTOR
        DO edge_index = start_index, end_index
!           bottom_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
          IF (bottom_level(edge_index) < 2 .OR. level > bottom_level(edge_index)-1) CYCLE ! nothing to diffuse
          fact(level)=c(edge_index,level)/b(edge_index,level+1)
          b(edge_index,level)=b(edge_index,level)-a(edge_index,level+1)*fact(level)
          c(edge_index,level) = 0.0_wp
          column_velocity(edge_index,level) = column_velocity(edge_index,level) - fact(level)*column_velocity(edge_index,level+1)
        END DO ! edge_index = start_index, end_index
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO edge_index = start_index, end_index
        IF (bottom_level(edge_index) < 2) CYCLE ! nothing to diffuse
        velocity(edge_index,1) = column_velocity(edge_index,1)/b(edge_index,1)
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO level = 2, max_end_level
        !$ACC LOOP GANG VECTOR
        DO edge_index = start_index, end_index
!           bottom_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
          IF (bottom_level(edge_index) < 2 .OR. level > bottom_level(edge_index)) CYCLE ! nothing to diffuse
          velocity(edge_index,level) = (column_velocity(edge_index,level) - &
            a(edge_index,level)*  velocity(edge_index,level-1)) / b(edge_index,level)    
        END DO ! edge_index = start_index, end_index
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)

    ELSE
      ! solve the tridiagonal matrix by eliminating a (the lower diagonal)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO level=2, max_end_level
        !$ACC LOOP GANG VECTOR
        DO edge_index = start_index, end_index
!           bottom_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
          IF (bottom_level(edge_index) < 2 .OR. level > bottom_level(edge_index)) CYCLE ! nothing to diffuse
          fact(level)=a(edge_index,level)/b(edge_index,level-1)
          b(edge_index,level)=b(edge_index,level)-c(edge_index,level-1)*fact(level)
          a(edge_index,level) = 0.0_wp
          column_velocity(edge_index,level) = column_velocity(edge_index,level) - fact(level)*column_velocity(edge_index,level-1)
        END DO ! edge_index = start_index, end_index
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)
      
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO edge_index = start_index, end_index
        IF (bottom_level(edge_index) < 2) CYCLE ! nothing to diffuse
        velocity(edge_index,bottom_level(edge_index)) = column_velocity(edge_index,bottom_level(edge_index))/ &
          & b(edge_index,bottom_level(edge_index))
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO level=max_end_level-1,1,-1
        !$ACC LOOP GANG VECTOR
        DO edge_index = start_index, end_index
!           bottom_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)
          IF (bottom_level(edge_index) < 2 .OR. level > bottom_level(edge_index)-1) CYCLE ! nothing to diffuse
          velocity(edge_index,level) = (column_velocity(edge_index,level) - &
            c(edge_index,level) * velocity(edge_index,level+1)) / b(edge_index,level)
            
        END DO ! edge_index = start_index, end_index
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)
      
    ENDIF

    !$ACC END DATA

  END SUBROUTINE velocity_diffusion_vertical_implicit_onBlock
#else
  SUBROUTINE velocity_diffusion_vertical_implicit_onBlock( &
    & patch_3d,                            &
    & velocity,                            &
    & a_v,                                 &
    & operators_coefficients,              &
    & start_index, end_index, edge_block, lacc)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)              :: velocity(:,:)   ! on edges, (nproma, levels)
    REAL(wp), INTENT(inout)              :: a_v(:,:)      ! on edges, (nproma, levels)
    TYPE(t_operator_coeff),INTENT(IN) ,TARGET :: operators_coefficients
    INTEGER , INTENT(in):: start_index, end_index, edge_block
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    !
!     REAL(wp) :: dt_inv
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
    REAL(wp) :: column_velocity(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)

    INTEGER :: bottom_level
    INTEGER :: edge_index, level
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
!     dt_inv = 1.0_wp/dtime

    !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(inv_prism_thickness, inv_prisms_center_distance) &
    !$ACC   PRIVATE(a, b, c, column_velocity, fact) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO edge_index = start_index, end_index
      bottom_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)

      IF (bottom_level < 2 ) CYCLE ! nothing to diffuse

      ! Note : the inv_prism_thick_e, inv_prism_center_dist_e should be updated in calculate_thickness
      DO level=1, bottom_level
        inv_prism_thickness(level)        = patch_3d%p_patch_1d(1)%inv_prism_thick_e(edge_index,level,edge_block)
        inv_prisms_center_distance(level) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_e(edge_index,level,edge_block)
        
        column_velocity(level) = velocity(edge_index,level)
      ENDDO

      !------------------------------------
      ! Fill triangular matrix
      ! b is diagonal, a is the upper diagonal, c is the lower
      !   top level
      a(1) = 0.0_wp
      c(1) = -a_v(edge_index,2) * inv_prism_thickness(1) * inv_prisms_center_distance(2)*dtime
      b(1) = 1.0_wp - c(1)
      DO level = 2, bottom_level-1
        a(level) = - a_v(edge_index,level)   * inv_prism_thickness(level) * inv_prisms_center_distance(level)*dtime
        c(level) = - a_v(edge_index,level+1) * inv_prism_thickness(level) * inv_prisms_center_distance(level+1)*dtime
        b(level) = 1.0_wp - a(level) - c(level)
      END DO
      ! bottom
      a(bottom_level) = -a_v(edge_index,bottom_level) *  &
        & inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level) * dtime
      b(bottom_level) = 1.0_wp - a(bottom_level)
      c(bottom_level) = 0.0_wp

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
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

  END SUBROUTINE velocity_diffusion_vertical_implicit_onBlock
#endif
  
END MODULE mo_ocean_velocity_diffusion
