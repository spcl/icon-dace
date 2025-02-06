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

! Contains the implementation of the limter of the ocean model

!----------------------------
#include "omp_definitions.inc"
#include "icon_definitions.inc"
!----------------------------
MODULE mo_ocean_limiter
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_math_types,                ONLY: t_cartesian_coordinates
  USE mo_math_constants,            ONLY: dbl_eps
  USE mo_impl_constants,            ONLY: sea_boundary, SEA
  USE mo_ocean_nml,                 ONLY: n_zlev
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma, p_test_run
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer
  USE mo_timer,                     ONLY: timer_start, timer_stop, timers_level, timer_adv_horz, timer_hflx_lim, &
    & timer_dif_horz, timer_extra10, timer_extra11, timer_extra12, timer_extra13, timer_extra15
  USE mo_ocean_tracer_transport_types,        ONLY: t_ocean_tracer
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish, message !, message_text, 
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_c1, sync_e, sync_patch_array, sync_patch_array_mult
  USE mo_mpi,                       ONLY: global_mpi_barrier
  USE mo_fortran_tools,             ONLY: set_acc_host_or_device
  
  IMPLICIT NONE
  
  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_ocean_limiter'
  
  CHARACTER(LEN=12)           :: str_module    = 'oceTracHorz '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  !
  ! PUBLIC INTERFACE
  !
  PUBLIC :: limiter_ocean_zalesak_horizontal, v_ppm_slimiter_mo_onBlock
    
  INTEGER, PARAMETER :: top=1
  
CONTAINS

  
  !-------------------------------------------------------------------------
  !>
  !! Flux limiter for horizontal advection
  !!
  !! Zalesak Flux-Limiter (Flux corrected transport)
  !! The corrected flux is a weighted average of the low order flux and the
  !! given high order flux. The high order flux is used to the greatest extent
  !! possible without introducing overshoots and undershoots.
  !! In vicinity of a lateral boundary only the low order flux is used: The criterion 
  !! is that at least one of the edges of the two neighboring cells of
  !! a central edges is a boundary edge.
  !! Note: This limiter is positive definite and almost monotone (but not strictly).
  !!
  !! Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !!
  !!  mpi note: computed on domain edges. Results is not synced.
  !!
!<Optimize:inUse>
 SUBROUTINE limiter_ocean_zalesak_horizontal( patch_3d,&
    & vert_velocity,          &
    & tracer,                 &
    & p_mass_flx_e,           &
    & flx_tracer_low,         &    
    & flx_tracer_high,        &
    & flx_tracer_final,       &
    & div_adv_flux_vert,      &   
    & operators_coefficients, &
    & h_old,                  &
    & h_new,                  &
    & lacc)                  
    
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp),INTENT(inout)              :: vert_velocity(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL       :: lacc

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (patch_3D%p_patch_2D(1)%cells%max_connectivity == 3) THEN
    
#ifdef __LVECTOR__
      CALL limiter_ocean_zalesak_horizontal_onTriangles_lvector( patch_3d,&
        & vert_velocity,          &
        & tracer,                 &
        & p_mass_flx_e,           &
        & flx_tracer_low,         &    
        & flx_tracer_high,        &
        & flx_tracer_final,       &
        & div_adv_flux_vert,      &   
        & operators_coefficients, &
        & h_old,                  &
        & h_new,                  &
        & lacc=lzacc)

#else     
      CALL limiter_ocean_zalesak_horizontal_onTriangles( patch_3d,&
        & vert_velocity,          &
        & tracer,                 &
        & p_mass_flx_e,           &
        & flx_tracer_low,         &    
        & flx_tracer_high,        &
        & flx_tracer_final,       &
        & div_adv_flux_vert,      &   
        & operators_coefficients, &
        & h_old,                  &
        & h_new,                  &
        & lacc=lzacc)
#endif
        
    ELSE
    
      CALL limiter_ocean_zalesak_horizontal_general( patch_3d,&
        & vert_velocity,          &
        & tracer,                 &
        & p_mass_flx_e,           &
        & flx_tracer_low,         &    
        & flx_tracer_high,        &
        & flx_tracer_final,       &
        & div_adv_flux_vert,      &   
        & operators_coefficients, &
        & h_old,                  &
        & h_new,                  &
        & lacc=lzacc)
        
    ENDIF
      
  END SUBROUTINE limiter_ocean_zalesak_horizontal
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  SUBROUTINE limiter_ocean_zalesak_horizontal_general( patch_3d,&
    & vert_velocity,          &
    & tracer,                 &
    & p_mass_flx_e,           &
    & flx_tracer_low,         &    
    & flx_tracer_high,        &
    & flx_tracer_final,       &
    & div_adv_flux_vert,      &   
    & operators_coefficients, &
    & h_old,                  &
    & h_new,                  &
    & lacc)                  
    
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp),INTENT(inout)              :: vert_velocity(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL       :: lacc
!     REAL(wp), INTENT(inout)             :: zlim(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    
    !Local variables
    !REAL(wp)              :: flx_tracer_high2  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp) :: z_mflx_anti(patch_3d%p_patch_2d(1)%cells%max_connectivity)
    REAL(wp) :: z_fluxdiv_c     !< flux divergence at cell center
    REAL(wp) :: z_anti          (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)    
    REAL(wp) :: z_tracer_new_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used
    REAL(wp) :: z_tracer_max    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local maximum of current tracer value and low order update
    REAL(wp) :: z_tracer_min    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local minimum of current tracer value and low order update
    REAL(wp) :: r_p             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< fraction which must multiply all in/out fluxes of cell jc to guarantee
    REAL(wp) :: r_m             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< no overshoot/undershoot
    REAL(wp) :: z_tracer_update_horz(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used    
    REAL(wp) :: r_frac          !< computed minimum fraction which must multiply< the flux at the edge
    REAL(wp) :: z_min, z_max    !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
    REAL(wp) :: p_p, p_m        !< sum of antidiffusive fluxes into and out of cell jc
    REAL(wp) :: prism_thick_old(n_zlev), inv_prism_thick_new(n_zlev)
    REAL(wp) :: delta_z_new, delta_z
    INTEGER, DIMENSION(:,:,:), POINTER ::  cellOfEdge_idx, cellOfEdge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    INTEGER :: start_level, end_level, edges_start_block, edges_end_block            
    INTEGER :: start_index, end_index, cells_start_block, cells_end_block 
    INTEGER :: edge_index, level, blockNo, jc,  cell_connect, sum_lsm_quad_edge, ctr
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc

    ! Pointers needed for GPU/OpenACC
    INTEGER, POINTER :: dolic_e(:,:)
    INTEGER, POINTER :: dolic_c(:,:)
    INTEGER, POINTER :: num_edges(:,:)
    INTEGER, POINTER :: edges_SeaBoundaryLevel(:,:,:)
    REAL(wp), POINTER :: div_coeff(:,:,:,:)
    REAL(wp), POINTER :: prism_thick_flat_sfc_c(:,:,:)
    REAL(wp), POINTER :: del_zlev_m(:)
    REAL(wp), POINTER :: inv_prism_thick_c(:,:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':limiter_ocean_zalesak_horizontal_general'
    !-------------------------------------------------------------------------

    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------
    start_level = 1
    end_level   = n_zlev
    edges_start_block = edges_in_domain%start_block
    edges_end_block   = edges_in_domain%end_block
    cells_start_block = cells_in_domain%start_block
    cells_end_block   = cells_in_domain%end_block

    cellOfEdge_idx  => patch_2d%edges%cell_idx
    cellOfEdge_blk  => patch_2d%edges%cell_blk
    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
    neighbor_cell_idx => patch_2d%cells%neighbor_idx
    neighbor_cell_blk => patch_2d%cells%neighbor_blk

    dolic_e => patch_3d%p_patch_1d(1)%dolic_e
    dolic_c => patch_3d%p_patch_1d(1)%dolic_c
    num_edges => patch_2d%cells%num_edges
    div_coeff => operators_coefficients%div_coeff
    prism_thick_flat_sfc_c => patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c
    del_zlev_m => patch_3d%p_patch_1d(1)%del_zlev_m
    inv_prism_thick_c => patch_3D%p_patch_1d(1)%inv_prism_thick_c
    edges_SeaBoundaryLevel => operators_coefficients%edges_SeaBoundaryLevel

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish(routine, 'OpenACC version currently not tested/validated')
#endif

    !$ACC DATA PRESENT(patch_3d%p_patch_2d(1)%nblks_e, patch_3d%p_patch_2d(1)%alloc_cell_blocks) &
    !$ACC   PRESENT(patch_3d%p_patch_2d(1)%cells%max_connectivity) IF(lzacc)

#ifdef NAGFOR
    z_tracer_max(:,:,:) = 0.0_wp
    z_tracer_min(:,:,:) = 0.0_wp
    r_m(:,:,:)          = 0.0_wp
    r_p(:,:,:)          = 0.0_wp
#endif
   
!ICON_OMP_PARALLEL
! !ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE        
!       DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block      
!         flux_div_vert(:,:,blockNo) = 0.0_wp     
!         CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)            
!         DO jc = start_index, end_index    
!           DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)        
!             ! positive vertical divergence in direction of w (upward positive)
!             flux_div_vert(jc,level,blockNo) = z_adv_flux_v(jc, level, blockNo) &
!             & - z_adv_flux_v(jc, level+1, blockNo)
!           ENDDO
!         END DO
!       END DO
! !ICON_OMP_END_DO      
! !       CALL sync_patch_array(sync_c, patch_2D, flux_div_vert)


!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE

    !$ACC DATA COPYIN(dolic_e, flx_tracer_high, flx_tracer_low) &
    !$ACC   COPY(z_anti) IF(lzacc)
    DO blockNo = edges_start_block, edges_end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      z_anti(:,:,blockNo)     = 0.0_wp
      !$ACC END KERNELS
      !$ACC WAIT(1)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO edge_index = start_index, end_index
        DO level = start_level, MIN(dolic_e(edge_index,blockNo), end_level)
          
          ! calculate antidiffusive flux for each edge
          z_anti(edge_index,level,blockNo) = flx_tracer_high(edge_index,level,blockNo)&
                                          &- flx_tracer_low(edge_index,level,blockNo)
        END DO  ! end loop over edges
      END DO  ! end loop over levels
      !$ACC END PARALLEL
      !$ACC WAIT(1)
    END DO  ! end loop over blocks
    !$ACC END DATA
!ICON_OMP_END_DO

    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, delta_z, delta_z_new, &
!ICON_OMP z_fluxdiv_c) ICON_OMP_DEFAULT_SCHEDULE
    !$ACC DATA COPYIN(div_adv_flux_vert, div_coeff, dolic_c, edge_of_cell_blk, edge_of_cell_idx) &
    !$ACC   COPYIN(flx_tracer_low, h_old, h_new, num_edges, prism_thick_flat_sfc_c, tracer) &
    !$ACC   COPY(z_tracer_max, z_tracer_min, z_tracer_new_low, z_tracer_update_horz) IF(lzacc)

    DO blockNo = cells_start_block, cells_end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      z_tracer_new_low(:,:,blockNo)    = 0.0_wp
      z_tracer_update_horz(:,:,blockNo)= 0.0_wp
      z_tracer_max(:,:,blockNo)        = 0.0_wp
      z_tracer_min(:,:,blockNo)        = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = start_index, end_index
        IF (dolic_c(jc,blockNo) < 1) CYCLE
        ! get prism thickness
        !inv_prism_thick_new(start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level) + h_new(jc,blockNo))
        !prism_thick_old(start_level)     = patch_3d%p_patch_1d(1)%del_zlev_m(start_level)           + h_old(jc,blockNo)
        !DO level = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
        !  prism_thick_old (level)    = patch_3d%p_patch_1d(1)%del_zlev_m(level)
        !  inv_prism_thick_new(level) = patch_3d%p_patch_1d(1)%inv_del_zlev_m(level)
        !ENDDO
        
        ! 3. Compute the complete (with horizontal and vertical divergence) updated low order solution z_tracer_new_low
        ! First at top level than in fluid interior
        !       
        level = start_level
        !  compute divergence of low order fluxes
        z_fluxdiv_c = 0
        !$ACC LOOP REDUCTION(+: z_fluxdiv_c)
        DO cell_connect = 1, num_edges(jc,blockNo)
          z_fluxdiv_c =  z_fluxdiv_c + &
            & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
            & div_coeff(jc,level,blockNo,cell_connect)
        ENDDO
       
        delta_z = prism_thick_flat_sfc_c(jc,level,blockNo)&
             &  + h_old(jc,blockNo)
        delta_z_new = prism_thick_flat_sfc_c(jc,level,blockNo)&
             &  + h_new(jc,blockNo)
             
        ! Low order flux at top level
        !z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
        !  & - dtime * (z_fluxdiv_c+flux_div_vert))/delta_z_new
        z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new

         z_tracer_update_horz(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
         & - dtime * (div_adv_flux_vert(jc,level,blockNo)))/delta_z_new

        !Fluid interior       
        DO level = start_level+1, MIN(dolic_c(jc,blockNo), end_level)       
          !  compute divergence of low order fluxes
          z_fluxdiv_c = 0
          !$ACC LOOP REDUCTION(+: z_fluxdiv_c)
          DO cell_connect = 1, num_edges(jc,blockNo)
            z_fluxdiv_c =  z_fluxdiv_c + &
              & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
              & div_coeff(jc,level,blockNo,cell_connect)
          ENDDO


          delta_z     = prism_thick_flat_sfc_c(jc,level,blockNo)
          delta_z_new = prism_thick_flat_sfc_c(jc,level,blockNo)
          !
          ! low order flux in flow interior
          !z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * prism_thick_old(level)     &
          !  & - dtime * (z_fluxdiv_c+flux_div_vert(jc,level,blockNo))) * inv_prism_thick_new(level)
          ! z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          !   & - dtime * (z_fluxdiv_c+flux_div_vert))/delta_z_new
          !z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          !    & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
              
          z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
            & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
            
          !z_tracer_update_horz(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          !  & - dtime * (div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
        ENDDO
      ENDDO
      !$ACC END PARALLEL
      
      ! precalculate local maximum/minimum of current tracer value and low order
      ! updated value
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      z_tracer_max(:,:,blockNo) =            &
        & MAX(          tracer(:,:,blockNo), &
        &     z_tracer_new_low(:,:,blockNo))
      z_tracer_min(:,:,blockNo) =            &
        & MIN(          tracer(:,:,blockNo), &
        &     z_tracer_new_low(:,:,blockNo))
      !$ACC END KERNELS

!      write(0,*) blockNo, ":", z_tracer_max(start_index:end_index,start_level:end_level,blockNo)
!      write(0,*) blockNo, ":", z_tracer_min(start_index:end_index,start_level:end_level,blockNo)
    ENDDO
    !$ACC WAIT(1)
    !$ACC END DATA
!ICON_OMP_END_DO

! DO level = 1, 4
! write(0,*)'tracer bounds',level, maxval(z_tracer_max(:,level,:)) ,minval(z_tracer_max(:,level,:)),&
! &maxval(z_tracer_min(:,level,:)) ,minval(z_tracer_min(:,level,:))
! END DO

!ICON_OMP_MASTER
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, z_tracer_max, z_tracer_min)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER
    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, inv_prism_thick_new, &
!ICON_OMP z_mflx_anti, z_max, z_min, cell_connect, p_p, p_m) ICON_OMP_DEFAULT_SCHEDULE
    !$ACC DATA COPYIN(del_zlev_m, div_coeff, dolic_c, edge_of_cell_blk, edge_of_cell_idx, h_new) &
    !$ACC   COPYIN(inv_prism_thick_c, neighbor_cell_blk, neighbor_cell_idx, num_edges) &
    !$ACC   COPYIN(z_anti, z_tracer_max, z_tracer_min, z_tracer_new_low) &
    !$ACC   COPY(inv_prism_thick_new, r_m, r_p, z_mflx_anti) IF(lzacc)

    DO blockNo = cells_start_block, cells_end_block

      ! this is only needed for the parallel test setups
      ! it will try  tocheck the uninitialized (land) parts
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      r_m(:,:,blockNo) = 0.0_wp
      r_p(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS
        
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = start_index, end_index
        
        ! get prism thickness
        inv_prism_thick_new(start_level) = 1.0_wp / (del_zlev_m(start_level)+ h_new(jc,blockNo))
        DO level = start_level+1, MIN(dolic_c(jc,blockNo), end_level)
          inv_prism_thick_new(level) = inv_prism_thick_c(jc,level,blockNo)
        ENDDO
        
        DO level = start_level, MIN(dolic_c(jc,blockNo), end_level)         
          ! 2. Define "antidiffusive" fluxes A(jc,level,blockNo,edge_index) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.          
          z_mflx_anti(:) = 0.0_wp
          z_max = z_tracer_max(jc,level,blockNo)
          z_min = z_tracer_min(jc,level,blockNo)
          p_p = 0.0_wp
          p_m = 0_wp
          !$ACC LOOP REDUCTION(+: p_m, p_p, z_max, z_min) PRIVATE(z_mflx_anti)
          DO cell_connect = 1, num_edges(jc,blockNo)
            IF (dolic_c(neighbor_cell_idx(jc,blockNo,cell_connect), neighbor_cell_blk(jc,blockNo,cell_connect)) >= level) THEN
              
              z_max = MAX(z_max, &
                & z_tracer_max(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))
              z_min = MIN(z_min, &
                & z_tracer_min(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))

              z_mflx_anti(cell_connect) =                                                        &
                & dtime * div_coeff(jc,level,blockNo,cell_connect) * inv_prism_thick_new(level)  &
                & * z_anti(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect))

              ! Sum of all incoming antidiffusive fluxes into cell jc
              ! outgoing fluxes carry a positive sign, incoming a negative
              p_p = p_p - MIN(0._wp, z_mflx_anti(cell_connect))
              ! Sum of all outgoing antidiffusive fluxes out of cell jc
              p_m = p_m + MAX(0._wp, z_mflx_anti(cell_connect))
            ENDIF
          ENDDO                
          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of tracer
          r_m(jc,level,blockNo) = (z_tracer_new_low(jc,level,blockNo) - z_min ) / (p_m + dbl_eps)!&
          !
          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of tracer
          r_p(jc,level,blockNo) = (z_max - z_tracer_new_low(jc,level,blockNo)) / (p_p + dbl_eps)!&
          !
          !update old tracer with low-order flux
          !!tracer(jc,level,blockNo)=z_tracer_new_low(jc,level,blockNo) 
          !tracer(jc,level,blockNo)=z_tracer_update_horz(jc,level,blockNo)         
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
    !$ACC WAIT(1)
    !$ACC END DATA
!ICON_OMP_END_DO

! DO level = 1, 10
! write(0,*)'updated:trac',level, maxval(tracer(:,level,:)) ,minval(tracer(:,level,:)),&
! &maxval(z_tracer_new_low(:,level,:)) ,minval(z_tracer_new_low(:,level,:))
! END DO

!ctr=0

    
!ICON_OMP_MASTER
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER

    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level, z_signum, r_frac) ICON_OMP_DEFAULT_SCHEDULE
    !$ACC DATA COPYIN(cellOfEdge_blk, cellOfEdge_idx, dolic_e, edges_SeaBoundaryLevel, flx_tracer_low, r_m, r_p) &
    !$ACC   COPYIN(z_anti) &
    !$ACC   COPY(flx_tracer_final) IF(lzacc)

    DO blockNo = edges_start_block, edges_end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      flx_tracer_final(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO edge_index = start_index, end_index
      
        DO level = start_level, MIN(dolic_e(edge_index,blockNo), end_level)
        
          IF(edges_SeaBoundaryLevel(edge_index,level,blockNo) > -2)THEN! edge < 2nd order boundary
          
            flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)
            
          ELSE!IF(sum_lsm_quad_edge==all_water_edges)THEN
          
            !z_anti>0 returns  1: here z_anti is outgoing, i.e. flux_high>flux_low
            !z_anti<0 returns -1: here z_anti is ingoing, i.e. flux_high<flux_low
            z_signum = SIGN(1._wp, z_anti(edge_index,level,blockNo))
                    
          ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
          ! but is computationally more efficient
          r_frac = 0.5_wp * (  &
            & (1._wp + z_signum) * & !<- active for z_signum=1
            & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)),  &
            &     r_p(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)))  &
            &+(1._wp - z_signum) * & !<- active for z_signum=-1
            & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)),  &
            &     r_p(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)))  )
          
          ! Limited flux
          flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)&
           & + MIN(1.0_wp,r_frac) *z_anti(edge_index,level,blockNo)      

          END IF
        END DO
      END DO
      !$ACC END PARALLEL
    ENDDO
    !$ACC WAIT(1)
    !$ACC END DATA
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
  !$ACC END DATA
  END SUBROUTINE limiter_ocean_zalesak_horizontal_general
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE limiter_ocean_zalesak_horizontal_onTriangles( patch_3d,&
    & vert_velocity,          &
    & tracer,                 &
    & p_mass_flx_e,           &
    & flx_tracer_low,         &    
    & flx_tracer_high,        &
    & flx_tracer_final,       &
    & div_adv_flux_vert,      &   
    & operators_coefficients, &
    & h_old,                  &
    & h_new,                  &
    & lacc)                  
    
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp),INTENT(inout)              :: vert_velocity(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL       :: lacc
!     REAL(wp), INTENT(inout)             :: zlim(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    
    !Local variables
    !REAL(wp)              :: flx_tracer_high2  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp) :: z_mflx_anti1, z_mflx_anti2, z_mflx_anti3
    REAL(wp) :: z_fluxdiv_c     !< flux divergence at cell center
    REAL(wp) :: z_anti          (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)    
    REAL(wp) :: z_tracer_new_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used
    REAL(wp) :: z_tracer_max    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local maximum of current tracer value and low order update
    REAL(wp) :: z_tracer_min    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local minimum of current tracer value and low order update
    REAL(wp) :: r_p             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< fraction which must multiply all in/out fluxes of cell jc to guarantee
    REAL(wp) :: r_m             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< no overshoot/undershoot
    REAL(wp) :: z_tracer_update_horz(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used    
    REAL(wp) :: r_frac          !< computed minimum fraction which must multiply< the flux at the edge
    REAL(wp) :: z_min, z_max    !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
    REAL(wp) :: p_p, p_m        !< sum of antidiffusive fluxes into and out of cell jc
    REAL(wp) :: inv_prism_thick_new
    REAL(wp) :: delta_z_new, delta_z
    INTEGER, DIMENSION(:,:,:), POINTER ::  cellOfEdge_idx, cellOfEdge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    INTEGER :: edge_blk1, edge_idx1, edge_blk2, edge_idx2, edge_blk3, edge_idx3
    INTEGER :: start_level !, end_level            
    INTEGER :: start_index, end_index
    INTEGER :: edges_start_block, edges_end_block, cells_start_block, cells_end_block
    INTEGER :: edge_index, level, blockNo, jc, sum_lsm_quad_edge, ctr
    INTEGER :: nidx1, nblk1, nidx2, nblk2, nidx3, nblk3
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc

    ! Pointers needed for GPU/OpenACC
    INTEGER, POINTER :: dolic_e(:,:)
    INTEGER, POINTER :: dolic_c(:,:)
    REAL(wp), POINTER :: div_coeff(:,:,:,:)
    REAL(wp), POINTER :: prism_thick_flat_sfc_c(:,:,:)
    REAL(wp), POINTER :: del_zlev_m(:)
    REAL(wp), POINTER :: inv_prism_thick_c(:,:,:)
    INTEGER, POINTER :: edges_SeaBoundaryLevel(:,:,:)
    !-------------------------------------------------------------------------
!     CALL message("","limiter_ocean_zalesak_horizontal_onTriangles is running...")

    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------
    start_level = 1
    edges_start_block = edges_in_domain%start_block
    edges_end_block   = edges_in_domain%end_block
    cells_start_block = cells_in_domain%start_block
    cells_end_block   = cells_in_domain%end_block

    cellOfEdge_idx  => patch_2d%edges%cell_idx
    cellOfEdge_blk  => patch_2d%edges%cell_blk
    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
    neighbor_cell_idx => patch_2d%cells%neighbor_idx
    neighbor_cell_blk => patch_2d%cells%neighbor_blk

    dolic_e => patch_3d%p_patch_1d(1)%dolic_e
    dolic_c => patch_3d%p_patch_1d(1)%dolic_c
    div_coeff => operators_coefficients%div_coeff
    prism_thick_flat_sfc_c => patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c
    del_zlev_m => patch_3d%p_patch_1d(1)%del_zlev_m
    inv_prism_thick_c => patch_3D%p_patch_1d(1)%inv_prism_thick_c
    edges_SeaBoundaryLevel => operators_coefficients%edges_SeaBoundaryLevel

    CALL set_acc_host_or_device(lzacc, lacc)

    ! All derived variables used to specify the size of OpenACC arrays in the declaration have to be present
    !$ACC DATA CREATE(r_m, r_p, z_anti, z_tracer_max, z_tracer_min, z_tracer_new_low, z_tracer_update_horz) IF(lzacc)

#ifdef NAGFOR
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    z_tracer_max(:,:,:) = 0.0_wp
    z_tracer_min(:,:,:) = 0.0_wp
    r_m(:,:,:)          = 0.0_wp
    r_p(:,:,:)          = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)
#endif
 
  IF (p_test_run) THEN
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    z_tracer_max(:,:,:) = 0.0_wp
    z_tracer_min(:,:,:) = 0.0_wp
    r_m(:,:,:)          = 0.0_wp
    r_p(:,:,:)          = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)
  ENDIF
 
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_start_block, edges_end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      z_anti(:,:,blockNo)     = 0.0_wp       
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO edge_index = start_index, end_index
        DO level = start_level, dolic_e(edge_index,blockNo)
          ! calculate antidiffusive flux for each edge
          z_anti(edge_index,level,blockNo) = flx_tracer_high(edge_index,level,blockNo)&
                                          &- flx_tracer_low(edge_index,level,blockNo)
        END DO  ! end loop over edges
      END DO  ! end loop over levels
      !$ACC END PARALLEL
    END DO  ! end loop over blocks
    !$ACC WAIT(1)
!ICON_OMP_END_DO

    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, delta_z, delta_z_new, &
!ICON_OMP z_fluxdiv_c, edge_blk1, edge_idx1, edge_blk2, edge_idx2, edge_blk3, edge_idx3) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_start_block, cells_end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      
!       z_tracer_new_low(:,:,blockNo)    = 0.0_wp
!       z_tracer_update_horz(:,:,blockNo)= 0.0_wp
!       z_tracer_max(:,:,blockNo)        = 0.0_wp
!       z_tracer_min(:,:,blockNo)        = 0.0_wp

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = start_index, end_index
        IF (dolic_c(jc,blockNo) < 1) CYCLE
        
        edge_blk1 = edge_of_cell_blk(jc,blockNo,1)
        edge_idx1 = edge_of_cell_idx(jc,blockNo,1)
        edge_blk2 = edge_of_cell_blk(jc,blockNo,2)
        edge_idx2 = edge_of_cell_idx(jc,blockNo,2)
        edge_blk3 = edge_of_cell_blk(jc,blockNo,3)
        edge_idx3 = edge_of_cell_idx(jc,blockNo,3)
        
        ! get prism thickness
        !inv_prism_thick_new(start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level) + h_new(jc,blockNo))
        !prism_thick_old(start_level)     = patch_3d%p_patch_1d(1)%del_zlev_m(start_level)           + h_old(jc,blockNo)
        !DO level = start_level+1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo), end_level)
        !  prism_thick_old (level)    = patch_3d%p_patch_1d(1)%del_zlev_m(level)
        !  inv_prism_thick_new(level) = patch_3d%p_patch_1d(1)%inv_del_zlev_m(level)
        !ENDDO
        
        ! 3. Compute the complete (with horizontal and vertical divergence) updated low order solution z_tracer_new_low
        ! First at top level than in fluid interior
        !       
        level = start_level
        !  compute divergence of low order fluxes
        z_fluxdiv_c = &
            & flx_tracer_low(edge_idx1,level,edge_blk1) * &
            & div_coeff(jc,level,blockNo,1) + & 
            & flx_tracer_low(edge_idx2,level,edge_blk2) * &
            & div_coeff(jc,level,blockNo,2) + & 
            & flx_tracer_low(edge_idx3,level,edge_blk3) * &
            & div_coeff(jc,level,blockNo,3) 
            
!         z_fluxdiv_c = 0
!         DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
!           z_fluxdiv_c =  z_fluxdiv_c + &
!             & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
!             & operators_coefficients%div_coeff(jc,level,blockNo,cell_connect)
!         ENDDO
       
        delta_z = prism_thick_flat_sfc_c(jc,level,blockNo)&
             &  + h_old(jc,blockNo)
        delta_z_new = prism_thick_flat_sfc_c(jc,level,blockNo)&
             &  + h_new(jc,blockNo)
             
        ! Low order flux at top level
        !z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
        !  & - dtime * (z_fluxdiv_c+flux_div_vert))/delta_z_new
        z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new

        z_tracer_update_horz(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
         & - dtime * (div_adv_flux_vert(jc,level,blockNo)))/delta_z_new

        z_tracer_max(jc,level,blockNo) =            &
          & MAX(          tracer(jc,level,blockNo), &
          &     z_tracer_new_low(jc,level,blockNo))
        z_tracer_min(jc,level,blockNo) =            &
          & MIN(          tracer(jc,level,blockNo), &
          &     z_tracer_new_low(jc,level,blockNo))
         
        !Fluid interior       
        DO level = start_level+1, dolic_c(jc,blockNo)       
          !  compute divergence of low order fluxes
          z_fluxdiv_c = &
            & flx_tracer_low(edge_idx1,level,edge_blk1) * &
            & div_coeff(jc,level,blockNo,1) + & 
            & flx_tracer_low(edge_idx2,level,edge_blk2) * &
            & div_coeff(jc,level,blockNo,2) + & 
            & flx_tracer_low(edge_idx3,level,edge_blk3) * &
            & div_coeff(jc,level,blockNo,3) 
         
!           z_fluxdiv_c = 0
!           DO cell_connect = 1, patch_2d%cells%num_edges(jc,blockNo)
!             z_fluxdiv_c =  z_fluxdiv_c + &
!               & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
!               & operators_coefficients%div_coeff(jc,level,blockNo,cell_connect)
!           ENDDO

          delta_z     = prism_thick_flat_sfc_c(jc,level,blockNo)
          delta_z_new = prism_thick_flat_sfc_c(jc,level,blockNo)
          !
          ! low order flux in flow interior
          !z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * prism_thick_old(level)     &
          !  & - dtime * (z_fluxdiv_c+flux_div_vert(jc,level,blockNo))) * inv_prism_thick_new(level)
          ! z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          !   & - dtime * (z_fluxdiv_c+flux_div_vert))/delta_z_new
          !z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          !    & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
              
          z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
            & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
            
          !z_tracer_update_horz(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
          !  & - dtime * (div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
          
          z_tracer_max(jc,level,blockNo) =            &
            & MAX(          tracer(jc,level,blockNo), &
            &     z_tracer_new_low(jc,level,blockNo))
          z_tracer_min(jc,level,blockNo) =            &
            & MIN(          tracer(jc,level,blockNo), &
            &     z_tracer_new_low(jc,level,blockNo))
         
        ENDDO
      ENDDO
      !$ACC END PARALLEL
      
      ! precalculate local maximum/minimum of current tracer value and low order
      ! updated value
!       z_tracer_max(:,:,blockNo) =            &
!         & MAX(          tracer(:,:,blockNo), &
!         &     z_tracer_new_low(:,:,blockNo))
!       z_tracer_min(:,:,blockNo) =            &
!         & MIN(          tracer(:,:,blockNo), &
!         &     z_tracer_new_low(:,:,blockNo))

!      write(0,*) blockNo, ":", z_tracer_max(start_index:end_index,start_level:end_level,blockNo)
!      write(0,*) blockNo, ":", z_tracer_min(start_index:end_index,start_level:end_level,blockNo)
    ENDDO
    !$ACC WAIT(1)
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

! DO level = 1, 4
! write(0,*)'tracer bounds',level, maxval(z_tracer_max(:,level,:)) ,minval(z_tracer_max(:,level,:)),&
! &maxval(z_tracer_min(:,level,:)) ,minval(z_tracer_min(:,level,:))
! END DO

! ! !ICON_OMP_MASTER
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, z_tracer_max, z_tracer_min)
! ! !ICON_OMP_END_MASTER
! ! !ICON_OMP_BARRIER

    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, inv_prism_thick_new, z_mflx_anti1, &
!ICON_OMP z_mflx_anti2, z_mflx_anti3, z_max, z_min, p_p, p_m, nidx1, nblk1, nidx2, nblk2, nidx3, nblk3, &
!ICON_OMP edge_blk1, edge_idx1, edge_blk2, edge_idx2, edge_blk3, edge_idx3) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_start_block, cells_end_block

      ! this is only needed for the parallel test setups
      ! it will try  tocheck the uninitialized (land) parts
!       r_m(:,:,blockNo) = 0.0_wp
!       r_p(:,:,blockNo) = 0.0_wp
        
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(edge_blk1, edge_idx1, edge_blk2, edge_idx2, edge_blk3, edge_idx3) &
      !$ACC   PRIVATE(inv_prism_thick_new, nblk1, nidx1, nblk2, nidx2, nblk3, nidx3) &
      !$ACC   PRIVATE(p_m, p_p, z_max, z_min, z_mflx_anti1, z_mflx_anti2, z_mflx_anti3)
      DO jc = start_index, end_index
            
        nidx1 = neighbor_cell_idx(jc,blockNo,1)
        nblk1 = neighbor_cell_blk(jc,blockNo,1)
        nidx2 = neighbor_cell_idx(jc,blockNo,2)
        nblk2 = neighbor_cell_blk(jc,blockNo,2)
        nidx3 = neighbor_cell_idx(jc,blockNo,3)
        nblk3 = neighbor_cell_blk(jc,blockNo,3)

        edge_blk1 = edge_of_cell_blk(jc,blockNo,1)
        edge_idx1 = edge_of_cell_idx(jc,blockNo,1)
        edge_blk2 = edge_of_cell_blk(jc,blockNo,2)
        edge_idx2 = edge_of_cell_idx(jc,blockNo,2)
        edge_blk3 = edge_of_cell_blk(jc,blockNo,3)
        edge_idx3 = edge_of_cell_idx(jc,blockNo,3)
        
        DO level = start_level, dolic_c(jc,blockNo)         
          ! 2. Define "antidiffusive" fluxes A(jc,level,blockNo,edge_index) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.
          
          
!           z_mflx_anti(:) = 0.0_wp

          IF (level == start_level) THEN
            inv_prism_thick_new = 1.0_wp / (del_zlev_m(start_level)+ h_new(jc,blockNo))
          ELSE
            inv_prism_thick_new = inv_prism_thick_c(jc,level,blockNo)
          ENDIF
          
          z_mflx_anti1 =                                                   &
            & dtime * div_coeff(jc,level,blockNo,1) * inv_prism_thick_new  &
            & * z_anti(edge_idx1,level,edge_blk1)
          z_mflx_anti2 =                                                   &
            & dtime * div_coeff(jc,level,blockNo,2) * inv_prism_thick_new  &
            & * z_anti(edge_idx2,level,edge_blk2)
          z_mflx_anti3 =                                                   &
            & dtime * div_coeff(jc,level,blockNo,3) * inv_prism_thick_new  &
            & * z_anti(edge_idx3,level,edge_blk3)
                  
          z_max = z_tracer_max(jc,level,blockNo)
          z_min = z_tracer_min(jc,level,blockNo)
          p_p = 0.0_wp
          p_m = 0_wp
          IF (dolic_c(nidx1, nblk1) >= level) THEN              
              z_max = MAX(z_max, &
                & z_tracer_max(nidx1,level,nblk1))
              z_min = MIN(z_min, &
                & z_tracer_min(nidx1,level,nblk1))

              ! Sum of all incoming antidiffusive fluxes into cell jc
              ! outgoing fluxes carry a positive sign, incoming a negative
              p_p = p_p - MIN(0._wp, z_mflx_anti1)
              ! Sum of all outgoing antidiffusive fluxes out of cell jc
              p_m = p_m + MAX(0._wp, z_mflx_anti1)
          ENDIF
          IF (dolic_c(nidx2, nblk2) >= level) THEN              
              z_max = MAX(z_max, &
                & z_tracer_max(nidx2,level,nblk2))
              z_min = MIN(z_min, &
                & z_tracer_min(nidx2,level,nblk2))

              ! Sum of all incoming antidiffusive fluxes into cell jc
              ! outgoing fluxes carry a positive sign, incoming a negative
              p_p = p_p - MIN(0._wp, z_mflx_anti2)
              ! Sum of all outgoing antidiffusive fluxes out of cell jc
              p_m = p_m + MAX(0._wp, z_mflx_anti2)
          ENDIF
          IF (dolic_c(nidx3, nblk3) >= level) THEN              
              z_max = MAX(z_max, &
                & z_tracer_max(nidx3,level,nblk3))
              z_min = MIN(z_min, &
                & z_tracer_min(nidx3,level,nblk3))

              ! Sum of all incoming antidiffusive fluxes into cell jc
              ! outgoing fluxes carry a positive sign, incoming a negative
              p_p = p_p - MIN(0._wp, z_mflx_anti3)
              ! Sum of all outgoing antidiffusive fluxes out of cell jc
              p_m = p_m + MAX(0._wp, z_mflx_anti3)
          ENDIF
                      
          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of tracer
          r_m(jc,level,blockNo) = (z_tracer_new_low(jc,level,blockNo) - z_min ) / (p_m + dbl_eps)!&
          !
          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of tracer
          r_p(jc,level,blockNo) = (z_max - z_tracer_new_low(jc,level,blockNo)) / (p_p + dbl_eps)!&
          !
          !update old tracer with low-order flux
          !!tracer(jc,level,blockNo)=z_tracer_new_low(jc,level,blockNo) 
          !tracer(jc,level,blockNo)=z_tracer_update_horz(jc,level,blockNo)         
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
    !$ACC WAIT(1)
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

! DO level = 1, 10
! write(0,*)'updated:trac',level, maxval(tracer(:,level,:)) ,minval(tracer(:,level,:)),&
! &maxval(z_tracer_new_low(:,level,:)) ,minval(z_tracer_new_low(:,level,:))
! END DO

!ctr=0

    
! ! !ICON_OMP_MASTER
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)
! ! !ICON_OMP_END_MASTER
! ! !ICON_OMP_BARRIER

    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
!ICON_OMP_DO_PARALLEL PRIVATE(start_index, end_index, edge_index, level, z_signum, r_frac) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_start_block, edges_end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      flx_tracer_final(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO edge_index = start_index, end_index
      
        DO level = start_level, dolic_e(edge_index,blockNo)
        
          IF( edges_SeaBoundaryLevel(edge_index,level,blockNo) > -2)THEN! edge < 2nd order boundary
          
            flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)
            
          ELSE!IF(sum_lsm_quad_edge==all_water_edges)THEN
          
            !z_anti>0 returns  1: here z_anti is outgoing, i.e. flux_high>flux_low
            !z_anti<0 returns -1: here z_anti is ingoing, i.e. flux_high<flux_low
            z_signum = SIGN(1._wp, z_anti(edge_index,level,blockNo))
                    
            ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
            ! but is computationally more efficient
            r_frac = 0.5_wp * (       &
              & (1._wp + z_signum) * & !<- active for z_signum=1
              & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)),  &
              &     r_p(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)))  &
              &+(1._wp - z_signum) * & !<- active for z_signum=-1
              & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)),  &
              &     r_p(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)))  )
            
            ! Limited flux
            flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)&
            & + MIN(1.0_wp,r_frac) *z_anti(edge_index,level,blockNo)      

          ENDIF
        END DO
      END DO
      !$ACC END PARALLEL
    ENDDO
    !$ACC WAIT(1)
!ICON_OMP_END_DO_PARALLEL
! !ICON_OMP_END_PARALLEL

  !$ACC END DATA
  END SUBROUTINE limiter_ocean_zalesak_horizontal_onTriangles
  !-------------------------------------------------------------------------
   
#ifdef __LVECTOR__
  !-------------------------------------------------------------------------
  ! #ifdef __LVECTOR__ is true to invoke this routine
  ! vectorized for nec
  !-------------------------------------------------------------------------
  SUBROUTINE limiter_ocean_zalesak_horizontal_onTriangles_lvector( patch_3d,&
    & vert_velocity,          &
    & tracer,                 &
    & p_mass_flx_e,           &
    & flx_tracer_low,         &    
    & flx_tracer_high,        &
    & flx_tracer_final,       &
    & div_adv_flux_vert,      &   
    & operators_coefficients, &
    & h_old,                  &
    & h_new,                  &
    & lacc)                  
    
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp),INTENT(inout)              :: vert_velocity(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e) 
    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), INTENT(inout)             :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(in)                :: h_old(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                :: h_new(1:nproma,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL :: lacc
!     REAL(wp), INTENT(inout)             :: zlim(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)

    
    !Local variables
    !REAL(wp)              :: flx_tracer_high2  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)     
    REAL(wp), DIMENSION(nproma) :: z_mflx_anti1, z_mflx_anti2, z_mflx_anti3
    REAL(wp) :: z_fluxdiv_c(nproma)     !< flux divergence at cell center
    REAL(wp) :: z_anti          (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)    
    REAL(wp) :: z_tracer_new_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used
    REAL(wp) :: z_tracer_max    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local maximum of current tracer value and low order update
    REAL(wp) :: z_tracer_min    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local minimum of current tracer value and low order update
    REAL(wp) :: r_p             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< fraction which must multiply all in/out fluxes of cell jc to guarantee
    REAL(wp) :: r_m             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< no overshoot/undershoot
    REAL(wp) :: z_tracer_update_horz(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used    
    REAL(wp) :: r_frac          !< computed minimum fraction which must multiply< the flux at the edge
    REAL(wp), DIMENSION(nproma) :: z_min, z_max    !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
    REAL(wp), DIMENSION(nproma) :: p_p, p_m        !< sum of antidiffusive fluxes into and out of cell jc
    REAL(wp) :: inv_prism_thick_new(nproma, n_zlev)
    REAL(wp) :: delta_z_new(nproma), delta_z(nproma)
    INTEGER, DIMENSION(:,:,:), POINTER ::  cellOfEdge_idx, cellOfEdge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edges_SeaBoundaryLevel
    REAL(wp), DIMENSION(:,:,:,:), POINTER :: div_coeff
    INTEGER :: start_level !, end_level            
    INTEGER :: start_index, end_index, nidx, nblk, max_dolic_c, max_dolic_e
    INTEGER :: edge_index, level, blockNo, jc,  cell_connect, sum_lsm_quad_edge, ctr
    INTEGER :: nidx1, nblk1, nidx2, nblk2, nidx3, nblk3
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc
    !-------------------------------------------------------------------------
!     CALL message("","limiter_ocean_zalesak_horizontal_onTriangles is running...")

    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------
    start_level = 1
!     end_level   = n_zlev
    cellOfEdge_idx  => patch_2d%edges%cell_idx
    cellOfEdge_blk  => patch_2d%edges%cell_blk
    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
    neighbor_cell_idx => patch_2d%cells%neighbor_idx
    neighbor_cell_blk => patch_2d%cells%neighbor_blk
    div_coeff => operators_coefficients%div_coeff
    edges_SeaBoundaryLevel => operators_coefficients%edges_SeaBoundaryLevel

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(patch_3d%p_patch_2d(1)%nblks_e, patch_3d%p_patch_2d(1)%alloc_cell_blocks) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%dolic_e, patch_3d%p_patch_1d(1)%dolic_c) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c, patch_3d%p_patch_1d(1)%inv_prism_thick_c) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%del_zlev_m, div_coeff) &
    !$ACC   COPYIN(cellOfEdge_idx, cellOfEdge_blk, edge_of_cell_idx, edge_of_cell_blk) &
    !$ACC   COPYIN(neighbor_cell_idx, neighbor_cell_blk, edges_SeaBoundaryLevel) &
    !$ACC   COPYIN(tracer, flx_tracer_low, flx_tracer_high, div_adv_flux_vert, h_old, h_new) &
    !$ACC   CREATE(z_mflx_anti1, z_mflx_anti2, z_mflx_anti3, z_fluxdiv_c, z_anti) &
    !$ACC   CREATE(z_tracer_new_low, z_tracer_max, z_tracer_min, r_p, r_m, z_tracer_update_horz) &
    !$ACC   CREATE(z_min, z_max, p_p, p_m, inv_prism_thick_new, delta_z_new, delta_z) &
    !$ACC   COPY(flx_tracer_final) IF(lzacc)

#ifdef NAGFOR
   !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
   z_tracer_max(:,:,:) = 0.0_wp
   z_tracer_min(:,:,:) = 0.0_wp
   r_m(:,:,:)          = 0.0_wp
   r_p(:,:,:)          = 0.0_wp
   !$ACC END KERNELS
   !$ACC WAIT(1)
#endif
 
  IF (p_test_run) THEN
    z_tracer_max(:,:,:) = 0.0_wp
    z_tracer_min(:,:,:) = 0.0_wp
    r_m(:,:,:)          = 0.0_wp
    r_p(:,:,:)          = 0.0_wp
  ENDIF
  
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      z_anti(:,:,blockNo)     = 0.0_wp       
      !$ACC END KERNELS
      !$ACC WAIT(1)

      max_dolic_e = -1
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) REDUCTION(MAX: max_dolic_e) IF(lzacc)
      DO edge_index = start_index, end_index
        max_dolic_e = MAX(max_dolic_e, patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo))
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO level = start_level, max_dolic_e
        DO edge_index = start_index, end_index
          IF (patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo) < level) CYCLE
          ! calculate antidiffusive flux for each edge
          z_anti(edge_index,level,blockNo) = flx_tracer_high(edge_index,level,blockNo)&
                                          &- flx_tracer_low(edge_index,level,blockNo)
        END DO  ! end loop over edges
      END DO  ! end loop over levels
      !$ACC END PARALLEL
      !$ACC WAIT(1)
    END DO  ! end loop over blocks
!ICON_OMP_END_DO

    
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, delta_z, delta_z_new, &
!ICON_OMP z_fluxdiv_c) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      level = start_level
     
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = start_index, end_index
        IF (patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) < level) CYCLE
        
        ! 3. Compute the complete (with horizontal and vertical divergence) updated low order solution z_tracer_new_low
        ! First at top level than in fluid interior
              
        !  compute divergence of low order fluxes
        z_fluxdiv_c(jc) = &
            & flx_tracer_low(edge_of_cell_idx(jc,blockNo,1),level,edge_of_cell_blk(jc,blockNo,1)) * &
            & div_coeff(jc,level,blockNo,1) + & 
            & flx_tracer_low(edge_of_cell_idx(jc,blockNo,2),level,edge_of_cell_blk(jc,blockNo,2)) * &
            & div_coeff(jc,level,blockNo,2) + & 
            & flx_tracer_low(edge_of_cell_idx(jc,blockNo,3),level,edge_of_cell_blk(jc,blockNo,3)) * &
            & div_coeff(jc,level,blockNo,3) 
            
       
        delta_z(jc) = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)&
             &  + h_old(jc,blockNo)
        delta_z_new(jc) = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)&
             &  + h_new(jc,blockNo)
             
        ! Low order flux at top level
        z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z(jc)                     &
          & - dtime * (z_fluxdiv_c(jc)+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new(jc)

         z_tracer_update_horz(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z(jc)                     &
         & - dtime * (div_adv_flux_vert(jc,level,blockNo)))/delta_z_new(jc)

        z_tracer_max(jc,level,blockNo) =            &
          & MAX(          tracer(jc,level,blockNo), &
          &     z_tracer_new_low(jc,level,blockNo))
        z_tracer_min(jc,level,blockNo) =            &
          & MIN(          tracer(jc,level,blockNo), &
          &     z_tracer_new_low(jc,level,blockNo))
      ENDDO      
      !$ACC END PARALLEL
    ENDDO
    !$ACC WAIT(1)
!ICON_OMP_END_DO
         
    !Fluid interior       
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, delta_z, delta_z_new, &
!ICON_OMP z_fluxdiv_c) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)

      max_dolic_c = -1
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) REDUCTION(MAX: max_dolic_c) IF(lzacc)
      DO jc = start_index, end_index
        max_dolic_c = MAX(max_dolic_c, patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo))
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO level = start_level+1, max_dolic_c
        DO jc = start_index, end_index
          IF (patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) < level) CYCLE
          !  compute divergence of low order fluxes
          z_fluxdiv_c(jc) = &
            & flx_tracer_low(edge_of_cell_idx(jc,blockNo,1),level,edge_of_cell_blk(jc,blockNo,1)) * &
            & div_coeff(jc,level,blockNo,1) + & 
            & flx_tracer_low(edge_of_cell_idx(jc,blockNo,2),level,edge_of_cell_blk(jc,blockNo,2)) * &
            & div_coeff(jc,level,blockNo,2) + & 
            & flx_tracer_low(edge_of_cell_idx(jc,blockNo,3),level,edge_of_cell_blk(jc,blockNo,3)) * &
            & div_coeff(jc,level,blockNo,3) 
         
          delta_z(jc)     = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
          delta_z_new(jc) = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
              
          z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z(jc)                     &
            & - dtime * (z_fluxdiv_c(jc)+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new(jc)
                      
          z_tracer_max(jc,level,blockNo) =            &
            & MAX(          tracer(jc,level,blockNo), &
            &     z_tracer_new_low(jc,level,blockNo))
          z_tracer_min(jc,level,blockNo) =            &
            & MIN(          tracer(jc,level,blockNo), &
            &     z_tracer_new_low(jc,level,blockNo))
         
        ENDDO
      ENDDO
      !$ACC END PARALLEL
      !$ACC WAIT(1)
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, z_tracer_max, z_tracer_min)

    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, inv_prism_thick_new, &
!ICON_OMP z_mflx_anti1, z_mflx_anti2, z_mflx_anti3, z_max, z_min, cell_connect, p_p, p_m, nidx, nblk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

      ! this is only needed for the parallel test setups
      ! it will try  tocheck the uninitialized (land) parts
!       r_m(:,:,blockNo) = 0.0_wp
!       r_p(:,:,blockNo) = 0.0_wp
        
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      inv_prism_thick_new(:,:) = patch_3D%p_patch_1d(1)%inv_prism_thick_c(:,:,blockNo)
      !$ACC END KERNELS

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        IF (patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) >= start_level) &              
          inv_prism_thick_new(jc, start_level) = 1.0_wp / (patch_3d%p_patch_1d(1)%del_zlev_m(start_level)+ h_new(jc,blockNo))
      ENDDO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      max_dolic_c = -1
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) REDUCTION(MAX: max_dolic_c) IF(lzacc)
      DO jc = start_index, end_index
        max_dolic_c = MAX(max_dolic_c, patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo))
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO level = start_level, max_dolic_c
        DO jc = start_index, end_index
          IF (patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo) < level) CYCLE              
          
          ! 2. Define "antidiffusive" fluxes A(jc,level,blockNo,edge_index) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.
                    
          z_mflx_anti1(jc) =                                                        &
            & dtime * div_coeff(jc,level,blockNo,1) * inv_prism_thick_new(jc, level)  &
            & * z_anti(edge_of_cell_idx(jc,blockNo,1),level,edge_of_cell_blk(jc,blockNo,1))
          z_mflx_anti2(jc) =                                                        &
            & dtime * div_coeff(jc,level,blockNo,2) * inv_prism_thick_new(jc, level)  &
            & * z_anti(edge_of_cell_idx(jc,blockNo,2),level,edge_of_cell_blk(jc,blockNo,2))
          z_mflx_anti3(jc) =                                                        &
            & dtime * div_coeff(jc,level,blockNo,3) * inv_prism_thick_new(jc, level)  &
            & * z_anti(edge_of_cell_idx(jc,blockNo,3),level,edge_of_cell_blk(jc,blockNo,3))
                  
          z_max(jc) = z_tracer_max(jc,level,blockNo)
          z_min(jc) = z_tracer_min(jc,level,blockNo)
          p_p(jc) = 0.0_wp
          p_m(jc) = 0_wp

          nidx1 = neighbor_cell_idx(jc,blockNo,1)
          nblk1 = neighbor_cell_blk(jc,blockNo,1)
          nidx2 = neighbor_cell_idx(jc,blockNo,2)
          nblk2 = neighbor_cell_blk(jc,blockNo,2)
          nidx3 = neighbor_cell_idx(jc,blockNo,3)
          nblk3 = neighbor_cell_blk(jc,blockNo,3)

          IF (patch_3d%p_patch_1d(1)%dolic_c(nidx1, nblk1) >= level) THEN
            z_max(jc) = MAX(z_max(jc), &
              & z_tracer_max(nidx1,level,nblk1))
            z_min(jc) = MIN(z_min(jc), &
              & z_tracer_min(nidx1,level,nblk1))

            ! Sum of all incoming antidiffusive fluxes into cell jc
            ! outgoing fluxes carry a positive sign, incoming a negative
            p_p(jc) = p_p(jc) - MIN(0._wp, z_mflx_anti1(jc))
            ! Sum of all outgoing antidiffusive fluxes out of cell jc
            p_m(jc) = p_m(jc) + MAX(0._wp, z_mflx_anti1(jc))
          ENDIF
          IF (patch_3d%p_patch_1d(1)%dolic_c(nidx2, nblk2) >= level) THEN
            z_max(jc) = MAX(z_max(jc), &
              & z_tracer_max(nidx2,level,nblk2))
            z_min(jc) = MIN(z_min(jc), &
              & z_tracer_min(nidx2,level,nblk2))

            ! Sum of all incoming antidiffusive fluxes into cell jc
            ! outgoing fluxes carry a positive sign, incoming a negative
            p_p(jc) = p_p(jc) - MIN(0._wp, z_mflx_anti2(jc))
            ! Sum of all outgoing antidiffusive fluxes out of cell jc
            p_m(jc) = p_m(jc) + MAX(0._wp, z_mflx_anti2(jc))
          ENDIF
          IF (patch_3d%p_patch_1d(1)%dolic_c(nidx3, nblk3) >= level) THEN
            z_max(jc) = MAX(z_max(jc), &
              & z_tracer_max(nidx3,level,nblk3))
            z_min(jc) = MIN(z_min(jc), &
              & z_tracer_min(nidx3,level,nblk3))

            ! Sum of all incoming antidiffusive fluxes into cell jc
            ! outgoing fluxes carry a positive sign, incoming a negative
            p_p(jc) = p_p(jc) - MIN(0._wp, z_mflx_anti3(jc))
            ! Sum of all outgoing antidiffusive fluxes out of cell jc
            p_m(jc) = p_m(jc) + MAX(0._wp, z_mflx_anti3(jc))
          ENDIF
           
          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of tracer
          r_m(jc,level,blockNo) = (z_tracer_new_low(jc,level,blockNo) - z_min(jc) ) / (p_m(jc) + dbl_eps)!&
          !
          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of tracer
          r_p(jc,level,blockNo) = (z_max(jc) - z_tracer_new_low(jc,level,blockNo)) / (p_p(jc) + dbl_eps)!&
          !
          !update old tracer with low-order flux
          !!tracer(jc,level,blockNo)=z_tracer_new_low(jc,level,blockNo) 
          !tracer(jc,level,blockNo)=z_tracer_update_horz(jc,level,blockNo)         
        ENDDO
      ENDDO
      !$ACC END PARALLEL
      !$ACC WAIT(1)
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
    
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)

    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
!ICON_OMP_DO_PARALLEL PRIVATE(start_index, end_index, edge_index, level, z_signum, r_frac) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      flx_tracer_final(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS
      !$ACC WAIT(1)

      max_dolic_e = -1
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) REDUCTION(MAX: max_dolic_e) IF(lzacc)
      DO edge_index = start_index, end_index
        max_dolic_e = MAX(max_dolic_e, patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo))
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_signum, r_frac)
      DO level = start_level, max_dolic_e
        DO edge_index = start_index, end_index
          IF (patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo) < level) CYCLE
 
          IF( edges_SeaBoundaryLevel(edge_index,level,blockNo) > -2)THEN! edge < 2nd order boundary
          
            flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)
            
          ELSE!IF(sum_lsm_quad_edge==all_water_edges)THEN

            !z_anti>0 returns  1: here z_anti is outgoing, i.e. flux_high>flux_low
            !z_anti<0 returns -1: here z_anti is ingoing, i.e. flux_high<flux_low
            z_signum = SIGN(1._wp, z_anti(edge_index,level,blockNo))
                    
            ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
            ! but is computationally more efficient
            r_frac = 0.5_wp * (      &
              & (1._wp + z_signum) * & !<- active for z_signum=1
              & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)),  &
              &     r_p(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)))  &
              &+(1._wp - z_signum) * & !<- active for z_signum=-1
              & MIN(r_m(cellOfEdge_idx(edge_index,blockNo,2),level,cellOfEdge_blk(edge_index,blockNo,2)),  &
              &     r_p(cellOfEdge_idx(edge_index,blockNo,1),level,cellOfEdge_blk(edge_index,blockNo,1)))  )
            
            ! Limited flux
            flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo)&
            & + MIN(1.0_wp,r_frac) * z_anti(edge_index,level,blockNo)   
          
          ENDIF
          
        END DO
      ENDDO
      !$ACC END PARALLEL
      !$ACC WAIT(1)
    ENDDO
!ICON_OMP_END_DO_PARALLEL
! !ICON_OMP_END_PARALLEL
     !$ACC END DATA
  END SUBROUTINE limiter_ocean_zalesak_horizontal_onTriangles_lvector
  !-------------------------------------------------------------------------
#endif

  !-------------------------------------------------------------------------
  !>
  !! Limiter for PPM (3rd order) vertical advection (monotone version)
  !!
  !! Removes over- and undershoots in first guess parabola by resetting the
  !! upper or lower interface values.
  !! Avoids non-physical over/undershoots in advected fields.
  !!
  !! Note that this limiter was coded assuming a pressure based vertical
  !! coordinate system. Nevertheless this limiter works for a height based
  !! vertical system, too. This is due to a 'wrong' computation of z_delta
  !! in the case of a height based coordinate system (i.e. z_delta is
  !! implicity multiplied by -1)
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  !!
  !! mpi parallelized, only cells_in_domain are computed, no sync
!<Optimize:inUse>
  SUBROUTINE v_ppm_slimiter_mo_onBlock( p_cc, p_face, p_slope, p_face_up, p_face_low, &
    & startIndex, endIndex, cells_noOfLevels, lacc)

    REAL(wp), INTENT(in)           :: p_cc(nproma,n_zlev)      !< advected cell centered variable
    REAL(wp), INTENT(in)           :: p_face(nproma,n_zlev+1)  !< reconstructed face values of the advected field
    REAL(wp), INTENT(in)           :: p_slope(nproma,n_zlev+1) !< monotonized slope
    REAL(wp), INTENT(inout)        :: p_face_up(nproma,n_zlev) !< final face value (upper face, height based)
    REAL(wp), INTENT(inout)        :: p_face_low(nproma,n_zlev)!< final face value (lower face, height based)
    INTEGER,  INTENT(in)           :: startIndex, endIndex
    INTEGER,  INTENT(in)           :: cells_noOfLevels(nproma)
    LOGICAL, INTENT(in), OPTIONAL :: lacc

    ! locals
    INTEGER :: nlev                      !< number of full levels
    INTEGER :: firstLevel                      !< vertical start thisLevel
    INTEGER :: jc, jk                   !< index of cell, vertical thisLevel
    INTEGER :: ikp1                      !< vertical thisLevel plus one
    REAL(wp) :: z_delta                   !< lower minus upper face value
    REAL(wp) :: z_a6i                     !< curvature of parabola
    TYPE(t_patch), POINTER :: patch_2D
    INTEGER :: kmax
    LOGICAL :: lzacc
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    firstLevel = 1
    nlev = n_zlev

#ifdef _OPENACC
    kmax = maxval(cells_noOfLevels)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = firstLevel, kmax
      DO jc = startIndex, endIndex
        IF (jk <= cells_noOfLevels(jc)) THEN
#else
    DO jc = startIndex, endIndex
      DO jk = firstLevel, cells_noOfLevels(jc)
#endif
        ! index of bottom half thisLevel
        ikp1 = jk + 1

        z_delta   = p_face(jc,ikp1) - p_face(jc,jk)
        z_a6i     = 6._wp * (p_cc(jc,jk)                           &
          & - 0.5_wp * (p_face(jc,jk) + p_face(jc,ikp1)))

        IF ( p_slope(jc,jk) == 0._wp) THEN
          p_face_up(jc,jk)  = p_cc(jc,jk)
          p_face_low(jc,jk) = p_cc(jc,jk)

        ELSE IF (z_delta * z_a6i > z_delta * z_delta) THEN
          p_face_up(jc,jk)  = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,ikp1)
          p_face_low(jc,jk) = p_face(jc,ikp1)

        ELSE IF (z_delta * z_a6i < -1._wp * (z_delta * z_delta)) THEN
          p_face_up(jc,jk)  = p_face(jc,jk)
          p_face_low(jc,jk) = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,jk)

        ELSE
          p_face_up(jc,jk)  = p_face(jc,jk)
          p_face_low(jc,jk) = p_face(jc,ikp1)
        ENDIF
#ifdef _OPENACC
        END IF   
#endif
      END DO
        
    END DO
#ifdef _OPENACC
    !$ACC END PARALLEL
    !$ACC WAIT(1)
#endif

  END SUBROUTINE v_ppm_slimiter_mo_onBlock
  !-------------------------------------------------------------------------
 
END MODULE mo_ocean_limiter
