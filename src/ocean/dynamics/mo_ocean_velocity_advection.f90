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

! Contains the implementation of velocity advection in vector invariant form
! that is used in the ocean model.

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_velocity_advection
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: finish
  USE mo_sync,                ONLY: sync_e, sync_c, sync_v, sync_patch_array
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_impl_constants,      ONLY: boundary, min_dolic
  USE mo_ocean_nml,           ONLY: n_zlev,NonlinearCoriolis_type,&
    & nonlinear_coriolis_primal_grid,nonlinear_coriolis_dual_grid, no_coriolis, &
    & VerticalAdvection_None, VerticalAdvection_RotationalForm, VerticalAdvection_MimeticRotationalForm, &
    & VerticalAdvection_DivergenceForm, HorizonatlVelocity_VerticalAdvection_form
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_ocean_types,         ONLY: t_hydro_ocean_diag
  USE mo_ocean_math_operators,ONLY: grad_fd_norm_oce_3d_onBlock, &
    & rot_vertex_ocean_3d, verticalDeriv_vec_midlevel_on_block,  &
    & verticalDiv_vector_onFullLevels_on_block
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_math_utilities,      ONLY: vector_product
  USE mo_scalar_product,      ONLY: map_cell2edges_3D, nonlinear_coriolis_3D,map_vec_prismtop2center_on_block, &
    & map_scalar_prismtop2center_onBlock, map_vector_center2prismtop_onBlock
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: veloc_adv_horz_mimetic
  PUBLIC  :: veloc_adv_vert_mimetic

  CHARACTER(len=12)  :: str_module = 'oceVelocAdv '  ! Output of module for 1 line debug
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_ocean_velocity_advection'
  INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug
! CHARACTER(len=12)  :: str_module = '__FILE__'


  INTEGER, PARAMETER, PRIVATE :: rotational_form = 0
  INTEGER, PARAMETER, PRIVATE :: divergence_form = 1
  INTEGER, PARAMETER, PRIVATE :: velocity_advection_form=0

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Computes horizontal advection of a (edge based) vector field.
  !! either by using rotational/vector-invariant form of velocity advection
  !! or the divergence form
  !!
  !!   mpi parallelized LL
!<Optimize:inUse>
  SUBROUTINE veloc_adv_horz_mimetic( patch_3D,        &
    & vn_old,          &
    & vn_new,          &
    & p_diag,          &
    & veloc_adv_horz_e,&
    & ocean_coefficients, &
    & lacc)
    !
    !
    TYPE(t_patch_3D ), POINTER, INTENT(IN) :: patch_3D
    REAL(wp), POINTER, INTENT(inout)     :: vn_old(:,:,:)
    REAL(wp), POINTER, INTENT(inout)     :: vn_new(:,:,:)
    TYPE(t_hydro_ocean_diag)             :: p_diag
    REAL(wp), POINTER, INTENT(inout)     :: veloc_adv_horz_e(:,:,:) ! out
    TYPE(t_operator_coeff), TARGET, INTENT(in) :: ocean_coefficients
    LOGICAL, INTENT(IN), OPTIONAL              :: lacc

    LOGICAL :: lzacc
    CHARACTER(len=*), PARAMETER :: routine = modname//':veloc_adv_horz_mimetic'
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (velocity_advection_form == rotational_form) THEN


      IF(NonlinearCoriolis_type==nonlinear_coriolis_dual_grid)THEN
        ! inUse
        CALL veloc_adv_horz_mimetic_rot( patch_3D,        &
          & vn_old,          &
          & p_diag,          &
          & veloc_adv_horz_e,&
          & ocean_coefficients,&
          & lacc=lzacc)

      ELSEIF(NonlinearCoriolis_type==nonlinear_coriolis_primal_grid)THEN

#ifdef _OPENACC
        IF (lzacc) CALL finish(routine, 'OpenACC version currently for nonlinear_coriolis_primal_grid not implemented')
#endif

        CALL veloc_adv_horz_mimetic_classicCgrid( patch_3D, &
          & vn_old,          &
          & p_diag,          &
          & veloc_adv_horz_e,&
          & ocean_coefficients)

      ELSEIF(NonlinearCoriolis_type==no_coriolis)THEN

#ifdef _OPENACC
        IF (lzacc) CALL finish(routine, 'OpenACC version currently for no_coriolis not implemented')
#endif

        CALL calculate_only_kineticGrad( patch_3D, &
          & vn_old,          &
          & p_diag,          &
          & veloc_adv_horz_e,&
          & ocean_coefficients)
      ENDIF


    ELSEIF (velocity_advection_form == divergence_form) THEN

#ifdef _OPENACC
      IF (lzacc) CALL finish(routine, 'OpenACC version currently for divergence_form not implemented')
#endif

      ! notInUse
      CALL veloc_adv_horz_mimetic_div( patch_3D,      &
        & vn_old,        &
        & p_diag,        &
        & ocean_coefficients,    &
        & veloc_adv_horz_e)
    ENDIF

  END SUBROUTINE veloc_adv_horz_mimetic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) vector field.
  !! either by using rotational/vector-invariant form of velocity advection
  !! or the divergence form
  !!
!<Optimize:inUse>
  SUBROUTINE veloc_adv_vert_mimetic( patch_3D, p_diag,ocean_coefficients, veloc_adv_vert_e, lacc)
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff), INTENT(in):: ocean_coefficients
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(:,:,:) ! (1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)
    LOGICAL, INTENT(IN), OPTIONAL              :: lacc

    LOGICAL :: lzacc
    CHARACTER(len=*), PARAMETER :: routine = modname//':veloc_adv_vert_mimetic'
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    SELECT CASE(HorizonatlVelocity_VerticalAdvection_form)
    CASE(VerticalAdvection_MimeticRotationalForm)
      CALL veloc_adv_vert_mimetic_rot( patch_3D, p_diag,ocean_coefficients, veloc_adv_vert_e, lacc=lzacc)

    CASE(VerticalAdvection_DivergenceForm)

#ifdef _OPENACC
      IF (lzacc) CALL finish(routine, 'OpenACC version currently for VerticalAdvection_DivergenceForm not implemented')
#endif

      CALL veloc_adv_vert_mimetic_div( patch_3D, p_diag,ocean_coefficients, veloc_adv_vert_e)

    CASE(VerticalAdvection_RotationalForm)

#ifdef _OPENACC
      IF (lzacc) CALL finish(routine, 'OpenACC version currently for VerticalAdvection_RotationalForm not implemented')
#endif

      CALL veloc_adv_vert_rot( patch_3D, p_diag,ocean_coefficients, veloc_adv_vert_e)

    CASE(VerticalAdvection_None)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      veloc_adv_vert_e(:,:,:) = 0.0_wp
      !$ACC END KERNELS
      !$ACC WAIT(1)

    CASE default
      CALL finish("veloc_adv_vert_mimetic","unknown HorizonatlVelocity_VerticalAdvection_form")
    END SELECT

  END SUBROUTINE veloc_adv_vert_mimetic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes horizontal advection of a (edge based) vector field.
  !!
  !! Computes rotational term of a vector field given by its components in
  !! the directions normal to triangle edges and the gradient of the kinetic energy
  !! which is calculated using the reconstructed velocity at cell centers. Both
  !! terms are combined and constitute the horizontal velocity advection.
  !!
  !!IMPORTANT: It is assumed that the reconstruction of the tangential velocity
  !!           has been done before.
  !1
  !! input:  lives on edges (velocity points)
  !! output: lives on edges (velocity points)
  !!
  !!
  !! veloc_adv_horz_e is on edges%in_domain
  !! p_diag%vort is on all vertices
!<Optimize:inUse>
  SUBROUTINE veloc_adv_horz_mimetic_rot( patch_3D,     &
    & vn,              &
    & p_diag,          &
    & veloc_adv_horz_e,&
    & ocean_coefficients,&
    & lacc)

    TYPE(t_patch_3D ), POINTER, INTENT(IN)   :: patch_3D
    REAL(wp), POINTER, INTENT(inout)  :: vn(:,:,:)
    TYPE(t_hydro_ocean_diag) :: p_diag
    REAL(wp), POINTER, INTENT(inout)  :: veloc_adv_horz_e(:,:,:) ! out
    TYPE(t_operator_coeff), INTENT(in), TARGET :: ocean_coefficients
    LOGICAL, INTENT(IN), OPTIONAL              :: lacc

    INTEGER :: jk, blockNo, cell_index,start_cell_index, end_cell_index, level,startLevel
    INTEGER :: start_edge_index, end_edge_index
    LOGICAL :: lzacc
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    !REAL(wp) :: z_vort_flx(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells
    TYPE(t_patch), POINTER         :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain
    all_cells       => patch_2D%cells%all

    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
    startLevel      =1
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    !calculate vorticity flux across dual edge
    ! p_diag%p_vn_dual and vn_old must have been synced
    CALL nonlinear_coriolis_3d( patch_3D, &
      & vn,              &
      & p_diag%p_vn_dual,&
      & p_diag%vort,     &
      & ocean_coefficients,      &
      & veloc_adv_horz_e, &
      & lacc=lzacc)

      ! this is calculated in mo_diagnostics, only used for output
!     p_diag%potential_vort_e=veloc_adv_horz_e
!     !Diagnostic potential vorticity at cells: !This work currently only for triangles
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
!
!       DO cell_index =  start_cell_index, end_cell_index
!         DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
!
!         p_diag%potential_vort_c(cell_index, level, blockNo)&
!         &=(p_diag%potential_vort_e(edge_of_cell_idx(cell_index,blockNo,1),level,edge_of_cell_blk(cell_index,blockNo,1))&
!         &+p_diag%potential_vort_e(edge_of_cell_idx(cell_index,blockNo,2),level,edge_of_cell_blk(cell_index,blockNo,2))&
!         &+p_diag%potential_vort_e(edge_of_cell_idx(cell_index,blockNo,3),level,edge_of_cell_blk(cell_index,blockNo,3)))/3.0_wp
!
!         END DO
!       END DO
!     END DO
    !-------------------------------------------------------------------------------
    ! IF(L_ENSTROPHY_DISSIPATION)THEN
    !  DO jk = start_level, elev
    !   write(*,*)'max/min vort flux: ',MAXVAL(z_vort_flx(:,jk,:)),&
    !                                         &MINVAL(z_vort_flx(:,jk,:))
    !   END DO
    ! z_vort_flx=laplacian4vortex_flux(patch_2D,z_vort_flx)
    ! ENDIF
    !-------------------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      !calculate gradient of kinetic energy
      CALL grad_fd_norm_oce_3d_onBlock ( &
        & p_diag%kin,                    &
        & patch_3D,                    &
        & ocean_coefficients%grad_coeff(:,:,blockNo), &
        & p_diag%grad(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo, lacc=lzacc)
      ! the result is on edges_in_domain


    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('HorzMimRot: kin energy'        ,p_diag%kin              ,str_module,idt_src, &
          patch_2D%cells%owned )
    CALL dbg_print('HorzMimRot: vorticity'         ,p_diag%vort             ,str_module,idt_src, &
          patch_2D%verts%owned )
    CALL dbg_print('HorzMimRot: grad kin en'       ,p_diag%grad             ,str_module,idt_src, &
          patch_2D%edges%owned )
    !---------------------------------------------------------------------
    !idt_src=2  ! output print level (1-5, fix)
    !CALL dbg_print('HorzMimRot: final Vel.Adv.'    ,veloc_adv_horz_e        ,str_module,idt_src, &
    !      patch_2D%edges%owned )
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_horz_mimetic_rot
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE calculate_only_kineticGrad( patch_3D,     &
    & vn,              &
    & p_diag,          &
    & veloc_adv_horz_e,&
    & ocean_coefficients)

    TYPE(t_patch_3D ),TARGET   :: patch_3D
    REAL(wp), POINTER, INTENT(inout)  :: vn(:,:,:)
    TYPE(t_hydro_ocean_diag) :: p_diag
    REAL(wp), POINTER, INTENT(inout)  :: veloc_adv_horz_e(:,:,:) ! out
    TYPE(t_operator_coeff), INTENT(in):: ocean_coefficients

    INTEGER :: blockNo, start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER         :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-----------------------------------------------------------------------

    veloc_adv_horz_e = 0.0_wp

!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      !calculate gradient of kinetic energy
      CALL grad_fd_norm_oce_3d_onBlock ( &
        & p_diag%kin,                    &
        & patch_3D,                    &
        & ocean_coefficients%grad_coeff(:,:,blockNo), &
        & p_diag%grad(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo)
      ! the result is on edges_in_domain


    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('HorzMimRot: kin energy'        ,p_diag%kin              ,str_module,idt_src, &
          patch_2D%cells%owned )
    CALL dbg_print('HorzMimRot: grad kin en'       ,p_diag%grad             ,str_module,idt_src, &
          patch_2D%edges%owned )

  END SUBROUTINE calculate_only_kineticGrad
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes horizontal advection of a (edge based) vector field.
  !!
  !! Computes rotational term of a vector field given by its components in
  !! the directions normal to triangle edges and the gradient of the kinetic energy
  !! which is calculated using the reconstructed velocity at cell centers. Both
  !! terms are combined and constitute the horizontal velocity advection.
  !!
  !!IMPORTANT: It is assumed that the reconstruction of the tangential velocity
  !!           has been done before.
  !1
  !! input:  lives on edges (velocity points)
  !! output: lives on edges (velocity points)
  !!
  !!
  !! veloc_adv_horz_e is on edges%in_domain
  !! p_diag%vort is on all vertices
!<Optimize:inUse>
  SUBROUTINE veloc_adv_horz_mimetic_classicCgrid( patch_3D,     &
    & vn,              &
    & p_diag,          &
    & veloc_adv_horz_e,&
    & ocean_coefficients)
    !
    !
    !  patch on which computation is performed
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    REAL(wp), INTENT(inout)  :: vn(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)
    !REAL(wp), INTENT(inout) :: vn_new(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_hydro_ocean_diag) :: p_diag
    REAL(wp), INTENT(inout)  :: veloc_adv_horz_e(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e) ! out
    !
    TYPE(t_operator_coeff), INTENT(in):: ocean_coefficients
    INTEGER :: jk, blockNo, je, jc
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: c1_idx, c1_blk, c2_idx, c2_blk
    INTEGER :: v1_idx, v1_blk, v2_idx, v2_blk
    REAL(wp) :: veloc_tangential
    INTEGER, POINTER :: edge_levels(:,:)
    !REAL(wp) :: z_vort_flx(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_edges, all_cells
    TYPE(t_patch), POINTER         :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain
    all_cells => patch_2D%cells%all
    edge_levels => patch_3D%p_patch_1D(1)%dolic_e
    !-----------------------------------------------------------------------
    CALL rot_vertex_ocean_3d( patch_3d, vn, p_diag%p_vn_dual, ocean_coefficients, p_diag%vort)
    ! sync not needed here, but used for example for the Leith
    CALL sync_patch_array(SYNC_V, patch_2D, p_diag%vort)
    !--------------------------------------------------------------
    !calculate nonlinear coriolis term by
    !1) projection cell reconstructed velocity vector in tangential direction
    !2) averaging the result from 1) from two adjecent cells to an edge
    !3) multiplying the result by the averaged vorticity
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_edge_index,end_edge_index, je, jk, c1_idx, c1_blk, c2_idx, c2_blk, &
!ICON_OMP  v1_idx, v1_blk, v2_idx, v2_blk, veloc_tangential) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO je = start_edge_index, end_edge_index
        !Neighbouring cells
        c1_idx = patch_2D%edges%cell_idx(je,blockNo,1)
        c1_blk = patch_2D%edges%cell_blk(je,blockNo,1)
        c2_idx = patch_2D%edges%cell_idx(je,blockNo,2)
        c2_blk = patch_2D%edges%cell_blk(je,blockNo,2)
        !Neighbouring verts
        v1_idx = patch_2D%edges%vertex_idx(je,blockNo,1)
        v1_blk = patch_2D%edges%vertex_blk(je,blockNo,1)
        v2_idx = patch_2D%edges%vertex_idx(je,blockNo,2)
        v2_blk = patch_2D%edges%vertex_blk(je,blockNo,2)

        DO jk = 1, edge_levels(je,blockNo)
          !calculation of tangential velocity
          veloc_tangential= dot_product( &
            &   p_diag%p_vn(c1_idx,jk,c1_blk)%x * ocean_coefficients%averageCellsToEdges(je,blockNo,1)  &
            &  +p_diag%p_vn(c2_idx,jk,c2_blk)%x * ocean_coefficients%averageCellsToEdges(je,blockNo,2), &
            &   patch_2D%edges%dual_cart_normal(je,blockNo)%x)

!           !This is an upwind version of the nonlinear coriolis.
!           !Not recommended just for testing purposes
!           IF(patch_2D%edges%system_orientation(je,blockNo)==1.0_wp)THEN
!           !dual normal vector points from vertex 1 to vertex 2
!             IF(veloc_tangential>0.0_wp)THEN
!              veloc_adv_horz_e(je,jk,blockNo)=veloc_tangential&
!              &*(patch_2d%edges%f_e(je,blockNo)+p_diag%vort(il_v1,jk,ib_v1))
!             ELSE
!               veloc_adv_horz_e(je,jk,blockNo)=veloc_tangential&
!               &*(patch_2d%edges%f_e(je,blockNo)+p_diag%vort(il_v2,jk,ib_v2))
!             ENDIF

!           ELSEIF(patch_2D%edges%system_orientation(je,blockNo)==-1.0_wp)THEN
            !!dual normal vector points from vertex 2 to vertex 1
!             IF(veloc_tangential>0.0)THEN
!               veloc_adv_horz_e(je,jk,blockNo)=veloc_tangential&
!               &*(patch_2d%edges%f_e(je,blockNo)+p_diag%vort(il_v2,jk,ib_v2))
!             ELSE
!               veloc_adv_horz_e(je,jk,blockNo)=veloc_tangential&
!               &*(patch_2d%edges%f_e(je,blockNo)+p_diag%vort(il_v1,jk,ib_v1))
!             ENDIF
!           ENDIF

            !calculation of nonlinear Coriolis
          veloc_adv_horz_e(je,jk,blockNo) = veloc_tangential &
            & * (patch_2d%edges%f_e(je,blockNo)  &
            & + 0.5_wp * (p_diag%vort(v1_idx,jk,v1_blk)+p_diag%vort(v2_idx,jk,v2_blk)))

        END DO
      END DO
    END DO
!ICON_OMP_END_DO

 !   DO blockNo = all_cells%start_block, all_cells%end_block
 !     CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
 !     DO jk = 1, n_zlev
 !       DO jc = start_cell_index, end_cell_index
 !         p_diag%kin(jc,jk,blockNo)= 0.5_wp*(&
 !         & p_diag%u(jc,jk,blockNo)*p_diag%u(jc,jk,blockNo)&
 !          &+p_diag%v(jc,jk,blockNo)*p_diag%v(jc,jk,blockNo))
 !       END DO
 !     END DO
 !   END DO
 !  CALL sync_patch_array(SYNC_C, patch_2D, p_diag%kin)



!ICON_OMP_DO PRIVATE(start_edge_index,end_edge_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      !calculate gradient of kinetic energy
      CALL grad_fd_norm_oce_3d_onBlock ( &
        & p_diag%kin,                    &
        & patch_3D,                    &
        & ocean_coefficients%grad_coeff(:,:,blockNo), &
        & p_diag%grad(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo)
    ! the result is on edges_in_domain

    END DO ! blocks
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
!     idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('advHorCgrid: kin energy'        ,p_diag%kin              ,str_module,idt_src, &
          patch_2D%cells%owned )
    CALL dbg_print('advHorCgrid: vorticity'         ,p_diag%vort             ,str_module,idt_src, &
          patch_2D%verts%owned )
    CALL dbg_print('advHorCgrid: f_e'  ,        patch_2d%edges%f_e        ,str_module,idt_src, &
          patch_2D%edges%owned )
    CALL dbg_print('advHorCgrid: grad kin en'       ,p_diag%grad             ,str_module,idt_src, &
          patch_2D%edges%owned )
    CALL dbg_print('advHorCgrid: veloc_adv_horz_e'  ,veloc_adv_horz_e        ,str_module,idt_src, &
          patch_2D%edges%owned )
    !---------------------------------------------------------------------
    !idt_src=2  ! output print level (1-5, fix)
    !CALL dbg_print('HorzMimRot: final Vel.Adv.'    ,veloc_adv_horz_e        ,str_module,idt_src, &
    !      patch_2D%edges%owned )
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_horz_mimetic_classicCgrid
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE veloc_adv_horz_mimetic_div( patch_3D,   &
    & vn,             &
    & p_diag,         &
    & ocean_coefficients,     &
    & veloc_adv_horz_e)
    !

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    REAL(wp), INTENT(inout)          :: vn(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e) ! dim: (nproma,n_zlev,nblks_e)
    TYPE(t_hydro_ocean_diag)         :: p_diag
    TYPE(t_operator_coeff),INTENT(in):: ocean_coefficients
    REAL(wp), INTENT(inout)          :: veloc_adv_horz_e(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e) ! out
    !
    !Local variables
    !
    INTEGER :: start_level, elev     ! vertical start and end level
    INTEGER :: jk, blockNo, je, jc
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_index_c, end_index_c
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: il_e1, ib_e1, il_e2, ib_e2, il_e3, ib_e3

    REAL(wp)                      :: z_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp)                      :: veloc_tangential(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp)                      :: div_veloc(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: u_v_cc_v(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: u_v_cc_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: u_v_cc_c(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: z_div_vec_c(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER         :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)

    ! #slo# set local variable to zero due to nag -nan compiler-option
    all_edges => patch_2D%edges%all
    all_cells => patch_2D%cells%all

    veloc_tangential(1:nproma,1:n_zlev,1:patch_2D%nblks_e) = 0.0_wp
    veloc_adv_horz_e(1:nproma,1:n_zlev,1:patch_2D%nblks_e) = 0.0_wp

        div_veloc(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)=0.0_wp

    start_level = 1
    elev = n_zlev

    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
      DO jk = start_level, elev
        DO je = start_edge_index, end_edge_index

          IF(patch_3D%lsm_e(je,jk,blockNo)< boundary)THEN
            !Neighbouring cells

            il_c1 = patch_2D%edges%cell_idx(je,blockNo,1)
            ib_c1 = patch_2D%edges%cell_blk(je,blockNo,1)
            il_c2 = patch_2D%edges%cell_idx(je,blockNo,2)
            ib_c2 = patch_2D%edges%cell_blk(je,blockNo,2)

                        !u_v_cc_e(je,jk,blockNo)%x= 0.5_wp*(p_diag%p_vn(il_c1,jk,ib_c1)%x + p_diag%p_vn(il_c2,jk,ib_c2)%x)
            u_v_cc_e(je,jk,blockNo)%x=&
                        &(  patch_2D%edges%edge_cell_length(je,blockNo,1)*p_diag%p_vn(il_c1,jk,ib_c1)%x&
                        & + patch_2D%edges%edge_cell_length(je,blockNo,2)*p_diag%p_vn(il_c2,jk,ib_c2)%x)&
                        &/(patch_2D%edges%dual_edge_length(je,blockNo))

                        z_e(je,jk,blockNo)=p_diag%vn_time_weighted(je,jk,blockNo)&
                        &*dot_product(u_v_cc_e(je,jk,blockNo)%x,patch_2D%edges%primal_cart_normal(je,blockNo)%x)

            veloc_tangential(je,jk,blockNo)=patch_2d%edges%f_e(je,blockNo)*&
            dot_product(u_v_cc_e(je,jk,blockNo)%x,patch_2D%edges%dual_cart_normal(je,blockNo)%x)

          ENDIF
        END DO
      END DO
    END DO

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index_c, end_index_c)
      DO jk = start_level, elev
        DO jc = start_index_c, end_index_c

         il_e1 = patch_2D%cells%edge_idx(jc,blockNo,1)
         ib_e1 = patch_2D%cells%edge_blk(jc,blockNo,1)

         il_e2 = patch_2D%cells%edge_idx(jc,blockNo,2)
         ib_e2 = patch_2D%cells%edge_blk(jc,blockNo,2)

         il_e3 = patch_2D%cells%edge_idx(jc,blockNo,3)
         ib_e3 = patch_2D%cells%edge_blk(jc,blockNo,3)

         div_veloc(jc,jk,blockNo)=  &
              & (z_e(il_e1,jk,ib_e1)* ocean_coefficients%div_coeff(jc,jk,blockNo,1)&
                          &*patch_3D%p_patch_1d(1)%prism_thick_e(il_e1,jk,ib_e1) + &
              & z_e(il_e2,jk,ib_e2) * ocean_coefficients%div_coeff(jc,jk,blockNo,2)&
                          &*patch_3D%p_patch_1d(1)%prism_thick_e(il_e2,jk,ib_e2) + &
              & z_e(il_e3,jk,ib_e3) * ocean_coefficients%div_coeff(jc,jk,blockNo,3)&
                          &*patch_3D%p_patch_1d(1)%prism_thick_e(il_e3,jk,ib_e3))&
                          &/patch_3D%p_patch_1d(1)%prism_thick_c(jc,jk,blockNo)


 !         z_div_vec_c(jc,jk,blockNo)%x =  &
!               & (u_v_cc_e(il_e1,jk,ib_e1)%x* ocean_coefficients%div_coeff(jc,jk,blockNo,1)&
!                         &*patch_3D%p_patch_1d(1)%prism_thick_e(il_e1,jk,ib_e1) + &
!               & u_v_cc_e(il_e2,jk,ib_e2)%x * ocean_coefficients%div_coeff(jc,jk,blockNo,2)&
!                         &*patch_3D%p_patch_1d(1)%prism_thick_e(il_e2,jk,ib_e2) + &
!               & u_v_cc_e(il_e3,jk,ib_e3)%x * ocean_coefficients%div_coeff(jc,jk,blockNo,3)&
!                         &*patch_3D%p_patch_1d(1)%prism_thick_e(il_e3,jk,ib_e3))&
!                         & / patch_3D%p_patch_1d(1)%prism_thick_c(jc,jk,blockNo)

        END DO
      END DO
    END DO



!     CALL map_cell2edges( patch_2D, z_div_vec_c, veloc_adv_horz_e, &
!       & opt_start_level=start_level, opt_elev=elev )
!    CALL map_cell2edges_3D( patch_3D, z_div_vec_c, veloc_adv_horz_e,ocean_coefficients)



DO blockNo = all_edges%start_block, all_edges%end_block
  CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
  DO jk = start_level, elev
    DO je = start_edge_index, end_edge_index

      IF(patch_3D%lsm_e(je,jk,blockNo)<= boundary)THEN
        !Neighbouring cells
        il_c1 = patch_2D%edges%cell_idx(je,blockNo,1)
        ib_c1 = patch_2D%edges%cell_blk(je,blockNo,1)
        il_c2 = patch_2D%edges%cell_idx(je,blockNo,2)
        ib_c2 = patch_2D%edges%cell_blk(je,blockNo,2)

!               veloc_adv_horz_e(je,jk,blockNo)=veloc_adv_horz_e(je,jk,blockNo)+veloc_tangential(je,jk,blockNo)
                veloc_adv_horz_e(je,jk,blockNo)=&
                &(  patch_2D%edges%edge_cell_length(je,blockNo,1)*div_veloc(il_c1,jk,ib_c1) &
                & + patch_2D%edges%edge_cell_length(je,blockNo,2)*div_veloc(il_c2,jk,ib_c2))&
                &/patch_2D%edges%dual_edge_length(je,blockNo)&
                !&/(patch_2D%edges%edge_cell_length(je,blockNo,1)+patch_2D%edges%edge_cell_length(je,blockNo,2))&
                &+veloc_tangential(je,jk,blockNo)
!               write(123,*)'dist', patch_2D%edges%edge_cell_length(je,blockNo,1), patch_2D%edges%edge_cell_length(je,blockNo,2),&
!               &patch_2D%edges%dual_edge_length(je,blockNo),&
!               &patch_2D%edges%edge_cell_length(je,blockNo,1)/patch_2D%edges%dual_edge_length(je,blockNo),&
!               &patch_2D%edges%edge_cell_length(je,blockNo,2)/patch_2D%edges%dual_edge_length(je,blockNo)
      ENDIF
    END DO
  END DO
END DO

DO jk=1,n_zlev
write(*,*)'ADV',jk,maxval(veloc_adv_horz_e(:,jk,:)),minval(veloc_adv_horz_e(:,jk,:)),&
&maxval(veloc_tangential(:,jk,:)),minval(veloc_tangential(:,jk,:))
ENDDO


    !calculates the curl. This is needed in Laplace-beltrami operator (velocity diffusion).
    !It is not needed for velocity advection.
    CALL rot_vertex_ocean_3D( patch_3D,              &
                              & p_diag%vn_time_weighted,&!vn,                  &
                              & p_diag%p_vn_dual,    &
                              & ocean_coefficients,          &
                              & p_diag%vort)
    CALL sync_patch_array(SYNC_V, patch_2D, p_diag%vort)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('HorzMimDiv: final Vel.Adv.'    ,veloc_adv_horz_e        ,str_module,idt_src, &
          patch_2D%edges%owned )
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('HorzMimDiv: vorticity'         ,p_diag%vort             ,str_module,idt_src, &
          patch_2D%verts%owned )
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_horz_mimetic_div
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) horizontal vector field that
  !! suits to rotational form of velocity equation.
  !!
  !! IMPORTANT: It is assumed that the velocity vector reconstruction from
  !! edges to cells has been done before.
  !!
  !! input:  lives on cells (velocity points)
  !! output: lives on edges (velocity points)
  !!
!<Optimize:inUse>
  SUBROUTINE veloc_adv_vert_rot( patch_3D, p_diag,ocean_coefficients, veloc_adv_vert_e)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(in) :: ocean_coefficients
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)

    !local variables
!     INTEGER :: start_level     ! vertical start and end level
    INTEGER :: jc, jk, blockNo
    INTEGER :: start_index, end_index
    INTEGER :: fin_level
    REAL(wp), POINTER :: inv_prism_center_distance(:,:)! ,prism_thick(:,:)
    TYPE(t_cartesian_coordinates) :: z_adv_u_fullLevels(nproma,n_zlev)
    TYPE(t_cartesian_coordinates) :: vn_halfLevels(nproma,n_zlev+1) !, z_adv_u_fullLevel(nproma,n_zlev)
    TYPE(t_cartesian_coordinates) :: z_adv_u_m(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: center_vertical_velocity(nproma,n_zlev) !, z_adv_u_fullLevel(nproma,n_zlev)
!     TYPE(t_cartesian_coordinates) :: vertDeriv_vec(nproma, n_zlev)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: patch_2D
    REAL(wp), POINTER             :: vertical_velocity(:,:,:)! , prism_center_distance(:,:,:)
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_cells => patch_2D%cells%all
!     prism_center_distance => patch_3D%p_patch_1D(1)%prism_center_dist_c(:,:,:)
   !-----------------------------------------------------------------------
!     start_level = 1

    vertical_velocity => p_diag%w

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index,jc, jk, center_vertical_velocity, &
!ICON_OMP vn_halfLevels, z_adv_u_fullLevels) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)

      ! map w on full levels
      CALL map_scalar_prismtop2center_onBlock(patch_3d, vertical_velocity(:,:,blockNo), center_vertical_velocity, &
        & blockNo, start_index, end_index)

      ! map p_vn on half levels
      CALL map_vector_center2prismtop_onBlock(patch_3d, p_diag%p_vn(:,:,blockNo), vn_halfLevels, &
        & blockNo, start_index, end_index)

      !vertical derivative at ocean interior Surface is handled below
      ! this does not include h
      CALL verticalDiv_vector_onFullLevels_on_block( patch_3d, vn_halfLevels, z_adv_u_fullLevels, &
        & 1, blockNo, start_index, end_index)

      ! multiply vertical velocity with vertical derivative of horizontal velocity
      DO jc = start_index, end_index
        DO jk = 1, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)
          z_adv_u_m(jc,jk,blockNo)%x =  &
            center_vertical_velocity(jc,jk) * z_adv_u_fullLevels(jc,jk)%x
        END DO
        DO jk = patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)+1, n_zlev
          z_adv_u_m(jc,jk,blockNo)%x = 0.0_wp
        END DO
      END DO

    END DO
!ICON_OMP_END_PARALLEL_DO

    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('vn 1%x(1)'    ,p_diag%p_vn(:,1,:)%x(1)        ,str_module,idt_src, &
          patch_2D%cells%owned )
    CALL dbg_print('vn 1%x(2)'    ,p_diag%p_vn(:,1,:)%x(2)        ,str_module,idt_src, &
          patch_2D%cells%owned )
    CALL dbg_print('vn 1%x(3)'    ,p_diag%p_vn(:,1,:)%x(3)        ,str_module,idt_src, &
          patch_2D%cells%owned )
    CALL dbg_print('vn 2%x(1)'    ,p_diag%p_vn(:,2,:)%x(1)        ,str_module,idt_src, &
          patch_2D%cells%owned )
    CALL dbg_print('vn 2%x(2)'    ,p_diag%p_vn(:,2,:)%x(2)        ,str_module,idt_src, &
          patch_2D%cells%owned )
    CALL dbg_print('vn 2%x(3)'    ,p_diag%p_vn(:,2,:)%x(3)        ,str_module,idt_src, &
          patch_2D%cells%owned )

    CALL dbg_print('VertAdvect x(1)'    ,z_adv_u_m(:,:,:)%x(1)        ,str_module,idt_src, &
          patch_2D%cells%owned )
    CALL dbg_print('VertAdvect x(2)'    ,z_adv_u_m(:,:,:)%x(2)        ,str_module,idt_src, &
          patch_2D%cells%owned )
    CALL dbg_print('VertAdvect x(3)'    ,z_adv_u_m(:,:,:)%x(3)        ,str_module,idt_src, &
          patch_2D%cells%owned )


    ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
    CALL map_cell2edges_3D( patch_3D, z_adv_u_m, veloc_adv_vert_e,ocean_coefficients)

!     CALL sync_patch_array(SYNC_E, patch_2D, veloc_adv_vert_e)
    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('VertMimRot: V.Adv. Final'    ,veloc_adv_vert_e         ,str_module,idt_src, &
          patch_2D%edges%owned )
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_vert_rot
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) horizontal vector field that
  !! suits to rotational form of velocity equation.
  !! The vertical derivative of the velocity vector at circumcenters that
  !! is reconstructed from edge data is calculated and then multiplied by
  !! the vertical velocity. The product is mapped from top of the computational
  !! prism to the middle (still at centers) via the transposed of vertical differentiation
  !! and then transformed to edges.
  !!
  !! IMPORTANT: It is assumed that the velocity vector reconstruction from
  !! edges to cells has been done before.
  !!
  !! input:  lives on cells (velocity points)
  !! output: lives on edges (velocity points)
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE veloc_adv_vert_mimetic_rot( patch_3D, p_diag,p_op_coeff, veloc_adv_vert_e, lacc)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)
    LOGICAL, INTENT(IN), OPTIONAL     :: lacc

    !local variables
    INTEGER :: start_level     ! vertical start and end level
    INTEGER :: jc, jk, blockNo
    INTEGER :: start_index, end_index
    INTEGER :: fin_level
    LOGICAL :: lzacc
    REAL(wp), POINTER :: inv_prism_center_distance(:,:)! ,prism_thick(:,:)
    TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1)
    TYPE(t_cartesian_coordinates) :: z_adv_u_m(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
!     TYPE(t_cartesian_coordinates) :: vertDeriv_vec(nproma, n_zlev)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: patch_2D
    REAL(wp), POINTER             :: vertical_velocity(:,:,:)
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_cells => patch_2D%cells%all
   !-----------------------------------------------------------------------
    start_level = 1

    CALL set_acc_host_or_device(lzacc, lacc)

    vertical_velocity => p_diag%w

    !$ACC DATA CREATE(z_adv_u_i, z_adv_u_m) IF(lzacc)

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(1) = 0.0_wp
    z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(2) = 0.0_wp
    z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(3) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index,jc, jk, fin_level,inv_prism_center_distance, &
!ICON_OMP z_adv_u_i) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)

      !vertical derivative at ocean interior Surface is handled below
      ! this does not include h
      CALL verticalDeriv_vec_midlevel_on_block( patch_3d, &
                                              & p_diag%p_vn(:,:,blockNo),  &
                                              & z_adv_u_i(:,:),&
                                              & start_level+1,             &
                                              & blockNo, start_index, end_index, lacc=lzacc)

      !Step 1: multiply vertical velocity with vertical derivative of horizontal velocity
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        fin_level = patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)

        IF(fin_level >= min_dolic) THEN

          !1a) ocean surface: vertical derivative times vertical velocity.
          !This form is consistent with energy conservation
!           z_adv_u_i(jc,start_level)%x =               &
!             & p_diag%w(jc,start_level,blockNo)*&  !/v_base%del_zlev_i(slev)
!             & (p_diag%p_vn(jc,start_level,blockNo)%x - p_diag%p_vn(jc,start_level+1,blockNo)%x)&!/del_zlev_i(slev)
!             & * inv_prism_center_distance(jc,start_level)
           z_adv_u_i(jc,start_level)%x =               &
             & -vertical_velocity(jc,start_level,blockNo) * p_diag%p_vn(jc,start_level,blockNo)%x &
             & * patch_3D%p_patch_1D(1)%constantPrismCenters_invZdistance(jc,start_level,blockNo)

          ! 1b) ocean interior
          DO jk = start_level+1, fin_level
            z_adv_u_i(jc,jk)%x =  vertical_velocity(jc,jk,blockNo) * z_adv_u_i(jc,jk)%x
          END DO

!           z_adv_u_i(jc,fin_level)%x   = 0.0_wp
          z_adv_u_i(jc,fin_level+1)%x = 0.0_wp

        ENDIF
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      ! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
      CALL map_vec_prismtop2center_on_block(patch_3d, z_adv_u_i, z_adv_u_m(:,:,blockNo), &
        & blockNo, start_index, end_index, lacc=lzacc)
    END DO

!ICON_OMP_END_PARALLEL_DO

    ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
    CALL map_cell2edges_3D( patch_3D, z_adv_u_m, veloc_adv_vert_e,p_op_coeff, lacc=lzacc)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('VertMimRot: V.Adv. Final'    ,veloc_adv_vert_e         ,str_module,idt_src, &
          patch_2D%edges%owned )
    !---------------------------------------------------------------------

    !$ACC END DATA
  END SUBROUTINE veloc_adv_vert_mimetic_rot
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) horizontal vector field that
  !! suits to rotational form of velocity equation.
  !! The vertical derivative of the velocity vector at circumcenters that
  !! is reconstructed from edge data is calculated and then multiplied by
  !! the vertical velocity. The product is mapped from top of the computational
  !! prism to the middle (still at centers) via the transposed of vertical differentiation
  !! and then transformed to edges.
  !!
  !! IMPORTANT: It is assumed that the velocity vector reconstruction from
  !! edges to cells has been done before.
  !!
  !! input:  lives on cells (velocity points)
  !! output: lives on edges (velocity points)
  !!
  !!
!   SUBROUTINE veloc_adv_vert_mimetic_rot_old( patch_3D, p_diag,ocean_coefficients, veloc_adv_vert_e)
!
!     TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
!     TYPE(t_hydro_ocean_diag)          :: p_diag
!     TYPE(t_operator_coeff),INTENT(in) :: ocean_coefficients
!     REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)
!
!     !local variables
!     INTEGER :: start_level     ! vertical start and end level
!     INTEGER :: jc, jk, blockNo
!     INTEGER :: start_index, end_index
!     INTEGER :: end_level
!     REAL(wp), POINTER :: prism_center_distance(:,:),prism_thick(:,:)
!     TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1)
!     TYPE(t_cartesian_coordinates) :: z_adv_u_m(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
!     TYPE(t_subset_range), POINTER :: all_cells
!     TYPE(t_patch), POINTER        :: patch_2D
!     REAL(wp), POINTER             :: vert_veloc(:,:,:)
!     !-----------------------------------------------------------------------
!     patch_2D   => patch_3D%p_patch_2D(1)
!     all_cells => patch_2D%cells%all
!    !-----------------------------------------------------------------------
!     start_level = 1
!
!     z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(1) = 0.0_wp
!     z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(2) = 0.0_wp
!     z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(3) = 0.0_wp
!
!     z_adv_u_i(1:nproma,1:n_zlev)%x(1) = 0.0_wp
!     z_adv_u_i(1:nproma,1:n_zlev)%x(2) = 0.0_wp
!     z_adv_u_i(1:nproma,1:n_zlev)%x(3) = 0.0_wp
!
!     vert_veloc =>p_diag%w
!
!
! !ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index,jc, jk, end_level,prism_center_distance, &
! !ICON_OMP z_adv_u_i) ICON_OMP_DEFAULT_SCHEDULE
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!
!       prism_center_distance => patch_3D%p_patch_1D(1)%prism_center_dist_c(:,:,blockNo)
!       prism_thick           => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,blockNo)
!
!
!       !Step 1: multiply vertical velocity with vertical derivative of horizontal velocity
!       DO jc = start_index, end_index
!         end_level = patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)
!
!         IF(end_level >= min_dolic) THEN
!
!           !1a) ocean surface: vertical derivative times vertical velocity.
!           !This form is consistent with energy conservation
!
!           z_adv_u_i(jc,start_level)%x =               &
!             & -vert_veloc(jc,start_level,blockNo)*p_diag%p_vn(jc,start_level+1,blockNo)%x &
!             & / prism_center_distance(jc,start_level)
!
!           ! 1b) ocean interior
!           DO jk = start_level+1, end_level-1
!             z_adv_u_i(jc,jk)%x                  &
!              & = vert_veloc(jc,jk,blockNo)* &
!              &   (p_diag%p_vn(jc,jk-1,blockNo)%x - p_diag%p_vn(jc,jk,blockNo)%x) &
!              &    / prism_center_distance(jc,jk)
!           END DO
!           z_adv_u_i(jc,end_level)%x = 0.0_wp
!
!         ENDIF
!       END DO
!       !----------------------------------------------------------------------------------
!
!       ! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
!        DO jc = start_index, end_index
!         end_level = patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)
!         DO jk = start_level,end_level-1
!           z_adv_u_m(jc,jk,blockNo)%x &
!           & = (prism_center_distance(jc,jk)   * z_adv_u_i(jc,jk)%x    &
!           & +  prism_center_distance(jc,jk+1) * z_adv_u_i(jc,jk+1)%x) &
!           & / (2.0_wp*prism_thick(jc,jk))!(prism_center_distance(jc,jk+1) + prism_center_distance(jc,jk))
!         END DO
!       END DO
!
!     END DO
! !ICON_OMP_END_PARALLEL_DO
!     ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
!     CALL map_cell2edges_3D( patch_3D, z_adv_u_m, veloc_adv_vert_e,ocean_coefficients)
!
!     !---------Debug Diagnostics-------------------------------------------
!     idt_src=3  ! output print level (1-5, fix)
!     CALL dbg_print('VertMimRot: V.Adv. Final'    ,veloc_adv_vert_e         ,str_module,idt_src, &
!           patch_2D%edges%owned )
!     !---------------------------------------------------------------------
!
!   END SUBROUTINE veloc_adv_vert_mimetic_rot_old
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) horizontal vector field that
  !! suits to rotational form of velocity equation. The rotational form excludes
  !! the flux form of the vertical advection, instead it implies the use of
  !! $w \partial_z Pv$. In this subroutine the vertical advection is discretized
  !! by transforming into flux-form minus a correction term,
  !! All calculations are carried out at cell centers and are
  !! mapped to edges at the end.
  !! The vertical derivative of the velocity vector at circumcenters that
  !! is reconstructed from edge data is calculated and then multiplied by
  !! the vertical velocity. The product is mapped from top of the computational
  !! prism to the middle (still at centers) via the transposed of vertical differentiation
  !! and then transformed to edges.
  !!
  !! IMPORTANT: It is assumed that the velocity vector reconstruction from
  !! edges to cells has been done before.
  !!
  !! input:  lives on cells (velocity points)
  !! output: lives on edges (velocity points)
  !!
  !!
  SUBROUTINE veloc_adv_vert_mim_rot_flux2( patch_3D, p_diag,ocean_coefficients, veloc_adv_vert_e)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(in) :: ocean_coefficients
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,patch_3D%p_patch_2D(1)%nblks_e)

    !local variables
    INTEGER :: start_level, elev     ! vertical start and end level
    INTEGER :: jc, jk, blockNo
    INTEGER :: start_index, end_index
    INTEGER :: fin_level
    REAL(wp), POINTER :: prism_center_distance(:)
    REAL(wp), POINTER :: del_zlev_m(:)
    REAL(wp)                      :: z_w_ave  (nproma,n_zlev,  patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                      :: z_w_diff (nproma,n_zlev-1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: z_adv_u_m(nproma,n_zlev,  patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_cells => patch_2D%cells%all
    !-----------------------------------------------------------------------

    start_level = 1
    elev = n_zlev

    z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(1) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(2) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(3) = 0.0_wp

    z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(1) = 0.0_wp
    z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(2) = 0.0_wp
    z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(3) = 0.0_wp

    z_w_ave (1:nproma,1:n_zlev,  1:patch_2D%alloc_cell_blocks) = 0.0_wp
    z_w_diff(1:nproma,1:n_zlev-1,1:patch_2D%alloc_cell_blocks) = 0.0_wp


    !Step 1: multiply vertical velocity with vertical derivative of horizontal velocity
    !This requires appropriate boundary conditions
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        fin_level = patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)!v_base%dolic_c(jc,blockNo)

        IF(fin_level>=min_dolic)THEN
          !del_zlev_m=>p_diag%inv_prism_thick_c(jc,:,blockNo)
          del_zlev_m => patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,:,blockNo)
          DO jk = start_level, fin_level-1
            z_w_ave(jc,jk,blockNo) = 0.5_wp*       (p_diag%w(jc,jk,blockNo)+p_diag%w(jc,jk+1,blockNo))
            z_w_diff(jc,jk,blockNo)=del_zlev_m(jk)*(p_diag%w(jc,jk,blockNo)-p_diag%w(jc,jk+1,blockNo))
             !&p_diag%inv_prism_thick_c(jc,jk,blockNo)!/v_base%del_zlev_m(jk)
          END DO
          z_w_ave(jc,fin_level,blockNo)= p_diag%w(jc,fin_level,blockNo)
        !ELSE
        !  z_w_ave(jc,:,blockNo)=0.0_wp
        ENDIF
      END DO
    END DO

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        fin_level = patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)!v_base%dolic_c(jc,blockNo)
        IF(fin_level>=min_dolic)THEN
          !prism_center_distance=>p_diag%inv_prism_center_dist_c(jc,:,blockNo)
          ! prism_center_distance => patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,:,blockNo)
          prism_center_distance => patch_3D%p_patch_1D(1)%constantPrismCenters_invZdistance(jc,:,blockNo)
          DO jk = start_level,fin_level-1
            !The last term is the correction term to the constructed flux form.
            !The result of this calculation lives on prism top or bottom.
            z_adv_u_i(jc,jk,blockNo)%x =                            &
              & (z_w_ave(jc,jk,blockNo)  *p_diag%p_vn(jc,jk,blockNo)%x   &
              & -z_w_ave(jc,jk+1,blockNo)*p_diag%p_vn(jc,jk+1,blockNo)%x)&
              &*prism_center_distance(jk)                                 &!!/v_base%prism_center_distance(jk)&
              & -z_w_diff(jc,jk,blockNo) * p_diag%p_vn(jc,jk,blockNo)%x
        END DO
        ENDIF
      END DO
    END DO
! write(*,*)'vert 2: A max/min vert adv:',jk,&
! & maxval(z_adv_u_i(:,jk,:)%x(1)), minval(z_adv_u_m(:,jk,:)%x(1)),&
! & maxval(z_adv_u_i(:,jk,:)%x(2)), minval(z_adv_u_m(:,jk,:)%x(2)),&
! & maxval(z_adv_u_i(:,jk,:)%x(3)), minval(z_adv_u_m(:,jk,:)%x(3))

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        fin_level = patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)!v_base%dolic_c(jc,blockNo)
        IF(fin_level>=min_dolic)THEN
          !prism_center_distance=>p_diag%prism_center_dist_c(jc,:,blockNo)
          ! prism_center_distance=>patch_3D%p_patch_1D(1)%prism_center_dist_c(jc,:,blockNo)
          prism_center_distance=>patch_3D%p_patch_1D(1)%constantPrismCenters_Zdistance(jc,:,blockNo)

          DO jk = start_level,fin_level-1!DO jk = start_level+1,fin_level-1
            !This seems to work well
            !z_adv_u_m(jc,jk,blockNo)%x &
            !& = 0.5_wp*(z_adv_u_i(jc,jk,blockNo)%x+z_adv_u_i(jc,jk+1,blockNo)%x)
            !Map back from cell top to cell center and weight according to thickness
            !  z_adv_u_m(jc,jk,blockNo)%x &
            !    & = (v_base%prism_center_distance(jk)  *z_adv_u_i(jc,jk,blockNo)%x&
            !    & +  v_base%prism_center_distance(jk+1)*z_adv_u_i(jc,jk+1,blockNo)%x) &
            !    & / (v_base%prism_center_distance(jk+1)+v_base%prism_center_distance(jk))
            z_adv_u_m(jc,jk,blockNo)%x &
              & = (prism_center_distance(jk)  *z_adv_u_i(jc,jk,blockNo)%x&
              & +  prism_center_distance(jk+1)*z_adv_u_i(jc,jk+1,blockNo)%x) &
              & / (prism_center_distance(jk+1)+prism_center_distance(jk))
          END DO
          ! 2c) ocean bottom
          !z_adv_u_m(jc,fin_level,blockNo)%x =0.0_wp!=  z_adv_u_i(jc,fin_level,blockNo)%x
        ENDIF
      END DO
    END DO

    ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
    CALL map_cell2edges_3D( patch_3D, z_adv_u_m, veloc_adv_vert_e,ocean_coefficients)


    CALL sync_patch_array(SYNC_E, patch_2D, veloc_adv_vert_e)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    !CALL dbg_print('VertMimRot2: z_adv_u_m%x(1)' ,z_adv_u_m%x(1)           ,str_module,idt_src)
    CALL dbg_print('VertMimRot2: VelAdv Final'   ,veloc_adv_vert_e         ,str_module,idt_src, &
          patch_2D%edges%owned )
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_vert_mim_rot_flux2
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) horizontal vector field that
  !! suits to rotational form of velocity equation. The rotational form excludes
  !! the flux form of the vertical advection, instead it implies the use of
  !! $w \partial_z Pv$. In this subroutine the vertical advection is discretized
  !! by transforming into flux-form minus a correction term,
  !! All calculations are carried out at cell centers and are
  !! mapped to edges at the end.
  !! The vertical derivative of the velocity vector at circumcenters that
  !! is reconstructed from edge data is calculated and then multiplied by
  !! the vertical velocity. The product is mapped from top of the computational
  !! prism to the middle (still at centers) via the transposed of vertical differentiation
  !! and then transformed to edges.
  !!
  !! IMPORTANT: It is assumed that the velocity vector reconstruction from
  !! edges to cells has been done before.
  !!
  !! input:  lives on cells (velocity points)
  !! output: lives on edges (velocity points)
  !!
  !!
  SUBROUTINE veloc_adv_vert_mimetic_rot_flux( patch_3D, p_diag,ocean_coefficients, veloc_adv_vert_e)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(in) :: ocean_coefficients
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,patch_3D%p_patch_2D(1)%nblks_e)

    !local variables
    INTEGER :: start_level, elev     ! vertical start and end level
    INTEGER :: jc, jk, blockNo
    INTEGER :: start_index, end_index
    INTEGER :: fin_level
    REAL(wp), POINTER :: del_zlev_m(:)
    REAL(wp)                      :: z_w_diff (nproma,n_zlev-1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_cells => patch_2D%cells%all
    !-----------------------------------------------------------------------
    start_level = 1
    elev = n_zlev

    z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(1) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(2) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(3) = 0.0_wp

    z_w_diff(1:nproma,1:n_zlev-1,1:patch_2D%alloc_cell_blocks) = 0.0_wp

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        fin_level = patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)!v_base%dolic_c(jc,blockNo)
        IF(fin_level>=min_dolic)THEN
          !del_zlev_m=>p_diag%inv_prism_thick_c(jc,:,blockNo)
          del_zlev_m=>patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,:,blockNo)

          DO jk = start_level, fin_level-1
            z_w_diff(jc,jk,blockNo)=del_zlev_m(jk)*(p_diag%w(jc,jk,blockNo)-p_diag%w(jc,jk+1,blockNo))
             !/v_base%del_zlev_m(jk)
          END DO

          !jk=start_level
          !z_adv_u_i(jc,jk,blockNo)%x = del_zlev_m(jk)*0.5_wp*&
          !& (p_diag%w(jc,jk,blockNo)  *(p_diag%p_vn(jc,jk,blockNo)%x + p_diag%p_vn(jc,jk,blockNo)%x) &
          !& -p_diag%w(jc,jk+1,blockNo)*(p_diag%p_vn(jc,jk,blockNo)%x + p_diag%p_vn(jc,jk+1,blockNo)%x)) &
          !& -z_w_diff(jc,jk,blockNo) * p_diag%p_vn(jc,jk,blockNo)%x

          DO jk = start_level+1,fin_level-1
            !The last term is the correction term to the constructed flux form.
            !The result of this calculation lives on prism top or bottom.
            z_adv_u_i(jc,jk,blockNo)%x = del_zlev_m(jk)*0.5_wp*&
              & (p_diag%w(jc,jk,blockNo)  *(p_diag%p_vn(jc,jk-1,blockNo)%x + p_diag%p_vn(jc,jk,blockNo)%x)    &
              & -p_diag%w(jc,jk+1,blockNo)*(p_diag%p_vn(jc,jk,blockNo)%x   + p_diag%p_vn(jc,jk+1,blockNo)%x)) &
              & -z_w_diff(jc,jk,blockNo) * p_diag%p_vn(jc,jk,blockNo)%x
          END DO
        ENDIF
      END DO
    END DO
    ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
    CALL map_cell2edges_3D( patch_3D, z_adv_u_i, veloc_adv_vert_e, ocean_coefficients)

    CALL sync_patch_array(SYNC_E, patch_2D, veloc_adv_vert_e)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=1  ! output print level (1-5, fix)
    !CALL dbg_print('VertMimRot2: z_adv_u_m%x(1)' ,z_adv_u_m%x(1)           ,str_module,idt_src)
    CALL dbg_print('VertMimRot: VelAdv Final'   ,veloc_adv_vert_e(:,:,:)         ,str_module,idt_src, &
          patch_2D%edges%owned )
    !---------------------------------------------------------------------
  END SUBROUTINE veloc_adv_vert_mimetic_rot_flux
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes vertical advection of a (edge based) horizontal vector field that
  !! suits to rotational form of velocity equation. The rotational form excludes
  !! the flux form of the vertical advection, instead it implies the use of
  !! $w \partial_z Pv$. In this subroutine the vertical advection is discretized
  !! by transforming into flux-form minus a correction term,
  !! All calculations are carried out at cell centers and are
  !! mapped to edges at the end.
  !! The vertical derivative of the velocity vector at circumcenters that
  !! is reconstructed from edge data is calculated and then multiplied by
  !! the vertical velocity. The product is mapped from top of the computational
  !! prism to the middle (still at centers) via the transposed of vertical differentiation
  !! and then transformed to edges.
  !!
  !! IMPORTANT: It is assumed that the velocity vector reconstruction from
  !! edges to cells has been done before.
  !!
  !! input:  lives on cells (velocity points)
  !! output: lives on edges (velocity points)
  !!
  !!
  SUBROUTINE veloc_adv_vert_mimetic_div( patch_3D, p_diag,ocean_coefficients, veloc_adv_vert_e)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(in) :: ocean_coefficients
    REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,patch_3D%p_patch_2D(1)%nblks_e)

    !local variables
    INTEGER :: start_level, elev     ! vertical start and end level
    INTEGER :: jc, jk, blockNo
    INTEGER :: start_index, end_index
    INTEGER :: fin_level
    REAL(wp), POINTER :: del_zlev_m(:)
    REAL(wp)                      :: z_w_diff (nproma,n_zlev-1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_cells => patch_2D%cells%all
    !-----------------------------------------------------------------------
    start_level = 1
    elev = n_zlev

    z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(1) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(2) = 0.0_wp
    z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(3) = 0.0_wp

    z_w_diff(1:nproma,1:n_zlev-1,1:patch_2D%alloc_cell_blocks) = 0.0_wp

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        fin_level = patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)!v_base%dolic_c(jc,blockNo)
        IF(fin_level>=min_dolic)THEN
          !del_zlev_m=>p_diag%inv_prism_thick_c(jc,:,blockNo)
          del_zlev_m=>patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,:,blockNo)

          DO jk = start_level, fin_level-1
            z_w_diff(jc,jk,blockNo)=del_zlev_m(jk)*(p_diag%w(jc,jk,blockNo)-p_diag%w(jc,jk+1,blockNo))
             !/v_base%del_zlev_m(jk)
          END DO

          jk=start_level
          z_adv_u_i(jc,jk,blockNo)%x = del_zlev_m(jk)*0.5_wp*&
          & (p_diag%w(jc,jk,blockNo)*(p_diag%p_vn(jc,jk,blockNo)%x + p_diag%p_vn(jc,jk,blockNo)%x) &
          & -p_diag%w(jc,jk+1,blockNo)*(p_diag%p_vn(jc,jk,blockNo)%x + p_diag%p_vn(jc,jk+1,blockNo)%x))

          DO jk = start_level+1,fin_level-1
            !The last term is the correction term to the constructed flux form.
            !The result of this calculation lives on prism top or bottom.
            z_adv_u_i(jc,jk,blockNo)%x = del_zlev_m(jk)*0.5_wp*&
              & (p_diag%w(jc,jk,blockNo)*(p_diag%p_vn(jc,jk-1,blockNo)%x + p_diag%p_vn(jc,jk,blockNo)%x)    &
              & -p_diag%w(jc,jk+1,blockNo)*(p_diag%p_vn(jc,jk,blockNo)%x + p_diag%p_vn(jc,jk+1,blockNo)%x))
          END DO
        ENDIF
      END DO
    END DO
    ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
    CALL map_cell2edges_3D( patch_3D, z_adv_u_i, veloc_adv_vert_e, ocean_coefficients)

    CALL sync_patch_array(SYNC_E, patch_2D, veloc_adv_vert_e)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    !CALL dbg_print('VertMimRot2: z_adv_u_m%x(1)' ,z_adv_u_m%x(1)           ,str_module,idt_src)
    CALL dbg_print('VertMimDiv: VelAdv Final'   ,veloc_adv_vert_e         ,str_module,idt_src, &
          patch_2D%edges%owned )
    !---------------------------------------------------------------------

  END SUBROUTINE veloc_adv_vert_mimetic_div
  !-------------------------------------------------------------------------

!   !-------------------------------------------------------------------------
!   !>
!   !! Computes vertical advection of a (edge based) horizontal vector field that
!   !! suits to divergence form of velocity equation.
!   !! The vertical derivative of the velocity vector at circumcenters that
!   !! is reconstructed from edge data is calculated and then multiplied by
!   !! the vertical velocity. The product is mapped from top of the computational
!   !! prism to the middle (still at centers) via the transposed of vertical differentiation
!   !! and then transformed to edges.
!   !!
!   !! IMPORTANT: It is assumed that the velocity vector reconstruction from
!   !! edges to cells has been done before.
!   !!
!   !! input:  lives on cells (velocity points)
!   !! output: lives on edges (velocity points)
!   !!
!   !!  mpi parallelized LL
!   !!
!   SUBROUTINE veloc_adv_vert_mimetic_div( patch_2D, p_diag,ocean_coefficients, veloc_adv_vert_e)
!
!     TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
!     TYPE(t_hydro_ocean_diag)          :: p_diag
!     TYPE(t_operator_coeff),INTENT(in) :: ocean_coefficients
!     REAL(wp), INTENT(inout)           :: veloc_adv_vert_e(1:nproma,1:n_zlev,patch_2D%nblks_e)
!
!     !local variables
!     INTEGER :: start_level, elev     ! vertical start and end level
!     INTEGER :: jc, jk, blockNo
!     INTEGER :: start_index, end_index
!     INTEGER :: fin_level
!     TYPE(t_cartesian_coordinates) :: z_adv_u_i(nproma,n_zlev+1,patch_2D%alloc_cell_blocks)
!     TYPE(t_cartesian_coordinates) :: z_adv_u_m(nproma,n_zlev,patch_2D%alloc_cell_blocks)
!     TYPE(t_subset_range), POINTER :: all_cells
!     !-----------------------------------------------------------------------
!     all_cells => patch_2D%cells%all
!     !-----------------------------------------------------------------------
!     start_level = 1
!     elev = n_zlev
!
!     z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(1) = 0.0_wp
!     z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(2) = 0.0_wp
!     z_adv_u_i(1:nproma,1:n_zlev+1,1:patch_2D%alloc_cell_blocks)%x(3) = 0.0_wp
!
!     z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(1) = 0.0_wp
!     z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(2) = 0.0_wp
!     z_adv_u_m(1:nproma,1:n_zlev,1:patch_2D%alloc_cell_blocks)%x(3) = 0.0_wp
!
!
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!       DO jc = start_index, end_index
!         fin_level = v_base%dolic_c(jc,blockNo)
!
!         IF(fin_level>=min_dolic)THEN
!           ! 1a) ocean surface  !Code below explodes: Use upper boundary condition for d_z u ?
!           z_adv_u_i(jc,start_level,blockNo)%x = p_diag%w(jc,start_level,blockNo)*p_diag%p_vn(jc,start_level,blockNo)%x
!
!           ! 1b) ocean interior
!           DO jk = start_level+1, fin_level-1
!             z_adv_u_i(jc,jk,blockNo)%x&
!               & = 0.5_wp*p_diag%w(jc,jk,blockNo)*(p_diag%p_vn(jc,jk-1,blockNo)%x + p_diag%p_vn(jc,jk,blockNo)%x)
!           END DO
!         ENDIF
!       END DO
!     END DO
!
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!       DO jc = start_index, end_index
!         fin_level = v_base%dolic_c(jc,blockNo)
!         IF(fin_level>=min_dolic)THEN
!           DO jk = start_level,fin_level-1
!             z_adv_u_m(jc,jk,blockNo)%x &
!               & = (z_adv_u_i(jc,jk,blockNo)%x+z_adv_u_i(jc,jk+1,blockNo)%x)/v_base%del_zlev_m(jk)
!           END DO
!         ENDIF
!       END DO
!       ! write(*,*)'B max/min vert adv:',jk, maxval(z_adv_u_m(:,jk,:)), minval(z_adv_u_m(:,jk,:)),&
!       ! & maxval(z_adv_v_m(:,jk,:)), minval(z_adv_v_m(:,jk,:))
!     END DO
!
!     ! ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
!     CALL map_cell2edges_3D( patch_2D, z_adv_u_m, veloc_adv_vert_e,ocean_coefficients)
!
!     CALL sync_patch_array(SYNC_E, patch_2D, veloc_adv_vert_e)
!
!     !---------Debug Diagnostics-------------------------------------------
!     idt_src=3  ! output print level (1-5, fix)
!     CALL dbg_print('VertMimDiv: z_adv_u_m%x(1)'  ,z_adv_u_m%x(1)           ,str_module,idt_src)
!     CALL dbg_print('VertMimDiv: VelAdv Final'    ,veloc_adv_vert_e         ,str_module,idt_src)
!     !---------------------------------------------------------------------
!
!   END SUBROUTINE veloc_adv_vert_mimetic_div
!   !-------------------------------------------------------------------------

!   !-------------------------------------------------------------------------
!   !>
!   !! Computes horizontal advection of a (edge based) vector field.
!   !!
!   !! Computes rotational term of a vector field given by its components in
!   !! the directions normal to triangle edges and the gradient of the kinetic energy
!   !! which is calculated using the reconstructed velocity at cell centers. Both
!   !! terms are combined and constitute the horizontal velocity advection.
!   !!
!   !!IMPORTANT: It is assumed that the reconstruction of the tangential velocity
!   !!           has been done before.
!   !1
!   !! input:  lives on edges (velocity points)
!   !! output: lives on edges (velocity points)
!   !!
!   !!  mpi parallelized LL
!   !!
!   SUBROUTINE veloc_adv_horz_rbf( patch_2D, vn, p_diag, grad_coeff, veloc_adv_horz_e, p_int)
!     !
!     !
!     !  patch on which computation is performed
!     !
!     TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
!
!     !
!     ! normal and tangential velocity  of which advection is computed
!     !
!     REAL(wp), INTENT(inout) :: vn(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
!     !
!     !diagnostic ocean state stores horizontally advected velocity
!     !
!     TYPE(t_hydro_ocean_diag) :: p_diag
!
!     REAL(wp), INTENT(in)    :: grad_coeff(:,:,:)
!     !
!     ! variable in which horizontally advected velocity is stored
!     !
!     REAL(wp), INTENT(inout) :: veloc_adv_horz_e(:,:,:)
!     !
!     !Interpolation necessary just for testing
!     TYPE(t_int_state),TARGET,INTENT(in)  :: p_int
!
!
!     INTEGER :: start_level, elev     ! vertical start and end level
!     INTEGER :: jk, blockNo, je, jc
!     INTEGER :: start_edge_index, end_edge_index
!     INTEGER :: start_index_c, end_index_c
!     INTEGER :: i_v1_idx, i_v1_blk, i_v2_idx, i_v2_blk
!     INTEGER :: jev, ile,ibe, i_v1_ctr, i_v2_ctr
!     REAL(wp) :: z_e  (nproma,n_zlev,patch_2D%nblks_e)
!     REAL(wp) :: z_vort_glb(nproma,n_zlev,patch_2D%nblks_v)
!     REAL(wp) :: z_grad_ekin_rbf(nproma,n_zlev,patch_2D%nblks_e)
!     REAL(wp) :: z_kin_rbf_e(nproma,n_zlev,patch_2D%nblks_e)
!     INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3
!     REAL(wp) :: z_vort_e(nproma,n_zlev,patch_2D%nblks_e)
!     REAL(wp) :: z_vort_flx_rbf(nproma,n_zlev,patch_2D%nblks_e)
!     REAL(wp) :: z_kin_e_rbf(nproma,n_zlev,patch_2D%nblks_e)
!     REAL(wp) :: z_weight_e1,z_weight_e2, z_weight_e3!, z_weight
!     TYPE(t_subset_range), POINTER :: all_edges, owned_edges, all_cells
!     !-----------------------------------------------------------------------
!     all_edges   => patch_2D%edges%all
!     owned_edges => patch_2D%edges%owned
!     all_cells   => patch_2D%cells%all
!
!     ! #slo# set local variable to zero due to nag -nan compiler-option
!     z_e             (:,:,:) = 0.0_wp
!     z_vort_glb      (:,:,:) = 0.0_wp
!     !z_vort_flx      (:,:,:) = 0.0_wp
!     veloc_adv_horz_e(:,:,:) = 0.0_wp
!     z_vort_e        (:,:,:) = 0.0_wp
!     z_vort_flx_rbf  (:,:,:) = 0.0_wp
!     z_kin_e_rbf     (:,:,:) = 0.0_wp
!     z_grad_ekin_rbf (:,:,:) = 0.0_wp
!
!     start_level = 1
!     elev = n_zlev
!
!     CALL rbf_vec_interpol_edge( vn,       &
!       & patch_2D,  &
!       & p_int,    &
!       & p_diag%vt,&
!       & opt_start_level=start_level,opt_elev=elev)
!
!     CALL sync_patch_array(SYNC_E, patch_2D, p_diag%v)
!
!
!     DO blockNo = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
!       DO jk = start_level, elev
!         DO je=start_edge_index, end_edge_index
!           IF ( v_base%lsm_e(je,jk,blockNo) == boundary ) THEN
!             p_diag%vt(je,jk,blockNo) = 0.0_wp
!             vn(je,jk,blockNo) = 0.0_wp
!           ENDIF
!         END DO
!       END DO
!     END DO
!
!     CALL rot_vertex_ocean_rbf(patch_2D,vn, p_diag%vt, p_diag%vort)
!     ! CALL verts2edges_scalar( p_diag%vort, patch_2D, p_int%v_1o2_e, &
!     !                          z_vort_e, opt_start_level=start_level,opt_elev=elev, opt_rlstart=3)
!     CALL sync_patch_array(SYNC_V, patch_2D, p_diag%vort)
!
!
!     DO blockNo = owned_edges%start_block, owned_edges%end_block
!       CALL get_index_range(owned_edges, blockNo, start_edge_index, end_edge_index)
!       DO jk = start_level, start_level
!         DO je=start_edge_index, end_edge_index
!           i_v1_idx = patch_2D%edges%vertex_idx(je,blockNo,1)
!           i_v1_blk = patch_2D%edges%vertex_blk(je,blockNo,1)
!           i_v2_idx = patch_2D%edges%vertex_idx(je,blockNo,2)
!           i_v2_blk = patch_2D%edges%vertex_blk(je,blockNo,2)
!           !count wet edges in vertex 1
!           i_v1_ctr = 0
!           DO jev = 1, patch_2D%verts%num_edges(i_v1_idx,i_v1_blk)
!             ile = patch_2D%verts%edge_idx(i_v1_idx,i_v1_blk,jev)
!             ibe = patch_2D%verts%edge_blk(i_v1_idx,i_v1_blk,jev)
!             IF ( v_base%lsm_e(ile,jk,ibe) == sea ) THEN
!               i_v1_ctr = i_v1_ctr +1
!             ENDIF
!           END DO
!           !count wet edges in vertex 2
!           i_v2_ctr = 0
!           DO jev = 1, patch_2D%verts%num_edges(i_v2_idx,i_v2_blk)
!             ile = patch_2D%verts%edge_idx(i_v2_idx,i_v2_blk,jev)
!             ibe = patch_2D%verts%edge_blk(i_v2_idx,i_v2_blk,jev)
!             IF ( v_base%lsm_e(ile,jk,ibe) == sea ) THEN
!               i_v2_ctr = i_v2_ctr +1
!             ENDIF
!           END DO
!           IF(   i_v1_ctr==patch_2D%verts%num_edges(i_v1_idx,i_v1_blk)&
!             & .AND.i_v2_ctr==patch_2D%verts%num_edges(i_v2_idx,i_v2_blk))THEN
!
!             z_vort_e(je,jk,blockNo) =&
!               & 0.5_wp*(p_diag%vort(i_v1_idx,jk,i_v1_blk)&
!               & +        p_diag%vort(i_v2_idx,jk,i_v2_blk))
!
!           ELSEIF(   i_v1_ctr==patch_2D%verts%num_edges(i_v1_idx,i_v1_blk)&
!             & .AND.i_v2_ctr <patch_2D%verts%num_edges(i_v2_idx,i_v2_blk))THEN
!
!             z_vort_e(je,jk,blockNo) = (REAL(i_v1_ctr,wp)*p_diag%vort(i_v1_idx,jk,i_v1_blk)&
!               & + REAL(i_v2_ctr,wp)*p_diag%vort(i_v2_idx,jk,i_v2_blk))/REAL(i_v1_ctr+i_v2_ctr,wp)
!
!           ELSEIF(   i_v1_ctr<patch_2D%verts%num_edges(i_v1_idx,i_v1_blk)&
!             & .AND.i_v2_ctr==patch_2D%verts%num_edges(i_v2_idx,i_v2_blk))THEN
!
!             z_vort_e(je,jk,blockNo) = (REAL(i_v1_ctr,wp)*p_diag%vort(i_v1_idx,jk,i_v1_blk)&
!               & + REAL(i_v2_ctr,wp)*p_diag%vort(i_v2_idx,jk,i_v2_blk))/REAL(i_v1_ctr+i_v2_ctr,wp)
!
!           ELSEIF(   i_v1_ctr<patch_2D%verts%num_edges(i_v1_idx,i_v1_blk)&
!             & .AND.i_v2_ctr<patch_2D%verts%num_edges(i_v2_idx,i_v2_blk))THEN
!
!             z_vort_e(je,jk,blockNo) = (REAL(i_v1_ctr,wp)*p_diag%vort(i_v1_idx,jk,i_v1_blk)&
!               & + REAL(i_v2_ctr,wp)*p_diag%vort(i_v2_idx,jk,i_v2_blk))/REAL(i_v1_ctr+i_v2_ctr,wp)
!           ELSE
!             z_vort_e(je,jk,blockNo) = 0.0_wp
!           ENDIF
!           !  IF(jk==1)THEN
!           !  IF(v_base%lsm_e(je,jk,blockNo)/=2)THEN
!           ! ! IF(v_base%lsm_e(je,jk,blockNo)==0)THEN
!           ! IF(i_v1_ctr==6.and.i_v2_ctr==6.and.v_base%lsm_e(je,jk,blockNo)==sea)THEN
!           ! ELSE
!           !  write(101,*)'vert ctr', jk,je,blockNo, i_v1_ctr, i_v2_ctr, v_base%lsm_e(je,jk,blockNo)!,&
!           ! ENDIF
!           ! ! ENDIF
!           !  ENDIF
!           !  ENDIF
!         END DO
!       END DO
!     ENDDO
!     CALL sync_patch_array(SYNC_E, patch_2D, z_vort_e)
!
!
!     DO blockNo = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
!       DO jk = start_level, elev
!         DO je=start_edge_index, end_edge_index
!           IF ( v_base%lsm_e(je,jk,blockNo) == sea ) THEN
!             z_vort_flx_rbf(je,jk,blockNo) =&
!               & p_diag%vt(je,jk,blockNo)*(patch_2D%edges%f_e(je,blockNo)+z_vort_e(je,jk,blockNo))
!             !          & p_diag%vt(je,jk,blockNo)*patch_2D%edges%f_e(je,blockNo)
!           ELSE
!             z_vort_flx_rbf(je,jk,blockNo) = 0.0_wp
!           ENDIF
!         END DO
!       END DO
!     ENDDO
!
!     CALL rbf_vec_interpol_cell( vn, patch_2D, p_int, p_diag%u,  &
!       & p_diag%v, opt_start_level=start_level, opt_elev=elev)
!     CALL sync_patch_array(SYNC_C, patch_2D, p_diag%v)
!
!     !write(*,*)'max/min vort flux:', MAXVAL(z_vort_flx_RBF(:,1,:)),MINVAL(z_vort_flx_RBF(:,1,:))
!     DO blockNo = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
!       DO jk = start_level, elev
!         DO je = start_edge_index, end_edge_index
!           ! calculate kinetic energy at edges from normal and tangential comp.
!           z_kin_rbf_e(je,jk,blockNo) =0.5_wp*(p_diag%vt(je,jk,blockNo)*p_diag%vt(je,jk,blockNo)&
!             & +    vn(je,jk,blockNo)*vn(je,jk,blockNo) )
!         ENDDO
!       ENDDO
!     ENDDO
!
!     !!$OMP END DO
!     !!$OMP END PARALLEL
!     ! Bilinear interpolation of kinetic energy from the edges to the cells
!     !    CALL edges2cells_scalar( z_kin_RBF_e,     &
!     !                           & patch_2D,         &
!     !                           & p_int%e_bln_c_s, &
!     !                           & p_diag%kin,      &
!     !                           & opt_start_level=start_level,opt_elev=elev)
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index_c, end_index_c)
!       DO jk = start_level, elev
!         DO jc = start_index_c, end_index_c
!           IF ( v_base%lsm_c(jc,jk,blockNo) > sea_boundary ) THEN
!             p_diag%kin(jc,jk,blockNo) = 0.0_wp
!           ELSE
!
!             ile1 = patch_2D%cells%edge_idx(jc,blockNo,1)
!             ibe1 = patch_2D%cells%edge_blk(jc,blockNo,1)
!             ile2 = patch_2D%cells%edge_idx(jc,blockNo,2)
!             ibe2 = patch_2D%cells%edge_blk(jc,blockNo,2)
!             ile3 = patch_2D%cells%edge_idx(jc,blockNo,3)
!             ibe3 = patch_2D%cells%edge_blk(jc,blockNo,3)
!             z_weight_e1 = 0.0_wp
!             z_weight_e2 = 0.0_wp
!             z_weight_e3 = 0.0_wp
!             IF(v_base%lsm_e(ile1,jk,ibe1)<= boundary)THEN
!               z_weight_e1 = patch_2D%edges%area_edge(ile1,ibe1)
!             ENDIF
!             IF(v_base%lsm_e(ile2,jk,ibe2)<= boundary)THEN
!               z_weight_e2 = patch_2D%edges%area_edge(ile2,ibe2)
!             ENDIF
!             IF(v_base%lsm_e(ile3,jk,ibe3)<= boundary)THEN
!               z_weight_e3 = patch_2D%edges%area_edge(ile3,ibe3)
!             ENDIF
!
!             !write(*,*)'weights',jc,jk,blockNo,z_weight_e1,z_weight_e2,z_weight_e3
!             p_diag%kin(jc,jk,blockNo) = (z_kin_rbf_e(ile1,jk,ibe1)*z_weight_e1&
!               & + z_kin_rbf_e(ile2,jk,ibe2)*z_weight_e2&
!               & + z_kin_rbf_e(ile3,jk,ibe3)*z_weight_e3)&
!               & /(z_weight_e1+z_weight_e2+z_weight_e3)
!             !       p_diag%kin(jc,jk,blockNo) = 0.5_wp*(p_diag%u(jc,jk,blockNo)*p_diag%u(jc,jk,blockNo)&
!             !                                     &+p_diag%v(jc,jk,blockNo)*p_diag%v(jc,jk,blockNo))
!           ENDIF
!         END DO
!       END DO
!     END DO
!
!    CALL grad_fd_norm_oce_3d( p_diag%kin, &
!       & patch_2D,    &
!       & grad_coeff, &
!       & z_grad_ekin_rbf)
!
! !     CALL grad_fd_norm_oce( p_diag%kin, &
! !       & patch_2D,    &
! !       & z_grad_ekin_rbf, opt_start_level=start_level,opt_elev=elev)
! !     CALL sync_patch_array(SYNC_C, patch_2D, z_grad_ekin_rbf)
!
!
!     !Add relative vorticity and gradient of kinetic energy to obtain complete horizontal advection
!     DO blockNo = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
!       DO jk = start_level, elev
!         DO je = start_edge_index, end_edge_index
!           IF ( v_base%lsm_e(je,jk,blockNo) <= sea_boundary ) THEN
!             veloc_adv_horz_e(je,jk,blockNo) =&
!               & z_vort_flx_rbf(je,jk,blockNo) + z_grad_ekin_rbf(je,jk,blockNo)
!           ELSE
!             veloc_adv_horz_e(je,jk,blockNo) = 0.0_wp
!           ENDIF
!         END DO
!       END DO
!     END DO
!
!     !---------Debug Diagnostics-------------------------------------------
!     idt_src=2  ! output print level (1-5, fix)
!     CALL dbg_print('HorzRBF: vorticity adv.'     ,veloc_adv_horz_e         ,str_module,idt_src)
!     idt_src=3  ! output print level (1-5, fix)
!     CALL dbg_print('HorzRBF: kin. energy'        ,p_diag%kin               ,str_module,idt_src)
!     CALL dbg_print('HorzRBF: vorticity'          ,p_diag%vort              ,str_module,idt_src)
!     idt_src=4  ! output print level (1-5, fix)
!     CALL dbg_print('HorzRBF: grad. kin. en.'     ,z_grad_ekin_rbf          ,str_module,idt_src)
!     CALL dbg_print('HorzRBF: vorticity_e'        ,z_vort_e                 ,str_module,idt_src)
!     CALL dbg_print('HorzRBF: vorticity flux'     ,z_vort_flx_rbf           ,str_module,idt_src)
!     !---------------------------------------------------------------------
!
!   END SUBROUTINE veloc_adv_horz_rbf
!   !-------------------------------------------------------------------------

!   !-------------------------------------------------------------------------
!   !>
!   !! Computes vertical advection of a (edge based) horizontal vector field.
!   !! The vertical derivative of the velocity vector at circumcenters that
!   !! is reconstructed from edge data is calculated and then multiplied by
!   !! the vertical velocity. The product is mapped from top of the computational
!   !! prism to the middle (still at centers) via the transposed of vertical differentiation
!   !! and then transformed to edges.
!   !!
!   !! IMPORTANT: It is assumed that the velocity vector reconstruction from
!   !! edges to cells has been done before.
!   !!
!   !! input:  lives on cells (velocity points)
!   !! output: lives on edges (velocity points)
!   !!
!   !!  mpi parallelized LL
!    !!
!   SUBROUTINE veloc_adv_vert_rbf( patch_2D, u_c, v_c, w_c, &
!     & top_bc_u_c, top_bc_v_c, &
!     & bot_bc_u_c,  bot_bc_v_c,&
!     & top_bc_w_c,  bot_bc_w_c,&
!     & veloc_adv_vert_e)
!     !
!     !  patch on which computation is performed
!     !
!     TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
!
!     !
!     ! Components of cell based variable which is vertically advected
!     REAL(wp), INTENT(in) :: u_c(:,:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
!     REAL(wp), INTENT(in) :: v_c(:,:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
!     REAL(wp), INTENT(in) :: w_c(:,:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
!     !
!     ! Top boundary condition for cell based variables
!     REAL(wp), INTENT(in) :: top_bc_u_c(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
!     REAL(wp), INTENT(in) :: top_bc_v_c(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
!     !
!     ! Bottom boundary condition for cell based variables
!     REAL(wp), INTENT(in) :: bot_bc_u_c(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
!     REAL(wp), INTENT(in) :: bot_bc_v_c(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
!     !
!     REAL(wp), INTENT(in) :: top_bc_w_c(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
!     REAL(wp), INTENT(in) :: bot_bc_w_c(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
!
!     ! variable in which horizontally advected velocity is stored
!     REAL(wp), INTENT(inout) :: veloc_adv_vert_e(:,:,:)
!
!     INTEGER :: start_level, elev     ! vertical start and end level
!     INTEGER :: jc, jk, blockNo, i_dolic
!     INTEGER :: start_index, end_index
!
!     TYPE(t_subset_range), POINTER :: all_cells
!
!     REAL(wp) :: z_adv_u_i(nproma,n_zlev+1,patch_2D%alloc_cell_blocks),  &
!       & z_adv_v_i(nproma,n_zlev+1,patch_2D%alloc_cell_blocks),  &
!       & z_adv_u_m(nproma,n_zlev,patch_2D%alloc_cell_blocks),  &
!       & z_adv_v_m(nproma,n_zlev,patch_2D%alloc_cell_blocks)
!     !-----------------------------------------------------------------------
!     all_cells => patch_2D%cells%all
!
!     ! #slo# set local variable to zero due to nag -nan compiler-option
!     z_adv_u_i(:,:,:) = 0.0_wp
!     z_adv_v_i(:,:,:) = 0.0_wp
!     z_adv_u_m(:,:,:) = 0.0_wp
!     z_adv_v_m(:,:,:) = 0.0_wp
!
!     start_level = 1
!     elev = n_zlev
!
!     !Step 1: multiply vertical velocity with vertical derivative of horizontal velocity
!     !This requires appropriate boundary conditions
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!       DO jk = start_level, elev
!         DO jc = start_index, end_index
!           !check if we have at least two layers of water
!           !  #slo# - 2011-04-01 - Is this really intended here
!           !  maybe this condition should be fulfilled everywhere
!           !  then it must be calculated in fill_vertical_domain
!           !  this condition could then be omitted here
!           IF (v_base%dolic_c(jc,blockNo) >= 2) THEN
!             !1b) ocean bottom
!             IF ( jk == v_base%dolic_c(jc,blockNo) ) THEN
!               ! u,v-component
!               z_adv_u_i(jc,jk+1,blockNo) = bot_bc_w_c(jc,blockNo)*bot_bc_u_c(jc,blockNo)
!               z_adv_v_i(jc,jk+1,blockNo) = bot_bc_w_c(jc,blockNo)*bot_bc_v_c(jc,blockNo)
!
!               !1c) ocean interior
!             ELSEIF( jk>start_level .AND.  jk < v_base%dolic_c(jc,blockNo))THEN
!               ! u,v-component
!               z_adv_u_i(jc,jk,blockNo)&
!                 & = w_c(jc,jk,blockNo) *( u_c(jc,jk-1,blockNo) - u_c(jc,jk,blockNo) )&
!                 & / v_base%prism_center_distance(jk)
!
!               z_adv_v_i(jc,jk,blockNo)&
!                 & = w_c(jc,jk,blockNo) *( v_c(jc,jk-1,blockNo) - v_c(jc,jk,blockNo) )&
!                 & / v_base%prism_center_distance(jk) !&
!               ! write(*,*)'vert adv:v: ',jk, jc,blockNo,w_c(jc,jk,blockNo) *( v_c(jc,jk,blockNo) - v_c(jc,jk-1,blockNo) ),&
!               ! &  v_c(jc,jk,blockNo), v_c(jc,jk-1,blockNo),v_base%prism_center_distance(jk-1)
!               ! write(*,*)'vert adv:u: ',jk, jc,blockNo,w_c(jc,jk,blockNo) *( u_c(jc,jk,blockNo) - u_c(jc,jk-1,blockNo) ),&
!               ! &  u_c(jc,jk,blockNo), u_c(jc,jk-1,blockNo)
!
!             ENDIF  ! jk-condition
!           ENDIF    ! at least 2 vertical layers
!         END DO
!       END DO
!       !  write(*,*)'A max/min vert adv:',jk, maxval(z_adv_u_i(:,jk,:)), minval(z_adv_u_i(:,jk,:)),&
!       !  & maxval(z_adv_v_i(:,jk,:)), minval(z_adv_v_i(:,jk,:))
!     END DO
!
!     ! Step 2: Map product of vertical velocity & vertical derivative from top of prism to mid position.
!     ! This mapping is the transposed of the vertical differencing.
!
!     !1) From surface down to one layer before bottom
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!       DO jk = start_level, elev-1
!         DO jc = start_index, end_index
!           !check if we are on land: To be replaced by 3D lsm
!           ! #slo# 2011-05-11 - replace by consistent formulation: vertical loop down to dolic
!           IF ( v_base%lsm_c(jc,jk,blockNo) <= sea_boundary ) THEN
!
!             z_adv_u_m(jc,jk,blockNo) &
!               & = (v_base%prism_center_distance(jk+1)*z_adv_u_i(jc,jk+1,blockNo)&
!               & +  v_base%prism_center_distance(jk)*z_adv_u_i(jc,jk,blockNo)) &
!               & / (v_base%prism_center_distance(jk+1)+v_base%prism_center_distance(jk))
!
!             z_adv_v_m(jc,jk,blockNo)&
!               & = (v_base%prism_center_distance(jk+1)*z_adv_v_i(jc,jk+1,blockNo)&
!               & +  v_base%prism_center_distance(jk)*z_adv_v_i(jc,jk,blockNo))&
!               & / (v_base%prism_center_distance(jk+1)+v_base%prism_center_distance(jk))
!           ENDIF
!         END DO
!       END DO
!       ! write(*,*)'B max/min vert adv:',jk, maxval(z_adv_u_m(:,jk,:)), minval(z_adv_u_m(:,jk,:)),&
!       ! & maxval(z_adv_v_m(:,jk,:)), minval(z_adv_v_m(:,jk,:))
!     END DO
!     !Bottom layer
!     !The value of v_base%prism_center_distance at the botom is 0.5*v_base%del_zlev_m
!     !The dimensioning of the firs arrays requires to seperate the vertical loop.
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!       DO jc = start_index, end_index
!         IF ( v_base%dolic_c(jc,blockNo)>0 ) THEN  ! wet points only
!           i_dolic = v_base%dolic_c(jc,blockNo)
!           z_adv_u_m(jc,i_dolic,blockNo) =0.0_wp!&
!           !& = (0.5_wp*v_base%del_zlev_m(elev)*z_adv_u_i(jc,elev+1,blockNo)&
!           !& +        v_base%prism_center_distance(elev)*z_adv_u_i(jc,elev,blockNo)) &
!           !& / (2.0_wp*v_base%del_zlev_m(elev))
!
!           z_adv_v_m(jc,i_dolic,blockNo)=0.0_wp!&
!           !& = (0.5_wp*v_base%del_zlev_m(elev)*z_adv_v_i(jc,elev+1,blockNo)&
!           !&   +        v_base%prism_center_distance(elev)*z_adv_v_i(jc,elev,blockNo))&
!           !& / (2.0_wp*v_base%del_zlev_m(elev))
!
!         END IF
!       END DO
!     END DO
!
!     ! Step 3: Map result of previous calculations from cell centers to edges (for all vertical layers)
! !     CALL primal_map_c2e( patch_2D,&
! !       & z_adv_u_m, z_adv_v_m,&
! !       & veloc_adv_vert_e )
!     ! result is synced in the called funtion
!
!   ! !---------Debug Diagnostics-------------------------------------------
!   !  NOT YET
!   ! idt_src=3  ! output print level (1-5, fix)
!   ! CALL dbg_print('VertRBF: kin. energy'        ,p_diag%kin               ,str_module,idt_src)
!   ! CALL dbg_print('VertRBF: vorticity'          ,p_diag%vort              ,str_module,idt_src)
!   ! idt_src=4  ! output print level (1-5, fix)
!   ! CALL dbg_print('VertRBF: grad. kin. en.'     ,z_grad_ekin_rbf          ,str_module,idt_src)
!   ! CALL dbg_print('VertRBF: vorticity_e'        ,z_vort_e                 ,str_module,idt_src)
!   ! CALL dbg_print('VertRBF: vorticity flux'     ,z_vort_flx_rbf           ,str_module,idt_src)
!   ! CALL dbg_print('VertRBF: vorticity adv.'     ,veloc_adv_horz_e         ,str_module,idt_src)
!   ! !---------------------------------------------------------------------
!
!   END SUBROUTINE veloc_adv_vert_rbf
!   !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! ! FUNCTION laplacian4vortex_flux(patch_2D, vort_flux_in) RESULT(vort_flux_out)
  ! ! TYPE(t_patch),TARGET, INTENT(in) :: patch_2D
  ! ! REAL(wp), INTENT(inout) :: vort_flux_in(:,:,:)
  ! ! REAL(wp) ::vort_flux_out(SIZE(vort_flux_in,1), SIZE(vort_flux_in,2), SIZE(vort_flux_in,3))
  ! !
  ! ! !
  ! ! !Local Variables
  ! ! !
  ! ! INTEGER :: start_level, elev
  ! ! INTEGER :: jk, blockNo, jc, je
  ! ! INTEGER :: i_startblk_c, i_endblk_c, start_index_c, end_index_c
  ! ! INTEGER :: i_startblk_e, i_endblk_e, start_edge_index, end_edge_index
  ! ! INTEGER :: rl_start_e, rl_end_e,rl_start_c, rl_end_c
  ! ! REAL(wp) :: z_tmp(nproma,n_zlev,patch_2D%nblks_e)
  ! ! TYPE(t_cartesian_coordinates)    :: z_pv_cc(nproma,n_zlev,patch_2D%alloc_cell_blocks)
  ! ! INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  ! ! INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
  ! ! TYPE(t_cartesian_coordinates) :: z_grad_u(nproma,n_zlev,patch_2D%nblks_e)
  ! ! TYPE(t_cartesian_coordinates) :: z_div_grad_u(nproma,n_zlev,patch_2D%alloc_cell_blocks)
  ! ! !-----------------------------------------------------------------------
  ! ! rl_start_e = 1
  ! ! rl_end_e   = min_rledge
  ! ! rl_start_c = 1
  ! ! rl_end_c   = min_rlcell
  ! !
  ! ! i_startblk_c = patch_2D%cells%start_blk(rl_start_c,1)
  ! ! i_endblk_c   = patch_2D%cells%end_blk(rl_end_c,1)
  ! ! i_startblk_e = patch_2D%edges%start_blk(rl_start_e,1)
  ! ! i_endblk_e   = patch_2D%edges%end_blk(rl_end_e,1)
  ! !
  ! ! start_level = 1
  ! ! elev = n_zlev
  ! !
  ! ! z_tmp = vort_flux_in
  ! !
  ! ! CALL map_edges2cell( patch_2D, z_tmp, z_pv_cc)
  ! ! DO blockNo = i_startblk_c, i_endblk_c
  ! !
  ! !   CALL get_indices_c(patch_2D, blockNo, i_startblk_c, i_endblk_c, &
  ! !                      start_index_c, end_index_c, rl_start_c, rl_end_c)
  ! !     DO jk = start_level, elev
  ! !       DO jc = start_index_c, end_index_c
  ! !           z_div_grad_u(jc,jk,blockNo)%x =  0.0_wp
  ! !       END DO
  ! !     END DO
  ! !   END DO
  ! ! DO blockNo = i_startblk_e, i_endblk_e
  ! !
  ! !   CALL get_indices_e( patch_2D, blockNo, i_startblk_e, i_endblk_e,&
  ! !                    &  start_edge_index, end_edge_index,&
  ! !                    &  rl_start_e, rl_end_e)
  ! !   DO jk = start_level, elev
  ! !     DO je = start_edge_index, end_edge_index
  ! !       z_grad_u(je,jk,blockNo)%x = 0.0_wp
  ! !     ENDDO
  ! !   END DO
  ! ! END DO
  ! !
  ! ! !Step 1: Calculate gradient of cell velocity vector.
  ! ! !Result is a gradient vector, located at edges
  ! ! !Step 2: Multiply each component of gradient vector with mixing coefficients
  ! ! DO blockNo = i_startblk_e, i_endblk_e
  ! !
  ! !   CALL get_indices_e( patch_2D, blockNo, i_startblk_e, i_endblk_e,&
  ! !                    &  start_edge_index, end_edge_index,&
  ! !                    &  rl_start_e, rl_end_e)
  ! !   DO jk = start_level, elev
  ! !     DO je = start_edge_index, end_edge_index
  ! !
  ! !       !Get indices of two adjacent triangles
  ! !       il_c1 = patch_2D%edges%cell_idx(je,blockNo,1)
  ! !       ib_c1 = patch_2D%edges%cell_blk(je,blockNo,1)
  ! !       il_c2 = patch_2D%edges%cell_idx(je,blockNo,2)
  ! !       ib_c2 = patch_2D%edges%cell_blk(je,blockNo,2)
  ! !
  ! !     IF ( v_base%lsm_e(je,jk,blockNo) <= sea_boundary ) THEN
  ! !       z_grad_u(je,jk,blockNo)%x = &
  ! !         &                  (z_pv_cc(il_c2,jk,ib_c2)%x &
  ! !         &                  - z_pv_cc(il_c1,jk,ib_c1)%x)              &
  ! !         &                  / patch_2D%edges%dual_edge_length(je,blockNo)
  ! !     ELSE
  ! !       z_grad_u(je,jk,blockNo)%x = 0.0_wp
  ! !     ENDIF
  ! !     ENDDO
  ! !   END DO
  ! ! END DO
  ! !
  ! ! !Step 2: Apply divergence to each component of mixing times gradient vector
  ! ! iidx => patch_2D%cells%edge_idx
  ! ! iblk => patch_2D%cells%edge_blk
  ! !
  ! ! DO blockNo = i_startblk_c, i_endblk_c
  ! !
  ! !   CALL get_indices_c(patch_2D, blockNo, i_startblk_c, i_endblk_c, &
  ! !                      start_index_c, end_index_c, rl_start_c, rl_end_c)
  ! !     DO jk = start_level, elev
  ! !       DO jc = start_index_c, end_index_c
  ! !
  ! !          IF ( v_base%lsm_c(jc,jk,blockNo) >= boundary ) THEN
  ! !            z_div_grad_u(jc,jk,blockNo)%x = 0.0_wp
  ! !          ELSE
  ! !           z_div_grad_u(jc,jk,blockNo)%x =  &
  ! !             z_grad_u(iidx(jc,blockNo,1),jk,iblk(jc,blockNo,1))%x * p_int_state(1)%geofac_div(jc,1,blockNo) + &
  ! !             z_grad_u(iidx(jc,blockNo,2),jk,iblk(jc,blockNo,2))%x * p_int_state(1)%geofac_div(jc,2,blockNo) + &
  ! !             z_grad_u(iidx(jc,blockNo,3),jk,iblk(jc,blockNo,3))%x * p_int_state(1)%geofac_div(jc,3,blockNo)
  ! !         ENDIF
  ! !       END DO
  ! !     END DO
  ! !   END DO
  ! !
  ! ! !Step 3: Map divergence back to edges
  ! ! CALL map_cell2edges( patch_2D, z_div_grad_u, z_tmp)
  ! ! vort_flux_out = vort_flux_in + 1000000.0_wp*z_tmp
  ! !
  ! !
  ! ! END FUNCTION laplacian4vortex_flux
  !-------------------------------------------------------------------------

END MODULE mo_ocean_velocity_advection
