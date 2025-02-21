! Contains the implementation of the tracer transport routines for the ICON ocean model
! using the z* vertical co-ordinate
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

MODULE mo_ocean_tracer_zstar
  !-------------------------------------------------------------------------
  USE mo_kind,                         ONLY: wp
  USE mo_math_constants,               ONLY: dbl_eps
  USE mo_ocean_nml,                    ONLY: n_zlev, &
    & l_with_vert_tracer_advection, l_with_vert_tracer_diffusion, &
    & GMRedi_configuration, Cartesian_Mixing, &
    & vert_mix_type
  USE mo_parallel_config,              ONLY: nproma
  USE mo_run_config,                   ONLY: dtime
  USE mo_model_domain,                 ONLY: t_patch, t_patch_3d
  USE mo_exception,                    ONLY: finish
  USE mo_ocean_math_operators,         ONLY: div_oce_3d
  USE mo_operator_ocean_coeff_3d,      ONLY: t_operator_coeff
  USE mo_grid_subset,                  ONLY: t_subset_range, get_index_range
  USE mo_sync,                         ONLY: sync_c, sync_c1, sync_patch_array, sync_patch_array_mult
  USE mo_ocean_tracer_transport_types, ONLY: t_ocean_transport_state, t_ocean_tracer, t_tracer_collection
  USE mo_scalar_product,               ONLY: map_edges2edges_sc_zstar
  USE mo_ocean_tracer_transport_vert,  ONLY: advect_flux_vertical
  USE mo_util_dbg_prnt,                ONLY: dbg_print
  USE mo_ocean_tracer_transport_horz,  ONLY: advect_horz, diffuse_horz
  USE mo_ocean_tracer_diffusion,       ONLY: tracer_diffusion_vertical_implicit
  USE mo_fortran_tools,                ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=15) :: str_module = 'oceTracer_zstar'  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug

  LOGICAL :: eliminate_upper_diag = .true.

  PUBLIC :: advect_ocean_tracers_zstar
  PUBLIC :: advect_individual_tracers_zstar
  PUBLIC :: upwind_zstar_hflux_oce
  PUBLIC :: limiter_ocean_zalesak_horz_zstar
  PUBLIC :: tracer_diffusion_vertical_implicit_zstar
  PUBLIC :: eliminate_upper_diag

  INTERFACE limiter_ocean_zalesak_horz_zstar



    MODULE PROCEDURE limiter_ocean_zalesak_horz_zstar_scalar

  END INTERFACE

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
  !!  Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !!
  !! Adapted for zstar
  !! FIXME: The limiter assumes no knowledge of eta for next time step which
  !! would be required if the formulation were to be correct
  !!
  SUBROUTINE limiter_ocean_zalesak_horz_zstar_vector( patch_3d,&
    & vert_velocity,          &
    & tracer,                 &
    & p_mass_flx_e,           &
    & flx_tracer_low,         &
    & flx_tracer_high,        &
    & div_adv_flux_vert,      &
    & stretch_c,              &
    & operators_coefficients, &
    & flx_tracer_final,       &
    & lacc )

    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp),INTENT(inout)              :: vert_velocity(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !
    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
    LOGICAL, INTENT(in), OPTIONAL       :: lacc

    !Local variables
    REAL(wp) :: z_mflx_anti
    REAL(wp) :: z_fluxdiv_c(nproma)     !< flux divergence at cell center
    REAL(wp) :: z_anti          (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)
    REAL(wp) :: z_tracer_new_low(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< new tracer field after transport, if low order fluxes are used
    REAL(wp) :: z_tracer_max    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local maximum of current tracer value and low order update
    REAL(wp) :: z_tracer_min    (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< local minimum of current tracer value and low order update
    REAL(wp) :: r_p             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< fraction which must multiply all in/out fluxes of cell jc to guarantee
    REAL(wp) :: r_m             (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)!< no overshoot/undershoot
    REAL(wp) :: r_frac          !< computed minimum fraction which must multiply< the flux at the edge
    REAL(wp) :: z_min(nproma), z_max(nproma)    !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_signum        !< sign of antidiffusive velocity
    REAL(wp) :: p_p(nproma), p_m(nproma)        !< sum of antidiffusive fluxes into and out of cell jc
    REAL(wp) :: inv_prism_thick_new
    REAL(wp) :: delta_z_new, delta_z
    INTEGER, DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  cellOfEdge_idx, cellOfEdge_blk
    INTEGER, DIMENSION(:,:,:), POINTER, CONTIGUOUS :: neighbor_cell_idx, neighbor_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER, CONTIGUOUS :: edge_of_cell_idx, edge_of_cell_blk
    INTEGER :: start_level, end_level, max_level
    INTEGER :: start_index, end_index
    INTEGER :: edge_index, level, blockNo, jc,  cell_connect, max_edges
    LOGICAL :: lzacc
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d

    ! Pointers needed for GPU/OpenACC
    INTEGER, POINTER :: dolic_c(:,:), dolic_e(:,:), cells_num_edges(:,:), edges_SeaBoundaryLevel(:,:,:)
    REAL(wp), POINTER :: inv_prism_thick_c(:,:,:), prism_thick_flat_sfc_c(:,:,:), div_coeff(:,:,:,:)
    !-------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    !-------------------------------------------------------------------------
    start_level = 1
    end_level   = n_zlev
    cellOfEdge_idx  => patch_2d%edges%cell_idx
    cellOfEdge_blk  => patch_2d%edges%cell_blk
    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
    neighbor_cell_idx => patch_2d%cells%neighbor_idx
    neighbor_cell_blk => patch_2d%cells%neighbor_blk

    dolic_c => patch_3d%p_patch_1d(1)%dolic_c
    dolic_e => patch_3d%p_patch_1d(1)%dolic_e
    inv_prism_thick_c => patch_3d%p_patch_1d(1)%inv_prism_thick_c
    prism_thick_flat_sfc_c => patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c
    cells_num_edges => patch_2d%cells%num_edges
    edges_SeaBoundaryLevel => operators_coefficients%edges_SeaBoundaryLevel
    div_coeff => operators_coefficients%div_coeff

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(patch_3d%p_patch_2d(1)%alloc_cell_blocks, patch_3d%p_patch_2d(1)%nblks_e) &
    !$ACC   CREATE(z_fluxdiv_c, z_anti, z_tracer_new_low, z_tracer_max, z_tracer_min) &
    !$ACC   CREATE(r_p, r_m, z_min, z_max, p_p, p_m) IF(lzacc)


    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL

!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level, max_level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)

      max_level = MIN(MAXVAL(dolic_e(start_index:end_index,blockNo)), end_level)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      z_anti(:,:,blockNo)     = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO level = start_level, max_level
        DO edge_index = start_index, end_index
          IF (level <= dolic_e(edge_index,blockNo)) THEN
            ! calculate antidiffusive flux for each edge
            z_anti(edge_index,level,blockNo) = flx_tracer_high(edge_index,level,blockNo)&
                                            &- flx_tracer_low(edge_index,level,blockNo)
          END IF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
      !$ACC END PARALLEL LOOP
    END DO  ! end loop over blocks
    !$ACC WAIT(1)
!ICON_OMP_END_DO


!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, delta_z, delta_z_new, max_level, max_edges, &
!ICON_OMP z_fluxdiv_c) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      z_tracer_new_low(:,:,blockNo)    = 0.0_wp
      z_tracer_max(:,:,blockNo)        = 0.0_wp
      z_tracer_min(:,:,blockNo)        = 0.0_wp
      !$ACC END KERNELS

      max_level = MIN(MAXVAL(dolic_c(start_index:end_index,blockNo)), end_level)
      max_edges = MAXVAL(cells_num_edges(start_index:end_index,blockNo))


      ! 3. Compute the complete (with horizontal and vertical divergence) updated low order solution z_tracer_new_low
      !  compute divergence of low order fluxes
      DO level = start_level, max_level
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        z_fluxdiv_c(:) = 0._wp
        !$ACC END KERNELS

        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO cell_connect = 1, max_edges
          DO jc = start_index, end_index
            IF (level <= dolic_c(jc,blockNo) .AND. cell_connect <= cells_num_edges(jc,blockNo)) THEN
              z_fluxdiv_c(jc) =  z_fluxdiv_c(jc) + &
                & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
                & div_coeff(jc,level,blockNo,cell_connect)
            ENDIF
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(delta_z, delta_z_new) ASYNC(1) IF(lzacc)
        DO jc = start_index, end_index
          IF (level <= dolic_c(jc,blockNo)) THEN
            delta_z     = stretch_c(jc, blockNo)*prism_thick_flat_sfc_c(jc,level,blockNo)
            delta_z_new = stretch_c(jc, blockNo)*prism_thick_flat_sfc_c(jc,level,blockNo)
            !

            z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
            & - dtime * (z_fluxdiv_c(jc)+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new
          ENDIF
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO

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

    ENDDO
    !$ACC WAIT(1)
!ICON_OMP_END_DO

!ICON_OMP_MASTER
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, z_tracer_max, z_tracer_min)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER
    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, inv_prism_thick_new, max_level, max_edges, &
!ICON_OMP z_mflx_anti, z_max, z_min, cell_connect, p_p, p_m) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)

      ! this is only needed for the parallel test setups
      ! it will try  tocheck the uninitialized (land) parts
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      r_m(:,:,blockNo) = 0.0_wp
      r_p(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      max_level = MIN(MAXVAL(dolic_c(start_index:end_index,blockNo)), end_level)
      max_edges = MAXVAL(cells_num_edges(start_index:end_index,blockNo))

      ! 2. Define "antidiffusive" fluxes A(jc,level,blockNo,edge_index) for each cell. It is the difference
      !    between the high order fluxes (given by the FFSL-scheme) and the low order
      !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
      !    - positive for outgoing fluxes
      !    - negative for incoming fluxes
      !    this sign convention is related to the definition of the divergence operator.
      DO level = start_level, max_level
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        p_p(:) = 0.0_wp
        p_m(:) = 0.0_wp

        z_max(:) = z_tracer_max(:,level,blockNo)
        z_min(:) = z_tracer_min(:,level,blockNo)
        !$ACC END KERNELS

        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) &
        !$ACC   PRIVATE(inv_prism_thick_new, p_p, p_m, z_max, z_min) ASYNC(1) IF(lzacc)
        DO cell_connect = 1, max_edges
          DO jc = start_index, end_index
            IF (level <= dolic_c(jc,blockNo) .AND. &
                cell_connect <= cells_num_edges(jc,blockNo) .AND. &
                level <= dolic_c( &
                    neighbor_cell_idx(jc,blockNo,cell_connect), &
                    neighbor_cell_blk(jc,blockNo,cell_connect))) THEN

              inv_prism_thick_new = inv_prism_thick_c(jc,level,blockNo) &
                & / stretch_c(jc, blockNo)

              z_max(jc) = MAX(z_max(jc), &
                & z_tracer_max(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))
              z_min(jc) = MIN(z_min(jc), &
                & z_tracer_min(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))

              z_mflx_anti = &
                & dtime * div_coeff(jc,level,blockNo,cell_connect) * inv_prism_thick_new &
                & * z_anti(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect))

              ! Sum of all incoming antidiffusive fluxes into cell jc
              ! outgoing fluxes carry a positive sign, incoming a negative
              p_p(jc) = p_p(jc) - MIN(0._wp, z_mflx_anti)
              ! Sum of all outgoing antidiffusive fluxes out of cell jc
              p_m(jc) = p_m(jc) + MAX(0._wp, z_mflx_anti)

            ENDIF
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = start_index, end_index
          IF (level <= dolic_c(jc,blockNo)) THEN
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
          ENDIF
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
    ENDDO
    !$ACC WAIT(1)
!ICON_OMP_END_DO


!ICON_OMP_MASTER
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER

    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level, z_signum, r_frac, max_level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)

      max_level = MIN(MAXVAL(dolic_e(start_index:end_index,blockNo)), end_level)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      flx_tracer_final(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO level = start_level, max_level
        DO edge_index = start_index, end_index
          IF (level <= dolic_e(edge_index,blockNo)) THEN
            IF (edges_SeaBoundaryLevel(edge_index,level,blockNo) > -2)THEN! edge < 2nd order boundary

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
              flx_tracer_final(edge_index,level,blockNo) = flx_tracer_low(edge_index,level,blockNo) &
                & + MIN(1.0_wp,r_frac) *z_anti(edge_index,level,blockNo)

            ENDIF
          ENDIF
        ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO
    !$ACC WAIT(1)
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

    !$ACC END DATA
  END SUBROUTINE limiter_ocean_zalesak_horz_zstar_vector
  !-------------------------------------------------------------------------


  SUBROUTINE limiter_ocean_zalesak_horz_zstar_scalar( patch_3d,&
    & vert_velocity,          &
    & tracer,                 &
    & p_mass_flx_e,           &
    & flx_tracer_low,         &
    & flx_tracer_high,        &
    & div_adv_flux_vert,      &
    & stretch_c,              &
    & operators_coefficients, &
    & flx_tracer_final,       &
    & lacc )

    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    REAL(wp),INTENT(inout)              :: vert_velocity(nproma,n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: tracer           (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)             :: p_mass_flx_e     (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_low   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_high  (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: flx_tracer_final (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(inout)             :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !
    TYPE(t_operator_coeff),INTENT(in)   :: operators_coefficients
    LOGICAL, INTENT(in), OPTIONAL       :: lacc

    !Local variables
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
    REAL(wp) :: inv_prism_thick_new
    REAL(wp) :: delta_z_new, delta_z
    INTEGER, DIMENSION(:,:,:), POINTER ::  cellOfEdge_idx, cellOfEdge_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    INTEGER :: start_level, end_level
    INTEGER :: start_index, end_index
    INTEGER :: edge_index, level, blockNo, jc,  cell_connect, sum_lsm_quad_edge
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc

    INTEGER  :: bt_lev
    TYPE(t_subset_range), POINTER :: all_cells

    ! Pointers needed for GPU/OpenACC
    INTEGER, POINTER :: dolic_c(:,:), dolic_e(:,:), cells_num_edges(:,:), edges_SeaBoundaryLevel(:,:,:)
    REAL(wp), POINTER :: inv_prism_thick_c(:,:,:), div_coeff(:,:,:,:)
    !-------------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    cells_in_domain => patch_2d%cells%in_domain
    all_cells       => patch_2d%cells%ALL
    !-------------------------------------------------------------------------
    start_level = 1
    end_level   = n_zlev
    cellOfEdge_idx  => patch_2d%edges%cell_idx
    cellOfEdge_blk  => patch_2d%edges%cell_blk
    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk
    neighbor_cell_idx => patch_2d%cells%neighbor_idx
    neighbor_cell_blk => patch_2d%cells%neighbor_blk

    dolic_c => patch_3d%p_patch_1d(1)%dolic_c
    dolic_e => patch_3d%p_patch_1d(1)%dolic_e
    inv_prism_thick_c => patch_3d%p_patch_1d(1)%inv_prism_thick_c
    cells_num_edges => patch_2d%cells%num_edges
    edges_SeaBoundaryLevel => operators_coefficients%edges_SeaBoundaryLevel
    div_coeff => operators_coefficients%div_coeff

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(patch_3d%p_patch_2d(1)%alloc_cell_blocks, patch_3d%p_patch_2d(1)%nblks_e) &
    !$ACC   PRESENT(patch_3d%p_patch_2d(1)%cells%max_connectivity) &
    !$ACC   CREATE(z_mflx_anti, z_anti, z_tracer_new_low, z_tracer_max, z_tracer_min) &
    !$ACC   CREATE(r_p, r_m, z_tracer_update_horz) IF(lzacc)


    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL

!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      z_anti(:,:,blockNo)     = 0.0_wp
      !$ACC END KERNELS

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
    END DO  ! end loop over blocks
    !$ACC WAIT(1)
!ICON_OMP_END_DO


!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, delta_z, delta_z_new, &
!ICON_OMP z_fluxdiv_c) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      z_tracer_new_low(:,:,blockNo)    = 0.0_wp
      z_tracer_update_horz(:,:,blockNo)= 0.0_wp
      z_tracer_max(:,:,blockNo)        = 0.0_wp
      z_tracer_min(:,:,blockNo)        = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(delta_z, delta_z_new, z_fluxdiv_c)
      DO jc = start_index, end_index
        IF (dolic_c(jc,blockNo) < 1) CYCLE

        ! 3. Compute the complete (with horizontal and vertical divergence) updated low order solution z_tracer_new_low
        DO level = start_level  , MIN(dolic_c(jc,blockNo), end_level)
          !  compute divergence of low order fluxes
          z_fluxdiv_c = 0
          DO cell_connect = 1, cells_num_edges(jc,blockNo)
            z_fluxdiv_c =  z_fluxdiv_c + &
              & flx_tracer_low(edge_of_cell_idx(jc,blockNo,cell_connect),level,edge_of_cell_blk(jc,blockNo,cell_connect)) * &
              & div_coeff(jc,level,blockNo,cell_connect)
          ENDDO

          delta_z     = stretch_c(jc, blockNo)*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
          delta_z_new = stretch_c(jc, blockNo)*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,blockNo)
          !

          z_tracer_new_low(jc,level,blockNo) = (tracer(jc,level,blockNo) * delta_z                     &
            & - dtime * (z_fluxdiv_c+div_adv_flux_vert(jc,level,blockNo)))/delta_z_new

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
    ENDDO
    !$ACC WAIT(1)
!ICON_OMP_END_DO

!ICON_OMP_MASTER
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, z_tracer_max, z_tracer_min)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER
    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, level, inv_prism_thick_new, &
!ICON_OMP z_mflx_anti, z_max, z_min, cell_connect, p_p, p_m) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

      ! this is only needed for the parallel test setups
      ! it will try  tocheck the uninitialized (land) parts
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      r_m(:,:,blockNo) = 0.0_wp
      r_p(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(inv_prism_thick_new, p_p, p_m, z_max, z_min)
      DO jc = start_index, end_index

        ! get prism thickness
        DO level = start_level, MIN(dolic_c(jc,blockNo), end_level)
          inv_prism_thick_new = inv_prism_thick_c(jc,level,blockNo) &
            & /stretch_c(jc, blockNo)

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
          DO cell_connect = 1, cells_num_edges(jc,blockNo)
            IF (dolic_c(neighbor_cell_idx(jc,blockNo,cell_connect), neighbor_cell_blk(jc,blockNo,cell_connect)) >= level) THEN

              z_max = MAX(z_max, &
                & z_tracer_max(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))
              z_min = MIN(z_min, &
                & z_tracer_min(neighbor_cell_idx(jc,blockNo,cell_connect),level,neighbor_cell_blk(jc,blockNo,cell_connect)))

              z_mflx_anti(cell_connect) =                                                        &
                & dtime * div_coeff(jc,level,blockNo,cell_connect) * inv_prism_thick_new  &
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
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
    !$ACC WAIT(1)
!ICON_OMP_END_DO


!ICON_OMP_MASTER
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(sync_c1, patch_2d, 2, r_m, r_p)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER

    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level, z_signum, r_frac) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      flx_tracer_final(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(z_signum, r_frac)
      DO edge_index = start_index, end_index

        DO level = start_level, MIN(dolic_e(edge_index,blockNo), end_level)

          IF (edges_SeaBoundaryLevel(edge_index,level,blockNo) > -2)THEN! edge < 2nd order boundary

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
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

  !$ACC END DATA
  END SUBROUTINE limiter_ocean_zalesak_horz_zstar_scalar
  !-------------------------------------------------------------------------





  SUBROUTINE upwind_zstar_hflux_oce( patch_3d, cell_value, edge_vn, edge_upwind_flux, opt_start_level, opt_end_level, lacc )

    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)              :: cell_value   (nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)      !< advected cell centered variable
    REAL(wp), INTENT(in)              :: edge_vn    (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)       !< normal velocity on edges
    REAL(wp), INTENT(inout)           :: edge_upwind_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level    ! optional vertical end level
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER, DIMENSION(:,:,:), POINTER :: idx, blk
    INTEGER  :: start_level, end_level
    INTEGER  :: start_index, end_index
    INTEGER  :: edge_index, level, blockNo         !< index of edge, vert level, block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc
    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    idx             => patch_3D%p_patch_2D(1)%edges%cell_idx
    blk             => patch_3D%p_patch_2D(1)%edges%cell_blk
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_start_level) ) THEN
      start_level = opt_start_level
    ELSE
      start_level = 1
    END IF
    IF ( PRESENT(opt_end_level) ) THEN
      end_level = opt_end_level
    ELSE
      end_level = n_zlev
    END IF

    CALL set_acc_host_or_device(lzacc, lacc)
    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
!ICON_OMP_PARALLEL PRIVATE(iilc, iibc)
    ! line and block indices of two neighboring cells
    iilc => patch_2d%edges%cell_idx
    iibc => patch_2d%edges%cell_blk

    ! loop through all patch edges (and blocks)
!ICON_OMP_DO PRIVATE(start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      edge_upwind_flux(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO edge_index = start_index, end_index
        DO level = start_level, MIN(patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo), end_level)
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          edge_upwind_flux(edge_index,level,blockNo) =  &
             0.5_wp * (        edge_vn(edge_index,level,blockNo)  *           &
               & ( cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)) + &
               &   cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2)) ) &
               &   - ABS( edge_vn(edge_index,level,blockNo) ) *               &
               & ( cell_value(iilc(edge_index,blockNo,2),level,iibc(edge_index,blockNo,2)) - &
               &   cell_value(iilc(edge_index,blockNo,1),level,iibc(edge_index,blockNo,1)) ) )

        END DO  ! end loop over edges
      END DO  ! end loop over levels
      !$ACC END PARALLEL
    END DO  ! end loop over blocks
    !$ACC WAIT(1)
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

  END SUBROUTINE upwind_zstar_hflux_oce
  !-----------------------------------------------------------------------


   !------------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_vertical_implicit_zstar( &
    & patch_3d,                  &
    & ocean_tracer,              &
    & a_v,     &
    & stretch_c, &
    & lacc )

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: ocean_tracer
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)
    REAL(wp), INTENT(in)                 :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) !
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    !
    INTEGER :: cell_block, start_index, end_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc

    !-----------------------------------------------------------------------
    cells_in_domain       =>  patch_3d%p_patch_2d(1)%cells%in_domain
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO cell_block = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, cell_block, start_index, end_index)

      CALL tracer_diffusion_vertical_implicit_zstar_onBlock( &
        & patch_3d,                  &
        & ocean_tracer,              &
        & a_v, stretch_c,            &
        & cell_block, start_index, end_index, &
        & lacc=lzacc)

    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE tracer_diffusion_vertical_implicit_zstar
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for scalar fields.
  !>
  !! The result ocean_tracer%concetration is calculated on domain_cells
  !-------------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_vertical_implicit_zstar_onBlock( &
    & patch_3d,                &
    & ocean_tracer,            &
    & a_v, stretch_c,          &
    & blockNo, start_index, end_index, &
    & lacc ) !,  &
    ! & diff_column)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: ocean_tracer
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)
    INTEGER, INTENT(in)                  :: blockNo, start_index, end_index
    REAL(wp), INTENT(in)                 :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL        :: lacc
    !
    !
    REAL(wp) :: inv_prism_thickness(nproma,1:n_zlev), inv_prisms_center_distance(nproma,1:n_zlev)
    REAL(wp) :: a(nproma,1:n_zlev), b(nproma,1:n_zlev), c(nproma,1:n_zlev)! , nb(1:n_zlev)
    REAL(wp) :: fact(nproma,1:n_zlev)
    REAL(wp) :: column_tracer(nproma,1:n_zlev)
    REAL(wp) :: dt_inv, diagonal_product
    REAL(wp), POINTER :: field_column(:,:,:)
    INTEGER :: cell_index, level, maxcell
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc

    REAL(wp) :: inv_str_c

    ! Pointers needed for GPU/OpenACC
    INTEGER, POINTER :: dolic_c(:,:)
    REAL(wp), POINTER :: inv_prism_thick_c(:,:,:), inv_prism_center_dist_c(:,:,:)

    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2d%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    dolic_c => patch_3d%p_patch_1d(1)%dolic_c
    inv_prism_thick_c => patch_3d%p_patch_1d(1)%inv_prism_thick_c
    inv_prism_center_dist_c => patch_3d%p_patch_1d(1)%inv_prism_center_dist_c
    !-----------------------------------------------------------------------
    dt_inv = 1.0_wp/dtime

    CALL set_acc_host_or_device(lzacc, lacc)



    !$ACC DATA PRESENT(patch_3d%p_patch_2d(1)%alloc_cell_blocks, patch_3d%p_patch_2d(1)%nblks_e) &
    !$ACC   CREATE(inv_prism_thickness, inv_prisms_center_distance, a, b, c, fact, column_tracer) &
    !$ACC   COPYIN(a_v, stretch_c, dolic_c, inv_prism_thick_c, inv_prism_center_dist_c) &
    !$ACC   COPY(field_column) IF(lzacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR PRIVATE(inv_str_c)
    DO cell_index = start_index, end_index
      IF (dolic_c(cell_index,blockNo) < 2) CYCLE ! nothing to diffuse

      inv_str_c    = 1._wp/stretch_c(cell_index, blockNo)

      DO level=1,dolic_c(cell_index,blockNo)
        inv_prism_thickness(cell_index,level)        = inv_str_c*inv_prism_thick_c(cell_index,level,blockNo)
        inv_prisms_center_distance(cell_index,level) = inv_str_c*inv_prism_center_dist_c(cell_index,level,blockNo)

        column_tracer(cell_index,level) = field_column(cell_index,level,blockNo)
      END DO
    END DO
    !$ACC END PARALLEL

    !------------------------------------
    ! Fill triangular matrix
    ! b is diagonal, a is the upper diagonal, c is the lower
    ! top level
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO cell_index = start_index, end_index
      a(cell_index,1) = 0.0_wp
      c(cell_index,1) = -a_v(cell_index,2,blockNo) * inv_prism_thickness(cell_index,1) &
        & * inv_prisms_center_distance(cell_index,2)*dtime
    END DO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO cell_index = start_index, end_index
      b(cell_index,1) = 1.0_wp - c(cell_index,1)
    END DO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO cell_index = start_index, end_index
      !$ACC LOOP SEQ
      DO level = 2, dolic_c(cell_index,blockNo)-1
        a(cell_index,level) = - a_v(cell_index,level,blockNo) * inv_prism_thickness(cell_index,level) &
          & * inv_prisms_center_distance(cell_index,level)*dtime
        c(cell_index,level) = - a_v(cell_index,level+1,blockNo) * inv_prism_thickness(cell_index,level) &
          & * inv_prisms_center_distance(cell_index,level+1)*dtime
        b(cell_index,level) = 1.0_wp - a(cell_index,level) - c(cell_index,level)
      END DO
    END DO

    ! bottom
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO cell_index = start_index, end_index
      level = dolic_c(cell_index,blockNo)
      IF (dolic_c(cell_index,blockNo) < 2) CYCLE ! nothing to diffuse

      a(cell_index,level) = - a_v(cell_index,level,blockNo) * inv_prism_thickness(cell_index,level) &
        & * inv_prisms_center_distance(cell_index,level)*dtime
      c(cell_index,level) = 0.0_wp
      b(cell_index,level) = 1.0_wp - a(cell_index,level)
    END DO

    IF (eliminate_upper_diag) THEN
      ! solve the tridiagonal matrix by eliminating c (the upper diagonal)

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index = start_index, end_index
        !$ACC LOOP SEQ
        DO level=dolic_c(cell_index,blockNo)-1,1,-1
          IF (dolic_c(cell_index,blockNo) < 2) CYCLE ! nothing to diffuse

          fact(cell_index,level) = c(cell_index,level)/b(cell_index,level+1)
          b(cell_index,level) = b(cell_index,level)-a(cell_index,level+1)*fact(cell_index,level)
          c(cell_index,level) = 0.0_wp
          column_tracer(cell_index,level) = column_tracer(cell_index,level) - &
            & fact(cell_index,level)*column_tracer(cell_index,level+1)
        END DO
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index = start_index, end_index
        IF (dolic_c(cell_index,blockNo) < 2) CYCLE ! nothing to diffuse

        field_column(cell_index,1,blockNo) = column_tracer(cell_index,1)/b(cell_index,1)
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index = start_index, end_index
        !$ACC LOOP SEQ
        DO level=2,dolic_c(cell_index,blockNo)
          IF (dolic_c(cell_index,blockNo) < 2) CYCLE ! nothing to diffuse

          field_column(cell_index,level,blockNo) = (column_tracer(cell_index,level) - &
            & a(cell_index,level)* field_column(cell_index,level-1,blockNo)) / b(cell_index,level)
        END DO
      END DO
    ELSE
      ! solve the tridiagonal matrix by eliminating a (the lower diagonal)

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index = start_index, end_index
        !$ACC LOOP SEQ
        DO level=2, dolic_c(cell_index,blockNo)
          IF (dolic_c(cell_index,blockNo) < 2) CYCLE ! nothing to diffuse
          fact(cell_index,level) = a(cell_index,level)/b(cell_index,level-1)
          b(cell_index,level) = b(cell_index,level)-c(cell_index,level-1)*fact(cell_index,level)
          a(cell_index,level) = 0.0_wp
          column_tracer(cell_index,level) = column_tracer(cell_index,level) - &
            & fact(cell_index,level)*column_tracer(cell_index,level-1)
        END DO
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index = start_index, end_index
        IF (dolic_c(cell_index,blockNo) < 2) CYCLE ! nothing to diffuse

        level = dolic_c(cell_index,blockNo)
        field_column(cell_index,level,blockNo) = column_tracer(cell_index,level)/b(cell_index,level)
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index = start_index, end_index
        !$ACC LOOP SEQ
        DO level=dolic_c(cell_index,blockNo)-1,1,-1
          IF (dolic_c(cell_index,blockNo) < 2) CYCLE ! nothing to diffuse

          field_column(cell_index,level,blockNo) = (column_tracer(cell_index,level) - &
            & c(cell_index,level)* field_column(cell_index,level+1,blockNo)) / b(cell_index,level)
        END DO
      END DO
    END IF
    !$ACC END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE tracer_diffusion_vertical_implicit_zstar_onBlock
  !------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !!Subroutine for all advection operations on a single tracer
  !>
  !-------------------------------------------------------------------------
  SUBROUTINE advect_individual_tracers_zstar( patch_3d, transport_state, &
    & operators_coefficients, stretch_e, stretch_c, stretch_c_new, old_tracer, new_tracer, lacc)


    TYPE(t_patch_3d), POINTER, INTENT(in)                :: patch_3d
    TYPE(t_ocean_transport_state), TARGET                :: transport_state
    TYPE(t_operator_coeff),   INTENT(in)                 :: operators_coefficients
    REAL(wp), INTENT(IN)               :: stretch_e(nproma, patch_3d%p_patch_2d(1)%nblks_e) !! stretch factor
    REAL(wp), INTENT(IN)               :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(IN)               :: stretch_c_new(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_ocean_tracer), TARGET       :: old_tracer
    TYPE(t_ocean_tracer)               :: new_tracer
    LOGICAL, INTENT(in), OPTIONAL      :: lacc

    TYPE(t_patch), POINTER :: patch_2d

    INTEGER  :: jb, jc, je, level
    REAL(wp) :: delta_t, delta_z,delta_z_new
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: top_bc(nproma)
    INTEGER  :: start_cell_index, end_cell_index
    REAL(wp) :: z_adv_flux_h (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_low (nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: z_adv_high(nproma, n_zlev, patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) :: temp(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL :: lzacc

    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    REAL(wp), POINTER :: old_tracer_concentration(:,:,:), new_tracer_concentration(:,:,:)

    CHARACTER(len=*), PARAMETER :: method_name = 'mo_ocean_tracer:advect_diffuse_tracer_zstar'

    CALL set_acc_host_or_device(lzacc, lacc)

    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain
    delta_t = dtime
    old_tracer_concentration => old_tracer%concentration
    new_tracer_concentration => new_tracer%concentration
    !---------------------------------------------------------------------

    !$ACC DATA COPYIN(old_tracer_concentration) &
    !$ACC   COPY(new_tracer_concentration) &
    !$ACC   CREATE(div_adv_flux_horz, div_adv_flux_vert, div_diff_flux_horz, top_bc) &
    !$ACC   CREATE(z_adv_flux_h, z_adv_low, z_adv_high) IF(lzacc)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('on entry: IndTrac: trac_old',old_tracer_concentration(:,:,:), &
      & str_module,idt_src, in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------

    ! these are probably not necessary
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    div_adv_flux_vert(:,:,:)  = 0.0_wp
    div_adv_flux_horz(:,:,:)  = 0.0_wp
    div_diff_flux_horz(:,:,:) = 0.0_wp
    !$ACC END KERNELS
    !$ACC WAIT(1)
    !---------------------------------------------------------------------
    IF ( l_with_vert_tracer_advection ) THEN

      CALL advect_flux_vertical( patch_3d,&
        & old_tracer_concentration, &
        & transport_state,                           &
        & operators_coefficients,                     &
        & div_adv_flux_vert,                          &
        & lacc=lzacc )

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('aft. AdvFluxVert:divfluxvert',div_adv_flux_vert          ,str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    ENDIF  ! l_with_vert_tracer_advection

    !---------------------------------------------------------------------
    !-Horizontal  advection
    !---------------------------------------------------------------------
    CALL upwind_zstar_hflux_oce( patch_3d,  &
      & old_tracer_concentration, &
      & transport_state%mass_flux_e,         &
      & z_adv_flux_h, lacc=lzacc )

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    z_adv_low = z_adv_flux_h
    !$ACC END KERNELS
    !$ACC WAIT(1)

    call map_edges2edges_sc_zstar( patch_3d, transport_state%vn, old_tracer_concentration, &
      & operators_coefficients, stretch_e, z_adv_flux_h, lacc=lzacc )

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    z_adv_high = z_adv_flux_h
    !$ACC END KERNELS
    !$ACC WAIT(1)

    CALL limiter_ocean_zalesak_horz_zstar( patch_3d,   &
      & transport_state%w,           &
      & old_tracer_concentration,              &
      & transport_state%mass_flux_e,           &
      & z_adv_low,                             &
      & z_adv_high,                            &
      & div_adv_flux_vert,                     &
      & stretch_c,                             &
      & operators_coefficients,                &
      & z_adv_flux_h,                          &
      & lacc=lzacc )

    !Calculate divergence of advective fluxes
    CALL div_oce_3d( z_adv_flux_h, patch_3D, operators_coefficients%div_coeff, &
      & div_adv_flux_horz, subset_range=cells_in_domain, lacc=lzacc )
    !---------------------------------------------------------------------

    !! horizontal diffusion, vertical is handled implicitely below
    !! Note that horizontal here implies constant z* lines
    !! Note that diffuse horz takes h_old and h_new but does not use them
    !! So we are passing a temp variable
    IF(GMRedi_configuration==Cartesian_Mixing)THEN
      !horizontal diffusion, vertical is handled implicitely below
      CALL diffuse_horz( patch_3d,         &
      & old_tracer_concentration,          &
      & transport_state,                   &
      & operators_coefficients,            &
      & old_tracer%hor_diffusion_coeff,    &
      & temp,                              &
      & temp,                              &
      & div_diff_flux_horz,                &
      & lacc=lzacc )
    ELSE
      CALL finish(method_name, "wrong GMredi call")
    ENDIF


    !Calculate preliminary tracer value out of horizontal advective and
    !diffusive fluxes and vertical advective fluxes, plus surface forcing.
    !Surface forcing applied as volume forcing at rhs, i.e.part of explicit term
    !in tracer (and also momentum) eqs. In this case, top boundary condition of
    !vertical Laplacians are homogeneous
    !ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, level, &
    !ICON_OMP delta_z, delta_z_new, top_bc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      IF (ASSOCIATED(old_tracer%top_bc)) THEN
        !$ACC DATA COPYIN(old_tracer%top_bc) IF(lzacc)
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        top_bc(:) = old_tracer%top_bc(:,jb)
        !$ACC END KERNELS
        !$ACC WAIT(1)
        !$ACC END DATA
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        top_bc(:) = 0.0_wp
        !$ACC END KERNELS
        !$ACC WAIT(1)
      ENDIF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(delta_z, delta_z_new)
      DO jc = start_cell_index, end_cell_index
        !! d_z*(coeff*w*C) = coeff*d_z(w*C) since coeff is constant for each column
        div_adv_flux_vert(jc, :, jb) = stretch_c(jc, jb)*div_adv_flux_vert(jc, :, jb)

        !! Apply boundary conditions
        DO level = 1, MIN(patch_3d%p_patch_1d(1)%dolic_c(jc,jb),1)  ! this at most should be 1
          delta_z     = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)*stretch_c(jc, jb)
          delta_z_new = patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,level,jb)*stretch_c_new(jc, jb)

          new_tracer_concentration(jc,level,jb)= &
            & (old_tracer_concentration(jc,level,jb) * delta_z &
            & - delta_t * (&
            &  div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
            & -div_diff_flux_horz(jc,level,jb)                               &
            & ) ) / delta_z_new

          new_tracer_concentration(jc,level,jb) =         &
            & ( new_tracer_concentration(jc,level,jb) +   &
            & (delta_t  / delta_z_new) * top_bc(jc))
        END DO

        DO level = 2, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)

          new_tracer_concentration(jc,level,jb) =                          &
            &  old_tracer_concentration(jc,level,jb)*(stretch_c(jc, jb)/stretch_c_new(jc, jb)) -         &
            &  (delta_t /  ( stretch_c_new(jc, jb)*patch_3d%p_patch_1D(1)%prism_thick_c(jc,level,jb) ) ) &
            & * (div_adv_flux_horz(jc,level,jb)  + div_adv_flux_vert(jc,level,jb)                        &
            & - div_diff_flux_horz(jc,level,jb) )

        ENDDO

      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)
    END DO
    !ICON_OMP_END_PARALLEL_DO


    !Vertical mixing: implicit and with coefficient a_v
    !that is the sum of PP-coeff and implicit part of Redi-scheme


    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('BefImplDiff: div_adv_flux_vert',div_adv_flux_vert, str_module,idt_src, in_subset=cells_in_domain)
    CALL dbg_print('BefImplDiff: trac_inter', new_tracer_concentration,  str_module,idt_src, in_subset=cells_in_domain)
    !---------------------------------------------------------------------

    !calculate vert diffusion impicit: result is stored in trac_out
    ! no sync because of columnwise computation
    IF ( l_with_vert_tracer_diffusion ) THEN

      !Vertical mixing: implicit and with coefficient a_v
      CALL tracer_diffusion_vertical_implicit_zstar( &
           & patch_3d,                  &
           & new_tracer,                &
           & old_tracer%ver_diffusion_coeff, &
           & stretch_c, &
           & lacc=lzacc)

    ENDIF!IF ( l_with_vert_tracer_diffusion )

    CALL sync_patch_array(sync_c, patch_2D, new_tracer_concentration)
    !$ACC END DATA
  END SUBROUTINE advect_individual_tracers_zstar


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  SUBROUTINE advect_ocean_tracers_zstar(old_tracers, new_tracers, transport_state, operators_coeff, &
      & stretch_e, stretch_c, stretch_c_new, lacc)
    TYPE(t_tracer_collection), INTENT(inout)      :: old_tracers
    TYPE(t_tracer_collection), INTENT(inout)      :: new_tracers
    TYPE(t_ocean_transport_state), TARGET         :: transport_state
    TYPE(t_operator_coeff), INTENT(in) :: operators_coeff
    REAL(wp), INTENT(IN)               :: stretch_e(nproma, transport_state%patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(IN)               :: stretch_c(nproma, transport_state%patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(IN)               :: stretch_c_new(nproma, transport_state%patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(in), OPTIONAL :: lacc


    !Local variables
    TYPE(t_patch_3d ), POINTER     :: patch_3d
    INTEGER :: tracer_index
    LOGICAL :: lzacc
    !-------------------------------------------------------------------------------
    patch_3d => transport_state%patch_3d

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(patch_3d%p_patch_2d(1)%alloc_cell_blocks, patch_3d%p_patch_2d(1)%nblks_e) IF(lzacc)

    DO tracer_index = 1, old_tracers%no_of_tracers

      IF ( old_tracers%tracer(tracer_index)%is_advected) THEN

        call advect_individual_tracers_zstar( patch_3d, transport_state, &
          & operators_coeff, stretch_e, stretch_c, stretch_c_new,        &
          & old_tracers%tracer(tracer_index), new_tracers%tracer(tracer_index), lacc=lzacc)

      ENDIF

    END DO
    !$ACC END DATA

  END SUBROUTINE advect_ocean_tracers_zstar
  !-------------------------------------------------------------------------


END MODULE mo_ocean_tracer_zstar
