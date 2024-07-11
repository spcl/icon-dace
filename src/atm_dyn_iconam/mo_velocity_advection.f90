!
! This module contains the subroutine calculating the velocity advection tendencies
! for the nonhydrostatic dynamical core. Separated from mo_solve_nonhydro in order
! to speed up compile time
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_velocity_advection

  USE mo_kind,                 ONLY: wp, vp
  USE mo_nonhydrostatic_config,ONLY: lextra_diffu
  USE mo_parallel_config,   ONLY: nproma
  USE mo_run_config,        ONLY: lvert_nest, timers_level
  USE mo_model_domain,      ONLY: t_patch
  USE mo_intp_data_strc,    ONLY: t_int_state
  USE mo_icon_interpolation_scalar, ONLY: cells2verts_scalar_ri
  USE mo_nonhydro_types,    ONLY: t_nh_metrics, t_nh_diag, t_nh_prog
  USE mo_math_divrot,       ONLY: rot_vertex_ri
  USE mo_vertical_grid,     ONLY: nrdmax
  USE mo_init_vgrid,        ONLY: nflatlev
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,    ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_timer,             ONLY: timer_solve_nh_veltend, timer_start, timer_stop

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: velocity_tendencies

#ifdef _CRAYFTN
#define __CRAY_FTN_VERSION (_RELEASE_MAJOR * 100 + _RELEASE_MINOR)
#endif

  CONTAINS


  !----------------------------------------------------------------------------
  !>
  !! velocity_tendencies
  !!
  !! Discretization of nonhydrostatic momentum equation similar to hydrostatic core
  !! In particular, the Lamb transformation is applied only to the horizontal
  !! equation of motion, whereas the vertical wind equation is discretized
  !! in advective form
  !!
  SUBROUTINE velocity_tendencies (p_prog, p_patch, p_int, p_metrics, p_diag, z_w_concorr_me, z_kin_hor_e, &
                                  z_vt_ie, ntnd, istep, lvn_only, dtime, dt_linintp_ubc, ldeepatmo)

    ! Passed variables
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_int_state), TARGET, INTENT(IN):: p_int
    TYPE(t_nh_prog), INTENT(INOUT)       :: p_prog
    TYPE(t_nh_metrics), INTENT(INOUT)    :: p_metrics
    TYPE(t_nh_diag), INTENT(INOUT)       :: p_diag

    ! Local variables from solve_nh that are passed for efficiency optimization
    REAL(vp), DIMENSION(:,:,:), INTENT(INOUT) :: z_w_concorr_me, z_kin_hor_e, z_vt_ie

    INTEGER, INTENT(IN)  :: ntnd     ! time level of ddt fields used to store tendencies
    INTEGER, INTENT(IN)  :: istep    ! 1: predictor step, 2: corrector step
    LOGICAL, INTENT(IN)  :: lvn_only ! true: compute only vn tendency
    REAL(wp),INTENT(IN)  :: dtime    ! time step
    REAL(wp),INTENT(IN)  :: dt_linintp_ubc  ! time shift for upper boundary condition
    LOGICAL, INTENT(IN)  :: ldeepatmo! true: deep-atmosphere modification

    ! Local variables
    INTEGER :: jb, jk, jc, je
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_startblk_2, i_endblk_2, i_startidx_2, i_endidx_2
    INTEGER :: rl_start, rl_end, rl_start_2, rl_end_2
    ! The data type vp (variable precision) is by default the same as wp but reduces
    ! to single precision when the __MIXED_PRECISION cpp flag is set at compile time
    REAL(vp):: z_w_concorr_mc(nproma,p_patch%nlev)
    REAL(vp):: z_w_con_c(nproma,p_patch%nlevp1)
    REAL(vp):: z_w_con_c_full(nproma,p_patch%nlev,p_patch%nblks_c)
    ! These fields in addition have reversed index order (vertical first) for optimization
#ifdef __LOOP_EXCHANGE
    REAL(vp):: z_v_grad_w(p_patch%nlev,nproma,p_patch%nblks_e)
    REAL(vp):: z_w_v(p_patch%nlevp1,nproma,p_patch%nblks_v)
    REAL(vp):: zeta(p_patch%nlev,nproma,p_patch%nblks_v)
    REAL(vp):: z_ekinh(p_patch%nlev,nproma,p_patch%nblks_c)
#else
    REAL(vp):: z_v_grad_w(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(vp):: z_w_v(nproma,p_patch%nlevp1,p_patch%nblks_v)
    REAL(vp):: zeta(nproma,p_patch%nlev,p_patch%nblks_v)
    REAL(vp):: z_ekinh(nproma,p_patch%nlev,p_patch%nblks_c)
#endif

    ! Pointers
    INTEGER, DIMENSION(:,:,:), POINTER, CONTIGUOUS :: &
      icidx, icblk, ieidx, ieblk, iqidx, iqblk, ividx, ivblk, incidx, incblk

    INTEGER  :: nlev, nlevp1          !< number of full and half levels
    ! Local control variable for vertical nesting
    LOGICAL :: l_vert_nested

    INTEGER :: jg

    ! Variables for conditional additional diffusion for vertical advection
    REAL(vp) :: cfl_w_limit, vcfl, maxvcfl, vcflmax(p_patch%nblks_c)
    REAL(wp) :: w_con_e, scalfac_exdiff, difcoef, max_vcfl_dyn
                
    INTEGER  :: ie, nrdmax_jg, nflatlev_jg, clip_count
    LOGICAL  :: levmask(p_patch%nblks_c,p_patch%nlev),levelmask(p_patch%nlev)
    LOGICAL  :: cfl_clipping(nproma,p_patch%nlevp1)   ! CFL > 0.85

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: z_w_concorr_mc,z_w_con_c,z_w_con_c_full
!DIR$ ATTRIBUTES ALIGN :64 :: z_v_grad_w,z_w_v,zeta,z_ekinh
!DIR$ ATTRIBUTES ALIGN :64 :: levmask,cfl_clipping
#endif
    !--------------------------------------------------------------------------

    IF (timers_level > 5) CALL timer_start(timer_solve_nh_veltend)

    IF ((lvert_nest) .AND. (p_patch%nshift > 0)) THEN  
      l_vert_nested = .TRUE.
    ELSE
      l_vert_nested = .FALSE.
    ENDIF

    !Get patch id
    jg = p_patch%id
    nrdmax_jg     = nrdmax(jg)
    nflatlev_jg   = nflatlev(jg)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Set pointers to neighbor cells/edges/vertices
    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    incidx => p_patch%cells%neighbor_idx
    incblk => p_patch%cells%neighbor_blk

    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk

    !$ACC DATA COPYIN(z_w_concorr_me, z_kin_hor_e, z_vt_ie) &
    !$ACC   CREATE(z_w_concorr_mc, z_w_con_c, cfl_clipping, z_w_con_c_full, z_v_grad_w, z_w_v, zeta, z_ekinh, levmask, levelmask) &
    !$ACC   PRESENT(p_diag, p_prog, p_int, p_metrics, p_patch) &
    !$ACC   PRESENT(iqidx, iqblk, ividx, icblk, icidx, ieidx, ieblk, incblk, ivblk, incidx)

    ! Limit on vertical CFL number for applying extra diffusion
    IF (lextra_diffu) THEN
      cfl_w_limit = 0.65_wp/dtime   ! this means 65% of the nominal CFL stability limit

      ! Scaling factor for extra diffusion
      scalfac_exdiff = 0.05_wp / ( dtime*(0.85_wp - cfl_w_limit*dtime) )
    ELSE
      cfl_w_limit = 0.85_wp/dtime   ! this means 65% of the nominal CFL stability limit
      scalfac_exdiff = 0._wp
    ENDIF

    ! Compute w at vertices
    IF (.NOT. lvn_only) CALL cells2verts_scalar_ri(p_prog%w, p_patch, &
      p_int%cells_aw_verts, z_w_v, opt_rlend=min_rlvert_int-1, opt_acc_async=.TRUE.)

    ! Compute vertical vorticity component at vertices
    CALL rot_vertex_ri (p_prog%vn, p_patch, p_int, zeta, &
      opt_rlend=min_rlvert_int-1, opt_acc_async=.TRUE.)

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk, rl_start_2, rl_end_2, i_startblk_2, i_endblk_2)

    IF (istep == 1) THEN ! Computations of velocity-derived quantities that come from solve_nh in istep=2

      rl_start = 5
      rl_end = min_rledge_int - 2

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR TILE(128, 1)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
!$NEC outerloop_unroll(4)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! RBF reconstruction of tangential wind component
            p_diag%vt(je,jk,jb) = &
              p_int%rbf_vec_coeff_e(1,je,jb) * p_prog%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) + &
              p_int%rbf_vec_coeff_e(2,je,jb) * p_prog%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) + &
              p_int%rbf_vec_coeff_e(3,je,jb) * p_prog%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) + &
              p_int%rbf_vec_coeff_e(4,je,jb) * p_prog%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

! Fusing the following loop with this one is a lot more efficient on the GPU
#ifndef _OPENACC
          ENDDO
        ENDDO
        ! Interpolate vn to interface levels and compute horizontal part of kinetic energy on edges
        DO jk = 2, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
#else
            IF (jk >= 2) THEN
#endif
            p_diag%vn_ie(je,jk,jb) =                                    &
              p_metrics%wgtfac_e(je,jk,jb)*p_prog%vn(je,jk,jb) +        &
             (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_prog%vn(je,jk-1,jb)
            z_kin_hor_e(je,jk,jb) = 0.5_wp*(p_prog%vn(je,jk,jb)**2 + p_diag%vt(je,jk,jb)**2)
#ifdef _OPENACC
            ENDIF
#endif
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        IF (.NOT. lvn_only) THEN ! Interpolate also vt to interface levels

!WS: this gang loop is independent of the previous one and could execute concurrently
!    but overlapping is MUCH SLOWER with PGI

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 2, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              z_vt_ie(je,jk,jb) =                                         &
                p_metrics%wgtfac_e(je,jk,jb)*p_diag%vt(je,jk,jb) +        &
               (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_diag%vt(je,jk-1,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        ! Compute contravariant correction for vertical velocity at interface levels
        ! (will be interpolated to cell centers below)

!WS: this gang loop is independent of the previous one and could execute concurrently
!    but overlapping is MUCH SLOWER with PGI
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = nflatlev_jg, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            z_w_concorr_me(je,jk,jb) =                              &
              p_prog%vn(je,jk,jb)*p_metrics%ddxn_z_full(je,jk,jb) + &
              p_diag%vt(je,jk,jb)*p_metrics%ddxt_z_full(je,jk,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        IF (.NOT. l_vert_nested) THEN

          ! Top and bottom levels
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            ! Quadratic extrapolation at the top turned out to cause numerical instability in pathological cases,
            ! thus we use a no-gradient condition in the upper half layer
            p_diag%vn_ie(je,1,jb) = p_prog%vn(je,1,jb)
            ! vt_ie(jk=1) is actually unused, but we need it for convenience of implementation
            z_vt_ie(je,1,jb) = p_diag%vt(je,1,jb)
            !
            z_kin_hor_e(je,1,jb) = 0.5_wp*(p_prog%vn(je,1,jb)**2 + p_diag%vt(je,1,jb)**2)
            p_diag%vn_ie(je,nlevp1,jb) =                           &
              p_metrics%wgtfacq_e(je,1,jb)*p_prog%vn(je,nlev,jb) +   &
              p_metrics%wgtfacq_e(je,2,jb)*p_prog%vn(je,nlev-1,jb) + &
              p_metrics%wgtfacq_e(je,3,jb)*p_prog%vn(je,nlev-2,jb)
          ENDDO
          !$ACC END PARALLEL

        ELSE

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_diag%vn_ie(je,1,jb) = p_diag%vn_ie_ubc(je,1,jb)+dt_linintp_ubc*p_diag%vn_ie_ubc(je,2,jb)
            ! vt_ie(jk=1) is actually unused, but we need it for convenience of implementation
            z_vt_ie(je,1,jb) = p_diag%vt(je,1,jb)
            !
            z_kin_hor_e(je,1,jb) = 0.5_wp*(p_prog%vn(je,1,jb)**2 + p_diag%vt(je,1,jb)**2)
            p_diag%vn_ie(je,nlevp1,jb) =                           &
              p_metrics%wgtfacq_e(je,1,jb)*p_prog%vn(je,nlev,jb) +   &
              p_metrics%wgtfacq_e(je,2,jb)*p_prog%vn(je,nlev-1,jb) + &
              p_metrics%wgtfacq_e(je,3,jb)*p_prog%vn(je,nlev-2,jb)
          ENDDO
          !$ACC END PARALLEL
        ENDIF

      ENDDO
!$OMP END DO

    ENDIF ! istep = 1

    rl_start = 7
    rl_end = min_rledge_int - 1

    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

    IF (.NOT. lvn_only) THEN
!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Compute v*grad w on edges (level nlevp1 is not needed because w(nlevp1) is diagnostic)
        ! Note: this implicitly includes a minus sign for the gradients, which is needed later on

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
            z_v_grad_w(jk,je,jb) = p_diag%vn_ie(je,jk,jb) * p_patch%edges%inv_dual_edge_length(je,jb)* &
             (p_prog%w(icidx(je,jb,1),jk,icblk(je,jb,1)) - p_prog%w(icidx(je,jb,2),jk,icblk(je,jb,2))) &
             + z_vt_ie(je,jk,jb) * p_patch%edges%inv_primal_edge_length(je,jb) *                       &
             p_patch%edges%tangent_orientation(je,jb) *                                                 &
             (z_w_v(jk,ividx(je,jb,1),ivblk(je,jb,1)) - z_w_v(jk,ividx(je,jb,2),ivblk(je,jb,2))) 
#else
!$NEC outerloop_unroll(2)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            z_v_grad_w(je,jk,jb) = p_diag%vn_ie(je,jk,jb) * p_patch%edges%inv_dual_edge_length(je,jb)* &
             (p_prog%w(icidx(je,jb,1),jk,icblk(je,jb,1)) - p_prog%w(icidx(je,jb,2),jk,icblk(je,jb,2))) &
             + z_vt_ie(je,jk,jb) * p_patch%edges%inv_primal_edge_length(je,jb) *                       &
             p_patch%edges%tangent_orientation(je,jb) *                                                 &
             (z_w_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) - z_w_v(ividx(je,jb,2),jk,ivblk(je,jb,2))) 
#endif
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!$OMP END DO

    ENDIF

    IF (.NOT. lvn_only .AND. ldeepatmo) THEN
!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
        ! -[ vn * ( dw/dn - vn / r + ft ) + vt * ( dw/dt - vt / r - fn ) ]
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
            z_v_grad_w(jk,je,jb) = z_v_grad_w(jk,je,jb) * p_metrics%deepatmo_gradh_ifc(jk)              &
             + p_diag%vn_ie(je,jk,jb) *                                                                 &
               ( p_diag%vn_ie(je,jk,jb) * p_metrics%deepatmo_invr_ifc(jk) - p_patch%edges%ft_e(je,jb) ) &
             + z_vt_ie(je,jk,jb) *                                                                      &
               ( z_vt_ie(je,jk,jb) * p_metrics%deepatmo_invr_ifc(jk) + p_patch%edges%fn_e(je,jb) )
#else
!$NEC outerloop_unroll(2)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            z_v_grad_w(je,jk,jb) = z_v_grad_w(je,jk,jb) * p_metrics%deepatmo_gradh_ifc(jk)              &
             + p_diag%vn_ie(je,jk,jb) *                                                                 &
               ( p_diag%vn_ie(je,jk,jb) * p_metrics%deepatmo_invr_ifc(jk) - p_patch%edges%ft_e(je,jb) ) &
             + z_vt_ie(je,jk,jb) *                                                                      &
               ( z_vt_ie(je,jk,jb) * p_metrics%deepatmo_invr_ifc(jk) + p_patch%edges%fn_e(je,jb) )
#endif
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!$OMP END DO
    ENDIF

    rl_start = 4
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    rl_start_2 = grf_bdywidth_c+1
    rl_end_2   = min_rlcell_int

    i_startblk_2 = p_patch%cells%start_block(rl_start_2)
    i_endblk_2   = p_patch%cells%end_block(rl_end_2)

!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, i_startidx_2, i_endidx_2, z_w_con_c, &
!$OMP            z_w_concorr_mc, difcoef, vcfl, maxvcfl, cfl_clipping, clip_count) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Interpolate horizontal kinetic energy to cell centers

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
!DIR$ IVDEP
        DO jk = 1, nlev
        z_ekinh(jk,jc,jb) =  &
#else
!$NEC outerloop_unroll(4)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
        z_ekinh(jc,jk,jb) =  &
#endif
          p_int%e_bln_c_s(jc,1,jb)*z_kin_hor_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
          p_int%e_bln_c_s(jc,2,jb)*z_kin_hor_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
          p_int%e_bln_c_s(jc,3,jb)*z_kin_hor_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      IF (istep == 1) THEN

        ! Interpolate contravariant correction to cell centers ...
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = nflatlev_jg, nlev
#else
!$NEC outerloop_unroll(4)
        DO jk = nflatlev_jg, nlev
          DO jc = i_startidx, i_endidx
#endif

            z_w_concorr_mc(jc,jk) =  &
              p_int%e_bln_c_s(jc,1,jb)*z_w_concorr_me(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
              p_int%e_bln_c_s(jc,2,jb)*z_w_concorr_me(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
              p_int%e_bln_c_s(jc,3,jb)*z_w_concorr_me(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))

          ENDDO
        ENDDO
        !$ACC END PARALLEL

        ! ... and to interface levels
        ! Remark: computation of w_concorr_c at nlevp1 is needed in solve_nh only
        ! because this serves solely for setting the lower boundary condition for w

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = nflatlev_jg+1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_diag%w_concorr_c(jc,jk,jb) =                                &
              p_metrics%wgtfac_c(jc,jk,jb)*z_w_concorr_mc(jc,jk) +        &
             (1._vp - p_metrics%wgtfac_c(jc,jk,jb))*z_w_concorr_mc(jc,jk-1)
          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ENDIF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_w_con_c(jc,jk) =  p_prog%w(jc,jk,jb)
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        z_w_con_c(jc,nlevp1) = 0.0_wp
      ENDDO
      !$ACC END PARALLEL

      ! Contravariant vertical velocity on w points and interpolation to full levels
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = nlev, nflatlev_jg+1, -1
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_w_con_c(jc,jk) = z_w_con_c(jc,jk) - p_diag%w_concorr_c(jc,jk,jb)
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      ! Search for grid points for which w_con is close to or above the CFL stability limit
      ! At these points, additional diffusion is applied in order to prevent numerical 
      ! instability if lextra_diffu = .TRUE.
      ! WS:  We split out levmask in order to collapse the subsequent loop, and avoid problems with two levels of REDUCTION
      !$ACC LOOP
      DO jk = MAX(3,nrdmax_jg-2), nlev-3
        levmask(jb,jk) = .FALSE.
      ENDDO
      !$ACC END PARALLEL

        ! Initialize the max CFL variable to zero for the following computation
        ! TODO: in OpenACC this logic can be moved to GPU completely
        maxvcfl = 0

! DA this kernel is ASYNC(1), so need to wait to retrieve the value
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) REDUCTION(MAX: maxvcfl)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(vcfl) REDUCTION(MAX: maxvcfl)
      DO jk = MAX(3,nrdmax_jg-2), nlev-3
#ifndef _OPENACC
        clip_count = 0
        DO jc = i_startidx, i_endidx
          cfl_clipping(jc,jk) = (ABS(z_w_con_c(jc,jk)) > cfl_w_limit*p_metrics%ddqz_z_half(jc,jk,jb))
          IF (cfl_clipping(jc,jk)) clip_count = clip_count+1
        ENDDO
        IF (clip_count == 0) CYCLE
        DO jc = i_startidx, i_endidx
#else
        DO jc = i_startidx, i_endidx
          cfl_clipping(jc,jk) = (ABS(z_w_con_c(jc,jk)) > cfl_w_limit*p_metrics%ddqz_z_half(jc,jk,jb))
#endif
          IF ( cfl_clipping(jc,jk) ) THEN
      ! WS:  setting levmask cannot create race conditions; following should be fine
            levmask(jb,jk) = .TRUE.
            vcfl = z_w_con_c(jc,jk)*dtime/p_metrics%ddqz_z_half(jc,jk,jb)
            maxvcfl = MAX( maxvcfl, ABS( vcfl ) )
            !
            ! limit w_con to 85% of the nominal CFL stability threshold
            IF (vcfl < -0.85_vp) THEN
              z_w_con_c(jc,jk)           = -0.85_vp*p_metrics%ddqz_z_half(jc,jk,jb)/dtime
            ELSE IF (vcfl > 0.85_vp) THEN
              z_w_con_c(jc,jk)           = 0.85_vp*p_metrics%ddqz_z_half(jc,jk,jb)/dtime
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(8)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          z_w_con_c_full(jc,jk,jb) = 0.5_vp*(z_w_con_c(jc,jk)+z_w_con_c(jc,jk+1))
        ENDDO
      ENDDO
      !$ACC END PARALLEL

! Get back the reduced max CFL
      !$ACC WAIT
      vcflmax(jb) = maxvcfl


      ! The remaining computations are not needed in vn_only mode and only on prognostic grid points
      IF (lvn_only) CYCLE
      IF (jb < i_startblk_2 .OR. jb > i_endblk_2) CYCLE

      CALL get_indices_c(p_patch, jb, i_startblk_2, i_endblk_2, &
                         i_startidx_2, i_endidx_2, rl_start_2, rl_end_2)

      ! Compute vertical derivative terms of vertical wind advection
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
!$NEC outerloop_unroll(8)
      DO jk = 2, nlev
!DIR$ IVDEP
        DO jc = i_startidx_2, i_endidx_2
          p_diag%ddt_w_adv_pc(jc,jk,jb,ntnd) =  - z_w_con_c(jc,jk)*                                 &
            (p_prog%w(jc,jk-1,jb)*p_metrics%coeff1_dwdz(jc,jk,jb) -                                 &
             p_prog%w(jc,jk+1,jb)*p_metrics%coeff2_dwdz(jc,jk,jb) +                                 &
             p_prog%w(jc,jk,jb)*(p_metrics%coeff2_dwdz(jc,jk,jb) - p_metrics%coeff1_dwdz(jc,jk,jb)) )
        ENDDO
      ENDDO

      ! Interpolate horizontal advection of w from edges to cells and add to advective tendency
      !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx_2, i_endidx_2
!DIR$ IVDEP
        DO jk = 2, nlev
          p_diag%ddt_w_adv_pc(jc,jk,jb,ntnd) = p_diag%ddt_w_adv_pc(jc,jk,jb,ntnd) + &
            p_int%e_bln_c_s(jc,1,jb)*z_v_grad_w(jk,ieidx(jc,jb,1),ieblk(jc,jb,1)) + &
            p_int%e_bln_c_s(jc,2,jb)*z_v_grad_w(jk,ieidx(jc,jb,2),ieblk(jc,jb,2)) + &
            p_int%e_bln_c_s(jc,3,jb)*z_v_grad_w(jk,ieidx(jc,jb,3),ieblk(jc,jb,3))
#else
!$NEC outerloop_unroll(4)
      DO jk = 2, nlev
        DO jc = i_startidx_2, i_endidx_2
          p_diag%ddt_w_adv_pc(jc,jk,jb,ntnd) = p_diag%ddt_w_adv_pc(jc,jk,jb,ntnd) + &
            p_int%e_bln_c_s(jc,1,jb)*z_v_grad_w(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
            p_int%e_bln_c_s(jc,2,jb)*z_v_grad_w(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
            p_int%e_bln_c_s(jc,3,jb)*z_v_grad_w(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))
#endif
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      IF (lextra_diffu) THEN

        ! Apply extra diffusion at grid points where w_con is close to or above the CFL stability limit
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG
        DO jk = MAX(3,nrdmax_jg-2), nlev-3
          IF (levmask(jb,jk)) THEN
            !$ACC LOOP VECTOR PRIVATE(difcoef)
            DO jc = i_startidx_2, i_endidx_2
              IF (cfl_clipping(jc,jk) .AND. p_patch%cells%decomp_info%owner_mask(jc,jb)) THEN
                difcoef = scalfac_exdiff * MIN(0.85_wp - cfl_w_limit*dtime,                       &
                  ABS(z_w_con_c(jc,jk))*dtime/p_metrics%ddqz_z_half(jc,jk,jb) - cfl_w_limit*dtime )

                ! nabla2 diffusion on w
                p_diag%ddt_w_adv_pc(jc,jk,jb,ntnd) = p_diag%ddt_w_adv_pc(jc,jk,jb,ntnd)  + &
                  difcoef * p_patch%cells%area(jc,jb) * (                                  &
                  p_prog%w(jc,jk,jb)                          *p_int%geofac_n2s(jc,1,jb) + &
                  p_prog%w(incidx(jc,jb,1),jk,incblk(jc,jb,1))*p_int%geofac_n2s(jc,2,jb) + &
                  p_prog%w(incidx(jc,jb,2),jk,incblk(jc,jb,2))*p_int%geofac_n2s(jc,3,jb) + &
                  p_prog%w(incidx(jc,jb,3),jk,incblk(jc,jb,3))*p_int%geofac_n2s(jc,4,jb)   )

              ENDIF
            ENDDO
          ENDIF
        ENDDO
        !$ACC END PARALLEL

      ENDIF

    ENDDO
!$OMP END DO


    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
!$OMP DO PRIVATE(jk)
    DO jk = MAX(3,nrdmax_jg-2), nlev-3
      levelmask(jk) = ANY(levmask(i_startblk:i_endblk,jk))
    ENDDO
!$OMP END DO
    !$ACC END PARALLEL

    rl_start = grf_bdywidth_e+1
    rl_end = min_rledge_int

    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx, ie, w_con_e, difcoef) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Sum up terms of horizontal wind advection: grad(Ekin_h) + vt*(f+relvort_e) + wcon_e*dv/dz
      IF (.NOT. ldeepatmo) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
#ifdef _CRAYFTN
!DIR$ PREFERVECTOR
#endif
          DO jk = 1, nlev
            p_diag%ddt_vn_apc_pc(je,jk,jb,ntnd) = - ( z_kin_hor_e(je,jk,jb) *                     &
             (p_metrics%coeff_gradekin(je,1,jb) - p_metrics%coeff_gradekin(je,2,jb)) +            &
              p_metrics%coeff_gradekin(je,2,jb)*z_ekinh(jk,icidx(je,jb,2),icblk(je,jb,2)) -       &
              p_metrics%coeff_gradekin(je,1,jb)*z_ekinh(jk,icidx(je,jb,1),icblk(je,jb,1)) +       &
              p_diag%vt(je,jk,jb) * ( p_patch%edges%f_e(je,jb) + 0.5_vp*                          &
             (zeta(jk,ividx(je,jb,1),ivblk(je,jb,1)) + zeta(jk,ividx(je,jb,2),ivblk(je,jb,2)))) + &
             (p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) +           &
              p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2)))*           &
             (p_diag%vn_ie(je,jk,jb) - p_diag%vn_ie(je,jk+1,jb))/p_metrics%ddqz_z_full_e(je,jk,jb))
          ENDDO
        ENDDO
#else
!$NEC outerloop_unroll(3)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            p_diag%ddt_vn_apc_pc(je,jk,jb,ntnd) = - ( z_kin_hor_e(je,jk,jb) *                     &
             (p_metrics%coeff_gradekin(je,1,jb) - p_metrics%coeff_gradekin(je,2,jb)) +            &
              p_metrics%coeff_gradekin(je,2,jb)*z_ekinh(icidx(je,jb,2),jk,icblk(je,jb,2)) -       &
              p_metrics%coeff_gradekin(je,1,jb)*z_ekinh(icidx(je,jb,1),jk,icblk(je,jb,1)) +       &
              p_diag%vt(je,jk,jb) * ( p_patch%edges%f_e(je,jb) + 0.5_vp*                          &
             (zeta(ividx(je,jb,1),jk,ivblk(je,jb,1)) + zeta(ividx(je,jb,2),jk,ivblk(je,jb,2)))) + &
             (p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) +           &
              p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2)))*           &
             (p_diag%vn_ie(je,jk,jb) - p_diag%vn_ie(je,jk+1,jb))/p_metrics%ddqz_z_full_e(je,jk,jb))
          ENDDO
        ENDDO
#endif
        !$ACC END PARALLEL

        ! If needed, compute separately the Coriolis effect: vt*f
        IF (p_diag%ddt_vn_adv_is_associated .OR. p_diag%ddt_vn_cor_is_associated) THEN
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR TILE(32, 4)
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              p_diag%ddt_vn_cor_pc(je,jk,jb,ntnd) = - p_diag%vt(je,jk,jb) * p_patch%edges%f_e(je,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        END IF

      ELSE
        ! Deep atmosphere: grad(Ekin_h) is multiplied by metrical modification factor to account 
        ! for spherical geometry, in addition metrical terms and contribution of vertical wind to 
        ! Coriolis acceleration have been added: wcon_e * dvn/dz -> wcon_e * ( dvn/dz + vn / r - ft ), 
        ! where r is radius and ft is tangential component of horizontal Coriolis parameter. 
        ! The vorticity 'zeta' has to be multiplied by 'deepatmo_gradh_mc(jk)', 
        ! because the subroutine 'rot_vertex_ri', which computes 'zeta', is itself not 
        ! modified for spherical geometry.

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
#ifdef _CRAYFTN
!DIR$ PREFERVECTOR
#endif
          DO jk = 1, nlev
            p_diag%ddt_vn_apc_pc(je,jk,jb,ntnd) = - (                                             & 
               (                                                                                  & 
                z_kin_hor_e(je,jk,jb) *                                                           &
                (p_metrics%coeff_gradekin(je,1,jb) - p_metrics%coeff_gradekin(je,2,jb)) +         &
                p_metrics%coeff_gradekin(je,2,jb)*z_ekinh(jk,icidx(je,jb,2),icblk(je,jb,2)) -     &
                p_metrics%coeff_gradekin(je,1,jb)*z_ekinh(jk,icidx(je,jb,1),icblk(je,jb,1))       &
               ) * p_metrics%deepatmo_gradh_mc(jk)                                                &
               + p_diag%vt(je,jk,jb) * ( p_patch%edges%f_e(je,jb) + 0.5_vp*                       &
               (zeta(jk,ividx(je,jb,1),ivblk(je,jb,1)) + zeta(jk,ividx(je,jb,2),ivblk(je,jb,2)))  &
               * p_metrics%deepatmo_gradh_mc(jk) )                                                &
               + (p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) +       &
               p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2))) *         & 
               (                                                                                  &
                (p_diag%vn_ie(je,jk,jb) - p_diag%vn_ie(je,jk+1,jb))/p_metrics%ddqz_z_full_e(je,jk,jb) &
                + p_prog%vn(je,jk,jb) * p_metrics%deepatmo_invr_mc(jk)                            & 
                - p_patch%edges%ft_e(je,jb)                                                       &
               )                                                                                  &
              )
          ENDDO
        ENDDO
#else
!$NEC outerloop_unroll(3)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            p_diag%ddt_vn_apc_pc(je,jk,jb,ntnd) = - (                                             &
               (                                                                                  &
                z_kin_hor_e(je,jk,jb) *                                                           &
                (p_metrics%coeff_gradekin(je,1,jb) - p_metrics%coeff_gradekin(je,2,jb)) +         &
                p_metrics%coeff_gradekin(je,2,jb)*z_ekinh(icidx(je,jb,2),jk,icblk(je,jb,2)) -     &
                p_metrics%coeff_gradekin(je,1,jb)*z_ekinh(icidx(je,jb,1),jk,icblk(je,jb,1))       &         
               ) * p_metrics%deepatmo_gradh_mc(jk)                                                &
               + p_diag%vt(je,jk,jb) * ( p_patch%edges%f_e(je,jb) + 0.5_vp*                       &
               (zeta(ividx(je,jb,1),jk,ivblk(je,jb,1)) + zeta(ividx(je,jb,2),jk,ivblk(je,jb,2)))  &
               * p_metrics%deepatmo_gradh_mc(jk) )                                                &
               + (p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) +       &
               p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2))) *         &
               (                                                                                  &
               (p_diag%vn_ie(je,jk,jb) - p_diag%vn_ie(je,jk+1,jb))/p_metrics%ddqz_z_full_e(je,jk,jb) & 
               + p_prog%vn(je,jk,jb) * p_metrics%deepatmo_invr_mc(jk)                             &
               - p_patch%edges%ft_e(je,jb)                                                        &
               )                                                                                  &
              )
          ENDDO
        ENDDO
#endif
        !$ACC END PARALLEL

        ! If needed, compute separately the Coriolis effect: vt*f
        IF (p_diag%ddt_vn_adv_is_associated .OR. p_diag%ddt_vn_cor_is_associated) THEN
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR TILE(32, 4)
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              p_diag%ddt_vn_cor_pc(je,jk,jb,ntnd) = - (                                        &
                 + p_diag%vt(je,jk,jb) * ( p_patch%edges%f_e(je,jb) )                          &
                 + (p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) +  &
                 p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2))) *    &
                 (                                                                             &
                 - p_patch%edges%ft_e(je,jb)                                                   &
                 )                                                                             &
                )
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        END IF

      ENDIF

      IF (lextra_diffu) THEN
        ! Search for grid points for which w_con is close to or above the CFL stability limit
        ! At these points, additional diffusion is applied in order to prevent numerical instability

        ie = 0
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG
        DO jk = MAX(3,nrdmax_jg-2), nlev-4
          IF (levelmask(jk) .OR. levelmask(jk+1)) THEN
            !$ACC LOOP VECTOR PRIVATE(difcoef, w_con_e)
            DO je = i_startidx, i_endidx
              w_con_e = p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) + &
                        p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2))
              IF (ABS(w_con_e) > cfl_w_limit*p_metrics%ddqz_z_full_e(je,jk,jb)) THEN
                difcoef = scalfac_exdiff * MIN(0.85_wp - cfl_w_limit*dtime,                &
                  ABS(w_con_e)*dtime/p_metrics%ddqz_z_full_e(je,jk,jb) - cfl_w_limit*dtime )

                p_diag%ddt_vn_apc_pc(je,jk,jb,ntnd) = p_diag%ddt_vn_apc_pc(je,jk,jb,ntnd) +             &
                  difcoef * p_patch%edges%area_edge(je,jb) * (                                          &
                  p_int%geofac_grdiv(je,1,jb)*p_prog%vn(je,jk,jb)                         +             &
                  p_int%geofac_grdiv(je,2,jb)*p_prog%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) +             &
                  p_int%geofac_grdiv(je,3,jb)*p_prog%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) +             &
                  p_int%geofac_grdiv(je,4,jb)*p_prog%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) +             &
                  p_int%geofac_grdiv(je,5,jb)*p_prog%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4)) +             &
                  p_patch%edges%tangent_orientation(je,jb)*p_patch%edges%inv_primal_edge_length(je,jb) * &
#ifdef __LOOP_EXCHANGE
                  (zeta(jk,ividx(je,jb,2),ivblk(je,jb,2)) - zeta(jk,ividx(je,jb,1),ivblk(je,jb,1))) )
#else
                  (zeta(ividx(je,jb,2),jk,ivblk(je,jb,2)) - zeta(ividx(je,jb,1),jk,ivblk(je,jb,1))) )
#endif
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        !$ACC END PARALLEL

      ENDIF

    ENDDO
!$OMP END DO 
!$OMP END PARALLEL

    !$ACC WAIT

    ! Save maximum vertical CFL number for substep number adaptation
    i_startblk = p_patch%cells%start_block(grf_bdywidth_c)
    i_endblk   = p_patch%cells%end_block(min_rlcell_int)

    ! DA: vcflmax already on cpu
    max_vcfl_dyn = MAX(p_diag%max_vcfl_dyn,MAXVAL(vcflmax(i_startblk:i_endblk)))

    ! DA: This is also on the CPU, should just stay there
    p_diag%max_vcfl_dyn = max_vcfl_dyn

    !$ACC END DATA

    IF (timers_level > 5) CALL timer_stop(timer_solve_nh_veltend)

  END SUBROUTINE velocity_tendencies

#ifdef _OPENACC
     SUBROUTINE h2d_velocity_tendencies( ntnd, p_prog, p_diag, z_w_concorr_me, z_kin_hor_e, z_vt_ie )
       INTEGER, INTENT(IN)                       :: ntnd
       TYPE(t_nh_prog), INTENT(INOUT)            :: p_prog
       TYPE(t_nh_diag), INTENT(INOUT)            :: p_diag
       REAL(vp), DIMENSION(:,:,:), INTENT(INOUT) :: z_w_concorr_me, z_kin_hor_e, z_vt_ie


       REAL(wp), DIMENSION(:,:,:),   POINTER  :: vn_tmp, w_tmp
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: vt_tmp, vn_ie_tmp
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: w_concorr_c_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_vn_apc_pc_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_vn_cor_pc_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_w_adv_pc_tmp

       vn_tmp              => p_prog%vn
       w_tmp               => p_prog%w
       vt_tmp              => p_diag%vt
       vn_ie_tmp           => p_diag%vn_ie
       w_concorr_c_tmp     => p_diag%w_concorr_c
       !$ACC UPDATE &
       !$ACC   DEVICE(vn_tmp, w_tmp, vt_tmp, vn_ie_tmp, w_concorr_c_tmp) &
       !$ACC   DEVICE(z_w_concorr_me, z_kin_hor_e, z_vt_ie) &
       !$ACC   ASYNC(1)

       ddt_vn_apc_pc_tmp   => p_diag%ddt_vn_apc_pc
       ddt_w_adv_pc_tmp    => p_diag%ddt_w_adv_pc
       !$ACC UPDATE DEVICE(ddt_vn_apc_pc_tmp(:,:,:,ntnd), ddt_w_adv_pc_tmp(:,:,:,ntnd)) ASYNC(1)
       IF (p_diag%ddt_vn_adv_is_associated .OR. p_diag%ddt_vn_cor_is_associated) THEN
          ddt_vn_cor_pc_tmp   => p_diag%ddt_vn_cor_pc
       !$ACC UPDATE DEVICE(ddt_vn_cor_pc_tmp(:,:,:,ntnd)) ASYNC(1)
       END IF

     END SUBROUTINE h2d_velocity_tendencies

     SUBROUTINE d2h_velocity_tendencies( istep, ntnd, p_diag, z_w_concorr_me, z_kin_hor_e, z_vt_ie )

       INTEGER, INTENT(IN)                       :: istep, ntnd
       TYPE(t_nh_diag), INTENT(INOUT)            :: p_diag
       REAL(vp), DIMENSION(:,:,:), INTENT(INOUT) :: z_w_concorr_me, z_kin_hor_e, z_vt_ie

       REAL(vp), DIMENSION(:,:,:),   POINTER  :: vt_tmp, vn_ie_tmp
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: w_concorr_c_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_vn_apc_pc_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_vn_cor_pc_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_w_adv_pc_tmp

       vt_tmp              => p_diag%vt
       vn_ie_tmp           => p_diag%vn_ie
       w_concorr_c_tmp     => p_diag%w_concorr_c
       !$ACC UPDATE &
       !$ACC   HOST(z_kin_hor_e, z_vt_ie, z_w_concorr_me, vt_tmp, vn_ie_tmp, w_concorr_c_tmp) &
       !$ACC   ASYNC(1) IF(istep==1)

       ddt_vn_apc_pc_tmp   => p_diag%ddt_vn_apc_pc
       ddt_w_adv_pc_tmp    => p_diag%ddt_w_adv_pc
       !$ACC UPDATE HOST(ddt_vn_apc_pc_tmp(:,:,:,ntnd), ddt_w_adv_pc_tmp(:,:,:,ntnd)) ASYNC(1)
       IF (p_diag%ddt_vn_adv_is_associated .OR. p_diag%ddt_vn_cor_is_associated) THEN
          ddt_vn_cor_pc_tmp   => p_diag%ddt_vn_cor_pc
       !$ACC UPDATE HOST(ddt_vn_cor_pc_tmp(:,:,:,ntnd)) ASYNC(1)
       END IF

       !$ACC WAIT(1)

     END SUBROUTINE d2h_velocity_tendencies
#endif

END MODULE mo_velocity_advection
