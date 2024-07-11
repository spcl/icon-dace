! Flux limiter for horizontal tracer transport
!
! This module contains flux limiters for horizontal
! tracer transport.
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
#define LAXFR_UPFLUX_MACRO(PPp_vn,PPp_psi_a,PPp_psi_b) (0.5_wp*((PPp_vn)*((PPp_psi_a)+(PPp_psi_b))-ABS(PPp_vn)*((PPp_psi_b)-(PPp_psi_a))))
#define LAXFR_UPFLUX_V_MACRO(PPp_w,PPp_psi_a,PPp_psi_b) (0.5_wp*((PPp_w)*((PPp_psi_a)+(PPp_psi_b))+ABS(PPp_w)*((PPp_psi_b)-(PPp_psi_a))))

#ifdef __INTEL_COMPILER
#define USE_LAXFR_MACROS
#define laxfr_upflux LAXFR_UPFLUX_MACRO
#endif
!----------------------------
MODULE mo_advection_hlimit

  USE mo_kind,                ONLY: wp, vp
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_fortran_tools,       ONLY: init
  USE mo_model_domain,        ONLY: t_patch, get_startrow_c
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_sync,                ONLY: SYNC_C1, sync_patch_array, &
    &                               sync_patch_array_mult
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_impl_constants,      ONLY: min_rledge_int, min_rlcell_int
#ifndef USE_LAXFR_MACROS
  USE mo_advection_utils,     ONLY: laxfr_upflux
#endif

  IMPLICIT NONE

  PRIVATE


  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: hflx_limiter_mo
  PUBLIC :: hflx_limiter_pd


CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Flux limiter for horizontal advection
  !!
  !! Zalesak Flux-Limiter (Flux corrected transport)
  !! The corrected flux is a weighted average of the low order flux and the
  !! given high order flux. The high order flux is used to the greatest extent
  !! possible without introducing overshoots and undershoots.
  !! Note: This limiter is positive definite and almost monotone (but not strictly).
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Schaer, C. and P.K. Smolarkiewicz (1996): A synchronous and iterative 
  !!   flux-correction formalism for coupled transport equations. J. comput. Phys., 
  !!   128, 101-120
  !!
  SUBROUTINE hflx_limiter_mo( ptr_patch, ptr_int, p_dtime, p_cc,            &
    &                         p_rhodz_now, p_rhodz_new, p_mass_flx_e,       &
    &                         p_mflx_tracer_h, slev, elev, opt_beta_fct,    &
    &                         opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_int_state), TARGET, INTENT(IN) :: & !< pointer to data structure for
      &  ptr_int                               !< interpolation

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable
      &  p_cc(:,:,:)                 !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< density times cell thickness at timestep n
      &  p_rhodz_now(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< density times cell thickness at timestep n+1
      &  p_rhodz_new(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(in) ::     &    !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)         !< (provided by dynamical core)
                                     !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::     &    !< time step
      &  p_dtime

    REAL(wp), INTENT(INOUT) ::  &    !< calculated horizontal tracer mass flux
      &  p_mflx_tracer_h(:,:,:)      !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN) ::      &    !< vertical start level
      &  slev

    INTEGER, INTENT(IN) ::      &    !< vertical end level
      &  elev

    REAL(wp), INTENT(IN), OPTIONAL ::  & !< factor for multiplicative spreading of range 
      &  opt_beta_fct                    !< of permissible values

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)


    REAL(wp) ::                 &    !< first order tracer mass flux
      &  z_mflx_low(nproma,slev:elev,ptr_patch%nblks_e)

    REAL(wp) ::                 &    !< antidiffusive tracer mass flux (F_H - F_L)
      &  z_anti(nproma,slev:elev,ptr_patch%nblks_e)

    REAL(vp) ::                 &    !< antidiffusive tracer mass flux (F_H - F_L)
      &  z_mflx_anti_1,         &    !< (units kg/kg)
      &  z_mflx_anti_2,         &
      &  z_mflx_anti_3

    REAL(vp) ::                 &    !< sum of incoming antidiffusive tracer mass fluxes
      &  z_mflx_anti_in (nproma,slev:elev,ptr_patch%nblks_c) !< (units kg/kg)

    REAL(vp) ::                 &    !< sum of outgoing antidiffusive tracer mass fluxes
      &  z_mflx_anti_out(nproma,slev:elev,ptr_patch%nblks_c) !< (units kg/kg)

    REAL(vp) ::                 &    !< flux divergence at cell center
      &  z_fluxdiv_c

    REAL(wp) ::                 &    !< new tracer field after hor. transport,
      &  z_tracer_new_low(nproma,slev:elev,ptr_patch%nblks_c) 
                                     !< if the low order fluxes are used

    REAL(vp) ::                 &    !< local maximum of current tracer value and low
      &  z_tracer_max(nproma,slev:elev,ptr_patch%nblks_c) !< order update

    REAL(vp) ::                 &    !< local minimum of current tracer value and low
      &  z_tracer_min(nproma,slev:elev,ptr_patch%nblks_c) !< order update

    ! remark: single precision would be sufficient for r_m and r_p, but SP-sync is not yet available
    REAL(wp) ::                 &    !< fraction which must multiply all in/out fluxes 
      &  r_p(nproma,slev:elev,ptr_patch%nblks_c),&   !< of cell jc to guarantee
      &  r_m(nproma,slev:elev,ptr_patch%nblks_c)     !< no overshoot/undershoot

    REAL(wp) :: r_frac !< computed minimum fraction which must multiply
                       !< the flux at the edge

    REAL(vp) :: z_min(nproma,slev:elev), & !< minimum/maximum value in cell and neighboring cells
      &         z_max(nproma,slev:elev) 
    REAL(wp) :: z_signum             !< sign of antidiffusive velocity
    REAL(wp) :: beta_fct             !< factor of allowed over-/undershooting in monotonous limiter
    REAL(wp) :: r_beta_fct           !< ... and its reverse value   

    INTEGER, CONTIGUOUS, POINTER :: &  !< Pointer to line and block indices of two
      &  iilc(:,:,:), iibc(:,:,:)      !< neighbor cells (array)
    INTEGER, CONTIGUOUS, POINTER :: &  !< Pointer to line and block indices of three
      &  iilnc(:,:,:), iibnc(:,:,:)    !< neighbor cells (array)
    INTEGER, CONTIGUOUS, POINTER :: &  !< Pointer to line and block indices (array)
      &  iidx(:,:,:), iblk(:,:,:)      !< of edges

    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend
    INTEGER  :: i_rlstart_e, i_rlend_e, i_rlstart_c, i_rlend_c
    INTEGER  :: je, jk, jb, jc         !< index of edge, vert level, block, cell

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: z_mflx_low,z_anti,z_mflx_anti_in,z_mflx_anti_out
!DIR$ ATTRIBUTES ALIGN :64 :: z_tracer_new_low,z_tracer_max,z_tracer_min
!DIR$ ATTRIBUTES ALIGN :64 :: r_p,r_m,z_min,z_max
#endif
  !-------------------------------------------------------------------------

    ! Set default values
    i_rlstart = grf_bdywidth_e
    i_rlend   = min_rledge_int - 1
    beta_fct  = 1._wp  ! the namelist default is 1.005, but it is passed to the limiter for the Miura3 scheme only

    ! Check for optional arguments
    IF (PRESENT(opt_rlstart)) i_rlstart = opt_rlstart
    IF (PRESENT(opt_rlend)) i_rlend = opt_rlend
    IF (PRESENT(opt_beta_fct)) beta_fct = opt_beta_fct

    r_beta_fct = 1._wp/beta_fct

    ! Set pointers to index-arrays

    ! line and block indices of two neighboring cells
    iilc => ptr_patch%edges%cell_idx
    iibc => ptr_patch%edges%cell_blk

    ! line and block indices of edges as seen from cells
    iidx => ptr_patch%cells%edge_idx
    iblk => ptr_patch%cells%edge_blk

    ! pointers to line and block indices of three neighbor cells
    iilnc => ptr_patch%cells%neighbor_idx
    iibnc => ptr_patch%cells%neighbor_blk

    !$ACC DATA CREATE(z_mflx_low, z_anti, z_mflx_anti_in, z_mflx_anti_out, r_m, r_p) &
    !$ACC   CREATE(z_tracer_new_low, z_tracer_max, z_tracer_min, z_min, z_max) &
    !$ACC   PRESENT(p_cc, p_mass_flx_e, p_rhodz_now, p_rhodz_new) PRESENT(p_mflx_tracer_h) &
    !$ACC   PRESENT(ptr_patch, ptr_int, iilc, iibc, iilnc, iibnc, iidx, iblk)

    IF (p_test_run) THEN
      !$ACC KERNELS PRESENT(r_p, r_m) ASYNC(1)
      r_p = 0._wp
      r_m = 0._wp
      !$ACC END KERNELS
    ENDIF

    !
    ! 1. Calculate low (first) order fluxes using the standard upwind scheme and the
    !    antidiffusive fluxes

    ! loop through all patch edges (and blocks)

    i_rlstart_e  = 5
    i_rlend_e    = min_rledge_int - 2
    i_startblk   = ptr_patch%edges%start_block(i_rlstart_e)
    i_endblk     = ptr_patch%edges%end_block(i_rlend_e)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart_e, i_rlend_e)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif
          !
          ! compute the first order upwind flux; notice
          ! that only the p_cc*p_vn value at cell edge is computed
          ! multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          !
          z_mflx_low(je,jk,jb) =  &
            &  laxfr_upflux(p_mass_flx_e(je,jk,jb),p_cc(iilc(je,jb,1),jk,iibc(je,jb,1)),p_cc(iilc(je,jb,2),jk,iibc(je,jb,2)))


          ! calculate antidiffusive flux for each edge
          ! only correct for i_rlend_e = min_rledge_int - 1, if p_mflx_tracer_h 
          ! is not synchronized. This is sufficient without iterative flux 
          ! correction which turned out to be overly expensive.
          z_anti(je,jk,jb)     = p_mflx_tracer_h(je,jk,jb) - z_mflx_low(je,jk,jb)


        END DO  ! end loop over edges
      END DO  ! end loop over levels
      !$ACC END PARALLEL

    END DO  ! end loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(i_rlstart_c,i_rlend_c,i_startblk,i_endblk)

    i_rlstart_c  = grf_bdywidth_c - 1
    i_rlend_c    = min_rlcell_int - 1
    i_startblk   = ptr_patch%cells%start_block(i_rlstart_c)
    i_endblk     = ptr_patch%cells%end_block(i_rlend_c)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_fluxdiv_c,z_mflx_anti_1,z_mflx_anti_2,z_mflx_anti_3 ) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,        &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_mflx_anti_1, z_mflx_anti_2, z_mflx_anti_3, z_fluxdiv_c)
#ifdef __LOOP_EXCHANGE
!DIR$ IVDEP,PREFERVECTOR
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          !
          ! 2. Define "antidiffusive" fluxes A(jc,jk,jb,je) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.

          z_mflx_anti_1 =                                                        &
            &     p_dtime * ptr_int%geofac_div(jc,1,jb) / p_rhodz_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,1),jk,iblk(jc,jb,1))

          z_mflx_anti_2 =                                                        &
            &     p_dtime * ptr_int%geofac_div(jc,2,jb) / p_rhodz_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,2),jk,iblk(jc,jb,2))

          z_mflx_anti_3 =                                                        &
            &     p_dtime * ptr_int%geofac_div(jc,3,jb) / p_rhodz_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,3),jk,iblk(jc,jb,3))

          ! Sum of all incoming antidiffusive fluxes into cell jc
          z_mflx_anti_in(jc,jk,jb) = -1._vp * (MIN(0._vp,z_mflx_anti_1) &
            &                                + MIN(0._vp,z_mflx_anti_2) &
            &                                + MIN(0._vp,z_mflx_anti_3) )

          ! Sum of all outgoing antidiffusive fluxes out of cell jc
          z_mflx_anti_out(jc,jk,jb) = MAX(0._vp,z_mflx_anti_1) &
            &                       + MAX(0._vp,z_mflx_anti_2) &
            &                       + MAX(0._vp,z_mflx_anti_3)

          !  compute also divergence of low order fluxes
          z_fluxdiv_c =  &
            & z_mflx_low(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            & z_mflx_low(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            & z_mflx_low(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)
!
! TODO:  the datum  z_mflx_low(iidx(jc,jb,3),jk,iblk(jc,jb,3)) yields differences later in z_tracer_new_low
!        The other entries do not cause a problem. 
!        Status 2015_09_07: problem still there in spite of corrections to mo_nonhydro_gpu_types,
!             both iidx(:,:,3) and iblk(:,:,3) possess problem.
!        Status 2015_09_22: this is related to the COLLAPSE directive mentioned above
!
          z_tracer_new_low(jc,jk,jb) =                        &
            &      ( p_cc(jc,jk,jb) * p_rhodz_now(jc,jk,jb)   &
            &      - p_dtime * z_fluxdiv_c )                  &
            &      / p_rhodz_new(jc,jk,jb)

          ! precalculate local maximum of current tracer value and low order
          ! updated value
          z_tracer_max(jc,jk,jb) = MAX(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))

          ! precalculate local minimum of current tracer value and low order
          ! updated value
          z_tracer_min(jc,jk,jb) = MIN(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      IF (ptr_patch%id > 1 .OR. l_limited_area) THEN

        ! Due to the lack of dynamic consistency between mass fluxes and cell mass changes
        ! in the boundary interpolation zone, the low-order advected tracer fields may be
        ! nonsense and therefore need artificial limitation

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          IF (ptr_patch%cells%refin_ctrl(jc,jb) == grf_bdywidth_c-1 .OR. &
              ptr_patch%cells%refin_ctrl(jc,jb) == grf_bdywidth_c) THEN
            DO jk = slev, elev
              z_tracer_new_low(jc,jk,jb) = MAX(0.9_wp*p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
              z_tracer_new_low(jc,jk,jb) = MIN(1.1_wp*p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
              z_tracer_max(jc,jk,jb) = MAX(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
              z_tracer_min(jc,jk,jb) = MIN(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
            ENDDO
          ENDIF
        ENDDO
        !$ACC END PARALLEL

      ENDIF

    ENDDO

!$OMP END DO

    ! Additional initialization of lateral boundary points is needed 
    ! for limited-area mode
    IF ( l_limited_area .AND. ptr_patch%id == 1 ) THEN

      i_startblk   = ptr_patch%cells%start_blk(1,1)
      i_endblk     = ptr_patch%cells%end_blk(grf_bdywidth_c-1,1)

      CALL init(r_m(:,:,i_startblk:i_endblk), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL init(r_p(:,:,i_startblk:i_endblk), lacc=.TRUE., opt_acc_async=.TRUE.)

!$OMP BARRIER

    ENDIF

    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.

    i_rlstart_c  = grf_bdywidth_c
    i_rlend_c    = min_rlcell_int
    i_startblk   = ptr_patch%cells%start_block(i_rlstart_c)
    i_endblk     = ptr_patch%cells%end_block(i_rlend_c)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_max,z_min) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          ! max value of cell and its neighbors
          ! also look back to previous time step
          z_max(jc,jk) = beta_fct * MAX( z_tracer_max(jc,jk,jb),               &
            &                 z_tracer_max(iilnc(jc,jb,1),jk,iibnc(jc,jb,1)),  &
            &                 z_tracer_max(iilnc(jc,jb,2),jk,iibnc(jc,jb,2)),  &
            &                 z_tracer_max(iilnc(jc,jb,3),jk,iibnc(jc,jb,3)) )

          ! min value of cell and its neighbors
          ! also look back to previous time step
          z_min(jc,jk) = r_beta_fct * MIN( z_tracer_min(jc,jk,jb),             &
            &                 z_tracer_min(iilnc(jc,jb,1),jk,iibnc(jc,jb,1)),  &
            &                 z_tracer_min(iilnc(jc,jb,2),jk,iibnc(jc,jb,2)),  &
            &                 z_tracer_min(iilnc(jc,jb,3),jk,iibnc(jc,jb,3)) )
        ENDDO
      ENDDO

      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of q
          r_m(jc,jk,jb) = (z_tracer_new_low(jc,jk,jb) - z_min(jc,jk))/ &
            &             (z_mflx_anti_out(jc,jk,jb) + dbl_eps)

          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of q
          r_p(jc,jk,jb) = (z_max(jc,jk) - z_tracer_new_low(jc,jk,jb))/ &
            &             (z_mflx_anti_in(jc,jk,jb) + dbl_eps)

        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Synchronize r_m and r_p and determine i_rlstart/i_rlend
    !
    CALL sync_patch_array_mult(SYNC_C1, ptr_patch, 2, r_m, r_p, opt_varname='r_m and r_p')

    !
    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !
    !    - at the end, compute new, limited fluxes which are then passed to 
    !      the main program. Note that p_mflx_tracer_h now denotes the 
    !      LIMITED flux.
    !

    i_startblk = ptr_patch%edges%start_block(i_rlstart)
    i_endblk   = ptr_patch%edges%end_block(i_rlend)



!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,r_frac,z_signum) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk,                &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      !
      ! compute final limited fluxes
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_signum, r_frac)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif

          z_signum = SIGN(1._wp,z_anti(je,jk,jb))

          ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
          ! but is computationally more efficient
          r_frac = 0.5_wp*( (1._wp+z_signum)*                &
             &     MIN(r_m(iilc(je,jb,1),jk,iibc(je,jb,1)),  &
             &         r_p(iilc(je,jb,2),jk,iibc(je,jb,2)))  &
             &     +  (1._wp-z_signum)*                      &
             &     MIN(r_m(iilc(je,jb,2),jk,iibc(je,jb,2)),  &
             &         r_p(iilc(je,jb,1),jk,iibc(je,jb,1)))  )

          ! Limited flux
          p_mflx_tracer_h(je,jk,jb) = z_mflx_low(je,jk,jb)               &
            &                       + MIN(1._wp,r_frac) * z_anti(je,jk,jb)

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE hflx_limiter_mo




  !-------------------------------------------------------------------------
  !>
  !! Positive definite flux limiter for horizontal advection
  !!
  !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
  !! Only outward fluxes are re-scaled, in order to maintain positive
  !! definiteness.
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
  !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
  !!   the cubed sphere grid. JCP, in press
  !!
  SUBROUTINE hflx_limiter_pd( ptr_patch, ptr_int, p_dtime, p_cc,        &
    &                         p_rhodz_now, p_mflx_tracer_h, slev, elev, &
    &                         opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_int_state), INTENT(IN) ::  &  !< pointer to data structure for
      &  ptr_int                         !< interpolation

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable
      &  p_cc(:,:,:)                 !< dim: (nproma,nlev,nblks_c)
                                     !< [kg kg^-1]

    REAL(wp), INTENT(IN) ::     &    !< density times cell thickness at timestep n
      &  p_rhodz_now(:,:,:)          !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &    !< time step [s]
      &  p_dtime

    REAL(wp), INTENT(INOUT) ::  &    !< calculated horizontal tracer mass flux
      &  p_mflx_tracer_h(:,:,:)      !< dim: (nproma,nlev,nblks_e)
                                     !< [kg m^-2 s^-1]

    INTEGER, INTENT(IN) ::      &    !< vertical start level
      &  slev

    INTEGER, INTENT(IN) ::      &    !< vertical end level
      &  elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    ! local variables
    !
    REAL(wp) :: &
      &  z_mflx1, z_mflx2, z_mflx3     !< tracer mass flux ( total mass crossing the edge )

    ! remark: single precision would be sufficient for r_m, but SP-sync is not yet available
    REAL(wp) ::                 &    !< fraction which must multiply all outgoing fluxes
      &  r_m(nproma,slev:elev,ptr_patch%nblks_c) !< of cell jc to guarantee
                                                      !< positive definiteness

    REAL(wp) :: z_signum                     !< sign of mass flux
                                             !< >0: out; <0: in
#if defined (__SX__) || defined ( _OPENACC )
    REAL(wp) :: p_m                          !< sum of fluxes out of cell jc
                                             !< [kg m^-3]
#else
    REAL(wp) :: p_m(nproma,slev:elev)
#endif
    INTEGER, CONTIGUOUS, POINTER :: &        !< Pointer to line and block indices of two
      &  iilc(:,:,:), iibc(:,:,:)            !< neighbor cells (array)

    INTEGER, CONTIGUOUS, POINTER :: &        !< Pointer to line and block indices (array)
      &  iidx(:,:,:), iblk(:,:,:)            !< of edges

    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_rlstart_c, i_rlend_c
    INTEGER  :: je, jk, jb, jc         !< index of edge, vert level, block, cell

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: p_m,r_m
#endif
  !-------------------------------------------------------------------------

    ! set default values
    i_rlstart = grf_bdywidth_e
    i_rlend   = min_rledge_int - 1

    ! Check for optional arguments
    IF (PRESENT(opt_rlstart)) i_rlstart = opt_rlstart
    IF (PRESENT(opt_rlend)) i_rlend = opt_rlend


    !
    ! Set pointers to index-arrays
    !
    ! line and block indices of two neighboring cells
    iilc => ptr_patch%edges%cell_idx
    iibc => ptr_patch%edges%cell_blk

    ! line and block indices of edges as seen from cells
    iidx => ptr_patch%cells%edge_idx
    iblk => ptr_patch%cells%edge_blk

    !$ACC ENTER DATA CREATE(r_m) ASYNC(1)
    !$ACC DATA PRESENT(p_cc, p_rhodz_now, p_mflx_tracer_h, r_m) &
    !$ACC   PRESENT(ptr_patch, ptr_int, iilc, iibc, iidx, iblk)

    IF (p_test_run) THEN
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
      r_m = 0._wp
      !$ACC END KERNELS
    ENDIF

!$OMP PARALLEL PRIVATE(i_rlstart_c,i_rlend_c,i_startblk,i_endblk)

    ! determine start row for cells from start row for edges.
    i_rlstart_c = get_startrow_c(startrow_e=i_rlstart)
    i_rlend_c   = min_rlcell_int

    ! Additional initialization of lateral boundary points is needed for limited-area mode
    IF ( l_limited_area .AND. ptr_patch%id == 1) THEN

      i_startblk   = ptr_patch%cells%start_block(1)
      i_endblk     = ptr_patch%cells%end_block(i_rlstart_c-1)

      CALL init(r_m(:,:,i_startblk:i_endblk), lacc=.TRUE.)
!$OMP BARRIER
    ENDIF

    i_startblk  = ptr_patch%cells%start_block(i_rlstart_c)
    i_endblk    = ptr_patch%cells%end_block(i_rlend_c)

    !
    ! 1. Reformulate all fluxes in terms of the total mass [kg m^-3]
    !    that crosses each of the CV-edges and store them in a cell-based structure.
    !
    !    z_mflx > 0: outward
    !    z_mflx < 0: inward
    !
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,p_m, &
!$OMP            z_mflx1,z_mflx2,z_mflx3) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,        &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_mflx1, z_mflx2, z_mflx3, p_m)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          z_mflx1 = ptr_int%geofac_div(jc,1,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,1),jk,iblk(jc,jb,1))
          z_mflx2 = ptr_int%geofac_div(jc,2,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,2),jk,iblk(jc,jb,2))
          z_mflx3 = ptr_int%geofac_div(jc,3,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,3),jk,iblk(jc,jb,3))

#if defined (__SX__) || defined ( _OPENACC )

          ! Sum of all outgoing fluxes out of cell jc
          p_m =  MAX(0._wp,z_mflx1) + MAX(0._wp,z_mflx2) + MAX(0._wp,z_mflx3)

          ! 2. fraction which must multiply all fluxes out of cell jc to guarantee no undershoot
          !    Nominator: maximum allowable decrease of \rho q
          r_m(jc,jk,jb) = MIN(1._wp, (p_cc(jc,jk,jb)*p_rhodz_now(jc,jk,jb)) / (p_m + dbl_eps) )

#else
          ! Sum of all outgoing fluxes out of cell jc
          p_m(jc,jk) = MAX(0._wp,z_mflx1) + &
            &          MAX(0._wp,z_mflx2) + &
            &          MAX(0._wp,z_mflx3)
        ENDDO
      ENDDO
      DO jk = slev, elev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          ! 2. fraction which must multiply all fluxes out of cell jc to guarantee
          !    no undershoot
          !    Nominator: maximum allowable decrease of \rho q
          r_m(jc,jk,jb) = MIN(1._wp, (p_cc(jc,jk,jb)*p_rhodz_now(jc,jk,jb)) &
            &                        /(p_m(jc,jk) + dbl_eps) )

#endif
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! synchronize r_m
    !
    IF(SIZE(r_m)/=0) CALL sync_patch_array(SYNC_C1,ptr_patch,r_m,opt_varname='r_m')

    !
    ! 3. Limit outward fluxes
    !    The inward ones remain untouched.
    !
    i_startblk = ptr_patch%edges%start_block(i_rlstart)
    i_endblk   = ptr_patch%edges%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,z_signum) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk,    &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_signum)
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif

          ! p_mflx_tracer_h > 0: flux directed from cell 1 -> 2
          ! p_mflx_tracer_h < 0: flux directed from cell 2 -> 1
          z_signum = SIGN(1._wp,p_mflx_tracer_h(je,jk,jb))

          p_mflx_tracer_h(je,jk,jb) = p_mflx_tracer_h(je,jk,jb) * 0.5_wp  &
            & *( (1._wp + z_signum) * r_m(iilc(je,jb,1),jk,iibc(je,jb,1)) &
            &   +(1._wp - z_signum) * r_m(iilc(je,jb,2),jk,iibc(je,jb,2)) )

        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT
    !$ACC END DATA
    !$ACC EXIT DATA DELETE(r_m)

  END SUBROUTINE hflx_limiter_pd

END MODULE mo_advection_hlimit

