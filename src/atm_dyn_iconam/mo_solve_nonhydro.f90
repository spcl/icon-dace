!
! This module contains the nonhydrostatic dynamical core for the triangular version
! Its routines were previously contained in mo_divergent_modes and mo_vector_operations
! but have been extracted for better memory efficiency
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

MODULE mo_solve_nonhydro

  USE mo_kind,                 ONLY: wp, vp
  USE mo_nonhydrostatic_config,ONLY: itime_scheme,iadv_rhotheta, igradp_method,             &
                                     kstart_moist, divdamp_order,                           &
                                     divdamp_fac, divdamp_fac2, divdamp_fac3, divdamp_fac4, &
                                     divdamp_z, divdamp_z2, divdamp_z3, divdamp_z4,         &
                                     divdamp_type, rayleigh_type, rhotheta_offctr,          &
                                     veladv_offctr, divdamp_fac_o2, kstart_dd3d, ndyn_substeps_var
  USE mo_dynamics_config,   ONLY: lmoist_thdyn, ldeepatmo
  USE mo_parallel_config,   ONLY: nproma, p_test_run, use_dycore_barrier, cpu_min_nproma
  USE mo_run_config,        ONLY: ltimer, timers_level, lvert_nest
  USE mo_model_domain,      ONLY: t_patch
  USE mo_grid_config,       ONLY: l_limited_area
  USE mo_gridref_config,    ONLY: grf_intmethod_e
  USE mo_interpol_config,   ONLY: nudge_max_coeff
  USE mo_intp_data_strc,    ONLY: t_int_state
  USE mo_intp,              ONLY: cells2verts_scalar
  USE mo_nonhydro_types,    ONLY: t_nh_state
  USE mo_physical_constants,ONLY: rd, cpd, cvd, cpv, grav, p0ref
  USE mo_math_gradients,    ONLY: grad_green_gauss_cell
  USE mo_velocity_advection,ONLY: velocity_tendencies
  USE mo_math_constants,    ONLY: dbl_eps
  USE mo_vertical_grid,     ONLY: nrdmax, nflat_gradp
  USE mo_init_vgrid,        ONLY: nflatlev
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,    ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int, &
    &                             min_rlcell, RAYLEIGH_CLASSIC, RAYLEIGH_KLEMP
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_sync,              ONLY: SYNC_E, SYNC_C, sync_patch_array,                             &
                                  sync_patch_array_mult, sync_patch_array_mult_mp
  USE mo_mpi,               ONLY: my_process_is_mpi_all_seq, work_mpi_barrier, i_am_accel_node
  USE mo_timer,             ONLY: timer_solve_nh, timer_barrier, timer_start, timer_stop,       &
                                  timer_solve_nh_cellcomp, timer_solve_nh_edgecomp,             &
                                  timer_solve_nh_vnupd, timer_solve_nh_vimpl, timer_solve_nh_exch
  USE mo_icon_comm_lib,     ONLY: icon_comm_sync
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_prepadv_types,     ONLY: t_prepare_adv
  USE mo_initicon_config,   ONLY: is_iau_active, iau_wgt_dyn
  USE mo_fortran_tools,     ONLY: init_zero_contiguous_dp, init_zero_contiguous_sp ! Import both for mixed prec.
#ifdef _OPENACC
  USE mo_mpi,               ONLY: my_process_is_work
#endif


  IMPLICIT NONE

  PRIVATE


  REAL(wp), PARAMETER :: rd_o_cvd = rd / cvd
  REAL(wp), PARAMETER :: cvd_o_rd = cvd / rd
  REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
  REAL(wp), PARAMETER :: grav_o_cpd = grav / cpd

  PUBLIC :: solve_nh

#ifdef _CRAYFTN
#define __CRAY_FTN_VERSION (_RELEASE_MAJOR * 100 + _RELEASE_MINOR)
#endif

  ! On the vectorizing DWD-NEC the diagnostics for the tendencies of the normal wind
  ! from terms xyz, ddt_vn_xyz, is disabled by default due to the fear that the
  ! conditional storage in conditionally allocated global fields is attempted even if
  ! the condition is not given and therefore the global field not allocated. If this
  ! happens, this would results in a corrupted memory.
  ! (Requested by G. Zaengl based on earlier problems with similar constructs.)
#ifndef __SX__
#define __ENABLE_DDT_VN_XYZ__
#endif

  CONTAINS


  !>
  !! solve_nh
  !!
  !! Main solver routine for nonhydrostatic dynamical core
  !!
  SUBROUTINE solve_nh (p_nh, p_patch, p_int, prep_adv, nnow, nnew, l_init, l_recompute, lsave_mflx, &
                       lprep_adv, lclean_mflx, idyn_timestep, jstep, dtime)

    TYPE(t_nh_state),    TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state),   TARGET, INTENT(IN)    :: p_int
    TYPE(t_patch),       TARGET, INTENT(INOUT) :: p_patch
    TYPE(t_prepare_adv), TARGET, INTENT(INOUT) :: prep_adv

    ! Initialization switch that has to be .TRUE. at the initial time step only (not for restart)
    LOGICAL,                   INTENT(IN)    :: l_init
    ! Switch to recompute velocity tendencies after a physics call irrespective of the time scheme option
    LOGICAL,                   INTENT(IN)    :: l_recompute
    ! Switch if mass flux needs to be saved for nest boundary interpolation tendency computation
    LOGICAL,                   INTENT(IN)    :: lsave_mflx
    ! Switch if preparations for tracer advection shall be computed
    LOGICAL,                   INTENT(IN)    :: lprep_adv
    ! Switch if mass fluxes computed for tracer advection need to be reinitialized
    LOGICAL,                   INTENT(IN)    :: lclean_mflx
    ! Counter of dynamics time step within a large time step (ranges from 1 to ndyn_substeps)
    INTEGER,                   INTENT(IN)    :: idyn_timestep
    ! Time step count since last boundary interpolation (ranges from 0 to 2*ndyn_substeps-1)
    INTEGER,                   INTENT(IN)    :: jstep
    ! Time levels
    INTEGER,                   INTENT(IN)    :: nnow, nnew
    ! Dynamics time step
    REAL(wp),                  INTENT(IN)    :: dtime

    ! Local variables
    INTEGER  :: jb, jk, jc, je, jks, jg
    INTEGER  :: nlev, nlevp1              !< number of full levels
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, ishift
    INTEGER  :: rl_start, rl_end, istep, ntl1, ntl2, nvar, nshift, nshift_total
    INTEGER  :: ic, ie, ilc0, ibc0, ikp1, ikp2

    REAL(wp) :: z_theta_v_fl_e  (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_vol_fl_e      (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_theta_v_e     (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_rho_e         (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_theta_v_v     (nproma,p_patch%nlev  ,p_patch%nblks_v), & ! used for iadv_rhotheta=1 only
                z_rho_v         (nproma,p_patch%nlev  ,p_patch%nblks_v)    ! used for iadv_rhotheta=1 only

    ! The data type vp (variable precision) is by default the same as wp but reduces
    ! to single precision when the __MIXED_PRECISION cpp flag is set at compile time
#ifdef __SWAPDIM
    REAL(vp) :: z_th_ddz_exner_c(nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_dexner_dz_c   (nproma,p_patch%nlev  ,p_patch%nblks_c,2), &
                z_vt_ie         (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_kin_hor_e     (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_exner_ex_pr   (nproma,p_patch%nlevp1,p_patch%nblks_c), & 
                z_gradh_exner   (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_rth_pr        (nproma,p_patch%nlev  ,p_patch%nblks_c,2), &
                z_grad_rth      (nproma,p_patch%nlev  ,p_patch%nblks_c,4), &
                z_w_concorr_me  (nproma,p_patch%nlev  ,p_patch%nblks_e)
#else
    REAL(vp) :: z_th_ddz_exner_c(nproma,p_patch%nlev,p_patch%nblks_c), &
                z_dexner_dz_c (2,nproma,p_patch%nlev,p_patch%nblks_c), &
                z_vt_ie         (nproma,p_patch%nlev,p_patch%nblks_e), &
                z_kin_hor_e     (nproma,p_patch%nlev,p_patch%nblks_e), &
                z_exner_ex_pr (nproma,p_patch%nlevp1,p_patch%nblks_c), & ! nlevp1 is intended here
                z_gradh_exner   (nproma,p_patch%nlev,p_patch%nblks_e), &
                z_rth_pr      (2,nproma,p_patch%nlev,p_patch%nblks_c), &
                z_grad_rth    (4,nproma,p_patch%nlev,p_patch%nblks_c), &
                z_w_concorr_me  (nproma,p_patch%nlev,p_patch%nblks_e)
#endif
    ! This field in addition has reversed index order (vertical first) for optimization
#ifdef __LOOP_EXCHANGE
    REAL(vp) :: z_graddiv_vn    (p_patch%nlev,nproma,p_patch%nblks_e)
#else
    REAL(vp) :: z_graddiv_vn    (nproma,p_patch%nlev,p_patch%nblks_e)
#endif

    REAL(wp) :: z_w_expl        (nproma,p_patch%nlevp1),          &
                z_vn_avg        (nproma,p_patch%nlev  ),          &
                z_mflx_top      (nproma,p_patch%nblks_c),         &
                z_contr_w_fl_l  (nproma,p_patch%nlevp1),          &
                z_contr_wq_fl_l (nproma,p_patch%nlevp1),          &
                z_rho_expl      (nproma,p_patch%nlev  ),          &
                z_exner_expl    (nproma,p_patch%nlev  )
    REAL(wp) :: z_theta_tavg_m1, z_theta_tavg, z_rho_tavg_m1, z_rho_tavg



    ! The data type vp (variable precision) is by default the same as wp but reduces
    ! to single precision when the __MIXED_PRECISION cpp flag is set at compile time

    ! TODO :  of these, fairly easy to scalarize:  z_theta_v_pr_ic
    REAL(vp) :: z_alpha         (nproma,p_patch%nlevp1),          &
                z_beta          (nproma,p_patch%nlev  ),          &
                z_beta_q        (nproma,p_patch%nlev  ),          &
                z_q             (nproma,p_patch%nlev  ),          &
                z_graddiv2_vn   (nproma,p_patch%nlev  ),          &
                z_theta_v_pr_ic (nproma,p_patch%nlevp1),          &
                z_exner_ic      (nproma,p_patch%nlevp1),          &
                z_w_concorr_mc  (nproma,p_patch%nlev  ),          &
                z_flxdiv_mass   (nproma,p_patch%nlev  ),          &
                z_flxdiv_theta  (nproma,p_patch%nlev  ),          &
                z_flxdiv_vol    (nproma,p_patch%nlev  ),          &
                z_hydro_corr    (nproma,p_patch%nblks_e)

    REAL(vp) :: z_a, z_b, z_c, z_g, z_gamma,      &
                z_w_backtraj, z_theta_v_pr_mc_m1, z_theta_v_pr_mc

#ifdef _OPENACC
    REAL(vp) :: z_w_concorr_mc_m0, z_w_concorr_mc_m1, z_w_concorr_mc_m2
#endif

    REAL(wp) :: z_theta1, z_theta2, wgt_nnow_vel, wgt_nnew_vel,     &
               dt_shift, wgt_nnow_rth, wgt_nnew_rth, dthalf,        &
               r_nsubsteps, r_dtimensubsteps, scal_divdamp_o2,      &
               alin, dz32, df32, dz42, df42, bqdr, aqdr,            &
               zf, dzlin, dzqdr
    ! time shifts for linear interpolation of nest UBC
    REAL(wp) :: dt_linintp_ubc, dt_linintp_ubc_nnow, dt_linintp_ubc_nnew
    REAL(wp) :: z_raylfac(nrdmax(p_patch%id))
    REAL(wp) :: z_ntdistv_bary_1, distv_bary_1, z_ntdistv_bary_2, distv_bary_2

    REAL(wp), DIMENSION(p_patch%nlev) :: scal_divdamp, bdy_divdamp, enh_divdamp_fac
    REAL(vp) :: z_dwdz_dd(nproma,kstart_dd3d(p_patch%id):p_patch%nlev,p_patch%nblks_c)

    ! Local variables for normal wind tendencies and differentials
    REAL(wp) :: z_ddt_vn_dyn, z_ddt_vn_apc, z_ddt_vn_cor, &
      &         z_ddt_vn_pgr, z_ddt_vn_ray, z_d_vn_dmp, z_d_vn_iau

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_theta_v_fl_e,z_theta_v_e,z_rho_e
!DIR$ ATTRIBUTES ALIGN : 64 :: z_theta_v_v,z_rho_v,z_dwdz_dd
!DIR$ ATTRIBUTES ALIGN : 64 :: z_th_ddz_exner_c,z_dexner_dz_c,z_vt_ie,z_kin_hor_e
!DIR$ ATTRIBUTES ALIGN : 64 :: z_exner_ex_pr,z_gradh_exner,z_rth_pr,z_grad_rth
!DIR$ ATTRIBUTES ALIGN : 64 :: z_w_concorr_me,z_graddiv_vn,z_w_expl
!DIR$ ATTRIBUTES ALIGN : 64 :: z_vn_avg,z_mflx_top,z_contr_w_fl_l,z_contr_wq_fl_l,z_rho_expl
!DIR$ ATTRIBUTES ALIGN : 64 :: z_exner_expl,z_alpha,z_beta,z_beta_q,z_q,z_graddiv2_vn
!DIR$ ATTRIBUTES ALIGN : 64 :: z_theta_v_pr_ic,z_exner_ic,z_w_concorr_mc
!DIR$ ATTRIBUTES ALIGN : 64 :: z_flxdiv_mass,z_flxdiv_theta,z_flxdiv_vol,z_hydro_corr
!DIR$ ATTRIBUTES ALIGN : 64 :: z_raylfac,scal_divdamp,bdy_divdamp,enh_divdamp_fac
#endif

    INTEGER :: nproma_gradp, nblks_gradp, npromz_gradp, nlen_gradp, jk_start
    LOGICAL :: lvn_only, lvn_pos, lapply_moist_thdyn

    ! Local variables to control vertical nesting
    LOGICAL :: l_vert_nested, l_child_vertnest

    ! Pointers
    INTEGER, POINTER, CONTIGUOUS :: &
      ! to cell indices
      icidx(:,:,:), icblk(:,:,:), &
      ! to edge indices
      ieidx(:,:,:), ieblk(:,:,:), &
      ! to vertex indices
      ividx(:,:,:), ivblk(:,:,:), &
      ! to vertical neighbor indices for pressure gradient computation
      ikidx(:,:,:,:),             &
      ! to quad edge indices
      iqidx(:,:,:), iqblk(:,:,:), &
      ! for igradp_method = 3
      iplev(:), ipeidx(:), ipeblk(:)

#ifdef __SX__
      REAL(wp) :: z_rho_tavg_m1_v(nproma), z_theta_tavg_m1_v(nproma)
      REAL(vp) :: z_theta_v_pr_mc_m1_v(nproma)
#endif
    !-------------------------------------------------------------------
    IF (use_dycore_barrier) THEN
      CALL timer_start(timer_barrier)
      CALL work_mpi_barrier()
      CALL timer_stop(timer_barrier)
    ENDIF
    !-------------------------------------------------------------------

    jg = p_patch%id

    IF (lvert_nest .AND. (p_patch%nshift_total > 0)) THEN
      l_vert_nested = .TRUE.
      nshift_total  = p_patch%nshift_total
    ELSE
      l_vert_nested = .FALSE.
      nshift_total  = 0
    ENDIF
    IF (lvert_nest .AND. p_patch%n_childdom > 0 .AND.              &
      (p_patch%nshift_child > 0 .OR. p_patch%nshift_total > 0)) THEN
      l_child_vertnest = .TRUE.
      nshift = p_patch%nshift_child + 1
    ELSE
      l_child_vertnest = .FALSE.
      nshift = 0
    ENDIF
    dthalf  = 0.5_wp*dtime

    IF (ltimer) CALL timer_start(timer_solve_nh)

    ! Inverse value of ndyn_substeps for tracer advection precomputations
    r_nsubsteps = 1._wp/REAL(ndyn_substeps_var(jg),wp)

    ! Inverse value of dtime * ndyn_substeps_var
    r_dtimensubsteps = 1._wp/(dtime*REAL(ndyn_substeps_var(jg),wp))

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Set pointers to neighbor cells
    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk

    ! Set pointers to neighbor edges
    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    ! Set pointers to vertices of an edge
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    ! Set pointer to vertical neighbor indices for pressure gradient
    ikidx => p_nh%metrics%vertidx_gradp

    ! Set pointers to quad edges
    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk

    ! DA: moved from below to here to get into the same ACC data section
    iplev  => p_nh%metrics%pg_vertidx
    ipeidx => p_nh%metrics%pg_edgeidx
    ipeblk => p_nh%metrics%pg_edgeblk

    
    ! Precompute Rayleigh damping factor
    DO jk = 2, nrdmax(jg)
       z_raylfac(jk) = 1.0_wp/(1.0_wp+dtime*p_nh%metrics%rayleigh_w(jk))
    ENDDO

    ! Fourth-order divergence damping
    !
    ! The divergence damping factor enh_divdamp_fac is defined as a profile in height z
    ! above sea level with 4 height sections:
    !
    ! enh_divdamp_fac(z) = divdamp_fac                                              !               z <= divdamp_z
    ! enh_divdamp_fac(z) = divdamp_fac  + (z-divdamp_z )* alin                      ! divdamp_z  <= z <= divdamp_z2
    ! enh_divdamp_fac(z) = divdamp_fac2 + (z-divdamp_z2)*(aqdr+(z-divdamp_z2)*bqdr) ! divdamp_z2 <= z <= divdamp_z4
    ! enh_divdamp_fac(z) = divdamp_fac4                                             ! divdamp_z4 <= z
    !
    alin = (divdamp_fac2-divdamp_fac)/(divdamp_z2-divdamp_z)
    !
    df32 = divdamp_fac3-divdamp_fac2; dz32 = divdamp_z3-divdamp_z2
    df42 = divdamp_fac4-divdamp_fac2; dz42 = divdamp_z4-divdamp_z2
    !
    bqdr = (df42*dz32-df32*dz42)/(dz32*dz42*(dz42-dz32))
    aqdr = df32/dz32-bqdr*dz32
    !
    DO jk = 1, nlev
      jks = jk + nshift_total
      zf = 0.5_wp*(vct_a(jks)+vct_a(jks+1))
      dzlin = MIN(divdamp_z2-divdamp_z ,MAX(0._wp,zf-divdamp_z ))
      dzqdr = MIN(divdamp_z4-divdamp_z2,MAX(0._wp,zf-divdamp_z2))
      !
      IF (divdamp_order == 24) THEN
        enh_divdamp_fac(jk) = MAX( 0._wp, divdamp_fac + dzlin*alin + dzqdr*(aqdr+dzqdr*bqdr) - 0.25_wp*divdamp_fac_o2 )
      ELSE
        enh_divdamp_fac(jk) =             divdamp_fac + dzlin*alin + dzqdr*(aqdr+dzqdr*bqdr)
      ENDIF
    ENDDO

    scal_divdamp(:) = - enh_divdamp_fac(:) * p_patch%geometry_info%mean_cell_area**2

    ! Time increment for backward-shifting of lateral boundary mass flux
    dt_shift = dtime*REAL(2*ndyn_substeps_var(jg)-1,wp)/2._wp    ! == dt_phy - 0.5*dtime

    ! Time increment for linear interpolation of nest UBC.
    ! The linear interpolation is of the form
    ! \phi(t) = \phi0 + (t-t0)*dphi/dt, with t=(jstep+0.5)*dtime, and t0=dt_phy
    !
    ! dt_linintp_ubc == (t-t0)
    dt_linintp_ubc      = jstep*dtime - dt_shift ! valid for center of current time step
    dt_linintp_ubc_nnow = dt_linintp_ubc - 0.5_wp*dtime
    dt_linintp_ubc_nnew = dt_linintp_ubc + 0.5_wp*dtime

    ! Coefficient for reduced fourth-order divergence damping along nest boundaries
    bdy_divdamp(:) = 0.75_wp/(nudge_max_coeff + dbl_eps)*ABS(scal_divdamp(:))

    !$ACC DATA CREATE(z_kin_hor_e, z_vt_ie, z_w_concorr_me, z_theta_v_fl_e) &
    !$ACC   CREATE(z_dexner_dz_c, z_exner_ex_pr, z_gradh_exner, z_rth_pr, z_grad_rth) &
    !$ACC   CREATE(z_theta_v_pr_ic, z_th_ddz_exner_c, z_w_concorr_mc) &
    !$ACC   CREATE(z_vn_avg, z_rho_e, z_theta_v_e, z_dwdz_dd, z_mflx_top) &
    !$ACC   CREATE(z_exner_ic, z_alpha, z_beta, z_beta_q, z_q, z_contr_w_fl_l, z_contr_wq_fl_l, z_exner_expl) &
    !$ACC   CREATE(z_flxdiv_mass, z_flxdiv_theta, z_flxdiv_vol, z_rho_expl, z_w_expl, z_vol_fl_e) &
    !$ACC   CREATE(z_rho_v, z_theta_v_v, z_graddiv_vn, z_hydro_corr, z_graddiv2_vn) &
    !$ACC   COPYIN(nflatlev, nflat_gradp, kstart_dd3d, kstart_moist, nrdmax) &
    !$ACC   COPYIN(z_raylfac, ndyn_substeps_var, scal_divdamp, bdy_divdamp) &
    !$ACC   PRESENT(prep_adv, p_int, p_patch, p_nh) &
    !$ACC   PRESENT(icidx, icblk, ividx, ivblk, ieidx, ieblk, ikidx, iqidx, iqblk) &
    !$ACC   PRESENT(ipeidx, ipeblk, iplev) &
    !$ACC   IF(i_am_accel_node)

    ! scaling factor for second-order divergence damping: divdamp_fac_o2*delta_x**2
    ! delta_x**2 is approximated by the mean cell area
    scal_divdamp_o2 = divdamp_fac_o2 * p_patch%geometry_info%mean_cell_area


    IF (p_test_run) THEN
      !$ACC KERNELS IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
      z_rho_e     = 0._wp
      z_theta_v_e = 0._wp
      z_dwdz_dd   = 0._wp
      z_graddiv_vn= 0._wp
      !$ACC END KERNELS
    ENDIF

    ! Set time levels of ddt_adv fields for call to velocity_tendencies
    IF (itime_scheme >= 4) THEN ! Velocity advection averaging nnow and nnew tendencies
      ntl1 = nnow
      ntl2 = nnew
    ELSE                        ! Velocity advection is taken at nnew only
      ntl1 = 1
      ntl2 = 1
    ENDIF

    ! Weighting coefficients for velocity advection if tendency averaging is used
    ! The off-centering specified here turned out to be beneficial to numerical
    ! stability in extreme situations
    wgt_nnow_vel = 0.5_wp - veladv_offctr ! default value for veladv_offctr is 0.25
    wgt_nnew_vel = 0.5_wp + veladv_offctr

    ! Weighting coefficients for rho and theta at interface levels in the corrector step
    ! This empirically determined weighting minimizes the vertical wind off-centering
    ! needed for numerical stability of vertical sound wave propagation
    wgt_nnew_rth = 0.5_wp + rhotheta_offctr ! default value for rhotheta_offctr is -0.1
    wgt_nnow_rth = 1._wp - wgt_nnew_rth

!$NEC sparse
    DO istep = 1, 2

      IF (istep == 1) THEN ! predictor step
        IF (itime_scheme >= 6 .OR. l_init .OR. l_recompute) THEN
          IF (itime_scheme < 6 .AND. .NOT. l_init) THEN
            lvn_only = .TRUE. ! Recompute only vn tendency
          ELSE
            lvn_only = .FALSE.
          ENDIF
          CALL velocity_tendencies(p_nh%prog(nnow),p_patch,p_int,p_nh%metrics,p_nh%diag,z_w_concorr_me, &
            z_kin_hor_e,z_vt_ie,ntl1,istep,lvn_only,dtime,dt_linintp_ubc_nnow,ldeepatmo)
        ENDIF
        nvar = nnow
      ELSE                 ! corrector step
        lvn_only = .FALSE.
        CALL velocity_tendencies(p_nh%prog(nnew),p_patch,p_int,p_nh%metrics,p_nh%diag,z_w_concorr_me, &
          z_kin_hor_e,z_vt_ie,ntl2,istep,lvn_only,dtime,dt_linintp_ubc_nnew,ldeepatmo)
        nvar = nnew
      ENDIF

      lapply_moist_thdyn = (istep == 2) .AND. lmoist_thdyn

      ! Preparations for igradp_method = 3/5 (reformulated extrapolation below the ground)
      IF (istep == 1 .AND. (igradp_method == 3 .OR. igradp_method == 5)) THEN

        nproma_gradp = cpu_min_nproma(nproma,256)
        nblks_gradp  = INT(p_nh%metrics%pg_listdim/nproma_gradp)
        npromz_gradp = MOD(p_nh%metrics%pg_listdim,nproma_gradp)
        IF (npromz_gradp > 0) THEN
          nblks_gradp = nblks_gradp + 1
        ELSE
          npromz_gradp = nproma_gradp
        ENDIF

      ENDIF

      IF (timers_level > 5) CALL timer_start(timer_solve_nh_cellcomp)

      ! Computations on mass points
!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

      rl_start = 3
      IF (istep == 1) THEN
        rl_end = min_rlcell_int - 1
      ELSE ! halo points are not needed in step 2
        rl_end = min_rlcell_int
      ENDIF

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

      ! initialize nest boundary points of z_rth_pr with zero
      IF (istep == 1 .AND. (jg > 1 .OR. l_limited_area)) THEN

#ifdef __SWAPDIM
#ifdef __MIXED_PRECISION
          CALL init_zero_contiguous_sp(z_rth_pr(1,1,1,1), nproma*nlev*i_startblk, opt_acc_async=.TRUE., &
            & lacc=i_am_accel_node)
          CALL init_zero_contiguous_sp(z_rth_pr(1,1,1,2), nproma*nlev*i_startblk, opt_acc_async=.TRUE., &
            & lacc=i_am_accel_node)
#else
          CALL init_zero_contiguous_dp(z_rth_pr(1,1,1,1), nproma*nlev*i_startblk, opt_acc_async=.TRUE., &
            & lacc=i_am_accel_node)
          CALL init_zero_contiguous_dp(z_rth_pr(1,1,1,2), nproma*nlev*i_startblk, opt_acc_async=.TRUE., &
            & lacc=i_am_accel_node)
#endif
#else
#ifdef __MIXED_PRECISION
          CALL init_zero_contiguous_sp(z_rth_pr(1,1,1,1), 2*nproma*nlev*i_startblk, opt_acc_async=.TRUE., &
            & lacc=i_am_accel_node)
#else
          CALL init_zero_contiguous_dp(z_rth_pr(1,1,1,1), 2*nproma*nlev*i_startblk, opt_acc_async=.TRUE., &
            & lacc=i_am_accel_node)
#endif
#endif

!$OMP BARRIER
      ENDIF

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,z_exner_ic,z_theta_v_pr_ic,z_w_backtraj,&
!$OMP            z_theta_v_pr_mc_m1,z_theta_v_pr_mc,z_rho_tavg_m1,z_rho_tavg, &
#ifdef __SX__
!$OMP            z_rho_tavg_m1_v,z_theta_tavg_m1_v,z_theta_v_pr_mc_m1_v, &
#endif
!$OMP            z_theta_tavg_m1,z_theta_tavg) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, rl_start, rl_end)

        IF (istep == 1) THEN ! to be executed in predictor step only

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              ! temporally extrapolated perturbation Exner pressure (used for horizontal gradients only)
              z_exner_ex_pr(jc,jk,jb) = (1._wp + p_nh%metrics%exner_exfac(jc,jk,jb)) *    &
                (p_nh%prog(nnow)%exner(jc,jk,jb) - p_nh%metrics%exner_ref_mc(jc,jk,jb)) - &
                 p_nh%metrics%exner_exfac(jc,jk,jb) * p_nh%diag%exner_pr(jc,jk,jb)

              ! non-extrapolated perturbation Exner pressure, saved in exner_pr for the next time step
              p_nh%diag%exner_pr(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb) - &
                                              p_nh%metrics%exner_ref_mc(jc,jk,jb)

            ENDDO
          ENDDO
          !$ACC END PARALLEL

          ! The purpose of the extra level of exner_pr is to simplify coding for
          ! igradp_method=4/5. It is multiplied with zero and thus actually not used
          !$ACC KERNELS IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          z_exner_ex_pr(:,nlevp1,jb) = 0._wp
          !$ACC END KERNELS


          IF (igradp_method <= 3) THEN
            ! Perturbation Exner pressure on bottom half level
!DIR$ IVDEP
            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR
            DO jc = i_startidx, i_endidx
              z_exner_ic(jc,nlevp1) =                                         &
                p_nh%metrics%wgtfacq_c(jc,1,jb)*z_exner_ex_pr(jc,nlev  ,jb) + &
                p_nh%metrics%wgtfacq_c(jc,2,jb)*z_exner_ex_pr(jc,nlev-1,jb) + &
                p_nh%metrics%wgtfacq_c(jc,3,jb)*z_exner_ex_pr(jc,nlev-2,jb)
            ENDDO
            !$ACC END PARALLEL

! WS: moved full z_exner_ic calculation here to avoid OpenACC dependency on jk+1 below
!     possibly GZ will want to consider the cache ramifications of this change for CPU
            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR TILE(32, 4)
            DO jk = nlev, MAX(2,nflatlev(jg)), -1
!DIR$ IVDEP
              DO jc = i_startidx, i_endidx
                ! Exner pressure on remaining half levels for metric correction term
                z_exner_ic(jc,jk) =                                                    &
                         p_nh%metrics%wgtfac_c(jc,jk,jb) *z_exner_ex_pr(jc,jk  ,jb) +  &
                  (1._vp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_exner_ex_pr(jc,jk-1,jb)
              ENDDO
            ENDDO
            !$ACC END PARALLEL

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR TILE(32, 4)
            DO jk = nlev, MAX(2,nflatlev(jg)), -1
!DIR$ IVDEP
              DO jc = i_startidx, i_endidx

                ! First vertical derivative of perturbation Exner pressure
#ifdef __SWAPDIM
                z_dexner_dz_c(jc,jk,jb,1) =                     &
#else
                z_dexner_dz_c(1,jc,jk,jb) =                     &
#endif
                  (z_exner_ic(jc,jk) - z_exner_ic(jc,jk+1)) *   &
                  p_nh%metrics%inv_ddqz_z_full(jc,jk,jb)
              ENDDO
            ENDDO
            !$ACC END PARALLEL

            IF (nflatlev(jg) == 1) THEN
              ! Perturbation Exner pressure on top half level
!DIR$ IVDEP
              !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
              !$ACC LOOP GANG VECTOR
              DO jc = i_startidx, i_endidx
                z_exner_ic(jc,1) =                                          &
                  p_nh%metrics%wgtfacq1_c(jc,1,jb)*z_exner_ex_pr(jc,1,jb) + &
                  p_nh%metrics%wgtfacq1_c(jc,2,jb)*z_exner_ex_pr(jc,2,jb) + &
                  p_nh%metrics%wgtfacq1_c(jc,3,jb)*z_exner_ex_pr(jc,3,jb)

                ! First vertical derivative of perturbation Exner pressure
#ifdef __SWAPDIM
                z_dexner_dz_c(jc,1,jb,1) =                    &
#else
                z_dexner_dz_c(1,jc,1,jb) =                    &
#endif
                  (z_exner_ic(jc,1) - z_exner_ic(jc,2)) *   &
                  p_nh%metrics%inv_ddqz_z_full(jc,1,jb)
              ENDDO
              !$ACC END PARALLEL
            ENDIF

          ENDIF

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
#ifdef __SWAPDIM
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            z_rth_pr(jc,1,jb,1) = p_nh%prog(nnow)%rho(jc,1,jb) - &
              p_nh%metrics%rho_ref_mc(jc,1,jb)
            z_rth_pr(jc,1,jb,2) = p_nh%prog(nnow)%theta_v(jc,1,jb) - &
              p_nh%metrics%theta_ref_mc(jc,1,jb)
          ENDDO
#else
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            z_rth_pr(1,jc,1,jb) =  p_nh%prog(nnow)%rho(jc,1,jb) - &
              p_nh%metrics%rho_ref_mc(jc,1,jb)
            z_rth_pr(2,jc,1,jb) =  p_nh%prog(nnow)%theta_v(jc,1,jb) - &
              p_nh%metrics%theta_ref_mc(jc,1,jb)
          ENDDO
#endif
          !$ACC END PARALLEL

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR TILE(32, 4)
          DO jk = 2, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              ! density at interface levels for vertical flux divergence computation
              p_nh%diag%rho_ic(jc,jk,jb) = p_nh%metrics%wgtfac_c(jc,jk,jb) *p_nh%prog(nnow)%rho(jc,jk  ,jb) + &
                                    (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*p_nh%prog(nnow)%rho(jc,jk-1,jb)

              ! perturbation density and virtual potential temperature at main levels for horizontal flux divergence term
              ! (needed in the predictor step only)
#ifdef __SWAPDIM
              z_rth_pr(jc,jk,jb,1) =  p_nh%prog(nnow)%rho(jc,jk,jb)     - p_nh%metrics%rho_ref_mc(jc,jk,jb)
              z_rth_pr(jc,jk,jb,2) =  p_nh%prog(nnow)%theta_v(jc,jk,jb) - p_nh%metrics%theta_ref_mc(jc,jk,jb)
#else
              z_rth_pr(1,jc,jk,jb) =  p_nh%prog(nnow)%rho(jc,jk,jb)     - p_nh%metrics%rho_ref_mc(jc,jk,jb)
              z_rth_pr(2,jc,jk,jb) =  p_nh%prog(nnow)%theta_v(jc,jk,jb) - p_nh%metrics%theta_ref_mc(jc,jk,jb)
#endif
#ifdef _OPENACC
            ENDDO
          ENDDO
          !$ACC END PARALLEL

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 2, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
#endif

              ! perturbation virtual potential temperature at interface levels
#ifdef __SWAPDIM
              z_theta_v_pr_ic(jc,jk) =                                           &
                       p_nh%metrics%wgtfac_c(jc,jk,jb) *z_rth_pr(jc,jk  ,jb,2) + &
                (1._vp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_rth_pr(jc,jk-1,jb,2)
#else
              z_theta_v_pr_ic(jc,jk) =                                           &
                       p_nh%metrics%wgtfac_c(jc,jk,jb) *z_rth_pr(2,jc,jk  ,jb) + &
                (1._vp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_rth_pr(2,jc,jk-1,jb)
#endif
              ! virtual potential temperature at interface levels
              p_nh%diag%theta_v_ic(jc,jk,jb) =                                                &
                       p_nh%metrics%wgtfac_c(jc,jk,jb) *p_nh%prog(nnow)%theta_v(jc,jk  ,jb) + &
                (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*p_nh%prog(nnow)%theta_v(jc,jk-1,jb)

              ! vertical pressure gradient * theta_v
              z_th_ddz_exner_c(jc,jk,jb) = p_nh%metrics%vwind_expl_wgt(jc,jb)* &
                p_nh%diag%theta_v_ic(jc,jk,jb) * (p_nh%diag%exner_pr(jc,jk-1,jb)-      &
                p_nh%diag%exner_pr(jc,jk,jb)) / p_nh%metrics%ddqz_z_half(jc,jk,jb) +   &
                z_theta_v_pr_ic(jc,jk)*p_nh%metrics%d_exner_dz_ref_ic(jc,jk,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL

        ELSE  ! istep = 2 - in this step, an upwind-biased discretization is used for rho_ic and theta_v_ic
          ! in order to reduce the numerical dispersion errors
#ifdef __SX__
          ! precompute values for jk = 1 which are previous values in first iteration of jk compute loop
          jk = 2
            DO jc = i_startidx, i_endidx
              z_rho_tavg_m1_v(jc) = wgt_nnow_rth*p_nh%prog(nnow)%rho(jc,jk-1,jb) + &
                              wgt_nnew_rth*p_nh%prog(nvar)%rho(jc,jk-1,jb)
              z_theta_tavg_m1_v(jc) = wgt_nnow_rth*p_nh%prog(nnow)%theta_v(jc,jk-1,jb) + &
                                wgt_nnew_rth*p_nh%prog(nvar)%theta_v(jc,jk-1,jb)
              z_theta_v_pr_mc_m1_v(jc)  = z_theta_tavg_m1_v(jc) - p_nh%metrics%theta_ref_mc(jc,jk-1,jb)
            ENDDO
#endif

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR TILE(128, 1) &
          !$ACC   PRIVATE(z_w_backtraj, z_rho_tavg_m1, z_theta_tavg_m1, z_rho_tavg) &
          !$ACC   PRIVATE(z_theta_tavg, z_theta_v_pr_mc_m1, z_theta_v_pr_mc)
          DO jk = 2, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              ! backward trajectory - use w(nnew) in order to be at the same time level as w_concorr
              z_w_backtraj = - (p_nh%prog(nnew)%w(jc,jk,jb) - p_nh%diag%w_concorr_c(jc,jk,jb)) * &
                dtime*0.5_wp/p_nh%metrics%ddqz_z_half(jc,jk,jb)

              ! temporally averaged density and virtual potential temperature depending on rhotheta_offctr
              ! (see pre-computation above)
#ifndef __SX__
              z_rho_tavg_m1 = wgt_nnow_rth*p_nh%prog(nnow)%rho(jc,jk-1,jb) + &
                              wgt_nnew_rth*p_nh%prog(nvar)%rho(jc,jk-1,jb)
              z_theta_tavg_m1 = wgt_nnow_rth*p_nh%prog(nnow)%theta_v(jc,jk-1,jb) + &
                                wgt_nnew_rth*p_nh%prog(nvar)%theta_v(jc,jk-1,jb)
#else
              z_rho_tavg_m1   = z_rho_tavg_m1_v(jc)
              z_theta_tavg_m1 = z_theta_tavg_m1_v(jc)
#endif

              z_rho_tavg = wgt_nnow_rth*p_nh%prog(nnow)%rho(jc,jk,jb) + &
                           wgt_nnew_rth*p_nh%prog(nvar)%rho(jc,jk,jb)
              z_theta_tavg = wgt_nnow_rth*p_nh%prog(nnow)%theta_v(jc,jk,jb) + &
                             wgt_nnew_rth*p_nh%prog(nvar)%theta_v(jc,jk,jb)

              ! density at interface levels for vertical flux divergence computation
              p_nh%diag%rho_ic(jc,jk,jb) = p_nh%metrics%wgtfac_c(jc,jk,jb) *z_rho_tavg    + &
                                    (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_rho_tavg_m1 + &
                z_w_backtraj*(z_rho_tavg_m1-z_rho_tavg)

              ! perturbation virtual potential temperature at main levels
#ifndef __SX__
              z_theta_v_pr_mc_m1  = z_theta_tavg_m1 - p_nh%metrics%theta_ref_mc(jc,jk-1,jb)
#else
              z_theta_v_pr_mc_m1 = z_theta_v_pr_mc_m1_v(jc)
#endif
              z_theta_v_pr_mc     = z_theta_tavg    - p_nh%metrics%theta_ref_mc(jc,jk,jb)

              ! perturbation virtual potential temperature at interface levels
              z_theta_v_pr_ic(jc,jk) =                                       &
                       p_nh%metrics%wgtfac_c(jc,jk,jb) *z_theta_v_pr_mc +    &
                (1._vp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_theta_v_pr_mc_m1

              ! virtual potential temperature at interface levels
              p_nh%diag%theta_v_ic(jc,jk,jb) = p_nh%metrics%wgtfac_c(jc,jk,jb) *z_theta_tavg    +  &
                                        (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_theta_tavg_m1 +  &
                z_w_backtraj*(z_theta_tavg_m1-z_theta_tavg)

              ! vertical pressure gradient * theta_v
              z_th_ddz_exner_c(jc,jk,jb) = p_nh%metrics%vwind_expl_wgt(jc,jb)* &
                p_nh%diag%theta_v_ic(jc,jk,jb) * (p_nh%diag%exner_pr(jc,jk-1,jb)-      &
                p_nh%diag%exner_pr(jc,jk,jb)) / p_nh%metrics%ddqz_z_half(jc,jk,jb) +   &
                z_theta_v_pr_ic(jc,jk)*p_nh%metrics%d_exner_dz_ref_ic(jc,jk,jb)

#ifdef __SX__
              ! save current values as previous values for next iteration
              z_rho_tavg_m1_v(jc) = z_rho_tavg
              z_theta_tavg_m1_v(jc) = z_theta_tavg
              z_theta_v_pr_mc_m1_v(jc) = z_theta_v_pr_mc
#endif

            ENDDO
          ENDDO
          !$ACC END PARALLEL

        ENDIF ! istep = 1/2


        IF (istep == 1) THEN

          ! Perturbation theta at top and surface levels
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
!DIR$ IVDEP
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            z_theta_v_pr_ic(jc,1)      = 0._wp
            z_theta_v_pr_ic(jc,nlevp1) =                                   &
#ifdef __SWAPDIM
              p_nh%metrics%wgtfacq_c(jc,1,jb)*z_rth_pr(jc,nlev  ,jb,2) +     &
              p_nh%metrics%wgtfacq_c(jc,2,jb)*z_rth_pr(jc,nlev-1,jb,2) +   &
              p_nh%metrics%wgtfacq_c(jc,3,jb)*z_rth_pr(jc,nlev-2,jb,2)
#else
              p_nh%metrics%wgtfacq_c(jc,1,jb)*z_rth_pr(2,jc,nlev  ,jb) +     &
              p_nh%metrics%wgtfacq_c(jc,2,jb)*z_rth_pr(2,jc,nlev-1,jb) +   &
              p_nh%metrics%wgtfacq_c(jc,3,jb)*z_rth_pr(2,jc,nlev-2,jb)
#endif
            p_nh%diag%theta_v_ic(jc,nlevp1,jb) =                                  &
              p_nh%metrics%theta_ref_ic(jc,nlevp1,jb) + z_theta_v_pr_ic(jc,nlevp1)
          ENDDO
          !$ACC END PARALLEL

          IF (igradp_method <= 3) THEN

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR TILE(32, 4)
            DO jk = nflat_gradp(jg), nlev
!DIR$ IVDEP
              DO jc = i_startidx, i_endidx
                ! Second vertical derivative of perturbation Exner pressure (hydrostatic approximation)
#ifdef __SWAPDIM
                z_dexner_dz_c(jc,jk,jb,2) = -0.5_vp *                              &
                  ((z_theta_v_pr_ic(jc,jk) - z_theta_v_pr_ic(jc,jk+1)) *           &
                  p_nh%metrics%d2dexdz2_fac1_mc(jc,jk,jb) + z_rth_pr(jc,jk,jb,2) * &
#else
                z_dexner_dz_c(2,jc,jk,jb) = -0.5_vp *                              &
                  ((z_theta_v_pr_ic(jc,jk) - z_theta_v_pr_ic(jc,jk+1)) *           &
                  p_nh%metrics%d2dexdz2_fac1_mc(jc,jk,jb) + z_rth_pr(2,jc,jk,jb) * &
#endif
                  p_nh%metrics%d2dexdz2_fac2_mc(jc,jk,jb))
              ENDDO
            ENDDO
            !$ACC END PARALLEL
          ENDIF

        ENDIF ! istep == 1

      ENDDO
!$OMP END DO NOWAIT

      IF (istep == 1) THEN
        ! Add computation of z_grad_rth (perturbation density and virtual potential temperature at main levels)
        ! at outer halo points: needed for correct calculation of the upwind gradients for Miura scheme
        rl_start = min_rlcell_int - 2
        rl_end   = min_rlcell_int - 2

        i_startblk = p_patch%cells%start_block(rl_start)
        i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
#ifdef __SWAPDIM
              z_rth_pr(jc,jk,jb,1) = p_nh%prog(nnow)%rho(jc,jk,jb)     - p_nh%metrics%rho_ref_mc(jc,jk,jb)
              z_rth_pr(jc,jk,jb,2) = p_nh%prog(nnow)%theta_v(jc,jk,jb) - p_nh%metrics%theta_ref_mc(jc,jk,jb)
#else
              z_rth_pr(1,jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb)     - p_nh%metrics%rho_ref_mc(jc,jk,jb)
              z_rth_pr(2,jc,jk,jb) = p_nh%prog(nnow)%theta_v(jc,jk,jb) - p_nh%metrics%theta_ref_mc(jc,jk,jb)
#endif
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO NOWAIT

      ENDIF
!$OMP END PARALLEL

      IF (timers_level > 5) THEN
        CALL timer_stop(timer_solve_nh_cellcomp)
        CALL timer_start(timer_solve_nh_vnupd)
      ENDIF

      ! Compute rho and theta at edges for horizontal flux divergence term
      IF (istep == 1) THEN
        IF (iadv_rhotheta == 1) THEN ! Simplified Miura scheme
          !DA: TODO: remove the wait after everything is async
          !$ACC WAIT
          ! Compute density and potential temperature at vertices
          CALL cells2verts_scalar(p_nh%prog(nnow)%rho,p_patch, p_int%cells_aw_verts, &
            z_rho_v, opt_rlend=min_rlvert_int-1)
          CALL cells2verts_scalar(p_nh%prog(nnow)%theta_v,p_patch, p_int%cells_aw_verts, &
            z_theta_v_v, opt_rlend=min_rlvert_int-1)

        ELSE IF (iadv_rhotheta == 2) THEN ! Miura second-order upwind scheme

          ! Compute Green-Gauss gradients for rho and theta
          CALL grad_green_gauss_cell(z_rth_pr, p_patch, p_int, z_grad_rth,    &
            opt_rlstart=3, opt_rlend=min_rlcell_int-1, opt_acc_async=.TRUE.)
        ENDIF
      ENDIF ! istep = 1

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
      IF (istep == 1) THEN

        ! Compute 'edge values' of density and virtual potential temperature for horizontal
        ! flux divergence term

        ! Initialize halo edges with zero in order to avoid access of uninitialized array elements
        i_startblk = p_patch%edges%start_block(min_rledge_int-2)
        i_endblk   = p_patch%edges%end_block(min_rledge_int-2)

        IF (i_endblk >= i_startblk) THEN
          CALL init_zero_contiguous_dp(z_rho_e    (1,1,i_startblk), nproma*nlev*(i_endblk-i_startblk+1), &
            & opt_acc_async=.TRUE., lacc=i_am_accel_node)
          CALL init_zero_contiguous_dp(z_theta_v_e(1,1,i_startblk), nproma*nlev*(i_endblk-i_startblk+1), &
            & opt_acc_async=.TRUE., lacc=i_am_accel_node)
!$OMP BARRIER
        ENDIF

        rl_start = 7
        rl_end   = min_rledge_int-1

        i_startblk = p_patch%edges%start_block(rl_start)
        i_endblk   = p_patch%edges%end_block  (rl_end)

        ! initialize also nest boundary points with zero
        IF (jg > 1 .OR. l_limited_area) THEN
          CALL init_zero_contiguous_dp(z_rho_e    (1,1,1), nproma*nlev*i_startblk, opt_acc_async=.TRUE., &
          & lacc=i_am_accel_node)
          CALL init_zero_contiguous_dp(z_theta_v_e(1,1,1), nproma*nlev*i_startblk, opt_acc_async=.TRUE., &
          & lacc=i_am_accel_node)
!$OMP BARRIER
        ENDIF

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,ilc0,ibc0,lvn_pos,&
!$OMP            z_ntdistv_bary_1,z_ntdistv_bary_2,distv_bary_1,distv_bary_2) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

          IF (iadv_rhotheta == 2) THEN
            !
            ! Operations from upwind_hflux_miura are inlined in order to process both
            ! fields in one step

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)

            ! Also the back-trajectory computation is inlined to improve efficiency
            !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) &
            !$ACC   PRIVATE(lvn_pos, ilc0, ibc0, z_ntdistv_bary_1, z_ntdistv_bary_2, distv_bary_1, distv_bary_2)
#ifdef __LOOP_EXCHANGE
            DO je = i_startidx, i_endidx
!DIR$ IVDEP, PREFERVECTOR
              DO jk = 1, nlev
#else
            DO jk = 1, nlev
              DO je = i_startidx, i_endidx
#endif
                lvn_pos = p_nh%prog(nnow)%vn(je,jk,jb) >= 0._wp

                ! line and block indices of upwind neighbor cell
                ilc0 = MERGE(p_patch%edges%cell_idx(je,jb,1),p_patch%edges%cell_idx(je,jb,2),lvn_pos)
                ibc0 = MERGE(p_patch%edges%cell_blk(je,jb,1),p_patch%edges%cell_blk(je,jb,2),lvn_pos)

                ! distances from upwind mass point to the end point of the backward trajectory
                ! in edge-normal and tangential directions
                ! (purpose of factor deepatmo_gradh_mc is to modify z_grad_rth below)
                z_ntdistv_bary_1 =  - ( p_nh%prog(nnow)%vn(je,jk,jb) * dthalf +                        &
                  MERGE(p_int%pos_on_tplane_e(je,jb,1,1), p_int%pos_on_tplane_e(je,jb,2,1),lvn_pos)) * &
                  p_nh%metrics%deepatmo_gradh_mc(jk)

                z_ntdistv_bary_2 =  - ( p_nh%diag%vt(je,jk,jb) * dthalf +                              &
                  MERGE(p_int%pos_on_tplane_e(je,jb,1,2), p_int%pos_on_tplane_e(je,jb,2,2),lvn_pos)) * &
                  p_nh%metrics%deepatmo_gradh_mc(jk)

                ! rotate distance vectors into local lat-lon coordinates:
                !
                ! component in longitudinal direction
                distv_bary_1 =                                                                     &
                      z_ntdistv_bary_1*MERGE(p_patch%edges%primal_normal_cell(je,jb,1)%v1,         &
                                             p_patch%edges%primal_normal_cell(je,jb,2)%v1,lvn_pos) &
                    + z_ntdistv_bary_2*MERGE(p_patch%edges%dual_normal_cell(je,jb,1)%v1,           &
                                             p_patch%edges%dual_normal_cell(je,jb,2)%v1,lvn_pos)

                ! component in latitudinal direction
                distv_bary_2 =                                                                     &
                      z_ntdistv_bary_1*MERGE(p_patch%edges%primal_normal_cell(je,jb,1)%v2,         &
                                             p_patch%edges%primal_normal_cell(je,jb,2)%v2,lvn_pos) &
                    + z_ntdistv_bary_2*MERGE(p_patch%edges%dual_normal_cell(je,jb,1)%v2,           &
                                             p_patch%edges%dual_normal_cell(je,jb,2)%v2,lvn_pos)


                ! Calculate "edge values" of rho and theta_v
                ! Note: z_rth_pr contains the perturbation values of rho and theta_v,
                ! and the corresponding gradients are stored in z_grad_rth.
#ifdef __SWAPDIM
                z_rho_e(je,jk,jb) =                                                     &
                  REAL(p_nh%metrics%rho_ref_me(je,jk,jb),wp) + z_rth_pr(ilc0,jk,ibc0,1) &
                  + distv_bary_1 * z_grad_rth(ilc0,jk,ibc0,1) &
                  + distv_bary_2 * z_grad_rth(ilc0,jk,ibc0,2)
                z_theta_v_e(je,jk,jb) =                                                   &
                  REAL(p_nh%metrics%theta_ref_me(je,jk,jb),wp) + z_rth_pr(ilc0,jk,ibc0,2) &
                  + distv_bary_1 * z_grad_rth(ilc0,jk,ibc0,3)                             &
                  + distv_bary_2 * z_grad_rth(ilc0,jk,ibc0,4)
#else
                z_rho_e(je,jk,jb) = REAL(p_nh%metrics%rho_ref_me(je,jk,jb),wp) &
                  +                      z_rth_pr(1,ilc0,jk,ibc0)              &
                  + distv_bary_1 * z_grad_rth(1,ilc0,jk,ibc0)                  &
                  + distv_bary_2 * z_grad_rth(2,ilc0,jk,ibc0)

                z_theta_v_e(je,jk,jb) = REAL(p_nh%metrics%theta_ref_me(je,jk,jb),wp) &
                  +                          z_rth_pr(2,ilc0,jk,ibc0)                &
                  + distv_bary_1 * z_grad_rth(3,ilc0,jk,ibc0)                        &
                  + distv_bary_2 * z_grad_rth(4,ilc0,jk,ibc0)
#endif
              ENDDO   ! loop over vertical levels
            ENDDO ! loop over edges
            !$ACC END PARALLEL

          ELSE

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
            DO je = i_startidx, i_endidx
!DIR$ IVDEP
              DO jk = 1, nlev
#else
            DO jk = 1, nlev
              DO je = i_startidx, i_endidx
#endif

                ! Compute upwind-biased values for rho and theta starting from centered differences
                ! Note: the length of the backward trajectory should be 0.5*dtime*(vn,vt) in order to arrive
                ! at a second-order accurate FV discretization, but twice the length is needed for numerical
                ! stability
                z_rho_e(je,jk,jb) =                                                                          &
                  p_int%c_lin_e(je,1,jb)*p_nh%prog(nnow)%rho(icidx(je,jb,1),jk,icblk(je,jb,1)) +             &
                  p_int%c_lin_e(je,2,jb)*p_nh%prog(nnow)%rho(icidx(je,jb,2),jk,icblk(je,jb,2)) -             &
                  dtime * (p_nh%prog(nnow)%vn(je,jk,jb)*p_patch%edges%inv_dual_edge_length(je,jb)*           &
                 (p_nh%prog(nnow)%rho(icidx(je,jb,2),jk,icblk(je,jb,2)) -                                    &
                  p_nh%prog(nnow)%rho(icidx(je,jb,1),jk,icblk(je,jb,1)) ) + p_nh%diag%vt(je,jk,jb) *         &
                  p_patch%edges%inv_primal_edge_length(je,jb) * p_patch%edges%tangent_orientation(je,jb) *   &
                 (z_rho_v(ividx(je,jb,2),jk,ivblk(je,jb,2)) - z_rho_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) ) )

                z_theta_v_e(je,jk,jb) =                                                                          &
                  p_int%c_lin_e(je,1,jb)*p_nh%prog(nnow)%theta_v(icidx(je,jb,1),jk,icblk(je,jb,1)) +             &
                  p_int%c_lin_e(je,2,jb)*p_nh%prog(nnow)%theta_v(icidx(je,jb,2),jk,icblk(je,jb,2)) -             &
                  dtime * (p_nh%prog(nnow)%vn(je,jk,jb)*p_patch%edges%inv_dual_edge_length(je,jb)*               &
                 (p_nh%prog(nnow)%theta_v(icidx(je,jb,2),jk,icblk(je,jb,2)) -                                    &
                  p_nh%prog(nnow)%theta_v(icidx(je,jb,1),jk,icblk(je,jb,1)) ) + p_nh%diag%vt(je,jk,jb) *         &
                  p_patch%edges%inv_primal_edge_length(je,jb) * p_patch%edges%tangent_orientation(je,jb) *       &
                 (z_theta_v_v(ividx(je,jb,2),jk,ivblk(je,jb,2)) - z_theta_v_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) ))

              ENDDO ! loop over edges
            ENDDO   ! loop over vertical levels
            !$ACC END PARALLEL
          ENDIF

        ENDDO
!$OMP END DO


      ELSE IF (istep == 2 .AND. divdamp_type >= 3) THEN ! apply div damping on 3D divergence

        ! add dw/dz contribution to divergence damping term

        rl_start = 7
        rl_end   = min_rledge_int-2

        i_startblk = p_patch%edges%start_block(rl_start)
        i_endblk   = p_patch%edges%end_block  (rl_end)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
!DIR$ IVDEP, PREFERVECTOR
            DO jk = kstart_dd3d(jg), nlev
              z_graddiv_vn(jk,je,jb) = z_graddiv_vn(jk,je,jb) +  p_nh%metrics%hmask_dd3d(je,jb)*            &
                p_nh%metrics%scalfac_dd3d(jk) * p_patch%edges%inv_dual_edge_length(je,jb)*                  &
                ( z_dwdz_dd(icidx(je,jb,2),jk,icblk(je,jb,2)) - z_dwdz_dd(icidx(je,jb,1),jk,icblk(je,jb,1)) )
#else
          DO jk = kstart_dd3d(jg), nlev
            DO je = i_startidx, i_endidx
              z_graddiv_vn(je,jk,jb) = z_graddiv_vn(je,jk,jb) +  p_nh%metrics%hmask_dd3d(je,jb)*            &
                p_nh%metrics%scalfac_dd3d(jk) * p_patch%edges%inv_dual_edge_length(je,jb)*                  &
                ( z_dwdz_dd(icidx(je,jb,2),jk,icblk(je,jb,2)) - z_dwdz_dd(icidx(je,jb,1),jk,icblk(je,jb,1)) )
#endif
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO

      ENDIF ! istep = 1/2

      ! Remaining computations at edge points

      rl_start = grf_bdywidth_e + 1   ! boundary update follows below
      rl_end   = min_rledge_int

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

      IF (istep == 1) THEN

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_theta1,z_theta2,ikp1,ikp2) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! Store values at nest interface levels; this is done here for the first sub-time step,
          ! the final averaging is done in mo_nh_nest_utilities:compute_tendencies
          IF (idyn_timestep == 1 .AND. l_child_vertnest) THEN
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
            DO je = i_startidx, i_endidx
              p_nh%diag%vn_ie_int(je,1,jb) = p_nh%diag%vn_ie(je,nshift,jb)
            ENDDO
            !$ACC END PARALLEL
          ENDIF

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
            DO jk = 1, nflatlev(jg)-1
#else
          DO jk = 1, nflatlev(jg)-1
            DO je = i_startidx, i_endidx
#endif
              ! horizontal gradient of Exner pressure where coordinate surfaces are flat
              z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)* &
               p_nh%metrics%deepatmo_gradh_mc(jk)*                                 &
               (z_exner_ex_pr(icidx(je,jb,2),jk,icblk(je,jb,2)) -                  &
                z_exner_ex_pr(icidx(je,jb,1),jk,icblk(je,jb,1)) )
            ENDDO
          ENDDO
          !$ACC END PARALLEL

          IF (igradp_method <= 3) THEN

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
            DO je = i_startidx, i_endidx
!DIR$ IVDEP
              DO jk = nflatlev(jg), nflat_gradp(jg)
#else
!$NEC outerloop_unroll(8)
            DO jk = nflatlev(jg), nflat_gradp(jg)
              DO je = i_startidx, i_endidx
#endif
                ! horizontal gradient of Exner pressure, including metric correction
                z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)*         &
                 p_nh%metrics%deepatmo_gradh_mc(jk)*                                         &
                 (z_exner_ex_pr(icidx(je,jb,2),jk,icblk(je,jb,2)) -                          &
                  z_exner_ex_pr(icidx(je,jb,1),jk,icblk(je,jb,1)) ) -                        &
                  p_nh%metrics%ddxn_z_full(je,jk,jb) *                                       &
#ifdef __SWAPDIM
                 (p_int%c_lin_e(je,1,jb)*z_dexner_dz_c(icidx(je,jb,1),jk,icblk(je,jb,1),1) + &
                  p_int%c_lin_e(je,2,jb)*z_dexner_dz_c(icidx(je,jb,2),jk,icblk(je,jb,2),1))
#else
                 (p_int%c_lin_e(je,1,jb)*z_dexner_dz_c(1,icidx(je,jb,1),jk,icblk(je,jb,1)) + &
                  p_int%c_lin_e(je,2,jb)*z_dexner_dz_c(1,icidx(je,jb,2),jk,icblk(je,jb,2)))
#endif
              ENDDO
            ENDDO
            !$ACC END PARALLEL

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
            DO je = i_startidx, i_endidx
!DIR$ IVDEP, PREFERVECTOR
              DO jk = nflat_gradp(jg)+1, nlev
#else
!$NEC outerloop_unroll(8)
            DO jk = nflat_gradp(jg)+1, nlev
              DO je = i_startidx, i_endidx
#endif
                ! horizontal gradient of Exner pressure, Taylor-expansion-based reconstruction
                z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)*          &
                  p_nh%metrics%deepatmo_gradh_mc(jk)*                                         &
                  (z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2)) +           &
                   p_nh%metrics%zdiff_gradp(2,je,jk,jb)*                                      &
#ifdef __SWAPDIM
                  (z_dexner_dz_c(icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2),1) +         &
                   p_nh%metrics%zdiff_gradp(2,je,jk,jb)*                                      &
                   z_dexner_dz_c(icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2),2)) -        &
                  (z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)) +           &
                   p_nh%metrics%zdiff_gradp(1,je,jk,jb)*                                      &
                  (z_dexner_dz_c(icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1),1) +         &
                   p_nh%metrics%zdiff_gradp(1,je,jk,jb)* &
                   z_dexner_dz_c(icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1),2))))
#else
                  (z_dexner_dz_c(1,icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2)) +         &
                   p_nh%metrics%zdiff_gradp(2,je,jk,jb)*                                      &
                   z_dexner_dz_c(2,icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2))) -        &
                  (z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)) +           &
                   p_nh%metrics%zdiff_gradp(1,je,jk,jb)*                                      &
                  (z_dexner_dz_c(1,icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)) +         &
                   p_nh%metrics%zdiff_gradp(1,je,jk,jb)*                                      &
                   z_dexner_dz_c(2,icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)))))
#endif
              ENDDO
            ENDDO
            !$ACC END PARALLEL

          ELSE IF (igradp_method == 4 .OR. igradp_method == 5) THEN

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
            DO je = i_startidx, i_endidx
              DO jk = nflatlev(jg), nlev
#else
            DO jk = nflatlev(jg), nlev
              DO je = i_startidx, i_endidx
#endif
                ! horizontal gradient of Exner pressure, cubic/quadratic interpolation
                z_gradh_exner(je,jk,jb) =  p_patch%edges%inv_dual_edge_length(je,jb)*   &
                  p_nh%metrics%deepatmo_gradh_mc(jk)*                                   &
                  (z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)-1,icblk(je,jb,2)) *   &
                   p_nh%metrics%coeff_gradp(5,je,jk,jb) +                               &
                   z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)  ,icblk(je,jb,2)) *   &
                   p_nh%metrics%coeff_gradp(6,je,jk,jb) +                               &
                   z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)+1,icblk(je,jb,2)) *   &
                   p_nh%metrics%coeff_gradp(7,je,jk,jb) +                               &
                   z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)+2,icblk(je,jb,2)) *   &
                   p_nh%metrics%coeff_gradp(8,je,jk,jb) -                               &
                  (z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)-1,icblk(je,jb,1)) *   &
                   p_nh%metrics%coeff_gradp(1,je,jk,jb) +                               &
                   z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)  ,icblk(je,jb,1)) *   &
                   p_nh%metrics%coeff_gradp(2,je,jk,jb) +                               &
                   z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)+1,icblk(je,jb,1)) *   &
                   p_nh%metrics%coeff_gradp(3,je,jk,jb) +                               &
                   z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)+2,icblk(je,jb,1)) *   &
                  p_nh%metrics%coeff_gradp(4,je,jk,jb)) )

              ENDDO
            ENDDO
            !$ACC END PARALLEL
          ENDIF

          ! compute hydrostatically approximated correction term that replaces downward extrapolation
          IF (igradp_method == 3) THEN

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR PRIVATE(z_theta1, z_theta2)
            DO je = i_startidx, i_endidx

              z_theta1 = &
                 p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb),icblk(je,jb,1)) +  &
                 p_nh%metrics%zdiff_gradp(1,je,nlev,jb)*                                       &
                (p_nh%diag%theta_v_ic(icidx(je,jb,1),ikidx(1,je,nlev,jb),  icblk(je,jb,1)) -   &
                 p_nh%diag%theta_v_ic(icidx(je,jb,1),ikidx(1,je,nlev,jb)+1,icblk(je,jb,1))) *  &
                 p_nh%metrics%inv_ddqz_z_full(icidx(je,jb,1),ikidx(1,je,nlev,jb),icblk(je,jb,1))

              z_theta2 = &
                 p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb),icblk(je,jb,2)) +  &
                 p_nh%metrics%zdiff_gradp(2,je,nlev,jb)*                                       &
                (p_nh%diag%theta_v_ic(icidx(je,jb,2),ikidx(2,je,nlev,jb),  icblk(je,jb,2)) -   &
                 p_nh%diag%theta_v_ic(icidx(je,jb,2),ikidx(2,je,nlev,jb)+1,icblk(je,jb,2))) *  &
                 p_nh%metrics%inv_ddqz_z_full(icidx(je,jb,2),ikidx(2,je,nlev,jb),icblk(je,jb,2))

              z_hydro_corr(je,jb) = grav_o_cpd*p_patch%edges%inv_dual_edge_length(je,jb)*    &
                (z_theta2-z_theta1)*4._wp/(z_theta1+z_theta2)**2

            ENDDO
            !$ACC END PARALLEL

          ELSE IF (igradp_method == 5) THEN

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR PRIVATE(ikp1, ikp2, z_theta1, z_theta2)
            DO je = i_startidx, i_endidx

              ikp1 = MIN(nlev,ikidx(1,je,nlev,jb)+2)
              ikp2 = MIN(nlev,ikidx(2,je,nlev,jb)+2)

              z_theta1 =                                                                       &
                p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb)-1,icblk(je,jb,1)) * &
                p_nh%metrics%coeff_gradp(1,je,nlev,jb) +                                         &
                p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb)  ,icblk(je,jb,1)) * &
                p_nh%metrics%coeff_gradp(2,je,nlev,jb) +                                         &
                p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb)+1,icblk(je,jb,1)) * &
                p_nh%metrics%coeff_gradp(3,je,nlev,jb) +                                         &
                p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikp1                 ,icblk(je,jb,1)) * &
                p_nh%metrics%coeff_gradp(4,je,nlev,jb)

              z_theta2 =                                                                       &
                p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb)-1,icblk(je,jb,2)) * &
                p_nh%metrics%coeff_gradp(5,je,nlev,jb) +                                         &
                p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb)  ,icblk(je,jb,2)) * &
                p_nh%metrics%coeff_gradp(6,je,nlev,jb) +                                         &
                p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb)+1,icblk(je,jb,2)) * &
                p_nh%metrics%coeff_gradp(7,je,nlev,jb) +                                         &
                p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikp2                 ,icblk(je,jb,2)) * &
                p_nh%metrics%coeff_gradp(8,je,nlev,jb)

              z_hydro_corr(je,jb) = grav_o_cpd*p_patch%edges%inv_dual_edge_length(je,jb)*    &
                (z_theta2-z_theta1)*4._wp/(z_theta1+z_theta2)**2

            ENDDO
            !$ACC END PARALLEL
          ENDIF

        ENDDO
!$OMP END DO

      ENDIF ! istep = 1


      IF (istep == 1 .AND. (igradp_method == 3 .OR. igradp_method == 5)) THEN

!$OMP DO PRIVATE(jb,je,ie,nlen_gradp,ishift) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, nblks_gradp
          IF (jb == nblks_gradp) THEN
            nlen_gradp = npromz_gradp
          ELSE
            nlen_gradp = nproma_gradp
          ENDIF
          ishift = (jb-1)*nproma_gradp
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
!$NEC ivdep
          !$ACC LOOP GANG VECTOR
          DO je = 1, nlen_gradp
            ie = ishift+je

            z_gradh_exner(ipeidx(ie),iplev(ie),ipeblk(ie))  =              &
              z_gradh_exner(ipeidx(ie),iplev(ie),ipeblk(ie)) +             &
              p_nh%metrics%pg_exdist(ie)*z_hydro_corr(ipeidx(ie),ipeblk(ie))

          ENDDO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO

      ENDIF

      ! Update horizontal velocity field: advection, Coriolis force, pressure-gradient term, and physics

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_graddiv2_vn,                                                   &
!$OMP            z_ddt_vn_dyn, z_ddt_vn_apc, z_ddt_vn_cor, z_ddt_vn_pgr, z_ddt_vn_ray, z_d_vn_dmp, z_d_vn_iau  &
!$OMP           ) ICON_OMP_DEFAULT_SCHEDULE

      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        IF ((itime_scheme >= 4) .AND. istep == 2) THEN ! use temporally averaged velocity advection terms

          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(z_ddt_vn_dyn, z_ddt_vn_apc, z_ddt_vn_cor, z_ddt_vn_pgr) TILE(32, 4)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              !
              z_ddt_vn_apc                      =  p_nh%diag%ddt_vn_apc_pc(je,jk,jb,ntl1)*wgt_nnow_vel  &
                &                                 +p_nh%diag%ddt_vn_apc_pc(je,jk,jb,ntl2)*wgt_nnew_vel
              z_ddt_vn_pgr                      = -cpd*z_theta_v_e(je,jk,jb)*z_gradh_exner(je,jk,jb)
              !
              z_ddt_vn_dyn                      =  z_ddt_vn_apc                   & ! advection plus Coriolis
                &                                 +z_ddt_vn_pgr                   & ! pressure gradient
                &                                 +p_nh%diag%ddt_vn_phy(je,jk,jb)   ! physics applied in dynamics
              !
              p_nh%prog(nnew)%vn(je,jk,jb)      =  p_nh%prog(nnow)%vn(je,jk,jb)   + dtime       * z_ddt_vn_dyn
              !
#ifdef __ENABLE_DDT_VN_XYZ__
              IF (p_nh%diag%ddt_vn_adv_is_associated .OR. p_nh%diag%ddt_vn_cor_is_associated) THEN
                z_ddt_vn_cor                    =  p_nh%diag%ddt_vn_cor_pc(je,jk,jb,ntl1)*wgt_nnow_vel  &
                  &                               +p_nh%diag%ddt_vn_cor_pc(je,jk,jb,ntl2)*wgt_nnew_vel
                !
                IF (p_nh%diag%ddt_vn_adv_is_associated) THEN
                  p_nh%diag%ddt_vn_adv(je,jk,jb)=  p_nh%diag%ddt_vn_adv(je,jk,jb) + r_nsubsteps *(z_ddt_vn_apc-z_ddt_vn_cor)
                END IF
                !
                IF (p_nh%diag%ddt_vn_cor_is_associated) THEN
                  p_nh%diag%ddt_vn_cor(je,jk,jb)=  p_nh%diag%ddt_vn_cor(je,jk,jb) + r_nsubsteps * z_ddt_vn_cor
                END IF
                !
              END IF
              !
              IF (p_nh%diag%ddt_vn_pgr_is_associated) THEN
                p_nh%diag%ddt_vn_pgr(je,jk,jb)  =  p_nh%diag%ddt_vn_pgr(je,jk,jb) + r_nsubsteps * z_ddt_vn_pgr
              END IF
              !
              IF (p_nh%diag%ddt_vn_phd_is_associated) THEN
                p_nh%diag%ddt_vn_phd(je,jk,jb)  =  p_nh%diag%ddt_vn_phd(je,jk,jb) + r_nsubsteps * p_nh%diag%ddt_vn_phy(je,jk,jb)
              END IF
              !
              IF (p_nh%diag%ddt_vn_dyn_is_associated) THEN
                p_nh%diag%ddt_vn_dyn(je,jk,jb)  =  p_nh%diag%ddt_vn_dyn(je,jk,jb) + r_nsubsteps * z_ddt_vn_dyn
              END IF
#endif
              !
            ENDDO
          ENDDO

        ELSE

          !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              !
              p_nh%prog(nnew)%vn(je,jk,jb)      =  p_nh%prog(nnow)%vn(je,jk,jb)   + dtime *                 &
                &                                ( p_nh%diag%ddt_vn_apc_pc(je,jk,jb,ntl1)                   &
                &                                 -cpd*z_theta_v_e(je,jk,jb)*z_gradh_exner(je,jk,jb)        &
                &                                 +p_nh%diag%ddt_vn_phy(je,jk,jb)                        )
              !
            ENDDO
          ENDDO
        ENDIF


        IF (istep == 2) THEN

          IF (divdamp_order == 4 .OR. divdamp_order == 24) THEN ! fourth-order divergence damping
          ! Compute gradient of divergence of gradient of divergence for fourth-order divergence damping
            !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
            DO je = i_startidx, i_endidx
!DIR$ IVDEP
              DO jk = 1, nlev
                z_graddiv2_vn(je,jk) = p_int%geofac_grdiv(je,1,jb)*z_graddiv_vn(jk,je,jb)      &
                  + p_int%geofac_grdiv(je,2,jb)*z_graddiv_vn(jk,iqidx(je,jb,1),iqblk(je,jb,1)) &
                  + p_int%geofac_grdiv(je,3,jb)*z_graddiv_vn(jk,iqidx(je,jb,2),iqblk(je,jb,2)) &
                  + p_int%geofac_grdiv(je,4,jb)*z_graddiv_vn(jk,iqidx(je,jb,3),iqblk(je,jb,3)) &
                  + p_int%geofac_grdiv(je,5,jb)*z_graddiv_vn(jk,iqidx(je,jb,4),iqblk(je,jb,4))
#else
!$NEC outerloop_unroll(6)
            DO jk = 1, nlev
              DO je = i_startidx, i_endidx
                z_graddiv2_vn(je,jk) = p_int%geofac_grdiv(je,1,jb)*z_graddiv_vn(je,jk,jb)      &
                  + p_int%geofac_grdiv(je,2,jb)*z_graddiv_vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
                  + p_int%geofac_grdiv(je,3,jb)*z_graddiv_vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
                  + p_int%geofac_grdiv(je,4,jb)*z_graddiv_vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
                  + p_int%geofac_grdiv(je,5,jb)*z_graddiv_vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))
#endif

              ENDDO
            ENDDO
          ENDIF

          ! apply divergence damping if diffusion is not called every sound-wave time step
          IF (divdamp_order == 2 .OR. (divdamp_order == 24 .AND. scal_divdamp_o2 > 1.e-6_wp) ) THEN ! 2nd-order divergence damping

            !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) PRIVATE(z_d_vn_dmp)
            DO jk = 1, nlev
!DIR$ IVDEP
              DO je = i_startidx, i_endidx
                !
#ifdef __LOOP_EXCHANGE
                z_d_vn_dmp = scal_divdamp_o2*z_graddiv_vn(jk,je,jb)
#else
                z_d_vn_dmp = scal_divdamp_o2*z_graddiv_vn(je,jk,jb)
#endif
                !
                p_nh%prog(nnew)%vn(je,jk,jb)      =  p_nh%prog(nnew)%vn(je,jk,jb)   + z_d_vn_dmp
                !
#ifdef __ENABLE_DDT_VN_XYZ__
                IF (p_nh%diag%ddt_vn_dmp_is_associated) THEN
                  p_nh%diag%ddt_vn_dmp(je,jk,jb)  =  p_nh%diag%ddt_vn_dmp(je,jk,jb) + z_d_vn_dmp * r_dtimensubsteps
                END IF
                !
                IF (p_nh%diag%ddt_vn_dyn_is_associated) THEN
                  p_nh%diag%ddt_vn_dyn(je,jk,jb)  =  p_nh%diag%ddt_vn_dyn(je,jk,jb) + z_d_vn_dmp * r_dtimensubsteps
                END IF
#endif
                !
              ENDDO
            ENDDO
          ENDIF
          IF (divdamp_order == 4 .OR. (divdamp_order == 24 .AND. divdamp_fac_o2 <= 4._wp*divdamp_fac) ) THEN
            IF (l_limited_area .OR. jg > 1) THEN
              ! fourth-order divergence damping with reduced damping coefficient along nest boundary
              ! (scal_divdamp is negative whereas bdy_divdamp is positive; decreasing the divergence
              ! damping along nest boundaries is beneficial because this reduces the interference
              ! with the increased diffusion applied in nh_diffusion)

              !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) PRIVATE(z_d_vn_dmp)
              DO jk = 1, nlev
!DIR$ IVDEP
!$NEC ivdep
                DO je = i_startidx, i_endidx
                  !
                  z_d_vn_dmp = (scal_divdamp(jk)+bdy_divdamp(jk)*p_int%nudgecoeff_e(je,jb))*z_graddiv2_vn(je,jk)
                  !
                  p_nh%prog(nnew)%vn(je,jk,jb)      =  p_nh%prog(nnew)%vn(je,jk,jb)   + z_d_vn_dmp
                  !
#ifdef __ENABLE_DDT_VN_XYZ__
                  IF (p_nh%diag%ddt_vn_dmp_is_associated) THEN
                    p_nh%diag%ddt_vn_dmp(je,jk,jb)  =  p_nh%diag%ddt_vn_dmp(je,jk,jb) + z_d_vn_dmp * r_dtimensubsteps
                  END IF
                  !
                  IF (p_nh%diag%ddt_vn_dyn_is_associated) THEN
                    p_nh%diag%ddt_vn_dyn(je,jk,jb)  =  p_nh%diag%ddt_vn_dyn(je,jk,jb) + z_d_vn_dmp * r_dtimensubsteps
                  END IF
#endif
                  !
                ENDDO
              ENDDO
            ELSE ! fourth-order divergence damping

              !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) PRIVATE(z_d_vn_dmp)
              DO jk = 1, nlev
!DIR$ IVDEP
                DO je = i_startidx, i_endidx
                  !
                  z_d_vn_dmp = scal_divdamp(jk)*z_graddiv2_vn(je,jk)
                  !
                  p_nh%prog(nnew)%vn(je,jk,jb)      =  p_nh%prog(nnew)%vn(je,jk,jb)   + z_d_vn_dmp
                  !
#ifdef __ENABLE_DDT_VN_XYZ__
                  IF (p_nh%diag%ddt_vn_dmp_is_associated) THEN
                    p_nh%diag%ddt_vn_dmp(je,jk,jb)  =  p_nh%diag%ddt_vn_dmp(je,jk,jb) + z_d_vn_dmp * r_dtimensubsteps
                  END IF
                  !
                  IF (p_nh%diag%ddt_vn_dyn_is_associated) THEN
                    p_nh%diag%ddt_vn_dyn(je,jk,jb)  =  p_nh%diag%ddt_vn_dyn(je,jk,jb) + z_d_vn_dmp * r_dtimensubsteps
                  END IF
#endif
                  !
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF

        IF (is_iau_active) THEN ! add analysis increment from data assimilation

          !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) PRIVATE(z_d_vn_iau)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              !
              z_d_vn_iau = iau_wgt_dyn*p_nh%diag%vn_incr(je,jk,jb)
              !
              p_nh%prog(nnew)%vn(je,jk,jb)        =  p_nh%prog(nnew)%vn(je,jk,jb)   + z_d_vn_iau
              !
#ifdef __ENABLE_DDT_VN_XYZ__
              IF (istep == 2) THEN
                IF (p_nh%diag%ddt_vn_iau_is_associated) THEN
                  p_nh%diag%ddt_vn_iau(je,jk,jb)  =  p_nh%diag%ddt_vn_iau(je,jk,jb) + z_d_vn_iau * r_dtimensubsteps
                END IF
                !
                IF (p_nh%diag%ddt_vn_dyn_is_associated) THEN
                  p_nh%diag%ddt_vn_dyn(je,jk,jb)  =  p_nh%diag%ddt_vn_dyn(je,jk,jb) + z_d_vn_iau * r_dtimensubsteps
                END IF
              END IF
#endif
              !
            ENDDO
          ENDDO
        ENDIF
        !$ACC END PARALLEL

        ! Classic Rayleigh damping mechanism for vn (requires reference state !!)
        !
        IF ( rayleigh_type == RAYLEIGH_CLASSIC ) THEN

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_ddt_vn_ray)
          DO jk = 1, nrdmax(jg)
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              !
              z_ddt_vn_ray = -p_nh%metrics%rayleigh_vn(jk) * (p_nh%prog(nnew)%vn(je,jk,jb) - p_nh%ref%vn_ref(je,jk,jb))
              !
              p_nh%prog(nnew)%vn(je,jk,jb)        =  p_nh%prog(nnew)%vn(je,jk,jb)   + z_ddt_vn_ray * dtime
              !
#ifdef __ENABLE_DDT_VN_XYZ__
              IF (istep == 2) THEN
                IF (p_nh%diag%ddt_vn_ray_is_associated) THEN
                  p_nh%diag%ddt_vn_ray(je,jk,jb)  =  p_nh%diag%ddt_vn_ray(je,jk,jb) + z_ddt_vn_ray * r_nsubsteps
                END IF
                !
                IF (p_nh%diag%ddt_vn_dyn_is_associated) THEN
                  p_nh%diag%ddt_vn_dyn(je,jk,jb)  =  p_nh%diag%ddt_vn_dyn(je,jk,jb) + z_ddt_vn_ray * r_nsubsteps
                END IF
              END IF
#endif
              !
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF
      ENDDO
!$OMP END DO

      ! Boundary update of horizontal velocity
      IF (istep == 1 .AND. (l_limited_area .OR. jg > 1)) THEN
        rl_start = 1
        rl_end   = grf_bdywidth_e

        i_startblk = p_patch%edges%start_block(rl_start)
        i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
            i_startidx, i_endidx, rl_start, rl_end)

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              !
              p_nh%prog(nnew)%vn(je,jk,jb)      =  p_nh%prog(nnow)%vn(je,jk,jb)   + p_nh%diag%grf_tend_vn(je,jk,jb) * dtime
              !
#ifdef __ENABLE_DDT_VN_XYZ__
              IF (p_nh%diag%ddt_vn_grf_is_associated) THEN
                p_nh%diag%ddt_vn_grf(je,jk,jb)  =  p_nh%diag%ddt_vn_grf(je,jk,jb) + p_nh%diag%grf_tend_vn(je,jk,jb) * r_nsubsteps
              END IF
              !
              IF (p_nh%diag%ddt_vn_dyn_is_associated) THEN
                p_nh%diag%ddt_vn_dyn(je,jk,jb)  =  p_nh%diag%ddt_vn_dyn(je,jk,jb) + p_nh%diag%grf_tend_vn(je,jk,jb) * r_nsubsteps
              END IF
#endif
              !
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO

      ENDIF

      ! Preparations for nest boundary interpolation of mass fluxes from parent domain
      IF (jg > 1 .AND. grf_intmethod_e == 6 .AND. jstep == 0 .AND. istep == 1) THEN

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG

!$OMP DO PRIVATE(ic,je,jb,jk) ICON_OMP_DEFAULT_SCHEDULE
        DO ic = 1, p_nh%metrics%bdy_mflx_e_dim
          je = p_nh%metrics%bdy_mflx_e_idx(ic)
          jb = p_nh%metrics%bdy_mflx_e_blk(ic)
!DIR$ IVDEP
          !$ACC LOOP VECTOR
          DO jk = 1, nlev
            p_nh%diag%grf_bdy_mflx(jk,ic,2) = p_nh%diag%grf_tend_mflx(je,jk,jb)
            p_nh%diag%grf_bdy_mflx(jk,ic,1) = prep_adv%mass_flx_me(je,jk,jb) - dt_shift*p_nh%diag%grf_bdy_mflx(jk,ic,2)
          ENDDO

        ENDDO
!$OMP END DO

        !$ACC END PARALLEL

      ENDIF

!$OMP END PARALLEL


      !-------------------------
      ! communication phase
      IF (timers_level > 5) THEN
        CALL timer_stop(timer_solve_nh_vnupd)
        CALL timer_start(timer_solve_nh_exch)
      ENDIF

      IF (istep == 1) THEN
        CALL sync_patch_array_mult(SYNC_E,p_patch,2,p_nh%prog(nnew)%vn,z_rho_e,opt_varname="vn_nnew and z_rho_e")
      ELSE
        CALL sync_patch_array(SYNC_E,p_patch,p_nh%prog(nnew)%vn,opt_varname="vn_nnew")
      ENDIF


      IF (timers_level > 5) THEN
        CALL timer_stop(timer_solve_nh_exch)
        CALL timer_start(timer_solve_nh_edgecomp)
      ENDIF
      ! end communication phase
      !-------------------------

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
      rl_start = 5
      rl_end   = min_rledge_int - 2

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_vn_avg) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        IF (istep == 1) THEN

          !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
!DIR$ IVDEP
            DO jk = 1, nlev
#else
!$NEC outerloop_unroll(8)
          DO jk = 1, nlev
!$NEC vovertake
            DO je = i_startidx, i_endidx
#endif
              ! Average normal wind components in order to get nearly second-order accurate divergence
              z_vn_avg(je,jk) = p_int%e_flx_avg(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)           &
                + p_int%e_flx_avg(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
                + p_int%e_flx_avg(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
                + p_int%e_flx_avg(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
                + p_int%e_flx_avg(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

              ! Compute gradient of divergence of vn for divergence damping
#ifdef __LOOP_EXCHANGE
              z_graddiv_vn(jk,je,jb) = p_int%geofac_grdiv(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)    &
#else
              z_graddiv_vn(je,jk,jb) = p_int%geofac_grdiv(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)    &
#endif
              + p_int%geofac_grdiv(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
                + p_int%geofac_grdiv(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
                + p_int%geofac_grdiv(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
                + p_int%geofac_grdiv(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

              ! RBF reconstruction of tangential wind component
              p_nh%diag%vt(je,jk,jb) = p_int%rbf_vec_coeff_e(1,je,jb)  &
                * p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
                + p_int%rbf_vec_coeff_e(2,je,jb)                       &
                * p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
                + p_int%rbf_vec_coeff_e(3,je,jb)                       &
                * p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
                + p_int%rbf_vec_coeff_e(4,je,jb)                       &
                * p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))
            ENDDO
          ENDDO

        ELSE IF (itime_scheme >= 5) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
!DIR$ IVDEP
            DO jk = 1, nlev
#else
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
#endif
              ! Average normal wind components in order to get nearly second-order accurate divergence
              z_vn_avg(je,jk) = p_int%e_flx_avg(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)           &
                + p_int%e_flx_avg(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
                + p_int%e_flx_avg(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
                + p_int%e_flx_avg(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
                + p_int%e_flx_avg(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

              ! RBF reconstruction of tangential wind component
              p_nh%diag%vt(je,jk,jb) = p_int%rbf_vec_coeff_e(1,je,jb)  &
                * p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
                + p_int%rbf_vec_coeff_e(2,je,jb)                       &
                * p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
                + p_int%rbf_vec_coeff_e(3,je,jb)                       &
                * p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
                + p_int%rbf_vec_coeff_e(4,je,jb)                       &
                * p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

            ENDDO
          ENDDO

        ELSE

          !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
!DIR$ IVDEP
            DO jk = 1, nlev
#else
!$NEC outerloop_unroll(8)
          DO jk = 1, nlev
!$NEC vovertake
            DO je = i_startidx, i_endidx
#endif
              ! Average normal wind components in order to get nearly second-order accurate divergence
              z_vn_avg(je,jk) = p_int%e_flx_avg(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)           &
                + p_int%e_flx_avg(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
                + p_int%e_flx_avg(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
                + p_int%e_flx_avg(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
                + p_int%e_flx_avg(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))
            ENDDO
          ENDDO
        ENDIF


        ! Compute fluxes at edges using averaged velocities
        !
        !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
        DO jk = 1,nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx

            p_nh%diag%mass_fl_e(je,jk,jb) = z_rho_e(je,jk,jb) *        &
              z_vn_avg(je,jk) * p_nh%metrics%ddqz_z_full_e(je,jk,jb)
            z_theta_v_fl_e(je,jk,jb) = p_nh%diag%mass_fl_e(je,jk,jb) * &
              z_theta_v_e(je,jk,jb)

          ENDDO
        ENDDO

        IF (lapply_moist_thdyn) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              z_vol_fl_e(je,jk,jb) = z_vn_avg(je,jk) * p_nh%metrics%ddqz_z_full_e(je,jk,jb)
            ENDDO
          ENDDO
        ENDIF

        IF (lsave_mflx .AND. istep == 2) THEN ! store mass flux for nest boundary interpolation
#ifndef _OPENACC
          DO je = i_startidx, i_endidx
            IF (p_patch%edges%refin_ctrl(je,jb) <= -4 .AND. p_patch%edges%refin_ctrl(je,jb) >= -6) THEN
!DIR$ IVDEP
              DO jk=1,nlev
                p_nh%diag%mass_fl_e_sv(je,jk,jb) = p_nh%diag%mass_fl_e(je,jk,jb)
              ENDDO
            ENDIF
          ENDDO
#else
            !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
            DO jk=1,nlev
              DO je = i_startidx, i_endidx
                IF (p_patch%edges%refin_ctrl(je,jb) <= -4 .AND. p_patch%edges%refin_ctrl(je,jb) >= -6) THEN
                  p_nh%diag%mass_fl_e_sv(je,jk,jb) = p_nh%diag%mass_fl_e(je,jk,jb)
                ENDIF
              ENDDO
            ENDDO
#endif
        ENDIF

        IF (lprep_adv .AND. istep == 2) THEN ! Preprations for tracer advection
          IF (lclean_mflx) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
            DO jk = 1, nlev
!$NEC ivdep
              DO je = i_startidx, i_endidx
                prep_adv%vn_traj(je,jk,jb)     = 0._wp
                prep_adv%mass_flx_me(je,jk,jb) = 0._wp
              ENDDO
            ENDDO
          ENDIF

          !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4)
          DO jk = 1, nlev
!$NEC ivdep
            DO je = i_startidx, i_endidx
              prep_adv%vn_traj(je,jk,jb)     = prep_adv%vn_traj(je,jk,jb)     + r_nsubsteps*z_vn_avg(je,jk)
              prep_adv%mass_flx_me(je,jk,jb) = prep_adv%mass_flx_me(je,jk,jb) + r_nsubsteps*p_nh%diag%mass_fl_e(je,jk,jb)
            ENDDO
          ENDDO
        ENDIF
        !$ACC END PARALLEL

        IF (istep == 1 .OR. itime_scheme >= 5) THEN
          ! Compute contravariant correction for vertical velocity at full levels

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = nflatlev(jg), nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              z_w_concorr_me(je,jk,jb) =                                          &
                p_nh%prog(nnew)%vn(je,jk,jb)*p_nh%metrics%ddxn_z_full(je,jk,jb) + &
                p_nh%diag%vt(je,jk,jb)      *p_nh%metrics%ddxt_z_full(je,jk,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        IF (istep == 1) THEN
          ! Interpolate vn to interface levels and compute horizontal part of kinetic energy on edges
          ! (needed in velocity tendencies called at istep=2)

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(3)
          DO jk = 2, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              p_nh%diag%vn_ie(je,jk,jb) =                                                    &
                           p_nh%metrics%wgtfac_e(je,jk,jb) *p_nh%prog(nnew)%vn(je,jk  ,jb) + &
                  (1._wp - p_nh%metrics%wgtfac_e(je,jk,jb))*p_nh%prog(nnew)%vn(je,jk-1,jb)
              z_vt_ie(je,jk,jb) =                                                      &
                           p_nh%metrics%wgtfac_e(je,jk,jb) *p_nh%diag%vt(je,jk  ,jb) + &
                  (1._wp - p_nh%metrics%wgtfac_e(je,jk,jb))*p_nh%diag%vt(je,jk-1,jb)
              z_kin_hor_e(je,jk,jb) = 0.5_wp*(p_nh%prog(nnew)%vn(je,jk,jb)**2 + p_nh%diag%vt(je,jk,jb)**2)
            ENDDO
          ENDDO
          !$ACC END PARALLEL

          IF (.NOT. l_vert_nested) THEN
            ! Top and bottom levels
!DIR$ IVDEP
            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR
            DO je = i_startidx, i_endidx
              ! Quadratic extrapolation at the top turned out to cause numerical instability in pathological cases,
              ! thus we use a no-gradient condition in the upper half layer
              p_nh%diag%vn_ie(je,1,jb) = p_nh%prog(nnew)%vn(je,1,jb)
              ! vt_ie(jk=1) is actually unused, but we need it for convenience of implementation
              z_vt_ie(je,1,jb) = p_nh%diag%vt(je,1,jb)
              !
              z_kin_hor_e(je,1,jb) = 0.5_wp*(p_nh%prog(nnew)%vn(je,1,jb)**2 + p_nh%diag%vt(je,1,jb)**2)
              p_nh%diag%vn_ie(je,nlevp1,jb) =                           &
                p_nh%metrics%wgtfacq_e(je,1,jb)*p_nh%prog(nnew)%vn(je,nlev,jb) +   &
                p_nh%metrics%wgtfacq_e(je,2,jb)*p_nh%prog(nnew)%vn(je,nlev-1,jb) + &
                p_nh%metrics%wgtfacq_e(je,3,jb)*p_nh%prog(nnew)%vn(je,nlev-2,jb)
            ENDDO
            !$ACC END PARALLEL
          ELSE
            ! vn_ie(jk=1) is interpolated horizontally from the parent domain, and linearly interpolated in time
            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              p_nh%diag%vn_ie(je,1,jb) = p_nh%diag%vn_ie_ubc(je,1,jb)+dt_linintp_ubc_nnew*p_nh%diag%vn_ie_ubc(je,2,jb)
              ! vt_ie(jk=1) is actually unused, but we need it for convenience of implementation
              z_vt_ie(je,1,jb) = p_nh%diag%vt(je,1,jb)
              !
              z_kin_hor_e(je,1,jb) = 0.5_wp*(p_nh%prog(nnew)%vn(je,1,jb)**2 + p_nh%diag%vt(je,1,jb)**2)
              p_nh%diag%vn_ie(je,nlevp1,jb) =                           &
                p_nh%metrics%wgtfacq_e(je,1,jb)*p_nh%prog(nnew)%vn(je,nlev,jb) +   &
                p_nh%metrics%wgtfacq_e(je,2,jb)*p_nh%prog(nnew)%vn(je,nlev-1,jb) + &
                p_nh%metrics%wgtfacq_e(je,3,jb)*p_nh%prog(nnew)%vn(je,nlev-2,jb)
            ENDDO
            !$ACC END PARALLEL
          ENDIF
        ENDIF

      ENDDO
!$OMP END DO

      ! Apply mass fluxes across lateral nest boundary interpolated from parent domain
      IF (jg > 1 .AND. grf_intmethod_e == 6) THEN

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        ! PGI 21.2 requires GANG-VECTOR on this level. (Having the jk as VECTOR crashes.)
        ! PRIVATE clause is required as je,jb are used in each vector thread.
        !$ACC LOOP GANG VECTOR PRIVATE(je, jb)

!$OMP DO PRIVATE(ic,je,jb,jk) ICON_OMP_DEFAULT_SCHEDULE
        DO ic = 1, p_nh%metrics%bdy_mflx_e_dim
          je = p_nh%metrics%bdy_mflx_e_idx(ic)
          jb = p_nh%metrics%bdy_mflx_e_blk(ic)

          ! This is needed for tracer mass consistency along the lateral boundaries
          IF (lprep_adv .AND. istep == 2) THEN ! subtract mass flux added previously...
            !$ACC LOOP SEQ
!$NEC ivdep
            DO jk = 1, nlev
              prep_adv%mass_flx_me(je,jk,jb) = prep_adv%mass_flx_me(je,jk,jb) - r_nsubsteps*p_nh%diag%mass_fl_e(je,jk,jb)
              prep_adv%vn_traj(je,jk,jb)     = prep_adv%vn_traj(je,jk,jb) - r_nsubsteps*p_nh%diag%mass_fl_e(je,jk,jb) / &
                (z_rho_e(je,jk,jb) * p_nh%metrics%ddqz_z_full_e(je,jk,jb))
            ENDDO
          ENDIF

!DIR$ IVDEP
          !$ACC LOOP SEQ
!$NEC ivdep
          DO jk = 1, nlev
            p_nh%diag%mass_fl_e(je,jk,jb) = p_nh%diag%grf_bdy_mflx(jk,ic,1) + &
              REAL(jstep,wp)*dtime*p_nh%diag%grf_bdy_mflx(jk,ic,2)
            z_theta_v_fl_e(je,jk,jb) = p_nh%diag%mass_fl_e(je,jk,jb) * z_theta_v_e(je,jk,jb)
          ENDDO

          IF (lapply_moist_thdyn) THEN
            !$ACC LOOP SEQ
            DO jk = 1, nlev
              z_vol_fl_e(je,jk,jb) = p_nh%diag%mass_fl_e(je,jk,jb) / z_rho_e(je,jk,jb)
            ENDDO
          ENDIF

          IF (lprep_adv .AND. istep == 2) THEN ! ... and add the corrected one again
            !$ACC LOOP SEQ
!$NEC ivdep
            DO jk = 1, nlev
              prep_adv%mass_flx_me(je,jk,jb) = prep_adv%mass_flx_me(je,jk,jb) + r_nsubsteps*p_nh%diag%mass_fl_e(je,jk,jb)
              prep_adv%vn_traj(je,jk,jb)     = prep_adv%vn_traj(je,jk,jb) + r_nsubsteps*p_nh%diag%mass_fl_e(je,jk,jb) / &
                (z_rho_e(je,jk,jb) * p_nh%metrics%ddqz_z_full_e(je,jk,jb))
            ENDDO
          ENDIF

        ENDDO
!$OMP END DO

      !$ACC END PARALLEL

      ENDIF


      ! It turned out that it is sufficient to compute the contravariant correction in the
      ! predictor step at time level n+1; repeating the calculation in the corrector step
      ! has negligible impact on the results except in very-high resolution runs with extremely steep mountains
      IF (istep == 1 .OR. itime_scheme >= 5) THEN

        rl_start = 3
        rl_end = min_rlcell_int - 1

        i_startblk = p_patch%cells%start_block(rl_start)
        i_endblk   = p_patch%cells%end_block(rl_end)

#ifdef _OPENACC
!
! This is one of the very few code divergences for OPENACC (see comment below)
!
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
            i_startidx, i_endidx, rl_start, rl_end)

          ! ... and to interface levels
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR TILE(32, 4) PRIVATE(z_w_concorr_mc_m1, z_w_concorr_mc_m0)
          DO jk = nflatlev(jg)+1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              ! COMMENT: this optimization yields drastically better performance in an OpenACC context
              ! Interpolate contravariant correction to cell centers...
              z_w_concorr_mc_m1 =  &
                p_int%e_bln_c_s(jc,1,jb)*z_w_concorr_me(ieidx(jc,jb,1),jk-1,ieblk(jc,jb,1)) + &
                p_int%e_bln_c_s(jc,2,jb)*z_w_concorr_me(ieidx(jc,jb,2),jk-1,ieblk(jc,jb,2)) + &
                p_int%e_bln_c_s(jc,3,jb)*z_w_concorr_me(ieidx(jc,jb,3),jk-1,ieblk(jc,jb,3))
              z_w_concorr_mc_m0 =  &
                p_int%e_bln_c_s(jc,1,jb)*z_w_concorr_me(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
                p_int%e_bln_c_s(jc,2,jb)*z_w_concorr_me(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
                p_int%e_bln_c_s(jc,3,jb)*z_w_concorr_me(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))
              p_nh%diag%w_concorr_c(jc,jk,jb) =                                &
                p_nh%metrics%wgtfac_c(jc,jk,jb)*z_w_concorr_mc_m0 +        &
                (1._vp - p_nh%metrics%wgtfac_c(jc,jk,jb))*z_w_concorr_mc_m1
            ENDDO
          ENDDO
          !$ACC END PARALLEL

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR PRIVATE(z_w_concorr_mc_m2, z_w_concorr_mc_m1, z_w_concorr_mc_m0)
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            ! Interpolate contravariant correction to cell centers...
            z_w_concorr_mc_m2 =  &
              p_int%e_bln_c_s(jc,1,jb)*z_w_concorr_me(ieidx(jc,jb,1),nlev-2,ieblk(jc,jb,1)) + &
              p_int%e_bln_c_s(jc,2,jb)*z_w_concorr_me(ieidx(jc,jb,2),nlev-2,ieblk(jc,jb,2)) + &
              p_int%e_bln_c_s(jc,3,jb)*z_w_concorr_me(ieidx(jc,jb,3),nlev-2,ieblk(jc,jb,3))

            z_w_concorr_mc_m1 =  &
              p_int%e_bln_c_s(jc,1,jb)*z_w_concorr_me(ieidx(jc,jb,1),nlev-1,ieblk(jc,jb,1)) + &
              p_int%e_bln_c_s(jc,2,jb)*z_w_concorr_me(ieidx(jc,jb,2),nlev-1,ieblk(jc,jb,2)) + &
              p_int%e_bln_c_s(jc,3,jb)*z_w_concorr_me(ieidx(jc,jb,3),nlev-1,ieblk(jc,jb,3))

            z_w_concorr_mc_m0   =  &
              p_int%e_bln_c_s(jc,1,jb)*z_w_concorr_me(ieidx(jc,jb,1),nlev,ieblk(jc,jb,1)) + &
              p_int%e_bln_c_s(jc,2,jb)*z_w_concorr_me(ieidx(jc,jb,2),nlev,ieblk(jc,jb,2)) + &
              p_int%e_bln_c_s(jc,3,jb)*z_w_concorr_me(ieidx(jc,jb,3),nlev,ieblk(jc,jb,3))

            p_nh%diag%w_concorr_c(jc,nlevp1,jb) =                         &
              p_nh%metrics%wgtfacq_c(jc,1,jb)*z_w_concorr_mc_m0 +         &
              p_nh%metrics%wgtfacq_c(jc,2,jb)*z_w_concorr_mc_m1 +       &
              p_nh%metrics%wgtfacq_c(jc,3,jb)*z_w_concorr_mc_m2
          ENDDO
          !$ACC END PARALLEL

        ENDDO
#else
!
! OMP-only code
!
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,z_w_concorr_mc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Interpolate contravariant correction to cell centers...
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
!DIR$ IVDEP
            DO jk = nflatlev(jg), nlev
#else
          DO jk = nflatlev(jg), nlev
            DO jc = i_startidx, i_endidx
#endif

              z_w_concorr_mc(jc,jk) =  &
                p_int%e_bln_c_s(jc,1,jb)*z_w_concorr_me(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
                p_int%e_bln_c_s(jc,2,jb)*z_w_concorr_me(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
                p_int%e_bln_c_s(jc,3,jb)*z_w_concorr_me(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))

            ENDDO
          ENDDO

          ! ... and to interface levels
          DO jk = nflatlev(jg)+1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              p_nh%diag%w_concorr_c(jc,jk,jb) =                                &
                p_nh%metrics%wgtfac_c(jc,jk,jb)*z_w_concorr_mc(jc,jk) +        &
               (1._vp - p_nh%metrics%wgtfac_c(jc,jk,jb))*z_w_concorr_mc(jc,jk-1)
            ENDDO
          ENDDO
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_nh%diag%w_concorr_c(jc,nlevp1,jb) =                         &
              p_nh%metrics%wgtfacq_c(jc,1,jb)*z_w_concorr_mc(jc,nlev) +   &
              p_nh%metrics%wgtfacq_c(jc,2,jb)*z_w_concorr_mc(jc,nlev-1) + &
              p_nh%metrics%wgtfacq_c(jc,3,jb)*z_w_concorr_mc(jc,nlev-2)
          ENDDO

        ENDDO
!$OMP END DO
#endif
      ENDIF
!$OMP END PARALLEL

      IF (timers_level > 5) THEN
        CALL timer_stop(timer_solve_nh_edgecomp)
        CALL timer_start(timer_solve_nh_vimpl)
      ENDIF


!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk,jk_start)

      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

      IF (l_vert_nested) THEN
        jk_start = 2
      ELSE
        jk_start = 1
      ENDIF

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,z_w_expl,z_contr_w_fl_l,z_contr_wq_fl_l,z_rho_expl,z_exner_expl, &
!$OMP   z_a,z_b,z_c,z_g,z_q,z_alpha,z_beta,z_beta_q,z_gamma,ic,z_flxdiv_mass,z_flxdiv_theta,z_flxdiv_vol  ) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! horizontal divergences of rho and rhotheta are inlined and processed in one step for efficiency

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
!DIR$ IVDEP, PREFERVECTOR
          DO jk = 1, nlev
#else
!$NEC outerloop_unroll(8)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
#endif
            z_flxdiv_mass(jc,jk) =  p_nh%metrics%deepatmo_divh_mc(jk) * (                         &
              p_nh%diag%mass_fl_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) + &
              p_nh%diag%mass_fl_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) + &
              p_nh%diag%mass_fl_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb) )

            z_flxdiv_theta(jc,jk) = p_nh%metrics%deepatmo_divh_mc(jk) * (                    &
              z_theta_v_fl_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) + &
              z_theta_v_fl_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) + &
              z_theta_v_fl_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb) )
          END DO
        END DO
        !$ACC END PARALLEL

        IF (lapply_moist_thdyn) THEN
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1, nlev
#else
!$NEC outerloop_unroll(8)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
#endif
              z_flxdiv_vol(jc,jk) = p_nh%metrics%deepatmo_divh_mc(jk) * (                  &
                z_vol_fl_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) + &
                z_vol_fl_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) + &
                z_vol_fl_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb)   )
            END DO
          END DO
          !$ACC END PARALLEL
        ENDIF

        ! upper boundary conditions for rho_ic and theta_v_ic in the case of vertical nesting
        !
        ! kept constant during predictor/corrector step, and linearly interpolated for 
        ! each dynamics substep. 
        ! Hence, copying them every dynamics substep during the predictor step (istep=1) is sufficient. 
        IF (l_vert_nested .AND. istep == 1) THEN
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            p_nh%diag%theta_v_ic(jc,1,jb) = p_nh%diag%theta_v_ic_ubc(jc,jb,1)  &
              &                           + dt_linintp_ubc * p_nh%diag%theta_v_ic_ubc(jc,jb,2)

            p_nh%diag%rho_ic(jc,1,jb) = p_nh%diag%rho_ic_ubc(jc,jb,1)  &
              &                       + dt_linintp_ubc * p_nh%diag%rho_ic_ubc(jc,jb,2)
 
            z_mflx_top(jc,jb) = p_nh%diag%mflx_ic_ubc(jc,jb,1)  &
              &               + dt_linintp_ubc * p_nh%diag%mflx_ic_ubc(jc,jb,2)

          ENDDO
          !$ACC END PARALLEL
        ENDIF

        ! Start of vertically implicit solver part for sound-wave terms;
        ! advective terms and gravity-wave terms are treated explicitly
        !
        IF (istep == 2 .AND. (itime_scheme >= 4)) THEN

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 2, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              ! explicit part for w - use temporally averaged advection terms for better numerical stability
              ! the explicit weight for the pressure-gradient term is already included in z_th_ddz_exner_c
              z_w_expl(jc,jk) = p_nh%prog(nnow)%w(jc,jk,jb) + dtime *   &
                (wgt_nnow_vel*p_nh%diag%ddt_w_adv_pc(jc,jk,jb,ntl1) +   &
                 wgt_nnew_vel*p_nh%diag%ddt_w_adv_pc(jc,jk,jb,ntl2)     &
                 -cpd*z_th_ddz_exner_c(jc,jk,jb) )

              ! contravariant vertical velocity times density for explicit part
              z_contr_wq_fl_l(jc,jk) = -p_nh%diag%w_concorr_c(jc,jk,jb) + &
                p_nh%metrics%vwind_expl_wgt(jc,jb) * p_nh%prog(nnow)%w(jc,jk,jb) 
              z_contr_w_fl_l(jc,jk) = p_nh%diag%rho_ic(jc,jk,jb) * z_contr_wq_fl_l(jc,jk)

            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ELSE

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 2, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              ! explicit part for w
              z_w_expl(jc,jk) = p_nh%prog(nnow)%w(jc,jk,jb) + dtime *                &
                (p_nh%diag%ddt_w_adv_pc(jc,jk,jb,ntl1)-cpd*z_th_ddz_exner_c(jc,jk,jb))

              ! contravariant vertical velocity times density for explicit part
              z_contr_w_fl_l(jc,jk) = p_nh%diag%rho_ic(jc,jk,jb)*(-p_nh%diag%w_concorr_c(jc,jk,jb) &
                + p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,jk,jb) )

            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        ! Solver coefficients
        IF (lapply_moist_thdyn) THEN
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              z_beta_q(jc,jk) = ( dtime * rd                    &
                       * p_nh%diag%chi_q(jc,jk,jb)              &
                       * p_nh%prog(nnow)%exner(jc,jk,jb)        &
                       * p_nh%metrics%inv_ddqz_z_full(jc,jk,jb) ) / cvd
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            z_beta(jc,jk)=dtime*rd*p_nh%prog(nnow)%exner(jc,jk,jb) /                 &
              (cvd*p_nh%prog(nnow)%rho(jc,jk,jb)*p_nh%prog(nnow)%theta_v(jc,jk,jb)) * &
              p_nh%metrics%inv_ddqz_z_full(jc,jk,jb)

            z_alpha(jc,jk)= p_nh%metrics%vwind_impl_wgt(jc,jb)*         &
              &  p_nh%diag%theta_v_ic(jc,jk,jb)*p_nh%diag%rho_ic(jc,jk,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL


        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          z_alpha(jc,nlevp1) = 0.0_wp
          !
          ! Note: z_q is used in the tridiagonal matrix solver for w below.
          !       z_q(1) is always zero, irrespective of w(1)=0 or w(1)/=0
          !       z_q(1)=0 is equivalent to cp(slev)=c(slev)/b(slev) in mo_math_utilities:tdma_solver_vec
          z_q(jc,1) = 0._vp
        ENDDO
        !$ACC END PARALLEL


        ! upper boundary condition for w
        ! interpolated from parent domain in case of vertical nesting, and
        ! rigid lid otherwise
        IF (.NOT. l_vert_nested) THEN
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            p_nh%prog(nnew)%w(jc,1,jb) = 0._wp
            z_contr_wq_fl_l(jc,1)      = 0._wp
            z_contr_w_fl_l(jc,1)       = 0._wp
          ENDDO
          !$ACC END PARALLEL
        ELSE  ! l_vert_nested
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            ! UBC for w: horizontally interpolated from the parent interface level,
            !            and linearly interpolated in time.
            p_nh%prog(nnew)%w(jc,1,jb) = p_nh%diag%w_ubc(jc,jb,1)  &
              &                        + dt_linintp_ubc_nnew * p_nh%diag%w_ubc(jc,jb,2)
            !
            z_contr_w_fl_l(jc,1) = z_mflx_top(jc,jb) * p_nh%metrics%vwind_expl_wgt(jc,jb)
            z_contr_wq_fl_l(jc,1) = z_contr_w_fl_l(jc,1) / p_nh%diag%rho_ic(jc,1,jb)
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        ! lower boundary condition for w, consistent with contravariant correction
        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,nlevp1,jb) = p_nh%diag%w_concorr_c(jc,nlevp1,jb)
          z_contr_wq_fl_l(jc,nlevp1)      = 0.0_wp
          z_contr_w_fl_l(jc,nlevp1)       = 0.0_wp
        ENDDO
        !$ACC END PARALLEL


        ! Explicit parts of density and Exner pressure
        !
        ! Top level first
        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_rho_expl(jc,1)=        p_nh%prog(nnow)%rho(jc,1,jb)     & 
            &        -dtime*p_nh%metrics%inv_ddqz_z_full(jc,1,jb)   &
            &                  *(z_flxdiv_mass(jc,1)                &
            &                  +z_contr_w_fl_l(jc,1   ) *           & 
            &                  p_nh%metrics%deepatmo_divzU_mc(1)    &  
            &                  -z_contr_w_fl_l(jc,2   ) *           &   
            &                  p_nh%metrics%deepatmo_divzL_mc(1) )

          z_exner_expl(jc,1)=     p_nh%diag%exner_pr(jc,1,jb)        &
            &      -z_beta (jc,1)*(z_flxdiv_theta(jc,1)              &
            & +p_nh%diag%theta_v_ic(jc,1,jb)*z_contr_w_fl_l(jc,1)  * &
            & p_nh%metrics%deepatmo_divzU_mc(1)                      &
            & -p_nh%diag%theta_v_ic(jc,2,jb)*z_contr_w_fl_l(jc,2)  * &
            & p_nh%metrics%deepatmo_divzL_mc(1) )                    & 
            & +dtime*p_nh%diag%ddt_exner_phy(jc,1,jb)
        ENDDO
        !$ACC END PARALLEL

        ! Other levels
        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 2, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            z_rho_expl(jc,jk)=       p_nh%prog(nnow)%rho(jc,jk  ,jb)   &
              &        -dtime*p_nh%metrics%inv_ddqz_z_full(jc,jk  ,jb) &
              &                     *(z_flxdiv_mass(jc,jk     )        &
              &                     +z_contr_w_fl_l(jc,jk     ) *      &
              &                     p_nh%metrics%deepatmo_divzU_mc(jk) &
              &                     -z_contr_w_fl_l(jc,jk+1   ) *      & 
              &                     p_nh%metrics%deepatmo_divzL_mc(jk) )

            z_exner_expl(jc,jk)=    p_nh%diag%exner_pr(jc,jk,jb) - z_beta(jc,jk) &
              &                             *(z_flxdiv_theta(jc,jk)              &
              &   +p_nh%diag%theta_v_ic(jc,jk  ,jb)*z_contr_w_fl_l(jc,jk  ) *    & 
              &   p_nh%metrics%deepatmo_divzU_mc(jk)                             &
              &   -p_nh%diag%theta_v_ic(jc,jk+1,jb)*z_contr_w_fl_l(jc,jk+1) *    & 
              &   p_nh%metrics%deepatmo_divzL_mc(jk) )                           &
              &   +dtime*p_nh%diag%ddt_exner_phy(jc,jk,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        IF (lapply_moist_thdyn) THEN
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx

              z_exner_expl(jc,jk) = z_exner_expl(jc,jk) + z_beta_q(jc,jk) *       &
                ( z_flxdiv_vol(jc,jk)                                             &
                  + z_contr_wq_fl_l(jc,jk  ) * p_nh%metrics%deepatmo_divzU_mc(jk) &
                  - z_contr_wq_fl_l(jc,jk+1) * p_nh%metrics%deepatmo_divzL_mc(jk) )

            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        IF (is_iau_active) THEN ! add analysis increments from data assimilation to density and exner pressure
          
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              z_rho_expl(jc,jk)   = z_rho_expl(jc,jk)   + iau_wgt_dyn*p_nh%diag%rho_incr(jc,jk,jb)
              z_exner_expl(jc,jk) = z_exner_expl(jc,jk) + iau_wgt_dyn*p_nh%diag%exner_incr(jc,jk,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        !
        ! Solve tridiagonal matrix for w
        !
! TODO: not parallelized
        IF (.NOT. lapply_moist_thdyn) THEN
        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP SEQ
        DO jk = 2, nlev
!DIR$ IVDEP
!$NEC ivdep
          !$ACC LOOP GANG VECTOR PRIVATE(z_gamma, z_a, z_c, z_b, z_g)
          DO jc = i_startidx, i_endidx
            z_gamma = dtime*cpd*p_nh%metrics%vwind_impl_wgt(jc,jb)*    &
              p_nh%diag%theta_v_ic(jc,jk,jb)/p_nh%metrics%ddqz_z_half(jc,jk,jb)
            z_a  = -z_gamma*z_beta(jc,jk-1)*z_alpha(jc,jk-1)*p_nh%metrics%deepatmo_divzU_mc(jk-1)
            z_c = -z_gamma*z_beta(jc,jk  )*z_alpha(jc,jk+1)*p_nh%metrics%deepatmo_divzL_mc(jk)
            z_b = 1.0_vp+z_gamma*z_alpha(jc,jk) &
              *(z_beta(jc,jk-1)*p_nh%metrics%deepatmo_divzL_mc(jk-1)+z_beta(jc,jk)*p_nh%metrics%deepatmo_divzU_mc(jk))
            z_g = 1.0_vp/(z_b+z_a*z_q(jc,jk-1))
            z_q(jc,jk) = - z_c*z_g
            p_nh%prog(nnew)%w(jc,jk,jb) = z_w_expl(jc,jk) - z_gamma  &
              &      *(z_exner_expl(jc,jk-1)-z_exner_expl(jc,jk))
            p_nh%prog(nnew)%w(jc,jk,jb) = (p_nh%prog(nnew)%w(jc,jk,jb)  &
              -z_a*p_nh%prog(nnew)%w(jc,jk-1,jb))*z_g
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        ELSE
        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP SEQ
        DO jk = 2, nlev
!$NEC ivdep
          !$ACC LOOP GANG VECTOR PRIVATE(z_gamma, z_a, z_c, z_b, z_g)
          DO jc = i_startidx, i_endidx
            z_gamma = dtime * cpd * p_nh%metrics%vwind_impl_wgt(jc,jb) * &
              p_nh%diag%theta_v_ic(jc,jk,jb) / p_nh%metrics%ddqz_z_half(jc,jk,jb)
            z_a = -z_gamma *                                               &
              ( z_beta  (jc,jk-1) * z_alpha(jc,jk-1) -                     &
                z_beta_q(jc,jk-1) * p_nh%metrics%vwind_impl_wgt(jc,jb) ) * &
              p_nh%metrics%deepatmo_divzU_mc(jk-1)
            z_c = -z_gamma *                                               &
              ( z_beta  (jc,jk  ) * z_alpha(jc,jk+1) -                     &
                z_beta_q(jc,jk  ) * p_nh%metrics%vwind_impl_wgt(jc,jb) ) * &
              p_nh%metrics%deepatmo_divzL_mc(jk  )
            z_b = 1.0_vp + z_gamma *                                         &
              ( z_alpha(jc,jk) *                                             &
                ( z_beta(jc,jk-1) * p_nh%metrics%deepatmo_divzL_mc(jk-1) +   &
                  z_beta(jc,jk  ) * p_nh%metrics%deepatmo_divzU_mc(jk  ) ) - &
                p_nh%metrics%vwind_impl_wgt(jc,jb) * &
                ( z_beta_q(jc,jk-1) * p_nh%metrics%deepatmo_divzL_mc(jk-1) + &
                  z_beta_q(jc,jk  ) * p_nh%metrics%deepatmo_divzU_mc(jk  ) ) )
            z_g = 1.0_vp / ( z_b + z_a * z_q(jc,jk-1) )
            z_q(jc,jk) = -z_c * z_g
            p_nh%prog(nnew)%w(jc,jk,jb) = z_w_expl(jc,jk) - z_gamma * &
              ( z_exner_expl(jc,jk-1) - z_exner_expl(jc,jk) )
            p_nh%prog(nnew)%w(jc,jk,jb) =           &
              ( p_nh%prog(nnew)%w(jc,jk,jb) -       &
                z_a * p_nh%prog(nnew)%w(jc,jk-1,jb) ) * z_g
          ENDDO
        ENDDO
        !$ACC END PARALLEL
        ENDIF

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP SEQ
        DO jk = nlev-1, 2, -1
!DIR$ IVDEP
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnew)%w(jc,jk,jb)&
              &             +p_nh%prog(nnew)%w(jc,jk+1,jb)*z_q(jc,jk)
          ENDDO
        ENDDO
        !$ACC END PARALLEL


        ! Rayleigh damping mechanism (Klemp,Dudhia,Hassiotis: MWR136,pp.3987-4004)
        !
        IF ( rayleigh_type == RAYLEIGH_KLEMP ) THEN

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 2, nrdmax(jg)
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              p_nh%prog(nnew)%w(jc,jk,jb) = z_raylfac(jk)*p_nh%prog(nnew)%w(jc,jk,jb) +    &
                                            (1._wp-z_raylfac(jk))*p_nh%prog(nnew)%w(jc,1,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ! Classic Rayleigh damping mechanism for w (requires reference state !!)
        !
        ELSE IF ( rayleigh_type == RAYLEIGH_CLASSIC ) THEN

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 2, nrdmax(jg)
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnew)%w(jc,jk,jb)       &
                &                         - dtime*p_nh%metrics%rayleigh_w(jk) &
                &                         * ( p_nh%prog(nnew)%w(jc,jk,jb)     &
                &                         - p_nh%ref%w_ref(jc,jk,jb) )
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        ! Results for thermodynamic variables

        IF (lapply_moist_thdyn) THEN
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR TILE(128, 1)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx

              z_exner_expl(jc,jk) = z_exner_expl(jc,jk) +                              &
                z_beta_q(jc,jk) * p_nh%metrics%vwind_impl_wgt(jc,jb) *                 &
                ( p_nh%prog(nnew)%w(jc,jk  ,jb) * p_nh%metrics%deepatmo_divzU_mc(jk) - &
                  p_nh%prog(nnew)%w(jc,jk+1,jb) * p_nh%metrics%deepatmo_divzL_mc(jk)   )

            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR TILE(128, 1)
!$NEC outerloop_unroll(8)
        DO jk = jk_start, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            ! density
            p_nh%prog(nnew)%rho(jc,jk,jb) = z_rho_expl(jc,jk)              &
              - p_nh%metrics%vwind_impl_wgt(jc,jb)*dtime                   &
              * p_nh%metrics%inv_ddqz_z_full(jc,jk,jb)                     &
              *(p_nh%diag%rho_ic(jc,jk  ,jb)*p_nh%prog(nnew)%w(jc,jk  ,jb) &
              * p_nh%metrics%deepatmo_divzU_mc(jk)                         &
              - p_nh%diag%rho_ic(jc,jk+1,jb)*p_nh%prog(nnew)%w(jc,jk+1,jb) &
              * p_nh%metrics%deepatmo_divzL_mc(jk) )

            ! exner
            p_nh%prog(nnew)%exner(jc,jk,jb) = z_exner_expl(jc,jk) &
              + p_nh%metrics%exner_ref_mc(jc,jk,jb)-z_beta(jc,jk) &
              *(z_alpha(jc,jk  )*p_nh%prog(nnew)%w(jc,jk  ,jb)    &
              * p_nh%metrics%deepatmo_divzU_mc(jk)                &
              - z_alpha(jc,jk+1)*p_nh%prog(nnew)%w(jc,jk+1,jb)    &
              * p_nh%metrics%deepatmo_divzL_mc(jk) )

            ! theta
            p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb)*p_nh%prog(nnow)%theta_v(jc,jk,jb) &
              *( (p_nh%prog(nnew)%exner(jc,jk,jb)/p_nh%prog(nnow)%exner(jc,jk,jb)-1.0_wp) * cvd_o_rd+1.0_wp   ) &
              / p_nh%prog(nnew)%rho(jc,jk,jb)

          ENDDO
        ENDDO
        !$ACC END PARALLEL

        ! Special treatment of uppermost layer in the case of vertical nesting
        IF (l_vert_nested) THEN
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            ! density
            p_nh%prog(nnew)%rho(jc,1,jb) = z_rho_expl(jc,1)           &
              - p_nh%metrics%vwind_impl_wgt(jc,jb)*dtime              &
              * p_nh%metrics%inv_ddqz_z_full(jc,1,jb)                 &
              *(z_mflx_top(jc,jb) * p_nh%metrics%deepatmo_divzU_mc(1) &
              - p_nh%diag%rho_ic(jc,2,jb)*p_nh%prog(nnew)%w(jc,2,jb)  &
              * p_nh%metrics%deepatmo_divzL_mc(1) )

            ! exner
            p_nh%prog(nnew)%exner(jc,1,jb) = z_exner_expl(jc,1)                  &
              + p_nh%metrics%exner_ref_mc(jc,1,jb)-z_beta(jc,1)                  &
              *(p_nh%metrics%vwind_impl_wgt(jc,jb)*p_nh%diag%theta_v_ic(jc,1,jb) &
              * z_mflx_top(jc,jb) * p_nh%metrics%deepatmo_divzU_mc(1)            &
              - z_alpha(jc,2)*p_nh%prog(nnew)%w(jc,2,jb)                         &
              * p_nh%metrics%deepatmo_divzL_mc(1) )

            ! theta
            p_nh%prog(nnew)%theta_v(jc,1,jb) = p_nh%prog(nnow)%rho(jc,1,jb)*p_nh%prog(nnow)%theta_v(jc,1,jb) &
              *( (p_nh%prog(nnew)%exner(jc,1,jb)/p_nh%prog(nnow)%exner(jc,1,jb)-1.0_wp) * cvd_o_rd+1.0_wp  ) &
              /p_nh%prog(nnew)%rho(jc,1,jb)

          ENDDO
          !$ACC END PARALLEL
        ENDIF


        ! compute dw/dz for divergence damping term
        IF (istep == 1 .AND. divdamp_type >= 3) THEN

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR TILE(32, 4)
          DO jk = kstart_dd3d(jg), nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              z_dwdz_dd(jc,jk,jb) = p_nh%metrics%inv_ddqz_z_full(jc,jk,jb) *          &
                ( (p_nh%prog(nnew)%w(jc,jk,jb)-p_nh%prog(nnew)%w(jc,jk+1,jb)) -       &
                (p_nh%diag%w_concorr_c(jc,jk,jb)-p_nh%diag%w_concorr_c(jc,jk+1,jb)) )
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        ! Preparations for tracer advection
        IF (lprep_adv .AND. istep == 2) THEN
          IF (lclean_mflx) THEN 
            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
            DO jk = 1, nlev
!$NEC ivdep
              DO jc = i_startidx, i_endidx
                prep_adv%mass_flx_ic(jc,jk,jb) = 0._wp
              ENDDO
            ENDDO
            !$ACC END PARALLEL
          ENDIF
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = jk_start, nlev
!$NEC ivdep
            DO jc = i_startidx, i_endidx
              prep_adv%mass_flx_ic(jc,jk,jb) = prep_adv%mass_flx_ic(jc,jk,jb) + r_nsubsteps * ( z_contr_w_fl_l(jc,jk) + &
                p_nh%diag%rho_ic(jc,jk,jb) * p_nh%metrics%vwind_impl_wgt(jc,jb) * p_nh%prog(nnew)%w(jc,jk,jb) )
            ENDDO
          ENDDO
          !$ACC END PARALLEL
          IF (l_vert_nested) THEN
            ! Use mass flux which has been interpolated to the upper nest boundary.
            ! This mass flux is also seen by the mass continuity equation (rho).
            ! Hence, by using the same mass flux for the tracer mass continuity equations,
            ! consistency with continuity (CWC) is ensured.
            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR
            DO jc = i_startidx, i_endidx
              prep_adv%mass_flx_ic(jc,1,jb) = prep_adv%mass_flx_ic(jc,1,jb) + &
                r_nsubsteps * z_mflx_top(jc,jb)
            ENDDO
            !$ACC END PARALLEL
          ENDIF
        ENDIF

        ! store dynamical part of exner time increment in exner_dyn_incr
        ! the conversion into a temperature tendency is done in the NWP interface
        IF (istep == 1 .AND. idyn_timestep == 1) THEN

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              p_nh%diag%exner_dyn_incr(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ELSE IF (istep == 2 .AND. idyn_timestep == ndyn_substeps_var(jg)) THEN
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              p_nh%diag%exner_dyn_incr(jc,jk,jb) = p_nh%prog(nnew)%exner(jc,jk,jb) - &
               (p_nh%diag%exner_dyn_incr(jc,jk,jb) + ndyn_substeps_var(jg)*dtime*p_nh%diag%ddt_exner_phy(jc,jk,jb))
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        IF (istep == 2 .AND. l_child_vertnest) THEN
          ! Store values at nest interface levels
!DIR$ IVDEP
          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx

            p_nh%diag%w_int(jc,jb,idyn_timestep) =  &
              0.5_wp*(p_nh%prog(nnow)%w(jc,nshift,jb) + p_nh%prog(nnew)%w(jc,nshift,jb))

            p_nh%diag%theta_v_ic_int(jc,jb,idyn_timestep) = p_nh%diag%theta_v_ic(jc,nshift,jb)

            p_nh%diag%rho_ic_int(jc,jb,idyn_timestep) =  p_nh%diag%rho_ic(jc,nshift,jb)

            p_nh%diag%mflx_ic_int(jc,jb,idyn_timestep) = p_nh%diag%rho_ic(jc,nshift,jb) * &
              (p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,nshift,jb) + &
              p_nh%metrics%vwind_impl_wgt(jc,jb)*p_nh%prog(nnew)%w(jc,nshift,jb))
          ENDDO
          !$ACC END PARALLEL
        ENDIF

      ENDDO
!$OMP END DO

      ! Boundary update in case of nesting
      IF (l_limited_area .OR. jg > 1) THEN

        rl_start = 1
        rl_end   = grf_bdywidth_c

        i_startblk = p_patch%cells%start_block(rl_start)
        i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! non-MPI-parallelized (serial) case
          IF (istep == 1 .AND. my_process_is_mpi_all_seq() ) THEN

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
            DO jk = 1, nlev
#if __INTEL_COMPILER != 1400 || __INTEL_COMPILER_UPDATE != 3
!DIR$ IVDEP
#endif
              DO jc = i_startidx, i_endidx

                p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb) + &
                  dtime*p_nh%diag%grf_tend_rho(jc,jk,jb)

                p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnow)%theta_v(jc,jk,jb) + &
                  dtime*p_nh%diag%grf_tend_thv(jc,jk,jb)

                ! Diagnose exner from rho*theta
                p_nh%prog(nnew)%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
                  p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)))

                p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnow)%w(jc,jk,jb) + &
                  dtime*p_nh%diag%grf_tend_w(jc,jk,jb)

              ENDDO
            ENDDO
            !$ACC END PARALLEL

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR
            DO jc = i_startidx, i_endidx
              p_nh%prog(nnew)%w(jc,nlevp1,jb) = p_nh%prog(nnow)%w(jc,nlevp1,jb) + &
                dtime*p_nh%diag%grf_tend_w(jc,nlevp1,jb)
            ENDDO
            !$ACC END PARALLEL

          ELSE IF (istep == 1 ) THEN

            ! In the MPI-parallelized case, only rho and w are updated here,
            ! and theta_v is preliminarily stored on exner in order to save
            ! halo communications

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
            DO jk = 1, nlev
#if __INTEL_COMPILER != 1400 || __INTEL_COMPILER_UPDATE != 3
!DIR$ IVDEP
#endif
              DO jc = i_startidx, i_endidx

                p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb) + &
                  dtime*p_nh%diag%grf_tend_rho(jc,jk,jb)

                ! *** Storing theta_v on exner is done to save MPI communications ***
                ! DO NOT TOUCH THIS!
                p_nh%prog(nnew)%exner(jc,jk,jb) = p_nh%prog(nnow)%theta_v(jc,jk,jb) + &
                  dtime*p_nh%diag%grf_tend_thv(jc,jk,jb)

                p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnow)%w(jc,jk,jb) + &
                  dtime*p_nh%diag%grf_tend_w(jc,jk,jb)

              ENDDO
            ENDDO
            !$ACC END PARALLEL

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR
            DO jc = i_startidx, i_endidx
              p_nh%prog(nnew)%w(jc,nlevp1,jb) = p_nh%prog(nnow)%w(jc,nlevp1,jb) + &
                dtime*p_nh%diag%grf_tend_w(jc,nlevp1,jb)
            ENDDO
            !$ACC END PARALLEL

          ENDIF

          ! compute dw/dz for divergence damping term
          IF (istep == 1 .AND. divdamp_type >= 3) THEN

            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR TILE(32, 4)
            DO jk = kstart_dd3d(jg), nlev
!DIR$ IVDEP
              DO jc = i_startidx, i_endidx
                z_dwdz_dd(jc,jk,jb) = p_nh%metrics%inv_ddqz_z_full(jc,jk,jb) *          &
                  ( (p_nh%prog(nnew)%w(jc,jk,jb)-p_nh%prog(nnew)%w(jc,jk+1,jb)) -       &
                  (p_nh%diag%w_concorr_c(jc,jk,jb)-p_nh%diag%w_concorr_c(jc,jk+1,jb)) )
              ENDDO
            ENDDO
            !$ACC END PARALLEL
          ENDIF

          ! Preparations for tracer advection
          !
          ! Note that the vertical mass flux at nest boundary points is required in case that 
          ! vertical tracer transport precedes horizontal tracer transport.
          IF (lprep_adv .AND. istep == 2) THEN
            IF (lclean_mflx) THEN
              !$ACC KERNELS IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
              prep_adv%mass_flx_ic(i_startidx:i_endidx,:,jb) = 0._wp
              !$ACC END KERNELS
            ENDIF
            !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
            DO jk = jk_start, nlev
!DIR$ IVDEP
!$NEC ivdep
              DO jc = i_startidx, i_endidx
                prep_adv%mass_flx_ic(jc,jk,jb) = prep_adv%mass_flx_ic(jc,jk,jb) + r_nsubsteps*p_nh%diag%rho_ic(jc,jk,jb)* &
                  (p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,jk,jb) +                                       &
                   p_nh%metrics%vwind_impl_wgt(jc,jb)*p_nh%prog(nnew)%w(jc,jk,jb) - p_nh%diag%w_concorr_c(jc,jk,jb) )
              ENDDO
            ENDDO
            !$ACC END PARALLEL
            IF (l_vert_nested) THEN
              !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
              !$ACC LOOP GANG VECTOR
              DO jc = i_startidx, i_endidx
                prep_adv%mass_flx_ic(jc,1,jb) = prep_adv%mass_flx_ic(jc,1,jb) + &
                  r_nsubsteps * (p_nh%diag%mflx_ic_ubc(jc,jb,1)                 &
                  + dt_linintp_ubc * p_nh%diag%mflx_ic_ubc(jc,jb,2))
              ENDDO
              !$ACC END PARALLEL
            ENDIF
          ENDIF

        ENDDO
!$OMP END DO

      ENDIF

!$OMP END PARALLEL

      !-------------------------
      ! communication phase

      IF (timers_level > 5) THEN
        CALL timer_stop(timer_solve_nh_vimpl)
        CALL timer_start(timer_solve_nh_exch)
      ENDIF

      IF (istep == 1) THEN
        IF (divdamp_type >= 3) THEN
          ! Synchronize w and vertical contribution to divergence damping
#ifdef __MIXED_PRECISION
          CALL sync_patch_array_mult_mp(SYNC_C,p_patch,1,1,p_nh%prog(nnew)%w,f3din1_sp=z_dwdz_dd, &
               &                        opt_varname="w_nnew and z_dwdz_dd")
#else
          CALL sync_patch_array_mult(SYNC_C,p_patch,2,p_nh%prog(nnew)%w,z_dwdz_dd, &
               &                     opt_varname="w_nnew and z_dwdz_dd")
#endif
        ELSE
          ! Only w needs to be synchronized
          CALL sync_patch_array(SYNC_C,p_patch,p_nh%prog(nnew)%w,opt_varname="w_nnew")
        ENDIF
      ELSE ! istep = 2: synchronize all prognostic variables
        CALL sync_patch_array_mult(SYNC_C,p_patch,3,p_nh%prog(nnew)%rho, &
          p_nh%prog(nnew)%exner,p_nh%prog(nnew)%w,opt_varname="rho, exner, w_nnew")
      ENDIF

      IF (timers_level > 5) CALL timer_stop(timer_solve_nh_exch)

      ! end communication phase
      !-------------------------

    ENDDO ! istep-loop


    ! The remaining computations are needed for MPI-parallelized applications only
    IF ( .NOT. my_process_is_mpi_all_seq() ) THEN

! OpenMP directives are commented for the NEC because the overhead is too large
#if !defined( __SX__ ) 
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
#endif
      IF (l_limited_area .OR. jg > 1) THEN

        ! Index list over halo points lying in the boundary interpolation zone
        ! Note: this list typically contains at most 10 grid points

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG
#ifndef __SX__
!$OMP DO PRIVATE(jb,ic,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
#endif
        DO ic = 1, p_nh%metrics%bdy_halo_c_dim

          jb = p_nh%metrics%bdy_halo_c_blk(ic)
          jc = p_nh%metrics%bdy_halo_c_idx(ic)
!DIR$ IVDEP
          !$ACC LOOP VECTOR
          DO jk = 1, nlev
            p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnew)%exner(jc,jk,jb)

            ! Diagnose exner from rho*theta
            p_nh%prog(nnew)%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
              p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)))

          ENDDO
        ENDDO
        !$ACC END PARALLEL
#ifndef __SX__
!$OMP END DO
#endif

        rl_start = 1
        rl_end   = grf_bdywidth_c

        i_startblk = p_patch%cells%start_block(rl_start)
        i_endblk   = p_patch%cells%end_block(rl_end)

#ifndef __SX__
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
#endif
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnew)%exner(jc,jk,jb)

              ! Diagnose exner from rhotheta
              p_nh%prog(nnew)%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
                p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)))

            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDDO
#ifndef __SX__
!$OMP END DO
#endif
      ENDIF

      rl_start = min_rlcell_int - 1
      rl_end   = min_rlcell

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

#ifndef __SX__
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
#endif
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)

#ifdef __LOOP_EXCHANGE
        !$ACC LOOP GANG
        DO jc = i_startidx, i_endidx
          IF (p_nh%metrics%mask_prog_halo_c(jc,jb)) THEN
!DIR$ IVDEP
            !$ACC LOOP VECTOR
            DO jk = 1, nlev
#else
        !$ACC LOOP GANG VECTOR TILE(32, 4)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            IF (p_nh%metrics%mask_prog_halo_c(jc,jb)) THEN
#endif
              p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb)*p_nh%prog(nnow)%theta_v(jc,jk,jb) &
                *( (p_nh%prog(nnew)%exner(jc,jk,jb)/p_nh%prog(nnow)%exner(jc,jk,jb)-1.0_wp) * cvd_o_rd+1.0_wp   ) &
                / p_nh%prog(nnew)%rho(jc,jk,jb)

#ifdef __LOOP_EXCHANGE
            ENDDO
          ENDIF
#else
            ENDIF
          ENDDO
#endif
        ENDDO
        !$ACC END PARALLEL

      ENDDO
#ifndef __SX__
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

    ENDIF  ! .NOT. my_process_is_mpi_all_seq()

    IF (ltimer) CALL timer_stop(timer_solve_nh)

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE solve_nh

#ifdef _OPENACC

     SUBROUTINE h2d_solve_nonhydro( nnow, jstep, jg, grf_intmethod_e, lprep_adv, l_vert_nested, is_iau_active, &
                                    p_nh, prep_adv )

       INTEGER, INTENT(IN)       :: nnow, jstep, jg, grf_intmethod_e
       LOGICAL, INTENT(IN)       :: l_vert_nested, lprep_adv, is_iau_active

       TYPE(t_nh_state),            INTENT(INOUT) :: p_nh
       TYPE(t_prepare_adv), TARGET, INTENT(INOUT) :: prep_adv

       REAL(wp), DIMENSION(:,:,:),   POINTER  :: exner_tmp, rho_tmp, theta_v_tmp, vn_tmp, w_tmp                 ! p_prog  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: vn_ie_ubc_tmp                                                 ! p_diag  WP 2D
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: w_ubc_tmp, mflx_ic_ubc_tmp, theta_v_ic_ubc_tmp, rho_ic_ubc_tmp ! p_diag  WP

       REAL(wp), DIMENSION(:,:,:),   POINTER  :: theta_v_ic_tmp, rho_ic_tmp                                     ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: mass_fl_e_tmp, exner_pr_tmp                                    ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: grf_bdy_mflx_tmp                                               ! p_diag  WP

       REAL(vp), DIMENSION(:,:,:),   POINTER  :: vt_tmp, vn_ie_tmp, w_concorr_c_tmp, ddt_exner_phy_tmp          ! p_diag  VP
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: exner_dyn_incr_tmp                                             ! p_diag  VP 
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: ddt_vn_phy_tmp                                                 ! p_diag  VP

       REAL(vp), DIMENSION(:,:,:),   POINTER  :: rho_incr_tmp, exner_incr_tmp                                   ! p_diag  VP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: vn_traj_tmp, mass_flx_me_tmp, mass_flx_ic_tmp                  ! prep_adv WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: vn_ref_tmp, w_ref_tmp                                          ! p_ref   WP

       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_vn_apc_pc_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_vn_cor_pc_tmp
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_w_adv_pc_tmp

       REAL(wp), DIMENSION(:,:,:),   POINTER  :: ddt_vn_dyn_tmp, ddt_vn_dmp_tmp, ddt_vn_adv_tmp, ddt_vn_cor_tmp ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: ddt_vn_pgr_tmp, ddt_vn_phd_tmp, ddt_vn_iau_tmp, ddt_vn_ray_tmp ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: ddt_vn_grf_tmp                                                 ! p_diag  WP

! p_patch:
!            p_patch%cells:   edge_idx/blk
!            p_patch%edges:   cell_idx/blk, vertex_idx/blk, quad_idx/blk, 
!                             primal/dual_normal_cell, inv_primal/dual_edge_length, tangent_orientation, refin_ctrl 

!
! p_nh%metrics:  vertidx_gradp, pg_vertidx, pg_edgeidx, pg_edgeblk,
!                bdy_halo_c_blk, bdy_halo_c_idx, bdy_mflx_e_blk, bdy_mflx_e_idx,
!                coeff_gradp, d_exner_dz_ref_ic, d2dexdz2_fac1_mc, 
!                ddqz_z_half, ddxn_z_full, ddxt_z_full, ddqz_z_full_e,
!                exner_exfac, exner_ref_mc, hmask_dd3d, inv_ddqz_z_full,
!                mask_prog_halo_c, nudge_e_blk, nudge_e_idx, pg_exdist,
!                rayleigh_vn, rayleigh_w, rho_ref_mc, rho_ref_me,
!                scalfac_dd3d, theta_ref_ic, theta_ref_mc, theta_ref_me,
!                vwind_expl_wgt, vwind_impl_wgt, 
!                wgtfac_c, wgtfac_e, wgtfacq_c, wgtfacq1_c, zdiff_gradp


! p_nh%prog(nnow)          All present (above)

       exner_tmp           => p_nh%prog(nnow)%exner 
       rho_tmp             => p_nh%prog(nnow)%rho
       theta_v_tmp         => p_nh%prog(nnow)%theta_v 
       vn_tmp              => p_nh%prog(nnow)%vn
       w_tmp               => p_nh%prog(nnow)%w
       !$ACC UPDATE DEVICE(exner_tmp, rho_tmp, theta_v_tmp, vn_tmp, w_tmp) ASYNC(1)

! p_nh%diag:

       rho_ic_tmp          => p_nh%diag%rho_ic
       theta_v_ic_tmp      => p_nh%diag%theta_v_ic
       !$ACC UPDATE DEVICE(rho_ic_tmp, theta_v_ic_tmp) ASYNC(1)

       vt_tmp              => p_nh%diag%vt
       vn_ie_tmp           => p_nh%diag%vn_ie
       w_concorr_c_tmp     => p_nh%diag%w_concorr_c
       !$ACC UPDATE DEVICE(vt_tmp, vn_ie_tmp, w_concorr_c_tmp) ASYNC(1)

       mass_fl_e_tmp       => p_nh%diag%mass_fl_e
       exner_pr_tmp        => p_nh%diag%exner_pr
       exner_dyn_incr_tmp  => p_nh%diag%exner_dyn_incr
       !$ACC UPDATE DEVICE(mass_fl_e_tmp, exner_pr_tmp, exner_dyn_incr_tmp) ASYNC(1)

! WS: I do not think these are necessary, but adding for completeness
       ddt_vn_apc_pc_tmp   => p_nh%diag%ddt_vn_apc_pc
       ddt_w_adv_pc_tmp    => p_nh%diag%ddt_w_adv_pc
       !$ACC UPDATE DEVICE(ddt_vn_apc_pc_tmp, ddt_w_adv_pc_tmp) ASYNC(1)
       IF (p_nh%diag%ddt_vn_adv_is_associated .OR. p_nh%diag%ddt_vn_cor_is_associated) THEN
          ddt_vn_cor_pc_tmp   => p_nh%diag%ddt_vn_cor_pc
          !$ACC UPDATE DEVICE(ddt_vn_cor_pc_tmp) ASYNC(1)
       END IF

! MAG: For completeness
       ddt_vn_dyn_tmp      => p_nh%diag%ddt_vn_dyn
       !$ACC UPDATE DEVICE(ddt_vn_dyn_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_dyn_is_associated)
       ddt_vn_dmp_tmp      => p_nh%diag%ddt_vn_dmp
       !$ACC UPDATE DEVICE(ddt_vn_dmp_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_dmp_is_associated)
       ddt_vn_adv_tmp      => p_nh%diag%ddt_vn_adv
       !$ACC UPDATE DEVICE(ddt_vn_adv_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_adv_is_associated)
       ddt_vn_cor_tmp      => p_nh%diag%ddt_vn_cor
       !$ACC UPDATE DEVICE(ddt_vn_cor_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_cor_is_associated)
       ddt_vn_pgr_tmp      => p_nh%diag%ddt_vn_pgr
       !$ACC UPDATE DEVICE(ddt_vn_pgr_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_pgr_is_associated)
       ddt_vn_phd_tmp      => p_nh%diag%ddt_vn_phd
       !$ACC UPDATE DEVICE(ddt_vn_phd_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_phd_is_associated)
       ddt_vn_iau_tmp      => p_nh%diag%ddt_vn_iau
       !$ACC UPDATE DEVICE(ddt_vn_iau_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_iau_is_associated)
       ddt_vn_ray_tmp      => p_nh%diag%ddt_vn_ray
       !$ACC UPDATE DEVICE(ddt_vn_ray_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_ray_is_associated)
       ddt_vn_grf_tmp      => p_nh%diag%ddt_vn_grf
       !$ACC UPDATE DEVICE(ddt_vn_grf_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_grf_is_associated)

       mflx_ic_ubc_tmp     => p_nh%diag%mflx_ic_ubc
       vn_ie_ubc_tmp       => p_nh%diag%vn_ie_ubc
       theta_v_ic_ubc_tmp  => p_nh%diag%theta_v_ic_ubc
       rho_ic_ubc_tmp      => p_nh%diag%rho_ic_ubc
       w_ubc_tmp           => p_nh%diag%w_ubc
       !$ACC UPDATE &
       !$ACC   DEVICE(mflx_ic_ubc_tmp, vn_ie_ubc_tmp) &
       !$ACC   DEVICE(theta_v_ic_ubc_tmp, rho_ic_ubc_tmp, w_ubc_tmp) &
       !$ACC   ASYNC(1) IF(l_vert_nested)

       ddt_exner_phy_tmp   => p_nh%diag%ddt_exner_phy
       ddt_vn_phy_tmp      => p_nh%diag%ddt_vn_phy
       !$ACC UPDATE DEVICE(ddt_exner_phy_tmp, ddt_vn_phy_tmp) ASYNC(1)

       rho_incr_tmp        => p_nh%diag%rho_incr
       exner_incr_tmp      => p_nh%diag%exner_incr
       !$ACC UPDATE DEVICE(rho_incr_tmp, exner_incr_tmp) ASYNC(1)

       grf_bdy_mflx_tmp   => p_nh%diag%grf_bdy_mflx
       !$ACC UPDATE DEVICE(grf_bdy_mflx_tmp) ASYNC(1) IF((jg > 1) .AND. (grf_intmethod_e == 6) .AND. (jstep == 0))

! prep_adv:

       vn_traj_tmp       => prep_adv%vn_traj
       mass_flx_me_tmp   => prep_adv%mass_flx_me
       mass_flx_ic_tmp   => prep_adv%mass_flx_ic
       !$ACC UPDATE DEVICE(vn_traj_tmp, mass_flx_me_tmp, mass_flx_ic_tmp) ASYNC(1) IF(lprep_adv)

! p_nh%ref:

       vn_ref_tmp          => p_nh%ref%vn_ref
       w_ref_tmp           => p_nh%ref%w_ref
       !$ACC UPDATE DEVICE(vn_ref_tmp, w_ref_tmp) ASYNC(1)

     END SUBROUTINE h2d_solve_nonhydro

     SUBROUTINE d2h_solve_nonhydro( nnew, jstep, jg, idyn_timestep, grf_intmethod_e, lsave_mflx, l_child_vertnest, &
          &                         lprep_adv, p_nh, prep_adv )

       INTEGER, INTENT(IN)       :: nnew, jstep, jg, idyn_timestep, grf_intmethod_e
       LOGICAL, INTENT(IN)       :: lsave_mflx, l_child_vertnest, lprep_adv

       TYPE(t_nh_state),            INTENT(INOUT) :: p_nh
       TYPE(t_prepare_adv), TARGET, INTENT(INOUT) :: prep_adv

       REAL(wp), DIMENSION(:,:,:),   POINTER  :: exner_tmp, rho_tmp, theta_v_tmp, vn_tmp, w_tmp                 ! p_prog  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: vn_ie_int_tmp                                                  ! p_diag  WP 2D
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: theta_v_ic_tmp, rho_ic_tmp, rho_ic_int_tmp, w_int_tmp          ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: theta_v_ic_int_tmp, grf_bdy_mflx_tmp                           ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: mass_fl_e_tmp,  mflx_ic_int_tmp, exner_pr_tmp                  ! p_diag  WP

       REAL(vp), DIMENSION(:,:,:),   POINTER  :: vt_tmp, vn_ie_tmp, w_concorr_c_tmp                             ! p_diag  VP
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: mass_fl_e_sv_tmp                                               ! p_diag  VP
       REAL(vp), DIMENSION(:,:,:),   POINTER  :: exner_dyn_incr_tmp                                             ! p_diag  VP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: vn_traj_tmp, mass_flx_me_tmp, mass_flx_ic_tmp                  ! prep_adv WP
       REAL(vp), DIMENSION(:,:,:,:), POINTER  :: ddt_vn_apc_pc_tmp, ddt_vn_cor_pc_tmp, ddt_w_adv_pc_tmp

       REAL(wp), DIMENSION(:,:,:),   POINTER  :: ddt_vn_dyn_tmp, ddt_vn_dmp_tmp, ddt_vn_adv_tmp, ddt_vn_cor_tmp ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: ddt_vn_pgr_tmp, ddt_vn_phd_tmp, ddt_vn_iau_tmp, ddt_vn_ray_tmp ! p_diag  WP
       REAL(wp), DIMENSION(:,:,:),   POINTER  :: ddt_vn_grf_tmp                                                 ! p_diag  WP

! The following code is necessary if the Dycore is to be run in isolation on the GPU
! Update all device output on host: the prognostic variables have shifted from nnow to nnew; diagnostics pointers set above

       exner_tmp           => p_nh%prog(nnew)%exner
       rho_tmp             => p_nh%prog(nnew)%rho
       theta_v_tmp         => p_nh%prog(nnew)%theta_v
       vn_tmp              => p_nh%prog(nnew)%vn
       w_tmp               => p_nh%prog(nnew)%w
       !$ACC UPDATE HOST(exner_tmp, rho_tmp, theta_v_tmp, vn_tmp, w_tmp) ASYNC(1)

       vt_tmp              => p_nh%diag%vt
       vn_ie_tmp           => p_nh%diag%vn_ie
       rho_ic_tmp          => p_nh%diag%rho_ic
       theta_v_ic_tmp      => p_nh%diag%theta_v_ic
       exner_pr_tmp        => p_nh%diag%exner_pr
       !$ACC UPDATE HOST(vt_tmp, vn_ie_tmp, rho_ic_tmp, theta_v_ic_tmp, exner_pr_tmp) ASYNC(1)

       w_concorr_c_tmp     => p_nh%diag%w_concorr_c
       mass_fl_e_tmp       => p_nh%diag%mass_fl_e
       exner_dyn_incr_tmp  => p_nh%diag%exner_dyn_incr
       !$ACC UPDATE HOST(w_concorr_c_tmp, mass_fl_e_tmp, exner_dyn_incr_tmp) ASYNC(1)

       ddt_vn_apc_pc_tmp   => p_nh%diag%ddt_vn_apc_pc
       ddt_w_adv_pc_tmp    => p_nh%diag%ddt_w_adv_pc
       !$ACC UPDATE HOST(ddt_vn_apc_pc_tmp, ddt_w_adv_pc_tmp) ASYNC(1)
       IF (p_nh%diag%ddt_vn_adv_is_associated .OR. p_nh%diag%ddt_vn_cor_is_associated) THEN
          ddt_vn_cor_pc_tmp   => p_nh%diag%ddt_vn_cor_pc
          !$ACC UPDATE HOST(ddt_vn_cor_pc_tmp) ASYNC(1)
       END IF

! MAG: For completeness
       ddt_vn_dyn_tmp      => p_nh%diag%ddt_vn_dyn
       !$ACC UPDATE HOST(ddt_vn_dyn_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_dyn_is_associated)
       ddt_vn_dmp_tmp      => p_nh%diag%ddt_vn_dmp
       !$ACC UPDATE HOST(ddt_vn_dmp_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_dmp_is_associated)
       ddt_vn_adv_tmp      => p_nh%diag%ddt_vn_adv
       !$ACC UPDATE HOST(ddt_vn_adv_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_adv_is_associated)
       ddt_vn_cor_tmp      => p_nh%diag%ddt_vn_cor
       !$ACC UPDATE HOST(ddt_vn_cor_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_cor_is_associated)
       ddt_vn_pgr_tmp      => p_nh%diag%ddt_vn_pgr
       !$ACC UPDATE HOST(ddt_vn_pgr_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_pgr_is_associated)
       ddt_vn_phd_tmp      => p_nh%diag%ddt_vn_phd
       !$ACC UPDATE HOST(ddt_vn_phd_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_phd_is_associated)
       ddt_vn_iau_tmp      => p_nh%diag%ddt_vn_iau
       !$ACC UPDATE HOST(ddt_vn_iau_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_iau_is_associated)
       ddt_vn_ray_tmp      => p_nh%diag%ddt_vn_ray
       !$ACC UPDATE HOST(ddt_vn_ray_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_ray_is_associated)
       ddt_vn_grf_tmp      => p_nh%diag%ddt_vn_grf
       !$ACC UPDATE HOST(ddt_vn_grf_tmp) ASYNC(1) IF(p_nh%diag%ddt_vn_grf_is_associated)

       mass_fl_e_sv_tmp    => p_nh%diag%mass_fl_e_sv
       !$ACC UPDATE HOST(mass_fl_e_sv_tmp) ASYNC(1) IF(lsave_mflx)

       w_int_tmp           => p_nh%diag%w_int
       mflx_ic_int_tmp     => p_nh%diag%mflx_ic_int
       theta_v_ic_int_tmp  => p_nh%diag%theta_v_ic_int
       rho_ic_int_tmp      => p_nh%diag%rho_ic_int
       !$ACC UPDATE HOST(w_int_tmp, mflx_ic_int_tmp, theta_v_ic_int_tmp, rho_ic_int_tmp) ASYNC(1) IF(l_child_vertnest)

      vn_ie_int_tmp      => p_nh%diag%vn_ie_int
      !$ACC UPDATE HOST(vn_ie_int_tmp) ASYNC(1) IF(idyn_timestep == 1 .AND. l_child_vertnest)

      grf_bdy_mflx_tmp    => p_nh%diag%grf_bdy_mflx
      !$ACC UPDATE HOST(grf_bdy_mflx_tmp) ASYNC(1) IF((jg > 1) .AND. (grf_intmethod_e == 6) .AND. (jstep == 0))

      vn_traj_tmp         => prep_adv%vn_traj
      mass_flx_me_tmp     => prep_adv%mass_flx_me
      mass_flx_ic_tmp     => prep_adv%mass_flx_ic
      !$ACC UPDATE HOST(vn_traj_tmp, mass_flx_me_tmp, mass_flx_ic_tmp) ASYNC(1) IF(lprep_adv)

      !$ACC WAIT(1)

     END SUBROUTINE d2h_solve_nonhydro

#endif

END MODULE mo_solve_nonhydro
