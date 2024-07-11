!
! This module contains the driver routines needed for initializing nested
! domains that are started sometime during the integration
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

MODULE mo_nh_init_nest_utils

  USE mo_kind,                  ONLY: wp
  USE mo_model_domain,          ONLY: t_patch, p_patch, p_patch_local_parent
  USE mo_nonhydro_types,        ONLY: t_nh_metrics, t_nh_prog, t_nh_diag
  USE mo_nonhydro_state,        ONLY: p_nh_state
  USE mo_initicon_types,        ONLY: t_initicon_state
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag
  USE mo_parallel_config,       ONLY: nproma, p_test_run
  USE mo_run_config,            ONLY: ltransport, msg_level, ntracer, iforcing
  USE mo_dynamics_config,       ONLY: nnow, nnow_rcf, nnew_rcf
  USE mo_physical_constants,    ONLY: rd, cvd_o_rd, p0ref, tmelt
  USE sfc_terra_data,           ONLY: crhosminf, cporv, cadp, csalb, ist_seawtr
  USE mo_impl_constants,        ONLY: min_rlcell, min_rlcell_int, &
    &                                 inwp, nclass_aero, ALB_SI_MISSVAL
  USE mo_grf_nudgintp,          ONLY: interpol_scal_nudging, interpol_vec_nudging
  USE mo_grf_bdyintp,           ONLY: interpol_scal_grf, interpol2_vec_grf
  USE mo_grid_config,           ONLY: lfeedback, ifeedback_type
  USE mo_exception,             ONLY: message, message_text
  USE mo_mpi,                   ONLY: my_process_is_mpi_parallel
  USE mo_communication,         ONLY: exchange_data, exchange_data_mult
  USE mo_sync,                  ONLY: sync_patch_array, sync_patch_array_mult, &
                                      SYNC_C, SYNC_E
  USE mo_intp_data_strc,        ONLY: t_int_state, p_int_state, p_int_state_local_parent
  USE mo_grf_intp_data_strc,    ONLY: t_gridref_single_state, t_gridref_state, &
                                      p_grf_state, p_grf_state_local_parent
  USE mo_loopindices,           ONLY: get_indices_c
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_fbk_start_c
  USE mo_nwp_lnd_types,         ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_lnd_nwp_config,        ONLY: ntiles_total, ntiles_water, nlev_soil, lseaice, itype_trvg, &
    &                                 llake, isub_lake, frlake_thrhld, frsea_thrhld, lprog_albsi, &
    &                                 itype_snowevap, dzsoil, frsi_min
  USE mo_initicon_config,       ONLY: icpl_da_sfcevap, icpl_da_skinc, icpl_da_sfcfric
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config, iprog_aero
  USE mo_interpol_config,       ONLY: nudge_zone_width
  USE mo_ext_data_types,        ONLY: t_external_data
  USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_cell
  USE mo_nwp_sfc_interp,        ONLY: smi_to_wsoil, wsoil_to_smi
  USE sfc_flake,                ONLY: flake_coldinit
  USE mo_input_instructions,    ONLY: t_readInstructionListPtr, kStateFailedFetch, &
    &                                 kInputSourceAnaI, kInputSourceFgAnaI, kInputSourceAna, &
    &                                 kInputSourceBoth

  IMPLICIT NONE

  PRIVATE


  REAL(wp), PARAMETER :: rd_o_cvd = 1._wp / cvd_o_rd

  PUBLIC :: initialize_nest, topo_blending_and_fbk, interpolate_scal_increments, &
            interpolate_vn_increments, interpolate_sfcana

  CONTAINS
  !-------------
  !>
  !! SUBROUTINE initialize_nest
  !!
  !! Driver routine for initializing a nested domain during runtime from the parent domain
  !!
  !! Note: the interpolation of the land variables does not yet include the
  !! land-water mask and adaptations to make optimal use of tiles
  !! The aggregate_landvars routine must be called before this routine
  !!
  SUBROUTINE initialize_nest(jg, jgc, ext_data, prm_diag, p_lnd_state)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'initialize_nest'

    INTEGER,                   INTENT(IN)    :: jg   ! parent (source) domain ID
    INTEGER,                   INTENT(IN)    :: jgc  ! child  (target) domain ID
    TYPE(t_external_data),     INTENT(IN)    :: ext_data(:)
    TYPE(t_nwp_phy_diag),      INTENT(INOUT) :: prm_diag(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)

    ! local pointers
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog
    TYPE(t_nh_prog),    POINTER     :: p_child_prog
    TYPE(t_nh_diag),    POINTER     :: p_parent_diag
    TYPE(t_nh_diag),    POINTER     :: p_child_diag
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog_rcf
    TYPE(t_nh_prog),    POINTER     :: p_child_prog_rcf
    TYPE(t_nh_metrics), POINTER     :: p_parent_metrics
    TYPE(t_nh_metrics), POINTER     :: p_child_metrics
    TYPE(t_lnd_prog),   POINTER     :: p_parent_lprog
    TYPE(t_lnd_prog),   POINTER     :: p_child_lprog
    TYPE(t_lnd_prog),   POINTER     :: p_child_lprog2
    TYPE(t_wtr_prog),   POINTER     :: p_parent_wprog
    TYPE(t_wtr_prog),   POINTER     :: p_child_wprog
    TYPE(t_lnd_diag),   POINTER     :: p_parent_ldiag
    TYPE(t_lnd_diag),   POINTER     :: p_child_ldiag
    TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
    TYPE(t_int_state), POINTER      :: p_int => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()

    TYPE(t_lnd_prog),   POINTER     :: lnd_prog

    ! local variables

    ! Indices
    INTEGER :: jb, jc, jk, jk1, jt, i_nchdom, i_chidx, i_startblk, i_endblk, &
               i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: nlev_c, nlev_p
    INTEGER :: nshift      ! difference between upper boundary of parent or feedback-parent
                           ! domain and upper boundary of child domain (in terms
                           ! of vertical levels)
    INTEGER :: num_lndvars, num_wtrvars, num_phdiagvars

    ! Local arrays for variables living on the local parent grid in the MPI case. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: vn_lp, w_lp, thv_pr_lp, rho_pr_lp, phdiag_lp, &
                                                 lndvars_lp, wtrvars_lp, aero_lp
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: tracer_lp

    ! Local arrays on the parent or child grid. These would not have to be allocatable,
    ! but the computational overhead does not matter for an initialization routine
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: thv_pr_par, rho_pr_par, lndvars_par, lndvars_chi, &
                                                 wtrvars_par, wtrvars_chi, phdiag_par, phdiag_chi
    REAL(wp), ALLOCATABLE :: tsfc_ref_p(:,:), tsfc_ref_c(:,:) ! Reference temperature at lowest level

    LOGICAL :: l_parallel, l_limit(ntracer)

    INTEGER :: i_count, ic, ist

    !-----------------------------------------------------------------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') 'Nest initialization, domain ',jg,' =>',jgc
      CALL message(routine,message_text)
    ENDIF

    l_parallel = my_process_is_mpi_parallel()

    p_parent_prog     => p_nh_state(jg)%prog(nnow(jg))
    p_child_prog      => p_nh_state(jgc)%prog(nnow(jgc))
    p_parent_diag     => p_nh_state(jg)%diag
    p_child_diag      => p_nh_state(jgc)%diag
    p_parent_prog_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
    p_child_prog_rcf  => p_nh_state(jgc)%prog(nnow_rcf(jgc))
    p_parent_metrics  => p_nh_state(jg)%metrics
    p_child_metrics   => p_nh_state(jgc)%metrics
    p_parent_lprog    => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
    p_parent_wprog    => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))
    p_child_lprog     => p_lnd_state(jgc)%prog_lnd(nnow_rcf(jgc))
    p_child_lprog2    => p_lnd_state(jgc)%prog_lnd(nnew_rcf(jgc))
    p_child_wprog     => p_lnd_state(jgc)%prog_wtr(nnow_rcf(jgc))
    p_parent_ldiag    => p_lnd_state(jg)%diag_lnd
    p_child_ldiag     => p_lnd_state(jgc)%diag_lnd
    p_pc              => p_patch(jgc)

    lnd_prog          => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))

    p_grf => p_grf_state_local_parent(jgc)
    p_int => p_int_state_local_parent(jgc)
    p_pp  => p_patch_local_parent(jgc)

    i_nchdom = MAX(1,p_patch(jg)%n_childdom)
    i_chidx  = p_pc%parent_child_index

    ! number of full levels of child domain
    nlev_c   = p_pc%nlev

    ! number of full levels of parent domain
    nlev_p   = p_patch(jg)%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift

    ! number of land and water variables to be interpolated
    ! Remark (GZ): the multi-layer snow variables are initialized afterwards in terra_multlay_init. This
    ! turned out to cause occasional conflicts with directly interpolating those variables here; thus
    ! the interpolation of the multi-layer snow fields has been completely removed from this routine
    num_lndvars = 2*nlev_soil+1+ &     ! multi-layer soil variables t_so and w_so (w_so_ice is initialized in terra_multlay_init)
                  5+12+1               ! single-layer prognostic variables + t_g, t_sk, freshsnow, t_seasfc, qv_s, plantevap, hsnow_max, 
                                       ! snow_age, t_avginc, t_sk, rh_avginc, t_wgt_avginc, vabs_avginc + aux variable for lake temp
    num_wtrvars  = 6                   ! water state fields + fr_seaice + alb_si
    num_phdiagvars = 32                ! number of physics diagnostic variables (copied from interpol_phys_grf)

    ALLOCATE(thv_pr_par  (nproma, nlev_p,      p_patch(jg)%nblks_c), &
             rho_pr_par  (nproma, nlev_p,      p_patch(jg)%nblks_c), &
             lndvars_par (nproma, num_lndvars, p_patch(jg)%nblks_c), &
             wtrvars_par (nproma, num_wtrvars, p_patch(jg)%nblks_c), &
             phdiag_par  (nproma, num_phdiagvars, p_patch(jg)%nblks_c), &
             lndvars_chi (nproma, num_lndvars, p_patch(jgc)%nblks_c),&
             wtrvars_chi (nproma, num_wtrvars, p_patch(jgc)%nblks_c),&
             phdiag_chi  (nproma, num_phdiagvars, p_patch(jgc)%nblks_c),&
             tsfc_ref_p  (nproma,              p_patch(jg)%nblks_c), &
             tsfc_ref_c  (nproma,              p_patch(jgc)%nblks_c) )


    ALLOCATE(vn_lp      (nproma, nlev_p,      p_pp%nblks_e),          &
             w_lp       (nproma, nlev_p+1,    p_pp%nblks_c),          &
             thv_pr_lp  (nproma, nlev_p,      p_pp%nblks_c),          &
             rho_pr_lp  (nproma, nlev_p,      p_pp%nblks_c),          &
             phdiag_lp  (nproma, num_phdiagvars,p_pp%nblks_c),        &
             tracer_lp  (nproma, nlev_p,      p_pp%nblks_c, ntracer), &
             aero_lp    (nproma, nclass_aero, p_pp%nblks_c),          &
             lndvars_lp (nproma, num_lndvars, p_pp%nblks_c),          &
             wtrvars_lp (nproma, num_wtrvars ,p_pp%nblks_c)           )

    IF (p_test_run) THEN
      lndvars_par  = 0._wp
      wtrvars_par  = 0._wp
      phdiag_par   = 0._wp
      lndvars_chi  = 0._wp
      wtrvars_chi  = 0._wp
      phdiag_chi   = 0._wp
    ENDIF



    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
      ! Step 0: get lake surface temperature
      !
      ! exclude nest boundary and halo points
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
      i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,i_count,ic,jk) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        ! the last element of lndvars contains the source data for lake sfc temperature;
        ! a) initialize with t_g
        lndvars_par(:,num_lndvars,jb) = lnd_prog%t_g(:,jb)

        i_count = ext_data(jg)%atm%list_lake%ncount(jb)
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data(jg)%atm%list_lake%idx(ic,jb)

          ! b) take lake sfc temperature where available
          lndvars_par(jc,num_lndvars,jb) = lnd_prog%t_g_t(jc,jb,isub_lake)

        ENDDO

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END IF



    ! Step 1: boundary interpolation
    !
    ! Note: the sepration between boundary interpolation and the interpolation
    ! of the fields in the prognostic region of the nest is needed because the
    ! existing interpolation routines are runtime-optimized for operations that
    ! are needed regularly. Boundary interpolation is executed from the parent grid
    ! to the child grid, whereas the so-called nudging interpolation works from
    ! the local parent grid to the child grid.


    ! Step 1a: fill buffer arrays

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_patch(jg)%cells%start_blk(grf_bdywidth_c+1,1)
    i_endblk   = p_patch(jg)%cells%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

      ! Compute perturbation quantities for rho and theta at parent level
      DO jk = 1, nlev_p
        DO jc = i_startidx, i_endidx
          rho_pr_par(jc,jk,jb) = &
            p_parent_prog%rho(jc,jk,jb) - p_parent_metrics%rho_ref_mc(jc,jk,jb)

          thv_pr_par(jc,jk,jb) = &
            p_parent_prog%theta_v(jc,jk,jb) - p_parent_metrics%theta_ref_mc(jc,jk,jb)
        ENDDO
      ENDDO

      ! Reference temperature at lowest level (taken as proxy for the surface)
      DO jc = i_startidx, i_endidx
        tsfc_ref_p(jc,jb) = p_parent_metrics%theta_ref_mc(jc,nlev_p,jb) * &
                            p_parent_metrics%exner_ref_mc(jc,nlev_p,jb)
      ENDDO

      IF (iforcing == inwp) THEN
        ! Collect diagnostic physics fields (the only really important ones are the precip fields)
        DO jc = i_startidx, i_endidx
          phdiag_par(jc,1,jb) = prm_diag(jg)%tot_prec(jc,jb)
          phdiag_par(jc,2,jb) = prm_diag(jg)%prec_gsp(jc,jb)
          phdiag_par(jc,3,jb) = prm_diag(jg)%prec_con(jc,jb)
          phdiag_par(jc,4,jb) = prm_diag(jg)%rain_gsp(jc,jb)
          phdiag_par(jc,5,jb) = prm_diag(jg)%snow_gsp(jc,jb)
          phdiag_par(jc,6,jb) = prm_diag(jg)%rain_con(jc,jb)
          phdiag_par(jc,7,jb) = prm_diag(jg)%snow_con(jc,jb)
          phdiag_par(jc,8,jb) = prm_diag(jg)%rain_gsp_rate(jc,jb)
          phdiag_par(jc,9,jb) = prm_diag(jg)%snow_gsp_rate(jc,jb)
          phdiag_par(jc,10,jb) = prm_diag(jg)%rain_con_rate(jc,jb)
          phdiag_par(jc,11,jb) = prm_diag(jg)%snow_con_rate(jc,jb)
          phdiag_par(jc,12,jb) = prm_diag(jg)%gz0(jc,jb)
          phdiag_par(jc,13,jb) = prm_diag(jg)%tcm(jc,jb)
          phdiag_par(jc,14,jb) = prm_diag(jg)%tch(jc,jb)
          phdiag_par(jc,15,jb) = prm_diag(jg)%tfm(jc,jb)
          phdiag_par(jc,16,jb) = prm_diag(jg)%tfh(jc,jb)
          phdiag_par(jc,17,jb) = prm_diag(jg)%tfv(jc,jb)
          phdiag_par(jc,18,jb) = prm_diag(jg)%t_2m(jc,jb)
          phdiag_par(jc,19,jb) = prm_diag(jg)%qv_2m(jc,jb)
          phdiag_par(jc,20,jb) = prm_diag(jg)%td_2m(jc,jb)
          phdiag_par(jc,21,jb) = prm_diag(jg)%rh_2m(jc,jb)
          phdiag_par(jc,22,jb) = prm_diag(jg)%u_10m(jc,jb)
          phdiag_par(jc,23,jb) = prm_diag(jg)%v_10m(jc,jb)
          IF (atm_phy_nwp_config(jg)%lhave_graupel) THEN
            phdiag_par(jc,24,jb) = prm_diag(jg)%graupel_gsp(jc,jb)
            phdiag_par(jc,25,jb) = prm_diag(jg)%graupel_gsp_rate(jc,jb)
          ELSE
            phdiag_par(jc,24,jb) = 0._wp
            phdiag_par(jc,25,jb) = 0._wp
          END IF
          IF (atm_phy_nwp_config(jg)%l2moment) THEN
            phdiag_par(jc,26,jb) = prm_diag(jg)%hail_gsp(jc,jb)
            phdiag_par(jc,27,jb) = prm_diag(jg)%hail_gsp_rate(jc,jb)
          ELSE
            phdiag_par(jc,26,jb) = 0._wp
            phdiag_par(jc,27,jb) = 0._wp
          ENDIF
          IF (ANY((/1,2,4,5,6,7,8/) == atm_phy_nwp_config(jg)%inwp_gscp)) THEN
            phdiag_par(jc,28,jb) = prm_diag(jg)%ice_gsp(jc,jb)
            phdiag_par(jc,29,jb) = prm_diag(jg)%ice_gsp_rate(jc,jb)
          ELSE
            phdiag_par(jc,28,jb) = 0._wp
            phdiag_par(jc,29,jb) = 0._wp
          END IF
          phdiag_par(jc,30,jb) = prm_diag(jg)%tot_prec_d(jc,jb)
          phdiag_par(jc,31,jb) = prm_diag(jg)%prec_gsp_d(jc,jb)
          phdiag_par(jc,32,jb) = prm_diag(jg)%prec_con_d(jc,jb)
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
        ! Collect soil variables
        DO jk = 1, nlev_soil
          jk1 = jk + nlev_soil
          DO jc = i_startidx, i_endidx
            lndvars_par(jc,jk,jb)  = p_parent_ldiag%w_so(jc,jk,jb)
            lndvars_par(jc,jk1,jb) = p_parent_ldiag%t_so(jc,jk,jb) - tsfc_ref_p(jc,jb)
          ENDDO
        ENDDO

        jk1 = 2*nlev_soil + 1
        DO jc = i_startidx, i_endidx
          lndvars_par(jc,jk1,jb)   = p_parent_ldiag%t_so(jc,nlev_soil+1,jb)-tsfc_ref_p(jc,jb)
          lndvars_par(jc,jk1+1,jb) = p_parent_lprog%t_g(jc,jb)    - tsfc_ref_p(jc,jb)
          IF (p_parent_ldiag%t_s(jc,jb) > 10._wp) THEN
            lndvars_par(jc,jk1+2,jb) = p_parent_ldiag%t_s(jc,jb)  - tsfc_ref_p(jc,jb)
          ELSE ! t_s has missing values over water - use t_g instead
            lndvars_par(jc,jk1+2,jb) = lndvars_par(jc,jk1+1,jb)
          ENDIF
          lndvars_par(jc,jk1+3,jb) = p_parent_ldiag%t_snow(jc,jb) - tsfc_ref_p(jc,jb)
          lndvars_par(jc,jk1+4,jb) = p_parent_ldiag%w_snow(jc,jb)
          lndvars_par(jc,jk1+5,jb) = p_parent_ldiag%rho_snow(jc,jb)
          lndvars_par(jc,jk1+6,jb) = p_parent_ldiag%w_i(jc,jb)
          lndvars_par(jc,jk1+7,jb) = p_parent_ldiag%freshsnow(jc,jb)
          lndvars_par(jc,jk1+8,jb) = MERGE(MAX(271._wp,p_parent_ldiag%t_so(jc,4,jb)), & ! fill t_seasfc with t_so where undefined
                                     p_parent_ldiag%t_seasfc(jc,jb),p_parent_ldiag%t_seasfc(jc,jb)<=0._wp)
          lndvars_par(jc,jk1+9,jb) = p_parent_ldiag%qv_s(jc,jb)
          IF (itype_trvg == 3) THEN
            lndvars_par(jc,jk1+10,jb) = p_parent_ldiag%plantevap(jc,jb)
          ELSE
            lndvars_par(jc,jk1+10,jb) = 0._wp
          ENDIF
          IF (itype_snowevap == 3) THEN
            lndvars_par(jc,jk1+11,jb) = p_parent_ldiag%hsnow_max(jc,jb)
            lndvars_par(jc,jk1+12,jb) = p_parent_ldiag%snow_age(jc,jb)
          ELSE
            lndvars_par(jc,jk1+11,jb) = 0._wp
            lndvars_par(jc,jk1+12,jb) = 0._wp
          ENDIF
          IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) THEN
            lndvars_par(jc,jk1+13,jb) = p_parent_diag%t2m_bias(jc,jb)
          ELSE IF (icpl_da_sfcevap >= 3) THEN
            lndvars_par(jc,jk1+13,jb) = p_parent_diag%t_avginc(jc,jb)
          ELSE
            lndvars_par(jc,jk1+13,jb) = 0._wp
          ENDIF
          lndvars_par(jc,jk1+14,jb) = p_parent_ldiag%t_sk(jc,jb)
          IF (icpl_da_sfcevap >= 2) THEN
            lndvars_par(jc,jk1+15,jb) = p_parent_diag%rh_avginc(jc,jb)
          ELSE
            lndvars_par(jc,jk1+15,jb) = 0._wp
          ENDIF
          IF (icpl_da_skinc >= 1) THEN
            lndvars_par(jc,jk1+16,jb) = p_parent_diag%t_wgt_avginc(jc,jb)
          ELSE
            lndvars_par(jc,jk1+16,jb) = 0._wp
          ENDIF
          IF (icpl_da_sfcfric >= 1) THEN
            lndvars_par(jc,jk1+17,jb) = p_parent_diag%vabs_avginc(jc,jb)
          ELSE
            lndvars_par(jc,jk1+17,jb) = 0._wp
          ENDIF
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lseaice) THEN
        DO jc = i_startidx, i_endidx
          IF (p_parent_wprog%t_ice(jc,jb) > 10._wp) THEN
            wtrvars_par(jc,1,jb) = p_parent_wprog%t_ice(jc,jb)
          ELSE
            wtrvars_par(jc,1,jb) = p_parent_lprog%t_g(jc,jb)
          ENDIF
          wtrvars_par(jc,2,jb) = p_parent_wprog%h_ice(jc,jb)
          wtrvars_par(jc,3,jb) = p_parent_wprog%t_snow_si(jc,jb)
          wtrvars_par(jc,4,jb) = p_parent_wprog%h_snow_si(jc,jb)
          wtrvars_par(jc,5,jb) = p_parent_ldiag%fr_seaice(jc,jb)
          IF (lprog_albsi) THEN
            wtrvars_par(jc,6,jb) = MAX(csalb(ist_seawtr),p_parent_wprog%alb_si(jc,jb))
          ELSE
            wtrvars_par(jc,6,jb) = ALB_SI_MISSVAL ! -1
          ENDIF
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Convert wsoil into SMI for interpolation
    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) &
      CALL wsoil_to_smi(p_patch(jg), ext_data(jg)%atm%list_land, &
        &               ext_data(jg)%atm%soiltyp, lndvars_par(:,1:nlev_soil,:))

    ! Step 1b: execute boundary interpolation

    CALL interpol2_vec_grf (p_pp=p_patch(jg), p_pc=p_pc, p_grf=p_grf_state(jg)%p_dom(i_chidx), nfields=1, lacc=.FALSE., &
      f3din1=p_parent_prog%vn, f3dout1=p_child_prog%vn)

    CALL interpol_scal_grf (p_pp=p_patch(jg), p_pc=p_pc, p_grf=p_grf_state(jg)%p_dom(i_chidx), nfields=3, lacc=.FALSE.,&
      f3din1=rho_pr_par,      f3dout1=p_child_prog%rho,     &
      f3din2=thv_pr_par,      f3dout2=p_child_prog%theta_v, &
      f3din3=p_parent_prog%w, f3dout3=p_child_prog%w        )

    IF (ltransport) THEN
      l_limit(:) = .TRUE. ! apply positive definite limiter on tracers

      CALL interpol_scal_grf ( p_pp=p_patch(jg), p_pc=p_pc, p_grf=p_grf_state(jg)%p_dom(i_chidx), nfields=ntracer, lacc=.FALSE.,  &
        f4din1=p_parent_prog_rcf%tracer, f4dout1=p_child_prog_rcf%tracer, llimit_nneg=l_limit)
    ENDIF

    IF (ltransport .AND. iprog_aero >= 1) THEN
      CALL interpol_scal_grf ( p_pp=p_patch(jg), p_pc=p_pc, p_grf=p_grf_state(jg)%p_dom(i_chidx), nfields=1, lacc=.FALSE.,  &
        f3din1=prm_diag(jg)%aerosol, f3dout1=prm_diag(jgc)%aerosol, llimit_nneg=(/.TRUE./), lnoshift=.TRUE.)
    ENDIF

    IF (iforcing == inwp) THEN
      CALL sync_patch_array(SYNC_C,p_patch(jg),phdiag_par)
      CALL interpol_scal_grf (p_pp=p_patch(jg), p_pc=p_pc, p_grf=p_grf_state(jg)%p_dom(i_chidx), nfields=1, lacc=.FALSE., &
        f3din1=phdiag_par, f3dout1=phdiag_chi, lnoshift=.TRUE.                 )
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
      CALL sync_patch_array(SYNC_C,p_patch(jg),lndvars_par)
      CALL interpol_scal_grf (p_pp=p_patch(jg), p_pc=p_pc, p_grf=p_grf_state(jg)%p_dom(i_chidx), nfields=1, lacc=.FALSE., &
        f3din1=lndvars_par, f3dout1=lndvars_chi, lnoshift=.TRUE.                 )
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lseaice) THEN
      CALL sync_patch_array(SYNC_C,p_patch(jg),wtrvars_par)
      CALL interpol_scal_grf (p_pp=p_patch(jg), p_pc=p_pc, p_grf=p_grf_state(jg)%p_dom(i_chidx), nfields=1, lacc=.FALSE., &
        f3din1=wtrvars_par, f3dout1=wtrvars_chi, lnoshift=.TRUE.                 )
    ENDIF

    ! Step 2: Interpolation of fields in the model interior

    ! Step 2a: Copy prognostic variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data_mult(p_pat=p_pp%comm_pat_glb_to_loc_c,          &
      &                     lacc=.FALSE.,                              &
      &                     nfields=3,       ndim2tot=3*nlev_p+1,      &
      &                     RECV1=rho_pr_lp, SEND1=rho_pr_par,         &
      &                     RECV2=thv_pr_lp, SEND2=thv_pr_par,         &
      &                     RECV3=w_lp,      SEND3=p_parent_prog%w     )

    CALL exchange_data(p_pat=p_pp%comm_pat_glb_to_loc_e, lacc=.FALSE., RECV=vn_lp, SEND=p_parent_prog%vn)

    IF (ltransport) THEN
      CALL exchange_data_mult(p_pat=p_pp%comm_pat_glb_to_loc_c,                    &
        &                     lacc=.FALSE.,                                        &
        &                     nfields=ntracer,  ndim2tot=ntracer*nlev_p,           &
        &                     RECV4D=tracer_lp, SEND4D=p_parent_prog_rcf%tracer    )
    ENDIF

    IF (ltransport .AND. iprog_aero >= 1) THEN
      CALL exchange_data(p_pat=p_pp%comm_pat_glb_to_loc_c, lacc=.FALSE., RECV=aero_lp, SEND=prm_diag(jg)%aerosol)
    ENDIF

    IF (iforcing == inwp) &
      &      CALL exchange_data(p_pat=p_pp%comm_pat_glb_to_loc_c, lacc=.FALSE., RECV=phdiag_lp, SEND=phdiag_par)

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) &
      &      CALL exchange_data(p_pat=p_pp%comm_pat_glb_to_loc_c, lacc=.FALSE., RECV=lndvars_lp, SEND=lndvars_par)

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lseaice) &
      &      CALL exchange_data(p_pat=p_pp%comm_pat_glb_to_loc_c, lacc=.FALSE., RECV=wtrvars_lp, SEND=wtrvars_par)

    ! Step 2b: Perform interpolation from local parent to child grid

    ! Note: the sync routines cannot be used for the synchronization on the local
    ! parent grid

    IF(l_parallel) CALL exchange_data(p_pat=p_pp%comm_pat_e, lacc=.FALSE., recv=vn_lp)
    CALL interpol_vec_nudging (p_pp, p_pc, p_int, p_grf%p_dom(i_chidx), &
                               nshift, 1, vn_lp, p_child_prog%vn, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E,p_pc,p_child_prog%vn)

    IF(l_parallel) CALL exchange_data_mult(p_pat=p_pp%comm_pat_c,                         &
                               lacc=.FALSE.,                                              &
                               nfields=3,  ndim2tot=3*nlev_p+1,                           &
                               recv1=thv_pr_lp, recv2=rho_pr_lp, recv3=w_lp)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), nshift, 3, 1,          &
                                f3din1=thv_pr_lp, f3dout1=p_child_prog%theta_v,           &
                                f3din2=rho_pr_lp, f3dout2=p_child_prog%rho,               &
                                f3din3=w_lp,      f3dout3=p_child_prog%w                  )
    CALL sync_patch_array_mult(SYNC_C,p_pc,3,p_child_prog%theta_v,p_child_prog%rho,  &
                               p_child_prog%w)

    IF (ltransport) THEN
      IF(l_parallel) CALL exchange_data_mult(p_pat=p_pp%comm_pat_c,   &
                                             lacc=.FALSE.,            &
                                             nfields=ntracer,         &
                                             ndim2tot=ntracer*nlev_p, &
                                             recv4d=tracer_lp)
      l_limit(:) = .TRUE. ! apply positive definite limiter
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx),          &
                                  nshift, ntracer, 1, f4din=tracer_lp,        &
                                  f4dout=p_child_prog_rcf%tracer, llimit_nneg=l_limit)
      CALL sync_patch_array_mult(SYNC_C,p_pc,ntracer,f4din=p_child_prog_rcf%tracer)
    ENDIF

    IF (ltransport .AND. iprog_aero >= 1) THEN
      IF(l_parallel) CALL exchange_data(p_pat=p_pp%comm_pat_c, lacc=.FALSE., recv=aero_lp)
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx),          &
         0, 1, 1, f3din1=aero_lp, f3dout1=prm_diag(jgc)%aerosol, llimit_nneg=(/.TRUE./))
      CALL sync_patch_array(SYNC_C,p_pc,prm_diag(jgc)%aerosol)
    ENDIF

    IF (iforcing == inwp) THEN
      IF(l_parallel) CALL exchange_data(p_pat=p_pp%comm_pat_c, lacc=.FALSE., recv=phdiag_lp)
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), 0,      &
                                  1, 1, f3din1=phdiag_lp, f3dout1=phdiag_chi, overshoot_fac=1.005_wp )
      CALL sync_patch_array(SYNC_C,p_pc,phdiag_chi)
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
      IF(l_parallel) CALL exchange_data(p_pat=p_pp%comm_pat_c, lacc=.FALSE., recv=lndvars_lp)
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), 0,      &
                                  1, 1, f3din1=lndvars_lp, f3dout1=lndvars_chi, overshoot_fac=1.005_wp )
      CALL sync_patch_array(SYNC_C,p_pc,lndvars_chi)
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lseaice) THEN
      IF(l_parallel) CALL exchange_data(p_pat=p_pp%comm_pat_c, lacc=.FALSE., recv=wtrvars_lp)
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), 0,      &
                                  1, 1, f3din1=wtrvars_lp, f3dout1=wtrvars_chi, overshoot_fac=1.005_wp )
      CALL sync_patch_array(SYNC_C,p_pc,wtrvars_chi)
    ENDIF

    ! Convert SMI back to wsoil
    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) &
      CALL smi_to_wsoil(p_patch(jgc), ext_data(jgc)%atm%list_land, &
        &               ext_data(jgc)%atm%soiltyp, lndvars_chi(:,1:nlev_soil,:))

    ! Step 3: Add reference state to thermodynamic variables and copy land fields
    ! from the container arrays to the prognostic variables (for the time being,
    ! all tiles are initialized with the interpolated aggregated variables)
    ! Note that nest boundary points have to be included here to obtain a proper
    ! initialization of all grid points

    i_nchdom = MAX(1,p_pc%n_childdom)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_pc%cells%start_blk(1,1)
    i_endblk   = p_pc%cells%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt,jk1,ist) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell)

      ! rho and theta currently contain perturbation fields; thus, the reference
      ! state still needs to be added. Then, rhotheta and exner are diagnosed
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
          p_child_prog%rho(jc,jk,jb) = &
            p_child_prog%rho(jc,jk,jb) + p_child_metrics%rho_ref_mc(jc,jk,jb)

          p_child_prog%theta_v(jc,jk,jb) = &
            p_child_prog%theta_v(jc,jk,jb) + p_child_metrics%theta_ref_mc(jc,jk,jb)

          p_child_prog%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd/p0ref* &
            p_child_prog%rho(jc,jk,jb) * p_child_prog%theta_v(jc,jk,jb)))

          ! exner_pr also needs to be initialized here
          p_child_diag%exner_pr(jc,jk,jb) = p_child_prog%exner(jc,jk,jb)
        ENDDO
      ENDDO

      ! Reference temperature at lowest level (taken as proxy for the surface)
      DO jc = i_startidx, i_endidx
        tsfc_ref_c(jc,jb) = p_child_metrics%theta_ref_mc(jc,nlev_c,jb) * &
                            p_child_metrics%exner_ref_mc(jc,nlev_c,jb)
      ENDDO

      IF (iforcing == inwp) THEN
        DO jc = i_startidx, i_endidx
          prm_diag(jgc)%tot_prec(jc,jb)       = MAX(0._wp,phdiag_chi(jc,1,jb))
          prm_diag(jgc)%prec_gsp(jc,jb)       = MAX(0._wp,phdiag_chi(jc,2,jb))
          prm_diag(jgc)%prec_con(jc,jb)       = MAX(0._wp,phdiag_chi(jc,3,jb))
          prm_diag(jgc)%rain_gsp(jc,jb)       = MAX(0._wp,phdiag_chi(jc,4,jb))
          prm_diag(jgc)%snow_gsp(jc,jb)       = MAX(0._wp,phdiag_chi(jc,5,jb))
          prm_diag(jgc)%rain_con(jc,jb)       = MAX(0._wp,phdiag_chi(jc,6,jb))
          prm_diag(jgc)%snow_con(jc,jb)       = MAX(0._wp,phdiag_chi(jc,7,jb))
          prm_diag(jgc)%rain_gsp_rate(jc,jb)  = phdiag_chi(jc,8,jb)
          prm_diag(jgc)%snow_gsp_rate(jc,jb)  = phdiag_chi(jc,9,jb)
          IF (atm_phy_nwp_config(jgc)%inwp_convection == 1) THEN
            prm_diag(jgc)%rain_con_rate(jc,jb)  = phdiag_chi(jc,10,jb)
            prm_diag(jgc)%snow_con_rate(jc,jb)  = phdiag_chi(jc,11,jb)
          ENDIF
          prm_diag(jgc)%gz0(jc,jb)            = phdiag_chi(jc,12,jb)
          prm_diag(jgc)%tcm(jc,jb)            = phdiag_chi(jc,13,jb)
          prm_diag(jgc)%tch(jc,jb)            = phdiag_chi(jc,14,jb)
          prm_diag(jgc)%tfm(jc,jb)            = phdiag_chi(jc,15,jb)
          prm_diag(jgc)%tfh(jc,jb)            = phdiag_chi(jc,16,jb)
          prm_diag(jgc)%tfv(jc,jb)            = phdiag_chi(jc,17,jb)
          prm_diag(jgc)%t_2m(jc,jb)           = phdiag_chi(jc,18,jb)
          prm_diag(jgc)%qv_2m(jc,jb)          = phdiag_chi(jc,19,jb)
          prm_diag(jgc)%td_2m(jc,jb)          = phdiag_chi(jc,20,jb)
          prm_diag(jgc)%rh_2m(jc,jb)          = phdiag_chi(jc,21,jb)
          prm_diag(jgc)%u_10m(jc,jb)          = phdiag_chi(jc,22,jb)
          prm_diag(jgc)%v_10m(jc,jb)          = phdiag_chi(jc,23,jb)
          IF (atm_phy_nwp_config(jgc)%lhave_graupel) THEN
            prm_diag(jgc)%graupel_gsp(jc,jb)      = MAX(0._wp,phdiag_chi(jc,24,jb))
            prm_diag(jgc)%graupel_gsp_rate(jc,jb) = phdiag_chi(jc,25,jb)
          END IF
          IF (atm_phy_nwp_config(jgc)%l2moment) THEN
            prm_diag(jgc)%hail_gsp(jc,jb)      = MAX(0._wp,phdiag_chi(jc,26,jb))
            prm_diag(jgc)%hail_gsp_rate(jc,jb) = phdiag_chi(jc,27,jb) 
          ENDIF
          IF (ANY((/1,2,4,5,6,7,8/) == atm_phy_nwp_config(jgc)%inwp_gscp)) THEN
            prm_diag(jgc)%ice_gsp(jc,jb)        = MAX(0._wp,phdiag_chi(jc,28,jb))
            prm_diag(jgc)%ice_gsp_rate(jc,jb)   = phdiag_chi(jc,29,jb)
          END IF
          prm_diag(jgc)%tot_prec_d(jc,jb)      = MAX(0._wp,phdiag_chi(jc,30,jb))
          prm_diag(jgc)%prec_gsp_d(jc,jb)      = MAX(0._wp,phdiag_chi(jc,31,jb))
          prm_diag(jgc)%prec_con_d(jc,jb)      = MAX(0._wp,phdiag_chi(jc,32,jb))
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jgc)%inwp_surface == 1) THEN
        ! Distribute soil variables
        DO jt = 1, ntiles_total
          DO jk = 1, nlev_soil
            jk1 = jk + nlev_soil
            DO jc = i_startidx, i_endidx
              p_child_lprog%w_so_t(jc,jk,jb,jt)  = MAX(0._wp,lndvars_chi(jc,jk,jb))
              p_child_lprog%t_so_t(jc,jk,jb,jt)  = lndvars_chi(jc,jk1,jb) + tsfc_ref_c(jc,jb)
              ! limit w_so_t to pore volume and dryness point - TERRA crashes otherwise
              ist = ext_data(jgc)%atm%soiltyp(jc,jb)
              IF (ist >= 3 .AND. ist <= 8) THEN ! soil types with non-zero water content
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*cadp(ist), &
                  MIN(dzsoil(jk)*cporv(ist),p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ENDIF
              p_child_lprog2%t_so_t(jc,jk,jb,jt) = p_child_lprog%t_so_t(jc,jk,jb,jt)
              p_child_lprog2%w_so_t(jc,jk,jb,jt) = p_child_lprog%w_so_t(jc,jk,jb,jt)
            ENDDO
          ENDDO
        ENDDO

        jk1 = 2*nlev_soil + 1
        DO jt = 1, ntiles_total
          DO jc = i_startidx, i_endidx
            p_child_lprog%t_so_t(jc,nlev_soil+1,jb,jt) = lndvars_chi(jc,jk1,jb) + &
              tsfc_ref_c(jc,jb)
            ! here we need to initialize t_g and t_g_t because t_g is copied to t_g_t
            ! in nwp_surface_init for land points
            p_child_lprog%t_g(jc,jb) = lndvars_chi(jc,jk1+1,jb) + tsfc_ref_c(jc,jb)
            p_child_lprog2%t_g(jc,jb) = p_child_lprog%t_g(jc,jb)
            p_child_lprog%t_g_t(jc,jb,jt) = p_child_lprog%t_g(jc,jb)
            p_child_lprog2%t_g_t(jc,jb,jt) = p_child_lprog%t_g(jc,jb)
            p_child_lprog%t_s_t(jc,jb,jt) = lndvars_chi(jc,jk1+2,jb) + tsfc_ref_c(jc,jb)
            p_child_lprog2%t_s_t(jc,jb,jt) = p_child_lprog%t_s_t(jc,jb,jt)
            p_child_lprog%t_snow_t(jc,jb,jt) = lndvars_chi(jc,jk1+3,jb) + tsfc_ref_c(jc,jb)
            p_child_lprog%w_snow_t(jc,jb,jt) = MAX(0._wp,lndvars_chi(jc,jk1+4,jb))
            p_child_lprog%rho_snow_t(jc,jb,jt) = lndvars_chi(jc,jk1+5,jb)
            IF (p_child_lprog%rho_snow_t(jc,jb,jt) < 0.75_wp*crhosminf .OR. &
                p_child_lprog%w_snow_t(jc,jb,jt) < 1.e-6_wp) p_child_lprog%rho_snow_t(jc,jb,jt) = 250._wp
            p_child_lprog%w_i_t(jc,jb,jt) = MAX(0._wp,lndvars_chi(jc,jk1+6,jb))
            p_child_ldiag%freshsnow_t(jc,jb,jt) = MAX(0._wp,MIN(1._wp,lndvars_chi(jc,jk1+7,jb)))
            p_child_ldiag%t_seasfc(jc,jb) = lndvars_chi(jc,jk1+8,jb)
            p_child_ldiag%qv_s(jc,jb)     = lndvars_chi(jc,jk1+9,jb)
            IF (itype_trvg == 3) THEN
              p_child_ldiag%plantevap(jc,jb)      = lndvars_chi(jc,jk1+10,jb)
              p_child_ldiag%plantevap_t(jc,jb,jt) = p_child_ldiag%plantevap(jc,jb)
            ENDIF
            IF (itype_snowevap == 3) THEN
              p_child_ldiag%hsnow_max(jc,jb) = lndvars_chi(jc,jk1+11,jb)
              p_child_ldiag%snow_age(jc,jb)  = lndvars_chi(jc,jk1+12,jb)
            ENDIF
            IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) THEN
              p_child_diag%t2m_bias(jc,jb) = lndvars_chi(jc,jk1+13,jb)
            ELSE IF (icpl_da_sfcevap >= 3) THEN
              p_child_diag%t_avginc(jc,jb) = lndvars_chi(jc,jk1+13,jb)
            ENDIF
            p_child_ldiag%t_sk(jc,jb)      = lndvars_chi(jc,jk1+14,jb)
            p_child_lprog%t_sk_t(jc,jb,jt) = p_child_ldiag%t_sk(jc,jb)
            p_child_lprog2%t_sk_t(jc,jb,jt) = p_child_ldiag%t_sk(jc,jb)
            IF (icpl_da_sfcevap >= 2) THEN
              p_child_diag%rh_avginc(jc,jb) = lndvars_chi(jc,jk1+15,jb)
            ENDIF
            IF (icpl_da_skinc >= 1) THEN
              p_child_diag%t_wgt_avginc(jc,jb) = lndvars_chi(jc,jk1+16,jb)
            ENDIF
            IF (icpl_da_sfcfric >= 1) THEN
              p_child_diag%vabs_avginc(jc,jb) = lndvars_chi(jc,jk1+17,jb)
            ENDIF
          ENDDO
        ENDDO
        DO jt = ntiles_total+1, ntiles_total+ntiles_water
          DO jc = i_startidx, i_endidx
            p_child_lprog%t_g_t(jc,jb,jt)  = p_child_lprog%t_g(jc,jb)
            p_child_lprog2%t_g_t(jc,jb,jt) = p_child_lprog%t_g(jc,jb)
            p_child_lprog%t_s_t(jc,jb,jt)  = p_child_ldiag%t_s(jc,jb)
            p_child_lprog2%t_s_t(jc,jb,jt) = p_child_ldiag%t_s(jc,jb)
            p_child_ldiag%qv_s_t(jc,jb,jt) = p_child_ldiag%qv_s(jc,jb)
          ENDDO
        ENDDO

      ENDIF


      IF (atm_phy_nwp_config(jgc)%inwp_surface == 1 .AND. lseaice) THEN
        DO jc = i_startidx, i_endidx
          p_child_wprog%t_ice(jc,jb)     = MIN(tmelt,wtrvars_chi(jc,1,jb))
          p_child_wprog%h_ice(jc,jb)     = MAX(0._wp,wtrvars_chi(jc,2,jb))
          p_child_wprog%t_snow_si(jc,jb) = MIN(tmelt,wtrvars_chi(jc,3,jb))
          p_child_wprog%h_snow_si(jc,jb) = MAX(0._wp,wtrvars_chi(jc,4,jb))
          p_child_ldiag%fr_seaice(jc,jb) = MAX(0._wp,MIN(1._wp,wtrvars_chi(jc,5,jb)))
          p_child_wprog%alb_si(jc,jb)    = wtrvars_chi(jc,6,jb)
          IF (p_child_ldiag%fr_seaice(jc,jb) < frsi_min )         p_child_ldiag%fr_seaice(jc,jb) = 0._wp
          IF (p_child_ldiag%fr_seaice(jc,jb) > (1._wp-frsi_min) ) p_child_ldiag%fr_seaice(jc,jb) = 1._wp
          IF (lprog_albsi) THEN
            IF (p_child_ldiag%fr_seaice(jc,jb) == 0._wp) THEN
              p_child_wprog%alb_si(jc,jb) = ALB_SI_MISSVAL ! -1
            ELSE
              p_child_wprog%alb_si(jc,jb) = MIN(0.8_wp, MAX(0.4_wp, p_child_wprog%alb_si(jc,jb)))
            ENDIF
          ENDIF
          IF (ext_data(jgc)%atm%fr_land(jc,jb) >= 1._wp-MAX(frlake_thrhld,frsea_thrhld)) THEN ! pure land point
            p_child_wprog%h_ice(jc,jb) = 0._wp
            p_child_wprog%h_snow_si(jc,jb) = 0._wp
            p_child_ldiag%fr_seaice(jc,jb) = 0._wp
          ENDIF
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jgc)%inwp_surface == 1 .AND. llake) THEN

        CALL flake_coldinit(                                        &
          &   nflkgb      = ext_data(jgc)%atm%list_lake%ncount(jb), &  ! in
          &   idx_lst_fp  = ext_data(jgc)%atm%list_lake%idx (:,jb), &  ! in
          &   depth_lk    = ext_data(jgc)%atm%depth_lk  (:,jb),     &  ! in
          ! here, a proper estimate of the lake surface temperature is required
          &   tskin       = lndvars_chi(:,num_lndvars,jb),          &  ! in
          &   t_snow_lk_p = p_child_wprog%t_snow_lk(:,jb),          &
          &   h_snow_lk_p = p_child_wprog%h_snow_lk(:,jb),          &
          &   t_ice_p     = p_child_wprog%t_ice    (:,jb),          &
          &   h_ice_p     = p_child_wprog%h_ice    (:,jb),          &
          &   t_mnw_lk_p  = p_child_wprog%t_mnw_lk (:,jb),          &
          &   t_wml_lk_p  = p_child_wprog%t_wml_lk (:,jb),          &
          &   t_bot_lk_p  = p_child_wprog%t_bot_lk (:,jb),          &
          &   c_t_lk_p    = p_child_wprog%c_t_lk   (:,jb),          &
          &   h_ml_lk_p   = p_child_wprog%h_ml_lk  (:,jb),          &
          &   t_b1_lk_p   = p_child_wprog%t_b1_lk  (:,jb),          &
          &   h_b1_lk_p   = p_child_wprog%h_b1_lk  (:,jb),          &
          &   t_g_lk_p    = p_child_lprog%t_g_t    (:,jb,isub_lake) )

      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL rbf_vec_interpol_cell(p_child_prog%vn, p_patch(jgc), p_int_state(jgc),&
                               p_child_diag%u, p_child_diag%v)

    CALL diagnose_pres_temp(p_child_metrics, p_child_prog, p_child_prog_rcf, p_child_diag, &
                            p_patch(jgc), opt_calc_temp=.TRUE., opt_calc_pres=.TRUE.,      &
                            lnd_prog=p_child_lprog)


    DEALLOCATE(thv_pr_par, rho_pr_par, lndvars_par, wtrvars_par, phdiag_par, lndvars_chi, &
      &       wtrvars_chi, phdiag_chi, tsfc_ref_p, tsfc_ref_c, vn_lp, w_lp, thv_pr_lp   , &
      &       rho_pr_lp, phdiag_lp, tracer_lp, lndvars_lp, wtrvars_lp, aero_lp)

  END SUBROUTINE initialize_nest


  !-------------
  !>
  !! SUBROUTINEs interpolate_scal|vn_increments
  !!
  !! Driver routines for interpolating data assimilation increments from a parent domain to a child domain
  !! Scalars and wind increments are processed separately because wind increments are filtered in
  !! SR create_dwdanainc_atm, and the parent-to-child interpolation is executed after filtering
  !!
  SUBROUTINE interpolate_scal_increments(initicon, jg, jgc)

    TYPE(t_initicon_state), TARGET, INTENT(INOUT) :: initicon(:)

    INTEGER, INTENT(IN) :: jg   ! parent (source) domain ID
    INTEGER, INTENT(IN) :: jgc  ! child  (target) domain ID

    ! local pointers

    TYPE(t_patch),         POINTER  :: p_pp, p_pc
    TYPE(t_gridref_state), POINTER  :: p_grf
    TYPE(t_int_state),     POINTER  :: p_int

    ! Local arrays for variables living on the local parent grid in the MPI case.
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: temp_lp, pres_lp, qv_lp

    INTEGER :: nlev_p, nshift, i_chidx

    LOGICAL :: l_parallel


    l_parallel = my_process_is_mpi_parallel()

    p_pp           => p_patch_local_parent(jgc)
    p_pc           => p_patch(jgc)
    p_grf          => p_grf_state_local_parent(jgc)
    p_int          => p_int_state_local_parent(jgc)

    i_chidx  = p_pc%parent_child_index

    ! number of full levels of parent domain
    nlev_p   = p_patch(jg)%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift

    IF (msg_level >= 7) THEN
      WRITE(message_text,'(a,i2,a,i2)') 'Interpolation of DA increments for scalar fields, domain ',jg,' =>',jgc
      CALL message('interpolate_increments',message_text)
    ENDIF

    ALLOCATE(temp_lp (nproma, nlev_p, p_pp%nblks_c),  &
             pres_lp (nproma, nlev_p, p_pp%nblks_c),  &
             qv_lp   (nproma, nlev_p, p_pp%nblks_c)   )


    ! Parent-to-child interpolation of atmospheric increments
    !
    ! Step 1: boundary interpolation
    !

    CALL interpol_scal_grf (p_pp=p_patch(jg), p_pc=p_pc, p_grf=p_grf_state(jg)%p_dom(i_chidx), nfields=3, lacc=.FALSE., &
      f3din1=initicon(jg)%atm_inc%temp, f3dout1=initicon(jgc)%atm_inc%temp, &
      f3din2=initicon(jg)%atm_inc%pres, f3dout2=initicon(jgc)%atm_inc%pres, &
      f3din3=initicon(jg)%atm_inc%qv,   f3dout3=initicon(jgc)%atm_inc%qv    )

    ! Step 2: Interpolation of fields in the model interior

    ! Step 2a: Copy prognostic variables from parent grid to fields on local parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)


    CALL exchange_data_mult(p_pat=p_pp%comm_pat_glb_to_loc_c,                &
      &                     lacc=.FALSE.,                                    &
      &                     nfields=3, ndim2tot=3*nlev_p,                    &
      &                     RECV1=temp_lp, SEND1=initicon(jg)%atm_inc%temp,  &
      &                     RECV2=pres_lp, SEND2=initicon(jg)%atm_inc%pres,  &
      &                     RECV3=qv_lp,   SEND3=initicon(jg)%atm_inc%qv     )

    ! Step 2b: Synchronize variables on local parent grids. Note that the sync routines cannot be used in this case
    IF (l_parallel) THEN
      CALL exchange_data_mult(p_pat=p_pp%comm_pat_c,                      &
                              lacc=.FALSE., nfields=3, ndim2tot=3*nlev_p, &
                              recv1=temp_lp, recv2=pres_lp, recv3=qv_lp)
    ENDIF

    ! Step 2c: Perform interpolation from local parent to child grid

    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), nshift, 3, 1,          &
                                f3din1=temp_lp, f3dout1=initicon(jgc)%atm_inc%temp,       &
                                f3din2=pres_lp, f3dout2=initicon(jgc)%atm_inc%pres,       &
                                f3din3=qv_lp,   f3dout3=initicon(jgc)%atm_inc%qv          )
    CALL sync_patch_array_mult(SYNC_C,p_pc,3,initicon(jgc)%atm_inc%temp,initicon(jgc)%atm_inc%pres,  &
                               initicon(jgc)%atm_inc%qv)

    DEALLOCATE(temp_lp, pres_lp, qv_lp)

  END SUBROUTINE interpolate_scal_increments

  SUBROUTINE interpolate_vn_increments(initicon, jg, jgc)

    TYPE(t_initicon_state), TARGET, INTENT(INOUT) :: initicon(:)

    INTEGER, INTENT(IN) :: jg   ! parent (source) domain ID
    INTEGER, INTENT(IN) :: jgc  ! child  (target) domain ID

    ! local pointers

    TYPE(t_patch),         POINTER  :: p_pp, p_pc
    TYPE(t_gridref_state), POINTER  :: p_grf
    TYPE(t_int_state),     POINTER  :: p_int

    ! Local arrays for variables living on the local parent grid in the MPI case.
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: vn_lp

    INTEGER :: nlev_p, nshift, i_chidx

    LOGICAL :: l_parallel


    IF (.NOT. my_process_is_mpi_parallel()) THEN
      l_parallel = .FALSE.
    ELSE
      l_parallel = .TRUE.
    ENDIF


    p_pp           => p_patch_local_parent(jgc)
    p_pc           => p_patch(jgc)
    p_grf          => p_grf_state_local_parent(jgc)
    p_int          => p_int_state_local_parent(jgc)

    i_chidx  = p_pc%parent_child_index

    ! number of full levels of parent domain
    nlev_p   = p_patch(jg)%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift

    IF (msg_level >= 7) THEN
      WRITE(message_text,'(a,i2,a,i2)') 'Interpolation of DA increments for wind, domain ',jg,' =>',jgc
      CALL message('interpolate_increments',message_text)
    ENDIF

    ALLOCATE(vn_lp(nproma, nlev_p, p_pp%nblks_e))


    ! Parent-to-child interpolation of atmospheric increments
    !
    ! Step 1: boundary interpolation
    !
    CALL interpol2_vec_grf (p_pp=p_patch(jg), p_pc=p_pc, p_grf=p_grf_state(jg)%p_dom(i_chidx), nfields=1, lacc=.FALSE., &
      f3din1=initicon(jg)%atm_inc%vn, f3dout1=initicon(jgc)%atm_inc%vn)

    ! Step 2: Interpolation of fields in the model interior

    ! Step 2a: Copy prognostic variables from parent grid to fields on local parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data(p_pat=p_pp%comm_pat_glb_to_loc_e, lacc=.FALSE., RECV=vn_lp, SEND=initicon(jg)%atm_inc%vn)


    ! Step 2b: Synchronize variables on local parent grids. Note that the sync routines cannot be used in this case
    IF (l_parallel) THEN
      CALL exchange_data(p_pat=p_pp%comm_pat_e, lacc=.FALSE., recv=vn_lp)
    ENDIF

    ! Step 2c: Perform interpolation from local parent to child grid
    CALL interpol_vec_nudging (p_pp, p_pc, p_int, p_grf%p_dom(i_chidx),   &
                               nshift, 1, vn_lp, initicon(jgc)%atm_inc%vn, lacc=.FALSE.)
    CALL sync_patch_array(SYNC_E,p_pc,initicon(jgc)%atm_inc%vn)

    DEALLOCATE(vn_lp)

  END SUBROUTINE interpolate_vn_increments


  !-------------
  !>
  !! SUBROUTINE interpolate_sfcana
  !!
  !! Driver routine for interpolating surface analysis data from a parent domain to a child domain
  !! The routine is supposed to work in combination with incremental analysis update only (MODE_IAU);
  !! it processes 
  !! * w_so      (increments) 
  !! * freshsnow (increments)
  !! * h_snow    (increments)
  !! * t_2m      (increments)
  !! * fr_ice    (full field)
  !! * sst       (full field)
  !!   i.e. t_so(0) over sea points only or t_seasfc
  !!
  SUBROUTINE interpolate_sfcana(initicon, inputInstructions, jg, jgc, p_lnd_state )

    TYPE(t_initicon_state),         INTENT(INOUT) :: initicon(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    INTEGER, INTENT(IN) :: jg   ! parent (source) domain ID
    INTEGER, INTENT(IN) :: jgc  ! child  (target) domain ID
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)

    ! local pointers
    TYPE(t_patch),         POINTER  :: p_pp, p_pc
    TYPE(t_gridref_state), POINTER  :: p_grf
    TYPE(t_int_state),     POINTER  :: p_int
    TYPE(t_lnd_diag),      POINTER  :: p_parent_ldiag
    TYPE(t_lnd_diag),      POINTER  :: p_child_ldiag

    ! local variables

    ! Indices
    INTEGER :: jb, jc, jk, jk1, i_chidx, i_startblk, i_endblk, &
               i_startidx, i_endidx, var_src

    INTEGER :: num_lndvars

    ! Local arrays for variables living on the local parent grid in the MPI case. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: lndvars_lp

    ! Local arrays on the parent or child grid. These would not have to be allocatable,
    ! but the computational overhead does not matter for an initialization routine
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: lndvars_par, lndvars_chi

    LOGICAL :: l_parallel

    !-----------------------------------------------------------------------


    IF (msg_level >= 7) THEN
      WRITE(message_text,'(a,i2,a,i2)') 'Interpolation of surface analysis data, domain ',jg,' =>',jgc
      CALL message('interpolate_sfcana',message_text)
    ENDIF

    l_parallel = my_process_is_mpi_parallel()

    p_parent_ldiag    => p_lnd_state(jg)%diag_lnd
    p_child_ldiag     => p_lnd_state(jgc)%diag_lnd

    p_pc              => p_patch(jgc)

    p_grf => p_grf_state_local_parent(jgc)
    p_int => p_int_state_local_parent(jgc)
    p_pp  => p_patch_local_parent(jgc)

    i_chidx  = p_pc%parent_child_index

    ! Number of variables to be interpolated
    num_lndvars = nlev_soil + 5

    ALLOCATE(lndvars_par (nproma, num_lndvars, p_patch(jg)%nblks_c), &
             lndvars_chi (nproma, num_lndvars, p_patch(jgc)%nblks_c) )

    ALLOCATE(lndvars_lp (nproma, num_lndvars, p_pp%nblks_c) )

    IF (p_test_run) THEN
      lndvars_par  = 0._wp
      lndvars_chi  = 0._wp
    ENDIF


    ! Step 1: boundary interpolation
    !
    ! Note: the sepration between boundary interpolation and the interpolation
    ! of the fields in the prognostic region of the nest is needed because the
    ! existing interpolation routines are runtime-optimized for operations that
    ! are needed regularly. Boundary interpolation is executed from the parent grid
    ! to the child grid, whereas the so-called nudging interpolation works from
    ! the local parent grid to the child grid.


    ! Step 1a: prepare boundary interpolation

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_patch(jg)%cells%start_block(grf_bdywidth_c+1)
    i_endblk   = p_patch(jg)%cells%end_block(min_rlcell)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

      DO jk = 1, nlev_soil
        DO jc = i_startidx, i_endidx
          lndvars_par(jc,jk,jb)  = initicon(jg)%sfc_inc%w_so(jc,jk,jb)
        ENDDO
      ENDDO

      jk1 = nlev_soil
      DO jc = i_startidx, i_endidx
        lndvars_par(jc,jk1+1,jb) = initicon(jg)%sfc%sst(jc,jb)
        lndvars_par(jc,jk1+2,jb) = p_parent_ldiag%fr_seaice(jc,jb)
        lndvars_par(jc,jk1+3,jb) = initicon(jg)%sfc_inc%h_snow(jc,jb)
        lndvars_par(jc,jk1+4,jb) = initicon(jg)%sfc_inc%freshsnow(jc,jb)
      ENDDO

      IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) THEN
        DO jc = i_startidx, i_endidx
          lndvars_par(jc,jk1+5,jb) = initicon(jg)%sfc_inc%t_2m(jc,jb)
        ENDDO
      ELSE
        DO jc = i_startidx, i_endidx
          lndvars_par(jc,jk1+5,jb) = 0._wp
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! Step 1b: execute boundary interpolation
    CALL sync_patch_array(SYNC_C,p_patch(jg),lndvars_par)
    CALL interpol_scal_grf (p_pp=p_patch(jg), p_pc=p_pc, p_grf=p_grf_state(jg)%p_dom(i_chidx), nfields=1, lacc=.FALSE., &
      f3din1=lndvars_par, f3dout1=lndvars_chi, lnoshift=.TRUE.                 )


    ! Step 2: Interpolation of fields in the model interior

    ! Step 2a: Copy variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data(p_pat=p_pp%comm_pat_glb_to_loc_c, lacc=.FALSE., RECV=lndvars_lp, SEND=lndvars_par)

    ! Step 2b: Perform interpolation from local parent to child grid

    ! Note: the sync routines cannot be used for the synchronization on the local
    ! parent grid

    IF(l_parallel) CALL exchange_data(p_pat=p_pp%comm_pat_c, lacc=.FALSE., RECV=lndvars_lp)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), 0, &
                                1, 1, f3din1=lndvars_lp, f3dout1=lndvars_chi, overshoot_fac=1.005_wp )
    CALL sync_patch_array(SYNC_C,p_pc,lndvars_chi)



    ! Step 3: Add reference state to thermodynamic variables and copy land fields
    ! from the container arrays to the prognostic variables (for the time being,
    ! all tiles are initialized with the interpolated aggregated variables)
    ! Note that nest boundary points have to be included here to obtain a proper
    ! initialization of all grid points


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_pc%cells%start_block(1)
    i_endblk   = p_pc%cells%end_block(min_rlcell)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, min_rlcell)

      DO jk = 1, nlev_soil
        DO jc = i_startidx, i_endidx
          initicon(jgc)%sfc_inc%w_so(jc,jk,jb) = lndvars_chi(jc,jk,jb)
        ENDDO
      ENDDO

      jk1 = nlev_soil

      DO jc = i_startidx, i_endidx
        initicon(jgc)%sfc%sst(jc,jb)           = lndvars_chi(jc,jk1+1,jb)
        p_child_ldiag%fr_seaice(jc,jb)         = lndvars_chi(jc,jk1+2,jb)
        initicon(jgc)%sfc_inc%h_snow(jc,jb)    = lndvars_chi(jc,jk1+3,jb)
        initicon(jgc)%sfc_inc%freshsnow(jc,jb) = lndvars_chi(jc,jk1+4,jb)

        ! set limits
        p_child_ldiag%fr_seaice(jc,jb) = MAX(0._wp,MIN(1._wp,p_child_ldiag%fr_seaice(jc,jb)))
      ENDDO

      IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) THEN
        DO jc = i_startidx, i_endidx
          initicon(jgc)%sfc_inc%t_2m(jc,jb) = lndvars_chi(jc,jk1+5,jb)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! set information about input source for interpolated fields
    !
    ! w_so (full field: fg, increment: ana(intp))
    var_src = inputInstructions(jg)%ptr%sourceOfVar('w_so')
    var_src = MERGE(kInputSourceFgAnaI, var_src, var_src == kInputSourceBoth)
    CALL inputInstructions(jgc)%ptr%setSource('w_so', var_src)
    !
    ! sst
    ! Branching is necessary here, as sst can either be read via the field t_seasfc of t_so(0)
    IF (inputInstructions(jg)%ptr%fetchStatus('t_seasfc', lIsFg=.FALSE.) == kStateFailedFetch) THEN
      ! since we cannot distinguish between t_so which is read from fg and 
      ! t_so(0)==sst which is read/interpolated from ana, we set
      ! full field: fg, ana(intp), increment: none
      var_src = inputInstructions(jg)%ptr%sourceOfVar('t_so')
      var_src = MERGE(kInputSourceFgAnaI, var_src, var_src == kInputSourceBoth)
      CALL inputInstructions(jgc)%ptr%setSource('t_so', var_src)
    ELSE
      ! full field: ana(intp), increment: none
      CALL inputInstructions(jgc)%ptr%setSource('t_seasfc', kInputSourceAnaI)
    ENDIF
    !
    ! fr_seaice (full field: ana(intp), increment: none)
    var_src = inputInstructions(jg)%ptr%sourceOfVar('fr_seaice')
    var_src = MERGE(kInputSourceAnaI, var_src, var_src == kInputSourceAna)
    CALL inputInstructions(jgc)%ptr%setSource('fr_seaice', var_src)
    !
    ! h_snow (full field: fg, increment: ana(intp))
    var_src = inputInstructions(jg)%ptr%sourceOfVar('h_snow')
    var_src = MERGE(kInputSourceFgAnaI, var_src, var_src == kInputSourceBoth)
    CALL inputInstructions(jgc)%ptr%setSource('h_snow', var_src)
    !
    ! freshsnow (full field: fg, increment: ana(intp))
    var_src = inputInstructions(jg)%ptr%sourceOfVar('freshsnow')
    var_src = MERGE(kInputSourceFgAnaI, var_src, var_src == kInputSourceBoth)
    CALL inputInstructions(jgc)%ptr%setSource('freshsnow', var_src)
    !
    ! t_2m (full field: none, increment: ana(intp))
    IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) THEN
      var_src = inputInstructions(jg)%ptr%sourceOfVar('t_2m')
      var_src = MERGE(kInputSourceAnaI, var_src, var_src == kInputSourceAna)
      CALL inputInstructions(jgc)%ptr%setSource('t_2m', var_src)
    ENDIF


    ! cleanup
    DEALLOCATE(lndvars_par, lndvars_chi, lndvars_lp)

  END SUBROUTINE interpolate_sfcana





  RECURSIVE SUBROUTINE topo_blending_and_fbk(jg, ext_data)

    INTEGER,               INTENT(IN)    :: jg
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    INTEGER :: jgc, jn


    ! Loop over nested domains
    DO jn = 1, p_patch(jg)%n_childdom

      jgc = p_patch(jg)%child_id(jn)

      CALL topography_blending(p_patch(jg), p_patch(jgc), p_grf_state(jg)%p_dom(jn), jn, &
               ext_data(jg)%atm%topography_c, ext_data(jgc)%atm%topography_c  )

      IF (p_patch(jgc)%n_childdom > 0) &
        CALL topo_blending_and_fbk(jgc, ext_data)

    ENDDO

    DO jn = 1, p_patch(jg)%n_childdom

      jgc = p_patch(jg)%child_id(jn)

      IF (lfeedback(jgc) .AND. ifeedback_type == 1) THEN
        CALL topography_feedback(p_patch(jg), jn, ext_data(jg)%atm%topography_c, ext_data(jgc)%atm%topography_c )
      ENDIF

    ENDDO

  END SUBROUTINE topo_blending_and_fbk

  !---------------------------------------------------------------------------
  !>
  !! Computes topography blending for nested domains initialized with real
  !! topography data
  !!
  SUBROUTINE topography_blending(p_pp, p_pc, p_grf, i_chidx, &
               topo_cp, topo_cc)

    ! patch at parent and child level
    TYPE(t_patch),                TARGET, INTENT(IN) :: p_pp
    TYPE(t_patch),                TARGET, INTENT(INOUT) :: p_pc
    ! grf state is needed at parent level only
    TYPE(t_gridref_single_state), TARGET, INTENT(IN) :: p_grf

    ! child domain index
    INTEGER, INTENT(IN) :: i_chidx

    ! topography on cells: p = parent level, c = child level
    REAL(wp), INTENT(IN),     DIMENSION(:,:) :: topo_cp
    REAL(wp), INTENT(INOUT),  DIMENSION(:,:) :: topo_cc

    ! auxiliary fields needed for SR calls (SR's expect 3D fields)
    REAL(wp), TARGET :: z_topo_cp(nproma,1,p_pp%nblks_c)
    REAL(wp)         :: z_topo_cc(nproma,1,p_pc%nblks_c)

    ! another auxiliary needed to handle MPI parallelization, and a related pointer
    REAL(wp), ALLOCATABLE, TARGET :: z_topo_clp (:,:,:)
    REAL(wp), POINTER             :: ptr_topo_cp(:,:,:)

    TYPE(t_gridref_single_state), POINTER :: ptr_grf
    TYPE(t_patch),                POINTER :: ptr_pp
    TYPE(t_int_state),            POINTER :: ptr_int

    INTEGER  :: jgc, jb, jc, nblks_c
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom, i_rlstart, i_rlend
    REAL(wp) :: wfac

    !-------------------------------------------------------------------------

    jgc = p_pc%id

    ptr_grf    => p_grf_state_local_parent(jgc)%p_dom(i_chidx)
    ptr_pp     => p_patch_local_parent(jgc)
    ptr_int    => p_int_state_local_parent(jgc)

    ALLOCATE(z_topo_clp(nproma,1,ptr_pp%nblks_c))

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(2(a,i2))') 'topography blending, domain ',&
        p_pp%id,'  => domain ',p_pc%id
      CALL message('', TRIM(message_text))
    ENDIF

    z_topo_cc(:,:,:) = 0._wp

    nblks_c = p_pp%nblks_c

    ! 1.(a) Copy input topography to auxiliary 3D field
    DO jb = 1, nblks_c
      CALL get_indices_c(p_pp, jb, 1, nblks_c, i_startidx, i_endidx, 1)

      DO jc = i_startidx, i_endidx
        z_topo_cp(jc,1,jb) = topo_cp(jc,jb)
      ENDDO
    ENDDO

    ! 1.(b) Copy this auxiliary field to the local parent in case of MPI parallelization

    CALL exchange_data(p_pat=ptr_pp%comm_pat_glb_to_loc_c, lacc=.FALSE., RECV=z_topo_clp, SEND=z_topo_cp)
    ptr_topo_cp => z_topo_clp

    IF (my_process_is_mpi_parallel()) THEN
      ! synchronization (CALL sync does not work on local parent)
      CALL exchange_data(p_pat=ptr_pp%comm_pat_c, lacc=.FALSE., RECV=ptr_topo_cp)
    END IF

    ! 1.(c) Interpolate coarse topography on fine mesh

    ! Comment: one could avoid the awkward procedure of interpolating boundary zone and
    ! prognostic zone separately by introducing appropriate interpolation routines and communication
    ! patterns (both of which are optimized for runtime tasks). However, those routines
    ! would not be needed anywhere else, and terrain blending is not runtime-critical

    ! Lateral boundary zone
    CALL interpol_scal_grf (p_pp=p_pp, p_pc=p_pc, p_grf=p_grf, nfields=1, lacc=.FALSE., &
      f3din1=z_topo_cp, f3dout1=z_topo_cc, lnoshift=.TRUE.)

    ! Prognostic part of the model domain
    ! Note: in contrast to boundary interpolation, nudging expects the input on the local parent grid
    CALL interpol_scal_nudging (ptr_pp, ptr_int, ptr_grf, 0, 1, 1,     &
                                f3din1=ptr_topo_cp, f3dout1=z_topo_cc  )

    ! 2. Apply terrain blending to cell points
    ! In the boundary interpolation zone (cell rows 1-4), the interpolated values
    ! are taken; in the subsequent 8 rows, linear blending is applied


    i_nchdom = MAX(1,p_pc%n_childdom)

    ! 2.(a) Pure interpolation in boundary interpolation zone
    i_rlstart  = 1
    i_rlend    = grf_bdywidth_c
    i_startblk = p_pc%cells%start_blk(i_rlstart,1)
    i_endblk   = p_pc%cells%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        topo_cc(jc,jb) = z_topo_cc(jc,1,jb)
      ENDDO
    ENDDO

    ! 2.(b) Blending in prognostic part of the nested domain
    i_rlstart  = grf_bdywidth_c + 1
    i_rlend    = min_rlcell
    i_startblk = p_pc%cells%start_blk(i_rlstart,1)
    i_endblk   = p_pc%cells%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        IF (p_pc%cells%refin_ctrl(jc,jb) > 0 .AND. &
            p_pc%cells%refin_ctrl(jc,jb) <= grf_bdywidth_c + nudge_zone_width) THEN

          wfac = REAL(p_pc%cells%refin_ctrl(jc,jb) - grf_bdywidth_c,wp)/ &
            REAL(nudge_zone_width+1,wp)
          topo_cc(jc,jb) = wfac*topo_cc(jc,jb) + (1._wp-wfac)*z_topo_cc(jc,1,jb)

        ENDIF ! on the interior grid points, topo_cc remains unchanged

      ENDDO
    ENDDO

    CALL sync_patch_array(SYNC_C,p_pc,topo_cc)

    DEALLOCATE(z_topo_clp)

  END SUBROUTINE topography_blending


  !---------------------------------------------------------------------------
  !>
  !! Computes topography feedback for nested domains initialized with real
  !! topography data
  !!
  SUBROUTINE topography_feedback(p_pp, i_chidx, topo_cp, topo_cc)

    ! patch at parent level
    TYPE(t_patch),                TARGET, INTENT(IN) :: p_pp

    ! child domain index
    INTEGER, INTENT(IN) :: i_chidx

    ! topography on cells: p = parent level, c = child level
    REAL(wp), INTENT(IN),    DIMENSION(:,:)         :: topo_cc
    REAL(wp), INTENT(INOUT), DIMENSION(:,:), TARGET :: topo_cp

    ! auxiliary field needed for SR calls (SR's expect 3D fields)
    REAL(wp) :: z_topo_cp(nproma,1,p_pp%nblks_c)

    ! another auxiliary to handle MPI parallelization, and related pointer
    REAL(wp), ALLOCATABLE, TARGET :: z_topo_clp(:,:)
    REAL(wp), POINTER             :: ptr_topo_cp(:,:)

    TYPE(t_gridref_state), POINTER :: ptr_grf
    TYPE(t_patch),         POINTER :: ptr_pp

    INTEGER  :: jgc, jb, jc
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom, i_rlstart, i_rlend

    INTEGER, POINTER, DIMENSION(:,:,:) :: iidx, iblk

    !-------------------------------------------------------------------------

    jgc = p_pp%child_id(i_chidx) ! child domain ID

    ALLOCATE(z_topo_clp(nproma,p_patch_local_parent(jgc)%nblks_c))
    ptr_pp      => p_patch_local_parent(jgc)
    ptr_grf     => p_grf_state_local_parent(jgc)
    ptr_topo_cp => z_topo_clp

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(2(a,i2))') 'topography feedback, domain ',&
        p_pp%child_id(i_chidx),'  => domain ',p_pp%id
        CALL message('', TRIM(message_text))
    ENDIF

    ! Start/End block in the parent domain
    i_rlstart  = grf_fbk_start_c
    i_rlend    = min_rlcell_int

    i_startblk = ptr_pp%cells%start_blk(i_rlstart,i_chidx)
    i_endblk   = ptr_pp%cells%end_blk(i_rlend,i_chidx)

    iidx => ptr_pp%cells%child_idx
    iblk => ptr_pp%cells%child_blk

    ! Perform feedback for cell points in the prognostic part of the nest overlap area
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        ptr_topo_cp(jc,jb) =                                                  &
          topo_cc(iidx(jc,jb,1),iblk(jc,jb,1))*ptr_grf%fbk_wgt_bln(jc,jb,1) + &
          topo_cc(iidx(jc,jb,2),iblk(jc,jb,2))*ptr_grf%fbk_wgt_bln(jc,jb,2) + &
          topo_cc(iidx(jc,jb,3),iblk(jc,jb,3))*ptr_grf%fbk_wgt_bln(jc,jb,3) + &
          topo_cc(iidx(jc,jb,4),iblk(jc,jb,4))*ptr_grf%fbk_wgt_bln(jc,jb,4)
      ENDDO
    ENDDO

    ! Attention: the feedback communication pattern is defined on the local parent (ptr_pp), ...
    CALL exchange_data(p_pat=ptr_pp%comm_pat_loc_to_glb_c_fbk, lacc=.FALSE., RECV=topo_cp, SEND=ptr_topo_cp, &
      &                l_recv_exists=.TRUE.)

    ! ... whereas the subsequent halo communication needs to be done on the parent grid (p_pp)
    CALL sync_patch_array(SYNC_C,p_pp,topo_cp)


    ! Copy fed-back topography to z_topo_cp to prepare interpolation to vertices
    i_nchdom = MAX(1,p_pp%n_childdom)
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_pp%cells%start_blk(i_rlstart,1)
    i_endblk   = p_pp%cells%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        z_topo_cp(jc,1,jb) = topo_cp(jc,jb)

      ENDDO
    ENDDO

    DEALLOCATE(z_topo_clp)

  END SUBROUTINE topography_feedback


END MODULE mo_nh_init_nest_utils

