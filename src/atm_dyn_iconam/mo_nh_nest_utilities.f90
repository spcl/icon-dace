!
!  This module contains the routines needed for nesting in the nonhydrostatic.
!  version.
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

MODULE mo_nh_nest_utilities
  !
  !
  USE mo_kind,                ONLY: wp, vp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_model_domain,        ONLY: t_patch, t_grid_cells, t_grid_edges, p_patch_local_parent, p_patch
  USE mo_grid_config,         ONLY: n_dom, n_dom_start
  USE mo_intp_data_strc,      ONLY: t_int_state, p_int_state, p_int_state_local_parent
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state, p_grf_state, p_grf_state_local_parent
  USE mo_gridref_config,      ONLY: grf_intmethod_c, grf_intmethod_e, grf_intmethod_ct, grf_scalfbk, grf_tracfbk
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf, interpol_vec_grf, interpol2_vec_grf
  USE mo_grf_nudgintp,        ONLY: interpol_scal_nudging, interpol_vec_nudging
  USE mo_grf_ubcintp,         ONLY: interpol_scal_ubc,interpol_vec_ubc
  USE mo_dynamics_config,     ONLY: nnow, nsav1, nnow_rcf
  USE mo_parallel_config,     ONLY: nproma, p_test_run, cpu_min_nproma
  USE mo_run_config,          ONLY: ltransport, msg_level, ntracer, lvert_nest, iqv, iqc, iforcing
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydro_state,      ONLY: p_nh_state
  USE mo_nwp_phy_state,       ONLY: prm_diag
  USE mo_prepadv_types,       ONLY: t_prepare_adv
  USE mo_nonhydrostatic_config,ONLY: ndyn_substeps_var
  USE mo_atm_phy_nwp_config,  ONLY: iprog_aero
  USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int, min_rledge_int, min_rlcell, min_rledge
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c, grf_bdyintp_end_c, grf_bdywidth_c, &
    &                               grf_bdywidth_e, grf_nudgintp_start_c, grf_nudgintp_start_e, &
    &                               grf_nudge_start_c
  USE mo_mpi,                 ONLY: my_process_is_mpi_parallel
#ifdef _OPENACC  
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif
  USE mo_communication,       ONLY: exchange_data, exchange_data_mult
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, sync_patch_array, sync_patch_array_mult
  USE mo_physical_constants,  ONLY: rd, cvd_o_rd, p0ref, vtmpc1, rd_o_cpd
  USE mo_limarea_config,      ONLY: latbc_config
  USE mo_initicon_types,      ONLY: t_pi_atm
  USE mo_async_latbc_types,   ONLY: t_latbc_state
  USE mo_advection_config,    ONLY: advection_config
  USE mo_nudging_config,      ONLY: nudging_config, indg_type, ithermdyn_type
  USE mtime,                  ONLY: datetimeToString, MAX_DATETIME_STR_LEN

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: compute_tendencies
  PUBLIC :: boundary_interpolation
  PUBLIC :: complete_nesting_setup
  PUBLIC :: save_progvars
  PUBLIC :: prep_bdy_nudging
  PUBLIC :: nest_boundary_nudging
  PUBLIC :: prep_rho_bdy_nudging
  PUBLIC :: density_boundary_nudging
  PUBLIC :: limarea_nudging_latbdy, limarea_nudging_upbdy
  PUBLIC :: intp_nestubc_nudging

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nh_nest_utilities'
CONTAINS


  !>
  !! Computes correction term needed to use perturbation density for
  !! lateral boundary nudging (for use with 1-way nesting).
  !!
  SUBROUTINE complete_nesting_setup(p_patch, p_patch_local_parent, p_grf_state_local_parent, p_nh_state)

    TYPE(t_patch),                 INTENT(IN)    :: p_patch(1:)
    TYPE(t_patch),         TARGET, INTENT(IN)    :: p_patch_local_parent(n_dom_start+1:)
    TYPE(t_gridref_state), TARGET, INTENT(IN)    :: p_grf_state_local_parent(n_dom_start+1:)
    TYPE(t_nh_state),              INTENT(INOUT) :: p_nh_state(1:)

    ! local
    TYPE(t_patch), POINTER     :: p_pp => NULL()   ! local parent patch of child domain
    REAL(vp),      ALLOCATABLE :: z_rho_ref(:,:,:)

    INTEGER :: jg, ji, jgc, jb, jc, jk, jks
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: nlev, nshift
    INTEGER :: nlev_c                  ! number of vertical levels for child domain
    INTEGER :: ist
    INTEGER,  POINTER :: iidx(:,:,:), iblk(:,:,:)
    REAL(wp), POINTER :: p_fbkwgt(:,:,:)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':complete_nesting_setup'

    DO jg = 1, n_dom-1

      nlev = p_patch(jg)%nlev

      ! skip the following, if no child domain exists
      IF (p_patch(jg)%n_childdom == 0) CYCLE

      ! compute correction term needed to use perturbation density for boundary nudging
      ! correction term lives on local-parent grid of domain jgc
      !
      DO ji = 1, p_patch(jg)%n_childdom

        jgc    =  p_patch(jg)%child_id(ji)

        nlev_c = p_patch(jgc)%nlev
        nshift = p_patch(jgc)%nshift

        p_fbkwgt => p_grf_state_local_parent(jgc)%fbk_wgt_bln
        p_pp     => p_patch_local_parent(jgc)

        iidx  => p_pp%cells%child_idx
        iblk  => p_pp%cells%child_blk

        i_startblk = p_pp%cells%start_blk(grf_nudgintp_start_c+1,ji)
        i_endblk   = p_pp%cells%end_blk(min_rlcell_int,ji)

        ALLOCATE(p_nh_state(jgc)%metrics%rho_ref_corr(nproma, nlev_c, p_pp%nblks_c), &
          &      z_rho_ref(nproma, nlev, p_pp%nblks_c), STAT=ist)
        IF (ist /= SUCCESS) THEN
          CALL finish(routine, 'allocation of rho_ref_corr and z_rho_ref failed')
        ENDIF
        z_rho_ref(:,:,:) = 0._vp
        p_nh_state(jgc)%metrics%rho_ref_corr(:,:,:) = 0._vp

        ! copy rho_ref of domain jg to local-parent grid of domain jgc
        CALL exchange_data(p_pp%comm_pat_glb_to_loc_c, RECV=z_rho_ref, &
          &                SEND=p_nh_state(jg)%metrics%rho_ref_mc)

        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            grf_nudgintp_start_c+1, min_rlcell_int)

          DO jk = 1, nlev_c
            jks = jk + nshift
            DO jc = i_startidx, i_endidx

              p_nh_state(jgc)%metrics%rho_ref_corr(jc,jk,jb) = - z_rho_ref(jc,jks,jb) + (           &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1)+ &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2)+ &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3)+ &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)  )

            ENDDO
          ENDDO
        ENDDO

        !
        ! rho_ref_corr is not defined with ADD_VAR, and needs to be explicitly copied to DEVICE
        !
        !$ACC ENTER DATA COPYIN(p_nh_state(jgc)%metrics%rho_ref_corr)

        DEALLOCATE(z_rho_ref, STAT=ist)
        IF (ist /= SUCCESS) THEN
          CALL finish(routine, 'deallocation of z_rho_ref failed')
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE complete_nesting_setup


  !-------------------------------------------------------------------------
  !
  !>
  !! Saves the dynamic prognostic variables needed afterwards for computing
  !! the lateral boundary tendencies for nested domains
  !!
  SUBROUTINE save_progvars (jg,p_nh_prog,p_nh_save)

    INTEGER,         INTENT(IN)    :: jg
    TYPE(t_nh_prog), INTENT(IN)    :: p_nh_prog
    TYPE(t_nh_prog), INTENT(INOUT) :: p_nh_save

    ! local variables
    !
    INTEGER :: ib, jb, ic, jc, ie, je, jk, jshift
    INTEGER :: nlev, nlevp1
    INTEGER :: nproma_bdyintp, nblks_bdyintp, npromz_bdyintp, nlen

    TYPE(t_gridref_state), POINTER :: p_grf

    !-----------------------------------------------------------------------


    p_grf  => p_grf_state(jg)

    ! number of vertical levels
    nlev   = p_patch(jg)%nlev
    nlevp1 = p_patch(jg)%nlevp1

    ! for dynamic nproma blocking
    nproma_bdyintp = cpu_min_nproma(nproma,256)

!$OMP PARALLEL PRIVATE(nblks_bdyintp,npromz_bdyintp)

    ! cell-based variables

    ! parameters for dynamic nproma blocking
    nblks_bdyintp  = INT(p_grf%npoints_bdyintp_src_c/nproma_bdyintp)
    npromz_bdyintp = MOD(p_grf%npoints_bdyintp_src_c,nproma_bdyintp)
    IF (npromz_bdyintp > 0) THEN
      nblks_bdyintp = nblks_bdyintp + 1
    ELSE
      npromz_bdyintp = nproma_bdyintp
    ENDIF

!$OMP DO PRIVATE(ib,jb,nlen,ic,jc,jk,jshift) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = 1, nblks_bdyintp
      IF (ib == nblks_bdyintp) THEN
        nlen = npromz_bdyintp
      ELSE
        nlen = nproma_bdyintp
      ENDIF
      jshift = (ib-1)*nproma_bdyintp

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG VECTOR PRIVATE(jc, jb)
      DO ic = jshift+1, jshift+nlen
        jc = p_grf%idxlist_bdyintp_src_c(ic)
        jb = p_grf%blklist_bdyintp_src_c(ic)
!DIR$ IVDEP
        DO jk = 1, nlev
#else
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(jc, jb)
      DO jk = 1, nlev
!$NEC ivdep
        DO ic = jshift+1, jshift+nlen
          jc = p_grf%idxlist_bdyintp_src_c(ic)
          jb = p_grf%blklist_bdyintp_src_c(ic)
#endif

          p_nh_save%w(jc,jk,jb)       = p_nh_prog%w(jc,jk,jb)
          p_nh_save%rho(jc,jk,jb)     = p_nh_prog%rho(jc,jk,jb)
          p_nh_save%theta_v(jc,jk,jb) = p_nh_prog%theta_v(jc,jk,jb)

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR PRIVATE(jc, jb)
      DO ic = jshift+1, jshift+nlen
        jc = p_grf%idxlist_bdyintp_src_c(ic)
        jb = p_grf%blklist_bdyintp_src_c(ic)
        p_nh_save%w(jc,nlevp1,jb)  = p_nh_prog%w(jc,nlevp1,jb)
      ENDDO
      !$ACC END PARALLEL

    ENDDO
    !$ACC WAIT(1)
!$OMP END DO

    ! edge-based variables

    ! parameters for dynamic nproma blocking
    nblks_bdyintp  = INT(p_grf%npoints_bdyintp_src_e/nproma_bdyintp)
    npromz_bdyintp = MOD(p_grf%npoints_bdyintp_src_e,nproma_bdyintp)
    IF (npromz_bdyintp > 0) THEN
      nblks_bdyintp = nblks_bdyintp + 1
    ELSE
      npromz_bdyintp = nproma_bdyintp
    ENDIF

!$OMP DO PRIVATE(ib,jb,nlen,ie,je,jk,jshift) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = 1, nblks_bdyintp
      IF (ib == nblks_bdyintp) THEN
        nlen = npromz_bdyintp
      ELSE
        nlen = nproma_bdyintp
      ENDIF
      jshift = (ib-1)*nproma_bdyintp

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
      DO ie = jshift+1, jshift+nlen
        je = p_grf%idxlist_bdyintp_src_e(ie)
        jb = p_grf%blklist_bdyintp_src_e(ie)
!DIR$ IVDEP
        DO jk = 1, nlev
#else
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(je, jb)
      DO jk = 1, nlev
!$NEC ivdep
        DO ie = jshift+1, jshift+nlen
          je = p_grf%idxlist_bdyintp_src_e(ie)
          jb = p_grf%blklist_bdyintp_src_e(ie)
#endif

          p_nh_save%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb)

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE save_progvars

  !-------------------------------------------------------------------------
  !
  !
  !
  !>
  !! Computes the time tendencies of the prognostic variables needed for
  !! interpolation to the lateral boundaries of the nested domains
  !! In addition, compute upper boundary condition for vertical nesting.
  !!
  SUBROUTINE compute_tendencies (jg,n_new,n_now,n_new_rcf,n_now_rcf,&
    &                            rdt,rdt_mflx)


    INTEGER,  INTENT(IN) :: jg  ! domain ID
    ! Time levels from which tendencies are computed
    INTEGER,  INTENT(IN) ::  n_new,n_now
    ! Time levels from which tracer-tendencies are computed
    INTEGER,  INTENT(IN) ::  n_new_rcf,n_now_rcf

    ! Inverse value of time step needed for computing the tendencies
    REAL(wp), INTENT(IN) ::  rdt
    ! Inverse value of time step needed for computing the mass flux tendencies
    REAL(wp), INTENT(IN) ::  rdt_mflx


    ! local variables
    !
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx,       &
      ib, jb, ic, jc, ie, je, jk, jt, js, nshift, jk_start, jshift
    INTEGER :: nlev, nlevp1           !< number of full and half levels
    INTEGER :: nproma_bdyintp, nblks_bdyintp, npromz_bdyintp, nlen, ntracer_bdyintp, nsubs

    REAL(wp) :: rnsubs        ! inverse: 1/nsubs
    REAL(wp) :: rdt_ubc       ! reciprocal time step for the computation of UBC time tendencies

    ! Switch to control if the child domain is vertically nested and therefore
    ! needs interpolation of upper boundary conditions
    LOGICAL :: l_child_vertnest

    TYPE(t_gridref_state), POINTER :: p_grf

    TYPE(t_nh_state), POINTER :: p_nh
    TYPE(t_nh_prog),  POINTER :: p_prog_now
    TYPE(t_nh_prog),  POINTER :: p_prog_new
    TYPE(t_nh_prog),  POINTER :: p_prog_now_rcf
    TYPE(t_nh_prog),  POINTER :: p_prog_new_rcf

    !-----------------------------------------------------------------------

    nsubs  = ndyn_substeps_var(jg)
    rnsubs = 1._wp/REAL(nsubs,wp)

    ! dt_ubc = dt * (nsubs-1)/nsubs
    ! i.e. we assume that the data for each substep are centered in time
    ! here we compute the reciprocal
    rdt_ubc = rdt*REAL(nsubs,wp)/REAL(MAX(1,nsubs-1),wp)


    p_grf          => p_grf_state(jg)
    p_nh           => p_nh_state(jg)
    p_prog_now     => p_nh%prog(n_now)
    p_prog_new     => p_nh%prog(n_new)
    p_prog_now_rcf => p_nh%prog(n_now_rcf)
    p_prog_new_rcf => p_nh%prog(n_new_rcf)

    ! number of vertical levels
    nlev   = p_patch(jg)%nlev
    nlevp1 = p_patch(jg)%nlevp1

    IF (advection_config(jg)%iadv_tke == 1) THEN
      ntracer_bdyintp = ntracer-1
    ELSE
      ntracer_bdyintp = ntracer
    ENDIF

    ! determine if upper boundary interpolation is needed
    IF (lvert_nest .AND. (p_patch(jg)%nshift_child > 0)) THEN
      l_child_vertnest = .TRUE.
      nshift = p_patch(jg)%nshift_child + 1
    ELSE
      l_child_vertnest = .FALSE.
      nshift = 1
    ENDIF

    jk_start = nshift ! start index for tendency computation

    ! for dynamic nproma blocking
    nproma_bdyintp = cpu_min_nproma(nproma,256)


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,nblks_bdyintp,npromz_bdyintp)

    IF (l_child_vertnest) THEN ! Compute upper boundary condition for nested domain

      ! cell-based variables
      i_startblk = p_patch(jg)%cells%start_block(grf_bdywidth_c+1)
      i_endblk   = p_patch(jg)%cells%end_block(min_rlcell_int)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,js) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          p_nh%diag%w_int          (jc,jb,nsubs+1) = 0._wp
          p_nh%diag%theta_v_ic_int (jc,jb,nsubs+1) = 0._wp
          p_nh%diag%rho_ic_int     (jc,jb,nsubs+1) = 0._wp
          p_nh%diag%mflx_ic_int    (jc,jb,nsubs+1) = 0._wp
        ENDDO
        !
        ! compute time averages over the nsubs dynamics substeps and store them at index nsubs+1
        !$ACC LOOP SEQ
        DO js = 1, nsubs
!DIR$ IVDEP
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
            p_nh%diag%w_int(jc,jb,nsubs+1)          = p_nh%diag%w_int(jc,jb,nsubs+1) +           &
              p_nh%diag%w_int(jc,jb,js)
            !
            p_nh%diag%theta_v_ic_int(jc,jb,nsubs+1) = p_nh%diag%theta_v_ic_int(jc,jb,nsubs+1) + &
              p_nh%diag%theta_v_ic_int(jc,jb,js)
            !
            p_nh%diag%rho_ic_int(jc,jb,nsubs+1)     = p_nh%diag%rho_ic_int(jc,jb,nsubs+1) +   &
              p_nh%diag%rho_ic_int(jc,jb,js)
            !
            p_nh%diag%mflx_ic_int(jc,jb,nsubs+1)    = p_nh%diag%mflx_ic_int(jc,jb,nsubs+1) +     &
              p_nh%diag%mflx_ic_int(jc,jb,js)
          ENDDO
        ENDDO

        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          p_nh%diag%w_int(jc,jb,nsubs+1)             = p_nh%diag%w_int(jc,jb,nsubs+1)          * rnsubs
          p_nh%diag%theta_v_ic_int(jc,jb,nsubs+1)    = p_nh%diag%theta_v_ic_int(jc,jb,nsubs+1) * rnsubs
          p_nh%diag%rho_ic_int(jc,jb,nsubs+1)        = p_nh%diag%rho_ic_int(jc,jb,nsubs+1)     * rnsubs
          p_nh%diag%mflx_ic_int(jc,jb,nsubs+1)       = p_nh%diag%mflx_ic_int(jc,jb,nsubs+1)    * rnsubs
        ENDDO

        ! Compute time tendencies to obtain second order in time accuracy for child nest UBC, 
        ! and store them at index nsubs+2.
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          p_nh%diag%w_int(jc,jb,nsubs+2)         = rdt_ubc                               &
            &   * (p_nh%diag%w_int(jc,jb,nsubs) - p_nh%diag%w_int(jc,jb,1))
          !
          p_nh%diag%theta_v_ic_int(jc,jb,nsubs+2)= rdt_ubc                               &
            &   * (p_nh%diag%theta_v_ic_int(jc,jb,nsubs) - p_nh%diag%theta_v_ic_int(jc,jb,1))
          !
          p_nh%diag%rho_ic_int(jc,jb,nsubs+2)    = rdt_ubc                               &
            &   * (p_nh%diag%rho_ic_int(jc,jb,nsubs) - p_nh%diag%rho_ic_int(jc,jb,1))
          !
          p_nh%diag%mflx_ic_int(jc,jb,nsubs+2)   = rdt_ubc                               &
            &   * (p_nh%diag%mflx_ic_int(jc,jb,nsubs) - p_nh%diag%mflx_ic_int(jc,jb,1))
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!$OMP END DO


      ! edge-based variables
      i_startblk = p_patch(jg)%edges%start_block(grf_bdywidth_e+1)
      i_endblk   = p_patch(jg)%edges%end_block(min_rledge_int-2)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, grf_bdywidth_e+1, min_rledge_int-2)

!DIR$ IVDEP
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR
        DO je = i_startidx, i_endidx
          p_nh%diag%vn_ie_int(je,2,jb) = rdt_ubc*(p_nh%diag%vn_ie(je,nshift,jb)-p_nh%diag%vn_ie_int(je,1,jb))
          p_nh%diag%vn_ie_int(je,1,jb) = 0.5_wp*(p_nh%diag%vn_ie_int(je,1,jb)+p_nh%diag%vn_ie(je,nshift,jb))
        ENDDO
        !$ACC END PARALLEL

      ENDDO
!$OMP END DO

    ENDIF ! l_child_vertnest


    !
    ! Computation of tendencies for lateral boundary interpolation
    !
    ! cell-based variables

    ! parameters for dynamic nproma blocking
    nblks_bdyintp  = INT(p_grf%npoints_bdyintp_src_c/nproma_bdyintp)
    npromz_bdyintp = MOD(p_grf%npoints_bdyintp_src_c,nproma_bdyintp)
    IF (npromz_bdyintp > 0) THEN
      nblks_bdyintp = nblks_bdyintp + 1
    ELSE
      npromz_bdyintp = nproma_bdyintp
    ENDIF

!$OMP DO PRIVATE(ib,jb,nlen,ic,jc,jk,jt,jshift) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = 1, nblks_bdyintp
      IF (ib == nblks_bdyintp) THEN
        nlen = npromz_bdyintp
      ELSE
        nlen = nproma_bdyintp
      ENDIF
      jshift = (ib-1)*nproma_bdyintp

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
      DO ic = jshift+1, jshift+nlen
        jc = p_grf%idxlist_bdyintp_src_c(ic)
        jb = p_grf%blklist_bdyintp_src_c(ic)
!DIR$ IVDEP
        DO jk = jk_start, nlev
#else
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(jc, jb)
      DO jk = jk_start, nlev
!$NEC ivdep
        DO ic = jshift+1, jshift+nlen
          jc = p_grf%idxlist_bdyintp_src_c(ic)
          jb = p_grf%blklist_bdyintp_src_c(ic)
#endif

          p_nh%diag%grf_tend_rho(jc,jk,jb) = &
            ( p_prog_new%rho(jc,jk,jb) - p_prog_now%rho(jc,jk,jb) )*rdt

          p_nh%diag%grf_tend_thv(jc,jk,jb) = &
            ( p_prog_new%theta_v(jc,jk,jb) - p_prog_now%theta_v(jc,jk,jb) )*rdt

          ! the div field carries perturbation density for use in SR boundary_interpolation
          p_nh%diag%div(jc,jk,jb) = &
            p_prog_now%rho(jc,jk,jb) - p_nh%metrics%rho_ref_mc(jc,jk,jb)

          ! the dpres_mc field carries perturbation potential temperature for use in SR boundary_interpolation
          p_nh%diag%dpres_mc(jc,jk,jb) = &
            p_prog_now%theta_v(jc,jk,jb) - p_nh%metrics%theta_ref_mc(jc,jk,jb)

          p_nh%diag%grf_tend_w(jc,jk,jb) = &
            ( p_prog_new%w(jc,jk,jb) - p_prog_now%w(jc,jk,jb) )*rdt
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR PRIVATE(jc, jb)
      DO ic = jshift+1, jshift+nlen
        jc = p_grf%idxlist_bdyintp_src_c(ic)
        jb = p_grf%blklist_bdyintp_src_c(ic)
        p_nh%diag%grf_tend_w(jc,nlevp1,jb) = &
          ( p_prog_new%w(jc,nlevp1,jb) - p_prog_now%w(jc,nlevp1,jb) )*rdt
      ENDDO
      !$ACC END PARALLEL

      IF (ltransport) THEN

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
        DO ic = jshift+1, jshift+nlen
          jc = p_grf%idxlist_bdyintp_src_c(ic)
          jb = p_grf%blklist_bdyintp_src_c(ic)
          DO jt = 1,ntracer_bdyintp
!DIR$ IVDEP
            DO jk = jk_start, nlev
#else
        !$ACC LOOP GANG VECTOR COLLAPSE(3) PRIVATE(jc, jb)
        DO jt = 1,ntracer_bdyintp
          DO jk = jk_start, nlev
!$NEC ivdep
            DO ic = jshift+1, jshift+nlen
              jc = p_grf%idxlist_bdyintp_src_c(ic)
              jb = p_grf%blklist_bdyintp_src_c(ic)
#endif

                p_nh%diag%grf_tend_tracer(jc,jk,jb,jt) =                 &
                  &            ( p_prog_new_rcf%tracer(jc,jk,jb,jt)      &
                  &            -  p_prog_now_rcf%tracer(jc,jk,jb,jt) )*rdt
            ENDDO
          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ENDIF

    ENDDO
!$OMP END DO

    ! edge-based variables

    ! parameters for dynamic nproma blocking
    nblks_bdyintp  = INT(p_grf%npoints_bdyintp_src_e/nproma_bdyintp)
    npromz_bdyintp = MOD(p_grf%npoints_bdyintp_src_e,nproma_bdyintp)
    IF (npromz_bdyintp > 0) THEN
      nblks_bdyintp = nblks_bdyintp + 1
    ELSE
      npromz_bdyintp = nproma_bdyintp
    ENDIF

!$OMP DO PRIVATE(ib,jb,nlen,ie,je,jk,jshift) ICON_OMP_DEFAULT_SCHEDULE
    DO ib = 1, nblks_bdyintp
      IF (ib == nblks_bdyintp) THEN
        nlen = npromz_bdyintp
      ELSE
        nlen = nproma_bdyintp
      ENDIF
      jshift = (ib-1)*nproma_bdyintp

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
      DO ie = jshift+1, jshift+nlen
        je = p_grf%idxlist_bdyintp_src_e(ie)
        jb = p_grf%blklist_bdyintp_src_e(ie)
!DIR$ IVDEP
        DO jk = jk_start, nlev
#else
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(je, jb)
      DO jk = jk_start, nlev
!$NEC ivdep
        DO ie = jshift+1, jshift+nlen
          je = p_grf%idxlist_bdyintp_src_e(ie)
          jb = p_grf%blklist_bdyintp_src_e(ie)
#endif
          p_nh%diag%grf_tend_vn(je,jk,jb) = &
            ( p_prog_new%vn(je,jk,jb) - p_prog_now%vn(je,jk,jb) )*rdt
          p_nh%diag%grf_tend_mflx(je,jk,jb) = &
            ( p_nh%diag%mass_fl_e(je,jk,jb) - p_nh%diag%mass_fl_e_sv(je,jk,jb) )*rdt_mflx
        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE compute_tendencies

  !-------------------------------------------------------------------------
  !
  !>
  !! Interpolates time tendencies of prognostic variables to the lateral boundary
  !! of a refined mesh.
  !! In addition, interpolates prognostic variables to child upper boundary  
  !! for vertical nesting.
  !!
  SUBROUTINE boundary_interpolation (jg,jgc,ntp_dyn,ntc_dyn,ntp_tr,ntc_tr, &
    p_patch, p_nh_state, prep_adv, p_grf_state)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':boundary_interpolation'

    INTEGER,  INTENT(IN)    :: jg, jgc      ! domain ID of parent and child grid

    ! Parent and child time levels for dynamical variables and tracers
    INTEGER,  INTENT(IN)    :: ntp_dyn, ntc_dyn, ntp_tr, ntc_tr

    TYPE(t_patch)        , INTENT(   IN), TARGET :: p_patch(:)
    TYPE(t_nh_state)     , INTENT(INOUT), TARGET :: p_nh_state(:)
    TYPE(t_prepare_adv)  , INTENT(INOUT), TARGET :: prep_adv(:)
    TYPE(t_gridref_state), INTENT(   IN), TARGET :: p_grf_state(:)

    ! local variables

    TYPE(t_nh_diag), POINTER           :: p_diagp => NULL()
    TYPE(t_nh_diag), POINTER           :: p_diagc  => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhp_dyn => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhc_dyn => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhp_tr  => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhc_tr  => NULL()
    TYPE(t_patch), POINTER             :: p_pp => NULL()
    TYPE(t_patch), POINTER             :: p_pc => NULL()
    TYPE(t_gridref_state), POINTER     :: p_grf => NULL()
    TYPE(t_grid_cells), POINTER        :: p_gcp => NULL()
    TYPE(t_prepare_adv), POINTER       :: prep_advp => NULL()
    TYPE(t_prepare_adv), POINTER       :: prep_advc => NULL()

    INTEGER :: i_startblk              ! start block
    INTEGER :: i_endblk                ! end index
    INTEGER :: i_startidx              ! start index
    INTEGER :: i_endidx                ! end index

    INTEGER :: jb, jc, jk, jt, ic      ! loop indices

    INTEGER :: nlev_c                  ! number of full levels (child domain)

    INTEGER :: i_chidx, i_sbc, i_ebc
    INTEGER :: ntracer_bdyintp, nsubs

    REAL(wp) :: aux3dp(nproma,ntracer+8,p_patch(jg)%nblks_c), &
      aux3dc(nproma,ntracer+8,p_patch(jgc)%nblks_c), &
      theta_prc(nproma,p_patch(jgc)%nlev,p_patch(jgc)%nblks_c), &
      rho_prc(nproma,p_patch(jgc)%nlev,p_patch(jgc)%nblks_c)

    ! Switch to control if the child domain is vertically nested and therefore
    ! needs interpolation of upper boundary conditions
    LOGICAL :: l_child_vertnest

    LOGICAL :: l_limit(2*ntracer)

#ifdef _OPENMP
    INTEGER :: num_threads_omp, omp_get_max_threads
#endif

    !-----------------------------------------------------------------------

    !$ACC DATA CREATE(aux3dp, aux3dc, theta_prc, rho_prc) IF(i_am_accel_node)
    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') '========= Interpolate:',jg,' =>',jgc
      CALL message(routine, message_text)
    ENDIF

#ifdef _OPENMP
    num_threads_omp = omp_get_max_threads()
#endif

    p_diagp       => p_nh_state(jg)%diag
    p_diagc       => p_nh_state(jgc)%diag
    p_nhp_dyn     => p_nh_state(jg)%prog(ntp_dyn)
    p_nhc_dyn     => p_nh_state(jgc)%prog(ntc_dyn)
    p_nhp_tr      => p_nh_state(jg)%prog(ntp_tr)
    p_nhc_tr      => p_nh_state(jgc)%prog(ntc_tr)
    p_grf         => p_grf_state(jg)
    p_pp          => p_patch(jg)
    p_pc          => p_patch(jgc)
    p_gcp         => p_patch(jg)%cells

    prep_advp     => prep_adv(jg)
    prep_advc     => prep_adv(jgc)

    i_chidx = p_patch(jgc)%parent_child_index

    nsubs    = ndyn_substeps_var(jg)

    ! number of vertical levels (child domain)
    nlev_c   = p_pc%nlev

    ! determine if upper boundary interpolation is needed
    IF (lvert_nest .AND. (p_pp%nshift_child > 0)) THEN
      l_child_vertnest = .TRUE.
    ELSE
      l_child_vertnest = .FALSE.
    ENDIF

    ! exclude TKE from boundary interpolation if it is only vertically advected
    IF (advection_config(jg)%iadv_tke == 1) THEN
      ntracer_bdyintp = ntracer-1
    ELSE
      ntracer_bdyintp = ntracer
    ENDIF

    ! Perform interpolations to upper nest boundary needed for vertical nesting
    IF (l_child_vertnest) THEN

      IF (p_test_run) THEN
        !$ACC KERNELS ASYNC(1) IF(i_am_accel_node)
        aux3dp = 0._wp
        aux3dc = 0._wp
        !$ACC END KERNELS
      ENDIF

      CALL sync_patch_array(SYNC_E,p_pp,p_diagp%vn_ie_int)

      CALL interpol_vec_ubc (p_pp, p_pc, p_grf%p_dom(i_chidx), &
        p_diagp%vn_ie_int, p_diagc%vn_ie_ubc)

      ! Start and end blocks for which interpolation is needed
      i_startblk = p_pp%cells%start_block(0)
      i_endblk   = p_pp%cells%end_block(min_rlcell_int)

      ! For back-copying at child level
      i_sbc      = p_pc%cells%start_block(grf_nudge_start_c-2)
      i_ebc      = p_pc%cells%end_block(min_rlcell_int)

      IF (ltransport) THEN

!$OMP PARALLEL DO PRIVATE(jb,i_startidx,i_endidx,jc,jt) ICON_OMP_DEFAULT_SCHEDULE
        DO jb =  i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            0, min_rlcell_int)
        
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
!$NEC ivdep
          DO jc = i_startidx, i_endidx
            aux3dp(jc,1:2,jb) = p_diagp%w_int         (jc,jb,nsubs+1:nsubs+2)
            aux3dp(jc,3:4,jb) = p_diagp%theta_v_ic_int(jc,jb,nsubs+1:nsubs+2)
            aux3dp(jc,5:6,jb) = p_diagp%rho_ic_int    (jc,jb,nsubs+1:nsubs+2)
            aux3dp(jc,7:8,jb) = p_diagp%mflx_ic_int   (jc,jb,nsubs+1:nsubs+2)
          ENDDO
          !$ACC LOOP SEQ
          DO jt = 1, ntracer
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO jc = i_startidx, i_endidx
              aux3dp(jc,8+jt,jb) = prep_advp%q_int(jc,jt,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL

        ENDDO
!$OMP END PARALLEL DO

        CALL sync_patch_array(SYNC_C,p_pp,aux3dp)

        CALL interpol_scal_ubc (p_pc, p_grf%p_dom(i_chidx),  &
          ntracer+8, aux3dp, aux3dc)

!$OMP PARALLEL DO PRIVATE(jb,i_startidx,i_endidx,jc,jt) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_sbc, i_ebc

          CALL get_indices_c(p_pc, jb, i_sbc, i_ebc, i_startidx, i_endidx, &
            grf_nudge_start_c-2, min_rlcell_int)

!$NEC ivdep
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
            p_diagc%w_ubc(jc,jb,1:2)          = aux3dc(jc,1:2,jb)
            p_diagc%theta_v_ic_ubc(jc,jb,1:2) = aux3dc(jc,3:4,jb)
            p_diagc%rho_ic_ubc(jc,jb,1:2)     = aux3dc(jc,5:6,jb)
            p_diagc%mflx_ic_ubc(jc,jb,1:2)    = aux3dc(jc,7:8,jb)
          ENDDO

          !$ACC LOOP SEQ
          DO jt = 1, ntracer
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO jc = i_startidx, i_endidx
              prep_advc%q_ubc(jc,jt,jb) = aux3dc(jc,jt+8,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL

        ENDDO
!$OMP END PARALLEL DO

      ELSE

!$OMP PARALLEL DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb =  i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            0, min_rlcell_int)
!$NEC ivdep
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            aux3dp(jc,1:2,jb) = p_diagp%w_int         (jc,jb,nsubs+1:nsubs+2)
            aux3dp(jc,3:4,jb) = p_diagp%theta_v_ic_int(jc,jb,nsubs+1:nsubs+2)
            aux3dp(jc,5:6,jb) = p_diagp%rho_ic_int    (jc,jb,nsubs+1:nsubs+2)
            aux3dp(jc,7:8,jb) = p_diagp%mflx_ic_int   (jc,jb,nsubs+1:nsubs+2)
          ENDDO
          !$ACC END PARALLEL

        ENDDO
!$OMP END PARALLEL DO

        CALL sync_patch_array(SYNC_C,p_pp,aux3dp)

        CALL interpol_scal_ubc(p_pc, p_grf%p_dom(i_chidx), 8, aux3dp, aux3dc)

!$OMP PARALLEL DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_sbc, i_ebc

          CALL get_indices_c(p_pc, jb, i_sbc, i_ebc, i_startidx, i_endidx, &
            grf_nudge_start_c-2, min_rlcell_int)

!$NEC ivdep
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            p_diagc%w_ubc(jc,jb,1:2)          = aux3dc(jc,1:2,jb)
            p_diagc%theta_v_ic_ubc(jc,jb,1:2) = aux3dc(jc,3:4,jb)
            p_diagc%rho_ic_ubc(jc,jb,1:2)     = aux3dc(jc,5:6,jb)
            p_diagc%mflx_ic_ubc(jc,jb,1:2)    = aux3dc(jc,7:8,jb)
          ENDDO
          !$ACC END PARALLEL

        ENDDO
!$OMP END PARALLEL DO

      ENDIF
    ENDIF  ! l_child_vertnest

    !
    ! Lateral boundary interpolation of cell-based dynamical variables
    !
    IF (grf_intmethod_c == 1) THEN ! tendency copying for all cell-based variables

      CALL exchange_data(p_pc%comm_pat_interpolation_c, &
        RECV=p_diagc%grf_tend_rho,     &
        SEND=p_diagp%grf_tend_rho)

      CALL exchange_data(p_pc%comm_pat_interpolation_c, &
        RECV=p_diagc%grf_tend_thv,     &
        SEND=p_diagp%grf_tend_thv)

      ! exchange_data should also work for w because it determines the
      ! vertical dimension with UBOUND
      CALL exchange_data(p_pc%comm_pat_interpolation_c, &
        RECV=p_diagc%grf_tend_w,       &
        SEND=p_diagp%grf_tend_w)

      ! grf_intmethod_c = 2, use gradient at cell center for interpolation
    ELSE IF (grf_intmethod_c == 2) THEN

      ! Interpolation of temporal tendencies, full w, perturbation density (stored in div)
      !  and perturbationvirtual potential temperature (stored in dpres_mc)
      CALL interpol_scal_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 6,         &
        p_diagp%grf_tend_rho, p_diagc%grf_tend_rho,  &
        p_diagp%grf_tend_thv, p_diagc%grf_tend_thv,  &
        p_diagp%grf_tend_w,   p_diagc%grf_tend_w,    &
        p_nhp_dyn%w,          p_nhc_dyn%w,           &
        p_nh_state(jg)%diag%div, rho_prc,            &
        p_nh_state(jg)%diag%dpres_mc, theta_prc      )

      ! Start and end blocks for which interpolation is needed
      i_startblk = p_pc%cells%start_block(1)
      i_endblk   = p_pc%cells%end_block(grf_bdywidth_c)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb =  i_startblk, i_endblk

        CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
          1, grf_bdywidth_c)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev_c
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_nhc_dyn%rho(jc,jk,jb) = rho_prc(jc,jk,jb) + &
              p_nh_state(jgc)%metrics%rho_ref_mc(jc,jk,jb)
            p_nhc_dyn%theta_v(jc,jk,jb) = theta_prc(jc,jk,jb) + &
              p_nh_state(jgc)%metrics%theta_ref_mc(jc,jk,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ENDDO
!$OMP END DO

      ! The following index list contains the halo points of the lateral boundary
      ! cells. These have to be copied as well in order for rho and theta to be
      ! synchronized.
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR PRIVATE(jc, jb)
!$OMP DO PRIVATE(ic,jb,jc,jk)
      DO ic = 1, p_nh_state(jgc)%metrics%bdy_halo_c_dim

        jb = p_nh_state(jgc)%metrics%bdy_halo_c_blk(ic)
        jc = p_nh_state(jgc)%metrics%bdy_halo_c_idx(ic)
!DIR$ IVDEP
        DO jk = 1, nlev_c
          p_nhc_dyn%rho(jc,jk,jb) = rho_prc(jc,jk,jb) + &
            p_nh_state(jgc)%metrics%rho_ref_mc(jc,jk,jb)
          p_nhc_dyn%theta_v(jc,jk,jb) = theta_prc(jc,jk,jb) + &
            p_nh_state(jgc)%metrics%theta_ref_mc(jc,jk,jb)

        ENDDO
      ENDDO
      !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL
    ENDIF

    ! Lateral boundary interpolation of cell based tracer variables
    IF (ltransport .AND. grf_intmethod_ct == 1) THEN

      ! Start and end blocks for which interpolation is needed
      i_startblk = p_gcp%start_block(grf_bdyintp_start_c)
      i_endblk   = p_gcp%end_block(grf_bdyintp_end_c)

      CALL exchange_data_mult(p_pc%comm_pat_interpolation_c, ntracer_bdyintp, ntracer_bdyintp*nlev_c, &
        RECV4D=p_diagc%grf_tend_tracer(:,:,:,1:ntracer_bdyintp),                                      &
        SEND4D=p_diagp%grf_tend_tracer(:,:,:,1:ntracer_bdyintp))


    ELSE IF (ltransport .AND. grf_intmethod_ct == 2) THEN

      ! Apply positive definite limiter on full tracer fields but not on tendencies
      l_limit(1:ntracer_bdyintp) = .FALSE.
      l_limit(ntracer_bdyintp+1:2*ntracer_bdyintp) = .TRUE.

      CALL interpol_scal_grf ( p_pp, p_pc, p_grf%p_dom(i_chidx), 2*ntracer_bdyintp, &
        f4din1 =  p_diagp%grf_tend_tracer(:,:,:,1:ntracer_bdyintp),                 &
        f4dout1 = p_diagc%grf_tend_tracer(:,:,:,1:ntracer_bdyintp),                 &
        f4din2  = p_nhp_tr%tracer(:,:,:,1:ntracer_bdyintp),                         &
        f4dout2 = p_nhc_tr%tracer(:,:,:,1:ntracer_bdyintp),                         &
        llimit_nneg=l_limit)

    ENDIF

    IF (ltransport .AND. iprog_aero >= 1) THEN
     CALL interpol_scal_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 1, prm_diag(jg)%aerosol,    &
                             prm_diag(jgc)%aerosol, llimit_nneg=(/.TRUE./), lnoshift=.TRUE.)
    ENDIF

    ! Lateral boundary interpolation of edge-based variables  (velocity components)
    IF (grf_intmethod_e == 2) THEN

      CALL interpol_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), p_diagp%grf_tend_vn, p_diagc%grf_tend_vn)

    ELSE IF (grf_intmethod_e == 4) THEN

      CALL interpol2_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 1, p_diagp%grf_tend_vn, p_diagc%grf_tend_vn)

    ELSE IF (grf_intmethod_e == 6) THEN

      CALL interpol2_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 3,                              &
        p_diagp%grf_tend_vn, p_diagc%grf_tend_vn, prep_advp%mass_flx_me, prep_advc%mass_flx_me, &
        p_diagp%grf_tend_mflx, p_diagc%grf_tend_mflx )

    ENDIF
    !$ACC WAIT
    !$ACC END DATA ! aux3dp, aux3dc, theta_prc, rho_prc

  END SUBROUTINE boundary_interpolation


  !>
  !! This routine prepares boundary nudging for use with 1-way nesting.
  !!
  !! The following steps are executed:
  !! 1. Mapping of parent grid prognostic variables to intermediate grid having 
  !!    the horizontal resolution of the parent grid, but sharing
  !!    the domain decomposition and vertical dimension with the child grid.
  !! 2. a) Interpolation/Averaging of child grid variables to intermediate grid. 
  !!    b) Computation of differences between mapped parent-grid values and averaged child grid
  !!    variables.
  !! 3. Interpolation of difference fields to the child grid
  !!
  SUBROUTINE prep_bdy_nudging(jgp, jg)

    CHARACTER(len=*), PARAMETER :: routine = modname//':prep_bdy_nudging'


    INTEGER, INTENT(IN) :: jg   ! child grid level
    INTEGER, INTENT(IN) :: jgp  ! parent grid level


    ! local variables
    !
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_child_prog  => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog_rcf => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_child_prog_rcf  => NULL()
    TYPE(t_nh_diag),    POINTER     :: p_diag        => NULL()
    TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
    TYPE(t_grid_edges), POINTER     :: p_gep => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
    TYPE(t_int_state), POINTER      :: p_int => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()

    ! Indices
    INTEGER :: jb, jc, jk, jt, je, js, i_chidx, i_startblk, i_endblk, &
      i_startidx, i_endidx, istartblk_c, istartblk_e
    INTEGER :: nlev_c, nlev_p
    INTEGER :: nshift      !< difference between upper boundary of parent or feedback-parent
                           !< domain and upper boundary of child domain (in terms
                           !< of vertical levels)
    INTEGER :: ntracer_nudge !< number of tracers to be nudged

    ! Local arrays for interpolated parent-grid values, and difference fields. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_thv, diff_thv
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_rho, diff_rho
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_vn , diff_vn
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_w  , diff_w
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: parent_tr , diff_tr

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk, ieidx, ieblk
    LOGICAL :: l_parallel
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt, p_fbkwgt_tr, p_fbkwgt_v
    !-----------------------------------------------------------------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') '1-way nesting: == Boundary nudging:',jg
      CALL message(routine, message_text)
    ENDIF

    l_parallel = my_process_is_mpi_parallel()

    IF (iforcing > 1) THEN  ! tracers represent moisture variables
      ntracer_nudge = 3     ! take only QV, QC and QI - nudging precip variables does not make much sense
    ELSE
      ntracer_nudge = ntracer
    ENDIF

    p_parent_prog     => p_nh_state(jgp)%prog(nsav1(jgp))
    p_child_prog      => p_nh_state(jg)%prog(nnow(jg))
    p_parent_prog_rcf => p_nh_state(jgp)%prog(nnow_rcf(jgp))
    p_child_prog_rcf  => p_nh_state(jg)%prog(nnow_rcf(jg))
    p_diag            => p_nh_state(jg)%diag
    p_pc              => p_patch(jg)

    p_grf => p_grf_state_local_parent(jg)
    p_int => p_int_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_gep => p_patch_local_parent(jg)%edges
    p_pp  => p_patch_local_parent(jg)

    i_chidx  = p_pc%parent_child_index

    ! number of full levels of child domain
    nlev_c   = p_pc%nlev

    ! number of full/half levels of parent domain
    nlev_p   = p_pp%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift
    js     = nshift

    ! Please note: In the parallel case
    ! - lower bound must be 1 due to synchronization calls
    ! - upper bound must be nblks_c/e to include halo cells/edges
    ! - this doesn't cost extra memory since p_patch_local_parent
    !   only includes the cells/edges really needed

    ! Value of i_startblk needed for subroutine call for parent-to-child interpolation
    istartblk_c = 1

    ALLOCATE(parent_thv  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c),  &
      diff_thv    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c),  &
      parent_rho  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c),  &
      diff_rho    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c),  &
      parent_w    (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c),  &
      diff_w      (nproma, nlev_c+1, p_patch_local_parent(jg)%nblks_c) )

    IF(ltransport) &
      ALLOCATE(parent_tr(nproma, nlev_p, p_patch_local_parent(jg)%nblks_c, ntracer_nudge),&
      &        diff_tr  (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c, ntracer_nudge) )

    ! Value of i_startblk needed for subroutine call for parent-to-child interpolation
    istartblk_e = 1

    ALLOCATE(parent_vn  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_e), &
      &      diff_vn    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_e)  )

    !$ACC DATA CREATE(parent_thv, diff_thv, parent_rho, diff_rho, parent_w, diff_w) &
    !$ACC   CREATE(parent_vn, diff_vn) IF(i_am_accel_node)
    !$ACC DATA CREATE(parent_tr, diff_tr) IF(i_am_accel_node .AND. ltransport)
    
    ! Set pointers to index and coefficient fields for cell-based variables
    iidx  => p_gcp%child_idx
    iblk  => p_gcp%child_blk
    ieidx => p_gep%child_idx
    ieblk => p_gep%child_blk

    IF (grf_scalfbk == 1) THEN
      p_fbkwgt    => p_grf%fbk_wgt_aw
    ELSE
      p_fbkwgt    => p_grf%fbk_wgt_bln
    ENDIF
    IF (grf_tracfbk == 1) THEN
      p_fbkwgt_tr => p_grf%fbk_wgt_aw
    ELSE
      p_fbkwgt_tr => p_grf%fbk_wgt_bln
    ENDIF
    p_fbkwgt_v  => p_grf%fbk_wgt_e

    ! 1st step: Copy prognostic variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 3, 3*nlev_p, &
      RECV1=parent_rho, SEND1=p_parent_prog%rho,         &
      RECV2=parent_thv, SEND2=p_parent_prog%theta_v,     &
      RECV3=parent_w,   SEND3=p_parent_prog%w            )

    CALL exchange_data(p_pp%comm_pat_glb_to_loc_e,     &
      RECV=parent_vn, SEND=p_parent_prog%vn)

    IF (ltransport) &
      CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, ntracer_nudge, ntracer_nudge*nlev_p, &
      RECV4D=parent_tr, SEND4D=p_parent_prog_rcf%tracer(:,:,:,1:ntracer_nudge))

    ! 2nd step: perform feedback from refined grid to intermediate grid and compute differences

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! a) cell-based variables
    ! Start/End block in the parent domain
    i_startblk = p_gcp%start_blk(grf_nudgintp_start_c+1,i_chidx)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
        grf_nudgintp_start_c+1, min_rlcell_int)

      ! initialize diff_w at surface with zero
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR
      DO jc = 1, nproma
        diff_w(jc,nlev_c+1,jb) = 0._wp
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
#endif

          diff_thv(jc,jk,jb) = parent_thv(jc,jk+js,jb) - (                           &
            p_child_prog%theta_v(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%theta_v(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%theta_v(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%theta_v(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)   )

          diff_rho(jc,jk,jb) = parent_rho(jc,jk+js,jb) - (                       &
            p_child_prog%rho(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%rho(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%rho(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%rho(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4))+ &
            p_nh_state(jg)%metrics%rho_ref_corr(jc,jk,jb)

          diff_w(jc,jk,jb) = parent_w(jc,jk+js,jb) - (                         &
            p_child_prog%w(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%w(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%w(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%w(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)   )

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      ! Tracers
      IF (ltransport) THEN

        DO jt = 1, ntracer_nudge

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1, nlev_c
#else
          DO jk = 1, nlev_c
            DO jc = i_startidx, i_endidx
#endif

              diff_tr(jc,jk,jb,jt) = parent_tr(jc,jk+js,jb,jt) - (                                &
                p_child_prog_rcf%tracer(iidx(jc,jb,1),jk,iblk(jc,jb,1),jt)*p_fbkwgt_tr(jc,jb,1) + &
                p_child_prog_rcf%tracer(iidx(jc,jb,2),jk,iblk(jc,jb,2),jt)*p_fbkwgt_tr(jc,jb,2) + &
                p_child_prog_rcf%tracer(iidx(jc,jb,3),jk,iblk(jc,jb,3),jt)*p_fbkwgt_tr(jc,jb,3) + &
                p_child_prog_rcf%tracer(iidx(jc,jb,4),jk,iblk(jc,jb,4),jt)*p_fbkwgt_tr(jc,jb,4)   )
            ENDDO
          ENDDO
          !$ACC END PARALLEL

        ENDDO

      ENDIF

    ENDDO
!$OMP END DO

    ! b) velocity
    ! Start/End block in the parent domain
    i_startblk = p_gep%start_blk(grf_nudgintp_start_e+2,i_chidx)
    i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
        grf_nudgintp_start_e+2, min_rledge_int)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
      DO jk = 1, nlev_c
        DO je = i_startidx, i_endidx
#endif

          diff_vn(je,jk,jb) = parent_vn(je,jk+js,jb) - (                            &
            p_child_prog%vn(ieidx(je,jb,1),jk,ieblk(je,jb,1))*p_fbkwgt_v(je,jb,1) + &
            p_child_prog%vn(ieidx(je,jb,2),jk,ieblk(je,jb,2))*p_fbkwgt_v(je,jb,2)   )

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! Interpolate differences to child grid; the differences are stored in the grf_tend fields

    ! interpol_vec/scal_nudging needs all fields with full boundaries, so the arrays
    ! have to be sync'd before calling these routines.
    ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
    ! since the arrays don't start with lower bound 1 in the non parallel case!

    ! Synchronization is needed after the interpolation step because the nudging tendencies are applied outside 
    ! the dynamical core. This is needed for the scalars for reasons of mass consistency, but is also done for the
    ! wind tendencies because this turns out to improve noise filtering

    IF(l_parallel) CALL exchange_data(p_pp%comm_pat_e, diff_vn)
    CALL interpol_vec_nudging (p_pp, p_pc, p_int, p_grf%p_dom(i_chidx),   &
      &                        0, istartblk_e, diff_vn,p_diag%grf_tend_vn )
    CALL sync_patch_array(SYNC_E,p_pc,p_diag%grf_tend_vn)

    IF(l_parallel) CALL exchange_data_mult(p_pp%comm_pat_c, 3, 3*nlev_c+1, &
      recv1=diff_thv, recv2=diff_rho, recv3=diff_w)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), 0, 3, istartblk_c,          &
      &                         f3din1=diff_thv, f3dout1=p_diag%grf_tend_thv,                  &
      &                         f3din2=diff_rho, f3dout2=p_diag%grf_tend_rho,                  &
      &                         f3din3=diff_w,   f3dout3=p_diag%grf_tend_w                     )
    CALL sync_patch_array_mult(SYNC_C,p_pc,3,p_diag%grf_tend_thv,p_diag%grf_tend_rho,  &
      p_diag%grf_tend_w)

    IF (ltransport) THEN
      IF(l_parallel) CALL exchange_data_mult(p_pp%comm_pat_c, ntracer_nudge, ntracer_nudge*nlev_c, recv4d=diff_tr)

      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx),                   &
        &                         0, ntracer_nudge, istartblk_c, f4din=diff_tr,        &
        &                         f4dout=p_diag%grf_tend_tracer(:,:,:,1:ntracer_nudge) )
      CALL sync_patch_array_mult(SYNC_C,p_pc,ntracer_nudge,f4din=p_diag%grf_tend_tracer(:,:,:,1:ntracer_nudge))
    ENDIF

    !$ACC WAIT
    !$ACC END DATA ! parent_thv, diff_thv, ...
    !$ACC END DATA ! parent_tr, diff_tr

    DEALLOCATE(parent_thv, diff_thv, parent_rho, diff_rho, parent_w, diff_w, parent_vn, diff_vn)

    IF(ltransport) DEALLOCATE(parent_tr, diff_tr)


  END SUBROUTINE prep_bdy_nudging

  !>
  !! This routine prepares boundary nudging for density only (for use with 2-way nesting)
  !!
  !! The following steps are executed:
  !! 1. Mapping of parent grid prognostic variables to intermediate grid having 
  !!    the horizontal resolution of the parent grid, but sharing
  !!    the domain decomposition and vertical dimension with the child grid.
  !! 2. a) Interpolation/Averaging of child grid variables to intermediate grid. 
  !!    b) Computation of differences between mapped parent-grid values and averaged child grid
  !!    variables.
  !! 3. Interpolation of difference fields to the child grid
  !!
  SUBROUTINE prep_rho_bdy_nudging(jgp, jg)


    INTEGER, INTENT(IN) :: jg   ! child grid level
    INTEGER, INTENT(IN) :: jgp  ! parent grid level


    ! local variables

    TYPE(t_nh_prog),    POINTER     :: p_parent_prog => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_child_prog  => NULL()
    TYPE(t_nh_diag),    POINTER     :: p_diag        => NULL()
    TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
    TYPE(t_int_state), POINTER      :: p_int => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()

    ! Indices
    INTEGER :: jb, jc, jk, js, i_chidx, i_startblk, i_endblk, &
      i_startidx, i_endidx, istartblk_c
    INTEGER :: nlev_c, nlev_p
    INTEGER :: nshift      !< difference between upper boundary of parent or feedback-parent
    !< domain and upper boundary of child domain (in terms
    !< of vertical levels)

    ! Local arrays for interpolated parent-grid values, and difference fields. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_rho, diff_rho

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
    LOGICAL :: l_parallel
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt
    !-----------------------------------------------------------------------


    l_parallel = my_process_is_mpi_parallel()

    p_parent_prog     => p_nh_state(jgp)%prog(nsav1(jgp))
    p_child_prog      => p_nh_state(jg)%prog(nnow(jg))
    p_diag            => p_nh_state(jg)%diag
    p_pc              => p_patch(jg)

    p_grf => p_grf_state_local_parent(jg)
    p_int => p_int_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)

    i_chidx  = p_pc%parent_child_index

    ! number of full levels of child domain
    nlev_c   = p_pc%nlev

    ! number of full/half levels of parent domain
    nlev_p   = p_pp%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift
    js     = nshift

    ! Please note: In the parallel case
    ! - lower bound must be 1 due to synchronization calls
    ! - upper bound must be nblks_c/e to include halo cells/edges
    ! - this doesn't cost extra memory since p_patch_local_parent
    !   only includes the cells/edges really needed

    ! Value of i_startblk needed for subroutine call for parent-to-child interpolation
    istartblk_c = 1

    ALLOCATE(parent_rho  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c), &
      &      diff_rho    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c)  )
    !$ACC DATA CREATE(parent_rho, diff_rho) IF(i_am_accel_node)

    ! Set pointers to index and coefficient fields for cell-based variables
    iidx  => p_gcp%child_idx
    iblk  => p_gcp%child_blk

    IF (grf_scalfbk == 1) THEN
      p_fbkwgt    => p_grf%fbk_wgt_aw
    ELSE
      p_fbkwgt    => p_grf%fbk_wgt_bln
    ENDIF

    ! 1st step: Copy prognostic variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data(p_pp%comm_pat_glb_to_loc_c,        &
      &                RECV=parent_rho, SEND=p_parent_prog%rho )

    ! 2nd step: perform feedback from refined grid to intermediate grid and compute differences

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! Start/End block in the parent domain
    i_startblk = p_gcp%start_blk(grf_nudgintp_start_c+1,i_chidx)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
        grf_nudgintp_start_c+1, min_rlcell_int)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
#endif

          diff_rho(jc,jk,jb) = parent_rho(jc,jk+js,jb) - (                       &
            p_child_prog%rho(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%rho(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%rho(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%rho(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4))+ &
            p_nh_state(jg)%metrics%rho_ref_corr(jc,jk,jb)

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL


    ! Interpolate differences to child grid; the differences are stored in the grf_tend fields

    ! interpol_scal_nudging needs all fields with full boundaries, so the arrays
    ! have to be sync'd before calling these routines.
    ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
    ! since the arrays don't start with lower bound 1 in the non parallel case!

    ! Synchronization is needed after the interpolation step for cell-based variables because for
    ! those, the nudging tendencies are applied outside the dynamical core for reasons of mass consistency

    IF(l_parallel) CALL exchange_data(p_pp%comm_pat_c, diff_rho)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), 0, 1, istartblk_c, &
      &                         f3din1=diff_rho, f3dout1=p_diag%grf_tend_rho                   )
    CALL sync_patch_array(SYNC_C,p_pc,p_diag%grf_tend_rho)

    !$ACC WAIT(1)
    !$ACC END DATA
    DEALLOCATE(parent_rho, diff_rho)


  END SUBROUTINE prep_rho_bdy_nudging



  !>
  !! This routine executes LATERAL boundary nudging for the limited-area mode.
  !!
  SUBROUTINE limarea_nudging_latbdy (p_patch, p_prog, ptr_tracer, p_metrics, p_diag, &
                                  p_int, p_latbc_const, p_latbc_old, p_latbc_new, lc1, lc2)

    TYPE(t_patch),      INTENT(IN)    :: p_patch
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog
    REAL(wp), CONTIGUOUS, INTENT(inout) :: ptr_tracer(:,:,:,:)
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag),    INTENT(INOUT) :: p_diag
    TYPE(t_int_state),  INTENT(IN)    :: p_int

    ! alternative input data, either for constant or time-dependent lateral boundary conditions
    TYPE(t_nh_prog),    INTENT(IN), OPTIONAL :: p_latbc_const
    TYPE(t_pi_atm),     INTENT(IN), OPTIONAL :: p_latbc_old, p_latbc_new
    REAL(wp),           INTENT(IN), OPTIONAL :: lc1     ! time interpolation weight
    REAL(wp),           INTENT(IN), OPTIONAL :: lc2     ! time interpolation weight
    
    ! local variables
    INTEGER  :: jb, jc, jk, je, ic, nlev
    INTEGER  :: jg, nshift
    REAL(wp) :: wfac_old, wfac_new, pres, temp, qv, tempv_inc, pres_inc
    REAL(wp) :: rho_tend, thv_tend, vn_tend, qv_tend
    REAL(wp) :: rd_o_cvd, rd_o_p0ref
    REAL(wp) :: zrho, ztheta_v, zexner
    !
    CHARACTER(len=*), PARAMETER :: routine = 'limarea_nudging_latbdy'


    ! domain id
    jg = p_patch%id

    ! number of full levels of child domain
    nlev = p_patch%nlev

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    IF (PRESENT(lc1) .AND. PRESENT(lc2)) THEN
      wfac_old = lc1
      wfac_new = lc2  
    ELSE
      wfac_old = -999._wp
      wfac_new = -999._wp
    ENDIF

    IF (nudging_config(jg)%ltype(indg_type%ubn)) THEN
      ! Upper boundary nudging is switched on
      nshift = nudging_config(jg)%ilev_end
    ELSE
      nshift = 0
    ENDIF

    IF (PRESENT(p_latbc_const) .AND. (PRESENT(p_latbc_old) .OR. PRESENT(p_latbc_new))) THEN

      CALL finish(TRIM(routine),'conflicting arguments')

    ELSE IF (PRESENT(p_latbc_const)) THEN ! Mode for constant lateral boundary data
      !$ACC DATA PRESENT(p_metrics, p_prog, p_latbc_const, p_int)

      ! compute differences between lateral boundary data and prognostic variables

!$OMP PARALLEL
!$OMP DO PRIVATE(jk,jc,jb,ic,rho_tend,thv_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
      DO ic = 1, p_metrics%nudge_c_dim
        jc = p_metrics%nudge_c_idx(ic)
        jb = p_metrics%nudge_c_blk(ic)
!DIR$ IVDEP
        DO jk = 1, nlev
#else
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(jc, jb, thv_tend, rho_tend, zrho, ztheta_v, zexner)
      DO jk = 1, nlev
!$NEC ivdep
        DO ic = 1, p_metrics%nudge_c_dim
          jc = p_metrics%nudge_c_idx(ic)
          jb = p_metrics%nudge_c_blk(ic)
#endif
          thv_tend = p_latbc_const%theta_v(jc,jk,jb) - p_prog%theta_v(jc,jk,jb)
          rho_tend = p_latbc_const%rho    (jc,jk,jb) - p_prog%rho    (jc,jk,jb)

#ifdef _OPENACC
            ! local variables zrho, ztheta_v and zexner are used because of GPU pgi compiler bug
            ! this section should be revisited once bug is resolved.
            ! BUG-FIX for nvfort 21.3.0 (Daint) fond and solved by Xavier Lapillonne
            zrho     = p_prog%rho(jc,jk,jb)     + p_int%nudgecoeff_c(jc,jb)*rho_tend
            ztheta_v = p_prog%theta_v(jc,jk,jb) + p_int%nudgecoeff_c(jc,jb)*thv_tend
            zexner   = EXP(rd_o_cvd*LOG(rd_o_p0ref*zrho*ztheta_v))
            p_prog%rho(jc,jk,jb)     = zrho
            p_prog%theta_v(jc,jk,jb) = ztheta_v
            p_prog%exner(jc,jk,jb)   = zexner
#else
            p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + p_int%nudgecoeff_c(jc,jb)*rho_tend
            p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + p_int%nudgecoeff_c(jc,jb)*thv_tend
            p_prog%exner(jc,jk,jb)   = EXP(rd_o_cvd*LOG(rd_o_p0ref*p_prog%rho(jc,jk,jb)*p_prog%theta_v(jc,jk,jb)))
#endif

        ENDDO
      ENDDO
      !$ACC END PARALLEL
!$OMP END DO

!$OMP DO PRIVATE(jk,je,jb,ic,vn_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
      DO ic = 1, p_metrics%nudge_e_dim
        je = p_metrics%nudge_e_idx(ic)
        jb = p_metrics%nudge_e_blk(ic)
!DIR$ IVDEP
        DO jk = 1, nlev
#else
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(jc, jb, vn_tend)
      DO jk = 1, nlev
!$NEC ivdep
        DO ic = 1, p_metrics%nudge_e_dim
          je = p_metrics%nudge_e_idx(ic)
          jb = p_metrics%nudge_e_blk(ic)
#endif
          vn_tend = p_latbc_const%vn(je,jk,jb) - p_prog%vn(je,jk,jb)

          p_prog%vn(je,jk,jb) = p_prog%vn(je,jk,jb) + p_int%nudgecoeff_e(je,jb)*vn_tend

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC WAIT
      !$ACC END DATA

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ELSE IF (PRESENT(p_latbc_old) .AND. PRESENT(p_latbc_new)) THEN ! Mode for time-dependent lateral boundary data

      !$ACC DATA PRESENT(p_metrics%nudge_c_idx, p_metrics%nudge_c_blk) &
      !$ACC   PRESENT(ptr_tracer, p_prog%theta_v, p_prog%rho, p_prog%exner, p_prog%vn) &
      !$ACC   PRESENT(p_diag%temp, p_diag%pres, p_diag%tempv) &
      !$ACC   PRESENT(p_latbc_old, p_latbc_new) &
      !$ACC   PRESENT(p_latbc_old%pres, p_latbc_old%temp, p_latbc_old%qv) &
      !$ACC   PRESENT(p_latbc_old%theta_v, p_latbc_old%rho, p_latbc_old%vn) &
      !$ACC   PRESENT(p_latbc_new%pres, p_latbc_new%temp, p_latbc_new%qv) &
      !$ACC   PRESENT(p_latbc_new%theta_v, p_latbc_new%rho, p_latbc_new%vn) &
      !$ACC   PRESENT(p_metrics%nudge_c_dim, p_int%nudgecoeff_c) &
      !$ACC   PRESENT(p_int, p_prog, p_diag, p_metrics) &
      !$ACC   PRESENT(p_metrics%nudge_e_idx, p_metrics%nudge_e_blk)


      ! check if interpolation weights are provided
      IF (.NOT. PRESENT(lc1) .OR. .NOT. PRESENT(lc2)) THEN
        CALL finish(routine,'missing interpolation weights wfac_old, wfac_new')
      ENDIF

      ! compute differences between lateral boundary data and prognostic variables


!$OMP PARALLEL
      IF (latbc_config%nudge_hydro_pres .AND. ltransport) THEN
!$OMP DO PRIVATE(jk,jc,jb,ic,pres,temp,qv,tempv_inc,pres_inc,rho_tend,thv_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
        DO ic = 1, p_metrics%nudge_c_dim
          jc = p_metrics%nudge_c_idx(ic)
          jb = p_metrics%nudge_c_blk(ic)
!DIR$ IVDEP
          DO jk = nshift+1, nlev
#else
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) &
        !$ACC   PRIVATE(jc, jb, pres, temp, qv, tempv_inc, pres_inc, thv_tend, rho_tend, zrho, ztheta_v, zexner)
        DO jk = nshift+1, nlev
!$NEC ivdep
          DO ic = 1, p_metrics%nudge_c_dim
            jc = p_metrics%nudge_c_idx(ic)
            jb = p_metrics%nudge_c_blk(ic)
#endif
            pres = wfac_old*p_latbc_old%pres(jc,jk,jb) + wfac_new*p_latbc_new%pres(jc,jk,jb)
            temp = wfac_old*p_latbc_old%temp(jc,jk,jb) + wfac_new*p_latbc_new%temp(jc,jk,jb)
            qv   = wfac_old*p_latbc_old%qv(jc,jk,jb)   + wfac_new*p_latbc_new%qv(jc,jk,jb)

            tempv_inc = (temp-p_diag%temp(jc,jk,jb))*(1._wp+vtmpc1*qv) + &
               (qv-ptr_tracer(jc,jk,jb,iqv))*vtmpc1*temp
            pres_inc  = pres-p_diag%pres(jc,jk,jb)

            thv_tend = tempv_inc/p_prog%exner(jc,jk,jb) - rd_o_cpd*p_prog%theta_v(jc,jk,jb)/pres*pres_inc
            rho_tend = ( pres_inc/p_diag%tempv(jc,jk,jb) - &
              tempv_inc*p_diag%pres(jc,jk,jb)/p_diag%tempv(jc,jk,jb)**2 )/rd

#ifdef _OPENACC
            ! local variables zrho, ztheta_v and zexner are used because of GPU pgi compiler bug
            ! this section should be revisited once bug is resolved.
            zrho     = p_prog%rho(jc,jk,jb)     + p_int%nudgecoeff_c(jc,jb)*rho_tend
            ztheta_v = p_prog%theta_v(jc,jk,jb) + p_int%nudgecoeff_c(jc,jb)*thv_tend
            zexner   = EXP(rd_o_cvd*LOG(rd_o_p0ref*zrho*ztheta_v))
            p_prog%rho(jc,jk,jb)     = zrho
            p_prog%theta_v(jc,jk,jb) = ztheta_v
            p_prog%exner(jc,jk,jb)   = zexner
#else
            p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + p_int%nudgecoeff_c(jc,jb)*rho_tend
            p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + p_int%nudgecoeff_c(jc,jb)*thv_tend
            p_prog%exner(jc,jk,jb)   = EXP(rd_o_cvd*LOG(rd_o_p0ref*p_prog%rho(jc,jk,jb)*p_prog%theta_v(jc,jk,jb)))
#endif
          ENDDO
        ENDDO
        !$ACC END PARALLEL
!$OMP END DO

      ELSE

!$OMP DO PRIVATE(jk,jc,jb,ic,rho_tend,thv_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
        DO ic = 1, p_metrics%nudge_c_dim
          jc = p_metrics%nudge_c_idx(ic)
          jb = p_metrics%nudge_c_blk(ic)
!DIR$ IVDEP
          DO jk = nshift+1, nlev
#else
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(jc, jb, thv_tend, rho_tend)
        DO jk = nshift+1, nlev
!$NEC ivdep
          DO ic = 1, p_metrics%nudge_c_dim
            jc = p_metrics%nudge_c_idx(ic)
            jb = p_metrics%nudge_c_blk(ic)
#endif
            thv_tend = wfac_old*p_latbc_old%theta_v(jc,jk,jb) + wfac_new*p_latbc_new%theta_v(jc,jk,jb) - p_prog%theta_v(jc,jk,jb)
            rho_tend = wfac_old*p_latbc_old%rho(jc,jk,jb) + wfac_new*p_latbc_new%rho(jc,jk,jb)- p_prog%rho(jc,jk,jb)

            p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + p_int%nudgecoeff_c(jc,jb)*rho_tend
            p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + p_int%nudgecoeff_c(jc,jb)*thv_tend
            p_prog%exner(jc,jk,jb)   = EXP(rd_o_cvd*LOG(rd_o_p0ref*p_prog%rho(jc,jk,jb)*p_prog%theta_v(jc,jk,jb)))

          ENDDO
        ENDDO
        !$ACC END PARALLEL
!$OMP END DO

      ENDIF

      IF (ltransport) THEN ! apply QV nudging in subsaturated regions
!$OMP DO PRIVATE(jk,jc,jb,ic,qv_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
        DO ic = 1, p_metrics%nudge_c_dim
          jc = p_metrics%nudge_c_idx(ic)
          jb = p_metrics%nudge_c_blk(ic)
!DIR$ IVDEP
          DO jk = 1, nlev
#else
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(jc, jb, qv_tend)
        DO jk = 1, nlev
!$NEC ivdep
          DO ic = 1, p_metrics%nudge_c_dim
            jc = p_metrics%nudge_c_idx(ic)
            jb = p_metrics%nudge_c_blk(ic)
#endif
            qv_tend = wfac_old*p_latbc_old%qv(jc,jk,jb) + wfac_new*p_latbc_new%qv(jc,jk,jb) - ptr_tracer(jc,jk,jb,iqv)

            ! Suppress positive nudging tendencies in saturated (=cloudy) regions in order to avoid runaway effects
            qv_tend = MERGE(MIN(0._wp,qv_tend), qv_tend, ptr_tracer(jc,jk,jb,iqc) > 1.e-10_wp)

            ! using a weaker nudging coefficient for QV than for thermodynamic variables turned out to have a slightly
            ! beneficial impact on forecast quality
            ptr_tracer(jc,jk,jb,iqv) = ptr_tracer(jc,jk,jb,iqv) + 0.5_wp*p_int%nudgecoeff_c(jc,jb)*qv_tend

          ENDDO
        ENDDO
        !$ACC END PARALLEL
!$OMP END DO

      ENDIF

!$OMP DO PRIVATE(jk,je,jb,ic,vn_tend) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
        DO ic = 1, p_metrics%nudge_e_dim
          je = p_metrics%nudge_e_idx(ic)
          jb = p_metrics%nudge_e_blk(ic)
!DIR$ IVDEP
          DO jk = nshift+1, nlev
#else
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(je, jb, vn_tend)
        DO jk = nshift+1, nlev
!$NEC ivdep
          DO ic = 1, p_metrics%nudge_e_dim
            je = p_metrics%nudge_e_idx(ic)
            jb = p_metrics%nudge_e_blk(ic)
#endif
            vn_tend = wfac_old*p_latbc_old%vn(je,jk,jb) + wfac_new*p_latbc_new%vn(je,jk,jb) - p_prog%vn(je,jk,jb)

            ! using a weaker nudging coefficient for vn than for thermodynamic variables turned out to have a
            ! beneficial impact on forecast quality
            p_prog%vn(je,jk,jb) = p_prog%vn(je,jk,jb) + 0.5_wp*p_int%nudgecoeff_e(je,jb)*vn_tend
          ENDDO
        ENDDO
        !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL
      !$ACC END DATA

    ELSE
      CALL finish(TRIM(routine), 'missing arguments')
    ENDIF

  END SUBROUTINE limarea_nudging_latbdy




  !>
  !! This routine executes UPPER boundary nudging for the limited-area mode.
  !!
  SUBROUTINE limarea_nudging_upbdy (p_patch, p_prog, ptr_tracer, p_metrics, p_diag, &
                                  p_int, p_latbc_const, p_latbc_old, p_latbc_new, lc1, lc2)

    TYPE(t_patch),      INTENT(IN)    :: p_patch
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog
    REAL(wp), CONTIGUOUS, INTENT(inout) :: ptr_tracer(:,:,:,:)
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag),    INTENT(INOUT) :: p_diag
    TYPE(t_int_state),  INTENT(IN)    :: p_int

    ! alternative input data, either for constant or time-dependent lateral boundary conditions
    TYPE(t_nh_prog),    INTENT(IN), OPTIONAL :: p_latbc_const
    TYPE(t_pi_atm),     INTENT(IN), OPTIONAL :: p_latbc_old, p_latbc_new
    REAL(wp),           INTENT(IN), OPTIONAL :: lc1     ! time interpolation weight
    REAL(wp),           INTENT(IN), OPTIONAL :: lc2     ! time interpolation weight

    ! local variables
    INTEGER  :: jg
    INTEGER  :: jb, jc, jk, je
    INTEGER  :: ke_nudge, i_startblk, i_endblk, i_startidx, i_endidx
    REAL(wp) :: wfac_old, wfac_new, pres, temp, qv, tempv_inc, pres_inc
    REAL(wp) :: rho_tend, thv_tend, vn_tend
    REAL(wp) :: rd_o_cvd, rd_o_p0ref, nudgecoeff
    REAL(wp) :: max_nudge_coeff_vn, max_nudge_coeff_thermdyn
    LOGICAL  :: bdymask(nproma)
    LOGICAL  :: lnudge_hydro_pres_ubn
    !
    CHARACTER(len=*), PARAMETER :: routine = 'limarea_nudging_upbdy'

    jg = p_patch%id

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    ke_nudge = nudging_config(jg)%ilev_end

    IF (PRESENT(lc1) .AND. PRESENT(lc2)) THEN
      wfac_old = lc1
      wfac_new = lc2  
    ELSE
      wfac_old = -999._wp
      wfac_new = -999._wp
   ENDIF
   
    !
    ! Check if hydrostatic or nonhydrostatic thermodynamic variables shall be used for computing nudging increments 
    lnudge_hydro_pres_ubn = nudging_config(jg)%thermdyn_type == ithermdyn_type%hydrostatic .AND. ltransport
    !
    ! Max. nudging coefficients (qv is not nudged in upper boundary zone)
    max_nudge_coeff_vn       = nudging_config(jg)%max_nudge_coeff_vn
    max_nudge_coeff_thermdyn = nudging_config(jg)%max_nudge_coeff_thermdyn

#ifdef _OPENACC
    CALL finish(TRIM(routine), 'UPPER boundary nudging for the limited-area mode is not tested on GPU')
#endif

    IF (PRESENT(p_latbc_const) .AND. (PRESENT(p_latbc_old) .OR. PRESENT(p_latbc_new))) THEN

      CALL finish(TRIM(routine), 'conflicting arguments')


    ELSE IF (PRESENT(p_latbc_const)) THEN ! Mode for constant lateral boundary data

      CALL finish(TRIM(routine), 'upper boundary nudging not implemented for constant forcing data')


    ELSE IF (PRESENT(p_latbc_old) .AND. PRESENT(p_latbc_new)) THEN ! Mode for time-dependent lateral boundary data


      !$ACC DATA PRESENT(p_diag%temp, p_diag%pres, p_diag%tempv, p_prog%rho) &
      !$ACC   PRESENT(ptr_tracer, p_prog%theta_v, p_prog%exner, p_prog%vn) &
      !$ACC   PRESENT(p_patch%cells%refin_ctrl, p_patch%edges%refin_ctrl) &
      !$ACC   CREATE(bdymask) &
      !$ACC   PRESENT(p_latbc_old, p_latbc_new) &
      !$ACC   PRESENT(p_metrics%nudgecoeff_vert) &
      !$ACC   PRESENT(p_int, p_int%nudgecoeff_c) &
      !$ACC   PRESENT(p_latbc_old%pres, p_latbc_old%temp, p_latbc_old%qv) &
      !$ACC   PRESENT(p_latbc_new%pres, p_latbc_new%temp, p_latbc_new%qv) &
      !$ACC   PRESENT(p_latbc_old%theta_v, p_latbc_old%rho, p_latbc_old%vn) &
      !$ACC   PRESENT(p_latbc_new%theta_v, p_latbc_new%rho, p_latbc_new%vn) &
      !$ACC   PRESENT(p_prog, p_diag, p_metrics, p_patch)

      IF (.NOT. PRESENT(lc1) .OR. .NOT. PRESENT(lc2)) THEN
        CALL finish(routine,'missing interpolation weights wfac_old, wfac_new')
      ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      i_startblk = p_patch%cells%start_block(grf_bdywidth_c+1)
      i_endblk   = p_patch%cells%end_block(min_rlcell)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,pres,temp,qv,tempv_inc,pres_inc,rho_tend,thv_tend,&
!$OMP            nudgecoeff,bdymask) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

        ! Exclude halo points of boundary interpolation zone (causes sync error otherwise)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          bdymask(jc) = p_patch%cells%refin_ctrl(jc,jb)>=1 .AND. p_patch%cells%refin_ctrl(jc,jb)<=grf_bdywidth_c
        ENDDO
        !$ACC END PARALLEL

        IF (lnudge_hydro_pres_ubn) THEN
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(nudgecoeff, pres, temp, qv, tempv_inc, pres_inc, thv_tend, rho_tend)
          DO jk = 1, ke_nudge
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              nudgecoeff = MERGE(0._wp, MAX(MERGE(p_int%nudgecoeff_c(jc,jb),0._wp,jg==1),  &
                &          max_nudge_coeff_thermdyn*p_metrics%nudgecoeff_vert(jk)), bdymask(jc))

              pres = wfac_old*p_latbc_old%pres(jc,jk,jb) + wfac_new*p_latbc_new%pres(jc,jk,jb)
              temp = wfac_old*p_latbc_old%temp(jc,jk,jb) + wfac_new*p_latbc_new%temp(jc,jk,jb)
              qv   = wfac_old*p_latbc_old%qv(jc,jk,jb)   + wfac_new*p_latbc_new%qv(jc,jk,jb)

              tempv_inc = (temp-p_diag%temp(jc,jk,jb))*(1._wp+vtmpc1*qv) + &
                 (qv-ptr_tracer(jc,jk,jb,iqv))*vtmpc1*temp
              pres_inc  = pres-p_diag%pres(jc,jk,jb)

              thv_tend = tempv_inc/p_prog%exner(jc,jk,jb) - rd_o_cpd*p_prog%theta_v(jc,jk,jb)/pres*pres_inc
              rho_tend = ( pres_inc/p_diag%tempv(jc,jk,jb) -                 &
                tempv_inc*p_diag%pres(jc,jk,jb)/p_diag%tempv(jc,jk,jb)**2 )/rd

              p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + nudgecoeff*rho_tend
              p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + nudgecoeff*thv_tend
              p_prog%exner(jc,jk,jb)   = MERGE(p_prog%exner(jc,jk,jb), &
                EXP(rd_o_cvd*LOG(rd_o_p0ref*p_prog%rho(jc,jk,jb)*p_prog%theta_v(jc,jk,jb))),bdymask(jc))

            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ELSE
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(nudgecoeff, thv_tend, rho_tend)
          DO jk = 1, ke_nudge
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              nudgecoeff = MERGE(0._wp, MAX(p_int%nudgecoeff_c(jc,jb),  &
                &          max_nudge_coeff_thermdyn*p_metrics%nudgecoeff_vert(jk)), bdymask(jc))

              thv_tend = wfac_old*p_latbc_old%theta_v(jc,jk,jb) + wfac_new*p_latbc_new%theta_v(jc,jk,jb) - p_prog%theta_v(jc,jk,jb)
              rho_tend = wfac_old*p_latbc_old%rho(jc,jk,jb) + wfac_new*p_latbc_new%rho(jc,jk,jb)- p_prog%rho(jc,jk,jb)

              p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + nudgecoeff*rho_tend
              p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + nudgecoeff*thv_tend
              p_prog%exner(jc,jk,jb)   = MERGE(p_prog%exner(jc,jk,jb), &
                EXP(rd_o_cvd*LOG(rd_o_p0ref*p_prog%rho(jc,jk,jb)*p_prog%theta_v(jc,jk,jb))),bdymask(jc))
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

      ENDDO
!$OMP END DO

      i_startblk = p_patch%edges%start_block(grf_bdywidth_e+1)
      i_endblk   = p_patch%edges%end_block(min_rledge)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk,vn_tend,nudgecoeff,bdymask) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, grf_bdywidth_e+1, min_rledge)

        ! Exclude halo points of boundary interpolation zone (causes sync error otherwise)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO je = i_startidx, i_endidx
          bdymask(je) = p_patch%edges%refin_ctrl(je,jb)>=1 .AND. p_patch%edges%refin_ctrl(je,jb)<=grf_bdywidth_e
        ENDDO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(nudgecoeff, vn_tend)
        DO jk = 1, ke_nudge
!DIR$ IVDEP
          DO je = i_startidx, i_endidx

            nudgecoeff = MERGE(0._wp, MAX(0.5_wp*p_int%nudgecoeff_e(je,jb),  &
              &          max_nudge_coeff_vn*p_metrics%nudgecoeff_vert(jk)), bdymask(je))

            vn_tend = wfac_old*p_latbc_old%vn(je,jk,jb) + wfac_new*p_latbc_new%vn(je,jk,jb) - p_prog%vn(je,jk,jb)

            p_prog%vn(je,jk,jb) = p_prog%vn(je,jk,jb) + nudgecoeff*vn_tend
          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ENDDO
      !$ACC WAIT(1)
      !$ACC END DATA
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    ELSE
      CALL finish(TRIM(routine), 'missing arguments')
    ENDIF

  END SUBROUTINE limarea_nudging_upbdy


  !>
  !! This routine interpolates lateral boundary data (or more generally forcing data) 
  !! from the base domain to a specific child domain (jg), and all childs contained therein.
  !!
  RECURSIVE SUBROUTINE intp_nestubc_nudging (p_patch, latbc_data, jg)

    TYPE(t_patch),       INTENT(INOUT) :: p_patch(:)
    TYPE(t_latbc_state), INTENT(INOUT) :: latbc_data        ! source data (for base domain)

    INTEGER, INTENT(IN) :: jg      ! domain ID of the target (child) domain. 
                                   ! i.e. the domain to which the data are interpolated

    INTEGER :: nlev, i_chidx, jn

    TYPE(t_patch),         POINTER  :: p_plp => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
    TYPE(t_int_state),     POINTER  :: p_int => NULL()

    REAL(wp), ALLOCATABLE :: pres_lp(:,:,:), temp_lp(:,:,:), qv_lp(:,:,:), vn_lp(:,:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: vDateTime_str_cur
    !
    CHARACTER(len=*), PARAMETER :: routine = 'intp_nestubc_nudging'


    ! only interpolate, if upper boundary nudging is activated for domain jg
    !
    IF (nudging_config(jg)%nudge_type==indg_type%ubn) THEN

      IF (msg_level >= 12) THEN
        CALL datetimeToString(latbc_data%vDateTime, vDateTime_str_cur)
        WRITE(message_text,'(a,i2,a,a)') 'latbc data interpolation to DOM', jg, &
          &                               ' at ', TRIM(vDateTime_str_cur) 
        CALL message(routine, message_text)
      ENDIF

      p_plp => p_patch_local_parent(jg)
      p_grf => p_grf_state_local_parent(jg)
      p_int => p_int_state_local_parent(jg)

      nlev = p_plp%nlev

      ! domain index of domain jg as seen from the parent domain
      i_chidx = p_patch(jg)%parent_child_index


      ! temporary local_parent arrays
      ALLOCATE(pres_lp(nproma, nlev, p_plp%nblks_c),  &
        &      temp_lp(nproma, nlev, p_plp%nblks_c),  & 
        &      qv_lp  (nproma, nlev, p_plp%nblks_c),  &
        &      vn_lp  (nproma, nlev, p_plp%nblks_e)   )


      CALL exchange_data_mult(p_plp%comm_pat_glb_to_loc_c, 3, 3*nlev,     &
        &                     RECV1=pres_lp, SEND1=latbc_data%atm%pres,   &
        &                     RECV2=temp_lp, SEND2=latbc_data%atm%temp,   &
        &                     RECV3=qv_lp,   SEND3=latbc_data%atm%qv )

      CALL exchange_data(p_plp%comm_pat_glb_to_loc_e, &
        &                RECV=vn_lp, SEND=latbc_data%atm%vn)


      ! parent to child interpolation
      !
      ! pres, temp, qv
      !
      CALL exchange_data_mult(p_plp%comm_pat_c, 3, 3*nlev,             &
        &                     recv1=pres_lp, recv2=temp_lp, recv3=qv_lp)
      !
      CALL interpol_scal_nudging (p_plp, p_int, p_grf%p_dom(i_chidx), 0, 3, 1,           &
        &                         f3din1=pres_lp, f3dout1=latbc_data%atm_child(jg)%pres, &
        &                         f3din2=temp_lp, f3dout2=latbc_data%atm_child(jg)%temp, &
        &                         f3din3=qv_lp,   f3dout3=latbc_data%atm_child(jg)%qv    )
      !
      CALL sync_patch_array_mult(SYNC_C, p_patch(jg), 3,        &
        &                        latbc_data%atm_child(jg)%pres, &
        &                        latbc_data%atm_child(jg)%temp, &
        &                        latbc_data%atm_child(jg)%qv)

      ! vn
      !
      CALL exchange_data(p_plp%comm_pat_e, vn_lp)
      !
      CALL interpol_vec_nudging (p_plp, p_patch(jg), p_int, p_grf%p_dom(i_chidx), &
        &                        0, 1, vn_lp, latbc_data%atm_child(jg)%vn )
      !
      CALL sync_patch_array(SYNC_E, p_patch(jg), latbc_data%atm_child(jg)%vn)

      DEALLOCATE(pres_lp, temp_lp, qv_lp, vn_lp)

    ENDIF


    ! in case that child domains exist for the current domain jg, repeat the interpolation 
    ! for each child domain.
    DO jn = 1, p_patch(jg)%n_childdom

       CALL intp_nestubc_nudging (p_patch    = p_patch(:),    &
        &                        latbc_data  = latbc_data,    &
        &                        jg          = p_patch(jg)%child_id(jn))
    ENDDO


  END SUBROUTINE intp_nestubc_nudging



  !>
  !! This routine executes boundary nudging for one-way nested domains
  !!
  !!
  SUBROUTINE nest_boundary_nudging(jg, nnew, nnew_rcf)


    INTEGER, INTENT(IN)  :: jg, nnew, nnew_rcf

    ! Pointers
    TYPE(t_nh_state),  POINTER ::  p_nh
    TYPE(t_int_state), POINTER ::  p_int

    ! Indices
    INTEGER :: jb, jc, je, jk, jt, ic

    INTEGER :: nlev          ! number of vertical full levels
    INTEGER :: ntracer_nudge !< number of tracers to be nudged

    REAL(wp) :: rd_o_cvd, rd_o_p0ref, upper_lim, lower_lim

    ! Set pointers
    p_nh  => p_nh_state(jg)
    p_int => p_int_state(jg)

    IF (iforcing > 1) THEN  ! tracers represent moisture variables
      ntracer_nudge = 3     ! take only QV, QC and QI - nudging precip variables does not make much sense
    ELSE
      ntracer_nudge = ntracer
    ENDIF

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    ! number of vertical levels
    nlev = p_patch(jg)%nlev

!$OMP PARALLEL

!$OMP DO PRIVATE(jk,jc,jb,ic,upper_lim,lower_lim) ICON_OMP_DEFAULT_SCHEDULE
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_nh%metrics%nudge_c_dim
      jc = p_nh%metrics%nudge_c_idx(ic)
      jb = p_nh%metrics%nudge_c_blk(ic)
!DIR$ IVDEP
      DO jk = 1, nlev
#else
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(jc, jb, lower_lim, upper_lim)
    DO jk = 1, nlev
!$NEC ivdep
      DO ic = 1, p_nh%metrics%nudge_c_dim
        jc = p_nh%metrics%nudge_c_idx(ic)
        jb = p_nh%metrics%nudge_c_blk(ic)
#endif
        upper_lim = 1.0025_wp*p_nh%prog(nnew)%rho(jc,jk,jb)
        lower_lim = 0.9975_wp*p_nh%prog(nnew)%rho(jc,jk,jb)
        p_nh%prog(nnew)%rho(jc,jk,jb) =                               &
          p_nh%prog(nnew)%rho(jc,jk,jb) + p_int%nudgecoeff_c(jc,jb)*  &
          p_nh%diag%grf_tend_rho(jc,jk,jb)
        p_nh%prog(nnew)%rho(jc,jk,jb) = MAX(lower_lim,MIN(upper_lim,p_nh%prog(nnew)%rho(jc,jk,jb)))

        upper_lim = 1.0025_wp*p_nh%prog(nnew)%theta_v(jc,jk,jb)
        lower_lim = 0.9975_wp*p_nh%prog(nnew)%theta_v(jc,jk,jb)
        p_nh%prog(nnew)%theta_v(jc,jk,jb) =                              &
          p_nh%prog(nnew)%theta_v(jc,jk,jb) + p_int%nudgecoeff_c(jc,jb)* &
          p_nh%diag%grf_tend_thv(jc,jk,jb)
        p_nh%prog(nnew)%theta_v(jc,jk,jb) = MAX(lower_lim,MIN(upper_lim,p_nh%prog(nnew)%theta_v(jc,jk,jb)))

        p_nh%prog(nnew)%exner(jc,jk,jb) =                                  &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)))

        p_nh%prog(nnew)%w(jc,jk,jb) =                               &
          p_nh%prog(nnew)%w(jc,jk,jb) + p_int%nudgecoeff_c(jc,jb)*  &
          p_nh%diag%grf_tend_w(jc,jk,jb)

      ENDDO
    ENDDO
    !$ACC END PARALLEL
!$OMP END DO

!$OMP DO PRIVATE(jb,jk,je,ic) ICON_OMP_DEFAULT_SCHEDULE
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_nh%metrics%nudge_e_dim
      je = p_nh%metrics%nudge_e_idx(ic)
      jb = p_nh%metrics%nudge_e_blk(ic)
!DIR$ IVDEP
      DO jk = 1, nlev
#else
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(je, jb)
    DO jk = 1, nlev
!$NEC ivdep
      DO ic = 1, p_nh%metrics%nudge_e_dim
        je = p_nh%metrics%nudge_e_idx(ic)
        jb = p_nh%metrics%nudge_e_blk(ic)
#endif
        p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)        &
          + p_int%nudgecoeff_e(je,jb)*p_nh%diag%grf_tend_vn(je,jk,jb)
      ENDDO
    ENDDO
    !$ACC END PARALLEL
!$OMP END DO

    IF (ltransport) THEN
!$OMP DO PRIVATE(jk,jc,jb,jt,ic,upper_lim,lower_lim) ICON_OMP_DEFAULT_SCHEDULE
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
#ifdef __LOOP_EXCHANGE
      DO ic = 1, p_nh%metrics%nudge_c_dim
        jc = p_nh%metrics%nudge_c_idx(ic)
        jb = p_nh%metrics%nudge_c_blk(ic)
        DO jt = 1, ntracer_nudge
!DIR$ IVDEP
          DO jk = 1, nlev
#else
      !$ACC LOOP GANG VECTOR COLLAPSE(3) PRIVATE(jc, jb, lower_lim, upper_lim)
      DO jt = 1, ntracer_nudge
        DO jk = 1, nlev
!$NEC ivdep
          DO ic = 1, p_nh%metrics%nudge_c_dim
            jc = p_nh%metrics%nudge_c_idx(ic)
            jb = p_nh%metrics%nudge_c_blk(ic)
#endif

            upper_lim = 1.01_wp*p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt)
            lower_lim = 0.99_wp*p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt)
            p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt) = p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt) + &
              p_int%nudgecoeff_c(jc,jb)*p_nh%diag%grf_tend_tracer(jc,jk,jb,jt)
            p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt) =  MAX(lower_lim,MIN(upper_lim,&
              p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt)))

          ENDDO
        ENDDO
      ENDDO
      !$ACC END PARALLEL
!$OMP END DO NOWAIT
    ENDIF

!$OMP END PARALLEL


  END SUBROUTINE nest_boundary_nudging

  !>
  !! This routine executes boundary nudging for density (for use with 2-way nesting)
  !!
  !!
  SUBROUTINE density_boundary_nudging(jg, nnew)


    INTEGER, INTENT(IN)  :: jg, nnew

    ! Pointers
    TYPE(t_nh_state),  POINTER ::  p_nh
    TYPE(t_int_state), POINTER ::  p_int

    ! Indices
    INTEGER :: jb, jc, jk, ic

    INTEGER :: nlev  ! number of vertical full levels

    REAL(wp) :: rd_o_cvd, rd_o_p0ref

    ! Set pointers
    p_nh  => p_nh_state(jg)
    p_int => p_int_state(jg)

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    ! number of vertical levels
    nlev = p_patch(jg)%nlev

!$OMP PARALLEL

!$OMP DO PRIVATE(jk,jc,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_nh%metrics%nudge_c_dim
      jc = p_nh%metrics%nudge_c_idx(ic)
      jb = p_nh%metrics%nudge_c_blk(ic)
!DIR$ IVDEP
      DO jk = 1, nlev
#else
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR PRIVATE(jc, jb) COLLAPSE(2)
    DO jk = 1, nlev
!$NEC ivdep
      DO ic = 1, p_nh%metrics%nudge_c_dim
        jc = p_nh%metrics%nudge_c_idx(ic)
        jb = p_nh%metrics%nudge_c_blk(ic)
#endif
        p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnew)%rho(jc,jk,jb) +  &
          MIN(0.333_wp,3._wp*p_int%nudgecoeff_c(jc,jb))*                 &
          p_nh%diag%grf_tend_rho(jc,jk,jb)

        p_nh%prog(nnew)%exner(jc,jk,jb) =                                  &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)))

      ENDDO
    ENDDO
    !$ACC END PARALLEL
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  END SUBROUTINE density_boundary_nudging


END MODULE mo_nh_nest_utilities

