!
! mo_nh_diffusion
!
! Diffusion in the nonhydrostatic model
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

MODULE mo_nh_diffusion

  USE mo_kind,                ONLY: wp, vp
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: l_limited_area, lfeedback
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_cell
  USE mo_interpol_config,     ONLY: nudge_max_coeff
  USE mo_intp,                ONLY: edges2cells_vector, cells2verts_scalar
  USE mo_nonhydrostatic_config, ONLY: l_zdiffu_t, ndyn_substeps
  USE mo_diffusion_config,    ONLY: diffusion_config
  USE mo_turbdiff_config,     ONLY: turbdiff_config
  USE mo_parallel_config,     ONLY: nproma, cpu_min_nproma
  USE mo_run_config,          ONLY: ltimer, iforcing, lvert_nest
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c
  USE mo_impl_constants    ,  ONLY: min_rledge_int, min_rlcell_int, min_rlvert_int, inwp, iaes
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_e, grf_bdywidth_c
  USE mo_math_laplace,        ONLY: nabla4_vec
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_gridref_config,      ONLY: denom_diffu_v
  USE mo_parallel_config,     ONLY: p_test_run
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, &
                                    sync_patch_array_mult, sync_patch_array_mult_mp
  USE mo_physical_constants,  ONLY: cvd_o_rd, grav
  USE mo_timer,               ONLY: timer_nh_hdiffusion, timer_start, timer_stop
  USE mo_fortran_tools,       ONLY: init
  USE mo_vertical_grid,       ONLY: nrdmax
  USE mo_mpi,                 ONLY: i_am_accel_node

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: diffusion

  ! On the vectorizing DWD-NEC the diagnostics for the tendencies of the normal wind
  ! from terms xyz, ddt_vn_xyz, is disabled by default due to the fear that the
  ! conditional storage in conditionally allocated global fields is attempted even if
  ! the condition is not given and therefore the global field not allocated. If this
  ! happens, this would results in a corrupted memory.
  ! (Requested by G. Zaengl based on earlier problems with similar constructs.)

  CONTAINS

  !>
  !! diffusion
  !!
  !! Computes the horizontal diffusion of velocity and temperature
  !!
  SUBROUTINE  diffusion(p_nh_prog,p_nh_diag,p_nh_metrics,p_patch,p_int,dtime,linit)

    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch    !< single patch
    TYPE(t_int_state),INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_nh_prog), INTENT(inout)    :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_diag), INTENT(inout)    :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    REAL(wp), INTENT(in)            :: dtime      !< time step
    LOGICAL,  INTENT(in)            :: linit      !< initial call or runtime call

    ! local variables - vp means variable precision depending on the __MIXED_PRECISION cpp flag
    REAL(vp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: z_temp
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_nabla2_e
    REAL(vp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: z_nabla2_c
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_nabla4_e
    REAL(vp), DIMENSION(nproma,p_patch%nlev) :: z_nabla4_e2

    REAL(wp):: diff_multfac_vn(p_patch%nlev), diff_multfac_w, diff_multfac_n2w(p_patch%nlev)
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc, je, ic, ishift, nshift, jk1
    INTEGER :: nlev, nlevp1              !< number of full and half levels

    ! start index levels and diffusion coefficient for boundary diffusion
    INTEGER :: start_bdydiff_e
    REAL(wp):: fac_bdydiff_v

    ! For Smagorinsky diffusion - vp means variable precision depending on the __MIXED_PRECISION cpp flag
    REAL(vp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: kh_smag_e
    REAL(vp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: kh_smag_ec
    REAL(vp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: u_vert
    REAL(vp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: v_vert
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: u_cell
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: v_cell
    REAL(vp), DIMENSION(nproma,p_patch%nlev) :: kh_c, div

    REAL(vp) :: dvt_norm, dvt_tang, vn_vert1, vn_vert2, vn_vert3, vn_vert4, vn_cell1, vn_cell2

    REAL(vp) :: smag_offset, nabv_tang, nabv_norm, rd_o_cvd, nudgezone_diff, bdy_diff, enh_diffu
    REAL(vp), DIMENSION(p_patch%nlev) :: smag_limit, diff_multfac_smag, enh_smag_fac, smag_blending
    INTEGER  :: nblks_zdiffu, nproma_zdiffu, npromz_zdiffu, nlen_zdiffu

    REAL(wp) :: alin, dz32, df32, dz42, df42, bqdr, aqdr, zf, dzlin, dzqdr

    ! Additional variables for 3D Smagorinsky coefficient
    REAL(wp):: z_w_v(nproma,p_patch%nlevp1,p_patch%nblks_v)
    REAL(wp), DIMENSION(nproma,p_patch%nlevp1) :: z_vn_ie, z_vt_ie
    REAL(wp) :: dvndz, dvtdz, dwdz, dthvdz, dwdn, dwdt, def_3d, def_2d, def_aux

    ! Variables for provisional fix against runaway cooling in local topography depressions
    INTEGER  :: icount(p_patch%nblks_c), iclist(2*nproma,p_patch%nblks_c), iklist(2*nproma,p_patch%nblks_c)
    REAL(wp) :: tdlist(2*nproma,p_patch%nblks_c), tdiff, trefdiff, thresh_tdiff, z_theta, fac2d

    INTEGER,  POINTER, CONTIGUOUS :: icidx(:,:,:), icblk(:,:,:), ieidx(:,:,:), &
      &                              ieblk(:,:,:), ividx(:,:,:), ivblk(:,:,:), &
      &                              iecidx(:,:,:), iecblk(:,:,:)
    INTEGER,  POINTER, CONTIGUOUS :: icell(:,:), ilev(:,:), iblk(:,:)
    REAL(wp), POINTER, CONTIGUOUS :: vcoef(:,:), geofac_n2s(:,:)

    LOGICAL :: ltemp_diffu, lfeedback_jg
    INTEGER :: diffu_type, discr_vn, discr_t
    INTEGER :: jg                 !< patch ID







    ! Variables for tendency diagnostics
    REAL(wp) :: z_d_vn_hdf
    REAL(wp) :: r_dtime

    !--------------------------------------------------------------------------

    ! get patch ID
    jg = p_patch%id

    ! prepare for tendency diagnostics
    r_dtime = 1._wp/dtime

    start_bdydiff_e = 5 ! refin_ctrl level at which boundary diffusion starts

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    nshift = p_patch%nshift_total

    ! Normalized diffusion coefficient for boundary diffusion
    fac_bdydiff_v = SQRT(REAL(ndyn_substeps,wp))/denom_diffu_v

    ! scaling factor for enhanced diffusion in nudging zone (if present, i.e. for
    ! limited-area runs and one-way nesting)
    nudgezone_diff = 0.04_wp/(nudge_max_coeff + dbl_eps)

    ! scaling factor for enhanced near-boundary diffusion for
    ! two-way nesting (used with Smagorinsky diffusion only; not needed otherwise)
    bdy_diff = 0.015_wp/(nudge_max_coeff + dbl_eps)

    ! threshold temperature deviation from neighboring grid points
    ! that activates extra diffusion against runaway cooling
    thresh_tdiff = - 5._wp

    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    icidx => p_patch%cells%neighbor_idx
    icblk => p_patch%cells%neighbor_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    rd_o_cvd = 1._wp/cvd_o_rd

    diffu_type  = diffusion_config(jg)%hdiff_order
    discr_vn    = diffusion_config(jg)%itype_vn_diffu
    discr_t     = diffusion_config(jg)%itype_t_diffu

    IF (linit) THEN ! enhanced diffusion at all levels for initial velocity filtering call
      diff_multfac_vn(:) = diffusion_config(jg)%k4/3._wp*diffusion_config(jg)%hdiff_efdt_ratio
      smag_offset        =  0.0_vp
      diffu_type = 5 ! always combine nabla4 background diffusion with Smagorinsky diffusion for initial filtering call
      smag_limit(:) = 0.125_wp-4._wp*diff_multfac_vn(:)
    ELSE ! combination with divergence damping inside the dynamical core
      IF (diffu_type == 4) THEN
        diff_multfac_vn(:) = MIN(1._wp/128._wp,diffusion_config(jg)%k4*REAL(ndyn_substeps,wp)/ &
                                 3._wp*p_nh_metrics%enhfac_diffu(:))
      ELSE ! For Smagorinsky diffusion, the Smagorinsky coefficient rather than the background
           ! diffusion coefficient is enhanced near the model top (see below)
        diff_multfac_vn(:) = MIN(1._wp/128._wp,diffusion_config(jg)%k4*REAL(ndyn_substeps,wp)/3._wp)
      ENDIF

      smag_offset   = REAL(0.25_wp*diffusion_config(jg)%k4*REAL(ndyn_substeps,wp),vp)
      smag_limit(:) = 0.125_wp-4._wp*diff_multfac_vn(:)
    ENDIF

    ! Multiplication factor for nabla4 diffusion on vertical wind speed
    diff_multfac_w = MIN(1._wp/48._wp,diffusion_config(jg)%k4w*REAL(ndyn_substeps,wp))

    ! Factor for additional nabla2 diffusion in upper damping zone
    diff_multfac_n2w(:) = 0._wp
    IF (nrdmax(jg) > 1) THEN ! seems to be redundant, but the NEC issues invalid operations otherwise
      DO jk = 2, nrdmax(jg)
        jk1 = jk + nshift
        diff_multfac_n2w(jk) = 1._wp/12._wp*((vct_a(jk1)-vct_a(nshift+nrdmax(jg)+1))/ &
                               (vct_a(2)-vct_a(nshift+nrdmax(jg)+1)))**4
      ENDDO
    ENDIF


    IF (diffu_type == 5) THEN

      ! temperature diffusion is used only in combination with Smagorinsky diffusion
      ltemp_diffu = diffusion_config(jg)%lhdiff_temp

      ! The Smagorinsky diffusion factor enh_divdamp_fac is defined as a profile in height z
      ! above sea level with 4 height sections:
      !
      ! enh_smag_fac(z) = hdiff_smag_fac                                                    !                  z <= hdiff_smag_z
      ! enh_smag_fac(z) = hdiff_smag_fac  + (z-hdiff_smag_z )* alin                         ! hdiff_smag_z  <= z <= hdiff_smag_z2
      ! enh_smag_fac(z) = hdiff_smag_fac2 + (z-hdiff_smag_z2)*(aqdr+(z-hdiff_smag_z2)*bqdr) ! hdiff_smag_z2 <= z <= hdiff_smag_z4
      ! enh_smag_fac(z) = hdiff_smag_fac4                                                   ! hdiff_smag_z4 <= z
      !
      alin = (diffusion_config(jg)%hdiff_smag_fac2-diffusion_config(jg)%hdiff_smag_fac)/ &
           & (diffusion_config(jg)%hdiff_smag_z2  -diffusion_config(jg)%hdiff_smag_z)
      !
      df32 = diffusion_config(jg)%hdiff_smag_fac3-diffusion_config(jg)%hdiff_smag_fac2
      df42 = diffusion_config(jg)%hdiff_smag_fac4-diffusion_config(jg)%hdiff_smag_fac2
      !
      dz32 = diffusion_config(jg)%hdiff_smag_z3-diffusion_config(jg)%hdiff_smag_z2
      dz42 = diffusion_config(jg)%hdiff_smag_z4-diffusion_config(jg)%hdiff_smag_z2
      !
      bqdr = (df42*dz32-df32*dz42)/(dz32*dz42*(dz42-dz32))
      aqdr =  df32/dz32-bqdr*dz32
      ! 
      DO jk = 1, nlev
        jk1 = jk + nshift
        !
        zf = 0.5_wp*(vct_a(jk1)+vct_a(jk1+1))
        dzlin = MIN( diffusion_config(jg)%hdiff_smag_z2-diffusion_config(jg)%hdiff_smag_z , &
             &  MAX( 0._wp,                          zf-diffusion_config(jg)%hdiff_smag_z ) )
        dzqdr = MIN( diffusion_config(jg)%hdiff_smag_z4-diffusion_config(jg)%hdiff_smag_z2, &
             &  MAX( 0._wp,                          zf-diffusion_config(jg)%hdiff_smag_z2) )
        !
        enh_smag_fac(jk) = REAL(diffusion_config(jg)%hdiff_smag_fac + dzlin*alin + dzqdr*(aqdr+dzqdr*bqdr),vp)
        !
        IF (zf < 100._wp) THEN
          smag_blending(jk) = 0._wp
        ELSE IF (zf < 500._wp) THEN
          smag_blending(jk) = 0.0025_wp*(zf-100._wp)
        ELSE
          smag_blending(jk) = 1._wp
        ENDIF
      ENDDO

      ! Smagorinsky coefficient is also enhanced in the six model levels beneath a vertical nest interface
      IF ((lvert_nest) .AND. (p_patch%nshift > 0)) THEN
        enh_smag_fac(1) = MAX(0.333_vp, enh_smag_fac(1))
        enh_smag_fac(2) = MAX(0.25_vp, enh_smag_fac(2))
        enh_smag_fac(3) = MAX(0.20_vp, enh_smag_fac(3))
        enh_smag_fac(4) = MAX(0.16_vp, enh_smag_fac(4))
        enh_smag_fac(5) = MAX(0.12_vp, enh_smag_fac(5))
        enh_smag_fac(6) = MAX(0.08_vp, enh_smag_fac(6))
      ENDIF

      ! empirically determined scaling factor
      diff_multfac_smag(:) = enh_smag_fac(:)*REAL(dtime,vp)

    ELSE
      ltemp_diffu = .FALSE.
    ENDIF

    !$ACC DATA CREATE(div, kh_c, kh_smag_e, kh_smag_ec, u_vert, v_vert, u_cell, v_cell, z_w_v, z_temp) &
    !$ACC   CREATE(z_nabla4_e, z_nabla4_e2, z_nabla2_e, z_nabla2_c, enh_diffu_3d, icount) &
    !$ACC   CREATE(z_vn_ie, z_vt_ie) &
    !$ACC   COPYIN(nrdmax, diff_multfac_vn, diff_multfac_n2w, diff_multfac_smag, smag_limit, smag_blending) &
    !$ACC   PRESENT(p_patch, p_int, p_nh_prog, p_nh_diag, p_nh_metrics) &
    !$ACC   PRESENT(ividx, ivblk, iecidx, iecblk, icidx, icblk, ieidx, ieblk) &
    !$ACC   IF(i_am_accel_node)

!!! Following variables may be present in certain situations, but we don't want it to fail in the general case.
!!! Should actually be in a separate data region with correct IF condition.
!!! !$ACC               p_nh_diag%div_ic, p_nh_diag%dwdx, p_nh_diag%dwdy, p_nh_diag%hdef_ic,                     &


    ! The diffusion is an intrinsic part of the NH solver, thus it is added to the timer
    IF (ltimer) CALL timer_start(timer_nh_hdiffusion)

    IF (diffu_type == 4) THEN

      CALL nabla4_vec( p_nh_prog%vn, p_patch, p_int, z_nabla4_e,  &
                       opt_rlstart=7,opt_nabla2=z_nabla2_e )

    ELSE IF ( diffu_type == 5 .AND. discr_vn == 1 .AND. .NOT. diffusion_config(jg)%lsmag_3d) THEN

      ! needs to be always initialized with OpenACC
      IF (p_test_run .OR. i_am_accel_node) THEN
        CALL init(u_vert, opt_acc_async=.TRUE.)
        CALL init(v_vert, opt_acc_async=.TRUE.)
      ENDIF

      !  RBF reconstruction of velocity at vertices
      CALL rbf_vec_interpol_vertex( p_nh_prog%vn, p_patch, p_int,             &
                                    u_vert, v_vert, opt_rlend=min_rlvert_int, &
                                    opt_acc_async=.TRUE. )

      rl_start = start_bdydiff_e
      rl_end   = min_rledge_int - 2

      CALL sync_patch_array_mult(SYNC_V,p_patch,2,u_vert,v_vert,                          &
                                 opt_varname="diffusion: u_vert and v_vert")

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,vn_vert1,vn_vert2,vn_vert3,vn_vert4, &
!$OMP            dvt_norm,dvt_tang), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Computation of wind field deformation

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
!$NEC outerloop_unroll(4)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

            vn_vert1 =  u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                        v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_vert(je,jb,1)%v2

            vn_vert2 =  u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                        v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_vert(je,jb,2)%v2

            dvt_tang =  p_patch%edges%tangent_orientation(je,jb)* (   &
                        u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                        p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                        v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                        p_patch%edges%dual_normal_vert(je,jb,2)%v2 - &
                        (u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                        p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                        v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                        p_patch%edges%dual_normal_vert(je,jb,1)%v2) )

            vn_vert3 =  u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                        p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                        v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                        p_patch%edges%primal_normal_vert(je,jb,3)%v2

            vn_vert4 =  u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                        p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                        v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                        p_patch%edges%primal_normal_vert(je,jb,4)%v2

            dvt_norm =  u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                        p_patch%edges%dual_normal_vert(je,jb,4)%v1 + &
                        v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                        p_patch%edges%dual_normal_vert(je,jb,4)%v2 - &
                        (u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                        p_patch%edges%dual_normal_vert(je,jb,3)%v1 + &
                        v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                        p_patch%edges%dual_normal_vert(je,jb,3)%v2)

            ! Smagorinsky diffusion coefficient
            kh_smag_e(je,jk,jb) = diff_multfac_smag(jk)*SQRT(                           (  &
              (vn_vert4-vn_vert3)*p_patch%edges%inv_vert_vert_length(je,jb)- &
              dvt_tang*p_patch%edges%inv_primal_edge_length(je,jb) )**2 + (         &
              (vn_vert2-vn_vert1)*p_patch%edges%tangent_orientation(je,jb)*   &
              p_patch%edges%inv_primal_edge_length(je,jb) +                                &
              dvt_norm*p_patch%edges%inv_vert_vert_length(je,jb))**2 )

            ! The factor of 4 comes from dividing by twice the "correct" length
            z_nabla2_e(je,jk,jb) = 4._wp * (                                      &
              (vn_vert4 + vn_vert3 - 2._wp*p_nh_prog%vn(je,jk,jb))  &
              *p_patch%edges%inv_vert_vert_length(je,jb)**2 +                     &
              (vn_vert2 + vn_vert1 - 2._wp*p_nh_prog%vn(je,jk,jb))  &
              *p_patch%edges%inv_primal_edge_length(je,jb)**2 )

            kh_smag_ec(je,jk,jb) = kh_smag_e(je,jk,jb)
            ! Subtract part of the fourth-order background diffusion coefficient
            kh_smag_e(je,jk,jb) = MAX(0._vp,kh_smag_e(je,jk,jb) - smag_offset)
            ! Limit diffusion coefficient to the theoretical CFL stability threshold
            kh_smag_e(je,jk,jb) = MIN(kh_smag_e(je,jk,jb),smag_limit(jk))
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

      ENDDO ! block jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ELSE IF ( diffu_type == 5 .AND. discr_vn == 1) THEN ! 3D Smagorinsky diffusion

      IF (p_test_run) THEN
        !$ACC KERNELS PRESENT(u_vert, v_vert, z_w_v) ASYNC(1) IF(i_am_accel_node)
        u_vert = 0._vp
        v_vert = 0._vp
        z_w_v  = 0._wp
        !$ACC END KERNELS
      ENDIF

      !  RBF reconstruction of velocity at vertices
      CALL rbf_vec_interpol_vertex( p_nh_prog%vn, p_patch, p_int,                 &
                                    u_vert, v_vert, opt_rlend=min_rlvert_int ,    &
                                    opt_acc_async=.TRUE. )

      rl_start = start_bdydiff_e
      rl_end   = min_rledge_int - 2

      CALL sync_patch_array_mult(SYNC_V,p_patch,2,u_vert,v_vert,                          &
                                 opt_varname="diffusion: u_vert and v_vert 2")

      CALL cells2verts_scalar(p_nh_prog%w, p_patch, p_int%cells_aw_verts, z_w_v, opt_rlend=min_rlvert_int)
      CALL sync_patch_array(SYNC_V,p_patch,z_w_v,opt_varname="diffusion: z_w_v")
      CALL sync_patch_array(SYNC_C,p_patch,p_nh_diag%theta_v_ic,opt_varname="diffusion: theta_v_ic")

      fac2d = 0.0625_wp ! Factor of the 2D deformation field which is used as minimum of the 3D def field

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,vn_vert1,vn_vert2,vn_vert3,vn_vert4,dvt_norm,dvt_tang, &
!$OMP            z_vn_ie,z_vt_ie,dvndz,dvtdz,dwdz,dthvdz,dwdn,dwdt,def_3d,def_2d,def_aux), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)


        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
        DO jk = 2, nlev
          DO je = i_startidx, i_endidx
            z_vn_ie(je,jk) = p_nh_metrics%wgtfac_e(je,jk,jb)*p_nh_prog%vn(je,jk,jb) +   &
             (1._wp - p_nh_metrics%wgtfac_e(je,jk,jb))*p_nh_prog%vn(je,jk-1,jb)
            z_vt_ie(je,jk) = p_nh_metrics%wgtfac_e(je,jk,jb)*p_nh_diag%vt(je,jk,jb) +   &
             (1._wp - p_nh_metrics%wgtfac_e(je,jk,jb))*p_nh_diag%vt(je,jk-1,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) IF(i_am_accel_node)
        DO je = i_startidx, i_endidx
          z_vn_ie(je,1) =                                            &
            p_nh_metrics%wgtfacq1_e(je,1,jb)*p_nh_prog%vn(je,1,jb) + &
            p_nh_metrics%wgtfacq1_e(je,2,jb)*p_nh_prog%vn(je,2,jb) + &
            p_nh_metrics%wgtfacq1_e(je,3,jb)*p_nh_prog%vn(je,3,jb)
          z_vn_ie(je,nlevp1) =                                           &
            p_nh_metrics%wgtfacq_e(je,1,jb)*p_nh_prog%vn(je,nlev,jb)   + &
            p_nh_metrics%wgtfacq_e(je,2,jb)*p_nh_prog%vn(je,nlev-1,jb) + &
            p_nh_metrics%wgtfacq_e(je,3,jb)*p_nh_prog%vn(je,nlev-2,jb)
          z_vt_ie(je,1) =                                            &
            p_nh_metrics%wgtfacq1_e(je,1,jb)*p_nh_diag%vt(je,1,jb) + &
            p_nh_metrics%wgtfacq1_e(je,2,jb)*p_nh_diag%vt(je,2,jb) + &
            p_nh_metrics%wgtfacq1_e(je,3,jb)*p_nh_diag%vt(je,3,jb)
          z_vt_ie(je,nlevp1) =                                           &
            p_nh_metrics%wgtfacq_e(je,1,jb)*p_nh_diag%vt(je,nlev,jb)   + &
            p_nh_metrics%wgtfacq_e(je,2,jb)*p_nh_diag%vt(je,nlev-1,jb) + &
            p_nh_metrics%wgtfacq_e(je,3,jb)*p_nh_diag%vt(je,nlev-2,jb)
        ENDDO
        !$ACC END PARALLEL LOOP

        ! Computation of wind field deformation

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
!$NEC outerloop_unroll(4)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

            vn_vert1 =  u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                        v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_vert(je,jb,1)%v2

            vn_vert2 =  u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                        v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_vert(je,jb,2)%v2

            dvt_tang =  p_patch%edges%tangent_orientation(je,jb)* (   &
                        u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                        p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                        v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                        p_patch%edges%dual_normal_vert(je,jb,2)%v2 - &
                        (u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                        p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                        v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                        p_patch%edges%dual_normal_vert(je,jb,1)%v2) )

            vn_vert3 =  u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                        p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                        v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                        p_patch%edges%primal_normal_vert(je,jb,3)%v2

            vn_vert4 =  u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                        p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                        v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                        p_patch%edges%primal_normal_vert(je,jb,4)%v2

            dvt_norm =  u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                        p_patch%edges%dual_normal_vert(je,jb,4)%v1 + &
                        v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                        p_patch%edges%dual_normal_vert(je,jb,4)%v2 - &
                        (u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                        p_patch%edges%dual_normal_vert(je,jb,3)%v1 + &
                        v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                        p_patch%edges%dual_normal_vert(je,jb,3)%v2)

            ! The factor of 4 comes from dividing by twice the "correct" length
            z_nabla2_e(je,jk,jb) = 4._wp * (                        &
              (vn_vert4 + vn_vert3 - 2._wp*p_nh_prog%vn(je,jk,jb))  &
              *p_patch%edges%inv_vert_vert_length(je,jb)**2 +       &
              (vn_vert2 + vn_vert1 - 2._wp*p_nh_prog%vn(je,jk,jb))  &
              *p_patch%edges%inv_primal_edge_length(je,jb)**2 )

            dvndz = (z_vn_ie(je,jk) - z_vn_ie(je,jk+1)) / p_nh_metrics%ddqz_z_full_e(je,jk,jb)
            dvtdz = (z_vt_ie(je,jk) - z_vt_ie(je,jk+1)) / p_nh_metrics%ddqz_z_full_e(je,jk,jb)

            dwdz =                                                                             &
              (p_int%c_lin_e(je,1,jb)*(p_nh_prog%w(iecidx(je,jb,1),jk,  iecblk(je,jb,1)) -     &
                                       p_nh_prog%w(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) ) +   &
               p_int%c_lin_e(je,2,jb)*(p_nh_prog%w(iecidx(je,jb,2),jk,  iecblk(je,jb,2)) -     &
                                       p_nh_prog%w(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) ) ) / &
               p_nh_metrics%ddqz_z_full_e(je,jk,jb)

            dthvdz =                                                                                    &
              (p_int%c_lin_e(je,1,jb)*(p_nh_diag%theta_v_ic(iecidx(je,jb,1),jk,  iecblk(je,jb,1)) -     &
                                       p_nh_diag%theta_v_ic(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) ) +   &
               p_int%c_lin_e(je,2,jb)*(p_nh_diag%theta_v_ic(iecidx(je,jb,2),jk,  iecblk(je,jb,2)) -     &
                                       p_nh_diag%theta_v_ic(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) ) ) / &
               p_nh_metrics%ddqz_z_full_e(je,jk,jb)

            dwdn = p_patch%edges%inv_dual_edge_length(je,jb)* (            &
              0.5_wp*(p_nh_prog%w(iecidx(je,jb,1),jk,  iecblk(je,jb,1)) +  &
                      p_nh_prog%w(iecidx(je,jb,1),jk+1,iecblk(je,jb,1))) - &
              0.5_wp*(p_nh_prog%w(iecidx(je,jb,2),jk,  iecblk(je,jb,2)) +  &
                      p_nh_prog%w(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)))   )

            dwdt = p_patch%edges%inv_primal_edge_length(je,jb)*p_patch%edges%tangent_orientation(je,jb) *( &
              0.5_wp*(z_w_v(ividx(je,jb,1),jk,ivblk(je,jb,1))+z_w_v(ividx(je,jb,1),jk+1,ivblk(je,jb,1))) - &
              0.5_wp*(z_w_v(ividx(je,jb,2),jk,ivblk(je,jb,2))+z_w_v(ividx(je,jb,2),jk+1,ivblk(je,jb,2)))   )

            def_aux = ( (vn_vert2-vn_vert1)*p_patch%edges%tangent_orientation(je,jb)*                             &
              p_patch%edges%inv_primal_edge_length(je,jb) + dvt_norm*p_patch%edges%inv_vert_vert_length(je,jb) )**2

            def_3d = 0.5_wp*(                                                                &
              ( (vn_vert4-vn_vert3)*p_patch%edges%inv_vert_vert_length(je,jb) )**2 +         &
              (dvt_tang*p_patch%edges%inv_primal_edge_length(je,jb))**2 + dwdz**2) +         &
              def_aux + (dvndz + dwdn)**2 + (dvtdz + dwdt)**2  - 3._wp*grav * dthvdz / (     &
              p_int%c_lin_e(je,1,jb)*p_nh_prog%theta_v(iecidx(je,jb,1),jk,iecblk(je,jb,1)) + &
              p_int%c_lin_e(je,2,jb)*p_nh_prog%theta_v(iecidx(je,jb,2),jk,iecblk(je,jb,2))   )

            def_2d = ((vn_vert4-vn_vert3)*p_patch%edges%inv_vert_vert_length(je,jb)- &
              dvt_tang*p_patch%edges%inv_primal_edge_length(je,jb) )**2 + def_aux

            ! 2D Smagorinsky diffusion coefficient (needed as input for turbulence)
            kh_smag_ec(je,jk,jb) = diff_multfac_smag(jk)*SQRT(def_2d) 

            ! 3D Smagorinsky diffusion coefficient with blending to 2D coefficient near the surface
            kh_smag_e(je,jk,jb) = diff_multfac_smag(jk)*SQRT( &
              MAX(smag_blending(jk)*def_3d + (1._wp-smag_blending(jk))*def_2d, fac2d*def_2d) )

            ! Subtract part of the fourth-order background diffusion coefficient
            kh_smag_e(je,jk,jb) = MAX(0._vp,kh_smag_e(je,jk,jb) - smag_offset)
            ! Limit diffusion coefficient to the theoretical CFL stability threshold
            kh_smag_e(je,jk,jb) = MIN(kh_smag_e(je,jk,jb),smag_limit(jk))

          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ELSE IF ( diffu_type == 5 .AND. discr_vn >= 2) THEN

      !  RBF reconstruction of velocity at vertices and cells
      CALL rbf_vec_interpol_vertex( p_nh_prog%vn, p_patch, p_int,                  &
                                    u_vert, v_vert, opt_rlend=min_rlvert_int-1,    &
                                    opt_acc_async=.TRUE. )

      ! DA: This wait ideally should be removed
      !$ACC WAIT

      IF (discr_vn == 2) THEN
        CALL rbf_vec_interpol_cell( p_nh_prog%vn, p_patch, p_int, &
                                    u_cell, v_cell, opt_rlend=min_rlcell_int-1 )
      ELSE
        CALL edges2cells_vector( p_nh_prog%vn, p_nh_diag%vt, p_patch, p_int, &
                                 u_cell, v_cell, opt_rlend=min_rlcell_int-1 )
      ENDIF

      IF (p_test_run) THEN
        !$ACC KERNELS ASYNC(1) IF(i_am_accel_node)
        z_nabla2_e = 0._wp
        !$ACC END KERNELS
      ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,rl_start,rl_end)

      rl_start = start_bdydiff_e
      rl_end   = min_rledge_int - 1

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,vn_vert1,vn_vert2,vn_cell1,vn_cell2,&
!$OMP             dvt_norm,dvt_tang), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Computation of wind field deformation

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
!$NEC outerloop_unroll(4)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

            vn_vert1 =        u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                              p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                              v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                              p_patch%edges%primal_normal_vert(je,jb,1)%v2

            vn_vert2 =        u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                              p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                              v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                              p_patch%edges%primal_normal_vert(je,jb,2)%v2

            dvt_tang =        p_patch%edges%tangent_orientation(je,jb)* (   &
                              u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                              p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                              v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                              p_patch%edges%dual_normal_vert(je,jb,2)%v2 - &
                             (u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                              p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                              v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                              p_patch%edges%dual_normal_vert(je,jb,1)%v2) )

            vn_cell1 =        u_cell(iecidx(je,jb,1),jk,iecblk(je,jb,1)) * &
                              p_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                              v_cell(iecidx(je,jb,1),jk,iecblk(je,jb,1)) * &
                              p_patch%edges%primal_normal_cell(je,jb,1)%v2

            vn_cell2 =        u_cell(iecidx(je,jb,2),jk,iecblk(je,jb,2)) * &
                              p_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                              v_cell(iecidx(je,jb,2),jk,iecblk(je,jb,2)) * &
                              p_patch%edges%primal_normal_cell(je,jb,2)%v2

            dvt_norm =        u_cell(iecidx(je,jb,2),jk,iecblk(je,jb,2)) * &
                              p_patch%edges%dual_normal_cell(je,jb,2)%v1 + &
                              v_cell(iecidx(je,jb,2),jk,iecblk(je,jb,2)) * &
                              p_patch%edges%dual_normal_cell(je,jb,2)%v2 - &
                             (u_cell(iecidx(je,jb,1),jk,iecblk(je,jb,1)) * &
                              p_patch%edges%dual_normal_cell(je,jb,1)%v1 + &
                              v_cell(iecidx(je,jb,1),jk,iecblk(je,jb,1)) * &
                              p_patch%edges%dual_normal_cell(je,jb,1)%v2)


            ! Smagorinsky diffusion coefficient
            kh_smag_e(je,jk,jb) = diff_multfac_smag(jk)*SQRT(                           (  &
              (vn_cell2-vn_cell1)*p_patch%edges%inv_dual_edge_length(je,jb)-               &
              dvt_tang*p_patch%edges%inv_primal_edge_length(je,jb) )**2 + (         &
              (vn_vert2-vn_vert1)*p_patch%edges%tangent_orientation(je,jb)*  &
              p_patch%edges%inv_primal_edge_length(je,jb) +                                &
              dvt_norm*p_patch%edges%inv_dual_edge_length(je,jb))**2 )

            ! The factor of 4 comes from dividing by twice the "correct" length
            z_nabla2_e(je,jk,jb) = 4._wp * (                                      &
              (vn_cell2 + vn_cell1 - 2._wp*p_nh_prog%vn(je,jk,jb))                &
              *p_patch%edges%inv_dual_edge_length(je,jb)**2 +                     &
              (vn_vert2 + vn_vert1 - 2._wp*p_nh_prog%vn(je,jk,jb))  &
              *p_patch%edges%inv_primal_edge_length(je,jb)**2 )

          ENDDO
        ENDDO

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

            kh_smag_ec(je,jk,jb) = kh_smag_e(je,jk,jb)
            ! Subtract part of the fourth-order background diffusion coefficient
            kh_smag_e(je,jk,jb) = MAX(0._vp,kh_smag_e(je,jk,jb) - smag_offset)
            ! Limit diffusion coefficient to the theoretical CFL stability threshold
            kh_smag_e(je,jk,jb) = MIN(kh_smag_e(je,jk,jb),smag_limit(jk))
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF

    ! Compute input quantities for turbulence scheme
    IF ( diffu_type == 5 .AND. (turbdiff_config(jg)%itype_sher >= 1 .OR. turbdiff_config(jg)%ltkeshs)) THEN

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx,kh_c,div), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

            kh_c(jc,jk) = (kh_smag_ec(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb) + &
                           kh_smag_ec(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb) + &
                           kh_smag_ec(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb))/ &
                          diff_multfac_smag(jk)

            div(jc,jk) = p_nh_prog%vn(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
                         p_nh_prog%vn(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
                         p_nh_prog%vn(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
        DO jk = 2, nlev ! levels 1 and nlevp1 are unused
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            p_nh_diag%div_ic(jc,jk,jb) = p_nh_metrics%wgtfac_c(jc,jk,jb)*div(jc,jk) + &
              (1._wp-p_nh_metrics%wgtfac_c(jc,jk,jb))*div(jc,jk-1)

            p_nh_diag%hdef_ic(jc,jk,jb) = (p_nh_metrics%wgtfac_c(jc,jk,jb)*kh_c(jc,jk) + &
              (1._wp-p_nh_metrics%wgtfac_c(jc,jk,jb))*kh_c(jc,jk-1))**2
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ENDIF

    ! Add fourth-order background diffusion
    !
    IF (diffu_type == 5 .AND. diffusion_config(jg)%hdiff_efdt_ratio >0._wp) THEN

      IF (discr_vn > 1) THEN
        CALL sync_patch_array(SYNC_E,p_patch,z_nabla2_e,      &
                              opt_varname="diffusion: nabla2_e")
      END IF

      ! Interpolate nabla2(v) to vertices in order to compute nabla2(nabla2(v))

      IF (p_test_run) THEN
        !$ACC KERNELS ASYNC(1) IF(i_am_accel_node)
        u_vert = 0._wp
        v_vert = 0._wp
        !$ACC END KERNELS
      ENDIF

      CALL rbf_vec_interpol_vertex( z_nabla2_e, p_patch, p_int, u_vert, v_vert, &
                                    opt_rlstart=4, opt_rlend=min_rlvert_int,    &
                                    opt_acc_async=.TRUE. )

      rl_start = grf_bdywidth_e+1
      rl_end   = min_rledge_int

      CALL sync_patch_array_mult(SYNC_V,p_patch,2,u_vert,v_vert,                          &
                                 opt_varname="diffusion: u_vert and v_vert 3")

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)
      lfeedback_jg = lfeedback(jg) ! for OpenACC

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,nabv_tang,nabv_norm,z_nabla4_e2,z_d_vn_hdf), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

         ! Compute nabla4(v)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) PRIVATE(nabv_tang, nabv_norm)
!$NEC outerloop_unroll(4)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

            nabv_tang = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                        v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_vert(je,jb,1)%v2 + &
                        u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                        v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_vert(je,jb,2)%v2

            nabv_norm = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                        p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                        v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                        p_patch%edges%primal_normal_vert(je,jb,3)%v2 + &
                        u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                        p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                        v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                        p_patch%edges%primal_normal_vert(je,jb,4)%v2

            ! The factor of 4 comes from dividing by twice the "correct" length
            z_nabla4_e2(je,jk) = 4._wp * (                          &
              (nabv_norm - 2._wp*z_nabla2_e(je,jk,jb))              &
              *p_patch%edges%inv_vert_vert_length(je,jb)**2 +       &
              (nabv_tang - 2._wp*z_nabla2_e(je,jk,jb))              &
              *p_patch%edges%inv_primal_edge_length(je,jb)**2 )

          ENDDO
        ENDDO

        ! Apply diffusion for the case of diffu_type = 5
        IF ( jg == 1 .AND. l_limited_area .OR. jg > 1 .AND. .NOT. lfeedback_jg) THEN
          !
          ! Domains with lateral boundary and nests without feedback
          !
          !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) PRIVATE(z_d_vn_hdf)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              !
              z_d_vn_hdf =   p_patch%edges%area_edge(je,jb)                              &
                &          * (  MAX(nudgezone_diff*p_int%nudgecoeff_e(je,jb),            &
                &                   REAL(kh_smag_e(je,jk,jb),wp)) * z_nabla2_e(je,jk,jb) &
                &             - p_patch%edges%area_edge(je,jb)                           &
                &             * diff_multfac_vn(jk) * z_nabla4_e2(je,jk)   )
              !
              p_nh_prog%vn(je,jk,jb)            =  p_nh_prog%vn(je,jk,jb)         + z_d_vn_hdf
              !
              IF (p_nh_diag%ddt_vn_hdf_is_associated) THEN
                p_nh_diag%ddt_vn_hdf(je,jk,jb)  =  p_nh_diag%ddt_vn_hdf(je,jk,jb) + z_d_vn_hdf * r_dtime
              END IF
              !
              IF (p_nh_diag%ddt_vn_dyn_is_associated) THEN
                p_nh_diag%ddt_vn_dyn(je,jk,jb)  =  p_nh_diag%ddt_vn_dyn(je,jk,jb) + z_d_vn_hdf * r_dtime
              END IF
              !
            ENDDO
          ENDDO

        ELSE IF (jg > 1) THEN
          !
          ! Nests with feedback
          !
          !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) PRIVATE(z_d_vn_hdf)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              !
              z_d_vn_hdf =   p_patch%edges%area_edge(je,jb)                                                       &
                &          * (  kh_smag_e(je,jk,jb) * z_nabla2_e(je,jk,jb)                                        &
                &             - p_patch%edges%area_edge(je,jb)                                                    &
                &             * MAX(diff_multfac_vn(jk),bdy_diff*p_int%nudgecoeff_e(je,jb)) * z_nabla4_e2(je,jk) )
              !
              p_nh_prog%vn(je,jk,jb)            =  p_nh_prog%vn(je,jk,jb)         + z_d_vn_hdf
              !
              IF (p_nh_diag%ddt_vn_hdf_is_associated) THEN
                p_nh_diag%ddt_vn_hdf(je,jk,jb)  =  p_nh_diag%ddt_vn_hdf(je,jk,jb) + z_d_vn_hdf * r_dtime
              END IF
              !
              IF (p_nh_diag%ddt_vn_dyn_is_associated) THEN
                p_nh_diag%ddt_vn_dyn(je,jk,jb)  =  p_nh_diag%ddt_vn_dyn(je,jk,jb) + z_d_vn_hdf * r_dtime
              END IF
              !
            ENDDO
          ENDDO

        ELSE
          !
          ! Global domains
          !
          !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) PRIVATE(z_d_vn_hdf)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              !
              z_d_vn_hdf =   p_patch%edges%area_edge(je,jb)                 &
                &          * (  kh_smag_e(je,jk,jb) * z_nabla2_e(je,jk,jb)  &
                &             - p_patch%edges%area_edge(je,jb)              &
                &             * diff_multfac_vn(jk) * z_nabla4_e2(je,jk)   )
              !
              p_nh_prog%vn(je,jk,jb)            =  p_nh_prog%vn(je,jk,jb)         + z_d_vn_hdf
              !
              IF (p_nh_diag%ddt_vn_hdf_is_associated) THEN
                p_nh_diag%ddt_vn_hdf(je,jk,jb)  =  p_nh_diag%ddt_vn_hdf(je,jk,jb) + z_d_vn_hdf * r_dtime
              END IF
              !
              IF (p_nh_diag%ddt_vn_dyn_is_associated) THEN
                p_nh_diag%ddt_vn_dyn(je,jk,jb)  =  p_nh_diag%ddt_vn_dyn(je,jk,jb) + z_d_vn_hdf * r_dtime
              END IF
              !
            ENDDO
          ENDDO
        ENDIF
        !$ACC END PARALLEL


      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF

    ! Apply diffusion for case diffu_type = 4

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,rl_start,rl_end)

    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

    IF (diffu_type == 4) THEN  ! fourth-order diffusion

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_d_vn_hdf) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRIVATE(z_d_vn_hdf) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
        DO jk = 1, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            !
            z_d_vn_hdf = - p_patch%edges%area_edge(je,jb)*p_patch%edges%area_edge(je,jb) &
              &          * diff_multfac_vn(jk) * z_nabla4_e(je,jk,jb)
            !
            p_nh_prog%vn(je,jk,jb)            =  p_nh_prog%vn(je,jk,jb)         + z_d_vn_hdf
            !
            IF (p_nh_diag%ddt_vn_hdf_is_associated) THEN
              p_nh_diag%ddt_vn_hdf(je,jk,jb)  =  p_nh_diag%ddt_vn_hdf(je,jk,jb) + z_d_vn_hdf * r_dtime
            END IF
            !
            IF (p_nh_diag%ddt_vn_dyn_is_associated) THEN
              p_nh_diag%ddt_vn_dyn(je,jk,jb)  =  p_nh_diag%ddt_vn_dyn(je,jk,jb) + z_d_vn_hdf * r_dtime
            END IF
            !
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
!$OMP END DO

    ENDIF

    IF (l_limited_area .OR. jg > 1) THEN

      ! Lateral boundary diffusion for vn
      i_startblk = p_patch%edges%start_block(start_bdydiff_e)
      i_endblk   = p_patch%edges%end_block(grf_bdywidth_e)

!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx,z_d_vn_hdf) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, start_bdydiff_e, grf_bdywidth_e)

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRIVATE(z_d_vn_hdf) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
        DO jk = 1, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            !
            z_d_vn_hdf = p_patch%edges%area_edge(je,jb) * fac_bdydiff_v * z_nabla2_e(je,jk,jb)
            !
            p_nh_prog%vn(je,jk,jb)            =  p_nh_prog%vn(je,jk,jb)         + z_d_vn_hdf
            !
            IF (p_nh_diag%ddt_vn_hdf_is_associated) THEN
              p_nh_diag%ddt_vn_hdf(je,jk,jb)  =  p_nh_diag%ddt_vn_hdf(je,jk,jb) + z_d_vn_hdf * r_dtime
            END IF
            !
            IF (p_nh_diag%ddt_vn_dyn_is_associated) THEN
              p_nh_diag%ddt_vn_dyn(je,jk,jb)  =  p_nh_diag%ddt_vn_dyn(je,jk,jb) + z_d_vn_hdf * r_dtime
            END IF
            !
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
!$OMP END DO

    ENDIF ! vn boundary diffusion


    IF (diffusion_config(jg)%lhdiff_w) THEN ! add diffusion on vertical wind speed
                     ! remark: the surface level (nlevp1) is excluded because w is diagnostic there

      IF (l_limited_area .AND. jg == 1) THEN
        rl_start = grf_bdywidth_c+1
      ELSE
        rl_start = grf_bdywidth_c
      ENDIF
      rl_end   = min_rlcell_int-1

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            z_nabla2_c(jc,jk,jb) =  &
              p_nh_prog%w(jc,jk,jb)                        *p_int%geofac_n2s(jc,1,jb) + &
              p_nh_prog%w(icidx(jc,jb,1),jk,icblk(jc,jb,1))*p_int%geofac_n2s(jc,2,jb) + &
              p_nh_prog%w(icidx(jc,jb,2),jk,icblk(jc,jb,2))*p_int%geofac_n2s(jc,3,jb) + &
              p_nh_prog%w(icidx(jc,jb,3),jk,icblk(jc,jb,3))*p_int%geofac_n2s(jc,4,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

        IF (turbdiff_config(jg)%itype_sher == 2) THEN ! compute horizontal gradients of w
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
          DO jk = 2, nlev
            DO jc = i_startidx, i_endidx
             p_nh_diag%dwdx(jc,jk,jb) =  p_int%geofac_grg(jc,1,jb,1)*p_nh_prog%w(jc,jk,jb) + &
               p_int%geofac_grg(jc,2,jb,1)*p_nh_prog%w(icidx(jc,jb,1),jk,icblk(jc,jb,1))   + &
               p_int%geofac_grg(jc,3,jb,1)*p_nh_prog%w(icidx(jc,jb,2),jk,icblk(jc,jb,2))   + &
               p_int%geofac_grg(jc,4,jb,1)*p_nh_prog%w(icidx(jc,jb,3),jk,icblk(jc,jb,3))

             p_nh_diag%dwdy(jc,jk,jb) =  p_int%geofac_grg(jc,1,jb,2)*p_nh_prog%w(jc,jk,jb) + &
               p_int%geofac_grg(jc,2,jb,2)*p_nh_prog%w(icidx(jc,jb,1),jk,icblk(jc,jb,1))   + &
               p_int%geofac_grg(jc,3,jb,2)*p_nh_prog%w(icidx(jc,jb,2),jk,icblk(jc,jb,2))   + &
               p_int%geofac_grg(jc,4,jb,2)*p_nh_prog%w(icidx(jc,jb,3),jk,icblk(jc,jb,3))

            ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
        ENDIF

      ENDDO
!$OMP END DO

      IF (l_limited_area .AND. jg == 1) THEN
        rl_start = 0
      ELSE
        rl_start = grf_bdywidth_c+1
      ENDIF
      rl_end   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            p_nh_prog%w(jc,jk,jb) = p_nh_prog%w(jc,jk,jb) - diff_multfac_w * p_patch%cells%area(jc,jb)**2 * &
             (z_nabla2_c(jc,jk,jb)                        *p_int%geofac_n2s(jc,1,jb) +                      &
              z_nabla2_c(icidx(jc,jb,1),jk,icblk(jc,jb,1))*p_int%geofac_n2s(jc,2,jb) +                      &
              z_nabla2_c(icidx(jc,jb,2),jk,icblk(jc,jb,2))*p_int%geofac_n2s(jc,3,jb) +                      &
              z_nabla2_c(icidx(jc,jb,3),jk,icblk(jc,jb,3))*p_int%geofac_n2s(jc,4,jb))
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

        ! add Smagorinsky diffusion on vertical wind speed
        IF (diffusion_config(jg)%lhdiff_smag_w) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
          DO jk = 2, nlev
            DO jc = i_startidx, i_endidx
              p_nh_prog%w(jc,jk,jb) = p_nh_prog%w(jc,jk,jb) + p_patch%cells%area(jc,jb)*z_nabla2_c(jc,jk,jb) * &
                MAX(diff_multfac_n2w(jk) ,                                                                     &
                 0.5_wp*(kh_smag_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)   +              &
                         kh_smag_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)   +              &
                         kh_smag_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb)   +              &
                         kh_smag_e(ieidx(jc,jb,1),jk-1,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb) +              &
                         kh_smag_e(ieidx(jc,jb,2),jk-1,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb) +              &
                         kh_smag_e(ieidx(jc,jb,3),jk-1,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb)) )
            ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
        ELSE
          ! Add nabla2 diffusion in upper damping layer (if present)
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
          DO jk = 2, nrdmax(jg)
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              p_nh_prog%w(jc,jk,jb) = p_nh_prog%w(jc,jk,jb) +                         &
                diff_multfac_n2w(jk) * p_patch%cells%area(jc,jb) * z_nabla2_c(jc,jk,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
        ENDIF

      ENDDO
!$OMP END DO

    ENDIF ! w diffusion

!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_E, p_patch, p_nh_prog%vn,opt_varname="diffusion: vn sync")

    IF (ltemp_diffu) THEN ! Smagorinsky temperature diffusion

      IF (l_zdiffu_t) THEN
        icell      => p_nh_metrics%zd_indlist
        iblk       => p_nh_metrics%zd_blklist
        ilev       => p_nh_metrics%zd_vertidx
   !     iedge      => p_nh_metrics%zd_edgeidx
   !     iedblk     => p_nh_metrics%zd_edgeblk
        vcoef      => p_nh_metrics%zd_intcoef
   !     blcoef     => p_nh_metrics%zd_e2cell
        geofac_n2s => p_nh_metrics%zd_geofac

        nproma_zdiffu = cpu_min_nproma(nproma,256)
        nblks_zdiffu = INT(p_nh_metrics%zd_listdim/nproma_zdiffu)
        npromz_zdiffu = MOD(p_nh_metrics%zd_listdim,nproma_zdiffu)
        IF (npromz_zdiffu > 0) THEN
          nblks_zdiffu = nblks_zdiffu + 1
        ELSE
          npromz_zdiffu = nproma_zdiffu
        ENDIF
      ENDIF

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

      ! Enhance Smagorinsky diffusion coefficient in the presence of excessive grid-point cold pools
      ! This is restricted to the two lowest model levels
      !
      rl_start = grf_bdywidth_c
      rl_end   = min_rlcell_int-1

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx,ic,tdiff,trefdiff), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ic = 0

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
        DO jk = nlev-1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            ! Perturbation potential temperature difference between local point and average of the three neighbors
            tdiff = p_nh_prog%theta_v(jc,jk,jb) -                          &
              (p_nh_prog%theta_v(icidx(jc,jb,1),jk,icblk(jc,jb,1)) +       &
               p_nh_prog%theta_v(icidx(jc,jb,2),jk,icblk(jc,jb,2)) +       &
               p_nh_prog%theta_v(icidx(jc,jb,3),jk,icblk(jc,jb,3)) ) / 3._wp
            trefdiff = p_nh_metrics%theta_ref_mc(jc,jk,jb) -                       &
              (p_nh_metrics%theta_ref_mc(icidx(jc,jb,1),jk,icblk(jc,jb,1)) +       &
               p_nh_metrics%theta_ref_mc(icidx(jc,jb,2),jk,icblk(jc,jb,2)) +       &
               p_nh_metrics%theta_ref_mc(icidx(jc,jb,3),jk,icblk(jc,jb,3)) ) / 3._wp

            ! Enahnced horizontal diffusion is applied if the theta perturbation is either
            ! - at least 5 K colder than the average of the neighbor points on valley points (determined by trefdiff < 0.) or
            ! - at least 7.5 K colder than the average of the neighbor points otherwise
            IF (tdiff-trefdiff < thresh_tdiff .AND. trefdiff < 0._wp .OR. tdiff-trefdiff < 1.5_wp*thresh_tdiff) THEN
              ic = ic+1
              iclist(ic,jb) = jc
              iklist(ic,jb) = jk
              tdlist(ic,jb) = thresh_tdiff - tdiff + trefdiff
            ENDIF
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
        icount(jb) = ic

      ENDDO
!$OMP END DO

      ! Enhance Smagorinsky coefficients at the three edges of the cells included in the list
      ! Attention: this operation is neither vectorizable nor OpenMP-parallelizable (race conditions!)

!$OMP MASTER
      DO jb = i_startblk,i_endblk

        IF (icount(jb) > 0) THEN
          DO ic = 1, icount(jb)
            jc = iclist(ic,jb)
            jk = iklist(ic,jb)
            enh_diffu = tdlist(ic,jb)*5.e-4_vp
            kh_smag_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) = MAX(enh_diffu,kh_smag_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)))
            kh_smag_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) = MAX(enh_diffu,kh_smag_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)))
            kh_smag_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) = MAX(enh_diffu,kh_smag_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)))
          ENDDO
        ENDIF

     ENDDO

!$OMP END MASTER
!$OMP BARRIER


      IF (discr_t == 1) THEN  ! use discretization K*nabla(theta)

        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int

        i_startblk = p_patch%cells%start_block(rl_start)
        i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! interpolated diffusion coefficient times nabla2(theta)
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              z_temp(jc,jk,jb) =  &
               (kh_smag_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)          + &
                kh_smag_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)          + &
                kh_smag_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb))         * &
               (p_nh_prog%theta_v(jc,jk,jb)                        *p_int%geofac_n2s(jc,1,jb) + &
                p_nh_prog%theta_v(icidx(jc,jb,1),jk,icblk(jc,jb,1))*p_int%geofac_n2s(jc,2,jb) + &
                p_nh_prog%theta_v(icidx(jc,jb,2),jk,icblk(jc,jb,2))*p_int%geofac_n2s(jc,3,jb) + &
                p_nh_prog%theta_v(icidx(jc,jb,3),jk,icblk(jc,jb,3))*p_int%geofac_n2s(jc,4,jb))
            ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
        ENDDO
!$OMP END DO

      ELSE IF (discr_t == 2) THEN ! use conservative discretization div(k*grad(theta))

        rl_start = grf_bdywidth_e
        rl_end   = min_rledge_int - 1

        i_startblk = p_patch%edges%start_block(rl_start)
        i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! compute kh_smag_e * grad(theta) (stored in z_nabla2_e for memory efficiency)
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              z_nabla2_e(je,jk,jb) = kh_smag_e(je,jk,jb) *              &
                p_patch%edges%inv_dual_edge_length(je,jb)*              &
               (p_nh_prog%theta_v(iecidx(je,jb,2),jk,iecblk(je,jb,2)) - &
                p_nh_prog%theta_v(iecidx(je,jb,1),jk,iecblk(je,jb,1)))
            ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
        ENDDO
!$OMP END DO

        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int

        i_startblk = p_patch%cells%start_block(rl_start)
        i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! now compute the divergence of the quantity above
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              z_temp(jc,jk,jb) =                                                         &
                z_nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
                z_nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
                z_nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL LOOP
        ENDDO
!$OMP END DO

      ENDIF


      IF (l_zdiffu_t) THEN ! Compute temperature diffusion truly horizontally over steep slopes
                           ! A conservative discretization is not possible here
!$OMP DO PRIVATE(jb,jc,ic,nlen_zdiffu,ishift) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, nblks_zdiffu
          IF (jb == nblks_zdiffu) THEN
            nlen_zdiffu = npromz_zdiffu
          ELSE
            nlen_zdiffu = nproma_zdiffu
          ENDIF
          ishift = (jb-1)*nproma_zdiffu
!$NEC ivdep
!DIR$ IVDEP
          DO jc = 1, nlen_zdiffu
            ic = ishift+jc
            z_temp(icell(1,ic),ilev(1,ic),iblk(1,ic)) =                                          &
              z_temp(icell(1,ic),ilev(1,ic),iblk(1,ic)) + p_nh_metrics%zd_diffcoef(ic)*          &
!              MAX(p_nh_metrics%zd_diffcoef(ic),        &
!              kh_smag_e(iedge(1,ic),ilev(1,ic),iedblk(1,ic))* blcoef(1,ic)  +                    &
!              kh_smag_e(iedge(2,ic),ilev(1,ic),iedblk(2,ic))* blcoef(2,ic)  +                    &
!              kh_smag_e(iedge(3,ic),ilev(1,ic),iedblk(3,ic))* blcoef(3,ic) ) *                   &
             (geofac_n2s(1,ic)*p_nh_prog%theta_v(icell(1,ic),ilev(1,ic),iblk(1,ic)) +            &
              geofac_n2s(2,ic)*(vcoef(1,ic)*p_nh_prog%theta_v(icell(2,ic),ilev(2,ic),iblk(2,ic))+&
              (1._wp-vcoef(1,ic))* p_nh_prog%theta_v(icell(2,ic),ilev(2,ic)+1,iblk(2,ic)))  +    &
              geofac_n2s(3,ic)*(vcoef(2,ic)*p_nh_prog%theta_v(icell(3,ic),ilev(3,ic),iblk(3,ic))+&
              (1._wp-vcoef(2,ic))*p_nh_prog%theta_v(icell(3,ic),ilev(3,ic)+1,iblk(3,ic)))  +     &
              geofac_n2s(4,ic)*(vcoef(3,ic)*p_nh_prog%theta_v(icell(4,ic),ilev(4,ic),iblk(4,ic))+&
              (1._wp-vcoef(3,ic))* p_nh_prog%theta_v(icell(4,ic),ilev(4,ic)+1,iblk(4,ic)))  )
          ENDDO
        ENDDO
        !$ACC END PARALLEL
!$OMP END DO

      ENDIF

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx,z_theta) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            z_theta = p_nh_prog%theta_v(jc,jk,jb)

            p_nh_prog%theta_v(jc,jk,jb) = p_nh_prog%theta_v(jc,jk,jb) + &
              p_patch%cells%area(jc,jb)*z_temp(jc,jk,jb)

            p_nh_prog%exner(jc,jk,jb) = p_nh_prog%exner(jc,jk,jb) *      &
              (1._wp+rd_o_cvd*(p_nh_prog%theta_v(jc,jk,jb)/z_theta-1._wp))

          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! This could be further optimized, but applications without physics are quite rare;
      IF ( linit .OR. (iforcing /= inwp .AND. iforcing /= iaes) ) THEN
        CALL sync_patch_array_mult(SYNC_C,p_patch,2,p_nh_prog%theta_v,p_nh_prog%exner,  &
                                   opt_varname="diffusion: theta and exner")
      ENDIF

    ENDIF ! temperature diffusion


    IF ( linit .OR. (iforcing /= inwp .AND. iforcing /= iaes) ) THEN
      IF (diffusion_config(jg)%lhdiff_w) THEN
        CALL sync_patch_array(SYNC_C,p_patch,p_nh_prog%w,"diffusion: w")
      END IF
    ENDIF

    IF (ltimer) CALL timer_stop(timer_nh_hdiffusion)

    !$ACC WAIT(1)
    !$ACC END DATA


  END SUBROUTINE diffusion


END MODULE mo_nh_diffusion
