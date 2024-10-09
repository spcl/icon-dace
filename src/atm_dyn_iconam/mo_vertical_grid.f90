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

! mo_vertical_grid provides all subroutines for handling of the vertical
! grid in the nonhydrostatic model branch.
! Routines are provided for
!  1) Creating the hybrid height based coordinate system.
!  2) (De)Allocations of related metric fields. Note, that they are defined
!     in mo_model_domain

!----------------------------
#include "omp_definitions.inc"
!----------------------------


MODULE mo_vertical_grid

  USE mo_kind,                  ONLY: wp, vp
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_model_domain,          ONLY: t_patch
  USE mo_ext_data_types,        ONLY: t_external_data
  USE mo_grid_config,           ONLY: n_dom, grid_sphere_radius
  USE mo_nonhydrostatic_config, ONLY: rayleigh_type, rayleigh_coeff, damp_height, &
    &                                 igradp_method, vwind_offctr,                &
    &                                 exner_expol, l_zdiffu_t, thslp_zdiffu,      &
    &                                 thhgtd_zdiffu, kstart_dd3d, divdamp_type,   &
    &                                 divdamp_trans_start, divdamp_trans_end
  USE mo_diffusion_config,      ONLY: diffusion_config
  USE mo_parallel_config,       ONLY: nproma, p_test_run
  USE mo_run_config,            ONLY: msg_level, iforcing
  USE mo_vertical_coord_table,  ONLY: vct_a
  USE mo_impl_constants,        ONLY: max_dom, RAYLEIGH_CLASSIC, &
    &                                 RAYLEIGH_KLEMP, min_rlcell_int, min_rlcell, min_rledge_int, &
    &                                 SUCCESS, iaes, ildf_echam
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e, grf_fbk_start_c, &
                                      grf_nudge_start_e, grf_nudgezone_width
  USE mo_physical_constants,    ONLY: grav, p0ref, rd, rd_o_cpd, cpd, p0sl_bg, rgrav
  USE mo_math_gradients,        ONLY: grad_fd_norm, grad_fd_tang
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp,                  ONLY: cells2edges_scalar, cells2verts_scalar, cell_avg
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_cell
  USE mo_math_constants,        ONLY: pi_2
  USE mo_loopindices,           ONLY: get_indices_e, get_indices_c
  USE mo_nonhydro_types,        ONLY: t_nh_state, t_nh_state_lists, t_nh_metrics
  USE mo_init_vgrid,            ONLY: nflatlev
  USE mo_util_vgrid_types,      ONLY: vgrid_buffer
  USE mo_sync,                  ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, global_sum_array, &
                                      sync_patch_array_mult, global_min, global_max
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_les_config,           ONLY: les_config
  USE mo_impl_constants,       ONLY: min_rlvert_int
  USE turb_data,               ONLY: akt
  USE mo_fortran_tools,        ONLY: init
  USE mo_util_string,          ONLY: int2string, real2string
  USE mo_mpi,                  ONLY: my_process_is_stdio
  USE mo_util_table,           ONLY: t_table, initialize_table, add_table_column, &
    &                                set_table_entry, print_table, finalize_table
  USE mo_nudging_config,       ONLY: nudging_config, indg_type, indg_profile
  USE mo_dynamics_config,      ONLY: ldeepatmo
  USE mo_aes_vdf_config,       ONLY: aes_vdf_config
  USE mo_turb_vdiff_params,    ONLY: VDIFF_TURB_3DSMAGORINSKY
  USE mo_var_list,             ONLY: t_var_list_ptr
  USE mo_nonhydro_state,       ONLY: new_zd_metrics  
  USE mo_deepatmo,             ONLY: deepatmo_htrafo
  IMPLICIT NONE

  PRIVATE


  ! Constants used for the computation of the background reference atmosphere of the nh-model
  !
  REAL(wp), PARAMETER :: h_scal_bg = 10000._wp    ! [m]      scale height
  REAL(wp), PARAMETER :: t0sl_bg   = 288.15_wp    ! [K]      sea level temperature
  REAL(wp), PARAMETER :: del_t_bg  = 75._wp       ! [K]      difference between sea level
  !                                                          temperature and asymptotic
  !                                                          stratospheric temperature

  INTEGER:: nrdmax(max_dom), nflat_gradp(max_dom)

  PUBLIC :: set_nh_metrics, nrdmax, nflat_gradp

  CHARACTER(*), PARAMETER :: modname = 'mo_vertical_grid'
  CONTAINS


  !----------------------------------------------------------------------------
  !>
  !! Initializes values for the fields in type nh_metrics of the NH state.
  !!
  SUBROUTINE set_nh_metrics(p_patch, p_nh, p_nh_lists, p_int, ext_data)

    CHARACTER(len=*), PARAMETER :: routine = modname//':set_nh_metrics'

    TYPE(t_patch),     TARGET, INTENT(INOUT) :: p_patch(n_dom)  !< patch
    TYPE(t_nh_state),          INTENT(INOUT) :: p_nh(n_dom)
    TYPE(t_int_state), TARGET, INTENT(   IN) :: p_int(n_dom)
    TYPE(t_nh_state_lists),    INTENT(INOUT) :: p_nh_lists(n_dom)
    TYPE(t_external_data),     INTENT(INOUT) :: ext_data(n_dom)

    INTEGER :: jg, jk, jk1, jk_start, jb, jc, je, jn, jgc, nlen, &
               nblks_c, npromz_c, nblks_e, npromz_e, nblks_v, ic
    INTEGER :: nlev, nlevp1              !< number of full levels

    ! Note: the present way of setting up coordinate surfaces will not work for vertical refinement
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk, icount_total
    INTEGER :: ica(max_dom)

    REAL(wp) :: z_diff, z_sin_diff, z_sin_diff_full, z_tanh_diff,      &
      &         z1, z2, z3, zf, z_help(nproma),                        &
      &         z_temp(nproma), z_aux1(nproma), z_aux2(nproma),        &
      &         z0, coef1, coef2, coef3, dn1, dn2, dn3, dn4, dn5, dn6, &
      &         zn_min, zn_avg, zn_sq, zn_rms
    REAL(wp) :: z_maxslope, z_maxhdiff, z_offctr
    REAL(wp), ALLOCATABLE :: z_ifv(:,:,:)
    REAL(wp), ALLOCATABLE :: z_me(:,:,:),z_maxslp(:,:,:),z_maxhgtd(:,:,:),z_shift(:,:,:), &
                             z_ddxt_z_half_e(:,:,:), z_ddxn_z_half_e(:,:,:)
    REAL(wp), ALLOCATABLE :: z_aux_c(:,:,:), z_aux_c2(:,:,:), z_aux_e(:,:,:)
    REAL(wp) :: extrapol_dist
    INTEGER,  ALLOCATABLE :: flat_idx(:,:), imask(:,:,:),icount(:)
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk, inidx, inblk
    LOGICAL :: l_found(nproma), lfound_all, is_les_phy
    INTEGER :: error_status

#ifdef INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: ica,z_help,z_temp,z_aux1,z_aux2,l_found
!DIR$ ATTRIBUTES ALIGN : 64 :: z_ifv,z_me,z_maxslp,z_maxhgtd,z_shift,z_ddxt_z_half_e
!DIR$ ATTRIBUTES ALIGN : 64 :: z_ddxn_z_half_e,z_aux_c,z_aux_c2,z_aux_e
!DIR$ ATTRIBUTES ALIGN : 64 :: flat_idx,imask,icount
#endif
    !------------------------------------------------------------------------

    DO jg = 1,n_dom

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c
      nblks_e   = p_patch(jg)%nblks_e
      npromz_e  = p_patch(jg)%npromz_e
      nblks_v   = p_patch(jg)%nblks_v

      ! number of vertical levels
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1


!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1,nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
          p_nh(jg)%metrics%z_mc(nlen+1:nproma,:,jb) = 0._wp
          p_nh(jg)%metrics%z_ifc(nlen+1:nproma,:,jb) = 0._wp
        ENDIF
        DO jk = 1, nlev+1
          ! geometric height of half levels
          !   The 3D coordinate field "z_ifc" exists already in a buffer
          !   variable of module "mo_util_vgrid_types". We move the data to its
          !   final place here:
          p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb) = vgrid_buffer(jg)%z_ifc(1:nlen,jk,jb)
        ENDDO
        DO jk = 1, nlev
          ! geometric height of full levels
          p_nh(jg)%metrics%z_mc(1:nlen,jk,jb) = 0.5_wp*(p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb) + &
            &                                           p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb))
        ENDDO
      ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1,nblks_c

        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
           ! To avoid uninitialized field elements going in the output routine
           p_nh(jg)%metrics%geopot(nlen+1:nproma,:,jb) = 0._wp
        ENDIF

        DO jk = 1, nlev
          ! geopotential on full levels
          p_nh(jg)%metrics%geopot(1:nlen,jk,jb) = grav*p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)
        ENDDO

        DO jk = 1, nlev
          ! functional determinant of the metrics (is positive), full levels
          p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk,jb) = &
          & p_nh(jg)%metrics%z_ifc(1:nlen,jk  ,jb)- &
          & p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb)
          ! difference of geopotential between the levels
          p_nh(jg)%metrics%dgeopot_mc (1:nlen,jk,jb) = grav * &
          &         ( p_nh(jg)%metrics%z_ifc(1:nlen,jk  ,jb) &
          &         - p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb) )
          ! inverse layer thickness (for runtime optimization)
          p_nh(jg)%metrics%inv_ddqz_z_full(1:nlen,jk,jb) = &
            1._wp/p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk,jb)
        ENDDO
        ! functional determinant of the metrics (is positive), half levels
        DO jk = 2, nlev
          p_nh(jg)%metrics%ddqz_z_half(1:nlen,jk,jb) = &
          & (p_nh(jg)%metrics%z_mc(1:nlen,jk-1,jb)     &
          & -p_nh(jg)%metrics%z_mc(1:nlen,jk  ,jb))
        ENDDO
        p_nh(jg)%metrics%ddqz_z_half(1:nlen,1,jb) =   &
        & 2.0_wp*(p_nh(jg)%metrics%z_ifc(1:nlen,1,jb) &
        &       - p_nh(jg)%metrics%z_mc (1:nlen,1,jb))
        p_nh(jg)%metrics%ddqz_z_half(1:nlen,nlevp1,jb) =    &
        & 2.0_wp*(p_nh(jg)%metrics%z_mc (1:nlen,nlev  ,jb)  &
        &       - p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb))
        ! coefficients for second-order acurate dw/dz term
        p_nh(jg)%metrics%coeff1_dwdz(:,1,jb) = 0._wp
        p_nh(jg)%metrics%coeff2_dwdz(:,1,jb) = 0._wp
        DO jk = 2, nlev
          p_nh(jg)%metrics%coeff1_dwdz(1:nlen,jk,jb) = &
            p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk,jb)/p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk-1,jb) / &
            ( p_nh(jg)%metrics%z_ifc(1:nlen,jk-1,jb) - p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb) )
          p_nh(jg)%metrics%coeff2_dwdz(1:nlen,jk,jb) = &
            p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk-1,jb)/p_nh(jg)%metrics%ddqz_z_full(1:nlen,jk,jb) / &
            ( p_nh(jg)%metrics%z_ifc(1:nlen,jk-1,jb) - p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb) )
        ENDDO


        ! surface geopotential
        !
        ext_data(jg)%atm%fis(1:nlen,jb) = grav*p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb)


        !-------------------------------------------------------------
        ! geopot above ground  - because physics needs positive values
        !-------------------------------------------------------------

          p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,nlevp1,jb) = 0._wp

          DO jk = nlev,1,-1
            ! geopotential (interfaces)
            p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,jk,jb)=             &
          &            p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,jk+1,jb) &
          &         +  grav *                                      &
          &          ( p_nh(jg)%metrics%z_ifc(1:nlen,jk  ,jb)      &
          &          - p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb))
          ENDDO

          ! geopotential above ground (at full levels)
          DO jk = 2,nlevp1
            p_nh(jg)%metrics%geopot_agl(1:nlen,jk-1,jb) =  0.5_wp *           &
           &              ( p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,jk-1  ,jb )&
           &              + p_nh(jg)%metrics%geopot_agl_ifc(1:nlen,jk    ,jb))
          ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! For the tangential slope we need temporarily also the height
      ! at the vertices
      ALLOCATE(z_ifv(nproma,nlevp1,nblks_v))
      ALLOCATE(z_ddxt_z_half_e(nproma,nlevp1,nblks_e))
      ALLOCATE(z_ddxn_z_half_e(nproma,nlevp1,nblks_e))
      ! Intermediate storage for fields that are optionally single precision (to circumvent duplicating subroutines)
      ALLOCATE(z_aux_e(nproma,nlev,nblks_e))

      z_ifv = 0._wp
      CALL cells2verts_scalar(p_nh(jg)%metrics%z_ifc, p_patch(jg), p_int(jg)%cells_aw_verts, z_ifv)

      ! Start index for slope computations
      i_startblk = p_patch(jg)%edges%start_block(2)

      z_ddxn_z_half_e = 0._wp
      z_ddxt_z_half_e = 0._wp
      p_nh(jg)%metrics%ddxt_z_full = 0._wp
      p_nh(jg)%metrics%ddxn_z_full = 0._wp

      ! slope of the terrain (tangential direction)
      CALL grad_fd_tang ( z_ifv,                  &
           &              p_patch(jg),            &
           &              z_ddxt_z_half_e, &
           &              1, nlevp1 )
      DEALLOCATE(z_ifv)

      ! slope of the terrain (normal direction)
      CALL grad_fd_norm ( p_nh(jg)%metrics%z_ifc, &
           &              p_patch(jg), &
           &              z_ddxn_z_half_e, &
           &              1, nlevp1 )

      CALL sync_patch_array_mult(SYNC_E,p_patch(jg),2,z_ddxt_z_half_e, &
        z_ddxn_z_half_e)
#ifndef __NO_ICON_LES__
      IF (atm_phy_nwp_config(jg)%is_les_phy) THEN
        ! remark: ddxt_z_half_e, ddxn_z_half_e in p_nh(jg)%metrics are optionally single precision
        p_nh(jg)%metrics%ddxt_z_half_e(:,:,:) = z_ddxt_z_half_e(:,:,:)
        p_nh(jg)%metrics%ddxn_z_half_e(:,:,:) = z_ddxn_z_half_e(:,:,:)
      ENDIF
#endif
      IF (aes_vdf_config(jg)%turb == VDIFF_TURB_3DSMAGORINSKY) THEN
        ! remark: ddxt_z_half_e, ddxn_z_half_e in p_nh(jg)%metrics are optionally single precision
        p_nh(jg)%metrics%ddxt_z_half_e(:,:,:) = z_ddxt_z_half_e(:,:,:)
        p_nh(jg)%metrics%ddxn_z_half_e(:,:,:) = z_ddxn_z_half_e(:,:,:)
      ENDIF

      ! vertically averaged metrics
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            p_nh(jg)%metrics%ddxn_z_full(je,jk,jb) = 0.5_wp * &
            & (z_ddxn_z_half_e(je,jk,jb) + z_ddxn_z_half_e(je,jk+1,jb))

            p_nh(jg)%metrics%ddxt_z_full(je,jk,jb) = 0.5_wp * &
            & (z_ddxt_z_half_e(je,jk,jb) + z_ddxt_z_half_e(je,jk+1,jb))
          ENDDO
        ENDDO


        ! Coefficients for improved calculation of kinetic energy gradient
        DO je = i_startidx, i_endidx
          p_nh(jg)%metrics%coeff_gradekin(je,1,jb) = p_patch(jg)%edges%edge_cell_length(je,jb,2) / &
            p_patch(jg)%edges%edge_cell_length(je,jb,1) * p_patch(jg)%edges%inv_dual_edge_length(je,jb)
          p_nh(jg)%metrics%coeff_gradekin(je,2,jb) = p_patch(jg)%edges%edge_cell_length(je,jb,1) / &
            p_patch(jg)%edges%edge_cell_length(je,jb,2) * p_patch(jg)%edges%inv_dual_edge_length(je,jb)
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      z_aux_e = 0._wp

      ! functional determinant at full level edges
      CALL cells2edges_scalar(p_nh(jg)%metrics%ddqz_z_full, &
           &                  p_patch(jg), p_int(jg)%c_lin_e, &
           &                  z_aux_e )

      ! remark: ddqz_z_full_e is optionally single precision
      p_nh(jg)%metrics%ddqz_z_full_e(:,:,:) = z_aux_e(:,:,:)
      CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%ddqz_z_full_e)

      DEALLOCATE(z_aux_e)

      ! slope angle and slope azimuth (used for slope-dependent radiation)
      ALLOCATE (z_aux_c(nproma,1,nblks_c), z_aux_c2(nproma,1,nblks_c))
      z_aux_c(:,:,:) = 0._wp ; z_aux_c2(:,:,:) = 0._wp
      
      CALL rbf_vec_interpol_cell(z_ddxn_z_half_e(:,nlevp1:nlevp1,:), p_patch(jg), p_int(jg), z_aux_c, z_aux_c2)

      i_startblk = p_patch(jg)%cells%start_block(2)

!$OMP PARALLEL DO PRIVATE(jb, i_startidx, i_endidx, jc)
      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch(jg), jb, i_startblk, nblks_c, i_startidx, i_endidx, 2)

! TEMPORARY WORKAROUND: Without novector, NEC 3.3.1 and 3.5.1 do not compile the following code correctly
!$NEC novector
        DO jc = i_startidx, i_endidx
          p_nh(jg)%metrics%slope_angle(jc,jb)   = ATAN(SQRT(z_aux_c(jc,1,jb)**2+z_aux_c2(jc,1,jb)**2))
          IF (z_aux_c(jc,1,jb) /= 0._wp .OR. z_aux_c2(jc,1,jb) /= 0._wp) THEN
            p_nh(jg)%metrics%slope_azimuth(jc,jb) = ATAN2(z_aux_c(jc,1,jb),z_aux_c2(jc,1,jb))
          ELSE
            p_nh(jg)%metrics%slope_azimuth(jc,jb) = 0._wp
          ENDIF
        ENDDO

      ENDDO
!$OMP END PARALLEL DO


      ! offcentering in vertical mass flux
      p_nh(jg)%metrics%vwind_impl_wgt(:,:)    = 0.5_wp + vwind_offctr
      p_nh(jg)%metrics%vwind_expl_wgt(:,:)    = 1.0_wp - p_nh(jg)%metrics%vwind_impl_wgt(:,:)

      ! Rayleigh damping properties

      p_nh(jg)%metrics%rayleigh_w(:)   = 0.0_wp
      p_nh(jg)%metrics%enhfac_diffu(:) = 1.0_wp

      ! Determine end index of damping layer
      nrdmax(jg) = 1
      DO jk = 2, nlevp1
        jk1 = jk + p_patch(jg)%nshift_total
        IF (vct_a(jk1) >= damp_height(jg))  nrdmax(jg) = jk
      ENDDO


      ! Rayleigh damping coefficient for vn and/or w
      DO jk = 1, nrdmax(jg)
        jk1 = jk + p_patch(jg)%nshift_total

        ! z - z_h  with z at half levels
        z_sin_diff  = MAX(0.0_wp,vct_a(jk1)-damp_height(jg))

        ! z_top - z
        z_tanh_diff = vct_a(1) - vct_a(jk1)

        ! z - z_h  with z at full levels
        z_sin_diff_full = MAX(0.0_wp,0.5_wp*(vct_a(jk1)+vct_a(jk1+1))-damp_height(jg))

        IF ( rayleigh_type == RAYLEIGH_CLASSIC ) THEN
          ! classical Rayleigh damping as known from idealized limited area test cases.
          ! Requires knowledge of a velocity reference state!
          !
          p_nh(jg)%metrics%rayleigh_vn(jk)= rayleigh_coeff(jg)*(SIN(pi_2*z_sin_diff_full/ &
              MAX(1.e-3_wp,vct_a(p_patch(jg)%nshift_total+1)-damp_height(jg))))**2

          p_nh(jg)%metrics%rayleigh_w(jk)= rayleigh_coeff(jg)*(SIN(pi_2*z_sin_diff/ &
              MAX(1.e-3_wp,vct_a(p_patch(jg)%nshift_total+1)-damp_height(jg))))**2

        ELSE IF ( rayleigh_type == RAYLEIGH_KLEMP ) THEN
          ! Rayleigh damping based on Klemp et al. (2008), MRW 136, pp 3987-4004
          ! No reference state needed, thus applicable to real cases!
          !
          p_nh(jg)%metrics%rayleigh_w(jk)= rayleigh_coeff(jg)*&
          (1._wp-TANH(3.8_wp*z_tanh_diff/MAX(1.e-6_wp,vct_a(1)-damp_height(jg))))

          p_nh(jg)%metrics%rayleigh_vn(jk)= 0._wp
        ENDIF
      ENDDO

      ! Enhancement coefficient for nabla4 background diffusion near model top
      DO jk = 1, nrdmax(jg)
        jk1 = jk + p_patch(jg)%nshift_total
!        z_diff = MAX(0.0_wp,0.5_wp*(vct_a(jk1)+vct_a(jk1+1))-damp_height(jg))
        z_diff = 0.5_wp*(vct_a(1)+vct_a(2))-0.5_wp*(vct_a(jk1)+vct_a(jk1+1))
        p_nh(jg)%metrics%enhfac_diffu(jk) = 1._wp + &
 !         (hdiff_efdt_ratio/hdiff_min_efdt_ratio-1._wp)*(SIN(pi_2*z_diff/ &
 !         MAX(1.e-3_wp,0.5_wp*(vct_a(1)+vct_a(2))-damp_height(jg))))**2
          (diffusion_config(jg)%hdiff_efdt_ratio/  &
           diffusion_config(jg)%hdiff_min_efdt_ratio-1._wp)* &
          (1._wp-TANH(3.8_wp*z_diff                          &
          /MAX(1.e-6_wp,0.5_wp*(vct_a(1)+vct_a(2))-damp_height(jg))))
      ENDDO

      IF (msg_level >= 10) THEN
        WRITE(message_text,'(a,i4,a,i4)') 'Domain', jg, &
          '; end index of Rayleigh damping layer for w: ', nrdmax(jg)
        CALL message(routine, message_text)
        WRITE(message_text,'(a)') &
          'Damping coefficient for w; diffusion enhancement coefficient:'
        CALL message(modname, message_text)
        DO jk = 1, nrdmax(jg)
          jk1 = jk + p_patch(jg)%nshift_total
          WRITE(message_text,'(a,i5,a,f8.1,2e13.5)') 'level',jk, &
               ', half-level height',vct_a(jk1), &
            p_nh(jg)%metrics%rayleigh_w(jk),p_nh(jg)%metrics%enhfac_diffu(jk)
          CALL message(modname, message_text)
        ENDDO
      ENDIF

      ! Scaling factor for 3D divergence damping terms, and start level from which they are > 0
      kstart_dd3d(jg) = 1
      p_nh(jg)%metrics%scalfac_dd3d(:) = 1._wp
      IF (divdamp_type == 32) THEN  ! 3D divergence damping transitioning into 2D divergence damping
        kstart_dd3d(jg) = 0
        DO jk = 1, nlev
          jk1 = jk + p_patch(jg)%nshift_total
          zf = 0.5_wp*(vct_a(jk1)+vct_a(jk1+1))
          IF (zf >= divdamp_trans_end) THEN
            p_nh(jg)%metrics%scalfac_dd3d(jk) = 0._wp
            kstart_dd3d(jg) = jk
          ELSE IF (zf >= divdamp_trans_start) THEN
            p_nh(jg)%metrics%scalfac_dd3d(jk) = (divdamp_trans_end-zf)/(divdamp_trans_end-divdamp_trans_start)
          ELSE
            p_nh(jg)%metrics%scalfac_dd3d(jk) = 1._wp
          ENDIF
        ENDDO
        kstart_dd3d(jg) = kstart_dd3d(jg) + 1
      ENDIF

      ! Horizontal mask field for 3D divergence damping term; 2D div damping is generally applied in the 
      ! immediate vicinity of nest boundaries
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 1)

        DO je = i_startidx, i_endidx
          IF (p_patch(jg)%edges%refin_ctrl(je,jb) <= 0 .OR.                                           &
              p_patch(jg)%edges%refin_ctrl(je,jb) >= grf_nudge_start_e+2*(grf_nudgezone_width-1) ) THEN
            p_nh(jg)%metrics%hmask_dd3d(je,jb) = 1._wp
          ELSE IF (p_patch(jg)%edges%refin_ctrl(je,jb) <= grf_nudge_start_e+grf_nudgezone_width-1) THEN
            p_nh(jg)%metrics%hmask_dd3d(je,jb) = 0._wp
          ELSE
            p_nh(jg)%metrics%hmask_dd3d(je,jb) = 1._wp/REAL(grf_nudgezone_width-1,wp) *              &
              REAL(p_patch(jg)%edges%refin_ctrl(je,jb) - (grf_nudge_start_e+grf_nudgezone_width-1), wp)
          ENDIF
        ENDDO
      ENDDO

      ! Compute variable Exner extrapolation factors and offcentering coefficients for the
      ! vertical wind solver to optimize numerical stability over steep orography
      iidx  => p_patch(jg)%cells%edge_idx
      iblk  => p_patch(jg)%cells%edge_blk
      inidx => p_int(jg)%rbf_c2grad_idx
      inblk => p_int(jg)%rbf_c2grad_blk

      i_startblk = p_patch(jg)%cells%start_block(2)


      ALLOCATE (z_maxslp(nproma,nlev,nblks_c), z_maxhgtd(nproma,nlev,nblks_c) )

!$OMP PARALLEL
      ! Initialization to ensure that values are properly set at lateral boundaries
      CALL init(p_nh(jg)%metrics%exner_exfac(:,:,:), exner_expol, lacc=.FALSE.)
      CALL init(z_maxslp(:,:,:), lacc=.FALSE.)
      CALL init(z_maxhgtd(:,:,:), lacc=.FALSE.)
      CALL init(z_aux_c(:,:,:), lacc=.FALSE.)
!$OMP BARRIER

!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jk1, jc, z_maxslope, z_offctr, z_diff, &
!$OMP            zn_min, zn_avg, zn_sq, zn_rms) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch(jg), jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

            z_maxslp(jc,jk,jb) =                                                     &
              MAX(ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,1),jk,iblk(jc,jb,1))), &
                  ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,2),jk,iblk(jc,jb,2))), &
                  ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,3),jk,iblk(jc,jb,3))) )

            z_maxhgtd(jc,jk,jb) =                                                  &
              MAX(ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,1),jk,iblk(jc,jb,1)) &
               * p_patch(jg)%edges%dual_edge_length(iidx(jc,jb,1),iblk(jc,jb,1))), &
                  ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,2),jk,iblk(jc,jb,2)) &
               * p_patch(jg)%edges%dual_edge_length(iidx(jc,jb,2),iblk(jc,jb,2))), &
                  ABS(p_nh(jg)%metrics%ddxn_z_full(iidx(jc,jb,3),jk,iblk(jc,jb,3)) &
               * p_patch(jg)%edges%dual_edge_length(iidx(jc,jb,3),iblk(jc,jb,3)) ) )

            ! Exner extrapolation reaches zero for a slope of 1/4 or a height difference of 500 m
            ! between adjacent grid points (empirically determined values)
            p_nh(jg)%metrics%exner_exfac(jc,jk,jb) = exner_expol * &
              MIN(1._wp-(4._wp*z_maxslp(jc,jk,jb))**2._wp,         &
                  1._wp-(0.002_wp*z_maxhgtd(jc,jk,jb))**2._wp)
            p_nh(jg)%metrics%exner_exfac(jc,jk,jb) =  &
              MAX(0._vp,p_nh(jg)%metrics%exner_exfac(jc,jk,jb))
            ! For extremely steep slopes, going a bit behind time level nnow turned out
            ! to further improve stability
            IF (z_maxslp(jc,jk,jb) > 1.5_wp) &
              p_nh(jg)%metrics%exner_exfac(jc,jk,jb) = &
                MAX(-1._wp/6._wp,1._wp/9._wp*(1.5_wp-z_maxslp(jc,jk,jb)))
          ENDDO
        ENDDO

        jk = nlevp1
        DO jc = i_startidx, i_endidx

          z_maxslope = MAX(ABS(z_ddxn_z_half_e(iidx(jc,jb,1),jk,iblk(jc,jb,1))),&
                           ABS(z_ddxn_z_half_e(iidx(jc,jb,2),jk,iblk(jc,jb,2))),&
                           ABS(z_ddxn_z_half_e(iidx(jc,jb,3),jk,iblk(jc,jb,3))),&
                           ABS(z_ddxt_z_half_e(iidx(jc,jb,1),jk,iblk(jc,jb,1))),&
                           ABS(z_ddxt_z_half_e(iidx(jc,jb,2),jk,iblk(jc,jb,2))),&
                           ABS(z_ddxt_z_half_e(iidx(jc,jb,3),jk,iblk(jc,jb,3))) )

          z_diff = MAX(ABS(z_ddxn_z_half_e(iidx(jc,jb,1),jk,iblk(jc,jb,1))                   &
                         * p_patch(jg)%edges%dual_edge_length(iidx(jc,jb,1),iblk(jc,jb,1))), &
                       ABS(z_ddxn_z_half_e(iidx(jc,jb,2),jk,iblk(jc,jb,2))                   &
                         * p_patch(jg)%edges%dual_edge_length(iidx(jc,jb,2),iblk(jc,jb,2))), &
                       ABS(z_ddxn_z_half_e(iidx(jc,jb,3),jk,iblk(jc,jb,3))                   &
                         * p_patch(jg)%edges%dual_edge_length(iidx(jc,jb,3),iblk(jc,jb,3)) ) )

          ! Empirically determined values to ensure stability over steep slopes
          z_offctr =   MAX(vwind_offctr, 0.425_wp*z_maxslope**0.75_wp, MIN(0.25_wp,2.5e-4_wp*(z_diff-250._wp)))
          z_offctr =   MIN(MAX(vwind_offctr,0.75_wp),z_offctr)

          p_nh(jg)%metrics%vwind_impl_wgt(jc,jb)    = 0.5_wp + z_offctr
          p_nh(jg)%metrics%vwind_expl_wgt(jc,jb)    = 1.0_wp - p_nh(jg)%metrics%vwind_impl_wgt(jc,jb)

        ENDDO

        ! Search for grid points with particularly strong squeezing of the low-level coordinate levels over
        ! local peaks - these may also need increased offcentering in order to avoid sound-wave instabilities
        DO jk = MAX(10,nlev-8),nlev
          jk1 = jk + p_patch(jg)%nshift_total
          DO jc = i_startidx, i_endidx
            ! squeezing ratio with respect to nominal layer thickness
            z_diff = (p_nh(jg)%metrics%z_ifc(jc,jk,jb)-p_nh(jg)%metrics%z_ifc(jc,jk+1,jb))/(vct_a(jk1)-vct_a(jk1+1))
            IF (z_diff < 0.6_wp) THEN
              p_nh(jg)%metrics%vwind_impl_wgt(jc,jb)    = MAX(p_nh(jg)%metrics%vwind_impl_wgt(jc,jb),1.2_wp-z_diff)
              p_nh(jg)%metrics%vwind_expl_wgt(jc,jb)    = 1._wp - p_nh(jg)%metrics%vwind_impl_wgt(jc,jb)
            ENDIF
          ENDDO
        ENDDO

        ! Compute mask field for mountain or upper slope grid points
        IF (ASSOCIATED(ext_data(jg)%atm%sso_stdh_raw)) THEN  ! only associated for NWP forcing
          DO jc = i_startidx,i_endidx
            zn_avg = 0._wp
            zn_sq  = 0._wp
            zn_min = p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb)
            DO ic = 2, 10
              zn_avg = zn_avg + p_nh(jg)%metrics%z_ifc(inidx(ic,jc,jb),nlevp1,inblk(ic,jc,jb))
              zn_sq  = zn_sq  + p_nh(jg)%metrics%z_ifc(inidx(ic,jc,jb),nlevp1,inblk(ic,jc,jb))**2
              zn_min = MIN(zn_min,p_nh(jg)%metrics%z_ifc(inidx(ic,jc,jb),nlevp1,inblk(ic,jc,jb)))
            ENDDO
            zn_avg = zn_avg/9._wp
            zn_sq  = zn_sq/9._wp
            zn_rms = SQRT(MAX(0._wp,zn_sq - zn_avg**2))
            IF (p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb)-zn_min > 200._wp .AND. p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb) > zn_avg) THEN
              p_nh(jg)%metrics%mask_mtnpoints(jc,jb) = MIN(1._wp,1.e5_wp*MAX(0._wp,MAX(zn_rms,                               &
                p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb)-zn_avg)-200._wp)/p_patch(jg)%geometry_info%mean_characteristic_length**2)
            ELSE
              p_nh(jg)%metrics%mask_mtnpoints(jc,jb) = 0._wp
            ENDIF
            ! Second mask field used for gust parameterization
            IF (p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb)-zn_min > 50._wp .AND. p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb) > zn_avg) THEN
              z_aux_c(jc,1,jb) = MIN(1.2_wp,70._wp*MAX(zn_rms,p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb)-zn_avg,  &
                ext_data(jg)%atm%sso_stdh_raw(jc,jb))/MAX(6750._wp,p_patch(jg)%geometry_info%mean_characteristic_length))
            ELSE
              z_aux_c(jc,1,jb) = 0._wp
            ENDIF
          ENDDO
        ENDIF  ! associated

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


      DEALLOCATE(z_ddxn_z_half_e,z_ddxt_z_half_e)

      ALLOCATE(z_aux_e(nproma,1,nblks_c))
      z_aux_e(:,1,:) = 0._wp

      CALL sync_patch_array(SYNC_C, p_patch(jg), z_aux_c)
      CALL cell_avg(z_aux_c, p_patch(jg), p_int(jg)%c_bln_avg, z_aux_e)
      p_nh(jg)%metrics%mask_mtnpoints_g(:,:) = MIN(1._wp,z_aux_e(:,1,:))

      CALL sync_patch_array(SYNC_C, p_patch(jg), p_nh(jg)%metrics%mask_mtnpoints)
      z_aux_c(:,1,:) = p_nh(jg)%metrics%mask_mtnpoints(:,:)
      CALL cell_avg(z_aux_c, p_patch(jg), p_int(jg)%c_bln_avg, z_aux_e)
      p_nh(jg)%metrics%mask_mtnpoints(:,:) = z_aux_e(:,1,:)

      DEALLOCATE(z_aux_c,z_aux_c2,z_aux_e)

      ! Index lists for boundary nudging (including halo cells so that no
      ! sync is needed afterwards)
      ic = 0
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jc = 1, nlen
          IF (p_int(jg)%nudgecoeff_c(jc,jb) > 1.e-10_wp) THEN
            ic = ic+1
          ENDIF
        ENDDO
      ENDDO
      p_nh(jg)%metrics%nudge_c_dim = ic
      !$ACC UPDATE DEVICE(p_nh(jg)%metrics%nudge_c_dim) ASYNC(1)

      IF ( ic == 0 ) THEN
         ALLOCATE(p_nh(jg)%metrics%nudge_c_idx(0:0),p_nh(jg)%metrics%nudge_c_blk(0:0))
      ELSE
         ALLOCATE(p_nh(jg)%metrics%nudge_c_idx(ic),p_nh(jg)%metrics%nudge_c_blk(ic))
      ENDIF
      ic = 0
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jc = 1, nlen
          IF (p_int(jg)%nudgecoeff_c(jc,jb) > 1.e-10_wp) THEN
            ic = ic+1
            p_nh(jg)%metrics%nudge_c_idx(ic) = jc
            p_nh(jg)%metrics%nudge_c_blk(ic) = jb
          ENDIF
        ENDDO
      ENDDO

      !$ACC ENTER DATA COPYIN(p_nh(jg)%metrics%nudge_c_idx, p_nh(jg)%metrics%nudge_c_blk)

      ic = 0
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO je = 1, nlen
          IF (p_int(jg)%nudgecoeff_e(je,jb) > 1.e-10_wp) THEN
            ic = ic+1
          ENDIF
        ENDDO
      ENDDO
      p_nh(jg)%metrics%nudge_e_dim = ic
      !$ACC UPDATE DEVICE(p_nh(jg)%metrics%nudge_e_dim) ASYNC(1)

      IF ( ic == 0 ) THEN
         ALLOCATE(p_nh(jg)%metrics%nudge_e_idx(0:0),p_nh(jg)%metrics%nudge_e_blk(0:0))
      ELSE
         ALLOCATE(p_nh(jg)%metrics%nudge_e_idx(ic),p_nh(jg)%metrics%nudge_e_blk(ic))
      ENDIF

      ic = 0
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO je = 1, nlen
          IF (p_int(jg)%nudgecoeff_e(je,jb) > 1.e-10_wp) THEN
            ic = ic+1
            p_nh(jg)%metrics%nudge_e_idx(ic) = je
            p_nh(jg)%metrics%nudge_e_blk(ic) = jb
          ENDIF
        ENDDO
      ENDDO

      !$ACC ENTER DATA COPYIN(p_nh(jg)%metrics%nudge_e_idx, p_nh(jg)%metrics%nudge_e_blk)

      ! Index lists needed to minimize the number of halo communications in solve_nh and feedback

      p_nh(jg)%metrics%mask_prog_halo_c(:,:) = .FALSE.

      i_startblk = p_patch(jg)%cells%start_block(min_rlcell_int-1)
      i_endblk   = p_patch(jg)%cells%end_block(min_rlcell)

      ic = 0
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, min_rlcell_int-1, min_rlcell)

        DO jc = i_startidx, i_endidx
          IF (p_patch(jg)%cells%refin_ctrl(jc,jb)>=1 .AND. &
              p_patch(jg)%cells%refin_ctrl(jc,jb)<=4) THEN
            ic = ic+1
            p_nh(jg)%metrics%bdy_halo_c_idx(ic) = jc
            p_nh(jg)%metrics%bdy_halo_c_blk(ic) = jb
          ELSE
            p_nh(jg)%metrics%mask_prog_halo_c(jc,jb) = .TRUE.
          ENDIF
        ENDDO
      ENDDO

      !$ACC ENTER DATA COPYIN(p_nh(jg)%metrics%bdy_halo_c_idx, p_nh(jg)%metrics%bdy_halo_c_blk)

      ! Index list for which interpolated mass fluxes along the lateral nest boundary need to be updated
      ! part 1: count nest boundary points of row 9
      i_startblk = p_patch(jg)%edges%start_block(grf_bdywidth_e)
      i_endblk   = p_patch(jg)%edges%end_block(grf_bdywidth_e)

      ic = 0
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, grf_bdywidth_e, grf_bdywidth_e)

        DO je = i_startidx, i_endidx
          ic = ic+1
        ENDDO
      ENDDO

      ! part 2: count halo points of levels 1 and 2 belonging to nest boundary points of row 9
      i_startblk = p_patch(jg)%edges%start_block(min_rledge_int-1)
      i_endblk   = p_patch(jg)%edges%end_block(min_rledge_int-2)

      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, min_rledge_int-1, min_rledge_int-2)

        DO je = i_startidx, i_endidx
          IF (p_patch(jg)%edges%refin_ctrl(je,jb) == grf_bdywidth_e) THEN
            ic = ic+1
          ENDIF
        ENDDO
      ENDDO
      p_nh(jg)%metrics%bdy_mflx_e_dim = ic
      !$ACC UPDATE DEVICE(p_nh(jg)%metrics%bdy_mflx_e_dim) ASYNC(1)

      ! Allocate index lists and storage field for boundary mass flux
      IF ( ic == 0 ) THEN
         ALLOCATE(p_nh(jg)%metrics%bdy_mflx_e_idx(0:0),p_nh(jg)%metrics%bdy_mflx_e_blk(0:0), &
                  p_nh(jg)%diag%grf_bdy_mflx(nlev,0:0,2))
      ELSE
         ALLOCATE(p_nh(jg)%metrics%bdy_mflx_e_idx(ic),p_nh(jg)%metrics%bdy_mflx_e_blk(ic), &
                  p_nh(jg)%diag%grf_bdy_mflx(nlev,ic,2))
      ENDIF
      p_nh(jg)%diag%grf_bdy_mflx(:,:,:) = 0._wp

      !$ACC ENTER DATA COPYIN(p_nh(jg)%diag%grf_bdy_mflx)

      ! part 3: fill index list with nest boundary points of row 9
      i_startblk = p_patch(jg)%edges%start_block(grf_bdywidth_e)
      i_endblk   = p_patch(jg)%edges%end_block(grf_bdywidth_e)

      ic = 0
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, grf_bdywidth_e, grf_bdywidth_e)

        DO je = i_startidx, i_endidx
          ic = ic+1
          p_nh(jg)%metrics%bdy_mflx_e_idx(ic) = je
          p_nh(jg)%metrics%bdy_mflx_e_blk(ic) = jb
        ENDDO
      ENDDO

      ! part 4: fill index list with halo points of levels 1 and 2 belonging to nest boundary points of row 9
      i_startblk = p_patch(jg)%edges%start_block(min_rledge_int-1)
      i_endblk   = p_patch(jg)%edges%end_block(min_rledge_int-2)

      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, min_rledge_int-1, min_rledge_int-2)

        DO je = i_startidx, i_endidx
          IF (p_patch(jg)%edges%refin_ctrl(je,jb) == grf_bdywidth_e) THEN
            ic = ic+1
            p_nh(jg)%metrics%bdy_mflx_e_idx(ic) = je
            p_nh(jg)%metrics%bdy_mflx_e_blk(ic) = jb
          ENDIF
        ENDDO
      ENDDO

      !$ACC ENTER DATA COPYIN(p_nh(jg)%metrics%bdy_mflx_e_idx, p_nh(jg)%metrics%bdy_mflx_e_blk)

      ! Index list for halo points belonging to the nest overlap zone
      i_startblk = p_patch(jg)%cells%start_block(min_rlcell_int-1)
      i_endblk   = p_patch(jg)%cells%end_block(min_rlcell)
      ica(:)     = 0

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, min_rlcell_int-1, min_rlcell)

        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          DO jc = i_startidx, i_endidx
            IF (p_patch(jg)%cells%refin_ctrl(jc,jb)<=grf_fbk_start_c .AND. &
                p_patch(jg)%cells%child_id(jc,jb)==jgc) THEN
              ica(jn) = ica(jn)+1
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ic = MAXVAL(ica)
      jn = MAX(1,p_patch(jg)%n_childdom)
      IF ( ic == 0 ) THEN
         ALLOCATE(p_nh(jg)%metrics%ovlp_halo_c_dim(jn),   &
                  p_nh(jg)%metrics%ovlp_halo_c_idx(0:0,jn),&
                  p_nh(jg)%metrics%ovlp_halo_c_blk(0:0,jn))
      ELSE
         ALLOCATE(p_nh(jg)%metrics%ovlp_halo_c_dim(jn),   &
                  p_nh(jg)%metrics%ovlp_halo_c_idx(ic,jn),&
                  p_nh(jg)%metrics%ovlp_halo_c_blk(ic,jn))
      ENDIF

      p_nh(jg)%metrics%ovlp_halo_c_dim(1:jn) = ica(1:jn)

      ica(:) = 0
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, min_rlcell_int-1, min_rlcell)

        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          DO jc = i_startidx, i_endidx
            IF (p_patch(jg)%cells%refin_ctrl(jc,jb)<=grf_fbk_start_c .AND. &
                p_patch(jg)%cells%child_id(jc,jb)==jgc) THEN
              ica(jn) = ica(jn)+1
              p_nh(jg)%metrics%ovlp_halo_c_idx(ica(jn),jn) = jc
              p_nh(jg)%metrics%ovlp_halo_c_blk(ica(jn),jn) = jb

            ENDIF
          ENDDO
        ENDDO
      ENDDO

      !$ACC ENTER DATA COPYIN(p_nh(jg)%metrics%ovlp_halo_c_dim, p_nh(jg)%metrics%ovlp_halo_c_idx) &
      !$ACC   COPYIN(p_nh(jg)%metrics%ovlp_halo_c_blk)

      IF (l_zdiffu_t) THEN
        CALL prepare_zdiffu(p_patch(jg), p_nh(jg), p_nh_lists(jg)%metrics_list, p_int(jg), z_maxslp, z_maxhgtd)
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jc, z1, z2, z3) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jk = 2, nlev
          p_nh(jg)%metrics%wgtfac_c(1:nlen,jk,jb) =  &
           (p_nh(jg)%metrics%z_ifc(1:nlen,jk-1,jb) - &
            p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb)) /  &
           (p_nh(jg)%metrics%z_ifc(1:nlen,jk-1,jb) - &
            p_nh(jg)%metrics%z_ifc(1:nlen,jk+1,jb))
        ENDDO
          p_nh(jg)%metrics%wgtfac_c(1:nlen,1,jb) = &
           (p_nh(jg)%metrics%z_ifc(1:nlen,2,jb) -  &
            p_nh(jg)%metrics%z_ifc(1:nlen,1,jb)) / &
           (p_nh(jg)%metrics%z_ifc(1:nlen,3,jb) -  &
            p_nh(jg)%metrics%z_ifc(1:nlen,1,jb))
          p_nh(jg)%metrics%wgtfac_c(1:nlen,nlevp1,jb) = &
           (p_nh(jg)%metrics%z_ifc(1:nlen,nlev,jb) -    &
            p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb)) / &
           (p_nh(jg)%metrics%z_ifc(1:nlen,nlev-1,jb) -  &
            p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb))
        DO jc = 1, nlen
          z1 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,nlev,jb) - &
                          p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb) )
          z2 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,nlev,jb) +     &
                          p_nh(jg)%metrics%z_ifc(jc,nlev-1,jb) ) - &
                          p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb)
          z3 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,nlev-1,jb) +   &
                          p_nh(jg)%metrics%z_ifc(jc,nlev-2,jb) ) - &
                          p_nh(jg)%metrics%z_ifc(jc,nlevp1,jb)
          p_nh(jg)%metrics%wgtfacq_c(jc,3,jb) = z1*z2/(z2-z3)/(z1-z3)
          p_nh(jg)%metrics%wgtfacq_c(jc,2,jb) = &
            (z1-p_nh(jg)%metrics%wgtfacq_c(jc,3,jb)*(z1-z3))/(z1-z2)
          p_nh(jg)%metrics%wgtfacq_c(jc,1,jb) = 1._wp - &
           (p_nh(jg)%metrics%wgtfacq_c(jc,2,jb) +       &
            p_nh(jg)%metrics%wgtfacq_c(jc,3,jb))
        ENDDO
        DO jc = 1, nlen
          z1 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,2,jb) -   &
                          p_nh(jg)%metrics%z_ifc(jc,1,jb) )
          z2 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,2,jb) +   &
                          p_nh(jg)%metrics%z_ifc(jc,3,jb) ) - &
                          p_nh(jg)%metrics%z_ifc(jc,1,jb)
          z3 = 0.5_wp * ( p_nh(jg)%metrics%z_ifc(jc,3,jb) +   &
                          p_nh(jg)%metrics%z_ifc(jc,4,jb) ) - &
                          p_nh(jg)%metrics%z_ifc(jc,1,jb)
          p_nh(jg)%metrics%wgtfacq1_c(jc,3,jb) = z1*z2/(z2-z3)/(z1-z3)
          p_nh(jg)%metrics%wgtfacq1_c(jc,2,jb) = &
            (z1-p_nh(jg)%metrics%wgtfacq1_c(jc,3,jb)*(z1-z3))/(z1-z2)
          p_nh(jg)%metrics%wgtfacq1_c(jc,1,jb) = 1._wp - &
           (p_nh(jg)%metrics%wgtfacq1_c(jc,2,jb) +       &
            p_nh(jg)%metrics%wgtfacq1_c(jc,3,jb))
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (msg_level >= 10) THEN
        z_offctr = MAXVAL(p_nh(jg)%metrics%vwind_impl_wgt,MASK=p_patch(jg)%cells%decomp_info%owner_mask(:,:))
        z_offctr = global_max(z_offctr) - 0.5_wp
        DO jk = 1, nlev
          WHERE(.NOT.p_patch(jg)%cells%decomp_info%owner_mask(:,:))
            z_maxslp(:,jk,:)  = -HUGE(0._wp)
            z_maxhgtd(:,jk,:) = -HUGE(0._wp)
          END WHERE
        ENDDO
        z_maxslope = MAXVAL(z_maxslp)
        z_maxslope = global_max(z_maxslope)
        z_maxhdiff = MAXVAL(z_maxhgtd)
        z_maxhdiff = global_max(z_maxhdiff)
        WRITE(message_text,'(a,f8.4)') 'Maximum vertical wind offcentering: ', z_offctr
        CALL message(routine, message_text)
        WRITE(message_text,'(a,f8.4)') 'Maximum slope: ', z_maxslope
        CALL message(routine, message_text)
        WRITE(message_text,'(a,f8.1)') 'Maximum height difference between adjacent points: ', &
          z_maxhdiff
        CALL message(routine, message_text)
      ENDIF

      DEALLOCATE (z_maxslp,z_maxhgtd)

      ALLOCATE(z_aux_c(nproma,nlevp1,nblks_c),z_aux_e(nproma,nlevp1,nblks_e))
      z_aux_e = 0._wp

      ! Interpolate weighting coefficients to edges
      z_aux_c(:,:,:) = p_nh(jg)%metrics%wgtfac_c(:,:,:) ! necessary because wgtfac* may be single precision

      CALL cells2edges_scalar(z_aux_c, p_patch(jg), p_int(jg)%c_lin_e, z_aux_e, 1, nlevp1)
      CALL sync_patch_array(SYNC_E,p_patch(jg),z_aux_e)

      p_nh(jg)%metrics%wgtfac_e(:,:,:) = z_aux_e(:,:,:)
      DEALLOCATE(z_aux_c,z_aux_e)

      ALLOCATE(z_aux_c(nproma,6,nblks_c),z_aux_e(nproma,6,nblks_e))
      z_aux_e = 0._wp

      z_aux_c(:,1:3,:) = p_nh(jg)%metrics%wgtfacq_c(:,1:3,:)
      z_aux_c(:,4:6,:) = p_nh(jg)%metrics%wgtfacq1_c(:,1:3,:)

      CALL cells2edges_scalar(z_aux_c, p_patch(jg), p_int(jg)%c_lin_e, z_aux_e, 1, 6)
      CALL sync_patch_array(SYNC_E,p_patch(jg),z_aux_e)

      p_nh(jg)%metrics%wgtfacq_e (:,1:3,:) = z_aux_e(:,1:3,:)
      p_nh(jg)%metrics%wgtfacq1_e(:,1:3,:) = z_aux_e(:,4:6,:)

      DEALLOCATE(z_aux_c,z_aux_e)


      ! Reference atmosphere fields
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, z_help, z_temp, z_aux1, z_aux2) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1,nblks_c
        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF

        ! Reference surface temperature
        p_nh(jg)%metrics%tsfc_ref(1:nlen,jb) = (t0sl_bg-del_t_bg)+del_t_bg*   &
          & EXP(-p_nh(jg)%metrics%z_ifc(1:nlen,nlevp1,jb)/h_scal_bg)

        DO jk = 1, nlev
          ! Reference pressure, full level mass points
          z_aux1(1:nlen) = p0sl_bg*EXP(-grav/rd*h_scal_bg/(t0sl_bg-del_t_bg) &
            &  *LOG((EXP(p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)/h_scal_bg)      &
            &  *(t0sl_bg-del_t_bg) +del_t_bg)/t0sl_bg))

          ! Reference Exner pressure, full level mass points
          p_nh(jg)%metrics%exner_ref_mc(1:nlen,jk,jb) = (z_aux1(1:nlen)/p0ref)**rd_o_cpd

          ! Reference temperature, full level mass points
          z_temp(1:nlen) = (t0sl_bg-del_t_bg)+del_t_bg*        &
            & EXP(-p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)/h_scal_bg)

          ! Reference density, full level mass points
          p_nh(jg)%metrics%rho_ref_mc(1:nlen,jk,jb) = z_aux1(1:nlen)/(rd*z_temp(1:nlen))

          ! Reference Potential temperature, full level mass points
          p_nh(jg)%metrics%theta_ref_mc(1:nlen,jk,jb) = z_temp(1:nlen) &
            & /p_nh(jg)%metrics%exner_ref_mc(1:nlen,jk,jb)
        ENDDO

        IF (igradp_method <= 3) THEN
          DO jk = 1, nlev
            ! First vertical derivative of reference Exner pressure, full level mass points,
            ! divided by theta_ref
            ! Note: for computational efficiency, this field is in addition divided by
            ! the vertical layer thickness
            p_nh(jg)%metrics%d2dexdz2_fac1_mc(1:nlen,jk,jb)   =             &
              & -grav/(cpd*p_nh(jg)%metrics%theta_ref_mc(1:nlen,jk,jb)**2)* &
              & p_nh(jg)%metrics%inv_ddqz_z_full(1:nlen,jk,jb)

            ! Vertical derivative of d_exner_dz/theta_ref, full level mass points
            p_nh(jg)%metrics%d2dexdz2_fac2_mc(1:nlen,jk,jb)   =                            &
              &  2._wp*grav/(cpd*p_nh(jg)%metrics%theta_ref_mc(1:nlen,jk,jb)**3)*(grav/cpd &
              & -del_t_bg/h_scal_bg*EXP(-p_nh(jg)%metrics%z_mc(1:nlen,jk,jb)/h_scal_bg))   &
              & /p_nh(jg)%metrics%exner_ref_mc(1:nlen,jk,jb)
          ENDDO
        ENDIF

        DO jk = 1, nlevp1
          ! Reference pressure, half level mass points
          z_aux1(1:nlen) = p0sl_bg*EXP(-grav/rd*h_scal_bg/(t0sl_bg-del_t_bg) &
            &  *LOG((EXP(p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb)/h_scal_bg)      &
            &  *(t0sl_bg-del_t_bg) +del_t_bg)/t0sl_bg))

          ! Reference Exner pressure, half level mass points
          z_help(1:nlen) = (z_aux1(1:nlen)/p0ref)**rd_o_cpd

          ! Reference temperature, half level mass points
          z_temp(1:nlen) = (t0sl_bg-del_t_bg)+del_t_bg*          &
            & EXP(-p_nh(jg)%metrics%z_ifc(1:nlen,jk,jb)/h_scal_bg)

          ! Reference density, half level mass points
          z_aux2(1:nlen) = z_aux1(1:nlen)/(rd*z_temp(1:nlen))

          ! Reference Potential temperature, half level mass points
          p_nh(jg)%metrics%theta_ref_ic(1:nlen,jk,jb) = z_temp(1:nlen)/z_help(1:nlen)

          ! First vertical derivative of reference Exner pressure, half level mass points
          p_nh(jg)%metrics%d_exner_dz_ref_ic(1:nlen,jk,jb)   =       &
            & -grav/cpd/p_nh(jg)%metrics%theta_ref_ic(1:nlen,jk,jb)

        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ALLOCATE(z_me(nproma,nlev,p_patch(jg)%nblks_e))
      IF (p_test_run) THEN
        z_me = 0._wp
      ENDIF

      ! Compute geometric height at edge points
      CALL cells2edges_scalar(p_nh(jg)%metrics%z_mc, p_patch(jg), &
             p_int(jg)%c_lin_e, z_me)

      CALL sync_patch_array(SYNC_E,p_patch(jg),z_me)

      i_startblk = p_patch(jg)%edges%start_block(2)

      ! Reference fields on edges
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, je, z_aux1, z_temp) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

            ! Reference pressure, full level edge points
            z_aux1(je) = p0sl_bg*EXP(-grav/rd*h_scal_bg/(t0sl_bg-del_t_bg) &
              &  *LOG((EXP(z_me(je,jk,jb)/h_scal_bg)      &
              &  *(t0sl_bg-del_t_bg) +del_t_bg)/t0sl_bg))

            ! Reference temperature, full level edge points
            z_temp(je) = (t0sl_bg-del_t_bg)+del_t_bg*        &
              & EXP(-z_me(je,jk,jb)/h_scal_bg)

            ! Reference density, full level edge points
            p_nh(jg)%metrics%rho_ref_me(je,jk,jb) = z_aux1(je)/(rd*z_temp(je))

            ! Reference Potential temperature, full level edge points
            p_nh(jg)%metrics%theta_ref_me(je,jk,jb) = z_temp(je)/((z_aux1(je)/p0ref)**rd_o_cpd)
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! Compute information needed for Taylor-expansion-based pressure gradient calculation
      IF (igradp_method == 1) THEN
        nflat_gradp(jg) = nlev
        DEALLOCATE(z_me)
        CYCLE
      ENDIF

      IF (nflatlev(jg) <= 2 .AND. igradp_method >= 4) THEN
        CALL finish(routine, 'flat_height must be at least 2 levels below top&
          & for igradp_method>3')
      ENDIF

      ALLOCATE(flat_idx(nproma,p_patch(jg)%nblks_e))
      flat_idx = nlev

      iidx => p_patch(jg)%edges%cell_idx
      iblk => p_patch(jg)%edges%cell_blk


!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

        IF (igradp_method <= 3) THEN

          p_nh(jg)%metrics%zdiff_gradp(:,:,:,jb) = 0._vp

          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              p_nh(jg)%metrics%vertidx_gradp(1:2,je,jk,jb) = jk
              p_nh(jg)%metrics%zdiff_gradp(1,je,jk,jb)     =  &
                z_me(je,jk,jb) - p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk,iblk(je,jb,1))
              p_nh(jg)%metrics%zdiff_gradp(2,je,jk,jb)     =  &
                z_me(je,jk,jb) - p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk,iblk(je,jb,2))
            ENDDO
          ENDDO
        ENDIF

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

            ! Compute the highest vertical index (counted from top to bottom) for which
            ! the edge point lies inside the cell box of the adjacent grid points
            IF (z_me(je,jk,jb) <= p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk,iblk(je,jb,1))   .AND. &
                z_me(je,jk,jb) >= p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk+1,iblk(je,jb,1)) .AND. &
                z_me(je,jk,jb) <= p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk,iblk(je,jb,2))   .AND. &
                z_me(je,jk,jb) >= p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk+1,iblk(je,jb,2))) THEN
              flat_idx(je,jb) = jk
            ENDIF

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! Compute global minimum of flat_idx

      ! Attention: Edges in the halo of flat_idx are set to "random" values.
      ! So we have to consider only inner edges (or to do a boundary exchange)
      ! when searching the minimum.
      ! Please note also that a patch may be completely empty on some processors,
      ! in this case MINVAL returns HUGE() so there is no special care needed.

      nflat_gradp(jg) = MINVAL(flat_idx(:,:), MASK=p_patch(jg)%edges%decomp_info%owner_mask(:,:))
      nflat_gradp(jg) = NINT(global_min(REAL(nflat_gradp(jg),wp)))

      ! Compute level indices of neighbor cells within which the local edge is located
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jk1, jk_start, je, l_found, lfound_all, &
!$OMP z0, z1, z2, z3, coef1, coef2, coef3, dn1, dn2, dn3, dn4, dn5, &
!$OMP dn6) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,nblks_e

        CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

        IF (igradp_method <= 3) THEN

          DO je = i_startidx, i_endidx

            jk_start = flat_idx(je,jb)
            DO jk = flat_idx(je,jb)+1,nlev
              DO jk1 = jk_start, nlev
                IF ( (jk1 == nlev) .OR. z_me(je,jk,jb) <=                         &
                    p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk1,iblk(je,jb,1)) .AND. &
                    z_me(je,jk,jb) >=                                             &
                    p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))) THEN

                  p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = jk1
                  p_nh(jg)%metrics%zdiff_gradp(1,je,jk,jb)   = z_me(je,jk,jb) -   &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1,iblk(je,jb,1))
                  jk_start = jk1
                  EXIT
                ENDIF
              ENDDO
            ENDDO
            jk_start = flat_idx(je,jb)
            DO jk = flat_idx(je,jb)+1,nlev
              DO jk1 = jk_start, nlev
                IF ( (jk1 == nlev) .OR. z_me(je,jk,jb) <=                         &
                    p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk1,iblk(je,jb,2)) .AND. &
                    z_me(je,jk,jb) >=                                             &
                    p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))) THEN

                  p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = jk1
                  p_nh(jg)%metrics%zdiff_gradp(2,je,jk,jb)   = z_me(je,jk,jb) -   &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1,iblk(je,jb,2))
                  jk_start = jk1
                  EXIT
                ENDIF
              ENDDO
            ENDDO

          ENDDO

        ELSE IF (igradp_method >= 4) THEN ! Coefficients for polynomial interpolation

          p_nh(jg)%metrics%coeff_gradp(:,:,:,jb) = 0._vp

          ! Workaround for MPI deadlock with cce 8.7.x
          IF (msg_level > 300) CALL message(routine, 'igradp_method>=4')

          jk_start = nflatlev(jg) - 1
          DO jk = nflatlev(jg),nlev
            l_found(:) = .FALSE.
            lfound_all = .FALSE.
            DO jk1 = jk_start, nlev - 2
              DO je = i_startidx, i_endidx
                IF (z_me(je,jk,jb) <=                                            &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1,iblk(je,jb,1)) .AND. &
                    z_me(je,jk,jb) >                                             &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))) THEN

                  ! cubic interpolation
                  p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = jk1

                  z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1-1,iblk(je,jb,1))
                  z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1  ,iblk(je,jb,1))
                  z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))
                  z3 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+2,iblk(je,jb,1))

                  coef1 =  z_me(je,jk,jb)-z0
                  coef2 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)
                  coef3 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)*(z_me(je,jk,jb)-z2)

                  dn1 = 1._wp/(z0-z1)
                  dn2 = 1._wp/(z1-z2)
                  dn3 = 1._wp/(z2-z3)
                  dn4 = 1._wp/(z0-z2)
                  dn5 = 1._wp/(z0-z3)
                  dn6 = 1._wp/(z1-z3)

                  p_nh(jg)%metrics%coeff_gradp(1,je,jk,jb) =            &
                    1._wp + coef1*dn1 + coef2*dn1*dn4 + coef3*dn1*dn4*dn5
                  p_nh(jg)%metrics%coeff_gradp(2,je,jk,jb) =                               &
                    -(coef1*dn1 + coef2*dn4*(dn1+dn2) + coef3*dn5*(dn1*dn4+dn2*dn4+dn2*dn6))
                  p_nh(jg)%metrics%coeff_gradp(3,je,jk,jb) =          &
                    coef2*dn2*dn4 + coef3*dn5*(dn2*dn4+dn2*dn6+dn3*dn6)
                  p_nh(jg)%metrics%coeff_gradp(4,je,jk,jb) = -coef3*dn3*dn5*dn6

                  l_found(je) = .TRUE.
                ELSE IF (z_me(je,jk,jb) <=                                        &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-1,iblk(je,jb,1))) THEN

                  ! quadratic interpolation/extrapolation
                  p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = nlev - 1

                  z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-2,iblk(je,jb,1))
                  z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-1,iblk(je,jb,1))
                  z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev  ,iblk(je,jb,1))

                  coef1 =  z_me(je,jk,jb)-z0
                  coef2 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)

                  dn1 = 1._wp/(z0-z1)
                  dn2 = 1._wp/(z1-z2)
                  dn4 = 1._wp/(z0-z2)

                  p_nh(jg)%metrics%coeff_gradp(1,je,jk,jb) = 1._wp + coef1*dn1 + coef2*dn1*dn4
                  p_nh(jg)%metrics%coeff_gradp(2,je,jk,jb) = -(coef1*dn1 + coef2*dn4*(dn1+dn2))
                  p_nh(jg)%metrics%coeff_gradp(3,je,jk,jb) = coef2*dn2*dn4
                  p_nh(jg)%metrics%coeff_gradp(4,je,jk,jb) = 0._wp

                  l_found(je) = .TRUE.
                ENDIF
              ENDDO
              IF (ALL(l_found(i_startidx:i_endidx))) THEN
                lfound_all = .TRUE.
                EXIT
              ENDIF
            ENDDO
            IF (lfound_all) THEN
              jk_start = MIN(nlev-2, &
                MINVAL(p_nh(jg)%metrics%vertidx_gradp(1,i_startidx:i_endidx,jk,jb)))
            ENDIF
          ENDDO

          jk_start = nflatlev(jg) - 1
          DO jk = nflatlev(jg),nlev
            l_found(:) = .FALSE.
            lfound_all = .FALSE.
            DO jk1 = jk_start, nlev - 2
              DO je = i_startidx, i_endidx
                IF (z_me(je,jk,jb) <=                                            &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1,iblk(je,jb,2)) .AND. &
                    z_me(je,jk,jb) >                                             &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))) THEN

                  ! cubic interpolation
                  p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = jk1

                  z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1-1,iblk(je,jb,2))
                  z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1  ,iblk(je,jb,2))
                  z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))
                  z3 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+2,iblk(je,jb,2))

                  coef1 =  z_me(je,jk,jb)-z0
                  coef2 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)
                  coef3 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)*(z_me(je,jk,jb)-z2)

                  dn1 = 1._wp/(z0-z1)
                  dn2 = 1._wp/(z1-z2)
                  dn3 = 1._wp/(z2-z3)
                  dn4 = 1._wp/(z0-z2)
                  dn5 = 1._wp/(z0-z3)
                  dn6 = 1._wp/(z1-z3)

                  p_nh(jg)%metrics%coeff_gradp(5,je,jk,jb) =            &
                    1._wp + coef1*dn1 + coef2*dn1*dn4 + coef3*dn1*dn4*dn5
                  p_nh(jg)%metrics%coeff_gradp(6,je,jk,jb) =                               &
                    -(coef1*dn1 + coef2*dn4*(dn1+dn2) + coef3*dn5*(dn1*dn4+dn2*dn4+dn2*dn6))
                  p_nh(jg)%metrics%coeff_gradp(7,je,jk,jb) =          &
                    coef2*dn2*dn4 + coef3*dn5*(dn2*dn4+dn2*dn6+dn3*dn6)
                  p_nh(jg)%metrics%coeff_gradp(8,je,jk,jb) = -coef3*dn3*dn5*dn6

                  l_found(je) = .TRUE.
                ELSE IF (z_me(je,jk,jb) <=                                        &
                    p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-1,iblk(je,jb,2))) THEN

                  ! quadratic interpolation/extrapolation
                  p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = nlev - 1

                  z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-2,iblk(je,jb,2))
                  z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-1,iblk(je,jb,2))
                  z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev  ,iblk(je,jb,2))

                  coef1 =  z_me(je,jk,jb)-z0
                  coef2 = (z_me(je,jk,jb)-z0)*(z_me(je,jk,jb)-z1)

                  dn1 = 1._wp/(z0-z1)
                  dn2 = 1._wp/(z1-z2)
                  dn4 = 1._wp/(z0-z2)

                  p_nh(jg)%metrics%coeff_gradp(5,je,jk,jb) = 1._wp + coef1*dn1 + coef2*dn1*dn4
                  p_nh(jg)%metrics%coeff_gradp(6,je,jk,jb) = -(coef1*dn1 + coef2*dn4*(dn1+dn2))
                  p_nh(jg)%metrics%coeff_gradp(7,je,jk,jb) = coef2*dn2*dn4
                  p_nh(jg)%metrics%coeff_gradp(8,je,jk,jb) = 0._wp

                  l_found(je) = .TRUE.
                ENDIF
              ENDDO
              IF (ALL(l_found(i_startidx:i_endidx))) THEN
                lfound_all = .TRUE.
                EXIT
              ENDIF
            ENDDO
            IF (lfound_all) THEN
              jk_start = MIN(nlev-2, &
                MINVAL(p_nh(jg)%metrics%vertidx_gradp(2,i_startidx:i_endidx,jk,jb)))
            ENDIF
          ENDDO

        ENDIF

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (igradp_method == 3 .OR. igradp_method == 5) THEN

        i_startblk = p_patch(jg)%edges%start_block(grf_bdywidth_e+1)

        ALLOCATE(imask(nproma,nlev,p_patch(jg)%nblks_e),icount(p_patch(jg)%nblks_e), &
                 z_shift(nproma,nlev,p_patch(jg)%nblks_e) )

        extrapol_dist = 5._wp ! maximum allowed extrapolation distance; may become a namelist variable later on

        ! Recompute indices and height differences if truly horizontal pressure gradient
        ! computation would intersect the ground
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jk1, jk_start, je, z_aux1, z_aux2, &
!$OMP z0, z1, z2, z3, coef1, coef2, coef3, dn1, dn2, dn3, dn4, dn5, dn6) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,nblks_e

          CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                             i_startidx, i_endidx, grf_bdywidth_e+1)

          imask(:,:,jb) = 0
          icount(jb)    = 0

          DO je = i_startidx, i_endidx

            z_aux1(je) = MAX(p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),nlevp1,iblk(je,jb,1)), &
                             p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),nlevp1,iblk(je,jb,2)))

            z_aux2(je) = z_aux1(je) - extrapol_dist ! allow for some limited downward extrapolation

            DO jk = flat_idx(je,jb)+1,nlev
              IF ( z_me(je,jk,jb) < z_aux2(je)) THEN

                ! Save information needed for index list setup
                IF (p_patch(jg)%edges%decomp_info%owner_mask(je,jb)) THEN
                  imask(je,jk,jb) = 1
                  icount(jb)      = icount(jb) + 1
                  z_shift(je,jk,jb) = z_me(je,jk,jb) - z_aux2(je)
                ENDIF
              ENDIF
            ENDDO

            IF (igradp_method == 3) THEN
              jk_start = flat_idx(je,jb)
              DO jk = flat_idx(je,jb)+1,nlev
                IF ( z_me(je,jk,jb) < z_aux2(je)) THEN
                  DO jk1 = jk_start, nlev
                    IF ( jk1 == nlev .OR. z_aux2(je) <=                             &
                      p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk1,iblk(je,jb,1)) .AND. &
                      z_aux2(je) >=                                                 &
                      p_nh(jg)%metrics%z_ifc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))) THEN

                      p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = jk1
                      p_nh(jg)%metrics%zdiff_gradp(1,je,jk,jb)   = z_aux2(je) -     &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1,iblk(je,jb,1))
                      jk_start = jk1
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
              jk_start = flat_idx(je,jb)
              DO jk = flat_idx(je,jb)+1,nlev
                IF ( z_me(je,jk,jb) < z_aux2(je)) THEN
                  DO jk1 = jk_start, nlev
                    IF ( jk1 == nlev .OR. z_aux2(je) <=                             &
                      p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk1,iblk(je,jb,2)) .AND. &
                      z_aux2(je) >=                                                 &
                      p_nh(jg)%metrics%z_ifc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))) THEN

                      p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = jk1
                      p_nh(jg)%metrics%zdiff_gradp(2,je,jk,jb)   = z_aux2(je) -     &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1,iblk(je,jb,2))
                      jk_start = jk1
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ELSE IF (igradp_method == 5) THEN
              jk_start = flat_idx(je,jb)
              DO jk = flat_idx(je,jb)+1,nlev
                IF ( z_me(je,jk,jb) < z_aux2(je)) THEN
                  DO jk1 = jk_start, nlev - 2
                    IF (z_aux2(je) <=                                                &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1,iblk(je,jb,1)) .AND. &
                        z_aux2(je) >                                                 &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))) THEN

                      ! cubic interpolation
                      p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = jk1

                      z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1-1,iblk(je,jb,1))
                      z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1  ,iblk(je,jb,1))
                      z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+1,iblk(je,jb,1))
                      z3 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),jk1+2,iblk(je,jb,1))

                      coef1 =  z_aux2(je)-z0
                      coef2 = (z_aux2(je)-z0)*(z_aux2(je)-z1)
                      coef3 = (z_aux2(je)-z0)*(z_aux2(je)-z1)*(z_aux2(je)-z2)

                      dn1 = 1._wp/(z0-z1)
                      dn2 = 1._wp/(z1-z2)
                      dn3 = 1._wp/(z2-z3)
                      dn4 = 1._wp/(z0-z2)
                      dn5 = 1._wp/(z0-z3)
                      dn6 = 1._wp/(z1-z3)

                      p_nh(jg)%metrics%coeff_gradp(1,je,jk,jb) =            &
                        1._wp + coef1*dn1 + coef2*dn1*dn4 + coef3*dn1*dn4*dn5
                      p_nh(jg)%metrics%coeff_gradp(2,je,jk,jb) =                               &
                        -(coef1*dn1 + coef2*dn4*(dn1+dn2) + coef3*dn5*(dn1*dn4+dn2*dn4+dn2*dn6))
                      p_nh(jg)%metrics%coeff_gradp(3,je,jk,jb) =          &
                        coef2*dn2*dn4 + coef3*dn5*(dn2*dn4+dn2*dn6+dn3*dn6)
                      p_nh(jg)%metrics%coeff_gradp(4,je,jk,jb) = -coef3*dn3*dn5*dn6

                      jk_start = jk1
                      EXIT

                    ELSE IF (z_aux2(je) <=                                        &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-1,iblk(je,jb,1))) THEN

                      ! quadratic interpolation/extrapolation
                      p_nh(jg)%metrics%vertidx_gradp(1,je,jk,jb) = nlev - 1

                      z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-2,iblk(je,jb,1))
                      z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev-1,iblk(je,jb,1))
                      z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,1),nlev  ,iblk(je,jb,1))

                      coef1 =  z_aux2(je)-z0
                      coef2 = (z_aux2(je)-z0)*(z_aux2(je)-z1)

                      dn1 = 1._wp/(z0-z1)
                      dn2 = 1._wp/(z1-z2)
                      dn4 = 1._wp/(z0-z2)

                      p_nh(jg)%metrics%coeff_gradp(1,je,jk,jb) = &
                        1._wp + coef1*dn1 + coef2*dn1*dn4
                      p_nh(jg)%metrics%coeff_gradp(2,je,jk,jb) = &
                        -(coef1*dn1 + coef2*dn4*(dn1+dn2))
                      p_nh(jg)%metrics%coeff_gradp(3,je,jk,jb) = coef2*dn2*dn4
                      p_nh(jg)%metrics%coeff_gradp(4,je,jk,jb) = 0._wp

                      jk_start = nlev - 2
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO

              jk_start = flat_idx(je,jb)
              DO jk = flat_idx(je,jb)+1,nlev
                IF ( z_me(je,jk,jb) < z_aux2(je)) THEN
                  DO jk1 = jk_start, nlev - 2
                    IF (z_aux2(je) <=                                                &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1,iblk(je,jb,2)) .AND. &
                        z_aux2(je) >                                                 &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))) THEN

                      ! cubic interpolation
                      p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = jk1

                      z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1-1,iblk(je,jb,2))
                      z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1  ,iblk(je,jb,2))
                      z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+1,iblk(je,jb,2))
                      z3 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),jk1+2,iblk(je,jb,2))

                      coef1 =  z_aux2(je)-z0
                      coef2 = (z_aux2(je)-z0)*(z_aux2(je)-z1)
                      coef3 = (z_aux2(je)-z0)*(z_aux2(je)-z1)*(z_aux2(je)-z2)

                      dn1 = 1._wp/(z0-z1)
                      dn2 = 1._wp/(z1-z2)
                      dn3 = 1._wp/(z2-z3)
                      dn4 = 1._wp/(z0-z2)
                      dn5 = 1._wp/(z0-z3)
                      dn6 = 1._wp/(z1-z3)

                      p_nh(jg)%metrics%coeff_gradp(5,je,jk,jb) =            &
                        1._wp + coef1*dn1 + coef2*dn1*dn4 + coef3*dn1*dn4*dn5
                      p_nh(jg)%metrics%coeff_gradp(6,je,jk,jb) =                               &
                        -(coef1*dn1 + coef2*dn4*(dn1+dn2) + coef3*dn5*(dn1*dn4+dn2*dn4+dn2*dn6))
                      p_nh(jg)%metrics%coeff_gradp(7,je,jk,jb) =          &
                        coef2*dn2*dn4 + coef3*dn5*(dn2*dn4+dn2*dn6+dn3*dn6)
                      p_nh(jg)%metrics%coeff_gradp(8,je,jk,jb) = -coef3*dn3*dn5*dn6

                      jk_start = jk1
                      EXIT

                    ELSE IF (z_aux2(je) <=                                        &
                        p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-1,iblk(je,jb,2))) THEN

                      ! quadratic interpolation/extrapolation
                      p_nh(jg)%metrics%vertidx_gradp(2,je,jk,jb) = nlev - 1

                      z0 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-2,iblk(je,jb,2))
                      z1 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev-1,iblk(je,jb,2))
                      z2 = p_nh(jg)%metrics%z_mc(iidx(je,jb,2),nlev  ,iblk(je,jb,2))

                      coef1 =  z_aux2(je)-z0
                      coef2 = (z_aux2(je)-z0)*(z_aux2(je)-z1)

                      dn1 = 1._wp/(z0-z1)
                      dn2 = 1._wp/(z1-z2)
                      dn4 = 1._wp/(z0-z2)

                      p_nh(jg)%metrics%coeff_gradp(5,je,jk,jb) = &
                        1._wp + coef1*dn1 + coef2*dn1*dn4
                      p_nh(jg)%metrics%coeff_gradp(6,je,jk,jb) = &
                        -(coef1*dn1 + coef2*dn4*(dn1+dn2))
                      p_nh(jg)%metrics%coeff_gradp(7,je,jk,jb) = coef2*dn2*dn4
                      p_nh(jg)%metrics%coeff_gradp(8,je,jk,jb) = 0._wp

                      jk_start = nlev - 2
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO

            ENDIF
          ENDDO
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        ! Generate index list for grid points requiring downward extrapolation of the pressure gradient
        icount_total = SUM(icount(i_startblk:nblks_e))
        p_nh(jg)%metrics%pg_listdim = icount_total
        !$ACC UPDATE DEVICE(p_nh(jg)%metrics%pg_listdim) ASYNC(1)
        ic = 0

        ALLOCATE (p_nh(jg)%metrics%pg_edgeidx(icount_total),&
                  p_nh(jg)%metrics%pg_edgeblk(icount_total),&
                  p_nh(jg)%metrics%pg_vertidx(icount_total),&
                  p_nh(jg)%metrics%pg_exdist (icount_total))

        DO jb = i_startblk,nblks_e

          CALL get_indices_e(p_patch(jg), jb, i_startblk, nblks_e, &
                             i_startidx, i_endidx, grf_bdywidth_e+1)

          DO je = i_startidx, i_endidx
            DO jk = flat_idx(je,jb)+1,nlev
              IF (imask(je,jk,jb) == 1) THEN
                ic = ic + 1
                p_nh(jg)%metrics%pg_edgeidx(ic) = je
                p_nh(jg)%metrics%pg_edgeblk(ic) = jb
                p_nh(jg)%metrics%pg_vertidx(ic) = jk
                p_nh(jg)%metrics%pg_exdist (ic) = z_shift(je,jk,jb)
              ENDIF
            ENDDO
          ENDDO
        ENDDO

        !$ACC ENTER DATA COPYIN(p_nh(jg)%metrics%pg_edgeidx, p_nh(jg)%metrics%pg_edgeblk) &
        !$ACC   COPYIN(p_nh(jg)%metrics%pg_vertidx, p_nh(jg)%metrics%pg_exdist)

        DEALLOCATE(imask,icount,z_shift)

      ENDIF

      IF (igradp_method <= 3) THEN
        CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%zdiff_gradp(1,:,:,:))
        CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%zdiff_gradp(2,:,:,:))
      ELSE
        DO ic = 1, 8
          CALL sync_patch_array(SYNC_E,p_patch(jg),p_nh(jg)%metrics%coeff_gradp(ic,:,:,:))
        ENDDO
      ENDIF
      DEALLOCATE(z_me,flat_idx)

      ! deep-atmosphere modifications
      IF (ldeepatmo) THEN
        ! (It seems that aes_vdf_config is only initialized if iforcing=iaes or ildf_echam)
        IF (iforcing == iaes .OR. iforcing == ildf_echam) THEN
          is_les_phy = aes_vdf_config(jg)%turb==VDIFF_TURB_3DSMAGORINSKY
        ELSE
          is_les_phy = .FALSE.
        ENDIF
#ifndef __NO_ICON_LES__
        is_les_phy = atm_phy_nwp_config(jg)%is_les_phy .OR. is_les_phy
#endif
        IF (.NOT. ALLOCATED(z_me)) THEN 
          ALLOCATE(z_me(nproma,nlev,nblks_e), STAT=error_status)
          IF (error_status /= SUCCESS) CALL finish(routine, 'Allocation of z_me failed')
        ENDIF
        CALL cells2edges_scalar(p_nh(jg)%metrics%z_mc, p_patch(jg), p_int(jg)%c_lin_e, z_me)
        CALL sync_patch_array(SYNC_E, p_patch(jg), z_me)
        CALL prepare_deepatmo_metrics(nblks_c      = nblks_c,                  &
          &                           nblks_e      = nblks_e,                  &
          &                           npromz_c     = npromz_c,                 &
          &                           npromz_e     = npromz_e,                 &
          &                           nlev         = nlev,                     &
          &                           nlevp1       = nlevp1,                   &
          &                           nshift_total = p_patch(jg)%nshift_total, &
          &                           is_les_phy   = is_les_phy,               &
          &                           p_metrics    = p_nh(jg)%metrics,         &
          &                           z_me         = z_me,                     &
          &                           ierror       = error_status)
        IF (error_status /= SUCCESS) CALL finish(routine, 'prepare_deepatmo_metrics failed')
        DEALLOCATE(z_me, STAT=error_status)
        IF (error_status /= SUCCESS) CALL finish(routine, 'Deallocation of z_me failed')
      ENDIF

    ENDDO  !jg

    ! Now vgrid_buffer(:)%z_ifc can be deallocated
    DO jg = 1,n_dom
      DEALLOCATE(vgrid_buffer(jg)%z_ifc, STAT=error_status)
      IF (error_status /= SUCCESS) CALL finish (routine, 'DEALLOCATE(vgrid_buffer(jg)%z_ifc) failed.')
    END DO

    !PREPARE LES, Anurag Dipankar MPIM (2013-04)
    DO jg = 1 , n_dom
#ifndef __NO_ICON_LES__
      IF(atm_phy_nwp_config(jg)%is_les_phy .OR. aes_vdf_config(jg)%turb == VDIFF_TURB_3DSMAGORINSKY) THEN
#else
      IF(aes_vdf_config(jg)%turb == VDIFF_TURB_3DSMAGORINSKY) THEN
#endif
        CALL prepare_les_model(p_patch(jg), p_nh(jg), p_int(jg), jg)
      END IF
    END DO

    ! Prepare vertically varying nudging (only for primary domain)
    CALL prepare_nudging(p_patch(1), p_nh(1))
    !
    DO jg = 2, n_dom
      IF (nudging_config(jg)%nudge_type==indg_type%ubn) THEN
        ! if upper boundary nudging is activated for domain jg, 
        ! copy nudging coefficients from DOM 1 to child domain.
        !
        DO jk = 1, p_patch(jg)%nlev
          jk1 = jk + p_patch(jg)%nshift_total
          p_nh(jg)%metrics%nudgecoeff_vert(jk) = p_nh(1)%metrics%nudgecoeff_vert(jk1)
        ENDDO
        !
        IF (msg_level >= nudging_config(jg)%msg_thr%high) THEN
          WRITE(message_text,'(a,i2)') 'Nudging coefficients copied to DOM: ', jg
          CALL message(routine, message_text)
        ENDIF
      ENDIF 
    ENDDO  ! jg

  END SUBROUTINE set_nh_metrics
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !>
  !! Computation of coefficients and index lists needed for truly horizontal
  !! temperature diffusion.
  !!
  SUBROUTINE prepare_zdiffu(p_patch, p_nh, p_nh_metrics_list, p_int, maxslp, maxhgtd)

    CHARACTER(len=*), PARAMETER :: routine = modname//':prepare_zdiffu'

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch
    TYPE(t_nh_state), INTENT(INOUT)      :: p_nh
    TYPE(t_var_list_ptr), INTENT(INOUT)      :: p_nh_metrics_list
    TYPE(t_int_state), TARGET,INTENT(IN) :: p_int
    REAL(wp),        INTENT(IN)        :: maxslp(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp),        INTENT(IN)        :: maxhgtd(nproma,p_patch%nlev,p_patch%nblks_c)

    INTEGER :: jk, jb, jc, nblks_c, jk1, jk_start, ji, ji1
    INTEGER :: nlev                  !< number of full levels
    INTEGER :: i_startidx, i_endidx, i_startblk

    INTEGER,  DIMENSION(nproma,p_patch%nblks_c) :: i_masklist, k_start, &
                                                   k_end

    INTEGER, DIMENSION(nproma*p_patch%nblks_c) :: i_indlist, i_blklist
    INTEGER :: i_listdim, numpoints, i_listreduce(p_patch%nblks_c)

    REAL(wp) :: max_nbhgt, z_vintcoeff(nproma,p_patch%nlev,p_patch%nblks_c,3), &
                z_maxslp_avg(nproma,p_patch%nlev,p_patch%nblks_c), z_maxhgtd_avg(nproma,p_patch%nlev,p_patch%nblks_c)

    INTEGER  :: nbidx(nproma,p_patch%nlev,p_patch%nblks_c,3)

    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: i_listreduce,i_masklist,k_start,k_end,i_indlist,i_blklist
!DIR$ ATTRIBUTES ALIGN : 64 :: z_vintcoeff,z_maxslp_avg,z_maxhgtd_avg,nbidx
#endif

    nblks_c    = p_patch%nblks_c

    ! number of vertical levels
    nlev   = p_patch%nlev

    ji = 0
    i_masklist(:,:) = 0

    iidx => p_patch%cells%neighbor_idx
    iblk => p_patch%cells%neighbor_blk

    ! Apply cell averaging to slope and height difference fields
    i_startblk = p_patch%cells%start_block(2)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, i_startidx, i_endidx, 2)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          z_maxslp_avg(jc,jk,jb) = maxslp(jc,jk,jb)                      *p_int%c_bln_avg(jc,1,jb) + &
                                   maxslp(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_int%c_bln_avg(jc,2,jb) + &
                                   maxslp(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_int%c_bln_avg(jc,3,jb) + &
                                   maxslp(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_int%c_bln_avg(jc,4,jb)

          z_maxhgtd_avg(jc,jk,jb) = maxhgtd(jc,jk,jb)                      *p_int%c_bln_avg(jc,1,jb) + &
                                    maxhgtd(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_int%c_bln_avg(jc,2,jb) + &
                                    maxhgtd(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_int%c_bln_avg(jc,3,jb) + &
                                    maxhgtd(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_int%c_bln_avg(jc,4,jb)

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Horizontal index list for truly horizontal diffusion
    i_startblk = p_patch%cells%start_block(grf_bdywidth_c+1)

    ! Attention: this loop is not suitable for OpenMP parallelization
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, grf_bdywidth_c+1)

      DO jc = i_startidx, i_endidx
        IF ((z_maxslp_avg(jc,nlev,jb)>=thslp_zdiffu .OR. z_maxhgtd_avg(jc,nlev,jb)>=thhgtd_zdiffu) &
            .AND. p_patch%cells%decomp_info%owner_mask(jc,jb)) THEN
          ji = ji+1
          i_blklist(ji) = jb
          i_indlist(ji) = jc
          i_masklist(jc,jb) = 1
        ENDIF
      ENDDO
    ENDDO

    i_listdim = ji
    i_listreduce(:) = 0

!$OMP PARALLEL
    ! Vertical start and end indices for each grid point
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc, ji, max_nbhgt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, grf_bdywidth_c+1)

      DO jc = i_startidx, i_endidx
        IF (i_masklist(jc,jb) == 1) THEN
          ! Criterion for setting the end index: Stencil points must not
          ! intersect the ground
          max_nbhgt = MAX(p_nh%metrics%z_mc(iidx(jc,jb,1),nlev,iblk(jc,jb,1)), &
                          p_nh%metrics%z_mc(iidx(jc,jb,2),nlev,iblk(jc,jb,2)), &
                          p_nh%metrics%z_mc(iidx(jc,jb,3),nlev,iblk(jc,jb,3)))
          DO jk = nlev, 1, -1
            IF( p_nh%metrics%z_mc(jc,jk,jb) >= max_nbhgt) THEN
              k_end(jc,jb) = jk
              EXIT
            ENDIF
          ENDDO
          DO jk = 1, nlev
            IF(z_maxslp_avg(jc,jk,jb) >= thslp_zdiffu .OR. z_maxhgtd_avg(jc,jk,jb)>=thhgtd_zdiffu) THEN
              k_start(jc,jb) = jk
              EXIT
            ENDIF
          ENDDO
          ! Reset mask list if vertical index range turns out to be empty
          IF (k_start(jc,jb) > k_end(jc,jb)) THEN
            i_masklist(jc,jb) = 0
            i_listreduce(jb) = i_listreduce(jb) + 1
            k_start(jc,jb) = nlev
          ENDIF
        ENDIF
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    i_listdim = i_listdim - SUM(i_listreduce)

    ! Recompute index list after removal of empty points
    ji = 0
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, grf_bdywidth_c+1)

      DO jc = i_startidx, i_endidx
        IF (i_masklist(jc,jb) == 1) THEN
          ji = ji+1
          i_blklist(ji) = jb
          i_indlist(ji) = jc
        ENDIF
      ENDDO
    ENDDO

!$OMP PARALLEL

    ! Vertical indices for neighbor cells and related vertical interpolation coefficients
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc, jk1, jk_start) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, grf_bdywidth_c+1)

      DO jc = i_startidx, i_endidx
        IF (i_masklist(jc,jb) == 1) THEN
          jk_start = nlev - 1
          DO jk = k_end(jc,jb), k_start(jc,jb), -1
            DO jk1 = jk_start, 1, -1
              IF ( p_nh%metrics%z_mc(jc,jk,jb) <=                        &
                p_nh%metrics%z_mc(iidx(jc,jb,1),jk1,iblk(jc,jb,1)) .AND. &
                p_nh%metrics%z_mc(jc,jk,jb) >=                           &
                p_nh%metrics%z_mc(iidx(jc,jb,1),jk1+1,iblk(jc,jb,1))) THEN

                ! upper point for linear vertical interpolation
                nbidx(jc,jk,jb,1) = jk1
                ! interpolation weight for upper point
                z_vintcoeff(jc,jk,jb,1) = (p_nh%metrics%z_mc(jc,jk,jb) -  &
                  p_nh%metrics%z_mc(iidx(jc,jb,1),jk1+1,iblk(jc,jb,1))) / &
                  (p_nh%metrics%z_mc(iidx(jc,jb,1),jk1,iblk(jc,jb,1)) -   &
                  p_nh%metrics%z_mc(iidx(jc,jb,1),jk1+1,iblk(jc,jb,1)) )
                jk_start = jk1
                EXIT
              ENDIF
            ENDDO
          ENDDO
          jk_start = nlev - 1
          DO jk = k_end(jc,jb), k_start(jc,jb), -1
            DO jk1 = jk_start, 1, -1
              IF ( p_nh%metrics%z_mc(jc,jk,jb) <=                        &
                p_nh%metrics%z_mc(iidx(jc,jb,2),jk1,iblk(jc,jb,2)) .AND. &
                p_nh%metrics%z_mc(jc,jk,jb) >=                           &
                p_nh%metrics%z_mc(iidx(jc,jb,2),jk1+1,iblk(jc,jb,2))) THEN

                ! upper point for linear vertical interpolation
                nbidx(jc,jk,jb,2) = jk1
                ! interpolation weight for upper point
                z_vintcoeff(jc,jk,jb,2) = (p_nh%metrics%z_mc(jc,jk,jb) -  &
                  p_nh%metrics%z_mc(iidx(jc,jb,2),jk1+1,iblk(jc,jb,2))) / &
                  (p_nh%metrics%z_mc(iidx(jc,jb,2),jk1,iblk(jc,jb,2)) -   &
                  p_nh%metrics%z_mc(iidx(jc,jb,2),jk1+1,iblk(jc,jb,2)) )
                jk_start = jk1
                EXIT
              ENDIF
            ENDDO
          ENDDO
          jk_start = nlev - 1
          DO jk = k_end(jc,jb), k_start(jc,jb), -1
            DO jk1 = jk_start, 1, -1
              IF ( p_nh%metrics%z_mc(jc,jk,jb) <=                        &
                p_nh%metrics%z_mc(iidx(jc,jb,3),jk1,iblk(jc,jb,3)) .AND. &
                p_nh%metrics%z_mc(jc,jk,jb) >=                           &
                p_nh%metrics%z_mc(iidx(jc,jb,3),jk1+1,iblk(jc,jb,3))) THEN

                ! upper point for linear vertical interpolation
                nbidx(jc,jk,jb,3) = jk1
                ! interpolation weight for upper point
                z_vintcoeff(jc,jk,jb,3) = (p_nh%metrics%z_mc(jc,jk,jb) -  &
                  p_nh%metrics%z_mc(iidx(jc,jb,3),jk1+1,iblk(jc,jb,3))) / &
                  (p_nh%metrics%z_mc(iidx(jc,jb,3),jk1,iblk(jc,jb,3)) -   &
                  p_nh%metrics%z_mc(iidx(jc,jb,3),jk1+1,iblk(jc,jb,3)) )
                jk_start = jk1
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! Compute total number of grid points contained in the 3D index
    ! lists to know how much memory to allocate
    numpoints = 0
    DO ji = 1, i_listdim
      jc = i_indlist(ji)
      jb = i_blklist(ji)
      numpoints = numpoints+(k_end(jc,jb)-k_start(jc,jb)+1)
    ENDDO

    CALL new_zd_metrics(p_nh%metrics, p_nh_metrics_list , numpoints)
    
    p_nh%metrics%zd_listdim = numpoints
    !$ACC UPDATE DEVICE(p_nh%metrics%zd_listdim) ASYNC(1)

    ! Fill index lists
    ji1 = 0
    DO ji = 1, i_listdim
      jc = i_indlist(ji)
      jb = i_blklist(ji)
      DO jk = k_start(jc,jb),k_end(jc,jb)
        ji1 = ji1 + 1
        p_nh%metrics%zd_indlist(1,ji1) = jc
        p_nh%metrics%zd_indlist(2,ji1) = p_patch%cells%neighbor_idx(jc,jb,1)
        p_nh%metrics%zd_indlist(3,ji1) = p_patch%cells%neighbor_idx(jc,jb,2)
        p_nh%metrics%zd_indlist(4,ji1) = p_patch%cells%neighbor_idx(jc,jb,3)
        p_nh%metrics%zd_blklist(1,ji1) = jb
        p_nh%metrics%zd_blklist(2,ji1) = p_patch%cells%neighbor_blk(jc,jb,1)
        p_nh%metrics%zd_blklist(3,ji1) = p_patch%cells%neighbor_blk(jc,jb,2)
        p_nh%metrics%zd_blklist(4,ji1) = p_patch%cells%neighbor_blk(jc,jb,3)
        p_nh%metrics%zd_edgeidx(1,ji1) = p_patch%cells%edge_idx(jc,jb,1)
        p_nh%metrics%zd_edgeidx(2,ji1) = p_patch%cells%edge_idx(jc,jb,2)
        p_nh%metrics%zd_edgeidx(3,ji1) = p_patch%cells%edge_idx(jc,jb,3)
        p_nh%metrics%zd_edgeblk(1,ji1) = p_patch%cells%edge_blk(jc,jb,1)
        p_nh%metrics%zd_edgeblk(2,ji1) = p_patch%cells%edge_blk(jc,jb,2)
        p_nh%metrics%zd_edgeblk(3,ji1) = p_patch%cells%edge_blk(jc,jb,3)
        p_nh%metrics%zd_e2cell(1,ji1)  = p_int%e_bln_c_s(jc,1,jb)
        p_nh%metrics%zd_e2cell(2,ji1)  = p_int%e_bln_c_s(jc,2,jb)
        p_nh%metrics%zd_e2cell(3,ji1)  = p_int%e_bln_c_s(jc,3,jb)
        p_nh%metrics%zd_geofac (1,ji1) = p_int%geofac_n2s(jc,1,jb)
        p_nh%metrics%zd_geofac (2,ji1) = p_int%geofac_n2s(jc,2,jb)
        p_nh%metrics%zd_geofac (3,ji1) = p_int%geofac_n2s(jc,3,jb)
        p_nh%metrics%zd_geofac (4,ji1) = p_int%geofac_n2s(jc,4,jb)
        p_nh%metrics%zd_vertidx(1,ji1) = jk
        p_nh%metrics%zd_vertidx(2,ji1) = nbidx(jc,jk,jb,1)
        p_nh%metrics%zd_vertidx(3,ji1) = nbidx(jc,jk,jb,2)
        p_nh%metrics%zd_vertidx(4,ji1) = nbidx(jc,jk,jb,3)
        p_nh%metrics%zd_intcoef(1,ji1) = z_vintcoeff(jc,jk,jb,1)
        p_nh%metrics%zd_intcoef(2,ji1) = z_vintcoeff(jc,jk,jb,2)
        p_nh%metrics%zd_intcoef(3,ji1) = z_vintcoeff(jc,jk,jb,3)

        p_nh%metrics%zd_diffcoef(ji1)   =                                   &
          MAX(0._wp,SQRT(MAX(0._wp,z_maxslp_avg(jc,jk,jb)-thslp_zdiffu))/250._wp, &
                    2.e-4_wp*SQRT(MAX(0._wp,z_maxhgtd_avg(jc,jk,jb)-thhgtd_zdiffu)) )
        p_nh%metrics%zd_diffcoef(ji1)   = MIN(1._wp/500._wp,p_nh%metrics%zd_diffcoef(ji1))

      ENDDO
    ENDDO

    numpoints = global_sum_array(p_nh%metrics%zd_listdim)

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i10)') 'Number of z-diffusion points: ', numpoints
      CALL message(routine, message_text)
    ENDIF

  END SUBROUTINE prepare_zdiffu
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !>
  !! Computation of coefficients for LES model in mo_sgs_turbulence
  !! Computes the square of the mixing length.
  !!
  !! lambda^2 = (Cs * Delta)^2 *(kappa*x_3)^2 / ((Cs * Delta)^2 + (kappa*x_3)^2) 
  !!          = (Cs * Delta * x_3)^2 / ( (Cs * Delta / kappa)^2 + x_3^2 )
  !!
  !! with   Cs    : Smagorinsky constant
  !!        Delta : filter/grid width
  !!        x_3   : vertical coordinate
  !!        kappa : von Karman constant
  !!
  !! See Dipankar et al. (2015) "LES using ... ICON"
  !!
  SUBROUTINE prepare_les_model(p_patch, p_nh, p_int, jg)

    CHARACTER(len=*), PARAMETER :: routine = modname//':prepare_les_model'

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch
    TYPE(t_nh_state), INTENT(INOUT)      :: p_nh
    TYPE(t_int_state), TARGET,INTENT(IN) :: p_int
    INTEGER, INTENT(IN)                  :: jg

    REAL(wp)  :: les_filter, z_mc, z_aux(nproma,p_patch%nlevp1,p_patch%nblks_c)
    REAL(wp)  :: smag_constant, max_turb_scale

    INTEGER :: jk, jb, jc, je, nblks_c, nblks_e, nlen, i_startidx, i_endidx, npromz_c
    INTEGER :: nlev, nlevp1, i_startblk
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_aux
#endif

    IF ( aes_vdf_config(jg)%turb == VDIFF_TURB_3DSMAGORINSKY ) THEN
      smag_constant  = aes_vdf_config(jg)%smag_constant
      max_turb_scale = aes_vdf_config(jg)%max_turb_scale
    ELSE
      smag_constant  = les_config(jg)%smag_constant
      max_turb_scale = les_config(jg)%max_turb_scale
    END IF 

    nlev = p_patch%nlev
    nlevp1 = nlev + 1
    nblks_c   = p_patch%nblks_c
    npromz_c  = p_patch%npromz_c
    nblks_e   = p_patch%nblks_e

    i_startblk = p_patch%edges%start_block(2)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,nblks_e
     CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)
      DO jk = 1 , nlev
       DO je = i_startidx, i_endidx
         p_nh%metrics%inv_ddqz_z_full_e(je,jk,jb) =  &
                1._wp / p_nh%metrics%ddqz_z_full_e(je,jk,jb)
       END DO
      END DO
    END DO
!$OMP END DO

!$OMP DO PRIVATE(jb,jc,jk,nlen,z_mc,les_filter) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1,nblks_c
      nlen = MERGE(nproma, npromz_c, jb /= nblks_c)
      DO jk = 1 , nlevp1
        DO jc = 1 , nlen
          z_mc  = p_nh%metrics%geopot_agl_ifc(jc,jk,jb) * rgrav

          les_filter = smag_constant * MIN( max_turb_scale, &
                      (p_nh%metrics%ddqz_z_half(jc,jk,jb)*p_patch%cells%area(jc,jb))**0.33333_wp )
          
          p_nh%metrics%mixing_length_sq(jc,jk,jb) = (les_filter*z_mc)**2    &
               / ((les_filter/akt)**2+z_mc**2)

          p_nh%metrics%inv_ddqz_z_half(jc,jk,jb) = 1._wp / p_nh%metrics%ddqz_z_half(jc,jk,jb)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF(p_test_run)THEN
!$OMP PARALLEL
      CALL init(p_nh%metrics%inv_ddqz_z_half_v(:,:,:), lacc=.FALSE.)
      CALL init(p_nh%metrics%inv_ddqz_z_half_e(:,:,:), lacc=.FALSE.)
      CALL init(p_nh%metrics%wgtfac_v(:,:,:), lacc=.FALSE.)
!$OMP END PARALLEL
    END IF

   CALL cells2verts_scalar(p_nh%metrics%inv_ddqz_z_half, p_patch, p_int%cells_aw_verts, &
                           p_nh%metrics%inv_ddqz_z_half_v, opt_rlend=min_rlvert_int)

   z_aux(:,:,:) = p_nh%metrics%wgtfac_c(:,:,:) ! needed because wgtfac_c may be single precision
   CALL cells2verts_scalar(z_aux, p_patch, p_int%cells_aw_verts,          &
                           p_nh%metrics%wgtfac_v, opt_rlend=min_rlvert_int)

   CALL sync_patch_array_mult(SYNC_V,p_patch,2,p_nh%metrics%wgtfac_v,       &
                                 p_nh%metrics%inv_ddqz_z_half_v)

   CALL cells2edges_scalar(p_nh%metrics%inv_ddqz_z_half, p_patch, p_int%c_lin_e, &
                           p_nh%metrics%inv_ddqz_z_half_e, opt_rlend=min_rledge_int)


  END SUBROUTINE prepare_les_model
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !>
  !! Computation of nudging coefficient for nudging types:
  !! - Upper boundary nudging
  !! - Global nudging
  !!
  SUBROUTINE prepare_nudging(p_patch, p_nh)

    ! In/out variables
    TYPE(t_patch),    INTENT(IN)    :: p_patch
    TYPE(t_nh_state), INTENT(INOUT) :: p_nh

    ! Local variables
    TYPE(t_table) :: table
    REAL(wp) :: scale_height, distance, distance_scaled
    REAL(wp) :: start_height, end_height
    REAL(wp) :: height, nudge_coeff
    INTEGER  :: jg, jk, irow
    INTEGER  :: istart, iend
    INTEGER  :: istat
    INTEGER  :: nexp
    CHARACTER(LEN=16), PARAMETER :: column_jk     = "Full level index"
    CHARACTER(LEN=10), PARAMETER :: column_height = "Height (m)"
    CHARACTER(LEN=35), PARAMETER :: column_coeff  = "Nudging coefficient/max_nudge_coeff"
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':prepare_nudging'

    !----------------------------------------------------

    ! Note: here, we compute the dimensionless vertical nudging coefficient profile.  
    ! The actual max. nudging coefficients 'max_nudge_coeff_vn', 'max_nudge_coeff_thermdyn', ...
    ! from the nudging namelist will be multiplied with the profile during runtime.

    jg = p_patch%id

    IF (.NOT. nudging_config(jg)%lconfigured) THEN
      ! ('lconfigured' should have been set to .true. 
      ! in 'src/configure_model/mo_nudging_config: configure_nudging'
      ! or in 'src/namelists/mo_nudging_nml: check_nudging')
      CALL finish(routine, "Configuration of nudging_config still pending. "// &
        & "Please, check the program sequence.")
    ELSEIF ((.NOT. nudging_config(jg)%lnudging) .OR. jg /= 1) THEN
      ! The following computations have to be done only, 
      ! if (upper boundary) nudging is switched on, 
      ! and only for the primary domain
      RETURN
    ENDIF
    

    ! Start and end indices for vertical loop
    istart = nudging_config(jg)%ilev_start
    iend   = nudging_config(jg)%ilev_end

    ! Start and end height of vertical nudging region
    start_height = nudging_config(jg)%nudge_start_height
    end_height   = nudging_config(jg)%nudge_end_height
    
    ! Scale height to control, how fast the nudging strength decreases 
    ! with increasing vertical distance from nudging end height
    scale_height = ABS(nudging_config(jg)%nudge_scale_height)
    
    ! Discriminate between the different profiles of the nudging strength/nudging coefficient
    SELECT CASE(nudging_config(jg)%nudge_profile)
    CASE(indg_profile%sqrddist)
      ! Inverse squared scaled vertical distance from nudging start height
      DO jk = istart, iend
        ! Distance from nudging start height
        distance = 0.5_wp*ABS(vct_a(jk)+vct_a(jk+1)) - start_height
        ! Scaled distance
        distance_scaled                  = distance / MAX(1.0e-6_wp, scale_height)
        p_nh%metrics%nudgecoeff_vert(jk) = distance_scaled**2
      ENDDO  !jk
    CASE(indg_profile%const) ! For global nudging
      ! Constant profile
      DO jk = istart, iend
        p_nh%metrics%nudgecoeff_vert(jk) = 1._wp
      ENDDO  !jk
    CASE(indg_profile%tanh) ! For global nudging
      ! Hyperbolic tangent profile, which decreases in magnitude from nudging end height downwards
      ! (compare the coefficient profile for the Rayleigh damping based on Klemp et al. (2008) above)
      DO jk = istart, iend
        ! Distance from nudging end height
        distance = ABS(end_height - 0.5_wp * (vct_a(jk) + vct_a(jk+1)))
        ! Scaled distance
        distance_scaled                  = distance / MAX(1.0e-6_wp, scale_height)
        p_nh%metrics%nudgecoeff_vert(jk) = 1._wp - TANH(3.8_wp * distance_scaled)
      ENDDO  !jk
    CASE(indg_profile%trapid) ! For global nudging
      ! Trapezoidal profile
      !
      !                          1
      ! ... = ---------------------------------------,         (I)
      !       [     2 * z - (a + b)      ]2*nexp
      !       [--------------------------]       + 1
      !       [ (1 - 2 / nexp) * (b - a) ]
      ! 
      ! where z, a and b denote nominal height, start height and end height, respectively, 
      ! and the larger the parameter nexp, the more the trapezoidal profile compares 
      ! to a rectangular profile (i.e. it is inversely proportional to the scale height).
      ! 
      ! A lower bound of 3 for nexp is given by the factor "1 / (1 - 2 / nexp)" in formular (I). 
      ! An upper bound of "iend - istart + 1" (if > 3) corresponds to a scale height 
      ! of the order of the grid layer thickness
      nexp = MIN( MAX(3, NINT(ABS(end_height - start_height)/scale_height)), MAX(3, iend-istart+1) )
      DO jk = istart, iend
        ! "Distance"
        distance = vct_a(jk) + vct_a(jk+1) - ABS(start_height + end_height)
        ! Scaled "distance"
        distance_scaled                  = distance / &
          & MAX( 1.0e-20_wp, (1._wp - 2._wp / REAL(nexp, wp)) * ABS(end_height - start_height) )
        p_nh%metrics%nudgecoeff_vert(jk) = 1._wp / ( distance_scaled**(2*nexp) + 1._wp )
      ENDDO  !jk      
    END SELECT


    ! Print some info
    IF (msg_level >= nudging_config(jg)%msg_thr%high .AND. my_process_is_stdio()) THEN 
      ! Print the vertical profile of the nudging coefficient (nudging strength)
      WRITE(0,*) routine, ': Vertical profile of the nudging coefficient ', &
        & '(only levels, where it is non-zero):'
      ! Set up table
      CALL initialize_table(table)
      ! Set up table columns
      CALL add_table_column(table, column_jk)
      CALL add_table_column(table, column_height)
      CALL add_table_column(table, column_coeff)
      ! Fill the table rows
      DO jk = istart, iend
        height      = 0.5_wp * (vct_a(jk) + vct_a(jk+1))
        nudge_coeff = p_nh%metrics%nudgecoeff_vert(jk)
        irow        = jk-istart+1  ! (row index should start with 1)
        CALL set_table_entry(table, irow, column_jk, TRIM(int2string(jk)))
        CALL set_table_entry(table, irow, column_height, TRIM(real2string(height)))
        CALL set_table_entry(table, irow, column_coeff, TRIM(real2string(nudge_coeff)))
      ENDDO  !jk
      ! Print table
      CALL print_table(table)
      ! Destruct table
      CALL finalize_table(table)
    ENDIF  !IF (msg_level >= nudging_config(jg)%msg_thr%high .AND. my_process_is_stdio())
    
  END SUBROUTINE prepare_nudging
  !----------------------------------------------------------------------------
  ! Procedures related to the deep-atmosphere (deepatmo) configuration:
  !----------------------------------------------------------------------------
  !>
  !! Computes metrical modification for the deep atmosphere
  !!
  SUBROUTINE prepare_deepatmo_metrics( nblks_c,       &  !in
    &                                  nblks_e,       &  !in
    &                                  npromz_c,      &  !in
    &                                  npromz_e,      &  !in
    &                                  nlev,          &  !in
    &                                  nlevp1,        &  !in
    &                                  nshift_total,  &  !in
    &                                  is_les_phy,    &  !in
    &                                  p_metrics,     &  !inout
    &                                  z_me,          &  !in
    &                                  ierror)           !optout

    ! In/out variables
    INTEGER,  INTENT(IN)              :: nblks_c       ! Number of cell blocks
    INTEGER,  INTENT(IN)              :: nblks_e       ! Number of edge blocks
    INTEGER,  INTENT(IN)              :: npromz_c      ! Length of last cell block
    INTEGER,  INTENT(IN)              :: npromz_e      ! Length of last edge block
    INTEGER,  INTENT(IN)              :: nlev          ! Number of full levels
    INTEGER,  INTENT(IN)              :: nlevp1        ! Number of half levels
    INTEGER,  INTENT(IN)              :: nshift_total  ! Shift for vertical nesting
    LOGICAL,  INTENT(IN)              :: is_les_phy    ! For recomputation of fields (see set_nh_metrics above)
    TYPE(t_nh_metrics), INTENT(INOUT) :: p_metrics     ! Metrics
    REAL(wp), INTENT(IN)              :: z_me(:,:,:)   ! Height of edge center
    INTEGER,  INTENT(OUT), OPTIONAL   :: ierror        ! Optional error flag

    ! Local variables
    REAL(wp), ALLOCATABLE :: zgpot_me(:,:,:), zgpot_mc(:,:,:), zgpot_ifc(:,:,:)
    REAL(wp) :: help, temp, aux1, aux2
    REAL(wp) :: height_mc, height_uifc, height_lifc
    REAL(wp) :: radial_distance_mc, radial_distance_uifc, radial_distance_lifc
    INTEGER  :: jb, jk, jc, je, jk_shift
    INTEGER  :: nlen
    INTEGER  :: istat
    LOGICAL  :: is_ierror_present

    !-------------------------------------------------------------------------------

    ! Note: Quantities related to damping, diffusion and the like 
    ! are not modified for the deep atmosphere

    is_ierror_present = PRESENT(ierror)

    ! Initialize error flag with error:
    IF (is_ierror_present) ierror = 1

    ! Auxiliary fields
    ALLOCATE( zgpot_me(nproma,nlev,nblks_e),    &
      &       zgpot_mc(nproma,nlev,nblks_c),    &
      &       zgpot_ifc(nproma,nlevp1,nblks_c), &
      &       STAT=istat                        )
    IF (istat /= SUCCESS) RETURN

    !-----------------------------------------------------
    !                Geopotential heights
    !-----------------------------------------------------

    ! Transform (geometric) heights of levels into geopotential heights of levels

    CALL deepatmo_htrafo( z_in                = z_me,                        & !in
      &                   z_out               = zgpot_me,                    & !out
      &                   nblks_nproma_npromz = [nblks_e, nproma, npromz_e], & !in
      &                   start_end_levels    = [1, nlev],                   & !in
      &                   radius              = grid_sphere_radius,          & !in
      &                   trafo_type          = 'z2zgpot',                   & !in
      &                   ierror              = istat,                       & !optout
      &                   lacc                = .FALSE.)                       !optin
    IF (istat /= SUCCESS) RETURN

    CALL deepatmo_htrafo(z_in=p_metrics%z_mc, z_out=zgpot_mc, nblks_nproma_npromz=[nblks_c,nproma,npromz_c],  &
      & start_end_levels=[1,nlev], radius=grid_sphere_radius, trafo_type='z2zgpot', ierror=istat, lacc=.FALSE.)
    IF (istat /= SUCCESS) RETURN

    CALL deepatmo_htrafo(z_in=p_metrics%z_ifc, z_out=zgpot_ifc, nblks_nproma_npromz=[nblks_c,nproma,npromz_c],  &
      & start_end_levels=[1,nlevp1], radius=grid_sphere_radius, trafo_type='z2zgpot', ierror=istat, lacc=.FALSE.)
    IF (istat /= SUCCESS) RETURN

    !-----------------------------------------------------
    !                 Modification factors 
    !       for the deep-atmosphere configuration
    !-----------------------------------------------------

    ! Note: For the deep-atmosphere-specific metrical modification factors 
    ! and other quantities to be 1d-fields (varying only in z-direction), 
    ! the terrain-dependence of coordinate surfaces below z = 'flat_height' 
    ! is neglected in the deep-atmosphere modifications, 
    ! otherwise the following fields would become 3d-fields which would be 
    ! too costly in terms of memory. 
    ! In addition the metrical modification factors, 
    ! e.g. for flux divergences, are relatively difficult to compute 
    ! in sperical geometry, if coordinate surfaces deviate from spherical shells, 
    ! and cell edges lose the center of Earth as curvature center, so that their 
    ! shape is no longer determined by being great circle sections.

    DO jk = 1, nlev

      ! Preparation:
      ! Account for differences in nested domains' vertical index
      jk_shift = jk + nshift_total
      ! Height of upper half level of full level 'jk'
      height_uifc = vct_a(jk_shift)
      ! Height of lower half level of full level 'jk'
      height_lifc = vct_a(jk_shift+1)
      ! Height of full level
      height_mc = 0.5_wp * ( height_lifc + height_uifc )
      ! Corresponding radial distances from center of Earth
      radial_distance_mc   = grid_sphere_radius + height_mc
      radial_distance_uifc = grid_sphere_radius + height_uifc
      radial_distance_lifc = grid_sphere_radius + height_lifc

      ! Compute full-level metrical modification factors ...
      ! ... for horizontal gradients
      p_metrics%deepatmo_gradh_mc(jk) = grid_sphere_radius / radial_distance_mc
      ! ... for divergence: 
      ! Horizontal part (= surface of side wall / cell volume * flux denisty over side wall) 
      !                    ----------------------------------
      ! (-> modification is necessary for underlined factor)  
      ! (There is almost no difference between the magnitude of 'deepatmo_divh_mc'  
      ! and 'deepatmo_gradh_mc', but nevertheless they are not identical)
      p_metrics%deepatmo_divh_mc(jk) = p_metrics%deepatmo_gradh_mc(jk) * ( 3._wp / 4._wp ) / & 
        & ( 1._wp - radial_distance_lifc * radial_distance_uifc / ( radial_distance_lifc + radial_distance_uifc )**2 )
      ! Vertical part
      ! 1) = surface of cell bottom / cell volume * flux density over cell bottom
      !      ------------------------------------
      p_metrics%deepatmo_divzL_mc(jk) = 3._wp / ( 1._wp + radial_distance_uifc / radial_distance_lifc & 
        &                                       + ( radial_distance_uifc / radial_distance_lifc )**2 )
      ! 2) = surface of cell lid / cell volume * flux density over cell lid
      !      ---------------------------------
      p_metrics%deepatmo_divzU_mc(jk) = 3._wp / ( 1._wp + radial_distance_lifc / radial_distance_uifc &
        &                                       + ( radial_distance_lifc / radial_distance_uifc )**2 )
      ! Full-level metrical modification factors ...
      ! ... for the volume of a cell
      ! (This is required e.g. for volume integrals 
      ! in 'src/atm_dyn_iconam/mo_nh_supervise/subervise_total_integrals_nh')
      p_metrics%deepatmo_vol_mc(jk) = ( radial_distance_lifc**2 + radial_distance_lifc * radial_distance_uifc &
        &                             + radial_distance_uifc**2 ) / ( 3._wp * grid_sphere_radius**2 )
      ! ...: inverse of radial distance from center of Earth
      p_metrics%deepatmo_invr_mc(jk) = 1._wp / radial_distance_mc

      ! Compute half-level metrical modification factors ...
      ! ... for horizontal gradients
      p_metrics%deepatmo_gradh_ifc(jk) = grid_sphere_radius / radial_distance_uifc
      ! ...: inverse of radial distance from center of Earth
      p_metrics%deepatmo_invr_ifc(jk) = 1._wp / radial_distance_uifc

    ENDDO  !jk

    ! Half-level metrical modification factors at mean sea level (jk=nlevp1)
    p_metrics%deepatmo_gradh_ifc(nlevp1) = 1._wp  ! = grid_sphere_radius / radial_distance_uifc = grid_sphere_radius
    p_metrics%deepatmo_invr_ifc(nlevp1)  = 1._wp / grid_sphere_radius

    !-----------------------------------------------------
    !                   Recomputations 
    !      for inclusion of deep-atmosphere effects      
    !-----------------------------------------------------

    ! The computations in set_nh_metrics above are so densely compressed for runtime-efficiency,
    ! that a minimally invasive incorporation of deep-atmosphere modifications is practically impossible. 
    ! We have therefore no choice but to recompute a number of quantities here, unfortunately.

!$OMP PARALLEL

    ! Metrical modification of the geopotential

!$OMP DO PRIVATE(jb, jk, jc, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_c
        
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
        p_metrics%geopot(nlen+1:nproma,1:nlev,jb)           = 0.0_wp
        p_metrics%geopot_agl(nlen+1:nproma,1:nlev,jb)       = 0.0_wp
        p_metrics%geopot_agl_ifc(nlen+1:nproma,1:nlevp1,jb) = 0.0_wp
        p_metrics%dgeopot_mc(nlen+1:nproma,1:nlev,jb)       = 0.0_wp
      ENDIF
        
      DO jk = 1, nlev   
        DO jc = 1, nlen 
          ! Geopotential on full levels: 
          ! Phi = g * z / ( 1 + z / a ), provided that 
          ! Phi(r=a)=Phi(z=0)=0 m2/s2 (a = radius of Earth)
          p_metrics%geopot(jc,jk,jb) = grav * zgpot_mc(jc,jk,jb)
        ENDDO  !jc
      ENDDO  !jk
        
      ! Geopot above ground
      DO jc = 1, nlen 
        p_metrics%geopot_agl_ifc(jc,nlevp1,jb) = 0._wp
      ENDDO

      DO jk = nlev, 1, -1 
        DO jc = 1, nlen              
          ! Geopotential difference between layer interfaces            
          p_metrics%dgeopot_mc(jc,jk,jb) = grav * ( zgpot_ifc(jc,jk,jb) - zgpot_ifc(jc,jk+1,jb) )       
          ! Geopotential (interfaces)
          p_metrics%geopot_agl_ifc(jc,jk,jb) = p_metrics%geopot_agl_ifc(jc,jk+1,jb) + p_metrics%dgeopot_mc(jc,jk,jb)              
          ! Unfortunately, gpot[(z1+z2)/2] /= [gpot(z1)+gpot(z2)]/2, 
          ! where gpot(z)=zgpot, so we have to compute it separately for full levels
          p_metrics%geopot_agl(jc,jk,jb) = p_metrics%geopot_agl_ifc(jc,jk+1,jb) &
            &                            + grav * ( zgpot_mc(jc,jk,jb) - zgpot_ifc(jc,jk+1,jb) )      
        ENDDO  !jc            
      ENDDO  !jk 
        
    ENDDO  !jb
!$OMP END DO

    ! Metrical modification of interface slope

    ! Note: numerous quantities are determined from the slope at (or close to) the ground. 
    ! A modification of them is not necessary, since the respective modification factors 
    ! are terrain-independent and would be equal to 1 everywhere at the ground 
    ! (or at least negligibly small)

!$OMP DO PRIVATE(jb, jk, je, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_e

      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSEIF (is_les_phy) THEN
        nlen = npromz_e
        p_metrics%ddxn_z_half_e(nlen+1:nproma,1:nlevp1,jb) = 0.0_vp
        p_metrics%ddxt_z_half_e(nlen+1:nproma,1:nlevp1,jb) = 0.0_vp
        p_metrics%ddxn_z_full(nlen+1:nproma,1:nlev,jb)     = 0.0_vp
        p_metrics%ddxt_z_full(nlen+1:nproma,1:nlev,jb)     = 0.0_vp
      ELSE
        nlen = npromz_e
        p_metrics%ddxn_z_full(nlen+1:nproma,1:nlev,jb)     = 0.0_vp
        p_metrics%ddxt_z_full(nlen+1:nproma,1:nlev,jb)     = 0.0_vp
      ENDIF

      IF (is_les_phy) THEN
        DO jk = 1, nlevp1
          DO je = 1, nlen
            p_metrics%ddxn_z_half_e(je,jk,jb) = p_metrics%deepatmo_gradh_ifc(jk) * p_metrics%ddxn_z_half_e(je,jk,jb)
            p_metrics%ddxt_z_half_e(je,jk,jb) = p_metrics%deepatmo_gradh_ifc(jk) * p_metrics%ddxt_z_half_e(je,jk,jb)
          ENDDO  !je
        ENDDO  !jk
      ENDIF  !IF (is_les_phy)

      DO jk = 1, nlev
        DO je = 1, nlen
          p_metrics%ddxn_z_full(je,jk,jb) = p_metrics%deepatmo_gradh_mc(jk) * p_metrics%ddxn_z_full(je,jk,jb)          
          p_metrics%ddxt_z_full(je,jk,jb) = p_metrics%deepatmo_gradh_mc(jk) * p_metrics%ddxt_z_full(je,jk,jb)
        ENDDO  !je
      ENDDO  !jk

    ENDDO  !jb
!$OMP END DO

    ! Metrical modification of reference atmosphere

    ! The deep-atmosphere case is more or less a copy of the shallow-atmosphere case 
    ! (for computational efficiency reasons), with the geometric height z 
    ! replaced by the geopotential height zgpot (we say that the reference temperature 
    ! has the same functional dependency on zgpot as it has on z in case of the 
    ! shallow atmosphere, this way the integral of -dp/dzgpot-rho*grav=0 in the 
    ! deep-atmosphere case is formally identical to the integral of -dp/dz-rho*grav=0 
    ! in case of the shallow atmospehre). 
    ! Note: the vertical derivative in 'd2dexdz2_fac1_mc', 'd2dexdz2_fac2_mc', 
    ! and 'd_exner_dz_ref_ic' is still with respect to z, not zgpot.

!$OMP DO PRIVATE(jb, jk, jc, nlen, help, temp, aux1, aux2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_c
        
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
        p_metrics%tsfc_ref(nlen+1:nproma,jb)                   = 0.0_vp
        p_metrics%exner_ref_mc(nlen+1:nproma,1:nlev,jb)        = 0.0_vp
        p_metrics%rho_ref_mc(nlen+1:nproma,1:nlev,jb)          = 0.0_vp
        p_metrics%theta_ref_mc(nlen+1:nproma,1:nlev,jb)        = 0.0_vp
        p_metrics%theta_ref_ic(nlen+1:nproma,1:nlevp1,jb)      = 0.0_vp
        p_metrics%d_exner_dz_ref_ic(nlen+1:nproma,1:nlevp1,jb) = 0.0_vp
        IF (igradp_method <= 3) THEN
          p_metrics%d2dexdz2_fac1_mc(nlen+1:nproma,1:nlev,jb)  = 0.0_vp
          p_metrics%d2dexdz2_fac2_mc(nlen+1:nproma,1:nlev,jb)  = 0.0_vp
        ENDIF
      ENDIF
        
      ! Reference surface temperature
      DO jc = 1, nlen
        p_metrics%tsfc_ref(jc,jb) = ( t0sl_bg - del_t_bg ) + del_t_bg * EXP( -zgpot_ifc(jc,nlevp1,jb) / h_scal_bg )
      ENDDO  !jc
        
      DO jk = 1, nlev
        DO jc = 1, nlen
          ! Reference pressure, full level mass points
          aux1 = p0sl_bg  * EXP( -grav / rd * h_scal_bg / ( t0sl_bg - del_t_bg ) * &
            &    LOG( ( EXP( zgpot_mc(jc,jk,jb) / h_scal_bg ) * ( t0sl_bg - del_t_bg ) + del_t_bg ) / t0sl_bg ) )      
          ! Reference Exner pressure, full level mass points
          p_metrics%exner_ref_mc(jc,jk,jb) = ( aux1 / p0ref )**rd_o_cpd            
          ! Reference temperature, full level mass points
          temp = ( t0sl_bg - del_t_bg ) + del_t_bg * EXP( -zgpot_mc(jc,jk,jb) / h_scal_bg )            
          ! Reference density, full level mass points
          p_metrics%rho_ref_mc(jc,jk,jb) = aux1 / ( rd * temp )            
          ! Reference potential temperature, full level mass points
          p_metrics%theta_ref_mc(jc,jk,jb) = temp / p_metrics%exner_ref_mc(jc,jk,jb)
        ENDDO  !jc
      ENDDO  !jk
        
      IF (igradp_method <= 3) THEN
        DO jk = 1, nlev
          DO jc = 1, nlen
            ! First vertical derivative of reference Exner pressure, full level mass points,
            ! divided by theta_ref
            ! Note: for computational efficiency, this field is in addition divided by
            ! the vertical layer thickness
            ! (Deep-atmosphere modification: the vertical derivative is with respect to z, 
            ! so we use dexner/dz=dexner/dzgpot*dzgpot/dz, with dzgpot/dz=(a/r)**2, 
            ! where a is the radius of Earth, and r=a+z. For the actual modification, 
            ! we can use that (a/r)**2=deepatmo_gradh_mc**2)
            p_metrics%d2dexdz2_fac1_mc(jc,jk,jb) = -grav / ( cpd * p_metrics%theta_ref_mc(jc,jk,jb)**2 ) * &
              &                                    p_metrics%inv_ddqz_z_full(jc,jk,jb) * p_metrics%deepatmo_gradh_mc(jk)**2              
            ! Vertical derivative of d_exner_dz/theta_ref, full level mass points
            ! (Deep-atmosphere modification: here we use that for a quantity X, 
            ! d**2X/dz**2=d**2X/dzgpot**2*dzgpot/dz*(a/r)**2 + dX/dzgpot*d(a/r)**2/dr 
            ! =d**2X/dzgpot**2*(a/r)**4 - 2*dX/dzgpot*(a/r)**2/r)
            p_metrics%d2dexdz2_fac2_mc(jc,jk,jb) = 2._wp * grav / ( cpd * p_metrics%theta_ref_mc(jc,jk,jb)**3 ) * ( grav / cpd &
              &                        - del_t_bg / h_scal_bg * EXP( -zgpot_mc(jc,jk,jb) / h_scal_bg ) )                       &
              &                        / p_metrics%exner_ref_mc(jc,jk,jb) * p_metrics%deepatmo_gradh_mc(jk)**4                 & 
              &                        - 2._wp * p_metrics%d2dexdz2_fac1_mc(jc,jk,jb) * p_metrics%ddqz_z_full(jc,jk,jb)        &
              &                        * p_metrics%deepatmo_gradh_mc(jk) / grid_sphere_radius
          ENDDO  !jc
        ENDDO  !jk
      ENDIF  !IF (igradp_method <= 3)
        
      DO jk = 1, nlevp1
        DO jc = 1, nlen
          ! Reference pressure, half level mass points
          aux1 = p0sl_bg * EXP( -grav / rd * h_scal_bg / ( t0sl_bg - del_t_bg ) * &
            &    LOG( ( EXP( zgpot_ifc(jc,jk,jb) / h_scal_bg ) * ( t0sl_bg - del_t_bg ) + del_t_bg ) / t0sl_bg ) )            
          ! Reference Exner pressure, half level mass points
          help = ( aux1 / p0ref )**rd_o_cpd            
          ! Reference temperature, half level mass points
          temp = ( t0sl_bg - del_t_bg ) + del_t_bg * EXP( -zgpot_ifc(jc,jk,jb) / h_scal_bg )            
          ! Reference density, half level mass points
          aux2 = aux1 / ( rd * temp )            
          ! Reference Potential temperature, half level mass points
          p_metrics%theta_ref_ic(jc,jk,jb) = temp / help           
          ! First vertical derivative of reference Exner pressure, half level mass points
          p_metrics%d_exner_dz_ref_ic(jc,jk,jb) = -grav / cpd / p_metrics%theta_ref_ic(jc,jk,jb) &
            &                                   * p_metrics%deepatmo_gradh_ifc(jk)**2
        ENDDO  !jc
      ENDDO  !jk
        
    ENDDO  !jb
!$OMP END DO

!$OMP DO PRIVATE(jb, jk, je, nlen, aux1, temp) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_e

      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
        p_metrics%rho_ref_me(nlen+1:nproma,1:nlev,jb)   = 0.0_vp
        p_metrics%theta_ref_me(nlen+1:nproma,1:nlev,jb) = 0.0_vp
      ENDIF

      DO jk = 1, nlev
        DO je = 1, nlen            
          ! Reference pressure, full level edge points
          aux1 = p0sl_bg * EXP( -grav / rd * h_scal_bg / ( t0sl_bg - del_t_bg ) * &
            &    LOG( ( EXP( zgpot_me(je,jk,jb) / h_scal_bg ) * ( t0sl_bg - del_t_bg ) + del_t_bg ) / t0sl_bg ) )              
          ! Reference temperature, full level edge points
          temp = ( t0sl_bg - del_t_bg ) + del_t_bg * EXP( -zgpot_me(je,jk,jb) / h_scal_bg )
          ! Reference density, full level edge points
          p_metrics%rho_ref_me(je,jk,jb) = aux1 / ( rd * temp )            
          ! Reference potential temperature, full level edge points
          p_metrics%theta_ref_me(je,jk,jb) = temp / ( ( aux1 / p0ref )**rd_o_cpd )
        ENDDO  !je
      ENDDO  !jk

    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    ! Clean-up
    DEALLOCATE(zgpot_me, zgpot_mc, zgpot_ifc, STAT=istat)
    IF (istat /= SUCCESS) RETURN

    IF (is_ierror_present) ierror = SUCCESS

  END SUBROUTINE prepare_deepatmo_metrics

END MODULE mo_vertical_grid

