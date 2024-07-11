! Horizontal transport and refraction of spectral surface wave energy
!
! Main routine for the horizontal transport of surface wave energy.
! Here, we integrate the spectral energy equation in time without
! sources and sinks, only taking into account advection and refraction.
!
! For the advection part, we make use of the horizontal transport scheme for
! tracers of the atmospheric model.
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
!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_advection_stepping

  USE mo_kind,                      ONLY: wp, vp
  USE mo_impl_constants,            ONLY: min_rlcell_int, min_rledge_int
  USE mo_impl_constants_grf,        ONLY: grf_bdywidth_c
  USE mo_loopindices,               ONLY: get_indices_c
  USE mo_model_domain,              ONLY: t_patch
  USE mo_parallel_config,           ONLY: nproma
  USE mo_grid_config,               ONLY: l_limited_area
  USE mo_run_config,                ONLY: ntracer, ltimer
  USE mo_interpol_config,           ONLY: llsq_lin_consv
  USE mo_intp_data_strc,            ONLY: t_int_state
  USE mo_wave_config,               ONLY: t_wave_config
  USE mo_wave_refraction,           ONLY: wave_refraction
  USE mo_wave_physics,              ONLY: set_energy2emin
  USE mo_energy_propagation_config, ONLY: t_energy_propagation_config
  USE mo_advection_traj,            ONLY: btraj_compute_o1, t_back_traj
  USE mo_advection_hflux,           ONLY: upwind_hflux_miura
  USE mo_fortran_tools,             ONLY: init
  USE mo_sync,                      ONLY: SYNC_C, sync_patch_array_mult
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_transport

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wave_step_advection

CONTAINS


  SUBROUTINE wave_step_advection( p_patch, p_int_state, wave_config, energy_propagation_config, &
    &                             p_dtime, wave_num_c, gv_c, bathymetry_c, geo_bath_grad_c, &
    &                             p_mflx_h, p_vn_traj, p_vt_traj, p_tracer_now, p_tracer_new )

    TYPE(t_patch), TARGET,            INTENT(IN):: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state),                INTENT(IN):: &  !< interpolation state
      &  p_int_state

    TYPE(t_wave_config), TARGET,      INTENT(IN):: &  !< wave config state
      &  wave_config

    TYPE(t_energy_propagation_config), TARGET, INTENT(IN):: & !< energy propagation config state
      &  energy_propagation_config

    REAL(wp),                         INTENT(IN):: & !< advective time step [s]
      &  p_dtime

    REAL(wp),                         INTENT(IN):: & !< wave number at cell center [1/m]
      &  wave_num_c(:,:,:)                           !< dim: (nproma,nblks_c,nfreqs)

    REAL(wp),                         INTENT(IN):: & !< group velocity at cell center [m/s]
      &  gv_c(:,:,:)                                 !< dim: (nproma,nblks_c,nfreqs)

    REAL(wp),                         INTENT(IN):: & !< bathymetry at cell center [m]
      &  bathymetry_c(:,:)                           !< dim: (nproma,nblks_c)

    REAL(vp),                         INTENT(IN):: & !< gradient of bathymetry [m/m]
      &  geo_bath_grad_c(:,:,:,:)                    !< dim: (2,nproma,nlev,nblks_c)

    REAL(wp),                         INTENT(IN):: & !< horizontal mass flux at edge midpoints
      &  p_mflx_h(:,:,:,:)                           !< WAVE: gv_n  [m/s]
                                                     !< dim: (nproma,nlev,nblks_e,ntracer)

    REAL(wp),                         INTENT(IN):: & !< edge-normal horizontal group velocity component at n+1/2
      &  p_vn_traj(:,:,:,:)                          !< for the calculation of backward trajectories
                                                     !< [m/s]
                                                     !< dim: (nproma,nlev,nblks_e,ntracer)

   REAL(wp),                          INTENT(IN):: & !< edge-tangential horizontal group velocity component at n+1/2
      &  p_vt_traj(:,:,:,:)                          !< for the calculation of backward trajectories
                                                     !< [m/s]
                                                     !< dim: (nproma,nlev,nblks_e,ntracer)

    REAL(wp), CONTIGUOUS,            INTENT(INOUT):: & !< spectral wave energy
      &  p_tracer_now(:,:,:,:)                       !< at current time level n (before transport)
                                                     !< [kg/kg]
                                                     !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), CONTIGUOUS,            INTENT(INOUT) :: & !< spectral wave energy
      &  p_tracer_new(:,:,:,:)                          !< at time level n+1 (after transport)
                                                        !< [kg/kg]
                                                        !< dim: (nproma,nlev,nblks_c,ntracer)

    ! local
    INTEGER :: jb, jk, jc, jt
    INTEGER :: i_startidx, i_endidx
    INTEGER :: i_rlstart_c, i_rlend_c
    INTEGER :: i_rlstart_e, i_rlend_e
    INTEGER :: i_startblk_c, i_endblk_c
    INTEGER :: i_rlstart_bdy, i_rlend_bdy
    INTEGER :: i_startblk_bdy, i_endblk_bdy

    TYPE(t_energy_propagation_config), POINTER :: &
      &  enprop_conf                                      !< convenience pointer to save fome paperwork

    INTEGER, CONTIGUOUS, POINTER ::   &   !< Pointer to line and block indices (array)
      &  iidx(:,:,:) => NULL(),       &   !< of edges
      &  iblk(:,:,:) => NULL()

    REAL(wp):: z_dthalf                   !< 0.5 * timestep

    REAL(wp)::  &                         !< horizontal fluxes of wave energy
      &  z_mflx_tracer_h(SIZE(p_tracer_now,1),SIZE(p_tracer_now,2),p_patch%nblks_e)

    REAL(wp):: z_fluxdiv_c                !< flux divergence at cell center

    TYPE(t_back_traj) :: btraj      !< backward trajectories for MIURA

    ! the following dummy fields are not needed for 2D transport of wave spectral energy.
    ! It is, however required by the horizontal transport routine, which is taken over from
    ! the atmosphere model.
    !
    ! dummy density weighted cell height [kg/m**2]
    ! set to 1 below
    REAL(wp)::  &
      &  z_rhodz(SIZE(p_tracer_now,1),SIZE(p_tracer_now,2),SIZE(p_tracer_now,3))
    !
    ! dummy lateral boundary tendencies of transported tracer quantity
    ! in preparation for WAVE-LAM
    ! set to 0 below
    REAL(wp), TARGET:: &
      &  z_grf_tend_tracer(SIZE(p_tracer_now,1),SIZE(p_tracer_now,2),SIZE(p_tracer_now,3),SIZE(p_tracer_now,4))
    !
    REAL(wp), POINTER, CONTIGUOUS:: p_grf_tend_tracer(:,:,:,:)
    !

    !-----------------------------------------------------------------------

    IF(ltimer) CALL timer_start(timer_transport)

    ! halo synchronization for spectral energy, before transport
    CALL sync_patch_array_mult(typ        = SYNC_C,               &
      &                        p_patch    = p_patch,              &
      &                        nfields    = SIZE(p_tracer_now,4), &
      &                        f4din      = p_tracer_now,         &
      &                        opt_varname='p_tracer_now')


    ! pointer to energy_propagation_config to save some paperwork
    enprop_conf => energy_propagation_config

    ! line and block indices of edges from a cell perspective
    iidx => p_patch%cells%edge_idx
    iblk => p_patch%cells%edge_blk

    p_grf_tend_tracer => z_grf_tend_tracer(:,:,:,:)


    !$OMP PARALLEL
    CALL init(init_var=z_rhodz, init_val=1._wp, lacc=.FALSE.)
    CALL init(init_var=z_grf_tend_tracer, lacc=.FALSE.)
    CALL init(init_var=z_mflx_tracer_h, lacc=.FALSE.)
    !$OMP END PARALLEL


    ! start and end cells for flux divergence calculation
    i_rlstart_c  = 1
    i_rlend_c    = min_rlcell_int
    i_startblk_c = p_patch%cells%start_block(i_rlstart_c)
    i_endblk_c   = p_patch%cells%end_block(i_rlend_c)

    ! start and end egdes for flux computation and backward trajectories
    i_rlstart_e  = 1
    i_rlend_e    = min_rledge_int - 1

    ! start and end cells for lateral boundary updates
    i_rlstart_bdy  = 1
    i_rlend_bdy    = grf_bdywidth_c
    i_startblk_bdy = p_patch%cells%start_block(i_rlstart_bdy)
    i_endblk_bdy   = p_patch%cells%end_block(i_rlend_bdy)


    ! initialize backward trajectory calculation
    CALL btraj%construct(nproma,p_patch%nlev,p_patch%nblks_e,2)

    z_dthalf = 0.5_wp * p_dtime

    ! there is only one dummy vertical level
    jk = 1

    TRACER_ADV: DO jt = 1, ntracer ! Tracer loop

      ! 1st order backward trajectory
      ! note, that the group velocity may depend on the frequency
      ! Hence, the computation of backward trajectories is required
      ! for every energy bin.
      !
      CALL btraj_compute_o1( btraj       = btraj,              & !inout
        &                  ptr_p         = p_patch,            & !in
        &                  ptr_int       = p_int_state,        & !in
        &                  p_vn          = p_vn_traj(:,:,:,jt),& !in
        &                  p_vt          = p_vt_traj(:,:,:,jt),& !in
        &                  p_dthalf      = z_dthalf,           & !in
        &                  opt_rlstart   = i_rlstart_e,        & !in
        &                  opt_rlend     = i_rlend_e,          & !in
        &                  opt_slev      = 1,                  & !in
        &                  opt_elev      = 1,                  & !in
        &                  opt_acc_async = .TRUE.              ) !in


      ! compute horizontal fluxes of wave energy
      !
      ! CALL MIURA with second order accurate reconstruction
      CALL upwind_hflux_miura(                                &
        &         p_patch         = p_patch,                  & !in
        &         p_cc            = p_tracer_now(:,:,:,jt),   & !in
        &         p_mass_flx_e    = p_mflx_h(:,:,:,jt),       & !in
        &         p_dtime         = p_dtime,                  & !in
        &         p_int           = p_int_state,              & !in
        &         btraj           = btraj,                    & !in
        &         p_igrad_c_miura = enprop_conf%igrad_c_miura,& !in
        &         p_itype_hlimit  = enprop_conf%itype_limit,  & !in
        &         p_out_e         = z_mflx_tracer_h(:,:,:),   & !inout
        &         opt_rhodz_now   = z_rhodz(:,:,:),           & !in
        &         opt_rhodz_new   = z_rhodz(:,:,:),           & !in
        &         opt_lconsv      = llsq_lin_consv,           & !in
        &         opt_rlstart_e   = i_rlstart_e,              & !in
        &         opt_rlend_e     = i_rlend_e,                & !in
        &         opt_slev        = 1,                        & !in
        &         opt_elev        = 1                         ) !in



      ! update wave energy, by computing the horizontal flux divergence
      !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,z_fluxdiv_c)
      DO jb = i_startblk_c, i_endblk_c

        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
                     i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

        ! compute horizontal flux divergences and update wave energy
        !
        DO jc = i_startidx, i_endidx
          z_fluxdiv_c =                                                                         &
            & z_mflx_tracer_h(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_int_state%geofac_div(jc,1,jb) + &
            & z_mflx_tracer_h(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_int_state%geofac_div(jc,2,jb) + &
            & z_mflx_tracer_h(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_int_state%geofac_div(jc,3,jb)

           ! update the wave energy, by applying the flux divergence
           !
           p_tracer_new(jc,jk,jb,jt) = p_tracer_now(jc,jk,jb,jt) - p_dtime * z_fluxdiv_c
        ENDDO  !jc
      ENDDO  !jb
!$OMP END DO

      ! update lateral boundary values of wave energy with interpolated time tendencies
      ! (limited area mode only)
      !
      IF (l_limited_area .OR. p_patch%id > 1) THEN

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
        DO jb = i_startblk_bdy, i_endblk_bdy

          CALL get_indices_c(p_patch, jb, i_startblk_bdy, i_endblk_bdy, &
                             i_startidx, i_endidx, i_rlstart_bdy, i_rlend_bdy)

          ! Tracer values are clipped here to avoid generation of negative values
          ! For mass conservation, a correction has to be applied in the
          ! feedback routine anyway
          DO jc = i_startidx, i_endidx
            p_tracer_new(jc,jk,jb,jt) =                            &
              &     MAX(0._wp, p_tracer_now(jc,jk,jb,jt)           &
              &   + p_dtime * p_grf_tend_tracer(jc,jk,jb,jt) )
          ENDDO
        ENDDO  !jb
!$OMP END DO NOWAIT

      ENDIF ! l_limited_area
!$OMP END PARALLEL

    END DO TRACER_ADV


    CALL btraj%destruct()


    ! Calculate wave refraction
    IF (enprop_conf%lgrid_refr) THEN
      CALL wave_refraction(p_patch     = p_patch,                  & !in
        &                  wave_config = wave_config,              & !in
        &                  dtime       = p_dtime,                  & !in
        &                  wave_num_c  = wave_num_c(:,:,:),        & !in
        &                  gv_c        = gv_c(:,:,:),              & !in
        &                  depth       = bathymetry_c(:,:),        & !in
        &                  depth_grad  = geo_bath_grad_c(:,:,:,:), & !in
        &                  tracer_now  = p_tracer_now(:,:,:,:),    & !in
        &                  tracer_new  = p_tracer_new(:,:,:,:))      !inout

      ! Set energy to absolute allowed minimum
      CALL set_energy2emin(p_patch     = p_patch,              & !in
        &                  wave_config = wave_config,          & !in
        &                  tracer      = p_tracer_new(:,:,:,:))  !inout
    END IF


    IF (ltimer) CALL timer_stop(timer_transport)

  END SUBROUTINE wave_step_advection

END MODULE mo_wave_advection_stepping

