!
!
! Prepares and postprocesses the fields for and from nwp physics
!
! Depending on the action item different sets of physics will be called:
! this is
! 1. condensation only so to apply saturation adjustment where needed
! 2. 'slow physics' means up to now the whole physical package beside
!     microphysics.
! 3. turbulence, microphysics and condensation are considered fast physical package
! 4. Updating the moist tracers in synchrone time intervalls to
!    advection and saturation adjustment
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

! Workaround note: With Cray Fortran 8.5.5, a segmentation fault occurred in
! the SUBROUTINE radheat (source line 1892) when accessing the dummy array
! "pqv".  the workaround here is to copy the "pqv" dummy array to a temporary
! pqv=prm_diag%tot_cld(:,:,jb,iqv)

MODULE mo_nh_interface_nwp

  USE mtime,                      ONLY: datetime
  USE mo_util_mtime,              ONLY: getElapsedSimTimeInSeconds
  USE mo_kind,                    ONLY: wp

  USE mo_timer
  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_impl_constants,          ONLY: itconv, itccov, itrad, itgscp,                        &
    &                                   itsatad, itturb, itsfc, itradheat,                    &
    &                                   itsso, itgwd, itfastphy, icosmo, igme, ivdiff,        &
    &                                   min_rlcell_int, min_rledge_int, min_rlcell, ismag, iprog
  USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,             ONLY: get_indices_c, get_indices_e
  USE mo_intp_rbf,                ONLY: rbf_vec_interpol_cell
  USE mo_model_domain,            ONLY: t_patch
  USE mo_intp_data_strc,          ONLY: t_int_state
  USE mo_nonhydro_types,          ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydrostatic_config,   ONLY: kstart_moist, ih_clch, ih_clcm, lcalc_dpsdt
  USE mo_nwp_lnd_types,           ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_ext_data_types,          ONLY: t_external_data
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag, t_nwp_phy_tend, t_nwp_phy_stochconv
  USE mo_coupling_config,         ONLY: is_coupled_to_ocean, is_coupled_to_waves, is_coupled_to_hydrodisc
  USE mo_parallel_config,         ONLY: nproma, p_test_run, use_physics_barrier
  USE mo_diffusion_config,        ONLY: diffusion_config
  USE mo_initicon_config,         ONLY: is_iau_active
  USE mo_run_config,              ONLY: ntracer, iqv, iqc, iqi, iqs, iqr, iqg, iqtke,  &
    &                                   msg_level, ltimer, timers_level, lart, ldass_lhn
  USE mo_grid_config,             ONLY: l_limited_area
  USE mo_physical_constants,      ONLY: rd, rd_o_cpd, vtmpc1, p0ref, rcvd, cvd, cvv, tmelt, grav

  USE mo_nh_diagnose_pres_temp,   ONLY: diagnose_pres_temp, diag_pres, diag_temp, calc_qsum
  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config, iprog_aero
  USE mo_iau,                     ONLY: iau_update_tracer
  USE mo_util_phys,               ONLY: tracer_add_phytend, inversion_height_index
  USE mo_lnd_nwp_config,          ONLY: ntiles_total, ntiles_water
  USE mo_cover_koe,               ONLY: cover_koe, cover_koe_config
  USE mo_satad,                   ONLY: satad_v_3D, satad_v_3D_gpu, latent_heat_sublimation
  USE mo_aerosol_util,            ONLY: prog_aerosol_2D
  USE mo_radiation,               ONLY: radheat, pre_radiation_nwp
  USE mo_radiation_config,        ONLY: irad_aero, iRadAeroTegen, iRadAeroCAMSclim, iRadAeroART
  USE mo_nwp_gw_interface,        ONLY: nwp_gwdrag
  USE mo_nwp_gscp_interface,      ONLY: nwp_microphysics
  USE mo_nwp_turbtrans_interface, ONLY: nwp_turbtrans
  USE mo_nwp_turbdiff_interface,  ONLY: nwp_turbdiff
  USE mo_nwp_sfc_interface,       ONLY: nwp_surface
  USE mo_nwp_conv_interface,      ONLY: nwp_convection
  USE mo_nwp_rad_interface,       ONLY: nwp_radiation
  USE mo_nwp_vdiff_interface,     ONLY: nwp_vdiff, nwp_vdiff_update_seaice_list, &
    &                                   nwp_vdiff_update_seaice
  USE mo_turb_vdiff_config,       ONLY: vdiff_config
  USE mo_ccycle_config,           ONLY: ccycle_config
  USE mo_nwp_ocean_coupling,      ONLY: nwp_couple_ocean
  USE mo_nwp_hydrodisc_coupling,  ONLY: nwp_couple_hydrodisc
  USE mo_sync,                    ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E,      &
                                        SYNC_C, SYNC_C1
  USE mo_mpi,                     ONLY: my_process_is_mpi_all_parallel, work_mpi_barrier
  USE mo_nwp_diagnosis,           ONLY: nwp_statistics, nwp_opt_diagnostics_2, &
                                    &   nwp_diag_output_1, nwp_diag_output_2




  USE mo_var_list,                ONLY: t_var_list_ptr

  USE mo_ls_forcing_nml,          ONLY: is_ls_forcing, is_nudging_uv, is_nudging_tq, is_sim_rad, &
    &                                   nudge_start_height, nudge_full_height, dt_relax
  USE mo_ls_forcing,              ONLY: apply_ls_forcing
  USE mo_les_turb_interface,      ONLY: les_turbulence
  USE mo_les_config,              ONLY: les_config

  USE mo_sim_rad,                 ONLY: sim_rad
  USE mo_advection_config,        ONLY: advection_config
  USE mo_o3_util,                 ONLY: calc_o3_gems
  USE mo_nh_supervise,            ONLY: compute_dpsdt

  USE mo_radar_data_state,        ONLY: radar_data, lhn_fields
  USE mo_latent_heat_nudging,     ONLY: organize_lhn
  USE mo_assimilation_config,     ONLY: assimilation_config
  USE mo_nwp_reff_interface,      ONLY: set_reff , combine_phases_radiation_reff
  USE mo_upatmo_impl_const,       ONLY: iUpatmoPrcStat, iUpatmoStat
  USE mo_upatmo_config,           ONLY: upatmo_config

  USE mo_nwp_upatmo_interface,    ONLY: nwp_upatmo_interface, nwp_upatmo_update

  USE mo_fortran_tools,           ONLY: set_acc_host_or_device, copy, init





  USE mo_nwp_sfc_utils,           ONLY: process_sst_and_seaice
  ! SPPT
  USE mo_sppt_state,              ONLY: sppt
  USE mo_sppt_config,             ONLY: sppt_config
  USE mo_sppt_util,               ONLY: construct_rn
  USE mo_sppt_core,               ONLY: calc_tend, pert_tend, apply_tend, save_state


  USE mo_nwp_tuning_config,       ONLY: tune_sc_eis
  USE mo_sbm_storage,             ONLY: t_sbm_storage, get_sbm_storage

  !$ser verbatim USE mo_ser_all,              ONLY: serialize_all

  IMPLICIT NONE

  PRIVATE


  REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
  REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd

  PUBLIC :: nwp_nh_interface

CONTAINS
  !
  !-----------------------------------------------------------------------
  !
  SUBROUTINE nwp_nh_interface(lcall_phy_jg, linit, lredgrid,       & !input
                            & dt_loc, dt_phy_jg,                   & !input
                            & mtime_datetime,                      & !input
                            & pt_patch, pt_int_state, p_metrics,   & !input
                            & pt_par_patch,                        & !input
                            & ext_data,                            & !input
                            & pt_prog,                             & !inout
                            & pt_prog_now_rcf,                     & !inout
                            & pt_prog_rcf,                         & !inout
                            & pt_diag ,                            & !inout
                            & prm_diag, prm_nwp_tend,              & !inout
                            & prm_nwp_stochconv, lnd_diag,         & !inout
                            & lnd_prog_now, lnd_prog_new,          & !inout
                            & wtr_prog_now, wtr_prog_new,          & !inout
                            & p_prog_list,                         & !in
                            & lacc                                 ) !in

    !>
    ! !INPUT PARAMETERS:

    LOGICAL, INTENT(IN)          ::   &             !< physics package time control (switches)
         &                          lcall_phy_jg(:) !< for domain jg
    LOGICAL, INTENT(IN)          :: linit           !< .TRUE. if initialization call (this switch is currently used
                                                    !  to call turbtran in addition to the slow-physics routines)
    LOGICAL, INTENT(IN)          :: lredgrid        !< use reduced grid for radiation
    REAL(wp),INTENT(in)          :: dt_loc          !< (advective) time step applicable to local grid level
    REAL(wp),INTENT(in)          :: dt_phy_jg(:)    !< time interval for all physics
                                                    !< packages on domain jg
    TYPE(datetime), POINTER                :: mtime_datetime !< date/time information (in)
    TYPE(t_patch),     TARGET,INTENT(inout):: pt_patch     !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in):: pt_par_patch !<grid/patch info (parent grid)

    TYPE(t_int_state),    TARGET,INTENT(in):: pt_int_state      !< interpolation state
    TYPE(t_nh_metrics)   ,       INTENT(in):: p_metrics
    TYPE(t_external_data),       INTENT(inout):: ext_data

    TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag     !<the diagnostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog_now_rcf !<old state for tke
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog_rcf !<the prognostic variables (with
                                                          !< red. calling frequency for tracers!
    TYPE(t_nwp_phy_diag),       INTENT(inout) :: prm_diag
    TYPE(t_nwp_phy_tend),TARGET,INTENT(inout) :: prm_nwp_tend
    TYPE(t_nwp_phy_stochconv),  INTENT(inout) :: prm_nwp_stochconv
    TYPE(t_lnd_prog),           INTENT(inout) :: lnd_prog_now, lnd_prog_new
    TYPE(t_wtr_prog),           INTENT(inout) :: wtr_prog_now, wtr_prog_new
    TYPE(t_lnd_diag),           INTENT(inout) :: lnd_diag

    TYPE(t_var_list_ptr), INTENT(inout) :: p_prog_list !current prognostic state list

    LOGICAL, INTENT(in), OPTIONAL :: lacc !flag to run on GPU

    ! !OUTPUT PARAMETERS:            !<variables induced by the whole physics
    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices

    ! Local scalars:

    INTEGER :: jc,jk,jb,jce,isubs!loop indices
    INTEGER :: jg,jgc            !domain id
    INTEGER :: convind           !help variable to circument compiler issue

    LOGICAL :: ltemp, lpres, ltemp_ifc, l_any_fastphys, l_any_slowphys
    LOGICAL :: lcall_lhn, lcall_lhn_v, lapply_lhn, lcall_lhn_c  !< switches for latent heat nudging
    LOGICAL :: lcompute_tt_lheat                                !< TRUE: store temperature tendency
                                                                ! due to grid scale microphysics
                                                                ! and satad for latent heat nudging

    LOGICAL :: l_any_upatmophys

    INTEGER,  POINTER ::  iidx(:,:,:), iblk(:,:,:)

    REAL(wp), TARGET :: &                                       !> temporal arrays for
      & z_ddt_u_tot (nproma,pt_patch%nlev,pt_patch%nblks_c),&
      & z_ddt_v_tot (nproma,pt_patch%nlev,pt_patch%nblks_c),&   !< hor. wind tendencies
      & z_ddt_temp  (nproma,pt_patch%nlev)                      !< Temperature tendency

    REAL(wp) :: z_exner_sv(nproma,pt_patch%nlev,pt_patch%nblks_c), z_tempv, sqrt_ri(nproma), n2, dvdz2, &
      zddt_u_raylfric(nproma,pt_patch%nlev), zddt_v_raylfric(nproma,pt_patch%nlev), convfac, wfac

    !< SBM microphysics:
    TYPE(t_sbm_storage), POINTER:: ptr_sbm_storage =>NULL()   ! pointer to SBM storage object

    !< vertical interfaces

    REAL(wp) :: zsct ! solar constant (at time of year)
    REAL(wp) :: zcosmu0 (nproma,pt_patch%nblks_c), cosmu0_slope(nproma,pt_patch%nblks_c)

    REAL(wp) :: z_qsum(nproma,pt_patch%nlev)       !< summand of virtual increment
    REAL(wp) :: z_ddt_alpha(nproma,pt_patch%nlev)  !< tendency of virtual increment

    ! auxiliaries for Rayleigh friction computation
    REAL(wp) :: vabs, rfric_fac, ustart, uoffset_q, ustart_q, max_relax

    ! Variables for LHN
    REAL(wp) :: dhumi_lhn,dhumi_lhn_tot

    ! communication ids, these do not need to be different variables,
    ! since they are not treated individualy
    INTEGER :: ddt_u_tot_comm, ddt_v_tot_comm, z_ddt_u_tot_comm, z_ddt_v_tot_comm, &
      & tracers_comm, tempv_comm, exner_pr_comm, w_comm

    INTEGER :: ntracer_sync
    LOGICAL :: lzacc ! non-optional version of lacc


    ! inversion height diagnostic for cover_koe with stratocumulus
    INTEGER :: kc_inversion(nproma), kc_entr_zone(nproma)
    LOGICAL :: lfound_inversion(nproma)

    ! Pointer to IDs of tracers which contain prognostic condensate.
    ! Required for computing the water loading term
    INTEGER, POINTER :: condensate_list(:)

    REAL(wp) :: p_sim_time      !< elapsed simulation time on this grid level

    LOGICAL :: lcalc_inv
    
    ! SCM Nudging
    REAL(wp) :: nudgecoeff

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (ltimer) CALL timer_start(timer_physics)

    ! calculate elapsed simulation time in seconds (local time for
    ! this domain!)
    p_sim_time = getElapsedSimTimeInSeconds(mtime_datetime)

    ! local variables related to the blocking

    jg        = pt_patch%id

    IF (pt_patch%n_childdom > 0) THEN
      jgc = pt_patch%child_id(jg)
    ELSE
      jgc = jg
    ENDIF

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !define pointers
    iidx  => pt_patch%edges%cell_idx
    iblk  => pt_patch%edges%cell_blk


    IF (lcall_phy_jg(itsatad) .OR. lcall_phy_jg(itgscp) .OR. &
        lcall_phy_jg(itturb)  .OR. lcall_phy_jg(itsfc)) THEN
      l_any_fastphys = .TRUE.
    ELSE
      l_any_fastphys = .FALSE.
    ENDIF


    IF (lcall_phy_jg(itrad) .OR.  lcall_phy_jg(itconv) .OR. lcall_phy_jg(itccov)  &
       .OR. lcall_phy_jg(itsso) .OR. lcall_phy_jg(itgwd)) THEN
      l_any_slowphys = .TRUE.
    ELSE
      l_any_slowphys = .FALSE.
    ENDIF

    ! upper-atmosphere physics
    IF (upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%enabled )) THEN
      l_any_upatmophys = (.NOT. upatmo_config(jg)%nwp_phy%isBeforeOpPhase(mtime_datetime)) .AND. &
        &                (.NOT. upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%afterActivePhase ))
    ELSE
      l_any_upatmophys = .FALSE.
    ENDIF

    ! condensate tracer IDs
    condensate_list => advection_config(jg)%trHydroMass%list


    ! Check whether latent heat nudging is active
    !
    IF (ldass_lhn .AND. assimilation_config(jg)%lvalid_data .AND. .NOT. linit) THEN
      !
      IF ( jg == jgc )  CALL assimilation_config(jg)%dass_g%reinitEvents()
      lcall_lhn   = assimilation_config(jg)%dass_lhn%isActive(mtime_datetime)
      lcall_lhn_v = assimilation_config(jg)%dass_lhn_verif%isActive(mtime_datetime)
      IF (msg_level >= 15) CALL assimilation_config(jg)%dass_g%printStatus(mtime_datetime)
      !
      lcompute_tt_lheat = lcall_lhn .OR. lcall_lhn_v
    ELSE
      lcompute_tt_lheat = .FALSE.
    ENDIF

    ! Inversion height is calculated only if the threshold is set to a non-default value
    lcalc_inv = tune_sc_eis < 1000._wp
    
    IF(sppt_config(jg)%lsppt .AND. .NOT. linit) THEN
      ! Construct field of random numbers for SPPT
      CALL construct_rn (pt_patch, mtime_datetime, sppt_config(jg), sppt(jg)%rn_3d, &
        &                sppt(jg)%rn_2d_now, sppt(jg)%rn_2d_new, lacc=lzacc)
    ENDIF ! end of lsppt

    !$ACC DATA CREATE(zddt_v_raylfric, zddt_u_raylfric, sqrt_ri, z_ddt_temp, z_ddt_alpha, z_ddt_v_tot) &
    !$ACC   CREATE(zcosmu0, z_ddt_u_tot, z_exner_sv, z_qsum, kc_inversion, kc_entr_zone, lfound_inversion) IF(lzacc)
    !$ACC DATA COPYIN(dt_phy_jg)

    IF ( lcall_phy_jg(itturb) .OR. lcall_phy_jg(itconv) .OR.           &
         lcall_phy_jg(itsso)  .OR. lcall_phy_jg(itgwd) .OR. linit ) THEN

      !-------------------------------------------------------------------------
      !>
      !!   Interpolation from v_n onto u,v =>  Reconstruct u and v
      !!   This is needed for turbulence, convection and SSO/GWdrag
      !!
      !-------------------------------------------------------------------------

      IF (msg_level >= 15) &
           & CALL message('mo_nh_interface_nwp:', 'reconstruct u/v')

      IF (timers_level > 3) CALL timer_start(timer_phys_u_v)

      CALL rbf_vec_interpol_cell(pt_prog%vn,            & !< normal wind comp.
        &                        pt_patch,              & !< patch
        &                        pt_int_state,          & !< interpolation state
        &                        pt_diag%u, pt_diag%v,  & !<  reconstr. u,v wind
        &                        opt_rlend=min_rlcell_int )

      IF (timers_level > 3) CALL timer_stop(timer_phys_u_v)

    ENDIF ! diagnose u/v

    IF (msg_level >= 18) THEN

      ! Diagnose temperature needed for debugging output
      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,    &
           &                              pt_diag, pt_patch,       &
           &                              opt_calc_temp=.TRUE.,    &
           &                              opt_calc_pres=.FALSE.,   &
           &                              opt_rlend=min_rlcell_int)

      ! Write extensive debugging output
      CALL nwp_diag_output_1(pt_patch, pt_diag, pt_prog_rcf)

    ENDIF


    !-------------------------------------------------------------------------
    !>  Update the slow-physics tendencies on the tracer fields,
    !!  afterwards perform saturation adjustment
    !-------------------------------------------------------------------------

    IF (msg_level >= 15) CALL message('mo_nh_interface_nwp:', 'satad')
    IF (timers_level > 2) CALL timer_start(timer_satad_v_3D)

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    ! Diagnose pressure and temperature in nest boundary zone
    ! (needed for the reduced radiation grid)
    rl_start = 1
    rl_end   = grf_bdywidth_c

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)


    IF (jg > 1 .OR. l_limited_area) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk,  &
                             i_startidx, i_endidx, rl_start, rl_end)

        CALL diag_temp (pt_prog, pt_prog_rcf, condensate_list, pt_diag,    &
                        jb, i_startidx, i_endidx, 1, kstart_moist(jg), nlev)
        
        CALL diag_pres (pt_prog, pt_diag, p_metrics, jb, i_startidx, i_endidx, 1, nlev)

      ENDDO
!$OMP END DO NOWAIT
    ENDIF

    ! Save Exner pressure field on halo points (prognostic points are treated below)
    rl_start = min_rlcell_int-1
    rl_end   = min_rlcell

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk,  &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      z_exner_sv(i_startidx:i_endidx,:,jb) = pt_prog%exner(i_startidx:i_endidx,:,jb)
      !$ACC END KERNELS
    ENDDO
!$OMP END DO NOWAIT


    ! computations on prognostic points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)


!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_qsum,z_tempv) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)


      IF (.NOT. linit) THEN

        IF (is_iau_active) THEN

          ! add analysis increments from data assimilation to qv (during IAU phase)
          CALL iau_update_tracer( pt_prog     = pt_prog,     & !in
           &                      p_metrics   = p_metrics,   & !in
           &                      pt_diag     = pt_diag,     & !inout
           &                      pt_prog_rcf = pt_prog_rcf, & !inout tracer
           &                      jg          = jg,          & !in
           &                      jb          = jb,          & !in
           &                      i_startidx  = i_startidx,  & !in
           &                      i_endidx    = i_endidx,    & !in
           &                      kend        = nlev,        & !in
           &                      lacc        = lzacc        ) !in
        ENDIF


        ! The provisional "new" tracer state, resulting from the advection
        ! step, still needs to be updated with the SLOW-physics tracer tendencies
        ! computed at the end of the last physics call for the then final
        ! "new" state. The corresponding update for the dynamics variables has
        ! already happened in the dynamical core.
        !
        CALL tracer_add_phytend( p_rho_now    = pt_prog%rho(:,:,jb),  & !in
          &                      prm_nwp_tend = prm_nwp_tend,         & !in
          &                      pdtime       = dt_phy_jg(itfastphy), & !in
          &                      prm_diag     = prm_diag,             & !inout phyfields
          &                      pt_prog_rcf  = pt_prog_rcf,          & !inout tracer
          &                      p_metrics    = p_metrics,            & !in
          &                      dt_loc       = dt_loc,               & !in
          &                      jg           = jg,                   & !in
          &                      jb           = jb,                   & !in
          &                      i_startidx   = i_startidx,           & !in
          &                      i_endidx     = i_endidx,             & !in
          &                      kend         = nlev,                 & !in
          &                      lacc         = lzacc                 ) !in


      ENDIF  ! linit


      IF (l_any_fastphys .OR. linit) THEN  ! diagnose temperature
        CALL diag_temp (pt_prog, pt_prog_rcf, condensate_list, pt_diag,    &
                        jb, i_startidx, i_endidx, 1, kstart_moist(jg), nlev)
      ENDIF


      IF(sppt_config(jg)%lsppt .AND. .NOT. linit) THEN
        ! Save prognostic/diagnostic variables for SPPT - Temperature and Tracer
        CALL save_state(jb, i_startidx, i_endidx, nlev,    &
          &             pt_diag%temp, pt_prog_rcf%tracer,  &
          &             sppt(jg), lacc=lzacc)
      ENDIF ! end of lsppt


      ! Save Exner pressure field (this is needed for a correction to reduce sound-wave generation by latent heating)
      !$ACC KERNELS ASYNC(1) IF(lzacc)
      z_exner_sv(i_startidx:i_endidx,:,jb) = pt_prog%exner(i_startidx:i_endidx,:,jb)
      !$ACC END KERNELS

      !!-------------------------------------------------------------------------
      !> Initial saturation adjustment (a second one follows at the end of the microphysics)
      !!-------------------------------------------------------------------------

      ! SBM microphysics
      ! store snapshots of qv and temp just before saturation adjustment
      IF (atm_phy_nwp_config(jg)%inwp_gscp == 8) THEN
        ptr_sbm_storage => get_sbm_storage(patch_id = jg)
!$OMP PARALLEL
        CALL copy(pt_prog_rcf%tracer(:,:,:,iqv), ptr_sbm_storage%qv_before_satad   )
        CALL copy(pt_diag%temp(:,:,:),           ptr_sbm_storage%temp_before_satad )
!$OMP END PARALLEL
      ENDIF


      IF (lcall_phy_jg(itsatad)) THEN

        ! initialize tt_lheat to be in used LHN
        IF (lcompute_tt_lheat) THEN
          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          prm_diag%tt_lheat (:,:,jb) = - pt_diag%temp   (:,:,jb)
          !$ACC END KERNELS
        ENDIF

        CALL satad_v_3D(                                &
           & maxiter  = 10                             ,& !> IN
           & tol      = 1.e-3_wp                       ,& !> IN
           & te       = pt_diag%temp       (:,:,jb)    ,& !> INOUT
           & qve      = pt_prog_rcf%tracer (:,:,jb,iqv),& !> INOUT
           & qce      = pt_prog_rcf%tracer (:,:,jb,iqc),& !> INOUT
           & rhotot   = pt_prog%rho        (:,:,jb)    ,& !> IN
           & idim     = nproma                         ,& !> IN
           & kdim     = nlev                           ,& !> IN
           & ilo      = i_startidx                     ,& !> IN
           & iup      = i_endidx                       ,& !> IN
           & klo      = kstart_moist(jg)               ,& !> IN
           & kup      = nlev                            & !> IN
           )
 
        CALL calc_qsum (pt_prog_rcf%tracer, z_qsum, condensate_list, jb, i_startidx, i_endidx, 1, kstart_moist(jg), nlev)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(z_tempv) COLLAPSE(2)
        DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            z_tempv                 = pt_diag%tempv(jc,jk,jb)

            pt_diag%tempv(jc,jk,jb) = pt_diag%temp(jc,jk,jb)                         &
              &                   * ( 1._wp +  vtmpc1                                &
              &                   * pt_prog_rcf%tracer(jc,jk,jb,iqv) - z_qsum(jc,jk) )

            pt_prog%exner(jc,jk,jb) = pt_prog%exner(jc,jk,jb) *        &
              & (1._wp+rd_o_cpd*(pt_diag%tempv(jc,jk,jb)/z_tempv-1._wp))

          ENDDO
        ENDDO
        !$ACC END PARALLEL

        ! initialize tt_lheat to be in used LHN
        IF (lcompute_tt_lheat) THEN
          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          prm_diag%tt_lheat (:,:,jb) = prm_diag%tt_lheat (:,:,jb) + pt_diag%temp   (:,:,jb)
          !$ACC END KERNELS
        ENDIF
      ENDIF

      IF (lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb) .OR. lcall_phy_jg(itsfc)) THEN
        ! diagnose pressure for subsequent fast-physics parameterizations
        CALL diag_pres (pt_prog, pt_diag, p_metrics, jb, i_startidx, i_endidx, 1, nlev)
      ENDIF


    ENDDO ! nblks


!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > 2) CALL timer_stop(timer_satad_v_3D)


    !!-------------------------------------------------------------------------
    !>  turbulent transfer and diffusion and microphysics
    !!
    !!  Because we consider the following physical processes as fast ones
    !!  we allow here the update of prognostic variables inside the subroutines
    !!  This means that the conversion back to the ICON-prognostic variables
    !!  has to be done afterwards
    !!-------------------------------------------------------------------------


    !For turbulence schemes NOT including the call to the surface scheme.
    !nwp_surface must even be called in inwp_surface = 0 because the
    !the lower boundary conditions for the turbulence scheme
    !are not set otherwise

    IF ( l_any_fastphys .AND. ( ANY( (/icosmo,igme,ismag,iprog/)==atm_phy_nwp_config(jg)%inwp_turb ) ) ) THEN

      IF (timers_level > 2) CALL timer_start(timer_nwp_surface)

      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "surface", .TRUE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)

       !> as pressure is needed only for an approximate adiabatic extrapolation
       !! of the temperature at the lowest model level towards ground level,
       !! a recalculation is not required

       CALL nwp_surface    (  dt_phy_jg(itfastphy),              & !>input
                             & pt_patch,                         & !>input
                             & ext_data,                         & !>input
                             & pt_prog, pt_prog_rcf,             & !>in/inout rcf=reduced calling freq.
                             & pt_diag, p_metrics,               & !>inout
                             & prm_diag,                         & !>inout
                             & lnd_prog_now, lnd_prog_new,       & !>inout
                             & wtr_prog_now, wtr_prog_new,       & !>inout
                             & lnd_diag,                         & !>input
                             & lacc=lzacc                         ) !>in

       !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "surface", .FALSE., opt_lupdate_cpu=.TRUE., opt_dt=mtime_datetime)
      IF (timers_level > 2) CALL timer_stop(timer_nwp_surface)
    END IF

    !Call to turbulent parameterization schemes
    IF (  lcall_phy_jg(itturb) ) THEN

      IF (timers_level > 1) CALL timer_start(timer_nwp_turbulence)

      SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)

      !Turbulence schemes NOT including the call to the surface scheme
      CASE(icosmo,igme)

      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "turbdiff", .TRUE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)
      ! compute turbulent diffusion (atmospheric column)
      CALL nwp_turbdiff   (  dt_phy_jg(itfastphy),              & !>in
                            & pt_patch, p_metrics,              & !>in
                            & ext_data,                         & !>in
                            & pt_prog,                          & !>in
                            & pt_prog_now_rcf, pt_prog_rcf,     & !>in/inout
                            & pt_diag,                          & !>inout
                            & prm_diag, prm_nwp_tend,           & !>inout
                            & wtr_prog_now,                     & !>in
                            & lnd_prog_now,                     & !>in
                            & lnd_diag,                         & !>in
                            & lacc=lzacc                        ) !>in

      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "turbdiff", .FALSE., opt_lupdate_cpu=.TRUE., opt_dt=mtime_datetime)

      CASE (ivdiff)

        !  The vdiff interface calls the land-surface scheme itself.
        CALL nwp_vdiff( &
            & mtime_datetime, dt_phy_jg(itfastphy), pt_patch, ccycle_config(jg), &
            & vdiff_config(jg), pt_prog, pt_prog_rcf, pt_diag, p_metrics, prm_diag, ext_data, &
            & lnd_diag, lnd_prog_new, wtr_prog_now, wtr_prog_new, prm_diag%nwp_vdiff_state, &
            & prm_nwp_tend, initialize=linit, lacc=lzacc &
          )

        IF ( is_coupled_to_hydrodisc() .AND. (.NOT. linit) ) THEN


          IF (ltimer) CALL timer_stop(timer_coupling)
        END IF

        IF (is_coupled_to_ocean()) THEN
          ! Sea-ice cover might change if ocean passed back new values.

          CALL nwp_vdiff_update_seaice ( &
              & pt_patch, .FALSE., lnd_diag%fr_seaice(:,:), ext_data%atm%list_sea, &
              & ext_data%atm%list_seaice, wtr_prog_new, lacc=lzacc &
            )

        ELSE
          CALL nwp_vdiff_update_seaice_list ( &
              & pt_patch, lnd_diag%fr_seaice(:,:), ext_data%atm%list_sea, &
              & ext_data%atm%list_seaice, lacc=lzacc &
            )
        END IF

      CASE(ismag,iprog)

        !----------------------------------------------------------------------------------
        !>  Additional syncs required for 3D turbulence.
        !----------------------------------------------------------------------------------
        IF(diffusion_config(jg)%lhdiff_w) THEN
          CALL sync_patch_array(SYNC_C, pt_patch, pt_prog%w)
        ENDIF

        ntracer_sync = iqc
        CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+5, pt_diag%temp, pt_diag%tempv, &
                                   pt_prog%exner, pt_diag%u, pt_diag%v, f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))

        CALL les_turbulence (  dt_phy_jg(itfastphy),             & !>in
                             & p_sim_time,                       & !>in
                             & pt_patch, p_metrics,              & !>in
                             & pt_int_state,                     & !>in
                             & pt_prog,                          & !>in
                             & pt_prog_now_rcf,                  & !>inout
                             & pt_prog_rcf,                      & !>inout
                             & pt_diag ,                         & !>inout
                             & prm_diag,prm_nwp_tend,            & !>inout
                             & lnd_prog_now,                     & !>in
                             & lnd_prog_new,                     & !>inout ONLY for idealized LES
                             & lnd_diag,                         & !>in
                             & lacc=lzacc                        ) !>in

      CASE DEFAULT

        CALL finish('mo_nh_interface_nwp:','unknown choice of turbulence scheme')

      END SELECT

      IF (timers_level > 1) CALL timer_stop(timer_nwp_turbulence)

    END IF

    !-------------------------------------------------------------------------
    !  prognostic microphysic and precipitation scheme
    !-------------------------------------------------------------------------
    IF ( lcall_phy_jg(itgscp)) THEN

      IF (msg_level >= 15) &
        & CALL message('mo_nh_interface_nwp:', 'microphysics')

      !> temperature and tracers have been updated by turbulence;
      !! an update of the pressure field is not needed because pressure
      !! is not needed at high accuracy in the microphysics scheme
      !! note: after the microphysics the second call to SATAD is within
      !!       the nwp_microphysics routine (first one is above)

      IF (timers_level > 1) CALL timer_start(timer_nwp_microphysics)

      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "microphysics", .TRUE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)
      CALL nwp_microphysics ( dt_phy_jg(itfastphy),             & !>input
                            & lcall_phy_jg(itsatad),            & !>input
                            & pt_patch, p_metrics,              & !>input
                            & pt_prog,                          & !>inout
                            & pt_prog_rcf%tracer,               & !>inout
                            & pt_prog_now_rcf%tke,              & !>in
                            & pt_diag ,                         & !>inout
                            & prm_diag, prm_nwp_tend,           & !>inout
                            & ext_data,                         & !>in
                            & lcompute_tt_lheat,                &
                            & lacc=lzacc ) !>in

      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "microphysics", .FALSE., opt_lupdate_cpu=.TRUE., opt_dt=mtime_datetime)

      IF (timers_level > 1) CALL timer_stop(timer_nwp_microphysics)

    ENDIF



    !!------------------------------------------------------------------
    !> Latent heat nudging (optional)
    !!------------------------------------------------------------------
    IF (ldass_lhn .AND. assimilation_config(jg)%lvalid_data .AND. .NOT. linit) THEN

      IF (msg_level >= 15) CALL message('mo_nh_interface_nwp:', 'applying LHN')
      IF (timers_level > 1) CALL timer_start(timer_datass)

      IF (lcall_lhn .OR. lcall_lhn_v) THEN
        !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "lhn", .TRUE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)
         CALL organize_lhn (   &
                               & dt_loc,                           & !>input
                               & p_sim_time,                       & ! in
                               & pt_patch, p_metrics,              & !>input
                               & pt_int_state,                     & !>input
                               & pt_prog_rcf,                      & !>inout
                               & pt_diag ,                         & !>inout
                               & prm_diag,                         & !>inout
                               & lhn_fields(jg),                   & !>inout
                               & radar_data(jg),                   &
                               & prm_nwp_tend,                     &
                               & mtime_datetime,                   &
                               & lcall_lhn, lcall_lhn_v,           &
                               & assimilation_config(jg)%lvalid_data)
        !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "lhn", .FALSE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)

        IF (msg_level >= 7 .AND. .NOT. assimilation_config(jg)%lvalid_data) THEN
          CALL message('mo_nh_interface_nwp:','LHN turned off due to lack of valid data')
        ENDIF
      ELSE IF (msg_level >= 15) THEN
        WRITE (message_text,'(a,f10.2,i5)') 'LHN not running because of specified times!',p_sim_time,jg
        CALL message('mo_nh_interface_nwp:', message_text)
      ENDIF


      lcall_lhn_c = assimilation_config(jgc)%dass_lhn%isActive(mtime_datetime)
      lapply_lhn  = (lcall_lhn .OR. lcall_lhn_c) .AND. assimilation_config(jg)%lvalid_data

      IF (lapply_lhn) THEN

        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int

        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,dhumi_lhn,dhumi_lhn_tot) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end )

          IF ( assimilation_config(jg)%lhn_updt_rule == 1 ) THEN  ! Humidity update of LHN goes to ice if qi>0 T<0

            ! Update in the two-moment scheme
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
                ! Temperature update
                pt_diag%temp(jc,jk,jb) = pt_diag%temp(jc,jk,jb) + lhn_fields(jg)%ttend_lhn(jc,jk,jb) * dt_loc

                ! The humididy update is made differntly because it afffects the ice nucleation

                ! Update directly ice if it is there some already there and it needs to grow, so no more ice is nucleated
                IF ( pt_prog_rcf%tracer(jc,jk,jb,iqi) > 1E-7 .AND. pt_diag%temp(jc,jk,jb) < 273.16 ) THEN
                  ! Calculate LH sublimation, add constant and cp/cv option
                  dhumi_lhn_tot = lhn_fields(jg)%qvtend_lhn(jc,jk,jb) * dt_loc
                  dhumi_lhn = MAX(dhumi_lhn_tot,-pt_prog_rcf%tracer(jc,jk,jb,iqi))
                  pt_prog_rcf%tracer(jc,jk,jb,iqi) = pt_prog_rcf%tracer(jc,jk,jb,iqi) + dhumi_lhn
                  pt_diag%temp(jc,jk,jb) = pt_diag%temp(jc,jk,jb) + &
                         & rcvd*latent_heat_sublimation(pt_diag%temp(jc,jk,jb))*dhumi_lhn
                  pt_prog_rcf%tracer(jc,jk,jb,iqv) = pt_prog_rcf%tracer(jc,jk,jb,iqv) + dhumi_lhn_tot - dhumi_lhn
                ELSE
                  pt_prog_rcf%tracer(jc,jk,jb,iqv) = pt_prog_rcf%tracer(jc,jk,jb,iqv) + lhn_fields(jg)%qvtend_lhn(jc,jk,jb) * dt_loc
                END IF
              ENDDO
            ENDDO

          ELSE ! Standard Update one-moment
            ! update prognostic variables
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
                pt_diag%temp(jc,jk,jb) = pt_diag%temp(jc,jk,jb) + lhn_fields(jg)%ttend_lhn(jc,jk,jb) * dt_loc
                pt_prog_rcf%tracer(jc,jk,jb,iqv) = pt_prog_rcf%tracer(jc,jk,jb,iqv) + lhn_fields(jg)%qvtend_lhn(jc,jk,jb) * dt_loc
              ENDDO
            ENDDO
            !$ACC END PARALLEL
          END IF
          !-------------------------------------------------------------------------
          !   call the saturation adjustment
          !-------------------------------------------------------------------------

          CALL satad_v_3D(                                &
               & maxiter  = 10                             ,& !> IN
               & tol      = 1.e-3_wp                       ,& !> IN
               & te       = pt_diag%temp       (:,:,jb)    ,& !> INOUT
               & qve      = pt_prog_rcf%tracer (:,:,jb,iqv),& !> INOUT
               & qce      = pt_prog_rcf%tracer (:,:,jb,iqc),& !> INOUT
               & rhotot   = pt_prog%rho        (:,:,jb)    ,& !> IN
               & idim     = nproma                         ,& !> IN
               & kdim     = nlev                           ,& !> IN
               & ilo      = i_startidx                     ,& !> IN
               & iup      = i_endidx                       ,& !> IN
               & klo      = kstart_moist(jg)               ,& !> IN
               & kup      = nlev                            ) !> IN

        ENDDO ! nblks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ENDIF

      IF (timers_level > 1) CALL timer_stop(timer_datass)

    ENDIF ! ldass_lhn




    IF (timers_level > 1) CALL timer_start(timer_fast_phys)

    ! Remark: in the (unusual) case that satad is used without any other physics,
    ! recalculation of the thermodynamic variables is duplicated here. However,
    ! this is the easiest way to combine minimization of halo communications
    ! with a failsafe flow control

    IF (msg_level >= 15) &
      & CALL message('mo_nh_interface_nwp:', 'recalculate thermodynamic variables')


    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx, z_qsum) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end )

      IF (lcall_phy_jg(itsatad) .OR. lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb)) THEN
        !-------------------------------------------------------------------------
        !>
        !! re-calculate scalar prognostic variables out of physics variables!
        !!
        !-------------------------------------------------------------------------

        CALL calc_qsum (pt_prog_rcf%tracer, z_qsum, condensate_list, jb, i_startidx, i_endidx, 1, kstart_moist(jg), nlev)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc =  i_startidx, i_endidx

            pt_diag%tempv(jc,jk,jb) =  pt_diag%temp(jc,jk,jb)          &
&                                  * ( 1._wp +  vtmpc1                 &
&                                  *  pt_prog_rcf%tracer(jc,jk,jb,iqv) &
&                                   - z_qsum(jc,jk) )

            pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
              &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

            pt_diag%exner_pr(jc,jk,jb) = pt_diag%exner_pr(jc,jk,jb) + &
              pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb)

            pt_prog%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) &
&                                       / pt_prog%exner(jc,jk,jb)

          ENDDO
        ENDDO
        !$ACC END PARALLEL

        ! compute dynamical temperature tendency from increments of Exner function and density
        ! the virtual increment is neglected here because this tendency is used only as
        ! input for the convection scheme, which is rather insensitive against this quantity
        IF ( lcall_phy_jg(itconv) ) THEN
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
            DO jc =  i_startidx, i_endidx

              pt_diag%ddt_temp_dyn(jc,jk,jb) = pt_diag%temp(jc,jk,jb)/dt_phy_jg(itfastphy) * &
                ( cpd_o_rd/pt_prog%exner(jc,jk,jb)*pt_diag%exner_dyn_incr(jc,jk,jb) -        &
                ( pt_diag%airmass_new(jc,jk,jb)-pt_diag%airmass_now(jc,jk,jb) ) /            &
                pt_diag%airmass_new(jc,jk,jb) )

            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

      ENDIF ! recalculation

      IF (lcall_phy_jg(itturb) .OR. linit .OR. l_any_slowphys) THEN
        ! rediagnose pressure
        CALL diag_pres (pt_prog, pt_diag, p_metrics, jb, i_startidx, i_endidx, 1, nlev)
      ENDIF

      IF (iprog_aero >= 1 .AND. .NOT. linit) THEN
        CALL prog_aerosol_2D (i_startidx, i_endidx, jg, nproma, nlev, dt_loc, iprog_aero,                &
          &                   prm_diag%aerosol(:,:,jb),prm_diag%aercl_ss(:,jb),prm_diag%aercl_or(:,jb),  &
          &                   prm_diag%aercl_bc(:,jb),prm_diag%aercl_su(:,jb),prm_diag%aercl_du(:,jb),   &
          &                   pt_prog%exner(:,:,jb),pt_diag%temp(:,:,jb),pt_prog_rcf%tracer(:,:,jb,iqv), &
          &                   prm_diag%cosmu0(:,jb),                                                     &
          &                   prm_diag%rain_gsp_rate(:,jb),prm_diag%snow_gsp_rate(:,jb),                 &
          &                   prm_diag%rain_con_rate(:,jb),prm_diag%snow_con_rate(:,jb),                 &
          &                   ext_data%atm%soiltyp(:,jb), ext_data%atm%plcov_t(:,jb,:),                  &
          &                   ext_data%atm%frac_t(:,jb,:),                                               &
          &                   lnd_prog_now%w_so_t(:,1,jb,:), lnd_prog_now%w_so_ice_t(:,1,jb,:),          &
          &                   lnd_diag%h_snow_t(:,jb,:), lnd_diag%t_seasfc(:,jb),                        &
          &                   ext_data%atm%lc_class_t(:,jb,:),                                           &
          &                   pt_prog%rho(:,nlev,jb), prm_diag%tcm_t(:,jb,:),                            &
          &                   pt_diag%u(:,nlev,jb), pt_diag%v(:,nlev,jb), prm_diag%sp_10m(:,jb),         &
          &                   ext_data%atm%emi_bc(:,jb), ext_data%atm%emi_oc(:,jb),                      &
          &                   ext_data%atm%emi_so2(:,jb), ext_data%atm%bcfire(:,jb),                     &
          &                   ext_data%atm%ocfire(:,jb), ext_data%atm%so2fire(:,jb),                     &
          &                   ext_data%atm%idx_lst_t(:,jb,:),                                            &
          &                   ext_data%atm%gp_count_t(jb,:), ext_data%atm%list_seawtr%ncount(jb),        &
          &                   ext_data%atm%list_seawtr%idx(:,jb))
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > 1) CALL timer_stop(timer_fast_phys)
    IF ( (lcall_phy_jg(itturb) .OR. linit) .AND. ( ANY((/icosmo,igme/)==atm_phy_nwp_config(jg)%inwp_turb) .OR. &
         (ANY((/ismag,iprog/)==atm_phy_nwp_config(jg)%inwp_turb) .AND. (les_config(jg)%isrfc_type==1)) ) ) THEN

      IF (timers_level > 1) CALL timer_start(timer_nwp_turbulence)

      ! compute turbulent transfer coefficients (atmosphere-surface interface)

      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "turbtrans", .TRUE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)
      CALL nwp_turbtrans  ( dt_phy_jg(itfastphy),             & !>in
                          & pt_patch, p_metrics,              & !>in
                          & ext_data,                         & !>in
                          & pt_prog,                          & !>in
                          & pt_prog_rcf,                      & !>inout
                          & pt_diag,                          & !>inout
                          & prm_diag,                         & !>inout
                          & prm_nwp_tend,                     & !>in
                          & wtr_prog_new,                     & !>in
                          & lnd_prog_new,                     & !>inout
                          & lnd_diag,                         & !>inout
                          & lacc=lzacc                         ) !>in
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "turbtrans", .FALSE., opt_lupdate_cpu=.TRUE., opt_dt=mtime_datetime)

      IF (timers_level > 1) CALL timer_stop(timer_nwp_turbulence)
 
    ENDIF !lcall(itturb)

    ! Calculate tendencies - temperatures and tracers
    IF (sppt_config(jg)%lsppt .AND. .NOT. linit) THEN

      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_block(rl_start)
      i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk,  &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Temperature
        CALL calc_tend(i_startidx, i_endidx, nlev, sppt(jg)%ddt_temp_fast(:,:,jb), pt_diag%temp(:,:,jb), sppt(jg)%temp_now(:,:,jb), dt_loc, lacc=lzacc)

        ! Wind components - use existing tendencies from turbulence
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            sppt(jg)%ddt_u_fast(jc,jk,jb) = prm_nwp_tend%ddt_u_turb(jc,jk,jb)
            sppt(jg)%ddt_v_fast(jc,jk,jb) = prm_nwp_tend%ddt_v_turb(jc,jk,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        ! Tracer
        CALL calc_tend(i_startidx, i_endidx, nlev, sppt(jg)%ddt_qv_fast(:,:,jb), pt_prog_rcf%tracer(:,:,jb,iqv), sppt(jg)%qv_now(:,:,jb), dt_loc, lacc=lzacc)
        CALL calc_tend(i_startidx, i_endidx, nlev, sppt(jg)%ddt_qi_fast(:,:,jb), pt_prog_rcf%tracer(:,:,jb,iqi), sppt(jg)%qi_now(:,:,jb), dt_loc, lacc=lzacc)
        CALL calc_tend(i_startidx, i_endidx, nlev, sppt(jg)%ddt_qr_fast(:,:,jb), pt_prog_rcf%tracer(:,:,jb,iqr), sppt(jg)%qr_now(:,:,jb), dt_loc, lacc=lzacc)
        CALL calc_tend(i_startidx, i_endidx, nlev, sppt(jg)%ddt_qs_fast(:,:,jb), pt_prog_rcf%tracer(:,:,jb,iqs), sppt(jg)%qs_now(:,:,jb), dt_loc, lacc=lzacc)
        CALL calc_tend(i_startidx, i_endidx, nlev, sppt(jg)%ddt_qc_fast(:,:,jb), pt_prog_rcf%tracer(:,:,jb,iqc), sppt(jg)%qc_now(:,:,jb), dt_loc, lacc=lzacc)

        IF ( iqg /= 0 ) THEN
          CALL calc_tend(i_startidx, i_endidx, nlev, sppt(jg)%ddt_qg_fast(:,:,jb), pt_prog_rcf%tracer(:,:,jb,iqg), sppt(jg)%qg_now(:,:,jb), dt_loc, lacc=lzacc)
        ENDIF

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      
    ENDIF ! end of lsppt

    !!-------------------------------------------------------------------------
    !!  slow physics part
    !!-------------------------------------------------------------------------

    IF (l_any_slowphys) THEN

      IF (msg_level >= 15) &
         CALL message('mo_nh_interface', 'diagnose pres/temp for slow physics')

      ! If slow physics is called without fast physics (which should happen
      ! at the initial time step only), temperature needs to be calculated
      ! Otherwise, temperature is up to date
      IF ( .NOT. (lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb))) THEN
        ltemp = .TRUE.
      ELSE
        ltemp = .FALSE.
      ENDIF


      ! Pressure has already been updated at the end of the fast physics part
      lpres = .FALSE.

      ! Temperature at interface levels is needed if irad_aero = 6/7/9
      IF ( lcall_phy_jg(itrad) .AND.  &
        & ( irad_aero == iRadAeroTegen .OR. irad_aero == iRadAeroCAMSclim .OR. irad_aero == iRadAeroART ) ) THEN
        ltemp_ifc = .TRUE.
      ELSE
        ltemp_ifc = .FALSE.
      ENDIF

      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,   &
        &                      pt_diag, pt_patch,                 &
        &                      opt_calc_temp     = ltemp,         &
        &                      opt_calc_pres     = lpres,         &
        &                      lnd_prog          = lnd_prog_new,  &
        &                      opt_calc_temp_ifc = ltemp_ifc,     &
        &                      opt_rlend         = min_rlcell_int)

    ENDIF


    !-------------------------------------------------------------------------
    !> Convection
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itconv)  ) THEN
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "convection", .TRUE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)

      IF (msg_level >= 15)  CALL message('mo_nh_interface', 'convection')

      IF (timers_level > 2) CALL timer_start(timer_nwp_convection)
      CALL nwp_convection (  dt_phy_jg(itconv),                 & !>input
                            & linit,                            & !>input
                            & pt_patch, p_metrics,              & !>input
                            & ext_data,                         & !>input
                            & pt_prog,                          & !>input
                            & pt_prog_rcf,                      & !>input
                            & mtime_datetime,                   & !>input
                            & pt_diag,                          & !>inout
                            & prm_diag,                         & !>inout
                            & prm_nwp_tend,                     & !>inout
                            & prm_nwp_stochconv,                & !>inout
                            & pt_int_state,                     & !>in
                            & lacc=lzacc                         ) !>in

      IF (timers_level > 2) CALL timer_stop(timer_nwp_convection)
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "convection", .FALSE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)

    ENDIF! convection


    !-------------------------------------------------------------------------
    !> Cloud cover
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itccov) ) THEN

      ! When using a reduced grid for radiation, part of the boundary points need
      ! to be included in order to compute spatial gradients for back-interpolation
      IF (lredgrid) THEN
        rl_start = grf_bdywidth_c-1
      ELSE
        rl_start = grf_bdywidth_c+1
      ENDIF
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_block(rl_start)
      i_endblk   = pt_patch%cells%end_block(rl_end)

      IF (msg_level >= 15) &
        &  CALL message('mo_nh_interface', 'cloud cover')

      IF (timers_level > 2) CALL timer_start(timer_cover_koe)

      !-------------------------------------------------------------------------
      !> Cloud water distribution: cloud cover, cloud water, cloud ice
      !  inwp_cldcover =
      !  (0) no clouds
      !  (1) diagnostic cloud cover
      !  (2) prognostic total water variance (not yet started)
      !  (3) clouds as in COSMO
      !  (4) clouds as in turbulence
      !  (5) grid-scale cloud cover [1 or 0]
      !-------------------------------------------------------------------------

      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "cover", .TRUE., opt_lupdate_cpu=.TRUE., opt_dt=mtime_datetime)
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        IF (lcalc_inv) THEN
          ! inversion height diagnostic for EIS-based stratocumulus parameterization in cover_koe
          ! ( for efficiency reasons this could be integrated in cover_koe and called with an index list )
          CALL inversion_height_index(                             &
             &  p_metrics%z_mc(:,:,jb),                          &
             &  p_metrics%z_ifc(:,nlev+1,jb),                    &
             &  pt_prog_rcf%tracer(:,:,jb,iqc),                  &
             &  pt_diag%temp(:,:,jb),                            &
             &  pt_diag%pres(:,:,jb),                            &
             &  i_startidx,i_endidx,kstart_moist(jg),nlev,nlev,  &
             &  kc_inversion(:),                                 &
             &  kc_entr_zone(:),                                 &
             &  lfound_inversion(:))
        ELSE
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            kc_inversion(jc)=0._wp
            kc_entr_zone(jc)=0._wp
            lfound_inversion(jc)=.FALSE.
          ENDDO
          !$ACC END PARALLEL
        ENDIF
        
        CALL cover_koe &
&             (kidia  = i_startidx ,   kfdia  = i_endidx  ,       & !! in:  horizonal begin, end indices
&              klon = nproma,  kstart = kstart_moist(jg)  ,       & !! in:  horiz. and vert. vector length
&              klev   = nlev                              ,       &
&              linit  = linit                             ,       &
&              dtime  = dt_phy_jg(itccov)                 ,       &
&              cover_koe_config = cover_koe_config(jg)    ,       & !! in:  physics config state
&              tt     = pt_diag%temp         (:,:,jb)     ,       & !! in:  temperature at full levels
&              pp     = pt_diag%pres         (:,:,jb)     ,       & !! in:  pressure at full levels
&              ps     = pt_diag%pres_sfc     (:,jb)       ,       & !! in:  surface pressure at full levels
&              t_g    = lnd_prog_new%t_g     (:,jb)       ,       & !! in:  surface temperature
&              pgeo   = p_metrics%geopot_agl (:,:,jb)     ,       & !! in:  geopotential height
&              deltaz = p_metrics%ddqz_z_full(:,:,jb)     ,       & !! in:  layer thickness
&              rho    = pt_prog%rho          (:,:,jb  )   ,       & !! in:  density
&              rcld   = prm_diag%rcld        (:,:,jb)     ,       & !! in:  standard deviation of saturation deficit
&              ldland = ext_data%atm%llsm_atm_c (:,jb)    ,       & !! in:  land/sea mask
&              ldcum  = prm_diag%locum       (:,jb)       ,       & !! in:  convection on/off
&              kcbot  = prm_diag%mbas_con    (:,jb)       ,       & !! in:  convective cloud base
&              kctop  = prm_diag%mtop_con    (:,jb)       ,       & !! in:  convective cloud top
&              ktype  = prm_diag%ktype       (:,jb)       ,       & !! in:  convection type
&              kcinv  = kc_inversion         (:)          ,       & !! in:  inversion height index
&              linversion = lfound_inversion (:)          ,       & !! in:  inversion height logical
&              peis     = prm_diag%conv_eis  (:,jb)       ,       & !! in:  estimated inversion strength
&              fac_ccqc = prm_diag%fac_ccqc  (:,jb)       ,       & !! in:  factor for CLC-QC relationship (for EPS perturbations) 
&              pmfude_rate = prm_diag%con_udd(:,:,jb,3)   ,       & !! in:  convective updraft detrainment rate
&              plu         = prm_diag%con_udd(:,:,jb,7)   ,       & !! in:  updraft condensate
&              pcore       = prm_diag%con_udd(:,:,jb,8)   ,       & !! in:  updraft core fraction
&              rhoc_tend= prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqc),& !! in:  convective rho_c tendency
&              qv     = pt_prog_rcf%tracer   (:,:,jb,iqv) ,       & !! in:  spec. humidity
&              qc     = pt_prog_rcf%tracer   (:,:,jb,iqc) ,       & !! in:  cloud water
&              qi     = pt_prog_rcf%tracer   (:,:,jb,iqi) ,       & !! in:  cloud ice
&              qs     = pt_prog_rcf%tracer   (:,:,jb,iqs) ,       & !! in:  snow
&              lacc=lzacc                                  ,       & !! in
&              ttend_clcov = prm_nwp_tend%ddt_temp_clcov(:,:,jb) ,& !! out: temp tendency from sgs condensation
&              cc_tot = prm_diag%clc         (:,:,jb)     ,       & !! out: cloud cover
&              qv_tot = prm_diag%tot_cld     (:,:,jb,iqv) ,       & !! out: qv       -"-
&              qc_tot = prm_diag%tot_cld     (:,:,jb,iqc) ,       & !! out: clw      -"-
&              qc_sgs = prm_diag%qc_sgs      (:,:,jb) ,           & !! inout: sgs clw from RH scheme
&              qi_tot = prm_diag%tot_cld     (:,:,jb,iqi)         ) !! out: ci       -"-

      ENDDO
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "cover", .FALSE., opt_lupdate_cpu=.TRUE., opt_dt=mtime_datetime)

      IF (timers_level > 2) CALL timer_stop(timer_cover_koe)

    ENDIF! cloud cover


    !-------------------------------------------------------------------------
    !> Effective Radius
    !-------------------------------------------------------------------------

    !! Call effective radius diagnostic calculation for every cloud cover time step to achieve consistency
    !! between qc_dia, qi_dia and reff's. This also updates clc_rad, the special clc_diagnostic
    !! for radiation and satellite operators (MEC and synsats). 
    !! clc_rad is a copy of clc, but modified in the presence of large hydrometeors in a way
    !! that it is set to 1.0 if qr, qs or qg are present. The latter is needed by the
    !! radiation schemes and RTTOV in order to correctly take into account the radiative
    !! effects of these grid-scale hydrometeors.

    IF ( lcall_phy_jg(itccov) .AND. atm_phy_nwp_config(jg)%icalc_reff > 0 ) THEN
      
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "set_reff", .TRUE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)
      IF (msg_level >= 15) CALL message('mo_nh_interface', 'effective radius')

      IF (timers_level > 10) CALL timer_start(timer_phys_reff)

      CALL set_reff( prm_diag, pt_patch, pt_prog, pt_diag, ext_data )

      IF (  atm_phy_nwp_config(jg)%icpl_rad_reff == 1 .AND. atm_phy_nwp_config(jg)%icalc_reff /= 101 ) THEN

        ! .. Copy clc to clc_rad, combine reff's and enforce clc_rad = 1.0 at points where qs, qg, qr are present:
        CALL combine_phases_radiation_reff( prm_diag, pt_patch, pt_prog )

      END IF

      IF (timers_level > 10) CALL timer_stop(timer_phys_reff)
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "set_reff", .FALSE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)
    END IF



    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itrad) ) THEN

      IF (ltimer) CALL timer_start(timer_nwp_radiation)
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "radiation", .TRUE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)
      CALL nwp_radiation (lredgrid,              & ! in
           &              p_sim_time,            & ! in
           &              mtime_datetime,        & ! in
           &              pt_patch,              & ! in
           &              pt_par_patch,          & ! in
           &              ext_data,              & ! in
           &              lnd_diag,              & ! in
           &              pt_prog,               & ! inout
           &              pt_diag,               & ! inout
           &              prm_diag,              & ! inout
           &              lnd_prog_new,          & ! in
           &              wtr_prog_new,          & ! in
           &              p_metrics%z_mc,        & ! in
           &              p_metrics%z_ifc,       & ! in
           &              p_metrics%ddqz_z_full, & ! in
           &              lacc=lzacc              ) ! in, optional
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "radiation", .FALSE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)
      IF (ltimer) CALL timer_stop(timer_nwp_radiation)

    ENDIF


    IF ( lcall_phy_jg(itradheat) ) THEN
      !$ACC DATA CREATE(cosmu0_slope) IF(lzacc)
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "radheat", .TRUE., opt_lupdate_cpu=.FALSE., opt_dt=mtime_datetime)

      IF (msg_level >= 15) &
&           CALL message('mo_nh_interface', 'radiative heating')


      IF (timers_level > 10) CALL timer_start(timer_pre_radiation_nwp)

      CALL pre_radiation_nwp (                      &
        & kbdim      = nproma,                      &
        & p_inc_rad  = dt_phy_jg(itfastphy),        &
        & p_sim_time = p_sim_time,                  &
        & pt_patch   = pt_patch,                    &
        & zsmu0      = zcosmu0,                     &
        & zsct       = zsct,                        &
        & slope_ang  = p_metrics%slope_angle,       &
        & slope_azi  = p_metrics%slope_azimuth,     &
        & horizon    = ext_data%atm%horizon,        &
        & cosmu0_slp = cosmu0_slope,                &
        & lacc=lzacc                                 )

      IF (timers_level > 10) CALL timer_stop(timer_pre_radiation_nwp)

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_block(rl_start)
      i_endblk   = pt_patch%cells%end_block(rl_end)

      IF (timers_level > 2) CALL timer_start(timer_radheat)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,isubs,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
!
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          zcosmu0 (jc,jb) &
            = 0.5_wp * (ABS(zcosmu0(jc,jb)) + zcosmu0(jc,jb))

          !calculate solar incoming flux at TOA
          prm_diag%flxdwswtoa(jc,jb) = zcosmu0(jc,jb) * zsct   ! zsct by pre_radiation
        ENDDO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = 1, nproma
          prm_diag%swflxsfc (jc,jb)=0._wp
          prm_diag%lwflxsfc (jc,jb)=0._wp
          prm_diag%swflxtoa (jc,jb)=0._wp
          prm_diag%lwflxtoa (jc,jb)=0._wp
        ENDDO
        !$ACC END PARALLEL

        IF (atm_phy_nwp_config(jg)%inwp_surface >= 1 .OR. is_coupled_to_ocean()) THEN

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO isubs = 1, ntiles_total+ntiles_water
            DO jc = 1, nproma
              prm_diag%swflxsfc_t (jc,jb,isubs)=0._wp
              prm_diag%lwflxsfc_t (jc,jb,isubs)=0._wp
            ENDDO
          ENDDO
          !$ACC END PARALLEL

          CALL radheat (                   &
          !
          ! input
          ! -----
          !
          & jcs=i_startidx                         ,&! in     start index of inner do loop
          & jce=i_endidx                           ,&! in     end index of inner do loop
          & jg=pt_patch%id                         ,&! in     patch ID
          & kbdim=nproma                           ,&! in     loop length and dimension size
          & klev=nlev                              ,&! in     vertical dimension size
          & klevp1=nlevp1                          ,&! in     vertical dimension size
          & ntiles=ntiles_total                    ,&! in     number of tiles of sfc flux fields
          & ntiles_wtr=ntiles_water                ,&! in     number of extra tiles for ocean and lakes
          & pmair=pt_diag%airmass_new(:,:,jb)      ,&! in     layer air mass             [kg/m2]
          & pqv=prm_diag%tot_cld(:,:,jb,iqv)       ,&! in     specific moisture           [kg/kg]
          & pcd=cvd                                ,&! in     specific heat of dry air  [J/kg/K]
          & pcv=cvv                                ,&! in     specific heat of vapor    [J/kg/K]
          & pi0=prm_diag%flxdwswtoa(:,jb)          ,&! in     solar incoming flux at TOA  [W/m2]
          & pemiss=prm_diag%lw_emiss(:,jb)         ,&! in     lw sfc emissivity
          & pqc=prm_diag%tot_cld    (:,:,jb,iqc)   ,&! in     specific cloud water        [kg/kg]
          & pqi=prm_diag%tot_cld    (:,:,jb,iqi)   ,&! in     specific cloud ice          [kg/kg]
          & ppres_ifc=pt_diag%pres_ifc(:,:,jb)     ,&! in     pressure at layer boundaries [Pa]
          & albedo=prm_diag%albdif(:,jb)           ,&! in     grid-box average shortwave albedo
          & albedo_t=prm_diag%albdif_t(:,jb,:)     ,&! in     tile-specific shortwave albedo
          & list_land_count  = ext_data%atm%list_land%ncount(jb),  &! in number of land points
          & list_land_idx    = ext_data%atm%list_land%idx(:,jb),   &! in index list of land points
          & list_seawtr_count= ext_data%atm%list_seawtr%ncount(jb),&! in number of water points
          & list_seawtr_idx  = ext_data%atm%list_seawtr%idx(:,jb), &! in index list of water points
          & list_seaice_count= ext_data%atm%list_seaice%ncount(jb),&! in number of seaice points
          & list_seaice_idx  = ext_data%atm%list_seaice%idx(:,jb), &! in index list of seaice points
          & list_lake_count  = ext_data%atm%list_lake%ncount(jb),  &! in number of (f)lake points
          & list_lake_idx    = ext_data%atm%list_lake%idx(:,jb),   &! in index list of (f)lake points
          & gp_count_t       = ext_data%atm%gp_count_t(jb,:),      &! in number of land points per tile
          & idx_lst_t        = ext_data%atm%idx_lst_t(:,jb,:),     &! in index list of land points per tile
          & cosmu0=zcosmu0(:,jb)                   ,&! in     cosine of solar zenith angle (w.r.t. plain surface)
          & cosmu0_slp=cosmu0_slope(:,jb)          ,&! in     slope-dependent cosine of solar zenith angle
          & skyview=ext_data%atm%skyview(:,jb)     ,&! in     skyview factor for islope_rad=2
          & opt_nh_corr=.TRUE.                     ,&! in     switch for NH mode
          & ptsfc=lnd_prog_new%t_g(:,jb)           ,&! in     surface temperature         [K]
          & ptsfc_t=lnd_prog_new%t_g_t(:,jb,:)     ,&! in     tile-specific surface temperature         [K]
          & ptsfctrad=prm_diag%tsfctrad(:,jb)      ,&! in     sfc temp. used for pflxlw   [K]
          & ptrmsw=prm_diag%trsolall (:,:,jb)      ,&! in     shortwave net tranmissivity []
          & pflxlw=prm_diag%lwflxall (:,:,jb)      ,&! in     longwave net flux           [W/m2]
          & lwflx_up_sfc_rs=prm_diag%lwflx_up_sfc_rs(:,jb), &! in longwave upward flux at surface [W/m2]
          & trsol_up_toa=prm_diag%trsol_up_toa(:,jb),   & ! in shortwave upward transm. at the top of the atmosphere
          & trsol_up_sfc=prm_diag%trsol_up_sfc(:,jb),   & ! in shortwave upward transm. at the surface
          & trsol_nir_sfc=prm_diag%trsol_nir_sfc(:,jb), & ! in near-infrared downward transm. at the surface
          & trsol_vis_sfc=prm_diag%trsol_vis_sfc(:,jb), & ! in visible downward transm. at the surface
          & trsol_par_sfc=prm_diag%trsol_par_sfc(:,jb), & ! in photosynthetically active downward transm. at the surface
          & trsol_dn_sfc_diff=prm_diag%trsol_dn_sfc_diff(:,jb),&! in shortwave diffuse downward transm. at the surface
          & trsol_clr_sfc=prm_diag%trsolclr_sfc(:,jb),  & ! in clear-sky net transmissivity at surface
          & use_trsolclr_sfc=.TRUE.                ,&     ! in use clear-sky surface transmissivity (optional)
          !
          ! output
          ! ------
          !
          & pdtdtradsw=prm_nwp_tend%ddt_temp_radsw(:,:,jb),&! out rad. heating by SW         [K/s]
          & pdtdtradlw=prm_nwp_tend%ddt_temp_radlw(:,:,jb),&! out rad. heating by LW         [K/s]
          & pflxsfcsw =prm_diag%swflxsfc (:,jb)        ,&   ! out shortwave surface net flux [W/m2]
          & pflxsfclw =prm_diag%lwflxsfc (:,jb)        ,&   ! out longwave surface net flux  [W/m2]
          & pflxsfcsw_t=prm_diag%swflxsfc_t (:,jb,:)   ,&   ! out tile-specific shortwave surface net flux [W/m2]
          & pflxsfclw_t=prm_diag%lwflxsfc_t (:,jb,:)   ,&   ! out tile-specific longwave surface net flux  [W/m2]
          & pflxtoasw =prm_diag%swflxtoa (:,jb)        ,&   ! out shortwave toa net flux     [W/m2]
          & pflxtoalw =prm_diag%lwflxtoa (:,jb)        ,&   ! out longwave  toa net flux     [W/m2]
          & lwflx_up_sfc=prm_diag%lwflx_up_sfc(:,jb)   ,&   ! out longwave upward flux at surface [W/m2]
          & swflx_up_toa=prm_diag%swflx_up_toa(:,jb)   ,&   ! out shortwave upward flux at the TOA [W/m2]
          & swflx_up_sfc=prm_diag%swflx_up_sfc(:,jb)   ,&   ! out shortwave upward flux at the surface [W/m2]
          & swflx_nir_sfc=prm_diag%swflx_nir_sfc(:,jb) ,&   ! out near-infrared downward flux at the surface [W/m2]
          & swflx_vis_sfc=prm_diag%swflx_vis_sfc(:,jb) ,&   ! out visible downward flux at the surface [W/m2]
          & swflx_par_sfc=prm_diag%swflx_par_sfc(:,jb) ,&   ! out PAR downward flux at the surface [W/m2]
          & swflx_clr_sfc=prm_diag%swflxclr_sfc(:,jb)  ,&   ! out clear-sky shortwave flux at the surface [W/m2]
          & swflx_dn_sfc_diff=prm_diag%swflx_dn_sfc_diff(:,jb), & ! out shortwave diffuse downward flux at the surface [W/m2]
          & lacc=lzacc                                          )

        ELSE
          CALL radheat (                   &
          !
          ! input
          ! -----
          !
          & jcs=i_startidx                         ,&! in     start index of inner do loop
          & jce=i_endidx                           ,&! in     end index of inner do loop
          & jg=pt_patch%id                         ,&! in     patch ID
          & kbdim=nproma                           ,&! in     loop length and dimension size
          & klev=nlev                              ,&! in     vertical dimension size
          & klevp1=nlevp1                          ,&! in     vertical dimension size
          & ntiles=1                               ,&! in     number of tiles of sfc flux fields
          & ntiles_wtr=0                           ,&! in     number of extra tiles for ocean and lakes
          & pmair=pt_diag%airmass_new(:,:,jb)      ,&! in     layer air mass             [kg/m2]
          & pqv=prm_diag%tot_cld(:,:,jb,iqv)       ,&! in     specific moisture           [kg/kg]
          & pcd=cvd                                ,&! in     specific heat of dry air  [J/kg/K]
          & pcv=cvv                                ,&! in     specific heat of vapor    [J/kg/K]
          & pi0=prm_diag%flxdwswtoa(:,jb)          ,&! in     solar incoming flux at TOA  [W/m2]
          & pemiss=prm_diag%lw_emiss(:,jb)         ,&! in     lw sfc emissivity
          & pqc=prm_diag%tot_cld    (:,:,jb,iqc)   ,&! in     specific cloud water        [kg/kg]
          & pqi=prm_diag%tot_cld    (:,:,jb,iqi)   ,&! in     specific cloud ice          [kg/kg]
          & ppres_ifc=pt_diag%pres_ifc(:,:,jb)     ,&! in     pressure at layer boundaries [Pa]
          & cosmu0=zcosmu0(:,jb)                   ,&! in     cosine of solar zenith angle
          & cosmu0_slp=cosmu0_slope(:,jb)          ,&! in     slope-dependent cosine of solar zenith angle
          & skyview=ext_data%atm%skyview(:,jb)     ,&! in     skyview factor for islope_rad=2
          & opt_nh_corr=.TRUE.                     ,&! in     switch for NH mode
          & ptsfc=lnd_prog_new%t_g(:,jb)           ,&! in     surface temperature         [K]
          & ptsfctrad=prm_diag%tsfctrad(:,jb)      ,&! in     sfc temp. used for pflxlw   [K]
          & ptrmsw=prm_diag%trsolall (:,:,jb)      ,&! in     shortwave net tranmissivity []
          & pflxlw=prm_diag%lwflxall (:,:,jb)      ,&! in     longwave net flux           [W/m2]
          & lwflx_up_sfc_rs=prm_diag%lwflx_up_sfc_rs(:,jb), &! in longwave upward flux at surface [W/m2]
          & trsol_up_toa=prm_diag%trsol_up_toa(:,jb),   & ! in shortwave upward transm. at the top of the atmosphere
          & trsol_up_sfc=prm_diag%trsol_up_sfc(:,jb),   & ! in shortwave upward transm. at the surface
          & trsol_nir_sfc=prm_diag%trsol_nir_sfc(:,jb), & ! in near-infrared downward transm. at the surface
          & trsol_vis_sfc=prm_diag%trsol_vis_sfc(:,jb), & ! in visible downward transm. at the surface
          & trsol_par_sfc=prm_diag%trsol_par_sfc(:,jb), & ! in photosynthetically active downward transm. at the surface
          & trsol_dn_sfc_diff=prm_diag%trsol_dn_sfc_diff(:,jb),&! in shortwave diffuse downward transm. at the surface
          !
          ! output
          ! ------
          !
          & pdtdtradsw=prm_nwp_tend%ddt_temp_radsw(:,:,jb),&! out    rad. heating by SW        [K/s]
          & pdtdtradlw=prm_nwp_tend%ddt_temp_radlw(:,:,jb),&! out    rad. heating by lw        [K/s]
          & pflxsfcsw =prm_diag%swflxsfc (:,jb)   ,&        ! out shortwave surface net flux [W/m2]
          & pflxsfclw =prm_diag%lwflxsfc (:,jb)   ,&        ! out longwave surface net flux  [W/m2]
          & pflxtoasw =prm_diag%swflxtoa (:,jb)   ,&        ! out shortwave toa net flux     [W/m2]
          & pflxtoalw =prm_diag%lwflxtoa (:,jb)   ,&        ! out longwave  toa net flux     [W/m2]
          & lwflx_up_sfc=prm_diag%lwflx_up_sfc(:,jb)   ,&   ! out longwave upward flux at surface [W/m2]
          & swflx_up_toa=prm_diag%swflx_up_toa(:,jb)   ,&   ! out shortwave upward flux at the TOA [W/m2]
          & swflx_up_sfc=prm_diag%swflx_up_sfc(:,jb)   ,&   ! out shortwave upward flux at the surface [W/m2]
          & swflx_nir_sfc=prm_diag%swflx_nir_sfc(:,jb) ,&   ! out near-infrared downward flux at the surface [W/m2]
          & swflx_vis_sfc=prm_diag%swflx_vis_sfc(:,jb) ,&   ! out visible downward flux at the surface [W/m2]
          & swflx_par_sfc=prm_diag%swflx_par_sfc(:,jb) ,&   ! out PAR downward flux at the surface [W/m2]
          & swflx_dn_sfc_diff=prm_diag%swflx_dn_sfc_diff(:,jb), & ! out shortwave diffuse downward flux at the surface [W/m2]
          & lacc=lzacc                                          )
        ENDIF

      ENDDO ! blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "radheat", .FALSE., opt_lupdate_cpu=.TRUE., opt_dt=mtime_datetime)

    !$ACC WAIT(1)
    !$ACC END DATA ! CREATE(cosmu0_slope)

      IF (timers_level > 2) CALL timer_stop(timer_radheat)

    ENDIF


    !-------------------------------------------------------------------------
    !> Gravity waves drag: orographic and non-orographic
    !-------------------------------------------------------------------------

    IF (lcall_phy_jg(itsso) .OR. lcall_phy_jg(itgwd)) THEN

      IF (msg_level >= 15) &
        &  CALL message('mo_nh_interface', 'gravity waves')

      IF (timers_level > 3) CALL timer_start(timer_sso)

      ! GZ: use fast-physics time step instead of dt_phy_jg(itsso) in order to avoid calling-frequency dependence of low-level blocking
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "gwdrag", .TRUE., opt_lupdate_cpu=.FALSE.)
      CALL nwp_gwdrag ( dt_loc,                    & !>input
        &               lcall_phy_jg(itsso),       & !>input
        &               dt_phy_jg(itgwd),          & !>input
        &               lcall_phy_jg(itgwd),       & !>input
        &               pt_patch, p_metrics,       & !>input
        &               ext_data,                  & !>input
        &               pt_diag,                   & !>inout
        &               prm_diag, prm_nwp_tend,    & !>inout
        &               lacc=lzacc                  ) !>in
      !$ser verbatim IF (.not. linit) CALL serialize_all(nproma, jg, "gwdrag", .FALSE., opt_lupdate_cpu=.FALSE.)

      IF (timers_level > 3) CALL timer_stop(timer_sso)
    ENDIF ! inwp_sso



    !-------------------------------------------------------------------------
    !> Waves coupling: if coupling time step
    !-------------------------------------------------------------------------

    IF ( is_coupled_to_waves() .AND. (.NOT. linit) ) THEN


    END IF


    !-------------------------------------------------------------------------
    !> Ocean coupling: if coupling time step (VDIFF calls this internally)
    !-------------------------------------------------------------------------

    IF ( is_coupled_to_ocean() .AND. (.NOT. linit) &
      & .AND. atm_phy_nwp_config(jg)%inwp_turb /= ivdiff) THEN

      IF (ltimer) CALL timer_start(timer_coupling)
      CALL nwp_couple_ocean( pt_patch, pt_diag, lnd_diag, wtr_prog_new, prm_diag, ext_data )

      !------------------------------------------------
      !  After receiving new sea-ice fraction (fr_seaice) and thickness (h_ice)
      !  update sea/seaice index lists
      !
      CALL process_sst_and_seaice( pt_patch, lnd_diag%fr_seaice, lnd_diag%t_seasfc, pt_diag%pres_sfc, &
        & ext_data, lnd_prog_now, lnd_prog_new, wtr_prog_now, wtr_prog_new, lnd_diag, wtr_prog_new%h_ice )

      IF (ltimer) CALL timer_stop(timer_coupling)

    END IF


    !-------------------------------------------------------------------------
    ! Anurag Dipankar MPIM (2013-May-29)
    ! Large-scale forcing is to be applied at the end of all physics so that
    ! the most updated variable is used. Ideally it should be "next" timestep
    ! variable. Also note that its not actually a part of physics (sub-grid
    ! activity). It is called here to take advantage of u,v.
    !
    ! These LS forcing act as slow process so the tendencies from them are
    ! accumulated with the slow physics tendencies next
    !
    !(2013-25-June) LS forcing is called every physics step
    !-------------------------------------------------------------------------
    IF(is_ls_forcing)THEN

      IF (msg_level >= 15) &
        &  CALL message('mo_nh_interface:', 'LS forcing')

      IF (timers_level > 3) CALL timer_start(timer_ls_forcing)

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_block(rl_start)
      i_endblk   = pt_patch%cells%end_block(rl_end)

      ! Call to apply_ls_forcing
      CALL apply_ls_forcing ( pt_patch,                         & !>in
        &                     p_metrics,                        & !>in
        &                     p_sim_time,                       & !>in
        &                     pt_prog,                          & !>in
        &                     pt_diag,                          & !>in
        &                     pt_prog_rcf%tracer(:,:,:,iqv),    & !>in
        &                     rl_start,                         & !>in
        &                     rl_end,                           & !>in
        &                     prm_nwp_tend%ddt_u_ls,            & !>out
        &                     prm_nwp_tend%ddt_v_ls,            & !>out
        &                     prm_nwp_tend%ddt_temp_ls,         & !>out
        &                     prm_nwp_tend%ddt_tracer_ls(:,iqv),& !>out
        &                     prm_nwp_tend%ddt_temp_subs_ls,    & !>out
        &                     prm_nwp_tend%ddt_qv_subs_ls,      & !>out
        &                     prm_nwp_tend%ddt_temp_adv_ls,     & !>out
        &                     prm_nwp_tend%ddt_qv_adv_ls,       & !>out
        &                     prm_nwp_tend%ddt_u_adv_ls,        & !>out
        &                     prm_nwp_tend%ddt_v_adv_ls,        & !>out
        &                     prm_nwp_tend%wsub,                & !>out
        &                     prm_nwp_tend%fc_sfc_lat_flx,      & !>out
        &                     prm_nwp_tend%fc_sfc_sens_flx,     & !>out
        &                     prm_nwp_tend%fc_ts,               & !>out
        &                     prm_nwp_tend%fc_tg,               & !>out
        &                     prm_nwp_tend%fc_qvs,              & !>out
        &                     prm_nwp_tend%fc_Ch,               & !>out
        &                     prm_nwp_tend%fc_Cq,               & !>out
        &                     prm_nwp_tend%fc_Cm,               & !>out
        &                     prm_nwp_tend%fc_ustar,            & !>out
        &                     prm_nwp_tend%temp_nudge,          & !>out
        &                     prm_nwp_tend%u_nudge,             & !>out
        &                     prm_nwp_tend%v_nudge,             & !>out
        &                     prm_nwp_tend%q_nudge              ) !>out

      !simplified raditation scheme for Stratocumulus cases
      IF (is_sim_rad) THEN
        DO jb = i_startblk, i_endblk
          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

          CALL  sim_rad(                          &
            & i_startidx,                         &  !>in
            & i_endidx,                           &  !>in
            & nproma,                             &  !>in
            & 1,                                  &  !>in
            & nlev,                               &  !>in
            & p_metrics%z_mc(:,:,jb),             &  !>in
            & p_metrics%z_ifc(:,:,jb),            &  !>in
            & pt_prog_rcf%tracer(:,:,jb,iqv),     &  !>in
            & pt_prog_rcf%tracer(:,:,jb,iqc),     &  !>in
            & pt_prog_rcf%tracer(:,:,jb,iqi),     &  !>in
            & pt_prog_rcf%rho(:,:,jb),            &  !>in
            & prm_nwp_tend%ddt_temp_sim_rad(:,:,jb)) !>inout
        ENDDO
      ENDIF

      IF (timers_level > 3) CALL timer_stop(timer_ls_forcing)

    ENDIF

    !-------------------------------------------------------------------------
    !  Upper-atmosphere physics: compute tendencies
    !-------------------------------------------------------------------------
    IF (l_any_upatmophys) THEN
      IF (upatmo_config(jg)%l_status( iUpatmoStat%timer )) CALL timer_start(timer_upatmo)
      ! This interface has to be called after all other slow physics.
      CALL nwp_upatmo_interface( dt_loc            = dt_loc,            & !in
        &                        mtime_datetime    = mtime_datetime,    & !in
        &                        p_patch           = pt_patch,          & !in
        &                        p_int_state       = pt_int_state,      & !in
        &                        p_metrics         = p_metrics,         & !in
        &                        p_prog            = pt_prog,           & !in
        &                        p_prog_rcf        = pt_prog_rcf,       & !in
        &                        p_diag            = pt_diag,           & !in
        &                        prm_nwp_diag      = prm_diag,          & !in
        &                        prm_nwp_tend      = prm_nwp_tend,      & !in
        &                        kstart_moist      = kstart_moist(jg)   ) !in
      IF (upatmo_config(jg)%l_status( iUpatmoStat%timer )) CALL timer_stop(timer_upatmo)

    ENDIF

    IF (timers_level > 2) CALL timer_start(timer_phys_acc)
    !-------------------------------------------------------------------------
    !>  accumulate tendencies of slow_physics: Not called when LS focing is ON
    !-------------------------------------------------------------------------
    IF( (l_any_slowphys .OR. lcall_phy_jg(itradheat)) .OR. is_ls_forcing) THEN
      ! needs to be always initialized with OpenACC
      IF (p_test_run .OR. lzacc) THEN
        CALL init(z_ddt_u_tot, opt_acc_async=.TRUE.)
        CALL init(z_ddt_v_tot, opt_acc_async=.TRUE.)
      ENDIF

      IF (timers_level > 10) CALL timer_start(timer_phys_acc_1)

      IF (msg_level >= 15) &
        &  CALL message('mo_nh_interface:', 'accumulate slow phys')

      ! Coefficients for extra Rayleigh friction
      ustart    = atm_phy_nwp_config(jg)%ustart_raylfric
      uoffset_q = ustart + 40._wp
      ustart_q  = ustart + 50._wp
      max_relax = -1._wp/atm_phy_nwp_config(jg)%efdt_min_raylfric
      !$ACC DATA COPYIN(ustart, uoffset_q, ustart_q, max_relax)

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_block(rl_start)
      i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_qsum,z_ddt_temp,z_ddt_alpha,vabs,nudgecoeff,&
!$OMP  rfric_fac,zddt_u_raylfric,zddt_v_raylfric,convfac,convind,sqrt_ri,n2,dvdz2,wfac) ICON_OMP_DEFAULT_SCHEDULE
!
      DO jb = i_startblk, i_endblk
!
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)


        ! artificial Rayleigh friction: active if GWD or SSO scheme is active
        IF (atm_phy_nwp_config(jg)%inwp_sso > 0 .OR. atm_phy_nwp_config(jg)%inwp_gwd > 0) THEN
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR PRIVATE(vabs, rfric_fac) COLLAPSE(2)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              vabs = SQRT(pt_diag%u(jc,jk,jb)**2 + pt_diag%v(jc,jk,jb)**2)
              rfric_fac = MAX(0._wp, 8.e-4_wp*(vabs-ustart))
              IF (vabs > ustart_q) THEN
                rfric_fac = MIN(1._wp,4.e-4_wp*(vabs-uoffset_q)**2)
              ENDIF
              zddt_u_raylfric(jc,jk) = max_relax*rfric_fac*pt_diag%u(jc,jk,jb)
              zddt_v_raylfric(jc,jk) = max_relax*rfric_fac*pt_diag%v(jc,jk,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ELSE
          !$ACC KERNELS ASYNC(1) IF(lzacc)
          zddt_u_raylfric(:,:) = 0._wp
          zddt_v_raylfric(:,:) = 0._wp
          !$ACC END KERNELS
        ENDIF


        ! Accumulate wind tendencies of slow physics
        ! Strictly spoken, this would not be necessary if only radiation was called
        ! in the current time step, but the radiation time step should be a multiple
        ! of the convection time step anyway in order to obtain up-to-date cloud cover fields
        IF (l_any_slowphys) THEN
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              z_ddt_u_tot(jc,jk,jb) = prm_nwp_tend%ddt_u_gwd(jc,jk,jb) + zddt_u_raylfric(jc,jk)  &
                &     + prm_nwp_tend%ddt_u_sso(jc,jk,jb)  + prm_nwp_tend%ddt_u_pconv(jc,jk,jb)
              z_ddt_v_tot(jc,jk,jb) = prm_nwp_tend%ddt_v_gwd(jc,jk,jb) + zddt_v_raylfric(jc,jk)  &
                &     + prm_nwp_tend%ddt_v_sso(jc,jk,jb)  + prm_nwp_tend%ddt_v_pconv(jc,jk,jb)
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ELSE IF (is_ls_forcing) THEN
          z_ddt_u_tot(i_startidx:i_endidx,:,jb) = 0._wp
          z_ddt_v_tot(i_startidx:i_endidx,:,jb) = 0._wp
        ENDIF


        ! SQRT of Richardson number between the two lowest model levels
        ! This is used below to reduce frictional heating near the surface under very stable conditions
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(n2, dvdz2)
        DO jc = i_startidx, i_endidx
          n2 = 2._wp*grav/(pt_prog%theta_v(jc,nlev,jb)+pt_prog%theta_v(jc,nlev-1,jb)) * MAX(1.e-4_wp,        &
               (pt_prog%theta_v(jc,nlev-1,jb)-pt_prog%theta_v(jc,nlev,jb))/p_metrics%ddqz_z_half(jc,nlev,jb) )
          dvdz2 = MAX(1.e-6_wp, ( (pt_diag%u(jc,nlev-1,jb)-pt_diag%u(jc,nlev,jb))**2 +                      &
                  (pt_diag%v(jc,nlev-1,jb)-pt_diag%v(jc,nlev,jb))**2 )/p_metrics%ddqz_z_half(jc,nlev,jb)**2 )
          sqrt_ri(jc) = MAX(1._wp, SQRT(n2/dvdz2))
        ENDDO
        !$ACC END PARALLEL


        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(wfac) COLLAPSE(2)
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            !
            ! heating related to momentum deposition by SSO, GWD and Rayleigh friction
            !
            wfac = MIN(1._wp, 0.004_wp*p_metrics%geopot_agl(jc,jk,jb)/grav)
            prm_nwp_tend%ddt_temp_drag(jc,jk,jb) = -rcvd*(pt_diag%u(jc,jk,jb)*             &
                                                   (prm_nwp_tend%ddt_u_sso(jc,jk,jb)+      &
                                                    prm_nwp_tend%ddt_u_gwd(jc,jk,jb)+      &
                                                    zddt_u_raylfric(jc,jk))                &
                                                   +      pt_diag%v(jc,jk,jb)*             &
                                                   (prm_nwp_tend%ddt_v_sso(jc,jk,jb)+      &
                                                    prm_nwp_tend%ddt_v_gwd(jc,jk,jb)+      &
                                                    zddt_v_raylfric(jc,jk)) ) /            &
                                                    ((1._wp-wfac)*sqrt_ri(jc) + wfac)
            !
            ! total slow physics heating rate
            !
            z_ddt_temp(jc,jk) = prm_nwp_tend%ddt_temp_radsw(jc,jk,jb) + prm_nwp_tend%ddt_temp_radlw(jc,jk,jb) &
              &              +  prm_nwp_tend%ddt_temp_drag (jc,jk,jb) + prm_nwp_tend%ddt_temp_pconv(jc,jk,jb) &
              &              +  prm_nwp_tend%ddt_temp_clcov(jc,jk,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL


        IF(sppt_config(jg)%lsppt .AND. .NOT. linit) THEN
          !
          ! Add tendencies from the fast physics to the ones derived from the slow pysics as well as perturbations
          !

          CALL pert_tend(jb, jg, i_startidx, i_endidx, pt_patch%nlev, sppt(jg),           &
                         prm_nwp_tend, pt_prog%rho(:,:,jb),                               &
                         z_ddt_temp(:,:), z_ddt_u_tot(:,:,jb), z_ddt_v_tot(:,:,jb), lacc=lzacc)
        ENDIF ! end of lsppt


        CALL calc_qsum (pt_prog_rcf%tracer, z_qsum, condensate_list, jb, i_startidx, i_endidx, 1, kstart_moist(jg), nlev)


        IF (kstart_moist(jg) > 1) THEN
          !$ACC KERNELS ASYNC(1) IF(lzacc)
          z_ddt_alpha(:,1:kstart_moist(jg)-1) = 0._wp
          !$ACC END KERNELS
        ENDIF

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            ! tendency of virtual increment
            ! tendencies of iqr,iqs are neglected (nonzero only for ldetrain_conv_prec=.TRUE.)
            z_ddt_alpha(jc,jk) = ( vtmpc1 * prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqv) &
             &                 - prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc)            &
             &                 - prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi) )          &
             &                 / pt_prog%rho(jc,jk,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        ! Convert temperature tendency into Exner function tendency
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            pt_diag%ddt_exner_phy(jc,jk,jb) = rd_o_cpd / pt_prog%theta_v(jc,jk,jb)           &
              &                             * (z_ddt_temp(jc,jk)                             &
              &                             *(1._wp + vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)&
              &                             - z_qsum(jc,jk))                                 &
              &                             + pt_diag%temp(jc,jk,jb) * z_ddt_alpha(jc,jk))

          ENDDO
        ENDDO
        !$ACC END PARALLEL


        !-------------------------------------------------------------------------
        !>  accumulate tendencies of slow_physics when LS forcing is ON
        !-------------------------------------------------------------------------
        IF (is_ls_forcing) THEN
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx

              ! add u/v/T forcing tendency
              z_ddt_u_tot(jc,jk,jb)   = z_ddt_u_tot(jc,jk,jb)       &
                &                     + prm_nwp_tend%ddt_u_ls(jk)

              z_ddt_v_tot(jc,jk,jb)   = z_ddt_v_tot(jc,jk,jb)       &
                &                     + prm_nwp_tend%ddt_v_ls(jk)

              z_ddt_temp(jc,jk)       = z_ddt_temp(jc,jk)           &
                                      + prm_nwp_tend%ddt_temp_ls(jk)

              ! simplified radiation scheme
              IF(is_sim_rad) THEN
                z_ddt_temp(jc,jk)     = z_ddt_temp(jc,jk)           &
                                      + prm_nwp_tend%ddt_temp_sim_rad(jc,jk,jb)
              ENDIF

              ! linear nudging profile between "start" and "full" heights - prevent sfc layer instability
              IF ( nudge_full_height == nudge_start_height ) THEN
                nudgecoeff = 1.0_wp
              ELSE
                nudgecoeff = ( p_metrics%geopot_agl(jc,jk,jb)/grav - nudge_start_height ) / &
                           & ( nudge_full_height                   - nudge_start_height )
                nudgecoeff = MAX( MIN( nudgecoeff, 1.0_wp ), 0.0_wp )
              END IF

              ! add u/v/T nudging
              IF ( is_nudging_uv ) THEN
                ! explicit:          (u,n+1 - u,n) / dt = (u,nudge - u,n)   / dt_relax

                ! implicit:          (u,n+1 - u,n) / dt = (u,nudge - u,n+1) / dt_relax

                ! analytic implicit: (u,n+1 - u,n) / dt = (u,nudge - u,n)   / dt_relax * exp(-dt/dt_relax)
                z_ddt_u_tot(jc,jk,jb) = z_ddt_u_tot(jc,jk,jb)       &
                  &  - ( pt_diag%u(jc,jk,jb) - prm_nwp_tend%u_nudge(jk) ) / dt_relax * exp(-dt_loc/dt_relax) &
                  &  * nudgecoeff

                z_ddt_v_tot(jc,jk,jb) = z_ddt_v_tot(jc,jk,jb)       &
                  &  - ( pt_diag%v(jc,jk,jb) - prm_nwp_tend%v_nudge(jk) ) / dt_relax * exp(-dt_loc/dt_relax) &
                  &  * nudgecoeff

              END IF

              ! attention: T nudging results in dt-step oscillation/instability!
              ! q nudging done in tracer_add_phytend in mo_util_phys

              IF ( is_nudging_tq ) THEN
                ! explicit:          (T,n+1 - T,n) / dt = (T,nudge - T,n)   / dt_relax

                ! implicit:          (T,n+1 - T,n) / dt = (T,nudge - T,n+1) / dt_relax

                ! analytic implicit: (T,n+1 - T,n) / dt = (T,nudge - T,n)   / dt_relax * exp(-dt/dt_relax)
                z_ddt_temp(jc,jk)     = z_ddt_temp(jc,jk)           &
                  &  - ( pt_diag%temp(jc,jk,jb) - prm_nwp_tend%temp_nudge(jk) ) / dt_relax &
                  &  * exp(-dt_loc/dt_relax) * nudgecoeff
              END IF

              ! Convert temperature tendency into Exner function tendency
              z_ddt_alpha(jc,jk) = z_ddt_alpha(jc,jk)                          &
                &                + vtmpc1 * prm_nwp_tend%ddt_tracer_ls(jk,iqv) &
                &                - prm_nwp_tend%ddt_tracer_ls(jk,iqc)          &
                &                - prm_nwp_tend%ddt_tracer_ls(jk,iqi)

              pt_diag%ddt_exner_phy(jc,jk,jb) = rd_o_cpd / pt_prog%theta_v(jc,jk,jb)           &
                &                             * (z_ddt_temp(jc,jk)                             &
                &                             *(1._wp + vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)&
                &                             - z_qsum(jc,jk))                                 &
                &                             + pt_diag%temp(jc,jk,jb) * z_ddt_alpha(jc,jk) )

            END DO  ! jc
          END DO  ! jk

        ENDIF ! END of LS forcing tendency accumulation

        IF (lcall_phy_jg(itconv)) THEN
!DIR$ IVDEP
          !$ACC DATA PRESENT(pt_diag%temp)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR PRIVATE(convfac, convind)
          DO jc = i_startidx, i_endidx
            convind = prm_diag%k950(jc,jb) ! using this varibale directly in the next row gave a memory "memory not mapped to object" error in PGI (GPU) 20.8
            ! rain-snow conversion factor to avoid 'snow showers' at temperatures when they don't occur in practice
            convfac = MIN(1._wp,MAX(0._wp,pt_diag%temp(jc,convind,jb)-tmelt)* &
              MAX(0._wp,prm_diag%t_2m(jc,jb)-(tmelt+1.5_wp)) )
            prm_diag%rain_con_rate(jc,jb) = prm_diag%rain_con_rate_3d(jc,nlevp1,jb) + &
              convfac*prm_diag%snow_con_rate_3d(jc,nlevp1,jb)
            prm_diag%snow_con_rate(jc,jb) = (1._wp-convfac)*prm_diag%snow_con_rate_3d(jc,nlevp1,jb)
          ENDDO
          !$ACC END PARALLEL
          !$ACC END DATA
        ENDIF

      ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      IF (timers_level > 10) CALL timer_stop(timer_phys_acc_1)

    !$ACC END DATA
    END IF  !END OF slow physics tendency accumulation


    ! Add only perturbed tendencies of tracers to the corresponding variable.
    ! Adding tendencies was taking care of during fast physics already.
    IF (sppt_config(jg)%lsppt .AND. .NOT. linit) THEN
      CALL apply_tend(pt_patch,sppt(jg), pt_prog_rcf, dt_loc, lacc=lzacc)
    ENDIF ! end of lsppt

    !--------------------------------------------------------
    ! Final section: Synchronization of updated prognostic variables,
    !                interpolation of u/v tendendies to edge points,
    !                and diagnostic computations
    !--------------------------------------------------------

    ! Synchronize tracers if any of the updating (fast-physics) processes was active.
    ! In addition, tempv and exner_pr needs to be synchronized.
    IF (advection_config(jg)%iadv_tke == 1) THEN
      ! TKE does not need to be synchronized if it is advected only vertically
      ntracer_sync = ntracer-1
    ELSE
      ntracer_sync = ntracer
    ENDIF

    IF (l_any_fastphys) THEN

      IF (timers_level > 10) CALL timer_start(timer_phys_sync_tracers)

      IF (diffusion_config(jg)%lhdiff_w .AND. iprog_aero >= 1) THEN
        CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+4, pt_diag%tempv, pt_prog%w, &
                                   pt_diag%exner_pr, prm_diag%aerosol,                         &
                                   f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))
      ELSE IF (diffusion_config(jg)%lhdiff_w) THEN
        CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+3, pt_diag%tempv, pt_prog%w, &
                                   pt_diag%exner_pr, f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))
      ELSE
        CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer_sync+2, pt_diag%tempv, &
                                   pt_diag%exner_pr, f4din=pt_prog_rcf%tracer(:,:,:,1:ntracer_sync))
      ENDIF

      IF (timers_level > 10) THEN
        CALL timer_stop(timer_phys_sync_tracers)
      ENDIF
    ELSE IF (linit .AND. advection_config(jg)%iadv_tke == 2) THEN
      CALL sync_patch_array(SYNC_C, pt_patch, pt_prog_rcf%tracer(:,:,:,iqtke))
    ENDIF


    !------------------------------------------------------------
    ! sync here the slowphys for aggregation
    !-------------------------------------------------------------------
    IF (use_physics_barrier) THEN
      CALL timer_start(timer_barrier)
      CALL work_mpi_barrier()
      CALL timer_stop(timer_barrier)
    ENDIF
    !-------------------------------------------------------------------
    IF (timers_level > 10) CALL timer_start(timer_phys_sync_ddt_u)

    IF ( (is_ls_forcing .OR. l_any_slowphys) .AND. lcall_phy_jg(itturb) ) THEN
      CALL sync_patch_array_mult(SYNC_C1, pt_patch, 4, z_ddt_u_tot, z_ddt_v_tot, &
                                 prm_nwp_tend%ddt_u_turb, prm_nwp_tend%ddt_v_turb)

    ELSE IF (lcall_phy_jg(itturb) ) THEN

      CALL sync_patch_array_mult(SYNC_C1, pt_patch, 2, prm_nwp_tend%ddt_u_turb, &
                                 prm_nwp_tend%ddt_v_turb)
    ENDIF

    ! DA: TODO: make kernels async in the interface and remove the wait
    !$ACC WAIT(1)

    IF (timers_level > 10) CALL timer_stop(timer_phys_sync_ddt_u)
    !------------------------------------------------------------


    !------------------------------------------------------------
    ! compute on the halos
    IF (timers_level > 10) CALL timer_start(timer_phys_acc_par)
    IF (l_any_fastphys) THEN
      IF (my_process_is_mpi_all_parallel() ) THEN

        rl_start = min_rlcell_int-1
        rl_end   = min_rlcell

        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end )

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc =  i_startidx, i_endidx

              IF (p_metrics%mask_prog_halo_c(jc,jb)) THEN
                pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
                  &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

                pt_prog%theta_v(jc,jk,jb) = pt_diag%tempv(jc,jk,jb) &
                  &                       / pt_prog%exner(jc,jk,jb)

              ENDIF

            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ENDIF ! my_process_is_mpi_all_parallel
    ENDIF ! fast-physics synchronization
    IF (timers_level > 10) CALL timer_stop(timer_phys_acc_par)

    !-------------------------------------------------------------------------
    !>
    !!    Interpolation from  u,v onto v_n
    !!      ddt_vn_phy  =interpol(ddt_u_tot)+interpol(ddt_v_tot)
    !!      Calculate normal velocity at edge midpoints
    !-------------------------------------------------------------------------

    IF (timers_level > 10)  CALL timer_start(timer_phys_acc_2)
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

    i_startblk = pt_patch%edges%start_block(rl_start)
    i_endblk   = pt_patch%edges%end_block(rl_end)


!$OMP DO PRIVATE(jb,jk,jce,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF ( (is_ls_forcing .OR. l_any_slowphys) .AND. lcall_phy_jg(itturb) ) THEN

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO jce = i_startidx, i_endidx

            pt_diag%ddt_vn_phy(jce,jk,jb) =   pt_int_state%c_lin_e(jce,1,jb)           &
&                                 * ( z_ddt_u_tot(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                                   + z_ddt_v_tot(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                                 * ( z_ddt_u_tot(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                                  +  z_ddt_v_tot(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

            pt_prog%vn(jce,jk,jb) = pt_prog%vn(jce,jk,jb) + dt_loc * (                 &
                                              pt_int_state%c_lin_e(jce,1,jb)           &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                       + prm_nwp_tend%ddt_v_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                      +  prm_nwp_tend%ddt_v_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                  *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 ) )

          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ELSE IF (lcall_phy_jg(itturb) ) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=8
        DO jk = 1, nlev
          DO jce = i_startidx, i_endidx

            pt_prog%vn(jce,jk,jb) = pt_prog%vn(jce,jk,jb) + dt_loc * (                 &
                                              pt_int_state%c_lin_e(jce,1,jb)           &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                       + prm_nwp_tend%ddt_v_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                      +  prm_nwp_tend%ddt_v_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                  *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 ) )

          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ENDIF

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    IF (timers_level > 10) CALL timer_stop(timer_phys_acc_2)

    !-------------------------------------------------------------------------
    !  Upper-atmosphere physics: add tendencies
    !-------------------------------------------------------------------------
    IF (l_any_upatmophys) THEN

      IF (upatmo_config(jg)%l_status( iUpatmoStat%timer )) CALL timer_start(timer_upatmo)
      ! This interface has to be called after all other tendencies have been accumulated.
      CALL nwp_upatmo_update( lslowphys         = l_any_slowphys,                  & !in
        &                     lradheat          = lcall_phy_jg(itradheat),         & !in
        &                     lturb             = lcall_phy_jg(itturb) .OR. linit, & !in
        &                     dt_loc            = dt_loc,                          & !in
        &                     p_patch           = pt_patch,                        & !inout
      	&                     p_prog_rcf        = pt_prog_rcf,                     & !inout
        &                     p_diag            = pt_diag                          ) !inout
      IF (upatmo_config(jg)%l_status( iUpatmoStat%timer )) CALL timer_stop(timer_upatmo)
    ENDIF

    IF (timers_level > 10) CALL timer_start(timer_phys_dpsdt)
    !
    ! dpsdt diagnostic
    IF (lcalc_dpsdt) THEN
      CALL compute_dpsdt (pt_patch      = pt_patch,  &
        &                 dt            = dt_loc,    &
        &                 pt_diag       = pt_diag,   &
        &                 lacc          = lzacc      )
    ENDIF
    IF (timers_level > 10) CALL timer_stop(timer_phys_dpsdt)


    IF (timers_level > 10) CALL timer_start(timer_phys_sync_vn)
    IF (lcall_phy_jg(itturb)) THEN
      CALL sync_patch_array(SYNC_E, pt_patch, pt_prog%vn)
    ENDIF
    IF (timers_level > 10) CALL timer_stop(timer_phys_sync_vn)
    IF (timers_level > 2) CALL timer_stop(timer_phys_acc)



    IF (lcall_phy_jg(itturb) .AND. msg_level >= 18) THEN ! extended diagnostic for turbulence quantities
      CALL nwp_diag_output_2(pt_patch, pt_prog_rcf, prm_nwp_tend)
    ENDIF

    CALL nwp_opt_diagnostics_2(pt_patch,             & !in
      &                        p_metrics,            & !in
      &                        pt_prog, pt_prog_rcf, & !in
      &                        pt_diag,              & !in
      &                        prm_diag,             & !inout
      &                        zcosmu0,              & !in
      &                        zsct,                 & !in
      &                        p_sim_time,           & !in
      &                        dt_phy_jg(itfastphy), & !in
      &                        lacc=lzacc             )

    ! time averages, accumulations and vertical integrals
    CALL nwp_statistics(lcall_phy_jg,                    & !in
                        & dt_phy_jg,p_sim_time,          & !in
                        & ext_data, kstart_moist(jg),    & !in
                        & ih_clch(jg), ih_clcm(jg),      & !in
                        & pt_patch, p_metrics,           & !in
                        & pt_prog, pt_prog_rcf,          & !in
                        & pt_diag,                       & !inout
                        & prm_diag, lnd_diag,            & !inout
                        & lacc=lzacc                      ) !in


    ! SBM microphysics
    ! store temperature snapshot at the end of a timestep
    IF (atm_phy_nwp_config(jg)%inwp_gscp == 8) THEN
      ptr_sbm_storage => get_sbm_storage(patch_id = jg)
!$OMP PARALLEL
      CALL copy(pt_prog_rcf%tracer(:,:,:,iqv), ptr_sbm_storage%qv_old )
      CALL copy(pt_diag%temp(:,:,:), ptr_sbm_storage%temp_old         )
!$OMP END PARALLEL
    ENDIF


    IF (ltimer) CALL timer_stop(timer_physics)

    !$ACC WAIT(1)
    !$ACC END DATA ! copyin
    !$ACC END DATA ! create

  END SUBROUTINE nwp_nh_interface

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------


END MODULE mo_nh_interface_nwp

