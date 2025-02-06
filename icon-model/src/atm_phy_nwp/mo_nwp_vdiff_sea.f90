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

! Sea and sea-ice handling for the NWP VDIFF interface.
!
! \note Do not use directly, use mo_nwp_vdiff_interface instead.

MODULE mo_nwp_vdiff_sea

  USE mo_albedo, ONLY: sfc_albedo_dir_rg, sfc_albedo_dir_taylor, sfc_albedo_dir_yang, &
      & sfc_albedo_whitecap
  USE mo_exception, ONLY: message, message_text
  USE mo_ext_data_init, ONLY: interpol_monthly_mean
  USE mo_ext_data_types, ONLY: t_external_data
  USE mo_fortran_tools, ONLY: copy, init, set_acc_host_or_device
  USE mo_idx_list, ONLY: t_idx_list_blocked
  USE mo_impl_constants, ONLY: SSTICE_ANA_CLINC, end_prog_cells, start_prog_cells
  USE mo_kind, ONLY: wp
  USE mo_lnd_nwp_config, ONLY: frsi_min, hice_min, hice_max, lprog_albsi, lseaice, sstice_mode
  USE mo_loopindices, ONLY: get_indices_c
  USE mo_master_config, ONLY: isRestart
  USE mo_model_domain, ONLY: t_patch
  USE mo_nwp_lnd_types, ONLY: t_lnd_diag, t_wtr_prog
  USE mo_nwp_ocean_coupling, ONLY: t_nwp_ocean_fields_rx, t_nwp_ocean_fields_tx, couple_ocean
  USE mo_nwp_vdiff_radfluxes, ONLY: t_nwp_vdiff_surface_rad_fluxes
  USE mo_nwp_vdiff_types, ONLY: t_nwp_vdiff_albedos, t_nwp_vdiff_sea_state, SFT_SICE, SFT_SWTR
  USE mo_parallel_config, ONLY: nproma
  USE mo_physical_constants, ONLY: als, alv, cpd, salinity_fac, stbo, tf_fresh => tmelt, tf_salt, &
      & zemiss_def
  USE mo_radiation_config, ONLY: albedo_whitecap, direct_albedo_water
  USE mo_run_config, ONLY: msg_level
  USE mo_satad, ONLY: sat_pres_ice, sat_pres_water, spec_humi
  USE mo_sync, ONLY: global_sum
  USE mo_timer, ONLY: ltimer, timer_coupling, timer_start, timer_stop
  USE mo_turb_vdiff, ONLY: vdiff_mixed_time_value, vdiff_surface_flux

  USE mtime, ONLY: datetime

  USE sfc_terra_data, ONLY: csalb, ist_seawtr, ist_seaice
  USE sfc_seaice, ONLY: alb_seaice_equil, seaice_init_nwp, seaice_timestep_nwp

#ifdef __NVCOMPILER
  USE mo_coupling_config, ONLY: is_coupled_to_ocean
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sea_model
  PUBLIC :: sea_model_init
  PUBLIC :: sea_model_update_sst
  PUBLIC :: sea_model_couple_ocean

  PUBLIC :: nwp_vdiff_update_seaice
  PUBLIC :: nwp_vdiff_update_seaice_list

  !
  ! Parameters
  !

  !> Ice emissivity. \todo Put this into `mo_physical_constants`.
  REAL(wp), PARAMETER :: lw_emissivity_ice = 0.975_wp

CONTAINS

  !>
  !! Sea and sea-ice model interface. There is no actual sea model, surface temperatures are
  !! obtained from external data (\c diag_lnd\%t_seasfc). The behavior mimics the one found in
  !! the nwp+terra scheme: Cells are considered ice-free once ice thickness drops below `hice_min`.
  !! It is not possible to re-freeze a cell, but the sea-ice fraction will be reset externally at
  !! midnight each day by `mo_td_ext_data::update_nwp_phy_bcs`. If the sea-ice scheme is switched
  !! off, every sea point is considered open water.
  !!
  SUBROUTINE sea_model ( &
        & iblk, ics, ice, dtime, steplen, ext_data, rain, snow, latent_hflx_ice_old, &
        & sensible_hflx_ice_old, flx_rad, sea_state, cos_zenith_angle, press_srf, wind_10m, &
        & prefactor_exchange, exchange_coeff_h_wtr, exchange_coeff_h_ice, t_acoef_wtr, &
        & t_bcoef_wtr, q_acoef_wtr, q_bcoef_wtr, t_acoef_ice, t_bcoef_ice, q_acoef_ice, &
        & q_bcoef_ice, prog_wtr_now, diag_lnd, t_wtr, t_ice, s_wtr, s_ice, qsat_wtr, qsat_ice, &
        & evapo_wtr, evapo_ice, latent_hflx_wtr, latent_hflx_ice, sensible_hflx_wtr, &
        & sensible_hflx_ice, conductive_hflx_ice, melt_potential_ice, alb, prog_wtr_new &
      )

    !> Current block index.
    INTEGER, INTENT(IN) :: iblk
    !> Cell index start.
    INTEGER, INTENT(IN) :: ics
    !> Cell index end.
    INTEGER, INTENT(IN) :: ice
    !> Length of time step [s].
    REAL(wp), INTENT(IN) :: dtime
    !> Length of time step in `prefactor_exchange` [s].
    REAL(wp), INTENT(IN) :: steplen
    !> External data.
    TYPE(t_external_data), TARGET, INTENT(IN) :: ext_data

    !> Rain at surface [kg/m**2/s].
    REAL(wp), INTENT(IN) :: rain(:)
    !> Snow at surface [kg/m**2/s].
    REAL(wp), INTENT(IN) :: snow(:)

    !> Latent heat flux into ice surface at time `t` [W/m**2].
    REAL(wp), INTENT(IN) :: latent_hflx_ice_old(:)
    !> Sensible heat flux into ice surface at time `t` [W/m**2].
    REAL(wp), INTENT(IN) :: sensible_hflx_ice_old(:)

    !> Radiation fluxes at surface [W/m**2].
    TYPE(t_nwp_vdiff_surface_rad_fluxes), INTENT(IN) :: flx_rad

    !> State of the model.
    TYPE(t_nwp_vdiff_sea_state), INTENT(INOUT) :: sea_state

    !> Cosine of the zenith angle [1].
    REAL(wp), INTENT(IN) :: cos_zenith_angle(:)
    !> Surface pressure [Pa].
    REAL(wp), INTENT(IN) :: press_srf(:)

    !> Wind speed at 10m [m/s].
    REAL(wp), INTENT(IN) :: wind_10m(:)

    !> Prefactor for exchange coefficients.
    REAL(wp), INTENT(IN) :: prefactor_exchange(:)
    !> Exchange coefficient for heat at water surface.
    REAL(wp), INTENT(IN) :: exchange_coeff_h_wtr(:)
    !> Exchange coefficient for heat at ice surface.
    REAL(wp), INTENT(IN) :: exchange_coeff_h_ice(:)

    !> Richtmyer E coefficients for heat at water surface [1].
    REAL(wp), INTENT(IN) :: t_acoef_wtr(:)
    !> Richtmyer F coefficients for heat at water surface [J/kg].
    REAL(wp), INTENT(IN) :: t_bcoef_wtr(:)
    !> Richtmyer E coefficients for specific humidity at water surface [1].
    REAL(wp), INTENT(IN) :: q_acoef_wtr(:)
    !> Richtmyer F coefficients for specific humidity at water surface [kg/kg].
    REAL(wp), INTENT(IN) :: q_bcoef_wtr(:)

    !> Richtmyer E coefficients for heat at ice surface [1].
    REAL(wp), INTENT(IN) :: t_acoef_ice(:)
    !> Richtmyer F coefficients for heat at ice surface [J/kg].
    REAL(wp), INTENT(IN) :: t_bcoef_ice(:)
    !> Richtmyer E coefficients for specific humidity at ice surface [1].
    REAL(wp), INTENT(IN) :: q_acoef_ice(:)
    !> Richtmyer F coefficients for specific humidity at ice surface [kg/kg].
    REAL(wp), INTENT(IN) :: q_bcoef_ice(:)

    !> Prognostic water variables at current time step.
    TYPE(t_wtr_prog), INTENT(IN) :: prog_wtr_now

    !> Diagnostic land variables. Reads `t_seasfc` and updates `fr_seaice`.
    TYPE(t_lnd_diag), INTENT(INOUT) :: diag_lnd

    !> Water temperature at time `t+1` [K].
    REAL(wp), INTENT(INOUT) :: t_wtr(:)
    !> Ice temperature at time `t+1` [K].
    REAL(wp), INTENT(INOUT) :: t_ice(:)

    !> Dry static energy over water at time `t+1` [J/kg].
    REAL(wp), INTENT(INOUT) :: s_wtr(:)
    !> Dry static energy over ice at time `t+1` [J/kg].
    REAL(wp), INTENT(INOUT) :: s_ice(:)

    !> Saturation specific humidity over water at time `t+1` [kg/kg].
    REAL(wp), INTENT(INOUT) :: qsat_wtr(:)
    !> Saturation specific humidity over ice at time `t+1` [kg/kg].
    REAL(wp), INTENT(INOUT) :: qsat_ice(:)

    !> Evaporation flux into water surface (mixed time) [kg/m**2/s].
    REAL(wp), INTENT(INOUT) :: evapo_wtr(:)
    !> Evaporation flux into ice surface (mixed time) [kg/m**2/s].
    REAL(wp), INTENT(INOUT) :: evapo_ice(:)

    !> Latent heat flux into water surface (mixed time) [W/m**2].
    REAL(wp), INTENT(INOUT) :: latent_hflx_wtr(:)
    !> Latent heat flux into ice surface (mixed time) [W/m**2].
    REAL(wp), INTENT(INOUT) :: latent_hflx_ice(:)

    !> Sensible heat flux into water surface (mixed time) [W/m**2].
    REAL(wp), INTENT(INOUT) :: sensible_hflx_wtr(:)
    !> Sensible heat flux into ice surface (mixed time) [W/m**2].
    REAL(wp), INTENT(INOUT) :: sensible_hflx_ice(:)
    !> Conductive heat flux at ice-ocean boundary [W/m**2].
    REAL(wp), OPTIONAL, INTENT(INOUT) :: conductive_hflx_ice(:)
    !> Melt potential at ice-atmosphere boundary [W/m**2].
    REAL(wp), OPTIONAL, INTENT(INOUT) :: melt_potential_ice(:)

    !> Albedos at time `t+1` sets sea water and ice indices [1].
    TYPE(t_nwp_vdiff_albedos), INTENT(INOUT) :: alb

    !> Prognostic water variables at time `t+1`.
    TYPE(t_wtr_prog), INTENT(INOUT) :: prog_wtr_new

    !
    ! Local variables
    !

    !> Mixed-time dry static energy over water [J/kg].
    REAL(wp) :: s_hat_wtr(SIZE(t_acoef_wtr))
    !> Mixed-time dry static energy over ice [J/kg].
    REAL(wp) :: s_hat_ice(SIZE(t_acoef_wtr))

    !> Mixed-time saturation specific humidity over water [kg/kg].
    REAL(wp) :: qsat_hat_wtr(SIZE(t_acoef_wtr))
    !> Mixed-time saturation specific humidity over ice [kg/kg].
    REAL(wp) :: qsat_hat_ice(SIZE(t_acoef_wtr))

    REAL(wp) :: t_ice_old !< Old ice temperature, sanitized [K].
    REAL(wp) :: qsat_ice_now !< Saturation specific humidity at time `t`.

    REAL(wp) :: wc_frac !< Whitecap fraction.
    REAL(wp) :: wc_alb !< Whitecap albedo.

    INTEGER :: ic !< Index list index.
    INTEGER :: jc !< Cell index.

    ! Temporaries for seaice scheme

    REAL(wp) :: qsen(SIZE(t_acoef_wtr)) !< Sensible heat flow [W/m**2].
    REAL(wp) :: qlat(SIZE(t_acoef_wtr)) !< Latent heat flow [W/m**2].
    REAL(wp) :: qlwrnet(SIZE(t_acoef_wtr)) !< Net long-wave radiation [W/m**2].
    REAL(wp) :: qsolnet(SIZE(t_acoef_wtr)) !< Net short-wave radiation [W/m**2].
    REAL(wp) :: snow_rate(SIZE(t_acoef_wtr)) !< Snow rate [kg/m**2/s].
    REAL(wp) :: rain_rate(SIZE(t_acoef_wtr)) !< Rain rate [kg/m**2/s].
    REAL(wp) :: tice_p(SIZE(t_acoef_wtr)) !< Ice temperature (time `t`) [K].
    REAL(wp) :: hice_p(SIZE(t_acoef_wtr)) !< Ice thickness (time `t`) [m].
    REAL(wp) :: tsnow_p(SIZE(t_acoef_wtr)) !< Snow temperature on ice (time `t`) [K].
    REAL(wp) :: hsnow_p(SIZE(t_acoef_wtr)) !< Snow thickness on ice (time `t`) [m].
    REAL(wp) :: albsi_p(SIZE(t_acoef_wtr)) !< Sea-ice albedo (time `t`) [1].
    REAL(wp) :: tice_n(SIZE(t_acoef_wtr)) !< Ice temperature (time `t+1`) [K].
    REAL(wp) :: hice_n(SIZE(t_acoef_wtr)) !< Ice thickness (time `t+1`) [m].
    REAL(wp) :: tsnow_n(SIZE(t_acoef_wtr)) !< Snow temperature on ice (time `t+1`) [K].
    REAL(wp) :: hsnow_n(SIZE(t_acoef_wtr)) !< Snow thickness on ice (time `t+1`) [m].
    REAL(wp) :: condhf(SIZE(t_acoef_wtr)) !< Conductive heat flux at ice-ocean boundary [W/m**2].
    REAL(wp) :: meltpot(SIZE(t_acoef_wtr)) !< Melt potential at ice-atmosphere boundary [W/m**2].
    REAL(wp) :: albsi_n(SIZE(t_acoef_wtr)) !< Sea-ice albedo (time `t+1`) [1].

    REAL(wp) :: alb_nir_dir, alb_nir_dif, alb_vis_dir, alb_vis_dif

    LOGICAL :: have_conductive_hflx_ice
    LOGICAL :: have_melt_potential_ice

    ! Asynchronous data regions are a too recent feature. We have to resort to unstructured ones.
    ! This crutch ensures that we don't forget to delete any variable.
#   define LIST_CREATE \
        s_hat_wtr, \
        s_hat_ice, \
        qsat_hat_wtr, \
        qsat_hat_ice, \
        qsen, \
        qlat, \
        qlwrnet, \
        qsolnet, \
        snow_rate, \
        rain_rate, \
        tice_p, \
        hice_p, \
        tsnow_p, \
        hsnow_p, \
        albsi_p, \
        tice_n, \
        hice_n, \
        tsnow_n, \
        hsnow_n, \
        condhf, \
        meltpot, \
        albsi_n
    !$ACC ENTER DATA ASYNC(1) CREATE(LIST_CREATE)

#ifdef __NVCOMPILER
    ! nvfortran does not understand passing a NULL pointer to an optional (Fortran 2008) :(
    have_conductive_hflx_ice = is_coupled_to_ocean()
    have_melt_potential_ice = is_coupled_to_ocean()
#else
    have_conductive_hflx_ice = PRESENT(conductive_hflx_ice)
    have_melt_potential_ice = PRESENT(melt_potential_ice)
#endif

    !$ACC DATA PRESENT(conductive_hflx_ice) IF(have_conductive_hflx_ice)
    !$ACC DATA PRESENT(melt_potential_ice) IF(have_melt_potential_ice)
    !$ACC DATA NO_CREATE(conductive_hflx_ice, melt_potential_ice)

    IF (lseaice) THEN

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR PRIVATE(jc, alb_nir_dir, alb_nir_dif, alb_vis_dir, alb_vis_dif)
        DO ic = 1, ext_data%atm%list_seaice%ncount(iblk)
          jc = ext_data%atm%list_seaice%idx(ic,iblk)

          alb_nir_dir = sea_state%alb_nir_dir(jc,iblk,SFT_SICE)
          alb_nir_dif = sea_state%alb_nir_dif(jc,iblk,SFT_SICE)
          alb_vis_dir = sea_state%alb_vis_dir(jc,iblk,SFT_SICE)
          alb_vis_dif = sea_state%alb_vis_dif(jc,iblk,SFT_SICE)

          qsen(ic) = sensible_hflx_ice_old(jc)
          qlat(ic) = als / alv * latent_hflx_ice_old(jc)
          qlwrnet(ic) = sea_state%lw_emissivity(jc,iblk,SFT_SICE) &
              & * (flx_rad%flx_lw_down(jc,iblk) - stbo * prog_wtr_now%t_ice(jc,iblk)**4)
          qsolnet(ic) = &
              & ((1._wp - alb_nir_dir) &
              &   + (alb_nir_dir - alb_nir_dif) * flx_rad%fr_nir_diffuse(jc,iblk) &
              & ) * flx_rad%flx_nir_down(jc,iblk) &
              & + ((1._wp - alb_vis_dir) &
              &   + (alb_vis_dir - alb_vis_dif) * flx_rad%fr_vis_diffuse(jc,iblk) &
              & ) * flx_rad%flx_vis_down(jc,iblk)
          snow_rate(ic) = snow(jc)
          rain_rate(ic) = rain(jc)
          tice_p(ic) = prog_wtr_now%t_ice(jc,iblk)
          hice_p(ic) = prog_wtr_now%h_ice(jc,iblk)
          tsnow_p(ic) = prog_wtr_now%t_snow_si(jc,iblk)
          hsnow_p(ic) = prog_wtr_now%h_snow_si(jc,iblk)
          albsi_p(ic) = prog_wtr_now%alb_si(jc,iblk)
        END DO
      !$ACC END PARALLEL

      CALL seaice_timestep_nwp ( &
          & dtime=dtime, &
          & nsigb=ext_data%atm%list_seaice%ncount(iblk), &
          & qsen=qsen(:), &
          & qlat=qlat(:), &
          & qlwrnet=qlwrnet(:), &
          & qsolnet=qsolnet(:), &
          & snow_rate=snow_rate(:), &
          & rain_rate=rain_rate(:), &
          & tice_p=tice_p(:), &
          & hice_p=hice_p(:), &
          & tsnow_p=tsnow_p(:), & ! DUMMY: not used yet
          & hsnow_p=hsnow_p(:), & ! DUMMY: not used yet
          & albsi_p=albsi_p(:), &
          & & ! outputs
          & tice_n=tice_n(:), &
          & hice_n=hice_n(:), &
          & tsnow_n=tsnow_n(:), & ! DUMMY: not used yet
          & hsnow_n=hsnow_n(:), & ! DUMMY: not used yet
          & condhf=condhf(:), &
          & meltpot=meltpot(:), &
          & albsi_n=albsi_n(:) &
        )
    END IF

    ! Initialize output with sane values.
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO ic = ics, ice
        t_ice(ic) = tf_fresh
        t_wtr(ic) = tf_salt

        IF (have_conductive_hflx_ice) &
            & conductive_hflx_ice(ic) = 0._wp
        IF (have_melt_potential_ice) &
            & melt_potential_ice(ic) = 0._wp

        prog_wtr_new%t_ice(ic,iblk) = tf_fresh
        prog_wtr_new%t_snow_si(ic,iblk) = tf_fresh
        prog_wtr_new%h_ice(ic,iblk) = 0._wp
        prog_wtr_new%h_snow_si(ic,iblk) = 0._wp
        prog_wtr_new%alb_si(ic,iblk) = csalb(ist_seaice)
      END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(jc)
      DO ic = 1, ext_data%atm%list_sea%ncount(iblk)
        jc = ext_data%atm%list_sea%idx(ic,iblk)

        t_wtr(jc) = MAX(diag_lnd%t_seasfc(jc,iblk), tf_salt)
      END DO
    !$ACC END PARALLEL

    IF (lseaice) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic = 1, ext_data%atm%list_seaice%ncount(iblk)
          jc = ext_data%atm%list_seaice%idx(ic,iblk)

          IF (hice_n(ic) >= hice_min) THEN
            t_ice(jc) = tice_n(ic)

            IF (have_conductive_hflx_ice) &
                & conductive_hflx_ice(jc) = condhf(ic)
            IF (have_melt_potential_ice) &
                & melt_potential_ice(jc) = meltpot(ic)

            prog_wtr_new%t_ice(jc,iblk) = t_ice(jc)
            prog_wtr_new%h_ice(jc,iblk) = hice_n(ic)
            prog_wtr_new%t_snow_si(jc,iblk) = tsnow_n(ic)
            prog_wtr_new%h_snow_si(jc,iblk) = hsnow_n(ic)

            prog_wtr_new%alb_si(jc,iblk) = MIN(MAX(0.001_wp, albsi_n(ic)), 0.999_wp)
          ELSE
            diag_lnd%fr_seaice(jc,iblk) = 0._wp
            ! Keep default sea-ice values.
          END IF

        END DO
      !$ACC END PARALLEL
    END IF

    ! We need a parallel region around these orphaned routines because we are in a parallel region
    ! ourselves.
    !$OMP PARALLEL
      CALL init (alb%alb_vis_dif(ics:ice,iblk,SFT_SWTR), csalb(ist_seawtr), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL init (alb%lw_emissivity(ics:ice,iblk,SFT_SWTR), zemiss_def, lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL init (alb%lw_emissivity(ics:ice,iblk,SFT_SICE), lw_emissivity_ice, lacc=.TRUE., opt_acc_async=.TRUE.)
    !$OMP END PARALLEL

    IF (lprog_albsi) THEN
      !$OMP PARALLEL
        CALL copy ( &
            & prog_wtr_new%alb_si(ics:ice,iblk), &
            & alb%alb_vis_dif(ics:ice,iblk,SFT_SICE), &
            & lacc=.TRUE., &
            & opt_acc_async=.TRUE. &
          )
      !$OMP END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO ic = ics, ice
          alb%alb_vis_dif(ic,iblk,SFT_SICE) = alb_seaice_equil(prog_wtr_new%t_ice(ic,iblk))
        END DO
      !$ACC END PARALLEL
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = ics, ice
        SELECT CASE (direct_albedo_water)
        CASE (1) ! Ritter-Geleyn
          alb%alb_vis_dir(ic,iblk,SFT_SWTR) = &
              & sfc_albedo_dir_rg (cos_zenith_angle(ic), alb%alb_vis_dif(ic,iblk,SFT_SWTR))
        CASE (2,4) ! (2) Yang, (4) Yang over sea water, RG over lakes.
          alb%alb_vis_dir(ic,iblk,SFT_SWTR) = &
              & sfc_albedo_dir_yang (cos_zenith_angle(ic), alb%alb_vis_dif(ic,iblk,SFT_SWTR))
        CASE (3) ! Taylor
          alb%alb_vis_dir(ic,iblk,SFT_SWTR) = &
              & sfc_albedo_dir_taylor (cos_zenith_angle(ic))
        END SELECT

        alb%alb_vis_dir(ic,iblk,SFT_SICE) = &
            & sfc_albedo_dir_rg (cos_zenith_angle(ic), alb%alb_vis_dif(ic,iblk,SFT_SICE))
      END DO

      IF (albedo_whitecap == 1) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wc_frac, wc_alb)
        DO ic = ics, ice
          CALL sfc_albedo_whitecap (wind_10m(ic), wc_frac, wc_alb)
          alb%alb_vis_dif(ic,iblk,SFT_SWTR) = &
              & (1._wp - wc_frac) * alb%alb_vis_dif(ic,iblk,SFT_SWTR) + wc_frac * wc_alb
        END DO
      END IF
    !$ACC END PARALLEL

    !$OMP PARALLEL
      ! No distinction between visible and near-IR albedo.
      CALL copy ( &
          & alb%alb_vis_dir(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & alb%alb_nir_dir(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )
      CALL copy ( &
          & alb%alb_vis_dif(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & alb%alb_nir_dif(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )

      ! Save state for next invocation.
      CALL copy ( &
          & alb%alb_nir_dir(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & sea_state%alb_nir_dir(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )
      CALL copy ( &
          & alb%alb_nir_dif(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & sea_state%alb_nir_dif(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )
      CALL copy ( &
          & alb%alb_vis_dir(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & sea_state%alb_vis_dir(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )
      CALL copy ( &
          & alb%alb_vis_dif(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & sea_state%alb_vis_dif(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )
      CALL copy ( &
          & alb%lw_emissivity(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & sea_state%lw_emissivity(ics:ice,iblk,SFT_SWTR:SFT_SICE), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )
    !$OMP END PARALLEL

    ! Compute saturation specific humidity and dry static energy over the sea surface. This is
    ! used to get the heat and humidity fluxes into the lowest grid box.
    ! Note: These should all be mixed-time values to be exact, but we do not have access to the
    ! old pressure. The sea surface temperature changes infrequently, so it can be considered
    ! constant.
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(t_ice_old, qsat_ice_now)
      DO ic = ics, ice

        ! Someone is messing with the old ice temperatures, setting them to zero for ice-free
        ! cells.
        t_ice_old = MERGE( &
            & tf_fresh, prog_wtr_now%t_ice(ic,iblk), prog_wtr_now%t_ice(ic,iblk) < 100._wp)

        qsat_ice_now = spec_humi(sat_pres_ice(t_ice_old),  press_srf(ic))

        qsat_wtr(ic) = salinity_fac * spec_humi(sat_pres_water(t_wtr(ic)), press_srf(ic))
        qsat_ice(ic) = spec_humi(sat_pres_ice(t_ice(ic)),  press_srf(ic))

        qsat_hat_wtr(ic) = qsat_wtr(ic)
        qsat_hat_ice(ic) = vdiff_mixed_time_value(qsat_ice(ic), qsat_ice_now)

        s_wtr(ic) = cpd * t_wtr(ic)
        s_ice(ic) = cpd * t_ice(ic)
        s_hat_wtr(ic) = s_wtr(ic)
        s_hat_ice(ic) = vdiff_mixed_time_value(s_ice(ic), cpd * t_ice_old)
      END DO
    !$ACC END PARALLEL

    ! Compute fluxes into the surface resulting from the new sea (ice) surface temperatures.
    ! Runs on async queue 1.
    CALL vdiff_surface_flux( &
        & ics, ice, steplen, prefactor_exchange(:), exchange_coeff_h_wtr(:), t_acoef_wtr(:), &
        & t_bcoef_wtr(:), s_hat_wtr(:), sensible_hflx_wtr(:) &
      )
    CALL vdiff_surface_flux( &
        & ics, ice, steplen, prefactor_exchange(:), exchange_coeff_h_wtr(:), q_acoef_wtr(:), &
        & q_bcoef_wtr(:), qsat_hat_wtr(:), evapo_wtr(:) &
      )
    CALL vdiff_surface_flux( &
        & ics, ice, steplen, prefactor_exchange(:), exchange_coeff_h_ice(:), t_acoef_ice(:), &
        & t_bcoef_ice(:), s_hat_ice(:), sensible_hflx_ice(:) &
      )
    CALL vdiff_surface_flux( &
        & ics, ice, steplen, prefactor_exchange(:), exchange_coeff_h_ice(:), q_acoef_ice(:), &
        & q_bcoef_ice(:), qsat_hat_ice(:), evapo_ice(:) &
      )

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO ic = ics, ice
        latent_hflx_wtr(ic) = alv * evapo_wtr(ic)
        latent_hflx_ice(ic) = alv * evapo_ice(ic)
      END DO
    !$ACC END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA ! NO_CREATE(conductive_hflx_ice, melt_potential_ice)
    !$ACC END DATA ! PRESENT(melt_potential_ice)
    !$ACC END DATA ! PRESENT(conductive_hflx_ice)
    !$ACC EXIT DATA DELETE(LIST_CREATE)
#   undef LIST_CREATE

  END SUBROUTINE sea_model


  !>
  !! Initialize the sea "model". This sets up the current albedos, which are needed to obtain net
  !! from downward fluxes in the first time step.
  SUBROUTINE sea_model_init (patch, ext_data, cos_zenith_angle, t_seasfc, sst_m, sea_state)

    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.
    TYPE(t_external_data), INTENT(IN) :: ext_data !< Patch external data.
    REAL(wp), INTENT(IN) :: cos_zenith_angle(:,:) !< Cosine of the zenith angle [1].
    REAL(wp), INTENT(IN) :: t_seasfc(:,:) !< Sea-surface temperature [K].
    REAL(wp), INTENT(IN) :: sst_m(:,:,:) !< Monthly mean SSTs [K].
    TYPE(t_nwp_vdiff_sea_state), INTENT(INOUT) :: sea_state !< sea model state structure.

    !> Climatological SST for experiment start date [K].
    REAL(wp) :: t_seasfc_clim(nproma, patch%nblks_c)

    INTEGER :: ic, ics, ice
    INTEGER :: i_blk, i_startblk, i_endblk

    IF (isRestart()) RETURN

    sea_state%alb_vis_dif(:,:,SFT_SWTR) = csalb(ist_seawtr)
    sea_state%alb_vis_dif(:,:,SFT_SICE) = csalb(ist_seaice)
    sea_state%alb_vis_dir(:,:,SFT_SWTR) = csalb(ist_seawtr)
    sea_state%alb_vis_dir(:,:,SFT_SICE) = csalb(ist_seaice)
    sea_state%lw_emissivity(:,:,SFT_SWTR) = zemiss_def
    sea_state%lw_emissivity(:,:,SFT_SICE) = lw_emissivity_ice

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    DO i_blk = i_startblk, i_endblk
      CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
          & end_prog_cells)
      DO ic = ics, ice
        SELECT CASE (direct_albedo_water)
        CASE (1) ! Ritter-Geleyn
          sea_state%alb_vis_dir(ic,i_blk,SFT_SWTR) = &
              & sfc_albedo_dir_rg( &
              &   cos_zenith_angle(ic,i_blk), sea_state%alb_vis_dif(ic,i_blk,SFT_SWTR))
        CASE (2,4) ! (2) Yang, (4) Yang over sea water, RG over lakes.
          sea_state%alb_vis_dir(ic,i_blk,SFT_SWTR) = &
              & sfc_albedo_dir_yang( &
              &   cos_zenith_angle(ic,i_blk), sea_state%alb_vis_dif(ic,i_blk,SFT_SWTR))
        CASE (3) ! Taylor
          sea_state%alb_vis_dir(ic,i_blk,SFT_SWTR) = &
              & sfc_albedo_dir_taylor(cos_zenith_angle(ic,i_blk))
        END SELECT

        sea_state%alb_vis_dir(ic,i_blk,SFT_SICE) = &
            & sfc_albedo_dir_rg( &
            &   cos_zenith_angle(ic,i_blk), sea_state%alb_vis_dif(ic,i_blk,SFT_SICE))
      END DO
    END DO

    sea_state%alb_nir_dif(:,:,SFT_SWTR) = sea_state%alb_vis_dif(:,:,SFT_SWTR)
    sea_state%alb_nir_dif(:,:,SFT_SICE) = sea_state%alb_vis_dif(:,:,SFT_SICE)
    sea_state%alb_nir_dir(:,:,SFT_SWTR) = sea_state%alb_vis_dir(:,:,SFT_SWTR)
    sea_state%alb_nir_dir(:,:,SFT_SICE) = sea_state%alb_vis_dir(:,:,SFT_SICE)

    IF (ASSOCIATED(sea_state%flx_co2_natural_sea)) THEN
      sea_state%flx_co2_natural_sea(:,:) = 0._wp
    END IF
    IF (ASSOCIATED(sea_state%ocean_u)) THEN
      sea_state%ocean_u(:,:) = 0._wp
      sea_state%ocean_v(:,:) = 0._wp
    END IF

    IF (sstice_mode == SSTICE_ANA_CLINC) THEN
      CALL interpol_monthly_mean( &
          & patch, sea_state%time_ref_t_seasfc, sst_m, t_seasfc_clim(:,:))

      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)

        sea_state%t_seasfc_offset(:,i_blk) = 0._wp
        ASSOCIATE (idx => ext_data%atm%list_sea%idx(1:ext_data%atm%list_sea%ncount(i_blk), i_blk))
          sea_state%t_seasfc_offset(idx,i_blk) = &
              & t_seasfc(idx,i_blk) - t_seasfc_clim(idx,i_blk)
        END ASSOCIATE
      END DO
    END IF

    CALL sea_state%h2d()

  END SUBROUTINE sea_model_init


  !> Update sea-surface temperature for the analysis + cimatological increments mode. For the
  !! other modes t_seasfc is kept current in mo_nwp_sfc_utils:update_sst_and_seaice.
  SUBROUTINE sea_model_update_sst (patch, sea_state, current_datetime, sst_m, t_seasfc)

    !> Current patch.
    TYPE(t_patch), INTENT(IN) :: patch
    !> Current state of the sea model.
    TYPE(t_nwp_vdiff_sea_state), INTENT(INOUT) :: sea_state
    !> Current time.
    TYPE(datetime), INTENT(IN) :: current_datetime
    !> Monthly SST means [K].
    REAL(wp), INTENT(IN) :: sst_m(:,:,:)
    !> Sea-surface temperature [K].
    REAL(wp), INTENT(INOUT) :: t_seasfc(:,:)

    INTEGER :: i_startblk, i_endblk, i_blk
    INTEGER :: ics, ice, ic

    IF (sstice_mode /= SSTICE_ANA_CLINC) RETURN
    IF (current_datetime%date%day == sea_state%time_last_update_t_seasfc%date%day) RETURN

    CALL interpol_monthly_mean (patch, current_datetime, sst_m, t_seasfc)

    !$ACC UPDATE DEVICE(t_seasfc) ASYNC(1)

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    !$OMP PARALLEL
      !$OMP DO PRIVATE(i_blk, ics, ice, ic)
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO ic = ics, ice
            t_seasfc(ic,i_blk) = t_seasfc(ic,i_blk) + sea_state%t_seasfc_offset(ic,i_blk)
          END DO
        !$ACC END PARALLEL
      END DO
    !$OMP END PARALLEL

    sea_state%time_last_update_t_seasfc = current_datetime

  END SUBROUTINE sea_model_update_sst


  SUBROUTINE sea_model_couple_ocean ( &
        & patch, list_sea, fr_sft, alb, flx_rad, pres_sfc, t_eff_sft, evapo_sft, flx_heat_latent_sft, &
        & flx_heat_sensible_sft, condhf_ice, meltpot_ice, co2_concentration_srf, rain_srf, snow_srf, &
        & umfl_sft, vmfl_sft, sp_10m, t_seasfc, fr_seaice, h_ice, sea_state &
      )

    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_sea !< Index list of sea points.
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: fr_sft(:,:,:) !< Fractions of each surface type.
    TYPE(t_nwp_vdiff_albedos), INTENT(IN) :: alb !< Albedo package.
    TYPE(t_nwp_vdiff_surface_rad_fluxes), INTENT(IN) :: flx_rad !< Surface flux package.

    !> Surface pressure [Pa] (nproma, nblks_c).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: pres_sfc(:,:)

    !> Effective temperature for radiation per surface type [K] (nproma, nblks_c, SFT_NUM).
    REAL(wp), INTENT(IN) :: t_eff_sft(:,:,:)

    !> Surface evaporation over surface types [kg/m**2/s] (nproma, nblks_c, SFT_NUM).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: evapo_sft(:,:,:)
    !> Latent surface flux over surface types [W/m**2] (nproma, nblks_c, SFT_NUM).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: flx_heat_latent_sft(:,:,:)
    !> Sensible surface flux over surface types [W/m**2] (nproma, nblks_c, SFT_NUM).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: flx_heat_sensible_sft(:,:,:)
    !> Conductive heat flux at ice bottom [W/m**2] (nproma, nblks_c).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: condhf_ice(:,:)
    !> Melt potential at ice top [W/m**2] (nproma, nblks_c).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: meltpot_ice(:,:)

    !> Surface CO2 concentration [kg(CO2)/kg(air)] (nproma, nblks_c).
    REAL(wp), TARGET, INTENT(IN) :: co2_concentration_srf(:,:)

    !> Surface rain [kg/m**2/s] (nproma, nblks_c).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: rain_srf(:,:)
    !> Surface snow [kg/m**2/s] (nproma, nblks_c).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: snow_srf(:,:)

    !> Zonal surface stress [N/m**2] (nproma, nblks_c, SFT_NUM).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: umfl_sft(:,:,:)
    !> Meridional surface stress [N/m**2] (nproma, nblks_c, SFT_NUM).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: vmfl_sft(:,:,:)
    !> 10m wind speed [m/s] (nproma, nblks_c).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: sp_10m(:,:)

    !> Sea surface temperature [K] (nproma, nblks_c).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: t_seasfc(:,:)
    !> Sea-ice fraction [m**2(ice)/m**2(ocean)] (nproma, nblks_c).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: fr_seaice(:,:)
    !> Sea-ice height [m] (nproma, nblks_c).
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: h_ice(:,:)
    !> Sea state.
    TYPE(t_nwp_vdiff_sea_state), INTENT(INOUT) :: sea_state

    ! Locals

    TYPE(t_nwp_ocean_fields_rx) :: rx
    TYPE(t_nwp_ocean_fields_tx) :: tx

    REAL(wp), TARGET :: lwflx(nproma, patch%nblks_c, SFT_SWTR:SFT_SICE) !< Net LW flux [W/m**2].
    REAL(wp), TARGET :: swflx(nproma, patch%nblks_c, SFT_SWTR:SFT_SICE) !< Net SW flux [W/m**2].
    REAL(wp), TARGET :: lhfl_s_i(nproma, patch%nblks_c) !< Latent heat flux over ice [W/m**2].

    REAL(wp) :: alb_nir_dir, alb_nir_dif, alb_vis_dir, alb_vis_dif, lw_emissivity

    INTEGER :: i_blk, i_startblk, i_endblk
    INTEGER :: ic, ics, ice
    INTEGER :: isft
    INTEGER :: jc

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    IF (ltimer) CALL timer_start(timer_coupling)

    !$OMP PARALLEL
      !$OMP DO PRIVATE(i_blk, ics, ice, ic, isft) &
      !$OMP   PRIVATE(alb_nir_dir, alb_nir_dif, alb_vis_dir, alb_vis_dif, lw_emissivity)
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)
        DO ic = ics, ice
          DO isft = SFT_SWTR, SFT_SICE
            alb_nir_dir = alb%alb_nir_dir(ic,i_blk,isft)
            alb_nir_dif = alb%alb_nir_dif(ic,i_blk,isft)
            alb_vis_dir = alb%alb_vis_dir(ic,i_blk,isft)
            alb_vis_dif = alb%alb_vis_dif(ic,i_blk,isft)
            lw_emissivity = alb%lw_emissivity(ic,i_blk,isft)

            lwflx(ic,i_blk,isft) = lw_emissivity &
                & * (flx_rad%flx_lw_down(ic,i_blk) - stbo * t_eff_sft(ic,i_blk,isft)**4)
            swflx(ic,i_blk,isft) = &
                & ((1._wp - alb_nir_dir) &
                &   + (alb_nir_dir - alb_nir_dif) * flx_rad%fr_nir_diffuse(ic,i_blk) &
                & ) * flx_rad%flx_nir_down(ic,i_blk) &
                & + ((1._wp - alb_vis_dir) &
                &   + (alb_vis_dir - alb_vis_dif) * flx_rad%fr_vis_diffuse(ic,i_blk) &
                & ) * flx_rad%flx_vis_down(ic,i_blk)
          END DO

          ! The amount of energy transferred to/from the ice is slightly higher than what ends up in
          ! the atmosphere because ice has negative latent heat.
          lhfl_s_i(ic,i_blk) = als / alv * flx_heat_latent_sft(ic,i_blk,SFT_SICE)
        END DO
      END DO
    !$OMP END PARALLEL

    tx%chfl_i => condhf_ice(:,:)
    tx%meltpot_i => meltpot_ice(:,:)
    tx%frac_w => fr_sft(:,:,SFT_SWTR)
    tx%frac_i => fr_sft(:,:,SFT_SICE)
    tx%lhfl_s_i => lhfl_s_i(:,:)
    tx%lhfl_s_w => flx_heat_latent_sft(:,:,SFT_SWTR)
    tx%lwflxsfc_i => lwflx(:,:,SFT_SICE)
    tx%lwflxsfc_w => lwflx(:,:,SFT_SWTR)
    tx%pres_sfc => pres_sfc(:,:)
    tx%q_co2 => co2_concentration_srf(:,:)
    tx%qhfl_s_i => evapo_sft(:,:,SFT_SICE)
    tx%qhfl_s_w => evapo_sft(:,:,SFT_SWTR)
    tx%rain_rate => rain_srf(:,:)
    tx%shfl_s_i => flx_heat_sensible_sft(:,:,SFT_SICE)
    tx%shfl_s_w => flx_heat_sensible_sft(:,:,SFT_SWTR)
    tx%snow_rate => snow_srf(:,:)
    tx%swflxsfc_i => swflx(:,:,SFT_SICE)
    tx%swflxsfc_w => swflx(:,:,SFT_SWTR)
    tx%umfl_s_i => umfl_sft(:,:,SFT_SICE)
    tx%umfl_s_w => umfl_sft(:,:,SFT_SWTR)
    tx%vmfl_s_i => vmfl_sft(:,:,SFT_SICE)
    tx%vmfl_s_w => vmfl_sft(:,:,SFT_SWTR)
    tx%sp_10m => sp_10m(:,:) ! Use grid-box average.

    rx%fr_seaice => fr_seaice(:,:)
    rx%h_ice => h_ice(:,:)
    rx%ocean_u => sea_state%ocean_u(:,:)
    rx%ocean_v => sea_state%ocean_v(:,:)
    rx%t_seasfc => t_seasfc(:,:)

    ! May not be associated.
    rx%flx_co2 => sea_state%flx_co2_natural_sea(:,:)

    CALL couple_ocean(patch, list_sea, tx, rx)

    ! Limit sea-ice height to allowed range.
    !$OMP PARALLEL
      !$OMP DO PRIVATE(i_blk, ic, jc)
      DO i_blk = i_startblk, i_endblk
        DO ic = 1, list_sea%ncount(i_blk)
          jc = list_sea%idx(ic,i_blk)

          IF (fr_seaice(jc,i_blk) >= frsi_min) THEN
            h_ice(jc,i_blk) = MIN(MAX(hice_min, h_ice(jc,i_blk)), hice_max)
          END IF
        END DO
      END DO
    !$OMP END PARALLEL

    IF (ltimer) CALL timer_stop(timer_coupling)

  END SUBROUTINE sea_model_couple_ocean
  !>
  !! Update the sea-ice fraction with new external data.
  SUBROUTINE nwp_vdiff_update_seaice ( &
        & patch, init_hice, fr_seaice, sea_list, seaice_list, prog_wtr, lacc &
      )

    !> Current patch.
    TYPE(t_patch), INTENT(IN) :: patch
    !> If .TRUE. initialize ice height to default initial value. Else leave height unchanged.
    LOGICAL :: init_hice
    !> New sea-ice fraction.
    REAL(wp), INTENT(IN) :: fr_seaice(:,:)
    !> Sea list.
    TYPE(t_idx_list_blocked), INTENT(IN) :: sea_list
    !> Sea-ice list.
    TYPE(t_idx_list_blocked), INTENT(INOUT) :: seaice_list
    !> Prognostic water variables.
    TYPE(t_wtr_prog), INTENT(INOUT) :: prog_wtr
    !> Use OpenACC.
    LOGICAL, OPTIONAL, INTENT(IN) :: lacc

    TYPE(t_idx_list_blocked) :: new_ice_list
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    CALL nwp_vdiff_update_seaice_list ( &
        & patch, fr_seaice(:,:), sea_list, seaice_list, new_ice_list=new_ice_list, &
        & prog_wtr=prog_wtr, lacc=lzacc &
      )

    CALL nwp_vdiff_update_seaice_vars ( &
        & patch, init_hice, fr_seaice(:,:), new_ice_list, prog_wtr, lacc=lzacc &
      )

  END SUBROUTINE nwp_vdiff_update_seaice

  !>
  !! Update dynamic sea-ice list with all cells with a sea-ice fraction of more than `frsi_min`.
  !! Optionally resets `h_ice` on newly ice-free points, and returns a list of new ice points.
  !!
  !! If `new_ice_list` is present, the sea-ice indices have to be sorted ascending. After the
  !! routine finishes, both lists are present on CPU and GPU.
  SUBROUTINE nwp_vdiff_update_seaice_list ( &
        & patch, fr_seaice, sea_list, seaice_list, new_ice_list, prog_wtr, lacc &
      )

    !> Current patch.
    TYPE(t_patch), INTENT(IN) :: patch
    !> Sea-ice fraction.
    REAL(wp), INTENT(IN) :: fr_seaice(:,:)
    !> Sea list.
    TYPE(t_idx_list_blocked), INTENT(IN) :: sea_list
    !> Sea-ice list.
    TYPE(t_idx_list_blocked), INTENT(INOUT) :: seaice_list
    !> List of new sea-ice points.
    TYPE(t_idx_list_blocked), OPTIONAL, INTENT(OUT) :: new_ice_list
    !> Prognostic water variables. Resets h_ice on ice-free points if present.
    TYPE(t_wtr_prog), OPTIONAL, INTENT(INOUT) :: prog_wtr
    !> Use OpenACC.
    LOGICAL, OPTIONAL, INTENT(IN) :: lacc

    INTEGER :: i_startblk, i_endblk
    INTEGER :: ic, jc, iblk, i_count, io
    INTEGER :: oldice_idx(nproma)
    LOGICAL :: newice_found
    LOGICAL :: lzacc

    INTEGER :: g_count_ice_prev
    INTEGER :: g_count_ice_now
    INTEGER :: g_count_ice_created

    LOGICAL :: have_new_ice_list
    LOGICAL :: have_prog_wtr

    INTEGER, PARAMETER :: MSG_LEVEL_DIAG = 12

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC ENTER DATA ASYNC(1) CREATE(oldice_idx) IF(lzacc)

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    have_new_ice_list = PRESENT(new_ice_list)
    have_prog_wtr = PRESENT(prog_wtr)

    g_count_ice_created = 0

    IF (msg_level >= MSG_LEVEL_DIAG) THEN
      g_count_ice_prev = seaice_list%get_sum_global(i_startblk, i_endblk)
    END IF

    IF (have_new_ice_list) THEN
      CALL new_ice_list%construct(nproma, patch%nblks_c, lopenacc=lzacc)
    END IF

    !$OMP PARALLEL
      !$OMP DO PRIVATE(iblk, ic, jc, io, i_count, oldice_idx, newice_found) &
      !$OMP   REDUCTION(+: g_count_ice_created)
      DO iblk = i_startblk, i_endblk
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          IF (have_prog_wtr) THEN
            !NEC$ ivdep
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc)
            DO ic = 1, seaice_list%ncount(iblk)
              jc = seaice_list%idx(ic, iblk)

              IF (fr_seaice(jc,iblk) < frsi_min) THEN
                prog_wtr%h_ice(jc,iblk) = 0._wp
              END IF
            END DO
          END IF

          IF (have_new_ice_list .OR. msg_level >= MSG_LEVEL_DIAG) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO ic = 1, seaice_list%ncount(iblk)
              oldice_idx(ic) = seaice_list%idx(ic,iblk)
            END DO
          END IF
        !$ACC END PARALLEL

        !$ACC SERIAL DEFAULT(PRESENT) PRIVATE(i_count, ic, io, jc, newice_found) REDUCTION(+: g_count_ice_created) ASYNC(1) IF(lzacc)
          i_count = seaice_list%ncount(iblk)
          seaice_list%ncount(iblk) = 0

          DO ic = 1, sea_list%ncount(iblk)
            jc = sea_list%idx(ic,iblk)

            IF (fr_seaice(jc,iblk) >= frsi_min) THEN
              seaice_list%ncount(iblk) = seaice_list%ncount(iblk) + 1
              seaice_list%idx(seaice_list%ncount(iblk),iblk) = jc
            END IF
          END DO

          IF (have_new_ice_list .OR. msg_level >= MSG_LEVEL_DIAG) THEN
            ! Find new ice points. Assumes that the old ice list is sorted (the new one is sorted by
            ! construction).

            io = 1
            DO ic = 1, seaice_list%ncount(iblk)
              jc = seaice_list%idx(ic,iblk)

              ! That's what you get when your language does not guarantee shortcut evaluation.
              DO WHILE (io <= i_count)
                IF (oldice_idx(io) >= jc) EXIT
                io = io + 1
              END DO

              IF (io > i_count) THEN
                ! Reached the end of the old ice list. Everything beyond must be new ice.
                newice_found = .TRUE.
              ELSE IF (oldice_idx(io) > jc) THEN
                ! The next index in the old ice list is larger. Current point cannot be in list.
                newice_found = .TRUE.
              ELSE
                ! oldice_idx(io) == jc, present in list.
                newice_found = .FALSE.
              END IF

              IF (newice_found) THEN
                g_count_ice_created = g_count_ice_created + 1
              END IF

              IF (newice_found .AND. have_new_ice_list) THEN
                new_ice_list%ncount(iblk) = new_ice_list%ncount(iblk) + 1
                new_ice_list%idx(new_ice_list%ncount(iblk),iblk) = jc
              END IF
            END DO
          END IF
        !$ACC END SERIAL
      END DO
    !$OMP END PARALLEL

    IF (msg_level >= MSG_LEVEL_DIAG) THEN
      g_count_ice_created = global_sum(g_count_ice_created)
      g_count_ice_now = seaice_list%get_sum_global(i_startblk, i_endblk)

      IF (g_count_ice_created == 0 .AND. g_count_ice_now == g_count_ice_prev) THEN
        WRITE (message_text, '(A,I8,A,I8)') &
            & 'Num cells: ', g_count_ice_prev, ' -> ', g_count_ice_now
      ELSE
        WRITE (message_text, '(A,I8,A,I8,A,A,I8,A,A,I8)') &
            & 'Num cells: ', g_count_ice_prev, ' -> ', g_count_ice_now, NEW_LINE('a'), &
            & '   Cells deleted: ', &
            &   g_count_ice_prev - (g_count_ice_now - g_count_ice_created), NEW_LINE('a'), &
            & '   Cells created: ', g_count_ice_created
      END IF

      CALL message ('nwp_vdiff_update_seaice_list', message_text)
    END IF

    ! Copy index lists to CPU. Depending on context, the CPU copy may be used.
    !$ACC UPDATE HOST(seaice_list%idx, seaice_list%ncount) ASYNC(1) IF(lzacc)
    !$ACC UPDATE HOST(new_ice_list%idx, new_ice_list%ncount) ASYNC(1) &
    !$ACC   IF(lzacc .AND. have_new_ice_list)

    !$ACC WAIT(1) IF(lzacc)
    !$ACC EXIT DATA DELETE(oldice_idx) IF(lzacc)

  END SUBROUTINE nwp_vdiff_update_seaice_list


  !>
  !! Initialize new ice points.
  SUBROUTINE nwp_vdiff_update_seaice_vars ( &
        & patch, init_hice, fr_seaice, new_ice_list, prog_wtr, lacc &
      )

    !> Current patch.
    TYPE(t_patch), INTENT(IN) :: patch
    !> If .TRUE. initialize ice height to default initial value. Else leave height unchanged.
    LOGICAL :: init_hice
    !> Sea-ice fraction.
    REAL(wp), INTENT(IN) :: fr_seaice(:,:)
    !> List of new ice points.
    TYPE(t_idx_list_blocked), INTENT(IN) :: new_ice_list
    !> Prognostic water variables.
    TYPE(t_wtr_prog), INTENT(INOUT) :: prog_wtr
    !> Use OpenACC.
    LOGICAL, OPTIONAL, INTENT(IN) :: lacc

    REAL(wp) :: frsi(nproma)

    REAL(wp) :: tice_p(nproma)
    REAL(wp) :: hice_p(nproma)
    REAL(wp) :: tsnow_p(nproma)
    REAL(wp) :: hsnow_p(nproma)
    REAL(wp) :: albsi_p(nproma)

    REAL(wp) :: tice_n(nproma)
    REAL(wp) :: hice_n(nproma)
    REAL(wp) :: tsnow_n(nproma)
    REAL(wp) :: hsnow_n(nproma)
    REAL(wp) :: albsi_n(nproma)

    INTEGER :: iblk, i_count
    INTEGER :: ic, jc

    LOGICAL :: lzacc

    IF (.NOT. lseaice) RETURN

    CALL set_acc_host_or_device(lzacc, lacc)

#   define LIST_CREATE \
        frsi, \
        tice_p, \
        hice_p, \
        tsnow_p, \
        hsnow_p, \
        albsi_p, \
        tice_n, \
        hice_n, \
        tsnow_n, \
        hsnow_n, \
        albsi_n
    !$ACC ENTER DATA ASYNC(1) CREATE(LIST_CREATE) IF(lzacc)

    !$OMP PARALLEL
      !$OMP DO PRIVATE(iblk, i_count, ic, jc, frsi, tice_p, hice_p, tsnow_p, hsnow_p, albsi_p) &
      !$OMP   PRIVATE(tice_n, hice_n, tsnow_n, hsnow_n, albsi_n)
      DO iblk = 1, patch%nblks_c
        i_count = new_ice_list%ncount(iblk)

        IF (i_count == 0) CYCLE

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic = 1, new_ice_list%ncount(iblk)
          jc = new_ice_list%idx(ic,iblk)

          frsi(ic) = fr_seaice(jc, iblk)
          tice_p(ic) = prog_wtr%t_ice(jc, iblk)
          hice_p(ic) = prog_wtr%h_ice(jc, iblk)
          tsnow_p(ic) = prog_wtr%t_snow_si(jc, iblk)
          hsnow_p(ic) = prog_wtr%h_snow_si(jc, iblk)
          albsi_p(ic) = prog_wtr%alb_si(jc, iblk)
        END DO
        !$ACC END PARALLEL

        !$ACC UPDATE HOST(frsi, tice_p, hice_p, tsnow_p, hsnow_p, albsi_p) ASYNC(1) IF(lzacc)
        !$ACC WAIT(1)

        tice_n(1:i_count) = tice_p(1:i_count)
        hice_n(1:i_count) = hice_p(1:i_count)
        tsnow_n(1:i_count) = tsnow_p(1:i_count)
        hsnow_n(1:i_count) = hsnow_p(1:i_count)
        albsi_n(1:i_count) = albsi_p(1:i_count)

        CALL seaice_init_nwp ( &
            & init_hice, i_count, frsi(:), &
            & tice_p(:), hice_p(:), tsnow_p(:), hsnow_p(:), albsi_p(:), &
            & tice_n(:), hice_n(:), tsnow_n(:), hsnow_n(:), albsi_n(:) &
          )

        !$ACC UPDATE DEVICE(tice_p, hice_p, tsnow_p, hsnow_p, albsi_p) ASYNC(1) IF(lzacc)

        !NEC$ ivdep
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR PRIVATE(jc)
          DO ic = 1, new_ice_list%ncount(iblk)
            jc = new_ice_list%idx(ic,iblk)

            prog_wtr%t_ice(jc, iblk) = tice_p(ic)
            prog_wtr%h_ice(jc, iblk) = hice_p(ic)
            prog_wtr%t_snow_si(jc, iblk) = tsnow_p(ic)
            prog_wtr%h_snow_si(jc, iblk) = hsnow_p(ic)
            prog_wtr%alb_si(jc, iblk) = albsi_p(ic)
          END DO
        !$ACC END PARALLEL
      END DO
    !$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(LIST_CREATE) IF(lzacc)

  END SUBROUTINE nwp_vdiff_update_seaice_vars

END MODULE mo_nwp_vdiff_sea
