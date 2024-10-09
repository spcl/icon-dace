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

! Interface to the VDIFF turbulence scheme and JSBACH land-surface scheme.

MODULE mo_nwp_vdiff_interface

  USE mo_bc_greenhouse_gases, ONLY: ghg_co2mmr
  USE mo_ccycle_config, ONLY: &
      & CCYCLE_MODE_NONE, CCYCLE_MODE_INTERACTIVE, CCYCLE_MODE_PRESCRIBED, CCYCLE_CO2CONC_CONST, &
      & CCYCLE_CO2CONC_FROMFILE, t_ccycle_config
  USE mo_convect_tables, ONLY: init_convect_tables
  USE mo_coupling_config, ONLY: is_coupled_to_ocean
  USE mo_aes_convect_tables, ONLY: init_aes_convect_tables => init_convect_tables
  USE mo_exception, ONLY: finish, message
  USE mo_ext_data_types, ONLY: t_external_data
  USE mo_fortran_tools, ONLY: assert_lacc_equals_i_am_accel_node, &
      & assert_acc_device_only, copy, init, if_associated
  USE mo_mpi, ONLY: i_am_accel_node
  USE mo_impl_constants, ONLY: end_prog_cells, start_prog_cells
  USE mo_kind, ONLY: wp
  USE mo_lnd_nwp_config, ONLY: frsea_thrhld
  USE mo_loopindices, ONLY: get_indices_c
  USE mo_master_config, ONLY: isRestart
  USE mo_model_domain, ONLY: t_patch
  USE mo_nonhydro_types, ONLY: t_nh_diag, t_nh_metrics, t_nh_prog
  USE mo_nwp_lnd_types, ONLY: t_lnd_diag, t_lnd_prog, t_wtr_prog
  USE mo_nwp_phy_types, ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_tuning_config, ONLY: tune_gust_factor
  USE mo_nwp_vdiff_radfluxes, ONLY: t_nwp_vdiff_surface_rad_fluxes
  USE mo_nwp_vdiff_sea, ONLY: sea_model, sea_model_init, sea_model_couple_ocean, &
      & sea_model_update_sst, nwp_vdiff_update_seaice, nwp_vdiff_update_seaice_list
  USE mo_nwp_vdiff_types, ONLY: &
      & SFT_LAND, SFT_LWTR, SFT_LICE, SFT_SWTR, SFT_SICE, SFT_NUM, SFT_CLASS, &
      & SFC_LAND, SFC_WATER, SFC_ICE, SFC_NUM, &
      & t_nwp_vdiff_albedos, t_nwp_vdiff_state
  USE mo_orbit, ONLY: orbit_vsop87
  USE mo_parallel_config, ONLY: nproma
  USE mo_physical_constants, ONLY: cpd, cvd, cvv, grav, rd, rdv, tf_fresh => tmelt, vtmpc1, &
      & zemiss_def, vmr_to_mmr_co2
  USE mo_run_config, ONLY: ico2, iqc, iqi, nqtendphy, iqt, iqv, ntracer
  USE mo_satad, ONLY: latent_heat_vaporization, sat_pres_water, spec_humi
  USE mo_turb_vdiff, ONLY: &
      & imh_vdiff => imh, imqv_vdiff => imqv, ih_vdiff => ih, iqc_vdiff => ixl, &
      & iqv_vdiff => iqv, nvar_vdiff, matrix_to_richtmyer_coeff, &
      & nmatrix_vdiff => nmatrix, vdiff_down, vdiff_init, &
      & vdiff_new_time_value, vdiff_get_richtmyer_coeff_momentum, vdiff_surface_flux, vdiff_up, &
      & vdiff_update_boundary, vdiff_get_tke
  USE mo_turb_vdiff_config, ONLY: t_vdiff_config
  USE mo_turb_vdiff_params, ONLY: vdiff_implfact => cvdifts, VDIFF_TURB_3DSMAGORINSKY
  USE mtime, ONLY: datetime, julianday, getJulianDayFromDatetime

#ifndef __NO_JSBACH__
  USE mo_jsb_interface, ONLY: &
      & jsbach_interface, jsbach_finish_timestep, jsbach_start_timestep, jsbach_get_var
  USE mo_jsb_model_init, ONLY: jsbach_init
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: nwp_vdiff
  PUBLIC :: nwp_vdiff_init
  PUBLIC :: nwp_vdiff_setup
  PUBLIC :: nwp_vdiff_update_seaice
  PUBLIC :: nwp_vdiff_update_seaice_list

  CHARACTER(len=*), PARAMETER :: module_name = 'mo_nwp_vdiff_interface'


  !>
  !! Compute a weighted average over a single dimension of a multidimensional array. By default,
  !! the last dimension is averaged. Optionally, the average can be performed for an arbitrary
  !! power `pow` of the array entries, \f$ \bar{x}(...) = (\sum_n c(...,n) x(...,n)^{\mathrm{pow}}
  !! )^{1/\mathrm{pow}} \f$.
  !!
  !! The 3D overload is OpenMP orphaned, and may be called in a parallel region to facilitate
  !! block-level worksharing.
  !!
  !! \warning The weight vector is assumed to sum to 1.
  INTERFACE weighted_average
    MODULE PROCEDURE weighted_average_2d
    MODULE PROCEDURE weighted_average_3d
  END INTERFACE weighted_average

CONTAINS

  !>
  !! Interface routine to the VDIFF turbulence scheme, coupled to the JSBACH LSS.
  !!
  !! This routine computes vertical diffusion of heat, moisture, and other tracers through the
  !! atmosphere. It also couples the land- and sea surfaces to the atmosphere dynamics (via
  !! boundary conditions).
  !!
  !! The routine takes the current state of the atmosphere (through the `nh_` arguments) and
  !! of land and sea (through `diag_`, `prog_`, and `mem`). It computes the turbulent diffusion
  !! coefficients, sets up the linear system to propagate heat, moisture, and tracer contents
  !! forward in time, and calls the LSS and a sea scheme to provide the surface boundary
  !! conditions. The routine then solves the linear system and computes tendencies. It also
  !! provides surface temperatures and albedos, as well as latent and sensible heat fluxes.
  !!
  SUBROUTINE nwp_vdiff ( &
        & datetime_now, delta_time, patch, ccycle_config, vdiff_config, nh_prog, nh_prog_rcf, &
        & nh_diag, nh_metrics, phy_diag, ext_data, diag_lnd, prog_lnd_new, prog_wtr_now, &
        & prog_wtr_new, mem, phy_tend, initialize, lacc &
      )

    TYPE(datetime), POINTER, INTENT(IN) :: datetime_now !< Current time.
    REAL(wp), INTENT(IN) :: delta_time !< Time interval.
    TYPE(t_patch), TARGET, INTENT(IN) :: patch !< Current patch.
    TYPE(t_ccycle_config), INTENT(IN) :: ccycle_config !< Carbon-cycle configuration.
    TYPE(t_vdiff_config), INTENT(IN) :: vdiff_config !< vdiff configuration.
    TYPE(t_nh_prog), INTENT(INOUT) :: nh_prog !< Prognostic variables on current patch.
    TYPE(t_nh_prog), INTENT(INOUT) :: nh_prog_rcf !< Prognostic variables (reduced calling freq).
    TYPE(t_nh_diag), INTENT(INOUT) :: nh_diag !< Diagnostic variables on current patch.
    TYPE(t_nh_metrics), INTENT(IN) :: nh_metrics !< Geometry of the patch.
    TYPE(t_nwp_phy_diag), INTENT(INOUT) :: phy_diag !< Diagnostic physics variables on current patch.
    TYPE(t_external_data), TARGET, INTENT(IN) :: ext_data !< External data for the patch.
    TYPE(t_lnd_diag), INTENT(INOUT) :: diag_lnd
    !< NWP LSS diagnostic land variables.
    TYPE(t_lnd_prog), INTENT(INOUT) :: prog_lnd_new
    !< NWP LSS prognostic land variables (time `t+1`).
    TYPE(t_wtr_prog), INTENT(IN) :: prog_wtr_now
    !< NWP LSS prognostic water variables (time `t`).
    TYPE(t_wtr_prog), INTENT(INOUT) :: prog_wtr_new
    !< NWP LSS prognostic water variables (time `t+1`).
    TYPE(t_nwp_vdiff_state), INTENT(INOUT) :: mem !< vdiff and jsbach state.
    TYPE(t_nwp_phy_tend), INTENT(INOUT) :: phy_tend
    !< NWP physics tendencies.
    LOGICAL, OPTIONAL, INTENT(IN) :: initialize
    !< Propagate land for 1s to get initial fluxes and temperatures. No update of atmospheric
    !! prognostic variables.
    LOGICAL, OPTIONAL, INTENT(IN) :: lacc

    !
    ! Local constants
    !
    INTEGER, PARAMETER :: n_essential_tracers = 3 !< Number of essential tracers.

    !
    ! Local variables
    !

    REAL(wp), TARGET :: zero2d(nproma, patch%nblks_c) !< All-zero 2D field.

    INTEGER :: ktrac !< number of additional tracers.

    !> Grid-box fraction of each surface class [1].
    REAL(wp) :: fr_sfc(nproma, patch%nblks_c, SFC_NUM)
    !> Grid-box fraction of each surface type [1].
    REAL(wp) :: fr_sft(nproma, patch%nblks_c, SFT_NUM)

    !> Temperature of each surface class [K].
    REAL(wp) :: temp_sfc(nproma, patch%nblks_c, SFC_NUM)

    REAL(wp), CONTIGUOUS, POINTER :: ocean_u(:,:) !< Ocean surface velocity (zonal).
    REAL(wp), CONTIGUOUS, POINTER :: ocean_v(:,:) !< Ocean surface velocity (meridional).

    !> Wind speed at lowest model level [m/s].
    REAL(wp) :: wind_lowest(nproma, patch%nblks_c)

    !> Cloud water content (liquid + ice).
    REAL(wp) :: cloud_water_total(nproma, patch%nlev, patch%nblks_c)
    !> Rain at surface [kg/m**2/s].
    REAL(wp), TARGET :: rain_srf(nproma, patch%nblks_c)
    !> Snow at surface [kg/m**2/s].
    REAL(wp), TARGET :: snow_srf(nproma, patch%nblks_c)

    !> Surface emission of tracers [kg/m**2/s].
    REAL(wp) :: tracer_srf_emission(nproma, ntracer-n_essential_tracers, patch%nblks_c)

    !> Convective velocity scale (grid-box mean) [m/s].
    REAL(wp) :: wstar(nproma, patch%nblks_c)
    !> Saturation specific humidity per surface class [kg/kg].
    REAL(wp) :: qsat_sfc(nproma, patch%nblks_c, SFC_NUM)
    !> Height of the top of the dry convective boundary layer [m].
    REAL(wp) :: height_top_dry_cbl(nproma, patch%nblks_c)
    !> Moist Richardson number [1].
    REAL(wp) :: ri_number(nproma, patch%nlev, patch%nblks_c)
    !> Moist surface Richardson number per surface class [1].
    REAL(wp) :: ri_number_sfc(nproma, patch%nblks_c, SFC_NUM)
    !> Mixing length at layer interfaces (excluding surface) [m].
    REAL(wp) :: mixing_length(nproma, patch%nlev, patch%nblks_c)

    !> Prefactor for exchange coefficients (effectively: `1.5 * air density * delta_t`)
    !! [kg*s/m**3].
    REAL(wp) :: prefactor_exchange(nproma, patch%nblks_c)
    !> Exchange coefficient for water variance [m**2/s].
    REAL(wp) :: exchange_coeff_water_var(nproma, patch%nlev, patch%nblks_c)
    !> Exchange coefficient for TTE [m**2/s].
    REAL(wp) :: exchange_coeff_tte(nproma, patch%nlev, patch%nblks_c)
    !> Exchange coefficient for temperature variance [m**2/s].
    REAL(wp) :: exchange_coeff_temp_var(nproma, patch%nlev, patch%nblks_c)

    !> Tridiagonal matrices for turbulent diffusion [1].
    REAL(wp) :: a_matrices(nproma, patch%nlev, 3, nmatrix_vdiff, patch%nblks_c)
    !> Tridiagonal matrices (bottom layer) for turbulent diffusion of heat and moisture per
    !! surface class [1]. Dimension 4 enumerates heat (1) and moisture (2).
    REAL(wp) :: a_matrices_btm(nproma, 3, SFC_NUM, 2, patch%nblks_c)
    !> Right-hand sides for turbulent diffusion [J/kg for heat, kg/kg for everything else].
    REAL(wp) :: b_rhs(nproma, patch%nlev, nvar_vdiff, patch%nblks_c)
    !> Right-hand sides (bottom layer) for turbulent diffusion of heat and moisture per surface
    !! class [J/kg for heat, kg/kg for moisture].
    REAL(wp) :: b_rhs_btm(nproma, SFC_NUM, 2, patch%nblks_c)

    !> Tendency for wind (zonal) [m/s**2]. Set by 3D Smagorinsky scheme.
    REAL(wp) :: ddt_u_smag(nproma, patch%nlev, patch%nblks_c)
    !> Tendency for wind (meridional) [m/s**2]. Set by 3D Smagorinsky scheme.
    REAL(wp) :: ddt_v_smag(nproma, patch%nlev, patch%nblks_c)
    !> Tendency for wind (vertical) [m/s**2]. Set by 3D Smagorinsky scheme.
    REAL(wp) :: ddt_w_smag(nproma, patch%nlev+1, patch%nblks_c)
    !> Tendency for temperature due to horizontal diffusion [K/s]. Set by 3D Smagorinsky scheme.
    REAL(wp) :: ddt_horiz_temp(nproma, patch%nlev, patch%nblks_c)
    !> Tendency for specific humidity due to horizontal diffusion [kg/kg/s]. Set by 3D Smagorinsky
    !! scheme.
    REAL(wp) :: ddt_horiz_qv(nproma, patch%nlev, patch%nblks_c)
    !> Tendency for cloud water due to horizontal diffusion [kg/kg/s]. Set by 3D Smagorinsky
    !! scheme.
    REAL(wp) :: ddt_horiz_qc(nproma, patch%nlev, patch%nblks_c)
    !> Tendency for cloud ice due to horizontal diffusion [kg/kg/s]. Set by 3D Smagorinsky scheme.
    REAL(wp) :: ddt_horiz_qi(nproma, patch%nlev, patch%nblks_c)

    !> Dry static energy at surface for each surface class [J/kg].
    REAL(wp) :: s_sfc(nproma, patch%nblks_c, SFC_NUM)
    !> Dry static energy in atmosphere [J/kg].
    REAL(wp) :: s_atm(nproma, patch%nlev, patch%nblks_c)

    !> Variance of virtual potential temperature (intermediate value) [K**2].
    REAL(wp) :: theta_v_var_intermediate(nproma, patch%nlev, patch%nblks_c)
    !> Total turbulent energy at intermediate time step [J/kg?].
    REAL(wp) :: total_turbulence_energy_intermediate(nproma, patch%nlev, patch%nblks_c)

    !> Exchange coefficient for heat [1]?
    REAL(wp) :: ch_sfc(nproma, patch%nblks_c, SFC_NUM)

    REAL(wp), TARGET :: bn_sfc(nproma, patch%nblks_c, SFC_NUM)
    REAL(wp), TARGET :: bhn_sfc(nproma, patch%nblks_c, SFC_NUM)
    REAL(wp) :: bm_sfc(nproma, patch%nblks_c, SFC_NUM)
    REAL(wp) :: bh_sfc(nproma, patch%nblks_c, SFC_NUM)

    !> CO2 concentration at surface [kg/kg].
    REAL(wp), TARGET :: co2_concentration_srf(nproma, patch%nblks_c)

    !> Ratio of bottom-layer to surface density times bottom-layer thickness Delta_z [m].
    REAL(wp) :: rho_ratio_delta_z(nproma, patch%nblks_c)

    !> Richtmyer E for heat per surface class [1].
    REAL(wp) :: t_acoef(nproma, SFC_NUM)
    !> Richtmyer F for heat per surface class [J/kg].
    REAL(wp) :: t_bcoef(nproma, SFC_NUM)
    !> Richtmyer E for humidity per surface class [1].
    REAL(wp) :: q_acoef(nproma, SFC_NUM)
    !> Richtmyer F for humidity per surface class [kg/kg].
    REAL(wp) :: q_bcoef(nproma, SFC_NUM)
    !> Richtmyer E for wind [1].
    REAL(wp) :: uv_acoef(nproma, SFC_NUM)
    !> Richtmyer F for zonal wind [m/s].
    REAL(wp) :: u_bcoef(nproma, SFC_NUM)
    !> Richtmyer F for meridional wind [m/s].
    REAL(wp) :: v_bcoef(nproma, SFC_NUM)

    !> Surface radiation fluxes.
    TYPE(t_nwp_vdiff_surface_rad_fluxes) :: flx_rad

    !> Drag coefficients [?].
    REAL(wp) :: drag_coef(nproma, SFC_NUM)

    !> Effective surface temperature for radiation [K].
    REAL(wp) :: t_eff_sft(nproma, patch%nblks_c, SFT_NUM)

    !> Saturation specific humidity over each surface type [kg/kg].
    REAL(wp) :: qsat_sft(nproma, patch%nblks_c, SFT_NUM)

    !> Dry static energy over each surface type [J/kg].
    REAL(wp) :: s_sft(nproma, patch%nblks_c, SFT_NUM)

    !> 2m temperature over each surface type [K].
    REAL(wp) :: t2m_sft(nproma, patch%nblks_c, SFT_NUM)
    !> 2m dew point over each surface type [K].
    REAL(wp) :: td2m_sft(nproma, patch%nblks_c, SFT_NUM)
    !> 2m relative humidity over each surface type [1].
    REAL(wp) :: rh2m_sft(nproma, patch%nblks_c, SFT_NUM)
    !> 2m specific humidity over each surface type [kg/kg].
    REAL(wp) :: qv2m_sft(nproma, patch%nblks_c, SFT_NUM)
    !> Zonal 10m wind over each surface type [m/s].
    REAL(wp) :: u10m_sft(nproma, patch%nblks_c, SFT_NUM)
    !> Meridional 10m wind over each surface type [m/s].
    REAL(wp) :: v10m_sft(nproma, patch%nblks_c, SFT_NUM)
    !> 10m wind speed over each surface type [m/s].
    REAL(wp) :: wind_10m_sft(nproma, patch%nblks_c, SFT_NUM)

    !> Effective specific humidity over each surface type [kg/kg].
    REAL(wp) :: qv_sft(nproma, patch%nblks_c, SFT_NUM)

    !> Evapotranspiration from each surface type [kg/m**2/s].
    REAL(wp), TARGET :: evapo_sft(nproma, patch%nblks_c, SFT_NUM)
    !> Potential evaporation from land [kg/m**2/s].
    REAL(wp) :: evapo_potential(nproma, patch%nblks_c)

    !> Ground heat flux [W/m**2].
    REAL(wp) :: flx_heat_ground(nproma, patch%nblks_c)
    !> Ground heat cpacity [J/K/m**2].
    REAL(wp) :: cap_heat_ground(nproma, patch%nblks_c)

    !> Latent heat flux over each surface type [W/m**2].
    REAL(wp), TARGET :: flx_heat_latent_sft(nproma, patch%nblks_c, SFT_NUM)
    !> Sensible heat flux over each surface type [W/m**2].
    REAL(wp), TARGET :: flx_heat_sensible_sft(nproma, patch%nblks_c, SFT_NUM)

    !> Zonal momentum flux over each surface type [N/m**2].
    REAL(wp), TARGET :: flx_mom_u_sft(nproma, patch%nblks_c, SFT_NUM)
    !> Meridional momentum flux over each surface type [N/m**2].
    REAL(wp), TARGET :: flx_mom_v_sft(nproma, patch%nblks_c, SFT_NUM)

    !> Heating rate used to melt snow on canopy [W/m**2].
    REAL(wp) :: Q_snowcanopymelt(nproma, patch%nblks_c)

    !> Albedos per radiation and surface type [1].
    TYPE(t_nwp_vdiff_albedos) :: alb

    !> Column integral of kinetic energy dissipation [W/m**2].
    REAL(wp) :: ekin_dissipation(nproma, patch%nblks_c)
    !> Tendency for zonal wind [m/s**2].
    REAL(wp) :: ddt_u(nproma, patch%nlev, patch%nblks_c)
    !> Tendency for meridional wind [m/s**2].
    REAL(wp) :: ddt_v(nproma, patch%nlev, patch%nblks_c)
    !> Heating rate in each layer [W/m**2]
    REAL(wp) :: ddt_Q(nproma, patch%nlev, patch%nblks_c)
    !> Tracer tendencies [kg/kg/s].
    REAL(wp) :: ddt_tracer(nproma, patch%nlev, patch%nblks_c, ntracer)
    !> Grid-box mean of roughness length for momentum [m].
    REAL(wp) :: z0m_gbm(nproma, patch%nblks_c)

    REAL(wp) :: flx_humidity(nproma) !< Block: Humidity flux [kg/m**2/s].
    REAL(wp) :: flx_sensible(nproma) !< Block: Sensible heat flux [W/m**2].
    REAL(wp) :: flx_mom_u(nproma) !< Block: Zonal momentum flux [N/m**2].
    REAL(wp) :: flx_mom_v(nproma) !< Block: Meridional momentum flux [N/m**2].

    REAL(wp) :: temp_srf_old(nproma, patch%nblks_c)

    INTEGER :: i_startblk !< Start block index.
    INTEGER :: i_endblk !< End block index.
    INTEGER :: i_blk !< Block index.
    INTEGER :: ics !< Start cell index.
    INTEGER :: ice !< End cell index.
    INTEGER :: ic !< Cell index.
    INTEGER :: isft !< Surface type index.
    INTEGER :: isfc !< Surface class index.
    INTEGER :: kl !< Level index.

    REAL(wp), CONTIGUOUS, POINTER :: b_neutral(:,:,:)
    REAL(wp), CONTIGUOUS, POINTER :: p_graupel_gsp_rate(:,:)
    REAL(wp), CONTIGUOUS, POINTER :: p_ice_gsp_rate(:,:)
    REAL(wp), CONTIGUOUS, POINTER :: p_hail_gsp_rate(:,:)

    REAL(wp), CONTIGUOUS, POINTER :: p_condhf_ice_blk(:)
    REAL(wp), CONTIGUOUS, POINTER :: p_meltpot_ice_blk(:)

    LOGICAL :: linit
    LOGICAL :: lis_coupled_to_ocean

    !
    ! Subroutine start
    !

    CALL assert_acc_device_only ('nwp_vdiff', lacc)
    CALL assert_lacc_equals_i_am_accel_node ('nwp_vdiff', lacc, i_am_accel_node)

    ! Asynchronous data regions are a too recent feature. We have to resort to unstructured ones.
    ! This crutch ensures that we don't forget to delete any variable. The compiler complains when
    ! the line gets too long, so we have to split it into chunks, too. Yay...
#   define LIST_CREATE1 \
        zero2d, \
        fr_sfc, \
        fr_sft, \
        temp_sfc, \
        ocean_u, \
        ocean_v, \
        wind_lowest, \
        tracer_srf_emission, \
        cloud_water_total, \
        rain_srf, \
        snow_srf, \
        wstar, \
        qsat_sfc, \
        height_top_dry_cbl, \
        ri_number, \
        ri_number_sfc
#   define LIST_CREATE2 \
        mixing_length, \
        prefactor_exchange, \
        exchange_coeff_water_var, \
        exchange_coeff_tte, \
        exchange_coeff_temp_var, \
        a_matrices, \
        a_matrices_btm, \
        b_rhs, \
        b_rhs_btm, \
        ddt_u_smag, \
        ddt_v_smag, \
        ddt_w_smag, \
        ddt_horiz_temp, \
        ddt_horiz_qv, \
        ddt_horiz_qc, \
        ddt_horiz_qi, \
        s_sfc, \
        s_atm, \
        theta_v_var_intermediate, \
        total_turbulence_energy_intermediate
#   define LIST_CREATE3 \
        ch_sfc, \
        bn_sfc, \
        bhn_sfc, \
        bm_sfc, \
        bh_sfc, \
        co2_concentration_srf, \
        rho_ratio_delta_z, \
        t_acoef, \
        t_bcoef, \
        q_acoef, \
        q_bcoef, \
        uv_acoef, \
        u_bcoef, \
        v_bcoef, \
        flx_rad, \
        drag_coef, \
        t_eff_sft, \
        qsat_sft, \
        s_sft, \
        t2m_sft
#   define LIST_CREATE4 \
        td2m_sft, \
        rh2m_sft, \
        qv2m_sft, \
        u10m_sft, \
        v10m_sft, \
        wind_10m_sft, \
        qv_sft, \
        evapo_sft, \
        evapo_potential, \
        flx_heat_ground, \
        cap_heat_ground, \
        flx_heat_latent_sft, \
        flx_heat_sensible_sft, \
        flx_mom_u_sft, \
        flx_mom_v_sft, \
        Q_snowcanopymelt, \
        alb, \
        ekin_dissipation, \
        ddt_u, \
        ddt_v, \
        ddt_Q
#   define LIST_CREATE5 \
        ddt_tracer, \
        z0m_gbm, \
        flx_humidity, \
        flx_sensible, \
        flx_mom_u, \
        flx_mom_v, \
        temp_srf_old
    !$ACC ENTER DATA ASYNC(1) &
    !$ACC   CREATE(LIST_CREATE1) &
    !$ACC   CREATE(LIST_CREATE2) &
    !$ACC   CREATE(LIST_CREATE3) &
    !$ACC   CREATE(LIST_CREATE4) &
    !$ACC   CREATE(LIST_CREATE5)

    IF (PRESENT(initialize)) THEN
      linit = initialize
    ELSE
      linit = .FALSE.
    END IF

    lis_coupled_to_ocean = is_coupled_to_ocean()

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    ! The NWP setup has v,c,i,r,s always enabled, in that order. `iqt` is the first tracer not
    ! related to moisture. See `mo_nml_crosscheck`.
    ! TODO: Lacking better information, we diffuse all extra tracers.
    ktrac = MAX(0, ntracer + 1 - iqt)

    !$OMP PARALLEL
      CALL init(zero2d(:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL init(tracer_srf_emission(:,:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
    !$OMP END PARALLEL

    CALL get_surface_type_fractions(patch, ext_data, mem, diag_lnd, fr_sfc, fr_sft)
    CALL get_surface_class_temperature(patch, fr_sft, mem%temp_sft, fr_sfc, temp_sfc)

    !$OMP PARALLEL
      CALL weighted_average(patch, fr_sfc(:,:,:), temp_sfc(:,:,:), temp_srf_old(:,:))

      !$OMP DO PRIVATE(i_blk, ics, ice, ic, kl)
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          ! Total cloud water: ice and liquid water.
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO kl = 1, patch%nlev
            DO ic = ics, ice
              cloud_water_total(ic,kl,i_blk) = &
                  & nh_prog_rcf%tracer(ic,kl,i_blk,iqc) + nh_prog_rcf%tracer(ic,kl,i_blk,iqi)
            END DO
          END DO

          ! Consider natural CO2 emissions to close the carbon cycle (if enabled).
          IF (ico2 > 0) THEN
            !$ACC LOOP GANG VECTOR
            DO ic = ics, ice
              tracer_srf_emission(ic,ico2 - iqt + 1,i_blk) = &
                  & mem%flx_co2_natural_land(ic,i_blk) * fr_sft(ic,i_blk,SFT_LAND) &
                  & + mem%sea_state%flx_co2_natural_sea(ic,i_blk) &
                  &   * (fr_sft(ic,i_blk,SFT_SWTR) + fr_sft(ic,i_blk,SFT_SICE))
            END DO
          END IF
        !$ACC END PARALLEL
      END DO
    !$OMP END PARALLEL

    ocean_u => if_associated(mem%sea_state%ocean_u, zero2d)
    ocean_v => if_associated(mem%sea_state%ocean_v, zero2d)

    p_graupel_gsp_rate => if_associated(phy_diag%graupel_gsp_rate, zero2d)
    p_ice_gsp_rate => if_associated(phy_diag%ice_gsp_rate, zero2d)
    p_hail_gsp_rate => if_associated(phy_diag%hail_gsp_rate, zero2d)

    !$ACC ENTER DATA ASYNC(1) &
    !$ACC   ATTACH(ocean_u, ocean_v, p_graupel_gsp_rate, p_ice_gsp_rate, p_hail_gsp_rate)

    !$OMP PARALLEL
      !$OMP DO PRIVATE(i_blk, ics, ice, ic)
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO ic = ics, ice
            rho_ratio_delta_z(ic,i_blk) = &
                & nh_prog%rho(ic,patch%nlev,i_blk) &
                & * (nh_metrics%z_ifc(ic,patch%nlev,i_blk) &
                &    - nh_metrics%z_ifc(ic,patch%nlev+1,i_blk) &
                & ) / (nh_diag%pres_sfc(ic,i_blk)/(rd*temp_srf_old(ic,i_blk) * ( &
                &      1._wp + vtmpc1*diag_lnd%qv_s(ic,i_blk) &
                &      - cloud_water_total(ic,patch%nlev,i_blk))))


            wind_lowest(ic,i_blk) = SQRT( &
                & nh_diag%u(ic,patch%nlev,i_blk)**2 + nh_diag%v(ic,patch%nlev,i_blk)**2)

            ! We treat everything that is not liquid water as snow.
            rain_srf(ic,i_blk) = phy_diag%rain_gsp_rate(ic,i_blk) &
                & + phy_diag%rain_con_rate(ic,i_blk)
            snow_srf(ic,i_blk) = phy_diag%snow_gsp_rate(ic,i_blk) &
                & + phy_diag%snow_con_rate(ic,i_blk) + p_hail_gsp_rate(ic,i_blk) &
                & + p_ice_gsp_rate(ic,i_blk) + p_graupel_gsp_rate(ic,i_blk)
          END DO
        !$ACC END PARALLEL
      END DO

      CALL flx_rad%init(patch, phy_diag)
      CALL alb%init(nproma, patch%nblks_c)
    !$OMP END PARALLEL

    CALL get_surface_co2_concentration(patch, ccycle_config, nh_prog_rcf%tracer(:,:,:,:), &
        & co2_concentration_srf)

    ! Routine queues on async queue 1. No need to wait.
    CALL vdiff_down( &
        & kbdim=nproma, &
        & nblks_c=patch%nblks_c, &
        & nblks_v=patch%nblks_v, &
        & nblks_e=patch%nblks_e, &
        & klev=patch%nlev, &
        & klevm1=patch%nlev-1, &
        & klevp1=patch%nlev+1, &
        & ktrac=ktrac, &
        & ksfc_type=SFC_NUM, &
        & idx_wtr=SFC_WATER, &
        & idx_ice=SFC_ICE, &
        & idx_lnd=SFC_LAND, &
        & pdtime=delta_time, &
        & pcoriol=patch%cells%f_c(:,:), &
        & patch=patch, &
        & pzf=nh_metrics%z_mc(:,:,:), &
        & pzh=nh_metrics%z_ifc(:,:,:), &
        & pgeom1=nh_metrics%geopot_agl(:,:,:), &
        & pfrc=fr_sfc(:,:,:), &
        & ptsfc_tile=temp_sfc(:,:,:), &
        & pocu=ocean_u(:,:), &
        & pocv=ocean_v(:,:), &
        & ppsfc=nh_diag%pres_sfc(:,:), &
        & pum1=nh_diag%u(:,:,:), &
        & pvm1=nh_diag%v(:,:,:), &
        & pwm1=nh_prog%w(:,:,:), &
        & ptm1=nh_diag%temp(:,:,:), &
        & pqm1=nh_prog_rcf%tracer(:,:,:,iqv), &
        & pxlm1=nh_prog_rcf%tracer(:,:,:,iqc), &
        & pxim1=nh_prog_rcf%tracer(:,:,:,iqi), &
        & pxm1=cloud_water_total(:,:,:), &
        & pxtm1=nh_prog_rcf%tracer(:,:,:,iqt:), &
        & pmair=nh_diag%airmass_new(:,:,:), &
        & rho=nh_prog%rho(:,:,:), &
        & paphm1=nh_diag%pres_ifc(:,:,:), &
        & papm1=nh_diag%pres(:,:,:), &
        & ptvm1=nh_diag%tempv(:,:,:), &
        & paclc=phy_diag%clc(:,:,:), &
        & pxt_emis=tracer_srf_emission(:,:,:), &
        & pthvvar=mem%theta_v_var(:,:,:),&
        & pxvar=mem%total_water_var(:,:,:), &
        & pz0m_tile=mem%z0m_sfc(:,:,:), &
        & ptottem1=mem%total_turbulence_energy(:,:,:), &
        & pcsat=mem%fact_qsat_srf(:,:), &
        & pcair=mem%fact_q_air(:,:), &
        & paz0lh=mem%z0h_land(:,:), &
        & vdiff_config=vdiff_config, &
        & & ! In/Outputs
        & pustar=mem%ustar(:,:), &
        & pwstar_tile=mem%wstar_sfc(:,:,:), &
        & & ! Outputs
        & pwstar=wstar(:,:), &
        & pqsat_tile=qsat_sfc(:,:,:), &
        & phdtcbl=height_top_dry_cbl(:,:), &
        & pri=ri_number(:,:,:), &
        & pri_tile=ri_number_sfc(:,:,:), &
        & pmixlen=mixing_length(:,:,:), &
        & pcfm=mem%exchange_coeff_m(:,:,:), &
        & pcfh=mem%exchange_coeff_h(:,:,:), &
        & pcfv=exchange_coeff_water_var(:,:,:), &
        & pcftotte=exchange_coeff_tte(:,:,:), &
        & pcfthv=exchange_coeff_temp_var(:,:,:), &
        & pcfm_tile=mem%exchange_coeff_m_sfc(:,:,:), &
        & pcfh_tile=mem%exchange_coeff_h_sfc(:,:,:), &
        & aa=a_matrices(:,:,:,:,:), &
        & aa_btm=a_matrices_btm(:,:,:,:,:), &
        & bb=b_rhs(:,:,:,:), &
        & bb_btm=b_rhs_btm(:,:,:,:), &
        & & ! <3D Smagorinsky outputs>
        & ddt_u=ddt_u_smag(:,:,:), &
        & ddt_v=ddt_v_smag(:,:,:), &
        & ddt_w=ddt_w_smag(:,:,:), &
        & ta_hori_tend=ddt_horiz_temp(:,:,:), &
        & qv_hori_tend=ddt_horiz_qv(:,:,:), &
        & ql_hori_tend=ddt_horiz_qc(:,:,:), &
        & qi_hori_tend=ddt_horiz_qi(:,:,:), &
        & & ! </3D Smagorinsky outputs>
        & pfactor_sfc=prefactor_exchange(:,:), &
        & pcpt_tile=s_sfc(:,:,:), &
        & pcptgz=s_atm(:,:,:), &
        & pzthvvar=theta_v_var_intermediate(:,:,:), &
        & pztottevn=total_turbulence_energy_intermediate(:,:,:), &
        & pch_tile=ch_sfc(:,:,:), &
        & pbn_tile=bn_sfc(:,:,:), &
        & pbhn_tile=bhn_sfc(:,:,:), &
        & pbm_tile=bm_sfc(:,:,:), &
        & pbh_tile=bh_sfc(:,:,:) &
      )

    ! vdiff_down does not initialize the INTENT(OUT) tendencies from the Smagorinsky if TTE is chosen :(
    IF (vdiff_config%turb /= VDIFF_TURB_3DSMAGORINSKY) THEN
      !$OMP PARALLEL
        CALL init(ddt_u_smag(:,:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
        CALL init(ddt_v_smag(:,:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
        CALL init(ddt_w_smag(:,:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
        CALL init(ddt_horiz_temp(:,:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
        CALL init(ddt_horiz_qv(:,:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
        CALL init(ddt_horiz_qc(:,:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
        CALL init(ddt_horiz_qi(:,:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
      !$OMP END PARALLEL
    END IF

    ! TODO: Call this once per day only.
    CALL update_earth_declination(datetime_now)

    CALL sea_model_update_sst( &
        & patch, mem%sea_state, datetime_now, ext_data%atm_td%sst_m(:,:,:), &
        & diag_lnd%t_seasfc(:,:) &
      )

#ifndef __NO_JSBACH__
    CALL jsbach_start_timestep(patch%id, datetime_now, delta_time)
#endif

    !$OMP PARALLEL PRIVATE(i_blk, ic, ics, ice, isfc, isft, kl) &
    !$OMP   PRIVATE(drag_coef, t_acoef, t_bcoef, q_acoef, q_bcoef, uv_acoef, u_bcoef, v_bcoef) &
    !$OMP   PRIVATE(flx_humidity, flx_sensible, flx_mom_u, flx_mom_v, b_neutral) &
    !$OMP   PRIVATE(p_condhf_ice_blk, p_meltpot_ice_blk)
      !$OMP DO
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO ic = ics, ice
            drag_coef(ic,SFC_LAND) = grav * prefactor_exchange(ic,i_blk) &
                & * mem%exchange_coeff_h_sfc(ic,i_blk,SFC_LAND)
            drag_coef(ic,SFC_WATER) = grav * prefactor_exchange(ic,i_blk) &
                & * mem%exchange_coeff_h_sfc(ic,i_blk,SFC_WATER)
            drag_coef(ic,SFC_ICE) = grav * prefactor_exchange(ic,i_blk) &
                & * mem%exchange_coeff_h_sfc(ic,i_blk,SFC_ICE)
            ch_sfc(ic,i_blk,SFC_LAND) = &
                & MERGE(ch_sfc(ic,i_blk,SFC_LAND), 1._wp, ext_data%atm%fr_land(ic,i_blk) > 0._wp)
          END DO
        !$ACC END PARALLEL

        CALL matrix_to_richtmyer_coeff( &
            & jcs=ics, &
            & kproma=ice, &
            & klev=patch%nlev, &
            & ksfc_type=SFC_NUM, &
            & idx_lnd=SFC_LAND, &
            & aa=a_matrices(:,:,:,imh_vdiff:imqv_vdiff,i_blk), &
            & bb=b_rhs(:,:,ih_vdiff:iqv_vdiff,i_blk), &
            & pdtime=delta_time, &
            & delz=rho_ratio_delta_z(:,i_blk), &
            & aa_btm=a_matrices_btm(:,:,:,:,i_blk), &
            & bb_btm=b_rhs_btm(:,:,:,i_blk), &
            & pen_h=t_acoef(:,:), &
            & pfn_h=t_bcoef(:,:), &
            & pen_qv=q_acoef(:,:), &
            & pfn_qv=q_bcoef(:,:), &
            & pcair=mem%fact_q_air(:,i_blk), &
            & pcsat=mem%fact_qsat_srf(:,i_blk) &
          )

        CALL vdiff_get_richtmyer_coeff_momentum( &
            & jcs=ics, &
            & kproma=ice, &
            & klev=patch%nlev, &
            & ksfc_type=SFC_NUM, &
            & aa=a_matrices(:,:,:,:,i_blk), &
            & bb=b_rhs(:,:,:,i_blk), &
            & pfactor_sfc=prefactor_exchange(:,i_blk), &
            & pcfm_tile=mem%exchange_coeff_m_sfc(:,i_blk,:), &
            & pmair=nh_diag%airmass_new(:,:,i_blk), &
            & pen_uv=uv_acoef(:,:), &
            & pfn_u=u_bcoef(:,:), &
            & pfn_v=v_bcoef(:,:) &
          )

        CALL vdiff_get_tke ( &
            & jcs=ics, &
            & kproma=ice, &
            & klev=patch%nlev, &
            & vdiff_config=vdiff_config, &
            & ptotte=mem%total_turbulence_energy(:,:,i_blk), &
            & pri=ri_number(:,:,i_blk), &
            & tke=nh_prog_rcf%tke(:,:,i_blk) &
          )

        ! Replace the vapor in the lowest atmospheric layer by total water (qv+qc) so the surface
        ! sees condensed water and gives consistent gradients, since any water would evaporate /
        ! condense instantly until the saturation specific humidity at the surface is reached. For
        ! consistency, dry static energy is reduced by the energy used to evaporate the liquid. The
        ! humidity flux will be applied to vapor and saturation adjustment takes care of the actual
        ! evaporation/condensation in the lowest level. We need to be careful to use the same
        ! latent heat that the saturation adjustment applies later.
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP SEQ
          DO isfc = 1, SFC_NUM
            !$ACC LOOP GANG VECTOR
            DO ic = ics, ice
              q_bcoef(ic,isfc) = q_bcoef(ic,isfc) + &
                  & vdiff_implfact * b_rhs(ic,patch%nlev,iqc_vdiff,i_blk)
              t_bcoef(ic,isfc) = t_bcoef(ic,isfc) - &
                  & latent_heat_vaporization(nh_diag%temp(ic,patch%nlev,i_blk)) &
                  & * vdiff_implfact * b_rhs(ic,patch%nlev,iqc_vdiff,i_blk)
            END DO
          END DO
        !$ACC END PARALLEL

#ifndef __NO_JSBACH__
        !$ACC WAIT

        CALL jsbach_interface( &
            & model_id=patch%id, &
            & iblk=i_blk, &
            & ics=ics, &
            & ice=ice, &
            & current_datetime=datetime_now,&
            & dtime=MERGE(1._wp, delta_time, linit), &
            & steplen=delta_time, &
            & t_air=nh_diag%temp(ics:ice, patch%nlev, i_blk), &
            & q_air=nh_prog_rcf%tracer(ics:ice, patch%nlev, i_blk, iqv), &
            & rain=rain_srf(ics:ice,i_blk), &
            & snow=snow_srf(ics:ice,i_blk), &
            & wind_air=wind_lowest(ics:ice,i_blk), &
            & wind_10m=phy_diag%sp_10m(ics:ice,i_blk), &
            & lw_srf_down=flx_rad%flx_lw_down(ics:ice,i_blk), &
            & swvis_srf_down=flx_rad%flx_vis_down(ics:ice,i_blk), &
            & swnir_srf_down=flx_rad%flx_nir_down(ics:ice,i_blk), &
            & swpar_srf_down=flx_rad%flx_par_down(ics:ice,i_blk), &
            & fract_par_diffuse=flx_rad%fr_par_diffuse(ics:ice,i_blk), &
            & press_srf=nh_diag%pres_sfc(ics:ice,i_blk), &
            & drag_srf=drag_coef(ics:ice,SFC_LAND), &
            & drag_wtr=drag_coef(ics:ice,SFC_WATER), &
            & drag_ice=drag_coef(ics:ice,SFC_ICE), &
            & t_acoef=t_acoef(ics:ice,SFC_LAND), &
            & t_bcoef=t_bcoef(ics:ice,SFC_LAND), &
            & q_acoef=q_acoef(ics:ice,SFC_LAND), &
            & q_bcoef=q_bcoef(ics:ice,SFC_LAND), &
            & t_acoef_wtr=t_acoef(ics:ice,SFC_WATER), &
            & t_bcoef_wtr=t_bcoef(ics:ice,SFC_WATER), &
            & q_acoef_wtr=q_acoef(ics:ice,SFC_WATER), &
            & q_bcoef_wtr=q_bcoef(ics:ice,SFC_WATER), &
            & t_acoef_ice=t_acoef(ics:ice,SFC_ICE), &
            & t_bcoef_ice=t_bcoef(ics:ice,SFC_ICE), &
            & q_acoef_ice=q_acoef(ics:ice,SFC_ICE), &
            & q_bcoef_ice=q_bcoef(ics:ice,SFC_ICE), &
            & pch=ch_sfc(ics:ice,i_blk,SFC_LAND), &
            & cos_zenith_angle=phy_diag%cosmu0(ics:ice,i_blk), &
            & CO2_air=co2_concentration_srf(ics:ice,i_blk), &
            & & ! Outputs
            & t_srf=mem%temp_sft(ics:ice,i_blk,SFT_LAND), &
            & t_lwtr=mem%temp_sft(ics:ice,i_blk,SFT_LWTR), &
            & t_lice=mem%temp_sft(ics:ice,i_blk,SFT_LICE), &
            & t_eff_srf=t_eff_sft(ics:ice,i_blk,SFT_LAND), &
            & qsat_srf=qsat_sft(ics:ice,i_blk,SFT_LAND), &
            & qsat_lwtr=qsat_sft(ics:ice,i_blk,SFT_LWTR), &
            & qsat_lice=qsat_sft(ics:ice,i_blk,SFT_LICE), &
            & s_srf=s_sft(ics:ice,i_blk,SFT_LAND), &
            & s_lwtr=s_sft(ics:ice,i_blk,SFT_LWTR), &
            & s_lice=s_sft(ics:ice,i_blk,SFT_LICE), &
            & fact_q_air=mem%fact_q_air(ics:ice,i_blk), &
            & fact_qsat_srf=mem%fact_qsat_srf(ics:ice,i_blk), &
            & evapotrans=evapo_sft(ics:ice,i_blk,SFT_LAND), &
            & evapopot=evapo_potential(ics:ice,i_blk), &
            & evapo_wtr=evapo_sft(ics:ice,i_blk,SFT_LWTR), &
            & evapo_ice=evapo_sft(ics:ice,i_blk,SFT_LICE), &
            & CO2_flux=mem%flx_co2_natural_land(ics:ice,i_blk), &
            & latent_hflx=flx_heat_latent_sft(ics:ice,i_blk,SFT_LAND), &
            & latent_hflx_wtr=flx_heat_latent_sft(ics:ice,i_blk,SFT_LWTR), &
            & latent_hflx_ice=flx_heat_latent_sft(ics:ice,i_blk,SFT_LICE), &
            & sensible_hflx=flx_heat_sensible_sft(ics:ice,i_blk,SFT_LAND), &
            & sensible_hflx_wtr=flx_heat_sensible_sft(ics:ice,i_blk,SFT_LWTR), &
            & sensible_hflx_ice=flx_heat_sensible_sft(ics:ice,i_blk,SFT_LICE), &
            & grnd_hflx=flx_heat_ground(ics:ice,i_blk), &
            & grnd_hcap=cap_heat_ground(ics:ice,i_blk), &
            & rough_h_srf=mem%z0h_land(ics:ice,i_blk), &
            & rough_m_srf=mem%z0m_sfc(ics:ice,i_blk,SFC_LAND), &
            & q_snocpymlt=Q_snowcanopymelt(ics:ice,i_blk), &
            & alb_vis_dir=alb%alb_vis_dir(ics:ice,i_blk,SFT_LAND), &
            & alb_vis_dif=alb%alb_vis_dif(ics:ice,i_blk,SFT_LAND), &
            & alb_nir_dir=alb%alb_nir_dir(ics:ice,i_blk,SFT_LAND), &
            & alb_nir_dif=alb%alb_nir_dif(ics:ice,i_blk,SFT_LAND), &
            & albedo_lwtr=alb%alb_vis_dif(ics:ice,i_blk,SFT_LWTR), &
            & albedo_lice=alb%alb_vis_dif(ics:ice,i_blk,SFT_LICE), &
            & ice_fract_lake=mem%fr_ice_on_lake(ics:ice,i_blk) &
          )
#endif

        ! condhf_ice and meltpot_ice are unallocated for uncoupled runs.
        IF (lis_coupled_to_ocean) THEN
          p_condhf_ice_blk => diag_lnd%condhf_ice(:,i_blk)
          p_meltpot_ice_blk => diag_lnd%meltpot_ice(:,i_blk)
        ELSE
          p_condhf_ice_blk => NULL()
          p_meltpot_ice_blk => NULL()
        END IF

        CALL sea_model( &
            & iblk=i_blk, &
            & ics=ics, &
            & ice=ice, &
            & dtime=MERGE(1._wp, delta_time, linit), &
            & steplen=delta_time, &
            & ext_data=ext_data, &
            & rain=rain_srf(:,i_blk), &
            & snow=snow_srf(:,i_blk), &
            & latent_hflx_ice_old=mem%flx_heat_latent_sft(:,i_blk,SFT_SICE), &
            & sensible_hflx_ice_old=mem%flx_heat_sensible_sft(:,i_blk,SFT_SICE), &
            & flx_rad=flx_rad, &
            & sea_state=mem%sea_state, &
            & cos_zenith_angle=phy_diag%cosmu0(:,i_blk), &
            & press_srf=nh_diag%pres_sfc(:,i_blk), &
            & wind_10m=phy_diag%sp_10m(:,i_blk), &
            & prefactor_exchange=prefactor_exchange(:,i_blk), &
            & exchange_coeff_h_wtr=mem%exchange_coeff_h_sfc(:,i_blk,SFC_WATER), &
            & exchange_coeff_h_ice=mem%exchange_coeff_h_sfc(:,i_blk,SFC_ICE), &
            & t_acoef_wtr=t_acoef(:,SFC_WATER), &
            & t_bcoef_wtr=t_bcoef(:,SFC_WATER), &
            & q_acoef_wtr=q_acoef(:,SFC_WATER), &
            & q_bcoef_wtr=q_bcoef(:,SFC_WATER), &
            & t_acoef_ice=t_acoef(:,SFC_ICE), &
            & t_bcoef_ice=t_bcoef(:,SFC_ICE), &
            & q_acoef_ice=q_acoef(:,SFC_ICE), &
            & q_bcoef_ice=q_bcoef(:,SFC_ICE), &
            & prog_wtr_now=prog_wtr_now, &
            & diag_lnd=diag_lnd, &
            & & ! Outputs
            & t_wtr=mem%temp_sft(:,i_blk,SFT_SWTR), &
            & t_ice=mem%temp_sft(:,i_blk,SFT_SICE), &
            & s_wtr=s_sft(:,i_blk,SFT_SWTR), &
            & s_ice=s_sft(:,i_blk,SFT_SICE), &
            & qsat_wtr=qsat_sft(:,i_blk,SFT_SWTR), &
            & qsat_ice=qsat_sft(:,i_blk,SFT_SICE), &
            & evapo_wtr=evapo_sft(:,i_blk,SFT_SWTR), &
            & evapo_ice=evapo_sft(:,i_blk,SFT_SICE), &
            & latent_hflx_wtr=flx_heat_latent_sft(:,i_blk,SFT_SWTR), &
            & latent_hflx_ice=flx_heat_latent_sft(:,i_blk,SFT_SICE), &
            & sensible_hflx_wtr=flx_heat_sensible_sft(:,i_blk,SFT_SWTR), &
            & sensible_hflx_ice=flx_heat_sensible_sft(:,i_blk,SFT_SICE), &
            & conductive_hflx_ice=p_condhf_ice_blk, &
            & melt_potential_ice=p_meltpot_ice_blk, &
            & alb=alb, &
            & prog_wtr_new=prog_wtr_new &
          )

        ! Diagnose surface stress (in N/m**2). This is a mixed-time flux.
        CALL get_surface_stress ( &
            & ics=ics, &
            & ice=ice, &
            & delta_time=delta_time, &
            & prefactor_exchange=prefactor_exchange(:,i_blk), &
            & exchange_coeff_m_sfc=mem%exchange_coeff_m_sfc(:,i_blk,:), &
            & uv_acoef=uv_acoef(:,:), &
            & u_bcoef=u_bcoef(:,:), &
            & v_bcoef=v_bcoef(:,:), &
            & ocean_u=ocean_u(:,i_blk), &
            & ocean_v=ocean_v(:,i_blk), &
            & zero=zero2d(:,1), &
            & umfl_sft=flx_mom_u_sft(:,i_blk,:), &
            & vmfl_sft=flx_mom_v_sft(:,i_blk,:) &
          )

        ! Including Q_snowcanopymelt here gives a slightly different flux from what JSBACH sees,
        ! but the energy has to come from somewhere. And the melted water goes into the ground
        ! reservoir, so it is definitely removing heat from the atmosphere.
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO ic = ics, ice
            flx_heat_sensible_sft(ic, i_blk, SFT_LAND) = &
                & flx_heat_sensible_sft(ic, i_blk, SFT_LAND) + Q_snowcanopymelt(ic,i_blk)
          END DO
        !$ACC END PARALLEL

        CALL weighted_average(fr_sft(:,i_blk,:), evapo_sft(:,i_blk,:), flx_humidity(:))
        CALL weighted_average(fr_sft(:,i_blk,:), flx_heat_sensible_sft(:,i_blk,:), flx_sensible(:))
        CALL weighted_average(fr_sft(:,i_blk,:), flx_mom_u_sft(:,i_blk,:), flx_mom_u(:))
        CALL weighted_average(fr_sft(:,i_blk,:), flx_mom_v_sft(:,i_blk,:), flx_mom_v(:))

        CALL vdiff_update_boundary( &
            & jcs=ics, &
            & kproma=ice, &
            & klev=patch%nlev, &
            & dtime=delta_time, &
            & pmair=nh_diag%airmass_new(:,patch%nlev,i_blk), &
            & shflx=flx_sensible(:), &
            & qflx=flx_humidity(:), &
            & uflx=flx_mom_u(:), &
            & vflx=flx_mom_v(:), &
            & aa=a_matrices(:,:,:,:,i_blk), &
            & aa_btm=a_matrices_btm(:,:,:,:,i_blk), &
            & s_btm=s_atm(:,patch%nlev,i_blk), &
            & q_btm=nh_prog_rcf%tracer(:,patch%nlev,i_blk,iqv), &
            & bb=b_rhs(:,:,:,i_blk) &
          )

        CALL vdiff_up( &
            & jcs=ics, &
            & kproma=ice, &
            & kbdim=nproma, &
            & klev=patch%nlev, &
            & klevm1=patch%nlev-1, &
            & ktrac=ktrac, &
            & ksfc_type=SFC_NUM, &
            & idx_wtr=SFC_WATER, &
            & pdtime=delta_time, &
            & pfrc=fr_sfc(:,i_blk,:), &
            & pcfm_tile=mem%exchange_coeff_m_sfc(:,i_blk,:), &
            & aa=a_matrices(:,:,:,:,i_blk), &
            & pcptgz=s_atm(:,:,i_blk), &
            & pum1=nh_diag%u(:,:,i_blk), &
            & pvm1=nh_diag%v(:,:,i_blk), &
            & ptm1=nh_diag%temp(:,:,i_blk), &
            & pmair=nh_diag%airmass_new(:,:,i_blk), &
            & pqm1=nh_prog_rcf%tracer(:,:,i_blk,iqv), &
            & pxlm1=nh_prog_rcf%tracer(:,:,i_blk,iqc), &
            & pxim1=nh_prog_rcf%tracer(:,:,i_blk,iqi), &
            & pxtm1=nh_prog_rcf%tracer(:,:,i_blk,iqt:), &
            & pgeom1=nh_metrics%geopot_agl(:,:,i_blk), &
            & pztottevn=total_turbulence_energy_intermediate(:,:,i_blk), &
            & vdiff_config=vdiff_config, &
            & bb=b_rhs(:,:,:,i_blk), &
            & pzthvvar=theta_v_var_intermediate(:,:,i_blk), &
            & & ! In/outputs
            & pxvar=mem%total_water_var(:,:,i_blk), &
            & pz0m_tile=mem%z0m_sfc(:,i_blk,:), &
            & pkedisp=ekin_dissipation(:,i_blk), &
            & pute_vdf=ddt_u(:,:,i_blk), &
            & pvte_vdf=ddt_v(:,:,i_blk), &
            & pq_vdf=ddt_Q(:,:,i_blk), &
            & pqte_vdf=ddt_tracer(:,:,i_blk,iqv), &
            & pxlte_vdf=ddt_tracer(:,:,i_blk,iqc), &
            & pxite_vdf=ddt_tracer(:,:,i_blk,iqi), &
            & pxtte_vdf=ddt_tracer(:,:,i_blk,iqt:), &
            & pz0m=z0m_gbm(:,i_blk), &
            & pthvvar=mem%theta_v_var(:,:,i_blk), &
            & ptotte=mem%total_turbulence_energy(:,:,i_blk) &
          )

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO ic = ics, ice
            qv_sft(ic,i_blk,SFT_LAND) = &
                & (1._wp - mem%fact_q_air(ic,i_blk)) * ( &
                &   nh_prog_rcf%tracer(ic,patch%nlev,i_blk,iqv) &
                &   + delta_time * ddt_tracer(ic,patch%nlev,i_blk,iqv)) &
                & + mem%fact_qsat_srf(ic,i_blk) * qsat_sft(ic,i_blk,SFT_LAND)
          END DO

          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO isft = SFT_LWTR, SFT_NUM
            DO ic = ics, ice
              qv_sft(ic,i_blk,isft) = qsat_sft(ic,i_blk,isft)
            END DO
          END DO

          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO kl = 1, patch%nlev
            DO ic = ics, ice
              s_atm(ic,kl,i_blk) = s_atm(ic,kl,i_blk) &
                  & + ddt_Q(ic,kl,i_blk) / nh_diag%airmass_new(ic,kl,i_blk) * delta_time
            END DO
          END DO
        !$ACC END PARALLEL

        ! Diagnose 2m temperature, dew point, and 10m wind.
        DO isft = 1, SFT_NUM
          isfc = SFT_CLASS(isft)

          SELECT CASE (isft)
          CASE (SFT_LAND)
            b_neutral => bhn_sfc
          CASE DEFAULT
            b_neutral => bn_sfc
          END SELECT

          !$ACC DATA ATTACH(b_neutral)
            CALL get_near_surface_temperature ( &
                & ics=ics, &
                & ice=ice, &
                & z_ref=2._wp, &
                & z_lowest=nh_metrics%z_mc(:,patch%nlev,i_blk), &
                & z_srf=nh_metrics%z_ifc(:,patch%nlev+1,i_blk), &
                & ri_number=ri_number_sfc(:,i_blk,isfc), &
                & b_neutral=b_neutral(:,i_blk,isfc), &
                & b_h_actual=bh_sfc(:,i_blk,isfc), &
                & s_lowest=s_atm(:,patch%nlev,i_blk), &
                & s_srf=s_sft(:,i_blk,isft), &
                & temp_ref=t2m_sft(:,i_blk,isft) &
              )

            CALL get_near_surface_dew_point ( &
                & ics=ics, &
                & ice=ice, &
                & z_ref=2._wp, &
                & temp_ref=t2m_sft(:,i_blk,isft), &
                & temp_lowest=nh_diag%temp(:,patch%nlev,i_blk), &
                & pres_lowest=nh_diag%pres(:,patch%nlev,i_blk), &
                & pres_srf=nh_diag%pres_sfc(:,i_blk), &
                & q_lowest=nh_prog_rcf%tracer(:,patch%nlev,i_blk,iqv), &
                & x_lowest=cloud_water_total(:,patch%nlev,i_blk), &
                & tdew_ref=td2m_sft(:,i_blk,isft), &
                & rh_ref=rh2m_sft(:,i_blk,isft), &
                & q_ref=qv2m_sft(:,i_blk,isft) &
              )

            CALL get_near_surface_wind ( &
                & ics=ics, &
                & ice=ice, &
                & z_ref=10._wp, &
                & z_lowest=nh_metrics%z_mc(:,patch%nlev,i_blk), &
                & z_srf=nh_metrics%z_ifc(:,patch%nlev+1,i_blk), &
                & u_lowest=nh_diag%u(:,patch%nlev,i_blk), &
                & v_lowest=nh_diag%v(:,patch%nlev,i_blk), &
                & ocean_u=ocean_u(:,i_blk), &
                & ocean_v=ocean_v(:,i_blk), &
                & ri_number=ri_number_sfc(:,i_blk,isfc), &
                & b_neutral=bn_sfc(:,i_blk,isfc), &
                & b_m_actual=bm_sfc(:,i_blk,isfc), &
                & u_ref=u10m_sft(:,i_blk,isft), &
                & v_ref=v10m_sft(:,i_blk,isft), &
                & sp_ref=wind_10m_sft(:,i_blk,isft) &
              )
          !$ACC END DATA
        END DO

      END DO
    !$OMP END PARALLEL

#ifndef __NO_JSBACH__
    CALL jsbach_finish_timestep (patch%id, datetime_now, delta_time)
#endif

    !$OMP PARALLEL
      ! For everything except land, the effective temperature for radiation is the surface
      ! temperature.
      CALL copy ( &
          & mem%temp_sft(:,:,SFT_LWTR:SFT_NUM), &
          & t_eff_sft(:,:,SFT_LWTR:SFT_NUM), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )

      ! JSBACH only provides a single albedo for lakes. Copy them to all flavors.
      CALL copy(alb%alb_vis_dif(:,:,SFT_LWTR), alb%alb_nir_dif(:,:,SFT_LWTR), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(alb%alb_vis_dif(:,:,SFT_LWTR), alb%alb_nir_dir(:,:,SFT_LWTR), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(alb%alb_vis_dif(:,:,SFT_LWTR), alb%alb_vis_dir(:,:,SFT_LWTR), lacc=.TRUE., opt_acc_async=.TRUE.)

      CALL copy(alb%alb_vis_dif(:,:,SFT_LICE), alb%alb_nir_dif(:,:,SFT_LICE), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(alb%alb_vis_dif(:,:,SFT_LICE), alb%alb_nir_dir(:,:,SFT_LICE), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(alb%alb_vis_dif(:,:,SFT_LICE), alb%alb_vis_dir(:,:,SFT_LICE), lacc=.TRUE., opt_acc_async=.TRUE.)

      ! For consistency, we set the long-wave emissivity of land and lake surfaces to the value
      ! used by JSBACH to deduce net long-wave radiation.
      CALL init(alb%lw_emissivity(:,:,SFT_LAND:SFT_LICE), zemiss_def, lacc=.TRUE., opt_acc_async=.TRUE.)

      ! Update grid-box temperature.
      CALL weighted_average ( &
          & patch, fr_sft(:,:,:), t_eff_sft(:,:,:), prog_lnd_new%t_g(:,:), pow=4._wp &
        )
      CALL weighted_average(patch, fr_sft(:,:,:), mem%temp_sft(:,:,:), diag_lnd%t_s(:,:))

      ! Update surface humidity.
      CALL weighted_average(patch, fr_sft(:,:,:), qv_sft(:,:,:), diag_lnd%qv_s(:,:))

      CALL weighted_average(patch, fr_sft(:,:,:), alb%alb_nir_dif(:,:,:), phy_diag%albnirdif(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), alb%alb_nir_dir(:,:,:), phy_diag%albnirdir(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), alb%alb_vis_dif(:,:,:), phy_diag%albvisdif(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), alb%alb_vis_dir(:,:,:), phy_diag%albvisdir(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), alb%lw_emissivity(:,:,:), phy_diag%lw_emiss(:,:))

      CALL flx_rad%get_albdif ( &
          & patch, phy_diag%albnirdif(:,:), phy_diag%albvisdif(:,:), phy_diag%albdif(:,:) &
        )

      ! Update tendencies.
      !$OMP DO PRIVATE(i_blk, ics, ice, ic, kl)
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO kl = 1, patch%nlev
            DO ic = ics, ice
              phy_tend%ddt_temp_turb(ic,kl,i_blk) = ddt_Q(ic,kl,i_blk) &
                  & / (cvd + (cvv - cvd) * (nh_prog_rcf%tracer(ic,kl,i_blk,iqv) &
                  &    + delta_time * ddt_tracer(ic,kl,i_blk,iqv))) &
                  & / nh_diag%airmass_new(ic,kl,i_blk) &
                  & + ddt_horiz_temp(ic,kl,i_blk)

              ddt_tracer(ic,kl,i_blk,iqv) = ddt_tracer(ic,kl,i_blk,iqv) + ddt_horiz_qv(ic,kl,i_blk)
              ddt_tracer(ic,kl,i_blk,iqc) = ddt_tracer(ic,kl,i_blk,iqc) + ddt_horiz_qc(ic,kl,i_blk)
              ddt_tracer(ic,kl,i_blk,iqi) = ddt_tracer(ic,kl,i_blk,iqi) + ddt_horiz_qi(ic,kl,i_blk)

              phy_tend%ddt_u_turb(ic,kl,i_blk) = ddt_u(ic,kl,i_blk) + ddt_u_smag(ic,kl,i_blk)
              phy_tend%ddt_v_turb(ic,kl,i_blk) = ddt_v(ic,kl,i_blk) + ddt_v_smag(ic,kl,i_blk)
            END DO
          END DO

          IF (ASSOCIATED(phy_tend%ddt_w_turb)) THEN
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
            DO kl = 1, patch%nlev
              DO ic = ics, ice
                phy_tend%ddt_w_turb(ic,kl,i_blk) = ddt_w_smag(ic,kl,i_blk)
              END DO
            END DO
          END IF
        !$ACC END PARALLEL
      END DO

      ! ddt_tracer_turb only stores qv, qc, qi.
      CALL copy ( &
          & ddt_tracer(:,:,:,1:nqtendphy), &
          & phy_tend%ddt_tracer_turb(:,:,:,1:nqtendphy), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )
      CALL copy ( &
          & flx_heat_sensible_sft(:,:,:), &
          & mem%flx_heat_sensible_sft(:,:,:), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )
      CALL copy ( &
          & flx_heat_latent_sft(:,:,:), &
          & mem%flx_heat_latent_sft(:,:,:), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )
      CALL copy ( &
        & evapo_sft(:,:,:), &
        & mem%flx_water_vapor_sft(:,:,:), &
        & opt_acc_async=.TRUE. &
      )

      CALL weighted_average ( &
          & patch, fr_sft(:,:,:), mem%flx_heat_latent_sft(:,:,:), phy_diag%lhfl_s(:,:) &
        )
      CALL weighted_average ( &
          & patch, fr_sft(:,:,:), mem%flx_heat_sensible_sft(:,:,:), phy_diag%shfl_s(:,:) &
        )
      CALL weighted_average ( &
          & patch, fr_sft(:,:,:), mem%flx_water_vapor_sft(:,:,:), phy_diag%qhfl_s(:,:) &
        )
      CALL weighted_average(patch, fr_sft(:,:,:), flx_mom_u_sft(:,:,:), phy_diag%umfl_s(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), flx_mom_v_sft(:,:,:), phy_diag%vmfl_s(:,:))

      CALL weighted_average(patch, fr_sft(:,:,:), t2m_sft(:,:,:), phy_diag%t_2m(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), td2m_sft(:,:,:), phy_diag%td_2m(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), rh2m_sft(:,:,:), phy_diag%rh_2m(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), qv2m_sft(:,:,:), phy_diag%qv_2m(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), u10m_sft(:,:,:), phy_diag%u_10m(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), v10m_sft(:,:,:), phy_diag%v_10m(:,:))
      CALL weighted_average(patch, fr_sft(:,:,:), wind_10m_sft(:,:,:), phy_diag%sp_10m(:,:))

      CALL get_gust_speed ( &
          & patch, phy_diag, nh_diag, nh_metrics, phy_diag%sp_10m, mem%ustar, phy_diag%dyn_gust &
        )

      ! These are not accessible from outside JSBACH. Fill with placeholders.
      CALL copy ( &
          & mem%flx_heat_latent_sft(:,:,SFT_LAND), &
          & phy_diag%lhfl_pl(:,1,:), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )
      CALL init (phy_diag%lhfl_pl(:,2:,:), lacc=.TRUE., opt_acc_async=.TRUE.)

      CALL copy ( &
          & mem%flx_heat_latent_sft(:,:,SFT_LAND), &
          & phy_diag%lhfl_bs(:,:), &
          & lacc=.TRUE., &
          & opt_acc_async=.TRUE. &
        )

      ! AES switches off theta_v and water variance. theta_v is replicated here for consistency.
      ! Water variance is diagnosed below.

      CALL init (mem%theta_v_var(:,:,:), lacc=.TRUE., opt_acc_async=.TRUE.)

      CALL copy (t2m_sft(:,:,SFT_LAND), phy_diag%t_2m_land(:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy (td2m_sft(:,:,SFT_LAND), phy_diag%td_2m_land(:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy (phy_diag%t_2m(:,:), phy_diag%t_tilemin_inst_2m(:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy (phy_diag%t_2m(:,:), phy_diag%t_tilemax_inst_2m(:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
    !$OMP END PARALLEL

    CALL update_nwp_tile_state (patch, delta_time, fr_sft, prog_lnd_new, diag_lnd, phy_diag)

    ! Update state with tendencies.
    IF (.NOT. linit) THEN
      CALL add_tendencies_to_state ( &
          & delta_time, patch, phy_tend, nh_prog, nh_prog_rcf, nh_diag, ddt_tracer &
        )

      CALL get_stddev_saturation_deficit ( &
          & patch=patch, &
          & z_mc=nh_metrics%z_mc(:,:,:), &
          & z_ifc=nh_metrics%z_ifc(:,:,:), &
          & tracer=nh_prog_rcf%tracer(:,:,:,:), &
          & mixing_length=mixing_length(:,:,:), &
          & total_turbulence_energy=mem%total_turbulence_energy(:,:,:), &
          & exchange_coeff_h=mem%exchange_coeff_h(:,:,:), &
          & total_water_var=mem%total_water_var(:,:,:), &
          & rcld=phy_diag%rcld(:,:,:) &
        )
    END IF

    !$ACC WAIT(1)
    !$ACC EXIT DATA &
    !$ACC   DETACH(ocean_u, ocean_v, p_graupel_gsp_rate, p_ice_gsp_rate, p_hail_gsp_rate)

    IF (is_coupled_to_ocean() .AND. .NOT. linit) THEN
      CALL sea_model_couple_ocean (&
          & patch=patch, &
          & list_sea=ext_data%atm%list_sea, &
          & fr_sft=fr_sft(:,:,:), &
          & alb=alb, &
          & flx_rad=flx_rad, &
          & pres_sfc=nh_diag%pres_sfc(:,:), &
          & t_eff_sft=t_eff_sft(:,:,:), &
          & evapo_sft=evapo_sft(:,:,:), &
          & flx_heat_latent_sft=flx_heat_latent_sft(:,:,:), &
          & flx_heat_sensible_sft=flx_heat_sensible_sft(:,:,:), &
          & condhf_ice=diag_lnd%condhf_ice(:,:), &
          & meltpot_ice=diag_lnd%meltpot_ice(:,:), &
          & co2_concentration_srf=co2_concentration_srf(:,:), &
          & rain_srf=rain_srf(:,:), &
          & snow_srf=snow_srf(:,:), &
          & umfl_sft=flx_mom_u_sft(:,:,:), &
          & vmfl_sft=flx_mom_v_sft(:,:,:), &
          & sp_10m=phy_diag%sp_10m(:,:), &
          & t_seasfc=diag_lnd%t_seasfc(:,:), &
          & fr_seaice=diag_lnd%fr_seaice(:,:), &
          & h_ice=prog_wtr_new%h_ice(:,:), &
          & sea_state=mem%sea_state &
        )
    END IF

    !$ACC WAIT(1)
    !$ACC EXIT DATA &
    !$ACC   DELETE(LIST_CREATE1) &
    !$ACC   DELETE(LIST_CREATE2) &
    !$ACC   DELETE(LIST_CREATE3) &
    !$ACC   DELETE(LIST_CREATE4) &
    !$ACC   DELETE(LIST_CREATE5)
#   undef LIST_CREATE5
#   undef LIST_CREATE4
#   undef LIST_CREATE3
#   undef LIST_CREATE2
#   undef LIST_CREATE1

  END SUBROUTINE nwp_vdiff


  !> Perform early setup tasks. This creates the variable lists for jsbach.
  SUBROUTINE nwp_vdiff_setup (patch, use_jsbach)

    TYPE(t_patch), INTENT(IN) :: patch(:)
    LOGICAL,       INTENT(IN) :: use_jsbach(:)

    INTEGER :: jg

    CALL vdiff_init(2, ntracer)

#ifndef __NO_JSBACH__
    DO jg = 1, SIZE(patch)
      IF (use_jsbach(jg)) CALL jsbach_init(patch(jg)%id)
    ENDDO
#endif

    ! TODO: This does double work, but we need some of the AES routines in VDIFF.
    CALL init_convect_tables
    CALL init_aes_convect_tables

    CALL setup_jsbach_init_vars

  END SUBROUTINE


  SUBROUTINE nwp_vdiff_init ( &
        & patch, vdiff_config, phy_diag, ext_data, lnd_prog_now, lnd_diag, vdiff_state &
      )

    TYPE(t_patch), INTENT(IN) :: patch !< Patch to initialize.
    TYPE(t_vdiff_config), INTENT(IN) :: vdiff_config !< VDIFF configuration.
    TYPE(t_nwp_phy_diag), INTENT(IN) :: phy_diag !< Physics variables
    TYPE(t_external_data), INTENT(IN) :: ext_data !< External parameters.
    TYPE(t_lnd_prog), INTENT(IN) :: lnd_prog_now !< Land prognostic variables.
    TYPE(t_lnd_diag), INTENT(IN) :: lnd_diag !< Land diagnostic variables.
    TYPE(t_nwp_vdiff_state), INTENT(INOUT) :: vdiff_state !< vdiff+jsbach state structure.

    IF (isRestart()) RETURN

    CALL message(module_name // ':nwp_vdiff_init', 'Cold-Initializing vdiff state')

    CALL sea_model_init( &
        & patch, ext_data, phy_diag%cosmu0(:,:), lnd_diag%t_seasfc(:,:), &
        & ext_data%atm_td%sst_m(:,:,:), vdiff_state%sea_state)

    vdiff_state%fact_q_air(:,:) = 0.5_wp
    vdiff_state%fact_qsat_srf(:,:) = 0.5_wp

    vdiff_state%flx_co2_natural_land(:,:) = 0._wp
    vdiff_state%flx_heat_latent_sft(:,:,:) = 0._wp
    vdiff_state%flx_heat_sensible_sft(:,:,:) = 0._wp
    vdiff_state%fr_ice_on_lake(:,:) = MERGE(1._wp, 0._wp, lnd_prog_now%t_g(:,:) < tf_fresh)
    vdiff_state%temp_sft(:,:,:) = SPREAD(lnd_prog_now%t_g(:,:), 3, SIZE(vdiff_state%temp_sft, 3))
    vdiff_state%theta_v_var(:,:,:) = 0._wp
    vdiff_state%total_turbulence_energy(:,:,:) = 1.e-4_wp
    vdiff_state%total_water_var(:,:,:) = 0._wp
    vdiff_state%ustar(:,:) = 1._wp
    vdiff_state%wstar_sfc(:,:,:) = 0._wp
    vdiff_state%z0h_land(:,:) = MAX(1e-3_wp, ext_data%atm%z0(:,:))
    vdiff_state%z0m_sfc(:,:,SFC_LAND) = MAX(1e-3_wp, ext_data%atm%z0(:,:))
    vdiff_state%z0m_sfc(:,:,SFC_WATER) = MAX(1e-3_wp, ext_data%atm%z0(:,:))
    vdiff_state%z0m_sfc(:,:,SFC_ICE) = vdiff_config%z0m_ice

    CALL vdiff_state%h2d()

  END SUBROUTINE nwp_vdiff_init


  !>
  !! Compute the fraction of each surface class recognized by vdiff.  Surface classes are land,
  !! open water, and ice. This routine needs lake and land area fractions for each grid box, as
  !! well as the fraction of sea and lake area covered by ice. Land and lake area fractions are
  !! obtained from external data in `ext_data`. The lake ice fraction is computed by JSBACH. Sea
  !! ice has to be provided by a sea ice scheme. Both are stored in `mem`.
  SUBROUTINE get_surface_type_fractions (patch, ext_data, mem, lnd_diag, fr_sfc, fr_sft)

    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.
    TYPE(t_external_data), INTENT(IN) :: ext_data !< External data for the patch.
    TYPE(t_nwp_vdiff_state), INTENT(IN) :: mem !< vdiff and jsbach state.
    TYPE(t_lnd_diag), INTENT(IN) :: lnd_diag
    !> Fraction of each surface class (nproma,nblks_c,SFC_NUM).
    REAL(wp), INTENT(OUT) :: fr_sfc(:,:,:)
    !> Fraction of each surface type (nproma,nblks_c,SFT_NUM).
    REAL(wp), INTENT(OUT) :: fr_sft(:,:,:)

    REAL(wp) :: fr_sea !< Sea fraction.
    REAL(wp) :: fr_total !< Total of fractions after correction.

    INTEGER :: i_startblk, i_endblk
    INTEGER :: ics, ice
    INTEGER :: ic, i_blk

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    !$OMP PARALLEL
    !$OMP DO PRIVATE(i_blk, ics, ice, ic, fr_sea, fr_total)
    DO i_blk = i_startblk, i_endblk
      CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
          & end_prog_cells)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR PRIVATE(fr_sea, fr_total)
        DO ic = ics, ice
          fr_sea = 1._wp - ext_data%atm%fr_land(ic,i_blk) - ext_data%atm%fr_lake(ic,i_blk)
          fr_total = 1._wp

          IF (fr_sea < frsea_thrhld) THEN
            fr_sea = 0._wp
            fr_total = ext_data%atm%fr_land(ic,i_blk) + ext_data%atm%fr_lake(ic,i_blk)
          END IF

          fr_sft(ic,i_blk,SFT_LAND) = ext_data%atm%fr_land(ic,i_blk) / fr_total
          fr_sft(ic,i_blk,SFT_LWTR) = ext_data%atm%fr_lake(ic,i_blk) / fr_total &
              & * (1._wp - mem%fr_ice_on_lake(ic,i_blk))
          fr_sft(ic,i_blk,SFT_SWTR) = fr_sea  / fr_total * (1._wp - lnd_diag%fr_seaice(ic,i_blk))
          fr_sft(ic,i_blk,SFT_LICE) = ext_data%atm%fr_lake(ic,i_blk) / fr_total &
              & * mem%fr_ice_on_lake(ic,i_blk)
          fr_sft(ic,i_blk,SFT_SICE) = fr_sea / fr_total * lnd_diag%fr_seaice(ic,i_blk)

          fr_sfc(ic,i_blk,SFC_LAND) = fr_sft(ic,i_blk,SFT_LAND)
          fr_sfc(ic,i_blk,SFC_WATER) = fr_sft(ic,i_blk,SFT_LWTR) + fr_sft(ic,i_blk,SFT_SWTR)
          fr_sfc(ic,i_blk,SFC_ICE) = fr_sft(ic,i_blk,SFT_LICE) + fr_sft(ic,i_blk,SFT_SICE)
        END DO
      !$ACC END PARALLEL
    END DO
    !$OMP END PARALLEL

  END SUBROUTINE get_surface_type_fractions


  !>
  !! Compute the temperature of each surface class recognized by vdiff. Surface classes are land,
  !! open water, and ice. This routine needs lake and land area fractions and temperatures for
  !! each grid box, as well as the fraction and temperature of sea and lake area covered by ice.
  !!
  !! \note The function takes water and ice temperatures as area-weighted averages of sea and lake
  !! values. Some kind of flux averaging would probably be more appropriate.
  !!
  SUBROUTINE get_surface_class_temperature (patch, fr_sft, temp_sft, fr_sfc, temp_sfc)

    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.
    !> Fraction of each surface type [1] (nproma,nblks_c,SFT_NUM).
    REAL(wp), INTENT(IN) :: fr_sft(:,:,:)
    !> Temperature of each surface type [K] (nproma,nblks_c,SFT_NUM).
    REAL(wp), INTENT(IN) :: temp_sft(:,:,:)
    !> Fraction of each surface class [1] (nproma,nblks_c,SFC_NUM).
    REAL(wp), INTENT(IN) :: fr_sfc(:,:,:)
    !> Temperature of each surface class [K] (nproma,nblks_c,SFC_NUM).
    REAL(wp), INTENT(OUT) :: temp_sfc(:,:,:)

    INTEGER :: i_startblk, i_endblk
    INTEGER :: ics, ice
    INTEGER :: ic, i_blk

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
          temp_sfc(ic,i_blk,SFC_LAND) = temp_sft(ic,i_blk,SFT_LAND)

          IF (fr_sfc(ic,i_blk,SFC_WATER) > 0._wp) THEN
            temp_sfc(ic,i_blk,SFC_WATER) = ( &
                &   temp_sft(ic,i_blk,SFT_LWTR) * fr_sft(ic,i_blk,SFT_LWTR) &
                &   + temp_sft(ic,i_blk,SFT_SWTR) * fr_sft(ic,i_blk,SFT_SWTR) &
                & ) / fr_sfc(ic,i_blk,SFC_WATER)
          ELSE
            ! Dummy value, should not have any significance.
            temp_sfc(ic,i_blk,SFC_WATER) = 0._wp
          END IF

          IF (fr_sfc(ic,i_blk,SFC_ICE) > 0._wp) THEN
            temp_sfc(ic,i_blk,SFC_ICE) = ( &
                &   temp_sft(ic,i_blk,SFT_LICE) * fr_sft(ic,i_blk,SFT_LICE) &
                &   + temp_sft(ic,i_blk,SFT_SICE) * fr_sft(ic,i_blk,SFT_SICE) &
                & ) / fr_sfc(ic,i_blk,SFC_ICE)
          ELSE
            ! Dummy value, should not have any significance.
            temp_sfc(ic,i_blk,SFC_ICE) = 0._wp
          END IF
        END DO
      !$ACC END PARALLEL
    END DO
    !$OMP END PARALLEL

  END SUBROUTINE get_surface_class_temperature


  !>
  !! Collect the surface CO2 concentration according to the carbon-cycle configuration.
  !! Time-dependent values are retrieved from `mo_bc_greenhouse_gases::ghg_co2mmr`.
  SUBROUTINE get_surface_co2_concentration (patch, ccycle_config, tracer, co2_concentration_srf)

    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.
    !> Carbon-cycle configuration.
    TYPE(t_ccycle_config), INTENT(IN) :: ccycle_config
    !> Tracer concentrations at current time step [kg/kg].
    REAL(wp), INTENT(IN) :: tracer(:,:,:,:)
    !> CO2 concentration at surface [kg/kg] (nproma, nblks_c).
    REAL(wp), INTENT(OUT) :: co2_concentration_srf(:,:)

    !> Routine name.
    CHARACTER(len=*), PARAMETER :: routine_name = &
        & TRIM(module_name) // ':get_surface_co2_concentration'

    !> CO2 volume mixing ration for 1990 [mol/mol].
    REAL(wp), PARAMETER :: CO2VMR_1990 = 348.0e-06_wp

    INTEGER :: i_startblk, i_endblk
    INTEGER :: ics, ice
    INTEGER :: ic, i_blk

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    SELECT CASE (ccycle_config%iccycle)
    CASE(CCYCLE_MODE_NONE)
      ! Should not make a difference, but according to a comment in the ECHAM interface, jsbach
      ! results change if this is changed.
      !$OMP PARALLEL
      !$OMP DO PRIVATE(i_blk, ics, ice, ic)
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO ic = ics, ice
            co2_concentration_srf(ic,i_blk) = CO2VMR_1990 * vmr_to_mmr_co2
          END DO
        !$ACC END PARALLEL
      END DO
      !$OMP END PARALLEL

    CASE(CCYCLE_MODE_INTERACTIVE)
      IF (ico2 <= 0) THEN
        CALL finish(TRIM(routine_name), 'Interactive carbon cycle activated but no CO2 tracer.')
      END IF

      !$OMP PARALLEL
      !$OMP DO PRIVATE(i_blk, ics, ice, ic)
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR
          DO ic = ics, ice
            co2_concentration_srf(ic,i_blk) = tracer(ic,patch%nlev,i_blk,ico2)
          END DO
        !$ACC END PARALLEL
      END DO
      !$OMP END PARALLEL

    CASE(CCYCLE_MODE_PRESCRIBED)
      SELECT CASE(ccycle_config%ico2conc)
      CASE(CCYCLE_CO2CONC_CONST)
        ! Constant concentration throughout the simulation.

        !$OMP PARALLEL
        !$OMP DO PRIVATE(i_blk, ics, ice, ic)
        DO i_blk = i_startblk, i_endblk
          CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
              & end_prog_cells)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR
            DO ic = ics, ice
              co2_concentration_srf(ic,i_blk) = ccycle_config%vmr_co2 * vmr_to_mmr_co2
            END DO
          !$ACC END PARALLEL
        END DO
        !$OMP END PARALLEL

      CASE(CCYCLE_CO2CONC_FROMFILE)
        ! Time-dependent concentration (location independent).

        !$OMP PARALLEL
        !$OMP DO PRIVATE(i_blk, ics, ice, ic)
        DO i_blk = i_startblk, i_endblk
          CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
              & end_prog_cells)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
            !$ACC LOOP GANG VECTOR
            DO ic = ics, ice
              co2_concentration_srf(ic,i_blk) = ghg_co2mmr
            END DO
          !$ACC END PARALLEL
        END DO
        !$OMP END PARALLEL

      END SELECT
    END SELECT

  END SUBROUTINE get_surface_co2_concentration


  !>
  !! Calculate a near-surface temperature from a dry static energy.
  !! Lifted from ECHAM physics' nsurf_diag, implements
  !! Geleyn, Tellus A 40, 347-351 (1988), doi: 10.3402/tellusa.v40i4.11805
  SUBROUTINE get_near_surface_temperature ( &
        & ics, ice, z_ref, z_lowest, z_srf, ri_number, b_neutral, b_h_actual, s_lowest, s_srf, &
        & temp_ref &
      )

    INTEGER, INTENT(IN) :: ics !< Starting cell index.
    INTEGER, INTENT(IN) :: ice !< End cell index.
    REAL(wp), INTENT(IN) :: z_ref !< Reference height above ground [m].
    REAL(wp), INTENT(IN) :: z_lowest(:) !< Height of the lowest level [m] (ics:ice).
    REAL(wp), INTENT(IN) :: z_srf(:) !< Height of the surface [m] (ics:ice).
    REAL(wp), INTENT(IN) :: ri_number(:) !< Richardson number [1] (ics:ice).
    REAL(wp), INTENT(IN) :: b_neutral(:) !< Neutral coefficient [1] (ics:ice).
    REAL(wp), INTENT(IN) :: b_h_actual(:) !< Coefficient for actual conditions [1] (ics:ice)
    REAL(wp), INTENT(IN) :: s_lowest(:) !< Dry static energy at lowest level [J/kg] (ics:ice).
    REAL(wp), INTENT(IN) :: s_srf(:) !< Dry static energy at surface [J/kg] (ics:ice).
    REAL(wp), INTENT(INOUT) :: temp_ref(:) !< Temperature at reference height [K] (ics:ice).

    REAL(wp) :: zrat
    REAL(wp) :: cbn
    REAL(wp) :: cbs
    REAL(wp) :: cbu
    REAL(wp) :: zred
    REAL(wp) :: s_ref
    REAL(wp) :: b_h

    INTEGER :: ic

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(zrat, cbn, cbs, cbu, zred, s_ref, b_h)
      DO ic = ics, ice
        IF (b_h_actual(ic) == 0._wp) THEN
          ! No heat flux from surface, probably because tile fraction is zero.
          ! No flux implies no gradient between surface and lowest layer, so either value can be
          ! used in principle. However, because tile fraction might be zero, the tile surface
          ! temperature may not be well-defined. Therefore, we take dry static energy in lowest
          ! layer, which should be sensible because it is tile independent.
          s_ref = s_lowest(ic)
        ELSE
          b_h = b_h_actual(ic)

          zrat = z_ref / (z_lowest(ic) - z_srf(ic))
          cbn = LOG(1._wp + (EXP(b_neutral(ic)) - 1._wp) * zrat )
          cbs = -(b_neutral(ic) - b_h) * zrat
          cbu = -LOG(1._wp + (EXP(b_neutral(ic) - b_h) - 1._wp) * zrat)

          zred = (cbn + MERGE(cbs, cbu, ri_number(ic) > 0._wp)) / b_h
          s_ref = s_srf(ic) + zred * (s_lowest(ic) - s_srf(ic))
        END IF

        temp_ref(ic) = (s_ref - z_ref*grav) / cpd
      END DO
    !$ACC END PARALLEL

  END SUBROUTINE get_near_surface_temperature


  !>
  !! Calculate the dew point at a reference height above ground. Requires temperature at that
  !! height.
  !!
  !! This routine inverts `spec_humi`, `sat_pres_water` and `sat_pres_ice` from `mo_satad` to
  !! compute the dew point from the specific humidity at the reference height. Specific humidity
  !! is computed by assuming that relative humidity is constant throughout the lowest atmospheric
  !! layer. The routine calculates pressure at the reference height, diagnoses the saturation
  !! specific humidity at that pressure and the temperature at reference height, and uses the
  !! relative humidity to get the specific humidity. It then runs that through the inverse of
  !! `spec_humi` to compute the vapor pressure, and through the inverse of `sat_pres_water` to
  !! compute the dew point. The output values are clamped assuming that there is no oversaturation
  !! present. This is true after saturation adjustment but the temperatures will be different, so
  !! the outputs are an approximation in case of oversaturated air.
  !!
  !! Nomenclature is that from `mo_satad`, except that `b3 = tf_fresh` (the melting point of
  !! water).
  SUBROUTINE get_near_surface_dew_point ( &
        & ics, ice, z_ref, temp_ref, temp_lowest, pres_lowest, pres_srf, q_lowest, x_lowest, &
        & tdew_ref, rh_ref, q_ref &
      )

    USE mo_lookup_tables_constants, ONLY: b1 => c1es, b2w => c3les, b4w => c4les

    INTEGER, INTENT(IN) :: ics !< Starting cell index.
    INTEGER, INTENT(IN) :: ice !< End cell index.
    REAL(wp), INTENT(IN) :: z_ref !< Reference height above ground [m].
    REAL(wp), INTENT(IN) :: temp_ref(:) !< Temperature at reference height [K] (ics:ice).
    REAL(wp), INTENT(IN) :: temp_lowest(:) !< Temperature at lowest level [K] (ics:ice).
    REAL(wp), INTENT(IN) :: pres_lowest(:) !< Pressure at lowest level [Pa] (ics:ice).
    REAL(wp), INTENT(IN) :: pres_srf(:) !< Pressure at surface [Pa] (ics:ice).
    REAL(wp), INTENT(IN) :: q_lowest(:) !< Specific humidity at lowest level [kg/kg] (ics:ice).
    REAL(wp), INTENT(IN) :: x_lowest(:) !< Total cloud water at lowest level [kg/kg] (ics:ice).
    REAL(wp), INTENT(INOUT) :: tdew_ref(:) !< Dew point at reference height [K] (ics:ice).

    !> Relative humidity at reference height [%] (ics:ice).
    REAL(wp), OPTIONAL, INTENT(INOUT) :: rh_ref(:)
    !> Specific humidity at reference height [kg/kg] (ics:ice).
    REAL(wp), OPTIONAL, INTENT(INOUT) :: q_ref(:)

    ! Constants
    REAL(wp), PARAMETER :: rh_min = 0.05_wp !< Minimum relative humidity.

    ! Locals
    INTEGER :: ic
    REAL(wp) :: qsat_lowest !< Saturation specific humidity in lowest layer [kg/kg].
    REAL(wp) :: qsat_ref !< Saturation specific humidity at reference height [kg/kg].
    REAL(wp) :: q !< Specific humidity at reference height [kg/kg].

    REAL(wp) :: rh !< Relative humidity in lowest layer [1].
    REAL(wp) :: pres_ref !< Pressure at reference height [Pa].

    REAL(wp) :: logfactor !< log(p_sat/b1)/b2

    LOGICAL :: have_rh_ref
    LOGICAL :: have_q_ref

    have_rh_ref = PRESENT(rh_ref)
    have_q_ref = PRESENT(q_ref)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(qsat_lowest, pres_ref, qsat_ref, q, logfactor)
      DO ic = ics, ice
        qsat_lowest = spec_humi(sat_pres_water(temp_lowest(ic)), pres_lowest(ic))

        rh  = MAX(rh_min, q_lowest(ic) / qsat_lowest)

        ! RW: I guess this is the leading term for the z-dependence of the pressure for the
        ! isothermal atmosphere, p(z_ref) = p(0) * exp(-g z_ref/(Rd Tv)). If one assumes
        ! linearity anyways, wouldn't it be more appropriate to do linear interpolation?
        pres_ref = pres_srf(ic) * &
            & (1._wp - z_ref*grav / &
            &   ( rd * temp_ref(ic) * (1._wp + vtmpc1 * q_lowest(ic) - x_lowest(ic))))

        qsat_ref = spec_humi(sat_pres_water(temp_ref(ic)), pres_ref)

        q = rh * qsat_ref

        logfactor = LOG(q * pres_ref / (rdv * b1 * (1._wp + vtmpc1 * q))) / b2w
        tdew_ref(ic) = MIN(temp_ref(ic), (tf_fresh - b4w * logfactor) / (1._wp - logfactor))

        IF (have_rh_ref) rh_ref(ic) = MIN(rh, 1._wp)
        IF (have_q_ref) q_ref(ic) = MIN(q, qsat_ref)
      END DO
    !$ACC END PARALLEL

  END SUBROUTINE get_near_surface_dew_point


  !>
  !! Calculate wind velocity and speed at a reference height above the surface.
  !! Lifted from ECHAM physics' nsurf_diag, implements
  !! Geleyn, Tellus A 40, 347-351 (1988), doi: 10.3402/tellusa.v40i4.11805
  SUBROUTINE get_near_surface_wind ( &
        & ics, ice, z_ref, z_lowest, z_srf, u_lowest, v_lowest, ocean_u, ocean_v, ri_number, &
        & b_neutral, b_m_actual, u_ref, v_ref, sp_ref &
      )

    INTEGER, INTENT(IN) :: ics !< Starting start index.
    INTEGER, INTENT(IN) :: ice !< Cell end index.
    REAL(wp), INTENT(IN) :: z_ref !< Reference height above ground [m].
    REAL(wp), INTENT(IN) :: z_lowest(:) !< Height of the lowest level [m] (ics:ice).
    REAL(wp), INTENT(IN) :: z_srf(:) !< Height of the surface [m] (ics:ice).
    REAL(wp), INTENT(IN) :: u_lowest(:) !< Zonal wind at lowest level [m/s] (ics:ice).
    REAL(wp), INTENT(IN) :: v_lowest(:) !< Meridional wind at lowest level [m/s] (ics:ice).
    REAL(wp), INTENT(IN) :: ocean_u(:) !< Zonal ocean surface velocity [m/s] (ics:ice).
    REAL(wp), INTENT(IN) :: ocean_v(:) !< Meridional ocean surface velocity [m/s] (ics:ice).
    REAL(wp), INTENT(IN) :: ri_number(:) !< Richardson number [1] (ics:ice).
    REAL(wp), INTENT(IN) :: b_neutral(:) !< Neutral coefficient [1] (ics:ice).
    REAL(wp), INTENT(IN) :: b_m_actual(:) !< Coefficient for actual conditions [1] (ics:ice)
    REAL(wp), INTENT(INOUT) :: u_ref(:) !< Zonal wind at reference height [m/s] (ics:ice).
    REAL(wp), INTENT(INOUT) :: v_ref(:) !< Meridional wind at reference height [m/s] (ics:ice).
    REAL(wp), INTENT(INOUT) :: sp_ref(:) !< Wind speed at reference height [m/s] (ics:ice).

    REAL(wp) :: zrat
    REAL(wp) :: cbn
    REAL(wp) :: cbs
    REAL(wp) :: cbu
    REAL(wp) :: zred
    REAL(wp) :: b_m

    INTEGER :: ic

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(zrat, cbn, cbs, cbu, zred, b_m)
      DO ic = ics, ice
        IF (b_m_actual(ic) == 0._wp) THEN
          ! No momentum flux from surface, probably because tile fraction is zero.
          ! Take wind in lowest layer.
          u_ref(ic) = u_lowest(ic)
          v_ref(ic) = v_lowest(ic)
        ELSE
          b_m = b_m_actual(ic)

          zrat = z_ref / (z_lowest(ic) - z_srf(ic))
          cbn = LOG(1._wp + (EXP(b_neutral(ic)) - 1._wp) * zrat)
          cbs = -(b_neutral(ic) - b_m) * zrat
          cbu = -LOG(1._wp + (EXP(b_neutral(ic) - b_m) - 1._wp) * zrat)
          zred = (cbn + MERGE(cbs, cbu, ri_number(ic) > 0._wp)) / b_m
          u_ref(ic) = zred * (u_lowest(ic) - ocean_u(ic)) + ocean_u(ic)
          v_ref(ic) = zred * (v_lowest(ic) - ocean_v(ic)) + ocean_v(ic)
        END IF

        sp_ref(ic) = SQRT(u_ref(ic)**2 + v_ref(ic)**2)
      END DO
    !$ACC END PARALLEL

  END SUBROUTINE get_near_surface_wind


  !>
  !! Calculate the gust speed at a reference level given the mean wind speed at the same level and
  !! the friction velocity. Compare `nwp_dyn_gust` in `mo_util_phys`. This routine is OpenMP
  !! orphaned.
  SUBROUTINE get_gust_speed (patch, phy_diag, nh_diag, nh_metrics, sp_ref, ustar, gust_ref)

    TYPE(t_patch), INTENT(IN) :: patch ! Current patch.
    TYPE(t_nwp_phy_diag), INTENT(IN) :: phy_diag !< Diagnostic physics fields.
    TYPE(t_nh_diag), INTENT(IN) :: nh_diag !< Diagnostic atmosphere fields.
    TYPE(t_nh_metrics), INTENT(IN) :: nh_metrics !< Metrics of the grid.
    REAL(wp), INTENT(IN) :: sp_ref(:,:) !< Wind speed at the reference level [m/s].
    REAL(wp), INTENT(IN) :: ustar(:,:) !< Friction velocity [m/s]
    REAL(wp), INTENT(INOUT) :: gust_ref(:,:) !< Gust speed at the reference level.

    !> Windspeed threshold above which gusts become stronger [m/s].
    REAL(wp), PARAMETER :: windspeed_threshold = 10._wp
    !> Coefficient determining gust increase if windspeed above threshold [ustar/(m/s)].
    !! Gets doubled in mountainous cells.
    REAL(wp), PARAMETER :: windspeed_coefficient = 0.2_wp
    !> Maximum for windspeed factor [ustar]. Gets doubled in mountainous cells.
    REAL(wp), PARAMETER :: windspeed_factor_max = 2._wp
    !> Coefficient for gust increase in mountainous cells [ustar].
    REAL(wp), PARAMETER :: mountain_factor = 2._wp

    INTEGER :: k_top !< Level index of top of SSO envelope layer.
    REAL(wp) :: sp_env !< Wind speed at the top of SSO envelope layer [m/s].
    REAL(wp) :: sp_lowest !< Wind speed at lowest model level [m/s].
    REAL(wp) :: uadd_sso !< Correction due to SSO [m/s].
    REAL(wp) :: windspeed_factor !< Additional gust coefficient [1].
    REAL(wp) :: mtnmask !< Mountain point mask.

    INTEGER :: i_startblk
    INTEGER :: i_endblk
    INTEGER :: i_blk
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: ic

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    !$OMP DO
    DO i_blk = i_startblk, i_endblk
      CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
          & end_prog_cells)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR PRIVATE(k_top, mtnmask, sp_env, sp_lowest, uadd_sso, windspeed_factor)
        DO ic = ics, ice
          IF (phy_diag%ktop_envel(ic,i_blk) < patch%nlev &
              & .AND. phy_diag%ktop_envel(ic,i_blk) > 1) THEN
            k_top = phy_diag%ktop_envel(ic,i_blk) - 1
          ELSE
            k_top = patch%nlev
          END IF

          mtnmask = nh_metrics%mask_mtnpoints_g(ic,i_blk)

          sp_env = SQRT(nh_diag%u(ic,k_top,i_blk)**2 + nh_diag%v(ic,k_top,i_blk)**2)
          sp_lowest = SQRT(nh_diag%u(ic,patch%nlev,i_blk)**2 + nh_diag%v(ic,patch%nlev,i_blk)**2)

          uadd_sso = MAX(0._wp, sp_env - sp_lowest)
          windspeed_factor = (1._wp + mtnmask)  &
              & * MAX(0._wp, MIN(windspeed_factor_max, &
              &     windspeed_coefficient*(sp_ref(ic,i_blk)-windspeed_threshold)))
          gust_ref(ic,i_blk) = sp_ref(ic,i_blk) + mtnmask*uadd_sso &
              & + (tune_gust_factor + windspeed_factor + mountain_factor*mtnmask) * ustar(ic,i_blk)
        END DO
      !$ACC END PARALLEL
    END DO

  END SUBROUTINE get_gust_speed


  !>
  !! Compute the standard deviation of the saturation deficit and total water variance.
  !!
  !! We use a simple diagnostic formula to get the total water variance, without advection or
  !! diffusion. The dissipation timescale is set to \f$ \tau = l_h/\sqrt{E} \f$, with mixing length
  !! \f$ l_h \f$ and specific total turbulent energy \f$ E \f$. It would be desirable for the
  !! timescale to be proportional to the horizontal resolution so the variance becomes smaller
  !! when resolution is increased.
  !!
  !! vdiff defines its total water variance on model levels, unlike TTE or virtual temperature
  !! variance, which are defined on interfaces. For consistency, we compute the water variance on
  !! the interior interfaces where we can make finite-difference approximations to dqt/dz. For the
  !! boundary interfaces we adopt the value of the closest interior interface. Finally, we
  !! interpolate the variance back to the model levels.
  !!
  !! In the standard deviation of the saturation deficit, temperature variations are neglected.
  !! Thus, it is just the square root of the total water variance, evaluated on the interfaces.
  SUBROUTINE get_stddev_saturation_deficit ( &
        & patch, z_mc, z_ifc, tracer, mixing_length, total_turbulence_energy, exchange_coeff_h, &
        & total_water_var, rcld &
      )

    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.
    REAL(wp), INTENT(IN) :: z_mc(:,:,:) !< Z coordinate of model level centers.
    REAL(wp), INTENT(IN) :: z_ifc(:,:,:) !< Z coordinate of level interfaces.
    REAL(wp), INTENT(IN) :: tracer(:,:,:,:) !< Tracer array [X/kg].
    REAL(wp), INTENT(IN) :: mixing_length(:,:,:) !< Turbulent mixing length [m].
    REAL(wp), INTENT(IN) :: total_turbulence_energy(:,:,:) !< Total turbulence energy [J/kg].
    REAL(wp), INTENT(IN) :: exchange_coeff_h(:,:,:)
    !< Turbulent diffusion coefficient for heat [m**2/s].
    REAL(wp), INTENT(INOUT) :: total_water_var(:,:,:) !< Total water variance [(kg/kg)**2].
    REAL(wp), INTENT(INOUT) :: rcld(:,:,:) !< Standard deviation of saturation deficit [kg/kg].

    !> Total water variance on interfaces.
    REAL(wp) :: qtvar(nproma, patch%nlev+1)

    !> Product of dissipation timescale and diffusion coefficient [m**2].
    REAL(wp) :: tau_k
    !> Difference in total water between levels [kg/kg].
    REAL(wp) :: delta_qt
    !> Height difference between levels [m].
    REAL(wp) :: delta_z
    !> Interpolation weight.
    REAL(wp) :: wg

    INTEGER :: i_startblk, i_endblk, i_blk
    INTEGER :: ics, ice, ic
    INTEGER :: kl

    !$ACC ENTER DATA CREATE(qtvar) ASYNC(1)

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    !$OMP PARALLEL
    !$OMP DO PRIVATE(i_blk, ics, ice, ic, kl, tau_k, delta_qt, delta_z, wg, qtvar)
    DO i_blk = i_startblk, i_endblk
      CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
          & end_prog_cells)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP SEQ
        DO kl = 2, patch%nlev
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(tau_k, delta_qt, delta_z)
          DO ic = ics, ice
            tau_k = mixing_length(ic,kl-1,i_blk) &
                & / SQRT(total_turbulence_energy(ic,kl-1,i_blk)) &
                & * exchange_coeff_h(ic,kl,i_blk)

            delta_qt = &
                & + tracer(ic,kl-1,i_blk,iqv) - tracer(ic,kl,i_blk,iqv) &
                & + tracer(ic,kl-1,i_blk,iqc) - tracer(ic,kl,i_blk,iqc) &
                & + tracer(ic,kl-1,i_blk,iqi) - tracer(ic,kl,i_blk,iqi)

            delta_z = z_mc(ic,kl-1,i_blk) - z_mc(ic,kl,i_blk)

            qtvar(ic,kl) = tau_k * (delta_qt / delta_z)**2
          END DO
        END DO

        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = ics, ice
          qtvar(ic,1) = qtvar(ic,2)
          qtvar(ic,patch%nlev+1) = qtvar(ic,patch%nlev)
        END DO

        !$ACC LOOP SEQ
        DO kl = 1, patch%nlev
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wg)
          DO ic = ics, ice
            rcld(ic,kl,i_blk) = SQRT(qtvar(ic,kl))

            wg = (z_mc(ic,kl,i_blk) - z_ifc(ic,kl+1,i_blk)) &
                & / ((z_ifc(ic,kl,i_blk) - z_ifc(ic,kl+1,i_blk)))

            total_water_var(ic,kl,i_blk) = wg * qtvar(ic,kl) + (1._wp - wg) * qtvar(ic,kl+1)
          END DO
        END DO
      !$ACC END PARALLEL
    END DO
    !$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(qtvar)

  END SUBROUTINE get_stddev_saturation_deficit


  !> Compute the surface stress (momentum flux) over each surface type.
  !!
  !! This routine accounts for movement of the ocean surface.
  !! Has to be called once for each block of cells.
  SUBROUTINE get_surface_stress ( &
        & ics, ice, delta_time, prefactor_exchange, exchange_coeff_m_sfc, uv_acoef, u_bcoef, &
        & v_bcoef, ocean_u, ocean_v, zero, umfl_sft, vmfl_sft &
      )

    INTEGER, INTENT(IN) :: ics !< Start cell index.
    INTEGER, INTENT(IN) :: ice !< End cell index.
    REAL(wp), INTENT(IN) :: delta_time !< Time step [s].

    !> Common prefactor for surface exchange coefficients [kg*s/m**3] (ics:ice).
    REAL(wp), INTENT(IN) :: prefactor_exchange(:)
    !> Exchange coefficient for momentum over surface classes [m/s] (ics:ice,SFC_NUM).
    REAL(wp), INTENT(IN) :: exchange_coeff_m_sfc(:,:)
    !> Richtmyer E for momentum over surface classes [1] (ics:ice,SFC_NUM).
    REAL(wp), INTENT(IN) :: uv_acoef(:,:)
    !> Richtmyer F for zonal momentum over surface classes [m/s] (ics:ice,SFC_NUM).
    REAL(wp), INTENT(IN) :: u_bcoef(:,:)
    !> Richtmyer F for meridional momentum over surface classes [m/s] (ics:ice,SFC_NUM).
    REAL(wp), INTENT(IN) :: v_bcoef(:,:)
    !> Zonal ocean velocity [m/s] (ics:ice).
    REAL(wp), TARGET, CONTIGUOUS, INTENT(IN) :: ocean_u(:)
    !> Meridional ocean velocity [m/s] (ics:ice).
    REAL(wp), TARGET, CONTIGUOUS, INTENT(IN) :: ocean_v(:)
    !> Zero field (for surface u and v over land, etc.) [1] (ics:ice).
    REAL(wp), TARGET, CONTIGUOUS, INTENT(IN) :: zero(:)

    !> Zonal momentum flux over surface types [N/m**2] (ics:ice,SFT_NUM).
    REAL(wp), INTENT(INOUT) :: umfl_sft(:,:)
    !> Meridional momentum flux over surface types [N/m**2] (ics:ice,SFT_NUM).
    REAL(wp), INTENT(INOUT) :: vmfl_sft(:,:)

    INTEGER :: isfc, isft

    REAL(wp), POINTER, CONTIGUOUS :: u_sfc(:)
    REAL(wp), POINTER, CONTIGUOUS :: v_sfc(:)

    DO isft = 1, SFT_NUM
      isfc = SFT_CLASS(isft)

      SELECT CASE (isft)
      CASE (SFT_SWTR)
        u_sfc => ocean_u
        v_sfc => ocean_v
      CASE DEFAULT
        u_sfc => zero
        v_sfc => zero
      END SELECT

      !$ACC ENTER DATA ASYNC(1) ATTACH(u_sfc, v_sfc)
        CALL vdiff_surface_flux ( &
            & jcs=ics, kproma=ice, dtime=delta_time, &
            & pfactor_sfc=prefactor_exchange(:), &
            & pcf=exchange_coeff_m_sfc(:,isfc), &
            & pen=uv_acoef(:,isfc), &
            & pfn=u_bcoef(:,isfc), &
            & x_sfc=u_sfc(:), &
            & flx=umfl_sft(:,isft) &
          )

        CALL vdiff_surface_flux ( &
            & jcs=ics, kproma=ice, dtime=delta_time, &
            & pfactor_sfc=prefactor_exchange(:), &
            & pcf=exchange_coeff_m_sfc(:,isfc), &
            & pen=uv_acoef(:,isfc), &
            & pfn=v_bcoef(:,isfc), &
            & x_sfc=v_sfc(:), &
            & flx=vmfl_sft(:,isft) &
          )
      !$ACC WAIT(1)
      !$ACC EXIT DATA DETACH(u_sfc, v_sfc)
    END DO

  END SUBROUTINE get_surface_stress


  !> Update tile-based variables in land and atmosphere state with results from JSBACH.
  !! Assumes that grid-scale variables are updated already. Grabs a few additional variables
  !! from inside JSBACH and updates the state with their contents (currently runoff and drainage)
  !! for a coupled HD model.
  SUBROUTINE update_nwp_tile_state (patch, delta_time, fr_sft, prog_lnd_new, diag_lnd, phy_diag)

    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.
    REAL(wp), INTENT(IN) :: delta_time !< Time step [s].
    REAL(wp), INTENT(IN) :: fr_sft(:,:,:) !< Surface-type fraction (nproma,nblk_c,SFT_NUM) [1].
    TYPE(t_lnd_prog), INTENT(INOUT) :: prog_lnd_new !< New prognostic land variables.
    TYPE(t_lnd_diag), INTENT(INOUT) :: diag_lnd !< Diagnostic land variables.
    TYPE(t_nwp_phy_diag), INTENT(INOUT) :: phy_diag !< Diagnostic atmosphere variables.

#ifndef __NO_JSBACH__
    REAL(wp), POINTER :: p_runoff(:,:)
    REAL(wp), POINTER :: p_drainage(:,:)

    REAL(wp) :: fr_land_lake
    INTEGER :: i_startblk, i_endblk, i_blk
    INTEGER :: ics, ice, ic

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    NULLIFY(p_runoff, p_drainage)

    CALL jsbach_get_var('hydro_runoff', patch%id, ptr2d=p_runoff, lacc=.TRUE.)
    CALL jsbach_get_var('hydro_drainage', patch%id, ptr2d=p_drainage, lacc=.TRUE.)
#endif

    !$OMP PARALLEL
      CALL copy(prog_lnd_new%t_g(:,:), prog_lnd_new%t_g_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(diag_lnd%t_s(:,:), prog_lnd_new%t_s_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(diag_lnd%t_s(:,:), prog_lnd_new%t_sk_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(diag_lnd%t_s(:,:), prog_lnd_new%t_so_t(:,1,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(diag_lnd%t_s(:,:), diag_lnd%t_sk(:,:), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(diag_lnd%t_s(:,:), diag_lnd%t_so(:,1,:), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(diag_lnd%qv_s(:,:), diag_lnd%qv_s_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)

      CALL copy(phy_diag%albnirdif(:,:), phy_diag%albnirdif_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(phy_diag%albvisdif(:,:), phy_diag%albvisdif_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(phy_diag%albdif(:,:), phy_diag%albdif_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)

      CALL copy(phy_diag%shfl_s(:,:), phy_diag%shfl_s_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(phy_diag%lhfl_s(:,:), phy_diag%lhfl_s_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(phy_diag%qhfl_s(:,:), phy_diag%qhfl_s_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(phy_diag%umfl_s(:,:), phy_diag%umfl_s_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy(phy_diag%vmfl_s(:,:), phy_diag%vmfl_s_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)

      CALL copy (phy_diag%u_10m(:,:), phy_diag%u_10m_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy (phy_diag%v_10m(:,:), phy_diag%v_10m_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)

      CALL copy (phy_diag%lhfl_pl(:,:,:), phy_diag%lhfl_pl_t(:,:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      CALL copy (phy_diag%lhfl_bs(:,:), phy_diag%lhfl_bs_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)

#ifndef __NO_JSBACH__
      !$ACC ENTER DATA ATTACH(p_runoff, p_drainage) ASYNC(1)
        CALL copy(p_runoff(:,:), diag_lnd%runoff_s_inst_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
        CALL copy(p_drainage(:,:), diag_lnd%runoff_g_inst_t(:,:,1), lacc=.TRUE., opt_acc_async=.TRUE.)
      !$ACC WAIT(1)
      !$ACC EXIT DATA DETACH(p_runoff, p_drainage)

      ! Convert runoff rate [kg/(m**2 s)] to the amount of runoff in the current time step
      ! [kg/m**2].

      !$OMP DO PRIVATE(i_blk, ics, ice, ic, fr_land_lake)
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c ( &
            & patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, end_prog_cells &
          )

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fr_land_lake)
          DO ic = ics, ice
            ! The runoff variables are given per m**2 of land+lake area. Convert to grid box values.
            fr_land_lake = fr_sft(ic,i_blk,SFT_LAND) + fr_sft(ic,i_blk,SFT_LWTR) &
                & + fr_sft(ic,i_blk,SFT_LICE)

            diag_lnd%runoff_s_inst_t(ic,i_blk,1) = delta_time * fr_land_lake &
                & * diag_lnd%runoff_s_inst_t(ic,i_blk,1)
            diag_lnd%runoff_g_inst_t(ic,i_blk,1) = delta_time * fr_land_lake &
                & * diag_lnd%runoff_g_inst_t(ic,i_blk,1)

            ! Accumulate runoff. For some reason this does not happen in mo_nwp_diagnosis.
            ! Unlike the tile-based runoff_{s,g}_t...
            diag_lnd%runoff_s(ic,i_blk) = diag_lnd%runoff_s(ic,i_blk) &
                & + diag_lnd%runoff_s_inst_t(ic,i_blk,1)
            diag_lnd%runoff_g(ic,i_blk) = diag_lnd%runoff_g(ic,i_blk) &
                & + diag_lnd%runoff_g_inst_t(ic,i_blk,1)
          END DO
        !$ACC END PARALLEL
      END DO
#endif
    !$OMP END PARALLEL

  END SUBROUTINE update_nwp_tile_state


  SUBROUTINE weighted_average_2d (w, x, xavg, pow)

    REAL(wp), INTENT(IN) :: w(:,:) !< Weights.
    REAL(wp), INTENT(IN) :: x(:,:) !< Values.
    REAL(wp), INTENT(OUT) :: xavg(:) !< Averaged value.
    REAL(wp), OPTIONAL, INTENT(IN) :: pow !< Power of the average (default: 1).

    INTEGER :: ic
    INTEGER :: k
    REAL(wp) :: ppow

    IF (PRESENT(pow)) THEN
      ! NVidia compiler doesn't like optionals.
      ppow = pow
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, SIZE(x, DIM=1)
          xavg(ic) = 0._wp
        END DO

        !$ACC LOOP SEQ
        DO k = 1, SIZE(x, DIM=2)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO ic = 1, SIZE(x, DIM=1)
            xavg(ic) = xavg(ic) + w(ic,k) * x(ic,k)**ppow
          END DO
        END DO

        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, SIZE(x, DIM=1)
          xavg(ic) = xavg(ic)**(1._wp/ppow)
        END DO
      !$ACC END PARALLEL

    ELSE

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, SIZE(x, DIM=1)
          xavg(ic) = 0._wp
        END DO

        !$ACC LOOP SEQ
        DO k = 1, SIZE(x, DIM=2)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO ic = 1, SIZE(x, DIM=1)
            xavg(ic) = xavg(ic) + w(ic,k) * x(ic,k)
          END DO
        END DO
      !$ACC END PARALLEL
    END IF

  END SUBROUTINE weighted_average_2d

  SUBROUTINE weighted_average_3d (patch, w, x, xavg, pow)

    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.
    REAL(wp), INTENT(IN) :: w(:,:,:) !< Weights.
    REAL(wp), INTENT(IN) :: x(:,:,:) !< Values.
    REAL(wp), INTENT(OUT) :: xavg(:,:) !< Averaged value.
    REAL(wp), OPTIONAL, INTENT(IN) :: pow !< Power of the average (default: 1).

    INTEGER :: i_startblk, i_endblk, i_blk
    INTEGER :: ics, ice, ic
    INTEGER :: k
    REAL(wp) :: ppow

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    IF (PRESENT(pow)) THEN
      ! NVidia compiler doesn't like optionals.
      ppow = pow

      !$OMP DO
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO ic = ics, ice
            xavg(ic,i_blk) = 0._wp
          END DO

          !$ACC LOOP SEQ
          DO k = 1, SIZE(x, DIM=3)
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO ic = ics, ice
              xavg(ic,i_blk) = xavg(ic,i_blk) + w(ic,i_blk,k) * x(ic,i_blk,k)**ppow
            END DO
          END DO

          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO ic = ics, ice
            xavg(ic,i_blk) = xavg(ic,i_blk)**(1._wp/ppow)
          END DO
        !$ACC END PARALLEL
      END DO

    ELSE

      !$OMP DO
      DO i_blk = i_startblk, i_endblk
        CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
            & end_prog_cells)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO ic = ics, ice
            xavg(ic,i_blk) = 0._wp
          END DO

          !$ACC LOOP SEQ
          DO k = 1, SIZE(x, DIM=3)
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO ic = ics, ice
              xavg(ic,i_blk) = xavg(ic,i_blk) + w(ic,i_blk,k) * x(ic,i_blk,k)
            END DO
          END DO
        !$ACC END PARALLEL
      END DO
    END IF

  END SUBROUTINE weighted_average_3d


  !>
  !! Update the state of the model with turbulence tendencies.
  SUBROUTINE add_tendencies_to_state ( &
        & dtime, patch, phy_tend, nh_prog, nh_prog_rcf, nh_diag, ddt_tracer &
      )

    REAL(wp), INTENT(IN) :: dtime !< Time step length [s].
    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.
    TYPE(t_nwp_phy_tend), INTENT(IN) :: phy_tend !< Physics tendencies.
    TYPE(t_nh_prog), INTENT(INOUT) :: nh_prog
    !< Prognostic variables of the current atmospheric state.
    TYPE(t_nh_prog), INTENT(INOUT) :: nh_prog_rcf
    !< Prognostic variables of the current atmospheric state (for tracers).
    TYPE(t_nh_diag), INTENT(INOUT) :: nh_diag
    !< Diagnostic variables of the current atmospheric state.
    REAL(wp), INTENT(IN) :: ddt_tracer(:,:,:,:) !< Extra tracer tendencies.

    LOGICAL :: nonnegative(ntracer)

    LOGICAL :: diffuse_tracer(ntracer)

    INTEGER :: i_startblk
    INTEGER :: i_endblk
    INTEGER :: i_blk

    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: ic
    INTEGER :: kl

    INTEGER :: it

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    diffuse_tracer(1:nqtendphy) = .TRUE.
    diffuse_tracer(nqtendphy+1:iqt-1) = .FALSE.
    diffuse_tracer(iqt:) = .TRUE.

    !$ACC ENTER DATA CREATE(nonnegative) ASYNC(1)

    !$ACC SERIAL DEFAULT(PRESENT) ASYNC(1)
      nonnegative(:) = .FALSE.
      IF (iqv > 0) nonnegative(iqv) = .TRUE.
      IF (iqc > 0) nonnegative(iqc) = .TRUE.
      IF (iqi > 0) nonnegative(iqi) = .TRUE.
    !$ACC END SERIAL

    !$OMP PARALLEL
    !$OMP DO PRIVATE(i_blk, ics, ice, ic, kl, it)
    DO i_blk = i_startblk, i_endblk
      CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
          & end_prog_cells)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO kl = 1, patch%nlev
          DO ic = ics, ice
            nh_diag%u(ic,kl,i_blk) = nh_diag%u(ic,kl,i_blk) &
                & + dtime * phy_tend%ddt_u_turb(ic,kl,i_blk)
            nh_diag%v(ic,kl,i_blk) = nh_diag%v(ic,kl,i_blk) &
                & + dtime * phy_tend%ddt_v_turb(ic,kl,i_blk)
            nh_diag%temp(ic,kl,i_blk) = nh_diag%temp(ic,kl,i_blk) &
                & + dtime * phy_tend%ddt_temp_turb(ic,kl,i_blk)
          END DO
        END DO
      !$ACC END PARALLEL

      IF (ASSOCIATED(phy_tend%ddt_w_turb)) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO kl = 1, patch%nlev
            DO ic = ics, ice
              nh_prog%w(ic,kl,i_blk) = nh_prog%w(ic,kl,i_blk) &
                & + phy_tend%ddt_w_turb(ic,kl,i_blk)
            END DO
          END DO
        !$ACC END PARALLEL
      END IF

      ! TODO: Check handling of other tracers!
      DO it = 1, ntracer
        IF (.NOT. diffuse_tracer(it)) CYCLE

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO kl = 1, patch%nlev
            DO ic = ics, ice
              nh_prog_rcf%tracer(ic,kl,i_blk,it) = nh_prog_rcf%tracer(ic,kl,i_blk,it) &
                  & + dtime * ddt_tracer(ic,kl,i_blk,it)
            END DO
          END DO

          IF (nonnegative(it)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
            DO kl = 1, patch%nlev
              DO ic = ics, ice
                nh_prog_rcf%tracer(ic,kl,i_blk,it) = &
                  & MAX(0._wp, nh_prog_rcf%tracer(ic,kl,i_blk,it))
              END DO
            END DO
          END IF
        !$ACC END PARALLEL
      END DO
    END DO
    !$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(nonnegative)

  END SUBROUTINE add_tendencies_to_state


  !>
  !! Call the orbit routine to update the current declination.
  SUBROUTINE update_earth_declination (dt)

    TYPE(datetime), TARGET, INTENT(IN) :: dt

    TYPE(julianday), TARGET :: jd

    REAL(wp) :: jd_real
    REAL(wp) :: rasc_sun
    REAL(wp) :: decl_sun
    REAL(wp) :: dist_sun

    REAL(wp), PARAMETER :: ms_per_day = 1000._wp * 60._wp * 60._wp * 24._wp

    CALL getJulianDayFromDatetime(dt, jd)

    jd_real = REAL(jd%day, wp) + jd%ms / ms_per_day

    CALL orbit_vsop87(jd_real, rasc_sun, decl_sun, dist_sun)

  END SUBROUTINE update_earth_declination

  !> Add all JSBACH restart variables to the jsb_init_vars group.
  !! This group is later used to create input instructions for initialization.
  !! We also use it to output the JSBACH state. Thus, we have to make all
  !! variables output variables.
  SUBROUTINE setup_jsbach_init_vars

    USE mo_var_groups, ONLY: var_groups_dyn
    USE mo_var_list_register, ONLY: t_vl_register_iter

    TYPE(t_vl_register_iter) :: iter
    INTEGER :: i
    INTEGER :: grp

    grp = var_groups_dyn%group_id('jsb_init_vars')

    DO WHILE(iter%next())
      IF (iter%cur%p%vlname(1:4) /= 'jsb_') CYCLE

      ASSOCIATE (vl => iter%cur%p)
        DO i = 1, vl%nvars
          IF (.NOT. vl%vl(i)%p%info%lrestart) CYCLE
            vl%vl(i)%p%info%in_group(grp) = .TRUE.
            vl%vl(i)%p%info%loutput = .TRUE.
        END DO
      END ASSOCIATE
    END DO

  END SUBROUTINE setup_jsbach_init_vars

END MODULE mo_nwp_vdiff_interface
