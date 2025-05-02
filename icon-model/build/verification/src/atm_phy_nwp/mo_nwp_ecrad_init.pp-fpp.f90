# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_nwp/mo_nwp_ecrad_init.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_nwp/mo_nwp_ecrad_init.f90"
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

! In this module the configuration state for the ecRad radiation code is being set up.
!
! - Setup information is stored inside the object ecrad_conf of the derived type t_ecrad_conf
!   containing all the configuration information needed to run the radiation scheme.
! - The intention is that this part of the configuration is fixed for a given model run.
! - ICON namelist settings are translated to ecRad conform settings, if unsupported values
!   are provided via namelist, the user gets an error. (These values should already be
!   checked by the nml_crosscheck)
! - Please note that only a subset of the available configuration options of ecRad is
!   filled by this routine. For a full list of ecRad settings, please have a look at
!   externals/ecrad/radiation/radiation_config.F90

!----------------------------

# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/include/omp_definitions.inc" 1
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

# 26 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_nwp/mo_nwp_ecrad_init.f90" 2
!----------------------------

MODULE mo_nwp_ecrad_init

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message
  USE mo_radiation_config,     ONLY: icld_overlap, irad_aero, ecrad_data_path,             &
                                 &   ecrad_isolver, ecrad_igas_model, isolrad,             &
                                 &   ecrad_llw_cloud_scat, ecrad_iliquid_scat,             &
                                 &   ecrad_iice_scat, ecrad_isnow_scat,                    &
                                 &   ecrad_irain_scat, ecrad_igraupel_scat,                &
                                 &   ecrad_use_general_cloud_optics,                       &
                                 &   iRadAeroConst, iRadAeroTegen, iRadAeroART,            &
                                 &   iRadAeroConstKinne, iRadAeroKinne, iRadAeroVolc,      &
                                 &   iRadAeroKinneVolc,  iRadAeroKinneVolcSP,              &
                                 &   iRadAeroKinneSP, iRadAeroNone,                        &
                                 &   iRadAeroCAMSclim, iRadAeroCAMStd

  USE mo_ecrad,                ONLY: t_ecrad_conf, ecrad_setup,                            &
                                 &   ISolverHomogeneous, ISolverMcICA, ISolverMcICAACC,    &
                                 &   ISolverSpartacus,                                     &
                                 &   ISolverTripleclouds, ISolverCloudless,                &
                                 &   IGasModelMonochromatic, IGasModelIFSRRTMG,            &
                                 &   IGasModelECCKD,                                       &
                                 &   ILiquidModelMonochromatic, ILiquidModelSOCRATES,      &
                                 &   ILiquidModelSlingo, IIceModelMonochromatic,           &
                                 &   IIceModelFu, IIceModelBaran,                          &
                                 &   IIceModelBaran2016, IIceModelBaran2017, IIceModelYi,  &
                                 &   IOverlapMaximumRandom, IOverlapExponentialRandom,     &
                                 &   IOverlapExponential,                                  &
                                 &   nweight_nir_ecrad, iband_nir_ecrad, weight_nir_ecrad, &
                                 &   nweight_vis_ecrad, iband_vis_ecrad, weight_vis_ecrad, &
                                 &   nweight_par_ecrad, iband_par_ecrad, weight_par_ecrad, &
                                 &   ecrad_hyd_list,                                       &
                                 &   ecrad_iqc, ecrad_iqi, ecrad_iqr, ecrad_iqs, ecrad_iqg



  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_init'


  PUBLIC :: setup_ecrad


CONTAINS


  !---------------------------------------------------------------------------------------
  SUBROUTINE setup_ecrad ( ecrad_conf )

    CHARACTER(len=*), PARAMETER :: routine = modname//'::setup_ecrad'

    TYPE(t_ecrad_conf),  INTENT(inout) :: &
      &  ecrad_conf                         !< ecRad configuration state

    ! Local variables
    REAL(wp)                  :: &
      &  wavelength_bound_sw(1), & !< Wavelength bound between VIS and NIR albedo (m)
      &  wavelength_bound_lw(0)    !< Wavelength bound LW emissivity (m)
    INTEGER                   :: &
      &  i_band_in_sw(2),        & !< The albedo band indices corresponding to each interval
      &  i_band_in_lw(1),        & !< The emissivity band indices corresponding to each interval
      &  cc_cloud                  !< counter for cloud_types

    CALL message(routine, 'Setup of ecRad radiation')

    !---------------------------------------------------------------------------------------
    ! Checks
    !---------------------------------------------------------------------------------------

    ! Compatibility check wp and JPRB. If this check fails, JPRB has to be adapted manually to wp.
    IF (PRECISION(wavelength_bound_sw) /= PRECISION(ecrad_conf%cloud_fraction_threshold)) &
      &  CALL finish(routine,'ICON working precision (wp) does not match ecRad precision.')
    IF (EPSILON(wavelength_bound_sw(1)) /= EPSILON(ecrad_conf%cloud_fraction_threshold)) &
      &  CALL finish(routine,'Smallest number in working precision (wp) is different from ecRad precision.')

    !---------------------------------------------------------------------------------------
    ! Configuration based on ICON namelist settings
    !---------------------------------------------------------------------------------------

    ! Directory with all input data required by ecRad
    ecrad_conf%directory_name = ecrad_data_path

    ! Overlap scheme
    SELECT CASE (icld_overlap)
      CASE (1)
        ecrad_conf%i_overlap_scheme = IOverlapMaximumRandom
      CASE (2)
        ecrad_conf%i_overlap_scheme = IOverlapExponentialRandom
      CASE (5)
        ecrad_conf%i_overlap_scheme = IOverlapExponential
      CASE DEFAULT
        CALL finish(routine, 'Only icld_overlap values of 1 (MAX-RAN), 2 (EXP-RAN) and 5 (EXP) are valid for ecRad')
    END SELECT

    ! Aerosol climatology
    SELECT CASE (irad_aero)
      CASE (iRadAeroNone) ! No aerosol
        ecrad_conf%use_aerosols = .false.
      CASE (iRadAeroConst, iRadAeroTegen, iRadAeroCAMSclim, iRadAeroCAMStd, iRadAeroART, iRadAeroConstKinne,  &
        &   iRadAeroKinne, iRadAeroVolc, iRadAeroKinneVolc,  iRadAeroKinneVolcSP, iRadAeroKinneSP)
        ecrad_conf%use_aerosols = .true.
      CASE DEFAULT
        CALL finish(routine, 'irad_aero not valid for ecRad')
    END SELECT

    ! LW scattering due to clouds
    ecrad_conf%do_lw_cloud_scattering = ecrad_llw_cloud_scat

    ! Select if generalized optics are used. This might be overruled by ecrad_igas_model later.
    ecrad_conf%use_general_cloud_optics = ecrad_use_general_cloud_optics

    ! Radiation solver
    SELECT CASE (ecrad_isolver)
      CASE(0)
        ecrad_conf%i_solver_sw   = ISolverMcICA        !< Short-wave solver
        ecrad_conf%i_solver_lw   = ISolverMcICA        !< Long-wave solver
        ecrad_conf%do_3d_effects = .false.             !< Do we include 3D effects?
      CASE(1)
        ecrad_conf%i_solver_sw   = ISolverTripleclouds !< Short-wave solver
        ecrad_conf%i_solver_lw   = ISolverTripleclouds !< Long-wave solver
        ecrad_conf%do_3d_effects = .false.             !< Do we include 3D effects?
      CASE(2)
        ecrad_conf%i_solver_sw   = ISolverMcICAACC     !< Short-wave solver
        ecrad_conf%i_solver_lw   = ISolverMcICAACC     !< Long-wave solver
        ecrad_conf%do_3d_effects = .false.             !< Do we include 3D effects?
      CASE(3)
        ecrad_conf%i_solver_sw   = ISolverSpartacus    !< Short-wave solver
        ecrad_conf%i_solver_lw   = ISolverSpartacus    !< Long-wave solver
        ecrad_conf%do_3d_effects = .true.              !< Do we include 3D effects?
      CASE DEFAULT
        CALL finish(routine, 'ecrad_isolver not valid for ecRad')
    END SELECT

    ! Gas model and spectral bands: RRTMG or ecckd
    SELECT CASE (ecrad_igas_model)
      CASE(0)
        ecrad_conf%i_gas_model_sw                  = IGasModelIFSRRTMG  !< Use RRTM gas model
        ecrad_conf%i_gas_model_lw                  = IGasModelIFSRRTMG  !< Use RRTM gas model
        ecrad_conf%do_cloud_aerosol_per_lw_g_point = .false.
        ecrad_conf%do_cloud_aerosol_per_sw_g_point = .false.
      CASE(1)
        ecrad_conf%i_gas_model_sw                  = IGasModelECCKD     !< Use ecckd gas model
        ecrad_conf%i_gas_model_lw                  = IGasModelECCKD     !< Use ecckd gas model
        ecrad_conf%do_cloud_aerosol_per_lw_g_point = .true.
        ecrad_conf%do_cloud_aerosol_per_sw_g_point = .true.
      CASE DEFAULT
        CALL finish(routine, 'ecrad_igas_model not valid for ecRad')
    END SELECT

    ! Generalized hydrometeors
    IF ( ecrad_conf%use_general_cloud_optics ) THEN
      ! Cont how many cloud types are.
      ecrad_conf%n_cloud_types   = 2
      IF ( ecrad_isnow_scat > -1 )    ecrad_conf%n_cloud_types = ecrad_conf%n_cloud_types + 1
      IF ( ecrad_irain_scat > -1 )    ecrad_conf%n_cloud_types = ecrad_conf%n_cloud_types + 1
      IF ( ecrad_igraupel_scat > -1 ) ecrad_conf%n_cloud_types = ecrad_conf%n_cloud_types + 1

      ALLOCATE (ecrad_hyd_list(ecrad_conf%n_cloud_types))

      ecrad_hyd_list(:)        = 0
      cc_cloud                 = 1
      ecrad_hyd_list(cc_cloud) = ecrad_iqc

      SELECT CASE (ecrad_iliquid_scat)
        CASE(0)
          ecrad_conf%cloud_type_name(cc_cloud) = "mie_droplet"
        CASE DEFAULT
          CALL finish(routine, 'ecrad_iliquid_scat not valid for ecRad and use_general_cloud_optics = T')
      END SELECT

      cc_cloud = cc_cloud + 1
      ecrad_hyd_list(cc_cloud) = ecrad_iqi

      SELECT CASE (ecrad_iice_scat)
        CASE(0)
          ecrad_conf%cloud_type_name(cc_cloud) = "fu-muskatel-rough_ice"
        CASE(10)
          ecrad_conf%cloud_type_name(cc_cloud) = "fu-muskatel_ice"
        CASE(11)
          ecrad_conf%cloud_type_name(cc_cloud) = "baum-general-habit-mixture_ice"
        CASE DEFAULT
          CALL finish(routine, 'ecrad_iice_scat not valid for ecRad and use_general_cloud_optics = T')
      END SELECT

      IF (ecrad_isnow_scat > -1) THEN ! Snow is optional for radiation
        cc_cloud = cc_cloud + 1
        ecrad_hyd_list(cc_cloud) = ecrad_iqs

        SELECT CASE (ecrad_isnow_scat)
          CASE(0)
            ecrad_conf%cloud_type_name(cc_cloud) = "fu-muskatel-rough_ice"
          CASE(10)
            ecrad_conf%cloud_type_name(cc_cloud) = "fu-muskatel_ice"
          CASE DEFAULT
            CALL finish(routine, 'ecrad_isnow_scat not valid for ecRad and use_general_cloud_optics = T')
        END SELECT
      ENDIF

      IF (ecrad_irain_scat > -1) THEN ! Rain is optional for radiation
        cc_cloud = cc_cloud + 1
        ecrad_hyd_list(cc_cloud) = ecrad_iqr

        SELECT CASE (ecrad_irain_scat)
          CASE(0)
            ecrad_conf%cloud_type_name(cc_cloud) = "mie_rain"
          CASE DEFAULT
            CALL finish(routine, 'ecrad_irain_scat not valid for ecRad and use_general_cloud_optics = T')
          END SELECT
      ENDIF

      IF (ecrad_igraupel_scat > -1) THEN ! Graupel is optional for radiation
        cc_cloud = cc_cloud + 1
        ecrad_hyd_list(cc_cloud) = ecrad_iqg

        SELECT CASE (ecrad_igraupel_scat)
          CASE(0)
            ecrad_conf%cloud_type_name(cc_cloud) = "fu-muskatel-rough_ice"
          CASE(10)
            ecrad_conf%cloud_type_name(cc_cloud) = "fu-muskatel_ice"
          CASE DEFAULT
            CALL finish(routine, 'ecrad_igraupel_scat not valid for ecRad and use_general_cloud_optics = T')
        END SELECT
      ENDIF
      
    ELSE ! .not.ecrad_conf%use_general_cloud_optics
      
      ! Liquid cloud particle scattering properties
      SELECT CASE (ecrad_iliquid_scat)
        CASE(0)
          ecrad_conf%i_liq_model = ILiquidModelSOCRATES
        CASE(1)
          ecrad_conf%i_liq_model = ILiquidModelSlingo
        CASE DEFAULT
          CALL finish(routine, 'ecrad_iliquid_scat not valid for ecRad and use_general_cloud_optics = F')
      END SELECT
    
      ! Ice cloud particle scattering properties
      SELECT CASE (ecrad_iice_scat)
        CASE(0)
          ecrad_conf%i_ice_model = IIceModelFu
        CASE(1)
          ecrad_conf%i_ice_model = IIceModelBaran2016
        CASE(2)
          ecrad_conf%i_ice_model = IIceModelYi
        CASE DEFAULT
          CALL finish(routine, 'ecrad_iice_scat not valid for ecRad and use_general_cloud_optics = F')
      END SELECT

    END IF ! ecrad_conf%use_general_cloud_optics

    IF (ecrad_conf%i_gas_model_sw == IGasModelIFSRRTMG .AND. ecrad_conf%i_gas_model_lw == IGasModelIFSRRTMG) THEN
      ecrad_conf%do_setup_ifsrrtm = .true.
    ELSE IF (ecrad_conf%i_gas_model_sw .NE. ecrad_conf%i_gas_model_lw) THEN
      CALL finish(routine, "Differing gas models for LW and SW are currently unsupported by ICON.")
    ELSE
      ecrad_conf%do_setup_ifsrrtm = .false.
    ENDIF

    IF (isolrad == 0 .OR. ecrad_igas_model == 1) THEN
      ! ecckd applies Coddington scaling by default
      ecrad_conf%use_spectral_solar_scaling  = .false.
    ELSE
      ecrad_conf%use_spectral_solar_scaling  = .true.
    ENDIF

    !---------------------------------------------------------------------------------------
    ! Currently hardcoded configuration
    !---------------------------------------------------------------------------------------
  
    ecrad_conf%do_lw                       = .true.       !< Do we compute longwave radiation?
    !
    ecrad_conf%do_sw                       = .true.       !< Do we compute shortwave radiation?
    !
    ecrad_conf%do_clear                    = .true.       !< Do we compute clear-sky fluxes?
    !
    ecrad_conf%do_sw_direct                = .true.       !< Do we compute solar direct fluxes?
    !
    ecrad_conf%do_lw_aerosol_scattering    = .false.      !< LW scattering due to aerosol
    !
    ecrad_conf%use_beta_overlap            = .false.      !< Use Shonk et al. (2010) "beta" overlap parameter
                                                          !< instead of "alpha" (Hogan and Illingworth, 2000)
    !
    ecrad_conf%iverbosesetup               = 0            !< Verbosity (0: none,     1: warning,  2: info,
    ecrad_conf%iverbose                    = 0            !<            3: progress, 4: detailed, 5: debug)
    !
    ecrad_conf%do_weighted_surface_mapping = .false.      !<Use spectral weighting method that was default before ecrad-1.5
    !
    ecrad_conf%do_surface_sw_spectral_flux = .true.       !< Save the surface downwelling shortwave fluxes in each band
                                                          !< Needed for photosynthetic active radiation
    !
    ecrad_conf%do_fu_lw_ice_optics_bug     = .false.      !< In the IFS environment there was a bug in the Fu longwave
                                                          !< ice optics producing better results than the fixed version
    !
    ecrad_conf%do_sw_delta_scaling_with_gases = .false.   !< Do SW delta-Eddington scaling on cloud-aerosol-gas mixture (.true.)
                                                          !< More correct approach of separately scaling the cloud and aerosol 
                                                          !< scattering properties before merging with gases (.false.)
    !
    ecrad_conf%cloud_fraction_threshold    = 1.0e-6_wp    !< Cloud is present in a layer if cloud fraction exceeds this value
    ecrad_conf%cloud_mixing_ratio_threshold= 1.0e-9_wp    !< ..and total cloud water mixing ratio exceeds this value
    !
    ecrad_conf%cloud_inhom_decorr_scaling  = 0.5_wp       !< Ratio of the overlap decorrelation length for cloud inhomogeneities
                                                          !< to the overlap decorrelation length for cloud boundaries.
                                                          !< Observations suggest this has a value of 0.5
    !
    ecrad_conf%min_gas_od_lw               = 1.0e-15_wp   !< Minimum gas optical depth, for stability (long-wave)
    !
    ecrad_conf%min_gas_od_sw               = 0.0_wp       !< Minimum gas optical depth, for stability (short-wave)
    !
    ecrad_conf%max_cloud_od                = 20.0_wp      !< Maximum total optical depth of a cloudy region, for stability

    ! Optical properties data is taken from aerosol_ifs_rrtm_46R1_with_NI_AM.nc :  
    IF (irad_aero == iRadAeroCAMSclim .OR. irad_aero == iRadAeroCAMStd) THEN
      ecrad_conf%use_aerosols              = .true.
      ecrad_conf%n_aerosol_types           = 11
      ecrad_conf%i_aerosol_type_map(1)     = -1           !< aermr01  Sea Salt Aerosol (0.03 - 0.5 um)   hydrophilic(1)
      ecrad_conf%i_aerosol_type_map(2)     = -2           !< aermr02  Sea Salt Aerosol (0.5 - 5 um)      hydrophilic(2)
      ecrad_conf%i_aerosol_type_map(3)     = -3           !< aermr03  Sea Salt Aerosol (5 - 20 um)       hydrophilic(3)
      ecrad_conf%i_aerosol_type_map(4)     =  1           !< aermr04  Dust Aerosol (0.03 - 0.55 um)      hydrophobic(1)
      ecrad_conf%i_aerosol_type_map(5)     =  2           !< aermr05  Dust Aerosol (0.55 - 0.9 um)       hydrophobic(2)
      ecrad_conf%i_aerosol_type_map(6)     =  3           !< aermr06  Dust Aerosol (0.9 - 20 um)         hydrophobic(3)
      ecrad_conf%i_aerosol_type_map(7)     = -4           !< aermr07  Hydrophilic Organic Matter Aerosol hydrophilic(4)
      ecrad_conf%i_aerosol_type_map(8)     = 10           !< aermr08  Hydrophobic Organic Matter Aerosol hydrophobic(10)
      ecrad_conf%i_aerosol_type_map(9)     = 11           !< aermr09  Hydrophilic Black Carbon Aerosol   hydrophobic(11)
      ecrad_conf%i_aerosol_type_map(10)    = 11           !< aermr10  Hydrophobic Black Carbon Aerosol   hydrophobic(11)
      ecrad_conf%i_aerosol_type_map(11)    = -5           !< aermr11  Sulphate Aerosol                   hydrophilic(5)
    ENDIF

    !---------------------------------------------------------------------------------------
    ! Call to ecRad setup routine. This also consolidates the configuration
    !---------------------------------------------------------------------------------------

    CALL ecrad_setup(ecrad_conf)

    !---------------------------------------------------------------------------------------
    ! Tell ecRad about the wavelength bounds of ICON data
    !---------------------------------------------------------------------------------------

    ! Set up the near-IR, visible, and photosynthetically active radiation wavelength bounds.
    ! Bounds for nir and vis are from mo_nwp_phy_types.
    CALL ecrad_conf%get_sw_weights(0.7e-6_wp, 5.0e-6_wp,       &
      &  nweight_nir_ecrad, iband_nir_ecrad, weight_nir_ecrad, &
      &  'near-IR radiation')
    CALL ecrad_conf%get_sw_weights(0.3e-6_wp, 0.7e-6_wp,       &
      &  nweight_vis_ecrad, iband_vis_ecrad, weight_vis_ecrad, &
      &  'visible radiation')
    CALL ecrad_conf%get_sw_weights(0.4e-6_wp, 0.7e-6_wp,       &
      &  nweight_par_ecrad, iband_par_ecrad, weight_par_ecrad, &
      &  'photosynthetically active radiation, PAR')

    !$ACC UPDATE DEVICE(iband_nir_ecrad, weight_nir_ecrad) ASYNC(1)
    !$ACC UPDATE DEVICE(iband_vis_ecrad, weight_vis_ecrad) ASYNC(1)
    !$ACC UPDATE DEVICE(iband_par_ecrad, weight_par_ecrad) ASYNC(1)

    ! ICON external parameters have SW albedo for two different wavelength bands, visible and near infrared. The following call to
    ! ecrad_conf%define_sw_albedo_intervals tells ecrad about the two bands and the wavelength bound which is at 700 nm (according
    ! to a comment in mo_nwp_phy_types).
    wavelength_bound_sw(1) = 0.7_wp * 1.e-6_wp !< 700 nm
    i_band_in_sw(1)        = 1
    i_band_in_sw(2)        = 2
    CALL ecrad_conf%define_sw_albedo_intervals(2, wavelength_bound_sw, i_band_in_sw, do_nearest=.false.)

    ! Similar to the short wave albedo bands, ecRad needs to know the number of longwave emissivity bands provided from ICON
    ! external data.
    i_band_in_lw           = 1
    CALL ecrad_conf%define_lw_emiss_intervals(1, wavelength_bound_lw, i_band_in_lw, do_nearest=.true.)

    !$ACC ENTER DATA COPYIN(ecrad_conf)
    !$ACC ENTER DATA COPYIN(ecrad_conf%cloud_optics) &
    !$ACC   COPYIN(ecrad_conf%i_albedo_from_band_sw) &
    !$ACC   COPYIN(ecrad_conf%i_band_from_reordered_g_lw) &
    !$ACC   COPYIN(ecrad_conf%i_band_from_reordered_g_sw) &
    !$ACC   COPYIN(ecrad_conf%i_band_from_g_lw) &
    !$ACC   COPYIN(ecrad_conf%i_g_from_reordered_g_lw) &
    !$ACC   COPYIN(ecrad_conf%i_emiss_from_band_lw) &
    !$ACC   COPYIN(ecrad_conf%pdf_sampler) &
    !$ACC   COPYIN(ecrad_conf%sw_albedo_weights)
    !$ACC ENTER DATA COPYIN(ecrad_conf%cloud_optics%liq_coeff_lw) &
    !$ACC   COPYIN(ecrad_conf%cloud_optics%liq_coeff_sw) &
    !$ACC   COPYIN(ecrad_conf%cloud_optics%ice_coeff_lw) &
    !$ACC   COPYIN(ecrad_conf%cloud_optics%ice_coeff_sw) &
    !$ACC   COPYIN(ecrad_conf%pdf_sampler%val)





  END SUBROUTINE setup_ecrad
  !---------------------------------------------------------------------------------------

# 569 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_nwp/mo_nwp_ecrad_init.f90"


END MODULE mo_nwp_ecrad_init
