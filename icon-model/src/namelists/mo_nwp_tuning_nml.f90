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

! Namelist for tuning and/or perturbing nwp physics
!
! These subroutines are called by read_atmo_namelists and do some
! nwp physics tuning

MODULE mo_nwp_tuning_nml

  USE mo_exception,           ONLY: finish
  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_nwp_tuning_config,   ONLY: config_tune_gkwake           => tune_gkwake,           &
    &                               config_tune_gkdrag           => tune_gkdrag,           &
    &                               config_tune_gkdrag_enh       => tune_gkdrag_enh,       &
    &                               config_tune_gfrcrit          => tune_gfrcrit,          &
    &                               config_tune_grcrit           => tune_grcrit,           &
    &                               config_tune_grcrit_enh       => tune_grcrit_enh,       &
    &                               config_tune_minsso           => tune_minsso,           &
    &                               config_tune_minsso_gwd       => tune_minsso_gwd,       &
    &                               config_tune_blockred         => tune_blockred,         &
    &                               config_tune_gfluxlaun        => tune_gfluxlaun,        &
    &                               config_tune_gcstar           => tune_gcstar,           &
    &                               config_tune_zceff_min        => tune_zceff_min,        &
    &                               config_tune_v0snow           => tune_v0snow,           &
    &                               config_tune_zvz0i            => tune_zvz0i,            &
    &                               config_tune_icesedi_exp      => tune_icesedi_exp,      &
    &                               config_tune_entrorg          => tune_entrorg,          &
    &                               config_tune_rprcon           => tune_rprcon,           &
    &                               config_tune_rdepths          => tune_rdepths,          &
    &                               config_tune_capdcfac_et      => tune_capdcfac_et,      &
    &                               config_tune_capdcfac_tr      => tune_capdcfac_tr,      &
    &                               config_tune_capethresh       => tune_capethresh,       &
    &                               config_tune_rhebc_land       => tune_rhebc_land,       &
    &                               config_tune_rhebc_ocean      => tune_rhebc_ocean,      &
    &                               config_tune_rcucov           => tune_rcucov,           &
    &                               config_tune_rhebc_land_trop  => tune_rhebc_land_trop,  &
    &                               config_tune_rhebc_ocean_trop => tune_rhebc_ocean_trop, &
    &                               config_tune_rcucov_trop      => tune_rcucov_trop,      &
    &                               config_tune_texc             => tune_texc,             &
    &                               config_tune_qexc             => tune_qexc,             &
    &                               config_tune_rcapqadv         => tune_rcapqadv,         &
    &                               config_tune_minsnowfrac      => tune_minsnowfrac,      &
    &                               config_tune_box_liq          => tune_box_liq,          &
    &                               config_tune_box_liq_asy      => tune_box_liq_asy,      &
    &                               config_tune_box_liq_sfc_fac  => tune_box_liq_sfc_fac,  &
    &                               config_allow_overcast        => allow_overcast,        &
    &                               config_tune_thicklayfac      => tune_thicklayfac,      &
    &                               config_tune_sgsclifac        => tune_sgsclifac,        &
    &                               config_icpl_turb_clc         => icpl_turb_clc,         &
    &                               config_tune_dust_abs         => tune_dust_abs,         &
    &                               config_tune_difrad_3dcont    => tune_difrad_3dcont,    &
    &                               config_tune_gust_factor      => tune_gust_factor,      &
    &                               config_itune_gust_diag       => itune_gust_diag,       &
    &                               config_tune_gustsso_lim      => tune_gustsso_lim,      &
    &                               config_tune_gustlim_agl      => tune_gustlim_agl,      &
    &                               config_tune_gustlim_fac      => tune_gustlim_fac,      &
    &                               config_itune_albedo          => itune_albedo,          &
    &                               config_tune_albedo_wso       => tune_albedo_wso,       &
    &                               config_itune_slopecorr       => itune_slopecorr,       &
    &                               config_itune_o3              => itune_o3,              &
    &                               config_lcalib_clcov          => lcalib_clcov,          &
    &                               config_max_calibfac_clcl     => max_calibfac_clcl,     &
    &                               config_max_freshsnow_inc     => max_freshsnow_inc,     &
    &                               config_tune_eiscrit          => tune_eiscrit,          &
    &                               config_tune_sc_eis           => tune_sc_eis,           &
    &                               config_tune_sc_invmin        => tune_sc_invmin,        &
    &                               config_tune_sc_invmax        => tune_sc_invmax,        &
    &                               config_tune_dursun_scaling   => tune_dursun_scaling,   &
    &                               config_tune_sbmccn           => tune_sbmccn,           &
    &                               config_tune_urbahf           => tune_urbahf,           &
    &                               config_tune_urbisa           => tune_urbisa
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_nwp_tuning_namelist


  !-----------------------------------!
  ! nwp_tuning_nml namelist variables !
  !-----------------------------------!

  REAL(wp) :: &                    !< low level wake drag constant
    &  tune_gkwake(max_dom)

  REAL(wp) :: &                    !< gravity wave drag constant; optional enhanced value for low latitudes
    &  tune_gkdrag(max_dom), tune_gkdrag_enh(max_dom)

  REAL(wp) :: &                    !< critical Froude number in SSO scheme
    &  tune_gfrcrit(max_dom)

  REAL(wp) :: &                    !< critical Richardson number in SSO scheme; optional enhanced value for low latitudes
    &  tune_grcrit(max_dom), tune_grcrit_enh(max_dom)

  REAL(wp) :: &                    !< minimum SSO standard deviation (m) for which SSO information is used
    &  tune_minsso(max_dom)

  REAL(wp) :: &                    !< minimum SSO standard deviation (m) above which GWD computation is activated
    &  tune_minsso_gwd(max_dom)    !  provided that the value is larger than tune_minsso

  REAL(wp) :: &                    !< multiple of SSO standard deviation above which blocking tendency is reduced
    &  tune_blockred(max_dom)

  REAL(wp) :: &                    !< total launch momentum flux in each azimuth (rho_o x F_o)
    &  tune_gfluxlaun

  REAL(wp) :: &                    !< constant in saturation wave spectrum (non-orographic GWD)
    &  tune_gcstar
  
  REAL(wp) :: &                    !< Minimum value for sticking efficiency
    &  tune_zceff_min

  REAL(wp) :: &                    !< factor in the terminal velocity for snow
    &  tune_v0snow

  REAL(wp) :: &                    !< Terminal fall velocity of ice 
    &  tune_zvz0i

  REAL(wp) :: &                    !< Exponent for density correction of cloud ice sedimentation
    &  tune_icesedi_exp

  REAL(wp) :: &                    !< sbm ccn factor
    &  tune_sbmccn

  REAL(wp) :: &                    !< Entrainment parameter for deep convection valid at dx=20 km 
    &  tune_entrorg

  REAL(wp) :: &                    !< Coefficient for conversion of cloud water into precipitation in convection scheme 
    &  tune_rprcon

  REAL(wp) :: &                    !< Maximum allowed shallow convection depth (Pa) 
    &  tune_rdepths

  REAL(wp) :: &                    !< Fraction of CAPE diurnal cycle correction applied in the extratropics
    &  tune_capdcfac_et            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< Fraction of CAPE diurnal cycle correction applied in the tropics
    &  tune_capdcfac_tr            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< CAPE threshold above which the convective adjustment time scale and entrainment
    &  tune_capethresh             !< are reduced for numerical stability [J/kg]

  REAL(wp) :: &                    !< RH threshold for onset of evaporation below cloud base over land
    &  tune_rhebc_land

  REAL(wp) :: &                    !< RH threshold for onset of evaporation below cloud base over sea
    &  tune_rhebc_ocean

  REAL(wp) :: &                    !< Convective area fraction
    &  tune_rcucov

  REAL(wp) :: &                    !< RH threshold for onset of evaporation below cloud base over tropical land
    &  tune_rhebc_land_trop        !  (relevant only if smaller than rhebc_land)

  REAL(wp) :: &                    !< RH threshold for onset of evaporation below cloud base over tropical sea
    &  tune_rhebc_ocean_trop       !  (relevant only if smaller than rhebc_ocean)

  REAL(wp) :: &                    !< Convective area fraction in the tropics
    &  tune_rcucov_trop            !  (relevant only if smaller than rcucov)

  REAL(wp) :: &                    !< Excess value for temperature used in test parcel ascent
    &  tune_texc

  REAL(wp) :: &                    !< Excess fraction of grid-scale QV used in test parcel ascent
    &  tune_qexc

  REAL(wp) :: &                    !< Factor for dynamic correction of cape closure
    &  tune_rcapqadv

  REAL(wp) :: &                    !< Minimum value to which the snow cover fraction is artificially reduced
    &  tune_minsnowfrac            !  in case of melting show (in case of idiag_snowfrac = 20)

  REAL(wp) :: &                    !< Box width for liquid clouds assumed in the cloud cover scheme
    &  tune_box_liq                ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Factor for increasing the box width in case of thick model layers
    &  tune_thicklayfac            ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Asymmetry factor liquid cloud parameterization
    &  tune_box_liq_asy            ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Tuning factor for box_liq reduction near the surface
    & tune_box_liq_sfc_fac(max_dom)! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Tuning factor for steeper dependence CLC(RH). This is an unphysical ad-hoc
    & allow_overcast               ! parameter to improve the cloud cover in the Mediterranean.
                                   ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Scaling factor for subgrid-scale contribution to diagnosed cloud ice
    &  tune_sgsclifac              ! (in case of inwp_cldcover = 1)

  INTEGER :: &                     !< Mode of coupling between turbulence and cloud cover
    &  icpl_turb_clc               ! 1: strong dependency of box width on rcld with upper and lower limit
                                   ! 2: weak dependency of box width on rcld with additive term and upper limit

  REAL(wp) :: &                    !< Tuning factor for enhanced LW absorption of mineral dust in the Saharan region
    &  tune_dust_abs               !

  REAL(wp) :: &                    !< Tuning factor for 3D contribution to diagnosed diffuse radiation
    &  tune_difrad_3dcont          !

  REAL(wp) :: &                    !< Tuning factor for gust parameterization
    &  tune_gust_factor            !

  REAL(wp) :: &                    !< Basic gust speed (m/s) at which the SSO correction starts to be reduced
    &  tune_gustsso_lim            !

  INTEGER :: &                     !< (MODIS) albedo tuning
    &  itune_albedo                ! 0: no tuning

  REAL(wp):: &                     !< bare soil albedo correction for soil types 3-6
    &  tune_albedo_wso(2)          ! tune_albedo_wso(1): albedo correction added over dry soil (w_so(1) < 0.001 m)
                                   ! tune_albedo_wso(2): albedo correction added over wet soil (w_so(1) > 0.002 m)

  INTEGER :: &                     !< slope-dependent tuning of parameters affecting stable PBLs
    &  itune_slopecorr             ! 1: slope-dependent reduction of rlam_heat and near-surface tkhmin

  INTEGER :: &                     !< type of artificial ozone tuning 
    &  itune_o3                    ! 0: no tuning
                                   ! 1: old tuning for RRTM radiation
                                   ! 2: (default) standard tuning for EcRad with RRTM gas optics
                                   ! 3: improved (for middle/upper stratosphere) tuning for EcRad with RRTM gas optics
                                   ! 4: provisional tuning for EcRad with EcCKD gas optics

  INTEGER :: &                     !< Type of gust tuning / SSO coupling
    &  itune_gust_diag             ! 1: use level above top of SSO envelope layer
                                   ! 2: use envelope top level, combined with adjusted tuning for MERIT/REMA orography
                                   ! 3: tuning for ICON-D2 with subgrid-scale condensation
                                   ! 4: tuning for ICON-D05 using time-averaged 10-m wind speed and additional limitations

  REAL(wp) :: &                    !< Height above ground up to which gust limitation is computed
    &  tune_gustlim_agl(max_dom)   !

  REAL(wp) :: &                    !< Tuning factor for gust limitation
    &  tune_gustlim_fac(max_dom)   !

  LOGICAL :: &                     ! cloud cover calibration over land points
    &  lcalib_clcov

  REAL(wp) :: &                    !< maximum calibration factor for low cloud cover (CLCL)
    &  max_calibfac_clcl

  REAL(wp) :: &                    !< maximum allowed positive freshsnow increment
       &  max_freshsnow_inc
  
  REAL(wp) :: &                    !< critical threshold for lower tropospheric stability (K)
       &  tune_eiscrit             !< to switch off conv param in stratocumulus regions

  REAL(wp) :: &                    !< critical threshold for lower tropospheric stability (K)
       &  tune_sc_eis              !< used for enhanced stratocumulus cloud cover

  REAL(wp) :: &                    !< minimum inversion height (m) used to define region with
       &  tune_sc_invmin           !< enhanced stratocumulus cloud cover

  REAL(wp) :: &                    !< maximum inversion height (m) used to define region with
       &  tune_sc_invmax           !< enhanced stratocumulus cloud cover

  REAL(wp) :: &                    !< scaling of direct solar rediation to tune sunshine duration
       &  tune_dursun_scaling      !< in corresponding diagnostic

  REAL(wp) :: &                    !< tuning of anthropogenic heat flux
       &  tune_urbahf(4)

  REAL(wp) :: &                    !< lower and upper bound for variable ISA paraeterization 
       &  tune_urbisa(2)           !< depending on smoothed urban fraction
  
  
  NAMELIST/nwp_tuning_nml/ tune_gkwake, tune_gkdrag, tune_gfluxlaun, tune_gcstar, &
    &                      tune_zceff_min, tune_v0snow, tune_zvz0i,               &
    &                      tune_entrorg, itune_albedo, tune_albedo_wso,           &
    &                      max_freshsnow_inc,                                     &
    &                      tune_capdcfac_et, tune_box_liq, tune_rhebc_land,       &
    &                      tune_rhebc_ocean, tune_rcucov, tune_texc,              &
    &                      tune_qexc, tune_minsnowfrac,tune_rhebc_land_trop,      &
    &                      tune_rhebc_ocean_trop, tune_rcucov_trop,               &
    &                      tune_dust_abs, tune_gfrcrit, tune_grcrit,              &
    &                      lcalib_clcov, tune_box_liq_asy, tune_capdcfac_tr,      &
    &                      tune_icesedi_exp, tune_rprcon, tune_gust_factor,       &
    &                      tune_rdepths, tune_thicklayfac, tune_sgsclifac,        &
    &                      icpl_turb_clc, tune_difrad_3dcont, max_calibfac_clcl,  &
    &                      tune_box_liq_sfc_fac, allow_overcast, tune_minsso,     &
    &                      tune_blockred, itune_gust_diag, tune_rcapqadv,         &
    &                      tune_gustsso_lim, tune_eiscrit, itune_o3,              &
    &                      tune_sc_eis, tune_sc_invmin, tune_sc_invmax,           &
    &                      tune_capethresh, tune_gkdrag_enh, tune_grcrit_enh,     &
    &                      tune_minsso_gwd, tune_dursun_scaling, tune_sbmccn,     &
    &                      itune_slopecorr, tune_gustlim_agl, tune_gustlim_fac,   &
    &                      tune_urbahf, tune_urbisa

CONTAINS


  !-------------------------------------------------------------------------
  !
  !! Read Namelist for NWP physics tuning. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP physics tuning
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  SUBROUTINE read_nwp_tuning_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit, jg

    REAL(wp) :: gkwake_def, gkdrag_def, gfrcrit_def, grcrit_def, minsso_def, minsso_gwd_def, blockred_def

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nwp_tuning_nml: read_tuning_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    ! Comment: In case we want to draw from a normal distribution, the namelist 
    ! parameters could be extended to arrays of size 2. The first value is the mean, 
    ! while the second one is the standard deviation. 

    ! SSO tuning
    gkwake_def  = 1.5_wp       ! original COSMO value 0.5
    gkdrag_def  = 0.075_wp     ! original COSMO value 0.075
    gfrcrit_def = 0.4_wp       ! original COSMO value 0.5
    grcrit_def  = 0.25_wp      ! original COSMO value 0.25
    minsso_def  = 10._wp       ! default 10 m (hardcoded value in original scheme)
    minsso_gwd_def = 0._wp     ! do not use enhanced threshold for activating GWD calculation
    blockred_def = 100._wp     ! effectively deactivates the blocking reduction

    tune_gkwake(:)  = gkwake_def
    tune_gkdrag(:)  = gkdrag_def
    tune_gkdrag_enh(:)  = gkdrag_def
    tune_gfrcrit(:) = gfrcrit_def
    tune_grcrit(:)  = grcrit_def
    tune_grcrit_enh(:)  = grcrit_def
    tune_minsso(:)  = minsso_def
    tune_minsso_gwd(:)  = minsso_gwd_def
    tune_blockred(:) = blockred_def

    !
    ! GWD tuning
    tune_gfluxlaun  = 2.50e-3_wp   ! original IFS value 3.75e-3
    tune_gcstar     = 1.0_wp
    !
    ! grid scale microphysics
    tune_zceff_min   = 0.01_wp
    tune_v0snow      = -1.0_wp      ! defaults are set in data_gscp depending on igscp
    tune_zvz0i       = 1.25_wp      ! original value of Heymsfield+Donner 1990: 3.29
    tune_icesedi_exp = 0.30_wp      ! exponent for density correction of cloud ice sedimentation
    tune_sbmccn      = 1.0_wp       ! [0-1] scaling factor to reduce the ccn concentration initial profile with respect to the polluted case

    !
    ! convection
    !

    !> entrainment parameter for deep convection:
    tune_entrorg     = 1.95e-3_wp   

    !> coefficient for conversion of cloud water into precipitation
    tune_rprcon      = 1.4e-3_wp

    !> maximum shallow convection depth (Pa)
    tune_rdepths     = 2.e4_wp

    !> fraction of CAPE diurnal cycle correction applied in the extratropics
    tune_capdcfac_et = 0.5_wp

    !> fraction of CAPE diurnal cycle correction applied in the tropics
    tune_capdcfac_tr = 0.5_wp

    !> CAPE threshold above which the convective adjustment time scale and entrainment
    !  are reduced for numerical stability [J/kg]
    tune_capethresh  = 7000._wp

    !> RH threshold for onset of evaporation below cloud base over land (original IFS value 0.7)
    tune_rhebc_land  = 0.75_wp

    !> RH threshold for onset of evaporation below cloud base over sea (original IFS value 0.9)
    tune_rhebc_ocean = 0.85_wp

    !> Convective area fraction used for computing evaporation below cloud base (original IFS value 0.05)
    tune_rcucov      = 0.05_wp

    !> Excess value for temperature used in test parcel ascent (K) (original IFS value 0.2 K)
    tune_texc        = 0.125_wp

    !> Excess fraction of grid-scale QV used in test parcel ascent
    !  (original IFS value 0.1 g/kg independent of grid-scale QV))
    tune_qexc        = 1.25e-2_wp

    ! Factor for dynamic correction of cape closure (turned off by default for the time being)
    tune_rcapqadv    = 0.0_wp

    ! The following switches allow separate tuning for evaporation below cloud base in the tropics
    tune_rhebc_land_trop  = 0.70_wp
    tune_rhebc_ocean_trop = 0.76_wp
    tune_rcucov_trop      = 0.03_wp

    !
    ! snow cover diagnosis
    tune_minsnowfrac = 0.20_wp     ! Minimum value to which the snow cover fraction is artificially reduced
                                   ! in case of melting show (in case of idiag_snowfrac = 20)
    !
    ! cloud cover
    tune_box_liq     = 0.05_wp     ! box width scale of liquid clouds
    tune_thicklayfac = 0.005_wp    ! factor [1/m] for increasing the box with for layer thicknesses exceeding 150 m
    tune_box_liq_asy = 3._wp       ! asymmetry factor for liquid cloud parameterization
    tune_box_liq_sfc_fac(:) = 1._wp   ! Tuning factor for box_liq reduction near the surface
    allow_overcast   = 1._wp       ! Tuning factor for steeper dependence CLC(RH)
    tune_sgsclifac   = 0._wp       ! Scaling factor for subgrid-scale contribution to diagnosed cloud ice
    lcalib_clcov     = .TRUE.      ! use calibration of layer-wise cloud cover diagnostics over land
    max_calibfac_clcl = 4._wp      ! maximum calibration factor for low cloud cover (CLCL); affects diagnostics only
    icpl_turb_clc    = 1           ! use strong dependency of box with on rcld (with factor 4) and upper and lower limit

    tune_gust_factor = 8.0_wp      ! tuning factor for gust parameterization
    itune_gust_diag  = 1           ! variant using level above SSO envelope
    tune_gustsso_lim = 100._wp     ! Basic gust speed at which the SSO correction starts to be reduced
    tune_gustlim_agl(:) = 1500._wp ! AGL height used for gust limitation in case of itune_gust_diag=4
    tune_gustlim_fac(:) = 0.0_wp   ! Corresponding tuning factor (0 means that limiting is deactivated)

    tune_dust_abs   = 0._wp        ! no tuning of LW absorption of mineral dust
    tune_difrad_3dcont = 0.5_wp    ! tuning factor for 3D contribution to diagnosed diffuse radiation (no impact on prognostic results!)
    itune_albedo    = 0            ! original (measured) albedo
    tune_albedo_wso = (/0._wp, 0._wp/) ! no bare soil albedo correction for soil types 3-6 (dry soil, wet soil)
    itune_o3        = 2            ! standard ozone tuning for EcRad
    itune_slopecorr  = 0           ! slope-dependent reduction of rlam_heat and near-surface tkhmin
    !
    ! IAU increment tuning
    max_freshsnow_inc = 0.025_wp   ! maximum allowed positive freshsnow increment

    !> critical stability threshold - if exceeded (stratocumulus regime)
    !> conv param is switched off to allow grid scale microphysics to create cloud
    !> If the default value of 1000 is kept, the atmosphere's EIS should never
    !> exceed the critical value, and this option remains effectively switched off
    tune_eiscrit     = 1000.0_wp

    !> critical stability threshold - used to identify region with stratocumulus
    !> cloud together with inversion height criterion, in order to enhance
    !> the diagnostic cloud cover in this region. The Sc cloud cover enhancement
    !> can be switched on/off independently of the EIS-based decision to use
    !> shallow convection in this region, hence two tuning parameters. With the
    !> default value of 1000 the enhancement is effectively turned off.
    tune_sc_eis      = 1000.0_wp

    !> min/max inversion height (m) to identify stratocumulus region for cloud cover
    !> enhancement. 
    tune_sc_invmin   = 200._wp
    tune_sc_invmax   = 1500._wp

    !> scaling of direct solar radiation in sunshine duration diagnostic
    tune_dursun_scaling = 1._wp

    !> Tuning parameters for TERRA-URB
    tune_urbahf = (/0._wp,2._wp,2._wp,50._wp/)   ! anthropogenic heat flux; base values and gradients for heating and cooling; upper limit
    tune_urbisa = (/0.6_wp,1._wp/)        ! lower and upper bound for variable ISA parameterization depending on smoothed urban fraction

    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, nwp_tuning_nml)   ! write defaults to temporary text file
    END IF

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('nwp_tuning_nml')
      READ(funit,NML=nwp_tuning_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('nwp_tuning_nml', STATUS=istat)

    SELECT CASE (istat)
    CASE (POSITIONED)

      ! Set array parameters to dummy values to determine which ones are actively set in the namelist
      tune_gkwake(:)  = -1._wp
      tune_gkdrag(:)  = -1._wp
      tune_gkdrag_enh(:)  = -1._wp
      tune_gfrcrit(:) = -1._wp
      tune_grcrit(:)  = -1._wp
      tune_grcrit_enh(:)  = -1._wp
      tune_minsso(:)  = -1._wp
      tune_minsso_gwd(:)  = -1._wp
      tune_blockred(:) = -1._wp

      READ (nnml, nwp_tuning_nml)    ! overwrite default settings

      ! Reset first element to default values if not specified in the namelist
      IF (tune_gkwake(1)  < 0._wp) tune_gkwake(1)  = gkwake_def
      IF (tune_gkdrag(1)  < 0._wp) tune_gkdrag(1)  = gkdrag_def
      IF (tune_gfrcrit(1) < 0._wp) tune_gfrcrit(1) = gfrcrit_def
      IF (tune_grcrit(1)  < 0._wp) tune_grcrit(1)  = grcrit_def
      IF (tune_minsso(1)  < 0._wp) tune_minsso(1)  = minsso_def
      IF (tune_minsso_gwd(1) < 0._wp) tune_minsso_gwd(1) = minsso_gwd_def
      IF (tune_blockred(1) < 0._wp) tune_blockred(1) = blockred_def

      ! Fill remaining array elements with entry of parent domain if not specified in the namelist
      DO jg = 2, max_dom
        IF (tune_gkwake(jg)  < 0._wp) tune_gkwake(jg)  = tune_gkwake(jg-1)
        IF (tune_gkdrag(jg)  < 0._wp) tune_gkdrag(jg)  = tune_gkdrag(jg-1)
        IF (tune_gfrcrit(jg) < 0._wp) tune_gfrcrit(jg) = tune_gfrcrit(jg-1)
        IF (tune_grcrit(jg)  < 0._wp) tune_grcrit(jg)  = tune_grcrit(jg-1)
        IF (tune_minsso(jg)  < 0._wp) tune_minsso(jg)  = tune_minsso(jg-1)
        IF (tune_minsso_gwd(jg)  < 0._wp) tune_minsso_gwd(jg)  = tune_minsso_gwd(jg-1)
        IF (tune_blockred(jg) < 0._wp) tune_blockred(jg) = tune_blockred(jg-1)
      ENDDO

      DO jg = 1, max_dom
        IF (tune_gkdrag_enh(jg)  < 0._wp) tune_gkdrag_enh(jg)  = tune_gkdrag(jg)
        IF (tune_grcrit_enh(jg)  < 0._wp) tune_grcrit_enh(jg)  = tune_grcrit(jg)
      ENDDO

      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, nwp_tuning_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------

    IF ( ANY (ABS(tune_albedo_wso) >= 0.3_wp)) THEN
      CALL finish(TRIM(routine), 'albedo correction tune_albedo_wso must be in the range -0.3 to 0.3')
    ENDIF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    config_tune_gkwake           = tune_gkwake
    config_tune_gkdrag           = tune_gkdrag
    config_tune_gkdrag_enh       = tune_gkdrag_enh
    config_tune_gfrcrit          = tune_gfrcrit
    config_tune_grcrit           = tune_grcrit
    config_tune_grcrit_enh       = tune_grcrit_enh
    config_tune_minsso           = tune_minsso
    config_tune_minsso_gwd       = tune_minsso_gwd
    config_tune_blockred         = tune_blockred
    config_tune_gfluxlaun        = tune_gfluxlaun
    config_tune_gcstar           = tune_gcstar
    config_tune_zceff_min        = tune_zceff_min
    config_tune_v0snow           = tune_v0snow
    config_tune_zvz0i            = tune_zvz0i
    config_tune_icesedi_exp      = tune_icesedi_exp
    config_tune_entrorg          = tune_entrorg
    config_tune_rprcon           = tune_rprcon
    config_tune_rdepths          = tune_rdepths
    config_tune_capdcfac_et      = tune_capdcfac_et
    config_tune_capdcfac_tr      = tune_capdcfac_tr
    config_tune_capethresh       = tune_capethresh
    config_tune_rhebc_land       = tune_rhebc_land
    config_tune_rhebc_ocean      = tune_rhebc_ocean
    config_tune_rcucov           = tune_rcucov
    config_tune_rhebc_land_trop  = tune_rhebc_land_trop
    config_tune_rhebc_ocean_trop = tune_rhebc_ocean_trop
    config_tune_rcucov_trop      = tune_rcucov_trop
    config_tune_texc             = tune_texc
    config_tune_qexc             = tune_qexc
    config_tune_rcapqadv         = tune_rcapqadv
    config_tune_minsnowfrac      = tune_minsnowfrac
    config_tune_box_liq          = tune_box_liq
    config_tune_box_liq_asy      = tune_box_liq_asy
    config_tune_box_liq_sfc_fac  = tune_box_liq_sfc_fac
    config_allow_overcast        = allow_overcast
    config_tune_thicklayfac      = tune_thicklayfac
    config_tune_sgsclifac        = tune_sgsclifac
    config_icpl_turb_clc         = icpl_turb_clc
    config_tune_dust_abs         = tune_dust_abs
    config_tune_difrad_3dcont    = tune_difrad_3dcont
    config_tune_gust_factor      = tune_gust_factor
    config_itune_gust_diag       = itune_gust_diag
    config_tune_gustsso_lim      = tune_gustsso_lim
    config_tune_gustlim_agl      = tune_gustlim_agl
    config_tune_gustlim_fac      = tune_gustlim_fac
    config_itune_albedo          = itune_albedo
    config_tune_albedo_wso       = tune_albedo_wso    
    config_itune_slopecorr       = itune_slopecorr
    config_itune_o3              = itune_o3
    config_lcalib_clcov          = lcalib_clcov
    config_max_calibfac_clcl     = max_calibfac_clcl
    config_max_freshsnow_inc     = max_freshsnow_inc
    config_tune_eiscrit          = tune_eiscrit
    config_tune_sc_eis           = tune_sc_eis
    config_tune_sc_invmin        = tune_sc_invmin
    config_tune_sc_invmax        = tune_sc_invmax
    config_tune_dursun_scaling   = tune_dursun_scaling
    config_tune_sbmccn           = tune_sbmccn
    config_tune_urbisa           = tune_urbisa
    config_tune_urbahf           = tune_urbahf

    !$ACC UPDATE DEVICE(config_tune_gust_factor, config_itune_gust_diag, config_tune_gustsso_lim) ASYNC(1)
    !$ACC UPDATE DEVICE(config_tune_gustlim_agl, config_tune_gustlim_fac, config_tune_albedo_wso) ASYNC(1)

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=nwp_tuning_nml)                    
      CALL store_and_close_namelist(funit, 'nwp_tuning_nml')             
    ENDIF

    !--------------------------------------------------------
    ! 7. write the contents of the namelist to an ASCII file
    !--------------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nwp_tuning_nml)


  END SUBROUTINE read_nwp_tuning_namelist

END MODULE mo_nwp_tuning_nml
