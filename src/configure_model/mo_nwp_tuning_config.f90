! @brief Tuning and/or perturbing nwp physics
!
! configuration setup for NWP physics tuning
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

MODULE mo_nwp_tuning_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tune_gkwake
  PUBLIC :: tune_gkdrag, tune_gkdrag_enh
  PUBLIC :: tune_gfrcrit
  PUBLIC :: tune_grcrit, tune_grcrit_enh
  PUBLIC :: tune_minsso, tune_minsso_gwd
  PUBLIC :: tune_blockred
  PUBLIC :: tune_gfluxlaun
  PUBLIC :: tune_gcstar
  PUBLIC :: tune_zceff_min
  PUBLIC :: tune_v0snow
  PUBLIC :: tune_zvz0i
  PUBLIC :: tune_icesedi_exp
  PUBLIC :: tune_entrorg
  PUBLIC :: tune_rprcon
  PUBLIC :: tune_rdepths
  PUBLIC :: tune_capdcfac_et
  PUBLIC :: tune_capdcfac_tr
  PUBLIC :: tune_capethresh
  PUBLIC :: tune_lowcapefac
  PUBLIC :: limit_negpblcape
  PUBLIC :: tune_rhebc_land
  PUBLIC :: tune_rhebc_ocean
  PUBLIC :: tune_rcucov
  PUBLIC :: tune_rhebc_land_trop
  PUBLIC :: tune_rhebc_ocean_trop
  PUBLIC :: tune_rcucov_trop
  PUBLIC :: tune_texc
  PUBLIC :: tune_qexc
  PUBLIC :: tune_rcapqadv
  PUBLIC :: tune_minsnowfrac
  PUBLIC :: tune_box_liq
  PUBLIC :: tune_box_liq_asy
  PUBLIC :: tune_box_liq_sfc_fac
  PUBLIC :: allow_overcast
  PUBLIC :: tune_thicklayfac
  PUBLIC :: tune_sgsclifac
  PUBLIC :: icpl_turb_clc
  PUBLIC :: tune_dust_abs
  PUBLIC :: tune_difrad_3dcont
  PUBLIC :: tune_gust_factor
  PUBLIC :: tune_gustsso_lim
  PUBLIC :: itune_gust_diag
  PUBLIC :: itune_albedo
  PUBLIC :: itune_slopecorr
  PUBLIC :: itune_o3
  PUBLIC :: lcalib_clcov
  PUBLIC :: max_calibfac_clcl
  PUBLIC :: max_freshsnow_inc
  PUBLIC :: tune_eiscrit
  PUBLIC :: tune_sc_eis  
  PUBLIC :: tune_sc_invmin
  PUBLIC :: tune_sc_invmax
  PUBLIC :: tune_dursun_scaling
  PUBLIC :: tune_sbmccn
  
  !!--------------------------------------------------------------------------
  !! Basic configuration setup for physics tuning
  !!--------------------------------------------------------------------------

!  TYPE :: t_nwp_tuning_config

    ! namelist variables
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

  REAL(wp) :: &                    !< [0-1] scaling factor to reduce the ccn concentration initial profile with respect to the polluted case
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

  REAL(wp) :: &                    !< Tuning factor for reducing the diurnal cycle correction in low-cape situations
    &  tune_lowcapefac = 1._wp     ! (relevant only if icapdcycl = 3; not a namelist variable)

  REAL(wp) :: &                    !< Minimum allowed negative PBL cape in diurnal cycle correction
    &  limit_negpblcape = 0._wp    ! (relevant only if icapdcycl = 3; not a namelist variable)

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

  REAL(wp) :: &                    !< Tuning factor for steeper dependence CLC(RH)
    & allow_overcast               ! (in case of inwp_cldcover = 1)

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
  !$ACC DECLARE CREATE(tune_gust_factor)

  INTEGER :: &                     !< Type of gust tuning / SSO coupling
    &  itune_gust_diag             ! 1: use level above top of SSO envelope layer
                                   ! 2: use envelope top level, combined with adjusted tuning
  !$ACC DECLARE CREATE(itune_gust_diag)

  REAL(wp) :: &                    !< Basic gust speed (m/s) at which the SSO correction starts to be reduced
    &  tune_gustsso_lim            !
  !$ACC DECLARE CREATE(tune_gustsso_lim)

  INTEGER :: &                     !< (MODIS) albedo tuning
    &  itune_albedo                ! 1: dimmed Sahara
                                   ! 2: dimmed Sahara and brighter Antarctica

  INTEGER :: &                     !< slope-dependent tuning of parameters affecting stable PBLs
    &  itune_slopecorr             ! 1: slope-dependent reduction of rlam_heat and near-surface tkhmin

  INTEGER :: &                     !< type of artificial ozone tuning 
    &  itune_o3                    ! 0: no tuning
                                   ! 1: old tuning for RRTM radiation
                                   ! 2: (default) standard tuning for EcRad with RRTM gas optics
                                   ! 3: improved (for middle/upper stratosphere) tuning for EcRad with RRTM gas optics
                                   ! 4: provisional tuning for EcRad with EcCKD gas optics

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
  
!  END TYPE t_nwp_tuning_config


!CONTAINS


END MODULE mo_nwp_tuning_config
