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

! Soil Vegetation Atmosphere Transfer (SVAT) scheme TERRA
! "Nihil in TERRA sine causa fit." (Cicero)!!

MODULE sfc_terra

!===============================================================================
!
! @par Description:
!   The module "sfc_terra.f90" performs calculations related to the
!   parameterization of soil processes. It contains the subroutine terra, which
!   is the combination of the former parts terra1 and terra2 of the COSMO-Model.
!
!   All parametric scalar and array data for this soil model routine are
!   defined in the data module sfc_terra_data.f90.
!   All global fields that are used by the soil model are passed through the
!   argument list.
!   All global scalar variables of the model that are used by the soil model
!   routine terra are imported by USE statements below.
!
! @par reference   This is an adaptation of subroutine terra_multlay in file
!  src_soil_multlay.f90 of the lm_f90 library (COSMO code). Equation numbers refer to
!  Doms, Foerstner, Heise, Herzog, Raschendorfer, Schrodin, Reinhardt, Vogel
!    (September 2005): "A Description of the Nonhydrostatic Regional Model LM",

!------------------------------------------------------------------------------
!
! Modules used:
#ifdef _OPENMP
  USE omp_lib,            ONLY: omp_get_thread_num
#endif


!------------------------------------------------------------------------------

USE sfc_terra_data  ! All variables from this data module are used by
                    ! this module. These variables start with letter "c"

!------------------------------------------------------------------------------


USE mo_mpi,                ONLY : get_my_global_mpi_id
!
USE mo_parallel_config,    ONLY: nproma
USE mo_kind,               ONLY: wp
USE mo_math_constants    , ONLY: pi
!
USE mo_physical_constants, ONLY: t0_melt => tmelt,& ! absolute zero for temperature
                                 r_d   => rd    , & ! gas constant for dry air
                                 rvd_m_o=>vtmpc1, & ! r_v/r_d - 1
                                 o_m_rdv        , & ! 1 - r_d/r_v
                                 rdv            , & ! r_d / r_v
                                 lh_v  => alv   , & ! latent heat of vapourization
                                 lh_s  => als   , & ! latent heat of sublimation
                                 lh_f  => alf   , & ! latent heat of fusion
                                 cp_d  => cpd   , & ! specific heat of dry air at constant press
                                 g     => grav  , & ! acceleration due to gravity
                                 b3    => t3    , & ! triple point of water at 611hPa
                                 sigma => stbo  , & ! Boltzmann-constant
                                 rho_w => rhoh2o, & ! density of liquid water (kg/m^3)
                                 rdocp => rd_o_cpd  ! r_d / cp_d
!
USE mo_lookup_tables_constants,ONLY: &
                                 b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b2i   => c3ies , & !!               -- " --
                                 b4w   => c4les , & !!               -- " --
                                 b4i   => c4ies     !!               -- " --

!
USE mo_lnd_nwp_config,     ONLY: lmulti_snow, l2lay_rho_snow,     &
  &                              itype_trvg, itype_evsl,          &
  &                              itype_root, itype_heatcond,      &
  &                              itype_hydbound,                  &
  &                              itype_canopy, tau_skin,          &
  &                              lterra_urb, itype_ahf, itype_eisa, &
  &                              lstomata,                        &
  &                              max_toplaydepth, itype_interception, &
  &                              cwimax_ml
!
!
USE mo_exception,          ONLY: finish
USE mo_run_config,         ONLY: msg_level
USE mo_fortran_tools,      ONLY: set_acc_host_or_device
USE mo_lnd_nwp_config,     ONLY: lcuda_graph_lnd

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE

!------------------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: terra

!------------------------------------------------------------------------------
! Public variables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Parameters and variables which are global in this module
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!+ Computation of the first part of the soil parameterization scheme
!------------------------------------------------------------------------------

  SUBROUTINE terra (                 &
                  nvec             , & ! array dimensions
                  ivstart          , & ! start index for computations in the parallel program
                  ivend            , & ! end index for computations in the parallel program
                  iblock           , & ! number of block
                  ke_soil, ke_snow , &
                  ke_soil_hy       , & ! number of active soil moisture layers
                  zmls             , & ! processing soil level structure
                  icant            , & ! canopy type
                  nclass_gscp      , & ! number of hydrometeor classes of grid scale microphysics
                  dt               , & ! time step
!
                  soiltyp_subs     , & ! type of the soil (keys 0-9)                     --
! for TERRA_URB
                  urb_isa          , & ! impervious surface area fraction of the urban canopy ( - )
                  urb_ai           , & ! surface area index of the urban canopy               ( - )
                  urb_h_bld        , & ! building height                                      ( m )
                  urb_hcap         , & ! volumetric heat capacity of urban material      (J/m**3/K)
                  urb_hcon         , & ! thermal conductivity of urban material             (W/m/K)
                  ahf              , & ! anthropogenic heat flux                           (W/m**2)
!
                  plcov            , & ! fraction of plant cover                         --
                  rootdp           , & ! depth of the roots                            ( m )
                  sai              , & ! surface area index                              --
                  tai              , & ! transpiration area index                        --
                  laifac           , & ! ratio between current LAI and laimax            --
                  eai              , & ! earth area (evaporative surface area) index     --
                  skinc            , & ! skin conductivity                             ( W/m**2/K )
!
                  heatcond_fac     , & ! tuning factor for soil thermal conductivity
                  heatcap_fac      , & ! tuning factor for soil heat capacity
                  hydiffu_fac      , & ! tuning factor for hydraulic diffusivity
!
                  rsmin2d          , & ! minimum stomatal resistance                   ( s/m )
                  r_bsmin          , & ! minimum bare soil evap resistance             ( s/m )
                  z0               , & ! vegetation roughness length                   ( m   )
!
                  u                , & ! zonal wind speed                              ( m/s )
                  v                , & ! meridional wind speed                         ( m/s )
                  t                , & ! temperature                                   (  k  )
                  qv               , & ! specific water vapor content                  (kg/kg)
                  qc               , & ! specific liquid-water content                 (kg/kg)
                  qi               , & ! specific frozen-water content                 (kg/kg)
                  ptot             , & ! full pressure                                 ( Pa  )
                  ps               , & ! surface pressure                              ( Pa  )
!
                  t_snow_now       , & ! temperature of the snow-surface               (  K  )
                  t_snow_new       , & ! temperature of the snow-surface               (  K  )
!
                  t_snow_mult_now  , & ! temperature of the snow-surface               (  K  )
                  t_snow_mult_new  , & ! temperature of the snow-surface               (  K  )
!
                  t_s_now          , & ! temperature of the ground surface             (  K  )
                  t_s_new          , & ! temperature of the ground surface             (  K  )
!
                  t_sk_now         , & ! skin temperature                              (  K  )
                  t_sk_new         , & ! skin temperature                              (  K  )
!
                  t_g              , & ! weighted surface temperature                  (  K  )
                  qv_s             , & ! specific humidity at the surface              (kg/kg)
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  w_snow_new       , & ! water content of snow                         (m H2O)
!
                  rho_snow_now     , & ! snow density                                  (kg/m**3)
                  rho_snow_new     , & ! snow density                                  (kg/m**3)
!
                  rho_snow_mult_now, & ! snow density                                  (kg/m**3)
                  rho_snow_mult_new, & ! snow density                                  (kg/m**3)
!
                  h_snow           , & ! snow depth                                   (  m  )
                  h_snow_gp        , & ! grid-point averaged snow depth               (  m  )
                  meltrate         , & ! snow melting rate                             (kg/(m**2*s))
                  tsnred           , & ! snow temperature offset for calculating evaporation  (K)
!
                  w_i_now          , & ! water content of interception store           (m H2O)
                  w_i_new          , & ! water content of interception store           (m H2O)
!
                  w_p_now          , & ! water content of pond interception water     (m H2O)
                  w_p_new          , & ! water content of pond interception water     (m H2O)
!
                  w_s_now          , & ! water content of interception snow           (m H2O)
                  w_s_new          , & ! water content of interception snow           (m H2O)
!
                  t_so_now         , & ! soil temperature (main level)                 (  K  )
                  t_so_new         , & ! soil temperature (main level)                 (  K  )
!
                  w_so_now         , & ! total soil water content (ice + liquid water) (m H20)
                  w_so_new         , & ! total soil water content (ice + liquid water) (m H20)
!
                  w_so_ice_now     , & ! soil ice content                              (m H20)
                  w_so_ice_new     , & ! soil ice content                              (m H20)
!
!                 t_2m             , & ! temperature in 2m                             (  K  )
                  u_10m            , & ! zonal wind in 10m                             ( m/s )
                  v_10m            , & ! meridional wind in 10m                        ( m/s )
                  freshsnow        , & ! indicator for age of snow in top of snow layer(  -  )
                  zf_snow          , & ! snow-cover fraction                           (  -  )
!
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wliq_snow_new    , & ! liquid water content in the snow              (m H2O)
!
                  wtot_snow_now    , & ! total (liquid + solid) water content of snow  (m H2O)
                  wtot_snow_new    , & ! total (liquid + solid) water content of snow  (m H2O)
!
                  dzh_snow_now     , & ! layer thickness between half levels in snow   (  m  )
                  dzh_snow_new     , & ! layer thickness between half levels in snow   (  m  )
!
                  prr_con          , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con          , & ! precipitation rate of snow, convective        (kg/m2*s)
                  conv_frac        , & ! convective area fraction as assumed in convection scheme
                  prr_gsp          , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp          , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                  pri_gsp          , & ! precipitation rate of ice, grid-scale        (kg/m2*s)
                  prg_gsp          , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
#ifdef TWOMOM_SB
                  prh_gsp          , & ! precipitation rate of hail, grid-scale        (kg/m2*s)
#endif
!
                  tch              , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm              , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv              , & ! laminar reduction factor for evaporation      ( -- )
                  tfvsn            , & ! reduction factor for snow evaporation from model-DA coupling     ( -- )
!
                  sobs             , & ! solar radiation at the ground                 ( W/m2)
                  thbs             , & ! thermal radiation at the ground               ( W/m2)
                  pabs             , & !!!! photosynthetic active radiation            ( W/m2)
!
                  runoff_s         , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g         , & ! soil water runoff; sum over forecast          (kg/m2)
                  resid_wso        , & ! soil water budget, residuum                   (kg/m2)
!
                  zshfl_s          , & ! sensible heat flux soil/air interface         (W/m2)
                  zlhfl_s          , & ! latent   heat flux soil/air interface         (W/m2)
                  zshfl_snow       , & ! sensible heat flux snow/air interface         (W/m2)
                  zlhfl_snow       , & ! latent   heat flux snow/air interface         (W/m2)
                  lhfl_bs          , & ! latent heat flux from bare soil evap.         (W/m2)
                  lhfl_pl          , & ! latent heat flux from plants                  (W/m2)
                  plevap           , & ! function of accumulated plant evaporation     (kg/m2)
                  rstom            , & ! stomatal resistance                           ( s/m )
                  zshfl_sfc        , & ! sensible heat flux surface interface          (W/m2)
                  zlhfl_sfc        , & ! latent   heat flux surface interface          (W/m2)
                  zqhfl_sfc        , & ! water vapor   flux surface interface          (kg/m2/s)
                  ldiff_qi         , & ! turbulent diffusion of frozen water is active
                  ldepo_qw         , & ! deposition of (frozen or liquid) cloud water required
                  lres_soilwatb    , & ! flag for computing the soil water budget
                  lacc             , & ! flag for activating OpenACC
                  opt_acc_async_queue) ! OpenACC stream
                                     

!-------------------------------------------------------------------------------
! Declarations
!-------------------------------------------------------------------------------


  INTEGER, INTENT(IN)  ::  &
                  icant,             & ! canopy type
                  nvec,              & ! array dimensions
                  ivstart,           & ! start index for computations in the parallel program
                  ivend,             & ! end index for computations in the parallel program
                  iblock,            & ! number of block
                  ke_soil, ke_snow,  &
                  ke_soil_hy           ! number of active soil moisture layers
  REAL    (KIND = wp), DIMENSION(ke_soil+1), INTENT(IN) :: &
                  zmls                 ! processing soil level structure
  INTEGER, INTENT(IN)  :: &
                  nclass_gscp          ! number of hydrometeor classes of grid scale microphysics
  REAL    (KIND = wp), INTENT(IN)  ::  &
                  dt                   ! time step

  INTEGER, DIMENSION(nvec), INTENT(IN) :: &
                  soiltyp_subs         ! type of the soil (keys 0-9)                     --

  REAL    (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
! for TERRA_URB
                  urb_isa          , & ! impervious surface area fraction of the urban canopy ( - )
                  urb_ai           , & ! surface area index of the urban canopy               ( - )
                  urb_h_bld        , & ! building height                                      ( m )
                  urb_hcap         , & ! volumetric heat capacity of urban material      (J/m**3/K)
                  urb_hcon         , & ! thermal conductivity of urban material             (W/m/K)
                  ahf              , & ! anthropogenic heat flux                           (W/m**2)
!
                  plcov            , & ! fraction of plant cover                         --
                  rootdp           , & ! depth of the roots                            ( m  )
                  sai              , & ! surface area index                              --
                  tai              , & ! transpiration area index                        --
                  laifac           , & ! ratio between current LAI and laimax
                  eai              , & ! earth area (evaporative surface area) index     --
                  skinc            , & ! skin conductivity                        ( W/m**2/K )
!
                  heatcond_fac     , & ! tuning factor for soil thermal conductivity
                  heatcap_fac      , & ! tuning factor for soil heat capacity
                  hydiffu_fac      , & ! tuning factor for hydraulic diffusivity
!
                  rsmin2d          , & ! minimum stomata resistance                    ( s/m )
                  r_bsmin          , & ! minimum bare soil evap resistance             ( s/m )
                  u                , & ! zonal wind speed                              ( m/s )
                  v                , & ! meridional wind speed                         ( m/s )
                  t                , & ! temperature                                   (  k  )
                  qv               , & ! specific water vapor content                  (kg/kg)
                  qc               , & ! specific liquid-water content                 (kg/kg)
                  qi               , & ! specific frozen-water content                 (kg/kg)
                  ptot             , & ! full pressure                                 ( Pa )
                  ps               , & ! surface pressure                              ( pa  )
                  h_snow_gp        , & ! grid-point averaged snow depth
                  tsnred           , & ! snow temperature offset for calculating evaporation (K)
                  u_10m            , & ! zonal wind in 10m                             ( m/s )
                  v_10m            , & ! meridional wind in 10m                        ( m/s )
                  prr_con          , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con          , & ! precipitation rate of snow, convective        (kg/m2*s)
                  conv_frac        , & ! convective area fraction
                  prr_gsp          , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp          , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                  pri_gsp          , & ! precipitation rate of ice, grid-scale         (kg/m2*s)
                  prg_gsp          , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
#ifdef TWOMOM_SB
                  prh_gsp          , & ! precipitation rate of hail, grid-scale        (kg/m2*s)
#endif
                  sobs             , & ! solar radiation at the ground                 ( W/m2)
                  thbs             , & ! thermal radiation at the ground               ( W/m2)
                  pabs                 !!!! photosynthetic active radiation            ( W/m2)

  REAL    (KIND = wp), DIMENSION(nvec), OPTIONAL, INTENT(IN) :: &
                  z0                   ! vegetation roughness length                    ( m )

  REAL    (KIND = wp), DIMENSION(nvec),           INTENT(IN) :: &
                  tfv              , & ! laminar reduction factor for evaporation      ( -- )
                  tfvsn                ! reduction factor for snow evaporation from model-DA coupling     ( -- )

  REAL    (KIND = wp), DIMENSION(nvec), OPTIONAL, INTENT(INOUT) :: &
                  plevap               ! function of accumulated plant evaporation     (kg/m2)

  REAL    (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
                  t_snow_now       , & ! temperature of the snow-surface (K)
                  t_s_now          , & ! temperature of the ground surface             (  K  )
                  t_sk_now         , & ! skin temperature                              (  K  )
                  t_g              , & ! weighted surface temperature                  (  K  )
                  qv_s             , & ! specific humidity at the surface              (kg/kg)
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  rho_snow_now     , & ! snow density                                  (kg/m**3)
                  h_snow           , & ! snow depth
                  w_i_now          , & ! water content of interception store           (m H2O)
                  w_p_now          , & ! water content of interception water pond      (m H2O)
                  w_s_now          , & ! water content of interception snow water      (m H2O)
                  freshsnow        , & ! indicator for age of snow in top of snow layer(  -  )
                  zf_snow          , & ! snow-cover fraction
                  tch              , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm              , & ! turbulent transfer coefficient for momentum   ( -- )
                  runoff_s         , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g         , & ! soil water runoff; sum over forecast          (kg/m2)
                  resid_wso            ! residuum of the budget of soil water content  (kg/m2)

  REAL    (KIND = wp), DIMENSION(nvec), INTENT(OUT) :: &
                  t_snow_new       , & !
                  t_s_new          , & ! temperature of the ground surface             (  K  )
                  t_sk_new         , & ! skin temperature                              (  K  )
                  w_snow_new       , & ! water content of snow                         (m H2O)
                  rho_snow_new     , & ! snow density                                  (kg/m**3)
                  meltrate         , & ! snow melting rate
                  w_i_new          , & ! water content of interception store           (m H2O)
                  w_p_new          , & ! water content of interception water pond      (m H2O)
                  w_s_new          , & ! water content of interception snow water
                  zshfl_s          , & ! sensible heat flux soil/air interface         (W/m2)
                  zlhfl_s          , & ! latent   heat flux soil/air interface         (W/m2)
                  zshfl_snow       , & ! sensible heat flux snow/air interface         (W/m2)
                  zlhfl_snow       , & ! latent   heat flux snow/air interface         (W/m2)
                  rstom            , & ! stomata resistance                            ( s/m )
                  lhfl_bs              ! latent heat flux from bare soil evap.         ( W/m2)


  REAL    (KIND = wp), DIMENSION(nvec,0:ke_snow), INTENT(INOUT) :: &
                  t_snow_mult_now      ! temperature of the snow-surface               (  K  )
  REAL    (KIND = wp), DIMENSION(nvec,0:ke_snow), INTENT(OUT) :: &
                  t_snow_mult_new      ! temperature of the snow-surface               (  K  )

  REAL    (KIND = wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
                  rho_snow_mult_now, & ! snow density                                  (kg/m**3)
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_now    , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow_now         ! layer thickness between half levels in snow   (  m  )

  REAL    (KIND = wp), DIMENSION(nvec,ke_snow), INTENT(OUT) :: &
                  rho_snow_mult_new, & ! snow density                                  (kg/m**3)
                  wliq_snow_new    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_new    , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow_new         ! layer thickness between half levels in snow   (  m  )

  REAL    (KIND = wp), DIMENSION(nvec,0:ke_soil+1), INTENT(INOUT) :: &
                  t_so_now             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = wp), DIMENSION(nvec,0:ke_soil+1), INTENT(OUT) :: &
                  t_so_new             ! soil temperature (main level)                 (  K  )

  REAL    (KIND = wp), DIMENSION(nvec,ke_soil+1), INTENT(INOUT) :: &
                  w_so_now         , & ! total soil water content (ice + liquid water) (m H20)
                  w_so_ice_now         ! soil ice content                              (m H20)
  REAL    (KIND = wp), DIMENSION(nvec,ke_soil+1), INTENT(OUT) :: &
                  w_so_new         , & ! total soil water content (ice + liquid water) (m H20)
                  w_so_ice_new         ! soil ice content                              (m H20)



!US why +1  REAL    (KIND = wp), DIMENSION(nvec,ke_soil+1), INTENT(OUT) :: &
!US is really +1 in ICON
  REAL    (KIND = wp), DIMENSION(nvec,ke_soil+1), INTENT(OUT) :: &
                  lhfl_pl          ! average latent heat flux from plants              ( W/m2)

  REAL    (KIND = wp), DIMENSION(nvec), OPTIONAL, INTENT(OUT) :: &
                  zshfl_sfc        , & ! sensible heat flux surface interface          (W/m2)
                  zlhfl_sfc        , & ! latent   heat flux surface interface          (W/m2)
                  zqhfl_sfc            ! latent   heat flux surface interface          (W/m2)

  LOGICAL, INTENT(IN) ::  &
                  ldiff_qi         , & ! turbulent diffusion of frozen water is active
                  ldepo_qw         , & ! deposition of (frozen or liquid) cloud water required
                  lres_soilwatb        ! calculation soil water budget desired

  LOGICAL, OPTIONAL, INTENT(IN) :: lacc
  INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

  LOGICAL :: lzacc
  INTEGER :: acc_async_queue

!--------------------------------------------------------------------------------
! TERRA Declarations
! ------------------

  INTEGER                  ::  &
    ierror                        ! error status variable

  CHARACTER (LEN=80)       ::  &
    yerror

! Local parameters:
! ----------------

  REAL(KIND=wp), PARAMETER ::  &
    !zepsi  = 1.0E-6_wp: is now used from sfc_terra_data: eps_soil, eps_temp, eps_div
    zalfa      =  1.0_wp, & ! degree of impliciteness (1: full implicit,
                            !    (0.5: Cranck-Nicholson)
    T_ref_ice  =  0.1_wp, & !degC Soil ice parameterization
    T_star_ice = 0.01_wp, & !degC according to K. Schaefer and Jafarov, E.,2016, doi:10.5194/bg-13-1991-2016
    b_clay     = -0.3_wp, &
    b_silt     = -0.5_wp, &
    b_sand     = -0.9_wp, &
    b_org      = -1.0_wp, &
  
    rho_i  = 910.0_wp       ! density of solid ice (soil model)  (kg/m**3)

! Local scalars:
! -------------

  INTEGER        ::  &

    ! Indices
    kso            , & ! loop index for soil moisture layers
    ksn            , & ! loop index for snow layers
    k              , & ! loop index for snow layers
    i,ic           , & ! loop index in x-direction
    jb             , & ! loop index for soil-type
    mstyp          , & ! soil type index
    msr_off        , & ! number of layers contributing to surface run off
    k10cm          , & ! index of half level closest to 0.1 m
    k100cm             ! index of half level closest to 1.0 m

  REAL(KIND=wp)  ::  &

    ! Timestep parameters
    zdt            , & ! integration time-step [s]
    zdtdrhw        , & ! timestep / density of water
    zrhwddt        , & ! density of water / timestep
    zroffdt        , & ! timestep for runoff calculation
    z1d2dt         , & ! 1./2*timestep

    ! Connection to the atmosphere
    zuv            , & ! wind velocity in lowest atmospheric layer
    ztvs           , & ! virtual temperature at the surface
    zplow          , & ! pressure of lowest atmospheric layer
    zqvlow         , & ! specific humidity of lowest atmospheric layer
    zrww           , & ! water covered part of grid element

    ! Snow parameters
    zsnow_rate     , & ! rate of snow fall            [kg/m**2 s]
    zrain_rate     , & ! rate of rain fall            [kg/m**2 s]
    zrime          , & ! ground riming rate
    zdsn_new       , & ! snow age refresh increment   [-]
    zdsn_old       , & ! snow age decay increment     [-]

    ! Multi-snow layer parameters
    zsn_porosity   , & ! snow porosity                    ( - )
    zp1            , &
    zfukt          , &
    zq0            , &
    zdzh_old       , &
    ze_in          , &
    zadd_dz        , &
    zeta           , &
    zdens_old      , &
    weight         , &

    ! Plant parameters
    zbeta          , & ! reduction factor for evaporation
    zevap          , & ! auxiliary variable
    zalpha             ! NP89 bare soil evaporation

  REAL(KIND=wp)  ::  &

    ! Local scalars for BATS-scheme
    zpsi0=0.01_wp  , & ! air entry potential at water saturation (m)
    zbf1           , & ! auxiliary variable
    zbf2           , & ! auxiliary variable
    zdmax          , & ! auxiliary variable
    zd             , & ! auxiliary variable
    zck            , & ! auxiliary variable
    z1             , & ! depth of half level closest to 0.1 m
    znull          , & ! depth of half level closest to 1.0 m
    zfqmax         , & ! maximum sustainable flux in uppermost soil layer
    zevapor        , & ! evaporation rate of bare soil in BATS-scheme
    zrla           , & ! atmospheric resistance
    zcatm          , & ! transfer function CA
    zpar           , & ! PAR (interim value)
    zf_wat         , & ! soil water function for stomatal resistance
    zf_tem         , & ! temperature function for stomatal resistance
    zf_sat         , & ! saturation deficit function for stomatal resistance
    zepsat         , & ! saturation vapour pressure at near surface temperature
    zepke          , & ! near surface water vapour pressure
    zrstom         , & ! stomata resistance
    zrveg          , & ! sum of atmospheric and stomata resistance
    zedrstom       , & ! 1./zrstom
    zustar         , & ! friction velocity
    ztrabpf        , & ! area average transpiration rate

    ! variables for non-uniform root-distribution
    ztrfr          , & ! fraction of total transpiration associated with each layer
    zrootfr        , & ! normalized root fraction in each layer (=1 for z_soil=0.0)
    zrootdz        , & ! normalized root fraction * min(layer thickness,rootlength-top_of_layer)

    ! Soil parameters
    zice           , & ! indicator of soil type ice
    zwqg           , & ! mean of fcap and pwp
    z4wdpv         , & ! 4*zwqg/porv
    zropart        , & ! fraction of layer filled by roots (Dickinson scheme)
    zrootfc        , & ! distributes transpiration to soil layers
    ze_sum             ! sum of all contributions to evapotranspiration

  REAL(KIND=wp)  ::  &

    ! Thermal parameters
    zgstr          , & ! downward longwave radiation
    zalas          , & ! heat conductivity of snow
    zrnet_snow     , & ! net radiation at snow surface
    zfak           , & ! utility variable for implicit snow temperature forecast
    ztsnow_im      , & ! utility variable for implicit snow temperature forecast
    ztsnew         , & ! preliminary value of snow surface temperature
    ztsnownew      , & ! preliminary value of snow surface temperature
    zfor_snow      , & ! total forcing at snow surface
    zfr_melt       , & ! melting snow fraction
    zdwgme         , & ! utility variable for snow melt determination
    zdelt_s        , & ! utility variable for snow melt determination
    ze_avail       , & ! utility variable for snow melt determination
    ze_total           ! utility variable for snow melt determination

  REAL    (KIND=wp) ::  &

    ! Hydraulic parameters
    zalf           , & ! utility variable
    zro_inf        , & ! surface runoff
    zdwsndtt       , & ! time rate of change of snow store
    zdwisndtt      , & ! time rate of change of snow store
    zdwpdtt        , & ! time rate of change of pond store
    zwsnstr        , & ! preliminary value of snow store
    zdwseps        , & ! artificial change of small amount of snow store
    zdwidtt        , & ! time rate of change of interception store
    zdwieps        , & ! artificial change of small amount of interception store
    zro_wi         , & ! surface runoff due to limited infiltration capacity
    zro_sfak       , & ! utility variable
    zro_gfak       , & ! utility variable
    zfmb_fak       , & ! utility variable
    zdwg           , & ! preliminary change of soil water content
    zwgn           , & ! preliminary change of soil water content
    zredfu         , & ! utility variable for runoff determination
    zro            , & ! utility variable for runoff determination
    zro2           , & ! utility variable for runoff determination
    zkorr          , & ! utility variable for runoff determination
    zfr_ice        , & ! reduction factor for water transport
    zfr_ice_free   , & ! reduction factor for water transport
    zwso_new       , & ! preliminary value of soil water content
    zw_ovpv        , & ! utility variable
    zfcorr_wi      , & ! utility variable
    zpercmax       , & ! utility variable

    ! Implicit solution of thermal and hydraulic equations
    zakb           , & ! utility variable
    zakb1          , & ! utility variable
    zakb2          , & ! utility variable
    zzz            , & ! utility variable
    z1dgam1        , & ! utility variable
    zredm          , & ! utility variable
    zredm05        , & ! utility variable
    zredp05        , & ! utility variable
    zgam2m05       , & ! utility variable
    zgam2p05           ! utility variable

  REAL    (KIND=wp) ::  &

    ! for FUNCTIONS below
    z2iw           , & ! dummy argument for Stmt. function
    z4iw           , & ! dummy argument for Stmt. function
    z234iw         , & ! dummy argument for Stmt. function
    zqs            , & ! saturation humidty at T_s and ps
    zdqs           , & ! derivative of zqs with respect to T_s
    zqsnow         , & ! saturation humidty at T_snow and ps
    zdqsnow        , & ! derivative of zqs with respect to T_snow

    !   Local scalars for transfer coefficient limitation calculations and for
    !   treatment of possible problems at lateral boundaries
    zch_soil  = 1.4E06_wp , & ! approximate average heat capacity of soil (J/m**3 K)
    zlim_dtdt = 2.5_wp    , & ! maximum allowed temperature increment (K) in
                              ! one time step in  uppermost soil layer
    zdT_snow                  ! maximum possible increment in snow layer before
                              ! melting sets in

  REAL    (KIND=wp) ::  &

    ! Local scalars for water content dependent freezing/melting
    zliquid        , & ! utility variable
    zxx            , & ! utility variable
    znen           , & ! utility variable
    zargu          , & ! utility variable

    ! Water transport
    zice_fr_ksom1  , & ! fractional ice content of actual layer - 1
    zice_fr_kso    , & ! fractional ice content of actual layer
    zice_fr_ksop1  , & ! fractional ice content of actual layer + 1
    zlw_fr_ksom1   , & ! fractional liquid water content of actual layer - 1
    zlw_fr_ksom1_new,& ! fractional liquid water content of actual layer -1
    zlw_fr_kso     , & ! fractional liquid water content of actual layer
    zlw_fr_kso_new , & ! fractional liquid water content of actual layer
    zlw_fr_ksop1   , & ! fractional liquid water content of actual layer + 1
    zlw_fr_ksom05  , & ! fractional liquid water content of actual layer - 1/2
    zdlw_fr_ksom05 , & ! hydraulic diffusivity coefficient at half level above
    zklw_fr_ksom05 , & ! hydraulic conductivity coefficient at half level above
    zklw_fr_kso_new, & ! hydraulic conductivity coefficient at main level
                       !    for actual runoff_g
    zdlw_fr_kso    , & ! hydraulic diff coefficient at main level
    zklw_fr_kso    , & ! hydraulic conductivity coefficient at main level
    zklw_fr_ksom1  , & ! hydraulic conductivity coefficient at main level above
    zlw_fr_ksop05  , & ! fractional liquid water content of actual layer + 1/2
    zdlw_fr_ksop05 , & ! hydraulic diffusivity coefficient at half level below
    zklw_fr_ksop05 , & ! hydraulic conductivity coefficient at half level below
    zinf           , & ! infiltration
    zcou_roffg     , & ! indicator to sum up runoffg in hydrological active layers

    ! Snow density
    ztau_snow      , & ! 'ageing constant' for snow density          (-)
    zrhosmax       , & ! temperature-dependent target density for snow ageing
    zrhosmin       , & ! wind-speed-dependent minimum density for new snow
    zgrfrac        , & ! fraction of graupel / convective snow
    zrho_snowe     , & ! updated density of existing snow ('ageing') kg/m**3
    zrho_snowf     , & ! density of fresh snow                       kg/m**3
    zrho_grauf     , & ! density of fresh graupel / convective snow  kg/m**3
    znorm          , & ! normalisation factor for weighted snow density mH2O

    ! Freezing/melting of soil water/ice:
    zdelwice       , & ! amount of melted soil ice/frozen soil water
    zdwi_scal      , & ! time scale parameter for freezing/melting soil water
    ztx            , & ! water content dependent freezing/melting temperature
    t_zw_up        , & ! temp -3 degC
    t_zw_low       , & ! temp -40 degC
    zd1, zd2, zd3, zd4 ! auxiliary variables

  REAL(KIND=wp)  ::  &
    ! for TERRA_URB
    zalpha_lnd     , & ! height-dependent conductivity factor for natural land
    zalpha_uf      , & ! height-dependent conductivity factor for buildings/street environment
    zisa_infil     , &
    ztalb          , & ! for modification of infrared albedo according to the building fraction
    ztemp

  INTEGER          :: m_limit        ! counter for application of limitation
  INTEGER          :: ke_soil_hy_m   ! number of active soil moisture layers for Mire
  CHARACTER(LEN=7) :: yhc            ! heating or cooling indicator

! Local (automatic) arrays:
! -------------------------

  INTEGER                        ::  &
    m_styp         (nvec)          , & ! soil type
    ke_soil_hy_b   (nvec)              ! lowest hydraulic active layer (array) 

  REAL(KIND = wp)                 :: &

    ! Two-time level variables exchange with interface
    h_snow_now     (nvec)          , & ! snow height  (m)
    h_snow_new     (nvec)          , & ! snow height  (m)

    ! Model geometry
    zdz_snow_fl    (nvec)          , & ! snow depth for snow temperature gradient

    ! Multi-layer snow model
    zhh_snow       (nvec,  ke_snow), & ! depth of the half level snow layers
    zhm_snow       (nvec,  ke_snow), & ! depth of snow main levels
    zdzh_snow      (nvec,  ke_snow), & ! layer thickness between half levels
    zdzm_snow      (nvec,  ke_snow), & ! distance between main levels
    zdtsnowdt_mult (nvec,0:ke_snow), & ! tendency of ztsnow

    ! External parameters
    zbwt           (nvec)          , & ! root depth (with artificial minimum value)
    zrock          (nvec)          , & ! ice/rock-indicator: 0 for ice and rock
    zsandf         (nvec)          , & ! mean fraction of sand (weight percent)
    zclayf         (nvec)          , & ! mean fraction of clay (weight percent)
    zsiltf         (nvec)          , & ! mean fraction of silt (weight percent)
    zb_por         (nvec)          , & ! pore size distribution index
    zpsis          (nvec)          , & ! air entry potential (m)
    zw_m_org       (nvec)          , & ! maximum of liquid water content organic
    zw_m_soil      (nvec)          , & ! maximum of liquid water content mineral soil
    zw_m_up        (nvec)          , & ! maximum of liquid water content at temp -3 degC
    zw_m_low       (nvec,ke_soil+1), & ! maximum of liquid water content at temp -40 degC
    zw_m           (nvec)              ! maximum of liquid water content  (m)

  REAL(KIND=wp)                  ::  &

    ! Connection to the atmosphere
    zrr            (nvec)          , & ! total rain rate including formation of dew and deposition of water droplets
    zrs            (nvec)          , & ! total snow rate including formation of rime and deposition of frozen particles
    zesoil         (nvec)          , & ! evaporation from bare soil
    zsfc_frac_bs   (nvec)          , & ! relative source surface of the bare soil
    zrhoch         (nvec)          , & ! transfer coefficient*rho*g
    zth_low        (nvec)          , & ! potential temperature of lowest layer
    zf_wi          (nvec)          , & ! surface fraction covered by interception water content
    ztmch          (nvec)          , & ! heat transfer coefficient*density*velocity
    zep_s          (nvec)          , & ! potential evaporation for t_s
    zep_snow       (nvec)          , & ! potential evaporation for t_snow
    zverbo         (nvec)          , & ! total evapotranspiration
    zversn         (nvec)          , & ! total evaporation from snow surface
    zthsoi         (nvec)          , & ! thermal flux at soil surface
    zthsnw         (nvec)          , & ! thermal flux at snow surface
    zfor_s         (nvec)          , & ! total forcing at soil surface
    zgsb           (nvec)          , & ! heat-flux through snow
    zrnet_s        (nvec)          , & ! net radiation
    zsprs          (nvec)          , & ! utility variable
    zahf           (nvec)          , & ! TERRA_URB: anthropogenic heat flux

    ! Tendencies
    zdwidt         (nvec)          , & ! tendency of water content of interception store
    zdwsndt        (nvec)          , & ! tendency of snow water content
    zdtsdt         (nvec)          , & ! tendency of zts
    zdtsnowdt      (nvec)          , & ! tendency of ztsnow
    zdwgdt         (nvec,  ke_soil)    ! tendency of water content [kg/(m**3 s)]

  REAL    (KIND=wp) ::  &

    !   Interception variables
    zwinstr        (nvec)          , & ! preliminary water content of interception store
    zinfmx         (nvec)          , & ! maximum infiltration rate
    zwimax         (nvec)          , & ! maximum water content of interception store
    zvers          (nvec)          , & ! water supply for infiltration
    zwisnstr       (nvec)          , & ! water content of snow interception store (t+1) (mH2O)
    zwpnstr        (nvec)          , & ! water content of pond store (t+1) (mH2O)
    zfd_wi         (nvec)          , & ! surface fraction of intercepted water
    zf_pd          (nvec)          , & ! surface fraction covered by pond
    zwisn          (nvec)          , & ! water content of snow interception store (mH2O)
    zwpn           (nvec)          , & ! water content of pond
    zewi           (nvec)          , & ! evaporation from interception store
    zepd           (nvec)          , & ! evaporation from pond
    zept           (nvec)          , & ! potential evaporation
    zesn           (nvec)          , & ! evaporation from snow
    zdrr           (nvec)          , & ! formation of dew
    zrrs           (nvec)              ! formation of rime

  REAL    (KIND=wp) ::  &
    zg1            (nvec)          , & ! heat flux between two uppermost soil layers (W/m**2)
    zlhfl          (nvec)          , & ! estimated       latent    heat flux surface (W/m**2)
    zshfl          (nvec)          , & ! estimated       sensible  heat flux surface (W/m**2)
    zthfl          (nvec)          , & ! estimated total turbulent heat flux surface (W/m**2)
    zradfl         (nvec)          , & ! total radiative flux at surface             (W/m**2)
    ze_melt        (nvec)          , & ! energy required for melting of snow         (J/m**2)
    zch_snow       (nvec)          , & ! heat capacity of snow * w_snow * Rho_w      (J/m**2 K)
    zeb1           (nvec)          , & ! estimated energy budget of first soil layer (W/m**2)
    ztchv          (nvec)          , & ! transfer coefficient multiplied by wind velocity
    ztchv_max      (nvec)          , & ! dito          , but upper limit
    zrho_atm       (nvec)          , & ! air density of lowest atmospheric layer     (kg/m**3)
    zdt_atm        (nvec)          , & ! temperature difference between lowest
    zdq_atm        (nvec)          , & ! dito for moisture (kg/kg)

    !additional variables for root distribution
    zroota         (nvec)          , & ! root density profile parameter (1/m)
    zwrootdz       (nvec, ke_soil) , & ! mean water content over root depth weighted by root density
    zrootdz_int    (nvec)          , & ! parameter needed to initialize the root density profile integral
    zwrootdz_int   (nvec)          , & ! parameter needed to initialize the root water content integral
    zqhfl_s        (nvec)          , & ! moisture flux at soil/air interface
    zqhfl_snow     (nvec)              ! moisture flux at snow/air interface

  REAL    (KIND=wp) ::  &

    ! Soil and plant parameters
    zroc     (nvec,ke_soil+1)      , & ! heat capacity
    zfcap    (nvec,ke_soil+1)      , & ! field capacity of soil
    zadp     (nvec,ke_soil+1)      , & ! air dryness point
    zporv    (nvec,ke_soil+1)      , & ! pore volume (fraction of volume)
    zdlam    (nvec)                , & ! heat conductivity parameter
    zdw      (nvec,ke_soil+1)      , & ! hydrological diff.parameter
    zdw1     (nvec,ke_soil+1)      , & ! hydrological diff.parameter
    zkw0     (nvec)                , & ! hydrological cond.parameter at surface
    zkw      (nvec,ke_soil+1)      , & ! hydrological cond.parameter
    zkwm     (nvec,ke_soil+1)      , & ! hydrological cond.parameter (main levels)
    zkw1     (nvec,ke_soil+1)      , & ! hydrological cond.parameter
    zik2     (nvec)                , & ! minimum infiltration rate
    zpwp     (nvec,ke_soil+1)      , & ! plant wilting point (fraction of volume)
    ztlpmwp  (nvec)                , & ! turgor-loss-point minus plant wilting point
    zedb     (nvec)                , & ! utility variable
    zaa      (nvec)                , & ! utility variable

    ! Hydraulic variables
    ztrang      (nvec,ke_soil)     , & ! transpiration contribution by the different soil layers
    ztrangs     (nvec)             , & ! total transpiration (transpiration from all soil layers)
    zwin        (nvec)             , & ! water content of interception store   (m H20)
    zwsnow      (nvec)             , & ! snow water equivalent                 (m H20)
    zwsnew      (nvec)             , & ! snow water equivalent                 (m H20)
    zdwsnm      (nvec)             , & ! utility variable for snow melt determination
    zw_fr       (nvec,ke_soil+1)   , & ! fractional total water content of the different soil layers
    zinfil      (nvec)             , & ! infiltration rate
    zlw_fr      (nvec,ke_soil+1)   , & ! fractional liquid water content of soil layers
    ziw_fr      (nvec,ke_soil+1)   , & ! fractional ice content of soil layers
    zwsnn       (nvec)             , & ! new value of zwsnow
    zflmg       (nvec,ke_soil+1)   , & ! flux of water at soil layer interfaces
    zrunoff_grav(nvec,ke_soil+1)       ! main level water gravitation
                                       !    flux (for runoff calc.)

  REAL    (KIND=wp) ::  &

    ! Local arrays for BATS-scheme
    zk0di       (nvec)             , & ! surface type dependent parameter
    zbedi       (nvec)             , & ! surface type dependent parameter
    zsnull      (nvec)             , & ! mean relative(zporv) water fraction of the
                                       !    so-called total active soil layer (above 1m)
    zs1         (nvec)             , & ! mean relative (zporv) water fraction  of
                                       !    layers above 0.1m
    zf_rad      (nvec)             , & ! radiation function for stomata resistance
    ztraleav    (nvec)             , & ! transpiration rate of dry leaves

    ! Root functions
    zwroot      (nvec)             , & ! mean water content over root depth
    zropartw    (nvec,ke_soil)         ! fraction of layer filled by roots * w_fr

  REAL    (KIND=wp) ::  &

    ! Thermal variables
    zts         (nvec)             , & ! soil surface temperature
    ztsk        (nvec)             , & ! skin temperature
    ztsnow      (nvec)             , & ! snow surface temperaure
    ztsnow_mult (nvec,0:ke_snow)   , & ! snow surface temperaure

    ! HEATCOND (soil moisture dependent thermal conductivity of the soil)
    zalamtmp    (nvec,ke_soil)     , & ! thermal conductivity
    zalam       (nvec,ke_soil)     , & ! thermal conductivity on half levels, meaning at the (lower) boundaries of soil layers 1 to ke_soil
    zrocg       (nvec,ke_soil+1)   , & ! total volumetric heat capacity of soil
    zrocg_soil  (nvec,ke_soil+1)   , & ! volumetric heat capacity of bare soil
    zrocs       (nvec)             , & ! heat capacity of snow
    ztsn        (nvec)             , & ! new value of zts
    ztsnown     (nvec)             , & ! new value of ztsnow
    ztsnown_mult(nvec,0:ke_snow)   , & ! new value of ztsnow
    znlgw1f     (ke_soil)              ! utility variable

  REAL(KIND=wp)                  ::  &

    !   Multi-snow layer parameters
    zqbase       (nvec)            , &
    zrefr        (nvec)            , & ! rate of liquid water refreezing
    zmelt        (nvec)            , & ! rate of snow melting
    ze_out       (nvec)            , &
    zrho_dry_old (nvec)            , &
    zp           (nvec,ke_snow)    , &
    zcounter     (nvec)            , &
    ze_rad       (nvec)            , &
    zswitch      (nvec)            , &
    tmp_num      (nvec)            , &
    sum_weight   (nvec)            , &
    t_new        (nvec,ke_snow)    , &
    rho_new      (nvec,ke_snow)    , &
    wl_new       (nvec,ke_snow)    , &
    z_old        (nvec,ke_snow)    , &
    dz_old       (nvec,ke_snow)    , &
    t_so_free_new(nvec,0:ke_soil+1), &
    t_so_snow_new(nvec,0:ke_soil+1), &
    sn_frac      (nvec)            , &
    zf_snow_lim  (nvec)            , &
    zdz_snow     (nvec)            , & ! snow depth (not for snow heat content but snow
                                       !             density calculation)
    zalas_mult    (nvec,  ke_snow) , & ! heat conductivity of snow
    ztsnownew_mult(nvec,0:ke_snow) , & ! preliminary value of snow surface temperature
    zextinct      (nvec,  ke_snow) , & ! solar radiation extinction coefficient in snow (1/m)
    zfor_snow_mult(nvec)               ! total forcing at snow surface

  REAL    (KIND=wp) ::  &
    ! Auxiliary variables
    hzalam      (nvec,ke_soil+1)   , & ! thermal conductivity on full levels, meaning at the centers of soil layers 1 to ke_soil+1
    zdqvtsnow   (nvec)             , & ! first derivative of saturation specific humidity
                                       !    with respect to t_snow
    zrho_snow   (nvec)             , & ! snow density used for computing heat capacity and conductivity
    zts_pm      (nvec)             , & ! indicator zts > < T_melt
    ztsk_pm     (nvec)             , & ! indicator ztsk > < T_melt
    ztfunc      (nvec)             , & ! smoothed transition function between T_melt and T_melt+2K
    ztsnow_pm   (nvec)             , & ! indicator ztsnow > < T_melt
    zeisa       (nvec)                 ! evaporation from impervious surfaces (puddles)

  REAL    (KIND=wp) ::  &
    ! Utility variables
    zaga        (nvec,0:ke_soil+ke_snow+1), &
    zagb        (nvec,0:ke_soil+ke_snow+1), &
    zagc        (nvec,0:ke_soil+ke_snow+1), &
    zagd        (nvec,0:ke_soil+ke_snow+1), &
    zage        (nvec,0:ke_soil+ke_snow+1)


  REAL    (KIND=wp) ::  &
    zhwso_budget(nvec)      ! utility variables for soil water budget

  LOGICAL     :: &
    limit_tch (nvec)         ! indicator for flux limitation problem


#ifdef __SX__
  REAL(wp) :: zfac(nvec),lhfl_pl_int(nvec)
#endif

  ! ground water as lower boundary of soil column
  REAL    (KIND=wp) ::  &
    zdelta_sm, zdhydcond_dlwfr

  ! HEATCOND
  REAL    (KIND=wp) ::          &
    zthetas, zlamli, zlamsat, zlams, rsandf, zlamq, zlam0, zrhod, zlamdry,  &
    zlamdry_soil, zsri, zKe, zthliq, zlamic, zlamdry_c1, zlamdry_c2,        &
    zlamspeat, zlamdrypeat, zkeai, zkebi, zkeal, zkebl !AYu 

  ! For performance improvement
  REAL    (KIND=wp) :: ln_2, ln_3, ln_10, ln_006

  INTEGER :: my_cart_id,        &
             itype_mire   = 0

  INTEGER :: my_thrd_id, mcid, mtid, mbid, mvid

  LOGICAL ::  ldebug = .FALSE.

!- End of header
!==============================================================================

CALL set_acc_host_or_device(lzacc, lacc)

IF(PRESENT(opt_acc_async_queue)) THEN
  acc_async_queue = opt_acc_async_queue
ELSE
  acc_async_queue = 1
ENDIF

my_cart_id = get_my_global_mpi_id()
#ifdef _OPENMP
my_thrd_id = omp_get_thread_num()
#endif
mbid = 714
mcid =   0
mtid =   0
mvid =   8

!------------------------------------------------------------------------------
! Begin Subroutine terra
!------------------------------------------------------------------------------
!==============================================================================
!  Computation of the diagnostic part I of the soil parameterization scheme
!  In this part, evaporation from the surface and transpiration is calculated.
!  A multi-layer soil water distribution is calculated by simple bulk
!  parameterisation or a Penman-Monteith-type version of evaporation
!  and transpiration which can be used alternatively.
!------------------------------------------------------------------------------

! Just do some checkout prints:
  IF (ldebug) THEN
    IF (iblock == mbid .AND. my_cart_id == mcid) THEN
#ifdef _OPENMP
     IF (my_thrd_id == mtid) THEN
#endif
      WRITE(*,'(A,3I5)'   ) 'SFC-DIAGNOSIS terra start (ke_soil, ke_snow, ke_soil_hy): ', &
                            ke_soil, ke_snow, ke_soil_hy
#ifdef _OPENMP
     ENDIF
#endif
    ENDIF
    DO i = ivstart, ivend
      IF (i== mvid .AND. iblock == mbid .AND. my_cart_id == mcid) THEN
#ifdef _OPENMP
       IF (my_thrd_id == mtid) THEN
#endif
        WRITE(*,'(A,2I5   )') ' SFC-DIAGNOSIS terra:  iblock = ', iblock, i
 
        WRITE(*,'(A       )') ' External Parameters:  '
        WRITE(*,'(A,I28   )') '   soiltyp          :  ', soiltyp_subs(i)
! for TERRA_URB
        WRITE(*,'(A,F28.16)') '   urb_isa          :  ', urb_isa     (i)
        WRITE(*,'(A,F28.16)') '   urb_ai           :  ', urb_ai      (i)
        WRITE(*,'(A,F28.16)') '   urb_h_bld        :  ', urb_h_bld   (i)
        WRITE(*,'(A,F28.16)') '   urb_hcap         :  ', urb_hcap    (i)
        WRITE(*,'(A,F28.16)') '   urb_hcon         :  ', urb_hcon    (i)
        WRITE(*,'(A,F28.16)') '   ahf              :  ', ahf         (i)
!
        WRITE(*,'(A,F28.16)') '   plcov            :  ', plcov       (i)
        WRITE(*,'(A,F28.16)') '   rootdp           :  ', rootdp      (i)
        WRITE(*,'(A,F28.16)') '   sai              :  ', sai         (i)
        WRITE(*,'(A,F28.16)') '   tai              :  ', tai         (i)
        WRITE(*,'(A,F28.16)') '   eai              :  ', eai         (i)
        WRITE(*,'(A,F28.16)') '   skinc            :  ', skinc       (i)
        WRITE(*,'(A,F28.16)') '   rsmin2d          :  ', rsmin2d     (i)
        WRITE(*,'(A,F28.16)') '   r_bsmin          :  ', r_bsmin     (i)
        WRITE(*,'(A       )') ' Other input parameters:'
        WRITE(*,'(A,F28.16)') '   u     ke         :  ', u           (i)
        WRITE(*,'(A,F28.16)') '   v     ke         :  ', v           (i)
        WRITE(*,'(A,F28.16)') '   t     ke         :  ', t           (i)
        WRITE(*,'(A,F28.16)') '   qv    ke         :  ', qv          (i)
      IF (ldepo_qw) THEN
        WRITE(*,'(A,F28.16)') '   qc    ke         :  ', qc          (i)
       IF (ldiff_qi) &
        WRITE(*,'(A,F28.16)') '   qi    ke         :  ', qi          (i)
      END IF
        WRITE(*,'(A,F28.16)') '   ptot  ke         :  ', ptot        (i)
        WRITE(*,'(A,F28.16)') '   ps               :  ', ps          (i)
        WRITE(*,'(A,F28.16)') '   h_snow_gp        :  ', h_snow_gp   (i)
        WRITE(*,'(A,F28.16)') '   u_10m            :  ', u_10m       (i)
        WRITE(*,'(A,F28.16)') '   v_10m            :  ', v_10m       (i)
        WRITE(*,'(A,F28.16)') '   prr_con          :  ', prr_con     (i)
        WRITE(*,'(A,F28.16)') '   prs_con          :  ', prs_con     (i)
        WRITE(*,'(A,F28.16)') '   conv_frac        :  ', conv_frac   (i)
        WRITE(*,'(A,F28.16)') '   prr_gsp          :  ', prr_gsp     (i)
        WRITE(*,'(A,F28.16)') '   prs_gsp          :  ', prs_gsp     (i)
        WRITE(*,'(A,F28.16)') '   prg_gsp          :  ', prg_gsp     (i)
#ifdef TWOMOM_SB
        WRITE(*,'(A,F28.16)') '   prh_gsp          :  ', prh_gsp     (i)
#endif
        WRITE(*,'(A,F28.16)') '   sobs             :  ', sobs        (i)
        WRITE(*,'(A,F28.16)') '   thbs             :  ', thbs        (i)
        WRITE(*,'(A,F28.16)') '   pabs             :  ', pabs        (i)
        WRITE(*,'(A       )') ' Soil initial parameters:'
        WRITE(*,'(A,F28.16)') '   t_snow_now       :  ', t_snow_now  (i)
        WRITE(*,'(A,F28.16)') '   t_s_now          :  ', t_s_now     (i)
        WRITE(*,'(A,F28.16)') '   t_g              :  ', t_g         (i)
DO k = 0, ke_soil+1
        WRITE(*,'(A,I1,A,F28.16)') '   t_so_now    (',k,')      :  ', t_so_now    (i,k)
ENDDO
DO k = 1, ke_soil+1
        WRITE(*,'(A,I1,A,F28.16)') '   w_so_now    (',k,')      :  ', w_so_now    (i,k)
ENDDO
DO k = 1, ke_soil+1
        WRITE(*,'(A,I1,A,F28.16)') '   w_so_ice_now   (',k,')      :  ',  w_so_ice_now (i,k)
ENDDO
        WRITE(*,'(A,F28.16)') '   qv_s (in)        :  ', qv_s        (i)
        WRITE(*,'(A,F28.16)') '   w_snow_now       :  ', w_snow_now  (i)
        WRITE(*,'(A,F28.16)') '   rho_snow_now     :  ', rho_snow_now(i)
        WRITE(*,'(A,F28.16)') '   h_snow (in)      :  ', h_snow      (i)
        WRITE(*,'(A,F28.16)') '   w_i_now          :  ', w_i_now     (i)
        WRITE(*,'(A,F28.16)') '   w_p_now          :  ', w_p_now     (i)
        WRITE(*,'(A,F28.16)') '   w_s_now          :  ', w_s_now     (i)
        WRITE(*,'(A,F28.16)') '   freshsnow (in)   :  ', freshsnow   (i)
        WRITE(*,'(A,F28.16)') '   zf_snow (in)     :  ', zf_snow     (i)
        WRITE(*,'(A,F28.16)') '   tch              :  ', tch         (i)
        WRITE(*,'(A,F28.16)') '   tcm              :  ', tcm         (i)
        WRITE(*,'(A,F28.16)') '   tfv              :  ', tfv         (i)
        WRITE(*,'(A,F28.16)') '   runoff_s (in)    :  ', runoff_s    (i)
        WRITE(*,'(A,F28.16)') '   runoff_g (in)    :  ', runoff_g    (i)

#ifdef _OPENMP
       ENDIF
#endif
      ENDIF
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section I.1: Initializations
!------------------------------------------------------------------------------

  ierror = 0
  yerror = '        '

  IF ( (itype_trvg == 3) .AND. ( (.NOT. PRESENT(z0)) .OR. (.NOT. PRESENT (plevap)) ) ) THEN
    IF (my_cart_id == 0) THEN
      PRINT *, ' *** ERROR:  Arguments (z0/plevap) not present in TERRA for itype_trvg=3! ***'
    ENDIF
    ierror = 1   ! that must not happen
    RETURN       ! whatever happens then
  ENDIF

#ifdef _OPENACC
  IF(lmulti_snow) CALL finish("sfc_terra:terra", "lmulti_snow has not been tested with OpenACC yet.")
#endif

  ! set number of active soil moisture layers for Mire
  ke_soil_hy_m            = 4

  ! time step for soil variables
  zdt      = dt

  ! time step for run-off computation
  zroffdt  = zdt

! Computation of derived constants
  zrhwddt  = rho_w/zdt       ! density of liquid water/timestep
  zdtdrhw  = zdt/rho_w       ! timestep/density of liquid water

  zdwi_scal = zdt/1800.0_wp  ! time scale parameter for freezing/melting soil water

! time constant for infiltration of water from interception store
  ctau_i   = MAX(ctau_i,zdt)

  ! Subroutine parameters IN
  !$ACC DATA &
  !$ACC   PRESENT(ivend) &
  !$ACC   PRESENT(zmls) &
  !$ACC   PRESENT(soiltyp_subs, urb_isa, urb_ai, urb_h_bld) &
  !$ACC   PRESENT(urb_hcap, urb_hcon, ahf) &
  !$ACC   PRESENT(plcov, rootdp, sai, eai, tai, laifac, skinc) &
  !$ACC   PRESENT(heatcond_fac, heatcap_fac, hydiffu_fac) &
  !$ACC   PRESENT(rsmin2d, r_bsmin, u, v, t, qv, qc, qi, ptot, ps, h_snow_gp, u_10m) &
  !$ACC   PRESENT(v_10m, prr_con, prs_con, conv_frac, prr_gsp, prs_gsp, pri_gsp) &
#ifdef TWOMOM_SB
  !$ACC   PRESENT(prh_gsp) &
#endif
  !$ACC   PRESENT(prg_gsp, sobs, thbs, pabs, zdzhs, tsnred) &

  ! Subroutine parameters INOUT
  !$ACC   PRESENT(t_snow_now, t_s_now, t_sk_now, t_g) &
  !$ACC   PRESENT(qv_s, w_snow_now) &
  !$ACC   PRESENT(rho_snow_now, h_snow, w_i_now, w_p_now, w_s_now) &
  !$ACC   PRESENT(freshsnow, zf_snow, tch, tcm, tfv, tfvsn, runoff_s) &
  !$ACC   PRESENT(runoff_g, t_snow_mult_now, rho_snow_mult_now) &
  !$ACC   PRESENT(wliq_snow_now, wtot_snow_now, dzh_snow_now) &
  !$ACC   PRESENT(t_so_now, w_so_now, w_so_ice_now) &

  ! Subroutine parameters OUT
  !$ACC   PRESENT(t_snow_new, t_s_new, t_sk_new) &
  !$ACC   PRESENT(w_snow_new, rho_snow_new) &
  !$ACC   PRESENT(meltrate, w_i_new, w_p_new, w_s_new, zshfl_s) &
  !$ACC   PRESENT(zlhfl_s, zshfl_snow, zlhfl_snow, rstom) &
  !$ACC   PRESENT(lhfl_bs, t_snow_mult_new, rho_snow_mult_new) &
  !$ACC   PRESENT(wliq_snow_new, wtot_snow_new, dzh_snow_new) &
  !$ACC   PRESENT(t_so_new, w_so_new, w_so_ice_new, lhfl_pl) &
  !$ACC   PRESENT(zshfl_sfc, zlhfl_sfc, zqhfl_sfc) &

  ! Local arrays
  !$ACC   CREATE(m_styp, ke_soil_hy_b) &
  !$ACC   CREATE(h_snow_now, h_snow_new, zzhls, zdzhs, zdzms) &
  !$ACC   CREATE(zdz_snow_fl, zhh_snow, zhm_snow, zdzh_snow) &
  !$ACC   CREATE(zdzm_snow, zdtsnowdt_mult, zbwt, zrock, zsandf) &
  !$ACC   CREATE(zclayf, zsiltf, zb_por, zpsis, zw_m) &
  !$ACC   CREATE(zrr, zrs, zesoil, zsfc_frac_bs) &
  !$ACC   CREATE(zw_m_org, zw_m_soil, zw_m_up, zw_m_low, zaa) &
  !$ACC   CREATE(zrhoch, zth_low, zf_wi, ztmch, zep_s, zep_snow) &
  !$ACC   CREATE(zverbo, zversn, zthsoi, zthsnw, zfor_s, zgsb) &
  !$ACC   CREATE(zrnet_s, zsprs, zahf, zdwidt, zdwsndt, zdtsdt) &
  !$ACC   CREATE(zdtsnowdt, zdwgdt, zwinstr, zinfmx, zwimax, zvers) &
  !$ACC   CREATE(zwisnstr, zwpnstr, zfd_wi, zf_pd, zwisn, zwpn) &
  !$ACC   CREATE(zewi, zepd, zept, zesn, zdrr, zrrs, zg1, zlhfl) &
  !$ACC   CREATE(zshfl, zthfl, zradfl, ze_melt, zch_snow, zeb1) &
  !$ACC   CREATE(ztchv, ztchv_max, zrho_atm, zdt_atm, zdq_atm) &
  !$ACC   CREATE(zroota, zwrootdz, zwrootdz_int, zrootdz_int) &
  !$ACC   CREATE(zqhfl_s, zqhfl_snow, zroc, zfcap, zadp, zporv) &
  !$ACC   CREATE(zdlam, zdw, zdw1, zkw0, zkw, zkwm, zkw1, zik2, zpwp, ztlpmwp) &
  !$ACC   CREATE(zedb, ztrang, ztrangs, zwin, zwsnow, zwsnew, zdwsnm) &
  !$ACC   CREATE(zw_fr, zinfil, zlw_fr, ziw_fr, zwsnn, zflmg) &
  !$ACC   CREATE(zrunoff_grav, zk0di, zbedi, zsnull, zs1, zf_rad) &
  !$ACC   CREATE(ztraleav, zwroot, zropartw, zts, ztsk, ztsnow) &
  !$ACC   CREATE(ztsnow_mult, zalamtmp, zalam, zrocg, zrocg_soil) &
  !$ACC   CREATE(zrocs, ztsn, ztsnown, ztsnown_mult, znlgw1f) &
  !$ACC   CREATE(zqbase, zrefr, zmelt, ze_out, zrho_dry_old) &
  !$ACC   CREATE(zp, zcounter, ze_rad, zswitch, tmp_num, sum_weight) &
  !$ACC   CREATE(t_new, rho_new, wl_new, dz_old, z_old) &
  !$ACC   CREATE(t_so_free_new, t_so_snow_new, sn_frac) &
  !$ACC   CREATE(zf_snow_lim, zdz_snow, zalas_mult, ztsnownew_mult) &
  !$ACC   CREATE(zextinct, zfor_snow_mult, hzalam, zdqvtsnow) &
  !$ACC   CREATE(zrho_snow, zts_pm, ztsk_pm, ztfunc, ztsnow_pm, zeisa) &
  !$ACC   CREATE(zagd, zage, limit_tch, zdz_snow) &
  !$ACC   CREATE(zaga, zagb, zagc) &

  ! Terra data module fields
  !$ACC   PRESENT(cporv, cfcap, cpwp, cadp, cik2, ckw0, ckw1, cdw0) &
  !$ACC   PRESENT(cdw1, crock, crhoc, cala1, cala0, ck0di, cbedi) &
  !$ACC   PRESENT(clgk0, csandf, cclayf) ASYNC(acc_async_queue)

  ! for the optional fields for itype_trvg
  !$ACC DATA PRESENT(plevap, z0) ASYNC(acc_async_queue) IF(PRESENT(plevap))

  ! for the soil water budget
  !$ACC DATA PRESENT(resid_wso) CREATE(zhwso_budget) ASYNC(acc_async_queue) IF(lres_soilwatb)

  IF (lres_soilwatb) THEN
    ! Calculation of soil water budget: store recent runoff fluxes
    ! the final runoff fluxes are subtracted at the end
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue)
    !$ACC LOOP GANG VECTOR
    DO i = ivstart, ivend
      zhwso_budget(i) = runoff_s(i) + runoff_g(i)
      ! in COSMO another procedure is used than in ICON to 
      !   initialize the values of runoff_s and runoff_g (COSMO:
      !   actual values, ICON: zeros) -> for ICON the following init is possible
      !   zhwso_budget(i) = 0.0_wp
    ENDDO
    !$ACC END PARALLEL
  ENDIF

  ln_10 = LOG(10.0_wp)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP SEQ
  DO kso = 0, UBOUND(zaga,2)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      zaga(i,kso) = 1.0_wp
      zagb(i,kso) = 1.0_wp
      zagc(i,kso) = 1.0_wp
    ENDDO
  ENDDO

! Temperaturedifference for liquid water content in frozen soil at -40 degC
!  J. Helmert: Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
!                                                    doi:10.5194/bg-13-1991-2016
  t_zw_up  = 270.15_wp     ! temp -3 degC
  t_zw_low = 233.15_wp     ! temp -40 degC

! Prepare basic surface properties (for land-points only)

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp)
  DO i = ivstart, ivend
    mstyp     = soiltyp_subs(i)        ! soil type
    m_styp(i) = mstyp                  ! array for soil type


    ! ensure that glaciers are covered with at least 1 m of snow
    IF (mstyp == 1) h_snow(i) = MAX(1._wp, h_snow(i))

#ifndef __SX__
    zdw       (i,:) = cdw0  (mstyp)
    zdw1      (i,:) = cdw1  (mstyp)
    zkw       (i,:) = ckw0  (mstyp)
    zkwm      (i,:) = ckw0  (mstyp)
    zkw1      (i,:) = ckw1  (mstyp)
    zporv     (i,:) = cporv (mstyp)              ! pore volume
    zpwp      (i,:) = cpwp  (mstyp)              ! plant wilting point
    zadp      (i,:) = cadp  (mstyp)              ! air dryness point
    zfcap     (i,:) = cfcap (mstyp)              ! field capacity
    zrocg     (i,:) = crhoc (mstyp)*heatcap_fac(i)   ! heat capacity
    zrocg_soil(i,:) = crhoc (mstyp)*heatcap_fac(i)   ! heat capacity
    zalam     (i,:) = cala0 (mstyp)              ! heat conductivity parameter
#endif
    zrock     (i)   = crock (mstyp)              ! EQ 0 for Ice and Rock EQ 1 else
    zik2      (i)   = cik2  (mstyp)
    zdlam     (i)   = cala1 (mstyp)-cala0(mstyp) ! heat conductivity parameter
    zbwt      (i)   = MAX(0.001_wp,rootdp(i))    ! Artificial minimum value for root depth
    zroota    (i)   = 3.0_wp/zbwt(i)             ! root density profile parameter (1/m)
                                                 ! zroota=0. creates the original TERRA_LM
                                                 ! version with constant root density

    ! New arrays for BATS-scheme
    zk0di     (i)   = ck0di(mstyp)               !
    zbedi     (i)   = cbedi(mstyp)               !

    meltrate(i)     = 0.0_wp
  ENDDO

#ifdef __SX__
!$NEC outerloop_unroll(8)
  DO kso = 1, ke_soil+1
    DO i = ivstart, ivend
    mstyp       = soiltyp_subs(i)        ! soil type
    zdw   (i,kso)  = cdw0  (mstyp)
    zdw1  (i,kso)  = cdw1  (mstyp)
    zkw   (i,kso)  = ckw0  (mstyp)
    zkwm  (i,kso)  = ckw0  (mstyp)
    zkw1  (i,kso)  = ckw1  (mstyp)
    zporv(i,kso)  = cporv(mstyp)              ! pore volume
    zpwp (i,kso)  = cpwp (mstyp)              ! plant wilting point
    zadp (i,kso)  = cadp (mstyp)              ! air dryness point
    zfcap(i,kso)  = cfcap(mstyp)              ! field capacity
    zrocg(i,kso)  = crhoc(mstyp)*heatcap_fac(i)      ! heat capacity
    zrocg_soil(i,kso) = crhoc(mstyp)*heatcap_fac(i)  ! heat capacity
   ENDDO
  ENDDO
!$NEC outerloop_unroll(7)
  DO kso = 1, ke_soil
    DO i = ivstart, ivend
    mstyp       = soiltyp_subs(i)        ! soil type
    zalam(i,kso)  = cala0(mstyp)              ! heat conductivity parameter
   ENDDO
  ENDDO
#endif

  ! Arrays for soil water freezing/melting
  zd = LOG((T_ref_ice-(t_zw_low-t0_melt))/T_star_ice)
  zd1 = EXP(b_sand*zd)
  zd2 = EXP(b_clay*zd)
  zd3 = EXP(b_silt*zd)
  zd4 = EXP(b_org*zd)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp)
  DO i = ivstart, ivend
    mstyp       = soiltyp_subs(i)        ! soil type
    zsandf(i)   = csandf(mstyp)
    zclayf(i)   = cclayf(mstyp)
    zsiltf(i)   = 100.0_wp -csandf(mstyp)-cclayf(mstyp) ! Residuum of sand and clay
    zpsis (i)   = -zpsi0 * EXP(ln_10*(1.88_wp-0.013_wp*zsandf(i)))
    zb_por(i)   = 2.91_wp + 0.159_wp*zclayf(i)
    zedb  (i)   = 1.0_wp/zb_por(i)
    zaa   (i)   = g*zpsis(i)/lh_f

    ! Liq. water content at -3 degC
    zw_m_up(i) = EXP(-zedb(i)*LOG((t_zw_up - t0_melt)/(t_zw_up*zaa(i))) ) ! Without zporv(i,kso)*zdzhs(kso)!!

    ! Determine liq. water content at -40 degC
    !  J. Helmert: Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
    !                                                    doi:10.5194/bg-13-1991-2016
    zw_m_soil(i) = 0.01_wp*(zsandf(i)*zd1 + zclayf(i)*zd2 + zsiltf(i)*zd3)
    zw_m_org(i) = zd4

    ! adaptive parameter tuning for near-surface hydraulic diffusivity
    zdw(i,1) = zdw(i,1)*hydiffu_fac(i)**2
    zdw(i,2) = zdw(i,2)*hydiffu_fac(i)
  ENDDO

  ! zkw0 will not be needed anymore if 'itype_interception = 2' is removed
  IF (itype_interception == 2) THEN
    DO i = ivstart, ivend
      !fc=2 1/m Exponential Ksat-profile decay parameter,see Decharme et al. (2006)
      zkw0   (i) = zkw   (i,1)*EXP(2.0_wp*rootdp(i))
    END DO
  ENDIF

  ! Set three-dimensional variables
!$NEC outerloop_unroll(7)
  !$ACC LOOP SEQ
  DO kso = 1, ke_soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zzz)
    DO i = ivstart, ivend
      !fc=2 1/m Exponential Ksat-profile decay parameter,see Decharme et al. (2006)
      zkw   (i,kso) = zkw   (i,kso)*EXP(-2._wp*(zzhls(kso)-rootdp(i)))
      zkwm  (i,kso) = zkwm  (i,kso)*EXP(-2._wp*(zmls(kso)-rootdp(i)))

      ! Scale soil heat capacity with organic fraction -> Chadburn et al., 2015
      IF (zmls(kso) < rootdp(i)) THEN
        zzz = plcov(i)*(rootdp(i)-zmls(kso))/rootdp(i)
        zrocg(i,kso)=(1.0_wp-zzz)*zrocg_soil(i,kso)+zzz*0.58E+06_wp
        !  J. Helmert: Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
        !  Organic fraction                                        doi:10.5194/bg-13-1991-2016
        zw_m_low(i,kso) = zporv(i,kso)*zdzhs(kso)*(zzz*zw_m_org(i) + (1.0_wp-zzz)*zw_m_soil(i))
      ELSE
        zw_m_low(i,kso) = zporv(i,kso)*zdzhs(kso)*zw_m_soil(i)
      END IF
    ENDDO
  ENDDO
  !$ACC END PARALLEL


  ! Determine constants clgk0 for BATS-scheme
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO jb       = 1, 10
    clgk0(jb) = LOG(MAX(eps_soil,ck0di(jb)/ckrdi))/ln_10
  END DO
  !$ACC END PARALLEL


  ! Determine the layer indices and total depth for water content averaging
  ! (in soil evaporation after Dickinson)
  k10cm  = 1
  k100cm = 1
  z1     = zzhls(1)
  znull  = zzhls(1)
  DO kso = 1, ke_soil
    IF (zmls(kso).le.0.1_wp) THEN
      z1      = zzhls(kso)
      k10cm   = kso
    END IF
    IF (zmls(kso).le.1.0_wp) THEN
      znull   = zzhls(kso)
      k100cm  = kso
    END IF
  END DO

  zpercmax = 2.E-3_wp

  ! ACCWA (nvhpc 23.3): NUM_GANGS(1024): The compiler was using a small number of
  ! gangs for this kernel for some reason, but by default if the loop bounds are
  ! present on the GPU only, the compiler would always use 1024 as some reasonable
  ! number
  !$ACC PARALLEL DEFAULT(PRESENT) NUM_GANGS(1024) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    ztrangs(i)            = 0.0_wp
    zsnull (i)            = 0.0_wp
    zs1    (i)            = 0.0_wp
    zwroot (i)            = 0.0_wp
    zdwidt (i)            = 0.0_wp         ! Initialisation of all
    zdwsndt(i)            = 0.0_wp         ! evaporation quantities
    zesoil (i)            = 0.0_wp         ! to be equal to zero
    zrr    (i)            = 0.0_wp         ! in first part formation of dew
    zrs    (i)            = 0.0_wp         ! in first part formation of rime
    zw_fr  (i,ke_soil+1)  = w_so_now(i,ke_soil+1)/zdzhs(ke_soil+1)
    lhfl_bs(i)            = 0.0_wp
    lhfl_pl(i,ke_soil+1)  = 0.0_wp
    rstom  (i)            = 0.0_wp
    zahf   (i)            = 0.0_wp         ! TERRA_URB: Anthropogenic heat flux
  ENDDO

  ! REORDER
!$NEC outerloop_unroll(7)
  !$ACC LOOP SEQ
  DO kso   = 1, ke_soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      zw_fr   (i,kso)     = w_so_now(i,kso)/zdzhs(kso)
      ztrang  (i,kso)     = 0.0_wp
      zropartw(i,kso)     = 0.0_wp
      lhfl_pl( i,kso)     = 0.0_wp
    ENDDO
  ENDDO


  ! Determine average soil water content over 10 and 100 cm, respectively.
  !$ACC LOOP SEQ
  DO kso   = 1, k100cm
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (kso.le.k10cm) zs1(i) = zs1(i) + w_so_now(i,kso)
      zsnull(i)   = zsnull(i) + w_so_now(i,kso)
    ENDDO
  ENDDO

  IF ( itype_mire == 1 ) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (soiltyp_subs(i) == 8) THEN
        ke_soil_hy_b(i) = 4 ! Set hydr. active layers for Mires 
      ELSE
        ke_soil_hy_b(i) = ke_soil_hy
      ENDIF
    END DO
    !WRITE(*,*) 'ITYPE_MIRE = ',itype_mire,' MINVAL(ke_soil_hy_b(nvec))= ', MINVAL(ke_soil_hy_b(:))
  END IF


  ! No soil moisture for Ice and Rock
!$NEC outerloop_unroll(8)
  !$ACC LOOP SEQ
  DO kso   = 1, ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (soiltyp_subs(i) < 3) THEN
        w_so_now(i,kso)         = 0._wp
        w_so_ice_now(i,kso)     = 0._wp
      END IF
    END DO
  END DO

!$NEC outerloop_unroll(8)
  !$ACC LOOP SEQ
  DO kso   = 1, ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      w_so_new(i,kso)         = w_so_now(i,kso)
      w_so_ice_new(i,kso)     = w_so_ice_now(i,kso)
    END DO
  END DO

  ! Decide which snow density is used for computing the heat capacity
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (l2lay_rho_snow) THEN
      zrho_snow(i) = rho_snow_mult_now(i,1)
    ELSE
      zrho_snow(i) = rho_snow_now(i)
    ENDIF
  END DO


  IF (itype_heatcond >= 2) THEN

    ! heat conductivity dependent on actual soil water content
    ! based on Peters-Lidard et al. (1998) and Johansen (1975),
    ! see also Block, Alexander (2007), Dissertation BTU Cottbus
    zlamli    = LOG(0.57_wp)       ! LOG(thermal conductivity of water)
    zlamic    = LOG(2.2_wp)        ! LOG(thermal conductivity of ice)
    zlamq     = LOG(7.7_wp)        ! LOG(thermal conductivity of quartz)
    ln_2      = LOG(2.0_wp)
    ln_3      = LOG(3.0_wp)
    ln_006    = LOG(0.06_wp)
    zlamspeat = LOG(0.25_wp)! AYu lambdas for peat (Lawrence & Slater, 2008)

    zlamdrypeat =   0.05_wp ! AYu lambdas for peat (Lawrence & Slater, 2008) 
    zkeai       = 0.6116_wp ! AYu 4 below for the new Kersten number parametrization 
                            ! for the peat (Russian construction database), ice and liquid water
    zkebi       = 1.4123_wp
    ! As long as there is soil type "peat" there is no need to make a correction 
    ! for the dependency of thermal conductivity on the water, therefore zkeal and zkebl are not needed 
    !    zkeal = 1.5858_wp     ! 
    !    zkebl = 1.6879_wp     ! 

    ! tuning constants for dry thermal conductivity formula
    IF (itype_evsl == 4 .OR. itype_evsl == 5) THEN
      zlamdry_c1 = 437.0_wp
      zlamdry_c2 = 0.901_wp
    ELSE
      zlamdry_c1 = 64.7_wp
      zlamdry_c2 = 0.947_wp
    ENDIF


!$NEC outerloop_unroll(8)
    !$ACC LOOP SEQ
    DO kso = 1, ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zthetas, zthliq, rsandf, zlams) &
      !$ACC   PRIVATE(zlam0, zlamsat, zrhod, zlamdry_soil, zzz, zlamdry, zsri) &
      !$ACC   PRIVATE(zKe, zxx)
      DO i = ivstart, ivend
        zthetas = zporv(i,kso)                                 ! porosity
        zthliq  = zthetas - w_so_ice_now(i,kso)/zdzhs(kso) ! unfrozen volume fraction

        rsandf = zsandf(i)/100._wp                     ! quartz content

        if (rsandf >= 0.2_wp)  zlam0 = ln_2      ! LOG(thermal conductivity non-quartz)
        if (rsandf <  0.2_wp)  zlam0 = ln_3

  ! saturated thermal conductivity

        IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) THEN
          zlams = zlamspeat ! AYu for peat
        ELSE
          zlams = zlamq*rsandf + zlam0*(1._wp-rsandf)  ! LOG(solids thermal conductivity)
        ENDIF

        zlamsat = EXP(zlams*(1.0_wp-zthetas) + zlamic*(zthetas-zthliq) + zthliq*zlamli)

  ! dry thermal conductivity

        zrhod   = 2700.0_wp*(1.0_wp-zthetas)       ! dry density
        zlamdry_soil = ( 0.135_wp*zrhod + zlamdry_c1 )                     &
                / ( 2700.0_wp - zlamdry_c2*zrhod )
        ! missing: crushed rock formulation for dry thermal conductivity (see PL98)
  ! Scale zlamdry with organic fraction
        IF(zmls(kso) < rootdp(i)) THEN
          zzz = plcov(i)*(rootdp(i)-zmls(kso))/rootdp(i)
          zlamdry = EXP(LOG(zlamdry_soil)*(1._wp-zzz)+ln_006*zzz) ! Chadburn et al.,2015, Dankers et al., 2011
        ELSE
          zzz = 0._wp
          zlamdry=zlamdry_soil
        END IF

        IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) zlamdry = zlamdrypeat !AYu for peat

  ! Kersten number

        zsri = MIN(1.0_wp, w_so_now(i,kso)/zdzhs(kso) / zthetas) ! degree of saturation

        IF ( t_so_now(i,kso) < t0_melt) THEN                         ! frozen
          zKe = zsri
          IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1)  zke = zkeai*zsri**zkebi !AYu for peat
        ELSE                                                         ! unfrozen
          zKe = 0.0_wp
          IF ( soiltyp_subs(i) == 3 .or. soiltyp_subs(i) == 4 ) THEN ! coarse soil
            IF ( zsri >= 0.05_wp ) THEN
              zKe = 0.7_wp*LOG(zsri)/ln_10 + 1.0_wp
            ENDIF
          ELSE                                                       ! fine soil (other)
            IF ( zsri >= 0.1_wp ) THEN
              zKe = LOG(zsri)/ln_10 + 1.0_wp
            ENDIF
          ENDIF
          IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) zKe = MAX(0.0_wp, zKe) ! AYu can be more than 1 for the peat
        ENDIF
        zKe = MAX(0.0_wp, (MIN(1.0_wp, zKe)) )

  ! thermal conductivity

        ! tuning factor to indirectly account for the impact of vegetation, which does not depend on soil moisture
        IF(itype_heatcond == 3 .AND. zmls(kso) < 0.075_wp) THEN
          zxx = 12.5_wp*(0.075_wp-zmls(kso))*zzz
        ELSE
          zxx = 0.0_wp
        ENDIF
        hzalam(i,kso) = heatcond_fac(i)*(zKe*(zlamsat - zlamdry) + zlamdry)*(1._wp-zxx) + zxx*0.06_wp

        ! heat conductivity is also artificially reduced on snow-free forest-covered tiles generated
        ! by the melting-rate parameterization
        IF (tsnred(i) < -1.0_wp .AND. z0(i) >= 0.4_wp) THEN
          zxx = MAX(0.0_wp,2.0_wp - ABS(tsnred(i)))
          hzalam(i,kso) = zxx*hzalam(i,kso) + (1.0_wp-zxx)*0.06_wp
        ENDIF
      ENDDO
    ENDDO


    IF (lterra_urb) THEN
      ! HW: Modification of the surface thermal conductivity:
      ! - according to the building materials (urb_hcon)
      ! - area index of buildings (urb_ai), height of building elements (urb_h_bld)
      !   and building fraction (urb_isa)
      !
      ! Because of the curvature of the surface, the uppermost soil layer heat
      ! transfer is larger compared to the thermal conductivity of a plan area.
      ! As a result, the effective thermal conductivity of the upper surface is increased.
      ! This is also the case in the layers beneath, in which the effective thermal conductivity
      ! decreases with depth.
      !
      ! This modification decreases with depth with respect to the
      ! natural soil below the buildings.

      !$ACC LOOP SEQ
      DO kso = 1, ke_soil+1
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zalpha_uf)
        DO i = ivstart, ivend

          zalpha_uf     = MAX(0.0_wp, MIN(zmls(kso)/urb_h_bld(i), 1.0_wp))

          hzalam(i,kso) =           urb_isa(i)  * ( (1.0_wp - zalpha_uf) * urb_hcon(i)*urb_ai(i)        &
                                                    +         zalpha_uf  * hzalam(i,kso)         )      &
                        + (1.0_wp - urb_isa(i)) * hzalam(i,kso)

        ENDDO
      ENDDO
    END IF


!$NEC nofuse
!$NEC outerloop_unroll(7)
    !$ACC LOOP SEQ
    DO kso = 1, ke_soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        ! interpolated thermal conductivity on half levels, meaning at the boundaries between soil layers
        zalam(i,kso) = hzalam(i,kso)*hzalam(i,kso+1)*(zmls(kso+1)-zmls(kso))                            &
                       / ( hzalam(i,kso)*(zmls(kso+1)-zzhls(kso))                                       &
                       +   hzalam(i,kso+1)*(zzhls(kso)-zmls(kso)) )
      ENDDO
    ENDDO

  ELSE

! heat conductivity based on assumption of a soil water content which is equal to the
! average between wilting point and field capacity
    !$ACC LOOP SEQ
!$NEC outerloop_unroll(7)
    DO kso = 1, ke_soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zwqg, z4wdpv)
      DO i = ivstart, ivend
        zwqg         = 0.5_wp*(zfcap(i,kso) + zpwp(i,kso))
        z4wdpv       = 4._wp*zwqg/zporv(i,kso)
        ! heat conductivity
        zalamtmp(i,kso) =              zdlam(i)                         &
                      * (0.25_wp + 0.30_wp*zdlam(i)           &
                      / (1._wp+0.75_wp*zdlam(i)))             &
                      * MIN (z4wdpv, 1.0_wp + (z4wdpv-1.0_wp)   &
                      *(1.0_wp+0.35_wp*zdlam(i))              &
                      /(1.0_wp+1.95_wp*zdlam(i)))
      ENDDO
    ENDDO

    !$ACC LOOP SEQ
    DO kso = 1, ke_soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        zalam(i,kso) = zalam(i,kso) + zalamtmp(i,kso)
        hzalam(i,kso) = zalam(i,kso)
      ENDDO
    ENDDO
  ENDIF


! Initialisations and conversion of tch to tmch
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    limit_tch(i) = .false.  !  preset problem indicator
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zuv, ztvs, zplow, zxx, zdT_snow)
  DO i = ivstart, ivend
    zuv        = SQRT ( u(i)**2 + v(i)**2 )
    ztvs       = t_g (i)*(1.0_wp + rvd_m_o*qv_s(i))

    !  'potential' temperature of lowest atmospheric layer
    zplow        =  ptot(i)
    zth_low (i)  =  t(i) * EXP(rdocp*LOG(ps(i)/zplow))

    zdt_atm (i)  =  zth_low(i)-t_g(i)
    zdq_atm (i)  =  qv(i)-qv_s(i)

    ! introduce an artificical upper boundary on transfer coefficients in cases
    ! where an extreme cooling/heating of topmost soil layer may be caused due
    ! to excessive sensible/latent heat fluxes (e.g. after creation of unbalanced
    ! structures in the data assimilation or following massive changes in radiative
    ! forcing due to infrequent radiation computations)

    ! estimate current energy budget of topmost soil layer

    ! heat flux between layers 1&2 based on current temperature profile
    zg1(i)= zalam(i,1)*(t_so_now(i,1)-t_so_now(i,2))/zdzms(2)

    ! estimates of sensible and latent heat flux
    zrho_atm(i)=ps(i)/(r_d*ztvs)
    zshfl(i) = tch(i)*zuv*zrho_atm(i)*cp_d*zdt_atm(i)
    zlhfl(i) = tch(i)*zuv*zrho_atm(i)*lh_v*zdq_atm(i)

    ! net radiative fluxes at surface
    zradfl(i) = sobs(i)+thbs(i)

    zxx = MIN(500.0_wp,200.0_wp+0.5_wp*ABS(zradfl(i)))

    IF (zshfl(i)*zlhfl(i) >= 0._wp) THEN
      zthfl(i) = zshfl(i) + zlhfl(i)
    ELSE IF (ABS(zshfl(i)) > ABS(zlhfl(i))) THEN
      zthfl(i) = zshfl(i) + SIGN(MIN(zxx,ABS(zlhfl(i))),zlhfl(i))
    ELSE
      zthfl(i) = zlhfl(i) + SIGN(MIN(zxx,ABS(zshfl(i))),zshfl(i))
    ENDIF

    IF (ABS(zthfl(i)) <= eps_soil) zthfl(i)=SIGN(eps_soil,zthfl(i))

    ! unconstrained estimated energy budget of topmost soil layer
    zeb1(i) = zthfl(i)+zradfl(i)-zg1(i)

    ! energy required to melt existing snow
    ze_melt(i)=w_snow_now(i)*rho_w*lh_f     ! (J/m**2)

    ! heat capacity of snow layer, limited to a snow depth of 1.5 m 
    ! for consistency with subsequent calculations
    zch_snow(i)=MIN(w_snow_now(i),1.5_wp*rho_snow_now(i)/rho_w)*rho_w*chc_i   ! (J/(m**2 K))

    ! constrain transfer coefficient, if energy budget  of topmost soil layer is:
    ! a) negative & surface layer is unstable (i.e   upward directed turbulent heat flux)
    ! b) positive & surface layer is stable   (i.e downward directed turbulent heat flux)

    IF (zeb1(i)<0.0_wp .AND. zthfl(i)<0.0_wp) THEN
      ! cooling of 1st soil layer&upward SHF+LHF

      ztchv_max(i) = ( zlim_dtdt*(zch_soil*zdzhs(1)+zch_snow(i))/zdt    &
                   + ABS(zg1(i)-zradfl(i)) ) / ABS(zthfl(i)) * tch(i)*zuv
    ELSEIF (zeb1(i)>0.0_wp .AND. zthfl(i)>0.0_wp) THEN
      ! heating of 1st soil layer & downward SHF+LHF
      !   Note: The heat capacity of snow is only relevant for the difference
      !         between the actual temperature and the melting point. The mean
      !         snow temperature is set to the average of t_snow & t_so(1)
      IF (lmulti_snow) THEN
        zdT_snow=MIN(0._wp, t0_melt-0.5_wp*(t_snow_mult_now(i,1)+t_so_now(i,1)))
      ELSE
        zdT_snow=MIN(0._wp,t0_melt-0.5_wp*(t_snow_now(i)+t_so_now(i,1)))
      ENDIF
      ztchv_max(i) = ( (zlim_dtdt*zch_soil*zdzhs(1)+zdT_snow*zch_snow(i)+ze_melt(i))/zdt  &
                   + ABS(zg1(i)-zradfl(i)) ) / zthfl(i) * tch(i)*zuv
    ELSE
      ! unlimited transfer coefficient
      ztchv_max(i) = HUGE(1._wp)
    ENDIF
                                                    ! required constraint as non-turbulent
                                                    ! energy budget components alone may
    ztchv_max(i) =MAX( ztchv_max(i), eps_soil)      ! exceed the upper limit in the energy
                                                    ! budget of the 1st soil layer

    ! Additional limitation for better numerical stability at long time steps
    ztchv_max(i) = MIN(ztchv_max(i),(4._wp*zlim_dtdt*(zch_soil*zdzhs(1)+zch_snow(i))/zdt &
                  +ABS(zradfl(i)))/ABS(zthfl(i))*tch(i)*zuv)

    ztchv(i)    = tch(i)*zuv  ! transfer coefficient * velocity

    LIM: IF (ztchv(i) > ztchv_max(i)) THEN
      tch(i)=ztchv_max(i)/MAX(zuv,1.E-06_wp)
      limit_tch(i) = .true.          ! set switch for later use
    END IF LIM

    ztmch(i) = tch(i)*zuv*g*zrho_atm(i)
  ENDDO


#ifndef _OPENACC
  IF (msg_level >= 20) THEN
    ! counter for limitation of transfer coefficients
    m_limit = COUNT( limit_tch(:) )

    ! In debugging mode and if transfer coefficient occured for at least one grid point
    IF (m_limit > 0) THEN
      WRITE(*,'(1X,A,/,2(1X,A,F10.2,A,/),1X,A,F10.2,/,1X,A,F10.3,/)')                  &
             'terra1: transfer coefficient had to be constrained',                     &
             'model time step                                 :', zdt     ,' seconds', &
             'max. temperature increment allowed per time step:',zlim_dtdt,' K',       &
             'upper soil model layer thickness                :', zdzhs(1)

      DO i = ivstart, ivend
        IF (limit_tch(i)) THEN
          zuv = SQRT (u(i)**2 + v(i)**2 )
          WRITE(*,*) 'TERRA flux limiter: TCH before and after, zshfl, zlhf, Tsfc, Tatm ', &
                      ztchv(i)/zuv, tch(i), zshfl(i), zlhfl(i), t_g(i), t(i)

          yhc = 'COOLING'
          IF (zeb1(i) > 0._wp) yhc='HEATING'
        END IF
      END DO
    ENDIF
  ENDIF
#endif

  ! Update indicator for age of snow in top of snow layer
  ! Note that cloud ice is deliberately disregarded here in order to avoid counting drifting snow as fresh snow
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zrain_rate, zsnow_rate, zdsn_new, zdsn_old, ztau_snow, zuv)
  DO i = ivstart, ivend
    IF (w_snow_now(i) <=0.0_wp) THEN
      ! if no snow exists, reinitialize age indicator
      freshsnow(i) = 1.0_wp
    ELSE
      IF ( nclass_gscp >= 2000 ) THEN
        ! only possible when running 2-moment microphysics
#ifdef TWOMOM_SB
        zsnow_rate = prs_gsp(i)+prs_con(i)+prg_gsp(i)+prh_gsp(i) ! [kg/m**2 s]
#endif
      ELSEIF ( nclass_gscp >= 6 ) THEN
        zsnow_rate = prs_gsp(i)+prs_con(i)+prg_gsp(i)            ! [kg/m**2 s]
      ELSE
        zsnow_rate = prs_gsp(i)+prs_con(i)                       ! [kg/m**2 s]
      ENDIF

      zrain_rate = prr_gsp(i)+prr_con(i)  ! [kg/m**2 s]

      ! temperature-dependent aging timescale: 7 days at freezing point, 28 days below -15 deg C
      ztau_snow = 86400._wp*MIN(28.0_wp,7._wp+1.4_wp*(t0_melt-MIN(t0_melt,t_snow_now(i))))

      ! wind-dependent snow aging: a thin snow cover tends to get broken under strong winds, which reduces the albedo
      ! an offset is added in order to ensure moderate aging for low snow depths
      zuv = MIN(300._wp, u_10m(i)**2 + v_10m(i)**2 + 12._wp )
      ztau_snow = MIN(ztau_snow,MAX(86400._wp,2.e8_wp*MAX(0.05_wp,h_snow_gp(i))/zuv))

      ! decay rate for fresh snow including contribution by rain (full aging after 10 mm of rain)
      zdsn_old   = zdt/ztau_snow + zdt*zrain_rate*0.1_wp

      ! linear growth rate equals 1.0 in 1 day for a temperature-dependent snow rate between
      ! 10 mmH2O (kg/m**2) per day (0.1) and 5 mmH2O (kg/m**2) per day (0.2)
      zdsn_new   = zdt*zsnow_rate*(0.1_wp + MIN(0.1_wp,0.02_wp*(t0_melt-t(i))))

      ! reduce decay rate, if new snow is falling and as function of snow
      ! age itself
      zdsn_old   = (zdsn_old - zdsn_new)*freshsnow(i)
      zdsn_old   = MAX(zdsn_old,0._wp)

      freshsnow(i) = freshsnow(i) + zdsn_new-zdsn_old

      freshsnow(i) = MIN(1._wp,MAX(0._wp,freshsnow(i)))

    END IF
  ENDDO

!------------------------------------------------------------------------------
! Section I.2: temperatures, water contents (in mH2O), surface pressure,
!------------------------------------------------------------------------------

  IF (lmulti_snow) THEN
    !$ACC LOOP SEQ
    DO ksn = 0,ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (w_snow_now(i) > 0.0_wp) THEN
          ! existence of snow
          t_snow_mult_now(i,ksn) = MIN (t0_melt - eps_temp, t_snow_mult_now(i,ksn) )
        ELSE IF (t_snow_mult_now(i,ke_snow) >= t0_melt) THEN
          ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
          t_snow_mult_now(i,ksn) = MAX (t0_melt + eps_temp, t_s_now(i) )
        ELSE
          ! no snow and  t_snow < t0_melt
          ! --> t_snow = t_s
          t_snow_mult_now(i,ksn) = MIN (t0_melt - eps_temp, t_s_now(i) )
        END IF
      ENDDO
    ENDDO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
       IF (w_snow_now(i) > 0.0_wp) THEN
         ! existence of snow
         ! --> no water in interception store and t_snow < t0_melt
         w_snow_now(i) = w_snow_now(i) + w_i_now(i)
         wtot_snow_now(i,1) = wtot_snow_now(i,1) + w_i_now(i)
         dzh_snow_now(i,1)  = dzh_snow_now(i,1)  + w_i_now(i) / rho_snow_mult_now(i,1)*rho_w
         w_i_now(i)         = 0.0_wp
         h_snow_now(i)      = h_snow(i)
       ELSE IF (t_snow_mult_now(i,ke_snow) >= t0_melt) THEN
         ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
         t_s_now   (i) = t_snow_mult_now(i,ke_snow)
         h_snow_now(i) = 0.0_wp
       ELSE
         ! no snow and  t_snow < t0_melt
         ! --> t_snow = t_s and no water w_i in interception store
         t_s_now   (i) = t_snow_mult_now(i,ke_snow)
         w_i_now(i)    = 0.0_wp
         h_snow_now(i) = 0.0_wp
       END IF
    ENDDO

  ELSE ! no multi-layer snow

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (w_snow_now(i) > 0.0_wp) THEN
        ! existence of snow
        ! --> no water in interception store and t_snow < t0_melt
        ! GZ: this effectively suppresses rime formation because deposition rates per time
        ! step are usually less than 1e-6 m (eps_temp, eps_soil)
   !!!  w_snow_now(i) = w_snow_now(i) + w_i_now(i)
   !!!  w_i_now(i) = 0.0_wp
        t_snow_now(i) = MIN (t0_melt - eps_temp, t_snow_now(i) )
      ELSE IF (t_snow_now(i) >= t0_melt) THEN
        ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
        t_s_now   (i) = MAX (t0_melt + eps_temp, t_s_now(i) )
        t_snow_now(i) = t_s_now(i)
      ELSE
        ! no snow and  t_snow < t0_melt
        ! --> t_snow = t_s and no water w_i in interception store
        t_s_now   (i) = MIN (t0_melt - eps_temp, t_s_now(i) )
        t_snow_now(i) = t_s_now(i)
   !!!  w_i_now(i) = 0.0_wp
      END IF
    ENDDO
  ENDIF  ! lmulti_snow

  ! Initializations for the next sections

  IF (lmulti_snow) THEN
    ! some preparations for ksn==0 and ksn==1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      ztsnow_mult   (i,0) = t_snow_mult_now(i,0)
      ztsnow_mult   (i,1) = t_snow_mult_now(i,1)

      zhh_snow(i,1) =  -h_snow_now(i) + dzh_snow_now(i,1)
      zhm_snow(i,1) = (-h_snow_now(i) + zhh_snow(i,1))/2._wp

      zdzh_snow(i,1) = dzh_snow_now(i,1)
      zextinct (i,1) = 0.13_wp*rho_snow_mult_now(i,1)+3.4_wp

      ! set ztsnow to ztsnow_mult(ksn=1)
      ztsnow   (i) = t_snow_mult_now(i,1)

      ! initialize zdz_snow (from Section I.3) to zdzh_snow(i,1)
      zdz_snow (i) = dzh_snow_now(i,1) ! zdzh_snow
      zwsnow   (i) = w_snow_now(i)
    ENDDO


    !$ACC LOOP SEQ
    DO  ksn = 2,ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        ztsnow_mult(i,ksn) = t_snow_mult_now(i,ksn)

        zhh_snow   (i,ksn) = zhh_snow(i,ksn-1) + dzh_snow_now(i,ksn)
        zhm_snow   (i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._wp

        zdzh_snow  (i,ksn) = dzh_snow_now(i,ksn)
        zextinct   (i,ksn) = 0.13_wp*rho_snow_mult_now(i,ksn)+3.4_wp

        ! build sum over all layers for zdzh_snow in zdz_snow
        zdz_snow   (i)     = zdz_snow(i) + zdzh_snow(i,ksn)
      ENDDO
    ENDDO

  ELSE  ! no lmulti_snow

    ! set ztsnow to t_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      ztsnow   (i) = t_snow_now(i)
      zwsnow   (i) = w_snow_now(i)
      zdz_snow (i) = zwsnow(i)*rho_w/rho_snow_now(i)
    ENDDO
  ENDIF  ! lmulti_snow

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(z2iw, z4iw, z234iw, zqs, zdqs, zqsnow, zdqsnow, zplow, zqvlow)
  DO i = ivstart, ivend
    ! ztsnow is now set according to lmulti_snow (see above);
    ! therefore we need not take care of it in this loop
    ! ztsnow   (i) = t_snow(i,nx)
    zts       (i) = t_s_now   (i)
    IF (itype_canopy == 1) THEN
      ztsk    (i) = t_s_now  (i)
    ELSE IF (itype_canopy == 2) THEN
      ztsk    (i) = t_sk_now  (i)
    END IF
    zts_pm    (i) = zsf_heav(zts   (i) - t0_melt)
    ztsk_pm   (i) = zsf_heav(ztsk  (i) - t0_melt)
    ztfunc    (i) = MAX(0.0_wp,1.0_wp - MAX(0.0_wp,0.5_wp*(ztsk(i)-t0_melt)))
    ztsnow_pm (i) = zsf_heav(ztsnow(i) - t0_melt)

    IF (itype_interception == 1) THEN
      zwin    (i) = w_i_now(i)
      w_p_now (i) =  0._wp
      w_p_new (i) =  0._wp
      w_s_now (i) =  0._wp
      w_s_new (i) =  0._wp
    ELSE IF (itype_interception == 2) THEN
      zwin    (i) = w_i_now   (i)
      zwpn    (i) = w_p_now   (i)
      zwisn   (i) = w_s_now   (i)
    END IF

    ! moisture and potential temperature of lowest atmospheric layer
    zplow       = ptot(i)
    zqvlow      =   qv(i)
    zth_low (i) =    t(i) * EXP(rdocp*LOG(ps(i)/zplow))

    ! density*transfer coefficient*wind velocity
    zrhoch(i)   = ztmch(i)*(1._wp/g) + eps_soil

    ! saturation specific humidity for t_s and t_snow and first derivative
    z2iw      = ztsk_pm(i)*b2w + (1._wp - ztsk_pm(i))*b2i
    z4iw      = ztsk_pm(i)*b4w + (1._wp - ztsk_pm(i))*b4i
    zqs       = zsf_psat_iw(ztsk(i), z2iw,z4iw)
    zqs       = zsf_qsat( zqs, ps(i) )
    zdqs      = zqvlow - zqs
    IF (ABS(zdqs).LT.0.01_wp*eps_soil) zdqs = 0.0_wp
    z2iw      = ztsnow_pm(i)*b2w + (1._wp - ztsnow_pm(i))*b2i
    z4iw      = ztsnow_pm(i)*b4w + (1._wp - ztsnow_pm(i))*b4i
    z234iw    = z2iw*(b3 - z4iw)
    zqsnow    = zsf_psat_iw(ztsnow(i)-MAX(0.0_wp,tsnred(i)),z2iw,z4iw)
    zqsnow    = zsf_qsat(zqsnow, ps(i))

    zdqvtsnow(i)= zsf_dqvdt_iw(ztsnow(i), zqsnow, z4iw,z234iw)
    zdqsnow     = zqvlow - zqsnow
    IF (ABS(zdqsnow).LT.0.01_wp*eps_soil) zdqsnow = 0.0_wp

    ! potential evaporation at T_snow and Ts
    zep_snow(i) = (1._wp-ztsnow_pm(i))*tfv(i)*zrhoch(i)*zdqsnow*MERGE(tfvsn(i), 1._wp, zdqsnow<0._wp)
    zep_s   (i) =                      tfv(i)*zrhoch(i)*zdqs
  ENDDO

!------------------------------------------------------------------------------
! Section I.3: heat conductivity, frozen fraction, snow and
!            water covered fraction, snow density and  height,
!            volumetric heat content
!------------------------------------------------------------------------------

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zzz, zrww)
  DO i = ivstart, ivend
     ! snow and water covered fraction
     !em        zrss = MAX( 0.01_wp, MIN(1.0_wp,zwsnow(i)/cf_snow) )
     zzz   = MAX( 0.25_wp*cf_w,0.4_wp*cwimax_ml*MAX(2.5_wp*plcov(i),tai(i)) )

     ! TERRA_URB: Puddles on the impervious surface area
     IF (lterra_urb .AND. (itype_eisa == 3)) THEN
       zzz = urb_isa(i) * cwisamax + (1.0_wp - urb_isa(i)) * zzz
     END IF

     zrww  = MAX( 0.01_wp, 1.0_wp - EXP(MAX( -5.0_wp, - zwin(i)/zzz) ) )
     !em        zf_snow(i) = zrss*zsf_heav(zwsnow(i) - eps_soil)

     IF (itype_interception == 1) THEN
       zf_wi  (i) = zrww*zsf_heav(zwin  (i) - 1.0E-4_wp*eps_soil)
     ELSE IF (itype_interception == 2) THEN
       zf_wi  (i) = plcov (i) !Fraction of interception store on grid box area scales with plcov
     END IF

     ! BR 7/2005 prognostic  snow density
     ! US DMironov: for the FLake Model, h_snow has to be a prognostic variable but
     !              here it is used only diagnostically. So the values are put
     !              to every time level

     ! zdz_snow has been computed in the initializations above depending
     ! on lmulti_snow
     ! zdz_snow(i)=zwsnow(i)*rho_w/rho_snow(i,nx)
     h_snow_new(i) = zdz_snow(i)
     h_snow_now(i) = zdz_snow(i)
     ! IF (.NOT. l2tls) h_snow(i,nold) = zdz_snow(i)

     ! constrain snow depth and consider this constraint for the computation
     ! of average snow density of snow layer
     zf_snow_lim(i) = MAX(0.01_wp,0.1_wp*freshsnow(i),zf_snow(i))
     zdz_snow (i) =  zdz_snow(i)/zf_snow_lim(i)
     zdz_snow (i) =  MAX(cdsmin,zdz_snow(i))

     ! limitation of snow depth to 1.5m for snow cover heat transfer
     zdz_snow_fl(i) = MIN(1.5_wp, zdz_snow(i))
     IF (.NOT. lmulti_snow) zrocs(i) = chc_i*zdz_snow_fl(i)*zrho_snow(i)
  ENDDO

!------------------------------------------------------------------------------
! Section I.4: Hydrology, 1.Section
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section I.4.1: Evaporation from interception store and from snow cover,
  !----------------------------------------------------------------------------
  ! Evaporation and transpiration are negative, dew and rime
  ! positive quantities, since positive sign indicates a flux
  ! directed towards the earth's surface!

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zzz)
  DO i = ivstart, ivend
    ! Evaporation from interception store if it contains water (wi>0) and
    ! if zep_s<0 indicates potential evaporation for temperature Ts
    ! amount of water evaporated is limited to total content of store
    zzz = (1.0_wp + 0.5_wp*ztfunc(i))/3.0_wp
    zdwidt(i) = zsf_heav(-zep_s(i)) * MAX(-zrhwddt*zwin(i),              &
      zzz*zf_wi(i)*zep_s(i), -MAX(300.0_wp,0.75_wp*zradfl(i))/lh_v)

    ! Evaporation of snow, if snow exists (wsnow>0) and if zep_snow<0
    ! indicates potential evaporation for temperature t_snow
    zdwsndt(i) = zsf_heav(-zep_snow(i))  &
                       * MAX(-zrhwddt*zwsnow(i), zf_snow(i)*zep_snow(i))

    ! Formation of dew or rime, if zep_s > 0 . distinction between
    ! dew or rime is only controlled by sign of surface temperature
    ! and not effected by presence of snow !
    zrr(i)=zsf_heav(zep_s   (i))*zep_s   (i)*        ztsk_pm(i)
    zrs(i)=zsf_heav(zep_snow(i))*zep_snow(i)*(1.0_wp-ztsk_pm(i))
  ENDDO

  IF (itype_interception == 2) THEN

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      zwimax(i) = MAX(1.E-6_wp,4.E-4_wp * sai(i)) ! Security min. value 1E-6 m
      zfd_wi(i)=MIN(1._wp,MAX(0.0_wp, EXP((2._wp/3._wp)*LOG(zwin(i)/zwimax(i))) ))
      zf_pd(i) = MIN(1._wp,MAX(0.0_wp, EXP((2._wp/3._wp)*LOG(zwpn(i)/zpercmax)) ))

      zewi(i)=zsf_heav(-zep_s(i)) * zf_wi(i)*zfd_wi(i)*zep_s(i) ! canopy covered part
      zepd(i)=zsf_heav(-zep_s(i)) * (1._wp - zf_wi(i))*zf_pd(i)*zep_s(i) ! bare soil part
      zept(i)=zsf_heav(-zep_s(i)) * zep_s(i) ! potential evaporation
      zesn(i)=zsf_heav(-zep_snow(i)) * zep_snow(i) ! Snow evaporation
      zdrr(i)=zsf_heav(zep_s   (i))*zep_s   (i)*        ztsk_pm(i)
      zrrs(i)=zsf_heav(zep_snow(i))*zep_snow(i)*(1.0_wp-ztsk_pm(i))
    END DO

  END IF

  !----------------------------------------------------------------------------
  ! Section I.4.2b: Bare soil evaporation, BATS version
  !----------------------------------------------------------------------------

  IF (itype_evsl.EQ.2) THEN

    ! Calculation of bare soil evaporation after Dickinson (1984)
    ! Determination of mean water content relative to volume of voids
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zice, zevap, zbeta, zbf1, zbf2, zdmax, zd) &
    !$ACC   PRIVATE(zck, zfqmax, zevapor)
    DO i = ivstart, ivend
      IF (zep_s(i) < 0.0_wp) THEN   ! upwards directed potential evaporation
        zsnull(i) = zsnull(i)/(znull*zporv(i,1))

        ! Treatment of ice (m_styp=1) and rocks (m_styp=2)
        zice   = zsf_heav(1.5_wp - REAL(m_styp(i),wp)) ! 1 only for ice
        zevap  = zrock(i) + zice                       ! 1 for all soil types
                                                       ! but rock and ice (=0)
        zbeta  = 0.0_wp

        IF (m_styp(i).ge.3) THEN ! Computations not for ice and rocks
          ! auxiliary quantities
          zbf1   = 5.5_wp - 0.8_wp* zbedi(i)*                &
                  (1.0_wp + 0.1_wp*(zbedi(i) - 4.0_wp)*  &
                   clgk0(m_styp(i)) )
          zbf2   = (zbedi(i) - 3.7_wp + 5.0_wp/zbedi(i))/  &
                  (5.0_wp + zbedi(i))
          zdmax  = zbedi(i)*cfinull*zk0di(i)/crhowm
          zs1(i)  = zs1(i)/(z1*zporv(i,1))
          zd     = 1.02_wp*zdmax*EXP( (zbedi(i)+2._wp)*LOG(zs1(i)) ) * &
                                     EXP( zbf1*LOG(zsnull(i)/zs1(i)) )
          zck    = (1.0_wp + 1550.0_wp*cdmin/zdmax)*zbf2
          ! maximum sustainable moisture flux in the uppermost surface
          ! layer in kg/(s*m**2)
          zfqmax = - rho_w*zck*zd*zsnull(i)/SQRT(znull*z1)
          zevapor= MAX(zep_s(i),zfqmax)
          IF(zw_fr(i,1)+zevapor*(1.0_wp - zf_wi(i))   &
                                 *(1.0_wp - zf_snow(i)) &
                                 *eai(i)/sai(i)* zdtdrhw/zdzhs(1) &
                                 .LE.zadp(i,1)) zevapor = 0._wp
          zbeta  = zevapor/MIN(zep_s(i),-eps_div)
        END IF ! Computations not for ice and rocks

        zbeta  = zbeta + (1.0_wp - zbeta)*zice
        ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other soil types
        zsfc_frac_bs(i)=eai(i)/sai(i)

        IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) THEN      ! AYu mire block
          zbeta           = 0.6_wp
          zsfc_frac_bs(i) = 1.0_wp
        ENDIF

        ! consideration of plant or snow/water cover
        IF (itype_interception == 1) THEN
          zesoil(i) = zevap*zbeta*zep_s(i)        & ! evaporation
          !!!              *(1.0_wp - zf_wi  (i)) & ! not water covered
                           *(1.0_wp - zf_snow(i)) & ! not snow covered
                           * zsfc_frac_bs(i)        ! relative source surface
                                                    ! of the bare soil

        ELSE IF (itype_interception == 2) THEN
          zesoil(i) = zevap*zbeta*zep_s(i)       & ! evaporation
                          *(1.0_wp - plcov(i))   & ! plant cover weighting
                          *(1.0_wp - zf_snow(i)) & ! not snow covered
                          *(1.0_wp - zf_pd(i))     ! not pond covered
        END IF ! interception

        lhfl_bs(i) = lh_v * zesoil(i)
        
      END IF  ! upwards directed potential evaporation
    END DO
  END IF ! BATS version

  !----------------------------------------------------------------------------
  ! Section I.4.2c: Bare soil evaporation, Noilhan and Planton, 1989
  !----------------------------------------------------------------------------

  IF (itype_evsl.EQ.3) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zice, zevap, zbeta, zalpha, z2iw, z4iw) &
    !$ACC   PRIVATE(zqs, zevapor)
    DO i = ivstart, ivend
      IF (zep_s(i) < 0.0_wp) THEN   ! upwards directed potential evaporation
        zsnull(i) = zsnull(i)/(znull*zporv(i,1))

        ! Treatment of ice (m_styp=1) and rocks (m_styp=2)
        zice   = zsf_heav(1.5_wp - REAL(m_styp(i),wp)) ! 1 only for ice
        zevap  = zrock(i) + zice                       ! 1 for all soil types
                                                       ! but rock and ice (=0)
        zbeta  = 0.0_wp

        IF (m_styp(i).ge.3) THEN ! Computations not for ice and rocks
          IF (zw_fr(i,1)> zfcap(i,1)) THEN
            zalpha = 1.0_wp
          ELSE
            zalpha = 0.5_wp * (1.0_wp - COS ( pi *                  &
                     (zw_fr(i,1) - zadp(i,1)) / ( zfcap(i,1) - zadp(i,1)) ) )
          ENDIF
          IF (itype_canopy == 1) THEN
            z2iw   = ztsnow_pm(i)*b2w + (1._wp - ztsnow_pm(i))*b2i
            z4iw   = ztsnow_pm(i)*b4w + (1._wp - ztsnow_pm(i))*b4i
            zqs    = zsf_qsat( zsf_psat_iw(zts(i), z2iw,z4iw), ps(i) )
          ELSE IF (itype_canopy == 2) THEN
            z2iw   = ztsk_pm(i)*b2w + (1.0_wp - ztsk_pm(i))*b2i
            z4iw   = ztsk_pm(i)*b4w + (1.0_wp - ztsk_pm(i))*b4i
            zqs    = zsf_qsat( zsf_psat_iw(ztsk(i), z2iw,z4iw), ps(i) )
          END IF
          zevapor= MIN(0.0_wp,zrhoch(i)*(qv(i)-zalpha*zqs))

          zbeta  = zevapor/MIN(zep_s(i),-eps_div)
        END IF ! Computations not for ice and rocks

        zbeta  = zbeta + (1.0_wp - zbeta)*zice
        ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other
        ! soil types
        ! consideration of plant or snow/water cover

        IF (itype_interception == 1) THEN
          zesoil(i) = zevap*zbeta*zep_s(i)       & ! evaporation
           !!!            *(1.0_wp - zf_wi  (i)) & ! not water covered
                          *(1.0_wp - zf_snow(i)) & ! not snow covered
                          * eai(i)/sai(i)          ! relative source surface
                                                   ! of the bare soil

        ELSE IF (itype_interception == 2) THEN
           zesoil(i) = zevap*zbeta*zep_s(i)      & ! evaporation
                          *(1.0_wp - plcov(i))   & ! plant cover weighting
                          *(1.0_wp - zf_snow(i)) & ! not snow covered
                          *(1.0_wp - zf_pd(i))     ! not pond covered
        END IF ! interception

        lhfl_bs(i) = lh_v * zesoil(i)

      END IF  ! upwards directed potential evaporation
    END DO
  END IF ! NP89

  !----------------------------------------------------------------------------
  ! Section I.4.2d: Bare soil evaporation, resistance version
  !----------------------------------------------------------------------------

  IF (itype_evsl == 4 .OR. itype_evsl == 5) THEN   ! Resistance version
    ! Calculation of bare soil evaporation using a resistance formulation.
    ! For a review see Schulz et al. (1998) and Schulz and Vogel (2020).
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zice, zevap, zbeta, zalpha)
    DO i = ivstart, ivend
      IF (zep_s(i) < 0.0_wp) THEN   ! upwards directed potential evaporation
        ! Treatment of ice (m_styp=1) and rocks (m_styp=2)
        zice  = zsf_heav(1.5_wp - REAL(m_styp(i),wp)) ! 1 only for ice
        zevap = zrock(i) + zice                       ! 1 for all, but rock
        zbeta = 0.0_wp

        IF (m_styp(i).GE.3) THEN ! Computations not for ice and rocks
           zalpha = MAX( 0.0_wp, MIN( 1.0_wp,                     &
                    (zw_fr(i,1) - zadp(i,1)) / (zfcap(i,1) - zadp(i,1)) ) )
           zalpha = r_bsmin(i) / (zalpha + eps_soil)
           zbeta  = 1.0_wp                                &
                  / (1.0_wp + zrhoch(i)*zalpha/zrho_atm(i))
        END IF ! Computations not for ice and rocks

        zbeta  = zbeta + (1.0_wp - zbeta)*zice
        ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other soil types
        zsfc_frac_bs(i)=eai(i)/sai(i)

        IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) THEN      ! AYu mire block
          zbeta           = 0.6_wp
          zsfc_frac_bs(i) = 1.0_wp
        ENDIF

        ! Consideration of plant or snow/water cover
        IF (itype_interception == 1) THEN          ! Interception
           zesoil(i) = zevap*zbeta*zep_s(i)      & ! evaporation
           !!!            *(1.0_wp - zf_wi  (i)) & ! not water covered
                     * (1.0_wp - zf_snow(i))     & ! not snow covered
                          * zsfc_frac_bs(i)        ! relative source surface
                                                   ! of the bare soil
        ELSE IF (itype_interception == 2) THEN
           zesoil(i) = zevap*zbeta*zep_s(i)      & ! evaporation
                     *(1.0_wp - plcov(i))        & ! plant cover weighting
                     *(1.0_wp - zf_snow(i))      & ! not snow covered
                     *(1.0_wp - zf_pd(i))          ! not pond covered
        END IF ! Interception

        lhfl_bs(i) = lh_v * zesoil(i)

      END IF  ! upwards directed potential evaporation
    END DO

  END IF      ! Resistance version

  !----------------------------------------------------------------------------
  ! Section I.4.3b: transpiration by plants, BATS version
  !----------------------------------------------------------------------------

  BATS: IF (itype_trvg == 2 .OR. itype_trvg == 3) THEN   ! BATS version
    ! This version is based on Dickinson's (1984) BATS scheme, simplified by
    ! neglecting the water and energy transports between the soil and the plant
    ! canopy. This leads to a Monteith combination formula for the computation
    ! of plant transpiration.
    ! Option 3 is an extended variant with an additional diagnostic variable for 
    ! accumulated plant evaporation since sunrise, allowing for a better representation of
    ! the diurnal cycle of plant evaporation (in particular trees)

    ! Root distribution

    IF (itype_root == 2) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        zrootdz_int (i)= 0.0_wp   !initialize the root density profile integral
        zwrootdz_int(i)= 0.0_wp   !initialize the root water content integral
      END DO


      !$ACC LOOP SEQ
      DO kso   = 1,ke_soil
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zrootfr, zrootdz)
        DO i = ivstart, ivend
          IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_wp) THEN
            ! consider the effect of root depth & root density
            zrootfr = EXP (-zroota(i)*zmls(kso)) ! root density
            zrootdz = zrootfr*MIN(zdzhs(kso),MAX(0.0_wp, zbwt(i)-(zmls(kso) &
                       -0.5_wp*zdzhs(kso)) ) )
            zrootdz_int(i)=zrootdz_int(i) + zrootdz
            ! The factor of 10 ensures that plants do not extract notable amounts of water from partly frozen soil
            zwrootdz(i,kso)=zrootdz*MAX(zpwp(i,kso),zw_fr(i,kso)-10.0_wp*w_so_ice_now(i,kso)/zdzhs(kso))
            zwrootdz_int(i)=zwrootdz_int(i) + zwrootdz(i,kso)
          END IF  ! negative potential evaporation only
        END DO
      END DO

      ! Compute root zone integrated average of liquid water content
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        zwrootdz_int(i)=zwrootdz_int(i)/MAX(zrootdz_int(i),eps_div)
      END DO

    ELSE   ! itype_root

      !$ACC LOOP SEQ
      DO kso   = 1,ke_soil
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zropart)
        DO i = ivstart, ivend
          IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_wp) THEN 
             ! upwards directed potential evaporation
              zropart  = MIN ( zdzhs(kso), MAX(0.0_wp,                     &
                               zbwt(i) - (zmls(kso) - 0.5_wp*zdzhs(kso))))
              zropartw(i,kso) = zropart*(zw_fr(i,kso)-w_so_ice_now(i,kso)/zdzhs(kso))
              zwroot(i) = zwroot(i) + zropartw(i,kso)/MAX(eps_div,zbwt(i))
          END IF  ! negative potential evaporation only
        END DO
      END DO

    ENDIF  ! itype_root


    ! Determination of the transfer functions CA, CF, and CV
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zuv, zcatm, zustar, zrla, zpar, zf_wat) &
    !$ACC   PRIVATE(zf_tem, z2iw, z4iw, zepsat, zepke, zf_sat, zedrstom) &
    !$ACC   PRIVATE(zrstom, zrveg, zxx, zzz)
    DO i = ivstart, ivend
      IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_wp) THEN
        ! upwards directed potential evaporation
        zuv        = SQRT (u(i)**2 + v(i)**2 )
        zcatm      = tch(i)*zuv           ! Function CA

        SELECT CASE ( icant )
          CASE (1)   ! Louis-transfer-scheme: additional laminar canopy resistance
            zustar = zuv*SQRT(tcm(i))
            zrla   = 1.0_wp/MAX(cdash*SQRT(zustar),eps_div)
          CASE (2)   ! Raschendorfer transfer scheme: laminar canopy resistance already considered
            zrla   = 0.0_wp
          CASE DEFAULT ! additional laminar canopy resistance for other schemes
            zrla   = 1.0_wp/MAX(tch(i)*zuv,eps_div)
        END SELECT

        ! to compute CV, first the stomatal resistance has to be determined
        ! this requires the determination of the F-functions:
        ! Radiation function
        IF (itype_trvg == 3) THEN
          ! modification depending on accumulated plant evaporation in order to reduce evaporation in the evening
          zxx      = 0.75_wp*MAX(eps_soil,ABS(plevap(i)))/MAX(0.2_wp,plcov(i))
          zzz      = MIN(3._wp,MAX(1._wp,zxx))
          ! stronger limitation for non-forest vegetation classes
          IF (z0(i) <= 0.4_wp) zzz = MIN(2._wp, zzz)
        ELSE
          zzz = 1._wp
          zxx = 1._wp
        ENDIF
        zpar     = pabs(i)/(cparcrit*zzz)   ! normalized PAR
        zf_rad(i)= MAX(0.0_wp,MIN(1.0_wp,zpar))
        ztlpmwp(i) = (zfcap(i,1) - zpwp(i,1))*(0.81_wp +       &
           0.121_wp*ATAN(-86400._wp*zep_s(i) - 4.75_wp))

        ! Soil water function
        IF (itype_root == 2) THEN
          zf_wat = MAX(0.0_wp,MIN(1.0_wp,(zwrootdz_int(i) - zpwp(i,1))/ztlpmwp(i)))
        ELSE
          zf_wat = MAX(0.0_wp,MIN(1.0_wp,(zwroot(i) - zpwp(i,1))/ztlpmwp(i)))
        ENDIF

        ! Temperature function
        ! T at lowest model level used (approximation of leaf height)
        zf_tem     = MAX(0.0_wp,MIN(1.0_wp,4.0_wp*     &
                     (t(i)-t0_melt)*(ctend-t(i))/(ctend-t0_melt)**2))

        ! Saturation deficit function (not used, see below)
        ! z2iw       = ztsk_pm(i)*b2w + (1.0_wp - ztsk_pm(i))*b2i
        ! z4iw       = ztsk_pm(i)*b4w + (1.0_wp - ztsk_pm(i))*b4i
        ! zepsat     = zsf_psat_iw(t(i),z2iw,z4iw)
        ! zepke      = qv(i) * ps(i) / (rdv + o_m_rdv*qv(i))
        ! zf_sat     = MAX(0.0_wp,MIN(1.0_wp,1.0_wp - (zepsat - zepke)/csatdef))

        ! zf_sat paralysed:
        zf_sat     = 1.0_wp

        IF (lstomata) THEN
          IF (itype_trvg == 3) THEN
            ! Modification of rsmin depending on accumulated plant evaporation; the z0 dependency
            ! is used to get a stronger effect for trees than for low vegetation

            zzz = MAX(0.5_wp+MIN(0.5_wp,1.0_wp-laifac(i)), EXP(SQRT(z0(i))*LOG(zxx)) )

          ELSE
            zzz = 1.0_wp
          ENDIF
          zedrstom   = 1.0_wp/crsmax + (1.0_wp/MAX(40._wp,zzz*rsmin2d(i)) - &
                       1.0_wp/crsmax)*zf_rad(i)*zf_wat*zf_tem*zf_sat
        ELSE
          zedrstom   = 1.0_wp/crsmax + (1.0_wp/crsmin -                     &
                       1.0_wp/crsmax)*zf_rad(i)*zf_wat*zf_tem*zf_sat
        END IF

        zrstom     = 1.0_wp/zedrstom              ! stomatal resistance
        rstom(i) = zrstom
        zrveg      = zrla + zrstom

        ! Transpiration rate of dry leaves:
        ztraleav(i)=zep_s(i)*tai(i)/(sai(i)+zrveg*zcatm)
      END IF  ! upwards directed potential evaporation only
    END DO

    ! Consideration of water and snow coverage, distribution to the different
    ! soil layers

    IF (itype_interception == 1) THEN
      !$ACC LOOP SEQ
      DO kso       = 1,ke_soil
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ztrabpf, ztrfr, zrootfc)
        DO i = ivstart, ivend
          IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_wp) THEN
            ! upwards potential evaporation

            ztrabpf  = ztraleav(i)*                  & ! plant covered part
            !!!           (1.0_wp - zf_wi(i))*       & ! not water covered
                          (1.0_wp - zf_snow(i))        ! not snow covered

            ! for root distribution
            IF (itype_root == 2) THEN
              ztrfr    = zwrootdz(i,kso)/(zrootdz_int(i)*zwrootdz_int(i))
              ztrang(i,kso) = ztrabpf*ztrfr
            ELSE
              zrootfc = zropartw(i,kso)/(zwroot(i) + eps_div)
              ztrang(i,kso) = ztrabpf*zrootfc/MAX(eps_div,zbwt(i))
            ENDIF

            ! Limit evaporation such that the soil water content does not fall beyond the wilting point
            IF(zw_fr(i,kso)+ztrang(i,kso)*zdtdrhw/zdzhs(kso) < zpwp(i,kso)) &
              ztrang(i,kso) = MIN(0._wp,(zpwp(i,kso)-zw_fr(i,kso))*zdzhs(kso)/zdtdrhw)

            IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) THEN    ! AYu mire block
              ztrang(i,kso) = 0.0_wp
            ENDIF

            lhfl_pl(i,kso)= lh_v * ztrang(i,kso)
            ztrangs(i)    = ztrangs(i) + ztrang(i,kso)
          END IF  ! upwards directed potential evaporation only
        END DO
      END DO          ! loop over soil layers
    ELSEIF (itype_interception == 2) THEN
      !$ACC LOOP SEQ
      DO kso = 1,ke_soil
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ztrabpf, ztrfr, zrootfc)
        DO i = ivstart, ivend
          IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_wp) THEN
            ! upwards potential evaporation

            ztrabpf  = ztraleav(i)*                   & ! plant covered part
                          (1.0_wp - zfd_wi(i))*       & ! not water covered
                          (1.0_wp - zf_snow(i))         ! not snow covered

            ! for root distribution
            IF (itype_root == 2) THEN
              ztrfr    = zwrootdz(i,kso)/(zrootdz_int(i)*zwrootdz_int(i))
              ztrang(i,kso) = ztrabpf*ztrfr
            ELSE
              zrootfc = zropartw(i,kso)/(zwroot(i) + eps_div)
              ztrang(i,kso) = ztrabpf*zrootfc/MAX(eps_div,zbwt(i))
              IF(zw_fr(i,kso)+ztrang(i,kso)*zdtdrhw/zdzhs(kso) &
                                .LT.zpwp(i,kso)) ztrang(i,kso) = 0._wp
            ENDIF
            lhfl_pl(i,kso)= lh_v * ztrang(i,kso)
            ztrangs(i)    = ztrangs(i) + ztrang(i,kso)
          END IF  ! upwards directed potential evaporation only
        END DO
      END DO          ! loop over soil layers
    END IF
  END IF BATS

  !----------------------------------------------------------------------------
  ! Section I.4.4: total evapotranspiration and
  !              associated fictitious soil humidity qv_s
  !----------------------------------------------------------------------------

#ifdef __SX__
  lhfl_pl_int(:) = 0._wp
  zfac(:)        = 1._wp
  DO kso = 1,ke_soil_hy
    DO i = ivstart, ivend
      lhfl_pl_int(i) = lhfl_pl_int(i) + lhfl_pl(i,kso)
    ENDDO
  ENDDO
#endif

  ! Ensure that the sum of the evaporation terms does not exceed the potential evaporation
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ze_sum, zxx, zzz)
  DO i = ivstart, ivend
    ze_sum = zdwsndt(i) + zdwidt(i) + zesoil(i) + ztrangs(i)
    zxx    = zf_snow(i)*zep_snow(i) + (1._wp-zf_snow(i))*zep_s(i) ! snow-weighted potential evaporation
    IF (zxx < 0._wp .AND. ze_sum < zxx) THEN
      zzz = zxx/ze_sum
      zdwsndt(i) = zdwsndt(i)*zzz
      zdwidt(i)  = zdwidt(i) *zzz
      zesoil(i)  = zesoil(i) *zzz
      ztrangs(i) = ztrangs(i)*zzz
      lhfl_bs(i) = lhfl_bs(i)*zzz
#ifdef __SX__
      zfac(i)    = zzz
#else
      ztrang(i,:)  = ztrang(i,:)*zzz
      lhfl_pl(i,:) = lhfl_pl(i,:)*zzz
#endif
    ENDIF
    IF (itype_trvg == 3) THEN
      ! accumulated plant evaporation since sunrise; an offset is subtracted to parameterize the
      ! amount of water that can be continuously supplied through the stem, and an increased recovery
      ! rate is assumed at night
      zzz = MAX(1._wp, 0.1_wp*(50._wp - pabs(i)))
      plevap(i) = MAX(-6._wp, MIN(0._wp, plevap(i) + zzz*dt/lh_v *           &
#ifdef __SX__
                 (lhfl_pl_int(i)+MAX(0.2_wp,plcov(i))*75._wp) ))
#else
                 (SUM(lhfl_pl(i,1:ke_soil_hy))+MAX(0.2_wp,plcov(i))*75._wp) ))
#endif
    ENDIF
    ! Negative values of tsnred indicate that snow is present on the corresponding snow tile
    ! and that the snow-free tile has been artificially generated by the melting-rate parameterization
    ! in this case, bare soil evaporation and, in the case of a long-lasting snow cover, plant evaporation, are turned off.
    IF ( (tsnred(i) < 0.0_wp) .AND. zep_s(i) < 0.0_wp) THEN
      zzz = MAX(0.0_wp,1.0_wp-ABS(tsnred(i)))
      zxx = MIN(1.0_wp,MAX(0.0_wp,2.0_wp-ABS(tsnred(i))))
      zesoil(i)    = zzz*zesoil(i)
      lhfl_bs(i)   = zzz*lhfl_bs(i)
      ztrangs(i)   = zxx*ztrangs(i)
#ifdef __SX__
      zfac(i)      = zxx*zfac(i)
#else
      ztrang(i,:)  = zxx*ztrang(i,:)
      lhfl_pl(i,:) = zxx*lhfl_pl(i,:)
#endif
    ENDIF
  ENDDO

#ifdef __SX__
!$NEC outerloop_unroll(7)
  DO kso = 1,ke_soil
    DO i = ivstart, ivend
      ztrang(i,kso)  = ztrang(i,kso)*zfac(i)
      lhfl_pl(i,kso) = lhfl_pl(i,kso)*zfac(i)
    ENDDO
  ENDDO
#endif

  IF (itype_interception == 1) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ze_sum)
    DO i = ivstart, ivend

      ze_sum = zdwsndt(i  )  & ! evaporation of snow
             + zdwidt (i  )  & ! evaporation from interception store
             + zesoil (i  )  & ! evaporation from bare soil
             + ztrangs(i  )  & ! transpiration from all soil layers
             + zrr    (i  )  & ! formation of dew
             + zrs    (i  )    ! formation of rime

      qv_s(i) = qv (i) - ze_sum /(zrhoch(i) + eps_div)
    END DO
  ELSE          IF (itype_interception == 2) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ze_sum)
    DO i = ivstart, ivend
      ze_sum = zesn   (i) &
             + zepd   (i) &
             + zewi   (i) &
             + zesoil (i) &
             + ztrangs(i) &
             + zdrr   (i) &
             + zrrs   (i)

      qv_s(i) = qv (i) - ze_sum /(zrhoch(i) + eps_div)
    END DO
  END IF
  !$ACC END PARALLEL


!------------------------------------------------------------------------------
! End of former module procedure terra1
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Former SUBROUTINE terra2 (yerror, ierror)
!------------------------------------------------------------------------------

!   In the prognostic part II the equation of heat conduction and water
!   transport is solved for a multi-layer soil using the same vertical grid
!   Freezing/melting of soil water/ice is accounted for (optionally). A
!   simple one-layer snow model provides the snow surface temperature and
!   the snow water equivalent.

!------------------------------------------------------------------------------
! Section II.1: Initializations
!------------------------------------------------------------------------------

  ! Computation of derived constants
  z1d2dt = 1.0_wp/zdt      ! 1./2*timestep

  ! Number of soil layers contributing to surface run-off
  ! (excess water removed from the uppermost layer constitutes the surface run-off)
  msr_off  = 1

  zdtdrhw = zdt/rho_w   ! timestep/density of liquid water

  ! time constant for infiltration of water from interception store
  ! must not be less than 2*time step
  ctau_i  = MAX( ctau_i, zdt )

  ! Utility variable to avoid IF-constructs
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO kso     = 1,ke_soil
    IF (kso == 1) THEN 
      znlgw1f(1)    = 1.0_wp
    ELSE
      znlgw1f(kso)  = 0.0_wp
    END IF
  END DO
  !$ACC END PARALLEL

  ! Initialisations
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc) &
  !$ACC   PRIVATE(zro_sfak, zro_gfak, zfmb_fak, zcou_roffg)
  IF (lmulti_snow) THEN
    !$ACC LOOP SEQ
    DO ksn = 0, ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        zdtsnowdt_mult(i,ksn)  = 0.0_wp
        zdtsdt(i)     = 0.0_wp
      END DO
    END DO
  ELSE
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      zdtsnowdt(i)  = 0.0_wp
      zdtsdt(i)     = 0.0_wp
    END DO
  END IF

!------------------------------------------------------------------------------
! Section II.2: Prepare basic surface properties and create some local
!               arrays of surface related quantities (for land-points only)
!               Initialise some fields
!------------------------------------------------------------------------------

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp)
  DO i = ivstart, ivend
!    mstyp      = soiltyp_subs(i)
!    m_styp (i) = mstyp
!    zsandf (i) = csandf(mstyp)
!    zclayf (i) = cclayf(mstyp)
     zgsb   (i) = 0.0_wp
     ztrangs(i) = 0.0_wp
  END DO


  !$ACC LOOP SEQ
  DO   kso = 1,ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      ziw_fr(i,kso) = w_so_ice_now(i,kso)/zdzhs(kso)   ! ice frac.
      zlw_fr(i,kso) = zw_fr(i,kso) - ziw_fr(i,kso)  ! liquid water frac.
      zroc(i,kso)   = zrocg(i,kso) + rho_w*zlw_fr(i,kso)*chc_w +          &
                                     rho_w*ziw_fr(i,kso)*chc_i
      IF (kso<=ke_soil) THEN
        ztrangs(i) = ztrangs(i) + ztrang(i,kso)
        zdwgdt(i,kso) = 0.0_wp
      END IF
      zflmg(i,kso)  = 0.0_wp
      zrunoff_grav(i,kso)  = 0.0_wp
    END DO
  END DO      !soil layers

  IF (lterra_urb) THEN
    ! HW: modification of soil heat capacity in urban areas:
    !     interpolation between buildings and non-buildings, 
    !     modification decreases with respect to the natural soil below. 
    !     Below urb_h_bld, everything is natural soil.

    !$ACC LOOP SEQ
    DO kso = 1, ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zalpha_uf)
      DO i = ivstart, ivend

        zalpha_uf   = MAX(0.0_wp, MIN(zmls(kso)/urb_h_bld(i), 1.0_wp))

        zroc(i,kso) =           urb_isa(i)  * ( (1.0_wp - zalpha_uf) * urb_hcap(i)*urb_ai(i)   &
                                                +         zalpha_uf  * zroc(i,kso)           ) &
                    + (1.0_wp - urb_isa(i)) * zroc(i,kso)

      ENDDO
    ENDDO
  END IF

!------------------------------------------------------------------------------
! Section II.3: Estimate thermal surface fluxes
!------------------------------------------------------------------------------

  IF (itype_interception == 1) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zrime, zfr_ice_free, zalf, zuv, zzz, zinf) &
    !$ACC   PRIVATE(zdwidtt, zdwieps, zro_wi, zro_inf, zdwsndtt, zwsnstr) &
    !$ACC   PRIVATE(zdwseps)
    DO i = ivstart, ivend
      ! store forcing terms due to evapotranspiration, formation of dew
      ! and rime for later use
      zverbo(i) = zdwidt(i) + zesoil(i) + ztrangs(i) +              &
                      (1._wp-zf_snow(i))*(zrr(i) + zrs(i))
      zversn(i) = zdwsndt(i) + zrs(i)                                 &
                                + zsf_heav (zwsnow(i) - eps_soil) * zrr(i)

      ! add grid scale and convective precipitation (and graupel, if present)
      ! to dew and rime
      zrr(i) = zrr(i) + prr_con(i) + prr_gsp(i)
      zrime  = zrs(i)
      zrs(i) = zrs(i) + prs_con(i) + prs_gsp(i) + pri_gsp(i)

      IF ( nclass_gscp >= 2000 ) THEN
        ! only possible when running 2-moment microphysics
#ifdef TWOMOM_SB
        zrs(i) = zrs(i) + prg_gsp(i) + prh_gsp(i)
#endif
      ELSEIF ( nclass_gscp >= 6 ) THEN
        zrs(i) = zrs(i) + prg_gsp(i)
      ENDIF

      ! add possible deposition of water droplets or ice particles:
      IF (ldepo_qw) THEN !cloud-particle deposition related to vertical diffsuion with zero-concentration condition
                      zrr(i) = zrr(i) + zrhoch(i)*qc(i) !dew -fall + liquid precipitation including deposition
        IF (ldiff_qi) zrs(i) = zrs(i) + zrhoch(i)*qi(i) !rime-fall + frozn. precipitation including deposition
      END IF  

      ! Decide whether riming is added to interception store or to snow cover
      IF (zrs(i) >= 1.05_wp*zrime .OR. zf_snow(i) >= 0.9_wp) zrime = 0._wp
      zrs(i) = zrs(i) - zrime

      !
      ! Infiltration and surface runoff
      !

      ! subtract evaporation from interception store to avoid negative
      ! values due to sum of evaporation+infiltration
      zwinstr(i) = zwin(i) + zdwidt(i)*zdtdrhw
      zwinstr(i) = MAX(0.0_wp,zwinstr(i))

      zwimax(i)   = cwimax_ml*(1.0_wp+ztfunc(i))*MAX(ztfunc(i), eps_soil, MAX(2.5_wp*plcov(i),tai(i)))

      ! TERRA_URB: Puddles on the impervious surface area
      IF (lterra_urb .AND. (itype_eisa == 3)) THEN
        zwimax(i) = urb_isa(i) * cwisamax + (1.0_wp - urb_isa(i)) * zwimax(i)
      END IF

      zalf   = SQRT(MAX(0.0_wp,1.0_wp - zwinstr(i)/zwimax(i)))

      ! water supply from interception store (if Ts above freezing)
      zuv    = SQRT ( u_10m(i)**2 + v_10m(i)**2 )
      zzz    = MAX(0.1_wp, 0.4_wp - 0.05_wp*zuv)
      zinf   = MAX(0._wp,zwinstr(i)-zzz*zwimax(i))*rho_w/ctau_i*        &
               (1._wp+0.75_wp*MAX(0._wp,zuv-1._wp))*(1._wp-ztfunc(i))**2

      ! possible contribution of rain to infiltration
      !    IF (zrr(i)-eps_soil > 0.0_wp) THEN
      !      zalf = MAX( zalf,                                                   &
      !             (zrhwddt*MAX(0.0_wp, zwimax(i)-zwinstr(i)) + zinf)/zrr(i) )
      !      zalf = MAX( 0.01_wp, MIN(1.0_wp, zalf) )

      ! if rain falls onto snow, all rain is considered for infiltration
      ! as no liquid water store is considered in the snowpack
      IF (zwsnow(i) > 0.0_wp) zalf = 0.0_wp
        ! Increase infiltration and reduce surface runoff (bugfix)
        ! this deactivates filling the interception store!!

      ! rain freezes on the snow surface
      IF (lmulti_snow .AND. zwsnow(i) > 0.0_wp) zalf = 1.0_wp

      ! interception store; convective precip is taken into account with a 
      ! fractional area passed from the convection scheme
      zdwidtt  = zalf*(zrr(i)+(conv_frac(i)-1._wp)*prr_con(i)) + zrime + zdwidt(i)-zinf
      zwinstr(i)  = zwin(i) + zdwidtt*zdtdrhw
      zwinstr(i)  = MAX(0.0_wp, zwinstr(i)) !avoid negative values (security)
      zdwieps  = 0.0_wp
      IF (zwinstr(i) > 0.0_wp .AND. zwinstr(i) < 1.0E-4_wp*eps_soil) THEN
        zdwieps    = zwinstr(i)*zrhwddt
        runoff_s(i)= runoff_s(i) + zdwieps*zroffdt
        zdwidtt    = zdwidtt - zdwieps
        zwinstr(i)    = 0.0_wp
      END IF

      ! add excess over zwimax(i) to infiltration
      zro_wi       = zrhwddt*MAX( 0.0_wp, zwinstr(i)-zwimax(i) )
      zdwidtt      = zdwidtt - zro_wi
      zdwidt(i)    = zdwidtt
      zinf         = zinf + zro_wi
      IF (zts(i) <= t0_melt) THEN
        ! add excess rime to snow
        zrs(i) = zrs(i) + zinf
        zinf = 0.0_wp
      ENDIF

      ! add rain contribution to water supply for infiltration
      ! surface runoff is evaluated after the calculation of infiltration
      zvers(i) = zinf + (1._wp - zalf)*zrr(i) + (1._wp-conv_frac(i))*zalf*prr_con(i)


      ! Infiltration
      !
      ! Avoid infiltration for rock, ice and snow-covered surfaces
      zinfil(i) = zvers(i)*zrock(i)*(1.0_wp - zf_snow(i))

      ! Avoid infiltration for urban impervious surface area fraction
      IF (lterra_urb .AND. ((itype_eisa == 2) .OR. (itype_eisa == 3))) THEN
        zinfil(i) = zinfil(i) * (1.0_wp - urb_isa(i))
      END IF


      ! Add difference to surface runoff
      zro_inf = zvers(i) - zinfil(i)

      ! Surface runoff, including the effect of the impervious surface area, if present
      !
      runoff_s(i) = runoff_s(i) + zro_inf*zroffdt


      ! change of snow water and interception water store
      ! (negligible residuals are added to the run-off)

      ! snow store
      zdwsndtt = zrs(i) + zdwsndt(i)
      zwsnstr  = zwsnow(i) + zdwsndtt*zdtdrhw
      zwsnstr  = MAX(0.0_wp, zwsnstr) ! avoid negative values (security)
      zdwseps  = 0.0_wp
      IF (zwsnstr > 0.0_wp .AND. zwsnstr < eps_soil) THEN
        ! shift marginal snow amounts to interception storage
!       IF (ztsnow_pm(i) > 0.0_wp) THEN
          zwinstr(i) = zwinstr(i) + zwsnstr
          zdwseps    = zwsnstr*zrhwddt
!         runoff_s(i) = runoff_s(i) + zdwseps*zroffdt   ! previous implementation
          zdwsndtt   = zdwsndtt - zdwseps
          zdwidt(i)  = zdwidt(i) + zdwseps
!       END IF
      END IF
      zdwsndt(i) = zdwsndtt
    END DO

  ELSEIF (itype_interception == 2) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zfr_ice_free, zfcorr_wi, zdwidtt, zro_wi) &
    !$ACC   PRIVATE(zuv, zdwisndtt, zro_inf, zdwpdtt, zdwsndtt, zwsnstr, zdwseps)
    DO i = ivstart, ivend
      ! Compute forcing terms after (???) interception, infiltration calculation --> below
      ! store forcing terms due to evapotranspiration, formation of dew
      ! and rime for later use
      zverbo(i) = zewi(i) + zesoil(i) + zepd(i) + ztrangs(i) +              &
                      (1._wp-zf_snow(i))*(zdrr(i) + zrrs(i))

      zversn(i) = zdwsndt(i) + zrrs(i) + zsf_heav (zwsnow(i) - eps_soil) * zdrr(i)

      ! ice free fraction of first soil layer scaled by pore volume
      ! is used as reduction factor for infiltration rate (previously zts_pm switch bobo)
      zfr_ice_free     = 1._wp-ziw_fr(i,1)/zporv(i,1)

      ! add grid scale and convective precipitation (and graupel, if present) to dew and rime
      zrr(i) = zdrr(i) + prr_con(i) + prr_gsp(i)

      IF ( nclass_gscp >= 2000 ) THEN
        ! only possible when running 2-moment microphysics
#ifdef TWOMOM_SB
        zrs(i) = zrrs(i) + prs_con(i) + prs_gsp(i) + pri_gsp(i) + prg_gsp(i) + prh_gsp(i)
#endif
      ELSEIF ( nclass_gscp >= 6 ) THEN
        zrs(i) = zrrs(i) + prs_con(i) + prs_gsp(i) + pri_gsp(i) + prg_gsp(i)
      ELSE
        zrs(i) = zrrs(i) + prs_con(i) + prs_gsp(i) + pri_gsp(i)
      ENDIF

      ! add possible deposition of water droplets or ice particles:
      IF (ldepo_qw) THEN !cloud-particle deposition related to vertical diffsuion with zero-concentration condition
                      zrr(i) = zrr(i) + zrhoch(i)*qc(i) !dew -fall + liquid precipitation including deposition
        IF (ldiff_qi) zrs(i) = zrs(i) + zrhoch(i)*qi(i) !rime-fall + frozn. precipitation including deposition
      END IF  

      ! Preliminary interception store budget
      ! According to Wang et al. (Evaluation of canopy interception schemes in land 
      !                           surface models, J. of Hydrology,2007)
      ! the water balance equation for the interception store is   
      !       zdwidtt  = Ic + Ew - Dr, 
      ! where Ic is the canopy interception rate, 
      !       Dr the canopy drip rate, and 
      !       Ew is the evaporation rate from wet foliage.
      ! This equation is independent of the approach used to model the canopy 
      ! hydrological processes.
      
      ! Comparable to the community land model (CLM), the precipitation intercepted 
      ! by vegetation canopy Ic is considered as an exponential function of canopy 
      ! density using a canopy cover fraction
      !      zf_wi  (j1,j2) = 1._wp - exp(-0.5_wp * plai(j1,j2)).  
      ! Therefore, similar to plai, the  canopy cover fraction shows an annual cycle.
      ! The canopy capacity zwimax(i) is considered linearly related to leaf area index
      ! (van Dijk and Bruijnzeel, 2001; Wang et al., 2007)
      ! The canopy dripping occurs only when canopy water storage winstr exceeds 
      ! the water holding capacity zwimax(i).
      ! This interception runoff is part of the soil infiltration.
      ! The melting of snow on the leafs, if T > T_melt is also considered, but from a 
      ! frozen interception store only evaporation is allowed.


      ! Interception of rain water
      IF (ztsnow(i).gt.t0_melt) THEN

        ! snow fall on leafs/soil with T>T_melt, snow water content increases
        ! interception store water content
        zfcorr_wi=1._wp ! Start value

        zdwidtt  = zf_wi(i)*zrr(i) + zf_wi(i)*zrs(i) + zewi(i)
        zwinstr(i)  = zwin (i) + zdtdrhw * zdwidtt
        zewi(i) = zewi(i) - MIN(0._wp,zwinstr(i) * zrhwddt) ! Partial evaporation


        ! Final calculation of the interception store budget
        zdwidtt  = zf_wi(i)*zrr(i) + zf_wi(i)*zrs(i) + zewi(i)
        zwinstr(i)  = zwin (i) + zdtdrhw * zdwidtt

      ELSE   IF (ztsnow(i).le.t0_melt) THEN!  T =< t0_melt

        zdwidtt = zewi(i) ! only evaporation from frozen interception store allowed
        zwinstr(i)  =  zwin (i) + zdtdrhw * zdwidtt
        zewi(i) = zewi(i) - MIN(0._wp,zwinstr(i) * zrhwddt) ! Partial evaporation


        ! Final calculation of the interception store budget
        zdwidtt  =  zewi(i) ! only evaporation from frozen interception store allowed
        zwinstr(i)  = zwin (i) + zdtdrhw * zdwidtt

      END IF ! Temperature

      zro_wi        = MAX( 0.0_wp, zwinstr(i)-zwimax(i) )   ! excess over zwimax -> runoff and infiltration
      zwinstr(i)  = MAX(0.0_wp, zwinstr(i))                 ! avoid negative values (security)
      zwinstr(i)  = MIN(zwinstr(i),zwimax(i))               ! correction of interc. store due to runoff
      zdwidt(i)   = zdwidtt - zro_wi

      zuv        = SQRT ( u(i)**2 + v(i)**2 )

      ! SNOW Interception Model (Roesch et al., 2001)
      ! Need forest_fraction
      !   zdwisndtt=(for_e(i)+for_d(i))*zrs(i) + (for_e(i)+for_d(i))*zesn(i) - zwisn(i)* &
      !   ((t(i,j,ke,nx)-270.15)/1.87E5_wp + SQRT(u(i,j,ke,nx)*u(i,j,ke,nx)+v(i,j,ke,nx)*v(i,j,ke,nx))/1.56E5_wp)
      ! Just for testing, using plcov
      zdwisndtt=plcov(i)*zrs(i) + plcov(i)*zesn(i) - zwisn(i)* &
                ((t(i)-270.15)/1.87E5_wp + zuv/1.56E5_wp)
      zwisnstr(i)  = zwisn (i) + zdtdrhw * zdwisndtt
      zwisnstr(i)  = MAX(0.0_wp, zwisnstr(i)) !avoid negative values (security)

      ! Calculation of infiltration
      ! add rain and interception excess to water supply for infiltration
      zvers(i) =  zrhwddt*zro_wi + (1._wp - zf_wi(i))*zrr(i) ! Test  Excess over zwimax will infiltrated

      ! maximum infiltration rate of the soil (rock/ice/water-exclusion
          !      zinfmx(i) = zrock(i)*zfr_ice_free *csvoro &
          !       *( cik1*MAX(0.5_wp,plcov(i))*MAX(0.0_wp,           &
          !          zporv(i,1)-zw_fr(i,1))/zporv(i,1) + zik2(i) )

      zinfmx(i) = zrock(i)*zfr_ice_free*csvoro*zkw0(i)*rho_w

      ! to avoid pore volume water excess of the uppermost layer by infiltration
      zinfmx(i) = MIN(zinfmx(i), (zporv(i,1) - zw_fr(i,1))*zdzhs(1)*zrhwddt)

      ! final infiltration rate limited by maximum value
      zinfil(i) = MIN(zinfmx(i),zvers(i)+ zwpn(i)*zrhwddt ) !GME

      ! surface run-off (residual of potential minus actual infiltration)
      zro_inf     = MAX(0._wp, zvers(i) - zinfil(i))

      ! Pond interception
      zdwpdtt  = zvers(i) - zinfil(i)  + zepd(i)
      zwpnstr(i)  = zwpn   (i) + zdtdrhw * zdwpdtt
      IF (zvers(i) - zinfil(i).ge.0._wp) then
        zepd(i) = zepd(i) - MIN(0._wp,zwpnstr(i) * zrhwddt) ! Partial evaporation
      ELSE
        zepd(i) = zepd(i) - MIN(0._wp,zwpn(i) * zrhwddt) ! Partial evaporation
      END IF

      ! Final calculation
      zdwpdtt  = zvers(i) - zinfil(i)  + zepd(i)
      zwpnstr(i)  = MAX(0._wp,zwpn (i) + zdtdrhw * zdwpdtt )
      zro_inf = MAX(0._wp, (zwpnstr(i)  - zpercmax)*zrhwddt )   ! Pond overflow?
      zwpnstr(i)  = MIN(zwpnstr(i),zpercmax) ! Correction of interc. store due to runoff

      ! Calculation of surface runoff
      ! change of snow water and interception water store
      ! (negligible residuals are added to the run-off)

      ! snow store
      zdwsndtt = zrs(i) + zdwsndt(i)
      zwsnstr  = w_s_now(i) + zdwsndtt*zdtdrhw
      zwsnstr  = MAX(0.0_wp, zwsnstr) ! avoid negative values (security)
      zdwseps  = 0.0_wp
      IF (zwsnstr > 0.0_wp .AND. zwsnstr < eps_soil) THEN
        zdwseps    = zwsnstr*zrhwddt
        runoff_s(i) = runoff_s(i) + zdwseps*zroffdt
        zdwsndtt   = zdwsndtt - zdwseps
      END IF
      zdwsndt(i) = zdwsndtt
      runoff_s(i)= runoff_s(i) + zro_wi*zroffdt
    END DO

  END IF

!------------------------------------------------------------------------------
! Section II.4: Soil water transport and runoff from soil layers
!------------------------------------------------------------------------------

! uppermost layer, kso = 1


  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zice_fr_kso, zice_fr_ksop1, zlw_fr_kso) &
  !$ACC   PRIVATE(zlw_fr_ksop1, zfr_ice, zredp05, zlw_fr_ksop05) &
  !$ACC   PRIVATE(zdlw_fr_ksop05, zklw_fr_ksop05, z1dgam1, zgam2p05)
  DO i = ivstart, ivend
    IF (soiltyp_subs(i)  >= 3) THEN
      ! sedimentation and capillary transport in soil
      ! Note: The fractional liquid water content (concentration)  of each layer
      !       is normalized by the ice free fraction of each layer in order to
      !       obtain a representative concentration of liquid water in the
      !       'active' part of each soil layer
      !       Hydraulic diffusivity and conductivity coefficients are multiplied
      !       by a reduction factor depending on the maximum ice fraction of the
      !       adjacent layers in order to avoid the transport of liquid water
      !       in to the frozen part of the adjacent layer
      !
      ! GZ, 2017-11-07:
      ! I wonder if this is sufficient to describe the reduction of water transport in frozen soils.
      ! Certainly, the hydraulic conductivity does not drop abruptly to zero once first ice crystals form in the soil,
      ! but if rain falls after a longer-term frost period without snow cover, ponds form even if the soil surface
      ! has already thawed. Should the hydraulic conductivity include a term like ...*(1-tuning_factor*ziw_fr/zporv)?
      !
      zice_fr_kso   = ziw_fr(i,1)
      zice_fr_ksop1 = ziw_fr(i,2)
      zlw_fr_kso    = zlw_fr(i,1)
      zlw_fr_ksop1  = zlw_fr(i,2)
  
      ! compute reduction factor for transport coefficients
      zfr_ice  = max (zice_fr_kso,zice_fr_ksop1)
      zredp05  = 1._wp-zfr_ice/MAX(zlw_fr_kso+zice_fr_kso,zlw_fr_ksop1+zice_fr_ksop1)
  
      ! interpolated scaled liquid water fraction at layer interface
      zlw_fr_ksop05  = 0.5_wp*(zdzhs(2)*zlw_fr_kso+zdzhs(1)*zlw_fr_ksop1) / zdzms(2)
      !!$          zdlw_fr_ksop05 = zredp05*zdw(i,1)*EXP(zdw1(i,1)*                        &
      !!$                           (zporv(i,1)-zlw_fr_ksop05)/(zporv(i,1)-zadp(i,1)) )
      !!$          zklw_fr_ksop05 = zredp05*zkw(i,1)*EXP(zkw1(i,1)*                        &
      !!$                           (zporv(i,1)-zlw_fr_ksop05)/(zporv(i,1)-zadp(i,1)) )
  
      zdlw_fr_ksop05= zredp05 * watrdiff_RT(zdw(i,1),zlw_fr_ksop05, zdw1(i,1),zporv(i,1),zadp(i,1))
      zklw_fr_ksop05= zredp05 * watrcon_RT (zkw(i,1),zlw_fr_ksop05, zkw1(i,1),zporv(i,1),zadp(i,1))
  
      ! coefficients for implicit flux computation
      z1dgam1     = zdt/zdzhs(1)
      zgam2p05    = zdlw_fr_ksop05/zdzms(2)
      zaga(i,1) = 0._wp
      zagb(i,1) = 1._wp+zalfa*zgam2p05*z1dgam1
      zagc(i,1) = -zalfa * zgam2p05*z1dgam1
      zagd(i,1) = zlw_fr(i,1) + zinfil(i)*z1dgam1/rho_w  &
                         -zklw_fr_ksop05*z1dgam1                     &
                         +(1._wp - zalfa)* zgam2p05*z1dgam1*(zlw_fr_ksop1 - zlw_fr_kso)  &
                         +                 zgam2p05*z1dgam1*(zice_fr_ksop1-zice_fr_kso)
  
      ! explicit part of soil surface water flux:
      zflmg (i,1) = - zinfil(i)! boundary value for soil water transport
    END IF
  END DO

! inner layers 2 <=kso<=ke_soil_hy-1

  !$ACC LOOP SEQ
!DIR$ ivdep, prefervector
  DO kso =2,ke_soil_hy-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zice_fr_ksom1, zice_fr_kso) &
    !$ACC   PRIVATE(zice_fr_ksop1, zlw_fr_ksom1, zlw_fr_kso, zlw_fr_ksop1) &
    !$ACC   PRIVATE(zlw_fr_ksom05, zlw_fr_ksop05, zfr_ice, zredm05, zredp05) &
    !$ACC   PRIVATE(zdlw_fr_ksom05, zdlw_fr_ksop05, zklw_fr_ksom05) &
    !$ACC   PRIVATE(zklw_fr_ksop05, z1dgam1, zgam2m05, zgam2p05)
    DO i = ivstart, ivend
      IF (soiltyp_subs(i) >= 3) THEN
        ! sedimentation and capillary transport in soil
        zice_fr_ksom1 = ziw_fr(i,kso-1)
        zice_fr_kso   = ziw_fr(i,kso  )
        zice_fr_ksop1 = ziw_fr(i,kso+1)
        zlw_fr_ksom1  = zlw_fr(i,kso-1)
        zlw_fr_kso    = zlw_fr(i,kso  )
        zlw_fr_ksop1  = zlw_fr(i,kso+1)
        ! interpolated scaled liquid water content at interface to layer
        ! above and below
        zlw_fr_ksom05 = 0.5_wp*(zdzhs(kso-1)*zlw_fr_kso+   &
                               zdzhs(kso)*zlw_fr_ksom1)/zdzms(kso)
        zlw_fr_ksop05 = 0.5_wp*(zdzhs(kso+1)*zlw_fr_kso+   &
                               zdzhs(kso)*zlw_fr_ksop1)/zdzms(kso+1)

        ! compute reduction factor for coefficients
        zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
        zredm05 = 1._wp-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)
        zfr_ice          = max (zice_fr_kso,zice_fr_ksop1)
        zredp05 = 1._wp-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksop1+zice_fr_ksop1)
             !!$    zdlw_fr_ksom05= zredm05*zdw(i,kso)*EXP( zdw1(i,kso)*   &
             !!$                       (zporv(i,kso)-zlw_fr_ksom05)/(zporv(i,kso)-zadp(i,kso)) )
             !!$    zdlw_fr_ksop05= zredp05*zdw(i,kso)*EXP( zdw1(i,kso)*   &
             !!$                       (zporv(i,kso)-zlw_fr_ksop05)/(zporv(i,kso)-zadp(i,kso)) )
             !!$    zklw_fr_ksom05= zredm05*zkw(i,kso)*EXP( zkw1(i,kso)*   &
             !!$                       (zporv(i,kso)-zlw_fr_ksom05)/(zporv(i,kso)-zadp(i,kso)) )
             !!$    zklw_fr_ksop05= zredp05*zkw(i,kso)*EXP( zkw1(i,kso)*   &
             !!$                       (zporv(i,kso)-zlw_fr_ksop05)/(zporv(i,kso)-zadp(i,kso)) )

        zdlw_fr_ksom05= zredm05*watrdiff_RT(zdw(i,kso),zlw_fr_ksom05,&
                              zdw1(i,kso),zporv(i,kso),zadp(i,kso))

        zdlw_fr_ksop05= zredp05*watrdiff_RT(zdw(i,kso),zlw_fr_ksop05,&
                              zdw1(i,kso),zporv(i,kso),zadp(i,kso))

        zklw_fr_ksom05= zredm05*watrcon_RT(zkw(i,kso-1),zlw_fr_ksom05,zkw1(i,kso),zporv(i,kso),zadp(i,kso))

        zklw_fr_ksop05= zredp05*watrcon_RT(zkw(i,kso),zlw_fr_ksop05,zkw1(i,kso),zporv(i,kso),zadp(i,kso))


        ! coefficients for implicit flux computation
        z1dgam1 = zdt/zdzhs(kso)
        zgam2m05  = zdlw_fr_ksom05/zdzms(kso)
        zgam2p05  = zdlw_fr_ksop05/zdzms(kso+1)
        zaga (i,kso) = -zalfa*zgam2m05*z1dgam1
        zagc (i,kso) = -zalfa*zgam2p05*z1dgam1
        zagb (i,kso) = 1._wp +zalfa*(zgam2m05+zgam2p05)*z1dgam1
        zagd (i,kso) = zlw_fr(i,kso)+                               &
                              z1dgam1*(-zklw_fr_ksop05+zklw_fr_ksom05)+ &
                              (1._wp-zalfa)*z1dgam1*                &
                              (zgam2p05*(zlw_fr_ksop1-zlw_fr_kso  )     &
                              -zgam2m05*(zlw_fr_kso  -zlw_fr_ksom1)   ) &
                             +z1dgam1*                                  &
                              (zgam2p05*(zice_fr_ksop1-zice_fr_kso  )   &
                              -zgam2m05*(zice_fr_kso-zice_fr_ksom1)   )

        !soil water flux, explicit part, for soil water flux investigations
        ! only)
        zflmg(i,kso) = rho_w &
                         *(zdlw_fr_ksom05*(zlw_fr_kso+zice_fr_kso-zlw_fr_ksom1-zice_fr_ksom1) &
                         / zdzms(kso) - zklw_fr_ksom05)

        IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) THEN
          IF(kso==ke_soil_hy_m-1) THEN
            zflmg(i,kso+1)=rho_w &
                    *(zdlw_fr_ksop05*(zlw_fr_ksop1+zice_fr_ksop1-zlw_fr_kso-zice_fr_kso) &
                    / zdzms(kso+1) - zklw_fr_ksop05)
          ENDIF
        ELSE
          IF(kso==ke_soil_hy-1) THEN
            zflmg(i,kso+1)=rho_w &
                    *(zdlw_fr_ksop05*(zlw_fr_ksop1+zice_fr_ksop1-zlw_fr_kso-zice_fr_kso) &
                    / zdzms(kso+1) - zklw_fr_ksop05)
          ENDIF
        ENDIF

      END IF
    END DO
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zice_fr_ksom1, zice_fr_kso, zlw_fr_ksom1) &
  !$ACC   PRIVATE(zlw_fr_kso, zlw_fr_ksom05, zfr_ice, zredm05, zdlw_fr_ksom05) &
  !$ACC   PRIVATE(z1dgam1, zgam2m05, zklw_fr_ksom05)
  DO i = ivstart, ivend
    IF (soiltyp_subs(i) >= 3) THEN
      IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) THEN
        ! lowest active hydrological layer ke_soil_hy-1
        zice_fr_ksom1 = ziw_fr(i,ke_soil_hy_m-1)
        zice_fr_kso   = ziw_fr(i,ke_soil_hy_m  ) 
        zlw_fr_ksom1  = zlw_fr(i,ke_soil_hy_m-1)
        zlw_fr_kso    = zlw_fr(i,ke_soil_hy_m  ) 
        zlw_fr_ksom05 = 0.5_wp*(zdzhs(ke_soil_hy_m-1)*zlw_fr_kso+ & 
                               zdzhs(ke_soil_hy_m)*zlw_fr_ksom1)/zdzms(ke_soil_hy_m)

        zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
        zredm05 = 1._wp-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)
          !!$   zdlw_fr_ksom05= zredm05*zdw(i,ke_soil_hy_m)*EXP( zdw1(i,ke_soil_hy_m)* &
          !!$                     (zporv(i,ke_soil_hy_m)-zlw_fr_ksom05)/(zporv(i,ke_soil_hy_m)-zadp(i,ke_soil_hy_m)) )

        zdlw_fr_ksom05= zredm05*watrdiff_RT(zdw(i,ke_soil_hy_m),zlw_fr_ksom05,&
                               zdw1(i,ke_soil_hy_m),zporv(i,ke_soil_hy_m),zadp(i,ke_soil_hy_m))


        z1dgam1 = zdt/zdzhs(ke_soil_hy_m)
        zgam2m05  = zdlw_fr_ksom05/zdzms(ke_soil_hy_m)
          !!$   zklw_fr_ksom05= zredm05*zkw(i,ke_soil_hy_m)*EXP( zkw1(i,ke_soil_hy_m)* &
          !!$               (zporv(i,ke_soil_hy_m)-zlw_fr_ksom05)/(zporv(i,ke_soil_hy_m)-zadp(i,ke_soil_hy_m)) )
          !!$

        zklw_fr_ksom05= zredm05*watrcon_RT(zkw(i,ke_soil_hy_m-1),zlw_fr_ksom05,&
                 zkw1(i,ke_soil_hy_m),zporv(i,ke_soil_hy_m),zadp(i,ke_soil_hy_m))

        zaga(i,ke_soil_hy_m) = -zalfa* zgam2m05*z1dgam1
        zagb(i,ke_soil_hy_m) = 1._wp+ zalfa*zgam2m05*z1dgam1
        zagc(i,ke_soil_hy_m) = 0.0_wp
        zagd(i,ke_soil_hy_m) = zlw_fr(i,ke_soil_hy_m)+z1dgam1*zklw_fr_ksom05         & 
                               +(1._wp-zalfa)*z1dgam1*                           & 
                                zgam2m05*(zlw_fr_ksom1  - zlw_fr_kso)            & 
                               +z1dgam1*                                         & 
                                zgam2m05*(zice_fr_ksom1-zice_fr_kso ) 
      ELSE
        ! lowest active hydrological layer ke_soil_hy-1
        zice_fr_ksom1 = ziw_fr(i,ke_soil_hy-1)
        zice_fr_kso   = ziw_fr(i,ke_soil_hy  )
        zlw_fr_ksom1  = zlw_fr(i,ke_soil_hy-1)
        zlw_fr_kso    = zlw_fr(i,ke_soil_hy  )
        zlw_fr_ksom05 = 0.5_wp*(zdzhs(ke_soil_hy-1)*zlw_fr_kso+ &
                              zdzhs(ke_soil_hy)*zlw_fr_ksom1)/zdzms(ke_soil_hy)

        zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
        zredm05 = 1._wp-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)
          !!$   zdlw_fr_ksom05= zredm05*zdw(i,ke_soil_hy)*EXP( zdw1(i,ke_soil_hy)* &
          !!$                     (zporv(i,ke_soil_hy)-zlw_fr_ksom05)/(zporv(i,ke_soil_hy)-zadp(i,ke_soil_hy)) )

        zdlw_fr_ksom05= zredm05*watrdiff_RT(zdw(i,ke_soil_hy),zlw_fr_ksom05,&
                              zdw1(i,ke_soil_hy),zporv(i,ke_soil_hy),zadp(i,ke_soil_hy))


        z1dgam1 = zdt/zdzhs(ke_soil_hy)
        zgam2m05  = zdlw_fr_ksom05/zdzms(ke_soil_hy)
          !!$   zklw_fr_ksom05= zredm05*zkw(i,ke_soil_hy)*EXP( zkw1(i,ke_soil_hy)* &
          !!$               (zporv(i,ke_soil_hy)-zlw_fr_ksom05)/(zporv(i,ke_soil_hy)-zadp(i,ke_soil_hy)) )
          !!$

        zklw_fr_ksom05= zredm05*watrcon_RT(zkw(i,ke_soil_hy-1),zlw_fr_ksom05,&
                zkw1(i,ke_soil_hy),zporv(i,ke_soil_hy),zadp(i,ke_soil_hy))

        zaga(i,ke_soil_hy) = -zalfa* zgam2m05*z1dgam1
        zagb(i,ke_soil_hy) = 1._wp+ zalfa*zgam2m05*z1dgam1
        zagc(i,ke_soil_hy) = 0.0_wp
        zagd(i,ke_soil_hy) = zlw_fr(i,ke_soil_hy)+z1dgam1*zklw_fr_ksom05         &
                               +(1._wp-zalfa)*z1dgam1*                           &
                                zgam2m05*(zlw_fr_ksom1  - zlw_fr_kso)            &
                               +z1dgam1*                                         &
                                zgam2m05*(zice_fr_ksom1-zice_fr_kso )
      END IF
    END IF
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (soiltyp_subs(i) >= 3) THEN
      ! generalized upper boundary condition
      zagc(i,1) = zagc(i,1)/zagb(i,1)
      zagd(i,1) = zagd(i,1)/zagb(i,1)
    END IF
  END DO

  !$ACC LOOP SEQ
  DO kso=2,ke_soil_hy-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zzz)
    DO i = ivstart, ivend
      IF (soiltyp_subs(i) >= 3) THEN
        zzz = 1._wp/(zagb(i,kso) - zaga(i,kso)*zagc(i,kso-1))
        zagc(i,kso) = zagc(i,kso) * zzz
        zagd(i,kso) = (zagd(i,kso) - zaga(i,kso)*zagd(i,kso-1)) * zzz
      END IF
    END DO
  END DO                ! soil layers


  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (soiltyp_subs(i) >= 3) THEN
      IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) THEN
        zage(i,ke_soil_hy_m) = (zagd(i,ke_soil_hy_m)-zaga(i,ke_soil_hy_m)*     &
                                zagd(i,ke_soil_hy_m-1))/                       &
                               (zagb(i,ke_soil_hy_m) - zaga(i,ke_soil_hy_m)*   &
                                zagc(i,ke_soil_hy_m-1))
      ELSE
        zage(i,ke_soil_hy) = (zagd(i,ke_soil_hy)-zaga(i,ke_soil_hy)*     &
                              zagd(i,ke_soil_hy-1))/                     &
                             (zagb(i,ke_soil_hy) - zaga(i,ke_soil_hy)*   &
                              zagc(i,ke_soil_hy-1))
      END IF
    END IF
  END DO

  IF (itype_mire == 1) THEN
!US index for k-loop depends on i! how to handle that???

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (soiltyp_subs(i) >= 3) THEN
        DO kso = ke_soil_hy_b(i)-1,1,-1
          zage(i,kso)     = zagd(i,kso) - zagc(i,kso)*zage(i,kso+1)

          ! compute implicit part of new liquid water content and add existing
          ! ice content
          w_so_new(i,kso) = zage(i,kso)*zdzhs(kso) + w_so_ice_now(i,kso)
        END DO
      END IF
    END DO
  ELSE
    !$ACC LOOP SEQ
    DO kso = ke_soil_hy-1,1,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (soiltyp_subs(i) >= 3) THEN
          zage(i,kso)     = zagd(i,kso) - zagc(i,kso)*zage(i,kso+1)

          ! compute implicit part of new liquid water content and add existing
          ! ice content
          w_so_new(i,kso) = zage(i,kso)*zdzhs(kso) + w_so_ice_now(i,kso)
        END IF
      END DO
    END DO                ! soil layers
  END IF

!lowest active hydrological level
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (soiltyp_subs(i) >= 3) THEN
      ! boundary values ensure that the calculation below leaves the climate
      ! layer water contents unchanged compute implicit part of new liquid
      ! water content and add existing ice content
      IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) THEN
       w_so_new(i,ke_soil_hy_m) = zage(i,ke_soil_hy_m)*zdzhs(ke_soil_hy_m) + &
                                 w_so_ice_now(i,ke_soil_hy_m)
      ELSE
        w_so_new(i,ke_soil_hy) = zage(i,ke_soil_hy)*zdzhs(ke_soil_hy) + &
                                 w_so_ice_now(i,ke_soil_hy)
      END IF
    END IF
  END DO

  ! to ensure vertical constant water concentration profile beginning at
  ! layer ke_soil_hy_b(i) for energetic treatment only
  ! soil water climate layer(s)

  ! J. Helmert: Activate the ground water table with water diffusion from layers beyond ke_soil_hy
  !             for mire points.
  !             This should avoid a dry falling of the upper mire layers due to strong bare soil evaporation!

  ! ground water as lower boundary of soil column

  IF (itype_hydbound == 3) THEN
    !$ACC LOOP SEQ
    DO kso = ke_soil_hy+1,ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (soiltyp_subs(i) >= 3) THEN
          w_so_new(i,kso) = zporv(i,kso)*zdzhs(kso)
        END IF
      END DO
    END DO
  ELSE
    !$ACC LOOP SEQ
    DO kso = ke_soil_hy+1,ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (soiltyp_subs(i) >= 3) THEN
          w_so_new(i,kso) = w_so_new(i,kso-1)*zdzhs(kso)/zdzhs(kso-1)
        END IF
      END DO
    END DO
  ENDIF

  IF (itype_mire == 1) THEN
    !$ACC LOOP SEQ
    DO kso = ke_soil_hy_m+1,ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (soiltyp_subs(i) == 8) THEN
          w_so_new(i,kso) = zporv(i,kso)*zdzhs(kso)
        END IF
      END DO
    END DO
  END IF

  ! combine implicit part of sedimentation and capillary flux with explicit part
  ! (for soil water flux investigations only)
  !$ACC LOOP SEQ
  DO kso = 2,ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR &
    !$ACC   PRIVATE(zice_fr_ksom1, zice_fr_kso, zlw_fr_ksom1_new) &
    !$ACC   PRIVATE(zlw_fr_kso_new, zlw_fr_ksom1, zlw_fr_kso, zfr_ice, zredm05) &
    !$ACC   PRIVATE(zlw_fr_ksom05, zdlw_fr_ksom05, zklw_fr_ksom05, zredm) &
    !$ACC   PRIVATE(zklw_fr_kso_new, zdelta_sm, zdlw_fr_kso, zklw_fr_kso) &
    !$ACC   PRIVATE(zklw_fr_ksom1, zdhydcond_dlwfr)
    DO i = ivstart, ivend
      IF (soiltyp_subs(i) >= 3) THEN
        zice_fr_ksom1 = ziw_fr(i,kso-1)
        zice_fr_kso   = ziw_fr(i,kso)
        zlw_fr_ksom1_new= w_so_new(i,kso-1)/zdzhs(kso-1) - zice_fr_ksom1
        zlw_fr_kso_new  = w_so_new(i,kso  )/zdzhs(kso  ) - zice_fr_kso
        zlw_fr_ksom1  = w_so_now(i,kso-1)/zdzhs(kso-1) - zice_fr_ksom1
        zlw_fr_kso    = w_so_now(i,kso  )/zdzhs(kso  ) - zice_fr_kso
        !... additionally for runoff_g at lower level of lowest active water
        ! layer calculated with (upstream) main level soil water content
        ! compute reduction factor for transport coefficients
        zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
        zredm05 = 1._wp-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)

        ! interpolated liquid water content at interface to layer above
        zlw_fr_ksom05 =0.5_wp*(zdzhs(kso)*zlw_fr_ksom1+zdzhs(kso-1)*zlw_fr_kso) /zdzms(kso)
          !!$    zdlw_fr_ksom05= zredm05*zdw(i,kso)*EXP(zdw1(i,kso) *  &
          !!$                    (zporv(i,kso)-zlw_fr_ksom05)/(zporv(i,kso)-zadp(i,kso)) )
          !!$    zklw_fr_ksom05= zredm05*zkw(i,kso) * EXP(zkw1(i,kso)* &
          !!$                    (zporv(i,kso)-zlw_fr_ksom05)/(zporv(i,kso)-zadp(i,kso)) )
          !!$

        zdlw_fr_ksom05= zredm05*watrdiff_RT(zdw(i,kso),zlw_fr_ksom05,&
                               zdw1(i,kso),zporv(i,kso),zadp(i,kso))

        zklw_fr_ksom05= zredm05*watrcon_RT(zkw(i,kso-1),zlw_fr_ksom05,zkw1(i,kso),&
                          zporv(i,kso),zadp(i,kso))

        IF (soiltyp_subs(i) == 8 .AND. itype_mire == 1) THEN
          IF (kso> ke_soil_hy_m) zdlw_fr_ksom05=0.0_wp   ! no flux gradient
                                                         ! contribution below 2.5m
          IF (kso> ke_soil_hy_m) zklw_fr_ksom05=0.0_wp   ! no gravitation flux below 2.5m
        ELSE

          IF (kso> ke_soil_hy) zdlw_fr_ksom05=0.0_wp   ! no flux gradient
                                                       ! contribution below 2.5m
          IF (kso> ke_soil_hy) zklw_fr_ksom05=0.0_wp   ! no gravitation flux below 2.5m
        END IF

        zflmg(i,kso) = (1._wp-zalfa) * zflmg(i,kso) + & ! explicit flux component
                                  zalfa  * rho_w *    & ! implicit flux component
                     (zdlw_fr_ksom05 * (zlw_fr_kso_new+zice_fr_kso-zlw_fr_ksom1_new-zice_fr_ksom1) &
                       /zdzms(kso) - zklw_fr_ksom05)
        zredm = 1._wp-zice_fr_kso/(zlw_fr_kso+zice_fr_kso)
          !!$    zklw_fr_kso_new = zredm*zkw(i,kso) * EXP(zkw1(i,kso)* &
          !!$                      (zporv(i,kso) - zlw_fr_kso_new)/(zporv(i,kso) - zadp(i,kso)) )

        zklw_fr_kso_new = zredm*watrcon_RT(zkwm(i,kso),&
             zlw_fr_kso_new,zkw1(i,kso),zporv(i,kso),zadp(i,kso))

        ! actual gravitation water flux
        IF (w_so_new(i,kso).LT.1.01_wp*zadp(i,kso)*zdzhs(kso)) THEN
          zklw_fr_kso_new=0._wp
        ENDIF
        zrunoff_grav(i,kso) =  - rho_w * zklw_fr_kso_new

        ! ground water as lower boundary of soil column
        IF ((kso == ke_soil_hy+1).and.(itype_hydbound == 3)) THEN
          zdelta_sm=( zlw_fr_kso_new - zlw_fr_ksom1_new )

              !!$   zdlw_fr_kso = zredm05*zdw(i,kso)*EXP(zdw1(i,kso) *  &
              !!$        (zporv(i,kso)-zlw_fr_kso_new)/(zporv(i,kso)-zadp(i,kso)) )
              !!$   zklw_fr_kso = zredm05*zkw(i,kso) * EXP(zkw1(i,kso)* &
              !!$        (zporv(i,kso)-zlw_fr_kso_new)/(zporv(i,kso)-zadp(i,kso)) )
              !!$   zklw_fr_ksom1 = zredm05*zkw(i,kso) * EXP(zkw1(i,kso)* &
              !!$        (zporv(i,kso)-zlw_fr_ksom1_new)/(zporv(i,kso)-zadp(i,kso)) )

          ! GZ, 2017-11-07: The variable staggering in the following expressions needs to be checked
          !   zredm05 is used for interface levels kso and kso-1; logically, "zredp05" would be 
          !   needed in the first two expressions in addition, the indexing of zkw1 would be 
          !   incorrect if this parameter is made level-dependent
          zdlw_fr_kso = zredm05*watrdiff_RT(zdw(i,kso),zlw_fr_kso_new,    &
                          zdw1(i,kso),zporv(i,kso),zadp(i,kso))
          zklw_fr_kso = zredm05*watrcon_RT(zkw(i,kso),                    &
                          zlw_fr_kso_new,zkw1(i,kso),zporv(i,kso),zadp(i,kso))
          zklw_fr_ksom1 = zredm05*watrcon_RT(zkw(i,kso-1),&
                          zlw_fr_ksom1_new,zkw1(i,kso),zporv(i,kso),zadp(i,kso))

          zdhydcond_dlwfr=( zklw_fr_kso - zklw_fr_ksom1 ) / zdelta_sm
          zrunoff_grav(i,ke_soil_hy)=zrunoff_grav(i,ke_soil_hy)+ zdhydcond_dlwfr / &
             (1.0_wp-EXP(-zdhydcond_dlwfr/zdlw_fr_kso*0.5_wp*zdzms(ke_soil_hy+1)))* zdelta_sm
        ENDIF
        IF ((kso == ke_soil_hy_m+1) .AND. (itype_mire == 1 .AND. soiltyp_subs(i) == 8)) THEN
          zdelta_sm=( zlw_fr_kso_new - zlw_fr_ksom1_new )

              !!$   zdlw_fr_kso = zredm05*zdw(i,kso)*EXP(zdw1(i,kso) *  &
              !!$        (zporv(i,kso)-zlw_fr_kso_new)/(zporv(i,kso)-zadp(i,kso)) )
              !!$   zklw_fr_kso = zredm05*zkw(i,kso) * EXP(zkw1(i,kso)* &
              !!$        (zporv(i,kso)-zlw_fr_kso_new)/(zporv(i,kso)-zadp(i,kso)) )
              !!$   zklw_fr_ksom1 = zredm05*zkw(i,kso) * EXP(zkw1(i,kso)* &
              !!$        (zporv(i,kso)-zlw_fr_ksom1_new)/(zporv(i,kso)-zadp(i,kso)) )

          ! GZ, 2017-11-07: The variable staggering in the following expressions needs to be checked
          !   zredm05 is used for interface levels kso and kso-1; logically, "zredp05" would be 
          !   needed in the first two expressions in addition, the indexing of zkw1 would be 
          !   incorrect if this parameter is made level-dependent
          zdlw_fr_kso = zredm05*watrdiff_RT(zdw(i,kso),zlw_fr_kso_new,    &
                          zdw1(i,kso),zporv(i,kso),zadp(i,kso))
          zklw_fr_kso = zredm05*watrcon_RT(zkw(i,kso),                    &
                          zlw_fr_kso_new,zkw1(i,kso),zporv(i,kso),zadp(i,kso))
          zklw_fr_ksom1 = zredm05*watrcon_RT(zkw(i,kso-1),&
                          zlw_fr_ksom1_new,zkw1(i,kso),zporv(i,kso),zadp(i,kso))

          zdhydcond_dlwfr=( zklw_fr_kso - zklw_fr_ksom1 ) / zdelta_sm
          zrunoff_grav(i,ke_soil_hy_m)=zrunoff_grav(i,ke_soil_hy_m)+ zdhydcond_dlwfr / &
             (1.0_wp-exp(-zdhydcond_dlwfr/zdlw_fr_kso*0.5_wp*zdzms(ke_soil_hy_m+1)))* zdelta_sm
        END IF
      END IF
    END DO
  END DO

  IF (itype_mire == 1) THEN
    !$ACC LOOP SEQ PRIVATE(zro_sfak, zro_gfak)
    DO  kso = 1,ke_soil
      ! utility variables used to avoid if-constructs in following loops
      zro_sfak = zsf_heav(0.5_wp + REAL(msr_off - kso,wp))      ! 1.0 for 'surface runoff'
      zro_gfak = 1._wp - zro_sfak                               ! 1.0 for 'ground runoff'

      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdwg, zredfu, zro, zwgn, zro2, zkorr, zfmb_fak)
      DO i = ivstart, ivend
        ! sedimentation and capillary transport in soil
        IF (soiltyp_subs(i) >= 3) THEN
          ! - Compute subsoil runoff (runoff_g) as drainage flux through bottom
          !   of layer ke_soil_hy (as suggested by the Rhone-Aggregation
          !   Experiment)
          ! - soil moisture gradient related flux is switched off below
          !   (i.e. only sedimentation flux allowed between ke_soil_hy and ke_soil_hy+1)
          IF (soiltyp_subs(i) == 8)  THEN
            zfmb_fak = MERGE(1.0_wp, 0.0_wp, kso==ke_soil_hy_m)
            zcou_roffg = MERGE(1.0_wp, 0.0_wp, kso<=ke_soil_hy_m)
          ELSE
            zfmb_fak = MERGE(1.0_wp, 0.0_wp, kso==ke_soil_hy)
            zcou_roffg = MERGE(1.0_wp, 0.0_wp, kso<=ke_soil_hy)
          END IF

          ! first runoff calculation without consideration of
          ! evapotranspiration
          !zdwg   =  zflmg(i,kso+1) - zflmg(i,kso)
          !zdwg calculated above by flux divergence has to be aequivalent with
          zdwg =  (w_so_new(i,kso)/zdzhs(kso)-zw_fr(i,kso))*zdzhs(kso) / zdtdrhw
          zdwg =  zdwg + zrunoff_grav(i,kso)*zfmb_fak
          zredfu =  MAX( 0.0_wp, MIN( 1.0_wp,(zw_fr(i,kso) -     &
                           zfcap(i,kso))/MAX(zporv(i,kso) - zfcap(i,kso),eps_div)) )
          zredfu = zsf_heav(zdwg)*zredfu
          zro    = zdwg*zredfu
          zdwg   = zdwg*(1._wp - zredfu)

          ! add evaporation (znlgw1f: first layer only)
          ! and transpiration (for each layer)
          zdwg   = zdwg + znlgw1f(kso) * zesoil(i) + ztrang (i,kso)
          zwgn   = zw_fr(i,kso) + zdtdrhw*zdwg/zdzhs(kso)
          zro2   = zrhwddt*zdzhs(kso)*MAX(0.0_wp, zwgn - zporv(i,kso))
          zkorr  = zrhwddt*zdzhs(kso)*MAX(0.0_wp, zadp(i,kso) - zwgn )
          zdwgdt(i,kso)= zdwg + zkorr - zro2
          zro    = zro      + zro2
          runoff_s(i) = runoff_s(i) + zro*zro_sfak*zroffdt
          ! only count runoff_g in the hydrological active layers
          runoff_g(i) = runoff_g(i) + zcou_roffg*zro*zro_gfak*zroffdt
          ! runoff_g reformulation due to drainage flux through bottom of layer ke_soil_hy
          !  zfmb_fak is only 1 at the last active horizont and 0 elsewhere
          !  only at this level the gravitational settling has to be counted for runoff (bulk view)
          runoff_g(i) = runoff_g(i) - (zrunoff_grav(i,kso) * zfmb_fak &
                                             + zkorr) * zroffdt
          !if (runoff_g(i) > 1000.0_wp) then
          !  print *, 'RUNOFF_G:  ', i, runoff_g(i), zrunoff_grav(i,kso), zro, zroffdt, zkorr
          !endif
        END IF
      END DO
    END DO         ! end loop over soil layers
  ELSE
    !$ACC LOOP SEQ PRIVATE(zro_sfak, zro_gfak, zfmb_fak, zcou_roffg)
!$NEC novector
    DO  kso = 1,ke_soil
      ! utility variables used to avoid if-constructs in following loops
      zro_sfak = zsf_heav(0.5_wp + REAL(msr_off - kso,wp))      ! 1.0 for 'surface runoff'
      zro_gfak = 1._wp - zro_sfak                               ! 1.0 for 'ground runoff'

      ! - Compute subsoil runoff (runoff_g) as drainage flux through bottom
      !   of layer ke_soil_hy (as suggested by the Rhone-Aggregation
      !   Experiment)
      ! - soil moisture gradient related flux is switched off below
      !   (i.e. only sedimentation flux allowed between ke_soil_hy and ke_soil_hy+1)
      zfmb_fak = MERGE(1.0_wp, 0.0_wp, kso==ke_soil_hy)
      zcou_roffg = MERGE(1.0_wp, 0.0_wp, kso<=ke_soil_hy)

      ! sedimentation and capillary transport in soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdwg, zredfu, zro, zwgn, zro2, zkorr)
      DO i = ivstart, ivend
        IF (soiltyp_subs(i) >= 3) THEN
          ! first runoff calculation without consideration of
          ! evapotranspiration
          !zdwg   =  zflmg(i,kso+1) - zflmg(i,kso)
          !zdwg calculated above by flux divergence has to be aequivalent with
          zdwg =  (w_so_new(i,kso)/zdzhs(kso)-zw_fr(i,kso))*zdzhs(kso) / zdtdrhw
          zdwg =  zdwg + zrunoff_grav(i,kso)*zfmb_fak
          zredfu =  MAX( 0.0_wp, MIN( 1.0_wp,(zw_fr(i,kso) -     &
                           zfcap(i,kso))/MAX(zporv(i,kso) - zfcap(i,kso),eps_div)) )
          zredfu = zsf_heav(zdwg)*zredfu
          zro    = zdwg*zredfu
          zdwg   = zdwg*(1._wp - zredfu)

          ! add evaporation (znlgw1f: first layer only)
          ! and transpiration (for each layer)
          zdwg   = zdwg + znlgw1f(kso) * zesoil(i) + ztrang (i,kso)
          zwgn   = zw_fr(i,kso) + zdtdrhw*zdwg/zdzhs(kso)
          zro2   = zrhwddt*zdzhs(kso)*MAX(0.0_wp, zwgn - zporv(i,kso))
          zkorr  = zrhwddt*zdzhs(kso)*MAX(0.0_wp, zadp(i,kso) - zwgn )
          zdwgdt(i,kso)= zdwg + zkorr - zro2
          zro    = zro      + zro2
          runoff_s(i) = runoff_s(i) + zro*zro_sfak*zroffdt
          ! only count runoff_g in the hydrological active layers
          runoff_g(i) = runoff_g(i) + zcou_roffg*zro*zro_gfak*zroffdt
          ! runoff_g reformulation due to drainage flux through bottom of layer ke_soil_hy
          !  zfmb_fak is only 1 at the last active horizont and 0 elsewhere
          !  only at this level the gravitational settling has to be counted for runoff (bulk view)
          runoff_g(i) = runoff_g(i) - (zrunoff_grav(i,kso) * zfmb_fak &
                                            + zkorr) * zroffdt
        END IF
      END DO
    END DO         ! end loop over soil layers
  ENDIF
  !$ACC END PARALLEL

!------------------------------------------------------------------------------
! Section II.5: Soil surface heat flux (thermal forcing)
!------------------------------------------------------------------------------

  IF (lmulti_snow) THEN

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      ! Estimate thermal surface fluxes:
      ! Estimate thermal surface fluxes over snow covered and snow free
      ! part of surface based on area mean values calculated in radiation
      ! code (positive = downward)

!     IF (lterra_urb) THEN
!       ! modification of the infrared albedo according to the buliding fraction
!       ztalb = sa_uf(i) * ctalb_bm * alb_red_uf(i) + (1.0_wp - sa_uf(i)) * Ctalb
!     ELSE
        ztalb = Ctalb
!     END IF

      zgstr =   sigma*(1._wp - ztalb) * ( (1._wp - zf_snow(i))* &
                zts(i) + zf_snow(i)*ztsnow_mult(i,1) )**4 + thbs(i)
      zthsnw(i) = - sigma*(1._wp - ztalb)*ztsnow_mult(i,1)**4 + zgstr
      zthsoi(i) = - sigma*(1._wp - ztalb)*zts(i)**4 + zgstr
      ! the estimation of the solar component would require the availability
      ! of the diffuse and direct components of the solar flux
      !
      ! Forcing for snow-free soil:
      ! (evaporation, transpiration, formation of dew and rime are already
      !  weighted by correspondind surface fraction)
      ! net radiation, sensible and latent heat flux

      zrnet_s(i) = sobs(i) + zthsoi(i)
      zshfl_s(i) = cp_d*zrhoch(i) * (zth_low(i) - zts(i))
      zlhfl_s(i) = (zts_pm(i)*lh_v + (1._wp-zts_pm(i))*lh_s)*zverbo(i) &
                   / MAX(eps_div,(1._wp - zf_snow(i)))  ! take out (1-f) scaling
      zqhfl_s(i) = zverbo(i)/ MAX(eps_div,(1._wp - zf_snow(i)))  ! take out (1-f) scaling
    END DO


    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      ! thawing of snow falling on soil with Ts > T0
      IF (ztsnow_pm(i)*zrs(i) > 0.0_wp) THEN
        ! snow fall on soil with T>T0, snow water content increases
        ! interception store water content
        zsprs  (i) = - lh_f*zrs(i)
        zdwidt (i) = zdwidt (i) + zrs(i)
        zdwsndt(i) = zdwsndt(i) - zrs(i)

        ! avoid overflow of interception store, add possible excess to
        ! surface run-off
        zwinstr(i)      = zwin(i) + zdwidt(i)*zdtdrhw
        IF (zwinstr(i) > zwimax(i)) THEN  ! overflow of interception store
          zro        = (zwinstr(i) - zwimax(i))*zrhwddt
          zdwidt(i)= zdwidt(i) - zro
          runoff_s(i)  = runoff_s(i) + zro*zroffdt
        ENDIF                       ! overflow of interception store

      ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
      ELSEIF ((1._wp-zts_pm(i))*zrr(i) > 0.0_wp) THEN
        zsprs  (i) = lh_f*zrr(i)
        ! keep freezing rain in interception storage rather than shifting it to snow
        ! zdwidt (i) = zdwidt (i) - zrr(i)
        ! zdwsndt(i) = zdwsndt(i) + zrr(i)
      ELSE
        zsprs  (i) = 0.0_wp
      END IF

      ! Influence of heatflux through snow on total forcing:
      zdwsndt(i) = zdwsndt(i)*zsf_heav(zdwsndt(i) - eps_soil/zdtdrhw)
      zwsnew(i)  = zwsnow(i) + zdwsndt(i)*zdtdrhw
      IF (zwsnew(i).GT.eps_soil) THEN

        zrho_snowf = crhosminf+(crhosmaxf-crhosminf)* (zth_low(i)-csnow_tmin) &
                                                /(t0_melt          -csnow_tmin)
        zrho_snowf = MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))

        IF(zdwsndt(i)-zrs(i)-zrr(i).GT.0.0_wp) THEN

          wtot_snow_now(i,1) = MAX(wtot_snow_now(i,1) + zdwsndt(i)*zdtdrhw, 0.0_wp)

          zhm_snow(i,1) = zhm_snow(i,1) - (zdwsndt(i)-zrs(i)-                             &
             zrr(i))*zdt/rho_i/2._wp- zrs(i)*zdt/zrho_snowf/2._wp- zrr(i)*zdt/rho_i/2._wp
          zdzh_snow(i,1) = zdzh_snow(i,1) + (zdwsndt(i)-zrs(i)-zrr(i))*zdt/rho_i +        &
             zrs(i)*zdt/zrho_snowf + zrr(i)*zdt/rho_i

          rho_snow_mult_now(i,1) = MAX(wtot_snow_now(i,1)*rho_w/zdzh_snow(i,1), 0.0_wp)
        ELSE

          wtot_snow_now(i,1) = MAX(wtot_snow_now(i,1) + (zrs(i)+zrr(i))*zdtdrhw, 0.0_wp)

          zhm_snow(i,1)  = zhm_snow(i,1) - zrs(i)*zdt/zrho_snowf/2._wp- &
                             zrr(i)*zdt/rho_i/2._wp
          zdzh_snow(i,1) = zdzh_snow(i,1) + zrs(i)*zdt/zrho_snowf + zrr(i)*zdt/rho_i

          IF (wtot_snow_now(i,1) .GT. 0._wp) THEN
            rho_snow_mult_now(i,1) = MAX(wtot_snow_now(i,1)*rho_w/zdzh_snow(i,1), 0.0_wp)

            wtot_snow_now(i,1) = MAX(wtot_snow_now(i,1) &
                                  + (zdwsndt(i)-zrs(i)-zrr(i))*zdtdrhw,0.0_wp)

            zhm_snow(i,1)  = zhm_snow(i,1) - (zdwsndt(i)-zrs(i)-zrr(i)) &
                                 *zdt/rho_snow_mult_now(i,1)/2._wp
            zdzh_snow(i,1) = zdzh_snow(i,1) + (zdwsndt(i)-zrs(i)-zrr(i)) &
                                 *zdt/rho_snow_mult_now(i,1)
          ELSE
            rho_snow_mult_now(i,1) = 0.0_wp
            zdzh_snow(i,1) = 0.0_wp
          END IF

        END IF
      END IF
      h_snow_now(i) = 0.0_wp
      sum_weight(i) = 0.0_wp
    END DO
 
    !$ACC LOOP SEQ
    DO ksn = 1,ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (zwsnew(i).GT.eps_soil) THEN
          h_snow_now(i) = h_snow_now(i) + zdzh_snow(i,ksn)
        END IF
      END DO
    END DO

    k = MIN(2,ke_snow-1)
    !$ACC LOOP SEQ
    DO ksn = 1,ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (zwsnew(i) .GT. eps_soil) THEN
          IF (zwsnow(i) .GT. eps_soil) THEN
            IF (ksn == 1) THEN ! Limit top layer to max_toplaydepth
              zhh_snow(i,ksn) = -MAX( h_snow_now(i)-max_toplaydepth, h_snow_now(i)/ke_snow*(ke_snow-ksn) )
            ELSE IF (ksn == 2 .AND. ke_snow > 2) THEN ! Limit second layer to 8*max_toplaydepth
              zhh_snow(i,ksn) = MIN( 8._wp*max_toplaydepth+zhh_snow(i,1), zhh_snow(i,1)/(ke_snow-1)*(ke_snow-ksn) )
            ELSE ! distribute the remaining snow equally among the layers
              zhh_snow(i,ksn) = zhh_snow(i,k)/(ke_snow-k)*(ke_snow-ksn)
            ENDIF
          ELSE ! a newly generated snow cover will not exceed max_toplaydepth
            zhh_snow(i,ksn) = -h_snow_now(i)/ke_snow*(ke_snow-ksn)
          END IF
        END IF
      END DO
    END DO


    !$ACC LOOP SEQ
    DO ksn = ke_snow,1,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF(zwsnew(i) .GT. eps_soil .AND. zwsnow(i) .GT. eps_soil) THEN
          dz_old(i,ksn) = zdzh_snow(i,ksn)
          z_old(i,ksn) = -sum_weight(i) - zdzh_snow(i,ksn)/2._wp
          sum_weight(i) = sum_weight(i) + zdzh_snow(i,ksn)
        END IF
      END DO
    END DO


    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF(zwsnew(i) .GT. eps_soil) THEN
        zhm_snow (i,1) = (-h_snow_now(i) + zhh_snow(i,1))/2._wp
        zdzh_snow(i,1) = zhh_snow(i,1) + h_snow_now(i)   !layer thickness betw. half levels of uppermost snow layer
        zdzm_snow(i,1) = zhm_snow(i,1) + h_snow_now(i)   !layer thickness between snow surface and 
                                                         !  main level of uppermost layer
        IF(zwsnow(i) .GT. eps_soil) THEN
          IF(dz_old(i,1).ne.0..and.rho_snow_mult_now(i,1).ne.0.) THEN
            wliq_snow_now(i,1) = wliq_snow_now(i,1)/dz_old(i,1)
          END IF
        END IF
      END IF
    END DO


    !$ACC LOOP SEQ
    DO ksn = 2,ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF(zwsnew(i) .GT. eps_soil) THEN
          zhm_snow(i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._wp
          zdzh_snow(i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
          zdzm_snow(i,ksn) = zhm_snow(i,ksn) - zhm_snow(i,ksn-1) ! layer thickness betw. main levels
          IF(zwsnow(i) .GT. eps_soil) THEN
            IF(dz_old(i,ksn).ne.0..and.rho_snow_mult_now(i,ksn).ne.0.) THEN
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn)/dz_old(i,ksn)
            END IF
          END IF
        END IF
      END DO
    END DO

    !$ACC LOOP SEQ
    DO ksn = ke_snow,1,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        t_new  (i,ksn) = 0.0_wp
        rho_new(i,ksn) = 0.0_wp
        wl_new (i,ksn) = 0.0_wp
      END DO

      !$ACC LOOP SEQ
      DO k = ke_snow,1,-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i = ivstart, ivend
          IF(zwsnew(i) .GT. eps_soil .AND. zwsnow(i) .GT. eps_soil) THEN

            weight = MIN(dz_old(i, k), &
                 &       z_old(i, k) + dz_old(i, k)*0.5_wp &
                 &         - zhm_snow(i, ksn) + zdzh_snow(i, ksn)*0.5_wp,&
                 &       zhm_snow(i, ksn) + zdzh_snow(i,ksn)*0.5_wp &
                 &         - z_old(i, k) + dz_old(i, k)*0.5_wp, &
                 &       zdzh_snow(i,ksn))

            weight = (weight + ABS(weight)) * 0.5_wp
            weight = weight / zdzh_snow(i,ksn)
            t_new  (i,ksn) = t_new  (i,ksn) + ztsnow_mult  (i,k      )*weight
            rho_new(i,ksn) = rho_new(i,ksn) + rho_snow_mult_now(i,k)*weight
            wl_new (i,ksn) = wl_new (i,ksn) + wliq_snow_now(i,k)*weight
          END IF
        END DO
      END DO
    END DO

    ! MIN(z_old(i, k) + dz_old(i, k)/2._wp, &
    ! &   zhm_snow(i, ksn) + zdzh_snow(i,ksn)/2._wp) &
    ! - MAX(z_old(i, k) - dz_old(i, k)/2._wp, &
    ! &     zhm_snow(i, ksn) - zdzh_snow(i, ksn)/2._wp), &

    !$ACC LOOP SEQ
    DO ksn = ke_snow,1,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF(zwsnew(i) .GT. eps_soil) THEN
          IF(zwsnow(i) .GT. eps_soil) THEN
            ztsnow_mult  (i,ksn      ) = t_new  (i,ksn)
            rho_snow_mult_now(i,ksn) = rho_new(i,ksn)
            wtot_snow_now    (i,ksn) = rho_new(i,ksn)*zdzh_snow(i,ksn)/rho_w
            wliq_snow_now    (i,ksn) = wl_new (i,ksn)*zdzh_snow(i,ksn)
          ELSE ! Remark: if there was now snow in the previous time step, snow depth 
               ! will not exceed the limit for equipartitioning
            ztsnow_mult  (i,ksn      ) = t_s_now(i)
            rho_snow_mult_now(i,ksn) = rho_snow_mult_now(i,1)
            wtot_snow_now    (i,ksn) = zwsnew(i)/ke_snow
          END IF
        END IF
      END DO
    END DO


    ! heat conductivity of snow as funtion of water content
    !$ACC LOOP SEQ
    DO ksn = 1, ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (zwsnew(i).GT.eps_soil) THEN
          zalas_mult(i,ksn) = 2.22_wp*EXP(1.88_wp*LOG(rho_snow_mult_now(i,ksn)/rho_i))
        END IF
      END DO
    END DO


    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (zwsnew(i).GT.eps_soil) THEN
        zgsb(i) = ((zalas_mult(i,ke_snow)*(-zhm_snow(i,ke_snow))+hzalam(i,1)*zdzms(1))/ &
                   (-zhm_snow(i,ke_snow)+zdzms(1)) * &
                   (ztsnow_mult(i,ke_snow) - t_so_now(i,1))/(-zhm_snow(i,ke_snow) &
                   +zdzms(1)))*zf_snow(i)

        !  ! GZ: use formulation of single-layer snow model, which is numerically more stable
        !  zgsb(i) = zalas_mult(i,ke_snow)*(ztsnow_mult(i,ke_snow) - t_so_now(i,1))/ &
        !  MAX(-zhm_snow(i,ke_snow),cdsmin)

      END IF

      ! total forcing for uppermost soil layer
        !<em new solution
        !        zfor_s(i) = ( zrnet_s(i) + zshfl_s(i) + zlhfl_s(i) + zsprs(i) ) * (1._wp - zf_snow(i)) &
        !                    + (1._wp-ztsnow_pm(i)) * zgsb(i)
      zfor_s(i) = ( zrnet_s(i) + zshfl_s(i) + zlhfl_s(i) + zsprs(i) ) * (1._wp - zf_snow(i))
        !em>

      IF(zwsnew(i) .GT. eps_soil) THEN
        IF(zextinct(i,1).gt.0.0_wp) THEN
          zrnet_snow = zthsnw(i)
        ELSE
          zrnet_snow = sobs(i) + zthsnw(i)
        END IF
      ELSE
        zrnet_snow = sobs(i) + zthsnw(i)
      END IF
      zshfl_snow(i) = zrhoch(i)*cp_d*(zth_low(i) - ztsnow_mult(i,1))
      zlhfl_snow(i) = lh_s*zversn(i)
      zqhfl_snow(i) = zversn(i)
      zfor_snow_mult(i)  = (zrnet_snow + zshfl_snow(i) + zlhfl_snow(i) + lh_f*zrr(i))*zf_snow(i)
    END DO

  !$ACC END PARALLEL
  ELSE  ! single-layer snow model

    ! TERRA_URB: Set anthropogenic heat flux
    IF (lterra_urb .AND. itype_ahf >= 1) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i = ivstart, ivend
        zahf(i) = urb_isa(i)*ahf(i)
      END DO
      !$ACC END PARALLEL
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP GANG VECTOR PRIVATE(ztalb, zgstr, zro, zalas, zzz, zd, zfr_melt, zw_ovpv)
    !$NEC sparse
    DO i = ivstart, ivend
      ! Estimate thermal surface fluxes:
      ! Estimate thermal surface fluxes over snow covered and snow free
      ! part of surface based on area mean values calculated in radiation
      ! code (positive = downward)

      ztalb = Ctalb

      zgstr     =   sigma*(1._wp - ztalb) * ( (1._wp - zf_snow(i))* &
                    ztsk(i) + zf_snow(i)*ztsnow(i) )**4 + thbs(i)
      zthsoi(i) = - sigma*(1._wp - ztalb)*ztsk(i)**4 + zgstr
      zthsnw(i) = - sigma*(1._wp - ztalb)*ztsnow(i)**4 + zgstr

      ! the estimation of the solar component would require the availability
      ! of the diffuse and direct components of the solar flux
      !
      ! Forcing for snow-free soil:
      ! (evaporation, transpiration, formation of dew and rime are already
      !  weighted by correspondind surface fraction)
      ! net radiation, sensible and latent heat flux

      zrnet_s(i) = sobs(i) + zthsoi(i)
      zshfl_s(i) = cp_d*zrhoch(i) * (zth_low(i) - ztsk(i))
      zlhfl_s(i) = (ztsk_pm(i)*lh_v + (1._wp-ztsk_pm(i))*lh_s)*zverbo(i) &
                   / MAX(eps_div,(1._wp - zf_snow(i)))  ! take out (1-f) scaling
      zqhfl_s(i) = zverbo(i)/ MAX(eps_div,(1._wp - zf_snow(i)))  ! take out (1-f) scaling
      zsprs  (i) = 0.0_wp

      ! thawing of snow falling on soil with Ts > T0
      !
      ! it turned out that ensuring (ztsnow_pm(i)*zrs(i) > 0.0_wp) is not sufficient, as in rare cases
      ! a floating point underflow might occur for zrs. By itself this is not a huge problem, however computing 
      ! 1/zrs(i) (as it is done below) results in a floating point overflow in these cases.
      IF (ztsnow_pm(i)*zrs(i) > eps_nounderflow) THEN
        ! snow fall on soil with T>T0, snow water content increases interception store water content
        ! melting rate is limited such that the two upper soil levels are not cooled significantly below the freezing point
        zzz = t0_melt-0.25_wp
        zd  = (t_so_now(i,1)-zzz)*zroc(i,1)*zdzhs(1) + MAX(0._wp,t_so_now(i,2)-zzz)*zroc(i,2)*zdzhs(2)
        zfr_melt = zd/MAX(zdt*lh_f*zrs(i), zd)
        zsprs  (i) = - lh_f*zrs(i)*zfr_melt
        zdwidt (i) = zdwidt (i) + zrs(i)*zfr_melt
        zdwsndt(i) = zdwsndt(i) - zrs(i)*zfr_melt

        ! avoid overflow of interception store, add possible excess to
        ! surface run-off
        zwinstr(i) = zwin(i) + zdwidt(i)*zdtdrhw
        IF (zwinstr(i) > zwimax(i)) THEN  ! overflow of interception store
          zro         = (zwinstr(i) - zwimax(i))*zrhwddt
          zdwidt(i)   = zdwidt(i) - zro
          zdwgdt(i,1) = zdwgdt(i,1) + zro
          ! check for pore volume overflow and add the remainder to surface runoff
          zw_ovpv = MAX(0.0_wp, (zw_fr(i,1)-zporv(i,1))*zdzhs(1)*zrhwddt + zdwgdt(i,1) )
          zdwgdt(i,1)= zdwgdt(i,1) - zw_ovpv
          runoff_s(i) = runoff_s(i) + zw_ovpv*zroffdt
        ENDIF                       ! overflow of interception store

      ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
      ELSEIF (zwsnow(i) == 0.0_wp .AND.                            &
             (1._wp-ztsnow_pm(i))*zrr(i) > 0.0_wp) THEN
        zsprs  (i) = MIN(lh_f*zrr(i),(t0_melt-t_s_now(i))*zroc(i,1)*zdzhs(1)/zdt)
        ! keep freezing rain in interception storage rather than shifting it to snow
        ! zdwidt (i) = zdwidt (i) - zrr(i)
        ! zdwsndt(i) = zdwsndt(i) + zrr(i)
      END IF

      ! Influence of heatflux through snow on total forcing:
      zwsnew(i)   = zwsnow(i) + zdwsndt(i)*zdtdrhw
      IF (zwsnew(i).GT.eps_soil) THEN
        !         heat conductivity of snow as funtion of water content
        ! BR      zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnow(i)))
        !
        ! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
        !
        zalas  = 2.22_wp*EXP(1.88_wp*LOG(zrho_snow(i)/rho_i))

        ! BR 11/2005 Use alternative formulation for heat conductivity by Sun et al., 1999
        !            The water vapour transport associated conductivity is not included.

        !        zalas   = 0.023_wp+(2.290_wp-0.023_wp)* &
        !                               (7.750E-05_wp*rho_snow(i,nx) + &
        !                                1.105E-06_wp*prho_snow(i,nx)**2)

        !          zgsb(i) = zalas*(ztsnow(i) - zts(i))/zdz_snow_fl(i)
        zgsb(i) = zalas*(ztsnow(i) - zts(i))/zdz_snow_fl(i)
      END IF

      ! Calculation of the surface energy balance

      ! total forcing for uppermost soil layer
      zfor_s(i) = ( zrnet_s(i) + zshfl_s(i) + zlhfl_s(i) )       &
                     * (1._wp - zf_snow(i)) + zsprs(i) + zahf(i) &
                  + zf_snow(i) * (1._wp-ztsnow_pm(i)) * zgsb(i)

    END DO
    !$ACC END PARALLEL

    IF (msg_level >= 20) THEN
      !$ACC UPDATE HOST(soiltyp_subs, zshfl_s, zlhfl_s, zrhoch, zth_low, t) ASYNC(acc_async_queue)
      !$ACC UPDATE HOST(zts_pm, zverbo, zf_snow, qv, qv_s, tch, tcm, zts) ASYNC(acc_async_queue)
      !$ACC WAIT(acc_async_queue)
      DO i = ivstart, ivend
        IF (soiltyp_subs(i) == 1) THEN  !1=glacier and Greenland
          IF ( ABS( zshfl_s(i) )  >  500.0_wp  .OR. &
               ABS( zlhfl_s(i) )  > 2000.0_wp ) THEN
            write(*,*) 'hello sfc_terra  2: ', zshfl_s(i), zrhoch(i),zth_low(i),t(i),zts(i), &
              '  ...LHF...  ',                 zlhfl_s(i), zts_pm(i),zverbo(i),zf_snow(i),qv(i),qv_s(i), &
              '  ...CH,CM...  ', tch(i), tcm(i)
          ENDIF
        ENDIF
      END DO
    ENDIF

  ENDIF ! lmulti_snow

!------------------------------------------------------------------------------
! Section II.6: Solution of the heat conduction equation, freezing/melting
!               of soil water/ice (optionally)
!------------------------------------------------------------------------------

  ! EM: If the single-layer snow model is used, nothing changes;
  ! if the multi-layer snow model is used and zf_snow(i) == 1,
  ! then the heat conduction equation for the whole column "soil + snow" is solved
  ! (see below, after the statement IF (lmulti_snow) THEN);
  ! if the multi-layer snow model is used, but zf_snow(i) < 1._wp,
  ! then the two partial temperature updates are computed: first, the heat conduction equation
  ! for the snow-free surface is solved, second, the heat conduction equation
  ! for the whole column "soil + snow" for snow-covered surface is solved.
  ! Then, the two updates are merged together.

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    sn_frac(i) = zf_snow(i)
    IF (zwsnow(i) < eps_soil .AND. zwsnew(i) >= eps_soil) THEN
      sn_frac(i) = 0.01_wp
    ELSEIF (zwsnow(i) >= eps_soil .AND. zwsnew(i) < eps_soil) THEN
      sn_frac(i) = 0._wp
    ENDIF
  END DO

  !$ACC LOOP SEQ
  DO kso = 1, ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      t_so_free_new(i,kso) = t_so_now(i,kso)
      t_so_snow_new(i,kso) = t_so_now(i,kso)
    END DO
  END DO        ! soil layers

  !$ACC LOOP SEQ
  DO kso = 2, ke_soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zakb1, zakb2)
    DO i = ivstart, ivend
      IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._wp) THEN
        ! for heat conductivity: zalam is now 3D
        zakb1 = zalam(i,kso-1)/zroc(i,kso)
        zakb2 = zalam(i,kso  )/zroc(i,kso)
        zaga(i,kso) = -zalfa*zdt*zakb1/(zdzhs(kso)*zdzms(kso))
        zagc(i,kso) = -zalfa*zdt*zakb2/(zdzhs(kso)*zdzms(kso+1))
        zagb(i,kso) = 1._wp - zaga(i,kso) - zagc(i,kso)
        zagd(i,kso) = t_so_now(i,kso) +                                     &
               (1._wp - zalfa)*( - zaga(i,kso)/zalfa*t_so_now(i,kso-1)+ &
               (zaga(i,kso)/zalfa + zagc(i,kso)/zalfa)*t_so_now(i,kso) -  &
                zagc(i,kso)/zalfa*t_so_now(i,kso+1)  )
      END IF
    END DO
  END DO        ! soil layers


  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zakb1, zakb2)
  DO i = ivstart, ivend
!    IF (.NOT. lmulti_snow .OR. (zf_snow(i) < 1._wp .AND. zwsnew(i).GT.eps_soil)) THEN
    IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._wp) THEN
      ! for heat conductivity: zalam is now 3D: here we need layer 1
      zakb1 = hzalam(i,1)/zroc(i,1)
      zakb2 =  zalam(i,1)/zroc(i,1)
      zaga(i,1) = -zalfa*zdt*zakb1/(zdzhs(1)*zdzms(1))
      zagc(i,1) = -zalfa*zdt*zakb2/(zdzhs(1)*zdzms(2))
      zagb(i,1) = 1._wp - zaga(i,1) - zagc(i,1)
      zagd(i,1) = t_so_now(i,1) + (1._wp - zalfa)* (                 &
                      - zaga(i,1)/zalfa * t_s_now(i) +                     &
                      (zaga(i,1) + zagc(i,1))/zalfa * t_so_now(i,1) -    &
                       zagc(i,1)/zalfa * t_so_now(i,2)   )
      zaga(i,0) = 0.0_wp
      zagb(i,0) = zalfa
      zagc(i,0) = -zalfa
      ! EM: In the case of multi-layer snow model, zfor_s(i) does not include the heat conductivity flux
      ! between soil and snow (zgsb). It will be accounted for at the next semi-step, see below.
      zagd(i,0)    = zdzms(1) * zfor_s(i)/hzalam(i,1)+(1._wp-zalfa)* &
                      (t_so_now(i,1) - t_s_now(i))
      zaga(i,ke_soil+1) = 0.0_wp
      zagb(i,ke_soil+1) = 1.0_wp
      zagc(i,ke_soil+1) = 0.0_wp
      zagd(i,ke_soil+1) = t_so_now(i,ke_soil+1)
    END IF
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._wp) THEN
      zagc(i,0) = zagc(i,0)/zagb(i,0)
      zagd(i,0) = zagd(i,0)/zagb(i,0)
    END IF
  END DO

  !$ACC LOOP SEQ
  DO kso = 1, ke_soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zzz)
    DO i = ivstart, ivend
      IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._wp) THEN
        zzz = 1._wp/(zagb(i,kso) - zaga(i,kso)*zagc(i,kso-1))
        zagc(i,kso) = zagc(i,kso) * zzz
        zagd(i,kso) = (zagd(i,kso) - zaga(i,kso)*zagd(i,kso-1)) * zzz
      END IF
    END DO
  END DO                ! soil layers

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._wp) THEN
      zage(i,ke_soil+1) = (zagd(i,ke_soil+1) - zaga(i,ke_soil+1)*       &
                             zagd(i,ke_soil))/                              &
                            (zagb(i,ke_soil+1) - zaga(i,ke_soil+1)*       &
                             zagc(i,ke_soil))
      t_so_free_new(i,ke_soil+1) = zage(i,ke_soil+1) ! climate value, unchanged
    END IF
  END DO

  !$ACC LOOP SEQ
  DO kso = ke_soil,0,-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._wp) THEN
        zage(i,kso)     = zagd(i,kso) - zagc(i,kso)*zage(i,kso+1)
        ! The surface temperature computed by t_so(i,0,nnew)=zage(i,0) is
        ! presently unused
        t_so_free_new(i,kso) = zage(i,kso)
      END IF
    END DO
  END DO                ! soil layers


  IF (lmulti_snow) THEN

    ! If there is snow, the solution of the heat conduction equation
    ! goes through the whole column "soil+snow"
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (sn_frac(i) > 0._wp) THEN
        ! Uppermost snow layer, Neumann boundary condition
        zrocs(i) = (wliq_snow_now(i,1)/wtot_snow_now(i,1)*chc_w + &
              (wtot_snow_now(i,1)-wliq_snow_now(i,1))/wtot_snow_now(i,1)*chc_i)*rho_snow_mult_now(i,1)
        zakb      = zalas_mult(i,1)/zrocs(i)
        zaga(i,1) = -zalfa*zdt*zakb/(zdzh_snow(i,1)*zdzm_snow(i,1))
        zagc(i,1) = -zalfa*zdt*zakb/(zdzh_snow(i,1)*zdzm_snow(i,2))
        zagb(i,1) = 1._wp - zaga(i,1) - zagc(i,1)
        zagd(i,1) = ztsnow_mult(i,1) + (1._wp - zalfa)* (-zaga(i,1)/zalfa * t_s_now(i) + &
                    (zaga(i,1) + zagc(i,1))/zalfa * t_so_now(i,1) - zagc(i,1)/zalfa * t_so_now(i,2) )
        zaga(i,0) = 0.0_wp
        zagb(i,0) = zalfa
        zagc(i,0) = -zalfa
        zagd(i,0) = zdzm_snow(i,1) * zfor_snow_mult(i)/zalas_mult(i,1)+(1._wp-zalfa)* &
                   (t_so_now(i,1) - t_s_now(i))
        zagc(i,0) = zagc(i,0)/zagb(i,0)
        zagd(i,0) = zagd(i,0)/zagb(i,0)

        ! Lowermost soil layer, Dirichlet boundary condition
        zaga(i,ke_snow + ke_soil+1) = 0.0_wp
        zagb(i,ke_snow + ke_soil+1) = 1.0_wp
        zagc(i,ke_snow + ke_soil+1) = 0.0_wp
        zagd(i,ke_snow + ke_soil+1) = t_so_now(i,ke_soil+1)

        ! Lowermost snow layer, special treatment
        zrocs(i) = (wliq_snow_now(i,ke_snow)/wtot_snow_now(i,ke_snow)*chc_w +  &
          & (wtot_snow_now(i,ke_snow)-wliq_snow_now(i,ke_snow))/                    &
          & wtot_snow_now(i,ke_snow)*chc_i)*rho_snow_mult_now(i,ke_snow)
        zakb = zalas_mult(i,ke_snow)/zrocs(i)
        zaga(i,ke_snow) = -zalfa*zdt*zakb/(zdzh_snow(i,ke_snow)*zdzm_snow(i,ke_snow))
        zakb = (zalas_mult(i,ke_snow)/zrocs(i)*(-zhm_snow(i,ke_snow))+hzalam(i,1)/ &
                zroc(i,1)*zmls(1))/(zmls(1)-zhm_snow(i,ke_snow))
        zagc(i,ke_snow) = -zalfa*zdt*zakb/(zdzh_snow(i,ke_snow)*(zmls(1)-zhm_snow(i,ke_snow)))
        zagb(i,ke_snow) = 1.0_wp - zaga(i,ke_snow) - zagc(i,ke_snow)
        zagd(i,ke_snow) = ztsnow_mult(i,ke_snow) + &
          &     (1._wp - zalfa)*( - zaga(i,ke_snow)/zalfa*ztsnow_mult(i,ke_snow-1) + &
          &     (zaga(i,ke_snow)/zalfa + zagc(i,ke_snow)/zalfa)*ztsnow_mult(i,ke_snow) - &
          &     zagc(i,ke_snow)/zalfa*t_so_now(i,1)  )

        ! Uppermost soil layer, special treatment
        zakb = (zalas_mult(i,ke_snow)/zrocs(i)*(-zhm_snow(i,ke_snow))+     &
                hzalam(i,1)/zroc(i,1)*zmls(1))/(zmls(1)-zhm_snow(i,ke_snow))
        zaga(i,ke_snow+1) = -zalfa*zdt*zakb/(zdzhs(1)*(zmls(1)-zhm_snow(i,ke_snow)))
        zakb = hzalam(i,1)/zroc(i,1)
        zagc(i,ke_snow+1) = -zalfa*zdt*zakb/(zdzhs(1)*zdzms(2))
        zagb(i,ke_snow+1) = 1._wp - zaga(i,ke_snow+1) - zagc(i,ke_snow+1)
        zagd(i,ke_snow+1) = t_so_now(i,1) + &
          &    (1._wp - zalfa)*( - zaga(i,ke_snow+1)/zalfa*ztsnow_mult(i,ke_snow) + &
          &    (zaga(i,ke_snow+1)/zalfa + zagc(i,ke_snow+1)/zalfa)*t_so_now(i,1)      - &
          &    zagc(i,ke_snow+1)/zalfa*t_so_now(i,2)  )
      END IF
    END DO

    ! Snow layers
    DO ksn = 2, ke_snow-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (sn_frac(i) > 0._wp) THEN
          zrocs(i) = (wliq_snow_now(i,ksn)/wtot_snow_now(i,ksn)*chc_w + &
            & (wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn))/ &
            & wtot_snow_now(i,ksn)*chc_i)*rho_snow_mult_now(i,ksn)
          zakb = zalas_mult(i,ksn)/zrocs(i)
          zaga(i,ksn) = -zalfa*zdt*zakb/(zdzh_snow(i,ksn)*zdzm_snow(i,ksn))
          zakb = (zalas_mult(i,ksn)/zrocs(i)*zdzm_snow(i,ksn)+zalas_mult(i,ksn+1)/zrocs(i)*zdzm_snow(i,ksn+1))/ &
                 (zdzm_snow(i,ksn)+zdzm_snow(i,ksn+1))
          zagc(i,ksn) = -zalfa*zdt*zakb/(zdzh_snow(i,ksn)*zdzm_snow(i,ksn+1))
          zagb(i,ksn) = 1._wp - zaga(i,ksn) - zagc(i,ksn)
          zagd(i,ksn) = ztsnow_mult(i,ksn) + &
            &     (1._wp - zalfa)*( - zaga(i,ksn)/zalfa*ztsnow_mult(i,ksn-1) + &
            &     (zaga(i,ksn)/zalfa + zagc(i,ksn)/zalfa)*ztsnow_mult(i,ksn)     - &
            &     zagc(i,ksn)/zalfa*ztsnow_mult(i,ksn+1)  )
        END IF
      END DO
    END DO                ! snow layers

    ! Soil layers
    DO ksn = ke_snow+2, ke_snow + ke_soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (sn_frac(i) > 0._wp) THEN
          kso = ksn - ke_snow
          zakb1 = zalam(i,kso-1)/zroc(i,kso)
          zakb2 = zalam(i,kso  )/zroc(i,kso)
          zaga(i,ksn) = -zalfa*zdt*zakb1/(zdzhs(kso)*zdzms(kso))
          zagc(i,ksn) = -zalfa*zdt*zakb2/(zdzhs(kso)*zdzms(kso+1))
          zagb(i,ksn) = 1._wp - zaga(i,ksn) - zagc(i,ksn)
          zagd(i,ksn) = t_so_now(i,kso) + &
            &    (1._wp - zalfa)*( - zaga(i,ksn)/zalfa*t_so_now(i,kso-1) +  &
            &    (zaga(i,ksn)/zalfa + zagc(i,ksn)/zalfa)*t_so_now(i,kso)     -  &
            &    zagc(i,ksn)/zalfa*t_so_now(i,kso+1)  )
        END IF
      END DO
    END DO                ! soil layers


    DO kso = 1, ke_snow + ke_soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (sn_frac(i) > 0._wp) THEN
          zzz = 1._wp/(zagb(i,kso) - zaga(i,kso)*zagc(i,kso-1))
          zagc(i,kso) = zagc(i,kso) * zzz
          zagd(i,kso) = (zagd(i,kso) - zaga(i,kso)*zagd(i,kso-1)) * zzz
        END IF
      END DO
    END DO                ! snow + soil layers

    ! Back substitution, lowermost soil layer
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (sn_frac(i) > 0._wp) THEN
        zage(i,ke_snow+ke_soil+1) =  &
          &    (zagd(i,ke_snow+ke_soil+1) - zaga(i,ke_snow+ke_soil+1) * zagd(i,ke_snow+ke_soil))/ &
          &    (zagb(i,ke_snow+ke_soil+1) - zaga(i,ke_snow+ke_soil+1) * zagc(i,ke_snow+ke_soil))
        t_so_snow_new(i,ke_soil+1) = zage(i,ke_snow+ke_soil+1) ! climate value, unchanged
      END IF
    END DO

    ! Back substitution, soil layers
    DO kso = ke_snow+ke_soil, ke_snow+1, -1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (sn_frac(i) > 0._wp) THEN
          zage(i,kso) = zagd(i,kso) - zagc(i,kso)*zage(i,kso+1)
          t_so_snow_new(i,kso-ke_snow) = zage(i,kso)
        END IF
      END DO
    END DO                ! soil layers

    ! Back substitution, snow layers
    DO ksn = ke_snow,1,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (sn_frac(i) > 0._wp) THEN
          zage(i,ksn) = zagd(i,ksn) - zagc(i,ksn)*zage(i,ksn+1)
          ztsnown_mult(i,ksn) = zage(i,ksn)
        END IF
      END DO
    END DO                ! snow layers

!in case of thin snowpack (less than zswitch), apply single-layer snow model
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (sn_frac(i) > 0._wp) THEN
        zgsb(i) = ((zalas_mult(i,ke_snow)*(-zhm_snow(i,ke_snow))+hzalam(i,1)*zdzms(1))/ &
                  (-zhm_snow(i,ke_snow)+zdzms(1)) * &
                  (ztsnown_mult(i,ke_snow) - t_so_now(i,1))/(-zhm_snow(i,ke_snow) &
                  +zdzms(1)))*zf_snow(i)
        zrocs(i) = (wliq_snow_now(i,1)/wtot_snow_now(i,1)*chc_w + &
          (wtot_snow_now(i,1)-wliq_snow_now(i,1))/wtot_snow_now(i,1)*chc_i)*rho_snow_mult_now(i,1)
        zswitch(i) = (-zfor_snow_mult(i)+zgsb(i))/50./zrocs(i)*zdt*ke_snow
        zswitch(i) = MAX(zswitch(i),1.E-03_wp)
      ELSE
        zswitch(i) = 0.0_wp
      END IF
    END DO


    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF(zwsnew(i) .LT. zswitch(i) .AND. sn_frac(i) > 0._wp) THEN
        ztsn  (i) = t_so_now(i,1)
        tmp_num(i) = ztsnow_mult(i,1) + zdt*2._wp*(zfor_snow_mult(i) - zgsb(i))  &
                           /zrocs(i)/(zswitch(i)/rho_snow_mult_now(i,1)*rho_w) &
                           &- ( ztsn(i) - zts(i) )
        zalas  = 2.22_wp*EXP(1.88_wp*LOG(rho_snow_mult_now(i,1)/rho_i))

        ztsnow_im    = - zrhoch(i) * (cp_d + zdqvtsnow(i) * lh_s) - zalas/h_snow_now(i)
        zfak  = MAX(eps_soil,1.0_wp - zdt*zalfa*ztsnow_im/zrocs(i)/h_snow_now(i))
        tmp_num(i) = ztsnow_mult(i,1) + (tmp_num(i)-ztsnow_mult(i,1))/zfak
      END IF
    END DO

!    DO ksn = 1, ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF(sn_frac(i) > 0._wp) THEN
!          IF(zwsnew(i) .LT. zswitch(i)) THEN
!            ztsnown_mult(i,ksn) = tmp_num(i)
!          END IF
          ztsnown_mult(i,0) = ztsnown_mult(i,1)
        END IF
      END DO
!    END DO


    DO ksn = 1, ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF(zwsnew(i) .GT. eps_soil .and. zwsnew(i) .LT. zswitch(i)) THEN

          IF((zfor_snow_mult(i)-zgsb(i))*zdt > zwsnew(i)*rho_w*lh_f) THEN
            ztsnown_mult(i,ksn) = t_so_now(i,0)
          ELSE IF(zfor_snow_mult(i)-zgsb(i) .GT. 0._wp) THEN
            ztsnown_mult(i,ksn) = ztsnow_mult(i,ksn) + &
              (zfor_snow_mult(i)-zgsb(i))*zdt/(chc_i*wtot_snow_now(i,ksn))/rho_w/ke_snow
          END IF
        END IF
      END DO
    END DO


    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF(zwsnew(i) .GT. eps_soil .and. zwsnew(i) .LT. zswitch(i)) THEN

        IF((zfor_snow_mult(i)-zgsb(i))*zdt > zwsnew(i)*rho_w*lh_f) THEN
          zdwsndt(i) = zdwsndt(i) - zwsnew(i)*rho_w/zdt
          zwsnew(i)  = 0._wp
          ztsnown_mult(i,0) = t_so_now(i,0)
        ELSE IF((zfor_snow_mult(i)-zgsb(i)) .GT. 0._wp) THEN
          ztsnown_mult(i,0) = ztsnown_mult(i,1)
        END IF

      END IF
    END DO

  END IF    ! lmulti_snow

  ! Combining the two partial updates of soil temperature
  IF(lmulti_snow) THEN
    DO kso = 1, ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        t_so_new(i,kso) = t_so_snow_new(i,kso)*sn_frac(i) + t_so_free_new(i,kso)*(1._wp - sn_frac(i))
      END DO
    END DO
  ELSE
    !$ACC LOOP SEQ
    DO kso = 1, ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        t_so_new(i,kso) = t_so_free_new(i,kso)
      END DO
    END DO
  END IF

! IF(lmelt) THEN ! + lmelt_var
    !$ACC LOOP SEQ
    DO kso = 1,ke_soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ztx, zxx, zliquid, znen, zfak, zdelwice) &
      !$ACC   PRIVATE(zwso_new, zargu)
      DO i = ivstart, ivend
        IF (soiltyp_subs(i) >= 3) THEN
          ztx      = t0_melt
          zw_m(i)  = zporv(i,kso)*zdzhs(kso)

          IF(t_so_new(i,kso).LT.(t0_melt-eps_temp)) THEN
!US         zxx    = g*zpsis(i)/lh_f
!US         zw_m(i) = zw_m(i)*EXP(-zedb(i)*LOG((t_so_new(i,kso) - t0_melt)/(t_so_new(i,kso)*zxx)) )

            IF (t_so_new(i,kso) < t_zw_low) THEN
              zw_m(i) = zw_m_low(i,kso)
            ELSE IF (t_so_new(i,kso) < t_zw_up) THEN ! Logarithmic Interpolation between -3 degC and -40 degC 
              zw_m(i) = zw_m_low(i,kso)*EXP((t_so_new(i,kso) - t_zw_low)*                          &
                (LOG(zporv(i,kso)*zdzhs(kso)*zw_m_up(i)) - LOG(zw_m_low(i,kso)))/(t_zw_up-t_zw_low))
            ELSE
              zw_m(i) = zw_m(i)*EXP(-zedb(i)*LOG((t_so_new(i,kso) - t0_melt)/(t_so_new(i,kso)*zaa(i))) )
            END IF

            zliquid= MAX(eps_div,w_so_now(i,kso) -  w_so_ice_now(i,kso))
            znen   = 1._wp-zaa(i)*EXP(zb_por(i)*LOG(zporv(i,kso)*zdzhs(kso)/zliquid))
            ztx    = t0_melt/znen
          ENDIF
          ztx      = MIN(t0_melt,ztx)
          zfak     = zroc(i,kso)*zdzhs(kso)/(lh_f*rho_w)
          zdelwice = - zfak*(t_so_new(i,kso)-ztx)
          zwso_new = w_so_now(i,kso) + zdt*zdwgdt(i,kso)/rho_w
          zargu = zwso_new - zw_m(i) - w_so_ice_now(i,kso)
          IF (t_so_new(i,kso) > t0_melt .AND. w_so_ice_now(i,kso) > 0.0_wp) THEN
            ! melting point adjustment (time scale 30 min)
            zdelwice = - MIN(w_so_ice_now(i,kso), zdwi_scal*(t_so_new(i,kso)-t0_melt)*zfak)
          ELSE IF (zdelwice < 0.0_wp) THEN ! this branch contains cases of melting and freezing
            zdelwice = - MIN( - zdelwice,MIN(-zargu,w_so_ice_now(i,kso)))
            ! limit latent heat consumption due to melting to half the temperature increase since last time step
            ! or 2.5 K within 30 min; the freezing rate is limited below
            zdelwice = - MIN( - zdelwice,MAX(2.5_wp*zdwi_scal,0.5_wp*(t_so_new(i,kso)-t_so_now(i,kso)))*zfak)
          ELSE
            zdelwice = MIN(zdelwice,MAX(zargu,0.0_wp))
          ENDIF
          IF (zdelwice > 0.0_wp) THEN
            ! limit latent heat release due to freezing to half the differene from the melting point
            zdelwice = MIN(zdelwice,0.5_wp*MAX(0.0_wp,(t0_melt-t_so_new(i,kso)))*zfak)
          ENDIF
          w_so_ice_new(i,kso) = w_so_ice_now(i,kso) + zdelwice
          t_so_new    (i,kso) = t_so_new    (i,kso) + zdelwice/zfak
        END IF
      END DO
    ENDDO
! END IF ! lmelt

!------------------------------------------------------------------------------
! Section II.7: Energy budget and temperature prediction at snow-surface
!------------------------------------------------------------------------------
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zrnet_snow, zfor_snow, zalas, ztsnow_im, zfak)
  DO i = ivstart, ivend
    ! next line has to be changed if a soil surface temperature is
    ! predicted by the heat conduction equation
    zdtsdt (i) = (t_so_new(i,1) - zts(i))*z1d2dt
    ztsn   (i) =  t_so_new(i,1)
    IF(.NOT. lmulti_snow)  &
      ztsnown(i) = ztsn(i)       ! default setting
    zwsnew(i)       = zwsnow(i) + zdwsndt(i)*zdtdrhw

    ! forcing contributions for snow formation of dew and rime are
    ! contained in ze_ges, heat fluxes must not be multiplied by
    ! snow covered fraction

    IF(.NOT. lmulti_snow) THEN

      zrnet_snow    = sobs(i) + zthsnw(i)
      zshfl_snow(i) = zrhoch(i)*cp_d*(zth_low(i) - ztsnow(i))
      zlhfl_snow(i) = lh_s*zversn(i)
      zqhfl_snow(i) = zversn(i)
      zfor_snow     = zrnet_snow + zshfl_snow(i) + zlhfl_snow(i)

      ! forecast of snow temperature Tsnow
      IF (ztsnow(i) < t0_melt .AND. zwsnew(i) > eps_soil) THEN
        ztsnown(i) = ztsnow(i) + zdt*2._wp*(zfor_snow - zgsb(i))  &
                       /zrocs(i) - ( ztsn(i) - zts(i) )

        ! implicit formulation
        ! BR        zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnew(i)))
        ! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
        !
        zalas  = 2.22_wp*EXP(1.88_wp*LOG(zrho_snow(i)/rho_i))

        ztsnow_im    = - zrhoch(i) * (cp_d + zdqvtsnow(i) * lh_s)       &
                                         - zalas/zdz_snow_fl(i)
        zfak  = MAX(eps_div,1.0_wp - zdt*zalfa*ztsnow_im/zrocs(i))
        ztsnown(i) = ztsnow(i) + (ztsnown(i)-ztsnow(i))/zfak
      END IF

      zdtsnowdt(i) = (ztsnown(i) - ztsnow(i))*z1d2dt
    ENDIF
  END DO
  !$ACC END PARALLEL

  IF (msg_level >= 20) THEN
    !$ACC UPDATE HOST(zshfl_s, zlhfl_s, zf_snow, zrhoch, t, zts, zts_pm, zverbo, qv, qv_s) ASYNC(acc_async_queue)
    !$ACC UPDATE HOST(zshfl_snow, zlhfl_snow, zth_low, ztsnow) ASYNC(acc_async_queue)
    !$ACC UPDATE HOST(zwsnow, zrr, zrs, zdwsndt) ASYNC(acc_async_queue)
    !$ACC WAIT(acc_async_queue)
    DO i = ivstart, ivend
!     IF (soiltyp_subs(i) == 1) THEN  !1=glacier and Greenland
      IF ( ABS( zshfl_s(i)    * (1.0_wp-zf_snow(i)) )  >  700.0_wp .OR. &
           ABS( zlhfl_s(i)    * (1.0_wp-zf_snow(i)) )  > 2000.0_wp        ) THEN
        write(*,*) 'soil soil: ', zshfl_s(i), zrhoch(i),zth_low(i),t(i),zts(i), &
          '  ...LHF...  ',        zlhfl_s(i), zts_pm(i),zverbo(i),zf_snow(i),qv(i),qv_s(i), &
          '  ...CH,CM...  ', tch(i), tcm(i)
      ENDIF
      IF ( ABS( zshfl_snow(i) )  >  500.0_wp  .OR. &
           ABS( zlhfl_snow(i) )  > 2000.0_wp ) THEN
        write(*,*) 'soil: ', zshfl_snow(i), zlhfl_snow(i), '....', &
        zth_low(i), ztsnow(i), '....', &
        zwsnow(i), zrr(i), zrs(i), zdwsndt(i)
      ENDIF
!     ENDIF
    END DO
  ENDIF

  IF (lmulti_snow) THEN
    DO ksn = 0, ke_snow
      DO i = ivstart, ivend
        IF (zwsnew(i) > eps_soil) THEN
          zdtsnowdt_mult(i,ksn) = (ztsnown_mult(i,ksn) - ztsnow_mult(i,ksn))*z1d2dt
        END IF
      END DO
    ENDDO
  ENDIF


!------------------------------------------------------------------------------
! Section II.8: Melting of snow ,infiltration and surface runoff of snow water
!------------------------------------------------------------------------------

! If the soil surface temperature predicted by the equation of heat conduction
! is used instead of using T_s = T_so(1), the following section has to be
! adjusted accordingly.
!
! Basically this snow model uses heat fluxes to either heat the uppermost soil
! layer, if the snow surface temperature exceeds t0_melt, or to heat the snow, if
! the temperature of the uppermost soil layer exceeds t0_melt. Melting is considered
! after this process, if the snow temperature equals t0_melt AND the temperature of
! the uppermost soil layer exceeds t0_melt. The excess heat (t_so(1)-t0_melt) is used for
! melting snow. In cases this melting may be postponed to the next time step.

  IF (.NOT. lmulti_snow) THEN

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP GANG VECTOR PRIVATE(ztsnownew, ze_avail, ze_total, zfr_melt) &
    !$ACC   PRIVATE(ztsnew, zdelt_s, zdwgme, zro, zredfu, zw_ovpv, zfr_ice_free)
    DO i = ivstart, ivend

      zwsnn  (i)  = zwsnow(i) + zdtdrhw*zdwsndt(i)
      zwsnew (i)  = zwsnn(i)
      ztsnownew     = ztsnown(i)
      ze_avail      = 0.0_wp
      ze_total      = 0.0_wp
      zfr_melt      = 0.0_wp

      IF (zwsnew(i) > eps_soil) THEN        ! points with snow cover only
        ! first case: T_snow > t0_melt: melting from above
        ! ----------
        IF (ztsnown(i) > t0_melt .AND. t_so_new(i,1) < t0_melt ) THEN
          ! Limit max w_snow in melting conditions for consistency with heat capacity calculation
          zdwsnm(i)    = MIN(1.5_wp*rho_snow_now(i)/rho_w,zwsnew(i))* &
                .5_wp*(ztsnown(i) - (t0_melt - eps_temp))/ &
               (.5_wp* (zts(i) - (t0_melt - eps_temp)) - lh_f/chc_i)
          zdwsnm(i)    = zdwsnm(i)*z1d2dt*rho_w
          zdwsndt(i)   = zdwsndt (i) + zdwsnm(i)
          meltrate(i)  = - zdwsnm(i)
          ztsnownew    = t0_melt - eps_temp
          zdtsnowdt(i) = zdtsnowdt(i) + (ztsnownew - ztsnown(i))*z1d2dt
          ! decide which parts of the meltwater are passed to w_so and runoff, respectively
          zfr_ice_free = 1._wp-ziw_fr(i,1)/zporv(i,1)
          zdwgme       = zfr_ice_free*meltrate(i)*zrock(i)       ! contribution to w_so
          zredfu       = MAX( 0.0_wp,  MIN( 1.0_wp, (zw_fr(i,1) -  &
                         zfcap(i,1))/MAX(zporv(i,1)-zfcap(i,1), eps_div)))
          zdwgdt(i,1)  = zdwgdt(i,1) + zdwgme*(1._wp - zredfu)
          zro          = meltrate(i) - zdwgme*(1._wp - zredfu)

          ! zro-, zdw_so_dt-correction in case of pore volume overshooting
          zw_ovpv = MAX(0._wp, (zw_fr(i,1)-zporv(i,1))*zdzhs(1)*zrhwddt + zdwgdt(i,1) )
          zro = zro + zw_ovpv
          zdwgdt(i,1)= zdwgdt(i,1) - zw_ovpv
          runoff_s (i)= runoff_s(i) + zro*zroffdt
        ENDIF ! melting from above

        IF (t_so_new(i,1) >= t0_melt) THEN
          !second case:  temperature of uppermost soil layer > t0_melt. First a
          !-----------   heat redistribution is performed. As a second step,
          !              melting of snow is considered.
          ! a) Heat redistribution
          ztsnew = t0_melt + eps_temp
          ztsnownew      = ztsnown(i) + zf_snow_lim(i)*(ztsn(i) - ztsnew) +  &
               2._wp*(ztsn(i) - ztsnew)*zroc(i,1)*zdzhs(1)/zrocs(i)
          zdtsdt(i)    = zdtsdt(i) + zf_snow_lim(i)*(ztsnew - t_so_new(i,1))*z1d2dt
          zdtsnowdt(i) = zdtsnowdt(i) + (ztsnownew - ztsnown(i))*z1d2dt
          ! b) Melting of snow (if possible)
          IF (ztsnownew > t0_melt) THEN
            ze_avail     = 0.5_wp*(ztsnownew - t0_melt)*zrocs(i)*zf_snow_lim(i)
            ze_total     = lh_f*zwsnew(i)*rho_w
            zfr_melt     = MIN(1.0_wp,ze_avail/ze_total)
            zdtsnowdt(i)= zdtsnowdt(i) + (t0_melt - ztsnownew)*z1d2dt
            zdelt_s      = MAX(0.0_wp,(ze_avail - ze_total)/(zroc(i,1)* &
                                                                        zdzhs(1)))
            zdtsdt(i)  = zdtsdt(i) + zf_snow_lim(i)*zdelt_s*z1d2dt

            ! melted snow is allowed to penetrate the soil (up to field
            ! capacity), if the soil type is neither ice nor rock (zrock = 0);
            ! else it contributes to surface run-off;
            ! fractional water content of the first soil layer determines
            ! a reduction factor which controls additional run-off
            zdwsnm(i)   = zfr_melt*zwsnew(i)*z1d2dt*rho_w  ! available water
            zdwsndt(i)  = zdwsndt (i) - zdwsnm(i)
            meltrate(i) = meltrate(i) + zdwsnm(i)
            zdwgme        = zdwsnm(i)*zrock(i)             ! contribution to w_so
            zro           = (1._wp - zrock(i))*zdwsnm(i)      ! surface runoff
            zredfu        = MAX( 0.0_wp,  MIN( 1.0_wp, (zw_fr(i,1) -  &
                            zfcap(i,1))/MAX(zporv(i,1)-zfcap(i,1), eps_div)))
            zdwgdt(i,1) = zdwgdt(i,1) + zdwgme*(1._wp - zredfu)
            zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                                   ! for this fraction

            ! zro-, zdw_so_dt-correction in case of pore volume overshooting
            zw_ovpv = MAX(0._wp, zw_fr(i,1)* zdzhs(1) * zrhwddt +  &
                       zdwgdt(i,1) - zporv(i,1) * zdzhs(1) * zrhwddt)
            zro = zro + zw_ovpv
            zdwgdt(i,1)= zdwgdt(i,1) - zw_ovpv


            IF (zfr_melt > 0.9999_wp) zdwsndt(i)= -zwsnow(i)*zrhwddt
            runoff_s (i)= runoff_s(i) + zro*zroffdt

          END IF   ! snow melting
        END IF     ! snow and/or soil temperatures
      END IF       ! points with snow cover only
    END DO
    !$ACC END PARALLEL

  ELSE       ! multi-layer snow scheme

    DO i = ivstart, ivend

      zwsnew(i) = zwsnow(i) + zdtdrhw*zdwsndt(i)
      zdwsnm(i) = 0.0_wp

      ze_out  (i) = 0.0_wp
      zqbase  (i) = 0.0_wp
      zcounter(i) = 0.0_wp

      ze_rad(i) = 0.0_wp
      IF(zextinct(i,1).gt.0.0_wp) ze_rad(i) = zf_snow(i) * sobs(i)

      ztsnownew_mult(i,0) = ztsnown_mult(i,0)
    END DO

    DO ksn = 1,ke_snow
      DO i = ivstart, ivend
        IF (zwsnew(i) > eps_soil) THEN        ! points with snow cover only

          zrefr(i) = 0.0_wp
          zmelt(i) = 0.0_wp
          ztsnownew_mult(i,ksn) = ztsnown_mult(i,ksn)

          IF(zdzh_snow(i,ksn) - wliq_snow_now(i,ksn).GT.eps_soil .OR. &
            wtot_snow_now(i,ksn) - wliq_snow_now(i,ksn).GT.eps_soil) THEN
            zrho_dry_old(i) = MAX(wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn), &
              &                     eps_soil)                                               &
              &                 *rho_w/(zdzh_snow(i,ksn) - wliq_snow_now(i,ksn))
          ELSE
            zrho_dry_old(i) = rho_w
          END IF

          ztsnownew_mult(i,ksn) = (ztsnown_mult(i,ksn)*wtot_snow_now(i,ksn) &
            &                       + t0_melt*zqbase(i)*zdt)/(zqbase(i)*zdt     &
            &                       + wtot_snow_now(i,ksn))

          IF(zextinct(i,ksn).eq.0.0_wp) THEN
            ze_in = ze_out(i)
          ELSE
            IF(ksn.eq.ke_snow) THEN     ! all the rest of radiation is absorbed by the lowermost snow layer
              ze_in = ze_out(i) + ze_rad(i)
            ELSE
              zcounter(i) = EXP (-zextinct(i,ksn)*zdzh_snow(i,ksn))
              ze_in = ze_out(i) + ze_rad(i) * (1._wp - zcounter(i))
              ze_rad(i) = ze_rad(i) * zcounter(i)
            END IF
          END IF

          ztsnownew_mult(i,ksn) = ztsnownew_mult(i,ksn) &
            &                       + ze_in*zdt/(chc_i*wtot_snow_now(i,ksn))/rho_w
          wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn) + zqbase(i)*zdt
          wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zqbase(i)*zdt

          zdzh_old = zdzh_snow(i,ksn)
          zdzh_snow(i,ksn) = zdzh_snow(i,ksn) + zqbase(i)*zdt

          rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*&
            &                                rho_w/zdzh_snow(i,ksn), &
            &                                0.0_wp)

          IF(ztsnownew_mult(i,ksn) .GT. t0_melt) THEN

            IF(wtot_snow_now(i,ksn) .LE. wliq_snow_now(i,ksn)) THEN
              ze_out(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
                &      *z1d2dt*rho_w
              zmelt(i) = 0.0_wp
            ELSEIF(chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt)/lh_f <= &
              wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn)) THEN
              zmelt(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
                &          *z1d2dt/lh_f
              ze_out(i) = 0.0_wp
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zmelt(i)*zdt
            ELSE
              zmelt(i) = (wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn))*z1d2dt
              ze_out(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
                &      *z1d2dt*rho_w - zmelt(i)*lh_f*rho_w
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zmelt(i)*zdt
            END IF
            ztsnownew_mult(i,ksn) = t0_melt

          ELSE
            ! T<0
            IF(wliq_snow_now(i,ksn) .GT. -chc_i*wtot_snow_now(i,ksn) &
              & *(ztsnownew_mult(i,ksn) - t0_melt)/lh_f) THEN
              zrefr(i) = -chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn) &
                &          - t0_melt)*z1d2dt/lh_f
              ztsnownew_mult(i,ksn)   = t0_melt
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) - zrefr(i)*zdt
            ELSE
              zrefr(i) = wliq_snow_now(i,ksn)*z1d2dt
              wliq_snow_now(i,ksn) = 0.0_wp
              ztsnownew_mult(i,ksn)   = ztsnownew_mult(i,ksn) + zrefr(i)*zdt*lh_f &
                &                         /(chc_i*wtot_snow_now(i,ksn))
            END IF
            ze_out(i) = 0.0_wp

          END IF

          zdtsnowdt_mult(i,ksn) = zdtsnowdt_mult(i,ksn) + &
                                    (ztsnownew_mult(i,ksn) - ztsnown_mult(i,ksn))*z1d2dt
          IF(wtot_snow_now(i,ksn) .LE. wliq_snow_now(i,ksn)) THEN
            zqbase(i)           = wliq_snow_now(i,ksn)*z1d2dt
            wliq_snow_now(i,ksn) = 0.0_wp
            wtot_snow_now(i,ksn) = 0.0_wp
            zdzh_snow(i,ksn)    = 0.0_wp
            rho_snow_mult_now(i,ksn)  = 0.0_wp
          ELSE
            IF(zrefr(i) .GT. 0.0_wp .OR. zmelt(i) .GT. 0.0_wp) THEN
              zadd_dz = 0.0_wp
              zadd_dz = MAX(zrefr(i),0._wp)*(-1.0_wp + 1.0_wp/rho_i*rho_w)*zdt
              zadd_dz = MAX(zmelt(i),0._wp)*(-1.0_wp/zrho_dry_old(i)*rho_w &
                &       + 1.0_wp)*zdt
              zdzh_snow(i,ksn)   = zdzh_snow(i,ksn) + zadd_dz
              rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*rho_w &
                &                            /zdzh_snow(i,ksn),0.0_wp)
              IF(wtot_snow_now(i,ksn) .LE. 0.0_wp) zdzh_snow(i,ksn) = 0.0_wp
              IF(rho_snow_mult_now(i,ksn) .GT. rho_w) THEN
                zdzh_snow(i,ksn)   = zdzh_snow(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
                rho_snow_mult_now(i,ksn) = rho_w
              END IF
            END IF

            zsn_porosity = 1._wp - (rho_snow_mult_now(i,ksn)/rho_w -  &
                           wliq_snow_now(i,ksn)/zdzh_snow(i,ksn))/rho_i*rho_w - &
                           wliq_snow_now(i,ksn)/zdzh_snow(i,ksn)
            zsn_porosity = MAX(zsn_porosity,cwhc + 0.1_wp)
            zp1 = zsn_porosity - cwhc

            IF (wliq_snow_now(i,ksn)/zdzh_snow(i,ksn) .GT. cwhc) THEN
              zfukt             = (wliq_snow_now(i,ksn)/zdzh_snow(i,ksn) - cwhc)/zp1
              zq0               = chcond * zfukt**3
              zqbase(i)       = MIN(zq0*zdt,wliq_snow_now(i,ksn))
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) - zqbase(i)
              wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn) - zqbase(i)

              zdzh_old = zdzh_snow(i,ksn)
              zdzh_snow(i,ksn) = zdzh_snow(i,ksn) - zqbase(i)
              zqbase(i)        = zqbase(i)*z1d2dt

              IF(zdzh_snow(i,ksn) .LT. eps_soil*0.01_wp) THEN
                wliq_snow_now(i,ksn) = 0.0_wp
                wtot_snow_now(i,ksn) = 0.0_wp
                zdzh_snow(i,ksn)     = 0.0_wp
                rho_snow_mult_now(i,ksn)  = 0.0_wp
              ELSE
                rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*rho_w &
                  &                            /zdzh_snow(i,ksn),0.0_wp)
                IF(wtot_snow_now(i,ksn) .LE. 0.0_wp) zdzh_snow(i,ksn) = 0.0_wp
                IF(rho_snow_mult_now(i,ksn) .GT. rho_w) THEN
                  zdzh_snow(i,ksn)   = zdzh_snow(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
                  rho_snow_mult_now(i,ksn) = rho_w
                END IF
              END IF
            ELSE
              zqbase(i) = 0.0_wp
            END IF
          END IF

        END IF       ! points with snow cover only
      END DO
    END DO        ! snow layers


    DO i = ivstart, ivend
      IF (zwsnew(i) > eps_soil) THEN        ! points with snow cover only
        zdwsnm(i) = zqbase(i)*rho_w       ! ksn == ke_snow
      END IF       ! points with snow cover only
    END DO

    DO i = ivstart, ivend
      IF (zwsnew(i) > eps_soil) THEN        ! points with snow cover only
        zdwsndt(i)  = zdwsndt(i) - zdwsnm(i)

        ! melted snow is allowed to penetrate the soil (up to field
        ! capacity), if the soil type is neither ice nor rock (zrock = 0);
        ! else it contributes to surface run-off;
        ! fractional water content of the first soil layer determines
        ! a reduction factor which controls additional run-off

        zdwgme        = zdwsnm(i)*zrock(i)             ! contribution to w_so
        zro           = (1._wp - zrock(i))*zdwsnm(i)      ! surface runoff
        zredfu        = MAX( 0.0_wp,  MIN( 1.0_wp, (zw_fr(i,1) -  &
                        zfcap(i,1))/MAX(zporv(i,1)-zfcap(i,1), eps_soil)))
        zdwgdt(i,1) = zdwgdt(i,1) + zdwgme*(1._wp - zredfu)
        zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                               ! for this fraction

        ! zro-, zdw_so_dt-correction in case of pore volume overshooting
        zw_ovpv = MAX(0._wp, zw_fr(i,1)* zdzhs(1) * zrhwddt +  &
                  zdwgdt(i,1) - zporv(i,1) * zdzhs(1) * zrhwddt)
        zro = zro + zw_ovpv
        zdwgdt(i,1)= zdwgdt(i,1) - zw_ovpv

        runoff_s(i) = runoff_s(i) + zro*zroffdt
      END IF       ! points with snow cover only
    END DO

    ! snow densification due to gravity and metamorphism
    DO ksn = 2, ke_snow
      zp(:,ksn) = 0.0_wp                         ! gravity, Pa
      DO k = ksn,1,-1
        DO i = ivstart, ivend
           zp(i,ksn) = zp(i,ksn) + rho_snow_mult_now(i,k)*g*zdzh_snow(i,ksn)
        END DO
      END DO
    END DO

    DO ksn = 2, ke_snow
      DO i = ivstart, ivend
        IF (zwsnew(i) > eps_soil) THEN        ! points with snow cover only
          IF(rho_snow_mult_now(i,ksn) .LT. 600._wp .AND. &
            rho_snow_mult_now(i,ksn) .NE. 0.0_wp) THEN
            zdens_old = rho_snow_mult_now(i,ksn)
            zeta =         &! compactive viscosity of snow
              ca2*EXP(19.3_wp*rho_snow_mult_now(i,ksn)/rho_i)* &
              EXP(67300._wp/8.31_wp/ztsnownew_mult(i,ksn))
            rho_snow_mult_now(i,ksn) = rho_snow_mult_now(i,ksn) + &
              zdt*rho_snow_mult_now(i,ksn)*(csigma+zp(i,ksn))/zeta
            rho_snow_mult_now(i,ksn) = MIN(rho_snow_mult_now(i,ksn),rho_i)
            zdzh_snow(i,ksn)   = zdzh_snow(i,ksn) * zdens_old/rho_snow_mult_now(i,ksn)
          END IF
        END IF       ! points with snow cover only
      END DO
    END DO


    DO i = ivstart, ivend
      IF (zwsnew(i) > eps_soil) THEN        ! points with snow cover only

        IF(ztsnownew_mult(i,0) .GT. t0_melt) THEN
          ztsnownew_mult(i,0) = t0_melt
          zdtsnowdt_mult(i,0) = zdtsnowdt_mult(i,0) +     &
                                  (ztsnownew_mult(i,0) - ztsnown_mult(i,0))*z1d2dt
        END IF
      END IF       ! points with snow cover only
    END DO

  END IF  ! multi-layer snow scheme

!------------------------------------------------------------------------------
! Section II.9: Final updating of prognostic values
!------------------------------------------------------------------------------

  IF (lmulti_snow) THEN

    ! First for ksn == 0
    DO i = ivstart, ivend
      t_snow_mult_new  (i,0) = t_snow_mult_now(i,0) + zdt*zdtsnowdt_mult(i,0)
      t_snow_new(i) = t_snow_mult_new (i,0)
    ENDDO


    DO ksn = 1, ke_snow
      DO i = ivstart, ivend
        t_snow_mult_new  (i,ksn) = t_snow_mult_now(i,ksn) + &
          &                              zdt*zdtsnowdt_mult(i,ksn)
        dzh_snow_new     (i,ksn) = zdzh_snow(i,ksn)
        wtot_snow_new    (i,ksn) = wtot_snow_now(i,ksn)
        rho_snow_mult_new(i,ksn) = rho_snow_mult_now(i,ksn)
        wliq_snow_new    (i,ksn) = wliq_snow_now(i,ksn)
      ENDDO
    ENDDO

  ELSE
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO i = ivstart, ivend
      t_snow_new(i)  = t_snow_now(i) + zdt*zdtsnowdt(i)
    ENDDO
    !$ACC END PARALLEL
  ENDIF

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    ! t_snow is computed above
    ! t_snow(i,nnew)  = t_snow(i,nx) + zdt*zdtsnowdt(i)
    t_so_new(i,1)  = t_so_now(i,1) + zdt*zdtsdt   (i)         ! (*)

    ! Next line has to be changed, if the soil surface temperature
    ! t_so(i,0,nnew) predicted by the heat conduction equation is used
    t_s_new   (i)    = t_so_new(i,1)
    t_so_new  (i,0)  = t_so_new(i,1)
    w_snow_new(i)  = w_snow_now(i) + zdt*zdwsndt  (i)/rho_w
    IF (itype_canopy == 1) THEN
      t_sk_new(i)    = t_s_new(i)
    ELSE IF (itype_canopy == 2) THEN

      ! Calculation of the skin temperature (snow free area)
      ! by Schulz and Vogel (2020), based on Viterbo and Beljaars (1995).
      ! A Newtonian relaxation approach is used to ensure numerical stability.

      IF (w_snow_now(i) > eps_soil .OR. w_snow_new(i) > eps_soil) THEN
        t_sk_new(i) = t_s_new(i) ! needs to be t_s rather than t_snow in order to obtain correct t_g afterwards
      ELSE
        t_sk_new(i) = t_sk_now(i) + 0.5_wp*(t_s_new(i) - t_s_now(i)) + zdt/tau_skin *                                     &
          ( (zrnet_s(i) + zshfl_s(i) + zlhfl_s(i) + zsprs(i) + zahf(i))/MAX(5.0_wp,skinc(i)) - (t_sk_now(i) - t_s_now(i)) )
      ENDIF

    END IF
    IF (itype_interception == 1) THEN
      w_i_new   (i)  = w_i_now(i) + zdt*zdwidt   (i)/rho_w
    ELSE IF (itype_interception == 2) THEN
      w_i_new (i)  = w_i_now(i) + zdt*zdwidt   (i)/rho_w
      w_p_new (i)  = zwpnstr(i)
      w_s_new (i)  = zwisnstr(i)
    END IF
  END DO

  ! Reset t_snow_new to t_so(0) if no snow was present at the beginning of the time step
  ! The heat balance calculation is incomplete in this case and sometimes yields unreasonable results
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (w_snow_now(i) < eps_soil .AND. w_snow_new(i) >= eps_soil) THEN
      t_snow_new(i) = MIN(t0_melt,t_so_new(i,0))
    ENDIF
  ENDDO

  IF (lmulti_snow) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (w_snow_now(i) < eps_soil .AND. w_snow_new(i) >= eps_soil) THEN
        t_snow_mult_new(i,:) = t_snow_new(i)
      ENDIF
    ENDDO
  ENDIF


  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    ! Next line has to be changed, if the soil surface temperature
    ! t_so(i,0,nnew) predicted by the heat conduction equation is used
    t_s_new   (i)    = t_so_new(i,1)
    t_so_new  (i,0)  = t_so_new(i,1)
    IF (itype_canopy == 1) THEN
      t_sk_new(i)    = t_s_new(i)
    END IF
    w_snow_new(i)  = w_snow_now(i) + zdt*zdwsndt  (i)/rho_w
    IF (itype_interception == 1) THEN
      w_i_new   (i)  = w_i_now(i) + zdt*zdwidt   (i)/rho_w
    ELSE IF (itype_interception == 2) THEN
      w_i_new   (i)  = w_i_now(i) + zdt*zdwidt   (i)/rho_w
      w_p_new (i)  = zwpnstr(i)
      w_s_new (i)  = zwisnstr(i)
    END IF
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (w_snow_new(i) <= eps_soil) THEN
      w_i_new(i)    = w_i_new(i) + w_snow_new(i)
      w_snow_new(i) = 0.0_wp
      t_snow_new(i) = t_so_new(i,0)
    ENDIF

    IF (w_i_new(i) <= 1.0E-4_wp*eps_soil) w_i_new(i) = 0.0_wp
  END DO
  !$ACC END PARALLEL
!<JH


  ! Eliminate snow for multi-layer snow model, if w_snow = 0
  IF (lmulti_snow) THEN
    DO ksn = 1, ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (w_snow_new(i) <= eps_soil) THEN
          t_snow_mult_new(i,ksn) = t_so_new(i,0)
          wliq_snow_new(i,ksn) = 0.0_wp
          wtot_snow_new(i,ksn) = 0.0_wp
          rho_snow_mult_new(i,ksn) = 0.0_wp
          dzh_snow_new(i,ksn) = 0.0_wp
        ENDIF
      END DO
    END DO
  ENDIF


  IF(.NOT. lmulti_snow) THEN

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP GANG VECTOR PRIVATE(zzz, ztau_snow, zrhosmax, zrho_snowe) &
    !$ACC   PRIVATE(zrho_snowf, zxx, zdwgme, zro, zredfu, zw_ovpv, zrhosmin)
    !$NEC sparse
    DO i = ivstart, ivend

      !BR 7/2005 Update snow density
      ! a) aging of existing snow
      !    temperature dependence of relaxation/ageing constant
      zzz       = (t_snow_new(i)-csnow_tmin)/(t0_melt-csnow_tmin)
      ztau_snow = crhosmint+(crhosmaxt-crhosmint)*zzz
      ztau_snow = MAX(0.05_wp,MIN(crhosmaxt,ztau_snow)) ! use 20 days in combination with 
                                                        ! temperature-dependent equilibrium density
      zrhosmax  = crhosmax_tmin+MAX(-0.25_wp,zzz)*(crhosmax_ml-crhosmax_tmin)
      zrho_snowe= MAX(rho_snow_now(i),zrhosmax+(rho_snow_now(i)-zrhosmax)* &
                  EXP(-ztau_snow*zdt/86400._wp) )

      ! b) density of fresh snow; the minimum density depends on wind speed to account for wind compression of powder snow
      zuv = SQRT ( u_10m(i)**2 + v_10m(i)**2 )
      zrhosmin = crhosminf*MIN(2._wp, 1._wp + 0.125_wp*MAX(0._wp, zuv-4._wp))
      zrho_snowf= zrhosmin+(crhosmaxf-zrhosmin)* ((zth_low(i)-csnow_tmin) / (t0_melt-csnow_tmin))**2
      zrho_snowf= MAX(zrhosmin,MIN(crhosmaxf,zrho_snowf))

      zrho_grauf= crhogminf+(crhogmaxf-crhogminf)* ((zth_low(i)-csnow_tmin)/(t0_melt-csnow_tmin))**2
      zrho_grauf= MAX(crhogminf,MIN(crhogmaxf,zrho_grauf))

      ! graupel fraction
      IF ( nclass_gscp >= 2000 ) THEN
        ! only possible when running 2-moment microphysics
#ifdef TWOMOM_SB
        zgrfrac = (prh_gsp(i)+prg_gsp(i)+prs_con(i)) / MAX(eps_soil,prs_gsp(i)+pri_gsp(i)+prs_con(i)+prg_gsp(i)+prh_gsp(i))
#endif
      ELSEIF ( nclass_gscp >= 6 ) THEN
        zgrfrac = (prg_gsp(i)+prs_con(i)) / MAX(eps_soil,prs_gsp(i)+pri_gsp(i)+prs_con(i)+prg_gsp(i))
      ELSE
        zgrfrac =             prs_con(i)  / MAX(eps_soil,prs_gsp(i)+pri_gsp(i)+prs_con(i))
      ENDIF

      ! c) new snow density is computed by adding depths of existing and new snow
      IF ( nclass_gscp >= 2000 ) THEN
       ! only possible when running the 2-moment microphysics
       !!$ UB: does that really make sense to integrate hail into snow density at the surface?
#ifdef TWOMOM_SB
        zzz = (prs_gsp(i)+pri_gsp(i)+prs_con(i)+prg_gsp(i)+prh_gsp(i))*zdtdrhw
#endif
      ELSEIF ( nclass_gscp >= 6 ) THEN
        zzz = (prs_gsp(i)+pri_gsp(i)+prs_con(i)+prg_gsp(i))*zdtdrhw
      ELSE
        zzz = (prs_gsp(i)+pri_gsp(i)+prs_con(i))*zdtdrhw
      ENDIF

      ! prevent accumulation of new snow if the air temperature is above 1 deg C with
      ! linear transition between 0.5 and 1 deg C
      IF (zzz > 0.5_wp*eps_soil .AND. zth_low(i) > t0_melt + 0.5_wp) THEN
        ! part of the new snow that accumulates on the ground
        zxx = MAX(0._wp, zzz*(t0_melt + 1._wp - zth_low(i))*2._wp)
        !
        ! the rest is transferred to the interception storage, soil moisture or runoff:
        w_i_new(i) = w_i_new(i) + zzz-zxx
        IF (w_i_new(i) > zwimax(i)) THEN  ! overflow of interception store
          zinf = w_i_new(i) - zwimax(i)
          w_i_new(i) = zwimax(i)
        ELSE
          zinf = 0.0_wp
        ENDIF
        zdwgme        = zinf*zrock(i)                    ! contribution to w_so
        zro           = (1._wp - zrock(i))*zinf      ! surface runoff
        zredfu        = MAX( 0.0_wp,  MIN( 1.0_wp, (zw_fr(i,1) -  &
                        zfcap(i,1))/MAX(zporv(i,1)-zfcap(i,1), eps_div)))
        zdwgdt(i,1) = zdwgdt(i,1) + zdwgme*(1._wp - zredfu)
        zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                               ! for this fraction

        ! zro-, zdw_so_dt-correction in case of pore volume overshooting
        zw_ovpv = MAX(0._wp, zw_fr(i,1)* zdzhs(1) * zrhwddt +  &
                  zdwgdt(i,1) - zporv(i,1) * zdzhs(1) * zrhwddt)
        zro = zro + zw_ovpv
        zdwgdt(i,1)= zdwgdt(i,1) - zw_ovpv

        runoff_s(i) = runoff_s(i) + zro*zroffdt

        ! correct SWE for immediately melted new snow
        zzz = zxx
        w_snow_new(i) = MIN(w_snow_new(i),w_snow_now(i)+zzz)
      ENDIF

      IF (soiltyp_subs(i) /= 1) THEN
        rho_snow_new(i)  = (w_snow_now(i)+zzz) / &
         ( MAX(w_snow_now(i),eps_soil)/zrho_snowe + zzz/( (1.0_wp-zgrfrac)*zrho_snowf + zgrfrac*zrho_grauf) )
      ELSE ! constant snow density over glaciers
        rho_snow_new(i) = rho_snow_now(i)
      ENDIF

      ! previous code based on weighted averaging of rho_snow
      !       znorm=MAX(w_snow_now(i)+(prs_gsp(i)+prs_con(i)+prg_gsp(i))      &
      !                 *zdtdrhw,eps_div)
      !       rho_snow_new(i)  = ( zrho_snowe*w_snow_now(i) + &
      !                           zrho_snowf*(prs_gsp(i)+prs_con(i)+prg_gsp(i)) &
      !                              *zdtdrhw )    /znorm
      !     ELSE
      !       znorm=MAX(w_snow_now(i)+(prs_gsp(i)+prs_con(i) )      &
      !                 *zdtdrhw,eps_div)
      !       rho_snow_new(i)  = ( zrho_snowe*w_snow_now(i) + &
      !                          zrho_snowf*(prs_gsp(i)+prs_con(i) ) &
      !                             *zdtdrhw )    /znorm
      !     ENDIF

      rho_snow_new(i) = MIN(crhosmax_ml,MAX(crhosmin_ml, rho_snow_new(i)))

      ! New calculation of snow height for single layer snow model
      zdz_snow(i)=w_snow_new(i)*rho_w/rho_snow_new(i)
      h_snow_new(i) = zdz_snow(i)

      ! Calculation of top-layer snow density for two-layer snow density scheme
      IF (l2lay_rho_snow) THEN
        zrho_snowe = MAX(rho_snow_mult_now(i,1),zrhosmax+(rho_snow_mult_now(i,1)-zrhosmax)* &
                     EXP(-ztau_snow*zdt/86400._wp) )
        zwsnow(i)  = MIN(max_toplaydepth,h_snow_gp(i))*rho_snow_mult_now(i,1)/rho_w
        rho_snow_mult_new(i,1) = (zwsnow(i)+zzz) / ( MAX(zwsnow(i),eps_div)/zrho_snowe + zzz/zrho_snowf )
        rho_snow_mult_new(i,1) = MIN(crhosmax_ml,MAX(crhosmin_ml, rho_snow_mult_new(i,1)))
        rho_snow_mult_new(i,2) = rho_snow_new(i)
      ENDIF

    END DO
    !$ACC END PARALLEL

  ELSE   ! new snow scheme

    DO i = ivstart, ivend
      h_snow_new(i) = 0.0_wp
      sum_weight(i) = 0.0_wp
    END DO

    DO ksn = 1,ke_snow
      DO i = ivstart, ivend
        IF(w_snow_new(i) .GT. eps_soil) THEN
          h_snow_new(i) = h_snow_new(i) + zdzh_snow(i,ksn)
        END IF
      END DO
    END DO

    k = MIN(2,ke_snow-1)
    DO ksn = 1,ke_snow
      DO i = ivstart, ivend
        IF (w_snow_new(i) .GT. eps_soil) THEN
          IF (ksn == 1) THEN ! Limit top layer to max_toplaydepth
            zhh_snow(i,ksn) = -MAX( h_snow_new(i)-max_toplaydepth, h_snow_new(i)/ke_snow*(ke_snow-ksn) )
          ELSE IF (ksn == 2 .AND. ke_snow > 2) THEN ! Limit second layer to 8*max_toplaydepth
            zhh_snow(i,ksn) = MIN( 8._wp*max_toplaydepth+zhh_snow(i,1), zhh_snow(i,1)/(ke_snow-1)*(ke_snow-ksn) )
          ELSE ! distribute the remaining snow equally among the layers
            zhh_snow(i,ksn) = zhh_snow(i,k)/(ke_snow-k)*(ke_snow-ksn)
          ENDIF
        ENDIF
      END DO
    END DO

    DO ksn = ke_snow,1,-1
      DO i = ivstart, ivend
        IF (w_snow_new(i) .GT. eps_soil) THEN
          dz_old(i,ksn) = dzh_snow_new(i,ksn)
          z_old(i,ksn) = -sum_weight(i) - dzh_snow_new(i,ksn)/2._wp
          sum_weight(i) = sum_weight(i) + dzh_snow_new(i,ksn)
        END IF
      END DO
    END DO

    DO i = ivstart, ivend
       IF(w_snow_new(i) .GT. eps_soil) THEN
         zhm_snow(i,1) = (-h_snow_new(i) + zhh_snow(i,1))/2._wp

         ! layer thickness betw. half levels of uppermost snow layer
         dzh_snow_new(i,1) = zhh_snow(i,1) + h_snow_new(i)      

         ! layer thickness between snow surface and main level of uppermost layer
         zdzm_snow(i,1     ) = zhm_snow(i,1) + h_snow_new(i)    

         IF(dz_old(i,1).ne.0..and.rho_snow_mult_new(i,1).ne.0.) THEN
           wliq_snow_new(i,1) = wliq_snow_new(i,1)/dz_old(i,1)
         END IF
       END IF
    END DO

    DO ksn = 2,ke_snow
      DO i = ivstart, ivend
        IF(w_snow_new(i) .GT. eps_soil) THEN
          zhm_snow(i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._wp
          dzh_snow_new(i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
          zdzm_snow(i,ksn     ) = zhm_snow(i,ksn) - zhm_snow(i,ksn-1) ! layer thickness betw. main levels
          IF(dz_old(i,ksn).ne.0..and.rho_snow_mult_new(i,ksn).ne.0.) THEN
            wliq_snow_new(i,ksn) = wliq_snow_new(i,ksn)/dz_old(i,ksn)
          END IF
        END IF
      END DO
    END DO

    DO ksn = ke_snow,1,-1
      DO i = ivstart, ivend
        t_new  (i,ksn) = 0.0_wp
        rho_new(i,ksn) = 0.0_wp
        wl_new (i,ksn) = 0.0_wp
      END DO

      DO k = ke_snow,1,-1
        DO i = ivstart, ivend
          IF(w_snow_new(i) .GT. eps_soil) THEN

            weight = MIN(&
                 dz_old(i,k), &
                 z_old(i,k) + dz_old(i,k)/2._wp &
                 - zhm_snow(i,ksn) + dzh_snow_new(i,ksn)/2._wp , &
                 zhm_snow(i,ksn) + dzh_snow_new(i,ksn)/2._wp &
                 - z_old(i,k) + dz_old(i,k)/2._wp, &
                 dzh_snow_new(i,ksn))

            weight = (weight + ABS(weight)) * 0.5_wp / dzh_snow_new(i,ksn)

            t_new  (i,ksn) = t_new  (i,ksn) + t_snow_mult_new  (i,k)*weight
            rho_new(i,ksn) = rho_new(i,ksn) + rho_snow_mult_new(i,k)*weight
            wl_new (i,ksn) = wl_new (i,ksn) + wliq_snow_new(i,k)*weight
          END IF
        END DO
      END DO
    END DO

    DO ksn = ke_snow,1,-1
      DO i = ivstart, ivend
        IF(w_snow_new(i) > eps_soil) THEN
          t_snow_mult_new  (i,ksn) = t_new  (i,ksn)
          rho_snow_mult_new(i,ksn) = rho_new(i,ksn)
          wtot_snow_new    (i,ksn) = rho_new(i,ksn)*dzh_snow_new(i,ksn)/rho_w
          wliq_snow_new    (i,ksn) = wl_new (i,ksn)*dzh_snow_new(i,ksn)
        END IF
      END DO
    END DO


    DO i = ivstart, ivend
      IF(w_snow_new(i) > eps_soil) THEN
        rho_snow_new(i) = w_snow_new(i)/h_snow_new(i)*rho_w
      ELSE !JH
        rho_snow_new(i) = 250._wp ! workaround need to be inspected!!
      END IF
      IF(w_snow_new(i) > eps_soil) THEN
        ! linear extrapolation from t_snow_mult_new(i,2) and t_snow_mult_new(i,1) to t_snow_mult_new(i,0)
        t_snow_mult_new(i,0) = (t_snow_mult_new(i,1)*(2._wp*dzh_snow_new(i,1)+dzh_snow_new(i,2))- &
                              t_snow_mult_new(i,2)*dzh_snow_new(i,1))/(dzh_snow_new(i,1)+dzh_snow_new(i,2))
        ! limiter to prevent unphysical values and/or numerical instabilities
        t_snow_mult_new(i,0) = MIN(273.15_wp,t_snow_mult_new(i,0),t_snow_mult_new(i,1)+5.0_wp)
        t_snow_mult_new(i,0) = MAX(t_snow_mult_new(i,0),t_snow_mult_new(i,1)-5.0_wp)
      END IF
      t_snow_new(i) = t_snow_mult_new (i,0)
    END DO

 ENDIF ! lmulti_snow


  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP SEQ
  DO kso = 1,ke_soil
  !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      w_so_new(i,kso) = w_so_now(i,kso) + zdt*zdwgdt(i,kso)/rho_w
    END DO
  END DO        ! soil layers

  ! Update of two-time level interface variables
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    h_snow(i) = h_snow_new(i)
  END DO
  !$ACC END PARALLEL

  ! computation of the weighted turbulent fluxes at the boundary surface-atmosphere
  IF(PRESENT(zshfl_sfc)) THEN
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO i = ivstart, ivend
      zshfl_sfc(i) = zshfl_s(i)*(1._wp - zf_snow(i)) + zshfl_snow(i)*zf_snow(i)
      zlhfl_sfc(i) = zlhfl_s(i)*(1._wp - zf_snow(i)) + zlhfl_snow(i)*zf_snow(i)
      zqhfl_sfc(i) = zqhfl_s(i)*(1._wp - zf_snow(i)) + zqhfl_snow(i)*zf_snow(i)

      !        zlhfl_s(i) = (zts_pm(i)*lh_v + (1._wp-zts_pm(i))*lh_s)*zverbo(i) &
      !                     / MAX(eps_div,(1._wp - zf_snow(i)))  ! take out (1-f) scaling
      !        zlhfl_snow(i) = lh_s*zversn(i)
      !      zlhfl_sfc(i) = zverbo(i) + zversn(i)*zf_snow(i)
    END DO
    !$ACC END PARALLEL
  END IF

  ! This block tests the residuum of water mass content in soil
  IF (lres_soilwatb) THEN
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue)
    !$ACC LOOP GANG VECTOR
    DO i = ivstart, ivend
      ! the rain rates 
      IF ( nclass_gscp >= 2000 ) THEN
        ! only possible when running 2-moment microphysics
#ifdef TWOMOM_SB
        zhwso_budget(i) = zhwso_budget(i) + &
          ( prr_con(i) + prs_con(i) + prr_gsp(i) + prs_gsp(i) + prg_gsp(i) + prh_gsp(i) )*zdt
#endif
      ELSEIF ( nclass_gscp >= 6 ) THEN
        zhwso_budget(i) = zhwso_budget(i) + &
          ( prr_con(i) + prs_con(i) + prr_gsp(i) + prs_gsp(i) + prg_gsp(i) )*zdt
      ELSE
        zhwso_budget(i) = zhwso_budget(i) + &
          ( prr_con(i) + prs_con(i) + prr_gsp(i) + prs_gsp(i) )*zdt
      ENDIF
      ! the evapotranspiration (not clear if zqhfl_sfc or zqhfl_s???????)
      zhwso_budget(i) = zhwso_budget(i) + zqhfl_sfc(i)*zdt
      ! surface + subsurface runoff (subtraction because water is lost)
      zhwso_budget(i) = zhwso_budget(i) - runoff_s(i) - runoff_g(i)
      ! snow + interception storage (convert m H2O in mm H2O)
      zhwso_budget(i) = zhwso_budget(i) - ( w_snow_new(i) - w_snow_now(i) )*1000.0_wp
      zhwso_budget(i) = zhwso_budget(i) - ( w_i_new(i)    - w_i_now(i)    )*1000.0_wp
      IF (itype_interception == 2) THEN
        zhwso_budget(i) = zhwso_budget(i) - ( w_p_new(i)    - w_p_now(i)    )*1000.0_wp
        zhwso_budget(i) = zhwso_budget(i) - ( w_s_new(i)    - w_s_now(i)    )*1000.0_wp
      ENDIF
      ! residuum = recharge - budget (RHS)
      resid_wso(i) = resid_wso(i) + ( ( SUM(w_so_new(i,1:ke_soil_hy)) - & 
        SUM(w_so_now(i,1:ke_soil_hy)) )*1000.0_wp - zhwso_budget(i) )
    END DO
    !$ACC END PARALLEL
  ENDIF

  IF (.NOT. lcuda_graph_lnd) THEN
    !$ACC WAIT(acc_async_queue)
  END IF

! for optional fields related to soil water budget
  !$ACC END DATA

! for optional fields plevap, z0
  !$ACC END DATA
!!!if (PRESENT(plevap))

! for general fields
  !$ACC END DATA

  IF (msg_level >= 19) THEN
    !$ACC UPDATE HOST(ivend) ASYNC(acc_async_queue)
    IF (.NOT. lcuda_graph_lnd) THEN
      !$ACC WAIT(acc_async_queue)
    END IF
    DO i = ivstart, ivend

     IF (ABS(t_s_now(i)-t_s_new(i)) > 15.0_wp .or. ABS(t_sk_now(i)-t_sk_new(i)) > 15.0_wp) THEN
!     the setting below was only useful in the frame of the COSMO model
!     IF (i== mvid .AND. iblock == mbid .AND. my_cart_id == mcid) THEN

        WRITE(*,'(A        )') '                                '
        WRITE(*,'(A,2I5)'  ) 'SFC-DIAGNOSIS terra output:  iblock = ', iblock, i

        WRITE(*,'(A        )') ' Temperatures and Humidities: '
        WRITE(*,'(A,2F28.16)') '   t_s      now/new :  ', t_s_now(i),       t_s_new(i)
        WRITE(*,'(A,2F28.16)') '   t_snow   now/new :  ', t_snow_now(i),    t_snow_new(i)
        WRITE(*,'(A, F28.16)') '   t_g              :  ', t_g(i)
        WRITE(*,'(A, F28.16)') '   qv_s (out)       :  ', qv_s(i)
        WRITE(*,'(A,2F28.16)') '   w_snow   now/new :  ', w_snow_now(i),    w_snow_new(i)
        WRITE(*,'(A,2F28.16)') '   rho_snow now/new :  ', rho_snow_now(i),  rho_snow_new(i)
        WRITE(*,'(A,2F28.16)') '   h_snow   now/new :  ', h_snow_now(i),    h_snow_new(i)
        WRITE(*,'(A, F28.16)') '   fresh_snow       :  ', freshsnow(i)
        WRITE(*,'(A, F28.16)') '   zf_snow (out)    :  ', zf_snow(i)
        WRITE(*,'(A, F28.16)') '   sn_frac (out)    :  ', sn_frac(i)
        WRITE(*,'(A,2F28.16)') '   w_i      now/new :  ', w_i_now(i),       w_i_new(i)
        WRITE(*,'(A,2F28.16)') '   w_p      now/new :  ', w_p_now(i),       w_p_new(i)
        WRITE(*,'(A,2F28.16)') '   w_s      now/new :  ', w_s_now(i),       w_s_new(i)
        DO k = 0, ke_soil+1
          WRITE(*,'(A,I1,A,2F28.16)') '   t_so    (',k,')      :  ', t_so_now (i,k), t_so_new (i,k)
        ENDDO
        DO k = 1, ke_soil+1
          WRITE(*,'(A,I1,A,2F28.16)') '   w_so    (',k,')      :  ', w_so_now (i,k), w_so_new (i,k)
        ENDDO
        DO k = 1, ke_soil+1
          WRITE(*,'(A,I1,A,2F28.16)') '   w_so_ice(',k,')      :  ', w_so_ice_now (i,k), w_so_ice_new (i,k)
        ENDDO
        WRITE(*,'(A        )') '                                '
        WRITE(*,'(A        )') ' Fluxes etc.            :  '
        WRITE(*,'(A, F28.16)') '   tcm              :  ', tcm(i)
        WRITE(*,'(A, F28.16)') '   tch              :  ', tch(i)
        WRITE(*,'(A, F28.16)') '   runoff_s (out)       :  ', runoff_s(i)
        WRITE(*,'(A, F28.16)') '   runoff_g (out)       :  ', runoff_g(i)
        WRITE(*,'(A, F28.16)') '   rstom            : ', rstom(i)
        IF (itype_trvg == 3) THEN
          WRITE(*,'(A, F28.16)') '   plevap               :  ', plevap(i)
        ENDIF
        WRITE(*,'(A,2F28.16)') '   zshfl/zlhfl (surface):  ', zshfl_sfc(i), zlhfl_sfc(i)
        WRITE(*,'(A, F28.16)') '   zqhfl (surface)      :  ', zqhfl_sfc(i)
        WRITE(*,'(A,2F28.16)') '   zshfl/zlhfl (soil)   :  ', zshfl_s(i), zlhfl_s(i)
        WRITE(*,'(A, F28.16)') '   zqhfl (soil)         :  ', zqhfl_s(i)
        WRITE(*,'(A,2F28.16)') '   zshfl/lhfl  (snow)   :  ', zshfl_snow(i), zlhfl_snow(i)
        WRITE(*,'(A, F28.16)') '   zqhfl (snow)         :  ', zqhfl_snow(i)
        WRITE(*,'(A, F28.16)') '   zgsb (heat flux throw snow) :  ', zgsb(i)
        WRITE(*,'(A, F28.16)') '   lhfl_bs              :  ', lhfl_bs(i)
        DO k = 1, ke_soil
          WRITE(*,'(A,I1,A, F28.16)') '   lhfl_pl (',k,')      :  ', lhfl_pl(i,k)
        ENDDO
        DO k = 1, ke_soil+1
          WRITE(*,'(A,I1,A, F28.16)') '   zrunoff_grav [kg/m2/s] (',k,')      :  ', zrunoff_grav(i,k)
        ENDDO

      ENDIF
    ENDDO
  ENDIF




!------------------------------------------------------------------------------
! End of module procedure terra
!------------------------------------------------------------------------------

END SUBROUTINE terra

!==============================================================================

REAL (KIND=wp) FUNCTION zsf_heav (zstx)
  ! Heaviside function

  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN)  :: zstx
  zsf_heav = 0.5_wp + SIGN (0.5_wp, zstx)

END FUNCTION zsf_heav

!==============================================================================

REAL (KIND=wp) FUNCTION zsf_psat_iw  (zstx, z2iw, z4iw)
  ! Saturation water vapour pressure over ice or water depending on temperature "zstx"

  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN)  :: zstx, z2iw, z4iw
  zsf_psat_iw   = b1*EXP(z2iw*(zstx - b3)/(zstx - z4iw))

END FUNCTION zsf_psat_iw

!==============================================================================

REAL (KIND=wp) FUNCTION zsf_qsat  (zspsatx, zspx)
  ! Specific humidity at saturation pressure (depending on the saturation water
  !  vapour pressure zspsatx" and the air pressure "zspx")

  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN)  :: zspsatx, zspx
  zsf_qsat      = rdv*zspsatx/(zspx-o_m_rdv*zspsatx)

END FUNCTION zsf_qsat

!==============================================================================

REAL (KIND=wp) FUNCTION zsf_dqvdt_iw (zstx, zsqsatx, z4iw, z234iw)
  ! First derivative of specific saturation humidity with respect to temperature
  ! (depending on temperature "zstx" and saturation specific humidity pressure
  !  "zsqsatx")

  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN)  :: zstx, zsqsatx, z4iw, z234iw
  zsf_dqvdt_iw  = z234iw*(1.0_wp+rvd_m_o*zsqsatx)*zsqsatx/(zstx-z4iw)**2

END FUNCTION zsf_dqvdt_iw

!==============================================================================


REAL (KIND=wp) FUNCTION watrcon_RT(ks,lw,kw1,pv,adp)

  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN)  :: ks,lw,kw1,pv,adp
  watrcon_RT=ks*EXP(kw1*MAX(0.0_wp,pv-lw)/(pv-adp))

END FUNCTION watrcon_RT

!==============================================================================

REAL (KIND=wp) FUNCTION watrdiff_RT(ds,lw,dw1,pv,adp)

  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN)  ::  ds,lw,dw1,pv,adp
  watrdiff_RT = ds*EXP(dw1*MAX(0.0_wp,pv-lw)/(pv-adp))

END FUNCTION watrdiff_RT

!==============================================================================

!------------------------------------------------------------------------------
! End of module sfc_terra
!------------------------------------------------------------------------------

END MODULE sfc_terra
