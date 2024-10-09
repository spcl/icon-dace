!> vegetation process memory (QUINCY)
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### definition and init of (memory) variables for the vegetation process
!>
MODULE mo_veg_memory_class
#ifndef __NO_QUINCY__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
  USE mo_util,                   ONLY: One_of, int2string
  USE mo_jsb_utils_iface,        ONLY: assign_if_present_allocatable

  USE mo_jsb_process_class,      ONLY: VEG_
  USE mo_lnd_bgcm_class,         ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, &
    &                                  LAST_ELEM_ID, get_name_for_element_id, get_short_name_for_element_id
  USE mo_lnd_bgcm_store_class,   ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_TOTAL_BIO_ID, VEG_BGCM_GROWTH_ID,               &
    &                                  VEG_BGCM_LITTERFALL_ID, VEG_BGCM_FRAC_ALLOC_ID,                            &
    &                                  VEG_BGCM_EXUDATION_ID, VEG_BGCM_ESTABLISHMENT_ID, VEG_BGCM_RESERVE_USE_ID, &
    &                                  VEG_BGCM_PP_FUEL_ID, VEG_BGCM_PP_PAPER_ID, VEG_BGCM_PP_FIBERBOARD_ID,      &
    &                                  VEG_BGCM_PP_OIRW_ID, VEG_BGCM_PP_PV_ID, VEG_BGCM_PP_SAWNWOOD_ID, VEG_BGCM_FPROD_DECAY_ID
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_pool_class,         ONLY: t_jsb_pool, ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
  PUBLIC :: t_veg_bgcm_main, t_veg_bgcm_with_components, t_veg_bgcm_with_elements

  INTEGER, PARAMETER :: max_no_of_vars = 1000

  ! ======================================================================================================= !
  !>Type definition for veg memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_veg_memory

    !------------------------------------------------------------------------------------------------------ !
    !> bgc_material structures
    !>
    TYPE(t_veg_bgcm_main),  POINTER :: veg_bgcm_main

    !------------------------------------------------------------------------------------------------------ !
    !> bgc_material diagnostics - summary of C, N, P, C13, C14, N15 across the particular bgcm
    !> used for output in IQ
    !> because the bgcm cannot be used with the 'aggregate()' routine yet @TODO
    !>
    TYPE(t_jsb_var_real2d) :: &
      & veg_pool_total_c          , &     !< sum of all veg_pool%*%carbon      [mol m-2]
      & veg_pool_total_n          , &     !< sum of all veg_pool%*%nitrogen    [mol m-2]
      & veg_pool_total_p          , &     !< sum of all veg_pool%*%phosphorus  [mol m-2]
      & veg_pool_total_c13        , &     !< sum of all veg_pool%*%carbon13    [mol m-2]
      & veg_pool_total_c14        , &     !< sum of all veg_pool%*%carbon14    [mol m-2]
      & veg_pool_total_n15                !< sum of all veg_pool%*%nitrogen15  [mol m-2]
    TYPE(t_jsb_var_real2d) :: &
      & veg_growth_total_c        , &     !< sum of all veg_growth%*%carbon      [mol m-2 timestep-1]
      & veg_growth_total_n        , &     !< sum of all veg_growth%*%nitrogen    [mol m-2 timestep-1]
      & veg_growth_total_p        , &     !< sum of all veg_growth%*%phosphorus  [mol m-2 timestep-1]
      & veg_growth_total_c13      , &     !< sum of all veg_growth%*%carbon13    [mol m-2 timestep-1]
      & veg_growth_total_c14      , &     !< sum of all veg_growth%*%carbon14    [mol m-2 timestep-1]
      & veg_growth_total_n15              !< sum of all veg_growth%*%nitrogen15  [mol m-2 timestep-1]
    TYPE(t_jsb_var_real2d) :: &
      & veg_litterfall_total_c    , &     !< sum of all veg_litterfall%*%carbon      [mol m-2 timestep-1]
      & veg_litterfall_total_n    , &     !< sum of all veg_litterfall%*%nitrogen    [mol m-2 timestep-1]
      & veg_litterfall_total_p    , &     !< sum of all veg_litterfall%*%phosphorus  [mol m-2 timestep-1]
      & veg_litterfall_total_c13  , &     !< sum of all veg_litterfall%*%carbon13    [mol m-2 timestep-1]
      & veg_litterfall_total_c14  , &     !< sum of all veg_litterfall%*%carbon14    [mol m-2 timestep-1]
      & veg_litterfall_total_n15          !< sum of all veg_litterfall%*%nitrogen15  [mol m-2 timestep-1]
    TYPE(t_jsb_var_real2d) :: &
      & veg_pool_leaf_c           , &     !< leaf carbon                  [mol m-2]
      & veg_pool_leaf_n           , &     !< leaf nitrogen                [mol m-2]
      & veg_pool_wood_c           , &     !< sapwood + heartwood carbon   [mol m-2]
      & veg_pool_wood_n           , &     !< sapwood + heartwood nitrogen [mol m-2]
      & veg_pool_fine_root_c      , &     !< fine root carbon             [mol m-2]
      & veg_pool_fine_root_n              !< fine root nitrogen           [mol m-2]
    TYPE(t_jsb_var_real2d) :: &
      & veg_products_total_c, &           !< total c in product pools     [mol m-2]
      & veg_products_total_n              !< total n in product pools     [mol m-2]
    !------------------------------------------------------------------------------------------------------ !
    !>basic plant 2D
    !>
    TYPE(t_jsb_var_real2d)          :: &
                                    height            , &            !< plant height [m]
                                    diameter          , &            !< plant diameter [m]
                                    dens_ind          , &            !< density of individuals [# m-2]
                                    lai               , &            !< actual leaf area index [m2 m-2]
                                    sai               , &            !< actual stem area index [m2 m-2]
                                    target_cn_leaf    , &            !< target foliar C:N ratio [mol C mol-1 N]
                                    target_np_leaf    , &            !< target foliar N:P ratio [mol N mol-1 P]
                                    mean_leaf_age     , &            !< average foliage age in days [days]
                                    fract_fpc         , &            !< foliage projected cover fraction of the tile
                                    blended_height    , &            !< foliage weighted height of the tile for turbulence calculations [m]
                                    cohort_age        , &            !< age of the cohort in years, used with veg_dynamics_scheme cohort [years]
                                    f_n_demand        , &            !< relative plant demand for N
                                    f_p_demand        , &            !< relative plant demand for P
                                    beta_sinklim_ps   , &            !< scaling factor to account for a direct sink limitation contraint on photosynthesis [unitless]
                                    dphi              , &            !< plant-soil water potential gradient [MPa]
                                    growth_req_n      , &            !< moles N required for a unit of C growth under current allocation fractions [mol N m-2 -1]
                                    growth_req_p      , &            !< moles P required for a unit of N growth under current allocation fractions [mol P m-2 -1]
                                    leaf_cn_direction , &            !< indicates direction of change in leaf CN ratio to maximise NPP [unitless]
                                    target_lai        , &            !< LAI implied from sap wood area [m2 m-2]
                                    kstar_labile      , &            !< labile pool turnover rate assuming only temperature and moisture constraints, no stoichiometry [day-1]
                                    root_limitation_state            !< dominating limiting factor for root growth [-1,1] (P = -1 | Water = 0 | N = 1) [unitless]
    ! used as logical: 0 = FALSE / 1 = TRUE
    TYPE(t_jsb_var_real2d)          :: &
                                    do_cohort_harvest                !< logical to identify whether the cohort may be harvested, in veg_dynamics_scheme cohort [logical]
    ! basic plant 3D
    TYPE(t_jsb_var_real3d)          :: &    ! DIMENSION(soil_layer_sb)
                                    root_fraction_sl          , &    !< fractions of root per layer [unitless]
                                    delta_root_fraction_sl    , &    !< change of fractions of root per layer [1/timestep]
                                    lai_cl                    , &    !< leaf area index per canopy layer [m2 m-2]
                                    cumm_lai_cl                      !< cummulative leaf area above the mid-point of the canopy layer [m2 m-2]

    !------------------------------------------------------------------------------------------------------ !
    !>target C:N and N:P ratios for various non-leaf tissues (mol/mol)
    !>
    TYPE(t_jsb_var_real2d)          :: &
                                    target_cn_fine_root       , &   !< target C:N ratio for fine roots [mol/mol]
                                    target_cn_coarse_root     , &   !< target C:N ratio for coarse roots [mol/mol]
                                    target_cn_sap_wood        , &   !< target C:N ratio for sap wood  [mol/mol]
                                    target_cn_fruit           , &   !< target C:N ratio for fruit pool [mol/mol]
                                    target_np_fine_root       , &   !< target N:P ratio for fine roots [mol/mol]
                                    target_np_coarse_root     , &   !< target N:P ratio for coarse roots [mol/mol]
                                    target_np_sap_wood        , &   !< target N:P ratio for sap wood pool [mol/mol]
                                    target_np_fruit                 !< target N:P ratio for fruit pool [mol/mol]

    !------------------------------------------------------------------------------------------------------ !
    !>canopy layer information
    !>
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    leaf_nitrogen_cl  , &           !< N in canopy layer [mmol m-2 LAI]
                                    fn_chl_cl         , &           !< fraction of N in chlorophyll, pigments only [unitless]
                                    fn_et_cl          , &           !< fraction of N in electron transport [unitless]
                                    fn_rub_cl         , &           !< fraction of N in Rubisco [unitless]
                                    fn_pepc_cl        , &           !< fraction of N in PEP carboylase, C4 plants only [unitless]
                                    fn_oth_cl         , &           !< fraction of N not associated with photosynthesis [unitless]
                                    fleaf_sunlit_cl                 !< fraction of layer which is sunlit [unitless]

    !------------------------------------------------------------------------------------------------------ !
    !>fluxes
    !>
    TYPE(t_jsb_var_real2d)          :: &
                                    maint_respiration_pot         , & !< potential maintenance respiration [micro-mol CO2 m-2 s-1]
                                    maint_respiration             , & !< maintenance respiration [micro-mol CO2 m-2 s-1]
                                    maint_respiration_c13         , & !< maintenance respiration [micro-mol 13CO2 m-2 s-1]
                                    maint_respiration_c14         , & !< maintenance respiration [micro-mol 14CO2 m-2 s-1]
                                    growth_respiration            , & !< growth respiration [micro-mol CO2 m-2 s-1]
                                    growth_respiration_c13        , & !< growth respiration [micro-mol 13CO2 m-2 s-1]
                                    growth_respiration_c14        , & !< growth respiration [micro-mol 14CO2 m-2 s-1]
                                    n_transform_respiration       , & !< respiration associated with N transformation [micro-mol CO2 m-2 s-1]
                                    n_fixation_respiration        , & !< respiration associated with N fixation [micro-mol CO2 m-2 s-1]
                                    n_processing_respiration      , & !< respiration associated with N processing [micro-mol CO2 m-2 s-1]
                                    n_processing_respiration_c13  , & !< respiration associated with N processing [micro-mol 13CO2 m-2 s-1]
                                    n_processing_respiration_c14  , & !< respiration associated with N processing [micro-mol 14CO2 m-2 s-1]
                                    npp                           , & !< net primary production [micro-mol CO2 m-2 s-1]
                                    npp_c13                       , & !< net primary production [micro-mol 13CO2 m-2 s-1]
                                    npp_c14                       , & !< net primary production [micro-mol 14CO2 m-2 s-1]
                                    uptake_nh4                    , & !< NH4 uptake from soil [micro-mol N m-2 s-1]
                                    uptake_nh4_n15                , & !< 15NH4 uptake from soil [micro-mol N m-2 s-1]
                                    uptake_no3                    , & !< NO3 uptake from soil [micro-mol N m-2 s-1]
                                    uptake_no3_n15                , & !< 15NO3 uptake from soil [micro-mol N m-2 s-1]
                                    uptake_po4                    , & !< PO4 uptake from soil [micro-mol P m-2 s-1]
                                    n_fixation                    , & !< symbiotic N fixation [micro-mol N m-2 s-1]
                                    n_fixation_n15                , & !< symbiotic N fixation [micro-mol 15N m-2 s-1]
                                    recycling_leaf_n              , & !< net flux of N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_leaf_n15            , & !< net flux of 15N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_leaf_p              , & !< net flux of P from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_n         , & !< net flux of N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_n15       , & !< net flux of 15N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_p         , & !< net flux of P from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_n        , & !< net flux of N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_n15      , & !< net flux of 15N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_p        , & !< net flux of P from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    net_growth                        !< npp minus turnover [micro-mol CO2 m-2 s-1]

    !------------------------------------------------------------------------------------------------------ !
    !>everything else
    !>
    TYPE(t_jsb_var_real2d)          :: &
                                    fmaint_rate_root          , & !< current N-specific maintentance respiration [~micro-mol C mol-1 N s-1]
                                    unit_npp                  , & !< incremental change in NPP [micro-mol CO2 m-2 s-1]
                                    unit_transpiration        , & !< incremental change in transpiration [kg m-2 s-1]
                                    unit_uptake_n_pot         , & !< incremental change in potential N uptake per unit root carbon [micro-mol N m-2 s-1]
                                    unit_uptake_p_pot         , & !< incremental change in potential P uptake per unit root carbon [micro-mol P m-2 s-1]
                                    unit_uptake_n_act         , & !< incremental change in actual N uptake per unit root carbon [micro-mol N m-2 s-1]
                                    unit_uptake_p_act         , & !< incremental change in actual P uptake per unit root carbon [micro-mol P m-2 s-1]
                                    delta_dens_ind            , & !< change in density of individuals [# m-2 timestep-1]
                                    cost_n_uptake_root        , & !< actual cost for N uptake by roots [mol C mol-1 N]
                                    mortality_rate            , & !< background_mort_rate_tree / grass [timestep-1]
                                    k1_opt                    , & !< coefficient for optimal leaf:root ratio [unitless]
                                    k2_opt                    , & !< coefficient for optimal leaf:root ratio [unitless]
                                    leaf2sapwood_mass_ratio   , & !< leaf to shoot mass ratio [unitless]
                                    leaf2root_mass_ratio          !< leaf to root mass ratio [unitless]

    !------------------------------------------------------------------------------------------------------ !
    !>daytime averages
    !>
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_daytime_cl          !< fraction of layer which is sunlit average of the previous day [unitless]
    TYPE(t_jsb_var_real2d)          :: &
                                    t_air_daytime             , &    !< air temperature average of the previous day [K]
                                    press_srf_daytime         , &    !< air pressure at surface average of the previous day [Pa]
                                    co2_mixing_ratio_daytime  , &    !< co2 mixing ratio of the air average of the previous day [co2 in ppm]
                                    ga_daytime                , &    !< aerodynamic conductance average of the previous day [m s-1]
                                    beta_sinklim_ps_daytime   , &    !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                                    t_jmax_opt_daytime               !< optimal temperature for electron transport of photosynthesis average of the previous day [deg C]

    !------------------------------------------------------------------------------------------------------ !
    !>accumlated over current daytime
    !>
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_daytime_dacc_cl     !< fraction of layer which is sunlit value accumulated over current day [unitless]
    TYPE(t_jsb_var_real2d)          :: &
                                    t_air_daytime_dacc            , & !< air temperature average of the previous day [K]
                                    press_srf_daytime_dacc        , & !< air pressure at surface average of the previous day [Pa]
                                    co2_mixing_ratio_daytime_dacc , & !< co2 mixing ratio of the air average of the previous day [co2 in ppm]
                                    ga_daytime_dacc               , & !< aerodynamic conductance average of the previous day [m s-1]
                                    beta_sinklim_ps_daytime_dacc  , & !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                                    t_jmax_opt_daytime_dacc           !< optimal temperature for electron transport of photosynthesis value accumulated over current day [deg C]

    !------------------------------------------------------------------------------------------------------ !
    !>moving averages
    !>
    ! 1.1 conditions for the labile pool dynamics
    TYPE(t_jsb_var_real2d)          :: &
                                    gpp_tlabile_mavg               , &  !< gpp  [micro-mol m-2 s-1]
                                    maint_respiration_tlabile_mavg , &  !< maintenance respiration [micro-mol m-2 s-1]
                                    growth_req_n_tlabile_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions time-averaged [mol N m-2 s-1]
                                    growth_req_p_tlabile_mavg           !< moles P required for a unit of N growth under current allocation fractions time-averaged [mol P m-2 s-1]
    ! 1.2 conditions for phenology calculations
    TYPE(t_jsb_var_real2d)          :: &
                                    t_air_tphen_mavg               , &  !< air temperature time-averaged [K]
                                    t_soil_srf_tphen_mavg
    ! 1.3 memory for plant uptake (average demand over a couple of days to get rid of diurnal C cycle)
    TYPE(t_jsb_var_real2d)          :: &
                                    npp_tuptake_mavg               , &  !< npp of the last week [micro-mol m-2 s-1]
                                    demand_uptake_n_tuptake_mavg   , &  !< scalar of Vmax for root N uptake time-averaged [unitless]
                                    demand_uptake_p_tuptake_mavg   , &  !< scalar of Vmax for root P uptake time-averaged [unitless]
                                    growth_req_n_tuptake_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions time-averaged [mol N m-2 s-1]
                                    growth_req_p_tuptake_mavg      , &  !< moles P required for a unit of N growth under current allocation fractions time-averaged [mol P m-2 s-1]
                                    exudation_c_tmyc_mavg               !< moving average of C allocation to mycorrhizal fungi [mol C]
    ! 1.4 conditions for photosynthesis for the optimisation of foliar N fractions
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_tfrac_mavg_cl          !< fraction of layer which is sunlit time-averaged [unitless]
    TYPE(t_jsb_var_real2d)          :: &
                                    t_air_tfrac_mavg               , &  !< air temperature time-averaged [K]
                                    press_srf_tfrac_mavg           , &  !< air pressure at surface time-averaged [Pa]
                                    co2_mixing_ratio_tfrac_mavg    , &  !< co2 mixing ratio of the air time-averaged [co2 in ppm]
                                    ga_tfrac_mavg                  , &  !< aerodynamic conductance time-averaged [m s-1]
                                    beta_sinklim_ps_tfrac_mavg          !< scaling factor to account for soil moisture constraints on photosynthesis time-averaged [unitless]
    ! 1.5 conditions for photosynthesis for the optimisation of foliar N
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_tcnl_mavg_cl           !< fraction of layer which is sunlit time-averaged [unitless]
    TYPE(t_jsb_var_real2d)          :: &
                                    t_air_tcnl_mavg                , &  !< air temperature time-averaged [K]
                                    press_srf_tcnl_mavg            , &  !< air pressure at surface time-averaged [Pa]
                                    co2_mixing_ratio_tcnl_mavg     , &  !< co2 mixing ratio of the air time-averaged [co2 in ppm]
                                    ga_tcnl_mavg                   , &  !< aerodynamic conductance time-averaged [m s-1]
                                    uptake_n_tcnl_mavg             , &  !< root N uptake time-averaged [micro-mol N m-2 s-1]
                                    growth_cn_tcnl_mavg            , &  !< new growth CN ratio time-averaged [mol C mol-1 N]
                                    growth_np_tcnl_mavg            , &  !< new growth NP ratio time-averaged [mol N mol-1 P]
                                    npp_tcnl_mavg                  , &  !< NPP time-averaged [micro-mol CO2 m-2 s-1]
                                    fmaint_rate_tcnl_mavg          , &  !< maintenance respiration rate per unit N for roots
                                    labile_nitrogen_tcnl_mavg      , &  !< labile nitrogen vegetation pool time-averaged [mol N m-2]
                                    labile_carbon_tcnl_mavg        , &  !< labile carbon vegetation pool time-averaged [mol C m-2]
                                    labile_phosphorus_tcnl_mavg         !< labile phosphorus vegetation pool time-averaged [mol P m-2]
    ! 1.6 conditions for photosynthesis for the optimisation of biomass allocation fractions
    TYPE(t_jsb_var_real2d)          :: &
                                    npp_talloc_mavg                , &  !< NPP [micro-mol m-2 s-1]
                                    n_fixation_talloc_mavg         , &  !< N fixation time-averaged [micro-mol N m-2 s-1]
                                    transpiration_talloc_mavg      , &  !< transpiration time-averaged [kg m-2 s-1]
                                    unit_npp_talloc_mavg           , &  !< rate of change in NPP with  increment in leaf mass [micro-mol C mol-1]
                                    unit_transpiration_talloc_mavg , &  !< incremental change in transpiration time-averaged [kg m-2 s-1]
                                    unit_uptake_n_talloc_mavg      , &  !< actual rate of N uptkae per unit root [micro-mol N mol-1]
                                    unit_uptake_p_talloc_mavg      , &  !< actual rate of P uptake per unit root [micro-mol P mol-1]
                                    dphi_talloc_mavg               , &  !< leaf-soil water potential gradient [MPa]
                                    growth_cn_talloc_mavg          , &  !< plant growth cn ratio, equivalent to 1/growth_req_n time-averaged [mol C mol-1 N]
                                    growth_cp_talloc_mavg          , &  !< plant growth cp ratio, equivalent to 1/(growthreq_n*growth_req_p) time-averaged [mol C mol-1 P]
                                    growth_np_talloc_mavg          , &  !< plant growth np ratio, equivalent to 1/(growth_req_p) time-averaged [mol N mol-1 P]
                                    labile_carbon_talloc_mavg      , &  !< labile carbon pool needed for allometry calculations time-averaged [mol C m-2]
                                    labile_nitrogen_talloc_mavg    , &  !< labile nitrogen pool neeeded for allometry calculations time-averaged [mol N m-2]
                                    labile_phosphorus_talloc_mavg  , &  !< labile phosphorus pool neeeded for allometry calculations time-averaged [mol P m-2]
                                    growth_p_limit_based_on_n_mavg , &  !< indicator of how much plant is P limited based on the N status (plant growth p:n ratio divided by labile p:n ratio)
                                    reserve_carbon_talloc_mavg     , &  !< reserve carbon needed for optimal allocation time-averaged [mol C m-2]
                                    reserve_nitrogen_talloc_mavg   , &  !< reserve nitrogen needed for optimal allocation time-averaged [mol N m-2]
                                    reserve_phosphorus_talloc_mavg , &  !< reserve phosphorus needed for optimal allocation time-averaged [mol P m-2]
                                    w_root_lim_talloc_mavg              !< root water limitation below which root allocation starts responding time-averaged
    ! 1.7 memory over root life time
    TYPE(t_jsb_var_real2d)          :: &
                                    leaf2root_troot_mavg           , &  !< leaf to root mass ratio time-averaged [unitless]
                                    unit_npp_troot_mavg            , &  !< rate of change in NPP with increment in leaf mass [micro-mol C mol-1]
                                    unit_uptake_n_pot_troot_mavg   , &  !< rate of N uptkae per unit root [micro-mol N mol-1]
                                    unit_uptake_p_pot_troot_mavg   , &  !< rate of P uptake per unit root [micro-mol P mol-1]
                                    fmaint_rate_troot_mavg              !< moving average of N-specific maintentance respiration [~micro-mol C mol-1 N s-1]
    ! 1.8 memory over timescale of vegetation dynamics
    TYPE(t_jsb_var_real2d)          :: &
                                    an_boc_tvegdyn_mavg            , &  !< net photosynthesis at bottom of canopy, lowest canopy layer existing time-averaged [micro-mol m-2 s-1]
                                    net_growth_tvegdyn_mavg        , &  !< net growth, npp minus turnover [micro-mol C m-2 s-1]
                                    lai_tvegdyn_mavg                    !< lai time-averaged [m2 m-2]
    ! 1.9 variables for respiration accclimation
    TYPE(t_jsb_var_real2d)          :: &
                                    beta_sinklim_ps_tacclim_mavg   , &  !< scaling factor to adjust beta sinklimitation to respond to labile carbon with delay [unitless]
                                    t_air_tacclim_mavg             , &  !< air temperature for respiration acclimation time-averaged [K]
                                    t_soil_root_tacclim_mavg            !< root temperature for respiration acclimation time-averaged [deg K]
    ! 1.10 certain variables, temperature, ...
    TYPE(t_jsb_var_real2d)          :: &
                                    t_air_week_mavg                , &  !< air temperature time-averaged [K]
                                    t_air_month_mavg               , &  !< air temperature time-averaged [K]
                                    t_jmax_opt_mavg                , &  !< optimal temperature for electron transport of photosynthesis time-averaged [deg C]
                                    t_soil_root                         !< Temperature of fine roots in the soil [K]

    ! 1.11 helper variables for L2A_ process
    !      these var are used at PFT tiles and aggregated to box tile
    !      in 'mo_atmland_interface:update_land2atm' their values are added to the "L2A_ counterpart variables"
    !      this is needed because the L2A_ var are not available at other tiles than the box tile (i.e., top tile)
    TYPE(t_jsb_var_real2d) ::       &
      & net_biosphere_production  , & !< Balance between carbon losses and gains [micro-mol CO2 m-2 s-1]
      & biological_n_fixation         !< = n_fixation (VEG_) + SUM(asymb_n_fixation (SB_), across soil_layers) [micro-mol N m-2 s-1]

  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory

  ! ======================================================================================================= !
  !>
  !> Type definition for veg_bgc_material_elements
  !>
  TYPE, EXTENDS(t_jsb_pool) :: t_veg_bgcm_with_elements
    TYPE(t_jsb_var_real2d), POINTER :: carbon
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen
    TYPE(t_jsb_var_real2d), POINTER :: phosphorus
    TYPE(t_jsb_var_real2d), POINTER :: carbon13
    TYPE(t_jsb_var_real2d), POINTER :: carbon14
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen15
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => veg_bgcm_with_elements_get_element_name_by_id
    PROCEDURE :: Init                       => veg_bgcm_with_elements_init
  END TYPE t_veg_bgcm_with_elements

  ! ======================================================================================================= !
  !>
  !> Type definition for veg_bgc_material_components
  !>
  TYPE, EXTENDS(t_jsb_pool) :: t_veg_bgcm_with_components
    TYPE(t_veg_bgcm_with_elements),    POINTER :: leaf
    TYPE(t_veg_bgcm_with_elements),    POINTER :: fine_root
    TYPE(t_veg_bgcm_with_elements),    POINTER :: coarse_root
    TYPE(t_veg_bgcm_with_elements),    POINTER :: sap_wood
    TYPE(t_veg_bgcm_with_elements),    POINTER :: heart_wood
    TYPE(t_veg_bgcm_with_elements),    POINTER :: labile
    TYPE(t_veg_bgcm_with_elements),    POINTER :: reserve
    TYPE(t_veg_bgcm_with_elements),    POINTER :: fruit
    TYPE(t_veg_bgcm_with_elements),    POINTER :: seed_bed
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => veg_bgcm_with_components_get_element_name_by_id
    PROCEDURE :: Init                       => veg_bgcm_with_components_init
  END TYPE t_veg_bgcm_with_components

  ! ======================================================================================================= !
  !>
  !> Type definition for veg_bgc_material_main
  !> the main bgc_material structure containing all VEG_ bgc_material structures
  !> main structure, added directly into the "TYPE(t_jsb_memory) :: t_veg_memory" as "CLASS(t_jsb_pool) :: bgc_material"
  !>
  TYPE, EXTENDS(t_jsb_pool) :: t_veg_bgcm_main
    !! pools
    TYPE(t_veg_bgcm_with_components),  POINTER :: vegbpool
    TYPE(t_veg_bgcm_with_elements),    POINTER :: vegbtotal_biomass
    ! ... and if with product pools
    TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_pp_fuel
    TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_pp_paper
    TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_pp_fiberboard
    TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_pp_oirw
    TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_pp_pv
    TYPE(t_veg_bgcm_with_elements),    POINTER :: veg_pp_sawnwood
    !! fluxes
    TYPE(t_veg_bgcm_with_components),  POINTER :: vegbgrowth
    TYPE(t_veg_bgcm_with_components),  POINTER :: vegblitterfall
    TYPE(t_veg_bgcm_with_components),  POINTER :: vegbfrac_alloc
    TYPE(t_veg_bgcm_with_elements),    POINTER :: vegbexudation
    TYPE(t_veg_bgcm_with_elements),    POINTER :: vegbestablishment
    TYPE(t_veg_bgcm_with_elements),    POINTER :: vegbreserve_use
    ! ... and if with product pools
    TYPE(t_veg_bgcm_with_elements),    POINTER :: fprod_decay
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => veg_bgcm_main_get_element_name_by_id
    PROCEDURE :: Init                       => veg_bgcm_main_init
  END TYPE t_veg_bgcm_main

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the VEG_ process
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,         ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,     ONLY: t_jsb_model
    USE mo_jsb_class,           ONLY: Get_model
    USE mo_quincy_output_class, ONLY: unitless
    dsl4jsb_Use_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_veg_memory),  INTENT(inout), TARGET :: mem             !< veg memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),                POINTER :: model                      !< model
    TYPE(t_jsb_grid),                 POINTER :: hgrid                      !< Horizontal grid
    TYPE(t_jsb_vgrid),                POINTER :: surface                    !< Vertical grid
    TYPE(t_jsb_vgrid),                POINTER :: vgrid_canopy_q_assimi      !< Vertical grid
    TYPE(t_jsb_vgrid),                POINTER :: vgrid_soil_sb              !< Vertical grid
    CHARACTER(len=30)                         :: unit_veg_pool              !< unit of element var
    CHARACTER(len=30)                         :: unit_veg_flux              !< unit of element var
    INTEGER                                   :: elem_idx_map(LAST_ELEM_ID) !< element mapper ID -> IND
    INTEGER                                   :: table                      !< ...
    TYPE(t_grib2)                             :: grib2_desc                 !< used for lai settings
    CHARACTER(len=*), PARAMETER               :: routine = TRIM(modname)//':Init_veg_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(model_id)
    table                 =  tables(1)
    hgrid                 => Get_grid(model%grid_id)
    surface               => Get_vgrid('surface')
    vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
    vgrid_soil_sb         => Get_vgrid('soil_layer_sb')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(VEG_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! init local variables
    unit_veg_pool = "mol m-2"
    unit_veg_flux = "mol m-2 timestep-1"

    ! ----------------------------------------------------------------------------------------------------- !
    ! add memory only for LAND & PFT tiles
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      ! --------------------------------------------------------------------------------------------------- !
      !>  1.0 bgc_material structures
      !>
      !>  @TODO replace "b" in bgcm object names by "_" (vegbpool -> veg_pool) !!
      !>
      ! copy elements_index_map from config to local
      elem_idx_map(:) = model%config%elements_index_map(:)
      ! --------------------------------------------------------------------------------------------------- !
      !> create and populate the main veg bgc_material structure mem%bgc_material that contains all other bgc_material
      !> note: no spaces in 'name' and/or 'shortname'
      !>
      ALLOCATE(t_veg_bgcm_main :: mem%bgc_material)
      CALL mem%bgc_material%Init(name = 'vegetation_bgc_material_main', shortname = 'veg_bgcm_main', &
        &                        elements_index_map = elem_idx_map(:), l_elements = .FALSE.)
      SELECT TYPE (selected => mem%bgc_material)
      CLASS IS (t_veg_bgcm_main)
        mem%veg_bgcm_main => selected
      END SELECT

      ! --------------------------------------------------------------------------------------------------- !
      !> Add the info to this memory that this memory contains bgc materials
      !>
      mem%has_bgc_materials = .TRUE.

      ! --------------------------------------------------------------------------------------------------- !
      !> create the bgc_material components structures
      !>
      ! each component structure is created within the main bgcm pool structure together with required bookkeeping information
      ! "shortname" is used for model output names
      ! veg_pool
      mem%veg_bgcm_main%vegbpool => Add_veg_bgc_material_with_components(mem%veg_bgcm_main, &
        & VEG_BGCM_POOL_ID, name = 'vegetation_bgcm_pool', shortname = 'veg_bgcm_pool',     &
        & unit = unit_veg_pool, elements_index_map = elem_idx_map(:))
      ! veg_growth
      mem%veg_bgcm_main%vegbgrowth => Add_veg_bgc_material_with_components(mem%veg_bgcm_main, &
        & VEG_BGCM_GROWTH_ID, name = 'vegetation_bgcm_growth', shortname = 'veg_bgcm_growth', &
        & unit = unit_veg_flux, elements_index_map = elem_idx_map(:))
      ! veg_litterfall
      mem%veg_bgcm_main%vegblitterfall => Add_veg_bgc_material_with_components(mem%veg_bgcm_main,         &
        & VEG_BGCM_LITTERFALL_ID, name = 'vegetation_bgcm_litterfall', shortname = 'veg_bgcm_litterfall', &
        & unit = unit_veg_flux, elements_index_map = elem_idx_map(:))
      ! veg_frac_alloc
      mem%veg_bgcm_main%vegbfrac_alloc => Add_veg_bgc_material_with_components(mem%veg_bgcm_main,         &
        & VEG_BGCM_FRAC_ALLOC_ID, name = 'vegetation_bgcm_frac_alloc', shortname = 'veg_bgcm_frac_alloc', &
        & unit = unit_veg_flux, elements_index_map = elem_idx_map(:))

      ! --------------------------------------------------------------------------------------------------- !
      !> create the bgc_material elements structures
      !>
      ! also here "shortname" is used for model output names
      ! veg_total_biomass
      mem%veg_bgcm_main%vegbtotal_biomass => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,    &
        & name = 'vegetation_bgcm_total_biomass', shortname = 'veg_bgcm_total_biomass',               &
        & elements_index_map = elem_idx_map(:), unit = unit_veg_pool, bgcm_id = VEG_BGCM_TOTAL_BIO_ID)
      ! veg_exudation
      mem%veg_bgcm_main%vegbexudation => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,        &
        & name = 'vegetation_bgcm_exudation', shortname = 'veg_bgcm_exudation',                       &
        & elements_index_map = elem_idx_map(:), unit = unit_veg_flux, bgcm_id = VEG_BGCM_EXUDATION_ID)
      ! veg_establishment
      mem%veg_bgcm_main%vegbestablishment => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,        &
        & name = 'vegetation_bgcm_establishment', shortname = 'veg_bgcm_establishment',                   &
        & elements_index_map = elem_idx_map(:), unit = unit_veg_flux, bgcm_id = VEG_BGCM_ESTABLISHMENT_ID)
      ! veg_reserve_use
      mem%veg_bgcm_main%vegbreserve_use => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,        &
        & name = 'vegetation_bgcm_reserve_use', shortname = 'veg_bgcm_reserve_use',                     &
        & elements_index_map = elem_idx_map(:), unit = unit_veg_flux, bgcm_id = VEG_BGCM_RESERVE_USE_ID)

      IF(dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
        ! 6 product pools represented as bgcms with elements
        ! the fuelwood product pool
        mem%veg_bgcm_main%veg_pp_fuel => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,         &
          & name = 'vegetation_bgcm_product_pool_fuelwood', shortname = 'veg_bgcm_pp_fuel',            &
          & elements_index_map = elem_idx_map(:), unit = unit_veg_pool, bgcm_id = VEG_BGCM_PP_FUEL_ID)
        ! the paper product pool
        mem%veg_bgcm_main%veg_pp_paper => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,         &
          & name = 'vegetation_bgcm_product_pool_paper', shortname = 'veg_bgcm_pp_paper',               &
          & elements_index_map = elem_idx_map(:), unit = unit_veg_pool, bgcm_id = VEG_BGCM_PP_PAPER_ID)
        ! the fiberboard product pool
        mem%veg_bgcm_main%veg_pp_fiberboard => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,         &
          & name = 'vegetation_bgcm_product_pool_fiberboard', shortname = 'veg_bgcm_pp_fiberboard',          &
          & elements_index_map = elem_idx_map(:), unit = unit_veg_pool, bgcm_id = VEG_BGCM_PP_FIBERBOARD_ID)
        ! the other industrial roundwood product pool
        mem%veg_bgcm_main%veg_pp_oirw => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,         &
          & name = 'vegetation_bgcm_product_pool_other_irw', shortname = 'veg_bgcm_pp_oirw',           &
          & elements_index_map = elem_idx_map(:), unit = unit_veg_pool, bgcm_id = VEG_BGCM_PP_OIRW_ID)
        ! the plywood and veneer product pool
        mem%veg_bgcm_main%veg_pp_pv => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,         &
          & name = 'vegetation_bgcm_product_pool_plywood_and_veneer', shortname = 'veg_bgcm_pp_pv',  &
          & elements_index_map = elem_idx_map(:), unit = unit_veg_pool, bgcm_id = VEG_BGCM_PP_PV_ID)
        ! the sawnwood product pool
        mem%veg_bgcm_main%veg_pp_sawnwood => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,         &
          & name = 'vegetation_bgcm_product_pool_sawnwood', shortname = 'veg_bgcm_pp_sawnwood',            &
          & elements_index_map = elem_idx_map(:), unit = unit_veg_pool, bgcm_id = VEG_BGCM_PP_SAWNWOOD_ID)
        ! and one flux represented as bgcm with elements: the product pool decay flux
        mem%veg_bgcm_main%fprod_decay => Add_veg_bgc_material_with_elements(mem%veg_bgcm_main,         &
          & name = 'vegetation_bgcm_flux_product_pool_decay', shortname = 'veg_bgcm_fprod_decay',      &
          & elements_index_map = elem_idx_map(:), unit = unit_veg_flux, bgcm_id = VEG_BGCM_FPROD_DECAY_ID)
      ENDIF

      ! create list of CHARACTERs with names veg_bgcm_main to the "lowest childs" of bgcm structures
      CALL mem%bgc_material%Set_paths()

      ! add variables to bgc_material
      CALL mem%Add_var(mem%veg_bgcm_main%vegbpool, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .TRUE.,                                                                           &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%veg_bgcm_main%vegbgrowth, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%veg_bgcm_main%vegblitterfall, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%veg_bgcm_main%vegbfrac_alloc, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%veg_bgcm_main%vegbtotal_biomass, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%veg_bgcm_main%vegbexudation, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%veg_bgcm_main%vegbestablishment, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%veg_bgcm_main%vegbreserve_use, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      ! Product pools and decay
      IF(dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
        CALL mem%Add_var(mem%veg_bgcm_main%veg_pp_fuel, &
          & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
          & prefix, suffix,                                                                              &
          & output_level = MEDIUM,                                                                       &
          & loutput = .TRUE.,                                                                            &
          & lrestart = .TRUE.,                                                                           &
          & initval_r = 0.0_wp)
        CALL mem%Add_var(mem%veg_bgcm_main%veg_pp_paper, &
          & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
          & prefix, suffix,                                                                              &
          & output_level = MEDIUM,                                                                       &
          & loutput = .TRUE.,                                                                            &
          & lrestart = .TRUE.,                                                                           &
          & initval_r = 0.0_wp)
        CALL mem%Add_var(mem%veg_bgcm_main%veg_pp_fiberboard, &
          & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
          & prefix, suffix,                                                                              &
          & output_level = MEDIUM,                                                                       &
          & loutput = .TRUE.,                                                                            &
          & lrestart = .TRUE.,                                                                           &
          & initval_r = 0.0_wp)
        CALL mem%Add_var(mem%veg_bgcm_main%veg_pp_oirw, &
          & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
          & prefix, suffix,                                                                              &
          & output_level = MEDIUM,                                                                       &
          & loutput = .TRUE.,                                                                            &
          & lrestart = .TRUE.,                                                                           &
          & initval_r = 0.0_wp)
        CALL mem%Add_var(mem%veg_bgcm_main%veg_pp_pv, &
          & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
          & prefix, suffix,                                                                              &
          & output_level = MEDIUM,                                                                       &
          & loutput = .TRUE.,                                                                            &
          & lrestart = .TRUE.,                                                                           &
          & initval_r = 0.0_wp)
        CALL mem%Add_var(mem%veg_bgcm_main%veg_pp_sawnwood, &
          & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
          & prefix, suffix,                                                                              &
          & output_level = MEDIUM,                                                                       &
          & loutput = .TRUE.,                                                                            &
          & lrestart = .TRUE.,                                                                           &
          & initval_r = 0.0_wp)
        CALL mem%Add_var(mem%veg_bgcm_main%fprod_decay, &
          & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
          & prefix, suffix,                                                                              &
          & output_level = MEDIUM,                                                                       &
          & loutput = .TRUE.,                                                                            &
          & lrestart = .FALSE.,                                                                          &
          & initval_r = 0.0_wp)
      ENDIF
      ! --------------------------------------------------------------------------------------------------- !
      !> 2.0 bgc_material diagnostics variables - used for model output with IQ
      !>
      CALL mem%Add_var('veg_pool_total_c', mem%veg_pool_total_c, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_c', 'mol m-2', 'total C, sum of all veg_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_total_n', mem%veg_pool_total_n, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_n', 'mol m-2', 'total N, sum of all veg_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_total_p', mem%veg_pool_total_p, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_p', 'mol m-2', 'total P, sum of all veg_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_total_c13', mem%veg_pool_total_c13, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_c13', 'mol m-2', 'total C13, sum of all veg_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_total_c14', mem%veg_pool_total_c14, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_c14', 'mol m-2', 'total C14, sum of all veg_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_total_n15', mem%veg_pool_total_n15, &
        & hgrid, surface, &
        & t_cf('veg_pool_total_n15', 'mol m-2', 'total N15, sum of all veg_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      ! diagnostic variables specific to QUINCY pools in ICON-Land
      CALL mem%Add_var('veg_growth_total_c', mem%veg_growth_total_c, &
        & hgrid, surface, &
        & t_cf('veg_growth_total_c', 'mol m-2 timestep-1', 'total C, sum of all veg_growth components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_growth_total_n', mem%veg_growth_total_n, &
        & hgrid, surface, &
        & t_cf('veg_growth_total_n', 'mol m-2 timestep-1', 'total N, sum of all veg_growth components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_growth_total_p', mem%veg_growth_total_p, &
        & hgrid, surface, &
        & t_cf('veg_growth_total_p', 'mol m-2 timestep-1', 'total P, sum of all veg_growth components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_growth_total_c13', mem%veg_growth_total_c13, &
        & hgrid, surface, &
        & t_cf('veg_growth_total_c13', 'mol m-2 timestep-1', 'total C13, sum of all veg_growth components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_growth_total_c14', mem%veg_growth_total_c14, &
        & hgrid, surface, &
        & t_cf('veg_growth_total_c14', 'mol m-2 timestep-1', 'total C14, sum of all veg_growth components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_growth_total_n15', mem%veg_growth_total_n15, &
        & hgrid, surface, &
        & t_cf('veg_growth_total_n15', 'mol m-2 timestep-1', 'total N15, sum of all veg_growth components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      ! diagnostic variables specific to QUINCY pools in ICON-Land
      CALL mem%Add_var('veg_litterfall_total_c', mem%veg_litterfall_total_c, &
        & hgrid, surface, &
        & t_cf('veg_litterfall_total_c', 'mol m-2 timestep-1', 'total C, sum of all veg_litterfall components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_litterfall_total_n', mem%veg_litterfall_total_n, &
        & hgrid, surface, &
        & t_cf('veg_litterfall_total_n', 'mol m-2 timestep-1', 'total N, sum of all veg_litterfall components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_litterfall_total_p', mem%veg_litterfall_total_p, &
        & hgrid, surface, &
        & t_cf('veg_litterfall_total_p', 'mol m-2 timestep-1', 'total P, sum of all veg_litterfall components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_litterfall_total_c13', mem%veg_litterfall_total_c13, &
        & hgrid, surface, &
        & t_cf('veg_litterfall_total_c13', 'mol m-2 timestep-1', 'total C13, sum of all veg_litterfall components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_litterfall_total_c14', mem%veg_litterfall_total_c14, &
        & hgrid, surface, &
        & t_cf('veg_litterfall_total_c14', 'mol m-2 timestep-1', 'total C14, sum of all veg_litterfall components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_litterfall_total_n15', mem%veg_litterfall_total_n15, &
        & hgrid, surface, &
        & t_cf('veg_litterfall_total_n15', 'mol m-2 timestep-1', 'total N15, sum of all veg_litterfall components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_leaf_c', mem%veg_pool_leaf_c, &
        & hgrid, surface, &
        & t_cf('veg_pool_leaf_c', 'mol m-2', 'vegetation pool leaf carbon'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_leaf_n', mem%veg_pool_leaf_n, &
        & hgrid, surface, &
        & t_cf('veg_pool_leaf_n', 'mol m-2', 'vegetation pool leaf nitrogen'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_wood_c', mem%veg_pool_wood_c, &
        & hgrid, surface, &
        & t_cf('veg_pool_wood_c', 'mol m-2', 'vegetation pool sapwood + heartwood carbon'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_wood_n', mem%veg_pool_wood_n, &
        & hgrid, surface, &
        & t_cf('veg_pool_wood_n', 'mol m-2', 'vegetation pool sapwood + heartwood nitrogen'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_fine_root_c', mem%veg_pool_fine_root_c, &
        & hgrid, surface, &
        & t_cf('veg_pool_fine_root_c', 'mol m-2', 'vegetation pool fine root carbon'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('veg_pool_fine_root_n', mem%veg_pool_fine_root_n, &
        & hgrid, surface, &
        & t_cf('veg_pool_fine_root_n', 'mol m-2', 'vegetation pool fine root nitrogen'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      IF(dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
        CALL mem%Add_var('veg_products_total_c', mem%veg_products_total_c, &
          & hgrid, surface, &
          & t_cf('veg_products_total_c', 'mol m-2', 'total c in product pools'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = MEDIUM, &
          & loutput = .TRUE., &
          & lrestart = .FALSE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('veg_products_total_n', mem%veg_products_total_n, &
          & hgrid, surface, &
          & t_cf('veg_products_total_n', 'mol m-2', 'total n in product pools'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = MEDIUM, &
          & loutput = .TRUE., &
          & lrestart = .FALSE., &
          & initval_r = 0.0_wp)
      END IF

      ! --------------------------------------------------------------------------------------------------- !
      !> 3.0 2D & 3D variables
      !>
      CALL mem%Add_var('height', mem%height, &
        & hgrid, surface, &
        & t_cf('height', 'm', 'plant height'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('diameter', mem%diameter, &
        & hgrid, surface, &
        & t_cf('diameter', 'm', 'plant diameter'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dens_ind', mem%dens_ind, &
        & hgrid, surface, &
        & t_cf('dens_ind', '# m-2', 'density of individuals'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      ! grib setting for: lai differentiation box tile and other tiles
      ! TODO this is done because ... (ask Reiner Schnur, MPI-M)
      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(2,0,28, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF

      CALL mem%Add_var('lai', mem%lai, &
        & hgrid, surface, &
        & t_cf('lai', 'm2 m-2', 'actual leaf area index'), &
        & t_grib1(table, 255, grib_bits), &
        & grib2_desc, &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sai', mem%sai, &
        & hgrid, surface, &
        & t_cf('sai', 'm2 m-2', 'actual stem area index'), &
        & t_grib1(table, 255, grib_bits), &
        & grib2_desc, &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_cn_leaf', mem%target_cn_leaf, &
        & hgrid, surface, &
        & t_cf('target_cn_leaf', '', 'target foliar C:N ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_np_leaf', mem%target_np_leaf, &
        & hgrid, surface, &
        & t_cf('target_np_leaf', '', 'target foliar N:P ratio '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('mean_leaf_age', mem%mean_leaf_age, &
        & hgrid, surface, &
        & t_cf('mean_leaf_age', '', 'average foliage age in days '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fract_fpc', mem%fract_fpc, &
        & hgrid, surface, &
        & t_cf('fract_fpc', '', 'foliage projected cover fraction'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('blended_height ', mem%blended_height, &
        & hgrid, surface, &
        & t_cf('blended_height', '', 'foliage weighted height'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('cohort_age', mem%cohort_age, &
        & hgrid, surface, &
        & t_cf('cohort_age', '', 'age of the cohort in years, used with veg_dynamics_scheme cohort'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('f_n_demand', mem%f_n_demand, &
        & hgrid, surface, &
        & t_cf('f_n_demand', '', 'relative plant demand for N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('f_p_demand', mem%f_p_demand, &
        & hgrid, surface, &
        & t_cf('f_p_demand', '', 'relative plant demand for P'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('beta_sinklim_ps', mem%beta_sinklim_ps, &
        & hgrid, surface, &
        & t_cf('beta_sinklim_ps', '', 'scaling factor to account for a direct sink limitation contraint on photosynthesis'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dphi', mem%dphi, &
        & hgrid, surface, &
        & t_cf('dphi', 'MPa', 'plant-soil water potential gradient'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_req_n', mem%growth_req_n, &
        & hgrid, surface, &
        & t_cf('growth_req_n', '', 'moles N required for a unit of C growth under current allocation fractions '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_req_p', mem%growth_req_p, &
        & hgrid, surface, &
        & t_cf('growth_req_p', '', 'moles P required for a unit of N growth under current allocation fractions '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaf_cn_direction', mem%leaf_cn_direction, &
        & hgrid, surface, &
        & t_cf('leaf_cn_direction', '', 'indicates direction of change in leaf CN ratio to maximise NPP'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_lai', mem%target_lai, &
        & hgrid, surface, &
        & t_cf('target_lai', 'm2 m-2', 'LAI implied from sap wood area'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('kstar_labile', mem%kstar_labile, &
        & hgrid, surface, &
        & t_cf('kstar_labile', '', 'labile pool turnover rate, only temperature and moisture constraints, no stoichiometry'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('root_limitation_state', mem%root_limitation_state, &
        & hgrid, surface, &
        & t_cf('root_limitation_state', '', 'root_limitation_state'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('do_cohort_harvest', mem%do_cohort_harvest, &
        & hgrid, surface, &
        & t_cf('do_cohort_harvest', '', 'logical: cohort harvest on/off, in veg_dynamics_scheme cohort'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('root_fraction_sl', mem%root_fraction_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('root_fraction_sl', '', 'fractions of root per layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('delta_root_fraction_sl', mem%delta_root_fraction_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('delta_root_fraction_sl', '', 'change in fractions of root per layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lai_cl', mem%lai_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('lai_cl', 'm2 m-2', 'leaf area index per canopy layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('cumm_lai_cl', mem%cumm_lai_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('cumm_lai_cl', 'm2 m-2', 'cummulative leaf area above the mid-point of the canopy layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_cn_fine_root', mem%target_cn_fine_root, &
        & hgrid, surface, &
        & t_cf('target_cn_fine_root', 'mol/mol', 'target C:N ratio for fine roots'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_cn_coarse_root', mem%target_cn_coarse_root, &
        & hgrid, surface, &
        & t_cf('target_cn_coarse_root', 'mol/mol', 'target C:N ratio for coarse roots'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_cn_sap_wood', mem%target_cn_sap_wood, &
        & hgrid, surface, &
        & t_cf('target_cn_sap_wood', 'mol/mol', 'target C:N ratio for sap wood '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_cn_fruit', mem%target_cn_fruit, &
        & hgrid, surface, &
        & t_cf('target_cn_fruit', 'mol/mol', 'target C:N ratio for fruit pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_np_fine_root', mem%target_np_fine_root, &
        & hgrid, surface, &
        & t_cf('target_np_fine_root', 'mol/mol', 'target N:P ratio for fine roots'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_np_coarse_root', mem%target_np_coarse_root, &
        & hgrid, surface, &
        & t_cf('target_np_coarse_root', 'mol/mol', 'target N:P ratio for coarse roots'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_np_sap_wood', mem%target_np_sap_wood, &
        & hgrid, surface, &
        & t_cf('target_np_sap_wood', 'mol/mol', 'target N:P ratio for sap wood pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('target_np_fruit', mem%target_np_fruit, &
        & hgrid, surface, &
        & t_cf('target_np_fruit', 'mol/mol', 'target N:P ratio for fruit pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaf_nitrogen_cl', mem%leaf_nitrogen_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('leaf_nitrogen_cl', 'mmol m-2 LAI', 'N in canopy layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fn_chl_cl', mem%fn_chl_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fn_chl_cl', unitless, 'fraction of N in chlorophyll, pigments only'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fn_et_cl', mem%fn_et_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fn_et_cl', unitless, 'fraction of N in electron transport '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fn_rub_cl', mem%fn_rub_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fn_rub_cl', unitless, 'fraction of N in Rubisco '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fn_pepc_cl', mem%fn_pepc_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fn_pepc_cl', unitless, 'fraction of N in PEP carboylase, C4 plants only'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fn_oth_cl', mem%fn_oth_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fn_oth_cl', unitless, 'fraction of N not associated with photosynthesis '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fleaf_sunlit_cl', mem%fleaf_sunlit_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fleaf_sunlit_cl', unitless, 'fraction of layer which is sunlit'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('maint_respiration_pot', mem%maint_respiration_pot, &
        & hgrid, surface, &
        & t_cf('maint_respiration_pot', 'micro-mol CO2 m-2 s-1', 'potential maintenance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('maint_respiration', mem%maint_respiration, &
        & hgrid, surface, &
        & t_cf('maint_respiration', 'micro-mol CO2 m-2 s-1', 'maintenance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('maint_respiration_c13', mem%maint_respiration_c13, &
        & hgrid, surface, &
        & t_cf('maint_respiration_c13', 'micro-mol 13CO2 m-2 s-1', 'maintenance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('maint_respiration_c14', mem%maint_respiration_c14, &
        & hgrid, surface, &
        & t_cf('maint_respiration_c14', 'micro-mol 14CO2 m-2 s-1', 'maintenance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_respiration', mem%growth_respiration, &
        & hgrid, surface, &
        & t_cf('growth_respiration', 'micro-mol CO2 m-2 s-1', 'growth respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_respiration_c13', mem%growth_respiration_c13, &
        & hgrid, surface, &
        & t_cf('growth_respiration_c13', 'micro-mol 13CO2 m-2 s-1', 'growth respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_respiration_c14', mem%growth_respiration_c14, &
        & hgrid, surface, &
        & t_cf('growth_respiration_c14', 'micro-mol 14CO2 m-2 s-1', 'growth respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n_transform_respiration', mem%n_transform_respiration, &
        & hgrid, surface, &
        & t_cf('n_transform_respiration', 'micro-mol CO2 m-2 s-1', 'respiration associated with N transformation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n_fixation_respiration', mem%n_fixation_respiration, &
        & hgrid, surface, &
        & t_cf('n_fixation_respiration', 'micro-mol CO2 m-2 s-1', 'respiration associated with N fixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n_processing_respiration', mem%n_processing_respiration, &
        & hgrid, surface, &
        & t_cf('n_processing_respiration', 'micro-mol CO2 m-2 s-1', 'respiration associated with N processing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n_processing_respiration_c13', mem%n_processing_respiration_c13, &
        & hgrid, surface, &
        & t_cf('n_processing_respiration_c13', 'micro-mol 13CO2 m-2 s-1', 'respiration associated with N processing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n_processing_respiration_c14', mem%n_processing_respiration_c14, &
        & hgrid, surface, &
        & t_cf('n_processing_respiration_c14', 'micro-mol 14CO2 m-2 s-1', 'respiration associated with N processing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('npp', mem%npp, &
        & hgrid, surface, &
        & t_cf('npp', 'micro-mol CO2 m-2 s-1', 'net primary production'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('npp_c13', mem%npp_c13, &
        & hgrid, surface, &
        & t_cf('npp_c13', 'micro-mol 13CO2 m-2 s-1', 'net primary production'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('npp_c14', mem%npp_c14, &
        & hgrid, surface, &
        & t_cf('npp_c14', 'micro-mol 14CO2 m-2 s-1', 'net primary production'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('uptake_nh4', mem%uptake_nh4, &
        & hgrid, surface, &
        & t_cf('uptake_nh4', 'micro-mol N m-2 s-1', 'NH4 uptake from soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('uptake_nh4_n15', mem%uptake_nh4_n15, &
        & hgrid, surface, &
        & t_cf('uptake_nh4_n15', 'micro-mol N m-2 s-1', '15NH4 uptake from soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('uptake_no3', mem%uptake_no3, &
        & hgrid, surface, &
        & t_cf('uptake_no3', 'micro-mol N m-2 s-1', 'NO3 uptake from soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('uptake_no3_n15', mem%uptake_no3_n15, &
        & hgrid, surface, &
        & t_cf('uptake_no3_n15', 'micro-mol N m-2 s-1', '15NO3 uptake from soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('uptake_po4', mem%uptake_po4, &
        & hgrid, surface, &
        & t_cf('uptake_po4', 'micro-mol P m-2 s-1', 'PO4 uptake from soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n_fixation', mem%n_fixation, &
        & hgrid, surface, &
        & t_cf('n_fixation', 'micro-mol N m-2 s-1', 'symbiontic N fixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = MEDIUM, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n_fixation_n15', mem%n_fixation_n15, &
        & hgrid, surface, &
        & t_cf('n_fixation_n15', 'micro-mol 15N m-2 s-1', 'symbiontic N fixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('recycling_leaf_n', mem%recycling_leaf_n, &
        & hgrid, surface, &
        & t_cf('recycling_leaf_n', 'mol m-2 timestep-1', 'net flux of N from leaf to labile (senescence and maintenance)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('recycling_leaf_n15', mem%recycling_leaf_n15, &
        & hgrid, surface, &
        & t_cf('recycling_leaf_n15', 'mol m-2 timestep-1', 'net flux of 15N from leaf to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('recycling_leaf_p', mem%recycling_leaf_p, &
        & hgrid, surface, &
        & t_cf('recycling_leaf_p', 'mol m-2 timestep-1', 'net flux of P from leaf to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('recycling_fine_root_n', mem%recycling_fine_root_n, &
        & hgrid, surface, &
        & t_cf('recycling_fine_root_n', 'mol m-2 timestep-1', 'net flux of N from fine roots to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('recycling_fine_root_n15', mem%recycling_fine_root_n15, &
        & hgrid, surface, &
        & t_cf('recycling_fine_root_n15', 'mol m-2 timestep-1', 'net flux of 15N from fine roots to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('recycling_fine_root_p', mem%recycling_fine_root_p, &
        & hgrid, surface, &
        & t_cf('recycling_fine_root_p', 'mol m-2 timestep-1', 'net flux of P from fine roots to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('recycling_heart_wood_n', mem%recycling_heart_wood_n, &
        & hgrid, surface, &
        & t_cf('recycling_heart_wood_n', 'mol m-2 timestep-1', 'net flux of N from sapwood to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('recycling_heart_wood_n15', mem%recycling_heart_wood_n15, &
        & hgrid, surface, &
        & t_cf('recycling_heart_wood_n15', 'mol m-2 timestep-1', 'net flux of 15N from sapwood to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('recycling_heart_wood_p', mem%recycling_heart_wood_p, &
        & hgrid, surface, &
        & t_cf('recycling_heart_wood_p', 'mol m-2 timestep-1', 'net flux of P from sapwood to labile, senescence & maint'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('net_growth', mem%net_growth, &
        & hgrid, surface, &
        & t_cf('net_growth', 'micro-mol CO2 m-2 s-1', 'npp minus turnover'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fmaint_rate_root', mem%fmaint_rate_root, &
        & hgrid, surface, &
        & t_cf('fmaint_rate_root', '~micro-mol C mol-1 N s-1', 'current N-specific maintentance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_npp', mem%unit_npp, &
        & hgrid, surface, &
        & t_cf('unit_npp', '', 'incremental change in NPP. Needed for optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_transpiration', mem%unit_transpiration, &
        & hgrid, surface, &
        & t_cf('unit_transpiration', '', 'incremental change in transpiration. Needed for optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_uptake_n_pot', mem%unit_uptake_n_pot, &
        & hgrid, surface, &
        & t_cf('unit_uptake_n_pot', '', 'incremental change in potential N uptake per unit root carbon, N uptake calculation '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_uptake_p_pot', mem%unit_uptake_p_pot, &
        & hgrid, surface, &
        & t_cf('unit_uptake_p_pot', '', 'incremental change in potential P uptake per unit root carbon, currently unused '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_uptake_n_act', mem%unit_uptake_n_act, &
        & hgrid, surface, &
        & t_cf('unit_uptake_n_act', '', 'incremental change in actual N uptake per unit root carbon, optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_uptake_p_act', mem%unit_uptake_p_act, &
        & hgrid, surface, &
        & t_cf('unit_uptake_p_act', '', 'incremental change in actual P uptake per unit root carbon, optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('delta_dens_ind', mem%delta_dens_ind, &
        & hgrid, surface, &
        & t_cf('delta_dens_ind', '# m-2 timestep-1', 'change in density of individuals'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('cost_n_uptake_root', mem%cost_n_uptake_root, &
        & hgrid, surface, &
        & t_cf('cost_n_uptake_root', 'mol C mol-1 N', 'actual cost for N uptake by roots'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('mortality_rate', mem%mortality_rate, &
        & hgrid, surface, &
        & t_cf('mortality_rate', 'timestep-1', 'background_mort_rate_tree / grass'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('k1_opt', mem%k1_opt, &
        & hgrid, surface, &
        & t_cf('k1_opt', unitless, 'coefficient for optimal leaf:root ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('k2_opt', mem%k2_opt, &
        & hgrid, surface, &
        & t_cf('k2_opt', unitless, 'coefficient for optimal leaf:root ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaf2sapwood_mass_ratio', mem%leaf2sapwood_mass_ratio, &
        & hgrid, surface, &
        & t_cf('leaf2sapwood_mass_ratio', unitless, 'leaf to shoot mass ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaf2root_mass_ratio', mem%leaf2root_mass_ratio, &
        & hgrid, surface, &
        & t_cf('leaf2root_mass_ratio', unitless, 'leaf to root mass ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fleaf_sunlit_daytime_cl', mem%fleaf_sunlit_daytime_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fleaf_sunlit_daytime_cl', unitless, 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_air_daytime', mem%t_air_daytime, &
        & hgrid, surface, &
        & t_cf('t_air_daytime', 'K', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('press_srf_daytime', mem%press_srf_daytime, &
        & hgrid, surface, &
        & t_cf('press_srf_daytime', 'Pa', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('co2_mixing_ratio_daytime', mem%co2_mixing_ratio_daytime, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio_daytime', 'co2 in ppm', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ga_daytime', mem%ga_daytime, &
        & hgrid, surface, &
        & t_cf('ga_daytime', 'co2 in ppm', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('beta_sinklim_ps_daytime', mem%beta_sinklim_ps_daytime, &
        & hgrid, surface, &
        & t_cf('beta_sinklim_ps_daytime', unitless, 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_jmax_opt_daytime', mem%t_jmax_opt_daytime, &
        & hgrid, surface, &
        & t_cf('t_jmax_opt_daytime', 'deg C', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fleaf_sunlit_daytime_dacc_cl', mem%fleaf_sunlit_daytime_dacc_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fleaf_sunlit_daytime_dacc_cl', unitless, 'value accumulated over current day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_air_daytime_dacc', mem%t_air_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('t_air_daytime_dacc', 'K', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('press_srf_daytime_dacc', mem%press_srf_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('press_srf_daytime_dacc', 'Pa', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('co2_mixing_ratio_daytime_dacc', mem%co2_mixing_ratio_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio_daytime_dacc', 'co2 in ppm', 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ga_daytime_dacc', mem%ga_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('ga_daytime_dacc', 'm s-1', 'aerodynamic conductance average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('beta_sinklim_ps_daytime_dacc', mem%beta_sinklim_ps_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('beta_sinklim_ps_daytime_dacc', unitless, 'average of the previous day'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_jmax_opt_daytime_dacc', mem%t_jmax_opt_daytime_dacc, &
        & hgrid, surface, &
        & t_cf('t_jmax_opt_daytime_dacc', 'deg C', 'value accumulated over current day '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('gpp_tlabile_mavg', mem%gpp_tlabile_mavg, &
        & hgrid, surface, &
        & t_cf('gpp_tlabile_mavg', 'micro-mol m-2 s-1', 'gpp '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('maint_respiration_tlabile_mavg', mem%maint_respiration_tlabile_mavg, &
        & hgrid, surface, &
        & t_cf('maint_respiration_tlabile_mavg', 'micro-mol m-2 s-1', 'maintenance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_req_n_tlabile_mavg', mem%growth_req_n_tlabile_mavg, &
        & hgrid, surface, &
        & t_cf('growth_req_n_tlabile_mavg', '', 'moles N required for a unit of C growth under current allocation fractions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_req_p_tlabile_mavg', mem%growth_req_p_tlabile_mavg, &
        & hgrid, surface, &
        & t_cf('growth_req_p_tlabile_mavg', '', 'moles P required for a unit of N growth under current allocation fractions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_air_tphen_mavg', mem%t_air_tphen_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_tphen_mavg', '', 'air temp mavg for phenology calc'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_soil_srf_tphen_mavg', mem%t_soil_srf_tphen_mavg, &
        & hgrid, surface, &
        & t_cf('t_soil_srf_tphen_mavg', '', 'soil temp mavg for phenology calc'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('npp_tuptake_mavg', mem%npp_tuptake_mavg, &
        & hgrid, surface, &
        & t_cf('npp_tuptake_mavg', 'micro-mol m-2 s-1', 'npp of the last week'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('demand_uptake_n_tuptake_mavg', mem%demand_uptake_n_tuptake_mavg, &
        & hgrid, surface, &
        & t_cf('demand_uptake_n_tuptake_mavg', '', 'scalar of Vmax for root N uptake '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('demand_uptake_p_tuptake_mavg', mem%demand_uptake_p_tuptake_mavg, &
        & hgrid, surface, &
        & t_cf('demand_uptake_p_tuptake_mavg', '', 'scalar of Vmax for root P uptake'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_req_n_tuptake_mavg', mem%growth_req_n_tuptake_mavg, &
        & hgrid, surface, &
        & t_cf('growth_req_n_tuptake_mavg', '', 'moles N required for a unit of C growth under current allocation fractions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_req_p_tuptake_mavg', mem%growth_req_p_tuptake_mavg, &
        & hgrid, surface, &
        & t_cf('growth_req_p_tuptake_mavg', '', 'moles P required for a unit of N growth under current allocation fractions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('exudation_c_tmyc_mavg', mem%exudation_c_tmyc_mavg, &
        & hgrid, surface, &
        & t_cf('exudation_c_tmyc_mavg', 'mol C', 'moving average of C allocation to mycorrhizae'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fleaf_sunlit_tfrac_mavg_cl', mem%fleaf_sunlit_tfrac_mavg_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fleaf_sunlit_tfrac_mavg_cl', '', 'fraction of sunlit leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_air_tfrac_mavg', mem%t_air_tfrac_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_tfrac_mavg', '', 'air temp'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('press_srf_tfrac_mavg', mem%press_srf_tfrac_mavg, &
        & hgrid, surface, &
        & t_cf('press_srf_tfrac_mavg', '', 'air pressure at surface'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('co2_mixing_ratio_tfrac_mavg', mem%co2_mixing_ratio_tfrac_mavg, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio_tfrac_mavg', '', 'co2 mixing ratio of the air'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ga_tfrac_mavg', mem%ga_tfrac_mavg, &
        & hgrid, surface, &
        & t_cf('ga_tfrac_mavg', '', 'aerodynamic conductance'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('beta_sinklim_ps_tfrac_mavg', mem%beta_sinklim_ps_tfrac_mavg, &
        & hgrid, surface, &
        & t_cf('beta_sinklim_ps_tfrac_mavg', '', 'scaling factor to account for soil moisture constraints on photosynthesis '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fleaf_sunlit_tcnl_mavg_cl', mem%fleaf_sunlit_tcnl_mavg_cl, &
        & hgrid, vgrid_canopy_q_assimi, &
        & t_cf('fleaf_sunlit_tcnl_mavg_cl', '', 'fraction of sunlit leaves'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_air_tcnl_mavg', mem%t_air_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_tcnl_mavg', '', 'air temp'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('press_srf_tcnl_mavg', mem%press_srf_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('press_srf_tcnl_mavg', '', 'air pressure at surface'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('co2_mixing_ratio_tcnl_mavg', mem%co2_mixing_ratio_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('co2_mixing_ratio_tcnl_mavg', '', 'co2 mixing ratio of the air'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ga_tcnl_mavg', mem%ga_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('ga_tcnl_mavg', '', 'aerodynamic conductance'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('uptake_n_tcnl_mavg', mem%uptake_n_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('uptake_n_tcnl_mavg', '', 'root N uptake'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_cn_tcnl_mavg', mem%growth_cn_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('growth_cn_tcnl_mavg', '', 'new growth CN ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_np_tcnl_mavg', mem%growth_np_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('growth_np_tcnl_mavg', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('npp_tcnl_mavg', mem%npp_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('npp_tcnl_mavg', '', 'NPP'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fmaint_rate_tcnl_mavg', mem%fmaint_rate_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('fmaint_rate_tcnl_mavg', '~micro-mol C mol-1 N s-1', 'N-specific maintentance respiration for optmal leaf N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('labile_nitrogen_tcnl_mavg', mem%labile_nitrogen_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('labile_nitrogen_tcnl_mavg', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('labile_carbon_tcnl_mavg', mem%labile_carbon_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('labile_carbon_tcnl_mavg', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('labile_phosphorus_tcnl_mavg', mem%labile_phosphorus_tcnl_mavg, &
        & hgrid, surface, &
        & t_cf('labile_phosphorus_tcnl_mavg', '', ' '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)


      CALL mem%Add_var('npp_talloc_mavg', mem%npp_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('npp_talloc_mavg', 'micro-mol m-2 s-1', 'NPP'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n_fixation_talloc_mavg', mem%n_fixation_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('n_fixation_talloc_mavg', '', 'N fixation averaged at the allocation timescale'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transpiration_talloc_mavg', mem%transpiration_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('transpiration_talloc_mavg', '', 'average transpiration at allocation time-scale'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_npp_talloc_mavg', mem%unit_npp_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('unit_npp_talloc_mavg', 'micro-mol C mol-1', 'rate of change in NPP with  increment in leaf mass'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_transpiration_talloc_mavg', mem%unit_transpiration_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('unit_transpiration_talloc_mavg', '', 'rate of tranpiration change per unit LAI change?'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_uptake_n_talloc_mavg', mem%unit_uptake_n_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('unit_uptake_n_talloc_mavg', 'micro-mol N mol-1', 'actual rate of N uptkae per unit root'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_uptake_p_talloc_mavg', mem%unit_uptake_p_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('unit_uptake_p_talloc_mavg', 'micro-mol P mol-1', 'actual rate of P uptake per unit root'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dphi_talloc_mavg', mem%dphi_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('dphi_talloc_mavg', 'MPa', 'leaf-soil water potential gradient'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_cn_talloc_mavg', mem%growth_cn_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('growth_cn_talloc_mavg', '', 'plant growth cn ratio, equivalent to 1/growth_req_n'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_cp_talloc_mavg', mem%growth_cp_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('growth_cp_talloc_mavg', '', 'plant growth cp ratio, equivalent to 1 / (growthreq_n * growth_req_p)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_np_talloc_mavg', mem%growth_np_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('growth_np_talloc_mavg', '', 'plant growth hp ratio, equivalent to 1/(growth_req_p)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('labile_carbon_talloc_mavg', mem%labile_carbon_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('labile_carbon_talloc_mavg', '', 'moving average of labile carbon pool, allometry calculations'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('labile_nitrogen_talloc_mavg', mem%labile_nitrogen_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('labile_nitrogen_talloc_mavg', '', 'moving average of labile nitrogen pool, allometry calculations'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('labile_phosphorus_talloc_mavg', mem%labile_phosphorus_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('labile_phosphorus_talloc_mavg', '', 'moving average of labile phosphorus pool, allometry calculations'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('growth_p_limit_based_on_n_mavg', mem%growth_p_limit_based_on_n_mavg, &
        & hgrid, surface, &
        & t_cf('growth_p_limit_based_on_n_mavg', '', 'moving average indicator of how much plant is P limited based on N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('reserve_carbon_talloc_mavg', mem%reserve_carbon_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('reserve_carbon_talloc_mavg', '', 'moving average of reserve carbon, optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('reserve_nitrogen_talloc_mavg', mem%reserve_nitrogen_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('reserve_nitrogen_talloc_mavg', '', 'moving average of reserve nitrogen, optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('reserve_phosphorus_talloc_mavg', mem%reserve_phosphorus_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('reserve_phosphorus_talloc_mavg', '', 'moving average of reserve phosphorus, optimal allocation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_root_lim_talloc_mavg', mem%w_root_lim_talloc_mavg, &
        & hgrid, surface, &
        & t_cf('w_root_lim_talloc_mavg', '', 'moving average of ...'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaf2root_troot_mavg', mem%leaf2root_troot_mavg, &
        & hgrid, surface, &
        & t_cf('leaf2root_troot_mavg', '', 'leaf to root mass ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_npp_troot_mavg', mem%unit_npp_troot_mavg, &
        & hgrid, surface, &
        & t_cf('unit_npp_troot_mavg', 'micro-mol C mol-1', 'rate of change in NPP with increment in leaf mass'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_uptake_n_pot_troot_mavg', mem%unit_uptake_n_pot_troot_mavg, &
        & hgrid, surface, &
        & t_cf('unit_uptake_n_pot_troot_mavg', 'micro-mol N mol-1', 'rate of N uptkae per unit root'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('unit_uptake_p_pot_troot_mavg', mem%unit_uptake_p_pot_troot_mavg, &
        & hgrid, surface, &
        & t_cf('unit_uptake_p_pot_troot_mavg', 'micro-mol P mol-1', 'rate of P uptake per unit root'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fmaint_rate_troot_mavg', mem%fmaint_rate_troot_mavg, &
        & hgrid, surface, &
        & t_cf('fmaint_rate_troot_mavg', '~micro-mol C mol-1 N s-1', 'moving average of N-specific maintentance respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('an_boc_tvegdyn_mavg', mem%an_boc_tvegdyn_mavg, &
        & hgrid, surface, &
        & t_cf('an_boc_tvegdyn_mavg', '', 'net photosynthesis at bottom of canopy, lowest canopy layer existing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('net_growth_tvegdyn_mavg', mem%net_growth_tvegdyn_mavg, &
        & hgrid, surface, &
        & t_cf('net_growth_tvegdyn_mavg', 'micro-mol C m-2 s-1', 'net growth, npp minus turnover'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lai_tvegdyn_mavg', mem%lai_tvegdyn_mavg, &
        & hgrid, surface, &
        & t_cf('lai_tvegdyn_mavg', 'm2 m-2', 'lai at vegdynamics timescale'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('beta_sinklim_ps_tacclim_mavg', mem%beta_sinklim_ps_tacclim_mavg, &
        & hgrid, surface, &
        & t_cf('beta_sinklim_ps_tacclim_mavg', '', 'scaling factor to adjust beta sinklimitation to labile carbon response'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_air_tacclim_mavg', mem%t_air_tacclim_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_tacclim_mavg', '', 'average air temperature for respiration acclimation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_soil_root', mem%t_soil_root, &
        & hgrid, surface, &
        & t_cf('t_soil_root', 'K', 'Temperature of fine roots in the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 273.15_wp)

      CALL mem%Add_var('t_soil_root_tacclim_mavg', mem%t_soil_root_tacclim_mavg, &
        & hgrid, surface, &
        & t_cf('t_soil_root_tacclim_mavg', '', 'average root temperature for respiration acclimation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_air_week_mavg', mem%t_air_week_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_week_mavg', '', 'avrg air temperature across 7 days'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_air_month_mavg', mem%t_air_month_mavg, &
        & hgrid, surface, &
        & t_cf('t_air_month_mavg', '', 'avrg air temperature across 30 days'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_jmax_opt_mavg', mem%t_jmax_opt_mavg, &
        & hgrid, surface, &
        & t_cf('t_jmax_opt_mavg', '', 'optimim temperature for electron transport of photosynthesis'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('net_biosphere_production', mem%net_biosphere_production, &
        & hgrid, surface, &
        & t_cf('net_biosphere_production', '[micro-mol CO2 m-2 s-1]', 'Balance between carbon losses and gains'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('biological_n_fixation', mem%biological_n_fixation, &
        & hgrid, surface, &
        & t_cf('biological_n_fixation', '[micro-mol N m-2 s-1]', 'symbiotic + asymbiotic N fixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

    ENDIF  ! IF One_of(LAND_TYPE  OR VEG_TYPE, lct_ids)
  END SUBROUTINE Init_veg_memory

  ! ======================================================================================================= !
  !>
  !> return name of element variable of bgc_material by its ENUM ID
  !> helper function for t_veg_bgcm_with_elements
  !> "bgcm_object%Get_element_name_by_id" => veg_bgcm_with_elements_get_element_name_by_id
  !>
  FUNCTION veg_bgcm_with_elements_get_element_name_by_id(this, id) RESULT(name)
    CLASS(t_veg_bgcm_with_elements), INTENT(in)  :: this
    INTEGER,                         INTENT(in)  :: id
    CHARACTER(LEN=:),                ALLOCATABLE :: name

    ! return name of element variable only if this bgc_material object contains elements
    IF (this%contains_elements) THEN
      name = get_name_for_element_id(id)
    ELSE
      name='no_elements'
    END IF
  END FUNCTION veg_bgcm_with_elements_get_element_name_by_id

  ! ======================================================================================================= !
  !>
  !> return name of element variable of bgc_material by its ENUM ID
  !> helper function for t_veg_bgcm_with_components
  !> "bgcm_object%Get_element_name_by_id" => veg_bgcm_with_components_get_element_name_by_id
  !>
  FUNCTION veg_bgcm_with_components_get_element_name_by_id(this, id) RESULT(name)
    CLASS(t_veg_bgcm_with_components), INTENT(in)  :: this
    INTEGER,                           INTENT(in)  :: id
    CHARACTER(LEN=:),                  ALLOCATABLE :: name

    ! return name of element variable only if this bgc_material object contains elements
    IF (this%contains_elements) THEN
      name = get_name_for_element_id(id)
    ELSE
      name='no_elements'
    END IF
  END FUNCTION veg_bgcm_with_components_get_element_name_by_id

  ! ======================================================================================================= !
  !>
  !> return name of element variable of bgc_material by its ENUM ID
  !> helper function for t_veg_bgcm_with_components
  !> "bgcm_object%Get_element_name_by_id" => veg_bgcm_main_get_element_name_by_id
  !>
  FUNCTION veg_bgcm_main_get_element_name_by_id(this, id) RESULT(name)
    CLASS(t_veg_bgcm_main), INTENT(in)  :: this
    INTEGER,                INTENT(in)  :: id
    CHARACTER(LEN=:),       ALLOCATABLE :: name

    ! return name of element variable only if this bgc_material object contains elements
    IF (this%contains_elements) THEN
      name = get_name_for_element_id(id)
    ELSE
      name='no_elements'
    END IF
  END FUNCTION veg_bgcm_main_get_element_name_by_id

  ! ======================================================================================================= !
  !>
  !> initializing the t_veg_bgcm_with_elements object
  !> "bgcm_object%Init()"
  !> by default it does contain element variables !
  !>
  SUBROUTINE veg_bgcm_with_elements_init(this, name, shortname, elements_index_map, &
    &                                    l_elements, element_list, element_unit)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_veg_bgcm_with_elements),    INTENT(inout) :: this                   ! object of TYPE t_sb_elements_pool
    CHARACTER(LEN=*),                   INTENT(in)    :: name                   ! pool name       (No spaces !)
    CHARACTER(LEN=*),                   INTENT(in)    :: shortname              ! pool short name (No spaces !)
    INTEGER,                            INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
    LOGICAL,          OPTIONAL,         INTENT(in)    :: l_elements             ! IF pool has elements
    INTEGER,          OPTIONAL,         INTENT(in)    :: element_list(:) ! @TODO remove | list of elements; used if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL,         INTENT(in)    :: element_unit           ! unit of elements
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: element_id
    CHARACTER(len=*), PARAMETER :: routine = modname//':veg_bgcm_with_elements_init'
    ! ----------------------------------------------------------------------------------------------------- !

    this%name              = name
    this%shortname         = shortname
    this%contains_elements = l_elements
    this%contains_bgcm     = .FALSE.

    IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')
    this%element_unit = element_unit

    ! add element variables
    DO element_id = 1, SIZE(elements_index_map)
      IF (elements_index_map(element_id) > 0) THEN
        CALL this%Add_element(get_name_for_element_id(element_id), get_short_name_for_element_id(element_id), element_id, dim=2)

        SELECT TYPE (selector => this%element_list(elements_index_map(element_id))%p)
        CLASS IS (t_jsb_var_real2d)
          SELECT CASE(element_id)
          CASE(ELEM_C_ID)
            this%carbon => selector
          CASE(ELEM_N_ID)
            this%nitrogen => selector
          CASE(ELEM_P_ID)
            this%phosphorus => selector
          CASE(ELEM_C13_ID)
            this%carbon13 => selector
          CASE(ELEM_C14_ID)
            this%carbon14 => selector
          CASE(ELEM_N15_ID)
            this%nitrogen15 => selector
          CASE DEFAULT
            CALL finish(routine, 'Invalid element id'//int2string(element_id))
          END SELECT
        END SELECT
      ENDIF
    ENDDO

  END SUBROUTINE veg_bgcm_with_elements_init

  ! ======================================================================================================= !
  !>
  !> initializing the t_veg_bgcm_with_components object
  !> "bgcm_object%Init()"
  !> by default it does not contain element variables (only t_veg_bgcm_with_elements do)
  !>
  SUBROUTINE veg_bgcm_with_components_init(this, name, shortname, elements_index_map, &
    &                                      l_elements, element_list, element_unit)
    CLASS(t_veg_bgcm_with_components),    INTENT(inout) :: this                   ! object of TYPE t_sb_components_pool
    CHARACTER(LEN=*),                     INTENT(in)    :: name                   ! pool name       (No spaces !)
    CHARACTER(LEN=*),                     INTENT(in)    :: shortname              ! pool short name (No spaces !)
    INTEGER,                              INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
    LOGICAL,          OPTIONAL,           INTENT(in)    :: l_elements             ! IF pool has elements
    INTEGER,          OPTIONAL,           INTENT(in)    :: element_list(:)        ! list of elements; used if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL,           INTENT(in)    :: element_unit           ! unit of elements
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':veg_bgcm_with_components_init'

    this%name              = name
    this%shortname         = shortname
    this%contains_elements = l_elements
    this%contains_bgcm     = .TRUE.

    IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')
    this%element_unit = element_unit
  END SUBROUTINE veg_bgcm_with_components_init

  ! ======================================================================================================= !
  !>
  !> initializing the t_veg_bgcm_main object
  !> "bgcm_object%Init()"
  !> by default it does not contain element variables (only t_veg_bgcm_with_elements do)
  !>
  SUBROUTINE veg_bgcm_main_init(this, name, shortname, elements_index_map, &
    &                           l_elements, element_list, element_unit)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_veg_bgcm_main),         INTENT(inout) :: this                   ! object of TYPE t_sb_components_pool
    CHARACTER(LEN=*),               INTENT(in)    :: name                   ! pool name       (No spaces !)
    CHARACTER(LEN=*),               INTENT(in)    :: shortname              ! pool short name (No spaces !)
    INTEGER,                        INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
    LOGICAL,          OPTIONAL,     INTENT(in)    :: l_elements             ! IF pool has elements
    INTEGER,          OPTIONAL,     INTENT(in)    :: element_list(:)        ! list of elements; used if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL,     INTENT(in)    :: element_unit           ! unit of elements
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':veg_bgcm_main_init'
    ! ----------------------------------------------------------------------------------------------------- !
    this%name              = name
    this%shortname         = shortname
    this%contains_elements = l_elements
    this%contains_bgcm     = .TRUE.
  END SUBROUTINE veg_bgcm_main_init

  ! ======================================================================================================= !
  !>
  !> create, init and append a veg bgcm with components
  !>
  FUNCTION Add_veg_bgc_material_with_components(pool, bgcm_id, name, shortname, unit, elements_index_map) &
      & RESULT(veg_bgcm_with_components)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_pool), TARGET,  INTENT(INOUT)  :: pool                  !< pool to which this pool with component should be added
    INTEGER,                    INTENT(in)     :: bgcm_id               !< id of this veg bgcm with components
    CHARACTER(LEN=*),           INTENT(in)     :: name                  !< name of this veg bgcm with components
    CHARACTER(LEN=*),           INTENT(in)     :: shortname             !< shortname of this veg bgcm with components
    CHARACTER(LEN=*),           INTENT(in)     :: unit                  !< unit of this veg bgcm with components
    INTEGER,                    INTENT(in)     :: elements_index_map(:) !< elements ID -> IND
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_veg_bgcm_with_components),  POINTER :: veg_bgcm_with_components !< the new veg bgcm with components
    ! ----------------------------------------------------------------------------------------------------- !
    ! bgcm_components init
    ALLOCATE(veg_bgcm_with_components)
    CALL veg_bgcm_with_components%Init(name, shortname, elements_index_map(:), l_elements=.FALSE., element_unit=unit)

    ! bgcm_elements allocate & init & add to bgcm_components
    !   "shortname" is used for model output
    !   "element_list" is optional and can be used to select specific elements (not needed when all elements may be added)

    ! leaf
    veg_bgcm_with_components%leaf => Add_veg_bgc_material_with_elements(veg_bgcm_with_components,                     &
      & name = 'leaf_bgcm', shortname = 'leaf', elements_index_map = elements_index_map(:), unit = unit)
    ! fine_root
    veg_bgcm_with_components%fine_root => Add_veg_bgc_material_with_elements(veg_bgcm_with_components,                &
      & name = 'fine_root_bgcm', shortname = 'fine_root', elements_index_map = elements_index_map(:), unit = unit)
    ! coarse_root
    veg_bgcm_with_components%coarse_root => Add_veg_bgc_material_with_elements(veg_bgcm_with_components,              &
      & name = 'coarse_root_bgcm', shortname = 'coarse_root', elements_index_map = elements_index_map(:), unit = unit)
    ! sap_wood
    veg_bgcm_with_components%sap_wood => Add_veg_bgc_material_with_elements(veg_bgcm_with_components,                 &
      & name = 'sap_wood_bgcm', shortname = 'sap_wood', elements_index_map = elements_index_map(:), unit = unit)
    ! heart_wood
    veg_bgcm_with_components%heart_wood => Add_veg_bgc_material_with_elements(veg_bgcm_with_components,               &
      & name = 'heart_wood_bgcm', shortname = 'heart_wood', elements_index_map = elements_index_map(:), unit = unit)
    ! labile
    veg_bgcm_with_components%labile => Add_veg_bgc_material_with_elements(veg_bgcm_with_components,                   &
      & name = 'labile_bgcm', shortname = 'labile', elements_index_map = elements_index_map(:), unit = unit)
    ! reserve
    veg_bgcm_with_components%reserve => Add_veg_bgc_material_with_elements(veg_bgcm_with_components,                  &
      & name = 'reserve_bgcm', shortname = 'reserve', elements_index_map = elements_index_map(:), unit = unit)
    ! fruit
    veg_bgcm_with_components%fruit => Add_veg_bgc_material_with_elements(veg_bgcm_with_components,                    &
      & name = 'fruit_bgcm', shortname = 'fruit', elements_index_map = elements_index_map(:), unit = unit)
    ! seed_bed
    veg_bgcm_with_components%seed_bed => Add_veg_bgc_material_with_elements(veg_bgcm_with_components,                 &
      & name = 'seed_bed_bgcm', shortname = 'seed_bed', elements_index_map = elements_index_map(:), unit = unit)

    ! Add the bookkeepoing for this veg bgcm with components
    CALL pool%Add_pool(veg_bgcm_with_components, bgcm_id)

  END FUNCTION Add_veg_bgc_material_with_components

  ! ======================================================================================================= !
  !>
  !> create, init and append a veg bgcm with elements
  !>
  FUNCTION Add_veg_bgc_material_with_elements(pool,  name, shortname, elements_index_map, unit, bgcm_id) &
      & RESULT(veg_bgcm_with_elements)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_pool), TARGET,  INTENT(INOUT) :: pool                   !< pool to which this pool with elements should be added
    CHARACTER(LEN=*),           INTENT(in)    :: name                   !< name of this veg bgcm with elements
    CHARACTER(LEN=*),           INTENT(in)    :: shortname              !< shortname of this veg bgcm with elements
    INTEGER,                    INTENT(in)    :: elements_index_map(:)  !< mapping elements ID -> IND
    CHARACTER(LEN=*),           INTENT(in)    :: unit                   !< unit of this veg bgcm with elements
    INTEGER,          OPTIONAL, INTENT(in)    :: bgcm_id                !< optional: id of this bgcm with elements
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_veg_bgcm_with_elements), POINTER   :: veg_bgcm_with_elements !< the new veg bgcm with elements
    ! ----------------------------------------------------------------------------------------------------- !

    ALLOCATE(veg_bgcm_with_elements)
    CALL veg_bgcm_with_elements%Init(name = name, shortname = shortname,  &
      & elements_index_map = elements_index_map(:), l_elements = .TRUE., element_unit = unit)
    CALL pool%Add_pool(veg_bgcm_with_elements, bgcm_id)

  END FUNCTION Add_veg_bgc_material_with_elements

#endif
END MODULE mo_veg_memory_class
