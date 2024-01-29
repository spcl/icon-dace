!> vegetation constants (QUINCY)
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
!>#### declare and define vegetation constants
!>
MODULE mo_veg_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp

  IMPLICIT NONE
  PUBLIC


  ! Respiration parameters
  REAL(wp), SAVE :: & 
      ! growth respiration per unit new biomass
      fresp_growth      , &   !< growth respiration per unit new biomass 
      ! N-specific foliar maintenance respiration rate (µmol CO2 / mol N / s)
      ! Taken from Sprugel et al. Respiration from the organ level to the stand. 
      ! In: Smith WK, Hinckley, T.M., editor. Physiological Ecology of Coniferous Forests.
      ! San Diego, CA, USA: Academic Press; 1996.
      ! Note converted from 7.4 e-7 to 1.25 e-6 to account for the different values of the Lloyd&Taylor equation
      ! such that the applied rates a comparable at 15-20°C
      fmaint_rate_base  , &   !< N-specific maintenance respiration rate for roots and leaves (µmol CO2 / mmol N / s) 
      fmaint_rate_wood  , &   !< N-specific maintenance respiration rate for wood (µmol CO2 / mmol N / s) 
      ! parameters in the Lloyd and Taylor temperature response function (Lloyd & Taylor, 1996, Funct. Ecol.) 
      t_maint_k1        , &   !< Coefficient in the equation for temperature sensitivity of foliar respiration
      t_maint_k2        , &   !< Coefficient in the equation for temperature sensitivity of foliar respiration
      t_maint_k3        , &   !< Coefficient in the equation for temperature sensitivity of foliar respiration
      f_resp_acclim     , &   !< respiration temperature acclimation factor 
      t_acclim_zero           !< base temperature for respiration acclimation 

  ! Parameters describing the PS ~ An relationship (describing photosynthetic activity )
  REAL(wp), SAVE :: &      
      k0_fn_pepc_C3     , &  !< fraction of N in PEP and PPDK (C3 plants) 
      k0_fn_pepc_C4     , &  !< fraction of N in PEP and PPDK (C4 plants) (Makino et al. 2003 / Tazoe et al. PCE 2006, 29, 691-700)
      ! in static leaf N mode: prescribed ratio of Jmax25/Vcmax25
      jmax2vcmax_C3     , &  !< text book; ratio of Jmax25/Vcmax25 for prescribed leaf N (C3 plants)
      jmax2vcmax_C4          !< calibrated to have Jmax=vcmax at 25 deg; ratio of Jmax25/Vcmax25 for prescribed leaf N (C4 plants)

  ! Labile pool parameters 
  REAL(wp), SAVE :: &
      ! default turnover time of the labile carbo-hydrates (in days)
      tau_labile,           & !< default turnover time of the labile carbo-hydrates (in days)
      ! temperature scalar on turnover time of labile pool (1-exp(-(lambda * T)^k))
      ! tuned to reduce growth between 5°C to zero at 0°C
      lambda_labile_temp,   & !< lambda in temperature response function of labile pool 
      k_labile_temp,        & !< k in temperature response function of labile pool 
      ! moisture scalar on turnover time of labile pool (1-exp(-(lambda * theta)^k))
      ! tuned to reduce growth from 20% saturation to zero at 0% saturation
      lambda_labile_theta,  & !< lambda in moisture response function of labile pool 
      k_labile_theta,       & !< k in moisture response function of labile pool 
      ! pull from reserve to labile pool due to phenology (1-exp(-(lambda * lai/target_lai)^k))
      lambda_phiphen,       & !< lambda in storage response function to phenology 
      k_phiphen,            & !< k in storage response function to phenology 
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      lambda_phimaint,      & !< lambda in storage response function to maintenance needs 
      k_phimaint,           & !< k in storage response function to maintenance needs
      lambda_phimaint_nut,  & !< lambda in storage response function to maintenance needs 
      k_phimaint_nut,       & !< k in storage response function to maintenance needs
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      lambda_phistore,      & !< lambda in labile response function to storage size
      k_phistore,           & !< k in labile response function to storage size 
      ! sink limitation signal strength to photosynthesis (exp(-(lambda * labile/target_labile)^k))
      knut_labile,          & !< rate at which nutrients can retrieved quicker for growth than carbon
      lambda_sinklim_ps,    & !< lambda in photosynthesis response function to labile pool size 
      k_sinklim_ps,         & !< k in photosynthesis response function to labile pool size 
      k_cnp_sinklim_ps,     & !< scaling factor to limit PS if labile CNP becomes too wide 
      beta_sinklim_ps_min,  & !< maximum downacclimation of photosynthetic capacity due to sink limitation 
      fstore_leaf_max,      & !< maximum reserve storage in leaves relative to leaf mass
      fstore_fine_root_max, & !< maximum reserve storage in fine roots relative to fine root mass
      fstore_sap_wood_max,  & !< maximum reserve storage in sap wood relative to sap wood mass
      k_phi_interact          !< value of reserve pull phi_X, at which the inverse reserve pull phi_Y is reduced

  ! Resorption parameters
  REAL(wp), SAVE :: & 
      resorp_fract_leaf,    & !< default fraction of nutrient resporption before leaf shedding
      resorp_fract_wood,    & !< default fraction of nutrient resporption before wood death
      root2leaf_cn,         & !< relative C:N of fine roots compared to leaves (mol C / mol N)
      wood2leaf_cn,         & !< relative C:N of sapwood biomass compared to leaves (mol C / mol N)
      root2leaf_np,         & !< relative N:P of fine roots compared to leaves (mol N / mol P)
      wood2leaf_np,         & !< relative N:P of sapwood biomass compared to leaves (mol N / mol P)
      tau_nutrient_recycling  !< time scale of enzyme, RNA etc. turnover (as in return to labile pool!, days) 

  ! Parameters concerning N uptake calculations
  REAL(wp), SAVE :: &
      kappa_f_demand_nc    , & !< fraction of maximum labile nitrogen concentration at which uptake is reduced to 50% 
      kappa_f_demand_pn    , & !< fraction of maximum labile phosphorus  concentration at which uptake is reduced to 50% 
      k_f_demand           , & !< k in nutrient uptake response funktion to labile nutrient concentration
      cost_n_fixation_max  , & !< maximum cost of N fixation at low / high temperatures (mol C / mol N)
      cost_n_fixation_min  , & !< mimumum cost of N fixation at 25°C (mol C / mol N)
      k1_temp_bnf          , & !< k1 .. k4 factors in temperature response function of nitrogenase activity
      k2_temp_bnf          , &
      k3_temp_bnf          , &
      k4_temp_bnf          , &
      transform_cost_nh4   , & !< N transformation costs of converting NH4 to organic N  (mol C / mol N)
      transform_cost_no3   , & !< N transformation costs of converting NO3 to organic N (mol C / mol N)
      transform_cost_org   , & !< N transformation costs of converting organic N to organic N (mol C / mol N)
      f_nacq_no3           , & !< fractions of NO3 and NH4 uptake based on transformation costs
      f_nacq_nh4           , &
      vmax_bnf             , & !< maximum rate of N fixation (µmol N / mol C fine root / s )
      k_bnf                , & !< half_saturating cost of the response of N fixation (mol C / mol N) 
      offset_bnf           , & !< factor for BNF response in resistance approach (mol C / mol N)
      eta_nfixation        , & !< fractionation of N15 due to N fixation (per mill), parameter value from Robinson, Tree 2001
      f_myc_exudation_max  , & !< maximum fraction of current growth allocated to mycorrhiza
      f_mycorrhization_min , & !< minimum root mycorrhization 
      f_mycorrhization_max     !< maximum root mycorrhization

  ! parameters to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp), SAVE :: &     
      k0_fn_chl_C3  ,   &  !< parameter to describe the varying fraction of N in chlorophyll through the canopy
      k0_fn_chl_C4  ,   &  !< parameter to describe the varying fraction of N in chlorophyll through the canopy
      k1_fn_chl_C3  ,   &  !< parameter to describe the varying fraction of N in chlorophyll through the canopy
      k1_fn_chl_C4  ,   &  !< parameter to describe the varying fraction of N in chlorophyll through the canopy
      kfn_chl              !< parameter to describe the varying fraction of N in chlorophyll through the canopy

  ! parameters used for phenology calculations
  REAL(wp), SAVE :: & 
      min_lai                , & !< minimum LAI beyond which all leaves are shed
      target_lai_max         , & !< maximum LAI target for reserve use calculations 
      max_leaf_shedding_rate , & !< maximum rate of leaf shedding (1/days) 
      k_leafon_canopy            !< scaling coefficient for leaf on in CANOPY-only model 

  ! vegetation dynamics
  REAL(wp), SAVE :: & 
      background_mort_rate_tree   , & !< background mortality rate for trees in static mortality case (1/year)
      background_mort_rate_grass  , & !< background mortality rate for grass in static mortality case (1/year)
      lambda_mort_light           , & !< parameter in the Weibull function controlling self-thinning mortality 
      k_mort_light                , & !< parameter in the Weibull function controlling self-thinning mortality 
      lambda_est_temp             , & !< parameter in the Weibull function controlling temperature-limited establishment 
      k_est_temp                  , & !< parameter in the Weibull function controlling temperature-limited establishment 
      lambda_est_moist            , & !< parameter in the Weibull function controlling moisture-limited establishment 
      k_est_moist                 , & !< parameter in the Weibull function controlling moisture-limited establishment 
      wood_extraction_rate        , & !< fraction of wood removed from site after harvesting
      min_greff                   , & !< minimum growth efficiency (mol C / m2 (LAI) / yr)
      k1_mort_greff               , & !< asymptotic growth efficiency mortality rate (1/yr)
      k2_mort_greff               , & !< scaling coefficient 
      k3_mort_greff               , & !< scaling coefficient in growth efficiency mortality below minimum growth efficiency
      max_crown_area              , & !< maximum crown area (m2)
      min_diameter                , & !< minimum diameter for crown area calculation (m)
      k_crown_area                , & !< scaling parameter in crownarea to diameter relationship
      k_rp                        , & !< scaling exponent in crownarea to diameter relationship
      k_fpc                       , & !< ... see ticket #475
      fpc_max                         !< maximum foliage projective cover

  ! other
  REAL(wp), SAVE :: & 
      ! carbon content of leaves per dry matter
      carbon_per_dryweight_leaf   , & !< carbon mass per unit dry weight of leaves (gC / g DW) 
      ! minimum permissable height of an individual
      min_height                  , & !< minimum permissable height of an individual
      ! stem mass to leaf mass ratio of grasses
      sm2lm_grass                 , & !< stem mass to leaf mass ratio of grasses
      ! Parameters controlling fruit allocation in: exp(-(lambda * reserve_use+offset)^k) 
      lambda_fruit_alloc          , & !< parameter in the Weibull function controlling fruit growth
      k4_fruit_alloc              , & !< parameter in the Weibull function controlling fruit growth
      k3_fruit_alloc              , & !< parameter in the Weibull function controlling fruit growth; ~µmol C / m2 / s
      ! minimum fraction of allocation going to fruit (flowers etc.)
      k1_fruit_alloc              , & !< minimum fraction of allocation going to fruit (flowers etc.)
      w_root_lim_max              , & !< value of the root water limitation below which root allocation starts responding
      leaf2root_min_ratio         , & !< minimum leaf to root ratio for dynamic allocation
      leaf2root_max_ratio             !< maximum leaf to root ratio for dynamic allocation

  ! parameters for changes in canopy N fractions
  REAL(wp), SAVE :: &     
      ! maximum change in photosynthetic N fractions
      delta_n_fraction    , & !< maximum fractional change for optimal adjustment of photosynthetic N fraction
      ! maximum change in CN ratio of canopy
      delta_n_leaf            !< maximum fractional change for optimal adjsutment of canopy N allocation


  ! some useful numbers for easy identification of photosynthetic pathways and growth forms
  ! need to be parameters (i.e. constant across model runtime)
      INTEGER, PARAMETER ::    &
      itree         = 1  ,     &
      igrass        = 2  


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_veg_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the vegetation process
  !! called by jsbach_setup_models()
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_veg_constants

    USE mo_jsb_math_constants,      ONLY: one_day, one_year
    USE mo_jsb_physical_constants,  ONLY: molar_mass_C, molar_mass_N, molar_mass_P

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut    
    ! 
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_veg_constants'


  ! Respiration parameters
      fresp_growth                = 0.25_wp
      fmaint_rate_base            = 1._wp 
      fmaint_rate_wood            = 0.25_wp * fmaint_rate_base
      t_maint_k1                  = 308.56_wp  
      t_maint_k2                  = 56.02_wp 
      t_maint_k3                  = 227.13_wp 
      f_resp_acclim               = -0.008_wp
      t_acclim_zero               = 283.15_wp 

  ! Parameters describing the PS ~ An relationship (describing photosynthetic activity )
      k0_fn_pepc_C3               = 0.0_wp 
      k0_fn_pepc_C4               = 0.045_wp 
      jmax2vcmax_C3               = 1.92_wp / 4._wp
      jmax2vcmax_C4               = 1.4_wp  

  ! Labile pool parameters 
      tau_labile                  = 5.0_wp  
      lambda_labile_temp          = 0.5_wp 
      k_labile_temp               = 2.0_wp 
      lambda_labile_theta         = 10.0_wp
      k_labile_theta              = 2.0_wp 
      lambda_phiphen              = 1.3_wp 
      k_phiphen                   = 8.0_wp 
      lambda_phimaint             = 4.0_wp 
      k_phimaint                  = 1.2_wp 
      lambda_phimaint_nut         = 1.6_wp 
      k_phimaint_nut              = 3._wp 
      lambda_phistore             = 2._wp
      k_phistore                  = 3._wp 
      knut_labile                 = 1.2_wp 
      lambda_sinklim_ps           = 0.1_wp
      k_sinklim_ps                = 2.0_wp 
      k_cnp_sinklim_ps            = 4.0_wp 
      beta_sinklim_ps_min         = 0.25_wp 
      fstore_leaf_max             = 0.02_wp
      fstore_fine_root_max        = 0.2_wp
      fstore_sap_wood_max         = 0.15_wp
      k_phi_interact              = 0.75_wp


  ! Resorption parameters
      resorp_fract_leaf           = 0.5_wp  
      resorp_fract_wood           = 0.2_wp 
      root2leaf_cn                = 0.85_wp 
      wood2leaf_cn                = 0.145_wp 
      root2leaf_np                = 1.0_wp 
      wood2leaf_np                = 1.0_wp 
      tau_nutrient_recycling      = 10._wp

  ! N uptake parameters
      kappa_f_demand_nc           = 0.75_wp 
      kappa_f_demand_pn           = 0.9_wp 
      k_f_demand                  = 2.0_wp
      cost_n_fixation_max         = 17.5_wp   / (molar_mass_C/molar_mass_N) 
      cost_n_fixation_min         = 8.0_wp    / (molar_mass_C/molar_mass_N) 
      k1_temp_bnf                 = 1.25_wp
      k2_temp_bnf                 = -3.62_wp
      k3_temp_bnf                 = 0.27_wp
      k4_temp_bnf                 = 50.3_wp
      transform_cost_nh4          = 1.7_wp    / (molar_mass_C/molar_mass_N)
      transform_cost_no3          = 2.3_wp    / (molar_mass_C/molar_mass_N)
      transform_cost_org          = 1.0_wp    / (molar_mass_C/molar_mass_N)
      f_nacq_nh4                  = transform_cost_no3 / (transform_cost_no3 + transform_cost_nh4)
      f_nacq_no3                  = transform_cost_nh4 / (transform_cost_no3 + transform_cost_nh4)
      vmax_bnf                    = 5.0_wp * 0.0225_wp / (molar_mass_N/molar_mass_C) / (one_day*one_year) * 1.e6_wp
      k_bnf                       = 50.0_wp   / (molar_mass_C/molar_Mass_N)
      offset_bnf                  = 6.0_wp * (f_nacq_nh4 * transform_cost_nh4 + f_nacq_no3 * transform_cost_no3)
      eta_nfixation               = 3.0_wp                   ! Robinson, Tree 2001
      f_myc_exudation_max         = 0.3_wp
      f_mycorrhization_min        = 0.01_wp
      f_mycorrhization_max        = 0.5_wp

  ! parameters to describe the varying fraction of N in chlorophyll through the canopy
      k0_fn_chl_C3                = 6.0_wp
      k0_fn_chl_C4                = 15.0_wp
      k1_fn_chl_C3                = 3.6_wp
      k1_fn_chl_C4                = 4.4_wp
      kfn_chl                     = 0.7_wp

  ! parameters used for phenology calculations
      min_lai                     = 0.1_wp
      target_lai_max              = 5.0_wp 
      max_leaf_shedding_rate      = 0.05_wp
      k_leafon_canopy             = 0.2_wp

  ! vegetation dynamics
      background_mort_rate_tree   = 0.01_wp 
      background_mort_rate_grass  = 0.05_wp 
      lambda_mort_light           = 0.15_wp 
      k_mort_light                = 4.0_wp  
      lambda_est_temp             = 0.075_wp  
      k_est_temp                  = 4.0_wp  
      lambda_est_moist            = 10._wp  
      k_est_moist                 = 2.0_wp  
      wood_extraction_rate        = 0.8_wp
      k1_mort_greff               = 0.05_wp 
      k2_mort_greff               = 0.3_wp 
      k3_mort_greff               = 25.0_wp 
      min_greff                   = 1.0_wp
      max_crown_area              = 15._wp
      min_diameter                = 0.01_wp
      k_crown_area                = 100._wp
      k_rp                        = 1.6_wp
      k_fpc                       = 0.5_wp
      fpc_max                     = 0.95_wp

  ! other
      carbon_per_dryweight_leaf   = 0.48_wp  
      min_height                  = 0.1_wp 
      sm2lm_grass                 = 0.05_wp 
      lambda_fruit_alloc          = 10.0_wp 
      k1_fruit_alloc              = 0.01_wp 
      k3_fruit_alloc              = 0.1_wp  
      k4_fruit_alloc              = 2.0_wp
      w_root_lim_max              = 0.8_wp
      leaf2root_min_ratio         = 0.2_wp
      leaf2root_max_ratio         = 15.0_wp

  ! parameters for changes in canopy N fractions
!@todo:should this be 1/20days?
      delta_n_fraction            = 0.048_wp        ! 0.001_wp (changed to account for multiplication with dtime/one_day in the code)
      delta_n_leaf                = 0.0048_wp       ! 0.0001 (changed to account for multiplication with dtime/one_day in the code)   

  END SUBROUTINE init_veg_constants

#endif
END MODULE mo_veg_constants
