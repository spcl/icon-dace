!> QUINCY soil-biogeochemistry constants
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
!>#### declare and define soil-biogeochemistry constants
!>
MODULE mo_sb_constants
#ifndef __NO_QUINCY__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! turnover times (mostly for simple_1d model)
  REAL(wp), SAVE :: &
      tau_soluable_litter, &    !< turnover time of soluable litter [s]
      tau_metabolic_litter, &   !< turnover time of metabolic litter pool [s]
      tau_structural_litter, &  !< turnover time of structural litter pool [s]
      tau_woody_litter, &       !< turnover time of woody litter [s]
      tau_fast, &               !< turnover time of fast SOM pool [s]
      tau_slow, &               !< turnover time of slow SOM pool [s]
      tau_biomineralisation     !< turnover time of slow SOM pool for biomineralisation [s]

  ! temperature and moisture rate modifiers
  REAL(wp), SAVE :: &
      t_ref_tresponse, &      !< reference temperature for temperature response function [K]
      ea_depolymerisation, &  !< activation energy for depolymerisation [J mol-1]
      ea_mic_uptake, &        !< activation energy for microbial uptake [J mol-1]
      ea_sorption, &          !< activation energy for sorption to mineral surfaces [J mol-1]
      ea_desorption, &        !< activation energy for desorption to mineral surfaces [J mol-1]
      ea_hsc, &               !< activation energy for half-saturation constants [J mol-1]
      t_ref_decomposition, &  !< temperature of peak decomposition rate [K]
      ea_decomposition, &     !< activation energy for decomposition [J mol-1]
      ed_decomposition, &     !< de-activation energy for decomposition [J mol-1]
      k1_afps, &              !< scaling factor in the air-filled pore space to modifier response
      k2_afps, &              !< scaling factor in the substrate modifier response
      kmr_depolymerisation, & !< scaling factor for the sensitivity of depolymerisation to oxygen limitation
      kmr_mic_uptake, &       !< scaling factor for the sensitivity of microbial uptake to oxygen limitation
      kmr_hsc, &              !< scaling factor for the sensitivity of half-saturation constant to moisture limitation
      kmr_sorption, &         !< scaling factor for the sensitivity of surface adsorption to moisture limitation
      kmr_decomposition, &    !< minimum water potential for decomposition
      kmr_asymb_bnf           !< minimum water potential for nfixation

  ! transport and sorption parameters
  REAL(wp), SAVE :: &
      k_diff_org, &                 !< diffusion coefficient for organic material due to bioturbation? [m2 (kg m-3) s-1]
      rho_bulk_org, &               !< bulk density of organic material [kg m-3]
      carbon_per_dryweight_SOM, &   !< C per unit dryweight SOM [gC / gOM]
      k_adsorpt_dom, &              !<
      k_desorpt_dom, &              !<
      k_adsorpt_det, &              !<
      k_desorpt_det, &              !<
      k_adsorpt_nh4, &              !< slow adsorption rate (from soluble Pi to sorbed Pi) of PO4, [kg soil mol-1 P s-1]
      k_desorpt_nh4, &              !< fast desorption rate (from NH4assoc to NH4solute) of nh4, [s-1]
      k_desorpt_po4, &              !< Slow desorption rate (between sorbed Pi and labile Pi) of PO4, [s-1], simple_1d model
      k_adsorpt_po4, &              !< Slow adsorption rate (between sorbed Pi and labile Pi) of PO4, [s-1], simple_1d model
      k_adsorpt_f_po4, &            !< fast adsorption rate (from soluble Pi to labile Pi) of PO4, [kg soil mol-1 P s-1]
      k_desorpt_f_po4               !< fast desorption rate (from labile Pi to soluble Pi) of PO4, [s-1]

  ! microbial parameters
  REAL(wp), SAVE :: &
      microbial_cn, &            !< molar ratio of C:N in microbes (testing only!)
      microbial_np, &            !< molar ratio of N:P in microbes (testing only!)
      microbial_nue, &           !< microbial nitrogen use efficiency
      microbial_pue, &           !< microbial phosphorus use efficiency
      microbial_cue, &           !< microbial carbon-use efficiency
      microbial_cue_max, &       !< maximum microbial carbon-use efficiency
      microbial_cue_min          !< minimum microbial carbon-use efficiency

  ! fast and slow parameters for simple_1d SB model
  REAL(wp), SAVE :: &
      k_fast_som_cn_max, &       !< maximum microbial C:N [mol/mol]
      k_fast_som_cn_min, &       !< minimum microbial C:N [mol/mol]
      f_fast_som_cn, &           !< slope of microbial C:N to soil NH4 [kg/mol]
      k_slow_som_cn, &           !< molar ratio of C:N in slow SOM
      k_fast_som_np, &           !< molar ratio of N:P in fast SOM at reference NH4 [mol/m3]
      kref_fast_som_np, &        !< reference NH4 for fast SOM N:P
      k_slow_som_np, &           !< molar ratio of N:P in slow SOM
      frac_woody2dec_litter, &   !< fraction of woody litter transformed into soluable or polymeric litter
      frac_litter2fast_som, &    !< fraction of litter transformed into fast SOM
      frac_fast2slow_som, &      !< fraction of fast SOM transformed into slow SOM
      frac_slow2fast_som         !< fraction of slow SOM transformed into fast SOM

  ! other parameters for initialisation
  REAL(wp), SAVE :: &
      k_init_soluable_litter_cn, &    !< initial soluable litter C:N [mol/mol]
      k_init_polymeric_litter_cn, &   !< initial polymeric C:N [mol/mol]
      k_init_woody_litter_cn, &       !< initial woody litter C:N [mol/mol]
      k_init_litter_np, &             !< initial litter N:P [mol/mol]
      frac_prod_not_root, &           !< fraction of production not roots
      min_fraction_default_sb         !< minimum fraction for various cases in the SB_ process [--]

  ! additional microbial parameters for JSM implementation
  REAL(wp), SAVE :: &
      fc_woody2polymeric, &             !< fraction of woody litter transformed into polymeric litter (JSM version)
      fc_soluable2dom, &                !< fraction of soluable litter transformed into dom pool (JSM version)
      k_mort_microbial, &               !< mortality rate of microbes [s-1]
      min_mic_biomass, &                !< minimum microbial biomass density (stopping mortality) [mol C/m3]
      min_enzyme_fraction, &            !< minimum enzymatic fraction [--]
      max_enzyme_fraction, &            !< maximum enzymatic fraction [--]
      init_enzyme_fraction, &           !< initial values for enzyme allocation etc.[--]
      vmax_mic_uptake, &                !< maximum uptake rates of microbes [s-1]
      km_mic_uptake, &                  !< half-saturation DOM density for microbial DOM consumption
      km_depolymerisation, &            !< half-saturation microbial biomass density for depolymerisation
      vmax_depolymerisation_residue, &  !< maximum depolymerisation rate of residue
      vmax_depolymerisation_litter, &   !< maximum depolymerisation rate of polymeric litter
      enzyme_frac_total, &              !< enzyme fraction of the total microbial biomass
      fc_p_res2dom, &                   !< fraction of residue P retranslocated to dom
      fc_n_res2dom, &                   !< fraction of residue n retranslocated to dom
      k_p_recyc_corr, &                 !< correction coefficient for microbial P recycling exponential decline curve
      k_n_recyc_corr                    !< correction coefficient for microbial P recycling exponential decline curve

  ! parameters controlling litter partitioning
  REAL(wp), SAVE :: &
      f_soluable_residue, &      !< fraction of microbial residue entering the DOM pool
      f_soluable_leaf, &         !< fraction of leaf litter that becomes DOM
      f_soluable_fine_root, &    !< fraction of fine root litter that becomes DOM
      f_soluable_fruit, &        !< fraction of fruit litter that becomes DOM
      f_soluable_seed_bed, &     !< fraction of seed bed litter that becomes DOM
      fc_soluable_max, &         !< maximum fraction of soluable litter formation
      k_fc_soluable, &           !< slope of soluable fraction with lignin to N ratio
      k_fn_soluable, &           !< proportionality factor controlling C:N of solubable vs. polymeric pool
      k_fp_soluable, &           !< proportionality factor controlling C:P of solubable vs. polymeric pool
      lc_leaf_max, &             !< maximum lignin content of leaves [1/mol]
      lc_leaf2sla, &             !< slope of lignin to sla relationship [1/m2]
      lc_fine_root, &            !< lignin content of fine roots [1/mol]
      lc_coarse_root, &          !< lignin content of coarse roots [1/mol]
      lc_woody_litter, &         !< lignin content of woody litter [1/mol]
      lc_fruit, &                !< lignin content of fruits [1/mol]
      lc_seed_bed                !< lignin content of seed bed [1/mol]

  ! mycorrhizal parameters
  REAL(wp), SAVE :: &
      mycorrhiza_cue, &          !< mycorrhizal carbon-use efficiency
      mycorrhiza_cn, &           !< molar ratio of C:N of mycorrhizal fungi
      mycorrhiza_np, &           !< molar ratio of N:P of mycorrhizal fungi
      f_maxresp_myc, &           !< maximum fraction of mycorrhizal carbon used for N transformation [1/day]
      f_surface_myc2fineroot, &  !< ratio of surface to mass relationship between mycorrhizal fungi and fine roots
      myc_root_coverage_max, &   !< maximum of root (tips) that can be covered by mycorrhizae
      vmax_uptake_norg, &        !< uptake rate of organic nitrogen by mycorrhiza
      km_up_org, &               !< half-saturation constant for organic uptake [mol N m-3]
      f_nin2norg, &              !< scaling factor to limit organic uptake according to accessibility to substrate
      f_norg_limit, &            !< smoothing factor to limit organic uptake
      vmax_priming, &            !< turnover time reduction of slow pool by mycorrhizae
      km_priming, &              !< half-saturation constant for turnover time reduction of slow pool by mycorrhizae [mol C m-3]
      resp_priming, &            !< respiration associated with priming
      organic_cn, &              !< molar ratio of C:N of mycorrhizal export to plants (amino acids)
      f_tau_myc                  !< increase of tau_slow in case of saprotrophic/priming macorrhizae to get similar tau_slow

  ! nutrient uptake kinetics & ECA approach: simple_1d model
  REAL(wp), SAVE :: &
      km1_up_nh4, &             !< low-affinity parameter for plant and mycorrhizal uptake [m3/mol]
      km1_up_no3, &             !< low-affinity parameter for plant and mycorrhizal uptake [m3/mol]
      km1_up_po4, &             !< low-affinity parameter for plant and mycorrhizal uptake [l/mol]  unit correct ?
      km2_up_nh4, &             !< high-affinity parameter for plant and mycorrhizal uptake [mol/m3]
      km2_up_no3, &             !< high-affinity parameter for plant and mycorrhizal uptake [mol/m3]
      km2_up_po4                !< high-affinity parameter for plant and mycorrhizal uptake [mol/l] unit correct ?

  ! nutrient uptake kinetics & ECA approach: additional jsm parameters
  REAL(wp), SAVE :: &
      vmax_mic_uptake_nh4, &      !< maximum NH4 uptake capacity of microbes [mol N / mol biomass C / s]
      km_mic_uptake_nh4, &        !< half-saturation concentration [mol N / m3]
      vmax_mic_uptake_no3, &      !< maximum NO3 uptake capacity of microbes [mol N / mol biomass C / s]
      km_mic_uptake_no3, &        !< half-saturation concentration [mol N / m3]
      f_mic_carrier_n, &          !< ratio between microbial N transport and biomass, []
      vmax_plant_uptake_nh4, &    !< maximum NH4 uptake capacity of fine roots [mol N / mol biomass C / s]
      km_plant_uptake_nh4, &      !< half-saturation concentration [mol N / m3]
      vmax_plant_uptake_no3, &    !< maximum NO3 uptake capacity of fine roots [mol N / mol biomass C / s]
      km_plant_uptake_no3, &      !< half-saturation concentration [mol N / m3]
      f_plant_carrier_n, &        !< ratio between plant N transport and fine root biomass, []
      vmax_nitrification_nh4, &   !< maximum nitrification capacity of nitrifiers [mol N / mol biomass C / s]
      vmax_denitrification_no3, & !< maximum nitrification capacity of denitrifiers [mol N / mol biomass C / s]
      km_nitrification_nh4, &     !< half-saturation concentration of nitrification [mol N / m3]
      km_denit_no3, &             !< half-saturation concentration of denitrification [mol N / m3]
      k_nit_carrier_n, &          !< ratio between N de/nitrifiers and microial biomass, [-]
      vmax_mic_uptake_po4, &      !< maximum PO4 uptake capacity of microbes [mol PO4 / mol biomass C / s]
      km_mic_uptake_po4, &        !< half-saturation concentration for decomposing microbe PO4 immobilization [mol PO4 m-3]
      vmax_plant_uptake_po4, &    !< maximum PO4 uptake capacity of fine roots [mol PO4 / mol biomass C / s]
      km_plant_uptake_po4, &      !< half-saturation concentration [mol PO4 / m3]
      k_mic_uptake_nh4, &         !< Reaction rate of microbe nh4 carrier enzyme [s-1], ECA specific
      k_plant_uptake_nh4, &       !< Reaction rate of plant nh4 carrier enzyme [s-1], ECA specific
      k_mic_uptake_no3, &         !< Reaction rate of microbe no3 carrier enzyme [s-1], ECA specific
      k_plant_uptake_no3, &       !< Reaction rate of plant no3 carrier enzyme [s-1], ECA specific
      k_mic_uptake_po4, &         !< Reaction rate of microbe PO4 carrier enzyme [s-1], ECA specific
      k_plant_uptake_po4, &       !< Reaction rate of plant PO4 carrier enzyme [s-1], ECA specific
      f_mic_carrier_po4, &        !< ratio between microbial P transport and microbial biomass, []
      f_plant_carrier_po4, &      !< ratio between plant P transport and fine root biomass, []
      f_nit_carrier_nh4, &        !<
      f_denit_carrier_no3, &      !<
      qmax_nh4_clay, &            !< Maximum NH4 sorption capacity of clay [mol P kg-1 clay], ECA and Langmuir
      km_adsorpt_mineral_nh4, &   !< Half-saturation concentration for NH4 adsorption to mineral soil [mol P kg-1 fine mineral soil], ECA and Langmuir
      km_adsorpt_OM_nh4, &        !< Half-saturation concentration for NH4 adsorption to OM, [mol N kg-1 OM], ECA and Langmuir
      qmax_po4_mineral, &         !< Maximum PO4 sorption capacity of mineral soil [mol P kg-1 soil], ECA and Langmuir
      km_adsorpt_mineral_po4, &   !< Half-saturation concentration for PO4 adsorption to mineral soil [mol P kg-1 fine mineral soil], ECA and Langmuir
      qmax_po4_clay, &            !< Maximum PO4 sorption capacity of clay [mol P kg-1 fine soil], ECA and Langmuir
      qmax_po4_silt, &            !< Maximum PO4 sorption capacity of silt [mol P kg-1 fine soil], ECA and Langmuir
      qmax_po4_sand, &            !< Maximum PO4 sorption capacity of sand [mol P kg-1 fine soil], ECA and Langmuir
      qmax_po4_OM, &              !< Maximum PO4 sorption capacity of organic matter [mol P kg-1 OM], ECA and Langmuir
      km_po4_ph, &                !< ph corrector for half-saturation concentration for PO4 adsorption [--], ECA and Langmuir
      km_adsorpt_clay_po4, &      !< Half-saturation concentration for PO4 adsorption to clay [mol P kg-1 fine mineral soil], ECA and Langmuir
      km_adsorpt_silt_po4, &      !< Half-saturation concentration for PO4 adsorption to silt [mol P kg-1 fine mineral soil], ECA and Langmuir
      km_adsorpt_sand_po4, &      !< Half-saturation concentration for PO4 adsorption to sand [mol P kg-1 fine mineral soil], ECA and Langmuir
      km_adsorpt_OM_po4, &        !< Half-saturation concentration for PO4 adsorption to OM, [mol P kg-1 OM], ECA and Langmuir
      k_deep_soil_moc_corr, &     !< numerical coefficient for correction of deep soil initialised with very low MOC, assuming the P sorption capacity of OM is one-tenth of mineral soil
      fact_repartition_om_init_cond1, &  !< re-partition factor when OM got over-saturation in SB_ init, soil condition 1 (see docu in sb_init)
      fact_repartition_om_init_cond2     !< re-partition factor when OM got over-saturation in SB_ init, soil condition 2 (see docu in sb_init)

  ! nitrification, denitrification parameters (sb_nloss_scheme = dynamic)
  REAL(wp), SAVE :: &
      ea_nitrification, &        !< activation energy for nitrification [J mol-1]
      ed_nitrification, &        !< deactivation energy for nitrification [J mol-1]
      t_opt_nitrification, &     !< optimum temperature for nitrification [K]
      ea_denitrification, &      !< activation energy for denitrification [J mol-1]
      ea_gasdiffusion, &         !< activation energy for gas diffusion [J mol-1]
      lambda_anvf, &             !< lambda in the Weibull function relating soil moisture to the anaerobic volume fraction
      k_anvf, &                  !< k in the Weibull function relating soil moisture to the anaerobic volume fraction
      vmax_nitrification, &      !< maximum nitrification rate [micro-mol N mol-1 N s-1]
      vmax_denitrification, &    !< maximum denitrification rate [micro-mol m-3 s-1]
      km_denitrification_c, &    !< half-saturation constant for carbon-limited (de-)nitrification [mol C m-3]
      km_denitrification_no3, &  !< half-saturation constant for nitrate-limited denitrifcation [mol N m-3]
      f_nit_noy, &               !< fraction of nitrification lost as NOy [--]
      f_nit_n2o, &               !< fraction of nitrification lost as N2O [--]
      f_denit_noy, &             !< fraction of denitrification lost as NOy [--]
      f_denit_n2o                !< fraction of denitrification lost as N2O [--]

  ! simple_1d model fixed loss and simple_1d leaching parameterisation
  REAL(wp), SAVE :: &
      floss_nmin, &              !< fractional loss of N per unit N mineralisation
      fnit_nmin, &               !< fractional nitrification per unit N mineralisation
      fleach_nh4, &              !< fraction of NH4 amenable to leaching
      fleach_po4, &              !< fraction of PO4 amenable to leaching
      km_asymb_bnf, &            !< half-saturation constant at which BNF w.r.t. soil N [mol N / m2]
      vmax_asymb_bnf             !< maximum fraction of N deficit in N decomposition that is closed by BNF []

  ! Isoptope stuff
  REAL(wp), SAVE :: &
      eta_mic_uptake_nh4, &      !< fractionation of N15 due to microbial NH4 assimilation [per mill]
      eta_mic_uptake_no3, &      !< fractionation of N15 due to microbial NO3 assimilation [per mill]
      eta_plant_uptake_nh4, &    !< fractionation of N15 due to plant NH4 assimilation [per mill]
      eta_plant_uptake_no3, &    !< fractionation of N15 due to plant NO3 assimilation [per mill]
      eta_nitrification, &       !< fractionation of N15 due to microbial nitrification to NO and N2O [per mill]
      eta_nitrate_production, &  !< fractionation of N15 due to microbial NO3 production during nitrification [per mill]
      eta_denitrification, &     !< fractionation of N15 due to microbial denitrification [per mill]
      eta_ammonification, &      !< fractionation of N15 due to microbial ammonification [per mill]
      eta_volatilisation_nh3     !< fractionation of N15 due to NH3 volatilisation [per mill]

  ! soil P fluxes
  REAL(wp), SAVE :: &
      k_weath_mineral, &                !< Weathering coefficient of P in mineral soil (composed of a background rate (frac_background_weath_mineral) and pH/root biomass regulated rate) [s-1]
      frac_background_weath_mineral, &  !< Background (i.e., minimum) rate of weathering coefficient of P in mineral soil 'k_weath_mineral' expressed as fraction of 'k_weath_mineral' [-]
      eq_sorp_mineral_po4, &        !< Slow Pi sorption equilibrium content of mineral soil [mol P kg-1 mineral soil]
      k_occlude_po4, &              !< Occlusion coefficient of sorbed PO4, [s-1]
      vmax_biochem_po4, &           !< Maximum PO4 biochemical mineralization rate [s-1], MM kinetics (potentially also for ECA)
      km_mic_biochem_po4, &         !< Half-saturation microbial biomass for PO4 biochemical mineralization [mol C m-3], MM kinetics (potentially also for ECA)
      km_dip_biochem_po4, &         !< Half-saturation solute P concentration for PO4 biochemical mineralization [mol C m-3], MM kinetics-alike
      km_pc_biochem_po4, &          !< Half-saturation SOM P:C ratio for PO4 biochemical mineralization, MM kinetics-alike
      km_rootc_biochem_po4, &       !< Half-saturation root C biomass for PO4 biochemical mineralization [mol C m-3], MM kinetics-alike
      km_rootc_weath_po4, &         !< Half-saturation root C biomass for PO4 weathering [mol C m-3], simple_1d
      km_root_enzyme_weath_po4, &   !< Half-saturation root enzyme C for PO4 weathering [mol C m-3], jsm
      biomin_threshold_pc_ratio_mic_res, &    !< Threshold P:C ratio for biomineralisation of microbial residue
      biomin_threshold_pc_ratio_min_assoc_om  !< Threshold P:C ratio for biomineralisation of mineral associated OM

  ! soil P pool
  REAL(wp), SAVE :: &
    k1_p_pool_depth_corr, &               !< denominator coefficient for P pool depth correction  [unitless ?]
    k2_p_pool_depth_corr, &               !< numerical coefficient for P pool depth correction  [unitless ?]
    k3_p_pool_depth_corr, &               !< numerical coefficient for P pool depth correction  [unitless ?]
    k4_p_pool_depth_corr, &               !< numerical coefficient for P pool depth correction  [unitless ?]
    k_p_pool_depth_corr_min_occlud, &     !< numerical coefficient of minimum occluded P for P pool depth correction  [unitless ?]
    p_labile_slow_global_avg, &           !< global average value of labile_P/(labile_P + slow_P)  [unitless ?]
    k_sb_pool_dom_assoc_c_corr            !< numerical coefficient for correction for sb_pool%dom_assoc%carbon init

  ! sb_pool default init values
  REAL(wp), SAVE :: &
    def_sb_pool_metabol_litter       , &      ! metabolic litter
    def_sb_pool_simple_struct_litter , &      ! structural litter, simple_1d soil model
    def_sb_pool_jsm_struct_litter    , &      ! structural litter, JSM
    def_sb_pool_woody_litter         , &      ! woody litter
    def_sb_pool_fast_som             , &      ! fast SOM
    def_sb_pool_slow_som             , &      ! slow SOM
    def_sb_pool_microbial_biomass    , &      ! microbial biomass
    def_sb_pool_microbial_necromass  , &      ! microbial necromass
    def_sb_pool_dom                           ! DOM

  ! constants used as factor for multiplication with namelist parameters (spq_config) when "parameter sensitivity" is enabled
  ! see sb_config & mo_qs_set_parameters:Modify_quincy_namelist_parameter
  REAL(wp) :: &
    f_psensi_soil_p_depth, &
    f_psensi_soil_p_labile, &
    f_psensi_soil_p_slow, &
    f_psensi_soil_p_occluded, &
    f_psensi_soil_p_primary, &
    f_psensi_qmax_org_fine_particle


   !! only temporarily, mol / m2 !! @TODO
   REAL(wp), SAVE :: &
      nh4_solute_prescribe, no3_solute_prescribe, po4_solute_prescribe, po4_solute_min

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from soil_biogeochemistry
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_sb_constants

    USE mo_jsb_math_constants,        ONLY: one_day, one_year
    USE mo_jsb_physical_constants,    ONLY: molar_mass_C, molar_mass_N, molar_mass_P

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_sb_constants'


    ! turnover times
    tau_soluable_litter   = 30._wp * one_day                ! 30 days
    tau_metabolic_litter  = 0.033_wp * one_year * one_day   ! Parton et al. 1993, GBC
    tau_structural_litter = 0.124_wp * one_year * one_day   ! Parton et al. 1993, GBC
    tau_woody_litter      = 2.5_wp * one_year * one_day     ! Thum et al. 2019, GMD
    tau_fast              = 2._wp * one_year * one_day      ! Thum et al. 2019, GMD
    tau_slow              = 100._wp * one_year * one_day    ! Thum et al. 2019, GMD
    tau_biomineralisation = 5._wp * one_year * one_day      ! Thum et al. 2019, GMD

    ! temperature and moisture rate modifiers
    t_ref_tresponse       = 293.15_wp                       ! Wang, G., et al. 2012., Soil Biology and Biochemistry 48, 28-38.
    ea_depolymerisation   = 53000._wp                       ! B. Ahrens, pers. comm. 2016
    ea_mic_uptake         = 47000._wp                       ! B. Ahrens, pers. comm. 2016
    ea_sorption           = 5000._wp                        ! B. Ahrens, pers. comm. 2016
    ea_desorption         = 20000._wp                       ! B. Ahrens, pers. comm. 2016
    ea_hsc                = 30000._wp                       ! B. Ahrens, pers. comm. 2016
    ea_decomposition      = 53000._wp                       ! Ahrens et al. 2015, SBB
    ed_decomposition      = 100000.0_wp                     ! Ahrens et al. 2015, SBB
    t_ref_decomposition   = 313.15_wp                       ! Thum et al. 2019, GMD
    k1_afps               = 4._wp / 3._wp                   ! B. Ahrens, pers. comm. 2016
    k2_afps               = 3._wp                           ! B. Ahrens, pers. comm. 2016
    kmr_depolymerisation  = 0.001_wp                        ! Thum et al. 2019, GMD ???
    kmr_mic_uptake        = 0.001_wp                        ! Thum et al. 2019, GMD ???
    kmr_hsc               = 0.001_wp                        ! Davidson et al. 2012, GCB
    kmr_sorption          = 0.001_wp                        ! Thum et al. 2019, GMD ???
    kmr_decomposition     = -2.0_wp                         ! Thum et al. 2019, GMD
    kmr_asymb_bnf         = -1.5_wp                         ! set to correspond to the permanent wilting point, SZ, 21/01/22

    ! transport and sorption parameters (all B. Ahrens, pers. comm. 2016)
    k_diff_org                = 0.15_wp / one_year / one_day                                ! Koven et al. 2013, BG [m2 (kg m-3) s-1]
    rho_bulk_org              = 150.3935_wp                                                 ! textbook knowledge, or from field measurement
    carbon_per_dryweight_SOM  = 0.4271_wp                                                   ! COMISSION for Hainich
    k_adsorpt_dom             = 59.9610_wp / one_year / one_day * molar_mass_C  / 1000._wp  ! conversion to m3 / mol C / s, Bernhard 20180913
    k_desorpt_dom             = 508.2151_wp / one_year / one_day / 1000._wp                 ! Bernhard 20180913
    k_adsorpt_det             = 0.31012_wp / one_year / one_day * molar_mass_C / 1000._wp   ! conversion to m3 / mol C / s, Bernhard 20180913
    k_desorpt_det             = 153.957_wp / one_year / one_day / 1000._wp                  ! Bernhard 20180913
    k_adsorpt_nh4             = 0.0_wp                                                      ! not activated yet
    k_desorpt_nh4             = 0.022_wp / 30.0_wp / one_day                                ! assumed to be the same as po4 desorption rate
    k_adsorpt_po4             = 0.022_wp * 8._wp / 9._wp / 30.0_wp / one_day                ! [s-1], simple_1d model, Yang et al. 2014, corrected for the global average,
                                                                                            ! assuming labile and slow P in equilibrium (Yang et al. 2013)
    k_desorpt_po4             = 0.022_wp / 30.0_wp / one_day                                ! [s-1], simple_1d model, Yang et al. 2014 * 100
    k_adsorpt_f_po4           = 12.9_wp / 3600._wp                                          ! [m3 mol-1 P s-1], based on Van der Zee et al. 1989
    k_desorpt_f_po4           = 0.14_wp / 3600._wp / 10._wp                                 ! [s-1], Van der Zee et al. 1989

    ! microbial parameters
    microbial_cn            = 7.6_wp * molar_mass_N / molar_mass_C   ! textbook knowledge
    microbial_np            = 5.6_wp * molar_mass_P / molar_mass_N   ! textbook knowledge
    microbial_nue           = 0.8_wp                                 ! Manzoni et al. 2008
    microbial_pue           = 0.8_wp                                 ! Manzoni et al. 2008
    microbial_cue           = 0.392_wp                               ! COMISSION paper
    microbial_cue_max       = 0.6_wp                                 ! Manzoni et al. 2008
    microbial_cue_min       = 0.3_wp                                 ! Manzoni et al. 2008

    ! fast and slow parameters for simple_1d SB model
    k_fast_som_cn_max       = 13._wp * molar_mass_N / molar_mass_C   ! Manzoni et al. 2008, slightly tighter than than (was 13)
    k_fast_som_cn_min       = 5._wp * molar_mass_N / molar_mass_C    ! Manzoni et al. 2008
    f_fast_som_cn           = 51000.0_wp                             ! calibrated to correspond to CENTURY (minimum reached at 2gN/m2)
    k_slow_som_cn           = 9._wp * molar_mass_N / molar_mass_C    ! textbook knowledge
    k_fast_som_np           = 14._wp * molar_mass_P / molar_mass_N   ! textbook knowledge
    kref_fast_som_np        = 0.1_wp                                 ! an average value chosen as reference level (SZ)
    k_slow_som_np           = 14._wp * molar_mass_P / molar_mass_N   ! textbook knowledge
    frac_woody2dec_litter   = 0.3_wp                                 ! ???
    frac_litter2fast_som    = 0.45_wp                                ! Parton et al. 1993, GBC
    frac_fast2slow_som      = 0.15_wp                                ! Parton et al. 1993, GBC
    frac_slow2fast_som      = 0.3_wp                                 ! Parton et al. 1993, GBC

    ! other parameters for initialisation
    k_init_soluable_litter_cn     = 25._wp * molar_mass_N / molar_mass_C    ! Thum et al. 2019, GMD ???
    k_init_polymeric_litter_cn    = 100._wp * molar_mass_N / molar_mass_C   ! Thum et al. 2019, GMD ???
    k_init_woody_litter_cn        = 260._wp * molar_mass_N / molar_mass_C   ! Thum et al. 2019, GMD ???
    k_init_litter_np              = 18._wp * molar_mass_P / molar_mass_N    ! Thum et al. 2019, GMD ???
    frac_prod_not_root            = 0.6_wp                                  ! ???
    min_fraction_default_sb       = 0.01_wp                                 ! useful assumption

    ! additional microbial parameters for JSM implementation
    fc_woody2polymeric            = 0.3_wp                                        ! Yu et al. 2020, GMD
    fc_soluable2dom               = 0.7_wp                                        ! Yu et al. 2020, GMD
    k_mort_microbial              = 7.48_wp * 1.e-8_wp                            ! COMISSION, Ahrens et al. 2015
    min_mic_biomass               = 0.01_wp                                       ! Yu et al. 2020, GMD
    min_enzyme_fraction           = 0.01_wp                                       ! Yu et al. 2020, GMD
    max_enzyme_fraction           = 0.6_wp                                        ! tuned, mostly against the very P-poor site LUE observations (in Germany)
    init_enzyme_fraction          = 0.5_wp                                        ! initial values for enzyme allocation etc.
    vmax_mic_uptake               = 3495.18053_wp / &                             ! ???
                                    one_year / one_day * 1.e6_wp * 10._wp
    km_mic_uptake                 = 1.02404788_wp / molar_mass_C * 1000._wp       ! COMISSION, Ahrens et al. 2015
    vmax_depolymerisation_residue = 0.2316605_wp / one_year / one_day             ! calibrated, from Yu et al. 2020
    vmax_depolymerisation_litter  = 0.3698_wp / one_year / one_day / 2._wp        ! calibrated, from Yu et al. 2020
    km_depolymerisation           = 0.08876_wp / 2._wp / molar_mass_C * 10._wp    ! calibrated, from Yu et al. 2020
    enzyme_frac_total             = 0.005_wp
    fc_p_res2dom                  = 0._wp                                         ! Yu et al. 2020, GMD
    fc_n_res2dom                  = 0.8_wp                                        ! Yu et al. 2020, GMD
    k_p_recyc_corr                = 0.3_wp                                        ! tuned [0,1], based on model-data comparison of 16 ecosystem traits (4 vegetation, 3 microbe, 9 soil) at 5 sites
    k_n_recyc_corr                = 0.1_wp                                        ! tuned [0,1], based on model-data comparison of 16 ecosystem traits (4 vegetation, 3 microbe, 9 soil) at 5 sites

    ! parameters controlling litter partitioning
    f_soluable_residue   = 0.172_wp                                   ! COMISSION, Ahrens et al. 2015
    f_soluable_leaf      = 0.5_wp                                     ! Not used
    f_soluable_fine_root = 0.5_wp                                     ! Not used
    f_soluable_fruit     = 0.5_wp                                     ! Not used
    f_soluable_seed_bed  = 0.5_wp                                     ! Not used
    fc_soluable_max      = 0.85_wp                                    ! Century, Parton et al. 1993
    k_fc_soluable        = 0.018_wp                                   ! Century, Parton et al. 1993
    k_fn_soluable        = 5.0_wp                                     ! Century, Parton et al. 1993
    k_fp_soluable        = 5.0_wp                                     ! Century, Parton et al. 1993
    lc_leaf_max          = 0.295_wp / molar_mass_C * molar_mass_N     ! White et al. 2000
    lc_leaf2sla          = -0.3712_wp / molar_mass_C * molar_mass_N   ! White et al. 2000
    lc_fine_root         = 0.22_wp / molar_mass_C * molar_mass_N      ! White et al. 2000
    lc_coarse_root       = 0.7_wp / molar_mass_C * molar_mass_N       ! assuming woody values
    lc_woody_litter      = 0.7_wp / molar_mass_C * molar_mass_N       ! White et al. 2000
    lc_fruit             = 0.22_wp / molar_mass_C * molar_mass_N      ! assuming fine roots values
    lc_seed_bed          = 0.22_wp / molar_mass_C * molar_mass_N      ! assuming fine roots values

    ! mycorrhizal parameters
    mycorrhiza_cue          = 0.7_wp                                  !
    mycorrhiza_cn           = 10.0_wp * molar_mass_N / molar_mass_C   !
    mycorrhiza_np           = 14.0_wp * molar_mass_P / molar_mass_N   !
    f_maxresp_myc           = 0.5_wp                                  !
    f_surface_myc2fineroot  = 40.0_wp                                 !
    myc_root_coverage_max   = 0.3_wp                                  !
    vmax_uptake_norg        = 1.0_wp                                  !
    km_up_org               = 400.0_wp                                !
    f_nin2norg              = 1.0_wp                                  !
    f_norg_limit            = 0.8_wp                                  !
    vmax_priming            = 20.0_wp                                 !
    km_priming              = 1.0_wp                                  !
    resp_priming            = 0.05_wp / one_day                       !
    organic_cn              = 3.0_wp * molar_mass_N / molar_mass_C    !
    f_tau_myc               = 5.0_wp                                  !

    ! nutrient uptake kinetics & ECA approach: simple_1d model
    km1_up_nh4 = 0.0416_wp                 ! Kronzucker et al. 1996, Plant Physiology
    km2_up_nh4 = 1.0_wp                    ! Kronzucker et al. 1996, Plant Physiology
    km1_up_no3 = 0.0416_wp                 ! Kronzucker et al. 1996, Plant Physiology
    km2_up_no3 = 1.0_wp                    ! Kronzucker et al. 1996, Plant Physiology
    km1_up_po4 = 0.0689_wp / 0.0003_wp     ! l/mol, corrected low affinity half-saturation, estimated based on Kavka and Polle 2016,
                                           ! assuming Km for low-affinity transporter is 300 micro-mol/l
    km2_up_po4 = 22.0_wp / 1.e6_wp         ! mol/l, Kavka and Polle 2016

    ! additional jsm parameters
    vmax_mic_uptake_nh4       = 106.56_wp / 3600._wp * molar_mass_C     ! micro-mol / mol / s from Zhu et al. 2017, and references therein
    km_mic_uptake_nh4         = 0.18_wp / molar_mass_N                  ! mol / m3 from Zhu et al. 2017, and references therein
    vmax_mic_uptake_no3       = 86.58_wp / 3600._wp * molar_mass_C      ! micro-mol / mol / s from Zhu et al. 2017, and references therein
    km_mic_uptake_no3         = 0.41_wp / molar_mass_N                  ! mol / m3 from Zhu et al. 2017, and references therein
    f_mic_carrier_n           = 0.5_wp / 10000._wp                      ! assuming 30% of the microbial enzymes are P transporter, based on Zhu et al. 2016
    vmax_plant_uptake_nh4     = 108.78_wp / 3600._wp * molar_mass_C     ! micro-mol / mol / s from Zhu et al. 2017, and references therein
    km_plant_uptake_nh4       = 0.896_wp / molar_mass_N                 ! mol / m3 from Zhu et al. 2017, and references therein
    vmax_plant_uptake_no3     = 18.2_wp / 3600._wp * molar_mass_C       ! micro-mol / mol / s from Zhu et al. 2017, and references therein
    km_plant_uptake_no3       = 0.672_wp / molar_mass_N                 ! mol / m3 from Zhu et al. 2017, and references therein
    f_plant_carrier_n         = 1.25_wp / 10000._wp                     ! assuming 30% of the microbial enzymes are P transporter, based on Zhu et al. 2016
    vmax_nitrification_nh4    = 636._wp / &                             ! micro-mol / mol / s , Kemp and Dodds 2002
                                molar_mass_N / one_day * molar_mass_C
    vmax_denitrification_no3  = 54._wp / &                              ! micro-mol / mol / s , Kemp and Dodds 2002
                                molar_mass_N / one_day * molar_mass_C
    km_nitrification_nh4      = 0.019_wp / molar_mass_N                 ! mol / m3 , Kemp and Dodds 2002
    km_denit_no3              = 0.3_wp / molar_mass_N                   ! mol / m3 , Kemp and Dodds 2002
    k_nit_carrier_n           = 0.001_wp / 1000._wp                     ! estimated
    vmax_mic_uptake_po4       = 15.72_wp / 3600._wp * molar_mass_C      ! micro-mol / mol / s, based on Zhu et al. 2016, and references therein
    km_mic_uptake_po4         = 0.02_wp / molar_mass_P                  ! mol / m3, Zhu et al. 2016
    vmax_plant_uptake_po4     = 0.0222 / 60._wp * molar_mass_C          ! micro-mol / mol / s, Kavka and Polle 2016
    km_plant_uptake_po4       = 0.067_wp / molar_mass_P                 ! mol / m3
    k_mic_uptake_nh4          = 0.32_wp * 24._wp / one_day              ! [s-1], Zhu et al. 2016; 0.32 h-1, Zhu et al. 2017
    k_plant_uptake_nh4        = 10.8_wp * 24._wp / one_day              ! [s-1], Zhu et al. 2016; 10.8 h-1, Zhu et al. 2017
    k_mic_uptake_no3          = 0.26_wp * 24._wp / one_day              ! [s-1], Zhu et al. 2016; 0.26 h-1, Zhu et al. 2017
    k_plant_uptake_no3        = 8.9_wp * 24._wp / one_day               ! [s-1], Zhu et al. 2016; 8.9 h-1, Zhu et al. 2017
    k_mic_uptake_po4          = 400._wp / one_day                       ! [s-1], Zhu et al. 2016, and references therein
    k_plant_uptake_po4        = 12._wp / one_day                        ! [s-1], Zhu et al. 2016, and references therein
    f_mic_carrier_po4         = 0.05_wp / 100._wp                       ! assuming 30% of the microbial enzymes are P transporter, based on Zhu et al. 2016
    f_plant_carrier_po4       = 1.25_wp / 10000._wp                     ! assuming 10% of the plant enzymes are P transporter, based on Zhu et al. 2016
    f_nit_carrier_nh4         = 0.05_wp / 10000._wp                     ! estimated
    f_denit_carrier_no3       = 0.05_wp / 10000._wp                     ! estimated
    qmax_nh4_clay             = 7.105 / molar_mass_N                    ! [mol N kg-1 clay], Nieder et al. 2011
    km_adsorpt_mineral_nh4    = 2.602 / 1000._wp / molar_mass_N         ! [mol N kg-1 fine mineral soil], based on 12.3 [mg N l-1 soil solution]  2.602_wp
                                                                        ! from literature data analysis, with .33 m wc and 1560 kg/m3 bulk density
    km_adsorpt_OM_nh4         = 26.02 / 1000._wp / molar_mass_N         ! [mol N kg-1 fine mineral soil], 10 times of mineral km
    qmax_po4_mineral          = 1.2_wp / molar_mass_P                   ! [mol P kg-1 soil], Olander and Vitousek 2005
    km_adsorpt_mineral_po4    = 0.070_wp / 1000._wp / molar_mass_P      ! [mol P kg-1 fine mineral soil], based on 0.332 [mg P l-1 soil solution]
                                                                        ! from literature data analysis, with .33 m wc and 1560 kg/m3 bulk density
    qmax_po4_clay             = 0.283159_wp / molar_mass_P              ! [mol P kg-1 clay], literature review, a analogous predictor of Feox+Alox																	
    qmax_po4_silt             = 0.283159_wp / molar_mass_P              ! [mol P kg-1 silt], literature review, clay+silt
    qmax_po4_sand             = 0.140475_wp / molar_mass_P              ! [mol P kg-1 sand], literature review
    qmax_po4_OM               = 0.140475_wp / molar_mass_P              ! [mol P kg-1 OM], literature review,
    km_po4_ph                 = 0.4_wp                                  ! corrector for km_po4 calculation [--]
    km_adsorpt_clay_po4       = 315.445_wp / 1000000._wp / molar_mass_P ! [mol P L-1 kg-1 clay], from literature data analysis
    km_adsorpt_silt_po4       = 2.176_wp / 1000._wp / molar_mass_P      ! [mol P L-1 kg-1 silt], from literature data analysis
    km_adsorpt_sand_po4       = 20._wp / 1000._wp / molar_mass_P        ! [mol P L-1 kg-1 sand], from literature data analysis																																	
    km_adsorpt_OM_po4         = 14._wp / 1000._wp / molar_mass_P        ! [mol P L-1 kg-1 OM], from literature data analysis
    k_deep_soil_moc_corr      = 0.01_wp                                 ! estimated, Yu et al. 2020, GMD
    fact_repartition_om_init_cond1 = 1.0_wp / 10.0_wp                   ! estimated, Yu et al. 2020, GMD
    fact_repartition_om_init_cond2 = 1.0_wp / 11.0_wp                   ! estimated, Yu et al. 2020, GMD

    ! nitrification, denitrification parameters (sb_nloss_scheme = dynamic)
    ea_nitrification          = 80000._wp                               ! Xu-Ri & Prentice, 2008, to fit their function
    ed_nitrification          = 200000._wp                              ! Xu-Ri & Prentice, 2008, to fit their function
    t_opt_nitrification       = 311.15_wp                               ! Xu-Ri & Prentice, 2008, to fit their function
    f_nit_noy                 = 0.02_wp                                 ! Xu-Ri & Prentice, 2008
    f_nit_n2o                 = 0.002_wp                                ! Xu-Ri & Prentice, 2008
    ea_denitrification        = 47000._wp                               ! Xu-Ri & Prentice, 2008, to fit their function
    vmax_nitrification        = 0.4_wp / one_day * 1.e6_wp              ! Xu-Ri & Prentice, 2008, 10% / day at 20degC
    vmax_denitrification      = 0.1_wp / one_day * 1.e6_wp              ! Xu-Ri & Prentice, 2008, 10% / day
    km_denitrification_c      = 20.0_wp                                 ! ???
    km_denitrification_no3    = 0.083_wp * 1000.0_wp / molar_mass_N     ! Xu-Ri & Prentice, 2008
    f_denit_noy               = 0.002_wp                                ! Xu-Ri & Prentice, 2008
    f_denit_n2o               = 0.02_wp                                 ! Xu-Ri & Prentice, 2008
    ea_gasdiffusion           = 47000._wp                               ! Xu-Ri & Prentice, 2008, to fit their function
    lambda_anvf               = 1.3_wp                                  ! fitted from Zaehle et al. 2011, NatGeo
    k_anvf                    = 3.0_wp                                  ! fitted from Zaehle et al. 2011, NatGeo

    ! simple_1d model fixed loss and simple_1d leaching parameterisation
    floss_nmin            = 0.05_wp                                 ! ???
    fnit_nmin             = 0.1_wp                                  ! ???
    ! simple_1d leaching
    fleach_nh4            = 0.1_wp                                  ! ???
    fleach_po4            = 0.001_wp                                ! ???
    ! asymbiotic BNF parameters
    km_asymb_bnf          = 0.05_wp                                 ! calibrated to match Davies-Barnard, 2020
    vmax_asymb_bnf        = 0.03_wp                                 ! calibrated to match Davies-Barnard, 2020

    ! Isoptope stuff
    eta_mic_uptake_nh4        = 17._wp                           ! Robinson 2001, TREE
    eta_mic_uptake_no3        = 13._wp                           ! Robinson 2001, TREE
    eta_plant_uptake_nh4      = 13.5_wp                          ! Robinson 2001, TREE
    eta_plant_uptake_no3      = 9.5_wp                           ! Robinson 2001, TREE
    eta_nitrification         = 47.5_wp                          ! Robinson 2001, TREE
    eta_nitrate_production    = 25._wp                           ! Robinson 2001, TREE
    eta_denitrification       = 31._wp                           ! Robinson 2001, TREE
    eta_ammonification        = 2.5_wp                           ! Robinson 2001, TREE
    eta_volatilisation_nh3    = 50._wp                           ! Robinson 2001, TREE

    ! soil P fluxes
    k_weath_mineral                        = 0.1_wp / 6.65_wp / &                          ! Wang et al. 2010 combining with BBR soil weight of 1m    [s-1]
                                             molar_mass_P / one_year / one_day
    frac_background_weath_mineral          = 0.1_wp                                        ! assumption based on expert knowledge
    eq_sorp_mineral_po4                    = 80._wp / molar_mass_P                         ! random guess based on init value of assoc_slow    [mol P m-3]
    k_occlude_po4                          = 1.0_wp / 1000000._wp / 30._wp / one_day       ! Yang et al. 2014            [s-1]
    vmax_biochem_po4                       = 0.25_wp / one_day * 5._wp                     ! estimated from Bunemann et al. 2016, LUS site, based on 50 mol/m3 mic   [s-1]
    km_mic_biochem_po4                     = 5._wp / molar_mass_C * 0.5_wp                 ! random number, equals to the usual microbe biomass (0.1 g C m-3)   [mol C m-3]
    km_dip_biochem_po4                     = 0.001_wp                                      ! estimated based on expert knowledge and model testing
    km_pc_biochem_po4                      = 0.01_wp                                       ! estimated based on expert knowledge and model testing
    km_rootc_biochem_po4                   = 20._wp                                        ! estimated based on expert knowledge and model testing
    km_rootc_weath_po4                     = 10._wp                                        ! estimated based on expert knowledge and model testing
    km_root_enzyme_weath_po4               = 1._wp / molar_mass_C / 10._wp                 ! estimated based on expert knowledge and model testing
    biomin_threshold_pc_ratio_mic_res      = 1._wp / 5000._wp                              ! estimated based on expert knowledge and model testing
    biomin_threshold_pc_ratio_min_assoc_om = 1._wp / 50000._wp                             ! estimated based on expert knowledge and model testing, 1/10 biomin_threshold_pc_ratio_mic_res

    ! soil P pool
    k1_p_pool_depth_corr              = 3.5_wp          ! calibrated
    k2_p_pool_depth_corr              = 2._wp           ! calibrated
    k3_p_pool_depth_corr              = 0.8_wp          ! calibrated
    k4_p_pool_depth_corr              = 0.5_wp          ! calibrated
    k_p_pool_depth_corr_min_occlud    = 0.1_wp          ! calibrated
    p_labile_slow_global_avg          = 9._wp / 17._wp  ! Yang et al. 2013
    k_sb_pool_dom_assoc_c_corr        = 5.0_wp          ! calibrated

    ! sb_pool default init values
    def_sb_pool_metabol_litter        = 100.0_wp        ! Thum et al. 2019, GMD
    def_sb_pool_simple_struct_litter  = 100.0_wp        ! Thum et al. 2019, GMD
    def_sb_pool_jsm_struct_litter     = 1000.0_wp       ! Yu et al. 2020, GMD
    def_sb_pool_woody_litter          = 1000.0_wp       ! Thum et al. 2019, GMD
    def_sb_pool_fast_som              = 1500.0_wp       ! Thum et al. 2019, GMD
    def_sb_pool_slow_som              = 10000.0_wp      ! Thum et al. 2019, GMD
    def_sb_pool_microbial_biomass     = 30.0_wp         ! Yu et al. 2020, GMD
    def_sb_pool_microbial_necromass   = 150.0_wp        ! Yu et al. 2020, GMD
    def_sb_pool_dom                   = 1.0_wp          ! Yu et al. 2020, GMD

    ! constants used as factor for multiplication with namelist parameters (spq_config) when "parameter sensitivity" is enabled
    f_psensi_soil_p_depth             = 1.0_wp
    f_psensi_soil_p_labile            = 1.0_wp
    f_psensi_soil_p_slow              = 1.0_wp
    f_psensi_soil_p_occluded          = 1.0_wp
    f_psensi_soil_p_primary           = 1.0_wp
    f_psensi_qmax_org_fine_particle   = 1.0_wp


    !! BLARPP only temporarily, mol / m2 !!  @TODO
    nh4_solute_prescribe = 0.2_wp         ! tuned, to have enough nh4 for microbial dynamics and a large enough associated-nh4 pool for initialization
    no3_solute_prescribe = 0.1_wp
    po4_solute_prescribe = 0.05_wp
    po4_solute_min       = 1.e-6_wp


  END SUBROUTINE init_sb_constants

#endif
END MODULE mo_sb_constants
