!> QUINCY soil-biogeochemistry process memory
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
!>#### definition and init of (memory) variables for the soil-biogeochemistry process
!>
MODULE mo_sb_memory_class
#ifndef __NO_QUINCY__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
  USE mo_util,                   ONLY: One_of, int2string
  USE mo_jsb_utils_iface,        ONLY: assign_if_present_allocatable

  USE mo_jsb_process_class,      ONLY: SB_
  USE mo_lnd_bgcm_class,         ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, &
    &                                  LAST_ELEM_ID, get_name_for_element_id, get_short_name_for_element_id
  USE mo_lnd_bgcm_store_class,   ONLY: SB_BGCM_POOL_ID, SB_BGCM_DELTA_POOLS_ID, SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, &
    &                                  SB_BGCM_TRANSPORT_ID, SB_BGCM_MYCO_EXPORT_ID
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_pool_class,         ONLY: t_jsb_pool, ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_memory, max_no_of_vars
  PUBLIC :: t_sb_bgcm_main, t_sb_bgcm_with_components, t_sb_bgcm_with_elements

  INTEGER, PARAMETER :: max_no_of_vars = 1000

  ! ======================================================================================================= !
  !>Type definition for sb memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_sb_memory

    !------------------------------------------------------------------------------------------------------ !
    !> bgc_material structures
    !>
    TYPE(t_sb_bgcm_main), POINTER :: sb_bgcm_main

    !------------------------------------------------------------------------------------------------------ !
    !> bgc_material diagnostics - summary of C, N, P, C13, C14, N15 across the particular bgcm
    !> used for output in IQ
    !> because the bgcm cannot be used with the 'aggregate()' routine yet @TODO
    !>
    TYPE(t_jsb_var_real3d)  :: &
      & sb_pool_total_c          , &        !< sum of all sb_pool%*%carbon      [mol m-3]
      & sb_pool_total_n          , &        !< sum of all sb_pool%*%nitrogen    [mol m-3]
      & sb_pool_total_p          , &        !< sum of all sb_pool%*%phosphorus  [mol m-3]
      & sb_pool_total_c13        , &        !< sum of all sb_pool%*%carbon13    [mol m-3]
      & sb_pool_total_c14        , &        !< sum of all sb_pool%*%carbon14    [mol m-3]
      & sb_pool_total_n15                   !< sum of all sb_pool%*%nitrogen15  [mol m-3]
    ! Sum of
    ! soil organic matter (som): dom & dom_assoc & fungi & mycorrhiza & microbial & residue & residue_assoc
    ! and below ground (bg) litter: bg parts of soluble_litter & polymeric_litter & woody_litter
    TYPE(t_jsb_var_real3d) ::  &
      & sb_pool_total_bg_soil_c,   & !< sum of som and below ground litter (layers 2:nsoil) of sb_pool%*%carbon     [mol m-3]
      & sb_pool_total_bg_soil_n,   & !< sum of som and below ground litter (layers 2:nsoil) of sb_pool%*%nitrogen   [mol m-3]
      & sb_pool_total_bg_soil_p,   & !< sum of som and below ground litter (layers 2:nsoil) of sb_pool%*%phosphorus [mol m-3]
      & sb_pool_total_bg_soil_c13, & !< sum of som and below ground litter (layers 2:nsoil) of sb_pool%*%carbon13   [mol m-3]
      & sb_pool_total_bg_soil_c14, & !< sum of som and below ground litter (layers 2:nsoil) of sb_pool%*%carbon14   [mol m-3]
      & sb_pool_total_bg_soil_n15    !< sum of som and below ground litter (layers 2:nsoil) of sb_pool%*%nitrogen15 [mol m-3]
    ! Sum of above ground (ag) litter: ag parts of soluble_litter & polymeric_litter & woody_litter
    TYPE(t_jsb_var_real2d) ::        &
      & sb_pool_total_ag_litter_c,   & !< sum of above ground litter (soil layer 1) of sb_pool%*%carbon     [mol m-2]
      & sb_pool_total_ag_litter_n,   & !< sum of above ground litter (soil layer 1) of sb_pool%*%nitrogen   [mol m-2]
      & sb_pool_total_ag_litter_p,   & !< sum of above ground litter (soil layer 1) of sb_pool%*%phosphorus [mol m-2]
      & sb_pool_total_ag_litter_c13, & !< sum of above ground litter (soil layer 1) of sb_pool%*%carbon13   [mol m-2]
      & sb_pool_total_ag_litter_c14, & !< sum of above ground litter (soil layer 1) of sb_pool%*%carbon14   [mol m-2]
      & sb_pool_total_ag_litter_n15    !< sum of above ground litter (soil layer 1) of sb_pool%*%nitrogen15 [mol m-2]
    ! various bgcm diagnostics
    TYPE(t_jsb_var_real3d)  :: &
      & sb_pool_woody_litter_c         !< woody litter or sb_pool [mol m-3]

    !------------------------------------------------------------------------------------------------------ !
    !>pools 3D
    !>
    TYPE(t_jsb_var_real3d)  :: &
                            nh4_solute                , &  !< ammonium in soil solution [mol m-3]
                            nh4_assoc                 , &  !< minerally associated ammonium [mol m-3]
                            no3_solute                , &  !< nitrate in soil solution [mol m-3]
                            noy                       , &  !< NO_y [mol m-3]
                            n2o                       , &  !< N2O [mol m-3]
                            n2                        , &  !< N2 [mol m-3]
                            nh4_n15_solute            , &  !< N15-ammonium in soil solution [mol m-3]
                            nh4_n15_assoc             , &  !< minerally associated N15-ammonium [mol m-3]
                            no3_n15_solute            , &  !< N15-nitrate in soil solution [mol m-3]
                            noy_n15                   , &  !< N15-NO_y [mol m-3]
                            n2o_n15                   , &  !< N15-N2O [mol m-3]
                            n2_n15                    , &  !< N15-N2 [mol m-3]
                            po4_assoc_fast            , &  !< fast minerally associated phosphate  [mol m-3]
                            po4_assoc_slow            , &  !< slow minerally associated phosphate  [mol m-3]
                            po4_occluded              , &  !< occluded phosphate [mol m-3]
                            po4_primary               , &  !< primary phosphate mineral [mol m-3]
                            po4_solute                , &  !< phophate in soil solution [mol m-3]
                            soil_litter_carbon_sl          !< (bulk) soil density of litter carbon [mol C m-3 soil]

    !------------------------------------------------------------------------------------------------------ !
    !>transport fluxes 3D
    !>
    TYPE(t_jsb_var_real3d)  :: &
                            transport_nh4_solute         , &  !< transport of ammonium in soil solution [mol m-3 timestep-1]
                            transport_nh4_assoc          , &  !< transport of minerally associated ammonium [mol m-3 timestep-1]
                            transport_no3_solute         , &  !< transport of nitrate in soil solution [mol m-3 timestep-1]
                            transport_noy                , &  !< transport of NO_y [mol m-3 timestep-1]
                            transport_n2o                , &  !< transport of N2O [mol m-3 timestep-1]
                            transport_n2                 , &  !< transport of N2 [mol m-3 timestep-1]
                            transport_nh4_n15_solute     , &  !< transport of N15-ammonium in soil solution [mol m-3 timestep-1]
                            transport_nh4_n15_assoc      , &  !< transport of minerally associated N15-ammonium [mol m-3 timestep-1]
                            transport_no3_n15_solute     , &  !< transport of N15-nitrate in soil solution [mol m-3 timestep-1]
                            transport_noy_n15            , &  !< transport of N15-NO_y [mol m-3 timestep-1]
                            transport_n2o_n15            , &  !< transport of N15-N2O [mol m-3 timestep-1]
                            transport_n2_n15             , &  !< transport of N15-N2 [mol m-3 timestep-1]
                            transport_po4_assoc_fast     , &  !< transport of fast minerally associated phosphate  [mol m-3 timestep-1]
                            transport_po4_assoc_slow     , &  !< transport of slow minerally associated phosphate  [mol m-3 timestep-1]
                            transport_po4_occluded       , &  !< transport of occluded phosphate [mol m-3 timestep-1]
                            transport_po4_primary        , &  !< transport of primary phosphate [mol m-3 timestep-1]
                            transport_po4_solute              !< transport of minerally associated phosphate  [mol m-3 timestep-1]

    !------------------------------------------------------------------------------------------------------ !
    !>lateral loss fluxes 3D
    !>
    !>   DOM indicate dissolved organic matter
    !>
    TYPE(t_jsb_var_real3d)  :: &
                            lateral_loss_dom_carbon_sl      , & !< land to groundwater flux of DOM-C [mol m-3 timestep-1]
                            lateral_loss_dom_nitrogen_sl    , & !< land to groundwater flux of DOM-N [mol m-3 timestep-1]
                            lateral_loss_dom_phosphorus_sl  , & !< land to groundwater flux of DOM-P [mol m-3 timestep-1]
                            lateral_loss_dom_carbon13_sl    , & !< land to groundwater flux of DOM-C13 [mol m-3 timestep-1]
                            lateral_loss_dom_carbon14_sl    , & !< land to groundwater flux of DOM-C14 [mol m-3 timestep-1]
                            lateral_loss_dom_nitrogen15_sl  , & !< land to groundwater flux of DOM-N15 [mol m-3 timestep-1]
                            lateral_loss_nh4_solute_sl      , & !< land to groundwater flux of NH4 [mol m-3 timestep-1]
                            lateral_loss_no3_solute_sl      , & !< land to groundwater flux of NO3 [mol m-3 timestep-1]
                            lateral_loss_po4_solute_sl      , & !< land to groundwater flux of PO4 [mol m-3 timestep-1]
                            lateral_loss_nh4_n15_solute_sl  , & !< land to groundwater flux of 15NH4 [mol m-3 timestep-1]
                            lateral_loss_no3_n15_solute_sl      !< land to groundwater flux of 15NO3 [mol m-3 timestep-1]
    !------------------------------------------------------------------------------------------------------ !
    !>response functions 3D
    !>  temperature response functions
    !>
    TYPE(t_jsb_var_real3d)  :: &
                            rtm_decomposition, &              !< rate modifier for decomposition processes [unitless]
                            rtm_depolymerisation, &           !< rate modifier for depolymerisation processes [unitless]
                            rtm_mic_uptake, &                 !< rate modifier for microbial uptake [unitless]
                            rtm_plant_uptake, &               !< rate modifier for plant uptake [unitless]
                            rtm_sorption, &                   !< rate modifier for sorption [unitless]
                            rtm_desorption, &                 !< rate modifier for desorption [unitless]
                            rtm_hsc, &                        !< rate modifier for half-saturation constants [unitless]
                            rtm_nitrification, &              !< rate modifier for nitrification [unitless]
                            rtm_denitrification, &            !< rate modifier for denitrification [unitless]
                            rtm_gasdiffusion, &               !< rate modifier for gasdiffusion [unitless]
                            rtm_asymb_bnf                     !< rate modifier for n fixation [unitless]

    !>  moisture response functions
    !>
    TYPE(t_jsb_var_real3d)  :: &
                            rmm_decomposition, &              !< rate modifier for decomposition processes [unitless]
                            rmm_depolymerisation, &           !< rate modifier for depolymerisation processes [unitless]
                            rmm_mic_uptake, &                 !< rate modifier for microbial uptake [unitless]
                            rmm_plant_uptake, &               !< rate modifier for plant uptake [unitless]
                            rmm_sorption, &                   !< rate modifier for sorption [unitless]
                            rmm_desorption, &                 !< rate modifier for desorption [unitless]
                            rmm_hsc, &                        !< rate modifier for half-saturation constants [unitless]
                            rmm_nitrification, &              !< rate modifier for nitrification [unitless]
                            rmm_gasdiffusion, &               !< rate modifier for gasdiffusion [unitless]
                            rmm_asymb_bnf                     !< rate modifier for n fixation [unitless]

    !------------------------------------------------------------------------------------------------------ !
    !>fluxes 3D | respiration, mineralisation, uptake, nitrification, denitrification, ...
    !>
    TYPE(t_jsb_var_real3d)  :: &
                            het_respiration, &                    !< heterotrophic respiration  [micro-mol m-3 s-1]
                            het_respiration_c13, &                !< c13 heterotrophic respiration   [micro-mol m-3 s-1]
                            het_respiration_c14, &                !< c14 heterotrophic respiration   [micro-mol m-3 s-1]
                            myc_respiration, &                    !< mycorrhizae respiration  [micro-mol m-3 s-1]
                            myc_respiration_c13, &                !< c13 mycorrhizae respiration   [micro-mol m-3 s-1]
                            myc_respiration_c14, &                !< c14 mycorrhizae respiration   [micro-mol m-3 s-1]
                            net_mineralisation_nh4, &             !< net mineralisation of N    [micro-mol m-3 s-1]
                            net_mineralisation_nh4_n15, &         !< net mineralisation of N15  [micro-mol m-3 s-1]
                            net_mineralisation_no3, &             !< net mineralisation of N    [micro-mol m-3 s-1]
                            net_mineralisation_no3_n15, &         !< net mineralisation of N15  [micro-mol m-3 s-1]
                            net_mineralisation_po4, &             !< net mineralisation of P  [micro-mol m-3 s-1]
                            plant_uptake_nh4_sl, &                !< plant uptake of N    [micro-mol m-3 s-1]
                            plant_uptake_nh4_n15_sl, &            !< plant uptake of N15  [micro-mol m-3 s-1]
                            plant_uptake_no3_sl, &                !< plant uptake of N    [micro-mol m-3 s-1]
                            plant_uptake_no3_n15_sl, &            !< plant uptake of N15  [micro-mol m-3 s-1]
                            plant_uptake_po4_sl, &                !< plant uptake of P  [micro-mol m-3 s-1]
                            microbial_uptake_nh4_sl, &            !< microbial uptake of N    [micro-mol m-3 s-1]
                            microbial_uptake_nh4_n15_sl, &        !< microbial uptake of N15  [micro-mol m-3 s-1]
                            microbial_uptake_no3_sl, &            !< microbial uptake of N    [micro-mol m-3 s-1]
                            microbial_uptake_no3_n15_sl, &        !< microbial uptake of N15  [micro-mol m-3 s-1]
                            microbial_uptake_po4_sl, &            !< microbial uptake of P  [micro-mol m-3 s-1]
                            mycorrhiza_uptake_nh4_sl, &           !< mycorrhiza uptake of N    [micro-mol m-3 s-1]
                            mycorrhiza_uptake_nh4_n15_sl, &       !< mycorrhiza uptake of N15  [micro-mol m-3 s-1]
                            mycorrhiza_uptake_no3_sl, &           !< mycorrhiza uptake of N    [micro-mol m-3 s-1]
                            mycorrhiza_uptake_no3_n15_sl, &       !< mycorrhiza uptake of N15  [micro-mol m-3 s-1]
                            mycorrhiza_uptake_norg_sl, &          !< mycorrhiza uptake of N    [micro-mol m-3 s-1]
                            mycorrhiza_uptake_norg_n15_sl, &      !< mycorrhiza uptake of N15  [micro-mol m-3 s-1]
                            mycorrhiza_uptake_po4_sl, &           !< mycorrhiza uptake of P  [micro-mol m-3 s-1]
                            myc_export_n_tlabile_mavg_sl, &       !< time averaged export from mycorrhizae to plants of N [mol m-3]
                            myc_export_p_tlabile_mavg_sl, &       !< time averaged export from mycorrhizae to plants of P [mol m-3]
                            myc_export_c_tmyc_mavg_sl, &          !< time averaged export from mycorrhizae to plants of C [mol m-3]
                            myc_export_n_tmyc_mavg_sl, &          !< time averaged export from mycorrhizae to plants of N [mol m-3]
                            volatilisation_nh4, &                 !< volatilisation of NH4  [micro-mol m-3 s-1]
                            volatilisation_nh4_n15, &             !< volatilisation of 15NH4  [micro-mol m-3 s-1]
                            nitrification_no3, &                  !< nitrification NO3 [micro-mol m-3 s-1]
                            nitrification_no3_n15, &              !< nitrification 15NO3 [micro-mol m-3 s-1]
                            nitrification_noy, &                  !< nitrification to NOy  [micro-mol m-3 s-1]
                            nitrification_noy_n15, &              !< nitrification to 15NOy  [micro-mol m-3 s-1]
                            nitrification_n2o, &                  !< nitrification to N2O [micro-mol m-3 s-1]
                            nitrification_n2o_n15, &              !< nitrification to 15N2O [micro-mol m-3 s-1]
                            denitrification_noy, &                !< nitrification NOy [micro-mol m-3 s-1]
                            denitrification_n2o_n15, &            !< nitrification to 15N2O  [micro-mol m-3 s-1]
                            denitrification_n2o, &                !< nitrification to N2O  [micro-mol m-3 s-1]
                            denitrification_noy_n15, &            !< nitrification 15NOy [micro-mol m-3 s-1]
                            denitrification_n2, &                 !< denitrification to N2  [micro-mol m-3 s-1]
                            denitrification_n2_n15, &             !< denitrification to 15N2  [micro-mol m-3 s-1]
                            asymb_n_fixation, &                   !< asymbiotic N fixation [micro-mol m-3 s-1]
                            asymb_n_fixation_n15, &               !< asymbiotic 15N fixation [micro-mol m-3 s-1]
                            weathering_po4, &                     !< phosphorus weathering rate [micro-mol P m-3 s-1]
                            occlusion_po4, &                      !< phosphorus occlusion rate, assoc_slow --> occluded [micro-mol P m-3 s-1]
                            slow_exchange_po4, &                  !< phosphorus slow exchange rate, assoc_fast <--> assoc_slow [micro-mol P m-3 s-1]
                            fast_exchange_po4, &                  !< phosphorus fast exchange rate, solute <--> assoc_fast [micro-mol P m-3 s-1]
                            fast_exchange_nh4, &                  !< nh4 fast exchange rate, solute <--> nh4_assoc [micro-mol N m-3 s-1]
                            fast_adsorpt_po4, &                   !< phosphorus fast exchange rate, solute <--> assoc_fast [micro-mol P m-3 s-1]
                            fast_adsorpt_nh4, &                   !< nh4 fast adsorption rate, solute --> nh4_assoc [micro-mol P m-3 s-1]
                            biochem_mineralisation_po4, &         !< phosphorus biochemical mineralization rate [micro-mol P m-3 s-1]
                            gross_mineralisation_po4              !< phosphorus gross (biological) mineralization rate [micro-mol P m-3 s-1]

    !------------------------------------------------------------------------------------------------------ !
    !>2D transport fluxes for communication with groundwater and atmosphere
    !>
    TYPE(t_jsb_var_real2d)  :: &
                            leaching_dom_carbon, &                !< land to groundwater flux of DOM-C [micro-mol m-2 s-1]
                            leaching_dom_nitrogen, &              !< land to groundwater flux of DOM-N [micro-mol m-2 s-1]
                            leaching_dom_phosphorus, &            !< land to groundwater flux of DOM-P [micro-mol m-2 s-1]
                            leaching_dom_carbon13, &              !< land to groundwater flux of DOM-C13 [micro-mol m-2 s-1]
                            leaching_dom_carbon14, &              !< land to groundwater flux of DOM-C14 [micro-mol m-2 s-1]
                            leaching_dom_nitrogen15, &            !< land to groundwater flux of DOM-N15 [micro-mol m-2 s-1]
                            leaching_nh4_solute, &                !< land to groundwater flux of NH4 [micro-mol m-2 s-1]
                            leaching_no3_solute, &                !< land to groundwater flux of NO3 [micro-mol m-2 s-1]
                            leaching_po4_solute, &                !< land to groundwater flux of PO4 [micro-mol m-2 s-1]
                            leaching_nh4_n15_solute, &            !< land to groundwater flux of 15NH4 [micro-mol m-2 s-1]
                            leaching_no3_n15_solute, &            !< land to groundwater flux of 15NO3 [micro-mol m-2 s-1]
                            emission_noy, &                       !< land to atmosphere flux of NOy [micro-mol m-2 s-1]
                            emission_n2o, &                       !< land to atmosphere flux of N2O [micro-mol m-2 s-1]
                            emission_n2, &                        !< land to atmosphere flux of N2 [micro-mol m-2 s-1]
                            emission_noy_n15, &                   !< land to atmosphere flux of 15NOy [micro-mol m-2 s-1]
                            emission_n2o_n15, &                   !< land to atmosphere flux of 15N2O [micro-mol m-2 s-1]
                            emission_n2_n15                       !< land to atmosphere flux of 15N2 [micro-mol m-2 s-1]

    !------------------------------------------------------------------------------------------------------ !
    !>3D ...
    !>
    TYPE(t_jsb_var_real3d)  :: &
                            qmax_org, &                           !< sorption capacity of dom of mineral soil [kg C kg-1 mineral soil]
                            qmax_po4, &                           !< sorption capacity of po4 of mineral soil [mol P m-3 soil]
                            qmax_nh4, &                           !< sorption capacity of nh4 of mineral soil [mol N m-3 soil]
                            qmax_fast_po4, &                      !< sorption capacity of po4 of the Fe/Al pool [mol P m-3 soil]
                            qmax_slow_po4, &                      !< sorption capacity of po4 of the clay/silt pool [mol P m-3 soil]
                            km_fast_po4, &                        !< half-saturation concentration of po4 for the Fe/Al pool [mol P m-3 soil]
                            km_slow_po4, &                        !< half-saturation concentration of po4 for the clay/silt pool [mol P m-3 soil]
                            km_adsorpt_po4_sl, &                  !< half saturation of po4 sorption isotherm [mol P m-3 soil]
                            km_adsorpt_nh4_sl, &                  !< half saturation of nh4 sorption isotherm [mol N m-3 soil]
                            k_bioturb, &                          !< diffusion rate resulting from bioturbation [m2 s-1]
                            bulk_dens_corr_sl, &                  !< organic fraction corrected bulk-density of the soil [kg m-3]
                            bulk_soil_carbon_sl, &                !< bulk soil density of carbon [kg m-3]
                            particle_fluxrate, &                  !< burial rate due to organic matter fall and decay [m s-1]
                            anaerobic_volume_fraction_sl, &       !< anaerobic volume fraction [unitless]
                            microbial_cue_eff, &                  !< instantaneous microbial CUE given NP constraints [unitless]
                            microbial_cue_eff_tmic_mavg, &        !< time-averaged microbial CUE given NP constraints [unitless]
                            microbial_nue_eff, &                  !< instantaneous microbial NUE given NP constraints [unitless]
                            microbial_nue_eff_tmic_mavg, &        !< time-averaged microbial NUE given NP constraints [unitless]
                            microbial_pue_eff, &                  !< instantaneous microbial PUE given NP constraints [unitless]
                            microbial_pue_eff_tmic_mavg, &        !< time-averaged microbial PUE given NP constraints [unitless]
                            fact_n_status_mic_c_growth, &         !< microbial N status, used to scale up/down certain microbial
                                                                  !< processes (enzyme allocation and nutrient recycling) [unitless]
                            fact_p_status_mic_c_growth, &         !< microbial P status, used to scale up/down certain microbial
                                                                  !< processes (enzyme allocation and nutrient recycling) [unitless]
                            enzyme_frac_AP, &                     !< fraction of Acid (alkaline) phosphatase of the total enzyme [unitless]
                            dom_cn, &                             !< C:N ratio of DOM pool [mol mol-1]
                            dom_cp, &                             !< C:P ratio of DOM pool [mol mol-1]
                            dom_cn_mavg, &                        !< time-averaged C:N ratio of DOM pool [mol mol-1]
                            dom_cp_mavg, &                        !< time-averaged C:P ratio of DOM pool [mol mol-1]
                            residue_som_c_form_mavg_sl, &         !< time-averaged residual som C formation rate [micro-mol m-3 s-1]
                            residue_som_c14_form_mavg_sl, &       !< time-averaged residual som C14 formation rate [micro-mol m-3 s-1]
                            residue_assoc_som_c_form_mavg_sl, &   !< time-averaged associated residual som C formation rate [micro-mol m-3 s-1]
                            residue_assoc_som_c14_form_mavg_sl, & !< time-averaged associated residual som C formation rate [micro-mol m-3 s-1]
                            assoc_dom_c_form_mavg_sl, &           !< time-averaged associated dom C formation rate [micro-mol m-3 s-1]
                            assoc_dom_c14_form_mavg_sl, &         !< time-averaged associated dom C formation rate [micro-mol m-3 s-1]
                            residue_som_c_loss_mavg_sl, &         !< time-averaged residual som C loss rate [micro-mol m-3 s-1]
                            residue_assoc_som_c_loss_mavg_sl, &   !< time-averaged residual associated som C loss rate [micro-mol m-3 s-1]
                            assoc_dom_c_loss_mavg_sl              !< time-averaged associated dom C loss rate [micro-mol m-3 s-1]



    !------------------------------------------------------------------------------------------------------ !
    !>3D soil properties
    !>
    TYPE(t_jsb_var_real3d)  :: &
      & ph_sl, &                   !< pH per soil layer [unitless]
      & vmax_weath_mineral_sl, &   !< weathering coefficient of mineral soil [mol P m-3 s-1]
      & Qmax_AlFe_cor, &           !< correction factor for phosphate Qmax due to the unknown Al/Fe [unitless]
      & enzyme_frac_poly, &        !< fraction of depolymerisation enzyme which depolymerises poly litter [unitless]
      & enzyme_frac_residue, &     !< fraction of depolymerisation enzyme which depolymerises microbial residue [unitless]
      & enzyme_frac_poly_c, &      !< fraction of depolymerisation enzyme which depolymerises poly litter, determined by C [unitless]
      & enzyme_frac_poly_n, &      !< fraction of depolymerisation enzyme which depolymerises poly litter, determined by N [unitless]
      & enzyme_frac_poly_p,&       !< fraction of depolymerisation enzyme which depolymerises poly litter, determined by P [unitless]
      & enzyme_frac_poly_c_mavg, & !< moving average of enzyme_frac_poly_c [unitless]
      & enzyme_frac_poly_n_mavg, & !< moving average of enzyme_frac_poly_n [unitless]
      & enzyme_frac_poly_p_mavg    !< moving average of enzyme_frac_poly_p [unitless]

    !------------------------------------------------------------------------------------------------------ !
    !>2D total fluxes for mass balance checks in soil_biogeochemistry model
    !>
    TYPE(t_jsb_var_real2d)  :: &
      & total_flux_carbon, &       !< carbon mass balance checks in soil_biogeochemistry model [mol C m-2 s-1]
      & total_flux_nitrogen, &     !< nitrogen mass balance checks in soil_biogeochemistry model [mol N m-2 s-1]
      & total_flux_phosphorus, &   !< phosphorus mass balance checks in soil_biogeochemistry model [mol P m-2 s-1]
      & total_flux_carbon13, &     !< carbon13 mass balance checks in soil_biogeochemistry model [mol C13 m-2 s-1]
      & total_flux_carbon14, &     !< carbon14 mass balance checks in soil_biogeochemistry model [mol C14 m-2 s-1]
      & total_flux_nitrogen15      !< nitrogen15 mass balance checks in soil_biogeochemistry model [mol N15 m-2 s-1]

    !------------------------------------------------------------------------------------------------------ !
    !> biosphere-level diagnostics
    !>
    TYPE(t_jsb_var_real3d)  :: &
      & total_soil_n, &            !< sb_pool_total_bg_soil_n + nh4_solute + nh4_assoc + no3_solute [mol m-3]
      & total_soil_inorg_n         !< nh4_solute + nh4_assoc + no3_solute [mol m-3]

    TYPE(t_jsb_var_real2d)  :: &
      & ecosystem_total_n_loss     !< emission_noy + emission_n2o + emission_n2 + leaching_dom_nitrogen + leaching_nh4_solute + leaching_no3_solute
                                   !<  + lateral_loss_dom_nitrogen_sl + lateral_loss_nh4_solute_sl + lateral_loss_no3_solute_sl [micro-mol m-2 s-1]

  CONTAINS
    PROCEDURE :: Init => Init_sb_memory
  END TYPE t_sb_memory

  ! ======================================================================================================= !
  !>
  !> Type definition for sb_bgc_material_elements
  !>
  TYPE, EXTENDS(t_jsb_pool) :: t_sb_bgcm_with_elements
    TYPE(t_jsb_var_real3d), POINTER :: carbon
    TYPE(t_jsb_var_real3d), POINTER :: nitrogen
    TYPE(t_jsb_var_real3d), POINTER :: phosphorus
    TYPE(t_jsb_var_real3d), POINTER :: carbon13
    TYPE(t_jsb_var_real3d), POINTER :: carbon14
    TYPE(t_jsb_var_real3d), POINTER :: nitrogen15
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => sb_bgcm_with_elements_get_element_name_by_id
    PROCEDURE :: Init                       => sb_bgcm_with_elements_init
  END TYPE t_sb_bgcm_with_elements

  ! ======================================================================================================= !
  !>
  !> Type definition for sb_bgc_material_components
  !>
  TYPE, EXTENDS(t_jsb_pool) :: t_sb_bgcm_with_components
    TYPE(t_sb_bgcm_with_elements), POINTER :: dom              !< elements of dissolved organic matter
    TYPE(t_sb_bgcm_with_elements), POINTER :: dom_assoc        !< elements of minerally associated dissolved organic matter
    TYPE(t_sb_bgcm_with_elements), POINTER :: soluable_litter  !< elements of metabolic litter (sugars, cellulose etc.)
    TYPE(t_sb_bgcm_with_elements), POINTER :: polymeric_litter !< elements of structural litter (lignified litter)
    TYPE(t_sb_bgcm_with_elements), POINTER :: woody_litter     !< elements of physically protected wood litter
    TYPE(t_sb_bgcm_with_elements), POINTER :: fungi            !< elements of saprophytic fungal biomass
    TYPE(t_sb_bgcm_with_elements), POINTER :: mycorrhiza       !< elements of mycorrhizal fungal biomass
    TYPE(t_sb_bgcm_with_elements), POINTER :: microbial        !< elements of microbial biomass (bacteria!)
    TYPE(t_sb_bgcm_with_elements), POINTER :: residue          !< elements of necromass of fungi, microbial and mycorrhiza
    TYPE(t_sb_bgcm_with_elements), POINTER :: residue_assoc    !< elements of minerally associated microbial necromass
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => sb_bgcm_with_components_get_element_name_by_id
    PROCEDURE :: Init                       => sb_bgcm_with_components_init
  END TYPE t_sb_bgcm_with_components

  ! ======================================================================================================= !
  !>
  !> Type definition for sb_bgc_material_main
  !> the main bgc_material structure containing all SB_ bgc_material structures
  !> main structure, added directly into the "TYPE(t_jsb_memory) :: t_sb_memory" as "CLASS(t_jsb_pool) :: bgc_material"
  !>
  TYPE, EXTENDS(t_jsb_pool) :: t_sb_bgcm_main
    !! pools
    TYPE(t_sb_bgcm_with_components), POINTER :: sbbpool
    TYPE(t_sb_bgcm_with_elements),   POINTER :: sbbdelta_pools
    ! !! fluxes
    TYPE(t_sb_bgcm_with_components), POINTER :: sbbformation
    TYPE(t_sb_bgcm_with_components), POINTER :: sbbloss
    TYPE(t_sb_bgcm_with_components), POINTER :: sbbtransport
    TYPE(t_sb_bgcm_with_elements),   POINTER :: sbbmyco_export
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => sb_bgcm_main_get_element_name_by_id
    PROCEDURE :: Init                       => sb_bgcm_main_init
  END TYPE t_sb_bgcm_main

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the SB_ process
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,         ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,     ONLY: t_jsb_model
    USE mo_jsb_class,           ONLY: Get_model
    USE mo_quincy_output_class, ONLY: unitless
    dsl4jsb_Use_config(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_sb_memory),   INTENT(inout), TARGET :: mem             !> SB_ memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),               POINTER :: model                       !< model
    TYPE(t_jsb_grid),                POINTER :: hgrid                       !< Horizontal grid
    TYPE(t_jsb_vgrid),               POINTER :: surface                     !< Vertical grid
    TYPE(t_jsb_vgrid),               POINTER :: vgrid_soil_sb               !< Vertical grid
    CHARACTER(len=30)                        :: unit_sb_pool                !< unit of element var
    CHARACTER(len=30)                        :: unit_sb_flux                !< unit of element var
    INTEGER                                  :: elem_idx_map(LAST_ELEM_ID)  !< element mapper ID -> IND
    INTEGER                                  :: table                       !< ...
    CHARACTER(len=*), PARAMETER              :: routine = TRIM(modname)//':Init_sb_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    table         =  tables(1)
    hgrid         => Get_grid(model%grid_id)
    surface       => Get_vgrid('surface')
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! init local variables
    unit_sb_pool = "mol m-3"
    unit_sb_flux = "mol m-3 timestep-1"

    ! ----------------------------------------------------------------------------------------------------- !
    ! add memory only for LAND & PFT tiles
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      ! --------------------------------------------------------------------------------------------------- !
      !>  1.0 bgc_material structures
      !>
      !>  @TODO replace "b" in bgcm object names by "_" (sbbpool -> sb_pool) !!
      !>
      ! copy elements_index_map from config to local
      elem_idx_map(:) = model%config%elements_index_map(:)
      ! --------------------------------------------------------------------------------------------------- !
      !> create and populate the main bgc_material structure mem%bgc_material that contains all other bgc_material
      !> note: no spaces in 'name' and/or 'shortname'
      !>
      ALLOCATE(t_sb_bgcm_main :: mem%bgc_material)
      CALL mem%bgc_material%Init(name = 'soil_biogeochemistry_bgc_material_main', shortname = 'sb_bgcm_main', &
        &                        elements_index_map = elem_idx_map(:), l_elements = .FALSE.)
      SELECT TYPE (selected => mem%bgc_material)
      CLASS IS (t_sb_bgcm_main)
        mem%sb_bgcm_main => selected
      END SELECT

      ! --------------------------------------------------------------------------------------------------- !
      !> Add the info to this memory that this memory contains bgc materials
      mem%has_bgc_materials = .TRUE.

      ! --------------------------------------------------------------------------------------------------- !
      !> create the bgc_material components structures
      !>
      ! "shortname" is used for model output names
      ! sb_pool
      mem%sb_bgcm_main%sbbpool => Add_sb_bgc_material_with_components(mem%sb_bgcm_main,         &
        & SB_BGCM_POOL_ID, name = 'soil_biogeochemistry_bgcm_pool', shortname = 'sb_bgcm_pool', &
        & unit = unit_sb_pool, elements_index_map = elem_idx_map(:))
      ! sb_formation
      mem%sb_bgcm_main%sbbformation => Add_sb_bgc_material_with_components(mem%sb_bgcm_main,                   &
        & SB_BGCM_FORMATION_ID, name = 'soil_biogeochemistry_bgcm_formation', shortname = 'sb_bgcm_formation', &
        & unit = unit_sb_flux, elements_index_map = elem_idx_map(:))
      ! sb_loss
      mem%sb_bgcm_main%sbbloss => Add_sb_bgc_material_with_components(mem%sb_bgcm_main,         &
        & SB_BGCM_LOSS_ID, name = 'soil_biogeochemistry_bgcm_loss', shortname = 'sb_bgcm_loss', &
        & unit = unit_sb_flux, elements_index_map = elem_idx_map(:))
      ! sb_transport
      mem%sb_bgcm_main%sbbtransport => Add_sb_bgc_material_with_components(mem%sb_bgcm_main,         &
        & SB_BGCM_TRANSPORT_ID, name = 'soil_biogeochemistry_bgcm_transport', shortname = 'sb_bgcm_transport', &
        & unit = unit_sb_flux, elements_index_map = elem_idx_map(:))

      ! --------------------------------------------------------------------------------------------------- !
      !> create the bgc_material elements structures
      !>
      ! "shortname" is used for model output names
      ! sb_delta_pools
      mem%sb_bgcm_main%sbbdelta_pools => Add_sb_bgc_material_with_elements(mem%sb_bgcm_main,            &
        & name = 'soil_biogeochemistry_bgcm_delta_pools', shortname = 'sb_bgcm_delta_pools',            &
        & elements_index_map = elem_idx_map(:), unit = unit_sb_pool, bgcm_id = SB_BGCM_DELTA_POOLS_ID)
      ! sb_myco_export
      mem%sb_bgcm_main%sbbmyco_export => Add_sb_bgc_material_with_elements(mem%sb_bgcm_main,        &
        & name = 'soil_biogeochemistry_bgcm_myco_export', shortname = 'sb_bgcm_myc_exp',            &
        & elements_index_map = elem_idx_map(:), unit = unit_sb_flux, bgcm_id = SB_BGCM_MYCO_EXPORT_ID)

      ! create list of CHARACTERs with names sb_bgcm_main to the "lowest childs" of bgcm structures
      CALL mem%bgc_material%Set_paths()

      ! add variables to bgc_material
      CALL mem%Add_var(mem%sb_bgcm_main%sbbpool, &
        & hgrid, surface, vgrid_soil_sb,                                                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .TRUE.,                                                                           &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%sb_bgcm_main%sbbformation, &
        & hgrid, surface, vgrid_soil_sb,                                                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%sb_bgcm_main%sbbloss, &
        & hgrid, surface, vgrid_soil_sb,                                                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%sb_bgcm_main%sbbtransport, &
        & hgrid, surface, vgrid_soil_sb,                                                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%sb_bgcm_main%sbbdelta_pools, &
        & hgrid, surface, vgrid_soil_sb,                                                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)
      CALL mem%Add_var(mem%sb_bgcm_main%sbbmyco_export, &
        & hgrid, surface, vgrid_soil_sb,                                                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = FULL,                                                                       &
        & loutput = .TRUE.,                                                                            &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)

      ! --------------------------------------------------------------------------------------------------- !
      !> 2.0 bgc_material diagnostics variables - used for model output with IQ
      !>
      CALL mem%Add_var('sb_pool_total_c', mem%sb_pool_total_c, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_c', 'mol m-3', 'total C, sum of all sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_n', mem%sb_pool_total_n, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_n', 'mol m-3', 'total N, sum of all sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_p', mem%sb_pool_total_p, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_p', 'mol m-3', 'total P, sum of all sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_c13', mem%sb_pool_total_c13, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_c13', 'mol m-3', 'total C13, sum of all sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_c14', mem%sb_pool_total_c14, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_c14', 'mol m-3', 'total C14, sum of all sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_n15', mem%sb_pool_total_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_n15', 'mol m-3', 'total N15, sum of all sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_bg_soil_c', mem%sb_pool_total_bg_soil_c, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_bg_soil_c', 'mol m-3', 'total C, sum of SOM and below ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_bg_soil_n', mem%sb_pool_total_bg_soil_n, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_bg_soil_n', 'mol m-3', 'total N, sum of SOM and below ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_bg_soil_p', mem%sb_pool_total_bg_soil_p, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_bg_soil_p', 'mol m-3', 'total P, sum of SOM and below ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_bg_soil_c13', mem%sb_pool_total_bg_soil_c13, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_bg_soil_c13', 'mol m-3', 'total C13, sum of SOM and below ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_bg_soil_c14', mem%sb_pool_total_bg_soil_c14, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_bg_soil_c14', 'mol m-3', 'total C14, sum of SOM and below ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_bg_soil_n15', mem%sb_pool_total_bg_soil_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_total_bg_soil_n15', 'mol m-3', 'total N15, sum of SOM and below ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_ag_litter_c', mem%sb_pool_total_ag_litter_c, &
        & hgrid, surface, &
        & t_cf('sb_pool_total_ag_litter_c', 'mol m-2', 'total C, sum of the above ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_ag_litter_n', mem%sb_pool_total_ag_litter_n, &
        & hgrid, surface, &
        & t_cf('sb_pool_total_ag_litter_n', 'mol m-2', 'total N, sum of the above ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_ag_litter_p', mem%sb_pool_total_ag_litter_p, &
        & hgrid, surface, &
        & t_cf('sb_pool_total_ag_litter_p', 'mol m-2', 'total P, sum of the above ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_ag_litter_c13', mem%sb_pool_total_ag_litter_c13, &
        & hgrid, surface, &
        & t_cf('sb_pool_total_ag_litter_c13', 'mol m-2', 'total C13, sum of the above ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_ag_litter_c14', mem%sb_pool_total_ag_litter_c14, &
        & hgrid, surface, &
        & t_cf('sb_pool_total_ag_litter_c14', 'mol m-2', 'total C14, sum of the above ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_total_ag_litter_n15', mem%sb_pool_total_ag_litter_n15, &
        & hgrid, surface, &
        & t_cf('sb_pool_total_ag_litter_n15', 'mol m-2', 'total N15, sum of the above ground litter in sb_pool components'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sb_pool_woody_litter_c', mem%sb_pool_woody_litter_c, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sb_pool_woody_litter_c', 'mol m-3', 'sb_pool woody litter'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      ! --------------------------------------------------------------------------------------------------- !
      !> 3.0 2D & 3D variables
      !>
      CALL mem%Add_var('nh4_solute', mem%nh4_solute, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('nh4_solute', 'mol m-3', 'ammonium in soil solution'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nh4_assoc', mem%nh4_assoc, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('nh4_assoc', 'mol m-3', 'minerally associated ammonium'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('no3_solute', mem%no3_solute, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('no3_solute', 'mol m-3', 'nitrate in soil solution'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('noy', mem%noy, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('noy', 'mol m-3', 'NO_y'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n2o', mem%n2o, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('n2o', 'mol m-3', 'N2O'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n2', mem%n2, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('n2', 'mol m-3', 'N2'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nh4_n15_solute', mem%nh4_n15_solute, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('nh4_n15_solute', 'mol m-3', 'N15-ammonium in soil solution'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nh4_n15_assoc', mem%nh4_n15_assoc, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('nh4_n15_assoc', 'mol m-3', 'minerally associated N15-ammonium'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('no3_n15_solute', mem%no3_n15_solute, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('no3_n15_solute', 'mol m-3', 'N15-nitrate in soil solution'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('noy_n15', mem%noy_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('noy_n15', 'mol m-3', 'N15-NO_y'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n2o_n15', mem%n2o_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('n2o_n15', 'mol m-3', 'N15-N2O'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('n2_n15', mem%n2_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('n2_n15', 'mol m-3', 'N15-N2'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('po4_assoc_fast', mem%po4_assoc_fast, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('po4_assoc_fast', 'mol m-3', 'fast minerally associated phosphate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('po4_assoc_slow', mem%po4_assoc_slow, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('po4_assoc_slow', 'mol m-3', 'slow minerally associated phosphate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('po4_occluded', mem%po4_occluded, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('po4_occluded', 'mol m-3', 'occluded phosphate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('po4_primary', mem%po4_primary, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('po4_primary', 'mol m-3', 'primary phosphate mineral'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('po4_solute', mem%po4_solute, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('po4_solute', 'mol m-3', 'phophate in soil solution'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_nh4_solute', mem%transport_nh4_solute, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_nh4_solute', 'mol m-3 timestep-1', 'transport of ammonium in soil solution'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_nh4_assoc', mem%transport_nh4_assoc, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_nh4_assoc', 'mol m-3 timestep-1', 'transport of minerally associated ammonium'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_no3_solute', mem%transport_no3_solute, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_no3_solute', 'mol m-3 timestep-1', 'transport of nitrate in soil solution'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_noy', mem%transport_noy, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_noy', 'mol m-3 timestep-1', 'transport of NO_y'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_n2o', mem%transport_n2o, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_n2o', 'mol m-3 timestep-1', 'transport of N2O'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_n2', mem%transport_n2, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_n2', 'mol m-3 timestep-1', 'transport of N2'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_nh4_n15_solute', mem%transport_nh4_n15_solute, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_nh4_n15_solute', 'mol m-3 timestep-1', 'transport of N15-ammonium in soil solution'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_nh4_n15_assoc', mem%transport_nh4_n15_assoc, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_nh4_n15_assoc', 'mol m-3 timestep-1', 'transport of minerally associated N15-ammonium'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_no3_n15_solute', mem%transport_no3_n15_solute, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_no3_n15_solute', 'mol m-3 timestep-1', 'transport of N15-nitrate in soil solution'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_noy_n15', mem%transport_noy_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_noy_n15', 'mol m-3 timestep-1', 'transport of N15-NO_y'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_n2o_n15', mem%transport_n2o_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_n2o_n15', 'mol m-3 timestep-1', 'transport of N15-N2O'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_n2_n15', mem%transport_n2_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_n2_n15', 'mol m-3 timestep-1', 'transport of N15-N2'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_po4_assoc_fast', mem%transport_po4_assoc_fast, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_po4_assoc_fast', 'mol m-3 timestep-1', 'transport of fast minerally associated phosphate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_po4_assoc_slow', mem%transport_po4_assoc_slow, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_po4_assoc_slow', 'mol m-3 timestep-1', 'transport of slow minerally associated phosphate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_po4_occluded', mem%transport_po4_occluded, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_po4_occluded', 'mol m-3 timestep-1', 'transport of occluded phosphate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_po4_primary', mem%transport_po4_primary, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_po4_primary', 'mol m-3 timestep-1', 'transport of primary phosphate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transport_po4_solute', mem%transport_po4_solute, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('transport_po4_solute', 'mol m-3 timestep-1', 'transport of minerally associated phosphate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_dom_carbon_sl', mem%lateral_loss_dom_carbon_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_dom_carbon_sl', 'mol m-3 timestep-1', 'lateral loss of dom carbon'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_dom_nitrogen_sl', mem%lateral_loss_dom_nitrogen_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_dom_nitrogen_sl', 'mol m-3 timestep-1', 'lateral loss of dom nitrogen'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_dom_phosphorus_sl', mem%lateral_loss_dom_phosphorus_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_dom_phosphorus_sl', 'mol m-3 timestep-1', 'lateral loss of dom phosphorus'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_dom_carbon13_sl', mem%lateral_loss_dom_carbon13_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_dom_carbon13_sl', 'mol m-3 timestep-1', 'lateral loss of dom C13'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_dom_carbon14_sl', mem%lateral_loss_dom_carbon14_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_dom_carbon14_sl', 'mol m-3 timestep-1', 'lateral loss of dom C14'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_dom_nitrogen15_sl', mem%lateral_loss_dom_nitrogen15_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_dom_nitrogen15_sl', 'mol m-3 timestep-1', 'lateral loss of dom N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_nh4_solute_sl', mem%lateral_loss_nh4_solute_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_nh4_solute_sl', 'mol m-3 timestep-1', 'lateral loss of NH4'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_no3_solute_sl', mem%lateral_loss_no3_solute_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_no3_solute_sl', 'mol m-3 timestep-1', 'lateral loss of NO3'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_po4_solute_sl', mem%lateral_loss_po4_solute_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_po4_solute_sl', 'mol m-3 timestep-1', 'lateral loss of PO4'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_nh4_n15_solute_sl', mem%lateral_loss_nh4_n15_solute_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_nh4_n15_solute_sl', 'mol m-3 timestep-1', 'lateral loss of 15NH4'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lateral_loss_no3_n15_solute_sl', mem%lateral_loss_no3_n15_solute_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('lateral_loss_no3_n15_solute_sl', 'mol m-3 timestep-1', 'lateral loss of 15NO3'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_decomposition', mem%rtm_decomposition, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_decomposition', unitless, 'rate modifier for decomposition processes'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_depolymerisation', mem%rtm_depolymerisation, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_depolymerisation', unitless, 'rate modifier for depolymerisation processes'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_mic_uptake', mem%rtm_mic_uptake, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_mic_uptake', unitless, 'rate modifier for microbial uptake'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_plant_uptake', mem%rtm_plant_uptake, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_plant_uptake', unitless, 'rate modifier for plant uptake'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_sorption', mem%rtm_sorption, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_sorption', unitless, 'rate modifier for sorption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_desorption', mem%rtm_desorption, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_desorption', unitless, 'rate modifier for desorption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_hsc', mem%rtm_hsc, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_hsc', unitless, 'rate modifier for half-saturation constants'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_nitrification', mem%rtm_nitrification, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_nitrification', unitless, 'rate modifier for nitrification'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_denitrification', mem%rtm_denitrification, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_denitrification', unitless, 'rate modifier for denitrification'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_gasdiffusion', mem%rtm_gasdiffusion, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_gasdiffusion', unitless, 'rate modifier for gasdiffusion'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rtm_asymb_bnf', mem%rtm_asymb_bnf, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rtm_asymb_bnf', unitless, 'rate modifier for nfixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rmm_decomposition', mem%rmm_decomposition, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rmm_decomposition', unitless, 'rate modifier for decomposition processes'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rmm_depolymerisation', mem%rmm_depolymerisation, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rmm_depolymerisation', unitless, 'rate modifier for depolymerisation processes'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rmm_mic_uptake', mem%rmm_mic_uptake, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rmm_mic_uptake', unitless, 'rate modifier for microbial uptake'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rmm_plant_uptake', mem%rmm_plant_uptake, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rmm_plant_uptake', unitless, 'rate modifier for plant uptake'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rmm_sorption', mem%rmm_sorption, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rmm_sorption', unitless, 'rate modifier for sorption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rmm_desorption', mem%rmm_desorption, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rmm_desorption', unitless, 'rate modifier for desorption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rmm_hsc', mem%rmm_hsc, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rmm_hsc', unitless, 'rate modifier for half-saturation constants'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rmm_nitrification', mem%rmm_nitrification, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rmm_nitrification', unitless, 'rate modifier for nitrification'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rmm_gasdiffusion', mem%rmm_gasdiffusion, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rmm_gasdiffusion', unitless, 'rate modifier for gasdiffusion'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('rmm_asymb_bnf', mem%rmm_asymb_bnf, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('rmm_asymb_bnf', unitless, 'rate modifier for nfixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('het_respiration', mem%het_respiration, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('het_respiration', 'micro-mol m-3 s-1', 'heterotrophic respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('het_respiration_c13', mem%het_respiration_c13, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('het_respiration_c13', 'micro-mol m-3 s-1', 'c13 heterotrophic respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('het_respiration_c14', mem%het_respiration_c14, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('het_respiration_c14', 'micro-mol m-3 s-1', 'c14 heterotrophic respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('myc_respiration', mem%myc_respiration, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('myc_respiration', 'micro-mol m-3 s-1', 'mycorrhizae respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('myc_respiration_c13', mem%myc_respiration_c13, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('myc_respiration_c13', 'micro-mol m-3 s-1', 'c13 mycorrhizae respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('myc_respiration_c14', mem%myc_respiration_c14, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('myc_respiration_c14', 'micro-mol m-3 s-1', 'c14 mycorrhizae respiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('net_mineralisation_nh4', mem%net_mineralisation_nh4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('net_mineralisation_nh4', 'micro-mol m-3 s-1', 'net mineralisation of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('net_mineralisation_nh4_n15', mem%net_mineralisation_nh4_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('net_mineralisation_nh4_n15', 'micro-mol m-3 s-1', 'net mineralisation of N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('net_mineralisation_no3', mem%net_mineralisation_no3, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('net_mineralisation_no3', 'micro-mol m-3 s-1', 'net mineralisation of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('net_mineralisation_no3_n15', mem%net_mineralisation_no3_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('net_mineralisation_no3_n15', 'micro-mol m-3 s-1', 'net mineralisation of N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('net_mineralisation_po4', mem%net_mineralisation_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('net_mineralisation_po4', 'micro-mol m-3 s-1', 'net mineralisation of P'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('plant_uptake_nh4_sl', mem%plant_uptake_nh4_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('plant_uptake_nh4_sl', 'micro-mol m-3 s-1', 'plant uptake of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('plant_uptake_nh4_n15_sl', mem%plant_uptake_nh4_n15_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('plant_uptake_nh4_n15_sl', 'micro-mol m-3 s-1', 'plant uptake of N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('plant_uptake_no3_sl', mem%plant_uptake_no3_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('plant_uptake_no3_sl', 'micro-mol m-3 s-1', 'plant uptake of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('plant_uptake_no3_n15_sl', mem%plant_uptake_no3_n15_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('plant_uptake_no3_n15_sl', 'micro-mol m-3 s-1', 'plant uptake of N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('plant_uptake_po4_sl', mem%plant_uptake_po4_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('plant_uptake_po4_sl', 'micro-mol m-3 s-1', 'plant uptake of P'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_uptake_nh4_sl', mem%microbial_uptake_nh4_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_uptake_nh4_sl', 'micro-mol m-3 s-1', 'microbial uptake of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_uptake_nh4_n15_sl', mem%microbial_uptake_nh4_n15_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_uptake_nh4_n15_sl', 'micro-mol m-3 s-1', 'microbial uptake of N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_uptake_no3_sl', mem%microbial_uptake_no3_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_uptake_no3_sl', 'micro-mol m-3 s-1', 'microbial uptake of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_uptake_no3_n15_sl', mem%microbial_uptake_no3_n15_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_uptake_no3_n15_sl', 'micro-mol m-3 s-1', 'microbial uptake of N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_uptake_po4_sl', mem%microbial_uptake_po4_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_uptake_po4_sl', 'micro-mol m-3 s-1', 'microbial uptake of P'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('mycorrhiza_uptake_nh4_sl', mem%mycorrhiza_uptake_nh4_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('mycorrhiza_uptake_nh4_sl', 'micro-mol m-3 s-1', 'mycorrhiza uptake of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('mycorrhiza_uptake_nh4_n15_sl', mem%mycorrhiza_uptake_nh4_n15_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('mycorrhiza_uptake_nh4_n15_sl', 'micro-mol m-3 s-1', 'mycorrhiza uptake of N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('mycorrhiza_uptake_no3_sl', mem%mycorrhiza_uptake_no3_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('mycorrhiza_uptake_no3_sl', 'micro-mol m-3 s-1', 'mycorrhiza uptake of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('mycorrhiza_uptake_no3_n15_sl', mem%mycorrhiza_uptake_no3_n15_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('mycorrhiza_uptake_no3_n15_sl', 'micro-mol m-3 s-1', 'mycorrhiza uptake of N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('mycorrhiza_uptake_norg_sl', mem%mycorrhiza_uptake_norg_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('mycorrhiza_uptake_norg_sl', 'micro-mol m-3 s-1', 'mycorrhiza uptake of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('mycorrhiza_uptake_norg_n15_sl', mem%mycorrhiza_uptake_norg_n15_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('mycorrhiza_uptake_norg_n15_sl', 'micro-mol m-3 s-1', 'mycorrhiza uptake of N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('mycorrhiza_uptake_po4_sl', mem%mycorrhiza_uptake_po4_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('mycorrhiza_uptake_po4_sl', 'micro-mol m-3 s-1', 'mycorrhiza uptake of P'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('myc_export_n_tlabile_mavg_sl', mem%myc_export_n_tlabile_mavg_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('myc_export_n_tlabile_mavg_sl', 'mol m-3', 'time averaged mycorrhizae export of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('myc_export_p_tlabile_mavg_sl', mem%myc_export_p_tlabile_mavg_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('myc_export_p_tlabile_mavg_sl', 'mol m-3', 'time averaged mycorrhizae export of P'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('myc_export_c_tmyc_mavg_sl', mem%myc_export_c_tmyc_mavg_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('myc_export_c_tmyc_mavg_sl', 'mol m-3', 'time averaged mycorrhizae export of C'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('myc_export_n_tmyc_mavg_sl', mem%myc_export_n_tmyc_mavg_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('myc_export_n_tmyc_mavg_sl', 'mol m-3', 'time averaged mycorrhizae export of N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('volatilisation_nh4', mem%volatilisation_nh4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('volatilisation_nh4', 'micro-mol m-3 s-1', 'volatilisation of NH4'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('volatilisation_nh4_n15', mem%volatilisation_nh4_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('volatilisation_nh4_n15', 'micro-mol m-3 s-1', 'volatilisation of 15NH4'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nitrification_no3', mem%nitrification_no3, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('nitrification_no3', 'micro-mol m-3 s-1', 'nitrification NO3'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nitrification_no3_n15', mem%nitrification_no3_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('nitrification_no3_n15', 'micro-mol m-3 s-1', 'nitrification 15NO3'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nitrification_noy', mem%nitrification_noy, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('nitrification_noy', 'micro-mol m-3 s-1', 'nitrification to NOy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nitrification_noy_n15', mem%nitrification_noy_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('nitrification_noy_n15', 'micro-mol m-3 s-1', 'nitrification to 15NOy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nitrification_n2o', mem%nitrification_n2o, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('nitrification_n2o', 'micro-mol m-3 s-1', 'nitrification to N2O'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nitrification_n2o_n15', mem%nitrification_n2o_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('nitrification_n2o_n15', 'micro-mol m-3 s-1', 'nitrification to 15N2O'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('denitrification_noy', mem%denitrification_noy, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('denitrification_noy', 'micro-mol m-3 s-1', 'nitrification NOy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('denitrification_n2o_n15', mem%denitrification_n2o_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('denitrification_n2o_n15', 'micro-mol m-3 s-1', 'nitrification to 15N2O'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('denitrification_n2o', mem%denitrification_n2o, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('denitrification_n2o', 'micro-mol m-3 s-1', 'nitrification to N2O'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('denitrification_noy_n15', mem%denitrification_noy_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('denitrification_noy_n15', 'micro-mol m-3 s-1', 'nitrification 15NOy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('denitrification_n2', mem%denitrification_n2, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('denitrification_n2', 'micro-mol m-3 s-1', 'denitrification to N2'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('denitrification_n2_n15', mem%denitrification_n2_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('denitrification_n2_n15', 'micro-mol m-3 s-1', 'denitrification to 15N2'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('asymb_n_fixation', mem%asymb_n_fixation, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('asymb_n_fixation', 'micro-mol m-3 s-1', 'asymboitic N fixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = MEDIUM, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('asymb_n_fixation_n15', mem%asymb_n_fixation_n15, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('asymb_n_fixation_n15', 'micro-mol m-3 s-1', 'asymboitic 15N fixation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('weathering_po4', mem%weathering_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('weathering_po4', 'micro-mol m-3 s-1', 'phosphorus weathering rate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('occlusion_po4', mem%occlusion_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('occlusion_po4', 'micro-mol m-3 s-1', 'phosphorus occlusion rate, assoc_slow --> occluded'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('slow_exchange_po4', mem%slow_exchange_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('slow_exchange_po4', 'micro-mol P m-3 s-1', 'phosphorus slow exchange rate, assoc_fast <--> assoc_slow'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fast_exchange_po4', mem%fast_exchange_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('fast_exchange_po4', 'micro-mol P m-3 s-1', 'phosphorus fast exchange rate, solute <--> assoc_fast'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fast_exchange_nh4', mem%fast_exchange_nh4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('fast_exchange_nh4', 'micro-mol N m-3 s-1', 'nh4 fast exchange rate, solute <--> nh4_assoc'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fast_adsorpt_po4', mem%fast_adsorpt_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('fast_adsorpt_po4', 'micro-mol P m-3 s-1', 'phosphorus fast exchange rate, solute <--> assoc_fast'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fast_adsorpt_nh4', mem%fast_adsorpt_nh4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('fast_adsorpt_nh4', 'micro-mol N m-3 s-1', 'nh4 fast adsorption rate, solute --> nh4_assoc'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('biochem_mineralisation_po4', mem%biochem_mineralisation_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('biochem_mineralisation_po4', 'micro-mol m-3 s-1', 'phosphorus biochemical mineralization rate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('gross_mineralisation_po4', mem%gross_mineralisation_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('gross_mineralisation_po4', 'micro-mol m-3 s-1', 'phosphorus gross (biological) mineralization rate'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_dom_carbon', mem%leaching_dom_carbon, &
        & hgrid, surface, &
        & t_cf('leaching_dom_carbon', 'micro-mol m-2 s-1', 'land to groundwater flux of DOM-C'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_dom_nitrogen', mem%leaching_dom_nitrogen, &
        & hgrid, surface, &
        & t_cf('leaching_dom_nitrogen', 'micro-mol m-2 s-1', 'land to groundwater flux of DOM-N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_dom_phosphorus', mem%leaching_dom_phosphorus, &
        & hgrid, surface, &
        & t_cf('leaching_dom_phosphorus', 'micro-mol m-2 s-1', 'land to groundwater flux of DOM-P'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_dom_carbon13', mem%leaching_dom_carbon13, &
        & hgrid, surface, &
        & t_cf('leaching_dom_carbon13', 'micro-mol m-2 s-1', 'land to groundwater flux of DOM-C13'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_dom_carbon14', mem%leaching_dom_carbon14, &
        & hgrid, surface, &
        & t_cf('leaching_dom_carbon14', 'micro-mol m-2 s-1', 'land to groundwater flux of DOM-C14'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_dom_nitrogen15', mem%leaching_dom_nitrogen15, &
        & hgrid, surface, &
        & t_cf('leaching_dom_nitrogen15', 'micro-mol m-2 s-1', 'land to groundwater flux of DOM-N15'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_nh4_solute', mem%leaching_nh4_solute, &
        & hgrid, surface, &
        & t_cf('leaching_nh4_solute', 'micro-mol m-2 s-1', 'land to groundwater flux of NH4'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_no3_solute', mem%leaching_no3_solute, &
        & hgrid, surface, &
        & t_cf('leaching_no3_solute', 'micro-mol m-2 s-1', 'land to groundwater flux of NO3'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_po4_solute', mem%leaching_po4_solute, &
        & hgrid, surface, &
        & t_cf('leaching_po4_solute', 'micro-mol m-2 s-1', 'land to groundwater flux of PO4'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_nh4_n15_solute', mem%leaching_nh4_n15_solute, &
        & hgrid, surface, &
        & t_cf('leaching_nh4_n15_solute', 'micro-mol m-2 s-1', 'land to groundwater flux of 15NH4'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('leaching_no3_n15_solute', mem%leaching_no3_n15_solute, &
        & hgrid, surface, &
        & t_cf('leaching_no3_n15_solute', 'micro-mol m-2 s-1', 'land to groundwater flux of 15NO3'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('emission_noy', mem%emission_noy, &
        & hgrid, surface, &
        & t_cf('emission_noy', 'micro-mol m-2 s-1', 'land to atmosphere flux of NOy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('emission_n2o', mem%emission_n2o, &
        & hgrid, surface, &
        & t_cf('emission_n2o', 'micro-mol m-2 s-1', 'land to atmosphere flux of N2O'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('emission_n2', mem%emission_n2, &
        & hgrid, surface, &
        & t_cf('emission_n2', 'micro-mol m-2 s-1', 'land to atmosphere flux of N2'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('emission_noy_n15', mem%emission_noy_n15, &
        & hgrid, surface, &
        & t_cf('emission_noy_n15', 'micro-mol m-2 s-1', 'land to atmosphere flux of 15NOy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('emission_n2o_n15', mem%emission_n2o_n15, &
        & hgrid, surface, &
        & t_cf('emission_n2o_n15', 'micro-mol m-2 s-1', 'land to atmosphere flux of 15N2O'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('emission_n2_n15', mem%emission_n2_n15, &
        & hgrid, surface, &
        & t_cf('emission_n2_n15', 'micro-mol m-2 s-1', 'land to atmosphere flux of 15N2'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_org', mem%qmax_org, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_org', 'kg C kg-1 mineral soil', 'sorption capacity of dom of mineral soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_po4', mem%qmax_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_po4', 'mol P m-3 soil', 'sorption capacity of po4 of mineral soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('km_adsorpt_po4_sl', mem%km_adsorpt_po4_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('km_adsorpt_po4_sl', 'mol P m-3 soil', 'half saturation of po4 sorption isotherm'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_nh4', mem%qmax_nh4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_nh4', 'mol N m-3 soil', 'sorption capacity of nh4 of mineral soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('km_adsorpt_nh4_sl', mem%km_adsorpt_nh4_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('km_adsorpt_nh4_sl', 'mol N m-3 soil', 'half saturation of nh4 sorption isotherm'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
			
      CALL mem%Add_var('qmax_fast_po4', mem%qmax_fast_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_fast_po4', 'mol P m-3 soil', 'sorption capacity of po4 of the clay(Fe/Al) and SOM pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_slow_po4', mem%qmax_slow_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_slow_po4', 'mol P m-3 soil', 'sorption capacity of po4 of the silt/sand pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('km_fast_po4', mem%km_fast_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('km_fast_po4', 'mol P m-3 soil', 'half-saturation concentration of po4 for the clay(Fe/Al) and SOM pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('km_slow_po4', mem%km_slow_po4, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('km_slow_po4', 'mol P m-3 soil', 'half-saturation concentration of po4 for the silt/sand pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('k_bioturb', mem%k_bioturb, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('k_bioturb', 'm2 s-1', 'diffusion rate resulting from bioturbation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('bulk_dens_corr_sl', mem%bulk_dens_corr_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('bulk_dens_corr_sl', 'kg m-3', 'organic fraction corrected bulk-density of the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('bulk_soil_carbon_sl', mem%bulk_soil_carbon_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('bulk_soil_carbon_sl', 'mol m-3', 'bulk soil density of non-litter carbon'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('soil_litter_carbon_sl', mem%soil_litter_carbon_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_litter_carbon_sl', 'mol m-3', 'soil density of litter carbon'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('particle_fluxrate', mem%particle_fluxrate, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('particle_fluxrate', 'm s-1', 'burial rate due to organic matter fall and decay'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('anaerobic_volume_fraction_sl', mem%anaerobic_volume_fraction_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('anaerobic_volume_fraction_sl', unitless, 'anaerobic volume fraction'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_cue_eff', mem%microbial_cue_eff, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_cue_eff', unitless, 'instantaneous microbial CUE given NP constraints'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_cue_eff_tmic_mavg', mem%microbial_cue_eff_tmic_mavg, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_cue_eff_tmic_mavg', unitless, 'time-averaged microbial CUE given NP constraints'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_nue_eff', mem%microbial_nue_eff, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_nue_eff', unitless, 'instantaneous microbial NUE given NP constraints'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_nue_eff_tmic_mavg', mem%microbial_nue_eff_tmic_mavg, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_nue_eff_tmic_mavg', unitless, 'time-averaged microbial NUE given NP constraints'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_pue_eff', mem%microbial_pue_eff, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_pue_eff', unitless, 'instantaneous microbial PUE given NP constraints'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('microbial_pue_eff_tmic_mavg', mem%microbial_pue_eff_tmic_mavg, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('microbial_pue_eff_tmic_mavg', unitless, 'time-averaged microbial CUE given NP constraints'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fact_n_status_mic_c_growth', mem%fact_n_status_mic_c_growth, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('fact_n_status_mic_c_growth', unitless, 'microbial N status, scaling up/down microbial processes'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('fact_p_status_mic_c_growth', mem%fact_p_status_mic_c_growth, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('fact_p_status_mic_c_growth', unitless, 'microbial P status, scaling up/down microbial processes'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('enzyme_frac_AP', mem%enzyme_frac_AP, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('enzyme_frac_AP', unitless, 'fraction of Acid (alkaline) phosphatase of the total enzyme'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dom_cn', mem%dom_cn, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('dom_cn', 'mol mol-1', 'C:N ratio of DOM pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dom_cp', mem%dom_cp, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('dom_cp', 'mol mol-1', 'C:P ratio of DOM pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dom_cn_mavg', mem%dom_cn_mavg, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('dom_cn_mavg', 'mol mol-1', 'time-averaged C:N ratio of DOM pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dom_cp_mavg', mem%dom_cp_mavg, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('dom_cp_mavg', 'mol mol-1', 'time-averaged C:P ratio of DOM pool'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      IF (model%config%flag_slow_sb_pool_spinup_accelerator) THEN
        CALL mem%Add_var('residue_som_c_form_mavg_sl', mem%residue_som_c_form_mavg_sl, &
          & hgrid, vgrid_soil_sb, &
          & t_cf('residue_som_c_form_mavg_sl', 'micro-mol m-3 s-1', 'time-averaged residual som C formation rate'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = FULL, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('residue_som_c14_form_mavg_sl', mem%residue_som_c14_form_mavg_sl, &
          & hgrid, vgrid_soil_sb, &
          & t_cf('residue_som_c14_form_mavg_sl', 'micro-mol m-3 s-1', 'time-averaged residual som C14 formation rate'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = FULL, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('residue_assoc_som_c_form_mavg_sl', mem%residue_assoc_som_c_form_mavg_sl, &
          & hgrid, vgrid_soil_sb, &
          & t_cf('residue_assoc_som_c_form_mavg_sl', 'micro-mol m-3 s-1', 'time-averaged assoc residual som C formation rate'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = FULL, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('residue_assoc_som_c14_form_mavg_sl', mem%residue_assoc_som_c14_form_mavg_sl, &
          & hgrid, vgrid_soil_sb, &
          & t_cf('residue_assoc_som_c14_form_mavg_sl', 'micro-mol m-3 s-1', 'time-averaged assoc resid som C14 formation rate'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = FULL, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('assoc_dom_c_form_mavg_sl', mem%assoc_dom_c_form_mavg_sl, &
          & hgrid, vgrid_soil_sb, &
          & t_cf('assoc_dom_c_form_mavg_sl', 'micro-mol m-3 s-1', 'time-averaged associated dom C formation rate'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = FULL, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('assoc_dom_c14_form_mavg_sl', mem%assoc_dom_c14_form_mavg_sl, &
          & hgrid, vgrid_soil_sb, &
          & t_cf('assoc_dom_c14_form_mavg_sl', 'micro-mol m-3 s-1', 'time-averaged associated dom C14 formation rate'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = FULL, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('residue_som_c_loss_mavg_sl', mem%residue_som_c_loss_mavg_sl, &
          & hgrid, vgrid_soil_sb, &
          & t_cf('residue_som_c_loss_mavg_sl', 'micro-mol m-3 s-1', 'time-averaged residual som C loss rate'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = FULL, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('residue_assoc_som_c_loss_mavg_sl', mem%residue_assoc_som_c_loss_mavg_sl, &
          & hgrid, vgrid_soil_sb, &
          & t_cf('residue_assoc_som_c_loss_mavg_sl', 'micro-mol m-3 s-1', 'time-averaged associated residual som C loss rate'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = FULL, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

        CALL mem%Add_var('assoc_dom_c_loss_mavg_sl', mem%assoc_dom_c_loss_mavg_sl, &
          & hgrid, vgrid_soil_sb, &
          & t_cf('assoc_dom_c_loss_mavg_sl', 'micro-mol m-3 s-1', 'time-averaged associated dom C loss rate'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & output_level = FULL, &
          & loutput = .FALSE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)
      END IF

      CALL mem%Add_var('ph_sl', mem%ph_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('ph_sl', unitless, 'pH per soil layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('vmax_weath_mineral_sl', mem%vmax_weath_mineral_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('vmax_weath_mineral_sl', 'mol P m-3 s-1', 'weathering coefficient of mineral soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('Qmax_AlFe_cor', mem%Qmax_AlFe_cor, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('Qmax_AlFe_cor', unitless, 'correction factor for phosphate Qmax due to the unknown Al/Fe'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('enzyme_frac_poly', mem%enzyme_frac_poly, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('enzyme_frac_poly', unitless, 'fraction of depolymerisation enzyme which depolymerises poly litter'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('enzyme_frac_residue', mem%enzyme_frac_residue, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('enzyme_frac_residue', unitless, 'fraction of depolymerisation enzyme which depolymerises microbial residue'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('enzyme_frac_poly_c', mem%enzyme_frac_poly_c, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('enzyme_frac_poly_c', unitless, 'fract of depolymerisation enzyme depol. microbial residue, determined by C'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('enzyme_frac_poly_n', mem%enzyme_frac_poly_n, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('enzyme_frac_poly_n', unitless, 'fract of depolymerisation enzyme depol. microbial residue, determined by N'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('enzyme_frac_poly_p', mem%enzyme_frac_poly_p, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('enzyme_frac_poly_p', unitless, 'fract of depolymerisation enzyme depol. microbial residue, determined by P'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('enzyme_frac_poly_c_mavg', mem%enzyme_frac_poly_c_mavg, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('enzyme_frac_poly_c_mavg', unitless, 'moving average of enzyme_frac_poly_c'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('enzyme_frac_poly_n_mavg', mem%enzyme_frac_poly_n_mavg, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('enzyme_frac_poly_n_mavg', unitless, 'moving average of enzyme_frac_poly_n'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('enzyme_frac_poly_p_mavg', mem%enzyme_frac_poly_p_mavg, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('enzyme_frac_poly_p_mavg', unitless, 'moving average of enzyme_frac_poly_p'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('total_flux_carbon', mem%total_flux_carbon, &
        & hgrid, surface, &
        & t_cf('total_flux_carbon', 'mol C m-2 s-1', 'carbon mass balance checks in soil_biogeochemistry model'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('total_flux_nitrogen', mem%total_flux_nitrogen, &
        & hgrid, surface, &
        & t_cf('total_flux_nitrogen', 'mol N m-2 s-1', 'nitrogen mass balance checks in soil_biogeochemistry model'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('total_flux_phosphorus', mem%total_flux_phosphorus, &
        & hgrid, surface, &
        & t_cf('total_flux_phosphorus', 'mol P m-2 s-1', 'phosphorus mass balance checks in soil_biogeochemistry model'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('total_flux_carbon13', mem%total_flux_carbon13, &
        & hgrid, surface, &
        & t_cf('total_flux_carbon13', 'mol C13 m-2 s-1', 'carbon13 mass balance checks in soil_biogeochemistry model'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('total_flux_carbon14', mem%total_flux_carbon14, &
        & hgrid, surface, &
        & t_cf('total_flux_carbon14', 'mol C14 m-2 s-1', 'carbon14 mass balance checks in soil_biogeochemistry model'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('total_flux_nitrogen15', mem%total_flux_nitrogen15, &
        & hgrid, surface, &
        & t_cf('total_flux_nitrogen15', 'mol N15 m-2 s-1', 'nitrogen15 mass balance checks in soil_biogeochemistry model'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('total_soil_n', mem%total_soil_n, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('total_soil_n', 'mol m-3', 'sb_pool_total_bg_soil_n + nh4_solute + nh4_assoc + no3_solute'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('total_soil_inorg_n', mem%total_soil_inorg_n, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('total_soil_inorg_n', 'mol m-3', 'nh4_solute + nh4_assoc + no3_solute'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ecosystem_total_n_loss', mem%ecosystem_total_n_loss, &
        & hgrid, surface, &
        & t_cf('ecosystem_total_n_loss', 'micro-mol m-2 s-1', 'emissions + leaching + lateral loss'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

    ENDIF  ! IF One_of(LAND_TYPE  OR VEG_TYPE, lct_ids)
  END SUBROUTINE Init_sb_memory

  ! ======================================================================================================= !
  !>
  !> return name of element variable of bgc_material by its ENUM ID
  !> helper function for t_sb_bgcm_with_elements
  !> "bgcm_object%Get_element_name_by_id" => sb_bgcm_with_elements_get_element_name_by_id
  !>
  FUNCTION sb_bgcm_with_elements_get_element_name_by_id(this, id) RESULT(name)
    CLASS(t_sb_bgcm_with_elements),     INTENT(in)  :: this
    INTEGER,                            INTENT(in)  :: id
    CHARACTER(LEN=:),                   ALLOCATABLE :: name

    ! return name of element variable only if this bgc_material object contains elements
    IF (this%contains_elements) THEN
      name = get_name_for_element_id(id)
    ELSE
      name='no_elements'
    END IF
  END FUNCTION sb_bgcm_with_elements_get_element_name_by_id

  ! ======================================================================================================= !
  !>
  !> return name of element variable of bgc_material by its ENUM ID
  !> helper function for t_sb_bgcm_with_components
  !> "bgcm_object%Get_element_name_by_id" => sb_bgcm_with_components_get_element_name_by_id
  !>
  FUNCTION sb_bgcm_with_components_get_element_name_by_id(this, id) RESULT(name)
    CLASS(t_sb_bgcm_with_components),     INTENT(in)  :: this
    INTEGER,                              INTENT(in)  :: id
    CHARACTER(LEN=:),                     ALLOCATABLE :: name

    ! return name of element variable only if this bgc_material object contains elements
    IF (this%contains_elements) THEN
      name = get_name_for_element_id(id)
    ELSE
      name='no_elements'
    END IF
  END FUNCTION sb_bgcm_with_components_get_element_name_by_id

  ! ======================================================================================================= !
  !>
  !> return name of element variable of bgc_material by its ENUM ID
  !> helper function for t_sb_bgcm_with_components
  !> "bgcm_object%Get_element_name_by_id" => sb_bgcm_main_get_element_name_by_id
  !>
  FUNCTION sb_bgcm_main_get_element_name_by_id(this, id) RESULT(name)
    CLASS(t_sb_bgcm_main),          INTENT(in)  :: this
    INTEGER,                        INTENT(in)  :: id
    CHARACTER(LEN=:),               ALLOCATABLE :: name

    ! return name of element variable only if this bgc_material object contains elements
    IF (this%contains_elements) THEN
      name = get_name_for_element_id(id)
    ELSE
      name='no_elements'
    END IF
  END FUNCTION sb_bgcm_main_get_element_name_by_id

  ! ======================================================================================================= !
  !>
  !> initializing the t_sb_bgcm_with_elements object
  !> "bgcm_object%Init()"
  !> by default it does contain element variables !
  !>
  SUBROUTINE sb_bgcm_with_elements_init(this, name, shortname, elements_index_map, &
    &                                   l_elements, element_list, element_unit)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_sb_bgcm_with_elements),     INTENT(inout) :: this                   ! bgcm object
    CHARACTER(LEN=*),                   INTENT(in)    :: name                   ! pool name       (No spaces !)
    CHARACTER(LEN=*),                   INTENT(in)    :: shortname              ! pool short name (No spaces !)
    INTEGER,                            INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
    LOGICAL,          OPTIONAL,         INTENT(in)    :: l_elements             ! IF pool has elements
    INTEGER,          OPTIONAL,         INTENT(in)    :: element_list(:) ! @TODO remove | list of elements; used if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL,         INTENT(in)    :: element_unit           ! unit of elements
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: element_id
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_bgcm_with_elements_init'
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
        CALL this%Add_element(get_name_for_element_id(element_id), get_short_name_for_element_id(element_id), element_id, dim=3)

        SELECT TYPE (selector => this%element_list(elements_index_map(element_id))%p)
        CLASS IS (t_jsb_var_real3d)
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

  END SUBROUTINE sb_bgcm_with_elements_init

  ! ======================================================================================================= !
  !>
  !> initializing the t_sb_bgcm_with_components object
  !> "bgcm_object%Init()"
  !> by default it does not contain element variables (only t_sb_bgcm_with_elements do)
  !>
  SUBROUTINE sb_bgcm_with_components_init(this, name, shortname, elements_index_map, &
    &                                     l_elements, element_list, element_unit)
    CLASS(t_sb_bgcm_with_components),     INTENT(inout) :: this                   ! bgcm object
    CHARACTER(LEN=*),                     INTENT(in)    :: name                   ! pool name       (No spaces !)
    CHARACTER(LEN=*),                     INTENT(in)    :: shortname              ! pool short name (No spaces !)
    INTEGER,                              INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
    LOGICAL,          OPTIONAL,           INTENT(in)    :: l_elements             ! IF pool has elements
    INTEGER,          OPTIONAL,           INTENT(in)    :: element_list(:)        ! list of elements; used if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL,           INTENT(in)    :: element_unit           ! unit of elements
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_bgcm_with_components_init'

    this%name              = name
    this%shortname         = shortname
    this%contains_elements = l_elements
    this%contains_bgcm     = .TRUE.

    IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')
    this%element_unit = element_unit
  END SUBROUTINE sb_bgcm_with_components_init

  ! ======================================================================================================= !
  !>
  !> initializing the t_sb_bgcm_main object
  !> "bgcm_object%Init()"
  !> by default it does not contain element variables (only t_sb_bgcm_with_elements do)
  !>
  SUBROUTINE sb_bgcm_main_init(this, name, shortname, elements_index_map, &
    &                          l_elements, element_list, element_unit)
    CLASS(t_sb_bgcm_main),          INTENT(inout) :: this                   ! bgcm object
    CHARACTER(LEN=*),               INTENT(in)    :: name                   ! pool name       (No spaces !)
    CHARACTER(LEN=*),               INTENT(in)    :: shortname              ! pool short name (No spaces !)
    INTEGER,                        INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
    LOGICAL,          OPTIONAL,     INTENT(in)    :: l_elements             ! IF pool has elements
    INTEGER,          OPTIONAL,     INTENT(in)    :: element_list(:)        ! list of elements; used if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL,     INTENT(in)    :: element_unit           ! unit of elements
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_bgcm_main_init'

    this%name              = name
    this%shortname         = shortname
    this%contains_elements = l_elements
    this%contains_bgcm     = .TRUE.
  END SUBROUTINE sb_bgcm_main_init

  ! ======================================================================================================= !
  !>
  !> create, init and append a sb bgcm with components
  !>
  FUNCTION Add_sb_bgc_material_with_components(pool, bgcm_id, name, shortname, unit, elements_index_map) &
      & RESULT(sb_bgcm_with_components)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_pool), TARGET,  INTENT(INOUT) :: pool                  !< pool to which this pool with component should be added
    INTEGER,                    INTENT(in)    :: bgcm_id               !< id of this sb bgcm with components
    CHARACTER(LEN=*),           INTENT(in)    :: name                  !< name of this sb bgcm with components
    CHARACTER(LEN=*),           INTENT(in)    :: shortname             !< shortname of this sb bgcm with components
    CHARACTER(LEN=*),           INTENT(in)    :: unit                  !< unit of this sb bgcm with components
    INTEGER,                    INTENT(in)    :: elements_index_map(:) !< elements ID -> IND
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_sb_bgcm_with_components), POINTER  :: sb_bgcm_with_components !< the sb bgcm with components
    ! ----------------------------------------------------------------------------------------------------- !
    ! bgcm_components init
    ALLOCATE(sb_bgcm_with_components)
    CALL sb_bgcm_with_components%Init(name, shortname, elements_index_map(:), l_elements = .FALSE., element_unit = unit)

    ! bgcm_elements allocate & init & add to bgcm_components
    !   "shortname" is used for model output
    !   "element_list" is optional and can be used to select specific elements (not needed when all elements may be added)

    ! dom
    sb_bgcm_with_components%dom => Add_sb_bgc_material_with_elements(sb_bgcm_with_components, &
      & name = 'dom_bgcm', shortname = 'dom', elements_index_map = elements_index_map(:), unit = unit)
    ! dom_assoc
    sb_bgcm_with_components%dom_assoc => Add_sb_bgc_material_with_elements(sb_bgcm_with_components, &
      & name = 'dom_assoc_bgcm', shortname = 'dom_assoc', elements_index_map = elements_index_map(:), unit = unit)
    ! soluable_litter
    sb_bgcm_with_components%soluable_litter => Add_sb_bgc_material_with_elements(sb_bgcm_with_components, &
      & name = 'soluable_litter_bgcm', shortname = 'sol_litter', elements_index_map = elements_index_map(:), unit = unit)
    ! polymeric_litter
    sb_bgcm_with_components%polymeric_litter => Add_sb_bgc_material_with_elements(sb_bgcm_with_components, &
      & name = 'polymeric_litter_bgcm', shortname = 'pol_litter', elements_index_map = elements_index_map(:), unit = unit)
    ! woody_litter
    sb_bgcm_with_components%woody_litter => Add_sb_bgc_material_with_elements(sb_bgcm_with_components, &
      & name = 'woody_litter_bgcm', shortname = 'woo_litter', elements_index_map = elements_index_map(:), unit = unit)
    ! fungi
    sb_bgcm_with_components%fungi => Add_sb_bgc_material_with_elements(sb_bgcm_with_components, &
      & name = 'fungi_bgcm', shortname = 'fungi', elements_index_map = elements_index_map(:), unit = unit)
    ! mycorrhiza
    sb_bgcm_with_components%mycorrhiza => Add_sb_bgc_material_with_elements(sb_bgcm_with_components, &
      & name = 'mycorrhiza_bgcm', shortname = 'mycorrhiza', elements_index_map = elements_index_map(:), unit = unit)
    ! microbial
    sb_bgcm_with_components%microbial => Add_sb_bgc_material_with_elements(sb_bgcm_with_components, &
      & name = 'microbial_bgcm', shortname = 'microbial', elements_index_map = elements_index_map(:), unit = unit)
    ! residue
    sb_bgcm_with_components%residue => Add_sb_bgc_material_with_elements(sb_bgcm_with_components, &
      & name = 'residue_bgcm', shortname = 'residue', elements_index_map = elements_index_map(:), unit = unit)
    ! residue_assoc
    sb_bgcm_with_components%residue_assoc => Add_sb_bgc_material_with_elements(sb_bgcm_with_components, &
      & name = 'residue_assoc_bgcm', shortname = 'residue_assoc', elements_index_map = elements_index_map(:), unit = unit)

    ! Add the bookkeepoing for this sb bgcm with components
    CALL pool%Add_pool(sb_bgcm_with_components, bgcm_id)
  END FUNCTION Add_sb_bgc_material_with_components

  ! ======================================================================================================= !
  !>
  !> create, init and append a sb bgcm with elements
  !>
  FUNCTION Add_sb_bgc_material_with_elements(pool,  name, shortname, elements_index_map, unit, bgcm_id) &
      & RESULT(sb_bgcm_with_elements)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_pool), TARGET,  INTENT(INOUT) :: pool                  !< pool to which this pool with elements should be added
    CHARACTER(LEN=*),           INTENT(in)    :: name                  !< name of this sb bgcm with elements
    CHARACTER(LEN=*),           INTENT(in)    :: shortname             !< shortname of this sb bgcm with elements
    INTEGER,                    INTENT(in)    :: elements_index_map(:) !< mapping elements ID -> IND
    CHARACTER(LEN=*),           INTENT(in)    :: unit                  !< unit of this sb bgcm with elements
    INTEGER,          OPTIONAL, INTENT(in)    :: bgcm_id               !< optional: id of this bgcm with elements
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_sb_bgcm_with_elements), POINTER   :: sb_bgcm_with_elements !< the new sb bgcm with elements
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(sb_bgcm_with_elements)

    CALL sb_bgcm_with_elements%Init(name = name, shortname = shortname,  &
      & elements_index_map = elements_index_map(:), l_elements = .TRUE., element_unit = unit)
    CALL pool%Add_pool(sb_bgcm_with_elements, bgcm_id)
  END FUNCTION Add_sb_bgc_material_with_elements

#endif
END MODULE mo_sb_memory_class
