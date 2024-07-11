!> QUINCY 'jena soil model' calculate soil-biogeochemical pools
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
!>#### main routines for JSM (jena soil model), calculating the soil biogeochemical pools
!>
MODULE mo_q_sb_jsm_main
#ifndef __NO_QUINCY__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message, message_text, finish

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID, SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, SB_BGCM_TRANSPORT_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_jsm

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sb_jsm_main'

CONTAINS

  ! ======================================================================================================= !
  !>update soil biogeochemical pools - JSM - Jena Soil Model
  !>
  SUBROUTINE update_sb_jsm(tile, options)
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: SB_, SPQ_, VEG_, A2L_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_jsb_math_constants,    ONLY: dtime, eps8
    USE mo_isotope_util,          ONLY: calc_fractionation
    USE mo_sb_constants
    USE mo_q_sb_litter_processes, ONLY: calc_litter_partitioning
    USE mo_q_sb_jsm_processes,    ONLY: calc_sinking_flux, calc_litter_turnover, &
      &                                 calc_microbial_growth, calc_sourcing_flux, &
      &                                 calc_microbial_turnover, calc_sorption_desorption_of_org_material, &
      &                                 calc_depolymerisation_SESAM, &
      &                                 calc_nitri_denitri_partitioning_rate, calc_asymb_bnf_fraction
    USE mo_q_sb_jsm_transport,    ONLY: calc_liquid_phase_transport_wrapper, &
      &                                 calc_bioturbation_transport_wrapper_jsm, &
      &                                 calc_bioturbation_rate
    USE mo_q_sb_ssm_main,         ONLY: calc_sb_inorganic_n15_fluxes
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model                !< the model
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store           !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_soil_sb        !< Vertical grid
    INTEGER                               :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    INTEGER                               :: ic                   !< looping over points
    INTEGER                               :: isoil                !< looping over soil layers
    REAL(wp)                              :: hlp1                 !< helper
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: arr_hlp1             !< helper, Dimension: (nc, nsoil_sb)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: arr_hlp2             !< helper, Dimension: (nc, nsoil_sb)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: arr_hlp3             !< helper, Dimension: (nc, nsoil_sb)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: fine_root_carbon_sl                      !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: km1_uptake_act                           !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: km2_uptake_act                           !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: km_array                                 !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: avail_soil_surface_sorption_capacity_sl  !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: graze_mineralisation_po4                 !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: graze_mineralisation_nh4                 !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: graze_mineralisation_nh4_n15             !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: unit_uptake_n_pot_sl                     !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: unit_uptake_p_pot_sl                     !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: growth_p_limit_smoothed_sl               !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: asymb_n_fixation_rel_rate_sl             !< asymbiotic N fixation rate (relativ to N demand)
    INTEGER                               :: iblk, ics, ice, nc             !< dimensions
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_jsm'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_formation_mt
    dsl4jsb_Def_mt2L3D :: sb_loss_mt
    dsl4jsb_Def_mt2L3D :: sb_transport_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: percolation_sl
    dsl4jsb_Real3D_onChunk :: frac_w_lat_loss_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl
    dsl4jsb_Real3D_onChunk :: volume_min_sl
    ! A2L_ 2D
    dsl4jsb_Real2D_onChunk :: p_deposition
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk :: growth_p_limit_based_on_n_mavg
    dsl4jsb_Real2D_onChunk :: unit_uptake_n_pot
    dsl4jsb_Real2D_onChunk :: unit_uptake_p_pot
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon
    dsl4jsb_Real2D_onChunk :: leaching_dom_nitrogen
    dsl4jsb_Real2D_onChunk :: leaching_dom_phosphorus
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon13
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon14
    dsl4jsb_Real2D_onChunk :: leaching_dom_nitrogen15
    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    dsl4jsb_Real2D_onChunk :: leaching_nh4_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_n15_solute
    dsl4jsb_Real2D_onChunk :: emission_noy
    dsl4jsb_Real2D_onChunk :: emission_n2o
    dsl4jsb_Real2D_onChunk :: emission_n2
    dsl4jsb_Real2D_onChunk :: emission_noy_n15
    dsl4jsb_Real2D_onChunk :: emission_n2o_n15
    dsl4jsb_Real2D_onChunk :: emission_n2_n15
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: rtm_decomposition
    dsl4jsb_Real3D_onChunk :: rtm_depolymerisation
    dsl4jsb_Real3D_onChunk :: rtm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rtm_hsc
    dsl4jsb_Real3D_onChunk :: rtm_plant_uptake
    dsl4jsb_Real3D_onChunk :: rmm_depolymerisation
    dsl4jsb_Real3D_onChunk :: rmm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_hsc
    dsl4jsb_Real3D_onChunk :: rmm_plant_uptake
    dsl4jsb_Real3D_onChunk :: rtm_nitrification
    dsl4jsb_Real3D_onChunk :: rmm_nitrification
    dsl4jsb_Real3D_onChunk :: rtm_denitrification
    dsl4jsb_Real3D_onChunk :: rtm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rmm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rtm_asymb_bnf
    dsl4jsb_Real3D_onChunk :: rmm_asymb_bnf
    dsl4jsb_Real3D_onChunk :: anaerobic_volume_fraction_sl
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: nitrification_no3
    dsl4jsb_Real3D_onChunk :: nitrification_noy
    dsl4jsb_Real3D_onChunk :: nitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_noy
    dsl4jsb_Real3D_onChunk :: denitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_n2
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: microbial_cue_eff
    dsl4jsb_Real3D_onChunk :: microbial_nue_eff
    dsl4jsb_Real3D_onChunk :: microbial_pue_eff
    dsl4jsb_Real3D_onChunk :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk :: weathering_po4
    dsl4jsb_Real3D_onChunk :: occlusion_po4
    dsl4jsb_Real3D_onChunk :: slow_exchange_po4
    dsl4jsb_Real3D_onChunk :: fast_adsorpt_po4
    dsl4jsb_Real3D_onChunk :: fast_adsorpt_nh4
    dsl4jsb_Real3D_onChunk :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: gross_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: nh4_solute
    dsl4jsb_Real3D_onChunk :: nh4_assoc
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: po4_primary
    dsl4jsb_Real3D_onChunk :: transport_noy
    dsl4jsb_Real3D_onChunk :: transport_n2o
    dsl4jsb_Real3D_onChunk :: transport_n2
    dsl4jsb_Real3D_onChunk :: noy
    dsl4jsb_Real3D_onChunk :: n2o
    dsl4jsb_Real3D_onChunk :: n2
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: no3_n15_solute
    dsl4jsb_Real3D_onChunk :: noy_n15
    dsl4jsb_Real3D_onChunk :: n2o_n15
    dsl4jsb_Real3D_onChunk :: n2_n15
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute
    dsl4jsb_Real3D_onChunk :: transport_nh4_assoc
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk :: transport_po4_primary
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute
    dsl4jsb_Real3D_onChunk :: transport_noy_n15
    dsl4jsb_Real3D_onChunk :: transport_n2o_n15
    dsl4jsb_Real3D_onChunk :: transport_n2_n15
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_carbon_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_nitrogen_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_phosphorus_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_carbon13_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_carbon14_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_nitrogen15_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_nh4_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_no3_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_po4_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_nh4_n15_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_no3_n15_solute_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: nitrification_no3_n15
    dsl4jsb_Real3D_onChunk :: nitrification_noy_n15
    dsl4jsb_Real3D_onChunk :: nitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_noy_n15
    dsl4jsb_Real3D_onChunk :: denitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_n2_n15
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation_n15
    dsl4jsb_Real3D_onChunk :: k_bioturb
    dsl4jsb_Real3D_onChunk :: qmax_org
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: qmax_nh4
    dsl4jsb_Real3D_onChunk :: qmax_fast_po4
    dsl4jsb_Real3D_onChunk :: km_fast_po4
    dsl4jsb_Real3D_onChunk :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onChunk :: km_adsorpt_nh4_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: vmax_weath_mineral_sl
    dsl4jsb_Real3D_onChunk :: particle_fluxrate
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly
    dsl4jsb_Real3D_onChunk :: enzyme_frac_residue
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_p
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_c_mavg
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_n_mavg
    dsl4jsb_Real3D_onChunk :: enzyme_frac_poly_p_mavg
    dsl4jsb_Real3D_onChunk :: fact_n_status_mic_c_growth
    dsl4jsb_Real3D_onChunk :: fact_p_status_mic_c_growth
    dsl4jsb_Real3D_onChunk :: enzyme_frac_AP
    ! ----------------------------------------------------------------------------------------------------- !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(arr_hlp1(nc, nsoil_sb))
    ALLOCATE(arr_hlp2(nc, nsoil_sb))
    ALLOCATE(arr_hlp3(nc, nsoil_sb))
    ALLOCATE(fine_root_carbon_sl(nc, nsoil_sb))
    ALLOCATE(km1_uptake_act(nc, nsoil_sb))
    ALLOCATE(km2_uptake_act(nc, nsoil_sb))
    ALLOCATE(km_array(nc, nsoil_sb))
    ALLOCATE(avail_soil_surface_sorption_capacity_sl(nc, nsoil_sb))
    ALLOCATE(graze_mineralisation_po4(nc, nsoil_sb))
    ALLOCATE(graze_mineralisation_nh4(nc, nsoil_sb))
    ALLOCATE(graze_mineralisation_nh4_n15(nc, nsoil_sb))
    ALLOCATE(unit_uptake_n_pot_sl(nc, nsoil_sb))
    ALLOCATE(unit_uptake_p_pot_sl(nc, nsoil_sb))
    ALLOCATE(growth_p_limit_smoothed_sl(nc, nsoil_sb))
    ALLOCATE(asymb_n_fixation_rel_rate_sl(nc, nsoil_sb))
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FORMATION_ID, sb_formation_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_LOSS_ID, sb_loss_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_TRANSPORT_ID, sb_transport_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_,   num_sl_above_bedrock)         ! in
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_,   soil_depth_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,   percolation_sl)               ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,   frac_w_lat_loss_sl)           ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,   w_soil_sl)                    ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,   bulk_dens_sl)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_,   volume_min_sl)                !
    ! A2L_ 2D
    dsl4jsb_Get_var2D_onChunk(A2L_,   p_deposition)                 ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_,   growth_p_limit_based_on_n_mavg)   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_,   unit_uptake_n_pot)                !
    dsl4jsb_Get_var2D_onChunk(VEG_,   unit_uptake_p_pot)                !
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_,   root_fraction_sl)             ! in
    ! SB_ 2D
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_dom_carbon)          ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_dom_nitrogen)        ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_dom_phosphorus)      ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_dom_carbon13)        ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_dom_carbon14)        ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_dom_nitrogen15)      ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_nh4_solute)          ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_no3_solute)          ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_po4_solute)          ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_nh4_n15_solute)      ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    leaching_no3_n15_solute)      ! out
    dsl4jsb_Get_var2D_onChunk(SB_,    emission_noy)                 !
    dsl4jsb_Get_var2D_onChunk(SB_,    emission_n2o)                 !
    dsl4jsb_Get_var2D_onChunk(SB_,    emission_n2)                  !
    dsl4jsb_Get_var2D_onChunk(SB_,    emission_noy_n15)             !
    dsl4jsb_Get_var2D_onChunk(SB_,    emission_n2o_n15)             !
    dsl4jsb_Get_var2D_onChunk(SB_,    emission_n2_n15)              !
    ! --------------------------
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_decomposition)                ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_depolymerisation)             ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_mic_uptake)                   ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_sorption)                     ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_desorption)                   ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_hsc)                          ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_plant_uptake)                 ! soil: in
    dsl4jsb_Get_var3D_onChunk(SB_,    rmm_depolymerisation)             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    rmm_mic_uptake)                   ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    rmm_sorption)                     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    rmm_desorption)                   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    rmm_hsc)                          ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    rmm_plant_uptake)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_nitrification)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    rmm_nitrification)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_denitrification)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_gasdiffusion)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    rmm_gasdiffusion)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    rtm_asymb_bnf)                    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    rmm_asymb_bnf)                    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    anaerobic_volume_fraction_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    het_respiration)                  ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    het_respiration_c13)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    het_respiration_c14)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    net_mineralisation_nh4)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    net_mineralisation_nh4_n15)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    net_mineralisation_no3)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    net_mineralisation_no3_n15)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    net_mineralisation_po4)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    nitrification_no3)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    nitrification_noy)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    nitrification_n2o)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    denitrification_noy)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    denitrification_n2o)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    denitrification_n2)               ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    volatilisation_nh4)               ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    asymb_n_fixation)                 ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    plant_uptake_nh4_sl)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    plant_uptake_nh4_n15_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    plant_uptake_no3_sl)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    plant_uptake_no3_n15_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    plant_uptake_po4_sl)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_uptake_nh4_sl)          !
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_uptake_nh4_n15_sl)      !
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_uptake_no3_sl)          !
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_uptake_no3_n15_sl)      !
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_uptake_po4_sl)          !
    dsl4jsb_Get_var3D_onChunk(SB_,    mycorrhiza_uptake_nh4_sl)         !
    dsl4jsb_Get_var3D_onChunk(SB_,    mycorrhiza_uptake_no3_sl)         !
    dsl4jsb_Get_var3D_onChunk(SB_,    mycorrhiza_uptake_norg_sl)        !
    dsl4jsb_Get_var3D_onChunk(SB_,    mycorrhiza_uptake_norg_n15_sl)    !
    dsl4jsb_Get_var3D_onChunk(SB_,    mycorrhiza_uptake_po4_sl)         !
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_cue_eff)                !
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_nue_eff)                !
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_pue_eff)                !
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_cue_eff_tmic_mavg)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_nue_eff_tmic_mavg)      !
    dsl4jsb_Get_var3D_onChunk(SB_,    microbial_pue_eff_tmic_mavg)      !
    dsl4jsb_Get_var3D_onChunk(SB_,    weathering_po4)                   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    occlusion_po4)                    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    slow_exchange_po4)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    fast_adsorpt_po4)                 ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    fast_adsorpt_nh4)                 ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    biochem_mineralisation_po4)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    gross_mineralisation_po4)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    nh4_solute)                       ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    nh4_assoc)                        ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    no3_solute)                       ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    po4_solute)                       ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    po4_assoc_fast)                   ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    po4_assoc_slow)                   ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    po4_occluded)                     ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    po4_primary)                      ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    nh4_n15_solute)                   ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    nh4_n15_assoc)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    no3_n15_solute)                   ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    noy_n15)                          !
    dsl4jsb_Get_var3D_onChunk(SB_,    n2o_n15)                          !
    dsl4jsb_Get_var3D_onChunk(SB_,    n2_n15)                           !
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_noy)                    !
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_n2o)                    !
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_n2)                     !
    dsl4jsb_Get_var3D_onChunk(SB_,    noy)                              !
    dsl4jsb_Get_var3D_onChunk(SB_,    n2o)                              !
    dsl4jsb_Get_var3D_onChunk(SB_,    n2)                               !
    dsl4jsb_Get_var3D_onChunk(SB_,    mycorrhiza_uptake_nh4_n15_sl)     !
    dsl4jsb_Get_var3D_onChunk(SB_,    mycorrhiza_uptake_no3_n15_sl)     !
    dsl4jsb_Get_var3D_onChunk(SB_,    nitrification_no3_n15)            !
    dsl4jsb_Get_var3D_onChunk(SB_,    nitrification_noy_n15)            !
    dsl4jsb_Get_var3D_onChunk(SB_,    nitrification_n2o_n15)            !
    dsl4jsb_Get_var3D_onChunk(SB_,    denitrification_noy_n15)          !
    dsl4jsb_Get_var3D_onChunk(SB_,    denitrification_n2o_n15)          !
    dsl4jsb_Get_var3D_onChunk(SB_,    denitrification_n2_n15)           !
    dsl4jsb_Get_var3D_onChunk(SB_,    volatilisation_nh4_n15)           !
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_noy_n15)                !
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_n2o_n15)                !
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_n2_n15)                 !
    dsl4jsb_Get_var3D_onChunk(SB_,    asymb_n_fixation_n15)             !
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_nh4_solute)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_nh4_assoc)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_no3_solute)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_po4_solute)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_po4_assoc_fast)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_po4_assoc_slow)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_po4_occluded)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_po4_primary)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_nh4_n15_solute)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_nh4_n15_assoc)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    transport_no3_n15_solute)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_dom_carbon_sl)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_dom_nitrogen_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_dom_phosphorus_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_dom_carbon13_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_dom_carbon14_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_dom_nitrogen15_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_nh4_solute_sl)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_no3_solute_sl)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_po4_solute_sl)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_nh4_n15_solute_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    lateral_loss_no3_n15_solute_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    k_bioturb)                        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    qmax_org)                         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    qmax_po4)                         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    qmax_nh4)                         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    qmax_fast_po4)                    !
    dsl4jsb_Get_var3D_onChunk(SB_,    km_fast_po4)                      !
    dsl4jsb_Get_var3D_onChunk(SB_,    km_adsorpt_po4_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    km_adsorpt_nh4_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    bulk_dens_corr_sl)                ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    vmax_weath_mineral_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(SB_,    particle_fluxrate)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_residue)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly_c)               ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly_n)               ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly_p)               ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly_c_mavg)          ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly_n_mavg)          ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_poly_p_mavg)          ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    fact_n_status_mic_c_growth)       ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    fact_p_status_mic_c_growth)       ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,    enzyme_frac_AP)                   ! inout
    !------------------------------------------------------------------------------------------------------ !

    !>1.0 calc fine root carbon, used for N & P calculations below, and smooth the plant growth P limitation
    !>
    DO isoil = 1,nsoil_sb
      WHERE(soil_depth_sl(:,isoil) > eps8)
        fine_root_carbon_sl(:,isoil) = root_fraction_sl(:,isoil) / soil_depth_sl(:,isoil) * veg_pool_mt(ix_fine_root, ixC, :)
        ! using LOG() here to smooth the response curve in case of extreme P limited condition,
        ! avoiding drastic changes in plant phosphatase productions as it is updated at each time step
        growth_p_limit_smoothed_sl(:,isoil) = MAX(0.0_wp, LOG(growth_p_limit_based_on_n_mavg(:)))
      ELSEWHERE
        fine_root_carbon_sl(:,isoil)        = 0.0_wp
        growth_p_limit_smoothed_sl(:,isoil) = 0.0_wp
      ENDWHERE
    ENDDO
    unit_uptake_n_pot(:) = 0.0_wp
    unit_uptake_p_pot(:) = 0.0_wp

    !>2.0 Inorganic N assimilation
    !>

    !>  2.1 PO4
    !>
    !! Calculate the P flux rates, and potential microbial N,P uptake for eca_none
    !!   NOTE: it is an ELEMENTAL routine, hence, not using a bgcm matrices as argument
    !!
    CALL calc_sourcing_flux( &
      & dsl4jsb_Config(SB_)%flag_sb_prescribe_po4, &  ! in
      & sb_pool_mt(ix_microbial, ixC, :, :), &
      & sb_pool_mt(ix_microbial, ixN, :, :), &
      & sb_pool_mt(ix_microbial, ixP, :, :), &
      & sb_pool_mt(ix_dom, ixC, :, :), &
      & sb_pool_mt(ix_dom, ixN, :, :), &
      & sb_pool_mt(ix_dom, ixP, :, :), &
      & rtm_mic_uptake(:,:), &
      & rmm_mic_uptake(:,:), &
      & rtm_hsc(:,:), &
      & rmm_hsc(:,:), &
      & rtm_sorption(:,:), &
      & rmm_sorption(:,:), &
      & rtm_desorption(:,:), &
      & rmm_desorption(:,:), &
      & sb_pool_mt(ix_dom_assoc, ixP, :, :), &
      & sb_pool_mt(ix_residue_assoc, ixP, :, :), &
      & sb_pool_mt(ix_residue, ixP, :, :), &
      & sb_pool_mt(ix_dom_assoc, ixC, :, :), &
      & sb_pool_mt(ix_residue_assoc, ixC, :, :), &
      & sb_pool_mt(ix_residue, ixC, :, :), &
      & po4_solute(:,:), &
      & nh4_solute(:,:), &
      & no3_solute(:,:), &
      & po4_assoc_fast(:,:), &
      & po4_assoc_slow(:,:), &
      & po4_primary(:,:), &
      & bulk_dens_corr_sl(:,:), &
      & volume_min_sl(:,:), &
      & fine_root_carbon_sl(:,:), &
      & vmax_weath_mineral_sl, &
      & enzyme_frac_AP(:,:), &
      & growth_p_limit_smoothed_sl(:,:), &            ! in
      & sb_loss_mt(ix_dom_assoc, ixP, :, :), &        ! inout
      & sb_loss_mt(ix_residue_assoc, ixP, :, :), &
      & sb_loss_mt(ix_residue, ixP, :, :), &
      & sb_loss_mt(ix_dom, ixP, :, :), &              ! inout
      & biochem_mineralisation_po4(:,:), &            ! out
      & weathering_po4(:,:), &
      & microbial_uptake_nh4_sl(:,:), &
      & microbial_uptake_no3_sl(:,:), &
      & microbial_uptake_po4_sl(:,:), &
      & gross_mineralisation_po4(:,:), &
      & occlusion_po4(:,:), &
      & slow_exchange_po4(:,:) )                      ! out

    !>3.0 Litter-related fluxes
    !>

    !>  3.1 decay of soluable litter
    !>
    CALL calc_litter_turnover( &
      & nc, &                                           ! in
      & nsoil_sb, &
      & tau_soluable_litter, &
      & fc_soluable2dom, &
      & rtm_depolymerisation(:,:), &
      & rmm_depolymerisation(:,:), &
      & sb_pool_mt(ix_soluable_litter, :, :, :), &      ! in
      & het_respiration(:,:), &                         ! inout
      & het_respiration_c13(:,:), &
      & het_respiration_c14(:,:), &
      & sb_loss_mt(ix_soluable_litter, :, :, :), &
      & sb_formation_mt(ix_dom, :, :, :) )              ! inout

    !>  3.2 decay of woody litter
    !>
    CALL calc_litter_turnover( &
      & nc, &                                           ! in
      & nsoil_sb, &
      & tau_woody_litter, &
      & fc_woody2polymeric, &
      & rtm_depolymerisation(:,:), &
      & rmm_depolymerisation(:,:), &
      & sb_pool_mt(ix_woody_litter, :, :, :), &         ! in
      & het_respiration(:,:), &                         ! inout
      & het_respiration_c13(:,:), &
      & het_respiration_c14(:,:), &
      & sb_loss_mt(ix_woody_litter,:,:,:), &
      & sb_formation_mt(ix_polymeric_litter,:,:,:) )    ! inout

    !>4.0 Update soil litter pools with litter fall
    !>
    CALL calc_litter_partitioning( &
      & nc, &                                         ! in
      & nsoil_sb, &
      & num_sl_above_bedrock(:), &
      & lctlib%sla, &
      & lctlib%growthform, &
      & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
      & soil_depth_sl(:,:), &
      & root_fraction_sl(:,:), &
      & veg_litterfall_mt(:,:,:), &                   ! in
      & sb_formation_mt(:,:,:,:) )                    ! inout

    !>5.0 vertical transport
    !>

    !>  5.1 liquid phase transport
    !>
    CALL calc_liquid_phase_transport_wrapper( &
      & nc, &                                   ! in
      & nsoil_sb, &
      & num_sl_above_bedrock(:), &
      & soil_depth_sl(:,:), &
      & rtm_gasdiffusion(:,:), &
      & rmm_gasdiffusion(:,:), &
      & percolation_sl(:,:), &
      & frac_w_lat_loss_sl(:,:), &
      & nh4_solute(:,:), &
      & no3_solute(:,:), &
      & po4_solute(:,:), &
      & nh4_n15_solute(:,:), &
      & no3_n15_solute(:,:), &
      & noy(:,:), &
      & noy_n15(:,:), &
      & n2o(:,:), &
      & n2o_n15(:,:), &
      & n2(:,:), &
      & n2_n15(:,:), &
      & sb_pool_mt(ix_dom, :, :, :), &          ! in
      & transport_nh4_solute(:,:), &            ! inout
      & transport_no3_solute(:,:), &
      & transport_po4_solute(:,:), &
      & transport_nh4_n15_solute(:,:), &
      & transport_no3_n15_solute(:,:), &
      & leaching_dom_carbon(:), &
      & leaching_dom_nitrogen(:), &
      & leaching_dom_phosphorus(:), &
      & leaching_dom_carbon13(:), &
      & leaching_dom_carbon14(:), &
      & leaching_dom_nitrogen15(:), &
      & leaching_nh4_solute(:), &
      & leaching_no3_solute(:), &
      & leaching_po4_solute(:), &
      & leaching_nh4_n15_solute(:), &
      & leaching_no3_n15_solute(:), &
      & sb_transport_mt(ix_dom, :, :, :), &     ! inout
      & lateral_loss_dom_carbon_sl(:,:), &      ! out
      & lateral_loss_dom_nitrogen_sl(:,:), &
      & lateral_loss_dom_phosphorus_sl(:,:), &
      & lateral_loss_dom_carbon13_sl(:,:), &
      & lateral_loss_dom_carbon14_sl(:,:), &
      & lateral_loss_dom_nitrogen15_sl(:,:), &
      & lateral_loss_nh4_solute_sl(:,:), &
      & lateral_loss_no3_solute_sl(:,:), &
      & lateral_loss_po4_solute_sl(:,:), &
      & lateral_loss_nh4_n15_solute_sl(:,:), &
      & lateral_loss_no3_n15_solute_sl(:,:), &
      & transport_noy(:,:), &
      & transport_noy_n15(:,:), &
      & transport_n2o(:,:), &
      & transport_n2o_n15(:,:), &
      & transport_n2(:,:), &
      & transport_n2_n15(:,:), &
      & emission_noy(:), &
      & emission_noy_n15(:), &
      & emission_n2o(:), &
      & emission_n2o_n15(:), &
      & emission_n2(:), &
      & emission_n2_n15(:) )                    ! out

    !! include the lateral loss of DOM into sb_loss_mt(ix_dom, :, :)
    sb_loss_mt(ix_dom, ixC, :, :)   = sb_loss_mt(ix_dom, ixC, :, :)   + lateral_loss_dom_carbon_sl(:,:)
    sb_loss_mt(ix_dom, ixC13, :, :) = sb_loss_mt(ix_dom, ixC13, :, :) + lateral_loss_dom_carbon13_sl(:,:)
    sb_loss_mt(ix_dom, ixC14, :, :) = sb_loss_mt(ix_dom, ixC14, :, :) + lateral_loss_dom_carbon14_sl(:,:)
    sb_loss_mt(ix_dom, ixN, :, :)   = sb_loss_mt(ix_dom, ixN, :, :)   + lateral_loss_dom_nitrogen_sl(:,:)
    sb_loss_mt(ix_dom, ixN15, :, :) = sb_loss_mt(ix_dom, ixN15, :, :) + lateral_loss_dom_nitrogen15_sl(:,:)
    sb_loss_mt(ix_dom, ixP, :, :)   = sb_loss_mt(ix_dom, ixP, :, :)   + lateral_loss_dom_phosphorus_sl(:,:)

    !>  5.2 transport due to bioturbation
    !>
    k_bioturb(:,:) = calc_bioturbation_rate(nc, nsoil_sb, num_sl_above_bedrock(:), soil_depth_sl(:,:), &
      &                                     bulk_dens_corr_sl(:,:), root_fraction_sl(:,:))
    ! bioturbation transport for JSM
    CALL calc_bioturbation_transport_wrapper_jsm( &
      & nc, &                                 ! in
      & nsoil_sb, &
      & num_sl_above_bedrock(:), &
      & soil_depth_sl(:,:), &
      & model%config%elements_index_map(:), &
      & model%config%is_element_used(:), &
      & k_bioturb(:,:), &
      & nh4_assoc(:,:), &
      & po4_assoc_fast(:,:), &
      & po4_assoc_slow(:,:), &
      & po4_occluded(:,:), &
      & nh4_n15_assoc(:,:), &
      & sb_pool_mt(:,:,:,:), &                ! in
      & transport_nh4_assoc(:,:), &           ! inout
      & transport_po4_assoc_fast(:,:), &
      & transport_po4_assoc_slow(:,:), &
      & transport_po4_occluded(:,:), &
      & transport_nh4_n15_assoc(:,:), &
      & sb_transport_mt(:,:,:,:) )            ! inout

    !>6.0 Microbial and plant P uptake
    !>

    !>  6.1 PO4
    !>
    arr_hlp1(:,:)        = po4_solute(:,:) + po4_assoc_fast(:,:)
    arr_hlp2(:,:)        = 1._wp
    km1_uptake_act(:,:)  = km1_up_po4 / (rtm_hsc(:,:) * rmm_hsc(:,:)) / 1000._wp / (w_soil_sl(:,:) / soil_depth_sl(:,:))
    km2_uptake_act(:,:)  = km2_up_po4 *  rtm_hsc(:,:) * rmm_hsc(:,:)  * 1000._wp * (w_soil_sl(:,:) / soil_depth_sl(:,:))
    km_array(:,:)        = km_adsorpt_po4_sl(:,:)
    avail_soil_surface_sorption_capacity_sl(:,:) = MAX(0._wp, qmax_po4(:,:) - po4_assoc_fast(:,:))
    CALL calc_sinking_flux(lctlib%vmax_uptake_p, "po4", TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &                       ! in
      &                    sb_pool_mt(ix_microbial, ixC, :, :), fine_root_carbon_sl(:,:), po4_solute(:,:), &                ! in
      &                    rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), &                          ! in
      &                    rtm_plant_uptake(:,:), rmm_plant_uptake(:,:), km1_uptake_act(:,:), km2_uptake_act(:,:), &        ! in
      &                    km_denitrification_c, km_denitrification_no3, rtm_nitrification(:,:), rmm_nitrification(:,:), &  ! in
      &                    arr_hlp2(:,:), km_array(:,:), avail_soil_surface_sorption_capacity_sl(:,:), &                    ! in
      &                    solute_available = arr_hlp1(:,:), &                        ! in
      &                    microbial_uptake = microbial_uptake_po4_sl(:,:), &         ! inout
      &                    plant_uptake     = plant_uptake_po4_sl(:,:), &             ! out
      &                    plant_uptake_pot = unit_uptake_p_pot_sl(:,:), &            ! out
      &                    fast_adsorpt     = fast_adsorpt_po4(:,:))                  ! OPTIONAL (out)
    !! calculate potential P uptake (micro-mol P / mol C fine-root / s), weighted by root distribution
    DO ic = 1,nc
      DO isoil = 1,nsoil_sb
        IF (fine_root_carbon_sl(ic,isoil) > eps8) THEN
          unit_uptake_p_pot(ic) = unit_uptake_p_pot(ic) + unit_uptake_p_pot_sl(ic,isoil) &
            &                     / fine_root_carbon_sl(ic,isoil) * root_fraction_sl(ic,isoil)
        ENDIF
      ENDDO
    ENDDO

    !>  6.2 NH4 - Microbial and plant N uptake, nitrification/denitrification
    !>
    arr_hlp1(:,:)        = nh4_solute(:,:) + nh4_assoc(:,:)
    arr_hlp2(:,:)        = 1._wp - anaerobic_volume_fraction_sl(:,:)
    km1_uptake_act(:,:)  = km1_up_nh4 / (rtm_hsc(:,:) * rmm_hsc(:,:)) / 1000._wp / (w_soil_sl(:,:) / soil_depth_sl(:,:))
    km2_uptake_act(:,:)  = km2_up_nh4 *  rtm_hsc(:,:) * rmm_hsc(:,:)  * 1000._wp * (w_soil_sl(:,:) / soil_depth_sl(:,:))
    km_array(:,:)        = km_adsorpt_nh4_sl(:,:)
    avail_soil_surface_sorption_capacity_sl(:,:) = MAX(0._wp, qmax_nh4(:,:) - nh4_assoc(:,:))
    CALL calc_sinking_flux(lctlib%vmax_uptake_n, "nh4", TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &                       ! in
      &                    sb_pool_mt(ix_microbial, ixC, :, :), fine_root_carbon_sl(:,:), nh4_solute(:,:), &                ! in
      &                    rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), &                          ! in
      &                    rtm_plant_uptake(:,:), rmm_plant_uptake(:,:), km1_uptake_act(:,:), km2_uptake_act(:,:), &        ! in
      &                    km_denitrification_c, km_denitrification_no3, rtm_nitrification(:,:), rmm_nitrification(:,:), &  ! in
      &                    arr_hlp2(:,:), km_array(:,:), avail_soil_surface_sorption_capacity_sl(:,:), &                    ! in
      &                    solute_available   = arr_hlp1(:,:), &                      ! in
      &                    vmax_nitri_denitri = vmax_nitrification, &                 ! OPTIONAL(in)
      &                    solute_alternative = no3_solute(:,:), &                    ! OPTIONAL(in)
      &                    microbial_uptake   = microbial_uptake_nh4_sl(:,:), &       ! inout
      &                    plant_uptake       = plant_uptake_nh4_sl(:,:), &           ! out
      &                    plant_uptake_pot   = unit_uptake_n_pot_sl(:,:), &          ! out
      &                    fast_adsorpt       = fast_adsorpt_nh4(:,:), &              ! OPTIONAL (out)
      &                    flux_nitri_denitri = nitrification_no3(:,:))               ! OPTIONAL (out)
    !! update potential N uptake due to nh4 uptake
    DO ic = 1,nc
      DO isoil = 1,nsoil_sb
        IF (fine_root_carbon_sl(ic,isoil) > eps8) THEN
          unit_uptake_n_pot(ic) = unit_uptake_n_pot(ic) + unit_uptake_n_pot_sl(ic,isoil) &
            &                     / fine_root_carbon_sl(ic,isoil) * root_fraction_sl(ic,isoil)
        ENDIF
      ENDDO
    ENDDO

    !>  6.3 NO3
    !>
    arr_hlp1(:,:)        = 1._wp
    km1_uptake_act(:,:)  = km1_up_no3 / (rtm_hsc(:,:) * rmm_hsc(:,:)) / 1000._wp / (w_soil_sl(:,:) / soil_depth_sl(:,:))
    km2_uptake_act(:,:)  = km2_up_no3 *  rtm_hsc(:,:) * rmm_hsc(:,:)  * 1000._wp * (w_soil_sl(:,:) / soil_depth_sl(:,:))
    km_array(:,:)        = 1._wp
    avail_soil_surface_sorption_capacity_sl(:,:) = 1._wp
    CALL calc_sinking_flux(lctlib%vmax_uptake_n, "no3", TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &                         ! in
      &                    sb_pool_mt(ix_microbial, ixC, :, :), fine_root_carbon_sl(:,:), no3_solute(:,:), &                  ! in
      &                    rtm_mic_uptake(:,:), rmm_mic_uptake(:,:), rtm_hsc(:,:), rmm_hsc(:,:), &                            ! in
      &                    rtm_plant_uptake(:,:), rmm_plant_uptake(:,:), km1_uptake_act(:,:), km2_uptake_act(:,:), &          ! in
      &                    km_denitrification_c, km_denitrification_no3, rtm_denitrification(:,:), arr_hlp1(:,:), &           ! in
      &                    anaerobic_volume_fraction_sl(:,:), km_array(:,:), avail_soil_surface_sorption_capacity_sl(:,:), &  ! in
      &                    solute_available   = no3_solute(:,:), &                    ! in
      &                    vmax_nitri_denitri = vmax_denitrification, &               ! OPTIONAL(in)
      &                    solute_alternative = nh4_solute(:,:), &                    ! OPTIONAL(in)
      &                    microbial_uptake   = microbial_uptake_no3_sl(:,:), &       ! inout
      &                    plant_uptake       = plant_uptake_no3_sl(:,:), &           ! out
      &                    plant_uptake_pot   = unit_uptake_n_pot_sl(:,:), &          ! out
      &                    flux_nitri_denitri = denitrification_n2(:,:))              ! OPTIONAL (out)
    !! update potential N uptake due to no3 uptake
    DO ic = 1,nc
      DO isoil = 1,nsoil_sb
        IF (fine_root_carbon_sl(ic,isoil) > eps8) THEN
          unit_uptake_n_pot(ic) = unit_uptake_n_pot(ic) + unit_uptake_n_pot_sl(ic,isoil) &
            &                     / fine_root_carbon_sl(ic,isoil) * root_fraction_sl(ic,isoil)
        ENDIF
      ENDDO
    ENDDO

    !>7.0 partition nitrification and denitrication rate
    !>
    CALL calc_nitri_denitri_partitioning_rate( &
      & nc, &
      & nsoil_sb, &
      & soil_depth_sl(:,:), &
      & TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme), &
      & rtm_denitrification(:,:), &
      & nitrification_no3(:,:), &
      & denitrification_n2(:,:), &
      & nitrification_noy(:,:), &
      & nitrification_n2o(:,:), &
      & denitrification_noy(:,:), &
      & denitrification_n2o(:,:), &
      & net_mineralisation_nh4(:,:), &
      & net_mineralisation_nh4_n15(:,:), &
      & nitrification_no3_n15(:,:), &
      & volatilisation_nh4(:,:), &
      & volatilisation_nh4_n15(:,:), &
      & emission_n2(:), &
      & emission_n2_n15(:) )

    !>  7.1 update n15 for inorganic fluxes (need to be done before microbial turnover and growth)
    !>
    ! here the function is called for JSM: flux rates in [micro-mol/m3/s] (unit conversion is applied within the routine)
    ! (simple_1d is using the same function but uses flux rates in [mol/m3/timestep])
    CALL calc_sb_inorganic_n15_fluxes( &
      & nc, &
      & nsoil_sb, &
      & num_sl_above_bedrock(:), &
      & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
      & nh4_solute(:,:), &
      & nh4_n15_solute(:,:), &
      & no3_solute(:,:), &
      & no3_n15_solute(:,:), &
      & plant_uptake_nh4_sl(:,:), &           ! input flux rates that get unit coverted
      & plant_uptake_no3_sl(:,:), &           ! input flux rates that get unit coverted
      & mycorrhiza_uptake_nh4_sl(:,:), &      ! input flux rates that get unit coverted
      & mycorrhiza_uptake_no3_sl(:,:), &      ! input flux rates that get unit coverted
      & microbial_uptake_nh4_sl(:,:), &       ! input flux rates that get unit coverted
      & microbial_uptake_no3_sl(:,:), &       ! input flux rates that get unit coverted
      & asymb_n_fixation(:,:), &              ! input flux rates that get unit coverted
      & nitrification_no3(:,:), &             ! input flux rates that get unit coverted
      & nitrification_noy(:,:), &             ! input flux rates that get unit coverted
      & nitrification_n2o(:,:), &             ! input flux rates that get unit coverted
      & denitrification_noy(:,:), &           ! input flux rates that get unit coverted
      & denitrification_n2o(:,:), &           ! input flux rates that get unit coverted
      & denitrification_n2(:,:), &            ! input flux rates that get unit coverted
      & transport_nh4_n15_solute(:,:), &      ! in
      & transport_no3_n15_solute(:,:), &      ! in
      & plant_uptake_nh4_n15_sl(:,:), &       ! output var that get unit coverted
      & plant_uptake_no3_n15_sl(:,:), &       ! output var that get unit coverted
      & mycorrhiza_uptake_nh4_n15_sl(:,:), &  ! output var that get unit coverted
      & mycorrhiza_uptake_no3_n15_sl(:,:), &  ! output var that get unit coverted
      & microbial_uptake_nh4_n15_sl(:,:), &   ! output var that get unit coverted
      & microbial_uptake_no3_n15_sl(:,:), &   ! output var that get unit coverted
      & asymb_n_fixation_n15(:,:), &          ! output var that get unit coverted
      & nitrification_no3_n15(:,:), &         ! output var that get unit coverted
      & nitrification_noy_n15(:,:), &         ! output var that get unit coverted
      & nitrification_n2o_n15(:,:), &         ! output var that get unit coverted
      & denitrification_noy_n15(:,:), &       ! output var that get unit coverted
      & denitrification_n2o_n15(:,:), &       ! output var that get unit coverted
      & denitrification_n2_n15(:,:) )         ! output var that get unit coverted

    !>8.0 co-metabolisation of polymeric litter & residues
    !>
    arr_hlp1(:,:) = sb_pool_mt(ix_microbial, ixC, :, :) * enzyme_frac_total * (1._wp - enzyme_frac_AP(:,:))
    CALL calc_depolymerisation_SESAM( &
      & nc, &                                 ! in
      & nsoil_sb, &
      & rtm_hsc(:,:), &
      & rmm_hsc(:,:), &
      & rtm_depolymerisation(:,:), &
      & rmm_depolymerisation(:,:), &
      & arr_hlp1(:,:), &
      & microbial_cue_eff_tmic_mavg(:,:), &
      & microbial_nue_eff_tmic_mavg(:,:), &
      & microbial_pue_eff_tmic_mavg(:,:), &
      & fact_n_status_mic_c_growth(:,:), &
      & fact_p_status_mic_c_growth(:,:), &
      & sb_pool_mt(:,:,:,:), &                ! in
      & enzyme_frac_poly(:,:), &              ! inout
      & enzyme_frac_residue(:,:), &
      & enzyme_frac_AP(:,:), &
      & enzyme_frac_poly_c(:,:), &
      & enzyme_frac_poly_n(:,:), &
      & enzyme_frac_poly_p(:,:), &
      & enzyme_frac_poly_c_mavg(:,:), &
      & enzyme_frac_poly_n_mavg(:,:), &
      & enzyme_frac_poly_p_mavg(:,:), &
      & sb_loss_mt(:,:,:,:), &
      & sb_formation_mt(:,:,:,:) )            ! inout

    !>9.0 microbial death
    !>  determine mortality rate, constrained by avoiding a complete death of the community
    !>
    CALL calc_microbial_turnover( &
      & nc, &                                 ! in
      & nsoil_sb, &
      & fact_n_status_mic_c_growth(:,:), &
      & fact_p_status_mic_c_growth(:,:), &
      & sb_pool_mt(ix_microbial,:,:,:), &     ! in
      & sb_loss_mt(ix_microbial,:,:,:), &     ! inout
      & sb_formation_mt(:,:,:,:) )            ! inout

    !>10.0 mineral stabilisation
    !>

    !>  10.1 calculates sorption and desorption of organic material
    !>
    CALL calc_sorption_desorption_of_org_material( &
      & nc, &                                 ! in
      & nsoil_sb, &
      & model%config%elements_index_map(:), &
      & model%config%is_element_used(:), &
      & qmax_org(:,:), &
      & rtm_sorption(:,:), &
      & rmm_sorption(:,:), &
      & rtm_desorption(:,:), &
      & rmm_desorption(:,:), &
      & sb_pool_mt(:,:,:,:), &                ! in
      & sb_formation_mt(:,:,:,:), &           ! inout
      & sb_loss_mt(:,:,:,:) )                 ! inout

    !>11.0 microbial growth
    !>

    !>  11.1 asymbiotic N fixation rate (relativ to N demand)
    !>
    asymb_n_fixation_rel_rate_sl(:,:) = calc_asymb_bnf_fraction(rtm_asymb_bnf(:,:), &
      &                                                         rmm_asymb_bnf(:,:), &
      &                                                         nh4_solute(:,:), &
      &                                                         no3_solute(:,:) )

    !>  11.2 actual microbial growth
    !>
    CALL calc_microbial_growth(nc, &
      & nsoil_sb, &                               ! in
      & TRIM(dsl4jsb_Config(SB_)%bnf_scheme), &
      & rtm_mic_uptake(:,:), &
      & rmm_mic_uptake(:,:), &
      & rtm_hsc(:,:), &
      & rmm_hsc(:,:), &
      & microbial_cue_eff_tmic_mavg(:,:), &
      & asymb_n_fixation_rel_rate_sl(:,:), &
      & sb_pool_mt(:,:,:,:), &                    ! in
      & microbial_uptake_nh4_sl(:,:), &           ! inout
      & microbial_uptake_nh4_n15_sl(:,:), &
      & microbial_uptake_no3_sl(:,:), &
      & microbial_uptake_no3_n15_sl(:,:), &
      & microbial_uptake_po4_sl(:,:), &
      & het_respiration(:,:), &
      & het_respiration_c13(:,:), &
      & het_respiration_c14(:,:), &
      & net_mineralisation_nh4(:,:), &
      & net_mineralisation_nh4_n15(:,:), &
      & net_mineralisation_no3(:,:), &
      & net_mineralisation_no3_n15(:,:), &
      & gross_mineralisation_po4(:,:), &
      & net_mineralisation_po4(:,:), &
      & microbial_cue_eff(:,:), &
      & microbial_nue_eff(:,:), &
      & microbial_pue_eff(:,:), &
      & asymb_n_fixation(:,:), &
      & fact_n_status_mic_c_growth(:,:), &
      & fact_p_status_mic_c_growth(:,:), &
      & sb_loss_mt(:,:,:,:), &
      & sb_formation_mt(:,:,:,:) )                ! inout

    !------------------------------------------------------------------------------------------------------ !
    ! de-allocate local allocatable variables
    DEALLOCATE(arr_hlp1)
    DEALLOCATE(arr_hlp2)
    DEALLOCATE(arr_hlp3)
    DEALLOCATE(fine_root_carbon_sl)
    DEALLOCATE(km1_uptake_act)
    DEALLOCATE(km2_uptake_act)
    DEALLOCATE(km_array)
    DEALLOCATE(avail_soil_surface_sorption_capacity_sl)
    DEALLOCATE(graze_mineralisation_po4)
    DEALLOCATE(graze_mineralisation_nh4)
    DEALLOCATE(graze_mineralisation_nh4_n15)
    DEALLOCATE(unit_uptake_n_pot_sl)
    DEALLOCATE(unit_uptake_p_pot_sl)
    DEALLOCATE(growth_p_limit_smoothed_sl)
    DEALLOCATE(asymb_n_fixation_rel_rate_sl)

  END SUBROUTINE update_sb_jsm

#endif
END MODULE mo_q_sb_jsm_main
