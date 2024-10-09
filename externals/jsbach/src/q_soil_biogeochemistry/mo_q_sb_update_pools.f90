!> QUINCY update soil-biogeochemical pools
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
!>#### routines for calculating the soil-biogeochemical pools using the fluxes that have been calculated
!>
MODULE mo_q_sb_update_pools
#ifndef __NO_QUINCY__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_jsb_math_constants,      ONLY: eps8, one_day, one_year
  USE mo_jsb_physical_constants,  ONLY: lambda_C14

  USE mo_sb_config_class,         ONLY: get_number_of_sb_compartments

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID, SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, &
    &                                   SB_BGCM_TRANSPORT_ID, SB_BGCM_DELTA_POOLS_ID, SB_BGCM_MYCO_EXPORT_ID
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_LITTERFALL_ID, VEG_BGCM_EXUDATION_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_pools

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sb_update_pools'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  !
  !-----------------------------------------------------------------------------------------------------
  !> update soil biogeochemical pools
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_pools(tile, options)

    USE mo_jsb_class,              ONLY: Get_model
    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,         ONLY: t_jsb_task_options
    USE mo_jsb_model_class,        ONLY: t_jsb_model
    USE mo_quincy_model_config,    ONLY: QLAND, QSOIL
    USE mo_jsb_lctlib_class,       ONLY: t_lctlib_element
    USE mo_jsb_task_class,         ONLY: t_jsb_task_options
    USE mo_jsb_process_class,      ONLY: A2L_, L2A_, SEB_, TURB_, SPQ_, HYDRO_, VEG_, HD_, Q_RAD_, Q_ASSIMI_, Q_PHENO_, SB_
    USE mo_jsb_grid_class,         ONLY: t_jsb_vgrid
    USE mo_jsb_grid,               ONLY: Get_vgrid
    USE mo_jsb_math_constants,     ONLY: eps1, eps8
    USE mo_sb_constants
    USE mo_q_sb_jsm_processes,     ONLY: calc_bulk_soil_carbon, calc_bulk_density_correction, calc_Psorption_parameter, &
                                         calc_qmax_bulk_density_correction
    USE mo_q_sb_jsm_transport,     ONLY: calc_particle_fluxrate, calc_particle_transport_wrapper
    USE mo_isotope_util,           ONLY: calc_mixing_ratio_N15N14
    USE mo_spq_util,               ONLY: calc_qmax_texture
    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)

    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4jsb_Use_memory(SPQ_)       !  USE mo_spq_memory_class,     ONLY: t_spq_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_spq_memory
    dsl4jsb_Use_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model                !< the model
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store           !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_soil_sb        !< Vertical grid
    INTEGER                               :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    INTEGER                               :: isoil                !< loop over soil layers
    REAL(wp), ALLOCATABLE, DIMENSION(:)   :: root_depth_avg       !< dim: nc
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: p_deposition_sl      !< ..
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: f_myc_sl             !< ..
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: hlp1, hlp2, hlp3     !< ..
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: km_adsorpt_po4_act, qmax_po4_mineral_act, &
                                             km_fast_po4_act, qmax_fast_po4_act, &
                                             km_slow_po4_act, qmax_slow_po4_act, &
                                             km_adsorpt_nh4_act, qmax_nh4_mineral_act, &
                                             k_desorpt_po4_act, k_adsorpt_po4_act, &
                                             rtm_sorption_act, rmm_sorption_act, &
                                             rtm_desorption_act, rmm_desorption_act, &
                                             partition_coef
    REAL(wp)                              :: dtime                          !< timestep length
    INTEGER                               :: iblk, ics, ice, nc, ic, is     !< grid dimensions
    INTEGER                               :: nr_elements                    !< number of elements
    INTEGER                               :: nr_sb_parts                    !< number of sb compartments
    INTEGER                               :: nr_veg_parts                   !< number of veg compartments
    REAL(wp), ALLOCATABLE                 :: sb_pool_prev_timestep_hlp_mt(:,:,:,:) !< dim: compartments, elements, nc, nsoil
    REAL(wp), ALLOCATABLE                 :: sb_particle_transport_hlp_mt(:,:,:,:) !< dim: compartments, elements, nc, nsoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt1L2D :: veg_exudation_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_formation_mt
    dsl4jsb_Def_mt2L3D :: sb_loss_mt
    dsl4jsb_Def_mt2L3D :: sb_transport_mt
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt
    dsl4jsb_Def_mt1L3D :: sb_delta_pools_mt
    ! ----------------------------------------------------------------------------------------------------- !
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    dsl4jsb_Def_memory(SPQ_)                             ! TYPE(t_spq_memory),     POINTER :: spq__mem
    dsl4jsb_Def_memory(VEG_)                              ! TYPE(t_spq_memory),     POINTER :: spq__mem
    dsl4jsb_Def_memory(A2L_)
    ! Declare pointers to variables in memory
    ! A2L
    dsl4jsb_Real2D_onChunk :: nhx_deposition
    dsl4jsb_Real2D_onChunk :: noy_deposition
    dsl4jsb_Real2D_onChunk :: nhx_n15_deposition
    dsl4jsb_Real2D_onChunk :: noy_n15_deposition
    dsl4jsb_Real2D_onChunk :: p_deposition
    dsl4jsb_Real2D_onChunk :: slow_sb_pool_accelerator_execute
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk :: f_n_demand
    dsl4jsb_Real2D_onChunk :: f_p_demand
    dsl4jsb_Real2D_onChunk :: net_biosphere_production
    dsl4jsb_Real2D_onChunk :: biological_n_fixation
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl
    dsl4jsb_Real3D_onChunk :: percolation_sl
    dsl4jsb_Real3D_onChunk :: qmax_org_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_po4_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_nh4_min_sl
    dsl4jsb_Real3D_onChunk :: qmax_po4_om_sl
    dsl4jsb_Real3D_onChunk :: sand_sl
    dsl4jsb_Real3D_onChunk :: silt_sl
    dsl4jsb_Real3D_onChunk :: clay_sl
    dsl4jsb_Real3D_onChunk :: volume_min_sl
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk :: total_flux_carbon
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen
    dsl4jsb_Real2D_onChunk :: total_flux_phosphorus
    dsl4jsb_Real2D_onChunk :: total_flux_carbon13
    dsl4jsb_Real2D_onChunk :: total_flux_carbon14
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen15
    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk :: leaching_nh4_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon
    dsl4jsb_Real2D_onChunk :: leaching_dom_nitrogen
    dsl4jsb_Real2D_onChunk :: leaching_dom_phosphorus
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon13
    dsl4jsb_Real2D_onChunk :: leaching_dom_carbon14
    dsl4jsb_Real2D_onChunk :: leaching_dom_nitrogen15
    dsl4jsb_Real2D_onChunk :: emission_noy
    dsl4jsb_Real2D_onChunk :: emission_noy_n15
    dsl4jsb_Real2D_onChunk :: emission_n2o
    dsl4jsb_Real2D_onChunk :: emission_n2o_n15
    dsl4jsb_Real2D_onChunk :: emission_n2
    dsl4jsb_Real2D_onChunk :: emission_n2_n15
    dsl4jsb_Real2D_onChunk :: ecosystem_total_n_loss
    dsl4jsb_Real2D_onChunk :: sb_pool_total_ag_litter_c
    dsl4jsb_Real2D_onChunk :: sb_pool_total_ag_litter_n
    dsl4jsb_Real2D_onChunk :: sb_pool_total_ag_litter_p
    dsl4jsb_Real2D_onChunk :: sb_pool_total_ag_litter_c13
    dsl4jsb_Real2D_onChunk :: sb_pool_total_ag_litter_c14
    dsl4jsb_Real2D_onChunk :: sb_pool_total_ag_litter_n15
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: myc_respiration
    dsl4jsb_Real3D_onChunk :: myc_respiration_c13
    dsl4jsb_Real3D_onChunk :: myc_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: nitrification_no3
    dsl4jsb_Real3D_onChunk :: nitrification_no3_n15
    dsl4jsb_Real3D_onChunk :: nitrification_noy
    dsl4jsb_Real3D_onChunk :: nitrification_noy_n15
    dsl4jsb_Real3D_onChunk :: nitrification_n2o
    dsl4jsb_Real3D_onChunk :: nitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_noy
    dsl4jsb_Real3D_onChunk :: denitrification_noy_n15
    dsl4jsb_Real3D_onChunk :: denitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_n2
    dsl4jsb_Real3D_onChunk :: denitrification_n2_n15
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation_n15
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
    dsl4jsb_Real3D_onChunk :: weathering_po4
    dsl4jsb_Real3D_onChunk :: occlusion_po4
    dsl4jsb_Real3D_onChunk :: slow_exchange_po4
    dsl4jsb_Real3D_onChunk :: fast_exchange_po4
    dsl4jsb_Real3D_onChunk :: fast_exchange_nh4
    dsl4jsb_Real3D_onChunk :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: gross_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: particle_fluxrate
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: transport_nh4_assoc
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk :: transport_po4_primary
    dsl4jsb_Real3D_onChunk :: transport_noy
    dsl4jsb_Real3D_onChunk :: transport_noy_n15
    dsl4jsb_Real3D_onChunk :: transport_n2o
    dsl4jsb_Real3D_onChunk :: transport_n2o_n15
    dsl4jsb_Real3D_onChunk :: transport_n2
    dsl4jsb_Real3D_onChunk :: transport_n2_n15
    dsl4jsb_Real3D_onChunk :: lateral_loss_nh4_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_nh4_n15_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_no3_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_no3_n15_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_po4_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_carbon_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_nitrogen_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_phosphorus_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_carbon13_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_carbon14_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_dom_nitrogen15_sl
    dsl4jsb_Real3D_onChunk :: nh4_solute
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: nh4_assoc
    dsl4jsb_Real3D_onChunk :: nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: no3_n15_solute
    dsl4jsb_Real3D_onChunk :: po4_solute
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: po4_primary
    dsl4jsb_Real3D_onChunk :: noy
    dsl4jsb_Real3D_onChunk :: noy_n15
    dsl4jsb_Real3D_onChunk :: n2o
    dsl4jsb_Real3D_onChunk :: n2o_n15
    dsl4jsb_Real3D_onChunk :: n2
    dsl4jsb_Real3D_onChunk :: n2_n15
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: residue_som_c_form_mavg_sl
    dsl4jsb_Real3D_onChunk :: residue_som_c14_form_mavg_sl
    dsl4jsb_Real3D_onChunk :: residue_assoc_som_c_form_mavg_sl
    dsl4jsb_Real3D_onChunk :: residue_assoc_som_c14_form_mavg_sl
    dsl4jsb_Real3D_onChunk :: assoc_dom_c_form_mavg_sl
    dsl4jsb_Real3D_onChunk :: assoc_dom_c14_form_mavg_sl
    dsl4jsb_Real3D_onChunk :: residue_som_c_loss_mavg_sl
    dsl4jsb_Real3D_onChunk :: residue_assoc_som_c_loss_mavg_sl
    dsl4jsb_Real3D_onChunk :: assoc_dom_c_loss_mavg_sl
    dsl4jsb_Real3D_onChunk :: qmax_org
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: qmax_fast_po4
    dsl4jsb_Real3D_onChunk :: qmax_slow_po4
    dsl4jsb_Real3D_onChunk :: km_fast_po4
    dsl4jsb_Real3D_onChunk :: km_slow_po4
    dsl4jsb_Real3D_onChunk :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onChunk :: qmax_nh4
    dsl4jsb_Real3D_onChunk :: km_adsorpt_nh4_sl
    dsl4jsb_Real3D_onChunk :: dom_cn
    dsl4jsb_Real3D_onChunk :: dom_cp
    dsl4jsb_Real3D_onChunk :: ph_sl
    dsl4jsb_Real3D_onChunk :: Qmax_AlFe_cor
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onChunk :: soil_litter_carbon_sl
    dsl4jsb_Real3D_onChunk :: sb_pool_total_c
    dsl4jsb_Real3D_onChunk :: sb_pool_total_n
    dsl4jsb_Real3D_onChunk :: sb_pool_total_p
    dsl4jsb_Real3D_onChunk :: sb_pool_total_c13
    dsl4jsb_Real3D_onChunk :: sb_pool_total_c14
    dsl4jsb_Real3D_onChunk :: sb_pool_total_n15
    dsl4jsb_Real3D_onChunk :: sb_pool_total_bg_soil_c
    dsl4jsb_Real3D_onChunk :: sb_pool_total_bg_soil_n
    dsl4jsb_Real3D_onChunk :: sb_pool_total_bg_soil_p
    dsl4jsb_Real3D_onChunk :: sb_pool_total_bg_soil_c13
    dsl4jsb_Real3D_onChunk :: sb_pool_total_bg_soil_c14
    dsl4jsb_Real3D_onChunk :: sb_pool_total_bg_soil_n15
    dsl4jsb_Real3D_onChunk :: sb_pool_woody_litter_c
    dsl4jsb_Real3D_onChunk :: total_soil_n
    dsl4jsb_Real3D_onChunk :: total_soil_inorg_n
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    dtime     = options%dtime
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    ! ---------------------------
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! Get process config
    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(A2L_)

    ! Allocate bgcm matrices
    nr_elements     = model%config%nr_of_used_elements
    nr_sb_parts     = get_number_of_sb_compartments()
    ALLOCATE(sb_pool_prev_timestep_hlp_mt(nr_sb_parts, nr_elements, nc, nsoil_sb))
    ALLOCATE(sb_particle_transport_hlp_mt(nr_sb_parts, nr_elements, nc, nsoil_sb))
    ! And get them from the store
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_EXUDATION_ID, veg_exudation_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FORMATION_ID, sb_formation_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_LOSS_ID, sb_loss_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_TRANSPORT_ID, sb_transport_mt)
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)
    dsl4jsb_Get_mt1L3D(SB_BGCM_DELTA_POOLS_ID, sb_delta_pools_mt)

    ! Allocate local allocatable variables (deallocate below at the end of this routine)
    ALLOCATE(p_deposition_sl(nc, nsoil_sb))
    ALLOCATE(f_myc_sl(nc, nsoil_sb))
    ALLOCATE(hlp1(nc, nsoil_sb))
    ALLOCATE(hlp2(nc, nsoil_sb))
    ALLOCATE(hlp3(nc, nsoil_sb))
    ALLOCATE(km_adsorpt_po4_act(nc, nsoil_sb))
    ALLOCATE(qmax_po4_mineral_act(nc, nsoil_sb))
    ALLOCATE(km_fast_po4_act(nc, nsoil_sb))
    ALLOCATE(qmax_fast_po4_act(nc, nsoil_sb))
    ALLOCATE(km_slow_po4_act(nc, nsoil_sb))
    ALLOCATE(qmax_slow_po4_act(nc, nsoil_sb))
    ALLOCATE(km_adsorpt_nh4_act(nc, nsoil_sb))
    ALLOCATE(qmax_nh4_mineral_act(nc, nsoil_sb))
    ALLOCATE(k_desorpt_po4_act(nc, nsoil_sb))
    ALLOCATE(k_adsorpt_po4_act(nc, nsoil_sb))
    ALLOCATE(rtm_sorption_act(nc, nsoil_sb))
    ALLOCATE(rmm_sorption_act(nc, nsoil_sb))
    ALLOCATE(rtm_desorption_act(nc, nsoil_sb))
    ALLOCATE(rmm_desorption_act(nc, nsoil_sb))
    ALLOCATE(partition_coef(nc, nsoil_sb))
    ALLOCATE(root_depth_avg(nc))
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_,      nhx_deposition)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      noy_deposition)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      nhx_n15_deposition)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      noy_n15_deposition)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      p_deposition)                ! in
    IF (model%config%flag_slow_sb_pool_spinup_accelerator) THEN
      dsl4jsb_Get_var2D_onChunk(A2L_,      slow_sb_pool_accelerator_execute)  !in
    ENDIF
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_,   f_n_demand)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_,   f_p_demand)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_,   net_biosphere_production)       ! out
    dsl4jsb_Get_var2D_onChunk(VEG_,   biological_n_fixation)          ! out
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_,   root_fraction_sl)                     ! in
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_,     num_sl_above_bedrock)        ! in
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_,     soil_depth_sl)               ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     bulk_dens_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     percolation_sl)              ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     qmax_org_min_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     qmax_po4_min_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     qmax_nh4_min_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     qmax_po4_om_sl)              ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     sand_sl)                     ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     silt_sl)                     ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     clay_sl)                     ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     volume_min_sl)               ! in
    ! SB_ 2D
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon)           ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_phosphorus)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon13)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon14)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen15)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_n15_solute)     ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_n15_solute)     ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_carbon)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_nitrogen)       ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_phosphorus)     ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_carbon13)       ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_carbon14)       ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_dom_nitrogen15)     ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy)                ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy_n15)            ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o)                ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o_n15)            ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2)                 ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2_n15)             ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       ecosystem_total_n_loss)      ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       sb_pool_total_ag_litter_c)   ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       sb_pool_total_ag_litter_n)   ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       sb_pool_total_ag_litter_p)   ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       sb_pool_total_ag_litter_c13) ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       sb_pool_total_ag_litter_c14) ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       sb_pool_total_ag_litter_n15) ! out
    ! SB_ 3D
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration_c13)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration_c14)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3_n15)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy_n15)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o_n15)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy_n15)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o_n15)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2)          ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2_n15)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4)          ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4_n15)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation)            ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation_n15)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_po4_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_n15_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_n15_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_n15_sl)! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_n15_sl)! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_sl)   ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_n15_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_po4_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       particle_fluxrate)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_assoc)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_assoc)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_fast)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_slow)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_occluded)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_primary)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy)               ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy_n15)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o)               ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o_n15)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2_n15)            ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_nh4_solute_sl)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_nh4_n15_solute_sl)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_no3_solute_sl)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_no3_n15_solute_sl)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_po4_solute_sl)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_dom_carbon_sl)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_dom_nitrogen_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_dom_phosphorus_sl)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_dom_carbon13_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_dom_carbon14_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_dom_nitrogen15_sl)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_assoc)                   ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_assoc)               ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded)                ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_primary)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       noy)                         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       noy_n15)                     ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o)                         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o_n15)                     ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2)                          ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2_n15)                      ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       weathering_po4)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       occlusion_po4)               ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       slow_exchange_po4)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       fast_exchange_po4)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       fast_exchange_nh4)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       biochem_mineralisation_po4)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       gross_mineralisation_po4)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)              ! in
    IF (model%config%flag_slow_sb_pool_spinup_accelerator) THEN
      dsl4jsb_Get_var3D_onChunk(SB_,       residue_som_c_form_mavg_sl)         ! in
      dsl4jsb_Get_var3D_onChunk(SB_,       residue_som_c14_form_mavg_sl)       ! in
      dsl4jsb_Get_var3D_onChunk(SB_,       residue_assoc_som_c_form_mavg_sl)   ! in
      dsl4jsb_Get_var3D_onChunk(SB_,       residue_assoc_som_c14_form_mavg_sl) ! in
      dsl4jsb_Get_var3D_onChunk(SB_,       assoc_dom_c_form_mavg_sl)           ! in
      dsl4jsb_Get_var3D_onChunk(SB_,       assoc_dom_c14_form_mavg_sl)         ! in
      dsl4jsb_Get_var3D_onChunk(SB_,       residue_som_c_loss_mavg_sl)         ! in
      dsl4jsb_Get_var3D_onChunk(SB_,       residue_assoc_som_c_loss_mavg_sl)   ! in
      dsl4jsb_Get_var3D_onChunk(SB_,       assoc_dom_c_loss_mavg_sl)           ! in
    END IF
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org)                    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_fast_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_slow_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       km_fast_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       km_slow_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       km_adsorpt_po4_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_nh4)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       km_adsorpt_nh4_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       dom_cn)                      ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       dom_cp)                      ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       ph_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       Qmax_AlFe_cor)
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_soil_carbon_sl)         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       soil_litter_carbon_sl)       ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_c)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_n)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_p)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_c13)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_c14)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_n15)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_bg_soil_c)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_bg_soil_n)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_bg_soil_p)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_bg_soil_c13)   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_bg_soil_c14)   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_total_bg_soil_n15)   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       sb_pool_woody_litter_c)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       total_soil_n)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       total_soil_inorg_n)          ! out
    !------------------------------------------------------------------------------------------------------ !

    !> 1.0 correct the plant uptake of nh4, no3 and po4 for soil model
    !!
    IF (model%config%qmodel_id==QSOIL) THEN
      DO ic = 1,nc
        DO is = 1,INT(num_sl_above_bedrock(ic))
          plant_uptake_nh4_sl(ic, is)     = plant_uptake_nh4_sl(ic, is)     * f_n_demand(ic)
          plant_uptake_nh4_n15_sl(ic, is) = plant_uptake_nh4_n15_sl(ic, is) * f_n_demand(ic)
          plant_uptake_no3_sl(ic, is)     = plant_uptake_no3_sl(ic, is)     * f_n_demand(ic)
          plant_uptake_no3_n15_sl(ic, is) = plant_uptake_no3_n15_sl(ic, is) * f_n_demand(ic)
          plant_uptake_po4_sl(ic, is)     = plant_uptake_po4_sl(ic, is)     * f_p_demand(ic)
        ENDDO
      ENDDO
    ENDIF

    !> 2.0 infrastructure for mass balance checks
    !!
    !! fluxes are calculated at the end of this routine (to include het_respiration)
    !!
    !! changes in pools per time step are calculated here AND at the end of this routine
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        sb_delta_pools_mt(ixC,ic, is) = SUM(sb_pool_mt(:, ixC, ic, is))
        sb_delta_pools_mt(ixN,ic, is) = SUM(sb_pool_mt(:, ixN, ic, is)) &
          & + nh4_solute(ic, is) + nh4_assoc(ic, is) + no3_solute(ic, is) + noy(ic, is) + n2o(ic, is) + n2(ic, is)
        sb_delta_pools_mt(ixP,ic, is) = SUM(sb_pool_mt(:, ixP, ic, is)) &
          & + po4_solute(ic, is) + po4_occluded(ic, is) + po4_assoc_slow(ic, is) + po4_assoc_fast(ic, is) + po4_primary(ic, is)
        sb_delta_pools_mt(ixC13,ic, is) = SUM(sb_pool_mt(:, ixC13, ic, is))
        ! for mass balance purposes, already reduce by radioactive decay, since this is not a reported flux
        sb_delta_pools_mt(ixC14,ic, is) = SUM(sb_pool_mt(:, ixC14, ic, is)) &
          & * (1._wp - lambda_C14 * dtime)
        sb_delta_pools_mt(ixN15,ic, is) = SUM(sb_pool_mt(:, ixN15, ic, is)) &
          & + nh4_n15_solute(ic, is) + nh4_n15_assoc(ic, is) + no3_n15_solute(ic, is) &
          & + noy_n15(ic, is) + n2o_n15(ic, is) + n2_n15(ic, is)
      ENDDO
    ENDDO

    !> 3.0 loss of C14 due to radioactive decay
    !!
    sb_loss_mt(:,ixC14,:,:) = sb_loss_mt(:,ixC14,:,:) + sb_pool_mt(:,ixC14,:,:) * lambda_C14 * dtime

    !> 4.0 update mycorrhiza with flux due to exudation from vegetation and N uptake from soil AND correct flux due to export to plants
    !!
    f_myc_sl(:,:)          = root_fraction_sl(:,:)

    f_myc_sl(:,:)          = root_fraction_sl(:,:)

    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        sb_formation_mt(ix_mycorrhiza, ixC, ic, is) = sb_formation_mt(ix_mycorrhiza, ixC, ic, is) &
          & + f_myc_sl(ic, is) * mycorrhiza_cue * veg_exudation_mt(ixC, ic) / soil_depth_sl(ic, is)

        sb_formation_mt(ix_mycorrhiza, ixN, ic, is) = sb_formation_mt(ix_mycorrhiza, ixN, ic, is) &
          & + f_myc_sl(ic, is) * veg_exudation_mt(ixN, ic) / soil_depth_sl(ic, is)                   &
          & + (mycorrhiza_uptake_nh4_sl(ic, is) + mycorrhiza_uptake_no3_sl(ic, is))                 &
          & * dtime / 1000000._wp

        sb_formation_mt(ix_mycorrhiza, ixP, ic, is) = sb_formation_mt(ix_mycorrhiza, ixP, ic, is) &
          & + f_myc_sl(ic, is) * veg_exudation_mt(ixP, ic) / soil_depth_sl(ic, is)                   &
          & + mycorrhiza_uptake_po4_sl(ic, is) * dtime / 1000000._wp

        sb_formation_mt(ix_mycorrhiza, ixC13, ic, is) = sb_formation_mt(ix_mycorrhiza, ixC13, ic, is) &
          & + f_myc_sl(ic, is) * mycorrhiza_cue * veg_exudation_mt(ixC13, ic) / soil_depth_sl(ic, is)

        sb_formation_mt(ix_mycorrhiza, ixC14, ic, is) = sb_formation_mt(ix_mycorrhiza, ixC14, ic, is) &
          & + f_myc_sl(ic, is) * mycorrhiza_cue * veg_exudation_mt(ixC14, ic) / soil_depth_sl(ic, is)

        sb_formation_mt(ix_mycorrhiza, ixN15, ic, is) = sb_formation_mt(ix_mycorrhiza, ixN15, ic, is) &
          & + f_myc_sl(ic, is) * veg_exudation_mt(ixN15, ic) / soil_depth_sl(ic, is)                     &
          & + (mycorrhiza_uptake_nh4_n15_sl(ic, is) + mycorrhiza_uptake_no3_n15_sl(ic, is))             &
          & * dtime / 1000000._wp

         het_respiration(ic, is)     = het_respiration(ic, is)     + f_myc_sl(ic, is) * (1._wp - mycorrhiza_cue) * &
          veg_exudation_mt(ixC,ic) / soil_depth_sl(ic, is) / dtime * 1000000._wp
         myc_respiration(ic, is)     = myc_respiration(ic, is)     + f_myc_sl(ic, is) * (1._wp - mycorrhiza_cue) * &
          veg_exudation_mt(ixC,ic) / soil_depth_sl(ic, is) / dtime * 1000000._wp
         het_respiration_c13(ic, is) = het_respiration_c13(ic, is) + f_myc_sl(ic, is) * (1._wp - mycorrhiza_cue) * &
          veg_exudation_mt(ixC13,ic) / soil_depth_sl(ic, is) / dtime * 1000000._wp
         myc_respiration_c13(ic, is) = myc_respiration_c13(ic, is) + f_myc_sl(ic, is) * (1._wp - mycorrhiza_cue) * &
          veg_exudation_mt(ixC13,ic) / soil_depth_sl(ic, is) / dtime * 1000000._wp
         het_respiration_c14(ic, is) = het_respiration_c14(ic, is) + f_myc_sl(ic, is) * (1._wp - mycorrhiza_cue) * &
          veg_exudation_mt(ixC14,ic) / soil_depth_sl(ic, is) / dtime * 1000000._wp
         myc_respiration_c14(ic, is) = myc_respiration_c14(ic, is) + f_myc_sl(ic, is) * (1._wp - mycorrhiza_cue) * &
          veg_exudation_mt(ixC14,ic) / soil_depth_sl(ic, is) / dtime * 1000000._wp
      ENDDO
    ENDDO

    sb_loss_mt(ix_mycorrhiza,:,:,:) = sb_loss_mt(ix_mycorrhiza,:,:,:) + sb_mycorrhiza_export_mt(:,:,:)

    !> 5.0 update of the organic pool: simple_1d
    !!
    !! update all components of the sb_pool applying the JSM equation for formation, loss and transport
    !!
    IF (TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_1d") THEN
      sb_pool_mt = sb_pool_mt + sb_formation_mt - sb_loss_mt + sb_transport_mt
    END IF

    !> 6.0 JSM:
    !!  a) update of the organic pool
    !!  b) update particle transport related diagnostics and state variables
    !!
    !! docu:
    !!  a) sb_last_timestep pool is used only for calc_particle_fluxrate() to calc (pseudo code):
    !!     particle_fluxrate = sb_pool - sb_last_timestep  (i.e., current timestep - previous timestep)
    !!  b) par_transport pool is used only as hlp to collect the output of calc_particle_transport()
    !!     within calc_particle_transport_wrapper() and apply these temporary values to sb_transport & sb_pool
    IF (TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "jsm") THEN
      !>  6.1 update all components of the sb_pool applying the JSM equation for formation, loss and transport
      !>
      sb_pool_prev_timestep_hlp_mt = sb_pool_mt
      sb_pool_mt                   = sb_pool_mt + sb_formation_mt - sb_loss_mt + sb_transport_mt
      !>  6.2 particle transport
      !>
      sb_particle_transport_hlp_mt = sb_transport_mt * 0._wp
      particle_fluxrate(:,:) = calc_particle_fluxrate(nc, nsoil_sb, ixC, soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), &
        &                                             sb_pool_mt, sb_pool_prev_timestep_hlp_mt)
      ! particle transport of organic and inorganic pools
      CALL calc_particle_transport_wrapper( &
        & nc, nsoil_sb, soil_depth_sl(:,:), &                                                 ! in
        & model%config%elements_index_map(:), &
        & model%config%is_element_used(:), &
        & particle_fluxrate(:,:), &
        & nh4_solute(:,:), nh4_assoc(:,:), no3_solute(:,:), po4_solute(:,:), &
        & po4_assoc_fast(:,:), po4_assoc_slow(:,:), po4_occluded(:,:), po4_primary(:,:), &
        & nh4_n15_solute(:,:), nh4_n15_assoc(:,:), no3_n15_solute(:,:), &
        & sb_pool_mt(:,:,:,:), &                                                              ! in
        & transport_nh4_solute(:,:), transport_nh4_assoc(:,:), &                              ! inout
        & transport_no3_solute(:,:), transport_po4_solute(:,:), &
        & transport_po4_assoc_fast(:,:), transport_po4_assoc_slow(:,:), &
        & transport_po4_occluded(:,:), transport_po4_primary(:,:), &
        & transport_nh4_n15_solute(:,:), transport_nh4_n15_assoc(:,:), &
        & transport_no3_n15_solute(:,:), &
        & sb_particle_transport_hlp_mt(:,:,:,:) )                                             ! inout
      ! ...
      sb_transport_mt(:,:,:,:) = sb_transport_mt(:, :, :, :) + sb_particle_transport_hlp_mt(:, :, :, :)
      sb_pool_mt(:,:,:,:)      = sb_pool_mt(:, :, :, :)      + sb_particle_transport_hlp_mt(:, :, :, :)
      dom_cn(:,:)              = sb_pool_mt(ix_dom, ixC, :, :) / sb_pool_mt(ix_dom, ixN, :, :)
      dom_cp(:,:)              = sb_pool_mt(ix_dom, ixC, :, :) / sb_pool_mt(ix_dom, ixP, :, :)
      ! unit conversion of inorganic transport fluxes to [mumol m3 s-1] for compatibility with the simple_1d soil model
      transport_nh4_solute(:,:)           = transport_nh4_solute(:,:)     * 1.e6_wp / dtime
      transport_nh4_assoc(:,:)            = transport_nh4_assoc(:,:)      * 1.e6_wp / dtime
      transport_no3_solute(:,:)           = transport_no3_solute(:,:)     * 1.e6_wp / dtime
      transport_po4_solute(:,:)           = transport_po4_solute(:,:)     * 1.e6_wp / dtime
      transport_po4_assoc_fast(:,:)       = transport_po4_assoc_fast(:,:) * 1.e6_wp / dtime
      transport_po4_assoc_slow(:,:)       = transport_po4_assoc_slow(:,:) * 1.e6_wp / dtime
      transport_po4_occluded(:,:)         = transport_po4_occluded(:,:)   * 1.e6_wp / dtime
      transport_po4_primary(:,:)          = transport_po4_primary(:,:)    * 1.e6_wp / dtime
      transport_nh4_n15_solute(:,:)       = transport_nh4_n15_solute(:,:) * 1.e6_wp / dtime
      transport_nh4_n15_assoc(:,:)        = transport_nh4_n15_assoc(:,:)  * 1.e6_wp / dtime
      transport_no3_n15_solute(:,:)       = transport_no3_n15_solute(:,:) * 1.e6_wp / dtime
      lateral_loss_dom_carbon_sl(:,:)     = lateral_loss_dom_carbon_sl(:,:)     * 1.e6_wp / dtime
      lateral_loss_dom_nitrogen_sl(:,:)   = lateral_loss_dom_nitrogen_sl(:,:)   * 1.e6_wp / dtime
      lateral_loss_dom_phosphorus_sl(:,:) = lateral_loss_dom_phosphorus_sl(:,:) * 1.e6_wp / dtime
      lateral_loss_dom_carbon13_sl(:,:)   = lateral_loss_dom_carbon13_sl(:,:)   * 1.e6_wp / dtime
      lateral_loss_dom_carbon14_sl(:,:)   = lateral_loss_dom_carbon14_sl(:,:)   * 1.e6_wp / dtime
      lateral_loss_dom_nitrogen15_sl(:,:) = lateral_loss_dom_nitrogen15_sl(:,:) * 1.e6_wp / dtime
      lateral_loss_nh4_solute_sl(:,:)     = lateral_loss_nh4_solute_sl(:,:)     * 1.e6_wp / dtime
      lateral_loss_no3_solute_sl(:,:)     = lateral_loss_no3_solute_sl(:,:)     * 1.e6_wp / dtime
      lateral_loss_po4_solute_sl(:,:)     = lateral_loss_po4_solute_sl(:,:)     * 1.e6_wp / dtime
      lateral_loss_nh4_n15_solute_sl(:,:) = lateral_loss_nh4_n15_solute_sl(:,:) * 1.e6_wp / dtime
      lateral_loss_no3_n15_solute_sl(:,:) = lateral_loss_no3_n15_solute_sl(:,:) * 1.e6_wp / dtime
      leaching_nh4_solute(:)              = leaching_nh4_solute(:)      * 1.e6_wp / dtime
      leaching_nh4_n15_solute(:)          = leaching_nh4_n15_solute(:)  * 1.e6_wp / dtime
      leaching_no3_solute(:)              = leaching_no3_solute(:)      * 1.e6_wp / dtime
      leaching_no3_n15_solute(:)          = leaching_no3_n15_solute(:)  * 1.e6_wp / dtime
      leaching_po4_solute(:)              = leaching_po4_solute(:)      * 1.e6_wp / dtime
      leaching_dom_carbon(:)              = leaching_dom_carbon(:)      * 1.e6_wp / dtime
      leaching_dom_nitrogen(:)            = leaching_dom_nitrogen(:)    * 1.e6_wp / dtime
      leaching_dom_phosphorus(:)          = leaching_dom_phosphorus(:)  * 1.e6_wp / dtime
      leaching_dom_carbon13(:)            = leaching_dom_carbon13(:)    * 1.e6_wp / dtime
      leaching_dom_carbon14(:)            = leaching_dom_carbon14(:)    * 1.e6_wp / dtime
      leaching_dom_nitrogen15(:)          = leaching_dom_nitrogen15(:)  * 1.e6_wp / dtime
    ENDIF

    !>7.0 update of inorganic pools
    !>
    !>  7.1 update of N pools
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        no3_solute(ic, is)     = no3_solute(ic, is)     + &
                              (net_mineralisation_no3(ic, is) + nitrification_no3(ic, is) + transport_no3_solute(ic, is) - &
                               lateral_loss_no3_solute_sl(ic, is) - &
                               plant_uptake_no3_sl(ic, is) - mycorrhiza_uptake_no3_sl(ic, is) - &
                               denitrification_noy(ic, is) - denitrification_n2o(ic, is) - denitrification_n2(ic, is)) &
                               * dtime / 1.e6_wp

        no3_n15_solute(ic, is) = no3_n15_solute(ic, is) + &
          (net_mineralisation_no3_n15(ic, is) + nitrification_no3_n15(ic, is) + transport_no3_n15_solute(ic, is) - &
          lateral_loss_no3_n15_solute_sl(ic, is) - &
          plant_uptake_no3_n15_sl(ic, is) - mycorrhiza_uptake_no3_n15_sl(ic, is) - &
          denitrification_noy_n15(ic, is) - denitrification_n2o_n15(ic, is) - denitrification_n2_n15(ic, is)) * dtime/1.e6_wp

        noy(ic, is)     = noy(ic, is)     + (nitrification_noy(ic, is) + denitrification_noy(ic, is) + &
                                                     transport_noy(ic, is)) * dtime/1.e6_wp
        noy_n15(ic, is) = noy_n15(ic, is) + (nitrification_noy_n15(ic, is) + denitrification_noy_n15(ic, is) + &
                                                     transport_noy_n15(ic, is)) * dtime/1.e6_wp
        n2o(ic, is)     = n2o(ic, is)     + (nitrification_n2o(ic, is) + denitrification_n2o(ic, is) + &
                                                     transport_n2o(ic, is)) * dtime/1.e6_wp
        n2o_n15(ic, is) = n2o_n15(ic, is) + (nitrification_n2o_n15(ic, is) + denitrification_n2o_n15(ic, is) + &
                                                     transport_n2o_n15(ic, is)) * dtime/1.e6_wp
        n2(ic, is)      = n2(ic, is)      + (denitrification_n2(ic, is) + transport_n2(ic, is)) * dtime/1.e6_wp
        n2_n15(ic, is)  = n2_n15(ic, is)  + (denitrification_n2_n15(ic, is) + transport_n2_n15(ic, is)) * dtime/1.e6_wp
      ENDDO
    ENDDO

    !> 8.0 update diagnostics and other state variables depending on organic matter content
    !! OC
    CALL calc_bulk_soil_carbon( &
      & nc, &                                        ! in
      & nsoil_sb, &
      & num_sl_above_bedrock(:), &
      & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
      & sb_pool_mt(:,:,:,:), &                       ! in
      & bulk_soil_carbon_sl, &                       ! inout
      & soil_litter_carbon_sl, &                     ! inout
      & volume_min_sl)                               ! inout
    bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:), soil_litter_carbon_sl(:,:), bulk_dens_sl(:,:))
    CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), &
      &                                    0._wp, qmax_org_min_sl(:,:), qmax_org(:,:))
    !! po4
    ! Correct the first layer Qmax_AlFe_cor with the woody litter
    Qmax_AlFe_cor(:,1) = Qmax_AlFe_cor(:,1) * 1._wp / ((1._wp - volume_min_sl(:,1))                                      &
      &                 * (bulk_soil_carbon_sl(:,1) + soil_litter_carbon_sl(:,1)  + sb_pool_mt(ix_woody_litter,ixC,:,1)) &
      &                 / (bulk_soil_carbon_sl(:,1) + soil_litter_carbon_sl(:,1)) + volume_min_sl(:,1))
    DO ic = 1,nc
      DO isoil = 1,INT(num_sl_above_bedrock(ic))
        CALL calc_Psorption_parameter( &
          & dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
          & clay_sl(ic, isoil), &
          & silt_sl(ic, isoil), &
          & rmm_sorption(ic, isoil) * 1000._wp, &
          & ph_sl(ic, isoil), &
          & bulk_soil_carbon_sl(ic, isoil), &
          & soil_litter_carbon_sl(ic, isoil), &
          & sb_pool_mt(ix_woody_litter,ixC,ic, isoil), &
          & bulk_dens_sl(ic, isoil), &
          & volume_min_sl(ic, isoil), &
          & po4_solute(ic, isoil), &
          & qmax_po4_min_sl(ic, isoil), &
          & qmax_po4_om_sl(ic, isoil), &
          & qmax_po4(ic, isoil), &
          & km_adsorpt_po4_sl(ic, isoil), &
          & qmax_fast_po4(ic, isoil), &
          & qmax_slow_po4(ic, isoil), &
          & km_fast_po4(ic, isoil), &
          & km_slow_po4(ic, isoil), &
          & Qmax_AlFe_cor(ic, isoil))
      ENDDO
    ENDDO
    !! nh4
    CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), &
      &                                    0.0_wp, qmax_nh4_min_sl(:,:), qmax_nh4(:,:))
    qmax_nh4(:,:) = qmax_nh4(:,:) * rmm_sorption(:,:)
    CALL calc_qmax_bulk_density_correction(bulk_soil_carbon_sl(:,:), volume_min_sl(:,:), bulk_dens_sl(:,:), &
      &                                    km_adsorpt_OM_nh4, km_adsorpt_mineral_nh4, km_adsorpt_nh4_sl(:,:))

    !> 9.0 atmospheric deposition
    !!
    no3_solute(:,1)      = no3_solute(:,1)     + noy_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime
    no3_n15_solute(:,1)  = no3_n15_solute(:,1) + noy_n15_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime

    !>10.0 update of mineral NH4 & P pools
    !>
    p_deposition_sl(:,:) = 0._wp
    p_deposition_sl(:,1) = p_deposition(:) / soil_depth_sl(:,1) ! micro-mol P / m3 / s
    !> 10.1 init vertical variables
    !>
    k_adsorpt_po4_act(:,:)    = rtm_sorption(:,:)   * rmm_sorption(:,:)   * k_adsorpt_po4
    k_desorpt_po4_act(:,:)    = rtm_desorption(:,:) * rmm_desorption(:,:) * k_desorpt_po4
    km_adsorpt_po4_act(:,:)   = km_adsorpt_po4_sl(:,:)
    qmax_po4_mineral_act(:,:) = qmax_po4(:,:)
    km_fast_po4_act(:,:)      = km_fast_po4(:,:)
    qmax_fast_po4_act(:,:)    = qmax_fast_po4(:,:)
    km_slow_po4_act(:,:)      = km_slow_po4(:,:)
    qmax_slow_po4_act(:,:)    = qmax_slow_po4(:,:)
    km_adsorpt_nh4_act(:,:)   = km_adsorpt_nh4_sl(:,:)
    qmax_nh4_mineral_act(:,:) = qmax_nh4(:,:)
    !> 10.2 update NH4 pools if not prescribed
    !>
    IF(.NOT. dsl4jsb_Config(SB_)%flag_sb_prescribe_nh4) THEN
      hlp1(:,:) = (net_mineralisation_nh4(:,:) + transport_nh4_solute(:,:) - &
                   plant_uptake_nh4_sl(:,:) - mycorrhiza_uptake_nh4_sl(:,:) - &
                   nitrification_no3(:,:) - nitrification_noy(:,:) - nitrification_n2o(:,:) - &
                   volatilisation_nh4(:,:) + transport_nh4_assoc(:,:) - &
                   lateral_loss_nh4_solute_sl(:,:) ) * dtime / 1.e6_wp
      hlp1(:,1) = hlp1(:,1) + nhx_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime
      hlp2(:,:) = -transport_nh4_assoc(:,:) * dtime / 1.e6_wp

      CALL calc_langmuir_kinetics(dtime, hlp1(:,:), nh4_solute(:,:),nh4_assoc(:,:), &
                                                   km_adsorpt_nh4_act(:,:), qmax_nh4_mineral_act(:,:), &
                                                   partition_coef(:,:))

      nh4_solute(:,:)        = nh4_solute(:,:) + partition_coef(:,:) * hlp1(:,:)
      nh4_assoc(:,:)         = nh4_assoc(:,:)  + (1._wp - partition_coef(:,:)) * hlp1(:,:)

      fast_exchange_nh4(:,:) = (partition_coef(:,:) * hlp1(:,:) + hlp2(:,:)) * 1.e6_wp / dtime
      hlp2(:,:) = (net_mineralisation_nh4_n15(:,:) + transport_nh4_n15_solute(:,:) - &
                   lateral_loss_nh4_n15_solute_sl(:,:) - &
                   plant_uptake_nh4_n15_sl(:,:) - mycorrhiza_uptake_nh4_n15_sl(:,:) - &
                   nitrification_no3_n15(:,:) - nitrification_noy_n15(:,:) - nitrification_n2o_n15(:,:) - &
                   volatilisation_nh4_n15(:,:) + transport_nh4_n15_assoc(:,:)) * dtime / 1.e6_wp
      hlp2(:,1) = hlp2(:,1) + nhx_n15_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime

      WHERE ((nh4_n15_solute(:,:) + nh4_n15_assoc(:,:)) < -hlp2(:,:))
        hlp2(:,:) = -(nh4_n15_solute(:,:) + nh4_n15_assoc(:,:))
      ENDWHERE
      WHERE (nh4_solute(:,:) + nh4_assoc(:,:) > eps8)
        partition_coef(:,:) = nh4_solute(:,:) / (nh4_solute(:,:) + nh4_assoc(:,:))
      ENDWHERE

      nh4_n15_solute(:,:) = partition_coef(:,:) * (nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + hlp2(:,:))
      nh4_n15_assoc(:,:)  = (1._wp - partition_coef(:,:)) * (nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + hlp2(:,:))
    ENDIF

    !> 10.3 update mineral P pools if not prescibed
    !!
    IF(.NOT. dsl4jsb_Config(SB_)%flag_sb_prescribe_po4) THEN
      CALL calc_p_inorg_dynamics(nc, nsoil_sb, dtime, &
                                       TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
                                       TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), &
                                       dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
                                       po4_solute(:,:), &
                                       po4_assoc_fast(:,:), &
                                       po4_assoc_slow(:,:), &
                                       qmax_po4_mineral_act(:,:), &
                                       km_adsorpt_po4_act(:,:), &
                                       qmax_fast_po4_act(:,:), &
                                       km_fast_po4_act(:,:), &
                                       qmax_slow_po4_act(:,:), &
                                       km_slow_po4_act(:,:), &
                                       k_adsorpt_po4_act(:,:), &
                                       k_desorpt_po4_act(:,:), &
                                       p_deposition_sl(:,:), &
                                       biochem_mineralisation_po4(:,:), &
                                       mycorrhiza_uptake_po4_sl(:,:), &
                                       plant_uptake_po4_sl(:,:), &
                                       microbial_uptake_po4_sl(:,:), &
                                       gross_mineralisation_po4(:,:), &
                                       net_mineralisation_po4(:,:), &
                                       transport_po4_solute(:,:), &
                                       lateral_loss_po4_solute_sl(:,:), &
                                       transport_po4_assoc_fast(:,:), &
                                       transport_po4_assoc_slow(:,:), &
                                       weathering_po4(:,:), &
                                       occlusion_po4(:,:), &
                                       slow_exchange_po4(:,:), &
                                       fast_exchange_po4(:,:))

      !! not update the primary and the occlusion pool during transient-forcing spinup (i.e., before the actual transient simulation)
      !!  -> not balancing the delta pool P and total flux P as we keep the primary P pool constant
      !!     despite a non-zero weathering flux from primary to solute P
      !!  -> this results in the "mass balance" not closed during spinup for both simple soil-model and jsm
      !!  -> To close the mass balance check, ONLY DURING SPIN-UP, we would need to add the weathering flux to the pool change
      !!     (as positive pool change).
      !! we do this because:
      !!   1) primary P will be depleted if the spin-up time is very long
      !!   2) we use the "soil P map" with the P status of about the year 1850 as model input
      !!      we may ensure that the soil P is not very different from this value after spinup, i.e., the transient simulation (usually starting 1900)
      IF (model%config%l_transient_spinup) THEN
        po4_occluded(:,:)   = po4_occluded(:,:)   +                         transport_po4_occluded(:,:)  * dtime / 1.e6_wp
        po4_primary(:,:)    = po4_primary(:,:)    +                         transport_po4_primary(:,:)   * dtime / 1.e6_wp
      ELSE
        po4_occluded(:,:)   = po4_occluded(:,:)   + (occlusion_po4(:,:)   + transport_po4_occluded(:,:)) * dtime / 1.e6_wp
        po4_primary(:,:)    = po4_primary(:,:)    + (-weathering_po4(:,:) + transport_po4_primary(:,:))  * dtime / 1.e6_wp
      ENDIF
    ENDIF

    !>11.0 Ensure that mineral nutrients are constant when prescibed
    !>

    !> 11.1 NH4
    !>
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_nh4)THEN
      hlp1(:,:) = nh4_solute_prescribe - nh4_solute(:,:)
      nh4_solute(:,:)     = nh4_solute(:,:) + hlp1(:,:)
      nh4_n15_solute(:,:) = nh4_n15_solute(:,:) + hlp1(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    ENDIF

    !> 11.2 NO3
    !>
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_no3)THEN
      hlp1(:,:) = no3_solute_prescribe - no3_solute(:,:)
      no3_solute(:,:)     = no3_solute(:,:) + hlp1(:,:)
      no3_n15_solute(:,:) = no3_n15_solute(:,:) + hlp1(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    ENDIF

    !> 11.3 PO4
    !>
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_po4)THEN
      po4_solute(:,:) = po4_solute_prescribe
    ENDIF


    !> 12.0 infrastructure for mass balance checks (see also section 0.9)
    !!
    !! pools: calc the difference to what was calculated at the beginning of this routine
    sb_delta_pools_mt(ixC,:,:) = SUM(sb_pool_mt(:,ixC,:,:), DIM=1) - sb_delta_pools_mt(ixC,:,:)
    sb_delta_pools_mt(ixN,:,:) = (SUM(sb_pool_mt(:,ixN,:,:), DIM=1)                           &
      & + nh4_solute(:,:) + nh4_assoc(:,:) + no3_solute(:,:) + noy(:,:) + n2o(:,:) + n2(:,:)) &
      & - sb_delta_pools_mt(ixN,:,:)
    sb_delta_pools_mt(ixP,:,:) = (SUM(sb_pool_mt(:,ixP,:,:), DIM=1)                                           &
      & + po4_solute(:,:) + po4_occluded(:,:) + po4_assoc_slow(:,:) + po4_assoc_fast(:,:) + po4_primary(:,:)) &
      & - sb_delta_pools_mt(ixP,:,:)
    sb_delta_pools_mt(ixC13,:,:) = SUM(sb_pool_mt(:,ixC13,:,:), DIM=1) - sb_delta_pools_mt(ixC13,:,:)
    sb_delta_pools_mt(ixC14,:,:) = SUM(sb_pool_mt(:,ixC14,:,:), DIM=1) - sb_delta_pools_mt(ixC14,:,:)

    sb_delta_pools_mt(ixN15,:,:) = (SUM(sb_pool_mt(:,ixN15,:,:), DIM=1)                                               &
      & + nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + no3_n15_solute(:,:) + noy_n15(:,:) + n2o_n15(:,:) + n2_n15(:,:)) &
      & - sb_delta_pools_mt(ixN15,:,:)

    !! fluxes
    total_flux_carbon(:) = SUM(veg_litterfall_mt(:,ixC,:), DIM=1) + veg_exudation_mt(ixC,:)               &
      & - SUM(sb_mycorrhiza_export_mt(ixC,:,:) * soil_depth_sl(:,:), DIM=2)                   &
      & - SUM(lateral_loss_dom_carbon_sl(:,:) * soil_depth_sl(:,:), DIM=2) * dtime / 1.e6_wp  &
      & - (leaching_dom_carbon(:) + SUM(het_respiration(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime / 1.e6_wp   ! DIM=2: sum only over dimension 2, i.e., soil layers

    total_flux_nitrogen(:) = SUM(veg_litterfall_mt(:,ixN,:), DIM=1)                            &
      & + SUM(asymb_n_fixation(:,:) * soil_depth_sl(:,:), DIM=2) *dtime/1.e6_wp                &
      & + (nhx_deposition(:) + noy_deposition(:)) * dtime/1.e6_wp                              &
      & + veg_exudation_mt(ixN, :)                                                             &
      & - SUM(sb_mycorrhiza_export_mt(ixN,:,:) * soil_depth_sl(:,:), DIM=2)                    &
      & - SUM(lateral_loss_dom_nitrogen_sl(:,:) * soil_depth_sl(:,:), DIM=2) * dtime / 1.e6_wp &
      & - SUM((lateral_loss_nh4_solute_sl(:,:) + lateral_loss_no3_solute_sl(:,:))              &
      &        * soil_depth_sl(:,:), DIM=2) * dtime / 1.e6_wp                                  &
      & - (leaching_nh4_solute(:) + leaching_no3_solute(:) + leaching_dom_nitrogen(:)          &
      & + emission_noy(:) + emission_n2o(:) + emission_n2(:)                                   &
      & + SUM(plant_uptake_nh4_sl(:,:) * soil_depth_sl(:,:), DIM=2)                            &
      & + SUM(plant_uptake_no3_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp

    total_flux_phosphorus(:) = SUM(veg_litterfall_mt(:,ixP,:), DIM=1) + p_deposition(:) * dtime/1.e6_wp &
      & + veg_exudation_mt(ixP,:)                                                                       &
      & - SUM(sb_mycorrhiza_export_mt(ixP,:,:) * soil_depth_sl(:,:), DIM=2)                             &
      & - SUM(lateral_loss_dom_phosphorus_sl(:,:) * soil_depth_sl(:,:), DIM=2) * dtime / 1.e6_wp        &
      & - SUM(lateral_loss_po4_solute_sl(:,:) * soil_depth_sl(:,:), DIM=2) * dtime / 1.e6_wp            &
      & - (leaching_po4_solute(:) + leaching_dom_phosphorus(:)                                          &
      & + SUM(plant_uptake_po4_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp

    total_flux_carbon13(:) = SUM(veg_litterfall_mt(:,ixC13,:), DIM=1) + veg_exudation_mt(ixC13,:)   &
      & - SUM(sb_mycorrhiza_export_mt(ixC13,:,:) * soil_depth_sl(:,:), DIM=2)                       &
      & - SUM(lateral_loss_dom_carbon13_sl(:,:) * soil_depth_sl(:,:), DIM=2) * dtime / 1.e6_wp      &
      & - (leaching_dom_carbon13(:) + SUM(het_respiration_c13(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime / 1.e6_wp
    total_flux_carbon14(:) = SUM(veg_litterfall_mt(:,ixC14,:), DIM=1) + veg_exudation_mt(ixC14,:)   &
      & - SUM(sb_mycorrhiza_export_mt(ixC14,:,:) * soil_depth_sl(:,:), DIM=2)                       &
      & - SUM(lateral_loss_dom_carbon14_sl(:,:) * soil_depth_sl(:,:), DIM=2) * dtime / 1.e6_wp      &
      & - (leaching_dom_carbon14(:) + SUM(het_respiration_c14(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime / 1.e6_wp

    total_flux_nitrogen15(:) = SUM(veg_litterfall_mt(:,ixN15,:), DIM=1)                          &
      & + SUM(asymb_n_fixation_n15(:,:) * soil_depth_sl(:,:), DIM=2) *dtime/1.e6_wp              &
      & + (nhx_deposition(:) + noy_deposition(:))                                                &
      &    / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) *dtime/1.e6_wp                  &
      & + veg_exudation_mt(ixN15, :)                                                             &
      & - SUM(sb_mycorrhiza_export_mt(ixN15,:,:) * soil_depth_sl(:,:), DIM=2)                    &
      & - SUM(lateral_loss_dom_nitrogen15_sl(:,:) * soil_depth_sl(:,:), DIM=2) * dtime / 1.e6_wp &
      & - SUM((lateral_loss_nh4_n15_solute_sl(:,:) + lateral_loss_no3_n15_solute_sl(:,:))        &
      &        * soil_depth_sl(:,:), DIM=2) * dtime / 1.e6_wp                                    &
      & - (leaching_nh4_n15_solute(:) + leaching_no3_n15_solute(:) + leaching_dom_nitrogen15(:)  &
      & + emission_noy_n15(:) + emission_n2o_n15(:) + emission_n2_n15(:)                         &
      & + SUM(plant_uptake_nh4_n15_sl(:,:) * soil_depth_sl(:,:), DIM=2)                          &
      & + SUM(plant_uptake_no3_n15_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp

    !> 13.0 bgc_material diagnostics - sum up sb_pool components
    !>      sum up: all, soil organic material, and litter for each particular element
    !>      or just pass the values from one/multiple variable to the diagnostics variable
    !>
    sb_pool_total_c(:,:)          = SUM(sb_pool_mt(:, ixC, :, :), DIM=1)
    sb_pool_total_n(:,:)          = SUM(sb_pool_mt(:, ixN, :, :), DIM=1)
    sb_pool_total_p(:,:)          = SUM(sb_pool_mt(:, ixP, :, :), DIM=1)
    sb_pool_total_c13(:,:)        = SUM(sb_pool_mt(:, ixC13, :, :), DIM=1)
    sb_pool_total_c14(:,:)        = SUM(sb_pool_mt(:, ixC14, :, :), DIM=1)
    sb_pool_total_n15(:,:)        = SUM(sb_pool_mt(:, ixN15, :, :), DIM=1)

    sb_pool_woody_litter_c(:,:)   = sb_pool_mt(ix_woody_litter, ixC, :, :)

    CALL calculate_sb_pool_bg_soil_for_element(sb_pool_mt, nsoil_sb, ixC, sb_pool_total_bg_soil_c)
    CALL calculate_sb_pool_bg_soil_for_element(sb_pool_mt, nsoil_sb, ixN, sb_pool_total_bg_soil_n)
    CALL calculate_sb_pool_bg_soil_for_element(sb_pool_mt, nsoil_sb, ixP, sb_pool_total_bg_soil_p)
    CALL calculate_sb_pool_bg_soil_for_element(sb_pool_mt, nsoil_sb, ixC13, sb_pool_total_bg_soil_c13)
    CALL calculate_sb_pool_bg_soil_for_element(sb_pool_mt, nsoil_sb, ixC14, sb_pool_total_bg_soil_c14)
    CALL calculate_sb_pool_bg_soil_for_element(sb_pool_mt, nsoil_sb, ixN15, sb_pool_total_bg_soil_n15)

    CALL calculate_sb_pool_ag_litter_for_element(sb_pool_mt, soil_depth_sl, ixC, sb_pool_total_ag_litter_c)
    CALL calculate_sb_pool_ag_litter_for_element(sb_pool_mt, soil_depth_sl, ixN, sb_pool_total_ag_litter_n)
    CALL calculate_sb_pool_ag_litter_for_element(sb_pool_mt, soil_depth_sl, ixP, sb_pool_total_ag_litter_p)
    CALL calculate_sb_pool_ag_litter_for_element(sb_pool_mt, soil_depth_sl, ixC13, sb_pool_total_ag_litter_c13)
    CALL calculate_sb_pool_ag_litter_for_element(sb_pool_mt, soil_depth_sl, ixC14, sb_pool_total_ag_litter_c14)
    CALL calculate_sb_pool_ag_litter_for_element(sb_pool_mt, soil_depth_sl, ixN15, sb_pool_total_ag_litter_n15)

    !> 14.0 biosphere-level diagnostics (SB_ process)
    !>

    !> 14.1 SB_ process
    !>
    ecosystem_total_n_loss(:) = &
      & emission_noy(:) + emission_n2o(:) + emission_n2(:) &
      & + leaching_dom_nitrogen(:) + leaching_nh4_solute(:) + leaching_no3_solute(:) &
      & + SUM(lateral_loss_dom_nitrogen_sl(:,:) * soil_depth_sl(:,:), DIM=2) &
      & + SUM(lateral_loss_nh4_solute_sl(:,:) * soil_depth_sl(:,:), DIM=2) &
      & + SUM(lateral_loss_no3_solute_sl(:,:) * soil_depth_sl(:,:), DIM=2)
    total_soil_n(:,:)         = sb_pool_total_bg_soil_n(:,:) + nh4_solute(:,:) + nh4_assoc(:,:) + no3_solute(:,:)
    total_soil_inorg_n(:,:)   = nh4_solute(:,:) + nh4_assoc(:,:) + no3_solute(:,:)

    !> 14.2 across multiple processes
    !>
    !>   see also 'update_veg_pools()'
    net_biosphere_production(:) = net_biosphere_production(:) &
      &          - SUM((het_respiration(:,:) + lateral_loss_dom_carbon_sl(:,:)) * soil_depth_sl(:,:), DIM = 2)
    ! biological N fixation (VEG_ and SB_ processes)
    biological_n_fixation(:)    = biological_n_fixation &
      &                                     + SUM(asymb_n_fixation(:,:) * soil_depth_sl(:,:), DIM = 2)

    !> 15.0 spin-up accelerator
    IF (model%config%flag_slow_sb_pool_spinup_accelerator) THEN
      SELECT CASE(model%config%qmodel_id)
      ! only execute spin-up accelerator when running in LAND or SOIL mode
      CASE(QLAND, QSOIL)
        DO ic = 1,nc
          ! Determine if running a spin-up with acceleration
          IF (slow_sb_pool_accelerator_execute(ic) == 1.0_wp) THEN
            IF (debug_on() .AND. iblk == 1 .AND. ic == 1) &
              & CALL message(TRIM(routine), 'Accelerate spin up for '//TRIM(tile%name)//' ...')

            ! fasten the spinup for residual pool
            CALL calc_accelerated_slow_soil_pool_spinup(nsoil_sb, dtime, soil_depth_sl(ic,:),&
                                                      & model%config%slow_sb_pool_spinup_accelerator_length, &
                                                      & residue_som_c_form_mavg_sl(ic,:), &
                                                      & residue_som_c14_form_mavg_sl(ic,:), &
                                                      & residue_som_c_loss_mavg_sl(ic,:), &
                                                      & sb_pool_mt(ix_residue,:,ic,:))
            IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "jsm") THEN
              ! fasten the spinup for the associated residual pool
              CALL calc_accelerated_slow_soil_pool_spinup(nsoil_sb, dtime, soil_depth_sl(ic,:), &
                                                        & model%config%slow_sb_pool_spinup_accelerator_length, &
                                                        & residue_assoc_som_c_form_mavg_sl(ic,:), &
                                                        & residue_assoc_som_c14_form_mavg_sl(ic,:), &
                                                        & residue_assoc_som_c_loss_mavg_sl(ic,:), &
                                                        & sb_pool_mt(ix_residue_assoc,:,ic,:))
              ! fasten the spinup for the associated dom pool
              CALL calc_accelerated_slow_soil_pool_spinup(nsoil_sb, dtime, soil_depth_sl(ic,:), &
                                                        & model%config%slow_sb_pool_spinup_accelerator_length, &
                                                        & assoc_dom_c_form_mavg_sl(ic,:), &
                                                        & assoc_dom_c14_form_mavg_sl(ic,:), &
                                                        & assoc_dom_c_loss_mavg_sl(ic,:), &
                                                        & sb_pool_mt(ix_dom_assoc,:,ic,:))
            END IF ! simple soil or jsm
          END IF
        END DO
      END SELECT
    END IF
    !------------------------------------------------------------------------------------------------------ !
    DEALLOCATE(sb_particle_transport_hlp_mt, sb_pool_prev_timestep_hlp_mt)

    ! Deallocate local allocatable variables
    DEALLOCATE(p_deposition_sl)
    DEALLOCATE(f_myc_sl)
    DEALLOCATE(hlp1)
    DEALLOCATE(hlp2)
    DEALLOCATE(hlp3)
    DEALLOCATE(km_adsorpt_po4_act)
    DEALLOCATE(qmax_po4_mineral_act)
    DEALLOCATE(km_fast_po4_act)
    DEALLOCATE(qmax_fast_po4_act)
    DEALLOCATE(km_slow_po4_act)
    DEALLOCATE(qmax_slow_po4_act)
    DEALLOCATE(k_desorpt_po4_act)
    DEALLOCATE(k_adsorpt_po4_act)
    DEALLOCATE(rtm_sorption_act)
    DEALLOCATE(rmm_sorption_act)
    DEALLOCATE(rtm_desorption_act)
    DEALLOCATE(rmm_desorption_act)
    DEALLOCATE(partition_coef)
    DEALLOCATE(root_depth_avg)

  END SUBROUTINE update_sb_pools


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate fast_exchange_po4 based on Langmuir equilibrium (Yang et al. 2014)
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_p_inorg_dynamics(nc, nsoil_sb, &
                                   dtime, &
                                   sb_model_scheme, &
                                   sb_adsorp_scheme, &
                                   flag_sb_double_langmuir, &
                                   po4_solute, &
                                   po4_assoc_fast, &
                                   po4_assoc_slow, &
                                   qmax_po4_mineral_act, &
                                   km_adsorpt_po4_act, &
                                   qmax_fast_po4_act, &
                                   km_fast_po4_act, &
                                   qmax_slow_po4_act, &
                                   km_slow_po4_act, &
                                   k_adsorpt_po4_act, &
                                   k_desorpt_po4_act, &
                                   p_deposition_sl, &
                                   biochem_mineralisation_po4, &
                                   mycorrhiza_uptake_po4_sl, &
                                   plant_uptake_po4_sl, &
                                   microbial_uptake_po4_sl, &
                                   gross_mineralisation_po4, &
                                   net_mineralisation_po4, &
                                   transport_po4_solute, &
                                   lateral_loss_po4_solute_sl, &
                                   transport_po4_assoc_fast, &
                                   transport_po4_assoc_slow, &
                                   weathering_po4, &
                                   occlusion_po4, &
                                   slow_exchange_po4, &
                                   fast_exchange_po4)

    USE mo_sb_constants
    USE mo_jsb_math_constants,  ONLY: eps4, eps8

    IMPLICIT NONE

    INTEGER,                                          INTENT(in)    :: nc, &                         !< dimensions
                                                                       nsoil_sb
    REAL(wp),                                         INTENT(in)    :: dtime                         !< timestep length
    CHARACTER(len=*),                                 INTENT(in)    :: sb_model_scheme               !< from sb config
    CHARACTER(len=*),                                 INTENT(in)    :: sb_adsorp_scheme              !< eca_full or eca_none
    LOGICAL ,                                         INTENT(in)    :: flag_sb_double_langmuir       !< T: double Langmuir; F: traditional Langmuir
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: qmax_po4_mineral_act, &       !< maximum sorption capacity of po4, [mol P/ m3]
                                                                       km_adsorpt_po4_act, &         !< half saturation of po4 sorption capacity, [mol P/ m3]
                                                                       qmax_fast_po4_act, &          !< maximum sorption capacity of fast_assoc_po4, [mol P/ m3]
                                                                       km_fast_po4_act, &            !< half saturation of po4 sorption in fast_assoc_po4, [mol P/ m3]
                                                                       qmax_slow_po4_act, &          !< maximum sorption capacity of slow_assoc_po4, [mol P/ m3]
                                                                       km_slow_po4_act, &            !< half saturation of po4 sorption in slow_assoc_po4, [mol P/ m3]
                                                                       k_adsorpt_po4_act, &          !< slow adsorption rate of po4, [1 / s]
                                                                       k_desorpt_po4_act, &          !< slow desorption rate of po4, [1 / s]
                                                                       p_deposition_sl, &            !< actual P deposition rate, [micro-mol P / m3 / s]
                                                                       biochem_mineralisation_po4, & !< biomineralisation rate [micro-mol P / m3 / s]
                                                                       transport_po4_solute, &       !< transport of soluble PO4 [micro-mol P / m3 / s]
                                                                       lateral_loss_po4_solute_sl, & !< lateral loss of soluble PO4 [micro-mol P / m3 / s]
                                                                       transport_po4_assoc_fast, &   !< transport of fast_assoc_po4 [micro-mol P / m3 / s]
                                                                       transport_po4_assoc_slow, &   !< transport of slow_assoc_po4 [micro-mol P / m3 / s]
                                                                       weathering_po4, &             !< phosphorus weathering rate [micro-mol P / m3 / s]
                                                                       mycorrhiza_uptake_po4_sl, &   !< mycorrhizal PO4 uptake [micro-mol P / m3 / s]
                                                                       plant_uptake_po4_sl, &        !< plant PO4 uptake [micro-mol P / m3 / s]
                                                                       microbial_uptake_po4_sl, &    !< microbial PO4 uptake [micro-mol P / m3 / s]
                                                                       gross_mineralisation_po4, &   !< gross mineralisation of PO4 [micro-mol P / m3 / s]
                                                                       net_mineralisation_po4        !< net mineralisation of PO4 [micro-mol P / m3 / s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: po4_solute, &                 !< PO4 in solution [mol / m3]
                                                                       po4_assoc_fast, &             !< fast minerally associated PO4 pool [mol/m3]
                                                                       po4_assoc_slow, &             !< slow minerally associated PO4 pool [mol/m3]
                                                                       occlusion_po4, &              !< phosphorus occlusion rate [micro-mol P / m3 / s]
                                                                       slow_exchange_po4             !< phosphorus desorption rate, assoc_fast --> assoc_slow [micro-mol P / m3 / s]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out)   :: fast_exchange_po4             !< phosphorus adsorption rate, solute --> assoc_fast [micro-mol P / m3 / s]

    REAL(wp), DIMENSION(nc, nsoil_sb)           :: hlp1, hlp2, hlp3, hlp4, &
                                                   partition_coef, & !< partition coefficient to soluble po4
                                                   partition_coef_1  !< partition coefficient to absorbed po4
    INTEGER                                     :: isoil, ichunk     !< loop over nsoil_sb and nc
    LOGICAL , DIMENSION(nc, nsoil_sb)           :: flag_arr          !< logical array helper
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_p_inorg_dynamics'

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 0.9 initilize local variables
    !>
    partition_coef(:,:)   = 0._wp
    partition_coef_1(:,:) = 0._wp
    flag_arr(:,:) = .TRUE.

    !>1.0 limit to possible assoc_slow
    !>
    hlp1(:,:) = k_desorpt_po4_act(:,:) * po4_assoc_slow(:,:) * dtime
    hlp2(:,:) = k_adsorpt_po4_act(:,:) * po4_assoc_fast(:,:) * dtime
    hlp3(:,:) = po4_assoc_slow(:,:) + hlp2(:,:) - hlp1(:,:) + (transport_po4_assoc_slow(:,:) - occlusion_po4(:,:)) * &
                dtime / 1000000._wp
    WHERE (hlp3(:,:) < eps8)
      hlp3(:,:) = po4_assoc_slow(:,:) + transport_po4_assoc_slow(:,:) * dtime / 1000000._wp + hlp2(:,:)
      hlp4(:,:) = hlp1(:,:) + occlusion_po4(:,:) * dtime / 1000000._wp
      WHERE (hlp4(:,:) < eps8)
        hlp1(:,:) = 0._wp
        occlusion_po4(:,:) = 0._wp
      ELSEWHERE
        hlp1(:,:) = hlp1(:,:) * hlp3(:,:) / hlp4(:,:)
        occlusion_po4(:,:) = occlusion_po4(:,:) * hlp3(:,:) / hlp4(:,:)
      ENDWHERE
      slow_exchange_po4(:,:) = (hlp2(:,:) - hlp1(:,:)) / dtime * 1000000._wp
    ENDWHERE

    !> 2.1 Langmuir equilibrium between soluble P and labile P, simple_1d and jsm model
    !> 2.1.1 traditional Langmuir isotherm
    !>
    IF (.NOT. flag_sb_double_langmuir) THEN ! if_iso2
      hlp1(:,:) = (weathering_po4(:,:) + biochem_mineralisation_po4(:,:) + net_mineralisation_po4(:,:) + &
                   transport_po4_solute(:,:) - lateral_loss_po4_solute_sl(:,:) + p_deposition_sl(:,:) - &
                   (plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:) + slow_exchange_po4(:,:)) + &
                   transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp
      hlp2(:,:) = (slow_exchange_po4(:,:) - transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp

      CALL calc_langmuir_kinetics(dtime, hlp1(:,:), &
                                  po4_solute(:,:),po4_assoc_fast(:,:), &
                                  km_adsorpt_po4_act(:,:), qmax_po4_mineral_act(:,:), &
                                  partition_coef(:,:))
    !> 2.1.2 double-surface Langmuir isotherm
    !>
    ELSE
      flag_arr(:,:) = ((qmax_slow_po4_act(:,:) - po4_assoc_slow(:,:))<eps4 .AND. po4_assoc_slow(:,:) > eps4)
      WHERE (flag_arr(:,:))
        slow_exchange_po4(:,:) = - k_desorpt_po4_act * po4_assoc_slow(:,:) * 1000000._wp
      ELSEWHERE
        slow_exchange_po4(:,:) = 0._wp
      ENDWHERE
      hlp1(:,:) = (weathering_po4(:,:) + biochem_mineralisation_po4(:,:) + net_mineralisation_po4(:,:) + &
                   transport_po4_solute(:,:) - lateral_loss_po4_solute_sl(:,:) + p_deposition_sl(:,:) - &
                   (plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:) + slow_exchange_po4(:,:)) + &
                   transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp
      hlp2(:,:) = ( - transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp

      DO ichunk = 1, nc
        DO isoil = 1, nsoil_sb
          IF (flag_arr(ichunk,isoil)) THEN
            CALL calc_langmuir_kinetics(dtime, hlp1(ichunk,isoil), &
                                        po4_solute(ichunk,isoil),po4_assoc_fast(ichunk,isoil), &
                                        km_fast_po4_act(ichunk,isoil), qmax_fast_po4_act(ichunk,isoil), &
                                        partition_coef(ichunk,isoil))
          ELSE
            CALL calc_langmuir_kinetics(dtime, hlp1(ichunk,isoil), &
                                        po4_solute(ichunk,isoil),po4_assoc_fast(ichunk,isoil), &
                                        km_adsorpt_po4_act(ichunk,isoil), qmax_po4_mineral_act(ichunk,isoil), &
                                        partition_coef(ichunk,isoil), &
                                        po4_assoc_slow(ichunk,isoil), &
                                        km_fast_po4_act(ichunk,isoil), qmax_fast_po4_act(ichunk,isoil), &
                                        km_slow_po4_act(ichunk,isoil), qmax_slow_po4_act(ichunk,isoil), &
                                        partition_coef_1(ichunk,isoil), slow_exchange_po4(ichunk,isoil))
          ENDIF
        ENDDO
      ENDDO
    ENDIF ! if_iso2

    !>3.0 update po4_solute, fast_assoc_po4, and slow_assoc_po4
    !>
    po4_solute(:,:)        = po4_solute(:,:) + partition_coef(:,:) * hlp1(:,:)
    po4_assoc_fast(:,:)    = po4_assoc_fast(:,:)  + (1._wp - partition_coef(:,:) - partition_coef_1(:,:)) * hlp1(:,:)
    fast_exchange_po4(:,:) = ((1._wp-partition_coef(:,:)-partition_coef_1(:,:)) * hlp1(:,:)+hlp2(:,:)) * 1.e6_wp / dtime
    po4_assoc_slow(:,:)    = po4_assoc_slow(:,:) + &
                            (slow_exchange_po4(:,:) - occlusion_po4(:,:) + transport_po4_assoc_slow(:,:)) * dtime / 1.e6_wp
  END SUBROUTINE calc_p_inorg_dynamics


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate Langmuir kinetics (modified from Yang et al. 2014)
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_langmuir_kinetics(dtime, solute_flux, solute, assoc_fast, &
                                              km, qmax, solute_partition_coef, &
                                              assoc_slow, km_fast, qmax_fast, km_slow, qmax_slow, &
                                              slow_partition_coef, slow_exchange)
    USE mo_sb_constants
    USE mo_jsb_math_constants,  ONLY: eps4, eps8

    IMPLICIT NONE
    ! 0.1 inout
    REAL(wp)          , INTENT(in)    :: dtime                  !< timestep length
    REAL(wp)          , INTENT(in)    :: solute_flux            !< total solute flux (mol/m3/dtime)
    REAL(wp)          , INTENT(in)    :: solute                 !< solute pool (mol/m3)
    REAL(wp)          , INTENT(in)    :: assoc_fast             !< adsorbed (fast) pool (mol/m3)
    REAL(wp)          , INTENT(in)    :: km                     !< half-saturation concentration of Langmuir isotherm (mol/m3)
    REAL(wp)          , INTENT(in)    :: qmax                   !< maximum sorption capacity of adsorbed pool (mol/m3)
    REAL(wp)          , INTENT(inout) :: solute_partition_coef  !< partition coefficient to soluble pool
    REAL(wp), OPTIONAL, INTENT(in)    :: assoc_slow             !< absorbed (slow) pool (mol/m3)
    REAL(wp), OPTIONAL, INTENT(in)    :: km_fast                !< half-saturation concentration of Langmuir isotherm (mol/m3)
    REAL(wp), OPTIONAL, INTENT(in)    :: qmax_fast              !< maximum sorption capacity of adsorbed pool (mol/m3)
    REAL(wp), OPTIONAL, INTENT(in)    :: km_slow                !< half-saturation concentration of Langmuir isotherm (mol/m3)
    REAL(wp), OPTIONAL, INTENT(in)    :: qmax_slow              !< maximum sorption capacity of adsorbed pool (mol/m3)
    REAL(wp), OPTIONAL, INTENT(inout) :: slow_partition_coef    !< partition coefficient to absorbed pool (po4 only)
    REAL(wp), OPTIONAL, INTENT(inout) :: slow_exchange          !< slow exchange rate between soluble and absorbed po4 pool [micro-mol P / m3 / s], double Langmuir only

    REAL(wp)                          :: hlp1, hlp2, hlp3, hlp4
    CHARACTER(len=*), PARAMETER       :: routine = TRIM(modname)//':calc_langmuir_kinetics'

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------
    !> calculate partition coefficient of soluble pool and absorbed pool (po4 double Langmuir only)
    !!


    !> 1.0 double Langmuir isotherm (po4 only)
    !!
    IF (PRESENT(assoc_slow)) THEN
      hlp1 = qmax_fast * km_fast / ((solute + km_fast)**2.0_wp )
      hlp2 = qmax_slow * km_slow / ((solute + km_slow)**2.0_wp )
      hlp3 = 1._wp
      IF ((qmax_fast * km_fast) < eps8 .AND. (qmax_slow * km_slow) < eps8) THEN
        hlp1 = 0._wp
        hlp2 = 0._wp
      ELSEIF ((qmax_fast * km_fast) < eps8) THEN
        hlp1 = 0._wp
      ELSEIF ((qmax_slow * km_slow) < eps8) THEN
        hlp2 = 0._wp
      ELSE
        IF ( assoc_fast < (-solute_flux) .OR. assoc_slow < (-solute_flux)) THEN
          hlp1 = assoc_fast / (assoc_fast + assoc_slow + solute)
          hlp2 = assoc_slow / (assoc_fast + assoc_slow + solute)
          hlp3 = solute     / (assoc_fast + assoc_slow + solute)
        ELSEIF ( (solute < po4_solute_min) .AND. solute_flux<0._wp ) THEN
          hlp3 = 0._wp
        ENDIF
      ENDIF
      solute_partition_coef = hlp3 / (hlp1 + hlp2 + hlp3)
      slow_partition_coef   = hlp2 / (hlp1 + hlp2 + hlp3)
      slow_exchange = slow_exchange + &
                      hlp2 / (hlp1 + hlp2 + hlp3) * solute_flux / dtime * 1000000._wp

    !> 1.1 check fast pool size
      hlp4 = assoc_fast + (1._wp - solute_partition_coef - slow_partition_coef) * solute_flux
      IF ((hlp4 > qmax) .AND. (qmax > eps8)) THEN
         IF (solute_flux < eps8) THEN
           solute_partition_coef = 0._wp
           slow_partition_coef   = 0._wp
         ELSEIF (solute_flux > eps8) THEN
           solute_partition_coef = solute_partition_coef / (solute_partition_coef + slow_partition_coef)
           slow_partition_coef   = slow_partition_coef   / (solute_partition_coef + slow_partition_coef)
         ENDIF
      ENDIF

    !> 2.0 traditional Langmuir isotherm
    ELSE
      IF((km > eps8) .AND. (qmax > eps8)) THEN
        !! calculate partition coefficient based on derivatives of solute
        hlp1 = qmax * km / ((solute + km) ** 2.0_wp )
        hlp3 = 1._wp
        solute_partition_coef = hlp3 / (hlp1 + hlp3)
        !! check for extreme low concentrations of soluble and adsorbed pools
        IF (assoc_fast < (-solute_flux)) THEN
          hlp1 = assoc_fast / (assoc_fast + solute)
          !! avoid assoc_fast to be lower than po4_solute_min
          IF ((assoc_fast + hlp1 * solute_flux) < po4_solute_min) THEN
            hlp1 = MAX(0._wp, (assoc_fast - po4_solute_min)) / (-solute_flux)
          ENDIF
          hlp3                  = 1._wp - hlp1
          solute_partition_coef = hlp3 / (hlp1 + hlp3)
        ELSEIF ((solute < po4_solute_min) .AND. solute_flux < 0._wp) THEN
          hlp3                  = 0._wp
          solute_partition_coef = hlp3 / (hlp1 + hlp3)
        ELSEIF (solute_flux < (-eps8)) THEN
          solute_partition_coef = MAX(hlp3 / (hlp1 + hlp3), solute / (solute + assoc_fast))
        ENDIF
      ELSE
        hlp1                  = 0._wp
        hlp3                  = 1._wp
        solute_partition_coef = hlp3 / (hlp1 + hlp3)
      ENDIF

      !> 2.1 check fast pool size against Qmax
      !!
      hlp4 = assoc_fast + (1._wp - solute_partition_coef) * solute_flux
      IF ((solute + assoc_fast) > eps8 ) THEN
         IF ((hlp4 > qmax) .AND. (qmax > eps8)) THEN
           IF (solute_flux < 0._wp) THEN
             solute_partition_coef = 1._wp - MIN(1._wp, MAX(0._wp, (qmax-po4_solute_min) / (-solute_flux)))
           ELSEIF (solute_flux > 0._wp) THEN
             solute_partition_coef = 1._wp
           ENDIF
         ENDIF
      ENDIF
    ENDIF  ! IF (PRESENT(assoc_slow))

  END SUBROUTINE calc_langmuir_kinetics

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to fasten the spinup by separately balancing slow soil pools
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_accelerated_slow_soil_pool_spinup(nsoil_sb, dtime, &
                                                  & soil_depth_sl, &
                                                  & spinup_accelerator_length, &
                                                  & formation_c_mavg_sl, &
                                                  & formation_c14_mavg_sl, &
                                                  & loss_c_mavg_sl, &
                                                  & pool_slow_sl)
    ! ----------------------------------------------------------------------------------------------------- !
    ! 0.1 inout
    INTEGER,  INTENT(in)    :: nsoil_sb                  !> soil dimension
    REAL(wp), INTENT(in)    :: dtime                     !> timestep length [s]
    REAL(wp), INTENT(in)    :: soil_depth_sl(:)          !> soil depth at each layer[m]
    INTEGER,  INTENT(in)    :: spinup_accelerator_length !> loop times of acceleration, read from namelist
    REAL(wp), INTENT(in)    :: formation_c_mavg_sl(:)    !> averaged slow soil C pool formation rate [mol m-3 timestep-1]
    REAL(wp), INTENT(in)    :: formation_c14_mavg_sl(:)  !> averaged slow soil C14 pool formation rate [mol m-3 timestep-1]
    REAL(wp), INTENT(in)    :: loss_c_mavg_sl(:)         !> averaged slow soil C pool loss rate [mol m-3timestep-1]
    REAL(wp), INTENT(inout) :: pool_slow_sl(:,:)         !> slow soil pool with all elements [mol m-3]
    ! ----------------------------------------------------------------------------------------------------- !
    ! 0.2 Local
    REAL, DIMENSION(nsoil_sb) :: c13c_slow_sl, nc_slow_sl, n15n_slow_sl,pn_slow_sl !> stiochiometry factors [mol / mol]
    REAL, DIMENSION(nsoil_sb) :: ratio_sl, ratio_c14_sl !> slow soil pool distributions across all soil layers[unitless]
    REAL                      :: tau_slow !> slow soil pool turnover rate [timestep]
    REAL                      :: pool_c_slow, pool_c14_slow !> the total slow soil pool [mol m-2]
    REAL                      :: formation_c_mavg, formation_c14_mavg, loss_c_mavg ! formation and loss rates of the total slow soil pool [mol m-2 timestep-1]
    INTEGER                   :: ilength, i_soil !counter for loops
    ! ------------------------------------------------------------------------------------------------------------
    !> accelerating the spin-up of slow pool
    !!

    !> 1.0 record current stoichiometry, apparent turnover rate of slow pool and slow pool distribution across all soil layers
    !!
    !> 1.1 record current stoichiometry
    DO i_soil = 1,nsoil_sb
      IF (pool_slow_sl(ixC,i_soil) > eps8) THEN
        c13c_slow_sl(i_soil) = pool_slow_sl(ixC13,i_soil) / pool_slow_sl(ixC,i_soil)
        nc_slow_sl(i_soil)   = pool_slow_sl(ixN,i_soil) / pool_slow_sl(ixC,i_soil)
        n15n_slow_sl(i_soil) = pool_slow_sl(ixN15,i_soil) / pool_slow_sl(ixN,i_soil)
        pn_slow_sl(i_soil)   = pool_slow_sl(ixP,i_soil) / pool_slow_sl(ixN,i_soil)
      ELSE
        c13c_slow_sl(i_soil) = 0._wp
        nc_slow_sl(i_soil)   = 0._wp
        n15n_slow_sl(i_soil) = 0._wp
        pn_slow_sl(i_soil)   = 0._wp
      END IF
    END DO

    !> 1.2 calculate the total turnover rate across all soil layers
    pool_c_slow = SUM(pool_slow_sl(ixC,:) * soil_depth_sl(:))
    pool_c14_slow = SUM(pool_slow_sl(ixC14,:) * soil_depth_sl(:))
    formation_c_mavg = SUM(formation_c_mavg_sl(:) * soil_depth_sl(:))
    formation_c14_mavg = SUM(formation_c14_mavg_sl(:) * soil_depth_sl(:))
    loss_c_mavg = SUM(loss_c_mavg_sl(:) * soil_depth_sl(:))

    !> 1.3 get the current slow soil pool profile
    IF ((pool_c_slow > eps8) .AND. (loss_c_mavg > eps8)) THEN
      tau_slow  = pool_c_slow / loss_c_mavg

      DO i_soil = 1,nsoil_sb
        ratio_sl(i_soil) = pool_slow_sl(ixC,i_soil) / pool_c_slow
        ratio_c14_sl(i_soil) = pool_slow_sl(ixC14,i_soil) / pool_c14_slow
      END DO

      !> 2.0 do the spin-up acceleration
      DO ilength = 1, spinup_accelerator_length
        pool_c_slow   = pool_c_slow + &
                            & ((-(pool_c_slow / tau_slow) + formation_c_mavg) / &
                            & dtime * one_day * one_year)
        !> C14 has the additional decay rate
        pool_c14_slow = pool_c14_slow + &
                            & ((-(pool_c14_slow / tau_slow) - (pool_c14_slow * lambda_C14 * dtime) + formation_c14_mavg) / &
                            & dtime * one_day * one_year)
      END DO

      !> 2.1 re-distribute the slow soil pool at each soil layer
      pool_slow_sl(ixC,:) = pool_c_slow * ratio_sl(:)
      pool_slow_sl(ixC14,:) = pool_c14_slow * ratio_c14_sl(:)

      !> 3.0 update other element pools
      pool_slow_sl(ixC13,:) = pool_slow_sl(ixC,:) * c13c_slow_sl(:)
      pool_slow_sl(ixN,:)   = pool_slow_sl(ixC,:) * nc_slow_sl(:)
      pool_slow_sl(ixN15,:) = pool_slow_sl(ixN,:) * n15n_slow_sl(:)
      pool_slow_sl(ixP,:)   = pool_slow_sl(ixN,:) * pn_slow_sl(:)
    END IF

  END SUBROUTINE calc_accelerated_slow_soil_pool_spinup

  ! ======================================================================================================= !
  !>
  !> Small helper routine to sum up total sb pool below ground soil (SOM and below ground litter) for a given element
  !>
  SUBROUTINE calculate_sb_pool_bg_soil_for_element(sb_pool_mt, nsoil_sb, ix_elem, sb_pool_total_bg_soil_elem)
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(IN)    :: sb_pool_mt(:,:,:,:)             !< sb pool for which the total BG soil is to be calculated
    INTEGER,  INTENT(IN)    :: nsoil_sb                        !< number of soil layers
    INTEGER,  INTENT(IN)    :: ix_elem                         !< index of the queried element in the sb pool
    REAL(wp), INTENT(INOUT) :: sb_pool_total_bg_soil_elem(:,:) !< total sb pool SOM + BG litter for given element
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER :: isoil
    ! ----------------------------------------------------------------------------------------------------- !
    sb_pool_total_bg_soil_elem(:,:) &
      & =  sb_pool_mt(ix_dom, ix_elem, :, :)       + sb_pool_mt(ix_dom_assoc, ix_elem, :, :)  &
      &  + sb_pool_mt(ix_fungi, ix_elem, :, :)     + sb_pool_mt(ix_mycorrhiza, ix_elem, :, :) &
      &  + sb_pool_mt(ix_microbial, ix_elem, :, :) + sb_pool_mt(ix_residue, ix_elem, :, :)    &
      &  + sb_pool_mt(ix_residue_assoc, ix_elem, :, :)

    DO isoil = 2,nsoil_sb
      sb_pool_total_bg_soil_elem(:,isoil) = sb_pool_total_bg_soil_elem(:,isoil) &
        & + sb_pool_mt(ix_soluable_litter, ix_elem, :, isoil)                   &
        & + sb_pool_mt(ix_polymeric_litter, ix_elem, :, isoil)                  &
        & + sb_pool_mt(ix_woody_litter, ix_elem, :, isoil)
    END DO

  END SUBROUTINE calculate_sb_pool_bg_soil_for_element

  ! ======================================================================================================= !
  !>
  !> Small helper routine to calculate total above ground litter in sb_pool for a given element
  !>
  SUBROUTINE calculate_sb_pool_ag_litter_for_element(sb_pool_mt, soil_depth_sl, ix_elem, sb_pool_total_ag_litter_elem)
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(IN)    :: sb_pool_mt(:,:,:,:)             !< sb pool for which the above ground litter is to be calculated
    REAL(wp), INTENT(IN)    :: soil_depth_sl(:,:)              !< soil depth at each layer[m]
    INTEGER,  INTENT(IN)    :: ix_elem                         !< index of the queried element in the sb pool
    REAL(wp), INTENT(INOUT) :: sb_pool_total_ag_litter_elem(:) !< total above ground litter for given element [mol m-2]
    ! ----------------------------------------------------------------------------------------------------- !
    sb_pool_total_ag_litter_elem(:) = (sb_pool_mt(ix_soluable_litter, ix_elem, :, 1)    &
      &                               + sb_pool_mt(ix_polymeric_litter, ix_elem, :, 1) &
      &                               + sb_pool_mt(ix_woody_litter, ix_elem, :, 1)) * soil_depth_sl(:,1)
  END SUBROUTINE calculate_sb_pool_ag_litter_for_element

#endif
END MODULE mo_q_sb_update_pools
