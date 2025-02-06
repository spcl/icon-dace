!> QUINCY 'simple soil model' calculate soil-biogeochemical pools
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
!>#### main routines for the SSM (simple soil model), calculating the soil biogeochemical pools
!>
MODULE mo_q_sb_ssm_main
#ifndef __NO_QUINCY__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message, message_text, finish

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID, SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, SB_BGCM_TRANSPORT_ID, &
    &                                   SB_BGCM_MYCO_EXPORT_ID

  USE mo_jsb_math_constants,      ONLY: one_day, one_year, eps8, eps4

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_simple_model
  PUBLIC :: calc_sb_inorganic_n15_fluxes  ! used by mo_q_sb_jsm_main

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sb_ssm_main'

CONTAINS

  ! ======================================================================================================= !
  !>update soil biogeochemical pools
  !>
  !>  routine to update soil litter and mineral pools (simple version for vegetation model testing)
  !>
  SUBROUTINE update_sb_simple_model(tile, options)
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: SB_, VEG_, SPQ_, A2L_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_veg_constants,         ONLY: ITREE
    USE mo_q_sb_jsm_transport,    ONLY: calc_bioturbation_rate
    USE mo_q_sb_litter_processes, ONLY: calc_litter_partitioning
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(SPQ_)
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
    INTEGER                               :: isoil                !< loop over nsoil_sb
    INTEGER                               :: ic, is               !< loop over dimensions
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: apparent_fast_som_cn      , &
                                             apparent_fast_som_np      , &
                                             fine_root_carbon_sl, &
                                             km1_uptake_nh4_act, &
                                             km2_uptake_nh4_act, &
                                             km1_uptake_no3_act, &
                                             km2_uptake_no3_act, &
                                             km1_uptake_po4_act, &
                                             km2_uptake_po4_act, &
                                             km_uptake_org_act, &
                                             km_mic_uptake_nh4_act, &
                                             km_mic_uptake_no3_act, &
                                             km_denitrification_c_act, &
                                             km_denitrification_no3_act, &
                                             p_deposition_act, &
                                             vmax_weath_mineral_act_sl, &
                                             k_adsorpt_po4_act, &
                                             k_desorpt_po4_act, &
                                             km_dip_biochem_po4_act, &
                                             qmax_po4_mineral_act
    REAL(wp)                              :: dtime
    INTEGER                               :: iblk, ics, ice, nc             !< dimensions
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_simple_model'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_formation_mt
    dsl4jsb_Def_mt2L3D :: sb_loss_mt
    dsl4jsb_Def_mt2L3D :: sb_transport_mt
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk :: unit_uptake_n_pot
    dsl4jsb_Real2D_onChunk :: unit_uptake_p_pot
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    ! A2L_ 2D
    dsl4jsb_Real2D_onChunk :: p_deposition
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl
    dsl4jsb_Real3D_onChunk :: percolation_sl
    dsl4jsb_Real3D_onChunk :: frac_w_lat_loss_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk :: leaching_nh4_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_n15_solute
    dsl4jsb_Real2D_onChunk :: emission_noy
    dsl4jsb_Real2D_onChunk :: emission_n2o
    dsl4jsb_Real2D_onChunk :: emission_n2
    dsl4jsb_Real2D_onChunk :: emission_noy_n15
    dsl4jsb_Real2D_onChunk :: emission_n2o_n15
    dsl4jsb_Real2D_onChunk :: emission_n2_n15
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: rtm_decomposition
    dsl4jsb_Real3D_onChunk :: rtm_plant_uptake
    dsl4jsb_Real3D_onChunk :: rtm_hsc
    dsl4jsb_Real3D_onChunk :: rtm_nitrification
    dsl4jsb_Real3D_onChunk :: rtm_denitrification
    dsl4jsb_Real3D_onChunk :: rtm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rtm_asymb_bnf
    dsl4jsb_Real3D_onChunk :: rmm_decomposition
    dsl4jsb_Real3D_onChunk :: rmm_plant_uptake
    dsl4jsb_Real3D_onChunk :: rmm_hsc
    dsl4jsb_Real3D_onChunk :: rmm_nitrification
    dsl4jsb_Real3D_onChunk :: rmm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_asymb_bnf
    dsl4jsb_Real3D_onChunk :: anaerobic_volume_fraction_sl
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: myc_respiration
    dsl4jsb_Real3D_onChunk :: myc_respiration_c13
    dsl4jsb_Real3D_onChunk :: myc_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: nitrification_no3
    dsl4jsb_Real3D_onChunk :: nitrification_noy
    dsl4jsb_Real3D_onChunk :: nitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_noy
    dsl4jsb_Real3D_onChunk :: denitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_n2
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: weathering_po4
    dsl4jsb_Real3D_onChunk :: occlusion_po4
    dsl4jsb_Real3D_onChunk :: slow_exchange_po4
    dsl4jsb_Real3D_onChunk :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: transport_nh4_assoc
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk :: transport_noy
    dsl4jsb_Real3D_onChunk :: transport_n2o
    dsl4jsb_Real3D_onChunk :: transport_n2
    dsl4jsb_Real3D_onChunk :: noy
    dsl4jsb_Real3D_onChunk :: n2o
    dsl4jsb_Real3D_onChunk :: n2
    dsl4jsb_Real3D_onChunk :: nh4_assoc
    dsl4jsb_Real3D_onChunk :: nh4_solute
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: po4_primary
    dsl4jsb_Real3D_onChunk :: k_bioturb
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation
    dsl4jsb_Real3D_onChunk :: qmax_po4
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: no3_n15_solute
    dsl4jsb_Real3D_onChunk :: noy_n15
    dsl4jsb_Real3D_onChunk :: n2o_n15
    dsl4jsb_Real3D_onChunk :: n2_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute
    dsl4jsb_Real3D_onChunk :: lateral_loss_nh4_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_no3_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_po4_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_nh4_n15_solute_sl
    dsl4jsb_Real3D_onChunk :: lateral_loss_no3_n15_solute_sl
    dsl4jsb_Real3D_onChunk :: transport_noy_n15
    dsl4jsb_Real3D_onChunk :: transport_n2o_n15
    dsl4jsb_Real3D_onChunk :: transport_n2_n15
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: nitrification_no3_n15
    dsl4jsb_Real3D_onChunk :: nitrification_noy_n15
    dsl4jsb_Real3D_onChunk :: nitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_noy_n15
    dsl4jsb_Real3D_onChunk :: denitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_n2_n15
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation_n15
    dsl4jsb_Real3D_onChunk :: vmax_weath_mineral_sl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    dtime     = options%dtime
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
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(apparent_fast_som_cn(nc, nsoil_sb))
    ALLOCATE(apparent_fast_som_np(nc, nsoil_sb))
    ALLOCATE(fine_root_carbon_sl(nc, nsoil_sb))
    ALLOCATE(km1_uptake_nh4_act(nc, nsoil_sb))
    ALLOCATE(km2_uptake_nh4_act(nc, nsoil_sb))
    ALLOCATE(km1_uptake_no3_act(nc, nsoil_sb))
    ALLOCATE(km2_uptake_no3_act(nc, nsoil_sb))
    ALLOCATE(km1_uptake_po4_act(nc, nsoil_sb))
    ALLOCATE(km2_uptake_po4_act(nc, nsoil_sb))
    ALLOCATE(km_uptake_org_act(nc, nsoil_sb))
    ALLOCATE(km_mic_uptake_nh4_act(nc, nsoil_sb))
    ALLOCATE(km_mic_uptake_no3_act(nc, nsoil_sb))
    ALLOCATE(km_denitrification_c_act(nc, nsoil_sb))
    ALLOCATE(km_denitrification_no3_act(nc, nsoil_sb))
    ALLOCATE(p_deposition_act(nc, nsoil_sb))
    ALLOCATE(vmax_weath_mineral_act_sl(nc, nsoil_sb))
    ALLOCATE(k_adsorpt_po4_act(nc, nsoil_sb))
    ALLOCATE(k_desorpt_po4_act(nc, nsoil_sb))
    ALLOCATE(km_dip_biochem_po4_act(nc, nsoil_sb))
    ALLOCATE(qmax_po4_mineral_act(nc, nsoil_sb))
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FORMATION_ID, sb_formation_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_LOSS_ID, sb_loss_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_TRANSPORT_ID, sb_transport_mt)
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_ 2D
    dsl4jsb_Get_var2D_onChunk(A2L_,      p_deposition)                ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_,      unit_uptake_n_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_,      unit_uptake_p_pot)
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_,      root_fraction_sl)
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_,      num_sl_above_bedrock)        ! in
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_,      soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SPQ_,      bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SPQ_,      percolation_sl)
    dsl4jsb_Get_var3D_onChunk(SPQ_,      frac_w_lat_loss_sl)
    dsl4jsb_Get_var3D_onChunk(SPQ_,      w_soil_sl)                   ! in
    ! SB_ 2D
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute)
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_n15_solute)
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute)
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute)
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_n15_solute)
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy)
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o)
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2)
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy_n15)
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o_n15)
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2_n15)
    ! SB_ 3D
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_decomposition)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_plant_uptake)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_denitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_asymb_bnf)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_decomposition)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_plant_uptake)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_asymb_bnf)
    dsl4jsb_Get_var3D_onChunk(SB_,       anaerobic_volume_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration)
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration_c13)
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration_c14)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2)
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       weathering_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       occlusion_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       slow_exchange_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       biochem_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       n2)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_fast)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_slow)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_occluded)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_nh4_solute_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_nh4_n15_solute_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_no3_solute_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_po4_solute_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_no3_n15_solute_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded)                ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_primary)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb)                   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       noy_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       n2_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation_n15)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       vmax_weath_mineral_sl)
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 provide soil profile values
    !>  calculate local variables
    !>
    CALL calc_local_soil_profile_values( &
      & nc, &                                   ! in
      & nsoil_sb, &
      & dtime, &
      & num_sl_above_bedrock(:), &
      & soil_depth_sl(:,:), &
      & bulk_dens_sl(:,:), &
      & w_soil_sl(:,:), &
      & nh4_solute(:,:), &
      & rtm_hsc(:,:), &
      & rmm_hsc(:,:), &
      & rtm_sorption(:,:), &
      & rtm_desorption(:,:), &
      & rmm_sorption(:,:), &
      & rmm_desorption(:,:), &
      & bulk_dens_corr_sl(:,:), &
      & p_deposition(:), &
      & vmax_weath_mineral_sl(:,:), &           ! in
      & km1_uptake_nh4_act(:,:), &              ! out
      & km2_uptake_nh4_act(:,:), &
      & km1_uptake_no3_act(:,:), &
      & km2_uptake_no3_act(:,:), &
      & km1_uptake_po4_act(:,:), &
      & km2_uptake_po4_act(:,:), &
      & km_uptake_org_act(:,:), &
      & km_mic_uptake_nh4_act, &
      & km_mic_uptake_no3_act, &
      & km_denitrification_c_act(:,:), &
      & km_denitrification_no3_act(:,:), &
      & apparent_fast_som_cn(:,:), &
      & apparent_fast_som_np(:,:), &
      & p_deposition_act(:,:), &
      & vmax_weath_mineral_act_sl(:,:), &
      & k_adsorpt_po4_act(:,:), &
      & k_desorpt_po4_act(:,:), &
      & km_dip_biochem_po4_act(:,:) )           ! out

    !>2.0 Calculate potential rates of uptake and losses
    !>

    !>  2.1 Potential plant nutrient uptake rates
    !>
    CALL calc_potential_plant_uptake_rate( &
      & nc, &                                     ! in
      & nsoil_sb, &
      & dtime, &
      & num_sl_above_bedrock(:), &
      & lctlib%vmax_uptake_n, &
      & lctlib%vmax_uptake_p, &
      & dsl4jsb_Config(SB_)%flag_mycorrhiza, &
      & soil_depth_sl(:,:), &
      & veg_pool_mt(ix_fine_root, ixC, :), &
      & sb_pool_mt(ix_mycorrhiza, ixC, :, :), &
      & root_fraction_sl(:,:), &
      & rtm_plant_uptake(:,:), &
      & rmm_plant_uptake(:,:), &
      & km1_uptake_nh4_act(:,:), &
      & km2_uptake_nh4_act(:,:), &
      & km1_uptake_no3_act(:,:), &
      & km2_uptake_no3_act(:,:), &
      & km1_uptake_po4_act(:,:), &
      & km2_uptake_po4_act(:,:), &
      & nh4_solute(:,:), &
      & no3_solute(:,:), &
      & po4_solute(:,:), &                        ! in
      & plant_uptake_nh4_sl(:,:), &               ! out
      & plant_uptake_no3_sl(:,:), &
      & plant_uptake_po4_sl(:,:), &
      & unit_uptake_n_pot(:), &
      & unit_uptake_p_pot(:) )                    ! out

    !>  2.2 Potential mycorrhiza nutrient uptake rates
    !>
    IF (dsl4jsb_Config(SB_)%flag_mycorrhiza) THEN
      CALL calc_potential_mycorrhiza_uptake_rate( &
        & nc, &                                     ! in
        & nsoil_sb, &
        & dtime, &
        & lctlib%vmax_uptake_n, &
        & lctlib%vmax_uptake_p, &
        & soil_depth_sl(:,:), &
        & veg_pool_mt(ix_fine_root, ixC, :), &
        & sb_pool_mt(ix_mycorrhiza, ixC, :, :), &
        & sb_pool_mt(ix_mycorrhiza, ixN, :, :), &
        & rtm_plant_uptake(:,:), &
        & rmm_plant_uptake(:,:), &
        & km1_uptake_nh4_act(:,:), &
        & km2_uptake_nh4_act(:,:), &
        & km1_uptake_no3_act(:,:), &
        & km2_uptake_no3_act(:,:), &
        & km1_uptake_po4_act(:,:), &
        & km2_uptake_po4_act(:,:), &
        & nh4_solute(:,:), &
        & no3_solute(:,:), &
        & po4_solute(:,:), &                        ! in
        & mycorrhiza_uptake_nh4_sl(:,:), &          ! out
        & mycorrhiza_uptake_no3_sl(:,:), &
        & mycorrhiza_uptake_po4_sl(:,:) )           ! out
    END IF

    !>  2.3 Calculate the potential decomposition rate of soil organic matter and litter pools
    !>
    CALL calc_potential_decomposition_rate( &
      & nc, &                                       ! in
      & nsoil_sb, &
      & dtime, &
      & soil_depth_sl(:,:), &
      & model%config%elements_index_map(:), &
      & model%config%is_element_used(:), &
      & dsl4jsb_Config(SB_)%flag_mycorrhiza_prim, &
      & dsl4jsb_Config(SB_)%flag_mycorrhiza_org, &
      & rtm_decomposition(:,:), &
      & rmm_decomposition(:,:), &
      & sb_pool_mt(:,:,:,:), &                      ! in
      & het_respiration(:,:), &                     ! inout
      & het_respiration_c13(:,:), &
      & het_respiration_c14(:,:), &
      & myc_respiration(:,:), &
      & myc_respiration_c13(:,:), &
      & myc_respiration_c14(:,:), &
      & sb_loss_mt(:,:,:,:) )                       ! inout

    !>  2.4 Loss pathways for the dynamic loss calculation scheme
    !>
    IF(TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme) == 'dynamic')THEN

      CALL calc_potential_nitrification_denitrification_rate( &
        & nc, &                                     ! in
        & nsoil_sb, &
        & dtime, &
        & num_sl_above_bedrock(:), &
        & rtm_nitrification(:,:), &
        & rmm_nitrification(:,:), &
        & rtm_denitrification(:,:), &
        & km_denitrification_c_act(:,:), &
        & km_denitrification_no3_act(:,:), &
        & anaerobic_volume_fraction_sl(:,:), &
        & nh4_solute(:,:), &
        & no3_solute(:,:), &
        & sb_pool_mt(ix_microbial, ixC, :, :), &    ! in
        & nitrification_no3(:,:), &                 ! out
        & denitrification_n2(:,:) )                 ! out

      CALL calc_inorganic_transport_rate( &
        & nc, &                                   ! in
        & nsoil_sb, &
        & dtime, &
        & num_sl_above_bedrock(:), &
        & soil_depth_sl(:,:), &
        & percolation_sl(:,:), &
        & frac_w_lat_loss_sl(:,:), &
        & rtm_gasdiffusion(:,:), &
        & rmm_gasdiffusion(:,:), &
        & nh4_solute(:,:), &
        & nh4_n15_solute(:,:), &
        & no3_solute(:,:), &
        & no3_n15_solute(:,:), &
        & po4_solute(:,:), &
        & noy(:,:), &
        & noy_n15(:,:), &
        & n2o(:,:), &
        & n2o_n15(:,:), &
        & n2(:,:), &
        & n2_n15(:,:), &                          ! in
        & transport_nh4_solute(:,:), &            ! inout
        & transport_nh4_n15_solute(:,:), &
        & transport_no3_solute(:,:), &
        & transport_no3_n15_solute(:,:), &
        & transport_po4_solute(:,:), &            ! inout
        & transport_noy(:,:), &                   ! out
        & transport_noy_n15(:,:), &
        & transport_n2o(:,:), &
        & transport_n2o_n15(:,:), &
        & transport_n2(:,:), &
        & transport_n2_n15(:,:), &
        & leaching_nh4_solute(:), &
        & leaching_nh4_n15_solute(:), &
        & leaching_no3_solute(:), &
        & leaching_no3_n15_solute(:), &
        & leaching_po4_solute(:), &
        & lateral_loss_nh4_solute_sl(:,:), &
        & lateral_loss_nh4_n15_solute_sl(:,:), &
        & lateral_loss_no3_solute_sl(:,:), &
        & lateral_loss_no3_n15_solute_sl(:,:), &
        & lateral_loss_po4_solute_sl(:,:), &
        & emission_noy(:), &
        & emission_noy_n15(:), &
        & emission_n2o(:), &
        & emission_n2o_n15(:), &
        & emission_n2(:), &
        & emission_n2_n15(:) )                    ! out

    ENDIF  ! IF(TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme)=='dynamic')

    !>  2.5 transport due to bioturbation
    !>
    ! calculate the bioturbation rate per soil layer (vertical transport rate [1/s] )
    !   from bulk density and root fraction and k_diff_org
    k_bioturb(:,:) = calc_bioturbation_rate(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                                     soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), root_fraction_sl(:,:))
    ! bioturbation transport for the 'simple_1d' soil model
    CALL calc_bioturbation_transport_wrapper_s1d( &
      & nc, &                                 ! in
      & nsoil_sb, &
      & dtime, &
      & num_sl_above_bedrock(:), &
      & soil_depth_sl(:,:), &
      & model%config%elements_index_map(:), &
      & model%config%is_element_used(:), &
      & k_bioturb(:,:), &
      & nh4_assoc, &
      & nh4_n15_assoc, &
      & po4_assoc_fast, &
      & po4_assoc_slow, &
      & po4_occluded, &
      & sb_pool_mt(:,:,:,:), &                ! in
      & transport_nh4_assoc(:,:), &           ! inout
      & transport_nh4_n15_assoc(:,:), &
      & transport_po4_assoc_fast(:,:), &
      & transport_po4_assoc_slow(:,:), &
      & transport_po4_occluded(:,:), &
      & sb_transport_mt(:,:,:,:) )            ! inout

    !>  2.6 Calculate the potential PO4 flux rate: biomineralisation, weathering, slow sorption, occlusion
    !>
    ! calculate the fine root C of each soil layer, fine_root_carbon_sl [mol C m-3]
    DO isoil = 1,nsoil_sb
      WHERE(soil_depth_sl(:,isoil) > eps8)
        fine_root_carbon_sl(:,isoil) = root_fraction_sl(:,isoil) / soil_depth_sl(:,isoil) * veg_pool_mt(ix_fine_root, ixC, :)
      ELSEWHERE
        fine_root_carbon_sl(:,isoil) = 0.0_wp
      ENDWHERE
    ENDDO
    ! calculate the potential PO4 flux rate
    CALL calc_potential_p_flux_rate( &
      & nc, &                                         ! in
      & nsoil_sb, &
      & dtime, &
      & num_sl_above_bedrock(:), &
      & dsl4jsb_Config(SB_)%flag_sb_prescribe_po4, &
      & rtm_decomposition(:,:), &
      & rmm_decomposition(:,:), &
      & km_dip_biochem_po4_act(:,:), &
      & bulk_dens_corr_sl(:,:), &
      & fine_root_carbon_sl(:,:), &
      & po4_solute(:,:), &
      & po4_assoc_fast(:,:), &
      & po4_assoc_slow(:,:), &
      & po4_occluded(:,:), &
      & po4_primary(:,:), &
      & vmax_weath_mineral_act_sl(:,:), &
      & k_adsorpt_po4_act(:,:), &
      & k_desorpt_po4_act(:,:), &
      & sb_pool_mt(ix_residue, ixP, :, :), &          ! in
      & weathering_po4(:,:), &                        ! out
      & occlusion_po4(:,:), &
      & slow_exchange_po4(:,:), &
      & biochem_mineralisation_po4(:,:) )             ! out

    !>3.0 Actual fluxes of in and out of SOM
    !>

    !>  3.1 Calculation of asymbiotic N fixation rate
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        asymb_n_fixation(ic, is) = calc_asymb_n_fixation( &
          &                         TRIM(dsl4jsb_Config(SB_)%bnf_scheme), &
          &                         sb_loss_mt(ix_soluable_litter, ixC, ic, is) + sb_loss_mt(ix_polymeric_litter, ixC, ic, is), &
          &                         sb_loss_mt(ix_soluable_litter, ixN, ic, is) + sb_loss_mt(ix_polymeric_litter, ixN, ic, is), &
          &                         apparent_fast_som_cn(ic, is), &
          &                         rtm_asymb_bnf(ic, is), &
          &                         rmm_asymb_bnf(ic, is), &
          &                         nh4_solute(ic, is), &
          &                         no3_solute(ic, is) )
      ENDDO
    ENDDO

    !>  3.2 Calculation of actual plant, mycorrhizal, microbial uptake rate as well as
    !>      nitrification, denitrification and biomineralisation rates
    !>
    CALL calc_actual_flux_rates( &
      & nc, &                                     ! in
      & nsoil_sb, &
      & num_sl_above_bedrock(:), &
      & TRIM(dsl4jsb_Config(SB_)%bnf_scheme), &
      & nh4_solute(:,:), &
      & no3_solute(:,:), &
      & po4_solute(:,:), &
      & po4_assoc_fast(:,:), &
      & transport_nh4_solute(:,:), &
      & transport_no3_solute(:,:), &
      & transport_po4_solute(:,:), &
      & transport_po4_assoc_fast(:,:), &
      & weathering_po4(:,:), &
      & p_deposition_act(:,:), &
      & slow_exchange_po4(:,:), &
      & apparent_fast_som_cn(:,:), &
      & km_mic_uptake_nh4_act(:,:), &
      & km_mic_uptake_no3_act(:,:), &             ! in
      & plant_uptake_nh4_sl(:,:), &               ! inout
      & plant_uptake_no3_sl(:,:), &
      & plant_uptake_po4_sl(:,:), &
      & mycorrhiza_uptake_nh4_sl(:,:), &
      & mycorrhiza_uptake_no3_sl(:,:), &
      & mycorrhiza_uptake_po4_sl(:,:), &
      & nitrification_no3(:,:), &
      & denitrification_n2(:,:), &
      & asymb_n_fixation(:,:), &
      & biochem_mineralisation_po4(:,:), &
      & sb_loss_mt(:,:,:,:), &                    ! inout
      & microbial_uptake_nh4_sl(:,:), &           ! out
      & microbial_uptake_no3_sl(:,:), &           ! out
      & microbial_uptake_po4_sl(:,:) )            ! out

    !>  3.3 ...
    !>
    IF (TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme) == 'dynamic') THEN
      CALL calc_partitioning_nitrification_denitrification_rate( &
        & nc, &
        & nsoil_sb, &
        & num_sl_above_bedrock(:), &
        & rtm_denitrification(:,:), &
        & nitrification_no3(:,:), &
        & denitrification_n2(:,:), &
        & nitrification_noy(:,:), &
        & nitrification_n2o(:,:), &
        & denitrification_noy(:,:), &
        & denitrification_n2o(:,:))
    ENDIF

    !>  3.4 update n15 for inorganic fluxes
    !>      here the function is called for the simple_1d soil-model: flux rates in [mol/m3/timestep]
    !>
    !>    NOTE: JSM is using the same function but uses flux rates in [micro-mol/m3/s]; unit conversion is applied within the routine
    !>
    CALL calc_sb_inorganic_n15_fluxes( &
      & nc, &                                         ! in
      & nsoil_sb, &
      & dtime, &
      & num_sl_above_bedrock(:), &
      & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
      & nh4_solute(:,:), &
      & nh4_n15_solute(:,:), &
      & no3_solute(:,:), &
      & no3_n15_solute(:,:), &
      & plant_uptake_nh4_sl(:,:), &
      & plant_uptake_no3_sl(:,:), &
      & mycorrhiza_uptake_nh4_sl(:,:), &
      & mycorrhiza_uptake_no3_sl(:,:), &
      & microbial_uptake_nh4_sl(:,:), &
      & microbial_uptake_no3_sl(:,:), &
      & asymb_n_fixation(:,:), &
      & nitrification_no3(:,:), &
      & nitrification_noy(:,:), &
      & nitrification_n2o(:,:), &
      & denitrification_noy(:,:), &
      & denitrification_n2o(:,:), &
      & denitrification_n2(:,:), &
      & transport_nh4_n15_solute(:,:), &
      & transport_no3_n15_solute(:,:), &              ! in
      & plant_uptake_nh4_n15_sl(:,:), &               ! out
      & plant_uptake_no3_n15_sl(:,:), &
      & mycorrhiza_uptake_nh4_n15_sl(:,:), &
      & mycorrhiza_uptake_no3_n15_sl(:,:), &
      & microbial_uptake_nh4_n15_sl(:,:), &
      & microbial_uptake_no3_n15_sl(:,:), &
      & asymb_n_fixation_n15(:,:), &
      & nitrification_no3_n15(:,:), &
      & nitrification_noy_n15(:,:), &
      & nitrification_n2o_n15(:,:), &
      & denitrification_noy_n15(:,:), &
      & denitrification_n2o_n15(:,:), &
      & denitrification_n2_n15(:,:) )                 ! out

    !>  3.5 formation of SOM and associated respiration, net mineralisation fluxes
    !>
    CALL calc_SOM_formation( &
      & nc, &                                 ! in
      & nsoil_sb, &
      & num_sl_above_bedrock(:), &
      & apparent_fast_som_cn(:,:), &
      & apparent_fast_som_np(:,:), &
      & microbial_uptake_nh4_sl(:,:), &
      & microbial_uptake_nh4_n15_sl(:,:), &
      & microbial_uptake_no3_sl(:,:), &
      & microbial_uptake_no3_n15_sl(:,:), &
      & asymb_n_fixation(:,:), &
      & asymb_n_fixation_n15(:,:), &
      & microbial_uptake_po4_sl(:,:), &       ! in
      & het_respiration(:,:), &               ! inout
      & het_respiration_c13(:,:), &
      & het_respiration_c14(:,:), &
      & net_mineralisation_nh4(:,:), &
      & net_mineralisation_no3(:,:), &
      & net_mineralisation_po4(:,:), &
      & net_mineralisation_nh4_n15(:,:), &
      & net_mineralisation_no3_n15(:,:), &
      & biochem_mineralisation_po4(:,:), &
      & sb_loss_mt(:,:,:,:), &
      & sb_formation_mt(:,:,:,:) )            ! inout

    !>  3.6 formation of metabolic and strucutral litter by woody liter decay and
    !>     associated respiration, net mineralisation fluxes
    !>
    IF (lctlib%growthform == ITREE) THEN
      CALL calc_litter_formation_from_woody_decay( &
        & nc, &                               ! in
        & nsoil_sb, &
        & sb_loss_mt(:,:,:,:), &              ! in
        & het_respiration(:,:), &             ! inout
        & het_respiration_c13(:,:), &
        & het_respiration_c14(:,:), &
        & net_mineralisation_nh4(:,:), &
        & net_mineralisation_po4(:,:), &
        & net_mineralisation_nh4_n15(:,:), &
        & sb_formation_mt(:,:,:,:) )          ! inout
    ENDIF

    !>4.0 Mycorrhiza dynamics
    !>
    CALL calc_mycorrhiza_dynamics( &
      & nc, &                                       ! in
      & nsoil_sb, &
      & dtime, &
      & soil_depth_sl(:,:), &
      & model%config%elements_index_map(:), &
      & model%config%is_element_used(:), &
      & lctlib%tau_mycorrhiza, &
      & dsl4jsb_Config(SB_)%flag_mycorrhiza_org, &
      & nh4_solute(:,:), &
      & no3_solute(:,:), &
      & km2_uptake_nh4_act(:,:), &
      & km2_uptake_no3_act(:,:), &
      & km_uptake_org_act(:,:), &
      & veg_pool_mt(ix_fine_root, ixC, :), &
      & sb_pool_mt(:,:,:,:), &                      ! in
      & mycorrhiza_uptake_nh4_sl(:,:), &            ! inout
      & mycorrhiza_uptake_nh4_n15_sl(:,:), &
      & mycorrhiza_uptake_no3_sl(:,:), &
      & mycorrhiza_uptake_no3_n15_sl(:,:), &
      & mycorrhiza_uptake_norg_sl(:,:), &
      & mycorrhiza_uptake_norg_n15_sl(:,:), &
      & mycorrhiza_uptake_po4_sl(:,:), &
      & het_respiration(:,:), &
      & het_respiration_c13(:,:), &
      & het_respiration_c14(:,:), &
      & myc_respiration(:,:), &
      & myc_respiration_c13(:,:), &
      & myc_respiration_c14(:,:), &
      & sb_formation_mt(:,:,:,:), &
      & sb_loss_mt(:,:,:,:), &
      & sb_mycorrhiza_export_mt(:,:,:) )            ! inout

    !>5.0 Update soil litter pools with litter fall
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

    !>6.0 for 'sb_model = fixed': estimate N losses
    !>
    IF (TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme) == 'fixed') THEN
      CALL calc_fixed_nloss_rate( &
        & nc, &
        & nsoil_sb, &
        & soil_depth_sl(:,:), &
        & net_mineralisation_nh4(:,:), &
        & net_mineralisation_nh4_n15(:,:), &
        & nitrification_no3(:,:), &
        & nitrification_no3_n15(:,:), &
        & volatilisation_nh4(:,:), &
        & volatilisation_nh4_n15(:,:), &
        & emission_n2(:), &
        & emission_n2_n15(:))
    ENDIF

    !>7.0 unit conversion of rates and fluxes
    !>
    !>    convert rates to [micro-mol/ m3 / s]
    !>    covert surface fluxes to [micro-mol / m2 / s]
    !>
    het_respiration(:,:)     = het_respiration(:,:) / dtime * 1000000._wp
    het_respiration_c13(:,:) = het_respiration_c13(:,:) / dtime * 1000000._wp
    het_respiration_c14(:,:) = het_respiration_c14(:,:) / dtime * 1000000._wp
    myc_respiration(:,:)     = myc_respiration(:,:) / dtime * 1000000._wp
    myc_respiration_c13(:,:) = myc_respiration_c13(:,:) / dtime * 1000000._wp
    myc_respiration_c14(:,:) = myc_respiration_c14(:,:) / dtime * 1000000._wp

    net_mineralisation_nh4(:,:)    = net_mineralisation_nh4(:,:) / dtime * 1000000._wp
    net_mineralisation_no3(:,:)    = net_mineralisation_no3(:,:) / dtime * 1000000._wp
    net_mineralisation_po4(:,:)    = net_mineralisation_po4(:,:) / dtime * 1000000._wp
    plant_uptake_nh4_sl(:,:)       = plant_uptake_nh4_sl(:,:) / dtime * 1000000._wp
    plant_uptake_no3_sl(:,:)       = plant_uptake_no3_sl(:,:) / dtime * 1000000._wp
    plant_uptake_po4_sl(:,:)       = plant_uptake_po4_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_nh4_sl(:,:)  = mycorrhiza_uptake_nh4_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_no3_sl(:,:)  = mycorrhiza_uptake_no3_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_norg_sl(:,:) = mycorrhiza_uptake_norg_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_po4_sl(:,:)  = mycorrhiza_uptake_po4_sl(:,:) / dtime * 1000000._wp
    microbial_uptake_nh4_sl(:,:)   = microbial_uptake_nh4_sl(:,:) / dtime * 1000000._wp
    microbial_uptake_no3_sl(:,:)   = microbial_uptake_no3_sl(:,:) / dtime * 1000000._wp
    microbial_uptake_po4_sl(:,:)   = microbial_uptake_po4_sl(:,:) / dtime * 1000000._wp
    asymb_n_fixation(:,:)          = asymb_n_fixation(:,:) / dtime * 1000000._wp

    volatilisation_nh4(:,:)  = volatilisation_nh4(:,:) / dtime * 1000000._wp
    nitrification_no3(:,:)   = nitrification_no3(:,:) / dtime * 1000000._wp
    nitrification_noy(:,:)   = nitrification_noy(:,:) / dtime * 1000000._wp
    nitrification_n2o(:,:)   = nitrification_n2o(:,:) / dtime * 1000000._wp
    denitrification_noy(:,:) = denitrification_noy(:,:) / dtime * 1000000._wp
    denitrification_n2o(:,:) = denitrification_n2o(:,:) / dtime * 1000000._wp
    denitrification_n2(:,:)  = denitrification_n2(:,:) / dtime * 1000000._wp

    transport_nh4_solute(:,:)     = transport_nh4_solute(:,:) / dtime * 1000000._wp
    transport_nh4_assoc(:,:)      = transport_nh4_assoc(:,:) / dtime * 1000000._wp
    transport_no3_solute(:,:)     = transport_no3_solute(:,:) / dtime * 1000000._wp
    transport_po4_solute(:,:)     = transport_po4_solute(:,:) / dtime * 1000000._wp
    transport_po4_assoc_fast(:,:) = transport_po4_assoc_fast(:,:) / dtime * 1000000._wp
    transport_po4_assoc_slow(:,:) = transport_po4_assoc_slow(:,:) / dtime * 1000000._wp
    transport_po4_occluded(:,:)   = transport_po4_occluded(:,:) / dtime * 1000000._wp
    transport_noy(:,:)            = transport_noy(:,:) / dtime * 1000000._wp
    transport_n2o(:,:)            = transport_n2o(:,:) / dtime * 1000000._wp
    transport_n2(:,:)             = transport_n2(:,:) / dtime * 1000000._wp

    lateral_loss_nh4_solute_sl(:,:) = lateral_loss_nh4_solute_sl(:,:) / dtime * 1000000._wp
    lateral_loss_no3_solute_sl(:,:) = lateral_loss_no3_solute_sl(:,:) / dtime * 1000000._wp
    lateral_loss_po4_solute_sl(:,:) = lateral_loss_po4_solute_sl(:,:) / dtime * 1000000._wp

    leaching_nh4_solute(:) = leaching_nh4_solute(:) / dtime * 1000000._wp
    leaching_no3_solute(:) = leaching_no3_solute(:) / dtime * 1000000._wp
    leaching_po4_solute(:) = leaching_po4_solute(:) / dtime * 1000000._wp
    emission_noy(:)        = emission_noy(:) / dtime * 1000000._wp
    emission_n2o(:)        = emission_n2o(:) / dtime * 1000000._wp
    emission_n2(:)         = emission_n2(:) / dtime * 1000000._wp

    net_mineralisation_nh4_n15(:,:)   = net_mineralisation_nh4_n15(:,:) / dtime * 1000000._wp
    net_mineralisation_no3_n15(:,:)   = net_mineralisation_no3_n15(:,:) / dtime * 1000000._wp
    plant_uptake_nh4_n15_sl(:,:)      = plant_uptake_nh4_n15_sl(:,:) / dtime * 1000000._wp
    plant_uptake_no3_n15_sl(:,:)      = plant_uptake_no3_n15_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_nh4_n15_sl(:,:) = mycorrhiza_uptake_nh4_n15_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_no3_n15_sl(:,:) = mycorrhiza_uptake_no3_n15_sl(:,:) / dtime * 1000000._wp
    microbial_uptake_nh4_n15_sl(:,:)  = microbial_uptake_nh4_n15_sl(:,:) / dtime * 1000000._wp
    microbial_uptake_no3_n15_sl(:,:)  = microbial_uptake_no3_n15_sl(:,:) / dtime * 1000000._wp
    asymb_n_fixation_n15(:,:)         = asymb_n_fixation_n15(:,:) / dtime * 1000000._wp

    volatilisation_nh4_n15(:,:)  = volatilisation_nh4_n15(:,:) / dtime * 1000000._wp
    nitrification_no3_n15(:,:)   = nitrification_no3_n15(:,:) / dtime * 1000000._wp
    nitrification_noy_n15(:,:)   = nitrification_noy_n15(:,:) / dtime * 1000000._wp
    nitrification_n2o_n15(:,:)   = nitrification_n2o_n15(:,:) / dtime * 1000000._wp
    denitrification_noy_n15(:,:) = denitrification_noy_n15(:,:) / dtime * 1000000._wp
    denitrification_n2o_n15(:,:) = denitrification_n2o_n15(:,:) / dtime * 1000000._wp
    denitrification_n2_n15(:,:)  = denitrification_n2_n15(:,:) / dtime * 1000000._wp

    transport_nh4_n15_solute(:,:) = transport_nh4_n15_solute(:,:) / dtime * 1000000._wp
    transport_nh4_n15_assoc(:,:)  = transport_nh4_n15_assoc(:,:) / dtime * 1000000._wp
    transport_no3_n15_solute(:,:) = transport_no3_n15_solute(:,:) / dtime * 1000000._wp
    transport_noy_n15(:,:)        = transport_noy_n15(:,:) / dtime * 1000000._wp
    transport_n2o_n15(:,:)        = transport_n2o_n15(:,:) / dtime * 1000000._wp
    transport_n2_n15(:,:)         = transport_n2_n15(:,:) / dtime * 1000000._wp

    leaching_nh4_n15_solute(:)          = leaching_nh4_n15_solute(:) / dtime * 1000000._wp
    leaching_no3_n15_solute(:)          = leaching_no3_n15_solute(:) / dtime * 1000000._wp
    lateral_loss_nh4_n15_solute_sl(:,:) = lateral_loss_nh4_n15_solute_sl(:,:) / dtime * 1000000._wp
    lateral_loss_no3_n15_solute_sl(:,:) = lateral_loss_no3_n15_solute_sl(:,:) / dtime * 1000000._wp

    emission_noy_n15(:)        = emission_noy_n15(:) / dtime * 1000000._wp
    emission_n2o_n15(:)        = emission_n2o_n15(:) / dtime * 1000000._wp
    emission_n2_n15(:)         = emission_n2_n15(:) / dtime * 1000000._wp

    weathering_po4(:,:)             = weathering_po4(:,:) / dtime * 1000000._wp
    occlusion_po4(:,:)              = occlusion_po4(:,:) / dtime * 1000000._wp
    slow_exchange_po4(:,:)          = slow_exchange_po4(:,:) / dtime * 1000000._wp
    biochem_mineralisation_po4(:,:) = biochem_mineralisation_po4(:,:) / dtime * 1000000._wp

    !------------------------------------------------------------------------------------------------------ !
    ! de-allocate local allocatable variables
    DEALLOCATE(apparent_fast_som_cn)
    DEALLOCATE(apparent_fast_som_np)
    DEALLOCATE(fine_root_carbon_sl)
    DEALLOCATE(km1_uptake_nh4_act)
    DEALLOCATE(km2_uptake_nh4_act)
    DEALLOCATE(km1_uptake_no3_act)
    DEALLOCATE(km2_uptake_no3_act)
    DEALLOCATE(km1_uptake_po4_act)
    DEALLOCATE(km2_uptake_po4_act)
    DEALLOCATE(km_uptake_org_act)
    DEALLOCATE(km_mic_uptake_nh4_act)
    DEALLOCATE(km_mic_uptake_no3_act)
    DEALLOCATE(km_denitrification_c_act)
    DEALLOCATE(km_denitrification_no3_act)
    DEALLOCATE(p_deposition_act)
    DEALLOCATE(vmax_weath_mineral_act_sl)
    DEALLOCATE(k_adsorpt_po4_act)
    DEALLOCATE(k_desorpt_po4_act)
    DEALLOCATE(km_dip_biochem_po4_act)
    DEALLOCATE(qmax_po4_mineral_act)

  END SUBROUTINE update_sb_simple_model

  ! ======================================================================================================= !
  !>provide soil profile values, i.e., local variables
  !>
  SUBROUTINE calc_local_soil_profile_values( &
    & nc, &
    & nsoil_sb, &
    & dtime, &
    & num_sl_above_bedrock, &
    & soil_depth_sl, &
    & bulk_dens_sl, &
    & w_soil_sl, &
    & nh4_solute, &
    & rtm_hsc, &
    & rmm_hsc, &
    & rtm_sorption, &
    & rtm_desorption, &
    & rmm_sorption, &
    & rmm_desorption, &
    & bulk_dens_corr_sl, &
    & p_deposition, &
    & vmax_weath_mineral_sl, &
    & km1_uptake_nh4_act, &
    & km2_uptake_nh4_act, &
    & km1_uptake_no3_act, &
    & km2_uptake_no3_act, &
    & km1_uptake_po4_act, &
    & km2_uptake_po4_act, &
    & km_uptake_org_act, &
    & km_mic_uptake_nh4_act, &
    & km_mic_uptake_no3_act, &
    & km_denitrification_c_act, &
    & km_denitrification_no3_act, &
    & apparent_fast_som_cn, &
    & apparent_fast_som_np, &
    & p_deposition_act, &
    & vmax_weath_mineral_act_sl, &
    & k_adsorpt_po4_act, &
    & k_desorpt_po4_act, &
    & km_dip_biochem_po4_act )

    USE mo_sb_constants,            ONLY: km_denitrification_c, km_denitrification_no3, km1_up_nh4, km1_up_no3, &
                                          km1_up_po4, km2_up_nh4, km2_up_no3, km2_up_po4, km_up_org, &
                                          k_fast_som_cn_min, k_fast_som_cn_max, km_mic_uptake_nh4, km_mic_uptake_no3, &
                                          f_fast_som_cn, k_fast_som_np, kref_fast_som_np, k_desorpt_po4, &
                                          k_adsorpt_po4, km_dip_biochem_po4
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                                     INTENT(in)  :: nc                           !< dimensions
    INTEGER,                                     INTENT(in)  :: nsoil_sb                     !< number of soil layers
    REAL(wp),                                    INTENT(in)  :: dtime                        !< timestep length
    REAL(wp), DIMENSION(nc),                     INTENT(in)  :: num_sl_above_bedrock         !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc),                     INTENT(in)  :: p_deposition                 !< p deposition, [micro-mol m-2 s-1]
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(in)  :: &
                                                                soil_depth_sl, &             !< soil depth of a layer [m]
                                                                bulk_dens_sl, &              !< bulk density of a layer [kg m-3]
                                                                w_soil_sl, &                 !< water content of a layer [m]
                                                                nh4_solute, &                !< NH4 content of a layer [mol m-3]
                                                                rtm_hsc, &                   !< temperature response scalar for half-saturation constants
                                                                rmm_hsc, &                   !< moisture response scalar for half-saturation constants
                                                                rtm_sorption, &              !< temperature response scalar for po4 sorption
                                                                rtm_desorption, &            !< temperature response scalar for po4 desorption
                                                                rmm_sorption, &              !< moisture response scalar for po4 sorption
                                                                rmm_desorption, &            !< moisture response scalar for po4 desorption
                                                                bulk_dens_corr_sl, &         !< bulk density of the soil layer, [kg/m3]
                                                                vmax_weath_mineral_sl        !< weathering coefficient of mineral soil [mol P m-3 s-1]
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(out) :: &
                                                                km1_uptake_nh4_act, &        !< actual km1 for NH4 uptake [m3 / mol]
                                                                km2_uptake_nh4_act, &        !< actual km2 for NH4 uptake [mol / m3]
                                                                km1_uptake_no3_act, &        !< actual km1 for NO3 uptake [m3 / mol]
                                                                km2_uptake_no3_act, &        !< actual km2 for NO3 uptake [mol / m3]
                                                                km1_uptake_po4_act, &        !< actual km1 for PO4 uptake [m3 / mol]
                                                                km2_uptake_po4_act, &        !< actual km2 for PO4 uptake [mol / m3]
                                                                km_uptake_org_act, &         !< actual km for organic uptake [mol / m3]
                                                                km_mic_uptake_nh4_act, &     !< actual km for microbial NH4 uptake [mol / m3]
                                                                km_mic_uptake_no3_act, &     !< actual km for microbial NO3 uptake [m3 / mol]
                                                                km_denitrification_c_act, &  !< actual km on C for denitrification [mol / m3]
                                                                km_denitrification_no3_act, & !< actual km on NO3 for denitrification [mol / m3]
                                                                apparent_fast_som_cn, &           !< apparent microbial molar carbon-nitrogen ratio
                                                                apparent_fast_som_np, &           !< apparent microbial molar carbon-nitrogen ratio
                                                                p_deposition_act, &          !< actual p deposition, [mol m-3 timestep-1]
                                                                vmax_weath_mineral_act_sl, & !< actual weathering rate for primary P [mol P kg-1 s-1]
                                                                k_adsorpt_po4_act, &         !< actual adsorption rate for fast_assoc_po4 [s-1]
                                                                k_desorpt_po4_act, &         !< actual desorption rate for slow_assoc_po4 [s-1]
                                                                km_dip_biochem_po4_act    !< actual km for soluble P control on biochem_mineralisation_po4
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                         :: ic, is
    CHARACTER(len=*), PARAMETER     :: routine = TRIM(modname)//':calc_local_soil_profile_values'
    ! ----------------------------------------------------------------------------------------------------- !

    !> init out var
    !>
    km1_uptake_nh4_act(:,:)         = 0.0_wp
    km2_uptake_nh4_act(:,:)         = 0.0_wp
    km1_uptake_no3_act(:,:)         = 0.0_wp
    km2_uptake_no3_act(:,:)         = 0.0_wp
    km1_uptake_po4_act(:,:)         = 0.0_wp
    km2_uptake_po4_act(:,:)         = 0.0_wp
    km_uptake_org_act(:,:)          = 0.0_wp
    km_mic_uptake_nh4_act(:,:)      = 0.0_wp
    km_mic_uptake_no3_act(:,:)      = 0.0_wp
    km_denitrification_c_act(:,:)   = 0.0_wp
    km_denitrification_no3_act(:,:) = 0.0_wp
    apparent_fast_som_cn(:,:)       = 0.0_wp
    apparent_fast_som_np(:,:)       = 0.0_wp
    p_deposition_act(:,:)           = 0.0_wp
    vmax_weath_mineral_act_sl(:,:)  = 0.0_wp
    k_adsorpt_po4_act(:,:)          = 0.0_wp
    k_desorpt_po4_act(:,:)          = 0.0_wp
    km_dip_biochem_po4_act(:,:)     = 0.0_wp


    !> run calc
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        km1_uptake_nh4_act(ic, is) = km1_up_nh4 / (rtm_hsc(ic, is) * rmm_hsc(ic, is))
        km1_uptake_no3_act(ic, is) = km1_up_no3 / (rtm_hsc(ic, is) * rmm_hsc(ic, is))
        km1_uptake_po4_act(ic, is) = km1_up_po4 / (rtm_hsc(ic, is) * rmm_hsc(ic, is)) &
          &                          / 1000._wp / (w_soil_sl(ic, is) / soil_depth_sl(ic, is))
        km2_uptake_nh4_act(ic, is) = km2_up_nh4 * rtm_hsc(ic, is) * rmm_hsc(ic, is)
        km2_uptake_no3_act(ic, is) = km2_up_no3 * rtm_hsc(ic, is) * rmm_hsc(ic, is)
        km2_uptake_po4_act(ic, is) = km2_up_po4 * rtm_hsc(ic, is) * rmm_hsc(ic, is) &
          &                          * 1000._wp * (w_soil_sl(ic, is) / soil_depth_sl(ic, is))
        km_uptake_org_act(ic, is)          = km_up_org  * rtm_hsc(ic, is) * rmm_hsc(ic, is)
        km_mic_uptake_nh4_act(ic, is)      = km_mic_uptake_nh4 / (rtm_hsc(ic, is) * rmm_hsc(ic, is))
        km_mic_uptake_no3_act(ic, is)      = km_mic_uptake_no3 / (rtm_hsc(ic, is) * rmm_hsc(ic, is))
        km_dip_biochem_po4_act(ic, is)     = km_dip_biochem_po4 * rtm_hsc(ic, is) * rmm_hsc(ic, is)
        km_denitrification_c_act(ic, is)   = km_denitrification_c
        km_denitrification_no3_act(ic, is) = km_denitrification_no3
        apparent_fast_som_cn(ic, is) = MAX(k_fast_som_cn_max - f_fast_som_cn * nh4_solute(ic, is) / bulk_dens_sl(ic, is), &
          &                         k_fast_som_cn_min)
        apparent_fast_som_np(ic, is) = k_fast_som_np &
          &                            * MAX(k_fast_som_cn_max - f_fast_som_cn * kref_fast_som_np / bulk_dens_sl(ic, is), &
          &                            k_fast_som_cn_min) / apparent_fast_som_cn(ic, is)
        vmax_weath_mineral_act_sl(ic, is) = rtm_hsc(ic, is) / rmm_hsc(ic, is) * vmax_weath_mineral_sl(ic, is)
        k_adsorpt_po4_act(ic, is) = rtm_sorption(ic, is) * k_adsorpt_po4
        k_desorpt_po4_act(ic, is) = rtm_desorption(ic, is) * k_desorpt_po4
        p_deposition_act(ic, is)  = 0.0_wp
        ! only for the 1st soil layer
        p_deposition_act(ic,1)  = p_deposition(ic) / 1.e6_wp / soil_depth_sl(ic,1) * dtime
      ENDDO
    ENDDO
  END SUBROUTINE calc_local_soil_profile_values

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates potential plant uptake of NH4, NO3 and PO4, given soil concentrations,
  !! and fine root mass per layer
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_plant_uptake_rate(nc, &                         ! in
                                              nsoil_sb, &
                                              dtime, &
                                              num_sl_above_bedrock, &
                                              lctlib_vmax_uptake_n, &
                                              lctlib_vmax_uptake_p, &
                                              flag_mycorrhiza, &
                                              soil_depth_sl, &
                                              veg_pool_fine_root_carbon, &
                                              sb_pool_mycorrhiza_carbon, &
                                              root_fraction_sl, &
                                              rtm_root_uptake_act, &
                                              rmm_root_uptake_act, &
                                              km1_uptake_nh4_act, &
                                              km2_uptake_nh4_act, &
                                              km1_uptake_no3_act, &
                                              km2_uptake_no3_act, &
                                              km1_uptake_po4_act, &
                                              km2_uptake_po4_act, &
                                              nh4_solute, &
                                              no3_solute, &
                                              po4_solute, &
                                              plant_uptake_nh4_sl, &        ! out
                                              plant_uptake_no3_sl, &
                                              plant_uptake_po4_sl, &
                                              unit_uptake_n_pot, &
                                              unit_uptake_p_pot  )

    USE mo_sb_constants,            ONLY: f_surface_myc2fineroot, myc_root_coverage_max

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                     INTENT(in)  :: nc, &                      !< dimensions
                                                                nsoil_sb
    REAL(wp),                                    INTENT(in)  :: dtime                      !< timestep length
    REAL(wp), DIMENSION(nc),                     INTENT(in)  :: num_sl_above_bedrock       !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp),                                    INTENT(in)  :: lctlib_vmax_uptake_n       !< lctlib parameter
    REAL(wp),                                    INTENT(in)  :: lctlib_vmax_uptake_p       !< lctlib parameter
    LOGICAL,                                     INTENT(in)  :: flag_mycorrhiza            !< from sb config
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(in)  :: soil_depth_sl              !< soil depth per layer [m]
    REAL(wp), DIMENSION(nc),                     INTENT(in)  :: veg_pool_fine_root_carbon  !< fine root carbon [mol / m2]
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(in)  :: sb_pool_mycorrhiza_carbon  !< mycorrhiza carbon [mol / m2]
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(in)  :: root_fraction_sl, &        !< fraction of fine roots per soil layer
                                                                rtm_root_uptake_act, &     !< temperature response of root uptake
                                                                rmm_root_uptake_act, &     !< moisture response of root uptake
                                                                km1_uptake_nh4_act, &      !< actual km1 for NH4 uptake [m3 / mol]
                                                                km2_uptake_nh4_act, &      !< actual km2 for NH4 uptake [mol / m3]
                                                                km1_uptake_no3_act, &      !< actual km1 for NO3 uptake [m3 / mol]
                                                                km2_uptake_no3_act, &      !< actual km2 for NO3 uptake [mol / m3]
                                                                km1_uptake_po4_act, &      !< actual km1 for PO4 uptake [m3 / mol]
                                                                km2_uptake_po4_act, &      !< actual km2 for PO4 uptake [mol / m3]
                                                                nh4_solute, &              !< NH4 solute concentration [mol / m3]
                                                                no3_solute, &              !< NO3 solute concentration [mol / m3]
                                                                po4_solute                 !< PO4 solute concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(out) :: &
                                                                plant_uptake_nh4_sl, &     !< potential plant uptake of NH4 [mol / m3 / timestep]
                                                                plant_uptake_no3_sl, &     !< potential plant uptake of NO3 [mol / m3 / timestep]
                                                                plant_uptake_po4_sl        !< potential plant uptake of PO4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc),                     INTENT(out) :: unit_uptake_n_pot, &       !< potential N uptake by fine roots per unit mass [micro-mol / mol / s]
                                                                unit_uptake_p_pot          !< potential P uptake by fine roots per unit mass [micro-mol / mol / s]
    ! ---------------------------
    ! 0.2 Local
    INTEGER                                          :: ic, is       !< loop over dimensions
    REAL(wp), DIMENSION(nc)                          :: myc_root_coverage_act, f_vmax_uptake_act
    REAL(wp), DIMENSION(nc,nsoil_sb)                 :: hlp1, hlp2   !< helpers
    CHARACTER(len=*), PARAMETER                      :: routine = TRIM(modname)//':calc_potential_plant_uptake_rate'



    !> 0.9 zero output variables
    !!
    plant_uptake_nh4_sl(:,:) = 0.0_wp
    plant_uptake_no3_sl(:,:) = 0.0_wp
    plant_uptake_po4_sl(:,:) = 0.0_wp
    unit_uptake_n_pot(:)     = 0.0_wp
    unit_uptake_p_pot(:)     = 0.0_wp


    !> 1.0 inorganic plant N uptake
    !!
    !! only root tips without mycorrhiza are allowed to take up N // max 30% of roots are infected! 'longer' myco if limited
    IF (flag_mycorrhiza) THEN
      DO ic = 1,nc
        IF (veg_pool_fine_root_carbon(ic) > eps8) THEN
          ! myc_root_coverage_act(ic) = myc_root_coverage_max
          ! dynamical root coverage to account for explorative mycorrhizae (idea for later...)
          myc_root_coverage_act(ic) = MIN(myc_root_coverage_max, &
            &                         SUM(sb_pool_mycorrhiza_carbon(ic,:) * soil_depth_sl(ic,:)) &
            &                         / veg_pool_fine_root_carbon(ic))
          f_vmax_uptake_act(ic)     = 1.0_wp / (1.0_wp - myc_root_coverage_max + myc_root_coverage_max * f_surface_myc2fineroot)
        ELSE
          myc_root_coverage_act(ic) = 0.0_wp
          f_vmax_uptake_act(ic)     = 1.0_wp
        ENDIF
      ENDDO
    ELSE
       f_vmax_uptake_act(:)     = 1.0_wp
       myc_root_coverage_act(:) = 0.0_wp
    ENDIF

    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        !> 2.1 inorganic plant N uptake
        !!
        !! potential N uptake rate (micro-mol N / mol C fine-root / s), weighted by root distribution
        hlp1(ic, is) = root_fraction_sl(ic, is) * rtm_root_uptake_act(ic, is) * rmm_root_uptake_act(ic, is) * &
                        (1.0_wp - myc_root_coverage_act(ic)) * lctlib_vmax_uptake_n * f_vmax_uptake_act(ic) * &
                        nh4_solute(ic, is) * ( km1_uptake_nh4_act(ic, is) + &
                        1._wp / ( km2_uptake_nh4_act(ic, is) + nh4_solute(ic, is) ))
        hlp2(ic, is) = root_fraction_sl(ic, is) * rtm_root_uptake_act(ic, is) * rmm_root_uptake_act(ic, is) * &
                        (1.0_wp - myc_root_coverage_act(ic)) * lctlib_vmax_uptake_n * f_vmax_uptake_act(ic) * &
                        no3_solute(ic, is) * ( km1_uptake_no3_act(ic, is) + &
                        1._wp / ( km2_uptake_no3_act(ic, is) + no3_solute(ic, is) ))
        unit_uptake_n_pot(ic) = unit_uptake_n_pot(ic) + hlp1(ic, is) + hlp2(ic, is)

        ! potential uptake rate in mol N / m3 / timestep
        plant_uptake_nh4_sl(ic, is) = hlp1(ic, is) * &
                                       veg_pool_fine_root_carbon(ic) / &
                                       soil_depth_sl(ic, is) * dtime / 1.e6_wp
        plant_uptake_no3_sl(ic, is) = hlp2(ic, is) * &
                                       veg_pool_fine_root_carbon(ic) / &
                                       soil_depth_sl(ic, is) * dtime / 1.e6_wp

        !> 2.2 inorganic plant P uptake
        !!
        !! potential P uptake rate (micro-mol P / mol C fine-root / s), weighted by root distribution
        hlp1(ic, is) = root_fraction_sl(ic, is)  * rtm_root_uptake_act(ic, is) * rmm_root_uptake_act(ic, is) * &
                        (1.0_wp - myc_root_coverage_act(ic)) * lctlib_vmax_uptake_p * f_vmax_uptake_act(ic) * &
                         po4_solute(ic, is) * ( km1_uptake_po4_act(ic, is) + &
                         1._wp / ( km2_uptake_po4_act(ic, is) + po4_solute(ic, is) ))
        unit_uptake_p_pot(ic) = unit_uptake_p_pot(ic) + hlp1(ic, is)
        ! potential uptake rate in mol P / m3 / timestep
        plant_uptake_po4_sl(ic, is) = hlp1(ic, is) * veg_pool_fine_root_carbon(ic) / soil_depth_sl(ic, is) * dtime / 1.e6_wp
      ENDDO
    ENDDO
  END SUBROUTINE calc_potential_plant_uptake_rate

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates potential mycorrhiza uptake of NH4, NO3 and PO4, given soil concentrations,
  !! and mycorrhiza mass per layer
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_mycorrhiza_uptake_rate(nc, &                         ! in
                                                   nsoil_sb, &
                                                   dtime, &
                                                   lctlib_vmax_uptake_n, &
                                                   lctlib_vmax_uptake_p, &
                                                   soil_depth_sl, &
                                                   veg_pool_fine_root_carbon, &
                                                   sb_pool_mycorrhiza_carbon, &
                                                   sb_pool_mycorrhiza_nitrogen, &
                                                   rtm_mycorrhiza_uptake_act, &
                                                   rmm_mycorrhiza_uptake_act, &
                                                   km1_uptake_nh4_act, &
                                                   km2_uptake_nh4_act, &
                                                   km1_uptake_no3_act, &
                                                   km2_uptake_no3_act, &
                                                   km1_uptake_po4_act, &
                                                   km2_uptake_po4_act, &
                                                   nh4_solute, &
                                                   no3_solute, &
                                                   po4_solute, &
                                                   mycorrhiza_uptake_nh4_sl, &   ! out
                                                   mycorrhiza_uptake_no3_sl, &
                                                   mycorrhiza_uptake_po4_sl )

    USE mo_veg_constants,           ONLY: transform_cost_nh4, transform_cost_no3
    USE mo_sb_constants,            ONLY: f_maxresp_myc, f_surface_myc2fineroot, myc_root_coverage_max

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)  :: nc                           !< dimensions
    INTEGER,                                          INTENT(in)  :: nsoil_sb                     !< number of soil layers
    REAL(wp),                                         INTENT(in)  :: dtime                        !< timestep length
    REAL(wp),                                         INTENT(in)  :: lctlib_vmax_uptake_n         !< lctlib parameter
    REAL(wp),                                         INTENT(in)  :: lctlib_vmax_uptake_p         !< lctlib parameter
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: soil_depth_sl                !< soil depth per layer [m]
    REAL(wp), DIMENSION(nc),                          INTENT(in)  :: veg_pool_fine_root_carbon    !< fine root carbon [mol / m2]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: sb_pool_mycorrhiza_carbon    !< mycorrhizal biomass [mol C / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: sb_pool_mycorrhiza_nitrogen  !< mycorrhizal biomass [mol C / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: rtm_mycorrhiza_uptake_act, & !< temperature response of mycorrhiza uptake
                                                                     rmm_mycorrhiza_uptake_act, & !< moisture response of mycorrhiza uptake
                                                                     km1_uptake_nh4_act, &        !< actual km1 for NH4 uptake [m3 / mol]
                                                                     km2_uptake_nh4_act, &        !< actual km2 for NH4 uptake [mol / m3]
                                                                     km1_uptake_no3_act, &        !< actual km1 for NO3 uptake [m3 / mol]
                                                                     km2_uptake_no3_act, &        !< actual km2 for NO3 uptake [mol / m3]
                                                                     km1_uptake_po4_act, &        !< actual km1 for PO4 uptake [m3 / mol]
                                                                     km2_uptake_po4_act, &        !< actual km2 for PO4 uptake [mol / m3]
                                                                     nh4_solute, &                !< NH4 solute concentration [mol / m3]
                                                                     no3_solute, &                !< NO3 solute concentration [mol / m3]
                                                                     po4_solute                   !< PO4 solute concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out) :: mycorrhiza_uptake_nh4_sl, &  !< potential mycorrhiza uptake of NH4 [mol / m3 / timestep]
                                                                     mycorrhiza_uptake_no3_sl, &  !< potential mycorrhiza uptake of NO3 [mol / m3 / timestep]
                                                                     mycorrhiza_uptake_po4_sl     !< potential mycorrhiza uptake of PO4 [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    INTEGER                                          :: isoil        !< loop over nsoil_sb
    REAL(wp), DIMENSION(nc)                          :: f_vmax_uptake_act
    REAL(wp), DIMENSION(nc,nsoil_sb)                 :: hlp1, hlp2   !< helpers
    CHARACTER(len=*), PARAMETER                      :: routine = TRIM(modname)//':calc_potential_mycorrhiza_uptake_rate'

    !> 0.9.1 init out arguments
    !>
    mycorrhiza_uptake_nh4_sl(:,:) = 0.0_wp
    mycorrhiza_uptake_no3_sl(:,:) = 0.0_wp
    mycorrhiza_uptake_po4_sl(:,:) = 0.0_wp

    !> 0.9.2 init local var
    !>
    f_vmax_uptake_act(:) = 0.0_wp

    !> 1.0 pot. N uptake by mycorrhizal fungi (mol / m3 / timestep)
    !>
    WHERE(veg_pool_fine_root_carbon(:) > eps8)
       f_vmax_uptake_act(:) =  1.0_wp / (1.0_wp - myc_root_coverage_max + myc_root_coverage_max * f_surface_myc2fineroot)
    ELSEWHERE
       f_vmax_uptake_act(:) = 1.0_wp
    ENDWHERE
    DO isoil = 1,nsoil_sb
      mycorrhiza_uptake_nh4_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          rtm_mycorrhiza_uptake_act(:,isoil) * rmm_mycorrhiza_uptake_act(:,isoil) * &
                                          f_surface_myc2fineroot * lctlib_vmax_uptake_n * f_vmax_uptake_act(:) * &
                                          nh4_solute(:,isoil) * ( km1_uptake_nh4_act(:,isoil) + &
                                          1._wp / ( km2_uptake_nh4_act(:,isoil) + nh4_solute(:,isoil) )) * &
                                          dtime / 1.e6_wp
      mycorrhiza_uptake_no3_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          rtm_mycorrhiza_uptake_act(:,isoil) * rmm_mycorrhiza_uptake_act(:,isoil) * &
                                          f_surface_myc2fineroot * lctlib_vmax_uptake_n * f_vmax_uptake_act(:) * &
                                          no3_solute(:,isoil) * ( km1_uptake_no3_act(:,isoil) + &
                                          1._wp / ( km2_uptake_no3_act(:,isoil) + no3_solute(:,isoil) )) * &
                                          dtime / 1.e6_wp
      !> 1.1 constrain potential N uptake by avaiable C for N transformation to organic
      !!
      hlp1(:,isoil) = mycorrhiza_uptake_nh4_sl(:,isoil) * transform_cost_nh4 + &
                      mycorrhiza_uptake_no3_sl(:,isoil) * transform_cost_no3
      WHERE(hlp1(:,isoil)  > (sb_pool_mycorrhiza_carbon(:,isoil) * f_maxresp_myc / one_day * dtime))
        WHERE(hlp1(:,isoil) > eps8)
          hlp2(:,isoil)                     = sb_pool_mycorrhiza_carbon(:,isoil) * f_maxresp_myc / one_day * dtime / &
                                              (hlp1(:,isoil))
          mycorrhiza_uptake_nh4_sl(:,isoil) = hlp2(:,isoil) * mycorrhiza_uptake_nh4_sl(:,isoil)
          mycorrhiza_uptake_no3_sl(:,isoil) = hlp2(:,isoil) * mycorrhiza_uptake_no3_sl(:,isoil)
        ELSEWHERE
          mycorrhiza_uptake_nh4_sl(:,isoil) = 0.0_wp
          mycorrhiza_uptake_no3_sl(:,isoil) = 0.0_wp
        ENDWHERE
      ENDWHERE
      !> 2.0 pot. P uptake by mycorrhizal fungi (mol / m3 / timestep)
      !!
      mycorrhiza_uptake_po4_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          rtm_mycorrhiza_uptake_act(:,isoil) * rmm_mycorrhiza_uptake_act(:,isoil) *  &
                                          f_surface_myc2fineroot * lctlib_vmax_uptake_p * f_vmax_uptake_act(:) * &
                                          po4_solute(:,isoil) * ( km1_uptake_po4_act(:,isoil) + &
                                          1._wp / ( km2_uptake_po4_act(:,isoil) + po4_solute(:,isoil) )) * &
                                          dtime / 1.e6_wp
    END DO
  END SUBROUTINE calc_potential_mycorrhiza_uptake_rate


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the potential rates of nitrification and denitrification
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_nitrification_denitrification_rate(nc, &                            ! in
                                                               nsoil_sb, &
                                                               dtime, &
                                                               num_sl_above_bedrock, &
                                                               rtm_nitrification, &
                                                               rmm_nitrification, &
                                                               rtm_denitrification, &
                                                               km_denitrification_c_act, &
                                                               km_denitrification_no3_act, &
                                                               anaerobic_volume_fraction, &
                                                               nh4_solute, &
                                                               no3_solute, &
                                                               sb_pool_microbial_carbon, &
                                                               nitrification_no3, &               ! out
                                                               denitrification_n2)

    USE mo_sb_constants,     ONLY: vmax_nitrification, vmax_denitrification

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)  :: nc                             !< dimensions
    INTEGER,                                          INTENT(in)  :: nsoil_sb                       !< number of soil layers
    REAL(wp),                                         INTENT(in)  :: dtime                          !< timestep length
    REAL(wp), DIMENSION(nc),                          INTENT(in)  :: num_sl_above_bedrock           !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: rtm_nitrification, &           !< actual temperature response for nitrification
                                                                     rmm_nitrification, &           !< actual moisture response for nitrification
                                                                     rtm_denitrification, &         !< actual temperature response for denitrification
                                                                     km_denitrification_c_act, &    !< actual km on C for denitrification
                                                                     km_denitrification_no3_act, &  !< actual km on NO3 for denitrification
                                                                     anaerobic_volume_fraction, &   !< actual anaerobic volume fraction
                                                                     nh4_solute, &                  !< NH4 in solution [mol / m3]
                                                                     no3_solute                     !< NO3 in solution [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: sb_pool_microbial_carbon       !< microbial C [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out) :: nitrification_no3, &           !< potential nitrification rate [mol / m3 / timestep]
                                                                     denitrification_n2             !< potential denitrification rate [mol / m3 / timestep]
    ! 0.2 Local
    INTEGER                     :: ic, is       !< loop over dimensions
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_nitrification_denitrification_rate'

    !> init out var
    !>
    nitrification_no3(:,:)  = 0.0_wp
    denitrification_n2(:,:) = 0.0_wp

    !> run calc
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        ! potential nitrification rate from NH4 to NO3 [mol/m3/timestep]
        nitrification_no3(ic,is)   = vmax_nitrification * rtm_nitrification(ic,is) * rmm_nitrification(ic,is) * &
                                   sb_pool_microbial_carbon(ic,is) / &
                                   (km_denitrification_c_act(ic,is) + sb_pool_microbial_carbon(ic,is))  * &
                                   (1._wp - anaerobic_volume_fraction(ic,is)) * nh4_solute(ic,is) * dtime / 1.e6_wp

        ! potential denitrification rate from NO3 to N2 [mol/m3/timestep]
        denitrification_n2(ic,is)  = vmax_denitrification * rtm_denitrification(ic,is) * &
                                   sb_pool_microbial_carbon(ic,is) / &
                                   (km_denitrification_c_act(ic,is) + sb_pool_microbial_carbon(ic,is))  * &
                                   anaerobic_volume_fraction(ic,is) * no3_solute(ic,is) / &
                                   (km_denitrification_no3_act(ic,is) + no3_solute(ic,is)) * dtime / 1.e6_wp
      ENDDO
    ENDDO
  END SUBROUTINE calc_potential_nitrification_denitrification_rate

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the leaching and gas transport within the soil column,
  !! given current soil N and P concentrations, soil temperature, moisture and water percolation rate
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_inorganic_transport_rate(nc, &                              ! in
                                           nsoil_sb, &
                                           dtime, &
                                           num_sl_above_bedrock, &
                                           soil_depth_sl, &
                                           percolation_sl, &
                                           frac_w_lat_loss_sl, &
                                           rtm_gasdiffusion, &
                                           rmm_gasdiffusion, &
                                           nh4_solute, &
                                           nh4_n15_solute, &
                                           no3_solute, &
                                           no3_n15_solute, &
                                           po4_solute, &
                                           noy, &
                                           noy_n15, &
                                           n2o, &
                                           n2o_n15, &
                                           n2, &
                                           n2_n15, &
                                           transport_nh4_solute, &            ! out
                                           transport_nh4_n15_solute, &
                                           transport_no3_solute, &
                                           transport_no3_n15_solute, &
                                           transport_po4_solute, &
                                           transport_noy, &
                                           transport_noy_n15, &
                                           transport_n2o, &
                                           transport_n2o_n15, &
                                           transport_n2, &
                                           transport_n2_n15, &
                                           leaching_nh4_solute, &
                                           leaching_nh4_n15_solute, &
                                           leaching_no3_solute, &
                                           leaching_no3_n15_solute, &
                                           leaching_po4_solute, &
                                           lateral_loss_nh4_solute_sl, &
                                           lateral_loss_nh4_n15_solute_sl, &
                                           lateral_loss_no3_solute_sl, &
                                           lateral_loss_no3_n15_solute_sl, &
                                           lateral_loss_po4_solute_sl, &
                                           emission_noy, &
                                           emission_noy_n15, &
                                           emission_n2o, &
                                           emission_n2o_n15, &
                                           emission_n2, &
                                           emission_n2_n15)

    USE mo_sb_constants,       ONLY: fleach_nh4, fleach_po4
    USE mo_q_sb_jsm_transport, ONLY: calc_liquid_phase_transport

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                           INTENT(in)    :: nc                                 !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                           !< number of soil layers
    REAL(wp),                          INTENT(in)    :: dtime                              !< timestep length
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock               !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: soil_depth_sl                      !< soil layer depth [m]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: percolation_sl, &                  !< fraction of water lost for 1d model [s-1]
                                                        frac_w_lat_loss_sl, &              !< constrained fraction of lateral (horizontal) water loss of 'w_soil_sl_old' (prev. timestep)
                                                        rtm_gasdiffusion, &                !< temperature modifier for gas diffusion
                                                        rmm_gasdiffusion, &                !< moisture modifier for gas diffusion
                                                        nh4_solute, &                      !< NH4 in solution [mol / m3]
                                                        nh4_n15_solute, &                  !< 15NH4 in solution [mol / m3]
                                                        no3_solute, &                      !< NO3 in solution [mol / m3]
                                                        no3_n15_solute, &                  !< 15NO3 in solution [mol / m3]
                                                        po4_solute, &                      !< PO4 in solution [mol / m3]
                                                        noy, &                             !< NOy concentration [mol / m3]
                                                        noy_n15, &                         !< 15NOy concentration [mol / m3]
                                                        n2o, &                             !< N2O concentration [mol / m3]
                                                        n2o_n15, &                         !< 15N2O concentration [mol / m3]
                                                        n2, &                              !< N2 concentration [mol / m3]
                                                        n2_n15                             !< 15N2 concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_nh4_solute, &            !< vertical transport of NH4 inside soil column [mol / m3 / timestep]
                                                        transport_nh4_n15_solute, &        !< vertical transport of 15NH4 inside soil column [mol / m3 / timestep]
                                                        transport_no3_solute, &            !< vertical transport of NO3 inside soil column [mol / m3 / timestep]
                                                        transport_no3_n15_solute, &        !< vertical transport of 15NO3 inside soil column [mol / m3 / timestep]
                                                        transport_po4_solute               !< vertical transport of PO4 inside soil column [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: transport_noy, &                   !< vertical transport of NOy [mol / m3 / timestep]
                                                        transport_noy_n15, &               !< vertical transport of 15NOy [mol / m3 / timestep]
                                                        transport_n2o, &                   !< vertical transport of N2O [mol / m3 / timestep]
                                                        transport_n2o_n15, &               !< vertical transport of 15N2O [mol / m3 / timestep]
                                                        transport_n2, &                    !< vertical transport of N2 [mol / m3 / timestep]
                                                        transport_n2_n15, &                !< vertical transport of 15N2 [mol / m3 / timestep]
                                                        lateral_loss_nh4_solute_sl, &      !< lateral loss of NH4 [mol/ m2 / timestep]
                                                        lateral_loss_nh4_n15_solute_sl, &  !< lateral loss of 15NH4 [mol / m2 / timestep]
                                                        lateral_loss_no3_solute_sl, &      !< lateral loss of NO3 [mol/ m2 / timestep]
                                                        lateral_loss_no3_n15_solute_sl, &  !< lateral loss of 15NO3 [mol / m2 / timestep]
                                                        lateral_loss_po4_solute_sl         !< lateral loss of PO4 [mol / m2 / timestep]
    REAL(wp), DIMENSION(nc),           INTENT(out)   :: leaching_nh4_solute, &             !< leaching of NH4 [mol / m2 / timestep]
                                                        leaching_nh4_n15_solute, &         !< leaching of 15NH4 [mol / m2 / timestep]
                                                        leaching_no3_solute, &             !< leaching of NO3 [mol / m2 / timestep]
                                                        leaching_no3_n15_solute, &         !< leaching of 15NO3 [mol / m2 / timestep]
                                                        leaching_po4_solute, &             !< leaching of PO4 [mol / m2 / timestep]
                                                        emission_noy, &                    !< soil efflux of NOy [mol / m2 / timestep]
                                                        emission_noy_n15, &                !< soil efflux of 15NOy [mol / m2 / timestep]
                                                        emission_n2o, &                    !< soil efflux of N2O [mol / m2 / timestep]
                                                        emission_n2o_n15, &                !< soil efflux of 15N2O [mol / m2 / timestep]
                                                        emission_n2, &                     !< soil efflux of N2 [mol / m2 / timestep]
                                                        emission_n2_n15                    !< soil efflux of 15N2 [mol / m2 / timestep]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_inorganic_transport_rate'

    !> init out var
    !>
    transport_noy(:,:)                  = 0.0_wp
    transport_noy_n15(:,:)              = 0.0_wp
    transport_n2o(:,:)                  = 0.0_wp
    transport_n2o_n15(:,:)              = 0.0_wp
    transport_n2(:,:)                   = 0.0_wp
    transport_n2_n15(:,:)               = 0.0_wp
    lateral_loss_nh4_solute_sl(:,:)     = 0.0_wp
    lateral_loss_nh4_n15_solute_sl(:,:) = 0.0_wp
    lateral_loss_no3_solute_sl(:,:)     = 0.0_wp
    lateral_loss_no3_n15_solute_sl(:,:) = 0.0_wp
    lateral_loss_po4_solute_sl(:,:)     = 0.0_wp
    leaching_nh4_solute(:)              = 0.0_wp
    leaching_nh4_n15_solute(:)          = 0.0_wp
    leaching_no3_solute(:)              = 0.0_wp
    leaching_no3_n15_solute(:)          = 0.0_wp
    leaching_po4_solute(:)              = 0.0_wp
    emission_noy(:)                     = 0.0_wp
    emission_noy_n15(:)                 = 0.0_wp
    emission_n2o(:)                     = 0.0_wp
    emission_n2o_n15(:)                 = 0.0_wp
    emission_n2(:)                      = 0.0_wp
    emission_n2_n15(:)                  = 0.0_wp


    !>1.0. Leaching rate
    !>
    ! For the 1D scheme with multiple soil hydrology layers, use adjection-diffusion scheme of JSM
    ! reduction of layer concentration in mol m-3 timestep-1, leaching rate in mol m-2 timestep-1
    CALL calc_liquid_phase_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
                                     percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
                                     nh4_solute(:,:), transport_nh4_solute(:,:), leaching_nh4_solute(:), &
                                     lateral_loss_nh4_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
                                     percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
                                     nh4_n15_solute(:,:), transport_nh4_n15_solute(:,:), leaching_nh4_n15_solute(:), &
                                     lateral_loss_nh4_n15_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
                                     percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
                                     no3_solute(:,:), transport_no3_solute(:,:), leaching_no3_solute(:), &
                                     lateral_loss_no3_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
                                     percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
                                     no3_n15_solute(:,:), transport_no3_n15_solute(:,:), leaching_no3_n15_solute(:), &
                                     lateral_loss_no3_n15_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
                                     percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
                                     po4_solute(:,:), transport_po4_solute(:,:), leaching_po4_solute(:), &
                                     lateral_loss_po4_solute_sl(:,:))

    !> 2.0 gaseous loss of NOy, N2O, and N2, for now following Xu-Ri and Prentice 2008,
    !!thus ignoring any actual gas-transport and instead assume that from each layer a certain fraction excapes directly to the atmosphere.
    !! [mol m-3 timestep-1]
    transport_noy(:,:)     = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * noy(:,:)     / one_day * dtime
    transport_noy_n15(:,:) = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * noy_n15(:,:) / one_day * dtime
    transport_n2o(:,:)     = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * n2o(:,:)     / one_day * dtime
    transport_n2o_n15(:,:) = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * n2o_n15(:,:) / one_day * dtime
    transport_n2(:,:)      = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * n2(:,:)      / one_day * dtime
    transport_n2_n15(:,:)  = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * n2_n15(:,:)  / one_day * dtime

    ! emission rates as positive flux to the atmosphere, integrated over the soil column [mol m-2 timestep-1]
    emission_noy(:)     = SUM(-transport_noy(:,:)     * soil_depth_sl(:,:), DIM=2)
    emission_n2o(:)     = SUM(-transport_n2o(:,:)     * soil_depth_sl(:,:), DIM=2)
    emission_n2(:)      = SUM(-transport_n2(:,:)      * soil_depth_sl(:,:), DIM=2)
    emission_noy_n15(:) = SUM(-transport_noy_n15(:,:) * soil_depth_sl(:,:), DIM=2)
    emission_n2o_n15(:) = SUM(-transport_n2o_n15(:,:) * soil_depth_sl(:,:), DIM=2)
    emission_n2_n15(:)  = SUM(-transport_n2_n15(:,:)  * soil_depth_sl(:,:), DIM=2)

  END SUBROUTINE calc_inorganic_transport_rate

  ! ======================================================================================================= !
  !>Calculate potential decomposition rate of SOM and litter pools according to:
  !>  `dX/dt = - f(T) * g(M) / tau_pool * X`
  !>
  SUBROUTINE calc_potential_decomposition_rate( &
    & nc, &
    & nsoil_sb, &
    & dtime, &
    & soil_depth_sl, &
    & elements_index_map, &
    & is_element_used, &
    & flag_mycorrhiza_prim, &
    & flag_mycorrhiza_org, &
    & rtm_decomposition, &
    & rmm_decomposition, &
    & sb_pool_mt, &
    & het_respiration, &
    & het_respiration_c13, &
    & het_respiration_c14, &
    & myc_respiration, &
    & myc_respiration_c13, &
    & myc_respiration_c14, &
    & sb_loss_mt )

    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_sb_constants,            ONLY: tau_metabolic_litter, tau_structural_litter, tau_woody_litter, tau_fast, tau_slow, &
      &                                   mycorrhiza_cn, vmax_priming, km_priming, resp_priming, f_tau_myc
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                            INTENT(in)    :: nc                       !< dimensions
    INTEGER,                            INTENT(in)    :: nsoil_sb                 !< number of soil layers
    REAL(wp),                           INTENT(in)    :: dtime                    !< timestep length
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: soil_depth_sl            !< soil layer depth [m]
    INTEGER,                            INTENT(in)    :: elements_index_map(:)    !< map bgcm element ID -> IDX
    LOGICAL,                            INTENT(in)    :: is_element_used(:)       !< is element in 'elements_index_map' used
    LOGICAL,                            INTENT(in)    :: flag_mycorrhiza_prim, &  !< from sb config
                                                         flag_mycorrhiza_org      !< from sb config
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: rtm_decomposition, &     !< temperature response of decomposition
                                                         rmm_decomposition        !< moisture response of decomposition
    REAL(wp),                           INTENT(in)    :: sb_pool_mt(:,:,:,:)      !< sb_pool [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: het_respiration, &       !< heterotrophic respiration [mol C / m3 / timestep]
                                                         het_respiration_c13, &   !< heterotrophic respiration [mol 13C / m3 / timestep]
                                                         het_respiration_c14, &   !< heterotrophic respiration [micro-mol 14C / m3 / timestep]
                                                         myc_respiration, &       !< mycorrhizae respiration [mol C / m3 / timestep]
                                                         myc_respiration_c13, &   !< mycorrhizae respiration [mol 13C / m3 / timestep]
                                                         myc_respiration_c14      !< mycorrhizae respiration [micro-mol 14C / m3 / timestep]
    REAL(wp),                           INTENT(inout) :: sb_loss_mt(:,:,:,:)      !< sb_loss flux [mol / m3  ??]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                           :: ielem                !< loop over bgcm elements
    INTEGER                           :: ix_elem              !< index of element in bgcm, used for looping
    REAL(wp), DIMENSION(nc,nsoil_sb)  :: tau_slow_myc, &
                                         km_priming_act, &
                                         resp_priming_act, &
                                         tau_slow_act
    REAL(wp), DIMENSION(nc,nsoil_sb)  :: hlp1, hlp2, hlp3
    REAL(wp), DIMENSION(nc,nsoil_sb)  :: hlp_cn
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_decomposition_rate'
    ! ----------------------------------------------------------------------------------------------------- !

    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm

        !>1.0 metabolic litter (= soluable pool)
        !>
        sb_loss_mt(ix_soluable_litter, ix_elem, :, :) = sb_loss_mt(ix_soluable_litter, ix_elem, :, :) &
          &                                             + rtm_decomposition(:,:) * rmm_decomposition(:,:) &
          &                                             / tau_metabolic_litter &
          &                                             * dtime * sb_pool_mt(ix_soluable_litter, ix_elem, :, :)

        !>2.0 structural litter (= polymeric pool)
        !>
        sb_loss_mt(ix_polymeric_litter, ix_elem, :, :) = sb_loss_mt(ix_polymeric_litter, ix_elem, :, :) &
          &                                              + rtm_decomposition(:,:) * rmm_decomposition(:,:) &
          &                                              / tau_structural_litter &
          &                                              * dtime * sb_pool_mt(ix_polymeric_litter, ix_elem, :, :)

        !>3.0 woody litter
        !>
        sb_loss_mt(ix_woody_litter, ix_elem, :, :) = sb_loss_mt(ix_woody_litter, ix_elem, :, :) &
          &                                          + rtm_decomposition(:,:) * rmm_decomposition(:,:) &
          &                                          / tau_woody_litter &
          &                                          * dtime * sb_pool_mt(ix_woody_litter, ix_elem, :, :)

        !>4.0 fast SOM pool (= microbial pool)
        !>
        sb_loss_mt(ix_microbial, ix_elem, :, :) = sb_loss_mt(ix_microbial, ix_elem, :, :) &
          &                                       + rtm_decomposition(:,:) * rmm_decomposition(:,:) &
          &                                       / tau_fast &
          &                                       * dtime * sb_pool_mt(ix_microbial, ix_elem, :, :)
      END IF  ! IF element is used
    END DO    ! loop over elements

    !>5.0 slow SOM pool (= residue pool)
    !>
    ! mycorrhiza_prim
    IF(flag_mycorrhiza_prim)THEN
      km_priming_act(:,:)   = km_priming
      resp_priming_act(:,:) = resp_priming
      WHERE(sb_pool_mt(ix_mycorrhiza, ixC, :, :) > eps8)
          tau_slow_act(:,:) = tau_slow / ((sb_pool_mt(ix_mycorrhiza, ixC, :, :) * vmax_priming) / &
                              (km_priming_act(:,:) + sb_pool_mt(ix_mycorrhiza, ixC, :, :)))
      ELSEWHERE
          hlp1(:,:)         = eps8
          tau_slow_act(:,:) = tau_slow /((hlp1(:,:) * vmax_priming) / (km_priming_act(:,:) + hlp1(:,:)))
      ENDWHERE

      sb_loss_mt(ix_mycorrhiza, ixC, :, :)   = sb_loss_mt(ix_mycorrhiza, ixC, :, :)   + resp_priming_act(:,:) &
        &                                      * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
      sb_loss_mt(ix_mycorrhiza, ixC13, :, :) = sb_loss_mt(ix_mycorrhiza, ixC13, :, :) + resp_priming_act(:,:) &
        &                                      * sb_pool_mt(ix_mycorrhiza, ixC13, :, :)
      sb_loss_mt(ix_mycorrhiza, ixC14, :, :) = sb_loss_mt(ix_mycorrhiza, ixC14, :, :) + resp_priming_act(:,:) &
        &                                      * sb_pool_mt(ix_mycorrhiza, ixC14, :, :)

      het_respiration(:,:)     = het_respiration(:,:)     + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
      het_respiration_c13(:,:) = het_respiration_c13(:,:) + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC13, :, :)
      het_respiration_c14(:,:) = het_respiration_c14(:,:) + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC14, :, :)

      myc_respiration(:,:)     = myc_respiration(:,:)     + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
      myc_respiration_c13(:,:) = myc_respiration_c13(:,:) + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC13, :, :)
      myc_respiration_c14(:,:) = myc_respiration_c14(:,:) + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC14, :, :)
    ! mycorrhiza_org
    ELSEIF(flag_mycorrhiza_org)THEN
      tau_slow_act(:,:) = tau_slow * f_tau_myc
    ! other
    ELSE
      tau_slow_act(:,:) = tau_slow
    END IF

    ! ...
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        sb_loss_mt(ix_residue, ix_elem, :, :) = sb_loss_mt(ix_residue, ix_elem, :, :) &
          &                                     + rtm_decomposition(:,:) * rmm_decomposition(:,:) &
          &                                     / tau_slow_act(:,:) &
          &                                     * dtime * sb_pool_mt(ix_residue, ix_elem, :, :)
      END IF  ! IF element is used
    END DO    ! loop over elements
  END SUBROUTINE calc_potential_decomposition_rate


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculate potential biomineralisation rate
  !!  dP/dt = - f(T)*g(M)/tau_pool * P
  !> Calculate weathering rate, slow adsorption rate and occlusion rate based on Yang et al. 2014
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_p_flux_rate(nc, &
                                        nsoil_sb, &
                                        dtime, &
                                        num_sl_above_bedrock, &
                                        flag_sb_prescribe_po4, &
                                        rtm_decomposition, &
                                        rmm_decomposition, &
                                        km_dip_biochem_po4_act, &
                                        bulk_dens_corr_sl, &
                                        fine_root_carbon_sl, &
                                        po4_solute, &
                                        po4_assoc_fast, &
                                        po4_assoc_slow, &
                                        po4_occluded, &
                                        po4_primary, &
                                        vmax_weath_mineral_act_sl, &
                                        k_adsorpt_po4_act, &
                                        k_desorpt_po4_act, &
                                        sb_pool_residue_phosphorus, &
                                        weathering_po4, &
                                        occlusion_po4, &
                                        slow_exchange_po4, &
                                        biochem_mineralisation_po4 )

    USE mo_sb_constants,            ONLY: tau_biomineralisation, k_occlude_po4, km_rootc_biochem_po4, km_rootc_weath_po4, &
      &                                   frac_background_weath_mineral

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                           INTENT(in)  :: nc                            !< dimensions
    INTEGER,                           INTENT(in)  :: nsoil_sb                      !< number of soil layers
    REAL(wp),                          INTENT(in)  :: dtime                         !< timestep length
    REAL(wp), DIMENSION(nc),           INTENT(in)  :: num_sl_above_bedrock          !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    LOGICAL,                           INTENT(in)  :: flag_sb_prescribe_po4         !< whether or not solute PO4 is prescribed
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)  :: rtm_decomposition, &          !< temperature response of decomposition
                                                      rmm_decomposition, &          !< moisture response of decomposition
                                                      km_dip_biochem_po4_act, &     !< actual half-saturation soluble P for biochem_mineralisation_po4 [mol P / m3]
                                                      bulk_dens_corr_sl, &          !< bulk density of the soil layer, [kg/m3]
                                                      fine_root_carbon_sl, &        !< carbon of fine roots per soil layer, [mol/m3]
                                                      po4_solute, &                 !< soluble PO4 pool [mol/m3]
                                                      po4_assoc_fast, &             !< fast minerally associated PO4 pool [mol/m3]
                                                      po4_assoc_slow, &             !< slow minerally associated PO4 pool [mol/m3]
                                                      po4_occluded, &               !< occluded PO4 pool [mol/m3]
                                                      po4_primary, &                !< primary PO4 pool [mol/m3]
                                                      vmax_weath_mineral_act_sl, &  !< actual weathering rate for primary P [mol P kg-1 s-1]
                                                      k_adsorpt_po4_act, &          !< actual adsorption rate for fast_assoc_po4 [s-1]
                                                      k_desorpt_po4_act, &          !< actual desorption rate for slow_assoc_po4 [s-1]
                                                      sb_pool_residue_phosphorus    !< slow SOM P pool [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out) :: weathering_po4, &             !< phosphorus weathering rate [mol P / m3 / timestep]
                                                      occlusion_po4, &              !< phosphorus occlusion rate, assoc_slow --> occluded [mol P / m3 / timestep]
                                                      slow_exchange_po4, &          !< phosphorus desorption rate, assoc_fast --> assoc_slow [mol P / m3 / timestep]
                                                      biochem_mineralisation_po4    !< potential biomineralisation rate [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    INTEGER                                          :: ic, is
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1, hlp2, hlp3
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_p_flux_rate'

    !> init out var
    !>
    weathering_po4(:,:)             = 0.0_wp
    occlusion_po4(:,:)              = 0.0_wp
    slow_exchange_po4(:,:)          = 0.0_wp
    biochem_mineralisation_po4(:,:) = 0.0_wp

    !> 1.0 Calculate P weathering rate from Primary P pool (mol /m3 / timestep)
    !! weathering rate is composed of a background rate (frac_background_weath_mineral) and pH/root biomass regulated rate
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        weathering_po4(ic, is) = vmax_weath_mineral_act_sl(ic, is) * po4_primary(ic, is) &
          &                   * (frac_background_weath_mineral + fine_root_carbon_sl(ic, is) &
          &                   / (km_rootc_weath_po4 + fine_root_carbon_sl(ic, is))) * dtime
        IF (weathering_po4(ic, is) > po4_primary(ic, is)) THEN
          weathering_po4(ic, is) = po4_primary(ic, is)
        ENDIF
        !> 1.1 Calculate P occlusion rate from occluded PO4 pool (mol /m3 / timestep)
        occlusion_po4(ic, is)     = po4_assoc_slow(ic, is) * k_occlude_po4 * dtime
        !> 1.2 Calculate P slow adsorption rate (mol /m3 / timestep)
        hlp1(ic, is)              = k_desorpt_po4_act(ic, is) * po4_assoc_slow(ic, is) * dtime
        hlp2(ic, is)              = k_adsorpt_po4_act(ic, is) * po4_assoc_fast(ic, is) * dtime
        slow_exchange_po4(ic, is) = hlp2(ic, is) - hlp1(ic, is)
      ENDDO
    ENDDO

    !> 2.0 Calculate potential P biomineralisation from slow organic P pool (mol /m3 / timestep)
    !!
    IF(.NOT. flag_sb_prescribe_po4) THEN
      DO ic = 1,nc
        DO is = 1,INT(num_sl_above_bedrock(ic))
          hlp1(ic, is) = km_dip_biochem_po4_act(ic, is) / (km_dip_biochem_po4_act(ic, is) + po4_solute(ic, is))
          IF (fine_root_carbon_sl(ic, is) > eps8) THEN
            hlp2(ic, is) = fine_root_carbon_sl(ic, is) / (km_rootc_biochem_po4 + fine_root_carbon_sl(ic, is))
          ELSE
            hlp2(ic, is) = 0.0_wp
          ENDIF
          biochem_mineralisation_po4(ic, is) = rtm_decomposition(ic, is) * rmm_decomposition(ic, is) / tau_biomineralisation * &
            &                                  sb_pool_residue_phosphorus(ic, is) * hlp1(ic, is) * hlp2(ic, is) * dtime
        ENDDO
      ENDDO
    ELSE
      biochem_mineralisation_po4(:,:) = 0.0_wp
    ENDIF

  END SUBROUTINE calc_potential_p_flux_rate

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates asymbiotic N fixation as a function of N deficit during litter decomposition
  !! and limiting functions from soil temperature, moisture and soil nutrient status
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_asymb_n_fixation(bnf_scheme, &
                                           sb_loss_litter_carbon, &
                                           sb_loss_litter_nitrogen, &
                                           apparent_fast_som_cn, &
                                           rtm_asymb_bnf, &
                                           rmm_asymb_bnf, &
                                           nh4_solute, &
                                           no3_solute) RESULT (asymb_n_fixation)

    USE mo_sb_constants,          ONLY: frac_litter2fast_som, microbial_nue, km_asymb_bnf, vmax_asymb_bnf

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*), INTENT(in)     :: bnf_scheme                    !< SB_ asymbiotic N fixation scheme
    REAL(wp), INTENT(in)             :: sb_loss_litter_carbon,  &     !< C flux from litter to fast pool [mol / m3 / timestep]
                                        sb_loss_litter_nitrogen,  &   !< N flux from litter to fast pool [mol / m3 / timestep]
                                        apparent_fast_som_cn, &       !< apparent microbial carbon-nitrogen ratio [mol/mol]
                                        rtm_asymb_bnf, &              !< temperature rate modifier []
                                        rmm_asymb_bnf, &              !< moisture rate modifier []
                                        nh4_solute, &                 !< NH4 in solution [mol / m3]
                                        no3_solute                    !< NO3 in solution [mol / m3]
    REAL(wp)                         :: asymb_n_fixation              !< asymbiotic N fixation [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                         :: f_nav,n_demand
    CHARACTER(len=*), PARAMETER      :: routine = TRIM(modname)//':calc_asymb_n_fixation'

    !> init out var
    !>
    asymb_n_fixation = 0.0_wp

    SELECT CASE(TRIM(bnf_scheme))
    ! for the SB_ asymbiotic N fixation no difference is made between "dynamic" and "fixed"
    !  ("fixed" is only available in the VEG_ symbiotic N fixation)
    CASE("dynamic", "fixed")
      !>1.0 N deficit of litter decomposition
      !>

      !>1.1 calculate potential transfer from litter decay to fast SOM pool, and associated potential microbial
      !>immobilisation rate as difference between maximal growth rate and N in litter
      !>
      n_demand = frac_litter2fast_som / apparent_fast_som_cn * sb_loss_litter_carbon - &
        &        microbial_nue * sb_loss_litter_nitrogen

      !>1.2 rate modifier for mineral soil nitrogen
      !>
      IF(nh4_solute + no3_solute > eps8) THEN
        f_nav = MIN(1._wp, 1.0_wp - (nh4_solute + no3_solute) / (nh4_solute + no3_solute + km_asymb_bnf))
      ELSE
        f_nav = 1._wp
      ENDIF

      asymb_n_fixation = vmax_asymb_bnf * rtm_asymb_bnf * rmm_asymb_bnf * f_nav * n_demand
    CASE("none", "unlimited")
      !>1.3 asymb_n_fixation is calculated during the decomposition process for these cases
      !>
      asymb_n_fixation = 0.0_wp
    END SELECT

  END FUNCTION calc_asymb_n_fixation

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the immobilisation rate of NH4 and PO4 of microbial biomass, given
  !! litter decay and stoichiometry, the co-occurring rates of plant, mycorrhizal uptake,
  !! nitrification, denitrification and transport, biomineralisation of PO4
  !! such that the newly formed fast SOM has the prescribed C:N:P stoichiometry; the carbon-
  !! use-efficiency of microbes is within the observed bound and mass-balance constraints are
  !! maintained. The result of this calculation are updates ("actual") rates of uptake, and processing
  !! as well as litter turnover rates.
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_actual_flux_rates( nc, &                         ! in
                                     nsoil_sb, &
                                     num_sl_above_bedrock, &
                                     bnf_scheme, &
                                     nh4_solute, &
                                     no3_solute, &
                                     po4_solute, &
                                     po4_assoc_fast, &
                                     transport_nh4_solute, &
                                     transport_no3_solute, &
                                     transport_po4_solute, &
                                     transport_po4_assoc_fast, &
                                     weathering_po4, &
                                     p_deposition_act, &
                                     slow_exchange_po4, &
                                     apparent_fast_som_cn , &
                                     km_mic_uptake_nh4_act, &
                                     km_mic_uptake_no3_act, &
                                     plant_uptake_nh4_sl, &                 ! inout
                                     plant_uptake_no3_sl, &
                                     plant_uptake_po4_sl, &
                                     mycorrhiza_uptake_nh4_sl, &
                                     mycorrhiza_uptake_no3_sl, &
                                     mycorrhiza_uptake_po4_sl, &
                                     nitrification_no3, &
                                     denitrification_n2, &
                                     asymb_n_fixation, &
                                     biochem_mineralisation_po4, &
                                     sb_loss_mt, &                          ! inout
                                     microbial_uptake_nh4_sl, &             ! out
                                     microbial_uptake_no3_sl, &             ! out
                                     microbial_uptake_po4_sl)               ! out

    USE mo_sb_constants,            ONLY: frac_litter2fast_som, k_fast_som_np, microbial_nue, microbial_pue, &
                                          vmax_mic_uptake_nh4, vmax_mic_uptake_no3

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                           INTENT(in)    :: nc, &                         !< dimensions
                                                        nsoil_sb
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock          !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    CHARACTER(len=*),                  INTENT(in)    :: bnf_scheme                    !< SB_ asymbiotic N fixation scheme
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: nh4_solute, &                 !< NH4 in solution [mol / m3]
                                                        no3_solute, &                 !< NO3 in solution [mol / m3]
                                                        po4_solute, &                 !< PO4 in solution [mol / m3]
                                                        po4_assoc_fast, &             !< labile mineral P [mol / m3]
                                                        transport_nh4_solute, &       !< transport of NH4 [mol / m3 / timestep]
                                                        transport_no3_solute, &       !< transport of NO3 [mol / m3 / timestep]
                                                        transport_po4_solute, &       !< transport of soluble PO4 [mol / m3 / timestep]
                                                        transport_po4_assoc_fast, &   !< transport of labile PO4 [mol / m3 / timestep]
                                                        weathering_po4, &             !< weathering of po4 [mol / m3 / timestep]
                                                        slow_exchange_po4, &          !< po4 desorption rate, assoc_fast --> assoc_slow [mol P / m3 / timestep]
                                                        p_deposition_act, &           !< actual P deposition rate, [mol / m3 / timestep]
                                                        apparent_fast_som_cn, &       !< apparent microbial carbon-nitrogen ratio [mol/mol]
                                                        km_mic_uptake_nh4_act, &      !< actual km for microbial NH4 uptake [mol / m3]
                                                        km_mic_uptake_no3_act         !< actual km for microbial NO3 uptake [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: plant_uptake_nh4_sl, &        !< plant NH4 uptake [mol / m3 / timestep]
                                                        plant_uptake_no3_sl, &        !< plant NO3 uptake [mol / m3 / timestep]
                                                        plant_uptake_po4_sl, &        !< plant PO4 uptake [mol / m3 / timestep]
                                                        mycorrhiza_uptake_nh4_sl, &   !< mycorrhizal NH4 uptake [mol / m3 / timestep]
                                                        mycorrhiza_uptake_no3_sl, &   !< mycorrhizal NO3 uptake [mol / m3 / timestep]
                                                        mycorrhiza_uptake_po4_sl, &   !< mycorrhizal PO4 uptake [mol / m3 / timestep]
                                                        nitrification_no3, &          !< nitrification rate [mol / m3 / timestep]
                                                        denitrification_n2, &         !< denitrification rate [mol / m3 / timestep]
                                                        asymb_n_fixation, &           !< asymbiotic N fixation (unlimited BNF only) [mol / m3 / timestep]
                                                        biochem_mineralisation_po4    !< biomineralisation rate [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_loss_mt(:,:,:,:)           !< bgcm sb_loss flux: decomposition of SOM and litter pools [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: microbial_uptake_nh4_sl, &    !< microbial immobilisation of NH4 [mol / m3 / timestep]
                                                        microbial_uptake_no3_sl, &    !< microbial immobilisation of NO3 [mol / m3 / timestep]
                                                        microbial_uptake_po4_sl       !< microbial immobilisation of PO4 [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    INTEGER                           :: ic, is      !< looping
    REAL(wp), DIMENSION(nc, nsoil_sb) :: hlp1, hlp2
    REAL(wp), DIMENSION(nc, nsoil_sb) :: transfer_c
    REAL(wp), DIMENSION(nc, nsoil_sb) :: f_uptake_nh4_sl
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_actual_flux_rates'

    !> init out var
    !>
    microbial_uptake_nh4_sl(:,:) = 0.0_wp
    microbial_uptake_no3_sl(:,:) = 0.0_wp
    microbial_uptake_po4_sl(:,:) = 0.0_wp

    !> 1.0 formation of fast SOM from litter
    !!

    !> 1.1 calculate potential transfer from litter decay to fast SOM pool, and associated potential microbial
    !!     immobilisation rate as difference between maximal growth rate and N in litter
    transfer_c(:,:) = frac_litter2fast_som * &
                      (sb_loss_mt(ix_soluable_litter, ixC, :, :) + &
                       sb_loss_mt(ix_polymeric_litter, ixC, :, :))

    !> 1.2 calculate partitioning of N immobilisation to NH4 and NO3 according to their abundance and the
    !!     affinities of the uptake channels to either species
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        IF (nh4_solute(ic, is) + no3_solute(ic, is) > eps8) THEN
          f_uptake_nh4_sl(ic, is) = vmax_mic_uptake_nh4 * nh4_solute(ic, is) &
            &                       / (nh4_solute(ic, is) + km_mic_uptake_nh4_act(ic, is)) &
            &                       / (vmax_mic_uptake_nh4 * nh4_solute(ic, is) / (nh4_solute(ic, is) &
            &                       + km_mic_uptake_nh4_act(ic, is)) + vmax_mic_uptake_no3 &
            &                       * no3_solute(ic, is) / (no3_solute(ic, is) + km_mic_uptake_no3_act(ic, is)))
        ELSE
          f_uptake_nh4_sl(ic, is) = 1.0_wp
        ENDIF

        microbial_uptake_nh4_sl(ic, is) = f_uptake_nh4_sl(ic, is) * &
                                                 (transfer_c(ic, is) / apparent_fast_som_cn(ic, is) - &
                                                  microbial_nue * (sb_loss_mt(ix_soluable_litter, ixN, ic, is) + &
                                                                  sb_loss_mt(ix_polymeric_litter, ixN, ic, is)) - &
                                                  asymb_n_fixation(ic, is))

        microbial_uptake_no3_sl(ic, is) = (1.0_wp - f_uptake_nh4_sl(ic, is)) * &
                                                 (transfer_c(ic, is) / apparent_fast_som_cn(ic, is) - &
                                                  microbial_nue * (sb_loss_mt(ix_soluable_litter, ixN, ic, is) + &
                                                                  sb_loss_mt(ix_polymeric_litter, ixN, ic, is)) - &
                                                  asymb_n_fixation(ic, is))

        !! limit uptake rate at low nitrogen assuming a Michaelis-Menten Kinetic
        hlp1(ic, is)                    = nh4_solute(ic, is) / ( km_mic_uptake_nh4_act(ic, is) + nh4_solute(ic, is) )
        microbial_uptake_nh4_sl(ic, is) = microbial_uptake_nh4_sl(ic, is) * MIN(hlp1(ic, is), 1.0_wp)

        hlp1(ic, is)                    = no3_solute(ic, is) / ( km_mic_uptake_no3_act(ic, is) + no3_solute(ic, is) )
        microbial_uptake_no3_sl(ic, is) = microbial_uptake_no3_sl(ic, is) * MIN(hlp1(ic, is), 1.0_wp)
      ENDDO
    ENDDO

    !>1.3  Limit microbial, mycorrhizal and plant uptake, as well as nitrification to available NH4 minus transport losses
    !>
    SELECT CASE(TRIM(bnf_scheme))
    CASE("unlimited")
      asymb_n_fixation(:,:) = MAX(microbial_uptake_nh4_sl(:,:) + plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:) + &
                              mycorrhiza_uptake_nh4_sl(:,:) - (nh4_solute(:,:) + transport_nh4_solute(:,:)), 0.0_wp)
    ! for the SB_ asymbiotic N fixation no difference is made between "dynamic" and "fixed"
    !  ("fixed" is only available in the VEG_ symbiotic N fixation)
    CASE("none", "dynamic", "fixed")
      DO ic = 1,nc
        DO is = 1,INT(num_sl_above_bedrock(ic))
          IF (microbial_uptake_nh4_sl(ic, is) >= 0.0_wp) THEN
            hlp1(ic, is) = microbial_uptake_nh4_sl(ic, is) + plant_uptake_nh4_sl(ic, is) + nitrification_no3(ic, is) + &
                        mycorrhiza_uptake_nh4_sl(ic, is)
            hlp2(ic, is) = nh4_solute(ic, is) + transport_nh4_solute(ic, is)
            IF (hlp2(ic, is) < hlp1(ic, is)) THEN
              IF (hlp1(ic, is) > eps8) THEN
                microbial_uptake_nh4_sl(ic, is)  = hlp2(ic, is) * microbial_uptake_nh4_sl(ic, is) / hlp1(ic, is)
                plant_uptake_nh4_sl(ic, is)      = hlp2(ic, is) * plant_uptake_nh4_sl(ic, is) / hlp1(ic, is)
                nitrification_no3(ic, is)        = hlp2(ic, is) * nitrification_no3(ic, is) / hlp1(ic, is)
                mycorrhiza_uptake_nh4_sl(ic, is) = hlp2(ic, is) * mycorrhiza_uptake_nh4_sl(ic, is) / hlp1(ic, is)
              ELSE
                microbial_uptake_nh4_sl(ic, is)  = 0.0_wp
                plant_uptake_nh4_sl(ic, is)      = 0.0_wp
                nitrification_no3(ic, is)        = 0.0_wp
                mycorrhiza_uptake_nh4_sl(ic, is) = 0.0_wp
              ENDIF
            ENDIF
          ELSE
            hlp1(ic, is) = plant_uptake_nh4_sl(ic, is) + nitrification_no3(ic, is) + mycorrhiza_uptake_nh4_sl(ic, is)
            hlp2(ic, is) = nh4_solute(ic, is) + transport_nh4_solute(ic, is) - microbial_uptake_nh4_sl(ic, is)
            IF (hlp2(ic, is) < hlp1(ic, is)) THEN
              IF (hlp1(ic, is) > eps8) THEN
                plant_uptake_nh4_sl(ic, is)      = hlp2(ic, is) * plant_uptake_nh4_sl(ic, is) / hlp1(ic, is)
                nitrification_no3(ic, is)        = hlp2(ic, is) * nitrification_no3(ic, is) / hlp1(ic, is)
                mycorrhiza_uptake_nh4_sl(ic, is) = hlp2(ic, is) * mycorrhiza_uptake_nh4_sl(ic, is) / hlp1(ic, is)
              ELSE
                plant_uptake_nh4_sl(ic, is)      = 0.0_wp
                nitrification_no3(ic, is)        = 0.0_wp
                mycorrhiza_uptake_nh4_sl(ic, is) = 0.0_wp
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    CASE DEFAULT
      CALL finish(TRIM(routine),'No such bnf_scheme available here.')
    END SELECT

    !> 1.4  limit plant uptake and denitrification to available NO3 minus transport losses
    !!
    DO ic = 1, nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        IF (microbial_uptake_no3_sl(ic,is) >= 0.0_wp) THEN
          hlp1(ic,is) = microbial_uptake_no3_sl(ic,is) + plant_uptake_no3_sl(ic,is) + &
                  denitrification_n2(ic,is) + mycorrhiza_uptake_no3_sl(ic,is)
          hlp2(ic,is) = no3_solute(ic,is) + transport_no3_solute(ic,is)
          IF (hlp2(ic,is) < hlp1(ic,is)) THEN
            IF (hlp1(ic,is) > eps8) THEN
              microbial_uptake_no3_sl(ic,is)  = hlp2(ic,is) * microbial_uptake_no3_sl(ic,is) / hlp1(ic,is)
              plant_uptake_no3_sl(ic,is)      = hlp2(ic,is) * plant_uptake_no3_sl(ic,is) / hlp1(ic,is)
              denitrification_n2(ic,is)       = hlp2(ic,is) * denitrification_n2(ic,is) / hlp1(ic,is)
              mycorrhiza_uptake_no3_sl(ic,is) = hlp2(ic,is) * mycorrhiza_uptake_no3_sl(ic,is) / hlp1(ic,is)
            ELSE
              microbial_uptake_no3_sl(ic,is)  = 0.0_wp
              plant_uptake_no3_sl(ic,is)      = 0.0_wp
              denitrification_n2(ic,is)       = 0.0_wp
              mycorrhiza_uptake_no3_sl(ic,is) = 0.0_wp
            ENDIF
          ENDIF
        ELSE
          hlp1(ic,is) = plant_uptake_no3_sl(ic,is) + denitrification_n2(ic,is) + mycorrhiza_uptake_no3_sl(ic,is)
          hlp2(ic,is) = no3_solute(ic,is) + transport_no3_solute(ic,is) - microbial_uptake_no3_sl(ic,is)
          IF (hlp2(ic,is) < hlp1(ic,is)) THEN
            IF (hlp1(ic,is) > eps8) THEN
              plant_uptake_no3_sl(ic,is)      = hlp2(ic,is) * plant_uptake_no3_sl(ic,is) / hlp1(ic,is)
              denitrification_n2(ic,is)       = hlp2(ic,is) * denitrification_n2(ic,is) / hlp1(ic,is)
              mycorrhiza_uptake_no3_sl(ic,is) = hlp2(ic,is) * mycorrhiza_uptake_no3_sl(ic,is) / hlp1(ic,is)
            ELSE
              plant_uptake_no3_sl(ic,is)      = 0.0_wp
              denitrification_n2(ic,is)       = 0.0_wp
              mycorrhiza_uptake_no3_sl(ic,is) = 0.0_wp
            ENDIF
          ENDIF
        ENDIF  !< microbial_uptake_no3_sl(ic,is) >= 0.0_wp
      ENDDO
    ENDDO

    !> 2.0 determine actual growth rate of fast SOM pool
    !!

    !> 2.1 Calculate actual litter loss rate, calculated such that litter decomposition
    !! matches the actual immobilisation flux and CN uptake from litter given actual CUE
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        hlp1(ic, is) = frac_litter2fast_som / apparent_fast_som_cn(ic, is) &
          &         * (sb_loss_mt(ix_soluable_litter, ixC, ic, is) + sb_loss_mt(ix_polymeric_litter, ixC, ic, is)) &
          &         - (asymb_n_fixation(ic, is) + microbial_nue &
          &         * (sb_loss_mt(ix_soluable_litter, ixN, ic, is) + sb_loss_mt(ix_polymeric_litter, ixN, ic, is)))
        IF ((microbial_uptake_nh4_sl(ic, is) + microbial_uptake_no3_sl(ic, is)) < hlp1(ic, is)) THEN
          IF (hlp1(ic, is) > eps8) THEN
            hlp2(ic, is) = (microbial_uptake_nh4_sl(ic, is) + microbial_uptake_no3_sl(ic, is)) / hlp1(ic, is)
          ELSE
            hlp2(ic, is) = 0.0_wp
          ENDIF
          sb_loss_mt(ix_soluable_litter, ixC, ic, is)    = hlp2(ic, is) * sb_loss_mt(ix_soluable_litter, ixC, ic, is)
          sb_loss_mt(ix_soluable_litter, ixN, ic, is)    = hlp2(ic, is) * sb_loss_mt(ix_soluable_litter, ixN, ic, is)
          sb_loss_mt(ix_soluable_litter, ixP, ic, is)    = hlp2(ic, is) * sb_loss_mt(ix_soluable_litter, ixP, ic, is)
          sb_loss_mt(ix_soluable_litter, ixC13, ic, is)  = hlp2(ic, is) * sb_loss_mt(ix_soluable_litter, ixC13, ic, is)
          sb_loss_mt(ix_soluable_litter, ixC14, ic, is)  = hlp2(ic, is) * sb_loss_mt(ix_soluable_litter, ixC14, ic, is)
          sb_loss_mt(ix_soluable_litter, ixN15, ic, is)  = hlp2(ic, is) * sb_loss_mt(ix_soluable_litter, ixN15, ic, is)

          sb_loss_mt(ix_polymeric_litter, ixC, ic, is)   = hlp2(ic, is) * sb_loss_mt(ix_polymeric_litter, ixC, ic, is)
          sb_loss_mt(ix_polymeric_litter, ixN, ic, is)   = hlp2(ic, is) * sb_loss_mt(ix_polymeric_litter, ixN, ic, is)
          sb_loss_mt(ix_polymeric_litter, ixP, ic, is)   = hlp2(ic, is) * sb_loss_mt(ix_polymeric_litter, ixP, ic, is)
          sb_loss_mt(ix_polymeric_litter, ixC13, ic, is) = hlp2(ic, is) * sb_loss_mt(ix_polymeric_litter, ixC13, ic, is)
          sb_loss_mt(ix_polymeric_litter, ixC14, ic, is) = hlp2(ic, is) * sb_loss_mt(ix_polymeric_litter, ixC14, ic, is)
          sb_loss_mt(ix_polymeric_litter, ixN15, ic, is) = hlp2(ic, is) * sb_loss_mt(ix_polymeric_litter, ixN15, ic, is)
        ENDIF
      ENDDO
    ENDDO

    !>3.0 calculate potential microbial P uptake from soluable pool (mol /m3 / timestep)
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        microbial_uptake_po4_sl(ic, is) = frac_litter2fast_som / apparent_fast_som_cn(ic, is) / k_fast_som_np &
          &                               * (sb_loss_mt(ix_soluable_litter, ixC, ic, is) &
          &                               + sb_loss_mt(ix_polymeric_litter, ixC, ic, is)) &
          &                               - microbial_pue &
          &                               * (sb_loss_mt(ix_soluable_litter, ixP, ic, is) &
          &                               + sb_loss_mt(ix_polymeric_litter, ixP, ic, is))

        !>  3.1 constrain P uptake by microbes and plants to the available sum of P from P_soluable,
        !>      biomineralisation and transport (mol /m3 / timestep)
        !>
        IF (microbial_uptake_po4_sl(ic, is) > 0.0_wp) THEN ! if microbes take up P, constrain both microbial and plant uptake to available P
          hlp1(ic, is) = microbial_uptake_po4_sl(ic, is) + plant_uptake_po4_sl(ic, is) + mycorrhiza_uptake_po4_sl(ic, is)
          hlp2(ic, is) = po4_solute(ic, is) + transport_po4_solute(ic, is)
          IF (hlp2(ic, is) < hlp1(ic, is)) THEN
            IF (hlp1(ic, is) > eps8) THEN
              microbial_uptake_po4_sl(ic, is)  = hlp2(ic, is) * microbial_uptake_po4_sl(ic, is) / hlp1(ic, is)
              plant_uptake_po4_sl(ic, is)      = hlp2(ic, is) * plant_uptake_po4_sl(ic, is) / hlp1(ic, is)
              mycorrhiza_uptake_po4_sl(ic, is) = hlp2(ic, is) * mycorrhiza_uptake_po4_sl(ic, is) / hlp1(ic, is)
            ELSE
              microbial_uptake_po4_sl(ic, is)  = 0.0_wp
              plant_uptake_po4_sl(ic, is)      = 0.0_wp
              mycorrhiza_uptake_po4_sl(ic, is) = 0.0_wp
            ENDIF
          ENDIF
        ELSE ! if microbes release P, constrain plant uptake to available P, note the negative sign of microbial_uptake_po4_sl(ic, is)
          hlp1(ic, is) = plant_uptake_po4_sl(ic, is) + mycorrhiza_uptake_po4_sl(ic, is)
          hlp2(ic, is) = po4_solute(ic, is) + transport_po4_solute(ic, is) - microbial_uptake_po4_sl(ic, is)
          IF (hlp2(ic, is) < hlp1(ic, is)) THEN
            IF (hlp1(ic, is) > eps8) THEN
              plant_uptake_po4_sl(ic, is)      = hlp2(ic, is) * plant_uptake_po4_sl(ic, is) / hlp1(ic, is)
              mycorrhiza_uptake_po4_sl(ic, is) = hlp2(ic, is) * mycorrhiza_uptake_po4_sl(ic, is) / hlp1(ic, is)
            ELSE
              plant_uptake_po4_sl(ic, is)      = 0.0_wp
              mycorrhiza_uptake_po4_sl(ic, is) = 0.0_wp
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE calc_actual_flux_rates

  ! ======================================================================================================= !
  !>For the sb_nloss_scheme = fixed: calculates the amount of N lost in proportion to
  !>  net mineralisation \n
  !>
  !>For the sb_nloss_scheme = dynamic: calculates the partitioning of nitrification and denitrification
  !>  to NO3, NOy, N2O and N2, and diagnoses the amount of leaching below the rooting zone
  !>
  SUBROUTINE calc_partitioning_nitrification_denitrification_rate( &
                                    nc, &                         ! in
                                    nsoil_sb, &
                                    num_sl_above_bedrock, &
                                    rtm_denitrification, &
                                    nitrification_no3, &
                                    denitrification_n2, &
                                    nitrification_noy, &
                                    nitrification_n2o, &
                                    denitrification_noy, &
                                    denitrification_n2o)

    USE mo_sb_constants,            ONLY: f_nit_noy, f_nit_n2o, f_denit_noy, f_denit_n2o
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                      !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                !< number of soil layers
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock    !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: rtm_denitrification     !< temperature modifier for denitrification
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: nitrification_no3, &    !< nitrification to NO3 [mol / m3 / timestep]
                                                        denitrification_n2      !< denitrification to N2 [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: nitrification_noy, &    !< nitrification to NOy [mol / m3 / timestep]
                                                        nitrification_n2o, &    !< nitrification to N2O [mol / m3 / timestep]
                                                        denitrification_noy, &  !< denitrification to NOy [mol / m3 / timestep]
                                                        denitrification_n2o     !< denitrification to N2O [mol / m3 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ic, is
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_partitioning_nitrification_denitrification_rate'
    ! ----------------------------------------------------------------------------------------------------- !

    !> init out var
    !>
    nitrification_noy(:,:)    = 0.0_wp
    nitrification_n2o(:,:)    = 0.0_wp
    denitrification_noy(:,:)  = 0.0_wp
    denitrification_n2o(:,:)  = 0.0_wp

    !>1.0 case of dynamic N loss calculations, following Zaehle et al. 2011
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        ! nitrification rate to produce NO3, NOy and N2O from NH4 [mol m-3 timestep-1]
        nitrification_noy(ic, is)   = f_nit_noy * nitrification_no3(ic, is)
        nitrification_n2o(ic, is)   = f_nit_n2o * nitrification_no3(ic, is)
        nitrification_no3(ic, is)   = nitrification_no3(ic, is) - nitrification_noy(ic, is) - nitrification_n2o(ic, is)
        ! denitrification rate to produce NOy, N2O and N2 from NO3 [mol m-3 timestep-1]
        denitrification_noy(ic, is) = f_denit_noy * rtm_denitrification(ic, is) * denitrification_n2(ic, is)
        denitrification_n2o(ic, is) = f_denit_n2o * rtm_denitrification(ic, is) * denitrification_n2(ic, is)
        denitrification_n2(ic, is)  = denitrification_n2(ic, is) - denitrification_noy(ic, is) - denitrification_n2o(ic, is)
      ENDDO
    ENDDO
  END SUBROUTINE calc_partitioning_nitrification_denitrification_rate


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to both: update_sb_simple_model and update_sb_jsm
  !   units for flux rates differ between the two models: \n
  !     in this routine they are assumed to be in [mol/m3/timestep] \n
  !     simple soil model: [mol/m3/timestep]
  !     JSM:               [micro-mol/m3/s] (unit conversion is applied within the routine)
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the N15 fluxes based on uptake or transformation of the inorganic NH4 / NO3 pools
  !! according to the framework proposed by Robinson, 2001, TREE. Includes a representation of the
  !! Rayleigh model to account for the change in fractionation with increasing consumption
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_sb_inorganic_n15_fluxes(nc, &
                                          nsoil_sb, &
                                          dtime, &
                                          num_sl_above_bedrock, &
                                          sb_model_scheme, &
                                          nh4_solute, &
                                          nh4_n15_solute, &
                                          no3_solute, &
                                          no3_n15_solute, &
                                          plant_uptake_nh4_sl, &
                                          plant_uptake_no3_sl, &
                                          mycorrhiza_uptake_nh4_sl, &
                                          mycorrhiza_uptake_no3_sl, &
                                          microbial_uptake_nh4_sl, &
                                          microbial_uptake_no3_sl, &
                                          asymb_n_fixation, &
                                          nitrification_no3, &
                                          nitrification_noy, &
                                          nitrification_n2o, &
                                          denitrification_noy, &
                                          denitrification_n2o, &
                                          denitrification_n2, &
                                          transport_nh4_n15_solute, &
                                          transport_no3_n15_solute, &
                                          plant_uptake_nh4_n15_sl, &
                                          plant_uptake_no3_n15_sl, &
                                          mycorrhiza_uptake_nh4_n15_sl, &
                                          mycorrhiza_uptake_no3_n15_sl, &
                                          microbial_uptake_nh4_n15_sl, &
                                          microbial_uptake_no3_n15_sl, &
                                          asymb_n_fixation_n15, &
                                          nitrification_no3_n15, &
                                          nitrification_noy_n15, &
                                          nitrification_n2o_n15, &
                                          denitrification_noy_n15, &
                                          denitrification_n2o_n15, &
                                          denitrification_n2_n15)

    USE mo_isotope_util,   ONLY: calc_fractionation, calc_mixing_ratio_N15N14
    USE mo_sb_constants,   ONLY: eta_plant_uptake_nh4, eta_plant_uptake_no3, eta_mic_uptake_nh4, eta_mic_uptake_no3, &
      &                          eta_nitrate_production, eta_nitrification, eta_denitrification
    USE mo_veg_constants,  ONLY: eta_nfixation
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                            INTENT(in)    :: nc, &                       !< dimensions
                                                         nsoil_sb
    REAL(wp),                           INTENT(in)    :: dtime                       !< timestep length
    REAL(wp), DIMENSION(nc),            INTENT(in)    :: num_sl_above_bedrock        !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    CHARACTER(len=*),                   INTENT(in)    :: sb_model_scheme             !< simple_1d / jsm from sb config
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: nh4_solute, &               !< soil solute NH4 concentration [mol m-3]
                                                         nh4_n15_solute, &           !< soil solute 15NH4 concentration [mol m-3]
                                                         no3_solute, &               !< soil solute NO3 concentration [mol m-3]
                                                         no3_n15_solute, &           !< soil solute 15NO3 concentration [mol m-3]
                                                         plant_uptake_nh4_sl, &      !< NH4 uptake by plants [mol / m3 / timestep]
                                                         plant_uptake_no3_sl, &      !< NO3 uptake by plants [mol / m3 / timestep]
                                                         mycorrhiza_uptake_nh4_sl, & !< NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                         mycorrhiza_uptake_no3_sl, & !< NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                         microbial_uptake_nh4_sl, &  !< NH4 uptake by microbes [mol / m3 / timestep]
                                                         microbial_uptake_no3_sl, &  !< NO3 uptake by microbes [mol / m3 / timestep]
                                                         asymb_n_fixation, &         !< asymbiotic N fixation [mol / m3 / timestep]
                                                         nitrification_no3, &        !< nitrification to NO3 [mol / m3 / timestep]
                                                         nitrification_noy, &        !< nitrification to NOy [mol / m3 / timestep]
                                                         nitrification_n2o, &        !< nitrification to N2O [mol / m3 / timestep]
                                                         denitrification_noy, &      !< denitrification to NOy [mol / m3 / timestep]
                                                         denitrification_n2o, &      !< denitrification to N2O [mol / m3 / timestep]
                                                         denitrification_n2, &       !< denitrification to N2 [mol / m3 / timestep]
                                                         transport_nh4_n15_solute, & !< transport of 15NH4 [mol / m3 / timestep]
                                                         transport_no3_n15_solute    !< transport of 15NO3 [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(out)   :: plant_uptake_nh4_n15_sl, &      !< 15NH4 uptake by plants [mol / m3 / timestep]
                                                         plant_uptake_no3_n15_sl, &      !< 15NO3 uptake by plants [mol / m3 / timestep]
                                                         mycorrhiza_uptake_nh4_n15_sl, & !< 15NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                         mycorrhiza_uptake_no3_n15_sl, & !< 15NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                         microbial_uptake_nh4_n15_sl, &  !< 15NH4 uptake by microbes [mol / m3 / timestep]
                                                         microbial_uptake_no3_n15_sl, &  !< 15NH4 uptake by microbes [mol / m3 / timestep]
                                                         asymb_n_fixation_n15, &         !< asymbiotic 15N fixation [mol / m3 / timestep]
                                                         nitrification_no3_n15, &        !< nitrification to 15NO3 [mol / m3 / timestep]
                                                         nitrification_noy_n15, &        !< nitrification to 15NOy [mol / m3 / timestep]
                                                         nitrification_n2o_n15, &        !< nitrification to 15N2O [mol / m3 / timestep]
                                                         denitrification_noy_n15, &      !< denitrification to 15NOy [mol / m3 / timestep]
                                                         denitrification_n2o_n15, &      !< denitrification to 15N2O [mol / m3 / timestep]
                                                         denitrification_n2_n15          !< denitrification to 15N2 [mol / m3 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: ic, is
    REAL(wp), DIMENSION(nc, nsoil_sb)         :: hlp1, hlp2
    REAL(wp), DIMENSION(nc, nsoil_sb)         :: &                ! flux rates with unit conversion for JSM
                                                 plant_uptake_nh4_local_sl, &
                                                 plant_uptake_no3_local_sl, &
                                                 mycorrhiza_uptake_nh4_local_sl, &
                                                 mycorrhiza_uptake_no3_local_sl, &
                                                 microbial_uptake_nh4_local_sl, &
                                                 microbial_uptake_no3_local_sl, &
                                                 asymb_n_fixation_local, &
                                                 nitrification_no3_local, &
                                                 nitrification_noy_local, &
                                                 nitrification_n2o_local, &
                                                 denitrification_noy_local, &
                                                 denitrification_n2o_local, &
                                                 denitrification_n2_local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sb_inorganic_n15_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !


    !> init out var
    !>
    plant_uptake_nh4_n15_sl(:,:)      = 0.0_wp
    plant_uptake_no3_n15_sl(:,:)      = 0.0_wp
    mycorrhiza_uptake_nh4_n15_sl(:,:) = 0.0_wp
    mycorrhiza_uptake_no3_n15_sl(:,:) = 0.0_wp
    microbial_uptake_nh4_n15_sl(:,:)  = 0.0_wp
    microbial_uptake_no3_n15_sl(:,:)  = 0.0_wp
    asymb_n_fixation_n15(:,:)         = 0.0_wp
    nitrification_no3_n15(:,:)        = 0.0_wp
    nitrification_noy_n15(:,:)        = 0.0_wp
    nitrification_n2o_n15(:,:)        = 0.0_wp
    denitrification_noy_n15(:,:)      = 0.0_wp
    denitrification_n2o_n15(:,:)      = 0.0_wp
    denitrification_n2_n15(:,:)       = 0.0_wp

    !>0.9 'INTENT(in) valiables': flux rates unit conversion when this routines is called from JSM
    !>  unit conversion: [micro-mol/m3/s] -> [mol/m3/timestep]
    !>
    !> @NOTE @TODO ? not all INTENT(in) variables are considered ?
    IF(sb_model_scheme == 'jsm') THEN
      plant_uptake_nh4_local_sl(:,:)      = plant_uptake_nh4_sl(:,:)      * dtime / 1.e6_wp
      plant_uptake_no3_local_sl(:,:)      = plant_uptake_no3_sl(:,:)      * dtime / 1.e6_wp
      mycorrhiza_uptake_nh4_local_sl(:,:) = mycorrhiza_uptake_nh4_sl(:,:) * dtime / 1.e6_wp
      mycorrhiza_uptake_no3_local_sl(:,:) = mycorrhiza_uptake_no3_sl(:,:) * dtime / 1.e6_wp
      microbial_uptake_nh4_local_sl(:,:)  = microbial_uptake_nh4_sl(:,:)  * dtime / 1.e6_wp
      microbial_uptake_no3_local_sl(:,:)  = microbial_uptake_no3_sl(:,:)  * dtime / 1.e6_wp
      asymb_n_fixation_local(:,:)         = asymb_n_fixation(:,:)         * dtime / 1.e6_wp
      nitrification_no3_local(:,:)        = nitrification_no3(:,:)        * dtime / 1.e6_wp
      nitrification_noy_local(:,:)        = nitrification_noy(:,:)        * dtime / 1.e6_wp
      nitrification_n2o_local(:,:)        = nitrification_n2o(:,:)        * dtime / 1.e6_wp
      denitrification_noy_local(:,:)      = denitrification_noy(:,:)      * dtime / 1.e6_wp
      denitrification_n2o_local(:,:)      = denitrification_n2o(:,:)      * dtime / 1.e6_wp
      denitrification_n2_local(:,:)       = denitrification_n2(:,:)       * dtime / 1.e6_wp
    ELSE
      plant_uptake_nh4_local_sl(:,:)      = plant_uptake_nh4_sl(:,:)
      plant_uptake_no3_local_sl(:,:)      = plant_uptake_no3_sl(:,:)
      mycorrhiza_uptake_nh4_local_sl(:,:) = mycorrhiza_uptake_nh4_sl(:,:)
      mycorrhiza_uptake_no3_local_sl(:,:) = mycorrhiza_uptake_no3_sl(:,:)
      microbial_uptake_nh4_local_sl(:,:)  = microbial_uptake_nh4_sl(:,:)
      microbial_uptake_no3_local_sl(:,:)  = microbial_uptake_no3_sl(:,:)
      asymb_n_fixation_local(:,:)         = asymb_n_fixation(:,:)
      nitrification_no3_local(:,:)        = nitrification_no3(:,:)
      nitrification_noy_local(:,:)        = nitrification_noy(:,:)
      nitrification_n2o_local(:,:)        = nitrification_n2o(:,:)
      denitrification_noy_local(:,:)      = denitrification_noy(:,:)
      denitrification_n2o_local(:,:)      = denitrification_n2o(:,:)
      denitrification_n2_local(:,:)       = denitrification_n2(:,:)
    ENDIF

    !> 1.0 calculate the actual discrimination for NH4,
    !!     which depends on the fraction of the source consumed f = (sum of fluxes)/stock = hlp1.
    !!     The discrimination is then reduced by a reduction factor r = - (1-f) ln (1-f) / f = hlp2
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        ! fraction of NH4 that is consumed
        IF (nh4_solute(ic, is) > eps8) THEN
          hlp1(ic, is) = MAX((plant_uptake_nh4_local_sl(ic, is) + mycorrhiza_uptake_nh4_local_sl(ic, is) + &
            &              microbial_uptake_nh4_local_sl(ic, is) + nitrification_no3_local(ic, is) + &
            &              nitrification_noy_local(ic, is) + nitrification_n2o_local(ic, is)) / nh4_solute(ic, is), eps8)
        ELSE
          hlp1(ic, is) = 1.0_wp
        ENDIF
        ! calculate reduction factor for discrimination
        IF (hlp1(ic, is) < (1.0_wp - eps8)) THEN
          hlp2(ic, is) = (hlp1(ic, is) - 1.0_wp) * LOG(1.0_wp - hlp1(ic, is)) / hlp1(ic, is)
        ELSE
          hlp2(ic, is) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO

    !> 1.1 15N fluxes involving ammonium
    !!
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        plant_uptake_nh4_n15_sl(ic, is)      = plant_uptake_nh4_local_sl(ic, is) * &
          &                                    calc_fractionation(nh4_solute(ic, is), nh4_n15_solute(ic, is), &
          &                                    eta_plant_uptake_nh4 * hlp2(ic, is))
        mycorrhiza_uptake_nh4_n15_sl(ic, is) = mycorrhiza_uptake_nh4_local_sl(ic, is) * &
          &                                    calc_fractionation(nh4_solute(ic, is), nh4_n15_solute(ic, is), &
          &                                    eta_plant_uptake_nh4 * hlp2(ic, is))
        microbial_uptake_nh4_n15_sl(ic, is)  = microbial_uptake_nh4_local_sl(ic, is) * &
          &                                    calc_fractionation(nh4_solute(ic, is), nh4_n15_solute(ic, is), &
          &                                    eta_mic_uptake_nh4 * hlp2(ic, is))
        nitrification_no3_n15(ic, is)        = nitrification_no3_local(ic, is) * &
          &                                    calc_fractionation(nh4_solute(ic, is), nh4_n15_solute(ic, is), &
          &                                    eta_nitrate_production * hlp2(ic, is))
        nitrification_noy_n15(ic, is)        = nitrification_noy_local(ic, is) * &
          &                                    calc_fractionation(nh4_solute(ic, is), nh4_n15_solute(ic, is), &
          &                                    eta_nitrification * hlp2(ic, is))
        nitrification_n2o_n15(ic, is)        = nitrification_n2o_local(ic, is) * &
          &                                    calc_fractionation(nh4_solute(ic, is), nh4_n15_solute(ic, is), &
          &                                    eta_nitrification * hlp2(ic, is))
      ENDDO
    ENDDO

    !> 1.2 contrain 15NH4 fluxes to available 15NH4
    !!
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        hlp1(ic, is) = plant_uptake_nh4_n15_sl(ic, is) + mycorrhiza_uptake_nh4_n15_sl(ic, is) &
          &            + microbial_uptake_nh4_n15_sl(ic, is) + nitrification_no3_n15(ic, is) &
          &            + nitrification_noy_n15(ic, is) + nitrification_n2o_n15(ic, is)
        hlp2(ic, is) = nh4_n15_solute(ic, is) + transport_nh4_n15_solute(ic, is)
        IF (hlp2(ic, is) < hlp1(ic, is)) THEN
           IF (hlp1(ic, is) > eps8) THEN
              plant_uptake_nh4_n15_sl(ic, is)      = hlp2(ic, is) * plant_uptake_nh4_n15_sl(ic, is) / hlp1(ic, is)
              mycorrhiza_uptake_nh4_n15_sl(ic, is) = hlp2(ic, is) * mycorrhiza_uptake_nh4_n15_sl(ic, is) / hlp1(ic, is)
              microbial_uptake_nh4_n15_sl(ic, is)  = hlp2(ic, is) * microbial_uptake_nh4_n15_sl(ic, is) / hlp1(ic, is)
              nitrification_no3_n15(ic, is)        = hlp2(ic, is) * nitrification_no3_n15(ic, is) / hlp1(ic, is)
              nitrification_noy_n15(ic, is)        = hlp2(ic, is) * nitrification_noy_n15(ic, is) / hlp1(ic, is)
              nitrification_n2o_n15(ic, is)        = hlp2(ic, is) * nitrification_n2o_n15(ic, is) / hlp1(ic, is)
           ELSE
              plant_uptake_nh4_n15_sl(ic, is)      = 0.0_wp
              mycorrhiza_uptake_nh4_n15_sl(ic, is) = 0.0_wp
              microbial_uptake_nh4_n15_sl(ic, is)  = 0.0_wp
              nitrification_no3_n15(ic, is)        = 0.0_wp
              nitrification_noy_n15(ic, is)        = 0.0_wp
              nitrification_n2o_n15(ic, is)        = 0.0_wp
           ENDIF
        ENDIF
      ENDDO
    ENDDO

    !> 2.0 calculate the actual discrimination for NO3,
    !!     which depends on the fraction of the source consumed f = (sum of fluxes)/stock = hlp1.
    !!     The discrimination is then reduced by a reduction factor r = - (1-f) ln (1-f) / f = hlp2
    ! fraction of NO3 that is consumed
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        IF (no3_solute(ic, is) > eps8) THEN
           hlp1(ic, is) = MAX((plant_uptake_no3_local_sl(ic, is) + mycorrhiza_uptake_no3_local_sl(ic, is) + &
                            microbial_uptake_no3_local_sl(ic, is) + denitrification_noy_local(ic, is) + &
                            denitrification_n2o_local(ic, is) + denitrification_n2_local(ic, is)) / no3_solute(ic, is), eps8)
        ELSE
           hlp1(ic, is) = 1.0_wp
        ENDIF
        IF (hlp1(ic, is) < (1.0_wp - eps8)) THEN
           hlp2(ic, is) = (hlp1(ic, is) - 1.0_wp) * LOG(1.0_wp - hlp1(ic, is)) / hlp1(ic, is)
        ELSE
           hlp2(ic, is) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO

    !> 2.1 N15 fluxes involving NO3
    !!
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        plant_uptake_no3_n15_sl(ic, is)      = plant_uptake_no3_local_sl(ic, is) * &
          &                                 calc_fractionation(no3_solute(ic, is), no3_n15_solute(ic, is), &
          &                                                    eta_plant_uptake_no3 * hlp2(ic, is))
        mycorrhiza_uptake_no3_n15_sl(ic, is) = mycorrhiza_uptake_no3_local_sl(ic, is) * &
          &                                 calc_fractionation(no3_solute(ic, is), no3_n15_solute(ic, is), &
          &                                                    eta_plant_uptake_no3 * hlp2(ic, is))
        microbial_uptake_no3_n15_sl(ic, is)  = microbial_uptake_no3_local_sl(ic, is) * &
          &                                 calc_fractionation(no3_solute(ic, is), no3_n15_solute(ic, is), &
          &                                                    eta_mic_uptake_no3 * hlp2(ic, is))
        denitrification_noy_n15(ic, is)      = denitrification_noy_local(ic, is) * &
          &                                 calc_fractionation(no3_solute(ic, is), no3_n15_solute(ic, is), &
          &                                                    eta_denitrification * hlp2(ic, is))
        denitrification_n2o_n15(ic, is)      = denitrification_n2o_local(ic, is) * &
          &                                 calc_fractionation(no3_solute(ic, is), no3_n15_solute(ic, is), &
          &                                                    eta_denitrification * hlp2(ic, is))
        denitrification_n2_n15(ic, is)       = denitrification_n2_local(ic, is) * &
          &                                 calc_fractionation(no3_solute(ic, is), no3_n15_solute(ic, is), &
          &                                                    eta_denitrification * hlp2(ic, is))
      ENDDO
    ENDDO

    !> 2.2 constrain 15NO3 fluxes to available 15NO3
    !!
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        hlp1(ic, is) = plant_uptake_no3_n15_sl(ic, is) + mycorrhiza_uptake_no3_n15_sl(ic, is) &
          &            + microbial_uptake_no3_n15_sl(ic, is) + denitrification_noy_n15(ic, is) &
          &            + denitrification_n2o_n15(ic, is) + denitrification_n2_n15(ic, is)
        hlp2(ic, is) = no3_n15_solute(ic, is) + transport_no3_n15_solute(ic, is)
        IF (hlp2(ic, is) < hlp1(ic, is)) THEN
           IF (hlp1(ic, is) > eps8) THEN
              plant_uptake_no3_n15_sl(ic, is)      = hlp2(ic, is) * plant_uptake_no3_n15_sl(ic, is) / hlp1(ic, is)
              mycorrhiza_uptake_no3_n15_sl(ic, is) = hlp2(ic, is) * mycorrhiza_uptake_no3_n15_sl(ic, is) / hlp1(ic, is)
              microbial_uptake_no3_n15_sl(ic, is)  = hlp2(ic, is) * microbial_uptake_no3_n15_sl(ic, is) / hlp1(ic, is)
              denitrification_noy_n15(ic, is)      = hlp2(ic, is) * denitrification_noy_n15(ic, is) / hlp1(ic, is)
              denitrification_n2o_n15(ic, is)      = hlp2(ic, is) * denitrification_n2o_n15(ic, is) / hlp1(ic, is)
              denitrification_n2_n15(ic, is)       = hlp2(ic, is) * denitrification_n2_n15(ic, is) / hlp1(ic, is)
           ELSE
              plant_uptake_no3_n15_sl(ic, is)      = 0.0_wp
              mycorrhiza_uptake_no3_n15_sl(ic, is) = 0.0_wp
              microbial_uptake_no3_n15_sl(ic, is)  = 0.0_wp
              denitrification_noy_n15(ic, is)      = 0.0_wp
              denitrification_n2o_n15(ic, is)      = 0.0_wp
              denitrification_n2_n15(ic, is)       = 0.0_wp
           ENDIF
        ENDIF
      ENDDO
    ENDDO

    !> 3.0 15N signal of asymbiotic fixation
    !!
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        asymb_n_fixation_n15(ic, is) = asymb_n_fixation_local(ic, is) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(-eta_nfixation))
      ENDDO
    ENDDO

    !> 4.0 unit conversion: [mol/m3/timestep] -> [micro-mol/m3/s]
    !!
    !! 'INTENT(out)' variables unit conversion when this routines is called from JSM
    !!
    IF(sb_model_scheme == 'jsm') THEN
      plant_uptake_nh4_n15_sl(:,:)      = plant_uptake_nh4_n15_sl(:,:)      / dtime * 1.e6_wp
      plant_uptake_no3_n15_sl(:,:)      = plant_uptake_no3_n15_sl(:,:)      / dtime * 1.e6_wp
      mycorrhiza_uptake_nh4_n15_sl(:,:) = mycorrhiza_uptake_nh4_n15_sl(:,:) / dtime * 1.e6_wp
      mycorrhiza_uptake_no3_n15_sl(:,:) = mycorrhiza_uptake_no3_n15_sl(:,:) / dtime * 1.e6_wp
      microbial_uptake_nh4_n15_sl(:,:)  = microbial_uptake_nh4_n15_sl(:,:)  / dtime * 1.e6_wp
      microbial_uptake_no3_n15_sl(:,:)  = microbial_uptake_no3_n15_sl(:,:)  / dtime * 1.e6_wp
      asymb_n_fixation_n15(:,:)         = asymb_n_fixation_n15(:,:)         / dtime * 1.e6_wp
      nitrification_no3_n15(:,:)        = nitrification_no3_n15(:,:)        / dtime * 1.e6_wp
      nitrification_noy_n15(:,:)        = nitrification_noy_n15(:,:)        / dtime * 1.e6_wp
      nitrification_n2o_n15(:,:)        = nitrification_n2o_n15(:,:)        / dtime * 1.e6_wp
      denitrification_noy_n15(:,:)      = denitrification_noy_n15(:,:)      / dtime * 1.e6_wp
      denitrification_n2o_n15(:,:)      = denitrification_n2o_n15(:,:)      / dtime * 1.e6_wp
      denitrification_n2_n15(:,:)       = denitrification_n2_n15(:,:)       / dtime * 1.e6_wp
    ENDIF

  END SUBROUTINE calc_sb_inorganic_n15_fluxes

  ! ======================================================================================================= !
  !>Calculates the formation of SOM resulting from the decomposition: dPrec/dt = CUE_Pdon * dPdon/dt
  !>
  !>  gross mineralisation follows the stoichiometric difference of the donor and recipient pool, taking account for CUE
  !>  for fast -> slow -> fast only causes gross mineralisation because of the C:N:P and CUE
  !>  for litter -> fast gross mineralisation is typcially less than microbial uptake from
  !>                mineral pools (or organic matter for P)
  !>
  SUBROUTINE calc_SOM_formation( &
    & nc, &
    & nsoil_sb, &
    & num_sl_above_bedrock, &
    & apparent_fast_som_cn, &
    & apparent_fast_som_np, &
    & microbial_uptake_nh4_sl, &
    & microbial_uptake_nh4_n15_sl, &
    & microbial_uptake_no3_sl, &
    & microbial_uptake_no3_n15_sl, &
    & asymb_n_fixation, &
    & asymb_n_fixation_n15, &
    & microbial_uptake_po4_sl, &
    & het_respiration, &
    & het_respiration_c13, &
    & het_respiration_c14, &
    & net_mineralisation_nh4, &
    & net_mineralisation_no3, &
    & net_mineralisation_po4, &
    & net_mineralisation_nh4_n15, &
    & net_mineralisation_no3_n15, &
    & biochem_mineralisation_po4, &
    & sb_loss_mt, &
    & sb_formation_mt )

    USE mo_isotope_util,        ONLY: calc_fractionation
    USE mo_sb_constants,        ONLY: frac_litter2fast_som, frac_fast2slow_som, frac_slow2fast_som, &
      &                               k_slow_som_cn, k_slow_som_np, microbial_nue, microbial_pue, &
      &                               eta_ammonification, eta_nitrification
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                              !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                        !< number of soil layers
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock            !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: apparent_fast_som_cn, &         !< apparent carbon-nitrogen ratio of microbial uptake [mol/mol]
                                                        apparent_fast_som_np, &         !< apparent carbon-nitrogen ratio of microbial uptake [mol/mol]
                                                        microbial_uptake_nh4_sl, &      !< microbial NH4 uptake [mol / m3 / timestep]
                                                        microbial_uptake_nh4_n15_sl, &  !< microbial 15NH4 uptake [mol / m3 / timestep]
                                                        microbial_uptake_no3_sl, &      !< microbial NO3 uptake [mol / m3 / timestep]
                                                        microbial_uptake_no3_n15_sl, &  !< microbial 15NO3 uptake [mol / m3 / timestep]
                                                        asymb_n_fixation, &             !< asymbiotic N fixation [mol / m3 / timestep]
                                                        asymb_n_fixation_n15, &         !< asymbiotic 15N fixation [mol / m3 / timestep]
                                                        microbial_uptake_po4_sl         !< microbial PO4 uptake [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: het_respiration, &              !< heterotrophic respiration [mol C / m3 / timestep]
                                                        het_respiration_c13, &          !< heterotrophic respiration [mol 13C / m3 / timestep]
                                                        het_respiration_c14, &          !< heterotrophic respiration [micro-mol 14C / m3 / timestep]
                                                        net_mineralisation_nh4, &       !< net mineralisation of NH4 [mol / m3 / timestep]
                                                        net_mineralisation_no3, &       !< net mineralisation of NO3 [mol / m3 / timestep]
                                                        net_mineralisation_po4, &       !< net mineralisation of PO4 [mol / m3 / timestep]
                                                        net_mineralisation_nh4_n15, &   !< net mineralisation of 15NH4 [mol / m3 / timestep]
                                                        net_mineralisation_no3_n15, &   !< net mineralisation of 15NO3 [mol / m3 / timestep]
                                                        biochem_mineralisation_po4      !< biomineralisation rate [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_loss_mt(:,:,:,:)             !< bgcm sb_loss flux: turnover of SOM and litter pools  [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)        !< bgcm sb_formation flux: formation of SOM and litter pools  [mol / m3 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                          :: ic, is       !< loop over dimensions
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp2
    CHARACTER(len=*), PARAMETER                      :: routine = TRIM(modname)//':calc_SOM_formation'
    ! ----------------------------------------------------------------------------------------------------- !

    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        !>1.0 formation of fast SOM from litter
        !>
        ! calculate 15N mineralised during litter decay, remainder enters fast pool
        hlp1(ic, is) = (1.0_wp - microbial_nue) * sb_loss_mt(ix_soluable_litter, ixN, ic, is) &
          &         * calc_fractionation(sb_loss_mt(ix_soluable_litter, ixN, ic, is), &
          &                              sb_loss_mt(ix_soluable_litter, ixN15, ic, is), &
          &                              eta_ammonification)
        hlp2(ic, is) = (1.0_wp - microbial_nue) * sb_loss_mt(ix_polymeric_litter, ixN, ic, is) &
          &         * calc_fractionation(sb_loss_mt(ix_polymeric_litter, ixN, ic, is), &
          &                              sb_loss_mt(ix_polymeric_litter, ixN15, ic, is), &
          &                              eta_ammonification)
        ! all elements
        sb_formation_mt(ix_microbial, ixC, ic, is)   = sb_formation_mt(ix_microbial, ixC, ic, is) + frac_litter2fast_som &
          &                                          * (sb_loss_mt(ix_soluable_litter, ixC, ic, is) &
          &                                          + sb_loss_mt(ix_polymeric_litter, ixC, ic, is))
        sb_formation_mt(ix_microbial, ixN, ic, is)   = sb_formation_mt(ix_microbial, ixN, ic, is) + microbial_nue &
          &                                          * (sb_loss_mt(ix_soluable_litter, ixN, ic, is) &
          &                                          + sb_loss_mt(ix_polymeric_litter, ixN, ic, is)) &
          &                                          + microbial_uptake_nh4_sl(ic, is) + microbial_uptake_no3_sl(ic, is) &
          &                                          + asymb_n_fixation(ic, is)
        sb_formation_mt(ix_microbial, ixP, ic, is)   = sb_formation_mt(ix_microbial, ixP, ic, is) + microbial_pue &
          &                                          * (sb_loss_mt(ix_soluable_litter, ixP, ic, is) &
          &                                          + sb_loss_mt(ix_polymeric_litter, ixP, ic, is)) &
          &                                          + microbial_uptake_po4_sl(ic, is)
        sb_formation_mt(ix_microbial, ixC13, ic, is) = sb_formation_mt(ix_microbial, ixC13, ic, is) + frac_litter2fast_som &
          &                                          * (sb_loss_mt(ix_soluable_litter, ixC13, ic, is) &
          &                                          + sb_loss_mt(ix_polymeric_litter, ixC13, ic, is))
        sb_formation_mt(ix_microbial, ixC14, ic, is) = sb_formation_mt(ix_microbial, ixC14, ic, is) + frac_litter2fast_som &
          &                                          * (sb_loss_mt(ix_soluable_litter, ixC14, ic, is) &
          &                                          + sb_loss_mt(ix_polymeric_litter, ixC14, ic, is))
        sb_formation_mt(ix_microbial, ixN15, ic, is) = sb_formation_mt(ix_microbial, ixN15, ic, is) &
          &                                          + sb_loss_mt(ix_soluable_litter, ixN15, ic, is)  - hlp1(ic, is) &
          &                                          + sb_loss_mt(ix_polymeric_litter, ixN15, ic, is) - hlp2(ic, is) &
          &                                          + microbial_uptake_nh4_n15_sl(ic, is) + microbial_uptake_no3_n15_sl(ic, is) &
          &                                          + asymb_n_fixation_n15(ic, is)

        het_respiration(ic, is)            = het_respiration(ic, is) + (1._wp - frac_litter2fast_som) * &
          &                                                      (sb_loss_mt(ix_soluable_litter, ixC, ic, is) + &
          &                                                       sb_loss_mt(ix_polymeric_litter, ixC, ic, is))
        het_respiration_c13(ic, is)        = het_respiration_c13(ic, is) + (1._wp - frac_litter2fast_som) * &
          &                                                          (sb_loss_mt(ix_soluable_litter, ixC13, ic, is) + &
          &                                                           sb_loss_mt(ix_polymeric_litter, ixC13, ic, is))
        het_respiration_c14(ic, is)        = het_respiration_c14(ic, is) + (1._wp - frac_litter2fast_som) * &
          &                                                          (sb_loss_mt(ix_soluable_litter, ixC14, ic, is) + &
          &                                                           sb_loss_mt(ix_polymeric_litter, ixC14, ic, is))
        net_mineralisation_nh4(ic, is)     = net_mineralisation_nh4(ic, is) - microbial_uptake_nh4_sl(ic, is) + &
          &                                                             (1.0_wp - microbial_nue) * &
          &                                                             (sb_loss_mt(ix_soluable_litter, ixN, ic, is) + &
          &                                                              sb_loss_mt(ix_polymeric_litter, ixN, ic, is))
        net_mineralisation_no3(ic, is)     = net_mineralisation_no3(ic, is) - microbial_uptake_no3_sl(ic, is)
        net_mineralisation_po4(ic, is)     = net_mineralisation_po4(ic, is) - microbial_uptake_po4_sl(ic, is) + &
          &                                                             (1.0_wp - microbial_pue) * &
          &                                                             (sb_loss_mt(ix_soluable_litter, ixP, ic, is) + &
          &                                                              sb_loss_mt(ix_polymeric_litter, ixP, ic, is))
        net_mineralisation_nh4_n15(ic, is) = net_mineralisation_nh4_n15(ic, is) - microbial_uptake_nh4_n15_sl(ic, is) &
          &                                  + hlp1(ic, is) + hlp2(ic, is)
        net_mineralisation_no3_n15(ic, is) = net_mineralisation_no3_n15(ic, is) - microbial_uptake_no3_n15_sl(ic, is)

        !>2.0 formation of slow SOM (= residue pool) from fast SOM decay
        !>
        sb_formation_mt(ix_residue, ixC, ic, is)   = sb_formation_mt(ix_residue, ixC, ic, is) &
          &                                        + frac_fast2slow_som &
          &                                        * sb_loss_mt(ix_microbial, ixC, ic, is)
        sb_formation_mt(ix_residue, ixN, ic, is)   = sb_formation_mt(ix_residue, ixN, ic, is) &
          &                                        + frac_fast2slow_som / k_slow_som_cn &
          &                                        * sb_loss_mt(ix_microbial, ixC, ic, is)
        sb_formation_mt(ix_residue, ixP, ic, is)   = sb_formation_mt(ix_residue, ixP, ic, is) &
          &                                        + frac_fast2slow_som / k_slow_som_cn / k_slow_som_np &
          &                                        * sb_loss_mt(ix_microbial, ixC, ic, is)
        sb_formation_mt(ix_residue, ixC13, ic, is) = sb_formation_mt(ix_residue, ixC13, ic, is) &
          &                                        + frac_fast2slow_som &
          &                                        * sb_loss_mt(ix_microbial, ixC13, ic, is)
        sb_formation_mt(ix_residue, ixC14, ic, is) = sb_formation_mt(ix_residue, ixC14, ic, is) &
          &                                        + frac_fast2slow_som &
          &                                        * sb_loss_mt(ix_microbial, ixC14, ic, is)
        ! store formation of slow 15N for later use
        hlp1(ic, is) = frac_fast2slow_som / k_slow_som_cn * sb_loss_mt(ix_microbial, ixC, ic, is) &
          &         * calc_fractionation(sb_loss_mt(ix_microbial, ixN, ic, is), &
          &                              sb_loss_mt(ix_microbial, ixN15, ic, is), &
          &                              0.0_wp)
        sb_formation_mt(ix_residue, ixN15, ic, is) = sb_formation_mt(ix_residue, ixN15, ic, is) + hlp1(ic, is)

        het_respiration(ic, is)        = het_respiration(ic, is)     + (1._wp - frac_fast2slow_som) &
          &                              * sb_loss_mt(ix_microbial, ixC, ic, is)
        het_respiration_c13(ic, is)    = het_respiration_c13(ic, is) + (1._wp - frac_fast2slow_som) &
          &                              * sb_loss_mt(ix_microbial, ixC13, ic, is)
        het_respiration_c14(ic, is)    = het_respiration_c14(ic, is) + (1._wp - frac_fast2slow_som) &
          &                              * sb_loss_mt(ix_microbial, ixC14, ic, is)
        net_mineralisation_nh4(ic, is) = net_mineralisation_nh4(ic, is) + sb_loss_mt(ix_microbial, ixN, ic, is) &
          &                           - frac_fast2slow_som / k_slow_som_cn * sb_loss_mt(ix_microbial, ixC, ic, is)
        net_mineralisation_po4(ic, is) = net_mineralisation_po4(ic, is) + sb_loss_mt(ix_microbial, ixP, ic, is) &
          &                           - frac_fast2slow_som / k_slow_som_cn / k_slow_som_np &
          &                           * sb_loss_mt(ix_microbial, ixC, ic, is)
        ! store net mineralisation of 15N for later use
        hlp2(ic, is) = (sb_loss_mt(ix_microbial, ixN, ic, is) - frac_fast2slow_som / k_slow_som_cn &
          &         * sb_loss_mt(ix_microbial, ixC, ic, is)) &
          &         * calc_fractionation(sb_loss_mt(ix_microbial, ixN, ic, is), &
          &                              sb_loss_mt(ix_microbial, ixN15, ic, is), &
          &                              eta_ammonification)
        net_mineralisation_nh4_n15(ic, is) = net_mineralisation_nh4_n15(ic, is) + hlp2(ic, is)
        net_mineralisation_no3_n15(ic, is) = net_mineralisation_no3_n15(ic, is)

        ! adjust fast 15N loss rate to account for discrimination during ammonification
        sb_loss_mt(ix_microbial, ixN15, ic, is) = sb_loss_mt(ix_microbial, ixN15, ic, is) &
          &                                     - (sb_loss_mt(ix_microbial, ixN15, ic, is) - hlp1(ic, is) - hlp2(ic, is))

        !>3.0 formation of fast SOM (= microbial pool) from slow SOM decay
        !>
        sb_formation_mt(ix_microbial, ixC, ic, is)   = sb_formation_mt(ix_microbial, ixC, ic, is) &
          &                                          + frac_slow2fast_som &
          &                                          * sb_loss_mt(ix_residue, ixC, ic, is)
        sb_formation_mt(ix_microbial, ixN, ic, is)   = sb_formation_mt(ix_microbial, ixN, ic, is) &
          &                                          + frac_slow2fast_som / apparent_fast_som_cn(ic, is) &
          &                                          * sb_loss_mt(ix_residue, ixC, ic, is)
        sb_formation_mt(ix_microbial, ixP, ic, is)   = sb_formation_mt(ix_microbial, ixP, ic, is) &
          &                                          + frac_slow2fast_som / apparent_fast_som_cn(ic, is) &
          &                                          / apparent_fast_som_np(ic, is) &
          &                                          * sb_loss_mt(ix_residue, ixC, ic, is)
        sb_formation_mt(ix_microbial, ixC13, ic, is) = sb_formation_mt(ix_microbial, ixC13, ic, is) &
          &                                          + frac_slow2fast_som &
          &                                          * sb_loss_mt(ix_residue, ixC13, ic, is)
        sb_formation_mt(ix_microbial, ixC14, ic, is) = sb_formation_mt(ix_microbial, ixC14, ic, is) &
          &                                          + frac_slow2fast_som &
          &                                          * sb_loss_mt(ix_residue, ixC14, ic, is)
        ! store formation of fast 15N for later use
        hlp1(ic, is) = frac_slow2fast_som / apparent_fast_som_cn(ic, is) * sb_loss_mt(ix_residue, ixC, ic, is) &
          &         * calc_fractionation(sb_loss_mt(ix_residue, ixN, ic, is), &
          &                              sb_loss_mt(ix_residue, ixN15, ic, is), &
          &                              0.0_wp)
        sb_formation_mt(ix_microbial, ixN15, ic, is) = sb_formation_mt(ix_microbial, ixN15, ic, is) + hlp1(ic, is)

        het_respiration(ic, is)        = het_respiration(ic, is)     + (1._wp - frac_slow2fast_som) &
          &                              * sb_loss_mt(ix_residue, ixC, ic, is)
        het_respiration_c13(ic, is)    = het_respiration_c13(ic, is) + (1._wp - frac_slow2fast_som) &
          &                              * sb_loss_mt(ix_residue, ixC13, ic, is)
        het_respiration_c14(ic, is)    = het_respiration_c14(ic, is) + (1._wp - frac_slow2fast_som) &
          &                              * sb_loss_mt(ix_residue, ixC14, ic, is)
        net_mineralisation_nh4(ic, is) = net_mineralisation_nh4(ic, is) + sb_loss_mt(ix_residue, ixN, ic, is) &
          &                           - frac_slow2fast_som / apparent_fast_som_cn(ic, is) * sb_loss_mt(ix_residue, ixC, ic, is)
        ! update the loss of residue phosphorus due to biochem_mineralisation_po4, restricted by the CNP ratio of slow SOM
        biochem_mineralisation_po4(ic, is) = MIN(biochem_mineralisation_po4(ic, is), sb_loss_mt(ix_residue, ixP, ic, is) &
          &                               - frac_slow2fast_som / apparent_fast_som_cn(ic, is) / apparent_fast_som_np(ic, is) &
          &                               * sb_loss_mt(ix_residue, ixC, ic, is))
        net_mineralisation_po4(ic, is)     = net_mineralisation_po4(ic, is) + sb_loss_mt(ix_residue, ixP, ic, is) &
          &                               - biochem_mineralisation_po4(ic, is) &
          &                               - frac_slow2fast_som / apparent_fast_som_cn(ic, is) / apparent_fast_som_np(ic, is) &
          &                               * sb_loss_mt(ix_residue, ixC, ic, is)
        ! store net mineralisation of 15N for later use
        hlp2(ic, is) = (sb_loss_mt(ix_residue, ixN, ic, is) - frac_slow2fast_som / apparent_fast_som_cn(ic, is) &
          &         * sb_loss_mt(ix_residue, ixC, ic, is)) &
          &         * calc_fractionation(sb_loss_mt(ix_residue, ixN, ic, is), &
          &                              sb_loss_mt(ix_residue, ixN15, ic, is), &
          &                              eta_ammonification)
        net_mineralisation_nh4_n15(ic, is) = net_mineralisation_nh4_n15(ic, is) + hlp2(ic, is)

        ! adjust slow 15N loss rate to account for discrimination during ammonification
        sb_loss_mt(ix_residue, ixN15, ic, is) = sb_loss_mt(ix_residue, ixN15, ic, is) &
          &                                   - (sb_loss_mt(ix_residue, ixN15, ic, is) - hlp1(ic, is) - hlp2(ic, is))
      ENDDO
    ENDDO
  END SUBROUTINE calc_SOM_formation

  ! ======================================================================================================= !
  !>Calculates the formation of metabolic and structural litter from woody litter decay.
  !>
  !>  Some C,N,P becomes lost during processing prescribed by the carbon and nutrient use
  !>  efficiencies. The remainder becomes incorporated into metabolic and structural litter
  !>  given the remaining litter C:N and the prescribed lignin content.
  !>
  SUBROUTINE calc_litter_formation_from_woody_decay( &
    & nc, &
    & nsoil_sb, &
    & sb_loss_mt, &
    & het_respiration, &
    & het_respiration_c13, &
    & het_respiration_c14, &
    & net_mineralisation_nh4, &
    & net_mineralisation_po4, &
    & net_mineralisation_nh4_n15, &
    & sb_formation_mt )

    USE mo_sb_constants,        ONLY: frac_woody2dec_litter, microbial_nue, microbial_pue, &
      &                               fc_soluable_max, k_fc_soluable, k_fn_soluable, k_fp_soluable, lc_woody_litter
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc, &                         !< dimensions
                                                        nsoil_sb
    REAL(wp),                          INTENT(in)    :: sb_loss_mt(:,:,:,:)           !< bgcm sb_loss flux: turnover of SOM and litter pools  [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: het_respiration, &            !< heterotrophic respiration [mol C / m3 / timestep]
                                                        het_respiration_c13, &        !< heterotrophic respiration [mol 13C / m3 / timestep]
                                                        het_respiration_c14, &        !< heterotrophic respiration [micro-mol 14C / m3 / timestep]
                                                        net_mineralisation_nh4, &     !< net mineralisation of NH4 [mol / m3 / timestep]
                                                        net_mineralisation_po4, &     !< net mineralisation of PO4 [mol / m3 / timestep]
                                                        net_mineralisation_nh4_n15    !< net mineralisation of 15NH4 [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)      !< bgcm sb_formation flux: formation of SOM and litter pools  [mol / m3 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_sb) :: fc_soluable_woody_litter, & !< fraction of non-respired woody litter becoming metabolic litter
                                         fn_soluable_woody_litter, &
                                         fp_soluable_woody_litter
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_formation_from_woody_decay'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 Calculate partitioning of decaying woody litter into metabolic and structural litter
    !>
    fc_soluable_woody_litter(:,:) = 0.0_wp
    fn_soluable_woody_litter(:,:) = 0.0_wp
    fp_soluable_woody_litter(:,:) = 0.0_wp
    WHERE(sb_loss_mt(ix_woody_litter, ixN, :, :) > eps8)
      fc_soluable_woody_litter(:,:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_woody_litter &
        &                             * sb_loss_mt(ix_woody_litter, ixC, :, :) / sb_loss_mt(ix_woody_litter, ixN, :, :))
      WHERE(fc_soluable_woody_litter(:,:) > eps8)
        fn_soluable_woody_litter(:,:) = 1._wp / (1._wp + (1._wp - fc_soluable_woody_litter(:,:)) &
          &                             / (k_fn_soluable * fc_soluable_woody_litter(:,:)))
        fp_soluable_woody_litter(:,:) = 1._wp / (1._wp + (1._wp - fc_soluable_woody_litter(:,:)) &
          &                             / (k_fp_soluable * fc_soluable_woody_litter(:,:)))
      ENDWHERE
    ENDWHERE

    !>2.0 formation of metabolic and structural litter
    !>     taking account of the processing losses of C, N and P
    !>
    sb_formation_mt(ix_soluable_litter, ixC, :, :)   = sb_formation_mt(ix_soluable_litter, ixC, :, :) &
      &                                                + fc_soluable_woody_litter(:,:) * frac_woody2dec_litter &
      &                                                * sb_loss_mt(ix_woody_litter, ixC, :, :)
    sb_formation_mt(ix_soluable_litter, ixN, :, :)   = sb_formation_mt(ix_soluable_litter, ixN, :, :) &
      &                                                + fn_soluable_woody_litter(:,:) * microbial_nue &
      &                                                * sb_loss_mt(ix_woody_litter, ixN, :, :)
    sb_formation_mt(ix_soluable_litter, ixP, :, :)   = sb_formation_mt(ix_soluable_litter, ixP, :, :) &
      &                                                + fp_soluable_woody_litter(:,:) * microbial_pue &
      &                                                * sb_loss_mt(ix_woody_litter, ixP, :, :)
    sb_formation_mt(ix_soluable_litter, ixC13, :, :) = sb_formation_mt(ix_soluable_litter, ixC13, :, :) &
      &                                                + fc_soluable_woody_litter(:,:) * frac_woody2dec_litter &
      &                                                * sb_loss_mt(ix_woody_litter, ixC13, :, :)
    sb_formation_mt(ix_soluable_litter, ixC14, :, :) = sb_formation_mt(ix_soluable_litter, ixC14, :, :) &
      &                                                + fc_soluable_woody_litter(:,:) * frac_woody2dec_litter &
      &                                                * sb_loss_mt(ix_woody_litter, ixC14, :, :)
    sb_formation_mt(ix_soluable_litter, ixN15, :, :) = sb_formation_mt(ix_soluable_litter, ixN15, :, :) &
      &                                                + fn_soluable_woody_litter(:,:) * microbial_nue &
      &                                                * sb_loss_mt(ix_woody_litter, ixN15, :, :)

    sb_formation_mt(ix_polymeric_litter, ixC, :, :)   = sb_formation_mt(ix_polymeric_litter, ixC, :, :) &
      &                                                 + (1.0_wp - fc_soluable_woody_litter(:,:)) * frac_woody2dec_litter &
      &                                                 * sb_loss_mt(ix_woody_litter, ixC, :, :)
    sb_formation_mt(ix_polymeric_litter, ixN, :, :)   = sb_formation_mt(ix_polymeric_litter, ixN, :, :) &
      &                                                 + (1.0_wp - fn_soluable_woody_litter(:,:)) * microbial_nue &
      &                                                 * sb_loss_mt(ix_woody_litter, ixN, :, :)
    sb_formation_mt(ix_polymeric_litter, ixP, :, :)   = sb_formation_mt(ix_polymeric_litter, ixP, :, :) &
      &                                                 + (1.0_wp - fp_soluable_woody_litter(:,:)) * microbial_pue &
      &                                                 * sb_loss_mt(ix_woody_litter, ixP, :, :)
    sb_formation_mt(ix_polymeric_litter, ixC13, :, :) = sb_formation_mt(ix_polymeric_litter, ixC13, :, :) &
      &                                                 + (1.0_wp - fc_soluable_woody_litter(:,:)) * frac_woody2dec_litter &
      &                                                 * sb_loss_mt(ix_woody_litter, ixC13, :, :)
    sb_formation_mt(ix_polymeric_litter, ixC14, :, :) = sb_formation_mt(ix_polymeric_litter, ixC14, :, :) &
      &                                                 + (1.0_wp - fc_soluable_woody_litter(:,:)) * frac_woody2dec_litter &
      &                                                 * sb_loss_mt(ix_woody_litter, ixC14, :, :)
    sb_formation_mt(ix_polymeric_litter, ixN15, :, :) = sb_formation_mt(ix_polymeric_litter, ixN15, :, :) &
      &                                                 + (1.0_wp - fn_soluable_woody_litter(:,:)) * microbial_nue &
      &                                                 * sb_loss_mt(ix_woody_litter, ixN15, :, :)

    !>3.0 record losses to heterotrophic respiration and mineralisation
    !>
    het_respiration(:,:)            = het_respiration(:,:) &
      &                               + (1.0_wp - frac_woody2dec_litter) &
      &                               * sb_loss_mt(ix_woody_litter, ixC, :, :)
    net_mineralisation_nh4(:,:)     = net_mineralisation_nh4(:,:) &
      &                               + (1.0_wp - microbial_nue) &
      &                               * sb_loss_mt(ix_woody_litter, ixN, :, :)
    net_mineralisation_po4(:,:)     = net_mineralisation_po4(:,:) &
      &                               + (1.0_wp - microbial_pue) &
      &                               * sb_loss_mt(ix_woody_litter, ixP, :, :)
    het_respiration_c13(:,:)        = het_respiration_c13(:,:) &
      &                               + (1.0_wp - frac_woody2dec_litter) &
      &                               * sb_loss_mt(ix_woody_litter, ixC13, :, :)
    het_respiration_c14(:,:)        = het_respiration_c14(:,:) &
      &                               + (1.0_wp - frac_woody2dec_litter) &
      &                               * sb_loss_mt(ix_woody_litter, ixC14, :, :)
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) &
      &                               + (1.0_wp - microbial_nue) &
      &                               * sb_loss_mt(ix_woody_litter, ixN15, :, :)
  END SUBROUTINE calc_litter_formation_from_woody_decay

  ! ======================================================================================================= !
  !>calculates matter transport by bioturbation, by calling calc_bioturbation_transport
  !>
  !> Input:  pools subject to bioturbation, percolation rate, and soil depth
  !>         This explicitly excludes mycorrhiza, as these die when eaten
  !>
  !> Output: vertical transport rate (mol/m2/timestep)
  !>
  SUBROUTINE calc_bioturbation_transport_wrapper_s1d( &
    & nc, &
    & nsoil_sb, &
    & dtime, &
    & num_sl_above_bedrock, &
    & soil_depth_sl, &
    & elements_index_map, &
    & is_element_used, &
    & k_bioturb, &
    & nh4_assoc, &
    & nh4_n15_assoc, &
    & po4_assoc_fast, &
    & po4_assoc_slow, &
    & po4_occluded, &
    & sb_pool_mt, &
    & transport_nh4_assoc, &
    & transport_nh4_n15_assoc, &
    & transport_po4_assoc_fast, &
    & transport_po4_assoc_slow, &
    & transport_po4_occluded, &
    & sb_transport_mt)

    USE mo_q_sb_jsm_transport,      ONLY: calc_bioturbation_transport
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                          !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                    !< number of soil layers
    REAL(wp),                          INTENT(in)    :: dtime                       !< timestep length
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock        !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: soil_depth_sl               !< depth of each soil layer [m]
    INTEGER,                           INTENT(in)    :: elements_index_map(:)       !< map bgcm element ID -> IDX
    LOGICAL,                           INTENT(in)    :: is_element_used(:)          !< is element in 'elements_index_map' used
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: k_bioturb                   !< diffusion factor for bioturbation
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: nh4_assoc                   !< adsorbed NH4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: nh4_n15_assoc               !< adsorbed 15NH4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: po4_assoc_fast              !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: po4_assoc_slow              !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: po4_occluded                !< occluded PO4 pool [mol/m3]
    REAL(wp),                          INTENT(in)    :: sb_pool_mt(:, :, :, :)      !< sb_pool
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_nh4_assoc         !< rate of associated NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_nh4_n15_assoc     !< rate of associated 15NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_po4_assoc_fast    !< rate of fast associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_po4_assoc_slow    !< rate of slow associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_po4_occluded      !< rate of occluded PO4 transport [mol/m3/timestep]
    REAL(wp),                          INTENT(inout) :: sb_transport_mt(:, :, :, :) !< sb_transport flux
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ielem             !< loop over bgcm elements
    INTEGER                     :: ix_elem           !< index of element in bgcm, used for looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_transport_wrapper_s1d'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 organic pools
    !>
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        !>  1.1 fast SOM pool
        !>
        CALL calc_bioturbation_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), &
          &                              sb_pool_mt(ix_soluable_litter, ix_elem, :, :), &
          &                              sb_transport_mt(ix_soluable_litter, ix_elem, :, :))
        CALL calc_bioturbation_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), &
          &                              sb_pool_mt(ix_polymeric_litter, ix_elem, :, :), &
          &                              sb_transport_mt(ix_polymeric_litter, ix_elem, :, :))
        CALL calc_bioturbation_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), &
          &                              sb_pool_mt(ix_microbial, ix_elem, :, :), &
          &                              sb_transport_mt(ix_microbial, ix_elem, :, :))
        !>  1.1 slow SOM pool
        !>
        CALL calc_bioturbation_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), &
          &                              sb_pool_mt(ix_residue, ix_elem, :, :), &
          &                              sb_transport_mt(ix_residue, ix_elem, :, :))
      END IF
    END DO

    !>2.0 inorganic pools
    !>
    CALL calc_bioturbation_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
      &                              k_bioturb(:,:),&
      &                              soil_depth_sl(:,:), &
      &                              nh4_assoc(:,:), transport_nh4_assoc(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
      &                              k_bioturb(:,:),&
      &                              soil_depth_sl(:,:), &
      &                              nh4_n15_assoc(:,:), transport_nh4_n15_assoc(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
      &                              k_bioturb(:,:),&
      &                              soil_depth_sl(:,:), &
      &                              po4_assoc_fast(:,:), transport_po4_assoc_fast(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
      &                              k_bioturb(:,:),&
      &                              soil_depth_sl(:,:), &
      &                              po4_assoc_slow(:,:), transport_po4_assoc_slow(:,:))
    CALL calc_bioturbation_transport(nc, nsoil_sb, dtime, num_sl_above_bedrock(:), &
      &                              k_bioturb(:,:),&
      &                              soil_depth_sl(:,:), &
      &                              po4_occluded(:,:), transport_po4_occluded(:,:))
  END SUBROUTINE calc_bioturbation_transport_wrapper_s1d

  ! ======================================================================================================= !
  !>For the sb_nloss_scheme = fixed: calculates the amount of N lost in proportion to net mineralisation
  !>
  SUBROUTINE calc_fixed_nloss_rate( &
    & nc, &
    & nsoil_sb, &
    & soil_depth_sl, &
    & net_mineralisation_nh4, &
    & net_mineralisation_nh4_n15, &
    & nitrification_no3, &
    & nitrification_no3_n15, &
    & volatilisation_nh4, &
    & volatilisation_nh4_n15, &
    & emission_n2, &
    & emission_n2_n15)

    USE mo_sb_constants,           ONLY: floss_nmin, fnit_nmin
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                                          INTENT(in)  :: nc                            !< dimensions
    INTEGER,                                          INTENT(in)  :: nsoil_sb                      !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: soil_depth_sl, &              !< soil depth per layer [m]
                                                                     net_mineralisation_nh4, &     !< net mineralisation rate of NH4 [mol / m3 / timestep]
                                                                     net_mineralisation_nh4_n15    !< net mineralisation rate of 15NH4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out) :: nitrification_no3, &          !< nitrification to NO3 [mol / m3 / timestep]
                                                                     nitrification_no3_n15, &      !< nitrification to 15NO3 [mol / m3 / timestep]
                                                                     volatilisation_nh4, &         !< volatilisation of NH4 [mol / m3 / timestep]
                                                                     volatilisation_nh4_n15        !< volatilisation of 15NH4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc),                          INTENT(out) :: emission_n2, &                !< soil efflux of N2 [mol / m2 / timestep]
                                                                     emission_n2_n15               !< soil efflux of 15N2 [mol / m2 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fixed_nloss_rate'
    ! ----------------------------------------------------------------------------------------------------- !

    !> init out var
    !>
    nitrification_no3(:,:)      = 0.0_wp
    nitrification_no3_n15(:,:)  = 0.0_wp
    volatilisation_nh4(:,:)     = 0.0_wp
    volatilisation_nh4_n15(:,:) = 0.0_wp
    emission_n2(:)              = 0.0_wp
    emission_n2_n15(:)          = 0.0_wp

    !>1.0 case of fixed N losses: N loss proportional to net mineralisation
    !>
    WHERE(net_mineralisation_nh4(:,:) > 0.0_wp)
      ! N lost to atmosphere during processing of net mineralisation [mol m-3 timestep-1]
      volatilisation_nh4(:,:)      = net_mineralisation_nh4(:,:)     * floss_nmin
      volatilisation_nh4_n15(:,:)  = net_mineralisation_nh4_n15(:,:) * floss_nmin
      ! NH4 transformed to NO3 [mol m-3 timestep-1]
      nitrification_no3(:,:)       = net_mineralisation_nh4(:,:)     * (1._wp - floss_nmin) * fnit_nmin
      nitrification_no3_n15(:,:)   = net_mineralisation_nh4_n15(:,:) * (1._wp - floss_nmin) * fnit_nmin
    ENDWHERE
    ! report volatilisation losses as N2 emission [mol m-2 timestep-1]
    emission_n2(:)     = SUM(volatilisation_nh4(:,:)     * soil_depth_sl(:,:), DIM=2)
    emission_n2_n15(:) = SUM(volatilisation_nh4_n15(:,:) * soil_depth_sl(:,:), DIM=2)
  END SUBROUTINE calc_fixed_nloss_rate


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to deal with mycorrhiza
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_mycorrhiza_dynamics( &
    & nc, &
    & nsoil_sb, &
    & dtime, &
    & soil_depth_sl, &
    & elements_index_map, &
    & is_element_used, &
    & lctlib_tau_mycorrhiza, &
    & flag_mycorrhiza_org, &
    & nh4_solute, &
    & no3_solute, &
    & km2_uptake_nh4_act, &
    & km2_uptake_no3_act, &
    & km_uptake_org_act, &
    & veg_pool_fine_root_carbon , &
    & sb_pool_mt, &
    & mycorrhiza_uptake_nh4_sl, &
    & mycorrhiza_uptake_nh4_n15_sl, &
    & mycorrhiza_uptake_no3_sl, &
    & mycorrhiza_uptake_no3_n15_sl, &
    & mycorrhiza_uptake_norg_sl, &
    & mycorrhiza_uptake_norg_n15_sl, &
    & mycorrhiza_uptake_po4_sl, &
    & het_respiration, &
    & het_respiration_c13, &
    & het_respiration_c14, &
    & myc_respiration, &
    & myc_respiration_c13, &
    & myc_respiration_c14, &
    & sb_formation_mt, &
    & sb_loss_mt, &
    & sb_mycorrhiza_export_mt)

    USE mo_veg_constants,           ONLY: transform_cost_nh4, transform_cost_no3
    USE mo_sb_constants,            ONLY: mycorrhiza_cn, mycorrhiza_np, mycorrhiza_cue, &
                                          tau_slow, vmax_uptake_norg, f_nin2norg, f_norg_limit, organic_cn
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                                !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                          !< number of soil layers
    REAL(wp),                          INTENT(in)    :: dtime                             !< timestep length
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: soil_depth_sl                     !< depth per soil layer
    INTEGER,                           INTENT(in)    :: elements_index_map(:)             !< map bgcm element ID -> IDX
    LOGICAL,                           INTENT(in)    :: is_element_used(:)                !< is element in 'elements_index_map' used
    REAL(wp),                          INTENT(in)    :: lctlib_tau_mycorrhiza             !< lctlib parameter
    LOGICAL,                           INTENT(in)    :: flag_mycorrhiza_org               !< from sb config
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: nh4_solute, &                     !< soil solute NH4 concentration [mol m-3]
                                                        no3_solute, &                     !< soil solute NO3 concentration [mol m-3]
                                                        km2_uptake_nh4_act, &             !< ..
                                                        km2_uptake_no3_act, &             !< ..
                                                        km_uptake_org_act                 !< ..
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: veg_pool_fine_root_carbon         !< fine root carbon [mol / m2]
    REAL(wp),                          INTENT(in)    :: sb_pool_mt(:,:,:,:)               !< bgcm sb_pool: soil biogeochemistry pools [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: mycorrhiza_uptake_nh4_sl, &       !< NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                        mycorrhiza_uptake_nh4_n15_sl, &   !< 15NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                        mycorrhiza_uptake_no3_sl, &       !< NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                        mycorrhiza_uptake_no3_n15_sl, &   !< 15NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                        mycorrhiza_uptake_norg_sl, &      !< N uptake by mycorrhiza from organic (residue) pool [mol / m3 / timestep]
                                                        mycorrhiza_uptake_norg_n15_sl, &  !< N uptake by mycorrhiza from organic (residue) pool [mol / m3 / timestep]
                                                        mycorrhiza_uptake_po4_sl, &       !< PO4 uptake by mycorrhiza [mol / m3 / timestep]
                                                        het_respiration, &                !< heterotrophic respiration [mol / m3 / timestep]
                                                        het_respiration_c13, &            !< heterotrophic respiration [mol / m3 / timestep]
                                                        het_respiration_c14, &            !< heterotrophic respiration [mol / m3 / timestep]
                                                        myc_respiration, &                !< mycorrhizae respiration [mol / m3 / timestep]
                                                        myc_respiration_c13, &            !< mycorrhizae respiration [mol / m3 / timestep]
                                                        myc_respiration_c14               !< mycorrhizae respiration [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)          !< bgcm sb_formation flux: formation rate  [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_loss_mt(:,:,:,:)               !< bgcm sb_loss flux: loss rate  [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_mycorrhiza_export_mt(:,:,:)    !< bgcm sb_mycorrhiza_export flux: export from mycorrhiza to plants [mol / m2 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: transform_respiration, &
                                                        fresp, &
                                                        fturn_mycorrhiza, &
                                                        fvmax, &
                                                        vmax_norg_up, &
                                                        k_norg_up, &
                                                        norg_act, &
                                                        fnorg, &
                                                        f_myc_opt
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1, hlp2, hlp3
    INTEGER                                          :: ielem               !< loop over bgcm elements
    INTEGER                                          :: ix_elem             !< index of element in bgcm, used for looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_mycorrhiza_dynamics'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 Mycorrhiza dynamics
    !>

    !>  1.1 uptake of nutrients by mycorrhiza and associated respiration losses to heterotrophic respiration
    !>
    WHERE(sb_pool_mt(ix_mycorrhiza, ixN, :, :) > eps8)
      hlp1(:,:) = MIN(1.0_wp, MAX(0.0_wp, &
        & (((sb_pool_mt(ix_mycorrhiza, ixC, :, :) / sb_pool_mt(ix_mycorrhiza, ixN, :, :)) &
        & - mycorrhiza_cn * f_norg_limit) &
        & / (mycorrhiza_cn * (1._wp - f_norg_limit)))))
      mycorrhiza_uptake_nh4_sl(:,:)     = mycorrhiza_uptake_nh4_sl(:,:)     * hlp1(:,:)
      mycorrhiza_uptake_nh4_n15_sl(:,:) = mycorrhiza_uptake_nh4_n15_sl(:,:) * hlp1(:,:)
      mycorrhiza_uptake_no3_sl(:,:)     = mycorrhiza_uptake_no3_sl(:,:)     * hlp1(:,:)
      mycorrhiza_uptake_no3_n15_sl(:,:) = mycorrhiza_uptake_no3_n15_sl(:,:) * hlp1(:,:)
    ENDWHERE

    WHERE(sb_pool_mt(ix_mycorrhiza, ixP, :, :) > eps8)
      hlp1(:,:) = MIN(1.0_wp,MAX(0.0_wp, &
        & (((sb_pool_mt(ix_mycorrhiza, ixN, :, :) / sb_pool_mt(ix_mycorrhiza, ixP, :, :)) &
        & - mycorrhiza_np * f_norg_limit) &
        & / (mycorrhiza_np * (1._wp - f_norg_limit)))))
      mycorrhiza_uptake_po4_sl(:,:) = mycorrhiza_uptake_po4_sl(:,:) * hlp1(:,:)
    ENDWHERE

    ! calculate respiration costs associated with N transformation (mol C / m3 / timestep)
    transform_respiration(:,:) = mycorrhiza_uptake_nh4_sl(:,:) * transform_cost_nh4 &
      &                          + mycorrhiza_uptake_no3_sl(:,:) * transform_cost_no3

    fresp(:,:) = 0.0_wp
    WHERE(sb_pool_mt(ix_mycorrhiza, ixC, :, :) > eps8)
      fresp(:,:) = transform_respiration(:,:) / sb_pool_mt(ix_mycorrhiza, ixC, :, :)
    ENDWHERE

    sb_loss_mt(ix_mycorrhiza, ixC, :, :)   = sb_loss_mt(ix_mycorrhiza, ixC, :, :)   + fresp(:,:) &
      &                                      * sb_pool_mt(ix_mycorrhiza, ixC, :, :)   !! respiration
    sb_loss_mt(ix_mycorrhiza, ixC13, :, :) = sb_loss_mt(ix_mycorrhiza, ixC13, :, :) + fresp(:,:) &
      &                                      * sb_pool_mt(ix_mycorrhiza, ixC13, :, :) !! respiration
    sb_loss_mt(ix_mycorrhiza, ixC14, :, :) = sb_loss_mt(ix_mycorrhiza, ixC14, :, :) + fresp(:,:) &
      &                                      * sb_pool_mt(ix_mycorrhiza, ixC14, :, :) !! respiration

    het_respiration(:,:)     = het_respiration(:,:)     + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC13, :, :)
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC14, :, :)

    myc_respiration(:,:)     = myc_respiration(:,:)     + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
    myc_respiration_c13(:,:) = myc_respiration_c13(:,:) + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC13, :, :)
    myc_respiration_c14(:,:) = myc_respiration_c14(:,:) + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC14, :, :)

    !>  1.2 organic N uptake by mycorrhiza and associated respiration losses to heterotrophic respiration
    !>
    IF (flag_mycorrhiza_org) THEN
      ! limit organic uptake when inorganic nitrogen is abundant
      fvmax(:,:) = MIN(1.0_wp, MAX(0.0_wp, &
        & 0.5_wp * (2.0_wp - nh4_solute(:,:) / (f_nin2norg * km2_uptake_nh4_act(:,:) + nh4_solute(:,:)) &
        & - no3_solute(:,:) / (f_nin2norg * km2_uptake_no3_act(:,:) + no3_solute(:,:)))))
      hlp1(:,:) = fvmax(:,:)

      ! limit organic uptake according to accessibility to substrate
      fvmax(:,:) = fvmax(:,:) * sb_pool_mt(ix_residue, ixN, :, :) / ( km_uptake_org_act(:,:) + sb_pool_mt(ix_residue, ixN, :, :))
      mycorrhiza_uptake_norg_sl(:,:) = fvmax(:,:) * vmax_uptake_norg * sb_pool_mt(ix_mycorrhiza, ixC, :, :) &
        &                              * dtime / 1000000._wp

      WHERE(sb_pool_mt(ix_residue, ixN, :, :) > eps8)
        mycorrhiza_uptake_norg_n15_sl(:,:) = mycorrhiza_uptake_norg_sl(:,:) &
          &                                  * sb_pool_mt(ix_residue, ixN15, :, :) / sb_pool_mt(ix_residue, ixN, :, :)
      ENDWHERE

      ! reduce organic uptake at low mycorrhiza CN (should be obsolete with above function)
      WHERE(sb_pool_mt(ix_mycorrhiza, ixN, :, :) > eps8)
        hlp1(:,:) = MIN(1.0_wp, MAX(0.0_wp, &
          & (((sb_pool_mt(ix_mycorrhiza, ixC, :, :) / sb_pool_mt(ix_mycorrhiza, ixN, :, :)) &
          & - mycorrhiza_cn * f_norg_limit) &
          & / (mycorrhiza_cn * (1._wp - f_norg_limit)))))
        mycorrhiza_uptake_norg_sl(:,:)     = mycorrhiza_uptake_norg_sl(:,:)     * hlp1(:,:)
        mycorrhiza_uptake_norg_n15_sl(:,:) = mycorrhiza_uptake_norg_n15_sl(:,:) * hlp1(:,:)
      ENDWHERE

      ! calculate N uptake fraction
      fnorg(:,:) = 0.0_wp
      WHERE(sb_pool_mt(ix_residue, ixN, :, :) > eps8)
        fnorg(:,:) = mycorrhiza_uptake_norg_sl(:,:) / sb_pool_mt(ix_residue, ixN, :, :)
      ENDWHERE

      sb_loss_mt(ix_residue, ixC, :, :)   = sb_loss_mt(ix_residue, ixC, :, :)   + fnorg(:,:) * sb_pool_mt(ix_residue, ixC, :, :)
      sb_loss_mt(ix_residue, ixN, :, :)   = sb_loss_mt(ix_residue, ixN, :, :)   + fnorg(:,:) * sb_pool_mt(ix_residue, ixN, :, :)
      sb_loss_mt(ix_residue, ixP, :, :)   = sb_loss_mt(ix_residue, ixP, :, :)   + fnorg(:,:) * sb_pool_mt(ix_residue, ixP, :, :)
      sb_loss_mt(ix_residue, ixC13, :, :) = sb_loss_mt(ix_residue, ixC13, :, :) + fnorg(:,:) * sb_pool_mt(ix_residue, ixC13, :, :)
      sb_loss_mt(ix_residue, ixC14, :, :) = sb_loss_mt(ix_residue, ixC14, :, :) + fnorg(:,:) * sb_pool_mt(ix_residue, ixC14, :, :)
      sb_loss_mt(ix_residue, ixN15, :, :) = sb_loss_mt(ix_residue, ixN15, :, :) + fnorg(:,:) * sb_pool_mt(ix_residue, ixN15, :, :)

      !new mycorrhiza formation out of uptake from slow N pool
      sb_formation_mt(ix_mycorrhiza, ixC, :, :)   = sb_formation_mt(ix_mycorrhiza, ixC, :, :) &
        &                                           + mycorrhiza_cue * fnorg(:,:) * sb_pool_mt(ix_residue, ixC, :, :)
      sb_formation_mt(ix_mycorrhiza, ixN, :, :)   = sb_formation_mt(ix_mycorrhiza, ixN, :, :) &
        &                                           + fnorg(:,:) * sb_pool_mt(ix_residue, ixN, :, :)
      sb_formation_mt(ix_mycorrhiza, ixP, :, :)   = sb_formation_mt(ix_mycorrhiza, ixP, :, :) &
        &                                           + fnorg(:,:) * sb_pool_mt(ix_residue, ixP, :, :)
      sb_formation_mt(ix_mycorrhiza, ixC13, :, :) = sb_formation_mt(ix_mycorrhiza, ixC13, :, :) &
        &                                           + mycorrhiza_cue * fnorg(:,:) * sb_pool_mt(ix_residue, ixC13, :, :)
      sb_formation_mt(ix_mycorrhiza, ixC14, :, :) = sb_formation_mt(ix_mycorrhiza, ixC14, :, :) &
        &                                           + mycorrhiza_cue * fnorg(:,:) * sb_pool_mt(ix_residue, ixC14, :, :)
      sb_formation_mt(ix_mycorrhiza, ixN15, :, :) = sb_formation_mt(ix_mycorrhiza, ixN15, :, :) &
        &                                           + fnorg(:,:) * sb_pool_mt(ix_residue, ixN15, :, :)

      het_respiration(:,:)     = het_respiration(:,:)     + (1.0_wp - mycorrhiza_cue) &
        &                        * fnorg(:,:) * sb_pool_mt(ix_residue, ixC, :, :)
      het_respiration_c13(:,:) = het_respiration_c13(:,:) + (1.0_wp - mycorrhiza_cue) &
        &                        * fnorg(:,:) * sb_pool_mt(ix_residue, ixC13, :, :)
      het_respiration_c14(:,:) = het_respiration_c14(:,:) + (1.0_wp - mycorrhiza_cue) &
        &                        * fnorg(:,:) * sb_pool_mt(ix_residue, ixC14, :, :)
    ENDIF

    !>  1.3 nutrient flux to vegetation (surplus uptake goes to plants) [mol (N/P) / m3 / timestep]
    !>
    ! maintain C:N of mycorrhizal fungi
    ! N supply for plants = N uptake by mycorrhiza - N demand by mycorrhiza
    ! N demand can be negative => N supply can be greater then current N uptake by plants!

    ! N demand of mycorrhiza in mol N / m3 / timestep
    hlp1(:,:) = MAX(0.0_wp, &
      &         sb_pool_mt(ix_mycorrhiza, ixC, :, :) / mycorrhiza_cn - sb_pool_mt(ix_mycorrhiza, ixN, :, :))
    ! N uptake incorporated into mycorrhiza (minimum of potential uptake and demand)
    hlp2(:,:) = MIN(hlp1(:,:), &
      &         mycorrhiza_uptake_nh4_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:) + mycorrhiza_uptake_norg_sl(:,:))
    ! remainder is exported to plants (mol N / m3 / timestep)
    sb_mycorrhiza_export_mt( ixN, :, :) = sb_mycorrhiza_export_mt( ixN, :, :) &
      &                                   + mycorrhiza_uptake_nh4_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:) &
      &                                   + mycorrhiza_uptake_norg_sl(:,:) - hlp2(:,:)
    WHERE(mycorrhiza_uptake_nh4_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:) + mycorrhiza_uptake_norg_sl(:,:) > eps8)
      sb_mycorrhiza_export_mt( ixN15, :, :) = sb_mycorrhiza_export_mt( ixN15, :, :) &
        &                                     + sb_mycorrhiza_export_mt( ixN, :, :) * (mycorrhiza_uptake_nh4_n15_sl(:,:) &
        &                                     + mycorrhiza_uptake_no3_n15_sl(:,:) &
        &                                     + mycorrhiza_uptake_norg_n15_sl(:,:)) / (mycorrhiza_uptake_nh4_sl(:,:) &
        &                                     + mycorrhiza_uptake_no3_sl(:,:) &
        &                                     + mycorrhiza_uptake_norg_sl(:,:))
    ENDWHERE
    ! assuming that mycorrhiza export organic N (C:N = 3:1)
    sb_mycorrhiza_export_mt( ixC, :, :) = sb_mycorrhiza_export_mt( ixC, :, :) &
      &                                   + sb_mycorrhiza_export_mt( ixN, :, :) * organic_cn
    WHERE(sb_pool_mt(ix_mycorrhiza, ixC, :, :) > eps8)
      sb_mycorrhiza_export_mt( ixC13, :, :) = sb_mycorrhiza_export_mt( ixC13, :, :) &
        &                                     + sb_mycorrhiza_export_mt( ixC, :, :) &
        &                                     * sb_pool_mt(ix_mycorrhiza, ixC13, :, :) / sb_pool_mt(ix_mycorrhiza, ixC, :, :)
      sb_mycorrhiza_export_mt( ixC14, :, :) = sb_mycorrhiza_export_mt( ixC14, :, :) &
        &                                     + sb_mycorrhiza_export_mt( ixC, :, :) &
        &                                     * sb_pool_mt(ix_mycorrhiza, ixC14, :, :) / sb_pool_mt(ix_mycorrhiza, ixC, :, :)
    ENDWHERE

    ! P demand of mycorrhiza in mol P / m3 / timestep
    hlp1(:,:) = MAX(sb_pool_mt(ix_mycorrhiza, ixN, :, :) / mycorrhiza_np &
      &         - sb_pool_mt(ix_mycorrhiza, ixP, :, :), 0.0_wp)
    ! incorporation of P into mycorrhiza (mol P / m3 / timestep), constrain flux to minimum of potential uptake and demand
    hlp2(:,:) = MIN(hlp1(:,:), mycorrhiza_uptake_po4_sl(:,:))
    ! remainder is exported to plants (mol P / m3 / timestep)
    hlp3(:,:) = mycorrhiza_uptake_po4_sl(:,:) - hlp2(:,:)
    sb_mycorrhiza_export_mt( ixP, :, :) = sb_mycorrhiza_export_mt( ixP, :, :) + hlp3(:,:)

    !>  1.4 Formation of litter from mycorrhiza turnover (mol / m3 / timestep)
    !>
    ! calc turnover
    fturn_mycorrhiza(:,:) = 1.0_wp / lctlib_tau_mycorrhiza / one_day / one_year * dtime
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        ! ...
        sb_loss_mt(ix_mycorrhiza, ixC, :, :) = sb_loss_mt(ix_mycorrhiza, ixC, :, :) &
          &                                    + fturn_mycorrhiza(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
        ! ...
        sb_formation_mt(ix_polymeric_litter, ixC, :, :) = sb_formation_mt(ix_polymeric_litter, ixC, :, :) &
          &                                               + fturn_mycorrhiza(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
      END IF
    END DO
  END SUBROUTINE calc_mycorrhiza_dynamics

#endif
END MODULE mo_q_sb_ssm_main
