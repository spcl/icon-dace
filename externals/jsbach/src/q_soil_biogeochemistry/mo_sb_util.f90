!> QUINCY helper routines for soil-beogeochemistry
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
!>#### various helper routines for the soil-beogeochemistry process
!>
MODULE mo_sb_util
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_jsb_math_constants,    ONLY: zero

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, SB_BGCM_TRANSPORT_ID, SB_BGCM_MYCO_EXPORT_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: reset_sb_fluxes, calculate_time_average_soilbiogeochemistry

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_util'

CONTAINS

  ! ======================================================================================================= !
  !>reset sb fluxes (to zero)
  !>
  SUBROUTINE reset_sb_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)    :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                !< the model
    TYPE(t_lnd_bgcm_store), POINTER   :: bgcm_store           !< the bgcm store of this tile
    TYPE(t_jsb_vgrid),      POINTER   :: vgrid_soil_sb        !< Vertical grid
    INTEGER                           :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    INTEGER                           :: iblk, ics, ice, nc   !< dimensions
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_sb_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt
    dsl4jsb_Def_mt2L3D :: sb_formation_mt
    dsl4jsb_Def_mt2L3D :: sb_loss_mt
    dsl4jsb_Def_mt2L3D :: sb_transport_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    IF (model%lctlib(tile%lcts(1)%lib_id)%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FORMATION_ID, sb_formation_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_LOSS_ID, sb_loss_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_TRANSPORT_ID, sb_transport_mt)
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 bgcm fluxes
    !>
    ! set the bgcm matrices to zero
    sb_formation_mt(:,:,:,:)       = zero
    sb_loss_mt(:,:,:,:)            = zero
    sb_transport_mt(:,:,:,:)       = zero
    sb_mycorrhiza_export_mt(:,:,:) = zero

    !>2.0 variables
    !>
    ! 2D
    dsl4jsb_var2D_onChunk(SB_, total_flux_carbon)              = zero
    dsl4jsb_var2D_onChunk(SB_, total_flux_nitrogen)            = zero
    dsl4jsb_var2D_onChunk(SB_, total_flux_phosphorus)          = zero
    dsl4jsb_var2D_onChunk(SB_, total_flux_carbon13)            = zero
    dsl4jsb_var2D_onChunk(SB_, total_flux_carbon14)            = zero
    dsl4jsb_var2D_onChunk(SB_, total_flux_nitrogen15)          = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_nh4_solute)            = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_no3_solute)            = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_po4_solute)            = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_carbon)            = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_nitrogen)          = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_phosphorus)        = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_carbon13)          = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_carbon14)          = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_nitrogen15)        = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_nh4_n15_solute)        = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_no3_n15_solute)        = zero
    dsl4jsb_var2D_onChunk(SB_, emission_noy)                   = zero
    dsl4jsb_var2D_onChunk(SB_, emission_n2o)                   = zero
    dsl4jsb_var2D_onChunk(SB_, emission_n2)                    = zero
    dsl4jsb_var2D_onChunk(SB_, emission_noy_n15)               = zero
    dsl4jsb_var2D_onChunk(SB_, emission_n2o_n15)               = zero
    dsl4jsb_var2D_onChunk(SB_, emission_n2_n15)                = zero
    ! 3D
    dsl4jsb_var3D_onChunk(SB_, het_respiration)                  = zero
    dsl4jsb_var3D_onChunk(SB_, het_respiration_c13)              = zero
    dsl4jsb_var3D_onChunk(SB_, het_respiration_c14)              = zero
    dsl4jsb_var3D_onChunk(SB_, myc_respiration)                  = zero
    dsl4jsb_var3D_onChunk(SB_, myc_respiration_c13)              = zero
    dsl4jsb_var3D_onChunk(SB_, myc_respiration_c14)              = zero
    dsl4jsb_var3D_onChunk(SB_, net_mineralisation_nh4)           = zero
    dsl4jsb_var3D_onChunk(SB_, net_mineralisation_nh4_n15)       = zero
    dsl4jsb_var3D_onChunk(SB_, net_mineralisation_no3)           = zero
    dsl4jsb_var3D_onChunk(SB_, net_mineralisation_no3_n15)       = zero
    dsl4jsb_var3D_onChunk(SB_, net_mineralisation_po4)           = zero
    dsl4jsb_var3D_onChunk(SB_, plant_uptake_nh4_sl)              = zero
    dsl4jsb_var3D_onChunk(SB_, plant_uptake_nh4_n15_sl)          = zero
    dsl4jsb_var3D_onChunk(SB_, plant_uptake_no3_sl)              = zero
    dsl4jsb_var3D_onChunk(SB_, plant_uptake_no3_n15_sl)          = zero
    dsl4jsb_var3D_onChunk(SB_, plant_uptake_po4_sl)              = zero
    dsl4jsb_var3D_onChunk(SB_, microbial_uptake_nh4_sl)          = zero
    dsl4jsb_var3D_onChunk(SB_, microbial_uptake_nh4_n15_sl)      = zero
    dsl4jsb_var3D_onChunk(SB_, microbial_uptake_no3_sl)          = zero
    dsl4jsb_var3D_onChunk(SB_, microbial_uptake_no3_n15_sl)      = zero
    dsl4jsb_var3D_onChunk(SB_, microbial_uptake_po4_sl)          = zero
    dsl4jsb_var3D_onChunk(SB_, weathering_po4)                   = zero
    dsl4jsb_var3D_onChunk(SB_, occlusion_po4)                    = zero
    dsl4jsb_var3D_onChunk(SB_, slow_exchange_po4)                = zero
    dsl4jsb_var3D_onChunk(SB_, fast_exchange_po4)                = zero
    dsl4jsb_var3D_onChunk(SB_, fast_exchange_nh4)                = zero
    dsl4jsb_var3D_onChunk(SB_, fast_adsorpt_po4)                 = zero
    dsl4jsb_var3D_onChunk(SB_, fast_adsorpt_nh4)                 = zero
    dsl4jsb_var3D_onChunk(SB_, biochem_mineralisation_po4)       = zero
    dsl4jsb_var3D_onChunk(SB_, gross_mineralisation_po4)         = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_nh4_sl)         = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_nh4_n15_sl)     = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_no3_sl)         = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_no3_n15_sl)     = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_norg_sl)        = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_norg_n15_sl)    = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_po4_sl)         = zero
    dsl4jsb_var3D_onChunk(SB_, transport_nh4_solute)             = zero
    dsl4jsb_var3D_onChunk(SB_, transport_nh4_assoc)              = zero
    dsl4jsb_var3D_onChunk(SB_, transport_no3_solute)             = zero
    dsl4jsb_var3D_onChunk(SB_, transport_po4_solute)             = zero
    dsl4jsb_var3D_onChunk(SB_, transport_po4_assoc_fast)         = zero
    dsl4jsb_var3D_onChunk(SB_, transport_po4_assoc_slow)         = zero
    dsl4jsb_var3D_onChunk(SB_, transport_po4_occluded)           = zero
    dsl4jsb_var3D_onChunk(SB_, transport_po4_primary)            = zero
    dsl4jsb_var3D_onChunk(SB_, transport_nh4_n15_solute)         = zero
    dsl4jsb_var3D_onChunk(SB_, transport_nh4_n15_assoc)          = zero
    dsl4jsb_var3D_onChunk(SB_, transport_no3_n15_solute)         = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_nh4_solute_sl)       = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_no3_solute_sl)       = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_po4_solute_sl)       = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_carbon_sl)       = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_nitrogen_sl)     = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_phosphorus_sl)   = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_carbon13_sl)     = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_carbon14_sl)     = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_nitrogen15_sl)   = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_nh4_n15_solute_sl)   = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_no3_n15_solute_sl)   = zero
    dsl4jsb_var3D_onChunk(SB_, transport_noy)                    = zero
    dsl4jsb_var3D_onChunk(SB_, transport_n2o)                    = zero
    dsl4jsb_var3D_onChunk(SB_, transport_n2)                     = zero
    dsl4jsb_var3D_onChunk(SB_, transport_noy_n15)                = zero
    dsl4jsb_var3D_onChunk(SB_, transport_n2o_n15)                = zero
    dsl4jsb_var3D_onChunk(SB_, transport_n2_n15)                 = zero
    dsl4jsb_var3D_onChunk(SB_, volatilisation_nh4)               = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_no3)                = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_noy)                = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_n2o)                = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_noy)              = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_n2o)              = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_n2)               = zero
    dsl4jsb_var3D_onChunk(SB_, asymb_n_fixation)                 = zero
    dsl4jsb_var3D_onChunk(SB_, volatilisation_nh4_n15)           = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_no3_n15)            = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_noy_n15)            = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_n2o_n15)            = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_noy_n15)          = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_n2o_n15)          = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_n2_n15)           = zero
    dsl4jsb_var3D_onChunk(SB_, asymb_n_fixation_n15)             = zero

  END SUBROUTINE reset_sb_fluxes

  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for SB_
  !>
  SUBROUTINE calculate_time_average_soilbiogeochemistry(tile, options)

    USE mo_jsb_control,             ONLY: debug_on
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_lnd_time_averages        ! e.g. calc_time_mavg, mavg_period_tphen, mavg_period_weekly
    dsl4jsb_Use_processes SB_
    dsl4jsb_Use_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                !< instance of the model
    TYPE(t_lnd_bgcm_store), POINTER   :: bgcm_store           !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER   :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER   :: vgrid_soil_sb        !< Vertical grid
    INTEGER                           :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    INTEGER                           :: iblk, ics, ice, nc   !< dimensions
    CHARACTER(len=*), PARAMETER       :: routine = modname//':calculate_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! 3D SB_
    dsl4jsb_Real3D_onChunk      :: microbial_cue_eff
    dsl4jsb_Real3D_onChunk      :: microbial_nue_eff
    dsl4jsb_Real3D_onChunk      :: microbial_pue_eff
    dsl4jsb_Real3D_onChunk      :: dom_cn
    dsl4jsb_Real3D_onChunk      :: dom_cp
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_p
    dsl4jsb_Real3D_onChunk      :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: myc_export_c_tmyc_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_n_tmyc_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_n_tlabile_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_p_tlabile_mavg_sl
    dsl4jsb_Real3D_onChunk      :: dom_cn_mavg
    dsl4jsb_Real3D_onChunk      :: dom_cp_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_c_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_n_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_p_mavg
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN !< only if the present tile is a pft
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! 3D SB_
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_cue_eff)                  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_nue_eff)                  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_pue_eff)                  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cn)                             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cp)                             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_c)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_n)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_p)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_cue_eff_tmic_mavg)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_nue_eff_tmic_mavg)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_pue_eff_tmic_mavg)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_c_tmyc_mavg_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_n_tmyc_mavg_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_n_tlabile_mavg_sl)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_p_tlabile_mavg_sl)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cn_mavg)                        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cp_mavg)                        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_c_mavg)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_n_mavg)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_p_mavg)            ! out
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 daytime averages - not used here
    !>

    !>1.0 moving averages
    !>
    ! docu:
    ! calc_time_mavg(current average, new value, length of avg_period,  !
    !                do_calc=LOGICAL, avg_period_unit='day')            ! OPTIONAL
    !                RETURN(new current average)
    ! the unit of the averaging period is 'day' by default, but can also be 'week' or 'year'

    !>  1.1 tlabile (averages at the timescale of the labile pool)
    !>
    myc_export_n_tlabile_mavg_sl(:,:) = calc_time_mavg(myc_export_n_tlabile_mavg_sl(:,:), &
      &                                                sb_mycorrhiza_export_mt(ixN,:,:), mavg_period_tlabile)
    myc_export_p_tlabile_mavg_sl(:,:) = calc_time_mavg(myc_export_p_tlabile_mavg_sl(:,:), &
      &                                                sb_mycorrhiza_export_mt(ixP,:,:), mavg_period_tlabile)

    !>  1.2 lctlib tau_mycorrhiza
    !>
    myc_export_c_tmyc_mavg_sl(:,:) = calc_time_mavg(myc_export_c_tmyc_mavg_sl(:,:),                &
      &                                             sb_mycorrhiza_export_mt(ixC,:,:),        &
                                                    lctlib%tau_mycorrhiza, avg_period_unit='year')
    myc_export_n_tmyc_mavg_sl(:,:) = calc_time_mavg(myc_export_n_tmyc_mavg_sl(:,:),                &
      &                                             sb_mycorrhiza_export_mt(ixN,:,:),        &
                                                    lctlib%tau_mycorrhiza, avg_period_unit='year')

    !>  1.3 tmic (microbial community acclimation)
    !>
    microbial_cue_eff_tmic_mavg(:,:) = calc_time_mavg(microbial_cue_eff_tmic_mavg(:,:), microbial_cue_eff(:,:), &
                                                      mavg_period_tmic)
    microbial_nue_eff_tmic_mavg(:,:) = calc_time_mavg(microbial_nue_eff_tmic_mavg(:,:), microbial_nue_eff(:,:), &
                                                      mavg_period_tmic)
    microbial_pue_eff_tmic_mavg(:,:) = calc_time_mavg(microbial_pue_eff_tmic_mavg(:,:), microbial_pue_eff(:,:), &
                                                      mavg_period_tmic)

    !>  1.4 tenzyme (memory time-scale for enzyme allocation)
    !>
    dom_cn_mavg(:,:)             = calc_time_mavg(dom_cn_mavg(:,:), dom_cn(:,:), &
                                                  mavg_period_tenzyme)
    dom_cp_mavg(:,:)             = calc_time_mavg(dom_cp_mavg(:,:), dom_cp(:,:), &
                                                  mavg_period_tenzyme)
    enzyme_frac_poly_c_mavg(:,:) = calc_time_mavg(enzyme_frac_poly_c_mavg(:,:), enzyme_frac_poly_c(:,:), &
                                                  mavg_period_tenzyme)
    enzyme_frac_poly_n_mavg(:,:) = calc_time_mavg(enzyme_frac_poly_n_mavg(:,:), enzyme_frac_poly_n(:,:), &
                                                  mavg_period_tenzyme)
    enzyme_frac_poly_p_mavg(:,:) = calc_time_mavg(enzyme_frac_poly_p_mavg(:,:), enzyme_frac_poly_p(:,:), &
                                                  mavg_period_tenzyme)

  END SUBROUTINE calculate_time_average_soilbiogeochemistry

#endif
END MODULE mo_sb_util
