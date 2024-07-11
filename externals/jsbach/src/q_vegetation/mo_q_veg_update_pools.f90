!> QUINCY update vegetation pools
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
!>#### routines for calculating the vegetation pools using the fluxes that have been calculated
!>
MODULE mo_q_veg_update_pools
#ifndef __NO_QUINCY__

  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_TOTAL_BIO_ID, VEG_BGCM_GROWTH_ID, &
    &                                   VEG_BGCM_LITTERFALL_ID, VEG_BGCM_EXUDATION_ID, VEG_BGCM_ESTABLISHMENT_ID, &
    &                                   VEG_BGCM_RESERVE_USE_ID
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID, SB_BGCM_MYCO_EXPORT_ID

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_veg_pools

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_update_pools'

CONTAINS

  ! ======================================================================================================= !
  !>update vegetation pools
  !>
  SUBROUTINE update_veg_pools(tile, options)
    USE mo_kind,                            ONLY: wp
    USE mo_jsb_impl_constants,              ONLY: test_false_true
    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_tile_class,                  ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,                  ONLY: t_jsb_task_options
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element
    USE mo_jsb_model_class,                 ONLY: t_jsb_model
    USE mo_quincy_model_config,             ONLY: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
    USE mo_jsb_process_class,               ONLY: SPQ_, VEG_, Q_RAD_, Q_ASSIMI_, Q_PHENO_, SB_
    USE mo_jsb_grid_class,                  ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                        ONLY: Get_vgrid
    USE mo_jsb_physical_constants,          ONLY: lambda_C14, rhoh2o, grav
    USE mo_jsb_math_constants,              ONLY: dtime, one_day, one_year, eps4, eps8
    USE mo_isotope_util,                    ONLY: calc_mixing_ratio_N15N14, calc_fractionation, calc_mixing_ratio_C14C
    USE mo_q_rad_parameters,                ONLY: rfr_ratio_toc
    USE mo_q_veg_canopy,                    ONLY: calc_canopy_layers
    USE mo_q_veg_growth,                    ONLY: calc_diameter_from_woody_biomass, calc_height_from_diameter
    USE mo_veg_constants,                   ONLY: itree,eta_nfixation, max_leaf_shedding_rate, k_leafon_canopy, fresp_growth
    USE mo_q_assimi_parameters,             ONLY: CiCa_default_C3, CiCa_default_C4
    USE mo_q_assimi_constants,              ONLY: ic3phot, ic4phot
    USE mo_atmland_constants,               ONLY: def_co2_mixing_ratio, def_co2_mixing_ratio_C13, def_co2_mixing_ratio_C14, &
                                                  def_co2_deltaC13, def_co2_deltaC14
    USE mo_q_assimi_process,                ONLY: discrimination_ps
    USE mo_q_pheno_constants,               ONLY: ievergreen
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(Q_ASSIMI_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(Q_RAD_)
    dsl4jsb_Use_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER :: model                  !< the model
    TYPE(t_lnd_bgcm_store),   POINTER :: bgcm_store             !< the bgcm store of this tile
    TYPE(t_lctlib_element),   POINTER :: lctlib                 !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),        POINTER :: vgrid_canopy_q_assimi  !< Vertical grid
    TYPE(t_jsb_vgrid),        POINTER :: vgrid_soil_sb          !< Vertical grid
    INTEGER                           :: ncanopy                !< number of canopy layers
    INTEGER                           :: nsoil_sb               !< number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(options%nc)   :: hlp1                   !< dummy helper variable
    REAL(wp)                          :: lctlib_g1              !< set to g1_medlyn or g1_bberry depending on canopy_conductance_scheme
    INTEGER                           :: iblk, ics, ice, nc     !< grid dimensions
    INTEGER                           :: ic, is, ix_c, ix_e     !< looping indices
    INTEGER                           :: nr_of_compartments     !< number of compartments in veg bgcms
    INTEGER                           :: nr_of_elements         !< number of elements in veg bgcms
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_pools'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L2D :: veg_growth_mt
    dsl4jsb_Def_mt1L2D :: veg_exudation_mt
    dsl4jsb_Def_mt1L2D :: veg_establishment_mt
    dsl4jsb_Def_mt1L2D :: veg_reserve_use_mt
    dsl4jsb_Def_mt1L2D :: veg_total_biomass_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_ASSIMI_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(Q_RAD_)
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk      :: gross_assimilation
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C13
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C14
    dsl4jsb_Real2D_onChunk      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soa_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tfrac_mavg
    ! Q_PHENO_ 2D
    dsl4jsb_Real2D_onChunk      :: growing_season
    dsl4jsb_Real2D_onChunk      :: lai_max
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk      :: w_soil_root_pot
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl
    ! Q_RAD_ 2D
    dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc_tvegdyn_mavg
    ! Q_RAD_ 3D
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: do_cohort_harvest
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: delta_dens_ind
    dsl4jsb_Real2D_onChunk      :: diameter
    dsl4jsb_Real2D_onChunk      :: height
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: mean_leaf_age
    dsl4jsb_Real2D_onChunk      :: dphi
    dsl4jsb_Real2D_onChunk      :: npp
    dsl4jsb_Real2D_onChunk      :: growth_respiration
    dsl4jsb_Real2D_onChunk      :: maint_respiration
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration
    dsl4jsb_Real2D_onChunk      :: net_growth
    dsl4jsb_Real2D_onChunk      :: npp_c13
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c13
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c13
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c13
    dsl4jsb_Real2D_onChunk      :: npp_c14
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c14
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c14
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c14
    dsl4jsb_Real2D_onChunk      :: uptake_nh4
    dsl4jsb_Real2D_onChunk      :: uptake_no3
    dsl4jsb_Real2D_onChunk      :: n_fixation
    dsl4jsb_Real2D_onChunk      :: uptake_nh4_n15
    dsl4jsb_Real2D_onChunk      :: uptake_no3_n15
    dsl4jsb_Real2D_onChunk      :: n_fixation_n15
    dsl4jsb_Real2D_onChunk      :: uptake_po4
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_p
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n15
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_p
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_p
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n15
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n15
    dsl4jsb_Real2D_onChunk      :: t_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onChunk      :: n_transform_respiration
    dsl4jsb_Real2D_onChunk      :: net_biosphere_production_l2aveghlp
    dsl4jsb_Real2D_onChunk      :: biological_n_fixation_l2aveghlp
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_p
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_n15
    dsl4jsb_Real2D_onChunk      :: veg_pool_leaf_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_leaf_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_wood_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_wood_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_fine_root_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_fine_root_n
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_c
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_n
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_p
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_n15
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_c
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_n
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_p
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_n15
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: root_fraction_sl
    dsl4jsb_Real3D_onChunk      :: delta_root_fraction_sl
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: cumm_lai_cl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(tile%owner_model_id)
    lctlib                => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
    vgrid_soil_sb         => Get_vgrid('soil_layer_sb')
    ncanopy               =  vgrid_canopy_q_assimi%n_levels
    nsoil_sb              =  vgrid_soil_sb%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_ASSIMI_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(Q_RAD_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_GROWTH_ID, veg_growth_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_EXUDATION_ID, veg_exudation_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_ESTABLISHMENT_ID, veg_establishment_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_RESERVE_USE_ID, veg_reserve_use_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_TOTAL_BIO_ID, veg_total_biomass_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_ASSIMI_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation)        ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C13)    ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C14)    ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_air_tfrac_mavg)       ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soa_tphen_mavg)       ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_ps_tfrac_mavg)   ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_tfrac_mavg)   ! in
    ! Q_PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)             ! in
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, lai_max)                    ! in
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_soil_root_pot)                ! in
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_, soil_depth_sl)                  ! in
    ! Q_RAD_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, rfr_ratio_boc_tvegdyn_mavg)     ! in
    ! Q_RAD_ 3D
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_shaded_tfrac_mavg_cl)      ! in
    ! VEG 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, do_cohort_harvest)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)                             ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)                             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                               ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                                  ! CASE "plant" & "land": out / other Cases: in
    dsl4jsb_Get_var2D_onChunk(VEG_, mean_leaf_age)                        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi)                                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, npp)                                  ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth)                           ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_c13)                              ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c13)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c13)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c13)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_c14)                              ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c14)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c14)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c14)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4_n15)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3_n15)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_n15)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_po4)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_p)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n15)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_p)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_p)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n15)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n15)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tfrac_mavg)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tacclim_mavg)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tfrac_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tfrac_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tfrac_mavg)                        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tfrac_mavg)           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)                      ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_transform_respiration)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production_l2aveghlp)   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, biological_n_fixation_l2aveghlp)      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_n)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_p)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c13)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c14)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_n15)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_leaf_c)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_leaf_n)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_wood_c)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_wood_n)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_fine_root_c)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_fine_root_n)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_c)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_n)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_p)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_c13)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_c14)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_n15)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_c)               ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_n)               ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_p)               ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_c13)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_c14)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_n15)             ! out
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)     ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)               ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)                      ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)                       ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)                     ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)                      ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)                      ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)               ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, delta_root_fraction_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                         ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)                    ! out
    ! ----------------------------------------------------------------------------------------------------- !

    ! TODO: I am sure there is a more elegant way to acess these?
    nr_of_compartments = SIZE(veg_pool_mt,1)
    nr_of_elements = SIZE(veg_pool_mt,2)

    !>0.9 set g1 according to canopy_conductance_scheme
    !>
    SELECT CASE(TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme))
      CASE ("medlyn")
        lctlib_g1 = lctlib%g1_medlyn
      CASE ("ballberry")
        lctlib_g1 = lctlib%g1_bberry
    END SELECT

    !>1.0 for plant and land model, update vegetation pools with fluxes
    !>
    SELECT CASE(model%config%qmodel_id)
    CASE(QPLANT, QLAND)

      DO ic = 1, nc

       !>  1.1 net balance of the labile carbon pool due to photosynthesis and respiration
       !>
       npp(ic)        = gross_assimilation(ic)     - &
                        growth_respiration(ic)     - maint_respiration(ic)     - n_processing_respiration(ic)
       npp_c13(ic)    = gross_assimilation_C13(ic) - &
                        growth_respiration_c13(ic) - maint_respiration_c13(ic) - n_processing_respiration_c13(ic)
       npp_c14(ic)    = gross_assimilation_C14(ic) - &
                        growth_respiration_c14(ic) - maint_respiration_c14(ic) - n_processing_respiration_c14(ic)
       net_growth(ic) = net_growth(ic) + npp(ic)

       veg_growth_mt(ix_labile, ixC, ic)   = veg_growth_mt(ix_labile, ixC, ic)   + npp(ic)     * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixC13, ic) = veg_growth_mt(ix_labile, ixC13, ic) + npp_c13(ic) * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixC14, ic) = veg_growth_mt(ix_labile, ixC14, ic) + npp_c14(ic) * dtime / 1000000.0_wp

       !>  1.2 Add nutrients from uptake calculations
       !>
       n_fixation_n15(ic)                  = n_fixation(ic) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(-eta_nfixation))
       veg_growth_mt(ix_labile, ixN, ic)   = veg_growth_mt(ix_labile, ixN, ic) + &
                                             ( uptake_nh4(ic) + uptake_no3(ic) + n_fixation(ic)) * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixN15, ic) = veg_growth_mt(ix_labile, ixN15, ic) + &
                                             ( uptake_nh4_n15(ic) + uptake_no3_n15(ic) + n_fixation_n15(ic)) * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixP, ic)   = veg_growth_mt(ix_labile, ixP, ic) + &
                                             uptake_po4(ic) * dtime / 1000000.0_wp

       !>  1.3 Add mycorrhiza export
       !>
       DO is = 1, nsoil_sb
        veg_growth_mt(ix_labile, ixC, ic)   = veg_growth_mt(ix_labile, ixC, ic)   + &
          &                                   sb_mycorrhiza_export_mt(ixC, ic, is)    * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixC13, ic) = veg_growth_mt(ix_labile, ixC13, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixC13, ic, is)  * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixC14, ic) = veg_growth_mt(ix_labile, ixC14, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixC14, ic, is)  * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixN, ic)   = veg_growth_mt(ix_labile, ixN, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixN, ic, is)    * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixN15, ic) = veg_growth_mt(ix_labile, ixN15, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixN15, ic, is)  * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixP, ic)   = veg_growth_mt(ix_labile, ixP, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixP, ic, is)    * soil_depth_sl(ic,is)
       END DO

      END DO

       !>  1.4 update labile pool and reserve
       !>
       veg_pool_mt(ix_labile, :, :)  = veg_pool_mt(ix_labile, :, :)  + veg_reserve_use_mt(:,:)
       veg_pool_mt(ix_reserve, :, :) = veg_pool_mt(ix_reserve, :, :) - veg_reserve_use_mt(:,:)

      DO ic = 1, nc
        !>  1.5 update fluxes from continuous nutrient turnover and resorption
       !>
       ! N
       veg_pool_mt(ix_labile, ixN, ic)     = veg_pool_mt(ix_labile, ixN, ic) &
         &                                   + recycling_leaf_n(ic) + recycling_fine_root_n(ic) + recycling_heart_wood_n(ic)
       veg_pool_mt(ix_leaf, ixN, ic)       = veg_pool_mt(ix_leaf, ixN, ic)       - recycling_leaf_n(ic)
       veg_pool_mt(ix_fine_root, ixN, ic)  = veg_pool_mt(ix_fine_root, ixN, ic)  - recycling_fine_root_n(ic)
       veg_pool_mt(ix_heart_wood, ixN, ic) = veg_pool_mt(ix_heart_wood, ixN, ic) - recycling_heart_wood_n(ic)
       ! P
       veg_pool_mt(ix_labile, ixP, ic)     = veg_pool_mt(ix_labile, ixP, ic) &
         &                                   + recycling_leaf_p(ic) + recycling_fine_root_p(ic) + recycling_heart_wood_p(ic)
       veg_pool_mt(ix_leaf, ixP, ic)       = veg_pool_mt(ix_leaf, ixP, ic)       - recycling_leaf_p(ic)
       veg_pool_mt(ix_fine_root, ixP, ic)  = veg_pool_mt(ix_fine_root, ixP, ic)  - recycling_fine_root_p(ic)
       veg_pool_mt(ix_heart_wood, ixP, ic) = veg_pool_mt(ix_heart_wood, ixP, ic) - recycling_heart_wood_p(ic)
       ! N15
       veg_pool_mt(ix_labile, ixN15, ic)     = veg_pool_mt(ix_labile, ixN15, ic) &
         &                                     + recycling_leaf_n15(ic) + recycling_fine_root_n15(ic) + recycling_heart_wood_n15(ic)
       veg_pool_mt(ix_leaf, ixN15, ic)       = veg_pool_mt(ix_leaf, ixN15, ic)       - recycling_leaf_n15(ic)
       veg_pool_mt(ix_fine_root, ixN15, ic)  = veg_pool_mt(ix_fine_root, ixN15, ic)  - recycling_fine_root_n15(ic)
       veg_pool_mt(ix_heart_wood, ixN15, ic) = veg_pool_mt(ix_heart_wood, ixN15, ic) - recycling_heart_wood_n15(ic)

       !>  1.6 update labile pool with fluxes to mycorrhiza
       !>
       veg_pool_mt(ix_labile, ixC, ic)   = veg_pool_mt(ix_labile, ixC, ic)   - veg_exudation_mt(ixC, ic)
       veg_pool_mt(ix_labile, ixN, ic)   = veg_pool_mt(ix_labile, ixN, ic)   - veg_exudation_mt(ixN, ic)
       veg_pool_mt(ix_labile, ixP, ic)   = veg_pool_mt(ix_labile, ixP, ic)   - veg_exudation_mt(ixP, ic)
       veg_pool_mt(ix_labile, ixC13, ic) = veg_pool_mt(ix_labile, ixC13, ic) - veg_exudation_mt(ixC13, ic)
       veg_pool_mt(ix_labile, ixC14, ic) = veg_pool_mt(ix_labile, ixC14, ic) - veg_exudation_mt(ixC14, ic)
       veg_pool_mt(ix_labile, ixN15, ic) = veg_pool_mt(ix_labile, ixN15, ic) - veg_exudation_mt(ixN15, ic)
      END DO


       !>  1.7 Update plant biogeochemical pools with growth, litterfall and establishment
       !>
      DO ic = 1, nc
       DO ix_e = 1, nr_of_elements
        DO ix_c = 1, nr_of_compartments
          veg_pool_mt(ix_c,ix_e,ic) = veg_pool_mt(ix_c,ix_e,ic) - veg_litterfall_mt(ix_c,ix_e,ic)
          veg_pool_mt(ix_c,ix_e,ic) = veg_pool_mt(ix_c,ix_e,ic) +     veg_growth_mt(ix_c,ix_e,ic)
        END DO
       END DO
      END DO
      DO ic = 1, nc
       DO ix_e = 1, nr_of_elements
        veg_pool_mt(ix_labile, ix_e, ic)   = veg_pool_mt(ix_labile, ix_e, ic)   + veg_establishment_mt(ix_e,ic)
        veg_pool_mt(ix_seed_bed, ix_e, ic) = veg_pool_mt(ix_seed_bed, ix_e, ic) - veg_establishment_mt(ix_e,ic)
       END DO
      END DO

       ! @TODO
       !>    1.7.1 quick-fix ensure veg_pool does not become negative
       !>          apply only with IQ as it changes results with QS
       !>
       !   ==> need additional debug code
       !       to output litterfall & growth fluxes + lat/lon/PFT for grid cells where veg_pool gets negative
#ifdef __QUINCY_STANDALONE__
#else
      DO ic = 1, nc
       DO ix_e = 1, nr_of_elements
        DO ix_c = 1, nr_of_compartments
         veg_pool_mt(ix_c,ix_e,ic) = MAX(0.0_wp, veg_pool_mt(ix_c,ix_e,ic))
        END DO
       END DO
      END DO
#endif

       !>  1.8 Radioactive decay of C14 in living plant tissue
       !>
      DO ic = 1, nc
       DO ix_c = 1, nr_of_compartments
        veg_pool_mt(ix_c, ixC14, ic) = veg_pool_mt(ix_c, ixC14, ic) * (1._wp - lambda_C14 * dtime )
       END DO
      END DO

    CASE (QCANOPY)

      DO ic = 1, nc

      ! diagnose required change in leaf C given phenology to attain prescribed LAI max
       IF (growing_season(ic) > test_false_true) THEN
         veg_growth_mt(ix_leaf, ixC, ic) = lai_max(ic) / lctlib%sla &
           &                              * (1._wp - exp(-k_leafon_canopy &
           &                              * (lai_max(ic) - lai(ic)) / lai_max(ic))) &
           &                              * dtime / one_day
       ELSE
         veg_litterfall_mt(ix_leaf, ixC, ic) = MIN(veg_pool_mt(ix_leaf, ixC, ic), &
           &                                  lai(ic)/lctlib%sla * max_leaf_shedding_rate * dtime / one_day &
           &                                  * lai_max(ic) &
           &                                  / MAX(eps4, lai(ic)))
       END IF

       ! calc leaf turnover for evergreen trees (needed, for e.g. mean_leaf_age, as evergreens are always within growing_season)
       IF (lctlib%phenology_type == ievergreen) THEN
         veg_litterfall_mt(ix_leaf, ixC, ic) = veg_pool_mt(ix_leaf, ixC, ic) / lctlib%tau_leaf / one_day / one_year * dtime
       ENDIF

       ! update leaf C:N:P pools
       veg_pool_mt(ix_leaf, ixC, ic) = veg_pool_mt(ix_leaf, ixC, ic) + veg_growth_mt(ix_leaf, ixC, ic) &
         &                            - veg_litterfall_mt(ix_leaf, ixC, ic)
       veg_pool_mt(ix_leaf, ixN, ic) = veg_pool_mt(ix_leaf, ixC, ic) / lctlib%cn_leaf
       veg_pool_mt(ix_leaf, ixP, ic) = veg_pool_mt(ix_leaf, ixN, ic) / lctlib%np_leaf

      END DO

    END SELECT

    !>2.0 Update plant diagnostics given the above fluxes
    !>
    lai(:)      = veg_pool_mt(ix_leaf, ixC, :) * lctlib%sla
    dens_ind(:) = dens_ind(:) + delta_dens_ind(:)
    DO ic = 1, nc
      IF (dens_ind(ic) > eps8 .AND. lctlib%growthform == itree) THEN
        diameter(ic)  = calc_diameter_from_woody_biomass(lctlib%allom_k1               , &
                                                         lctlib%allom_k2               , &
                                                         lctlib%wood_density           , &
                                                         veg_pool_mt(ix_sap_wood, ixC, ic)   , &
                                                         veg_pool_mt(ix_heart_wood, ixC, ic) , &
                                                         dens_ind(ic))
        height(ic)    = calc_height_from_diameter(lctlib%allom_k1, lctlib%allom_k2, diameter(ic))
      END IF

      IF (veg_pool_mt(ix_leaf, ixC, ic) > eps8) THEN
        mean_leaf_age(ic) = (mean_leaf_age(ic) + dtime / one_day) &
          &                * (1._wp - veg_growth_mt(ix_leaf, ixC, ic) / veg_pool_mt(ix_leaf, ixC, ic)) &
          &                - mean_leaf_age(ic) * veg_litterfall_mt(ix_leaf, ixC, ic) / veg_pool_mt(ix_leaf, ixC, ic)
      ELSE
        mean_leaf_age(ic) = 0.0_wp
      END IF

      dphi(ic) = w_soil_root_pot(ic) - lctlib%phi_leaf_min - grav * rhoh2o * height(ic) * 1.e-6_wp
    END DO

    !>  2.1 implied change of root fraction given root growth
    !>
    IF (dsl4jsb_Config(VEG_)%flag_dynamic_roots) THEN
      DO is = 1, nsoil_sb
        DO ic = 1, nc
          root_fraction_sl(ic,is) = root_fraction_sl(ic,is) + delta_root_fraction_sl(ic,is)
        END DO
      END DO
    END IF

    !>3.0 stand/cohort harvest
    !>
    ! in case of complete harvest, set all pools to zero
    ! (to account for wood extraction, which would otherwise leave some wood, and reset the key plant diagnostics)
    IF(model%config%l_do_stand_replacing_harvest .OR. do_cohort_harvest(ics) > test_false_true) THEN
      ! differ between cohort & population scheme
      IF (do_cohort_harvest(ics) > test_false_true) THEN    ! -> cohort scheme
        dens_ind(:)                     = 1.0_wp
        veg_pool_mt(ix_reserve, ixC, :) = 1.0_wp
      ELSE ! model%config%l_do_stand_replacing_harvest -> population scheme
        dens_ind(:)                     = 0.1_wp
        veg_pool_mt(ix_reserve, ixC, :) = 0.1_wp
      ENDIF
      IF(dsl4jsb_Config(SB_)%flag_mycorrhiza) THEN !save some additional C in case plant has to feed mycorrhizae
        veg_pool_mt(ix_reserve, ixC, :) = veg_pool_mt(ix_reserve, ixC, :) * 2.0_wp
      ENDIF
      ! set zero
      veg_pool_mt(ix_leaf, :, :)         = 0.0_wp
      veg_pool_mt(ix_fine_root, :, :)    = 0.0_wp
      veg_pool_mt(ix_coarse_root, :, :)  = 0.0_wp
      veg_pool_mt(ix_sap_wood, :, :)     = 0.0_wp
      veg_pool_mt(ix_heart_wood, :, :)   = 0.0_wp
      veg_pool_mt(ix_fruit, :, :)        = 0.0_wp
      veg_pool_mt(ix_labile, :, :)       = 0.0_wp
      height(:)                          = 0.0_wp
      diameter(:)                        = 0.0_wp
      lai(:)                             = 0.0_wp
      mean_leaf_age(:)                   = 0.0_wp
      ! assume that some (10% ?) mycorrhizae will survive the harvest - most will die shortly after harvest caused by turnover
      sb_pool_mt(ix_mycorrhiza, :, :, :) = sb_pool_mt(ix_mycorrhiza, :, :, :) * 0.1_wp
      ! init
      dphi(:)                       = w_soil_root_pot(:) - lctlib%phi_leaf_min
      rfr_ratio_boc_tvegdyn_mavg(:) = rfr_ratio_toc
      ! pool%reserve% C N P
      veg_pool_mt(ix_reserve, ixN, :) = veg_pool_mt(ix_reserve, ixC, :) / (1._wp + fresp_growth) / lctlib%cn_leaf
      veg_pool_mt(ix_reserve, ixP, :) = veg_pool_mt(ix_reserve, ixN, :) / lctlib%np_leaf
      ! pool%reserve%carbon13
      IF(lctlib%ps_pathway == ic3phot)THEN
         hlp1(:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
      ELSE
         hlp1(:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
      ENDIF
      veg_pool_mt(ix_reserve, ixC13, :) = calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13, hlp1(:)) &
        &                                 * veg_pool_mt(ix_reserve, ixC, :)
      ! pool%reserve%carbon14
      IF(lctlib%ps_pathway == ic3phot)THEN
         hlp1(:) = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
      ELSE
         hlp1(:) = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
      ENDIF
      veg_pool_mt(ix_reserve, ixC14, :) = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14-hlp1(:)) &
        &                                 * veg_pool_mt(ix_reserve, ixC, :)
      ! pool%reserve%nitrogen15
      veg_pool_mt(ix_reserve, ixN15, :) = veg_pool_mt(ix_reserve, ixN, :) / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    ENDIF

    !>4.0 Given current vegetation pools, update canopy layers
    !>
    CALL calc_canopy_layers( nc                                                       , & ! in
                             ncanopy                                                  , &
                             vgrid_canopy_q_assimi%dz(:)                              , &
                             vgrid_canopy_q_assimi%lbounds(:)                         , &
                             vgrid_canopy_q_assimi%ubounds(:)                         , &
                             lctlib%ps_pathway                                        , &
                             lctlib%k0_fn_struc                                       , &
                             lctlib%fn_oth_min                                        , &
                             lctlib%sla                                               , &
                             lctlib%np_leaf                                           , &
                             lctlib%gmin                                              , &
                             lctlib%g0                                                , &
                             lctlib_g1                                                , &
                             lctlib%t_jmax_omega                                      , &
                             dsl4jsb_Config(Q_ASSIMI_)%flag_optimal_Nfraction         , &
                             TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme), & !   Q_ASSIMI_ config (medlyn/ballberry)
                             veg_pool_mt(ix_leaf, ixN, :)                             , &
                             lai(:)                                                   , &
                             ppfd_sunlit_tfrac_mavg_cl(:,:)                           , &
                             ppfd_shaded_tfrac_mavg_cl(:,:)                           , &
                             fleaf_sunlit_tfrac_mavg_cl(:,:)                          , &
                             fn_rub_cl(:,:)                                           , &
                             fn_et_cl(:,:)                                            , &
                             fn_pepc_cl(:,:)                                          , &
                             fn_chl_cl(:,:)                                           , &
                             fn_oth_cl(:,:)                                           , & ! in
                             lai_cl(:,:)                                              , & ! out
                             cumm_lai_cl(:,:)                                         , & ! out
                             leaf_nitrogen_cl(:,:)                                    , & ! out
                             t_air            = t_air_tfrac_mavg(:), &                    ! optional in
                             t_acclim         = t_air_tacclim_mavg(:), &
                             press_srf        = press_srf_tfrac_mavg(:), &
                             co2_mixing_ratio = co2_mixing_ratio_tfrac_mavg(:), &
                             ga               = ga_tfrac_mavg(:), &
                             beta_air         = beta_air_tfrac_mavg(:), &
                             beta_soa         = beta_soa_tphen_mavg(:), &
                             beta_soil_ps     = beta_soil_ps_tfrac_mavg(:), &
                             beta_sinklim_ps  = beta_sinklim_ps_tfrac_mavg(:), &
                             beta_soil_gs     = beta_soil_gs_tfrac_mavg(:), &
                             t_jmax_opt       = t_jmax_opt_mavg(:) )                      ! optional in

    !> 5.0 analysis of mass conservation: sum compartments of veg_pool into veg_total_biomass
    !>
    veg_total_biomass_mt(:,:) = SUM(veg_pool_mt(:, :, :), DIM = 1)

    !> 6.0 bgc_material diagnostics
    !>     sum up veg_pool, veg_growth and veg_litterfall across components for each particular element
    !>     or simply add value (of one/two/few pool variables) to a memory variable for aggregation and output
    !>
    ! pools
    veg_pool_total_c(:)         = SUM(veg_pool_mt(:, ixC, :), DIM = 1)
    veg_pool_total_n(:)         = SUM(veg_pool_mt(:, ixN, :), DIM = 1)
    veg_pool_total_p(:)         = SUM(veg_pool_mt(:, ixP, :), DIM = 1)
    veg_pool_total_c13(:)       = SUM(veg_pool_mt(:, ixC13, :), DIM = 1)
    veg_pool_total_c14(:)       = SUM(veg_pool_mt(:, ixC14, :), DIM = 1)
    veg_pool_total_n15(:)       = SUM(veg_pool_mt(:, ixN15, :), DIM = 1)
    veg_pool_leaf_c(:)          = veg_pool_mt(ix_leaf, ixC, :)
    veg_pool_leaf_n(:)          = veg_pool_mt(ix_leaf, ixN, :)
    veg_pool_wood_c(:)          = veg_pool_mt(ix_sap_wood, ixC, :) + veg_pool_mt(ix_heart_wood, ixC, :)
    veg_pool_wood_n(:)          = veg_pool_mt(ix_sap_wood, ixN, :) + veg_pool_mt(ix_heart_wood, ixN, :)
    veg_pool_fine_root_c(:)     = veg_pool_mt(ix_fine_root, ixC, :)
    veg_pool_fine_root_n(:)     = veg_pool_mt(ix_fine_root, ixN, :)
    ! fluxes
    veg_growth_total_c(:)       = SUM(veg_growth_mt(:, ixC, :), DIM = 1)
    veg_growth_total_n(:)       = SUM(veg_growth_mt(:, ixN, :), DIM = 1)
    veg_growth_total_p(:)       = SUM(veg_growth_mt(:, ixP, :), DIM = 1)
    veg_growth_total_c13(:)     = SUM(veg_growth_mt(:, ixC13, :), DIM = 1)
    veg_growth_total_c14(:)     = SUM(veg_growth_mt(:, ixC14, :), DIM = 1)
    veg_growth_total_n15(:)     = SUM(veg_growth_mt(:, ixN15, :), DIM = 1)
    veg_litterfall_total_c(:)   = SUM(veg_litterfall_mt(:, ixC, :), DIM = 1)
    veg_litterfall_total_n(:)   = SUM(veg_litterfall_mt(:, ixN, :), DIM = 1)
    veg_litterfall_total_p(:)   = SUM(veg_litterfall_mt(:, ixP, :), DIM = 1)
    veg_litterfall_total_c13(:) = SUM(veg_litterfall_mt(:, ixC13, :), DIM = 1)
    veg_litterfall_total_c14(:) = SUM(veg_litterfall_mt(:, ixC14, :), DIM = 1)
    veg_litterfall_total_n15(:) = SUM(veg_litterfall_mt(:, ixN15, :), DIM = 1)

    !> 6.0 biosphere-level diagnostics across multiple processes
    !>
    !>   see also 'update_sb_pools()'
    !>
    !>   TODO calculation may include 'fFire' and 'fLUC' (land-use change) once available
    !>
    ! could also be calc as: 'NPP - het_resp - fFire - fLUC'
    ! (n_processing_respiration = n_fixation_respiration + n_transform_respiration)
    net_biosphere_production_l2aveghlp(:) = net_biosphere_production_l2aveghlp(:) &
      & + gross_assimilation(:) - maint_respiration(:) - growth_respiration(:) - n_processing_respiration(:)
    ! biological N fixation (VEG_ and SB_ processes)
    biological_n_fixation_l2aveghlp(:)    = biological_n_fixation_l2aveghlp(:) + n_fixation(:)

  END SUBROUTINE update_veg_pools

#endif
END MODULE mo_q_veg_update_pools
