!> QUINCY soil-biogeochemistry variables init
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
!>#### initialization of soil-biogeochemistry memory variables using, e.g., ic & bc input files
!>
MODULE mo_sb_init
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_jsb_math_constants,    ONLY: eps4, eps8

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store, get_bgcm_idx
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID

  USE mo_jsb_process_class,       ONLY: SB_, SPQ_, VEG_

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sb_init

  TYPE t_sb_init_vars
    REAL(wp), POINTER :: &
      & usda_taxonomy_class     (:,:) => NULL(), &
      & nwrb_taxonomy_class     (:,:) => NULL(), &
      & qmax_org_fine_particle  (:,:) => NULL(), &
      & soil_ph                 (:,:) => NULL(), &
      & soil_p_labile           (:,:) => NULL(), &
      & soil_p_slow             (:,:) => NULL(), &
      & soil_p_occluded         (:,:) => NULL(), &
      & soil_p_primary          (:,:) => NULL()
    INTEGER, POINTER :: &
      & usda_taxonomy_class_int (:,:) => NULL(), &
      & nwrb_taxonomy_class_int (:,:) => NULL()
  END TYPE t_sb_init_vars

  TYPE(t_sb_init_vars) :: sb_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_init'

CONTAINS

  ! ======================================================================================================= !
  !> Intialize SB_ process
  !>
  SUBROUTINE sb_init(tile)

    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_init'

#ifdef __QUINCY_STANDALONE__
    CALL sb_qs_read_init_vars(tile)
    CALL sb_init_ic_bc(tile)
    CALL sb_finalize_init_vars()
#else
    ! call only at root tile, i.e., the only tile without associated parent tile (to avoid unnessary i/o)
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL sb_read_init_vars(tile)
    END IF

    CALL sb_init_ic_bc(tile)

    IF (tile%Is_last_process_tile(SB_)) THEN
      CALL sb_finalize_init_vars()
    END IF
#endif
  END SUBROUTINE sb_init

  ! ======================================================================================================= !
  !> Intialize SB_ process from ic and bc input files
  !>
  SUBROUTINE sb_init_ic_bc(tile)
    USE mo_jsb_class,                     ONLY: Get_model
    USE mo_jsb_tile_class,                ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,               ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,              ONLY: t_lctlib_element
    USE mo_jsb_grid_class,                ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,                      ONLY: Get_grid, Get_vgrid
    USE mo_q_assimi_constants,            ONLY: IC3PHOT
    USE mo_q_assimi_process,              ONLY: discrimination_ps
    USE mo_q_assimi_parameters,           ONLY: CiCa_default_C3, CiCa_default_C4
    USE mo_atmland_constants,             ONLY: def_co2_mixing_ratio, def_co2_mixing_ratio_C13, def_co2_mixing_ratio_C14, &
      &                                         def_co2_deltaC13, def_co2_deltaC14
    USE mo_sb_constants,                  ONLY: microbial_cue_max, microbial_nue, microbial_pue, microbial_cn, microbial_np, &
      &                                         k_weath_mineral
    USE mo_isotope_util,                  ONLY: calc_fractionation, calc_mixing_ratio_C14C
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER :: model                 !< the model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store            !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER :: lctlib                !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_grid),       POINTER :: hgrid                 !< Horizontal grid
    TYPE(t_jsb_vgrid),      POINTER :: vgrid_soil_sb         !< Vertical grid
    INTEGER                         :: nsoil_sb              !< number of soil layers as used/defined by the SB_ process
    REAL(wp)                        :: hlp1                  !< helper var
    REAL(wp)                        :: hlp_frac_c13          !< helper var
    REAL(wp)                        :: hlp_frac_c14          !< helper var
    CHARACTER(len=3)                :: site_ID_spp1685       !< QUINCY standalone, SPP site-set, sideID
    INTEGER                         :: iblk, ic, is          !< loop over dimensions
    INTEGER                         :: nblks                 !< number of blocks
    INTEGER                         :: nproma                !< number of grid points per iblk
    INTEGER                         :: sb_pool_idx           !< index of the sb pool within the bgcm store
    REAL(wp),           POINTER :: sb_pool_mt_domain(:,:,:,:,:) !< dim: compartments, elements, nc, nsoil, nblks
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init_ic_bc'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! VEG_ 3D
    dsl4jsb_Real3D_onDomain :: root_fraction_sl
    ! SPQ_ 2D
    dsl4jsb_Real2D_onDomain :: num_sl_above_bedrock
    ! SPQ_ 3D
    dsl4jsb_Real3D_onDomain :: bulk_dens_sl
    dsl4jsb_Real3D_onDomain :: w_soil_sl
    dsl4jsb_Real3D_onDomain :: soil_depth_sl             ! (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Real3D_onDomain :: soil_lay_width_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onDomain :: qmax_org_min_sl
    dsl4jsb_Real3D_onDomain :: qmax_po4_min_sl
    dsl4jsb_Real3D_onDomain :: qmax_nh4_min_sl
    dsl4jsb_Real3D_onDomain :: qmax_po4_om_sl
    dsl4jsb_Real3D_onDomain :: silt_sl
    dsl4jsb_Real3D_onDomain :: clay_sl
    dsl4jsb_Real3D_onDomain :: volume_min_sl
    ! SB_ 3D
    dsl4jsb_Real3D_onDomain :: nh4_assoc
    dsl4jsb_Real3D_onDomain :: nh4_solute
    dsl4jsb_Real3D_onDomain :: no3_solute
    dsl4jsb_Real3D_onDomain :: po4_solute
    dsl4jsb_Real3D_onDomain :: nh4_n15_assoc
    dsl4jsb_Real3D_onDomain :: nh4_n15_solute
    dsl4jsb_Real3D_onDomain :: no3_n15_solute
    dsl4jsb_Real3D_onDomain :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onDomain :: soil_litter_carbon_sl
    dsl4jsb_Real3D_onDomain :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onDomain :: qmax_org
    dsl4jsb_Real3D_onDomain :: qmax_po4
    dsl4jsb_Real3D_onDomain :: qmax_nh4
    dsl4jsb_Real3D_onDomain :: qmax_fast_po4
    dsl4jsb_Real3D_onDomain :: qmax_slow_po4
    dsl4jsb_Real3D_onDomain :: km_fast_po4
    dsl4jsb_Real3D_onDomain :: km_slow_po4
    dsl4jsb_Real3D_onDomain :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onDomain :: km_adsorpt_nh4_sl
    dsl4jsb_Real3D_onDomain :: k_bioturb
    dsl4jsb_Real3D_onDomain :: po4_assoc_fast
    dsl4jsb_Real3D_onDomain :: po4_assoc_slow
    dsl4jsb_Real3D_onDomain :: po4_occluded
    dsl4jsb_Real3D_onDomain :: ph_sl
    dsl4jsb_Real3D_onDomain :: Qmax_AlFe_cor
    dsl4jsb_Real3D_onDomain :: po4_primary
    dsl4jsb_Real3D_onDomain :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: vmax_weath_mineral_sl
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly
    dsl4jsb_Real3D_onDomain :: enzyme_frac_residue
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_p
    dsl4jsb_Real3D_onDomain :: dom_cn
    dsl4jsb_Real3D_onDomain :: dom_cp
    dsl4jsb_Real3D_onDomain :: fact_n_status_mic_c_growth
    dsl4jsb_Real3D_onDomain :: fact_p_status_mic_c_growth
    dsl4jsb_Real3D_onDomain :: enzyme_frac_AP
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN                !< run this init only if the present tile is a pft
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    hgrid         => Get_grid(model%grid_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    nblks         =  hgrid%nblks
    nproma        =  hgrid%nproma
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    sb_pool_idx = get_bgcm_idx(bgcm_store, SB_BGCM_POOL_ID, tile%name, routine)
    sb_pool_mt_domain => bgcm_store%store_2l_3d_bgcms(bgcm_store%idx_in_store(sb_pool_idx))%mt_2l_3d_bgcm(:,:,:,:,:)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onDomain(VEG_,   root_fraction_sl)
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onDomain(SPQ_,   num_sl_above_bedrock)     ! in
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onDomain(SPQ_,   bulk_dens_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   w_soil_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   soil_depth_sl)            ! in  (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   soil_lay_width_sl)        ! in
    dsl4jsb_Get_var3D_onDomain(SPQ_,   soil_lay_depth_center_sl) ! in
    dsl4jsb_Get_var3D_onDomain(SPQ_,   soil_lay_depth_ubound_sl) ! in
    dsl4jsb_Get_var3D_onDomain(SPQ_,   qmax_org_min_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   qmax_po4_min_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   qmax_nh4_min_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   qmax_po4_om_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   silt_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   clay_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_,   volume_min_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_assoc)
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    no3_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_n15_assoc)
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_n15_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    no3_n15_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    bulk_soil_carbon_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    soil_litter_carbon_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    bulk_dens_corr_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_org)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_nh4)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_fast_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_slow_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_fast_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_slow_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_adsorpt_po4_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_adsorpt_nh4_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    k_bioturb)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_assoc_fast)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_assoc_slow)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_occluded)
    dsl4jsb_Get_var3D_onDomain(SB_,    ph_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    Qmax_AlFe_cor)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_primary)
    dsl4jsb_Get_var3D_onDomain(SB_,    microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_,    microbial_nue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_,    microbial_pue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_,    vmax_weath_mineral_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_residue)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly_c)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly_n)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly_p)
    dsl4jsb_Get_var3D_onDomain(SB_,    dom_cn)
    dsl4jsb_Get_var3D_onDomain(SB_,    dom_cp)
    dsl4jsb_Get_var3D_onDomain(SB_,    fact_n_status_mic_c_growth)
    dsl4jsb_Get_var3D_onDomain(SB_,    fact_p_status_mic_c_growth)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_AP)
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 init bgcm with zero
    !>  NOTE: reading the restart files is done after init
    !>  NOTE: not using a "Get" function in process init for bgcm
    !>
    sb_pool_mt_domain(:,:,:,:,:) = 0.0_wp

    !>1.0 all soil models
    !>  simple_sm & jsm
    !>

    !>  1.1 calc hlp_frac_c13 & hlp_frac_c14
    !>
    ! C13
    IF (lctlib%ps_pathway == IC3PHOT) THEN
      hlp1 = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, &
        &                      def_co2_mixing_ratio, &
        &                      'C13', &
        &                      lctlib%ps_pathway)
    ELSE
      hlp1 = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, &
        &                      def_co2_mixing_ratio, &
        &                      'C13', &
        &                      lctlib%ps_pathway)
    ENDIF
    hlp_frac_c13 = calc_fractionation(def_co2_mixing_ratio, &
      &                               def_co2_mixing_ratio_C13, &
      &                               hlp1 )
    ! C14
    IF (lctlib%ps_pathway == IC3PHOT) THEN
      hlp1 = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, &
        &                      def_co2_mixing_ratio, &
        &                      'C14', &
        &                      lctlib%ps_pathway)
    ELSE
      hlp1 = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, &
        &                      def_co2_mixing_ratio, &
        &                      'C14', &
        &                      lctlib%ps_pathway)
    ENDIF
    hlp_frac_c14 = calc_mixing_ratio_C14C(def_co2_deltaC13, &
      &                                   def_co2_deltaC14 - hlp1)

    !>  1.2 init site-specific soil properties
    !>
    ph_sl(:,:,:)              = SPREAD(sb_init_vars%soil_ph(:,:), DIM = 2, ncopies = nsoil_sb)
    Qmax_AlFe_cor(:,:,:)      = 1._wp                           !< ... docu ...

    !>  1.3 init SB_ variables from constants
    !>
    microbial_cue_eff_tmic_mavg(:,:,:)  = microbial_cue_max
    microbial_nue_eff_tmic_mavg(:,:,:)  = microbial_nue
    microbial_pue_eff_tmic_mavg(:,:,:)  = microbial_pue
    vmax_weath_mineral_sl(:,:,:)        = k_weath_mineral
    dom_cn(:,:,:)                       = microbial_cn
    dom_cp(:,:,:)                       = microbial_cn * microbial_np
    fact_n_status_mic_c_growth(:,:,:)   = 1._wp                           !< ... docu ...
    fact_p_status_mic_c_growth(:,:,:)   = 1._wp                           !< ... docu ...

#ifdef __QUINCY_STANDALONE__
    !------------------------------------------------------------------------------------------------------ !
    !>  1.4 SPP site-set: identify site based on longitude
    !>  five sites from the "SPP 1685 Project: Forest Strategies for limited Phosphorus Resources"
    !>
    site_ID_spp1685 = "OTH"   ! default: OTH = other
    IF (model%config%flag_spp1685) THEN
      IF (ABS(model%config%lon - 9.75_wp) < eps8) THEN
        site_ID_spp1685 = "BBR"
      ELSEIF (ABS(model%config%lon - 10.25_wp) < eps8) THEN
        site_ID_spp1685 = "LUE"
      ELSEIF (ABS(model%config%lon - 10.75_wp) < eps8) THEN
        site_ID_spp1685 = "VES"
      ELSEIF (ABS(model%config%lon - 12.75_wp) < eps8) THEN
        site_ID_spp1685 = "MIT"
      ELSEIF (ABS(model%config%lon - 7.75_wp)  < eps8) THEN
        site_ID_spp1685 = "CON"
      ENDIF
    ENDIF
    !------------------------------------------------------------------------------------------------------ !
#endif

    SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
    !>2.0 ssm - a simple four pool soil biogeochemical model
    !>
    CASE("simple_1d")
      CALL sb_init_simple_sm( &
        & nproma, &                                   ! in
        & nsoil_sb, &
        & nblks, &
        & num_sl_above_bedrock(:,:), &
        & lctlib%growthform, &
        & lctlib%k_som_fast_init, &
        & lctlib%k_som_slow_init, &
        & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
        & dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
        & SPREAD(sb_init_vars%qmax_org_fine_particle(:,:), DIM = 2, ncopies = nsoil_sb), &
        & SPREAD(sb_init_vars%soil_p_labile(:,:), DIM = 2, ncopies = nsoil_sb), &
        & SPREAD(sb_init_vars%soil_p_slow(:,:), DIM = 2, ncopies = nsoil_sb), &
        & SPREAD(sb_init_vars%soil_p_occluded(:,:), DIM = 2, ncopies = nsoil_sb), &
        & SPREAD(sb_init_vars%soil_p_primary(:,:), DIM = 2, ncopies = nsoil_sb), &
        & dsl4jsb_Config(SB_)%soil_p_depth, &
        & soil_depth_sl(:,:,:), &
        & soil_lay_width_sl(:,:,:), &
        & soil_lay_depth_ubound_sl(:,:,:), &
        & soil_lay_depth_center_sl(:,:,:), &
        & hlp_frac_c13, &
        & hlp_frac_c14, &
        & root_fraction_sl(:,:,:), &
        & clay_sl(:,:,:), &
        & silt_sl(:,:,:), &
        & w_soil_sl(:,:,:), &
        & bulk_dens_sl(:,:,:), &                      ! in
        & ph_sl(:,:,:), &                             ! inout
        & bulk_soil_carbon_sl(:,:,:), &
        & soil_litter_carbon_sl(:,:,:), &
        & po4_solute(:,:,:), &
        & volume_min_sl(:,:,:), &
        & qmax_po4(:,:,:), &
        & qmax_po4_min_sl(:,:,:), &
        & qmax_po4_om_sl(:,:,:), &
        & qmax_fast_po4(:,:,:), &
        & qmax_slow_po4(:,:,:), &
        & km_fast_po4(:,:,:), &
        & km_slow_po4(:,:,:), &
        & km_adsorpt_po4_sl(:,:,:), &
        & po4_assoc_fast(:,:,:), &
        & po4_assoc_slow(:,:,:), &
        & po4_occluded(:,:,:), &
        & po4_primary(:,:,:), &
        & Qmax_AlFe_cor(:,:,:), &
        & sb_pool_mt_domain(:,:,:,:,:), &                    ! inout
        & nh4_solute(:,:,:), &                        ! out
        & no3_solute(:,:,:), &
        & nh4_n15_solute(:,:,:), &
        & no3_n15_solute(:,:,:), &
        & bulk_dens_corr_sl(:,:,:), &
        & qmax_org(:,:,:), &
        & qmax_org_min_sl(:,:,:), &
        & qmax_nh4(:,:,:), &
        & qmax_nh4_min_sl(:,:,:), &
        & km_adsorpt_nh4_sl(:,:,:), &
        & nh4_assoc(:,:,:), &
        & nh4_n15_assoc(:,:,:) &                      ! out
#ifdef __QUINCY_STANDALONE__
        & , &
        & model%config%flag_spp1685, &                ! optional in
        & site_ID_spp1685 &                           ! optional in
#endif
        & )
    !>3.0 jsm - Jena Soil Model (described in https://doi.org/10.5194/gmd-13-783-2020)
    !>
    CASE("jsm")
      CALL sb_init_jsm( &
        & nproma, &                                   ! in
        & nsoil_sb, &
        & nblks, &
        & num_sl_above_bedrock(:,:), &
        & lctlib%growthform, &
        & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
        & dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
        & SPREAD(sb_init_vars%qmax_org_fine_particle(:,:), DIM = 2, ncopies = nsoil_sb), &
        & SPREAD(sb_init_vars%soil_p_labile(:,:), DIM = 2, ncopies = nsoil_sb), &
        & SPREAD(sb_init_vars%soil_p_slow(:,:), DIM = 2, ncopies = nsoil_sb), &
        & SPREAD(sb_init_vars%soil_p_occluded(:,:), DIM = 2, ncopies = nsoil_sb), &
        & SPREAD(sb_init_vars%soil_p_primary(:,:), DIM = 2, ncopies = nsoil_sb), &
        & dsl4jsb_Config(SB_)%soil_p_depth, &
        & soil_depth_sl(:,:,:), &
        & soil_lay_width_sl(:,:,:), &
        & soil_lay_depth_ubound_sl(:,:,:), &
        & soil_lay_depth_center_sl(:,:,:), &
        & hlp_frac_c13, &
        & hlp_frac_c14, &
        & root_fraction_sl(:,:,:), &
        & clay_sl(:,:,:), &
        & silt_sl(:,:,:), &
        & w_soil_sl(:,:,:), &
        & bulk_dens_sl(:,:,:), &                      ! in
        & ph_sl(:,:,:), &                             ! inout
        & qmax_po4_min_sl(:,:,:), &
        & qmax_po4_om_sl(:,:,:), &
        & volume_min_sl(:,:,:), &
        & po4_solute(:,:,:), &
        & bulk_soil_carbon_sl(:,:,:), &
        & soil_litter_carbon_sl(:,:,:), &
        & qmax_po4(:,:,:), &
        & qmax_nh4(:,:,:), &
        & qmax_fast_po4(:,:,:), &
        & qmax_slow_po4(:,:,:), &
        & km_fast_po4(:,:,:), &
        & km_slow_po4(:,:,:), &
        & km_adsorpt_po4_sl(:,:,:), &
        & po4_assoc_fast(:,:,:), &
        & po4_assoc_slow(:,:,:), &
        & po4_occluded(:,:,:), &
        & Qmax_AlFe_cor(:,:,:), &
        & po4_primary(:,:,:), &
        & sb_pool_mt_domain(:,:,:,:,:), &                    ! inout
        & qmax_org_min_sl(:,:,:), &                   ! out
        & qmax_nh4_min_sl(:,:,:), &
        & nh4_assoc(:,:,:), &
        & nh4_solute(:,:,:), &
        & no3_solute(:,:,:), &
        & nh4_n15_assoc(:,:,:), &
        & nh4_n15_solute(:,:,:), &
        & no3_n15_solute(:,:,:), &
        & bulk_dens_corr_sl(:,:,:), &
        & qmax_org(:,:,:), &
        & km_adsorpt_nh4_sl(:,:,:), &
        & k_bioturb(:,:,:), &
        & enzyme_frac_poly(:,:,:), &
        & enzyme_frac_residue(:,:,:), &
        & enzyme_frac_poly_c(:,:,:), &
        & enzyme_frac_poly_n(:,:,:), &
        & enzyme_frac_poly_p(:,:,:), &
        & enzyme_frac_AP(:,:,:) &                     ! out
#ifdef __QUINCY_STANDALONE__
        & , &
        & model%config%flag_spp1685, &                ! optional in
        & site_ID_spp1685 &                           ! optional in
#endif
        & )
    CASE DEFAULT
      CALL finish(routine, 'Invalid setting of sb_model_scheme (simple_sm / jsm).')
    END SELECT  ! simple_sm | jsm

    !> 3.0 set sb_pool bgcm to zero for tiles such as bare soil and urban area
    !>
    IF (lctlib%BareSoilFlag) THEN
      sb_pool_mt_domain(:,:,:,:,:) = 0.0_wp
    END IF

  END SUBROUTINE sb_init_ic_bc

  ! ======================================================================================================= !
  !> de-allocate SPQ_ init vars
  !>
  SUBROUTINE sb_finalize_init_vars
    DEALLOCATE( &
      & sb_init_vars%usda_taxonomy_class      , &
      & sb_init_vars%nwrb_taxonomy_class      , &
      & sb_init_vars%qmax_org_fine_particle   , &
      & sb_init_vars%soil_ph                  , &
      & sb_init_vars%soil_p_labile            , &
      & sb_init_vars%soil_p_slow              , &
      & sb_init_vars%soil_p_occluded          , &
      & sb_init_vars%soil_p_primary           , &
      & sb_init_vars%usda_taxonomy_class_int  , &
      & sb_init_vars%nwrb_taxonomy_class_int)
  END SUBROUTINE sb_finalize_init_vars

#ifdef __QUINCY_STANDALONE__
  SUBROUTINE sb_qs_read_init_vars(tile)
    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
    USE mo_jsb_io,             ONLY: missval

    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SB_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    INTEGER                     :: nsoil_sb
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_qs_read_init_vars'

    IF (debug_on()) CALL message(routine, 'Reading/setting SB_ init vars')

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels

    ALLOCATE( &
      & sb_init_vars%usda_taxonomy_class     (nproma, nblks),       &
      & sb_init_vars%nwrb_taxonomy_class     (nproma, nblks),       &
      & sb_init_vars%qmax_org_fine_particle  (nproma, nblks),       &
      & sb_init_vars%soil_ph                 (nproma, nblks),       &
      & sb_init_vars%soil_p_labile           (nproma, nblks),       &
      & sb_init_vars%soil_p_slow             (nproma, nblks),       &
      & sb_init_vars%soil_p_occluded         (nproma, nblks),       &
      & sb_init_vars%soil_p_primary          (nproma, nblks),       &
      & sb_init_vars%usda_taxonomy_class_int (nproma, nblks),       &
      & sb_init_vars%nwrb_taxonomy_class_int (nproma, nblks)        &
      & )

    sb_init_vars%usda_taxonomy_class     (:,:) = missval
    sb_init_vars%nwrb_taxonomy_class     (:,:) = missval
    sb_init_vars%qmax_org_fine_particle  (:,:) = missval
    sb_init_vars%soil_ph                 (:,:) = missval
    sb_init_vars%soil_p_labile           (:,:) = missval
    sb_init_vars%soil_p_slow             (:,:) = missval
    sb_init_vars%soil_p_occluded         (:,:) = missval
    sb_init_vars%soil_p_primary          (:,:) = missval
    sb_init_vars%usda_taxonomy_class_int (:,:) = -9999
    sb_init_vars%nwrb_taxonomy_class_int (:,:) = -9999

    dsl4jsb_Get_config(SB_)
    ! ----------------------------------------------------------------------------------------------------- !

    !> USDA taxonomy soil class
    !>
    sb_init_vars%usda_taxonomy_class_int(:,:) = dsl4jsb_Config(SB_)%usda_taxonomy_class

    !> NWRB taxonomy soil class
    !>
    sb_init_vars%nwrb_taxonomy_class_int(:,:) = dsl4jsb_Config(SB_)%nwrb_taxonomy_class

    !> qmax organic fine particle
    !>
    sb_init_vars%qmax_org_fine_particle(:,:) = dsl4jsb_Config(SB_)%qmax_org_fine_particle

    !> soil ph
    !>
    sb_init_vars%soil_ph(:,:) = dsl4jsb_Config(SB_)%soil_ph

    !> P labile
    !>
    sb_init_vars%soil_p_labile(:,:) = dsl4jsb_Config(SB_)%soil_p_labile

    !> P slow
    !>
    sb_init_vars%soil_p_slow(:,:) = dsl4jsb_Config(SB_)%soil_p_slow

    !> P occluded
    !>
    sb_init_vars%soil_p_occluded(:,:) = dsl4jsb_Config(SB_)%soil_p_occluded

    !> P primary (apatite)
    !>
    sb_init_vars%soil_p_primary(:,:) = dsl4jsb_Config(SB_)%soil_p_primary
  END SUBROUTINE sb_qs_read_init_vars
#else
  ! ======================================================================================================= !
  !> read SB_ init vars from input file
  !>
  SUBROUTINE sb_read_init_vars(tile)
    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_jsb_io,             ONLY: missval

    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SB_)

    REAL(wp), POINTER           :: ptr_2D(:, :)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb
    TYPE(t_input_file)          :: input_file
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    INTEGER                     :: nsoil_sb
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_read_init_vars'

    IF (debug_on()) CALL message(routine, 'Reading/setting SB_ init vars')

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels

    ALLOCATE( &
      & sb_init_vars%usda_taxonomy_class     (nproma, nblks),       &
      & sb_init_vars%nwrb_taxonomy_class     (nproma, nblks),       &
      & sb_init_vars%qmax_org_fine_particle  (nproma, nblks),       &
      & sb_init_vars%soil_ph                 (nproma, nblks),       &
      & sb_init_vars%soil_p_labile           (nproma, nblks),       &
      & sb_init_vars%soil_p_slow             (nproma, nblks),       &
      & sb_init_vars%soil_p_occluded         (nproma, nblks),       &
      & sb_init_vars%soil_p_primary          (nproma, nblks),       &
      & sb_init_vars%usda_taxonomy_class_int (nproma, nblks),       &
      & sb_init_vars%nwrb_taxonomy_class_int (nproma, nblks)        &
      & )

    sb_init_vars%usda_taxonomy_class     (:,:) = missval
    sb_init_vars%nwrb_taxonomy_class     (:,:) = missval
    sb_init_vars%qmax_org_fine_particle  (:,:) = missval
    sb_init_vars%soil_ph                 (:,:) = missval
    sb_init_vars%soil_p_labile           (:,:) = missval
    sb_init_vars%soil_p_slow             (:,:) = missval
    sb_init_vars%soil_p_occluded         (:,:) = missval
    sb_init_vars%soil_p_primary          (:,:) = missval
    sb_init_vars%usda_taxonomy_class_int (:,:) = -9999
    sb_init_vars%nwrb_taxonomy_class_int (:,:) = -9999

    dsl4jsb_Get_config(SB_)
    ! ----------------------------------------------------------------------------------------------------- !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SB_)%bc_quincy_soil_filename), model%grid_id)


    !> USDA taxonomy soil class
    !>
    ptr_2D => input_file%Read_2d(       &
      & variable_name='TAXOUSDA_10km',  &
      & fill_array = sb_init_vars%usda_taxonomy_class)
    sb_init_vars%usda_taxonomy_class(:,:)     = MERGE(ptr_2D, 1.0_wp, ptr_2D > 1.0_wp)
    sb_init_vars%usda_taxonomy_class_int(:,:) = INT(sb_init_vars%usda_taxonomy_class(:,:))

    !> NWRB taxonomy soil class
    !>
    ptr_2D => input_file%Read_2d(  &
      & variable_name='TAXNWRB',   &
      & fill_array = sb_init_vars%nwrb_taxonomy_class)
    sb_init_vars%nwrb_taxonomy_class(:,:)     = MERGE(ptr_2D, 1.0_wp, ptr_2D > 1.0_wp)
    sb_init_vars%nwrb_taxonomy_class_int(:,:) = INT(sb_init_vars%nwrb_taxonomy_class(:,:))

    !> qmax organic fine particle
    !>
    ptr_2D => input_file%Read_2d(  &
      & variable_name='qmax_org',  &
      & fill_array = sb_init_vars%qmax_org_fine_particle)
    sb_init_vars%qmax_org_fine_particle(:,:) = MERGE(ptr_2D, 6.5537_wp, ptr_2D > 0.0_wp)

    !> soil ph
    !>
    ptr_2D => input_file%Read_2d(   &
      & variable_name='PHIHOX',     &
      & fill_array = sb_init_vars%soil_ph)
    sb_init_vars%soil_ph(:,:) = MERGE(ptr_2D, 65.0_wp, ptr_2D > 0.0_wp)
    ! divide by 10 to get the actual values within the correct range
    sb_init_vars%soil_ph(:,:) = sb_init_vars%soil_ph(:,:) / 10.0_wp

    !> P labile
    !>
    ptr_2D => input_file%Read_2d(           &
      & variable_name='Labile_Inorganic_P', &
      & fill_array = sb_init_vars%soil_p_labile)
    sb_init_vars%soil_p_labile(:,:) = MERGE(ptr_2D, 50._wp, ptr_2D > 0.0_wp)

    !> P slow
    !>
    ptr_2D => input_file%Read_2d(             &
      & variable_name='Seconday_Mineral__P',  &
      & fill_array = sb_init_vars%soil_p_slow)
    sb_init_vars%soil_p_slow(:,:) = MERGE(ptr_2D, 70._wp, ptr_2D > 0.0_wp)

    !> P occluded
    !>
    ptr_2D => input_file%Read_2d(     &
      & variable_name='Occluded_P',   &
      & fill_array = sb_init_vars%soil_p_occluded)
    sb_init_vars%soil_p_occluded(:,:) = MERGE(ptr_2D, 150._wp, ptr_2D > 0.0_wp)

    !> P primary (apatite)
    !>
    ptr_2D => input_file%Read_2d(  &
      & variable_name='Apatite_P',  &
      & fill_array = sb_init_vars%soil_p_primary)
    sb_init_vars%soil_p_primary(:,:) = MERGE(ptr_2D, 100._wp, ptr_2D > 0.0_wp)

    CALL input_file%Close()
  END SUBROUTINE sb_read_init_vars
#endif

  ! ======================================================================================================= !
  !>init SB_ with the simple soil model
  !>
  SUBROUTINE sb_init_simple_sm( &
    & nproma, &
    & nsoil_sb, &
    & nblks, &
    & num_sl_above_bedrock, &
    & lctlib_growthform, &
    & lctlib_k_som_fast_init, &
    & lctlib_k_som_slow_init, &
    & sb_model_scheme, &
    & flag_sb_double_langmuir, &
    & qmax_org_fine_particle, &
    & soil_p_labile, &
    & soil_p_slow, &
    & soil_p_occluded, &
    & soil_p_primary, &
    & soil_p_depth, &
    & soil_depth_sl, &
    & soil_lay_width_sl, &
    & soil_lay_depth_ubound_sl, &
    & soil_lay_depth_center_sl, &
    & hlp_frac_c13, &
    & hlp_frac_c14, &
    & root_fraction_sl, &
    & clay_sl, &
    & silt_sl, &
    & w_soil_sl, &
    & bulk_dens_sl, &
    & ph_sl, &
    & bulk_soil_carbon_sl, &
    & soil_litter_carbon_sl, &
    & po4_solute, &
    & volume_min_sl, &
    & qmax_po4, &
    & qmax_po4_min_sl, &
    & qmax_po4_om_sl, &
    & qmax_fast_po4, &
    & qmax_slow_po4, &
    & km_fast_po4, &
    & km_slow_po4, &
    & km_adsorpt_po4_sl, &
    & po4_assoc_fast, &
    & po4_assoc_slow, &
    & po4_occluded, &
    & po4_primary, &
    & Qmax_AlFe_cor, &
    & sb_pool_mt_domain, &
    & nh4_solute, &
    & no3_solute, &
    & nh4_n15_solute, &
    & no3_n15_solute, &
    & bulk_dens_corr_sl, &
    & qmax_org, &
    & qmax_org_min_sl, &
    & qmax_nh4, &
    & qmax_nh4_min_sl, &
    & km_adsorpt_nh4_sl, &
    & nh4_assoc, &
    & nh4_n15_assoc, &
    & flag_spp1685, &   ! optional in
    & site_ID_spp1685 ) ! optional in

    USE mo_q_sb_jsm_processes,            ONLY: calc_bulk_soil_carbon, calc_bulk_density_correction, &
      &                                         calc_qmax_bulk_density_correction, calc_Psorption_parameter, &
      &                                         calc_fast_po4
    USE mo_spq_util,                      ONLY: calc_qmax_texture
    USE mo_isotope_util,                  ONLY: calc_mixing_ratio_N15N14
    USE mo_veg_constants,                 ONLY: ITREE
    USE mo_sb_constants,                  ONLY: def_sb_pool_metabol_litter, def_sb_pool_simple_struct_litter, &
      &                                         def_sb_pool_woody_litter, def_sb_pool_fast_som, def_sb_pool_slow_som, &
      &                                         nh4_solute_prescribe, no3_solute_prescribe, po4_solute_prescribe, &
      &                                         p_labile_slow_global_avg, frac_prod_not_root, qmax_nh4_clay, &
      &                                         km_adsorpt_OM_nh4, km_adsorpt_mineral_nh4, &
      &                                         k_init_soluable_litter_cn, k_init_litter_np, k_init_polymeric_litter_cn, &
      &                                         k_init_woody_litter_cn, k_fast_som_np, k_fast_som_cn_min, k_fast_som_cn_max, &
      &                                         k_slow_som_cn, k_slow_som_np, &
      &                                         reference_depth_som_init, k_som_init, &
      &                                         sb_pool_total_som_at_ref_depth_init, fract_som_fast
    USE mo_jsb_physical_constants,        ONLY: molar_mass_C
#ifdef __QUINCY_STANDALONE__
    USE mo_qs_process_init_util,          ONLY: init_sb_soil_properties_spp1685_sites, &
      &                                         init_sb_soil_p_pools_spp1685_sites
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                    INTENT(in)    :: nproma
    INTEGER,                    INTENT(in)    :: nsoil_sb
    INTEGER,                    INTENT(in)    :: nblks
    REAL(wp),                   INTENT(in)    :: num_sl_above_bedrock(:,:)
    INTEGER,                    INTENT(in)    :: lctlib_growthform
    REAL(wp),                   INTENT(in)    :: lctlib_k_som_fast_init
    REAL(wp),                   INTENT(in)    :: lctlib_k_som_slow_init
    CHARACTER(len=*),           INTENT(in)    :: sb_model_scheme                  !< SB_ config
    LOGICAL,                    INTENT(in)    :: flag_sb_double_langmuir          !< SB_ config
    REAL(wp),                   INTENT(in)    :: qmax_org_fine_particle(:,:,:)    !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_labile(:,:,:)             !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_slow(:,:,:)               !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_occluded(:,:,:)           !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_primary(:,:,:)            !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_depth                     !< SB_ config
    REAL(wp),                   INTENT(in)    :: soil_depth_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_width_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_depth_ubound_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_depth_center_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: hlp_frac_c13
    REAL(wp),                   INTENT(in)    :: hlp_frac_c14
    REAL(wp),                   INTENT(in)    :: root_fraction_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: clay_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: silt_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: w_soil_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: bulk_dens_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: ph_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: bulk_soil_carbon_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: soil_litter_carbon_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_solute(:,:,:)
    REAL(wp),                   INTENT(inout) :: volume_min_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4_min_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4_om_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_fast_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_slow_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_fast_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_slow_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_adsorpt_po4_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_assoc_fast(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_assoc_slow(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_occluded(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_primary(:,:,:)
    REAL(wp),                   INTENT(inout) :: Qmax_AlFe_cor(:,:,:)
    REAL(wp),                   INTENT(inout) :: sb_pool_mt_domain(:,:,:,:,:)          !< bgc_material sb_pool
    REAL(wp),                   INTENT(out)   :: nh4_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: no3_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_n15_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: no3_n15_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: bulk_dens_corr_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_org(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_org_min_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_nh4(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_nh4_min_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: km_adsorpt_nh4_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_assoc(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_n15_assoc(:,:,:)
    LOGICAL,          OPTIONAL, INTENT(in)    :: flag_spp1685
    CHARACTER(len=3), OPTIONAL, INTENT(in)    :: site_ID_spp1685
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), ALLOCATABLE       :: hlp1(:,:)
    REAL(wp), ALLOCATABLE       :: hlp2(:,:)
    REAL(wp), ALLOCATABLE       :: weight_fast_sl(:,:,:)
    REAL(wp), ALLOCATABLE       :: weight_slow_sl(:,:,:)
    REAL(wp), ALLOCATABLE       :: init_soluable_litter(:,:)
    REAL(wp), ALLOCATABLE       :: init_polymeric_litter(:,:)
    REAL(wp), ALLOCATABLE       :: init_fast(:,:)
    REAL(wp), ALLOCATABLE       :: init_slow(:,:)
    REAL(wp), ALLOCATABLE       :: p_labile(:,:,:)
    REAL(wp), ALLOCATABLE       :: p_labile_slow_site(:,:,:)
    INTEGER                     :: iblk, is, ic
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init_simple_sm'
    IF (debug_on()) CALL message(TRIM(routine), 'Starting ...')
    ! ----------------------------------------------------------------------------------------------------- !

    !> init out var
    !>
    nh4_solute(:,:,:)         = 0.0_wp
    no3_solute(:,:,:)         = 0.0_wp
    nh4_n15_solute(:,:,:)     = 0.0_wp
    no3_n15_solute(:,:,:)     = 0.0_wp
    bulk_dens_corr_sl(:,:,:)  = 0.0_wp
    qmax_org(:,:,:)           = 0.0_wp
    qmax_org_min_sl(:,:,:)    = 0.0_wp
    qmax_nh4(:,:,:)           = 0.0_wp
    qmax_nh4_min_sl(:,:,:)    = 0.0_wp
    km_adsorpt_nh4_sl(:,:,:)  = 0.0_wp
    nh4_assoc(:,:,:)          = 0.0_wp
    nh4_n15_assoc(:,:,:)      = 0.0_wp

    !>0.9 allocate local var
    !>
    ALLOCATE(hlp1(nproma, nblks))
    ALLOCATE(hlp2(nproma, nblks))
    ALLOCATE(weight_fast_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(weight_slow_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(init_soluable_litter(nproma, nblks))
    ALLOCATE(init_polymeric_litter(nproma, nblks))
    ALLOCATE(init_fast(nproma, nblks))
    ALLOCATE(init_slow(nproma, nblks))
    ALLOCATE(p_labile(nproma, nsoil_sb, nblks))
    ALLOCATE(p_labile_slow_site(nproma, nsoil_sb, nblks))

    !>1.0 init
    !>
    hlp1(:,:)             = 0.0_wp
    hlp2(:,:)             = 0.0_wp
    weight_fast_sl(:,:,:) = 0.0_wp
    weight_slow_sl(:,:,:) = 0.0_wp
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,NINT(num_sl_above_bedrock(ic, iblk))
          weight_fast_sl(ic, is, iblk) = EXP(-(lctlib_k_som_fast_init * soil_lay_depth_center_sl(ic, is, iblk)) ** k_som_init)
          weight_slow_sl(ic, is, iblk) = EXP(-(lctlib_k_som_slow_init * soil_lay_depth_center_sl(ic, is, iblk)) ** k_som_init)
          ! only count up to reference depth for scaling
          IF (soil_lay_depth_center_sl(ic, is, iblk) < reference_depth_som_init) THEN
            hlp1(ic, iblk) = hlp1(ic, iblk) + weight_fast_sl(ic, is, iblk) * soil_lay_width_sl(ic, is, iblk)
            hlp2(ic, iblk) = hlp2(ic, iblk) + weight_slow_sl(ic, is, iblk) * soil_lay_width_sl(ic, is, iblk)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! calculate top-soil SOM concentration and convert from default gC/m2 to mol C / m3
    WHERE(hlp1 > eps8)
      init_soluable_litter(:,:)   = def_sb_pool_metabol_litter * fract_som_fast / molar_mass_C / hlp1(:,:)
      init_polymeric_litter(:,:)  = def_sb_pool_simple_struct_litter * fract_som_fast / molar_mass_C / hlp1(:,:)
      init_fast(:,:)              = sb_pool_total_som_at_ref_depth_init * fract_som_fast / molar_mass_C / hlp1(:,:)
    END WHERE
    WHERE(hlp2 > eps8)
      init_slow(:,:) = sb_pool_total_som_at_ref_depth_init * (1._wp - fract_som_fast) / molar_mass_C / hlp2(:,:)
    END WHERE

    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,NINT(num_sl_above_bedrock(ic, iblk))
          sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk)  = init_soluable_litter(ic, iblk) &
            &                                                         * weight_fast_sl(ic, is, iblk)
          sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = init_polymeric_litter(ic, iblk) &
            &                                                         * weight_fast_sl(ic, is, iblk)
          sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk)        = init_fast(ic, iblk) * weight_fast_sl(ic, is, iblk)
          sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk)          = init_slow(ic, iblk) * weight_slow_sl(ic, is, iblk)
        ENDDO
      ENDDO
    ENDDO
    ! re-calculate the first layer of the wood litter pool for trees only
    IF (lctlib_growthform == ITREE) THEN
      sb_pool_mt_domain(ix_woody_litter, ixC, :, 1, :) = &
        &   sb_pool_total_som_at_ref_depth_init / molar_mass_C / soil_lay_width_sl(:,1,:)
    ENDIF

    ! NOTE: this differs on purpose from the calculation in mo_q_sb_update_pools for the prescribed option
    !       this is because here I'm trying to ensure that the models are initialised with the same amount of nutrient
    !       whereas in the prescribed option the target is unlimiting concentrations across the profile. The latter leads
    !       to more nutrient in the soil column in 1d compared to 0d.
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,NINT(num_sl_above_bedrock(ic, iblk))
          nh4_solute(ic, is, iblk) = nh4_solute_prescribe &
            &                           / soil_lay_depth_ubound_sl(ic, NINT(num_sl_above_bedrock(ic, iblk)), iblk)
          no3_solute(ic, is, iblk) = no3_solute_prescribe &
            &                           / soil_lay_depth_ubound_sl(ic, NINT(num_sl_above_bedrock(ic, iblk)), iblk)
        ENDDO
      ENDDO
    ENDDO

    ! inital bulk density and qmax_po4 based on current soil carbon content
    DO iblk = 1,nblks
      CALL calc_bulk_soil_carbon( &
        & nproma, &                           ! in
        & nsoil_sb, &
        & num_sl_above_bedrock(:,iblk), &
        & TRIM(sb_model_scheme), &
        & sb_pool_mt_domain(:,:,:,:,iblk), &  ! in
        & bulk_soil_carbon_sl(:,:,iblk), &    ! inout
        & soil_litter_carbon_sl(:,:,iblk), &
        & volume_min_sl(:,:,iblk) )           ! inout

      bulk_dens_corr_sl(:,:,iblk) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:,iblk), &
        &                                                        soil_litter_carbon_sl(:,:,iblk), &
        &                                                        bulk_dens_sl(:,:,iblk))
      qmax_org_min_sl(:,:,iblk) = calc_qmax_texture(qmax_org_fine_particle(:,:,iblk), &
        &                                           silt_sl(:,:,iblk), &
        &                                           clay_sl(:,:,iblk))

      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &  ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & 0._wp, &
        & qmax_org_min_sl(:,:,iblk), &      ! in
        & qmax_org(:,:,iblk) )              ! out
    END DO

#ifdef __QUINCY_STANDALONE__
    ! specific init values for sites of the SPP site-set
    ! for: ph_sl(:,:,:) and Qmax_AlFe_cor(:,:,:)
    IF (flag_spp1685) THEN
      DO iblk = 1,nblks
        CALL init_sb_soil_properties_spp1685_sites( &
          & nproma, &                   ! in
          & nsoil_sb, &
          & site_ID_spp1685, &          ! in
          & ph_sl(:,:,iblk), &          ! out
          & Qmax_AlFe_cor(:,:,iblk) )   ! out
      END DO
    ENDIF
#endif

    ! Correct the 1st soil layer Qmax_AlFe_cor(:,:,:) with the woody litter
    Qmax_AlFe_cor(:, 1, :) = Qmax_AlFe_cor(:, 1, :) * 1._wp / ((1._wp - volume_min_sl(:, 1, :)) &
      &                      * (bulk_soil_carbon_sl(:, 1, :) + soil_litter_carbon_sl(:, 1, :)  &
      &                      + sb_pool_mt_domain(ix_woody_litter, ixC, :,  1, :)) &
      &                      / (bulk_soil_carbon_sl(:, 1, :) + soil_litter_carbon_sl(:, 1, :)) &
      &                      + volume_min_sl(:, 1, :))
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
          CALL calc_Psorption_parameter( &
            & flag_sb_double_langmuir, &          ! in
            & clay_sl(ic, is, iblk), &
            & silt_sl(ic, is, iblk), &
            & w_soil_sl(ic, is, iblk) / soil_depth_sl(ic, is, iblk) * 1000._wp, &
            & ph_sl(ic, is, iblk), &
            & bulk_soil_carbon_sl(ic, is, iblk), &
            & soil_litter_carbon_sl(ic, is, iblk), &
            & sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk), &
            & bulk_dens_sl(ic, is, iblk), &
            & volume_min_sl(ic, is, iblk), &
            & po4_solute_prescribe, &             ! in - po4_solute(ic, is, iblk)
            & qmax_po4_min_sl(ic, is, iblk), &        ! inout
            & qmax_po4_om_sl(ic, is, iblk), &
            & qmax_po4(ic, is, iblk), &
            & km_adsorpt_po4_sl(ic, is, iblk), &
            & qmax_fast_po4(ic, is, iblk), &
            & qmax_slow_po4(ic, is, iblk), &
            & km_fast_po4(ic, is, iblk), &
            & km_slow_po4(ic, is, iblk), &
            & Qmax_AlFe_cor(ic, is, iblk) )           ! inout
        ENDDO
      ENDDO
    ENDDO
    ! calculate the partition coefficient between labile and slow P pool based on
    ! their Qmax (double-Langmuir) or global average values (single-Langmuir)
    IF (flag_sb_double_langmuir) THEN
      p_labile_slow_site(:,:,:) = eps4
      ! if 'qmax_po4(:,:,:) > eps8'
      DO iblk = 1,nblks
        DO is = 1,nsoil_sb
          DO ic = 1,nproma
            IF (qmax_po4(ic, is, iblk) > eps8) THEN
              p_labile_slow_site(ic, is, iblk) = MAX(qmax_fast_po4(ic, is, iblk) / qmax_po4(ic, is, iblk), eps4)
            END IF
          END DO
        END DO
      END DO
    ELSE
      p_labile_slow_site(:,:,:) = p_labile_slow_global_avg
    ENDIF

    ! nh4 adsorption parameter intialisation
    qmax_nh4_min_sl(:,:,:) = calc_qmax_texture(qmax_nh4_clay, &
      &                                        0.0_wp, &
      &                                        clay_sl(:,:,:))
    DO iblk = 1,nblks
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &  ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & 0.0_wp, &
        & qmax_nh4_min_sl(:,:,iblk), &      ! in
        & qmax_nh4(:,:,iblk) )              ! out
    ENDDO
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
          qmax_nh4(ic, is, iblk) = qmax_nh4(ic, is, iblk) * w_soil_sl(ic, is, iblk) / soil_depth_sl(ic, is, iblk)
        ENDDO
      ENDDO
    ENDDO
    DO iblk = 1,nblks
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &  ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & km_adsorpt_OM_nh4, &
        & km_adsorpt_mineral_nh4, &         ! in
        & km_adsorpt_nh4_sl(:,:,iblk) )     ! out
    ENDDO
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
          nh4_assoc(ic, is, iblk)     = qmax_nh4(ic, is, iblk) * nh4_solute(ic, is, iblk) &
              &                           / (nh4_solute(ic, is, iblk) + km_adsorpt_nh4_sl(ic, is, iblk))
          nh4_n15_assoc(ic, is, iblk) = nh4_assoc(ic, is, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! metabolic litter all other elements
          sb_pool_mt_domain(ix_soluable_litter, ixN, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) / k_init_soluable_litter_cn
          sb_pool_mt_domain(ix_soluable_litter, ixP, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_soluable_litter, ixN, ic, is, iblk) / k_init_litter_np
          sb_pool_mt_domain(ix_soluable_litter, ixC13, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) * hlp_frac_c13
          sb_pool_mt_domain(ix_soluable_litter, ixC14, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) * hlp_frac_c14
          sb_pool_mt_domain(ix_soluable_litter, ixN15, ic, is, iblk) = sb_pool_mt_domain(ix_soluable_litter, ixN, ic, is, iblk) &
            &                                                          / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! non-woody litter all other elements
          sb_pool_mt_domain(ix_polymeric_litter, ixN, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) / k_init_polymeric_litter_cn
          sb_pool_mt_domain(ix_polymeric_litter, ixP, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_polymeric_litter, ixN, ic, is, iblk) / k_init_litter_np
          sb_pool_mt_domain(ix_polymeric_litter, ixC13, ic, is, iblk) = sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
            &                                                           * hlp_frac_c13
          sb_pool_mt_domain(ix_polymeric_litter, ixC14, ic, is, iblk) = sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
            &                                                           * hlp_frac_c14
          sb_pool_mt_domain(ix_polymeric_litter, ixN15, ic, is, iblk) = sb_pool_mt_domain(ix_polymeric_litter, ixN, ic, is, iblk) &
            &                                                           / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! woody litter all other elements
          sb_pool_mt_domain(ix_woody_litter, ixN, ic, is, iblk) &
            &  = sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk) / k_init_woody_litter_cn
          sb_pool_mt_domain(ix_woody_litter, ixP, ic, is, iblk)   = sb_pool_mt_domain(ix_woody_litter, ixN, ic, is, iblk) &
            &                                                       / k_init_litter_np
          sb_pool_mt_domain(ix_woody_litter, ixC13, ic, is, iblk) = sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk) &
            &                                                       * hlp_frac_c13
          sb_pool_mt_domain(ix_woody_litter, ixC14, ic, is, iblk) = sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk) &
            &                                                       * hlp_frac_c14
          sb_pool_mt_domain(ix_woody_litter, ixN15, ic, is, iblk) = sb_pool_mt_domain(ix_woody_litter, ixN, ic, is, iblk) &
            &                                                       / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! fast pool (= microbial) all other elements
          sb_pool_mt_domain(ix_microbial, ixN, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) / (k_fast_som_cn_min + k_fast_som_cn_max) * 2.0_wp
          sb_pool_mt_domain(ix_microbial, ixP, ic, is, iblk)   = sb_pool_mt_domain(ix_microbial, ixN, ic, is, iblk) / k_fast_som_np
          sb_pool_mt_domain(ix_microbial, ixC13, ic, is, iblk) = sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) * hlp_frac_c13
          sb_pool_mt_domain(ix_microbial, ixC14, ic, is, iblk) = sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) * hlp_frac_c14
          sb_pool_mt_domain(ix_microbial, ixN15, ic, is, iblk) = sb_pool_mt_domain(ix_microbial, ixN, ic, is, iblk) &
            &                                                    / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! ..
          sb_pool_mt_domain(ix_residue, ixN, ic, is, iblk)   = sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) / k_slow_som_cn
          sb_pool_mt_domain(ix_residue, ixP, ic, is, iblk)   = sb_pool_mt_domain(ix_residue, ixN, ic, is, iblk) / k_slow_som_np
          sb_pool_mt_domain(ix_residue, ixC13, ic, is, iblk) = sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) * hlp_frac_c13
          sb_pool_mt_domain(ix_residue, ixC14, ic, is, iblk) = sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) * hlp_frac_c14
          sb_pool_mt_domain(ix_residue, ixN15, ic, is, iblk) = sb_pool_mt_domain(ix_residue, ixN, ic, is, iblk) &
            &                                                  / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! ..
          nh4_n15_solute(ic, is, iblk) = nh4_solute(ic, is, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          no3_n15_solute(ic, is, iblk) = no3_solute(ic, is, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        ENDDO
      ENDDO
    ENDDO

    !> Soil P pools intialisation
    !>
#ifdef __QUINCY_STANDALONE__
    ! SPP site-set - initialise mineral soil P pools
    IF (flag_spp1685) THEN
      DO iblk = 1,nblks
        ! for SPP sites, read in initial state from site-specific input data
        CALL init_sb_soil_p_pools_spp1685_sites( &
          & nproma, &
          & nsoil_sb, &
          & soil_lay_depth_center_sl(:,:,iblk), &
          & soil_p_labile(:,:,iblk), &
          & soil_p_slow(:,:,iblk), &
          & soil_p_occluded(:,:,iblk), &
          & soil_p_primary(:,:,iblk), &
          & site_ID_spp1685, &
          & p_labile_slow_site(:,:,iblk), &
          & p_labile(:,:,iblk), &
          & po4_assoc_slow(:,:,iblk), &
          & po4_occluded(:,:,iblk), &
          & po4_primary(:,:,iblk) )
      END DO
      DO iblk = 1,nblks
        DO ic = 1,nproma
          DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
            IF (flag_sb_double_langmuir) THEN
              po4_assoc_fast(ic, is, iblk) = calc_fast_po4(p_labile(ic, is, iblk), &
                &                            qmax_fast_po4(ic, is, iblk), &
                &                            km_fast_po4(ic, is, iblk))
            ELSE
              po4_assoc_fast(ic, is, iblk) = calc_fast_po4(p_labile(ic, is, iblk), &
                &                            qmax_po4(ic, is, iblk), &
                &                            km_adsorpt_po4_sl(ic, is, iblk) )
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      po4_solute(:,:,:) = p_labile(:,:,:) - po4_assoc_fast(:,:,:)
    ELSE
#endif

    ! initialise mineral soil P pools
    DO iblk = 1,nblks
      CALL init_soil_p_pools( &
        & nproma, &                         ! in
        & nsoil_sb, &
        & num_sl_above_bedrock(:, iblk), &
        & soil_lay_depth_center_sl(:,:,iblk), &
        & flag_sb_double_langmuir, &
        & soil_depth_sl(:,:,iblk), &
        & bulk_dens_corr_sl(:,:,iblk), &
        & root_fraction_sl(:,:,iblk), &
        & soil_p_labile(:,:,iblk), &
        & soil_p_slow(:,:,iblk), &
        & soil_p_occluded(:,:,iblk), &
        & soil_p_primary(:,:,iblk), &
        & soil_p_depth, &
        & p_labile_slow_site(:,:,iblk), &   ! in
        & po4_solute(:,:,iblk), &           ! inout
        & po4_assoc_fast(:,:,iblk), &
        & po4_assoc_slow(:,:,iblk), &
        & po4_occluded(:,:,iblk), &
        & po4_primary(:,:,iblk), &
        & qmax_po4(:,:,iblk), &
        & km_adsorpt_po4_sl(:,:,iblk), &
        & qmax_fast_po4(:,:,iblk), &
        & qmax_slow_po4(:,:,iblk), &
        & km_fast_po4(:,:,iblk), &
        & km_slow_po4(:,:,iblk) )           ! inout
    END DO

#ifdef __QUINCY_STANDALONE__
    END IF  ! flag_spp1685
#endif

    ! de-allocate local allocatable var
    DEALLOCATE(p_labile_slow_site)
  END SUBROUTINE sb_init_simple_sm

  ! ======================================================================================================= !
  !>init SB_ with the jsm - jena soil model
  !>
  SUBROUTINE sb_init_jsm( &
    & nproma, &
    & nsoil_sb, &
    & nblks, &
    & num_sl_above_bedrock, &
    & lctlib_growthform, &
    & sb_model_scheme, &
    & flag_sb_double_langmuir, &
    & qmax_org_fine_particle, &
    & soil_p_labile, &
    & soil_p_slow, &
    & soil_p_occluded, &
    & soil_p_primary, &
    & soil_p_depth, &
    & soil_depth_sl, &
    & soil_lay_width_sl, &
    & soil_lay_depth_ubound_sl, &
    & soil_lay_depth_center_sl, &
    & hlp_frac_c13, &
    & hlp_frac_c14, &
    & root_fraction_sl, &
    & clay_sl, &
    & silt_sl, &
    & w_soil_sl, &
    & bulk_dens_sl, &
    & ph_sl, &
    & qmax_po4_min_sl, &
    & qmax_po4_om_sl, &
    & volume_min_sl, &
    & po4_solute, &
    & bulk_soil_carbon_sl, &
    & soil_litter_carbon_sl, &
    & qmax_po4, &
    & qmax_nh4, &
    & qmax_fast_po4, &
    & qmax_slow_po4, &
    & km_fast_po4, &
    & km_slow_po4, &
    & km_adsorpt_po4_sl, &
    & po4_assoc_fast, &
    & po4_assoc_slow, &
    & po4_occluded, &
    & Qmax_AlFe_cor, &
    & po4_primary, &
    & sb_pool_mt_domain, &
    & qmax_org_min_sl, &
    & qmax_nh4_min_sl, &
    & nh4_assoc, &
    & nh4_solute, &
    & no3_solute, &
    & nh4_n15_assoc, &
    & nh4_n15_solute, &
    & no3_n15_solute, &
    & bulk_dens_corr_sl, &
    & qmax_org, &
    & km_adsorpt_nh4_sl, &
    & k_bioturb, &
    & enzyme_frac_poly, &
    & enzyme_frac_residue, &
    & enzyme_frac_poly_c, &
    & enzyme_frac_poly_n, &
    & enzyme_frac_poly_p, &
    & enzyme_frac_AP, &
    & flag_spp1685, &           ! optional in
    & site_ID_spp1685 )         ! optional in

    USE mo_q_sb_jsm_processes,            ONLY: calc_bulk_soil_carbon, calc_bulk_density_correction, &
      &                                         calc_qmax_bulk_density_correction, calc_Psorption_parameter, &
      &                                         calc_fast_po4
    USE mo_q_sb_jsm_transport,            ONLY: calc_bioturbation_rate
    USE mo_spq_util,                      ONLY: calc_qmax_texture
    USE mo_isotope_util,                  ONLY: calc_mixing_ratio_N15N14
    USE mo_veg_constants,                 ONLY: ITREE
    USE mo_sb_constants                   ! ..
    USE mo_jsb_physical_constants,        ONLY: molar_mass_C
#ifdef __QUINCY_STANDALONE__
    USE mo_qs_process_init_util,          ONLY: init_sb_soil_properties_spp1685_sites, &
      &                                         init_sb_soil_p_pools_spp1685_sites
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                    INTENT(in)    :: nproma
    INTEGER,                    INTENT(in)    :: nsoil_sb
    INTEGER,                    INTENT(in)    :: nblks
    REAL(wp),                   INTENT(in)    :: num_sl_above_bedrock(:,:)
    INTEGER,                    INTENT(in)    :: lctlib_growthform
    CHARACTER(len=*),           INTENT(in)    :: sb_model_scheme                  !< SB_ config
    LOGICAL,                    INTENT(in)    :: flag_sb_double_langmuir          !< SB_ config
    REAL(wp),                   INTENT(in)    :: qmax_org_fine_particle(:,:,:)    !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_labile(:,:,:)             !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_slow(:,:,:)               !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_occluded(:,:,:)           !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_primary(:,:,:)            !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_depth                     !< SB_ config
    REAL(wp),                   INTENT(in)    :: soil_depth_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_width_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_depth_ubound_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_depth_center_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: hlp_frac_c13
    REAL(wp),                   INTENT(in)    :: hlp_frac_c14
    REAL(wp),                   INTENT(in)    :: root_fraction_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: clay_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: silt_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: w_soil_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: bulk_dens_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: ph_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4_min_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4_om_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: volume_min_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_solute(:,:,:)
    REAL(wp),                   INTENT(inout) :: bulk_soil_carbon_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: soil_litter_carbon_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_nh4(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_fast_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_slow_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_fast_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_slow_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_adsorpt_po4_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_assoc_fast(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_assoc_slow(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_occluded(:,:,:)
    REAL(wp),                   INTENT(inout) :: Qmax_AlFe_cor(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_primary(:,:,:)
    REAL(wp),                   INTENT(inout) :: sb_pool_mt_domain(:,:,:,:,:)          !< bgc_material sb_pool
    REAL(wp),                   INTENT(out)   :: qmax_org_min_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_nh4_min_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_assoc(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: no3_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_n15_assoc(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_n15_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: no3_n15_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: bulk_dens_corr_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_org(:,:,:)
    REAL(wp),                   INTENT(out)   :: km_adsorpt_nh4_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: k_bioturb(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_poly(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_residue(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_poly_c(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_poly_n(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_poly_p(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_AP(:,:,:)
    LOGICAL,          OPTIONAL, INTENT(in)    :: flag_spp1685
    CHARACTER(len=3), OPTIONAL, INTENT(in)    :: site_ID_spp1685
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), ALLOCATABLE       :: p_labile(:,:,:)
    REAL(wp), ALLOCATABLE       :: p_labile_slow_site(:,:,:)
    REAL(wp), ALLOCATABLE       :: arr_hlp7(:,:,:)
    REAL(wp), ALLOCATABLE       :: arr_hlp8(:,:,:)
    REAL(wp), ALLOCATABLE       :: arr_hlp9(:,:,:)
    REAL(wp)                    :: frac_mineral_assoc_soil_c  !< fraction of mineral-associated C in the soil
    INTEGER                     :: iblk, is, ic
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init_jsm'
    IF (debug_on()) CALL message(TRIM(routine), 'Starting ...')
    ! ----------------------------------------------------------------------------------------------------- !

    !> init out var
    !>
    qmax_org_min_sl(:,:,:)      = 0.0_wp
    qmax_nh4_min_sl(:,:,:)      = 0.0_wp
    nh4_assoc(:,:,:)            = 0.0_wp
    nh4_solute(:,:,:)           = 0.0_wp
    no3_solute(:,:,:)           = 0.0_wp
    nh4_n15_assoc(:,:,:)        = 0.0_wp
    nh4_n15_solute(:,:,:)       = 0.0_wp
    no3_n15_solute(:,:,:)       = 0.0_wp
    bulk_dens_corr_sl(:,:,:)    = 0.0_wp
    qmax_org(:,:,:)             = 0.0_wp
    km_adsorpt_nh4_sl(:,:,:)    = 0.0_wp
    k_bioturb(:,:,:)            = 0.0_wp
    enzyme_frac_poly(:,:,:)     = 0.0_wp
    enzyme_frac_residue(:,:,:)  = 0.0_wp
    enzyme_frac_poly_c(:,:,:)   = 0.0_wp
    enzyme_frac_poly_n(:,:,:)   = 0.0_wp
    enzyme_frac_poly_p(:,:,:)   = 0.0_wp
    enzyme_frac_AP(:,:,:)       = 0.0_wp

    !>0.9 allocate local var
    !>
    ALLOCATE(p_labile(nproma, nsoil_sb, nblks))
    ALLOCATE(p_labile_slow_site(nproma, nsoil_sb, nblks))
    ALLOCATE(arr_hlp7(nproma, nsoil_sb, nblks))
    ALLOCATE(arr_hlp8(nproma, nsoil_sb, nblks))
    ALLOCATE(arr_hlp9(nproma, nsoil_sb, nblks))

    !>1.0 init
    !>

    ! initialisation of litter pools
    sb_pool_mt_domain(ix_soluable_litter, ixC, :, 1, :) &
      & = def_sb_pool_metabol_litter /  molar_mass_C / soil_lay_width_sl(:,1,:) &
      &                                             * (frac_prod_not_root + (1.0_wp - frac_prod_not_root) &
      &                                             * root_fraction_sl(:,1,:))
    sb_pool_mt_domain(ix_polymeric_litter, ixC, :, 1, :) &
      & = def_sb_pool_jsm_struct_litter / molar_mass_C / soil_lay_width_sl(:,1,:) &
      &                                             * (frac_prod_not_root + (1.0_wp - frac_prod_not_root) &
      &                                             * root_fraction_sl(:,1,:))
    IF (lctlib_growthform == ITREE) THEN
      sb_pool_mt_domain(ix_woody_litter, ixC, :, 1, :)  = def_sb_pool_woody_litter / molar_mass_C / soil_lay_width_sl(:,1,:)
    ENDIF
    DO is = 2,nsoil_sb
      sb_pool_mt_domain(ix_soluable_litter, ixC, :, is, :)  = def_sb_pool_metabol_litter /  molar_mass_C &
        &                                                 / soil_lay_width_sl(:, is, :) &
        &                                                 * (1.0_wp - frac_prod_not_root) * root_fraction_sl(:, is, :)
      sb_pool_mt_domain(ix_polymeric_litter, ixC, :, is, :) = def_sb_pool_jsm_struct_litter / molar_mass_C &
        &                                                 / soil_lay_width_sl(:, is, :) &
        &                                                 * (1.0_wp - frac_prod_not_root) * root_fraction_sl(:, is, :)
    ENDDO

    ! for each soil layer
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
          sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) = MAX(def_sb_pool_microbial_biomass &
            &                                          / soil_lay_width_sl(ic, is, iblk) &
            &                                          * root_fraction_sl(ic, is, iblk), min_mic_biomass)
          sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk) = def_sb_pool_dom / soil_lay_width_sl(ic, is, iblk) &
            &                                          * root_fraction_sl(ic, is, iblk)
          sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk) &
            & = k_sb_pool_dom_assoc_c_corr * sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk)
          sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) = def_sb_pool_microbial_necromass / soil_lay_width_sl(ic, is, iblk) &
            &                                              * root_fraction_sl(ic, is, iblk)
          ! calc fraction of mineral-associated C in the soil
          ! @TODO make the '7.0_wp' a parameter
          ! @TODO add docu for this calculation
          frac_mineral_assoc_soil_c  = EXP(-is / 7._wp)
          ! deep soil is initialised with very low MOC, corrected by k_deep_soil_moc_corr
          IF (is > 10._wp) THEN
            frac_mineral_assoc_soil_c = frac_mineral_assoc_soil_c * k_deep_soil_moc_corr
          ENDIF
          sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk) = qmax_org_fine_particle(ic, is, iblk) &
            &                                                      * (clay_sl(ic, is, iblk) + silt_sl(ic, is, iblk)) &
            &                                                      * bulk_dens_sl(ic, is, iblk) * frac_mineral_assoc_soil_c
          ! To avoid over-saturation of OM in soil layer
          ! calculate the fast C and stable C content (saturation at 1._wp) for each layer
          IF (is == 1) THEN
            arr_hlp7(ic, is, iblk) = (sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) &
              &                     + sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk) + &
              &                     sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk)) * &
              &                     molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
            arr_hlp8(ic, is, iblk) = (sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) &
              &                     + sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)) * &
              &                     molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
          ELSE
            arr_hlp7(ic, is, iblk) = (sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk)   &
              &                     + sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
              &                     + sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk)        &
              &                     + sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk)              &
              &                     + sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk))       &
              &                     * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
            arr_hlp8(ic, is, iblk) = (sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) &
              &                     + sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)) * &
              &                     molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
          ENDIF
          ! check if soil is over-saturated with fast C
          IF (arr_hlp7(ic, is, iblk) > 1._wp ) THEN
            sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk)  = sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk)  &
              &                                                         * 1._wp / arr_hlp7(ic, is, iblk)
            sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
              &                                                         * 1._wp / arr_hlp7(ic, is, iblk)
            sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk)        = sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk)        &
              &                                                         * 1._wp / arr_hlp7(ic, is, iblk)
            sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk)              = sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk)              &
              &                                                         * 1._wp / arr_hlp7(ic, is, iblk)
            sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk)        = sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk)        &
              &                                                         * 1._wp / arr_hlp7(ic, is, iblk)
            sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk)          = 0._wp
            sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)    = 0._wp
          ! check if soil is over-saturated with fast and stable C
          ELSE IF ((arr_hlp7(ic, is, iblk) + arr_hlp8(ic, is, iblk)) > 1._wp) THEN
            arr_hlp9(ic, is, iblk) = (1._wp - arr_hlp7(ic, is, iblk)) * fact_repartition_om_init_cond1 * &
              &                     qmax_org_fine_particle(ic, is, iblk) * &
              &                     (clay_sl(ic, is, iblk) + silt_sl(ic, is, iblk)) * bulk_dens_sl(ic, is, iblk)
            arr_hlp8(ic, is, iblk) = (1._wp - arr_hlp7(ic, is, iblk)) * (1._wp - fact_repartition_om_init_cond1) * &
              &                     rho_bulk_org * 1000._wp * carbon_per_dryweight_SOM / molar_mass_C
            IF (arr_hlp8(ic, is, iblk) > arr_hlp9(ic, is, iblk)) THEN
              sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) = &
                &                                                 sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) &
                &                                                 + (arr_hlp8(ic, is, iblk) - arr_hlp9(ic, is, iblk)) &
                &                                                 * fact_repartition_om_init_cond2
              sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = &
                &                                                 sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
                &                                                 + (arr_hlp8(ic, is, iblk) - arr_hlp9(ic, is, iblk)) &
                &                                                 * (1._wp - fact_repartition_om_init_cond2)
              arr_hlp8(ic, is, iblk) = arr_hlp9(ic, is, iblk)
            ENDIF
            sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk) = calc_fast_po4(arr_hlp8(ic, is, iblk), &
              &                                                      arr_hlp9(ic, is, iblk), &
              &                                                      k_desorpt_det / k_adsorpt_det)
            sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) &
              & = arr_hlp8(ic, is, iblk) - sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)
          ! check if the available mineral surface is enough for sorbed OM
          ELSE
            arr_hlp9(ic, is, iblk) = (1._wp - arr_hlp7(ic, is, iblk) - arr_hlp8(ic, is, iblk)) &
              &                     * qmax_org_fine_particle(ic, is, iblk) &
              &                     * (clay_sl(ic, is, iblk) + silt_sl(ic, is, iblk)) * bulk_dens_sl(ic, is, iblk)
            arr_hlp8(ic, is, iblk) &
              & = (sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) + sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk))
            IF (arr_hlp8(ic, is, iblk) > arr_hlp9(ic, is, iblk)) THEN
              sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) = &
                &                                               sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) &
                &                                               + (arr_hlp8(ic, is, iblk) - arr_hlp9(ic, is, iblk)) &
                &                                               * fact_repartition_om_init_cond2
              sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = &
                &                                               sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
                &                                               + (arr_hlp8(ic, is, iblk) - arr_hlp9(ic, is, iblk)) &
                &                                               * (1._wp - fact_repartition_om_init_cond2)
              arr_hlp8(ic, is, iblk) = arr_hlp9(ic, is, iblk)
            ENDIF
            sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk) = calc_fast_po4(arr_hlp8(ic, is, iblk), &
              &                                                      arr_hlp9(ic, is, iblk), &
              &                                                      k_desorpt_det / k_adsorpt_det)
            sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) &
              & = arr_hlp8(ic, is, iblk) - sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! initial values of soluble nh4 and no3
    DO iblk = 1,nblks
      DO ic = 1,nproma
        nh4_solute(ic, :, iblk)      = nh4_solute_prescribe / soil_lay_depth_ubound_sl(ic, nsoil_sb, iblk)
        nh4_n15_solute(ic, :, iblk)  = nh4_solute(ic, :, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        no3_solute(ic, :, iblk)      = no3_solute_prescribe / soil_lay_depth_ubound_sl(ic, nsoil_sb, iblk)
        no3_n15_solute(ic, :, iblk)  = no3_solute(ic, :, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
      END DO
    END DO
    ! inital bulk density, Qmax_org & Qmax_po4 based on current soil carbon content
    DO iblk = 1,nblks
      CALL calc_bulk_soil_carbon( &
        & nproma, &                           ! in
        & nsoil_sb, &
        & num_sl_above_bedrock(:,iblk), &
        & TRIM(sb_model_scheme), &
        & sb_pool_mt_domain(:,:,:,:,iblk), &         ! in
        & bulk_soil_carbon_sl(:,:,iblk), &    ! inout
        & soil_litter_carbon_sl(:,:,iblk), &
        & volume_min_sl(:,:,iblk) )           ! inout
      bulk_dens_corr_sl(:,:,iblk) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:,iblk), &
        &                                                        soil_litter_carbon_sl(:,:,iblk), &
        &                                                        bulk_dens_sl(:,:,iblk))
      !> po4/OC adsorption parameter intialisation
      !! OC adsorption
      qmax_org_min_sl(:,:,iblk) = calc_qmax_texture(qmax_org_fine_particle(:,:,iblk), &
        &                                           silt_sl(:,:,iblk), &
        &                                           clay_sl(:,:,iblk))
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &    ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & 0._wp, &
        & qmax_org_min_sl(:,:,iblk), &        ! in
        & qmax_org(:,:,iblk) )                ! out
    END DO

    !! po4 adsorption

#ifdef __QUINCY_STANDALONE__
    ! specific init values for sites of the SPP site-set
    ! for: ph_sl(:,:,:) and Qmax_AlFe_cor(:,:,:)
    IF (flag_spp1685) THEN
      DO iblk = 1,nblks
        CALL init_sb_soil_properties_spp1685_sites( &
          & nproma, &                   ! in
          & nsoil_sb, &
          & site_ID_spp1685, &          ! in
          & ph_sl(:,:,iblk), &          ! out
          & Qmax_AlFe_cor(:,:,iblk) )   ! out
      END DO
    ENDIF
#endif

    ! Correct the first layer Qmax_AlFe_cor(:,:,:) with the woody litter
    Qmax_AlFe_cor(:,1,:) = Qmax_AlFe_cor(:,1,:) * 1._wp / ((1._wp - volume_min_sl(:,1,:)) &
      &                    * (bulk_soil_carbon_sl(:,1,:) + soil_litter_carbon_sl(:,1,:) &
      &                    + sb_pool_mt_domain(ix_woody_litter, ixC, :, 1, :)) &
      &                    / (bulk_soil_carbon_sl(:,1,:) + soil_litter_carbon_sl(:,1,:)) &
      &                    + volume_min_sl(:,1,:))
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
          CALL calc_Psorption_parameter( &
            & flag_sb_double_langmuir, &          ! in
            & clay_sl(ic, is, iblk), &
            & silt_sl(ic, is, iblk), &
            & w_soil_sl(ic, is, iblk) / soil_depth_sl(ic, is, iblk) * 1000._wp, &
            & ph_sl(ic, is, iblk), &
            & bulk_soil_carbon_sl(ic, is, iblk), &
            & soil_litter_carbon_sl(ic, is, iblk), &
            & sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk), &
            & bulk_dens_sl(ic, is, iblk), &
            & volume_min_sl(ic, is, iblk), &
            & po4_solute_prescribe, &             ! in - po4_solute(ic, is, iblk)
            & qmax_po4_min_sl(ic, is, iblk), &        ! inout
            & qmax_po4_om_sl(ic, is, iblk), &
            & qmax_po4(ic, is, iblk), &
            & km_adsorpt_po4_sl(ic, is, iblk), &
            & qmax_fast_po4(ic, is, iblk), &
            & qmax_slow_po4(ic, is, iblk), &
            & km_fast_po4(ic, is, iblk), &
            & km_slow_po4(ic, is, iblk), &
            & Qmax_AlFe_cor(ic, is, iblk) )           ! inout
        ENDDO
      ENDDO
    ENDDO

    ! calculate the partition coefficient between labile and slow P pool based on
    ! their Qmax (double-Langmuir) or global average values (single-Langmuir)
    IF (flag_sb_double_langmuir) THEN
      p_labile_slow_site(:,:,:) = eps4
      ! if 'qmax_po4(:,:,:) > eps8'
      DO iblk = 1,nblks
        DO is = 1,nsoil_sb
          DO ic = 1,nproma
            IF (qmax_po4(ic, is, iblk) > eps8) THEN
              p_labile_slow_site(ic, is, iblk) = MAX(qmax_fast_po4(ic, is, iblk) / qmax_po4(ic, is, iblk), eps4)
            END IF
          END DO
        END DO
      END DO
    ELSE
      p_labile_slow_site(:,:,:) = p_labile_slow_global_avg
    ENDIF

    !> nh4 adsorption parameter intialisation
    !!
    qmax_nh4_min_sl(:,:,:) = calc_qmax_texture(qmax_nh4_clay, &
      &                                        0.0_wp, &
      &                                        clay_sl(:,:,:))
    DO iblk = 1,nblks
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &   ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & 0.0_wp, &
        & qmax_nh4_min_sl(:,:,iblk), &       ! in
        & qmax_nh4(:,:,iblk) )               ! out
    END DO
    qmax_nh4(:,:,:) = qmax_nh4(:,:,:) * w_soil_sl(:,:,:) / soil_depth_sl(:,:,:)
    DO iblk = 1,nblks
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &   ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & km_adsorpt_OM_nh4, &
        & km_adsorpt_mineral_nh4, &          ! in
        & km_adsorpt_nh4_sl(:,:,iblk) )      ! out
    END DO
    nh4_assoc(:,:,:)     = qmax_nh4(:,:,:) * nh4_solute (:,:,:) / (nh4_solute(:,:,:) + km_adsorpt_nh4_sl(:,:,:))
    nh4_n15_assoc(:,:,:) = nh4_assoc(:,:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

    !> Soil P pools intialisation
    !>
#ifdef __QUINCY_STANDALONE__
    ! SPP site-set - initialise mineral soil P pools
    IF (flag_spp1685) THEN
      DO iblk = 1,nblks
        ! for SPP sites, read in initial state from site-specific input data
        CALL init_sb_soil_p_pools_spp1685_sites( &
          & nproma, &
          & nsoil_sb, &
          & soil_lay_depth_center_sl(:,:,iblk), &
          & soil_p_labile(:,:,iblk), &
          & soil_p_slow(:,:,iblk), &
          & soil_p_occluded(:,:,iblk), &
          & soil_p_primary(:,:,iblk), &
          & site_ID_spp1685, &
          & p_labile_slow_site(:,:,iblk), &
          & p_labile(:,:,iblk), &
          & po4_assoc_slow(:,:,iblk), &
          & po4_occluded(:,:,iblk), &
          & po4_primary(:,:,iblk) )
      END DO
      ! calculate the soluble and adsorbed P pool based on Langmuir isotherm
      DO iblk = 1,nblks
        DO ic = 1,nproma
          DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
            IF (flag_sb_double_langmuir) THEN
              po4_assoc_fast(ic, is, iblk) = calc_fast_po4(p_labile(ic, is, iblk), &
                &                            qmax_fast_po4(ic, is, iblk), &
                &                            km_fast_po4(ic, is, iblk))
            ELSE
              po4_assoc_fast(ic, is, iblk) = calc_fast_po4(p_labile(ic, is, iblk), &
                &                            qmax_po4(ic, is, iblk), &
                &                            km_adsorpt_po4_sl(ic, is, iblk))
            ENDIF
            po4_solute(ic, is, iblk) = p_labile(ic, is, iblk) - po4_assoc_fast(ic, is, iblk)
          ENDDO
        ENDDO
      ENDDO
    ELSE
#endif

    ! initialise mineral soil P pools
    DO iblk = 1,nblks
      CALL init_soil_p_pools( &
        & nproma, &                               ! in
        & nsoil_sb, &
        & num_sl_above_bedrock(:, iblk), &
        & soil_lay_depth_center_sl(:,:,iblk), &
        & flag_sb_double_langmuir, &
        & soil_depth_sl(:,:,iblk), &
        & bulk_dens_corr_sl(:,:,iblk), &
        & root_fraction_sl(:,:,iblk), &
        & soil_p_labile(:,:,iblk), &
        & soil_p_slow(:,:,iblk), &
        & soil_p_occluded(:,:,iblk), &
        & soil_p_primary(:,:,iblk), &
        & soil_p_depth, &
        & p_labile_slow_site(:,:,iblk), &         ! in
        & po4_solute(:,:,iblk), &                 ! inout
        & po4_assoc_fast(:,:,iblk), &
        & po4_assoc_slow(:,:,iblk), &
        & po4_occluded(:,:,iblk), &
        & po4_primary(:,:,iblk), &
        & qmax_po4(:,:,iblk), &
        & km_adsorpt_po4_sl(:,:,iblk), &
        & qmax_fast_po4(:,:,iblk), &
        & qmax_slow_po4(:,:,iblk), &
        & km_fast_po4(:,:,iblk), &
        & km_slow_po4(:,:,iblk) )                 ! inout
    END DO

#ifdef __QUINCY_STANDALONE__
    END IF  ! flag_spp1685
#endif

    ! initialise diffusion by bioturbation
    DO iblk = 1,nblks
      k_bioturb(:,:,iblk) = calc_bioturbation_rate(nproma, &
        &                                          nsoil_sb, &
        &                                          num_sl_above_bedrock(:,iblk), &
        &                                          soil_depth_sl(:,:,iblk), &
        &                                          bulk_dens_corr_sl(:,:,iblk), &
        &                                          root_fraction_sl(:,:,iblk))
    END DO

    ! initialisation of fractional enzyme allocation (range [0,1])
    enzyme_frac_poly(:,:,:)     = init_enzyme_fraction
    enzyme_frac_residue(:,:,:)  = init_enzyme_fraction
    enzyme_frac_poly_c(:,:,:)   = init_enzyme_fraction
    enzyme_frac_poly_n(:,:,:)   = init_enzyme_fraction
    enzyme_frac_poly_p(:,:,:)   = init_enzyme_fraction
    enzyme_frac_AP(:,:,:)       = min_enzyme_fraction

    ! isotopic signals in all pools
    ! microbial, DOM, aDOM, residue, and aRes pools
    sb_pool_mt_domain(ix_soluable_litter, ixN, :, :, :) &
      & = sb_pool_mt_domain(ix_soluable_litter, ixC, :, :, :) / k_init_soluable_litter_cn
    sb_pool_mt_domain(ix_soluable_litter, ixP, :, :, :) &
      &  = sb_pool_mt_domain(ix_soluable_litter, ixN, :, :, :) / k_init_litter_np
    sb_pool_mt_domain(ix_soluable_litter, ixC13, :, :, :) = sb_pool_mt_domain(ix_soluable_litter, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_soluable_litter, ixC14, :, :, :) = sb_pool_mt_domain(ix_soluable_litter, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_soluable_litter, ixN15, :, :, :) = sb_pool_mt_domain(ix_soluable_litter, ixN, :, :, :) / &
      &                                              (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_polymeric_litter, ixN, :, :, :) &
      & = sb_pool_mt_domain(ix_polymeric_litter, ixC, :, :, :) / k_init_polymeric_litter_cn
    sb_pool_mt_domain(ix_polymeric_litter, ixP, :, :, :) &
      & = sb_pool_mt_domain(ix_polymeric_litter, ixN, :, :, :) / k_init_litter_np
    sb_pool_mt_domain(ix_polymeric_litter, ixC13, :, :, :) = sb_pool_mt_domain(ix_polymeric_litter, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_polymeric_litter, ixC14, :, :, :) = sb_pool_mt_domain(ix_polymeric_litter, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_polymeric_litter, ixN15, :, :, :) = sb_pool_mt_domain(ix_polymeric_litter, ixN, :, :, :) / &
      &                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_woody_litter, ixN, :, :, :) &
      & = sb_pool_mt_domain(ix_woody_litter, ixC, :, :, :) / k_init_woody_litter_cn
    sb_pool_mt_domain(ix_woody_litter, ixP, :, :, :)   = sb_pool_mt_domain(ix_woody_litter, ixN, :, :, :) / k_init_litter_np
    sb_pool_mt_domain(ix_woody_litter, ixC13, :, :, :) = sb_pool_mt_domain(ix_woody_litter, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_woody_litter, ixC14, :, :, :) = sb_pool_mt_domain(ix_woody_litter, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_woody_litter, ixN15, :, :, :) = sb_pool_mt_domain(ix_woody_litter, ixN, :, :, :) / &
      &                                           (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_microbial, ixN, :, :, :)   = sb_pool_mt_domain(ix_microbial, ixC, :, :, :) / microbial_cn
    sb_pool_mt_domain(ix_microbial, ixP, :, :, :)   = sb_pool_mt_domain(ix_microbial, ixN, :, :, :) / microbial_np
    sb_pool_mt_domain(ix_microbial, ixC13, :, :, :) = sb_pool_mt_domain(ix_microbial, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_microbial, ixC14, :, :, :) = sb_pool_mt_domain(ix_microbial, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_microbial, ixN15, :, :, :) = sb_pool_mt_domain(ix_microbial, ixN, :, :, :) / &
      &                                        (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_residue, ixN, :, :, :)   = sb_pool_mt_domain(ix_residue, ixC, :, :, :) / microbial_cn
    sb_pool_mt_domain(ix_residue, ixP, :, :, :)   = sb_pool_mt_domain(ix_residue, ixN, :, :, :) / k_slow_som_np
    sb_pool_mt_domain(ix_residue, ixC13, :, :, :) = sb_pool_mt_domain(ix_residue, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_residue, ixC14, :, :, :) = sb_pool_mt_domain(ix_residue, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_residue, ixN15, :, :, :) = sb_pool_mt_domain(ix_residue, ixN, :, :, :) / &
      &                                      (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_residue_assoc, ixN, :, :, :)   = sb_pool_mt_domain(ix_residue_assoc, ixC, :, :, :) / microbial_cn
    sb_pool_mt_domain(ix_residue_assoc, ixP, :, :, :)   = sb_pool_mt_domain(ix_residue_assoc, ixN, :, :, :) / k_slow_som_np
    sb_pool_mt_domain(ix_residue_assoc, ixC13, :, :, :) = sb_pool_mt_domain(ix_residue_assoc, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_residue_assoc, ixC14, :, :, :) = sb_pool_mt_domain(ix_residue_assoc, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_residue_assoc, ixN15, :, :, :) = sb_pool_mt_domain(ix_residue_assoc, ixN, :, :, :) / &
      &                                            (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_dom, ixN, :, :, :)   = sb_pool_mt_domain(ix_dom, ixC, :, :, :) / microbial_cn
    sb_pool_mt_domain(ix_dom, ixP, :, :, :)   = sb_pool_mt_domain(ix_dom, ixN, :, :, :) / microbial_np
    sb_pool_mt_domain(ix_dom, ixC13, :, :, :) = sb_pool_mt_domain(ix_dom, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_dom, ixC14, :, :, :) = sb_pool_mt_domain(ix_dom, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_dom, ixN15, :, :, :) = sb_pool_mt_domain(ix_dom, ixN, :, :, :) / &
      &                                  (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_dom_assoc, ixN, :, :, :)   = sb_pool_mt_domain(ix_dom_assoc, ixC, :, :, :) / microbial_cn
    sb_pool_mt_domain(ix_dom_assoc, ixP, :, :, :)   = sb_pool_mt_domain(ix_dom_assoc, ixN, :, :, :) / microbial_np
    sb_pool_mt_domain(ix_dom_assoc, ixC13, :, :, :) = sb_pool_mt_domain(ix_dom_assoc, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_dom_assoc, ixC14, :, :, :) = sb_pool_mt_domain(ix_dom_assoc, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_dom_assoc, ixN15, :, :, :) = sb_pool_mt_domain(ix_dom_assoc, ixN, :, :, :) / &
      &                                        (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    ! fungi and mycorrhiza
    sb_pool_mt_domain(ix_fungi, ixC13, :, :, :) = sb_pool_mt_domain(ix_fungi, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_fungi, ixC14, :, :, :) = sb_pool_mt_domain(ix_fungi, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_fungi, ixN15, :, :, :) = sb_pool_mt_domain(ix_fungi, ixN, :, :, :) / &
      &                                    (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_mycorrhiza, ixC13, :, :, :) = sb_pool_mt_domain(ix_mycorrhiza, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_mycorrhiza, ixC14, :, :, :) = sb_pool_mt_domain(ix_mycorrhiza, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_mycorrhiza, ixN15, :, :, :) = sb_pool_mt_domain(ix_mycorrhiza, ixN, :, :, :) / &
      &                                         (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

    ! de-allocate local allocatable var
    DEALLOCATE(p_labile)
    DEALLOCATE(p_labile_slow_site)
    DEALLOCATE(arr_hlp7)
    DEALLOCATE(arr_hlp8)
    DEALLOCATE(arr_hlp9)
  END SUBROUTINE sb_init_jsm

  ! ======================================================================================================= !
  !>Intialize soil P pools for all layers from the soil P maps (Yang et al. 2013)
  !> assuming that the soil P pool sizes are correlated with soil weight
  !>
  SUBROUTINE init_soil_p_pools(nc, &
                               nsoil_sb, &
                               num_sl_above_bedrock, &
                               soil_lay_depth_center_sl, &
                               flag_sb_double_langmuir, &
                               soil_depth_sl, bulk_dens_sl, root_fraction_sl,&
                               soil_p_labile, soil_p_slow, soil_p_occluded, soil_p_primary, &
                               soil_p_depth, p_labile_slow_ratio, &
                               po4_solute, po4_assoc_fast, po4_assoc_slow, po4_occluded, po4_primary, &
                               qmax_po4,km_po4, qmax_fast_po4, km_fast_po4, qmax_slow_po4, km_slow_po4)

    USE mo_jsb_physical_constants, ONLY: molar_mass_P
    USE mo_sb_constants           ! e.g.: k1_p_pool_depth_corr, k2_p_pool_depth_corr, k3_p_pool_depth_corr, &
                                  !       k4_p_pool_depth_corr, k_p_pool_depth_corr_min_occlud, p_labile_slow_global_avg
    USE mo_q_sb_jsm_processes,     ONLY: calc_fast_po4

    INTEGER,                          INTENT(in)     :: nc                        !< dimensions
    INTEGER                         , INTENT(in)     :: nsoil_sb                  !< number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(nc)         , INTENT(in)     :: num_sl_above_bedrock      !< number of soil layers with thickness > eps8
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)     :: soil_lay_depth_center_sl  !< depth at the center of each soil layer
    LOGICAL                         , INTENT(in)     :: flag_sb_double_langmuir   !< T: double Langmuir; F: traditional Langmuir
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)     :: soil_depth_sl, &      !< soil depth per soil layer, [m]
                                                        root_fraction_sl, &   !<
                                                        bulk_dens_sl, &       !< corrected soil bulk density per soil layer, [kg m-3]
                                                        p_labile_slow_ratio   !< ratio between labile and labile+slow P pools
    REAL(wp)                        , INTENT(in)     :: soil_p_depth          !< soil depth of the soil P maps
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(in)     :: soil_p_labile, &      !< total labile P pool (solute + assoc_fast) of the soil P maps, [g m-2]
                                                        soil_p_slow, &        !< total slow P pool (secondary P) of the soil P maps, [g m-2]
                                                        soil_p_occluded, &    !< total occluded P pool of the soil P maps, [g m-2]
                                                        soil_p_primary        !< total primary P pool of the soil P maps, [g m-2]
    REAL(wp), DIMENSION(nc,nsoil_sb), INTENT(inout)  :: po4_solute, &         !< soluble Pi pool per soil layer, [mol m-3]
                                                        po4_assoc_fast, &     !< fast minerally associated PO4 pool per soil layer, [mol m-3]
                                                        po4_assoc_slow, &     !< slow minerally associated PO4 pool per soil layer, [mol m-3]
                                                        po4_occluded, &       !< occluded PO4 pool per soil layer, [mol m-3]
                                                        po4_primary, &        !< primary PO4 pool per soil layer, [mol m-3]
                                                        qmax_po4, &           !< maximum PO4 sorption capacity per soil layer, [mol m-3]
                                                        km_po4, &             !< half-saturation PO4 concentration for sorption per soil layer, [mol m-3]
                                                        qmax_fast_po4, &      !< maximum PO4 sorption capacity for po4_assoc_fast, [mol m-3]
                                                        km_fast_po4, &        !< half-saturation PO4 concentration for po4_assoc_fast, [mol m-3]
                                                        qmax_slow_po4, &      !< maximum PO4 sorption capacity for po4_assoc_slow, [mol m-3]
                                                        km_slow_po4           !< half-saturation PO4 concentration for po4_assoc_slow, [mol m-3]
    INTEGER                                          :: ic, is
    REAL(wp), DIMENSION(nc)                          :: rhlp1
    REAL(wp), DIMENSION(nc)                          :: soil_p_weight         !< total soil weight up to the soil_p_depth [kg]
    REAL(wp), DIMENSION(nc)                          :: soil_p_layer          !< up to soil_p_layer keep the P init map info
    REAL(wp), DIMENSION(nc)                          :: hlp1, hlp2
    REAL(wp), DIMENSION(nc)                          :: hlp3, hlp4
    LOGICAL,  DIMENSION(nc)                          :: aflag
    REAL(wp), DIMENSION(nc,nsoil_sb)                 :: soil_weight_sl, &                 !< soil weight per soil layer
                                                        arr_hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_p_pools'


    !> 0.9 init output variables
    !!
    po4_solute(:,:)     = 0._wp
    po4_assoc_fast(:,:) = 0._wp
    po4_assoc_slow(:,:) = 0._wp
    po4_occluded(:,:)   = 0._wp
    po4_primary(:,:)    = 0._wp

    !> 1.0 calculate the total soil weight up to the soil_p_depth,
    !! and the soil weight per soil layer
    hlp1(:)  = 0._wp
    hlp2(:)  = 0._wp
    aflag(:) = .TRUE.
    !> NOTE at the moment the soil_p_* values are identical across soil layers, hence, using the value from sl 1
    rhlp1(:) = soil_p_slow(:,1) + soil_p_labile(:,1) + soil_p_occluded(:,1) + soil_p_primary(:,1)
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        soil_weight_sl(ic,is) = soil_depth_sl(ic,is) * bulk_dens_sl(ic,is)
        IF (aflag(ic)) THEN
          hlp1(ic) = hlp1(ic) + soil_depth_sl(ic,is)
          hlp2(ic) = hlp2(ic) + soil_weight_sl(ic,is)
          IF (hlp1(ic) > soil_p_depth) THEN
            soil_p_weight(ic) = hlp2(ic) - soil_weight_sl(ic,is) * (hlp1(ic) - soil_p_depth) / soil_depth_sl(ic,is)
            aflag(ic)         = .FALSE.
            soil_p_layer(ic)  = is
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    hlp4(:) = rhlp1(:) / soil_p_weight(:) ! the P density from init map, [g P / kg soil]

    !> 2.0 calculate the initial slow P, occluded P, primary P pool, and p_labile
    !! based on the assumption that soil P pool are linearly corrected with soil weight
    !! the depth correction of initial P pool distributions are fitted against the observed pattern of a young cambisol soil
    !! with basalt as parent material (BBR, Lang et al. 2017) that,
    !! the proportion of: labile+sorbed P decreases with increasing depth, and primary P increases with increasing depth (Yu et al., 2020)
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        ! primary P, corrected with depth, [g P]
        hlp1(ic) = MIN(k2_p_pool_depth_corr * soil_lay_depth_center_sl(ic,is) ** k4_p_pool_depth_corr, &
                      (k1_p_pool_depth_corr - k_p_pool_depth_corr_min_occlud)) / k1_p_pool_depth_corr
        hlp2(ic) = k3_p_pool_depth_corr ** (k2_p_pool_depth_corr * soil_lay_depth_center_sl(ic,is)) / k1_p_pool_depth_corr
        ! occluded P, corrected with depth
        hlp3(ic) = MAX((1.0_wp - hlp1(ic) - hlp2(ic)), k_p_pool_depth_corr_min_occlud / k1_p_pool_depth_corr)

        IF ((hlp1(ic) + hlp2(ic) + hlp3(ic)) > 1._wp) THEN
          hlp1(ic) = hlp1(ic) / (hlp1(ic) + hlp2(ic) + hlp3(ic))
          hlp2(ic) = hlp2(ic) / (hlp1(ic) + hlp2(ic) + hlp3(ic))
          hlp3(ic) = hlp3(ic) / (hlp1(ic) + hlp2(ic) + hlp3(ic))
        ENDIF
        IF (soil_p_layer(ic) > is) THEN
          po4_assoc_slow(ic,is) = bulk_dens_sl(ic,is) * hlp4(ic) * (soil_p_slow(ic,is) + soil_p_labile(ic,is)) &
            &                     / rhlp1(ic) / molar_mass_P * (1._wp - p_labile_slow_ratio(ic,is))
          po4_occluded(ic,is)   = bulk_dens_sl(ic,is) * hlp4(ic) * soil_p_occluded(ic,is) &
            &                     / rhlp1(ic) / molar_mass_P
          po4_primary(ic,is)    = bulk_dens_sl(ic,is) * hlp4(ic) * soil_p_primary(ic,is) &
            &                     / rhlp1(ic) / molar_mass_P
          arr_hlp1(ic,is)       = bulk_dens_sl(ic,is) * hlp4(ic) * (soil_p_slow(ic,is) + soil_p_labile(ic,is)) &
            &                     / rhlp1(ic) / molar_mass_P * p_labile_slow_ratio(ic,is)
        ELSE
          po4_assoc_slow(ic,is) = bulk_dens_sl(ic,is) * hlp4(ic) * hlp2(ic) / molar_mass_P &
            &                     * (1._wp - p_labile_slow_ratio(ic,is))
          po4_occluded(ic,is)   = bulk_dens_sl(ic,is) * hlp4(ic) * hlp3(ic) / molar_mass_P
          po4_primary(ic,is)    = bulk_dens_sl(ic,is) * hlp4(ic) * hlp1(ic) / molar_mass_P
          arr_hlp1(ic,is)       = bulk_dens_sl(ic,is) * hlp4(ic) * hlp2(ic) / molar_mass_P &
            &                     * p_labile_slow_ratio(ic,is)
        ENDIF
      ENDDO
    ENDDO

    !> 2.1 calculate the initial solute P and fast P,
    !! based on the assumption that labile P is the sum of both,
    !! and they follow the Langmuir equilibrium
    !! po4_assoc_fast = (p_labile_sl + Smax + K - sqrt((p_labile_sl + Smax + K)**2._wp - 4*p_labile_sl*Smax)) / 2._wp
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        IF (flag_sb_double_langmuir) THEN
          po4_assoc_fast(ic, is) = calc_fast_po4(arr_hlp1(ic, is), qmax_fast_po4(ic, is), km_fast_po4(ic, is))
          po4_solute(ic, is)     = arr_hlp1(ic, is) - po4_assoc_fast(ic, is)
        ELSE
          po4_assoc_fast(ic, is) = calc_fast_po4(arr_hlp1(ic, is), qmax_po4(ic, is), km_po4(ic, is))
          po4_solute(ic, is)     = arr_hlp1(ic, is) - po4_assoc_fast(ic, is)
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE init_soil_p_pools

#endif
END MODULE mo_sb_init
