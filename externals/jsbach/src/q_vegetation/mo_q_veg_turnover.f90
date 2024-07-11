!> QUINCY calculate vegetation turnover and litterfall
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
!>#### calculate vegetation turnover and litterfall from tissue senesence, including nutrient recycling
!>
!> also calculate turnover of sapwood to heartwood and fruits to seed_bed
!>
MODULE mo_q_veg_turnover
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message
  USE mo_jsb_impl_constants,  ONLY: test_false_true
  USE mo_jsb_math_constants,  ONLY: dtime, one_day, one_year, eps8, eps4

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, VEG_BGCM_GROWTH_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_veg_turnover

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_turnover'

CONTAINS

  ! ======================================================================================================= !
  !>Calculates litter fall resulting from tissue senesence, including nutrient recycling,
  !> as well as turnover of sapwood to heartwood and fruits to seed_bed
  !>
  SUBROUTINE update_veg_turnover(tile, options)
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: VEG_, Q_PHENO_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(Q_PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER :: model                  !< the model
    TYPE(t_lnd_bgcm_store),   POINTER :: bgcm_store             !< the bgcm store of this tile
    TYPE(t_lctlib_element),   POINTER :: lctlib                 !< land-cover-type library - parameter across pft's
    INTEGER                           :: iblk, ics, ice, nc     !< dimensions
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_turnover'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L2D :: veg_growth_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(Q_PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_PHENO_ 2D
    dsl4jsb_Real2D_onChunk      :: growing_season
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: target_lai
    dsl4jsb_Real2D_onChunk      :: target_cn_fine_root
    dsl4jsb_Real2D_onChunk      :: target_np_fine_root
    dsl4jsb_Real2D_onChunk      :: target_cn_leaf
    dsl4jsb_Real2D_onChunk      :: target_np_leaf
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_p
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_p
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n15
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n15
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n15
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_p
    dsl4jsb_Real2D_onChunk      :: net_growth
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(Q_PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_GROWTH_ID, veg_growth_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)           ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, target_lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_fine_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_fine_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth)
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 Calculate nutrient flux from enzyme, RNA and other turnover in leaves and roots
    !>
    CALL calc_nutrient_recycling( &
      & nc                                    , & ! in
      & target_cn_fine_root(:)                , &
      & target_np_fine_root(:)                , &
      & target_cn_leaf(:)                     , &
      & target_np_leaf(:)                     , &
      & growing_season(:)                     , &
      & veg_pool_mt(:,:,:)                    , & ! in
      & recycling_leaf_n(:)                   , & ! inout
      & recycling_fine_root_n(:)              , &
      & recycling_leaf_p(:)                   , &
      & recycling_fine_root_p(:)              , &
      & recycling_leaf_n15(:)                 , &
      & recycling_fine_root_n15(:)            )   ! inout

    !>2.0 Calculate recycling and retranslocation
    !>
    CALL calc_turnover( &
      & nc                                    , & ! in
      & model%config%elements_index_map(:)    , &
      & model%config%is_element_used(:)       , &
      & lctlib%phenology_type                 , &
      & lctlib%growthform                     , &
      & lctlib%frac_sapwood_branch            , &
      & lctlib%tau_leaf                       , &
      & lctlib%tau_fine_root                  , &
      & lctlib%tau_coarse_root                , &
      & lctlib%tau_branch                     , &
      & lctlib%tau_sap_wood                   , &
      & lctlib%tau_fruit                      , &
      & lctlib%tau_seed_litter                , &
      & growing_season(:)                     , & ! in
      & lai(:)                                , & ! in
      & target_lai(:)                         , & ! in
      & veg_pool_mt(:,:,:)                    , & ! in
      & recycling_leaf_n(:)                   , & ! inout
      & recycling_leaf_n15(:)                 , &
      & recycling_leaf_p(:)                   , &
      & recycling_heart_wood_n(:)             , &
      & recycling_heart_wood_n15(:)           , &
      & recycling_heart_wood_p(:)             , &
      & veg_growth_mt(:,:,:)                  , & ! inout
      & veg_litterfall_mt(:,:,:)              , & ! inout
      & net_growth(:)                          )  ! out

  END SUBROUTINE update_veg_turnover

  ! ======================================================================================================= !
  !>Calculate nutrient flux from enzyme turnover and other biochemical turnover in leaves and roots
  !>  at the time scale of tau_nutrient_recycling
  !>
  !>  Calculations only during growing season
  !>
  !>  Input: plant status and leaf and root C:N:P target ratios
  !>  Output: net flux between labile and leaves/roots (+ve: flux to labile, -ve: flux from labile)
  !>
  SUBROUTINE calc_nutrient_recycling( &
    & nc                                , &
    & target_cn_fine_root               , &
    & target_np_fine_root               , &
    & target_cn_leaf                    , &
    & target_np_leaf                    , &
    & growing_season                    , &
    & veg_pool_mt                       , &
    & recycling_leaf_n                  , &
    & recycling_fine_root_n             , &
    & recycling_leaf_p                  , &
    & recycling_fine_root_p             , &
    & recycling_leaf_n15                , &
    & recycling_fine_root_n15           )

    USE mo_veg_constants,         ONLY: tau_nutrient_recycling, tau_labile
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                          INTENT(in)    :: nc                         !< dimensions
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: target_cn_fine_root        !< target C:N of fine_roots
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: target_np_fine_root        !< target N:P of fine_roots
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: target_cn_leaf             !< ..
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: target_np_leaf             !< ..
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: growing_season             !< whether the plant is in the growing season
    REAL(wp),                         INTENT(in)    :: veg_pool_mt(:,:,:)         !< bgcm veg_pool
    REAL(wp), DIMENSION(nc),          INTENT(inout) :: recycling_leaf_n           !< ..
    REAL(wp), DIMENSION(nc),          INTENT(inout) :: recycling_fine_root_n      !< ..
    REAL(wp), DIMENSION(nc),          INTENT(inout) :: recycling_leaf_p           !< ..
    REAL(wp), DIMENSION(nc),          INTENT(inout) :: recycling_fine_root_p      !< ..
    REAL(wp), DIMENSION(nc),          INTENT(inout) :: recycling_leaf_n15         !< ..
    REAL(wp), DIMENSION(nc),          INTENT(inout) :: recycling_fine_root_n15    !< ..
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc)     :: hlp1                     !< helper variable
    REAL(wp), DIMENSION(nc)     :: recycling_leaf           !< more helpers1
    REAL(wp), DIMENSION(nc)     :: recycling_fine_root      !< more helpers2
    INTEGER                     :: ic                       !< loop over chunk of points
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_nutrient_recycling'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 replace overturning N from leaves and roots by new N from the labile pool using the target C:N
    !>
    DO ic = 1,nc
      IF (growing_season(ic) > test_false_true) THEN
        recycling_leaf(ic)      = (veg_pool_mt(ix_leaf, ixN, ic) &
          &                       - veg_pool_mt(ix_leaf, ixC, ic) / target_cn_leaf(ic)) / tau_nutrient_recycling &
          &                       / one_day * dtime
        recycling_fine_root(ic) = (veg_pool_mt(ix_fine_root, ixN, ic) &
          &                       - veg_pool_mt(ix_fine_root, ixC, ic) / target_cn_fine_root(ic)) / tau_nutrient_recycling &
          &                       / one_day * dtime
        !>  1.1 if N is addded to leaves or fine root (negative flux), make sure there is enough labile N to support this
        !>
        IF ((-1.0_wp * (recycling_leaf(ic) + recycling_fine_root(ic))) > eps8) THEN
          hlp1(ic) = -1.0_wp * (1.0_wp / tau_labile / one_day * veg_pool_mt(ix_labile, ixN, ic) * dtime) &
            &        / (recycling_leaf(ic) + recycling_fine_root(ic))
          IF (hlp1(ic) < 1.0_wp) THEN
            recycling_leaf(ic)       = recycling_leaf(ic)      * hlp1(ic)
            recycling_fine_root(ic)  = recycling_fine_root(ic) * hlp1(ic)
          END IF
        END IF
        recycling_leaf_n(ic)      = recycling_leaf_n(ic)      + recycling_leaf(ic)
        recycling_fine_root_n(ic) = recycling_fine_root_n(ic) + recycling_fine_root(ic)
        !>  1.2 calculate signature of N transfer from labile to leaf or inverse
        !>
        IF (recycling_leaf(ic) >= 0.0_wp) THEN  ! for all positive values of recycling_fine_root(ic)
          IF (veg_pool_mt(ix_leaf, ixN, ic) > eps8) THEN
            recycling_leaf_n15(ic) = recycling_leaf_n15(ic) &
              &                      + recycling_leaf(ic) / veg_pool_mt(ix_leaf, ixN, ic) * veg_pool_mt(ix_leaf, ixN15, ic)
          END IF
        ! for all negative values of recycling_fine_root(ic)
        ELSE
          IF (veg_pool_mt(ix_labile, ixN, ic) > eps8) THEN
            recycling_leaf_n15(ic) = recycling_leaf_n15(ic) &
              &                      + recycling_leaf(ic) / veg_pool_mt(ix_labile, ixN, ic) * veg_pool_mt(ix_labile, ixN15, ic)
          END IF
        END IF
        IF (recycling_fine_root(ic) >= 0.0_wp) THEN  ! for all positive values of recycling_fine_root(ic)
          IF (veg_pool_mt(ix_fine_root, ixN, ic) > eps8) THEN
            recycling_fine_root_n15(ic) = recycling_fine_root_n15(ic) &
              &                           + recycling_fine_root(ic) / veg_pool_mt(ix_fine_root, ixN, ic) &
              &                           * veg_pool_mt(ix_fine_root, ixN15, ic)
          END IF
        ! for all negative values of recycling_fine_root(ic)
        ELSE
          IF (veg_pool_mt(ix_labile, ixN, ic) > eps8) THEN
            recycling_fine_root_n15(ic) = recycling_fine_root_n15(ic) &
              &                           + recycling_fine_root(ic) / veg_pool_mt(ix_labile, ixN, ic) &
              &                           * veg_pool_mt(ix_labile, ixN15, ic)
          END IF
        END IF

        !>2.0 replace overturning P from leaves and roots by new N from the labile pool using the target N:P
        !>
        recycling_leaf(ic)      = (veg_pool_mt(ix_leaf, ixP, ic) - veg_pool_mt(ix_leaf, ixN, ic) / target_np_leaf(ic)) &
          &                       / tau_nutrient_recycling / one_day * dtime
        recycling_fine_root(ic) = (veg_pool_mt(ix_fine_root, ixP, ic) &
         &                       - veg_pool_mt(ix_fine_root, ixN, ic) / target_np_fine_root(ic)) &
         &                       / tau_nutrient_recycling / one_day * dtime
        !>  2.1 if P is addded to leaves or fine root, make sure there is enough labile P to support this
        !>
        IF ((-1.0_wp * (recycling_leaf(ic) + recycling_fine_root(ic))) > eps8) THEN
          hlp1(ic) = -1.0_wp * (1.0_wp/tau_labile/one_day * veg_pool_mt(ix_labile, ixP, ic) * dtime) &
            &        / (recycling_leaf(ic) + recycling_fine_root(ic))
          IF (hlp1(ic) < 1.0_wp) THEN
            recycling_leaf(ic)       = recycling_leaf(ic)      * hlp1(ic)
            recycling_fine_root(ic)  = recycling_fine_root(ic) * hlp1(ic)
          END IF
        END IF
        recycling_leaf_p(ic)        = recycling_leaf_p(ic)        + recycling_leaf(ic)
        recycling_fine_root_p(ic)   = recycling_fine_root_p(ic)   + recycling_fine_root(ic)
      END IF  ! growing_season(:) > test_false_true
    END DO    ! loop over '1,nc'
  END SUBROUTINE calc_nutrient_recycling

  ! ======================================================================================================= !
  !>calculate turnover: recycling and retranslocation
  !>
  SUBROUTINE calc_turnover( &
    & nc                          , &
    & elements_index_map          , &
    & is_element_used             , &
    & lctlib_phenology_type       , &
    & lctlib_growthform           , &
    & lctlib_frac_sapwood_branch  , &
    & lctlib_tau_leaf             , &
    & lctlib_tau_fine_root        , &
    & lctlib_tau_coarse_root      , &
    & lctlib_tau_branch           , &
    & lctlib_tau_sap_wood         , &
    & lctlib_tau_fruit            , &
    & lctlib_tau_seed_litter      , &
    & growing_season              , &
    & lai                         , &
    & target_lai                  , &
    & veg_pool_mt                 , &
    & recycling_leaf_n            , &
    & recycling_leaf_n15          , &
    & recycling_leaf_p            , &
    & recycling_heart_wood_n      , &
    & recycling_heart_wood_n15    , &
    & recycling_heart_wood_p      , &
    & veg_growth_mt               , &
    & veg_litterfall_mt           , &
    & net_growth                  )

    USE mo_veg_constants,         ONLY: max_leaf_shedding_rate, igrass, ITREE
    USE mo_q_pheno_constants,     ONLY: ievergreen, iraingreen, isummergreen, iperennial
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                  INTENT(in)    :: nc                           !< dimensions
    INTEGER,                  INTENT(in)    :: elements_index_map(:)        !< map bgcm element ID -> IDX
    LOGICAL,                  INTENT(in)    :: is_element_used(:)           !< is element in 'elements_index_map' used
    INTEGER,                  INTENT(in)    :: lctlib_phenology_type        !< lctlib parameter
    INTEGER,                  INTENT(in)    :: lctlib_growthform            !< lctlib parameter
    REAL(wp),                 INTENT(in)    :: lctlib_frac_sapwood_branch   !< lctlib parameter
    REAL(wp),                 INTENT(in)    :: lctlib_tau_leaf              !< lctlib parameter
    REAL(wp),                 INTENT(in)    :: lctlib_tau_fine_root         !< lctlib parameter
    REAL(wp),                 INTENT(in)    :: lctlib_tau_coarse_root       !< lctlib parameter
    REAL(wp),                 INTENT(in)    :: lctlib_tau_branch            !< lctlib parameter
    REAL(wp),                 INTENT(in)    :: lctlib_tau_sap_wood          !< lctlib parameter
    REAL(wp),                 INTENT(in)    :: lctlib_tau_fruit             !< lctlib parameter
    REAL(wp),                 INTENT(in)    :: lctlib_tau_seed_litter       !< lctlib parameter
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: growing_season               !< whether the plant is in the growing season
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: lai                          !< the plant's LAI
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: target_lai                   !< the plant's annual maximum LAI
    REAL(wp),                 INTENT(in)    :: veg_pool_mt(:,:,:)           !< bgcm veg_pool: current state of vegetation pools
    REAL(wp), DIMENSION(nc),  INTENT(inout) :: recycling_leaf_n             !< foliar N returning to labile pool
    REAL(wp), DIMENSION(nc),  INTENT(inout) :: recycling_leaf_n15           !< foliar N15 returning to labile pool
    REAL(wp), DIMENSION(nc),  INTENT(inout) :: recycling_leaf_p             !< foliar P returning to labile pool
    REAL(wp), DIMENSION(nc),  INTENT(inout) :: recycling_heart_wood_n       !< woody N returning to labile pool
    REAL(wp), DIMENSION(nc),  INTENT(inout) :: recycling_heart_wood_n15     !< woody N15 returning to labile pool
    REAL(wp), DIMENSION(nc),  INTENT(inout) :: recycling_heart_wood_p       !< woody P returning to labile pool
    REAL(wp),                 INTENT(inout) :: veg_growth_mt(:,:,:)         !< bgcm flux: current growth rate of vegetation pools
    REAL(wp),                 INTENT(inout) :: veg_litterfall_mt(:,:,:)     !< bgcm flux: current litter fall from vegetation pools
    REAL(wp), DIMENSION(nc),  INTENT(out)   :: net_growth                   !< current net growth [mumol C / m2 / s]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ielem                            !< loop over bgcm elements
    INTEGER                     :: ix_elem                          !< index of element in bgcm, used for looping
    REAL(wp), DIMENSION(nc)     :: fturn_leaf                       !< fraction of leaves turning over
    REAL(wp), DIMENSION(nc)     :: fturn                            !< other tissue turnover fraction
    REAL(wp), DIMENSION(nc)     :: retranslocation                  !< retranslocation before shedding
    REAL(wp), DIMENSION(nc)     :: nitrogen_recycling_fract         !< recycling fraction of N
    REAL(wp), DIMENSION(nc)     :: phosphorus_recycling_fract       !< recycling fraction of P
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_turnover'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 init output variable
    !>
    net_growth(:) = 0.0_wp

    !>1.0 foliage turnover, including recycling of nutrients
    !>
    ! flat rate for evergreens, dual rates for deciduous (senescence after the growing season)
    SELECT CASE (lctlib_phenology_type)
    CASE (ievergreen)
      fturn_leaf(:) = 1.0_wp / lctlib_tau_leaf / one_day / one_year * dtime
    CASE (iraingreen, isummergreen)
      WHERE(growing_season(:) > test_false_true)
        fturn_leaf(:) = 0.0_wp              ! herbivory
      ELSEWHERE
        WHERE(lai(:) > eps4)
          fturn_leaf(:) = MAX(eps4, max_leaf_shedding_rate * dtime / one_day * target_lai(:) / lai(:))
        ELSEWHERE
          fturn_leaf(:) = 1.0_wp
        END WHERE
      END WHERE
    CASE (iperennial)
      WHERE(growing_season(:) > test_false_true)
        fturn_leaf(:) = 1.0_wp / lctlib_tau_leaf / one_day / one_year * dtime
      ELSEWHERE
        WHERE(lai(:) > eps4)
          fturn_leaf(:) = MAX(1.0_wp / lctlib_tau_leaf / one_day / one_year * dtime, &
            &             max_leaf_shedding_rate * dtime / one_day * target_lai(:) / lai(:))
        ELSEWHERE
          fturn_leaf(:) = 1.0_wp
        END WHERE
      END WHERE
    END SELECT
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        veg_litterfall_mt(ix_leaf, ix_elem, :) = veg_litterfall_mt(ix_leaf, ix_elem, :) &
          &                                      + veg_pool_mt(ix_leaf, ix_elem, :) * fturn_leaf(:)
      END IF
    END DO
    net_growth(:) = net_growth(:) - veg_pool_mt(ix_leaf, ixC, :) * fturn_leaf(:)

    !>  1.1 Nutrient recycling from foliage turnover
    !>
    nitrogen_recycling_fract(:)   = calc_nitrogen_recycling_fract("leaf")
    phosphorus_recycling_fract(:) = calc_phosphorus_recycling_fract("leaf")
    ! N
    retranslocation(:)                 = nitrogen_recycling_fract(:) * veg_pool_mt(ix_leaf, ixN, :) * fturn_leaf(:)
    recycling_leaf_n(:)                = recycling_leaf_n(:) + retranslocation(:)
    veg_litterfall_mt(ix_leaf, ixN, :) = veg_litterfall_mt(ix_leaf, ixN, :) - retranslocation(:)
    ! N15
    retranslocation(:)                   = nitrogen_recycling_fract(:) * veg_pool_mt(ix_leaf, ixN15, :) * fturn_leaf(:)
    recycling_leaf_n15(:)                = recycling_leaf_n15(:) + retranslocation(:)
    veg_litterfall_mt(ix_leaf, ixN15, :) = veg_litterfall_mt(ix_leaf, ixN15, :) - retranslocation(:)
    ! P
    retranslocation(:)                 = phosphorus_recycling_fract(:) * veg_pool_mt(ix_leaf, ixP, :) * fturn_leaf(:)
    recycling_leaf_p(:)                = recycling_leaf_p(:) + retranslocation(:)
    veg_litterfall_mt(ix_leaf, ixP, :) = veg_litterfall_mt(ix_leaf, ixP, :) - retranslocation(:)

    !>2.0 fine root turnover to litter
    !>
    SELECT CASE (lctlib_growthform)
    CASE (ITREE)
      fturn(:) = 1.0_wp / lctlib_tau_fine_root / one_day / one_year * dtime
    CASE (igrass)
      WHERE(growing_season(:) > test_false_true)
        fturn(:) = 1.0_wp / lctlib_tau_fine_root / one_day / one_year * dtime
      ELSEWHERE
        fturn(:) = fturn_leaf(:)
      END WHERE
    END SELECT
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        veg_litterfall_mt(ix_fine_root, ix_elem, :) = veg_litterfall_mt(ix_fine_root, ix_elem, :) &
          &                                           + veg_pool_mt(ix_fine_root, ix_elem, :) * fturn(:)
      END IF
    END DO
    net_growth(:) = net_growth(:) - veg_pool_mt(ix_fine_root, ixC, :) * fturn(:)

    !>3.0 coarse root turnover to litter
    !>
    fturn(:) = 1.0_wp / lctlib_tau_coarse_root / one_day / one_year * dtime
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        veg_litterfall_mt(ix_coarse_root, ix_elem, :) = veg_litterfall_mt(ix_coarse_root, ix_elem, :) &
          &                                             + veg_pool_mt(ix_coarse_root, ix_elem, :) * fturn(:)
      END IF
    END DO
    net_growth(:) = net_growth(:) - veg_pool_mt(ix_coarse_root, ixC, :) * fturn(:)

    !>4.0 sapwood turnover
    !>

    !>  4.1 litter production from halms (grasses), or branches (trees)
    !>
    SELECT CASE (lctlib_growthform)
    CASE (ITREE)
      fturn(:) = 1.0_wp / lctlib_tau_branch / one_day / one_year * dtime
    CASE (igrass) ! for grases follows leaf phenology
      fturn(:) = fturn_leaf(:)
    ENDSELECT
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        veg_litterfall_mt(ix_sap_wood, ix_elem, :) = veg_litterfall_mt(ix_sap_wood, ix_elem, :) &
          &                                          + lctlib_frac_sapwood_branch &
          &                                          * veg_pool_mt(ix_sap_wood, ix_elem, :) * fturn(:)
      END IF
    END DO
    net_growth(:) = net_growth(:) - lctlib_frac_sapwood_branch * veg_pool_mt(ix_sap_wood, ixC, :) * fturn(:)

    !>  4.2 heart wood formation, including nutrient recycling
    !>
    IF (lctlib_growthform == ITREE) THEN
      fturn(:) = (1.0_wp - lctlib_frac_sapwood_branch) / lctlib_tau_sap_wood / one_day / one_year * dtime
      ! loop over bgcm elements
      DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
        IF (is_element_used(ielem)) THEN
          ix_elem = elements_index_map(ielem)    ! get element index in bgcm
          veg_growth_mt(ix_heart_wood, ix_elem, :) = veg_growth_mt(ix_heart_wood, ix_elem, :) &
            &                                        + veg_pool_mt(ix_sap_wood, ix_elem, :) * fturn(:)
          veg_growth_mt(ix_sap_wood, ix_elem, :)   = veg_growth_mt(ix_sap_wood, ix_elem, :) &
            &                                        - veg_pool_mt(ix_sap_wood, ix_elem, :) * fturn(:)
        END IF
      END DO
      ! Nutrient recycling
      nitrogen_recycling_fract(:)   = calc_nitrogen_recycling_fract("wood")
      phosphorus_recycling_fract(:) = calc_phosphorus_recycling_fract("wood")
      recycling_heart_wood_n(:)     = recycling_heart_wood_n(:)   &
        &                             + veg_pool_mt(ix_sap_wood, ixN, :) * fturn(:) * nitrogen_recycling_fract(:)
      recycling_heart_wood_n15(:)   = recycling_heart_wood_n15(:) &
        &                             + veg_pool_mt(ix_sap_wood, ixN15, :) * fturn(:) * nitrogen_recycling_fract(:)
      recycling_heart_wood_p(:)     = recycling_heart_wood_p(:)   &
        &                             + veg_pool_mt(ix_sap_wood, ixP, :) * fturn(:) * phosphorus_recycling_fract(:)
    END IF

    !>5.0 fruit turnover to seed_bed
    !>
    fturn(:) = 1.0_wp / lctlib_tau_fruit / one_day / one_year * dtime
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        veg_growth_mt(ix_seed_bed, ix_elem, :) = veg_growth_mt(ix_seed_bed, ix_elem, :) &
          &                                      + veg_pool_mt(ix_fruit, ix_elem, :) * fturn(:)
        veg_growth_mt(ix_fruit, ix_elem, :)    = veg_growth_mt(ix_fruit, ix_elem, :) &
          &                                      - veg_pool_mt(ix_fruit, ix_elem, :) * fturn(:)
      END IF
    END DO
    net_growth(:) = net_growth(:) - veg_pool_mt(ix_fruit, ixC, :) * fturn(:)

    !>6.0 seed-bed turnover to litter
    !>
    fturn(:) = 1.0_wp / lctlib_tau_seed_litter / one_day / one_year * dtime
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        veg_litterfall_mt(ix_seed_bed, ix_elem, :) = veg_litterfall_mt(ix_seed_bed, ix_elem, :) &
          &                                          + veg_pool_mt(ix_seed_bed, ix_elem, :) * fturn(:)
      END IF
    END DO
  END SUBROUTINE calc_turnover

  ! ======================================================================================================= !
  !>calculate fractional N recycling from leaves
  !>
  ! @TODO should be made dynamically dependended on N stress later
  ! @TODO veg_turnover:calc_nitrogen_recycling_fract implement functionality other than "set to default"
  ELEMENTAL FUNCTION calc_nitrogen_recycling_fract(char_bgcm_compartment) RESULT (recycling_fract)
    USE mo_veg_constants,         ONLY: resorp_fract_leaf, resorp_fract_wood
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(4), INTENT(in) :: char_bgcm_compartment
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: recycling_fract
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_nitrogen_recycling_fract'
    ! ----------------------------------------------------------------------------------------------------- !

    SELECT CASE(TRIM(char_bgcm_compartment))
    CASE("leaf")
      recycling_fract = resorp_fract_leaf
    CASE("wood")
      recycling_fract = resorp_fract_wood
    ENDSELECT
  END FUNCTION calc_nitrogen_recycling_fract

  ! ======================================================================================================= !
  !>calculate fractional P recycling from leaves
  !>
  ! @TODO should be made dynamically dependended on P stress later
  ! @TODO veg_turnover:calc_phosphorus_recycling_fract implement functionality other than "set to default"
  ELEMENTAL FUNCTION calc_phosphorus_recycling_fract(char_bgcm_compartment) RESULT (recycling_fract)
    USE mo_veg_constants,         ONLY: resorp_fract_leaf, resorp_fract_wood
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(4), INTENT(in) :: char_bgcm_compartment
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: recycling_fract
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_phosphorus_recycling_fract'
    ! ----------------------------------------------------------------------------------------------------- !

    SELECT CASE(TRIM(char_bgcm_compartment))
    CASE("leaf")
      recycling_fract = resorp_fract_leaf
    CASE("wood")
      recycling_fract = resorp_fract_wood
    ENDSELECT
  END FUNCTION calc_phosphorus_recycling_fract

#endif
END MODULE mo_q_veg_turnover
